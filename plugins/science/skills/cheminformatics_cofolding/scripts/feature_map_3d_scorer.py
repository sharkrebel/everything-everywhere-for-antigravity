# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "gemmi",
#     "rdkit",
# ]
# ///

"""Calculates 3D Pharmacophore Feature Map (FeatMaps) overlap scores between aligned ligand poses.

Useful for evaluating whether an AlphaFold 3 predicted ligand pose reproduces
the key 3D pharmacophoric interactions (H-bond donors, acceptors, aromatics,
hydrophobes) of an experimental crystal reference ligand.
"""

import argparse
import json
import os
import tempfile
from typing import Dict, Optional

from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, ChemicalFeatures, rdMolAlign
from rdkit.Chem.FeatMaps import FeatMaps

try:
  import gemmi
except ImportError:
  gemmi = None


def load_molecule_3d(
    path: str, smiles: Optional[str] = None, chain: str = "B"
) -> Chem.Mol:
  """Loads 3D molecule from SDF, MOL2, PDB, or mmCIF."""
  if path.endswith(".sdf") or path.endswith(".mol"):
    suppl = Chem.SDMolSupplier(path, removeHs=False)
    if suppl and len(suppl) > 0 and suppl[0] is not None:
      return suppl[0]
    raise ValueError(f"Failed to parse SDF from {path}")
  elif path.endswith(".pdb"):
    mol = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False)
    if mol is not None:
      return mol
    raise ValueError(f"Failed to parse PDB from {path}")
  elif path.endswith(".cif") or path.endswith(".mmcif"):
    if gemmi is None:
      raise ImportError(
          "gemmi is required to parse .cif/.mmcif files in"
          " feature_map_3d_scorer.py."
      )
    doc = gemmi.cif.read(path)
    st = gemmi.make_structure_from_block(doc.sole_block())
    # Extract ligand chain
    new_st = gemmi.Structure()
    new_model = gemmi.Model("1")
    for ch in st[0]:
      if ch.name == chain:
        new_model.add_chain(ch)
    new_st.add_model(new_model)
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_f:
      tmp_pdb = tmp_f.name
    try:
      new_st.write_minimal_pdb(tmp_pdb)
      mol = Chem.MolFromPDBFile(tmp_pdb, removeHs=False, sanitize=False)
    finally:
      if os.path.exists(tmp_pdb):
        os.remove(tmp_pdb)
    if smiles:
      ref_mol = Chem.MolFromSmiles(smiles)
      if ref_mol:
        mol = AllChem.AssignBondOrdersFromTemplate(ref_mol, mol)
    return mol
  else:
    raise ValueError(f"Unsupported file format for {path}")


def compute_3d_feature_map_score(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    fdef_file: Optional[str] = None,
    align: bool = False,
) -> Dict[str, any]:
  """Calculates 3D pharmacophore feature map overlap score between two 3D conformers."""
  fdef_path = fdef_file or os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
  if not os.path.exists(fdef_path):
    raise FileNotFoundError(f"Feature definition file not found: {fdef_path}")

  feat_factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)

  # Optional alignment if poses are not already in common reference frame
  if align:
    pyO3A = rdMolAlign.GetO3A(mol1, mol2)
    rmsd = pyO3A.Align()
  else:
    rmsd = None

  feats1 = feat_factory.GetFeaturesForMol(mol1)
  feats2 = feat_factory.GetFeaturesForMol(mol2)

  fparams = {f.GetFamily(): FeatMaps.FeatMapParams() for f in feats1 + feats2}

  # FeatMap for target mol1
  fmap1 = FeatMaps.FeatMap(
      feats=feats1, weights=[1.0] * len(feats1), params=fparams
  )
  score1_2 = fmap1.ScoreFeats(feats2)
  self_score1 = fmap1.ScoreFeats(feats1)

  # FeatMap for ref mol2
  fmap2 = FeatMaps.FeatMap(
      feats=feats2, weights=[1.0] * len(feats2), params=fparams
  )
  _ = fmap2.ScoreFeats(feats1)
  self_score2 = fmap2.ScoreFeats(feats2)

  # Normalized overlap score (0.0 to 1.0)
  norm_score = (
      float(score1_2 / max(self_score1, self_score2))
      if max(self_score1, self_score2) > 0
      else 0.0
  )

  family_counts1 = {}
  for f in feats1:
    family_counts1[f.GetFamily()] = family_counts1.get(f.GetFamily(), 0) + 1
  family_counts2 = {}
  for f in feats2:
    family_counts2[f.GetFamily()] = family_counts2.get(f.GetFamily(), 0) + 1

  return {
      "featmap_score": round(score1_2, 3),
      "normalized_similarity": round(norm_score, 4),
      "mol1_features_count": len(feats1),
      "mol2_features_count": len(feats2),
      "mol1_feature_families": family_counts1,
      "mol2_feature_families": family_counts2,
      "aligned_rmsd": round(rmsd, 3) if rmsd is not None else "pre-aligned",
  }


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Calculate 3D Pharmacophore Feature Map (FeatMaps) overlap score."
      )
  )
  parser.add_argument(
      "--mol1",
      required=True,
      help="Path to first 3D structure (e.g. AF3 extracted ligand SDF or CIF).",
  )
  parser.add_argument(
      "--mol2",
      required=True,
      help="Path to second 3D structure (e.g. crystal reference SDF or PDB).",
  )
  parser.add_argument(
      "--align",
      action="store_true",
      help="Perform Open3DAlign (O3A) before computing FeatMap score.",
  )
  parser.add_argument(
      "--smiles1",
      default=None,
      help="Optional SMILES for mol1 bond order assignment if CIF.",
  )
  parser.add_argument(
      "--smiles2",
      default=None,
      help="Optional SMILES for mol2 bond order assignment if CIF.",
  )
  parser.add_argument(
      "--output-json",
      default=None,
      help="Optional path to output JSON results.",
  )
  args = parser.parse_args()

  m1 = load_molecule_3d(args.mol1, smiles=args.smiles1)
  m2 = load_molecule_3d(args.mol2, smiles=args.smiles2)

  res = compute_3d_feature_map_score(m1, m2, align=args.align)
  res["mol1_path"] = args.mol1
  res["mol2_path"] = args.mol2

  print(json.dumps(res, indent=2))

  if args.output_json:
    with open(args.output_json, "w") as f:
      json.dump(res, f, indent=2)
    print(f"Results written to {args.output_json}")


if __name__ == "__main__":
  main()
