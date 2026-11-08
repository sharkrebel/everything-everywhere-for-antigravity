# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "pandas",
#     "scipy",
#     "scikit-learn",
#     "matplotlib",
# ]
# ///

"""Computes 2D physicochemical properties, complexity indices, and structural filter alerts."""

import argparse
import json
import os
import sys
from typing import Any, Dict, List, Optional
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, FilterCatalog, GraphDescriptors, Lipinski, QED, rdMolDescriptors


def get_filter_catalog(
    catalogs: Optional[List[str]] = None,
) -> FilterCatalog.FilterCatalog:
  """Builds a unified FilterCatalog for PAINS, BRENK, and NIH alerts."""
  params = FilterCatalog.FilterCatalogParams()
  catalog_map = {
      "pains": FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS,
      "brenk": FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK,
      "nih": FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH,
      "zinc": FilterCatalog.FilterCatalogParams.FilterCatalogs.ZINC,
  }
  selected = catalogs or ["pains", "brenk", "nih"]
  for cat in selected:
    cat_key = cat.lower()
    if cat_key in catalog_map:
      params.AddCatalog(catalog_map[cat_key])
  return FilterCatalog.FilterCatalog(params)


def compute_properties(
    mol: Chem.Mol, catalog: Optional[FilterCatalog.FilterCatalog] = None
) -> Dict[str, Any]:
  """Computes physicochemical descriptors, complexity metrics, and filter catalog alerts."""
  mw = Descriptors.MolWt(mol)
  logp = Crippen.MolLogP(mol)
  tpsa = rdMolDescriptors.CalcTPSA(mol)
  hbd = Lipinski.NumHDonors(mol)
  hba = Lipinski.NumHAcceptors(mol)
  rotb = Lipinski.NumRotatableBonds(mol)
  fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
  heavy_atoms = mol.GetNumHeavyAtoms()
  rings = Lipinski.RingCount(mol)
  aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
  bertz = GraphDescriptors.BertzCT(mol)
  qed_score = QED.qed(mol)

  # Lipinski Ro5 rule violations
  lipinski_violations = sum([
      1 if mw > 500.0 else 0,
      1 if logp > 5.0 else 0,
      1 if hbd > 5 else 0,
      1 if hba > 10 else 0,
  ])

  # Veber rule compliance
  veber_pass = (rotb <= 10) and (tpsa <= 140.0)

  # Structural alerts
  if catalog is None:
    catalog = get_filter_catalog(["pains", "brenk", "nih"])
  matches = catalog.GetMatches(mol)
  alerts = [m.GetDescription() for m in matches]

  return {
      "molecular_weight": round(mw, 3),
      "logp": round(logp, 2),
      "tpsa": round(tpsa, 1),
      "hbd": hbd,
      "hba": hba,
      "rotatable_bonds": rotb,
      "fraction_csp3": round(fsp3, 3),
      "bertz_complexity": round(bertz, 1),
      "heavy_atom_count": heavy_atoms,
      "ring_count": rings,
      "aromatic_ring_count": aromatic_rings,
      "qed": round(qed_score, 4),
      "lipinski_violations": lipinski_violations,
      "lipinski_pass": lipinski_violations <= 1,
      "veber_pass": veber_pass,
      "structural_alerts_count": len(alerts),
      "alerts": alerts,
  }


def main():
  parser = argparse.ArgumentParser(
      description="Compute molecular descriptors and filter alerts."
  )
  parser.add_argument(
      "--smiles", type=str, required=True, help="Input SMILES string."
  )
  parser.add_argument(
      "--json", action="store_true", help="Output in JSON format."
  )
  parser.add_argument(
      "--output",
      "-o",
      type=str,
      default=None,
      help="Output file path to save descriptors as JSON.",
  )
  args = parser.parse_args()

  smiles = args.smiles
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    sys.stderr.write(f"Error: Invalid SMILES string '{smiles}'\n")
    sys.exit(1)

  props = compute_properties(mol)
  if args.output:
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
      json.dump(props, f, indent=2)
  if args.json:
    print(json.dumps(props, indent=2))
  else:
    print(f"Properties for SMILES: {smiles}")
    for k, v in props.items():
      print(f"  {k:<25}: {v}")


if __name__ == "__main__":
  main()
