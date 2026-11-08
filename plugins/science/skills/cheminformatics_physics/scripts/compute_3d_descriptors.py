# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "pandas",
#     "numpy",
# ]
# ///

"""3D conformer ensemble generation and 3D shape/inertia descriptors (PMI, Rg, spherocity, asphericity)."""

import argparse
import gzip
import os
import sys
import time
from typing import Any, Dict, List, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors3D, Lipinski

# Ensure line-buffered stdout for real-time progress logging
if hasattr(sys.stdout, "reconfigure"):
  sys.stdout.reconfigure(line_buffering=True)


def generate_macrocycle_conformers(
    mol: Chem.Mol,
    num_confs: int = 20,
    random_seed: int = 42,
) -> Tuple[Chem.Mol, List[int]]:
  """Generates 3D conformers using ETKDGv3 with small-ring and macrocycle torsion rules."""
  mol_3d = Chem.AddHs(mol)
  params = AllChem.ETKDGv3()
  params.randomSeed = random_seed
  params.useMacrocycleTorsions = True  # Critical for macrocycles
  params.useSmallRingTorsions = True
  params.numThreads = 0  # Use all available CPU cores

  cids = list(
      AllChem.EmbedMultipleConfs(mol_3d, numConfs=num_confs, params=params)
  )
  if not cids:
    params.useRandomCoords = True
    cids = list(
        AllChem.EmbedMultipleConfs(mol_3d, numConfs=num_confs, params=params)
    )

  # Optimize with MMFF94s (using all CPU cores)
  if AllChem.MMFFHasAllMoleculeParams(mol_3d):
    AllChem.MMFFOptimizeMoleculeConfs(
        mol_3d, mmffVariant="MMFF94s", numThreads=0
    )
  else:
    AllChem.UFFOptimizeMoleculeConfs(mol_3d, numThreads=0)

  return mol_3d, cids


def compute_3d_shape_descriptors(
    mol_3d: Chem.Mol, conf_id: int = -1
) -> Dict[str, float]:
  """Calculates 3D shape, inertia, and sphericity descriptors."""
  try:
    return {
        "pmi1": round(float(Descriptors3D.PMI1(mol_3d, confId=conf_id)), 3),
        "pmi2": round(float(Descriptors3D.PMI2(mol_3d, confId=conf_id)), 3),
        "pmi3": round(float(Descriptors3D.PMI3(mol_3d, confId=conf_id)), 3),
        "radius_of_gyration": round(
            float(Descriptors3D.RadiusOfGyration(mol_3d, confId=conf_id)), 3
        ),
        "spherocity_index": round(
            float(Descriptors3D.SpherocityIndex(mol_3d, confId=conf_id)), 3
        ),
        "asphericity": round(
            float(Descriptors3D.Asphericity(mol_3d, confId=conf_id)), 3
        ),
        "inertial_shape_factor": round(
            float(Descriptors3D.InertialShapeFactor(mol_3d, confId=conf_id)), 3
        ),
    }
  except Exception:
    return {
        "pmi1": 0.0,
        "pmi2": 0.0,
        "pmi3": 0.0,
        "radius_of_gyration": 0.0,
        "spherocity_index": 0.0,
        "asphericity": 0.0,
        "inertial_shape_factor": 0.0,
    }


def _has_3d_coordinates(mol: Chem.Mol) -> bool:
  """Returns True if `mol` carries a 3D conformer."""
  if mol.GetNumConformers() == 0:
    return False
  conf = mol.GetConformer(0)
  if conf.Is3D():
    return True
  return any(
      abs(conf.GetAtomPosition(i).z) > 1e-4 for i in range(mol.GetNumAtoms())
  )


def _ring_and_rotatable_metadata(heavy_mol: Chem.Mol) -> Dict[str, Any]:
  """Computes macrocycle and rotatable-bond metadata on a heavy-atom graph."""
  atom_rings = heavy_mol.GetRingInfo().AtomRings()
  max_ring = max((len(r) for r in atom_rings), default=0)
  rotb = Lipinski.NumRotatableBonds(heavy_mol)
  return {
      "is_macrocycle": max_ring >= 12,
      "max_ring_size": max_ring,
      "rot_entropy_penalty_kcal": round(rotb * 0.65, 2),
  }


def compute_3d_descriptors_from_structure_file(
    input_path: str,
    num_confs: int = 15,
) -> pd.DataFrame:
  """Computes 3D shape descriptors from an SDF, SDF.GZ, or MOL file."""
  if not os.path.exists(input_path):
    raise FileNotFoundError(f"Input structure file not found: {input_path}")

  lower = input_path.lower()
  if lower.endswith(".sdf.gz"):
    with gzip.open(input_path, "rb") as gz:
      mols = list(Chem.ForwardSDMolSupplier(gz, removeHs=False))
  elif lower.endswith(".sdf"):
    mols = list(Chem.SDMolSupplier(input_path, removeHs=False))
  elif lower.endswith((".mol", ".mdl")):
    mols = [Chem.MolFromMolFile(input_path, removeHs=False)]
  else:
    raise ValueError(
        f"Unsupported structure format for {input_path}; expected .sdf,"
        " .sdf.gz, .mol, or .csv."
    )

  if not any(m is not None for m in mols):
    raise ValueError(f"No valid molecules could be parsed from {input_path}.")

  rows: List[Dict[str, Any]] = []
  for idx, mol in enumerate(mols, 1):
    if mol is None:
      continue
    name = (
        mol.GetProp("_Name").strip()
        if mol.HasProp("_Name") and mol.GetProp("_Name").strip()
        else f"mol_{idx}"
    )
    heavy_mol = Chem.RemoveHs(mol)
    smi = Chem.MolToSmiles(heavy_mol)
    meta = _ring_and_rotatable_metadata(heavy_mol)

    if _has_3d_coordinates(mol):
      mol_3d = Chem.AddHs(mol, addCoords=True)
      conf_id = mol_3d.GetConformer(0).GetId()
    else:
      mol_3d, cids = generate_macrocycle_conformers(
          heavy_mol, num_confs=num_confs
      )
      conf_id = cids[0] if cids else -1

    desc = compute_3d_shape_descriptors(mol_3d, conf_id=conf_id)
    row: Dict[str, Any] = {"compound_id": name, "SMILES": smi}
    row.update(desc)
    row.update(meta)
    rows.append(row)

  return pd.DataFrame(rows)


def main():
  parser = argparse.ArgumentParser(
      description="3D conformer and shape descriptor calculation."
  )
  parser.add_argument(
      "--smiles", type=str, default=None, help="Single input SMILES."
  )
  parser.add_argument(
      "--input_csv",
      "--input-csv",
      dest="input_csv",
      type=str,
      default=None,
      help="Path to input CSV for batch processing.",
  )
  parser.add_argument(
      "--input",
      "--input-sdf",
      "--input_sdf",
      dest="input_path",
      type=str,
      default=None,
      help="Path to input 3D SDF/SDF.GZ/MOL file or batch CSV.",
  )
  parser.add_argument(
      "--smiles_col",
      "--smiles-col",
      dest="smiles_col",
      type=str,
      default="SMILES",
      help="Column name for SMILES in CSV.",
  )
  parser.add_argument(
      "--output_csv",
      "--output-csv",
      "--output",
      dest="output_csv",
      type=str,
      default=None,
      help="Path to save enriched output CSV.",
  )
  parser.add_argument(
      "--num_confs",
      "--num-confs",
      dest="num_confs",
      type=int,
      default=15,
      help="Number of conformers per molecule.",
  )
  args = parser.parse_args()

  if (
      args.input_path
      and args.input_path.lower().endswith(".csv")
      and not args.input_csv
  ):
    args.input_csv = args.input_path
    args.input_path = None

  if args.input_path:
    start_time = time.time()
    df_out = compute_3d_descriptors_from_structure_file(
        args.input_path, num_confs=args.num_confs
    )
    stem = (
        args.input_path[:-7]
        if args.input_path.lower().endswith(".sdf.gz")
        else os.path.splitext(args.input_path)[0]
    )
    out_path = (
        args.output_csv
        if args.output_csv
        else f"{stem}_with_3d_descriptors.csv"
    )
    out_dir = os.path.dirname(os.path.abspath(out_path))
    if out_dir:
      os.makedirs(out_dir, exist_ok=True)
    df_out.to_csv(out_path, index=False)
    total_time = time.time() - start_time
    print(
        f"Successfully computed 3D descriptors for {len(df_out)}"
        f" structure(s) from {args.input_path} in {total_time:.2f}s.",
        flush=True,
    )
    print(f"Enriched CSV saved to: {os.path.abspath(out_path)}", flush=True)
  elif args.input_csv:
    df = pd.read_csv(args.input_csv)
    s_col = args.smiles_col if args.smiles_col in df.columns else None
    if not s_col:
      for c in ["SMILES", "smiles", "canonical_smiles", "Smiles"]:
        if c in df.columns:
          s_col = c
          break
    if not s_col:
      raise ValueError(
          f"SMILES column '{args.smiles_col}' not found in {args.input_csv}"
      )

    print(
        f"Loaded {len(df)} compounds from {args.input_csv}. Computing"
        " 3D shape descriptors...",
        flush=True,
    )
    start_time = time.time()
    results = []

    for idx, smi in enumerate(df[s_col], 1):
      t0 = time.time()
      mol = Chem.MolFromSmiles(str(smi))
      if mol is None:
        results.append({
            "is_macrocycle": False,
            "max_ring_size": 0,
            "rot_entropy_penalty_kcal": 0.0,
            "pmi1": 0.0,
            "pmi2": 0.0,
            "pmi3": 0.0,
            "radius_of_gyration": 0.0,
            "spherocity_index": 0.0,
            "asphericity": 0.0,
            "inertial_shape_factor": 0.0,
        })
        continue

      ri = mol.GetRingInfo()
      atom_rings = ri.AtomRings()
      max_ring = max([len(r) for r in atom_rings]) if atom_rings else 0
      is_macro = max_ring >= 12
      rotb = Lipinski.NumRotatableBonds(mol)
      rot_penalty = round(rotb * 0.65, 2)  # ~0.65 kcal/mol per rotatable bond

      mol_3d, cids = generate_macrocycle_conformers(
          mol, num_confs=args.num_confs
      )
      best_cid = cids[0] if cids else -1
      desc = compute_3d_shape_descriptors(mol_3d, conf_id=best_cid)
      desc["is_macrocycle"] = is_macro
      desc["max_ring_size"] = max_ring
      desc["rot_entropy_penalty_kcal"] = rot_penalty
      results.append(desc)

      elapsed = time.time() - t0
      print(
          f"  [{idx}/{len(df)}] Processed in {elapsed:.2f}s (Macrocycle:"
          f" {is_macro}, Ring: {max_ring}, RotB: {rotb})",
          flush=True,
      )

    desc_df = pd.DataFrame(results)
    for col in desc_df.columns:
      df[col] = desc_df[col]

    out_path = (
        args.output_csv
        if args.output_csv
        else os.path.splitext(args.input_csv)[0] + "_with_3d_descriptors.csv"
    )
    out_dir = os.path.dirname(os.path.abspath(out_path))
    if out_dir:
      os.makedirs(out_dir, exist_ok=True)
    df.to_csv(out_path, index=False)
    total_time = time.time() - start_time
    print(
        f"\nSuccessfully computed 3D descriptors for {len(df)}"
        f" compounds in {total_time:.2f}s.",
        flush=True,
    )
    print(f"Enriched CSV saved to: {os.path.abspath(out_path)}", flush=True)
  elif args.smiles:
    mol = Chem.MolFromSmiles(args.smiles)
    if mol is None:
      sys.stderr.write("Invalid SMILES\n")
      sys.exit(1)

    mol_3d, cids = generate_macrocycle_conformers(mol, num_confs=args.num_confs)
    print(
        f"Generated {len(cids)} conformers with useMacrocycleTorsions=True.",
        flush=True,
    )

    if cids:
      props = compute_3d_shape_descriptors(mol_3d, conf_id=cids[0])
      # Report the same fields the batch path reports, so the two modes
      # agree and the reference guide's field list holds in both.
      atom_rings = mol.GetRingInfo().AtomRings()
      max_ring = max([len(r) for r in atom_rings]) if atom_rings else 0
      is_macro = max_ring >= 12
      rotb = Lipinski.NumRotatableBonds(mol)
      props["is_macrocycle"] = is_macro
      props["max_ring_size"] = max_ring
      props["rot_entropy_penalty_kcal"] = round(rotb * 0.65, 2)
      print(f"3D Descriptors for Conformer {cids[0]}:", flush=True)
      for k, v in props.items():
        print(f"  {k:<22}: {v}", flush=True)

      if args.output_csv:
        out_dir = os.path.dirname(os.path.abspath(args.output_csv))
        if out_dir:
          os.makedirs(out_dir, exist_ok=True)
        row = {"SMILES": args.smiles}
        row.update(props)
        pd.DataFrame([row]).to_csv(args.output_csv, index=False)
        print(
            f"Descriptors saved to: {os.path.abspath(args.output_csv)}",
            flush=True,
        )
    else:
      sys.stderr.write(
          "ERROR: no conformer could be embedded for this molecule; "
          "no descriptors computed.\n"
      )
      sys.exit(1)
  else:
    print(
        "Please provide --input (SDF/MOL/CSV), --input_csv, or --smiles.",
        flush=True,
    )
    sys.exit(1)


if __name__ == "__main__":
  main()
