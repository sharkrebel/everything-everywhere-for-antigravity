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

"""Extracts and analyzes ring systems and frequency distributions across chemical libraries."""

import argparse
from collections import Counter
import json
import os
from typing import List, Tuple
import pandas as pd
from rdkit import Chem


def extract_ring_systems(mol: Chem.Mol) -> List[str]:
  """Extracts all connected ring systems (including fused and bridged systems) from a molecule."""
  ring_info = mol.GetRingInfo()
  atom_rings = ring_info.AtomRings()
  if not atom_rings:
    return []

  # Find connected components of ring atoms
  ring_atoms_set = set()
  for r in atom_rings:
    ring_atoms_set.update(r)

  # Group connected ring atoms
  visited = set()
  ring_systems = []

  for start_atom_idx in ring_atoms_set:
    if start_atom_idx in visited:
      continue
    component = set()
    queue = [start_atom_idx]
    visited.add(start_atom_idx)

    while queue:
      curr = queue.pop(0)
      component.add(curr)
      curr_atom = mol.GetAtomWithIdx(curr)
      for neighbor in curr_atom.GetNeighbors():
        n_idx = neighbor.GetIdx()
        if n_idx in ring_atoms_set and n_idx not in visited:
          visited.add(n_idx)
          queue.append(n_idx)

    # Extract submolecule for this ring system
    submol = Chem.PathToSubmol(
        mol,
        [
            b.GetIdx()
            for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in component
            and b.GetEndAtomIdx() in component
        ],
    )
    if submol and submol.GetNumAtoms() > 0:
      try:
        Chem.SanitizeMol(submol)
        ring_systems.append(Chem.MolToSmiles(submol))
      except Exception:
        continue

  return ring_systems


def analyze_library_rings(
    smiles_list: List[str],
) -> List[Tuple[str, int, float]]:
  """Analyzes ring system frequency across a library."""
  ring_counter = Counter()
  total_mols = len(smiles_list)

  for s in smiles_list:
    mol = Chem.MolFromSmiles(s)
    if mol is None:
      continue
    rings = extract_ring_systems(mol)
    ring_counter.update(rings)

  results = []
  for ring_smi, count in ring_counter.most_common():
    pct = (count / total_mols) * 100.0
    results.append((ring_smi, count, round(pct, 2)))

  return results


def main():
  parser = argparse.ArgumentParser(
      description="Ring system extraction and frequency analyzer."
  )
  parser.add_argument("--smiles", nargs="+", help="List of SMILES to analyze.")
  parser.add_argument(
      "--input_csv", type=str, help="Input CSV file containing a SMILES column."
  )
  parser.add_argument(
      "--smiles_col",
      type=str,
      default="SMILES",
      help="Name of the SMILES column when --input_csv is used.",
  )
  parser.add_argument(
      "--output",
      "-o",
      type=str,
      default=None,
      help="Output path to save full ring distribution (.csv or .json).",
  )
  args = parser.parse_args()

  smiles = []
  if args.smiles:
    smiles.extend(args.smiles)
  if args.input_csv:
    df = pd.read_csv(args.input_csv)
    col = (
        args.smiles_col
        if args.smiles_col in df.columns
        else next(
            (c for c in df.columns if c.lower() == "smiles"), df.columns[0]
        )
    )
    smiles.extend(df[col].dropna().astype(str).tolist())

  if not smiles:
    parser.error("Either --smiles or --input_csv is required.")

  rings = analyze_library_rings(smiles)
  if args.output:
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    if args.output.endswith(".json"):
      with open(args.output, "w", encoding="utf-8") as f:
        json.dump(
            [
                {"Ring_SMILES": r, "Count": c, "Frequency_Pct": p}
                for r, c, p in rings
            ],
            f,
            indent=2,
        )
    else:
      df_out = pd.DataFrame(
          rings, columns=["Ring_SMILES", "Count", "Frequency_Pct"]
      )
      df_out.to_csv(args.output, index=False)
    print(
        f"Saved {len(rings)} unique ring systems across {len(smiles)} molecules"
        f" to {args.output}"
    )
  else:
    print(
        f"Extracted {len(rings)} unique ring systems across {len(smiles)}"
        " molecules:"
    )
    for i, (ring_smi, count, pct) in enumerate(rings):
      print(f"  Ring {i+1:2d}: {ring_smi:<25} Count: {count:2d} ({pct:.1f}%)")


if __name__ == "__main__":
  main()
