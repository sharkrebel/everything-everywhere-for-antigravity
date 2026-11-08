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

"""Performs MaxMin diversity selection on chemical libraries."""

import argparse
from typing import List
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.SimDivFilters import rdSimDivPickers


def pick_diverse_subset(
    smiles_list: List[str],
    num_to_pick: int,
    seed: int = 42,
    fp_radius: int = 2,
    fp_size: int = 2048,
) -> List[int]:
  """Selects a maximally diverse subset using the C++ accelerated LazyBitVectorPick algorithm.

  Args:
      smiles_list: List of SMILES strings.
      num_to_pick: Number of diverse compounds to select.
      seed: Random seed for initial pick.
      fp_radius: Morgan fingerprint radius.
      fp_size: Fingerprint bit size.

  Returns:
      List of selected indices in original smiles_list.
  """
  mols = []
  valid_indices = []
  for idx, s in enumerate(smiles_list):
    m = Chem.MolFromSmiles(str(s)) if pd.notna(s) else None
    if m is not None:
      mols.append(m)
      valid_indices.append(idx)

  if not mols:
    return []

  if num_to_pick > len(mols):
    num_to_pick = len(mols)

  gen = rdFingerprintGenerator.GetMorganGenerator(
      radius=fp_radius, fpSize=fp_size
  )
  fps = [gen.GetFingerprint(m) for m in mols]

  picker = rdSimDivPickers.MaxMinPicker()
  # LazyBitVectorPick runs directly on bit vectors with C++ Tanimoto acceleration
  picked_local_indices = list(
      picker.LazyBitVectorPick(fps, len(fps), num_to_pick, seed=seed)
  )
  return [valid_indices[i] for i in picked_local_indices]


def pick_diverse_dataframe(
    df: pd.DataFrame,
    num_to_pick: int,
    smiles_col: str = "SMILES",
    seed: int = 42,
) -> pd.DataFrame:
  """In-Memory DataFrame Diversity Selection: Returns a DataFrame containing only the diverse subset."""
  picked_indices = pick_diverse_subset(
      df[smiles_col].tolist(), num_to_pick=num_to_pick, seed=seed
  )
  return df.iloc[picked_indices].copy()


def main():
  parser = argparse.ArgumentParser(
      description="MaxMin Diversity Picker for compounds."
  )
  parser.add_argument(
      "--input_csv", type=str, required=True, help="Input CSV path."
  )
  parser.add_argument(
      "--output_csv",
      type=str,
      default=None,
      help="Optional output CSV path for selected compounds.",
  )
  parser.add_argument(
      "--smiles_col", type=str, default="SMILES", help="SMILES column name."
  )
  parser.add_argument(
      "--num_picks", type=int, default=4, help="Number of compounds to select."
  )
  args = parser.parse_args()

  df = pd.read_csv(args.input_csv)
  diverse_df = pick_diverse_dataframe(
      df, num_to_pick=args.num_picks, smiles_col=args.smiles_col
  )
  print(
      f"Selected {len(diverse_df)} diverse compounds from {len(df)} total"
      " entries."
  )

  if args.output_csv:
    diverse_df.to_csv(args.output_csv, index=False)
    print(f"Saved selected diverse compounds to {args.output_csv}")


if __name__ == "__main__":
  main()
