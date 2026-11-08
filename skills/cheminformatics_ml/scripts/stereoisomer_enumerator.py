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

"""Enumerates undefined stereocenters and double-bond stereoisomers."""

import argparse
from typing import List
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions


def enumerate_isomers(
    smiles: str,
    only_unassigned: bool = True,
    max_isomers: int = 32,
) -> List[str]:
  """Generates all possible stereoisomers for undefined chiral centers and double bonds."""
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    raise ValueError(f"Invalid SMILES: {smiles}")

  opts = StereoEnumerationOptions(
      onlyUnassigned=only_unassigned,
      maxIsomers=max_isomers,
      unique=True,
  )

  isomers = list(EnumerateStereoisomers(mol, options=opts))
  return [Chem.MolToSmiles(iso, isomericSmiles=True) for iso in isomers]


def enumerate_isomers_df(
    smiles: str,
    only_unassigned: bool = True,
    max_isomers: int = 32,
) -> pd.DataFrame:
  """In-Memory DataFrame Stereoisomer Enumeration: Returns isomers with 3D chiral descriptors."""
  isomers = enumerate_isomers(
      smiles, only_unassigned=only_unassigned, max_isomers=max_isomers
  )
  records = []
  for idx, iso_smi in enumerate(isomers):
    mol = Chem.MolFromSmiles(iso_smi)
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    records.append({
        "Isomer_ID": f"Stereoisomer_{idx+1}",
        "Isomeric_SMILES": iso_smi,
        "Chiral_Centers": str(centers),
        "MW": round(Descriptors.ExactMolWt(mol), 1),
        "cLogP": round(Crippen.MolLogP(mol), 2),
    })
  return pd.DataFrame(records)


def main():
  parser = argparse.ArgumentParser(
      description="Enumerate stereoisomers for undefined chiral centers."
  )
  parser.add_argument(
      "--smiles",
      type=str,
      required=True,
      help="Input SMILES with undefined stereo.",
  )
  parser.add_argument(
      "--max_isomers",
      type=int,
      default=32,
      help="Max stereoisomers to generate.",
  )
  parser.add_argument(
      "--output_csv",
      "-o",
      type=str,
      default=None,
      help="Output CSV path to save stereoisomers.",
  )
  args = parser.parse_args()

  df = enumerate_isomers_df(args.smiles, max_isomers=args.max_isomers)
  if args.output_csv:
    df.to_csv(args.output_csv, index=False)
    print(f"Saved {len(df)} stereoisomers to {args.output_csv}")

  print(f"Enumerated {len(df)} stereoisomers for '{args.smiles}':")
  for _, row in df.iterrows():
    print(
        f"  {row['Isomer_ID']}: {row['Isomeric_SMILES']:30s} | Centers:"
        f" {row['Chiral_Centers']}"
    )


if __name__ == "__main__":
  main()
