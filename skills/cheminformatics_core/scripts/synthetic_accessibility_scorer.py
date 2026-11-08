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

"""Calculates Ertl & Schuffenhauer Synthetic Accessibility Score (SA_Score 1-10)."""

import argparse
import json
import os
import sys
from typing import Dict, List, Union
import pandas as pd
from rdkit import Chem, RDConfig

# Load sascorer from RDKit contrib directory
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer


def calculate_sa_score(mol_or_smiles: Union[str, Chem.Mol]) -> float:
  """Calculates Synthetic Accessibility Score (1 = very easy to synthesize, 10 = extremely difficult).

  Accepts either a SMILES string or an already-parsed RDKit Mol, so callers
  that have a Mol in hand (e.g. the shared descriptor engine) do not have to
  round-trip through SMILES.
  """
  if isinstance(mol_or_smiles, Chem.Mol):
    mol = mol_or_smiles
  else:
    mol = Chem.MolFromSmiles(mol_or_smiles)
    if mol is None:
      raise ValueError(f"Invalid SMILES: {mol_or_smiles}")
  if mol is None:
    raise ValueError("calculate_sa_score received a null molecule")
  return float(sascorer.calculateScore(mol))


def score_dataset(smiles_list: List[str]) -> List[Dict[str, any]]:
  """Calculates SA_Scores for a list of SMILES."""
  results = []
  for s in smiles_list:
    try:
      score = calculate_sa_score(s)
      results.append({
          "smiles": s,
          "sa_score": round(score, 3),
          "synthesizability": (
              "Easy"
              if score < 4.0
              else ("Moderate" if score < 6.0 else "Difficult")
          ),
      })
    except Exception as e:
      results.append({"smiles": s, "error": str(e)})
  return results


def main():
  parser = argparse.ArgumentParser(
      description="Synthetic Accessibility (SA_Score) calculator."
  )
  parser.add_argument(
      "--smiles", nargs="+", required=True, help="List of SMILES to score."
  )
  parser.add_argument(
      "--output",
      "-o",
      type=str,
      default=None,
      help="Output file path to save SA scores (.json or .csv).",
  )
  args = parser.parse_args()

  scores = score_dataset(args.smiles)
  if args.output:
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    if args.output.endswith(".csv"):
      pd.DataFrame(scores).to_csv(args.output, index=False)
    else:
      with open(args.output, "w", encoding="utf-8") as f:
        json.dump(scores, f, indent=2)
  print(f"Synthetic Accessibility Scores for {len(scores)} compounds:")
  for r in scores:
    if "sa_score" in r:
      print(
          f"  {r['smiles']:<35} -> SA Score: {r['sa_score']:5.2f}"
          f" [{r['synthesizability']}]"
      )


if __name__ == "__main__":
  main()
