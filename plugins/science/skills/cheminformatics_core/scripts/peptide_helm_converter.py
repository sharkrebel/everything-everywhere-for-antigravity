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

"""Biopolymer & Sequence conversions: Peptides, FASTA, HELM to RDKit Mol."""

import argparse
import json
import os
from typing import Dict
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


def convert_sequence_to_mol(sequence_str: str) -> Chem.Mol:
  """Converts 1-letter amino acid peptide sequence to RDKit Mol."""
  mol = Chem.MolFromSequence(sequence_str)
  if mol is None:
    raise ValueError(f"Failed to parse peptide sequence: {sequence_str}")
  return mol


def convert_helm_to_mol(helm_str: str) -> Chem.Mol:
  """Converts HELM string (Hierarchical Editing Language for Macromolecules) to RDKit Mol."""
  mol = Chem.MolFromHELM(helm_str)
  if mol is None:
    raise ValueError(f"Failed to parse HELM string: {helm_str}")
  return mol


def analyze_peptide(sequence_str: str) -> Dict[str, any]:
  """Generates molecular properties for a peptide sequence."""
  mol = convert_sequence_to_mol(sequence_str)
  helm_notation = Chem.MolToHELM(mol)
  smiles = Chem.MolToSmiles(mol)
  mw = Descriptors.ExactMolWt(mol)
  formula = rdMolDescriptors.CalcMolFormula(mol)

  return {
      "sequence": sequence_str,
      "length": len(sequence_str),
      "exact_mw": round(mw, 3),
      "formula": formula,
      "helm": helm_notation,
      "smiles": smiles,
  }


def main():
  parser = argparse.ArgumentParser(
      description="Peptide sequence, FASTA, and HELM analyzer."
  )
  parser.add_argument(
      "--sequence", type=str, default=None, help="1-letter amino acid sequence."
  )
  parser.add_argument(
      "--helm",
      type=str,
      default=None,
      help="HELM macromolecule notation string.",
  )
  parser.add_argument(
      "--output",
      "-o",
      type=str,
      default=None,
      help="Output file path to save peptide/HELM analysis as JSON.",
  )
  args = parser.parse_args()

  if args.helm:
    mol = convert_helm_to_mol(args.helm)
    info = {
        "helm": args.helm,
        "exact_mw": round(Descriptors.ExactMolWt(mol), 3),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "smiles": Chem.MolToSmiles(mol),
    }
    if args.output:
      os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
      with open(args.output, "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2)
    print(f"HELM Analysis for '{args.helm}':")
    print(f"  Exact MW:     {info['exact_mw']} g/mol")
    print(f"  Formula:      {info['formula']}")
    print(f"  SMILES:       {info['smiles']}")
  elif args.sequence:
    info = analyze_peptide(args.sequence)
    if args.output:
      os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
      with open(args.output, "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2)
    print(f"Peptide Analysis for sequence '{info['sequence']}':")
    print(f"  Length:       {info['length']} amino acids")
    print(f"  Exact MW:     {info['exact_mw']} g/mol")
    print(f"  Formula:      {info['formula']}")
    print(f"  HELM:         {info['helm']}")
    print(f"  SMILES:       {info['smiles']}")
  else:
    parser.error("Either --sequence or --helm must be provided.")


if __name__ == "__main__":
  main()
