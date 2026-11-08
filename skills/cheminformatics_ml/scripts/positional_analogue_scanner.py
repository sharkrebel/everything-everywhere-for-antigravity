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

"""Performs Positional Analogue Scanning (PAS) for systematic substituent exploration."""

import argparse
from typing import List, Optional
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors

SCAN_PRESETS = {
    "fluorine": "[c;H1:1]>>[c:1]F",
    "chlorine": "[c;H1:1]>>[c:1]Cl",
    "methyl": "[c;H1:1]>>[c:1]C",
    "hydroxyl": "[c;H1:1]>>[c:1]O",
    "methoxy": "[c;H1:1]>>[c:1]OC",
    "trifluoromethyl": "[c;H1:1]>>[c:1]C(F)(F)F",
    "nitrogen_walk": "[c;H1:1]>>[n:1]",
}


def scan_positional_analogues(
    smiles: str,
    scan_type: str = "fluorine",
    custom_smarts: Optional[str] = None,
) -> List[str]:
  """Applies systematic positional scanning across all eligible aromatic C-H positions."""
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    raise ValueError(f"Invalid SMILES string: {smiles}")

  rxn_smarts = custom_smarts or SCAN_PRESETS.get(scan_type.lower())
  if not rxn_smarts:
    raise ValueError(
        f"Unknown scan type '{scan_type}'. Available:"
        f" {list(SCAN_PRESETS.keys())}"
    )

  rxn = AllChem.ReactionFromSmarts(rxn_smarts)
  analogues = set()

  for outcome in rxn.RunReactants((mol,)):
    prod = outcome[0]
    try:
      Chem.SanitizeMol(prod)
      prod_smi = Chem.MolToSmiles(prod)
      if prod_smi != smiles:
        analogues.add(prod_smi)
    except Exception:
      continue

  return sorted(list(analogues))


def scan_positional_analogues_df(
    smiles: str,
    scan_type: str = "fluorine",
    custom_smarts: Optional[str] = None,
) -> pd.DataFrame:
  """In-Memory DataFrame Positional Scanning: Computes analogues with property deltas (MW, cLogP, TPSA)."""
  parent_mol = Chem.MolFromSmiles(smiles)
  if parent_mol is None:
    raise ValueError(f"Invalid SMILES: {smiles}")

  p_mw = Descriptors.ExactMolWt(parent_mol)
  p_clogp = Crippen.MolLogP(parent_mol)
  _ = rdMolDescriptors.CalcTPSA(parent_mol)

  analogues = scan_positional_analogues(
      smiles, scan_type=scan_type, custom_smarts=custom_smarts
  )
  records = []

  for idx, a_smi in enumerate(analogues):
    a_mol = Chem.MolFromSmiles(a_smi)
    mw = round(Descriptors.ExactMolWt(a_mol), 1)
    clogp = round(Crippen.MolLogP(a_mol), 2)
    tpsa = round(rdMolDescriptors.CalcTPSA(a_mol), 1)

    records.append({
        "Scan_Type": scan_type,
        "Analogue_ID": f"{scan_type.capitalize()}_Scan_{idx+1}",
        "Analogue_SMILES": a_smi,
        "MW": mw,
        "cLogP": clogp,
        "TPSA": tpsa,
        "Delta_MW": round(mw - p_mw, 1),
        "Delta_cLogP": round(clogp - p_clogp, 2),
    })

  return pd.DataFrame(records)


def main():
  parser = argparse.ArgumentParser(
      description="Generate systematic positional analogues (PAS)."
  )
  parser.add_argument(
      "--smiles", type=str, required=True, help="Input parent SMILES string."
  )
  parser.add_argument(
      "--scan_type",
      type=str,
      default="fluorine",
      choices=[
          "fluorine",
          "chlorine",
          "methyl",
          "hydroxyl",
          "methoxy",
          "trifluoromethyl",
          "nitrogen_walk",
      ],
      help="Substituent walk type.",
  )
  parser.add_argument(
      "--output_csv",
      "-o",
      type=str,
      default=None,
      help="Output CSV path to save analogue library.",
  )
  args = parser.parse_args()

  df = scan_positional_analogues_df(args.smiles, scan_type=args.scan_type)
  if args.output_csv:
    df.to_csv(args.output_csv, index=False)
    print(f"Saved {len(df)} positional analogues to {args.output_csv}")

  print(
      f"Positional Analogue Scan ({args.scan_type}) for '{args.smiles}'"
      f" ({len(df)} analogues generated):"
  )
  for _, row in df.iterrows():
    print(
        f"  {row['Analogue_ID']}: {row['Analogue_SMILES']:30s} | MW:"
        f" {row['MW']} (Δ {row['Delta_MW']:+0.1f}) | cLogP: {row['cLogP']} (Δ"
        f" {row['Delta_cLogP']:+0.2f})"
    )


if __name__ == "__main__":
  main()
