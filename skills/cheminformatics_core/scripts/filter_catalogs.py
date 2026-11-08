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

"""Medicinal Chemistry Structural Alert Filtering & Multi-Parameter Optimization (MPO).

Provides high-speed C++ screening of:
1. Structural Alerts: PAINS (A/B/C), Brenk, NIH, and ZINC filter catalogs.
2. Multi-Parameter Optimization (MPO): Pfizer CNS MPO score (0.0 to 6.0).
3. Drug-Likeness Rules: Lipinski Rule of 5, Fragment Rule of 3, and Veber rules.
"""

import argparse
import json
import os
import sys
from typing import Dict, List, Optional, Tuple
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, PandasTools, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

# Canonical CNS MPO implementation (single source of truth).
_script_dir = os.path.dirname(os.path.abspath(__file__))
if _script_dir not in sys.path:
  sys.path.insert(0, _script_dir)

from descriptors import calculate_cns_mpo  # noqa: E402


def get_filter_catalog(
    *,
    include_pains: bool = True,
    include_brenk: bool = True,
    include_nih: bool = True,
    include_zinc: bool = False,
) -> FilterCatalog:
  """Builds an RDKit FilterCatalog with the specified standard alert catalogs."""
  params = FilterCatalogParams()
  if include_pains:
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
  if include_brenk:
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
  if include_nih:
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
  if include_zinc:
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.ZINC)
  return FilterCatalog(params)


def calculate_rule_of_5(mol: Chem.Mol) -> Tuple[int, List[str]]:
  """Evaluates Lipinski's Rule of 5: returns (violation_count, list_of_violations)."""
  if mol is None:
    return 4, ["Invalid Molecule"]
  violations = []
  mw = Descriptors.MolWt(mol)
  clogp = Crippen.MolLogP(mol)
  hbd = Lipinski.NumHDonors(mol)
  hba = Lipinski.NumHAcceptors(mol)

  if mw > 500.0:
    violations.append(f"MW > 500 ({mw:.1f})")
  if clogp > 5.0:
    violations.append(f"cLogP > 5.0 ({clogp:.2f})")
  if hbd > 5:
    violations.append(f"HBD > 5 ({hbd})")
  if hba > 10:
    violations.append(f"HBA > 10 ({hba})")

  return len(violations), violations


def calculate_rule_of_3(mol: Chem.Mol) -> Tuple[bool, List[str]]:
  """Evaluates Congreve's Rule of 3 for Fragment screening (MW<=300, cLogP<=3, HBD<=3, HBA<=3, RotB<=3)."""
  if mol is None:
    return False, ["Invalid Molecule"]
  violations = []
  mw = Descriptors.MolWt(mol)
  clogp = Crippen.MolLogP(mol)
  hbd = Lipinski.NumHDonors(mol)
  hba = Lipinski.NumHAcceptors(mol)
  rotb = Lipinski.NumRotatableBonds(mol)

  if mw > 300.0:
    violations.append(f"MW > 300 ({mw:.1f})")
  if clogp > 3.0:
    violations.append(f"cLogP > 3.0 ({clogp:.2f})")
  if hbd > 3:
    violations.append(f"HBD > 3 ({hbd})")
  if hba > 3:
    violations.append(f"HBA > 3 ({hba})")
  if rotb > 3:
    violations.append(f"RotB > 3 ({rotb})")

  return len(violations) == 0, violations


def calculate_veber_rules(mol: Chem.Mol) -> Tuple[bool, List[str]]:
  """Evaluates Veber's Oral Bioavailability Rules (RotB <= 10, TPSA <= 140 A^2)."""
  if mol is None:
    return False, ["Invalid Molecule"]
  violations = []
  rotb = Lipinski.NumRotatableBonds(mol)
  tpsa = rdMolDescriptors.CalcTPSA(mol)

  if rotb > 10:
    violations.append(f"RotB > 10 ({rotb})")
  if tpsa > 140.0:
    violations.append(f"TPSA > 140 ({tpsa:.1f})")

  return len(violations) == 0, violations


def screen_structural_alerts(
    mol: Chem.Mol,
    pains: bool = True,
    brenk: bool = True,
    nih: bool = True,
) -> Dict[str, List[str]]:
  """Screens a molecule for PAINS, Brenk, and NIH structural alerts."""
  if mol is None:
    return {"pains": ["Invalid Molecule"], "brenk": [], "nih": []}

  results = {"pains": [], "brenk": [], "nih": []}

  if pains:
    p_params = FilterCatalogParams()
    p_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    p_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    p_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    cat = FilterCatalog(p_params)
    results["pains"] = [m.GetDescription() for m in cat.GetMatches(mol)]

  if brenk:
    b_params = FilterCatalogParams()
    b_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    cat = FilterCatalog(b_params)
    results["brenk"] = [m.GetDescription() for m in cat.GetMatches(mol)]

  if nih:
    n_params = FilterCatalogParams()
    n_params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
    cat = FilterCatalog(n_params)
    results["nih"] = [m.GetDescription() for m in cat.GetMatches(mol)]

  return results


def screen_chemical_dataset(
    df: pd.DataFrame,
    smiles_col: Optional[str] = None,
) -> pd.DataFrame:
  """Annotates a DataFrame with full MPO, Rule of 5/3, Veber, and PAINS/Brenk alert profiles."""
  df = df.copy()
  if not smiles_col:
    cands = [
        c for c in df.columns if "smi" in c.lower() or "structure" in c.lower()
    ]
    smiles_col = cands[0] if cands else df.columns[0]

  # Precompute RDKit Mols if not present
  mols = [
      Chem.MolFromSmiles(str(s)) if pd.notna(s) else None
      for s in df[smiles_col]
  ]

  cns_mpos = []
  ro5_viols = []
  ro3_passes = []
  veber_passes = []
  pains_alerts = []
  brenk_alerts = []
  nih_alerts = []

  pains_cat = get_filter_catalog(
      include_pains=True, include_brenk=False, include_nih=False
  )
  brenk_cat = get_filter_catalog(
      include_pains=False, include_brenk=True, include_nih=False
  )
  nih_cat = get_filter_catalog(
      include_pains=False, include_brenk=False, include_nih=True
  )

  for m in mols:
    if m is None:
      cns_mpos.append(0.0)
      ro5_viols.append(4)
      ro3_passes.append(False)
      veber_passes.append(False)
      pains_alerts.append("Invalid")
      brenk_alerts.append("Invalid")
      nih_alerts.append("Invalid")
      continue

    cns_mpos.append(calculate_cns_mpo(m))
    v_count, _ = calculate_rule_of_5(m)
    ro5_viols.append(v_count)
    r3_pass, _ = calculate_rule_of_3(m)
    ro3_passes.append(r3_pass)
    veber_pass, _ = calculate_veber_rules(m)
    veber_passes.append(veber_pass)

    p_matches = pains_cat.GetMatches(m)
    pains_alerts.append(
        ", ".join([x.GetDescription() for x in p_matches])
        if p_matches
        else "Clean"
    )

    b_matches = brenk_cat.GetMatches(m)
    brenk_alerts.append(
        ", ".join([x.GetDescription() for x in b_matches])
        if b_matches
        else "Clean"
    )

    n_matches = nih_cat.GetMatches(m)
    nih_alerts.append(
        ", ".join([x.GetDescription() for x in n_matches])
        if n_matches
        else "Clean"
    )

  df["CNS_MPO"] = cns_mpos
  df["Ro5_Violations"] = ro5_viols
  df["Ro3_Pass"] = ro3_passes
  df["Veber_Pass"] = veber_passes
  df["PAINS_Status"] = pains_alerts
  df["Brenk_Status"] = brenk_alerts
  df["NIH_Status"] = nih_alerts

  return df


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Screen molecules for PAINS/Brenk/NIH alerts, CNS MPO, and"
          " Drug-Likeness rules."
      )
  )
  parser.add_argument(
      "--smiles", type=str, default=None, help="Single SMILES string to screen."
  )
  parser.add_argument(
      "--input", type=str, default=None, help="Input CSV or SDF dataset file."
  )
  parser.add_argument(
      "--smiles_col",
      type=str,
      default=None,
      help="Name of SMILES column in CSV.",
  )
  parser.add_argument(
      "--output",
      type=str,
      default=None,
      help="Output CSV path for annotated dataset.",
  )
  parser.add_argument(
      "--json", action="store_true", help="Output results in JSON format."
  )
  args = parser.parse_args()

  if args.smiles:
    mol = Chem.MolFromSmiles(args.smiles)
    if mol is None:
      print(f"Error: Invalid SMILES '{args.smiles}'", file=sys.stderr)
      sys.exit(1)

    cns_mpo = calculate_cns_mpo(mol)
    ro5_v_count, ro5_viols = calculate_rule_of_5(mol)
    ro3_pass, ro3_viols = calculate_rule_of_3(mol)
    veber_pass, veber_viols = calculate_veber_rules(mol)
    alerts = screen_structural_alerts(mol)

    report = {
        "SMILES": args.smiles,
        "CNS_MPO_Score": cns_mpo,
        "Lipinski_Ro5_Violations": ro5_v_count,
        "Lipinski_Violations_Detail": ro5_viols,
        "Fragment_Ro3_Pass": ro3_pass,
        "Veber_Oral_Bioavailability_Pass": veber_pass,
        "PAINS_Alerts": alerts["pains"] if alerts["pains"] else "Clean",
        "Brenk_Alerts": alerts["brenk"] if alerts["brenk"] else "Clean",
        "NIH_Alerts": alerts["nih"] if alerts["nih"] else "Clean",
    }

    if args.json:
      print(json.dumps(report, indent=2))
    else:
      print("=" * 60)
      print("  MEDICINAL CHEMISTRY & STRUCTURAL ALERT REPORT")
      print("=" * 60)
      print(f"  SMILES:             {args.smiles}")
      print(
          f"  CNS MPO Score:      {cns_mpo} / 6.0"
          f" ({'High CNS Probability' if cns_mpo >= 4.0 else 'Moderate/Low CNS'})"
      )
      print(
          f"  Lipinski Ro5:       {ro5_v_count} violations"
          f" {'(PASS)' if ro5_v_count <= 1 else '(FAIL)'}"
      )
      if ro5_viols:
        print(f"    - Details:        {', '.join(ro5_viols)}")
      print(f"  Fragment Ro3:       {'PASS' if ro3_pass else 'FAIL'}")
      print(f"  Veber Oral Rule:    {'PASS' if veber_pass else 'FAIL'}")
      print(
          "  PAINS Status:      "
          f" {', '.join(alerts['pains']) if alerts['pains'] else 'CLEAN (No PAINS alerts)'}"
      )
      print(
          "  Brenk Status:      "
          f" {', '.join(alerts['brenk']) if alerts['brenk'] else 'CLEAN'}"
      )
      print(
          "  NIH Status:        "
          f" {', '.join(alerts['nih']) if alerts['nih'] else 'CLEAN'}"
      )
      print("=" * 60)

  elif args.input and os.path.exists(args.input):
    if args.input.endswith(".csv"):
      df = pd.read_csv(args.input)
    elif args.input.endswith(".sdf"):
      df = PandasTools.LoadSDF(args.input)
    else:
      print("Error: Unsupported format. Provide .csv or .sdf", file=sys.stderr)
      sys.exit(1)

    annotated_df = screen_chemical_dataset(df, smiles_col=args.smiles_col)

    pains_clean = (annotated_df["PAINS_Status"] == "Clean").sum()
    ro5_pass = (annotated_df["Ro5_Violations"] <= 1).sum()
    mean_cns = annotated_df["CNS_MPO"].mean()

    print("=" * 60)
    print(f"  DATASET FILTERING REPORT ({len(df)} compounds)")
    print("=" * 60)
    print(
        f"  Lipinski Ro5 Pass Rate:    {ro5_pass}/{len(df)}"
        f" ({ro5_pass/len(df)*100:.1f}%)"
    )
    print(
        f"  PAINS Clean Rate:          {pains_clean}/{len(df)}"
        f" ({pains_clean/len(df)*100:.1f}%)"
    )
    print(f"  Mean Pfizer CNS MPO:       {mean_cns:.2f} / 6.0")
    print("=" * 60)

    if args.output:
      annotated_df.to_csv(args.output, index=False)
      print(f"Saved annotated dataset to: {os.path.abspath(args.output)}")
  else:
    print("Error: Please provide --smiles or --input.", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
  main()
