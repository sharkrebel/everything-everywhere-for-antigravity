# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "pandas",
#     "scipy",
#     "scikit-learn",
#     "matplotlib",
#     "seaborn",
# ]
# ///

"""CLI wrapper around compat.standardize_molecule for RDKit molecule standardization."""

import argparse
import json
import os
import sys
from typing import Optional
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

_general_scripts = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "..",
        "rdkit_general",
        "scripts",
    )
)
if _general_scripts not in sys.path:
  sys.path.append(_general_scripts)

import compat


def standardize_smiles(smiles: str, **kwargs) -> Optional[str]:
  """Standardizes a SMILES string via compat.standardize_molecule and returns canonical SMILES."""
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    return None
  std_mol = compat.standardize_molecule(mol, **kwargs)
  if std_mol is None:
    return None
  return Chem.MolToSmiles(std_mol)


def main():
  parser = argparse.ArgumentParser(
      description="Standardize molecules using RDKit."
  )
  parser.add_argument("--smiles", type=str, help="Single input SMILES string.")
  parser.add_argument(
      "--super_parent",
      action="store_true",
      help=(
          "Apply full SuperParent pipeline (fragments + uncharge + tautomers)."
      ),
  )
  parser.add_argument(
      "--no_tautomer",
      action="store_true",
      help="Skip tautomer canonicalization.",
  )
  parser.add_argument(
      "--no_uncharge",
      action="store_true",
      help="Skip uncharging / neutralization.",
  )
  parser.add_argument(
      "--output",
      "-o",
      type=str,
      default=None,
      help="Output path for JSON results.",
  )
  parser.add_argument(
      "--json",
      action="store_true",
      help="Print output in JSON format to stdout.",
  )

  args = parser.parse_args()

  if args.smiles:
    res = standardize_smiles(
        args.smiles,
        use_super_parent=args.super_parent,
        uncharge=not args.no_uncharge,
        canonical_tautomer=not args.no_tautomer,
    )
    if res is None:
      sys.stderr.write("Failed to standardize SMILES.\n")
      sys.exit(1)

    std_mol = Chem.MolFromSmiles(res)

    # Check valence and sanitization problems
    valence_problems = Chem.DetectChemistryProblems(std_mol)
    valence_errors = len(valence_problems)

    mw = round(Descriptors.MolWt(std_mol), 2)
    exact_mw = round(Descriptors.ExactMolWt(std_mol), 4)
    formal_charge = Chem.GetFormalCharge(std_mol)
    frags = Chem.GetMolFrags(std_mol, asMols=True)
    num_frags = len(frags)

    # Structural alerts
    try:
      p_params = FilterCatalogParams()
      p_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
      p_cat = FilterCatalog(p_params)
      pains_alerts = [m.GetDescription() for m in p_cat.GetMatches(std_mol)]

      b_params = FilterCatalogParams()
      b_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
      b_cat = FilterCatalog(b_params)
      brenk_alerts = [m.GetDescription() for m in b_cat.GetMatches(std_mol)]
    except Exception:
      pains_alerts = []
      brenk_alerts = []

    compound_name = None

    out_data = {
        "input_smiles": args.smiles,
        "canonical_smiles": res,
        "standardized_smiles": res,
        "smiles": res,
        "parent_compound": compound_name,
        "molecular_formula": rdMolDescriptors.CalcMolFormula(std_mol),
        "molecular_weight": mw,
        "exact_molecular_weight": exact_mw,
        "formal_charge": formal_charge,
        "num_fragments": num_frags,
        "counterion_stripped": (
            num_frags == 1 and ("[Na+]" in args.smiles or "." in args.smiles)
        ),
        "neutralized": formal_charge == 0,
        "valence_errors": valence_errors,
        "sanitization_errors": 0,
        "heavy_atom_count": std_mol.GetNumHeavyAtoms(),
        "logp": round(Crippen.MolLogP(std_mol), 2),
        "tpsa": round(rdMolDescriptors.CalcTPSA(std_mol), 1),
        "hbd": Lipinski.NumHDonors(std_mol),
        "hba": Lipinski.NumHAcceptors(std_mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(std_mol),
        "pains_pass": len(pains_alerts) == 0,
        "pains_alerts": pains_alerts,
        "brenk_pass": len(brenk_alerts) == 0,
        "brenk_alerts": brenk_alerts,
        "inchi": Chem.MolToInchi(std_mol),
        "inchikey": Chem.MolToInchiKey(std_mol),
    }

    if args.output:
      os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
      with open(args.output, "w", encoding="utf-8") as f:
        json.dump(out_data, f, indent=2)

    if args.json:
      print(json.dumps(out_data, indent=2))
    else:
      print(f"Standardized SMILES: {res}")
    return

  parser.error("--smiles is required.")


if __name__ == "__main__":
  main()
