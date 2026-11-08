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

"""Arbitrary Reaction SMARTS Construction, Validation, Transformation, and Search Template Extraction Engine."""

import argparse
import json
import os
import sys
from typing import Any, Dict, List
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors


def validate_and_parse_reaction_smarts(reaction_smarts: str) -> Dict[str, Any]:
  """Validates an arbitrary reaction SMARTS, checks atom mapping balance, and extracts search templates."""
  rxn = AllChem.ReactionFromSmarts(reaction_smarts)
  if rxn is None:
    return {
        "valid": False,
        "error": f"Failed to parse Reaction SMARTS: {reaction_smarts}",
    }

  n_warn, n_err = rxn.Validate()
  num_reactants = rxn.GetNumReactantTemplates()
  num_products = rxn.GetNumProductTemplates()

  # Extract clean unmapped query SMARTS for each reactant template for database/supplier searching
  reactant_templates = []
  for i in range(num_reactants):
    tmpl = rxn.GetReactantTemplate(i)
    clean_tmpl = Chem.Mol(tmpl)
    for atom in clean_tmpl.GetAtoms():
      atom.SetAtomMapNum(0)
    reactant_templates.append({
        "index": i + 1,
        "mapped_smarts": Chem.MolToSmarts(tmpl),
        "search_query_smarts": Chem.MolToSmarts(clean_tmpl),
    })

  return {
      "valid": n_err == 0,
      "reaction_smarts": reaction_smarts,
      "num_reactants": num_reactants,
      "num_products": num_products,
      "validation_warnings": n_warn,
      "validation_errors": n_err,
      "reactant_templates": reactant_templates,
  }


def run_reaction_smarts(
    reaction_smarts: str,
    reactant_smiles: List[str],
    sanitize_products: bool = True,
) -> List[List[str]]:
  """Executes an arbitrary reaction SMARTS transformation on input reactant SMILES strings."""
  rxn = AllChem.ReactionFromSmarts(reaction_smarts)
  if rxn is None:
    raise ValueError(f"Invalid Reaction SMARTS: {reaction_smarts}")

  if rxn.GetNumReactantTemplates() != len(reactant_smiles):
    raise ValueError(
        f"Reaction expects {rxn.GetNumReactantTemplates()} reactant(s), "
        f"but {len(reactant_smiles)} were provided."
    )

  reactant_mols = []
  for i, smi in enumerate(reactant_smiles):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
      raise ValueError(f"Invalid reactant SMILES at position {i}: {smi}")
    reactant_mols.append(mol)

  raw_products = rxn.RunReactants(tuple(reactant_mols))
  unique_product_sets = []
  seen_smiles = set()

  for prod_tuple in raw_products:
    prod_smiles_list = []
    valid_set = True
    for prod_mol in prod_tuple:
      if sanitize_products:
        try:
          Chem.SanitizeMol(prod_mol)
        except Exception as e:
          sys.stderr.write(f"Product sanitization failed: {e}\n")
          valid_set = False
          break
      prod_smiles = Chem.MolToSmiles(prod_mol)
      prod_smiles_list.append(prod_smiles)

    if valid_set:
      key = tuple(prod_smiles_list)
      if key not in seen_smiles:
        seen_smiles.add(key)
        unique_product_sets.append(prod_smiles_list)

  return unique_product_sets


REACTIVE_HANDLES = [
    {
        "name": "Boronic acid pinacol ester (Bpin)",
        "smarts": "[#5]1-O-C(C)(C)-C(C)(C)-O-1",
        "reactions": [
            {
                "reaction_name": "Suzuki-Miyaura Cross-Coupling (Recommended)",
                "description": (
                    "Forms C-C biaryl/heteroaryl bonds with aryl/heteroaryl"
                    " halides. High chemoselectivity; tolerates free amino"
                    " groups."
                ),
                "reaction_smarts": (
                    "[c:1][B]1[O][C]([C])([C])[C]([C])([C])[O]1.[c:2][Br,I,Cl]>>[c:1][c:2]"
                ),
                "partner_smarts": "c[Br]",
                "partner_type": "Aryl/Heteroaryl Halides (c[Br], c[I], c[Cl])",
            },
            {
                "reaction_name": "Chan-Lam C-N / C-O Coupling",
                "description": (
                    "Copper-catalyzed coupling with amines, phenols, or amides"
                    " at room temperature."
                ),
                "reaction_smarts": (
                    "[c:1][B]1[O][C]([C])([C])[C]([C])([C])[O]1.[N;H2,H1:2]>>[c:1][N:2]"
                ),
                "partner_smarts": "[N;H2,H1]",
                "partner_type": "Amines / Heterocycles",
            },
        ],
    },
    {
        "name": "Boronic acid",
        "smarts": "[#5](-[OH])-[OH]",
        "reactions": [{
            "reaction_name": "Suzuki-Miyaura Cross-Coupling",
            "description": "Couples boronic acids with aryl halides.",
            "reaction_smarts": "[c:1][B]([OH])[OH].[c:2][Br,I,Cl]>>[c:1][c:2]",
            "partner_smarts": "c[Br]",
            "partner_type": "Aryl Halides",
        }],
    },
    {
        "name": "Primary aromatic amine",
        "smarts": "c[NH2]",
        "reactions": [
            {
                "reaction_name": "Buchwald-Hartwig Amination",
                "description": (
                    "Palladium-catalyzed C-N bond formation with aryl halides."
                ),
                "reaction_smarts": "[c:1][NH2].[c:2][Br,I,Cl]>>[c:1][NH][c:2]",
                "partner_smarts": "c[Br]",
                "partner_type": "Aryl Halides",
            },
            {
                "reaction_name": "Amide Coupling / Acylation",
                "description": "Forms aromatic carboxamides or sulfonamides.",
                "reaction_smarts": (
                    "[c:1][NH2].[C:2](=[O:3])[OH]>>[c:1][NH][C:2](=[O:3])"
                ),
                "partner_smarts": "C(=O)O",
                "partner_type": "Carboxylic Acids",
            },
        ],
    },
    {
        "name": "Carboxylic acid",
        "smarts": "C(=O)[OH]",
        "reactions": [{
            "reaction_name": "Amide Bond Formation",
            "description": (
                "Coupling with primary or secondary amines via HATU/EDC/HOBt."
            ),
            "reaction_smarts": (
                "[C:1](=[O:2])[OH].[N;H2,H1:3]>>[C:1](=[O:2])[N:3]"
            ),
            "partner_smarts": "[N;H2,H1]",
            "partner_type": "Amines",
        }],
    },
    {
        "name": "Aryl halide (Br/I)",
        "smarts": "c[Br,I]",
        "reactions": [{
            "reaction_name": "Suzuki-Miyaura Cross-Coupling",
            "description": (
                "Coupling with boronic acids or boronate pinacol esters."
            ),
            "reaction_smarts": "[c:1][Br,I].[c:2][B]>>[c:1][c:2]",
            "partner_smarts": "B(O)O",
            "partner_type": "Boronic Acids / Esters",
        }],
    },
]


def analyze_core_molecule(smiles: str) -> Dict[str, Any]:
  """Analyzes a starting core molecule, perceives reactive functional handles, and suggests compatible reactions."""
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    raise ValueError(f"Invalid SMILES string: {smiles}")

  can_smi = Chem.MolToSmiles(mol)
  formula = rdMolDescriptors.CalcMolFormula(mol)
  mw = round(Descriptors.ExactMolWt(mol), 2)
  clogp = round(Crippen.MolLogP(mol), 2)
  tpsa = round(rdMolDescriptors.CalcTPSA(mol), 1)

  detected_handles = []
  for h in REACTIVE_HANDLES:
    p = Chem.MolFromSmarts(h["smarts"])
    if p and mol.HasSubstructMatch(p):
      matches = mol.GetSubstructMatches(p)
      detected_handles.append({
          "handle_name": h["name"],
          "matches_count": len(matches),
          "reactions": h["reactions"],
      })

  return {
      "smiles": can_smi,
      "formula": formula,
      "mw": mw,
      "clogp": clogp,
      "tpsa": tpsa,
      "handles": detected_handles,
  }


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Construct, validate, analyze, and execute Reaction SMARTS"
          " transformations."
      )
  )
  parser.add_argument(
      "--analyze",
      type=str,
      default=None,
      help="Analyze starting material SMILES and suggest compatible reactions.",
  )
  parser.add_argument(
      "--smarts", type=str, default=None, help="Reaction SMARTS string."
  )
  parser.add_argument(
      "--reactants", nargs="*", default=None, help="Reactant SMILES strings."
  )
  parser.add_argument(
      "--validate",
      action="store_true",
      help="Validate SMARTS and extract unmapped search queries.",
  )
  parser.add_argument(
      "--output",
      "-o",
      type=str,
      default=None,
      help="Output file path to save reaction analysis or products as JSON.",
  )
  args = parser.parse_args()

  if args.analyze:
    analysis = analyze_core_molecule(args.analyze)
    if args.output:
      os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
      with open(args.output, "w", encoding="utf-8") as f:
        json.dump(analysis, f, indent=2)
    print("=" * 80)
    print(f"  CORE MOLECULE ANALYSIS: {analysis['smiles']}")
    print(
        f"  Formula: {analysis['formula']} | MW: {analysis['mw']} | cLogP:"
        f" {analysis['clogp']} | TPSA: {analysis['tpsa']}"
    )
    print("=" * 80)
    if not analysis["handles"]:
      print("No common reactive functional handles detected.")
    else:
      print(
          f"\nDetected {len(analysis['handles'])} reactive functional"
          " handle(s):"
      )
      for idx, h in enumerate(analysis["handles"], 1):
        print(
            f"\n[{idx}] Handle: {h['handle_name']} ({h['matches_count']}"
            " site(s))"
        )
        for rxn in h["reactions"]:
          print(f"    * Reaction: {rxn['reaction_name']}")
          print(f"      Description: {rxn['description']}")
          print(f"      Reaction SMARTS: {rxn['reaction_smarts']}")
          print(f"      Compatible Partners: {rxn['partner_type']}")
          print("      Discovery Command:")
          print("        python3 scripts/find_coupling_partners.py \\")
          print(f"          --core \"{analysis['smiles']}\" \\")
          print(f"          --reaction_smarts \"{rxn['reaction_smarts']}\" \\")
          print(f"          --partner_smarts \"{rxn['partner_smarts']}\" \\")
          print("          --all --output_csv scratch/derivatives.csv")
    return

  if not args.smarts:
    parser.error("Either --analyze or --smarts must be specified.")

  info = validate_and_parse_reaction_smarts(args.smarts)
  if not info["valid"]:
    print(
        f"Error: {info.get('error', 'Reaction SMARTS validation failed')}",
        file=sys.stderr,
    )
    sys.exit(1)

  if args.validate or not args.reactants:
    print(f"Reaction SMARTS: {args.smarts}")
    print(
        f"  Status: {'VALID' if info['valid'] else 'INVALID'}"
        f" ({info['validation_warnings']} warnings, {info['validation_errors']}"
        " errors)"
    )
    print(
        f"  Expects: {info['num_reactants']} reactant(s) -> yields"
        f" {info['num_products']} product(s)"
    )
    print("\n  Reactant Search Templates (Use to query suppliers / databases):")
    for tmpl in info.get("reactant_templates", []):
      print(
          f"    Reactant {tmpl['index']} Query SMARTS:"
          f" {tmpl['search_query_smarts']}"
      )

  products = []
  if args.reactants:
    products = run_reaction_smarts(args.smarts, args.reactants)
    print(
        f"\nExecution on reactants {args.reactants} yielded {len(products)}"
        " distinct product set(s):"
    )
    for i, prod_set in enumerate(products):
      print(f"  Outcome {i+1}: {' + '.join(prod_set)}")

  if args.output:
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
      json.dump({"validation": info, "products": products}, f, indent=2)


if __name__ == "__main__":
  main()
