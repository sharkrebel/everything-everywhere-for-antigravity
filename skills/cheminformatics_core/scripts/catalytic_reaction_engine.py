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

"""Organometallic Catalysis Engine & Privileged Ligand Screening Suite.

Provides curated, verified Reaction SMARTS for named transition-metal catalyzed
cross-coupling reactions, C-H functionalization, and olefin metathesis:
1. Suzuki-Miyaura (Boronic Acid + Aryl Halide -> Biaryl)
2. Buchwald-Hartwig (Amine + Aryl Halide -> Arylamine)
3. Heck Reaction (Alkene + Aryl Halide -> Substituted Alkene)
4. Sonogashira (Terminal Alkyne + Aryl Halide -> Internal Alkyne)
5. Negishi Coupling (Organozinc + Aryl Halide -> Coupled Product)
6. Stille Coupling (Organotin + Aryl Halide -> Coupled Product)
7. C-H Borylation (Arene + B2pin2 -> Arylboronate)
8. Olefin Metathesis (Ring-Closing & Cross Metathesis)

Includes a catalog of privileged phosphine, bisphosphine, and NHC ligands with
precomputed steric and electronic descriptors to accelerate catalyst selection.
"""

import argparse
import json
from typing import Any, Dict, List, Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Curated Transition-Metal Catalyzed Reaction SMARTS
NAMED_CATALYTIC_REACTIONS: Dict[str, Dict[str, str]] = {
    "suzuki_miyaura": {
        "name": "Suzuki-Miyaura Cross-Coupling",
        "description": (
            "Pd-catalyzed coupling of boronic acids/esters with aryl/heteroaryl"
            " halides/triflates."
        ),
        "smarts": (
            "[#6:1][B;$(B([OH,O])[OH,O])].[#6:2][Cl,Br,I,$(OS(=O)(=O)[#6])]>>[#6:1]-[#6:2]"
        ),
        "catalyst_class": "Pd(0) / Phosphines (SPhos, XPhos, RuPhos, PPh3)",
    },
    "buchwald_hartwig": {
        "name": "Buchwald-Hartwig Amination",
        "description": (
            "Pd-catalyzed C-N bond formation between aryl halides and"
            " primary/secondary amines."
        ),
        "smarts": (
            "[#6:1][Cl,Br,I,$(OS(=O)(=O)[#6])].[N;H2,H1;!$(N[C,S]=O):2][#6:3]>>[#6:1]-[N:2]-[#6:3]"
        ),
        "catalyst_class": (
            "Pd(0) / Dialkylbiaryl Phosphines (BrettPhos, RuPhos, BINAP)"
        ),
    },
    "heck": {
        "name": "Mizoroki-Heck Alkene Coupling",
        "description": (
            "Pd-catalyzed coupling of aryl halides with activated alkenes."
        ),
        "smarts": "[#6:1][Cl,Br,I].[C:2]=[C;H2,H1:3]>>[#6:1]-[C:2]=[C:3]",
        "catalyst_class": "Pd(0) / Bulky Phosphines (P(o-tol)3, P(t-Bu)3)",
    },
    "sonogashira": {
        "name": "Sonogashira Alkyne Coupling",
        "description": (
            "Pd/Cu co-catalyzed coupling of terminal alkynes with aryl/alkenyl"
            " halides."
        ),
        "smarts": "[#6:1][Cl,Br,I].[C;H1:2]#[C:3]>>[#6:1]-[C:2]#[C:3]",
        "catalyst_class": "PdCl2(PPh3)2 / CuI",
    },
    "negishi": {
        "name": "Negishi Cross-Coupling",
        "description": (
            "Pd/Ni-catalyzed coupling of organozinc halides with organohalides."
        ),
        "smarts": "[#6:1][Zn][Cl,Br,I].[#6:2][Cl,Br,I]>>[#6:1]-[#6:2]",
        "catalyst_class": "Pd(0) / PEPPSI (NHC) or Ni(II)",
    },
    "stille": {
        "name": "Stille Cross-Coupling",
        "description": (
            "Pd-catalyzed coupling of organotin reagents with aryl/alkenyl"
            " halides."
        ),
        "smarts": "[#6:1][Sn]([#6])([#6])[#6].[#6:2][Cl,Br,I]>>[#6:1]-[#6:2]",
        "catalyst_class": "Pd(PPh3)4 / Pd2(dba)3",
    },
    "ch_borylation": {
        "name": "Iridium-Catalyzed C-H Borylation",
        "description": (
            "Ir-catalyzed direct C-H borylation of arenes and heteroarenes with"
            " B2pin2."
        ),
        "smarts": (
            "[c;H1:1].[B]([O:2])([O:3])[B]([O:4])[O:5]>>[c:1][B]([O:2])[O:3]"
        ),
        "catalyst_class": (
            "[Ir(cod)(OMe)]2 / dtbpy (4,4'-di-tert-butylbipyridine)"
        ),
    },
    "olefin_metathesis": {
        "name": "Olefin Metathesis (Cross & Ring-Closing)",
        "description": (
            "Ru-catalyzed alkene redistribution via ruthenacyclobutane"
            " intermediates."
        ),
        "smarts": "[C:1]=[C;H2,H1:2].[C:3]=[C;H2,H1:4]>>[C:1]=[C:3]",
        "catalyst_class": "Grubbs 2nd Gen / Hoveyda-Grubbs 2nd Gen",
    },
}

# Catalog of Privileged Catalytic Ligands with Precomputed Descriptors
PRIVILEGED_LIGANDS: List[Dict[str, Any]] = [
    {
        "name": "SPhos",
        "cas": "657408-07-6",
        "smiles": "COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1",
        "class": "Dialkylbiaryl Phosphine",
        "percent_v_bur": 36.8,
        "cone_angle_deg": 195.0,
        "bite_angle_deg": None,
        "application": (
            "Sterically hindered Suzuki-Miyaura coupling & aryl chlorides"
        ),
    },
    {
        "name": "XPhos",
        "cas": "564483-18-7",
        "smiles": "CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2CCCCC2)c(C(C)C)c1",
        "class": "Dialkylbiaryl Phosphine",
        "percent_v_bur": 38.2,
        "cone_angle_deg": 204.0,
        "bite_angle_deg": None,
        "application": (
            "General Suzuki, Buchwald amination & challenging cross-couplings"
        ),
    },
    {
        "name": "RuPhos",
        "cas": "787618-22-8",
        "smiles": "CC(C)Oc1cccc(OC(C)C)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1",
        "class": "Dialkylbiaryl Phosphine",
        "percent_v_bur": 37.4,
        "cone_angle_deg": 198.0,
        "bite_angle_deg": None,
        "application": "Buchwald-Hartwig amination of secondary amines",
    },
    {
        "name": "BrettPhos",
        "cas": "1070663-78-3",
        "smiles": "COc1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C(C)(C)C)c(C(C)C)c1OC",
        "class": "Dialkylbiaryl Phosphine",
        "percent_v_bur": 44.5,
        "cone_angle_deg": 221.0,
        "bite_angle_deg": None,
        "application": "Primary amine amination & C-O/C-F functionalization",
    },
    {
        "name": "BINAP",
        "cas": "98327-87-8",
        "smiles": (
            "c1ccc(P(c2ccccc2)c2ccc3ccccc3c2-c2c(P(c3ccccc3)c3ccccc3)ccc3ccccc23)cc1"
        ),
        "class": "Chiral Bisphosphine",
        "percent_v_bur": 42.1,
        "cone_angle_deg": 218.0,
        "bite_angle_deg": 92.7,
        "application": (
            "Asymmetric hydrogenation & enantioselective cross-coupling"
        ),
    },
    {
        "name": "Xantphos",
        "cas": "161265-03-8",
        "smiles": (
            "CC1(C)c2cccc(P(c3ccccc3)c3ccccc3)c2Oc2c(P(c3ccccc3)c3ccccc3)cccc21"
        ),
        "class": "Wide-Bite-Angle Bisphosphine",
        "percent_v_bur": 46.0,
        "cone_angle_deg": 228.0,
        "bite_angle_deg": 111.4,
        "application": "Hydroformylation, carbonylation & borylation",
    },
    {
        "name": "IPr",
        "cas": "244187-81-3",
        "smiles": "CC(C)c1cccc(C(C)C)c1N1C=CN(c2c(C(C)C)cccc2C(C)C)C1=[C]",
        "class": "N-Heterocyclic Carbene (NHC)",
        "percent_v_bur": 51.2,
        "cone_angle_deg": 235.0,
        "bite_angle_deg": None,
        "application": "Buchwald, Negishi & metathesis catalyst stabilization",
    },
    {
        "name": "dtbpy",
        "cas": "72914-19-3",
        "smiles": "CC(C)(C)c1ccn(-c2cc(C(C)(C)C)ccn2)cc1",
        "class": "Bidentate Bipyridine",
        "percent_v_bur": 28.5,
        "cone_angle_deg": None,
        "bite_angle_deg": 75.0,
        "application": "Ir-catalyzed C-H borylation & photoredox catalysis",
    },
]


def run_catalytic_reaction(
    reaction_name: str,
    reactants: List[str],
) -> List[str]:
  """Executes a named organometallic catalytic transformation on input reactant SMILES."""
  rxn_key = reaction_name.lower().replace("-", "_").replace(" ", "_")
  if rxn_key not in NAMED_CATALYTIC_REACTIONS:
    raise ValueError(
        f"Unknown reaction '{reaction_name}'. Supported:"
        f" {list(NAMED_CATALYTIC_REACTIONS.keys())}"
    )

  rxn_info = NAMED_CATALYTIC_REACTIONS[rxn_key]
  rxn_smarts = rxn_info["smarts"]
  rxn = AllChem.ReactionFromSmarts(rxn_smarts)

  mol_reactants = [Chem.MolFromSmiles(r) for r in reactants]
  if any(m is None for m in mol_reactants):
    raise ValueError("Could not parse one or more reactant SMILES.")

  products_tuple = rxn.RunReactants(tuple(mol_reactants))
  unique_products = set()

  for prod_set in products_tuple:
    for p in prod_set:
      try:
        Chem.SanitizeMol(p)
        unique_products.add(Chem.MolToSmiles(p))
      except Exception:
        pass

  return sorted(list(unique_products))


def list_privileged_ligands(
    *,
    min_vbur: Optional[float] = None,
    max_vbur: Optional[float] = None,
    ligand_class: Optional[str] = None,
) -> pd.DataFrame:
  """Returns a filtered DataFrame of privileged catalytic ligands."""
  df = pd.DataFrame(PRIVILEGED_LIGANDS)
  if min_vbur is not None:
    df = df[df["percent_v_bur"] >= min_vbur]
  if max_vbur is not None:
    df = df[df["percent_v_bur"] <= max_vbur]
  if ligand_class is not None:
    df = df[df["class"].str.contains(ligand_class, case=False)]
  return df


def main():
  parser = argparse.ArgumentParser(
      description="Organometallic Catalysis & Ligand Screening Suite."
  )
  parser.add_argument(
      "--reaction",
      type=str,
      choices=list(NAMED_CATALYTIC_REACTIONS.keys()),
      help="Named catalytic reaction.",
  )
  parser.add_argument(
      "--reactants", type=str, nargs="+", help="Reactant SMILES strings."
  )
  parser.add_argument(
      "--list_ligands",
      action="store_true",
      help="Print privileged ligand catalog.",
  )
  parser.add_argument(
      "--filter_vbur",
      type=float,
      nargs=2,
      metavar=("MIN", "MAX"),
      help="Filter ligands by %%V_bur range.",
  )
  parser.add_argument("--json", action="store_true", help="Output JSON format.")
  parser.add_argument(
      "--output",
      "--output-json",
      dest="output_path",
      type=str,
      default=None,
      help="Save JSON output to a file.",
  )
  args = parser.parse_args()

  if args.list_ligands or args.filter_vbur:
    min_v = args.filter_vbur[0] if args.filter_vbur else None
    max_v = args.filter_vbur[1] if args.filter_vbur else None
    df = list_privileged_ligands(min_vbur=min_v, max_vbur=max_v)
    if args.output_path:
      df.to_json(args.output_path, orient="records", indent=2)
      print(f"[+] Saved privileged ligand catalog JSON to: {args.output_path}")
    if args.json:
      print(df.to_json(orient="records", indent=2))
    else:
      print("\n=== Privileged Catalytic Ligand Catalog ===")
      print(
          df[[
              "name",
              "cas",
              "class",
              "percent_v_bur",
              "cone_angle_deg",
              "bite_angle_deg",
              "application",
          ]].to_string(index=False)
      )
    return

  if args.reaction and args.reactants:
    prods = run_catalytic_reaction(args.reaction, args.reactants)
    rxn_info = NAMED_CATALYTIC_REACTIONS[args.reaction]
    payload = {
        "reaction": rxn_info["name"],
        "catalyst_class": rxn_info["catalyst_class"],
        "reactants": args.reactants,
        "products": prods,
    }
    if args.output_path:
      with open(args.output_path, "w") as f:
        json.dump(payload, f, indent=2)
      print(f"[+] Saved catalytic reaction JSON to: {args.output_path}")
    print(f"\n=== Catalytic Transformation: {rxn_info['name']} ===")
    print(f"Catalyst System: {rxn_info['catalyst_class']}")
    print(f"Reactants:       {' + '.join(args.reactants)}")
    print(f"Products ({len(prods)}):")
    for i, p in enumerate(prods, 1):
      print(f"  {i}. {p}")
    return

  parser.error("Must provide --catalyst_smiles or --reaction_smarts.")


if __name__ == "__main__":
  main()
