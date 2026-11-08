# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "pandas",
#     "scipy",
#     "scikit-learn",
#     "matplotlib",
#     "polite-http",
# ]
# ///

"""Retrosynthetic Analysis, Strategic 1-Step Disconnections, RECAP Tree, and BRICS Synthon Sourcing."""

import argparse
import concurrent.futures
import os
from typing import Any, Dict, List
import pandas as pd
from polite_http import http_client
from rdkit import Chem
from rdkit.Chem import AllChem, BRICS, Descriptors, Recap

_PUBCHEM_CLIENT = http_client.HttpClient(
    "https://pubchem.ncbi.nlm.nih.gov/rest/", qps=5.0
)

COMMERCIAL_VENDORS = [
    "Enamine",
    "Chemspace",
    "MolPort",
    "Sigma-Aldrich",
    "WuXi AppTec",
    "Selleck Chemicals",
    "MedChemExpress (MCE)",
    "Tocris Bioscience",
    "ChemDiv",
    "Life Chemicals",
    "Specs",
    "TargetMol",
    "BOC Sciences",
    "Combi-Blocks",
    "Fluorochem",
    "Oakwood Chemical",
    "Thermo Fisher Scientific",
    "Alfa Aesar",
    "Ambeed",
    "BLDpharm",
]

STRATEGIC_DISCONNECTIONS = [
    (
        "Amide Coupling",
        "[C:1](=[O:2])!@[N;!$(N=O):3]>>[C:1](=[O:2])O.[N:3]",
        "Carboxylic Acid + Amine",
    ),
    (
        "Suzuki-Miyaura Biaryl",
        "[c:1]!@[c:2]>>[c:1]B(O)O.[c:2]Br",
        "Aryl Boronic Acid + Aryl Bromide",
    ),
    (
        "Buchwald-Hartwig C-N",
        "[c:1]!@[N;!$(NC=O);!$(N=O):2]>>[c:1]Br.[N:2]",
        "Aryl Bromide + Amine",
    ),
    (
        "Sulfonamide Formation",
        "[S:1](=[O:2])(=[O:3])!@[N:4]>>[S:1](=[O:2])(=[O:3])Cl.[N:4]",
        "Sulfonyl Chloride + Amine",
    ),
    (
        "Esterification",
        "[C:1](=[O:2])!@[O:3][C:4]>>[C:1](=[O:2])O.[O:3][C:4]",
        "Carboxylic Acid + Alcohol",
    ),
    (
        "Ether / SNAr Alkylation",
        "[c:1]!@[O:2][C:3]>>[c:1][O:2].[C:3]Br",
        "Phenol + Alkyl Bromide",
    ),
    (
        "Reductive Amination",
        "[c:1][C;H2:2]!@[N;!$(NC=O);!$(N=O):3]>>[c:1][C:2]=O.[N:3]",
        "Aldehyde + Amine",
    ),
    (
        "Urea Formation",
        "[N:1]!@[C:2](=[O:3])!@[N:4]>>[N:1][C:2](=[O:3])Cl.[N:4]",
        "Carbamoyl Chloride + Amine",
    ),
]


def is_valid_cas(cas_str: str) -> bool:
  """Validates format and Mod-10 checksum of a CAS Registry Number."""
  if not isinstance(cas_str, str):
    return False
  parts = cas_str.strip().split("-")
  if len(parts) != 3:
    return False
  p1, p2, p3 = parts
  if not (
      p1.isdigit()
      and p2.isdigit()
      and p3.isdigit()
      and 2 <= len(p1) <= 7
      and len(p2) == 2
      and len(p3) == 1
  ):
    return False
  digits = p1 + p2
  return sum(
      int(d) * (len(digits) - i) for i, d in enumerate(digits)
  ) % 10 == int(p3)


def batch_resolve_precursor_sourcing(
    smiles_list: List[str],
) -> Dict[str, Dict[str, Any]]:
  """Parallel batch resolution of precursor commercial availability, CAS numbers, and vendors."""
  results = {
      smi: {"cas": "N/A", "vendors": [], "available": False, "cid": None}
      for smi in smiles_list
  }
  unique_smiles = [s for s in set(smiles_list) if s]

  def resolve_single_smiles(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
      return smi, None, None, []

    inchikey = Chem.MolToInchiKey(mol)
    cid = None
    try:
      data = _PUBCHEM_CLIENT.fetch_json(
          f"pug/compound/inchikey/{inchikey}/property/CanonicalSMILES/JSON",
          timeout=8.0,
      )
      props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
      cid = props.get("CID")
    except Exception:
      return smi, None, None, []

    if not cid:
      return smi, None, None, []

    # Fetch Synonyms (CAS) and Vendors (XRefs)
    cas_num = None
    vendors = []

    try:
      data_syns = _PUBCHEM_CLIENT.fetch_json(
          f"pug/compound/cid/{cid}/synonyms/JSON",
          timeout=8.0,
      )
      syns = (
          data_syns.get("InformationList", {})
          .get("Information", [{}])[0]
          .get("Synonym", [])
      )
      valid_cas = [s for s in syns if is_valid_cas(s)]
      if valid_cas:
        cas_num = valid_cas[0]
    except Exception:
      pass

    try:
      data_xrefs = _PUBCHEM_CLIENT.fetch_json(
          f"pug/compound/cid/{cid}/xrefs/SourceName/JSON",
          timeout=8.0,
      )
      sources = (
          data_xrefs.get("InformationList", {})
          .get("Information", [{}])[0]
          .get("SourceName", [])
      )
      matched = set()
      for s in sources:
        for v in COMMERCIAL_VENDORS:
          if v.lower() in s.lower():
            matched.add(v)
      vendors = sorted(list(matched))
    except Exception:
      pass

    return smi, cid, cas_num, vendors

  with concurrent.futures.ThreadPoolExecutor(
      max_workers=min(12, len(unique_smiles) or 1)
  ) as executor:
    for smi, cid, cas_num, vendors in executor.map(
        resolve_single_smiles, unique_smiles
    ):
      if cid:
        results[smi] = {
            "cid": cid,
            "cas": cas_num or "N/A",
            "vendors": vendors,
            "available": len(vendors) > 0,
        }

  return results


def disconnect_target_molecule(
    target_smiles: str,
    check_availability: bool = True,
) -> pd.DataFrame:
  """Performs strategic 1-step retrosynthetic disconnections and resolves commercial precursor availability."""
  mol = Chem.MolFromSmiles(target_smiles)
  if mol is None:
    raise ValueError(f"Invalid target SMILES: {target_smiles}")

  disconnection_records = []
  all_precursors = set()

  for rxn_name, rxn_smarts, precursor_types in STRATEGIC_DISCONNECTIONS:
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    outcomes = rxn.RunReactants((mol,))
    if not outcomes:
      continue

    seen_pairs = set()
    for out in outcomes:
      if len(out) < 2:
        continue
      mol_a, mol_b = out[0], out[1]
      try:
        Chem.SanitizeMol(mol_a)
        Chem.SanitizeMol(mol_b)
        smi_a = Chem.MolToSmiles(mol_a)
        smi_b = Chem.MolToSmiles(mol_b)
      except Exception:
        continue

      pair_key = tuple(sorted([smi_a, smi_b]))
      if pair_key in seen_pairs:
        continue
      seen_pairs.add(pair_key)

      all_precursors.add(smi_a)
      all_precursors.add(smi_b)

      disconnection_records.append({
          "Disconnection_Reaction": rxn_name,
          "Precursor_Classes": precursor_types,
          "Precursor_A_SMILES": smi_a,
          "Precursor_A_MW": round(Descriptors.ExactMolWt(mol_a), 1),
          "Precursor_B_SMILES": smi_b,
          "Precursor_B_MW": round(Descriptors.ExactMolWt(mol_b), 1),
      })

  if not disconnection_records:
    return pd.DataFrame()

  df = pd.DataFrame(disconnection_records)

  if check_availability:
    print(
        f"Resolving commercial availability for {len(all_precursors)} unique"
        " precursor building blocks..."
    )
    sourcing_info = batch_resolve_precursor_sourcing(list(all_precursors))

    df["Precursor_A_CAS"] = df["Precursor_A_SMILES"].apply(
        lambda s: sourcing_info.get(s, {}).get("cas", "N/A")
    )
    df["Precursor_A_Vendors"] = df["Precursor_A_SMILES"].apply(
        lambda s: ", ".join(sourcing_info.get(s, {}).get("vendors", []))
    )
    df["Precursor_A_Purchasable"] = df["Precursor_A_SMILES"].apply(
        lambda s: sourcing_info.get(s, {}).get("available", False)
    )

    df["Precursor_B_CAS"] = df["Precursor_B_SMILES"].apply(
        lambda s: sourcing_info.get(s, {}).get("cas", "N/A")
    )
    df["Precursor_B_Vendors"] = df["Precursor_B_SMILES"].apply(
        lambda s: ", ".join(sourcing_info.get(s, {}).get("vendors", []))
    )
    df["Precursor_B_Purchasable"] = df["Precursor_B_SMILES"].apply(
        lambda s: sourcing_info.get(s, {}).get("available", False)
    )

    df["Both_Precursors_Purchasable"] = (
        df["Precursor_A_Purchasable"] & df["Precursor_B_Purchasable"]
    )
    df.sort_values(
        by=["Both_Precursors_Purchasable"], ascending=False, inplace=True
    )

  return df


def recap_decompose_molecule(target_smiles: str) -> List[str]:
  """Hierarchical RECAP decomposition into building block leaves."""
  mol = Chem.MolFromSmiles(target_smiles)
  if mol is None:
    raise ValueError(f"Invalid SMILES: {target_smiles}")
  tree = Recap.RecapDecompose(mol)
  return list(tree.GetLeaves().keys()) if tree else []


def brics_decompose_molecule(target_smiles: str) -> List[str]:
  """BRICS decomposition into 16-motif labeled synthons."""
  mol = Chem.MolFromSmiles(target_smiles)
  if mol is None:
    raise ValueError(f"Invalid SMILES: {target_smiles}")
  return list(BRICS.BRICSDecompose(mol))


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Retrosynthetic Analysis, Strategic Disconnections, RECAP & BRICS."
      )
  )
  parser.add_argument(
      "--smiles", type=str, required=True, help="Target molecule SMILES string."
  )
  parser.add_argument(
      "--recap", action="store_true", help="Perform RECAP tree decomposition."
  )
  parser.add_argument(
      "--brics",
      action="store_true",
      help="Perform BRICS synthon decomposition.",
  )
  parser.add_argument(
      "--no_sourcing",
      action="store_true",
      help="Skip commercial vendor checking.",
  )
  parser.add_argument(
      "--output_csv",
      type=str,
      default=None,
      help="Path to save disconnection table CSV.",
  )
  args = parser.parse_args()

  print("\n" + "=" * 75)
  print(f"       RETROSYNTHETIC ANALYSIS REPORT FOR TARGET: {args.smiles}")
  print("=" * 75)

  if args.recap:
    leaves = recap_decompose_molecule(args.smiles)
    print(
        f"\nRECAP Hierarchical Decomposition ({len(leaves)} Building Block"
        " Leaves):"
    )
    for i, l in enumerate(leaves):
      print(f"  Leaf {i+1}: {l}")

  if args.brics:
    synthons = brics_decompose_molecule(args.smiles)
    print(f"\nBRICS Synthon Decomposition ({len(synthons)} Labeled Synthons):")
    for i, s in enumerate(synthons):
      print(f"  Synthon {i+1}: {s}")

  df = disconnect_target_molecule(
      args.smiles, check_availability=not args.no_sourcing
  )

  if df.empty:
    print(
        "\nNo standard 1-step strategic disconnections found for this target"
        " structure."
    )
    return

  print(f"\nStrategic 1-Step Disconnections ({len(df)} pathways identified):")
  cols_to_show = [
      "Disconnection_Reaction",
      "Precursor_A_SMILES",
      "Precursor_A_CAS",
      "Precursor_B_SMILES",
      "Precursor_B_CAS",
  ]
  if "Both_Precursors_Purchasable" in df.columns:
    cols_to_show.append("Both_Precursors_Purchasable")

  print(df[cols_to_show].to_string(index=False))

  if args.output_csv:
    df.to_csv(args.output_csv, index=False)
    print(
        "\nSaved full retrosynthetic disconnection plan to:"
        f" {os.path.abspath(args.output_csv)}"
    )


if __name__ == "__main__":
  main()
