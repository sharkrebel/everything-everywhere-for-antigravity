# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "pandas",
#     "polite-http",
# ]
# ///

"""Programmatic verification of chemical nomenclature, stereochemistry, and topological descriptors.

Verifies that any chemical names, stereocenters ((R)/(S)/(E)/(Z)), macrocyclic
ring sizes, and topological attributes cited in reports, synthesis proposals, or
compound series strictly match the underlying molecular graph in RDKit.

Prevents hallucinated nomenclature (e.g. misstating ring sizes, misassigning
chiral descriptors, or proposing incorrect systematic IUPAC names).
"""

import argparse
import json
import os
import sys
from typing import Any, Dict, List, Optional, Tuple
import urllib.parse

import pandas as pd
from polite_http import http_client
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

OPSIN_BASE_URL = "https://opsin.ch.cam.ac.uk"
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov"

# Rate-limited, retrying clients. PubChem asks for <= 5 requests/second.
_OPSIN_CLIENT = http_client.HttpClient(OPSIN_BASE_URL, qps=3.0)
_PUBCHEM_CLIENT = http_client.HttpClient(PUBCHEM_BASE_URL, qps=5.0)


def analyze_molecular_topology(mol: Chem.Mol) -> Dict[str, Any]:
  """Analyzes rings, stereochemistry, and topological features of a molecule."""
  if mol is None:
    return {"valid": False, "error": "Mol is None"}

  # 1. Ring perception
  ring_info = mol.GetRingInfo()
  atom_rings = list(ring_info.AtomRings())
  ring_sizes = [len(r) for r in atom_rings]
  macrocyclic_rings = [len(r) for r in atom_rings if len(r) >= 12]

  ring_details = []
  for idx, r in enumerate(atom_rings):
    atoms = [mol.GetAtomWithIdx(i) for i in r]
    symbols = [a.GetSymbol() for a in atoms]
    heteroatoms = [s for s in symbols if s != "C"]
    is_aromatic = all(a.GetIsAromatic() for a in atoms)
    ring_details.append({
        "ring_index": idx + 1,
        "size": len(r),
        "is_macrocycle": len(r) >= 12,
        "is_aromatic": is_aromatic,
        "heteroatom_count": len(heteroatoms),
        "heteroatoms": heteroatoms,
        "atom_indices": list(r),
    })

  # 2. Chiral centers
  chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
  chiral_details = [
      {"atom_idx": idx, "descriptor": desc, "assigned": desc in ["R", "S"]}
      for idx, desc in chiral_centers
  ]

  # 3. Double bond stereochemistry
  double_bonds = []
  for b in mol.GetBonds():
    if b.GetBondType() == Chem.BondType.DOUBLE:
      stereo = str(b.GetStereo())
      double_bonds.append({
          "begin_idx": b.GetBeginAtomIdx(),
          "end_idx": b.GetEndAtomIdx(),
          "stereo": stereo,
      })

  # 4. Core descriptors
  formula = rdMolDescriptors.CalcMolFormula(mol)
  mw = Descriptors.MolWt(mol)
  heavy_atoms = mol.GetNumHeavyAtoms()
  rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
  inchikey = Chem.MolToInchiKey(mol)
  canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

  return {
      "valid": True,
      "canonical_smiles": canonical_smiles,
      "inchikey": inchikey,
      "formula": formula,
      "molecular_weight": round(mw, 3),
      "heavy_atoms": heavy_atoms,
      "rotatable_bonds": rotatable_bonds,
      "ring_count": len(atom_rings),
      "ring_sizes": ring_sizes,
      "macrocycle_sizes": macrocyclic_rings,
      "rings": ring_details,
      "chiral_centers": chiral_details,
      "double_bonds": double_bonds,
  }


def verify_ring_size(mol: Chem.Mol, claimed_size: int) -> Tuple[bool, str]:
  """Verifies whether a claimed ring size (e.g. 14-membered ring) is present."""
  if mol is None:
    return False, "Molecule is None."
  atom_rings = list(mol.GetRingInfo().AtomRings())
  sizes = [len(r) for r in atom_rings]
  if claimed_size in sizes:
    matching_count = sizes.count(claimed_size)
    return (
        True,
        (
            f"Verified: Molecule contains {matching_count} ring(s) of size"
            f" {claimed_size}."
        ),
    )
  return (
      False,
      (
          f"MISMATCH: Claimed {claimed_size}-membered ring, but detected ring"
          f" sizes are: {sorted(sizes)}."
      ),
  )


def verify_stereocenters(
    mol: Chem.Mol, claimed_chirality: List[str]
) -> Tuple[bool, str]:
  """Verifies that claimed stereochemical designations (e.g.

  ['(R)']) match detected centers.
  """
  if mol is None:
    return False, "Molecule is None."
  detected = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
  detected_descs = [desc for _, desc in detected]

  missing = []
  for c in claimed_chirality:
    clean_c = c.strip("() ")
    if clean_c.upper() in (
        "R/S",
        "S/R",
        "R,S",
        "S,R",
        "RS",
        "SR",
        "±",
        "+/-",
        "RAC",
    ):
      if detected:
        continue
      missing.append(c)
    elif clean_c not in detected_descs:
      missing.append(c)

  if missing:
    return (
        False,
        (
            f"MISMATCH: Expected stereocenters {missing} not found. Detected"
            f" chiral centers: {detected}"
        ),
    )
  return True, f"Verified: Stereocenters match {detected}."


def resolve_name_with_opsin(name: str) -> Optional[Dict[str, Any]]:
  """Resolves an IUPAC / systematic name to structure via the OPSIN web API."""
  encoded_name = urllib.parse.quote(name.strip())
  url = f"{OPSIN_BASE_URL}/opsin/{encoded_name}.json"
  try:
    data = _OPSIN_CLIENT.fetch_json(url)
  except http_client.HttpError as e:
    # HTTP 400/404 is expected when OPSIN cannot parse a non-systematic name.
    if e.status_code not in (400, 404):
      sys.stderr.write(
          f"[verify_chemical_nomenclature] Warning: OPSIN HTTP {e.status_code}"
          f" for {name!r}: {e}\n"
      )
    return None
  except Exception as e:
    sys.stderr.write(
        "[verify_chemical_nomenclature] Warning: OPSIN lookup failed for"
        f" {name!r}: {e}\n"
    )
    return None
  if not isinstance(data, dict):
    return None
  ikey = (
      data.get("stdinchikey") or data.get("stdinchiKey") or data.get("inchikey")
  )
  if not ikey:
    return None
  return {
      "smiles": data.get("smiles"),
      "inchikey": ikey,
      "resolved_by": "OPSIN",
  }


def resolve_name_with_pubchem(name: str) -> Optional[Dict[str, Any]]:
  """Resolves a chemical name to structure via PubChem PUG REST API."""
  encoded_name = urllib.parse.quote(name.strip())
  url = (
      f"{PUBCHEM_BASE_URL}/rest/pug/compound/name/{encoded_name}"
      "/property/InChIKey,IsomericSMILES/JSON"
  )
  try:
    data = _PUBCHEM_CLIENT.fetch_json(url)
  except http_client.HttpError as e:
    if e.status_code not in (400, 404):
      sys.stderr.write(
          "[verify_chemical_nomenclature] Warning: PubChem HTTP"
          f" {e.status_code} for {name!r}: {e}\n"
      )
    return None
  except Exception as e:
    sys.stderr.write(
        "[verify_chemical_nomenclature] Warning: PubChem lookup failed for"
        f" {name!r}: {e}\n"
    )
    return None
  props_list = (
      data.get("PropertyTable", {}).get("Properties", [])
      if isinstance(data, dict)
      else []
  )
  if not props_list:
    return None
  props = props_list[0]
  return {
      "smiles": (
          props.get("IsomericSMILES")
          or props.get("SMILES")
          or props.get("CanonicalSMILES")
          or props.get("ConnectivitySMILES")
      ),
      "inchikey": props.get("InChIKey"),
      "resolved_by": "PubChem",
  }


def _has_assigned_stereo(mol: Optional[Chem.Mol]) -> bool:
  """Returns whether the molecule has any specified stereocenters or stereo bonds."""
  if mol is None:
    return False
  return any(
      s.specified == Chem.StereoSpecified.Specified
      for s in Chem.FindPotentialStereo(mol)
  )


def verify_iupac_name(
    mol: Chem.Mol, claimed_name: str
) -> Tuple[bool, str, Dict[str, Any]]:
  """Verifies that a claimed systematic IUPAC name resolves round-trip to the exact molecular structure."""
  if mol is None:
    return False, "Molecule is None.", {}

  mol_inchikey = Chem.MolToInchiKey(mol)
  res = resolve_name_with_opsin(claimed_name) or resolve_name_with_pubchem(
      claimed_name
  )

  if res is None:
    return (
        False,
        (
            f"UNRESOLVED: The claimed name '{claimed_name}' could not be parsed"
            " by OPSIN or PubChem."
        ),
        {},
    )

  resolved_inchikey = res.get("inchikey")
  if resolved_inchikey == mol_inchikey:
    return (
        True,
        (
            f"Verified: Name '{claimed_name}' resolves to exact structure"
            f" ({res['resolved_by']})."
        ),
        res,
    )

  # Check skeleton match (first 14 chars of InChIKey)
  if (
      resolved_inchikey
      and mol_inchikey
      and resolved_inchikey[:14] == mol_inchikey[:14]
  ):
    resolved_smiles = res.get("smiles")
    resolved_mol = (
        Chem.MolFromSmiles(resolved_smiles) if resolved_smiles else None
    )
    if (
        _has_assigned_stereo(resolved_mol)
        and resolved_inchikey[15:23] != mol_inchikey[15:23]
    ):
      return (
          False,
          (
              "MISMATCH: Constitutional skeleton matches, but stereochemistry"
              f" differs between claimed name '{claimed_name}'"
              f" ({resolved_inchikey}) and target molecule ({mol_inchikey})."
          ),
          res,
      )
    return (
        True,
        (
            "Verified (Constitutional skeleton matches; stereoisomer/tautomer"
            f" variant): '{claimed_name}'."
        ),
        res,
    )

  return (
      False,
      (
          f"MISMATCH: Claimed name '{claimed_name}' resolved to InChIKey"
          f" {resolved_inchikey}, but target molecule InChIKey is"
          f" {mol_inchikey}."
      ),
      res,
  )


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Programmatic Chemical Nomenclature & Topological Verification"
          " Utility."
      )
  )
  parser.add_argument(
      "--smiles",
      type=str,
      nargs="+",
      action="extend",
      help="One or more SMILES strings of target molecule(s).",
  )
  parser.add_argument(
      "--input_csv",
      type=str,
      help="Input CSV file to verify batch of molecules.",
  )
  parser.add_argument(
      "--smiles_col",
      type=str,
      default="SMILES",
      help="SMILES column name in CSV.",
  )
  parser.add_argument(
      "--claimed_name",
      type=str,
      help="Claimed IUPAC or systematic name to verify.",
  )
  parser.add_argument(
      "--claimed_ring_size",
      type=int,
      help="Claimed ring size (e.g. 14 for 14-membered macrocycle).",
  )
  parser.add_argument(
      "--claimed_chirality",
      type=str,
      nargs="+",
      help="Claimed chiral center descriptors (e.g. 'R' 'S').",
  )
  parser.add_argument(
      "--json", action="store_true", help="Output results as JSON."
  )
  parser.add_argument(
      "--output",
      "-o",
      type=str,
      default=None,
      help="Output path to save verification report (.json or .csv).",
  )
  args = parser.parse_args()

  if args.smiles:
    all_outputs = []
    for smi in args.smiles:
      mol = Chem.MolFromSmiles(smi)
      if mol is None:
        print(f"Error: Invalid SMILES: {smi}", file=sys.stderr)
        sys.exit(1)

      topology = analyze_molecular_topology(mol)
      verifications = {}

      if args.claimed_ring_size is not None:
        ok, msg = verify_ring_size(mol, args.claimed_ring_size)
        verifications["ring_size"] = {"valid": ok, "message": msg}

      if args.claimed_chirality:
        ok, msg = verify_stereocenters(mol, args.claimed_chirality)
        verifications["chirality"] = {"valid": ok, "message": msg}

      if args.claimed_name:
        ok, msg, res = verify_iupac_name(mol, args.claimed_name)
        verifications["iupac_name"] = {
            "valid": ok,
            "message": msg,
            "resolution": res,
        }

      all_outputs.append({
          "smiles": smi,
          "topology": topology,
          "verifications": verifications,
      })

    output_data = all_outputs[0] if len(all_outputs) == 1 else all_outputs

    if args.output:
      os.makedirs(
          os.path.dirname(os.path.abspath(args.output)) or ".", exist_ok=True
      )
      with open(args.output, "w", encoding="utf-8") as f:
        json.dump(output_data, f, indent=2)

    if args.json:
      print(json.dumps(output_data, indent=2))
    else:
      for item in all_outputs:
        topology = item["topology"]
        verifications = item["verifications"]
        print("=== MOLECULAR TOPOLOGY & NOMENCLATURE REPORT ===")
        if len(all_outputs) > 1:
          print(f"SMILES:           {item['smiles']}")
        print(f"Formula:          {topology['formula']}")
        print(f"Molecular Weight: {topology['molecular_weight']}")
        print(f"InChIKey:         {topology['inchikey']}")
        print(f"Ring Sizes:       {topology['ring_sizes']}")
        print(f"Macrocycles:      {topology['macrocycle_sizes'] or 'None'}")
        print(f"Chiral Centers:   {topology['chiral_centers'] or 'None'}")

        if verifications:
          print("\n--- Verification Findings ---")
          for k, v in verifications.items():
            status = "[PASS]" if v["valid"] else "[FAIL]"
            print(f"{status} {k.upper()}: {v['message']}")

  elif args.input_csv:
    df = pd.read_csv(args.input_csv)
    col = args.smiles_col
    if col not in df.columns:
      for cand in (
          "canonical_smiles",
          "smiles",
          "SMILES",
          "Clean_SMILES",
          "Structure",
      ):
        if cand in df.columns:
          col = cand
          break
    if col not in df.columns:
      print(
          f"Error: Column {col} not found in {args.input_csv}. Available:"
          f" {list(df.columns)}",
          file=sys.stderr,
      )
      sys.exit(1)

    id_col = next(
        (
            c
            for c in (
                "compound_id",
                "paper_compound_id",
                "patent_example_id",
                "compound_name",
                "chembl_id",
                "Name",
                "ID",
            )
            if c in df.columns
        ),
        None,
    )

    records = []
    for idx, row in df.iterrows():
      smi = row[col]
      mol = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None
      top = analyze_molecular_topology(mol) if mol else {}
      records.append({
          "Index": idx,
          "Compound_ID": (
              str(row[id_col])
              if id_col and pd.notna(row.get(id_col))
              else f"Row_{idx}"
          ),
          "SMILES": smi,
          "Formula": top.get("formula"),
          "MW": top.get("molecular_weight"),
          "Ring_Sizes": str(top.get("ring_sizes")),
          "Macrocycle_Sizes": str(top.get("macrocycle_sizes")),
          "Chiral_Centers": str(
              [c["descriptor"] for c in top.get("chiral_centers", [])]
          ),
      })

    res_df = pd.DataFrame(records)
    if args.output:
      os.makedirs(
          os.path.dirname(os.path.abspath(args.output)) or ".", exist_ok=True
      )
      if args.output.endswith(".json"):
        res_df.to_json(args.output, orient="records", indent=2)
      else:
        res_df.to_csv(args.output, index=False)
    if args.json:
      print(res_df.to_json(orient="records", indent=2))
    else:
      print(res_df.to_string(index=False))
  else:
    parser.error("Either --smiles or --input_csv is required.")


if __name__ == "__main__":
  main()
