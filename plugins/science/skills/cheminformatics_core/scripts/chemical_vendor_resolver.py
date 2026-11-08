# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "polite-http",
#     "rdkit",
#     "numpy",
#     "pandas",
# ]
# ///

"""Unified Chemical Vendor, CAS Registry, and Commercial Availability Resolver.

Consolidates CAS validation, InChIKey round-trip structure verification,
and commercial vendor lookup (Enamine, Chemspace, MolPort, Mcule, PubChem).
"""

import argparse
import json
import os
from typing import Any, Dict, Optional
import urllib.parse

from polite_http import http_client
from rdkit import Chem

_PUBCHEM_CLIENT = http_client.HttpClient(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug", qps=5.0
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
]


def is_valid_cas(cas_str: str) -> bool:
  """Validates format and Mod-10 checksum of a CAS Registry Number.

  Format: [2-7 digits]-[2 digits]-[1 digit]
  Checksum: (sum(d_i * weight_i)) mod 10 == check_digit
  """
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
  check_digit = int(p3)
  total = sum(int(d) * (len(digits) - i) for i, d in enumerate(digits))
  return total % 10 == check_digit


def verify_cas_roundtrip(smiles: str, cas_number: str) -> Dict[str, Any]:
  """Round-trip verification: ensures a CAS number maps back to the expected SMILES structure via InChIKey."""
  if not is_valid_cas(cas_number):
    return {
        "verified": False,
        "reason": "INVALID_CAS_CHECKSUM",
        "expected_inchikey": None,
        "resolved_inchikey": None,
    }

  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    return {
        "verified": False,
        "reason": "INVALID_SMILES",
        "expected_inchikey": None,
        "resolved_inchikey": None,
    }

  expected_ik = Chem.MolToInchiKey(mol)
  enc_cas = urllib.parse.quote(cas_number.strip(), safe="")
  path = f"compound/name/{enc_cas}/property/InChIKey/JSON"

  resolved_ik = None
  try:
    data = _PUBCHEM_CLIENT.fetch_json(path, timeout=10)
    props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
    resolved_ik = props.get("InChIKey")
  except Exception:
    pass

  if not resolved_ik:
    return {
        "verified": False,
        "reason": "CAS_NOT_RESOLVED_IN_PUBCHEM",
        "expected_inchikey": expected_ik,
        "resolved_inchikey": None,
    }

  # Match primary connectivity layer (first 14 chars of InChIKey)
  verified = expected_ik[:14] == resolved_ik[:14]
  return {
      "verified": verified,
      "reason": "MATCH" if verified else "STEREOCHEMISTRY_OR_SKELETON_MISMATCH",
      "expected_inchikey": expected_ik,
      "resolved_inchikey": resolved_ik,
  }


def resolve_compound_suppliers(identifier: str) -> Dict[str, Any]:
  """Resolves chemical structure, CID, CAS numbers, and commercial supplier linkouts via PubChem PUG-REST."""
  mol = Chem.MolFromSmiles(identifier)
  is_smiles = mol is not None
  inchikey = Chem.MolToInchiKey(mol) if mol else None

  if is_smiles and inchikey:
    path = f"compound/inchikey/{inchikey}/property/CanonicalSMILES,IsomericSMILES,IUPACName,MolecularWeight,InChIKey/JSON"
  else:
    encoded_id = urllib.parse.quote(identifier, safe="")
    path = f"compound/name/{encoded_id}/property/CanonicalSMILES,IsomericSMILES,IUPACName,MolecularWeight,InChIKey/JSON"

  try:
    data = _PUBCHEM_CLIENT.fetch_json(path, timeout=12)
    props = data["PropertyTable"]["Properties"][0]
    cid = props["CID"]
    canonical_smiles = (
        props.get("ConnectivitySMILES")
        or props.get("CanonicalSMILES")
        or props.get("SMILES")
    )
    inchikey = props.get("InChIKey", "")
    iupac_name = props.get("IUPACName", "")
  except Exception as e:
    return {
        "query": identifier,
        "error": f"Failed to resolve compound: {str(e)}",
    }

  # Fetch synonyms / CAS
  cas_numbers = []
  syn_data = None
  try:
    syn_data = _PUBCHEM_CLIENT.fetch_json(
        f"compound/cid/{cid}/synonyms/JSON", timeout=10
    )
  except http_client.HttpError as e:
    # HTTP 404 is expected when a CID has no indexed synonyms.
    if e.status_code != 404:
      sys.stderr.write(
          f"[chemical_vendor_resolver] Warning: HTTP {e.status_code} fetching"
          f" synonyms for CID {cid}: {e}\n"
      )
  except Exception as e:
    sys.stderr.write(
        "[chemical_vendor_resolver] Warning: failed to fetch synonyms for CID"
        f" {cid}: {e}\n"
    )
  if isinstance(syn_data, dict):
    all_syns = (
        syn_data.get("InformationList", {})
        .get("Information", [{}])[0]
        .get("Synonym", [])
    )
    for s in all_syns:
      if is_valid_cas(s) and s not in cas_numbers:
        cas_numbers.append(s)

  primary_cas = cas_numbers[0] if cas_numbers else "NOT_FOUND"

  # Query vendor sources
  matched_vendors = []
  src_data = None
  try:
    src_data = _PUBCHEM_CLIENT.fetch_json(
        f"compound/cid/{cid}/xrefs/SourceName/JSON", timeout=10
    )
  except http_client.HttpError as e:
    # HTTP 404 is expected when a CID has no external source records.
    if e.status_code != 404:
      sys.stderr.write(
          f"[chemical_vendor_resolver] Warning: HTTP {e.status_code} fetching"
          f" vendor sources for CID {cid}: {e}\n"
      )
  except Exception as e:
    sys.stderr.write(
        "[chemical_vendor_resolver] Warning: failed to fetch vendor sources"
        f" for CID {cid}: {e}\n"
    )
  if isinstance(src_data, dict):
    source_names = (
        src_data.get("InformationList", {})
        .get("Information", [{}])[0]
        .get("SourceName", [])
    )
    source_names_lower = [str(s).lower() for s in source_names]
    for vendor in COMMERCIAL_VENDORS:
      if any(vendor.lower() in s for s in source_names_lower):
        matched_vendors.append(vendor)

  return {
      "query": identifier,
      "cid": cid,
      "canonical_smiles": canonical_smiles,
      "inchikey": inchikey,
      "iupac_name": iupac_name,
      "primary_cas": primary_cas,
      "all_cas_numbers": cas_numbers,
      "available_commercial_vendors": matched_vendors,
      "enamine_search_url": (
          f"https://enamine.net/hit-locator?search={urllib.parse.quote(canonical_smiles)}"
      ),
      "chemspace_search_url": (
          f"https://chem-space.com/search?query={urllib.parse.quote(canonical_smiles)}"
      ),
      "molport_search_url": (
          f"https://www.molport.com/shop/find-compounds?search={urllib.parse.quote(canonical_smiles)}"
      ),
  }


def _write_or_print(res: Dict[str, Any], output_path: Optional[str]) -> None:
  """Writes full result dictionary to JSON file and prints a concise summary to stdout to avoid context window bloat."""
  if output_path:
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
      json.dump(res, f, indent=2)
  vendors = (
      res.get("available_commercial_vendors")
      or res.get("vendors")
      or res.get("suppliers")
      or []
  )
  summary_id = res.get("primary_cas") or res.get("cas") or res.get("query", "")
  print(f"=== Chemical Vendor & CAS Resolution Summary ({summary_id}) ===")
  if "cid" in res:
    print(f"PubChem CID:     {res.get('cid')}")
  if "canonical_smiles" in res:
    print(f"SMILES:          {res.get('canonical_smiles')}")
  if "primary_cas" in res:
    print(
        f"Primary CAS:     {res.get('primary_cas')}"
        f" (status={res.get('status')})"
    )
  if vendors:
    print(f"Vendors Found:   {len(vendors)} (top 5 shown below)")
    for v in vendors[:5]:
      if isinstance(v, dict):
        v_name = v.get("vendor_name") or v.get("name", "Unknown")
        cat_id = v.get("catalog_id") or v.get("registry_id", "N/A")
        print(f"  - {v_name}: {cat_id}")
      else:
        print(f"  - {v}")
  if output_path:
    print(f"Full JSON Payload Saved: {os.path.abspath(output_path)}")
  elif len(json.dumps(res)) <= 1200:
    print(json.dumps(res, indent=2))
  else:
    print(
        "(Full vendor payload omitted from stdout; pass --output-json <path> to"
        " save full catalog JSON)"
    )


def main():
  parser = argparse.ArgumentParser(
      description="Unified Chemical Vendor and CAS Resolver."
  )
  parser.add_argument("--cas", default=None, help="CAS Registry Number.")
  parser.add_argument(
      "--smiles",
      default="",
      help="Expected SMILES for round-trip verification.",
  )
  parser.add_argument(
      "--query", default=None, help="SMILES, chemical name, or CID."
  )
  parser.add_argument(
      "--output",
      "--output_json",
      "--output-json",
      dest="output",
      default=None,
      help="Optional output JSON file path.",
  )
  subparsers = parser.add_subparsers(dest="subcommand", required=False)

  # Subcommand: resolve-cas
  cas_p = subparsers.add_parser(
      "resolve-cas", help="Verify and resolve CAS Registry Number."
  )
  cas_p.add_argument("--cas", required=True, help="CAS Registry Number.")
  cas_p.add_argument(
      "--smiles",
      default="",
      help="Expected SMILES for round-trip verification.",
  )
  cas_p.add_argument(
      "--output",
      "--output_json",
      "--output-json",
      dest="output",
      default=None,
      help="Optional output JSON file path.",
  )

  # Subcommand: lookup-suppliers
  sup_p = subparsers.add_parser(
      "lookup-suppliers",
      help="Look up commercial chemical vendors and catalogs.",
  )
  sup_p.add_argument(
      "--query", required=True, help="SMILES, chemical name, or CID."
  )
  sup_p.add_argument(
      "--output",
      "--output_json",
      "--output-json",
      dest="output",
      default=None,
      help="Optional output JSON file path.",
  )

  args = parser.parse_args()
  mode = args.subcommand
  if not mode:
    if args.cas:
      mode = "resolve-cas"
    elif args.query:
      mode = "lookup-suppliers"
    else:
      parser.error(
          "Specify a subcommand (resolve-cas / lookup-suppliers) or pass --cas"
          " or --query."
      )

  if mode == "resolve-cas":
    if args.smiles:
      res = verify_cas_roundtrip(args.smiles, args.cas)
    else:
      valid = is_valid_cas(args.cas)
      res = {"cas": args.cas, "valid_format_and_checksum": valid}
    _write_or_print(res, args.output)

  elif mode == "lookup-suppliers":
    res = resolve_compound_suppliers(args.query)
    _write_or_print(res, args.output)


if __name__ == "__main__":
  main()
