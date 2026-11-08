# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "rdkit>=2024.3.1",
#   "beautifulsoup4>=4.12.0",
#   "pandas>=2.2.0",
#   "polite-http",
# ]
# ///

"""NIST Chemistry WebBook Thermochemical & Physical Property Client.

Queries the official NIST Chemistry WebBook
(https://webbook.nist.gov/chemistry/) by CAS Registry Number, InChIKey, SMILES,
or chemical name using polite_http.http_client.

Extracts:
- Gas-phase thermochemistry (Delta_f H_gas, S_gas, Cp)
- Condensed-phase thermochemistry (Delta_f H_liquid, Delta_f H_solid)
- Phase change properties (T_boil, T_fus, T_crit, P_crit, Delta_vap H, Antoine
  A/B/C)
- Ion energetics (Ionization Energy IE, Proton Affinity PA, Gas Basicity)
- Official NIST 2D MOL structure and validated RDKit canonical SMILES
"""

import argparse
import json
import os
import re
import sys
from typing import Any, Dict, List, Optional
import urllib.parse

from bs4 import BeautifulSoup
import pandas as pd
from polite_http import http_client
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

USER_AGENT = ""


class NistWebBookClient:
  """Rate-limited client for the NIST Chemistry WebBook."""

  BASE_URL = "https://webbook.nist.gov/"

  def __init__(self, qps: float = 2.0, timeout: float = 15.0):
    self._client = http_client.HttpClient(
        self.BASE_URL,
        qps=qps,
        timeout=timeout,
        user_agent=USER_AGENT or None,
    )

  def _get_text(self, path: str) -> Optional[str]:
    try:
      return self._client.fetch_text(path)
    except Exception as e:
      print(f"Error: polite_http request failed ({e})", file=sys.stderr)
      return None

  def resolve_nist_id(
      self,
      cas: Optional[str] = None,
      inchikey: Optional[str] = None,
      smiles: Optional[str] = None,
      name: Optional[str] = None,
  ) -> Optional[str]:
    """Resolves query identifiers to a canonical NIST WebBook ID (e.g., C64175)."""
    if smiles and not inchikey:
      mol = Chem.MolFromSmiles(smiles)
      if mol:
        inchikey = Chem.MolToInchiKey(mol)

    if cas:
      clean_cas = re.sub(r"[^0-9]", "", cas)
      return f"C{clean_cas}"

    query_params = {"Units": "SI"}
    if inchikey:
      query_params["InChI"] = inchikey
    elif name:
      query_params["Name"] = name
    else:
      raise ValueError("Must provide --cas, --inchikey, --smiles, or --name.")

    search_path = f"cgi/cbook.cgi?{urllib.parse.urlencode(query_params)}"
    html = self._get_text(search_path)
    if not html:
      return None

    # Check if search directly landed on a compound page containing a NIST ID link
    m = re.search(r"ID=(C[0-9]+)", html)
    if m:
      return m.group(1)
    return None

  def fetch_mol_structure(self, nist_id: str) -> Optional[str]:
    """Downloads the official 2D MOL structure block from NIST."""
    path = f"cgi/cbook.cgi?Str2File={urllib.parse.quote(nist_id)}"
    mol_block = self._get_text(path)
    if mol_block and ("M  END" in mol_block or "V2000" in mol_block):
      return mol_block
    return None

  def fetch_thermochemical_data(
      self,
      cas: Optional[str] = None,
      inchikey: Optional[str] = None,
      smiles: Optional[str] = None,
      name: Optional[str] = None,
  ) -> Dict[str, Any]:
    """Fetches and parses comprehensive thermochemical and phase-change tables from NIST WebBook."""
    nist_id = self.resolve_nist_id(
        cas=cas, inchikey=inchikey, smiles=smiles, name=name
    )
    if not nist_id:
      raise ValueError(
          f"Could not resolve compound in NIST Chemistry WebBook (cas={cas},"
          f" inchikey={inchikey}, name={name})."
      )

    path = f"cgi/cbook.cgi?ID={urllib.parse.quote(nist_id)}&Units=SI&Mask=1FFF"
    html = self._get_text(path)
    if not html:
      raise RuntimeError(
          f"Failed to retrieve NIST WebBook page for ID={nist_id}"
      )

    soup = BeautifulSoup(html, "html.parser")
    title = soup.find("h1", id="Top")
    compound_name = title.get_text(strip=True) if title else (name or nist_id)

    # Extract metadata bullets (Formula, MW, CAS, InChIKey)
    formula = None
    mw = None
    cas_found = cas
    for li in soup.find_all("li"):
      txt = li.get_text(" ", strip=True)
      if txt.startswith("Formula:"):
        formula = txt.split("Formula:", 1)[1].strip()
      elif txt.startswith("Molecular weight:"):
        try:
          mw = float(txt.split("Molecular weight:", 1)[1].strip().split()[0])
        except Exception:
          pass
      elif txt.startswith("CAS Registry Number:"):
        cas_found = txt.split("CAS Registry Number:", 1)[1].strip()

    # Parse all data tables
    thermochem_entries: List[Dict[str, str]] = []
    antoine_entries: List[Dict[str, str]] = []

    for table in soup.find_all("table", class_="data"):
      caption_el = table.find("caption") or table.find_previous(["h2", "h3"])
      section_name = (
          caption_el.get_text(" ", strip=True) if caption_el else "Data"
      )

      headers = [th.get_text(" ", strip=True) for th in table.find_all("th")]
      for row in table.find_all("tr"):
        cols = [td.get_text(" ", strip=True) for td in row.find_all("td")]
        if not cols:
          continue
        if "Quantity" in headers and len(cols) >= 3:
          thermochem_entries.append({
              "section": section_name,
              "quantity": cols[0],
              "value": cols[1],
              "units": cols[2],
              "method": cols[3] if len(cols) > 3 else "",
              "reference": cols[4] if len(cols) > 4 else "",
          })
        elif (
            "A" in headers
            and "B" in headers
            and "C" in headers
            and len(cols) >= 4
        ):
          antoine_entries.append({
              "temperature_range_K": cols[0],
              "A": cols[1],
              "B": cols[2],
              "C": cols[3],
              "reference": cols[4] if len(cols) > 4 else "",
          })

    # Download NIST 2D MOL block and convert to RDKit canonical SMILES
    mol_block = self.fetch_mol_structure(nist_id)
    canonical_smiles = smiles
    if mol_block:
      mol = Chem.MolFromMolBlock(mol_block)
      if mol:
        canonical_smiles = Chem.MolToSmiles(mol)
        if not inchikey:
          inchikey = Chem.MolToInchiKey(mol)
        if not formula:
          formula = rdMolDescriptors.CalcMolFormula(mol)
        if not mw:
          mw = round(Descriptors.MolWt(mol), 3)

    # Summarize key properties for fast programmatic access
    summary: Dict[str, Any] = {}
    for entry in thermochem_entries:
      q = entry["quantity"]
      val_str = entry["value"].split("±")[0].strip()
      try:
        val_num = float(val_str)
      except Exception:
        continue
      if "f H° gas" in q and "dfH_gas_kJ_mol" not in summary:
        summary["dfH_gas_kJ_mol"] = val_num
      elif "f H° liquid" in q and "dfH_liquid_kJ_mol" not in summary:
        summary["dfH_liquid_kJ_mol"] = val_num
      elif "f H° solid" in q and "dfH_solid_kJ_mol" not in summary:
        summary["dfH_solid_kJ_mol"] = val_num
      elif "S° gas" in q and "S_gas_J_mol_K" not in summary:
        summary["S_gas_J_mol_K"] = val_num
      elif "T boil" in q and "T_boil_K" not in summary:
        summary["T_boil_K"] = val_num
      elif "T fus" in q and "T_fus_K" not in summary:
        summary["T_fus_K"] = val_num
      elif "vap H°" in q and "dvapH_kJ_mol" not in summary:
        summary["dvapH_kJ_mol"] = val_num
      elif "IE (evaluated)" in q or (
          q == "IE" and "ionization_energy_eV" not in summary
      ):
        summary["ionization_energy_eV"] = val_num

    return {
        "nist_id": nist_id,
        "compound_name": compound_name,
        "cas": cas_found,
        "formula": formula,
        "molecular_weight": mw,
        "inchikey": inchikey,
        "canonical_smiles": canonical_smiles,
        "summary_properties": summary,
        "thermochemical_table": thermochem_entries,
        "antoine_parameters": antoine_entries,
        "mol_block": mol_block,
        "webbook_url": (
            f"https://webbook.nist.gov/cgi/cbook.cgi?ID={nist_id}&Units=SI"
        ),
    }


def main() -> None:
  parser = argparse.ArgumentParser(
      description=(
          "Query NIST Chemistry WebBook for experimental thermochemical,"
          " phase-change, and physical properties."
      )
  )
  parser.add_argument(
      "--cas",
      type=str,
      default=None,
      help="CAS Registry Number (e.g., 64-17-5)",
  )
  parser.add_argument(
      "--inchikey", type=str, default=None, help="Standard InChIKey"
  )
  parser.add_argument(
      "--smiles",
      type=str,
      default=None,
      help="SMILES string (automatically converted to InChIKey)",
  )
  parser.add_argument(
      "--name", type=str, default=None, help="Chemical name (e.g., ethanol)"
  )
  parser.add_argument(
      "--output-json",
      "--output_json",
      dest="output_json",
      type=str,
      default=None,
      help="Save full JSON report to path",
  )
  parser.add_argument(
      "--output-csv",
      "--output_csv",
      dest="output_csv",
      type=str,
      default=None,
      help="Save thermochemical table to CSV path",
  )
  parser.add_argument(
      "--output-mol",
      "--output_mol",
      dest="output_mol",
      type=str,
      default=None,
      help="Save official NIST 2D MOL file to path",
  )

  args = parser.parse_args()
  client = NistWebBookClient()
  data = client.fetch_thermochemical_data(
      cas=args.cas,
      inchikey=args.inchikey,
      smiles=args.smiles,
      name=args.name,
  )

  print(
      f"=== NIST Chemistry WebBook: {data['compound_name']}"
      f" ({data['nist_id']}) ==="
  )
  print(f"CAS Registry No:   {data['cas']}")
  print(
      f"Molecular Formula: {data['formula']} (MW: {data['molecular_weight']})"
  )
  print(f"Canonical SMILES:  {data['canonical_smiles']}")
  print(f"NIST WebBook URL:  {data['webbook_url']}")
  print("--- Key Evaluated Thermochemical Properties ---")
  for k, v in data["summary_properties"].items():
    print(f"  {k:24s}: {v}")
  print(f"Total Thermochemical Entries: {len(data['thermochemical_table'])}")
  if data["antoine_parameters"]:
    print(
        f"Antoine Vapor Pressure Equations: {len(data['antoine_parameters'])}"
        " temperature ranges available"
    )

  if args.output_mol and data.get("mol_block"):
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_mol)), exist_ok=True
    )
    with open(args.output_mol, "w", encoding="utf-8") as f:
      f.write(data["mol_block"])
    print(f"[+] Saved official NIST 2D MOL structure to: {args.output_mol}")

  if args.output_csv and data["thermochemical_table"]:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True
    )
    df = pd.DataFrame(data["thermochemical_table"])
    df.to_csv(args.output_csv, index=False)
    print(
        f"[+] Saved thermochemical CSV table ({len(df)} rows) to:"
        f" {args.output_csv}"
    )

  if args.output_json:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_json)), exist_ok=True
    )
    out_copy = dict(data)
    out_copy.pop("mol_block", None)  # Keep JSON compact
    with open(args.output_json, "w", encoding="utf-8") as f:
      json.dump(out_copy, f, indent=2)
    print(f"[+] Saved full NIST JSON report to: {args.output_json}")


if __name__ == "__main__":
  main()
