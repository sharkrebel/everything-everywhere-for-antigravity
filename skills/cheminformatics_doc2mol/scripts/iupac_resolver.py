# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pandas",
#     "polite-http",
#     "rdkit",
# ]
# ///

"""Multi-Tier Chemical Nomenclature & IUPAC Resolver.

Translates systematic IUPAC names, trivial trade names, developmental codes,
and CAS numbers into validated, standardized RDKit Mol objects.
Cascade: Local/Online OPSIN -> PubChem PUG REST -> NCI CIR.
"""

import argparse
from collections.abc import Mapping
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Optional, Tuple
import urllib.parse

import pandas as pd
from polite_http import http_client
from rdkit import Chem

_OPSIN_CLIENT = http_client.HttpClient(
    "https://opsin.ch.cam.ac.uk/opsin/", qps=5.0, max_retries=2, timeout=8.0
)
_PUBCHEM_CLIENT = http_client.HttpClient(
    "https://pubchem.ncbi.nlm.nih.gov/rest/",
    qps=5.0,
    max_retries=2,
    timeout=8.0,
)
_CIR_CLIENT = http_client.HttpClient(
    "https://cactus.nci.nih.gov/chemical/structure/",
    qps=5.0,
    max_retries=2,
    backoff_base=2.0,
    backoff_max=10.0,
    timeout=8.0,
)

# Import text cleaner from same package
_cur_dir = os.path.dirname(os.path.abspath(__file__))
if _cur_dir not in sys.path:
  sys.path.insert(0, _cur_dir)
from chemical_text_sanitizer import clean_chemical_text
from chemical_text_sanitizer import normalize_smiles_string
from chemical_text_sanitizer import repair_chemical_ocr_text


def extract_pubchem_smiles(prop: Mapping[str, Any]) -> Optional[str]:
  """Selects the preferred SMILES from a PubChem PropertyTable record.

  Prefers stereo-preserving keys (`SMILES`, legacy `IsomericSMILES`) before
  falling back to connectivity-only keys (`ConnectivitySMILES`, legacy
  `CanonicalSMILES`).

  Args:
    prop: A single compound mapping from `PropertyTable.Properties`.

  Returns:
    Non-empty SMILES string if present, otherwise `None`.
  """
  for key in (
      "SMILES",
      "IsomericSMILES",
      "ConnectivitySMILES",
      "CanonicalSMILES",
  ):
    val = prop.get(key)
    if isinstance(val, str) and val.strip():
      return val.strip()
  return None


class IUPACResolver:
  """Multi-tier cascade chemical nomenclature resolver."""

  def __init__(self, local_opsin_jar: Optional[str] = None):
    self.local_opsin_jar = local_opsin_jar or self._find_opsin_jar()
    self.java_bin = shutil.which("java")
    if self.local_opsin_jar and not self.java_bin:
      sys.stderr.write(
          f"[OPSIN] Warning: local OPSIN JAR found at {self.local_opsin_jar},"
          " but 'java' binary is not on PATH; falling back to online OPSIN"
          " API.\n"
      )

  def _find_opsin_jar(self) -> Optional[str]:
    """Checks potential local paths for opsin jar."""
    candidates = [
        os.path.expanduser("~/.local/share/opsin/opsin-cli.jar"),
        "/usr/local/share/java/opsin.jar",
        "/usr/share/java/opsin.jar",
    ]
    for c in candidates:
      if c and os.path.exists(c):
        return c
    return None

  def resolve_with_opsin_local(self, name: str) -> Optional[str]:
    """Converts IUPAC via local Java OPSIN JAR if available."""
    if not self.local_opsin_jar or not self.java_bin:
      return None

    inp_path = ""
    out_path = ""
    try:
      with tempfile.NamedTemporaryFile(
          suffix=".txt", mode="w", delete=False
      ) as inp_f:
        inp_f.write(name)
        inp_path = inp_f.name

      with tempfile.NamedTemporaryFile(
          suffix=".txt", mode="r", delete=False
      ) as out_f:
        out_path = out_f.name

      cmd = [self.java_bin, "-jar", self.local_opsin_jar, inp_path, out_path]
      _ = subprocess.run(cmd, capture_output=True, timeout=10, check=False)
      smi = ""
      if os.path.exists(out_path):
        with open(out_path, "r", encoding="utf-8") as f:
          smi = f.read().strip()

      if smi and Chem.MolFromSmiles(smi) is not None:
        return smi
    except (OSError, subprocess.SubprocessError) as e:
      sys.stderr.write(f"[OPSIN_LOCAL] Subprocess failed for {name!r}: {e}\n")
    finally:
      for path in (inp_path, out_path):
        if path and os.path.exists(path):
          os.remove(path)
    return None

  def resolve_with_opsin_online(self, name: str) -> Optional[str]:
    """Converts IUPAC via Cambridge OPSIN REST API."""
    encoded = urllib.parse.quote(name)
    try:
      data = _OPSIN_CLIENT.fetch_json(f"{encoded}.json", timeout=8.0)
    except http_client.HttpError as e:
      # HTTP 400/404 is expected when OPSIN does not recognize a trivial/trade
      # name; return None quietly so the cascade falls through to PubChem.
      if e.status_code not in (400, 404):
        sys.stderr.write(
            f"[OPSIN] HTTP {e.status_code} resolving {name!r}: {e}\n"
        )
      return None
    except Exception as e:
      sys.stderr.write(f"[OPSIN] Resolution failed for {name!r}: {e}\n")
      return None
    smi = data.get("smiles") if isinstance(data, dict) else None
    if smi and Chem.MolFromSmiles(smi) is not None:
      return smi
    return None

  def resolve_with_pubchem(self, name: str) -> Optional[Tuple[str, str]]:
    """Resolves trivial names, trade names, or IUPAC via PubChem PUG REST."""
    encoded = urllib.parse.quote(name)
    try:
      data = _PUBCHEM_CLIENT.fetch_json(
          f"pug/compound/name/{encoded}/property/SMILES,ConnectivitySMILES,InChIKey/JSON",
          timeout=8.0,
      )
    except http_client.HttpError as e:
      # HTTP 404 is expected when PubChem has no record for the name; fall
      # through quietly to the next resolver tier.
      if e.status_code not in (400, 404):
        sys.stderr.write(
            f"[PubChem] HTTP {e.status_code} resolving {name!r}: {e}\n"
        )
      return None
    except Exception as e:
      sys.stderr.write(f"[PubChem] Resolution failed for {name!r}: {e}\n")
      return None
    props = (
        data.get("PropertyTable", {}).get("Properties", [])
        if isinstance(data, dict)
        else []
    )
    if props and isinstance(props[0], Mapping):
      smi = extract_pubchem_smiles(props[0])
      inchikey = props[0].get("InChIKey", "")
      if smi and Chem.MolFromSmiles(smi) is not None:
        return smi, inchikey
    return None

  def resolve_with_cir(self, name: str) -> Optional[str]:
    """Resolves via NIH Chemical Identifier Resolver (CIR)."""
    try:
      encoded = urllib.parse.quote(name)
      resp = _CIR_CLIENT.fetch(f"{encoded}/smiles", timeout=8.0)
      smi = resp.text.strip()
      if smi and Chem.MolFromSmiles(smi) is not None:
        return smi
    except Exception as e:
      sys.stderr.write(
          f"[CIR] Resolution failed for {name!r} via NCI CACTUS: {e}\n"
      )
    return None

  def resolve(self, name_or_iupac: str) -> Optional[Chem.Mol]:
    """Master cascade: cleans name, evaluates engines, and returns sanitized RDKit Mol."""
    cleaned = clean_chemical_text(name_or_iupac)
    if not cleaned:
      return None

    # Step 0: Check if input is already a valid SMILES
    smi_try = normalize_smiles_string(cleaned)
    mol_direct = Chem.MolFromSmiles(smi_try)
    if mol_direct is not None:
      mol_direct.SetProp("_Name", cleaned)
      mol_direct.SetProp("SourceText", cleaned)
      mol_direct.SetProp("ResolutionEngine", "DIRECT_SMILES")
      return mol_direct

    resolved_smi = None
    engine_used = None
    inchikey = ""

    # Step 1: Local OPSIN
    smi = self.resolve_with_opsin_local(cleaned)
    if smi:
      resolved_smi = smi
      engine_used = "OPSIN_LOCAL"

    # Step 2: Online Cambridge OPSIN
    if not resolved_smi:
      smi = self.resolve_with_opsin_online(cleaned)
      if smi:
        resolved_smi = smi
        engine_used = "OPSIN_CAMBRIDGE_ONLINE"

    # Step 3: PubChem PUG REST
    if not resolved_smi:
      res_pc = self.resolve_with_pubchem(cleaned)
      if res_pc:
        resolved_smi, inchikey = res_pc
        engine_used = "PUBCHEM_PUG_REST"

    # Step 4: NIH CIR
    if not resolved_smi:
      smi = self.resolve_with_cir(cleaned)
      if smi:
        resolved_smi = smi
        engine_used = "NIH_CIR"

    # Step 5: Fallback OCR stem repair pass if initial cascade was unsuccessful
    if not resolved_smi:
      ocr_repaired = repair_chemical_ocr_text(cleaned)
      if ocr_repaired and ocr_repaired != cleaned:
        smi = self.resolve_with_opsin_local(ocr_repaired)
        if not smi:
          smi = self.resolve_with_opsin_online(ocr_repaired)
        if not smi:
          res_pc = self.resolve_with_pubchem(ocr_repaired)
          if res_pc:
            smi, inchikey = res_pc
        if not smi:
          smi = self.resolve_with_cir(ocr_repaired)
        if smi:
          resolved_smi = smi
          engine_used = "OCR_STEM_REPAIR"
          cleaned = ocr_repaired

    if not resolved_smi:
      return None

    # Validate in RDKit
    mol = Chem.MolFromSmiles(resolved_smi)
    if mol is None:
      return None

    Chem.SanitizeMol(mol)
    mol.SetProp("_Name", cleaned)
    mol.SetProp("SourceText", cleaned)
    mol.SetProp("ResolutionEngine", engine_used or "UNKNOWN")
    mol.SetProp("CanonicalSMILES", Chem.MolToSmiles(mol))
    if inchikey:
      mol.SetProp("InChIKey", inchikey)
    else:
      try:
        mol.SetProp("InChIKey", Chem.MolToInchiKey(mol))
      except Exception:
        pass

    return mol


def main():
  parser = argparse.ArgumentParser(
      description="Multi-tier IUPAC & Chemical Name Resolver."
  )
  parser.add_argument(
      "--name", type=str, help="Single IUPAC name or chemical identifier."
  )
  parser.add_argument(
      "--input_file",
      "--input-file",
      "--input-csv",
      "--input_csv",
      "--input",
      type=str,
      dest="input_file",
      help="Path to text/CSV file with names.",
  )
  parser.add_argument(
      "--name-col",
      "--name_col",
      type=str,
      default=None,
      help="Column name for chemical names in CSV.",
  )
  parser.add_argument(
      "--output_sdf",
      "--output-sdf",
      type=str,
      dest="output_sdf",
      help="Path to write resolved molecules as SDF.",
  )
  parser.add_argument(
      "--output_csv",
      "--output-csv",
      type=str,
      dest="output_csv",
      help="Path to write resolution summary CSV.",
  )
  args = parser.parse_args()

  resolver = IUPACResolver()

  items_to_process = []
  if args.name:
    items_to_process.append({"compound_id": "CMPD_1", "name": args.name})
  elif args.input_file:
    if args.input_file.endswith(".csv"):
      df = pd.read_csv(args.input_file)
      col = args.name_col
      if not col:
        candidates = [
            c
            for c in df.columns
            if any(
                k in c.lower()
                for k in ["name", "iupac", "chemical", "compound"]
            )
        ]
        col = candidates[0] if candidates else df.columns[0]
      id_col = "compound_id" if "compound_id" in df.columns else df.columns[0]
      for _, row in df.iterrows():
        items_to_process.append(
            {"compound_id": str(row.get(id_col, "")), "name": str(row[col])}
        )
    else:
      with open(args.input_file, "r") as f:
        for i, line in enumerate(f):
          line = line.strip()
          if line and not line.startswith("#"):
            items_to_process.append(
                {"compound_id": f"CMPD_{i+1}", "name": line}
            )
  else:
    parser.error("Either --name or --input_file is required.")

  print(f"Resolving {len(items_to_process)} chemical identifiers...")
  results = []
  mols = []

  for item in items_to_process:
    name = item["name"]
    cid = item["compound_id"]
    mol = resolver.resolve(name)
    if mol is not None:
      smi = Chem.MolToSmiles(mol)
      engine = mol.GetProp("ResolutionEngine")
      print(f"✅ Resolved '{name}' via {engine}: {smi}")
      results.append({
          "compound_id": cid,
          "chemical_name": name,
          "input_name": name,
          "smiles": smi,
          "resolved_smiles": smi,
          "engine": engine,
          "inchikey": (
              mol.GetProp("InChIKey") if mol.HasProp("InChIKey") else ""
          ),
          "status": "SUCCESS",
      })
      mols.append(mol)
    else:
      print(f"❌ Could not resolve '{name}'")
      results.append({
          "compound_id": cid,
          "chemical_name": name,
          "input_name": name,
          "smiles": "",
          "resolved_smiles": "",
          "engine": "NONE",
          "inchikey": "",
          "status": "FAILED",
      })

  if args.output_sdf and mols:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_sdf)), exist_ok=True
    )
    with Chem.SDWriter(args.output_sdf) as writer:
      for m in mols:
        writer.write(m)
    print(f"Saved {len(mols)} molecules to SDF: {args.output_sdf}")

  if args.output_csv and results:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True
    )
    pd.DataFrame(results).to_csv(args.output_csv, index=False)
    print(f"Saved summary table to CSV: {args.output_csv}")


if __name__ == "__main__":
  main()
