# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "pandas",
#     "polite-http",
# ]
# ///

"""High-Throughput Commercial Reactant Discovery, In-Silico Reactions & 2D Image Generator.

Streams commercial coupling partners from PubChem PUG-REST, executes in-silico
reactions, profiles drug-likeness (MW, cLogP, TPSA, PAINS/Brenk), resolves CAS
numbers, provides direct vendor linkouts, and renders 2D reaction schemes and
product grids.
"""

import argparse
import io
import json
import os
import re
import subprocess
import sys
import time
from typing import Dict, List, Optional
import urllib.parse

import pandas as pd
from polite_http import http_client
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, FilterCatalog, Lipinski, QED, rdMolDescriptors

_PUBCHEM_CLIENT = http_client.HttpClient(
    "https://pubchem.ncbi.nlm.nih.gov/rest/", qps=5.0
)

# Line-buffered stdout for real-time background task logging
if hasattr(sys.stdout, "reconfigure"):
  sys.stdout.reconfigure(line_buffering=True)

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

from skill_paths import add_skill_scripts, skill_dir

add_skill_scripts("cheminformatics_core")

from synthetic_accessibility_scorer import calculate_sa_score


def _resolve_grid_viewer():
  """Locates the authoritative interactive_grid_viewer.py in rdkit_general.

  Returns:
    Absolute path to the viewer script, or None if it cannot be located.
  """
  candidate = os.path.join(
      skill_dir("rdkit_general", "scripts"), "interactive_grid_viewer.py"
  )
  return candidate if os.path.exists(candidate) else None


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
  total = sum(int(d) * (len(digits) - i) for i, d in enumerate(digits))
  return total % 10 == int(p3)


def get_filter_catalog() -> FilterCatalog.FilterCatalog:
  """Builds an RDKit FilterCatalog with PAINS and Brenk alerts."""
  params = FilterCatalog.FilterCatalogParams()
  params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
  params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
  params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
  params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
  return FilterCatalog.FilterCatalog(params)


def search_pubchem_cids(
    smarts_query: str, max_cids: Optional[int] = None
) -> List[int]:
  """Queries PubChem fastsubstructure endpoint returning all matching CIDs without artificial limits."""
  clean = (
      smarts_query.replace(":1", "")
      .replace(":2", "")
      .replace(":3", "")
      .replace(":4", "")
  )
  enc = urllib.parse.quote(clean, safe="")
  cids: List[int] = []
  try:
    data = _PUBCHEM_CLIENT.fetch_json(
        f"pug/compound/fastsubstructure/smarts/{enc}/cids/JSON",
        timeout=45.0,
    )
    cids = data.get("IdentifierList", {}).get("CID", [])
  except Exception:
    pass

  if max_cids and max_cids > 0:
    cids = cids[:max_cids]
  return cids


def fetch_pubchem_properties_chunk(cids: List[int]) -> pd.DataFrame:
  """Bulk downloads properties for up to 10,000 CIDs in a single HTTP POST request."""
  if not cids:
    return pd.DataFrame()

  cid_str = ",".join(str(c) for c in cids)
  props = "CanonicalSMILES,MolecularWeight,InChIKey,IUPACName,Title,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount"
  post_data = urllib.parse.urlencode({"cid": cid_str}).encode("utf-8")
  try:
    resp = _PUBCHEM_CLIENT.fetch(
        f"pug/compound/cid/property/{props}/CSV",
        method="POST",
        data=post_data,
        timeout=45.0,
    )
    return pd.read_csv(io.BytesIO(resp.data))
  except Exception:
    return pd.DataFrame()


def fetch_pubchem_cas_chunk(
    cids: List[int], batch_size: int = 100
) -> Dict[int, str]:
  """Bulk fetches verified CAS numbers for CIDs via PUG-REST synonyms."""
  cas_map: Dict[int, str] = {}
  if not cids:
    return cas_map

  for i in range(0, len(cids), batch_size):
    sub_cids = cids[i : i + batch_size]
    cid_str = ",".join(str(c) for c in sub_cids)
    post_data = urllib.parse.urlencode({"cid": cid_str}).encode("utf-8")
    try:
      data = _PUBCHEM_CLIENT.fetch_json(
          "pug/compound/cid/synonyms/JSON",
          method="POST",
          data=post_data,
          timeout=15.0,
      )
      for info in data.get("InformationList", {}).get("Information", []):
        cid = info.get("CID")
        for s in info.get("Synonym", []):
          if is_valid_cas(s):
            cas_map[cid] = s
            break
    except Exception:
      pass
  return cas_map


def dataframe_to_markdown_table(df: pd.DataFrame) -> str:
  """Renders a pandas DataFrame as a GitHub-flavored Markdown table."""
  if df.empty:
    return ""
  cols = list(df.columns)
  lines = [
      "| " + " | ".join(str(c) for c in cols) + " |",
      "| " + " | ".join("---" for _ in cols) + " |",
  ]
  for _, row in df.iterrows():
    lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
  return "\n".join(lines)


def discover_and_enumerate_coupling(
    core_smiles: str,
    reaction_smarts: str,
    target_partner_smarts: Optional[str] = None,
    max_partners: int = 0,
    batch_size: int = 5000,
    output_csv: Optional[str] = None,
    output_json: Optional[str] = None,
) -> pd.DataFrame:
  """High-throughput scalable pipeline: Fastsubstructure CID query -> Chunked POST -> RDKit Stream -> CSV/JSON."""
  core_mol = Chem.MolFromSmiles(core_smiles)
  if core_mol is None:
    raise ValueError(f"Invalid core SMILES: {core_smiles}")

  rxn = AllChem.ReactionFromSmarts(reaction_smarts)
  if rxn is None:
    raise ValueError(f"Invalid Reaction SMARTS: {reaction_smarts}")

  if not target_partner_smarts:
    if rxn.GetNumReactantTemplates() < 2:
      raise ValueError("Reaction requires at least 2 reactant templates.")
    tmpl = rxn.GetReactantTemplate(1)
    clean_tmpl = Chem.Mol(tmpl)
    for a in clean_tmpl.GetAtoms():
      a.SetAtomMapNum(0)
    search_query = Chem.MolToSmarts(clean_tmpl)
    # PubChem fastsubstructure endpoint does not support comma lists like [Br,I,Cl]
    if "," in search_query and "[" in search_query:
      search_query = re.sub(r"\[([A-Za-z]+),[^\]]*\]", r"[\1]", search_query)
  else:
    search_query = target_partner_smarts

  catalog = get_filter_catalog()

  is_unlimited = max_partners is None or max_partners <= 0
  print(
      "1. Querying PubChem Database for commercial partners matching:"
      f" '{search_query}'...",
      flush=True,
  )
  all_cids = search_pubchem_cids(
      search_query, max_cids=None if is_unlimited else max_partners
  )

  if not all_cids:
    print("   No matching compound CIDs found in PubChem.", flush=True)
    return pd.DataFrame()

  total_cids = len(all_cids)
  print(
      f"   Successfully matched {total_cids:,} compound CIDs in PubChem"
      " database.",
      flush=True,
  )
  if is_unlimited:
    print(
        "   Exhaustive mode: Retrieving, reacting, and profiling ALL"
        f" {total_cids:,} compounds...",
        flush=True,
    )
  else:
    print(
        "   Retrieving, reacting, and profiling top"
        f" {min(max_partners, total_cids):,} compounds...",
        flush=True,
    )

  target_cids = all_cids if is_unlimited else all_cids[:max_partners]
  chunks = [
      target_cids[i : i + batch_size]
      for i in range(0, len(target_cids), batch_size)
  ]
  total_processed = 0
  first_chunk = True
  display_records = []
  all_records = []

  if output_csv:
    out_dir = os.path.dirname(os.path.abspath(output_csv))
    if out_dir:
      os.makedirs(out_dir, exist_ok=True)

  print(
      f"2. Processing {len(chunks)} batch(es) of properties and in-silico"
      " reactions...",
      flush=True,
  )
  start_time = time.time()

  for idx, chunk in enumerate(chunks, 1):
    chunk_start = time.time()
    df_prop = fetch_pubchem_properties_chunk(chunk)
    if df_prop.empty:
      continue

    # Bulk fetch CAS records
    cas_map = fetch_pubchem_cas_chunk(chunk)

    records = []
    for _, r in df_prop.iterrows():
      smi = r.get("ConnectivitySMILES") or r.get("CanonicalSMILES")
      if not isinstance(smi, str):
        continue
      mol = Chem.MolFromSmiles(smi)
      if not mol:
        continue

      outcomes = rxn.RunReactants((core_mol, mol))
      if not outcomes or not outcomes[0]:
        continue
      prod_mol = outcomes[0][0]
      try:
        Chem.SanitizeMol(prod_mol)
      except Exception:
        continue

      # Ensure no residual boron remains in the product (monofunctional coupling completion)
      if any(atom.GetAtomicNum() == 5 for atom in prod_mol.GetAtoms()):
        continue

      cid = int(r["CID"])
      name = str(r.get("Title") or r.get("IUPACName") or f"CID_{cid}")
      prod_smi = Chem.MolToSmiles(prod_mol)
      mw = round(Descriptors.ExactMolWt(prod_mol), 1)
      clogp = round(Crippen.MolLogP(prod_mol), 2)
      tpsa = round(rdMolDescriptors.CalcTPSA(prod_mol), 1)
      hbd = Lipinski.NumHDonors(prod_mol)
      hba = Lipinski.NumHAcceptors(prod_mol)
      rotb = Lipinski.NumRotatableBonds(prod_mol)
      fsp3 = round(rdMolDescriptors.CalcFractionCSP3(prod_mol), 2)
      qed_val = round(QED.qed(prod_mol), 2)
      sa_val = round(calculate_sa_score(prod_mol), 2)

      matches = catalog.GetMatches(prod_mol)
      pains_list = [
          m.GetDescription() for m in matches if "PAINS" in m.GetDescription()
      ]
      brenk_list = [
          m.GetDescription() for m in matches if "Brenk" in m.GetDescription()
      ]

      ro5_violations = sum([
          1 if mw > 500.0 else 0,
          1 if clogp > 5.0 else 0,
          1 if hbd > 5 else 0,
          1 if hba > 10 else 0,
      ])

      rec = {
          "Reactant_CID": cid,
          "Reactant_Name": name,
          "Reactant_SMILES": smi,
          "CAS": cas_map.get(cid, None),
          "Where_to_Purchase": (
              f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=Chemical-Vendors"
          ),
          "Product_SMILES": prod_smi,
          "Product_MW": mw,
          "Product_cLogP": clogp,
          "Product_TPSA": tpsa,
          "Product_HBD": hbd,
          "Product_HBA": hba,
          "Product_RotB": rotb,
          "Product_Fsp3": fsp3,
          "Product_QED": qed_val,
          "Product_SA_Score": sa_val,
          "PAINS_Alerts": "; ".join(pains_list) if pains_list else "Zero PAINS",
          "Brenk_Alerts": "; ".join(brenk_list) if brenk_list else "Zero Brenk",
          "Ro5_Violations": ro5_violations,
      }
      records.append(rec)

    if not records:
      continue

    all_records.extend(records)
    chunk_df = pd.DataFrame(records)
    total_processed += len(chunk_df)

    if output_csv:
      chunk_df.to_csv(
          output_csv,
          mode="w" if first_chunk else "a",
          header=first_chunk,
          index=False,
      )
      first_chunk = False

    if len(display_records) < 25:
      needed = 25 - len(display_records)
      display_records.extend(records[:needed])

    chunk_elapsed = time.time() - chunk_start
    print(
        f"   [Batch {idx}/{len(chunks)}] Processed {len(chunk_df):,} products"
        f" ({chunk_elapsed:.2f}s) | Total streamed: {total_processed:,}",
        flush=True,
    )

  elapsed = time.time() - start_time
  rate = total_processed / elapsed if elapsed > 0 else 0
  print(
      f"\nCompleted reaction transformation: {total_processed:,} products in"
      f" {elapsed:.2f}s ({rate:,.1f} rxns/sec).",
      flush=True,
  )

  if output_json and all_records:
    out_jdir = os.path.dirname(os.path.abspath(output_json))
    if out_jdir:
      os.makedirs(out_jdir, exist_ok=True)
    with open(output_json, "w", encoding="utf-8") as f:
      json.dump(all_records, f, indent=2)
    print(
        f"JSON dataset ({len(all_records):,} compounds) saved to:"
        f" {os.path.abspath(output_json)}",
        flush=True,
    )

  display_df = pd.DataFrame(display_records)
  return display_df


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Find real commercial coupling partners in PubChem at scale and"
          " generate 2D depictions."
      )
  )
  parser.add_argument(
      "--core", type=str, required=True, help="Core starting material SMILES."
  )
  parser.add_argument(
      "--reaction_smarts",
      type=str,
      required=True,
      help="Reaction SMARTS string.",
  )
  parser.add_argument(
      "--partner_smarts",
      type=str,
      default=None,
      help="Explicit partner query SMARTS (optional).",
  )
  parser.add_argument(
      "--batch_size",
      type=int,
      default=5000,
      help="Batch size for PubChem property POST requests.",
  )
  parser.add_argument(
      "--all",
      action="store_true",
      default=False,
      help="Exhaustive unconstrained retrieval.",
  )
  parser.add_argument(
      "--max_partners",
      type=int,
      default=0,
      help="Maximum number of partners to retrieve (0 for all).",
  )
  parser.add_argument(
      "--output", type=str, default=None, help="Output file path (CSV or JSON)."
  )
  parser.add_argument(
      "--output_json", type=str, default=None, help="Save output to JSON file."
  )
  parser.add_argument(
      "--output_csv", type=str, default=None, help="Save output to CSV file."
  )
  parser.add_argument(
      "--output_html",
      type=str,
      default=None,
      help="Generate interactive 2D HTML chemical structure grid.",
  )
  parser.add_argument(
      "--output_scheme_png",
      type=str,
      default=None,
      help="Generate 2D reaction scheme PNG image.",
  )
  parser.add_argument(
      "--output_products_png",
      type=str,
      default=None,
      help="Generate 2D chemical structure grid PNG of top products.",
  )
  args = parser.parse_args()

  out_csv = args.output_csv
  out_json = args.output_json
  if args.output:
    if args.output.endswith(".json"):
      out_json = args.output
    else:
      out_csv = args.output

  max_p = args.max_partners
  if args.all:
    max_p = 0
  elif max_p == 0 and not args.all:
    # Default to 500 partners for responsive retrieval
    max_p = 500

  df = discover_and_enumerate_coupling(
      core_smiles=args.core,
      reaction_smarts=args.reaction_smarts,
      target_partner_smarts=args.partner_smarts,
      max_partners=max_p,
      batch_size=args.batch_size,
      output_csv=out_csv,
      output_json=out_json,
  )

  # 2D Reaction Scheme Generation
  scheme_png = args.output_scheme_png
  if not scheme_png and args.output_csv:
    scheme_png = os.path.splitext(args.output_csv)[0] + "_reaction_scheme.png"
  if scheme_png:
    try:
      draw_script = os.path.join(
          skill_dir("rdkit_general", "scripts"), "draw_molecules.py"
      )
      if os.path.exists(draw_script):
        cmd = [
            sys.executable,
            draw_script,
            "--reaction_smarts",
            args.reaction_smarts,
            "--output",
            os.path.abspath(scheme_png),
        ]
        subprocess.run(cmd, check=True)
        print(
            "\nRendered 2D reaction scheme image:"
            f" {os.path.abspath(scheme_png)}",
            flush=True,
        )
        print(
            "  Markdown embed: ![Reaction"
            f" Scheme]({os.path.abspath(scheme_png)})",
            flush=True,
        )
    except Exception as e:
      print(
          f"Warning: Could not render reaction scheme PNG: {e}",
          flush=True,
          file=sys.stderr,
      )

  # 2D Top Products Grid Generation
  products_png = args.output_products_png
  if not products_png and args.output_csv:
    products_png = os.path.splitext(args.output_csv)[0] + "_top_products.png"
  if products_png and not df.empty:
    try:
      draw_script = os.path.join(
          skill_dir("rdkit_general", "scripts"), "draw_molecules.py"
      )
      if os.path.exists(draw_script):
        top_smiles = df["Product_SMILES"].head(12).tolist()
        top_names = df["Reactant_Name"].head(12).tolist()
        cmd = [
            sys.executable,
            draw_script,
            "--smiles",
            *top_smiles,
            "--legends",
            *top_names,
            "--align_core",
            args.core,
            "--output",
            os.path.abspath(products_png),
        ]
        subprocess.run(cmd, check=True)
        print(
            "Rendered 2D chemical structure grid of top products:"
            f" {os.path.abspath(products_png)}",
            flush=True,
        )
        print(
            "  Markdown embed: ![Top Products"
            f" Grid]({os.path.abspath(products_png)})",
            flush=True,
        )
    except Exception as e:
      print(
          f"Warning: Could not render products grid PNG: {e}",
          flush=True,
          file=sys.stderr,
      )

  print("\n" + "=" * 80, flush=True)
  print(
      "  REPRESENTATIVE REAL COMMERCIAL COUPLING PARTNERS & REACTION PRODUCTS",
      flush=True,
  )
  print("=" * 80, flush=True)
  if df.empty:
    print(
        "No commercial coupling partners found matching reaction criteria.",
        flush=True,
    )
  else:
    show_cols = [
        "Reactant_Name",
        "Reactant_SMILES",
        "CAS",
        "Where_to_Purchase",
        "Product_MW",
        "Product_cLogP",
        "Product_TPSA",
        "Product_QED",
        "PAINS_Alerts",
    ]
    available_cols = [c for c in show_cols if c in df.columns]
    print(dataframe_to_markdown_table(df[available_cols].head(10)), flush=True)

  if args.output_csv:
    print(
        f"\nFull dataset saved to: {os.path.abspath(args.output_csv)}",
        flush=True,
    )

  if args.output_html and args.output_csv:
    try:
      viewer_script = _resolve_grid_viewer()
      if viewer_script:
        cmd = [
            sys.executable,
            viewer_script,
            "--input",
            os.path.abspath(args.output_csv),
            "--output",
            os.path.abspath(args.output_html),
            "--align_core",
            args.core,
        ]
        subprocess.run(cmd, check=True)
        print(
            "Interactive 2D HTML grid saved to:"
            f" {os.path.abspath(args.output_html)}",
            flush=True,
        )
    except Exception as e:
      print(
          f"Warning: Could not generate HTML grid: {e}",
          flush=True,
          file=sys.stderr,
      )


if __name__ == "__main__":
  main()
