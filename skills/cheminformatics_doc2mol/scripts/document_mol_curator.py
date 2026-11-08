# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy",
#     "pandas",
#     "rdkit",
# ]
# ///

"""Chemical Document Curation, Deduplication, and Multi-Format Exporter.

Curates extracted chemical molecules from literature, patents, Markush tables,
and synthesis procedures:
1. Canonicalizes structures and deduplicates by standard InChIKey.
2. Aggregates extraction provenance (source paragraphs, table rows, citations).
3. Computes medicinal chemistry & ADMET physicochemical properties.
4. Exports structured SDF with metadata tags, clean CSV, and interactive HTML
reports.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from typing import Any, Dict, List

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import Lipinski

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

from draw_molecules import compute_clean_2d_coords


def normalize_patent_token(token: str) -> str:
  """Normalizes PCT and national patent identifiers to canonical standard (e.g.

  WO-YYYY-NNNNNN).
  """
  if not token or not isinstance(token, str):
    return ""
  token = token.strip()
  # Match PCT application formats: PCT/US2021/045123 or PCT/EP21/099881
  pct_match = re.match(r"^PCT/([A-Z]{2})(\d{2,4})/(\d+)", token, re.IGNORECASE)
  if pct_match:
    cc, yr, num = pct_match.groups()
    yr = f"20{yr}" if len(yr) == 2 else yr
    return f"WO-{yr}-{num.zfill(6)}"

  # Match standard publication numbers with optional kind codes: WO-2021-188392-A1 or WO-2021-188392
  wo_match = re.match(
      r"^(?:WO)[-_ ]?(\d{4})[-_ ]?(\d+)(?:[-_ ]?[A-Z0-9]+)?$",
      token,
      re.IGNORECASE,
  )
  if wo_match:
    return f"WO-{wo_match.group(1)}-{wo_match.group(2).zfill(6)}"

  # Match national patents: US-11234567-B2 or US11234567
  us_match = re.match(
      r"^US[-_ ]?(\d+)(?:[-_ ]?([A-Z0-9]+))?$", token, re.IGNORECASE
  )
  if us_match:
    num = us_match.group(1)
    kind = f"-{us_match.group(2)}" if us_match.group(2) else ""
    return f"US-{num}{kind}"

  return token


def parse_stoichiometry_multiplier(formulation_str: Any) -> float:
  """Extracts exact quantitative stoichiometry multiplier from salt description."""
  if formulation_str is None or pd.isna(formulation_str):
    return 1.0
  if isinstance(formulation_str, (int, float)):
    return float(formulation_str)
  text = str(formulation_str).lower()
  if "hemi" in text or "half" in text or "0.5" in text or "semi" in text:
    return 0.5
  if "sesqui" in text or "1.5" in text:
    return 1.5
  if "bis" in text or "di" in text or "2.0" in text:
    return 2.0
  if "tri" in text or "tris" in text or "3.0" in text:
    return 3.0
  if "tetra" in text or "4.0" in text:
    return 4.0
  return 1.0


def curate_and_deduplicate(
    records: List[Dict[str, Any]],
    id_key: str = "compound_id",
    smiles_key: str = "smiles",
) -> List[Dict[str, Any]]:
  """Curate a list of compound records, computing canonical SMILES and InChIKeys,

  and merging provenance for duplicate chemical structures.
  """
  by_inchikey: Dict[str, Dict[str, Any]] = {}

  for rec in records:
    raw_smiles = rec.get(smiles_key) or rec.get("canonical_smiles") or ""
    if not raw_smiles:
      continue

    mol = Chem.MolFromSmiles(raw_smiles)
    if mol is None:
      cid = rec.get(id_key, "unknown")
      sys.stderr.write(
          f"Warning: dropping invalid SMILES {raw_smiles!r} (compound {cid})\n"
      )
      continue

    # Strip salts/fragments if compound is a salt, or keep largest fragment
    # For document curation, keep the largest organic fragment if multi-component
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
      mol = max(frags, key=lambda f: f.GetNumHeavyAtoms())

    canon_smi = Chem.MolToSmiles(mol, isomericSmiles=True)
    inchikey = Chem.MolToInchiKey(mol)

    if not inchikey or "*" in canon_smi:
      inchikey = f"DUMMY_{hashlib.sha256(canon_smi.encode()).hexdigest()[:24]}"

    if inchikey in by_inchikey:
      # Merge provenance
      existing = by_inchikey[inchikey]
      src = rec.get("source", "unknown")
      if src not in existing["sources"]:
        existing["sources"].append(src)
      prov = (
          rec.get("provenance")
          or rec.get("document_id")
          or rec.get("patent_id")
      )
      if prov and prov not in existing["provenances"]:
        existing["provenances"].append(str(prov))
      # If current record has an informative name, append alias
      name = rec.get("compound_name") or rec.get(id_key)
      if name and name not in existing["aliases"]:
        existing["aliases"].append(str(name))
      for k, v in rec.items():
        if v is None or pd.isna(v) or str(v).strip() in ("", "nan", "none"):
          continue
        if k in (
            "pIC50",
            "Normalized_pIC50",
            "IC50_nM",
            "pKi",
            "Ki_nM",
            "pActivity",
            "Activity",
            "activity_type",
            "activity_value",
            "activity_bin",
            "activity_range",
            "potency_tier",
            "patent_example_id",
        ) or str(k).startswith(("patent_assay_", "patent_table_col_")):
          existing["activity_data"].setdefault(k, v)
    else:
      # Compute 2D coordinates for visual depiction
      compute_clean_2d_coords(mol)

      # Compute medicinal chemistry properties
      mw = Descriptors.MolWt(mol)
      logp = Descriptors.MolLogP(mol)
      tpsa = Descriptors.TPSA(mol)
      hbd = Lipinski.NumHDonors(mol)
      hba = Lipinski.NumHAcceptors(mol)
      rotb = Descriptors.NumRotatableBonds(mol)
      fcsp3 = Descriptors.FractionCSP3(mol)
      ring_count = Lipinski.RingCount(mol)

      # Lipinski Rule of 5 violations
      lipinski_violations = sum([
          mw > 500,
          logp > 5.0,
          hbd > 5,
          hba > 10,
      ])

      cid = rec.get(id_key) or f"DOC_MOL_{len(by_inchikey)+1:03d}"
      cname = rec.get("compound_name") or str(cid)

      # Formulation & stoichiometry multiplier
      form_val = (
          rec.get("Formulation")
          or rec.get("formulation")
          or rec.get("formulation_type")
          or rec.get("salt")
          or rec.get("Salt")
      )
      stoich_mult = (
          parse_stoichiometry_multiplier(form_val)
          if form_val is not None
          else float(rec.get("stoichiometry_multiplier", 1.0))
      )

      # Patent identifier
      patent_val = (
          rec.get("Patent")
          or rec.get("patent")
          or rec.get("Patent_ID")
          or rec.get("patent_id")
      )
      norm_patent = (
          normalize_patent_token(str(patent_val)) if patent_val else ""
      )

      # Bioactivity endpoints (quantitative and discrete/binned) if present
      act_data = {}
      for k, v in rec.items():
        if v is None or pd.isna(v) or str(v).strip() in ("", "nan", "none"):
          continue
        if k in (
            "pIC50",
            "Normalized_pIC50",
            "IC50_nM",
            "pKi",
            "Ki_nM",
            "pActivity",
            "Activity",
            "activity_type",
            "activity_value",
            "activity_bin",
            "activity_range",
            "potency_tier",
            "patent_example_id",
        ) or str(k).startswith(("patent_assay_", "patent_table_col_")):
          act_data[k] = v

      sources = [rec.get("source", "extracted")] if rec.get("source") else []
      prov = (
          rec.get("provenance")
          or rec.get("document_id")
          or rec.get("patent_id")
      )
      provenances = [str(prov)] if prov else []

      by_inchikey[inchikey] = {
          "compound_id": str(cid),
          "compound_name": str(cname),
          "canonical_smiles": canon_smi,
          "Clean_SMILES": canon_smi,
          "stoichiometry_multiplier": float(stoich_mult),
          "Patent_ID": norm_patent,
          "activity_data": act_data,
          "inchikey": inchikey,
          "mol": mol,
          "mw": round(mw, 2),
          "logp": round(logp, 2),
          "tpsa": round(tpsa, 2),
          "hbd": int(hbd),
          "hba": int(hba),
          "rotb": int(rotb),
          "fcsp3": round(fcsp3, 2),
          "ring_count": int(ring_count),
          "lipinski_violations": int(lipinski_violations),
          "sources": sources,
          "provenances": provenances,
          "aliases": [str(cname)],
      }

  return list(by_inchikey.values())


def export_to_sdf(
    curated_records: List[Dict[str, Any]], output_sdf: str
) -> None:
  """Export curated molecules to an SD file with structured property tags."""
  os.makedirs(os.path.dirname(os.path.abspath(output_sdf)), exist_ok=True)
  with Chem.SDWriter(output_sdf) as writer:

    for rec in curated_records:
      mol = rec["mol"]
      mol.SetProp("_Name", rec["compound_id"])
      mol.SetProp("CompoundName", rec["compound_name"])
      mol.SetProp("CanonicalSMILES", rec["canonical_smiles"])
      mol.SetProp("Clean_SMILES", rec["Clean_SMILES"])
      mol.SetProp("InChIKey", rec["inchikey"])
      mol.SetProp("MolecularWeight", str(rec["mw"]))
      mol.SetProp("LogP", str(rec["logp"]))
      mol.SetProp("TPSA", str(rec["tpsa"]))
      mol.SetProp("HBD", str(rec["hbd"]))
      mol.SetProp("HBA", str(rec["hba"]))
      mol.SetProp("RotatableBonds", str(rec["rotb"]))
      mol.SetProp("LipinskiViolations", str(rec["lipinski_violations"]))
      if rec.get("stoichiometry_multiplier") is not None:
        mol.SetProp(
            "StoichiometryMultiplier", str(rec["stoichiometry_multiplier"])
        )
      if rec.get("Patent_ID"):
        mol.SetProp("PatentID", str(rec["Patent_ID"]))
      for k, v in rec.get("activity_data", {}).items():
        mol.SetProp(k, str(v))
      mol.SetProp("Sources", "; ".join(rec["sources"]))
      mol.SetProp("Provenance", "; ".join(rec["provenances"]))
      mol.SetProp("Aliases", "; ".join(rec["aliases"]))
      writer.write(mol)


def export_to_csv(
    curated_records: List[Dict[str, Any]], output_csv: str
) -> None:
  """Export curated molecules to a tabular CSV file."""
  os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
  rows = []
  for rec in curated_records:
    row = {
        "compound_id": rec["compound_id"],
        "compound_name": rec["compound_name"],
        "canonical_smiles": rec["canonical_smiles"],
        "Clean_SMILES": rec.get("Clean_SMILES", rec["canonical_smiles"]),
        "stoichiometry_multiplier": float(
            rec.get("stoichiometry_multiplier", 1.0)
        ),
        "Patent_ID": rec.get("Patent_ID", ""),
        "patent_id": rec.get("Patent_ID", ""),
        "inchikey": rec["inchikey"],
        "mw": rec["mw"],
        "logp": rec["logp"],
        "tpsa": rec["tpsa"],
        "hbd": rec["hbd"],
        "hba": rec["hba"],
        "rotb": rec["rotb"],
        "fcsp3": rec["fcsp3"],
        "ring_count": rec["ring_count"],
        "lipinski_violations": rec["lipinski_violations"],
        "sources": "; ".join(rec["sources"]),
        "provenance": "; ".join(rec["provenances"]),
        "aliases": "; ".join(rec["aliases"]),
    }
    for k, v in rec.get("activity_data", {}).items():
      if k not in row:
        row[k] = v
    rows.append(row)

  df = pd.DataFrame(rows)
  df.to_csv(output_csv, index=False)


def export_to_html_report(
    curated_records: List[Dict[str, Any]],
    output_html: str,
    title: str = "Document Extracted Chemical Library",
) -> None:
  """Generate an interactive, modern HTML report with 2D chemical structure depictions."""
  os.makedirs(os.path.dirname(os.path.abspath(output_html)), exist_ok=True)

  # Generate individual SVGs for each molecule
  card_htmls = []
  for rec in curated_records:
    mol = rec["mol"]
    svg_drawer = Draw.rdMolDraw2D.MolDraw2DSVG(260, 200)
    svg_drawer.DrawMolecule(mol)
    svg_drawer.FinishDrawing()
    raw_svg = svg_drawer.GetDrawingText()
    # Clean up XML headers in inline SVG
    clean_svg = raw_svg.replace(
        '<?xml version="1.0" encoding="iso-8859-1"?>', ""
    )

    lip_badge = (
        '<span class="badge badge-success">Ro5 Clean</span>'
        if rec["lipinski_violations"] == 0
        else (
            '<span class="badge badge-warning">Ro5 Viol:'
            f' {rec["lipinski_violations"]}</span>'
        )
    )

    act_data = rec.get("activity_data", {}) or {}
    act_badge = ""
    act_val_disp = (
        act_data.get("activity_value")
        or act_data.get("activity_bin")
        or act_data.get("pIC50")
        or act_data.get("pKi")
        or act_data.get("pActivity")
        or act_data.get("IC50_nM")
        or act_data.get("Activity")
    )
    if act_val_disp is not None:
      act_rng = act_data.get("activity_range")
      val_txt = str(act_val_disp)
      if act_rng and str(act_rng) not in val_txt:
        val_txt = f"{val_txt} ({act_rng})"
      act_badge = (
          '<span class="badge badge-success" style="margin-left:4px;">Potency:'
          f" {val_txt}</span>"
      )

    act_rows_html = ""
    if act_data:
      act_items = []
      for ak, av in act_data.items():
        if ak.endswith("_range") and ak[:-6] in act_data:
          continue
        if ak == "activity_range" and (
            "activity_value" in act_data or "activity_bin" in act_data
        ):
          continue
        disp_ak = str(ak).replace("patent_assay_", "")
        av_s = str(av)
        rng_v = act_data.get(f"{ak}_range")
        if rng_v and str(rng_v) not in av_s:
          av_s = f"{av_s} ({rng_v})"
        act_items.append(
            f"<div><small>{disp_ak}:</small> <strong>{av_s}</strong></div>"
        )
      if act_items:
        act_rows_html = (
            '<div class="prop-grid" style="background:#eef7ff;">'
            + "".join(act_items[:4])
            + "</div>"
        )

    card = f"""
        <div class="chem-card">
          <div class="card-header">
            <h3 class="compound-title">{rec['compound_id']}</h3>
            <div>{lip_badge}{act_badge}</div>
          </div>
          <div class="svg-container">
            {clean_svg}
          </div>
          <div class="card-body">
            <p class="compound-name"><strong>Name:</strong> {rec['compound_name']}</p>
            <div class="prop-grid">
              <div><small>MW:</small> <strong>{rec['mw']}</strong></div>
              <div><small>LogP:</small> <strong>{rec['logp']}</strong></div>
              <div><small>TPSA:</small> <strong>{rec['tpsa']}</strong></div>
              <div><small>HBD/A:</small> <strong>{rec['hbd']}/{rec['hba']}</strong></div>
            </div>
            {act_rows_html}
            <div class="provenance-section">
              <small class="text-muted">Sources: {', '.join(rec['sources'])}</small><br>
              <small class="text-muted">InChIKey: {rec['inchikey']}</small>
            </div>
            <div class="smiles-box">
              <code>{rec['canonical_smiles']}</code>
            </div>
          </div>
        </div>
        """
    card_htmls.append(card)

  cards_joined = "\n".join(card_htmls)

  html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>{title}</title>
  <style>
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      background-color: #f8f9fa;
      color: #212529;
      margin: 0;
      padding: 24px;
    }}
    .container {{
      max-width: 1200px;
      margin: 0 auto;
    }}
    .header {{
      background: white;
      padding: 24px;
      border-radius: 8px;
      box-shadow: 0 1px 3px rgba(0,0,0,0.1);
      margin-bottom: 24px;
    }}
    .header h1 {{ margin: 0 0 8px 0; color: #1a73e8; }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(320px, 1fr));
      gap: 20px;
    }}
    .chem-card {{
      background: white;
      border: 1px solid #e0e0e0;
      border-radius: 8px;
      overflow: hidden;
      box-shadow: 0 1px 2px rgba(0,0,0,0.05);
      transition: transform 0.15s ease-in-out;
    }}
    .chem-card:hover {{
      transform: translateY(-2px);
      box-shadow: 0 4px 10px rgba(0,0,0,0.1);
    }}
    .card-header {{
      padding: 12px 16px;
      background: #f1f3f4;
      display: flex;
      justify-content: space-between;
      align-items: center;
      border-bottom: 1px solid #e0e0e0;
    }}
    .compound-title {{ margin: 0; font-size: 16px; }}
    .svg-container {{
      display: flex;
      justify-content: center;
      padding: 12px;
      background: #ffffff;
      border-bottom: 1px solid #f1f3f4;
    }}
    .card-body {{ padding: 14px 16px; }}
    .compound-name {{ margin: 0 0 10px 0; font-size: 14px; word-break: break-word; }}
    .prop-grid {{
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 6px;
      background: #f8f9fa;
      padding: 8px 12px;
      border-radius: 4px;
      font-size: 12px;
      margin-bottom: 10px;
    }}
    .badge {{
      font-size: 11px;
      padding: 2px 8px;
      border-radius: 12px;
      font-weight: bold;
    }}
    .badge-success {{ background: #e6f4ea; color: #137333; }}
    .badge-warning {{ background: #fef7e0; color: #b06000; }}
    .provenance-section {{ font-size: 11px; margin-bottom: 8px; }}
    .text-muted {{ color: #5f6368; }}
    .smiles-box {{
      background: #202124;
      color: #8ab4f8;
      padding: 6px 10px;
      border-radius: 4px;
      font-size: 11px;
      overflow-x: auto;
      white-space: nowrap;
    }}
  </style>
</head>
<body>
  <div class="container">
    <div class="header">
      <h1>{title}</h1>
      <p>Curated Chemical Structures Extracted from Literature & Document Sources: <strong>{len(curated_records)} unique molecules</strong></p>
    </div>
    <div class="grid">
      {cards_joined}
    </div>
  </div>
</body>
</html>
"""
  with open(output_html, "w", encoding="utf-8") as f:
    f.write(html_template)


def main() -> None:
  parser = argparse.ArgumentParser(
      description="Curate, deduplicate, and export chemical document molecules."
  )
  parser.add_argument(
      "--input",
      "--input-csv",
      "--input_csv",
      "--input-file",
      "--input_file",
      dest="input",
      required=True,
      help="Input CSV or JSON file with smiles column",
  )
  parser.add_argument(
      "--smiles-col",
      "--smiles_col",
      "--smiles",
      dest="smiles_col",
      default="smiles",
      help="Name of SMILES column",
  )
  parser.add_argument(
      "--id-col",
      "--id_col",
      "--id",
      dest="id_col",
      default="compound_id",
      help="Name of compound identifier column",
  )
  parser.add_argument(
      "--output-sdf",
      "--output_sdf",
      dest="output_sdf",
      help="Path to write curated SDF file",
  )
  parser.add_argument(
      "--output-csv",
      "--output_csv",
      dest="output_csv",
      help="Path to write curated CSV file",
  )
  parser.add_argument(
      "--output-html",
      "--output_html",
      dest="output_html",
      help="Path to write interactive HTML report",
  )
  parser.add_argument(
      "--title",
      default="Extracted Document Chemical Library",
      help="Report title",
  )

  args = parser.parse_args()

  # Load input records
  if args.input.endswith(".json"):
    with open(args.input, "r", encoding="utf-8") as f:
      raw_records = json.load(f)
  elif args.input.endswith(".csv"):
    df = pd.read_csv(args.input)
    raw_records = df.to_dict(orient="records")
  else:
    print(
        f"Error: Unsupported input format {args.input}. Expected .csv or .json",
        file=sys.stderr,
    )
    sys.exit(1)

  print(f"[*] Loaded {len(raw_records)} raw records from: {args.input}")
  curated = curate_and_deduplicate(
      raw_records, id_key=args.id_col, smiles_key=args.smiles_col
  )
  print(
      f"[+] Deduplication complete: {len(curated)} unique chemical structures"
      " identified."
  )

  if args.output_csv:
    export_to_csv(curated, args.output_csv)
    print(f"[+] Exported CSV to: {args.output_csv}")

  if args.output_sdf:
    export_to_sdf(curated, args.output_sdf)
    print(f"[+] Exported SDF to: {args.output_sdf}")

  if args.output_html:
    export_to_html_report(curated, args.output_html, title=args.title)
    print(f"[+] Exported HTML visual report to: {args.output_html}")


if __name__ == "__main__":
  main()
