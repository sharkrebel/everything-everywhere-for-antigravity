# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "rdkit",
#   "pandas",
#   "numpy",
#   "polite-http",
# ]
# ///

"""Academic Paper & Bioassay Literature Harvester (paper_fetcher.py).

Dedicated CLI and library engine for harvesting small-molecule structures,
paper compound keys (`paper_compound_id`, e.g. '1a', '2b', 'Compound-101'),
and wide pivoted bioactivity endpoints (`Ki_<target>`, `Ki_<target>_<mutant>`,
`IC50_<target>`, `pchembl_value`) from academic publications by:
  - PubMed IDs (`--pmids <PMID_1> <PMID_2>` or `PMID:<PMID>`)
  - ChEMBL Document IDs (`--chembl-docs CHEMBL<DOC_ID>`)
  - DOIs (`--dois 10.1021/<DOI_SUFFIX>`)

Uses rate-limited retries via `polite_http.http_client.HttpClient` and automatic
fallback to NCBI E-utilities (`elink.fcgi?dbfrom=pubmed&db=pcassay` +
`db=pccompound`) and PubChem BioAssay concise CSVs whenever EMBL-EBI ChEMBL is
unavailable or does not index a paper.
"""

from __future__ import annotations

import argparse
import io
import json
import os
import re
import sys
from typing import Any, Optional, Sequence
import urllib.parse

import numpy as np
import pandas as pd
from polite_http import http_client
from rdkit import Chem

_cur_dir = os.path.dirname(os.path.abspath(__file__))
if _cur_dir not in sys.path:
  sys.path.insert(0, _cur_dir)

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
from iupac_resolver import extract_pubchem_smiles

CHEMBL_BASE_URL = "https://www.ebi.ac.uk"
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov"
EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov"
_DEFAULT_TIMEOUT_SEC = 12.0

_MUTATION_TOKEN_RE = re.compile(r"\b([A-Z]\d{3,4}[A-Z])\b")
# Matches single-letter amino-acid point-mutation codes (e.g. t790m, l858r, g12c)
# in lowercased assay descriptions to distinguish mutant vs wild-type target assays.
_AMINO_ACID_POINT_MUTATION_RE = re.compile(
    r"\b[acdefghiklmnpqrstvwy]\d{2,4}[acdefghiklmnpqrstvwy]\b"
)

_CHEMBL_CLIENT = http_client.HttpClient(
    CHEMBL_BASE_URL, qps=5.0, timeout=_DEFAULT_TIMEOUT_SEC, max_retries=3
)
_PUBCHEM_CLIENT = http_client.HttpClient(
    PUBCHEM_BASE_URL, qps=5.0, timeout=_DEFAULT_TIMEOUT_SEC, max_retries=3
)
_EUTILS_CLIENT = http_client.HttpClient(
    EUTILS_BASE_URL, qps=3.0, timeout=_DEFAULT_TIMEOUT_SEC, max_retries=3
)


def _build_endpoint_column_key(
    std_type: object, target_name: object, assay_desc: object = None
) -> Optional[str]:
  """Builds a normalized wide-table column name for a (standard_type, target, mutant) bioassay endpoint."""
  if not std_type or pd.isna(std_type):
    return None
  st = re.sub(r"[^A-Za-z0-9]+", "_", str(std_type).strip()).strip("_")
  tgt = (
      re.sub(r"[^A-Za-z0-9]+", "_", str(target_name).strip()).strip("_")
      if target_name and pd.notna(target_name)
      else "assay"
  )
  if not st:
    return None
  col_key = f"{st}_{tgt}"
  if assay_desc and pd.notna(assay_desc):
    muts = _MUTATION_TOKEN_RE.findall(str(assay_desc).upper())
    if muts:
      uniq_muts = "_".join(dict.fromkeys(muts))
      if uniq_muts not in col_key:
        col_key = f"{col_key}_{uniq_muts}"
  return col_key


def _clean_calendar_year(raw_val: Any) -> Optional[int]:
  """Extracts a valid calendar year in [1800, 2100], rejecting booleans."""
  if raw_val is None or isinstance(raw_val, (bool, np.bool_)):
    return None
  try:
    if pd.isna(raw_val):
      return None
  except (TypeError, ValueError):
    pass
  try:
    fval = float(raw_val)
    if not np.isfinite(fval):
      return None
    yr = int(round(fval))
    return yr if 1800 <= yr <= 2100 else None
  except (TypeError, ValueError):
    m_yr = re.match(r"^\s*(18\d{2}|19\d{2}|20\d{2})\b", str(raw_val))
    return int(m_yr.group(1)) if m_yr else None


def _fetch_pubmed_publication_info(
    pmid_clean: str,
) -> tuple[Optional[int], Optional[str]]:
  """Resolves publication year and ISO date from NCBI PubMed esummary."""
  if not pmid_clean or not pmid_clean.isdigit():
    return None, None
  esum_url = (
      f"{EUTILS_BASE_URL}/entrez/eutils/esummary.fcgi"
      f"?db=pubmed&id={pmid_clean}&retmode=json"
  )
  try:
    esum_json = _EUTILS_CLIENT.fetch_json(esum_url) or {}
    entry = (esum_json.get("result") or {}).get(pmid_clean) or {}
    pub_date: Optional[str] = None
    doc_year: Optional[int] = None
    for date_key in ("sortpubdate", "epubdate", "pubdate"):
      raw_dt = str(entry.get(date_key) or "").strip()
      if not raw_dt:
        continue
      m_iso = re.search(r"\b(19\d{2}|20\d{2})[-/](\d{2})[-/](\d{2})\b", raw_dt)
      if m_iso and pub_date is None:
        pub_date = f"{m_iso.group(1)}-{m_iso.group(2)}-{m_iso.group(3)}"
      if doc_year is None:
        m_yr = re.search(r"\b(18\d{2}|19\d{2}|20\d{2})\b", raw_dt)
        if m_yr:
          doc_year = int(m_yr.group(1))
    return doc_year, pub_date
  except (
      http_client.HttpError,
      ValueError,
      KeyError,
      TypeError,
      OSError,
  ):
    return None, None


def _fetch_pubchem_pmid_fallback(pmid: str) -> list[dict[str, Any]]:
  """Fallback harvester via NCBI E-utilities (elink) and PubChem PUG REST when EMBL-EBI ChEMBL is unavailable for a PMID."""
  records: list[dict[str, Any]] = []
  pmid_clean = re.sub(r"(?i)^pmid:\s*", "", str(pmid).strip())
  if not pmid_clean.isdigit():
    return records

  print(
      "  [+] Querying NCBI E-utilities (elink) & PubChem BioAssay fallback for"
      f" PMID:{pmid_clean}...",
      flush=True,
  )
  cids: list[int] = []
  try:
    elink_c_url = (
        f"{EUTILS_BASE_URL}/entrez/eutils/elink.fcgi"
        f"?dbfrom=pubmed&db=pccompound&id={pmid_clean}&retmode=json"
    )
    elink_c_json = _EUTILS_CLIENT.fetch_json(elink_c_url)
    for ls in (elink_c_json or {}).get("linksets", []):
      for ldb in ls.get("linksetdbs", []):
        for lk in ldb.get("links", []):
          c_val = int(lk["id"] if isinstance(lk, dict) else lk)
          if c_val not in cids:
            cids.append(c_val)
  except Exception:
    pass

  if not cids:
    try:
      cid_url = f"{PUBCHEM_BASE_URL}/rest/pug/compound/xref/PubMedID/{pmid_clean}/cids/JSON"
      cid_json = _PUBCHEM_CLIENT.fetch_json(cid_url)
      if isinstance(cid_json, dict):
        for c_val in cid_json.get("IdentifierList", {}).get("CID", []):
          if int(c_val) not in cids:
            cids.append(int(c_val))
    except Exception:
      pass

  aids: list[int] = []
  try:
    elink_a_url = (
        f"{EUTILS_BASE_URL}/entrez/eutils/elink.fcgi"
        f"?dbfrom=pubmed&db=pcassay&id={pmid_clean}&retmode=json"
    )
    elink_a_json = _EUTILS_CLIENT.fetch_json(elink_a_url)
    for ls in (elink_a_json or {}).get("linksets", []):
      for ldb in ls.get("linksetdbs", []):
        if ldb.get("linkname") in (None, "pubmed_pcassay"):
          for lk in ldb.get("links", []):
            a_val = int(lk["id"] if isinstance(lk, dict) else lk)
            if a_val not in aids:
              aids.append(a_val)
  except Exception:
    pass

  aids = sorted(aids)
  assay_rows_by_cid: dict[int, list[dict[str, Any]]] = {}
  for batch_start in range(0, min(len(aids), 100), 25):
    batch_aids = ",".join(str(a) for a in aids[batch_start : batch_start + 25])
    try:
      csv_url = (
          f"{PUBCHEM_BASE_URL}/rest/pug/assay/aid/{batch_aids}/concise/CSV"
      )
      csv_text = _PUBCHEM_CLIENT.fetch_text(csv_url)
      if not csv_text or not csv_text.strip():
        continue
      adf = pd.read_csv(io.StringIO(csv_text))
      if "PubMed ID" in adf.columns:
        has_pmid_match = (
            adf["PubMed ID"].astype(str).str.contains(pmid_clean, na=False)
        )
        if has_pmid_match.any():
          adf = adf[has_pmid_match]
      for _, row in adf.iterrows():
        cid_val = row.get("CID")
        if pd.isna(cid_val):
          continue
        cid_int = int(float(cid_val))
        if cid_int not in cids:
          cids.append(cid_int)
        raw_act_name = row.get("Activity Name")
        act_type = (
            str(raw_act_name).strip()
            if pd.notna(raw_act_name) and str(raw_act_name).strip()
            else "IC50"
        )
        act_uM = row.get("Activity Value [uM]")
        aid_cur = row.get("AID", batch_aids)
        target_nm = str(
            row.get("Target Name")
            or row.get("Target Accession")
            or f"AID_{aid_cur}"
        ).strip()
        raw_assay_name = row.get("Assay Name")
        assay_desc = (
            str(raw_assay_name).strip()
            if pd.notna(raw_assay_name) and str(raw_assay_name).strip()
            else f"PubChem AID {aid_cur} (PMID {pmid_clean})"
        )
        std_val = None
        std_units = None
        pchembl_val = None
        if pd.notna(act_uM):
          try:
            uM_f = float(act_uM)
            if uM_f > 0:
              if act_type.upper() in ("KI", "IC50", "KD", "EC50", "KIS"):
                std_val = round(uM_f * 1000.0, 4)
                std_units = "nM"
                pchembl_val = round(-np.log10(uM_f * 1e-6), 3)
              else:
                std_val = round(uM_f, 4)
          except Exception:
            pass
        assay_rows_by_cid.setdefault(cid_int, []).append({
            "standard_type": act_type,
            "standard_value": std_val,
            "standard_units": std_units,
            "pchembl_value": pchembl_val,
            "target_name": target_nm,
            "assay_description": assay_desc,
        })
    except Exception:
      continue

  if not cids:
    return records

  doc_year, pub_date = _fetch_pubmed_publication_info(pmid_clean)

  # Fetch structures and titles in batches of 100 CIDs
  cid_props: dict[int, dict[str, Any]] = {}
  for idx in range(0, min(len(cids), 500), 100):
    batch = [str(c) for c in cids[idx : idx + 100]]
    prop_url = (
        f"{PUBCHEM_BASE_URL}/rest/pug/compound/cid/{','.join(batch)}"
        "/property/SMILES,ConnectivitySMILES,Title,IUPACName/JSON"
    )
    try:
      p_json = _PUBCHEM_CLIENT.fetch_json(prop_url)
      if isinstance(p_json, dict):
        plist = p_json.get("PropertyTable", {}).get("Properties", [])
        for p in plist:
          c_id = int(p.get("CID"))
          cid_props[c_id] = p
    except Exception:
      continue

  for cid_int, p in cid_props.items():
    smi = extract_pubchem_smiles(p)
    if not smi:
      continue
    title = str(p.get("Title") or f"CID_{cid_int}").strip()
    iupac = str(p.get("IUPACName") or "").strip()
    a_list = assay_rows_by_cid.get(cid_int) or [{
        "standard_type": None,
        "standard_value": None,
        "standard_units": None,
        "pchembl_value": None,
        "target_name": None,
        "assay_description": f"PubChem Xref PMID:{pmid_clean}",
    }]
    for a_item in a_list:
      rec: dict[str, Any] = {
          "source": f"PubChem_PMID_{pmid_clean}",
          "compound_id": f"CID{cid_int}",
          "paper_compound_id": title if len(title) <= 20 else f"CID{cid_int}",
          "compound_name": title[:28] if title else f"CID{cid_int}",
          "iupac_name": iupac,
          "canonical_smiles": smi,
          "target_name": a_item["target_name"],
          "target_chembl_id": None,
          "standard_type": a_item["standard_type"],
          "standard_relation": "=",
          "standard_value": a_item["standard_value"],
          "standard_units": a_item["standard_units"],
          "pchembl_value": a_item["pchembl_value"],
          "assay_description": a_item["assay_description"],
          "document_chembl_id": f"PMID:{pmid_clean}",
          "document_year": doc_year,
          "publication_date": pub_date,
      }
      col_key = _build_endpoint_column_key(
          a_item["standard_type"],
          a_item["target_name"],
          a_item["assay_description"],
      )
      if col_key and a_item["standard_value"] is not None:
        rec[col_key] = a_item["standard_value"]
        rec["all_endpoints"] = json.dumps({col_key: a_item["standard_value"]})
      records.append(rec)
  print(
      f"  -> Retrieved {len(records)} PubChem fallback records across"
      f" {len(cid_props)} compounds for PMID:{pmid_clean}.",
      flush=True,
  )
  return records


def fetch_chembl_document_activities(
    doc_chembl_ids: Sequence[str],
    max_records: int = 2000,
) -> pd.DataFrame:
  """Fetches bioactivity records, compound keys, and chemical structures from ChEMBL with retry/backoff and PubChem PMID fallback."""
  all_records: list[dict[str, Any]] = []

  expanded_ids: list[str] = []
  for raw_doc in doc_chembl_ids:
    for part in str(raw_doc).split(","):
      if part.strip():
        expanded_ids.append(part.strip())

  for doc_idx, raw_id in enumerate(expanded_ids):
    token = raw_id.strip()
    pmid_token: Optional[str] = None
    resolved_docs: list[str] = []
    doc_years: dict[str, int] = {}

    if re.match(r"(?i)^pmid:\s*\d+$", token) or token.isdigit():
      pmid_token = re.sub(r"(?i)^pmid:\s*", "", token)
      doc_lookup_url = (
          f"{CHEMBL_BASE_URL}/chembl/api/data/document.json"
          f"?pubmed_id={pmid_token}&limit=20"
      )
      try:
        doc_json = _CHEMBL_CLIENT.fetch_json(doc_lookup_url)
        for d_entry in (doc_json or {}).get("documents") or []:
          dc_id = d_entry.get("document_chembl_id")
          if dc_id:
            dc_norm = str(dc_id).strip().upper()
            resolved_docs.append(dc_norm)
            yr_val = _clean_calendar_year(d_entry.get("year"))
            if yr_val is not None:
              doc_years[dc_norm] = yr_val
      except Exception as e:
        print(
            f"      [!] Could not resolve PMID:{pmid_token} on ChEMBL ({e});"
            " activating PubChem fallback.",
            flush=True,
        )
      if not resolved_docs:
        for fb_rec in _fetch_pubchem_pmid_fallback(pmid_token):
          fb_rec["_doc_order"] = doc_idx
          all_records.append(fb_rec)
        continue
    elif token.startswith("10.") and "/" in token:
      doc_lookup_url = (
          f"{CHEMBL_BASE_URL}/chembl/api/data/document.json"
          f"?doi={token}&limit=20"
      )
      try:
        doc_json = _CHEMBL_CLIENT.fetch_json(doc_lookup_url)
        for d_entry in (doc_json or {}).get("documents") or []:
          dc_id = d_entry.get("document_chembl_id")
          if dc_id:
            dc_norm = str(dc_id).strip().upper()
            resolved_docs.append(dc_norm)
            yr_val = _clean_calendar_year(d_entry.get("year"))
            if yr_val is not None:
              doc_years[dc_norm] = yr_val
      except Exception as e:
        print(
            f"      [!] Could not resolve DOI {token} on ChEMBL ({e});"
            " resolving via NCBI PubMed...",
            flush=True,
        )
      if not resolved_docs:
        try:
          es_url = (
              f"{EUTILS_BASE_URL}/entrez/eutils/esearch.fcgi"
              f"?db=pubmed&term={urllib.parse.quote(token + '[aid]')}&retmode=json"
          )
          es_json = _EUTILS_CLIENT.fetch_json(es_url)
          id_list = (es_json or {}).get("esearchresult", {}).get("idlist", [])
          if id_list:
            pmid_token = str(id_list[0]).strip()
            for fb_rec in _fetch_pubchem_pmid_fallback(pmid_token):
              fb_rec["_doc_order"] = doc_idx
              all_records.append(fb_rec)
        except Exception:
          pass
        continue
    else:
      resolved_docs.append(token.upper())

    for doc_clean in resolved_docs:
      print(
          f"  [+] Harvesting ChEMBL document activities for {doc_clean}...",
          flush=True,
      )
      if doc_clean not in doc_years:
        try:
          d_meta = _CHEMBL_CLIENT.fetch_json(
              f"{CHEMBL_BASE_URL}/chembl/api/data/document/{doc_clean}.json"
          )
          yr_meta = (
              _clean_calendar_year(d_meta.get("year"))
              if isinstance(d_meta, dict)
              else None
          )
          if yr_meta is not None:
            doc_years[doc_clean] = yr_meta
          elif isinstance(d_meta, dict) and d_meta.get("pubmed_id"):
            fb_yr, _ = _fetch_pubmed_publication_info(str(d_meta["pubmed_id"]))
            if fb_yr is not None:
              doc_years[doc_clean] = fb_yr
        except (
            http_client.HttpError,
            ValueError,
            KeyError,
            TypeError,
            OSError,
        ):
          pass
      doc_records_before = len(all_records)
      # 1. Fetch compound_record.json to preserve paper compound labels (e.g. '1a', '2b')
      mol_to_key: dict[str, str] = {}
      mol_to_name: dict[str, str] = {}
      cr_offset = 0
      while cr_offset < max_records:
        cr_url = (
            f"{CHEMBL_BASE_URL}/chembl/api/data/compound_record.json"
            f"?document_chembl_id={doc_clean}&limit=1000&offset={cr_offset}"
        )
        try:
          cr_data = _CHEMBL_CLIENT.fetch_json(cr_url)
          c_recs = (cr_data or {}).get("compound_records") or []
          if not c_recs:
            break
          for cr in c_recs:
            mid = cr.get("molecule_chembl_id")
            if not mid:
              continue
            ckey = cr.get("compound_key")
            cname = cr.get("compound_name")
            if ckey and mid not in mol_to_key:
              mol_to_key[mid] = str(ckey).strip()
            if cname and mid not in mol_to_name:
              mol_to_name[mid] = str(cname).strip()
          cr_offset += len(c_recs)
          if len(c_recs) < 1000:
            break
        except Exception:
          break

      # 2. Fetch activity.json with pagination
      offset = 0
      limit = 100
      while offset < max_records:
        url = (
            f"{CHEMBL_BASE_URL}/chembl/api/data/activity.json"
            f"?document_chembl_id={doc_clean}&limit={limit}&offset={offset}"
        )
        try:
          data = _CHEMBL_CLIENT.fetch_json(url) or {}
          acts = data.get("activities", [])
          if not acts:
            break
          for a in acts:
            mid = a.get("molecule_chembl_id")
            std_type = a.get("standard_type")
            raw_std_val = a.get("standard_value")
            try:
              std_val = float(raw_std_val) if raw_std_val is not None else None
            except (TypeError, ValueError):
              std_val = raw_std_val
            raw_pchembl = a.get("pchembl_value")
            try:
              pchembl_val = (
                  float(raw_pchembl) if raw_pchembl is not None else None
              )
            except (TypeError, ValueError):
              pchembl_val = raw_pchembl
            target_name = a.get("target_pref_name")
            assay_desc = a.get("assay_description")
            raw_cname = mol_to_name.get(mid, "")
            raw_ckey = mol_to_key.get(mid, "")
            if raw_cname and len(raw_cname) <= 28:
              short_cname = raw_cname
            elif raw_ckey and len(raw_ckey) <= 20:
              is_code = bool(
                  re.match(
                      r"^(?:cmpd|compound|example|ex\b|chembl|[a-z]{2,}[-_]?\d)",
                      raw_ckey,
                      re.IGNORECASE,
                  )
              )
              short_cname = raw_ckey if is_code else f"Cmpd {raw_ckey}"
            else:
              short_cname = mid or raw_cname[:24]
            act_yr = _clean_calendar_year(a.get("document_year"))
            if act_yr is None:
              act_yr = doc_years.get(doc_clean)
            rec: dict[str, Any] = {
                "source": f"ChEMBL_{doc_clean}",
                "compound_id": mid,
                "paper_compound_id": raw_ckey,
                "compound_name": short_cname,
                "iupac_name": raw_cname if len(raw_cname) > 28 else "",
                "canonical_smiles": a.get("canonical_smiles"),
                "target_name": target_name,
                "target_chembl_id": a.get("target_chembl_id"),
                "standard_type": std_type,
                "standard_relation": a.get("standard_relation"),
                "standard_value": std_val,
                "standard_units": a.get("standard_units"),
                "pchembl_value": pchembl_val,
                "assay_description": assay_desc,
                "document_chembl_id": doc_clean,
                "document_year": act_yr,
                "_doc_order": doc_idx,
            }
            col_key = _build_endpoint_column_key(
                std_type, target_name, assay_desc
            )
            if col_key and std_val is not None:
              rec[col_key] = std_val
              rec["all_endpoints"] = json.dumps({col_key: std_val})
            all_records.append(rec)
          offset += len(acts)
          total_count = data.get("page_meta", {}).get("total_count", offset)
          if offset >= total_count:
            break
        except Exception as e:
          print(
              f"      Warning: ChEMBL API request failed at offset {offset}:"
              f" {e}",
              flush=True,
          )
          break

      if len(all_records) == doc_records_before and pmid_token:
        for fb_rec in _fetch_pubchem_pmid_fallback(pmid_token):
          fb_rec["_doc_order"] = doc_idx
          all_records.append(fb_rec)

  raw_df = pd.DataFrame(all_records)
  if raw_df.empty or "canonical_smiles" not in raw_df.columns:
    print(
        f"  -> Successfully retrieved {len(raw_df)} raw ChEMBL/PubChem"
        " activity records.",
        flush=True,
    )
    return raw_df

  work = raw_df.copy()
  for str_col in ("paper_compound_id", "compound_name", "iupac_name"):
    if str_col in work.columns:
      work[str_col] = work[str_col].replace("", pd.NA)
  work["_pchembl_num"] = pd.to_numeric(
      work.get("pchembl_value"), errors="coerce"
  )

  def _is_mutant_desc(desc_s: str) -> bool:
    s = str(desc_s).lower()
    if any(k in s for k in ("mutant", "mutation", "variant")):
      return True
    return bool(_AMINO_ACID_POINT_MUTATION_RE.search(s))

  work["_is_wt"] = (
      work.get("assay_description", pd.Series("", index=work.index))
      .astype(str)
      .apply(lambda s: 0 if _is_mutant_desc(s) else 1)
  )
  if "_doc_order" not in work.columns:
    work["_doc_order"] = 0
  work = work.sort_values(
      by=["_doc_order", "_is_wt", "_pchembl_num"],
      ascending=[True, False, False],
      na_position="last",
      kind="stable",
  )
  doc_col = (
      work["document_chembl_id"].fillna(work.get("source", ""))
      if "document_chembl_id" in work.columns
      else work.get("source", "")
  )
  smi_key = work["canonical_smiles"].fillna("").astype(str).str.strip()
  work["_doc_cmpd_key"] = [
      f"{d}||{s}" if s else f"__row_{i}"
      for i, (d, s) in enumerate(zip(doc_col, smi_key))
  ]
  if "all_endpoints" in work.columns:

    def _merge_ep_json(vals: Sequence[object]) -> Optional[str]:
      merged_ep: dict[str, Any] = {}
      for v in vals:
        if pd.notna(v) and str(v).strip():
          try:
            parsed = json.loads(str(v)) if isinstance(v, str) else v
            if isinstance(parsed, dict):
              for ek, ev in parsed.items():
                if ek not in merged_ep and ev is not None:
                  merged_ep[ek] = ev
          except Exception:
            pass
      return json.dumps(merged_ep) if merged_ep else None

    work["all_endpoints"] = work.groupby("_doc_cmpd_key", sort=False)[
        "all_endpoints"
    ].transform(_merge_ep_json)

  wide_df = (
      work.groupby("_doc_cmpd_key", sort=False, as_index=False)
      .first()
      .drop(
          columns=["_pchembl_num", "_is_wt", "_doc_cmpd_key", "_doc_order"],
          errors="ignore",
      )
  )
  if (
      "canonical_smiles" in wide_df.columns
      and "inchikey" not in wide_df.columns
  ):

    def _smi_to_ik(smi_val: Any) -> Optional[str]:
      if not pd.notna(smi_val) or not str(smi_val).strip():
        return None
      mol = Chem.MolFromSmiles(str(smi_val).strip())
      return Chem.MolToInchiKey(mol) if mol is not None else None

    wide_df["inchikey"] = wide_df["canonical_smiles"].apply(_smi_to_ik)
  print(
      f"  -> Successfully retrieved {len(raw_df)} raw ChEMBL/PubChem activity"
      f" records ({len(wide_df)} unique document-compound records with wide"
      " endpoint columns).",
      flush=True,
  )
  return wide_df


def main() -> None:
  parser = argparse.ArgumentParser(
      description=(
          "Academic Paper & Bioassay Literature Harvester: retrieves chemical"
          " structures, paper compound IDs, and pivoted bioassay endpoints"
          " from ChEMBL and PubMed/PubChem by PMID, ChEMBL Document ID, or DOI."
      )
  )
  parser.add_argument(
      "documents",
      nargs="*",
      help=(
          "Positional document identifiers (PMIDs, ChEMBL document IDs"
          " like CHEMBL<DOC_ID>, or DOIs like 10.1021/...)."
      ),
  )
  parser.add_argument(
      "--pmids",
      "--pmid",
      dest="pmids",
      nargs="*",
      default=None,
      help="One or more PubMed IDs (e.g. <PMID> or PMID:<PMID>).",
  )
  parser.add_argument(
      "--chembl-docs",
      "--chembl_docs",
      "--doc-ids",
      "--doc_ids",
      dest="chembl_docs",
      nargs="*",
      default=None,
      help="One or more ChEMBL document IDs (e.g. CHEMBL<DOC_ID>).",
  )
  parser.add_argument(
      "--dois",
      "--doi",
      dest="dois",
      nargs="*",
      default=None,
      help="One or more DOIs (e.g. 10.1021/<DOI_SUFFIX>).",
  )
  parser.add_argument(
      "--output-csv",
      "--output_csv",
      "--output",
      dest="output_csv",
      default="academic_papers.csv",
      help="Output CSV path for harvested academic literature bioassays.",
  )
  parser.add_argument(
      "--output-sdf",
      "--output_sdf",
      dest="output_sdf",
      default=None,
      help="Optional output SDF path.",
  )

  args = parser.parse_args()
  doc_tokens: list[str] = []
  for group in (args.documents, args.pmids, args.chembl_docs, args.dois):
    if group:
      for item in group:
        for part in str(item).split(","):
          if part.strip():
            doc_tokens.append(part.strip())

  if not doc_tokens:
    parser.error(
        "Provide at least one document identifier via --pmids, --chembl-docs,"
        " --dois, or positional arguments."
    )

  df = fetch_chembl_document_activities(doc_tokens)
  if df.empty:
    print(
        "[!] WARNING: 0 bioactivity records found in ChEMBL/PubChem for"
        f" {doc_tokens}. Many journal articles (e.g., ACS Med. Chem. Lett.)"
        " are not deposited in ChEMBL by PMID/DOI; extract structures/tables"
        " from the paper full text or query ChEMBL by target instead.",
        flush=True,
    )
  os.makedirs(
      os.path.dirname(os.path.abspath(args.output_csv)) or ".", exist_ok=True
  )
  df.to_csv(args.output_csv, index=False)
  print(
      f"[+] Saved academic paper bioassay CSV ({len(df)} records,"
      f" {len(df.columns)} columns) to: {args.output_csv}"
  )

  if args.output_sdf and not df.empty and "canonical_smiles" in df.columns:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_sdf)) or ".", exist_ok=True
    )
    with Chem.SDWriter(args.output_sdf) as writer:
      seen_smi: set[str] = set()
      for _, row in df.iterrows():
        smi = str(row.get("canonical_smiles") or "").strip()
        if not smi or smi in seen_smi:
          continue
        seen_smi.add(smi)
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
          compute_clean_2d_coords(mol)
          mol.SetProp(
              "_Name",
              str(
                  row.get("paper_compound_id")
                  or row.get("compound_id")
                  or "Compound"
              ),
          )
          for k, v in row.items():
            if pd.notna(v):
              mol.SetProp(str(k), str(v))
          writer.write(mol)
    print(f"[+] Saved academic paper SDF library to: {args.output_sdf}")

  curate_script = os.path.abspath(
      os.path.join(
          os.path.dirname(os.path.abspath(__file__)),
          "..",
          "..",
          "cheminformatics_core",
          "scripts",
          "curate_congeneric_corpus.py",
      )
  )
  print(
      "[NEXT STEP] Merge academic paper and patent CSVs, align columns"
      " (--col-map), compute RDKit descriptors, and plot chemical space with:\n"
      f"            uv run {curate_script} --input-csv {args.output_csv}"
      " <triaged_patents.csv> --output-csv master_corpus.csv",
      flush=True,
  )


if __name__ == "__main__":
  main()
