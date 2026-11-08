# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "rdkit>=2024.3.1",
#   "pandas>=2.2.0",
#   "beautifulsoup4>=4.12.0",
#   "polite-http",
# ]
# ///
"""Public Patent Chemical Retrieval Utility (Zero-Auth Open-Access).

Queries public chemical databases (EMBL-EBI ChEMBL and NCBI PubChem PUG REST)
by patent number to extract small molecules, chemical structures (SMILES),
and annotated biological activities without requiring private API keys
or proprietary Google-internal services.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import io
import math
import os
import re
import sys
from typing import Any, Dict, List, Optional, Set, Tuple
import urllib.parse
import zlib

from bs4 import BeautifulSoup
import pandas as pd
from polite_http import http_client
from rdkit import Chem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

_CHEMBL_CLIENT = http_client.HttpClient("https://www.ebi.ac.uk/", qps=5.0)
_CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"
_PUBCHEM_CLIENT = http_client.HttpClient(
    "https://pubchem.ncbi.nlm.nih.gov/", qps=3.0
)
_OPSIN_CLIENT = http_client.HttpClient("https://opsin.ch.cam.ac.uk/", qps=5.0)
_GOOGLE_PATENTS_CLIENT = http_client.HttpClient(
    "https://patents.google.com/", qps=5.0
)
_CLIENTS: Dict[str, http_client.HttpClient] = {
    "https://www.ebi.ac.uk/": _CHEMBL_CLIENT,
    "https://pubchem.ncbi.nlm.nih.gov/": _PUBCHEM_CLIENT,
    "https://opsin.ch.cam.ac.uk/": _OPSIN_CLIENT,
    "https://patents.google.com/": _GOOGLE_PATENTS_CLIENT,
}

_PREP_OR_INTERMEDIATE_HEADER_RE = re.compile(
    r"\b(?:prep(?:arations?|\b|\.)|interm(?:ediates?|\b|\.)|"
    r"int\.|starting\s+materials?)"
)


def _expand_wipo_year(yr_str: str) -> str:
  """Expands a 2- or 4-digit WIPO/PCT year using the 1978 PCT pivot."""
  s = yr_str.strip()
  if len(s) == 2 and s.isdigit():
    yr2 = int(s)
    return str(1900 + yr2 if yr2 >= 78 else 2000 + yr2)
  return s


def normalize_patent_token(token: str) -> str:
  """Normalizes PCT and national patent identifiers to canonical standard (e.g.

  WO-YYYY-NNNNNN).
  """
  if not token or not isinstance(token, str):
    return ""
  token = token.strip()
  # Match PCT application formats: PCT/US2021/045123, PCT/EP98/01234
  pct_match = re.match(
      r"^PCT/[A-Z]{2}((?:19|20)\d{2}|\d{2})/(\d+)", token, re.IGNORECASE
  )
  if pct_match:
    yr = _expand_wipo_year(pct_match.group(1))
    return f"WO-{yr}-{pct_match.group(2).zfill(6)}"

  # Match punctuated WIPO numbers: WO-2021-188392-A1, WO 99/12345, WO99/12345A1
  wo_punct_match = re.match(
      r"^WO[-_ /]?((?:19|20)\d{2}|\d{2})[-_ /]+(\d+)(?:[-_ /]?[A-Z0-9]+)?$",
      token,
      re.IGNORECASE,
  )
  if wo_punct_match:
    yr = _expand_wipo_year(wo_punct_match.group(1))
    return f"WO-{yr}-{wo_punct_match.group(2).zfill(6)}"

  # Match unpunctuated 4-digit-year WIPO numbers before 2-digit-year tokens
  wo_4yr_match = re.match(
      r"^WO((?:19|20)\d{2})(\d{4,7})(?:[-_ ]?[A-Z0-9]+)?$",
      token,
      re.IGNORECASE,
  )
  if wo_4yr_match:
    return f"WO-{wo_4yr_match.group(1)}-{wo_4yr_match.group(2).zfill(6)}"

  wo_2yr_match = re.match(
      r"^WO(\d{2})(\d{4,6})(?:[-_ ]?[A-Z0-9]+)?$",
      token,
      re.IGNORECASE,
  )
  if wo_2yr_match:
    yr = _expand_wipo_year(wo_2yr_match.group(1))
    return f"WO-{yr}-{wo_2yr_match.group(2).zfill(6)}"

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
  if (
      formulation_str is None
      or str(formulation_str).strip() == ""
      or str(formulation_str).lower() == "nan"
  ):
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


def normalize_patent_id(raw_patent: str) -> Dict[str, Any]:
  """Normalize a raw patent identifier into standard database lookup variants.

  Handles US, EP, and PCT (WO) applications with year, serial, and kind
  tokenization.
  """
  clean = raw_patent.strip().upper()

  # Pre-process PCT application tokens like PCT/US2021/045123 or PCT/EP98/01234
  m_pct = re.match(
      r"^PCT/[A-Z]{2}((?:19|20)\d{2}|\d{2})/(\d+)", clean, re.IGNORECASE
  )
  if m_pct:
    yr = _expand_wipo_year(m_pct.group(1))
    clean = f"WO{yr}{m_pct.group(2).zfill(6)}"

  # Extract 2-letter authority code
  m_auth = re.match(r"^(WO|EP|US|GB|DE|JP|FR|CA|CN|AU)", clean)
  country = m_auth.group(1) if m_auth else "US"
  rest = clean[len(country) :].strip(" -/.,")

  # Detect kind code at the end (e.g. A1, A2, B1, B2, C1)
  m_kind = re.search(r"([A-Z]\d?)$", rest)
  if m_kind:
    kind = m_kind.group(1)
    body = rest[: -len(kind)].strip(" -/.,")
  else:
    kind = ""
    body = rest

  # Extract pure digits
  digits = re.sub(r"[\s\-\/\.,]", "", body)

  year = None
  serial = digits
  pct_variants = []

  if country == "WO":
    if "/" in body or "-" in body or " " in body:
      parts = [p for p in re.split(r"[\s\-\/]+", body) if p]
      if len(parts) >= 2 and parts[0].isdigit() and len(parts[0]) in (2, 4):
        year = _expand_wipo_year(parts[0])
        serial = "".join(parts[1:])
    elif len(digits) >= 8 and (
        digits.startswith("19") or digits.startswith("20")
    ):
      year = digits[:4]
      serial = digits[4:]
    elif 6 <= len(digits) <= 8 and digits.isdigit():
      year = _expand_wipo_year(digits[:2])
      serial = digits[2:]

    if year:
      serial_6 = serial.zfill(6)
      pct_variants = [
          f"WO-{year}-{serial_6}",
          f"WO-{year}{serial_6}",
          f"WO {year}/{serial_6}",
          f"WO/{year}/{serial_6}",
          f"WO{year}{serial_6}",
      ]
      if kind:
        pct_variants.extend([
            f"WO-{year}{serial_6}-{kind}",
            f"WO-{year}-{serial_6}-{kind}",
            f"WO {year}/{serial_6} {kind}",
            f"WO{year}{serial_6}{kind}",
        ])
      yr2 = year[-2:]
      if int(year) < 2004:
        pct_variants.extend([
            f"WO{yr2}/{serial}",
            f"WO{yr2}/{serial_6}",
            f"WO-{yr2}-{serial}",
            f"WO{yr2}{serial}",
        ])

  canon_digits = (
      f"{year}{serial.zfill(6)}" if (country == "WO" and year) else digits
  )
  compact = f"{country}{canon_digits}{kind}"
  chembl_hyphen = (
      f"WO-{year}-{serial.zfill(6)}-{kind}"
      if (country == "WO" and year and kind)
      else (
          f"WO-{year}-{serial.zfill(6)}"
          if (country == "WO" and year)
          else (f"{country}-{digits}-{kind}" if kind else f"{country}-{digits}")
      )
  )
  pubchem_term = f"{country}{canon_digits}"

  # Build prioritized list of search terms for PubChem and ChEMBL
  search_terms = []
  if pct_variants:
    search_terms.extend(pct_variants)
  search_terms.extend([
      chembl_hyphen,
      compact,
      pubchem_term,
      f"{country}{digits}{kind}",
      f"{country}{digits}",
      f"{country}-{digits}",
      f"{country} {digits}",
      canon_digits,
      digits,
  ])
  seen = set()
  uniq_search_terms = []
  for t in search_terms:
    if t not in seen:
      seen.add(t)
      uniq_search_terms.append(t)

  return {
      "raw": raw_patent,
      "country": country,
      "number": digits,
      "year": year,
      "kind": kind,
      "compact": compact,
      "chembl_variant": chembl_hyphen,
      "pubchem_term": pubchem_term,
      "pct_variants": pct_variants,
      "search_terms": uniq_search_terms,
  }


def _get_client_for_url(url: str) -> http_client.HttpClient:
  parsed = urllib.parse.urlparse(url)
  origin = f"{parsed.scheme}://{parsed.netloc}/"
  if origin not in _CLIENTS:
    qps = 3.0 if "ncbi.nlm.nih.gov" in (parsed.netloc or "") else 5.0
    _CLIENTS[origin] = http_client.HttpClient(origin, qps=qps)
  return _CLIENTS[origin]


def _warn(msg: str) -> None:
  """Emits a diagnostic to stderr so silent degradation is always visible."""
  print(f"[patent_fetcher] WARNING: {msg}", file=sys.stderr, flush=True)


def _http_get_json(
    url: str, timeout: int = 10, max_retries: int = 2
) -> Optional[Dict[str, Any]]:
  """Execute a fast HTTP GET request returning parsed JSON or None, returning None immediately on 404/414."""
  client = _get_client_for_url(url)
  try:
    return client.fetch_json(url, timeout=timeout, max_retries=max_retries)
  except http_client.HttpError as e:
    # 400/404/414 mean "no such record", which is an expected, quiet outcome.
    # Anything else (429 throttling, 5xx, TLS) is a real failure and must not
    # be mistaken by the caller for an empty result.
    if e.status_code not in (400, 404, 414):
      _warn(f"HTTP {e.status_code} fetching JSON from {url}")
    return None
  except (
      OSError,
      EOFError,
      ValueError,
      TypeError,
      LookupError,
      zlib.error,
  ) as e:
    _warn(f"{type(e).__name__} fetching JSON from {url}: {e}")
    return None


def _http_get_text(
    url: str, timeout: int = 10, max_retries: int = 1
) -> Optional[str]:
  """Execute a fast HTTP GET request returning raw text or None, returning None immediately on 404/414."""
  client = _get_client_for_url(url)
  try:
    return client.fetch_text(url, timeout=timeout, max_retries=max_retries)
  except http_client.HttpError as e:
    if e.status_code not in (400, 404, 414):
      _warn(f"HTTP {e.status_code} fetching text from {url}")
    return None
  except (
      OSError,
      EOFError,
      ValueError,
      TypeError,
      LookupError,
      zlib.error,
  ) as e:
    _warn(f"{type(e).__name__} fetching text from {url}: {e}")
    return None


def fetch_patent_chembl(
    patent_dict: Dict[str, Any], limit: Optional[int] = None
) -> List[Dict[str, Any]]:
  """Retrieve chemical structures, metadata, and experimental bioactivities from EMBL-EBI ChEMBL."""
  candidates = []
  if patent_dict.get("search_terms"):
    candidates.extend(patent_dict["search_terms"])
  candidates.extend([
      patent_dict["chembl_variant"],
      patent_dict["compact"],
      f"{patent_dict['country']}-{patent_dict['number']}",
  ])
  seen = set()
  uniq_candidates = []
  for c in candidates:
    if c not in seen:
      seen.add(c)
      uniq_candidates.append(c)

  doc_info = None
  if uniq_candidates:
    url = f"{_CHEMBL_API_BASE}/document.json?patent_id__in={urllib.parse.quote(','.join(uniq_candidates))}"
    data = _http_get_json(url, timeout=4, max_retries=1)
    if data and data.get("documents"):
      doc_info = data["documents"][0]

  if not doc_info:
    return []

  doc_chembl_id = doc_info.get("document_chembl_id")
  doc_title = doc_info.get("title", "")
  raw_doc_yr = doc_info.get("year")
  doc_year: Optional[int] = None
  if raw_doc_yr is not None and not isinstance(raw_doc_yr, bool):
    try:
      yr_int = int(round(float(raw_doc_yr)))
      if 1800 <= yr_int <= 2100:
        doc_year = yr_int
    except (TypeError, ValueError, OverflowError):
      pass

  # Fetch all experimental bioactivity records linked to this document with full pagination
  act_map: Dict[str, Dict[str, Any]] = {}
  act_offset = 0
  while True:
    act_url = f"{_CHEMBL_API_BASE}/activity.json?document_chembl_id={doc_chembl_id}&limit=1000&offset={act_offset}"
    act_data = _http_get_json(act_url)
    batch_acts = (act_data or {}).get("activities") or []
    if not batch_acts:
      break
    for act in batch_acts:
      mol_id = act.get("molecule_chembl_id")
      if not mol_id:
        continue
      pchembl = act.get("pchembl_value")
      std_val = act.get("standard_value")
      std_type = act.get("standard_type")
      std_units = act.get("standard_units")
      target_name = act.get("target_pref_name")
      existing = act_map.setdefault(mol_id, {"all_endpoints": {}})
      if std_type and std_val is not None:
        col_key = f"{std_type}_{target_name or 'assay'}".replace(
            " ", "_"
        ).replace("/", "_")
        existing["all_endpoints"][col_key] = std_val
      # Prioritize records with valid pChEMBL or IC50/Ki for primary columns
      if not existing.get("activity_type") or (
          pchembl and not existing.get("pchembl_value")
      ):
        existing.update({
            "activity_type": std_type,
            "activity_value": std_val,
            "activity_units": std_units,
            "pchembl_value": pchembl,
            "target_name": target_name,
        })
    if len(batch_acts) < 1000:
      break
    act_offset += 1000

  # Fetch compound records linked to this document to obtain genuine molecule_chembl_ids (with pagination)
  compound_records = []
  offset = 0
  while True:
    page_limit = (
        1000
        if (limit is None or limit - len(compound_records) > 1000)
        else max(1, limit - len(compound_records))
    )
    cr_url = f"{_CHEMBL_API_BASE}/compound_record.json?document_chembl_id={doc_chembl_id}&limit={page_limit}&offset={offset}"
    cr_data = _http_get_json(cr_url)
    if not cr_data or not cr_data.get("compound_records"):
      break
    batch_records = cr_data["compound_records"]
    compound_records.extend(batch_records)
    if limit is not None and len(compound_records) >= limit:
      compound_records = compound_records[:limit]
      break
    if len(batch_records) < page_limit:
      break
    offset += len(batch_records)

  if not compound_records:
    return []

  mol_id_to_key = {}
  for cr in compound_records:
    m_id = cr.get("molecule_chembl_id")
    if m_id and m_id not in mol_id_to_key:
      mol_id_to_key[m_id] = (
          cr.get("compound_name") or cr.get("compound_key") or m_id
      )

  mids = list(mol_id_to_key.keys())
  if not mids:
    return []

  # If document-level activities were empty, query activities by molecule_chembl_id__in
  if not act_map:
    for i in range(0, len(mids), 25):
      batch_ids = ",".join(mids[i : i + 25])
      m_act_url = f"{_CHEMBL_API_BASE}/activity.json?molecule_chembl_id__in={batch_ids}&limit=500"
      m_act_data = _http_get_json(m_act_url)
      if m_act_data and m_act_data.get("activities"):
        for act in m_act_data["activities"]:
          mol_id = act.get("molecule_chembl_id")
          if not mol_id:
            continue
          pchembl = act.get("pchembl_value")
          std_val = act.get("standard_value")
          std_type = act.get("standard_type")
          std_units = act.get("standard_units")
          target_name = act.get("target_pref_name")
          existing = act_map.setdefault(mol_id, {"all_endpoints": {}})
          if std_type and std_val is not None:
            col_key = f"{std_type}_{target_name or 'assay'}".replace(
                " ", "_"
            ).replace("/", "_")
            existing["all_endpoints"][col_key] = std_val
          if not existing.get("activity_type") or (
              pchembl and not existing.get("pchembl_value")
          ):
            existing.update({
                "activity_type": std_type,
                "activity_value": std_val,
                "activity_units": std_units,
                "pchembl_value": pchembl,
                "target_name": target_name,
            })

  # Fetch SMILES and properties from molecule.json in batches
  results = []
  for i in range(0, len(mids), 50):
    batch_ids = ",".join(mids[i : i + 50])
    mol_url = f"{_CHEMBL_API_BASE}/molecule.json?molecule_chembl_id__in={batch_ids}&limit=100"
    mol_data = _http_get_json(mol_url)
    if not mol_data or not mol_data.get("molecules"):
      continue

    for m in mol_data["molecules"]:
      m_id = m.get("molecule_chembl_id")
      structs = m.get("molecule_structures")
      if not structs:
        continue
      smiles = structs.get("canonical_smiles")
      if not smiles:
        continue

      pref_name = m.get("pref_name") or mol_id_to_key.get(m_id) or m_id
      act_info = act_map.get(m_id, {})

      mol_entry = {
          "compound_id": pref_name,
          "chembl_id": m_id,
          "smiles": smiles,
          "source": "ChEMBL",
          "patent_id": patent_dict["compact"],
          "document_chembl_id": doc_chembl_id,
          "document_title": doc_title,
          "document_year": doc_year,
          "activity_type": act_info.get("activity_type"),
          "activity_value": act_info.get("activity_value"),
          "activity_units": act_info.get("activity_units"),
          "pchembl_value": act_info.get("pchembl_value"),
          "target_name": act_info.get("target_name"),
      }
      for k, v in act_info.get("all_endpoints", {}).items():
        mol_entry[k] = v
      results.append(mol_entry)

  return results


# Canonical organic elements & bioisosteres (H, B, C, N, O, F, Si, P, S, Cl, Se, Br, I)
ORGANIC_ELEMENTS: Set[int] = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 34, 35, 53}

# MedChem Metallodrugs, Coordination Complexes & Radiopharmaceuticals:
# Pt (cisplatin/oxaliplatin), Au (auranofin), Ru (NAMI-A/KP1019), Fe (ferrocenes/ferroquine),
# Cu (casiopeínas), Zn (zinc complexes), Gd (MRI contrast), Lu/Y/Tc/Ga/In (radiopharmaceuticals),
# Ti, V, Mn, Co, Ni, As (Trisenox), Bi (bismuth subsalicylate), Ag, etc.
METALLODRUG_ELEMENTS: Set[int] = {
    3,
    12,
    13,
    20,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    32,
    33,
    38,
    39,
    40,
    42,
    43,
    44,
    45,
    46,
    47,
    48,
    49,
    50,
    51,
    56,
    57,
    58,
    62,
    63,
    64,
    65,
    67,
    68,
    70,
    71,
    73,
    74,
    75,
    76,
    77,
    78,
    79,
    80,
    83,
}

# By default, allow organic small molecules as well as metallodrugs and coordination complexes:
DEFAULT_ALLOWED_ELEMENTS: Set[int] = ORGANIC_ELEMENTS | METALLODRUG_ELEMENTS


def fetch_patent_pubchem(
    patent_dict: Dict[str, Any],
    limit: Optional[int] = None,
    allowed_elements: Optional[Set[int]] = None,
) -> List[Dict[str, Any]]:
  """Retrieve chemical structures from NCBI PubChem using authoritative PUG REST xref/PatentID and E-utilities."""
  cids: List[str] = []

  # 1. Primary authoritative lookup via PubChem PUG REST xref/PatentID
  xref_candidates = [
      patent_dict.get("compact", ""),
      patent_dict.get("pubchem_term", ""),
      patent_dict.get("chembl_variant", ""),
  ]
  for st in patent_dict.get("search_terms", []):
    if not st.isdigit() and st not in xref_candidates:
      xref_candidates.append(st.replace(" ", "").replace("-", ""))

  seen_xref = set()
  for cand in xref_candidates:
    if not cand or cand in seen_xref or cand.isdigit():
      continue
    seen_xref.add(cand)
    xref_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/PatentID/{urllib.parse.quote(cand)}/cids/JSON"
    data = _http_get_json(xref_url)
    if data and data.get("IdentifierList", {}).get("CID"):
      found = [str(c) for c in data["IdentifierList"]["CID"]]
      if len(found) > len(cids):
        cids = found
      if cids:
        break

  # 2. Fallback to Entrez E-utilities [Patent] tag search only if xref returned no CIDs (excluding bare digits)
  if not cids:
    query_terms = [
        f"{patent_dict['pubchem_term']}[Patent]",
        f"{patent_dict['compact']}[Patent]",
    ]
    for st in patent_dict.get("search_terms", []):
      if not st.isdigit():
        query_terms.append(f"{st}[Patent]")

    seen = set()
    retmax_val = limit if limit is not None else 100000
    for term in query_terms:
      if term in seen:
        continue
      seen.add(term)
      esearch_url = (
          "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
          f"db=pccompound&term={urllib.parse.quote(term)}&retmode=json&retmax={retmax_val}"
      )
      data = _http_get_json(esearch_url)
      if data and data.get("esearchresult", {}).get("idlist"):
        found = data["esearchresult"]["idlist"]
        if len(found) > len(cids):
          cids = found
        if limit is not None and len(cids) >= limit:
          break

  if not cids:
    return []

  # Sort CIDs descending so newly patented example molecules come before common
  # building blocks, while retaining enough xref CIDs (at least 2500) so
  # lower-numbered canonical stereoisomer CIDs are not dropped.
  cids_sorted = sorted(
      [c for c in cids if str(c).isdigit()], key=int, reverse=True
  ) + [c for c in cids if not str(c).isdigit()]
  if limit is not None:
    cids_sorted = cids_sorted[: max(limit * 6, 2500)]

  print(
      f"[patent_fetcher] Fetching PubChem properties for {len(cids_sorted)}"
      " CIDs...",
      file=sys.stderr,
      flush=True,
  )

  batches = [cids_sorted[i : i + 150] for i in range(0, len(cids_sorted), 150)]
  results = []

  def _fetch_prop_batch(batch_cids: List[str]) -> List[Dict[str, Any]]:
    cids_str = ",".join(batch_cids)
    prop_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids_str}/"
        "property/SMILES,ConnectivitySMILES,MolecularFormula,IUPACName,Title/JSON"
    )
    prop_data = _http_get_json(prop_url, timeout=12)
    if not prop_data or not prop_data.get("PropertyTable", {}).get(
        "Properties"
    ):
      return []
    batch_res = []
    for p in prop_data["PropertyTable"]["Properties"]:
      cid = p.get("CID")
      smiles = p.get("SMILES") or p.get("ConnectivitySMILES")
      if not smiles:
        continue
      mol = Chem.MolFromSmiles(smiles)
      if not mol:
        continue
      frags = Chem.GetMolFrags(mol, asMols=True)
      if len(frags) > 1:
        mol = max(frags, key=lambda f: f.GetNumHeavyAtoms())
        smiles = Chem.MolToSmiles(mol)
      if mol.GetNumHeavyAtoms() < 12:
        continue
      cur_allowed = (
          allowed_elements
          if allowed_elements is not None
          else DEFAULT_ALLOWED_ELEMENTS
      )
      if cur_allowed is not None:
        if any(
            atom.GetAtomicNum() not in cur_allowed for atom in mol.GetAtoms()
        ):
          continue
      raw_iupac = p.get("IUPACName") or ""
      # PubChem returns a semicolon-joined string for multi-component and
      # parent/child CIDs. The structure above was reduced to its largest
      # fragment, so carrying the joined name forward labels the row with a
      # different compound's name. Keep one name, preserve the original, and
      # record that the choice was ambiguous rather than hiding it.
      name_parts = [s.strip() for s in raw_iupac.split(";") if s.strip()]
      name_ambiguous = len(name_parts) > 1
      if name_parts:
        iupac = max(name_parts, key=len)
      else:
        iupac = f"CID_{cid}"
      # Title is PubChem's preferred common name and is served by the same
      # request, so compound_name need not duplicate the IUPAC string.
      title = (p.get("Title") or "").strip() or iupac
      batch_res.append({
          "compound_id": f"CID_{cid}",
          "pubchem_cid": cid,
          "smiles": smiles,
          "canonical_smiles": Chem.MolToSmiles(mol),
          "inchikey": Chem.MolToInchiKey(mol),
          "iupac_name": iupac,
          "compound_name": title,
          "pubchem_name_raw": raw_iupac,
          "name_is_ambiguous": name_ambiguous,
          "molecular_formula": p.get("MolecularFormula", ""),
          "source": "PubChem",
          "patent_id": patent_dict["compact"],
      })
    return batch_res

  with concurrent.futures.ThreadPoolExecutor(max_workers=6) as executor:
    futures = [executor.submit(_fetch_prop_batch, b) for b in batches]
    for idx, fut in enumerate(concurrent.futures.as_completed(futures), 1):
      try:
        results.extend(fut.result())
      except Exception:
        pass
      if idx % 5 == 0 or idx == len(batches):
        print(
            f"[patent_fetcher] PubChem property batches: {idx}/{len(batches)}"
            f" ({len(results)} valid organic compounds)",
            file=sys.stderr,
            flush=True,
        )

  results.sort(
      key=lambda r: int(r.get("pubchem_cid") or 0)
      if str(r.get("pubchem_cid") or "").isdigit()
      else 0,
      reverse=True,
  )
  return results


_OPSIN_CACHE: Dict[str, Optional[str]] = {}
_CHEMBL_IK_CACHE: Dict[str, Dict[str, Any]] = {}
_PUBCHEM_BIOASSAY_CACHE: Dict[str, Dict[str, Any]] = {}


def enrich_records_with_pubchem_bioassays(
    records: List[Dict[str, Any]],
) -> None:
  """Enriches PubChem patent records with quantitative IC50/Ki/Kd/EC50 bioactivity from PubChem BioAssay summaries."""
  has_discrete_bins = any(r.get("activity_bin") is not None for r in records)
  cid_to_rec: Dict[str, Dict[str, Any]] = {}
  for r in records:
    if (
        r.get("pchembl_value") is not None
        or r.get("activity_value") is not None
        or r.get("activity_bin") is not None
    ):
      continue
    c = str(r.get("pubchem_cid", "")).strip()
    if not c:
      continue
    if c in _PUBCHEM_BIOASSAY_CACHE:
      cached = _PUBCHEM_BIOASSAY_CACHE[c]
      if cached:
        r.update(cached)
      continue
    cid_to_rec[c] = r

  cids = list(cid_to_rec.keys())
  if not cids:
    return

  for c in cids:
    _PUBCHEM_BIOASSAY_CACHE.setdefault(c, {})

  batches = [cids[i : i + 20] for i in range(0, len(cids), 20)]
  print(
      f"[patent_fetcher] Enriching {len(cids)} compounds with PubChem BioAssays"
      f" ({len(batches)} batches)...",
      file=sys.stderr,
      flush=True,
  )

  def _fetch_bioassay_batch(batch: List[str]) -> Optional[str]:
    cids_str = ",".join(batch)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids_str}/assaysummary/CSV"
    return _http_get_text(url, timeout=12, max_retries=0)

  with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
    futures = [executor.submit(_fetch_bioassay_batch, b) for b in batches]
    for idx, fut in enumerate(concurrent.futures.as_completed(futures), 1):
      try:
        text = fut.result()
      except Exception:
        text = None
      if text:
        reader = csv.DictReader(io.StringIO(text))
        for row in reader:
          cid = str(row.get("CID", "")).strip()
          val_str = str(row.get("Activity Value [uM]", "")).strip()
          act_name = str(row.get("Activity Name", "")).strip()
          assay_name = str(row.get("Assay Name", "")).strip()
          outcome = str(row.get("Activity Outcome", "")).strip().lower()
          if not cid or not val_str or cid not in cid_to_rec:
            continue
          if outcome in ("inactive", "inconclusive", "unspecified"):
            continue
          try:
            val_um = float(val_str)
            if val_um <= 0:
              continue
            p_val = round(6.0 - math.log10(val_um), 3)
            val_nm = round(val_um * 1000.0, 4)
          except Exception:
            continue
          # ChEMBL defines pChEMBL only for quantitative potencies <= 100 uM
          # (pChEMBL >= 4.0).
          if p_val < 4.0:
            continue
          # When the patent already uses discrete potency bins, skip secondary
          # patent-scraped numeric approximations ("from US Patent ...") so bin
          # upper bounds do not masquerade as continuous IC50s.
          if has_discrete_bins and re.search(
              r"\bfrom\s+(?:[A-Z]{2}\s+)?Patent\b",
              assay_name,
              re.IGNORECASE,
          ):
            continue
          rec = cid_to_rec[cid]
          existing_p = rec.get("pchembl_value")
          is_preferred = act_name in ("IC50", "Ki", "Kd", "EC50", "Potency")

          should_update = existing_p is None or (
              is_preferred and p_val > float(existing_p)
          )
          if should_update:
            upd = {
                "activity_type": act_name or "Potency",
                "activity_value": val_nm,
                "activity_units": "nM",
                "pchembl_value": p_val,
                "target_name": assay_name[:100],
            }
            rec.update(upd)
            _PUBCHEM_BIOASSAY_CACHE[cid] = upd
      if idx % 10 == 0 or idx == len(batches):
        print(
            f"[patent_fetcher] PubChem BioAssay batches: {idx}/{len(batches)}"
            " complete",
            file=sys.stderr,
            flush=True,
        )


def enrich_records_with_chembl_by_inchikey(
    records: List[Dict[str, Any]],
    fetch_activities: bool = True,
) -> None:
  """Cross-reference validated RDKit molecules against ChEMBL using standard InChIKeys to attach ChEMBL IDs and (optionally) bioactivity data."""
  if not records:
    return

  inchikey_to_records: Dict[str, List[Dict[str, Any]]] = {}
  for r in records:
    mol = r.get("rdkit_mol")
    if mol is None and r.get("canonical_smiles"):
      mol = Chem.MolFromSmiles(str(r["canonical_smiles"]))
      r["rdkit_mol"] = mol
    if mol is None:
      continue
    try:
      ik = Chem.MolToInchiKey(mol)
      if ik:
        r["inchikey"] = ik
        needs_act = (
            fetch_activities
            and r.get("pchembl_value") is None
            and r.get("activity_value") is None
            and r.get("activity_bin") is None
        )
        if ik in _CHEMBL_IK_CACHE:
          cached_ik = _CHEMBL_IK_CACHE[ik]
          if cached_ik.get("chembl_id") and not r.get("chembl_id"):
            r["chembl_id"] = cached_ik["chembl_id"]
          if cached_ik.get("pref_name") and str(
              r.get("compound_id") or ""
          ).startswith("CID_"):
            r["compound_id"] = cached_ik["pref_name"]
          act_info = cached_ik.get("act_info")
          if act_info and needs_act:
            r["activity_type"] = act_info.get("activity_type")
            r["activity_value"] = act_info.get("activity_value")
            r["activity_units"] = act_info.get("activity_units")
            r["pchembl_value"] = act_info.get("pchembl_value")
            r["target_name"] = act_info.get("target_name")
            for k, v in act_info.get("all_endpoints", {}).items():
              r[k] = v
          if not (needs_act and "act_info" not in cached_ik):
            continue
        if not r.get("chembl_id") or needs_act:
          inchikey_to_records.setdefault(ik, []).append(r)
    except Exception:
      pass

  if not inchikey_to_records:
    return

  all_iks = list(inchikey_to_records.keys())
  for ik in all_iks:
    _CHEMBL_IK_CACHE.setdefault(ik, {})
  chembl_id_to_ik: Dict[str, str] = {}
  iks_to_lookup = [
      ik for ik in all_iks if not _CHEMBL_IK_CACHE.get(ik, {}).get("chembl_id")
  ]
  for ik in all_iks:
    cached_cid = _CHEMBL_IK_CACHE.get(ik, {}).get("chembl_id")
    if cached_cid:
      chembl_id_to_ik[cached_cid] = ik

  batches = [
      iks_to_lookup[i : i + 25] for i in range(0, len(iks_to_lookup), 25)
  ]
  if batches:
    print(
        f"[patent_fetcher] Cross-referencing {len(iks_to_lookup)} InChIKeys"
        f" with ChEMBL ({len(batches)} batches)...",
        file=sys.stderr,
        flush=True,
    )
    circuit_open = False
    consecutive_none = 0

    def _fetch_ik_batch(batch: List[str]) -> Optional[Dict[str, Any]]:
      if circuit_open:
        return None
      url = (
          f"{_CHEMBL_API_BASE}/molecule.json?"
          f"molecule_structures__standard_inchi_key__in={','.join(batch)}&limit=100"
      )
      return _http_get_json(url, timeout=6, max_retries=1)

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
      futures = [executor.submit(_fetch_ik_batch, b) for b in batches]
      completed_count = 0
      for fut in concurrent.futures.as_completed(futures):
        completed_count += 1
        if completed_count % 5 == 0 or completed_count == len(batches):
          print(
              "[patent_fetcher] ChEMBL InChIKey lookup progress:"
              f" {completed_count}/{len(batches)} batches complete",
              file=sys.stderr,
              flush=True,
          )
        try:
          data = fut.result()
        except Exception:
          data = None
        if data is None:
          consecutive_none += 1
          if consecutive_none >= 3 and not circuit_open:
            circuit_open = True
            print(
                "[patent_fetcher] ChEMBL InChIKey API returned consecutive"
                " errors/timeouts; circuit breaker skipping remaining batches.",
                file=sys.stderr,
                flush=True,
            )
          continue
        consecutive_none = 0
        if not data.get("molecules"):
          continue
        for m in data["molecules"]:
          m_id = m.get("molecule_chembl_id")
          structs = m.get("molecule_structures") or {}
          ik = structs.get("standard_inchi_key")
          pref_name = m.get("pref_name")
          if m_id and ik and ik in inchikey_to_records:
            chembl_id_to_ik[m_id] = ik
            _CHEMBL_IK_CACHE[ik]["chembl_id"] = m_id
            if pref_name:
              _CHEMBL_IK_CACHE[ik]["pref_name"] = pref_name
            for r in inchikey_to_records[ik]:
              if not r.get("chembl_id"):
                r["chembl_id"] = m_id
              if pref_name and str(r.get("compound_id") or "").startswith(
                  "CID_"
              ):
                r["compound_id"] = pref_name

  if not fetch_activities or not chembl_id_to_ik:
    return

  mids = [
      m_id
      for m_id, ik in chembl_id_to_ik.items()
      if any(
          r.get("pchembl_value") is None
          and r.get("activity_value") is None
          and r.get("activity_bin") is None
          for r in inchikey_to_records.get(ik, [])
      )
  ]
  if not mids:
    return

  for i in range(0, len(mids), 40):
    batch_mids = ",".join(mids[i : i + 40])
    act_url = f"{_CHEMBL_API_BASE}/activity.json?molecule_chembl_id__in={batch_mids}&limit=500"
    act_data = _http_get_json(act_url, timeout=10)
    if not act_data or not act_data.get("activities"):
      continue

    best_acts: Dict[str, Dict[str, Any]] = {}
    for act in act_data["activities"]:
      m_id = act.get("molecule_chembl_id")
      if not m_id:
        continue
      pchembl = act.get("pchembl_value")
      std_val = act.get("standard_value")
      std_type = act.get("standard_type")
      std_units = act.get("standard_units")
      target_name = act.get("target_pref_name")
      tgt_clean = str(target_name or "").strip().lower()
      if tgt_clean in (
          "",
          "unchecked",
          "unspecified",
          "unknown",
          "no relevant target",
          "admet",
      ):
        continue
      try:
        if std_val is not None and float(std_val) <= 0:
          continue
      except (TypeError, ValueError):
        continue
      existing = best_acts.setdefault(m_id, {"all_endpoints": {}})
      raw_doc_year = act.get("document_year")
      if raw_doc_year is not None and not isinstance(raw_doc_year, bool):
        try:
          yr_int = int(float(raw_doc_year))
          if 1800 <= yr_int <= 2100:
            prev_yr = existing.get("document_year")
            if prev_yr is None or yr_int < prev_yr:
              existing["document_year"] = yr_int
        except (TypeError, ValueError):
          pass
      if std_type and std_val is not None:
        col_key = f"{std_type}_{target_name or 'assay'}".replace(
            " ", "_"
        ).replace("/", "_")
        existing["all_endpoints"][col_key] = std_val
      if not existing.get("activity_type") or (
          pchembl and not existing.get("pchembl_value")
      ):
        existing.update({
            "activity_type": std_type,
            "activity_value": std_val,
            "activity_units": std_units,
            "pchembl_value": pchembl,
            "target_name": target_name,
        })

    for m_id, act_info in best_acts.items():
      ik = chembl_id_to_ik.get(m_id)
      if not ik:
        continue
      _CHEMBL_IK_CACHE[ik]["act_info"] = act_info
      for r in inchikey_to_records.get(ik, []):
        if act_info.get("document_year") is not None and (
            r.get("document_year") is None
            or isinstance(r.get("document_year"), bool)
        ):
          r["document_year"] = act_info["document_year"]
        if (
            r.get("pchembl_value") is None
            and r.get("activity_value") is None
            and r.get("activity_bin") is None
        ):
          r["activity_type"] = act_info.get("activity_type")
          r["activity_value"] = act_info.get("activity_value")
          r["activity_units"] = act_info.get("activity_units")
          r["pchembl_value"] = act_info.get("pchembl_value")
          r["target_name"] = act_info.get("target_name")
          for k, v in act_info.get("all_endpoints", {}).items():
            r[k] = v


_ASSAY_HEADER_NON_TARGET_WORDS = frozenset({
    "wt",
    "wild",
    "type",
    "wild-type",
    "ki",
    "kd",
    "ic50",
    "ec50",
    "gi50",
    "ed50",
    "pm",
    "nm",
    "um",
    "mm",
    "cell",
    "cells",
    "cellular",
    "enzymatic",
    "enzyme",
    "biochemical",
    "assay",
    "assays",
    "activity",
    "inhibition",
    "potency",
    "mean",
    "avg",
    "average",
    "value",
    "values",
    "col",
    "column",
    "table",
    "example",
    "compound",
    "human",
    "recombinant",
})


def _extract_target_from_assay_header(hdr_lbl: str) -> Optional[str]:
  """Extracts the primary biological target token from a table assay column header."""
  cleaned = re.sub(r"\([^)]*\)", " ", str(hdr_lbl or ""))
  for tok in re.findall(r"[A-Za-z0-9]+(?:-[A-Za-z0-9]+)?", cleaned):
    if (
        len(tok) >= 3
        and not tok.isdigit()
        and tok.lower() not in _ASSAY_HEADER_NON_TARGET_WORDS
    ):
      return tok
  return None


_INACTIVE_OR_MISSING_TOKENS = frozenset({
    "",
    "na",
    "nd",
    "nt",
    "n/a",
    "n.d.",
    "n.t.",
    "none",
    "nan",
    "null",
    "-",
    "—",
    "–",
    "not determined",
    "not tested",
})


def _is_discrete_potency_token(val_str: Any) -> bool:
  """Returns True if val_str is a discrete/categorical potency grade or range bin."""
  if val_str is None:
    return False
  s = str(val_str).strip()
  if not s or s.lower() in _INACTIVE_OR_MISSING_TOKENS:
    return False
  # Strip optional trailing range annotation like "A (<= 100 nM)"
  s_head = re.sub(r"\s*\([^)]*\)\s*$", "", s).strip()
  if re.match(r"^[A-F][+-]?$", s_head):
    return True
  if re.match(r"^(?:\+{1,6}|\*{1,6}|active|inactive)$", s_head, re.IGNORECASE):
    return True
  if re.match(
      r"^\d+(?:\.\d+)?\s*(?:pM|nM|uM|µM|μM|mM)?\s*(?:[-–—]|to)\s*\d+(?:\.\d+)?\s*(?:pM|nM|uM|µM|μM|mM)$",
      s,
      re.IGNORECASE,
  ):
    return True
  return False


def _extract_patent_dates_from_soup(
    soup: Any,
    patent_dict: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
  """Extracts publication, filing, and priority dates from Google Patents HTML."""
  pub_date: Optional[str] = None
  filing_date: Optional[str] = None
  priority_date: Optional[str] = None

  def _norm_iso_date(raw_val: Optional[str]) -> Optional[str]:
    if not raw_val:
      return None
    m_dt = re.search(r"\b(19\d{2}|20\d{2})[-/](\d{2})[-/](\d{2})\b", raw_val)
    if m_dt:
      return f"{m_dt.group(1)}-{m_dt.group(2)}-{m_dt.group(3)}"
    return None

  if soup is not None:
    for prop_name, target_attr in (
        ("publicationDate", "pub"),
        ("filingDate", "filing"),
        ("priorityDate", "priority"),
    ):
      for tag_name in ("time", "dd", "td", "span"):
        el = soup.find(tag_name, {"itemprop": prop_name})
        if el is None:
          continue
        cand_dt = _norm_iso_date(el.get("datetime")) or _norm_iso_date(
            el.get_text(" ", strip=True)
        )
        if cand_dt:
          if target_attr == "pub" and pub_date is None:
            pub_date = cand_dt
          elif target_attr == "filing" and filing_date is None:
            filing_date = cand_dt
          elif target_attr == "priority" and priority_date is None:
            priority_date = cand_dt
          break

    for meta in soup.find_all(
        "meta", {"name": re.compile(r"^DC\.date$", re.I)}
    ):
      scheme = str(meta.get("scheme") or "").strip().lower()
      content_dt = _norm_iso_date(meta.get("content"))
      if not content_dt:
        continue
      if "submitted" in scheme and filing_date is None:
        filing_date = content_dt
      elif ("issue" in scheme or not scheme) and pub_date is None:
        pub_date = content_dt

  doc_year: Optional[int] = None
  if pub_date:
    yr_val = int(pub_date[:4])
    if 1800 <= yr_val <= 2100:
      doc_year = yr_val

  if doc_year is None and patent_dict:
    raw_yr = patent_dict.get("year")
    if raw_yr and str(raw_yr).isdigit():
      yr_val = int(str(raw_yr))
      if 1800 <= yr_val <= 2100:
        doc_year = yr_val
    if doc_year is None:
      compact_id = str(patent_dict.get("compact") or "")
      m_us_a = re.match(r"^US(20\d{2})\d{6,8}A\d$", compact_id, re.IGNORECASE)
      if m_us_a:
        doc_year = int(m_us_a.group(1))

  if doc_year is None:
    for date_cand in (filing_date, priority_date):
      if date_cand:
        yr_val = int(date_cand[:4])
        if 1800 <= yr_val <= 2100:
          doc_year = yr_val
          break

  return {
      "publication_date": pub_date,
      "filing_date": filing_date,
      "priority_date": priority_date,
      "document_year": doc_year,
  }


def _apply_patent_dates_to_records(
    records: List[Dict[str, Any]],
    date_info: Dict[str, Any],
) -> None:
  """Backfills missing, boolean, or sentinel date fields across patent records."""
  if not date_info:
    return
  for rec in records:
    for key in (
        "publication_date",
        "filing_date",
        "priority_date",
        "document_year",
    ):
      new_val = date_info.get(key)
      if new_val is None:
        continue
      cur_val = rec.get(key)
      if key == "document_year":
        valid_cur = False
        if cur_val is not None and not isinstance(cur_val, bool):
          try:
            yr_i = int(round(float(cur_val)))
            valid_cur = 1800 <= yr_i <= 2100
          except (TypeError, ValueError, OverflowError):
            valid_cur = False
        if not valid_cur:
          rec[key] = new_val
      else:
        cur_str = "" if cur_val is None else str(cur_val).strip()
        if (
            cur_val is None
            or isinstance(cur_val, bool)
            or not re.match(r"^(19\d{2}|20\d{2})-\d{2}-\d{2}$", cur_str)
        ):
          rec[key] = new_val


def _extract_patent_bin_legend(html: str) -> Dict[str, str]:
  """Extracts discrete potency bin definitions (e.g., A <= 100 nM, B 100-250 nM) from patent specification text."""
  legend: Dict[str, str] = {}
  if not html:
    return legend
  text = re.sub(r"<[^>]+>", " ", html)
  text = re.sub(r"\s+", " ", text)

  def _norm_bound(op_phrase: str, val_s: str, unit_s: str) -> str:
    u = unit_s.replace("μ", "u").replace("µ", "u")
    op_low = op_phrase.lower()
    if "less than" in op_low and "equal" in op_low:
      return f"<= {val_s} {u}"
    if "less than" in op_low or "<" in op_low:
      return f"< {val_s} {u}"
    if "greater than" in op_low and "equal" in op_low:
      return f">= {val_s} {u}"
    if "greater than" in op_low or ">" in op_low:
      return f"> {val_s} {u}"
    return f"{val_s} {u}"

  # 1. Dual-bounded prose: "greater than X <unit> [but/and] less than or equal
  # to Y <unit> are represented as 'B'"
  for m in re.finditer(
      r"greater\s+than(?:\s+or\s+equal\s+to)?\s+(\d+(?:\.\d+)?)\s*(pM|nM|uM|µM|μM|mM)\s+"
      r"(?:but|and)\s+less\s+than(?:\s+or\s+equal\s+to)?\s+(\d+(?:\.\d+)?)\s*(pM|nM|uM|µM|μM|mM)\s+"
      r"(?:are|is)\s+(?:represented|designated|indicated|assigned|scored)\s+(?:as|by)\s+"
      r"[\"'“”‘’]?([A-F]|\+{1,5}|13)[\"'“”‘’]?",
      text,
      re.IGNORECASE,
  ):
    v1, u1, v2, u2, raw_code = m.groups()
    code = "B" if raw_code == "13" else raw_code.upper()
    u1_n = u1.replace("μ", "u").replace("µ", "u")
    u2_n = u2.replace("μ", "u").replace("µ", "u")
    rng = (
        f"{v1}-{v2} {u2_n}"
        if u1_n.lower() == u2_n.lower()
        else f"{v1} {u1_n} - {v2} {u2_n}"
    )
    legend[code] = rng

  # 2. Single-bounded prose: "(less|greater) than X <unit> are represented as A"
  for m in re.finditer(
      r"(?<!but\s)(?<!and\s)((?:less|greater)\s+than(?:\s+or\s+equal\s+to)?)\s+"
      r"(\d+(?:\.\d+)?)\s*(pM|nM|uM|µM|μM|mM)\s+"
      r"(?:are|is)\s+(?:represented|designated|indicated|assigned|scored)\s+(?:as|by)\s+"
      r"[\"'“”‘’]?([A-F]|\+{1,5})[\"'“”‘’]?",
      text,
      re.IGNORECASE,
  ):
    op_p, val_s, unit_s, raw_code = m.groups()
    code = raw_code.upper()
    legend.setdefault(code, _norm_bound(op_p, val_s, unit_s))

  # 3. Symbolic footnote definitions: "A: <= 100 nM", "'A' means IC50 < 100 nM"
  for m in re.finditer(
      r"[\"'“”‘’]?\b([A-F])\b[\"'“”‘’]?\s*(?:[:=]|means|represents|is)\s*(?:IC\s*50|Ki|EC\s*50|Kd)?\s*"
      r"([<>≤≥]=?\s*\d+(?:\.\d+)?\s*(?:pM|nM|uM|µM|μM|mM)|\d+(?:\.\d+)?\s*(?:pM|nM|uM|µM|μM|mM)?\s*[-–]\s*\d+(?:\.\d+)?\s*(?:pM|nM|uM|µM|μM|mM))",
      text,
      re.IGNORECASE,
  ):
    code = m.group(1).upper()
    bound_s = (
        m.group(2)
        .replace("μ", "u")
        .replace("µ", "u")
        .replace("≤", "<=")
        .replace("≥", ">=")
    )
    legend.setdefault(code, re.sub(r"\s+", " ", bound_s).strip())

  return legend


def _parse_p(val_str: str) -> Tuple[Optional[float], Optional[float]]:
  """Parses a raw potency string with units into (nM_value, pActivity)."""
  if not val_str or str(val_str).strip().upper() in ("NA", "ND", "N/A", "-"):
    return None, None
  s = str(val_str).strip()
  # Do not misparse range strings ("100-250 nM") as a single scalar 100 nM
  if re.match(
      r"^\d+(?:\.\d+)?\s*(?:pM|nM|uM|µM|μM|mM)?\s*(?:[-–—]|to)\s*\d+(?:\.\d+)?\s*(?:pM|nM|uM|µM|μM|mM)$",
      s,
      re.IGNORECASE,
  ):
    return None, None
  m = re.search(
      r"([0-9]+(?:\.[0-9]+)?)\s*(pM|nM|uM|µM|μM|mM)?",
      s,
      re.IGNORECASE,
  )
  if not m:
    return None, None
  num = float(m.group(1))
  unit = (m.group(2) or "nM").lower()
  factor = {
      "pm": 1e-12,
      "nm": 1e-9,
      "um": 1e-6,
      "µm": 1e-6,
      "μm": 1e-6,
      "mm": 1e-3,
  }.get(unit, 1e-9)
  val_m = num * factor
  if val_m <= 0:
    return None, None
  return val_m * 1e9, round(-math.log10(val_m), 3)


def _inject_assay_columns(
    target_dict: Dict[str, Any],
    assays: List[str],
    header_labels: Optional[List[str]] = None,
    bin_legend: Optional[Dict[str, str]] = None,
) -> None:
  """Injects 1-based patent_table_col_{k} and header-named patent_assay_* columns while preserving column indices across empty cells."""
  pat_prefix = re.sub(
      r"[^A-Za-z0-9]+",
      "",
      str(target_dict.get("patent_id") or "PAT"),
  )

  for idx_a, a_val in enumerate(assays, 1):
    if not a_val or not str(a_val).strip():
      continue
    clean_val = str(a_val).strip()
    target_dict[f"patent_table_col_{idx_a}"] = a_val
    target_dict[f"{pat_prefix}_col_{idx_a}"] = a_val
    nm_v, p_v = _parse_p(a_val)
    if nm_v is not None:
      target_dict[f"patent_table_col_{idx_a}_nM"] = round(nm_v, 4)
      target_dict[f"patent_table_col_{idx_a}_pActivity"] = p_v
      target_dict[f"{pat_prefix}_col_{idx_a}_nM"] = round(nm_v, 4)
    elif _is_discrete_potency_token(clean_val):
      target_dict.setdefault("activity_bin", clean_val)
      if bin_legend and clean_val.upper() in bin_legend:
        rng_lbl = bin_legend[clean_val.upper()]
        target_dict.setdefault("activity_range", rng_lbl)
        target_dict[f"patent_table_col_{idx_a}_range"] = rng_lbl
    if header_labels and idx_a - 1 < len(header_labels):
      hdr_lbl = header_labels[idx_a - 1]
      norm_hdr_lbl = hdr_lbl.replace("µ", "u").replace("μ", "u")
      clean_hdr = re.sub(r"[^a-zA-Z0-9]+", "_", norm_hdr_lbl).strip("_").lower()
      if clean_hdr and not clean_hdr.startswith("col_"):
        col_key = f"patent_assay_{clean_hdr}"
        target_dict[col_key] = a_val
        if nm_v is not None:
          target_dict[f"{col_key}_nM"] = round(nm_v, 4)
          target_dict[f"{col_key}_pActivity"] = p_v
        elif (
            _is_discrete_potency_token(clean_val)
            and bin_legend
            and clean_val.upper() in bin_legend
        ):
          target_dict[f"{col_key}_range"] = bin_legend[clean_val.upper()]
      if not target_dict.get("target_name"):
        tgt_tok = _extract_target_from_assay_header(hdr_lbl)
        if tgt_tok:
          target_dict["target_name"] = tgt_tok


def _cell_text_no_sup(cell_tag: Any) -> str:
  """Extracts text from a table cell after stripping <sup> footnote tags."""
  raw_c = str(cell_tag)
  if "<sup" in raw_c.lower():
    raw_c = re.sub(
        r"<sup[^>]*>.*?</sup>", "", raw_c, flags=re.IGNORECASE | re.DOTALL
    )
    return BeautifulSoup(raw_c, "html.parser").get_text(" ", strip=True)
  return cell_tag.get_text(" ", strip=True)


def _looks_like_table_iupac(txt: str) -> bool:
  """Returns True if a table cell string resembles an IUPAC chemical name."""
  s = txt.strip()
  if len(s) < 20 or "δ" in s or "ppm" in s.lower() or "mhz" in s.lower():
    return False
  if re.search(r"\bJ\s*=\s*\d", s) or re.match(r"^\d+\s*%", s):
    return False
  low = s.lower()
  return "-" in s and any(
      k in low
      for k in (
          "amino",
          "fluoro",
          "chloro",
          "bromo",
          "methyl",
          "ethyl",
          "methoxy",
          "oxy",
          "oxo",
          "pyrazol",
          "pyridin",
          "pyrimidin",
          "phenyl",
          "benzo",
          "carbonitrile",
          "carboxamide",
          "yl)",
          "yl]",
      )
  )


def _extract_fulltext_patent_tables(
    all_tables: List[Any],
) -> Tuple[
    Dict[str, List[str]],
    Dict[str, List[Tuple[str, str]]],
    Dict[str, str],
    Dict[str, float],
]:
  """Extracts Example assay values, column headers, IUPAC names, and MS [M+H]+ from HTML <table> elements."""
  example_to_assays: Dict[str, List[str]] = {}
  example_to_assay_cols: Dict[str, List[Tuple[str, str]]] = {}
  example_to_iupac: Dict[str, str] = {}
  example_to_mh_mass: Dict[str, float] = {}

  carried_header_labels: List[str] = []
  carried_col1_is_example: bool = False
  carried_is_synth_table: bool = False

  for tbl in all_tables:
    rows = tbl.find_all("tr")
    col1_is_example_header = False
    header_labels: List[str] = []
    is_synth_table = False
    prev_h_tag = tbl.find_previous("heading")
    prev_h_txt = prev_h_tag.get_text(" ", strip=True) if prev_h_tag else ""
    is_prep_section_table = bool(
        re.match(
            r"^\s*(?:Preparations?|Intermediates?)\b",
            prev_h_txt,
            re.IGNORECASE,
        )
        and not re.search(r"\bExamples?\b", prev_h_txt, re.IGNORECASE)
    )

    # Determine how many leading rows precede the first Example data row
    # (up to 8 rows for multi-row headers).
    first_data_r = len(rows)
    for r_idx, hr in enumerate(rows[:8]):
      hr_tags = hr.find_all(["td", "th"])
      h_cells_id = [_cell_text_no_sup(c) for c in hr_tags]
      if len(h_cells_id) >= 2 and re.match(
          r"^(?:Example|Ex\.?)?\s*\d+[a-zA-Z]?$",
          h_cells_id[0].strip(),
          re.IGNORECASE,
      ):
        first_data_r = r_idx
        break
      if len(h_cells_id) >= 3 and re.match(
          r"^(?:Example|Ex\.?)\s*\d+[a-zA-Z]?$",
          h_cells_id[1].strip(),
          re.IGNORECASE,
      ):
        first_data_r = r_idx
        break

    header_row_slices = [
        [c.get_text(" ", strip=True) for c in hr.find_all(["td", "th"])]
        for hr in rows[: min(first_data_r, 8)]
    ]
    header_row_slices = [
        rc
        for rc in header_row_slices
        if len(rc) >= 2 and any(x.strip() for x in rc)
    ]
    if header_row_slices:
      max_h_cols = max(len(rc) for rc in header_row_slices)
      merged_hcols: List[str] = []
      for c_i in range(max_h_cols):
        parts = [
            rc[c_i].strip()
            for rc in header_row_slices
            if c_i < len(rc)
            and rc[c_i].strip()
            and not re.match(r"^TABLE\s+\d+$", rc[c_i].strip(), re.IGNORECASE)
        ]
        merged_hcols.append(" ".join(parts).strip())
      hcells_low = [c.lower() for c in merged_hcols]
      joined_hlow = " ".join(hcells_low)
      # Skip Preparation / synthetic intermediate tables and bibliographic
      # metadata tables.
      if (
          is_prep_section_table
          or any(
              _PREP_OR_INTERMEDIATE_HEADER_RE.search(hc)
              for hc in hcells_low[:2]
          )
          or any(
              k in joined_hlow
              for k in (
                  "application number",
                  "priority date",
                  "publication number",
                  "filing date",
              )
          )
      ):
        carried_header_labels = []
        carried_col1_is_example = False
        carried_is_synth_table = False
        continue

      has_bio_hdr = any(
          any(
              k in hc
              for k in (
                  "ic50",
                  "ic 50",
                  "ec50",
                  "ec 50",
                  "ki",
                  "kd",
                  "ed50",
                  "gi50",
                  "elisa",
                  "assay",
                  "enzyme",
                  "inhibition",
                  "activity",
              )
          )
          for hc in hcells_low
      )
      has_synth_hdr = (not has_bio_hdr) and any(
          any(
              k in hc
              for k in ("name", "iupac", "nmr", "yield", "form", "data", "ms")
          )
          for hc in hcells_low
      )
      if (
          has_bio_hdr
          or has_synth_hdr
          or any(
              any(k in hc for k in ("example", "ex.", "compound", "cpd", "no."))
              for hc in hcells_low
          )
      ):
        header_labels = merged_hcols
        is_synth_table = has_synth_hdr
        if (
            len(hcells_low) >= 2
            and any(
                k in hcells_low[1]
                for k in ("example", "ex.", "compound", "cpd", "no.", "id")
            )
            and not any(k in hcells_low[0] for k in ("example", "ex."))
        ):
          col1_is_example_header = True
        carried_header_labels = list(header_labels)
        carried_col1_is_example = col1_is_example_header
        carried_is_synth_table = is_synth_table

    if is_prep_section_table:
      carried_header_labels = []
      carried_col1_is_example = False
      carried_is_synth_table = False
      continue

    if not header_labels and carried_header_labels and first_data_r == 0:
      header_labels = list(carried_header_labels)
      col1_is_example_header = carried_col1_is_example
      is_synth_table = carried_is_synth_table

    for r in rows:
      r_tags = r.find_all(["td", "th"])
      cells = [c.get_text(" ", strip=True) for c in r_tags]
      if len(cells) < 2:
        continue
      c0_nosup = _cell_text_no_sup(r_tags[0]).strip()
      ex_id = None
      val_start = 1
      m0 = re.match(
          r"^(?:Example|Ex\.?)?\s*(\d+[a-zA-Z]?)$",
          c0_nosup,
          re.IGNORECASE,
      )
      if m0:
        ex_id = m0.group(1)
        val_start = 1
      elif len(cells) >= 3:
        c1_nosup = _cell_text_no_sup(r_tags[1]).strip()
        m1_pat = (
            r"^(?:Example|Ex\.?)?\s*(\d+[a-zA-Z]?)$"
            if col1_is_example_header
            else r"^(?:Example|Ex\.?)\s*(\d+[a-zA-Z]?)$"
        )
        m1 = re.match(m1_pat, c1_nosup, re.IGNORECASE)
        if m1:
          ex_id = m1.group(1)
          val_start = 2
      if not ex_id:
        continue

      # If this table (or continuation table) is a synthesis/characterization
      # table, extract IUPAC name and MS [M+H]+ without polluting
      # example_to_assays.
      row_has_table_iupac = any(
          _looks_like_table_iupac(c_txt) for c_txt in cells[val_start:]
      )
      if is_synth_table or (not header_labels and row_has_table_iupac):
        for c_pos in range(val_start, len(cells)):
          c_txt = cells[c_pos].strip()
          if not c_txt:
            continue
          hdr_low = (
              header_labels[c_pos].lower() if c_pos < len(header_labels) else ""
          )
          if (
              (
                  "name" in hdr_low
                  or "iupac" in hdr_low
                  or _looks_like_table_iupac(c_txt)
              )
              and _looks_like_table_iupac(c_txt)
              and ex_id not in example_to_iupac
          ):
            cleaned_tbl_iup = re.sub(r"-\s+", "-", c_txt)
            cleaned_tbl_iup = re.sub(
                r"\s*\((?:Compound\s+)?\d+[a-zA-Z]?\)\s*$", "", cleaned_tbl_iup
            ).strip()
            example_to_iupac[ex_id] = cleaned_tbl_iup
          m_ms = re.findall(
              r"\b(\d{3,4}(?:\.\d+)?)\s*\[M\s*\+?\s*(?:H|1)\]",
              c_txt,
              re.IGNORECASE,
          )
          if m_ms and ex_id not in example_to_mh_mass:
            try:
              example_to_mh_mass[ex_id] = float(m_ms[-1])
            except Exception:
              pass
        continue

      vals = []
      named_vals: List[Tuple[str, str]] = []
      i = val_start
      col_idx = 0
      n_hdr_cols = max(0, len(header_labels) - val_start)
      row_has_unit_cells = any(
          c.strip().lower() in ("nm", "um", "µm", "μm", "pm", "mm", "%")
          for c in cells[val_start:]
      )
      is_paired_unit_table = row_has_unit_cells or (
          n_hdr_cols >= 1 and (len(cells) - val_start) == 2 * n_hdr_cols
      )
      while i < len(cells):
        c_val = cells[i].strip()
        has_separate_unit = i + 1 < len(cells) and cells[
            i + 1
        ].strip().lower() in ("nm", "um", "µm", "μm", "pm", "mm", "%")
        has_empty_unit_pair = (
            not c_val
            and is_paired_unit_table
            and i + 1 < len(cells)
            and not cells[i + 1].strip()
        )
        hdr_pos = (
            val_start + col_idx
            if (
                has_separate_unit
                or has_empty_unit_pair
                or len(header_labels) < len(cells)
            )
            else i
        )
        hdr_lbl = (
            header_labels[hdr_pos].strip()
            if hdr_pos < len(header_labels) and header_labels[hdr_pos].strip()
            else f"col_{col_idx + 1}"
        )
        if has_separate_unit:
          c_val = f"{c_val} {cells[i + 1].strip()}"
          i += 2
        elif has_empty_unit_pair:
          i += 2
        else:
          i += 1
        col_idx += 1
        vals.append(c_val)
        if c_val:
          named_vals.append((hdr_lbl, c_val))
      if any(v.strip() for v in vals):
        example_to_assays.setdefault(ex_id, []).extend(vals)
        example_to_assay_cols.setdefault(ex_id, []).extend(named_vals)

  return (
      example_to_assays,
      example_to_assay_cols,
      example_to_iupac,
      example_to_mh_mass,
  )


def enrich_records_with_fulltext_patent_html(
    patent_dict: Dict[str, Any], records: List[Dict[str, Any]]
) -> None:
  """Fetches Google Patents full-text HTML to extract Table 1 bioactivities and resolve Example IUPAC headings."""
  candidates = [patent_dict.get("compact", "")]
  for st in patent_dict.get("search_terms", []):
    clean_st = st.replace("-", "").replace(" ", "").replace("/", "")
    if clean_st and clean_st not in candidates:
      candidates.append(clean_st)

  html = None
  used_id = ""
  tried = []
  fetched_but_no_tables = []
  for cand in candidates:
    if not cand or cand.isdigit():
      continue
    tried.append(cand)
    url = f"https://patents.google.com/patent/{urllib.parse.quote(cand)}/en"
    page = _http_get_text(url, timeout=12)
    if not page:
      continue
    if (
        "<table" in page
        or "<heading" in page
        or 'class="description-paragraph"' in page
    ):
      html = page
      used_id = cand
      break
    # The page came back but carries neither tables nor description headings.
    # Remember it so a successful fetch is never silently discarded by a
    # subsequent failing candidate.
    fetched_but_no_tables.append(cand)
    if html is None:
      html = page
      used_id = cand

  if not html:
    fallback_dates = _extract_patent_dates_from_soup(None, patent_dict)
    patent_dict["_patent_dates"] = fallback_dates
    _apply_patent_dates_to_records(records, fallback_dates)
    _warn(
        f"{patent_dict.get('compact', '?')}: no full-text HTML retrieved from"
        f" Google Patents (tried {tried or 'no usable identifiers'}). No"
        " example-table compounds will be extracted for this patent."
    )
    return
  if fetched_but_no_tables and used_id in fetched_but_no_tables:
    _warn(
        f"{patent_dict.get('compact', '?')}: full-text HTML for {used_id}"
        " contains no <table> or <heading> elements. Example tables in this"
        " patent are most likely images, which this extractor cannot read."
    )

  soup = BeautifulSoup(html, "html.parser")
  patent_dates = _extract_patent_dates_from_soup(soup, patent_dict)
  patent_dict["_patent_dates"] = patent_dates
  bin_legend = _extract_patent_bin_legend(html)

  # 1. Extract Table rows indexed by Example number (with multi-row column
  # headers & continuation table carry-over).
  _all_tables = soup.find_all("table")
  _table_count = len(_all_tables)
  (
      example_to_assays,
      example_to_assay_cols,
      example_to_iupac,
      example_to_mh_mass,
  ) = _extract_fulltext_patent_tables(_all_tables)
  example_to_peak: Dict[str, int] = {}
  heading_iupac_ex_ids: Set[str] = set()

  # 2. Extract Example -> IUPAC, Peak label, LCMS [M+H]+, and inline Ki/IC50
  # from <heading> slices (with fallback to <div class="description-paragraph">
  # / <p> Example headers for WIPO WO patents).
  heading_matches = list(
      re.finditer(
          r"<heading[^>]*>(.*?)</heading>", html, re.IGNORECASE | re.DOTALL
      )
  )
  _heading_count = len(heading_matches)
  headings = [
      BeautifulSoup(hm.group(1), "html.parser").get_text(" ", strip=True)
      for hm in heading_matches
  ]
  if not any(re.match(r"^\s*Examples?\b", h, re.IGNORECASE) for h in headings):
    para_matches = list(
        re.finditer(
            r"<(?:div[^>]*class=\"description-paragraph\"[^>]*|p[^>]*?)>\s*(?:<(?:b|strong)>\s*)?(Examples?\s+\d+[a-zA-Z]?\s*[:.;].*?)</(?:div|p)>",
            html,
            re.IGNORECASE | re.DOTALL,
        )
    )
    if para_matches:
      heading_matches = para_matches
      _heading_count = len(heading_matches)
      headings = [
          BeautifulSoup(hm.group(1), "html.parser").get_text(" ", strip=True)
          for hm in heading_matches
      ]

  # Note: _is_iupac_heading and _clean_iupac_heading rely on Google Patents'
  # HTML heading conventions for USPTO/EPO/WIPO full-text specifications.
  # End-to-end extraction is exercised in tests/patent_fetcher_test.py and
  # TEST.md (Scenario 4) so markup drift is caught early.
  def _is_iupac_heading(txt: str) -> bool:
    if len(txt) <= 18:
      return False
    low = txt.lower().strip()
    if low.startswith((
        "example",
        "peak",
        "atropisomer",
        "step",
        "general procedure",
        "general scheme",
        "biological",
        "assay",
        "table",
        "field of",
        "background",
        "summary",
        "detailed description",
        "claims",
        "select starting",
    )):
      return False
    return True

  def _clean_iupac_heading(txt: str) -> str:
    cleaned = re.sub(
        r"^(?:Preparation of|Synthesis of)\s+", "", txt, flags=re.IGNORECASE
    ).strip()
    cleaned = re.sub(
        r"\s*\((?:Compound\s+)?\d+[a-zA-Z]?\)\s*$", "", cleaned
    ).strip()
    cleaned = re.sub(r"\s+H\s*\d+\s*$", "", cleaned).strip()
    cleaned = re.sub(r"\.\s*$", "", cleaned).strip()
    return cleaned

  def _extract_final_step_lcms(body_text: str) -> Optional[float]:
    step_splits = re.split(
        r"\bStep\s+\d+\s*[:.]", body_text, flags=re.IGNORECASE
    )
    target_txt = step_splits[-1] if len(step_splits) > 1 else body_text
    m_lcms = re.findall(
        r"(?:LCMS|LC-MS|MS|APCI|ESI|HRMS)[^.;]{0,55}?\b(\d{3,4}(?:\.\d+)?)\s*\[M\s*\+?\s*(?:H|1)\]",
        target_txt,
        re.IGNORECASE,
    )
    if m_lcms:
      try:
        return float(m_lcms[-1])
      except Exception:
        return None
    return None

  example_to_prose: Dict[str, str] = {}
  example_in_group: Set[str] = set()
  current_group_iupac = ""
  current_group_iupac_pair: Optional[Tuple[str, str]] = None
  current_group_tokens: Set[str] = set()
  for i, h in enumerate(headings):
    hm = heading_matches[i]
    body_start = hm.end()
    body_end = (
        heading_matches[i + 1].start()
        if i + 1 < len(heading_matches)
        else len(html)
    )
    if re.match(
        r"^(?:Biological|General Procedure|General Scheme|Claims|Assay)",
        h,
        re.IGNORECASE,
    ):
      current_group_iupac = ""
      current_group_iupac_pair = None
      current_group_tokens = set()
      continue

    if (
        re.match(r"^(?:Preparation of|Synthesis of)\s+", h, re.IGNORECASE)
        and len(h) > 25
    ):
      cleaned_prep = _clean_iupac_heading(h)
      current_group_tokens = set()
      # Check for combined heading: "Preparation of <IUPAC_A> (Example N) and <IUPAC_B>"
      m_comb = re.match(
          r"^(.+?)\s*\(\s*Examples?\s+(\d+[a-zA-Z]?)\s*\)\s+and\s+(.+)$",
          cleaned_prep,
          re.IGNORECASE,
      )
      if m_comb:
        ex_comb = m_comb.group(2).strip()
        iup_a = _clean_iupac_heading(m_comb.group(1))
        iup_b = _clean_iupac_heading(m_comb.group(3))
        if iup_a and ex_comb not in example_to_iupac:
          example_to_iupac[ex_comb] = iup_a
          heading_iupac_ex_ids.add(ex_comb)
        comb_body = BeautifulSoup(
            html[body_start:body_end], "html.parser"
        ).get_text(" ", strip=True)
        mh_comb = _extract_final_step_lcms(comb_body)
        if mh_comb is not None and ex_comb not in example_to_mh_mass:
          example_to_mh_mass[ex_comb] = mh_comb
        current_group_iupac = iup_b
        current_group_iupac_pair = None
        continue

      # Check for slash-separated dual IUPAC heading: "<IUPAC_1>/<IUPAC_2>"
      slash_parts = [
          p.strip() for p in cleaned_prep.split("/") if len(p.strip()) > 25
      ]
      if len(slash_parts) == 2 and "/" in cleaned_prep:
        current_group_iupac = slash_parts[0]
        current_group_iupac_pair = (slash_parts[0], slash_parts[1])
      else:
        current_group_iupac = cleaned_prep
        current_group_iupac_pair = None
      continue

    if re.match(r"^\s*Examples?\b", h, re.IGNORECASE):
      # Extract only Example tokens from the Example header prefix (before any colon/IUPAC)
      ex_prefix = re.split(r"[:;]|Preparation of|Synthesis of", h, maxsplit=1)[
          0
      ]
      ex_tokens = re.findall(r"\b(\d+[a-zA-Z]?)\b", ex_prefix)
      if not ex_tokens:
        continue
      if len(ex_tokens) > 1:
        current_group_tokens = set(ex_tokens)
        example_in_group.update(ex_tokens)
      elif ex_tokens[0] in current_group_tokens:
        example_in_group.add(ex_tokens[0])
      already_defined_single = (
          len(ex_tokens) == 1
          and ex_tokens[0] in example_to_iupac
          and ex_tokens[0] not in current_group_tokens
      )
      iupac = ""
      m_self = re.search(
          r"(?:Preparation of|Synthesis of|[:])\s*(.+)", h, re.IGNORECASE
      )
      if m_self and _is_iupac_heading(m_self.group(1)):
        self_rhs = m_self.group(1).strip()
        # If a WIPO paragraph heading wraps across two
        # <div class="description-paragraph"> / <p> tags ending with a trailing
        # hyphen or unclosed bracket, append the continuation paragraph.
        if (
            self_rhs.endswith("-")
            or self_rhs.count("[") > self_rhs.count("]")
            or self_rhs.count("(") > self_rhs.count(")")
            or self_rhs.count("{") > self_rhs.count("}")
        ):
          m_next_p = re.match(
              r"\s*<(?:div[^>]*class=\"description-paragraph\"[^>]*|p[^>]*?)>(.*?)</(?:div|p)>",
              html[hm.end() : body_end],
              re.IGNORECASE | re.DOTALL,
          )
          if m_next_p and "<img" not in m_next_p.group(1).lower():
            cont_txt = BeautifulSoup(m_next_p.group(1), "html.parser").get_text(
                " ", strip=True
            )
            if (
                cont_txt
                and len(cont_txt) < 100
                and not re.match(
                    r"^(?:Examples?|Step|To\s+a|A\s+mixture|A\s+solution|The\s+title|1\s*H\s+NMR|LCMS|MS)\b",
                    cont_txt,
                    re.IGNORECASE,
                )
            ):
              self_rhs = f"{self_rhs}{cont_txt}"
              body_start = hm.end() + m_next_p.end()
        # For single-paragraph WIPO entries ("Example N: <IUPAC>. To a
        # solution..." or "... <Target> Ki: 0.01 uM"), trim at the sentence
        # boundary preceding experimental prose or inline bioactivity.
        self_iup = re.split(
            r"\.\s+(?:To\s+a|A\s+mixture|A\s+solution|The\s+title|Prepared\s+according|Step\s+\d|1\s*H\s+NMR|LCMS|MS\b|[A-Za-z0-9_-]+\s+(?:Ki|IC\s*50))|\s+[A-Za-z0-9_]+(?:-\s*[A-Za-z0-9_]+)?\s+(?:Ki|IC\s*50|EC\s*50)\s*[:=]",
            self_rhs,
            maxsplit=1,
            flags=re.IGNORECASE,
        )[0]
        iupac = _clean_iupac_heading(self_iup)
      elif (
          i + 1 < len(headings)
          and _is_iupac_heading(headings[i + 1])
          and not re.match(
              r"^(?:Preparation of|Synthesis of)\s+",
              headings[i + 1],
              re.IGNORECASE,
          )
      ):
        # Forward heading convention: <heading>Example 1</heading> ->
        # <heading><IUPAC></heading>. In multi-step syntheses where a
        # <chemistry> reaction scheme or Step 1 block precedes the first IUPAC
        # sub-heading (or the first sub-heading is a synthetic intermediate),
        # collect all IUPAC sub-headings before the next Example/Table and take
        # the final product.
        gap_before_first = html[hm.end() : heading_matches[i + 1].start()]
        is_multistep_subheadings = (
            "<chemistry" in gap_before_first.lower()
            or "step " in gap_before_first.lower()
            or bool(
                re.search(
                    r"(?:dioxaborolan|boronic|stannane|silyl|^tert-butyl\b|carbamate\b)",
                    headings[i + 1],
                    re.IGNORECASE,
                )
            )
        )
        chosen_idx = i + 1
        if is_multistep_subheadings:
          sub_indices = [i + 1]
          k_sub = i + 2
          while k_sub < len(headings):
            hk = headings[k_sub]
            if re.match(
                r"^\s*(?:Examples?|Biological|General Procedure|General"
                r" Scheme|Claims|Assay|Table)\b",
                hk,
                re.IGNORECASE,
            ) or re.match(
                r"^(?:Preparation of|Synthesis of)\s+", hk, re.IGNORECASE
            ):
              break
            seg_between = html[
                heading_matches[k_sub - 1]
                .end() : heading_matches[k_sub]
                .start()
            ]
            if "<table" in seg_between.lower():
              break
            if _is_iupac_heading(hk):
              sub_indices.append(k_sub)
            k_sub += 1
          if len(sub_indices) > 1:
            non_interm = [
                idx_s
                for idx_s in sub_indices
                if not re.search(
                    r"(?:dioxaborolan|boronic|stannane|silyl|^tert-butyl\b|carbamate\b)",
                    headings[idx_s],
                    re.IGNORECASE,
                )
            ]
            chosen_idx = non_interm[-1] if non_interm else sub_indices[-1]
        iupac = _clean_iupac_heading(headings[chosen_idx])
        body_start = heading_matches[chosen_idx].end()
        body_end = (
            heading_matches[chosen_idx + 1].start()
            if chosen_idx + 1 < len(heading_matches)
            else len(html)
        )
      elif (
          i + 1 < len(headings)
          and re.match(
              r"^(?:Preparation of|Synthesis of)\s+",
              headings[i + 1],
              re.IGNORECASE,
          )
          and len(
              BeautifulSoup(
                  html[hm.end() : heading_matches[i + 1].start()], "html.parser"
              ).get_text(strip=True)
          )
          < 20
      ):
        # Immediate sub-heading: <heading>Example N</heading> (0-char gap) -> <heading>Preparation of <IUPAC></heading>
        iupac = _clean_iupac_heading(headings[i + 1])
        body_start = heading_matches[i + 1].end()
        body_end = (
            heading_matches[i + 2].start()
            if i + 2 < len(heading_matches)
            else len(html)
        )
      elif current_group_iupac:
        # Preceding group-header convention: <heading>Preparation of <IUPAC></heading> -> <heading>Example N</heading>
        iupac = current_group_iupac

      if len(ex_tokens) == 2 and current_group_iupac_pair is not None:
        if ex_tokens[0] not in example_to_iupac:
          example_to_iupac[ex_tokens[0]] = current_group_iupac_pair[0]
          heading_iupac_ex_ids.add(ex_tokens[0])
        if ex_tokens[1] not in example_to_iupac:
          example_to_iupac[ex_tokens[1]] = current_group_iupac_pair[1]
          heading_iupac_ex_ids.add(ex_tokens[1])
        current_group_iupac = current_group_iupac_pair[0]
        current_group_iupac_pair = None
      elif len(ex_tokens) > 1 and iupac:
        current_group_iupac = iupac

      # Check if immediately followed by <heading>[Atropisomer] Peak 1</heading> or <heading>[Atropisomer] Peak 2</heading>
      if len(ex_tokens) == 1 and i + 1 < len(headings):
        m_pk = re.match(
            r"^(?:Atropisomer\s+|Enantiomer\s+|Isomer\s+)?Peak\s*([12])\b",
            headings[i + 1],
            re.IGNORECASE,
        )
        if m_pk:
          if not already_defined_single:
            example_to_peak[ex_tokens[0]] = int(m_pk.group(1))
          body_start = heading_matches[i + 1].end()
          body_end = (
              heading_matches[i + 2].start()
              if i + 2 < len(heading_matches)
              else len(html)
          )

      body_txt = (
          h
          + " "
          + BeautifulSoup(html[body_start:body_end], "html.parser").get_text(
              " ", strip=True
          )
      )
      if not already_defined_single:
        for ex in ex_tokens:
          if len(ex_tokens) == 1 or ex not in example_to_prose:
            example_to_prose[ex] = body_txt
      # Extract LCMS [M+H]+ m/z from the final step reported in prose
      mh_val = _extract_final_step_lcms(body_txt)
      if mh_val is not None and not already_defined_single:
        for ex in ex_tokens:
          if len(ex_tokens) == 1 or ex not in example_to_mh_mass:
            example_to_mh_mass[ex] = mh_val

      # Extract inline prose Ki / IC50 / EC50 / Kd bioactivities when not already in Table
      for km in re.finditer(
          r"(?:([A-Za-z0-9_]+(?:-\s*[A-Za-z0-9_]+)?)\s+)?(Ki|IC\s*50|EC\s*50|Kd)\s*[:=]\s*([<>~]?\s*[0-9]+(?:\.[0-9]+)?\s*(?:pM|nM|uM|µM|μM|mM))",
          body_txt,
          re.IGNORECASE,
      ):
        raw_tgt_tok = re.sub(r"-\s+", "-", (km.group(1) or "").strip())
        act_kind_tok = re.sub(r"\s+", "", km.group(2) or "Ki")
        val_with_unit = km.group(3).strip()
        hdr_prose = (
            f"{raw_tgt_tok} {act_kind_tok}".strip()
            if raw_tgt_tok
            and raw_tgt_tok.lower()
            not in (
                "the",
                "a",
                "an",
                "with",
                "had",
                "showed",
                "exhibited",
                "of",
            )
            else act_kind_tok
        )
        for ex in ex_tokens:
          if ex not in example_to_assays:
            example_to_assays.setdefault(ex, []).append(val_with_unit)
            example_to_assay_cols.setdefault(ex, []).append(
                (hdr_prose, val_with_unit)
            )

      if iupac:
        for ex in ex_tokens:
          if ex not in example_to_iupac:
            example_to_iupac[ex] = iupac
            heading_iupac_ex_ids.add(ex)

  if not example_to_assays and not example_to_iupac:
    _warn(
        f"{patent_dict.get('compact', '?')}: parsed"
        f" {_table_count} <table> and {_heading_count} <heading> elements but"
        " recovered 0 example rows and 0 example headings. The example"
        " identifiers in this patent do not match the expected"
        " 'Example N' / 'N' cell format, so no example-table compounds will"
        " be extracted."
    )
    return

  print(
      f"[patent_fetcher] Full-text HTML ({used_id}): found"
      f" {len(example_to_assays)} assay rows & {len(example_to_iupac)}"
      " Example IUPAC headings",
      file=sys.stderr,
      flush=True,
  )

  def _parse_pactivity(val_str: str) -> Tuple[Optional[float], Optional[float]]:
    if not val_str:
      return None, None
    s_clean = val_str.strip()
    if (
        s_clean.upper() in _INACTIVE_OR_MISSING_TOKENS
        or s_clean in bin_legend
        or _is_discrete_potency_token(s_clean)
    ):
      return None, None
    if re.search(r"\d\s*[-–—]\s*\d|\bto\b", s_clean, re.IGNORECASE):
      return None, None
    m = re.search(
        r"([0-9]+(?:\.[0-9]+)?)\s*(pM|nM|uM|µM|μM|mM)?",
        s_clean,
        re.IGNORECASE,
    )
    if not m:
      return None, None
    num = float(m.group(1))
    unit = (m.group(2) or "nM").lower()
    factor = {
        "pm": 1e-12,
        "nm": 1e-9,
        "um": 1e-6,
        "µm": 1e-6,
        "μm": 1e-6,
        "mm": 1e-3,
    }.get(unit, 1e-9)
    val_m = num * factor
    if val_m <= 0:
      return None, None
    return round(num * (factor / 1e-9), 4), round(-math.log10(val_m), 3)

  def norm_iupac(t: str) -> str:
    return re.sub(r"[^a-zA-Z0-9]", "", str(t)).lower()

  # Index existing records by explicit PubChem Example title ("US..., Ex 343"),
  # normalized IUPAC name, full InChIKey, and 14-char connectivity InChIKey.
  pubchem_ex_to_recs: Dict[str, List[Dict[str, Any]]] = {}
  iupac_to_recs: Dict[str, List[Dict[str, Any]]] = {}
  ik_to_recs: Dict[str, List[Dict[str, Any]]] = {}
  conn_to_recs: Dict[str, List[Dict[str, Any]]] = {}
  norm_iupac_to_smi: Dict[str, str] = {}
  for r in records:
    cname = str(r.get("compound_name") or "").strip()
    if cname and len(cname) < 45:
      m_pex = re.search(
          r"\b(?:Example|Ex\.?|Compound|Cmpd\.?)\s+(\d+[a-zA-Z]?)\b",
          cname,
          re.IGNORECASE,
      )
      if m_pex:
        ex_k = m_pex.group(1)
        r["_pubchem_ex_id"] = ex_k
        pubchem_ex_to_recs.setdefault(ex_k, []).append(r)
    if r.get("iupac_name"):
      key = norm_iupac(r["iupac_name"])
      if len(key) >= 15:
        iupac_to_recs.setdefault(key, []).append(r)
        if key not in norm_iupac_to_smi and r.get("canonical_smiles"):
          norm_iupac_to_smi[key] = str(r["canonical_smiles"])
    mol_r = r.get("rdkit_mol")
    if mol_r is not None and not r.get("inchikey"):
      try:
        r["inchikey"] = Chem.MolToInchiKey(mol_r)
      except Exception:
        pass
    if r.get("inchikey"):
      ik_str = str(r["inchikey"])
      ik_to_recs.setdefault(ik_str, []).append(r)
      conn_block = ik_str.split("-")[0]
      if conn_block:
        conn_to_recs.setdefault(conn_block, []).append(r)

  # When PubChem indexes patent Example numbers in Title for most records, a few
  # named drug/INN records (e.g. Lirafugratinib, CID 155309155) have their INN
  # as Title and store "US11780845, Example 343" in PubChem Synonyms.
  if len(pubchem_ex_to_recs) >= 10:
    short_named_recs = [
        r
        for r in records
        if not r.get("_pubchem_ex_id")
        and r.get("pubchem_cid")
        and len(str(r.get("compound_name") or "")) < 35
    ][:300]
    if short_named_recs:
      cid_to_rec_map = {str(r["pubchem_cid"]): r for r in short_named_recs}
      pat_base_num = re.sub(
          r"[A-Z]\d?$", "", str(patent_dict.get("compact") or used_id)
      )
      cids_list = list(cid_to_rec_map.keys())
      for start_idx in range(0, len(cids_list), 100):
        chunk_cids = ",".join(cids_list[start_idx : start_idx + 100])
        syn_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{chunk_cids}/synonyms/JSON"
        syn_data = _http_get_json(syn_url, timeout=10)
        if not syn_data:
          continue
        info_list = syn_data.get("InformationList", {}).get("Information", [])
        for info_item in info_list:
          cid_str = str(info_item.get("CID", ""))
          rec_target = cid_to_rec_map.get(cid_str)
          if rec_target is None:
            continue
          for syn_str in info_item.get("Synonym", []):
            if pat_base_num and pat_base_num.lower() not in syn_str.lower():
              continue
            m_syn_ex = re.search(
                r"\b(?:Example|Ex\.?|Compound|Cmpd\.?)\s+(\d+[a-zA-Z]?)\b",
                syn_str,
                re.IGNORECASE,
            )
            if m_syn_ex:
              ex_syn = m_syn_ex.group(1)
              if not rec_target.get("_pubchem_ex_id"):
                rec_target["_pubchem_ex_id"] = ex_syn
              if rec_target not in pubchem_ex_to_recs.get(ex_syn, []):
                pubchem_ex_to_recs.setdefault(ex_syn, []).append(rec_target)

  # Detect if synthetic <heading>Example N</heading> numbering is decoupled from
  # Table 1 compound numbering (e.g., US11780845B2 has 54 synthetic procedure
  # Examples 1..54 while Table 1 numbers compounds 1..1496 across Tables 2..63).
  if len(pubchem_ex_to_recs) >= 20 and heading_iupac_ex_ids:
    matched_same = 0
    matched_diff = 0
    for ex_h in heading_iupac_ex_ids:
      nk_h = norm_iupac(example_to_iupac.get(ex_h, ""))
      if nk_h and nk_h in iupac_to_recs:
        for cand_h in iupac_to_recs[nk_h]:
          pex_h = cand_h.get("_pubchem_ex_id")
          if pex_h:
            if pex_h == ex_h:
              matched_same += 1
            else:
              matched_diff += 1
    if matched_diff >= 3 and matched_diff > matched_same:
      for ex_h in list(heading_iupac_ex_ids):
        if ex_h in pubchem_ex_to_recs:
          example_to_iupac.pop(ex_h, None)
          example_to_mh_mass.pop(ex_h, None)

  # 3. Resolve remaining IUPAC names via OPSIN in parallel for examples not directly in PubChem
  unmatched_iupacs = {}
  for ex_id, iup in example_to_iupac.items():
    if ex_id not in pubchem_ex_to_recs and norm_iupac(iup) not in iupac_to_recs:
      unmatched_iupacs[ex_id] = iup

  _VALID_HYPHEN_PREFIXES = frozenset({
      "tert",
      "sec",
      "iso",
      "neo",
      "cis",
      "trans",
      "endo",
      "exo",
      "syn",
      "anti",
      "meso",
      "gem",
      "vic",
      "sym",
      "unsym",
      "ortho",
      "meta",
      "para",
      "alpha",
      "beta",
      "gamma",
      "delta",
      "epsilon",
      "zeta",
      "eta",
      "theta",
      "kappa",
      "lambda",
      "mu",
      "omega",
      "catena",
      "closo",
      "nido",
      "arachno",
      "hypho",
      "klado",
      "retro",
      "nor",
      "homo",
      "seco",
      "abeo",
      "friedo",
      "cyclo",
      "spiro",
      "bicyclo",
      "tricyclo",
      "tetracyclo",
      "poly",
      "copoly",
      "mono",
      "bis",
      "tris",
  })

  def _repair_iupac_ocr_typos(raw_name: str) -> List[str]:
    """Generates candidate IUPAC strings when an initial string fails OPSIN due to OCR locant/hyphen/brace typos."""
    candidates: List[str] = []
    # 0. Replace WIPO '⁇' OCR placeholders for curly braces '{' and '}'
    c_br = raw_name.strip()
    c_br = re.sub(r"\s+H\s*\d+\s*$", "", c_br)
    if "⁇" in c_br:
      parts_q = c_br.split("⁇")
      rebuilt_q = []
      open_b = True
      for idx_q, pq in enumerate(parts_q):
        rebuilt_q.append(pq)
        if idx_q < len(parts_q) - 1:
          rebuilt_q.append("{" if open_b else "}")
          open_b = not open_b
      c_br = "".join(rebuilt_q)

    # 1. General character-level & glyph-fusion OCR confusables
    # Opening '[' misread as 'f' after numeric locant ("3-f2-" -> "3-[2-")
    c_br = re.sub(r"\b(\d+)-f\s*(?=[\d\(])", r"\1-[", c_br)
    # 'y' misread as 'v' before 'l' ("vl" -> "yl") or after 'x' ("xv" -> "xy")
    c_br = re.sub(r"vl", "yl", c_br)
    c_br = re.sub(r"xv", "xy", c_br)
    # In unclosed '(' expressions, 'l)' before an alkoxy suffix fused into 'π' or 'n'
    if c_br.count("(") > c_br.count(")"):
      c_br = re.sub(
          r"(phenv|pheny|[a-z]yl)[πn](?=(?:meth|eth|prop|but)?oxy)",
          r"\1l)",
          c_br,
          flags=re.IGNORECASE,
      )
      c_br = re.sub(r"vl", "yl", c_br)
    # Greek OCR confusables & 'l)' fused into 'O' after 'y' ("yOmethyl" -> "yl)methyl")
    c_br = c_br.replace("θ", "e").replace("β", "e").replace("π", "n")
    c_br = re.sub(r"yO(?=[a-z])", "yl)", c_br)
    # 'l' misread as '1' after 'y' ("y1" -> "yl" when not part of a closing bracket)
    if c_br.count("[") == c_br.count("]"):
      c_br = re.sub(r"y1(?!\d)", "yl", c_br)
    # Closing '}' or ']' misread as 'j', 'J', '1', or 'l' after '-yl' / '-oxy' when braces/brackets are unclosed
    if c_br.count("{") > c_br.count("}"):
      c_br = re.sub(r"(yl|oxy)[jJ1l](?=[-\)])", r"\1}", c_br)
    if c_br.count("[") > c_br.count("]"):
      c_br = re.sub(
          r"(yl|oxy)[jJ1l](?!ate\b|ic\b)(?=[A-Za-z\-\)]|$)", r"\1]", c_br
      )
    c_br = re.sub(r"y1(?!\d)", "yl", c_br)
    # Sans-serif glyph confusables: 'ri' <-> 'h', 'chior' -> 'chlor', 'pvr'/'dvr' -> 'pyr', 'f luor' -> 'fluor'
    c_br = re.sub(r"(py|pipe)h(d|m)", r"\1ri\2", c_br, flags=re.IGNORECASE)
    c_br = re.sub(r"prien", "phen", c_br, flags=re.IGNORECASE)
    c_br = re.sub(r"chior", "chlor", c_br, flags=re.IGNORECASE)
    c_br = re.sub(r"f\s+luor", "fluor", c_br, flags=re.IGNORECASE)
    c_br = re.sub(r"[pd]vr", "pyr", c_br, flags=re.IGNORECASE)
    # Triple repeated consonant ("pyrrrol" -> "pyrrol")
    c_br = re.sub(r"([a-z])\1\1", r"\1\1", c_br, flags=re.IGNORECASE)
    # Locant punctuation & fusion typos ("N.N-" -> "N,N-", "2,-(" -> "2'-(", "[3,3 l ]" -> "[3,3']", "1,2,33-" -> "1,2,3-")
    c_br = re.sub(r"\b(N|\d)\.(\d|N)-", r"\1,\2-", c_br)
    c_br = re.sub(r"(\b\d),\s*-\s*\(", r"\1'-(", c_br)
    c_br = re.sub(r"\[(\d,\d)\s*l\s*\]", r"[\1']", c_br)
    c_br = re.sub(r"\b(\d,\d,)(\d)\2-", r"\1\2-", c_br)
    # 'l-1' stroke sequence fused into 'M' before 'H-' ("methyMH-" -> "methyl-1H-")
    c_br = re.sub(r"([a-z])MH-", r"\1l-1H-", c_br)
    # Dropped 'l' at end of '-thyl' / '-phenyl' / '-azol' suffix before punctuation or locant
    c_br = re.sub(
        r"((?:m|e|naph)th|phen)y(?=[-,\]\)\}])",
        r"\1yl",
        c_br,
        flags=re.IGNORECASE,
    )
    c_br = re.sub(
        r"([a-z]{3,}z)o(?=-\d+-yl\b)",
        r"\1ol",
        c_br,
        flags=re.IGNORECASE,
    )
    # General IUPAC vowel elision: drop terminal 'e' of '-ine' parent ring before '-N-amine/-ol/-one'
    c_br = re.sub(
        r"([a-z]{4,}in)e-(\d+-(?:amine|ol|one)\b)",
        r"\1-\2",
        c_br,
        flags=re.IGNORECASE,
    )
    c_br = re.sub(
        r"\b(mono|di|tri)\s+(hydro[a-z]+|hydrate|[a-z]+acetate|[a-z]+sulfonate)\b",
        r"\1\2",
        c_br,
        flags=re.IGNORECASE,
    )
    c_br = re.sub(
        r"(pyr|ind)?azo-\s*(?=\[\d)", r"\1azolo", c_br, flags=re.IGNORECASE
    )
    c_br = re.sub(r"\b([1234])\s+H\b", r"\1H", c_br)
    c_br = re.sub(r"\((?:\d+-)?(?:exo|endo)\)-", "", c_br, flags=re.IGNORECASE)

    # 2. Line-wrap spaces around hyphens, locant commas, and closing brackets/parens/braces
    c_ws = re.sub(r"\s*-\s*", "-", c_br)
    c_ws = re.sub(r"(?<=\d)\s*,\s*(?=\d)", ",", c_ws)
    c_ws = re.sub(r"([\)\]\}])\s+(?=[a-z\(])", r"\1", c_ws)
    c_ws = re.sub(r"([\[\(\{])\s+", r"\1", c_ws)
    c_ws = re.sub(r"\s+([\)\]\}])", r"\1", c_ws)
    # 3. Digit + capital 'O' in numeric locants (e.g., "1O-" -> "10-", "2O," -> "20,")
    c_loc = re.sub(r"(?<=\d)O(?=[-,\]\)])", "0", c_ws)
    # 4. Lowercase 'l' or uppercase 'I' misread for digit '1' in locants (e.g., "l-methyl", "l0-oxa", "-i-yl")
    c_loc = re.sub(r"(?:^|(?<=[\(\[,-]))[lI](?=\d*[-,\]\)])", "1", c_loc)
    c_loc = re.sub(r"-i-yl\b", "-1-yl", c_loc)
    # Missing '1-' locant after '[(1R)-' or '[(1S)-' before '('
    c_loc = re.sub(r"(\[\(1[RSrs]\))-\((?=\d)", r"\1-1-(", c_loc)
    # Duplicated outer locant around stereodescriptor ("1-(R)-1-<ring>-2-yl..." -> "1-[(2R)-<ring>-2-yl...]")
    c_loc = re.sub(
        r"\b(\d+)-\(([RSrs])\)-\1-([a-z]+-(\d+)-yl[a-z]*)",
        r"\1-[(\4\2)-\3]",
        c_loc,
    )
    # 5. Unbalanced parenthesis/bracket/brace from OCR or truncated heading
    if c_loc.count("(") == c_loc.count(")") + 1:
      c_loc = re.sub(r"(\([^\(\)]+?-yl)-([a-z]{4,}\b)", r"\1)\2", c_loc)

    def _strip_soft_hyphen(m: re.Match[str]) -> str:
      left, right = m.group(1), m.group(2)
      if left.lower() in _VALID_HYPHEN_PREFIXES:
        return m.group(0)
      return f"{left}{right}"

    c_dehyphen = re.sub(
        r"\b([a-z]{2,})-([a-z]{2,})\b", _strip_soft_hyphen, c_loc
    )
    extra_cands: List[str] = []
    if c_dehyphen != c_loc:
      extra_cands.append(c_dehyphen)
    for base_c in ([c_loc, c_dehyphen] if c_dehyphen != c_loc else [c_loc]):
      if (
          base_c.count("{") == base_c.count("}") + 1
          and base_c.count(")") == base_c.count("(") + 1
      ):
        for m_p in re.finditer(r"\)(?=[a-z\-])", base_c):
          extra_cands.append(base_c[: m_p.start()] + "}" + base_c[m_p.end() :])
      if (
          base_c.count("[") == base_c.count("]") + 1
          and base_c.count(")") == base_c.count("(") + 1
      ):
        for m_p in re.finditer(r"\)(?=[a-z\-])", base_c):
          extra_cands.append(base_c[: m_p.start()] + "]" + base_c[m_p.end() :])
      if (
          base_c.count(")") == base_c.count("(") + 1
          and base_c.count("[") == base_c.count("]")
          and base_c.count("{") == base_c.count("}")
      ):
        for m_p in re.finditer(r"yl\)(?=[a-z])", base_c):
          extra_cands.append(base_c[: m_p.start() + 2] + base_c[m_p.end() :])
      if base_c.count("(") == base_c.count(")") + 1:
        c_extra_p = base_c.replace("{[(", "{[", 1)
        if c_extra_p != base_c:
          extra_cands.append(c_extra_p)
        extra_cands.append(base_c + ")")
      if base_c.count("[") == base_c.count("]") + 1:
        extra_cands.append(base_c + "]")
    for c in [c_loc] + extra_cands + [c_ws, c_br]:
      if c and c != raw_name and c not in candidates:
        candidates.append(c)
    return candidates

  def _resolve_opsin(item: Tuple[str, str]) -> Tuple[str, Optional[str]]:
    ex_id, name = item
    norm_k = norm_iupac(name)
    if norm_k in _OPSIN_CACHE:
      return ex_id, _OPSIN_CACHE[norm_k]
    # Try verbatim IUPAC name first
    opsin_url = (
        f"https://opsin.ch.cam.ac.uk/opsin/{urllib.parse.quote(name)}.json"
    )
    data = _http_get_json(opsin_url, timeout=8)
    if data and data.get("smiles"):
      _OPSIN_CACHE[norm_k] = data["smiles"]
      return ex_id, data["smiles"]
    # Only if verbatim parsing fails, check for potential OCR locant/hyphen typos
    for cand_name in _repair_iupac_ocr_typos(name):
      cand_url = f"https://opsin.ch.cam.ac.uk/opsin/{urllib.parse.quote(cand_name)}.json"
      cand_data = _http_get_json(cand_url, timeout=8)
      if cand_data and cand_data.get("smiles"):
        _warn(
            f"Resolved IUPAC for Example {ex_id} after repairing OCR typo:"
            f" '{name}' -> '{cand_name}'"
        )
        _OPSIN_CACHE[norm_k] = cand_data["smiles"]
        return ex_id, cand_data["smiles"]
    _OPSIN_CACHE[norm_k] = None
    return ex_id, None

  resolved_smiles: Dict[str, str] = {}
  if unmatched_iupacs:
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
      for ex_id, smi in executor.map(_resolve_opsin, unmatched_iupacs.items()):
        if smi:
          resolved_smiles[ex_id] = smi
          nk_res = norm_iupac(example_to_iupac.get(ex_id, ""))
          if nk_res and nk_res not in norm_iupac_to_smi:
            norm_iupac_to_smi[nk_res] = smi
    opsin_failures = len(unmatched_iupacs) - len(resolved_smiles)
    if opsin_failures:
      _warn(
          f"{patent_dict.get('compact', '?')}: OPSIN could not resolve"
          f" {opsin_failures} of {len(unmatched_iupacs)} example names to a"
          " structure. Those examples are retained with assay data but no"
          " structure."
      )

  for ex_id_k, iup_k in example_to_iupac.items():
    nk_k = norm_iupac(iup_k)
    if ex_id_k not in resolved_smiles and nk_k in norm_iupac_to_smi:
      resolved_smiles[ex_id_k] = norm_iupac_to_smi[nk_k]

  _fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

  def _mh_compatible(mol: Any, ex_id: str) -> bool:
    if mol is None or ex_id not in example_to_mh_mass:
      return True
    exp_mh = example_to_mh_mass[ex_id]
    calc_mh = Descriptors.ExactMolWt(mol) + 1.0078
    return abs(calc_mh - exp_mh) <= 1.6

  # Detect copy-pasted duplicate headings across distinct non-peak, non-grouped Examples
  # where one sibling Example's LCMS [M+H]+ agrees with the heading's exact mass while
  # another sibling's LCMS [M+H]+ conflicts with the copy-pasted heading and can be
  # repaired from an explicit alkylating/coupling reagent in the experimental prose.
  non_peak_by_norm: Dict[str, List[str]] = {}
  for ex_id_np, iup_np in list(example_to_iupac.items()):
    if ex_id_np in example_to_peak or ex_id_np in example_in_group:
      continue
    nk_np = norm_iupac(iup_np)
    if nk_np:
      non_peak_by_norm.setdefault(nk_np, []).append(ex_id_np)
  for nk_dup, dup_ex_list in non_peak_by_norm.items():
    if len(dup_ex_list) < 2:
      continue
    smi_dup = norm_iupac_to_smi.get(nk_dup) or resolved_smiles.get(
        dup_ex_list[0]
    )
    mol_dup = Chem.MolFromSmiles(smi_dup) if smi_dup else None
    if mol_dup is None:
      continue
    has_matching_sib = any(
        ex_d in example_to_mh_mass and _mh_compatible(mol_dup, ex_d)
        for ex_d in dup_ex_list
    )
    if not has_matching_sib:
      continue
    for ex_d in dup_ex_list:
      if ex_d in example_to_mh_mass and not _mh_compatible(mol_dup, ex_d):
        orig_iup_d = example_to_iupac.get(ex_d, "")
        prose_d = example_to_prose.get(ex_d, "")
        repaired_ok = False
        reagent_stems = re.findall(
            r"\b(?:bromo|chloro|iodo)([a-z]{4,})ane\b",
            prose_d,
            re.IGNORECASE,
        )
        for r_stem in dict.fromkeys(s.lower() for s in reagent_stems):
          new_yl = f"{r_stem}yl"
          for sub_m in re.finditer(
              r"\b([a-z]{4,}(?:-\d+)?-yl)\b", orig_iup_d, re.IGNORECASE
          ):
            old_yl = sub_m.group(1)
            if old_yl.lower() == new_yl:
              continue
            cand_iup_d = (
                orig_iup_d[: sub_m.start(1)]
                + new_yl
                + orig_iup_d[sub_m.end(1) :]
            )
            _, cand_smi_d = _resolve_opsin((ex_d, cand_iup_d))
            cand_mol_d = Chem.MolFromSmiles(cand_smi_d) if cand_smi_d else None
            if cand_mol_d is not None and _mh_compatible(cand_mol_d, ex_d):
              _warn(
                  f"Repaired copy-pasted IUPAC heading on Example {ex_d} using"
                  f" prose reagent & LCMS ({example_to_mh_mass[ex_d]}):"
                  f" '{orig_iup_d}' -> '{cand_iup_d}'"
              )
              example_to_iupac[ex_d] = cand_iup_d
              resolved_smiles[ex_d] = cand_smi_d
              norm_iupac_to_smi[norm_iupac(cand_iup_d)] = cand_smi_d
              repaired_ok = True
              break
          if repaired_ok:
            break

  # Build the Tier A/B (IUPAC/OPSIN-verified) applicability domain anchor set
  # so Tier C [M+H]+ mass matching cannot assign off-scaffold prior-art,
  # synthetic intermediates, or co-therapy comparator agents.
  tier_ab_scaffolds = set()
  tier_ab_fps = []
  tier_ab_has_aryl_br_i = False
  for ex_id_ab, iup_ab in example_to_iupac.items():
    nk_ab = norm_iupac(iup_ab) if iup_ab else ""
    mols_ab = []
    if nk_ab and nk_ab in iupac_to_recs:
      mols_ab.extend(
          r.get("rdkit_mol")
          for r in iupac_to_recs[nk_ab]
          if r.get("rdkit_mol") is not None
      )
    if ex_id_ab in resolved_smiles:
      m_ab = Chem.MolFromSmiles(resolved_smiles[ex_id_ab])
      if m_ab is not None:
        mols_ab.append(m_ab)
    for m_ab in mols_ab:
      try:
        if _ARYL_BR_I_SMARTS is not None and m_ab.HasSubstructMatch(
            _ARYL_BR_I_SMARTS
        ):
          tier_ab_has_aryl_br_i = True
        scaf_s = MurckoScaffold.MurckoScaffoldSmiles(mol=m_ab)
        if scaf_s:
          tier_ab_scaffolds.add(scaf_s)
        tier_ab_fps.append(_fp_gen.GetFingerprint(m_ab))
      except Exception:
        pass

  def _passes_tier_ab_domain(mol: Any, min_sim: float = 0.40) -> bool:
    if mol is None or not tier_ab_fps:
      return True
    try:
      if _BORON_SI_SN_SMARTS is not None and mol.HasSubstructMatch(
          _BORON_SI_SN_SMARTS
      ):
        return False
      if (
          not tier_ab_has_aryl_br_i
          and _ARYL_BR_I_SMARTS is not None
          and mol.HasSubstructMatch(_ARYL_BR_I_SMARTS)
      ):
        return False
      scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
      if scaf and scaf in tier_ab_scaffolds:
        return True
      fp = _fp_gen.GetFingerprint(mol)
      sims = DataStructs.BulkTanimotoSimilarity(fp, tier_ab_fps)
      return bool(sims and max(sims) >= min_sim)
    except Exception:
      return True

  def _desalt_organic_mol(mol: Any) -> Any:
    if mol is None:
      return None
    try:
      frags = Chem.GetMolFrags(mol, asMols=True)
      if len(frags) <= 1:
        return mol
      org = [
          f for f in frags if any(a.GetAtomicNum() == 6 for a in f.GetAtoms())
      ]
      return max(org or frags, key=lambda f: f.GetNumHeavyAtoms())
    except Exception:
      return mol

  def _inject_assay_columns(
      target_dict: Dict[str, Any], ex_id: str, assays: List[str]
  ) -> None:
    pat_prefix = re.sub(
        r"[^A-Za-z0-9]+",
        "",
        str(patent_dict.get("compact") or used_id or "PAT"),
    )
    for idx_a, a_val in enumerate(assays, 1):
      if not a_val or not str(a_val).strip():
        continue
      a_clean = str(a_val).strip()
      if a_clean == "13" and "B" in bin_legend:
        a_clean = "B"
      target_dict[f"patent_table_col_{idx_a}"] = a_clean
      target_dict[f"{pat_prefix}_col_{idx_a}"] = a_clean
      nm_v, p_v = _parse_pactivity(a_clean)
      if nm_v is not None:
        target_dict[f"patent_table_col_{idx_a}_nM"] = round(nm_v, 4)
        target_dict[f"patent_table_col_{idx_a}_pActivity"] = p_v
        target_dict[f"{pat_prefix}_col_{idx_a}_nM"] = round(nm_v, 4)
      elif a_clean in bin_legend or _is_discrete_potency_token(a_clean):
        b_rng = bin_legend.get(a_clean)
        if not target_dict.get("activity_bin"):
          target_dict["activity_bin"] = a_clean
          if b_rng:
            target_dict["activity_range"] = b_rng
        if b_rng:
          target_dict[f"patent_table_col_{idx_a}_range"] = b_rng
    for hdr_lbl, a_val in example_to_assay_cols.get(ex_id, []):
      a_clean = str(a_val).strip()
      if a_clean == "13" and "B" in bin_legend:
        a_clean = "B"
      norm_hdr_lbl = hdr_lbl.replace("µ", "u").replace("μ", "u")
      clean_hdr = re.sub(r"[^a-zA-Z0-9]+", "_", norm_hdr_lbl).strip("_").lower()
      if clean_hdr and not clean_hdr.startswith("col_"):
        col_key = f"patent_assay_{clean_hdr}"
        target_dict[col_key] = a_clean
        nm_v, p_v = _parse_pactivity(a_clean)
        if nm_v is not None:
          target_dict[f"{col_key}_nM"] = round(nm_v, 4)
          target_dict[f"{col_key}_pActivity"] = p_v
        elif a_clean in bin_legend or _is_discrete_potency_token(a_clean):
          b_rng = bin_legend.get(a_clean)
          if not target_dict.get("activity_bin"):
            target_dict["activity_bin"] = a_clean
            if b_rng:
              target_dict["activity_range"] = b_rng
          if b_rng:
            target_dict[f"{col_key}_range"] = b_rng
        if not target_dict.get("target_name"):
          inferred_tgt = _extract_target_from_assay_header(hdr_lbl)
          if inferred_tgt:
            target_dict["target_name"] = inferred_tgt

  assigned_record_ids = set()
  enriched_count = 0
  added_count = 0
  structureless_count = 0

  all_ex_ids = sorted(
      set(list(example_to_assays.keys()) + list(example_to_iupac.keys())),
      key=lambda x: (
          int(re.findall(r"\d+", x)[0]) if re.findall(r"\d+", x) else 9999,
          x,
      ),
  )

  for ex_id in all_ex_ids:
    assays = example_to_assays.get(ex_id, [])
    primary_nm, primary_p = None, None
    primary_bin: Optional[str] = None
    primary_bin_range: Optional[str] = None
    primary_bin_hdr: Optional[str] = None
    for a_str in assays:
      if not a_str or not str(a_str).strip():
        continue
      a_clean = str(a_str).strip()
      if a_clean == "13" and "B" in bin_legend:
        a_clean = "B"
      nm_val, p_val = _parse_pactivity(a_clean)
      if p_val is not None:
        primary_nm, primary_p = nm_val, p_val
        break
      if primary_bin is None and (
          a_clean in bin_legend or _is_discrete_potency_token(a_clean)
      ):
        primary_bin = a_clean
        primary_bin_range = bin_legend.get(a_clean)
    if primary_bin is not None:
      for hdr_lbl, a_val in example_to_assay_cols.get(ex_id, []):
        a_c = (
            "B"
            if (str(a_val).strip() == "13" and "B" in bin_legend)
            else str(a_val).strip()
        )
        if a_c == primary_bin and not hdr_lbl.startswith("col_"):
          primary_bin_hdr = hdr_lbl
          break

    iup = example_to_iupac.get(ex_id, "")
    norm_key = norm_iupac(iup) if iup else ""
    peak_label = (
        f"Peak {example_to_peak[ex_id]}" if ex_id in example_to_peak else ""
    )
    ex_smi = resolved_smiles.get(ex_id) or (
        norm_iupac_to_smi.get(norm_key, "") if norm_key else ""
    )
    if not ex_smi and ex_id in pubchem_ex_to_recs and pubchem_ex_to_recs[ex_id]:
      unclaimed_pex = [
          r
          for r in pubchem_ex_to_recs[ex_id]
          if id(r) not in assigned_record_ids
      ]
      if unclaimed_pex:
        ex_smi = str(unclaimed_pex[0].get("canonical_smiles") or "")
        if not iup:
          iup = str(unclaimed_pex[0].get("iupac_name") or "")

    def _can_claim_cand(cand: Dict[str, Any], cur_ex: str) -> bool:
      if id(cand) in assigned_record_ids:
        return False
      cand_pex = cand.get("_pubchem_ex_id")
      if cand_pex and cand_pex != cur_ex:
        return False
      return True

    target_rec = None
    # Tier 0: Direct PubChem Example title/synonym match ("US..., Example 343")
    if ex_id in pubchem_ex_to_recs:
      for cand in pubchem_ex_to_recs[ex_id]:
        if id(cand) not in assigned_record_ids:
          target_rec = cand
          break

    # Tier A: Direct normalized IUPAC match against PubChem records
    if target_rec is None and norm_key and norm_key in iupac_to_recs:
      for cand in iupac_to_recs[norm_key]:
        if _can_claim_cand(cand, ex_id) and _mh_compatible(
            cand.get("rdkit_mol"), ex_id
        ):
          target_rec = cand
          break
      if target_rec is None:
        for cand in iupac_to_recs[norm_key]:
          if _can_claim_cand(cand, ex_id):
            target_rec = cand
            break

    # Tier B: Match via OPSIN InChIKey (exact or 14-char connectivity layer)
    if target_rec is None and ex_smi:
      m = _desalt_organic_mol(Chem.MolFromSmiles(ex_smi))
      if m:
        ik = Chem.MolToInchiKey(m)
        if ik in ik_to_recs:
          for cand in ik_to_recs[ik]:
            if _can_claim_cand(cand, ex_id):
              target_rec = cand
              break
        if target_rec is None and ik:
          conn_block = ik.split("-")[0]
          if conn_block in conn_to_recs:
            for cand in conn_to_recs[conn_block]:
              if _can_claim_cand(cand, ex_id):
                target_rec = cand
                break

    # Tier C: Match unlinked PubChem record by experimental LCMS [M+H]+ mass ONLY when no IUPAC/OPSIN structure
    # is resolved for this Example (or its chiral separation group), gated by Tier A/B Murcko scaffold + Tanimoto (0.40).
    if target_rec is None and not ex_smi and ex_id in example_to_mh_mass:
      exp_mh = example_to_mh_mass[ex_id]
      mh_cands = [
          r
          for r in records
          if _can_claim_cand(r, ex_id)
          and not r.get("patent_example_id")
          and r.get("rdkit_mol") is not None
          and abs((Descriptors.ExactMolWt(r["rdkit_mol"]) + 1.0078) - exp_mh)
          <= 0.65
          and _passes_tier_ab_domain(r["rdkit_mol"], min_sim=0.40)
      ]
      if mh_cands:
        # Rank by closest exact mass match, breaking ties by stereochemical completeness
        mh_cands.sort(
            key=lambda r: (
                abs((Descriptors.ExactMolWt(r["rdkit_mol"]) + 1.0078) - exp_mh),
                -_stereo_richness(r.get("canonical_smiles")),
            )
        )
        target_rec = mh_cands[0]

    bin_display = (
        f"{primary_bin} ({primary_bin_range})"
        if (primary_bin and primary_bin_range)
        else primary_bin
    )
    bin_type = (
        f"{primary_bin_hdr} (Patent Table Bin)"
        if primary_bin_hdr
        else "Potency Bin (Patent Table)"
    )

    if target_rec is not None:
      assigned_record_ids.add(id(target_rec))
      if (
          norm_key
          and norm_key not in norm_iupac_to_smi
          and target_rec.get("canonical_smiles")
      ):
        norm_iupac_to_smi[norm_key] = str(target_rec["canonical_smiles"])
      pat_code = str(
          target_rec.get("patent_id") or patent_dict["compact"]
      ).strip()
      target_rec["patent_example_id"] = (
          f"{pat_code}:Example_{ex_id}" if pat_code else f"Example_{ex_id}"
      )
      if peak_label:
        target_rec["chiral_peak"] = peak_label
      target_rec["all_sources"] = ",".join(
          sorted({str(target_rec.get("source", "")), "Patent_FullText_Table"})
      )
      if primary_p is not None and (
          target_rec.get("pchembl_value") is None
          or primary_p > float(target_rec.get("pchembl_value", 0))
      ):
        target_rec["activity_type"] = "Ki/IC50 (Patent Table)"
        target_rec["activity_value"] = primary_nm
        target_rec["activity_units"] = "nM"
        target_rec["pchembl_value"] = primary_p
      elif primary_bin is not None:
        target_rec["activity_bin"] = primary_bin
        if primary_bin_range:
          target_rec["activity_range"] = primary_bin_range
        if (
            target_rec.get("pchembl_value") is None
            and target_rec.get("activity_value") is None
        ):
          target_rec["activity_type"] = bin_type
          target_rec["activity_value"] = bin_display
          target_rec["activity_units"] = "bin"
      _inject_assay_columns(target_rec, ex_id, assays)
      enriched_count += 1
    elif ex_smi or assays:
      smi = ex_smi
      mol = _desalt_organic_mol(Chem.MolFromSmiles(smi)) if smi else None
      can_smi = Chem.MolToSmiles(mol, isomericSmiles=True) if mol else None
      ik = Chem.MolToInchiKey(mol) if mol else None
      if mol is None:
        structureless_count += 1
      max_ring = (
          max([len(r) for r in mol.GetRingInfo().AtomRings()])
          if (mol and mol.GetRingInfo().AtomRings())
          else (0 if mol is not None else None)
      )
      pat_code = str(patent_dict["compact"]).strip()
      qual_ex_id = (
          f"{pat_code}:Example_{ex_id}" if pat_code else f"Example_{ex_id}"
      )
      new_rec = {
          "compound_id": qual_ex_id,
          "patent_example_id": qual_ex_id,
          "chiral_peak": peak_label,
          "smiles": can_smi,
          "canonical_smiles": can_smi,
          "iupac_name": iup,
          "inchikey": ik,
          "rdkit_mol": mol,
          "max_ring_size": max_ring,
          "is_macrocycle": (max_ring >= 12) if max_ring is not None else None,
          "source": "Patent_FullText_Table",
          "patent_id": patent_dict["compact"],
          "activity_type": (
              "Ki/IC50 (Patent Table)"
              if primary_p is not None
              else (bin_type if primary_bin is not None else None)
          ),
          "activity_value": (
              primary_nm if primary_nm is not None else bin_display
          ),
          "activity_units": (
              "nM"
              if primary_nm is not None
              else ("bin" if primary_bin is not None else None)
          ),
          "pchembl_value": primary_p,
          "activity_bin": primary_bin,
          "activity_range": primary_bin_range,
      }
      _inject_assay_columns(new_rec, ex_id, assays)
      records.append(new_rec)
      assigned_record_ids.add(id(new_rec))
      added_count += 1
      if primary_p is not None or primary_bin is not None:
        enriched_count += 1

  _apply_patent_dates_to_records(records, patent_dates)
  print(
      f"[patent_fetcher] Full-text HTML enrichment complete: {enriched_count}"
      f" compounds with Table/inline bioactivity ({added_count} newly added"
      f" structures, of which {structureless_count} have assay data but no"
      " resolved structure)",
      file=sys.stderr,
      flush=True,
  )


_BOC_SMARTS = Chem.MolFromSmarts("CC(C)(C)OC(=O)N")
_CBZ_SMARTS = Chem.MolFromSmarts("c1ccccc1COC(=O)N")
_BORON_SI_SN_SMARTS = Chem.MolFromSmarts("[#5,#14,#50]")
_ARYL_BR_I_SMARTS = Chem.MolFromSmarts("c[Br,I]")


def filter_records_by_exemplified_domain(
    records: List[Dict[str, Any]],
    min_tanimoto: float = 0.40,
    keep_all_refs: bool = False,
    query_smarts: Optional[str] = None,
    query_scaffold: Optional[str] = None,
    prune_unexemplified_unassayed: bool = False,
) -> List[Dict[str, Any]]:
  """Filters un-exemplified associative records (PubChem text-mined & OCSR concepts)

  against the patent's exemplified compounds using Murcko scaffold overlap or
  Morgan ECFP4 Tanimoto similarity >= min_tanimoto (default 0.40, matching
  predict_potency.py), strips synthetic intermediates (boronic esters,
  silyl/stannyl,
  Boc/Cbz, aryl-Br/I), and optionally prunes unassayed un-exemplified PubChem
  bloat
  when full-text HTML Examples exist.
  """
  fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

  smarts_pat = Chem.MolFromSmarts(query_smarts) if query_smarts else None
  scaf_pat = Chem.MolFromSmiles(query_scaffold) if query_scaffold else None
  if query_scaffold and scaf_pat is None:
    scaf_pat = Chem.MolFromSmarts(query_scaffold)

  exemplified_scaffolds = set()
  exemplified_fps = []
  exemplified_conn_blocks: Set[str] = set()
  exemplified_has_boc = False
  exemplified_has_aryl_halide = False
  has_html_examples = False
  for r in records:
    if r.get("patent_example_id"):
      has_html_examples = True
    if r.get("patent_example_id") or r.get("source") == "ChEMBL":
      m = r.get("rdkit_mol")
      if m is None and r.get("canonical_smiles"):
        m = Chem.MolFromSmiles(str(r["canonical_smiles"]))
      if m is not None:
        try:
          ik = r.get("inchikey") or Chem.MolToInchiKey(m)
          if ik:
            exemplified_conn_blocks.add(str(ik).split("-")[0])
          if _BOC_SMARTS is not None and m.HasSubstructMatch(_BOC_SMARTS):
            exemplified_has_boc = True
          if _ARYL_BR_I_SMARTS is not None and m.HasSubstructMatch(
              _ARYL_BR_I_SMARTS
          ):
            exemplified_has_aryl_halide = True
          scaf_s = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
          if scaf_s:
            exemplified_scaffolds.add(scaf_s)
          exemplified_fps.append(fp_gen.GetFingerprint(m))
        except Exception:
          pass

  def _is_unexemplified_intermediate(mol_obj: Any) -> bool:
    if mol_obj is None:
      return False
    if _BORON_SI_SN_SMARTS is not None and mol_obj.HasSubstructMatch(
        _BORON_SI_SN_SMARTS
    ):
      return True
    if (
        not exemplified_has_boc
        and _BOC_SMARTS is not None
        and mol_obj.HasSubstructMatch(_BOC_SMARTS)
    ):
      return True
    if (
        not exemplified_has_boc
        and _CBZ_SMARTS is not None
        and mol_obj.HasSubstructMatch(_CBZ_SMARTS)
    ):
      return True
    if (
        not exemplified_has_aryl_halide
        and _ARYL_BR_I_SMARTS is not None
        and mol_obj.HasSubstructMatch(_ARYL_BR_I_SMARTS)
    ):
      return True
    return False

  kept: List[Dict[str, Any]] = []
  pruned_count = 0
  for r in records:
    m = r.get("rdkit_mol")
    if m is None and r.get("canonical_smiles"):
      m = Chem.MolFromSmiles(str(r["canonical_smiles"]))
      r["rdkit_mol"] = m

    # Explicit user SMARTS / scaffold filter applies to all structural rows
    if m is not None:
      if smarts_pat is not None and not m.HasSubstructMatch(smarts_pat):
        pruned_count += 1
        continue
      if scaf_pat is not None and not m.HasSubstructMatch(scaf_pat):
        pruned_count += 1
        continue

    # Always keep exemplified rows (or assay-only rows with an Example ID) and ChEMBL rows
    if (
        keep_all_refs
        or r.get("patent_example_id")
        or r.get("source") == "ChEMBL"
    ):
      kept.append(r)
      continue

    # Strip obvious synthetic intermediates from un-exemplified associative/OCSR records
    if _is_unexemplified_intermediate(m):
      pruned_count += 1
      continue

    # When full-text HTML examples were extracted, prune un-exemplified PubChem
    # patent-association rows that either duplicate an already-exemplified structure
    # or carry no verified target bioactivity.
    if (
        prune_unexemplified_unassayed
        and has_html_examples
        and r.get("source") == "PubChem"
    ):
      r_conn_block = ""
      if m is not None:
        try:
          r_ik = r.get("inchikey") or Chem.MolToInchiKey(m)
          if r_ik:
            r_conn_block = str(r_ik).split("-")[0]
        except Exception:
          pass
      if (
          (r_conn_block and r_conn_block in exemplified_conn_blocks)
          or (
              r.get("pchembl_value") is None
              and r.get("activity_value") is None
              and r.get("activity_bin") is None
          )
          or str(r.get("target_name") or "").strip().lower()
          in (
              "",
              "nan",
              "none",
              "unchecked",
              "unspecified",
              "unknown",
              "no relevant target",
              "admet",
          )
      ):
        pruned_count += 1
        continue

    # If no exemplified compounds exist in the patent to define an applicability domain, keep
    if not exemplified_fps or m is None:
      kept.append(r)
      continue

    try:
      scaf_s = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
      if scaf_s and scaf_s in exemplified_scaffolds:
        r["max_tanimoto_to_exemplified"] = 1.0
        kept.append(r)
        continue
      fp = fp_gen.GetFingerprint(m)
      sims = DataStructs.BulkTanimotoSimilarity(fp, exemplified_fps)
      max_sim = max(sims) if sims else 0.0
      r["max_tanimoto_to_exemplified"] = round(float(max_sim), 3)
      if max_sim >= min_tanimoto:
        kept.append(r)
      else:
        pruned_count += 1
    except Exception:
      kept.append(r)

  if pruned_count > 0:
    print(
        "[patent_fetcher] Applicability domain filter (Murcko scaffold or"
        f" Tanimoto >= {min_tanimoto:.2f}): retained {len(kept)} compounds,"
        f" pruned {pruned_count} off-scaffold associative/concept records.",
        file=sys.stderr,
        flush=True,
    )
  return kept


def fetch_patent_compounds(
    patent_id: str,
    limit: Optional[int] = None,
    preferred_source: str = "all",
    min_tanimoto: float = 0.40,
    keep_all_pubchem_refs: bool = False,
    query_smarts: Optional[str] = None,
    query_scaffold: Optional[str] = None,
    allowed_elements: Optional[Set[int]] = None,
    skip_external_bioassays: bool = False,
) -> List[Dict[str, Any]]:
  """Master retrieval function combining ChEMBL, PubChem, and full-text Google Patents HTML tables."""
  p_dict = normalize_patent_id(patent_id)
  all_records: List[Dict[str, Any]] = []

  if preferred_source in ("chembl", "all"):
    chembl_mols = fetch_patent_chembl(p_dict, limit=limit)
    all_records.extend(chembl_mols)

  if preferred_source in ("pubchem", "all") or not all_records:
    pubchem_mols = fetch_patent_pubchem(
        p_dict, limit=limit, allowed_elements=allowed_elements
    )
    seen_smiles = {r["smiles"] for r in all_records}
    for pm in pubchem_mols:
      if pm["smiles"] not in seen_smiles:
        all_records.append(pm)
        seen_smiles.add(pm["smiles"])

  validated_records = []
  for rec in all_records:
    smi = rec.get("smiles", "")
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
      rec["rdkit_mol"] = mol
      rec["canonical_smiles"] = Chem.MolToSmiles(mol)
      ring_info = mol.GetRingInfo()
      ring_sizes = [len(r) for r in ring_info.AtomRings()]
      max_ring = max(ring_sizes) if ring_sizes else 0
      rec["max_ring_size"] = max_ring
      rec["is_macrocycle"] = max_ring >= 12
      validated_records.append(rec)

  enrich_records_with_fulltext_patent_html(p_dict, validated_records)
  validated_records = filter_records_by_exemplified_domain(
      validated_records,
      min_tanimoto=min_tanimoto,
      keep_all_refs=keep_all_pubchem_refs,
      query_smarts=query_smarts,
      query_scaffold=query_scaffold,
      prune_unexemplified_unassayed=False,
  )
  if limit is not None and len(validated_records) > limit:
    exemplified_recs = [
        r
        for r in validated_records
        if r.get("patent_example_id") or r.get("source") == "ChEMBL"
    ]
    assoc_recs = [
        r
        for r in validated_records
        if not (r.get("patent_example_id") or r.get("source") == "ChEMBL")
    ]
    assoc_budget = max(0, limit - len(exemplified_recs))
    validated_records = exemplified_recs + assoc_recs[:assoc_budget]
  html_ex_recs = [r for r in validated_records if r.get("patent_example_id")]
  has_html_assays = bool(html_ex_recs) and any(
      r.get("pchembl_value") is not None
      or r.get("activity_value") is not None
      or r.get("activity_bin") is not None
      for r in html_ex_recs
  )
  skip_ext = has_html_assays or skip_external_bioassays
  enrich_records_with_chembl_by_inchikey(
      validated_records, fetch_activities=not skip_ext
  )
  if not skip_ext:
    enrich_records_with_pubchem_bioassays(validated_records)
  validated_records = filter_records_by_exemplified_domain(
      validated_records,
      min_tanimoto=min_tanimoto,
      keep_all_refs=keep_all_pubchem_refs,
      query_smarts=query_smarts,
      query_scaffold=query_scaffold,
      prune_unexemplified_unassayed=True,
  )
  merged_dates = dict(p_dict.get("_patent_dates") or {})
  if not merged_dates:
    merged_dates = _extract_patent_dates_from_soup(None, p_dict)
  if merged_dates.get("document_year") is None:
    sibling_years = []
    for rec in validated_records:
      ry = rec.get("document_year")
      if ry is not None and not isinstance(ry, bool):
        try:
          ry_int = int(float(ry))
          if 1800 <= ry_int <= 2100:
            sibling_years.append(ry_int)
        except (TypeError, ValueError):
          pass
    if sibling_years:
      merged_dates["document_year"] = min(sibling_years)
  _apply_patent_dates_to_records(validated_records, merged_dates)

  for rec in validated_records:
    has_act = (
        (rec.get("pchembl_value") is not None)
        or (rec.get("activity_value") is not None)
        or (rec.get("activity_bin") is not None)
    )
    rec["dataset_role"] = "labeled_patent" if has_act else "prospective_patent"
    rec.setdefault("all_sources", str(rec.get("source", "")))
  return validated_records


# Google Patents' OCSR writes unresolved R-groups and attachment points as
# dummy atoms, but it also emits [Y], [U], [V] and [W], which RDKit parses
# happily as yttrium, uranium, vanadium and tungsten. A generic structure will
# therefore sail through MolFromSmiles unless it is caught as a string first.
_PSEUDO_ATOM_TOKENS = ("[Y]", "[U]", "[V]", "[W]")


def _is_generic_structure(smiles: str) -> bool:
  """Reports whether a SMILES carries Markush placeholders rather than atoms."""
  if "*" in smiles:
    return True
  return any(tok in smiles for tok in _PSEUDO_ATOM_TOKENS)


def _stereo_richness(smiles: Optional[str]) -> int:
  """Counts declared stereo markers, used to pick between two reads of a compound."""
  if not smiles:
    return 0
  return smiles.count("@") + smiles.count("/") + smiles.count("\\")


def _concept_span(li: Any, prop: str) -> Optional[str]:
  """Returns the text of the first child span carrying the given itemprop."""
  el = li.find("span", {"itemprop": prop})
  if el is None:
    return None
  return el.get_text(strip=True) or None


def fetch_patent_concepts(
    patent_id: str, min_heavy_atoms: int = 6
) -> List[Dict[str, Any]]:
  """Harvests the structures Google Patents' own OCSR read out of the figures.

  Patent structures are drawn, not named, so neither the example tables nor
  OPSIN can recover most of them. Google already runs OCSR over the figures and
  publishes the result as the page's concept index. Those records are almost
  all dot-separated mixtures, because each figure is a multi-step reaction
  scheme rather than a single compound, so every fragment is emitted as its own
  candidate row and left for a downstream substructure filter to triage.

  Args:
    patent_id: A raw patent identifier in any form normalize_patent_id accepts.
    min_heavy_atoms: Fragments smaller than this are dropped as reagents,
      counterions and solvent.

  Returns:
    A list of records in the same shape the rest of this module emits, tagged
    source="Patent_Concept_OCSR", plus any generic structures tagged
    source="Patent_Markush". Concept records are not indexed to an example
    number by Google, so patent_example_id is always None.
  """
  p_dict = normalize_patent_id(patent_id)
  used_id = p_dict["compact"]
  url = f"https://patents.google.com/patent/{used_id}/en"
  html = _http_get_text(url, timeout=20)
  if not html:
    _warn(f"Concept harvest failed for {patent_id} ({url})")
    return []

  soup = BeautifulSoup(html, "html.parser")
  concept_dates = _extract_patent_dates_from_soup(soup, p_dict)
  matches = soup.find_all("li", {"itemprop": "match"})
  if not matches:
    _warn(
        f"No concept records found on the Google Patents page for {used_id}."
        " The page markup may have changed."
    )
    return []

  records: List[Dict[str, Any]] = []
  # Fragments are deduplicated on the skeleton (connectivity) block of the
  # InChIKey so that the same compound read out of two different schemes, or
  # read once with and once without stereochemistry, does not appear twice.
  by_skeleton: Dict[str, Dict[str, Any]] = {}
  n_chemical = 0
  n_markush = 0
  n_unparsed = 0
  n_too_small = 0
  n_fragments_seen = 0

  for li in matches:
    if _concept_span(li, "domain") != "Chemical compound":
      continue
    smiles = _concept_span(li, "smiles")
    if not smiles:
      continue
    n_chemical += 1

    count_txt = _concept_span(li, "count")
    try:
      concept_count = int(count_txt) if count_txt else None
    except ValueError:
      concept_count = None
    sections = ";".join(
        sorted({
            s.get_text(strip=True)
            for s in li.find_all("span", {"itemprop": "sections"})
            if s.get_text(strip=True)
        })
    )
    concept_ik = _concept_span(li, "inchi_key")

    if _is_generic_structure(smiles):
      # Kept, not dropped: a general formula is the one place the patent states
      # its claimed scope. It is stored as a raw string because numbered dummy
      # atoms make it unparseable as a molecule, and it is deliberately NOT
      # converted into a substructure filter -- the OCSR read of a general
      # formula image is too garbled to match its own exemplified compounds.
      n_markush += 1
      records.append({
          "compound_id": f"Markush_{used_id}_{n_markush}",
          "patent_example_id": None,
          "smiles": None,
          "canonical_smiles": None,
          "iupac_name": None,
          "inchikey": concept_ik,
          "rdkit_mol": None,
          "max_ring_size": 0,
          "is_macrocycle": False,
          "source": "Patent_Markush",
          "all_sources": "Patent_Markush",
          "patent_id": used_id,
          "activity_type": None,
          "activity_value": None,
          "activity_units": None,
          "pchembl_value": None,
          "markush_smiles": smiles,
          "concept_count": concept_count,
          "concept_sections": sections,
      })
      continue

    fragments = [f.strip() for f in smiles.split(".") if f.strip()]
    for frag in fragments:
      n_fragments_seen += 1
      mol = Chem.MolFromSmiles(frag)
      if mol is None:
        n_unparsed += 1
        continue
      if mol.GetNumHeavyAtoms() < min_heavy_atoms:
        n_too_small += 1
        continue
      can_smi = Chem.MolToSmiles(mol)
      ik = Chem.MolToInchiKey(mol)
      skeleton = ik.split("-")[0] if ik else can_smi
      existing = by_skeleton.get(skeleton)
      if existing is not None:
        # The same fragment appearing in several schemes is a stronger signal
        # that it is a real product, so keep the largest occurrence count.
        if concept_count is not None and (
            existing.get("concept_count") is None
            or concept_count > existing["concept_count"]
        ):
          existing["concept_count"] = concept_count
        # Deduplicating on the skeleton means a flat read and a
        # stereochemically defined read of the same compound collide. Keep the
        # richer one: OCSR frequently reads the wedge bonds in one depiction
        # and misses them in another, and first-seen would otherwise decide it.
        if _stereo_richness(can_smi) > _stereo_richness(
            existing["canonical_smiles"]
        ):
          existing["smiles"] = can_smi
          existing["canonical_smiles"] = can_smi
          existing["inchikey"] = ik
          existing["rdkit_mol"] = mol
        continue
      ring_sizes = [len(r) for r in mol.GetRingInfo().AtomRings()]
      max_ring = max(ring_sizes) if ring_sizes else 0
      rec = {
          "compound_id": f"Concept_{used_id}_{len(by_skeleton) + 1}",
          "patent_example_id": None,
          "smiles": can_smi,
          "canonical_smiles": can_smi,
          "iupac_name": None,
          "inchikey": ik,
          "rdkit_mol": mol,
          "max_ring_size": max_ring,
          "is_macrocycle": max_ring >= 12,
          "source": "Patent_Concept_OCSR",
          "all_sources": "Patent_Concept_OCSR",
          "patent_id": used_id,
          "activity_type": None,
          "activity_value": None,
          "activity_units": None,
          "pchembl_value": None,
          "concept_count": concept_count,
          "concept_sections": sections,
          "concept_parent_inchikey": concept_ik,
          "concept_parent_fragments": len(fragments),
      }
      by_skeleton[skeleton] = rec
      records.append(rec)

  _apply_patent_dates_to_records(records, concept_dates)
  n_macro = sum(1 for r in records if r.get("is_macrocycle"))
  print(
      f"[patent_fetcher] Concept harvest for {used_id}: {n_chemical} chemical"
      f" concept records -> {n_fragments_seen} fragments -> {len(by_skeleton)}"
      f" unique structures ({n_macro} macrocyclic); {n_markush} generic"
      f" structures retained separately; dropped {n_unparsed} unparseable and"
      f" {n_too_small} below {min_heavy_atoms} heavy atoms.",
      file=sys.stderr,
      flush=True,
  )
  if records:
    print(
        "[patent_fetcher] Concept records are not indexed to example numbers"
        " by Google, so they carry no patent_example_id and no activity."
        " Triage them with a substructure filter before treating them as a"
        " series.",
        file=sys.stderr,
        flush=True,
    )
  return records


_RELATED_STOP_WORDS = frozenset({
    "compounds",
    "compound",
    "derivatives",
    "derivative",
    "treatment",
    "treating",
    "diseases",
    "disease",
    "useful",
    "inhibitors",
    "inhibitor",
    "protein",
    "kinase",
    "kinases",
    "substituted",
    "pure",
    "enantiomerically",
    "optically",
    "novel",
    "heterocyclic",
    "heteroaryl",
    "heteroaromatic",
    "bicyclic",
    "agents",
    "agent",
    "method",
    "methods",
    "thereof",
    "pharmaceutical",
    "pharmaceutically",
    "acceptable",
    "compositions",
    "composition",
    "for",
    "the",
    "and",
    "with",
    "from",
    "that",
    "this",
    "which",
    "wherein",
    "having",
    "formula",
    "selected",
    "group",
    "consisting",
    "alkyl",
    "aryl",
    "hydrogen",
    "halogen",
    "methyl",
    "ethyl",
    "phenyl",
    "example",
    "examples",
    "table",
    "preparation",
    "synthesis",
    "step",
    "solution",
    "mixture",
    "added",
    "acid",
    "water",
    "ester",
    "amine",
    "amino",
    "human",
    "activity",
    "assay",
    "values",
    "data",
    "patent",
    "invention",
    "relates",
    "mediated",
    "abnormal",
    "cell",
    "cells",
    "growth",
    "cancer",
    "cancers",
    "anticancer",
    "hyperproliferative",
    "proliferative",
    "disorder",
    "disorders",
    "are",
    "provided",
    "well",
    "their",
    "use",
    "uses",
    "such",
    "preferred",
    "potent",
    "inhibit",
    "inhibiting",
    "present",
    "defined",
    "herein",
    "further",
    "comprising",
    "salts",
    "salt",
    "ameliorating",
    "mammals",
    "mammal",
    "administering",
    "effective",
    "amount",
    "least",
    "one",
    "two",
    "three",
    "also",
    "directed",
    "includes",
    "including",
})


def _expand_iupac_and_plural_tokens(tokens: Set[str]) -> Set[str]:
  """Expands tokens with general singular/plural and IUPAC elided-vowel variants."""
  expanded = set(tokens)
  for tok in list(tokens):
    if len(tok) < 4:
      continue
    if tok.endswith("s") and len(tok) >= 5:
      expanded.add(tok[:-1])
    else:
      expanded.add(f"{tok}s")
    if not tok.endswith(("e", "s", "y", "i", "o", "u")):
      expanded.add(f"{tok}e")
      expanded.add(f"{tok}es")
    for prefix in (
        "amino",
        "alkoxy",
        "hydroxy",
        "methyl",
        "fluoro",
        "chloro",
        "bromo",
        "cyano",
        "oxo",
    ):
      if tok.startswith(prefix) and len(tok) - len(prefix) >= 5:
        stem = tok[len(prefix) :]
        if stem not in _RELATED_STOP_WORDS:
          expanded.add(stem)
          if not stem.endswith("e"):
            expanded.add(f"{stem}e")
            expanded.add(f"{stem}es")
  return expanded - _RELATED_STOP_WORDS


def discover_related_patents(
    patent_id: str, max_related: int = 5
) -> List[Dict[str, Any]]:
  """Discovers linked continuation, DOCDB/PCT family, parent, and follow-on patents.

  Scrapes the Google Patents specification metadata for:
  1. Direct continuations, divisionals, child/parent applications, and DOCDB/PCT
     family members (US B1/B2 and WO A1/A2).
  2. Backward and forward citations that share applicant/assignee tokens and
     chemical/target title tokens with the seed patent's title, abstract,
     example headings, and assay headers, excluding solid-form/polymorph/salt/
     formulation process patents.

  Args:
    patent_id: A raw patent identifier in any form normalize_patent_id accepts.
    max_related: Maximum number of related patents to return per seed patent.

  Returns:
    A list of metadata dicts (`patent_id`, `seed_patent`, `relationship`,
    `title`, `score`) ordered by relevance.
  """
  p_dict = normalize_patent_id(patent_id)
  used_id = p_dict["compact"]
  url = f"https://patents.google.com/patent/{used_id}/en"
  html = _http_get_text(url, timeout=10)
  if not html:
    _warn(f"Patent family lookup failed for {patent_id} ({url})")
    return []

  soup = BeautifulSoup(html, "html.parser")
  seed_title_el = soup.find("meta", {"name": "DC.title"})
  seed_title = (seed_title_el.get("content") if seed_title_el else "") or ""
  seed_desc_el = soup.find("meta", {"name": "DC.description"})
  seed_desc = (seed_desc_el.get("content") if seed_desc_el else "") or ""
  title_tokens_raw = {
      w.lower()
      for w in re.findall(r"[A-Za-z]{3,}", f"{seed_title} {seed_desc}")
      if w.lower() not in _RELATED_STOP_WORDS and len(w) >= 3
  }
  title_tokens = _expand_iupac_and_plural_tokens(title_tokens_raw)
  seed_tokens = set(title_tokens)

  # Extract dominant chemical ring stems from Example IUPAC headings and
  # biological target tokens from assay headers in the specification.
  heading_els = soup.find_all("heading") + soup.find_all(
      "p", {"id": re.compile(r"^h-", re.I)}
  )
  headings_txt = " ".join(h.get_text(" ", strip=True) for h in heading_els)
  for m_ring_am in re.finditer(
      r"\b([a-z]{5,})-(?:\d+-)?(?:yl)?amin[eo]?\b", headings_txt, re.IGNORECASE
  ):
    ring_stem = m_ring_am.group(1).lower()
    if ring_stem not in _RELATED_STOP_WORDS:
      seed_tokens.update((
          ring_stem,
          f"{ring_stem}e",
          f"{ring_stem}es",
          f"amino{ring_stem}",
          f"amino{ring_stem}e",
          f"amino{ring_stem}es",
      ))
  yl_counts: Dict[str, int] = {}
  for m_yl in re.finditer(
      r"\b([a-z]{5,})-(?:\d+-)?yl\b", headings_txt, re.IGNORECASE
  ):
    stem = m_yl.group(1).lower()
    if stem not in _RELATED_STOP_WORDS:
      yl_counts[stem] = yl_counts.get(stem, 0) + 1
  for stem, cnt in yl_counts.items():
    if cnt >= 5:
      seed_tokens.update((stem, f"{stem}e", f"{stem}es"))
  for m_hdr in re.finditer(
      r"\b([A-Za-z][A-Za-z0-9-]{2,10})\s+(?:WT\s+)?(?:Ki|IC\s*50|EC\s*50)\b",
      html,
      re.IGNORECASE,
  ):
    t_tok = m_hdr.group(1).lower().replace("-", "")
    if t_tok not in _RELATED_STOP_WORDS and len(t_tok) >= 3:
      seed_tokens.add(t_tok)

  _corp_stop = {
      "inc",
      "corp",
      "corporation",
      "ltd",
      "limited",
      "llc",
      "company",
      "pharmaceuticals",
      "products",
      "therapeutics",
      "sciences",
  }
  assignees = set()
  for dd in soup.find_all("dd", {"itemprop": re.compile(r"^assignee", re.I)}):
    txt = dd.get_text(" ", strip=True)
    for tok in re.findall(r"[A-Za-z]{4,}", txt):
      if tok.lower() not in _corp_stop:
        assignees.add(tok.lower())
  for mc in soup.find_all("meta", {"name": "DC.contributor"}):
    c_val = (mc.get("content") or "").strip()
    if any(
        k in c_val.lower()
        for k in (
            "inc",
            "corp",
            "ltd",
            "llc",
            "pharma",
            "therapeutics",
            "biosciences",
        )
    ):
      for tok in re.findall(r"[A-Za-z]{4,}", c_val):
        if tok.lower() not in _corp_stop:
          assignees.add(tok.lower())

  results: List[Dict[str, Any]] = []
  seen_base = {re.sub(r"[A-Z]\d?$", "", used_id)}
  seen_ids = {used_id}
  family_base_ids: Set[str] = set()
  family_dates: Set[str] = set()

  # 1. Continuation / Divisional / Child / DOCDB Family
  for tr in soup.find_all(
      ["tr", "dd", "span", "div"],
      {
          "itemprop": re.compile(
              r"applications|afterApplications|childApplications|"
              r"parentApplications|continuation|docdbFamily",
              re.IGNORECASE,
          )
      },
  ):
    ip = tr.get("itemprop", "")
    txt = " ".join(tr.get_text(" ", strip=True).split())
    for m in re.finditer(r"\b((?:US\d{7,8}B[12])|(?:WO\d{10}A[123]))\b", txt):
      cand = m.group(1)
      base_cand = re.sub(r"[A-Z]\d?$", "", cand)
      family_base_ids.add(base_cand)
      if tr.name == "tr":
        row_dates = re.findall(r"\b(\d{4}-\d{2}-\d{2})\b", txt)
        if row_dates:
          family_dates.add(row_dates[0])
      if cand in seen_ids or base_cand in seen_base:
        continue
      seen_ids.add(cand)
      seen_base.add(base_cand)
      rel_type = (
          "Continuation / Divisional"
          if "docdb" not in ip.lower()
          else "DOCDB Family / PCT"
      )
      clean_title = re.sub(r"^.*?\(\s*en\s*\)\s*\*?\s*", "", txt).strip()
      results.append({
          "patent_id": cand,
          "seed_patent": used_id,
          "relationship": rel_type,
          "title": (clean_title or seed_title or "Family Counterpart")[:95],
          "score": 150 if cand[:2].upper() == used_id[:2].upper() else 140,
      })

  # Expand assignees and co-filed priority dates from citation rows that either
  # belong to the seed's DOCDB family/priority date or already match seed_tokens.
  all_ref_rows = soup.find_all(
      "tr",
      {
          "itemprop": re.compile(
              r"backwardReferences|forwardReferences|similarDocuments", re.I
          )
      },
  )
  co_filed_dates: Set[str] = set(family_dates)
  for tr in all_ref_rows:
    txt = " ".join(tr.get_text(" ", strip=True).split())
    m_dates = re.findall(r"\b(\d{4}-\d{2}-\d{2})\b", txt)
    if not m_dates:
      continue
    txt_low = txt.lower()
    m = re.search(r"\b(US\d{7,8}B[12]|WO\d{10}A[123])\b", txt)
    if m and re.sub(r"[A-Z]\d?$", "", m.group(1)) in family_base_ids:
      co_filed_dates.add(m_dates[0])
      continue
    if assignees and any(a in txt_low for a in assignees):
      t_part = re.sub(
          r"^.*?\d{4}-\d{2}-\d{2}\s+\d{4}-\d{2}-\d{2}\s+", "", txt_low
      )
      t_words = (
          set(re.findall(r"[a-z]{3,}", t_part))
          - assignees
          - _RELATED_STOP_WORDS
          - _corp_stop
      )
      if seed_tokens & t_words:
        co_filed_dates.add(m_dates[0])

  for tr in all_ref_rows:
    txt = " ".join(tr.get_text(" ", strip=True).split())
    m_dates = re.findall(r"\b(\d{4}-\d{2}-\d{2})\b", txt)
    if m_dates and m_dates[0] in co_filed_dates:
      m_ass = re.search(
          r"\d{4}-\d{2}-\d{2}\s+\d{4}-\d{2}-\d{2}\s+([A-Z][A-Za-z]+)", txt
      )
      if m_ass and m_ass.group(1).lower() not in _corp_stop:
        assignees.add(m_ass.group(1).lower())

  seed_tokens -= assignees | _corp_stop
  title_tokens -= assignees | _corp_stop

  # 2. Backward & Forward References scored by seed_title/abstract token overlap & assignee
  for tr in all_ref_rows:
    ip = tr.get("itemprop", "")
    txt = " ".join(tr.get_text(" ", strip=True).split())
    m = re.search(r"\b(US\d{7,8}B[12]|WO\d{10}A[123])\b", txt)
    if not m:
      continue
    cand = m.group(1)
    base_cand = re.sub(r"[A-Z]\d?$", "", cand)
    if cand in seen_ids or base_cand in seen_base:
      continue
    txt_low = txt.lower()
    title_part = re.sub(
        r"^.*?\d{4}-\d{2}-\d{2}\s+\d{4}-\d{2}-\d{2}\s+", "", txt_low
    )
    # Exclude polymorph, crystalline solid form, salt, formulation, combination, method-of-use, and biotransformation process patents
    if re.search(
        r"\b(?:polymorphs?|crystalline|solid\s+forms?|hydrate|maleate|salts?|formulations?|combinations?|biotransformation|enantiomer\s+separation|method\s+of\s+treating|for\s+use\s+in\s+the\s+treatment|use\s+of)\b",
        title_part,
    ):
      continue
    m_dates = re.findall(r"\b(\d{4}-\d{2}-\d{2})\b", txt)
    is_co_filed = bool(m_dates and m_dates[0] in co_filed_dates)
    same_ass = bool(assignees and any(a in txt_low for a in assignees))
    if assignees and not same_ass:
      continue
    title_words = {
        w
        for w in re.findall(r"[a-z]{3,}", title_part)
        if w not in _RELATED_STOP_WORDS
        and w not in assignees
        and w not in _corp_stop
    }
    shared = len(seed_tokens & title_words) + (1 if is_co_filed else 0)
    is_forward = "forward" in ip.lower()
    if not same_ass and shared < 2:
      continue
    if (
        same_ass
        and shared < 1
        and not is_forward
        and "similar" not in ip.lower()
    ):
      continue
    seen_ids.add(cand)
    seen_base.add(base_cand)
    rel_type = (
        "Parent / Prior Series (Cited)"
        if "backward" in ip.lower()
        else "Follow-On Series (Citing)"
    )
    clean_title = re.sub(r"^.*?\(\s*en\s*\)\s*\*?\s*", "", txt).strip()
    is_us_grant = bool(re.match(r"^US\d+B[12]$", cand))
    has_title_overlap = bool(title_tokens & title_words)
    results.append({
        "patent_id": cand,
        "seed_patent": used_id,
        "relationship": rel_type,
        "title": clean_title,
        "score": (
            (60 if same_ass else 45)
            + shared * 15
            + (15 if is_us_grant else 0)
            + (10 if has_title_overlap else 0)
        ),
    })

  results.sort(key=lambda r: (-int(r["score"]), str(r["patent_id"])))
  return results[:max_related]


def _parse_allowed_elements(
    spec_str: Optional[str], strict_organic: bool
) -> Optional[Set[int]]:
  """Parses a comma-separated element specification into atomic numbers."""
  if strict_organic:
    return set(ORGANIC_ELEMENTS)
  if not spec_str:
    return set(DEFAULT_ALLOWED_ELEMENTS)
  if spec_str.strip().lower() in ("all", "any", "none", "*"):
    return None
  pt = Chem.GetPeriodicTable()
  el_set = set()
  for tok in spec_str.split(","):
    t = tok.strip()
    if not t:
      continue
    if t.isdigit():
      el_set.add(int(t))
    else:
      try:
        el_set.add(pt.GetAtomicNumber(t.capitalize()))
      except (RuntimeError, ValueError, TypeError):
        pass
  return el_set or set(DEFAULT_ALLOWED_ELEMENTS)


def _build_arg_parser() -> argparse.ArgumentParser:
  """Constructs the CLI argument parser for patent_fetcher."""
  parser = argparse.ArgumentParser(
      description=(
          "Extract chemical structures associated with patents and their linked"
          " patent families via Google Patents, OPSIN, PubChem & ChEMBL."
      )
  )
  parser.add_argument(
      "positional_patents",
      nargs="*",
      default=[],
      help="Optional positional patent IDs (e.g. US-20130089624-A1)",
  )
  parser.add_argument(
      "--patent-id",
      "--patent_id",
      "--patent",
      dest="patent_id",
      required=False,
      default="",
      help=(
          "Patent ID (e.g. US-20130089624-A1 or US 8,497,274, or"
          " comma-separated list)"
      ),
  )
  parser.add_argument(
      "--include-family",
      action="store_true",
      default=True,
      help=(
          "Automatically discover and ingest continuation, DOCDB/PCT family,"
          " and same-assignee parent/follow-on patents (enabled by default)"
      ),
  )
  parser.add_argument(
      "--no-include-family",
      "--single-patent",
      dest="include_family",
      action="store_false",
      help="Disable automatic discovery of related family and series patents",
  )
  parser.add_argument(
      "--max-related-patents",
      type=int,
      default=3,
      help=(
          "Maximum additional related US/WO patents to harvest per seed patent"
          " (default: 3)"
      ),
  )
  parser.add_argument(
      "--include-concepts",
      "--include_concepts",
      dest="include_concepts",
      action="store_true",
      default=False,
      help=(
          "Harvest the structures Google Patents' own OCSR read out of the"
          " figures. Patent compounds are drawn rather than named, so this is"
          " usually the only route to the exemplified series. Rows are tagged"
          " Patent_Concept_OCSR, carry no example number, and need a"
          " substructure filter downstream."
      ),
  )
  parser.add_argument(
      "--concept-min-heavy-atoms",
      "--concept_min_heavy_atoms",
      dest="concept_min_heavy_atoms",
      type=int,
      default=6,
      help=(
          "Minimum heavy atoms for a concept fragment to be kept"
          " (default: 6). Only applies with --include-concepts."
      ),
  )
  parser.add_argument(
      "--limit",
      type=int,
      default=None,
      help=(
          "Maximum number of compounds to fetch (default: None, fetches all"
          " available compounds)"
      ),
  )
  parser.add_argument(
      "--source",
      choices=["all", "chembl", "pubchem"],
      default="all",
      help="Data source",
  )
  parser.add_argument(
      "--min-tanimoto",
      "--min_tanimoto",
      dest="min_tanimoto",
      type=float,
      default=0.40,
      help=(
          "Minimum Morgan ECFP4 Tanimoto similarity to the patent's exemplified"
          " series (unless sharing a Murcko scaffold) required to keep"
          " un-exemplified PubChem/OCSR records (default: 0.40)."
      ),
  )
  parser.add_argument(
      "--keep-all-pubchem-refs",
      "--keep_all_pubchem_refs",
      dest="keep_all_pubchem_refs",
      action="store_true",
      default=False,
      help=(
          "Retain all un-exemplified PubChem/OCSR records without applicability"
          " domain filtering."
      ),
  )
  parser.add_argument(
      "--query-smarts",
      "--query_smarts",
      dest="query_smarts",
      default=None,
      help="Optional SMARTS pattern to filter retrieved compounds.",
  )
  parser.add_argument(
      "--query-scaffold",
      "--query_scaffold",
      dest="query_scaffold",
      default=None,
      help="Optional SMILES/SMARTS scaffold to filter retrieved compounds.",
  )
  parser.add_argument(
      "--output-csv",
      "--output_csv",
      "--output",
      dest="output_csv",
      help="Optional output CSV path",
  )
  parser.add_argument(
      "--output-sdf",
      "--output_sdf",
      dest="output_sdf",
      help="Optional output SDF path",
  )
  parser.add_argument(
      "--allowed-elements",
      "--allowed_elements",
      dest="allowed_elements",
      default=None,
      help=(
          "Optional comma-separated list of allowed chemical element symbols or"
          " atomic numbers (e.g. 'C,H,N,O,F,Pt,Au' or 'all' to disable"
          " filtering). By default, all organic elements, bioisosteres (Si,"
          " Se), and metallodrugs are allowed."
      ),
  )
  parser.add_argument(
      "--strict-organic-only",
      "--strict_organic_only",
      dest="strict_organic_only",
      action="store_true",
      default=False,
      help=(
          "Restrict allowed elements strictly to canonical non-metal organic"
          " elements (exclude metallodrugs)."
      ),
  )
  return parser


def _run_cli(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
  """Executes the patent retrieval, deduplication, and export workflow."""
  eff_allowed_elements = _parse_allowed_elements(
      args.allowed_elements, args.strict_organic_only
  )

  raw_tokens: List[str] = []
  if args.patent_id:
    raw_tokens.extend(p.strip() for p in args.patent_id.split(",") if p.strip())
  for pos_p in args.positional_patents:
    raw_tokens.extend(p.strip() for p in str(pos_p).split(",") if p.strip())
  if not raw_tokens:
    parser.error(
        "Provide at least one patent ID via --patent-id or positional"
        " arguments."
    )

  _curate_script_path = os.path.abspath(
      os.path.join(
          os.path.dirname(__file__),
          "..",
          "..",
          "cheminformatics_core",
          "scripts",
          "curate_congeneric_corpus.py",
      )
  )

  seed_patent_ids: List[str] = []
  for tok in raw_tokens:
    comp = normalize_patent_id(tok)["compact"]
    if comp not in seed_patent_ids:
      seed_patent_ids.append(comp)

  patent_ids = list(seed_patent_ids)
  patent_relationship_map: Dict[str, str] = {
      pid: "Seed Patent" for pid in seed_patent_ids
  }
  patent_family_group: Dict[str, str] = {pid: pid for pid in seed_patent_ids}
  if args.include_family:
    related_rows: List[Dict[str, Any]] = []
    for pid in seed_patent_ids:
      rels = discover_related_patents(
          pid, max_related=max(args.max_related_patents + 2, 5)
      )
      related_rows.extend(rels)
      added_for_seed = 0
      seed_grp = patent_family_group.get(pid, pid)
      for r in rels:
        r_pid = str(r["patent_id"])
        rel_str = str(r.get("relationship") or "")
        is_same_doc_family = ("Continuation" in rel_str) or ("DOCDB" in rel_str)
        if is_same_doc_family:
          if r_pid in patent_family_group:
            old_grp = patent_family_group[r_pid]
            for k_p, v_g in list(patent_family_group.items()):
              if v_g == old_grp:
                patent_family_group[k_p] = seed_grp
          else:
            patent_family_group[r_pid] = seed_grp
        else:
          patent_family_group.setdefault(r_pid, r_pid)
        if r_pid not in patent_relationship_map:
          patent_relationship_map[r_pid] = f"{r['relationship']} (from {pid})"
        if (
            any(r_pid.startswith(pfx) for pfx in ("US", "WO", "EP"))
            and r_pid not in patent_ids
            and added_for_seed < args.max_related_patents
        ):
          patent_ids.append(r_pid)
          added_for_seed += 1
    if related_rows:
      print(
          "\n### Linked Patent Family, Parent & Follow-On Series Discovered\n"
      )
      print(
          "| Patent ID | Seed Patent | Relationship | Harvested | Title /"
          " Assignee Summary |"
      )
      print("| :--- | :--- | :--- | :--- | :--- |")
      seen_tbl = set()
      for r in related_rows:
        key_tbl = (r["patent_id"], r["seed_patent"])
        if key_tbl in seen_tbl:
          continue
        seen_tbl.add(key_tbl)
        harv = "Yes" if r["patent_id"] in patent_ids else "Reference Only"
        print(
            f"| `{r['patent_id']}` | `{r['seed_patent']}` |"
            f" {r['relationship']} | {harv} | {r['title']} |"
        )
      print("", flush=True)

  def _patent_harvest_priority(idx_pid: Tuple[int, str]) -> Tuple[int, int]:
    orig_idx, p_str = idx_pid
    if re.match(r"^US\d+B\d*$", p_str, re.IGNORECASE):
      return (0, orig_idx)
    if p_str.upper().startswith("US"):
      return (1, orig_idx)
    if p_str.upper().startswith("EP"):
      return (2, orig_idx)
    return (3, orig_idx)

  patent_ids = [
      p for _, p in sorted(enumerate(patent_ids), key=_patent_harvest_priority)
  ]

  all_fetched_records: List[Dict[str, Any]] = []
  seen_keys = set()
  by_structure: Dict[str, Dict[str, Any]] = {}
  structureless_seen = 0
  concept_merged = 0
  any_html_assays_found = False
  for pid in patent_ids:
    eff_limit = args.limit
    print(
        f"[*] Querying public patent repositories for: {pid}"
        f" ({patent_relationship_map.get(pid, 'Related Patent')})"
    )
    recs = fetch_patent_compounds(
        pid,
        limit=eff_limit,
        preferred_source=args.source,
        min_tanimoto=args.min_tanimoto,
        keep_all_pubchem_refs=args.keep_all_pubchem_refs,
        query_smarts=args.query_smarts,
        query_scaffold=args.query_scaffold,
        allowed_elements=eff_allowed_elements,
        skip_external_bioassays=any_html_assays_found,
    )
    if any(
        r.get("patent_example_id")
        and (
            r.get("pchembl_value") is not None
            or r.get("activity_value") is not None
            or r.get("activity_bin") is not None
        )
        for r in recs
    ):
      any_html_assays_found = True
    for r in recs:
      r["patent_relationship"] = patent_relationship_map.get(
          pid, "Related Patent"
      )
      ex_id = r.get("patent_example_id") or r.get("compound_id") or ""
      structure_key = r.get("inchikey") or r.get("canonical_smiles")
      if not structure_key:
        # Assay-only rows resolve to no structure at all, so there is nothing
        # to compare them on. Give each one a unique key instead of letting
        # every structureless row collapse onto a shared empty one.
        structureless_seen += 1
        structure_key = f"__no_structure__:{structureless_seen}"
      else:
        by_structure.setdefault(
            f"{structure_key.split('-')[0]}|{r.get('patent_id')}", r
        )
      peak_key = str(r.get("chiral_peak") or "")
      key = (
          (structure_key, str(ex_id), peak_key, str(r.get("patent_id")))
          if (ex_id or peak_key)
          else (structure_key, str(r.get("patent_id")))
      )
      if key not in seen_keys:
        seen_keys.add(key)
        all_fetched_records.append(r)

    if args.include_concepts:
      for c in fetch_patent_concepts(
          pid, min_heavy_atoms=args.concept_min_heavy_atoms
      ):
        ik = c.get("inchikey") or c.get("canonical_smiles")
        if not ik:
          # Markush rows have no structure by definition and must not be
          # folded into anything.
          all_fetched_records.append(c)
          continue
        skel = f"{ik.split('-')[0]}|{c.get('patent_id')}"
        prior = by_structure.get(skel)
        if prior is not None and prior is not c:
          prior["also_patent_concept"] = True
          prior["all_sources"] = ",".join(
              sorted(
                  (
                      set(str(prior.get("all_sources") or "").split(","))
                      | {"Patent_Concept_OCSR"}
                  )
                  - {""}
              )
          )
          prior["concept_count"] = c.get("concept_count")
          prior["concept_sections"] = c.get("concept_sections")
          concept_merged += 1
          continue
        by_structure[skel] = c
        all_fetched_records.append(c)

  # Cross-family backfill: if a continuation/divisional patent produced a
  # structureless Example_X row whose Example_X was resolved with a valid
  # chemical structure in a parent/sibling patent IN THE SAME PATENT FAMILY,
  # merge any missing assay endpoints onto the structured row and omit structureless table rows.
  structured_by_ex: Dict[str, Dict[str, Any]] = {}
  has_any_structured = False
  for r in all_fetched_records:
    if r.get("canonical_smiles"):
      has_any_structured = True
      ex_id = str(r.get("patent_example_id") or "").strip()
      if ex_id:
        structured_by_ex.setdefault(ex_id, r)
        r_pat = str(r.get("patent_id") or "").strip()
        fam_grp = patent_family_group.get(r_pat, r_pat)
        bare_ex = ex_id.split(":", 1)[-1]
        if fam_grp and bare_ex:
          structured_by_ex.setdefault(f"{fam_grp}::{bare_ex}", r)

  if has_any_structured:
    cleaned_fetched: List[Dict[str, Any]] = []
    omitted_structureless = 0
    for r in all_fetched_records:
      if (
          not r.get("canonical_smiles")
          and r.get("source") == "Patent_FullText_Table"
      ):
        ex_id = str(r.get("patent_example_id") or "").strip()
        r_pat = str(r.get("patent_id") or "").strip()
        fam_grp = patent_family_group.get(r_pat, r_pat)
        bare_ex = ex_id.split(":", 1)[-1]
        target_r = structured_by_ex.get(ex_id) or (
            structured_by_ex.get(f"{fam_grp}::{bare_ex}")
            if (fam_grp and bare_ex)
            else None
        )
        if target_r is not None:
          for k_col, v_col in r.items():
            if v_col is not None and target_r.get(k_col) is None:
              target_r[k_col] = v_col
        omitted_structureless += 1
        continue
      cleaned_fetched.append(r)
    if omitted_structureless > 0:
      print(
          "[*] Cross-family deduplication: merged/omitted"
          f" {omitted_structureless} unresolved structureless table rows."
      )
    all_fetched_records = cleaned_fetched

  records = filter_records_by_exemplified_domain(
      all_fetched_records,
      min_tanimoto=args.min_tanimoto,
      keep_all_refs=args.keep_all_pubchem_refs,
      query_smarts=args.query_smarts,
      query_scaffold=args.query_scaffold,
      prune_unexemplified_unassayed=True,
  )
  if concept_merged:
    print(
        f"[*] {concept_merged} concept structures were already present from"
        " ChEMBL/PubChem; those rows are flagged also_patent_concept=True"
        " rather than duplicated."
    )

  macro_count = sum(1 for r in records if r.get("is_macrocycle"))
  act_count = sum(
      1
      for r in records
      if r.get("pchembl_value") is not None
      or r.get("activity_value") is not None
      or r.get("activity_bin") is not None
  )
  print(
      f"[+] Successfully extracted {len(records)} validated molecules"
      f" ({macro_count} macrocycles, {act_count} with bioactivity data)."
  )

  # Provenance breakdown. PubChem rows are compounds ASSOCIATED with the
  # patent document -- cited prior art, comparators and reagents indexed from
  # the experimental section -- not the exemplified series. Reporting only a
  # combined total lets that distinction disappear.
  by_source: Dict[str, int] = {}
  for r in records:
    by_source[str(r.get("source", "?"))] = (
        by_source.get(str(r.get("source", "?")), 0) + 1
    )
  print("[*] Provenance breakdown:")
  for src, n in sorted(by_source.items(), key=lambda kv: -kv[1]):
    print(f"      {src:28s} {n}")
  assoc = by_source.get("PubChem", 0)
  concepts = by_source.get("Patent_Concept_OCSR", 0)
  # Count every row the patent's own example table contributed to, not just
  # rows whose only source is that table. Most exemplified compounds are also
  # returned by PubChem, so they carry source="PubChem" and would otherwise be
  # miscounted as associative.
  exemplified = sum(
      1
      for r in records
      if "Patent_FullText_Table" in str(r.get("all_sources", ""))
  )
  paired = sum(
      1
      for r in records
      if r.get("patent_example_id")
      and (
          r.get("pchembl_value") is not None
          or r.get("activity_value") is not None
          or r.get("activity_bin") is not None
      )
  )
  if exemplified:
    print(
        f"[*] {exemplified} rows carry data from the patent's own example"
        f" tables, of which {paired} pair an example number with a potency."
    )
  if assoc:
    print(
        f"[!] {assoc} of {len(records)} rows came from PubChem's patent"
        " association index. These are compounds the patent DOCUMENT"
        " references (prior art, comparators, reagents, intermediates), not"
        " the exemplified series. Filter structurally before treating this"
        " as a congeneric series -- see"
        " cheminformatics_core/scripts/"
        "library_search.py, which takes --scaffold_survey to find the cores"
        " present and then one or more --query_smarts to select them."
    )
  if concepts:
    print(
        f"[!] {concepts} rows are individual fragments of Google's OCSR read"
        " of the patent figures. Each figure is a whole reaction scheme, so"
        " starting materials, reagents and intermediates are in there"
        " alongside the products. This set is a candidate pool, not a series."
    )
  if exemplified == 0 and concepts == 0:
    print(
        "[!] ZERO compounds were recovered from the patent's own example"
        " tables. Everything above is associative. Treat this result as"
        " unverified for series membership. If the patent draws its compounds"
        " rather than naming them, which is the common case, retry with"
        " --include-concepts."
    )

  if not records:
    print(
        "[-] No chemical structures found for this patent in public databases."
    )
    sys.exit(0)

  # Print summary of first few
  for i, r in enumerate(records[:5]):
    smi = r.get("canonical_smiles") or ""
    shown = f"{smi[:60]}..." if smi else "(no structure resolved)"
    print(f"  [{i+1}] {r['compound_id']} ({r['source']}): {shown}")

  if args.output_csv:
    csv_rows = []
    core_cols = {
        "patent_id",
        "source",
        "document_id",
        "document_year",
        "publication_date",
        "filing_date",
        "priority_date",
        "compound_id",
        "patent_example_id",
        "chiral_peak",
        "chembl_id",
        "pubchem_cid",
        "compound_name",
        "iupac_name",
        "canonical_smiles",
        "inchikey",
        "dataset_role",
        "is_macrocycle",
        "max_ring_size",
        "activity_type",
        "activity_value",
        "activity_units",
        "pchembl_value",
        "activity_bin",
        "activity_range",
        "target_name",
        "stoichiometry_multiplier",
        "molecular_formula",
        "patent_relationship",
        "all_sources",
        "max_tanimoto_to_exemplified",
        "also_patent_concept",
        "concept_count",
        "concept_sections",
    }
    for r in records:
      has_label = (
          r.get("pchembl_value") is not None
          or r.get("activity_value") is not None
          or r.get("activity_bin") is not None
      )
      row_dict = {
          "patent_id": r.get("patent_id", ""),
          "source": r.get("source", ""),
          "document_id": r.get("document_id", r.get("patent_id", "")),
          "document_year": r.get("document_year"),
          "publication_date": r.get("publication_date"),
          "filing_date": r.get("filing_date"),
          "priority_date": r.get("priority_date"),
          "compound_id": r.get("compound_id", ""),
          "chiral_peak": r.get("chiral_peak", ""),
          "chembl_id": r.get("chembl_id", ""),
          "pubchem_cid": r.get("pubchem_cid", ""),
          "compound_name": r.get("compound_name", r.get("iupac_name", "")),
          "canonical_smiles": r.get("canonical_smiles", ""),
          "inchikey": r.get("inchikey", ""),
          "dataset_role": (
              "labeled_patent" if has_label else "prospective_patent"
          ),
          "is_macrocycle": r.get("is_macrocycle", False),
          "max_ring_size": r.get("max_ring_size", 0),
          "activity_type": r.get("activity_type"),
          "activity_value": r.get("activity_value"),
          "activity_units": r.get("activity_units"),
          "pchembl_value": r.get("pchembl_value"),
          "activity_bin": r.get("activity_bin"),
          "activity_range": r.get("activity_range"),
          "target_name": r.get("target_name"),
          "stoichiometry_multiplier": 1.0,
          "molecular_formula": r.get("molecular_formula"),
      }
      for k, v in r.items():
        if (
            k not in row_dict
            and k not in ("rdkit_mol", "smiles")
            and not k.startswith("_")
            and not isinstance(v, (dict, list))
        ):
          row_dict[k] = v
      csv_rows.append(row_dict)
    df = pd.DataFrame(csv_rows)
    # Sparsity-based pruning of dynamic ChEMBL all_endpoints columns (< max(3, 2% of rows))
    min_non_null = max(3, int(0.02 * len(df))) if len(df) >= 5 else 1
    keep_cols = []
    for col in df.columns:
      if (
          col in core_cols
          or col.startswith("patent_table_col_")
          or col.startswith("patent_assay_")
          or int(df[col].notna().sum()) >= min_non_null
      ):
        keep_cols.append(col)
    df = df[keep_cols]
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True
    )
    df.to_csv(args.output_csv, index=False)
    print(
        f"[+] Saved CSV table to: {args.output_csv} (columns:"
        f" {len(df.columns)})"
    )
    print(
        "[NEXT STEP] Triage scaffolds with `library_search.py"
        " --scaffold_survey`, then harvest ChEMBL/PubMed literature & curate"
        f" the master series with:\n            uv run {_curate_script_path}"
        f" --chembl-docs <CHEMBL_IDs_or_PMIDs> --input-csv {args.output_csv}",
        flush=True,
    )

  if args.output_sdf:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_sdf)), exist_ok=True
    )
    with Chem.SDWriter(args.output_sdf) as writer:
      skipped_no_mol = 0
      for r in records:
        m = r.get("rdkit_mol")
        if m is None:
          # Assay-only example rows carry no structure and cannot be written to
          # an SDF. Skipping them silently would misreport the export size.
          skipped_no_mol += 1
          continue
        m.SetProp("_Name", str(r["compound_id"]))
        m.SetProp("PatentID", str(r["patent_id"]))
        m.SetProp("Source", str(r["source"]))
        m.SetProp("DocumentTitle", str(r.get("document_title", "")))
        writer.write(m)
    print(f"[+] Saved SDF molecular library to: {args.output_sdf}")
    if skipped_no_mol:
      print(
          f"[!] {skipped_no_mol} records had no resolved structure and were"
          " omitted from the SDF (they remain in the CSV)."
      )


def main() -> None:
  parser = _build_arg_parser()
  args = parser.parse_args()
  _run_cli(args, parser)


if __name__ == "__main__":
  main()
