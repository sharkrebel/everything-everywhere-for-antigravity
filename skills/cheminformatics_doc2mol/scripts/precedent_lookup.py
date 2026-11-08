# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "rdkit>=2024.3.1",
#   "pandas>=2.2.0",
#   "numpy>=1.26.0",
#   "drfp>=0.3.7",
#   "pyarrow>=14.0.0",
#   "ord-schema>=0.3.8",
#   "polite-http",
# ]
# ///

"""Reaction Precedent Lookup Engine.

Searches reaction precedent databases (ORD/USPTO corpora via the pre-indexed
ORD v2 API, PubChem Patent cross-references, and full-text Google Patents
specifications) via local datasets or online APIs, computes vectorized batch
DRFP reaction fingerprint similarity, and clusters reaction condition
archetypes with empirical yield and reaction time statistics.
"""

from __future__ import annotations

import argparse
import base64
from collections.abc import Mapping, Sequence, Set as AbstractSet
import concurrent.futures
import dataclasses
import html as html_lib
from html import parser as html_parser
import http.client
import json
import math
import os
import re
import sys
import threading
from typing import Any, Final
import urllib.parse
import zlib

import chemical_text_sanitizer
from drfp import DrfpEncoder
from drfp import fingerprint as drfp_fingerprint
from google.protobuf import message as protobuf_message
import numpy as np
from ord_schema.proto import reaction_pb2
import pandas as pd
from polite_http import http_client
import procedure_extractor
import pyarrow as pa
from rdkit import Chem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import rdFMCS

RDLogger.DisableLog("rdApp.*")

_ORD_CLIENT = http_client.HttpClient(
    "https://open-reaction-database.org",
    qps=5.0,
    timeout=12.0,
)
_PUBCHEM_CLIENT = http_client.HttpClient(
    "https://pubchem.ncbi.nlm.nih.gov/",
    qps=5.0,
    timeout=10.0,
    max_retries=1,
)
_GOOGLE_PATENTS_CLIENT = http_client.HttpClient(
    "https://patents.google.com/",
    qps=8.0,
    timeout=8.0,
    max_retries=1,
)

_ORD_BATCH_FETCH_MAX: Final[int] = 100  # Max IDs per /api/reactions request.
_ORD_SIMILARITY_THRESHOLD: Final[float] = 0.65  # Default Morgan cutoff.
_ORD_SMALL_MOL_SIMILARITY_THRESHOLD: Final[float] = 0.80  # Cutoff (<8 atoms).
_ORD_SMALL_MOL_HEAVY_ATOMS: Final[int] = 8  # Heavy-atom cutoff for anchors.
_ORD_SEARCH_TIMEOUT_SECS: Final[float] = 6.0  # Per-pass /api/search timeout.
_ORD_ORDER_BY = [{"key": "reaction_id"}]
_PUBCHEM_PATENT_XREF_TIMEOUT_SECS: Final[float] = 25.0
_GOOGLE_PATENTS_HTML_TIMEOUT_SECS: Final[float] = 8.0
_MAX_PATENT_XREF_REACTANTS: Final[int] = 2
_MAX_PATENT_HTML_WORKERS: Final[int] = 6
_MAX_PATENT_PROC_PARSE_CANDIDATES: Final[int] = 8
_MAX_BATCH_STEP_WORKERS: Final[int] = 3
_MIN_ANALOG_PRODUCT_MORGAN_SIM: Final[float] = 0.35
_UNSPECIFIED_STEREO_PRODUCT_SIM: Final[float] = 0.90
_CONFLICTING_STEREO_PRODUCT_SIM_CAP: Final[float] = 0.80
# Clears default min_similarity=0.35 when the exact query product is synthesized
# from unattested starting materials, while staying bounded below 1-of-N DRFP.
_PRODUCT_ONLY_PRECEDENT_WEIGHT: Final[float] = 0.36
_PRODUCT_ONLY_DRFP_CAP_RATIO: Final[float] = 0.85
_REACTANT_ONLY_PATENT_WEIGHT: Final[float] = 0.50
_MCS_TIMEOUT_SECS: Final[int] = 1

_MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=2048,
)
_CHIRAL_MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=2048,
    includeChirality=True,
)
_MOL_SANITIZE_ERRORS: Final = (
    Chem.rdchem.MolSanitizeException,
    ValueError,
    RuntimeError,
)
_HTTP_FETCH_ERRORS = (
    http_client.HttpError,
    http.client.HTTPException,
    OSError,
    EOFError,
    zlib.error,
    LookupError,
    ValueError,
)
_DRFP_ENCODE_ERRORS = (
    drfp_fingerprint.NoReactionError,
    ValueError,
    RuntimeError,
    TypeError,
)

_PATENT_SEP_TRANSLATION = str.maketrans("", "", "-_ /")
_PCT_TOKEN_RE = re.compile(r"^PCT/[A-Z]{2}(\d{2,4})/(\d+)", re.IGNORECASE)
_WO_HYPHEN_RE = re.compile(
    r"^WO[-_ /]+(\d{2}|\d{4})[-_ /]+(\d{1,6})(?:[-_ /]*([A-Z]\d?))?$",
    re.IGNORECASE,
)
_WO_4DIGIT_RE = re.compile(r"^WO((?:19|20)\d{2})(\d{5,6})([A-Z]\d?)?$")
_WO_2DIGIT_RE = re.compile(r"^WO(\d{2})(\d{5,6})([A-Z]\d?)?$")
_US_PREGRANT_KIND_RE = re.compile(r"^US(20\d{2})(\d{1,6})([AP]\d)$")
_US_PREGRANT_STEM_RE = re.compile(r"^US(20\d{2})(\d{6})$")
_PATENT_KIND_SUFFIX_RE = re.compile(r"^([A-Z]{2,4}\d+)[A-Z][A-Z0-9]?$")
_PATENT_YEAR_RE = re.compile(r"^[A-Z]{2,4}((?:19|20)\d{2})\d{4,}")
_TEMP_RANGE_C_RE = re.compile(
    r"(?<![\d.\-–—])(-?\d{1,3}(?:\.\d+)?)\s*(?:°\s*C\b|℃)?\s*"
    r"(?:[\-–—~]|to\b)\s*(-?\d{1,3}(?:\.\d+)?)\s*(?:°\s*C\b|℃)",
    re.IGNORECASE,
)
_TEMP_C_RE = re.compile(
    r"(?<![\d.\-–—])(-?\d{1,3}(?:\.\d+)?)\s*(?:°\s*C\b|℃)",
    re.IGNORECASE,
)
_ROOM_TEMP_RE = re.compile(
    r"\b(?:room\s+temperature|r\.t\.|ambient\s+temperature)\b",
    re.IGNORECASE,
)
_TIME_RE = re.compile(
    r"(?<![A-Za-z\d.])"
    r"(?:(\d{1,3}(?:\.\d+)?)\s*(?:[\-–—~]|(?i:to)\b)\s*)?"
    r"(\d{1,3}(?:\.\d+)?)"
    r"(?:\s+(?i:(hr|hrs|hours|min|mins|minutes))\b|"
    r"\s*(h)\b(?!\s*[\-–]?\s*(?i:nmr)\b|[\-–][A-Za-z]))"
)
_OVERNIGHT_RE = re.compile(r"\bovernight\b", re.IGNORECASE)
_EXPERIMENTAL_HINT_RE = re.compile(
    r"\b(?:mmol|mol\b|mg\b|g\b|mL\b|yield|afford|isolat|stir|heat|"
    r"room\s+temperature|r\.t\.|overnight|example|embodiment|preparation|"
    r"synthesis|step\s+\d)\b|°\s*C|℃",
    re.IGNORECASE,
)
_EXAMPLE_ID_RE = re.compile(
    r"\b((?:Example|Embodiment|Preparation|Step|Table)\s+\d+[A-Za-z]?)\b",
    re.IGNORECASE,
)
_CATALYST_RE = re.compile(
    r"\bcatalyst\s+([A-Za-z0-9\-\(\)\[\],\s]{3,45}?)"
    r"(?:\s*\(|\s*,|\s+was\b|\s+and\b|\s+in\b)",
    re.IGNORECASE,
)
_EXPLICIT_YIELD_PCT_RE = re.compile(
    r"(?:yield[^0-9%]{0,20}([0-9]+(?:\.[0-9]+)?)\s*%|"
    r"([0-9]+(?:\.[0-9]+)?)\s*%\s*(?:isolated\s+|isolation\s+)?yield)",
    re.IGNORECASE,
)
_PRODUCT_CUE_RE = re.compile(
    r"\b(?:afford\w*|giv(?:e|es|en|ing)|gave|yield\w*|obtain\w*|provid\w*|"
    r"furnish\w*|produc\w*|isolat\w*|synthesi[sz]\w*\s+(?:of\s+)?|"
    r"prepar\w*\s+(?:of\s+)?|title\s+compound)\b",
    re.IGNORECASE,
)
_PRODUCT_AFTER_CUE_RE = re.compile(
    r"^[^.;]{0,65}\b(?:yield\b|was\s+(?:obtained|isolated|afforded|prepared|"
    r"synthesized)\b)",
    re.IGNORECASE,
)
_REACTANT_CHARGE_PREFIX_RE = re.compile(
    r"\b(?:to\s+(?:a\s+)?(?:stirred\s+)?(?:solution|mixture|suspension|slurry)"
    r"|(?:solution|mixture|suspension|slurry)\s+of|into|charged\s+with|using"
    r"|adding|charging|mixing|dissolving|weighing|taking|treating"
    r"|reacting)\s*$",
    re.IGNORECASE,
)
_REACTANT_CHARGE_QTY_RE = re.compile(
    r"^\s*(?:\d+(?:\.\d+)?\s*(?:mg|g|mL)\s*,?\s*)?"
    r"\(?[^);.]{0,35}\b(?:mmol|mol|equiv|eq)\b",
    re.IGNORECASE,
)
_REACTANT_ACTION_SUFFIX_RE = re.compile(
    r"^[^.;]{0,80}\b(?:was|were)\s+"
    r"(?:added|dissolved|stirred|heated|cooled|treated|charged|mixed)\b",
    re.IGNORECASE,
)
_HEADING_LINE_RE = re.compile(
    r"^(?:Example|Embodiment|Preparation|Step|Synthesis)\b",
    re.IGNORECASE,
)
_HTML_TAG_STRIP_RE = re.compile(r"<[^>]+>")
_HEAD_END_RE = re.compile(r"</head\s*>", re.IGNORECASE)
_DESC_SECTION_START_RE = re.compile(
    r"<section\b[^>]*\bitemprop=['\"]description['\"][^>]*>",
    re.IGNORECASE,
)
_SECTION_TAG_RE = re.compile(r"<(/?)section\b[^>]*>", re.IGNORECASE)
_CLAUSE_SPLIT_RE = re.compile(r"[.;]\s+")
_BLOCK_LEAF_TAGS: Final[frozenset[str]] = frozenset({
    "p",
    "tr",
    "heading",
    "h1",
    "h2",
    "h3",
    "h4",
    "h5",
    "h6",
})
_SPACE_BOUNDARY_TAGS: Final[frozenset[str]] = frozenset({
    "td",
    "th",
    "br",
    "li",
    "dd",
    "dt",
    "para-num",
})

_CACHE_LOCK = threading.Lock()
_ERROR_LOCK = threading.Lock()
_PUBCHEM_INFO_BY_SMILES: dict[str, dict[str, Any]] = {}
_PATENT_IDS_BY_CID: dict[int, tuple[str, ...]] = {}
_PATENT_XREF_ERROR_BY_CID: dict[int, str] = {}
_HTML_BY_PATENT_ID: dict[str, str] = {}
_HTML_ERROR_BY_PATENT_ID: dict[str, str] = {}
_RESOLVED_SMILES_BY_PRODUCT_NAME: dict[str, str] = {}

_CHEMICAL_NAME_MARKER_RE = re.compile(r"[\d\-\(\)\[\]]")
_GENERIC_PLACEHOLDER_NAME_RE = re.compile(
    r"^(?:compound|intermediate|example|embodiment|preparation|step|formula|"
    r"product|title\s+compound|isomer|enantiomer|diastereomer|fraction|peak)"
    r"\b[\s\dA-Za-z\-\(\)\[\]]*$",
    re.IGNORECASE,
)
_NON_ALNUM_RE = re.compile(r"[^a-z0-9]+")
_SALT_WORD_STRIP_RE = re.compile(
    r"\b(?:mono|di|tri|hemi|bis|tris)?(?:hydrochloride|hydrobromide|"
    r"hydroiodide|trifluoroacetate|trifluoroacetic\s+acid|formate|acetate|"
    r"mesylate|methanesulfonate|tosylate|p\-toluenesulfonate|besylate|"
    r"benzenesulfonate|tartrate|citrate|fumarate|maleate|sulfate|bisulfate|"
    r"phosphate|hexafluorophosphate|tetrafluoroborate|nitrate|sodium|"
    r"potassium|hydrate)\b",
    re.IGNORECASE,
)
_STEM_REQUIRED_ELEMENTS: Final[tuple[tuple[re.Pattern[str], str], ...]] = (
    (re.compile(r"fluoro", re.IGNORECASE), "F"),
    (re.compile(r"chloro", re.IGNORECASE), "Cl"),
    (re.compile(r"bromo", re.IGNORECASE), "Br"),
    (re.compile(r"iodo", re.IGNORECASE), "I"),
    (re.compile(r"(?:cyano|nitril|carbonitril|nitro)", re.IGNORECASE), "N"),
    (re.compile(r"(?:sulfonyl|sulfonamido|thio|mercapto)", re.IGNORECASE), "S"),
    (re.compile(r"phospho", re.IGNORECASE), "P"),
)
_MS_HEADER_WINDOW_RE = re.compile(
    r"(?:\bm\s*/\s*z\b|\b(?:LC[\-\s]?MS|HRMS|LRMS|ESI(?:\-MS)?|APCI)\b)"
    r"([^.;\n]{1,80})",
    re.IGNORECASE,
)
_MS_STANDALONE_MZ_RE = re.compile(
    r"(?<![A-Za-z\d.])(\d{2,4}(?:\.\d{1,4})?)"
    r"(?!\s*(?:MHz|Hz|mg|g|mL|mmol|mol|%|°\s*C|℃|min|h|hr|hrs)\b)",
    re.IGNORECASE,
)
_MS_ADDUCT_AFTER_MZ_RE = re.compile(
    r"(?<![A-Za-z\d.])(\d{2,4}(?:\.\d{1,4})?)\s*"
    r"(?:\[|\()\s*M\s*[\+\-]\s*(?:H|Na|2H|NH4)?[^\)\]]*(?:\]|\))",
    re.IGNORECASE,
)
_MS_ADDUCT_BEFORE_MZ_RE = re.compile(
    r"(?:\[|\()\s*M\s*[\+\-]\s*(?:H|Na|2H|NH4)?[^\)\]]*(?:\]|\))"
    r"[^.;\n\d]{0,25}(\d{2,4}(?:\.\d{1,4})?)",
    re.IGNORECASE,
)


# ==============================================================================
# Chemical Data Models
# ==============================================================================


@dataclasses.dataclass
class PrecedentRecord:
  """Historical reaction precedent with reaction conditions and citation."""

  reaction_id: str
  reaction_smiles: str
  reactants: list[str]
  product: str
  reagents: str
  solvent: str
  temperature_c: float | None
  time_hours: float | None
  yield_percent: float | None
  citation: str
  example_id: str
  similarity_score: float = 0.0
  title: str = ""
  url: str = ""
  notes: str = ""
  source: str = ""

  def to_dict(self) -> dict[str, Any]:
    return dataclasses.asdict(self)


@dataclasses.dataclass
class ConditionCluster:
  """Aggregated statistical condition archetype across retrieved precedents."""

  cluster_name: str
  reagents: str
  solvent: str
  temp_bracket: str
  precedent_count: int
  median_yield: float | None
  yield_range: tuple[float, float] | None
  median_time_hours: float | None
  time_range_hours: tuple[float, float] | None
  temperature_range_c: tuple[float, float] | None
  citations: list[str]

  def to_dict(self) -> dict[str, Any]:
    return dataclasses.asdict(self)


@dataclasses.dataclass
class PatentCandidate:
  """Candidate patent ranked by empirical IDF or synthesis co-occurrence."""

  patent_id: str
  idf_score: float
  idf_coverage: float
  matched_reactant_indices: list[int]
  matched_product_indices: list[int]
  xhr_title: str = ""
  xhr_snippet: str = ""


# ==============================================================================
# SMILES Canonicalization & Reaction Normalization
# ==============================================================================


def _canonicalize_valid_smiles(smi: str) -> str:
  """Returns canonical unmapped RDKit SMILES, or '' if smi cannot be parsed."""
  if not smi or not smi.strip():
    return ""
  without_cx = smi.split("|", 1)[0].strip()
  if not without_cx:
    return ""
  token = without_cx.split()[0]
  mol = Chem.MolFromSmiles(token)
  if mol is None:
    return ""
  if ":" in token:
    for atom in mol.GetAtoms():
      if atom.GetAtomMapNum():
        atom.SetAtomMapNum(0)
  return Chem.MolToSmiles(mol)


def canonicalize_smiles(smi: str) -> str:
  """Converts SMILES string to canonical unmapped RDKit SMILES."""
  if not smi or not smi.strip():
    return ""
  return _canonicalize_valid_smiles(smi) or smi.strip()


def _split_smiles_fragments(smi: str) -> list[str]:
  """Strips CXSMILES '|...|' extensions before splitting SMILES fragments."""
  if not smi or not smi.strip():
    return []
  without_cx = smi.split("|", 1)[0].strip()
  if not without_cx:
    return []
  core = without_cx.split()[0]
  return [s.strip() for s in core.split(".") if s.strip()]


def _split_reaction_sides(rxn_smiles: str) -> tuple[list[str], list[str]]:
  """Extracts reactant and product tokens from 2- or 3-part reaction SMILES."""
  if not rxn_smiles or not rxn_smiles.strip():
    return [], []
  without_cx = rxn_smiles.split("|", 1)[0].strip()
  compact = re.sub(r"\s*(>>|>|\.)\s*", r"\1", without_cx)
  if not compact:
    return [], []
  core = compact.split()[0]
  if ">>" in core:
    lhs, rhs = core.split(">>", 1)
  elif core.count(">") == 2:
    lhs, _, rhs = core.split(">")
  else:
    return [s for s in core.split(".") if s], []
  reacts = [s for s in lhs.split(".") if s]
  prods = [s for s in rhs.split(".") if s]
  return reacts, prods


def normalize_reaction_smiles_for_drfp(rxn_smiles: str) -> str:
  """Normalizes reaction SMILES to unmapped canonical 'reactants>>products'."""
  if not rxn_smiles or ">" not in rxn_smiles:
    return ""
  lhs_parts, rhs_parts = _split_reaction_sides(rxn_smiles)
  if not lhs_parts or not rhs_parts:
    return ""
  lhs_clean = [_canonicalize_valid_smiles(s) for s in lhs_parts]
  rhs_clean = [_canonicalize_valid_smiles(s) for s in rhs_parts]
  if not all(lhs_clean) or not all(rhs_clean):
    return ""
  return f"{'.'.join(sorted(lhs_clean))}>>{'.'.join(sorted(rhs_clean))}"


def _molecular_complexity_score(smi: str) -> float:
  """Computes Bertz topological complexity for closed-shell organic SMILES."""
  if not smi or "*" in smi:
    return -1.0
  clean_smi = _canonicalize_valid_smiles(smi)
  if not clean_smi:
    return -1.0
  mol = Chem.MolFromSmiles(clean_smi)
  if mol is None or mol.GetNumHeavyAtoms() < 2:
    return -1.0
  for atom in mol.GetAtoms():
    z = atom.GetAtomicNum()
    if z == 0 or atom.GetNumRadicalElectrons() > 0:
      return -1.0
    if (21 <= z <= 30) or (39 <= z <= 48) or (72 <= z <= 80):
      return -1.0
  return max(float(Descriptors.BertzCT(mol)), float(mol.GetNumHeavyAtoms()))


def _is_searchable_organic_smiles(smi: str) -> bool:
  """Returns True for closed-shell, non-dummy, non-transition-metal SMILES."""
  return _molecular_complexity_score(smi) > 0.0


def _rank_organic_smiles(raw_smiles: Sequence[str]) -> list[str]:
  """Deduplicates, canonicalizes, and ranks organic SMILES by complexity."""
  unique_clean: list[str] = []
  seen: set[str] = set()
  for item in raw_smiles:
    for part in _split_smiles_fragments(str(item)):
      smi = _canonicalize_valid_smiles(part)
      if smi and smi not in seen and _is_searchable_organic_smiles(smi):
        seen.add(smi)
        unique_clean.append(smi)
  return sorted(unique_clean, key=_molecular_complexity_score, reverse=True)


# ==============================================================================
# Database Ingestion & Live Pre-Indexed ORD v2 API Client
# ==============================================================================


def _ord_time_to_hours(time_msg: Any) -> float | None:
  """Converts an ORD Time message into hours."""
  if not time_msg or not time_msg.HasField("value"):
    return None
  val = float(time_msg.value)
  unit = int(time_msg.units)
  if unit == reaction_pb2.Time.MINUTE:
    return round(val / 60.0, 4)
  if unit == reaction_pb2.Time.SECOND:
    return round(val / 3600.0, 4)
  if unit == reaction_pb2.Time.DAY:
    return val * 24.0
  return val


def _ord_temperature_to_celsius(temp_msg: Any) -> float | None:
  """Converts an ORD ReactionTemperature message into degrees Celsius."""
  if not temp_msg:
    return None
  if temp_msg.HasField("setpoint") and temp_msg.setpoint.HasField("value"):
    val = float(temp_msg.setpoint.value)
    unit = int(temp_msg.setpoint.units)
    if unit == reaction_pb2.Temperature.FAHRENHEIT:
      return round((val - 32.0) * 5.0 / 9.0, 2)
    if unit == reaction_pb2.Temperature.KELVIN:
      return round(val - 273.15, 2)
    return val
  if temp_msg.HasField("control"):
    ctrl_type = int(temp_msg.control.type)
    if (
        ctrl_type
        == reaction_pb2.TemperatureConditions.TemperatureControl.AMBIENT
    ):
      return 25.0
    if (
        ctrl_type
        == reaction_pb2.TemperatureConditions.TemperatureControl.ICE_BATH
    ):
      return 0.0
    if (
        ctrl_type
        == reaction_pb2.TemperatureConditions.TemperatureControl.DRY_ICE_BATH
    ):
      return -78.0
  return None


def parse_ord_protobuf(raw_bytes: bytes) -> dict[str, Any] | None:
  """Parses raw ORD reaction protobuf bytes into a standard dictionary."""
  rxn = reaction_pb2.Reaction()
  try:
    rxn.ParseFromString(raw_bytes)
  except (protobuf_message.DecodeError, TypeError) as exc:
    sys.stderr.write(f"[ORD] Failed to parse Reaction protobuf: {exc}\n")
    return None

  raw_rxn_smi = ""
  for ident in rxn.identifiers:
    if ident.type == reaction_pb2.ReactionIdentifier.REACTION_SMILES:
      raw_rxn_smi = ident.value
      break
    if (
        ident.type == reaction_pb2.ReactionIdentifier.REACTION_CXSMILES
        and not raw_rxn_smi
    ):
      raw_rxn_smi = ident.value

  reagents_list: list[str] = []
  solvents_list: list[str] = []
  reactants_list: list[str] = []
  reactant_smiles: list[str] = []
  unspecified_inputs: list[str] = []

  for inp in rxn.inputs.values():
    for comp in inp.components:
      comp_smi = ""
      for ident in comp.identifiers:
        if ident.type == reaction_pb2.CompoundIdentifier.SMILES:
          comp_smi = _canonicalize_valid_smiles(ident.value)
          break
      name = comp_smi or (
          comp.identifiers[0].value if comp.identifiers else "Component"
      )
      role = int(comp.reaction_role)
      if role == reaction_pb2.ReactionRole.REACTANT:
        reactants_list.append(name)
        if comp_smi:
          reactant_smiles.append(comp_smi)
      elif role == reaction_pb2.ReactionRole.UNSPECIFIED and comp_smi:
        unspecified_inputs.append(comp_smi)
      elif role in (
          reaction_pb2.ReactionRole.REAGENT,
          reaction_pb2.ReactionRole.CATALYST,
      ):
        reagents_list.append(name)
      elif role == reaction_pb2.ReactionRole.SOLVENT:
        solvents_list.append(name)

  if not reactants_list and unspecified_inputs:
    reactants_list = list(unspecified_inputs)
    reactant_smiles = list(unspecified_inputs)

  products_list: list[str] = []
  desired_products: list[str] = []
  desired_yield: float | None = None
  fallback_yield: float | None = None
  time_h: float | None = None

  if rxn.outcomes:
    out = rxn.outcomes[0]
    if out.HasField("reaction_time"):
      time_h = _ord_time_to_hours(out.reaction_time)
    for p in out.products:
      p_smi = ""
      for ident in p.identifiers:
        if ident.type == reaction_pb2.CompoundIdentifier.SMILES:
          p_smi = _canonicalize_valid_smiles(ident.value)
          break
      is_desired = bool(p.is_desired_product)
      if p_smi:
        products_list.append(p_smi)
        if is_desired:
          desired_products.append(p_smi)
      prod_yield: float | None = None
      for m in p.measurements:
        if m.HasField("percentage") and m.percentage.HasField("value"):
          if m.type == reaction_pb2.ProductMeasurement.YIELD:
            prod_yield = float(m.percentage.value)
            break
          if (
              m.type == reaction_pb2.ProductMeasurement.UNSPECIFIED
              and prod_yield is None
          ):
            prod_yield = float(m.percentage.value)
      if prod_yield is not None:
        if is_desired and desired_yield is None:
          desired_yield = prod_yield
        if fallback_yield is None:
          fallback_yield = prod_yield

  yield_val = desired_yield if desired_yield is not None else fallback_yield
  chosen_products = desired_products or products_list
  prod_smi = ".".join(chosen_products)

  if (
      reactant_smiles
      and len(reactant_smiles) == len(reactants_list)
      and chosen_products
  ):
    rxn_smi = f"{'.'.join(reactant_smiles)}>>{prod_smi}"
  elif raw_rxn_smi:
    raw_core = raw_rxn_smi.split("|", 1)[0].strip().split()[0]
    if ":" in raw_core or (raw_core.count(">") == 2 and ">>" not in raw_core):
      rxn_smi = normalize_reaction_smiles_for_drfp(raw_rxn_smi) or raw_core
    else:
      rxn_smi = raw_core
    if not reactants_list or not prod_smi:
      parsed_r, parsed_p = _split_reaction_sides(rxn_smi)
      if not reactants_list:
        reactants_list = parsed_r
      if not prod_smi and parsed_p:
        prod_smi = ".".join(parsed_p)
  elif reactant_smiles and chosen_products:
    rxn_smi = f"{'.'.join(reactant_smiles)}>>{prod_smi}"
  else:
    rxn_smi = ""

  temp_c = (
      _ord_temperature_to_celsius(rxn.conditions.temperature)
      if rxn.HasField("conditions") and rxn.conditions.HasField("temperature")
      else None
  )

  cit = ""
  has_patent = False
  if rxn.HasField("provenance"):
    has_patent = bool(rxn.provenance.patent)
    cit = rxn.provenance.patent or rxn.provenance.doi or rxn.reaction_id
  else:
    cit = rxn.reaction_id

  reagents_str = (
      " + ".join(dict.fromkeys(reagents_list))
      if reagents_list
      else "Unspecified"
  )
  solvent_str = (
      " + ".join(dict.fromkeys(solvents_list))
      if solvents_list
      else "Unspecified"
  )
  ord_url = (
      f"https://open-reaction-database.org/id/{rxn.reaction_id}"
      if rxn.reaction_id
      else ""
  )
  notes_str = (
      rxn.notes.procedure_details.strip()
      if rxn.HasField("notes") and rxn.notes.procedure_details
      else ""
  )

  return {
      "reaction_id": rxn.reaction_id,
      "reaction_smiles": rxn_smi,
      "reactants": reactants_list,
      "product": prod_smi,
      "reagents": reagents_str,
      "solvent": solvent_str,
      "temperature_c": temp_c,
      "time_hours": time_h,
      "yield_percent": yield_val,
      "citation": cit or "ORD",
      "example_id": "",
      "notes": notes_str,
      "url": ord_url,
      "source": "ORD/USPTO" if has_patent else "ORD",
  }


def _build_ord_component_clause(
    collection_path: str,
    smiles: str,
    *,
    op: str,
    threshold: float = _ORD_SIMILARITY_THRESHOLD,
) -> dict[str, Any]:
  """Builds an ORD v2 '/api/search' exists clause on inputs or outcomes."""
  inner: dict[str, Any] = {"op": op, "path": "smiles", "smiles": smiles}
  if op == "similarity":
    inner["threshold"] = threshold
  return {"op": "exists", "path": collection_path, "where": inner}


def _merge_ord_candidate_ids(
    and_ids: Sequence[str],
    exact_ids: Sequence[str],
    sim_ids: Sequence[str],
    budget: int = _ORD_BATCH_FETCH_MAX,
) -> list[str]:
  """Combines multi-clause AND hits with a balanced exact/similarity split."""
  merged: list[str] = list(dict.fromkeys(and_ids))
  seen = set(merged)
  remaining = max(0, budget - len(merged))
  if remaining == 0:
    return merged[:budget]

  exact_unique = [rid for rid in exact_ids if rid not in seen]
  sim_unique = [rid for rid in sim_ids if rid not in seen]
  half = remaining // 2
  exact_take = exact_unique[:half]
  sim_take = sim_unique[: remaining - len(exact_take)]
  for rid in exact_take + sim_take:
    if rid not in seen:
      seen.add(rid)
      merged.append(rid)

  if len(merged) < budget:
    for rid in exact_unique[len(exact_take) :] + sim_unique[len(sim_take) :]:
      if len(merged) >= budget:
        break
      if rid not in seen:
        seen.add(rid)
        merged.append(rid)
  return merged[:budget]


def query_ord_api(
    query_smiles: str = "",
    *,
    reactants: Sequence[str] = (),
    product: str = "",
    limit: int = 50,
    timeout_seconds: int = 12,
    errors: list[str] | None = None,
) -> list[dict[str, Any]]:
  """Queries the pre-indexed ORD v2 API (/api/search and /api/reactions)."""
  print("Querying Open Reaction Database (ORD) via pre-indexed v2 API...")

  raw_reactants: list[str] = []
  for r in reactants:
    raw_reactants.extend(_split_smiles_fragments(r))
  raw_products = _split_smiles_fragments(product)
  if query_smiles and (not raw_reactants and not raw_products):
    raw_reactants, raw_products = _split_reaction_sides(query_smiles)

  ranked_reactants = _rank_organic_smiles(raw_reactants)
  ranked_products = _rank_organic_smiles(raw_products)

  if not ranked_reactants and not ranked_products:
    return []

  fetch_cap = min(max(limit * 2, 40), _ORD_BATCH_FETCH_MAX)
  search_timeout = min(float(timeout_seconds), _ORD_SEARCH_TIMEOUT_SECS)

  def _fetch_ord_ids(
      where_clause: dict[str, Any], *, max_ids: int = fetch_cap
  ) -> list[str]:
    payload = {
        "query": {
            "where": where_clause,
            "limit": max_ids,
            "order_by": _ORD_ORDER_BY,
        }
    }
    try:
      resp = _ORD_CLIENT.fetch_json(
          "api/search",
          method="POST",
          json_body=payload,
          timeout=search_timeout,
          max_retries=1,
      )
    except _HTTP_FETCH_ERRORS as exc:
      _record_backend_error(
          f"[ORD] /api/search request failed (where={where_clause}): {exc}",
          errors,
      )
      return []
    if not isinstance(resp, Mapping):
      return []
    return [str(rid) for rid in resp.get("reaction_ids", ()) if rid]

  # Build multi-component AND clause (when both reactants and product exist),
  # single-anchor exact parent hash, and single-anchor Morgan similarity passes.
  anchor_path = "outcomes.products" if ranked_products else "inputs.components"
  anchor_smi = ranked_products[0] if ranked_products else ranked_reactants[0]
  exact_where = _build_ord_component_clause(
      anchor_path, anchor_smi, op="same_parent"
  )

  anchor_mol = Chem.MolFromSmiles(anchor_smi)
  is_small_anchor = (
      anchor_mol is not None
      and anchor_mol.GetNumHeavyAtoms() < _ORD_SMALL_MOL_HEAVY_ATOMS
  )
  sim_threshold = (
      _ORD_SMALL_MOL_SIMILARITY_THRESHOLD
      if is_small_anchor
      else _ORD_SIMILARITY_THRESHOLD
  )
  sim_where = _build_ord_component_clause(
      anchor_path,
      anchor_smi,
      op="similarity",
      threshold=sim_threshold,
  )

  and_where: dict[str, Any] | None = None
  if ranked_products and ranked_reactants:
    and_clauses = [
        _build_ord_component_clause(
            "inputs.components", r_smi, op="same_parent"
        )
        for r_smi in ranked_reactants[:2]
    ] + [
        _build_ord_component_clause(
            "outcomes.products", ranked_products[0], op="same_parent"
        )
    ]
    and_where = {"op": "and", "clauses": and_clauses}

  with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    and_fut = (
        pool.submit(_fetch_ord_ids, and_where, max_ids=_ORD_BATCH_FETCH_MAX)
        if and_where is not None
        else None
    )
    exact_fut = pool.submit(_fetch_ord_ids, exact_where)
    sim_fut = pool.submit(_fetch_ord_ids, sim_where)
    and_ids = and_fut.result() if and_fut is not None else []
    exact_ids = exact_fut.result()
    sim_ids = sim_fut.result()

  # If a 3-clause reactant+reactant+product AND returned empty, try the top
  # reactant + product 2-clause AND before falling back to single-anchor hits.
  if not and_ids and ranked_products and len(ranked_reactants) >= 2:
    pair_where = {
        "op": "and",
        "clauses": [
            _build_ord_component_clause(
                "inputs.components", ranked_reactants[0], op="same_parent"
            ),
            _build_ord_component_clause(
                "outcomes.products", ranked_products[0], op="same_parent"
            ),
        ],
    }
    and_ids = _fetch_ord_ids(pair_where, max_ids=_ORD_BATCH_FETCH_MAX)

  reaction_ids = _merge_ord_candidate_ids(
      and_ids, exact_ids, sim_ids, budget=_ORD_BATCH_FETCH_MAX
  )

  # Fallback to reactant similarity only when product search returned 0 hits.
  if not reaction_ids and ranked_products and ranked_reactants:
    react_sim_where = _build_ord_component_clause(
        "inputs.components",
        ranked_reactants[0],
        op="similarity",
        threshold=_ORD_SIMILARITY_THRESHOLD,
    )
    reaction_ids = list(dict.fromkeys(_fetch_ord_ids(react_sim_where)))

  if not reaction_ids:
    print("ORD v2 search returned 0 candidate reactions.")
    return []

  batch_ids = reaction_ids[:_ORD_BATCH_FETCH_MAX]
  qs = urllib.parse.urlencode([("id", rid) for rid in batch_ids])
  try:
    raw_records = _ORD_CLIENT.fetch_json(
        f"api/reactions?{qs}",
        timeout=float(timeout_seconds),
        max_retries=1,
    )
  except _HTTP_FETCH_ERRORS as exc:
    _record_backend_error(
        f"[ORD] /api/reactions fetch failed for {len(batch_ids)} reaction"
        f" ID(s): {exc}",
        errors,
    )
    return []

  records: list[dict[str, Any]] = []
  if isinstance(raw_records, Sequence):
    for item in raw_records:
      if not isinstance(item, Mapping) or "binpb" not in item:
        continue
      try:
        pb_bytes = base64.b64decode(item["binpb"], validate=True)
      except (ValueError, TypeError) as exc:
        rxn_id = item.get("reaction_id", "<unknown>")
        _record_backend_error(
            f"[ORD] Invalid base64 protobuf payload for reaction '{rxn_id}':"
            f" {exc}",
            errors,
        )
        continue
      rec = parse_ord_protobuf(pb_bytes)
      if rec and rec.get("reaction_smiles"):
        records.append(rec)

  print(f"Retrieved {len(records)} verified ORD reaction record(s).")
  return records


def _tag_local_records(
    records: list[dict[str, Any]],
    precedent_db_path: str,
) -> list[dict[str, Any]]:
  """Sets a default local source label on records loaded from a local file."""
  default_source = f"local:{os.path.basename(precedent_db_path)}"
  for rec in records:
    if isinstance(rec, dict) and not rec.get("source"):
      rec["source"] = default_source
  return records


def load_precedent_database(
    precedent_db_path: str | None = None,
) -> list[dict[str, Any]]:
  """Loads precedent records from a local .parquet, .jsonl, .json, or .csv."""
  if not precedent_db_path:
    return []

  if not os.path.exists(precedent_db_path):
    print(
        f"Warning: Precedent database path '{precedent_db_path}' does not"
        " exist."
    )
    return []

  try:
    if ".parquet" in precedent_db_path:
      df = pd.read_parquet(precedent_db_path)
      if "reaction" in df.columns:
        records = []
        for raw_bytes in df["reaction"]:
          rec = parse_ord_protobuf(raw_bytes)
          if rec and rec.get("reaction_smiles"):
            records.append(rec)
        return _tag_local_records(records, precedent_db_path)
      return _tag_local_records(
          df.astype(object).where(pd.notna(df), None).to_dict(orient="records"),
          precedent_db_path,
      )
    if precedent_db_path.endswith(".jsonl"):
      records = []
      with open(precedent_db_path, "r", encoding="utf-8") as f:
        for line in f:
          stripped = line.strip()
          if stripped:
            records.append(json.loads(stripped))
      return _tag_local_records(records, precedent_db_path)
    if precedent_db_path.endswith(".json"):
      with open(precedent_db_path, "r", encoding="utf-8") as f:
        data = json.load(f)
      if isinstance(data, list):
        return _tag_local_records(data, precedent_db_path)
      if isinstance(data, dict) and "reactions" in data:
        return _tag_local_records(data["reactions"], precedent_db_path)
      return []
    if precedent_db_path.endswith(".csv"):
      df = pd.read_csv(precedent_db_path)
      return _tag_local_records(
          df.astype(object).where(pd.notna(df), None).to_dict(orient="records"),
          precedent_db_path,
      )
    print(
        f"Warning: Unsupported format for '{precedent_db_path}'. Supported:"
        " .parquet, .jsonl, .json, .csv"
    )
    return []
  except (OSError, ValueError, KeyError, pa.ArrowException) as exc:
    sys.stderr.write(
        f"[PrecedentDB] Failed to load precedent database '{precedent_db_path}'"
        f" ({type(exc).__name__}): {exc}\n"
    )
    return []


# ==============================================================================
# Fast Patent Precedent Lookup (Batched PubChem + Google Patents XHR & HTML)
# ==============================================================================


def _record_backend_error(
    msg: str,
    errors: list[str] | None = None,
) -> None:
  """Logs a backend error to stderr and appends it thread-safely to errors."""
  sys.stderr.write(f"{msg}\n")
  if errors is not None:
    with _ERROR_LOCK:
      if msg not in errors:
        errors.append(msg)


def _is_http_404(exc: BaseException) -> bool:
  """Returns True when exc represents an HTTP 404 Not Found response."""
  return isinstance(exc, http_client.HttpError) and exc.status_code == 404


def normalize_google_patent_id(raw_id: str) -> str:
  """Normalizes PCT, legacy WO, and US pre-grant IDs for Google Patents."""
  if not raw_id:
    return ""
  token = raw_id.strip()
  pct_match = _PCT_TOKEN_RE.match(token)
  if pct_match:
    yr, serial = pct_match.groups()
    if len(yr) == 2:
      yr = f"19{yr}" if int(yr) >= 78 else f"20{yr}"
    return f"WO{yr}{serial.zfill(6)}"

  wo_hyphen_match = _WO_HYPHEN_RE.match(token)
  if wo_hyphen_match:
    yr, serial, kind = wo_hyphen_match.groups()
    if len(yr) == 2:
      yr = f"19{yr}" if int(yr) >= 78 else f"20{yr}"
    return f"WO{yr}{serial.zfill(6)}{(kind or '').upper()}"

  cleaned = token.translate(_PATENT_SEP_TRANSLATION).strip().upper()
  wo4_match = _WO_4DIGIT_RE.match(cleaned)
  if wo4_match:
    yr, serial, kind = wo4_match.groups()
    return f"WO{yr}{serial.zfill(6)}{kind or ''}"

  wo2_match = _WO_2DIGIT_RE.match(cleaned)
  if wo2_match:
    yr2, serial, kind = wo2_match.groups()
    full_yr = f"19{yr2}" if int(yr2) >= 78 else f"20{yr2}"
    return f"WO{full_yr}{serial.zfill(6)}{kind or ''}"

  us_pg_match = _US_PREGRANT_KIND_RE.match(cleaned)
  if us_pg_match:
    yr, serial, kind = us_pg_match.groups()
    return f"US{yr}{serial.zfill(7)}{kind}"

  us_stem_match = _US_PREGRANT_STEM_RE.match(cleaned)
  if us_stem_match:
    yr, serial = us_stem_match.groups()
    return f"US{yr}0{serial}"

  return cleaned


def _patent_doc_stem(norm_id: str) -> str:
  """Strips trailing patent kind code (such as B2 or A1) from a patent ID."""
  m = _PATENT_KIND_SUFFIX_RE.match(norm_id)
  return m.group(1) if m else norm_id


def resolve_pubchem_cid_for_smiles(
    smi: str,
    *,
    errors: list[str] | None = None,
) -> int | None:
  """Resolves an organic SMILES string to its primary PubChem CID."""
  clean_smi = _canonicalize_valid_smiles(smi)
  if not _is_searchable_organic_smiles(clean_smi):
    return None
  with _CACHE_LOCK:
    cached = _PUBCHEM_INFO_BY_SMILES.get(clean_smi)
  if cached is not None:
    return cached.get("cid")

  encoded = urllib.parse.quote(clean_smi, safe="")
  url = f"rest/pug/compound/smiles/{encoded}/cids/JSON"
  try:
    data = _PUBCHEM_CLIENT.fetch_json(url, timeout=6.0, max_retries=1)
    ident_list = (
        data.get("IdentifierList") if isinstance(data, Mapping) else None
    )
    cids = ident_list.get("CID") if isinstance(ident_list, Mapping) else None
    if (
        isinstance(cids, Sequence)
        and not isinstance(cids, (str, bytes))
        and cids
        and int(cids[0]) > 0
    ):
      cid = int(cids[0])
      with _CACHE_LOCK:
        _PUBCHEM_INFO_BY_SMILES.setdefault(
            clean_smi,
            {"cid": cid, "names": [], "names_loaded": False},
        )
      return cid
  except (*_HTTP_FETCH_ERRORS, TypeError) as exc:
    if not _is_http_404(exc):
      _record_backend_error(
          f"[PubChem] CID lookup failed for SMILES '{clean_smi}': {exc}",
          errors,
      )
    return None
  return None


def _extract_pubchem_names(row: Mapping[str, Any]) -> list[str]:
  """Extracts non-placeholder Title and IUPACName from a PropertyTable row."""
  names: list[str] = []
  for field in ("Title", "IUPACName"):
    val = str(row.get(field) or "").strip()
    if val and not val.startswith("CID ") and val not in names:
      names.append(val)
  return names


def _batch_resolve_pubchem_compounds(
    smiles_list: Sequence[str],
    *,
    errors: list[str] | None = None,
) -> dict[str, dict[str, Any]]:
  """Resolves CIDs, Titles, and IUPAC names in batched PubChem requests."""
  unique_smiles: list[str] = []
  seen_smiles: set[str] = set()
  smiles_by_inchikey: dict[str, list[str]] = {}
  info_by_smiles: dict[str, dict[str, Any]] = {}

  for raw_smi in smiles_list:
    smi = _canonicalize_valid_smiles(raw_smi)
    if not _is_searchable_organic_smiles(smi) or smi in seen_smiles:
      continue
    seen_smiles.add(smi)
    with _CACHE_LOCK:
      cached = _PUBCHEM_INFO_BY_SMILES.get(smi)
    if cached is not None and cached.get("names_loaded", True):
      info_by_smiles[smi] = cached
      continue
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
      continue
    unique_smiles.append(smi)
    ik = Chem.MolToInchiKey(mol)
    if ik:
      smiles_by_inchikey.setdefault(ik, []).append(smi)

  if smiles_by_inchikey:
    ik_csv = ",".join(smiles_by_inchikey.keys())
    prop_url = (
        f"rest/pug/compound/inchikey/{ik_csv}/property/"
        "InChIKey,Title,IUPACName/JSON"
    )
    try:
      data = _PUBCHEM_CLIENT.fetch_json(prop_url, timeout=8.0, max_retries=1)
      prop_table = (
          data.get("PropertyTable") if isinstance(data, Mapping) else None
      )
      raw_props = (
          prop_table.get("Properties")
          if isinstance(prop_table, Mapping)
          else None
      )
      rows = (
          [r for r in raw_props if isinstance(r, Mapping)]
          if isinstance(raw_props, Sequence)
          else []
      )
      props = sorted(
          rows,
          key=lambda r: int(r.get("CID") or 0)
          if int(r.get("CID") or 0) > 0
          else 10**12,
      )
      for row in props:
        ik = str(row.get("InChIKey") or "")
        cid = int(row.get("CID") or 0)
        if not ik or cid <= 0:
          continue
        entry = {
            "cid": cid,
            "names": _extract_pubchem_names(row),
            "names_loaded": True,
        }
        for matched_smi in smiles_by_inchikey.get(ik, ()):
          if matched_smi not in info_by_smiles:
            info_by_smiles[matched_smi] = entry
            with _CACHE_LOCK:
              _PUBCHEM_INFO_BY_SMILES[matched_smi] = entry
    except (*_HTTP_FETCH_ERRORS, TypeError) as exc:
      if not _is_http_404(exc):
        _record_backend_error(
            "[PubChem] Batch InChIKey lookup failed; falling back to"
            f" per-SMILES resolution: {exc}",
            errors,
        )

  unresolved = [s for s in unique_smiles if s not in info_by_smiles]
  if unresolved:
    smiles_by_cid: dict[int, list[str]] = {}
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=min(4, len(unresolved))
    ) as pool:
      futs = {
          pool.submit(
              resolve_pubchem_cid_for_smiles,
              s,
              errors=errors,
          ): s
          for s in unresolved
      }
      for fut, smi in futs.items():
        cid = fut.result()
        if cid:
          entry = {"cid": cid, "names": [], "names_loaded": False}
          info_by_smiles[smi] = entry
          smiles_by_cid.setdefault(cid, []).append(smi)
    if smiles_by_cid:
      cid_csv = ",".join(str(c) for c in smiles_by_cid)
      cid_prop_url = (
          f"rest/pug/compound/cid/{cid_csv}/property/Title,IUPACName/JSON"
      )
      try:
        cid_data = _PUBCHEM_CLIENT.fetch_json(
            cid_prop_url,
            timeout=6.0,
            max_retries=1,
        )
        prop_table = (
            cid_data.get("PropertyTable")
            if isinstance(cid_data, Mapping)
            else None
        )
        raw_props = (
            prop_table.get("Properties")
            if isinstance(prop_table, Mapping)
            else None
        )
        for row in raw_props if isinstance(raw_props, Sequence) else ():
          if not isinstance(row, Mapping):
            continue
          cid = int(row.get("CID") or 0)
          if cid not in smiles_by_cid:
            continue
          hydrated = {
              "cid": cid,
              "names": _extract_pubchem_names(row),
              "names_loaded": True,
          }
          for matched_smi in smiles_by_cid[cid]:
            info_by_smiles[matched_smi] = hydrated
            with _CACHE_LOCK:
              _PUBCHEM_INFO_BY_SMILES[matched_smi] = hydrated
      except (*_HTTP_FETCH_ERRORS, TypeError) as exc:
        if not _is_http_404(exc):
          _record_backend_error(
              "[PubChem] Batch CID property hydration failed for"
              f" {len(smiles_by_cid)} CID(s): {exc}",
              errors,
          )

  return info_by_smiles


def fetch_pubchem_patent_ids_for_cid(
    cid: int,
    *,
    errors: list[str] | None = None,
) -> list[str]:
  """Fetches PatentID cross-references for a PubChem CID via TXT stream."""
  if not cid or cid <= 0:
    return []
  with _CACHE_LOCK:
    cached = _PATENT_IDS_BY_CID.get(cid)
    cached_err = _PATENT_XREF_ERROR_BY_CID.get(cid)
  if cached is not None:
    if not cached and cached_err:
      _record_backend_error(cached_err, errors)
    return list(cached)

  txt_url = f"rest/pug/compound/cid/{cid}/xrefs/PatentID/TXT"
  try:
    text = _PUBCHEM_CLIENT.fetch_text(
        txt_url,
        timeout=_PUBCHEM_PATENT_XREF_TIMEOUT_SECS,
        max_retries=0,
    )
    pids = tuple(line.strip() for line in text.splitlines() if line.strip())
    with _CACHE_LOCK:
      _PATENT_IDS_BY_CID[cid] = pids
      _PATENT_XREF_ERROR_BY_CID.pop(cid, None)
    return list(pids)
  except _HTTP_FETCH_ERRORS as exc:
    err_msg = (
        f"[PubChem] Patent cross-reference lookup failed for CID {cid}: {exc}"
    )
    is_404 = _is_http_404(exc)
    with _CACHE_LOCK:
      _PATENT_IDS_BY_CID.setdefault(cid, ())
      if is_404:
        _PATENT_XREF_ERROR_BY_CID.pop(cid, None)
      else:
        _PATENT_XREF_ERROR_BY_CID.setdefault(cid, err_msg)
    if not is_404:
      _record_backend_error(err_msg, errors)
    return []


def _patent_publication_year(norm_id: str) -> int:
  """Extracts a 4-digit publication year from a patent ID when present."""
  yr_match = _PATENT_YEAR_RE.match(norm_id)
  if not yr_match:
    return 0
  yr = int(yr_match.group(1))
  return yr if 1900 <= yr <= 2099 else 0


def _empirical_component_idf(patent_set: set[str], pool_size: int) -> float:
  """Computes Robertson-Sparck Jones empirical IDF over the patent pool."""
  if not patent_set or pool_size <= 0:
    return 0.0
  return math.log(1.0 + float(pool_size) / float(len(patent_set)))


def select_candidate_patent_ids(
    reactant_patent_sets: Sequence[set[str]],
    product_patent_sets: Sequence[set[str]],
    max_patents: int = 8,
) -> list[PatentCandidate]:
  """Ranks candidate patents by empirical PubChem IDF co-occurrence."""
  best_id_by_stem: dict[str, str] = {}

  def _to_stem_set(raw_set: set[str]) -> set[str]:
    stems: set[str] = set()
    for raw_pid in raw_set:
      norm_id = normalize_google_patent_id(raw_pid)
      if not norm_id:
        continue
      stem = _patent_doc_stem(norm_id)
      stems.add(stem)
      prev = best_id_by_stem.get(stem)
      if prev is None or (len(norm_id), norm_id) > (len(prev), prev):
        best_id_by_stem[stem] = norm_id
    return stems

  stem_r_sets = [_to_stem_set(s) for s in reactant_patent_sets]
  stem_p_sets = [_to_stem_set(s) for s in product_patent_sets]

  pool_size = sum(len(s) for s in stem_r_sets) + sum(
      len(s) for s in stem_p_sets
  )
  if pool_size <= 0:
    return []

  r_idfs = [_empirical_component_idf(s, pool_size) for s in stem_r_sets]
  p_idfs = [_empirical_component_idf(s, pool_size) for s in stem_p_sets]
  total_idf = sum(r_idfs) + sum(p_idfs)

  # Order non-empty component sets by ascending cardinality (descending IDF).
  # By the pigeonhole principle, every multi-component co-occurrence appears in
  # at least one of the K-1 smallest sets, avoiding iteration over singleton
  # entries of massive bulk-chemical patent sets.
  indexed_sets = sorted(
      [s for s in stem_r_sets + stem_p_sets if s],
      key=len,
  )

  seed_pool: set[str] = set()
  for idx, pat_set in enumerate(indexed_sets):
    if idx < len(indexed_sets) - 1 or len(pat_set) <= max_patents * 16:
      seed_pool |= pat_set
    elif len(seed_pool) < max_patents:
      recent_slice = sorted(
          pat_set,
          key=lambda st: (
              _patent_publication_year(best_id_by_stem.get(st, st)),
              best_id_by_stem.get(st, st),
          ),
          reverse=True,
      )[: max_patents * 8]
      seed_pool |= set(recent_slice)

  scored_candidates: list[PatentCandidate] = []
  for stem in seed_pool:
    r_hits = [i for i, s in enumerate(stem_r_sets) if stem in s]
    p_hits = [j for j, s in enumerate(stem_p_sets) if stem in s]
    w_r = sum(r_idfs[i] for i in r_hits)
    w_p = sum(p_idfs[j] for j in p_hits)
    idf_cov = min(1.0, (w_r + w_p) / total_idf)
    bridge_bonus = min(w_r, w_p) if (r_hits and p_hits) else 0.0
    idf_score = (w_r + w_p) + bridge_bonus
    scored_candidates.append(
        PatentCandidate(
            patent_id=best_id_by_stem.get(stem, stem),
            idf_score=idf_score,
            idf_coverage=idf_cov,
            matched_reactant_indices=r_hits,
            matched_product_indices=p_hits,
        )
    )

  scored_candidates.sort(
      key=lambda c: (
          round(c.idf_score, 6),
          _patent_publication_year(c.patent_id),
          c.patent_id,
      ),
      reverse=True,
  )
  return scored_candidates[:max_patents]


def _pick_searchable_name(info: Mapping[str, Any] | None) -> str:
  """Selects a concise chemical name for Google Patents full-text search."""
  if not info:
    return ""
  for name in info.get("names", ()):
    cleaned = str(name).strip()
    if 3 <= len(cleaned) <= 90:
      return cleaned
  return ""


def _search_google_patents_xhr(
    *,
    reactant_names: Sequence[str],
    product_names: Sequence[str],
    keywords: Sequence[str] | None,
    max_patents: int,
    errors: list[str] | None = None,
) -> list[PatentCandidate]:
  """Queries Google Patents' XHR endpoint for synthesis precedents."""
  top_prod_idx = next(
      (i for i, n in enumerate(product_names) if n),
      None,
  )
  top_react_idx = next(
      (i for i, n in enumerate(reactant_names) if n),
      None,
  )
  kw_clean = [k.strip() for k in keywords or () if k and k.strip()]

  queries: list[tuple[str, list[int], list[int], float]] = []
  if top_prod_idx is not None and top_react_idx is not None:
    p_name = product_names[top_prod_idx]
    r_name = reactant_names[top_react_idx]
    kw_part = f' "{kw_clean[0]}"' if kw_clean else ""
    queries.append((
        f'"{p_name}" "{r_name}"{kw_part} ("mmol" OR "yield")',
        [top_react_idx],
        [top_prod_idx],
        1.0,
    ))
    if kw_clean:
      queries.append((
          f'"{p_name}" "{r_name}" ("mmol" OR "yield")',
          [top_react_idx],
          [top_prod_idx],
          0.95,
      ))
  elif top_prod_idx is not None:
    p_name = product_names[top_prod_idx]
    kw_part = " ".join(f'"{k}"' for k in kw_clean[:2])
    kw_clause = f" {kw_part}" if kw_part else ""
    queries.append((
        f'"{p_name}"{kw_clause} ("mmol" OR "yield")',
        [],
        [top_prod_idx],
        0.75,
    ))
  elif top_react_idx is not None:
    r_name = reactant_names[top_react_idx]
    kw_part = " ".join(f'"{k}"' for k in kw_clean[:2])
    kw_clause = f" {kw_part}" if kw_part else ""
    queries.append((
        f'"{r_name}"{kw_clause} ("mmol" OR "yield")',
        [top_react_idx],
        [],
        0.65,
    ))

  if not queries:
    return []

  candidates: list[PatentCandidate] = []
  seen_stems: set[str] = set()
  num_fetch = min(max(max_patents * 2, 10), 20)

  for q_str, r_indices, p_indices, base_cov in queries:
    if len(candidates) >= max_patents:
      break
    inner_qs = f"q={urllib.parse.quote(q_str)}&num={num_fetch}&language=ENGLISH"
    xhr_path = f"xhr/query?{urllib.parse.urlencode({'url': inner_qs})}"
    try:
      xhr_data = _GOOGLE_PATENTS_CLIENT.fetch_json(
          xhr_path,
          timeout=6.0,
          max_retries=1,
      )
    except _HTTP_FETCH_ERRORS as exc:
      if not _is_http_404(exc):
        _record_backend_error(
            f"[Google Patents] XHR search failed for query '{q_str}': {exc}",
            errors,
        )
      continue
    if not isinstance(xhr_data, Mapping):
      continue
    results_obj = xhr_data.get("results")
    if not isinstance(results_obj, Mapping):
      continue
    clusters = results_obj.get("cluster")
    if (
        not isinstance(clusters, Sequence)
        or not clusters
        or not isinstance(clusters[0], Mapping)
    ):
      continue
    items = clusters[0].get("result")
    if not isinstance(items, Sequence):
      continue
    for rank_idx, item in enumerate(items):
      pat = item.get("patent") if isinstance(item, Mapping) else None
      if not isinstance(pat, Mapping):
        continue
      norm_id = normalize_google_patent_id(
          str(pat.get("publication_number") or "")
      )
      stem = _patent_doc_stem(norm_id)
      if not norm_id or stem in seen_stems:
        continue
      seen_stems.add(stem)

      raw_title = html_lib.unescape(
          _HTML_TAG_STRIP_RE.sub("", str(pat.get("title") or "")).strip()
      )
      raw_snippet = html_lib.unescape(
          _HTML_TAG_STRIP_RE.sub(" ", str(pat.get("snippet") or "")).strip()
      )
      combined_low = f"{raw_title} {raw_snippet}".lower()

      matched_r = list(r_indices)
      for idx, r_name in enumerate(reactant_names):
        if r_name and idx not in matched_r and r_name.lower() in combined_low:
          matched_r.append(idx)
      matched_p = list(p_indices)
      for idx, p_name in enumerate(product_names):
        if p_name and idx not in matched_p and p_name.lower() in combined_low:
          matched_p.append(idx)

      rank_decay = max(0.0, 1.0 - 0.03 * rank_idx)
      candidates.append(
          PatentCandidate(
              patent_id=norm_id,
              idf_score=round(4.0 * base_cov * rank_decay, 4),
              idf_coverage=base_cov,
              matched_reactant_indices=matched_r,
              matched_product_indices=matched_p,
              xhr_title=raw_title,
              xhr_snippet=raw_snippet,
          )
      )

  return candidates[:max_patents]


def _extract_temperature_and_time(
    text: str,
) -> tuple[float | None, float | None]:
  """Extracts reaction temperature (°C) and duration (h) from procedure text."""
  text = text.replace("\u2212", "-")
  temp_c: float | None = None
  tr_match = _TEMP_RANGE_C_RE.search(text)
  if tr_match:
    t_low = float(tr_match.group(1))
    t_high = float(tr_match.group(2))
    if -120.0 <= t_low <= 400.0 and -120.0 <= t_high <= 400.0:
      temp_c = round((t_low + t_high) / 2.0, 2)
  if temp_c is None:
    t_match = _TEMP_C_RE.search(text)
    if t_match:
      t_num = float(t_match.group(1))
      if -120.0 <= t_num <= 400.0:
        temp_c = t_num
    elif _ROOM_TEMP_RE.search(text):
      temp_c = 25.0

  time_h: float | None = None
  tm_match = _TIME_RE.search(text)
  if tm_match:
    range_start, end_val, word_unit, h_unit = tm_match.groups()
    val = (
        (float(range_start) + float(end_val)) / 2.0
        if range_start is not None
        else float(end_val)
    )
    unit = (word_unit or h_unit or "h").lower()
    time_h = round(val / 60.0, 2) if unit.startswith("min") else round(val, 4)
  elif _OVERNIGHT_RE.search(text):
    time_h = 16.0

  return temp_c, time_h


def _score_experimental_paragraph(
    para: str,
    extractor: procedure_extractor.ProcedureExtractor,
    kw_list: Sequence[str],
    *,
    already_cleaned: bool = False,
    has_target_name: bool | None = None,
) -> int:
  """Counts extracted synthesis procedure attributes via ProcedureExtractor."""
  cleaned = (
      para
      if already_cleaned
      else chemical_text_sanitizer.clean_chemical_text(para)
  )
  score = 0
  yield_info = extractor.extract_yield(cleaned)
  if yield_info.get("yield_percent") or yield_info.get("yield_mass"):
    score += 1
  if yield_info.get("yield_moles"):
    score += 1
  if extractor.extract_solvents(cleaned):
    score += 1
  if extractor.extract_reagents(cleaned):
    score += 1
  target_found = (
      has_target_name
      if has_target_name is not None
      else bool(extractor.extract_target_product_name(cleaned))
  )
  if target_found:
    score += 1
  temp_c, time_h = _extract_temperature_and_time(cleaned)
  if temp_c is not None or time_h is not None:
    score += 1
  p_low = cleaned.lower()
  score += sum(kw in p_low for kw in kw_list)
  return score


def _slice_patent_sections(raw_html: str) -> tuple[str, str]:
  """Slices <head> and <section itemprop='description'> before DOM parsing."""
  head_match = _HEAD_END_RE.search(raw_html)
  head_html = (
      raw_html[: head_match.end()]
      if head_match is not None
      else raw_html[:24000]
  )

  desc_match = _DESC_SECTION_START_RE.search(raw_html)
  if desc_match is None:
    return head_html, raw_html

  sec_start = desc_match.start()
  depth = 1
  sec_end = len(raw_html)
  for tag_m in _SECTION_TAG_RE.finditer(raw_html, desc_match.end()):
    if tag_m.group(1):
      depth -= 1
      if depth <= 0:
        sec_end = tag_m.end()
        break
    elif not tag_m.group(0).endswith("/>"):
      depth += 1
  return head_html, raw_html[sec_start:sec_end]


class _HeadMetaParser(html_parser.HTMLParser):
  """Extracts Dublin Core metadata (<meta name='DC.*'>) from <head> HTML."""

  def __init__(self) -> None:
    super().__init__()
    self.title: str = ""
    self.contributors: list[str] = []
    self.pub_date: str = ""

  def handle_starttag(
      self,
      tag: str,
      attrs: list[tuple[str, str | None]],
  ) -> None:
    if tag != "meta":
      return
    attr_by_name = {k.lower(): v or "" for k, v in attrs}
    name = attr_by_name.get("name", "").lower()
    content = " ".join(attr_by_name.get("content", "").split())
    if not content:
      return
    if name == "dc.title" and not self.title:
      self.title = content
    elif name == "dc.contributor":
      self.contributors.append(content)
    elif name == "dc.date" and not self.pub_date:
      self.pub_date = content


@dataclasses.dataclass
class _BlockBuffer:
  """Accumulates text chunks inside an open HTML block tag."""

  tag: str
  depth: int
  parts: list[str] = dataclasses.field(default_factory=list)


@dataclasses.dataclass
class _DivBuffer:
  """Accumulates text chunks inside a fallback <div> element."""

  depth: int
  parts: list[str] = dataclasses.field(default_factory=list)


class _DescriptionBlockParser(html_parser.HTMLParser):
  """Extracts Google Patents description paragraphs without third-party DOM."""

  def __init__(self) -> None:
    super().__init__()
    self.blocks: list[str] = []
    self.fallback_divs: list[str] = []
    self.all_text_parts: list[str] = []
    self._active_blocks: list[_BlockBuffer] = []
    self._active_divs: list[_DivBuffer] = []
    self._skip_depth: int = 0

  def _append_space(self) -> None:
    self.all_text_parts.append(" ")
    for entry in self._active_blocks:
      entry.parts.append(" ")
    for div_entry in self._active_divs:
      div_entry.parts.append(" ")

  def _flush_tag(self, target_tag: str) -> None:
    kept: list[_BlockBuffer] = []
    for entry in self._active_blocks:
      if entry.tag == target_tag:
        joined = " ".join("".join(entry.parts).split())
        if joined:
          self.blocks.append(joined)
      else:
        kept.append(entry)
    self._active_blocks = kept

  def handle_starttag(
      self,
      tag: str,
      attrs: list[tuple[str, str | None]],
  ) -> None:
    if tag in ("script", "style", "noscript"):
      self._skip_depth += 1
      return
    if self._skip_depth > 0:
      return
    if tag in _BLOCK_LEAF_TAGS:
      self._append_space()
      self._flush_tag(tag)
      self._active_blocks.append(_BlockBuffer(tag=tag, depth=1))
    elif tag == "div":
      self._append_space()
      for entry in self._active_blocks:
        if entry.tag == "div":
          entry.depth += 1
      for div_entry in self._active_divs:
        div_entry.depth += 1
      self._active_divs.append(_DivBuffer(depth=1))
      attr_by_name = {k.lower(): v or "" for k, v in attrs}
      classes = set(attr_by_name.get("class", "").split())
      if "description-paragraph" in classes or "description-line" in classes:
        self._active_blocks.append(_BlockBuffer(tag="div", depth=1))
    elif tag in _SPACE_BOUNDARY_TAGS:
      self._append_space()

  def handle_endtag(self, tag: str) -> None:
    if tag in ("script", "style", "noscript"):
      self._skip_depth = max(0, self._skip_depth - 1)
      return
    if self._skip_depth > 0:
      return
    if tag in _BLOCK_LEAF_TAGS:
      self._flush_tag(tag)
      self._append_space()
    elif tag == "table":
      self._flush_tag("tr")
      self._append_space()
    elif tag == "div":
      kept: list[_BlockBuffer] = []
      for entry in self._active_blocks:
        if entry.tag == "div":
          entry.depth -= 1
          if entry.depth <= 0:
            joined = " ".join("".join(entry.parts).split())
            if joined:
              self.blocks.append(joined)
          else:
            kept.append(entry)
        else:
          kept.append(entry)
      self._active_blocks = kept

      kept_divs: list[_DivBuffer] = []
      for div_entry in self._active_divs:
        div_entry.depth -= 1
        if div_entry.depth <= 0:
          joined = " ".join("".join(div_entry.parts).split())
          if joined:
            self.fallback_divs.append(joined)
        else:
          kept_divs.append(div_entry)
      self._active_divs = kept_divs
      self._append_space()
    elif tag in _SPACE_BOUNDARY_TAGS:
      self._append_space()

  def handle_data(self, data: str) -> None:
    if self._skip_depth > 0 or not data:
      return
    self.all_text_parts.append(data)
    for entry in self._active_blocks:
      entry.parts.append(data)
    for div_entry in self._active_divs:
      div_entry.parts.append(data)

  def finish(self) -> None:
    self.close()
    for entry in self._active_blocks:
      joined = " ".join("".join(entry.parts).split())
      if joined:
        self.blocks.append(joined)
    self._active_blocks.clear()
    for div_entry in self._active_divs:
      joined = " ".join("".join(div_entry.parts).split())
      if joined:
        self.fallback_divs.append(joined)
    self._active_divs.clear()


def _compile_synonym_patterns(
    synonym_groups: Sequence[Sequence[str]],
) -> list[tuple[re.Pattern[str], ...]]:
  """Compiles chemical-token-bounded regexes for each component's synonyms."""
  compiled: list[tuple[re.Pattern[str], ...]] = []
  for syns in synonym_groups:
    pats: list[re.Pattern[str]] = []
    seen: set[str] = set()
    for raw_syn in syns:
      syn = chemical_text_sanitizer.clean_chemical_text(raw_syn).strip()
      key = syn.lower()
      if len(syn) < 3 or key in seen:
        continue
      seen.add(key)
      pats.append(
          re.compile(
              rf"(?<![A-Za-z0-9\-]){re.escape(syn)}(?![A-Za-z0-9\-])",
              re.IGNORECASE,
          )
      )
    compiled.append(tuple(pats))
  return compiled


def _is_product_synthesized_in_text(
    text: str,
    patterns: Sequence[re.Pattern[str]],
) -> bool:
  """Returns True if a product synonym appears without reactant-only cues."""
  for pat in patterns:
    for m in pat.finditer(text):
      before = text[max(0, m.start() - 80) : m.start()]
      if _REACTANT_CHARGE_PREFIX_RE.search(before):
        continue
      after = text[m.end() : m.end() + 90]
      clause_before = _CLAUSE_SPLIT_RE.split(before)[-1]
      clause_after = _CLAUSE_SPLIT_RE.split(after)[0]
      if _PRODUCT_CUE_RE.search(clause_before) or _PRODUCT_AFTER_CUE_RE.search(
          clause_after
      ):
        return True
      if not _REACTANT_CHARGE_QTY_RE.search(
          after
      ) and not _REACTANT_ACTION_SUFFIX_RE.search(after):
        return True
  return False


def _extract_uncovered_star_fragments(
    mol: Chem.Mol,
    covered_atoms: AbstractSet[int],
) -> list[str]:
  """Extracts uncovered molecular fragments with [*] capping severed bonds."""
  if not covered_atoms or len(covered_atoms) >= mol.GetNumAtoms():
    return []
  work_mol = Chem.Mol(mol)
  Chem.KekulizeIfPossible(work_mol, clearAromaticFlags=True)
  cut_bonds = [
      b.GetIdx()
      for b in work_mol.GetBonds()
      if (b.GetBeginAtomIdx() in covered_atoms)
      != (b.GetEndAtomIdx() in covered_atoms)
  ]
  if cut_bonds:
    frag_mol = Chem.FragmentOnBonds(
        work_mol,
        cut_bonds,
        addDummies=True,
        dummyLabels=[(0, 0)] * len(cut_bonds),
    )
  else:
    frag_mol = work_mol
  frag_atom_indices: list[tuple[int, ...]] = []
  frags = Chem.GetMolFrags(
      frag_mol,
      asMols=True,
      sanitizeFrags=False,
      fragsMolAtomMapping=frag_atom_indices,
  )
  out: list[str] = []
  n_orig = work_mol.GetNumAtoms()
  for fmol, orig_idxs in zip(frags, frag_atom_indices, strict=True):
    real_idxs = [idx for idx in orig_idxs if idx < n_orig]
    if real_idxs and all(idx not in covered_atoms for idx in real_idxs):
      for atom in fmol.GetAtoms():
        if atom.GetAtomicNum() == 0:
          atom.SetIsotope(0)
      try:
        Chem.SanitizeMol(fmol)
      except _MOL_SANITIZE_ERRORS:
        continue
      smi = Chem.MolToSmiles(fmol, canonical=True)
      if smi:
        out.append(smi)
  return out


def _best_disjoint_mcs_match(
    r_mol: Chem.Mol,
    p_mol: Chem.Mol,
    covered_p_atoms: AbstractSet[int],
) -> tuple[set[int], set[int]]:
  """Finds uncovered MCS atom sets in (r_mol, p_mol) minimizing cut bonds."""

  def _cut_count(mol: Chem.Mol, match_set: AbstractSet[int]) -> int:
    return sum(
        (b.GetBeginAtomIdx() in match_set) != (b.GetEndAtomIdx() in match_set)
        for b in mol.GetBonds()
    )

  keep_bonds = [
      b.GetIdx()
      for b in p_mol.GetBonds()
      if b.GetBeginAtomIdx() not in covered_p_atoms
      and b.GetEndAtomIdx() not in covered_p_atoms
  ]
  if r_mol.GetNumAtoms() == 1 or (covered_p_atoms and not keep_bonds):
    atom_pairs: list[tuple[int, int]] = []
    for r_idx in range(r_mol.GetNumAtoms()):
      r_z = r_mol.GetAtomWithIdx(r_idx).GetAtomicNum()
      for p_idx in range(p_mol.GetNumAtoms()):
        if (
            p_idx not in covered_p_atoms
            and p_mol.GetAtomWithIdx(p_idx).GetAtomicNum() == r_z
        ):
          atom_pairs.append((r_idx, p_idx))
    if not atom_pairs:
      return set(), set()
    best_r_idx, best_p_idx = min(
        atom_pairs,
        key=lambda rp: (
            _cut_count(p_mol, covered_p_atoms | {rp[1]}),
            _cut_count(r_mol, {rp[0]}),
        ),
    )
    return {best_r_idx}, {best_p_idx}

  sub_to_orig: dict[int, int] | None = None
  if covered_p_atoms:
    sub_map: dict[int, int] = {}
    p_target = Chem.PathToSubmol(p_mol, keep_bonds, atomMap=sub_map)
    Chem.FastFindRings(p_target)
    sub_to_orig = {sub_idx: orig_idx for orig_idx, sub_idx in sub_map.items()}
  else:
    p_target = p_mol

  mcs = rdFMCS.FindMCS(
      [r_mol, p_target],
      ringMatchesRingOnly=True,
      completeRingsOnly=False,
      bondCompare=rdFMCS.BondCompare.CompareOrderExact,
      atomCompare=rdFMCS.AtomCompare.CompareElements,
      timeout=_MCS_TIMEOUT_SECS,
  )
  min_expected = min(r_mol.GetNumAtoms(), p_target.GetNumAtoms())
  if mcs.numAtoms < min_expected:
    mcs_any = rdFMCS.FindMCS(
        [r_mol, p_target],
        ringMatchesRingOnly=False,
        completeRingsOnly=False,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        timeout=_MCS_TIMEOUT_SECS,
    )
    if mcs_any.numAtoms > mcs.numAtoms:
      mcs = mcs_any

  if not mcs.smartsString or mcs.numAtoms == 0:
    return set(), set()
  q_mcs = Chem.MolFromSmarts(mcs.smartsString)
  if q_mcs is None:
    return set(), set()

  if sub_to_orig is not None:
    p_candidates = [
        {sub_to_orig[idx] for idx in m}
        for m in p_target.GetSubstructMatches(q_mcs)
    ]
  else:
    p_candidates = [set(m) for m in p_mol.GetSubstructMatches(q_mcs)]
  r_candidates = [set(m) for m in r_mol.GetSubstructMatches(q_mcs)]
  if not p_candidates or not r_candidates:
    return set(), set()

  best_r = min(r_candidates, key=lambda ms: _cut_count(r_mol, ms))
  best_p = min(
      p_candidates, key=lambda ms: _cut_count(p_mol, covered_p_atoms | ms)
  )
  return best_r, best_p


def _deduce_minimal_coreactants(
    attested_reactants: Sequence[str],
    attested_products: Sequence[str],
    *,
    include_leaving_groups: bool = False,
) -> tuple[list[str], list[str], bool]:
  """Deduces [*]-capped missing co-reactant fragments from uncovered product.

  Args:
    attested_reactants: Canonical SMILES of reactants attested in the patent.
    attested_products: Canonical SMILES of products attested in the patent.
    include_leaving_groups: Whether to also return [*]-capped leaving groups
      cleaved from attested reactants.

  Returns:
    A 3-tuple `(coreactants, leaving_groups, matched_product)` where
    `matched_product` is True iff at least one attested reactant mapped onto
    `attested_products`.
  """
  p_mol = Chem.MolFromSmiles(".".join(attested_products))
  if p_mol is None or p_mol.GetNumAtoms() == 0:
    return [], [], False

  parsed_r = [
      m
      for r_smi in attested_reactants
      if (m := Chem.MolFromSmiles(r_smi)) is not None and m.GetNumAtoms() > 0
  ]
  parsed_r.sort(key=lambda m: m.GetNumAtoms(), reverse=True)

  covered_p: set[int] = set()
  leaving_groups: list[str] = []
  for r_mol in parsed_r:
    r_match, p_match = _best_disjoint_mcs_match(r_mol, p_mol, covered_p)
    if not p_match:
      continue
    covered_p.update(p_match)
    if include_leaving_groups and r_match:
      leaving_groups.extend(_extract_uncovered_star_fragments(r_mol, r_match))

  return (
      _extract_uncovered_star_fragments(p_mol, covered_p),
      leaving_groups,
      bool(covered_p),
  )


def _score_partial_drfp_with_rgroups(
    *,
    query_rxn_smi: str,
    query_reactants: Sequence[str],
    query_products: Sequence[str],
    attested_reactants: Sequence[str],
    attested_products: Sequence[str],
) -> float:
  """Scores an incomplete reaction via DRFP after [*]-capped R1 completion."""
  deduced_r, deduced_extra, matched_product = _deduce_minimal_coreactants(
      attested_reactants,
      attested_products,
      include_leaving_groups=len(query_products) > len(attested_products),
  )
  scoring_lhs = list(attested_reactants) + deduced_r
  scoring_rhs = list(attested_products) + deduced_extra
  sim = calculate_reaction_similarity(
      query_rxn_smi,
      list(query_reactants),
      f"{'.'.join(scoring_lhs)}>>{'.'.join(scoring_rhs)}",
      scoring_lhs,
  )
  if not matched_product:
    total_q = max(1, len(query_reactants) + len(query_products))
    sim *= (len(attested_reactants) + len(attested_products)) / total_q
  return sim


def _score_patent_precedent_similarity(
    *,
    query_rxn_smi: str,
    query_reactants: Sequence[str],
    query_products: Sequence[str],
    prec_rxn_smi: str,
    attested_reactants: Sequence[str],
    attested_products: Sequence[str],
    prod_sim: float,
) -> float:
  """Computes DRFP-convergent similarity for complete or partial patent hits."""
  if not query_reactants:
    return round(prod_sim, 3) if attested_products else 0.0

  if (
      attested_reactants
      and len(attested_reactants) >= len(query_reactants)
      and len(attested_products) >= len(query_products)
  ):
    return round(
        calculate_reaction_similarity(
            query_rxn_smi,
            list(query_reactants),
            prec_rxn_smi,
            attested_reactants,
        ),
        3,
    )

  if attested_reactants and attested_products:
    return round(
        _score_partial_drfp_with_rgroups(
            query_rxn_smi=query_rxn_smi,
            query_reactants=query_reactants,
            query_products=query_products,
            attested_reactants=attested_reactants,
            attested_products=attested_products,
        ),
        3,
    )

  if attested_products:
    sim = _PRODUCT_ONLY_PRECEDENT_WEIGHT * prod_sim
    single_r_sims = [
        _score_partial_drfp_with_rgroups(
            query_rxn_smi=query_rxn_smi,
            query_reactants=query_reactants,
            query_products=query_products,
            attested_reactants=[r_smi],
            attested_products=attested_products,
        )
        for r_smi in query_reactants
    ]
    if single_r_sims:
      sim = min(sim, _PRODUCT_ONLY_DRFP_CAP_RATIO * min(single_r_sims))
    return round(sim, 3)

  if attested_reactants:
    r_sim = calculate_reaction_similarity(
        query_rxn_smi,
        list(query_reactants),
        prec_rxn_smi,
        attested_reactants,
    )
    weight = _REACTANT_ONLY_PATENT_WEIGHT if query_products else 1.0
    return round(weight * r_sim, 3)

  return 0.0


@dataclasses.dataclass(frozen=True)
class _ScoredParagraph:
  """Experimental patent paragraph scored during the single-pass scan."""

  raw_para: str
  cleaned_para: str
  proc_score: int
  rank_score: int
  r_hits: list[int]
  p_hits: list[int]
  prod_name: str
  para_idx: int


def _is_resolvable_chemical_product_name(prod_name: str) -> bool:
  """Returns True if prod_name is a chemical name rather than a placeholder."""
  cleaned = prod_name.strip()
  if not cleaned or not _CHEMICAL_NAME_MARKER_RE.search(cleaned):
    return False
  return not _GENERIC_PLACEHOLDER_NAME_RE.match(cleaned)


def _alnum_4grams(text: str) -> set[str]:
  """Extracts lowercase 4-character alphanumeric substrings."""
  alnum = _NON_ALNUM_RE.sub("", text.lower())
  if len(alnum) < 4:
    return set()
  return {alnum[k : k + 4] for k in range(len(alnum) - 3)}


def _extract_ms_mz_values(text: str) -> list[float]:
  """Extracts observed mass-spectrometry m/z values from procedure text."""
  values: list[float] = []
  for m in _MS_HEADER_WINDOW_RE.finditer(text):
    for num_m in _MS_STANDALONE_MZ_RE.finditer(m.group(1)):
      try:
        values.append(float(num_m.group(1)))
      except ValueError:
        continue
  for pat in (_MS_ADDUCT_AFTER_MZ_RE, _MS_ADDUCT_BEFORE_MZ_RE):
    for m in pat.finditer(text):
      try:
        values.append(float(m.group(1)))
      except ValueError:
        continue
  return values


def _query_product_ms_and_elements(
    query_products: Sequence[str],
) -> tuple[list[float], set[str]]:
  """Computes expected ESI/APCI ion m/z values and element symbols."""
  expected_ions: list[float] = []
  elements: set[str] = set()
  for q_smi in query_products:
    for frag in _split_smiles_fragments(q_smi):
      canon = _canonicalize_valid_smiles(frag)
      mol = Chem.MolFromSmiles(canon) if canon else None
      if mol is None:
        continue
      mw = float(Descriptors.ExactMolWt(mol))
      expected_ions.extend((mw, mw + 1.0078, mw - 1.0078, mw + 22.9898))
      elements.update(atom.GetSymbol() for atom in mol.GetAtoms())
  return expected_ions, elements


def parse_google_patent_precedent(
    candidate: PatentCandidate,
    query_rxn_smi: str,
    query_reactants: Sequence[str],
    query_products: Sequence[str],
    keywords: Sequence[str] | None = None,
    extractor: procedure_extractor.ProcedureExtractor | None = None,
    *,
    reactant_synonyms: Sequence[Sequence[str]] = (),
    product_synonyms: Sequence[Sequence[str]] = (),
    errors: list[str] | None = None,
) -> dict[str, Any] | None:
  """Fetches Google Patents HTML and extracts experimental precedent details."""
  norm_id = normalize_google_patent_id(candidate.patent_id)
  if not norm_id:
    return None

  with _CACHE_LOCK:
    raw_html = _HTML_BY_PATENT_ID.get(norm_id)
    cached_err = _HTML_ERROR_BY_PATENT_ID.get(norm_id, "")
  if raw_html is None:
    try:
      raw_html = _GOOGLE_PATENTS_CLIENT.fetch_text(
          f"patent/{norm_id}/en",
          timeout=_GOOGLE_PATENTS_HTML_TIMEOUT_SECS,
          max_retries=1,
      )
    except _HTTP_FETCH_ERRORS as exc:
      err_msg = (
          f"[Google Patents] Failed to fetch specification for '{norm_id}':"
          f" {exc}"
      )
      is_404 = _is_http_404(exc)
      with _CACHE_LOCK:
        _HTML_BY_PATENT_ID.setdefault(norm_id, "")
        if is_404:
          _HTML_ERROR_BY_PATENT_ID.pop(norm_id, None)
        else:
          _HTML_ERROR_BY_PATENT_ID.setdefault(norm_id, err_msg)
      if not is_404:
        _record_backend_error(err_msg, errors)
      return None
    with _CACHE_LOCK:
      _HTML_BY_PATENT_ID[norm_id] = raw_html
      _HTML_ERROR_BY_PATENT_ID.pop(norm_id, None)
  elif not raw_html:
    if cached_err:
      _record_backend_error(cached_err, errors)
    return None
  if not raw_html:
    return None

  proc_extractor = extractor or procedure_extractor.ProcedureExtractor()
  head_html, desc_html = _slice_patent_sections(raw_html)

  head_parser = _HeadMetaParser()
  head_parser.feed(head_html)
  head_parser.close()
  meta_title = head_parser.title or candidate.xhr_title or norm_id
  contribs = head_parser.contributors
  pub_date = head_parser.pub_date

  desc_parser = _DescriptionBlockParser()
  desc_parser.feed(desc_html)
  desc_parser.finish()
  raw_blocks = desc_parser.blocks or desc_parser.fallback_divs

  paragraphs: list[str] = []
  pending_heading = ""
  for block in raw_blocks:
    if len(block) < 40:
      pending_heading = block if _HEADING_LINE_RE.match(block) else ""
      continue
    if pending_heading:
      sep = " " if pending_heading.endswith((".", ":", ";")) else ". "
      merged = f"{pending_heading}{sep}{block}".strip()
    else:
      merged = block
    pending_heading = ""
    if 40 <= len(merged) <= 2400:
      paragraphs.append(merged)

  if not paragraphs:
    raw_text = " ".join("".join(desc_parser.all_text_parts).split())
    paragraphs = [
        raw_text[i : i + 900]
        for i in range(0, min(len(raw_text), 18000), 750)
        if raw_text[i : i + 900].strip()
    ]
  if candidate.xhr_snippet and 40 <= len(candidate.xhr_snippet) <= 2400:
    paragraphs.append(candidate.xhr_snippet)

  r_patterns = _compile_synonym_patterns(reactant_synonyms)
  p_patterns = _compile_synonym_patterns(product_synonyms)
  has_any_syns = any(bool(p) for p in (*r_patterns, *p_patterns))
  kw_list = [k.strip().lower() for k in keywords or () if k and k.strip()]

  best_para = ""
  best_cleaned = ""
  best_proc_score = -1
  best_rank_score = -1
  best_r_hits: list[int] = []
  best_p_hits: list[int] = []
  scored_paragraphs: list[_ScoredParagraph] = []

  for para_idx, para in enumerate(paragraphs):
    if not _EXPERIMENTAL_HINT_RE.search(para) and not any(
        kw in para.lower() for kw in kw_list
    ):
      continue
    cleaned_p = chemical_text_sanitizer.clean_chemical_text(para)
    para_r_hits = [
        i
        for i, pats in enumerate(r_patterns)
        if any(pat.search(cleaned_p) for pat in pats)
    ]
    para_p_hits = [
        j
        for j, pats in enumerate(p_patterns)
        if _is_product_synthesized_in_text(cleaned_p, pats)
    ]
    prod_name = (
        proc_extractor.extract_target_product_name(cleaned_p) or ""
    ).strip()
    proc_score = _score_experimental_paragraph(
        cleaned_p,
        proc_extractor,
        kw_list,
        already_cleaned=True,
        has_target_name=bool(prod_name),
    )
    syn_bonus = (
        8 * len(para_p_hits)
        + 5 * len(para_r_hits)
        + (6 if para_r_hits and para_p_hits else 0)
    )
    rank_score = proc_score + syn_bonus
    if proc_score >= 2:
      scored_paragraphs.append(
          _ScoredParagraph(
              raw_para=para,
              cleaned_para=cleaned_p,
              proc_score=proc_score,
              rank_score=rank_score,
              r_hits=para_r_hits,
              p_hits=para_p_hits,
              prod_name=prod_name,
              para_idx=para_idx,
          )
      )
    if rank_score > best_rank_score:
      best_rank_score = rank_score
      best_proc_score = proc_score
      best_para = para
      best_cleaned = cleaned_p
      best_r_hits = para_r_hits
      best_p_hits = para_p_hits

  if best_rank_score < 0 and paragraphs:
    best_para = paragraphs[0]
    best_cleaned = chemical_text_sanitizer.clean_chemical_text(best_para)
    fallback_prod_name = (
        proc_extractor.extract_target_product_name(best_cleaned) or ""
    ).strip()
    best_proc_score = _score_experimental_paragraph(
        best_cleaned,
        proc_extractor,
        kw_list,
        already_cleaned=True,
        has_target_name=bool(fallback_prod_name),
    )
    best_r_hits = [
        i
        for i, pats in enumerate(r_patterns)
        if any(pat.search(best_cleaned) for pat in pats)
    ]
    best_p_hits = [
        j
        for j, pats in enumerate(p_patterns)
        if _is_product_synthesized_in_text(best_cleaned, pats)
    ]
    if best_proc_score >= 2:
      scored_paragraphs.append(
          _ScoredParagraph(
              raw_para=best_para,
              cleaned_para=best_cleaned,
              proc_score=best_proc_score,
              rank_score=best_proc_score,
              r_hits=best_r_hits,
              p_hits=best_p_hits,
              prod_name=fallback_prod_name,
              para_idx=0,
          )
      )

  if not best_para:
    return None

  def _compute_effective_r_indices(
      r_hits: Sequence[int],
      p_hits: Sequence[int],
      *,
      has_attested_product: bool = False,
  ) -> list[int]:
    out: list[int] = []
    for i in range(len(query_reactants)):
      if i < len(r_patterns) and r_patterns[i]:
        matched_r = i in r_hits
      else:
        matched_r = i in candidate.matched_reactant_indices and (
            bool(r_hits or p_hits) or has_attested_product
        )
      if matched_r:
        out.append(i)
    return out

  if has_any_syns:
    effective_r_indices = _compute_effective_r_indices(best_r_hits, best_p_hits)
    effective_p_indices: list[int] = []
    for j in range(len(query_products)):
      if j < len(p_patterns) and p_patterns[j]:
        matched_p = j in best_p_hits
      else:
        matched_p = j in candidate.matched_product_indices and bool(best_r_hits)
      if matched_p:
        effective_p_indices.append(j)
    if (
        not effective_p_indices
        and query_products
        and len(query_reactants) >= 2
        and len(effective_r_indices) == len(query_reactants)
    ):
      raw_doc_text = f"{meta_title} {' '.join(desc_parser.all_text_parts)}"
      for j in range(len(query_products)):
        pats = p_patterns[j] if j < len(p_patterns) else ()
        if j in candidate.matched_product_indices or any(
            pat.search(raw_doc_text) for pat in pats
        ):
          effective_p_indices.append(j)
  else:
    effective_r_indices = [
        i
        for i in candidate.matched_reactant_indices
        if 0 <= i < len(query_reactants)
    ]
    effective_p_indices = (
        [
            j
            for j in candidate.matched_product_indices
            if 0 <= j < len(query_products)
        ]
        if effective_r_indices
        else []
    )

  resolved_prod_smi = ""
  if query_products and not effective_p_indices and scored_paragraphs:
    require_cand_r_hits = bool(
        has_any_syns and query_reactants and any(r_patterns)
    )
    fallback_gate_ok = (
        not has_any_syns
        or bool(effective_r_indices)
        or (not query_reactants and bool(candidate.matched_product_indices))
    )
    eligible_cands = [
        sp
        for sp in scored_paragraphs
        if _is_resolvable_chemical_product_name(sp.prod_name)
        and (bool(sp.r_hits) if require_cand_r_hits else fallback_gate_ok)
    ]
    if eligible_cands:
      expected_ions, query_elements = _query_product_ms_and_elements(
          query_products
      )
      query_syn_4grams: set[str] = set()
      for syn_group in product_synonyms:
        for syn in syn_group:
          query_syn_4grams.update(_alnum_4grams(syn))

      def _cand_prerank_key(
          sp: _ScoredParagraph,
      ) -> tuple[int, int, int, int, int, int]:
        obs_mzs = _extract_ms_mz_values(sp.cleaned_para)
        ms_ion_match = (
            int(
                any(
                    abs(obs - ion) <= 1.6
                    for obs in obs_mzs
                    for ion in expected_ions
                )
            )
            if expected_ions and obs_mzs
            else 0
        )
        stripped_name = _SALT_WORD_STRIP_RE.sub(" ", sp.prod_name)
        stem_compat = (
            int(
                all(
                    elem in query_elements
                    for pat, elem in _STEM_REQUIRED_ELEMENTS
                    if pat.search(stripped_name)
                )
            )
            if query_elements
            else 1
        )
        syn_4gram_overlap = (
            len(_alnum_4grams(sp.prod_name) & query_syn_4grams)
            if query_syn_4grams
            else 0
        )
        return (
            ms_ion_match,
            stem_compat,
            syn_4gram_overlap,
            len(sp.r_hits),
            sp.rank_score,
            -sp.para_idx,
        )

      ranked_cands = sorted(
          eligible_cands,
          key=_cand_prerank_key,
          reverse=True,
      )
      use_shared_name_cache = (
          type(proc_extractor) is procedure_extractor.ProcedureExtractor
      )
      local_resolved_by_name: dict[str, str] = {}
      parse_calls_used = 0
      win_cand: _ScoredParagraph | None = None
      win_prod_smi = ""
      win_key = (-1.0, -1)

      for cand in ranked_cands:
        norm_name_key = " ".join(cand.prod_name.lower().split())
        if norm_name_key in local_resolved_by_name:
          cand_smi = local_resolved_by_name[norm_name_key]
        else:
          cached_smi: str | None = None
          if use_shared_name_cache:
            with _CACHE_LOCK:
              cached_smi = _RESOLVED_SMILES_BY_PRODUCT_NAME.get(norm_name_key)
          if cached_smi is not None:
            cand_smi = cached_smi
            local_resolved_by_name[norm_name_key] = cand_smi
          else:
            if parse_calls_used >= _MAX_PATENT_PROC_PARSE_CANDIDATES:
              continue
            parse_calls_used += 1
            cand_ex_m = _EXAMPLE_ID_RE.search(cand.raw_para)
            cand_ex_id = cand_ex_m.group(1) if cand_ex_m else ""
            parsed_proc = proc_extractor.parse_procedure(
                cand.raw_para,
                procedure_id=cand_ex_id or norm_id,
            )
            raw_cand_smi = str(
                parsed_proc.get("target_product_smiles") or ""
            ).strip()
            cand_smi = _canonicalize_valid_smiles(raw_cand_smi)
            local_resolved_by_name[norm_name_key] = cand_smi
            if use_shared_name_cache:
              with _CACHE_LOCK:
                _RESOLVED_SMILES_BY_PRODUCT_NAME[norm_name_key] = cand_smi

        if not cand_smi:
          continue
        cand_prod_sim = score_product_similarity(query_products, [cand_smi])
        cand_key = (cand_prod_sim, cand.rank_score)
        if cand_key > win_key:
          win_key = cand_key
          win_cand = cand
          win_prod_smi = cand_smi
        if cand_prod_sim >= 0.999:
          break

      if win_cand is not None and win_key[0] >= _MIN_ANALOG_PRODUCT_MORGAN_SIM:
        best_para = win_cand.raw_para
        best_cleaned = win_cand.cleaned_para
        best_proc_score = win_cand.proc_score
        best_r_hits = win_cand.r_hits
        best_p_hits = win_cand.p_hits
        if has_any_syns:
          effective_r_indices = _compute_effective_r_indices(
              best_r_hits,
              best_p_hits,
              has_attested_product=True,
          )
        resolved_prod_smi = win_prod_smi

  if best_proc_score < 2 and not (effective_r_indices and effective_p_indices):
    return None

  ex_match = _EXAMPLE_ID_RE.search(best_para)
  example_id = ex_match.group(1) if ex_match else ""
  cleaned_para = best_cleaned

  if resolved_prod_smi:
    attested_products = [resolved_prod_smi]
  else:
    attested_products = [
        query_products[j]
        for j in effective_p_indices
        if 0 <= j < len(query_products)
    ]

  if (
      query_products
      and not attested_products
      and (has_any_syns or not effective_r_indices)
  ):
    return None
  if not query_products and not effective_r_indices:
    return None

  yield_info = proc_extractor.extract_yield(cleaned_para)
  solvents_list = proc_extractor.extract_solvents(cleaned_para)
  reagents_list = list(proc_extractor.extract_reagents(cleaned_para))

  yield_val: float | None = None
  raw_yield_str = str(yield_info.get("yield_percent") or "").strip()
  if raw_yield_str:
    is_non_yield_pct = bool(
        re.search(
            rf"\b{re.escape(raw_yield_str)}\s*(?:mol\s*%|wt\s*\.?\s*%|"
            r"%\s*(?:mol|equiv|eq|by\s+(?:weight|mass|volume|mol)|"
            r"of\s+(?!the\s+theor|theor)|Pd|Pt|Rh|Ru|Ni|Cu|Fe|Ir|Au|Ag)\b)",
            cleaned_para,
            re.IGNORECASE,
        )
    )
    if is_non_yield_pct:
      explicit_m = _EXPLICIT_YIELD_PCT_RE.search(cleaned_para)
      raw_yield_str = (
          (explicit_m.group(1) or explicit_m.group(2) or "").strip()
          if explicit_m
          else ""
      )
    if raw_yield_str:
      try:
        y_num = float(raw_yield_str)
        if 0.0 < y_num <= 100.0:
          yield_val = y_num
      except (TypeError, ValueError):
        pass

  temp_c, time_h = _extract_temperature_and_time(cleaned_para)
  solvent_str = (
      " / ".join(solvents_list) if solvents_list else "Patent conditions"
  )

  if not reagents_list:
    cat_match = _CATALYST_RE.search(cleaned_para)
    if cat_match:
      reagents_list.append(cat_match.group(1).strip())
  reagents_str = (
      ", ".join(reagents_list)
      if reagents_list
      else "Patent experimental procedure"
  )

  attested_reactants = [
      query_reactants[i]
      for i in effective_r_indices
      if 0 <= i < len(query_reactants)
  ]
  prec_prod_str = ".".join(attested_products)
  prec_rxn_smi = f"{'.'.join(attested_reactants)}>>{prec_prod_str}"

  prod_sim = (
      score_product_similarity(query_products, attested_products)
      if query_products and attested_products
      else 0.0
  )
  if (
      resolved_prod_smi
      and query_products
      and prod_sim < _MIN_ANALOG_PRODUCT_MORGAN_SIM
  ):
    return None

  similarity_score = _score_patent_precedent_similarity(
      query_rxn_smi=query_rxn_smi,
      query_reactants=query_reactants,
      query_products=query_products,
      prec_rxn_smi=prec_rxn_smi,
      attested_reactants=attested_reactants,
      attested_products=attested_products,
      prod_sim=prod_sim,
  )

  contrib_str = ", ".join(contribs[:4]) if contribs else "Google Patents"
  date_sfx = f" ({pub_date[:4]})" if pub_date else ""
  citation_str = f"{contrib_str} — {norm_id}{date_sfx}"
  patent_url = f"https://patents.google.com/patent/{norm_id}/en"
  notes_snippet = cleaned_para[:320] + (
      "..." if len(cleaned_para) > 320 else ""
  )

  return {
      "reaction_id": norm_id,
      "reaction_smiles": prec_rxn_smi,
      "reactants": attested_reactants,
      "product": prec_prod_str,
      "reagents": reagents_str,
      "solvent": solvent_str,
      "temperature_c": temp_c,
      "time_hours": time_h,
      "yield_percent": yield_val,
      "citation": citation_str,
      "example_id": example_id,
      "similarity_score": similarity_score,
      "title": f"[Patent {norm_id}] {meta_title}",
      "url": patent_url,
      "notes": notes_snippet,
      "source": "Google Patents / PubChem",
  }


def _split_unique_canonical_smiles(items: Sequence[str]) -> list[str]:
  """Splits dot-delimited SMILES entries into unique canonical components."""
  out: list[str] = []
  for item in items:
    for tok in _split_smiles_fragments(item):
      canon = _canonicalize_valid_smiles(tok)
      if canon:
        out.append(canon)
  return list(dict.fromkeys(out))


def query_patent_precedents(
    query_reactants: Sequence[str],
    query_product: str = "",
    query_reaction_smiles: str = "",
    keywords: Sequence[str] | None = None,
    max_patents: int = 8,
    *,
    extractor: procedure_extractor.ProcedureExtractor | None = None,
    errors: list[str] | None = None,
) -> list[dict[str, Any]]:
  """Queries patent precedents concurrently via PubChem and Google Patents."""
  print(
      "Querying PubChem & Google Patents for experimental reaction"
      " precedents..."
  )
  reactants_clean = _split_unique_canonical_smiles(query_reactants)
  products_clean = (
      _split_unique_canonical_smiles([query_product]) if query_product else []
  )

  ranked_reactants = _rank_organic_smiles(reactants_clean)
  ranked_products = _rank_organic_smiles(products_clean)

  if not ranked_reactants and not ranked_products:
    return []

  # 1. Batched PubChem call for CIDs, Titles, and IUPAC names
  info_by_smiles = _batch_resolve_pubchem_compounds(
      ranked_reactants + ranked_products,
      errors=errors,
  )
  reactant_names = [
      _pick_searchable_name(info_by_smiles.get(s)) for s in ranked_reactants
  ]
  product_names = [
      _pick_searchable_name(info_by_smiles.get(s)) for s in ranked_products
  ]

  # Select only the highest-complexity components for PubChem PatentID TXT sets
  # (skipping low-complexity bulk reagents that carry >100k patent IDs).
  xref_reactants = ranked_reactants[:_MAX_PATENT_XREF_REACTANTS]
  xref_products = ranked_products[:1]

  # 2. Run Google Patents XHR search and PubChem PatentID lookups concurrently.
  with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:

    def _submit_cid_xref(
        smi: str,
    ) -> concurrent.futures.Future[list[str]] | None:
      info = info_by_smiles.get(smi)
      if not info or not info.get("cid"):
        return None
      return pool.submit(
          fetch_pubchem_patent_ids_for_cid,
          int(info["cid"]),
          errors=errors,
      )

    xhr_fut = pool.submit(
        _search_google_patents_xhr,
        reactant_names=reactant_names,
        product_names=product_names,
        keywords=keywords,
        max_patents=max_patents,
        errors=errors,
    )
    r_xref_futs = [_submit_cid_xref(s) for s in xref_reactants]
    p_xref_futs = [_submit_cid_xref(s) for s in xref_products]

    xhr_candidates = xhr_fut.result()
    reactant_patent_sets = [
        set(fut.result()) if fut is not None else set() for fut in r_xref_futs
    ]
    product_patent_sets = [
        set(fut.result()) if fut is not None else set() for fut in p_xref_futs
    ]

  idf_candidates = select_candidate_patent_ids(
      reactant_patent_sets=reactant_patent_sets,
      product_patent_sets=product_patent_sets,
      max_patents=max_patents,
  )

  # Merge IDF bridge and Google Patents XHR candidates by document stem.
  merged_by_stem = {
      _patent_doc_stem(cand.patent_id): cand for cand in idf_candidates
  }

  for xhr_cand in xhr_candidates:
    stem = _patent_doc_stem(xhr_cand.patent_id)
    existing = merged_by_stem.get(stem)
    if existing is None:
      merged_by_stem[stem] = xhr_cand
    else:
      preferred_id = (
          xhr_cand.patent_id
          if len(xhr_cand.patent_id) >= len(existing.patent_id)
          else existing.patent_id
      )
      merged_r = sorted(
          set(existing.matched_reactant_indices)
          | set(xhr_cand.matched_reactant_indices)
      )
      merged_p = sorted(
          set(existing.matched_product_indices)
          | set(xhr_cand.matched_product_indices)
      )
      merged_by_stem[stem] = dataclasses.replace(
          existing,
          patent_id=preferred_id,
          idf_score=max(existing.idf_score, xhr_cand.idf_score) + 0.5,
          idf_coverage=max(existing.idf_coverage, xhr_cand.idf_coverage),
          matched_reactant_indices=merged_r,
          matched_product_indices=merged_p,
          xhr_title=existing.xhr_title or xhr_cand.xhr_title,
          xhr_snippet=existing.xhr_snippet or xhr_cand.xhr_snippet,
      )

  candidates = sorted(
      merged_by_stem.values(),
      key=lambda c: (
          bool(c.matched_reactant_indices and c.matched_product_indices),
          round(c.idf_score, 6),
          _patent_publication_year(c.patent_id),
      ),
      reverse=True,
  )[:max_patents]

  if not candidates:
    print("No candidate patent precedents found for query structures.")
    return []

  rxn_smi = query_reaction_smiles or (
      f"{'.'.join(reactants_clean)}>>{'.'.join(products_clean)}"
  )
  shared_extractor = extractor or procedure_extractor.ProcedureExtractor()
  q_reacts_for_parse = ranked_reactants or reactants_clean
  q_prods_for_parse = ranked_products or products_clean
  r_syns_for_parse = [
      tuple(info_by_smiles.get(s, {}).get("names", ()))
      for s in q_reacts_for_parse
  ]
  p_syns_for_parse = [
      tuple(info_by_smiles.get(s, {}).get("names", ()))
      for s in q_prods_for_parse
  ]

  # 3. Fetch & parse candidate Google Patents specifications concurrently
  patent_records: list[dict[str, Any]] = []
  num_workers = min(_MAX_PATENT_HTML_WORKERS, len(candidates))
  with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as pool:
    parse_futs = [
        pool.submit(
            parse_google_patent_precedent,
            candidate=cand,
            query_rxn_smi=rxn_smi,
            query_reactants=q_reacts_for_parse,
            query_products=q_prods_for_parse,
            keywords=keywords,
            extractor=shared_extractor,
            reactant_synonyms=r_syns_for_parse,
            product_synonyms=p_syns_for_parse,
            errors=errors,
        )
        for cand in candidates
    ]
    for fut in parse_futs:
      rec = fut.result()
      if rec is not None:
        patent_records.append(rec)

  print(
      f"Retrieved {len(patent_records)} verified Google Patents reaction"
      " precedent(s)."
  )
  return patent_records


# ==============================================================================
# Vectorized Fingerprint & Similarity Calculation
# ==============================================================================


def compute_drfp_fingerprint(rxn_smiles: str) -> np.ndarray | None:
  """Computes 2048-bit DRFP fingerprint on normalized unmapped SMILES."""
  norm_smi = normalize_reaction_smiles_for_drfp(rxn_smiles)
  if not norm_smi:
    return None
  try:
    fp = DrfpEncoder.encode(norm_smi)
    if fp:
      return np.asarray(fp[0], dtype=bool)
  except _DRFP_ENCODE_ERRORS as exc:
    sys.stderr.write(
        f"[DRFP] Failed to encode reaction SMILES '{norm_smi}': {exc}\n"
    )
    return None
  return None


def _mol_has_defined_stereo(mol: Chem.Mol) -> bool:
  """Returns True if mol specifies tetrahedral chirality or alkene geometry."""
  for atom in mol.GetAtoms():
    if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
      return True
  for bond in mol.GetBonds():
    if bond.GetStereo() > Chem.rdchem.BondStereo.STEREOANY:
      return True
  return False


def compute_morgan_reactants_fp(
    reactants: Sequence[str],
    *,
    include_chirality: bool = False,
) -> Any:
  """Computes combined 2048-bit Morgan fingerprint across all reactants."""
  gen = _CHIRAL_MORGAN_GEN if include_chirality else _MORGAN_GEN
  bits = None
  for smi in reactants:
    clean_smi = canonicalize_smiles(smi)
    mol = Chem.MolFromSmiles(clean_smi) if clean_smi else None
    if mol is not None:
      fp = gen.GetFingerprint(mol)
      bits = fp if bits is None else (bits | fp)
  return bits


@dataclasses.dataclass(frozen=True)
class _StereoQuery:
  """Precomputed Morgan and constitutional stereo state for a SMILES group."""

  fp: Any
  has_stereo: bool
  flat_canon: str
  mol: Chem.Mol
  isomeric_parts: tuple[str, ...]


def _prepare_stereo_query(smiles_list: Sequence[str]) -> _StereoQuery | None:
  """Builds a query-conditioned stereochemical descriptor from SMILES."""
  unique_parts: list[str] = []
  seen: set[str] = set()
  mols: list[Chem.Mol] = []
  for smi in smiles_list:
    for frag in _split_smiles_fragments(smi):
      canon = _canonicalize_valid_smiles(frag)
      if not canon or canon in seen:
        continue
      mol = Chem.MolFromSmiles(canon)
      if mol is not None:
        seen.add(canon)
        unique_parts.append(canon)
        mols.append(mol)
  if not mols:
    return None
  has_stereo = any(_mol_has_defined_stereo(m) for m in mols)
  combined_mol = Chem.MolFromSmiles(".".join(unique_parts))
  if combined_mol is None:
    return None
  fp = compute_morgan_reactants_fp(
      unique_parts,
      include_chirality=has_stereo,
  )
  if fp is None:
    return None
  return _StereoQuery(
      fp=fp,
      has_stereo=has_stereo,
      flat_canon=Chem.MolToSmiles(combined_mol, isomericSmiles=False),
      mol=combined_mol,
      isomeric_parts=tuple(unique_parts),
  )


def _score_against_stereo_query(
    query: _StereoQuery,
    cand_smiles: Sequence[str],
) -> float:
  """Scores candidate SMILES against a stereo-conditioned query."""
  c_unique: list[str] = []
  seen: set[str] = set()
  for smi in cand_smiles:
    for frag in _split_smiles_fragments(smi):
      canon = _canonicalize_valid_smiles(frag)
      if canon and canon not in seen:
        seen.add(canon)
        c_unique.append(canon)
  if not c_unique:
    return 0.0

  c_fp = compute_morgan_reactants_fp(
      c_unique,
      include_chirality=query.has_stereo,
  )
  if c_fp is None:
    return 0.0
  base_sim = float(DataStructs.TanimotoSimilarity(query.fp, c_fp))
  if not query.has_stereo:
    return base_sim

  c_mol = Chem.MolFromSmiles(".".join(c_unique))
  if (
      c_mol is None
      or c_mol.GetNumAtoms() != query.mol.GetNumAtoms()
      or Chem.MolToSmiles(c_mol, isomericSmiles=False) != query.flat_canon
  ):
    return base_sim

  strict_params = Chem.SubstructMatchParameters()
  strict_params.useChirality = True
  strict_params.specifiedStereoQueryMatchesUnspecified = False
  if c_mol.HasSubstructMatch(query.mol, strict_params):
    return 1.0

  relaxed_params = Chem.SubstructMatchParameters()
  relaxed_params.useChirality = True
  relaxed_params.specifiedStereoQueryMatchesUnspecified = True
  if c_mol.HasSubstructMatch(query.mol, relaxed_params):
    return _UNSPECIFIED_STEREO_PRODUCT_SIM

  return min(base_sim, _CONFLICTING_STEREO_PRODUCT_SIM_CAP)


def _score_candidate_parts_against_stereo_query(
    query: _StereoQuery,
    cand_smiles: Sequence[str],
) -> float:
  """Scores candidate SMILES as a whole and per component for single targets."""
  c_parts: list[str] = []
  for smi in cand_smiles:
    c_parts.extend(_split_smiles_fragments(smi))
  if not c_parts:
    return 0.0

  best_sim = _score_against_stereo_query(query, c_parts)
  if len(query.isomeric_parts) == 1 and len(c_parts) > 1:
    part_sims: list[float] = []
    has_exact_iso = False
    has_conflict_iso = False
    for part in c_parts:
      p_sim = _score_against_stereo_query(query, [part])
      part_sims.append(p_sim)
      if query.has_stereo:
        p_canon = _canonicalize_valid_smiles(part)
        p_mol = Chem.MolFromSmiles(p_canon) if p_canon else None
        if (
            p_mol is not None
            and p_mol.GetNumAtoms() == query.mol.GetNumAtoms()
            and Chem.MolToSmiles(p_mol, isomericSmiles=False)
            == query.flat_canon
        ):
          if p_sim >= 0.999:
            has_exact_iso = True
          elif p_sim < _UNSPECIFIED_STEREO_PRODUCT_SIM:
            has_conflict_iso = True
    if part_sims:
      max_part_sim = max(part_sims)
      if has_exact_iso and has_conflict_iso:
        max_part_sim = min(max_part_sim, _UNSPECIFIED_STEREO_PRODUCT_SIM)
      if max_part_sim > best_sim:
        best_sim = max_part_sim
  return best_sim


def score_product_similarity(
    query_products: Sequence[str],
    cand_products: Sequence[str],
) -> float:
  """Scores product similarity with query-conditioned stereochemistry."""
  q_query = _prepare_stereo_query(query_products)
  if q_query is None:
    return 0.0
  return _score_candidate_parts_against_stereo_query(q_query, cand_products)


def calculate_reaction_similarity(
    query_rxn_smiles: str,
    query_reactants: Sequence[str],
    prec_rxn_smiles: str,
    prec_reactants: Sequence[str],
    *,
    query_product: str = "",
    prec_product: str = "",
) -> float:
  """Calculates reaction similarity via DRFP or reactant/product Morgan FPs."""
  sims = calculate_batch_reaction_similarities(
      query_rxn_smiles=query_rxn_smiles,
      query_reactants=query_reactants,
      candidate_rxn_smiles=[prec_rxn_smiles],
      candidate_reactants=[prec_reactants],
      query_product=query_product,
      candidate_products=[prec_product] if prec_product else (),
  )
  return sims[0] if sims else 0.0


def calculate_batch_reaction_similarities(
    query_rxn_smiles: str,
    query_reactants: Sequence[str],
    candidate_rxn_smiles: Sequence[str],
    candidate_reactants: Sequence[Sequence[str]],
    *,
    query_product: str = "",
    candidate_products: Sequence[str] = (),
) -> list[float]:
  """Scores candidate reactions in one batch DRFP pass with Morgan fallback."""
  n_cands = len(candidate_rxn_smiles)
  if n_cands == 0:
    return []

  scores = [0.0] * n_cands
  norm_query = normalize_reaction_smiles_for_drfp(query_rxn_smiles)

  drfp_indices: list[int] = []
  drfp_smiles: list[str] = []
  if norm_query:
    for idx, cand_smi in enumerate(candidate_rxn_smiles):
      norm_cand = normalize_reaction_smiles_for_drfp(cand_smi)
      if norm_cand:
        drfp_indices.append(idx)
        drfp_smiles.append(norm_cand)

  drfp_succeeded: set[int] = set()
  if norm_query and drfp_smiles:
    try:
      encoded = DrfpEncoder.encode([norm_query] + drfp_smiles)
      fps = np.asarray(encoded, dtype=bool)
      q_fp = fps[0]
      cand_fps = fps[1:]
      inter = np.logical_and(cand_fps, q_fp).sum(axis=1)
      union = np.logical_or(cand_fps, q_fp).sum(axis=1)
      batch_sims = np.where(union > 0, inter / union, 0.0)
      for out_idx, sim_val in zip(drfp_indices, batch_sims, strict=True):
        scores[out_idx] = float(sim_val)
        drfp_succeeded.add(out_idx)
    except _DRFP_ENCODE_ERRORS as exc:
      sys.stderr.write(
          f"[DRFP] Batch encoding failed for {len(drfp_smiles)} candidate(s);"
          f" falling back to reactant Morgan similarity: {exc}\n"
      )

  fallback_indices = [i for i in range(n_cands) if i not in drfp_succeeded]
  if fallback_indices:
    q_lhs, q_rhs = _split_reaction_sides(query_rxn_smiles)
    eff_q_reactants = list(query_reactants) if query_reactants else q_lhs
    q_react_query = _prepare_stereo_query(eff_q_reactants)
    if q_react_query is not None:
      for i in fallback_indices:
        scores[i] = _score_against_stereo_query(
            q_react_query,
            candidate_reactants[i],
        )
    elif not eff_q_reactants:
      q_prod_parts = _split_smiles_fragments(query_product) or q_rhs
      q_prod_query = _prepare_stereo_query(q_prod_parts)
      if q_prod_query is not None:
        for i in fallback_indices:
          raw_cand_prod = (
              candidate_products[i] if i < len(candidate_products) else ""
          )
          c_parts = (
              _split_smiles_fragments(raw_cand_prod)
              or _split_reaction_sides(candidate_rxn_smiles[i])[1]
          )
          if not c_parts:
            continue
          scores[i] = _score_candidate_parts_against_stereo_query(
              q_prod_query,
              c_parts,
          )

  return scores


# ==============================================================================
# Precedent Search, Condition Clustering & Batch Multi-Step Execution
# ==============================================================================


def _to_optional_float(value: Any) -> float | None:
  """Converts finite numeric values to float, or returns None for NaN/empty."""
  if value is None or isinstance(value, bool):
    return None
  try:
    num = float(value)
  except (TypeError, ValueError):
    return None
  return num if math.isfinite(num) else None


def _clean_text_field(value: Any, default: str = "") -> str:
  """Normalizes optional string fields while rejecting NaN float values."""
  if value is None:
    return default
  if isinstance(value, float) and math.isnan(value):
    return default
  text = str(value).strip()
  return text if text else default


def search_precedents(
    query_reactants: Sequence[str],
    query_product: str | None = None,
    query_reaction_smiles: str | None = None,
    precedent_db_path: str | None = None,
    use_api: str | None = None,
    min_similarity: float = 0.35,
    limit: int = 50,
    errors: list[str] | None = None,
    include_patents: bool = False,
    keywords: Sequence[str] | None = None,
    max_patents: int = 8,
) -> list[PrecedentRecord]:
  """Searches precedent DB or live ORD v2 and Patent APIs concurrently."""
  parsed_r, parsed_p = _split_reaction_sides(query_reaction_smiles or "")
  raw_r_tokens: list[str] = []
  for r in query_reactants:
    raw_r_tokens.extend(_split_smiles_fragments(r))
  eff_r_tokens = raw_r_tokens if raw_r_tokens else parsed_r
  reactants = [
      c for tok in eff_r_tokens if (c := _canonicalize_valid_smiles(tok))
  ]
  if eff_r_tokens and not reactants:
    return []

  raw_p_tokens = _split_smiles_fragments(query_product or "")
  eff_p_tokens = raw_p_tokens if raw_p_tokens else parsed_p
  prod_parts = [
      c for tok in eff_p_tokens if (c := _canonicalize_valid_smiles(tok))
  ]
  if eff_p_tokens and not prod_parts and not reactants:
    return []
  product = ".".join(prod_parts)
  if not reactants and not product:
    return []

  canon_parsed_r = [c for r in parsed_r if (c := _canonicalize_valid_smiles(r))]
  canon_parsed_p = ".".join(
      c for p in parsed_p if (c := _canonicalize_valid_smiles(p))
  )
  if (
      query_reaction_smiles
      and sorted(reactants) == sorted(canon_parsed_r)
      and product == canon_parsed_p
  ):
    query_rxn_smi = query_reaction_smiles
  else:
    query_rxn_smi = f"{'.'.join(sorted(reactants))}>>{product}"

  records: list[dict[str, Any]] = []
  api_mode = (use_api or "").lower()
  want_ord = api_mode in ("ord", "all")
  want_patents = api_mode in ("patents", "all") or include_patents

  # Run ORD v2 and Patent lookup concurrently when both sources are requested.
  if want_ord and want_patents:
    ord_errors: list[str] | None = [] if errors is not None else None
    pat_errors: list[str] | None = [] if errors is not None else None
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
      ord_fut = pool.submit(
          query_ord_api,
          query_rxn_smi,
          reactants=reactants,
          product=product,
          limit=limit,
          errors=ord_errors,
      )
      pat_fut = pool.submit(
          query_patent_precedents,
          query_reactants=reactants,
          query_product=product,
          query_reaction_smiles=query_rxn_smi,
          keywords=keywords,
          max_patents=max_patents,
          errors=pat_errors,
      )
      records.extend(ord_fut.result())
      records.extend(pat_fut.result())
    if errors is not None:
      if ord_errors:
        errors.extend(ord_errors)
      if pat_errors:
        errors.extend(pat_errors)
  elif want_ord:
    records.extend(
        query_ord_api(
            query_rxn_smi,
            reactants=reactants,
            product=product,
            limit=limit,
            errors=errors,
        )
    )
  elif want_patents:
    records.extend(
        query_patent_precedents(
            query_reactants=reactants,
            query_product=product,
            query_reaction_smiles=query_rxn_smi,
            keywords=keywords,
            max_patents=max_patents,
            errors=errors,
        )
    )

  if precedent_db_path:
    records.extend(load_precedent_database(precedent_db_path))

  if not records:
    return []

  # Prepare normalized candidates and batch-compute DRFP/Morgan similarities.
  prepared_rxn_smiles: list[str] = []
  prepared_reactants: list[list[str]] = []
  prepared_products: list[str] = []
  unscored_indices: list[int] = []
  precomputed_sims: dict[int, float] = {}

  for idx, rec in enumerate(records):
    prec_rxn_smi = _clean_text_field(
        rec.get("reaction_smiles") or rec.get("rxn_smiles")
    )
    raw_prec_reactants = rec.get("reactants", [])
    if isinstance(raw_prec_reactants, str):
      prec_reactants = _split_smiles_fragments(raw_prec_reactants)
    elif isinstance(raw_prec_reactants, Sequence):
      prec_reactants = [
          _clean_text_field(s)
          for s in raw_prec_reactants
          if _clean_text_field(s)
      ]
    else:
      prec_reactants = []

    prec_product = _clean_text_field(rec.get("product"))
    if not prec_rxn_smi:
      prec_rxn_smi = f"{'.'.join(sorted(prec_reactants))}>>{prec_product}"

    prepared_rxn_smiles.append(prec_rxn_smi)
    prepared_reactants.append(prec_reactants)
    prepared_products.append(prec_product)

    existing_sim = _to_optional_float(rec.get("similarity_score"))
    if existing_sim is not None:
      precomputed_sims[idx] = existing_sim
    else:
      unscored_indices.append(idx)

  if unscored_indices:
    batch_sims = calculate_batch_reaction_similarities(
        query_rxn_smiles=query_rxn_smi,
        query_reactants=reactants,
        candidate_rxn_smiles=[prepared_rxn_smiles[i] for i in unscored_indices],
        candidate_reactants=[prepared_reactants[i] for i in unscored_indices],
        query_product=product,
        candidate_products=[prepared_products[i] for i in unscored_indices],
    )
    for rec_idx, sim_val in zip(unscored_indices, batch_sims, strict=True):
      precomputed_sims[rec_idx] = sim_val

  matched: list[PrecedentRecord] = []
  for idx, rec in enumerate(records):
    sim = precomputed_sims.get(idx, 0.0)
    if sim < min_similarity:
      continue

    base_str = _clean_text_field(rec.get("base"))
    cat_str = _clean_text_field(rec.get("catalyst"))
    if base_str and cat_str:
      reagents = f"{base_str} + {cat_str}"
    else:
      reagents = (
          _clean_text_field(rec.get("reagents"))
          or base_str
          or cat_str
          or "Unspecified"
      )

    temp_c = _to_optional_float(rec.get("temperature_c"))
    raw_time = (
        rec.get("time_hours")
        if rec.get("time_hours") is not None
        else rec.get("time_h")
    )
    time_h = _to_optional_float(raw_time)
    raw_yld = (
        rec.get("yield_percent")
        if rec.get("yield_percent") is not None
        else rec.get("yield")
    )
    yld = _to_optional_float(raw_yld)

    cit = (
        _clean_text_field(rec.get("citation"))
        or _clean_text_field(rec.get("patent_or_ord_id"))
        or _clean_text_field(rec.get("patent_id"))
        or _clean_text_field(rec.get("ord_id"))
        or "Unknown"
    )
    ex_id = _clean_text_field(rec.get("example_id")) or _clean_text_field(
        rec.get("procedure_id")
    )
    rxn_id = _clean_text_field(
        rec.get("reaction_id")
    ) or f"{cit}_{ex_id}".strip("_")

    matched.append(
        PrecedentRecord(
            reaction_id=rxn_id,
            reaction_smiles=prepared_rxn_smiles[idx],
            reactants=prepared_reactants[idx],
            product=prepared_products[idx],
            reagents=reagents,
            solvent=_clean_text_field(rec.get("solvent"), "Unspecified"),
            temperature_c=temp_c,
            time_hours=time_h,
            yield_percent=yld,
            citation=cit,
            example_id=ex_id,
            similarity_score=round(sim, 3),
            title=_clean_text_field(rec.get("title")),
            url=_clean_text_field(rec.get("url")),
            notes=_clean_text_field(rec.get("notes")),
            source=_clean_text_field(rec.get("source"), ""),
        )
    )

  def _rank_key(p: PrecedentRecord) -> tuple[float, int, float, int]:
    yld = p.yield_percent
    yld_tier = 1 if yld is None else (2 if yld > 0.0 else 0)
    completeness = (
        int(
            bool(p.reagents)
            and p.reagents not in (
                "Unspecified",
                "Patent experimental procedure",
            )
        )
        + int(
            bool(p.solvent)
            and p.solvent not in ("Unspecified", "Patent conditions")
        )
        + int(p.temperature_c is not None)
    )
    return (p.similarity_score, yld_tier, yld or 0.0, completeness)

  matched.sort(key=_rank_key, reverse=True)
  return matched[:limit]


def _temperature_bracket(temp_c: float | None) -> str:
  """Classifies reaction temperature (°C) into a standard condition bracket."""
  if temp_c is None:
    return "Ambient/Unspecified"
  if temp_c <= 30.0:
    return "RT (<=30°C)"
  if temp_c <= 70.0:
    return "Warm (31-70°C)"
  return "Heated (>70°C)"


def cluster_reaction_conditions(
    precedents: Sequence[PrecedentRecord],
) -> list[ConditionCluster]:
  """Clusters precedents into condition archetypes by reagents/solvent/temp."""
  if not precedents:
    return []

  buckets: dict[str, list[PrecedentRecord]] = {}
  for p in precedents:
    t_str = _temperature_bracket(p.temperature_c)
    key = f"{p.reagents} in {p.solvent} at {t_str}"
    buckets.setdefault(key, []).append(p)

  clusters: list[ConditionCluster] = []
  for cluster_name, p_list in buckets.items():
    valid_yields = [
        p.yield_percent for p in p_list if p.yield_percent is not None
    ]
    valid_times = [p.time_hours for p in p_list if p.time_hours is not None]
    valid_temps = [
        p.temperature_c for p in p_list if p.temperature_c is not None
    ]

    med_yield = float(np.median(valid_yields)) if valid_yields else None
    yield_rng = (
        (float(np.min(valid_yields)), float(np.max(valid_yields)))
        if valid_yields
        else None
    )

    med_time = float(np.median(valid_times)) if valid_times else None
    time_rng = (
        (float(np.min(valid_times)), float(np.max(valid_times)))
        if valid_times
        else None
    )

    temp_rng = (
        (float(np.min(valid_temps)), float(np.max(valid_temps)))
        if valid_temps
        else None
    )

    citations = list(
        dict.fromkeys(
            f"{p.citation} ({p.example_id})" if p.example_id else p.citation
            for p in p_list
        )
    )

    first_p = p_list[0]
    clusters.append(
        ConditionCluster(
            cluster_name=cluster_name,
            reagents=first_p.reagents,
            solvent=first_p.solvent,
            temp_bracket=_temperature_bracket(first_p.temperature_c),
            precedent_count=len(p_list),
            median_yield=med_yield,
            yield_range=yield_rng,
            median_time_hours=med_time,
            time_range_hours=time_rng,
            temperature_range_c=temp_rng,
            citations=citations,
        )
    )

  clusters.sort(
      key=lambda c: (
          c.precedent_count,
          c.median_yield if c.median_yield is not None else 0.0,
      ),
      reverse=True,
  )
  return clusters


def lookup_reaction_precedents(
    reactants: Sequence[str] = (),
    product: str | None = None,
    reaction_smiles: str | None = None,
    precedent_db_path: str | None = None,
    use_api: str | None = None,
    min_similarity: float = 0.35,
    limit: int = 50,
    include_patents: bool = False,
    keywords: Sequence[str] | None = None,
    max_patents: int = 8,
) -> dict[str, Any]:
  """Retrieves reaction precedents and clusters condition archetypes."""
  parsed_r, parsed_p = _split_reaction_sides(reaction_smiles or "")
  resolved_reactants = list(reactants) if reactants else parsed_r
  prod_frags = _split_smiles_fragments(product or "")
  if prod_frags:
    resolved_product: str | None = ".".join(prod_frags)
  elif parsed_p:
    resolved_product = ".".join(parsed_p)
  else:
    resolved_product = None

  errors: list[str] = []
  precedents = search_precedents(
      query_reactants=resolved_reactants,
      query_product=resolved_product,
      query_reaction_smiles=reaction_smiles,
      precedent_db_path=precedent_db_path,
      use_api=use_api,
      min_similarity=min_similarity,
      limit=limit,
      errors=errors,
      include_patents=include_patents,
      keywords=keywords,
      max_patents=max_patents,
  )

  clusters = cluster_reaction_conditions(precedents)
  if reaction_smiles and not reactants and not prod_frags:
    query_smi = reaction_smiles
  else:
    query_smi = (
        f"{'.'.join(sorted(resolved_reactants))}>>{resolved_product or ''}"
    )

  return {
      "query": {
          "reaction_smiles": query_smi,
          "reactants": resolved_reactants,
          "product": resolved_product,
      },
      "precedent_count": len(precedents),
      "clusters": [c.to_dict() for c in clusters],
      "precedents": [p.to_dict() for p in precedents],
      "errors": errors,
  }


def _resolve_cli_source_flags(
    api: str | None,
    include_patents_flag: bool | None,
) -> tuple[str | None, bool]:
  """Resolves CLI source flags into (use_api, include_patents)."""
  if include_patents_flag is not None:
    if not include_patents_flag:
      if api == "all":
        return "ord", False
      if api == "patents":
        return None, False
      return api, False
    return api, True
  return api, api in ("patents", "all")


class _BatchValidationError(ValueError):
  """Raised when a batch precedent lookup payload or step schema is invalid."""


@dataclasses.dataclass(frozen=True, kw_only=True)
class _PreparedBatchStep:
  """Validated single-step parameters resolved against batch defaults."""

  reactants: list[str]
  product: str | None
  reaction_smiles: str | None
  precedent_db_path: str | None
  use_api: str | None
  min_similarity: float
  limit: int
  include_patents: bool
  keywords: list[str] | None
  max_patents: int
  metadata: dict[str, Any]


def _extract_batch_steps(raw_batch: Any) -> list[Any]:
  """Validates top-level --batch-json structure and returns its step list."""
  if isinstance(raw_batch, Mapping):
    if "steps" not in raw_batch:
      raise _BatchValidationError(
          "Top-level JSON object must contain a 'steps' list."
      )
    raw_steps = raw_batch["steps"]
  else:
    raw_steps = raw_batch

  if not isinstance(raw_steps, Sequence) or isinstance(
      raw_steps, (str, bytes, bytearray)
  ):
    raise _BatchValidationError(
        "Expected a non-empty JSON list of step objects or {'steps': [...]}."
    )
  if not raw_steps:
    raise _BatchValidationError("Step list must contain at least one step.")
  return list(raw_steps)


def _validate_query_smiles_inputs(
    reactants: Sequence[str],
    product: str | None,
    reaction_smiles: str | None,
    *,
    prefix: str = "",
) -> None:
  """Validates that CLI/batch query inputs contain parseable RDKit SMILES."""
  if not reactants and product is None and reaction_smiles is None:
    raise _BatchValidationError(
        f"{prefix}Must provide non-empty 'reaction_smiles', 'reactants', or"
        " 'product'."
    )
  for r in reactants:
    r_frags = _split_smiles_fragments(r)
    if not r_frags or any(
        not _canonicalize_valid_smiles(tok) for tok in r_frags
    ):
      raise _BatchValidationError(f"{prefix}Invalid reactant SMILES: {r!r}.")
  if product is not None:
    p_frags = _split_smiles_fragments(product)
    if not p_frags or any(
        not _canonicalize_valid_smiles(tok) for tok in p_frags
    ):
      raise _BatchValidationError(
          f"{prefix}Invalid product SMILES: {product!r}."
      )
  if reaction_smiles is not None:
    parsed_r, parsed_p = _split_reaction_sides(reaction_smiles)
    if not parsed_r and not parsed_p and not reactants and product is None:
      raise _BatchValidationError(
          f"{prefix}'reaction_smiles' must contain at least one reactant or"
          " product SMILES."
      )
    for tok in parsed_r + parsed_p:
      if not _canonicalize_valid_smiles(tok):
        raise _BatchValidationError(
            f"{prefix}Invalid SMILES token {tok!r} in 'reaction_smiles':"
            f" {reaction_smiles!r}."
        )


def _prepare_batch_step(
    step_index: int,
    step_spec: Any,
    *,
    precedent_db_path: str | None,
    use_api: str | None,
    min_similarity: float,
    limit: int,
    include_patents: bool,
    default_keywords: list[str] | None,
    max_patents: int,
) -> _PreparedBatchStep:
  """Validates a batch step mapping and resolves overrides against defaults."""
  if not isinstance(step_spec, Mapping):
    raise _BatchValidationError(
        f"Step {step_index} must be a JSON object/mapping, got"
        f" {type(step_spec).__name__}."
    )

  raw_rxn = step_spec.get("reaction_smiles")
  if raw_rxn is not None and not isinstance(raw_rxn, str):
    raise _BatchValidationError(
        f"Step {step_index} 'reaction_smiles' must be a string."
    )
  rxn_smiles = raw_rxn.strip() if isinstance(raw_rxn, str) else ""
  resolved_rxn_smiles = rxn_smiles or None

  raw_reactants = step_spec.get("reactants")
  if raw_reactants is None:
    step_reactants: list[str] = []
  elif isinstance(raw_reactants, str):
    if raw_reactants.strip():
      step_reactants = [raw_reactants.strip()]
    else:
      step_reactants = []
  elif isinstance(raw_reactants, Sequence) and not isinstance(
      raw_reactants, (bytes, bytearray)
  ):
    step_reactants = []
    for item in raw_reactants:
      if not isinstance(item, str):
        raise _BatchValidationError(
            f"Step {step_index} 'reactants' entries must be strings."
        )
      if item.strip():
        step_reactants.append(item.strip())
  else:
    raise _BatchValidationError(
        f"Step {step_index} 'reactants' must be a SMILES string or list of"
        " SMILES strings."
    )

  raw_product = step_spec.get("product")
  if raw_product is not None and not isinstance(raw_product, str):
    raise _BatchValidationError(
        f"Step {step_index} 'product' must be a string."
    )
  step_product = (
      raw_product.strip() or None if isinstance(raw_product, str) else None
  )

  _validate_query_smiles_inputs(
      step_reactants,
      step_product,
      resolved_rxn_smiles,
      prefix=f"Step {step_index}: ",
  )
  split_reactants: list[str] = []
  for r in step_reactants:
    split_reactants.extend(_split_smiles_fragments(r))
  step_reactants = split_reactants

  raw_db = step_spec.get("precedent_db")
  if raw_db is None:
    step_db = precedent_db_path
  elif isinstance(raw_db, str) and raw_db.strip():
    step_db = raw_db.strip()
  else:
    raise _BatchValidationError(
        f"Step {step_index} 'precedent_db' must be a non-empty string."
    )

  raw_api = step_spec.get("api")
  if raw_api is not None and (
      not isinstance(raw_api, str) or raw_api not in ("ord", "patents", "all")
  ):
    raise _BatchValidationError(
        f"Step {step_index} 'api' must be 'ord', 'patents', or 'all', got"
        f" {raw_api!r}."
    )

  raw_inc = step_spec.get("include_patents")
  if raw_inc is not None and not isinstance(raw_inc, bool):
    raise _BatchValidationError(
        f"Step {step_index} 'include_patents' must be a boolean, got"
        f" {raw_inc!r}."
    )

  if raw_api is not None or raw_inc is not None:
    step_use_api, step_include_patents = _resolve_cli_source_flags(
        raw_api if raw_api is not None else use_api,
        raw_inc,
    )
  else:
    step_use_api, step_include_patents = use_api, include_patents

  raw_min_sim = step_spec.get("min_similarity")
  if raw_min_sim is None:
    step_min_sim = float(min_similarity)
  elif isinstance(raw_min_sim, bool) or not isinstance(
      raw_min_sim, (int, float)
  ):
    raise _BatchValidationError(
        f"Step {step_index} 'min_similarity' must be a number between 0.0 and"
        f" 1.0, got {raw_min_sim!r}."
    )
  else:
    step_min_sim = float(raw_min_sim)
    if math.isnan(step_min_sim) or not (0.0 <= step_min_sim <= 1.0):
      raise _BatchValidationError(
          f"Step {step_index} 'min_similarity' must be between 0.0 and 1.0,"
          f" got {raw_min_sim!r}."
      )

  raw_limit = step_spec.get("limit")
  if raw_limit is None:
    step_limit = int(limit)
  elif (
      isinstance(raw_limit, bool)
      or not isinstance(raw_limit, int)
      or raw_limit < 1
  ):
    raise _BatchValidationError(
        f"Step {step_index} 'limit' must be a positive integer, got"
        f" {raw_limit!r}."
    )
  else:
    step_limit = raw_limit

  raw_max_patents = step_spec.get("max_patents")
  if raw_max_patents is None:
    step_max_patents = int(max_patents)
  elif (
      isinstance(raw_max_patents, bool)
      or not isinstance(raw_max_patents, int)
      or raw_max_patents < 1
  ):
    raise _BatchValidationError(
        f"Step {step_index} 'max_patents' must be a positive integer, got"
        f" {raw_max_patents!r}."
    )
  else:
    step_max_patents = raw_max_patents

  raw_kw = step_spec.get("keywords")
  if raw_kw is None:
    step_kw = list(default_keywords) if default_keywords is not None else None
  elif isinstance(raw_kw, str):
    step_kw = [k.strip() for k in raw_kw.split(",") if k.strip()]
  elif isinstance(raw_kw, Sequence) and not isinstance(
      raw_kw, (bytes, bytearray)
  ):
    step_kw = []
    for kw_item in raw_kw:
      if not isinstance(kw_item, str):
        raise _BatchValidationError(
            f"Step {step_index} 'keywords' entries must be strings."
        )
      if kw_item.strip():
        step_kw.append(kw_item.strip())
  else:
    raise _BatchValidationError(
        f"Step {step_index} 'keywords' must be a string or list of strings."
    )

  metadata = {
      meta_key: step_spec[meta_key]
      for meta_key in ("proposal_id", "step_id", "step_number", "title")
      if meta_key in step_spec
  }

  return _PreparedBatchStep(
      reactants=step_reactants,
      product=step_product,
      reaction_smiles=resolved_rxn_smiles,
      precedent_db_path=step_db,
      use_api=step_use_api,
      min_similarity=step_min_sim,
      limit=step_limit,
      include_patents=step_include_patents,
      keywords=step_kw,
      max_patents=step_max_patents,
      metadata=metadata,
  )


def lookup_reaction_precedents_batch(
    steps: Sequence[Mapping[str, Any]],
    *,
    precedent_db_path: str | None = None,
    use_api: str | None = None,
    min_similarity: float = 0.35,
    limit: int = 25,
    include_patents: bool = False,
    keywords: Sequence[str] | None = None,
    max_patents: int = 8,
) -> dict[str, Any]:
  """Looks up precedents across multiple pathway steps with shared caches."""
  if not isinstance(steps, Sequence) or isinstance(
      steps, (str, bytes, bytearray)
  ):
    raise _BatchValidationError("steps must be a sequence of step mappings.")
  if not steps:
    return {"step_count": 0, "steps": []}

  default_kw = [str(k) for k in keywords if k] if keywords else None
  prepared_steps = [
      _prepare_batch_step(
          idx,
          step_spec,
          precedent_db_path=precedent_db_path,
          use_api=use_api,
          min_similarity=min_similarity,
          limit=limit,
          include_patents=include_patents,
          default_keywords=default_kw,
          max_patents=max_patents,
      )
      for idx, step_spec in enumerate(steps, start=1)
  ]

  def _run_step(prepared: _PreparedBatchStep) -> dict[str, Any]:
    res = lookup_reaction_precedents(
        reactants=prepared.reactants,
        product=prepared.product,
        reaction_smiles=prepared.reaction_smiles,
        precedent_db_path=prepared.precedent_db_path,
        use_api=prepared.use_api,
        min_similarity=prepared.min_similarity,
        limit=prepared.limit,
        include_patents=prepared.include_patents,
        keywords=prepared.keywords,
        max_patents=prepared.max_patents,
    )
    return {**prepared.metadata, **res}

  workers = min(_MAX_BATCH_STEP_WORKERS, len(prepared_steps))
  with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
    step_results = list(pool.map(_run_step, prepared_steps))

  return {
      "step_count": len(step_results),
      "steps": step_results,
  }


# ==============================================================================
# Formatted Terminal & Markdown Reporting
# ==============================================================================


def print_precedent_report(results: Mapping[str, Any]) -> None:
  """Prints a terminal report of precedent results and condition clusters."""
  query = results["query"]
  print("\n" + "=" * 80)
  print(" REACTION PRECEDENT LOOKUP REPORT")
  print("=" * 80)
  print(f"Query Reaction:   {query['reaction_smiles']}")
  print(f"Precedents Found: {results['precedent_count']} matching records\n")

  api_errors = results.get("errors") or []
  if api_errors:
    print("Backend / API Errors:")
    for err in api_errors:
      print(f"  - {err}")
    print()

  if results["precedent_count"] == 0:
    print("No matching precedents found.")
    print(
        "Note: Provide a database via --precedent-db <path.parquet|jsonl|csv>"
        " or query live APIs via --api all."
    )
    print("=" * 80 + "\n")
    return

  print("-" * 80)
  print(" RANKED REACTION CONDITION ARCHETYPES")
  print("-" * 80)

  for idx, c in enumerate(results["clusters"], start=1):
    print(f"\n{idx}. {c['cluster_name']}")
    print(f"   Reagents:       {c['reagents']} | Solvent: {c['solvent']}")

    yield_str = (
        f"{c['median_yield']:.1f}%" if c["median_yield"] is not None else "N/A"
    )
    if c.get("yield_range"):
      yield_str += (
          f" (range: {c['yield_range'][0]:.0f}-{c['yield_range'][1]:.0f}%)"
      )

    time_str = (
        f"{c['median_time_hours']:.1f} h"
        if c["median_time_hours"] is not None
        else "N/A"
    )
    if c.get("time_range_hours"):
      time_str += (
          f" ({c['time_range_hours'][0]:.1f}-{c['time_range_hours'][1]:.1f} h)"
      )

    print(
        f"   Performance:    Median Yield: {yield_str} | Typical Time:"
        f" {time_str}"
    )
    cit_preview = ", ".join(c["citations"][:4])
    if len(c["citations"]) > 4:
      cit_preview += f" (+{len(c['citations']) - 4} more)"
    print(f"   Citations (N={c['precedent_count']}): {cit_preview}")

  print("\n" + "=" * 80 + "\n")


def _format_step_markdown_section(
    results: Mapping[str, Any],
    *,
    nested: bool = False,
) -> str:
  """Formats a single reaction's precedent archetypes into Markdown."""
  heading = (
      "### Ranked Reaction Condition Archetypes\n"
      if nested
      else "## Ranked Reaction Condition Archetypes\n"
  )
  lines = [
      f"**Query Reaction:** `{results['query']['reaction_smiles']}`\n",
      f"**Precedents Retrieved:** {results['precedent_count']}\n",
      heading,
      (
          "| Rank | Conditions | Median Yield | Typical Time | Support (N) |"
          " Citations |"
      ),
      "| :--- | :--- | :--- | :--- | :--- | :--- |",
  ]
  for idx, c in enumerate(results["clusters"], start=1):
    y_str = (
        f"{c['median_yield']:.1f}%" if c["median_yield"] is not None else "N/A"
    )
    t_str = (
        f"{c['median_time_hours']:.1f} h"
        if c["median_time_hours"] is not None
        else "N/A"
    )
    cits = ", ".join(c["citations"])
    lines.append(
        f"| {idx} | {c['cluster_name']} | {y_str} | {t_str} |"
        f" {c['precedent_count']} | {cits} |"
    )
  return "\n".join(lines) + "\n"


def _write_markdown_report(
    results: Mapping[str, Any],
    output_path: str,
) -> None:
  """Writes single-step or multi-step precedent results to a Markdown file."""
  os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
  with open(output_path, "w", encoding="utf-8") as f:
    f.write("# Reaction Precedent Lookup Report\n\n")
    if "steps" in results:
      for idx, step_res in enumerate(results["steps"], start=1):
        label = (
            step_res.get("title") or step_res.get("step_id") or f"Step {idx}"
        )
        f.write(f"## {label}\n\n")
        f.write(_format_step_markdown_section(step_res, nested=True))
        f.write("\n")
    else:
      f.write(_format_step_markdown_section(results))
  print(f"Exported Markdown report to {output_path}")


# ==============================================================================
# Command-Line Interface
# ==============================================================================


def main() -> None:
  parser = argparse.ArgumentParser(
      description=(
          "Reaction Precedent Lookup Engine: queries ORD v2 and"
          " PubChem/Google Patents concurrently, computes vectorized batch"
          " DRFP fingerprint similarity, and clusters condition archetypes."
      )
  )
  parser.add_argument(
      "--reactants",
      nargs="*",
      default=None,
      help=(
          "One or more reactant SMILES strings (e.g. 'CC(O)c1cccnc1'"
          " 'Fc1ccccc1[N+](=O)[O-]')."
      ),
  )
  parser.add_argument(
      "--product",
      default=None,
      help=(
          "Target product SMILES string (pass with --reactants for a specific"
          " transformation, or alone to look up existing syntheses of the"
          " target)."
      ),
  )
  parser.add_argument(
      "--reaction-smiles",
      default=None,
      help="Reaction SMILES string in 'reactant.reactant>>product' format.",
  )
  parser.add_argument(
      "--batch-json",
      default=None,
      help=(
          "Path to JSON list (or {'steps': [...]}) of multi-step reactions to"
          " evaluate in a single process with shared in-memory PubChem and"
          " Google Patents HTML caches."
      ),
  )
  parser.add_argument(
      "--precedent-db",
      default=None,
      help=(
          "Path to local reaction precedent database (.parquet, .jsonl,"
          " .json, .csv)."
      ),
  )
  parser.add_argument(
      "--api",
      default=None,
      choices=["ord", "patents", "all"],
      help=(
          "Query public reaction database APIs directly ('ord' for ORD/USPTO"
          " v2 API, 'patents' for PubChem + Google Patents, or 'all' for"
          " concurrent execution of both)."
      ),
  )
  parser.add_argument(
      "--include-patents",
      dest="include_patents",
      action="store_true",
      default=None,
      help=(
          "Query PubChem + full-text Google Patents specifications (enabled by"
          " default when --api is 'all' or 'patents')."
      ),
  )
  parser.add_argument(
      "--no-include-patents",
      dest="include_patents",
      action="store_false",
      help="Disable PubChem / Google Patents precedent search.",
  )
  parser.add_argument(
      "--keywords",
      nargs="*",
      default=None,
      help=(
          "Optional reaction, catalyst, or enzyme keywords to prioritize"
          " experimental paragraphs in Google Patents (e.g. 'amination'"
          " 'P411')."
      ),
  )
  parser.add_argument(
      "--max-patents",
      type=int,
      default=8,
      help=(
          "Maximum number of candidate Google Patents specifications to"
          " inspect (default: 8)."
      ),
  )
  parser.add_argument(
      "--min-similarity",
      type=float,
      default=0.35,
      help=(
          "Minimum reaction fingerprint Tanimoto similarity threshold (default:"
          " 0.35)."
      ),
  )
  parser.add_argument(
      "--limit",
      type=int,
      default=50,
      help=(
          "Maximum number of individual precedent records to retrieve (default:"
          " 50)."
      ),
  )
  parser.add_argument(
      "--output-json",
      default=None,
      help="Destination path for structured JSON results export.",
  )
  parser.add_argument(
      "--output-markdown",
      default=None,
      help="Destination path for markdown summary report.",
  )

  args = parser.parse_args()
  effective_api, include_patents = _resolve_cli_source_flags(
      args.api, args.include_patents
  )

  if args.batch_json is not None:
    try:
      with open(args.batch_json, "r", encoding="utf-8") as f:
        raw_batch = json.load(f)
      step_list = _extract_batch_steps(raw_batch)
      results: dict[str, Any] = lookup_reaction_precedents_batch(
          steps=step_list,
          precedent_db_path=args.precedent_db,
          use_api=effective_api,
          min_similarity=args.min_similarity,
          limit=args.limit,
          include_patents=include_patents,
          keywords=args.keywords,
          max_patents=args.max_patents,
      )
    except (OSError, json.JSONDecodeError, _BatchValidationError) as exc:
      parser.error(f"Invalid --batch-json ({args.batch_json}): {exc}")
    for step_res in results["steps"]:
      print_precedent_report(step_res)
  else:
    reactants = [r.strip() for r in args.reactants or () if r and r.strip()]
    product = (
        args.product.strip() if args.product and args.product.strip() else None
    )
    rxn_smiles = (
        args.reaction_smiles.strip()
        if args.reaction_smiles and args.reaction_smiles.strip()
        else None
    )
    if not reactants and product is None and rxn_smiles is None:
      parser.error(
          "Must provide --reactants, --product, --reaction-smiles, or"
          " --batch-json."
      )
    try:
      _validate_query_smiles_inputs(reactants, product, rxn_smiles)
    except _BatchValidationError as exc:
      parser.error(str(exc))
    results = lookup_reaction_precedents(
        reactants=reactants,
        product=product,
        reaction_smiles=rxn_smiles,
        precedent_db_path=args.precedent_db,
        use_api=effective_api,
        min_similarity=args.min_similarity,
        limit=args.limit,
        include_patents=include_patents,
        keywords=args.keywords,
        max_patents=args.max_patents,
    )
    print_precedent_report(results)

  if args.output_json:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_json)), exist_ok=True
    )
    with open(args.output_json, "w", encoding="utf-8") as f:
      json.dump(results, f, indent=2)
    print(f"Exported JSON results to {args.output_json}")

  if args.output_markdown:
    _write_markdown_report(results, args.output_markdown)

  step_results = results.get("steps", [results])
  if any(s.get("errors") and s["precedent_count"] == 0 for s in step_results):
    sys.exit(1)


if __name__ == "__main__":
  main()
