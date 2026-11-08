# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "pandas",
#     "scipy",
#     "scikit-learn",
#     "matplotlib",
# ]
# ///

"""agent UI-Compatible Interactive 2D Chemical Grid & SAR Scatter Plot Viewer (Zero-CDN, Self-Contained, CSP-Compliant).

Features:
- Dual-View Architecture: Responsive 2D Grid View + Interactive SAR Scatter Plot
(MW vs. cLogP with Ro5 Sweet Spot).
- Property Badges: Rule of 5 (Ro5), PAINS structural alerts, and Pfizer CNS MPO
badges on cards and modals.
- Scaffold-Aligned 2D Depictions (--align_core): Aligns all derivative
depictions to a common core scaffold.
- Instant Modal Pop-Out Zoom: Viewport-centered, zero-lag (<16ms) interactive
inspection with clipboard copy.
- Zero External CDN / 100% CSP-Compliant: Standalone single-file HTML keeping
total size strictly under 5MB.
- Instant Client-Side CSV Export: Download filtered compound sets directly from
the browser.
"""

import argparse
import base64
from collections.abc import Iterable
import datetime
import gzip
import html
import json
import math
import os
import re
import sys
import tempfile
from typing import Any, Final, Optional
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import PandasTools
from rdkit.Chem import QED
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

# Canonical CNS MPO implementation lives in cheminformatics_core.
_script_dir = os.path.dirname(os.path.abspath(__file__))
if _script_dir not in sys.path:
  sys.path.insert(0, _script_dir)

from draw_molecules import compute_clean_2d_coords
from skill_paths import add_skill_scripts

add_skill_scripts("cheminformatics_core")

from descriptors import calculate_cns_mpo

_METAL_ATOMIC_NUMS: Final[frozenset[int]] = frozenset({
    3,
    4,
    11,
    12,
    13,
    19,
    20,
    *range(21, 32),
    37,
    38,
    *range(39, 51),
    55,
    56,
    *range(57, 85),
    87,
    88,
    *range(89, 117),
})
_ANON_PATENT_COL_RE: Final[re.Pattern[str]] = re.compile(
    r"^(?:patent_table|pat|[a-z]{2,4}\d{5,}[a-z0-9]*)_col_\d+"
    r"(?:_(?:nm|pactivity|range))?$",
    re.I,
)


def _format_col_label(col_s: str) -> str:
  """Formats snake_case column names into space-separated UI labels."""
  cleaned = re.sub(r"[_\s]+", " ", str(col_s)).strip()
  return cleaned or str(col_s).strip()


def _parse_grid_smiles(raw_val: Any) -> Optional[Chem.Mol]:
  """Parses SMILES into an RDKit Mol, repairing dative and metal valence."""
  if not pd.notna(raw_val):
    return None
  smi_str = str(raw_val).strip()
  if not smi_str or smi_str.lower() == "nan":
    return None
  with rdBase.BlockLogs():
    mol = Chem.MolFromSmiles(smi_str)
    if mol is not None:
      return mol
    ps_raw = Chem.SmilesParserParams()
    ps_raw.sanitize = False
    ps_raw.removeHs = False
    raw = Chem.MolFromSmiles(smi_str, ps_raw)
    if raw is None:
      return None
    raw.UpdatePropertyCache(strict=False)
    for b in list(raw.GetBonds()):
      if b.GetBondType() != Chem.BondType.DATIVE:
        continue
      da, aa = b.GetBeginAtom(), b.GetEndAtom()
      if (
          da.GetAtomicNum() not in _METAL_ATOMIC_NUMS
          and aa.GetAtomicNum() not in _METAL_ATOMIC_NUMS
          and da.GetFormalCharge() == 0
          and aa.GetFormalCharge() == 0
      ):
        dh, ah = da.GetTotalNumHs(), aa.GetTotalNumHs()
        b.SetBondType(Chem.BondType.SINGLE)
        da.SetFormalCharge(1)
        aa.SetFormalCharge(-1)
        da.SetNoImplicit(True)
        da.SetNumExplicitHs(dh)
        aa.SetNoImplicit(True)
        aa.SetNumExplicitHs(ah)
    problems = Chem.DetectChemistryProblems(raw)
    if not problems:
      try:
        Chem.SanitizeMol(raw)
      except (RuntimeError, ValueError):
        return None
    elif all(
        p.GetType() == "AtomValenceException"
        and p.GetAtomIdx() >= 0
        and raw.GetAtomWithIdx(p.GetAtomIdx()).GetAtomicNum()
        in _METAL_ATOMIC_NUMS
        for p in problems
    ):
      try:
        Chem.SanitizeMol(
            raw,
            sanitizeOps=(
                Chem.SanitizeFlags.SANITIZE_ALL
                ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
            ),
        )
        raw.UpdatePropertyCache(strict=False)
      except (RuntimeError, ValueError):
        return None
    else:
      return None
    out = Chem.RemoveHs(raw, sanitize=False)
    out.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(out)
    Chem.AssignStereochemistry(out, cleanIt=False, force=True)
    return out


def _safe_qed(mol: Chem.Mol) -> float:
  """Computes rounded QED drug-likeness, returning NaN on RDKit failure."""
  with rdBase.BlockLogs():
    try:
      return round(QED.qed(mol), 2)
    except (RuntimeError, ValueError, ZeroDivisionError):
      return np.nan


def get_pains_catalog() -> FilterCatalog:
  """Builds an RDKit FilterCatalog with PAINS (A/B/C) and Brenk alerts."""
  params = FilterCatalogParams()
  params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
  params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
  params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
  params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
  return FilterCatalog(params)


def generate_svg_string(
    mol: Chem.Mol,
    width: int = 250,
    height: int = 180,
    core_mol: Optional[Chem.Mol] = None,
) -> str:
  """Generates crisp standalone SVG string with transparent background and optional core alignment."""
  if mol is None:
    return (
        '<div style="color:#d93025; padding:40px 0; text-align:center;'
        ' font-size:12px;">Invalid Structure</div>'
    )
  try:
    compute_clean_2d_coords(mol, core_mol=core_mol)

    d2d = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = d2d.drawOptions()
    opts.clearBackground = False  # Transparent background
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    svg = d2d.GetDrawingText()

    # Clean XML header, comments, and unneeded namespaces
    svg = re.sub(r"<\?xml.*?\?>", "", svg)
    svg = re.sub(r"<!--.*?-->", "", svg, flags=re.DOTALL)
    svg = re.sub(r"xmlns:rdkit=[\x27\x22].*?[\x27\x22]", "", svg)
    svg = re.sub(r"xmlns:xlink=[\x27\x22].*?[\x27\x22]", "", svg)
    svg = re.sub(r"xml:space=[\x27\x22].*?[\x27\x22]", "", svg)
    svg = re.sub(r"version=[\x27\x22].*?[\x27\x22]", "", svg)
    svg = re.sub(r"baseProfile=[\x27\x22].*?[\x27\x22]", "", svg)
    svg = re.sub(r" class=[\x27\x22][^\x27\x22]*[\x27\x22]", "", svg)
    svg = (
        svg.replace("fill-rule:evenodd;", "")
        .replace("stroke-linecap:butt;", "")
        .replace("stroke-linejoin:miter;", "")
        .replace("stroke-opacity:1", "")
    )
    svg = re.sub(r"(\d+\.\d)\d+", r"\1", svg)
    svg = re.sub(r"\s+", " ", svg)
    svg = re.sub(r">\s+<", "><", svg)
    if "<svg" in svg:
      svg = svg[svg.find("<svg") :]
    return svg.strip()
  except Exception as e:
    return (
        '<div style="color:#d93025; padding:20px; font-size:11px;">Render'
        f" error: {html.escape(str(e))}</div>"
    )


_INACTIVE_POTENCY_TOKENS = frozenset({
    "",
    "-",
    "--",
    "—",
    "NA",
    "N/A",
    "ND",
    "N.D.",
    "NT",
    "N.T.",
    "NONE",
    "NULL",
    "NAN",
    "INACTIVE",
    "NOT DETERMINED",
    "NOT TESTED",
})

_LETTER_BIN_TIERS = {
    "A": 4.0,
    "B": 3.0,
    "C": 2.0,
    "D": 1.0,
    "E": 0.5,
    "F": 0.25,
}

_UNIT_FACTORS: Final[dict[str, float]] = {
    "fm": 1e-15,
    "pm": 1e-12,
    "nm": 1e-9,
    "um": 1e-6,
    "µm": 1e-6,
    "μm": 1e-6,
    "mm": 1e-3,
    "m": 1.0,
}
_UNIT_CANONICAL_LABELS: Final[dict[str, str]] = {
    "fm": "fM",
    "pm": "pM",
    "nm": "nM",
    "um": "uM",
    "µm": "uM",
    "μm": "uM",
    "mm": "mM",
    "m": "M",
}

_CONC_UNIT_PAT: Final[str] = (
    r"(?:f|p|n|u|µ|μ|m)?mol\s*"
    r"(?:/\s*l|[\s._·]l(?:\s*(?:\^?\s*[-−]1|⁻¹)|_1)?|l[-−]1)|"
    r"molar|fm|pm|nm|um|µm|μm|mm|m"
)
_NORM_UNIT_PAT: Final[str] = (
    r"(?:f|p|n|u|µ|μ|m)?mol_l(?:_[1¹])?|molar|fm|pm|nm|um|µm|μm|mm"
)
_RAW_CONC_PAT: Final[str] = r"(?:ic|ec|ac|gi)_?50|ki|kd|mic"
_LOG_POTENCY_PAT: Final[str] = (
    r"p(?:ic|ec|ac|gi)_?50|pki|pkd|pmic|pchembl(?:_value)?|pactivity|pact|"
    rf"log(?:_?10)?_(?:1_)?(?:{_RAW_CONC_PAT})"
)
_POTENCY_META_TOKENS: Final[frozenset[str]] = frozenset({
    "unit",
    "units",
    "qualifier",
    "operator",
    "relation",
    "flag",
    "comment",
    "range",
})


def _norm_col(col_name: str) -> str:
  """Normalizes column separators to single lowercase underscores."""
  return re.sub(r"[\W_]+", "_", col_name.lower()).strip("_")


def _normalize_unit_token(unit_str: str) -> str:
  """Canonicalizes short molar symbols and SI mol/L spellings."""
  s = (
      unit_str.strip()
      .lower()
      .replace("µ", "μ")
      .replace("·", ".")
      .replace("−", "-")
      .replace("⁻¹", "-1")
  )
  if s == "molar":
    return "m"
  m_mol = re.match(
      r"^(f|p|n|u|μ|m)?mol(?:\s*(?:/|[\s._])\s*l(?:\s*(?:\^?\s*-1|_1))?|l-1)$",
      s,
  )
  if m_mol:
    return f"{m_mol.group(1) or ''}m"
  return s


def _normalize_potency_col_tokens(col_name: str) -> str:
  """Collapses split `(p)IC_50` / `IC-50` / `IC 50` tokens into `(p)ic50`."""
  return re.sub(
      r"(?<![a-z0-9])(p?(?:ic|ec|ac|gi))[\s_-]?50(?![a-z0-9])",
      r"\g<1>50",
      col_name,
      flags=re.IGNORECASE,
  )


def _conc_col_molar_unit(
    col_name: str, unit_hint: str = ""
) -> tuple[float, str]:
  """Returns `(molar_factor, canonical_label)` for a concentration column."""
  if unit_hint.strip():
    hint_key = _normalize_unit_token(unit_hint)
    if hint_key in _UNIT_FACTORS:
      return _UNIT_FACTORS[hint_key], _UNIT_CANONICAL_LABELS[hint_key]
  norm = _normalize_potency_col_tokens(col_name)
  all_units: list[str] = []
  readout_units: list[str] = []
  for m in re.finditer(
      rf"(?:^|[^a-z0-9µμ]|(?:ic50|ki|kd|ec50|ac50|gi50|mic))"
      rf"({_CONC_UNIT_PAT})(?=[^a-z0-9]|$)",
      norm,
      re.IGNORECASE,
  ):
    raw_u = m.group(1)
    u_start = m.start(1)
    prefix = norm[:u_start]
    if len(raw_u) == 1:
      prev_ch = norm[u_start - 1] if u_start > 0 else ""
      has_pot_pre = bool(
          re.search(
              r"(?:^|[^a-z0-9])(?:p?(?:ic|ec|ac|gi)50|ki|kd|mic|"
              r"potency|activity|affinity|conc\w*|value)(?:[^a-z0-9]|$)",
              prefix,
              re.IGNORECASE,
          )
      )
      if (
          raw_u != "M"
          or not prev_ch
          or prev_ch.isalnum()
          or prev_ch in "µμ"
          or not (
              prev_ch in "([{" or not norm[m.end(1) :].strip() or has_pot_pre
          )
      ):
        continue
    unit = _normalize_unit_token(raw_u)
    if unit not in _UNIT_FACTORS:
      continue
    all_units.append(unit)
    if not re.search(
        r"(?<![a-z0-9<>≤≥.=])\d+(?:\.\d+)?[\s_]*$", prefix, re.IGNORECASE
    ):
      readout_units.append(unit)
  chosen = (
      readout_units[-1]
      if readout_units
      else (all_units[-1] if all_units else None)
  )
  if chosen:
    return _UNIT_FACTORS[chosen], _UNIT_CANONICAL_LABELS[chosen]
  return 1e-9, "nM"


_NUM_LIT_PAT: Final[str] = r"(?:\d+(?:\.\d+)?|\.\d+)(?:e[+-]?\d+)?"


def _parse_conc_string_molar(
    s: str, default_molar: float | None = None
) -> tuple[float | None, str]:
  """Parses a concentration range or inequality into `(molar_value, op)`."""
  s_num = re.sub(r"(?<=\d),(?=\d{3}\b)", "", s.strip())
  m_rng = re.search(
      rf"({_NUM_LIT_PAT})\s*({_CONC_UNIT_PAT})?\s*(?:[-–—]|to)\s*"
      rf"({_NUM_LIT_PAT})\s*({_CONC_UNIT_PAT})(?![a-z0-9])",
      s_num,
      re.IGNORECASE,
  )
  if m_rng and not (len(m_rng.group(4)) == 1 and m_rng.group(4) != "M"):
    u2 = _UNIT_FACTORS.get(_normalize_unit_token(m_rng.group(4)), 1e-9)
    raw_u1 = m_rng.group(2)
    u1 = (
        _UNIT_FACTORS.get(_normalize_unit_token(raw_u1), u2)
        if raw_u1 and not (len(raw_u1) == 1 and raw_u1 != "M")
        else u2
    )
    m1, m2 = float(m_rng.group(1)) * u1, float(m_rng.group(3)) * u2
    return (math.sqrt(m1 * m2), "") if (m1 > 0 and m2 > 0) else (None, "")
  m_ineq = re.match(
      rf"^([<>≤≥~=]*)\s*({_NUM_LIT_PAT})\s*({_CONC_UNIT_PAT})(?![a-z0-9])",
      s_num,
      re.IGNORECASE,
  )
  if m_ineq and not (len(m_ineq.group(3)) == 1 and m_ineq.group(3) != "M"):
    unit_key = _normalize_unit_token(m_ineq.group(3))
    v_m = float(m_ineq.group(2)) * _UNIT_FACTORS.get(unit_key, 1e-9)
    return (v_m, m_ineq.group(1)) if v_m > 0 else (None, "")
  if default_molar is not None:
    m_bare = re.match(rf"^([<>≤≥~=]+)\s*({_NUM_LIT_PAT})\s*$", s_num, re.I)
    if m_bare:
      v_m = float(m_bare.group(2)) * default_molar
      return (v_m, m_bare.group(1)) if v_m > 0 else (None, "")
  return None, ""


def _parse_discrete_potency_tier(val: Any) -> Optional[float]:
  """Maps discrete potency bins to ordinal scores (higher = more potent)."""
  if val is None or pd.isna(val):
    return None
  s = str(val).strip()
  if not s or s.upper() in _INACTIVE_POTENCY_TOKENS:
    return None
  # 1. Letter bin with optional range/label, e.g. "A", "A (<= 100 nM)", "Bin C"
  m_let = re.match(
      r"^(?:(?:Grade|Bin|Category|Class|Tier)\s+)?([A-F])(?:\s*\(.*\))?$",
      s,
      re.IGNORECASE,
  )
  if m_let:
    return _LETTER_BIN_TIERS.get(m_let.group(1).upper())
  # 2. Plus-sign or star potency grades: ++++, +++, ++, +, ****, ***, **, *
  if re.match(r"^(\+{1,6}|\*{1,6})$", s):
    return float(len(s))
  # 3. Concentration range or inequality (e.g. "<= 100 nM", "100-250 nM")
  v_m, op = _parse_conc_string_molar(s)
  if v_m is not None:
    base_p = -math.log10(v_m)
    if "<" in op or "≤" in op:
      return round(base_p + 0.15, 3)
    if ">" in op or "≥" in op:
      return round(base_p - 0.15, 3)
    return round(base_p, 3)
  return None


def _unit_molar_factor(col_name: str, unit_hint: str = "") -> float | None:
  """Resolves molar conversion factor from an explicit unit or column suffix."""
  if unit_hint.strip():
    return _UNIT_FACTORS.get(_normalize_unit_token(unit_hint))
  return _conc_col_molar_unit(col_name)[0]


def _is_potency_meta_col(col_name: str) -> bool:
  """Returns True if col_name is a qualifier/unit/range/type metadata column."""
  norm = _norm_col(col_name)
  if set(norm.split("_")) & _POTENCY_META_TOKENS:
    return True
  return norm == "type" or (
      norm.endswith("_type") and not norm.endswith(("wild_type", "cell_type"))
  )


def _is_ratio_or_percent_col(col_name: str) -> bool:
  """Returns True if col_name is a percentage or fold/selectivity ratio."""
  if "%" in col_name:
    return True
  return bool(
      re.search(
          r"(?:^|_)(?:fold|ratio|selectivity|shift|pct|percent|emax)(?:_|$)",
          _norm_col(col_name),
      )
  )


def _is_potency_col_name(col_name: str) -> bool:
  """Returns True if col_name names a potency or bioactivity metric."""
  if _is_potency_meta_col(col_name) or _is_ratio_or_percent_col(col_name):
    return False
  return bool(
      re.search(
          rf"(?:^|_)(?:{_LOG_POTENCY_PAT}|{_RAW_CONC_PAT}|potency|"
          rf"activity(?:_value|_bin)?|standard_value)"
          rf"(?:_|$|(?:{_NORM_UNIT_PAT})(?:_|$))",
          _norm_col(col_name),
      )
  )


def _is_raw_conc_column(col_name: str) -> bool:
  """Returns True if higher values in col_name mean weaker potency."""
  if _is_ratio_or_percent_col(col_name):
    return False
  norm = _norm_col(col_name)
  if re.search(
      rf"(?:^|_)(?:{_LOG_POTENCY_PAT})(?:_|$|(?:{_NORM_UNIT_PAT})(?:_|$))",
      norm,
  ):
    return False
  if re.search(
      rf"(?:^|_)(?:{_RAW_CONC_PAT})(?:_|$|(?:{_NORM_UNIT_PAT})(?:_|$))",
      norm,
  ) or norm in ("activity_value", "standard_value"):
    return True
  if re.search(r"(?:^|_)(?:inhib\w*|viab\w*)(?:_|$)", norm):
    return False
  m_tail = re.search(rf"_({_NORM_UNIT_PAT})$", norm)
  if m_tail:
    prefix = norm[: m_tail.start(1)]
    if not re.search(r"(?<![a-z0-9])\d+_$", prefix):
      return True
  return False


def _resolve_act_companions(
    act_col: str, columns: Iterable[Any]
) -> tuple[str, str, tuple[str, ...]]:
  """Returns companion (type_col, units_col, range_cols) scoped to act_col."""
  col_by_norm: dict[str, str] = {}
  for c in columns:
    col_by_norm.setdefault(_norm_col(str(c)), str(c))
  norm_act = _norm_col(act_col)
  if norm_act in ("activity_value", "activity_bin", "activity"):
    type_norms = ("activity_type", f"{norm_act}_type")
    unit_norms = (
        "activity_units",
        "activity_unit",
        f"{norm_act}_units",
        f"{norm_act}_unit",
    )
    range_norms = (f"{norm_act}_range", "activity_range")
  elif norm_act == "standard_value":
    type_norms = ("standard_type", f"{norm_act}_type")
    unit_norms = (
        "standard_units",
        "standard_unit",
        f"{norm_act}_units",
        f"{norm_act}_unit",
    )
    range_norms = (f"{norm_act}_range", "standard_range")
  else:
    type_norms = (f"{norm_act}_type",)
    unit_norms = (f"{norm_act}_units", f"{norm_act}_unit")
    range_norms = (f"{norm_act}_range",)
  type_col = next((col_by_norm[k] for k in type_norms if k in col_by_norm), "")
  units_col = next((col_by_norm[k] for k in unit_norms if k in col_by_norm), "")
  range_cols = tuple(
      dict.fromkeys(col_by_norm[k] for k in range_norms if k in col_by_norm)
  )
  return (type_col, units_col, range_cols)


def _series_potency_scores(
    series: pd.Series,
    col_name: str = "",
    units_series: pd.Series | None = None,
) -> pd.Series:
  """Converts numeric or discrete potencies to ordinal potency scores."""
  num_s = pd.to_numeric(series, errors="coerce")
  finite_mask = num_s.notna() & np.isfinite(num_s)
  is_conc = _is_raw_conc_column(col_name)
  out = pd.Series(np.nan, index=series.index, dtype=float)
  valid_num = (finite_mask & (num_s > 0)) if is_conc else finite_mask
  if valid_num.any():
    if is_conc:
      if units_series is not None:
        factors = units_series.loc[valid_num].apply(
            lambda u: _unit_molar_factor(col_name, _as_clean_str(u))
        )
      else:
        factors = pd.Series(
            _unit_molar_factor(col_name), index=num_s.loc[valid_num].index
        )
      has_factor = factors.notna()
      if has_factor.any():
        idx_f = factors.index[has_factor]
        molar = (num_s.loc[idx_f] * factors.loc[idx_f].astype(float)).clip(
            lower=1e-30
        )
        out.loc[idx_f] = -np.log10(molar)
    else:
      out.loc[valid_num] = num_s.loc[valid_num]
  rem_mask = out.isna() & (~finite_mask) & series.notna()
  if rem_mask.any():
    out.loc[rem_mask] = pd.to_numeric(
        series.loc[rem_mask].apply(_parse_discrete_potency_tier),
        errors="coerce",
    )
  return out


def find_activity_column(
    df: pd.DataFrame, user_col: Optional[str] = None
) -> Optional[str]:
  """Auto-detects the primary bioactivity or potency column."""
  if user_col and user_col in df.columns:
    return user_col

  ordered_cols = sorted(
      df.columns, key=lambda c: bool(_ANON_PATENT_COL_RE.match(str(c)))
  )
  numeric_patterns = (
      # Priority 1: Prefix negative-log molar potency/affinity
      rf"^(?:{_LOG_POTENCY_PAT}|potency)(?:_|$|(?:{_NORM_UNIT_PAT})(?:_|$))",
      # Priority 2: Prefix raw concentration / standard activity metrics
      (
          rf"^(?:{_RAW_CONC_PAT})(?:_|$|(?:{_NORM_UNIT_PAT})(?:_|$))"
          r"|^(?:activity_value|standard_value)$"
      ),
      # Priority 3: Substring negative-log molar potency columns
      rf"(?:^|_)(?:{_LOG_POTENCY_PAT})(?:_|$|(?:{_NORM_UNIT_PAT})(?:_|$))",
      # Priority 4: Substring raw concentration potency columns
      rf"(?:^|_)(?:{_RAW_CONC_PAT})(?:_|$|(?:{_NORM_UNIT_PAT})(?:_|$))",
  )
  for pat in numeric_patterns:
    for col in ordered_cols:
      if (
          not _is_potency_meta_col(col)
          and not _is_ratio_or_percent_col(col)
          and re.search(pat, _norm_col(col))
      ):
        num_s = pd.to_numeric(df[col], errors="coerce")
        if (num_s.notna() & np.isfinite(num_s)).any():
          return col

  # Priority 5: Discrete / binned / inequality potency columns
  for col in ordered_cols:
    if (
        _is_potency_col_name(col)
        and df[col].apply(_parse_discrete_potency_tier).notna().any()
    ):
      return col

  return None


def _as_float(value: Any, default: float = 0.0) -> float:
  """Coerces a cell to float, mapping None, NaN and unparseable values to default.

  The `value or default` idiom is unsafe here: NaN is truthy, so it survives
  the `or` and then raises in float()/int().

  Args:
    value: Raw cell value from the input table.
    default: Value to substitute when the cell is missing or unparseable.

  Returns:
    A finite float.
  """
  try:
    result = float(value)
  except (TypeError, ValueError):
    return default
  return default if math.isnan(result) or math.isinf(result) else result


def _as_int(value: Any, default: int = 0) -> int:
  """Coerces a cell to int, mapping None, NaN and unparseable values to default.

  Args:
    value: Raw cell value from the input table.
    default: Value to substitute when the cell is missing or unparseable.

  Returns:
    An int.
  """
  fval = _as_opt_float(value)
  return int(round(fval)) if fval is not None else default


def _as_opt_float(value: Any) -> Optional[float]:
  """Coerces a cell to float, mapping None, NaN and unparseable values to None (never 0.0)."""
  if value is None or pd.isna(value):
    return None
  try:
    result = float(value)
  except (TypeError, ValueError):
    return None
  return None if math.isnan(result) or math.isinf(result) else result


def _as_opt_int(value: Any) -> Optional[int]:
  """Coerces a cell to int, mapping None, NaN and unparseable values to None (never 0)."""
  fval = _as_opt_float(value)
  return int(round(fval)) if fval is not None else None


def _as_clean_str(value: Any) -> str:
  """Returns stripped cell text, mapping None and NaN to empty string."""
  if value is None or pd.isna(value):
    return ""
  s = str(value).strip()
  return "" if s.lower() == "nan" else s


def create_molecular_grid(
    df: pd.DataFrame,
    smiles_col: Optional[str] = None,
    id_col: Optional[str] = None,
    activity_col: Optional[str] = None,
    filter_active: Optional[bool] = None,
    output_html: str = "molecular_grid.html",
    title: str = "Interactive Molecular Grid & SAR Explorer",
    align_core: Optional[str] = None,
    max_cards: int = 250,
) -> str:
  """Renders a 100% self-contained, CSP-compliant chemical grid HTML with SAR scatter plot and property badges."""
  df = df.copy()

  # Identify SMILES column first and build ROMol so structureless rows cannot
  # crowd out valid structures or become MW=0 / cLogP=0 ghost cards.
  if not smiles_col:
    priority_smiles = (
        "product_smiles",
        "canonical_smiles",
        "smiles",
        "clean_smiles",
        "mol_smiles",
        "reactant_smiles",
    )
    smi_col_by_norm: dict[str, str] = {}
    for c in df.columns:
      smi_col_by_norm.setdefault(_norm_col(c), c)
    matched_smi = [
        smi_col_by_norm[p] for p in priority_smiles if p in smi_col_by_norm
    ]
    if matched_smi:
      smiles_col = matched_smi[0]
    else:
      tok_smi = [
          c
          for c in df.columns
          if re.search(r"(?:^|[\W_])smiles(?:$|[\W_])", c, re.I)
      ]
      if tok_smi:
        smiles_col = tok_smi[0]
      elif "ROMol" in df.columns:
        df["SMILES"] = df["ROMol"].apply(
            lambda m: Chem.MolToSmiles(m) if m is not None else ""
        )
        smiles_col = "SMILES"
      else:
        smi_sub = [c for c in df.columns if "smi" in c.lower()]
        if smi_sub:
          smiles_col = smi_sub[0]
        elif "structure" in smi_col_by_norm:
          smiles_col = smi_col_by_norm["structure"]
        elif len(df.columns) > 0:
          smiles_col = df.columns[0]
        else:
          smiles_col = "SMILES"
        if smiles_col not in df.columns:
          df[smiles_col] = pd.Series(dtype=object)

  if "ROMol" not in df.columns:
    df["ROMol"] = df[smiles_col].apply(_parse_grid_smiles)

  valid_mol_mask = df["ROMol"].notna()
  if not df.empty and not valid_mol_mask.any():
    raise ValueError(
        f"No valid chemical structures parsed from column '{smiles_col}'."
    )
  dropped_no_mol = int((~valid_mol_mask).sum())
  if dropped_no_mol > 0:
    print(
        f"Omitted {dropped_no_mol} structureless (missing/NaN SMILES) rows"
        " from the 2D molecular grid viewer."
    )
    df = df[valid_mol_mask].copy()

  total_records = len(df)

  # Detect primary activity/potency column (numeric or discrete/binned)
  act_col = find_activity_column(df, activity_col)
  act_type_col, act_units_col, act_range_cols = (
      _resolve_act_companions(act_col, df.columns) if act_col else ("", "", ())
  )
  should_filter = False
  if act_col:
    units_s = df[act_units_col] if act_units_col in df.columns else None
    act_scores = _series_potency_scores(df[act_col], act_col, units_s)
    df["_potency_sort_score"] = act_scores
    num_act = pd.to_numeric(df[act_col], errors="coerce")
    finite_act = num_act.notna() & np.isfinite(num_act)
    valid_num_act = (
        (finite_act & (num_act > 0))
        if _is_raw_conc_column(act_col)
        else finite_act
    )
    has_activity = act_scores.notna() | valid_num_act
    num_active = int(has_activity.sum())

    # Default behavior: If activity data exists and dataset contains uncharacterized/untested rows,
    # filter to show only active compounds unless filter_active is explicitly False.
    should_filter = (filter_active is True) or (
        filter_active is None and num_active > 0 and num_active < len(df)
    )
    if should_filter and num_active > 0:
      print(
          f"Defaulting to active compounds: filtered from {len(df)} to"
          f" {num_active} compounds with valid activity in '{act_col}'."
      )
      df = df[has_activity].copy()
    if num_active > 0:
      df = df.sort_values(
          by="_potency_sort_score",
          ascending=False,
          na_position="last",
          kind="stable",
      )

  # Automatically identify ID/Name column
  if not id_col:
    priority_id_cols = (
        "compound_id",
        "mol_id",
        "id",
        "name",
        "compound_name",
        "compound",
    )
    id_col_by_norm: dict[str, str] = {}
    for c in df.columns:
      if c != smiles_col:
        id_col_by_norm.setdefault(_norm_col(c), c)
    matched_priority = [
        id_col_by_norm[p] for p in priority_id_cols if p in id_col_by_norm
    ]
    if matched_priority:
      id_col = matched_priority[0]
    else:
      id_candidates = [
          c
          for c in df.columns
          if c != smiles_col
          and (
              "id" in c.lower()
              or "name" in c.lower()
              or "compound" in c.lower()
          )
      ]
      title_candidates = [
          c for c in df.columns if c != smiles_col and "title" in c.lower()
      ]
      if id_candidates:
        id_col = id_candidates[0]
      elif title_candidates:
        id_col = title_candidates[0]
      else:
        id_col = smiles_col

  # Parse align core if given
  core_mol = None
  if align_core:
    core_mol = _parse_grid_smiles(align_core) or Chem.MolFromSmarts(align_core)
    if core_mol is not None:
      compute_clean_2d_coords(core_mol)

  # Filter catalog for PAINS / Brenk
  pains_catalog = get_pains_catalog()

  # Compute or backfill missing/NaN baseline properties in C++ (never fallback to 0.0 for missing mols)
  def _fill_prop(col_name: str, fn, decimals: Optional[int] = None) -> None:
    if col_name not in df.columns:
      df[col_name] = np.nan
    df[col_name] = pd.to_numeric(df[col_name], errors="coerce")
    miss = df[col_name].isna() & df["ROMol"].notna()
    if miss.any():
      df.loc[miss, col_name] = df.loc[miss, "ROMol"].apply(fn)
    if decimals is not None:
      df[col_name] = df[col_name].round(decimals)

  _fill_prop("MW", lambda m: round(Descriptors.MolWt(m), 1), decimals=1)
  _fill_prop("cLogP", lambda m: round(Crippen.MolLogP(m), 2), decimals=2)
  _fill_prop(
      "TPSA", lambda m: round(rdMolDescriptors.CalcTPSA(m), 1), decimals=1
  )
  _fill_prop("HBD", lambda m: int(Lipinski.NumHDonors(m)))
  _fill_prop("HBA", lambda m: int(Lipinski.NumHAcceptors(m)))
  _fill_prop("RotB", lambda m: int(Lipinski.NumRotatableBonds(m)))
  _fill_prop(
      "Fsp3",
      lambda m: round(rdMolDescriptors.CalcFractionCSP3(m), 2),
      decimals=2,
  )
  _fill_prop("QED", _safe_qed, decimals=2)
  _fill_prop("CNS_MPO", calculate_cns_mpo, decimals=2)

  # Compute Ro5 Violations and PAINS Status if missing
  def check_ro5(row):
    mw_v = _as_opt_float(row.get("MW"))
    lp_v = _as_opt_float(row.get("cLogP"))
    hbd_v = _as_opt_int(row.get("HBD"))
    hba_v = _as_opt_int(row.get("HBA"))
    if mw_v is None and lp_v is None:
      return None
    v = 0
    if mw_v is not None and mw_v > 500:
      v += 1
    if lp_v is not None and lp_v > 5.0:
      v += 1
    if hbd_v is not None and hbd_v > 5:
      v += 1
    if hba_v is not None and hba_v > 10:
      v += 1
    return v

  if "Ro5_Violations" not in df.columns:
    df["Ro5_Violations"] = (
        pd.Series(dtype=object) if df.empty else df.apply(check_ro5, axis=1)
    )
  else:
    miss_ro5 = df["Ro5_Violations"].isna()
    if miss_ro5.any():
      df.loc[miss_ro5, "Ro5_Violations"] = df.loc[miss_ro5].apply(
          check_ro5, axis=1
      )

  def check_pains(m):
    if not m:
      return "Clean"
    matches = pains_catalog.GetMatches(m)
    return (
        "; ".join([e.GetDescription() for e in matches]) if matches else "Clean"
    )

  if "PAINS_Status" not in df.columns:
    df["PAINS_Status"] = df["ROMol"].apply(check_pains)

  # Display safety limit: Stay well under the 20MB IDE webview display limit
  display_df = df
  if max_cards and max_cards > 0 and len(df) > max_cards:
    display_df = df.iloc[:max_cards]
    print(
        f"Showing top {len(display_df)} of {len(df)} compounds; pass"
        " --max_cards to increase."
    )

  # Identify numeric and discrete-potency properties for dynamic SAR plotting
  numeric_cols = []
  discrete_potency_cols = []
  for c in display_df.columns:
    if c in ("ROMol", "svg_data", "_potency_sort_score"):
      continue
    num_s = pd.to_numeric(display_df[c], errors="coerce")
    if (num_s.notna() & np.isfinite(num_s)).any():
      numeric_cols.append(c)
    elif c == act_col or _is_potency_col_name(c):
      tier_s = display_df[c].apply(_parse_discrete_potency_tier)
      if tier_s.notna().any():
        discrete_potency_cols.append(c)

  # Only embed columns that carry a value for at least one displayed compound.
  candidate_cols = [
      c
      for c in display_df.columns
      if c not in ("ROMol", "svg_data", "_potency_sort_score")
  ]
  embed_cols = [c for c in candidate_cols if display_df[c].notna().any()]
  dropped_cols = len(candidate_cols) - len(embed_cols)
  if dropped_cols > 0 and not display_df.empty:
    print(
        f"Omitted {dropped_cols} columns that were empty for all"
        f" {len(display_df)} displayed compounds."
    )

  modal_skip_keys = [
      "ROMol",
      "svg_data",
      "MW",
      "cLogP",
      "TPSA",
      "HBD",
      "HBA",
      "RotB",
      "Fsp3",
      "QED",
      "CNS_MPO",
      "PAINS_Status",
      "Ro5_Violations",
  ]
  for dynamic_col in (smiles_col, id_col):
    if dynamic_col and dynamic_col not in modal_skip_keys:
      modal_skip_keys.append(dynamic_col)

  compounds_data = []
  for idx, (_, row) in enumerate(display_df.iterrows()):
    mol = row.get("ROMol")
    smi = str(row.get(smiles_col, ""))
    name = _as_clean_str(row.get(id_col)) or f"Compound_{idx+1}"
    svg = generate_svg_string(mol, width=250, height=180, core_mol=core_mol)

    meta_dict = {}
    props_dict = {}
    for c in embed_cols:
      val = row.get(c)
      if pd.isna(val):
        continue
      text = str(val).strip()
      if not text:
        continue
      meta_dict[c] = text
      if c in numeric_cols:
        fval = _as_opt_float(val)
        if fval is not None:
          props_dict[c] = round(fval, 4)
        elif _is_raw_conc_column(c):
          c_units_col = (
              act_units_col
              if c == act_col
              else _resolve_act_companions(c, display_df.columns)[1]
          )
          c_unit_hint = (
              _as_clean_str(row.get(c_units_col))
              if c_units_col in display_df.columns
              else ""
          )
          col_molar = _unit_molar_factor(c, c_unit_hint)
          if col_molar is not None:
            v_m, op_c = _parse_conc_string_molar(text, default_molar=col_molar)
            if v_m is not None and not any(
                q in op_c for q in ("<", "≤", ">", "≥")
            ):
              props_dict[c] = round(v_m / col_molar, 6)
        else:
          m_approx = re.match(r"^[~=]+\s*([+-]?\d+(?:\.\d+)?)\s*$", text)
          if m_approx:
            props_dict[c] = round(float(m_approx.group(1)), 4)
      elif c in discrete_potency_cols:
        tval = _parse_discrete_potency_tier(val)
        if tval is not None:
          props_dict[c] = tval

    # Determine primary potency score and human-readable badge display
    potency_score: Optional[float] = None
    potency_display: Optional[str] = None
    if act_col and pd.notna(row.get(act_col)):
      raw_act = _as_clean_str(row.get(act_col))
      if raw_act:
        f_act = _as_opt_float(raw_act)
        if f_act is not None:
          if _is_raw_conc_column(act_col):
            if f_act > 0:
              unit_s = (
                  _as_clean_str(row.get(act_units_col))
                  if act_units_col in display_df.columns
                  else ""
              )
              factor = _unit_molar_factor(act_col, unit_s)
              if factor is not None:
                potency_score = round(
                    -math.log10(max(f_act * factor, 1e-30)), 3
                )
              type_s = (
                  _as_clean_str(row.get(act_type_col))
                  if act_type_col in display_df.columns
                  else ""
              )
              lbl_s = _format_col_label(type_s or act_col)
              has_unit_in_lbl = bool(
                  unit_s
                  and re.search(
                      rf"(?:^|[\s(\[-]){re.escape(unit_s)}(?:$|[\s)\]-])",
                      lbl_s,
                      flags=re.I,
                  )
              )
              suffix_s = (
                  f" {unit_s}" if (unit_s and not has_unit_in_lbl) else ""
              )
              potency_display = f"{lbl_s}: {raw_act}{suffix_s}".strip()
          else:
            potency_score = round(f_act, 3)
            potency_display = f"{_format_col_label(act_col)}: {round(f_act, 2)}"
        else:
          potency_score = _parse_discrete_potency_tier(raw_act)
          if potency_score is not None:
            rng_s = next(
                (
                    _as_clean_str(row.get(rc))
                    for rc in act_range_cols
                    if rc in display_df.columns and _as_clean_str(row.get(rc))
                ),
                "",
            )
            val_with_rng = (
                f"{raw_act} ({rng_s})"
                if (rng_s and rng_s not in raw_act)
                else raw_act
            )
            type_s = (
                _as_clean_str(row.get(act_type_col))
                if act_type_col in display_df.columns
                else ""
            )
            lbl_s = _format_col_label(type_s or act_col)
            potency_display = f"{lbl_s}: {val_with_rng}".strip()

    pains_raw = row.get("PAINS_Status")
    pains_val = (
        str(pains_raw).strip()
        if pd.notna(pains_raw) and str(pains_raw).strip().lower() != "nan"
        else "Clean"
    )
    ro5_v = _as_opt_int(row.get("Ro5_Violations"))
    cns_mpo_val = _as_opt_float(row.get("CNS_MPO"))

    for k_prop, v_prop in [
        ("MW", _as_opt_float(row.get("MW"))),
        ("cLogP", _as_opt_float(row.get("cLogP"))),
        ("TPSA", _as_opt_float(row.get("TPSA"))),
        ("HBD", _as_opt_int(row.get("HBD"))),
        ("HBA", _as_opt_int(row.get("HBA"))),
        ("RotB", _as_opt_int(row.get("RotB"))),
        ("Fsp3", _as_opt_float(row.get("Fsp3"))),
        ("QED", _as_opt_float(row.get("QED"))),
        (
            "CNS_MPO",
            cns_mpo_val
            if (cns_mpo_val is not None and cns_mpo_val >= 0)
            else None,
        ),
        ("Ro5_Violations", ro5_v),
    ]:
      if v_prop is not None:
        props_dict.setdefault(k_prop, v_prop)

    compounds_data.append({
        "idx": idx,
        "id": name,
        "smiles": smi,
        "mw": _as_opt_float(row.get("MW")),
        "clogp": _as_opt_float(row.get("cLogP")),
        "tpsa": _as_opt_float(row.get("TPSA")),
        "hbd": _as_opt_int(row.get("HBD")),
        "hba": _as_opt_int(row.get("HBA")),
        "rotb": _as_opt_int(row.get("RotB")),
        "fsp3": _as_opt_float(row.get("Fsp3")),
        "qed": _as_opt_float(row.get("QED")),
        "cns_mpo": cns_mpo_val,
        "ro5_violations": ro5_v,
        "pains": pains_val,
        "potency_score": potency_score,
        "potency_display": potency_display,
        "svg": svg,
        "props": props_dict,
        "meta": meta_dict,
    })

  def _js_json(obj: Any) -> str:
    return json.dumps(obj, allow_nan=False).replace("</", "<\\/")

  json_bytes = json.dumps(compounds_data, allow_nan=False).encode("utf-8")
  gzip_b64_payload = base64.b64encode(
      gzip.compress(json_bytes, compresslevel=9)
  ).decode("ascii")

  header_subtitle = f"Showing {len(compounds_data)} compounds"
  if total_records > len(compounds_data):
    if act_col and should_filter:
      header_subtitle = (
          f"Showing {len(compounds_data)} active compounds with measured"
          f" {act_col} (from {total_records} total in dataset)"
      )
    else:
      header_subtitle = (
          f"Showing top {len(compounds_data)} of {total_records} compounds"
          f" (Full {total_records} in dataset)"
      )

  active_col_hint_js = _js_json(act_col or "")
  smiles_col_hint_js = _js_json(smiles_col or "")

  # Self-contained HTML template with Grid View, SAR Scatter View, and Property Badges
  html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{html.escape(title)}</title>
    <style>
        :where(:root) {{
            /* Generative UI Design System Core Semantic Tokens (0-specificity fallback so host :root tokens win) */
            --background: #f8fafc;
            --card: #ffffff;
            --border: #e2e8f0;
            --card-border: #e2e8f0;
            --foreground: #0f172a;
            --muted-foreground: #64748b;
            --muted: #f1f5f9;
            --primary: #4f46e5;
            --primary-foreground: #ffffff;
            --accent: #6366f1;
            --accent-hover: #4338ca;
            
            /* Status Badges */
            --badge-success-bg: #ecfdf5;
            --badge-success-txt: #059669;
            --badge-warning-bg: #fffbeb;
            --badge-warning-txt: #d97706;
            --badge-danger-bg: #fff1f2;
            --badge-danger-txt: #e11d48;
            --badge-info-bg: #e0e7ff;
            --badge-info-txt: #4338ca;

            /* Compatibility Aliases */
            --bg-page: var(--background);
            --bg-card: var(--card);
            --border-card: var(--card-border, var(--border));
            --text-primary: var(--foreground);
            --text-secondary: var(--muted-foreground);
            --badge-green-bg: var(--badge-success-bg);
            --badge-green-txt: var(--badge-success-txt);
            --badge-amber-bg: var(--badge-warning-bg);
            --badge-amber-txt: var(--badge-warning-txt);
            --badge-red-bg: var(--badge-danger-bg);
            --badge-red-txt: var(--badge-danger-txt);
            --badge-blue-bg: var(--badge-info-bg);
            --badge-blue-txt: var(--badge-info-txt);
            --font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Inter", "Google Sans", Helvetica, Arial, sans-serif;
        }}
        :where(html.dark) {{
            --background: #0f172a;
            --card: #1e293b;
            --border: #334155;
            --card-border: #334155;
            --foreground: #f8fafc;
            --muted-foreground: #94a3b8;
            --muted: #1e293b;
            --primary: #6366f1;
            --primary-foreground: #ffffff;
            --accent: #818cf8;
            --accent-hover: #a5b4fc;

            --badge-success-bg: rgba(16, 185, 129, 0.15);
            --badge-success-txt: #34d399;
            --badge-warning-bg: rgba(245, 158, 11, 0.15);
            --badge-warning-txt: #fbbf24;
            --badge-danger-bg: rgba(244, 63, 94, 0.15);
            --badge-danger-txt: #fb7185;
            --badge-info-bg: rgba(99, 102, 241, 0.15);
            --badge-info-txt: #a5b4fc;
        }}
        @media (prefers-color-scheme: dark) {{
            :where(:root:not(.light)) {{
                --background: #0f172a;
                --card: #1e293b;
                --border: #334155;
                --card-border: #334155;
                --foreground: #f8fafc;
                --muted-foreground: #94a3b8;
                --muted: #1e293b;
                --primary: #6366f1;
                --primary-foreground: #ffffff;
                --accent: #818cf8;
                --accent-hover: #a5b4fc;

                --badge-success-bg: rgba(16, 185, 129, 0.15);
                --badge-success-txt: #34d399;
                --badge-warning-bg: rgba(245, 158, 11, 0.15);
                --badge-warning-txt: #fbbf24;
                --badge-danger-bg: rgba(244, 63, 94, 0.15);
                --badge-danger-txt: #fb7185;
                --badge-info-bg: rgba(99, 102, 241, 0.15);
                --badge-info-txt: #a5b4fc;
            }}
        }}
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: var(--font-family);
            background: var(--bg-page);
            color: var(--text-primary);
            padding: 16px;
        }}
        .header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 12px;
            margin-bottom: 16px;
            padding-bottom: 12px;
            border-bottom: 1px solid var(--border-card);
        }}
        .title-group h1 {{
            font-size: 18px;
            font-weight: 600;
            color: var(--text-primary);
        }}
        .title-group .count {{
            font-size: 12px;
            color: var(--text-secondary);
        }}
        .controls {{
            display: flex;
            gap: 8px;
            align-items: center;
            flex-wrap: wrap;
        }}
        .view-toggle {{
            display: inline-flex;
            background: var(--bg-card);
            border: 1px solid var(--border-card);
            border-radius: 8px;
            padding: 3px;
            gap: 2px;
        }}
        .toggle-btn {{
            background: none;
            border: none;
            padding: 5px 12px;
            border-radius: 6px;
            font-size: 12px;
            font-weight: 500;
            color: var(--text-secondary);
            cursor: pointer;
            transition: background 0.15s, color 0.15s;
        }}
        .toggle-btn.active {{
            background: var(--primary);
            color: #ffffff;
        }}
        input.search {{
            padding: 6px 12px;
            border: 1px solid var(--border-card);
            border-radius: 8px;
            background: var(--bg-card);
            color: var(--text-primary);
            font-size: 13px;
            outline: none;
            width: min(260px, 100%);
            transition: border-color 0.15s;
        }}
        input.search:focus {{
            border-color: var(--primary);
        }}
        select.sort {{
            padding: 6px 10px;
            border: 1px solid var(--border-card);
            border-radius: 8px;
            background: var(--bg-card);
            color: var(--text-primary);
            font-size: 13px;
            outline: none;
            cursor: pointer;
        }}
        .action-btn {{
            padding: 6px 12px;
            border: 1px solid var(--border-card);
            border-radius: 8px;
            background: var(--bg-card);
            color: var(--text-primary);
            font-size: 12px;
            font-weight: 500;
            cursor: pointer;
            transition: background 0.15s, border-color 0.15s;
        }}
        .action-btn:hover {{
            background: rgba(127,127,127,0.1);
            border-color: var(--primary);
        }}

        /* Grid View */
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(260px, 1fr));
            gap: 14px;
        }}
        .card {{
            background: var(--bg-card);
            border: 1px solid var(--border-card);
            border-radius: 12px;
            padding: 12px;
            display: flex;
            flex-direction: column;
            box-shadow: 0 1px 3px 0 rgba(0,0,0,0.05), 0 1px 2px -1px rgba(0,0,0,0.05);
            transition: transform 0.15s ease, box-shadow 0.15s ease, border-color 0.15s ease;
            min-width: 0;
            overflow: hidden;
            contain: content;
            cursor: pointer;
            position: relative;
        }}
        .card:hover {{
            box-shadow: 0 10px 15px -3px rgba(0,0,0,0.08), 0 4px 6px -4px rgba(0,0,0,0.05);
            transform: translateY(-2px);
            border-color: var(--primary);
        }}
        .card-header {{
            font-size: 13px;
            font-weight: 600;
            color: var(--primary);
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            margin-bottom: 6px;
            width: 100%;
            min-width: 0;
        }}
        .badge-row {{
            display: flex;
            gap: 4px;
            flex-wrap: wrap;
            margin-bottom: 6px;
        }}
        .badge {{
            font-size: 10px;
            font-weight: 500;
            padding: 2px 8px;
            border-radius: 9999px;
            letter-spacing: 0.2px;
            display: inline-block;
            vertical-align: middle;
            max-width: 100%;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }}
        .badge-green {{ background: var(--badge-green-bg); color: var(--badge-green-txt); }}
        .badge-amber {{ background: var(--badge-amber-bg); color: var(--badge-amber-txt); }}
        .badge-red {{ background: var(--badge-red-bg); color: var(--badge-red-txt); }}
        .badge-blue {{ background: var(--badge-blue-bg); color: var(--badge-blue-txt); }}

        .svg-container {{
            background: #ffffff !important;
            border: 1px solid rgba(127,127,127,0.18);
            border-radius: 6px;
            display: flex;
            align-items: center;
            justify-content: center;
            height: 165px;
            margin-bottom: 8px;
            overflow: hidden;
            width: 100%;
            box-sizing: border-box;
            transition: border-color 0.15s, box-shadow 0.15s;
        }}
        .card:hover .svg-container {{
            border-color: var(--primary);
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
        }}
        .svg-container svg {{
            width: 100% !important;
            height: 100% !important;
            max-height: 155px;
            display: block;
        }}
        .props-grid {{
            display: grid;
            grid-template-columns: repeat(3, minmax(0, 1fr));
            gap: 6px;
            font-size: 11px;
            color: var(--text-secondary);
            margin-bottom: 6px;
            background: rgba(127,127,127,0.04);
            padding: 6px 8px;
            border-radius: 6px;
        }}
        .prop-item {{
            display: flex;
            flex-direction: column;
            align-items: flex-start;
            min-width: 0;
            line-height: 1.25;
        }}
        .prop-label {{
            font-size: 9.5px;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.03em;
            color: var(--text-secondary);
        }}
        .prop-val {{
            font-size: 11.5px;
            font-weight: 600;
            color: var(--text-primary);
            font-variant-numeric: tabular-nums;
            white-space: nowrap;
        }}
        .smiles-wrapper {{
            display: flex;
            align-items: center;
            width: 100%;
            min-width: 0;
            overflow: hidden;
        }}
        .smiles-row {{
            font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
            font-size: 10.5px;
            color: var(--text-secondary);
            background: rgba(127,127,127,0.07);
            padding: 3px 6px;
            border-radius: 6px;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            cursor: pointer;
            width: 100%;
            box-sizing: border-box;
            border: 1px solid transparent;
            transition: border-color 0.15s, background 0.15s, color 0.15s;
        }}
        .smiles-row:hover {{
            background: rgba(79,70,229,0.1);
            color: var(--primary);
            border-color: var(--primary);
        }}
        .meta-list {{
            margin-top: 6px;
            font-size: 10.5px;
            color: var(--text-secondary);
            border-top: 1px dashed rgba(127,127,127,0.15);
            padding-top: 5px;
            display: flex;
            flex-direction: column;
            gap: 2px;
            overflow: hidden;
            width: 100%;
        }}
        .meta-row {{
            display: flex;
            align-items: baseline;
            justify-content: flex-start;
            gap: 6px;
            min-width: 0;
            width: 100%;
            line-height: 1.3;
        }}
        .meta-key {{
            font-weight: 600;
            color: var(--text-primary);
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            max-width: 68%;
            flex-shrink: 1;
        }}
        .meta-val {{
            color: var(--text-secondary);
            font-variant-numeric: tabular-nums;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            flex-shrink: 0;
        }}
        .zoom-hint {{
            position: absolute;
            top: 8px;
            right: 8px;
            font-size: 10px;
            background: rgba(0,0,0,0.65);
            color: #ffffff;
            padding: 2px 7px;
            border-radius: 10px;
            opacity: 0;
            transition: opacity 0.15s;
            pointer-events: none;
        }}
        .card:hover .zoom-hint {{
            opacity: 1;
        }}

        /* Scatter Plot View */
        .scatter-view {{
            background: var(--bg-card);
            border: 1px solid var(--border-card);
            border-radius: 12px;
            padding: 14px 16px;
            display: flex;
            flex-direction: column;
            align-items: center;
            position: relative;
        }}
        .scatter-svg {{
            width: 100%;
            max-width: 680px;
            max-height: min(340px, 40vh);
            height: auto;
            display: block;
            margin: 0 auto;
        }}
        .scatter-dot {{
            cursor: pointer;
            transition: r 0.15s, stroke-width 0.15s;
        }}
        .scatter-dot:hover {{
            r: 9;
            stroke: var(--primary);
            stroke-width: 3;
        }}
        .scatter-tooltip {{
            position: absolute;
            background: rgba(15,23,42,0.94);
            color: #ffffff;
            padding: 8px 12px;
            border-radius: 8px;
            font-size: 11.5px;
            pointer-events: none;
            opacity: 0;
            transition: opacity 0.15s;
            z-index: 1000;
            line-height: 1.4;
            max-width: 240px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.2);
        }}

        /* Instant Modal Pop-Out Zoom */
        .modal-overlay {{
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: rgba(0, 0, 0, 0.75);
            display: none;
            justify-content: center;
            align-items: center;
            z-index: 99999999;
            padding: 16px;
            box-sizing: border-box;
            overflow-y: auto;
            backdrop-filter: blur(4px);
        }}
        .modal-overlay.active {{
            display: flex !important;
        }}
        .modal-card {{
            background: var(--bg-card);
            border: 1px solid var(--border-card);
            border-radius: 16px;
            max-width: 860px;
            width: 95%;
            max-height: 90vh;
            display: flex;
            flex-direction: column;
            overflow: hidden;
            box-shadow: 0 20px 25px -5px rgba(0,0,0,0.2), 0 8px 10px -6px rgba(0,0,0,0.2);
            position: relative;
            z-index: 100000000;
            margin: auto;
        }}
        .modal-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 14px 20px;
            border-bottom: 1px solid var(--border-card);
            background: rgba(127,127,127,0.04);
        }}
        .modal-header h2 {{
            font-size: 16px;
            font-weight: 600;
            color: var(--primary);
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            max-width: 680px;
        }}
        .modal-close {{
            background: none;
            border: none;
            font-size: 24px;
            line-height: 1;
            color: var(--text-secondary);
            cursor: pointer;
            padding: 4px 8px;
            border-radius: 6px;
        }}
        .modal-close:hover {{
            background: rgba(127,127,127,0.15);
            color: var(--text-primary);
        }}
        .modal-body {{
            padding: 20px;
            overflow-y: auto;
            display: grid;
            grid-template-columns: 1.15fr 1fr;
            gap: 20px;
            align-items: start;
            max-height: calc(90vh - 130px);
        }}
        @media (max-width: 720px) {{
            .modal-body {{ grid-template-columns: 1fr; }}
        }}
        .modal-svg-view {{
            background: #ffffff !important;
            border: 1px solid rgba(127,127,127,0.22);
            border-radius: 12px;
            padding: 16px;
            display: flex;
            align-items: center;
            justify-content: center;
            min-height: 290px;
            width: 100%;
            overflow: hidden;
            box-sizing: border-box;
            box-shadow: inset 0 1px 3px rgba(0,0,0,0.05);
        }}
        .modal-svg-view svg {{
            width: 100% !important;
            height: auto !important;
            max-height: 280px;
            display: block;
        }}
        .modal-details {{
            display: flex;
            flex-direction: column;
            gap: 12px;
            width: 100%;
            min-width: 0;
        }}
        .modal-details > div {{
            min-width: 0;
        }}
        .modal-smiles-box {{
            background: rgba(127,127,127,0.07);
            border: 1px solid rgba(127,127,127,0.18);
            border-radius: 8px;
            padding: 8px 10px;
            font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
            font-size: 11px;
            word-break: break-all;
            position: relative;
            color: var(--text-primary);
        }}
        .copy-btn {{
            margin-top: 6px;
            background: var(--primary);
            color: #ffffff;
            border: none;
            border-radius: 6px;
            padding: 5px 12px;
            font-size: 11.5px;
            font-weight: 500;
            cursor: pointer;
            display: inline-flex;
            align-items: center;
            gap: 4px;
            transition: background 0.15s;
        }}
        .copy-btn:hover {{
            background: var(--accent-hover);
        }}
        .modal-props-table {{
            width: 100%;
            table-layout: fixed;
            border-collapse: collapse;
            font-size: 12px;
        }}
        .modal-props-table tr {{
            border-bottom: 1px solid rgba(127,127,127,0.1);
        }}
        .modal-props-table td {{
            padding: 6px 4px;
            vertical-align: top;
        }}
        .modal-props-table td.label {{
            color: var(--text-secondary);
            font-weight: 500;
            width: 46%;
            padding-right: 10px;
            overflow-wrap: anywhere;
            word-break: break-word;
        }}
        .modal-props-table td.val {{
            width: 54%;
            font-weight: 600;
            color: var(--text-primary);
            font-variant-numeric: tabular-nums;
            overflow-wrap: anywhere;
        }}
        .modal-footer {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 12px 20px;
            border-top: 1px solid var(--border-card);
            background: rgba(127,127,127,0.04);
        }}
        .nav-btn {{
            background: var(--bg-card);
            border: 1px solid var(--border-card);
            border-radius: 6px;
            padding: 6px 14px;
            font-size: 12px;
            font-weight: 500;
            color: var(--text-primary);
            cursor: pointer;
            transition: background 0.15s;
        }}
        .nav-btn:hover {{
            background: rgba(127,127,127,0.12);
        }}
        .nav-btn:disabled {{
            opacity: 0.35;
            cursor: not-allowed;
        }}

        .copy-toast {{
            position: fixed;
            bottom: 20px;
            right: 20px;
            background: #202124;
            color: #ffffff;
            padding: 8px 16px;
            border-radius: 6px;
            font-size: 12px;
            opacity: 0;
            transition: opacity 0.2s;
            pointer-events: none;
            z-index: 2000000;
        }}
    </style>
</head>
<body>
    <div class="header">
        <div class="title-group">
            <h1>{html.escape(title)}</h1>
            <span class="count" id="itemCount">{html.escape(header_subtitle)}</span>
        </div>
        <div class="controls">
            <div class="view-toggle">
                <button class="toggle-btn active" id="btnGridView" onclick="switchView('grid')">Grid View</button>
                <button class="toggle-btn" id="btnScatterView" onclick="switchView('scatter')">SAR Scatter Plot</button>
            </div>
            <input type="text" class="search" id="searchInput" placeholder="Filter: name, SMILES, mw&lt;450, pki&gt;9..." title="Smart filter: multi-term AND, OR, -exclude, phrases (&quot;ro5 pass&quot;), fields &amp; comparisons (mw&lt;450, clogp&lt;=3, mpo&gt;=4, pki&gt;9, ic50&lt;10, ro5=0, mw:300..450, pains:clean, series:macrocyclic)" oninput="filterGrid()">
            <select class="sort" id="sortSelect" onchange="sortGrid()">
                <option value="default">Default Order</option>
                <option value="potency_desc">Potency (High &rarr; Low)</option>
                <option value="potency_asc">Potency (Low &rarr; High)</option>
                <option value="mw_asc">MW (Low &rarr; High)</option>
                <option value="mw_desc">MW (High &rarr; Low)</option>
                <option value="clogp_asc">cLogP (Low &rarr; High)</option>
                <option value="clogp_desc">cLogP (High &rarr; Low)</option>
                <option value="cns_mpo_desc">CNS MPO (High &rarr; Low)</option>
                <option value="qed_desc">QED (High &rarr; Low)</option>
                <option value="tpsa_asc">TPSA (Low &rarr; High)</option>
            </select>
            <button class="action-btn" onclick="exportFilteredCSV()">Export CSV</button>
        </div>
    </div>

    <!-- 1. Grid View -->
    <div class="grid" id="gridContainer"></div>

    <!-- 2. SAR Scatter Plot View -->
    <div class="scatter-view" id="scatterContainer" style="display: none;">
        <div class="scatter-toolbar" style="width: 100%; display: flex; flex-wrap: wrap; justify-content: space-between; align-items: center; gap: 12px; margin-bottom: 14px; padding: 8px 14px; background: rgba(127,127,127,0.06); border-radius: 8px; border: 1px solid var(--border-card);">
            <div id="scatterTitle" style="font-size: 13px; font-weight: 600; color: var(--text-primary);">
                Structure-Activity Relationship (SAR) Landscape
            </div>
            <div style="display: flex; gap: 12px; align-items: center; flex-wrap: wrap;">
                <div style="display: flex; align-items: center; gap: 6px;">
                    <label style="font-size: 11px; font-weight: 600; color: var(--text-secondary); text-transform: uppercase;">X-Axis:</label>
                    <select class="sort" id="scatterXSelect" onchange="updateScatterPlot()" style="font-size: 12px; padding: 4px 8px;"></select>
                </div>
                <div style="display: flex; align-items: center; gap: 6px;">
                    <label style="font-size: 11px; font-weight: 600; color: var(--text-secondary); text-transform: uppercase;">Y-Axis:</label>
                    <select class="sort" id="scatterYSelect" onchange="updateScatterPlot()" style="font-size: 12px; padding: 4px 8px;"></select>
                </div>
                <div style="display: flex; align-items: center; gap: 6px;">
                    <label style="font-size: 11px; font-weight: 600; color: var(--text-secondary); text-transform: uppercase;">Color:</label>
                    <select class="sort" id="scatterColorSelect" onchange="updateScatterPlot()" style="font-size: 12px; padding: 4px 8px;"></select>
                </div>
            </div>
        </div>
        <svg class="scatter-svg" id="scatterSvg" viewBox="0 0 800 480"></svg>
        <div class="scatter-tooltip" id="scatterTooltip"></div>
    </div>

    <!-- Modal Pop-Out Zoom View -->
    <div class="modal-overlay" id="modalOverlay" onclick="closeModalOnBackdrop(event)">
        <div class="modal-card">
            <div class="modal-header">
                <h2 id="modalTitle">Compound Name</h2>
                <button class="modal-close" onclick="closeModal()">&times;</button>
            </div>
            <div class="modal-body">
                <div class="modal-svg-view" id="modalSvg"></div>
                <div class="modal-details">
                    <div class="badge-row" id="modalBadges"></div>
                    <div>
                        <div style="font-size:11px; font-weight:600; color:var(--text-secondary); margin-bottom:4px;">CANONICAL SMILES</div>
                        <div class="modal-smiles-box" id="modalSmiles"></div>
                        <button class="copy-btn" id="modalCopyBtn" onclick="copyCurrentModalSmiles()">Copy SMILES</button>
                    </div>
                    <div>
                        <div style="font-size:11px; font-weight:600; color:var(--text-secondary); margin-bottom:4px; margin-top:8px;">PROPERTIES & METADATA</div>
                        <table class="modal-props-table" id="modalPropsTable"></table>
                    </div>
                </div>
            </div>
            <div class="modal-footer">
                <button class="nav-btn" id="prevBtn" onclick="navigateModal(-1)">&larr; Previous</button>
                <span style="font-size:12px; color:var(--text-secondary);" id="modalIndexText">1 of 10</span>
                <button class="nav-btn" id="nextBtn" onclick="navigateModal(1)">Next &rarr;</button>
            </div>
        </div>
    </div>

    <div class="copy-toast" id="copyToast">SMILES copied to clipboard!</div>

    <script>
        let rawData = [];
        let currentData = [];
        let currentModalIndex = 0;
        let currentView = 'grid';

        function showCopyToast(msg = "SMILES copied to clipboard!") {{
            const toast = document.getElementById("copyToast");
            toast.textContent = msg;
            toast.style.opacity = "1";
            setTimeout(() => {{ toast.style.opacity = "0"; }}, 1500);
        }}

        function escapeHtmlAttr(str) {{
            return String(str ?? '')
                .replace(/&/g, '&amp;')
                .replace(/"/g, '&quot;')
                .replace(/'/g, '&#39;')
                .replace(/</g, '&lt;')
                .replace(/>/g, '&gt;');
        }}

        function copySmiles(smiles, e) {{
            if (e) e.stopPropagation();
            if (navigator.clipboard && navigator.clipboard.writeText) {{
                navigator.clipboard.writeText(smiles);
            }}
            showCopyToast();
        }}

        function copySmilesByIdx(dataIndex, e) {{
            if (e) e.stopPropagation();
            const item = currentData.find(d => d.idx === dataIndex) || rawData.find(d => d.idx === dataIndex);
            if (item && item.smiles) {{
                copySmiles(item.smiles);
            }}
        }}

        function copyCurrentModalSmiles() {{
            const item = currentData[currentModalIndex];
            if (item) copySmiles(item.smiles);
        }}

        function switchView(viewName) {{
            currentView = viewName;
            const btnGrid = document.getElementById("btnGridView");
            const btnScatter = document.getElementById("btnScatterView");
            const gridEl = document.getElementById("gridContainer");
            const scatterEl = document.getElementById("scatterContainer");

            if (viewName === 'grid') {{
                btnGrid.classList.add("active");
                btnScatter.classList.remove("active");
                gridEl.style.display = "grid";
                scatterEl.style.display = "none";
                renderCards(currentData);
            }} else {{
                btnScatter.classList.add("active");
                btnGrid.classList.remove("active");
                gridEl.style.display = "none";
                scatterEl.style.display = "flex";
                renderScatterPlot(currentData);
            }}
        }}

        const _DISC_TIER_RE = /^(?:(?:Grade|Bin|Category|Class|Tier)\\s+)?[A-F](?:\\s*\\(.*\\))?$|^(\\+{{1,6}}|\\*{{1,6}})$/i;
        const _DISC_TIER_CELL_RE = /^(?:(?:Grade|Bin|Category|Class|Tier)\\s+)?([A-F])(?:\\s*\\(.*\\))?$/i;
        const _DISC_POTENCY_DISP_RE = /:\\s*(?:(?:(?:Grade|Bin|Category|Class|Tier)\\s+)?[A-F](?:\\s*\\(|$)|\\++|\\*+)/i;
        const _hasOwn = (obj, k) => Boolean(obj) && Object.prototype.hasOwnProperty.call(obj, k);

        function getBadgesHtml(item) {{
            let badges = '';
            // Primary Potency Badge (continuous or discrete bin/range)
            if (item.potency_display) {{
                const pScore = item.potency_score;
                let potCls = 'badge-blue';
                if (pScore !== null && pScore !== undefined && !isNaN(pScore)) {{
                    const isDiscBin = _DISC_POTENCY_DISP_RE.test(item.potency_display);
                    if (isDiscBin ? pScore >= 3.5 : pScore >= 7.0) potCls = 'badge-green';
                    else if (isDiscBin ? pScore < 1.5 : pScore < 5.5) potCls = 'badge-amber';
                }}
                badges += `<span class="badge ${{potCls}}" title="${{escapeHtmlAttr(item.potency_display)}}">${{escapeHtmlAttr(item.potency_display)}}</span>`;
            }}

            // Ro5 Badge
            if (item.ro5_violations === 0) {{
                badges += '<span class="badge badge-green">Ro5 Pass</span>';
            }} else if (item.ro5_violations === 1) {{
                badges += '<span class="badge badge-amber">Ro5: 1 viol</span>';
            }} else if (item.ro5_violations !== null && item.ro5_violations !== undefined) {{
                badges += `<span class="badge badge-red">Ro5: ${{item.ro5_violations}} viol</span>`;
            }}

            // PAINS Badge
            if (item.pains && item.pains.toLowerCase().includes("clean")) {{
                badges += '<span class="badge badge-green">PAINS Clean</span>';
            }} else if (item.pains) {{
                badges += `<span class="badge badge-red" title="${{escapeHtmlAttr(item.pains)}}">PAINS Alert</span>`;
            }}

            // CNS MPO Badge
            if (item.cns_mpo !== null && item.cns_mpo !== undefined && !isNaN(item.cns_mpo) && item.cns_mpo >= 0) {{
                const cls = item.cns_mpo >= 4.0 ? 'badge-blue' : (item.cns_mpo >= 3.0 ? 'badge-amber' : 'badge-red');
                badges += `<span class="badge ${{cls}}">MPO: ${{item.cns_mpo}}</span>`;
            }}

            return badges;
        }}

        function openModal(dataIndex) {{
            currentModalIndex = currentData.findIndex(item => item.idx === dataIndex);
            if (currentModalIndex === -1) currentModalIndex = 0;
            renderModalContent();
            const overlay = document.getElementById("modalOverlay");
            overlay.classList.add("active");
            const card = document.querySelector(".modal-card");
            if (card) {{
                card.scrollIntoView({{ block: "center", inline: "center", behavior: "auto" }});
            }}
        }}

        function closeModal() {{
            document.getElementById("modalOverlay").classList.remove("active");
        }}

        function closeModalOnBackdrop(e) {{
            if (e.target === document.getElementById("modalOverlay")) {{
                closeModal();
            }}
        }}

        function navigateModal(direction) {{
            const nextIdx = currentModalIndex + direction;
            if (nextIdx >= 0 && nextIdx < currentData.length) {{
                currentModalIndex = nextIdx;
                renderModalContent();
            }}
        }}

        const _SMILES_COL_HINT = {smiles_col_hint_js};
        const _MODAL_SKIP_KEYS = new Set({_js_json(modal_skip_keys)});
        const _isModalSkipKey = k => _MODAL_SKIP_KEYS.has(k);
        function _isCardSkipKey(k) {{
            return (
                _MODAL_SKIP_KEYS.has(k) ||
                _isSmilesColKey(k) ||
                _isInchikeyColKey(k) ||
                /(?:^|[\\s_-])(?:inchi|smarts)$/i.test(k) ||
                /(?:InChI|SMARTS)$/.test(k)
            );
        }}

        const _ANON_PATENT_COL_RE = /^(?:patent_table|pat|[a-z]{{2,4}}\\d{{5,}}[a-z0-9]*)_col_(\\d+)(?:_(?:nm|pactivity|range))?$/i;
        const _NAMED_ASSAY_KEY_RE = /(?:^|_)p?(?:(?:ic|ec|ac|gi|ed|cc)_?50|ki|kd|mic)(?:_|$)/i;
        const _LOG_POTENCY_KEY_RE = /(?:_pactivity$|(?:^|_)p(?:(?:ic|ec|ac|gi|ed|cc)_?50|ki|kd)(?:_|$))/i;
        const _META_UNIT_SUFFIX_RE = /(?:^|_)(?:units?|type|qualifier|operator|relation|comment)$/i;
        const _BIO_TYPE_TOKEN_RE = /(?:wild|cell)_type|(?:sub|geno|pheno)type/gi;
        const _POTENCY_CELL_JS_RE = /^([<>~=\\u2248\\u2264\\u2265]*)\\s*(\\d+(?:\\.\\d+)?)\\s*(pm|nm|um|\\u03bcm|\\u00b5m|mm)?$/i;
        const _LOG_CELL_JS_RE = /^([<>~=\\u2248\\u2264\\u2265]*)\\s*(-?\\d+(?:\\.\\d+)?)$/;

        function _isNamedAssayMetaKey(k) {{
            const kl = String(k).toLowerCase();
            if (kl.endsWith('_range') || _ANON_PATENT_COL_RE.test(kl)) return false;
            if (kl.startsWith('patent_assay_')) return kl.length > 'patent_assay_'.length;
            if (_META_UNIT_SUFFIX_RE.test(kl.replace(_BIO_TYPE_TOKEN_RE, ''))) return false;
            return _NAMED_ASSAY_KEY_RE.test(kl);
        }}

        function _canonPotencyQual(rawQual) {{
            const q = String(rawQual || '').trim();
            if (q.startsWith('<') || q.startsWith('\\u2264')) return '<';
            if (q.startsWith('>') || q.startsWith('\\u2265')) return '>';
            return '';
        }}

        function _defaultUnitForKey(k) {{
            const kl = String(k || '').toLowerCase().replace(/_range$/, '');
            if (kl.endsWith('_pm')) return 'pm';
            if (kl.endsWith('_um') || kl.endsWith('_\\u03bcm') || kl.endsWith('_\\u00b5m')) return 'um';
            if (kl.endsWith('_mm')) return 'mm';
            return 'nm';
        }}

        function _cellPotencySig(val, key = '') {{
            if (!_hasMetaVal(val)) return null;
            const s = String(val).trim();
            if (!s || /^(?:nan|none|nd|n\\/a|-)$/i.test(s)) return null;
            const kl = String(key || '').toLowerCase();
            if (_LOG_POTENCY_KEY_RE.test(kl)) {{
                const lm = _LOG_CELL_JS_RE.exec(s);
                if (lm && Number.isFinite(Number(lm[2]))) {{
                    return {{ pAct: Math.round(Number(lm[2]) * 1000) / 1000, qual: _canonPotencyQual(lm[1]), text: s.toLowerCase() }};
                }}
                return {{ pAct: null, qual: '', text: s.toLowerCase().replace(/\\s+/g, ' ') }};
            }}
            const m = _POTENCY_CELL_JS_RE.exec(s);
            if (m) {{
                const num = Number(m[2]);
                if (Number.isFinite(num) && num > 0) {{
                    const unit = (m[3] || _defaultUnitForKey(kl)).toLowerCase();
                    const offset = unit === 'pm' ? 12.0 : ((unit === 'um' || unit === '\\u03bcm' || unit === '\\u00b5m') ? 6.0 : (unit === 'mm' ? 3.0 : 9.0));
                    const pAct = Math.round((offset - Math.log10(num)) * 1000) / 1000;
                    return {{ pAct, qual: _canonPotencyQual(m[1]), text: s.toLowerCase() }};
                }}
            }}
            return {{ pAct: null, qual: '', text: s.toLowerCase().replace(/\\s+/g, ' ') }};
        }}

        function _potencySigsMatch(a, b) {{
            if (!a || !b) return false;
            if (a.pAct !== null && b.pAct !== null) {{
                return a.qual === b.qual && Math.abs(a.pAct - b.pAct) <= 0.02;
            }}
            if (a.pAct === null && b.pAct === null) {{
                return a.text === b.text;
            }}
            return false;
        }}

        function _authoritativeSigsFromKeys(meta, keys) {{
            const rawSigs = [];
            const derivedSigs = [];
            for (const k of keys) {{
                const sig = _cellPotencySig(_hasOwn(meta, k) ? meta[k] : null, k);
                if (!sig) continue;
                const kl = String(k).toLowerCase();
                if (!/_(?:nm|pactivity)$/i.test(kl) && !_LOG_POTENCY_KEY_RE.test(kl)) {{
                    rawSigs.push(sig);
                }} else {{
                    derivedSigs.push(sig);
                }}
            }}
            return rawSigs.length > 0 ? rawSigs : derivedSigs;
        }}

        function _findMetaKeyCI(meta, targetKey) {{
            const tl = targetKey.toLowerCase();
            return Object.keys(meta).find(mk => mk.toLowerCase() === tl && _hasOwnMetaVal(meta, mk)) || null;
        }}

        function _areAnonColsCovered(meta) {{
            const anonKeysByIdxPrefix = new Map();
            const namedKeysByStem = new Map();
            let hasPatentAssay = false;
            for (const mk of Object.keys(meta)) {{
                if (!_hasOwnMetaVal(meta, mk)) continue;
                const m = _ANON_PATENT_COL_RE.exec(mk);
                if (m) {{
                    if (/_range$/i.test(mk)) continue;
                    const idx = Number(m[1]);
                    const pfx = mk.toLowerCase().replace(/_(?:nm|pactivity)$/i, '');
                    if (!anonKeysByIdxPrefix.has(idx)) anonKeysByIdxPrefix.set(idx, new Map());
                    const byPfx = anonKeysByIdxPrefix.get(idx);
                    if (!byPfx.has(pfx)) byPfx.set(pfx, []);
                    byPfx.get(pfx).push(mk);
                    continue;
                }}
                if (_isNamedAssayMetaKey(mk)) {{
                    if (/^patent_assay_/i.test(mk)) hasPatentAssay = true;
                    const stem = mk.toLowerCase()
                        .replace(/^patent_assay_/i, '')
                        .replace(/_(?:nm|pactivity|range)$/i, '')
                        .replace(/(^|_)p((?:ic|ec|ac|gi|ed|cc)_?50|ki|kd)(?=_|$)/gi, '$1$2')
                        .split('_').filter(Boolean).sort().join('_');
                    if (!namedKeysByStem.has(stem)) namedKeysByStem.set(stem, []);
                    namedKeysByStem.get(stem).push(mk);
                }}
            }}
            if (anonKeysByIdxPrefix.size === 0) return false;
            const anonEntries = [];
            for (const [idx, byPfx] of anonKeysByIdxPrefix.entries()) {{
                const sigs = [];
                for (const pfxKeys of byPfx.values()) {{
                    sigs.push(..._authoritativeSigsFromKeys(meta, pfxKeys));
                }}
                if (sigs.length > 0) anonEntries.push([idx, sigs]);
            }}
            if (anonEntries.length === 0) return false;
            const namedEntries = [];
            for (const [stem, skeys] of namedKeysByStem.entries()) {{
                const sigs = _authoritativeSigsFromKeys(meta, skeys);
                if (sigs.length > 0) namedEntries.push([stem, sigs]);
            }}
            const maxAnonIdx = Math.max(...anonEntries.map(e => e[0]));
            if (namedEntries.length < anonEntries.length || (!hasPatentAssay && namedEntries.length < maxAnonIdx)) {{
                return false;
            }}
            const adj = anonEntries.map(([, asigs]) => {{
                const matches = [];
                for (let j = 0; j < namedEntries.length; j++) {{
                    const nsigs = namedEntries[j][1];
                    if (asigs.every(asig => nsigs.some(nsig => _potencySigsMatch(asig, nsig)))) {{
                        matches.push(j);
                    }}
                }}
                return matches;
            }});
            if (adj.some(m => m.length === 0)) return false;
            const matchRight = new Array(namedEntries.length).fill(-1);
            function _dfs(i, seen) {{
                for (const j of adj[i]) {{
                    if (seen[j]) continue;
                    seen[j] = true;
                    if (matchRight[j] === -1 || _dfs(matchRight[j], seen)) {{
                        matchRight[j] = i;
                        return true;
                    }}
                }}
                return false;
            }}
            for (let i = 0; i < anonEntries.length; i++) {{
                if (!_dfs(i, new Array(namedEntries.length).fill(false))) return false;
            }}
            return true;
        }}

        function _hasMetaKeyCI(meta, targetKey) {{
            const tl = targetKey.toLowerCase();
            return Object.keys(meta).some(mk => mk.toLowerCase() === tl && _hasOwnMetaVal(meta, mk));
        }}

        function _metaPriorityScore(k) {{
            const kl = k.toLowerCase();
            if (kl === 'patent_example_id' || kl === 'paper_compound_id') return 1;
            if (kl.startsWith('patent_assay_') && !kl.endsWith('_range')) return 2;
            if (_isNamedAssayMetaKey(kl) || /^(?:pchembl_value|activity_value|activity_bin)$/i.test(kl)) return 3;
            if (kl === 'target_name' || kl === 'chiral_peak') return 4;
            if (kl.startsWith('patent_table_col_') && !kl.endsWith('_range')) return 5;
            if (kl === 'patent_id' || kl === 'document_id' || kl === 'source') return 20;
            return 10;
        }}

        function _hasMetaVal(v) {{
            return v !== null && v !== undefined && v !== '';
        }}

        function _hasOwnMetaVal(meta, k) {{
            return _hasOwn(meta, k) && _hasMetaVal(meta[k]);
        }}

        function _formatMetaKey(k) {{
            const raw = String(k).trim();
            return raw.replace(/[_\\s]+/g, ' ').trim() || raw;
        }}

        function _formatMetaVal(v, k = '') {{
            if (v === null || v === undefined) return '';
            const s = String(v).trim();
            if (k && /(?:^|[\\s_-])(?:id|cas|cid|sid|code|number|no\\.?|num|idx)$/i.test(String(k))) {{
                return s;
            }}
            if (/^-?\\d+\\.\\d{{3,}}$/.test(s)) {{
                const num = Number(s);
                if (num === 0 || Math.abs(num) >= 0.1) {{
                    return parseFloat(num.toFixed(2)).toString();
                }}
                return parseFloat(num.toPrecision(2)).toString();
            }}
            return s;
        }}

        function _normFieldKey(s) {{
            return String(s ?? '').toLowerCase().replace(/[^a-z0-9]+/g, '');
        }}

        function _isSmilesColKey(k) {{
            if (_SMILES_COL_HINT && k === _SMILES_COL_HINT) return true;
            const nk = _normFieldKey(k);
            return nk === 'smi' || nk.endsWith('smiles');
        }}

        function _isInchikeyColKey(k) {{
            return _normFieldKey(k).endsWith('inchikey');
        }}

        function _isFormulaColKey(k) {{
            const nk = _normFieldKey(k);
            return nk === 'molecularformula' || nk === 'formula' || nk === 'molformula' || nk === 'chemicalformula';
        }}

        const _normMetaKey = k => String(k).toLowerCase().replace(/[\\W_]+/g, '_').replace(/^_+|_+$/g, '');

        function _isVisibleMetaKey(meta, k, skipOrFn = _isModalSkipKey, anonCovered = null) {{
            const isSkipped = typeof skipOrFn === 'function'
                ? skipOrFn(k)
                : Boolean(skipOrFn && typeof skipOrFn.has === 'function' && skipOrFn.has(k));
            if (isSkipped || !_hasOwnMetaVal(meta, k)) return false;
            const nk = _normMetaKey(k);
            if (nk.endsWith('_range')) {{
                const baseNorm = nk.slice(0, -6);
                const hasBase = Object.keys(meta).some(mk => mk !== k && _hasOwnMetaVal(meta, mk) && (
                    _normMetaKey(mk) === baseNorm ||
                    (nk === 'activity_range' && (_normMetaKey(mk) === 'activity_value' || _normMetaKey(mk) === 'activity_bin')) ||
                    (nk === 'standard_range' && _normMetaKey(mk) === 'standard_value')
                ));
                if (hasBase) return false;
            }}
            const anonMatch = _ANON_PATENT_COL_RE.exec(k);
            if (anonMatch) {{
                if (anonCovered === null ? _areAnonColsCovered(meta) : anonCovered) return false;
                if (!/^patent_table_col_/i.test(k)) {{
                    const colIdx = anonMatch[1];
                    const kl = k.toLowerCase();
                    if (kl.endsWith('_range')) {{
                        const canonRng = _findMetaKeyCI(meta, `patent_table_col_${{colIdx}}_range`);
                        if (canonRng && String(meta[canonRng]).trim().toLowerCase() === String(meta[k]).trim().toLowerCase()) {{
                            return false;
                        }}
                    }} else if (/_(?:nm|pactivity)$/i.test(kl)) {{
                        const sfx = kl.endsWith('_nm') ? '_nM' : '_pActivity';
                        const sameUnitCanon = _findMetaKeyCI(meta, `patent_table_col_${{colIdx}}${{sfx}}`);
                        const rawCanon = _findMetaKeyCI(meta, `patent_table_col_${{colIdx}}`);
                        const refKey = sameUnitCanon || rawCanon;
                        if (refKey && _potencySigsMatch(_cellPotencySig(meta[k], k), _cellPotencySig(meta[refKey], refKey))) {{
                            return false;
                        }}
                    }} else {{
                        const canonKeys = ['', '_nM', '_pActivity']
                            .map(sfx => _findMetaKeyCI(meta, `patent_table_col_${{colIdx}}${{sfx}}`))
                            .filter(Boolean);
                        const canonSigs = _authoritativeSigsFromKeys(meta, canonKeys);
                        const kSig = _cellPotencySig(meta[k], k);
                        if (kSig && canonSigs.some(cs => _potencySigsMatch(kSig, cs))) {{
                            return false;
                        }}
                    }}
                }}
            }}
            if (k.length > 3 && /_nm$/i.test(k) && _hasMetaKeyCI(meta, k.slice(0, -3))) return false;
            if (k.length > 10 && /_pactivity$/i.test(k)) {{
                const stem = k.slice(0, -10);
                if (_hasMetaKeyCI(meta, stem) || _hasMetaKeyCI(meta, `${{stem}}_nM`)) return false;
            }}
            return true;
        }}

        function _visibleMetaKeys(meta, skipOrFn = _isModalSkipKey) {{
            const anonCovered = _areAnonColsCovered(meta);
            return Object.keys(meta)
                .filter(k => _isVisibleMetaKey(meta, k, skipOrFn, anonCovered))
                .sort((a, b) => _metaPriorityScore(a) - _metaPriorityScore(b));
        }}

        function _formatMetaEntryVal(meta, k) {{
            let valStr = _formatMetaVal(_hasOwn(meta, k) ? meta[k] : null, k);
            const nk = _normMetaKey(k);
            const targetRng = nk + '_range';
            const hasActBin = Object.keys(meta).some(mk => _normMetaKey(mk) === 'activity_bin' && _hasOwnMetaVal(meta, mk));
            const allowActRng = nk === 'activity_bin' || (nk === 'activity_value' && !hasActBin);
            const rngKey = Object.keys(meta).find(mk => _hasOwnMetaVal(meta, mk) && (
                _normMetaKey(mk) === targetRng ||
                (allowActRng && _normMetaKey(mk) === 'activity_range') ||
                (nk === 'standard_value' && _normMetaKey(mk) === 'standard_range')
            ));
            const rngVal = rngKey ? meta[rngKey] : null;
            if (rngVal && !valStr.includes(String(rngVal))) {{
                valStr = `${{valStr}} (${{rngVal}})`;
            }}
            return valStr;
        }}
        function renderModalContent() {{
            const item = currentData[currentModalIndex];
            if (!item) return;

            document.getElementById("modalTitle").textContent = item.id;
            document.getElementById("modalSvg").innerHTML = item.svg;
            document.getElementById("modalSmiles").textContent = item.smiles;
            document.getElementById("modalIndexText").textContent = `${{currentModalIndex + 1}} of ${{currentData.length}}`;
            document.getElementById("modalBadges").innerHTML = getBadgesHtml(item);

            document.getElementById("prevBtn").disabled = currentModalIndex === 0;
            document.getElementById("nextBtn").disabled = currentModalIndex === currentData.length - 1;

            const fmtVal = (v, suffix = '') => (v !== null && v !== undefined && !isNaN(v)) ? `${{v}}${{suffix}}` : 'N/A';

            let propsHtml = '<colgroup><col style="width:46%"><col style="width:54%"></colgroup>';
            if (item.potency_display) {{
                propsHtml += `<tr><td class="label">Primary Bioactivity:</td><td class="val">${{escapeHtmlAttr(item.potency_display)}}</td></tr>`;
            }}
            propsHtml += `
                <tr><td class="label">Molecular Weight:</td><td class="val">${{fmtVal(item.mw, ' g/mol')}}</td></tr>
                <tr><td class="label">cLogP (Crippen):</td><td class="val">${{fmtVal(item.clogp)}}</td></tr>
                <tr><td class="label">TPSA:</td><td class="val">${{fmtVal(item.tpsa, ' Å²')}}</td></tr>
                <tr><td class="label">H-Bond Donors (HBD):</td><td class="val">${{fmtVal(item.hbd)}}</td></tr>
                <tr><td class="label">H-Bond Acceptors (HBA):</td><td class="val">${{fmtVal(item.hba)}}</td></tr>
                <tr><td class="label">Rotatable Bonds:</td><td class="val">${{fmtVal(item.rotb)}}</td></tr>
                <tr><td class="label">Fraction Csp3:</td><td class="val">${{fmtVal(item.fsp3)}}</td></tr>
                <tr><td class="label">QED Drug-Likeness:</td><td class="val">${{fmtVal(item.qed)}}</td></tr>
                <tr><td class="label">Pfizer CNS MPO:</td><td class="val">${{fmtVal(item.cns_mpo, ' / 6.0')}}</td></tr>
                <tr><td class="label">PAINS / Brenk Alerts:</td><td class="val">${{escapeHtmlAttr(item.pains || 'Clean')}}</td></tr>
            `;

            if (item.meta) {{
                for (const k of _visibleMetaKeys(item.meta)) {{
                    const dispK = _formatMetaKey(k);
                    const valStr = _formatMetaEntryVal(item.meta, k);
                    propsHtml += `<tr><td class="label" title="${{escapeHtmlAttr(k)}}">${{escapeHtmlAttr(dispK)}}:</td><td class="val">${{escapeHtmlAttr(valStr)}}</td></tr>`;
                }}
            }}
            document.getElementById("modalPropsTable").innerHTML = propsHtml;
        }}

        // Keyboard Navigation
        document.addEventListener("keydown", (e) => {{
            const overlay = document.getElementById("modalOverlay");
            if (!overlay || !overlay.classList.contains("active")) return;
            if (e.key === "Escape") closeModal();
            else if (e.key === "ArrowLeft") navigateModal(-1);
            else if (e.key === "ArrowRight") navigateModal(1);
        }});

        function renderCards(data) {{
            const container = document.getElementById("gridContainer");
            const countEl = document.getElementById("itemCount");
            countEl.textContent = `Showing ${{data.length}} of ${{rawData.length}} compounds`;
            
            if (data.length === 0) {{
                container.innerHTML = '<div style="grid-column: 1/-1; text-align:center; padding: 40px; color: var(--text-secondary);">No matching compounds found.</div>';
                return;
            }}

            const fmtCardVal = (v, decimals = 1) => {{
                if (v === null || v === undefined || isNaN(v)) return 'N/A';
                const num = Number(v);
                return Number.isInteger(num) ? String(num) : parseFloat(num.toFixed(decimals)).toString();
            }};

            container.innerHTML = data.map(item => {{
                let extraMeta = '';
                if (item.meta) {{
                    const keys = _visibleMetaKeys(item.meta, _isCardSkipKey);
                    if (keys.length > 0) {{
                        extraMeta = '<div class="meta-list">' + keys.slice(0, 3).map(k => {{
                            const dispK = _formatMetaKey(k);
                            const valStr = _formatMetaEntryVal(item.meta, k);
                            return `<div class="meta-row"><span class="meta-key" title="${{escapeHtmlAttr(k)}}">${{escapeHtmlAttr(dispK)}}:</span><span class="meta-val" title="${{escapeHtmlAttr(valStr)}}">${{escapeHtmlAttr(valStr)}}</span></div>`;
                        }}).join('') + '</div>';
                    }}
                }}

                const badgesHtml = getBadgesHtml(item);
                const fullTitle = (item.meta && _hasOwnMetaVal(item.meta, 'compound_name') && item.meta.compound_name !== item.id)
                    ? `${{item.id}} — ${{item.meta.compound_name}}`
                    : item.id;

                return `
                <div class="card" onclick="openModal(${{item.idx}})">
                    <div class="zoom-hint">Click to Zoom</div>
                    <div class="card-header" title="${{escapeHtmlAttr(fullTitle)}}">${{escapeHtmlAttr(item.id)}}</div>
                    <div class="badge-row">${{badgesHtml}}</div>
                    <div class="svg-container">${{item.svg}}</div>
                    <div class="props-grid">
                        <div class="prop-item"><span class="prop-label">MW</span><span class="prop-val">${{fmtCardVal(item.mw, 1)}}</span></div>
                        <div class="prop-item"><span class="prop-label">cLogP</span><span class="prop-val">${{fmtCardVal(item.clogp, 2)}}</span></div>
                        <div class="prop-item"><span class="prop-label">TPSA</span><span class="prop-val">${{fmtCardVal(item.tpsa, 1)}}</span></div>
                    </div>
                    <div class="smiles-wrapper">
                        <div class="smiles-row" title="Click to copy SMILES" onclick="copySmilesByIdx(${{item.idx}}, event)">${{escapeHtmlAttr(item.smiles)}}</div>
                    </div>
                    ${{extraMeta}}
                </div>
                `;
            }}).join('');
        }}

        let scatterInitialized = false;
        const _ACTIVE_COL_HINT = {active_col_hint_js};

        function initScatterControls() {{
            if (scatterInitialized || rawData.length === 0) return;
            const xSel = document.getElementById("scatterXSelect");
            const ySel = document.getElementById("scatterYSelect");
            const cSel = document.getElementById("scatterColorSelect");
            if (!xSel || !ySel || !cSel) return;

            const propKeys = new Set();
            rawData.forEach(item => {{
                if (item.props) {{
                    Object.keys(item.props).forEach(k => {{
                        if (typeof item.props[k] === 'number' && !isNaN(item.props[k])) {{
                            propKeys.add(k);
                        }}
                    }});
                }}
            }});

            const keysArr = Array.from(propKeys);
            const priorityOrder = [_ACTIVE_COL_HINT, "pIC50", "pchembl_value", "pKi", "IC50_nM", "activity_value", "activity_bin", "cLogP", "MW", "TPSA", "CNS_MPO", "QED", "RotB", "Fsp3", "HBD", "HBA", "Ro5_Violations"].filter(Boolean);
            keysArr.sort((a, b) => {{
                const ia = priorityOrder.indexOf(a);
                const ib = priorityOrder.indexOf(b);
                if (ia !== -1 && ib !== -1) return ia - ib;
                if (ia !== -1) return -1;
                if (ib !== -1) return 1;
                return a.localeCompare(b);
            }});

            [xSel, ySel, cSel].forEach(sel => {{
                sel.innerHTML = keysArr.map(k => `<option value="${{escapeHtmlAttr(k)}}">${{escapeHtmlAttr(k)}}</option>`).join('');
            }});

            const _isMetaAxis = k => /(?:^|[\\W_])(?:units?|qualifier|operator|relation|flag|comment|range)(?:$|[\\W_])/i.test(k);
            let potencyKey = (_ACTIVE_COL_HINT && keysArr.includes(_ACTIVE_COL_HINT))
                ? _ACTIVE_COL_HINT
                : keysArr.find(k => !_isMetaAxis(k) && /(?:^|[\\W_])(?:p?(?:ic|ec|ac|gi)_?50|p?ki|p?kd|p?mic|pchembl|pactivity|potency|activity)/i.test(k));
            const clogpKey = keysArr.find(k => /^cLogP$|^LogP$/i.test(k)) || "cLogP";
            const mwKey = keysArr.find(k => /^MW$/i.test(k)) || "MW";
            const colorKey = keysArr.find(k => /CNS_MPO|QED|Ro5_Violations/i.test(k)) || "Ro5_Violations";

            if (potencyKey) {{
                ySel.value = potencyKey;
                const xCandidates = [clogpKey, mwKey, ...keysArr].filter(k => k !== potencyKey && keysArr.includes(k));
                xSel.value = xCandidates[0] || keysArr[0];
            }} else {{
                xSel.value = keysArr.includes(mwKey) ? mwKey : keysArr[0];
                ySel.value = keysArr.includes(clogpKey) ? clogpKey : (keysArr[1] || keysArr[0]);
            }}
            cSel.value = keysArr.includes(colorKey) ? colorKey : (keysArr[2] || keysArr[0]);

            scatterInitialized = true;
        }}

        function updateScatterPlot() {{
            renderScatterPlot(currentData);
        }}

        const _SCALAR_ITEM_FIELDS = new Set(['mw', 'clogp', 'tpsa', 'hbd', 'hba', 'rotb', 'fsp3', 'qed', 'cns_mpo', 'ro5_violations']);

        function getCompoundValue(item, propName) {{
            if (!item || !propName) return null;
            if (_hasOwn(item.props, propName) && item.props[propName] !== undefined && item.props[propName] !== null) {{
                return item.props[propName];
            }}
            const lk = propName.toLowerCase();
            if (_SCALAR_ITEM_FIELDS.has(lk) && item[lk] !== undefined && item[lk] !== null && typeof item[lk] === 'number') {{
                return item[lk];
            }}
            if (_hasOwn(item.meta, propName) && item.meta[propName] !== undefined) {{
                const s = String(item.meta[propName]).trim();
                const v = /^[+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:e[+-]?\\d+)?$/i.test(s) ? Number(s) : NaN;
                return Number.isFinite(v) ? v : null;
            }}
            return null;
        }}

        function _discreteTierLabels(propName, data) {{
            if (!propName) return null;
            const labels = new Map();
            for (const d of data) {{
                if (!_hasOwnMetaVal(d.meta, propName)) continue;
                const raw = String(d.meta[propName]).trim();
                if (!_DISC_TIER_RE.test(raw)) continue;
                const num = getCompoundValue(d, propName);
                if (num !== null && !isNaN(num) && !labels.has(num)) {{
                    labels.set(num, raw.replace(/\\s*\\(.*\\)$/, '').trim());
                }}
            }}
            return labels.size > 0 ? labels : null;
        }}

        function _formatAxisTick(val, tierLabels) {{
            if (tierLabels && tierLabels.has(val)) {{
                return tierLabels.get(val);
            }}
            return Math.abs(val) >= 100 ? val.toFixed(0) : (Math.abs(val) >= 10 ? val.toFixed(1) : val.toFixed(2));
        }}

        function renderScatterPlot(data) {{
            initScatterControls();
            const countEl = document.getElementById("itemCount");
            if (countEl) countEl.textContent = `Showing ${{data.length}} of ${{rawData.length}} compounds`;
            const svg = document.getElementById("scatterSvg");
            const xSel = document.getElementById("scatterXSelect");
            const ySel = document.getElementById("scatterYSelect");
            const cSel = document.getElementById("scatterColorSelect");
            const titleEl = document.getElementById("scatterTitle");
            if (!svg || !xSel || !ySel) return;

            const xProp = xSel.value || "MW";
            const yProp = ySel.value || "cLogP";
            const cProp = cSel ? cSel.value : "Ro5_Violations";

            if (titleEl) {{
                titleEl.textContent = `Structure-Activity Relationship (SAR) Landscape: ${{yProp}} vs. ${{xProp}}`;
            }}

            const discX = _discreteTierLabels(xProp, data);
            const discY = _discreteTierLabels(yProp, data);

            const validPoints = data.map(item => {{
                const vx = getCompoundValue(item, xProp);
                const vy = getCompoundValue(item, yProp);
                const vc = getCompoundValue(item, cProp);
                return {{ item, x: vx, y: vy, c: vc }};
            }}).filter(p => p.x !== null && !isNaN(p.x) && p.y !== null && !isNaN(p.y));

            if (validPoints.length === 0) {{
                const emptyMsg = data.length === 0
                    ? "No matching compounds found."
                    : "No compounds with numerical values for selected axes";
                svg.innerHTML = `<text x="400" y="240" text-anchor="middle" fill="var(--text-secondary)" font-size="14">${{emptyMsg}}</text>`;
                return;
            }}

            const padL = 65, padR = 40, padT = 30, padB = 55;
            const w = 800, h = 480;
            const pw = w - padL - padR, ph = h - padT - padB;

            let minX = Math.min(...validPoints.map(p => p.x));
            let maxX = Math.max(...validPoints.map(p => p.x));
            let minY = Math.min(...validPoints.map(p => p.y));
            let maxY = Math.max(...validPoints.map(p => p.y));

            if (minX === maxX) {{ minX -= 1; maxX += 1; }}
            if (minY === maxY) {{ minY -= 1; maxY += 1; }}

            const spanX = maxX - minX;
            const spanY = maxY - minY;
            const plotMinX = minX - 0.05 * spanX;
            const plotMaxX = maxX + 0.05 * spanX;
            const plotMinY = minY - 0.05 * spanY;
            const plotMaxY = maxY + 0.05 * spanY;

            const xRange = plotMaxX - plotMinX;
            const yRange = plotMaxY - plotMinY;

            const scaleX = (val) => padL + ((val - plotMinX) / xRange) * pw;
            const scaleY = (val) => padT + ph - ((val - plotMinY) / yRange) * ph;

            let s = '';

            // Gridlines & Ticks
            s += '<g stroke="rgba(127,127,127,0.15)" stroke-width="1">';
            const xTicks = discX
                ? Array.from(discX.keys()).sort((a, b) => a - b)
                : [0, 1, 2, 3, 4, 5].map(i => plotMinX + (i / 5.0) * xRange);
            xTicks.forEach(tickXVal => {{
                if (tickXVal < plotMinX || tickXVal > plotMaxX) return;
                const px = scaleX(tickXVal);
                s += `<line x1="${{px}}" y1="${{padT}}" x2="${{px}}" y2="${{padT + ph}}" />`;
                const labelX = escapeHtmlAttr(_formatAxisTick(tickXVal, discX));
                s += `<text x="${{px}}" y="${{padT + ph + 18}}" font-size="11" fill="var(--text-secondary)" text-anchor="middle">${{labelX}}</text>`;
            }});
            const yTicks = discY
                ? Array.from(discY.keys()).sort((a, b) => a - b)
                : [0, 1, 2, 3, 4, 5].map(i => plotMinY + (i / 5.0) * yRange);
            yTicks.forEach(tickYVal => {{
                if (tickYVal < plotMinY || tickYVal > plotMaxY) return;
                const py = scaleY(tickYVal);
                s += `<line x1="${{padL}}" y1="${{py}}" x2="${{padL + pw}}" y2="${{py}}" />`;
                const labelY = escapeHtmlAttr(_formatAxisTick(tickYVal, discY));
                s += `<text x="${{padL - 10}}" y="${{py + 4}}" font-size="11" fill="var(--text-secondary)" text-anchor="end">${{labelY}}</text>`;
            }});
            s += '</g>';

            // Ro5 Sweet Spot Box only if MW vs cLogP is selected
            const isMwX = /^MW$/i.test(xProp);
            const isLogpY = /^cLogP$|^LogP$/i.test(yProp);
            if (isMwX && isLogpY) {{
                const sweetX1 = Math.max(scaleX(Math.max(plotMinX, 0)), padL);
                const sweetX2 = Math.min(scaleX(500), padL + pw);
                const sweetY1 = Math.max(scaleY(5.0), padT);
                const sweetY2 = Math.min(scaleY(0.0), padT + ph);
                if (sweetX2 > sweetX1 && sweetY2 > sweetY1) {{
                    s += `<rect x="${{sweetX1}}" y="${{sweetY1}}" width="${{sweetX2 - sweetX1}}" height="${{sweetY2 - sweetY1}}" fill="rgba(16,185,129,0.08)" stroke="rgba(16,185,129,0.35)" stroke-width="1.5" rx="6" />`;
                    s += `<text x="${{sweetX1 + 10}}" y="${{sweetY1 + 18}}" font-size="11" font-weight="600" fill="var(--badge-green-txt)">Lipinski Ro5 Sweet Spot</text>`;
                }}
            }}

            // Axes
            s += `<line x1="${{padL}}" y1="${{padT + ph}}" x2="${{padL + pw}}" y2="${{padT + ph}}" stroke="var(--border-card)" stroke-width="1.5" />`;
            s += `<line x1="${{padL}}" y1="${{padT}}" x2="${{padL}}" y2="${{padT + ph}}" stroke="var(--border-card)" stroke-width="1.5" />`;
            s += `<text x="${{padL + pw / 2}}" y="${{h - 12}}" font-size="12" font-weight="600" fill="var(--text-primary)" text-anchor="middle">${{escapeHtmlAttr(xProp)}}</text>`;
            s += `<text x="18" y="${{padT + ph / 2}}" font-size="12" font-weight="600" fill="var(--text-primary)" text-anchor="middle" transform="rotate(-90 18 ${{padT + ph/2}})">${{escapeHtmlAttr(yProp)}}</text>`;

            // Min/Max for color interpolation
            const cVals = validPoints.map(p => p.c).filter(v => v !== null && !isNaN(v));
            const minC = cVals.length > 0 ? Math.min(...cVals) : 0;
            const maxC = cVals.length > 0 ? Math.max(...cVals) : 1;
            const spanC = (maxC - minC) || 1.0;

            // Data Points
            validPoints.forEach(p => {{
                const cx = scaleX(p.x);
                const cy = scaleY(p.y);
                let color = "#4f46e5";

                if (/ro5_violations/i.test(cProp)) {{
                    if (p.c === 0) color = "#10b981";
                    else if (p.c === 1) color = "#f59e0b";
                    else color = "#f43f5e";
                }} else if (p.c !== null && !isNaN(p.c)) {{
                    const frac = Math.max(0, Math.min(1, (p.c - minC) / spanC));
                    const r = Math.round(79 + frac * (16 - 79));
                    const g = Math.round(70 + frac * (185 - 70));
                    const b = Math.round(229 + frac * (129 - 229));
                    color = `rgb(${{r}},${{g}},${{b}})`;
                }}

                s += `<circle class="scatter-dot" cx="${{cx}}" cy="${{cy}}" r="6" fill="${{color}}" stroke="#ffffff" stroke-width="1.5"
                    data-idx="${{p.item.idx}}"
                    onclick="openModal(${{p.item.idx}})"
                    onmousemove="showScatterTooltip(event, ${{p.item.idx}})"
                    onmouseleave="hideScatterTooltip()" />`;
            }});

            svg.innerHTML = s;
        }}

        function showScatterTooltip(e, idx) {{
            const item = currentData.find(d => d.idx === idx);
            if (!item) return;
            const xName = document.getElementById("scatterXSelect")?.value || "MW";
            const yName = document.getElementById("scatterYSelect")?.value || "cLogP";
            const cName = document.getElementById("scatterColorSelect")?.value || "";
            const xVal = getCompoundValue(item, xName);
            const yVal = getCompoundValue(item, yName);
            const cVal = cName ? getCompoundValue(item, cName) : null;
            const tooltip = document.getElementById("scatterTooltip");
            const rect = document.getElementById("scatterContainer").getBoundingClientRect();

            const formatPropVal = (pName, v) => {{
                if (_hasOwnMetaVal(item.meta, pName) && isNaN(Number(item.meta[pName]))) {{
                    return escapeHtmlAttr(item.meta[pName]);
                }}
                return (v !== null && v !== undefined && !isNaN(v)) ? (Math.abs(v) >= 100 ? v.toFixed(1) : v.toFixed(3)) : 'N/A';
            }};
            let colorRow = '';
            if (cName && cName !== xName && cName !== yName && cVal !== null && cVal !== undefined) {{
                colorRow = `<br><b>${{escapeHtmlAttr(cName)}}:</b> ${{formatPropVal(cName, cVal)}}`;
            }}

            tooltip.innerHTML = `<b>${{escapeHtmlAttr(item.id)}}</b><br><b>${{escapeHtmlAttr(xName)}}:</b> ${{formatPropVal(xName, xVal)}}<br><b>${{escapeHtmlAttr(yName)}}:</b> ${{formatPropVal(yName, yVal)}}${{colorRow}}`;
            tooltip.style.left = `${{e.clientX - rect.left + 14}}px`;
            tooltip.style.top = `${{e.clientY - rect.top + 14}}px`;
            tooltip.style.opacity = "1";
        }}

        function hideScatterTooltip() {{
            const tooltip = document.getElementById("scatterTooltip");
            tooltip.style.opacity = "0";
        }}

        const _UNSCOPED_META_SKIP = new Set([
            'ROMol', 'svg_data', '_potency_sort_score', 'all_endpoints',
            'MW', 'cLogP', 'TPSA', 'HBD', 'HBA', 'RotB', 'Fsp3', 'QED',
            'CNS_MPO', 'Ro5_Violations'
        ]);
        const _ID_LIKE_META_RE = /(?:^|[^a-z0-9])(?:id|cid|cas|year|example|peak|lot|batch)(?:[^a-z0-9]|$)/i;
        const _PLAIN_NUM_RE = /^[+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:e[+-]?\\d+)?$/i;
        const _INCHIKEY_QUERY_RE = /^[A-Z]{{6,14}}(?:-[A-Z]{{1,10}}(?:-[A-Z])?)?$/i;
        const _OP_LIST = ['<=', '>=', '!=', '==', '=', '<', '>', ':'];
        const _CONC_FACTOR_NM = Object.assign(Object.create(null), {{
            fm: 1e-6, pm: 1e-3, nm: 1, um: 1e3, 'µm': 1e3, 'μm': 1e3, mm: 1e6, m: 1e9
        }});

        function _normConcUnitJs(u) {{
            if (u === null || u === undefined) return null;
            const raw = String(u).trim();
            if (!raw) return null;
            if (raw.length === 1) return raw === 'M' ? 'm' : null;
            const s = raw.toLowerCase()
                .replace(/µ/g, 'μ')
                .replace(/·/g, '.')
                .replace(/−/g, '-')
                .replace(/⁻¹/g, '-1');
            if (s === 'molar') return 'm';
            if (_CONC_FACTOR_NM[s] !== undefined) return s;
            const mm = s.match(/^(f|p|n|u|μ|m)?mol(?:\\s*(?:\\/|[\\s._])\\s*l(?:\\s*(?:\\^?\\s*-1|_1))?|l-1)$/);
            if (mm) return (mm[1] || '') + 'm';
            return null;
        }}

        function _parseCensoredCellJs(sv) {{
            const m = String(sv ?? '').trim().match(/^(<=|>=|[<>≤≥])\\s*([+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:e[+-]?\\d+)?)(?:\\s*(.+))?$/);
            if (!m) return null;
            const rawBound = Number(m[2]);
            if (isNaN(rawBound)) return null;
            const unitRaw = (m[3] || '').trim();
            if (!unitRaw) return {{ op: m[1], bound: rawBound, concUnit: null }};
            const uLow = unitRaw.toLowerCase();
            if (uLow === 'da' || uLow === 'g/mol') return {{ op: m[1], bound: rawBound, concUnit: null }};
            const cu = _normConcUnitJs(unitRaw);
            return cu !== null ? {{ op: m[1], bound: rawBound, concUnit: cu }} : null;
        }}

        function _isConcStringCellJs(sv) {{
            if (_DISC_TIER_CELL_RE.test(sv)) return false;
            const s = String(sv ?? '').trim().replace(/(\\d),(?=\\d{{3}}\\b)/g, '$1');
            const mSingle = s.match(/^[<>≤≥~=]*\\s*[+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:e[+-]?\\d+)?\\s*(.+)$/);
            if (mSingle && _normConcUnitJs(mSingle[1]) !== null) return true;
            const mRng = s.match(/^(?:\\d+\\.?\\d*|\\.\\d+)(?:e[+-]?\\d+)?(?:\\s*((?:f|p|n|u|µ|μ|m)?mol(?:\\s*(?:\\/|\\s|[._·])\\s*l(?:\\s*(?:\\^?\\s*[-−]1|⁻¹)|_1)?|l[-−]1)|molar|fm|pm|nm|um|µm|μm|mm|M))?\\s*(?:[-–—]|\\bto\\b)\\s*(?:\\d+\\.?\\d*|\\.\\d+)(?:e[+-]?\\d+)?\\s*(.+)$/i);
            return Boolean(mRng && _normConcUnitJs(mRng[2]) !== null && (!mRng[1] || _normConcUnitJs(mRng[1]) !== null));
        }}

        function _normalizePotencyColTokensJs(s) {{
            return String(s ?? '').replace(
                /(^|[^a-z0-9])(p?(?:ic|ec|ac|gi))[\\s_-]?50(?![a-z0-9])/gi,
                (_, pre, tok) => pre + tok + '50'
            );
        }}

        function _splitFieldTokens(s) {{
            const prepped = _normalizePotencyColTokensJs(String(s ?? '').toLowerCase());
            return prepped.split(/[^a-z0-9]+/).filter(Boolean);
        }}

        function _colUnitFactorNm(colName) {{
            const norm = _normalizePotencyColTokensJs(colName);
            const re = /(?:^|[^a-z0-9µμ]|(?:ic50|ki|kd|ec50|ac50|gi50|mic))((?:f|p|n|u|µ|μ|m)?mol(?:\\s*(?:\\/|\\s|[._·])\\s*l(?:\\s*(?:\\^?\\s*[-−]1|⁻¹)|_1)?|l[-−]1)|molar|fm|pm|nm|um|µm|μm|mm|M)(?=[^a-z0-9]|$)/gi;
            const allUnits = [];
            const readoutUnits = [];
            let m;
            while ((m = re.exec(norm)) !== null) {{
                const rawTok = m[1];
                const uStart = m.index + m[0].length - rawTok.length;
                const prefix = norm.slice(0, uStart);
                if (rawTok.length === 1) {{
                    const prevCh = uStart > 0 ? norm[uStart - 1] : '';
                    const hasPotPre = /(?:^|[^a-z0-9])(?:p?(?:ic|ec|ac|gi)50|ki|kd|mic|potency|activity|affinity|conc\\w*|value)(?:[^a-z0-9]|$)/i.test(prefix);
                    if (
                        rawTok !== 'M' ||
                        !prevCh ||
                        /[a-z0-9µμ]/i.test(prevCh) ||
                        !('([{{'.includes(prevCh) || !norm.slice(uStart + 1).trim() || hasPotPre)
                    ) continue;
                }}
                const unit = _normConcUnitJs(rawTok);
                if (!unit) continue;
                allUnits.push(unit);
                if (!/(?:^|[^a-z0-9<>≤≥.=])\\d+(?:\\.\\d+)?[\\s_]*$/i.test(prefix)) {{
                    readoutUnits.push(unit);
                }}
            }}
            const chosen = readoutUnits.length > 0
                ? readoutUnits[readoutUnits.length - 1]
                : (allUnits.length > 0 ? allUnits[allUnits.length - 1] : null);
            return chosen ? (_CONC_FACTOR_NM[chosen] || 1) : 1;
        }}

        function _isFieldIdent(s) {{
            return Boolean(
                s &&
                /^[a-zA-Z0-9][a-zA-Z0-9_.'-]*$/.test(s) &&
                /[a-zA-Z]/.test(s) &&
                !/^[+-]?(?:\\d+\\.?\\d*|\\.\\d+)e[+-]?\\d+$/i.test(s)
            );
        }}

        let _schemaCacheRef = null;
        let _schemaCacheLen = -1;
        let _schemaCache = null;

        function _getSchemaIndex() {{
            if (_schemaCache && _schemaCacheRef === rawData && _schemaCacheLen === rawData.length) {{
                return _schemaCache;
            }}
            const metaCols = new Set();
            const propCols = new Set();
            const numericMetaSeen = new Set();
            const censoredMetaSeen = new Set();
            const concStrMetaSeen = new Set();
            const tierMetaSeen = new Set();
            const boolStats = new Map();
            let hasContinuousPotency = false;

            for (const item of rawData) {{
                const dispIsDisc = item.potency_display && _DISC_POTENCY_DISP_RE.test(item.potency_display);
                if (!dispIsDisc && typeof item.potency_score === 'number' && !isNaN(item.potency_score)) {{
                    hasContinuousPotency = true;
                }}
                if (item.props) {{
                    for (const k of Object.keys(item.props)) {{
                        const v = item.props[k];
                        if (typeof v === 'number' && !isNaN(v)) propCols.add(k);
                    }}
                }}
                if (item.meta) {{
                    for (const k of Object.keys(item.meta)) {{
                        if (k === 'ROMol' || k === 'svg_data' || k === '_potency_sort_score' || k === 'all_endpoints') continue;
                        metaCols.add(k);
                        const rawV = item.meta[k];
                        if (rawV === null || rawV === undefined || rawV === '') continue;
                        const sv = String(rawV).trim();
                        if (!sv) continue;
                        if (_PLAIN_NUM_RE.test(sv)) numericMetaSeen.add(k);
                        if (_parseCensoredCellJs(sv) !== null) censoredMetaSeen.add(k);
                        if (_isConcStringCellJs(sv)) concStrMetaSeen.add(k);
                        if (_DISC_TIER_CELL_RE.test(sv) || /^(\\+{{1,6}}|\\*{{1,6}})$/.test(sv)) tierMetaSeen.add(k);
                        if (/^(?:is|has)_[a-z0-9_]+$/i.test(k) && k.toLowerCase() !== 'name_is_ambiguous') {{
                            let st = boolStats.get(k);
                            if (!st) {{ st = {{ total: 0, ok: true }}; boolStats.set(k, st); }}
                            st.total += 1;
                            if (!/^(?:true|false|0|1|yes|no)$/i.test(sv)) st.ok = false;
                        }}
                    }}
                }}
            }}

            const allCols = Array.from(new Set([...metaCols, ...propCols]));
            const colByNormMeta = new Map();
            for (const c of allCols) {{
                const nm = _normMetaKey(c);
                if (nm && !colByNormMeta.has(nm)) colByNormMeta.set(nm, c);
            }}
            const colCompanionUnits = new Map();
            const molarCompanionCols = new Set();
            for (const c of allCols) {{
                const nm = _normMetaKey(c);
                const unitNorms = (nm === 'activity_value' || nm === 'activity_bin' || nm === 'activity')
                    ? ['activity_units', 'activity_unit', nm + '_units', nm + '_unit']
                    : (nm === 'standard_value'
                        ? ['standard_units', 'standard_unit', nm + '_units', nm + '_unit']
                        : [nm + '_units', nm + '_unit']);
                for (const un of unitNorms) {{
                    if (colByNormMeta.has(un)) {{
                        const uCol = colByNormMeta.get(un);
                        colCompanionUnits.set(c, uCol);
                        if (rawData.some(it => _hasOwnMetaVal(it.meta, uCol) && _normConcUnitJs(it.meta[uCol]) !== null)) {{
                            molarCompanionCols.add(c);
                        }}
                        break;
                    }}
                }}
            }}
            const exactByNorm = new Map();
            const addExact = (normKey, col) => {{
                if (!normKey) return;
                let arr = exactByNorm.get(normKey);
                if (!arr) {{ arr = []; exactByNorm.set(normKey, arr); }}
                if (!arr.includes(col)) arr.push(col);
            }};

            const boolTagsByCol = new Map();
            const boolColsByTag = new Map();
            for (const [col, st] of boolStats.entries()) {{
                if (st.total > 0 && st.ok) {{
                    const tag = col.replace(/^(?:is|has)_/i, '').toLowerCase();
                    boolTagsByCol.set(col, tag);
                    const normTag = _normFieldKey(tag);
                    let arr = boolColsByTag.get(normTag);
                    if (!arr) {{ arr = []; boolColsByTag.set(normTag, arr); }}
                    arr.push(col);
                }}
            }}

            const tierCols = new Set();
            const pScaleCols = new Set();
            const concCols = new Set();
            const concAsPCols = new Set();
            const concColFactorNm = new Map();
            const assayGroups = Object.assign(Object.create(null), {{
                pki: [], pic50: [], pkd: [], pec50: [], pchembl: [],
                ic50: [], ki: [], kd: [], ec50: []
            }});

            const candidateScaleCols = new Set([...propCols, ...censoredMetaSeen]);
            for (const k of candidateScaleCols) {{
                if (tierMetaSeen.has(k) && !numericMetaSeen.has(k) && !concStrMetaSeen.has(k) && !censoredMetaSeen.has(k)) {{
                    tierCols.add(k);
                    continue;
                }}
                const normK = _normalizePotencyColTokensJs(k.toLowerCase());
                if (normK.includes('%') || /(?:^|[^a-z0-9])(?:fold|ratio|selectivity|shift|pct|percent|emax)(?:[^a-z0-9]|$)/i.test(normK)) continue;
                const hasConcSibling = ['_fM', '_fm', '_nM', '_nm', '_uM', '_um', '_pM', '_pm', '_mM', '_mm'].some(sfx => candidateScaleCols.has(k + sfx));
                const hasPlainOrCensNum = !metaCols.has(k) || numericMetaSeen.has(k) || (censoredMetaSeen.has(k) && !propCols.has(k));
                const explicitPName = /(?:^|[^a-z0-9])(?:p(?:ic|ec|ac|gi)50|pki|pkd|pmic|pchembl|pactivity|pact)(?:[^a-z0-9]|$)|(?:^|[^a-z0-9])log(?:10)?[^a-z0-9]*(?:1[^a-z0-9]+)?(?:(?:ic|ec|ac|gi)50|ki|kd|mic)(?:[^a-z0-9]|$)/i.test(normK);
                const isConcAsP = !explicitPName && !numericMetaSeen.has(k) && propCols.has(k) && concStrMetaSeen.has(k) && /(?:^|[^a-z0-9])(?:ic50|ki|kd|ec50|ac50|gi50|mic|potency|activity)(?:[^a-z0-9]|$)/i.test(normK);
                if (explicitPName || isConcAsP) {{
                    pScaleCols.add(k);
                    if (isConcAsP) concAsPCols.add(k);
                    for (const [pKey, baseKey] of [['pki', 'ki'], ['pic50', 'ic50'], ['pkd', 'kd'], ['pec50', 'ec50']]) {{
                        if (new RegExp(`(?:^|[^a-z0-9])(?:${{pKey}}|${{baseKey}})(?:[^a-z0-9]|$)`, 'i').test(normK)) {{
                            assayGroups[pKey].push(k);
                        }}
                    }}
                    if (/(?:^|[^a-z0-9])pchembl/i.test(normK)) assayGroups.pchembl.push(k);
                }} else if (!hasConcSibling && hasPlainOrCensNum) {{
                    const normMetaK = _normMetaKey(k);
                    const concUnitSfx = '(?:(?:f|p|n|u|µ|μ|m)?mol[ ._·/]*l(?:[ ._^-]*1|⁻¹)?|molar|fm|pm|nm|um|µm|μm|mm)';
                    for (const cKey of ['ic50', 'ki', 'kd', 'ec50']) {{
                        if (new RegExp(`(?:^|[^a-z0-9])${{cKey}}(?:[^a-z0-9]|$|${{concUnitSfx}}(?:[^a-z0-9]|$))`, 'i').test(normK)) {{
                            concCols.add(k);
                            concColFactorNm.set(k, _colUnitFactorNm(k));
                            assayGroups[cKey].push(k);
                        }}
                    }}
                    if (!concCols.has(k) && (
                        normMetaK === 'activity_value' ||
                        normMetaK === 'standard_value' ||
                        new RegExp(`(?:^|[^a-z0-9])(?:ac50|gi50|mic)(?:[^a-z0-9]|$|${{concUnitSfx}}(?:[^a-z0-9]|$))`, 'i').test(normK)
                    )) {{
                        concCols.add(k);
                        concColFactorNm.set(k, _colUnitFactorNm(k));
                    }}
                }}
            }}

            const colTokensMap = new Map();
            const numericOnlyCols = new Set();
            for (const k of allCols) {{
                addExact(_normFieldKey(k), k);
                addExact(_normFieldKey(_formatMetaKey(k)), k);
                colTokensMap.set(k, _splitFieldTokens(k));
                if (propCols.has(k) && numericMetaSeen.has(k) && !_ID_LIKE_META_RE.test(k) && !boolTagsByCol.has(k)) {{
                    numericOnlyCols.add(k);
                }}
            }}

            _schemaCacheRef = rawData;
            _schemaCacheLen = rawData.length;
            _schemaCache = {{
                allCols, propCols, metaCols, exactByNorm,
                boolTagsByCol, boolColsByTag, colCompanionUnits, molarCompanionCols,
                tierCols, pScaleCols, concCols, concAsPCols, concColFactorNm, assayGroups,
                colTokensMap, numericOnlyCols, censoredMetaSeen, hasContinuousPotency
            }};
            return _schemaCache;
        }}

        const _BUILTIN_FIELDS = (() => {{
            const map = Object.create(null);
            const add = (names, spec) => {{ for (const n of names) map[n] = spec; }};
            add(['id', 'name', 'compound', 'title'], {{ kind: 'id' }});
            add(['smiles', 'canonicalsmiles', 'smi', 'structure'], {{ kind: 'smiles' }});
            add(['mw', 'molwt', 'molecularweight'], {{ kind: 'item_num', key: 'mw', unitType: 'mw' }});
            add(['clogp', 'logp'], {{ kind: 'item_num', key: 'clogp' }});
            add(['tpsa', 'psa'], {{ kind: 'item_num', key: 'tpsa' }});
            add(['hbd', 'donors', 'hbonddonors'], {{ kind: 'item_num', key: 'hbd' }});
            add(['hba', 'acceptors', 'hbondacceptors'], {{ kind: 'item_num', key: 'hba' }});
            add(['rotb', 'rotbonds', 'rotatablebonds'], {{ kind: 'item_num', key: 'rotb' }});
            add(['fsp3', 'csp3', 'fractioncsp3'], {{ kind: 'item_num', key: 'fsp3' }});
            add(['qed'], {{ kind: 'item_num', key: 'qed' }});
            add(['mpo', 'cnsmpo'], {{ kind: 'item_num', key: 'cns_mpo' }});
            add(['ro5', 'ro5violations', 'lipinski', 'violations'], {{ kind: 'ro5' }});
            add(['pains', 'painsstatus', 'brenk', 'alert', 'alerts'], {{ kind: 'pains' }});
            add(['potency'], {{ kind: 'potency' }});
            return map;
        }})();

        function _resolveFieldSpec(rawField, schema) {{
            if (!_isFieldIdent(rawField)) return null;
            const norm = _normFieldKey(rawField);
            if (!norm) return null;

            if (norm === 'inchikey' || norm === 'standardinchikey') {{
                const ikCols = schema.allCols.filter(_isInchikeyColKey);
                return {{
                    kind: 'cols',
                    cols: ikCols.length > 0 ? ikCols : ['inchikey', 'calculated_inchikey'],
                    unitType: null,
                    fieldHasUnit: false
                }};
            }}

            if (norm === 'pact' || norm === 'pactivity') {{
                if (schema.exactByNorm.has(norm)) {{
                    const cols = schema.exactByNorm.get(norm);
                    const acceptsConc = cols.some(c => schema.concCols.has(c) || schema.concAsPCols.has(c) || schema.pScaleCols.has(c));
                    return {{ kind: 'cols', cols, unitType: acceptsConc ? 'conc' : null, wantsConc: false, wantsPScale: true, fieldHasUnit: false }};
                }}
                return _BUILTIN_FIELDS.potency;
            }}

            if (_BUILTIN_FIELDS[norm]) return _BUILTIN_FIELDS[norm];

            const fieldHasUnitToken = /(?:^|[^a-z0-9µμ]|(?:ic50|ki|kd|ec50|ac50|gi50|mic))(?:(?:f|p|n|u|µ|μ|m)?mol[._-]*l(?:[._-]*1)?|molar|fm|pm|nm|um|µm|μm|mm)(?:[^a-z0-9]|$)|_m$/i.test(rawField);

            if (schema.assayGroups[norm]) {{
                const cols = schema.assayGroups[norm];
                const isConc = (norm === 'ic50' || norm === 'ki' || norm === 'kd' || norm === 'ec50');
                const pFallbackCols = isConc ? (schema.assayGroups['p' + norm] || []) : [];
                const textCols = schema.exactByNorm.get(norm) || [];
                if (cols.length > 0 || pFallbackCols.length > 0 || textCols.length === 0) {{
                    return {{
                        kind: 'assay',
                        assayType: norm,
                        cols,
                        pFallbackCols,
                        textCols,
                        unitType: isConc ? 'conc' : null,
                        fieldHasUnit: false
                    }};
                }}
            }}

            if (schema.exactByNorm.has(norm)) {{
                const cols = schema.exactByNorm.get(norm);
                const acceptsConc = cols.some(c => schema.concCols.has(c) || schema.concAsPCols.has(c) || schema.pScaleCols.has(c));
                const wantsConcExact = cols.some(c => schema.concCols.has(c) || schema.concAsPCols.has(c));
                return {{
                    kind: 'cols',
                    cols,
                    unitType: acceptsConc ? 'conc' : null,
                    wantsConc: wantsConcExact,
                    wantsPScale: !wantsConcExact && cols.some(c => schema.pScaleCols.has(c)),
                    fieldHasUnit: fieldHasUnitToken
                }};
            }}

            const qTokens = _splitFieldTokens(rawField);
            if (qTokens.length === 0) return null;
            if (qTokens.length === 1 && qTokens[0].length <= 2 && !schema.boolColsByTag.has(norm)) {{
                return null;
            }}

            const wantsFold = qTokens.some(t => t === 'fold' || t === 'ratio' || t === 'selectivity' || t === 'shift');
            const wantsConc = !wantsFold && qTokens.some(t => t === 'ic50' || t === 'ki' || t === 'kd' || t === 'ec50' || t === 'ac50' || t === 'gi50' || t === 'mic');
            const wantsPScale = !wantsFold && qTokens.some(t =>
                t === 'pki' || t === 'pic50' || t === 'pkd' || t === 'pec50' ||
                t === 'pchembl' || t === 'pact' || t === 'pactivity'
            );

            const matchedCols = [];
            if (schema.boolColsByTag.has(norm)) {{
                for (const bc of schema.boolColsByTag.get(norm)) matchedCols.push(bc);
            }}
            for (const col of schema.allCols) {{
                if (wantsConc && !schema.concCols.has(col) && !schema.concAsPCols.has(col) && !schema.tierCols.has(col)) continue;
                if (wantsPScale && !schema.pScaleCols.has(col)) continue;
                const cToks = schema.colTokensMap.get(col) || [];
                if (qTokens.every(qt => cToks.includes(qt))) {{
                    if (!matchedCols.includes(col)) matchedCols.push(col);
                }}
            }}

            if (matchedCols.length > 0) {{
                const acceptsConc = wantsConc || matchedCols.some(c => schema.concCols.has(c) || schema.concAsPCols.has(c) || schema.pScaleCols.has(c));
                return {{
                    kind: 'cols',
                    cols: matchedCols,
                    boolCols: schema.boolColsByTag.get(norm) || null,
                    unitType: acceptsConc ? 'conc' : null,
                    wantsConc,
                    wantsPScale,
                    fieldHasUnit: fieldHasUnitToken
                }};
            }}
            return null;
        }}

        function _scanQueryLexemes(input) {{
            const out = [];
            const n = input.length;
            let i = 0;
            while (i < n) {{
                while (i < n && /\\s/.test(input[i])) i++;
                if (i >= n) break;
                const start = i;
                let text = '';
                let unquotedChars = 0;
                let firstCharQuoted = false;
                let hasQuoted = false;
                while (i < n && !/\\s/.test(input[i])) {{
                    const ch = input[i];
                    if (ch === '"') {{
                        if (i === start) firstCharQuoted = true;
                        hasQuoted = true;
                        const q = ch;
                        i++;
                        while (i < n && input[i] !== q) {{
                            if (input[i] === '\\\\' && i + 1 < n) {{
                                const nxt = input[i + 1];
                                if (nxt === '\\\\' || (nxt === q && input.indexOf(q, i + 2) !== -1)) i++;
                            }}
                            text += input[i];
                            i++;
                        }}
                        if (i < n && input[i] === q) i++;
                    }} else {{
                        text += ch;
                        unquotedChars++;
                        i++;
                    }}
                }}
                const wholeQuoted = firstCharQuoted && unquotedChars === 0;
                if (text.length > 0 || wholeQuoted) {{
                    out.push({{ text, wholeQuoted, hasQuoted }});
                }}
            }}
            return out;
        }}

        function _applyUnitMultiplier(val, unit, unitType) {{
            if (!unit) return val;
            const rawU = String(unit).trim();
            if (!rawU) return val;
            const uLow = rawU.toLowerCase();
            if (uLow === 'da' || uLow === 'g/mol') return unitType === 'conc' ? null : val;
            const normU = _normConcUnitJs(rawU);
            if (normU !== null) return unitType === 'conc' ? val * _CONC_FACTOR_NM[normU] : null;
            return null;
        }}

        function _parseNumericWithUnit(rhsStr, unitType) {{
            const s = String(rhsStr ?? '').trim();
            const m = s.match(/^([+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:e[+-]?\\d+)?)(?:\\s*(.+))?$/);
            if (!m) return null;
            const rawVal = Number(m[1]);
            if (isNaN(rawVal)) return null;
            const scaled = _applyUnitMultiplier(rawVal, m[2], unitType);
            if (scaled === null) return null;
            return {{ val: scaled, hasUnit: Boolean(m[2]), rawNum: m[1] }};
        }}

        function _parseRangeWithUnit(rhsStr, unitType) {{
            const parts = String(rhsStr ?? '').trim().split('..');
            if (parts.length !== 2) return null;
            const m1 = parts[0].trim().match(/^([+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:e[+-]?\\d+)?)(?:\\s*(.+))?$/);
            const m2 = parts[1].trim().match(/^([+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:e[+-]?\\d+)?)(?:\\s*(.+))?$/);
            if (!m1 || !m2) return null;
            const u2 = m2[2] || '';
            const u1 = m1[2] || u2;
            const v1 = _applyUnitMultiplier(Number(m1[1]), u1, unitType);
            const v2 = _applyUnitMultiplier(Number(m2[1]), u2, unitType);
            if (v1 === null || v2 === null || isNaN(v1) || isNaN(v2)) return null;
            return {{ bounds: [Math.min(v1, v2), Math.max(v1, v2)], hasUnit: Boolean(u1 || u2) }};
        }}

        function _parseQueryBranches(rawQuery, schema) {{
            const lexemes = _scanQueryLexemes(rawQuery);
            if (lexemes.length === 0) return [];

            const rawBranches = [[]];
            for (const lex of lexemes) {{
                if (!lex.wholeQuoted && (lex.text.toUpperCase() === 'OR' || lex.text === '||')) {{
                    if (rawBranches[rawBranches.length - 1].length > 0) rawBranches.push([]);
                    continue;
                }}
                rawBranches[rawBranches.length - 1].push(lex);
            }}

            const isConj = (lx) => !lx.wholeQuoted && (lx.text.toUpperCase() === 'AND' || lx.text === '&&');
            const canStitchUnit = (rhs) => _PLAIN_NUM_RE.test(rhs) || /\\.\\.[+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:e[+-]?\\d+)?$/i.test(rhs);

            const branches = [];
            for (const bLex of rawBranches) {{
                if (bLex.length === 0) continue;
                const clauses = [];
                let i = 0;
                while (i < bLex.length) {{
                    const cur = bLex[i];
                    if (isConj(cur)) {{ i++; continue; }}
                    if (cur.wholeQuoted) {{
                        if (cur.text) clauses.push({{ type: 'text', raw: cur.text, lower: cur.text.toLowerCase(), wholeQuoted: true, negate: false }});
                        i++;
                        continue;
                    }}
                    if (cur.text === '-' || cur.text === '!' || _OP_LIST.includes(cur.text)) {{
                        i++;
                        continue;
                    }}

                    let rawTok = cur.text;
                    let negate = false;
                    if ((rawTok.startsWith('-') && rawTok.length > 1 && !/^-\\d/.test(rawTok) && !/^-\\.\\d/.test(rawTok)) ||
                        (rawTok.startsWith('!') && rawTok.length > 1 && !rawTok.startsWith('!='))) {{
                        negate = true;
                        rawTok = rawTok.slice(1);
                    }}

                    if (i + 1 < bLex.length && !bLex[i + 1].wholeQuoted && !isConj(bLex[i + 1])) {{
                        const lw0 = rawTok.toLowerCase();
                        const lw1 = bLex[i + 1].text.toLowerCase();
                        if ((lw0 === 'pains' || lw0 === 'brenk') && (lw1 === 'clean' || lw1 === 'pass' || lw1 === 'alert' || lw1 === 'alerts' || lw1 === 'fail')) {{
                            clauses.push({{ type: 'field', spec: _BUILTIN_FIELDS.pains, rawField: lw0, op: ':', rhs: lw1, rhsQuoted: false, negate }});
                            i += 2;
                            continue;
                        }}
                        if ((lw0 === 'ro5' || lw0 === 'lipinski') && (lw1 === 'pass' || lw1 === 'clean' || lw1 === 'fail' || lw1 === 'viol' || lw1 === 'alert')) {{
                            clauses.push({{ type: 'field', spec: _BUILTIN_FIELDS.ro5, rawField: lw0, op: ':', rhs: lw1, rhsQuoted: false, negate }});
                            i += 2;
                            continue;
                        }}
                        if (lw0 === 'col' && /^\\d+$/.test(lw1) && !(i + 2 < bLex.length && !bLex[i + 2].wholeQuoted && _OP_LIST.some(op => bLex[i + 2].text === op || bLex[i + 2].text.startsWith(op)))) {{
                            clauses.push({{ type: 'text', raw: `${{rawTok}} ${{bLex[i + 1].text}}`, lower: `col ${{lw1}}`, wholeQuoted: false, negate }});
                            i += 2;
                            continue;
                        }}
                    }}

                    let candLhs = null, candOp = null, candRhs = '', candQuoted = false, advance = 1;
                    let opPos = -1;
                    for (const op of _OP_LIST) {{
                        const idx = rawTok.indexOf(op);
                        if (idx > 0 && (opPos === -1 || idx < opPos || (idx === opPos && op.length > candOp.length))) {{
                            opPos = idx;
                            candOp = op;
                        }}
                    }}
                    if (opPos > 0) {{
                        const looksSmilesDoubleBond = !cur.hasQuoted && candOp === '=' && (rawTok.indexOf('=', opPos + 1) !== -1 || /[()\\[\\]#\\/\\\\@]/.test(rawTok));
                        if (!looksSmilesDoubleBond) {{
                            candLhs = rawTok.slice(0, opPos);
                            candRhs = rawTok.slice(opPos + candOp.length);
                            candQuoted = cur.hasQuoted;
                            if (candRhs === '' && i + 1 < bLex.length && !isConj(bLex[i + 1])) {{
                                const specTry = _resolveFieldSpec(candLhs, schema);
                                if (specTry || (candOp !== ':' && _isFieldIdent(candLhs))) {{
                                    candRhs = bLex[i + 1].text;
                                    candQuoted = bLex[i + 1].wholeQuoted || bLex[i + 1].hasQuoted;
                                    advance = 2;
                                }}
                            }}
                        }} else {{
                            candOp = null;
                        }}
                    }} else if (i + 1 < bLex.length && !bLex[i + 1].wholeQuoted && !isConj(bLex[i + 1]) && _isFieldIdent(rawTok)) {{
                        const nextTxt = bLex[i + 1].text;
                        const matchedOp = _OP_LIST.find(op => nextTxt === op || nextTxt.startsWith(op));
                        if (matchedOp) {{
                            candLhs = rawTok;
                            candOp = matchedOp;
                            if (nextTxt === matchedOp) {{
                                advance = 2;
                                if (i + 2 < bLex.length && !isConj(bLex[i + 2])) {{
                                    candRhs = bLex[i + 2].text;
                                    candQuoted = bLex[i + 2].wholeQuoted || bLex[i + 2].hasQuoted;
                                    advance = 3;
                                }}
                            }} else {{
                                candRhs = nextTxt.slice(matchedOp.length);
                                candQuoted = bLex[i + 1].hasQuoted;
                                advance = 2;
                            }}
                        }}
                    }}

                    if (candLhs && candOp && _isFieldIdent(candLhs)) {{
                        const spec = _resolveFieldSpec(candLhs, schema);
                        if (candRhs === '') {{
                            if (spec || candOp === ':') {{
                                i += advance;
                                continue;
                            }}
                        }} else {{
                            let rhsFull = candRhs;
                            let advFull = advance;
                            const isUnitLex = (t) => _normConcUnitJs(t) !== null || /^(?:da|g\\/mol)$/i.test(t);
                            if (!candQuoted && i + advFull < bLex.length && !bLex[i + advFull].wholeQuoted &&
                                isUnitLex(bLex[i + advFull].text) && canStitchUnit(rhsFull)) {{
                                rhsFull += ' ' + bLex[i + advFull].text;
                                advFull += 1;
                            }}
                            const uType = spec ? spec.unitType : (/(?:ic50|ki|kd|ec50)/i.test(candLhs) ? 'conc' : null);
                            const rhsLooksNumOrRange = (
                                _parseNumericWithUnit(rhsFull, uType) !== null ||
                                _parseRangeWithUnit(rhsFull, uType) !== null
                            );
                            const validUnknownNumCmp = !spec && candOp !== ':' && candOp !== '=' && candOp !== '==' && rhsLooksNumOrRange && _isFieldIdent(candLhs);
                            if (spec || validUnknownNumCmp) {{
                                clauses.push({{ type: 'field', spec, rawField: candLhs, op: candOp, rhs: rhsFull, rhsQuoted: candQuoted, negate }});
                                i += advFull;
                                continue;
                            }}
                        }}
                    }}

                    clauses.push({{ type: 'text', raw: rawTok, lower: rawTok.toLowerCase(), wholeQuoted: false, negate }});
                    i++;
                }}
                if (clauses.length > 0) branches.push(clauses);
            }}
            return branches;
        }}

        function _compareNum(val, op, target, rhsRaw = '', tol = 1e-6) {{
            if (val === null || val === undefined || typeof val !== 'number' || isNaN(val)) return false;
            if (op === '<') return val < target;
            if (op === '<=') return val <= target + tol;
            if (op === '>') return val > target;
            if (op === '>=') return val >= target - tol;
            if (op === ':' && rhsRaw && /^[+-]?\\d+\\.\\d+$/.test(rhsRaw)) {{
                if (Math.abs(val - target) <= tol) return true;
                const normRhs = rhsRaw.startsWith('+') ? rhsRaw.slice(1) : rhsRaw;
                return String(val).startsWith(normRhs);
            }}
            if (op === '=' || op === '==' || op === ':') return Math.abs(val - target) <= tol;
            if (op === '!=') return Math.abs(val - target) > tol;
            return false;
        }}

        function _matchNumOrRange(val, op, numInfo, rangeInfo, rhsRaw = '') {{
            if (val === null || val === undefined || typeof val !== 'number' || isNaN(val)) return false;
            if (rangeInfo) {{
                const inRange = val >= rangeInfo.bounds[0] - 1e-6 && val <= rangeInfo.bounds[1] + 1e-6;
                return op === '!=' ? !inRange : inRange;
            }}
            return numInfo !== null && _compareNum(val, op, numInfo.val, rhsRaw);
        }}

        function _evalConcOnPScale(pVal, op, numInfo, rangeInfo) {{
            if (pVal === null || pVal === undefined || typeof pVal !== 'number' || isNaN(pVal)) return false;
            const tol = 1.5e-3;
            if (rangeInfo) {{
                const [loNm, hiNm] = rangeInfo.bounds;
                if (hiNm <= 0) return op === '!=';
                const pMin = 9 - Math.log10(hiNm);
                const pMax = loNm > 0 ? (9 - Math.log10(loNm)) : Infinity;
                const inRange = pVal >= pMin - tol && pVal <= pMax + tol;
                return op === '!=' ? !inRange : inRange;
            }}
            if (!numInfo) return false;
            const cNm = numInfo.val;
            if (cNm <= 0) return op === '>' || op === '>=' || op === '!=';
            const pTarget = 9 - Math.log10(cNm);
            if (op === '<') return pVal > pTarget;
            if (op === '<=') return pVal >= pTarget - tol;
            if (op === '>') return pVal < pTarget;
            if (op === '>=') return pVal <= pTarget + tol;
            if (op === '!=') return Math.abs(pVal - pTarget) > tol;
            return Math.abs(pVal - pTarget) <= tol;
        }}

        function _parseBoolLiteral(s) {{
            const lw = String(s ?? '').trim().toLowerCase();
            if (lw === 'true' || lw === 'yes') return true;
            if (lw === 'false' || lw === 'no') return false;
            return null;
        }}

        function _evalCensoredInterval(isUpper, cStrict, boundVal, effOp, numInfo, rangeInfo) {{
            const tol = 1e-6;
            if (rangeInfo) {{
                if (effOp !== '!=') return false;
                const [lo, hi] = rangeInfo.bounds;
                return isUpper
                    ? (cStrict ? boundVal <= lo + tol : boundVal < lo - tol)
                    : (cStrict ? boundVal >= hi - tol : boundVal > hi + tol);
            }}
            if (!numInfo) return false;
            const target = numInfo.val;
            if (isUpper) {{
                if (effOp === '<') return cStrict ? boundVal <= target + tol : boundVal < target - tol;
                if (effOp === '<=') return boundVal <= target + tol;
                if (effOp === '!=') return cStrict ? boundVal <= target + tol : boundVal < target - tol;
                return false;
            }}
            if (effOp === '>') return cStrict ? boundVal >= target - tol : boundVal > target + tol;
            if (effOp === '>=') return boundVal >= target - tol;
            if (effOp === '!=') return cStrict ? boundVal >= target - tol : boundVal > target + tol;
            return false;
        }}

        function _evalFieldClause(item, clause, schema) {{
            const {{ spec, op, rhs, rhsQuoted }} = clause;
            if (!spec) return false;
            const rhsTrim = String(rhs ?? '').trim();
            const rhsLower = rhsTrim.toLowerCase();

            const rangeInfo = _parseRangeWithUnit(rhsTrim, spec.unitType);
            const numInfo = rangeInfo ? null : _parseNumericWithUnit(rhsTrim, spec.unitType);
            const hasExplicitUnit = Boolean((rangeInfo && rangeInfo.hasUnit) || (numInfo && numInfo.hasUnit));
            const isIneqOp = (op === '<' || op === '<=' || op === '>' || op === '>=');

            if (spec.kind === 'id') {{
                const candidates = [item.id];
                const seenCols = new Set();
                const pushCol = (ck) => {{
                    if (!_hasOwn(item.meta, ck) || seenCols.has(ck)) return;
                    seenCols.add(ck);
                    const mv = item.meta[ck];
                    if (mv !== null && mv !== undefined && mv !== '') candidates.push(String(mv));
                }};
                const exactCols = schema.exactByNorm.get(_normFieldKey(clause.rawField)) || [];
                for (const ck of exactCols) pushCol(ck);
                if (item.meta) {{
                    for (const ck of Object.keys(item.meta)) {{
                        if (/(?:^|[\\W_])(?:id|name|title|compound|synonyms?)(?:$|[\\W_])/i.test(ck) && !/(?:^|[\\W_])ambiguous(?:$|[\\W_])/i.test(ck)) {{
                            pushCol(ck);
                        }}
                    }}
                }}
                if (op === '=' || op === '==') return candidates.some(c => String(c ?? '').toLowerCase() === rhsLower);
                if (op === '!=') return candidates.every(c => String(c ?? '').toLowerCase() !== rhsLower);
                return candidates.some(c => _textMatchesToken(String(c ?? '').toLowerCase(), rhsLower, Boolean(rhsQuoted)));
            }}

            if (spec.kind === 'smiles') {{
                const smi = String(item.smiles ?? '');
                const matched = rhsTrim.length === 1 && /[A-Z]/.test(rhsTrim)
                    ? new RegExp(rhsTrim.replace(/[.*+?^${{}}()|[\\]\\\\]/g, '\\\\$&') + '(?![a-z])').test(smi)
                    : smi.includes(rhsTrim);
                return op === '!=' ? !matched : matched;
            }}

            if (spec.kind === 'item_num') {{
                return _matchNumOrRange(item[spec.key], op, numInfo, rangeInfo, rhsTrim);
            }}

            if (spec.kind === 'ro5') {{
                const v = item.ro5_violations;
                if (rangeInfo || numInfo) return _matchNumOrRange(v, op, numInfo, rangeInfo, rhsTrim);
                if (rhsLower === 'pass' || rhsLower === 'clean' || rhsLower === 'ok' || rhsLower === 'true') {{
                    return op === '!=' ? (v !== 0) : (v === 0);
                }}
                if (rhsLower === 'fail' || rhsLower === 'alert' || rhsLower === 'viol' || rhsLower === 'false') {{
                    return op === '!=' ? (v === 0) : (v !== null && v !== undefined && v > 0);
                }}
                const badgeTxt = v === 0 ? 'ro5 pass' : (v != null ? `ro5: ${{v}} viol` : '');
                const matched = badgeTxt.includes(rhsLower);
                return op === '!=' ? !matched : matched;
            }}

            if (spec.kind === 'pains') {{
                const lowP = String(item.pains ?? 'Clean').toLowerCase();
                const isClean = lowP.includes('clean');
                if (rhsLower === 'clean' || rhsLower === 'pass' || rhsLower === 'none' || rhsLower === 'false') {{
                    return op === '!=' ? !isClean : isClean;
                }}
                if (rhsLower === 'alert' || rhsLower === 'alerts' || rhsLower === 'fail' || rhsLower === 'true') {{
                    return op === '!=' ? isClean : !isClean;
                }}
                const matched = lowP.includes(rhsLower);
                return op === '!=' ? !matched : matched;
            }}

            if (spec.kind === 'potency') {{
                const pVals = [];
                const dispIsDisc = item.potency_display && _DISC_POTENCY_DISP_RE.test(item.potency_display);
                if (schema.hasContinuousPotency) {{
                    if (!dispIsDisc && typeof item.potency_score === 'number' && !isNaN(item.potency_score)) {{
                        pVals.push(item.potency_score);
                    }}
                }} else if (item.props) {{
                    for (const pk of schema.pScaleCols) {{
                        if (typeof item.props[pk] === 'number' && !isNaN(item.props[pk])) pVals.push(item.props[pk]);
                    }}
                }}
                if (rangeInfo || numInfo) {{
                    if (op === '!=') return pVals.length > 0 && pVals.every(v => _matchNumOrRange(v, '!=', numInfo, rangeInfo, rhsTrim));
                    return pVals.some(v => _matchNumOrRange(v, op, numInfo, rangeInfo, rhsTrim));
                }}
                const disp = String(item.potency_display ?? '').toLowerCase();
                const ok = _textMatchesToken(disp, rhsLower, Boolean(rhsQuoted));
                return op === '!=' ? !ok : ok;
            }}

            let cols = spec.cols || [];
            if (spec.kind === 'assay') {{
                if (rangeInfo || numInfo) {{
                    if (cols.length === 0 && spec.pFallbackCols && spec.pFallbackCols.length > 0) {{
                        cols = spec.pFallbackCols;
                    }}
                }} else {{
                    cols = Array.from(new Set([...cols, ...(spec.pFallbackCols || []), ...(spec.textCols || [])]));
                }}
            }} else if (spec.kind === 'cols' && (rangeInfo || numInfo) && !hasExplicitUnit && !spec.wantsConc && !spec.wantsPScale && !spec.fieldHasUnit) {{
                const isTruePCol = (c) => schema.pScaleCols.has(c) && !schema.concAsPCols.has(c);
                const isConcSemCol = (c) => schema.concCols.has(c) || schema.concAsPCols.has(c);
                if (cols.some(isTruePCol) && cols.some(isConcSemCol)) {{
                    const preferConc = Boolean(_ACTIVE_COL_HINT && cols.includes(_ACTIVE_COL_HINT) && isConcSemCol(_ACTIVE_COL_HINT));
                    cols = cols.filter(preferConc ? isConcSemCol : isTruePCol);
                }}
            }}
            const boolLit = _parseBoolLiteral(rhsTrim);
            if (boolLit !== null && spec.boolCols && spec.boolCols.length > 0) {{
                cols = spec.boolCols;
            }}

            const evalColNum = (col, effOp) => {{
                if (schema.tierCols.has(col)) return null;
                const resolveColFactorNm = () => {{
                    const uCol = schema.colCompanionUnits.get(col);
                    if (uCol && _hasOwnMetaVal(item.meta, uCol)) {{
                        const cu = _normConcUnitJs(item.meta[uCol]);
                        if (cu !== null) return _CONC_FACTOR_NM[cu];
                        return (hasExplicitUnit || schema.molarCompanionCols.has(col)) ? null : 1;
                    }}
                    return schema.concColFactorNm.get(col) || 1;
                }};
                const useConcOnP = spec.unitType === 'conc' && (hasExplicitUnit || !spec.wantsPScale) && (
                    schema.concAsPCols.has(col) || (
                        schema.pScaleCols.has(col) && (
                            spec.kind === 'assay' ||
                            Boolean(spec.wantsConc || spec.fieldHasUnit || hasExplicitUnit)
                        )
                    )
                );
                if (_hasOwn(item.meta, col) && item.meta[col] !== undefined && item.meta[col] !== null) {{
                    const cens = _parseCensoredCellJs(item.meta[col]);
                    if (cens) {{
                        const cOpRaw = cens.op;
                        let isUpper = (cOpRaw === '<' || cOpRaw === '<=' || cOpRaw === '≤');
                        const cStrict = (cOpRaw === '<' || cOpRaw === '>');
                        const rawBound = cens.bound;
                        let boundScaled = rawBound;
                        if (cens.concUnit !== null) {{
                            const boundNm = rawBound * _CONC_FACTOR_NM[cens.concUnit];
                            if (spec.unitType === 'conc' && (useConcOnP || schema.concCols.has(col))) {{
                                if (hasExplicitUnit || !spec.fieldHasUnit) {{
                                    boundScaled = boundNm;
                                }} else {{
                                    const cf = resolveColFactorNm();
                                    if (cf === null) return null;
                                    boundScaled = boundNm / cf;
                                }}
                            }} else if (schema.pScaleCols.has(col)) {{
                                if (boundNm <= 0) return false;
                                boundScaled = 9 - Math.log10(boundNm);
                                isUpper = !isUpper;
                            }}
                        }} else if (useConcOnP && !schema.concAsPCols.has(col)) {{
                            boundScaled = Math.pow(10, 9 - rawBound);
                            isUpper = !isUpper;
                        }} else if (spec.unitType === 'conc' && schema.concCols.has(col)) {{
                            const cf = resolveColFactorNm();
                            if (cf === null) return null;
                            if (hasExplicitUnit || !spec.fieldHasUnit) {{
                                boundScaled = rawBound * cf;
                            }}
                        }}
                        return _evalCensoredInterval(isUpper, cStrict, boundScaled, effOp, numInfo, rangeInfo);
                    }}
                }}
                const nv = getCompoundValue(item, col);
                if (nv === null || isNaN(nv)) return null;
                if (useConcOnP) {{
                    return _evalConcOnPScale(nv, effOp, numInfo, rangeInfo);
                }}
                let scaledVal = nv;
                if (spec.unitType === 'conc' && schema.concCols.has(col)) {{
                    const cf = resolveColFactorNm();
                    if (cf === null) return null;
                    if (hasExplicitUnit || !spec.fieldHasUnit) {{
                        scaledVal = nv * cf;
                    }}
                }}
                return _matchNumOrRange(scaledVal, effOp, numInfo, rangeInfo, rhsTrim);
            }};

            const evalColPositive = (col, effOp) => {{
                if (boolLit !== null && schema.boolTagsByCol.has(col)) {{
                    const mv = _hasOwn(item.meta, col) ? item.meta[col] : undefined;
                    if (mv === null || mv === undefined || mv === '') return false;
                    return /^(?:true|1|yes)$/i.test(String(mv).trim()) === boolLit;
                }}
                if (rangeInfo) {{
                    return evalColNum(col, effOp) === true;
                }}
                if (numInfo !== null && (isIneqOp || (!schema.tierCols.has(col) && ((_hasOwn(item.props, col) && item.props[col] !== undefined) || schema.censoredMetaSeen.has(col))))) {{
                    const res = evalColNum(col, effOp);
                    if (res !== null) return res;
                    if (isIneqOp) return false;
                }}
                if (isIneqOp) return false;
                const mv = _hasOwn(item.meta, col) ? item.meta[col] : undefined;
                if (mv === null || mv === undefined || mv === '') return false;
                const rawMv = String(mv).trim();
                if (schema.tierCols.has(col) && /^[a-f]$/i.test(rhsTrim)) {{
                    const mTier = rawMv.match(_DISC_TIER_CELL_RE);
                    return Boolean(mTier && mTier[1].toLowerCase() === rhsLower);
                }}
                const sv = rawMv.toLowerCase();
                const fv = _formatMetaVal(mv, col).toLowerCase();
                if (effOp === '=' || effOp === '==') return sv === rhsLower || fv === rhsLower;
                return _textMatchesToken(sv, rhsLower, Boolean(rhsQuoted)) || _textMatchesToken(fv, rhsLower, Boolean(rhsQuoted));
            }};

            if (spec.kind === 'assay' && cols.length === 0) {{
                const disp = String(item.potency_display ?? '');
                const dispNorm = _normalizePotencyColTokensJs(disp.toLowerCase());
                const baseAssay = spec.assayType.replace(/^p/, '');
                const hasAssayTok = new RegExp(`(?:^|[^a-z0-9])p?${{baseAssay}}(?:[^a-z0-9]|$)`, 'i').test(dispNorm);
                if (disp && !_DISC_POTENCY_DISP_RE.test(disp) && hasAssayTok) {{
                    if (rangeInfo || numInfo) {{
                        return spec.unitType === 'conc'
                            ? _evalConcOnPScale(item.potency_score, op, numInfo, rangeInfo)
                            : _matchNumOrRange(item.potency_score, op, numInfo, rangeInfo, rhsTrim);
                    }}
                }}
                return false;
            }}

            if (op === '!=') {{
                if (boolLit !== null && cols.some(c => schema.boolTagsByCol.has(c))) {{
                    return cols.some(col => {{
                        if (!schema.boolTagsByCol.has(col)) return false;
                        const mv = _hasOwn(item.meta, col) ? item.meta[col] : undefined;
                        if (mv === null || mv === undefined || mv === '') return false;
                        return /^(?:true|1|yes)$/i.test(String(mv).trim()) !== boolLit;
                    }});
                }}
                if (rangeInfo || numInfo) {{
                    const evals = cols.filter(c => !schema.tierCols.has(c)).map(c => evalColNum(c, '!=')).filter(r => r !== null);
                    if (evals.length > 0) return evals.every(Boolean);
                }}
                return !cols.some(col => evalColPositive(col, '='));
            }}
            return cols.some(col => evalColPositive(col, op));
        }}

        function _badgeWordMatch(word, queryLower) {{
            if (queryLower === word) return true;
            return queryLower.length >= 4 && word.startsWith(queryLower);
        }}

        function _textMatchesToken(haystackLower, queryLower, wholeQuoted) {{
            if (!haystackLower) return false;
            if (wholeQuoted || queryLower.length >= 3) return haystackLower.includes(queryLower);
            const esc = queryLower.replace(/[.*+?^${{}}()|[\\]\\\\]/g, '\\\\$&');
            const isDigitQuery = /^\\d/.test(queryLower);
            const lead = isDigitQuery ? '(?:^|[^a-z0-9.]|(?:^|\\\\D)\\\\.)' : '(?:^|[^a-z0-9])';
            const trail = isDigitQuery ? '(?:$|[^a-z0-9.]|\\\\.(?!\\\\d))' : '(?:$|[^a-z0-9])';
            const pat = queryLower.length === 1
                ? `${{lead}}${{esc}}${{trail}}`
                : `${{lead}}${{esc}}`;
            return new RegExp(pat).test(haystackLower);
        }}

        function _matchesUnscopedToken(item, clause, schema) {{
            const {{ raw, lower, wholeQuoted }} = clause;
            if (!raw) return true;

            if (!wholeQuoted && (raw === 'Cl' || raw === 'Br')) {{
                if (String(item.smiles ?? '').includes(raw)) return true;
                if (item.meta) {{
                    for (const fk of Object.keys(item.meta)) {{
                        if (_isFormulaColKey(fk) && String(item.meta[fk] ?? '').includes(raw)) return true;
                    }}
                }}
                return false;
            }}

            const idLow = String(item.id ?? '').toLowerCase();
            if ((wholeQuoted || lower.length >= 2) ? idLow.includes(lower) : _textMatchesToken(idLow, lower, false)) return true;

            const smi = String(item.smiles ?? '');
            if (smi) {{
                if (/[A-Z]/.test(raw)) {{
                    if (raw.length === 1) {{
                        const esc = raw.replace(/[.*+?^${{}}()|[\\]\\\\]/g, '\\\\$&');
                        if (new RegExp(esc + '(?![a-z])').test(smi)) return true;
                    }} else if (smi.includes(raw)) {{
                        return true;
                    }}
                }} else if ((!/^[+-]?\\d+(?:\\.\\d+)?$/.test(raw) && /[0-9=#@\\[\\]()\\/\\\\+-]/.test(raw)) || (raw.length >= 3 && /^[cnospb]+$/.test(raw))) {{
                    if (smi.includes(raw)) return true;
                }}
            }}

            if (item.potency_display && (wholeQuoted || (!/^\\d+$/.test(raw) && (raw.length > 1 || raw === lower)))) {{
                if (_textMatchesToken(String(item.potency_display).toLowerCase(), lower, wholeQuoted)) return true;
            }}

            if (item.ro5_violations === 0) {{
                if (wholeQuoted ? 'ro5 pass'.includes(lower) : (lower === 'ro5 pass' || lower === 'ro5' || _badgeWordMatch('pass', lower))) {{
                    return true;
                }}
            }} else if (item.ro5_violations !== null && item.ro5_violations !== undefined && item.ro5_violations > 0) {{
                const ro5Badge = `ro5: ${{item.ro5_violations}} viol`;
                if (wholeQuoted ? ro5Badge.includes(lower) : (ro5Badge === lower || _badgeWordMatch('viol', lower))) {{
                    return true;
                }}
            }}

            if (item.pains) {{
                const painsRaw = String(item.pains).toLowerCase();
                const isClean = painsRaw.includes('clean');
                if (isClean) {{
                    if (wholeQuoted ? 'pains clean'.includes(lower) : (lower === 'pains clean' || _badgeWordMatch('clean', lower))) {{
                        return true;
                    }}
                }} else {{
                    if (((wholeQuoted || !/^\\d+$/.test(raw)) && _textMatchesToken(painsRaw, lower, wholeQuoted)) || lower === 'pains alert' || _badgeWordMatch('pains', lower) || _badgeWordMatch('alert', lower)) {{
                        return true;
                    }}
                }}
            }}
            if (wholeQuoted && item.cns_mpo !== null && item.cns_mpo !== undefined && !isNaN(item.cns_mpo) && item.cns_mpo >= 0) {{
                if (`mpo: ${{item.cns_mpo}}`.includes(lower)) return true;
            }}

            for (const [col, tag] of schema.boolTagsByCol.entries()) {{
                const mv = _hasOwn(item.meta, col) ? item.meta[col] : undefined;
                if (mv !== null && mv !== undefined && /^(?:true|1|yes)$/i.test(String(mv).trim())) {{
                    const tagSpaced = tag.replace(/_/g, ' ');
                    if (tag === lower || tagSpaced === lower || (lower.length >= 4 && (tag.includes(lower) || tagSpaced.includes(lower)))) {{
                        return true;
                    }}
                }}
            }}

            if (item.meta) {{
                const allowInchikey = raw.length >= 6 && _INCHIKEY_QUERY_RE.test(raw);
                const anonCovered = _areAnonColsCovered(item.meta);
                const anonColIdxMatch = /^col[\\s_]+(\\d+)$/i.exec(lower);
                for (const k of Object.keys(item.meta)) {{
                    if (_UNSCOPED_META_SKIP.has(k) || _isSmilesColKey(k) || schema.boolTagsByCol.has(k)) continue;
                    const anonKeyMatch = _ANON_PATENT_COL_RE.exec(k);
                    if (anonKeyMatch) {{
                        if (!_isVisibleMetaKey(item.meta, k, _isModalSkipKey, anonCovered)) continue;
                        if (anonColIdxMatch && anonKeyMatch[1] === anonColIdxMatch[1]) return true;
                    }}
                    if (_isInchikeyColKey(k)) {{
                        if (allowInchikey && String(item.meta[k] ?? '').toLowerCase().includes(lower)) return true;
                        continue;
                    }}
                    if (_isFormulaColKey(k)) {{
                        if (/[A-Z]/.test(raw) && raw.length >= 2 && String(item.meta[k] ?? '').includes(raw)) return true;
                        continue;
                    }}
                    const mv = item.meta[k];
                    if (mv === null || mv === undefined || mv === '') continue;
                    const sv = String(mv).trim();
                    if (!sv) continue;
                    const svLow = sv.toLowerCase();
                    if (svLow === 'true' || svLow === 'false') continue;
                    if (k === 'PAINS_Status' && (svLow === 'clean' || (!wholeQuoted && /^\\d+$/.test(raw)))) continue;
                    if (schema.numericOnlyCols.has(k)) {{
                        if (wholeQuoted || lower.length >= 3) {{
                            const fmtLow = _formatMetaVal(mv, k).toLowerCase();
                            if (svLow === lower || fmtLow === lower) {{
                                return true;
                            }}
                            if (_PLAIN_NUM_RE.test(lower) && !/^[+-]?0{{2,}}$/.test(lower) && _PLAIN_NUM_RE.test(svLow) && Number(svLow) === Number(lower)) {{
                                return true;
                            }}
                            if (wholeQuoted || !_PLAIN_NUM_RE.test(lower)) {{
                                const entryLow = _formatMetaEntryVal(item.meta, k).toLowerCase();
                                if (entryLow !== fmtLow && _textMatchesToken(entryLow, lower, wholeQuoted)) {{
                                    return true;
                                }}
                            }}
                        }}
                        continue;
                    }}
                    if (!wholeQuoted && raw.length === 1 && /[A-Z]/.test(raw)) continue;
                    if (
                        _textMatchesToken(svLow, lower, wholeQuoted) ||
                        _textMatchesToken(_formatMetaVal(mv, k).toLowerCase(), lower, wholeQuoted) ||
                        _textMatchesToken(_formatMetaEntryVal(item.meta, k).toLowerCase(), lower, wholeQuoted)
                    ) {{
                        return true;
                    }}
                }}
            }}
            return false;
        }}

        function _evalQueryBranches(item, branches, schema) {{
            if (branches.length === 0) return true;
            return branches.some(clauses => clauses.every(clause => {{
                const matched = clause.type === 'field'
                    ? _evalFieldClause(item, clause, schema)
                    : _matchesUnscopedToken(item, clause, schema);
                return clause.negate ? !matched : matched;
            }}));
        }}

        function filterGrid() {{
            const searchEl = document.getElementById("searchInput");
            const rawQuery = searchEl ? searchEl.value.trim() : "";
            if (!rawQuery) {{
                currentData = [...rawData];
            }} else {{
                const schema = _getSchemaIndex();
                const branches = _parseQueryBranches(rawQuery, schema);
                currentData = branches.length === 0
                    ? [...rawData]
                    : rawData.filter(item => _evalQueryBranches(item, branches, schema));
            }}
            sortGrid(true);
        }}

        function _cmpNum(va, vb, asc = true) {{
            const na = (va === null || va === undefined || isNaN(va));
            const nb = (vb === null || vb === undefined || isNaN(vb));
            if (na && nb) return 0;
            if (na) return 1;
            if (nb) return -1;
            return asc ? (va - vb) : (vb - va);
        }}

        const _SORT_SPECS = Object.assign(Object.create(null), {{
            potency_desc: ['potency_score', false],
            potency_asc: ['potency_score', true],
            mw_asc: ['mw', true],
            mw_desc: ['mw', false],
            clogp_asc: ['clogp', true],
            clogp_desc: ['clogp', false],
            cns_mpo_desc: ['cns_mpo', false],
            qed_desc: ['qed', false],
            tpsa_asc: ['tpsa', true],
        }});

        function sortGrid(doRender = true) {{
            const sortEl = document.getElementById("sortSelect");
            const sortMode = sortEl ? sortEl.value : "default";
            const spec = _SORT_SPECS[sortMode];
            if (sortMode === "default") {{
                currentData.sort((a, b) => _cmpNum(a.idx, b.idx, true));
            }} else if (spec) {{
                currentData.sort((a, b) => _cmpNum(a[spec[0]], b[spec[0]], spec[1]) || (a.idx - b.idx));
            }}
            if (doRender) {{
                if (currentView === 'grid') renderCards(currentData);
                else renderScatterPlot(currentData);
            }}
        }}

        function exportFilteredCSV() {{
            if (!currentData || currentData.length === 0) {{
                alert("No compounds to export.");
                return;
            }}
            const headers = ["ID", "SMILES", "Potency", "MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "Fsp3", "QED", "CNS_MPO", "Ro5_Violations", "PAINS_Status"];
            const rows = [headers.join(",")];
            currentData.forEach(d => {{
                const vals = [
                    `"${{d.id.replace(/"/g, '""')}}"`,
                    `"${{d.smiles.replace(/"/g, '""')}}"`,
                    `"${{(d.potency_display || '').replace(/"/g, '""')}}"`,
                    d.mw ?? '',
                    d.clogp ?? '',
                    d.tpsa ?? '',
                    d.hbd ?? '',
                    d.hba ?? '',
                    d.rotb ?? '',
                    d.fsp3 ?? '',
                    d.qed ?? '',
                    d.cns_mpo ?? '',
                    d.ro5_violations ?? '',
                    `"${{(d.pains || '').replace(/"/g, '""')}}"`
                ];
                rows.push(vals.join(","));
            }});
            const blob = new Blob([rows.join("\\n")], {{ type: "text/csv;charset=utf-8;" }});
            const url = URL.createObjectURL(blob);
            const a = document.createElement("a");
            a.href = url;
            a.download = "compounds_export.csv";
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
            showCopyToast("Exported filtered CSV file!");
        }}

        // Initial Render (Decompress lossless gzipped JSON payload via non-blocking pipeThrough)
        (async function initPayload() {{
            try {{
                const b64 = "{gzip_b64_payload}";
                const binStr = atob(b64);
                const bytes = new Uint8Array(binStr.length);
                for (let i = 0; i < binStr.length; i++) {{
                    bytes[i] = binStr.charCodeAt(i);
                }}
                const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
                const jsonText = await new Response(stream).text();
                rawData = JSON.parse(jsonText);
                filterGrid();
            }} catch (err) {{
                console.error("Failed to initialize compound payload:", err);
            }}
        }})();
    </script>
</body>
</html>"""

  def _write_text_atomic(target_path: str, text_content: str) -> None:
    target_dir = os.path.dirname(os.path.abspath(target_path))
    os.makedirs(target_dir, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(
        prefix=".tmp_", suffix=os.path.basename(target_path), dir=target_dir
    )
    try:
      with os.fdopen(fd, "w", encoding="utf-8") as tf:
        tf.write(text_content)
      os.replace(tmp_path, target_path)
    except BaseException:
      if os.path.exists(tmp_path):
        try:
          os.remove(tmp_path)
        except OSError:
          pass
      raise

  abs_out = os.path.abspath(output_html)
  _write_text_atomic(abs_out, html_content)

  now_iso = (
      datetime.datetime.now(datetime.timezone.utc)
      .isoformat()
      .replace("+00:00", "Z")
  )
  summary_txt = (
      f"{title}: Interactive 2D chemical structure grid and live SAR"
      f" scatter explorer ({len(compounds_data)} compounds)."
  )
  artifact_meta = {
      "summary": summary_txt,
      "updatedAt": now_iso,
      "userFacing": True,
      "requestFeedback": False,
  }
  meta_json_str = json.dumps(artifact_meta, indent=2)
  try:
    _write_text_atomic(abs_out + ".metadata.json", meta_json_str)
  except OSError:
    pass

  # If executed inside an explicit agent artifact workspace (via output_html or
  # cwd matching .../brain/<uuid>), mirror both the HTML and .metadata.json to
  # the root artifact directory so the agent UI renders it as a top-level
  # interactive HTML artifact rather than raw source code.
  brain_root = None
  for candidate_path in (abs_out, os.path.abspath(os.getcwd())):
    m_brain = re.search(r"(.*/brain/[0-9a-fA-F-]{36})(?:/|$)", candidate_path)
    if m_brain and os.path.isdir(m_brain.group(1)):
      brain_root = m_brain.group(1)
      break

  canonical_artifact_path = abs_out
  if brain_root:
    root_html = os.path.join(brain_root, os.path.basename(abs_out))
    canonical_artifact_path = root_html
    if os.path.abspath(root_html) != abs_out:
      try:
        _write_text_atomic(root_html, html_content)
        _write_text_atomic(root_html + ".metadata.json", meta_json_str)
      except OSError:
        canonical_artifact_path = abs_out

  base_name = os.path.basename(canonical_artifact_path)
  print(
      "Generated agent UI-compatible interactive grid & SAR explorer"
      f" ({len(compounds_data)} cards, {os.path.getsize(abs_out)/1024/1024:.2f}"
      f" MB) at: {canonical_artifact_path}\n[IMPORTANT ARTIFACT LINKING]: In"
      " Markdown reports and chat responses, ALWAYS link this interactive HTML"
      " viewer using its root artifact path:\n "
      f" [{base_name}](file://{canonical_artifact_path})\nNEVER link a source"
      " code workspace or external directory path for .html viewers"
      " (non-artifact paths open raw HTML source code instead of the"
      " interactive viewer)."
  )
  return canonical_artifact_path


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Generate a agent UI-compatible interactive 2D chemical structure"
          " grid & SAR scatter HTML."
      )
  )
  parser.add_argument(
      "--input",
      "--input_csv",
      "--input-csv",
      dest="input",
      type=str,
      required=True,
      help="Input CSV/SDF file path.",
  )
  parser.add_argument(
      "--smiles_col", type=str, default=None, help="Name of SMILES column."
  )
  parser.add_argument(
      "--id_col",
      "--name_col",
      type=str,
      default=None,
      dest="id_col",
      help="Name of ID/Compound Name column.",
  )
  parser.add_argument(
      "--align_core",
      type=str,
      default=None,
      help="Core SMARTS or SMILES to align all 2D depictions to.",
  )
  parser.add_argument(
      "--output",
      "--output_html",
      "--output-html",
      dest="output",
      type=str,
      default="molecular_grid.html",
      help="Output HTML file path.",
  )
  parser.add_argument(
      "--title",
      type=str,
      default="Interactive Molecular Grid & SAR Explorer",
      help="Page title.",
  )
  parser.add_argument(
      "--max_cards",
      type=int,
      default=250,
      help=(
          "Maximum number of cards to embed in interactive HTML (default: 250)."
      ),
  )
  parser.add_argument(
      "--activity_col",
      type=str,
      default=None,
      help=(
          "Name of activity/potency column (e.g. pIC50, pKi). Auto-detected if"
          " not specified."
      ),
  )
  parser.add_argument(
      "--only_active",
      "--filter_active",
      action="store_true",
      default=None,
      dest="filter_active",
      help=(
          "Filter dataset to only compounds with valid, non-null activity"
          " (default: True when activity column is present)."
      ),
  )
  parser.add_argument(
      "--all_compounds",
      action="store_true",
      default=False,
      help=(
          "Include all compounds including uncharacterized/untested records"
          " (overrides default active filtering)."
      ),
  )
  args = parser.parse_args()

  if args.input:
    if not os.path.exists(args.input):
      print(
          f"ERROR: input file not found: {args.input}",
          file=sys.stderr,
      )
      sys.exit(1)
    if args.input.endswith(".csv"):
      df = pd.read_csv(args.input)
    elif args.input.endswith(".sdf"):
      df = PandasTools.LoadSDF(args.input)
    else:
      print(
          f"ERROR: unsupported input format: {args.input}. Use .csv or .sdf",
          file=sys.stderr,
      )
      sys.exit(1)
  else:
    parser.error("--input is required.")

  filter_active = None
  if args.all_compounds:
    filter_active = False
  elif args.filter_active:
    filter_active = True

  try:
    create_molecular_grid(
        df,
        smiles_col=args.smiles_col,
        id_col=args.id_col,
        activity_col=args.activity_col,
        filter_active=filter_active,
        output_html=args.output,
        title=args.title,
        align_core=args.align_core,
        max_cards=args.max_cards,
    )
  except ValueError as exc:
    print(f"ERROR: {exc}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
  main()
