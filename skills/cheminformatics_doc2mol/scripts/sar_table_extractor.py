# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pandas",
# ]
# ///

"""Tabular SAR Bioactivity Extractor & Inequality Normalizer.

Parses experimental biological assay tables from patents and medicinal chemistry
papers, extracting quantitative concentrations with inequality operators (<, >,
<=, >=, =), standardizing units (pM, nM, uM, mM) to molarity (M), and computing
logarithmic potency scores (pIC50, pKi, pEC50, pActivity).
"""

from __future__ import annotations

import argparse
import math
import os
import re
from typing import Any, Dict, List, Optional

import pandas as pd

UNIT_FACTORS_TO_MOLAR = {
    "pm": 1e-12,
    "picomolar": 1e-12,
    "nm": 1e-9,
    "nanomolar": 1e-9,
    "um": 1e-6,
    "µm": 1e-6,
    "micromolar": 1e-6,
    "mm": 1e-3,
    "millimolar": 1e-3,
    "m": 1.0,
    "molar": 1.0,
    "mol/l": 1.0,
}


_LETTER_TIER_MAP = {
    "A": 4.0,
    "B": 3.0,
    "C": 2.0,
    "D": 1.0,
    "E": 0.5,
    "F": 0.25,
}


def parse_bioactivity_entry(raw_val: Any) -> Dict[str, Any]:
  """Parses a single raw bioactivity string into operator, numeric scalar, standardized molarity, and discrete potency bin/range."""
  empty_res = {
      "operator": "=",
      "numeric": None,
      "unit": "",
      "value_M": None,
      "pActivity": None,
      "activity_bin": None,
      "activity_range": None,
      "potency_tier": None,
  }
  if raw_val is None or pd.isna(raw_val):
    return dict(empty_res)

  s = str(raw_val).strip()
  if not s or s.lower() in ("nan", "none", "nd", "n/a", "-", "nt"):
    return dict(empty_res)

  # 1. Check letter-binned format with optional parenthetical concentration
  # range: e.g. 'A', 'A (<= 100 nM)', 'B (100-250 nM)'
  m_let = re.match(r"^([A-F])(?:\s*\(\s*(.+?)\s*\))?$", s, re.IGNORECASE)
  if m_let:
    letter = m_let.group(1).upper()
    rng_str = (m_let.group(2) or "").strip() or None
    tier = _LETTER_TIER_MAP.get(letter)
    res = dict(empty_res)
    res["activity_bin"] = letter
    res["activity_range"] = rng_str
    res["potency_tier"] = tier
    if rng_str:
      sub_parsed = parse_bioactivity_entry(rng_str)
      for k in ("operator", "numeric", "unit", "value_M", "pActivity"):
        val_k = sub_parsed.get(k)
        if val_k is not None and bool(str(val_k)):
          res[k] = val_k
    return res

  # 2. Check plus/star discrete potency bins: e.g. '++++', '+++', '++', '+'
  if re.match(r"^(\+{1,6}|\*{1,6})$", s):
    res = dict(empty_res)
    res["activity_bin"] = s
    res["potency_tier"] = float(min(len(s), 6))
    return res

  # 3. Check explicit concentration range: e.g. '100-250 nM', '250 nM - 1 uM'
  m_rng = re.match(
      r"^([0-9]+(?:\.[0-9]+)?)\s*([A-Za-zµμ/]+)?\s*(?:[-–—]|to)\s*([0-9]+(?:\.[0-9]+)?)\s*([A-Za-zµμ/]+)$",
      s,
      re.IGNORECASE,
  )
  if m_rng:
    v1 = float(m_rng.group(1))
    u1_raw = (m_rng.group(2) or "").strip()
    v2 = float(m_rng.group(3))
    u2_raw = (m_rng.group(4) or "").strip()
    u2 = u2_raw.replace("μ", "u").replace("µ", "u")
    u1 = (u1_raw or u2_raw).replace("μ", "u").replace("µ", "u")
    f1 = UNIT_FACTORS_TO_MOLAR.get(u1.lower(), 1e-9)
    f2 = UNIT_FACTORS_TO_MOLAR.get(u2.lower(), 1e-9)
    m1_val = v1 * f1
    m2_val = v2 * f2
    geom_m = (
        math.sqrt(m1_val * m2_val)
        if (m1_val > 0 and m2_val > 0)
        else 0.5 * (m1_val + m2_val)
    )
    p_act = round(-math.log10(geom_m), 3) if geom_m > 0 else None
    nm_hi = max(m1_val, m2_val) * 1e9
    tier = (
        4.0
        if nm_hi <= 100.0
        else (3.0 if nm_hi <= 250.0 else (2.0 if nm_hi <= 1000.0 else 1.0))
    )
    return {
        "operator": "range",
        "numeric": round(geom_m / f2, 4) if f2 > 0 else None,
        "unit": u2_raw,
        "value_M": geom_m,
        "pActivity": p_act,
        "activity_bin": None,
        "activity_range": s,
        "potency_tier": tier,
    }

  op = "="
  s_rem = s
  # Extract operator
  m_op = re.match(r"^([<>=≤≥~]+)\s*", s_rem)
  if m_op:
    raw_op = m_op.group(1)
    s_rem = s_rem[len(raw_op) :].strip()
    if "<=" in raw_op or "≤" in raw_op:
      op = "<="
    elif ">=" in raw_op or "≥" in raw_op:
      op = ">="
    elif "<" in raw_op:
      op = "<"
    elif ">" in raw_op:
      op = ">"
    elif "~" in raw_op or "=" in raw_op:
      op = "="

  # Extract numeric scalar and unit suffix
  m_num = re.match(
      r"^([0-9]+(?:\.[0-9]+)?(?:[eE][-+]?[0-9]+)?)\s*([A-Za-zµμ/]+)?", s_rem
  )
  if not m_num:
    res = dict(empty_res)
    res["operator"] = op
    return res

  num = float(m_num.group(1))
  unit = (m_num.group(2) or "").strip()
  unit_norm = unit.replace("μ", "u").replace("µ", "u").lower()

  factor = UNIT_FACTORS_TO_MOLAR.get(unit_norm, 1.0)
  val_m = num * factor

  p_act = None
  if val_m > 0:
    try:
      p_act = round(-math.log10(val_m), 3)
    except (ValueError, OverflowError):
      p_act = None

  tier = None
  if unit_norm in UNIT_FACTORS_TO_MOLAR and val_m > 0:
    nm_val = val_m * 1e9
    if op in (">", ">="):
      tier = 1.0 if nm_val >= 1000.0 else 2.0
    elif nm_val <= 100.0:
      tier = 4.0
    elif nm_val <= 250.0:
      tier = 3.0
    elif nm_val <= 1000.0:
      tier = 2.0
    else:
      tier = 1.0

  return {
      "operator": op,
      "numeric": num,
      "unit": unit,
      "value_M": val_m,
      "pActivity": p_act,
      "activity_bin": None,
      "activity_range": s if op != "=" else None,
      "potency_tier": tier,
  }


def extract_sar_table(
    input_csv: str,
    activity_cols: Optional[List[str]] = None,
    output_csv: Optional[str] = None,
) -> pd.DataFrame:
  """Extracts and normalizes bioactivity measurements across multiple assay columns."""
  if (
      isinstance(activity_cols, str)
      and output_csv is None
      and activity_cols.endswith(".csv")
  ):
    output_csv = activity_cols
    activity_cols = None
  df = pd.read_csv(input_csv)

  cols_to_process = [c for c in activity_cols or [] if c in df.columns]
  if not cols_to_process:
    # Auto-detect potential activity columns if none match explicitly
    cols_to_process = [
        c
        for c in df.columns
        if any(
            k in c.lower()
            for k in [
                "ic50",
                "ki",
                "ec50",
                "kd",
                "ed50",
                "potency",
                "activity",
                "patent_assay_",
                "patent_table_col_",
            ]
        )
        and not c.lower().endswith((
            "_range",
            "_units",
            "_unit",
            "_type",
            "_operator",
            "_nm",
            "_pactivity",
            "_numeric",
            "_value_m",
            "_potency_tier",
            "_bin",
        ))
        and c.lower()
        not in (
            "activity_type",
            "activity_units",
            "activity_range",
            "pchembl_value",
            "pactivity",
            "pic50",
            "pki",
            "pec50",
            "pkd",
            "potency_tier",
        )
    ]

  for col in cols_to_process:
    col_lower = col.lower()
    parsed = [parse_bioactivity_entry(v) for v in df[col]]

    df[f"{col}_operator"] = [p["operator"] for p in parsed]
    df[f"{col}_numeric"] = [p["numeric"] for p in parsed]
    df[f"{col}_unit"] = [p["unit"] for p in parsed]
    df[f"{col}_value_M"] = [p["value_M"] for p in parsed]
    df[f"{col}_pActivity"] = [p["pActivity"] for p in parsed]
    df[f"{col}_potency_tier"] = [p["potency_tier"] for p in parsed]
    if any(p["activity_bin"] is not None for p in parsed):
      df[f"{col}_bin"] = [p["activity_bin"] for p in parsed]
      if "activity_bin" not in df.columns:
        df["activity_bin"] = df[f"{col}_bin"]
    if any(p["activity_range"] is not None for p in parsed):
      df[f"{col}_range"] = [p["activity_range"] for p in parsed]
      if "activity_range" not in df.columns:
        df["activity_range"] = df[f"{col}_range"]

    # Add standard convenience aliases
    if "ic50" in col_lower:
      df["pIC50"] = df[f"{col}_pActivity"]
      df["IC50_operator"] = df[f"{col}_operator"]
    elif "ki" in col_lower:
      df["pKi"] = df[f"{col}_pActivity"]
      df["Ki_operator"] = df[f"{col}_operator"]
    elif "ec50" in col_lower:
      df["pEC50"] = df[f"{col}_pActivity"]
      df["EC50_operator"] = df[f"{col}_operator"]
    else:
      df[f"p{col}"] = df[f"{col}_pActivity"]

  # Global pActivity / potency_tier convenience columns across processed columns
  if cols_to_process:
    primary = cols_to_process[0]
    df["pActivity"] = df[f"{primary}_pActivity"]
    df["operator"] = df[f"{primary}_operator"]
    df["potency_tier"] = df[f"{primary}_potency_tier"]
    for extra_col in cols_to_process[1:]:
      df["pActivity"] = df["pActivity"].combine_first(
          df[f"{extra_col}_pActivity"]
      )
      df["potency_tier"] = df["potency_tier"].combine_first(
          df[f"{extra_col}_potency_tier"]
      )

  if output_csv:
    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    df.to_csv(output_csv, index=False)
    print(f"Exported {len(df)} normalized SAR records to {output_csv}")

  return df


def main():
  parser = argparse.ArgumentParser(
      description="SAR Bioactivity Extractor & Inequality Normalizer."
  )
  parser.add_argument(
      "--input-csv",
      "--input_csv",
      "--input",
      type=str,
      required=True,
      help="Input CSV containing SAR table.",
  )
  parser.add_argument(
      "--output-csv",
      "--output_csv",
      "--output",
      type=str,
      required=True,
      help="Path for output normalized CSV.",
  )
  parser.add_argument(
      "--activity-cols",
      "--activity_cols",
      nargs="+",
      default=["IC50_raw", "Ki_raw"],
      help="Activity column names to parse.",
  )
  args = parser.parse_args()

  extract_sar_table(
      input_csv=args.input_csv,
      activity_cols=args.activity_cols,
      output_csv=args.output_csv,
  )


if __name__ == "__main__":
  main()
