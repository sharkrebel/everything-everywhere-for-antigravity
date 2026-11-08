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

"""Performs Bemis-Murcko scaffold splitting for machine learning dataset partitioning."""

import argparse
import calendar
from collections import defaultdict
import datetime
import math
import re
import sys
from typing import Any, Dict, Final, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def group_by_scaffold(smiles_list: List[str]) -> Dict[str, List[int]]:
  """Groups dataset row indices by their Bemis-Murcko scaffold.

  Args:
      smiles_list: List of SMILES strings.

  Returns:
      Mapping of scaffold SMILES to the row indices sharing that scaffold.
      Unparseable entries are collected under the "INVALID_SMILES" key.
  """
  scaffold_groups = defaultdict(list)
  for idx, smi in enumerate(smiles_list):
    mol = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None
    if mol is None:
      scaffold = "INVALID_SMILES"
    else:
      scaffold = MurckoScaffold.MurckoScaffoldSmiles(
          mol=mol, includeChirality=False
      )
    scaffold_groups[scaffold].append(idx)
  return dict(scaffold_groups)


def scaffold_split(
    smiles_list: List[str],
    train_ratio: float = 0.8,
    val_ratio: float = 0.1,
    test_ratio: float = 0.1,
    strategy: str = "balanced",
) -> Tuple[List[int], List[int], List[int]]:
  """Partitions molecules into Train, Val, and Test sets by Bemis-Murcko scaffold.

  A scaffold group is never split across partitions, so there is zero scaffold
  overlap between train and test regardless of strategy.

  Two assignment strategies are available:

  - ``balanced`` (default, deterministic): scaffold groups are visited largest
    first and each is placed in whichever split is currently furthest below its
    quota, measured as absolute remaining capacity. This keeps every split
    representative of the dataset's scaffold size distribution.
  - ``greedy`` (legacy): fills train to its cutoff, then val, then dumps the
    remainder into test. Because groups are visited largest first, every
    populous scaffold lands in train and the test set degenerates into the
    tail of singleton scaffolds. Retained only for reproducing older splits.

  Both strategies are fully deterministic: ties are broken by scaffold SMILES,
  so the same input always yields the same partition.

  Args:
      smiles_list: List of SMILES strings.
      train_ratio: Fraction for training set.
      val_ratio: Fraction for validation set.
      test_ratio: Fraction for test set.
      strategy: Either "balanced" or "greedy".

  Returns:
      Tuple of (train_indices, val_indices, test_indices).

  Raises:
      ValueError: If strategy is unrecognized, any split ratio is negative, or
        no split ratio is positive.
  """
  if strategy not in ("balanced", "greedy"):
    raise ValueError(f"Unknown split strategy: {strategy!r}")

  ratio_items = [
      ("train", train_ratio),
      ("val", val_ratio),
      ("test", test_ratio),
  ]
  if any(r < 0 for _, r in ratio_items):
    raise ValueError(
        "Split ratios must be non-negative; got"
        f" train_ratio={train_ratio}, val_ratio={val_ratio},"
        f" test_ratio={test_ratio}."
    )
  active_names = [name for name, r in ratio_items if r > 0]
  if not active_names:
    raise ValueError(
        "At least one split ratio must be positive; got"
        f" train_ratio={train_ratio}, val_ratio={val_ratio},"
        f" test_ratio={test_ratio}."
    )
  if not smiles_list:
    return [], [], []

  scaffold_groups = group_by_scaffold(smiles_list)

  # Deterministic ordering: largest groups first, ties broken by scaffold SMILES.
  ordered = sorted(scaffold_groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))

  total = len(smiles_list)
  train_idx: List[int] = []
  val_idx: List[int] = []
  test_idx: List[int] = []
  ratio_sum = sum(r for _, r in ratio_items if r > 0)

  if strategy == "greedy":
    norm = {name: r / ratio_sum if r > 0 else 0.0 for name, r in ratio_items}
    train_cutoff = norm["train"] * total
    val_cutoff = (norm["train"] + norm["val"]) * total
    last_active = active_names[-1]
    buckets = {"train": train_idx, "val": val_idx, "test": test_idx}
    for _, group in ordered:
      if norm["train"] > 0 and (
          last_active == "train" or len(train_idx) + len(group) <= train_cutoff
      ):
        train_idx.extend(group)
      elif norm["val"] > 0 and (
          last_active == "val"
          or len(train_idx) + len(val_idx) + len(group) <= val_cutoff
      ):
        val_idx.extend(group)
      else:
        buckets[last_active].extend(group)
    return train_idx, val_idx, test_idx

  # Balanced: place each group in the split where it fits best proportionally.
  # The score is the fill ratio *after* hypothetically adding the group, so a
  # large scaffold family is not dumped into a small quota. Scoring the ratio
  # before assignment (or using absolute remaining capacity) lets the largest
  # quota monopolise every group and starves the test split.
  assignment: Dict[str, List[List[int]]] = {"train": [], "val": [], "test": []}
  quotas = {
      name: max((r / ratio_sum) * total, 1e-9)
      for name, r in ratio_items
      if r > 0
  }

  def filled(split: str) -> int:
    return sum(len(g) for g in assignment[split])

  for _, group in ordered:
    # Enumerate index is the final tie-break, keeping the result deterministic.
    target = min(
        enumerate(active_names),
        key=lambda inm: (
            (filled(inm[1]) + len(group)) / quotas[inm[1]],
            inm[0],
        ),
    )[1]
    assignment[target].append(group)

  # Repair: an active split can end up empty when every scaffold family is
  # larger than its quota. Donate the group that best fits, taken from the
  # active split that can most afford to lose one.
  if len(ordered) >= 2 and len(active_names) >= 2:
    for split in ("test", "val", "train"):
      if split not in quotas or assignment[split]:
        continue
      donor = max(
          (n for n in active_names if len(assignment[n]) >= 2),
          key=lambda n: (len(assignment[n]), -active_names.index(n)),
          default=None,
      )
      if donor is None:
        continue
      best = min(
          range(len(assignment[donor])),
          key=lambda i: (
              abs(len(assignment[donor][i]) - quotas[split]),
              -len(assignment[donor][i]),
          ),
      )
      assignment[split].append(assignment[donor].pop(best))
  elif len(active_names) >= 2:
    print(
        f"Warning: only {len(ordered)} distinct scaffold(s) present; a "
        f"{train_ratio:.0%}/{val_ratio:.0%}/{test_ratio:.0%} split with zero "
        "scaffold overlap is not achievable.",
        file=sys.stderr,
    )

  for group in assignment["train"]:
    train_idx.extend(group)
  for group in assignment["val"]:
    val_idx.extend(group)
  for group in assignment["test"]:
    test_idx.extend(group)

  if (val_ratio > 0 and not val_idx) or (test_ratio > 0 and not test_idx):
    print(
        "Warning: scaffold split produced an empty evaluation partition. The "
        "dataset has too few distinct scaffolds for the requested ratios.",
        file=sys.stderr,
    )

  return train_idx, val_idx, test_idx


def summarize_split(
    smiles_list: List[str],
    train_idx: List[int],
    val_idx: List[int],
    test_idx: List[int],
) -> Dict[str, Dict[str, Any]]:
  """Reports per-split compound counts and scaffold size distribution.

  Args:
      smiles_list: The SMILES list that was split.
      train_idx: Training row indices.
      val_idx: Validation row indices.
      test_idx: Test row indices.

  Returns:
      Mapping of split name to its compound count, distinct scaffold count,
      and the largest and median scaffold group size within that split.
  """
  scaffold_of: Dict[int, str] = {}
  for scaffold, idxs in group_by_scaffold(smiles_list).items():
    for i in idxs:
      scaffold_of[i] = scaffold

  out: Dict[str, Dict[str, Any]] = {}
  for name, idxs in (
      ("train", train_idx),
      ("val", val_idx),
      ("test", test_idx),
  ):
    sizes: Dict[str, int] = defaultdict(int)
    for i in idxs:
      sizes[scaffold_of[i]] += 1
    counts = sorted(sizes.values(), reverse=True)
    median = counts[len(counts) // 2] if counts else 0
    out[name] = {
        "compounds": len(idxs),
        "scaffolds": len(counts),
        "largest_scaffold": counts[0] if counts else 0,
        "median_scaffold": median,
    }
  return out


_INVALID_TIME_STRINGS: Final[frozenset[str]] = frozenset(
    {"", "nan", "none", "null", "na", "n/a", "true", "false"}
)

_PUBLIC_CALENDAR_COLUMNS: Final[Tuple[str, ...]] = (
    "publication_date",
    "document_year",
    "publication_year",
    "chembl_year",
    "year",
    "Year",
    "date",
    "Date",
)

_FALLBACK_CALENDAR_COLUMNS: Final[Tuple[str, ...]] = (
    "filing_date",
    "priority_date",
)

_CALENDAR_TIME_COLUMNS: Final[Tuple[str, ...]] = (
    *_PUBLIC_CALENDAR_COLUMNS,
    *_FALLBACK_CALENDAR_COLUMNS,
)

_PUBLIC_SUBYEAR_COLUMNS: Final[Tuple[str, ...]] = (
    "publication_date",
    "date",
    "Date",
)

_DOCUMENT_ID_COLUMNS: Final[Tuple[str, ...]] = (
    "document_chembl_id",
    "document_id",
    "patent_id",
    "pubmed_id",
    "pmid",
    "doi",
)

_COMPOUND_ID_COLUMNS: Final[Tuple[str, ...]] = (
    "chembl_id",
    "molecule_chembl_id",
    "inchikey",
    "standard_inchi_key",
)

_ISO_DATE_RE: Final[re.Pattern[str]] = re.compile(
    r"^(\d{4})[-/](\d{1,2})(?:[-/](\d{1,2}))?(?:[T\s].*)?$"
)
_COMPACT_DATE_RE: Final[re.Pattern[str]] = re.compile(
    r"^(19\d{2}|20\d{2})(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])$"
)
_WIPO_YEAR_RE: Final[re.Pattern[str]] = re.compile(
    r"\bWO[-\s/]?((?:19|20)\d{2}|\d{2})[-\s/]?\d{4,7}(?:[-\s/]?[A-Z]\d?)?\b",
    re.IGNORECASE,
)
_US_PUB_YEAR_RE: Final[re.Pattern[str]] = re.compile(
    r"\bUS[-\s/]?(20\d{2})[-\s/]?\d{6,8}(?:[-\s/]?A\d)?\b(?![-\s/]?B\d)",
    re.IGNORECASE,
)


def is_calendar_col_name(col: Optional[str]) -> bool:
  """Returns True when `col` denotes a calendar date or publication year."""
  if col is None:
    return True
  lowered = str(col).strip().lower()
  if lowered in {c.lower() for c in _CALENDAR_TIME_COLUMNS}:
    return True
  return any(tok in lowered for tok in ("year", "date"))


_is_calendar_col_name = is_calendar_col_name


def _is_coarse_calendar_year(val: float) -> bool:
  """Returns True if `val` is an exact integer year in [1800, 2100]."""
  return (
      math.isfinite(val)
      and 1800.0 <= val <= 2100.0
      and abs(val - round(val)) < 1e-7
  )


def _is_subyear_calendar_float(val: float) -> bool:
  """Returns True if `val` encodes a sub-year calendar date in [1800, 2101)."""
  return (
      math.isfinite(val)
      and 1800.0 <= val < 2101.0
      and abs(val - round(val)) >= 1e-7
  )


def _date_to_fractional_year(year: int, month: int, day: int) -> float:
  """Converts a calendar date into a monotonic fractional year.

  Uses `(yday - 0.5) / days_in_year` so every calendar date has a non-zero
  fractional part in `[0.00136, 0.99864]`, keeping exact integers (`float(Y)`)
  bijectively reserved for coarse year-only values.

  Args:
      year: Four-digit calendar year.
      month: Calendar month (1..12).
      day: Day of month (1..31).

  Returns:
      Monotonic fractional year in `(year, year + 1)`.
  """
  dt = datetime.date(year, month, day)
  days_in_year = 366 if calendar.isleap(year) else 365
  yday = dt.timetuple().tm_yday
  return float(year) + (float(yday) - 0.5) / float(days_in_year)


def _parse_single_temporal_token(val: object) -> Tuple[float, bool]:
  """Parses a single scalar timestamp/year token, rejecting booleans.

  Args:
      val: Scalar timestamp, ISO date string, integer year, or object.

  Returns:
      Tuple of (numeric_value_or_nan, is_subyear_date).
  """
  if val is None or isinstance(val, (bool, np.bool_)):
    return np.nan, False
  try:
    if pd.isna(val):
      return np.nan, False
  except (TypeError, ValueError):
    pass

  if isinstance(val, datetime.datetime):
    return _date_to_fractional_year(val.year, val.month, val.day), True
  if isinstance(val, datetime.date):
    return _date_to_fractional_year(val.year, val.month, val.day), True

  if isinstance(val, (int, float, np.integer, np.floating)):
    fval = float(val)
    if not math.isfinite(fval):
      return np.nan, False
    if fval.is_integer() and 19000101 <= int(fval) <= 20991231:
      m_compact = _COMPACT_DATE_RE.match(str(int(fval)))
      if m_compact:
        try:
          return (
              _date_to_fractional_year(
                  int(m_compact.group(1)),
                  int(m_compact.group(2)),
                  int(m_compact.group(3)),
              ),
              True,
          )
        except ValueError:
          return np.nan, False
    return fval, _is_subyear_calendar_float(fval)

  text = str(val).strip()
  if text.lower() in _INVALID_TIME_STRINGS:
    return np.nan, False

  m_iso = _ISO_DATE_RE.match(text)
  if m_iso:
    year = int(m_iso.group(1))
    month = int(m_iso.group(2))
    day = int(m_iso.group(3)) if m_iso.group(3) is not None else 1
    try:
      return _date_to_fractional_year(year, month, day), True
    except ValueError:
      return np.nan, False

  m_compact = _COMPACT_DATE_RE.match(text)
  if m_compact:
    try:
      return (
          _date_to_fractional_year(
              int(m_compact.group(1)),
              int(m_compact.group(2)),
              int(m_compact.group(3)),
          ),
          True,
      )
    except ValueError:
      return np.nan, False

  try:
    fval = float(text)
  except ValueError:
    return np.nan, False
  if not math.isfinite(fval):
    return np.nan, False
  return fval, _is_subyear_calendar_float(fval)


def parse_temporal_series(
    time_series: Any,
    calendar_bounds: bool = False,
) -> Tuple[pd.Series, bool]:
  """Parses a time series into finite floats, rejecting boolean placeholders.

  Boolean values (`True`, `False`, `np.bool_`, `"True"`, `"False"`) are never
  treated as numeric `0`/`1`, preventing `False < True` from masquerading as a
  valid chronological split.

  Args:
      time_series: Sequence or Series of years, ISO dates, or numeric indices.
      calendar_bounds: If True, reject numeric values outside [1800, 2101) (e.g.
        0/1/-1/9999 placeholders in calendar year/date columns).

  Returns:
      Tuple of (float64 Series preserving input index, has_subyear_dates).
  """
  if isinstance(time_series, pd.Series):
    orig_index = time_series.index
    s = time_series
  else:
    s = pd.Series(time_series)
    orig_index = s.index

  if pd.api.types.is_bool_dtype(s.dtype):
    return pd.Series(np.nan, index=orig_index, dtype=float), False

  values: List[float] = []
  has_subyear = False
  for item in s.tolist():
    parsed, is_subyear = _parse_single_temporal_token(item)
    if calendar_bounds and not math.isnan(parsed):
      if parsed < 1800.0 or parsed >= 2101.0:
        parsed = np.nan
        is_subyear = False
    if is_subyear:
      has_subyear = True
    values.append(parsed)
  return pd.Series(values, index=orig_index, dtype=float), has_subyear


def parse_temporal_cutoff(
    cutoff: Union[int, float, str],
    has_subyear_dates: bool = False,
) -> float:
  """Parses an explicit `--time-cutoff` value, rejecting booleans and NaNs.

  When `has_subyear_dates` is True and `cutoff` is an integer year (e.g. `2014`
  or `"2014"`), the cutoff is placed at December 31 of that calendar year so
  `time <= 2014` includes all dates within 2014.

  Args:
      cutoff: Explicit numeric year or ISO date cutoff.
      has_subyear_dates: Whether the target series contains sub-year dates.

  Returns:
      Parsed float64 cutoff threshold.
  """
  if cutoff is None or isinstance(cutoff, (bool, np.bool_)):
    raise ValueError(
        f"Invalid temporal cutoff {cutoff!r}: boolean or None values are not"
        " valid time cutoffs."
    )
  if isinstance(cutoff, str) and cutoff.strip().lower() in ("true", "false"):
    raise ValueError(
        f"Invalid temporal cutoff {cutoff!r}: boolean strings are not valid"
        " time cutoffs."
    )

  is_plain_year = False
  if isinstance(cutoff, (int, np.integer)):
    is_plain_year = True
  elif isinstance(cutoff, (float, np.floating)) and float(cutoff).is_integer():
    is_plain_year = True
  elif isinstance(cutoff, str) and re.fullmatch(r"\s*\d{4}\s*", cutoff):
    is_plain_year = True

  parsed, is_subyear = _parse_single_temporal_token(cutoff)
  if math.isnan(parsed):
    raise ValueError(f"Could not parse temporal cutoff {cutoff!r}.")

  if has_subyear_dates and is_plain_year and not is_subyear:
    yr = int(round(parsed))
    if 1800 <= yr <= 2100:
      return _date_to_fractional_year(yr, 12, 31)
    return float(yr) + 1.0 - 1e-6
  return parsed


def format_temporal_value(val: float, has_subyear_dates: bool = False) -> str:
  """Formats a parsed temporal float for human-readable logging and JSON."""
  if not math.isfinite(val):
    return "NaN"
  if _is_subyear_calendar_float(val) or (
      has_subyear_dates
      and 1800.0 <= val < 2101.0
      and abs(val - round(val)) >= 1e-7
  ):
    year = int(math.floor(val + 1e-9))
    days_in_year = 366 if calendar.isleap(year) else 365
    yday = min(
        days_in_year,
        max(1, int(round((val - year) * days_in_year + 0.5))),
    )
    dt = datetime.date(year, 1, 1) + datetime.timedelta(days=yday - 1)
    return dt.isoformat()
  if abs(val - round(val)) < 1e-6:
    return str(int(round(val)))
  return f"{val:.4f}"


def extract_identifier_year(value: object) -> Optional[float]:
  """Extracts the earliest unambiguous publication year from patent IDs.

  Only extracts years from identifiers that canonically encode a publication
  year: WIPO PCT (`WO-YYYY-NNNNNN` / `WO99/12345`) and US pre-grant publications
  (`USYYYYNNNNNNNA1`). Never guesses from granted US patent numbers.

  Args:
      value: Raw patent or document identifier token(s).

  Returns:
      Earliest extracted four-digit year as a float, or None if absent.
  """
  if value is None or isinstance(value, (bool, np.bool_)):
    return None
  try:
    if pd.isna(value):
      return None
  except (TypeError, ValueError):
    pass
  text = str(value).strip()
  if not text or text.lower() in _INVALID_TIME_STRINGS:
    return None

  years: List[int] = []
  for match in _WIPO_YEAR_RE.finditer(text):
    raw_yr = match.group(1)
    if len(raw_yr) == 4:
      yr = int(raw_yr)
    else:
      yr2 = int(raw_yr)
      yr = (1900 + yr2) if yr2 >= 78 else (2000 + yr2)
    if 1978 <= yr <= 2100:
      years.append(yr)

  for match in _US_PUB_YEAR_RE.finditer(text):
    yr = int(match.group(1))
    if 2001 <= yr <= 2100:
      years.append(yr)

  return float(min(years)) if years else None


def backfill_temporal_series(
    df: pd.DataFrame,
    time_col: Optional[str] = None,
) -> Tuple[pd.Series, bool, Optional[str]]:
  """Sanitizes and backfills true chronological dates across a DataFrame.

  1. Rejects boolean values (`False`/`True`) and placeholder strings -> `NaN`.
  2. Backfills missing entries in `time_col` from companion date columns
     (`publication_date`, `document_year`, `priority_date`, `filing_date`, etc.)
     and clamps to any strictly earlier calendar disclosure year.
  3. Propagates known document dates onto undated rows of the same document
     (`document_chembl_id`, `document_id`, `patent_id`, `pmid`, `doi`), then
     extracts unambiguous WIPO/US-A1 patent publication years when needed.
  4. Propagates each compound's earliest disclosure date across rows sharing
     `chembl_id` or `inchikey` without letting coarse year-only values erase
     sub-year precision within the same calendar year.

  Args:
      df: Input DataFrame.
      time_col: Preferred time column name (if None or missing, auto-discovers
        standard calendar columns preferring ISO `publication_date`).

  Returns:
      Tuple of (backfilled float64 Series preserving `df.index`,
      has_subyear_dates, resolved_col_name).
  """
  orig_index = df.index
  n_rows = len(df)
  if n_rows == 0:
    return pd.Series([], index=orig_index, dtype=float), False, time_col

  work_df = df.reset_index(drop=True)
  resolved_col: Optional[str] = None
  explicit_time_col = bool(time_col and time_col in work_df.columns)
  if explicit_time_col:
    resolved_col = time_col
  else:
    for cand in _CALENDAR_TIME_COLUMNS:
      if cand in work_df.columns:
        resolved_col = cand
        break

  is_calendar = _is_calendar_col_name(resolved_col)
  if resolved_col is not None:
    out, has_subyear = parse_temporal_series(
        work_df[resolved_col], calendar_bounds=is_calendar
    )
  else:
    out = pd.Series(np.nan, index=range(n_rows), dtype=float)
    has_subyear = False

  # If a custom non-calendar numeric column is already complete, return it as-is
  if not is_calendar and out.notna().sum() >= 2 and not out.isna().any():
    out.index = orig_index
    return out, has_subyear, resolved_col

  # Step 1: Backfill missing entries, clamp across public disclosure columns,
  # and refine coarse Y.0 entries from matching-year public sub-year dates.
  if is_calendar:
    for sec_col in _PUBLIC_CALENDAR_COLUMNS:
      if sec_col == resolved_col or sec_col not in work_df.columns:
        continue
      sec_vals, _ = parse_temporal_series(
          work_df[sec_col], calendar_bounds=True
      )
      if not sec_vals.notna().any():
        continue
      allow_subyear_refine = sec_col in _PUBLIC_SUBYEAR_COLUMNS
      updated_any = False
      for i in range(n_rows):
        s_val = float(sec_vals.iloc[i])
        if math.isnan(s_val):
          continue
        c_val = float(out.iloc[i])
        if math.isnan(c_val):
          out.iloc[i] = s_val
          updated_any = True
        elif int(math.floor(s_val + 1e-9)) < int(math.floor(c_val + 1e-9)):
          out.iloc[i] = s_val
          updated_any = True
        elif (
            allow_subyear_refine
            and _is_coarse_calendar_year(c_val)
            and _is_subyear_calendar_float(s_val)
            and int(round(c_val)) == int(math.floor(s_val + 1e-9))
        ):
          out.iloc[i] = s_val
          updated_any = True
      if updated_any and resolved_col is None:
        resolved_col = sec_col

  def _prefer_earlier_date(cand: float, existing: float) -> float:
    cand_yr = int(math.floor(cand + 1e-9))
    exist_yr = int(math.floor(existing + 1e-9))
    if cand_yr < exist_yr:
      return cand
    if exist_yr < cand_yr:
      return existing
    # Same calendar year: prefer fine-grained sub-year date over coarse integer
    cand_sub = _is_subyear_calendar_float(cand)
    exist_sub = _is_subyear_calendar_float(existing)
    if cand_sub and exist_sub:
      return min(cand, existing)
    if cand_sub:
      return cand
    return existing

  def _propagate_by_columns(
      series: pd.Series,
      id_cols: Tuple[str, ...],
      overwrite_with_earlier: bool = False,
  ) -> pd.Series:
    updated = series.copy()
    for id_col in id_cols:
      if id_col not in work_df.columns:
        continue
      if not overwrite_with_earlier and not updated.isna().any():
        continue
      token_min_date: Dict[str, float] = {}
      raw_ids = work_df[id_col].tolist()
      for idx, raw_id in enumerate(raw_ids):
        t_val = updated.iloc[idx]
        if (
            pd.isna(t_val)
            or raw_id is None
            or isinstance(raw_id, (bool, np.bool_))
        ):
          continue
        try:
          if pd.isna(raw_id):
            continue
        except (TypeError, ValueError):
          pass
        for tok in re.split(r"[,;|]+", str(raw_id)):
          tok = tok.strip()
          if not tok or tok.lower() in _INVALID_TIME_STRINGS:
            continue
          if tok not in token_min_date:
            token_min_date[tok] = float(t_val)
          else:
            token_min_date[tok] = _prefer_earlier_date(
                float(t_val), token_min_date[tok]
            )
      if not token_min_date:
        continue
      for idx, raw_id in enumerate(raw_ids):
        cur_val = updated.iloc[idx]
        if (
            (not overwrite_with_earlier and pd.notna(cur_val))
            or raw_id is None
            or isinstance(raw_id, (bool, np.bool_))
        ):
          continue
        try:
          if pd.isna(raw_id):
            continue
        except (TypeError, ValueError):
          pass
        matched = [
            token_min_date[tok.strip()]
            for tok in re.split(r"[,;|]+", str(raw_id))
            if tok.strip() in token_min_date
        ]
        if matched:
          earliest = matched[0]
          for m_val in matched[1:]:
            earliest = _prefer_earlier_date(m_val, earliest)
          if pd.isna(cur_val):
            updated.iloc[idx] = earliest
          else:
            updated.iloc[idx] = _prefer_earlier_date(earliest, float(cur_val))
    return updated

  # Step 2: Propagate known document dates onto undated rows of same document
  if is_calendar and out.isna().any() and out.notna().any():
    out = _propagate_by_columns(out, _DOCUMENT_ID_COLUMNS)

  # Step 3: Extract unambiguous WIPO/US-A1 years for remaining undated rows
  if is_calendar and out.isna().any():
    id_source_cols = [
        c
        for c in ("patent_id", "document_id", "source_id", "compound_id")
        if c in work_df.columns
    ]
    extracted_any = False
    for idx in out.index[out.isna()]:
      inferred_years = [
          yr
          for c in id_source_cols
          if (yr := extract_identifier_year(work_df[c].iloc[idx])) is not None
      ]
      if inferred_years:
        out.iloc[idx] = min(inferred_years)
        extracted_any = True
        if resolved_col is None:
          resolved_col = "patent_id"
    if extracted_any and out.isna().any():
      out = _propagate_by_columns(out, _DOCUMENT_ID_COLUMNS)

  # Step 3b: Fallback to non-public filing_date / priority_date ONLY when public
  # disclosure dates and unambiguous patent publication years are absent.
  if is_calendar and out.isna().any():
    filled_fallback = False
    for fb_col in _FALLBACK_CALENDAR_COLUMNS:
      if fb_col == resolved_col or fb_col not in work_df.columns:
        continue
      fb_vals, _ = parse_temporal_series(work_df[fb_col], calendar_bounds=True)
      if not fb_vals.notna().any():
        continue
      for i in range(n_rows):
        if math.isnan(float(out.iloc[i])) and not math.isnan(
            float(fb_vals.iloc[i])
        ):
          out.iloc[i] = float(fb_vals.iloc[i])
          filled_fallback = True
          if resolved_col is None:
            resolved_col = fb_col
    if filled_fallback and out.isna().any():
      out = _propagate_by_columns(out, _DOCUMENT_ID_COLUMNS)

  # Step 4: Assign each compound its earliest disclosure date across rows
  if is_calendar and out.notna().any():
    out = _propagate_by_columns(
        out, _COMPOUND_ID_COLUMNS, overwrite_with_earlier=True
    )

  has_subyear = any(
      _is_subyear_calendar_float(float(v)) for v in out.dropna().tolist()
  )
  out.index = orig_index
  if resolved_col is not None:
    out.name = resolved_col
  return out, has_subyear, resolved_col


def resolve_temporal_cutoff(
    numeric_s: pd.Series,
    train_ratio: float = 0.8,
    val_ratio: float = 0.0,
    test_ratio: float = 0.2,
) -> Tuple[float, Optional[float]]:
  """Selects strict chronological cutoff(s) with zero timestamp overlap.

  Every observation with an identical timestamp stays on the same side of the
  cutoff (`train <= c_train`, `c_train < val <= c_val`, `test > c_end`),
  guaranteeing `max(t_train) < min(t_test)`. When coarse year-only values and
  sub-year ISO dates coexist within the same calendar year `Y`, candidate
  cutoffs inside year `Y` are snapped to Dec 31 of year `Y` so observations of
  unknown month/day in year `Y` are never split against `Y-MM-DD` observations.

  Args:
      numeric_s: Parsed float64 Series of timestamps.
      train_ratio: Target training fraction.
      val_ratio: Target validation fraction (0.0 for 2-way split).
      test_ratio: Target test fraction.

  Returns:
      Tuple of (train_cutoff, val_cutoff_or_None).
  """
  ratio_items = [
      ("train", train_ratio),
      ("val", val_ratio),
      ("test", test_ratio),
  ]
  if any(r < 0 for _, r in ratio_items):
    raise ValueError(
        "Split ratios must be non-negative; got"
        f" train_ratio={train_ratio}, val_ratio={val_ratio},"
        f" test_ratio={test_ratio}."
    )
  if train_ratio <= 0 or (val_ratio <= 0 and test_ratio <= 0):
    raise ValueError(
        "Temporal split requires positive train_ratio and at least one positive"
        f" evaluation ratio; got train_ratio={train_ratio},"
        f" val_ratio={val_ratio}, test_ratio={test_ratio}."
    )

  valid_vals = np.sort(numeric_s.dropna().to_numpy(dtype=float))
  unique_vals = np.unique(valid_vals)
  n_active = sum(1 for _, r in ratio_items if r > 0)
  if len(unique_vals) < n_active:
    raise ValueError(
        f"Temporal split requires at least {n_active} distinct valid"
        f" non-boolean timestamps to form {n_active} non-overlapping"
        f" chronological partitions, but found {len(unique_vals)}."
    )

  coarse_years = {
      int(round(v)) for v in unique_vals if _is_coarse_calendar_year(float(v))
  }
  subyear_years = {
      int(math.floor(v + 1e-9))
      for v in unique_vals
      if _is_subyear_calendar_float(float(v))
  }
  ambiguous_years = coarse_years & subyear_years

  raw_candidates: List[float] = []
  max_val = float(valid_vals[-1])
  for u in unique_vals[:-1]:
    u_f = float(u)
    yr_u = int(math.floor(u_f + 1e-9))
    if yr_u in ambiguous_years:
      u_f = _date_to_fractional_year(yr_u, 12, 31)
    if u_f < max_val - 1e-9:
      raw_candidates.append(u_f)

  candidates = sorted(set(raw_candidates))
  if len(candidates) < (n_active - 1):
    raise ValueError(
        "Cannot form non-overlapping chronological partitions: coarse year-only"
        f" and sub-year dates coexist in year(s) {sorted(ambiguous_years)} with"
        " no observations in distinct calendar years."
    )

  ratio_sum = train_ratio + val_ratio + test_ratio
  norm_train = train_ratio / ratio_sum
  norm_val = val_ratio / ratio_sum
  n_valid = float(len(valid_vals))

  if val_ratio <= 0:
    best_c = min(
        candidates,
        key=lambda c: (
            abs(
                float(np.searchsorted(valid_vals, c, side="right")) / n_valid
                - norm_train
            ),
            -float(c),
        ),
    )
    return float(best_c), None

  # 3-way chronological split: choose c_tr < c_va from valid candidates
  target_train = norm_train
  target_train_val = norm_train + norm_val
  best_pair: Optional[Tuple[float, float]] = None
  best_cost: Optional[Tuple[float, float, float]] = None
  m = len(candidates)
  for i in range(0, m - 1):
    c_tr = float(candidates[i])
    cnt_tr = int(np.searchsorted(valid_vals, c_tr, side="right"))
    if cnt_tr <= 0:
      continue
    frac_tr = float(cnt_tr) / n_valid
    for j in range(i + 1, m):
      c_va = float(candidates[j])
      cnt_tr_va = int(np.searchsorted(valid_vals, c_va, side="right"))
      if cnt_tr_va <= cnt_tr or cnt_tr_va >= len(valid_vals):
        continue
      frac_tr_va = float(cnt_tr_va) / n_valid
      cost = (
          abs(frac_tr - target_train) + abs(frac_tr_va - target_train_val),
          abs(frac_tr - target_train),
          -c_tr,
      )
      if best_cost is None or cost < best_cost:
        best_cost = cost
        best_pair = (c_tr, c_va)
  if best_pair is None:
    raise ValueError(
        "Could not resolve 3 non-empty non-overlapping chronological partitions"
        " from the available timestamps."
    )
  return best_pair


def temporal_split(
    time_series: Any,
    train_ratio: float = 0.8,
    val_ratio: float = 0.0,
    test_ratio: float = 0.2,
    cutoff: Optional[Union[int, float, str]] = None,
    has_subyear_dates: Optional[bool] = None,
) -> Tuple[List[int], List[int], List[int]]:
  """Partitions molecules into chronological sets using a strict `--time-cutoff`.

  Enforces zero timestamp overlap (`max(t_train) < min(t_test)`) and rejects
  boolean (`False`/`True`) or single-timestamp columns so degenerate metadata
  can never masquerade as a prospective Sheridan temporal split.

  Args:
      time_series: Series or sequence of time values (publication year, ISO
        date, or chronological index). Booleans are rejected.
      train_ratio: Target fraction for training set when `cutoff` is omitted.
      val_ratio: Target fraction for validation set when `cutoff` is omitted.
      test_ratio: Target fraction for test set when `cutoff` is omitted.
      cutoff: Optional explicit cutoff value (e.g. `2013` or `"2014-06-30"`).
        Samples with `time <= cutoff` are assigned to train and `time > cutoff`
        to test. When omitted, a strict cutoff is resolved from `time_series` so
        tied timestamps are never split across partitions.
      has_subyear_dates: Optional override indicating `time_series` contains
        sub-year dates so plain-year `cutoff` values resolve to Dec 31.

  Returns:
      Tuple of (train_indices, val_indices, test_indices).
  """
  if cutoff is not None and val_ratio > 0:
    raise ValueError(
        "Explicit `cutoff` / `--time-cutoff` defines a 2-way train/test"
        f" boundary (`val_ratio` must be 0.0, got {val_ratio}); omit"
        " `--time-cutoff` for automatic 3-way train/val/test cutoffs."
    )

  series_name = getattr(time_series, "name", None)
  use_cal_bounds = (
      _is_calendar_col_name(str(series_name))
      if series_name is not None and str(series_name).strip()
      else False
  )
  raw_s = pd.Series(time_series).reset_index(drop=True)
  if raw_s.empty:
    return [], [], []

  numeric_s, detected_subyear = parse_temporal_series(
      raw_s, calendar_bounds=use_cal_bounds
  )
  effective_subyear = (
      detected_subyear
      if has_subyear_dates is None
      else bool(has_subyear_dates or detected_subyear)
  )
  valid_count = int(numeric_s.notna().sum())
  unique_count = int(numeric_s.dropna().nunique())
  if unique_count < 2:
    raise ValueError(
        "Temporal split requires at least 2 distinct valid non-boolean time"
        f" values (found {unique_count} distinct across"
        f" {valid_count}/{len(raw_s)} non-null rows; boolean True/False and"
        " out-of-range placeholders are rejected)."
    )

  if cutoff is not None:
    c_train = parse_temporal_cutoff(cutoff, has_subyear_dates=effective_subyear)
    if _is_subyear_calendar_float(c_train):
      cut_yr = int(math.floor(c_train + 1e-9))
      if c_train < _date_to_fractional_year(cut_yr, 12, 31) - 1e-7 and any(
          _is_coarse_calendar_year(float(v)) and int(round(float(v))) == cut_yr
          for v in numeric_s.dropna().tolist()
      ):
        raise ValueError(
            f"Mid-year temporal cutoff {cutoff!r} cannot disambiguate coarse"
            f" year-only {cut_yr} entries lacking month/day precision; supply"
            f" sub-year ISO dates or use a year-end boundary ('{cut_yr - 1}'"
            f" or '{cut_yr}')."
        )
    c_val: Optional[float] = None
  else:
    c_train, c_val = resolve_temporal_cutoff(
        numeric_s,
        train_ratio=train_ratio,
        val_ratio=val_ratio,
        test_ratio=test_ratio,
    )

  # Order valid indices chronologically (ties broken by original row index).
  valid_pairs = sorted(
      (
          (i, float(numeric_s.iloc[i]))
          for i in range(len(numeric_s))
          if pd.notna(numeric_s.iloc[i])
      ),
      key=lambda pair: (pair[1], pair[0]),
  )
  nan_idx = [i for i in range(len(numeric_s)) if pd.isna(numeric_s.iloc[i])]

  dated_train_idx: List[int] = []
  val_idx: List[int] = []
  test_idx: List[int] = []
  for idx, t_val in valid_pairs:
    if t_val <= c_train:
      dated_train_idx.append(idx)
    elif c_val is not None and t_val <= c_val:
      val_idx.append(idx)
    else:
      test_idx.append(idx)

  if not dated_train_idx or not test_idx:
    min_fmt = format_temporal_value(
        float(numeric_s.min()), has_subyear_dates=effective_subyear
    )
    max_fmt = format_temporal_value(
        float(numeric_s.max()), has_subyear_dates=effective_subyear
    )
    cut_fmt = format_temporal_value(
        c_train, has_subyear_dates=effective_subyear
    )
    raise ValueError(
        f"Temporal cutoff {cut_fmt} produces an empty partition (dated train="
        f"{len(dated_train_idx)}, test={len(test_idx)}); dataset time range is"
        f" [{min_fmt}, {max_fmt}]."
    )

  if nan_idx:
    print(
        f"WARNING: {len(nan_idx)} of {len(numeric_s)} rows remain undated after"
        " temporal parsing/backfill; placing undated rows at the start of the"
        " training split only (never in validation/test).",
        file=sys.stderr,
    )
  train_idx = nan_idx + dated_train_idx
  return train_idx, val_idx, test_idx


def temporal_split_df(
    df: pd.DataFrame,
    time_col: Optional[str] = None,
    train_ratio: float = 0.8,
    val_ratio: float = 0.0,
    test_ratio: float = 0.2,
    cutoff: Optional[Union[int, float, str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
  """DataFrame-level temporal split with automatic true-date backfilling."""
  backfilled_s, has_subyear, _ = backfill_temporal_series(df, time_col=time_col)
  train_idx, val_idx, test_idx = temporal_split(
      backfilled_s,
      train_ratio=train_ratio,
      val_ratio=val_ratio,
      test_ratio=test_ratio,
      cutoff=cutoff,
      has_subyear_dates=has_subyear,
  )
  return (
      df.iloc[train_idx].copy(),
      df.iloc[val_idx].copy(),
      df.iloc[test_idx].copy(),
  )


def main():
  parser = argparse.ArgumentParser(
      description="Split chemical datasets by scaffold or chronological time."
  )
  parser.add_argument(
      "--input_csv",
      "--input-csv",
      dest="input_csv",
      type=str,
      required=True,
      help="Input CSV path.",
  )
  parser.add_argument(
      "--smiles_col",
      "--smiles-col",
      dest="smiles_col",
      type=str,
      default="SMILES",
      help="SMILES column name.",
  )
  parser.add_argument(
      "--split_method",
      "--split-method",
      dest="split_method",
      type=str,
      choices=["scaffold", "time"],
      default="scaffold",
      help="Split strategy.",
  )
  parser.add_argument(
      "--time_col",
      "--time-col",
      dest="time_col",
      type=str,
      default=None,
      help=(
          "Column containing time/date for temporal split (auto-discovers"
          " publication_date, document_year, priority_date, filing_date when"
          " omitted)."
      ),
  )
  parser.add_argument(
      "--time_cutoff",
      "--time-cutoff",
      dest="time_cutoff",
      type=str,
      default=None,
      help=(
          "Explicit cutoff value (e.g. 2013 or 2014-06-30) for strict"
          " non-overlapping temporal split (train <= cutoff, test > cutoff)."
      ),
  )
  parser.add_argument(
      "--train_ratio",
      "--train-ratio",
      dest="train_ratio",
      type=float,
      default=0.8,
      help="Train set ratio.",
  )
  parser.add_argument(
      "--val_ratio",
      "--val-ratio",
      dest="val_ratio",
      type=float,
      default=None,
      help=(
          "Validation set ratio (defaults to 0.1 for scaffold split, 0.0 for"
          " temporal split)."
      ),
  )
  parser.add_argument(
      "--test_ratio",
      "--test-ratio",
      dest="test_ratio",
      type=float,
      default=None,
      help=(
          "Test set ratio (defaults to 0.1 for scaffold split, 0.2 for"
          " temporal split)."
      ),
  )
  parser.add_argument(
      "--strategy",
      type=str,
      choices=["balanced", "greedy"],
      default="balanced",
      help=(
          "Scaffold assignment strategy. 'balanced' (default) keeps each split"
          " representative; 'greedy' reproduces the legacy train-first fill."
      ),
  )
  parser.add_argument(
      "--output_csv",
      "--output-csv",
      "--output",
      dest="output_csv",
      type=str,
      default=None,
      help="Output CSV path to save dataset with assigned 'Split' column.",
  )
  args = parser.parse_args()

  df = pd.read_csv(args.input_csv)
  scaffold_stats = None
  if args.split_method == "time":
    val_ratio = 0.0 if args.val_ratio is None else args.val_ratio
    test_ratio = 0.2 if args.test_ratio is None else args.test_ratio
    if args.time_col is not None and args.time_col not in df.columns:
      print(
          f"Error: Time column '{args.time_col}' not found in"
          f" {args.input_csv}. Available columns: {list(df.columns)}",
          file=sys.stderr,
      )
      sys.exit(1)
    backfilled_s, has_subyear, resolved_col = backfill_temporal_series(
        df, time_col=args.time_col
    )
    if resolved_col is None:
      print(
          "Error: No temporal column or backfillable date metadata exists in"
          f" {args.input_csv}",
          file=sys.stderr,
      )
      sys.exit(1)
    try:
      if args.time_cutoff is not None:
        c_train = parse_temporal_cutoff(
            args.time_cutoff, has_subyear_dates=has_subyear
        )
        c_val = None
      else:
        c_train, c_val = resolve_temporal_cutoff(
            backfilled_s,
            train_ratio=args.train_ratio,
            val_ratio=val_ratio,
            test_ratio=test_ratio,
        )
      train_i, val_i, test_i = temporal_split(
          backfilled_s,
          train_ratio=args.train_ratio,
          val_ratio=val_ratio,
          test_ratio=test_ratio,
          cutoff=args.time_cutoff,
          has_subyear_dates=has_subyear,
      )
    except ValueError as exc:
      print(f"Error: {exc}", file=sys.stderr)
      sys.exit(1)

    # Persist backfilled dates onto output copies when missing/placeholder
    work_df = df.copy()
    if _is_calendar_col_name(resolved_col) or resolved_col == "patent_id":
      out_year_col = (
          resolved_col
          if resolved_col in work_df.columns and resolved_col != "patent_id"
          else "document_year"
      )
      if out_year_col not in work_df.columns:
        work_df[out_year_col] = pd.Series(
            np.nan, index=work_df.index, dtype=object
        )
      cur_parsed, _ = parse_temporal_series(
          work_df[out_year_col], calendar_bounds=True
      )
      fill_mask = cur_parsed.isna() & backfilled_s.notna()
      is_year_only_col = "year" in out_year_col.lower()
      if fill_mask.any():
        work_df[out_year_col] = work_df[out_year_col].astype(object)
        for idx in work_df.index[fill_mask]:
          val_f = float(backfilled_s.loc[idx])
          work_df.at[idx, out_year_col] = (
              int(math.floor(val_f + 1e-9))
              if is_year_only_col
              else format_temporal_value(val_f, has_subyear)
          )
      if out_year_col != "document_year" and "document_year" in work_df.columns:
        doc_yr_parsed, _ = parse_temporal_series(
            work_df["document_year"], calendar_bounds=True
        )
        yr_fill = doc_yr_parsed.isna() & backfilled_s.notna()
        if yr_fill.any():
          work_df["document_year"] = work_df["document_year"].astype(object)
          for idx in work_df.index[yr_fill]:
            work_df.at[idx, "document_year"] = int(
                math.floor(float(backfilled_s.loc[idx]) + 1e-9)
            )

    train_df = work_df.iloc[train_i].copy()
    val_df = work_df.iloc[val_i].copy()
    test_df = work_df.iloc[test_i].copy()
    cut_fmt = format_temporal_value(c_train, has_subyear_dates=has_subyear)
    tr_times = backfilled_s.iloc[train_i].dropna()
    te_times = backfilled_s.iloc[test_i].dropna()
    n_undated_tr = len(train_i) - len(tr_times)
    tr_range = (
        f"[{format_temporal_value(float(tr_times.min()), has_subyear)},"
        f" {format_temporal_value(float(tr_times.max()), has_subyear)}]"
    )
    te_range = (
        f"[{format_temporal_value(float(te_times.min()), has_subyear)},"
        f" {format_temporal_value(float(te_times.max()), has_subyear)}]"
    )
    if c_val is None:
      cutoff_desc = f"--time-cutoff {cut_fmt}"
    else:
      val_cut_fmt = format_temporal_value(c_val, has_subyear_dates=has_subyear)
      cutoff_desc = f"--time-cutoff {cut_fmt}, val_cutoff {val_cut_fmt}"
    undated_note = (
        f" ({n_undated_tr} undated in train)" if n_undated_tr > 0 else ""
    )
    print(
        f"Temporal Split Results for {len(df)} compounds (source:"
        f" '{resolved_col}', {cutoff_desc}):"
    )
    print(
        f"  Train Time Range: {tr_range}{undated_note} | Test Time Range:"
        f" {te_range}"
    )
  else:
    val_ratio = 0.1 if args.val_ratio is None else args.val_ratio
    test_ratio = 0.1 if args.test_ratio is None else args.test_ratio
    s_col = (
        args.smiles_col
        if args.smiles_col in df.columns
        else (
            "canonical_smiles"
            if "canonical_smiles" in df.columns
            else args.smiles_col
        )
    )
    if s_col not in df.columns:
      print(
          f"Error: SMILES column '{args.smiles_col}' not found in "
          f"{args.input_csv}. Available columns: {list(df.columns)}",
          file=sys.stderr,
      )
      sys.exit(1)
    smiles = df[s_col].tolist()
    n_unparseable = sum(
        1
        for s in smiles
        if not pd.notna(s) or Chem.MolFromSmiles(str(s)) is None
    )
    if n_unparseable:
      print(
          f"WARNING: {n_unparseable} of {len(smiles)} rows have "
          f"unparseable SMILES in column '{s_col}'. They have no "
          "Bemis-Murcko scaffold and are grouped together under "
          "'INVALID_SMILES', so they will all land in the same "
          "partition. Clean the input if this matters.",
          file=sys.stderr,
      )
    train_i, val_i, test_i = scaffold_split(
        smiles,
        args.train_ratio,
        val_ratio,
        test_ratio,
        strategy=args.strategy,
    )
    train_df = df.iloc[train_i].copy()
    val_df = df.iloc[val_i].copy()
    test_df = df.iloc[test_i].copy()
    print(
        f"Scaffold Split Results for {len(df)} compounds (strategy:"
        f" {args.strategy}):"
    )
    scaffold_stats = summarize_split(smiles, train_i, val_i, test_i)

  print(
      f"  Train Set: {len(train_df):5d} compounds"
      f" ({len(train_df)/len(df)*100:.1f}%)"
  )
  print(
      f"  Val Set:   {len(val_df):5d} compounds"
      f" ({len(val_df)/len(df)*100:.1f}%)"
  )
  print(
      f"  Test Set:  {len(test_df):5d} compounds"
      f" ({len(test_df)/len(df)*100:.1f}%)"
  )

  if scaffold_stats:
    print(
        "\n  Scaffold distribution (a healthy split keeps populous scaffolds in"
        " every partition):"
    )
    print(f"    {'split':6s} {'scaffolds':>10s} {'largest':>9s} {'median':>7s}")
    for name in ("train", "val", "test"):
      s = scaffold_stats[name]
      print(
          f"    {name:6s} {s['scaffolds']:10d} {s['largest_scaffold']:9d}"
          f" {s['median_scaffold']:7d}"
      )

  if args.output_csv:
    train_df = train_df.copy()
    val_df = val_df.copy()
    test_df = test_df.copy()
    train_df["Split"] = "train"
    val_df["Split"] = "val"
    test_df["Split"] = "test"
    combined_df = pd.concat([train_df, val_df, test_df], ignore_index=True)
    combined_df.to_csv(args.output_csv, index=False)
    print(f"Saved split dataset to {args.output_csv}")


if __name__ == "__main__":
  main()
