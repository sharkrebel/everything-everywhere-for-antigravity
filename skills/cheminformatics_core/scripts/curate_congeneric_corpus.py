# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "pandas",
#     "numpy",
#     "matplotlib",
#     "scipy",
#     "seaborn",
#     "polite-http",
# ]
# ///

"""Congeneric Series & Literature Corpus Aggregation Engine.

Collapses multi-step manual scraping, ChEMBL pagination, patent table
ingestion, descriptor calculation, and property visualization into a
single-call reproducible CLI.

What this script does:
1. Ingests ChEMBL document activities and targets with automated pagination,
   retry/backoff, and automatic NCBI PubChem BioAssay (elink) fallback when
   ChEMBL/EBI is unreachable.
2. Ingests patent tables and chemical series CSVs, then concatenates them.
3. Canonicalizes SMILES via RDKit and deduplicates compounds while preserving
   distinct stereoisomers (`@`/`@@`), separated chiral peaks (`Peak 1`/`Peak
   2`), and single-stereocenter enantiomer pairs across patent examples (e.g.,
   active `(R)` vs inactive `(S)` controls with distinct potencies), merging
   provenance (`all_sources`, `patent_example_id`, `all_endpoints`) across
   colliding records.
4. Computes canonical RDKit descriptors (including `is_macrocycle`,
   `MaxRingSize`, `CNS_MPO`, `SA_Score`, and `calculated_inchikey`) via
   `descriptors.compute_standard_descriptors`.
5. Renders a configurable multi-panel chemical space distribution figure
   (`--panels` supporting `scatter`, `hist`, `box`, `kde`, `bar` templates over
   any descriptor or potency feature) colored consistently across series or
   categorical cohorts (`--group-by`) using the publication palette with
   automatic integer tick detection. Supports `--plot-only` for instant
   re-plotting of existing curated CSVs.
6. Exports the aggregated table as CSV and, optionally, as SDF and interactive
   HTML (`--output-html`).
"""

from __future__ import annotations

import argparse
from collections.abc import Hashable, Sequence
import hashlib
import io
import json
import os
import re
import sys
import textwrap
import time
from typing import Any, Final, Optional, Tuple, Union

import matplotlib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy import stats as scipy_stats

matplotlib.use("Agg")
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

_script_dir = os.path.dirname(os.path.abspath(__file__))
_general_scripts = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "..",
        "rdkit_general",
        "scripts",
    )
)
if _script_dir not in sys.path:
  sys.path.insert(0, _script_dir)
if _general_scripts not in sys.path:
  sys.path.append(_general_scripts)

from compat import (
    PUBLICATION_CYCLE,
    PUBLICATION_PALETTE,
    autodetect_integer_axes,
    set_safe_plot_style,
)
from draw_molecules import compute_clean_2d_coords, render_grid_molecules
from interactive_grid_viewer import create_molecular_grid
from skill_paths import add_skill_scripts

add_skill_scripts("cheminformatics_doc2mol")

import descriptors
from paper_fetcher import (
    _build_endpoint_column_key,
    fetch_chembl_document_activities,
)


def _resolve_grouping(
    df: pd.DataFrame, group_by: Optional[str] = None
) -> Tuple[Optional[str], list[tuple[str, pd.Series, str]]]:
  """Resolves the categorical grouping column and assigns Generative UI cycle colors.

  Args:
    df: Input DataFrame containing molecular properties and metadata.
    group_by: Optional column name to group by, or 'none' to disable grouping.
      If omitted, auto-detects from standard categorical/series columns.

  Returns:
    Tuple of (resolved_column_name, list of (label, boolean_mask, color_hex)
    tuples).
  """
  if group_by and group_by.lower() == "none":
    full_mask = pd.Series(True, index=df.index)
    return None, [("All Compounds", full_mask, PUBLICATION_CYCLE[0])]

  is_single_patent_only = (
      "patent_id" in df.columns
      and df["patent_id"].dropna().nunique() <= 1
      and (
          "source" not in df.columns
          or not df["source"]
          .dropna()
          .astype(str)
          .str.contains(r"chembl|literature|paper|journal|jmc", case=False)
          .any()
      )
  )
  default_cands = [
      "series",
      "series_name",
      "cohort",
      "Matched_Core",
      "scaffold_label",
  ]
  if is_single_patent_only:
    default_cands.extend(["activity_bin", "source"])
  else:
    default_cands.extend(["source", "activity_bin"])

  candidates = [group_by] if group_by else default_cands
  chosen_col = None
  for col in candidates:
    if col and col in df.columns:
      n_uniq = df[col].dropna().nunique()
      if 2 <= n_uniq <= 8 or (group_by and n_uniq >= 1):
        chosen_col = col
        break

  if not chosen_col:
    full_mask = pd.Series(True, index=df.index)
    return None, [("All Compounds", full_mask, PUBLICATION_CYCLE[0])]

  vc = df[chosen_col].value_counts(dropna=True)
  ordered_vals = list(vc.index)
  if chosen_col == "activity_bin":
    ordered_vals.sort(
        key=lambda v: (
            0
            if re.match(r"^[A-F](?:\s*\(.*\))?$", str(v).strip(), re.I)
            else 1,
            str(v).strip().upper(),
        )
    )
  else:
    ordered_vals.sort(
        key=lambda v: (
            1
            if str(v).strip().lower() in ("other", "unclassified", "unknown")
            else 0
        )
    )
  groups = []
  for idx, val in enumerate(ordered_vals):
    mask = df[chosen_col] == val
    color = PUBLICATION_CYCLE[idx % len(PUBLICATION_CYCLE)]
    groups.append((str(val), mask, color))
  return chosen_col, groups


FEATURE_METADATA: dict[str, dict[str, object]] = {
    "MW": {
        "label": "Molecular Weight (Da)",
        "short_title": "Molecular Weight (MW)",
        "threshold": (500.0, "Ro5 MW (500)", "muted", "--"),
        "show_mean": False,
    },
    "cLogP": {
        "label": "cLogP (Lipophilicity)",
        "short_title": "Lipophilicity (cLogP)",
        "threshold": (5.0, "Ro5 cLogP (5.0)", "muted", ":"),
        "show_mean": False,
    },
    "TPSA": {
        "label": "Topological Polar Surface Area (Å²)",
        "short_title": "Polar Surface Area (TPSA)",
        "threshold": (140.0, "Oral TPSA (140 Å²)", "muted", "--"),
        "show_mean": True,
    },
    "RotB": {
        "label": "Rotatable Bond Count (Flexibility)",
        "short_title": "Rotatable Bonds (Conformational Rigidity)",
        "threshold": None,
        "show_mean": False,
    },
    "HBD": {
        "label": "Hydrogen Bond Donors (HBD)",
        "short_title": "H-Bond Donors (HBD)",
        "threshold": (5.0, "Ro5 HBD (5)", "muted", "--"),
        "show_mean": False,
    },
    "HBA": {
        "label": "Hydrogen Bond Acceptors (HBA)",
        "short_title": "H-Bond Acceptors (HBA)",
        "threshold": (10.0, "Ro5 HBA (10)", "muted", "--"),
        "show_mean": False,
    },
    "Fsp3": {
        "label": "Fraction sp³ Carbons (Fsp3)",
        "short_title": "3D Character (Fsp3)",
        "threshold": (0.42, "Drug-like Fsp3 (0.42)", "success", "--"),
        "show_mean": True,
    },
    "QED": {
        "label": "Quantitative Estimate of Drug-likeness (QED)",
        "short_title": "Drug-likeness (QED)",
        "threshold": (0.67, "High QED (>= 0.67)", "success", "--"),
        "show_mean": True,
    },
    "SA_Score": {
        "label": "Synthetic Accessibility Score (1=Easy, 10=Hard)",
        "short_title": "Synthetic Accessibility (SA Score)",
        "threshold": (6.0, "Tractable Synthesis (<= 6.0)", "warning", "--"),
        "show_mean": True,
    },
    "CNS_MPO": {
        "label": "Pfizer CNS MPO Score (0.0 - 6.0)",
        "short_title": "Central Nervous System Penetrance (CNS MPO)",
        "threshold": (
            4.0,
            "CNS Drug-like Desirability (>= 4.0)",
            "success",
            "--",
        ),
        "show_mean": True,
    },
    "HeavyAtoms": {
        "label": "Heavy Atom Count",
        "short_title": "Heavy Atom Count",
        "threshold": None,
        "show_mean": False,
    },
    "Rings": {
        "label": "Ring Count",
        "short_title": "Total Ring Count",
        "threshold": None,
        "show_mean": False,
    },
    "AromaticRings": {
        "label": "Aromatic Ring Count",
        "short_title": "Aromatic Ring Count",
        "threshold": (3.0, "Aromatic Ring Limit (3)", "warning", "--"),
        "show_mean": False,
    },
    "AliphaticRings": {
        "label": "Aliphatic Ring Count",
        "short_title": "Aliphatic Ring Count",
        "threshold": None,
        "show_mean": False,
    },
    "MaxRingSize": {
        "label": "Largest Ring Size (Atoms)",
        "short_title": "Largest Ring Size",
        "threshold": None,
        "show_mean": False,
    },
}

_POTENCY_CANDIDATES = (
    "pchembl_value",
    "pActivity",
    "pactivity",
    "pic50",
    "pki",
    "standard_value",
)


def _human_assay_label(col_name: str) -> str:
  """Formats an assay column identifier into a concise human-readable title."""
  s = str(col_name).strip()
  if s in ("activity_bin", "activity_value"):
    return s
  s = re.sub(r"^patent_assay_", "", s, flags=re.IGNORECASE)
  tokens = [t for t in s.split("_") if t]
  upper_acronyms = {
      "ic50": "IC50",
      "ec50": "EC50",
      "ki": "Ki",
      "kd": "Kd",
      "wt": "WT",
      "nm": "nM",
      "um": "uM",
      "elisa": "ELISA",
      "htrf": "HTRF",
      "trfret": "TR-FRET",
  }
  formatted: list[str] = []
  for tok in tokens:
    low = tok.lower()
    if low in upper_acronyms:
      formatted.append(upper_acronyms[low])
    elif re.match(r"^[a-z]{2,5}\d+[a-z0-9]*$", low):
      formatted.append(tok.upper())
    else:
      formatted.append(tok.capitalize())
  return " ".join(formatted) if formatted else str(col_name)


def _count_discrete_potency_values(df: pd.DataFrame, col: str) -> int:
  """Counts non-empty categorical potency values in a candidate column."""
  if col not in df.columns:
    return 0
  s = (
      df[col]
      .dropna()
      .astype(str)
      .str.strip()
      .replace({"nan": "", "none": "", "nd": "", "n/a": "", "-": ""})
  )
  return int(s.ne("").sum())


def _is_discrete_potency_col(df: pd.DataFrame, col: Optional[str]) -> bool:
  """Returns True if col holds categorical/binned potency labels rather than continuous floats."""
  if not col or col not in df.columns:
    return False
  n_num = int(pd.to_numeric(df[col], errors="coerce").notna().sum())
  n_str = _count_discrete_potency_values(df, col)
  if n_str == 0 or n_num > n_str * 0.5:
    return False
  if col == "activity_bin":
    return True
  sample = (
      df[col].dropna().astype(str).str.strip().loc[lambda s: s.ne("")].head(20)
  )
  disc_pat = (
      r"^(?:[A-F][+-]?(?:\s*\(.*\))?|\+{1,6}|\*{1,6}|[<>≤≥~].+|"
      r"\d+(?:\.\d+)?\s*(?:pM|nM|uM|µM|μM|mM)?\s*(?:[-–—]|to)\s*\d+.*)$"
  )
  return bool(
      not sample.empty and sample.str.match(disc_pat, case=False).mean() >= 0.5
  )


def _resolve_feature_col(
    df: pd.DataFrame, feature: str
) -> Tuple[Optional[str], dict[str, Any]]:
  """Resolves a feature name or shorthand alias (e.g.

  'potency') to a DataFrame column and metadata.
  """
  feat_clean = feature.strip()
  if (
      feat_clean.lower() in ("potency", "activity")
      or feat_clean in _POTENCY_CANDIDATES
  ):
    candidates = (
        [feat_clean] if feat_clean in df.columns else list(_POTENCY_CANDIDATES)
    )
    resolved = None
    n_quant = 0
    for cand in candidates:
      if cand in df.columns:
        cnt_q = int(pd.to_numeric(df[cand], errors="coerce").notna().sum())
        if cnt_q > n_quant:
          resolved = cand
          n_quant = cnt_q

    disc_candidates = [
        c
        for c in df.columns
        if str(c).startswith("patent_assay_") and not str(c).endswith("_range")
    ] + ["activity_bin", "activity_value"]
    best_disc: Optional[str] = None
    n_disc = 0
    for disc_cand in disc_candidates:
      if _is_discrete_potency_col(df, disc_cand):
        cnt_d = _count_discrete_potency_values(df, disc_cand)
        if cnt_d > n_disc:
          best_disc = disc_cand
          n_disc = cnt_d

    # Prefer discrete potency bins when 'potency'/'activity' is requested (or no
    # quantitative column exists) and discrete bins cover more compounds than a
    # sparse/single-row quantitative column.
    prefer_discrete = best_disc is not None and (
        resolved is None
        or (
            feat_clean.lower() in ("potency", "activity")
            and n_disc > n_quant
            and n_quant < max(5, int(0.10 * max(len(df), 1)))
        )
    )
    if prefer_discrete and best_disc is not None:
      human_lbl = _human_assay_label(best_disc)
      return best_disc, {
          "label": f"Discrete Potency Bin ({human_lbl})",
          "short_title": f"Discrete Bioactivity ({human_lbl})",
          "threshold": None,
          "show_mean": False,
          "show_median": False,
          "is_discrete_potency": True,
          "empty_msg": "No quantitative or discrete potency records annotated",
      }
    col_label = resolved or "pchembl_value"
    meta: dict[str, Any] = {
        "label": f"Potency ({col_label} / -log M)",
        "short_title": "Bioactivity Distribution",
        "threshold": (9.0, "Sub-nM (pAct >= 9.0)", "danger", ":"),
        "show_mean": False,
        "show_median": True,
        "empty_msg": "No quantitative potency records annotated",
    }
    return resolved, meta

  resolved = None
  if feat_clean in df.columns:
    resolved = feat_clean
  else:
    col_map = {c.lower(): c for c in df.columns}
    if feat_clean.lower() in col_map:
      resolved = col_map[feat_clean.lower()]

  canon_key = resolved or feat_clean
  if (
      resolved
      and canon_key not in FEATURE_METADATA
      and _is_discrete_potency_col(df, resolved)
  ):
    human_lbl = _human_assay_label(resolved)
    return resolved, {
        "label": f"Discrete Potency Bin ({human_lbl})",
        "short_title": human_lbl,
        "threshold": None,
        "show_mean": False,
        "show_median": False,
        "is_discrete_potency": True,
        "empty_msg": f"No records annotated for {resolved}",
    }

  base_meta = FEATURE_METADATA.get(
      canon_key,
      {
          "label": _human_assay_label(canon_key),
          "short_title": _human_assay_label(canon_key),
          "threshold": None,
          "show_mean": False,
      },
  )
  return resolved, dict(base_meta)


def _expand_panel_specs(
    panels: Sequence[Union[str, dict[str, object]]],
) -> list[Union[str, dict[str, object]]]:
  """Expands comma- or semicolon-separated panel spec strings into individual panel specs."""
  expanded: list[Union[str, dict[str, object]]] = []
  for item in panels:
    if isinstance(item, dict):
      expanded.append(item)
      continue
    s = str(item).strip()
    if not s:
      continue
    if ("," in s or ";" in s) and ":" in s:
      tokens = [t.strip() for t in re.split(r"[;,]", s) if t.strip()]
      if len(tokens) > 1 and all(
          ":" in t or re.search(r"_vs_|\s+vs\s+", t, re.I) for t in tokens
      ):
        expanded.extend(tokens)
        continue
    expanded.append(s)
  return expanded


def _parse_panel_spec(spec: Union[str, dict[str, object]]) -> dict[str, object]:
  """Parses a panel specification string or dictionary into a normalized configuration dict."""
  supported_types = {"scatter", "hist", "box", "kde", "bar"}
  if isinstance(spec, dict):
    ptype = str(spec.get("type", "hist")).lower()
    if ptype not in supported_types:
      raise ValueError(
          f"Unsupported plot type '{ptype}' in panel spec. Must be one of"
          f" {sorted(supported_types)}."
      )
    return dict(spec, type=ptype)

  raw = str(spec).strip()
  parts = [p.strip() for p in raw.split(":") if p.strip()]
  if not parts:
    raise ValueError(f"Empty panel specification: '{spec}'")

  first_lower = parts[0].lower()
  if first_lower in supported_types:
    ptype = first_lower
    rest = parts[1:]
  elif len(parts) >= 2 and parts[1].lower() in supported_types:
    ptype = parts[1].lower()
    rest = [parts[0]] + parts[2:]
  elif len(parts) >= 3 and parts[2].lower() in supported_types:
    ptype = parts[2].lower()
    rest = parts[:2] + parts[3:]
  else:
    ptype = (
        "scatter"
        if (len(parts) >= 2 or re.search(r"_vs_|\s+vs\s+", parts[0], re.I))
        else "hist"
    )
    rest = parts

  if (
      ptype == "scatter"
      and len(rest) >= 1
      and re.search(r"_vs_|\s+vs\s+", rest[0], re.I)
  ):
    vs_parts = [
        vp.strip()
        for vp in re.split(r"_vs_|\s+vs\s+", rest[0], maxsplit=1, flags=re.I)
        if vp.strip()
    ]
    if len(vs_parts) == 2:
      rest = vs_parts + rest[1:]

  if ptype == "scatter":
    if len(rest) < 2:
      raise ValueError(
          f"Scatter panel spec '{spec}' requires two features (e.g."
          " 'scatter:MW:cLogP')."
      )
    out: dict[str, object] = {"type": "scatter", "x": rest[0], "y": rest[1]}
    if len(rest) >= 3:
      out["title"] = rest[2]
    return out
  else:
    if not rest:
      raise ValueError(
          f"Panel spec '{spec}' requires a feature name (e.g."
          f" '{ptype}:SA_Score')."
      )
    out = {"type": ptype, "x": rest[0]}
    if len(rest) >= 2:
      out["title"] = rest[1]
    return out


def _format_missing_panels_help(
    df: pd.DataFrame,
    output_csv: str = "curated_corpus.csv",
    output_plot: str = "chemical_space_distributions.png",
) -> str:
  """Builds a dataset-aware helper error message when --panels is omitted."""
  avail_desc = [c for c in FEATURE_METADATA if c in df.columns]
  bio_items: list[Tuple[int, str]] = []
  skipped_constant: list[str] = []
  for c in df.columns:
    cs = str(c)
    if cs.endswith("_range") or cs == "activity_type":
      continue
    if (
        cs in _POTENCY_CANDIDATES
        or cs.startswith(("patent_assay_", "patent_table_col_"))
        or cs in ("activity_bin", "activity_value")
        or cs.endswith("_nM")
        or re.search(
            r"(?:^|_)(p?(?:Ki|IC50|Kd|EC50|Kis)|Ratio)(?:_|$)",
            cs,
            re.IGNORECASE,
        )
    ):
      num_s = pd.to_numeric(df[c], errors="coerce")
      n_num = int(num_s.notna().sum())
      n_uniq = int(num_s.dropna().nunique())
      if _is_discrete_potency_col(df, cs):
        n_disc = _count_discrete_potency_values(df, cs)
        top_bins = (
            df[c]
            .dropna()
            .astype(str)
            .str.strip()
            .value_counts()
            .head(5)
            .to_dict()
        )
        bins_str = ", ".join(f"{k}:{v}" for k, v in top_bins.items())
        bio_items.append((
            n_disc,
            (
                f"    - {cs} (n={n_disc}, discrete bins: {bins_str}) -> use"
                f" 'bar:{cs}'"
            ),
        ))
      elif n_num > 0:
        if n_uniq < 2 and len(df) > 1:
          const_val = num_s.dropna().iloc[0]
          skipped_constant.append(f"{cs} (constant={const_val:g}, n={n_num})")
          continue
        bio_items.append((
            n_num,
            (
                f"    - {cs} (n={n_num}, continuous numeric) -> use 'hist:{cs}'"
                f" or 'box:{cs}'"
            ),
        ))

  bio_items.sort(key=lambda kv: -kv[0])
  if any(cnt >= 2 for cnt, _ in bio_items):
    bio_items = [item for item in bio_items if item[0] >= 2]
  bio_lines = [line for _, line in bio_items[:12]]
  bio_block = (
      "\n".join(bio_lines)
      if bio_lines
      else "    - (no bioactivity columns detected)"
  )
  const_note = (
      "  [WARNING] Omitted constant/single-valued numeric column(s) from panel"
      f" suggestions: {', '.join(skipped_constant)}\n"
      if skipped_constant
      else ""
  )
  desc_str = ", ".join(avail_desc) if avail_desc else "MW, cLogP, TPSA, ..."
  return (
      "\n[ERROR] Plotting requires explicit --panels (no default panels are"
      " assumed so plots match the dataset and scientific question).\n"
      f"  Curated dataset saved to: {output_csv} (N={len(df)} compounds).\n"
      f"  Available Physicochemical Descriptors:\n    {desc_str}\n"
      f"  Available Bioactivity / Assay Columns in this dataset:\n{bio_block}\n"
      f"{const_note}"
      "  Panel Syntax: 'type:feature[:feature2][:Custom Title]' where type is"
      " one of: scatter, hist, box, kde, bar.\n"
      "  Re-run instantly in --plot-only mode with your chosen panels:\n"
      f"    uv run {os.path.abspath(__file__)} --plot-only"
      f" --input-csv {output_csv} --output-plot {output_plot}"
      " --panels 'scatter:MW:cLogP' 'hist:TPSA' '<panel3>' '<panel4>'\n"
  )


def _draw_threshold_line(
    ax: plt.Axes,
    threshold: Any,
    axis: str = "x",
) -> None:
  """Draws a vertical or horizontal reference threshold line using Generative UI palette tokens."""
  if not isinstance(threshold, tuple) or len(threshold) != 4:
    return
  val, label, color_key, style = threshold
  color_hex = PUBLICATION_PALETTE.get(color_key, color_key)
  if axis == "x":
    ax.axvline(
        val,
        color=color_hex,
        linestyle=style,
        linewidth=1.4,
        alpha=0.75,
        label=label,
    )
  else:
    ax.axhline(
        val,
        color=color_hex,
        linestyle=style,
        linewidth=1.4,
        alpha=0.75,
        label=label,
    )


def _add_reference_legend_with_headroom(
    ax: plt.Axes, is_zero_bottom: bool = False
) -> None:
  """Draws a compact legend for threshold/reference lines only and adds top y-headroom."""
  handles, labels = ax.get_legend_handles_labels()
  if not handles:
    return
  ymin, ymax = ax.get_ylim()
  span = max(ymax - ymin, 1e-3)
  headroom = (0.14 + 0.06 * min(len(handles), 3)) * span
  if is_zero_bottom:
    ax.set_ylim(0.0, ymax + headroom)
  else:
    ax.set_ylim(ymin, ymax + headroom)
  ax.legend(
      handles,
      labels,
      loc="upper right",
      frameon=True,
      framealpha=0.92,
      fontsize=8.5,
      ncol=min(len(handles), 2),
  )


def _format_box_xtick_label(
    name: str, count: int, mean_val: Optional[float] = None, max_width: int = 15
) -> str:
  """Formats and wraps a cohort x-axis tick label so adjacent boxplots never collide."""
  clean = re.sub(r"_+", " ", str(name).strip())
  # Split parenthetical qualifiers onto a new line first for cleaner wrapping
  if " (" in clean and clean.endswith(")"):
    head, tail = clean.split(" (", 1)
    wrapped_head = textwrap.fill(head.strip(), width=max_width)
    wrapped_tail = textwrap.fill("(" + tail.strip(), width=max_width + 3)
    base = f"{wrapped_head}\n{wrapped_tail}"
  else:
    base = textwrap.fill(clean, width=max_width)

  if mean_val is not None and not np.isnan(mean_val):
    return f"{base}\n(n={count}, μ={mean_val:.2f})"
  return f"{base}\n(n={count})"


def _draw_empty_panel_message(ax: plt.Axes, msg: str) -> None:
  """Renders a centered placeholder message on an empty subplot panel."""
  ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes)


def _render_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    spec: dict[str, object],
    groups: list[tuple[str, pd.Series, str]],
    panel_letter: str,
) -> None:
  """Renders a single configured panel template onto the provided Matplotlib Axes."""
  ptype = str(spec.get("type", "hist")).lower()
  x_req = str(spec.get("x", ""))
  x_col, x_meta = _resolve_feature_col(df, x_req)

  raw_title = str(spec.get("title") or x_meta.get("short_title") or x_req)
  if not raw_title.startswith(f"{panel_letter}. "):
    title_str = f"{panel_letter}. {raw_title}"
  else:
    title_str = raw_title
  ax.set_title(title_str, fontsize=12, fontweight="bold", pad=8)

  if ptype == "scatter":
    y_req = str(spec.get("y", ""))
    y_col, y_meta = _resolve_feature_col(df, y_req)
    if not spec.get("title"):
      ax.set_title(
          f"{panel_letter}. {x_meta.get('short_title', x_req)} vs"
          f" {y_meta.get('short_title', y_req)}",
          fontsize=12,
          fontweight="bold",
          pad=8,
      )
    if not x_col or not y_col:
      _draw_empty_panel_message(
          ax, f"Columns '{x_req}' / '{y_req}' not available"
      )
      return

    has_points = False
    for _, mask, color in groups:
      sub = df.loc[mask].copy()
      sub[x_col] = pd.to_numeric(sub[x_col], errors="coerce")
      sub[y_col] = pd.to_numeric(sub[y_col], errors="coerce")
      sub = sub.dropna(subset=[x_col, y_col])
      if not sub.empty:
        has_points = True
        ax.scatter(
            sub[x_col],
            sub[y_col],
            color=color,
            alpha=0.75,
            edgecolors=PUBLICATION_PALETTE["foreground"]
            if len(groups) > 1
            else "none",
            linewidths=0.4 if len(groups) > 1 else 0.0,
            s=42,
        )

    if not has_points:
      _draw_empty_panel_message(
          ax, f"No numeric records for {x_col} vs {y_col}"
      )
      return

    hl_list = spec.get("_highlight_compounds")
    if isinstance(hl_list, (list, tuple)) and hl_list:
      hl_markers = ["*", "D", "^", "s", "P", "X"]
      hl_colors = [
          "#E11D48",
          "#16A34A",
          "#D97706",
          "#7C3AED",
          "#0284C7",
          "#DC2626",
      ]
      id_cols = [
          c
          for c in (
              "compound_name",
              "pref_name",
              "paper_compound_id",
              "compound_id",
              "patent_example_id",
              "chembl_id",
          )
          if c in df.columns
      ]
      seen_hl_indices: set[Any] = set()
      marker_idx = 0
      for h_query in hl_list:
        hq = str(h_query).strip()
        if not hq:
          continue
        hq_norm = re.sub(r"[\s_-]+", "_", hq.lower()).strip("_")
        matched_row = None
        matched_idx = None
        for r_idx, r_cand in df.iterrows():
          for ic in id_cols:
            v_str = str(r_cand.get(ic) or "").strip()
            if not v_str or v_str.lower() in ("nan", "none"):
              continue
            v_norm = re.sub(r"[\s_-]+", "_", v_str.lower()).strip("_")
            if hq_norm == v_norm or bool(
                re.search(
                    rf"(?:^|[:,;(\s]){re.escape(hq_norm)}(?:$|[,;)\s])",
                    v_norm,
                )
            ):
              xv = pd.to_numeric(r_cand.get(x_col), errors="coerce")
              yv = pd.to_numeric(r_cand.get(y_col), errors="coerce")
              if pd.notna(xv) and pd.notna(yv):
                matched_row = (float(xv), float(yv))
                matched_idx = r_idx
                break
          if matched_row is not None:
            break
        if matched_row is not None and matched_idx not in seen_hl_indices:
          seen_hl_indices.add(matched_idx)
          ax.scatter(
              [matched_row[0]],
              [matched_row[1]],
              color=hl_colors[marker_idx % len(hl_colors)],
              marker=hl_markers[marker_idx % len(hl_markers)],
              s=180,
              edgecolors="#0F172A",
              linewidths=1.1,
              zorder=10,
              label=hq,
          )
          marker_idx += 1

    _draw_threshold_line(ax, x_meta.get("threshold"), axis="x")
    _draw_threshold_line(ax, y_meta.get("threshold"), axis="y")
    ax.set_xlabel(
        str(spec.get("xlabel") or x_meta.get("label") or x_col),
        fontsize=11,
        fontweight="bold",
    )
    ax.set_ylabel(
        str(spec.get("ylabel") or y_meta.get("label") or y_col),
        fontsize=11,
        fontweight="bold",
    )
    ax.grid(True, linestyle="--", alpha=0.5)
    autodetect_integer_axes(
        ax,
        x_data=pd.to_numeric(df[x_col], errors="coerce"),
        y_data=pd.to_numeric(df[y_col], errors="coerce"),
        is_histogram=False,
    )
    _add_reference_legend_with_headroom(ax, is_zero_bottom=False)
    return

  if ptype in ("hist", "kde", "box") and x_meta.get("is_discrete_potency"):
    ptype = "bar"

  if ptype == "bar":
    if not x_col or x_col not in df.columns:
      _draw_empty_panel_message(ax, f"Column '{x_req}' not available")
      return
    vc_cats = df[x_col].dropna().astype(str).value_counts()
    if x_meta.get("is_discrete_potency") or x_col in (
        "activity_bin",
        "activity_value",
    ):
      top_cats = sorted(
          vc_cats.head(8).index.tolist(),
          key=lambda s: (
              0
              if re.match(r"^[A-F](?:\s*\(.*\))?$", s.strip(), re.I)
              else (1 if re.match(r"^(\+{1,6}|\*{1,6})$", s.strip()) else 2),
              s.strip().upper()
              if not re.match(r"^(\+{1,6}|\*{1,6})$", s.strip())
              else -len(s.strip()),
          ),
      )
    else:
      top_cats = vc_cats.head(8).index.tolist()
    if not top_cats:
      _draw_empty_panel_message(ax, f"No records for '{x_col}'")
      return
    disp_cats = list(top_cats)
    rng_col = None
    if f"{x_col}_range" in df.columns:
      rng_col = f"{x_col}_range"
    elif (
        x_col in ("activity_bin", "activity_value")
        and "activity_range" in df.columns
    ):
      rng_col = "activity_range"
    elif x_meta.get("is_discrete_potency") and "activity_range" in df.columns:
      rng_col = "activity_range"
    if rng_col:
      bin_to_rng: dict[str, str] = {}
      for b_v, r_v in zip(df[x_col], df[rng_col]):
        if pd.notna(b_v) and pd.notna(r_v) and str(r_v).strip():
          bin_to_rng.setdefault(str(b_v).strip(), str(r_v).strip())
      disp_cats = [
          f"{c}\n({bin_to_rng[c]})"
          if c in bin_to_rng and bin_to_rng[c] not in c
          else c
          for c in top_cats
      ]
    x_pos = np.arange(len(top_cats))
    bottoms = np.zeros(len(top_cats))
    for _, mask, color in groups:
      sub_s = df.loc[mask, x_col].dropna().astype(str)
      counts = np.array(
          [int((sub_s == c).sum()) for c in top_cats], dtype=float
      )
      if counts.sum() > 0:
        ax.bar(
            x_pos,
            counts,
            bottom=bottoms,
            color=color,
            edgecolor=PUBLICATION_PALETTE["foreground"],
            linewidth=0.6,
            width=0.65,
        )
        bottoms += counts
    ax.set_xticks(x_pos)
    ax.set_xticklabels(
        disp_cats,
        rotation=15 if any("\n" in dc for dc in disp_cats) else 25,
        ha="center" if any("\n" in dc for dc in disp_cats) else "right",
        fontsize=8.5,
    )
    ax.set_xlabel(
        str(spec.get("xlabel") or x_meta.get("label") or x_col),
        fontsize=11,
        fontweight="bold",
    )
    ax.set_ylabel("Compound Count", fontsize=11, fontweight="bold")
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    autodetect_integer_axes(ax, y_data=bottoms, is_histogram=True)
    _add_reference_legend_with_headroom(ax, is_zero_bottom=True)
    return

  # Numeric 1D distributions: hist, box, kde
  if not x_col or x_col not in df.columns:
    empty_msg = str(
        x_meta.get("empty_msg") or f"{x_req} not calculated or annotated"
    )
    _draw_empty_panel_message(ax, empty_msg)
    return

  series_list: list[pd.Series] = []
  colors_list: list[str] = []
  raw_names: list[str] = []
  show_mean = bool(x_meta.get("show_mean", False))

  for label, mask, color in groups:
    vals = pd.to_numeric(df.loc[mask, x_col], errors="coerce").dropna()
    if not vals.empty:
      series_list.append(vals)
      raw_names.append(label)
      colors_list.append(color)

  if not series_list:
    empty_msg = str(
        x_meta.get("empty_msg") or f"No numeric records annotated for {x_col}"
    )
    _draw_empty_panel_message(ax, empty_msg)
    return

  all_vals = pd.concat(series_list)

  if ptype == "hist":
    arr = all_vals.to_numpy(dtype=float)
    is_discrete_int = (
        len(arr) > 0
        and np.all(np.isclose(arr, np.round(arr)))
        and (np.max(arr) - np.min(arr)) <= 35
    )
    if is_discrete_int:
      min_i = int(np.min(arr))
      max_i = int(np.max(arr))
      bins: Union[int, np.ndarray] = np.arange(min_i - 0.5, max_i + 1.5, 1.0)
    elif x_col == "CNS_MPO":
      bins = np.linspace(0.5, 6.0, 12)
    else:
      bins = (
          np.linspace(float(all_vals.min()), float(all_vals.max()), 16)
          if all_vals.nunique() > 1
          else 10
      )
    ax.hist(
        series_list,
        bins=bins,
        color=colors_list,
        stacked=True,
        edgecolor=PUBLICATION_PALETTE["foreground"],
        linewidth=0.6,
        alpha=0.85,
        rwidth=0.85,
    )
    if x_meta.get("show_median"):
      ax.axvline(
          all_vals.median(),
          color=PUBLICATION_PALETTE["foreground"],
          linestyle="--",
          linewidth=1.2,
          label=f"Median ({all_vals.median():.2f})",
      )
    _draw_threshold_line(ax, x_meta.get("threshold"), axis="x")
    ax.set_xlabel(
        str(spec.get("xlabel") or x_meta.get("label") or x_col),
        fontsize=11,
        fontweight="bold",
    )
    ylabel = (
        "Annotated Compound Count"
        if x_meta.get("show_median")
        else "Compound Count"
    )
    ax.set_ylabel(
        str(spec.get("ylabel") or ylabel), fontsize=11, fontweight="bold"
    )
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    autodetect_integer_axes(ax, x_data=all_vals, is_histogram=True)
    _add_reference_legend_with_headroom(ax, is_zero_bottom=True)
    return

  if ptype == "box":
    positions = np.arange(1, len(series_list) + 1)
    bp = ax.boxplot(
        [s.to_numpy(dtype=float) for s in series_list],
        positions=positions,
        widths=0.48,
        patch_artist=True,
        showfliers=False,
    )
    for patch, color in zip(bp["boxes"], colors_list):
      patch.set_facecolor(color)
      patch.set_alpha(0.35)
      patch.set_edgecolor(PUBLICATION_PALETTE["foreground"])
      patch.set_linewidth(1.1)
    for median_line in bp["medians"]:
      median_line.set_color(PUBLICATION_PALETTE["foreground"])
      median_line.set_linewidth(1.8)

    rng = np.random.default_rng(42)
    for pos, vals, color in zip(positions, series_list, colors_list):
      jitter = rng.uniform(-0.14, 0.14, size=len(vals))
      ax.scatter(
          np.full(len(vals), pos) + jitter,
          vals,
          color=color,
          alpha=0.8,
          edgecolors=PUBLICATION_PALETTE["foreground"],
          linewidths=0.4,
          s=26,
          zorder=3,
      )
    _draw_threshold_line(ax, x_meta.get("threshold"), axis="y")
    ax.set_xticks(positions)
    max_w = (
        18 if len(series_list) <= 2 else (14 if len(series_list) == 3 else 12)
    )
    xtick_labels = [
        _format_box_xtick_label(
            name,
            len(s),
            mean_val=float(s.mean()) if show_mean else None,
            max_width=max_w,
        )
        for name, s in zip(raw_names, series_list)
    ]
    needs_rotation = len(series_list) >= 4 or any(
        len(str(n)) > 20 for n in raw_names
    )
    ax.set_xticklabels(
        xtick_labels,
        fontsize=8.0 if len(series_list) >= 4 else 8.5,
        rotation=15 if needs_rotation else 0,
        ha="right" if needs_rotation else "center",
    )
    ax.set_ylabel(
        str(spec.get("ylabel") or x_meta.get("label") or x_col),
        fontsize=11,
        fontweight="bold",
    )
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    autodetect_integer_axes(ax, y_data=all_vals, is_histogram=False)
    _add_reference_legend_with_headroom(ax, is_zero_bottom=False)
    return

  if ptype == "kde":
    vmin = float(all_vals.min())
    vmax = float(all_vals.max())
    span = max(vmax - vmin, 1e-3)
    grid_min = max(0.0, vmin - 0.15 * span) if vmin >= 0 else vmin - 0.15 * span
    grid_max = vmax + 0.15 * span
    grid = np.linspace(grid_min, grid_max, 200)

    for vals, color in zip(series_list, colors_list):
      arr = vals.to_numpy(dtype=float)
      if len(arr) >= 3 and np.std(arr) > 1e-6:
        kde = scipy_stats.gaussian_kde(arr)
        density = kde(grid)
        ax.plot(grid, density, color=color, linewidth=2.0)
        ax.fill_between(grid, density, color=color, alpha=0.2)
      else:
        ax.hist(
            arr,
            bins=10,
            density=True,
            histtype="step",
            color=color,
            linewidth=1.8,
        )
    if x_meta.get("show_median"):
      ax.axvline(
          all_vals.median(),
          color=PUBLICATION_PALETTE["foreground"],
          linestyle="--",
          linewidth=1.2,
          label=f"Median ({all_vals.median():.2f})",
      )
    _draw_threshold_line(ax, x_meta.get("threshold"), axis="x")
    ax.set_xlabel(
        str(spec.get("xlabel") or x_meta.get("label") or x_col),
        fontsize=11,
        fontweight="bold",
    )
    ax.set_ylabel("Density", fontsize=11, fontweight="bold")
    ax.set_ylim(bottom=0)
    ax.grid(True, linestyle="--", alpha=0.5)
    autodetect_integer_axes(ax, x_data=all_vals, is_histogram=False)
    _add_reference_legend_with_headroom(ax, is_zero_bottom=True)
    return


def plot_chemical_space_distributions(
    df: pd.DataFrame,
    output_png: str,
    title: str = "Congeneric Chemical Space & Property Distributions",
    group_by: Optional[str] = None,
    panels: Optional[Sequence[Union[str, dict[str, object]]]] = None,
    highlight_compounds: Optional[Sequence[str]] = None,
    output_csv: str = "curated_corpus.csv",
) -> None:
  """Renders a configurable multi-panel publication-grade chemical space distribution figure."""
  if not panels:
    raise ValueError(
        _format_missing_panels_help(
            df, output_csv=output_csv, output_plot=output_png
        )
    )
  set_safe_plot_style()
  _, groups = _resolve_grouping(df, group_by=group_by)

  panel_specs = [_parse_panel_spec(p) for p in _expand_panel_specs(panels)]
  if highlight_compounds:
    hl_list = [str(x).strip() for x in highlight_compounds if str(x).strip()]
    for spec in panel_specs:
      if not spec.get("_highlight_compounds"):
        spec["_highlight_compounds"] = hl_list
  n_panels = len(panel_specs)
  if n_panels == 1:
    nrows, ncols, figsize = 1, 1, (8.0, 6.2)
  elif n_panels == 2:
    nrows, ncols, figsize = 1, 2, (13.5, 6.2)
  elif n_panels <= 4:
    nrows, ncols, figsize = 2, 2, (13.5, 10.8)
  else:
    nrows, ncols, figsize = 2, 3, (19.0, 11.2)

  fig, axes = plt.subplots(nrows, ncols, figsize=figsize, dpi=300)
  axes_flat = np.atleast_1d(axes).flatten()

  for idx, spec in enumerate(panel_specs):
    if idx >= len(axes_flat):
      break
    letter = chr(65 + idx)
    _render_panel(axes_flat[idx], df, spec, groups, letter)

  for idx in range(n_panels, len(axes_flat)):
    axes_flat[idx].set_visible(False)

  plt.suptitle(title, fontsize=14, fontweight="bold", y=0.985)

  cohort_handles = [
      Patch(
          facecolor=color,
          edgecolor=PUBLICATION_PALETTE["foreground"],
          linewidth=0.7,
          alpha=0.85,
          label=f"{label} (n={int(mask.sum())})",
      )
      for label, mask, color in groups
      if int(mask.sum()) > 0
  ]
  if cohort_handles:
    ncol_leg = min(len(cohort_handles), 3)
    n_leg_rows = (len(cohort_handles) + ncol_leg - 1) // ncol_leg
    top_rect = 0.925 if n_leg_rows == 1 else 0.895
    fig.legend(
        handles=cohort_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
        ncol=ncol_leg,
        frameon=True,
        framealpha=0.95,
        fontsize=9.5,
        title="Series / Structural Cohorts",
        title_fontsize=9.5,
    )
    plt.tight_layout(rect=[0.0, 0.0, 1.0, top_rect])
  else:
    plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])

  os.makedirs(
      os.path.dirname(os.path.abspath(output_png)) or ".", exist_ok=True
  )
  plt.savefig(output_png)
  plt.close()
  print(
      "Saved chemical space distribution plot to:"
      f" {os.path.abspath(output_png)}",
      flush=True,
  )
  print(
      "\n![Chemical Space & Property"
      f" Distributions]({os.path.abspath(output_png)})\n",
      flush=True,
  )


def _is_patent_fulltext_source(value: object) -> bool:
  """Reports whether a provenance label denotes a patent's own full-text table."""
  text = str(value).lower()
  return "patent" in text and "full" in text and "text" in text


def _combine_provenance(values: Sequence[object]) -> str:
  """Joins provenance labels into a sorted, deduplicated, comma-separated string."""
  seen: set[str] = set()
  for v in values:
    if pd.notna(v) and str(v).strip() != "":
      for part in str(v).split(","):
        p_clean = part.strip()
        if p_clean and p_clean.lower() not in ("nan", "none"):
          seen.add(p_clean)
  return ",".join(sorted(seen))


_POTENCY_CELL_RE = re.compile(
    r"^([<>~=≈≤≥]*)\s*(\d+(?:\.\d+)?)\s*(pm|nm|um|μm|µm|mm)?$",
    re.IGNORECASE,
)
_LOG_CELL_RE: Final[re.Pattern[str]] = re.compile(
    r"^([<>~=≈≤≥]*)\s*(-?\d+(?:\.\d+)?)$"
)


def _canon_potency_qual(raw_qual: str) -> str:
  """Normalizes inequality prefixes into '<', '>', or ''."""
  q = raw_qual.strip()
  if q.startswith(("<", "≤")):
    return "<"
  if q.startswith((">", "≥")):
    return ">"
  return ""


def _default_unit_for_col(col_name: str) -> str:
  """Infers default molar concentration unit from column name suffix."""
  cl = col_name.strip().lower()
  if cl.endswith("_range"):
    cl = cl[:-6]
  if cl.endswith("_pm"):
    return "pm"
  if cl.endswith(("_um", "_μm", "_µm")):
    return "um"
  if cl.endswith("_mm"):
    return "mm"
  return "nm"


def _parse_potency_cell_with_qual(
    val: Any, default_unit: str = "nm"
) -> tuple[float | None, str]:
  """Parses a raw potency cell into (-log10(M), canonical_qualifier)."""
  if val is None or pd.isna(val):
    return (None, "")
  s = str(val).strip()
  if not s or s.lower() in ("nan", "none", "nd", "n/a", "-"):
    return (None, "")
  m = _POTENCY_CELL_RE.match(s)
  if not m:
    return (None, "")
  try:
    num = float(m.group(2))
  except ValueError:
    return (None, "")
  if not np.isfinite(num) or num <= 0:
    return (None, "")
  qual = _canon_potency_qual(m.group(1) or "")
  unit = (m.group(3) or default_unit).lower()
  if unit == "pm":
    return (round(12.0 - float(np.log10(num)), 3), qual)
  if unit in ("um", "μm", "µm"):
    return (round(6.0 - float(np.log10(num)), 3), qual)
  if unit == "mm":
    return (round(3.0 - float(np.log10(num)), 3), qual)
  return (round(9.0 - float(np.log10(num)), 3), qual)


def _parse_potency_cell_to_pactivity(val: object) -> Optional[float]:
  """Parses a raw potency cell (e.g.

  '<0.200', '0.5 nM', '0.003 uM') into -log10(M).
  """
  return _parse_potency_cell_with_qual(val, default_unit="nm")[0]


def _is_nonexact_raw_potency_cell(val: object) -> bool:
  """Returns True if `val` carries censorship, approx, range, or bin text."""
  if val is None or pd.isna(val):
    return False
  s = str(val).strip()
  if not s or s.lower() in ("nan", "none", "nd", "n/a", "-"):
    return False
  m = _POTENCY_CELL_RE.match(s)
  if not m:
    return True
  raw_q = (m.group(1) or "").strip()
  return bool(raw_q and raw_q != "=")


_POTENCY_METRIC: Final[str] = r"(?:(?:ic|ec|ac|gi|ed|cc)_?50|ki|kd)"
_ASSAY_SUFFIX_RE: Final[re.Pattern[str]] = re.compile(
    r"_(?:nM|pActivity|range)$", re.IGNORECASE
)
_ANON_PATENT_COL_RE: Final[re.Pattern[str]] = re.compile(
    r"^(?:patent_table|pat|[A-Za-z]{2,4}\d{5,}[A-Za-z0-9]*)_col_(\d+)"
    r"(?:_(?:nM|pActivity|range))?$",
    re.IGNORECASE,
)
_GENERIC_PATENT_COL_RE: Final[re.Pattern[str]] = re.compile(
    r"^patent_table_col_(\d+)(?:_nM)?$"
)
_NON_ALNUM_RE: Final[re.Pattern[str]] = re.compile(r"[^A-Za-z0-9]+")
_NAMED_POTENCY_COL_RE: Final[re.Pattern[str]] = re.compile(
    rf"(?:^|_)p?(?:{_POTENCY_METRIC}|mic)(?:_|$)", re.IGNORECASE
)
_LOG_POTENCY_COL_RE: Final[re.Pattern[str]] = re.compile(
    rf"(?:_pactivity$|(?:^|_)p{_POTENCY_METRIC}(?:_|$))", re.IGNORECASE
)
_BIO_TYPE_TOKEN_RE: Final[re.Pattern[str]] = re.compile(
    r"(?:wild|cell)_type|(?:sub|geno|pheno)type", re.IGNORECASE
)
_META_UNIT_SUFFIX_RE: Final[re.Pattern[str]] = re.compile(
    r"(?:^|_)(?:units?|type|qualifier|operator|relation|comment)$",
    re.IGNORECASE,
)
_POTENCY_STEM_NORM_RE: Final[re.Pattern[str]] = re.compile(
    rf"(^|_)p({_POTENCY_METRIC})(?=_|$)", re.IGNORECASE
)
_PACTIVITY_MATCH_TOL: Final[float] = 0.02


def _cell_potency_sig(
    val: object, col_name: str = ""
) -> tuple[float | None, str, str] | None:
  """Returns a comparable (pActivity, qual, lower_text) potency signature."""
  if val is None or pd.isna(val):
    return None
  s = str(val).strip()
  if not s or s.lower() in ("nan", "none", "nd", "n/a", "-"):
    return None
  norm_s = re.sub(r"\s+", " ", s.lower())
  cl = col_name.lower()
  if _LOG_POTENCY_COL_RE.search(cl):
    m_log = _LOG_CELL_RE.match(s)
    if m_log:
      f = float(m_log.group(2))
      if np.isfinite(f):
        return (round(f, 3), _canon_potency_qual(m_log.group(1) or ""), norm_s)
    return (None, "", norm_s)
  pa, qual = _parse_potency_cell_with_qual(
      s, default_unit=_default_unit_for_col(cl)
  )
  if pa is not None:
    return (pa, qual, norm_s)
  return (None, "", norm_s)


_DOC_INVALID_STRINGS = frozenset(
    {"", "nan", "none", "null", "na", "n/a", "true", "false"}
)
_DOC_WIPO_RE = re.compile(
    r"\bWO[-\s/]?((?:19|20)\d{2}|\d{2})[-\s/]?\d{4,7}(?:[-\s/]?[A-Z]\d?)?\b",
    re.IGNORECASE,
)
_DOC_US_PUB_RE = re.compile(
    r"\bUS[-\s/]?(20\d{2})[-\s/]?\d{6,8}(?:[-\s/]?A\d)?\b(?![-\s/]?B\d)",
    re.IGNORECASE,
)


def _parse_valid_calendar_year(val: object) -> Optional[int]:
  """Extracts a valid calendar year in [1800, 2100], rejecting booleans."""
  if val is None or isinstance(val, (bool, np.bool_)):
    return None
  try:
    if pd.isna(val):
      return None
  except (TypeError, ValueError):
    pass
  if isinstance(val, (int, float, np.integer, np.floating)):
    fval = float(val)
    if not np.isfinite(fval):
      return None
    ival = int(round(fval))
    if abs(fval - ival) < 1e-6:
      if 1800 <= ival <= 2100:
        return ival
      if 19000101 <= ival <= 20991231:
        return int(str(ival)[:4])
    return None
  text = str(val).strip()
  if text.lower() in _DOC_INVALID_STRINGS:
    return None
  m_iso = re.match(
      r"^(18\d{2}|19\d{2}|20\d{2})(?:[-/]\d{1,2}(?:[-/]\d{1,2})?)?(?:[T\s].*)?$",
      text,
  )
  if m_iso:
    return int(m_iso.group(1))
  if re.fullmatch(
      r"(19\d{2}|20\d{2})(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])", text
  ):
    return int(text[:4])
  try:
    fval = float(text)
    if np.isfinite(fval) and abs(fval - round(fval)) < 1e-6:
      ival = int(round(fval))
      if 1800 <= ival <= 2100:
        return ival
  except ValueError:
    pass
  return None


def _parse_valid_iso_date(val: object) -> Optional[str]:
  """Extracts a normalized YYYY-MM-DD date string, rejecting booleans."""
  if val is None or isinstance(val, (bool, np.bool_)):
    return None
  try:
    if pd.isna(val):
      return None
  except (TypeError, ValueError):
    pass
  text = str(val).strip()
  if text.lower() in _DOC_INVALID_STRINGS:
    return None
  m_dt = re.match(
      r"^(18\d{2}|19\d{2}|20\d{2})[-/](\d{1,2})[-/](\d{1,2})(?:[T\s].*)?$",
      text,
  )
  if m_dt:
    return f"{int(m_dt.group(1)):04d}-{int(m_dt.group(2)):02d}-{int(m_dt.group(3)):02d}"
  return None


def _extract_unambiguous_patent_year(val: object) -> Optional[int]:
  """Extracts earliest publication year from WIPO PCT or US pre-grant IDs."""
  if val is None or isinstance(val, (bool, np.bool_)):
    return None
  try:
    if pd.isna(val):
      return None
  except (TypeError, ValueError):
    pass
  text = str(val).strip()
  if not text or text.lower() in _DOC_INVALID_STRINGS:
    return None
  years: list[int] = []
  for match in _DOC_WIPO_RE.finditer(text):
    raw_yr = match.group(1)
    if len(raw_yr) == 4:
      yr = int(raw_yr)
    else:
      yr2 = int(raw_yr)
      yr = (1900 + yr2) if yr2 >= 78 else (2000 + yr2)
    if 1978 <= yr <= 2100:
      years.append(yr)
  for match in _DOC_US_PUB_RE.finditer(text):
    yr = int(match.group(1))
    if 2001 <= yr <= 2100:
      years.append(yr)
  return min(years) if years else None


def backfill_document_dates(df: pd.DataFrame) -> pd.DataFrame:
  """Sanitizes boolean placeholders and backfills true chronological dates."""
  if df.empty:
    return df
  work = df.copy()
  year_like_cols = (
      "document_year",
      "publication_year",
      "chembl_year",
      "year",
      "Year",
  )
  date_like_cols = (
      "publication_date",
      "priority_date",
      "filing_date",
      "date",
      "Date",
  )
  has_any_date_col = any(
      c in work.columns for c in (*year_like_cols, *date_like_cols)
  )

  # Extract any calendar years embedded in date_like_cols before strict ISO
  # sanitization so year-only tokens still backfill `document_year` without
  # leaving non-ISO strings in `publication_date`.
  raw_date_col_years: dict[str, list[Optional[int]]] = {}
  for dcol in date_like_cols:
    if dcol in work.columns:
      raw_vals = work[dcol].tolist()
      raw_date_col_years[dcol] = [
          _parse_valid_calendar_year(v) for v in raw_vals
      ]
      cleaned_d = [
          iso if (iso := _parse_valid_iso_date(v)) is not None else np.nan
          for v in raw_vals
      ]
      work[dcol] = pd.Series(cleaned_d, index=work.index, dtype=object)

  # Sanitize boolean False/True or invalid calendar sentinels in year columns
  for ycol in year_like_cols:
    if ycol in work.columns:
      cleaned_y = [_parse_valid_calendar_year(v) for v in work[ycol].tolist()]
      work[ycol] = pd.Series(cleaned_y, index=work.index, dtype="Int64")

  years: list[int | None] = (
      [int(v) if pd.notna(v) else None for v in work["document_year"].tolist()]
      if "document_year" in work.columns
      else [None] * len(work)
  )

  # Backfill/clamp across public disclosure date/year columns only (never clamp
  # public publication years down to non-public filing_date / priority_date).
  for sec_col in (
      "publication_date",
      "publication_year",
      "chembl_year",
      "year",
      "Year",
      "date",
      "Date",
  ):
    if sec_col not in work.columns:
      continue
    sec_yrs = (
        raw_date_col_years[sec_col]
        if sec_col in raw_date_col_years
        else [_parse_valid_calendar_year(v) for v in work[sec_col].tolist()]
    )
    for i, cur in enumerate(years):
      cand_yr = sec_yrs[i]
      if cand_yr is not None and (cur is None or cand_yr < cur):
        years[i] = cand_yr

  def _propagate_years(
      yr_list: list[int | None],
      id_cols: tuple[str, ...],
      overwrite_with_earlier: bool = False,
  ) -> list[int | None]:
    out_yr = list(yr_list)
    for id_col in id_cols:
      if id_col not in work.columns:
        continue
      if not overwrite_with_earlier and all(y is not None for y in out_yr):
        continue
      tok_to_min: dict[str, int] = {}
      col_vals = work[id_col].tolist()
      for i, raw_id in enumerate(col_vals):
        yv = out_yr[i]
        if yv is None or raw_id is None or isinstance(raw_id, (bool, np.bool_)):
          continue
        try:
          if pd.isna(raw_id):
            continue
        except (TypeError, ValueError):
          pass
        for tok in re.split(r"[,;|]+", str(raw_id)):
          tok = tok.strip()
          if not tok or tok.lower() in _DOC_INVALID_STRINGS:
            continue
          if tok not in tok_to_min or yv < tok_to_min[tok]:
            tok_to_min[tok] = yv
      if not tok_to_min:
        continue
      for i, raw_id in enumerate(col_vals):
        cur_yv = out_yr[i]
        if (
            (not overwrite_with_earlier and cur_yv is not None)
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
            tok_to_min[tok.strip()]
            for tok in re.split(r"[,;|]+", str(raw_id))
            if tok.strip() in tok_to_min
        ]
        if matched:
          earliest = min(matched)
          if cur_yv is None or earliest < cur_yv:
            out_yr[i] = earliest
    return out_yr

  doc_id_cols = ("document_chembl_id", "document_id", "patent_id")
  cmpd_id_cols = ("chembl_id", "inchikey")

  if any(y is not None for y in years) and any(y is None for y in years):
    years = _propagate_years(years, doc_id_cols)

  # Fallback to unambiguous WIPO/US-A1 patent IDs when explicit public dates are absent
  if any(y is None for y in years):
    pat_cols = [
        c for c in ("patent_id", "document_id", "source") if c in work.columns
    ]
    inferred = False
    for i, cur in enumerate(years):
      if cur is not None:
        continue
      cands = [
          yr
          for c in pat_cols
          if (yr := _extract_unambiguous_patent_year(work[c].iloc[i]))
          is not None
      ]
      if cands:
        years[i] = min(cands)
        inferred = True
    if inferred and any(y is None for y in years):
      years = _propagate_years(years, doc_id_cols)

  # Non-public fallback (filing_date before priority_date) only when public dates are absent
  if any(y is None for y in years):
    filled_fallback = False
    for fb_col in ("filing_date", "priority_date"):
      if fb_col not in work.columns:
        continue
      fb_yrs = raw_date_col_years.get(fb_col, [])
      for i, cur in enumerate(years):
        if cur is None and i < len(fb_yrs) and fb_yrs[i] is not None:
          years[i] = fb_yrs[i]
          filled_fallback = True
    if filled_fallback and any(y is None for y in years):
      years = _propagate_years(years, doc_id_cols)

  pre_cmpd_years = list(years)
  if any(y is not None for y in years):
    years = _propagate_years(years, cmpd_id_cols, overwrite_with_earlier=True)

  if has_any_date_col or any(y is not None for y in years):
    work["document_year"] = pd.Series(years, index=work.index, dtype="Int64")
    for dcol in date_like_cols:
      if dcol in work.columns:
        reconciled = []
        for i, dval in enumerate(work[dcol].tolist()):
          iso = _parse_valid_iso_date(dval)
          yr_i = years[i]
          pre_yr = pre_cmpd_years[i]
          if (
              iso is not None
              and yr_i is not None
              and pre_yr is not None
              and yr_i < pre_yr
              and int(iso[:4]) > yr_i
          ):
            reconciled.append(np.nan)
          else:
            reconciled.append(iso if iso is not None else np.nan)
        work[dcol] = pd.Series(reconciled, index=work.index, dtype=object)
  return work


_backfill_document_dates = backfill_document_dates


def deduplicate_by_smiles(df: pd.DataFrame, smiles_col: str) -> pd.DataFrame:
  """Deduplicates rows on the SMILES string while preserving provenance.

  Collisions are resolved with an explicit, stable ranking: a row carrying a
  bioactivity wins over one that does not, and a row sourced from a patent's
  own full-text table wins over any other source at equal activity status.
  Before any row is dropped, the provenance of every row sharing a structure is
  merged into companion columns ('all_sources', and 'all_dataset_roles' when
  'dataset_role' exists) and 'patent_example_id' / 'patent_id' are combined
  so no contributing source or patent example identifier is lost.

  Rows whose SMILES is missing, empty, or whitespace-only carry no structure,
  so they are used to backfill matching patent_example_id / compound_id rows
  before being dropped.

  Args:
    df: Concatenated input table.
    smiles_col: Name of the column holding the SMILES strings.

  Returns:
    A new DataFrame, in input order, holding one row per distinct non-empty
    SMILES string.
  """
  work = _backfill_document_dates(df)

  def _row_source_priority(idx_pos: int) -> int:
    if "source" not in work.columns:
      return 1
    src_val = str(work["source"].iloc[idx_pos] or "").strip().lower()
    if any(
        k in src_val
        for k in ("chembl", "literature", "paper", "journal", "jmc")
    ):
      return 0
    if _is_patent_fulltext_source(src_val):
      return 1
    return 2

  def _canon_and_flat_smi(val: object) -> Tuple[str, str, bool, str]:
    if pd.isna(val):
      return "", "", False, ""
    s = str(val).strip()
    if not s or s.lower() in ("nan", "none"):
      return "", "", False, ""
    m = Chem.MolFromSmiles(s)
    if m is None:
      return s, s, ("@" in s or "/" in s or "\\" in s), ""
    iso_s = Chem.MolToSmiles(m, isomericSmiles=True)
    flat_s = Chem.MolToSmiles(m, isomericSmiles=False)
    try:
      ik_s = Chem.MolToInchiKey(m) or ""
    except Exception:
      ik_s = ""
    return iso_s, flat_s, (iso_s != flat_s), ik_s

  parsed_smi = [_canon_and_flat_smi(v) for v in work[smiles_col]]
  raw_iso_keys = [t[0] for t in parsed_smi]
  raw_flat_keys = [t[1] for t in parsed_smi]
  raw_has_stereo = [t[2] for t in parsed_smi]
  raw_iks = [t[3] for t in parsed_smi]

  # Unify tautomer and protonation-state SMILES that share the identical
  # Standard InChIKey connectivity + stereochemistry blocks (blocks 0 and 1 of
  # 'CONNECTIVITY-STEREO-PROTONATION', e.g. annular 1H/2H tautomers or
  # protonated cation '-O' vs neutral parent '-N'), and unify _flat_smiles_key
  # across identical connectivity blocks (block 0).
  ik_to_best_iso: dict[str, Tuple[Tuple[int, int], str]] = {}
  conn_to_best_flat: dict[str, Tuple[Tuple[int, int], str]] = {}
  for pos_i, (iso_k, flat_k, ik_k) in enumerate(
      zip(raw_iso_keys, raw_flat_keys, raw_iks)
  ):
    if not ik_k:
      continue
    ik_blocks = ik_k.split("-")
    conn_block = ik_blocks[0]
    stereo_block = "-".join(ik_blocks[:2])
    is_neutral = len(ik_blocks) == 3 and ik_blocks[2] == "N"
    prio = (0 if is_neutral else 1, _row_source_priority(pos_i))
    if (
        stereo_block not in ik_to_best_iso
        or prio < ik_to_best_iso[stereo_block][0]
    ):
      ik_to_best_iso[stereo_block] = (prio, iso_k)
    if (
        conn_block not in conn_to_best_flat
        or prio < conn_to_best_flat[conn_block][0]
    ):
      conn_to_best_flat[conn_block] = (prio, flat_k)

  work["_smiles_key"] = [
      ik_to_best_iso["-".join(ik_k.split("-")[:2])][1]
      if (ik_k and "-".join(ik_k.split("-")[:2]) in ik_to_best_iso)
      else iso_k
      for iso_k, ik_k in zip(raw_iso_keys, raw_iks)
  ]
  work["_flat_smiles_key"] = [
      conn_to_best_flat[ik_k.split("-")[0]][1]
      if (ik_k and ik_k.split("-")[0] in conn_to_best_flat)
      else flat_k
      for flat_k, ik_k in zip(raw_flat_keys, raw_iks)
  ]
  work["_has_defined_stereo"] = raw_has_stereo
  work["_has_structure"] = work["_smiles_key"] != ""
  work["_is_neutral_ik"] = [
      bool(ik_k and ik_k.endswith("-N")) for ik_k in raw_iks
  ]

  # Qualify patent_example_id with patent_id (e.g. '<PATENT_ID>:Example_N')
  # whenever patent_id is available and patent_example_id is not yet qualified.
  if "patent_example_id" in work.columns and "patent_id" in work.columns:
    qualified_ex = []
    for pat_v, ex_v in zip(work["patent_id"], work["patent_example_id"]):
      ex_s = str(ex_v).strip() if pd.notna(ex_v) else ""
      pat_s = str(pat_v).strip() if pd.notna(pat_v) else ""
      if (
          ex_s
          and ex_s.lower() not in ("nan", "none")
          and ":" not in ex_s
          and pat_s
          and pat_s.lower() not in ("nan", "none")
          and "," not in pat_s
      ):
        qualified_ex.append(f"{pat_s}:{ex_s}")
      else:
        qualified_ex.append(ex_v)
    work["patent_example_id"] = qualified_ex

  if "chiral_peak" in work.columns:
    work["chiral_peak"] = work["chiral_peak"].astype(object)
  if "compound_id" in work.columns:
    if "chembl_id" not in work.columns:
      work["chembl_id"] = pd.Series(index=work.index, dtype=object)
    else:
      work["chembl_id"] = work["chembl_id"].astype(object)
    for idx_c, row_c in work.iterrows():
      cid_str = str(row_c.get("compound_id") or "").strip()
      cur_chembl = row_c.get("chembl_id")
      if re.match(r"^CHEMBL\d+$", cid_str, flags=re.IGNORECASE) and (
          pd.isna(cur_chembl)
          or str(cur_chembl).strip().lower() in ("", "nan", "none")
      ):
        work.at[idx_c, "chembl_id"] = cid_str.upper()

  # Backfill missing pchembl_value on patent/assay rows before ranking or peak-matching
  if "pchembl_value" not in work.columns:
    work["pchembl_value"] = np.nan
  else:
    work["pchembl_value"] = pd.to_numeric(
        work["pchembl_value"], errors="coerce"
    )

  pat_cols = [
      c
      for c in work.columns
      if str(c).startswith(("patent_table_col_", "patent_assay_"))
      and not str(c).lower().endswith("_range")
  ]
  for idx_w, row_w in work[work["pchembl_value"].isna()].iterrows():
    std_v = row_w.get("standard_value")
    std_u = str(row_w.get("standard_units") or "nM").strip()
    candidates = [(row_w.get(pc), str(pc)) for pc in pat_cols]
    if pd.notna(std_v):
      candidates.insert(0, (f"{std_v} {std_u}", ""))
    for raw_v, col_n in candidates:
      sig = _cell_potency_sig(raw_v, col_n)
      if sig is not None and sig[0] is not None:
        work.at[idx_w, "pchembl_value"] = sig[0]
        break

  # Reconcile chiral_peak ("Peak 1" / "Peak 2") and unassigned stereocenters
  # within each flat 2D constitution (_flat_smiles_key):
  # 1. If a flat constitution has stereospecific @/@@ isomers in the corpus,
  #    map Peak 1 / Peak 2 groups (ranked by potency) onto the corresponding
  #    stereospecific isomers (ranked by potency), and map unannotated flat (NaN)
  #    rows onto the primary stereoisomer when matching in potency or unassayed
  #    (unless the unannotated row is a distinct racemic patent Example alongside
  #    Peak 1 / Peak 2 patent Examples).
  # 2. If a flat constitution has both Peak 1 and Peak 2 but fewer than 2
  #    stereospecific isomers in the corpus, generate the two tetrahedral
  #    enantiomers (@ and @@) when an unassigned chiral center is present so each
  #    peak resolves to a distinct stereoisomeric SMILES; otherwise keep
  #    Peak 1 and Peak 2 as distinct _dedup_key entries while merging NaN rows
  #    into the primary peak sharing that _smiles_key.
  cp_series = (
      work["chiral_peak"]
      .astype(str)
      .where(work["chiral_peak"].notna(), "")
      .str.strip()
      if "chiral_peak" in work.columns
      else pd.Series("", index=work.index)
  )
  cp_series = cp_series.where(~cp_series.str.lower().isin(["nan", "none"]), "")
  work["_clean_peak"] = cp_series
  work["_dedup_key"] = work["_smiles_key"]

  def _two_enantiomer_smiles(
      flat_smi: str, known_iso: Optional[str] = None
  ) -> Optional[Tuple[str, str]]:
    base_m = Chem.MolFromSmiles(known_iso if known_iso else flat_smi)
    if base_m is None:
      return None
    centers = Chem.FindMolChiralCenters(base_m, includeUnassigned=True)
    if not centers:
      return None
    # Prefer inverting a defined stereocenter on known_iso, or assigning the first unassigned center
    defined_c = [idx for idx, tag in centers if tag != "?"]
    unassigned_c = [idx for idx, tag in centers if tag == "?"]
    c_idx = defined_c[0] if defined_c else unassigned_c[0]
    m1 = Chem.Mol(base_m)
    m2 = Chem.Mol(base_m)
    m1.GetAtomWithIdx(c_idx).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    m2.GetAtomWithIdx(c_idx).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
    s1 = Chem.MolToSmiles(m1, isomericSmiles=True)
    s2 = Chem.MolToSmiles(m2, isomericSmiles=True)
    if s1 == s2:
      return None
    if known_iso == s2:
      return s2, s1
    return s1, s2

  for flat_k, grp in work[work["_has_structure"]].groupby(
      "_flat_smiles_key", sort=False
  ):
    if len(grp) <= 1:
      continue
    def_rows = grp[grp["_has_defined_stereo"]]
    undef_rows = grp[~grp["_has_defined_stereo"]]

    # Rank distinct defined stereoisomers by max potency descending
    stereo_ranks: list[Tuple[float, str]] = []
    for iso_s, iso_grp in def_rows.groupby("_smiles_key", sort=False):
      max_p = iso_grp["pchembl_value"].max()
      stereo_ranks.append(
          (float(max_p) if pd.notna(max_p) else -1.0, str(iso_s))
      )
    stereo_ranks.sort(key=lambda x: x[0], reverse=True)
    ordered_isos = [x[1] for x in stereo_ranks]

    # Distinct non-empty chiral peaks in this flat constitution
    non_empty_peaks = [p for p in grp["_clean_peak"].unique() if p]
    if len(non_empty_peaks) >= 2:
      # Rank peaks by max potency descending
      peak_ranks: list[Tuple[float, str]] = []
      for pk in non_empty_peaks:
        pk_max = grp[grp["_clean_peak"] == pk]["pchembl_value"].max()
        peak_ranks.append(
            (float(pk_max) if pd.notna(pk_max) else -1.0, str(pk))
        )
      peak_ranks.sort(key=lambda x: x[0], reverse=True)
      ordered_peaks = [x[1] for x in peak_ranks]

      if len(ordered_isos) < 2:
        pair = _two_enantiomer_smiles(
            str(flat_k), ordered_isos[0] if ordered_isos else None
        )
        if pair is not None:
          ordered_isos = [pair[0], pair[1]]

      peak_rows_have_pat_ex = False
      if "patent_example_id" in grp.columns:
        pk_sub = grp[grp["_clean_peak"] != ""]
        peak_rows_have_pat_ex = any(
            pd.notna(v) and str(v).strip().lower() not in ("", "nan", "none")
            for v in pk_sub["patent_example_id"]
        )

      if len(ordered_isos) >= 2:
        peak_to_iso = {
            ordered_peaks[i]: ordered_isos[min(i, len(ordered_isos) - 1)]
            for i in range(len(ordered_peaks))
        }
        for idx_g, row_g in grp.iterrows():
          pk = row_g["_clean_peak"]
          if pk in peak_to_iso:
            target_iso = peak_to_iso[pk]
            work.at[idx_g, "_smiles_key"] = target_iso
            work.at[idx_g, "_dedup_key"] = target_iso
          elif not row_g["_has_defined_stereo"]:
            row_ex = (
                str(row_g.get("patent_example_id") or "").strip()
                if pd.notna(row_g.get("patent_example_id"))
                else ""
            )
            if (
                peak_rows_have_pat_ex
                and row_ex
                and row_ex.lower() not in ("nan", "none")
            ):
              # Distinct racemic patent Example alongside Peak 1 / Peak 2
              # patent Examples: preserve its flat racemic SMILES and do not
              # merge into Peak 1 or Peak 2.
              work.at[idx_g, "_dedup_key"] = f"{row_g['_smiles_key']}::racemate"
              continue
            # Unannotated flat row (chiral_peak=NaN): merge into closest/primary stereoisomer
            row_p = row_g["pchembl_value"]
            if pd.notna(row_p) and len(stereo_ranks) >= 2:
              best_iso = min(
                  stereo_ranks, key=lambda sr: abs(sr[0] - float(row_p))
              )[1]
            else:
              best_iso = ordered_isos[0]
            work.at[idx_g, "_smiles_key"] = best_iso
            work.at[idx_g, "_dedup_key"] = best_iso
      else:
        # True 2D-identical atropisomers without tetrahedral stereocenters:
        # keep Peak 1 and Peak 2 separate, and merge any NaN row into the primary peak.
        primary_pk = ordered_peaks[0]
        for idx_g, row_g in grp.iterrows():
          pk = row_g["_clean_peak"] or primary_pk
          work.at[idx_g, "_dedup_key"] = f"{row_g['_smiles_key']}::peak={pk}"
    elif not undef_rows.empty:

      def _extract_id_tokens(sub_df: pd.DataFrame, col_name: str) -> set[str]:
        tokens: set[str] = set()
        if col_name not in sub_df.columns:
          return tokens
        for val_tok in sub_df[col_name]:
          if pd.notna(val_tok):
            for part_tok in str(val_tok).split(","):
              clean_tok = part_tok.strip()
              if clean_tok and clean_tok.lower() not in ("nan", "none"):
                tokens.add(clean_tok)
        return tokens

      def_pat_ex_set = _extract_id_tokens(def_rows, "patent_example_id")
      def_paper_id_set = _extract_id_tokens(def_rows, "paper_compound_id")

      has_distinct_racemic_example = False
      for _, row_u in undef_rows.iterrows():
        if str(row_u.get("_clean_peak") or "").strip():
          continue
        u_pat_ex = _extract_id_tokens(
            pd.DataFrame([row_u]), "patent_example_id"
        )
        if def_pat_ex_set and u_pat_ex and u_pat_ex.isdisjoint(def_pat_ex_set):
          has_distinct_racemic_example = True
          break
        u_paper_id = _extract_id_tokens(
            pd.DataFrame([row_u]), "paper_compound_id"
        )
        if (
            def_paper_id_set
            and u_paper_id
            and u_paper_id.isdisjoint(def_paper_id_set)
        ):
          has_distinct_racemic_example = True
          break

      if has_distinct_racemic_example:
        # Distinct racemic patent Example or paper compound alongside an
        # explicitly stereospecific isomer: preserve its flat racemic SMILES
        # and do not collapse into or invert the stereospecific isomer.
        for idx_u, row_u in undef_rows.iterrows():
          if not str(row_u.get("_clean_peak") or "").strip():
            work.at[idx_u, "_dedup_key"] = f"{row_u['_smiles_key']}::racemate"
          elif ordered_isos:
            work.at[idx_u, "_smiles_key"] = ordered_isos[0]
            work.at[idx_u, "_dedup_key"] = ordered_isos[0]
      else:
        p_vals_grp = grp["pchembl_value"].dropna()
        p_span = (
            float(p_vals_grp.max() - p_vals_grp.min())
            if len(p_vals_grp) >= 2
            else 0.0
        )
        if len(ordered_isos) < 2 and p_span >= 0.5:
          pair = _two_enantiomer_smiles(
              str(flat_k), ordered_isos[0] if ordered_isos else None
          )
          if pair is not None:
            ordered_isos = [pair[0], pair[1]]
            stereo_ranks = [
                (float(p_vals_grp.max()), pair[0]),
                (float(p_vals_grp.min()), pair[1]),
            ]
        if len(ordered_isos) >= 1:
          # Upgrade unannotated/single-peak flat rows to the matching defined/inferred stereoisomer
          for idx_u, row_u in undef_rows.iterrows():
            row_p = row_u["pchembl_value"]
            if pd.notna(row_p) and len(stereo_ranks) >= 2:
              best_idx = min(
                  range(len(stereo_ranks)),
                  key=lambda i: abs(stereo_ranks[i][0] - float(row_p)),
              )
              best_iso = stereo_ranks[best_idx][1]
              if (
                  "chiral_peak" in work.columns
                  and not str(work.at[idx_u, "_clean_peak"]).strip()
              ):
                work.at[idx_u, "chiral_peak"] = (
                    "Enantiomer 1 (High Potency)"
                    if best_idx == 0
                    else "Enantiomer 2 (Low Potency)"
                )
            else:
              best_iso = ordered_isos[0]
            work.at[idx_u, "_smiles_key"] = best_iso
            work.at[idx_u, "_dedup_key"] = best_iso

  work.loc[work["_has_structure"], smiles_col] = work.loc[
      work["_has_structure"], "_smiles_key"
  ]
  if "canonical_smiles" in work.columns:
    work.loc[work["_has_structure"], "canonical_smiles"] = work.loc[
        work["_has_structure"], "_smiles_key"
    ]
  if "inchikey" in work.columns:

    def _smi_to_ik(s_val: object) -> Optional[str]:
      if pd.isna(s_val) or not str(s_val).strip():
        return None
      m_ik = Chem.MolFromSmiles(str(s_val).strip())
      if m_ik is None:
        return None
      try:
        return Chem.MolToInchiKey(m_ik)
      except Exception:
        return None

    work.loc[work["_has_structure"], "inchikey"] = [
        _smi_to_ik(s) for s in work.loc[work["_has_structure"], "_smiles_key"]
    ]

  for col, merged_col in (
      ("source", "all_sources"),
      ("dataset_role", "all_dataset_roles"),
  ):
    if col not in work.columns:
      continue
    merged = work.groupby("_dedup_key")[col].transform(_combine_provenance)
    # Structureless rows share no structure, so they only describe themselves.
    work[merged_col] = merged.where(
        work["_has_structure"],
        work[col].apply(lambda v: _combine_provenance([v])),
    )

  def _compute_source_rank(v: object) -> int:
    s_low = str(v or "").strip().lower()
    if any(
        k in s_low for k in ("chembl", "literature", "paper", "journal", "jmc")
    ):
      return 0
    if _is_patent_fulltext_source(v):
      return 1
    return 2

  has_act_mask = work["pchembl_value"].notna()
  for disc_act_col in ["activity_bin", "activity_value", "potency_tier"] + [
      c
      for c in work.columns
      if (
          (str(c).startswith("patent_assay_") and not str(c).endswith("_range"))
          or re.match(
              r"^(?:p?(?:Ki|IC50|Kd|EC50)|pActivity)(?:_|$)", str(c), re.I
          )
      )
  ]:
    if disc_act_col in work.columns:
      has_act_mask = has_act_mask | (
          work[disc_act_col].notna()
          & work[disc_act_col].astype(str).str.strip().str.lower().ne("")
          & work[disc_act_col].astype(str).str.strip().str.lower().ne("nan")
          & work[disc_act_col].astype(str).str.strip().str.lower().ne("none")
      )
  work["_has_act"] = has_act_mask
  work["_source_rank"] = (
      work["source"].apply(_compute_source_rank)
      if "source" in work.columns
      else 1
  )

  # Compute per-row source document size before provenance columns are merged
  # across duplicate SMILES so larger primary discovery documents break ties.
  doc_series = pd.Series("", index=work.index, dtype=str)
  for c in ("document_chembl_id", "patent_id", "source"):
    if c in work.columns:
      s = work[c].fillna("").astype(str).str.strip()
      doc_series = doc_series.where(doc_series != "", s)
  doc_counts = doc_series[doc_series != ""].value_counts().to_dict()
  work["_doc_size"] = doc_series.map(lambda k: doc_counts.get(k, 0) if k else 0)

  # Combine patent_example_id, patent_id, and chiral_peak across identical _dedup_key rows
  # so duplicate SMILES from multiple examples retain full provenance when
  # collapsed into a single row.
  for id_col_merge in ("patent_example_id", "patent_id", "chiral_peak"):
    if id_col_merge in work.columns:
      merged_ids = work.groupby("_dedup_key")[id_col_merge].transform(
          _combine_provenance
      )
      work[id_col_merge] = merged_ids.where(
          work["_has_structure"], work[id_col_merge]
      )

  # Pivot any raw (standard_type, target_name, assay_description) activity rows
  # into wide endpoint columns before collapsing duplicate SMILES.
  if "standard_type" in work.columns and "standard_value" in work.columns:
    t_col = (
        work["target_name"]
        if "target_name" in work.columns
        else pd.Series(None, index=work.index)
    )
    d_col = (
        work["assay_description"]
        if "assay_description" in work.columns
        else pd.Series(None, index=work.index)
    )
    for idx, (st_v, sv_v, tg_v, ad_v) in enumerate(
        zip(work["standard_type"], work["standard_value"], t_col, d_col)
    ):
      if pd.notna(st_v) and pd.notna(sv_v):
        ckey = _build_endpoint_column_key(st_v, tg_v, ad_v)
        if ckey:
          if ckey not in work.columns:
            work[ckey] = pd.Series(index=work.index, dtype=object)
          row_label = work.index[idx]
          if pd.isna(work.at[row_label, ckey]):
            try:
              sv_clean: object = float(sv_v)
            except (ValueError, TypeError):
              sv_clean = str(sv_v).strip()
            if work[ckey].dtype != object and isinstance(sv_clean, str):
              work[ckey] = work[ckey].astype(object)
            work.at[row_label, ckey] = sv_clean

  work["_orig_idx"] = work.index
  ordered = work.sort_values(
      by=["_has_act", "_source_rank", "_is_neutral_ik", "_doc_size"],
      ascending=[False, True, False, False],
      kind="stable",
  )
  struct_part = ordered[ordered["_has_structure"]].copy()
  for str_col in (
      "paper_compound_id",
      "compound_name",
      "patent_id",
      "patent_example_id",
      "document_chembl_id",
  ):
    if str_col in struct_part.columns:
      struct_part[str_col] = struct_part[str_col].replace("", np.nan)

  if not struct_part.empty:
    if "all_endpoints" in struct_part.columns:

      def _merge_endpoints_json(vals: Sequence[object]) -> Optional[str]:
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

      struct_part["all_endpoints"] = struct_part.groupby(
          "_dedup_key", sort=False
      )["all_endpoints"].transform(_merge_endpoints_json)

    def _pick_best_human_id(vals: Sequence[object]) -> object:
      non_null = [
          v
          for v in vals
          if pd.notna(v) and str(v).strip().lower() not in ("", "nan", "none")
      ]
      if not non_null:
        return np.nan
      first_val = non_null[0]
      first_str = str(first_val).strip()
      if re.match(
          r"^(?:CHEMBL\d+|CID_\d+|(?:[A-Z]{2}\d+[A-Z0-9]*:)?(?:Example|Compound)_[A-Za-z0-9_.-]+)$",
          first_str,
          re.I,
      ):
        for v in non_null[1:]:
          s = str(v).strip()
          if re.match(r"^[A-Za-z][A-Za-z0-9 -]{3,27}$", s) and not re.match(
              r"^(?:CHEMBL\d+|CID_|US\d+|WO\d+|EP\d+|CN\d+|JP\d+|Example|Compound)",
              s,
              re.I,
          ):
            return v
      return first_val

    for id_c in ("compound_id", "compound_name"):
      if id_c in struct_part.columns:
        struct_part[id_c] = struct_part.groupby("_dedup_key", sort=False)[
            id_c
        ].transform(_pick_best_human_id)

    def _min_valid_year(vals: Sequence[object]) -> object:
      valid_yrs = [
          yr for v in vals if (yr := _parse_valid_calendar_year(v)) is not None
      ]
      return min(valid_yrs) if valid_yrs else pd.NA

    def _min_valid_iso_date(vals: Sequence[object]) -> object:
      valid_dts = [
          dt for v in vals if (dt := _parse_valid_iso_date(v)) is not None
      ]
      return min(valid_dts) if valid_dts else np.nan

    for y_col in ("document_year", "publication_year", "chembl_year", "year"):
      if y_col in struct_part.columns:
        struct_part[y_col] = (
            struct_part.groupby("_dedup_key", sort=False)[y_col]
            .transform(_min_valid_year)
            .astype("Int64")
        )
    dedup_grp_sizes = (
        struct_part.groupby("_dedup_key", sort=False)["_dedup_key"]
        .transform("size")
        .tolist()
    )
    for d_col in ("publication_date", "filing_date", "priority_date"):
      if d_col in struct_part.columns:
        struct_part[d_col] = struct_part.groupby("_dedup_key", sort=False)[
            d_col
        ].transform(_min_valid_iso_date)
        if "document_year" in struct_part.columns:
          reconciled_col = []
          for g_sz, iso_v, yr_v in zip(
              dedup_grp_sizes,
              struct_part[d_col].tolist(),
              struct_part["document_year"].tolist(),
              strict=True,
          ):
            if (
                g_sz > 1
                and isinstance(iso_v, str)
                and len(iso_v) >= 4
                and pd.notna(yr_v)
                and int(iso_v[:4]) > int(yr_v)
            ):
              reconciled_col.append(np.nan)
            else:
              reconciled_col.append(iso_v)
          struct_part[d_col] = pd.Series(
              reconciled_col, index=struct_part.index, dtype=object
          )

    structured = (
        struct_part.groupby("_dedup_key", sort=False, as_index=False)
        .first()
        .set_index("_orig_idx")
    )
  else:
    structured = struct_part.drop(columns=["_orig_idx"], errors="ignore")

  structureless = ordered[~ordered["_has_structure"]].drop(
      columns=["_orig_idx"], errors="ignore"
  )

  if not structured.empty and not structureless.empty:
    # Backfill any missing assay / potency values on structured rows that share
    # the exact patent-qualified patent_example_id or compound_id (or bare
    # Example_N ONLY when one of the rows is unprefixed, never across two
    # different patent_id prefixes).
    for id_key in ("patent_example_id", "compound_id"):
      if id_key in structured.columns and id_key in structureless.columns:
        sless_by_exact: dict[str, pd.Series] = {}
        sless_by_unprefixed: dict[str, pd.Series] = {}
        sless_by_bare_any: dict[str, pd.Series] = {}
        for _, s_row in structureless.dropna(subset=[id_key]).iterrows():
          for raw_tok in str(s_row.get(id_key, "")).split(","):
            tok = raw_tok.strip()
            if not tok:
              continue
            sless_by_exact.setdefault(tok, s_row)
            if ":" in tok:
              bare_tok = tok.split(":", 1)[1].strip()
              if bare_tok:
                sless_by_bare_any.setdefault(bare_tok, s_row)
            else:
              sless_by_unprefixed.setdefault(tok, s_row)
              sless_by_bare_any.setdefault(tok, s_row)
        if sless_by_exact:
          for idx_s, row_s in structured.iterrows():
            kv = row_s.get(id_key)
            if pd.notna(kv):
              for raw_tok in str(kv).split(","):
                tok = raw_tok.strip()
                if not tok:
                  continue
                if ":" in tok:
                  # Both sides have patent prefixes unless sless_by_unprefixed matched
                  cand_row = sless_by_exact.get(tok) or sless_by_unprefixed.get(
                      tok.split(":", 1)[1].strip()
                  )
                else:
                  cand_row = sless_by_exact.get(tok) or sless_by_bare_any.get(
                      tok
                  )
                if cand_row is not None:
                  for col in structured.columns:
                    if col in (
                        "_smiles_key",
                        "_dedup_key",
                        "_has_structure",
                        "_has_act",
                        "_source_rank",
                        "_is_neutral_ik",
                    ):
                      continue
                    if pd.isna(structured.at[idx_s, col]) and pd.notna(
                        cand_row.get(col)
                    ):
                      val_to_set = cand_row.get(col)
                      try:
                        structured.at[idx_s, col] = val_to_set
                      except (TypeError, ValueError):
                        structured[col] = structured[col].astype(object)
                        structured.at[idx_s, col] = val_to_set
    deduped = structured.sort_index()
  else:
    deduped = pd.concat([structured, structureless]).sort_index()
  if "compound_name" in deduped.columns:
    for idx_d, row_d in deduped.iterrows():
      cname_val = (
          str(row_d.get("compound_name", "")).strip()
          if pd.notna(row_d.get("compound_name"))
          else ""
      )
      if len(cname_val) > 28:
        if (
            "iupac_name" not in deduped.columns
            or pd.isna(row_d.get("iupac_name"))
            or not str(row_d.get("iupac_name", "")).strip()
        ):
          deduped.at[idx_d, "iupac_name"] = cname_val
        short_alt = ""
        p_id = (
            str(row_d.get("paper_compound_id", "")).strip()
            if pd.notna(row_d.get("paper_compound_id"))
            else ""
        )
        if p_id and p_id.lower() not in ("nan", "none") and len(p_id) <= 20:
          is_code = bool(
              re.match(
                  r"^(?:cmpd|compound|example|ex\b|chembl|[a-z]{2,}[-_]?\d)",
                  p_id,
                  re.IGNORECASE,
              )
          )
          short_alt = p_id if is_code else f"Cmpd {p_id}"
        if not short_alt:
          for cand_col in (
              "patent_example_id",
              "compound_id",
              "chembl_id",
              "molecule_chembl_id",
          ):
            if cand_col in deduped.columns and pd.notna(row_d.get(cand_col)):
              cv = str(row_d.get(cand_col, "")).strip()
              if cv and cv.lower() not in ("nan", "none") and len(cv) <= 24:
                short_alt = cv
                break
        if short_alt:
          deduped.at[idx_d, "compound_name"] = short_alt
  return deduped.drop(
      columns=[
          "_smiles_key",
          "_flat_smiles_key",
          "_has_defined_stereo",
          "_clean_peak",
          "_dedup_key",
          "_has_structure",
          "_has_act",
          "_source_rank",
          "_is_neutral_ik",
          "_doc_size",
      ],
      errors="ignore",
  ).reset_index(drop=True)


def _print_dedup_summary(
    before: pd.DataFrame, after: pd.DataFrame, smiles_col: str
) -> None:
  """Prints input, removal, missing-structure, and per-source dedup counts."""
  key = (
      before[smiles_col]
      .astype(str)
      .where(before[smiles_col].notna(), "")
      .str.strip()
  )
  n_missing = int((key == "").sum())
  print(f"      Rows in: {len(before)}")
  print(f"      Rows removed by SMILES dedup: {len(before) - len(after)}")
  print(
      f"      Rows with missing/empty SMILES: {n_missing} (merged onto"
      " matching Example IDs / omitted)"
  )
  if "source" not in before.columns:
    print("      Per-source counts unavailable: no 'source' column.")
    return
  counts_before = before["source"].fillna("(unlabelled)").value_counts()
  counts_after = (
      after["source"].fillna("(unlabelled)").value_counts()
      if "source" in after.columns
      else pd.Series(dtype=int)
  )
  print("      Per-source row counts (before -> after):")
  for src in sorted(counts_before.index.astype(str)):
    print(
        f"        {src}: {int(counts_before.get(src, 0))} ->"
        f" {int(counts_after.get(src, 0))}"
    )


def _featurize_dataframe(
    master_df: pd.DataFrame, s_col: str
) -> Tuple[pd.DataFrame, list[Optional[Chem.Mol]]]:
  """Computes canonical RDKit descriptors via descriptors.compute_standard_descriptors for rows."""
  props_records = []
  mols = []
  for smi in master_df[s_col]:
    if not isinstance(smi, str) or not smi.strip():
      mols.append(None)
      props_records.append({})
      continue
    m = Chem.MolFromSmiles(smi)
    if m is None:
      mols.append(None)
      props_records.append({})
      continue

    mols.append(m)
    try:
      rec = dict(descriptors.compute_standard_descriptors(m, style="short"))
    except Exception:
      rec = {}
    try:
      rec["calculated_inchikey"] = Chem.MolToInchiKey(m)
    except Exception:
      rec["calculated_inchikey"] = None
    props_records.append(rec)

  props_df = pd.DataFrame(props_records)
  for c in props_df.columns:
    master_df[c] = props_df[c]
  return master_df, mols


def _print_cohort_summary_table(
    df: pd.DataFrame,
    group_by: Optional[str] = None,
    panels: Optional[Sequence[Union[str, dict[str, object]]]] = None,
) -> None:
  """Prints a structured Markdown summary table of descriptors & potency across cohorts to stdout."""
  _, groups = _resolve_grouping(df, group_by)
  if len(groups) > 1:
    groups = [
        ("All Compounds", pd.Series(True, index=df.index), "#000000")
    ] + list(groups)

  disc_bin_col = None
  n_disc_total = 0
  for dc in ["activity_bin", "activity_value"] + [
      c
      for c in df.columns
      if str(c).startswith("patent_assay_") and not str(c).endswith("_range")
  ]:
    if _is_discrete_potency_col(df, dc):
      cnt_d = _count_discrete_potency_values(df, dc)
      if cnt_d > n_disc_total:
        disc_bin_col = dc
        n_disc_total = cnt_d

  pot_col = None
  for cand in ("pchembl_value", "pActivity", "pKi", "pIC50", "potency"):
    if cand in df.columns:
      n_num = int(pd.to_numeric(df[cand], errors="coerce").notna().sum())
      if n_num >= 2 or (n_num >= 1 and disc_bin_col is None):
        pot_col = cand
        break

  per_endpoint_pcols = [
      c
      for c in df.columns
      if c != pot_col
      and re.match(r"^p(Ki|IC50|Kd|EC50)_", str(c))
      and int(pd.to_numeric(df[c], errors="coerce").notna().sum()) >= 5
  ]

  panels_str = " ".join(str(p) for p in panels or [])
  base_metric_cands = [
      "MW",
      "cLogP",
      "TPSA",
      "HBD",
      "HBA",
      "RotB",
      "Fsp3",
      "QED",
      "SA_Score",
  ]
  if "CNS_MPO" in panels_str:
    base_metric_cands.append("CNS_MPO")

  metrics = [c for c in base_metric_cands if c in df.columns]
  if pot_col:
    metrics = [pot_col] + per_endpoint_pcols + metrics
  elif per_endpoint_pcols:
    metrics = per_endpoint_pcols + metrics

  display_labels = {
      m: (
          "pchembl_value (Pooled)"
          if (m == "pchembl_value" and len(per_endpoint_pcols) >= 2)
          else m
      )
      for m in metrics
  }

  col_headers = ["Property"]
  cohort_stats: list[dict[str, str]] = []
  for label, mask, _ in groups:
    sub = df[mask]
    n_tot = len(sub)
    has_num_act = (
        pd.to_numeric(sub[pot_col], errors="coerce").notna()
        if pot_col
        else pd.Series(False, index=sub.index)
    )
    if disc_bin_col:
      has_disc_act = sub[disc_bin_col].notna() & sub[disc_bin_col].astype(
          str
      ).str.strip().str.lower().replace({"nan": "", "none": ""}).ne("")
    else:
      has_disc_act = pd.Series(False, index=sub.index)
    n_act = int((has_num_act | has_disc_act).sum())
    col_headers.append(f"{label} (N={n_tot}, Assayed={n_act})")
    stat_map: dict[str, str] = {"Assayed Compounds": f"{n_act} / {n_tot}"}
    if disc_bin_col:
      vc_disc = (
          sub.loc[has_disc_act, disc_bin_col]
          .astype(str)
          .str.strip()
          .value_counts()
      )
      if not vc_disc.empty:
        sorted_bins = sorted(
            vc_disc.index.tolist(),
            key=lambda b: (
                0 if re.match(r"^[A-F](?:\s*\(.*\))?$", b, re.I) else 1,
                b.upper(),
            ),
        )
        stat_map["Discrete Potency Bins"] = ", ".join(
            f"{b}: {int(vc_disc[b])}" for b in sorted_bins[:6]
        )
      else:
        stat_map["Discrete Potency Bins"] = "—"
    for m in metrics:
      s = pd.to_numeric(sub[m], errors="coerce").dropna()
      disp_m = display_labels.get(m, m)
      if len(s) == 0:
        stat_map[disp_m] = "—"
      else:
        q25, med, q75 = np.percentile(s, [25, 50, 75])
        stat_map[disp_m] = (
            f"{s.mean():.2f} ± {s.std():.2f} (med {med:.2f}, IQR"
            f" {q75 - q25:.2f})"
        )
    cohort_stats.append(stat_map)

  print("\n### Cohort Descriptor & Potency Summary (Mean ± SD | Median, IQR)\n")
  print("| " + " | ".join(col_headers) + " |")
  print("| " + " | ".join(":---" for _ in col_headers) + " |")
  all_row_keys = ["Assayed Compounds"]
  if disc_bin_col:
    all_row_keys.append("Discrete Potency Bins")
  all_row_keys.extend(display_labels.get(m, m) for m in metrics)
  for rk in all_row_keys:
    row_cells = [f"**{rk}**"] + [cs.get(rk, "—") for cs in cohort_stats]
    print("| " + " | ".join(row_cells) + " |")
  print("", flush=True)


def _is_purely_numeric_series(s: pd.Series) -> bool:
  """Returns True when every non-empty value in `s` parses as a finite float."""
  has_val = False
  for v in s:
    if v is None or pd.isna(v):
      continue
    sv = str(v).strip()
    if not sv or sv.lower() in ("nan", "none"):
      continue
    has_val = True
    try:
      fv = float(sv)
    except ValueError:
      return False
    if not np.isfinite(fv):
      return False
  return has_val


def _resolve_col_candidate(df: pd.DataFrame, cand: str) -> str | None:
  """Resolves a column name candidate against df.columns (exact, with/without _nM, or case-insensitive)."""
  c = cand.strip()
  if not c:
    return None
  if c in df.columns:
    return c
  if f"{c}_nM" in df.columns:
    return f"{c}_nM"
  if (
      c.lower().endswith("_nm")
      and c[:-3] in df.columns
      and _is_purely_numeric_series(df[c[:-3]])
  ):
    return c[:-3]
  low_map = {str(col).lower(): str(col) for col in df.columns}
  if c.lower() in low_map:
    return low_map[c.lower()]
  if f"{c.lower()}_nm" in low_map:
    return low_map[f"{c.lower()}_nm"]
  if (
      c.lower().endswith("_nm")
      and c.lower()[:-3] in low_map
      and _is_purely_numeric_series(df[low_map[c.lower()[:-3]]])
  ):
    return low_map[c.lower()[:-3]]
  return None


def _is_named_assay_col(col: object) -> bool:
  """Returns True when col is a named potency/assay column."""
  cl = str(col).strip().lower()
  if cl.endswith("_range") or _ANON_PATENT_COL_RE.match(cl):
    return False
  if cl.startswith("patent_assay_"):
    return len(cl) > len("patent_assay_")
  if _META_UNIT_SUFFIX_RE.search(_BIO_TYPE_TOKEN_RE.sub("", cl)):
    return False
  return bool(_NAMED_POTENCY_COL_RE.search(cl))


def _is_raw_patent_source(name: str) -> bool:
  """Returns True for raw patent_assay_* or positional patent column names."""
  nl = name.strip().lower()
  return nl.startswith(("patent_assay_", "patent_table_col_")) or bool(
      _ANON_PATENT_COL_RE.match(nl)
  )


def _potency_sigs_match(
    sig_a: tuple[float | None, str, str],
    sig_b: tuple[float | None, str, str],
) -> bool:
  """Returns True when two potency signatures represent the same value and qualifier."""
  pa, qa, sa = sig_a
  pb, qb, sb = sig_b
  if pa is not None and pb is not None:
    return qa == qb and abs(pa - pb) <= _PACTIVITY_MATCH_TOL
  if pa is None and pb is None:
    return sa == sb
  return False


def _authoritative_sigs_for_cols(
    cdf: pd.DataFrame, idx_r: Hashable, cols: Sequence[str]
) -> list[tuple[float | None, str, str]]:
  """Returns row potency signatures preferring unsuffixed raw columns over derived floats."""
  raw_sigs: list[tuple[float | None, str, str]] = []
  derived_sigs: list[tuple[float | None, str, str]] = []
  for col in cols:
    if col not in cdf.columns or col.lower().endswith("_range"):
      continue
    sig = _cell_potency_sig(cdf.at[idx_r, col], col)
    if sig is None:
      continue
    cl = col.lower()
    if not _ASSAY_SUFFIX_RE.search(cl) and not _LOG_POTENCY_COL_RE.search(cl):
      raw_sigs.append(sig)
    else:
      derived_sigs.append(sig)
  return raw_sigs if raw_sigs else derived_sigs


def _authoritative_pos_sigs_for_cols(
    cdf: pd.DataFrame, idx_r: Hashable, cols: Sequence[str]
) -> list[tuple[float | None, str, str]]:
  """Returns authoritative signatures per positional prefix on row `idx_r`."""
  by_prefix: dict[str, list[str]] = {}
  for col in cols:
    if col.lower().endswith("_range"):
      continue
    pfx = _ASSAY_SUFFIX_RE.sub("", col).lower()
    by_prefix.setdefault(pfx, []).append(col)
  out: list[tuple[float | None, str, str]] = []
  for pfx_cols in by_prefix.values():
    out.extend(_authoritative_sigs_for_cols(cdf, idx_r, pfx_cols))
  return out


def _has_injective_stem_match(
    anon_sig_groups: Sequence[Sequence[tuple[float | None, str, str]]],
    named_sig_groups: Sequence[Sequence[tuple[float | None, str, str]]],
) -> bool:
  """Returns True iff every positional group matches a distinct named stem."""
  if len(named_sig_groups) < len(anon_sig_groups):
    return False
  adj: list[list[int]] = []
  for asigs in anon_sig_groups:
    matches = [
        j
        for j, nsigs in enumerate(named_sig_groups)
        if all(any(_potency_sigs_match(ps, ns) for ns in nsigs) for ps in asigs)
    ]
    if not matches:
      return False
    adj.append(matches)
  match_right = [-1] * len(named_sig_groups)

  def _dfs(i: int, seen: list[bool]) -> bool:
    for j in adj[i]:
      if seen[j]:
        continue
      seen[j] = True
      if match_right[j] == -1 or _dfs(match_right[j], seen):
        match_right[j] = i
        return True
    return False

  for i in range(len(anon_sig_groups)):
    if not _dfs(i, [False] * len(named_sig_groups)):
      return False
  return True


def _has_nonempty_cells(s: pd.Series) -> bool:
  """Returns True if `s` has at least one non-empty value."""
  for v in s:
    if v is None or pd.isna(v):
      continue
    sv = str(v).strip()
    if sv and sv.lower() not in ("nan", "none"):
      return True
  return False


def _disambiguate_and_rename_input_columns(
    cdf: pd.DataFrame,
    source_path: str = "",
    rename_cols: Sequence[str] | None = None,
) -> pd.DataFrame:
  """Qualifies patent_table_col_*, coalesces --col-map, and prunes duplicates."""
  cdf = cdf.copy()
  orig_cols_by_pos: dict[int, list[str]] = {}
  for c in cdf.columns:
    m_orig = _ANON_PATENT_COL_RE.match(str(c))
    if not m_orig or str(c).lower().endswith("_range"):
      continue
    orig_cols_by_pos.setdefault(int(m_orig.group(1)), []).append(str(c))
  orig_sigs_by_row: dict[
      Hashable, dict[int, list[tuple[float | None, str, str]]]
  ] = {}
  for idx_r in cdf.index:
    for c_idx, pcols in orig_cols_by_pos.items():
      psigs = _authoritative_pos_sigs_for_cols(cdf, idx_r, pcols)
      if psigs:
        orig_sigs_by_row.setdefault(idx_r, {})[c_idx] = psigs

  # 1. Qualify generic patent_table_col_{k}_nM columns by patent_id when present
  generic_pcols = [
      c for c in cdf.columns if _GENERIC_PATENT_COL_RE.match(str(c))
  ]
  if generic_pcols and "patent_id" in cdf.columns:
    numeric_gcols = {
        gcol
        for gcol in generic_pcols
        if str(gcol).endswith("_nM")
        or pd.to_numeric(cdf[gcol], errors="coerce").notna().any()
    }
    for idx_r, row_r in cdf.iterrows():
      pid = str(row_r.get("patent_id") or "").strip().upper()
      if not pid or pid in ("NAN", "NONE"):
        continue
      pid_tokens = [
          t
          for tok in re.split(r"[,;\s]+", pid)
          if (t := _NON_ALNUM_RE.sub("", tok))
      ]
      if not pid_tokens:
        continue
      doc_tok = _NON_ALNUM_RE.sub(
          "", str(row_r.get("document_id") or "").strip().upper()
      )
      safe_pid = doc_tok if doc_tok in pid_tokens else pid_tokens[0]
      for gcol in generic_pcols:
        val_nm = row_r.get(gcol)
        if pd.isna(val_nm):
          continue
        col_num = int(_GENERIC_PATENT_COL_RE.match(str(gcol)).group(1))
        if (
            len(pid_tokens) > 1
            and doc_tok not in pid_tokens
            and any(
                f"{t}_col_{col_num}" in cdf.columns
                and pd.notna(cdf.at[idx_r, f"{t}_col_{col_num}"])
                for t in pid_tokens
            )
        ):
          continue
        q_cols = (
            (f"{safe_pid}_col_{col_num}_nM", f"{safe_pid}_col_{col_num}")
            if gcol in numeric_gcols
            else (f"{safe_pid}_col_{col_num}",)
        )
        for q_col in q_cols:
          if q_col not in cdf.columns:
            cdf[q_col] = pd.Series(index=cdf.index, dtype=object)
          if pd.isna(cdf.at[idx_r, q_col]):
            cdf.at[idx_r, q_col] = val_nm

  # 2. Apply --rename-cols / --col-map rules with multi-column coalescing
  coalesce_groups: dict[str, list[str]] = {}
  touched_anon_cols: set[str] = set()
  if rename_cols:
    base_fn = os.path.basename(source_path).lower() if source_path else ""
    for raw_rule in rename_cols:
      for part in str(raw_rule).split(","):
        rule = part.strip()
        if "=" not in rule:
          continue
        lhs, rhs = [x.strip() for x in rule.split("=", 1)]
        has_file_prefix = ":" in lhs
        if has_file_prefix:
          f_pat, lhs = [x.strip() for x in lhs.split(":", 1)]
          if not source_path or not (
              f_pat.lower() in base_fn or f_pat.lower() in source_path.lower()
          ):
            continue

        if "|" in rhs:
          for src_cand in rhs.split("|"):
            resolved = _resolve_col_candidate(cdf, src_cand)
            if resolved and resolved not in coalesce_groups.setdefault(lhs, []):
              coalesce_groups[lhs].append(resolved)
          continue

        res_lhs = _resolve_col_candidate(cdf, lhs)
        res_rhs = _resolve_col_candidate(cdf, rhs)
        if has_file_prefix or (
            _is_raw_patent_source(lhs) and not _is_raw_patent_source(rhs)
        ):
          if res_lhs is None:
            continue
          target_col, src_col = rhs, res_lhs
        elif res_rhs is not None:
          target_col, src_col = lhs, res_rhs
        else:
          continue
        if src_col not in coalesce_groups.setdefault(target_col, []):
          coalesce_groups[target_col].append(src_col)

  protected_targets: set[str] = set(coalesce_groups.keys())
  if coalesce_groups:
    merged_by_target: dict[str, pd.Series] = {}
    all_drop_cols: list[str] = []
    consumed_pos_targets: list[tuple[int, str, str]] = []
    for target_col, src_cols in coalesce_groups.items():
      if not src_cols:
        continue
      chain_cols = (
          [target_col] + [c for c in src_cols if c != target_col]
          if target_col in cdf.columns
          else list(src_cols)
      )
      merged_s = cdf[chain_cols[0]].copy()
      for next_c in chain_cols[1:]:
        merged_s = merged_s.combine_first(cdf[next_c])
      merged_by_target[target_col] = merged_s

      m_tgt_sfx = re.search(r"_(?:nM|pActivity)$", target_col, flags=re.I)
      tgt_stem = target_col[: m_tgt_sfx.start()] if m_tgt_sfx else target_col
      tgt_rng_col = f"{tgt_stem}_range"
      raw_by_src: dict[str, str] = {}
      rng_by_src: dict[str, str] = {}
      if (
          target_col in cdf.columns
          and tgt_stem != target_col
          and tgt_stem in cdf.columns
      ):
        raw_by_src[target_col] = tgt_stem
      if target_col in cdf.columns and tgt_rng_col in cdf.columns:
        rng_by_src[target_col] = tgt_rng_col

      for c in src_cols:
        if m_pos := _ANON_PATENT_COL_RE.match(str(c)):
          consumed_pos_targets.append(
              (int(m_pos.group(1)), target_col, tgt_stem)
          )
        if c != target_col and c in cdf.columns:
          c_stem = _ASSAY_SUFFIX_RE.sub("", str(c)).lower()
          for sib in cdf.columns:
            sib_s = str(sib)
            if sib_s == target_col or sib_s in coalesce_groups:
              continue
            if _ASSAY_SUFFIX_RE.sub("", sib_s).lower() != c_stem:
              continue
            sib_low = sib_s.lower()
            if sib_low.endswith("_range"):
              rng_by_src.setdefault(str(c), sib_s)
            elif not _ASSAY_SUFFIX_RE.search(
                sib_low
            ) and not _LOG_POTENCY_COL_RE.search(sib_low):
              raw_by_src.setdefault(str(c), sib_s)
            if not _ANON_PATENT_COL_RE.match(sib_s) or sib_s == str(c):
              all_drop_cols.append(sib_s)

      if (
          m_tgt_sfx is not None
          and raw_by_src
          and tgt_stem not in coalesce_groups
      ):
        raw_acc = (
            cdf[tgt_stem].copy()
            if tgt_stem in cdf.columns
            else pd.Series(index=cdf.index, dtype=object)
        )
        claimed = pd.Series(False, index=cdf.index)
        for col_i in chain_cols:
          win_mask = ~claimed & cdf[col_i].notna()
          raw_i = raw_by_src.get(col_i)
          if raw_i is not None:
            raw_acc = raw_acc.where(~win_mask, cdf[raw_i])
          claimed = claimed | cdf[col_i].notna()
        unclaimed = ~claimed
        if unclaimed.any():
          for col_i in chain_cols:
            raw_i = raw_by_src.get(col_i)
            if raw_i is not None:
              fill_mask = unclaimed & raw_acc.isna() & cdf[raw_i].notna()
              raw_acc = raw_acc.where(~fill_mask, cdf[raw_i])
        if any(_is_nonexact_raw_potency_cell(v) for v in raw_acc):
          if tgt_stem in merged_by_target:
            merged_by_target[tgt_stem] = merged_by_target[
                tgt_stem
            ].combine_first(raw_acc)
          else:
            merged_by_target[tgt_stem] = raw_acc
          protected_targets.add(tgt_stem)

      if rng_by_src and tgt_rng_col not in coalesce_groups:
        rng_acc = (
            cdf[tgt_rng_col].copy()
            if tgt_rng_col in cdf.columns
            else pd.Series(index=cdf.index, dtype=object)
        )
        claimed_rng = pd.Series(False, index=cdf.index)
        for col_i in chain_cols:
          win_mask = ~claimed_rng & cdf[col_i].notna()
          rng_i = rng_by_src.get(col_i)
          if rng_i is not None:
            rng_acc = rng_acc.where(~win_mask, cdf[rng_i])
          claimed_rng = claimed_rng | cdf[col_i].notna()
        unclaimed_rng = ~claimed_rng
        if unclaimed_rng.any():
          for col_i in chain_cols:
            rng_i = rng_by_src.get(col_i)
            if rng_i is not None:
              fill_mask = unclaimed_rng & rng_acc.isna() & cdf[rng_i].notna()
              rng_acc = rng_acc.where(~fill_mask, cdf[rng_i])
        if _has_nonempty_cells(rng_acc):
          if tgt_rng_col in merged_by_target:
            merged_by_target[tgt_rng_col] = merged_by_target[
                tgt_rng_col
            ].combine_first(rng_acc)
          else:
            merged_by_target[tgt_rng_col] = rng_acc
          protected_targets.add(tgt_rng_col)

    if merged_by_target:
      cdf = cdf.assign(**merged_by_target)
    if all_drop_cols:
      drop_unique = [
          c for c in dict.fromkeys(all_drop_cols) if c not in protected_targets
      ]
      if drop_unique:
        cdf = cdf.drop(columns=drop_unique, errors="ignore")
    for col_idx, tgt_name, tgt_stem in consumed_pos_targets:
      sib_pos = [
          str(c)
          for c in cdf.columns
          if c not in protected_targets
          and (m := _ANON_PATENT_COL_RE.match(str(c)))
          and int(m.group(1)) == col_idx
      ]
      tgt_check_cols = (
          [tgt_stem, tgt_name]
          if tgt_stem != tgt_name and tgt_stem in cdf.columns
          else [tgt_name]
      )
      for idx_r in cdf.index:
        tgt_sigs = _authoritative_sigs_for_cols(cdf, idx_r, tgt_check_cols)
        if not tgt_sigs:
          continue
        row_sigs = _authoritative_pos_sigs_for_cols(cdf, idx_r, sib_pos)
        if row_sigs and all(
            any(_potency_sigs_match(s, ts) for ts in tgt_sigs) for s in row_sigs
        ):
          for sc in sib_pos:
            cdf.at[idx_r, sc] = np.nan
            touched_anon_cols.add(sc)

  # 3. Prune anonymous positional cells on rows fully covered by named assays
  anon_by_idx: dict[int, list[str]] = {}
  for c in cdf.columns:
    if c not in protected_targets and (m := _ANON_PATENT_COL_RE.match(str(c))):
      anon_by_idx.setdefault(int(m.group(1)), []).append(str(c))

  if anon_by_idx:
    named_by_stem: dict[str, list[str]] = {}
    for c in cdf.columns:
      if _is_named_assay_col(c):
        raw_stem = re.sub(r"^patent_assay_", "", str(c), flags=re.I).lower()
        norm_stem = _POTENCY_STEM_NORM_RE.sub(
            r"\1\2", _ASSAY_SUFFIX_RE.sub("", raw_stem)
        )
        stem_k = "_".join(sorted(t for t in norm_stem.split("_") if t))
        named_by_stem.setdefault(stem_k, []).append(str(c))
    if named_by_stem:
      for idx_r in cdf.index:
        orig_row_pos = orig_sigs_by_row.get(idx_r, {})
        if not orig_row_pos:
          continue
        named_sig_groups = [
            stem_sigs
            for scols in named_by_stem.values()
            if (stem_sigs := _authoritative_sigs_for_cols(cdf, idx_r, scols))
        ]
        if len(named_sig_groups) < len(orig_row_pos):
          continue
        req_by_idx: dict[int, list[tuple[float | None, str, str]]] = {
            c_idx: list(psigs) for c_idx, psigs in orig_row_pos.items()
        }
        active_anon_cols: list[str] = []
        for c_idx, idx_cols in anon_by_idx.items():
          idx_sigs = _authoritative_pos_sigs_for_cols(cdf, idx_r, idx_cols)
          if not idx_sigs:
            continue
          active_anon_cols.extend(idx_cols)
          cur_list = req_by_idx.setdefault(c_idx, [])
          for isig in idx_sigs:
            if not any(_potency_sigs_match(isig, esig) for esig in cur_list):
              cur_list.append(isig)
        if active_anon_cols and _has_injective_stem_match(
            list(req_by_idx.values()), named_sig_groups
        ):
          for ac in active_anon_cols:
            cdf.at[idx_r, ac] = np.nan
            touched_anon_cols.add(ac)

  if touched_anon_cols:
    empty_anon = [
        c
        for c in cdf.columns
        if c in touched_anon_cols
        and cdf[c]
        .dropna()
        .astype(str)
        .str.strip()
        .str.lower()
        .isin(("", "nan", "none"))
        .all()
    ]
    if empty_anon:
      cdf = cdf.drop(columns=empty_anon, errors="ignore")

  return cdf


def _load_any_chemical_input(
    path: str, rename_cols: Optional[Sequence[str]] = None
) -> list[pd.DataFrame]:
  """Loads CSV, JSON, or SDF file (or directory of such files) into DataFrames."""
  loaded: list[pd.DataFrame] = []
  if os.path.isdir(path):
    for entry in sorted(os.listdir(path)):
      full_p = os.path.join(path, entry)
      if os.path.isfile(full_p) and entry.lower().endswith(
          (".csv", ".json", ".sdf")
      ):
        loaded.extend(_load_any_chemical_input(full_p, rename_cols=rename_cols))
    return loaded

  if not os.path.exists(path):
    print(f"  Warning: File not found: {path}")
    return loaded

  low = path.lower()
  print(f"  [+] Ingesting input source: {path}...", flush=True)
  if low.endswith(".csv"):
    cdf = pd.read_csv(path)
    cdf = _disambiguate_and_rename_input_columns(
        cdf, path, rename_cols=rename_cols
    )
    cdf = descriptors.normalize_descriptor_columns(cdf, target_style="short")
    loaded.append(cdf)
  elif low.endswith(".json"):
    with open(path, "r", encoding="utf-8") as fh:
      jdata = json.load(fh)
    if isinstance(jdata, dict):
      for k in ("compounds", "records", "molecules", "data", "activities"):
        if isinstance(jdata.get(k), list):
          jdata = jdata[k]
          break
    if isinstance(jdata, list):
      cdf = pd.DataFrame(jdata)
      if "molecule_chembl_id" in cdf.columns and "chembl_id" not in cdf.columns:
        cdf["chembl_id"] = cdf["molecule_chembl_id"]
      if "compound_id" not in cdf.columns and "chembl_id" in cdf.columns:
        cdf["compound_id"] = cdf["chembl_id"]
      if "target_pref_name" in cdf.columns and "target_name" not in cdf.columns:
        cdf["target_name"] = cdf["target_pref_name"]
      if "source" not in cdf.columns and "document_chembl_id" in cdf.columns:
        cdf["source"] = "ChEMBL"
      if "pchembl_value" in cdf.columns:
        cdf["pchembl_value"] = pd.to_numeric(
            cdf["pchembl_value"], errors="coerce"
        )
      if "standard_value" in cdf.columns:
        num_sv = pd.to_numeric(cdf["standard_value"], errors="coerce")
        if num_sv.notna().any():
          cdf["standard_value"] = cdf["standard_value"].where(
              num_sv.isna(), num_sv
          )
      cdf = _disambiguate_and_rename_input_columns(
          cdf, path, rename_cols=rename_cols
      )
      cdf = descriptors.normalize_descriptor_columns(cdf, target_style="short")
      loaded.append(cdf)
  elif low.endswith(".sdf"):
    suppl = Chem.SDMolSupplier(path, removeHs=False)
    sdf_rows = []
    for idx, m in enumerate(suppl):
      if m is None:
        continue
      row_d = m.GetPropsAsDict()
      row_d.setdefault(
          "compound_id",
          m.GetProp("_Name")
          if m.HasProp("_Name") and m.GetProp("_Name")
          else f"SDF_{idx+1}",
      )
      row_d["canonical_smiles"] = Chem.MolToSmiles(m, isomericSmiles=True)
      sdf_rows.append(row_d)
    if sdf_rows:
      cdf = pd.DataFrame(sdf_rows)
      cdf = _disambiguate_and_rename_input_columns(
          cdf, path, rename_cols=rename_cols
      )
      cdf = descriptors.normalize_descriptor_columns(cdf, target_style="short")
      loaded.append(cdf)
  return loaded


def _apply_corpus_structural_filters_and_series(
    df: pd.DataFrame,
    s_col: str,
    min_tanimoto: float = 0.40,
    query_smarts: Optional[str] = None,
    query_scaffold: Optional[str] = None,
    series_smarts: Optional[str] = None,
) -> pd.DataFrame:
  """Applies Murcko scaffold + Morgan ECFP4 Tanimoto (0.40) applicability domain

  filtering to un-exemplified associative records, optional SMARTS/scaffold
  filters,
  and optional --series-smarts ('Label1=SMARTS1,Label2=SMARTS2') cohort tagging.
  """
  fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
  smarts_pat = Chem.MolFromSmarts(query_smarts) if query_smarts else None
  scaf_pat = Chem.MolFromSmiles(query_scaffold) if query_scaffold else None
  if query_scaffold and scaf_pat is None:
    scaf_pat = Chem.MolFromSmarts(query_scaffold)

  synth_reagent_pat = Chem.MolFromSmarts("[#5,#14,#50]")
  boc_pat = Chem.MolFromSmarts("CC(C)(C)OC(=O)N")
  cbz_pat = Chem.MolFromSmarts("OC(=O)NCc1ccccc1")
  aryl_br_i_pat = Chem.MolFromSmarts("[c,n][Br,I]")

  parsed_series: list[tuple[str, str, Any]] = []
  if series_smarts:
    for item in re.split(r"[;,](?=[A-Za-z0-9_ -]+=)", str(series_smarts)):
      if "=" in item:
        lbl, smt = item.split("=", 1)
        smt_str = smt.strip()
        smt_relaxed = re.sub(r"\[([nN])[hH]\d*\]", r"\1", smt_str)
        pat = (
            Chem.MolFromSmarts(smt_relaxed)
            or Chem.MolFromSmiles(smt_relaxed)
            or Chem.MolFromSmarts(smt_str)
            or Chem.MolFromSmiles(smt_str)
        )
        if pat is not None or smt_str:
          parsed_series.append((lbl.strip(), smt_str, pat))

  # Identify anchor/exemplified compounds (ChEMBL, Literature, or Patent_FullText_Table / patent_example_id)
  anchor_scaffolds = set()
  anchor_fps = []
  anchor_ik_to_indices: dict[str, list[Any]] = {}
  anchors_have_boc = False
  anchors_have_aryl_br_i = False
  pat_has_examples: set[str] = set()
  pat_parent_map: dict[str, set[str]] = {}
  has_unexemplified_assoc = False
  for idx, row in df.iterrows():
    src = str(row.get("source") or "")
    all_src = str(row.get("all_sources") or src)
    ex_id = row.get("patent_example_id")
    is_assoc = (
        src in ("PubChem", "Patent_Concept_OCSR")
        and "Patent_FullText_Table" not in all_src
        and (pd.isna(ex_id) or not str(ex_id).strip())
    )
    rel_str = str(row.get("patent_relationship") or "")
    if rel_str and rel_str.lower() not in ("", "nan", "none"):
      parent_matches = re.findall(
          r"from\s+([A-Z]{2}\d+[A-Z0-9]*)", rel_str, flags=re.IGNORECASE
      )
      if parent_matches:
        for pat_tok in str(row.get("patent_id") or "").split(","):
          pt = pat_tok.strip().upper()
          if pt and pt not in ("NAN", "NONE"):
            for pm in parent_matches:
              pat_parent_map.setdefault(pt, set()).add(pm.strip().upper())
    if is_assoc:
      has_unexemplified_assoc = True
    else:
      for pat_tok in str(row.get("patent_id") or "").split(","):
        pt = pat_tok.strip().upper()
        if pt and pt not in ("NAN", "NONE"):
          pat_has_examples.add(pt)
      raw_s = row.get(s_col)
      smi = (
          str(raw_s).strip()
          if pd.notna(raw_s)
          and str(raw_s).strip().lower() not in ("", "nan", "none")
          else ""
      )
      m = Chem.MolFromSmiles(smi) if smi else None
      if m is not None:
        try:
          ik_full = Chem.MolToInchiKey(m)
          if ik_full:
            anchor_ik_to_indices.setdefault(ik_full, []).append(idx)
        except Exception:
          pass
        if (boc_pat is not None and m.HasSubstructMatch(boc_pat)) or (
            cbz_pat is not None and m.HasSubstructMatch(cbz_pat)
        ):
          anchors_have_boc = True
        if aryl_br_i_pat is not None and m.HasSubstructMatch(aryl_br_i_pat):
          anchors_have_aryl_br_i = True
        try:
          sc = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
          if sc:
            anchor_scaffolds.add(sc)
          anchor_fps.append(fp_gen.GetFingerprint(m))
        except Exception:
          pass

  for child_pt, parent_pts in pat_parent_map.items():
    if any(parent_pt in pat_has_examples for parent_pt in parent_pts):
      pat_has_examples.add(child_pt)
  # When a patent family harvest includes both a granted US patent (with
  # extracted HTML Examples) and a table-less PCT/seed sibling (e.g.
  # WO2020231990A1 alongside US11780845B2), propagate pat_has_examples to
  # co-harvested Seed / DOCDB Family patents so unexemplified PubChem
  # associative rows from the WO sibling are pruned.
  if pat_has_examples and "patent_id" in df.columns:
    for _, row_p in df.iterrows():
      rel_p = str(row_p.get("patent_relationship") or "").lower()
      if (
          not rel_p
          or "seed" in rel_p
          or "docdb" in rel_p
          or "continuation" in rel_p
      ):
        for pat_tok in str(row_p.get("patent_id") or "").split(","):
          pt = pat_tok.strip().upper()
          if pt and pt not in ("NAN", "NONE"):
            pat_has_examples.add(pt)

  keep_indices = []
  assigned_series = []
  pruned_count = 0
  for idx, row in df.iterrows():
    raw_s = row.get(s_col)
    smi = (
        str(raw_s).strip()
        if pd.notna(raw_s)
        and str(raw_s).strip().lower() not in ("", "nan", "none")
        else ""
    )
    m = Chem.MolFromSmiles(smi) if smi else None
    if m is not None:
      if smarts_pat is not None and not m.HasSubstructMatch(smarts_pat):
        pruned_count += 1
        continue
      if scaf_pat is not None and not m.HasSubstructMatch(scaf_pat):
        pruned_count += 1
        continue

    src = str(row.get("source") or "")
    all_src = str(row.get("all_sources") or src)
    ex_id = row.get("patent_example_id")
    is_assoc = (
        src in ("PubChem", "Patent_Concept_OCSR")
        and "Patent_FullText_Table" not in all_src
        and (pd.isna(ex_id) or not str(ex_id).strip())
    )
    if is_assoc and has_unexemplified_assoc and anchor_fps:
      # When a patent (or its parent patent) already has extracted Examples in
      # the corpus, any remaining PubChem patent-association row with no
      # patent_example_id is not an Example of that patent. Backfill its
      # chembl_id onto any matching exemplified row and prune the associative row.
      row_pats = [
          pt.strip().upper()
          for pt in str(row.get("patent_id") or "").split(",")
          if pt.strip() and pt.strip().upper() not in ("NAN", "NONE")
      ]
      row_has_chembl = pd.notna(row.get("chembl_id")) and str(
          row.get("chembl_id")
      ).strip().lower() not in ("", "nan", "none")
      row_ik = ""
      if m is not None:
        try:
          row_ik = Chem.MolToInchiKey(m) or ""
        except Exception:
          pass
      if (
          row_has_chembl
          and row_ik
          and row_ik in anchor_ik_to_indices
          and "chembl_id" in df.columns
      ):
        if df["chembl_id"].dtype != object:
          df["chembl_id"] = df["chembl_id"].astype(object)
        for a_idx in anchor_ik_to_indices[row_ik]:
          cur_cid = df.at[a_idx, "chembl_id"]
          if pd.isna(cur_cid) or str(cur_cid).strip().lower() in (
              "",
              "nan",
              "none",
          ):
            df.at[a_idx, "chembl_id"] = row.get("chembl_id")
      if (
          src == "PubChem"
          and row_pats
          and all(pt in pat_has_examples for pt in row_pats)
      ):
        pruned_count += 1
        continue

      if m is not None:
        # Prune synthetic intermediates (boronic/silyl/stannyl, Boc/Cbz, aryl-Br/I)
        if synth_reagent_pat is not None and m.HasSubstructMatch(
            synth_reagent_pat
        ):
          pruned_count += 1
          continue
        if not anchors_have_boc and (
            (boc_pat is not None and m.HasSubstructMatch(boc_pat))
            or (cbz_pat is not None and m.HasSubstructMatch(cbz_pat))
        ):
          pruned_count += 1
          continue
        if (
            not anchors_have_aryl_br_i
            and aryl_br_i_pat is not None
            and m.HasSubstructMatch(aryl_br_i_pat)
        ):
          pruned_count += 1
          continue

      if m is not None and min_tanimoto > 0:
        try:
          sc = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
          if not (sc and sc in anchor_scaffolds):
            fp = fp_gen.GetFingerprint(m)
            sims = DataStructs.BulkTanimotoSimilarity(fp, anchor_fps)
            if not (sims and max(sims) >= min_tanimoto):
              pruned_count += 1
              continue
        except Exception:
          pass

    keep_indices.append(idx)
    if parsed_series:
      s_prev = row.get("series")
      mc_prev = row.get("Matched_Core")
      mc_str = ""
      if pd.notna(mc_prev):
        cand_mc = str(mc_prev).strip()
        if cand_mc.lower() not in ("", "nan", "none"):
          mc_str = cand_mc
      matched_lbl: Optional[str] = None
      for lbl, smt_str, pat in parsed_series:
        if mc_str and mc_str == smt_str:
          matched_lbl = lbl
          break
        if m is not None and pat is not None and m.HasSubstructMatch(pat):
          matched_lbl = lbl
          break
      if matched_lbl is None:
        if (
            pd.notna(s_prev)
            and str(s_prev).strip().lower() not in ("", "nan", "none")
            and not re.search(r"[=()\[\]1-9]", str(s_prev))
        ):
          matched_lbl = str(s_prev).strip()
        else:
          matched_lbl = "Other"
      assigned_series.append(matched_lbl)

  out_df = df.loc[keep_indices].copy().reset_index(drop=True)
  if parsed_series and len(assigned_series) == len(out_df):
    out_df["series"] = assigned_series
    out_df["Matched_Core"] = assigned_series
  if pruned_count > 0:
    print(
        "  [+] Applicability domain / structure filter (Murcko or Tanimoto >="
        f" {min_tanimoto:.2f}): retained {len(out_df)} compounds, pruned"
        f" {pruned_count} off-scaffold records.",
        flush=True,
    )
  return out_df


_PIVOTED_ENDPOINT_PREFIX_RE = re.compile(
    r"^(p?(?:Ki|IC50|Kd|EC50|Kis)|Inhibition|Ratio|Selectivity_ratio|AUC|fAUC|CL|Vdss|T1_2|Cmax|fCmax|Tmax|Cave|F|Fu|Stability|TGI)_"
)


def _augment_pactivity_columns(
    df: pd.DataFrame, min_endpoint_fraction: float = 0.05
) -> pd.DataFrame:
  """Automatically computes -log10(M) (pKi_*, pIC50_*, pKd_*, pEC50_*) columns for pivoted or coalesced nM potency columns, backfills missing pchembl_value, and folds sparse (<5%) pivoted endpoint columns into all_endpoints."""
  pat_prefix = re.compile(r"^(Ki|IC50|Kd|EC50)(?:_(.+))?$")
  pat_suffix = re.compile(r"^(.+?)_(Ki|IC50|Kd|EC50)(?:_nM)?$", re.IGNORECASE)
  new_cols: dict[str, pd.Series] = {}
  for col in list(df.columns):
    col_str = str(col)
    m_pre = pat_prefix.match(col_str)
    m_suf = pat_suffix.match(col_str) if not m_pre else None
    if not m_pre and not m_suf:
      continue
    if m_pre:
      prefix, suffix = m_pre.group(1), m_pre.group(2)
      if suffix and suffix.lower().endswith("_nm"):
        suffix = suffix[:-3]
      p_col = f"p{prefix}_{suffix}" if suffix else f"p{prefix}"
    else:
      target_stem, kind = m_suf.group(1), m_suf.group(2)
      kind_canon = {"ki": "Ki", "ic50": "IC50", "kd": "Kd", "ec50": "EC50"}.get(
          kind.lower(), kind
      )
      p_col = f"p{kind_canon}_{target_stem}"
    if p_col in df.columns:
      continue
    num_s = pd.to_numeric(df[col], errors="coerce")
    valid = num_s > 0
    if valid.any():
      p_vals = pd.Series(np.nan, index=df.index, dtype=float)
      p_vals.loc[valid] = np.round(9.0 - np.log10(num_s.loc[valid]), 3)
      new_cols[p_col] = p_vals
  if new_cols:
    df = pd.concat([df.copy(), pd.DataFrame(new_cols, index=df.index)], axis=1)

  # Backfill missing pchembl_value from primary pKi_*/pIC50_* or *_nM / patent potency columns
  if "pchembl_value" not in df.columns:
    df["pchembl_value"] = np.nan
  else:
    df["pchembl_value"] = pd.to_numeric(df["pchembl_value"], errors="coerce")

  p_act_cols = [
      c for c in df.columns if re.match(r"^p(Ki|IC50|Kd|EC50)(?:_.+)?$", str(c))
  ]
  pat_cols = [
      c
      for c in df.columns
      if (
          str(c).startswith(("patent_table_col_", "patent_assay_"))
          or str(c).endswith("_nM")
      )
      and not str(c).lower().endswith("_range")
  ]
  for idx_r, row_r in df[df["pchembl_value"].isna()].iterrows():
    filled = False
    for pc in p_act_cols:
      val_p = pd.to_numeric(row_r.get(pc), errors="coerce")
      if pd.notna(val_p):
        df.at[idx_r, "pchembl_value"] = round(float(val_p), 3)
        filled = True
        break
    if not filled:
      for ptc in pat_cols:
        sig = _cell_potency_sig(row_r.get(ptc), str(ptc))
        if sig is not None and sig[0] is not None:
          df.at[idx_r, "pchembl_value"] = sig[0]
          break

  # Prune sparse pivoted endpoint columns (< min_endpoint_fraction non-null rows)
  # after preserving their values in the per-row 'all_endpoints' JSON column.
  if len(df) >= 20 and min_endpoint_fraction > 0:
    min_req = max(5, int(np.ceil(min_endpoint_fraction * len(df))))
    sparse_cols = [
        c
        for c in df.columns
        if _PIVOTED_ENDPOINT_PREFIX_RE.match(str(c))
        and int(df[c].notna().sum()) < min_req
    ]
    if sparse_cols:
      if "all_endpoints" not in df.columns:
        df["all_endpoints"] = None
      raw_sparse = [
          c
          for c in sparse_cols
          if not str(c).startswith(("pKi_", "pIC50_", "pKd_", "pEC50_"))
      ]
      if raw_sparse:
        updated_ep: list[Optional[str]] = []
        for _, r in df.iterrows():
          cur_ep = r.get("all_endpoints")
          ep_dict: dict[str, Any] = {}
          if pd.notna(cur_ep) and str(cur_ep).strip():
            try:
              parsed = (
                  json.loads(str(cur_ep)) if isinstance(cur_ep, str) else cur_ep
              )
              if isinstance(parsed, dict):
                ep_dict.update(parsed)
            except Exception:
              pass
          for sc in raw_sparse:
            v_sc = r.get(sc)
            if pd.notna(v_sc) and sc not in ep_dict:
              ep_dict[sc] = v_sc
          updated_ep.append(json.dumps(ep_dict) if ep_dict else None)
        df["all_endpoints"] = updated_ep
      df = df.drop(columns=sparse_cols, errors="ignore")
  return df


def _print_discovered_endpoint_columns_hint(
    df: pd.DataFrame, output_csv: str
) -> None:
  """Prints discovered bioactivity/assay columns and a CLI hint for coalescing columns via --col-map without custom scripts."""
  endpoint_cols: list[tuple[str, int]] = []
  for c in df.columns:
    cs = str(c)
    if cs.endswith("_range") or cs == "activity_type":
      continue
    if (
        _PIVOTED_ENDPOINT_PREFIX_RE.match(cs)
        or cs.startswith(("patent_assay_", "patent_table_col_"))
        or cs in ("activity_bin", "activity_value")
        or cs.endswith("_nM")
        or re.search(
            r"(?:^|_)(Ki|IC50|Kd|EC50|Ratio)(?:_|$)", cs, re.IGNORECASE
        )
    ):
      cnt = int(pd.to_numeric(df[c], errors="coerce").notna().sum())
      if cnt == 0:
        disc_s = (
            df[c]
            .dropna()
            .astype(str)
            .str.strip()
            .loc[
                lambda s: s.str.match(
                    r"^(?:[A-F](?:\s*\(.*\))?|\+{1,6}|\*{1,6}|[<>≤≥~].+)$",
                    case=False,
                )
            ]
        )
        cnt = int(len(disc_s))
      if cnt > 0:
        endpoint_cols.append((cs, cnt))
  if not endpoint_cols:
    return
  summary_items = ", ".join(f"{c} (n={cnt})" for c, cnt in endpoint_cols[:10])
  print(
      f"  [+] Active Bioactivity / Assay Columns: {summary_items}", flush=True
  )
  if len(endpoint_cols) >= 2:
    c0, c1 = endpoint_cols[0][0], endpoint_cols[1][0]
    print(
        "  [TIP] To coalesce heterogeneous endpoint columns across sources"
        " or export split subsets, re-run curate_congeneric_corpus.py"
        " directly (do NOT write a custom Python script):\n"
        f"        uv run {os.path.abspath(__file__)} --plot-only --input-csv"
        f" {output_csv} --output-csv {output_csv} --col-map"
        f" 'Unified_Ki_nM={c0}|{c1}' --split-by series",
        flush=True,
    )


def _export_split_subsets(
    df: pd.DataFrame,
    mols: list[Optional[Chem.Mol]],
    split_by: str,
    output_csv: str,
    output_sdf: Optional[str] = None,
) -> None:
  """Exports per-group CSV and SDF files for each cohort in split_by."""
  if split_by not in df.columns:
    print(
        f"  [!] Warning: --split-by column '{split_by}' not found in dataframe."
    )
    return

  base_csv_dir = os.path.dirname(os.path.abspath(output_csv)) or "."
  base_csv_stem = os.path.splitext(os.path.basename(output_csv))[0]
  base_sdf_dir = (
      os.path.dirname(os.path.abspath(output_sdf))
      if output_sdf
      else base_csv_dir
  )
  base_sdf_stem = (
      os.path.splitext(os.path.basename(output_sdf))[0]
      if output_sdf
      else base_csv_stem
  )

  for grp_val, sub_df in df.groupby(
      df[split_by].fillna("unlabelled"), sort=False
  ):
    if sub_df.empty:
      continue
    safe_tag = (
        re.sub(r"[^A-Za-z0-9_.-]+", "_", str(grp_val)).strip("_") or "group"
    )
    sub_csv_path = os.path.join(base_csv_dir, f"{base_csv_stem}_{safe_tag}.csv")
    sub_df.to_csv(sub_csv_path, index=False)
    print(
        f"  [+] Saved split subset CSV ({split_by}={grp_val}, N={len(sub_df)})"
        f" to: {sub_csv_path}"
    )
    if output_sdf:
      sub_sdf_path = os.path.join(
          base_sdf_dir, f"{base_sdf_stem}_{safe_tag}.sdf"
      )
      with Chem.SDWriter(sub_sdf_path) as writer:
        for idx_row in sub_df.index:
          mol = mols[idx_row] if idx_row < len(mols) else None
          if mol is not None:
            try:
              compute_clean_2d_coords(mol)
              row_data = df.loc[idx_row]
              for k, v in row_data.items():
                if pd.notna(v):
                  mol.SetProp(str(k), str(v))
              writer.write(mol)
            except Exception:
              pass
      print(
          f"  [+] Saved split subset SDF ({split_by}={grp_val},"
          f" N={len(sub_df)}) to: {sub_sdf_path}"
      )


def _export_html_and_structures_png(
    master_df: pd.DataFrame,
    s_col: Optional[str],
    group_by: Optional[str] = None,
    plot_title: Optional[str] = None,
    output_html: Optional[str] = None,
    output_structures_png: Optional[str] = None,
    highlight_compounds: Optional[Sequence[str]] = None,
) -> None:
  """Exports the interactive 2D HTML compound viewer and representative structures PNG grid."""
  if not s_col:
    return
  if output_html:
    try:
      create_molecular_grid(
          master_df,
          smiles_col=s_col,
          filter_active=False,
          max_cards=len(master_df) if len(master_df) <= 1000 else 500,
          output_html=output_html,
          title=plot_title or "Interactive Congeneric Structure & SAR Explorer",
      )
      print(f"  [+] Saved interactive HTML structure viewer to: {output_html}")
    except Exception as e:
      print(f"  [!] Warning: Could not generate interactive HTML viewer: {e}")

  if output_structures_png:
    try:
      _, grp_list = _resolve_grouping(master_df, group_by)
      rep_smiles: list[str] = []
      rep_legends: list[str] = []
      per_grp = max(2, 6 // max(len(grp_list), 1))
      has_dense_pchembl = (
          "pchembl_value" in master_df.columns
          and int(
              pd.to_numeric(master_df["pchembl_value"], errors="coerce")
              .notna()
              .sum()
          )
          >= 2
      )
      hl_tokens = [
          re.sub(r"[\s_-]+", "_", str(t).strip().lower()).strip("_")
          for t in highlight_compounds or []
          if str(t).strip()
      ]
      bin_order_map = {
          "A": 0,
          "B": 1,
          "C": 2,
          "D": 3,
          "E": 4,
          "++++": 0,
          "+++": 1,
          "++": 2,
          "+": 3,
      }
      for g_label, g_mask, _ in grp_list:
        g_df = master_df[g_mask].copy()
        if hl_tokens:

          def _is_hl_row(row_s: pd.Series) -> int:
            for c_name in (
                "compound_name",
                "pref_name",
                "compound_id",
                "paper_compound_id",
                "patent_example_id",
                "chembl_id",
                "Compound_Name",
            ):
              raw_v = str(row_s.get(c_name, "")).strip()
              if not raw_v or raw_v.lower() in ("nan", "none"):
                continue
              val_norm = re.sub(r"[\s_-]+", "_", raw_v.lower()).strip("_")
              if any(
                  tok == val_norm
                  or bool(
                      re.search(
                          rf"(?:^|[:,;(\s]){re.escape(tok)}(?:$|[,;)\s])",
                          val_norm,
                      )
                  )
                  for tok in hl_tokens
              ):
                return 0
            return 1

          g_df["_sort_hl"] = g_df.apply(_is_hl_row, axis=1)
        else:
          g_df["_sort_hl"] = 1

        if has_dense_pchembl and "pchembl_value" in g_df.columns:
          g_df["_sort_p"] = -pd.to_numeric(
              g_df["pchembl_value"], errors="coerce"
          ).fillna(-1.0)
        else:
          g_df["_sort_p"] = 0.0

        if "activity_bin" in g_df.columns:
          g_df["_sort_bin"] = (
              g_df["activity_bin"]
              .astype(str)
              .str.strip()
              .map(lambda b: bin_order_map.get(b, 99))
          )
        else:
          g_df["_sort_bin"] = 99

        if "patent_example_id" in g_df.columns:

          def _ex_num(val_e: object) -> int:
            m_ex = re.search(r"\d+", str(val_e))
            return int(m_ex.group(0)) if m_ex else 999999

          g_df["_sort_ex"] = g_df["patent_example_id"].map(_ex_num)
        else:
          g_df["_sort_ex"] = 999999

        g_df = g_df.sort_values(
            ["_sort_hl", "_sort_p", "_sort_bin", "_sort_ex"],
            ascending=[True, True, True, True],
        )
        n_hl_in_grp = int((g_df["_sort_hl"] == 0).sum())
        take_n = max(per_grp, n_hl_in_grp)
        for _, r in g_df.head(take_n).iterrows():
          smi_val = str(r.get(s_col, "")).strip()
          if smi_val and smi_val not in rep_smiles:
            raw_cid = str(
                r.get("compound_id") or r.get("patent_example_id") or "Compound"
            ).strip()
            raw_cname = str(
                r.get("compound_name") or r.get("pref_name") or ""
            ).strip()
            raw_ex = str(r.get("patent_example_id") or "").split(",")[0].strip()
            if (
                raw_cname
                and raw_cname.lower() not in ("nan", "none")
                and re.match(r"^[A-Za-z][A-Za-z0-9 -]{3,27}$", raw_cname)
                and not re.match(
                    r"^(?:CHEMBL\d+|CID_|US\d+|WO\d+|EP\d+|CN\d+|JP\d+|Example|Compound)",
                    raw_cname,
                    re.I,
                )
            ):
              cid = raw_cname
            elif (
                re.match(r"^(?:CID_\d+|CHEMBL\d+)$", raw_cid, re.I)
                and raw_ex
                and raw_ex.lower() not in ("nan", "none")
            ):
              cid = raw_ex
            else:
              cid = raw_cid.split(",")[0].strip()
            p_val = r.get("pchembl_value")
            p_txt = ""
            if pd.notna(p_val) and has_dense_pchembl:
              p_txt = f" | pAct={float(p_val):.2f}"
            elif pd.notna(r.get("activity_bin")) and str(
                r.get("activity_bin")
            ).strip() not in ("", "nan", "None"):
              b_str = str(r.get("activity_bin")).strip()
              r_str = str(r.get("activity_range", "")).strip()
              if r_str and r_str not in ("", "nan", "None"):
                p_txt = f" | Bin {b_str} ({r_str})"
              else:
                p_txt = f" | Bin {b_str}"
            rep_smiles.append(smi_val)
            rep_legends.append(f"{cid} ({g_label}){p_txt}")
      if rep_smiles:
        render_grid_molecules(
            rep_smiles[:6],
            output_structures_png,
            legends=rep_legends[:6],
            mols_per_row=3,
            sub_img_size=(400, 340),
            transparent=False,
        )
        print(
            "  [+] Rendered opaque-white 2D representative structures grid to:"
            f" {os.path.abspath(output_structures_png)}"
        )
    except Exception as e:
      print(
          f"  [!] Warning: Could not render representative structures PNG: {e}"
      )


def curate_congeneric_corpus(
    input_csvs: Optional[Sequence[str]] = None,
    chembl_docs: Optional[Sequence[str]] = None,
    output_csv: str = "master_corpus.csv",
    output_sdf: Optional[str] = None,
    output_plot: Optional[str] = None,
    output_html: Optional[str] = None,
    output_structures_png: Optional[str] = None,
    smiles_col: Optional[str] = None,
    group_by: Optional[str] = None,
    split_by: Optional[str] = None,
    panels: Optional[Sequence[Union[str, dict[str, object]]]] = None,
    plot_title: Optional[str] = None,
    plot_only: bool = False,
    min_tanimoto: float = 0.40,
    query_smarts: Optional[str] = None,
    query_scaffold: Optional[str] = None,
    series_smarts: Optional[str] = None,
    highlight_compounds: Optional[Sequence[str]] = None,
    rename_cols: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
  """Unified execution pipeline: ingests, featurizes, annotates, and plots congeneric libraries."""
  if panels and not output_plot:
    output_plot = "chemical_space_distributions.png"

  dfs = []

  if chembl_docs and not plot_only:
    chembl_df = fetch_chembl_document_activities(chembl_docs)
    if not chembl_df.empty:
      chembl_df = _disambiguate_and_rename_input_columns(
          chembl_df, source_path="chembl", rename_cols=rename_cols
      )
      dfs.append(chembl_df)

  if input_csvs:
    expanded_paths: list[str] = []
    for raw_item in input_csvs:
      for part in str(raw_item).split(","):
        if part.strip():
          expanded_paths.append(part.strip())
    for path_item in expanded_paths:
      dfs.extend(_load_any_chemical_input(path_item, rename_cols=rename_cols))

  if not dfs:
    raise ValueError(
        "No input data sources provided (provide --inputs / --input-csv or"
        " --chembl-docs)."
    )

  master_df = pd.concat(dfs, ignore_index=True)

  if "canonical_smiles" in master_df.columns and "SMILES" in master_df.columns:
    master_df["canonical_smiles"] = master_df["canonical_smiles"].fillna(
        master_df["SMILES"]
    )
    master_df["SMILES"] = master_df["SMILES"].fillna(
        master_df["canonical_smiles"]
    )
  elif (
      "SMILES" in master_df.columns
      and "canonical_smiles" not in master_df.columns
  ):
    master_df["canonical_smiles"] = master_df["SMILES"]
  elif (
      "canonical_smiles" in master_df.columns
      and "SMILES" not in master_df.columns
  ):
    master_df["SMILES"] = master_df["canonical_smiles"]

  s_col = None
  for cand in [
      smiles_col,
      "canonical_smiles",
      "SMILES",
      "smiles",
      "Canonical_SMILES",
  ]:
    if cand and cand in master_df.columns:
      s_col = cand
      break

  if plot_only:
    mols_plot: list[Optional[Chem.Mol]] = []
    if (
        "MW" not in master_df.columns or "cLogP" not in master_df.columns
    ) and s_col:
      print(
          "  [+] Descriptor columns missing in --plot-only mode; featurizing"
          f" {len(master_df)} rows...",
          flush=True,
      )
      master_df, mols_plot = _featurize_dataframe(master_df, s_col)
    elif s_col and (output_sdf or split_by):
      mols_plot = [
          Chem.MolFromSmiles(str(s)) if pd.notna(s) else None
          for s in master_df[s_col]
      ]
    master_df = _backfill_document_dates(master_df)
    master_df = _augment_pactivity_columns(master_df)
    eff_csv_plot = (
        output_csv
        if (output_csv != "master_corpus.csv" or rename_cols)
        else (input_csvs[0] if input_csvs else output_csv)
    )
    if rename_cols or (
        output_csv != "master_corpus.csv"
        and input_csvs
        and len(input_csvs) == 1
        and os.path.abspath(output_csv) != os.path.abspath(input_csvs[0])
    ):
      os.makedirs(
          os.path.dirname(os.path.abspath(eff_csv_plot)) or ".", exist_ok=True
      )
      master_df.to_csv(eff_csv_plot, index=False)
      print(f"  [+] Updated master congeneric corpus CSV: {eff_csv_plot}")
    if output_sdf and mols_plot:
      os.makedirs(
          os.path.dirname(os.path.abspath(output_sdf)) or ".", exist_ok=True
      )
      with Chem.SDWriter(output_sdf) as writer:
        for mol, (_, row) in zip(mols_plot, master_df.iterrows()):
          if mol is not None:
            try:
              compute_clean_2d_coords(mol)
              for k, v in row.items():
                if pd.notna(v):
                  mol.SetProp(str(k), str(v))
              writer.write(mol)
            except (RuntimeError, ValueError, TypeError, KeyError) as exc:
              print(
                  f"  [!] Warning: Skipped molecule in SDF export: {exc}",
                  file=sys.stderr,
              )
      print(f"  [+] Updated master SDF library: {output_sdf}")
    if split_by:
      _export_split_subsets(
          master_df, mols_plot, split_by, eff_csv_plot, output_sdf
      )
    _export_html_and_structures_png(
        master_df,
        s_col=s_col,
        group_by=group_by,
        plot_title=plot_title,
        output_html=output_html,
        output_structures_png=output_structures_png,
        highlight_compounds=highlight_compounds,
    )
    _print_discovered_endpoint_columns_hint(master_df, eff_csv_plot)
    _print_cohort_summary_table(master_df, group_by, panels=panels)
    if (
        output_plot
        or panels
        or (
            not rename_cols
            and not split_by
            and not output_sdf
            and not output_html
            and not output_structures_png
        )
    ):
      plot_chemical_space_distributions(
          master_df,
          output_plot or "chemical_space_distributions.png",
          title=plot_title
          or "Congeneric Chemical Space & Property Distributions",
          group_by=group_by,
          panels=panels,
          highlight_compounds=highlight_compounds,
          output_csv=eff_csv_plot,
      )
    return master_df

  if not s_col:
    raise ValueError("Could not find valid SMILES column in merged datasets.")

  master_df = _apply_corpus_structural_filters_and_series(
      master_df,
      s_col=s_col,
      min_tanimoto=min_tanimoto,
      query_smarts=query_smarts,
      query_scaffold=query_scaffold,
      series_smarts=series_smarts,
  )

  print(
      f"  [+] Deduplicating {len(master_df)} raw entries by SMILES string...",
      flush=True,
  )
  raw_df = master_df
  master_df = deduplicate_by_smiles(raw_df, s_col)
  if "series" in master_df.columns and "Matched_Core" in master_df.columns:
    master_df["series"] = (
        master_df["series"]
        .replace({"": np.nan, "nan": np.nan, "None": np.nan})
        .fillna(master_df["Matched_Core"])
        .fillna("Other")
    )
    master_df["Matched_Core"] = (
        master_df["Matched_Core"]
        .replace({"": np.nan, "nan": np.nan, "None": np.nan})
        .fillna(master_df["series"])
    )
  elif "series" in master_df.columns:
    master_df["series"] = (
        master_df["series"]
        .replace({"": np.nan, "nan": np.nan, "None": np.nan})
        .fillna("Other")
    )
  master_df = _disambiguate_and_rename_input_columns(master_df)
  _print_dedup_summary(raw_df, master_df, s_col)
  master_df = master_df.reset_index(drop=True)

  print(
      f"  [+] Featurizing {len(master_df)} retained rows with canonical"
      " RDKit descriptors and CNS MPO...",
      flush=True,
  )
  master_df, mols = _featurize_dataframe(master_df, s_col)
  master_df = _augment_pactivity_columns(master_df)
  if "canonical_smiles" in master_df.columns:
    master_df["SMILES"] = master_df["canonical_smiles"]
    if "Clean_SMILES" in master_df.columns:
      master_df["Clean_SMILES"] = master_df["canonical_smiles"]

  os.makedirs(
      os.path.dirname(os.path.abspath(output_csv)) or ".", exist_ok=True
  )
  master_df.to_csv(output_csv, index=False)
  print(f"  [+] Saved master congeneric corpus CSV to: {output_csv}")

  if output_sdf:
    os.makedirs(
        os.path.dirname(os.path.abspath(output_sdf)) or ".", exist_ok=True
    )
    with Chem.SDWriter(output_sdf) as writer:
      for idx, (mol, (_, row)) in enumerate(zip(mols, master_df.iterrows())):
        if mol is not None:
          try:
            compute_clean_2d_coords(mol)
            for k, v in row.items():
              if pd.notna(v):
                mol.SetProp(str(k), str(v))
            writer.write(mol)
          except (RuntimeError, ValueError, TypeError, KeyError) as exc:
            print(
                f"  [!] Warning: Skipped molecule {idx} in SDF export: {exc}",
                file=sys.stderr,
            )
    print(f"  [+] Saved master SDF library to: {output_sdf}")

  if split_by:
    _export_split_subsets(master_df, mols, split_by, output_csv, output_sdf)

  _export_html_and_structures_png(
      master_df,
      s_col=s_col,
      group_by=group_by,
      plot_title=plot_title,
      output_html=output_html,
      output_structures_png=output_structures_png,
      highlight_compounds=highlight_compounds,
  )

  _print_discovered_endpoint_columns_hint(master_df, output_csv)
  _print_cohort_summary_table(master_df, group_by, panels=panels)

  if output_plot:
    plot_chemical_space_distributions(
        master_df,
        output_plot,
        title=plot_title
        or "Congeneric Chemical Space & Property Distributions",
        group_by=group_by,
        panels=panels,
        highlight_compounds=highlight_compounds,
        output_csv=output_csv,
    )

  return master_df


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Congeneric Series & Literature Corpus Aggregation Engine:"
          " concatenates inputs, deduplicates by SMILES string, computes RDKit"
          " descriptors, and renders configurable multi-panel property"
          " distributions."
      )
  )
  parser.add_argument(
      "--chembl-docs",
      "--chembl_docs",
      "--pmids",
      "--pmid",
      dest="chembl_docs",
      nargs="*",
      default=None,
      help=(
          "One or more ChEMBL document IDs (CHEMBL...), PubMed IDs"
          " (PMID:<ID> or <PMID>), or DOIs (10.1021/...) to harvest"
          " with built-in EBI retry and PubChem BioAssay fallback."
      ),
  )
  parser.add_argument(
      "--inputs",
      "--input-csv",
      "--input_csv",
      "--input",
      dest="input_csv",
      nargs="*",
      default=None,
      help=(
          "Input CSV, JSON, SDF files or directories (e.g. patent tables,"
          " ChEMBL JSONs, SDFs)."
      ),
  )
  parser.add_argument(
      "--output-csv",
      "--output_csv",
      "--output",
      dest="output_csv",
      default="master_corpus.csv",
      help="Destination path for curated master CSV.",
  )
  parser.add_argument(
      "--output-sdf",
      "--output_sdf",
      dest="output_sdf",
      default=None,
      help="Destination path for curated master SDF.",
  )
  parser.add_argument(
      "--output-plot",
      "--output_plot",
      dest="output_plot",
      default=None,
      help=(
          "Path for saving chemical space distribution graphic (requires"
          " --panels or --panels-json)."
      ),
  )
  parser.add_argument(
      "--output-html",
      "--output_html",
      dest="output_html",
      default=None,
      help="Optional path to auto-generate interactive 2D HTML card viewer.",
  )
  parser.add_argument(
      "--output-structures-png",
      "--output_structures_png",
      dest="output_structures_png",
      default=None,
      help=(
          "Optional path to auto-render an opaque-white 2D representative"
          " structures grid PNG."
      ),
  )
  parser.add_argument(
      "--min-tanimoto",
      "--min_tanimoto",
      dest="min_tanimoto",
      type=float,
      default=0.40,
      help=(
          "Minimum Morgan ECFP4 Tanimoto similarity to exemplified/ChEMBL"
          " compounds (unless sharing a Murcko scaffold) for un-exemplified"
          " associative records (default: 0.40)."
      ),
  )
  parser.add_argument(
      "--query-smarts",
      "--query_smarts",
      dest="query_smarts",
      default=None,
      help="Optional SMARTS pattern to filter the merged corpus.",
  )
  parser.add_argument(
      "--query-scaffold",
      "--query_scaffold",
      dest="query_scaffold",
      default=None,
      help="Optional SMILES/SMARTS scaffold to filter the merged corpus.",
  )
  parser.add_argument(
      "--series-smarts",
      "--series_smarts",
      dest="series_smarts",
      default=None,
      help=(
          "Optional comma-separated 'Label=SMARTS' definitions to classify"
          " compounds into series cohorts (e.g. 'SeriesA=...,SeriesB=...')."
      ),
  )
  parser.add_argument(
      "--smiles-col",
      "--smiles_col",
      dest="smiles_col",
      default=None,
      help="Column name containing ligand SMILES.",
  )
  parser.add_argument(
      "--group-by",
      "--group_by",
      dest="group_by",
      default=None,
      help=(
          "Column name to group/color plots by (e.g. series, source,"
          " Matched_Core). Auto-detected if omitted; pass 'none' for"
          " single-series plotting."
      ),
  )
  parser.add_argument(
      "--split-by",
      "--split_by",
      dest="split_by",
      default=None,
      help=(
          "Optional column name to export per-group subset CSV and SDF files"
          " (e.g. series, source, Matched_Core)."
      ),
  )
  parser.add_argument(
      "--panels",
      nargs="*",
      default=None,
      help=(
          "Required when generating --output-plot. Custom panel specs using"
          " 'type:feature[:feature2][:title]' shorthand. Supported types:"
          " scatter, hist, box, kde, bar. Example: --panels"
          " 'scatter:MW:cLogP' 'box:SA_Score' 'kde:TPSA' 'bar:potency'"
      ),
  )
  parser.add_argument(
      "--panels-json",
      "--panels_json",
      dest="panels_json",
      default=None,
      help="JSON list of panel configuration dicts (overrides --panels).",
  )
  parser.add_argument(
      "--highlight",
      "--highlight-compounds",
      "--highlight_compounds",
      dest="highlight_compounds",
      nargs="*",
      default=None,
      help=(
          "Optional compound IDs, names, or paper_compound_id labels (comma-"
          " or space-separated) to highlight with labeled markers on scatter"
          " panels."
      ),
  )
  parser.add_argument(
      "--plot-title",
      "--plot_title",
      dest="plot_title",
      default=None,
      help="Custom main figure title (suptitle).",
  )
  parser.add_argument(
      "--plot-only",
      "--plot_only",
      dest="plot_only",
      action="store_true",
      help=(
          "Re-plot an existing curated CSV without re-running ChEMBL"
          " pagination, SMILES deduplication, or CSV/SDF export."
      ),
  )
  parser.add_argument(
      "--rename-cols",
      "--rename_cols",
      "--col-map",
      "--col_map",
      dest="rename_cols",
      nargs="*",
      default=None,
      help=(
          "Optional column rename/coalesce rules on input files"
          " ('new=old1|old2', 'new=old', 'raw_patent_col=new', or"
          " 'file.csv:old=new') to align per-file endpoint columns before"
          " merging."
      ),
  )

  args = parser.parse_args()

  panels_arg: Optional[Sequence[Union[str, dict[str, object]]]] = args.panels
  if args.panels_json:
    panels_arg = json.loads(args.panels_json)

  hl_arg: list[str] = []
  if args.highlight_compounds:
    for item in args.highlight_compounds:
      for part in str(item).split(","):
        if part.strip():
          hl_arg.append(part.strip())

  try:
    curate_congeneric_corpus(
        input_csvs=args.input_csv,
        chembl_docs=args.chembl_docs,
        output_csv=args.output_csv,
        output_sdf=args.output_sdf,
        output_plot=args.output_plot,
        output_html=args.output_html,
        output_structures_png=args.output_structures_png,
        smiles_col=args.smiles_col,
        group_by=args.group_by,
        split_by=args.split_by,
        panels=panels_arg,
        plot_title=args.plot_title,
        plot_only=args.plot_only,
        min_tanimoto=args.min_tanimoto,
        query_smarts=args.query_smarts,
        query_scaffold=args.query_scaffold,
        series_smarts=args.series_smarts,
        highlight_compounds=hl_arg if hl_arg else None,
        rename_cols=args.rename_cols,
    )
  except ValueError as exc:
    if "Plotting requires explicit --panels" in str(exc):
      print(f"{exc}", file=sys.stderr)
      sys.exit(2)
    raise


if __name__ == "__main__":
  main()
