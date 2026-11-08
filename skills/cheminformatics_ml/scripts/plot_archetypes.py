# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "matplotlib>=3.7.0",
#     "seaborn>=0.12.0",
#     "pandas>=2.0.0",
#     "numpy>=1.24.0",
#     "scipy>=1.10.0",
#     "scikit-learn>=1.2.0",
# ]
# ///
"""Reusable multi-panel plotting routines with strict anti-overlap and complete evidence formatting."""

import argparse
import os
import re
import sys
import textwrap
from typing import Any, List, Optional, Sequence, Tuple, Union
import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, r2_score

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

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

from compat import (
    PUBLICATION_CYCLE,
    PUBLICATION_PALETTE,
    apply_publication_spines,
    set_safe_plot_style,
)


def short_plot_label(raw_id: Any, max_len: int = 16) -> str:
  """Condenses a compound identifier to a short, non-overlapping plot label.

  Strips parenthetical secondary codes from trivial names (e.g.
  '<LeadName> (<CodeID>)' -> '<LeadName>', 'Cmpd 1a (CHEMBL...)' ->
  'Cmpd 1a') and replaces or truncates long IUPAC/systematic names so point
  annotations and bar tick labels never overflow plot boundaries.
  """
  if raw_id is None or (isinstance(raw_id, float) and np.isnan(raw_id)):
    return "Cmpd"
  s = str(raw_id).strip()
  if not s or s.lower() in ("nan", "none", "null"):
    return "Cmpd"

  # 1. Trivial prefix with parenthetical code: '<LeadName> (<CodeID>)' or 'Cmpd 1a (CHEMBL...)'
  m_triv = re.match(r"^([A-Za-z][A-Za-z0-9_.\s-]{1,15})\s*\(([^)]+)\)$", s)
  if m_triv:
    prefix, paren = m_triv.group(1).strip(), m_triv.group(2).strip()
    if not re.search(r"\b(rac|rel|cis|trans)\b", prefix, re.I):
      if prefix.upper().startswith("CHEMBL") and len(paren) <= 10:
        return paren if paren.lower().startswith("cmpd") else f"Cmpd {paren}"
      return prefix

  # 2. Check for embedded short identifier inside a long string
  if len(s) > max_len:
    m_ex = re.search(
        r"\b(Example[_\s-]?\d+[a-z]?|Cmpd[_\s-]?[0-9]+[a-z]?|CHEMBL\d+)\b",
        s,
        re.I,
    )
    if m_ex:
      return m_ex.group(1).strip()
    return s[: max_len - 1].rstrip(" -_,([") + "…"
  return s


def _short_metric_label(col_name: str) -> str:
  """Shortens verbose Delta_* column names for inline plot callout boxes."""
  c = str(col_name).strip()
  replacements = {
      "Delta_pActivity": "ΔpAct",
      "Delta_Activity": "ΔAct",
      "Delta_CNS_MPO": "ΔMPO",
      "Delta_cLogP": "ΔcLogP",
      "Delta_TPSA": "ΔTPSA",
      "Delta_RotB": "ΔRotB",
      "Delta_MW": "ΔMW",
      "Delta_LipE": "ΔLipE",
  }
  if c in replacements:
    return replacements[c]
  if c.startswith("Delta_"):
    return "Δ" + c[6:14]
  return c[:12]


def _resolve_metric_column(
    df: pd.DataFrame,
    explicit: Optional[str],
    candidates: list[str],
    what: str,
) -> str:
  """Resolves which benchmark metric column to plot.

  `congeneric_series_ml.py` emits two different column schemas depending on
  whether a held-out test split was supplied, so a hard-coded column name
  silently works for one schema and raises `KeyError` for the other. Pass an
  explicit name to override; otherwise the first candidate actually present in
  the frame wins.

  Raises:
    KeyError: if none of the candidates are present, with the frame's real
      columns listed. Failing loudly is deliberate: plotting the wrong metric
      under a correct-looking axis label is worse than not plotting at all.
  """
  if explicit is not None:
    if explicit not in df.columns:
      raise KeyError(
          f"Requested {what} column {explicit!r} is not in the benchmark frame."
          f" Available columns: {sorted(df.columns)}"
      )
    return explicit
  for cand in candidates:
    if cand in df.columns:
      return cand
  raise KeyError(
      f"No recognised {what} column found in the benchmark frame. Looked for"
      f" {candidates}; available columns: {sorted(df.columns)}"
  )


def plot_directional_barh(
    ax: plt.Axes,
    categories: Sequence[str],
    values: Sequence[float],
    colors: Optional[Sequence[str]] = None,
    fmt: str = "%+.2f",
    pad_ratio: float = 0.35,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
) -> Any:
  """Plots a horizontal bar chart with guaranteed anti-collision text placement for negative & positive values.

  Solves the classic Matplotlib bug where negative bar labels collide with
  y-axis category tick labels.
  Enforces:
  1. Symmetric auto-padded axis limits (-max_val * (1 + pad_ratio), +max_val *
  (1 + pad_ratio)).
  2. Directional label placement (left of bar tip for negative, right for
  positive).
  3. Explicit y-axis tick padding (pad=8) to ensure category names never touch
  bar tips.
  """
  vals = np.array(values, dtype=float)
  cats = list(categories)
  max_val = (
      float(np.max(np.abs(vals)))
      if len(vals) > 0 and np.max(np.abs(vals)) > 0
      else 1.0
  )

  ax.set_xlim(-max_val * (1.0 + pad_ratio), max_val * (1.0 + pad_ratio))
  ax.axvline(0, color="#64748B", linestyle="--", linewidth=0.8)

  if colors is None:
    colors = [
        PUBLICATION_PALETTE["success"]
        if v >= 0
        else PUBLICATION_PALETTE["danger"]
        for v in vals
    ]

  bars = ax.barh(
      cats,
      vals,
      color=colors,
      edgecolor="#0F172A",
      height=0.65,
      linewidth=0.8,
      alpha=0.9,
  )

  for b in bars:
    w = b.get_width()
    offset = max_val * 0.03
    if w >= 0:
      ax.text(
          w + offset,
          b.get_y() + b.get_height() / 2,
          fmt % w,
          va="center",
          ha="left",
          fontsize=8.5,
          fontweight="bold",
      )
    else:
      ax.text(
          w - offset,
          b.get_y() + b.get_height() / 2,
          fmt % w,
          va="center",
          ha="right",
          fontsize=8.5,
          fontweight="bold",
      )

  ax.tick_params(axis="y", pad=8)
  if title:
    ax.set_title(title, fontweight="bold", pad=10)
  if xlabel:
    ax.set_xlabel(xlabel)
  ax.grid(axis="x", color="#F1F5F9", linewidth=0.8)
  apply_publication_spines(ax)
  return bars


def plot_screening_attrition_and_properties(
    stages: List[str],
    counts: List[int],
    deck_df: pd.DataFrame,
    hits_df: pd.DataFrame,
    output_png: str,
    property_col: str = "MW",
    qed_col: str = "QED",
    mpo_col: str = "CNS_MPO",
) -> str:
  """Archetype 1: Virtual screening funnel, property shift, MPO distribution."""
  abs_out = os.path.abspath(output_png)
  os.makedirs(os.path.dirname(abs_out), exist_ok=True)

  fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
      2, 2, figsize=(14, 10), dpi=300, layout="constrained"
  )
  set_safe_plot_style(font_scale=1.05)

  # 1. Attrition Funnel
  ax1.barh(
      stages[::-1],
      counts[::-1],
      color=PUBLICATION_PALETTE["primary"],
      edgecolor="#0F172A",
      height=0.6,
  )
  ax1.set_title("A. Multi-Tier Screening Attrition", fontweight="bold", pad=10)
  ax1.set_xlabel("Compound Count")
  apply_publication_spines(ax1)

  # 2. Property Distribution
  if property_col in deck_df.columns and property_col in hits_df.columns:
    sns.kdeplot(
        data=deck_df[property_col].dropna(),
        ax=ax2,
        label=f"Deck (N={len(deck_df)})",
        color=PUBLICATION_PALETTE["muted"],
        fill=True,
        alpha=0.35,
    )
    sns.kdeplot(
        data=hits_df[property_col].dropna(),
        ax=ax2,
        label=f"Hits (N={len(hits_df)})",
        color=PUBLICATION_PALETTE["success"],
        fill=True,
        alpha=0.5,
    )
    ax2.set_title(
        f"B. {property_col} Distribution Shift", fontweight="bold", pad=10
    )
    ax2.set_xlabel(property_col)
    ax2.legend(framealpha=0.9, loc="upper right")
    apply_publication_spines(ax2)

  # 3. QED vs MPO
  if qed_col in hits_df.columns and mpo_col in hits_df.columns:
    ax3.scatter(
        hits_df[qed_col],
        hits_df[mpo_col],
        color=PUBLICATION_PALETTE["purple"],
        s=60,
        edgecolors="#0F172A",
        alpha=0.85,
    )
    ax3.set_title(
        "C. Hit Drug-Likeness & MPO Distribution", fontweight="bold", pad=10
    )
    ax3.set_xlabel(qed_col)
    ax3.set_ylabel(mpo_col)
    apply_publication_spines(ax3)

  # 4. Summary Card
  ax4.axis("off")
  summary_text = (
      "Screening Campaign Summary:\n\n"
      f"• Total Library Deck: {counts[0]:,} molecules\n"
      f"• Hits Passing Filters: {len(hits_df):,} molecules\n"
      f"• Overall Hit Rate: {len(hits_df)/max(1, counts[0])*100:.2f}%\n"
      f"• Destination Output: {os.path.basename(output_png)}"
  )
  ax4.text(
      0.1,
      0.85,
      summary_text,
      fontsize=11,
      verticalalignment="top",
      bbox=dict(
          boxstyle="round,pad=0.6,rounding_size=0.3",
          facecolor=PUBLICATION_PALETTE["card"],
          edgecolor=PUBLICATION_PALETTE["border"],
          alpha=0.95,
      ),
  )
  ax4.set_title("D. Campaign Attrition Statistics", fontweight="bold", pad=10)

  fig.savefig(abs_out, dpi=300)
  plt.close(fig)
  print(f"Saved Archetype 1 plot: {abs_out}")
  return abs_out


def plot_ml_qsar_parity(
    y_exp: np.ndarray,
    y_pred: np.ndarray,
    model_name: str,
    feature_importance_df: Optional[pd.DataFrame],
    output_png: str,
    target_label: str = "pIC50",
    anchor_annotations: Optional[List[Tuple[float, float, str]]] = None,
    split_labels: Optional[Union[np.ndarray, List[Any], pd.Series]] = None,
) -> str:
  """Archetype 3: QSAR ML Parity plot, residuals, and feature attributions with anti-overlap formatting."""
  abs_out = os.path.abspath(output_png)
  os.makedirs(os.path.dirname(abs_out), exist_ok=True)

  fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
      2, 2, figsize=(15, 11), dpi=300, layout="constrained"
  )
  set_safe_plot_style(font_scale=1.05)

  r2 = r2_score(y_exp, y_pred) if len(y_exp) >= 2 else 0.0
  pr, _ = (
      pearsonr(y_exp, y_pred)
      if len(y_exp) >= 3 and np.std(y_exp) > 1e-9 and np.std(y_pred) > 1e-9
      else (0.0, 0.0)
  )
  rmse = np.sqrt(mean_squared_error(y_exp, y_pred)) if len(y_exp) >= 1 else 0.0

  has_split = False
  is_tr = None
  is_te = None
  if split_labels is not None and len(split_labels) == len(y_exp):
    split_str = np.array([str(s).strip() for s in split_labels])
    is_te = np.array([
        bool(re.search(r"test|holdout|val|external", s, re.IGNORECASE))
        for s in split_str
    ])
    is_tr = np.array(
        [bool(re.search(r"train|cv", s, re.IGNORECASE)) for s in split_str]
    )
    if np.any(is_te) and np.any(is_tr):
      has_split = True

  # 1. Parity Plot with clean Corner Stats Box
  if has_split:
    ax1.scatter(
        y_exp[is_tr],
        y_pred[is_tr],
        color=PUBLICATION_PALETTE["primary"],
        s=70,
        edgecolors="#0F172A",
        lw=0.8,
        zorder=3,
        alpha=0.85,
        label=f"Train CV (N={np.sum(is_tr)})",
    )
    ax1.scatter(
        y_exp[is_te],
        y_pred[is_te],
        color=PUBLICATION_PALETTE["danger"],
        marker="D",
        s=85,
        edgecolors="#FB7185",
        lw=0.8,
        zorder=4,
        alpha=1.0,
        label=f"Holdout Test (N={np.sum(is_te)})",
    )

    r2_tr = r2_score(y_exp[is_tr], y_pred[is_tr]) if np.sum(is_tr) >= 2 else 0.0
    pr_tr, _ = (
        pearsonr(y_exp[is_tr], y_pred[is_tr])
        if np.sum(is_tr) >= 3
        and np.std(y_exp[is_tr]) > 1e-9
        and np.std(y_pred[is_tr]) > 1e-9
        else (0.0, 0.0)
    )
    r2_te = r2_score(y_exp[is_te], y_pred[is_te]) if np.sum(is_te) >= 2 else 0.0
    pr_te, _ = (
        pearsonr(y_exp[is_te], y_pred[is_te])
        if np.sum(is_te) >= 3
        and np.std(y_exp[is_te]) > 1e-9
        and np.std(y_pred[is_te]) > 1e-9
        else (0.0, 0.0)
    )
    rmse_te = (
        np.sqrt(mean_squared_error(y_exp[is_te], y_pred[is_te]))
        if np.sum(is_te) >= 1
        else 0.0
    )
    stats_str = (
        f"Train CV $Q^2 = {r2_tr:.2f}$\n"
        f"Train CV $r = {pr_tr:.2f}$\n"
        f"Test $R^2 = {r2_te:.2f}$\n"
        f"Test RMSE = {rmse_te:.2f}\n"
        f"Test $r = {pr_te:.2f}"
    )
  else:
    ax1.scatter(
        y_exp,
        y_pred,
        color=PUBLICATION_PALETTE["primary"],
        s=75,
        edgecolors="#0F172A",
        lw=0.8,
        zorder=3,
        label="Evaluated Congeners",
    )
    stats_str = (
        f"$R^2 = {r2:.2f}$\nPearson $r = {pr:.2f}$\nRMSE = {rmse:.2f}\n$N ="
        f" {len(y_exp)}$"
    )

  min_v = min(min(y_exp), min(y_pred)) - 0.5
  max_v = max(max(y_exp), max(y_pred)) + 0.5
  ax1.plot(
      [min_v, max_v],
      [min_v, max_v],
      color="#64748B",
      linestyle="--",
      lw=1.5,
      alpha=0.8,
      label="1:1 Parity Line",
  )

  ax1.text(
      0.05,
      0.95,
      stats_str,
      transform=ax1.transAxes,
      fontsize=9.5,
      verticalalignment="top",
      bbox=dict(
          boxstyle="round,pad=0.4,rounding_size=0.3",
          fc=PUBLICATION_PALETTE["card"],
          ec=PUBLICATION_PALETTE["border"],
          alpha=0.95,
      ),
  )

  # Selective Non-overlapping Annotations
  if anchor_annotations:
    for x_a, y_a, lbl in anchor_annotations[:3]:
      ax1.annotate(
          short_plot_label(lbl, max_len=18),
          xy=(x_a, y_a),
          xytext=(15, 15),
          textcoords="offset points",
          fontweight="bold",
          fontsize=9,
          color=PUBLICATION_PALETTE["primary_dark"],
          arrowprops=dict(
              arrowstyle="->", color=PUBLICATION_PALETTE["primary_dark"], lw=1.0
          ),
          bbox=dict(
              boxstyle="round,pad=0.2,rounding_size=0.2",
              fc="white",
              ec=PUBLICATION_PALETTE["primary_light"],
              alpha=0.9,
              lw=0.8,
          ),
      )

  ax1.set_title(f"A. {model_name} Parity Plot", fontweight="bold", pad=10)
  ax1.set_xlabel(f"Experimental {target_label}")
  ax1.set_ylabel(f"Predicted {target_label}")
  ax1.set_xlim(min_v, max_v)
  ax1.set_ylim(min_v, max_v)
  ax1.legend(loc="lower right", framealpha=0.9, fontsize=9)
  apply_publication_spines(ax1)

  # 2. Residuals vs Fitted
  residuals = y_exp - y_pred
  if has_split:
    ax2.scatter(
        y_pred[is_tr],
        residuals[is_tr],
        color=PUBLICATION_PALETTE["primary"],
        s=65,
        edgecolors="#0F172A",
        lw=0.6,
        zorder=3,
        alpha=0.8,
        label="Train CV",
    )
    ax2.scatter(
        y_pred[is_te],
        residuals[is_te],
        color=PUBLICATION_PALETTE["danger"],
        marker="D",
        s=75,
        edgecolors="#FB7185",
        lw=0.8,
        zorder=4,
        alpha=1.0,
        label="Holdout Test",
    )
    ax2.legend(loc="lower right", framealpha=0.85, fontsize=8)
  else:
    ax2.scatter(
        y_pred,
        residuals,
        color=PUBLICATION_PALETTE["cyan"],
        s=70,
        edgecolors="#0F172A",
        zorder=3,
    )
  ax2.axhline(0, color="#64748B", linestyle="--", lw=1.2)
  ax2.set_title("B. Residuals vs Predicted Values", fontweight="bold", pad=10)
  ax2.set_xlabel(f"Predicted {target_label}")
  ax2.set_ylabel("Residual (Exp - Pred)")
  apply_publication_spines(ax2)

  # 3. Residual Error Histogram
  if has_split and np.sum(is_te) >= 3:
    sns.histplot(
        residuals[is_tr],
        ax=ax3,
        kde=True,
        color=PUBLICATION_PALETTE["purple"],
        bins=12,
        label="Train CV",
        stat="density",
        alpha=0.6,
    )
    sns.histplot(
        residuals[is_te],
        ax=ax3,
        kde=True,
        color=PUBLICATION_PALETTE["danger"],
        bins=8,
        label="Holdout Test",
        stat="density",
        alpha=0.6,
    )
    ax3.legend(loc="upper right", framealpha=0.85, fontsize=8)
  else:
    sns.histplot(
        residuals,
        ax=ax3,
        kde=True,
        color=PUBLICATION_PALETTE["purple"],
        bins=12,
    )
  ax3.set_title("C. Prediction Error Distribution", fontweight="bold", pad=10)
  ax3.set_xlabel("Residual Error")
  apply_publication_spines(ax3)

  # 4. Standardized Feature Attribution with Auto-Expanded Limits
  if feature_importance_df is not None and not feature_importance_df.empty:
    # Filter out anonymous fingerprint hash bits (e.g. FP_42, Bit_842) - only display interpretable descriptors
    interpretable_df = feature_importance_df[
        ~feature_importance_df["feature"]
        .astype(str)
        .str.contains(r"^(?:FP_\d+|Bit_\d+|fp_\d+|bit_\d+)$", regex=True)
    ]
    if not interpretable_df.empty:
      top = interpretable_df.sort_values(by="importance", ascending=True).tail(
          10
      )
      colors = [
          PUBLICATION_PALETTE["success"]
          if v >= 0
          else PUBLICATION_PALETTE["danger"]
          for v in top["importance"]
      ]
      bars = ax4.barh(
          top["feature"],
          top["importance"],
          color=colors,
          edgecolor="#0F172A",
          linewidth=0.8,
      )
      ax4.axvline(0, color="#64748B", lw=1.0)

      # Auto-expand xlim to prevent text label clipping
      max_w = max(abs(top["importance"])) if len(top) > 0 else 1.0
      ax4.set_xlim(-max_w * 1.35, max_w * 1.35)
      for bar in bars:
        w = bar.get_width()
        offset = max_w * 0.04 if w >= 0 else -max_w * 0.04
        ha = "left" if w >= 0 else "right"
        ax4.text(
            w + offset,
            bar.get_y() + bar.get_height() / 2,
            f"{w:.2f}",
            va="center",
            ha=ha,
            fontsize=8.5,
            fontweight="bold",
        )

      ax4.tick_params(axis="y", pad=8)
      ax4.set_title(
          "D. Standardized Feature Importance (β)", fontweight="bold", pad=10
      )
      ax4.set_xlabel("Standardized Coefficient (β)")
      apply_publication_spines(ax4)
    else:
      ax4.axis("off")
  else:
    ax4.axis("off")

  fig.savefig(abs_out, dpi=300)
  plt.close(fig)
  print(f"Saved Archetype 3 plot: {abs_out}")
  return abs_out


def plot_representation_benchmark_comparison(
    benchmark_df: pd.DataFrame,
    output_png: str,
    feature_col: str = "Feature_Set",
    model_col: str = "Model",
    r2_col: Optional[str] = None,
    r_col: Optional[str] = None,
    target_col: Optional[str] = "Target",
) -> str:
  """Plots multi-representation comparison bar chart (Descriptors vs Fingerprints vs Physics vs Hybrid).

  Supports multi-target grouped comparisons (e.g. WT vs Resistant Mutants vs
  Cellular) using hue faceting.

  The benchmark engine emits two column schemas: `Train_CV_R2`/`Train_CV_r`
  when a held-out test split is supplied, and `Q2_LOO`/`Pearson_r` otherwise.
  When `r2_col`/`r_col` are not given explicitly, whichever schema is actually
  present is resolved automatically, so a held-out-test benchmark no longer
  raises `KeyError: 'Q2_LOO'`.
  """
  abs_out = os.path.abspath(output_png)
  os.makedirs(os.path.dirname(abs_out), exist_ok=True)

  r2_col = _resolve_metric_column(
      benchmark_df, r2_col, ["Q2_LOO", "Train_CV_R2", "CV_R2", "Test_R2"], "R²"
  )
  r_col = _resolve_metric_column(
      benchmark_df,
      r_col,
      ["Pearson_r", "Train_CV_r", "Test_Pearson_r"],
      "Pearson r",
  )
  r2_label = {
      "Q2_LOO": r"$Q^2_{\rm LOO}$",
      "CV_R2": r"$R^2_{\rm CV}$",
      "Train_CV_R2": r"$R^2_{\rm CV}$",
      "Test_R2": r"$R^2_{\rm test}$",
  }[r2_col]

  has_multi_target = (
      target_col is not None
      and target_col in benchmark_df.columns
      and benchmark_df[target_col].nunique() > 1
  )

  fig, (ax1, ax2) = plt.subplots(
      1, 2, figsize=(16, 6.5), dpi=300, layout="constrained"
  )
  set_safe_plot_style(font_scale=1.05)

  if has_multi_target:
    best_df = (
        benchmark_df.sort_values(by=r2_col, ascending=False)
        .groupby([feature_col, target_col])
        .first()
        .reset_index()
    )
    # Grouped bar plots
    sns.barplot(
        data=best_df,
        x=feature_col,
        y=r2_col,
        hue=target_col,
        ax=ax1,
        palette=PUBLICATION_CYCLE,
        edgecolor="#0F172A",
        linewidth=0.8,
    )
    sns.barplot(
        data=best_df,
        x=feature_col,
        y=r_col,
        hue=target_col,
        ax=ax2,
        palette=PUBLICATION_CYCLE,
        edgecolor="#0F172A",
        linewidth=0.8,
    )

    for ax in (ax1, ax2):
      for container in ax.containers:
        ax.bar_label(
            container, fmt="%.2f", fontsize=8, padding=2, fontweight="bold"
        )
      ax.legend(
          title=target_col, loc="upper right", framealpha=0.9, fontsize=8.5
      )

    wrapped_x = [
        textwrap.fill(s, width=13) for s in best_df[feature_col].unique()
    ]
    ax1.set_xticks(range(len(wrapped_x)))
    ax1.set_xticklabels(wrapped_x, rotation=15, ha="right", fontsize=9.5)
    ax2.set_xticks(range(len(wrapped_x)))
    ax2.set_xticklabels(wrapped_x, rotation=15, ha="right", fontsize=9.5)

  else:
    # Single target / grouped by feature set
    best_df = (
        benchmark_df.sort_values(by=r2_col, ascending=False)
        .groupby(feature_col)
        .first()
        .reset_index()
    )
    x = np.arange(len(best_df))
    width = 0.4

    # Panel A: Q2 / R2 Comparison
    bars1 = ax1.bar(
        x,
        best_df[r2_col],
        width=width,
        color=PUBLICATION_PALETTE["primary"],
        edgecolor="#0F172A",
        label=f"{r2_label} Score",
    )
    for b in bars1:
      h = b.get_height()
      ax1.text(
          b.get_x() + b.get_width() / 2,
          h + 0.02,
          f"{h:.2f}",
          ha="center",
          va="bottom",
          fontweight="bold",
          fontsize=9,
      )

    # Panel B: Pearson Correlation Comparison
    bars2 = ax2.bar(
        x,
        best_df[r_col],
        width=width,
        color=PUBLICATION_PALETTE["success"],
        edgecolor="#0F172A",
        label="Pearson $r$",
    )
    for b in bars2:
      h = b.get_height()
      ax2.text(
          b.get_x() + b.get_width() / 2,
          h + 0.02,
          f"{h:.2f}",
          ha="center",
          va="bottom",
          fontweight="bold",
          fontsize=9,
      )

    wrapped_x = [textwrap.fill(s, width=14) for s in best_df[feature_col]]
    ax1.set_xticks(x)
    ax1.set_xticklabels(wrapped_x, rotation=20, ha="right", fontsize=9.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(wrapped_x, rotation=20, ha="right", fontsize=9.5)

  # Common titles and bounds
  min_q2 = min(0.0, float(best_df[r2_col].min()) - 0.1)
  max_q2 = max(1.0, float(best_df[r2_col].max()) + 0.15)
  ax1.set_title(
      f"A. Cross-Validated Predictive Power ({r2_label})",
      fontweight="bold",
      pad=10,
  )
  ax1.set_ylabel(f"{r2_label} Score")
  ax1.set_ylim(min_q2, max_q2)
  apply_publication_spines(ax1)

  min_r = min(0.0, float(best_df[r_col].min()) - 0.1)
  max_r = max(1.0, float(best_df[r_col].max()) + 0.15)
  ax2.set_title(
      "B. Linear Correlation to Experiment (Pearson $r$)",
      fontweight="bold",
      pad=10,
  )
  ax2.set_ylabel("Pearson $r$")
  ax2.set_ylim(min_r, max_r)
  apply_publication_spines(ax2)

  fig.savefig(abs_out, dpi=300)
  plt.close(fig)
  print(f"Saved benchmark comparison plot: {abs_out}")
  return abs_out


def plot_multi_target_parity(
    df: pd.DataFrame,
    target_col: str,
    exp_col: str,
    pred_col: str,
    output_png: str,
    model_name: str = "Hybrid ML",
) -> str:
  """Plots multi-target parity regressions in a clean multi-panel grid (1 panel per target/assay).

  Guarantees zero cherry-picking by displaying all evaluated endpoints
  side-by-side with
  verbatim stats cards.
  """
  abs_out = os.path.abspath(output_png)
  os.makedirs(os.path.dirname(abs_out), exist_ok=True)

  targets = sorted(df[target_col].unique())
  num_targets = len(targets)
  if num_targets == 0:
    raise ValueError(
        "No targets found in DataFrame for multi-target parity plotting."
    )

  ncols = min(3, num_targets)
  nrows = int(np.ceil(num_targets / ncols))

  fig, axes = plt.subplots(
      nrows,
      ncols,
      figsize=(6.0 * ncols, 5.5 * nrows),
      dpi=300,
      layout="constrained",
      squeeze=False,
  )
  set_safe_plot_style(font_scale=1.05)

  palette = PUBLICATION_CYCLE * (1 + num_targets // len(PUBLICATION_CYCLE))

  for i, tgt in enumerate(targets):
    r_idx = i // ncols
    c_idx = i % ncols
    ax = axes[r_idx, c_idx]
    sub_df = df[df[target_col] == tgt].dropna(subset=[exp_col, pred_col])

    if len(sub_df) < 2:
      ax.text(
          0.5,
          0.5,
          f"Insufficient data for {tgt} (N={len(sub_df)})",
          ha="center",
          va="center",
      )
      continue

    model_label = (
        sub_df["Best_Model"].iloc[0]
        if "Best_Model" in sub_df.columns
        else model_name
    )
    fset_label = (
        f" [{sub_df['Best_Feature_Set'].iloc[0]}]"
        if "Best_Feature_Set" in sub_df.columns
        else ""
    )
    ax.set_title(
        f"{tgt}: {model_label}{fset_label}",
        fontweight="bold",
        pad=10,
        fontsize=10,
    )
    ax.set_xlabel(f"Experimental {tgt}")
    ax.set_ylabel(f"Predicted {tgt}")

    split_col = None
    for cand in [
        "Split",
        "split",
        "Set",
        "set",
        "Partition",
        "partition",
        "Fold",
        "fold",
        "subset",
        "Subset",
    ]:
      if cand in sub_df.columns:
        split_col = cand
        break

    has_test_split = False
    sub_tr = None
    sub_te = None
    if split_col is not None:
      sub_tr = sub_df[
          sub_df[split_col]
          .astype(str)
          .str.contains(r"train|cv", case=False, na=False)
      ]
      sub_te = sub_df[
          sub_df[split_col]
          .astype(str)
          .str.contains(r"test|holdout|val|external", case=False, na=False)
      ]
      has_test_split = len(sub_tr) > 0 and len(sub_te) > 0

    color = palette[i % len(palette)]

    if has_test_split:

      y_tr_t = sub_tr[exp_col].values
      y_tr_p = sub_tr[pred_col].values
      y_te_t = sub_te[exp_col].values
      y_te_p = sub_te[pred_col].values

      r2_tr = r2_score(y_tr_t, y_tr_p) if len(y_tr_t) >= 2 else 0.0
      pr_tr, _ = (
          pearsonr(y_tr_t, y_tr_p)
          if len(y_tr_t) >= 3
          and np.std(y_tr_t) > 1e-9
          and np.std(y_tr_p) > 1e-9
          else (0.0, 0.0)
      )

      r2_te = r2_score(y_te_t, y_te_p) if len(y_te_t) >= 2 else 0.0
      rmse_te = (
          np.sqrt(mean_squared_error(y_te_t, y_te_p))
          if len(y_te_t) >= 1
          else 0.0
      )
      pr_te, _ = (
          pearsonr(y_te_t, y_te_p)
          if len(y_te_t) >= 3
          and np.std(y_te_t) > 1e-9
          and np.std(y_te_p) > 1e-9
          else (0.0, 0.0)
      )

      ax.scatter(
          y_tr_t,
          y_tr_p,
          color=color,
          s=65,
          edgecolors="#0F172A",
          lw=0.8,
          alpha=0.75,
          label=f"Train CV (N={len(y_tr_t)})",
          zorder=3,
      )
      ax.scatter(
          y_te_t,
          y_te_p,
          color=PUBLICATION_PALETTE["danger"],
          marker="D",
          s=85,
          edgecolors="#0F172A",
          lw=1.0,
          alpha=1.0,
          label=f"Holdout Test (N={len(y_te_t)})",
          zorder=4,
      )

      min_v = min(np.min(sub_df[exp_col]), np.min(sub_df[pred_col])) - 0.4
      max_v = max(np.max(sub_df[exp_col]), np.max(sub_df[pred_col])) + 0.4
      ax.plot(
          [min_v, max_v],
          [min_v, max_v],
          color="#64748B",
          linestyle="--",
          lw=1.2,
          alpha=0.8,
          label="1:1 Parity",
      )

      stats_str = (
          f"Train CV $Q^2 = {r2_tr:.2f}$\n"
          f"Train CV $r = {pr_tr:.2f}$\n"
          f"Test $R^2 = {r2_te:.2f}$\n"
          f"Test RMSE = {rmse_te:.2f}\n"
          f"Test $r = {pr_te:.2f}$"
      )
      ax.legend(loc="lower right", fontsize=8, framealpha=0.85)
      print(
          f"Parity stats for {tgt} ({model_label}{fset_label}): Train CV"
          f" Q2={r2_tr:.2f}, Train CV r={pr_tr:.2f} | Test R2={r2_te:.2f}, Test"
          f" RMSE={rmse_te:.2f}, Test r={pr_te:.2f}"
      )
    else:
      y_t = sub_df[exp_col].values
      y_p = sub_df[pred_col].values

      r2 = r2_score(y_t, y_p)
      pr, _ = (
          pearsonr(y_t, y_p)
          if len(y_t) >= 3 and np.std(y_t) > 1e-9 and np.std(y_p) > 1e-9
          else (0.0, 0.0)
      )
      sr, _ = (
          spearmanr(y_t, y_p)
          if len(y_t) >= 3 and np.std(y_t) > 1e-9 and np.std(y_p) > 1e-9
          else (0.0, 0.0)
      )
      rmse = np.sqrt(mean_squared_error(y_t, y_p))

      ax.scatter(
          y_t,
          y_p,
          color=color,
          s=70,
          edgecolors="#0F172A",
          lw=0.8,
          zorder=3,
          alpha=0.9,
      )

      min_v = min(np.min(y_t), np.min(y_p)) - 0.4
      max_v = max(np.max(y_t), np.max(y_p)) + 0.4
      ax.plot(
          [min_v, max_v],
          [min_v, max_v],
          color="#64748B",
          linestyle="--",
          lw=1.2,
          alpha=0.8,
          label="1:1 Parity",
      )

      stats_str = (
          f"$R^2 = {r2:.2f}$\n"
          f"Pearson $r = {pr:.2f}$\n"
          f"Spearman $\\rho = {sr:.2f}$\n"
          f"RMSE = {rmse:.2f}\n"
          f"$N = {len(y_t)}$"
      )
      print(
          f"Parity stats for {tgt} ({model_label}{fset_label}): R2={r2:.2f},"
          f" Pearson r={pr:.2f}, Spearman rho={sr:.2f}, RMSE={rmse:.2f},"
          f" N={len(y_t)}"
      )

    ax.text(
        0.05,
        0.95,
        stats_str,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(
            boxstyle="round,pad=0.4,rounding_size=0.3",
            fc=PUBLICATION_PALETTE["card"],
            ec=PUBLICATION_PALETTE["border"],
            alpha=0.95,
        ),
    )

    ax.set_xlim(min_v, max_v)
    ax.set_ylim(min_v, max_v)
    apply_publication_spines(ax)

  # Hide any unused subplots
  for j in range(num_targets, nrows * ncols):
    r_idx = j // ncols
    c_idx = j % ncols
    axes[r_idx, c_idx].axis("off")

  fig.savefig(abs_out, dpi=300)
  plt.close(fig)
  print(f"Saved multi-target parity plot: {abs_out}")
  return abs_out


def _detect_col(
    df: pd.DataFrame,
    patterns: Tuple[str, ...],
    exclude: str = "",
) -> Optional[str]:
  """Resolves an optional DataFrame column name by substring match."""
  for col in df.columns:
    if col == exclude:
      continue
    col_lower = str(col).lower()
    if any(pat in col_lower for pat in patterns):
      return col
  return None


def plot_physics_component_correlations(
    df: pd.DataFrame,
    target_col: str,
    docking_col: Optional[str] = None,
    strain_col: Optional[str] = None,
    rotb_col: Optional[str] = None,
    eff_energy_col: Optional[str] = None,
    output_png: str = "physics_correlations.png",
    target_label: str = "Experimental pKi",
) -> str:
  """Archetype 4: 4-Panel physical component correlation regressions with verified stats cards."""
  abs_out = os.path.abspath(output_png)
  os.makedirs(os.path.dirname(abs_out), exist_ok=True)
  if not docking_col:
    docking_col = _detect_col(
        df, ("dock", "vina", "smina", "affinity"), exclude=target_col
    )
  if not strain_col:
    strain_col = _detect_col(df, ("strain",), exclude=target_col)
  if not rotb_col:
    rotb_col = _detect_col(df, ("rotb", "rotatable"), exclude=target_col)
  if not eff_energy_col:
    eff_energy_col = _detect_col(
        df, ("eff", "free_energy", "combined"), exclude=target_col
    )

  fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
      2, 2, figsize=(14, 11), dpi=300, layout="constrained"
  )
  set_safe_plot_style(font_scale=1.05)

  def _add_stats_card(ax, pr, sr, n):
    card_text = (
        f"Pearson $r = {pr:.2f}$\nSpearman $\\rho = {sr:.2f}$\n$N = {n}$"
    )
    ax.text(
        0.05,
        0.95,
        card_text,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(
            boxstyle="round,pad=0.4,rounding_size=0.3",
            fc=PUBLICATION_PALETTE["card"],
            ec=PUBLICATION_PALETTE["border"],
            alpha=0.95,
        ),
    )
    apply_publication_spines(ax)

  # 1. Docking Affinity
  if docking_col and docking_col in df.columns:
    valid1 = df[[docking_col, target_col]].dropna()
    sns.regplot(
        data=valid1,
        x=docking_col,
        y=target_col,
        ax=ax1,
        color=PUBLICATION_PALETTE["primary"],
        scatter_kws={"s": 65, "edgecolor": "#0F172A"},
    )
    pr1, _ = pearsonr(valid1[docking_col], valid1[target_col])
    sr1, _ = spearmanr(valid1[docking_col], valid1[target_col])
    ax1.set_title("A. Docking Affinity vs Potency", fontweight="bold", pad=10)
    ax1.set_xlabel(f"Docking Energy ({docking_col})")
    ax1.set_ylabel(target_label)
    _add_stats_card(ax1, pr1, sr1, len(valid1))

  # 2. Bound Strain Energy
  if strain_col and strain_col in df.columns:
    valid2 = df[[strain_col, target_col]].dropna()
    sns.regplot(
        data=valid2,
        x=strain_col,
        y=target_col,
        ax=ax2,
        color=PUBLICATION_PALETTE["warning"],
        scatter_kws={"s": 65, "edgecolor": "#0F172A"},
    )
    pr2, _ = pearsonr(valid2[strain_col], valid2[target_col])
    sr2, _ = spearmanr(valid2[strain_col], valid2[target_col])
    ax2.set_title(
        "B. Bound Conformational Strain vs Potency", fontweight="bold", pad=10
    )
    ax2.set_xlabel("Conformational Strain ΔEstrain (kcal/mol)")
    ax2.set_ylabel(target_label)
    _add_stats_card(ax2, pr2, sr2, len(valid2))

  # 3. Rotatable Bonds / Entropy Penalty
  if rotb_col and rotb_col in df.columns:
    valid3 = df[[rotb_col, target_col]].dropna()
    sns.regplot(
        data=valid3,
        x=rotb_col,
        y=target_col,
        ax=ax3,
        color=PUBLICATION_PALETTE["purple"],
        scatter_kws={"s": 65, "edgecolor": "#0F172A"},
    )
    pr3, _ = pearsonr(valid3[rotb_col], valid3[target_col])
    sr3, _ = spearmanr(valid3[rotb_col], valid3[target_col])
    ax3.set_title(
        "C. Rotatable Bond Entropy vs Potency", fontweight="bold", pad=10
    )
    ax3.set_xlabel(f"Rotatable Bond Count ({rotb_col})")
    ax3.set_ylabel(target_label)
    _add_stats_card(ax3, pr3, sr3, len(valid3))

  # 4. Effective Combined Free Energy
  if eff_energy_col and eff_energy_col in df.columns:
    valid4 = df[[eff_energy_col, target_col]].dropna()
    sns.regplot(
        data=valid4,
        x=eff_energy_col,
        y=target_col,
        ax=ax4,
        color=PUBLICATION_PALETTE["success"],
        scatter_kws={"s": 65, "edgecolor": "#0F172A"},
    )
    pr4, _ = pearsonr(valid4[eff_energy_col], valid4[target_col])
    sr4, _ = spearmanr(valid4[eff_energy_col], valid4[target_col])
    ax4.set_title(
        "D. Effective Physical Free Energy vs Potency",
        fontweight="bold",
        pad=10,
    )
    ax4.set_xlabel("Effective Binding Free Energy ΔG_eff")
    ax4.set_ylabel(target_label)
    _add_stats_card(ax4, pr4, sr4, len(valid4))

  print(
      f"\n[Physics Component Correlations against {target_label}"
      f" ({target_col})]"
  )
  if docking_col and docking_col in df.columns:
    print(
        f"  Panel A ({docking_col}): Pearson r = {pr1:.2f}, Spearman rho ="
        f" {sr1:.2f}, N = {len(valid1)}"
    )
  if strain_col and strain_col in df.columns:
    print(
        f"  Panel B ({strain_col}): Pearson r = {pr2:.2f}, Spearman rho ="
        f" {sr2:.2f}, N = {len(valid2)}"
    )
  if rotb_col and rotb_col in df.columns:
    print(
        f"  Panel C ({rotb_col}): Pearson r = {pr3:.2f}, Spearman rho ="
        f" {sr3:.2f}, N = {len(valid3)}"
    )
  if eff_energy_col and eff_energy_col in df.columns:
    print(
        f"  Panel D ({eff_energy_col}): Pearson r = {pr4:.2f}, Spearman rho ="
        f" {sr4:.2f}, N = {len(valid4)}"
    )

  fig.savefig(abs_out, dpi=300)
  plt.close(fig)
  print(f"Saved physics component correlations plot: {abs_out}\n")
  return abs_out


def plot_mmpa_sar_landscape(
    mmpa_df: Optional[pd.DataFrame] = None,
    sali_df: Optional[pd.DataFrame] = None,
    output_png: str = "mmpa_sar_landscape.png",
    target_label: str = "pIC50",
    top_k_transforms: int = 12,
    top_k_cliffs: int = 5,
    paired_shifts_df: Optional[pd.DataFrame] = None,
    paired_shifts_label: str = "Scaffold / Topological Matched Pairs",
    paired_x_col: Optional[str] = None,
    paired_y_col: Optional[str] = None,
) -> str:
  """Archetype 5: Matched Molecular Pair (MMPA) Transformations & SALI Activity Cliffs Landscape.

  Renders a publication-quality multi-panel visualization:
  - Panel A: Top recurring chemical transformations ranked by raw standard score
  (|mean ΔActivity| / SD) with standard deviations and pair counts.
  - Panel B: Structure-Activity Landscape Index (SALI) activity cliffs
  (ΔActivity vs Tanimoto similarity) with annotated cliff leads.
  - Panels C & D (optional): When paired_shifts_df is provided, renders
  matched-pair property shifts (using agent-selected paired_x_col and
  paired_y_col)
  and top pair comparisons.
  """
  abs_out = os.path.abspath(output_png)
  os.makedirs(os.path.dirname(abs_out) or ".", exist_ok=True)

  has_mmpa = mmpa_df is not None and not mmpa_df.empty
  has_sali = sali_df is not None and not sali_df.empty
  has_paired = paired_shifts_df is not None and not paired_shifts_df.empty

  if not has_mmpa and not has_sali and not has_paired:
    fig, ax = plt.subplots(figsize=(8, 4), dpi=300)
    ax.text(
        0.5,
        0.5,
        "No matched pairs or activity cliffs found in dataset.",
        ha="center",
        va="center",
        fontsize=12,
    )
    ax.axis("off")
    fig.savefig(abs_out, dpi=300)
    plt.close(fig)
    print(f"Saved empty MMPA placeholder plot: {abs_out}")
    return abs_out

  set_safe_plot_style(font_scale=1.05)

  ax1 = ax2 = ax3 = ax4 = None
  if has_paired and (has_mmpa or has_sali):
    fig, axes = plt.subplots(
        2, 2, figsize=(16.5, 12.5), dpi=300, layout="constrained"
    )
    ax1, ax2 = axes[0, 0], axes[0, 1]
    ax3, ax4 = axes[1, 0], axes[1, 1]
  elif has_mmpa and has_sali:
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(16, 7), dpi=300, layout="constrained"
    )
  elif has_paired:
    fig, (ax3, ax4) = plt.subplots(
        1, 2, figsize=(16, 7), dpi=300, layout="constrained"
    )
  elif has_mmpa:
    fig, ax1 = plt.subplots(
        1, 1, figsize=(9.5, 7), dpi=300, layout="constrained"
    )
  else:
    fig, ax2 = plt.subplots(
        1, 1, figsize=(9.5, 7), dpi=300, layout="constrained"
    )

  # --- Panel A: Top Matched Molecular Pair Transformations ---
  if ax1 is not None and has_mmpa:
    delta_col = next(
        (
            c
            for c in [
                "Delta_Activity",
                "delta_activity",
                "Delta_pActivity",
                "delta_pactivity",
            ]
            if c in mmpa_df.columns
        ),
        None,
    )
    trans_col = next(
        (
            c
            for c in ["Transformation", "transformation", "Transform"]
            if c in mmpa_df.columns
        ),
        None,
    )

    if delta_col and trans_col:
      valid_mmp = mmpa_df.dropna(subset=[delta_col, trans_col]).copy()
      valid_mmp["Clean_Trans"] = (
          valid_mmp[trans_col]
          .astype(str)
          .str.replace(r"\[\*:1\]", "R-", regex=True)
          .str.replace(r"\[\*:2\]", "R'-", regex=True)
      )

      grp = (
          valid_mmp.groupby("Clean_Trans")[delta_col]
          .agg(
              mean_delta="mean",
              std_delta=lambda s: float(s.std()) if len(s) > 1 else 0.0,
              count="count",
          )
          .reset_index()
      )

      grp["abs_mean"] = grp["mean_delta"].abs()
      rec_grp = grp[(grp["count"] >= 2) & (grp["std_delta"] > 0)].copy()
      if not rec_grp.empty:
        rec_grp["std_score"] = rec_grp["abs_mean"] / rec_grp["std_delta"]
        top_grp = rec_grp.sort_values(
            by=["std_score", "count"], ascending=[False, False]
        ).head(top_k_transforms)
      else:
        top_grp = grp.sort_values(
            by=["count", "abs_mean"], ascending=[False, False]
        ).head(top_k_transforms)
      top_grp = top_grp.sort_values(by="mean_delta", ascending=True)

      y_pos = np.arange(len(top_grp))
      colors = [
          PUBLICATION_PALETTE["success"]
          if v >= 0
          else PUBLICATION_PALETTE["danger"]
          for v in top_grp["mean_delta"]
      ]
      xerr = (
          top_grp["std_delta"].values
          if (top_grp["std_delta"] > 0).any()
          else None
      )

      bars = ax1.barh(
          y_pos,
          top_grp["mean_delta"],
          xerr=xerr,
          color=colors,
          edgecolor="#0F172A",
          linewidth=0.8,
          alpha=0.9,
          capsize=3,
          height=0.6,
      )
      ax1.axvline(0, color="#64748B", linestyle="--", linewidth=1.0)
      ax1.set_yticks(y_pos)
      ax1.set_yticklabels(top_grp["Clean_Trans"], fontsize=9)

      max_abs_x = max(
          0.5,
          float((top_grp["abs_mean"] + top_grp["std_delta"]).max()) * 1.75,
      )
      ax1.set_xlim(-max_abs_x * 0.85, max_abs_x)

      for bar, row in zip(bars, top_grp.itertuples()):
        val = float(row.mean_delta)
        std_v = float(getattr(row, "std_delta", 0.0) or 0.0)
        cnt = row.count
        if val >= 0:
          x_txt = val + std_v + max_abs_x * 0.03
          ha = "left"
        else:
          # Place negative bar labels on the empty positive side of x=0
          x_txt = max_abs_x * 0.03
          ha = "left"
        if std_v > 0:
          z_val = val / std_v
          lbl_txt = f"{val:+.2f} ± {std_v:.2f} (z={z_val:+.2f}, N={cnt})"
        else:
          lbl_txt = f"{val:+.2f} (N={cnt})"
        ax1.text(
            x_txt,
            bar.get_y() + bar.get_height() / 2,
            lbl_txt,
            va="center",
            ha=ha,
            fontsize=8.0,
            fontweight="bold",
        )

      ax1.set_title(
          f"A. Top Recurring MMP Shifts Ranked by |Mean Δ{target_label} / SD|",
          fontweight="bold",
          pad=10,
      )
      ax1.set_xlabel(
          f"Mean Shift in {target_label} (Substituent 2 - Substituent 1)"
      )
      ax1.set_ylabel("Chemical Transformation (R1 -> R2)")
      apply_publication_spines(ax1)
    else:
      ax1.text(
          0.5,
          0.5,
          "MMPA data missing Transformation or Delta column",
          ha="center",
          va="center",
      )
      ax1.axis("off")

  # --- Panel B: Structure-Activity Landscape Index (SALI) Activity Cliffs ---
  if ax2 is not None and has_sali:
    sim_col = next(
        (
            c
            for c in ["Similarity", "similarity", "Tanimoto"]
            if c in sali_df.columns
        ),
        None,
    )
    act_delta_col = next(
        (
            c
            for c in [
                "Delta_Activity",
                "delta_activity",
                "Abs_Delta_Activity",
                "abs_delta_activity",
            ]
            if c in sali_df.columns
        ),
        None,
    )
    sali_score_col = next(
        (c for c in ["SALI", "sali"] if c in sali_df.columns), None
    )

    if sim_col and act_delta_col and sali_score_col:
      valid_sali = sali_df.dropna(
          subset=[sim_col, act_delta_col, sali_score_col]
      ).copy()
      valid_sali["Abs_Delta"] = valid_sali[act_delta_col].abs()

      sc = ax2.scatter(
          valid_sali[sim_col],
          valid_sali["Abs_Delta"],
          c=valid_sali[sali_score_col],
          cmap="plasma",
          s=65,
          edgecolors="#0F172A",
          linewidth=0.7,
          alpha=0.85,
          zorder=3,
      )
      cbar = fig.colorbar(sc, ax=ax2, pad=0.02)
      cbar.set_label("SALI Index", fontweight="bold", fontsize=9.5)

      sali_thresh = max(
          2.0, float(np.percentile(valid_sali[sali_score_col], 85))
      )
      ax2.set_title(
          "B. Structure-Activity Landscape (SALI Activity Cliffs)",
          fontweight="bold",
          pad=10,
      )
      ax2.set_xlabel("Chemical Similarity (Tanimoto ECFP4)")
      ax2.set_ylabel(f"|Δ{target_label}|")

      y_min_b = float(valid_sali["Abs_Delta"].min())
      y_max_b = float(valid_sali["Abs_Delta"].max())
      y_span_b = max(0.5, y_max_b - y_min_b)
      ax2.set_ylim(
          bottom=max(0.0, y_min_b - 0.08 * y_span_b),
          top=y_max_b + 0.25 * y_span_b,
      )

      top_cliffs = valid_sali.sort_values(
          by=sali_score_col, ascending=False
      ).head(top_k_cliffs * 4)
      selected_cliffs = []
      placed_coords = []
      max_b_ann = min(top_k_cliffs, 4)
      for r in top_cliffs.itertuples():
        if len(selected_cliffs) >= max_b_ann:
          break
        sim_val = float(getattr(r, sim_col))
        delta_val = float(getattr(r, "Abs_Delta"))
        if any(
            abs(sim_val - ps) < 0.05 and abs(delta_val - pd_v) < 0.10 * y_span_b
            for ps, pd_v in placed_coords
        ):
          continue
        placed_coords.append((sim_val, delta_val))
        selected_cliffs.append(r)

      # Sort selected cliffs by delta_val descending so vertical callout slots never cross or collide
      selected_cliffs.sort(
          key=lambda row_t: float(getattr(row_t, "Abs_Delta")), reverse=True
      )
      placed_box_ys: list[float] = []
      min_y_sep = 0.14 * y_span_b
      for idx_b, r in enumerate(selected_cliffs):
        id1 = short_plot_label(
            getattr(r, "Compound_1", getattr(r, "idx_i", "Mol1"))
        )
        id2 = short_plot_label(
            getattr(r, "Compound_2", getattr(r, "idx_j", "Mol2"))
        )
        sim_val = float(getattr(r, sim_col))
        delta_val = float(getattr(r, "Abs_Delta"))
        sali_val = float(getattr(r, sali_score_col))

        target_y = delta_val + (
            0.06 * y_span_b if idx_b == 0 else -0.04 * y_span_b
        )
        for prev_y in placed_box_ys:
          if target_y > prev_y - min_y_sep:
            target_y = prev_y - min_y_sep
        target_y = max(y_min_b + 0.06 * y_span_b, target_y)
        placed_box_ys.append(target_y)

        target_x = sim_val - 0.06 if sim_val > 0.65 else sim_val + 0.06
        ha_align = "right" if sim_val > 0.65 else "left"

        ax2.annotate(
            f"{id1} / {id2}\nSALI={sali_val:.1f}, Δ={delta_val:.2f}",
            xy=(sim_val, delta_val),
            xytext=(target_x, target_y),
            textcoords="data",
            ha=ha_align,
            va="center",
            fontsize=8,
            fontweight="bold",
            bbox=dict(
                boxstyle="round,pad=0.25,rounding_size=0.2",
                fc="white",
                ec=PUBLICATION_PALETTE["border"],
                alpha=0.95,
                lw=0.8,
            ),
            arrowprops=dict(arrowstyle="->", color="#334155", lw=0.8),
            zorder=5,
        )

      stats_str = (
          f"Total Evaluated Pairs: {len(valid_sali):,}\nMax SALI:"
          f" {valid_sali[sali_score_col].max():.2f}\nMean Similarity:"
          f" {valid_sali[sim_col].mean():.2f}\nCliffs (SALI >="
          f" {sali_thresh:.1f}):"
          f" {(valid_sali[sali_score_col] >= sali_thresh).sum()}"
      )
      ax2.text(
          0.96,
          0.96,
          stats_str,
          transform=ax2.transAxes,
          fontsize=8.5,
          horizontalalignment="right",
          verticalalignment="top",
          bbox=dict(
              boxstyle="round,pad=0.4,rounding_size=0.3",
              fc=PUBLICATION_PALETTE["card"],
              ec=PUBLICATION_PALETTE["border"],
              alpha=0.95,
          ),
      )
      apply_publication_spines(ax2)
    else:
      ax2.text(
          0.5,
          0.5,
          "SALI data missing Similarity, Delta, or SALI score columns",
          ha="center",
          va="center",
      )
      ax2.axis("off")
  elif ax2 is not None:
    ax2.axis("off")

  # --- Panels C & D: Generic Multi-Parameter Paired Property Shifts ---
  if has_paired and ax3 is not None and ax4 is not None:
    pdf = paired_shifts_df.copy()
    delta_cols = [
        c
        for c in pdf.columns
        if c.startswith("Delta_")
        and pd.api.types.is_numeric_dtype(pdf[c])
        and pdf[c].notna().any()
    ]
    sim_col_p = next(
        (
            c
            for c in ["Topological_Similarity", "Similarity", "Tanimoto"]
            if c in pdf.columns and pd.api.types.is_numeric_dtype(pdf[c])
        ),
        None,
    )
    id_cols_p = [
        c
        for c in pdf.columns
        if not c.startswith("SMILES")
        and not c.startswith("Delta_")
        and c != sim_col_p
        and pdf[c].dtype == object
    ]
    id_a = id_cols_p[0] if len(id_cols_p) >= 1 else pdf.columns[0]
    id_b = (
        id_cols_p[1]
        if len(id_cols_p) >= 2
        else pdf.columns[min(1, len(pdf.columns) - 1)]
    )

    def _resolve_paired_col(req_col: Optional[str]) -> Optional[str]:
      if not req_col:
        return None
      for cand in (req_col, f"Delta_{req_col}"):
        if cand in pdf.columns and pd.api.types.is_numeric_dtype(pdf[cand]):
          return cand
      return None

    req_y = _resolve_paired_col(paired_y_col)
    req_x = _resolve_paired_col(paired_x_col)

    if len(delta_cols) >= 2 or (req_x and req_y):
      # Choose y_col (agent-selected paired_y_col or activity/MPO shift)
      y_col = req_y or next(
          (
              c
              for c in ["Delta_pActivity", "Delta_Activity", "Delta_CNS_MPO"]
              if c in delta_cols and pdf[c].notna().sum() >= 3
          ),
          delta_cols[0],
      )
      # Choose x_col (agent-selected paired_x_col or continuous physicochemical/MPO shift)
      preferred_x = [
          "Delta_CNS_MPO",
          "Delta_cLogP",
          "Delta_TPSA",
          "Delta_MW",
          "Delta_RotB",
          "Delta_HBD",
          "Delta_HBA",
      ] + delta_cols
      x_col = (
          req_x
          if (req_x and req_x != y_col)
          else next(
              (
                  c
                  for c in preferred_x
                  if c in pdf.columns
                  and c != y_col
                  and pdf[c].notna().sum() >= 3
              ),
              next((c for c in delta_cols if c != y_col), delta_cols[-1]),
          )
      )
      sub_p = pdf.dropna(subset=[x_col, y_col]).copy()
      c_vals = (
          sub_p[sim_col_p]
          if sim_col_p and sim_col_p in sub_p.columns
          else np.ones(len(sub_p))
      )
      sc3 = ax3.scatter(
          sub_p[x_col],
          sub_p[y_col],
          c=c_vals,
          cmap="viridis",
          s=75,
          edgecolors="#0F172A",
          linewidth=0.8,
          alpha=0.9,
          zorder=3,
      )
      if sim_col_p and sim_col_p in sub_p.columns:
        cbar3 = fig.colorbar(sc3, ax=ax3, pad=0.02)
        cbar3.set_label(
            sim_col_p.replace("_", " "), fontweight="bold", fontsize=9.5
        )
      ax3.axhline(0, color="#64748B", linestyle="--", linewidth=1.0)
      ax3.axvline(0, color="#64748B", linestyle=":", linewidth=1.0)

      ymin_p, ymax_p = float(sub_p[y_col].min()), float(sub_p[y_col].max())
      yspan_p = max(0.5, ymax_p - ymin_p)
      xmin_p, xmax_p = float(sub_p[x_col].min()), float(sub_p[x_col].max())
      xspan_p = max(1.0, xmax_p - xmin_p)
      ax3.set_ylim(ymin_p - 0.12 * yspan_p, ymax_p + 0.28 * yspan_p)
      ax3.set_xlim(xmin_p - 0.15 * xspan_p, xmax_p + 0.25 * xspan_p)

      y_short = _short_metric_label(y_col)
      x_short = _short_metric_label(x_col)
      if sim_col_p and sim_col_p in sub_p.columns:
        sub_p["_rank_score"] = sub_p[y_col].abs() / (
            1.05 - sub_p[sim_col_p].clip(0.0, 0.99)
        )
      else:
        sub_p["_rank_score"] = sub_p[y_col].abs()
      ranked_p = sub_p.sort_values(by="_rank_score", ascending=False)

      top_p_rows = []
      seen_c_cmpds: set[str] = set()
      for r in ranked_p.itertuples():
        if len(top_p_rows) >= 3:
          break
        ca = short_plot_label(getattr(r, id_a, "PairA"))
        cb = short_plot_label(getattr(r, id_b, "PairB"))
        if ca in seen_c_cmpds or cb in seen_c_cmpds:
          continue
        seen_c_cmpds.update((ca, cb))
        top_p_rows.append(r)

      top_p_rows.sort(
          key=lambda row_t: float(getattr(row_t, y_col, 0.0)), reverse=True
      )
      placed_c_ys: list[float] = []
      min_c_sep = 0.16 * yspan_p
      for idx, r in enumerate(top_p_rows):
        ca = short_plot_label(getattr(r, id_a, "PairA"))
        cb = short_plot_label(getattr(r, id_b, "PairB"))
        xv = float(getattr(r, x_col, 0.0))
        yv = float(getattr(r, y_col, 0.0))
        target_cy = yv + (0.10 * yspan_p if yv >= 0 else -0.10 * yspan_p)
        for prev_cy in placed_c_ys:
          if abs(target_cy - prev_cy) < min_c_sep:
            target_cy = prev_cy - min_c_sep
        placed_c_ys.append(target_cy)
        target_cx = xv + (
            0.10 * xspan_p
            if xv < (xmin_p + 0.50 * xspan_p)
            else -0.10 * xspan_p
        )
        ha_p = "left" if xv < (xmin_p + 0.50 * xspan_p) else "right"
        ax3.annotate(
            f"{ca} -> {cb}\n{y_short}={yv:+.2f}, {x_short}={xv:+.1f}",
            xy=(xv, yv),
            xytext=(target_cx, target_cy),
            textcoords="data",
            ha=ha_p,
            va="center",
            fontsize=8,
            fontweight="bold",
            bbox=dict(
                boxstyle="round,pad=0.25,rounding_size=0.2",
                fc="white",
                ec=PUBLICATION_PALETTE["border"],
                alpha=0.95,
                lw=0.8,
            ),
            arrowprops=dict(arrowstyle="->", color="#334155", lw=0.8),
            zorder=5,
        )

      ax3.set_title(
          f"C. {paired_shifts_label} ({y_col.replace('_', ' ')} vs"
          f" {x_col.replace('_', ' ')})",
          fontweight="bold",
          pad=10,
      )
      ax3.set_xlabel(x_col.replace("_", " "))
      ax3.set_ylabel(y_col.replace("_", " "))
      apply_publication_spines(ax3)

      # Panel D: Top Paired Shifts Bar Chart (deduplicated per compound so 10 distinct pairs are shown)
      dedup_idx = []
      cmpd_counts: dict[str, int] = {}
      for idx_r, r in ranked_p.iterrows():
        if len(dedup_idx) >= 10:
          break
        ca = short_plot_label(r.get(id_a, "A"))
        cb = short_plot_label(r.get(id_b, "B"))
        if cmpd_counts.get(ca, 0) >= 1 or cmpd_counts.get(cb, 0) >= 1:
          continue
        cmpd_counts[ca] = cmpd_counts.get(ca, 0) + 1
        cmpd_counts[cb] = cmpd_counts.get(cb, 0) + 1
        dedup_idx.append(idx_r)
      top_b = sub_p.loc[dedup_idx] if len(dedup_idx) >= 4 else ranked_p.head(10)
      top_b = top_b.sort_values(by=y_col, ascending=True)
      p_labels = [
          f"{short_plot_label(getattr(r, id_a, 'A'))} ->"
          f" {short_plot_label(getattr(r, id_b, 'B'))}"
          for r in top_b.itertuples()
      ]
      y_pos4 = np.arange(len(top_b))
      colors4 = [
          PUBLICATION_PALETTE["primary"]
          if v >= 0
          else PUBLICATION_PALETTE["warning"]
          for v in top_b[y_col]
      ]
      bars4 = ax4.barh(
          y_pos4,
          top_b[y_col],
          color=colors4,
          edgecolor="#0F172A",
          linewidth=0.8,
          alpha=0.9,
          height=0.6,
      )
      ax4.axvline(0, color="#64748B", linestyle="--", linewidth=1.0)
      ax4.set_yticks(y_pos4)
      ax4.set_yticklabels(p_labels, fontsize=8.5)
      max_abs_p = max(0.5, float(top_b[y_col].abs().max()) * 1.65)
      ax4.set_xlim(-max_abs_p * 0.75, max_abs_p)
      for bar, r in zip(bars4, top_b.itertuples()):
        val = float(getattr(r, y_col))
        xval = float(getattr(r, x_col, 0.0))
        if val >= 0:
          x_txt = val + max_abs_p * 0.03
          ha_b = "left"
        else:
          # Place negative bar labels on the empty positive side of x=0
          x_txt = max_abs_p * 0.03
          ha_b = "left"
        ax4.text(
            x_txt,
            bar.get_y() + bar.get_height() / 2,
            f"{val:+.2f} ({x_short}={xval:+.1f})",
            va="center",
            ha=ha_b,
            fontsize=8,
            fontweight="bold",
        )
      ax4.set_title(
          f"D. Top {paired_shifts_label} Ranked by |{y_col.replace('_', ' ')}|",
          fontweight="bold",
          pad=10,
      )
      ax4.set_xlabel(
          f"{y_col.replace('_', ' ')} ({id_b.replace('_', ' ')} -"
          f" {id_a.replace('_', ' ')})"
      )
      ax4.set_ylabel(
          f"Compound Pair ({id_a.replace('_', ' ')} ->"
          f" {id_b.replace('_', ' ')})"
      )
      apply_publication_spines(ax4)
    else:
      ax3.axis("off")
      ax4.axis("off")

  fig.savefig(abs_out, dpi=300)
  plt.close(fig)
  print(f"Saved MMPA and SALI SAR Landscape plot: {abs_out}")
  return abs_out


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Multi-panel publication plot archetypes for QSAR and physics"
          " workflows."
      )
  )
  parser.add_argument(
      "--mode",
      type=str,
      choices=["attrition", "parity", "physics_correlations", "mmpa"],
      help="Plot archetype to render.",
  )
  parser.add_argument("--input_csv", type=str, help="Input CSV dataset path.")
  parser.add_argument(
      "--mmpa_csv",
      type=str,
      default=None,
      help="MMPA pairs CSV path (for --mode mmpa).",
  )
  parser.add_argument(
      "--sali_csv",
      type=str,
      default=None,
      help="SALI cliffs CSV path (for --mode mmpa).",
  )
  parser.add_argument(
      "--paired_csv",
      type=str,
      default=None,
      help="Topological/scaffold matched pairs CSV path (for --mode mmpa).",
  )
  parser.add_argument(
      "--paired_x_col",
      "--x_col",
      dest="paired_x_col",
      type=str,
      default=None,
      help=(
          "Agent-selected X-axis property shift column for Panel C (e.g.,"
          " Delta_CNS_MPO, Delta_cLogP, Delta_TPSA, Delta_MW)."
      ),
  )
  parser.add_argument(
      "--paired_y_col",
      "--y_col",
      dest="paired_y_col",
      type=str,
      default=None,
      help="Optional Y-axis potency/MPO shift column for Panel C and D.",
  )
  parser.add_argument(
      "--output_png",
      "-o",
      type=str,
      default="plot_output.png",
      help="Output PNG path.",
  )
  parser.add_argument(
      "--target_col",
      type=str,
      default="pIC50",
      help="Experimental potency column.",
  )
  parser.add_argument(
      "--pred_col",
      type=str,
      default="pred_pIC50",
      help="Predicted potency column.",
  )
  parser.add_argument(
      "--model_name",
      type=str,
      default="QSAR Model",
      help="Model label used in the parity panel title.",
  )
  parser.add_argument(
      "--split_col",
      type=str,
      default=None,
      help=(
          "Column indicating data split (train/test/CV/holdout). Auto-detected"
          " if omitted."
      ),
  )
  parser.add_argument(
      "--hits_csv",
      type=str,
      default=None,
      help="Optional filtered hits CSV path (for --mode attrition).",
  )
  parser.add_argument(
      "--docking_col",
      type=str,
      default=None,
      help=(
          "Docking score column (auto-detected for --mode"
          " physics_correlations)."
      ),
  )
  parser.add_argument(
      "--strain_col",
      type=str,
      default=None,
      help=(
          "Strain energy column (auto-detected for --mode"
          " physics_correlations)."
      ),
  )
  parser.add_argument(
      "--rotb_col",
      type=str,
      default=None,
      help=(
          "Rotatable bonds column (auto-detected for --mode"
          " physics_correlations)."
      ),
  )
  parser.add_argument(
      "--eff_energy_col",
      type=str,
      default=None,
      help=(
          "Effective free energy column (auto-detected for --mode"
          " physics_correlations)."
      ),
  )
  args = parser.parse_args()

  set_safe_plot_style()

  if args.mode == "mmpa":
    mmpa_df = (
        pd.read_csv(args.mmpa_csv)
        if args.mmpa_csv and os.path.exists(args.mmpa_csv)
        else None
    )
    sali_df = (
        pd.read_csv(args.sali_csv)
        if args.sali_csv and os.path.exists(args.sali_csv)
        else None
    )
    paired_df = (
        pd.read_csv(args.paired_csv)
        if args.paired_csv and os.path.exists(args.paired_csv)
        else None
    )
    plot_mmpa_sar_landscape(
        mmpa_df=mmpa_df,
        sali_df=sali_df,
        output_png=args.output_png,
        target_label=args.target_col,
        paired_shifts_df=paired_df,
        paired_x_col=args.paired_x_col,
        paired_y_col=args.paired_y_col,
    )
  elif args.mode and args.input_csv:
    df = pd.read_csv(args.input_csv)
    if args.mode == "attrition":
      if args.hits_csv and os.path.exists(args.hits_csv):
        hits_df = pd.read_csv(args.hits_csv)
      else:
        hit_col = _detect_col(df, ("hit", "pass"))
        hits_df = df[df[hit_col].astype(bool)].copy() if hit_col else df.copy()
      stages = ["Screening Deck", "Selected Hits"]
      counts = [len(df), len(hits_df)]
      plot_screening_attrition_and_properties(
          stages=stages,
          counts=counts,
          deck_df=df,
          hits_df=hits_df,
          output_png=args.output_png,
      )
    elif args.mode == "parity":
      missing = [
          c for c in (args.target_col, args.pred_col) if c not in df.columns
      ]
      if missing:
        print(
            f"Error: column(s) {missing} not found in {args.input_csv}. "
            f"Available: {list(df.columns)}",
            file=sys.stderr,
        )
        sys.exit(1)
      valid = df.dropna(subset=[args.target_col, args.pred_col])
      split_col = args.split_col
      if split_col is None:
        for cand in [
            "Split",
            "split",
            "Set",
            "set",
            "Partition",
            "partition",
            "Fold",
            "fold",
            "subset",
            "Subset",
        ]:
          if cand in valid.columns:
            split_col = cand
            break
      split_labels = (
          valid[split_col].to_numpy()
          if split_col and split_col in valid.columns
          else None
      )
      plot_ml_qsar_parity(
          y_exp=valid[args.target_col].to_numpy(),
          y_pred=valid[args.pred_col].to_numpy(),
          model_name=args.model_name,
          feature_importance_df=None,
          output_png=args.output_png,
          target_label=args.target_col,
          split_labels=split_labels,
      )
    elif args.mode == "physics_correlations":
      plot_physics_component_correlations(
          df=df,
          target_col=args.target_col,
          docking_col=args.docking_col,
          strain_col=args.strain_col,
          rotb_col=args.rotb_col,
          eff_energy_col=args.eff_energy_col,
          output_png=args.output_png,
          target_label=args.target_col,
      )
  else:
    print(
        "Plot Archetypes anti-overlap and multi-endpoint library loaded"
        " successfully."
    )


if __name__ == "__main__":
  main()
