# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "rdkit>=2024.3.1",
#   "numpy>=1.26.0",
#   "pandas>=2.1.0",
#   "matplotlib>=3.8.0",
# ]
# ///
"""Multi-Parameter Optimization (MPO), LipE/LE & Pareto Radar Plotter.

Computes exact medicinal chemistry efficiency indices (Ligand Efficiency LE,
Lipophilic Efficiency LipE/LLE, LELP, CNS MPO, QED, SA Score), evaluates
Derringer-Suich Target Product Profile (TTP) desirability functions, extracts
exact K-objective non-dominated Pareto frontiers (NSGA-II front ranking), and
renders a 2-panel publication figure:
  - Panel A: pIC50 vs. cLogP optimization landscape with diagonal Iso-LipE
    contour lines (LipE = 3, 4, 5, 6, 7) and highlighted Pareto-optimal front.
  - Panel B: Normalized Multi-Axis Radar / Spider chart comparing the top
    Pareto-optimal leads against the Target Product Profile (TTP) polygon.
"""

import argparse
import json
import math
import os
import re
import sys
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
from rdkit.Chem import QED
from rdkit.Chem import RDConfig
from rdkit.Chem import rdMolDescriptors

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # type: ignore[import-not-found]


def _compute_cns_mpo_approx(
    clogp: float, tpsa: float, mw: float, hbd: int
) -> float:
  """Computes Wager CNS MPO desirability score (0-6 scale, using cLogP for cLogD)."""

  def ramp(
      val: float, low: float, high: float, decreasing: bool = True
  ) -> float:
    if decreasing:
      if val <= low:
        return 1.0
      if val >= high:
        return 0.0
      return (high - val) / (high - low)
    else:
      if val <= low:
        return 0.0
      if val >= high:
        return 1.0
      return (val - low) / (high - low)

  s_clogp = ramp(clogp, 3.0, 5.0, decreasing=True)
  s_clogd = ramp(clogp, 2.0, 4.0, decreasing=True)
  s_mw = ramp(mw, 360.0, 500.0, decreasing=True)
  if 40.0 <= tpsa <= 90.0:
    s_tpsa = 1.0
  elif 20.0 < tpsa < 40.0:
    s_tpsa = (tpsa - 20.0) / 20.0
  elif 90.0 < tpsa < 120.0:
    s_tpsa = (120.0 - tpsa) / 30.0
  else:
    s_tpsa = 0.0
  s_hbd = ramp(float(hbd), 0.5, 3.5, decreasing=True)
  # Neutral basic pKa proxy component (~1.0 when <= 8.0)
  s_pka = 1.0
  return float(s_clogp + s_clogd + s_mw + s_tpsa + s_hbd + s_pka)


def enrich_medchem_metrics(
    df: pd.DataFrame,
    smiles_col: str = "SMILES",
    potency_col: str = "pIC50",
    ic50_nm_col: str | None = None,
) -> pd.DataFrame:
  """Computes RDKit properties and LE, LipE (LLE), and LELP for each compound."""
  df = df.copy()
  if potency_col not in df.columns:
    if ic50_nm_col and ic50_nm_col in df.columns:
      df[potency_col] = 9.0 - np.log10(
          np.clip(pd.to_numeric(df[ic50_nm_col], errors="coerce"), 1e-6, None)
      )
    elif "IC50_nM" in df.columns:
      df[potency_col] = 9.0 - np.log10(
          np.clip(pd.to_numeric(df["IC50_nM"], errors="coerce"), 1e-6, None)
      )
    else:
      raise ValueError(
          f"Potency column '{potency_col}' (or 'IC50_nM') not found in input"
          f" columns: {list(df.columns)}"
      )

  records = []
  for idx, row in df.iterrows():
    smi = str(row.get(smiles_col, ""))
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
      continue
    pic50 = float(row[potency_col])
    n_heavy = max(1, mol.GetNumHeavyAtoms())
    mw = float(Descriptors.MolWt(mol))
    clogp = float(Crippen.MolLogP(mol))
    tpsa = float(rdMolDescriptors.CalcTPSA(mol))
    hbd = int(rdMolDescriptors.CalcNumHBD(mol))
    hba = int(rdMolDescriptors.CalcNumHBA(mol))
    fsp3 = float(rdMolDescriptors.CalcFractionCSP3(mol))
    qed_val = float(QED.qed(mol))
    sa_val = float(sascorer.calculateScore(mol))
    cns_mpo = _compute_cns_mpo_approx(clogp, tpsa, mw, hbd)

    # Thermodynamic Ligand Efficiency: LE = 1.372 * pIC50 / N_heavy (kcal/mol/HA at 298.15 K)
    le = (1.372 * pic50) / n_heavy
    # Lipophilic Efficiency (LipE / LLE) = pIC50 - cLogP
    lipe = pic50 - clogp
    # LELP = cLogP / LE
    lelp = clogp / le if abs(le) > 1e-6 else 0.0

    rec = dict(row)
    rec.update({
        "HeavyAtoms": n_heavy,
        "MW": round(mw, 2),
        "cLogP": round(clogp, 3),
        "TPSA": round(tpsa, 2),
        "HBD": hbd,
        "HBA": hba,
        "Fsp3": round(fsp3, 3),
        "QED": round(qed_val, 3),
        "SA_Score": round(sa_val, 3),
        "CNS_MPO": round(cns_mpo, 3),
        "LE": round(le, 3),
        "LipE": round(lipe, 3),
        "LELP": round(lelp, 3),
    })
    records.append(rec)

  return pd.DataFrame(records)


def score_objective_desirability(
    values: np.ndarray, mode: str, low: float, high: float
) -> np.ndarray:
  """Computes normalized [0.01, 1.0] desirability score for an objective."""
  v = np.asarray(values, dtype=float)
  if mode == "max":
    d = (v - low) / max(1e-9, (high - low))
  elif mode == "min":
    d = (high - v) / max(1e-9, (high - low))
  elif mode == "window":
    mid = 0.5 * (low + high)
    half_w = max(1e-9, 0.5 * (high - low))
    d = 1.0 - np.abs(v - mid) / (1.5 * half_w)
  else:
    raise ValueError(f"Unsupported objective mode: {mode!r}")
  return np.clip(d, 0.01, 1.0)


def compute_pareto_ranks(
    obj_matrix: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
  """Computes exact non-dominated Pareto boolean mask and NSGA-II front ranks.

  Args:
    obj_matrix: (N, K) array where HIGHER values are better across all K
      columns.

  Returns:
    Tuple of (is_pareto_front_1, pareto_rank_1_indexed).
  """
  n = obj_matrix.shape[0]
  ranks = np.zeros(n, dtype=int)
  remaining = set(range(n))
  current_rank = 1

  while remaining:
    rem_list = sorted(remaining)
    front = []
    for i in rem_list:
      dominated = False
      for j in rem_list:
        if i == j:
          continue
        if np.all(obj_matrix[j] >= obj_matrix[i]) and np.any(
            obj_matrix[j] > obj_matrix[i]
        ):
          dominated = True
          break
      if not dominated:
        front.append(i)
    for idx in front:
      ranks[idx] = current_rank
      remaining.remove(idx)
    current_rank += 1

  return (ranks == 1), ranks


def parse_objectives(
    spec_str: str | None,
) -> list[dict[str, Any]]:
  """Parses '--objectives' specification or returns standard MedChem TTP."""
  if not spec_str:
    return [
        {"col": "pIC50", "mode": "max", "low": 5.5, "high": 8.5, "weight": 1.5},
        {"col": "LipE", "mode": "max", "low": 2.5, "high": 6.0, "weight": 1.5},
        {"col": "LE", "mode": "max", "low": 0.22, "high": 0.40, "weight": 1.0},
        {
            "col": "CNS_MPO",
            "mode": "max",
            "low": 3.0,
            "high": 5.5,
            "weight": 1.0,
        },
        {"col": "QED", "mode": "max", "low": 0.35, "high": 0.80, "weight": 1.0},
        {
            "col": "SA_Score",
            "mode": "min",
            "low": 2.0,
            "high": 4.5,
            "weight": 1.0,
        },
    ]
  objs = []
  for item in spec_str.split(","):
    parts = [p.strip() for p in item.split(":")]
    if len(parts) < 4:
      raise ValueError(
          f"Invalid objective spec '{item}'. Expected"
          " 'col:mode:low:high[:weight]'."
      )
    objs.append({
        "col": parts[0],
        "mode": parts[1].lower(),
        "low": float(parts[2]),
        "high": float(parts[3]),
        "weight": float(parts[4]) if len(parts) >= 5 else 1.0,
    })
  return objs


def render_mpo_pareto_radar(
    df: pd.DataFrame,
    objectives: list[dict[str, Any]],
    id_col: str,
    potency_col: str,
    output_plot: str,
    top_n_radar: int = 4,
) -> None:
  """Renders the 2-panel Iso-LipE optimization landscape + Spider/Radar overlay."""
  fig = plt.figure(figsize=(13.5, 6.0), dpi=220)
  ax_scatter = fig.add_subplot(1, 2, 1)
  ax_radar = fig.add_subplot(1, 2, 2, polar=True)

  # --- Panel A: pIC50 vs cLogP with Iso-LipE Contours & Pareto Front ---
  x = df["cLogP"].to_numpy(dtype=float)
  y = df[potency_col].to_numpy(dtype=float)
  mpo_d = df["MPO_Desirability"].to_numpy(dtype=float)
  is_p = df["is_pareto_efficient"].to_numpy(dtype=bool)

  x_min, x_max = min(float(np.min(x)) - 0.6, 0.5), max(
      float(np.max(x)) + 0.6, 5.5
  )
  y_min, y_max = min(float(np.min(y)) - 0.6, 4.5), max(
      float(np.max(y)) + 0.6, 9.5
  )
  x_grid = np.linspace(x_min, x_max, 100)

  for lipe_val in [2, 3, 4, 5, 6, 7]:
    y_iso = x_grid + lipe_val
    ax_scatter.plot(
        x_grid,
        y_iso,
        linestyle="--",
        linewidth=0.9,
        color="#94a3b8" if lipe_val < 5 else "#10b981",
        alpha=0.65,
    )
    lbl_x = x_min + 0.35
    lbl_y = lbl_x + lipe_val
    if y_min + 0.2 < lbl_y < y_max - 0.2:
      ax_scatter.text(
          lbl_x,
          lbl_y + 0.08,
          f"LipE = {lipe_val}",
          fontsize=7.5,
          color="#475569" if lipe_val < 5 else "#047857",
          rotation=32,
      )

  sc = ax_scatter.scatter(
      x[~is_p],
      y[~is_p],
      c=mpo_d[~is_p],
      cmap="viridis",
      vmin=0.0,
      vmax=1.0,
      s=55,
      alpha=0.82,
      edgecolors="#334155",
      linewidths=0.6,
      label="Series Analogues",
  )
  ax_scatter.scatter(
      x[is_p],
      y[is_p],
      c=mpo_d[is_p],
      cmap="viridis",
      vmin=0.0,
      vmax=1.0,
      s=115,
      marker="D",
      edgecolors="#dc2626",
      linewidths=1.8,
      label="Pareto-Optimal Frontier",
      zorder=5,
  )

  # Connect 2D (cLogP min, pIC50 max) frontier points
  p_idx = np.where(is_p)[0]
  if len(p_idx) > 1:
    order = p_idx[np.argsort(x[p_idx])]
    ax_scatter.plot(
        x[order],
        y[order],
        color="#dc2626",
        linestyle="-",
        linewidth=1.4,
        alpha=0.75,
        zorder=4,
    )

  def _short_row_label(row_s: pd.Series, fallback: str) -> str:
    raw_v = (
        str(row_s.get(id_col, "")).strip()
        if id_col in row_s and pd.notna(row_s.get(id_col))
        else ""
    )
    if raw_v.lower() in ("nan", "none", "null"):
      raw_v = ""
    m_triv = re.match(
        r"^([A-Za-z][A-Za-z0-9_.\s-]{1,15})\s*\(([^)]+)\)$", raw_v
    )
    if m_triv and not re.search(
        r"\b(rac|rel|cis|trans)\b", m_triv.group(1), re.I
    ):
      return m_triv.group(1).strip()
    if raw_v and len(raw_v) <= 20 and raw_v.count("(") < 2:
      return raw_v
    for alt_c in (
        "paper_compound_id",
        "patent_example_id",
        "compound_id",
        "chembl_id",
        "molecule_chembl_id",
    ):
      if alt_c in row_s and pd.notna(row_s.get(alt_c)):
        av = str(row_s.get(alt_c, "")).strip()
        if av and av.lower() not in ("nan", "none") and len(av) <= 20:
          is_code = bool(
              re.match(
                  r"^(?:cmpd|compound|example|ex\b|chembl|[a-z]{2,}[-_]?\d)",
                  av,
                  re.IGNORECASE,
              )
          )
          return (
              av if (alt_c != "paper_compound_id" or is_code) else f"Cmpd {av}"
          )
    return (
        (raw_v[:15].rstrip(" -_,([") + "…")
        if len(raw_v) > 16
        else (raw_v or fallback)
    )

  for i in p_idx:
    cid = _short_row_label(df.iloc[i], f"#{i + 1}")
    ax_scatter.annotate(
        f"{cid}\n(LipE={df.iloc[i]['LipE']:.1f})",
        xy=(x[i], y[i]),
        xytext=(6, 6),
        textcoords="offset points",
        fontsize=7.5,
        fontweight="bold",
        color="#991b1b",
    )

  cbar = fig.colorbar(sc, ax=ax_scatter, pad=0.02)
  cbar.set_label("Derringer-Suich MPO Desirability (D)", fontsize=9)

  ax_scatter.set_xlim(x_min, x_max)
  ax_scatter.set_ylim(y_min, y_max)
  ax_scatter.set_xlabel(
      "Lipophilicity (RDKit cLogP)", fontsize=10, fontweight="bold"
  )
  ax_scatter.set_ylabel(
      f"Potency ({potency_col})", fontsize=10, fontweight="bold"
  )
  ax_scatter.set_title(
      "A. Potency vs. Lipophilicity & Iso-LipE Landscape",
      fontsize=11,
      fontweight="bold",
  )
  ax_scatter.legend(loc="lower right", frameon=True, fontsize=8.5)
  ax_scatter.spines["top"].set_visible(False)
  ax_scatter.spines["right"].set_visible(False)

  # --- Panel B: Multi-Axis Normalized Radar / Spider Plot ---
  labels = [o["col"] for o in objectives]
  k = len(labels)
  angles = np.linspace(0, 2 * math.pi, k, endpoint=False).tolist()
  angles += angles[:1]

  # TTP Threshold polygon (d = 0.70 target ring)
  ttp_ring = [0.70] * (k + 1)
  ax_radar.plot(
      angles,
      ttp_ring,
      linestyle="--",
      linewidth=1.5,
      color="#059669",
      label="Target Product Profile (d=0.70)",
  )
  ax_radar.fill(angles, ttp_ring, color="#10b981", alpha=0.08)

  top_df = df.sort_values(
      by=["pareto_rank", "MPO_Desirability"], ascending=[True, False]
  ).head(top_n_radar)
  palette = ["#1d4ed8", "#dc2626", "#d97706", "#7c3aed", "#0d9488"]

  for idx_row, (_, r) in enumerate(top_df.iterrows()):
    vals = [float(r[f"d_{o['col']}"]) for o in objectives]
    vals += vals[:1]
    cid = _short_row_label(r, f"Lead {idx_row + 1}")
    color = palette[idx_row % len(palette)]
    ax_radar.plot(
        angles,
        vals,
        linewidth=2.0,
        color=color,
        marker="o",
        markersize=4.5,
        label=(
            f"{cid} (D={r['MPO_Desirability']:.2f},"
            f" Rank={int(r['pareto_rank'])})"
        ),
    )
    ax_radar.fill(angles, vals, color=color, alpha=0.10)

  ax_radar.set_theta_offset(math.pi / 2)
  ax_radar.set_theta_direction(-1)
  ax_radar.set_xticks(angles[:-1])
  ax_radar.set_xticklabels(labels, fontsize=9, fontweight="bold")
  ax_radar.set_ylim(0.0, 1.05)
  ax_radar.set_yticks([0.25, 0.50, 0.75, 1.00])
  ax_radar.set_yticklabels(
      ["0.25", "0.50", "0.75", "1.00"], fontsize=7.5, color="#64748b"
  )
  ax_radar.set_title(
      "B. Multi-Parameter Optimization (MPO) Spider Overlay",
      fontsize=11,
      fontweight="bold",
      pad=18,
  )
  ax_radar.legend(
      loc="upper right", bbox_to_anchor=(1.28, 1.12), frameon=True, fontsize=8
  )

  fig.tight_layout()
  os.makedirs(os.path.dirname(os.path.abspath(output_plot)), exist_ok=True)
  fig.savefig(output_plot, dpi=220, bbox_inches="tight")
  plt.close(fig)


def run_mpo_analysis(
    input_csv: str,
    output_csv: str,
    output_json: str | None = None,
    output_plot: str | None = None,
    smiles_col: str = "SMILES",
    id_col: str = "Compound_ID",
    potency_col: str = "pIC50",
    ic50_nm_col: str | None = None,
    objectives_spec: str | None = None,
    top_n_radar: int = 4,
) -> dict[str, Any]:
  """Runs full MPO desirability, Pareto frontier ranking, and radar plotting."""
  raw_df = pd.read_csv(input_csv)
  if id_col not in raw_df.columns:
    for alt_id in ["ID", "Name", "compound_id", "mol_id", "ChEMBL_ID"]:
      if alt_id in raw_df.columns:
        id_col = alt_id
        break
    else:
      raw_df[id_col] = [f"CMPD_{i + 1:03d}" for i in range(len(raw_df))]

  df = enrich_medchem_metrics(
      raw_df,
      smiles_col=smiles_col,
      potency_col=potency_col,
      ic50_nm_col=ic50_nm_col,
  )
  objectives = parse_objectives(objectives_spec)

  desirability_cols = []
  raw_obj_cols = []
  weights = []
  for o in objectives:
    col = o["col"]
    if col not in df.columns:
      raise ValueError(
          f"Objective column '{col}' not found in enriched dataframe columns:"
          f" {list(df.columns)}"
      )
    col_vals = df[col].to_numpy(dtype=float)
    d_col = f"d_{col}"
    df[d_col] = np.round(
        score_objective_desirability(col_vals, o["mode"], o["low"], o["high"]),
        4,
    )
    desirability_cols.append(d_col)
    if o["mode"] == "max":
      raw_obj_cols.append(col_vals)
    elif o["mode"] == "min":
      raw_obj_cols.append(-col_vals)
    else:
      mid = 0.5 * (o["low"] + o["high"])
      raw_obj_cols.append(-np.abs(col_vals - mid))
    weights.append(float(o["weight"]))

  d_mat = df[desirability_cols].to_numpy(dtype=float)
  raw_obj_mat = np.column_stack(raw_obj_cols)
  w_arr = np.asarray(weights, dtype=float)
  # Derringer-Suich weighted geometric mean desirability
  df["MPO_Desirability"] = np.round(
      np.exp(np.sum(w_arr * np.log(d_mat), axis=1) / np.sum(w_arr)), 4
  )

  is_pareto, ranks = compute_pareto_ranks(raw_obj_mat)
  df["is_pareto_efficient"] = is_pareto
  df["pareto_rank"] = ranks

  df = df.sort_values(
      by=["pareto_rank", "MPO_Desirability"], ascending=[True, False]
  ).reset_index(drop=True)

  os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
  df.to_csv(output_csv, index=False)

  if output_plot:
    render_mpo_pareto_radar(
        df,
        objectives=objectives,
        id_col=id_col,
        potency_col=potency_col,
        output_plot=output_plot,
        top_n_radar=top_n_radar,
    )

  pareto_df = df[df["is_pareto_efficient"]]
  summary = {
      "total_compounds_evaluated": int(len(df)),
      "pareto_frontier_count": int(len(pareto_df)),
      "pareto_frontier_compounds": (
          pareto_df[[
              id_col,
              smiles_col,
              potency_col,
              "cLogP",
              "LipE",
              "LE",
              "CNS_MPO",
              "MPO_Desirability",
          ]].to_dict(orient="records")
      ),
      "objectives_applied": objectives,
      "output_csv": os.path.abspath(output_csv),
      "output_plot": os.path.abspath(output_plot) if output_plot else None,
  }
  if output_json:
    os.makedirs(os.path.dirname(os.path.abspath(output_json)), exist_ok=True)
    with open(output_json, "w", encoding="utf-8") as f:
      json.dump(summary, f, indent=2)

  return summary


def main() -> None:
  parser = argparse.ArgumentParser(
      description=(
          "Compute Medicinal Chemistry Efficiency Metrics (LE, LipE, LELP),"
          " Derringer-Suich MPO Desirability, and Non-Dominated Pareto"
          " Frontiers with Iso-LipE & Radar visualization."
      )
  )
  parser.add_argument(
      "--input", "-i", required=True, help="Input CSV file path."
  )
  parser.add_argument(
      "--output",
      "-o",
      required=True,
      help="Output enriched & Pareto-ranked CSV file path.",
  )
  parser.add_argument(
      "--output-json",
      default=None,
      help="Optional output JSON summary file path.",
  )
  parser.add_argument(
      "--output-plot",
      default=None,
      help=(
          "Optional path (.png or .svg) for the 2-panel Iso-LipE & Radar chart."
      ),
  )
  parser.add_argument(
      "--smiles-col", default="SMILES", help="Column name containing SMILES."
  )
  parser.add_argument(
      "--id-col", default="Compound_ID", help="Column name for compound IDs."
  )
  parser.add_argument(
      "--potency-col",
      default="pIC50",
      help="Column name for log potency (e.g. pIC50 or pKi).",
  )
  parser.add_argument(
      "--ic50-nm-col",
      default=None,
      help=(
          "Optional column name for nanomolar IC50/Ki if pIC50 is not"
          " precomputed."
      ),
  )
  parser.add_argument(
      "--objectives",
      default=None,
      help=(
          "Comma-separated objective specs 'col:mode:low:high[:weight]' (e.g."
          " 'pIC50:max:5.5:8.5:1.5,LipE:max:2.5:6.0:1.5,SA_Score:min:2.0:4.5:1.0')."
      ),
  )
  parser.add_argument(
      "--top-n-radar",
      type=int,
      default=4,
      help="Number of top Pareto leads to overlay on the Spider/Radar plot.",
  )
  args = parser.parse_args()

  try:
    summary = run_mpo_analysis(
        input_csv=args.input,
        output_csv=args.output,
        output_json=args.output_json,
        output_plot=args.output_plot,
        smiles_col=args.smiles_col,
        id_col=args.id_col,
        potency_col=args.potency_col,
        ic50_nm_col=args.ic50_nm_col,
        objectives_spec=args.objectives,
        top_n_radar=args.top_n_radar,
    )
    print(
        f"[SUCCESS] Evaluated {summary['total_compounds_evaluated']} compounds;"
        f" identified {summary['pareto_frontier_count']} Pareto-optimal lead(s)"
        f" -> {args.output}"
    )
  except Exception as exc:
    print(f"[ERROR] {exc}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
  main()
