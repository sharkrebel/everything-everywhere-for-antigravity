# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "pandas",
# ]
# ///

"""Multi-Objective Optimization using Pareto Frontier analysis on chemical datasets.

Identifies non-dominated chemical leads balancing competing criteria
(e.g. maximizing potency, minimizing molecular weight, optimizing cLogP,
minimizing SA_Score).
Computes Pareto ranks (Rank 1 = frontier) and optionally checks k-NN Morgan
applicability domain.
"""

import argparse
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


def calculate_pareto_ranks(
    df: pd.DataFrame, objectives: List[Tuple[str, bool]]
) -> np.ndarray:
  """Calculates Pareto ranks (1 = non-dominated frontier, 2 = second layer, etc.).

  Args:
      df: DataFrame containing the data.
      objectives: List of tuples (column_name, is_maximize).

  Returns:
      1D array of integer ranks for each row.
  """
  n = len(df)
  if n == 0:
    return np.array([], dtype=int)

  # Extract score matrix and invert minimized objectives so all are maximized
  scores = np.zeros((n, len(objectives)), dtype=float)
  for idx, (col, is_max) in enumerate(objectives):
    val = (
        pd.to_numeric(df[col], errors="coerce")
        .fillna(-np.inf if is_max else np.inf)
        .values
    )
    scores[:, idx] = val if is_max else -val

  ranks = np.zeros(n, dtype=int)
  current_rank = 1
  remaining_indices = np.arange(n)

  while len(remaining_indices) > 0:
    sub_scores = scores[remaining_indices]
    num_sub = len(remaining_indices)
    is_efficient = np.ones(num_sub, dtype=bool)

    for i in range(num_sub):
      if not is_efficient[i]:
        continue
      # A point j dominates i if j >= i everywhere and j > i somewhere
      for j in range(num_sub):
        if i != j and is_efficient[j]:
          if np.all(sub_scores[j] >= sub_scores[i]) and np.any(
              sub_scores[j] > sub_scores[i]
          ):
            is_efficient[i] = False
            break

    frontier_indices = remaining_indices[is_efficient]
    ranks[frontier_indices] = current_rank
    remaining_indices = remaining_indices[~is_efficient]
    current_rank += 1

  return ranks


def evaluate_applicability_domain(
    train_smiles: List[str],
    query_smiles: List[str],
    k: int = 3,
    percentile_threshold: float = 95.0,
) -> List[Dict[str, any]]:
  """Evaluates whether query molecules reside within the Morgan k-NN applicability domain."""
  gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
  train_fps = [
      gen.GetFingerprint(Chem.MolFromSmiles(s))
      for s in train_smiles
      if Chem.MolFromSmiles(s)
  ]
  query_fps = [
      gen.GetFingerprint(Chem.MolFromSmiles(s))
      if Chem.MolFromSmiles(s)
      else None
      for s in query_smiles
  ]

  if not train_fps:
    return [{"in_domain": True, "avg_knn_dist": 0.0} for _ in query_smiles]

  # Training internal k-NN distance threshold
  train_knn_dists = []
  for i, fp in enumerate(train_fps):
    sims = sorted(
        DataStructs.BulkTanimotoSimilarity(fp, train_fps), reverse=True
    )
    # Skip self-similarity at index 0
    knn_sims = sims[1 : k + 1] if len(sims) > 1 else [1.0]
    train_knn_dists.append(1.0 - np.mean(knn_sims))
  cutoff = float(np.percentile(train_knn_dists, percentile_threshold))

  results = []
  for q_fp in query_fps:
    if q_fp is None:
      results.append(
          {"in_domain": False, "avg_knn_dist": 1.0, "cutoff": cutoff}
      )
      continue
    sims = sorted(
        DataStructs.BulkTanimotoSimilarity(q_fp, train_fps), reverse=True
    )[:k]
    dist = float(1.0 - np.mean(sims)) if sims else 1.0
    results.append({
        "in_domain": bool(dist <= cutoff),
        "avg_knn_dist": round(dist, 3),
        "cutoff": round(cutoff, 3),
    })
  return results


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Multi-Objective Pareto Frontier Optimization on chemical datasets."
      )
  )
  parser.add_argument(
      "--input-csv", required=True, help="Input CSV containing compound data."
  )
  parser.add_argument(
      "--smiles-col", default="SMILES", help="Name of SMILES column."
  )
  parser.add_argument(
      "--maximize",
      nargs="+",
      default=[],
      help="Columns to maximize (e.g. --maximize pIC50 Solubility).",
  )
  parser.add_argument(
      "--minimize",
      nargs="+",
      default=[],
      help="Columns to minimize (e.g. --minimize MW cLogP SA_Score).",
  )
  parser.add_argument(
      "--training-csv",
      default=None,
      help="Optional training CSV for applicability domain checking.",
  )
  parser.add_argument(
      "--output-csv",
      default="pareto_optimized.csv",
      help="Output CSV with Pareto ranks.",
  )
  args = parser.parse_args()

  df = pd.read_csv(args.input_csv)
  objectives = []
  for col in args.maximize:
    if col not in df.columns:
      raise ValueError(f"Maximize column '{col}' not found in {args.input_csv}")
    objectives.append((col, True))
  for col in args.minimize:
    if col not in df.columns:
      raise ValueError(f"Minimize column '{col}' not found in {args.input_csv}")
    objectives.append((col, False))

  if not objectives:
    raise ValueError(
        "Must specify at least one --maximize or --minimize objective column."
    )

  print(
      f"Calculating Pareto frontiers across {len(df)} compounds for objectives:"
  )
  for col, is_max in objectives:
    print(f"  - {'Maximize' if is_max else 'Minimize'}: {col}")

  ranks = calculate_pareto_ranks(df, objectives)
  df["Pareto_Rank"] = ranks
  df["Is_Pareto_Front"] = ranks == 1

  if args.training_csv and args.smiles_col in df.columns:
    train_df = pd.read_csv(args.training_csv)
    if args.smiles_col in train_df.columns:
      print(
          f"Evaluating Morgan k-NN applicability domain against {len(train_df)}"
          " training compounds..."
      )
      ad_res = evaluate_applicability_domain(
          train_df[args.smiles_col].dropna().tolist(),
          df[args.smiles_col].fillna("").tolist(),
      )
      df["In_Applicability_Domain"] = [r["in_domain"] for r in ad_res]
      df["AD_Avg_kNN_Dist"] = [r["avg_knn_dist"] for r in ad_res]

  df.to_csv(args.output_csv, index=False)
  front_count = int((ranks == 1).sum())
  print(
      f"Pareto analysis complete! Identified {front_count} non-dominated"
      " compounds on Rank 1 frontier."
  )
  print(f"Results written to {args.output_csv}")


if __name__ == "__main__":
  main()
