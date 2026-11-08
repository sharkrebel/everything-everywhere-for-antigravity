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

"""Performs Butina clustering of compounds based on Morgan fingerprint Tanimoto similarity."""

import argparse
from typing import List, Tuple
import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina


def cluster_smiles(
    smiles_list: List[str],
    cutoff_distance: float = 0.4,
    fp_radius: int = 2,
    fp_size: int = 2048,
) -> Tuple[List[Tuple[int, ...]], List[int]]:
  """Clusters SMILES strings using Butina algorithm on Morgan fingerprints.

  Args:
      smiles_list: List of SMILES strings.
      cutoff_distance: Distance cutoff (1 - TanimotoSimilarity). E.g., 0.4 means
        >= 0.6 similarity.
      fp_radius: Morgan radius (default 2 for ECFP4).
      fp_size: Bit vector length.

  Returns:
      Tuple of (list of clusters with member indices, list of centroid indices).
  """
  mols = []
  valid_map = {}
  for idx, s in enumerate(smiles_list):
    m = Chem.MolFromSmiles(str(s)) if pd.notna(s) else None
    if m is not None:
      valid_map[len(mols)] = idx
      mols.append(m)

  if not mols:
    return [], []

  # Generate Morgan fingerprints
  mfpgen = rdFingerprintGenerator.GetMorganGenerator(
      radius=fp_radius, fpSize=fp_size
  )
  fps = [mfpgen.GetFingerprint(m) for m in mols]

  # Compute 1D lower-triangular distance list (1 - Tanimoto)
  n = len(fps)
  dists = []
  for i in range(1, n):
    sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
    dists.extend([1.0 - x for x in sims])

  # Butina clustering
  raw_clusters = Butina.ClusterData(
      dists, n, distThresh=cutoff_distance, isDistData=True
  )

  # Map back to original indices
  clusters = [
      tuple(valid_map[local_idx] for local_idx in c) for c in raw_clusters
  ]
  centroids = [valid_map[c[0]] for c in raw_clusters]
  return clusters, centroids


def cluster_dataframe(
    df: pd.DataFrame,
    smiles_col: str = "SMILES",
    cutoff_distance: float = 0.4,
    fp_radius: int = 2,
) -> pd.DataFrame:
  """In-Memory DataFrame Clustering: Adds 'Cluster_ID', 'Is_Centroid', and 'Cluster_Size' columns."""
  df = df.copy()
  clusters, centroids = cluster_smiles(
      df[smiles_col].tolist(),
      cutoff_distance=cutoff_distance,
      fp_radius=fp_radius,
  )

  df["Cluster_ID"] = -1
  df["Is_Centroid"] = False
  df["Cluster_Size"] = 0

  centroid_set = set(centroids)
  for c_id, members in enumerate(clusters):
    size = len(members)
    for m_idx in members:
      df.at[m_idx, "Cluster_ID"] = c_id + 1
      df.at[m_idx, "Cluster_Size"] = size
      if m_idx in centroid_set:
        df.at[m_idx, "Is_Centroid"] = True

  return df


def main():
  parser = argparse.ArgumentParser(
      description="Cluster molecules by Morgan fingerprint similarity."
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
      help="SMILES column.",
  )
  parser.add_argument(
      "--cutoff",
      "--cutoff_distance",
      "--cutoff-distance",
      dest="cutoff",
      type=float,
      default=0.4,
      help="Distance cutoff (default: 0.4).",
  )
  parser.add_argument(
      "--output_csv",
      "--output-csv",
      dest="output_csv",
      type=str,
      default=None,
      help="Output CSV path to save clustered compounds.",
  )
  args = parser.parse_args()

  df = pd.read_csv(args.input_csv)
  clustered_df = cluster_dataframe(
      df, smiles_col=args.smiles_col, cutoff_distance=args.cutoff
  )
  num_clusters = (
      int(clustered_df["Cluster_ID"].max()) if len(clustered_df) > 0 else 0
  )
  print(
      f"Clustered {len(df)} compounds into {num_clusters} clusters"
      f" (cutoff={args.cutoff})."
  )
  if args.output_csv:
    clustered_df.to_csv(args.output_csv, index=False)
    print(f"Saved clustered results to {args.output_csv}")


if __name__ == "__main__":
  main()
