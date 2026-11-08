# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy",
#     "rdkit",
# ]
# ///

"""Torsion Fingerprint Deviation (TFD) matrix calculation and conformer clustering.

Quantifies conformational diversity across 3D conformer ensembles using internal
dihedral angles (rdkit.Chem.TorsionFingerprints) rather than Cartesian
coordinates. Normalizes distances to [0, 1] and performs Taylor-Butina
clustering to extract structurally non-redundant centroid conformers.
"""

import argparse
import json
import sys
from typing import List, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit.ML.Cluster import Butina


def compute_tfd_matrix_and_clusters(
    mol: Chem.Mol,
    threshold: float = 0.10,
) -> Tuple[List[float], List[Tuple[int, ...]]]:
  """Calculates pairwise 1D condensed TFD distance matrix and Taylor-Butina clusters."""
  num_confs = mol.GetNumConformers()
  if num_confs < 2:
    return [], [(0,)]

  # Compute pairwise lower-triangular TFD distance matrix
  tfd_matrix = TorsionFingerprints.GetTFDMatrix(mol)

  # Butina clustering on condensed distance matrix
  clusters = Butina.ClusterData(
      tfd_matrix, num_confs, distThresh=threshold, isDistData=True
  )
  return list(tfd_matrix), list(clusters)


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Compute Torsion Fingerprint Deviation (TFD) and cluster conformers."
      )
  )
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument("--sdf", help="Input multi-conformer SDF file.")
  group.add_argument(
      "--smiles", help="SMILES string to generate conformers for."
  )

  parser.add_argument(
      "--num-confs",
      "--num_confs",
      dest="num_confs",
      type=int,
      default=25,
      help="Number of conformers to generate if using --smiles (default: 25).",
  )
  parser.add_argument(
      "--tfd-cutoff",
      "--tfd_cutoff",
      "--tfd-thresh",
      "--tfd_thresh",
      dest="tfd_cutoff",
      type=float,
      default=0.10,
      help="TFD distance threshold for Butina clustering (default: 0.10).",
  )
  parser.add_argument(
      "--output-sdf",
      "--output_sdf",
      dest="output_sdf",
      default=None,
      help=(
          "Output SDF containing clustered conformers (or only centroids if"
          " --centroids-only)."
      ),
  )
  parser.add_argument(
      "--centroids-only",
      "--centroids_only",
      dest="centroids_only",
      action="store_true",
      help="Only write cluster centroid conformers to output SDF.",
  )
  parser.add_argument(
      "--output-json",
      "--output_json",
      "--output-csv",
      "--output_csv",
      dest="output_json",
      default=None,
      help="Output JSON/CSV report with TFD statistics.",
  )
  args = parser.parse_args()

  if args.sdf:
    suppl = Chem.SDMolSupplier(args.sdf, removeHs=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
      raise ValueError(f"No valid molecules found in {args.sdf}")
    # Merge individual molecules into one multi-conformer Mol
    base_mol = Chem.Mol(mols[0])
    base_mol.RemoveAllConformers()
    for m in mols:
      for conf in m.GetConformers():
        base_mol.AddConformer(Chem.Conformer(conf), assignId=True)
    mol = base_mol
    print(f"Loaded {mol.GetNumConformers()} conformers from {args.sdf}.")
  else:
    print(
        f"Generating {args.num_confs} conformers for {args.smiles} using"
        " ETKDGv3..."
    )
    parsed = Chem.MolFromSmiles(args.smiles)
    if parsed is None:
      # Chem.AddHs(None) raises an opaque Boost.Python.ArgumentError that
      # tells the user nothing about which input was bad.
      print(f"Error: could not parse SMILES: {args.smiles}", file=sys.stderr)
      sys.exit(1)
    mol = Chem.AddHs(parsed)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMultipleConfs(mol, numConfs=args.num_confs, params=params)
    AllChem.MMFFOptimizeMoleculeConfs(mol)
    print(f"Generated and minimized {mol.GetNumConformers()} conformers.")

  num_confs = mol.GetNumConformers()
  tfd_mat, clusters = compute_tfd_matrix_and_clusters(
      mol, threshold=args.tfd_cutoff
  )

  # Compute summary metrics
  mean_tfd = float(np.mean(tfd_mat)) if tfd_mat else 0.0
  max_tfd = float(np.max(tfd_mat)) if tfd_mat else 0.0
  min_tfd = float(np.min(tfd_mat)) if tfd_mat else 0.0

  cluster_info = []
  centroid_indices = []
  for c_idx, members in enumerate(clusters):
    centroid_id = members[0]
    centroid_indices.append(centroid_id)
    cluster_info.append({
        "cluster_id": c_idx + 1,
        "centroid_conformer_id": int(centroid_id),
        "size": len(members),
        "member_conformer_ids": [int(m) for m in members],
    })

  results = {
      "conformer_count": num_confs,
      "tfd_cutoff": args.tfd_cutoff,
      "cluster_count": len(clusters),
      "mean_tfd": round(mean_tfd, 4),
      "max_tfd": round(max_tfd, 4),
      "min_tfd": round(min_tfd, 4),
      "clusters": cluster_info,
      "centroid_conformer_ids": [int(c) for c in centroid_indices],
  }

  print("\nTFD Conformer Clustering Summary:")
  print(f"  Total conformers: {num_confs}")
  print(f"  TFD Cutoff: {args.tfd_cutoff}")
  print(f"  Unique clusters: {len(clusters)}")
  print(
      f"  Pairwise TFD Mean: {mean_tfd:.3f} | Max: {max_tfd:.3f} | Min:"
      f" {min_tfd:.3f}"
  )
  for c in cluster_info[:5]:
    print(
        f"    Cluster {c['cluster_id']}: Centroid Conf"
        f" {c['centroid_conformer_id']} ({c['size']} conformers)"
    )
  if len(cluster_info) > 5:
    print(
        f"    (Showing top 5 of {len(cluster_info)} clusters; pass"
        " --output-json to save all clusters)"
    )

  if args.output_sdf:
    with Chem.SDWriter(args.output_sdf) as writer:
      # Map conformer ID to cluster info
      conf_to_cluster = {}
      for c in cluster_info:
        for m in c["member_conformer_ids"]:
          conf_to_cluster[m] = (
              c["cluster_id"],
              m == c["centroid_conformer_id"],
          )

      for cid in range(num_confs):
        is_centroid = conf_to_cluster.get(cid, (1, True))[1]
        if args.centroids_only and not is_centroid:
          continue
        conf_mol = Chem.Mol(mol)
        # Retain only current conformer
        for other_cid in reversed(range(num_confs)):
          if other_cid != cid:
            conf_mol.RemoveConformer(other_cid)
        cl_id, is_cent = conf_to_cluster.get(cid, (1, True))
        conf_mol.SetIntProp("TFD_Cluster_ID", cl_id)
        conf_mol.SetBoolProp("Is_Conformer_Centroid", is_cent)
        writer.write(conf_mol)
    print(f"Saved conformers to {args.output_sdf}")

  if args.output_json:
    with open(args.output_json, "w") as f:
      json.dump(results, f, indent=2)
    print(f"Saved JSON report to {args.output_json}")


if __name__ == "__main__":
  main()
