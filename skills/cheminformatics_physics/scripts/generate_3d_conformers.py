# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "pandas",
#     "scipy",
#     "scikit-learn",
#     "matplotlib",
#     "seaborn",
# ]
# ///

"""Generates 3D conformers with ETKDGv3, energy minimization, and RMSD pruning."""

import argparse
import os
import sys
from typing import List, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem

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

from compat import get_etkdg_params


def generate_conformers(
    mol: Chem.Mol,
    num_confs: int = 50,
    prune_rms_thresh: float = 0.5,
    energy_window: float = 10.0,
    force_field: str = "MMFF94s",
    random_seed: int = 42,
) -> Tuple[Chem.Mol, List[Tuple[int, float]]]:
  """Generates energy-minimized 3D conformers pruned by RMSD.

  Args:
      mol: Input RDKit Mol object.
      num_confs: Number of initial conformers to embed.
      prune_rms_thresh: RMSD threshold for initial pruning during embedding.
      energy_window: Max energy difference (kcal/mol) from global minimum to
        keep.
      force_field: 'MMFF94s', 'MMFF94', or 'UFF'.
      random_seed: Random seed for reproducibility.

  Returns:
      Tuple of (Mol with conformers, list of (conf_id,
      relative_energy_kcal_mol)).
  """
  mol_3d = Chem.AddHs(mol)

  # 1. Setup ETKDGv3 parameters
  params = get_etkdg_params(
      random_seed=random_seed,
      prune_rms_thresh=prune_rms_thresh,
      use_small_ring_torsions=True,
  )

  # 2. Embed multiple conformers with multi-tier fallback cascade
  cids = list(
      AllChem.EmbedMultipleConfs(mol_3d, numConfs=num_confs, params=params)
  )
  if not cids:
    # Tier 2: ETKDG with random coordinates
    params.useRandomCoords = True
    cids = list(
        AllChem.EmbedMultipleConfs(mol_3d, numConfs=num_confs, params=params)
    )

  if not cids:
    # Tier 3: Standard Distance Geometry without ETKDG torsional bounds
    p_dg = AllChem.EmbedParameters()
    p_dg.randomSeed = random_seed
    p_dg.useRandomCoords = True
    p_dg.ignoreSmoothingFailures = True
    p_dg.enforceChirality = False
    p_dg.maxIterations = 100
    p_dg.pruneRmsThresh = prune_rms_thresh
    cids = list(
        AllChem.EmbedMultipleConfs(mol_3d, numConfs=num_confs, params=p_dg)
    )

  if not cids:
    # Tier 4: Embed heavy atoms first, then add coordinates for hydrogens
    mol_heavy = Chem.Mol(mol)
    p_heavy = AllChem.EmbedParameters()
    p_heavy.randomSeed = random_seed
    p_heavy.useRandomCoords = True
    p_heavy.ignoreSmoothingFailures = True
    p_heavy.enforceChirality = False
    heavy_cids = list(
        AllChem.EmbedMultipleConfs(
            mol_heavy, numConfs=num_confs, params=p_heavy
        )
    )
    if heavy_cids:
      mol_3d = Chem.Mol(mol_heavy)
      mol_3d = Chem.AddHs(mol_3d, addCoords=True)
      cids = [conf.GetId() for conf in mol_3d.GetConformers()]

  if not cids:
    raise RuntimeError("Conformer embedding failed for molecule.")

  # 3. Energy minimization
  energies = []
  if force_field.upper().startswith("MMFF"):
    variant = "MMFF94s" if "s" in force_field.lower() else "MMFF94"
    # Check MMFF parameterization using MMFFHasAllMoleculeParams
    if AllChem.MMFFHasAllMoleculeParams(mol_3d):
      results = AllChem.MMFFOptimizeMoleculeConfs(
          mol_3d, mmffVariant=variant, maxIters=500
      )
      energies = [res[1] for res in results]
    else:
      sys.stderr.write("MMFF parameters missing, falling back to UFF.\n")
      results = AllChem.UFFOptimizeMoleculeConfs(mol_3d, maxIters=500)
      energies = [res[1] for res in results]
  else:
    results = AllChem.UFFOptimizeMoleculeConfs(mol_3d, maxIters=500)
    energies = [res[1] for res in results]

  # 4. Sort and calculate relative energies
  conf_energies = list(zip(cids, energies))
  conf_energies.sort(key=lambda x: x[1])

  min_energy = conf_energies[0][1]
  filtered_confs = []
  for cid, e in conf_energies:
    delta_e = e - min_energy
    if delta_e <= energy_window:
      filtered_confs.append((cid, delta_e))

  return mol_3d, filtered_confs


def main():
  parser = argparse.ArgumentParser(
      description="Generate 3D conformers with ETKDGv3 and MMFF/UFF."
  )
  parser.add_argument(
      "--smiles",
      "--input_smiles",
      dest="smiles",
      type=str,
      required=True,
      help="Input SMILES.",
  )
  parser.add_argument(
      "--num_confs",
      type=int,
      default=20,
      help="Number of conformers to generate.",
  )
  parser.add_argument(
      "--output_sdf",
      "--output",
      dest="output_sdf",
      type=str,
      help="Path to save output SDF.",
  )
  parser.add_argument(
      "--energy_window",
      type=float,
      default=10.0,
      help="Energy window (kcal/mol).",
  )
  parser.add_argument(
      "--force_field",
      "--force-field",
      dest="force_field",
      type=str,
      default="MMFF94s",
      choices=["MMFF94s", "MMFF94", "UFF"],
      help=(
          "Force field for conformer minimization (default: MMFF94s). Falls"
          " back to UFF automatically if MMFF parameters are unavailable for"
          " the molecule."
      ),
  )
  parser.add_argument(
      "--prune_rmsd",
      "--prune-rmsd",
      dest="prune_rmsd",
      type=float,
      default=0.5,
      help=(
          "RMSD threshold (Angstrom) for pruning redundant conformers during"
          " embedding (default: 0.5). Raise it for a coarser, more diverse"
          " ensemble."
      ),
  )
  parser.add_argument(
      "--top_k",
      "--top-k",
      dest="top_k",
      type=int,
      default=10,
      help=(
          "Maximum number of conformers to display in console output"
          " (default: 10; all conformers are saved to --output_sdf)."
      ),
  )
  args = parser.parse_args()

  mol = Chem.MolFromSmiles(args.smiles)
  if mol is None:
    sys.stderr.write("Invalid SMILES\n")
    sys.exit(1)

  mol_3d, confs = generate_conformers(
      mol,
      num_confs=args.num_confs,
      energy_window=args.energy_window,
      force_field=args.force_field,
      prune_rms_thresh=args.prune_rmsd,
  )
  print(
      f"Generated {len(confs)} valid conformers within {args.energy_window}"
      " kcal/mol window:"
  )
  for cid, delta_e in confs[: args.top_k]:
    print(f"  Conformer ID {cid:2d}: delta E = {delta_e:.3f} kcal/mol")
  if len(confs) > args.top_k:
    print(
        f"  ... showing top {args.top_k} of {len(confs)} conformers; pass"
        " --output_sdf to save the full ensemble."
    )

  if args.output_sdf:
    with Chem.SDWriter(args.output_sdf) as writer:
      for cid, delta_e in confs:
        mol_3d.SetProp("RelativeEnergy_kcal_mol", f"{delta_e:.3f}")
        writer.write(mol_3d, confId=cid)
    print(f"Saved conformers to {args.output_sdf}")


if __name__ == "__main__":
  main()
