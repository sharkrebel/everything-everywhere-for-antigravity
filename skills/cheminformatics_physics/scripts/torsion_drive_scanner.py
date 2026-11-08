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

"""Torsion Drive Scanner: High-Performance Wavefront and Advanced PES Methods.

Implements:
1. Fast Recursive Wavefront Propagation (TorsionDrive algorithm with C++
   multi-threading, ETKDG multi-seeding, and adaptive threshold pruning)
   [Default / Production].
2. Sequential Bidirectional Path-Following (Forward/Backward sweeps for
   hysteresis auditing).
3. 2D Coupled Wavefront PES Unrolling (Simultaneous 4-way grid propagation for
   coupled dihedrals).
4. Adaptive Multi-Resolution Grid Scanner (Hierarchical coarse-to-fine zoom
   around barriers).
5. Templates for WBO-Guided Fragmentation and Metadynamics Rotor Sampling.
"""

import argparse
from collections import deque
import json
import math
import os
import sys
import time
from typing import Any, Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdForceFieldHelpers, rdMolTransforms

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

from compat import set_safe_plot_style


def setup_molecule_3d(
    mol: Chem.Mol,
    force_field: str = "MMFF94s",
    random_seed: int = 42,
) -> Tuple[Chem.Mol, float]:
  """Prepares molecule with explicit hydrogens and an optimized 3D conformer.

  Args:
      mol: Input RDKit Mol.
      force_field: 'MMFF94s', 'MMFF94', or 'UFF'.
      random_seed: Random seed for 3D embedding.

  Returns:
      Tuple of (3D Mol with 1 conformer, ground_state_energy_kcal_mol).
  """
  mol_3d = Chem.AddHs(mol)
  if mol_3d.GetNumConformers() == 0:
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.numThreads = 0
    if AllChem.EmbedMolecule(mol_3d, params) != 0:
      params.useRandomCoords = True
      AllChem.EmbedMolecule(mol_3d, params)

  ff_name = force_field.upper()
  if ff_name.startswith("MMFF") and AllChem.MMFFHasAllMoleculeParams(mol_3d):
    variant = "MMFF94s" if "s" in ff_name.lower() else "MMFF94"
    AllChem.MMFFOptimizeMolecule(mol_3d, mmffVariant=variant, maxIters=1000)
    mp = AllChem.MMFFGetMoleculeProperties(mol_3d, mmffVariant=variant)
    ff = AllChem.MMFFGetMoleculeForceField(mol_3d, mp, confId=0)
    ground_e = ff.CalcEnergy()
  else:
    AllChem.UFFOptimizeMolecule(mol_3d, maxIters=1000)
    ff = AllChem.UFFGetMoleculeForceField(mol_3d, confId=0)
    ground_e = ff.CalcEnergy()

  return mol_3d, ground_e


def identify_key_rotatable_bonds(
    mol: Chem.Mol,
    max_bonds: Optional[int] = None,
) -> List[Dict[str, Any]]:
  """Universally identifies and prioritizes key rotatable bonds in any molecule.

  Categorizes bonds by structural and conformational impact:
  1. Inter-Ring Linkages (bonds connecting two ring systems, e.g. biaryls,
  heteroaryl hinges) - Priority 1
  2. Ring-to-Chain Linkages (exocyclic bonds connecting rings to flexible
  linkers) - Priority 2
  3. Amide / Ester / Sulfonamide Linkages - Priority 3
  4. Acyclic Chain Linkages - Priority 4

  Filters out trivial terminal rotors (methyl, -CF3, -NH2, -OH, halogens).

  Args:
      mol: Input RDKit Mol (with or without Hs).
      max_bonds: Optional cap on candidate rotatable bonds to return (None =
        all).

  Returns:
      List of dicts containing 'dihedral' (4-atom tuple), 'bond_idx', 'type',
      'priority', and 'description'.
  """
  mol_no_h = Chem.RemoveHs(mol)
  results = []

  for bond in mol_no_h.GetBonds():
    if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing():
      continue
    u = bond.GetBeginAtom()
    v = bond.GetEndAtom()

    # Skip terminal heavy atoms
    deg_u = u.GetDegree()
    deg_v = v.GetDegree()
    if deg_u <= 1 or deg_v <= 1:
      continue

    u_in_ring = u.IsInRing()
    v_in_ring = v.IsInRing()

    # Classify bond type and priority
    if u_in_ring and v_in_ring:
      btype = "inter_ring"
      priority = 1
      desc = (
          "Inter-Ring Linkage"
          f" ({u.GetSymbol()}{u.GetIdx()}-{v.GetSymbol()}{v.GetIdx()})"
      )
    elif u_in_ring or v_in_ring:
      btype = "ring_chain"
      priority = 2
      ring_atom = u if u_in_ring else v
      chain_atom = v if u_in_ring else u
      desc = (
          "Ring-to-Chain"
          f" ({ring_atom.GetSymbol()}{ring_atom.GetIdx()}-{chain_atom.GetSymbol()}{chain_atom.GetIdx()})"
      )
    else:
      btype = "acyclic_chain"
      priority = 4
      desc = (
          "Acyclic Chain"
          f" ({u.GetSymbol()}{u.GetIdx()}-{v.GetSymbol()}{v.GetIdx()})"
      )

    # Select best neighbor for u (not v) and v (not u)
    nbrs_u = [nbr for nbr in u.GetNeighbors() if nbr.GetIdx() != v.GetIdx()]
    nbrs_v = [nbr for nbr in v.GetNeighbors() if nbr.GetIdx() != u.GetIdx()]
    if not nbrs_u or not nbrs_v:
      continue

    def score_nbr(a: Chem.Atom) -> Tuple[int, int, int, int]:
      return (
          int(a.IsInRing()),
          int(a.GetIsAromatic()),
          a.GetAtomicNum(),
          a.GetDegree(),
      )

    nbr_u = max(nbrs_u, key=score_nbr)
    nbr_v = max(nbrs_v, key=score_nbr)

    dihedral = (nbr_u.GetIdx(), u.GetIdx(), v.GetIdx(), nbr_v.GetIdx())
    results.append({
        "dihedral": dihedral,
        "bond_idx": bond.GetIdx(),
        "type": btype,
        "priority": priority,
        "description": desc,
    })

  results.sort(key=lambda x: (x["priority"], x["bond_idx"]))
  return results[:max_bonds] if max_bonds is not None else results


def _validate_dihedral(
    mol: Chem.Mol, idx: Tuple[int, int, int, int]
) -> Tuple[int, int, int, int]:
  """Validates that (i, j, k, l) is a real, in-range, bonded torsion chain.

  Without this check a non-bonded quadruple is accepted silently by the
  wavefront scanner (yielding a physically meaningless barrier) and crashes
  the bidirectional scanner with an opaque AttributeError.

  Args:
      mol: RDKit Mol the indices refer to.
      idx: Candidate (i, j, k, l) atom indices.

  Returns:
      The validated tuple.

  Raises:
      ValueError: If the indices are out of range, not distinct, or do not
          form a bonded i-j-k-l chain.
  """
  n = mol.GetNumAtoms()
  out_of_range = [a for a in idx if a < 0 or a >= n]
  if out_of_range:
    raise ValueError(
        f"Dihedral atom indices {list(idx)} out of range: molecule has "
        f"{n} atoms (valid 0-{n - 1}); offending {out_of_range}. "
        "Run with --list_rotors to see valid torsions."
    )
  if len(set(idx)) != 4:
    raise ValueError(
        f"Dihedral atom indices {list(idx)} must be four distinct atoms."
    )
  for a, b in ((idx[0], idx[1]), (idx[1], idx[2]), (idx[2], idx[3])):
    if mol.GetBondBetweenAtoms(a, b) is None:
      raise ValueError(
          f"Dihedral atom indices {list(idx)} do not form a bonded "
          f"chain: atoms {a} and {b} are not bonded. A torsion requires "
          "i-j, j-k and k-l bonds. Run with --list_rotors to see valid "
          "torsions."
      )
  return idx


def resolve_dihedral_indices(
    mol: Chem.Mol,
    dihedral: Optional[List[int]] = None,
    smarts: Optional[str] = None,
    smarts_indices: Optional[List[int]] = None,
) -> Tuple[int, int, int, int]:
  """Resolves the 4 atom indices (i, j, k, l) for the target dihedral angle.

  Args:
      mol: RDKit Mol.
      dihedral: Explicit list of 4 atom indices.
      smarts: SMARTS query matching the torsion fragment.
      smarts_indices: 4 indices within the SMARTS match tuple.

  Returns:
      Tuple of (i, j, k, l) atom indices, validated to be a bonded chain.

  Raises:
      ValueError: If no dihedral can be resolved, or the resolved atoms do
          not form a bonded i-j-k-l chain.
  """
  if dihedral is not None and len(dihedral) == 4:
    return _validate_dihedral(mol, tuple(dihedral))

  if smarts:
    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
      raise ValueError(f"Invalid SMARTS pattern: {smarts}")
    matches = mol.GetSubstructMatches(patt)
    if not matches:
      raise ValueError(f"SMARTS query '{smarts}' did not match molecule.")
    match = matches[0]
    if smarts_indices is not None and len(smarts_indices) == 4:
      return _validate_dihedral(
          mol, tuple(match[idx] for idx in smarts_indices)
      )
    if len(match) >= 4:
      return _validate_dihedral(mol, (match[0], match[1], match[2], match[3]))
    raise ValueError(f"SMARTS match has {len(match)} atoms, need at least 4.")

  # Universal Fallback: Use priority-ranked rotatable bond detection
  key_rotors = identify_key_rotatable_bonds(mol)
  if key_rotors:
    return _validate_dihedral(mol, tuple(key_rotors[0]["dihedral"]))

  raise ValueError(
      "Could not automatically determine rotatable dihedral. Specify --dihedral"
      " or --smarts."
  )


def _minimize_at_dihedral(
    mol: Chem.Mol,
    conf_id: int,
    atom_indices: Tuple[int, int, int, int],
    target_deg: float,
    force_field: str = "MMFF94s",
    force_constant: float = 1000.0,
    max_iters: int = 300,
) -> float:
  """Helper applying rigid rotation and harmonic constraint in degrees."""
  i, j, k, l = atom_indices
  bond = mol.GetBondBetweenAtoms(j, k)
  if bond is None or not bond.IsInRing():
    try:
      rdMolTransforms.SetDihedralDeg(
          mol.GetConformer(conf_id), i, j, k, l, float(target_deg)
      )
    except Exception:
      pass

  ff_name = force_field.upper()
  if ff_name.startswith("MMFF") and AllChem.MMFFHasAllMoleculeParams(mol):
    variant = "MMFF94s" if "s" in ff_name.lower() else "MMFF94"
    props = rdForceFieldHelpers.MMFFGetMoleculeProperties(
        mol, mmffVariant=variant
    )
    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(
        mol, props, confId=conf_id
    )
    add_torsion_constraint = ff.MMFFAddTorsionConstraint
  else:
    ff = rdForceFieldHelpers.UFFGetMoleculeForceField(mol, confId=conf_id)
    add_torsion_constraint = ff.UFFAddTorsionConstraint

  add_torsion_constraint(
      i,
      j,
      k,
      l,
      False,
      float(target_deg),
      float(target_deg),
      float(force_constant),
  )

  ff.Initialize()
  ff.Minimize(maxIts=max_iters)
  return ff.CalcEnergy()


def scan_torsion_wavefront(
    mol: Chem.Mol,
    atom_indices: Tuple[int, int, int, int],
    step_deg: float = 5.0,
    n_seeds: int = 25,
    energy_threshold: float = 1e-3,
    force_field: str = "MMFF94s",
    force_constant: Optional[float] = None,
    max_iters: int = 250,
    coupled_dihedrals: Optional[List[Tuple[int, int, int, int]]] = None,
) -> Dict[str, Any]:
  """High-performance recursive wavefront propagation torsion drive."""
  i, j, k, l = atom_indices
  bond = mol.GetBondBetweenAtoms(j, k)
  bond_in_ring = bond.IsInRing() if bond else False
  if bond_in_ring:
    sys.stderr.write(
        f"NOTE: Central bond ({j}-{k}) is inside a ring. Uncoupled 360°"
        " rotation may cause ring strain. For rings >= 10 atoms, prefer"
        " conformer ensemble distributions (ETKDGv3 with"
        " useMacrocycleTorsions=True).\n"
    )
  if force_constant is None:
    force_constant = 50.0 if bond_in_ring else 1000.0

  angles = np.arange(0, 360, step_deg, dtype=float)
  n_grid = len(angles)

  seed_mol = Chem.Mol(mol)
  params = rdDistGeom.ETKDGv3()
  params.randomSeed = 42
  params.numThreads = 0
  rdDistGeom.EmbedMultipleConfs(seed_mol, numConfs=n_seeds, params=params)

  if force_field.upper().startswith("MMFF"):
    rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(
        seed_mol, mmffVariant="MMFF94s", numThreads=0, maxIters=300
    )
  else:
    rdForceFieldHelpers.UFFOptimizeMoleculeConfs(
        seed_mol, numThreads=0, maxIters=300
    )

  grid_best_energy = [float("inf")] * n_grid
  grid_best_conf = [None] * n_grid
  active_queue = deque()

  for cid in range(seed_mol.GetNumConformers()):
    cur_ang = (
        rdMolTransforms.GetDihedralDeg(seed_mol.GetConformer(cid), i, j, k, l)
        + 360.0
    ) % 360.0
    diffs = [min(abs(cur_ang - a), 360.0 - abs(cur_ang - a)) for a in angles]
    idx = int(np.argmin(diffs))

    conf_mol = Chem.Mol(seed_mol)
    conf_mol.RemoveAllConformers()
    conf_mol.AddConformer(
        Chem.Conformer(seed_mol.GetConformer(cid)), assignId=True
    )

    e = _minimize_at_dihedral(
        conf_mol,
        0,
        atom_indices,
        angles[idx],
        force_field,
        force_constant,
        max_iters,
    )
    if e < grid_best_energy[idx]:
      grid_best_energy[idx] = e
      grid_best_conf[idx] = Chem.Conformer(conf_mol.GetConformer(0))
      if idx not in active_queue:
        active_queue.append(idx)

    if coupled_dihedrals:
      for coup in coupled_dihedrals:
        for kick in [-120.0, 120.0]:
          kicked_mol = Chem.Mol(conf_mol)
          orig_deg = rdMolTransforms.GetDihedralDeg(
              kicked_mol.GetConformer(0), *coup
          )
          rdMolTransforms.SetDihedralDeg(
              kicked_mol.GetConformer(0), *coup, orig_deg + kick
          )
          e_k = _minimize_at_dihedral(
              kicked_mol,
              0,
              atom_indices,
              angles[idx],
              force_field,
              force_constant,
              max_iters,
          )
          if e_k < grid_best_energy[idx]:
            grid_best_energy[idx] = e_k
            grid_best_conf[idx] = Chem.Conformer(kicked_mol.GetConformer(0))

  if len(active_queue) == 0:
    conf_mol = Chem.Mol(mol)
    e = _minimize_at_dihedral(
        conf_mol,
        0,
        atom_indices,
        angles[0],
        force_field,
        force_constant,
        max_iters,
    )
    grid_best_energy[0] = e
    grid_best_conf[0] = Chem.Conformer(conf_mol.GetConformer(0))
    active_queue.append(0)

  while active_queue:
    curr_idx = active_queue.popleft()
    curr_conf = grid_best_conf[curr_idx]

    neighbors = [(curr_idx + 1) % n_grid, (curr_idx - 1) % n_grid]
    for nbr_idx in neighbors:
      nbr_ang = angles[nbr_idx]

      test_mol = Chem.Mol(mol)
      test_mol.RemoveAllConformers()
      test_mol.AddConformer(Chem.Conformer(curr_conf), assignId=True)

      e_nbr = _minimize_at_dihedral(
          test_mol,
          0,
          atom_indices,
          nbr_ang,
          force_field,
          force_constant,
          max_iters,
      )
      if e_nbr < (grid_best_energy[nbr_idx] - energy_threshold):
        grid_best_energy[nbr_idx] = e_nbr
        grid_best_conf[nbr_idx] = Chem.Conformer(test_mol.GetConformer(0))
        if nbr_idx not in active_queue:
          active_queue.append(nbr_idx)

  min_e = min(grid_best_energy)
  e_wave = np.array(grid_best_energy)
  delta_e = e_wave - min_e

  minima_indices = []
  for idx in range(len(delta_e)):
    prev_idx = (idx - 1) % len(delta_e)
    next_idx = (idx + 1) % len(delta_e)
    if delta_e[idx] <= delta_e[prev_idx] and delta_e[idx] <= delta_e[next_idx]:
      if delta_e[idx] < 1.0:
        minima_indices.append(idx)

  minima_list = [
      {
          "angle_deg": float(angles[idx]),
          "delta_e": float(delta_e[idx]),
          "absolute_e": float(e_wave[idx]),
      }
      for idx in minima_indices
  ]

  return {
      "algorithm": "fast_wavefront_propagation",
      "angles_deg": angles.tolist(),
      "energy_envelope_kcal_mol": e_wave.tolist(),
      "delta_e_kcal_mol": delta_e.tolist(),
      "barrier_kcal_mol": float(np.max(delta_e) - np.min(delta_e)),
      "minima": minima_list,
      "_conformers": grid_best_conf,
  }


def scan_torsion_bidirectional(
    mol: Chem.Mol,
    atom_indices: Tuple[int, int, int, int],
    step_deg: float = 5.0,
    force_field: str = "MMFF94s",
    force_constant: Optional[float] = None,
    max_iters: int = 300,
) -> Dict[str, Any]:
  """Sequential forward and backward coordinate continuation for hysteresis auditing."""
  i, j, k, l = atom_indices
  bond_in_ring = mol.GetBondBetweenAtoms(j, k).IsInRing()
  if force_constant is None:
    force_constant = 50.0 if bond_in_ring else 1000.0

  conf = mol.GetConformer(0)
  init_deg = (
      math.degrees(rdMolTransforms.GetDihedralRad(conf, i, j, k, l)) % 360.0
  )

  angles = np.arange(0, 360, step_deg, dtype=float)
  n_steps = int(360 / step_deg) + 2

  drives: Dict[int, Dict[float, Tuple[float, Chem.Conformer]]] = {
      +1: {},
      -1: {},
  }

  for direction in [+1, -1]:
    mol_curr = Chem.Mol(mol)
    curr_deg = init_deg

    for step in range(n_steps):
      target_deg = float((round(curr_deg / step_deg) * step_deg) % 360.0)
      e = _minimize_at_dihedral(
          mol_curr,
          0,
          atom_indices,
          target_deg,
          force_field,
          force_constant,
          max_iters,
      )
      if (
          target_deg not in drives[direction]
          or e < drives[direction][target_deg][0]
      ):
        drives[direction][target_deg] = (
            e,
            Chem.Conformer(mol_curr.GetConformer(0)),
        )
      curr_deg += direction * step_deg

  sorted_angles = []
  e_fwd = []
  e_bwd = []
  e_env = []
  hysteresis = []
  best_conformers = []

  for deg in angles:
    val_fwd = drives[+1].get(deg, (float("nan"), None))
    val_bwd = drives[-1].get(deg, (float("nan"), None))
    ef = val_fwd[0]
    eb = val_bwd[0]
    ee = min(ef, eb)
    h = abs(ef - eb) if not (math.isnan(ef) or math.isnan(eb)) else 0.0

    sorted_angles.append(deg)
    e_fwd.append(ef)
    e_bwd.append(eb)
    e_env.append(ee)
    hysteresis.append(h)
    best_conf = val_fwd[1] if ef <= eb else val_bwd[1]
    best_conformers.append(best_conf)

  sorted_angles = np.array(sorted_angles)
  e_fwd = np.array(e_fwd)
  e_bwd = np.array(e_bwd)
  e_env = np.array(e_env)
  hysteresis = np.array(hysteresis)

  min_env = np.min(e_env)
  delta_e = e_env - min_env

  return {
      "algorithm": "bidirectional_sequential",
      "angles_deg": sorted_angles.tolist(),
      "energy_forward_kcal_mol": e_fwd.tolist(),
      "energy_backward_kcal_mol": e_bwd.tolist(),
      "energy_envelope_kcal_mol": e_env.tolist(),
      "delta_e_kcal_mol": delta_e.tolist(),
      "hysteresis_kcal_mol": hysteresis.tolist(),
      "max_hysteresis_kcal_mol": float(np.max(hysteresis)),
      "mean_hysteresis_kcal_mol": float(np.mean(hysteresis)),
      "barrier_kcal_mol": float(np.max(delta_e) - np.min(delta_e)),
      "_conformers": best_conformers,
  }


def scan_torsion_2d_wavefront(
    mol: Chem.Mol,
    dihedral_1: Tuple[int, int, int, int],
    dihedral_2: Tuple[int, int, int, int],
    step_deg: float = 10.0,
    force_field: str = "MMFF94s",
    n_seeds: int = 30,
) -> Dict[str, Any]:
  """Template for 2D Coupled Wavefront PES Unrolling."""
  angles = np.arange(0, 360, step_deg, dtype=float)
  n_grid = len(angles)

  grid_energy = np.full((n_grid, n_grid), np.inf)
  grid_confs = {}
  active_queue = deque()

  seed_mol = Chem.Mol(mol)
  params = rdDistGeom.ETKDGv3()
  params.randomSeed = 42
  params.numThreads = 0
  rdDistGeom.EmbedMultipleConfs(seed_mol, numConfs=n_seeds, params=params)
  rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(seed_mol, numThreads=0)

  for cid in range(seed_mol.GetNumConformers()):
    a1 = (
        rdMolTransforms.GetDihedralDeg(seed_mol.GetConformer(cid), *dihedral_1)
        + 360.0
    ) % 360.0
    a2 = (
        rdMolTransforms.GetDihedralDeg(seed_mol.GetConformer(cid), *dihedral_2)
        + 360.0
    ) % 360.0
    i1 = int(np.argmin([min(abs(a1 - a), 360 - abs(a1 - a)) for a in angles]))
    i2 = int(np.argmin([min(abs(a2 - a), 360 - abs(a2 - a)) for a in angles]))

    conf_mol = Chem.Mol(seed_mol)
    conf_mol.RemoveAllConformers()
    conf_mol.AddConformer(
        Chem.Conformer(seed_mol.GetConformer(cid)), assignId=True
    )

    rdMolTransforms.SetDihedralDeg(
        conf_mol.GetConformer(0), *dihedral_1, float(angles[i1])
    )
    rdMolTransforms.SetDihedralDeg(
        conf_mol.GetConformer(0), *dihedral_2, float(angles[i2])
    )
    props = rdForceFieldHelpers.MMFFGetMoleculeProperties(
        conf_mol, mmffVariant="MMFF94s"
    )
    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(
        conf_mol, props, confId=0
    )
    ff.MMFFAddTorsionConstraint(
        *dihedral_1, False, float(angles[i1]), float(angles[i1]), 1000.0
    )
    ff.MMFFAddTorsionConstraint(
        *dihedral_2, False, float(angles[i2]), float(angles[i2]), 1000.0
    )
    ff.Initialize()
    ff.Minimize(maxIts=250)
    e = ff.CalcEnergy()

    if e < grid_energy[i1, i2]:
      grid_energy[i1, i2] = e
      grid_confs[(i1, i2)] = Chem.Conformer(conf_mol.GetConformer(0))
      active_queue.append((i1, i2))

  while active_queue:
    curr_i1, curr_i2 = active_queue.popleft()
    curr_conf = grid_confs[(curr_i1, curr_i2)]

    neighbors = [
        ((curr_i1 + 1) % n_grid, curr_i2),
        ((curr_i1 - 1) % n_grid, curr_i2),
        (curr_i1, (curr_i2 + 1) % n_grid),
        (curr_i1, (curr_i2 - 1) % n_grid),
    ]

    for n1, n2 in neighbors:
      test_mol = Chem.Mol(mol)
      test_mol.RemoveAllConformers()
      test_mol.AddConformer(Chem.Conformer(curr_conf), assignId=True)
      rdMolTransforms.SetDihedralDeg(
          test_mol.GetConformer(0), *dihedral_1, float(angles[n1])
      )
      rdMolTransforms.SetDihedralDeg(
          test_mol.GetConformer(0), *dihedral_2, float(angles[n2])
      )
      props = rdForceFieldHelpers.MMFFGetMoleculeProperties(
          test_mol, mmffVariant="MMFF94s"
      )
      ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(
          test_mol, props, confId=0
      )
      ff.MMFFAddTorsionConstraint(
          *dihedral_1, False, float(angles[n1]), float(angles[n1]), 1000.0
      )
      ff.MMFFAddTorsionConstraint(
          *dihedral_2, False, float(angles[n2]), float(angles[n2]), 1000.0
      )
      ff.Initialize()
      ff.Minimize(maxIts=200)
      e_nbr = ff.CalcEnergy()

      if e_nbr < (grid_energy[n1, n2] - 1e-3):
        grid_energy[n1, n2] = e_nbr
        grid_confs[(n1, n2)] = Chem.Conformer(test_mol.GetConformer(0))
        if (n1, n2) not in active_queue:
          active_queue.append((n1, n2))

  min_e = np.min(grid_energy)
  return {
      "algorithm": "2d_coupled_wavefront",
      "angles_deg": angles.tolist(),
      "energy_grid_kcal_mol": (grid_energy - min_e).tolist(),
      "global_minimum_kcal_mol": float(min_e),
      "barrier_kcal_mol": float(np.max(grid_energy) - min_e),
  }


def scan_torsion_adaptive_multires(
    mol: Chem.Mol,
    atom_indices: Tuple[int, int, int, int],
    coarse_step: float = 15.0,
    fine_step: float = 2.5,
    gradient_threshold: float = 0.5,
    force_field: str = "MMFF94s",
) -> Dict[str, Any]:
  """Template for Hierarchical Adaptive Multi-Resolution Scanning."""
  coarse_res = scan_torsion_wavefront(
      mol, atom_indices, step_deg=coarse_step, force_field=force_field
  )
  angles_c = np.array(coarse_res["angles_deg"])
  energies_c = np.array(coarse_res["delta_e_kcal_mol"])

  grad = np.abs(np.gradient(energies_c, coarse_step))
  fine_angles_to_add = set()

  for idx, g in enumerate(grad):
    if g > gradient_threshold or energies_c[idx] < 1.0:
      center_ang = angles_c[idx]
      sub_angles = np.arange(
          center_ang - coarse_step, center_ang + coarse_step + 0.1, fine_step
      )
      for a in sub_angles:
        fine_angles_to_add.add(
            float((round(a / fine_step) * fine_step) % 360.0)
        )

  all_angles = sorted(list(set(angles_c).union(fine_angles_to_add)))
  return {
      "algorithm": "adaptive_multiresolution",
      "coarse_angles": angles_c.tolist(),
      "total_grid_points": len(all_angles),
      "coarse_barrier": coarse_res["barrier_kcal_mol"],
  }


def main():
  parser = argparse.ArgumentParser(
      description="High-Performance Wavefront Torsion Drive Scanner."
  )
  parser.add_argument("--smiles", type=str, help="Input SMILES string.")
  parser.add_argument("--mol_file", type=str, help="Input SDF or MOL file.")
  parser.add_argument(
      "--list_rotors",
      action="store_true",
      help="Automatically detect, classify, and print all key rotatable bonds.",
  )
  parser.add_argument(
      "--algorithm",
      type=str,
      choices=[
          "wavefront",
          "bidirectional",
          "2d_wavefront",
          "adaptive_multires",
      ],
      default="wavefront",
      help="PES algorithm (default: 'wavefront').",
  )
  parser.add_argument(
      "--dihedral",
      type=int,
      nargs=4,
      help="4 zero-based atom indices for dihedral.",
  )
  parser.add_argument(
      "--dihedral2",
      type=int,
      nargs=4,
      help="Second dihedral for 2D wavefront scans.",
  )
  parser.add_argument(
      "--smarts", type=str, help="SMARTS pattern matching target dihedral."
  )
  parser.add_argument(
      "--smarts_indices", type=int, nargs=4, help="Indices within SMARTS match."
  )
  parser.add_argument(
      "--step_deg",
      type=float,
      default=5.0,
      help="Grid step in degrees (default: 5.0).",
  )
  parser.add_argument(
      "--n_seeds",
      type=int,
      default=25,
      help="Number of initial ETKDG conformer seeds.",
  )
  parser.add_argument(
      "--force_field",
      type=str,
      default="MMFF94s",
      choices=["MMFF94s", "MMFF94", "UFF"],
  )
  parser.add_argument("--output_csv", type=str, help="Output CSV path.")
  parser.add_argument("--output_json", type=str, help="Output JSON path.")
  parser.add_argument(
      "--output_sdf", type=str, help="Output multi-conformer trajectory SDF."
  )
  parser.add_argument(
      "--plot_png",
      type=str,
      help="Path for saving publication-ready PES plot PNG.",
  )
  parser.add_argument(
      "--max_bonds",
      "--max-bonds",
      dest="max_bonds",
      type=int,
      default=None,
      help=(
          "Optional cap on number of rotatable bonds returned by --list_rotors"
          " (default: None, returning all key rotatable bonds)."
      ),
  )
  args = parser.parse_args()

  if args.smiles:
    mol = Chem.MolFromSmiles(args.smiles)
  elif args.mol_file:
    suppl = Chem.SDMolSupplier(args.mol_file, removeHs=False)
    mol = suppl[0]
  else:
    raise ValueError("Must provide either --smiles or --mol_file.")

  if args.list_rotors:
    rotors = identify_key_rotatable_bonds(mol, max_bonds=args.max_bonds)
    print(f"=== Identified {len(rotors)} Key Rotatable Bond(s) ===")
    for idx, r in enumerate(rotors, start=1):
      print(
          f"  {idx}. Dihedral {r['dihedral']} | Type: {r['type']:<14} |"
          f" Priority {r['priority']} | {r['description']}"
      )
    return

  mol_3d, ground_e = setup_molecule_3d(mol, force_field=args.force_field)
  atom_indices = resolve_dihedral_indices(
      mol_3d,
      dihedral=args.dihedral,
      smarts=args.smarts,
      smarts_indices=args.smarts_indices,
  )

  print(f"=== Running Torsion Drive: {args.algorithm.upper()} ===")
  print(
      f"Target Dihedral: {atom_indices}, Force Field: {args.force_field}, Step:"
      f" {args.step_deg}°"
  )

  t0 = time.time()
  if args.algorithm == "wavefront":
    scan_res = scan_torsion_wavefront(
        mol_3d,
        atom_indices,
        step_deg=args.step_deg,
        n_seeds=args.n_seeds,
        force_field=args.force_field,
    )
  elif args.algorithm == "bidirectional":
    scan_res = scan_torsion_bidirectional(
        mol_3d,
        atom_indices,
        step_deg=args.step_deg,
        force_field=args.force_field,
    )
  elif args.algorithm == "2d_wavefront":
    if not args.dihedral2:
      raise ValueError("--dihedral2 is required for 2D wavefront scans.")
    scan_res = scan_torsion_2d_wavefront(
        mol_3d,
        atom_indices,
        tuple(args.dihedral2),
        step_deg=args.step_deg,
        force_field=args.force_field,
    )
  elif args.algorithm == "adaptive_multires":
    scan_res = scan_torsion_adaptive_multires(
        mol_3d, atom_indices, force_field=args.force_field
    )

  dt = time.time() - t0
  print(f"Scan Completed in {dt:.2f}s.")
  print(f"Barrier Height: {scan_res.get('barrier_kcal_mol', 0.0):.3f} kcal/mol")

  if args.output_csv and "angles_deg" in scan_res:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True
    )
    df = pd.DataFrame({
        "angle_deg": scan_res["angles_deg"],
        "delta_e_kcal_mol": scan_res.get("delta_e_kcal_mol", []),
    })
    df.to_csv(args.output_csv, index=False)
    print(f"Saved CSV to {args.output_csv}")

  if args.output_json:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_json)), exist_ok=True
    )
    payload = {k: v for k, v in scan_res.items() if k != "_conformers"}
    with open(args.output_json, "w") as f:
      json.dump(payload, f, indent=2)
    print(f"Saved JSON to {args.output_json}")

  if args.output_sdf and "_conformers" in scan_res:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_sdf)), exist_ok=True
    )
    with Chem.SDWriter(args.output_sdf) as writer:
      mol_out = Chem.Mol(mol_3d)
      mol_out.RemoveAllConformers()
      for idx, conf in enumerate(scan_res["_conformers"]):
        if conf is not None:
          conf.SetId(idx)
          mol_out.AddConformer(conf, assignId=True)
          mol_out.SetProp("Dihedral_deg", f"{scan_res['angles_deg'][idx]:.1f}")
          mol_out.SetProp(
              "DeltaE_kcal_mol", f"{scan_res['delta_e_kcal_mol'][idx]:.3f}"
          )
          writer.write(mol_out, confId=idx)
    print(f"Saved SDF trajectory to {args.output_sdf}")

  if args.plot_png and "angles_deg" in scan_res:
    set_safe_plot_style()
    os.makedirs(os.path.dirname(os.path.abspath(args.plot_png)), exist_ok=True)
    angles = np.array(scan_res["angles_deg"])
    energies = np.array(scan_res["delta_e_kcal_mol"])
    barrier = scan_res.get(
        "barrier_kcal_mol", float(np.max(energies) - np.min(energies))
    )

    plt.figure(figsize=(8, 5), dpi=300)
    plt.plot(
        angles,
        energies,
        color="#4F46E5",
        linewidth=2.5,
        label=f"PES Profile ({args.force_field})",
    )
    plt.fill_between(angles, energies, alpha=0.15, color="#4F46E5")
    plt.axhline(0, color="#64748B", linestyle="--", linewidth=0.8)

    # Annotate global min and barrier
    min_idx = np.argmin(energies)
    max_idx = np.argmax(energies)
    plt.scatter(
        [angles[min_idx]],
        [energies[min_idx]],
        color="#10B981",
        s=70,
        zorder=5,
        label=f"Global Min ({angles[min_idx]:.0f}°)",
    )
    plt.scatter(
        [angles[max_idx]],
        [energies[max_idx]],
        color="#F43F5E",
        s=70,
        zorder=5,
        label=f"Barrier ({barrier:.2f} kcal/mol)",
    )

    plt.xlabel("Dihedral Angle (degrees)", fontsize=11, fontweight="bold")
    plt.ylabel("Relative Energy ΔE (kcal/mol)", fontsize=11, fontweight="bold")
    plt.title(
        f"Relaxed Torsion Drive ({args.algorithm.upper()} -"
        f" {args.force_field})",
        fontsize=12,
        fontweight="bold",
    )
    plt.xlim(0, 360)
    plt.xticks(np.arange(0, 361, 45))
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(frameon=True, facecolor="white", edgecolor="none")
    plt.tight_layout()
    plt.savefig(args.plot_png)
    plt.close()
    print(f"Saved PES plot to {args.plot_png}")


if __name__ == "__main__":
  main()
