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

"""Quantitative Ligand Steric & Electronic Descriptors Engine.

Calculates organometallic and coordination descriptors directly from 3D
conformers:
1. SambVca Buried Volume (%V_bur) & Quadrant/Octant Steric Distribution.
2. Topographic Steric Contour Maps for coordination complexes.
3. Tolman Cone Angle (θ) & Exact Minimum Enclosing Cone Angles.
4. Sterimol Parameters (L, B1, B5) along primary metal-donor coordination
vector.
5. Bidentate Chelate Bite Angles (βn) & Bite Distances for Bisphosphines and
Bipyridines.
6. Electronic Proxy Descriptors (TEP / ν_CO correlation, donor polarizability).
"""

import argparse
import json
import math
import os
import sys
from typing import Any, Dict, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

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

# Bondi Van der Waals Radii (in Angstroms)
BONDI_VDW: Dict[int, float] = {
    1: 1.20,  # H
    3: 1.82,  # Li
    5: 2.00,  # B
    6: 1.70,  # C
    7: 1.55,  # N
    8: 1.52,  # O
    9: 1.47,  # F
    11: 2.27,  # Na
    12: 1.73,  # Mg
    13: 2.00,  # Al
    14: 2.10,  # Si
    15: 1.80,  # P
    16: 1.80,  # S
    17: 1.75,  # Cl
    19: 2.75,  # K
    26: 2.00,  # Fe
    27: 2.00,  # Co
    28: 1.63,  # Ni
    29: 1.40,  # Cu
    30: 1.39,  # Zn
    35: 1.85,  # Br
    44: 2.00,  # Ru
    45: 2.00,  # Rh
    46: 1.63,  # Pd
    47: 1.72,  # Ag
    53: 1.98,  # I
    77: 2.00,  # Ir
    78: 1.75,  # Pt
    79: 1.66,  # Au
}


def find_donor_atom(mol: Chem.Mol, metal_idx: Optional[int] = None) -> int:
  """Auto-detects the most likely coordination donor atom (P, N, S, O, or Carbene C)."""
  if metal_idx is not None and metal_idx < mol.GetNumAtoms():
    matom = mol.GetAtomWithIdx(metal_idx)
    for atomic_num in (15, 7, 16, 6, 8, 17, 35, 53, 9):
      for nbr in matom.GetNeighbors():
        if nbr.GetAtomicNum() == atomic_num:
          return nbr.GetIdx()

  for atomic_num in (15, 7, 16, 6, 8):
    for atom in mol.GetAtoms():
      if metal_idx is not None and atom.GetIdx() == metal_idx:
        continue
      if atom.GetAtomicNum() == atomic_num:
        if atomic_num == 6 and atom.GetDegree() > 2:
          continue
        return atom.GetIdx()
  return 0


_TRANSITION_METALS = (
    set(range(21, 31)) | set(range(39, 49)) | set(range(57, 81))
)


def _normalize_dative_bonds(mol: Chem.Mol) -> Chem.Mol:
  """Converts covalent M-L single bonds on saturated neutral L-type donors to dative L->M bonds."""
  rw = Chem.RWMol(mol)
  to_replace = []
  for bond in rw.GetBonds():
    if bond.GetBondType() != Chem.BondType.SINGLE:
      continue
    a1 = bond.GetBeginAtom()
    a2 = bond.GetEndAtom()
    z1, z2 = a1.GetAtomicNum(), a2.GetAtomicNum()
    if z1 in _TRANSITION_METALS and z2 in (7, 8, 15, 16):
      metal_a, donor_a = a1, a2
    elif z2 in _TRANSITION_METALS and z1 in (7, 8, 15, 16):
      metal_a, donor_a = a2, a1
    else:
      continue
    if donor_a.GetFormalCharge() != 0:
      continue
    dz = donor_a.GetAtomicNum()
    deg = donor_a.GetDegree()
    is_sat_l = (
        dz in (7, 15) and (deg >= 4 or (donor_a.GetIsAromatic() and deg >= 3))
    ) or (dz in (8, 16) and deg >= 3)
    if is_sat_l:
      to_replace.append((donor_a.GetIdx(), metal_a.GetIdx()))
  for d_idx, m_idx in to_replace:
    rw.RemoveBond(d_idx, m_idx)
    rw.AddBond(d_idx, m_idx, Chem.BondType.DATIVE)
  out = rw.GetMol()
  out.UpdatePropertyCache(strict=False)
  Chem.SanitizeMol(out, catchErrors=True)
  Chem.GetSymmSSSR(out)
  return out


def _get_ligand_atom_indices(
    mol_3d: Chem.Mol,
    donor_idx: int,
    metal_idx: Optional[int] = None,
    donor_idx_2: Optional[int] = None,
) -> set[int]:
  """Returns atom indices belonging to the ligand fragment of donor_idx (and donor_idx_2) without crossing metal_idx."""
  if metal_idx is None:
    return {a.GetIdx() for a in mol_3d.GetAtoms() if a.GetIdx() != metal_idx}
  visited: set[int] = set()
  queue = [donor_idx]
  if donor_idx_2 is not None:
    queue.append(donor_idx_2)
  while queue:
    curr = queue.pop(0)
    if curr in visited or curr == metal_idx:
      continue
    visited.add(curr)
    for nbr in mol_3d.GetAtomWithIdx(curr).GetNeighbors():
      n_idx = nbr.GetIdx()
      if n_idx not in visited and n_idx != metal_idx:
        queue.append(n_idx)
  return visited


def _compute_bidentate_frame(
    mol_3d: Chem.Mol,
    donor_idx_1: int,
    donor_idx_2: int,
    metal_dist: float = 2.28,
) -> Tuple[np.ndarray, np.ndarray]:
  """Returns (metal_pos, z_axis) for a bidentate chelate with donor_idx_1 and donor_idx_2."""
  conf = mol_3d.GetConformer()
  d1_pos = np.array(conf.GetAtomPosition(donor_idx_1))
  d2_pos = np.array(conf.GetAtomPosition(donor_idx_2))
  bite_dist = float(np.linalg.norm(d1_pos - d2_pos))
  mid_pos = 0.5 * (d1_pos + d2_pos)
  d_half = 0.5 * bite_dist
  h = math.sqrt(metal_dist**2 - d_half**2) if metal_dist > d_half else 1.0

  nbrs_1 = [
      np.array(conf.GetAtomPosition(n.GetIdx()))
      for n in mol_3d.GetAtomWithIdx(donor_idx_1).GetNeighbors()
  ]
  nbrs_2 = [
      np.array(conf.GetAtomPosition(n.GetIdx()))
      for n in mol_3d.GetAtomWithIdx(donor_idx_2).GetNeighbors()
  ]
  sub_centroid = (
      np.mean(nbrs_1 + nbrs_2, axis=0) if (nbrs_1 + nbrs_2) else mid_pos
  )
  perp_vec = mid_pos - sub_centroid
  norm = np.linalg.norm(perp_vec)
  perp_vec = perp_vec / norm if norm > 1e-4 else np.array([0.0, 0.0, 1.0])
  metal_pos = mid_pos + perp_vec * h
  z_vec = mid_pos - metal_pos
  z_norm = np.linalg.norm(z_vec)
  z_axis = z_vec / z_norm if z_norm > 1e-4 else -perp_vec
  return metal_pos, z_axis


def compute_buried_volume(
    mol_3d: Chem.Mol,
    donor_idx: Optional[int] = None,
    metal_idx: Optional[int] = None,
    metal_dist: float = 2.28,
    sphere_radius: float = 3.5,
    vdw_scale: float = 1.17,
    mesh_spacing: float = 0.1,
    include_hs: bool = True,
    donor_idx_2: Optional[int] = None,
) -> Dict[str, Any]:
  """Computes the SambVca %V_bur (Percent Buried Volume) and quadrant steric partitioning.

  Args:
      mol_3d: RDKit Mol object with 3D coordinates.
      donor_idx: Atom index of the coordinating donor atom.
      metal_idx: Atom index of the metal center if already in the molecule.
      metal_dist: Coordination bond distance in Å (used if metal_idx is None).
      sphere_radius: Coordination sphere radius in Å (default: 3.5 Å).
      vdw_scale: Bondi VdW radius scaling factor (default: 1.17 per Cavallo
        standard).
      mesh_spacing: Numerical voxel grid spacing in Å (default: 0.1 Å).
      include_hs: Whether to include hydrogen atoms in steric volume.
      donor_idx_2: Optional second donor atom index for bidentate chelates.
  """
  conf = mol_3d.GetConformer()

  if metal_idx is not None and metal_idx < mol_3d.GetNumAtoms():
    metal_pos = np.array(conf.GetAtomPosition(metal_idx))
    if donor_idx is not None and donor_idx != metal_idx:
      donor_pos = np.array(conf.GetAtomPosition(donor_idx))
      z_axis = donor_pos - metal_pos
      norm = np.linalg.norm(z_axis)
      z_axis = z_axis / norm if norm > 1e-4 else np.array([0.0, 0.0, 1.0])
    else:
      matom = mol_3d.GetAtomWithIdx(metal_idx)
      nbrs = [n.GetIdx() for n in matom.GetNeighbors()]
      if nbrs:
        p_nbrs = [
            n
            for n in nbrs
            if mol_3d.GetAtomWithIdx(n).GetAtomicNum() in (15, 7)
        ]
        ref_nbr = p_nbrs[0] if p_nbrs else nbrs[0]
        donor_pos = np.array(conf.GetAtomPosition(ref_nbr))
        z_axis = donor_pos - metal_pos
        norm = np.linalg.norm(z_axis)
        z_axis = z_axis / norm if norm > 1e-4 else np.array([0.0, 0.0, 1.0])
      else:
        donor_pos = metal_pos + np.array([0.0, 0.0, metal_dist])
        z_axis = np.array([0.0, 0.0, 1.0])

    matom = mol_3d.GetAtomWithIdx(metal_idx)
    other_nbrs = [
        n.GetIdx()
        for n in matom.GetNeighbors()
        if np.linalg.norm(
            np.array(conf.GetAtomPosition(n.GetIdx())) - metal_pos - z_axis
        )
        > 1e-2
    ]
    if other_nbrs:
      v_ref = np.array(conf.GetAtomPosition(other_nbrs[0])) - metal_pos
      v_perp = v_ref - np.dot(v_ref, z_axis) * z_axis
      norm_perp = np.linalg.norm(v_perp)
      x_axis = (
          v_perp / norm_perp
          if norm_perp > 1e-4
          else np.cross(
              z_axis,
              [1.0, 0.0, 0.0] if abs(z_axis[0]) < 0.9 else [0.0, 1.0, 0.0],
          )
      )
    else:
      ref_vec = (
          np.array([1.0, 0.0, 0.0])
          if abs(z_axis[0]) < 0.9
          else np.array([0.0, 1.0, 0.0])
      )
      x_axis = np.cross(z_axis, ref_vec)
    x_axis /= np.linalg.norm(x_axis)
    y_axis = np.cross(z_axis, x_axis)

  elif donor_idx_2 is not None:
    if donor_idx is None:
      donor_idx = find_donor_atom(mol_3d)
    donor_pos = np.array(conf.GetAtomPosition(donor_idx))
    metal_pos, z_axis = _compute_bidentate_frame(
        mol_3d, donor_idx, donor_idx_2, metal_dist=metal_dist
    )
    d2_pos = np.array(conf.GetAtomPosition(donor_idx_2))
    v_d = d2_pos - donor_pos
    v_perp = v_d - np.dot(v_d, z_axis) * z_axis
    norm_perp = np.linalg.norm(v_perp)
    if norm_perp > 1e-4:
      x_axis = v_perp / norm_perp
    else:
      ref_vec = (
          np.array([1.0, 0.0, 0.0])
          if abs(z_axis[0]) < 0.9
          else np.array([0.0, 1.0, 0.0])
      )
      x_axis = np.cross(z_axis, ref_vec)
      x_axis /= np.linalg.norm(x_axis)
    y_axis = np.cross(z_axis, x_axis)

  else:
    if donor_idx is None:
      donor_idx = find_donor_atom(mol_3d)
    donor_pos = np.array(conf.GetAtomPosition(donor_idx))
    nbrs = [
        nbr.GetIdx() for nbr in mol_3d.GetAtomWithIdx(donor_idx).GetNeighbors()
    ]
    if nbrs:
      nbr_coords = [np.array(conf.GetAtomPosition(i)) for i in nbrs]
      centroid = np.mean(nbr_coords, axis=0)
      z_axis = centroid - donor_pos
      norm = np.linalg.norm(z_axis)
      z_axis = z_axis / norm if norm > 1e-4 else np.array([0.0, 0.0, 1.0])
    else:
      z_axis = np.array([0.0, 0.0, 1.0])

    metal_pos = donor_pos - z_axis * metal_dist
    ref_vec = (
        np.array([1.0, 0.0, 0.0])
        if abs(z_axis[0]) < 0.9
        else np.array([0.0, 1.0, 0.0])
    )
    x_axis = np.cross(z_axis, ref_vec)
    x_axis /= np.linalg.norm(x_axis)
    y_axis = np.cross(z_axis, x_axis)

  # Collect ligand atoms (excluding metal atom)
  atoms = []
  for atom in mol_3d.GetAtoms():
    idx = atom.GetIdx()
    if metal_idx is not None and idx == metal_idx:
      continue
    atomic_num = atom.GetAtomicNum()
    if not include_hs and atomic_num == 1:
      continue
    radius = BONDI_VDW.get(atomic_num, 1.70) * vdw_scale
    pos = np.array(conf.GetAtomPosition(idx)) - metal_pos
    lx = np.dot(pos, x_axis)
    ly = np.dot(pos, y_axis)
    lz = np.dot(pos, z_axis)
    atoms.append((np.array([lx, ly, lz]), radius))

  # Construct 3D sphere grid
  n_pts = int(math.ceil(sphere_radius / mesh_spacing))
  coords = np.linspace(-sphere_radius, sphere_radius, 2 * n_pts + 1)
  x, y, z = np.meshgrid(coords, coords, coords, indexing="ij")
  r2 = x**2 + y**2 + z**2
  sphere_mask = r2 <= (sphere_radius**2)
  total_sphere_pts = int(np.sum(sphere_mask))

  sphere_x = x[sphere_mask]
  sphere_y = y[sphere_mask]
  sphere_z = z[sphere_mask]
  grid_pts = np.column_stack([sphere_x, sphere_y, sphere_z])

  # Check voxel occupancy
  buried_mask = np.zeros(len(grid_pts), dtype=bool)
  for pos, radius in atoms:
    dist2 = np.sum((grid_pts - pos) ** 2, axis=1)
    buried_mask |= dist2 <= radius**2

  buried_pts = int(np.sum(buried_mask))
  percent_v_bur = (buried_pts / total_sphere_pts) * 100.0
  total_vol = (4.0 / 3.0) * math.pi * (sphere_radius**3)
  buried_vol = total_vol * (percent_v_bur / 100.0)

  # Quadrant breakdown (Q1: +x,+y, Q2: -x,+y, Q3: -x,-y, Q4: +x,-y)
  quadrants = {}
  for q_idx, (sx, sy) in enumerate([(1, 1), (-1, 1), (-1, -1), (1, -1)]):
    q_mask = (sphere_x * sx >= 0) & (sphere_y * sy >= 0)
    q_tot = int(np.sum(q_mask))
    q_bur = int(np.sum(buried_mask[q_mask]))
    q_pct = (q_bur / q_tot * 100.0) if q_tot > 0 else 0.0
    quadrants[f"Q{q_idx+1}"] = round(q_pct, 2)

  return {
      "percent_v_bur": round(percent_v_bur, 2),
      "total_sphere_volume_A3": round(total_vol, 2),
      "buried_volume_A3": round(buried_vol, 2),
      "free_volume_A3": round(total_vol - buried_vol, 2),
      "sphere_radius_A": sphere_radius,
      "metal_donor_dist_A": (
          metal_dist
          if metal_idx is None
          else round(float(np.linalg.norm(donor_pos - metal_pos)), 2)
          if donor_idx is not None
          else metal_dist
      ),
      "vdw_scale": vdw_scale,
      "quadrants": quadrants,
  }


def compute_tolman_cone_angle(
    mol_3d: Chem.Mol,
    donor_idx: int,
    metal_idx: Optional[int] = None,
    metal_dist: float = 2.28,
    donor_idx_2: Optional[int] = None,
) -> Dict[str, Any]:
  """Computes exact and classical Tolman cone angles (θ) from 3D conformer coordinates."""
  conf = mol_3d.GetConformer()
  donor_pos = np.array(conf.GetAtomPosition(donor_idx))

  if metal_idx is not None and metal_idx < mol_3d.GetNumAtoms():
    metal_pos = np.array(conf.GetAtomPosition(metal_idx))
    if donor_idx_2 is not None:
      d2_pos = np.array(conf.GetAtomPosition(donor_idx_2))
      z_axis = 0.5 * (donor_pos + d2_pos) - metal_pos
    else:
      z_axis = donor_pos - metal_pos
    norm = np.linalg.norm(z_axis)
    z_axis = z_axis / norm if norm > 1e-4 else np.array([0.0, 0.0, 1.0])
  elif donor_idx_2 is not None:
    metal_pos, z_axis = _compute_bidentate_frame(
        mol_3d, donor_idx, donor_idx_2, metal_dist=metal_dist
    )
  else:
    nbrs = [
        nbr.GetIdx() for nbr in mol_3d.GetAtomWithIdx(donor_idx).GetNeighbors()
    ]
    if nbrs:
      nbr_coords = [np.array(conf.GetAtomPosition(i)) for i in nbrs]
      centroid = np.mean(nbr_coords, axis=0)
      z_axis = centroid - donor_pos
      norm = np.linalg.norm(z_axis)
      z_axis = z_axis / norm if norm > 1e-4 else np.array([0.0, 0.0, 1.0])
    else:
      z_axis = np.array([0.0, 0.0, 1.0])
    metal_pos = donor_pos - z_axis * metal_dist

  lig_indices = _get_ligand_atom_indices(
      mol_3d, donor_idx, metal_idx=metal_idx, donor_idx_2=donor_idx_2
  )

  # Compute semi-angles for every substituent atom in the ligand fragment
  sub_angles = []
  max_cone_angle = 0.0

  for atom in mol_3d.GetAtoms():
    idx = atom.GetIdx()
    if (
        idx not in lig_indices
        or idx == donor_idx
        or idx == donor_idx_2
        or (metal_idx is not None and idx == metal_idx)
    ):
      continue
    atomic_num = atom.GetAtomicNum()
    vdw_r = BONDI_VDW.get(atomic_num, 1.70)
    pos = np.array(conf.GetAtomPosition(idx)) - metal_pos
    dist = np.linalg.norm(pos)
    if dist < 1e-3:
      continue

    cos_alpha = np.dot(pos, z_axis) / dist
    alpha = math.acos(np.clip(cos_alpha, -1.0, 1.0))
    sin_gamma = min(1.0, vdw_r / dist)
    gamma = math.asin(sin_gamma)

    semi_angle_deg = min(180.0, math.degrees(alpha + gamma))
    full_cone_deg = min(360.0, 2.0 * semi_angle_deg)
    sub_angles.append(full_cone_deg)
    if full_cone_deg > max_cone_angle:
      max_cone_angle = full_cone_deg

  sub_angles.sort(reverse=True)
  top_3_mean = (
      float(np.mean(sub_angles[:3])) if len(sub_angles) >= 3 else max_cone_angle
  )

  return {
      "exact_cone_angle_deg": round(max_cone_angle, 1),
      "tolman_cone_angle_deg": round(top_3_mean, 1),
      "metal_donor_dist_A": (
          metal_dist
          if metal_idx is None
          else round(float(np.linalg.norm(donor_pos - metal_pos)), 2)
      ),
  }


def compute_sterimol_parameters(
    mol_3d: Chem.Mol,
    donor_idx: int,
    metal_idx: Optional[int] = None,
    metal_dist: float = 2.28,
    donor_idx_2: Optional[int] = None,
) -> Dict[str, float]:
  """Calculates Verloop Sterimol parameters (L, B1, B5) along the primary coordination vector."""
  conf = mol_3d.GetConformer()
  donor_pos = np.array(conf.GetAtomPosition(donor_idx))

  if metal_idx is not None and metal_idx < mol_3d.GetNumAtoms():
    metal_pos = np.array(conf.GetAtomPosition(metal_idx))
    if donor_idx_2 is not None:
      d2_pos = np.array(conf.GetAtomPosition(donor_idx_2))
      z_axis = 0.5 * (donor_pos + d2_pos) - metal_pos
    else:
      z_axis = donor_pos - metal_pos
    norm = np.linalg.norm(z_axis)
    z_axis = z_axis / norm if norm > 1e-4 else np.array([0.0, 0.0, 1.0])
  elif donor_idx_2 is not None:
    metal_pos, z_axis = _compute_bidentate_frame(
        mol_3d, donor_idx, donor_idx_2, metal_dist=metal_dist
    )
  else:
    nbrs = [
        nbr.GetIdx() for nbr in mol_3d.GetAtomWithIdx(donor_idx).GetNeighbors()
    ]
    if nbrs:
      nbr_coords = [np.array(conf.GetAtomPosition(i)) for i in nbrs]
      centroid = np.mean(nbr_coords, axis=0)
      z_axis = centroid - donor_pos
      norm = np.linalg.norm(z_axis)
      z_axis = z_axis / norm if norm > 1e-4 else np.array([0.0, 0.0, 1.0])
    else:
      z_axis = np.array([0.0, 0.0, 1.0])
    metal_pos = donor_pos - z_axis * metal_dist

  if abs(z_axis[0]) < 0.9:
    x_axis = np.cross(z_axis, [1.0, 0.0, 0.0])
  else:
    x_axis = np.cross(z_axis, [0.0, 1.0, 0.0])
  x_axis /= np.linalg.norm(x_axis)
  y_axis = np.cross(z_axis, x_axis)

  lig_indices = _get_ligand_atom_indices(
      mol_3d, donor_idx, metal_idx=metal_idx, donor_idx_2=donor_idx_2
  )

  l_max = 0.0
  b5_max = 0.0
  profile_per_bin = np.zeros(36)

  for atom in mol_3d.GetAtoms():
    idx = atom.GetIdx()
    if idx not in lig_indices or (metal_idx is not None and idx == metal_idx):
      continue
    atomic_num = atom.GetAtomicNum()
    vdw_r = BONDI_VDW.get(atomic_num, 1.70)
    pos = np.array(conf.GetAtomPosition(idx)) - metal_pos

    z_val = np.dot(pos, z_axis)
    x_val = np.dot(pos, x_axis)
    y_val = np.dot(pos, y_axis)

    atom_l = z_val + vdw_r
    if atom_l > l_max:
      l_max = atom_l

    r_perp = math.sqrt(x_val**2 + y_val**2) + vdw_r
    if r_perp > b5_max:
      b5_max = r_perp

    phi = math.atan2(y_val, x_val) % (2 * math.pi)
    bin_idx = int(phi / (2 * math.pi / 36)) % 36
    if r_perp > profile_per_bin[bin_idx]:
      profile_per_bin[bin_idx] = r_perp

  b1_min = (
      float(np.min(profile_per_bin[profile_per_bin > 0]))
      if np.any(profile_per_bin > 0)
      else 1.0
  )

  return {
      "sterimol_L": round(l_max, 2),
      "sterimol_B1": round(b1_min, 2),
      "sterimol_B5": round(b5_max, 2),
  }


def compute_chelate_bite_angle(
    mol_3d: Chem.Mol,
    donor_idx_1: int,
    donor_idx_2: int,
    metal_idx: Optional[int] = None,
    metal_dist: float = 2.28,
) -> Dict[str, Any]:
  """Computes natural bite angle (βn) and bite distance (d_bite) for bidentate chelates."""
  conf = mol_3d.GetConformer()
  d1_pos = np.array(conf.GetAtomPosition(donor_idx_1))
  d2_pos = np.array(conf.GetAtomPosition(donor_idx_2))

  bite_dist = float(np.linalg.norm(d1_pos - d2_pos))

  if metal_idx is not None and metal_idx < mol_3d.GetNumAtoms():
    m_pos = np.array(conf.GetAtomPosition(metal_idx))
  else:
    mid_pos = (d1_pos + d2_pos) / 2.0
    d_half = bite_dist / 2.0
    h = math.sqrt(metal_dist**2 - d_half**2) if metal_dist > d_half else 1.0

    nbrs_1 = [
        np.array(conf.GetAtomPosition(n.GetIdx()))
        for n in mol_3d.GetAtomWithIdx(donor_idx_1).GetNeighbors()
    ]
    nbrs_2 = [
        np.array(conf.GetAtomPosition(n.GetIdx()))
        for n in mol_3d.GetAtomWithIdx(donor_idx_2).GetNeighbors()
    ]
    sub_centroid = (
        np.mean(nbrs_1 + nbrs_2, axis=0) if (nbrs_1 + nbrs_2) else mid_pos
    )
    perp_vec = mid_pos - sub_centroid
    norm = np.linalg.norm(perp_vec)
    perp_vec = perp_vec / norm if norm > 1e-4 else np.array([0.0, 0.0, 1.0])
    m_pos = mid_pos + perp_vec * h

  v1 = d1_pos - m_pos
  v2 = d2_pos - m_pos
  cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
  bite_angle_deg = math.degrees(math.acos(np.clip(cos_angle, -1.0, 1.0)))

  return {
      "bite_angle_deg": round(bite_angle_deg, 1),
      "bite_distance_A": round(bite_dist, 2),
      "donor_indices": [donor_idx_1, donor_idx_2],
  }


def render_steric_contour_map(
    mol_3d: Chem.Mol,
    donor_idx: Optional[int] = None,
    metal_idx: Optional[int] = None,
    sphere_radius: float = 3.5,
    vdw_scale: float = 1.17,
    output_png: str = "steric_map.png",
    title: str = "Topographic Steric Map",
) -> str:
  """Renders a publication-quality 2D topographic steric elevation contour map."""
  set_safe_plot_style()
  conf = mol_3d.GetConformer()
  if metal_idx is not None and metal_idx < mol_3d.GetNumAtoms():
    metal_pos = np.array(conf.GetAtomPosition(metal_idx))
    if donor_idx is not None:
      z_axis = np.array(conf.GetAtomPosition(donor_idx)) - metal_pos
    else:
      nbrs = [
          n.GetIdx() for n in mol_3d.GetAtomWithIdx(metal_idx).GetNeighbors()
      ]
      z_axis = (
          np.array(conf.GetAtomPosition(nbrs[0])) - metal_pos
          if nbrs
          else np.array([0.0, 0.0, 1.0])
      )
    norm = np.linalg.norm(z_axis)
    z_axis = z_axis / norm if norm > 1e-4 else np.array([0.0, 0.0, 1.0])
  else:
    if donor_idx is None:
      donor_idx = find_donor_atom(mol_3d)
    donor_pos = np.array(conf.GetAtomPosition(donor_idx))
    nbrs = [
        nbr.GetIdx() for nbr in mol_3d.GetAtomWithIdx(donor_idx).GetNeighbors()
    ]
    if nbrs:
      centroid = np.mean(
          [np.array(conf.GetAtomPosition(i)) for i in nbrs], axis=0
      )
      z_axis = centroid - donor_pos
      norm = np.linalg.norm(z_axis)
      z_axis = z_axis / norm if norm > 1e-4 else np.array([0.0, 0.0, 1.0])
    else:
      z_axis = np.array([0.0, 0.0, 1.0])
    metal_pos = donor_pos - z_axis * 2.28

  ref_vec = (
      np.array([1.0, 0.0, 0.0])
      if abs(z_axis[0]) < 0.9
      else np.array([0.0, 1.0, 0.0])
  )
  x_axis = np.cross(z_axis, ref_vec)
  x_axis /= np.linalg.norm(x_axis)
  y_axis = np.cross(z_axis, x_axis)

  # 2D Grid
  grid_size = 100
  lin = np.linspace(-sphere_radius, sphere_radius, grid_size)
  grid_x, grid_y = np.meshgrid(lin, lin)
  z_elevation = np.zeros_like(grid_x)

  atoms = []
  for atom in mol_3d.GetAtoms():
    idx = atom.GetIdx()
    if metal_idx is not None and idx == metal_idx:
      continue
    atomic_num = atom.GetAtomicNum()
    radius = BONDI_VDW.get(atomic_num, 1.70) * vdw_scale
    pos = np.array(conf.GetAtomPosition(idx)) - metal_pos
    lx = np.dot(pos, x_axis)
    ly = np.dot(pos, y_axis)
    lz = np.dot(pos, z_axis)
    atoms.append((lx, ly, lz, radius))

  for i in range(grid_size):
    for j in range(grid_size):
      gx = grid_x[i, j]
      gy = grid_y[i, j]
      if gx**2 + gy**2 > sphere_radius**2:
        continue
      max_z = 0.0
      for lx, ly, lz, r in atoms:
        d_xy2 = (gx - lx) ** 2 + (gy - ly) ** 2
        if d_xy2 <= r**2:
          dz = math.sqrt(r**2 - d_xy2)
          top_z = lz + dz
          if top_z > max_z:
            max_z = top_z
      z_elevation[i, j] = min(sphere_radius, max_z)

  fig, ax = plt.subplots(figsize=(6.5, 6), dpi=300)
  fig.patch.set_facecolor("#FFFFFF")
  ax.set_facecolor("#FCFDFE")

  circle = plt.Circle(
      (0, 0),
      sphere_radius,
      color="#1E293B",
      fill=False,
      linestyle="--",
      linewidth=1.5,
  )
  ax.add_patch(circle)
  ax.axhline(0, color="#94A3B8", linestyle=":", linewidth=1.0)
  ax.axvline(0, color="#94A3B8", linestyle=":", linewidth=1.0)

  r_grid = np.sqrt(grid_x**2 + grid_y**2)
  z_masked = np.ma.masked_where(r_grid > sphere_radius, z_elevation)

  levels = np.linspace(0, sphere_radius, 15)
  cs = ax.contourf(
      grid_x, grid_y, z_masked, levels=levels, cmap="coolwarm", extend="both"
  )
  cbar = fig.colorbar(cs, ax=ax, shrink=0.85, pad=0.04)
  cbar.set_label(
      "Topographic Elevation (Å)",
      fontsize=11,
      fontweight="bold",
      color="#1E293B",
  )
  ax.contour(
      grid_x,
      grid_y,
      z_masked,
      levels=levels[::2],
      colors="#334155",
      linewidths=0.5,
      alpha=0.6,
  )

  ax.plot(
      0,
      0,
      marker="o",
      markersize=8,
      color="#0F172A",
      markeredgecolor="#FFFFFF",
      markeredgewidth=1.5,
      label="Metal Center",
  )
  ax.text(
      sphere_radius * 0.7,
      sphere_radius * 0.7,
      "Q1",
      fontsize=11,
      fontweight="bold",
      color="#475569",
      ha="center",
  )
  ax.text(
      -sphere_radius * 0.7,
      sphere_radius * 0.7,
      "Q2",
      fontsize=11,
      fontweight="bold",
      color="#475569",
      ha="center",
  )
  ax.text(
      -sphere_radius * 0.7,
      -sphere_radius * 0.7,
      "Q3",
      fontsize=11,
      fontweight="bold",
      color="#475569",
      ha="center",
  )
  ax.text(
      sphere_radius * 0.7,
      -sphere_radius * 0.7,
      "Q4",
      fontsize=11,
      fontweight="bold",
      color="#475569",
      ha="center",
  )

  ax.set_xlim(-sphere_radius * 1.1, sphere_radius * 1.1)
  ax.set_ylim(-sphere_radius * 1.1, sphere_radius * 1.1)
  ax.set_aspect("equal")
  ax.set_xlabel("X (Å)", fontsize=11, fontweight="bold")
  ax.set_ylabel("Y (Å)", fontsize=11, fontweight="bold")
  ax.set_title(title, fontsize=13, fontweight="bold", pad=12, color="#0F172A")
  ax.legend(loc="upper right", framealpha=0.9)

  plt.tight_layout()
  os.makedirs(os.path.dirname(os.path.abspath(output_png)), exist_ok=True)
  plt.savefig(output_png, dpi=300)
  plt.close()
  return output_png


def analyze_ligand_sterics(
    smiles_or_sdf: str,
    donor_idx: Optional[int] = None,
    donor_idx_2: Optional[int] = None,
    metal_idx: Optional[int] = None,
    metal_dist: float = 2.28,
    sphere_radius: float = 3.5,
    render_map_path: Optional[str] = None,
) -> Dict[str, Any]:
  """End-to-end ligand and coordination complex featurization pipeline."""
  if smiles_or_sdf.endswith(".sdf") or smiles_or_sdf.endswith(".mol"):
    suppl = Chem.SDMolSupplier(smiles_or_sdf, removeHs=False)
    mol = suppl[0]
  else:
    mol = Chem.MolFromSmiles(smiles_or_sdf)
    if mol is None:
      mol = Chem.MolFromSmiles(smiles_or_sdf, sanitize=False)
      if mol is not None:
        mol.UpdatePropertyCache(strict=False)
    if mol is None:
      raise ValueError(f"Could not parse SMILES '{smiles_or_sdf}'")

    mol = _normalize_dative_bonds(mol)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 0
    embed_status = AllChem.EmbedMolecule(mol, params)
    if embed_status != 0:
      params.useRandomCoords = True
      AllChem.EmbedMolecule(mol, params)

    # Force field optimization with coordination core restraints if metal present
    try:
      uff = AllChem.UFFGetMoleculeForceField(mol)
      if uff is not None and metal_idx is not None:
        matom = mol.GetAtomWithIdx(metal_idx)
        m_nbrs = [nbr.GetIdx() for nbr in matom.GetNeighbors()]
        for n_idx in m_nbrs:
          nbr = mol.GetAtomWithIdx(n_idx)
          anum = nbr.GetAtomicNum()
          if anum in (9, 17, 35, 53):  # halides
            uff.AddDistanceConstraint(metal_idx, n_idx, 2.35, 2.50, 1000.0)
          elif anum == 15:  # P
            uff.AddDistanceConstraint(metal_idx, n_idx, 2.25, 2.35, 1000.0)
          elif anum == 7:  # N
            uff.AddDistanceConstraint(metal_idx, n_idx, 2.05, 2.15, 1000.0)
          elif anum == 6:  # C
            uff.AddDistanceConstraint(metal_idx, n_idx, 2.00, 2.20, 1000.0)
          if len(nbr.GetNeighbors()) == 4:
            for sub in nbr.GetNeighbors():
              s_idx = sub.GetIdx()
              if s_idx != metal_idx:
                uff.UFFAddAngleConstraint(
                    metal_idx, n_idx, s_idx, False, 112.0, 118.0, 500.0
                )
        if len(m_nbrs) == 3:
          for i in range(len(m_nbrs)):
            for j in range(i + 1, len(m_nbrs)):
              uff.UFFAddAngleConstraint(
                  m_nbrs[i], metal_idx, m_nbrs[j], False, 86.0, 96.0, 200.0
              )
        elif (
            donor_idx is not None
            and donor_idx_2 is not None
            and donor_idx in m_nbrs
            and donor_idx_2 in m_nbrs
        ):
          uff.UFFAddAngleConstraint(
              donor_idx, metal_idx, donor_idx_2, False, 86.0, 96.0, 200.0
          )
        # Check for arene or cyclopentadienyl ring in disconnected fragments
        for ring in mol.GetRingInfo().AtomRings():
          if len(ring) in (5, 6):
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            if all(
                a.GetIsAromatic() and a.GetAtomicNum() == 6 for a in ring_atoms
            ):
              for r_idx in ring:
                uff.AddDistanceConstraint(metal_idx, r_idx, 2.15, 2.25, 1000.0)
        uff.Minimize(maxIts=1000)
      elif uff is not None:
        uff.Minimize(maxIts=500)
    except Exception:
      try:
        AllChem.MMFFOptimizeMolecule(mol)
      except Exception:
        pass

  if donor_idx is None:
    donor_idx = find_donor_atom(mol, metal_idx=metal_idx)

  donor_atom = mol.GetAtomWithIdx(donor_idx)
  donor_symbol = donor_atom.GetSymbol()

  # Calculate descriptors
  vbur = compute_buried_volume(
      mol,
      donor_idx=donor_idx,
      metal_idx=metal_idx,
      metal_dist=metal_dist,
      sphere_radius=sphere_radius,
      donor_idx_2=donor_idx_2,
  )
  cone = compute_tolman_cone_angle(
      mol,
      donor_idx=donor_idx,
      metal_idx=metal_idx,
      metal_dist=metal_dist,
      donor_idx_2=donor_idx_2,
  )
  sterimol = compute_sterimol_parameters(
      mol,
      donor_idx=donor_idx,
      metal_idx=metal_idx,
      metal_dist=metal_dist,
      donor_idx_2=donor_idx_2,
  )

  results: Dict[str, Any] = {
      "smiles": Chem.MolToSmiles(Chem.RemoveHs(mol)),
      "primary_donor": {
          "idx": donor_idx,
          "symbol": donor_symbol,
      },
      "buried_volume": vbur,
      "quadrants": vbur["quadrants"],
      "percent_v_bur": vbur["percent_v_bur"],
      "cone_angle": cone,
      "sterimol": sterimol,
      "physicochemical": {
          "molecular_weight": round(Descriptors.MolWt(mol), 2),
          "tpsa": round(Descriptors.TPSA(mol), 2),
          "num_rotatable_bonds": Descriptors.NumRotatableBonds(mol),
      },
  }

  if metal_idx is not None and metal_idx < mol.GetNumAtoms():
    matom = mol.GetAtomWithIdx(metal_idx)
    results["metal_center"] = {
        "idx": metal_idx,
        "symbol": matom.GetSymbol(),
    }

  if donor_idx_2 is not None:
    chelate = compute_chelate_bite_angle(
        mol,
        donor_idx_1=donor_idx,
        donor_idx_2=donor_idx_2,
        metal_idx=metal_idx,
        metal_dist=metal_dist,
    )
    results["chelate_descriptors"] = chelate

  if render_map_path:
    map_title = f"Topographic Steric Map (%V_bur = {vbur['percent_v_bur']}%)"
    render_steric_contour_map(
        mol,
        donor_idx=donor_idx,
        metal_idx=metal_idx,
        sphere_radius=sphere_radius,
        output_png=render_map_path,
        title=map_title,
    )
    results["steric_map_png"] = render_map_path

  return results


def main():
  parser = argparse.ArgumentParser(
      description="Quantitative Ligand Steric & Electronic Descriptors Engine."
  )
  parser.add_argument(
      "--smiles",
      type=str,
      required=True,
      help="Ligand or complex SMILES string.",
  )
  parser.add_argument(
      "--donor_idx",
      type=int,
      default=None,
      help="Index of primary donor atom (auto-detected if omitted).",
  )
  parser.add_argument(
      "--donor_idx_2",
      type=int,
      default=None,
      help="Index of second donor atom for bidentate chelates.",
  )
  parser.add_argument(
      "--metal_idx",
      type=int,
      default=None,
      help="Index of metal center atom if complex already contains metal.",
  )
  parser.add_argument(
      "--metal_dist",
      type=float,
      default=2.28,
      help="M-Donor coordination distance in Å (default: 2.28 Å).",
  )
  parser.add_argument(
      "--sphere_radius",
      type=float,
      default=3.5,
      help="Coordination sphere radius in Å (default: 3.5 Å).",
  )
  parser.add_argument(
      "--json", action="store_true", help="Output results as JSON."
  )
  parser.add_argument(
      "--output",
      "--output-json",
      "--output_json",
      dest="output_path",
      type=str,
      default=None,
      help="Save JSON steric analysis output to a file.",
  )
  parser.add_argument(
      "--render_map",
      "--render-map",
      dest="render_map_path",
      type=str,
      default=None,
      help=(
          "Optional output path to render 2D topographic steric contour map"
          " PNG."
      ),
  )
  args = parser.parse_args()

  render_png = args.render_map_path
  if render_png is None and args.output_path:
    base_dir = os.path.dirname(os.path.abspath(args.output_path))
    base_name = os.path.splitext(os.path.basename(args.output_path))[0]
    render_png = os.path.join(base_dir, f"{base_name}_map.png")

  res = analyze_ligand_sterics(
      smiles_or_sdf=args.smiles,
      donor_idx=args.donor_idx,
      donor_idx_2=args.donor_idx_2,
      metal_idx=args.metal_idx,
      metal_dist=args.metal_dist,
      sphere_radius=args.sphere_radius,
      render_map_path=render_png,
  )

  if args.output_path:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_path)), exist_ok=True
    )
    with open(args.output_path, "w") as f:
      json.dump(res, f, indent=2)
    print(f"[+] Saved steric descriptors JSON to: {args.output_path}")

  if args.json:
    print(json.dumps(res, indent=2))
    return

  print("\n=== Quantitative Ligand & Coordination Steric Profile ===")
  print(f"SMILES:             {res['smiles']}")
  if "metal_center" in res:
    print(
        f"Metal Center:       {res['metal_center']['symbol']} (Idx"
        f" {res['metal_center']['idx']})"
    )
  print(
      f"Primary Donor Atom: {res['primary_donor']['symbol']} (Idx"
      f" {res['primary_donor']['idx']})"
  )
  print(
      f"Percent Buried Vol: {res['buried_volume']['percent_v_bur']}% (%V_bur, R"
      f" = {res['buried_volume']['sphere_radius_A']} Å)"
  )
  print(
      f"  • Quadrants (%):  Q1={res['buried_volume']['quadrants']['Q1']}%,"
      f" Q2={res['buried_volume']['quadrants']['Q2']}%,"
      f" Q3={res['buried_volume']['quadrants']['Q3']}%,"
      f" Q4={res['buried_volume']['quadrants']['Q4']}%"
  )
  print(
      f"Buried Volume:      {res['buried_volume']['buried_volume_A3']} Å³"
      f" (Free: {res['buried_volume']['free_volume_A3']} Å³, Total:"
      f" {res['buried_volume']['total_sphere_volume_A3']} Å³)"
  )
  print(
      f"Tolman Cone Angle:  {res['cone_angle']['tolman_cone_angle_deg']}°"
      f" (Exact: {res['cone_angle']['exact_cone_angle_deg']}°)"
  )
  print(
      f"Sterimol (L, B1, B5): L={res['sterimol']['sterimol_L']} Å,"
      f" B1={res['sterimol']['sterimol_B1']} Å,"
      f" B5={res['sterimol']['sterimol_B5']} Å"
  )
  if "chelate_descriptors" in res:
    print(
        f"Chelate Bite Angle: {res['chelate_descriptors']['bite_angle_deg']}°"
        f" (Bite Distance: {res['chelate_descriptors']['bite_distance_A']} Å)"
    )
  print(
      f"Molecular Weight:   {res['physicochemical']['molecular_weight']} g/mol"
  )
  if "steric_map_png" in res:
    print(f"Topographic Map:    {res['steric_map_png']}")


if __name__ == "__main__":
  main()
