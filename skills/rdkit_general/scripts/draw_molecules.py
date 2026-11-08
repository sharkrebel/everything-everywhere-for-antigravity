# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "pandas",
#     "numpy",
# ]
# ///

"""High-Quality 2D Chemical Structure & Reaction Scheme Renderer (PNG / SVG).

Supports rendering:
1. Single molecules with optional core alignment and substructure highlights.
2. Aligned molecular grids with legends, auto-pagination, and core coordinates.
3. Batch rendering directly from CSV files (e.g. top N hits or leads).
4. Chemical reaction schemes from reaction SMARTS.
"""

import argparse
import math
import os
import sys
from typing import Final, List, Optional

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import AllChem, Draw, rdDepictor
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D

# Transition Metals (d-block and f-block)
TRANSITION_METALS = (
    set(range(21, 31))
    | set(range(39, 49))
    | set(range(57, 81))
    | set(range(89, 113))
)


def sanitize_legend(text: Optional[str]) -> str:
  """Normalizes Unicode characters (e.g.

  Greek eta/alpha/beta) to ASCII equivalents

  to prevent font rendering dropouts in 2D depiction engines.
  """
  if not text:
    return ""
  text = str(text)
  replacements = {
      "\u03b7": "eta",
      "\u0397": "Eta",
      "\u03b1": "alpha",
      "\u0391": "Alpha",
      "\u03b2": "beta",
      "\u0392": "Beta",
      "\u03b3": "gamma",
      "\u0393": "Gamma",
      "\u03b4": "delta",
      "\u0394": "Delta",
      "\u03bc": "mu",
      "\u039c": "Mu",
  }
  for k, v in replacements.items():
    text = text.replace(k, v)
  return text


def prepare_organometallic_depiction(mol: Chem.Mol):
  """Detects and lays out organometallic complexes (sandwich and half-sandwich/piano-stool)

  with clean centroid coordination bonds and natural 2D ligand coordinates.
  Eliminates disconnected floating fragment artifacts from database SMILES.
  Returns: (processed_mol, set_of_centroid_atom_indices)
  """
  if mol is None:
    return mol, set()

  metal_atoms = [
      a for a in mol.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS
  ]
  if not metal_atoms:
    return mol, set()

  frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
  if len(frags) <= 1:
    return mol, set()

  metal_frags = []
  arene_frags = []
  other_frags = []

  for f in frags:
    Chem.FastFindRings(f)
    if any(a.GetAtomicNum() in TRANSITION_METALS for a in f.GetAtoms()):
      metal_frags.append(f)
    else:
      ring_info = f.GetRingInfo().AtomRings()
      has_pi_ring = any(
          len(r) in (5, 6, 7, 8)
          and all(f.GetAtomWithIdx(i).GetIsAromatic() for i in r)
          for r in ring_info
      )
      if has_pi_ring:
        arene_frags.append(f)
      else:
        other_frags.append(f)

  # Case A: Sandwich complex (1 metal, 2 pi-rings, e.g. ferrocene)
  if len(metal_frags) == 1 and len(arene_frags) == 2 and not other_frags:
    m_frag = metal_frags[0]
    m_atom = [
        a for a in m_frag.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS
    ][0]
    metal_sym = m_atom.GetSymbol()

    rw = Chem.RWMol()
    metal_idx = rw.AddAtom(Chem.Atom(metal_sym))
    c_top = rw.AddAtom(Chem.Atom(0))
    c_bot = rw.AddAtom(Chem.Atom(0))
    rw.AddBond(c_top, metal_idx, Chem.BondType.ZERO)
    rw.AddBond(c_bot, metal_idx, Chem.BondType.ZERO)

    conf = Chem.Conformer()
    conf.SetAtomPosition(metal_idx, Point3D(0.0, 0.0, 0.0))
    conf.SetAtomPosition(c_top, Point3D(0.0, 2.2, 0.0))
    conf.SetAtomPosition(c_bot, Point3D(0.0, -2.2, 0.0))

    for ring_frag, y_offset in [(arene_frags[0], 2.2), (arene_frags[1], -2.2)]:
      rdDepictor.Compute2DCoords(ring_frag)
      r_conf = ring_frag.GetConformer()
      ar_rings = [
          r for r in ring_frag.GetRingInfo().AtomRings() if len(r) in (5, 6)
      ]
      ar_ring = (
          ar_rings[0] if ar_rings else list(range(ring_frag.GetNumAtoms()))
      )
      c_orig = np.mean(
          [np.array(r_conf.GetAtomPosition(i))[:2] for i in ar_ring], axis=0
      )

      r_map = {}
      for a in ring_frag.GetAtoms():
        idx = rw.AddAtom(Chem.Atom(a.GetSymbol()))
        rw.GetAtomWithIdx(idx).SetIsAromatic(a.GetIsAromatic())
        r_map[a.GetIdx()] = idx
      for b in ring_frag.GetBonds():
        rw.AddBond(
            r_map[b.GetBeginAtomIdx()],
            r_map[b.GetEndAtomIdx()],
            b.GetBondType(),
        )

      for a in ring_frag.GetAtoms():
        orig = np.array(r_conf.GetAtomPosition(a.GetIdx()))[:2] - c_orig
        conf.SetAtomPosition(
            r_map[a.GetIdx()], Point3D(orig[0], orig[1] + y_offset, 0.0)
        )

    out_mol = rw.GetMol()
    out_mol.AddConformer(conf, assignId=True)
    return out_mol, {c_top, c_bot}

  # Case B: Half-sandwich / Piano-stool (1 metal fragment, >= 1 arene ring)
  if len(metal_frags) == 1 and len(arene_frags) >= 1:
    m_frag = metal_frags[0]
    arene_frag = arene_frags[0]

    m_atom_candidates = [
        a for a in m_frag.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS
    ]
    if not m_atom_candidates:
      return mol, set()
    m_atom = m_atom_candidates[0]
    metal_sym = m_atom.GetSymbol()
    metal_local_idx = m_atom.GetIdx()

    rdDepictor.Compute2DCoords(arene_frag)
    a_conf = arene_frag.GetConformer()
    coord_rings = [
        r for r in arene_frag.GetRingInfo().AtomRings() if len(r) in (5, 6)
    ]
    coord_ring = (
        coord_rings[0]
        if coord_rings
        else [a.GetIdx() for a in arene_frag.GetAtoms() if a.GetIsAromatic()]
    )
    c_2d = np.mean(
        [np.array(a_conf.GetAtomPosition(idx))[:2] for idx in coord_ring],
        axis=0,
    )

    target_c = np.array([0.0, 2.6])
    arene_shift = target_c - c_2d

    rw = Chem.RWMol()
    metal_idx = rw.AddAtom(Chem.Atom(metal_sym))
    centroid_idx = rw.AddAtom(Chem.Atom(0))
    rw.AddBond(centroid_idx, metal_idx, Chem.BondType.ZERO)

    arene_map = {}
    for a in arene_frag.GetAtoms():
      idx = rw.AddAtom(Chem.Atom(a.GetSymbol()))
      rw.GetAtomWithIdx(idx).SetIsAromatic(a.GetIsAromatic())
      arene_map[a.GetIdx()] = idx
    for b in arene_frag.GetBonds():
      rw.AddBond(
          arene_map[b.GetBeginAtomIdx()],
          arene_map[b.GetEndAtomIdx()],
          b.GetBondType(),
      )

    donor_neighbors = [n.GetIdx() for n in m_atom.GetNeighbors()]
    metal_rw = Chem.RWMol(m_frag)
    for n_idx in donor_neighbors:
      metal_rw.RemoveBond(metal_local_idx, n_idx)
    metal_rw.RemoveAtom(metal_local_idx)

    lig_frags = Chem.GetMolFrags(
        metal_rw.GetMol(), asMols=True, sanitizeFrags=False
    )

    n_legs = len(lig_frags)
    if n_legs == 2:
      leg_angles = [225, 315]
    elif n_legs == 3:
      leg_angles = [210, 270, 330]
    elif n_legs == 1:
      leg_angles = [270]
    else:
      leg_angles = np.linspace(205, 335, n_legs)

    conf = Chem.Conformer()
    conf.SetAtomPosition(centroid_idx, Point3D(target_c[0], target_c[1], 0.0))
    conf.SetAtomPosition(metal_idx, Point3D(0.0, 0.0, 0.0))

    for a in arene_frag.GetAtoms():
      pos = np.array(a_conf.GetAtomPosition(a.GetIdx()))[:2] + arene_shift
      conf.SetAtomPosition(arene_map[a.GetIdx()], Point3D(pos[0], pos[1], 0.0))

    for lig, leg_ang in zip(lig_frags, leg_angles):
      rdDepictor.Compute2DCoords(lig)
      l_conf = lig.GetConformer()

      donor_in_lig = 0
      for target_sym in ("P", "N", "S", "Cl", "Br", "I", "O"):
        matches = [
            a.GetIdx() for a in lig.GetAtoms() if a.GetSymbol() == target_sym
        ]
        if matches:
          donor_in_lig = matches[0]
          break

      d_len = 2.4 if lig.GetNumAtoms() > 3 else 1.9
      d_pos = np.array([
          d_len * np.cos(np.deg2rad(leg_ang)),
          d_len * np.sin(np.deg2rad(leg_ang)),
      ])

      l_map = {}
      for a in lig.GetAtoms():
        idx = rw.AddAtom(Chem.Atom(a.GetSymbol()))
        rw.GetAtomWithIdx(idx).SetIsAromatic(a.GetIsAromatic())
        l_map[a.GetIdx()] = idx
      for b in lig.GetBonds():
        rw.AddBond(
            l_map[b.GetBeginAtomIdx()],
            l_map[b.GetEndAtomIdx()],
            b.GetBondType(),
        )

      donor_sym = lig.GetAtomWithIdx(donor_in_lig).GetSymbol()
      btype = (
          Chem.BondType.DATIVE
          if donor_sym in ("P", "N", "O", "S")
          else Chem.BondType.SINGLE
      )
      rw.AddBond(l_map[donor_in_lig], metal_idx, btype)

      l_orig = np.array(l_conf.GetAtomPosition(donor_in_lig))[:2]
      if lig.GetNumAtoms() > 1:
        other_atoms = [
            a.GetIdx() for a in lig.GetAtoms() if a.GetIdx() != donor_in_lig
        ]
        lig_center = np.mean(
            [np.array(l_conf.GetAtomPosition(idx))[:2] for idx in other_atoms],
            axis=0,
        )
        v = lig_center - l_orig
        rot = np.deg2rad(leg_ang) - np.arctan2(v[1], v[0])
        c, s = np.cos(rot), np.sin(rot)
        R = np.array([[c, -s], [s, c]])
      else:
        R = np.eye(2)

      for a in lig.GetAtoms():
        orig = np.array(l_conf.GetAtomPosition(a.GetIdx()))[:2] - l_orig
        new_pos = np.dot(R, orig) + d_pos
        conf.SetAtomPosition(
            l_map[a.GetIdx()], Point3D(new_pos[0], new_pos[1], 0.0)
        )

    out_mol = rw.GetMol()
    out_mol.AddConformer(conf, assignId=True)
    return out_mol, {centroid_idx}

  return mol, set()


def is_ring_convex(conf: Chem.Conformer, ring: tuple[int, ...]) -> bool:
  """Returns True if the 2D polygon for `ring` is strictly convex."""
  pts = np.array(
      [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in ring],
      dtype=float,
  )
  crosses = []
  for k in range(len(ring)):
    v1 = pts[k] - pts[k - 1]
    v2 = pts[(k + 1) % len(ring)] - pts[k]
    crosses.append(v1[0] * v2[1] - v1[1] * v2[0])
  return min(crosses) * max(crosses) > 1e-5


def _find_alignment_pattern(
    mol: Chem.Mol, core_mol: Chem.Mol
) -> Optional[Chem.Mol]:
  """Returns `core_mol` if matched, else the >=5-atom MCS query molecule."""
  if mol.HasSubstructMatch(core_mol):
    return core_mol
  mcs = rdFMCS.FindMCS(
      [mol, core_mol],
      ringMatchesRingOnly=False,
      completeRingsOnly=False,
      bondCompare=rdFMCS.BondCompare.CompareOrderExact,
      timeout=1,
  )
  if mcs.numAtoms >= 5 and mcs.smartsString:
    return Chem.MolFromSmarts(mcs.smartsString)
  return None


def regularize_nonconvex_small_rings(mol: Chem.Mol) -> None:
  """Regularizes non-convex small rings (<=7 atoms) via RDKit submolecules.

  Extracts each non-convex, non-cage small-ring component and its 1-bond
  exocyclic attachment stubs into an isolated submolecule, computes ideal 2D
  coordinates, fits the ring component onto the existing conformer via 2D
  orthogonal Procrustes (Kabsch SVD), and propagates non-macrocycle exocyclic
  branches via the rigid SO(2) motion of their attachment bonds.

  Args:
    mol: RDKit molecule with an existing 2D conformer modified in place.
  """
  if mol.GetNumConformers() == 0 or mol.GetConformer().Is3D():
    return
  Chem.GetSymmSSSR(mol)
  conf = mol.GetConformer()
  atom_rings = mol.GetRingInfo().AtomRings()
  small_rings = [set(r) for r in atom_rings if len(r) <= 7]
  macro_rings = [set(r) for r in atom_rings if len(r) >= 8]
  macro_atoms = set().union(*macro_rings) if macro_rings else set()

  for r in atom_rings:
    if len(r) > 7 or is_ring_convex(conf, r):
      continue
    r_set = set(r)
    if any(len(r_set & sr) >= 3 for sr in small_rings if sr != r_set):
      continue
    comp = set(r_set)
    changed = True
    while changed:
      changed = False
      for sr in small_rings:
        if not sr.issubset(comp) and (comp & sr):
          comp |= sr
          changed = True

    branches: list[tuple[int, int, list[int]]] = []
    stub_atoms: set[int] = set()
    for a_idx in list(comp):
      for nb in mol.GetAtomWithIdx(a_idx).GetNeighbors():
        nb_idx = nb.GetIdx()
        if nb_idx in comp or mol.GetBondBetweenAtoms(a_idx, nb_idx).IsInRing():
          continue
        visited = set(comp)
        queue = [nb_idx]
        branch = []
        in_macro = False
        while queue:
          curr = queue.pop(0)
          if curr in visited:
            continue
          visited.add(curr)
          if curr in macro_atoms:
            in_macro = True
            break
          branch.append(curr)
          for n2 in mol.GetAtomWithIdx(curr).GetNeighbors():
            if n2.GetIdx() not in visited:
              queue.append(n2.GetIdx())
        if not in_macro and branch:
          branches.append((a_idx, nb_idx, branch))
          stub_atoms.add(nb_idx)

    sub_atoms = comp | stub_atoms
    bonds = [
        b.GetIdx()
        for b in mol.GetBonds()
        if (b.GetBeginAtomIdx() in comp and b.GetEndAtomIdx() in sub_atoms)
        or (b.GetEndAtomIdx() in comp and b.GetBeginAtomIdx() in sub_atoms)
    ]
    amap: dict[int, int] = {}
    submol = Chem.PathToSubmol(mol, bonds, atomMap=amap)
    for u_idx in stub_atoms:
      if u_idx in amap:
        submol.GetAtomWithIdx(amap[u_idx]).SetIsAromatic(False)
    Chem.FastFindRings(submol)
    rdDepictor.Compute2DCoords(submol)
    sub_conf = submol.GetConformer()

    comp_indices = [i for i in comp if i in amap]
    sub_comp_indices = [amap[i] for i in comp_indices]
    x_orig_comp = np.array(
        [
            [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
            for i in comp_indices
        ],
        dtype=float,
    )
    x_sub_comp = np.array(
        [
            [sub_conf.GetAtomPosition(i).x, sub_conf.GetAtomPosition(i).y]
            for i in sub_comp_indices
        ],
        dtype=float,
    )
    c_orig = x_orig_comp.mean(axis=0)
    c_sub = x_sub_comp.mean(axis=0)
    u_mat, _, vt_mat = np.linalg.svd(
        (x_sub_comp - c_sub).T @ (x_orig_comp - c_orig)
    )
    rot = u_mat @ vt_mat

    old_pos = {
        i: np.array(
            [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y], dtype=float
        )
        for i in range(mol.GetNumAtoms())
    }
    new_pos: dict[int, np.ndarray] = {}
    for orig_i, sub_i in amap.items():
      pt_sub = np.array(
          [
              sub_conf.GetAtomPosition(sub_i).x,
              sub_conf.GetAtomPosition(sub_i).y,
          ],
          dtype=float,
      )
      new_pos[orig_i] = (pt_sub - c_sub) @ rot + c_orig

    for orig_i in comp_indices:
      pt = new_pos[orig_i]
      conf.SetAtomPosition(orig_i, Point3D(float(pt[0]), float(pt[1]), 0.0))

    if np.linalg.det(rot) < 0:
      for b in mol.GetBonds():
        if b.GetBeginAtomIdx() in comp:
          if (
              mol.GetAtomWithIdx(b.GetBeginAtomIdx()).GetChiralTag()
              != Chem.ChiralType.CHI_UNSPECIFIED
          ):
            b.SetBondDir(Chem.BondDir.NONE)
          elif b.GetBondDir() == Chem.BondDir.BEGINWEDGE:
            b.SetBondDir(Chem.BondDir.BEGINDASH)
          elif b.GetBondDir() == Chem.BondDir.BEGINDASH:
            b.SetBondDir(Chem.BondDir.BEGINWEDGE)

    for a_idx, nb_idx, branch in branches:
      v_old = old_pos[nb_idx] - old_pos[a_idx]
      v_new = new_pos[nb_idx] - new_pos[a_idx]
      n_old = float(np.linalg.norm(v_old))
      n_new = float(np.linalg.norm(v_new))
      if n_old > 1e-6 and n_new > 1e-6:
        u_old = v_old / n_old
        u_new = v_new / n_new
        cos_t = float(np.clip(np.dot(u_old, u_new), -1.0, 1.0))
        sin_t = float(u_old[0] * u_new[1] - u_old[1] * u_new[0])
        r_branch = np.array([[cos_t, -sin_t], [sin_t, cos_t]], dtype=float)
      else:
        r_branch = np.eye(2)
      p_nb_new = new_pos[nb_idx]
      p_nb_old = old_pos[nb_idx]
      for w_idx in branch:
        pt_w = p_nb_new + (old_pos[w_idx] - p_nb_old) @ r_branch.T
        conf.SetAtomPosition(
            w_idx, Point3D(float(pt_w[0]), float(pt_w[1]), 0.0)
        )


# Standard RDKit 2D bond length is 1.50 A; geminal 1,3-neighbors sit >= 2.12 A
# apart at 90 deg and 2.60 A at 120 deg, so non-bonded pairs below 0.90 A or
# covalent bonds compressed below 1.00 A indicate 2D steric/fold congestion.
_MIN_NONBONDED_2D_DIST: Final[float] = 0.90
_MIN_BONDED_2D_DIST: Final[float] = 1.00
_COORDGEN_FALLBACK_2D_DIST: Final[float] = 0.85
_BOND_COMPRESSION_WEIGHT: Final[float] = 0.85
WEDGE_BOND_DIRS: Final[frozenset[Chem.BondDir]] = frozenset({
    Chem.BondDir.BEGINWEDGE,
    Chem.BondDir.BEGINDASH,
})
_CIS_STEREO: Final[frozenset[Chem.BondStereo]] = frozenset({
    Chem.BondStereo.STEREOZ,
    Chem.BondStereo.STEREOCIS,
})
_TRANS_STEREO: Final[frozenset[Chem.BondStereo]] = frozenset({
    Chem.BondStereo.STEREOE,
    Chem.BondStereo.STEREOTRANS,
})


def depiction_clearance(mol: Chem.Mol) -> tuple[float, float]:
  """Returns `(min_nonbonded_dist, min_bonded_dist)` in 2D conformer 0."""
  if mol.GetNumAtoms() < 2 or mol.GetNumConformers() == 0:
    return (math.inf, math.inf)
  xy = mol.GetConformer().GetPositions()[:, :2]
  dists = np.hypot(
      xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1]
  )
  np.fill_diagonal(dists, math.inf)
  bonded_dists: list[float] = []
  for bond in mol.GetBonds():
    i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    bonded_dists.append(float(dists[i, j]))
    dists[i, j] = dists[j, i] = math.inf
  return (float(np.min(dists)), min(bonded_dists, default=math.inf))


def has_inverted_double_bond_2d(mol: Chem.Mol) -> bool:
  """Returns True if any E/Z double bond has inverted or collinear 2D stereo."""
  if mol.GetNumConformers() == 0:
    return False
  xy = mol.GetConformer().GetPositions()[:, :2]
  for bond in mol.GetBonds():
    stereo = bond.GetStereo()
    if stereo not in _CIS_STEREO and stereo not in _TRANS_STEREO:
      continue
    stereo_atoms = list(bond.GetStereoAtoms())
    if len(stereo_atoms) != 2:
      continue
    u_idx, v_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    s0_idx, s1_idx = stereo_atoms[0], stereo_atoms[1]
    if mol.GetBondBetweenAtoms(u_idx, s0_idx) is None:
      s0_idx, s1_idx = s1_idx, s0_idx
    vec_b = xy[v_idx] - xy[u_idx]
    if float(np.hypot(vec_b[0], vec_b[1])) <= 1e-4:
      return True
    vec_0 = xy[s0_idx] - xy[u_idx]
    vec_1 = xy[s1_idx] - xy[v_idx]
    cross_prod = (vec_b[0] * vec_0[1] - vec_b[1] * vec_0[0]) * (
        vec_b[0] * vec_1[1] - vec_b[1] * vec_1[0]
    )
    if stereo in _CIS_STEREO and cross_prod <= 1e-4:
      return True
    if stereo in _TRANS_STEREO and cross_prod >= -1e-4:
      return True
  return False


def _compute_sampled_2d_coords(mol: Chem.Mol) -> None:
  """Computes flip-sampled 2D coords, handling multicomponent salts safely."""
  frag_atoms_list = Chem.GetMolFrags(mol)
  if len(frag_atoms_list) <= 1:
    if mol.GetNumBonds() > 0:
      rdDepictor.Compute2DCoords(
          mol,
          nFlipsPerSample=2,
          nSample=100,
          sampleSeed=42,
          permuteDeg4Nodes=True,
      )
    else:
      rdDepictor.Compute2DCoords(mol)
    return
  frag_mols = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
  conf = Chem.Conformer(mol.GetNumAtoms())
  conf.Set3D(False)
  x_cursor = 0.0
  for frag_atoms, frag_mol in zip(frag_atoms_list, frag_mols, strict=True):
    Chem.GetSymmSSSR(frag_mol)
    if frag_mol.GetNumBonds() > 0:
      rdDepictor.Compute2DCoords(
          frag_mol,
          nFlipsPerSample=2,
          nSample=100,
          sampleSeed=42,
          permuteDeg4Nodes=True,
      )
    else:
      rdDepictor.Compute2DCoords(frag_mol)
    frag_xy = frag_mol.GetConformer().GetPositions()[:, :2].copy()
    frag_xy[:, 1] -= float(np.mean(frag_xy[:, 1]))
    frag_xy[:, 0] += x_cursor - float(np.min(frag_xy[:, 0]))
    x_cursor = float(np.max(frag_xy[:, 0])) + 2.5
    for local_i, orig_i in enumerate(frag_atoms):
      conf.SetAtomPosition(
          orig_i,
          Point3D(float(frag_xy[local_i, 0]), float(frag_xy[local_i, 1]), 0.0),
      )
  mol.RemoveAllConformers()
  mol.AddConformer(conf, assignId=True)


def resolve_2d_clashes(mol: Chem.Mol, *, used_template: bool = False) -> bool:
  """Unfolds congested 2D depictions when non-bonded atoms or bonds collide.

  Args:
    mol: RDKit molecule with an existing 2D conformer modified in place.
    used_template: Whether `mol` was laid out with template constraints (adds an
      unconstrained native 2D retry before flip sampling).

  Returns:
    True if `mol`'s 2D conformer was updated to a higher-clearance layout.
  """
  min_nb, min_b = depiction_clearance(mol)
  if min_nb >= _MIN_NONBONDED_2D_DIST and min_b >= _MIN_BONDED_2D_DIST:
    return False

  small_rings = [
      ring for ring in mol.GetRingInfo().AtomRings() if len(ring) <= 7
  ]
  orig_conf = mol.GetConformer()
  orig_nonconvex = sum(
      not is_ring_convex(conf=orig_conf, ring=ring) for ring in small_rings
  )
  best_score = min(min_nb, _BOND_COMPRESSION_WEIGHT * min_b)
  best_nb, best_b = min_nb, min_b
  best_mol: Chem.Mol | None = None

  modes = (*(("default",) if used_template else ()), "sampled", "coordgen")
  prev_coordgen = rdDepictor.GetPreferCoordGen()
  try:
    with rdBase.BlockLogs():
      rdDepictor.SetPreferCoordGen(False)
      for mode in modes:
        if best_nb >= _MIN_NONBONDED_2D_DIST and best_b >= _MIN_BONDED_2D_DIST:
          break
        if mode == "coordgen" and best_score >= _COORDGEN_FALLBACK_2D_DIST:
          break
        cand = Chem.Mol(mol)
        for bond in cand.GetBonds():
          if bond.GetBondDir() in WEDGE_BOND_DIRS:
            bond.SetBondDir(Chem.BondDir.NONE)
        if mode == "default":
          rdDepictor.Compute2DCoords(cand)
        elif mode == "sampled":
          _compute_sampled_2d_coords(cand)
        else:
          try:
            rdDepictor.SetPreferCoordGen(True)
            rdDepictor.Compute2DCoords(cand)
          finally:
            rdDepictor.SetPreferCoordGen(False)
          # CoordGen emits ~1.0 A bonds; rescale to RDKit's 1.5 A standard so
          # clearance thresholds, ring templates, and core alignment match.
          rdDepictor.NormalizeDepiction(cand, canonicalize=0)
        regularize_nonconvex_small_rings(cand)
        if has_inverted_double_bond_2d(cand):
          continue
        cand_conf = cand.GetConformer()
        if (
            sum(
                not is_ring_convex(conf=cand_conf, ring=ring)
                for ring in small_rings
            )
            > orig_nonconvex
        ):
          continue
        cand_nb, cand_b = depiction_clearance(cand)
        cand_score = min(cand_nb, _BOND_COMPRESSION_WEIGHT * cand_b)
        if cand_score > best_score + 1e-4:
          best_score, best_nb, best_b, best_mol = (
              cand_score,
              cand_nb,
              cand_b,
              cand,
          )
  finally:
    rdDepictor.SetPreferCoordGen(prev_coordgen)

  if best_mol is None:
    return False
  mol.RemoveAllConformers()
  mol.AddConformer(best_mol.GetConformer(), assignId=True)
  for dst_b, src_b in zip(mol.GetBonds(), best_mol.GetBonds(), strict=True):
    dst_b.SetBondDir(src_b.GetBondDir())
  return True


def _align_rigid_to_core(
    mol: Chem.Mol, core_mol: Chem.Mol, *, prefer_largest_anchor: bool = False
) -> None:
  """Rigid-body aligns `mol` onto `core_mol` via Kabsch or ring+stub anchor."""
  ref_patt = _find_alignment_pattern(mol, core_mol)
  if ref_patt is None:
    return
  mol_matches = mol.GetSubstructMatches(ref_patt, uniquify=False)[:8]
  core_matches = core_mol.GetSubstructMatches(ref_patt, uniquify=False)[:8]
  conf_m = mol.GetConformer()
  conf_c = core_mol.GetConformer()
  if mol_matches and core_matches and len(mol_matches[0]) >= 4:
    pts_m_all = conf_m.GetPositions()[:, :2]
    pts_c_all = conf_c.GetPositions()[:, :2]
    Chem.GetSymmSSSR(ref_patt)
    blocks: list[set[int]] = []
    for ring_atoms in ref_patt.GetRingInfo().AtomRings():
      merged = set(ring_atoms)
      kept: list[set[int]] = []
      for bset in blocks:
        if merged & bset:
          merged |= bset
        else:
          kept.append(bset)
      blocks = kept + [merged]
    anchors: list[list[int]] = []
    for bset in blocks:
      nbs = set(bset)
      for q_i in bset:
        for nb in ref_patt.GetAtomWithIdx(q_i).GetNeighbors():
          nbs.add(nb.GetIdx())
      if len(nbs) >= 4:
        anchors.append(sorted(nbs))

    def _fit_candidates(pm: np.ndarray, pc: np.ndarray):
      cm, cc = pm.mean(axis=0), pc.mean(axis=0)
      u_mat, _, vt_mat = np.linalg.svd((pm - cm).T @ (pc - cc))
      vt_pos = vt_mat.copy()
      if np.linalg.det(u_mat @ vt_pos) < 0:
        vt_pos[-1, :] *= -1.0
      vt_neg = vt_pos.copy()
      vt_neg[-1, :] *= -1.0
      for r_mat, is_proper in ((u_mat @ vt_pos, True), (u_mat @ vt_neg, False)):
        diff = (pm - cm) @ r_mat - (pc - cc)
        rmsd = float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))
        yield rmsd, (pts_m_all - cm) @ r_mat + cc, is_proper

    best_full = None
    best_anchor = None
    for mol_match in mol_matches:
      pm_sub = pts_m_all[list(mol_match)]
      for core_match in core_matches:
        if len(mol_match) != len(core_match):
          continue
        for rmsd_f, placed_f, _ in _fit_candidates(
            pm_sub, pts_c_all[list(core_match)]
        ):
          if best_full is None or rmsd_f < best_full[0]:
            best_full = (rmsd_f, placed_f)
        for anc in anchors:
          pm_a = pts_m_all[[mol_match[q] for q in anc]]
          pc_a = pts_c_all[[core_match[q] for q in anc]]
          for rmsd_a, placed_a, is_proper in _fit_candidates(pm_a, pc_a):
            if rmsd_a <= 0.50:
              rank = (
                  len(anc) if prefer_largest_anchor else 0,
                  float(np.ptp(placed_a[:, 0]))
                  - 0.5 * float(np.ptp(placed_a[:, 1]))
                  - 0.5 * rmsd_a
                  + (0.05 if is_proper else 0.0),
              )
              if best_anchor is None or rank > best_anchor[0]:
                best_anchor = (rank, placed_a)

    chosen_xy = None
    if best_full is not None and best_full[0] <= 0.75:
      chosen_xy = best_full[1]
    elif best_anchor is not None:
      chosen_xy = best_anchor[1]
    if chosen_xy is not None:
      for idx_atom in range(mol.GetNumAtoms()):
        conf_m.SetAtomPosition(
            idx_atom,
            Point3D(
                float(chosen_xy[idx_atom, 0]),
                float(chosen_xy[idx_atom, 1]),
                0.0,
            ),
        )
      for bond in mol.GetBonds():
        if bond.GetBondDir() in WEDGE_BOND_DIRS:
          bond.SetBondDir(Chem.BondDir.NONE)
      # Leave wedging to MolDraw2D/PrepareMolForDrawing, which adds explicit
      # chiral Hs at ring junctions first; pre-set wedges cause it to skip
      # those centers and leave the explicit H unwedged.
      return

  params = rdDepictor.ConstrainedDepictionParams()
  params.alignOnly = True
  params.acceptFailure = True
  params.adjustMolBlockWedging = True
  rdDepictor.GenerateDepictionMatching2DStructure(
      mol, core_mol, -1, ref_patt, params
  )


def compute_clean_2d_coords(
    mol: Chem.Mol | None, core_mol: Chem.Mol | None = None
) -> None:
  """Computes 2D coords with regular small rings, clash relief, and alignment.

  Args:
    mol: RDKit molecule modified in place with clean 2D conformer coordinates;
      no-op when `None`.
    core_mol: Optional scaffold or core template molecule for 2D orientation
      alignment (populated with 2D coordinates in place if conformerless).
  """
  if mol is None:
    return
  for bond in mol.GetBonds():
    if bond.GetBondDir() in WEDGE_BOND_DIRS:
      bond.SetBondDir(Chem.BondDir.NONE)
  prev_coordgen = rdDepictor.GetPreferCoordGen()
  try:
    rdDepictor.SetPreferCoordGen(False)
    Chem.GetSymmSSSR(mol)
    atom_rings = mol.GetRingInfo().AtomRings()
    has_macrocycle = any(len(ring) >= 8 for ring in atom_rings)
    core_has_macrocycle = False

    if core_mol is not None:
      Chem.GetSymmSSSR(core_mol)
      if core_mol.GetNumConformers() == 0:
        compute_clean_2d_coords(core_mol)
      core_has_macrocycle = any(
          len(ring) >= 8 for ring in core_mol.GetRingInfo().AtomRings()
      )

    # 1. When neither molecule has a macrocycle, use template-constrained
    # depiction (exact core or MCS) so rotatable open-chain bonds share pose.
    depicted = False
    if core_mol is not None and not has_macrocycle and not core_has_macrocycle:
      try:
        ref_patt = _find_alignment_pattern(mol, core_mol)
        if ref_patt is not None:
          rdDepictor.GenerateDepictionMatching2DStructure(
              mol, core_mol, refPatt=ref_patt, acceptFailure=True
          )
          depicted = not has_inverted_double_bond_2d(mol)
      except (ValueError, RuntimeError, KeyError):
        pass

    if not depicted:
      rdDepictor.Compute2DCoords(mol)

    # 2. Regularize any non-convex small ring (<= 7 atoms) and relieve clashes.
    regularize_nonconvex_small_rings(mol)
    resolved_clash = resolve_2d_clashes(mol, used_template=depicted)

    # 3. For macrocycles or unfolded poses, rigid-body align to core_mol without
    # distorting ring or rotatable-bond geometry.
    if core_mol is not None and (
        has_macrocycle or core_has_macrocycle or resolved_clash or not depicted
    ):
      try:
        _align_rigid_to_core(
            mol=mol,
            core_mol=core_mol,
            prefer_largest_anchor=not (has_macrocycle or core_has_macrocycle),
        )
      except (ValueError, RuntimeError, KeyError, np.linalg.LinAlgError):
        pass
  finally:
    rdDepictor.SetPreferCoordGen(prev_coordgen)


def render_single_molecule(
    smiles: str,
    output_path: str,
    width: int = 450,
    height: int = 400,
    highlight_smarts: Optional[str] = None,
    align_core: Optional[str] = None,
    legend: str = "",
    transparent: bool = False,
    dark_mode: bool = False,
) -> str:
  """Renders a single molecule to high-resolution PNG or vector SVG with optional core alignment.

  By default (transparent=False, dark_mode=False), renders an opaque white
  (#FFFFFF) background so black (#000000) bonds and atom labels remain legible
  in both light and dark UIs. Pass dark_mode=True (with or without
  transparent=True) to switch RDKit's drawing palette to light strokes via
  rdMolDraw2D.SetDarkMode.
  """
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    raise ValueError(f"Invalid SMILES: {smiles}")

  mol, centroids = prepare_organometallic_depiction(mol)
  core_mol = None
  if align_core and not centroids:
    core_mol = Chem.MolFromSmiles(align_core) or Chem.MolFromSmarts(align_core)
    if core_mol is not None:
      compute_clean_2d_coords(core_mol)
  if not centroids:
    compute_clean_2d_coords(mol, core_mol=core_mol)
  elif mol.GetNumConformers() == 0:
    AllChem.Compute2DCoords(mol)

  highlight_atoms = None
  if highlight_smarts:
    patt = Chem.MolFromSmarts(highlight_smarts)
    if patt and mol.HasSubstructMatch(patt):
      highlight_atoms = list(mol.GetSubstructMatch(patt))

  clean_legend = sanitize_legend(legend)

  ext = os.path.splitext(output_path)[1].lower()
  if ext == ".svg":
    d2d = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = d2d.drawOptions()
    if dark_mode:
      rdMolDraw2D.SetDarkMode(opts)
    opts.addAtomIndices = False
    opts.bondLineWidth = 2.0
    opts.clearBackground = not transparent
    for c in centroids:
      opts.atomLabels[c] = ""
    if highlight_atoms:
      d2d.DrawMolecule(mol, highlightAtoms=highlight_atoms, legend=clean_legend)
    else:
      d2d.DrawMolecule(mol, legend=clean_legend)
    d2d.FinishDrawing()
    svg_text = d2d.GetDrawingText()
    with open(output_path, "w", encoding="utf-8") as f:
      f.write(svg_text)
  else:
    d2d = rdMolDraw2D.MolDraw2DCairo(width, height)
    opts = d2d.drawOptions()
    if dark_mode:
      rdMolDraw2D.SetDarkMode(opts)
    opts.addAtomIndices = False
    opts.bondLineWidth = 2.0
    opts.clearBackground = not transparent
    for c in centroids:
      opts.atomLabels[c] = ""
    if highlight_atoms:
      d2d.DrawMolecule(mol, highlightAtoms=highlight_atoms, legend=clean_legend)
    else:
      d2d.DrawMolecule(mol, legend=clean_legend)
    d2d.FinishDrawing()
    with open(output_path, "wb") as f:
      f.write(d2d.GetDrawingText())

  return os.path.abspath(output_path)


def render_grid_molecules(
    smiles_list: List[str],
    output_path: str,
    legends: Optional[List[str]] = None,
    mols_per_row: int = 4,
    sub_img_size: tuple = (300, 300),
    highlight_smarts: Optional[str] = None,
    align_core: Optional[str] = None,
    transparent: bool = False,
    dark_mode: bool = False,
) -> str:
  """Renders multiple molecules into an aligned 2D image grid with shared core coordinate alignment."""
  mols = []
  valid_legends = []
  all_centroids = []

  core_mol = None
  if align_core:
    core_mol = Chem.MolFromSmiles(align_core) or Chem.MolFromSmarts(align_core)
    if core_mol is not None:
      compute_clean_2d_coords(core_mol)

  for idx, smi in enumerate(smiles_list):
    m = Chem.MolFromSmiles(str(smi))
    if m is not None:
      m, centroids = prepare_organometallic_depiction(m)
      if not centroids:
        compute_clean_2d_coords(m, core_mol=core_mol)
      elif m.GetNumConformers() == 0:
        AllChem.Compute2DCoords(m)
      mols.append(m)
      all_centroids.append(centroids)
      if legends and idx < len(legends):
        valid_legends.append(sanitize_legend(legends[idx]))
      else:
        valid_legends.append(f"Mol_{idx+1}")

  if not mols:
    raise ValueError("No valid molecules to render.")

  highlight_atom_lists = None
  if highlight_smarts:
    patt = Chem.MolFromSmarts(highlight_smarts)
    if patt:
      highlight_atom_lists = [
          list(m.GetSubstructMatch(patt)) if m.HasSubstructMatch(patt) else []
          for m in mols
      ]

  # Configure DrawOptions to hide centroid labels across grid
  dopts = rdMolDraw2D.MolDrawOptions()
  if dark_mode:
    rdMolDraw2D.SetDarkMode(dopts)
  dopts.bondLineWidth = 2.0
  dopts.clearBackground = not transparent
  for m, cents in zip(mols, all_centroids):
    for c in cents:
      dopts.atomLabels[c] = ""

  ext = os.path.splitext(output_path)[1].lower()
  if ext == ".svg":
    svg_text = Draw.MolsToGridImage(
        mols,
        molsPerRow=mols_per_row,
        subImgSize=sub_img_size,
        legends=valid_legends,
        highlightAtomLists=highlight_atom_lists,
        drawOptions=dopts,
        useSVG=True,
    )
    with open(output_path, "w", encoding="utf-8") as f:
      f.write(svg_text)
  else:
    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=mols_per_row,
        subImgSize=sub_img_size,
        legends=valid_legends,
        highlightAtomLists=highlight_atom_lists,
        drawOptions=dopts,
        useSVG=False,
    )
    img.save(output_path)

  return os.path.abspath(output_path)


def render_reaction_scheme(
    reaction_smarts: str,
    output_path: str,
    width: int = 800,
    height: int = 300,
    transparent: bool = False,
    dark_mode: bool = False,
) -> str:
  """Renders a chemical reaction scheme to PNG or SVG."""
  rxn = AllChem.ReactionFromSmarts(reaction_smarts)
  if rxn is None:
    raise ValueError(f"Invalid Reaction SMARTS: {reaction_smarts}")

  ext = os.path.splitext(output_path)[1].lower()
  is_svg = ext == ".svg"
  d2d = (
      rdMolDraw2D.MolDraw2DSVG(width, height)
      if is_svg
      else rdMolDraw2D.MolDraw2DCairo(width, height)
  )
  opts = d2d.drawOptions()
  if dark_mode:
    rdMolDraw2D.SetDarkMode(opts)
  opts.clearBackground = not transparent
  d2d.DrawReaction(rxn)
  d2d.FinishDrawing()
  if is_svg:
    with open(output_path, "w", encoding="utf-8") as f:
      f.write(d2d.GetDrawingText())
  else:
    with open(output_path, "wb") as f:
      f.write(d2d.GetDrawingText())

  return os.path.abspath(output_path)


def _find_matching_row(
    df: object,
    query_lower: str,
    search_cols: list[str],
    matched_indices: list[object],
) -> object:
  """Returns the first unmatched row index matching `query_lower`, or None."""
  for idx_r, row_r in df.iterrows():
    if idx_r in matched_indices:
      continue
    for sc in search_cols:
      cell_s = str(row_r.get(sc, "")).strip()
      if not cell_s or cell_s.lower() in ("nan", "none"):
        continue
      tokens = [t.strip().lower() for t in cell_s.split(",")]
      bare_tokens = [t.split(":", 1)[1].strip() for t in tokens if ":" in t]
      if (
          cell_s.lower() == query_lower
          or query_lower in tokens
          or query_lower in bare_tokens
          or f"cmpd {query_lower}" == cell_s.lower()
      ):
        return idx_r
  return None


def main():
  if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
  parser = argparse.ArgumentParser(
      description=(
          "Render high-quality 2D chemical structure depictions and reaction"
          " schemes."
      )
  )
  parser.add_argument(
      "--smiles", nargs="*", default=None, help="One or more SMILES strings."
  )
  parser.add_argument(
      "--input_csv",
      "--input-csv",
      "--input",
      dest="input_csv",
      type=str,
      default=None,
      help="Path to input CSV containing molecules.",
  )
  parser.add_argument(
      "--smiles_col",
      "--smiles-col",
      dest="smiles_col",
      type=str,
      default="SMILES",
      help="Column name for SMILES in CSV.",
  )
  parser.add_argument(
      "--name_col",
      "--name-col",
      "--id_col",
      "--id-col",
      dest="name_col",
      type=str,
      default=None,
      help="Column name for molecule labels/names in CSV.",
  )
  parser.add_argument(
      "--ids",
      "--compound_ids",
      "--compound-ids",
      dest="ids",
      nargs="*",
      default=None,
      help=(
          "Optional comma- or space-separated compound IDs, names,"
          " paper_compound_id, or patent_example_id values to select from"
          " --input_csv in order."
      ),
  )
  parser.add_argument(
      "--legend_cols",
      "--legend-cols",
      dest="legend_cols",
      nargs="*",
      default=None,
      help=(
          "Optional comma- or space-separated column names from --input_csv"
          " to format into multi-field legends (e.g."
          " 'compound_name,pchembl_value,CNS_MPO')."
      ),
  )
  parser.add_argument(
      "--top_n",
      "--top-n",
      dest="top_n",
      type=int,
      default=12,
      help="Maximum number of molecules to render from CSV (default: 12).",
  )
  parser.add_argument(
      "--reaction_smarts",
      "--reaction-smarts",
      dest="reaction_smarts",
      type=str,
      default=None,
      help="Reaction SMARTS string to render as a reaction scheme.",
  )
  parser.add_argument(
      "--output",
      "--output_path",
      "--output-path",
      "--out",
      dest="output",
      type=str,
      required=True,
      help="Output image file path (.png or .svg).",
  )
  parser.add_argument(
      "--legends", nargs="*", default=None, help="Labels/legends for molecules."
  )
  parser.add_argument(
      "--highlight_smarts",
      "--highlight-smarts",
      dest="highlight_smarts",
      type=str,
      default=None,
      help="SMARTS pattern to highlight on the structures.",
  )
  parser.add_argument(
      "--align_core",
      "--align-core",
      dest="align_core",
      type=str,
      default=None,
      help=(
          "SMARTS/SMILES of core scaffold to lock and align 2D coordinates"
          " across all structures."
      ),
  )
  parser.add_argument(
      "--mols_per_row",
      "--mols-per-row",
      dest="mols_per_row",
      type=int,
      default=4,
      help="Number of molecules per row in grid mode.",
  )
  parser.add_argument(
      "--width",
      type=int,
      default=500,
      help="Per-molecule panel width in pixels (default: 500).",
  )
  parser.add_argument(
      "--height",
      type=int,
      default=420,
      help="Per-molecule panel height in pixels (default: 420).",
  )
  parser.add_argument(
      "--transparent",
      action="store_true",
      help=(
          "Render with a transparent background instead of the default opaque"
          " white (#FFFFFF). Keep default opaque white for dark-mode reports"
          " so black (#000000) bonds remain legible, or pair with --dark-mode."
      ),
  )
  parser.add_argument(
      "--dark-mode",
      "--dark_mode",
      dest="dark_mode",
      action="store_true",
      help=(
          "Enable RDKit's dark-mode drawing palette (rdMolDraw2D.SetDarkMode)"
          " so carbon bonds and atom labels render in light colors suitable"
          " for dark backgrounds."
      ),
  )
  args = parser.parse_args()

  if args.transparent and not args.dark_mode:
    print(
        "Warning: --transparent without --dark-mode renders black (#000000)"
        " bonds on a transparent canvas. Ensure the target HTML container uses"
        " a white surface (.chem-svg-well / background: #ffffff) or pass"
        " --dark-mode for dark backgrounds.",
        file=sys.stderr,
    )

  def _grid_cell_size(
      n_mols: int, cols: int, w: int, h: int
  ) -> tuple[int, int]:
    cols = max(1, cols)
    rows = max(1, (n_mols + cols - 1) // cols)
    # Treat w <= 800 as per-cell width; if caller passed a full-canvas width (>800), divide by cols/rows with a 400x340 floor
    cell_w = w if w <= 800 else max(400, w // cols)
    cell_h = h if h <= 800 else max(340, h // rows)
    return (max(360, cell_w), max(300, cell_h))

  out_dir = os.path.dirname(os.path.abspath(args.output))
  if out_dir:
    os.makedirs(out_dir, exist_ok=True)

  if args.reaction_smarts:
    out = render_reaction_scheme(
        args.reaction_smarts,
        args.output,
        width=args.width if args.width > 500 else 900,
        height=args.height if args.height != 420 else 350,
        transparent=args.transparent,
        dark_mode=args.dark_mode,
    )
    print(f"Rendered reaction scheme to: {out}", flush=True)
  elif args.input_csv:
    df = pd.read_csv(args.input_csv)
    s_col = args.smiles_col if args.smiles_col in df.columns else None
    if not s_col:
      for c in [
          "SMILES",
          "smiles",
          "Product_SMILES",
          "Reactant_SMILES",
          "canonical_smiles",
      ]:
        if c in df.columns:
          s_col = c
          break
    if not s_col:
      raise ValueError(
          f"SMILES column '{args.smiles_col}' not found in {args.input_csv}."
      )

    # Match requested compound IDs (--ids) in order across common identifier
    # columns, supporting comma-separated alias lists or 'source:id' tokens.
    if args.ids:
      req_ids = [
          p.strip()
          for item in args.ids
          for p in str(item).split(",")
          if p.strip()
      ]
      known_id_cols = {
          "compound_id",
          "compound_name",
          "pref_name",
          "paper_compound_id",
          "patent_example_id",
          "molecule_chembl_id",
          "name",
          "id",
      }
      search_cols = ([args.name_col] if args.name_col in df.columns else []) + [
          c
          for c in df.columns
          if c != args.name_col and str(c).lower() in known_id_cols
      ]
      matched_indices: list[object] = []
      for qid in req_ids:
        if (
            found_idx := _find_matching_row(
                df, qid.lower(), search_cols, matched_indices
            )
        ) is not None:
          matched_indices.append(found_idx)
      if not matched_indices:
        sys.stderr.write(
            f"Error: None of the requested --ids {req_ids} matched any row "
            f"in {args.input_csv}.\n"
        )
        sys.exit(1)
      df_sub = df.loc[matched_indices].copy()
    else:
      df_sub = df.head(args.top_n)

    smiles_list = df_sub[s_col].tolist()
    # Build per-molecule legends from explicit --legends, multi-field
    # --legend_cols ('<ID> | col=val'), or the primary identifier column.
    legends = []
    if args.legends:
      legends = list(args.legends)
    elif args.legend_cols:
      l_cols = [
          p.strip()
          for item in args.legend_cols
          for p in str(item).split(",")
          if p.strip() and p.strip() in df_sub.columns
      ]
      if l_cols:
        for _, row_l in df_sub.iterrows():
          parts_l = []
          for idx_c, lc in enumerate(l_cols):
            val_l = row_l.get(lc)
            if pd.isna(val_l) or str(val_l).strip().lower() in (
                "",
                "nan",
                "none",
            ):
              continue
            if idx_c == 0:
              parts_l.append(str(val_l).strip())
            else:
              try:
                fval = float(val_l)
                parts_l.append(f"{lc}={fval:.2f}")
              except (ValueError, TypeError):
                parts_l.append(f"{lc}={str(val_l).strip()}")
          legends.append(" | ".join(parts_l))
    if not legends:
      if args.name_col and args.name_col in df.columns:
        legends = df_sub[args.name_col].astype(str).tolist()
      else:
        for c in [
            "compound_name",
            "Name",
            "name",
            "Reactant_Name",
            "ID",
            "id",
            "compound_id",
            "Compound_ID",
            "molecule_chembl_id",
        ]:
          if c in df.columns:
            legends = df_sub[c].astype(str).tolist()
            break

    out = render_grid_molecules(
        smiles_list,
        args.output,
        legends=legends if legends else None,
        mols_per_row=args.mols_per_row,
        sub_img_size=_grid_cell_size(
            len(smiles_list), args.mols_per_row, args.width, args.height
        ),
        highlight_smarts=args.highlight_smarts,
        align_core=args.align_core,
        transparent=args.transparent,
        dark_mode=args.dark_mode,
    )
    print(
        f"Rendered grid of {len(smiles_list)} molecules from {args.input_csv}"
        f" to: {out}",
        flush=True,
    )
  elif args.smiles:
    if len(args.smiles) == 1:
      legend = args.legends[0] if args.legends else ""
      out = render_single_molecule(
          args.smiles[0],
          args.output,
          width=args.width,
          height=args.height,
          highlight_smarts=args.highlight_smarts,
          align_core=args.align_core,
          legend=legend,
          transparent=args.transparent,
          dark_mode=args.dark_mode,
      )
      print(f"Rendered molecule to: {out}", flush=True)
    else:
      out = render_grid_molecules(
          args.smiles,
          args.output,
          legends=args.legends,
          mols_per_row=args.mols_per_row,
          sub_img_size=_grid_cell_size(
              len(args.smiles), args.mols_per_row, args.width, args.height
          ),
          highlight_smarts=args.highlight_smarts,
          align_core=args.align_core,
          transparent=args.transparent,
          dark_mode=args.dark_mode,
      )
      print(
          f"Rendered grid of {len(args.smiles)} molecules to: {out}", flush=True
      )
  else:
    print(
        "Error: Please provide --input_csv, --smiles, or --reaction_smarts.",
        file=sys.stderr,
        flush=True,
    )
    sys.exit(1)


if __name__ == "__main__":
  main()
