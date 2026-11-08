# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "rdkit>=2024.3.1",
#   "numpy>=1.26.0",
#   "pandas",
#   "networkx>=3.0",
#   "pyyaml>=6.0",
#   "matplotlib>=3.8.0",
#   "pillow>=10.0.0",
# ]
# ///

import argparse
from collections.abc import Sequence
import glob
import itertools
import json
import math
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any, Dict, Final, List, Optional, Set, Tuple

import matplotlib
import networkx as nx
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D, Point3D
import yaml

matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.getcwd()
CACHE_DIR = os.path.join(BASE_DIR, "mol_cache")
TEX_DIR = os.path.join(BASE_DIR, "tex")
IMG_DIR = os.path.join(BASE_DIR, "images")

if SCRIPT_DIR not in sys.path:
  sys.path.insert(0, SCRIPT_DIR)
if BASE_DIR not in sys.path:
  sys.path.insert(0, BASE_DIR)

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

from draw_molecules import regularize_nonconvex_small_rings
from electron_graph import (
    CurvedArrow,
    D_F_BLOCK_METALS,
    MolecularGraph,
    ReactionVerifierEngine,
    VALENCE_ELECTRONS,
)
import reaction_specs_100
import reactions_data
from rotation_solver import solve_reaction_rotations


def ensure_dirs(base_dir: str | None = None, out_dir: str | None = None):
  """Ensure runtime working directories exist on-demand without top-level import side effects."""
  global BASE_DIR, CACHE_DIR, TEX_DIR, IMG_DIR
  if base_dir:
    BASE_DIR = os.path.abspath(base_dir)
  CACHE_DIR = os.path.join(BASE_DIR, "mol_cache")
  TEX_DIR = os.path.join(BASE_DIR, "tex")
  IMG_DIR = (
      os.path.abspath(out_dir) if out_dir else os.path.join(BASE_DIR, "images")
  )
  os.makedirs(CACHE_DIR, exist_ok=True)
  os.makedirs(TEX_DIR, exist_ok=True)
  os.makedirs(IMG_DIR, exist_ok=True)


def _atomic_write_bytes(dst_path: str, payload: bytes) -> None:
  """Atomically writes `payload` to `dst_path` via a temp file."""
  dst_dir = os.path.dirname(os.path.abspath(dst_path))
  base = os.path.basename(dst_path)
  os.makedirs(dst_dir, exist_ok=True)
  with tempfile.NamedTemporaryFile(
      dir=dst_dir, prefix=f".{base}.", suffix=".tmp", delete=False
  ) as tmp_f:
    tmp_path = tmp_f.name
    tmp_f.write(payload)
  try:
    os.replace(tmp_path, dst_path)
  except OSError:
    if os.path.exists(tmp_path):
      os.unlink(tmp_path)
    raise


def _resolve_mol_atom_token(mol, token, map_to_idx=None):
  """Resolves an integer atom index or ':N' SMILES atom-map token on an RDKit Mol."""
  if isinstance(token, str) and token.startswith(":"):
    want = int(token[1:])
    if map_to_idx and want in map_to_idx:
      return map_to_idx[want]
    hits = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum() == want]
    if not hits:
      raise ValueError(
          f"Atom map label '{token}' matches no atom. Write the atom as "
          f"[X:{want}] in the SMILES to tag it."
      )
    if len(hits) > 1:
      raise ValueError(
          f"Atom map label '{token}' is ambiguous: matches atoms {hits}."
      )
    return hits[0]
  return int(token)


def _parse_smiles_preserving_mapped_hs(smi: str) -> Chem.Mol | None:
  """Parses `smi` while retaining atom-mapped explicit hydrogens."""
  if "[H:" in smi or "[2H:" in smi or "[3H:" in smi:
    ps = Chem.SmilesParserParams()
    ps.removeHs = False
    mol = Chem.MolFromSmiles(smi, ps)
    if mol is None:
      return None
    rm_params = Chem.RemoveHsParameters()
    rm_params.removeDegreeZero = False
    rm_params.removeOnlyHNeighbors = False
    rm_params.removeIsotopes = False
    rm_params.removeMapped = False
    rm_params.removeDefiningBondStereo = False
    rm_params.showWarnings = False
    return Chem.RemoveHs(mol, rm_params)
  return Chem.MolFromSmiles(smi)


def latex_escape_label(text: str) -> str:
  """Escapes plain-text TikZ labels while preserving explicit `$...$` math."""
  if not text:
    return ""
  parts = re.split(r"(\$[^$]*\$)", str(text))
  out = []
  for i, part in enumerate(parts):
    if i % 2 == 1:
      out.append(part)
    else:
      s = part.replace("\\\\", "\x00NL\x00").replace("\n", "\x00NL\x00")
      s = s.replace("&", r"\&")
      s = s.replace("%", r"\%")
      s = s.replace("#", r"\#")
      s = s.replace("_", r"\_")
      s = s.replace("·", r"$\cdot$")
      s = s.replace("±", r"$\pm$")
      s = s.replace("<-->", r"$\leftrightarrow$")
      s = s.replace("<->", r"$\leftrightarrow$")
      s = s.replace("-->", r"$\rightarrow$")
      s = s.replace("->", r"$\rightarrow$")
      s = s.replace("↔", r"$\leftrightarrow$")
      s = s.replace("→", r"$\rightarrow$")
      s = s.replace("\x00NL\x00", r" \\ ")
      out.append(s)
  return "".join(out)


def _reflect_pts_across_bond(
    pts: dict[int, tuple[float, float]],
    u_idx: int,
    v_idx: int,
    moving_atoms: set[int],
) -> dict[int, tuple[float, float]]:
  """Reflects `moving_atoms` across bond `u_idx-v_idx` (180-deg 3D rotamer)."""
  ux, uy = pts[u_idx]
  vx, vy = pts[v_idx]
  dx, dy = vx - ux, vy - uy
  bond_len = math.hypot(dx, dy)
  if bond_len < 1e-8:
    return dict(pts)
  tx, ty = dx / bond_len, dy / bond_len
  out = dict(pts)
  for idx in moving_atoms:
    px, py = pts[idx][0] - ux, pts[idx][1] - uy
    proj = px * tx + py * ty
    perp_x = px - proj * tx
    perp_y = py - proj * ty
    out[idx] = (ux + proj * tx - perp_x, uy + proj * ty - perp_y)
  return out


def _collect_bond_subtree(
    mol: Chem.Mol, u_block: int, v_start: int, n_atoms_limit: int
) -> set[int]:
  """Returns connected atoms on the `v_start` side of `u_block-v_start`."""
  comp: set[int] = set()
  stack = [v_start]
  while stack:
    curr = stack.pop()
    if curr in comp or curr == u_block or curr >= n_atoms_limit:
      continue
    comp.add(curr)
    for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
      stack.append(nbr.GetIdx())
  return comp


_ATROPISOMER_STEREOS: Final = tuple(
    st
    for st in (
        getattr(Chem.BondStereo, "STEREOATROPCW", None),
        getattr(Chem.BondStereo, "STEREOATROPCCW", None),
    )
    if st is not None
)


def _clear_stale_single_bond_wedges(mol: Chem.Mol) -> None:
  """Clears cached wedge/dash flags so RDKit recomputes CIP-faithful wedges."""
  for b in mol.GetBonds():
    if b.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH):
      b.SetBondDir(Chem.BondDir.NONE)


def _is_planar_achiral_mol(mol: Chem.Mol) -> bool:
  """Returns True if `mol` lacks tetrahedral, wedge, or axial stereocenters."""
  if any(
      a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
      for a in mol.GetAtoms()
  ):
    return False
  for b in mol.GetBonds():
    if b.GetBondDir() in (
        Chem.BondDir.BEGINWEDGE,
        Chem.BondDir.BEGINDASH,
        Chem.BondDir.UNKNOWN,
    ):
      return False
    if b.GetStereo() in _ATROPISOMER_STEREOS:
      return False
  return True


def _get_conf_pts_2d(
    mol: Chem.Mol, n_atoms: int | None = None
) -> dict[int, tuple[float, float]]:
  """Extracts 2D `(x, y)` coordinates from `mol`'s conformer."""
  conf = mol.GetConformer()
  limit = mol.GetNumAtoms() if n_atoms is None else n_atoms
  return {
      i: (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)
      for i in range(limit)
  }


def _set_conf_pts_2d(
    mol: Chem.Mol, pts: dict[int, tuple[float, float]]
) -> None:
  """Writes 2D `(x, y)` coordinates into `mol` and clears cached wedges."""
  conf = mol.GetConformer()
  for k_idx, (px_val, py_val) in pts.items():
    conf.SetAtomPosition(k_idx, Point3D(float(px_val), float(py_val), 0.0))
  _clear_stale_single_bond_wedges(mol)


def _resolve_intramol_heavy_pairs(
    mol: Chem.Mol,
    intramol_targets: list[tuple[Any, ...]] | None,
    map_to_idx: dict[int, int],
    mol_name: str,
    h_to_parent: dict[int, int] | None = None,
    n_heavy_limit: int | None = None,
) -> tuple[list[tuple[int, int]], dict[int, int]]:
  """Resolves intramolecular `(src, tgt[, conn])` tokens to heavy-atom pairs.

  Args:
      mol: Target RDKit molecule.
      intramol_targets: Intramolecular `(src, tgt[, conn])` token tuples.
      map_to_idx: Mapping from SMILES atom-map integer to 0-based atom index.
      mol_name: Molecule identifier for diagnostic warnings.
      h_to_parent: Optional mapping from added explicit H index to parent atom.
      n_heavy_limit: Optional upper bound on original heavy-atom indices.

  Returns:
      Tuple `(heavy_pairs, h_attackers)` of deduplicated `(src_heavy,
      tgt_heavy)` index pairs and explicit-H attacker indices.
  """
  h_parent_map = h_to_parent or {}
  limit = n_heavy_limit if n_heavy_limit is not None else mol.GetNumAtoms()
  h_attackers: dict[int, int] = {}
  heavy_pairs: list[tuple[int, int]] = []

  def _endpoint_indices(tok: Any) -> list[int]:
    if isinstance(tok, (list, tuple)):
      return [_resolve_mol_atom_token(mol, sub_t, map_to_idx) for sub_t in tok]
    return [_resolve_mol_atom_token(mol, tok, map_to_idx)]

  def _to_heavy(idx: int) -> int | None:
    if idx in h_parent_map:
      return h_parent_map[idx]
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() == 1:
      nbrs = [
          nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1
      ]
      return nbrs[0] if nbrs else None
    return idx if idx < limit else None

  for item in intramol_targets or []:
    src_tok, tgt_tok = item[0], item[1]
    conn_tok = item[2] if len(item) >= 3 else None
    try:
      eff_src_tok = conn_tok if conn_tok is not None else src_tok
      s_raw = _endpoint_indices(eff_src_tok)
      t_raw = _endpoint_indices(tgt_tok)
      if len(s_raw) == 1 and len(t_raw) == 1 and t_raw[0] in h_parent_map:
        if s_raw[0] < limit:
          h_attackers[t_raw[0]] = s_raw[0]
      s_cands = [h for i in s_raw if (h := _to_heavy(i)) is not None]
      t_cands = [h for i in t_raw if (h := _to_heavy(i)) is not None]
      best_pair = None
      best_plen = 10**9
      for si in s_cands:
        for ti in t_cands:
          if si == ti:
            continue
          plen = len(Chem.rdmolops.GetShortestPath(mol, si, ti))
          if 0 < plen < best_plen:
            best_plen = plen
            best_pair = (si, ti)
      if best_pair is not None and best_pair not in heavy_pairs:
        heavy_pairs.append(best_pair)
    except (KeyError, ValueError, IndexError, TypeError) as exc:
      print(
          f"Warning: Could not resolve intramolecular target ({src_tok},"
          f" {tgt_tok}) on {mol_name}: {exc}",
          file=sys.stderr,
      )
  return heavy_pairs, h_attackers


def _fold_macrocycle_rotamers(
    mol: Chem.Mol,
    pts: dict[int, tuple[float, float]],
    s_heavy: int,
    t_heavy: int,
    path: list[int],
    n_atoms_limit: int,
) -> dict[int, tuple[float, float]] | None:
  """Folds acyclic single bonds into a collision-free C-shaped reactive bay.

  Args:
      mol: Target RDKit molecule.
      pts: Mapping from atom index to 2D `(x, y)` coordinates.
      s_heavy: Nucleophilic heavy-atom index.
      t_heavy: Electrophilic heavy-atom index.
      path: Shortest atom-index path connecting `s_heavy` and `t_heavy`.
      n_atoms_limit: Upper bound on atom indices included in subtree flips.

  Returns:
      Updated 2D coordinate mapping when a stereoretentive 180-deg rotamer
      improves pocket geometry without steric clashes, or `None` if unchanged.
  """
  rot_bonds: list[tuple[int, int, set[int]]] = []
  for pi in range(len(path) - 1):
    u_b, v_b = path[pi], path[pi + 1]
    b_obj = mol.GetBondBetweenAtoms(u_b, v_b)
    if (
        b_obj is not None
        and not b_obj.IsInRing()
        and b_obj.GetBondType() == Chem.BondType.SINGLE
        and b_obj.GetStereo() not in _ATROPISOMER_STEREOS
        and mol.GetAtomWithIdx(u_b).GetDegree() >= 2
        and mol.GetAtomWithIdx(v_b).GetDegree() >= 2
    ):
      # Skip linear sp centers (alkynes/nitriles/allenes) where flip is a no-op.
      end_bonds = (
          *mol.GetAtomWithIdx(u_b).GetBonds(),
          *mol.GetAtomWithIdx(v_b).GetBonds(),
      )
      if any(nb_b.GetBondTypeAsDouble() >= 2.5 for nb_b in end_bonds):
        continue
      rot_bonds.append(
          (u_b, v_b, _collect_bond_subtree(mol, u_b, v_b, n_atoms_limit))
      )
  if not rot_bonds:
    return None

  active_atoms = [i for i in range(n_atoms_limit) if i in pts]
  is_heavy = {i: mol.GetAtomWithIdx(i).GetAtomicNum() > 1 for i in active_atoms}
  bonds = [
      (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
      for b in mol.GetBonds()
      if b.GetBeginAtomIdx() in pts and b.GetEndAtomIdx() in pts
  ]

  def _src_lone_pair_dot(
      p_map: dict[int, tuple[float, float]], d_st_val: float
  ) -> float:
    if d_st_val < 1e-6:
      return 0.0
    ps, pt = p_map[s_heavy], p_map[t_heavy]
    s_nbs = [
        nb.GetIdx()
        for nb in mol.GetAtomWithIdx(s_heavy).GetNeighbors()
        if nb.GetIdx() in p_map
    ]
    vx = -sum(p_map[nb][0] - ps[0] for nb in s_nbs)
    vy = -sum(p_map[nb][1] - ps[1] for nb in s_nbs)
    ns = math.hypot(vx, vy)
    if ns < 0.25:
      return 0.0
    proj = (vx / ns) * (pt[0] - ps[0]) + (vy / ns) * (pt[1] - ps[1])
    return proj / d_st_val

  ps0, pt0 = pts[s_heavy], pts[t_heavy]
  d_st0 = math.hypot(pt0[0] - ps0[0], pt0[1] - ps0[1])
  has_chord_cross = any(
      u1 not in (s_heavy, t_heavy)
      and v1 not in (s_heavy, t_heavy)
      and intersect_segments(ps0, pt0, pts[u1], pts[v1])
      for u1, v1 in bonds
  )
  dot_s0 = _src_lone_pair_dot(pts, d_st0)
  needs_fold = (
      d_st0 < 2.35
      or has_chord_cross
      or (len(path) >= 6 and (d_st0 > 4.8 or dot_s0 < -0.25))
  )
  if not needs_fold:
    return None

  base_pair_dists: dict[tuple[int, int], float] = {}
  for idx_i, i in enumerate(active_atoms):
    for j in active_atoms[idx_i + 1 :]:
      if mol.GetBondBetweenAtoms(i, j) is None:
        base_pair_dists[(i, j)] = math.hypot(
            pts[i][0] - pts[j][0], pts[i][1] - pts[j][1]
        )

  base_bond_cross: set[tuple[int, int]] = set()
  for bi, (u1, v1) in enumerate(bonds):
    for bj in range(bi + 1, len(bonds)):
      u2, v2 = bonds[bj]
      if (
          u1 not in (u2, v2)
          and v1 not in (u2, v2)
          and intersect_segments(pts[u1], pts[v1], pts[u2], pts[v2])
      ):
        base_bond_cross.add((bi, bj))

  def _score_rotamer(
      cand: dict[int, tuple[float, float]], n_flips: int
  ) -> float | None:
    min_changed_nb = 2.5
    for (i, j), d0 in base_pair_dists.items():
      d = math.hypot(cand[i][0] - cand[j][0], cand[i][1] - cand[j][1])
      if abs(d - d0) > 1e-4:
        min_steric = 1.15 if (is_heavy[i] and is_heavy[j]) else 0.85
        if d < min_steric:
          return None
        if is_heavy[i] and is_heavy[j] and d < min_changed_nb:
          min_changed_nb = d
    for bi, (u1, v1) in enumerate(bonds):
      for bj in range(bi + 1, len(bonds)):
        if (bi, bj) in base_bond_cross:
          continue
        u2, v2 = bonds[bj]
        if u1 in (u2, v2) or v1 in (u2, v2):
          continue
        if intersect_segments(cand[u1], cand[v1], cand[u2], cand[v2]):
          return None
    ps, pt = cand[s_heavy], cand[t_heavy]
    d_st = math.hypot(pt[0] - ps[0], pt[1] - ps[1])
    if d_st < 2.25:
      return None
    chord_cross = 0
    min_chord_d = 999.0
    for u1, v1 in bonds:
      if u1 in (s_heavy, t_heavy) or v1 in (s_heavy, t_heavy):
        continue
      if intersect_segments(ps, pt, cand[u1], cand[v1]):
        chord_cross += 1
      min_chord_d = min(
          min_chord_d,
          dist_pt_to_seg(cand[u1], ps, pt),
          dist_pt_to_seg(cand[v1], ps, pt),
          dist_pt_to_seg(ps, cand[u1], cand[v1]),
          dist_pt_to_seg(pt, cand[u1], cand[v1]),
      )
    dot_s = _src_lone_pair_dot(cand, d_st)
    return (
        -60.0 * chord_cross
        - 1.8 * abs(d_st - 3.2)
        + 1.5 * dot_s
        + 1.2 * min(min_chord_d, 1.5)
        + 0.8 * min(min_changed_nb, 2.0)
        - 0.25 * n_flips
    )

  def _apply_mask(mask: tuple[int, ...]) -> dict[int, tuple[float, float]]:
    cur = dict(pts)
    for bit, (u_b, v_b, mset) in zip(mask, rot_bonds, strict=True):
      if bit:
        cur = _reflect_pts_across_bond(cur, u_b, v_b, mset)
    return cur

  orig_sc = _score_rotamer(pts, 0)
  k_rot = len(rot_bonds)
  best_sc: float | None = orig_sc
  best_mask: tuple[int, ...] = (0,) * k_rot
  best_cand: dict[int, tuple[float, float]] | None = None

  if k_rot <= 8:
    candidate_masks = [
        m for m in itertools.product((0, 1), repeat=k_rot) if any(m)
    ]
  else:
    candidate_masks = []
    for r_bits in (1, 2):
      for idxs in itertools.combinations(range(k_rot), r_bits):
        m_arr = [0] * k_rot
        for idx_bit in idxs:
          m_arr[idx_bit] = 1
        candidate_masks.append(tuple(m_arr))

  for mask in candidate_masks:
    cand = _apply_mask(mask)
    sc = _score_rotamer(cand, sum(mask))
    if sc is not None and (best_sc is None or sc > best_sc):
      best_sc = sc
      best_mask = mask
      best_cand = cand

  if k_rot > 8 and best_cand is not None:
    for _ in range(2):
      improved = False
      for i in range(k_rot):
        trial = list(best_mask)
        trial[i] = 1 - trial[i]
        t_mask = tuple(trial)
        cand = _apply_mask(t_mask)
        sc = _score_rotamer(cand, sum(t_mask))
        if sc is not None and (best_sc is None or sc > best_sc + 1e-6):
          best_sc = sc
          best_mask = t_mask
          best_cand = cand
          improved = True
      if not improved:
        break

  if best_cand is not None and (
      orig_sc is None or (best_sc is not None and best_sc > orig_sc + 0.75)
  ):
    return best_cand
  return None


def _is_stereo_sensitive_pivot(mol: Chem.Mol, atom_idx: int) -> bool:
  """Returns True if rotating one bond angle at `atom_idx` alters stereochem."""
  atom = mol.GetAtomWithIdx(atom_idx)
  if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
    return True
  for b_inc in atom.GetBonds():
    if b_inc.GetStereo() != Chem.BondStereo.STEREONONE:
      return True
    if b_inc.GetBondType() == Chem.BondType.DOUBLE:
      return True
  return False


def _unfold_crowded_intramolecular_pockets(
    mol: Chem.Mol,
    base_pts: dict[int, tuple[float, float]],
    heavy_pairs: list[tuple[int, int]],
    n_heavy_orig: int,
) -> dict[int, tuple[float, float]]:
  """Folds macrocyclic linkers or unfolds crowded (<2.35 A) reactive pockets."""
  pts = dict(base_pts)
  folded_bridge_bonds: set[tuple[int, int]] = set()
  for s_heavy, t_heavy in heavy_pairs:
    if s_heavy not in pts or t_heavy not in pts:
      continue
    path = list(Chem.rdmolops.GetShortestPath(mol, s_heavy, t_heavy))
    path_edges = {
        (min(path[pi], path[pi + 1]), max(path[pi], path[pi + 1]))
        for pi in range(len(path) - 1)
    }
    half_edges = max(1, len(path_edges) // 2)
    if len(path) >= 5 and len(path_edges & folded_bridge_bonds) < half_edges:
      rot_pts = _fold_macrocycle_rotamers(
          mol, pts, s_heavy, t_heavy, path, n_heavy_orig
      )
      if rot_pts is not None:
        pts = rot_pts
        folded_bridge_bonds.update(path_edges)
        continue

    dx0 = pts[s_heavy][0] - pts[t_heavy][0]
    dy0 = pts[s_heavy][1] - pts[t_heavy][1]
    if math.hypot(dx0, dy0) >= 2.35:
      continue
    rot_bonds = []
    for pi in range(len(path) - 1):
      u_b, v_b = path[pi], path[pi + 1]
      b_obj = mol.GetBondBetweenAtoms(u_b, v_b)
      if (
          b_obj is not None
          and not b_obj.IsInRing()
          and b_obj.GetBondType() == Chem.BondType.SINGLE
          and b_obj.GetStereo() not in _ATROPISOMER_STEREOS
      ):
        if not _is_stereo_sensitive_pivot(mol, u_b):
          rot_bonds.append((u_b, v_b))
        elif not _is_stereo_sensitive_pivot(mol, v_b):
          rot_bonds.append((v_b, u_b))
    if not rot_bonds:
      continue

    u_piv, v_start = rot_bonds[len(rot_bonds) // 2]
    comp_v = _collect_bond_subtree(mol, u_piv, v_start, n_heavy_orig)
    comp_rest = {i for i in range(n_heavy_orig) if i in pts} - comp_v - {u_piv}
    ux, uy = pts[u_piv]
    best_pts = dict(pts)
    best_rot_score = -1e9
    for r_deg in range(-160, 170, 10):
      r_rad = math.radians(r_deg)
      cr, sr = math.cos(r_rad), math.sin(r_rad)
      cand_pts = dict(pts)
      for c_idx in comp_v:
        ox, oy = pts[c_idx][0] - ux, pts[c_idx][1] - uy
        cand_pts[c_idx] = (
            ux + ox * cr - oy * sr,
            uy + ox * sr + oy * cr,
        )
      d_st = math.hypot(
          cand_pts[s_heavy][0] - cand_pts[t_heavy][0],
          cand_pts[s_heavy][1] - cand_pts[t_heavy][1],
      )
      if d_st < 2.40 or d_st > 3.35:
        continue
      min_nb = min(
          (
              math.hypot(
                  cand_pts[a][0] - cand_pts[b][0],
                  cand_pts[a][1] - cand_pts[b][1],
              )
              for a in comp_v
              for b in comp_rest
          ),
          default=2.0,
      )
      if min_nb < 1.25:
        continue
      r_score = min_nb - 1.2 * abs(d_st - 2.75) - 0.002 * abs(r_deg)
      if r_score > best_rot_score:
        best_rot_score = r_score
        best_pts = cand_pts
    pts = best_pts
  return pts


def _try_swap_exterior_leaf_substituent(
    mol: Chem.Mol,
    coord_map: Dict[int, Tuple[float, float]],
    parent_idx: int,
    h_idx: int,
) -> bool:
  """Swaps a terminal leaf substituent to the interior bisector on a degree-3 center.

  When `parent_idx` has 2 chain neighbors and 1 terminal heteroatom/leaf
  substituent occupying the wide exterior bisector, moving that leaf substituent
  to the interior bisector frees the exterior bisector for the reactive explicit
  hydrogen `h_idx` so curved arrows have an unobstructed approach path.
  """
  px, py = coord_map[parent_idx]
  heavy_nbrs = [
      nbr.GetIdx()
      for nbr in mol.GetAtomWithIdx(parent_idx).GetNeighbors()
      if nbr.GetIdx() in coord_map and nbr.GetAtomicNum() > 1
  ]
  if len(heavy_nbrs) != 3:
    return False

  leaf_subs = [
      nidx
      for nidx in heavy_nbrs
      if sum(
          1
          for nb2 in mol.GetAtomWithIdx(nidx).GetNeighbors()
          if nb2.GetAtomicNum() > 1
      )
      == 1
  ]
  chain_nbrs = [nidx for nidx in heavy_nbrs if nidx not in leaf_subs]
  if len(leaf_subs) != 1 or len(chain_nbrs) != 2:
    return False

  sub_idx = leaf_subs[0]
  c1, c2 = chain_nbrs
  a_c1 = math.atan2(coord_map[c1][1] - py, coord_map[c1][0] - px)
  a_c2 = math.atan2(coord_map[c2][1] - py, coord_map[c2][0] - px)
  a_sub = math.atan2(coord_map[sub_idx][1] - py, coord_map[sub_idx][0] - px)

  # Interior bisector of (c1, c2)
  ux = math.cos(a_c1) + math.cos(a_c2)
  uy = math.sin(a_c1) + math.sin(a_c2)
  if math.hypot(ux, uy) <= 1e-3:
    return False

  ang_in = math.atan2(uy, ux)
  diff_sub_in = abs((a_sub - ang_in + math.pi) % (2 * math.pi) - math.pi)
  if diff_sub_in <= math.radians(120.0):
    return False

  r_sub = math.hypot(coord_map[sub_idx][0] - px, coord_map[sub_idx][1] - py)
  new_sub_x = px + r_sub * math.cos(ang_in)
  new_sub_y = py + r_sub * math.sin(ang_in)
  old_sub_x, old_sub_y = coord_map[sub_idx]
  coord_map[sub_idx] = (new_sub_x, new_sub_y)
  # Shift any already-placed H on sub_idx by the same delta
  for sub_nb in mol.GetAtomWithIdx(sub_idx).GetNeighbors():
    s_nb_idx = sub_nb.GetIdx()
    if s_nb_idx != parent_idx and s_nb_idx in coord_map:
      ox, oy = coord_map[s_nb_idx]
      coord_map[s_nb_idx] = (
          new_sub_x + (ox - old_sub_x),
          new_sub_y + (oy - old_sub_y),
      )
  coord_map[h_idx] = (
      px + 1.12 * math.cos(a_sub),
      py + 1.12 * math.sin(a_sub),
  )
  return True


def _place_added_explicit_hydrogens(
    mol: Chem.Mol,
    base_pts: Dict[int, Tuple[float, float]],
    added_h_pairs: List[Tuple[int, int]],
    h_attackers: Dict[int, int],
) -> Dict[int, Tuple[float, float]]:
  """Computes 2D coordinates for added explicit hydrogens to avoid bond overlaps."""
  coord_map = dict(base_pts)
  for parent_idx, h_idx in added_h_pairs:
    if _try_swap_exterior_leaf_substituent(mol, coord_map, parent_idx, h_idx):
      continue

    px, py = coord_map[parent_idx]
    nbr_angles = [
        math.atan2(
            coord_map[nbr.GetIdx()][1] - py, coord_map[nbr.GetIdx()][0] - px
        )
        for nbr in mol.GetAtomWithIdx(parent_idx).GetNeighbors()
        if nbr.GetIdx() in coord_map
    ]
    atk_idx = h_attackers.get(h_idx)
    target_ang = None
    if atk_idx is not None and atk_idx in coord_map:
      target_ang = math.atan2(
          coord_map[atk_idx][1] - py, coord_map[atk_idx][0] - px
      )

    # Scan candidate angles in 5-degree increments to maximize angular clearance
    # from existing bonds and align toward any intramolecular attacking atom.
    best_ang = 0.0
    best_score = -1e9
    for deg in range(0, 360, 5):
      rad = math.radians(deg)
      min_ang_diff = (
          min(
              abs((rad - na + math.pi) % (2 * math.pi) - math.pi)
              for na in nbr_angles
          )
          if nbr_angles
          else math.pi
      )
      if min_ang_diff < math.radians(30.0):
        continue
      hx = px + 1.12 * math.cos(rad)
      hy = py + 1.12 * math.sin(rad)
      min_dist = min(
          (
              math.hypot(hx - ox, hy - oy)
              for oi, (ox, oy) in coord_map.items()
              if oi != parent_idx
          ),
          default=2.0,
      )
      score = min_dist + 0.6 * min_ang_diff
      if target_ang is not None:
        score += 3.5 * math.cos(rad - target_ang)
      if score > best_score:
        best_score = score
        best_ang = rad
    coord_map[h_idx] = (
        px + 1.12 * math.cos(best_ang),
        py + 1.12 * math.sin(best_ang),
    )
  return coord_map


def _align_geminal_bond_shift_substituents(
    mol: Chem.Mol,
    bond_shift_pairs: list[tuple[Any, Any]] | None,
    map_to_idx: dict[int, int],
) -> None:
  """Places acyclic bond-shift substituents adjacent to their target bond."""
  if not bond_shift_pairs or mol.GetNumConformers() == 0:
    return
  conf = mol.GetConformer()
  n_atoms = mol.GetNumAtoms()
  for b1_raw, b2_raw in bond_shift_pairs:
    if not (
        isinstance(b1_raw, (list, tuple))
        and isinstance(b2_raw, (list, tuple))
        and len(b1_raw) == 2
        and len(b2_raw) == 2
    ):
      continue
    try:
      b1 = (
          _resolve_mol_atom_token(mol, b1_raw[0], map_to_idx),
          _resolve_mol_atom_token(mol, b1_raw[1], map_to_idx),
      )
      b2 = (
          _resolve_mol_atom_token(mol, b2_raw[0], map_to_idx),
          _resolve_mol_atom_token(mol, b2_raw[1], map_to_idx),
      )
    except (TypeError, ValueError):
      continue
    if (
        mol.GetBondBetweenAtoms(*b1) is None
        or mol.GetBondBetweenAtoms(*b2) is None
    ):
      continue
    shared = set(b1) & set(b2)
    if len(shared) != 1:
      continue
    c_idx = next(iter(shared))
    if mol.GetAtomWithIdx(c_idx).GetDegree() < 4:
      continue
    a_idx = b1[0] if b1[1] == c_idx else b1[1]
    b_idx = b2[0] if b2[1] == c_idx else b2[1]
    pc = conf.GetAtomPosition(c_idx)
    pa = conf.GetAtomPosition(a_idx)
    pb = conf.GetAtomPosition(b_idx)
    ang_a = math.atan2(pa.y - pc.y, pa.x - pc.x)
    ang_b = math.atan2(pb.y - pc.y, pb.x - pc.x)
    other_nbrs = [
        nb.GetIdx()
        for nb in mol.GetAtomWithIdx(c_idx).GetNeighbors()
        if nb.GetIdx() not in (a_idx, b_idx)
    ]
    other_angs = [
        math.atan2(
            conf.GetAtomPosition(ni).y - pc.y,
            conf.GetAtomPosition(ni).x - pc.x,
        )
        for ni in other_nbrs
    ]
    d_ccw = (ang_b - ang_a) % (2.0 * math.pi)
    d_cw = (ang_a - ang_b) % (2.0 * math.pi)
    has_ccw = any(
        0.05 < ((oa - ang_a) % (2.0 * math.pi)) < d_ccw - 0.05
        for oa in other_angs
    )
    has_cw = any(
        0.05 < ((ang_a - oa) % (2.0 * math.pi)) < d_cw - 0.05
        for oa in other_angs
    )
    if not (has_ccw and has_cw):
      continue

    def _can_swap(nb_i: int, piv: int = c_idx) -> bool:
      bd = mol.GetBondBetweenAtoms(piv, nb_i)
      if (
          bd is None
          or bd.IsInRing()
          or bd.GetBondType() != Chem.BondType.SINGLE
          or bd.GetStereo() in _ATROPISOMER_STEREOS
      ):
        return False
      return not any(
          b_n.GetBondType() == Chem.BondType.DOUBLE
          for b_n in mol.GetAtomWithIdx(nb_i).GetBonds()
      )

    x_shift = next((x for x in (a_idx, b_idx) if _can_swap(x)), None)
    x_other = next((x for x in other_nbrs if _can_swap(x)), None)
    if x_shift is None or x_other is None:
      continue
    sub_s = _collect_bond_subtree(mol, c_idx, x_shift, n_atoms)
    sub_o = _collect_bond_subtree(mol, c_idx, x_other, n_atoms)
    if sub_s & sub_o:
      continue
    ps_pt = conf.GetAtomPosition(x_shift)
    po_pt = conf.GetAtomPosition(x_other)
    dtheta = math.atan2(po_pt.y - pc.y, po_pt.x - pc.x) - math.atan2(
        ps_pt.y - pc.y, ps_pt.x - pc.x
    )
    cand = _get_conf_pts_2d(mol, n_atoms)
    for atoms_sub, ang_rot in ((sub_s, dtheta), (sub_o, -dtheta)):
      ca, sa = math.cos(ang_rot), math.sin(ang_rot)
      for idx_s in atoms_sub:
        dx_s, dy_s = cand[idx_s][0] - pc.x, cand[idx_s][1] - pc.y
        cand[idx_s] = (
            pc.x + ca * dx_s - sa * dy_s,
            pc.y + sa * dx_s + ca * dy_s,
        )
    edges = [
        (b_e.GetBeginAtomIdx(), b_e.GetEndAtomIdx()) for b_e in mol.GetBonds()
    ]
    clash = False
    for i_e, (u1, v1) in enumerate(edges):
      for u2, v2 in edges[i_e + 1 :]:
        if (
            u1 not in (u2, v2)
            and v1 not in (u2, v2)
            and intersect_segments(cand[u1], cand[v1], cand[u2], cand[v2])
        ):
          clash = True
          break
      if clash:
        break
    if clash:
      continue
    adj = set(edges) | {(v_e, u_e) for u_e, v_e in edges}
    moved = sub_s | sub_o
    too_close = False
    for i_a in moved:
      for j_a in range(n_atoms):
        if (
            i_a != j_a
            and (i_a, j_a) not in adj
            and math.hypot(
                cand[i_a][0] - cand[j_a][0], cand[i_a][1] - cand[j_a][1]
            )
            < 0.72
        ):
          too_close = True
          break
      if too_close:
        break
    if too_close:
      continue
    _set_conf_pts_2d(mol, {idx_u: cand[idx_u] for idx_u in moved})


class MoleculeCanvas:

  def __init__(
      self,
      smi,
      name,
      rot_deg=0.0,
      add_h=None,
      coord_map=None,
      w_cm=3.6,
      canvas_pos=(0.0, 0.0),
      w_px=350,
      h_px=300,
      atom_labels=None,
      intramol_targets=None,
      kekule_bonds=None,
      ref_canvas=None,
      flip_y: bool = False,
      bond_shift_pairs=None,
      is_product: bool = False,
  ):
    self.smi = smi
    self.name = name
    self.rot_deg = float(rot_deg)
    self.flip_y = bool(flip_y)
    self.is_product = bool(is_product)
    self.w_cm = float(w_cm)
    self.canvas_pos = np.array(canvas_pos, dtype=float)
    self.w_px = w_px
    self.h_px = h_px
    self.atom_labels = dict(atom_labels) if atom_labels else None

    mol = _parse_smiles_preserving_mapped_hs(smi)
    if mol is None:
      raise ValueError(f"Invalid SMILES: {smi}")
    self.map_to_idx = {
        a.GetAtomMapNum(): a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomMapNum() > 0
    }

    added_h_pairs = []
    if add_h is not None:
      # First compute 2D coords on the heavy-atom scaffold so adding H never distorts the scaffold
      rdDepictor.Compute2DCoords(mol)
      regularize_nonconvex_small_rings(mol)
      n_heavy_orig = mol.GetNumAtoms()
      base_pts = _get_conf_pts_2d(mol, n_heavy_orig)
      rw = Chem.RWMol(mol)
      add_h_list = [add_h] if isinstance(add_h, (int, str)) else list(add_h)
      for h_target in add_h_list:
        h_idx_target = _resolve_mol_atom_token(mol, h_target, self.map_to_idx)
        target_atom = rw.GetAtomWithIdx(h_idx_target)
        num_h = target_atom.GetTotalNumHs()
        if num_h > 0:
          target_atom.SetNumExplicitHs(max(0, num_h - 1))
          target_atom.SetNoImplicit(True)
        h_idx = rw.AddAtom(Chem.Atom(1))
        rw.AddBond(h_idx_target, h_idx, Chem.BondType.SINGLE)
        added_h_pairs.append((h_idx_target, h_idx))
      mol = rw.GetMol()
      Chem.SanitizeMol(mol)

      if coord_map is None:
        # Resolve any intramolecular attack targets (e.g. base atom attacking an added H)
        h_to_parent = {h_i: p_i for p_i, h_i in added_h_pairs}
        heavy_pairs, h_attackers = _resolve_intramol_heavy_pairs(
            mol,
            intramol_targets,
            self.map_to_idx,
            name,
            h_to_parent=h_to_parent,
            n_heavy_limit=n_heavy_orig,
        )
        base_pts = _unfold_crowded_intramolecular_pockets(
            mol, base_pts, heavy_pairs, n_heavy_orig
        )
        coord_map = _place_added_explicit_hydrogens(
            mol, base_pts, added_h_pairs, h_attackers
        )

    if coord_map is not None:
      pt_map = {}
      for k, v in coord_map.items():
        k_idx = _resolve_mol_atom_token(mol, k, self.map_to_idx)
        if isinstance(v, (tuple, list)):
          pt_map[k_idx] = Point2D(float(v[0]), float(v[1]))
        else:
          pt_map[k_idx] = v
      rdDepictor.Compute2DCoords(mol, coordMap=pt_map)
      if len(pt_map) == mol.GetNumAtoms():
        _set_conf_pts_2d(
            mol, {k_idx: (pt.x, pt.y) for k_idx, pt in pt_map.items()}
        )
    else:
      rdDepictor.Compute2DCoords(mol)
      regularize_nonconvex_small_rings(mol)
      frags = Chem.GetMolFrags(mol)
      if len(frags) > 1:
        conf = mol.GetConformer()
        curr_top_y = None
        for frag in frags:
          xs = [conf.GetAtomPosition(i).x for i in frag]
          ys = [conf.GetAtomPosition(i).y for i in frag]
          mean_x = 0.5 * (min(xs) + max(xs))
          min_y, max_y = min(ys), max(ys)
          if curr_top_y is None:
            for i in frag:
              p = conf.GetAtomPosition(i)
              conf.SetAtomPosition(i, Point3D(p.x - mean_x, p.y, 0.0))
            curr_top_y = max_y
          else:
            dy = (curr_top_y + 2.4) - min_y
            for i in frag:
              p = conf.GetAtomPosition(i)
              conf.SetAtomPosition(i, Point3D(p.x - mean_x, p.y + dy, 0.0))
            curr_top_y = max_y + dy
      elif intramol_targets:
        n_atoms_mol = mol.GetNumAtoms()
        base_pts = _get_conf_pts_2d(mol, n_atoms_mol)
        heavy_pairs, _ = _resolve_intramol_heavy_pairs(
            mol, intramol_targets, self.map_to_idx, name
        )
        new_pts = _unfold_crowded_intramolecular_pockets(
            mol, base_pts, heavy_pairs, n_atoms_mol
        )
        if new_pts != base_pts:
          _set_conf_pts_2d(mol, new_pts)
      if bond_shift_pairs:
        _align_geminal_bond_shift_substituents(
            mol, bond_shift_pairs, self.map_to_idx
        )

    if ref_canvas is not None and coord_map is None:
      try:
        r_mol = ref_canvas.mol
        pair_candidates: List[Tuple[List[int], List[int]]] = []
        if isinstance(kekule_bonds, list):
          for atom_pairs, _ in kekule_bonds:
            valid_pairs = [
                (r_i, p_i)
                for r_i, p_i in atom_pairs
                if r_i < r_mol.GetNumAtoms() and p_i < mol.GetNumAtoms()
            ]
            if len(valid_pairs) >= 3:
              pair_candidates.append((
                  [r_i for r_i, _ in valid_pairs],
                  [p_i for _, p_i in valid_pairs],
              ))
        if not pair_candidates:
          mcs = rdFMCS.FindMCS([r_mol, mol], timeout=2)
          mcs_any = rdFMCS.FindMCS(
              [r_mol, mol],
              bondCompare=rdFMCS.BondCompare.CompareAny,
              ringMatchesRingOnly=True,
              timeout=2,
          )
          if mcs_any.numAtoms > mcs.numAtoms:
            mcs = mcs_any
          if mcs.numAtoms >= 3:
            q_mcs = Chem.MolFromSmarts(mcs.smartsString)
            r_matches = r_mol.GetSubstructMatches(q_mcs, uniquify=False)
            p_matches = mol.GetSubstructMatches(q_mcs, uniquify=False)
            if r_matches and p_matches:
              match_r0 = list(r_matches[0])
              pair_candidates = [(match_r0, list(mp)) for mp in p_matches]

        if pair_candidates:
          r_conf = r_mol.GetConformer()
          p_conf = mol.GetConformer()
          has_chiral = any(
              a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
              for a in mol.GetAtoms()
          )
          best_r = None
          best_c_p = None
          best_c_r = None
          best_score = 1e9
          for match_r, match_p in pair_candidates:
            x_r = np.array(
                [
                    [
                        r_conf.GetAtomPosition(idx).x,
                        r_conf.GetAtomPosition(idx).y,
                    ]
                    for idx in match_r
                ],
                dtype=float,
            )
            c_r = np.mean(x_r, axis=0)
            x_r_c = x_r - c_r
            x_p = np.array(
                [
                    [
                        p_conf.GetAtomPosition(idx).x,
                        p_conf.GetAtomPosition(idx).y,
                    ]
                    for idx in match_p
                ],
                dtype=float,
            )
            c_p = np.mean(x_p, axis=0)
            x_p_c = x_p - c_p
            u_mat, _, vt_mat = np.linalg.svd(x_p_c.T @ x_r_c)
            vt_rot = vt_mat.copy()
            if np.linalg.det(u_mat @ vt_rot) < 0:
              vt_rot[-1, :] *= -1
            r_rot = u_mat @ vt_rot
            rmsd_rot = float(
                np.sqrt(np.mean(np.sum((x_p_c @ r_rot - x_r_c) ** 2, axis=1)))
            )
            if rmsd_rot < best_score:
              best_score = rmsd_rot
              best_r = r_rot
              best_c_p = c_p
              best_c_r = c_r

            if not has_chiral:
              r_o2 = u_mat @ vt_mat
              rmsd_o2 = (
                  float(
                      np.sqrt(
                          np.mean(np.sum((x_p_c @ r_o2 - x_r_c) ** 2, axis=1))
                      )
                  )
                  + 0.02
              )
              if rmsd_o2 < best_score:
                best_score = rmsd_o2
                best_r = r_o2
                best_c_p = c_p
                best_c_r = c_r

          if (
              best_r is not None
              and best_c_p is not None
              and best_c_r is not None
          ):
            for idx_a in range(mol.GetNumAtoms()):
              pt_a = p_conf.GetAtomPosition(idx_a)
              v_a = np.array([pt_a.x, pt_a.y], dtype=float) - best_c_p
              v_new = v_a @ best_r + best_c_r
              p_conf.SetAtomPosition(
                  idx_a, Point3D(float(v_new[0]), float(v_new[1]), 0.0)
              )
            self.rot_deg = float(ref_canvas.rot_deg)
      except (RuntimeError, ValueError, np.linalg.LinAlgError) as exc:
        raise ValueError(
            f"Failed to align product 2D coordinates for {name!r} ({smi!r}):"
            f" {exc}"
        ) from exc

    if isinstance(kekule_bonds, list) and kekule_bonds:
      if ref_canvas is not None and len(kekule_bonds) > 1:
        r_conf = ref_canvas.mol.GetConformer()
        p_conf = mol.GetConformer()
        best_kb = kekule_bonds[0][1]
        best_dist = 1e9
        for atom_pairs, kb_cand in kekule_bonds:
          dists = []
          for r_idx, p_idx in atom_pairs:
            if (
                r_idx < ref_canvas.mol.GetNumAtoms()
                and p_idx < mol.GetNumAtoms()
            ):
              pr = r_conf.GetAtomPosition(r_idx)
              pp = p_conf.GetAtomPosition(p_idx)
              dists.append((pr.x - pp.x) ** 2 + (pr.y - pp.y) ** 2)
          if dists:
            mean_d = float(np.mean(dists))
            if mean_d < best_dist:
              best_dist = mean_d
              best_kb = kb_cand
        kekule_bonds = best_kb
      else:
        kekule_bonds = kekule_bonds[0][1]

    if kekule_bonds:
      for b in mol.GetBonds():
        u_i, v_i = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        edge_key = (min(u_i, v_i), max(u_i, v_i))
        if edge_key in kekule_bonds:
          b.SetBondType(kekule_bonds[edge_key])
        b.SetIsAromatic(False)
      for a in mol.GetAtoms():
        a.SetIsAromatic(False)
    else:
      try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
      except (Chem.KekulizeException, ValueError) as exc:
        raise ValueError(
            f"Failed to kekulize molecule {name!r} ({smi!r}): {exc}"
        ) from exc

    if any(
        b.GetBondType() != Chem.BondType.SINGLE
        and b.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH)
        for b in mol.GetBonds()
    ):
      _clear_stale_single_bond_wedges(mol)

    if self.flip_y and ref_canvas is None:
      if not _is_planar_achiral_mol(mol):
        raise ValueError(
            f"Cannot mirror chiral or wedged molecule {name!r} ({smi!r}) with"
            " flip_y=True."
        )
      if mol.GetNumConformers() > 0:
        conf_f = mol.GetConformer()
        for idx_f in range(mol.GetNumAtoms()):
          pf = conf_f.GetAtomPosition(idx_f)
          conf_f.SetAtomPosition(idx_f, Point3D(float(pf.x), float(-pf.y), 0.0))
        self.rot_deg = (-self.rot_deg) % 360.0

    d2d = rdMolDraw2D.MolDraw2DCairo(w_px, h_px)
    opts = d2d.drawOptions()
    opts.clearBackground = True
    opts.bondLineWidth = 2.0
    if self.rot_deg != 0.0:
      opts.rotate = self.rot_deg
    svg_atom_labels = {}
    if self.atom_labels:
      for k, v in self.atom_labels.items():
        k_idx = _resolve_mol_atom_token(mol, k, self.map_to_idx)
        opts.atomLabels[k_idx] = str(v)
        svg_atom_labels[k_idx] = str(v)
    # Strip atom-map numbers on the rendered copy so ':1' tags do not clutter atom symbols
    draw_mol = Chem.Mol(mol)
    for a in draw_mol.GetAtoms():
      if a.GetAtomMapNum() != 0:
        a.SetAtomMapNum(0)
    # Native proportions - no font size restrictions
    d2d.DrawMolecule(draw_mol)
    d2d.FinishDrawing()

    png_path = os.path.join(CACHE_DIR, f"{name}.png")
    _atomic_write_bytes(png_path, d2d.GetDrawingText())
    self.png_path = png_path
    self.mol = mol

    scale = self.w_cm / float(w_px)
    self.scale = scale
    self.h_cm = float(h_px) * scale
    raw_positions = {}
    for i in range(mol.GetNumAtoms()):
      pt = d2d.GetDrawCoords(i)
      cx = self.canvas_pos[0] + (pt.x - w_px / 2.0) * scale
      cy = self.canvas_pos[1] - (pt.y - h_px / 2.0) * scale
      raw_positions[i] = np.array([cx, cy])

    bond_lengths_cm = [
        float(
            np.linalg.norm(
                raw_positions[b.GetBeginAtomIdx()]
                - raw_positions[b.GetEndAtomIdx()]
            )
        )
        for b in mol.GetBonds()
    ]
    median_bond_cm = (
        float(np.median(bond_lengths_cm)) if bond_lengths_cm else 0.80
    )
    glyph_scale = min(1.0, max(0.52, median_bond_cm / 0.82))

    self.atoms = {}
    for i in range(mol.GetNumAtoms()):
      sym = mol.GetAtomWithIdx(i).GetSymbol()
      base_r = (
          0.25
          if (len(sym) > 1 or mol.GetAtomWithIdx(i).GetFormalCharge() != 0)
          else (0.12 if sym == "C" else (0.16 if sym == "H" else 0.19))
      )
      self.atoms[i] = {
          "idx": i,
          "sym": sym,
          "pos": raw_positions[i],
          "glyph_r": base_r * glyph_scale,
      }

    labeled_idxs = set(svg_atom_labels.keys())
    svg_bond_mids = {}
    if any(b.GetBondTypeAsDouble() >= 1.5 for b in mol.GetBonds()):
      svg_bond_mids = _extract_rdkit_svg_bond_offsets_cm(
          draw_mol,
          w_px,
          h_px,
          scale,
          self.canvas_pos,
          rot_deg=self.rot_deg,
          atom_labels_idx=svg_atom_labels,
      )

    self.bonds = []
    center_of_mass = np.mean([a["pos"] for a in self.atoms.values()], axis=0)
    for b in mol.GetBonds():
      a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
      p1 = self.atoms[a1]["pos"]
      p2 = self.atoms[a2]["pos"]
      v = p2 - p1
      l = np.linalg.norm(v)
      if l == 0:
        continue
      u = v / l
      n = np.array([-u[1], u[0]])
      mid = (p1 + p2) / 2.0
      if np.dot(mid + 0.1 * n - center_of_mass, n) < 0:
        n = -n
      b_idx = b.GetIdx()
      if b.GetBondTypeAsDouble() >= 1.5 and b_idx in svg_bond_mids:
        raw_offs = sorted(
            float(np.dot(m_pt - mid, n)) for m_pt in svg_bond_mids[b_idx]
        )
        stroke_offsets_cm = []
        for off_val in raw_offs:
          if (
              not stroke_offsets_cm
              or abs(off_val - stroke_offsets_cm[-1]) > 0.005
          ):
            stroke_offsets_cm.append(off_val)
      else:
        stroke_offsets_cm = _compute_bond_stroke_offsets_cm(
            mol, b, self.atoms, n, l, labeled_idxs
        )
      stroke_segments = [
          (p1 + off * n, p2 + off * n) for off in stroke_offsets_cm
      ]
      self.bonds.append({
          "a1": a1,
          "a2": a2,
          "p1": p1,
          "p2": p2,
          "mid": mid,
          "tangent": u,
          "norm": n,
          "order": b.GetBondTypeAsDouble(),
          "stroke_offsets_cm": stroke_offsets_cm,
          "stroke_segments": stroke_segments,
      })

  def find_bond(self, a1, a2):
    for b in self.bonds:
      if (b["a1"] == a1 and b["a2"] == a2) or (b["a1"] == a2 and b["a2"] == a1):
        return b
    return None


def _extract_rdkit_svg_bond_offsets_cm(
    draw_mol: Chem.Mol,
    w_px: int,
    h_px: int,
    scale: float,
    canvas_pos: np.ndarray | Tuple[float, float],
    rot_deg: float = 0.0,
    atom_labels_idx: Optional[Dict[int, str]] = None,
) -> Dict[int, List[np.ndarray]]:
  """Extracts drawn bond stroke midpoints (cm) by bond index via MolDraw2DSVG.

  Args:
      draw_mol: Kekulized RDKit molecule with atom-map numbers cleared.
      w_px: Canvas width in pixels.
      h_px: Canvas height in pixels.
      scale: Centimeters-per-pixel scale factor (`w_cm / w_px`).
      canvas_pos: `(x_cm, y_cm)` center of the molecule canvas.
      rot_deg: Drawing rotation angle in degrees.
      atom_labels_idx: Optional mapping from 0-based atom index to custom text
        label.

  Returns:
      Mapping from 0-based RDKit bond index to a list of `(x_cm, y_cm)` midpoint
      vectors for each parallel bond stroke rendered by RDKit.
  """
  d2d_svg = rdMolDraw2D.MolDraw2DSVG(w_px, h_px)
  opts = d2d_svg.drawOptions()
  opts.clearBackground = True
  opts.bondLineWidth = 2.0
  if rot_deg != 0.0:
    opts.rotate = rot_deg
  if atom_labels_idx:
    for k_idx, v in atom_labels_idx.items():
      opts.atomLabels[k_idx] = str(v)
  d2d_svg.DrawMolecule(draw_mol)
  d2d_svg.FinishDrawing()
  svg_text = d2d_svg.GetDrawingText()

  cx0, cy0 = canvas_pos
  bond_mids = {}
  for path_match in re.finditer(r"<path\b[^>]*>", svg_text):
    tag = path_match.group(0)
    m_cls = re.search(r"\bbond-(\d+)\b", tag)
    if not m_cls:
      continue
    m_d = re.search(r"\bd=['\"]([^'\"]+)['\"]", tag)
    if not m_d:
      continue
    nums = [
        float(x) for x in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", m_d.group(1))
    ]
    if len(nums) >= 4:
      mx = 0.5 * (nums[0] + nums[2])
      my = 0.5 * (nums[1] + nums[3])
      mid_cm = np.array([
          cx0 + (mx - w_px / 2.0) * scale,
          cy0 - (my - h_px / 2.0) * scale,
      ])
      bond_mids.setdefault(int(m_cls.group(1)), []).append(mid_cm)
  return bond_mids


def _compute_bond_stroke_offsets_cm(
    mol: Chem.Mol,
    bond: Chem.Bond,
    atoms_dict: Dict[int, Dict[str, Any]],
    norm_vec: np.ndarray,
    bond_len_cm: float,
    labeled_idxs: Optional[Set[int]] = None,
) -> List[float]:
  """Computes signed normal offsets (cm) of RDKit's parallel bond strokes.

  Args:
      mol: Kekulized RDKit molecule.
      bond: RDKit bond whose stroke offsets are being computed.
      atoms_dict: Mapping from atom index to canvas atom metadata (`pos`, etc.).
      norm_vec: Unit normal vector of the bond on the canvas.
      bond_len_cm: Drawn bond length in centimeters.
      labeled_idxs: Optional set of atom indices carrying custom `atomLabels`.

  Returns:
      Sorted list of signed perpendicular offsets in centimeters along
      `norm_vec`
      for each parallel stroke of `bond`.
  """
  order = bond.GetBondTypeAsDouble()
  mult_off = 0.15 * float(bond_len_cm)
  if order < 1.5:
    return [0.0]
  if order >= 2.5:
    return [-mult_off, 0.0, mult_off]
  # Cumulenes, allenes, ketenes, and CO2 (sp centers with adjacent =/# bonds)
  for end_atom in (bond.GetBeginAtom(), bond.GetEndAtom()):
    for nb_b in end_atom.GetBonds():
      if nb_b.GetIdx() != bond.GetIdx() and nb_b.GetBondTypeAsDouble() >= 1.5:
        return [-0.5 * mult_off, 0.5 * mult_off]
  a1 = bond.GetBeginAtomIdx()
  a2 = bond.GetEndAtomIdx()
  mid = (atoms_dict[a1]["pos"] + atoms_dict[a2]["pos"]) / 2.0
  if bond.IsInRing():
    for ring in mol.GetRingInfo().AtomRings():
      if a1 in ring and a2 in ring:
        rc = np.mean([atoms_dict[k]["pos"] for k in ring], axis=0)
        sign = 1.0 if float(np.dot(rc - mid, norm_vec)) >= 0.0 else -1.0
        return sorted([0.0, sign * mult_off])
  labeled = labeled_idxs or set()
  skeletal_c_1nbr = []
  endpoint_nbrs = []
  for atom_obj, other_idx in (
      (bond.GetBeginAtom(), a2),
      (bond.GetEndAtom(), a1),
  ):
    idx = atom_obj.GetIdx()
    has_label = (
        atom_obj.GetSymbol() != "C"
        or atom_obj.GetFormalCharge() != 0
        or atom_obj.GetNumExplicitHs() > 0
        or atom_obj.GetIsotope() != 0
        or atom_obj.GetNumRadicalElectrons() > 0
        or idx in labeled
    )
    nbrs = [
        nb.GetIdx()
        for nb in atom_obj.GetNeighbors()
        if nb.GetIdx() != other_idx
    ]
    endpoint_nbrs.append((idx, nbrs, has_label))
    if not has_label and len(nbrs) == 1:
      skeletal_c_1nbr.append((idx, nbrs[0]))
  if skeletal_c_1nbr:
    c_idx, nb_idx = skeletal_c_1nbr[0]
    v_nb = atoms_dict[nb_idx]["pos"] - atoms_dict[c_idx]["pos"]
    sign = 1.0 if float(np.dot(v_nb, norm_vec)) >= 0.0 else -1.0
    return sorted([0.0, sign * mult_off])
  if any(not nbrs for _, nbrs, _ in endpoint_nbrs):
    return [-0.5 * mult_off, 0.5 * mult_off]
  net_proj = 0.0
  for idx, nbrs, _ in endpoint_nbrs:
    for nb_idx in nbrs:
      net_proj += float(
          np.dot(atoms_dict[nb_idx]["pos"] - atoms_dict[idx]["pos"], norm_vec)
      )
  sign = 1.0 if net_proj >= 0.0 else -1.0
  return sorted([0.0, sign * mult_off])


def bond_outer_offset_cm(
    bond: Optional[Dict[str, Any]], normal_dir: np.ndarray
) -> float:
  """Returns outermost parallel stroke offset (cm) from bond centerline along normal_dir.

  Args:
      bond: Bond geometry dictionary from `MoleculeCanvas.bonds`.
      normal_dir: Unit normal vector pointing toward the arrow's side of the
        bond.

  Returns:
      Non-negative perpendicular offset in centimeters from the bond centerline
      to the outermost drawn bond stroke on the `normal_dir` side.
  """
  if not bond or "stroke_offsets_cm" not in bond:
    return 0.0
  sign = 1.0 if float(np.dot(normal_dir, bond["norm"])) >= 0.0 else -1.0
  return max(0.0, max(sign * off for off in bond["stroke_offsets_cm"]))


def bezier_turns_ccw(
    p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray
) -> bool:
  """Returns True when the cubic Bézier curve turns CCW at its arrowhead tip (t=0.90).

  Args:
      p0: Start point in centimeters.
      p1: First control point in centimeters.
      p2: Second control point in centimeters.
      p3: End point in centimeters.

  Returns:
      True if the curve's signed curvature at `t=0.90` is counterclockwise (or,
      when the tip segment is nearly collinear, if the midpoint bows below the
      chord).
  """
  t = 0.90
  b1_t = (
      3.0 * (1.0 - t) ** 2 * (p1 - p0)
      + 6.0 * (1.0 - t) * t * (p2 - p1)
      + 3.0 * t**2 * (p3 - p2)
  )
  b2_t = 6.0 * (1.0 - t) * (p2 - 2.0 * p1 + p0) + 6.0 * t * (p3 - 2.0 * p2 + p1)
  cross_tip = float(b1_t[0] * b2_t[1] - b1_t[1] * b2_t[0])
  if abs(cross_tip) > 1e-5:
    return cross_tip > 0.0
  v_chord = p3 - p0
  mid_pt = 0.125 * p0 + 0.375 * p1 + 0.375 * p2 + 0.125 * p3
  v_mid = mid_pt - p0
  cross_z = float(v_chord[0] * v_mid[1] - v_chord[1] * v_mid[0])
  return cross_z < -1e-6


def ccw(A, B, C):
  return (C[1] - A[1]) * (B[0] - A[0]) > (B[1] - A[1]) * (C[0] - A[0])


def intersect_segments(A, B, C, D):
  return ccw(A, C, D) != ccw(B, C, D) and ccw(A, B, C) != ccw(A, B, D)


def sample_bezier_pts(
    p0, p1, p2, p3, n: int = 25, t_min: float = 0.05, t_max: float = 0.95
):
  t = np.linspace(t_min, t_max, n)[:, None]
  return (
      (1 - t) ** 3 * p0
      + 3 * (1 - t) ** 2 * t * p1
      + 3 * (1 - t) * t**2 * p2
      + t**3 * p3
  )


def dist_pt_to_seg(p, a, b) -> float:
  vx, vy = float(b[0] - a[0]), float(b[1] - a[1])
  l2 = vx * vx + vy * vy
  if l2 < 1e-8:
    return math.hypot(float(p[0] - a[0]), float(p[1] - a[1]))
  t = max(
      0.0,
      min(1.0, (float(p[0] - a[0]) * vx + float(p[1] - a[1]) * vy) / l2),
  )
  return math.hypot(
      float(p[0] - (a[0] + t * vx)), float(p[1] - (a[1] + t * vy))
  )


_ARROW_TYPES: Final = (
    "attack",
    "bond_cleave",
    "bond_to_bond",
    "bond_to_atom",
    "lone_pair_to_bond",
)


def _endpoint_y_cm(canvas: Any, tok: Any) -> float:
  """Returns canvas y-coordinate (cm) of an atom or bond-midpoint token."""
  if canvas is None or tok is None:
    return 0.0
  try:
    if isinstance(tok, (list, tuple)) and len(tok) == 2:
      a1 = resolve_atom_token(canvas, tok[0])
      a2 = resolve_atom_token(canvas, tok[1])
      return 0.5 * float(
          canvas.atoms[a1]["pos"][1] + canvas.atoms[a2]["pos"][1]
      )
    a_idx = resolve_atom_token(canvas, tok)
    return float(canvas.atoms[a_idx]["pos"][1])
  except (KeyError, ValueError, TypeError, IndexError):
    return 0.0


def _arrow_route_sort_key(
    a_spec: dict[str, Any], mol_map: dict[str, Any]
) -> tuple[int, float]:
  """Orders tight local shifts first, attacks top-to-bottom, cleavages last."""
  atype = a_spec["type"]
  if atype == "bond_cleave":
    m = mol_map.get(a_spec.get("mol", "r1"))
    return (3, -_endpoint_y_cm(m, a_spec.get("bond")))
  if atype == "lone_pair_to_bond":
    m = mol_map.get(a_spec.get("mol", "r1"))
    y_sum = _endpoint_y_cm(m, a_spec.get("atom")) + _endpoint_y_cm(
        m, a_spec.get("bond")
    )
    return (0, -y_sum)
  if atype == "bond_to_bond":
    sm = a_spec.get("src_mol") or a_spec.get("mol", "r1")
    tm = a_spec.get("tgt_mol") or a_spec.get("mol", sm)
    b1 = a_spec.get("b1") or a_spec.get("src_bond")
    b2 = a_spec.get("b2") or a_spec.get("tgt_bond")
    y_sum = _endpoint_y_cm(mol_map.get(sm), b1) + _endpoint_y_cm(
        mol_map.get(tm), b2
    )
    return (0 if sm == tm else 2, -y_sum)
  if atype == "attack":
    sm, s_tok = a_spec["src"]
    tm, t_tok = a_spec["tgt"]
    y_sum = _endpoint_y_cm(mol_map.get(sm), s_tok) + _endpoint_y_cm(
        mol_map.get(tm), t_tok
    )
    return (1 if sm == tm else 2, -y_sum)
  if atype == "bond_to_atom":
    sm = a_spec.get("src_mol", "r1")
    tm = a_spec.get("tgt_mol", "r2")
    y_sum = _endpoint_y_cm(mol_map.get(sm), a_spec.get("bond")) + (
        _endpoint_y_cm(mol_map.get(tm), a_spec.get("tgt_atom"))
    )
    return (1 if sm == tm else 2, -y_sum)
  return (2, 0.0)


def _cross2d(u_arr: np.ndarray, v_arr: np.ndarray) -> np.ndarray:
  """Computes 2D cross product `u_x * v_y - u_y * v_x` over trailing axis."""
  return u_arr[..., 0] * v_arr[..., 1] - u_arr[..., 1] * v_arr[..., 0]


def _polyline_min_dist_cm(
    poly1: np.ndarray, poly2: np.ndarray, fast_cutoff_cm: float = 0.0
) -> float:
  """Returns minimum distance (cm) between 2D polylines (`0.0` if crossing)."""
  if len(poly1) < 2 or len(poly2) < 2:
    return float("inf")
  if fast_cutoff_cm > 0.0:
    min1 = np.min(poly1, axis=0)
    max1 = np.max(poly1, axis=0)
    min2 = np.min(poly2, axis=0)
    max2 = np.max(poly2, axis=0)
    dx_b = max(0.0, float(min1[0] - max2[0]), float(min2[0] - max1[0]))
    dy_b = max(0.0, float(min1[1] - max2[1]), float(min2[1] - max1[1]))
    if dx_b * dx_b + dy_b * dy_b >= fast_cutoff_cm * fast_cutoff_cm:
      return math.hypot(dx_b, dy_b)

  a = poly1[:-1]
  b = poly1[1:]
  c = poly2[:-1]
  d = poly2[1:]

  da = d[None, :, :] - a[:, None, :]
  ca = c[None, :, :] - a[:, None, :]
  db = d[None, :, :] - b[:, None, :]
  cb = c[None, :, :] - b[:, None, :]
  ba = (b - a)[:, None, :]
  c1 = _cross2d(da, ca)
  c2 = _cross2d(db, cb)
  c3 = _cross2d(ba, ca)
  c4 = _cross2d(ba, da)
  if bool(np.any(((c1 > 0.0) != (c2 > 0.0)) & ((c3 > 0.0) != (c4 > 0.0)))):
    return 0.0

  v = d - c
  v_l2 = np.maximum(np.sum(v * v, axis=-1), 1e-12)
  t1 = np.clip(
      np.sum((poly1[:, None, :] - c[None, :, :]) * v[None, :, :], axis=-1)
      / v_l2[None, :],
      0.0,
      1.0,
  )
  proj1 = c[None, :, :] + t1[:, :, None] * v[None, :, :]
  d1_sq = float(np.min(np.sum((poly1[:, None, :] - proj1) ** 2, axis=-1)))

  u = b - a
  u_l2 = np.maximum(np.sum(u * u, axis=-1), 1e-12)
  t2 = np.clip(
      np.sum((poly2[:, None, :] - a[None, :, :]) * u[None, :, :], axis=-1)
      / u_l2[None, :],
      0.0,
      1.0,
  )
  proj2 = a[None, :, :] + t2[:, :, None] * u[None, :, :]
  d2_sq = float(np.min(np.sum((poly2[:, None, :] - proj2) ** 2, axis=-1)))
  return math.sqrt(min(d1_sq, d2_sq))


def arrowhead_wing_segments(
    p0: np.ndarray,
    p1: np.ndarray,
    p2: np.ndarray,
    p3: np.ndarray,
    electrons: int = 2,
) -> list[tuple[np.ndarray, np.ndarray]]:
  """Returns segments approximating the rendered Stealth/harpoon arrowhead."""
  poly_full = sample_bezier_pts(p0, p1, p2, p3, n=32, t_min=0.0, t_max=1.0)
  arc_l = float(np.sum(np.linalg.norm(poly_full[1:] - poly_full[:-1], axis=1)))
  if electrons == 1:
    head_l = 0.185 if arc_l < 0.45 else 0.245
    half_w = 0.088 if arc_l < 0.45 else 0.112
  else:
    head_l = 0.135 if arc_l < 0.48 else 0.200
    half_w = 0.049 if arc_l < 0.48 else 0.065

  tangent = np.array([1.0, 0.0])
  norm_t = 0.0
  for cand_t in (poly_full[-1] - poly_full[-3], p3 - p2, p3 - p1, p3 - p0):
    norm_c = float(np.linalg.norm(cand_t))
    if norm_c >= 1e-6:
      tangent = cand_t
      norm_t = norm_c
      break
  u_t = tangent / max(norm_t, 1e-6)
  u_n = np.array([-u_t[1], u_t[0]])

  p_tip = p3
  p_base = p_tip - head_l * u_t
  wing_sweep = 0.12 * head_l
  p_left = (p_base - wing_sweep * u_t) + half_w * u_n
  p_right = (p_base - wing_sweep * u_t) - half_w * u_n
  p_notch = p_base + 0.20 * head_l * u_t

  if electrons == 1:
    barb_pt = p_right if bezier_turns_ccw(p0, p1, p2, p3) else p_left
    return [(p_tip, barb_pt), (barb_pt, p_notch)]
  return [
      (p_tip, p_left),
      (p_tip, p_right),
      (p_left, p_notch),
      (p_right, p_notch),
  ]


def _wing_bond_penalty(
    wing_segs: list[tuple[np.ndarray, np.ndarray]],
    b_segs: list[tuple[np.ndarray, np.ndarray]],
    min_dist_cm: float = 0.032,
) -> int:
  """Returns crossing penalty count between Stealth wings and bond strokes."""
  for w1, w2 in wing_segs:
    for bp1, bp2 in b_segs:
      if intersect_segments(w1, w2, bp1, bp2):
        return 4
      if (
          dist_pt_to_seg(w1, bp1, bp2) < min_dist_cm
          or dist_pt_to_seg(w2, bp1, bp2) < min_dist_cm
          or dist_pt_to_seg(bp1, w1, w2) < min_dist_cm
          or dist_pt_to_seg(bp2, w1, w2) < min_dist_cm
      ):
        return 2
  return 0


def _wings_hit_polyline(
    wing_segs: list[tuple[np.ndarray, np.ndarray]],
    poly: np.ndarray,
    min_cm: float = 0.045,
) -> bool:
  """Returns True if any arrowhead wing segment is within `min_cm` of `poly`."""
  for w1, w2 in wing_segs:
    if _polyline_min_dist_cm(np.vstack([w1, w2]), poly) < min_cm:
      return True
  return False


def _prior_arrow_penalty(
    poly: np.ndarray,
    wing_segs: list[tuple[np.ndarray, np.ndarray]],
    prior_arrows: list[np.ndarray] | None,
    *,
    ctrl_pts: (
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None
    ) = None,
) -> tuple[int, float]:
  """Scores candidate polyline and head against already-routed arrows."""
  if not prior_arrows:
    return 0, 0.0
  poly_chk = (
      sample_bezier_pts(*ctrl_pts, n=50, t_min=0.0, t_max=1.0)
      if ctrl_pts is not None
      else poly
  )
  cross_cnt = 0
  soft_pen = 0.0
  for prev_poly in prior_arrows:
    d_arr = _polyline_min_dist_cm(poly_chk, prev_poly, fast_cutoff_cm=0.26)
    if d_arr <= 1e-6:
      cross_cnt += 8
    elif d_arr < 0.060:
      cross_cnt += 5
    elif d_arr < 0.115:
      soft_pen += 260.0 * (0.115 - d_arr) / 0.055
    if d_arr < 0.25 and _wings_hit_polyline(wing_segs, prev_poly, min_cm=0.045):
      cross_cnt += 4
  return cross_cnt, soft_pen


def _signs_for_side(n_base: np.ndarray, side: Any) -> list[float]:
  """Resolves optional directional `side` hint to ordered normal signs."""
  if side is None:
    return [1.0, -1.0]
  if side in ("up", "top", 1, 1.0):
    return [1.0 if n_base[1] > 0 else -1.0]
  if side in ("down", "bottom", -1, -1.0):
    return [-1.0 if n_base[1] > 0 else 1.0]
  if side == "left":
    return [1.0 if n_base[0] < 0 else -1.0]
  if side == "right":
    return [1.0 if n_base[0] > 0 else -1.0]
  return [1.0, -1.0]


def _scene_obstacle_penalty(
    poly: np.ndarray,
    wing_segs: list[tuple[np.ndarray, np.ndarray]],
    scene_canvases: list[Any],
    *,
    owner_mol: Any = None,
    exclude_bonds: tuple[tuple[int, int], ...] = (),
    endpoint_atom_idx: int | None = None,
    owner_glyph_pad: float = 0.005,
    other_glyph_pad: float = 0.005,
    glyph_weight: int = 8,
    min_bond_cm: float = 0.08,
    bond_weight: int = 5,
    end_bond_cm: float | None = None,
    end_bond_weight: int = 0,
    trim_frac: float = 0.0,
    wing_cm: float = 0.028,
) -> int:
  """Scores `poly` and `wing_segs` against atom glyphs and bond strokes."""
  excl_sets = [frozenset(bp) for bp in exclude_bonds]
  cross = 0
  for c_obj in scene_canvases:
    is_owner = c_obj is owner_mol
    for a_i, ainfo in c_obj.atoms.items():
      pad = (
          0.0
          if (is_owner and a_i == endpoint_atom_idx)
          else (owner_glyph_pad if is_owner else other_glyph_pad)
      )
      d_glyph = float(np.min(np.linalg.norm(poly - ainfo["pos"], axis=1)))
      if d_glyph < ainfo["glyph_r"] + pad:
        cross += glyph_weight
    for ob in c_obj.bonds:
      is_excl = is_owner and frozenset((ob["a1"], ob["a2"])) in excl_sets
      if is_excl and end_bond_cm is None:
        continue
      ob_segs = ob["stroke_segments"]
      if not is_excl:
        cross += _wing_bond_penalty(wing_segs, ob_segs, min_dist_cm=wing_cm)
      req_cm = (
          end_bond_cm if (is_excl and end_bond_cm is not None) else min_bond_cm
      )
      w_pen = end_bond_weight if is_excl else bond_weight
      for s1_pt, s2_pt in ob_segs:
        bp1 = s1_pt + trim_frac * (s2_pt - s1_pt) if trim_frac > 0.0 else s1_pt
        bp2 = s2_pt + trim_frac * (s1_pt - s2_pt) if trim_frac > 0.0 else s2_pt
        for s_i in range(len(poly) - 1):
          if (
              intersect_segments(poly[s_i], poly[s_i + 1], bp1, bp2)
              or dist_pt_to_seg(poly[s_i], bp1, bp2) < req_cm
          ):
            cross += w_pen
            break
  return cross


def _open_sector_bisector_degs(nbr_angles: list[float]) -> list[float]:
  """Returns open angular-sector bisectors (degrees in [0, 360))."""
  if not nbr_angles:
    return []
  if len(nbr_angles) == 1:
    a0 = nbr_angles[0]
    return [
        math.degrees(a0 + math.pi) % 360.0,
        math.degrees(a0 + 2.0 * math.pi / 3.0) % 360.0,
        math.degrees(a0 - 2.0 * math.pi / 3.0) % 360.0,
    ]
  angs = sorted(((a + math.pi) % (2.0 * math.pi) - math.pi) for a in nbr_angles)
  out_degs: list[float] = []
  for idx_a, a_curr in enumerate(angs):
    a_next = (
        angs[idx_a + 1] if idx_a + 1 < len(angs) else angs[0] + 2.0 * math.pi
    )
    gap = a_next - a_curr
    if gap >= math.radians(42.0):
      mid = a_curr + 0.5 * gap
      out_degs.append(math.degrees(mid) % 360.0)
      if gap >= math.radians(95.0):
        out_degs.append(math.degrees(mid - math.radians(22.0)) % 360.0)
        out_degs.append(math.degrees(mid + math.radians(22.0)) % 360.0)
  return out_degs


# Half-extent of `\node[fill=white, inner sep=2pt] {\huge $+$}` (~0.17 cm glyph
# half-size + 2pt = 0.07 cm padding) and clearance box around `+` separators.
_PLUS_GLYPH_RADIUS_CM: Final[float] = 0.24
_PLUS_BOX_HALF_X_CM: Final[float] = 0.35
_PLUS_BOX_HALF_Y_CM: Final[float] = 0.34


def hits_plus_sign(
    poly: np.ndarray,
    plus_xs: Sequence[float],
    wing_segs: Sequence[tuple[np.ndarray, np.ndarray]] = (),
) -> bool:
  """Returns True if a curve or arrowhead wing enters a `+` separator box."""
  if not plus_xs:
    return False
  pts = (
      np.vstack([poly, *(np.linspace(w1, w2, 9) for w1, w2 in wing_segs)])
      if wing_segs
      else poly
  )
  return any(
      np.any(
          (np.abs(pts[:, 0] - px) < _PLUS_BOX_HALF_X_CM)
          & (np.abs(pts[:, 1]) < _PLUS_BOX_HALF_Y_CM)
      )
      for px in plus_xs
  )


def audit_arrow_clearances(
    arrows,
    canvases,
    min_clearance_mm=0.5,
    min_arrow_clearance_mm=0.5,
    *,
    plus_xs: Sequence[float] = (),
):
  """Measures arrow clearances against bond strokes, glyphs, and other arrows.

  Samples every Bézier curve that was emitted and compares it against every
  bond stroke segment, every atom glyph disk and `+` separator on the canvas,
  and every other emitted curved arrow in the scene, so bond grazes, glyph
  collisions, and arrow-to-arrow crossings are all audited quantitatively.

  Distances are computed in canvas centimetres and reported in millimetres.
  The audit measures only: it does not move arrows and does not raise.

  Args:
      arrows: Sequence of `(label, (p0, p1, p2, p3))` pairs whose control points
        are canvas positions in cm.
      canvases: MoleculeCanvas objects making up the scene.
      min_clearance_mm: Bond clearance at or above which an arrow counts as
        clear. A 1.0 pt arrow stroke is about 0.35 mm wide, so 0.5 mm of
        measured centre-to-segment gap leaves visible whitespace.
      min_arrow_clearance_mm: Centerline-to-centerline clearance (mm) between
        distinct curved arrows required for `clear` when two or more arrows are
        present (0.0 indicates a direct curve-curve intersection).
      plus_xs: Optional x-coordinates (cm) of `+` separator glyphs at `y = 0`.

  Returns:
      Dict with `all_clear`, `min_clearance_mm`, `min_arrow_clearance_mm`,
      scene-wide `min_bond_clearance_mm`, `min_glyph_clearance_mm`, and
      `closest_arrow_pair_mm` (or `inf` when fewer than 2 arrows exist), plus
      an `arrows` list of per-arrow `label`, `bond_clearance_mm`,
      `glyph_clearance_mm`, `arrow_clearance_mm`, and `clear`.
  """
  segments = []
  for c in canvases:
    for b in c.bonds:
      segments.extend(b["stroke_segments"])
  glyphs = [
      (a["pos"], a["glyph_r"]) for c in canvases for a in c.atoms.values()
  ]
  plus_centers = [np.array([float(px), 0.0]) for px in plus_xs]
  sampled_polys = [
      sample_bezier_pts(p0, p1, p2, p3, n=100) for _, (p0, p1, p2, p3) in arrows
  ]
  full_arrow_polys = [
      sample_bezier_pts(p0, p1, p2, p3, n=100, t_min=0.0, t_max=1.0)
      for _, (p0, p1, p2, p3) in arrows
  ]

  report = {
      "all_clear": True,
      "min_clearance_mm": float(min_clearance_mm),
      "min_arrow_clearance_mm": float(min_arrow_clearance_mm),
      "min_bond_clearance_mm": float("inf"),
      "min_glyph_clearance_mm": float("inf"),
      "closest_arrow_pair_mm": float("inf"),
      "arrows": [],
  }

  for idx_a, (label, _) in enumerate(arrows):
    poly = sampled_polys[idx_a]
    bond_mm = 10.0 * min(
        (min(dist_pt_to_seg(pt, s1, s2) for pt in poly) for s1, s2 in segments),
        default=float("inf"),
    )
    # Negative clearance means the curve enters an atom or `+` glyph.
    glyph_mm = 10.0 * min(
        (
            float(np.min(np.linalg.norm(poly - pos, axis=1))) - r
            for pos, r in glyphs
        ),
        default=float("inf"),
    )
    if plus_centers:
      plus_mm = 10.0 * min(
          float(np.min(np.linalg.norm(full_arrow_polys[idx_a] - pc, axis=1)))
          - _PLUS_GLYPH_RADIUS_CM
          for pc in plus_centers
      )
      glyph_mm = min(glyph_mm, plus_mm)
    arrow_mm = 10.0 * min(
        (
            _polyline_min_dist_cm(
                full_arrow_polys[idx_a], full_arrow_polys[idx_b]
            )
            for idx_b in range(len(full_arrow_polys))
            if idx_b != idx_a
        ),
        default=float("inf"),
    )
    clear = (
        bond_mm >= min_clearance_mm
        and glyph_mm >= 0.0
        and arrow_mm >= min_arrow_clearance_mm
    )

    report["all_clear"] = report["all_clear"] and clear
    report["min_bond_clearance_mm"] = min(
        report["min_bond_clearance_mm"], bond_mm
    )
    report["min_glyph_clearance_mm"] = min(
        report["min_glyph_clearance_mm"], glyph_mm
    )
    report["closest_arrow_pair_mm"] = min(
        report["closest_arrow_pair_mm"], arrow_mm
    )
    report["arrows"].append({
        "label": label,
        "bond_clearance_mm": bond_mm,
        "glyph_clearance_mm": glyph_mm,
        "arrow_clearance_mm": arrow_mm,
        "clear": clear,
    })

  return report


def report_arrow_audit(rid, audit, n_specified, *, strict: bool = False):
  """Prints the clearance audit for one reaction and fails on missing arrows.

  Args:
      rid: Reaction id, which attributes the lines when renders run in parallel.
      audit: Report returned by `audit_arrow_clearances`.
      n_specified: Number of arrows the spec asked for.
      strict: When True, also raises `ValueError` if `not audit["all_clear"]`.

  Raises:
      ValueError: If fewer arrows were emitted than specified in the spec, or
        if `strict` is True and any arrow fails the clearance audit.
  """
  n_measured = len(audit["arrows"])
  if n_measured < n_specified:
    msg = (
        f"[AUDIT] {rid}: {n_specified - n_measured} of {n_specified} specified"
        " arrows produced no geometry and are missing from the diagram."
    )
    print(msg)
    raise ValueError(msg)
  for a in audit["arrows"]:
    if a["clear"]:
      continue
    arr_str = (
        f", nearest arrow by {a['arrow_clearance_mm']:.2f} mm"
        if math.isfinite(a["arrow_clearance_mm"])
        else ""
    )
    print(
        f"[AUDIT] {rid}: {a['label']} clears the nearest bond by"
        f" {a['bond_clearance_mm']:.2f} mm and the nearest atom label by"
        f" {a['glyph_clearance_mm']:.2f} mm{arr_str}; wanted at least"
        f" {audit['min_clearance_mm']:.2f} mm bond, 0.00 mm glyph, and"
        f" {audit['min_arrow_clearance_mm']:.2f} mm arrow-to-arrow."
    )
  if strict and not audit["all_clear"]:
    pair_str = (
        f", min_arrow={audit['closest_arrow_pair_mm']:.2f} mm"
        if math.isfinite(audit["closest_arrow_pair_mm"])
        else ""
    )
    raise ValueError(
        f"[AUDIT] {rid}: arrow clearance audit failed"
        f" (min_bond={audit['min_bond_clearance_mm']:.2f} mm,"
        f" min_glyph={audit['min_glyph_clearance_mm']:.2f} mm{pair_str})."
    )
  if n_measured and audit["all_clear"]:
    pair_str = (
        f", {audit['closest_arrow_pair_mm']:.2f} mm between arrows"
        if math.isfinite(audit["closest_arrow_pair_mm"])
        else ""
    )
    print(
        f"[AUDIT] {rid}: {n_measured} arrows clear, closest approach"
        f" {audit['min_bond_clearance_mm']:.2f} mm to a bond and"
        f" {audit['min_glyph_clearance_mm']:.2f} mm to an atom label{pair_str}."
    )


def resolve_atom_token(canvas, token):
  """Resolves a single atom reference to a zero-based atom index.

  Accepts either a plain integer index or an atom-map label of the form
  ``":1"``, matching an atom written as ``[NH2:1]`` in the input SMILES.
  Map labels are the documented syntax for hand-written mechanism JSON,
  because they survive edits to the rest of the structure.

  Args:
      canvas: MoleculeCanvas whose `mol` carries the atom map numbers.
      token: Integer index, numeric string, or ``":N"`` map label.

  Returns:
      The zero-based atom index.

  Raises:
      ValueError: If a map label matches no atom, or matches more than one.
  """
  if isinstance(token, str) and token.startswith(":"):
    want = int(token[1:])
    hits = [
        a.GetIdx() for a in canvas.mol.GetAtoms() if a.GetAtomMapNum() == want
    ]
    if not hits:
      raise ValueError(
          f"Atom map label '{token}' matches no atom. Write the atom as "
          f"[X:{want}] in the SMILES to tag it."
      )
    if len(hits) > 1:
      raise ValueError(
          f"Atom map label '{token}' is ambiguous: matches atoms {hits}."
      )
    return hits[0]
  return int(token)


def resolve_arrow_endpoint(canvas, idx):
  """Resolves an arrow endpoint that may be an atom index or a bond pair.

  Curved-arrow specs use a bare atom index when the electrons come from (or
  land on) a single centre - a lone pair or an empty orbital - and a
  two-element ``(a1, a2)`` bond pair when a sigma or pi bond is itself the
  electron source, e.g. the C-Mg bond of a Grignard, the B-H bond of
  borohydride, or an enolate C=C. In the bond case the arrow tail belongs at
  the bond midpoint.

  Either form may use atom-map labels (``":1"``) instead of raw indices.

  Args:
      canvas: MoleculeCanvas holding laid-out 2D atom positions.
      idx: Atom reference, or a two-element sequence naming a bond.

  Returns:
      Tuple of (position, clearance_radius, neighbour_angles, atom_ids) where
      neighbour_angles are the bearings of substituents to steer clear of and
      atom_ids are the atoms the endpoint is incident to.

  Raises:
      ValueError: If a bond pair does not name exactly two bonded atoms.
  """
  if isinstance(idx, (tuple, list)):
    if len(idx) != 2:
      raise ValueError(
          f"Bond endpoint must name exactly 2 atoms, got {list(idx)}."
      )
    a1 = resolve_atom_token(canvas, idx[0])
    a2 = resolve_atom_token(canvas, idx[1])
    if canvas.mol.GetBondBetweenAtoms(a1, a2) is None:
      raise ValueError(
          f"Bond endpoint ({a1}, {a2}) is not a bond in"
          f" {canvas.mol.GetNumAtoms()}-atom molecule."
      )
    pos = (canvas.atoms[a1]["pos"] + canvas.atoms[a2]["pos"]) / 2.0
    angles = []
    for end, other in ((a1, a2), (a2, a1)):
      for nbr in canvas.mol.GetAtomWithIdx(end).GetNeighbors():
        if nbr.GetIdx() == other:
          continue
        v = canvas.atoms[nbr.GetIdx()]["pos"] - pos
        angles.append(math.atan2(v[1], v[0]))
    # Small fixed clearance: the tail starts just off the bond midpoint,
    # not outside an atom glyph.
    return pos, 0.05, angles, (a1, a2)

  aidx = resolve_atom_token(canvas, idx)
  pos = canvas.atoms[aidx]["pos"]
  angles = []
  for nbr in canvas.mol.GetAtomWithIdx(aidx).GetNeighbors():
    v = canvas.atoms[nbr.GetIdx()]["pos"] - pos
    angles.append(math.atan2(v[1], v[0]))
  return pos, canvas.atoms[aidx]["glyph_r"], angles, (aidx,)


def _ring_exterior_normal(
    canvas, a1: int, a2: int, n_default: np.ndarray
) -> Optional[np.ndarray]:
  """Returns outward-pointing unit normal for a peripheral ring bond, else None."""
  if not hasattr(canvas, "mol") or canvas.mol is None:
    return None
  bond = canvas.mol.GetBondBetweenAtoms(a1, a2)
  if bond is None or not bond.IsInRing():
    return None
  mid = 0.5 * (canvas.atoms[a1]["pos"] + canvas.atoms[a2]["pos"])
  ring_signs = set()
  for ring in canvas.mol.GetRingInfo().AtomRings():
    if a1 in ring and a2 in ring:
      rc = np.mean([canvas.atoms[k]["pos"] for k in ring], axis=0)
      dot_rc = float(np.dot(rc - mid, n_default))
      if abs(dot_rc) > 1e-4:
        ring_signs.add(1 if dot_rc > 0.0 else -1)
  if len(ring_signs) != 1:
    return None
  rc_sign = next(iter(ring_signs))
  return -n_default if rc_sign > 0 else n_default


def compute_attack_bezier(
    src_mol,
    src_idx,
    tgt_mol,
    tgt_idx,
    all_bonds=None,
    plus_x=None,
    side=None,
    lift=None,
    canvases=None,
    prior_arrows=None,
    electrons=2,
):
  """Computes an obstacle-avoiding cubic Bézier for nucleophilic/bond attack.

  Args:
      src_mol: Source `MoleculeCanvas`.
      src_idx: Source atom token/index or `(a1, a2)` bond token/pair.
      tgt_mol: Target `MoleculeCanvas`.
      tgt_idx: Target atom token/index or `(a1, a2)` bond token/pair.
      all_bonds: Optional prebuilt bond dicts with `"mol"` references; populated
        automatically from `canvases` when omitted.
      plus_x: Optional x-coordinate(s) of inter-molecule `+` glyphs to clear.
      side: Optional trajectory preference (`"over"`, `"under"`, `"interior"`).
      lift: Optional target handle-length scale factor.
      canvases: Optional scene `MoleculeCanvas` sequence for obstacle checks.
      prior_arrows: Optional list of already-routed arrow polylines to avoid.
      electrons: Electron count (`2` for standard pair, `1` for fishhook).

  Returns:
      Tuple `(p0, p1, p2, p3)` of 2D control points in cm, or `None` if no
      feasible geometry exists.
  """
  p_src, r_src, src_nbr_angles, src_atom_ids = resolve_arrow_endpoint(
      src_mol, src_idx
  )
  p_tgt, r_tgt, nbr_angles, tgt_atom_ids = resolve_arrow_endpoint(
      tgt_mol, tgt_idx
  )

  b_src = None
  u_bond = None
  n_bond = None
  n_src_ring_ext = None
  if len(src_atom_ids) == 2:
    b_src = src_mol.find_bond(src_atom_ids[0], src_atom_ids[1])
    if b_src:
      u_bond = b_src["tangent"]
      n_bond = b_src["norm"]
      n_src_ring_ext = _ring_exterior_normal(
          src_mol, src_atom_ids[0], src_atom_ids[1], n_bond
      )

  b_tgt = None
  u_tgt_bond = None
  n_tgt_bond = None
  n_tgt_ring_ext = None
  if len(tgt_atom_ids) == 2:
    b_tgt = tgt_mol.find_bond(tgt_atom_ids[0], tgt_atom_ids[1])
    if b_tgt:
      u_tgt_bond = b_tgt["tangent"]
      n_tgt_bond = b_tgt["norm"]
      n_tgt_ring_ext = _ring_exterior_normal(
          tgt_mol, tgt_atom_ids[0], tgt_atom_ids[1], n_tgt_bond
      )

  v_chord = p_tgt - p_src
  chord_len = float(np.linalg.norm(v_chord))
  n_chord = (
      np.array([-v_chord[1] / chord_len, v_chord[0] / chord_len])
      if chord_len > 1e-4
      else np.array([0.0, 1.0])
  )
  chord_deg = math.degrees(math.atan2(v_chord[1], v_chord[0])) % 360.0
  rev_deg = (chord_deg + 180.0) % 360.0

  cand_deps = []
  dx = p_tgt[0] - p_src[0]

  if dx < 0:  # R2 -> R1
    d_range = (
        np.linspace(15.0, 165.0, 14)
        if n_bond is not None
        else np.linspace(95.0, 165.0, 10)
    )
    a_range = (
        np.linspace(15.0, 165.0, 14)
        if n_tgt_bond is not None
        else np.linspace(15.0, 85.0, 10)
    )
  else:  # R1 -> R2
    d_range = (
        np.linspace(15.0, 165.0, 14)
        if n_bond is not None
        else np.linspace(15.0, 85.0, 10)
    )
    a_range = (
        np.linspace(15.0, 165.0, 14)
        if n_tgt_bond is not None
        else np.linspace(95.0, 165.0, 10)
    )
  for d in d_range:
    for a in a_range:
      cand_deps.append((float(d), float(a), "over"))

  if dx < 0:
    d_range_u = (
        np.linspace(-165.0, -15.0, 14)
        if n_bond is not None
        else np.linspace(-165.0, -95.0, 10)
    )
    a_range_u = (
        np.linspace(-165.0, -15.0, 14)
        if n_tgt_bond is not None
        else np.linspace(-85.0, -15.0, 10)
    )
  else:
    d_range_u = (
        np.linspace(-165.0, -15.0, 14)
        if n_bond is not None
        else np.linspace(-85.0, -15.0, 10)
    )
    a_range_u = (
        np.linspace(-165.0, -15.0, 14)
        if n_tgt_bond is not None
        else np.linspace(-165.0, -95.0, 10)
    )
  for d in d_range_u:
    for a in a_range_u:
      cand_deps.append((float(d), float(a), "under"))

  # Direct / interior C-curve trajectories (both ends bow to the same normal side of chord)
  c_mags = (
      (22.0, 32.0, 42.0, 52.0, 66.0, 82.0)
      if src_mol is tgt_mol
      else (22.0, 32.0, 42.0, 52.0)
  )
  for side_sign in (1.0, -1.0):
    for d_mag in c_mags:
      for a_mag in c_mags:
        d_off = side_sign * d_mag
        a_off = -side_sign * a_mag
        cand_deps.append(
            ((chord_deg + d_off) % 360.0, (rev_deg + a_off) % 360.0, "interior")
        )

  # Open-sector bisector trajectories avoid dead zones on vertical/tilted chords
  open_src_degs = _open_sector_bisector_degs(src_nbr_angles)
  open_tgt_degs = _open_sector_bisector_degs(nbr_angles)
  for d_open in open_src_degs:
    for a_open in open_tgt_degs:
      mean_sin = 0.5 * (
          math.sin(math.radians(d_open)) + math.sin(math.radians(a_open))
      )
      m_tag = (
          "interior"
          if (src_mol is tgt_mol and side not in ("over", "under"))
          else ("over" if mean_sin >= 0.0 else "under")
      )
      cand_deps.append((d_open, a_open, m_tag))

  # Full-circle rescue fan so sterically crowded or exocyclic pockets always
  # retain a collision-free escape sector when standard over/under arcs jam.
  full_fan = np.linspace(-165.0, 165.0, 12)
  for d in full_fan:
    for a in full_fan:
      cand_deps.append((float(d), float(a), "rescue"))

  scene_canvases = (
      list(canvases)
      if canvases is not None
      else [m for m in dict.fromkeys((src_mol, tgt_mol)) if m is not None]
  )
  eff_bonds = all_bonds
  if eff_bonds is None:
    eff_bonds = []
    for c_item in scene_canvases:
      for b_item in c_item.bonds:
        eff_bonds.append(dict(b_item, mol=c_item))

  min_scene_y = None
  for m_obj in scene_canvases:
    for ainfo in m_obj.atoms.values():
      ay = float(ainfo["pos"][1])
      min_scene_y = ay if min_scene_y is None else min(min_scene_y, ay)
  y_floor = min(
      -1.85, (min_scene_y if min_scene_y is not None else -1.0) - 0.85
  )

  best_bezier = None
  best_score = -1e9
  target_lift = (
      lift
      if lift is not None
      else (
          0.65
          if (
              plus_x is not None
              and (not isinstance(plus_x, (list, tuple)) or len(plus_x) > 0)
          )
          else 0.45
      )
  )

  for d_deg, a_deg, mode in cand_deps:
    if mode == "rescue" and best_score > -75.0:
      break
    if side is not None and mode not in (side, "rescue"):
      continue

    d_rad = math.radians(d_deg)
    a_rad = math.radians(a_deg)
    u_dep = np.array([math.cos(d_rad), math.sin(d_rad)])
    u_app = np.array([math.cos(a_rad), math.sin(a_rad)])

    # Penalize arrows running nearly parallel to a source/target bond instead of
    # hard-discarding moderate angles, so wrap-around approaches stay available.
    par_pen = 0.0
    if side is not None and mode == "rescue":
      par_pen += 500.0
    elif mode == "rescue":
      par_pen += 90.0
    if u_bond is not None:
      dot_src = abs(float(np.dot(u_dep, u_bond)))
      if dot_src > 0.96:
        continue
      if dot_src > 0.65:
        par_pen += 120.0 * (dot_src - 0.65)
    if u_tgt_bond is not None:
      dot_tgt = abs(float(np.dot(u_app, u_tgt_bond)))
      if dot_tgt > 0.96:
        continue
      if dot_tgt > 0.65:
        par_pen += 120.0 * (dot_tgt - 0.65)

    clr_dep = (
        min([
            abs((d_rad - na + math.pi) % (2 * math.pi) - math.pi)
            for na in src_nbr_angles
        ])
        if src_nbr_angles
        else 1.0
    )
    clr_app = (
        min([
            abs((a_rad - na + math.pi) % (2 * math.pi) - math.pi)
            for na in nbr_angles
        ])
        if nbr_angles
        else 1.0
    )
    if clr_dep < math.radians(20.0) or clr_app < math.radians(20.0):
      continue

    off_src = r_src + 0.06
    off_tgt = r_tgt + 0.05
    sum_off = off_src + off_tgt
    if chord_len > 1e-4 and sum_off > 0.68 * chord_len:
      scale_off = (0.68 * chord_len) / sum_off
      off_src = max(r_src + 0.015, off_src * scale_off)
      off_tgt = max(r_tgt + 0.015, off_tgt * scale_off)

    if n_bond is not None:
      if n_src_ring_ext is not None:
        dep_n = n_src_ring_ext
        if float(np.dot(u_dep, dep_n)) < -0.35:
          par_pen += 2500.0
      else:
        sign_n = 1.0 if np.dot(u_dep, n_bond) >= 0 else -1.0
        dep_n = sign_n * n_bond
      outer_off = bond_outer_offset_cm(b_src, dep_n)
      p0 = (
          p_src + (outer_off + min(0.115, max(0.065, 0.21 * chord_len))) * dep_n
      )
    else:
      p0 = p_src + off_src * u_dep
    if n_tgt_bond is not None:
      if n_tgt_ring_ext is not None:
        app_n = n_tgt_ring_ext
        if float(np.dot(u_app, app_n)) < -0.35:
          par_pen += 2500.0
      else:
        sign_tn = 1.0 if np.dot(u_app, n_tgt_bond) >= 0 else -1.0
        app_n = sign_tn * n_tgt_bond
      outer_tgt_off = bond_outer_offset_cm(b_tgt, app_n)
      p3 = p_tgt + (outer_tgt_off + min(0.105, 0.20 * chord_len)) * app_n
    else:
      p3 = p_tgt + off_tgt * u_app
    if (
        mode != "rescue"
        and chord_len > 1e-4
        and float(np.dot(p3 - p0, v_chord)) < 0.18 * (chord_len**2)
    ):
      continue
    dist = float(np.linalg.norm(p3 - p0))

    c_symmetry = float(np.dot(u_dep, n_chord) * np.dot(u_app, n_chord))

    lift_scales = (
        [0.24, 0.30, 0.35, 0.45, 0.55, 0.70, 0.85, 1.15]
        if mode != "interior"
        else [0.24, 0.32, 0.42, 0.52, 0.65]
    )
    max_l1 = max(3.5, 0.85 * dist)
    max_l2 = max(3.2, 0.80 * dist)
    for lift_scale in lift_scales:
      min_handle = max(0.14, 0.32 * dist) if dist < 1.6 else 0.8
      handle_l1 = min(max_l1, max(min_handle, dist * lift_scale))
      handle_l2 = min(max_l2, max(min_handle, dist * (lift_scale * 0.95)))
      p1 = p0 + handle_l1 * u_dep
      p2 = p3 + handle_l2 * u_app

      n_cross = 0
      poly = sample_bezier_pts(p0, p1, p2, p3, n=45)
      wing_segs = arrowhead_wing_segments(p0, p1, p2, p3, electrons=electrons)
      if float(np.min(poly[:, 1])) < y_floor:
        n_cross += 4

      for mol_obj in scene_canvases:
        for ainfo in mol_obj.atoms.values():
          d_min = float(np.min(np.linalg.norm(poly - ainfo["pos"], axis=1)))
          if d_min < ainfo["glyph_r"] + 0.005:
            n_cross += 6

      if plus_x is not None:
        p_list = plus_x if isinstance(plus_x, (list, tuple)) else [plus_x]
        for px in p_list:
          if (p0[0] - px) * (p3[0] - px) < 0:
            idx_near = np.where(np.abs(poly[:, 0] - px) < 0.35)[0]
            if len(idx_near) > 0:
              y_near = poly[idx_near, 1]
              if np.any(np.abs(y_near) < 0.30):
                n_cross += 10

      min_noninc_d = 0.20
      for b in eff_bonds:
        is_incident_src = (
            b["a1"] in src_atom_ids or b["a2"] in src_atom_ids
        ) and b["mol"] is src_mol
        is_incident_tgt = (
            b["a1"] in tgt_atom_ids or b["a2"] in tgt_atom_ids
        ) and b["mol"] is tgt_mol
        is_src_bond = (
            len(src_atom_ids) == 2
            and b["mol"] is src_mol
            and {b["a1"], b["a2"]} == set(src_atom_ids)
        )
        is_tgt_bond = (
            len(tgt_atom_ids) == 2
            and b["mol"] is tgt_mol
            and {b["a1"], b["a2"]} == set(tgt_atom_ids)
        )
        b_segs = b["stroke_segments"]

        if not (is_src_bond or is_tgt_bond):
          n_cross += _wing_bond_penalty(wing_segs, b_segs, min_dist_cm=0.032)

        if not (is_incident_src or is_incident_tgt):
          for bp1, bp2 in b_segs:
            for s in range(len(poly) - 1):
              if intersect_segments(poly[s], poly[s + 1], bp1, bp2):
                n_cross += 5
                min_noninc_d = 0.0
                break
              d_seg = dist_pt_to_seg(poly[s], bp1, bp2)
              if d_seg < min_noninc_d:
                min_noninc_d = d_seg
              if d_seg < 0.11:
                n_cross += 2
                break
        elif is_src_bond or is_tgt_bond:
          for bp1, bp2 in b_segs:
            for s in range(len(poly) - 1):
              if intersect_segments(poly[s], poly[s + 1], bp1, bp2):
                n_cross += 5
                break
              if dist_pt_to_seg(poly[s], bp1, bp2) < 0.065:
                n_cross += 3
                break
        else:
          ref_pos = p_src if is_incident_src else p_tgt
          for bp1, bp2 in b_segs:
            for s in range(len(poly) - 1):
              d_inc = dist_pt_to_seg(poly[s], bp1, bp2)
              if (
                  intersect_segments(poly[s], poly[s + 1], bp1, bp2)
                  or d_inc < 0.055
              ):
                n_cross += 5
                break
              if np.linalg.norm(poly[s] - ref_pos) > 0.22 and d_inc < 0.10:
                n_cross += 3
                break

      arr_cross, arr_soft_pen = _prior_arrow_penalty(
          poly, wing_segs, prior_arrows, ctrl_pts=(p0, p1, p2, p3)
      )
      n_cross += arr_cross

      dpoly = poly[1:] - poly[:-1]
      turns = dpoly[:-1, 0] * dpoly[1:, 1] - dpoly[:-1, 1] * dpoly[1:, 0]
      has_hook = (
          c_symmetry < 0.04
          or float(np.dot(p2 - p1, p3 - p0)) < 0.0
          or (float(np.min(turns)) < -1e-5 and float(np.max(turns)) > 1e-5)
      )
      ang_dep = abs(
          (d_rad - math.radians(chord_deg) + math.pi) % (2.0 * math.pi)
          - math.pi
      )
      ang_app = abs(
          (a_rad - math.radians(rev_deg) + math.pi) % (2.0 * math.pi) - math.pi
      )
      bow_pen = (
          3.5 * abs(ang_dep - ang_app)
          + 12.0 * max(0.0, ang_dep - 1.35)
          + 12.0 * max(0.0, ang_app - 1.35)
      )

      score = (
          -1000.0 * n_cross
          - par_pen
          - arr_soft_pen
          - bow_pen
          + 5.0 * min(clr_dep, 1.15)
          + 5.0 * min(clr_app, 1.15)
          + 20.0 * min(min_noninc_d, 0.15)
          - 0.8 * abs(lift_scale - target_lift)
      )
      if has_hook:
        score -= 180.0
      if mode == "under":
        score -= 6.0
      elif mode == "interior":
        score += 6.0
      if score > best_score:
        best_score = score
        best_bezier = (p0, p1, p2, p3)

  return best_bezier


def compute_bond_cleave_bezier(
    mol,
    a1_idx,
    a2_idx,
    tgt_idx,
    incoming_from_right=True,
    attack_pos=None,
    side=None,
    attack_dir=None,
    canvases=None,
    prior_arrows=None,
    electrons=2,
):
  """Computes a bond-cleavage cubic Bézier landing on `tgt_idx`.

  Args:
      mol: `MoleculeCanvas` owning the cleaved bond.
      a1_idx: First bond atom token or index.
      a2_idx: Second bond atom token or index.
      tgt_idx: Target atom token or index receiving the electron(s).
      incoming_from_right: Fallback hemisphere preference when no driver exists.
      attack_pos: Optional 2D position of an incoming attacking center.
      side: Optional normal-side override (`"up"`, `"down"`, `"left"`,
        `"right"`).
      attack_dir: Optional incoming tangent vector of a driving arrow.
      canvases: Optional scene `MoleculeCanvas` sequence for obstacle checks.
      prior_arrows: Optional list of already-routed arrow polylines to avoid.
      electrons: Electron count (`2` for heterolysis, `1` for homolysis).

  Returns:
      Tuple `(p0, p1, p2, p3)` of 2D control points in cm, or `None` if bond
      `a1_idx-a2_idx` does not exist on `mol`.
  """
  a1_idx = resolve_atom_token(mol, a1_idx)
  a2_idx = resolve_atom_token(mol, a2_idx)
  tgt_idx = resolve_atom_token(mol, tgt_idx)
  b = mol.find_bond(a1_idx, a2_idx)
  if not b:
    return None
  p_c = (
      mol.atoms[a1_idx]["pos"]
      if a1_idx != tgt_idx
      else mol.atoms[a2_idx]["pos"]
  )
  p_tgt = mol.atoms[tgt_idx]["pos"]
  r_tgt = mol.atoms[tgt_idx]["glyph_r"]

  v = p_tgt - p_c
  bond_len = float(np.linalg.norm(v))
  if bond_len < 1e-6:
    return None
  t = v / bond_len
  n_base = np.array([-t[1], t[0]])
  signs = _signs_for_side(n_base, side)

  scene_canvases = list(canvases) if canvases is not None else [mol]
  pref_f_start = 0.56 if electrons == 1 else 0.28
  f_start_cands = (0.56, 0.46, 0.32) if electrons == 1 else (0.28, 0.36, 0.52)

  cand = []
  com_x = float(np.mean([a["pos"][0] for a in mol.atoms.values()]))
  tgt_has_h_label = mol.atoms[tgt_idx]["sym"] != "C" and (
      mol.mol.GetAtomWithIdx(tgt_idx).GetTotalNumHs() > 0
      or mol.mol.GetAtomWithIdx(tgt_idx).GetFormalCharge() != 0
  )
  b_segs = b["stroke_segments"]
  r_clear = max(r_tgt + 0.02, min(r_tgt + 0.035, 0.38 * bond_len))

  for sign in signs:
    n = sign * n_base
    outer_off = bond_outer_offset_cm(b, n)
    n0 = outer_off + min(0.105, 0.15 * bond_len)
    for f_start in f_start_cands:
      p0 = p_c + f_start * v + n0 * n
      for r_back in (0.35, 0.18, 0.0, 0.58, 0.80):
        rim_n_frac = math.sqrt(max(0.0, 1.01 - min(r_back, 0.98) ** 2))
        n3 = max(rim_n_frac * r_clear, outer_off + 0.075)
        p3 = p_tgt - (r_back * r_clear) * t + n3 * n
        base_arch = max(n3 + 0.075, n0 + 0.085, 0.42 * bond_len)
        f1 = f_start + 0.24 * (1.0 - f_start)
        f2 = f_start + 0.58 * (1.0 - f_start)
        f3 = 1.0 - (r_back * r_clear) / bond_len
        u_end = (p3 - p_tgt) / max(float(np.linalg.norm(p3 - p_tgt)), 1e-6)
        u_p2 = 0.60 * u_end + 0.40 * n
        u_p2 = u_p2 / max(float(np.linalg.norm(u_p2)), 1e-6)
        for arch_scale in (1.0, 1.25, 1.55):
          h_arch = arch_scale * base_arch
          p1 = p_c + f1 * v + h_arch * n
          if f2 <= f3 - 0.14:
            p2 = p_c + f2 * v + (h_arch + 0.01) * n
          else:
            p2 = p3 + (0.22 * arch_scale) * u_p2

          poly = sample_bezier_pts(p0, p1, p2, p3, n=30)
          wing_segs = arrowhead_wing_segments(
              p0, p1, p2, p3, electrons=electrons
          )

          min_d_bond = float("inf")
          for pt in poly[3:]:
            for s1, s2 in b_segs:
              min_d_bond = min(min_d_bond, dist_pt_to_seg(pt, s1, s2))
          cross = 20 if min_d_bond < 0.065 else 0
          cross += _scene_obstacle_penalty(
              poly,
              wing_segs,
              scene_canvases,
              owner_mol=mol,
              exclude_bonds=((a1_idx, a2_idx),),
              endpoint_atom_idx=tgt_idx,
          )
          arr_cross, arr_soft_pen = _prior_arrow_penalty(
              poly, wing_segs, prior_arrows, ctrl_pts=(p0, p1, p2, p3)
          )
          cross += arr_cross

          score = (
              -100.0 * cross
              - arr_soft_pen
              + min(min_d_bond, 0.20) * 20.0
              - 5.0 * abs(f_start - pref_f_start)
              - 1.8 * abs(r_back - 0.35)
              - 2.2 * (arch_scale - 1.0)
          )
          # Left-side HN+/HO+ labels place the 'H' prefix toward -x; prefer +x
          # normal side when `side` is unspecified to avoid crossing 'H'.
          if side is None:
            if tgt_has_h_label and p_tgt[0] < com_x:
              score += 18.0 * float(n[0])
            if attack_dir is not None and np.linalg.norm(attack_dir) > 1e-6:
              u_att = attack_dir / float(np.linalg.norm(attack_dir))
              score += 28.0 if float(np.dot(n, u_att)) < 0.0 else -28.0
            elif attack_pos is not None:
              v_att = attack_pos - p_c
              if np.dot(n, v_att) < 0:
                score += 15.0
            elif incoming_from_right:
              if abs(t[0]) < 0.35:
                score += 10.0 if n[0] < 0 else -10.0
              else:
                score += -n[0] + 0.6 * n[1]
            else:
              score += n[0] + 0.6 * n[1]

          cand.append((score, (p0, p1, p2, p3)))

  cand.sort(key=lambda x: x[0], reverse=True)
  return cand[0][1]


def compute_bond_to_bond_bezier(
    mol,
    b1_pair,
    b2_pair,
    *,
    canvases=None,
    allow_incipient: bool = False,
    tgt_mol=None,
    all_bonds=None,
    plus_x=None,
    side=None,
    lift=None,
    prior_arrows=None,
    electrons=2,
):
  """Computes a cubic Bézier shifting electrons from `b1_pair` to `b2_pair`.

  Args:
      mol: `MoleculeCanvas` owning the source bond.
      b1_pair: Source `(a1, a2)` bond atom tokens or indices.
      b2_pair: Target `(a1, a2)` bond atom tokens or indices.
      canvases: Optional scene `MoleculeCanvas` sequence for obstacle checks.
      allow_incipient: When True and `b2_pair` is an incipient (not yet formed)
        bond, routes from `b1_pair` to the incipient target atom.
      tgt_mol: Optional `MoleculeCanvas` owning the target bond when cross-mol.
      all_bonds: Optional prebuilt bond dicts forwarded for incipient routing.
      plus_x: Optional `+` glyph x-coordinates forwarded for incipient routing.
      side: Optional trajectory side hint forwarded for incipient routing.
      lift: Optional lift scale forwarded for incipient routing.
      prior_arrows: Optional list of already-routed arrow polylines to avoid.
      electrons: Electron count (`2` for pair shift, `1` for radical shift).

  Returns:
      Tuple `(p0, p1, p2, p3)` of 2D control points in cm, or `None` if source
      or target bond cannot be resolved.
  """

  def _resolve_b2b_ep(pref_canvas, raw_tok):
    try:
      return pref_canvas, resolve_atom_token(pref_canvas, raw_tok)
    except (ValueError, KeyError, TypeError, IndexError):
      if canvases and isinstance(raw_tok, str) and raw_tok.startswith(":"):
        for cand_c in canvases:
          if cand_c is pref_canvas or cand_c.is_product:
            continue
          try:
            resolved_i = resolve_atom_token(cand_c, raw_tok)
          except (ValueError, KeyError, TypeError, IndexError):
            continue
          if resolved_i in cand_c.atoms:
            return cand_c, resolved_i
      raise

  c1_a, i1_a = _resolve_b2b_ep(mol, b1_pair[0])
  c1_b, i1_b = _resolve_b2b_ep(mol, b1_pair[1])
  tgt_pref = tgt_mol if tgt_mol is not None else c1_a
  c2_a, i2_a = _resolve_b2b_ep(
      c1_a if b2_pair[0] in b1_pair else tgt_pref, b2_pair[0]
  )
  c2_b, i2_b = _resolve_b2b_ep(
      c1_a if b2_pair[1] in b1_pair else tgt_pref, b2_pair[1]
  )
  if c1_a is not c1_b:
    return None
  mol = c1_a
  b1_pair = (i1_a, i1_b)
  b1 = mol.find_bond(*b1_pair)
  if c2_a is not mol or c2_b is not mol:
    if not b1:
      return None
    if c2_a is not mol and c2_b is mol:
      other_c, other_i = c2_a, i2_a
    elif c2_b is not mol and c2_a is mol:
      other_c, other_i = c2_b, i2_b
    elif c2_a is c2_b:
      other_b = c2_a.find_bond(i2_a, i2_b)
      other_c = c2_a
      other_i = [i2_a, i2_b] if other_b else i2_a
    else:
      return None
    return compute_attack_bezier(
        mol,
        list(b1_pair),
        other_c,
        other_i,
        all_bonds=all_bonds,
        plus_x=plus_x,
        side=side,
        lift=lift,
        canvases=canvases,
        prior_arrows=prior_arrows,
        electrons=electrons,
    )
  b2_pair = (i2_a, i2_b)
  b2 = mol.find_bond(*b2_pair)
  if not b1 or not b2:
    if allow_incipient and b1 and not b2:
      shared_inc = list(set(b1_pair).intersection(set(b2_pair)))
      if shared_inc:
        tgt_inc = b2_pair[0] if b2_pair[1] == shared_inc[0] else b2_pair[1]
      else:
        tgt_inc = min(
            b2_pair,
            key=lambda idx: float(
                np.linalg.norm(mol.atoms[idx]["pos"] - b1["mid"])
            ),
        )
      return compute_attack_bezier(
          mol,
          list(b1_pair),
          mol,
          tgt_inc,
          all_bonds=all_bonds,
          plus_x=plus_x,
          side=side,
          lift=lift,
          canvases=canvases,
          prior_arrows=prior_arrows,
          electrons=electrons,
      )
    return None

  scene_canvases = list(canvases) if canvases is not None else [mol]
  shared = list(set(b1_pair).intersection(set(b2_pair)))
  if not shared:
    cand = []
    for sign1, sign2 in [(1.0, 1.0), (-1.0, -1.0), (1.0, -1.0), (-1.0, 1.0)]:
      n1 = sign1 * b1["norm"]
      n2 = sign2 * b2["norm"]
      off1 = bond_outer_offset_cm(b1, n1)
      off2 = bond_outer_offset_cm(b2, n2)
      n_lift_vec = n1 + n2
      n_lift_norm = float(np.linalg.norm(n_lift_vec))
      n_lift = n_lift_vec / n_lift_norm if n_lift_norm > 1e-4 else n1
      for f0 in [0.14, 0.10, 0.20]:
        for f3 in [0.16, 0.12, 0.22]:
          for lift_norm in [0.25, 0.18, 0.35, 0.55, 0.80, 1.15]:
            p0 = b1["mid"] + (f0 + off1) * n1
            p3 = b2["mid"] + (f3 + off2) * n2
            mid_ctrl = (p0 + p3) / 2.0 + lift_norm * n_lift
            p1 = p0 + 0.5 * (mid_ctrl - p0)
            p2 = p3 + 0.5 * (mid_ctrl - p3)
            poly = sample_bezier_pts(p0, p1, p2, p3, n=25)
            wing_segs = arrowhead_wing_segments(
                p0, p1, p2, p3, electrons=electrons
            )
            cross = _scene_obstacle_penalty(
                poly,
                wing_segs,
                scene_canvases,
                owner_mol=mol,
                exclude_bonds=(b1_pair, b2_pair),
                owner_glyph_pad=0.055,
                glyph_weight=2,
                bond_weight=1,
                end_bond_cm=0.065,
                end_bond_weight=2,
                trim_frac=0.10,
                wing_cm=0.028,
            )
            arr_cross, arr_soft_pen = _prior_arrow_penalty(
                poly, wing_segs, prior_arrows, ctrl_pts=(p0, p1, p2, p3)
            )
            cross += arr_cross
            cand.append((
                cross,
                arr_soft_pen + lift_norm,
                -(sign1 + sign2),
                (p0, p1, p2, p3),
            ))
    cand.sort(key=lambda x: (x[0], x[1], x[2]))
    return cand[0][3]

  c_idx = shared[0]
  c_pos = mol.atoms[c_idx]["pos"]

  a_idx = b1_pair[0] if b1_pair[1] == c_idx else b1_pair[1]
  b_idx = b2_pair[0] if b2_pair[1] == c_idx else b2_pair[1]
  a_pos = mol.atoms[a_idx]["pos"]
  b_pos = mol.atoms[b_idx]["pos"]

  v_a = a_pos - c_pos
  v_b = b_pos - c_pos
  l_a = np.linalg.norm(v_a)
  l_b = np.linalg.norm(v_b)
  if l_a == 0 or l_b == 0:
    return None
  u_a = v_a / l_a
  u_b = v_b / l_b

  ang_a = math.atan2(u_a[1], u_a[0])
  ang_b = math.atan2(u_b[1], u_b[0])

  other_nbr_angs = []
  for nbr in mol.mol.GetAtomWithIdx(c_idx).GetNeighbors():
    n_i = nbr.GetIdx()
    if n_i not in (a_idx, b_idx) and n_i in mol.atoms:
      v_n = mol.atoms[n_i]["pos"] - c_pos
      if np.linalg.norm(v_n) > 0:
        other_nbr_angs.append(math.atan2(v_n[1], v_n[0]))

  diff_ccw = (ang_b - ang_a) % (2 * math.pi)
  diff_cw = (ang_a - ang_b) % (2 * math.pi)

  def count_in_sector(start_ang, span_ang, sign):
    cnt = 0
    for oa in other_nbr_angs:
      rel = ((oa - start_ang) * sign) % (2 * math.pi)
      if 0.05 < rel < span_ang - 0.05:
        cnt += 1
    return cnt

  cnt_ccw = count_in_sector(ang_a, diff_ccw, 1.0)
  cnt_cw = count_in_sector(ang_a, diff_cw, -1.0)

  if cnt_ccw == 0 and cnt_cw > 0:
    chosen_dir = 1.0
  elif cnt_cw == 0 and cnt_ccw > 0:
    chosen_dir = -1.0
  else:
    chosen_dir = 1.0 if diff_ccw <= diff_cw else -1.0

  b1_obj = mol.mol.GetBondBetweenAtoms(b1_pair[0], b1_pair[1])
  b2_obj = mol.mol.GetBondBetweenAtoms(b2_pair[0], b2_pair[1])
  both_in_ring = (b1_obj is not None and b1_obj.IsInRing()) and (
      b2_obj is not None and b2_obj.IsInRing()
  )
  ring_center = None
  if both_in_ring:
    ring_info = mol.mol.GetRingInfo()
    for r in ring_info.AtomRings():
      if c_idx in r and a_idx in r and b_idx in r:
        ring_center = np.mean([mol.atoms[a]["pos"] for a in r], axis=0)
        break

  cand = []
  for dir_try, dir_pen in [(chosen_dir, 0.0), (-chosen_dir, 2.5)]:
    span_try = diff_ccw if dir_try > 0 else diff_cw
    n_a = (
        np.array([-u_a[1], u_a[0]])
        if dir_try > 0
        else np.array([u_a[1], -u_a[0]])
    )
    n_b = (
        np.array([u_b[1], -u_b[0]])
        if dir_try > 0
        else np.array([-u_b[1], u_b[0]])
    )
    if both_in_ring and ring_center is not None and dir_pen == 0.0:
      v_rc = ring_center - c_pos
      if np.linalg.norm(v_rc) > 1e-6:
        u_rc = v_rc / np.linalg.norm(v_rc)
        if np.dot(n_a, u_rc) < 0:
          n_a = -n_a
        if np.dot(n_b, u_rc) < 0:
          n_b = -n_b

    outer_a = bond_outer_offset_cm(b1, n_a)
    outer_b = bond_outer_offset_cm(b2, n_b)
    for extra_off in (0.0, 0.045):
      off_a = (outer_a + 0.085 if outer_a > 0 else 0.10) + extra_off
      off_b = (outer_b + 0.085 if outer_b > 0 else 0.10) + extra_off

      for f_start in (0.45, 0.55, 0.62):
        for f_end in (0.40, 0.50, 0.58):
          p0 = c_pos + f_start * v_a + off_a * n_a
          p3 = c_pos + f_end * v_b + off_b * n_b

          r_start = float(np.linalg.norm(p0 - c_pos))
          r_end = float(np.linalg.norm(p3 - c_pos))
          if r_start < 1e-6 or r_end < 1e-6:
            continue
          r_avg = (r_start + r_end) / 2.0

          u_rad_0 = (p0 - c_pos) / r_start
          t_dep = dir_try * np.array([-u_rad_0[1], u_rad_0[0]])

          u_rad_3 = (p3 - c_pos) / r_end
          t_arr = dir_try * np.array([-u_rad_3[1], u_rad_3[0]])

          theta = min(span_try, math.pi * 0.85)
          k_base = (4.0 / 3.0) * math.tan(theta / 4.0) * r_avg

          for k_mult in (0.85, 1.0, 1.25, 1.55):
            k = max(0.12, min(k_base * k_mult, 0.55))
            p1 = p0 + k * t_dep
            p2 = p3 - k * t_arr

            poly = sample_bezier_pts(p0, p1, p2, p3, n=30)
            wing_segs = arrowhead_wing_segments(
                p0, p1, p2, p3, electrons=electrons
            )
            cross = _scene_obstacle_penalty(
                poly,
                wing_segs,
                scene_canvases,
                owner_mol=mol,
                exclude_bonds=(b1_pair, b2_pair),
                owner_glyph_pad=0.0,
                other_glyph_pad=0.0,
                glyph_weight=3,
                bond_weight=4,
                wing_cm=0.025,
            )
            min_b_clr = 0.052
            for b_obj in (b1, b2):
              hit_b = False
              for bp1, bp2 in b_obj["stroke_segments"]:
                for s in range(len(poly) - 1):
                  if (
                      intersect_segments(poly[s], poly[s + 1], bp1, bp2)
                      or dist_pt_to_seg(poly[s], bp1, bp2) < min_b_clr
                      or dist_pt_to_seg(poly[s + 1], bp1, bp2) < min_b_clr
                  ):
                    hit_b = True
                    break
                if hit_b:
                  cross += 4
                  break

            arr_cross, arr_soft_pen = _prior_arrow_penalty(
                poly, wing_segs, prior_arrows, ctrl_pts=(p0, p1, p2, p3)
            )
            cross += arr_cross
            d_c = float(np.min(np.linalg.norm(poly - c_pos, axis=1)))
            cand.append((
                cross,
                dir_pen
                + arr_soft_pen
                + 15.0 * extra_off
                + 0.8 * abs(k_mult - 1.0)
                + abs(d_c - 0.35),
                (p0, p1, p2, p3),
            ))

  cand.sort(key=lambda x: (x[0], x[1]))
  return cand[0][2]


def compute_lone_pair_to_bond_bezier(
    mol,
    atom_idx,
    bond_pair,
    side=None,
    canvases=None,
    allow_incipient: bool = False,
    prior_arrows=None,
    electrons=2,
    all_bonds=None,
    plus_x=None,
):
  """Computes a lone-pair-to-bond cubic Bézier from `atom_idx` onto `bond_pair`.

  Args:
      mol: `MoleculeCanvas` owning both `atom_idx` and `bond_pair`.
      atom_idx: Donor atom token or index.
      bond_pair: Target `(a1, a2)` bond atom tokens or indices.
      side: Optional normal-side override (`"up"`, `"down"`, `"left"`,
        `"right"`).
      canvases: Optional scene `MoleculeCanvas` sequence for obstacle checks.
      allow_incipient: When True and `bond_pair` does not yet exist on `mol`,
        falls back to an incipient atom-to-atom trajectory.
      prior_arrows: Optional list of already-routed arrow polylines to avoid.
      electrons: Electron count (`2` for lone pair, `1` for unpaired radical).
      all_bonds: Optional prebuilt bond dicts forwarded to incipient attacks.
      plus_x: Optional `+` separator x-coordinate(s) forwarded to incipient
        intermolecular attacks.

  Returns:
      Tuple `(p0, p1, p2, p3)` of 2D control points in cm, or `None` if
      `atom_idx` or `bond_pair` does not exist on `mol`.
  """
  atom_idx = resolve_atom_token(mol, atom_idx)
  try:
    bond_pair = (
        resolve_atom_token(mol, bond_pair[0]),
        resolve_atom_token(mol, bond_pair[1]),
    )
  except (ValueError, KeyError, TypeError, IndexError):
    if allow_incipient and canvases and atom_idx in mol.atoms:
      for raw_ep in bond_pair:
        if not (isinstance(raw_ep, str) and raw_ep.startswith(":")):
          continue
        try:
          if resolve_atom_token(mol, raw_ep) == atom_idx:
            continue
        except (ValueError, KeyError, TypeError, IndexError):
          pass
        for cand_canvas in canvases:
          if cand_canvas is mol or cand_canvas.is_product:
            continue
          try:
            other_idx = resolve_atom_token(cand_canvas, raw_ep)
          except (ValueError, KeyError, TypeError, IndexError):
            continue
          if other_idx in cand_canvas.atoms:
            return compute_attack_bezier(
                mol,
                atom_idx,
                cand_canvas,
                other_idx,
                all_bonds=all_bonds,
                plus_x=plus_x,
                side=side,
                canvases=canvases,
                prior_arrows=prior_arrows,
                electrons=electrons,
            )
    raise
  b = mol.find_bond(*bond_pair)
  if not b or atom_idx not in mol.atoms:
    if allow_incipient and not b and atom_idx in mol.atoms:
      if atom_idx in bond_pair:
        tgt_inc = bond_pair[0] if bond_pair[1] == atom_idx else bond_pair[1]
      else:
        tgt_inc = min(
            bond_pair,
            key=lambda idx: float(
                np.linalg.norm(
                    mol.atoms[idx]["pos"] - mol.atoms[atom_idx]["pos"]
                )
            ),
        )
      return compute_attack_bezier(
          mol,
          atom_idx,
          mol,
          tgt_inc,
          all_bonds=all_bonds,
          side=side,
          canvases=canvases,
          prior_arrows=prior_arrows,
          electrons=electrons,
      )
    return None

  p_atom = mol.atoms[atom_idx]["pos"]
  r_atom = mol.atoms[atom_idx]["glyph_r"]
  other_idx = (
      bond_pair[1]
      if bond_pair[0] == atom_idx
      else (bond_pair[0] if bond_pair[1] == atom_idx else None)
  )
  p_other = mol.atoms[other_idx]["pos"] if other_idx in mol.atoms else b["mid"]

  v_ab = p_other - p_atom
  l_ab = float(np.linalg.norm(v_ab))
  if other_idx is not None:
    u_along = v_ab / l_ab if l_ab > 1e-6 else b["tangent"]
    n_base = np.array([-u_along[1], u_along[0]])
  else:
    u_along = (
        b["tangent"]
        if float(np.dot(v_ab, b["tangent"])) >= 0.0
        else -b["tangent"]
    )
    n_base = b["norm"]

  signs = _signs_for_side(n_base, side)
  n_ring_ext = _ring_exterior_normal(mol, bond_pair[0], bond_pair[1], n_base)
  sub_vecs: list[np.ndarray] = []
  for end_i in (
      (atom_idx, other_idx) if other_idx is not None else (atom_idx,)
  ):
    if end_i is None or end_i not in mol.atoms:
      continue
    p_end = mol.atoms[end_i]["pos"]
    for nbr in mol.mol.GetAtomWithIdx(end_i).GetNeighbors():
      nb_i = nbr.GetIdx()
      if nb_i not in bond_pair and nb_i in mol.atoms:
        v_nb = mol.atoms[nb_i]["pos"] - p_end
        if np.linalg.norm(v_nb) > 1e-6:
          sub_vecs.append(v_nb / float(np.linalg.norm(v_nb)))

  scene_canvases = list(canvases) if canvases is not None else [mol]
  b_segs = b["stroke_segments"]
  cand = []

  for sign in signs:
    n = sign * n_base
    outer_off = bond_outer_offset_cm(b, n)
    n3_off = 0.14 if outer_off == 0 else outer_off + 0.095

    hemi_bonus = 0.0
    if n_ring_ext is not None:
      hemi_bonus += 18.0 if float(np.dot(n, n_ring_ext)) > 0.0 else -18.0
    for u_sub in sub_vecs:
      dot_s = float(np.dot(n, u_sub))
      if dot_s > 0.20:
        hemi_bonus -= 12.0 * dot_s
      elif dot_s < -0.20:
        hemi_bonus += 4.0 * (-dot_s)

    for f_land in (0.50, 0.43, 0.57):
      p_land_base = (
          p_atom + f_land * v_ab if other_idx is not None else b["mid"]
      )
      p3 = p_land_base + n3_off * n
      for dep_tilt in (0.0, -0.22, 0.22):
        dep_dir = n + dep_tilt * (-u_along)
        dep_dir = dep_dir / max(float(np.linalg.norm(dep_dir)), 1e-6)
        p0 = p_atom + (r_atom + 0.055) * dep_dir
        for arch_h in (0.20, 0.28, 0.38):
          p1 = p0 + arch_h * n + 0.09 * u_along
          p2 = p3 + (0.90 * arch_h) * n - 0.07 * u_along
          poly = sample_bezier_pts(p0, p1, p2, p3, n=30)
          wing_segs = arrowhead_wing_segments(
              p0, p1, p2, p3, electrons=electrons
          )

          min_d_tgt = float("inf")
          for pt in poly[3:-3]:
            for s1, s2 in b_segs:
              min_d_tgt = min(min_d_tgt, dist_pt_to_seg(pt, s1, s2))
          cross = 12 if min_d_tgt < 0.065 else 0
          cross += _scene_obstacle_penalty(
              poly,
              wing_segs,
              scene_canvases,
              owner_mol=mol,
              exclude_bonds=(bond_pair,),
              endpoint_atom_idx=atom_idx,
              min_bond_cm=0.085,
              wing_cm=0.030,
          )
          arr_cross, arr_soft_pen = _prior_arrow_penalty(
              poly, wing_segs, prior_arrows, ctrl_pts=(p0, p1, p2, p3)
          )
          cross += arr_cross

          score = (
              -100.0 * cross
              - arr_soft_pen
              + hemi_bonus
              + 15.0 * min(min_d_tgt, 0.18)
              - 4.0 * abs(f_land - 0.50)
              - 2.0 * abs(dep_tilt)
              - 3.0 * (arch_h - 0.20)
          )
          cand.append((score, (p0, p1, p2, p3)))

  cand.sort(key=lambda x: x[0], reverse=True)
  return cand[0][1]


def compute_bond_to_atom_bezier(
    src_mol,
    bond_pair,
    tgt_mol,
    tgt_atom,
    *,
    lift=None,
    canvases=None,
    prior_arrows=None,
    electrons=2,
    plus_x=None,
):
  """Computes a cubic Bézier from bond `bond_pair` onto atom `tgt_atom`.

  Args:
      src_mol: `MoleculeCanvas` owning `bond_pair`.
      bond_pair: Source `(a1, a2)` bond atom tokens or indices.
      tgt_mol: `MoleculeCanvas` owning `tgt_atom`.
      tgt_atom: Target atom token or index.
      lift: Optional vertical control-handle lift in cm.
      canvases: Optional scene `MoleculeCanvas` sequence for obstacle checks.
      prior_arrows: Optional list of already-routed arrow polylines to avoid.
      electrons: Electron count (`2` for pair, `1` for radical).
      plus_x: Optional x-coordinate(s) of inter-molecule `+` glyphs to clear.

  Returns:
      Tuple `(p0, p1, p2, p3)` of 2D control points in cm, or `None` if the
      source bond or target atom cannot be resolved.
  """
  a1 = resolve_atom_token(src_mol, bond_pair[0])
  a2 = resolve_atom_token(src_mol, bond_pair[1])
  b = src_mol.find_bond(a1, a2)
  tgt_a_idx = resolve_atom_token(tgt_mol, tgt_atom)
  tgt_a = tgt_mol.atoms.get(tgt_a_idx)
  if not b or not tgt_a:
    return None

  scene_canvases = (
      list(canvases)
      if canvases is not None
      else [m for m in dict.fromkeys((src_mol, tgt_mol)) if m is not None]
  )
  p_list = (
      (plus_x if isinstance(plus_x, (list, tuple)) else [plus_x])
      if plus_x is not None
      else ()
  )
  lift_cands = (
      [float(lift)] if lift is not None else [0.8, -0.8, 1.1, -1.1, 1.45, -1.45]
  )
  best_bta = None
  best_bta_sc = -1e9
  for ctrl_lift in lift_cands:
    n = b["norm"]
    if abs(n[1]) < 0.25:
      if float(np.dot(tgt_a["pos"] - b["mid"], n)) < 0.0:
        n = -n
    elif (ctrl_lift > 0 and n[1] < 0) or (ctrl_lift < 0 and n[1] > 0):
      n = -n
    outer_off = bond_outer_offset_cm(b, n)
    p0_c = b["mid"] + (0.14 if outer_off == 0 else outer_off + 0.095) * n
    dx_c = tgt_a["pos"][0] - p0_c[0]
    offset_x = (
        -(tgt_a["glyph_r"] + 0.06) if dx_c > 0 else (tgt_a["glyph_r"] + 0.06)
    )
    p3_c = tgt_a["pos"] + np.array([offset_x, 0.0])
    dist_c = float(np.linalg.norm(p3_c - p0_c))
    step_x = dist_c * 0.25 * (1 if dx_c > 0 else -1)
    p1_c = p0_c + np.array([step_x, ctrl_lift])
    p2_c = p3_c + np.array([-step_x, ctrl_lift * 0.8])
    poly_bta = sample_bezier_pts(p0_c, p1_c, p2_c, p3_c, n=35)
    wings_bta = arrowhead_wing_segments(
        p0_c, p1_c, p2_c, p3_c, electrons=electrons
    )
    bta_cross = _scene_obstacle_penalty(
        poly_bta,
        wings_bta,
        scene_canvases,
        owner_mol=src_mol,
        exclude_bonds=((a1, a2),),
        glyph_weight=4,
        bond_weight=3,
        wing_cm=0.030,
    )
    if p_list and hits_plus_sign(poly_bta, p_list, wings_bta):
      bta_cross += 6
    arr_cr, bta_soft = _prior_arrow_penalty(
        poly_bta, wings_bta, prior_arrows, ctrl_pts=(p0_c, p1_c, p2_c, p3_c)
    )
    bta_cross += arr_cr
    bta_sc = (
        -1000.0 * bta_cross
        - bta_soft
        - (4.0 if ctrl_lift < 0 else 0.0)
        - 2.0 * abs(abs(ctrl_lift) - 0.8)
    )
    if bta_sc > best_bta_sc:
      best_bta_sc = bta_sc
      best_bta = (p0_c, p1_c, p2_c, p3_c)
  return best_bta


def _extract_molecules(rxn_def, key_prefix: str):
  """Extracts molecule definitions from rxn_def for 'r' (reactants) or 'p' (products).

  Returns:
      List of (key, config_dict).
  """
  plural_key = "reactants" if key_prefix == "r" else "products"
  if plural_key in rxn_def:
    raw = rxn_def[plural_key]
    if isinstance(raw, list):
      items = []
      for idx, item in enumerate(raw):
        if item is not None:
          k = item.get("key", f"{key_prefix}{idx + 1}")
          items.append((k, dict(item)))
      return items
    elif isinstance(raw, dict):
      return [(k, dict(v)) for k, v in raw.items() if v is not None]
    else:
      raise ValueError(
          f"Invalid {plural_key} format: expected list or dict, got {type(raw)}"
      )

  # Search for key_prefix + digits, e.g. r1, r2, ... or p1, p2, ...
  keys = [
      k
      for k in rxn_def.keys()
      if re.match(rf"^{key_prefix}\d+$", k) and rxn_def[k] is not None
  ]
  keys.sort(key=lambda k: int(k[len(key_prefix) :]))
  return [(k, dict(rxn_def[k])) for k in keys]


def _canon_smi_no_iso(s: str) -> str:
  m = _parse_smiles_preserving_mapped_hs(s)
  if m is None:
    return s
  for atom in m.GetAtoms():
    atom.SetAtomMapNum(0)
  try:
    m = Chem.RemoveHs(m)
  except Exception:
    pass
  return Chem.MolToSmiles(m, isomericSmiles=False)


def build_graph_and_arrows_from_spec(rxn_def, spec):
  """Builds a MolecularGraph and list of CurvedArrow objects from rxn_def and spec."""
  rid = rxn_def.get("id", "reaction")
  add_h_spec = spec.get("add_h", {})
  reactant_items = _extract_molecules(rxn_def, "r")
  mg = MolecularGraph(name=rid)
  r_mols = {}
  r_map_to_idx = {}

  for k, r_cfg in reactant_items:
    smi = r_cfg["smi"]
    m = _parse_smiles_preserving_mapped_hs(smi)
    if m is None:
      raise ValueError(f"{rid}: Invalid reactant SMILES for '{k}': {smi}")
    map_to_idx = {
        a.GetAtomMapNum(): a.GetIdx()
        for a in m.GetAtoms()
        if a.GetAtomMapNum() > 0
    }

    if k in add_h_spec:
      add_h_val = add_h_spec[k]
      add_h_list = (
          [add_h_val] if isinstance(add_h_val, (int, str)) else list(add_h_val)
      )
      rw = Chem.RWMol(m)
      for atom_tok in add_h_list:
        atom_idx = _resolve_mol_atom_token(m, atom_tok, map_to_idx)
        ta = rw.GetAtomWithIdx(int(atom_idx))
        num_h = ta.GetTotalNumHs()
        if num_h > 0:
          ta.SetNumExplicitHs(max(0, num_h - 1))
          ta.SetNoImplicit(True)
        h_idx = rw.AddAtom(Chem.Atom(1))
        rw.AddBond(int(atom_idx), h_idx, Chem.BondType.SINGLE)
      m = rw.GetMol()
      Chem.SanitizeMol(m)

    r_mols[k] = m
    r_map_to_idx[k] = map_to_idx

    arom_bonds = {
        (
            min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
        )
        for b in m.GetBonds()
        if b.GetIsAromatic()
    }
    selected_kekule = None
    if arom_bonds:
      net_arom_delta = {}

      def _resolve_edge(u_tok, v_tok, mol_obj=m, map_dict=map_to_idx):
        ui = int(_resolve_mol_atom_token(mol_obj, u_tok, map_dict))
        vi = int(_resolve_mol_atom_token(mol_obj, v_tok, map_dict))
        return (min(ui, vi), max(ui, vi))

      for a_spec in spec.get("arrows", []):
        atype = a_spec.get("type")
        try:
          if atype == "attack":
            src = a_spec.get("src")
            if (
                isinstance(src, (list, tuple))
                and len(src) == 2
                and src[0] == k
                and isinstance(src[1], (list, tuple))
            ):
              e = _resolve_edge(src[1][0], src[1][1])
              if e in arom_bonds:
                net_arom_delta[e] = net_arom_delta.get(e, 0) - 1
          elif atype == "bond_cleave" and a_spec.get("mol", "r1") == k:
            u_t, v_t = a_spec["bond"]
            e = _resolve_edge(u_t, v_t)
            if e in arom_bonds:
              net_arom_delta[e] = net_arom_delta.get(e, 0) - 1
          elif atype == "lone_pair_to_bond" and a_spec.get("mol", "r1") == k:
            u_t, v_t = a_spec["bond"]
            e = _resolve_edge(u_t, v_t)
            if e in arom_bonds:
              net_arom_delta[e] = net_arom_delta.get(e, 0) + 1
          elif atype == "bond_to_bond":
            b2b_src_k = a_spec.get("src_mol") or a_spec.get("mol", "r1")
            b2b_tgt_k = a_spec.get("tgt_mol") or a_spec.get("mol", b2b_src_k)
            b1 = a_spec.get("b1") or a_spec.get("src_bond")
            b2 = a_spec.get("b2") or a_spec.get("tgt_bond")
            if b1 and b2b_src_k == k:
              e1 = _resolve_edge(b1[0], b1[1])
              if e1 in arom_bonds:
                net_arom_delta[e1] = net_arom_delta.get(e1, 0) - 1
            if b2 and b2b_src_k == k and b2b_tgt_k == k:
              e2 = _resolve_edge(b2[0], b2[1])
              if e2 in arom_bonds:
                net_arom_delta[e2] = net_arom_delta.get(e2, 0) + 1
          elif atype == "bond_to_atom" and a_spec.get("src_mol", "r1") == k:
            u_t, v_t = a_spec["bond"]
            e = _resolve_edge(u_t, v_t)
            if e in arom_bonds:
              net_arom_delta[e] = net_arom_delta.get(e, 0) - 1
        except (KeyError, IndexError, TypeError, ValueError):
          pass

      if any(v != 0 for v in net_arom_delta.values()):
        try:
          base_q = [a.GetFormalCharge() for a in m.GetAtoms()]
          non_arom_bt = {
              b.GetIdx(): b.GetBondType()
              for b in m.GetBonds()
              if not b.GetIsAromatic()
          }
          suppl = Chem.ResonanceMolSupplier(m, Chem.KEKULE_ALL)
          for cand in suppl:
            if cand is None:
              continue
            if [a.GetFormalCharge() for a in cand.GetAtoms()] != base_q:
              continue
            if any(
                cand.GetBondWithIdx(i).GetBondType() != bt
                for i, bt in non_arom_bt.items()
            ):
              continue
            ok = True
            for edge_uv, delta_bo in net_arom_delta.items():
              b_cand = cand.GetBondBetweenAtoms(edge_uv[0], edge_uv[1])
              if b_cand is None:
                ok = False
                break
              if delta_bo < 0 and b_cand.GetBondType() != Chem.BondType.DOUBLE:
                ok = False
                break
              if delta_bo > 0 and b_cand.GetBondType() != Chem.BondType.SINGLE:
                ok = False
                break
            if ok:
              selected_kekule = Chem.Mol(cand)
              for b_c in selected_kekule.GetBonds():
                b_c.SetIsAromatic(False)
              for a_c in selected_kekule.GetAtoms():
                a_c.SetIsAromatic(False)
              break
        except (Chem.KekulizeException, RuntimeError, ValueError) as exc:
          raise ValueError(
              f"{rid}: Failed to enumerate Kekulé resonance contributors for"
              f" reactant {k!r} ({smi!r}): {exc}"
          ) from exc

    if selected_kekule is not None:
      m = selected_kekule
      r_mols[k] = m
    else:
      try:
        Chem.Kekulize(m, clearAromaticFlags=True)
      except (Chem.KekulizeException, ValueError) as exc:
        raise ValueError(
            f"{rid}: Failed to kekulize reactant {k!r} ({smi!r}): {exc}"
        ) from exc

    if not hasattr(mg, "r_kekule_map"):
      mg.r_kekule_map = {}
    mg.r_kekule_map[k] = {
        (
            min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
        ): b.GetBondType()
        for b in m.GetBonds()
    }

    for atom in m.GetAtoms():
      aid = f"{k}_{atom.GetIdx()}"
      elem = atom.GetSymbol()
      q = atom.GetFormalCharge()
      rad = 0 if elem in D_F_BLOCK_METALS else atom.GetNumRadicalElectrons()
      impl_h = atom.GetTotalNumHs()
      v_free = VALENCE_ELECTRONS.get(elem, 4) if elem != "*" else 0
      exp_val = 0
      for b in atom.GetBonds():
        if b.GetBondType() == Chem.BondType.DATIVE:
          if b.GetBeginAtomIdx() == atom.GetIdx():
            exp_val += 2
        else:
          exp_val += int(round(b.GetBondTypeAsDouble()))
      rem_e = v_free - q - rad - (exp_val + impl_h) if elem != "*" else 0
      lp = max(0, rem_e // 2) if elem != "*" else 0
      unpaired_d = (
          (rem_e % 2) if (elem in D_F_BLOCK_METALS and rem_e > 0) else 0
      )
      mg.add_atom(
          atom_id=aid,
          element=elem,
          lone_pairs=lp,
          radicals=rad,
          unpaired_d=unpaired_d,
          formal_charge=q,
          implicit_h=impl_h,
      )

    for bond in m.GetBonds():
      u = f"{k}_{bond.GetBeginAtomIdx()}"
      v = f"{k}_{bond.GetEndAtomIdx()}"
      is_dat = bond.GetBondType() == Chem.BondType.DATIVE
      bo = 1 if is_dat else int(round(bond.GetBondTypeAsDouble()))
      st = bond.GetStereo()
      stereo_str = "NONE"
      stereo_atoms = None
      if bo == 2 and st in (
          Chem.BondStereo.STEREOE,
          Chem.BondStereo.STEREOZ,
          Chem.BondStereo.STEREOCIS,
          Chem.BondStereo.STEREOTRANS,
      ):
        stereo_str = (
            "E"
            if st in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS)
            else "Z"
        )
        sa = list(bond.GetStereoAtoms())
        if len(sa) == 2:
          stereo_atoms = (f"{k}_{sa[0]}", f"{k}_{sa[1]}")
      mg.add_bond(
          u,
          v,
          order=max(1, bo),
          stereo=stereo_str,
          stereo_atoms=stereo_atoms,
          is_dative=is_dat,
          dative_begin=u if is_dat else None,
          dative_end=v if is_dat else None,
      )

  curved_arrows = []

  def to_aid(m_key, a_tok):
    if (
        isinstance(a_tok, (tuple, list))
        and len(a_tok) == 2
        and isinstance(a_tok[0], str)
        and a_tok[0] in r_mols
    ):
      m_key, a_tok = a_tok[0], a_tok[1]
    if isinstance(a_tok, str) and a_tok.startswith(":") and a_tok[1:].isdigit():
      want_map = int(a_tok[1:])
      if want_map not in (r_map_to_idx.get(m_key) or {}):
        alt_keys = [k2 for k2, m2 in r_map_to_idx.items() if want_map in m2]
        if len(alt_keys) == 1:
          m_key = alt_keys[0]
    resolved = _resolve_mol_atom_token(
        r_mols[m_key], a_tok, r_map_to_idx.get(m_key)
    )
    return f"{m_key}_{resolved}"

  for a in spec.get("arrows", []):
    atype = a.get("type")
    electrons = int(a.get("electrons", 2))
    stereo_policy = a.get("stereo_policy", "auto")
    if atype == "attack":
      src = a["src"]
      tgt = a["tgt"]
      tgt_mol, tgt_a = tgt
      tgt_id = to_aid(tgt_mol, tgt_a)
      if isinstance(src[1], (tuple, list)):
        src_mol = src[0]
        u, v = src[1]
        src_id = (to_aid(src_mol, u), to_aid(src_mol, v))
        conn_atom = a.get("conn", u)
        conn_id = to_aid(src_mol, conn_atom)
        curved_arrows.append(
            CurvedArrow(
                "bond",
                src_id,
                "bond",
                (conn_id, tgt_id),
                electrons=electrons,
                label=a.get("desc", ""),
                stereo_policy=stereo_policy,
            )
        )
      else:
        src_mol, src_a = src
        src_id = to_aid(src_mol, src_a)
        src_node = mg.atoms.get(src_id)
        src_type = (
            "radical"
            if (electrons == 1 and src_node and src_node.radicals > 0)
            else a.get("source_type", "lone_pair")
        )
        curved_arrows.append(
            CurvedArrow(
                src_type,
                src_id,
                "bond",
                (src_id, tgt_id),
                electrons=electrons,
                label=a.get("desc", ""),
                stereo_policy=stereo_policy,
            )
        )
    elif atype == "bond_cleave":
      mol = a.get("mol", "r1")
      u, v = a["bond"]
      tgt = a.get("tgt", v)
      src_id = (to_aid(mol, u), to_aid(mol, v))
      tgt_id = to_aid(mol, tgt)
      curved_arrows.append(
          CurvedArrow(
              "bond",
              src_id,
              "atom",
              tgt_id,
              electrons=electrons,
              label=a.get("desc", ""),
              stereo_policy=stereo_policy,
          )
      )
    elif atype == "lone_pair_to_bond":
      mol = a.get("mol", "r1")
      atom = a["atom"]
      u, v = a["bond"]
      src_id = to_aid(mol, atom)
      sink_id = (to_aid(mol, u), to_aid(mol, v))
      src_node = mg.atoms.get(src_id)
      src_type = (
          "radical"
          if (electrons == 1 and src_node and src_node.radicals > 0)
          else a.get("source_type", "lone_pair")
      )
      curved_arrows.append(
          CurvedArrow(
              src_type,
              src_id,
              "bond",
              sink_id,
              electrons=electrons,
              label=a.get("desc", ""),
              stereo_policy=stereo_policy,
          )
      )
    elif atype == "bond_to_bond":
      src_mol = a.get("src_mol") or a.get("mol", "r1")
      tgt_mol = a.get("tgt_mol") or a.get("mol", src_mol)
      u1, v1 = a.get("b1") or a["src_bond"]
      u2, v2 = a.get("b2") or a["tgt_bond"]
      src_id = (to_aid(src_mol, u1), to_aid(src_mol, v1))
      u2_mol = src_mol if (tgt_mol != src_mol and u2 in (u1, v1)) else tgt_mol
      v2_mol = src_mol if (tgt_mol != src_mol and v2 in (u1, v1)) else tgt_mol
      sink_id = (to_aid(u2_mol, u2), to_aid(v2_mol, v2))
      curved_arrows.append(
          CurvedArrow(
              "bond",
              src_id,
              "bond",
              sink_id,
              electrons=electrons,
              label=a.get("desc", ""),
              stereo_policy=stereo_policy,
          )
      )
    elif atype == "bond_to_atom":
      src_mol = a.get("src_mol", "r1")
      tgt_mol = a.get("tgt_mol", "r2")
      u, v = a["bond"]
      tgt_atom = a["tgt_atom"]
      src_id = (to_aid(src_mol, u), to_aid(src_mol, v))
      tgt_id = to_aid(tgt_mol, tgt_atom)
      curved_arrows.append(
          CurvedArrow(
              "bond",
              src_id,
              "atom",
              tgt_id,
              electrons=electrons,
              label=a.get("desc", ""),
              stereo_policy=stereo_policy,
          )
      )

  return mg, curved_arrows


_AA_PATTERNS = [
    ("HIS", r"(?i)\b(his\d+[a-z]?|his|histidine)\b"),
    ("ARG", r"(?i)\b(arg\d+[a-z]?|arg|arginine)\b"),
    ("LYS", r"(?i)\b(lys\d+[a-z]?|lys|lysine)\b"),
    ("ASP", r"(?i)\b(asp\d+[a-z]?|asp|aspartate|aspartic)\b"),
    ("GLU", r"(?i)\b(glu\d+[a-z]?|glu|glutamate|glutamic)\b"),
    ("ASN", r"(?i)\b(asn\d+[a-z]?|asn|asparagine)\b"),
    ("GLN", r"(?i)\b(gln\d+[a-z]?|gln|glutamine)\b"),
    ("CYS", r"(?i)\b(cys\d+[a-z]?|cys|cysteine)\b"),
    ("SER", r"(?i)\b(ser\d+[a-z]?|ser|serine)\b"),
    ("THR", r"(?i)\b(thr\d+[a-z]?|thr|threonine)\b"),
    ("TYR", r"(?i)\b(tyr\d+[a-z]?|tyr|tyrosine)\b"),
    ("TRP", r"(?i)\b(trp\d+[a-z]?|trp|tryptophan)\b"),
    ("PHE", r"(?i)\b(phe\d+[a-z]?|phe|phenylalanine)\b"),
    ("MET", r"(?i)\b(met\d+[a-z]?|met|methionine)\b"),
    ("PRO", r"(?i)\b(pro\d+[a-z]?|pro|proline)\b"),
]

_COFACTOR_ANCHOR_LABELS = frozenset({
    "P",
    "PO3",
    "PO4",
    "PP",
    "TRNA",
    "RIB",
    "AD",
    "AMP",
    "ADP",
    "ATP",
    "NAD",
    "NADP",
    "FAD",
    "FMN",
    "COA",
    "SAM",
    "PLP",
    "TPP",
    "BIOTIN",
    "HEME",
    "FE",
    "MO",
    "MG",
    "ZN",
    "MN",
    "CU",
    "CO",
    "NI",
    "CA",
    "W",
    "V",
})


def _validate_amino_acid_sidechain(rid: str, mk: str, mcfg: dict) -> None:
  """Validates that any molecule representing a canonical amino acid residue has its expected chemical side chain."""
  mname = str(mcfg.get("name", ""))
  msmi = str(mcfg.get("smi", ""))
  if not msmi:
    return

  atom_labels = mcfg.get("atom_labels") or {}
  lbl_values = [str(v).strip() for v in atom_labels.values() if str(v).strip()]
  atom_lbl_text = " ".join(lbl_values)

  matched_aas = []
  for aa_code, pat in _AA_PATTERNS:
    if re.search(pat, atom_lbl_text):
      matched_aas.append(aa_code)

  if not matched_aas:
    if lbl_values and all(
        re.sub(r"[^A-Za-z0-9]", "", v).upper() in _COFACTOR_ANCHOR_LABELS
        for v in lbl_values
    ):
      return
    for aa_code, pat in _AA_PATTERNS:
      if re.search(pat, mname):
        matched_aas.append(aa_code)

  matched_aas = list(dict.fromkeys(matched_aas))
  if len(matched_aas) != 1:
    return

  aa = matched_aas[0]
  mol = Chem.MolFromSmiles(msmi)
  if mol is None:
    return

  rings = mol.GetRingInfo().AtomRings()
  elems = [a.GetSymbol() for a in mol.GetAtoms() if a.GetSymbol() != "*"]

  def _has_c3n2_ring() -> bool:
    for r in rings:
      if len(r) == 5:
        re_syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in r]
        if re_syms.count("N") == 2 and re_syms.count("C") == 3:
          return True
    return False

  def _has_guanidinium() -> bool:
    for a in mol.GetAtoms():
      if a.GetSymbol() == "C":
        n_nbs = [nb for nb in a.GetNeighbors() if nb.GetSymbol() == "N"]
        if len(n_nbs) >= 3:
          return True
    return False

  def _has_carboxyl_or_acyl_o2() -> bool:
    for a in mol.GetAtoms():
      if a.GetSymbol() == "C":
        o_nbs = [nb for nb in a.GetNeighbors() if nb.GetSymbol() == "O"]
        if len(o_nbs) >= 2:
          return True
    return False

  def _has_amide_cn_o() -> bool:
    for a in mol.GetAtoms():
      if a.GetSymbol() == "C":
        nb_syms = {nb.GetSymbol() for nb in a.GetNeighbors()}
        if "O" in nb_syms and "N" in nb_syms:
          return True
    return False

  def _has_phenol_ring() -> bool:
    for r in rings:
      if len(r) == 6:
        re_syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in r]
        if re_syms.count("C") == 6:
          for i in r:
            for nb in mol.GetAtomWithIdx(i).GetNeighbors():
              if nb.GetIdx() not in r and nb.GetSymbol() == "O":
                return True
    return False

  if aa == "HIS":
    if not _has_c3n2_ring():
      raise ValueError(
          f"[Amino Acid Side-Chain Error] {rid} ({mk} {mname!r}): Molecule is"
          f" labeled as a Histidine (His) residue, but SMILES {msmi!r} does not"
          " contain a 5-membered C3N2 imidazole/imidazolium ring. Use"
          " '*c1c[nH+]c[nH]1' for protonated His-H+ (general acid) and"
          " '*c1cnc[nH]1' for neutral His (general base), not a primary"
          " alkylamine/alkylammonium ('*[NH3+]' / '*N') which represents"
          " Lysine."
      )
  elif aa == "ARG":
    if not _has_guanidinium():
      raise ValueError(
          f"[Amino Acid Side-Chain Error] {rid} ({mk} {mname!r}): Molecule is"
          f" labeled as an Arginine (Arg) residue, but SMILES {msmi!r} does not"
          " contain a guanidinium/guanidino group (C bonded to 3 N atoms). Use"
          " '*NC(=[NH2+])N' for protonated Arg-H+ and '*NC(=N)N' for neutral"
          " Arg."
      )
  elif aa == "LYS":
    if "N" not in elems or _has_c3n2_ring() or _has_guanidinium():
      raise ValueError(
          f"[Amino Acid Side-Chain Error] {rid} ({mk} {mname!r}): Molecule is"
          f" labeled as a Lysine (Lys) residue, but SMILES {msmi!r} does not"
          " contain an amine/nitrogen side chain ('N' required without"
          " imidazole/guanidinium). Use '*CC[NH3+]' or '*[NH3+]' for protonated"
          " Lys-H+ and '*CCN' or '*N' for neutral Lys."
      )
  elif aa in ("ASP", "GLU"):
    aa_full = "Aspartate (Asp)" if aa == "ASP" else "Glutamate (Glu)"
    if not _has_carboxyl_or_acyl_o2():
      raise ValueError(
          f"[Amino Acid Side-Chain Error] {rid} ({mk} {mname!r}): Molecule is"
          f" labeled as an {aa_full} residue, but SMILES {msmi!r} does not"
          " contain a carboxyl/carboxylate group (C bonded to two O atoms). Use"
          " '*CC(=O)[O-]' or '*C(=O)[O-]' for the carboxylate base and"
          " '*CC(=O)O' or '*C(=O)O' for the carboxylic acid, not a bare"
          " alkoxide/hydroxyl ('*[O-]' / '*O') or amine."
      )
  elif aa in ("ASN", "GLN"):
    aa_full = "Asparagine (Asn)" if aa == "ASN" else "Glutamine (Gln)"
    if not _has_amide_cn_o():
      raise ValueError(
          f"[Amino Acid Side-Chain Error] {rid} ({mk} {mname!r}): Molecule is"
          f" labeled as an {aa_full} residue, but SMILES {msmi!r} does not"
          " contain a carboxamide group (C bonded to N and O). Use '*CC(=O)N'"
          " or '*C(=O)N'."
      )
  elif aa == "CYS":
    if "S" not in elems:
      raise ValueError(
          f"[Amino Acid Side-Chain Error] {rid} ({mk} {mname!r}): Molecule is"
          f" labeled as a Cysteine (Cys) residue, but SMILES {msmi!r} does not"
          " contain a sulfur ('S') atom. Use '*CS' or '*S' for thiol and"
          " '*C[S-]' or '*[S-]' for thiolate."
      )
  elif aa in ("SER", "THR"):
    aa_full = "Serine (Ser)" if aa == "SER" else "Threonine (Thr)"
    has_sp3_or_dummy_o = False
    for a in mol.GetAtoms():
      if a.GetSymbol() == "O":
        for nb in a.GetNeighbors():
          if nb.GetSymbol() == "*" or (
              nb.GetSymbol() == "C"
              and sum(1 for x in nb.GetNeighbors() if x.GetSymbol() == "O") == 1
              and not nb.GetIsAromatic()
          ):
            has_sp3_or_dummy_o = True
            break
    if not has_sp3_or_dummy_o or _has_phenol_ring() or "S" in elems:
      raise ValueError(
          f"[Amino Acid Side-Chain Error] {rid} ({mk} {mname!r}): Molecule is"
          f" labeled as a {aa_full} residue, but SMILES {msmi!r} does not"
          " represent an aliphatic hydroxyl/alkoxide/ester side chain. Use"
          " '*CO' or '*O' (Ser) / '*C(C)O' (Thr) for alcohol and '*C[O-]' or"
          " '*[O-]' for alkoxide."
      )
  elif aa == "TYR":
    if not _has_phenol_ring():
      raise ValueError(
          f"[Amino Acid Side-Chain Error] {rid} ({mk} {mname!r}): Molecule is"
          f" labeled as a Tyrosine (Tyr) residue, but SMILES {msmi!r} does not"
          " contain a 6-membered phenolic/phenolate ring. Use '*c1ccc(O)cc1'"
          " for neutral Tyr-OH and '*c1ccc([O-])cc1' for tyrosinate Tyr-O-,"
          " not a bare aliphatic alcohol ('*O' / '*[O-]')."
      )
  elif aa == "TRP":
    has_indole = (
        any(len(r) == 5 for r in rings)
        and any(len(r) == 6 for r in rings)
        and "N" in elems
    )
    if not has_indole:
      raise ValueError(
          f"[Amino Acid Side-Chain Error] {rid} ({mk} {mname!r}): Molecule is"
          f" labeled as a Tryptophan (Trp) residue, but SMILES {msmi!r} does"
          " not contain a bicyclic indole ring ('*c1c[nH]c2ccccc12')."
      )
  elif aa == "PHE":
    has_phenyl = any(
        len(r) == 6 and all(mol.GetAtomWithIdx(i).GetSymbol() == "C" for i in r)
        for r in rings
    )
    if not has_phenyl or _has_phenol_ring():
      raise ValueError(
          f"[Amino Acid Side-Chain Error] {rid} ({mk} {mname!r}): Molecule is"
          f" labeled as a Phenylalanine (Phe) residue, but SMILES {msmi!r}"
          " does not contain an unsubstituted phenyl ring ('*Cc1ccccc1')."
      )
  elif aa == "MET":
    has_thioether = False
    for a in mol.GetAtoms():
      if a.GetSymbol() == "S":
        c_nbs = [nb for nb in a.GetNeighbors() if nb.GetSymbol() in ("C", "*")]
        if len(c_nbs) >= 2:
          has_thioether = True
    if not has_thioether:
      raise ValueError(
          f"[Amino Acid Side-Chain Error] {rid} ({mk} {mname!r}): Molecule is"
          f" labeled as a Methionine (Met) residue, but SMILES {msmi!r} does"
          " not contain a thioether/sulfonium group ('*CCSC' or '*CSC')."
      )
  elif aa == "PRO":
    has_pyrrolidine = any(
        len(r) == 5
        and [mol.GetAtomWithIdx(i).GetSymbol() for i in r].count("N") == 1
        and [mol.GetAtomWithIdx(i).GetSymbol() for i in r].count("C") == 4
        for r in rings
    )
    if not has_pyrrolidine:
      raise ValueError(
          f"[Amino Acid Side-Chain Error] {rid} ({mk} {mname!r}): Molecule is"
          f" labeled as a Proline (Pro) residue, but SMILES {msmi!r} does not"
          " contain a 5-membered pyrrolidine/prolyl ring."
      )


def _strip_explicit_hydrogens(
    mol: Chem.Mol, clear_atom_maps: bool = False
) -> Tuple[Chem.Mol, List[int]]:
  """Returns a copy of `mol` without explicit hydrogens and heavy->orig indices."""
  heavy_to_orig = [
      i
      for i in range(mol.GetNumAtoms())
      if mol.GetAtomWithIdx(i).GetAtomicNum() != 1
  ]
  rw_mol = Chem.RWMol(mol)
  for idx in range(rw_mol.GetNumAtoms() - 1, -1, -1):
    if rw_mol.GetAtomWithIdx(idx).GetAtomicNum() == 1:
      rw_mol.RemoveAtom(idx)
  if clear_atom_maps:
    for atom in rw_mol.GetAtoms():
      atom.SetAtomMapNum(0)
  return rw_mol.GetMol(), heavy_to_orig


def _parse_graph_atom_origin(aid_str: str) -> Tuple[str, Optional[int]]:
  """Parses a graph atom ID (e.g. 'r1_3') into (reactant_key, atom_idx)."""
  if "_" not in aid_str:
    return ("", None)
  mol_key, idx_str = aid_str.rsplit("_", 1)
  return (mol_key, int(idx_str) if idx_str.isdigit() else None)


def _build_product_kekule_candidate(
    match: Tuple[int, ...],
    p_mol: Chem.Mol,
    p_heavy: Chem.Mol,
    p_heavy_to_orig: List[int],
    prod_heavy_kek: Chem.Mol,
    prod_heavy_origins: List[Tuple[str, Optional[int]]],
    r_partner_key: str,
) -> Tuple[List[Tuple[int, int]], Dict[Tuple[int, int], Chem.BondType]]:
  """Extracts (reactant-product atom pairs, Kekule bond map) for one match."""
  atom_pairs = []
  for u_h, prod_h_idx in enumerate(match):
    rk_s, r_idx = prod_heavy_origins[prod_h_idx]
    if rk_s == r_partner_key and r_idx is not None:
      atom_pairs.append((r_idx, p_heavy_to_orig[u_h]))

  kek_bonds = {}
  for bond in p_heavy.GetBonds():
    u_h, v_h = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    u_orig, v_orig = p_heavy_to_orig[u_h], p_heavy_to_orig[v_h]
    prod_bond = prod_heavy_kek.GetBondBetweenAtoms(match[u_h], match[v_h])
    if prod_bond is not None:
      edge = (min(u_orig, v_orig), max(u_orig, v_orig))
      kek_bonds[edge] = prod_bond.GetBondType()

  for bond in p_mol.GetBonds():
    edge = (
        min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
        max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
    )
    kek_bonds.setdefault(edge, Chem.BondType.SINGLE)

  return atom_pairs, kek_bonds


def _extract_product_kekule_map(
    prod_g: MolecularGraph,
    product_items: List[Tuple[str, Dict[str, Any]]],
    rid: str,
) -> Tuple[
    Dict[
        str,
        List[
            Tuple[List[Tuple[int, int]], Dict[Tuple[int, int], Chem.BondType]]
        ],
    ],
    Dict[str, str],
]:
  """Extracts post-arrow Kekule bond orders and reactant partners for products."""
  if not product_items:
    return {}, {}
  try:
    prod_aids = list(prod_g.atoms.keys())
    prod_kek_orig = prod_g.to_rdkit(sanitize=False)
    prod_heavy_kek, heavy_indices = _strip_explicit_hydrogens(prod_kek_orig)
    prod_heavy_origins = [
        _parse_graph_atom_origin(prod_aids[i]) for i in heavy_indices
    ]
    prod_heavy_arom = Chem.Mol(prod_heavy_kek)
    Chem.SanitizeMol(prod_heavy_arom)
    atom_to_comp_size: Dict[int, int] = {}
    for frag in Chem.GetMolFrags(prod_heavy_arom):
      for f_idx in frag:
        atom_to_comp_size[f_idx] = len(frag)

    p_kekule_map = {}
    p_ref_map: Dict[str, str] = {}
    for p_k, p_cfg in product_items:
      p_smi = p_cfg.get("smi")
      if not p_smi:
        continue
      p_mol = _parse_smiles_preserving_mapped_hs(p_smi)
      if p_mol is None:
        continue
      p_heavy, p_heavy_to_orig = _strip_explicit_hydrogens(
          p_mol, clear_atom_maps=True
      )
      Chem.SanitizeMol(p_heavy)
      matches = prod_heavy_arom.GetSubstructMatches(p_heavy, uniquify=False)
      if not matches:
        continue
      p_frags = Chem.GetMolFrags(p_heavy)
      exact_matches = [
          m_tup
          for m_tup in matches
          if all(
              atom_to_comp_size.get(m_tup[pf[0]]) == len(pf) for pf in p_frags
          )
      ]
      if exact_matches:
        matches = exact_matches
      origin_counts: Dict[str, int] = {}
      for prod_h_idx in matches[0]:
        rk_s, r_idx_s = prod_heavy_origins[prod_h_idx]
        if rk_s and r_idx_s is not None:
          origin_counts[rk_s] = origin_counts.get(rk_s, 0) + 1
      r_partner_key = (
          max(origin_counts, key=origin_counts.get)
          if origin_counts
          else f"r{p_k[1:]}"
      )
      if origin_counts.get(r_partner_key, 0) >= 3:
        p_ref_map[p_k] = r_partner_key
      p_kekule_map[p_k] = [
          _build_product_kekule_candidate(
              match,
              p_mol,
              p_heavy,
              p_heavy_to_orig,
              prod_heavy_kek,
              prod_heavy_origins,
              r_partner_key,
          )
          for match in matches
      ]
    return p_kekule_map, p_ref_map
  except (Chem.MolSanitizeException, RuntimeError, ValueError) as exc:
    raise ValueError(
        f"{rid}: Failed to extract post-arrow Kekulé bond orders for products:"
        f" {exc}"
    ) from exc


def verify_mechanism_with_electron_graph(rxn_def, spec):
  """Formally verifies the elementary reaction using ElectronGraph prior to rendering.

  Enforces:
    - Topological locality & non-teleportation
    - Lewis octet rule and non-hypervalency
    - Non-negative lone pairs and radicals
    - Tier 5 elementary-step single-proton transfer invariant
    - Stereochemical consistency (cis/Z and trans/E SMILES vs name)
    - Target product matching (when products are specified)
  """
  rid = rxn_def.get("id", "reaction")
  for side_prefix in ("r", "p"):
    for mk, mcfg in _extract_molecules(rxn_def, side_prefix):
      mname = str(mcfg.get("name", ""))
      msmi = str(mcfg.get("smi", ""))
      is_cis_named = bool(
          re.search(r"(?i)((?<!s-)\bcis\s*-|\(\d*z\)\s*-)", mname)
      )
      is_trans_named = bool(
          re.search(r"(?i)((?<!s-)\btrans\s*-|\(\d*e\)\s*-)", mname)
      )
      if is_cis_named or is_trans_named:
        mol_chk = Chem.MolFromSmiles(msmi)
        if mol_chk is not None:
          acyclic_cc_dbls = [
              b
              for b in mol_chk.GetBonds()
              if b.GetBondType() == Chem.BondType.DOUBLE
              and not b.IsInRing()
              and b.GetBeginAtom().GetAtomicNum() == 6
              and b.GetEndAtom().GetAtomicNum() == 6
          ]
          if acyclic_cc_dbls and ("/" not in msmi and "\\" not in msmi):
            raise ValueError(
                f"[Stereochemistry Error] {rid} ({mk} {mname!r}): Name"
                " specifies cis/(Z) or trans/(E) stereochemistry, but SMILES"
                f" {msmi!r} omits directional bond stereo ('/' and '\\\\')."
                " Unspecified acyclic C=C bonds default to s-trans zig-zag."
            )
          if acyclic_cc_dbls:
            rdDepictor.Compute2DCoords(mol_chk)
            conf_chk = mol_chk.GetConformer()
            for b_cc in acyclic_cc_dbls:
              u_i = b_cc.GetBeginAtomIdx()
              v_i = b_cc.GetEndAtomIdx()
              het_u = [
                  nb.GetIdx()
                  for nb in mol_chk.GetAtomWithIdx(u_i).GetNeighbors()
                  if nb.GetIdx() != v_i and nb.GetAtomicNum() not in (1, 6)
              ]
              het_v = [
                  nb.GetIdx()
                  for nb in mol_chk.GetAtomWithIdx(v_i).GetNeighbors()
                  if nb.GetIdx() != u_i and nb.GetAtomicNum() not in (1, 6)
              ]
              if len(het_u) == 1 and len(het_v) == 1:
                pu = conf_chk.GetAtomPosition(u_i)
                pv = conf_chk.GetAtomPosition(v_i)
                psu = conf_chk.GetAtomPosition(het_u[0])
                psv = conf_chk.GetAtomPosition(het_v[0])
                nx_uv = -(pv.y - pu.y)
                ny_uv = pv.x - pu.x
                side_u = (psu.x - pu.x) * nx_uv + (psu.y - pu.y) * ny_uv
                side_v = (psv.x - pv.x) * nx_uv + (psv.y - pv.y) * ny_uv
                if is_cis_named and (side_u * side_v < 0):
                  raise ValueError(
                      f"[Stereochemistry Error] {rid} ({mk} {mname!r}): Name"
                      f" specifies cis/(Z) stereochemistry, but SMILES {msmi!r}"
                      " places the C=C heteroatom substituents on opposite"
                      " (trans/E) faces."
                  )
                if is_trans_named and (side_u * side_v > 0):
                  raise ValueError(
                      f"[Stereochemistry Error] {rid} ({mk} {mname!r}): Name"
                      " specifies trans/(E) stereochemistry, but SMILES"
                      f" {msmi!r} places the C=C heteroatom substituents on the"
                      " same (cis/Z) face."
                  )

      # Comprehensive Amino Acid Side-Chain Chemical Identity Validator
      _validate_amino_acid_sidechain(rid, mk, mcfg)

  allow_multiproton = bool(
      (spec or {}).get("allow_multiproton") or rxn_def.get("allow_multiproton")
  )
  mg, arrows = build_graph_and_arrows_from_spec(rxn_def, spec)
  valid, prod_g, audit = ReactionVerifierEngine.apply_and_verify_step(
      mg, arrows, step_name=rid, allow_multiproton=allow_multiproton
  )
  if not valid:
    err_msg = (
        f"[ELECTRON_GRAPH ERROR] Verification failed for {rid}:\n"
        + "\n".join(f"  - {e}" for e in audit.get("errors", []))
    )
    raise ValueError(err_msg)

  product_items = _extract_molecules(rxn_def, "p")
  if product_items:
    target_smis = [p["smi"] for _, p in product_items if p.get("smi")]
    if target_smis:
      target_smi = ".".join(target_smis)
      skip_prod_match = bool(
          (spec or {}).get("skip_product_match")
          or rxn_def.get("skip_product_match")
          or (spec or {}).get("partial_products")
          or rxn_def.get("partial_products")
      )
      if not skip_prod_match:
        try:
          num_prod_components = nx.number_connected_components(prod_g.graph)
          target_mg = MolecularGraph.from_smiles(target_smi)
          r_counts = mg.nuclei_multiset()
          p_counts = target_mg.nuclei_multiset()
          has_excess_atoms = any(
              p_counts[k] > r_counts.get(k, 0) for k in p_counts
          )
          if has_excess_atoms or (
              len(product_items) >= num_prod_components and r_counts != p_counts
          ):
            raise ValueError(
                "[Stoichiometry Error] Stoichiometry mismatch between"
                f" reactants and target products for {rid}:\n"
                f"  Reactant atoms: {dict(r_counts)}\n"
                f"  Target product atoms: {dict(p_counts)}\n"
                f"  Target SMILES: {target_smi}"
            )
          prod_rads = sum(
              a.radicals + a.unpaired_d for a in prod_g.atoms.values()
          )
          target_rads = sum(
              a.radicals + a.unpaired_d for a in target_mg.atoms.values()
          )
          if (
              len(product_items) >= num_prod_components
              and prod_rads != target_rads
          ):
            raise ValueError(
                f"[Stoichiometry Error] Radical count mismatch for {rid}:"
                f" generated product has {prod_rads} radical electron(s), but"
                f" target product SMILES ({target_smi}) has {target_rads}."
            )
        except ValueError:
          raise
        except Exception:
          pass

      match, msg = ReactionVerifierEngine.verify_product_match(
          prod_g, target_smi
      )
      if not match:
        if "specifies @/@@ chirality on homotopic/achiral" in msg:
          raise ValueError(f"[ELECTRON_GRAPH ERROR] {msg}")
        prod_smi = prod_g.to_smiles()

        def _canon_smi_iso(s: str) -> str:
          m_iso = Chem.MolFromSmiles(s)
          if m_iso is None:
            return s
          try:
            m_iso = Chem.RemoveHs(m_iso)
          except (
              Chem.rdchem.MolSanitizeException,
              ValueError,
              RuntimeError,
          ) as exc:
            raise ValueError(
                "[ELECTRON_GRAPH ERROR] Failed to remove explicit hydrogens"
                f" when canonicalizing stereoisomer {s!r} for {rid}: {exc}"
            ) from exc
          return Chem.MolToSmiles(m_iso, canonical=True, isomericSmiles=True)

        has_explicit_policy = any(
            arr.stereo_policy and arr.stereo_policy != "auto" for arr in arrows
        )
        if not skip_prod_match and (
            has_explicit_policy
            or ("@" in prod_smi and any("@" in s for s in target_smis if s))
        ):
          prod_parts_iso = [_canon_smi_iso(s) for s in prod_smi.split(".") if s]
          target_parts_iso = [
              _canon_smi_iso(s)
              for s in target_smis
              if s and (has_explicit_policy or "@" in s)
          ]
          missing_iso = [
              tp for tp in target_parts_iso if tp not in prod_parts_iso
          ]
          if missing_iso:
            raise ValueError(
                "[ELECTRON_GRAPH ERROR] Target stereoisomer mismatch for"
                f" {rid}: {msg} (Generated: {prod_smi}, Expected: {target_smi})"
            )

        # Compare components canonically (non-isomeric)
        prod_parts_no_iso = [
            _canon_smi_no_iso(s) for s in prod_smi.split(".") if s
        ]
        target_parts_no_iso = [_canon_smi_no_iso(s) for s in target_smis if s]
        missing_no_iso = [
            tp for tp in target_parts_no_iso if tp not in prod_parts_no_iso
        ]

        if missing_no_iso and not skip_prod_match:
          err_msg = (
              f"[ELECTRON_GRAPH ERROR] Target product mismatch for {rid}:\n"
              f"  Missing: {missing_no_iso}\n"
              f"  Generated: {prod_smi}\n"
              f"  Expected: {target_smi}"
          )
          raise ValueError(err_msg)

  audit["r_kekule_map"] = dict(getattr(mg, "r_kekule_map", {}))
  p_kekule_map, p_ref_map = _extract_product_kekule_map(
      prod_g, product_items, rid
  )
  audit["p_kekule_map"] = p_kekule_map
  audit["p_ref_map"] = p_ref_map
  return audit


def _sample_bezier_pts_2d(p0, p1, p2, p3, n=100):
  t = np.linspace(0, 1, n)[:, None]
  return (
      (1 - t) ** 3 * p0
      + 3 * (1 - t) ** 2 * t * p1
      + 3 * (1 - t) * t**2 * p2
      + t**3 * p3
  )


def _draw_matplotlib_stealth_arrow(
    ax, p0, p1, p2, p3, color="#E11D48", base_lw=1.05, electrons=2
):
  """Renders a publication-grade, TikZ Stealth-equivalent curved electron-pushing arrow.

  Arrowhead dimensions and shaft line widths are scaled dynamically to the curve
  arc length so that short bonds (e.g., bond cleavage) never appear squeezed or
  crammed, and long trajectories retain a slender, elegant dart profile.
  """
  poly = _sample_bezier_pts_2d(p0, p1, p2, p3, n=100)
  diffs = poly[1:] - poly[:-1]
  seg_lens = np.linalg.norm(diffs, axis=1)
  arc_l = float(np.sum(seg_lens))

  # Scale arrowhead and stroke based on arc length (canvas coordinates in cm)
  # TikZ reference: chemarrow (2.0x1.3mm), fishhookarrow (harpoon 2.45x2.25mm)
  if electrons == 1:
    if arc_l < 0.45:
      head_l = min(0.16, max(0.11, arc_l * 0.34))
      head_w = head_l * 0.95
      lw = base_lw * 0.92
    else:
      head_l = min(0.22, max(0.15, arc_l * 0.25))
      head_w = head_l * 0.92
      lw = base_lw
  elif arc_l < 0.48:
    head_l = min(0.12, max(0.07, arc_l * 0.25))
    head_w = head_l * 0.65
    lw = base_lw * 0.88
  else:
    head_l = min(0.17, max(0.11, arc_l * 0.20))
    head_w = head_l * 0.62
    lw = base_lw

  # Compute smooth tangent at arrowhead tip
  tangent = poly[-1] - poly[-6] if len(poly) >= 6 else (p3 - p2)
  norm_t = np.linalg.norm(tangent)
  if norm_t < 1e-6:
    tangent = p3 - p2
    norm_t = np.linalg.norm(tangent)
  u_t = tangent / max(norm_t, 1e-6)
  u_n = np.array([-u_t[1], u_t[0]])

  p_tip = p3
  p_base = p_tip - head_l * u_t
  wing_sweep = 0.12 * head_l
  p_left = (p_base - wing_sweep * u_t) + 0.5 * head_w * u_n
  p_right = (p_base - wing_sweep * u_t) - 0.5 * head_w * u_n
  p_notch = p_base + 0.20 * head_l * u_t

  # Terminate shaft seamlessly at the indented notch of the stealth head
  dists = np.linalg.norm(poly - p_notch, axis=1)
  idx_stop = int(np.argmin(dists))
  if idx_stop < 5:
    idx_stop = 5
  poly_shaft = np.vstack([poly[:idx_stop], p_notch])

  # Draw shaft
  ax.plot(
      poly_shaft[:, 0],
      poly_shaft[:, 1],
      color=color,
      lw=lw,
      solid_capstyle="round",
      zorder=6,
  )

  # Draw stealth arrowhead polygon (place single fishhook barb on convex side)
  if electrons == 1:
    barb_pt = p_right if bezier_turns_ccw(p0, p1, p2, p3) else p_left
    arrow_poly = [p_tip, barb_pt, p_notch]
  else:
    arrow_poly = [p_tip, p_left, p_notch, p_right]

  stealth_patch = patches.Polygon(
      arrow_poly,
      closed=True,
      facecolor=color,
      edgecolor=color,
      lw=lw * 0.3,
      joinstyle="miter",
      zorder=7,
  )
  ax.add_patch(stealth_patch)


def _draw_matplotlib_reaction_arrow(
    ax, start_x, end_x, y=0.0, color="#0f172a", lw=1.2
):
  """Renders a crisp, publication-grade main transformation arrow with sleek stealth tip."""
  head_l = 0.25
  head_w = 0.16
  shaft_end_x = end_x - 0.20

  ax.plot(
      [start_x, shaft_end_x],
      [y, y],
      color=color,
      lw=lw,
      solid_capstyle="butt",
      zorder=4,
  )

  p_tip = np.array([end_x, y])
  p_base = np.array([end_x - head_l, y])
  p_left = np.array([end_x - head_l - 0.03, y + head_w / 2])
  p_right = np.array([end_x - head_l - 0.03, y - head_w / 2])
  p_notch = np.array([end_x - head_l + 0.04, y])

  arrow_poly = [p_tip, p_left, p_notch, p_right]
  head_patch = patches.Polygon(
      arrow_poly,
      closed=True,
      facecolor=color,
      edgecolor=color,
      lw=lw * 0.3,
      joinstyle="miter",
      zorder=5,
  )
  ax.add_patch(head_patch)


def _render_mechanism_matplotlib(
    rid,
    title,
    reactant_canvases,
    product_canvases,
    plus_xs,
    product_plus_xs,
    rxn_arrow_start_x,
    rxn_arrow_end_x,
    cond,
    audit_arrows,
    png_final,
):
  """High-resolution fallback renderer for reaction mechanisms with curved electron-pushing arrows when pdflatex is unavailable."""
  xs = [float(rxn_arrow_start_x), float(rxn_arrow_end_x)]
  ys = [0.0]

  all_canvases = []
  for k, r_cfg, canvas, center_x in reactant_canvases:
    aspect = canvas.h_px / max(canvas.w_px, 1)
    h_cm = canvas.w_cm * aspect
    all_canvases.append(
        (canvas, center_x, canvas.w_cm, h_cm, r_cfg.get("name", ""))
    )
    xs.extend([center_x - canvas.w_cm / 2, center_x + canvas.w_cm / 2])
    ys.extend([-h_cm / 2, h_cm / 2])

  for k, p_cfg, canvas, center_x in product_canvases:
    aspect = canvas.h_px / max(canvas.w_px, 1)
    h_cm = canvas.w_cm * aspect
    all_canvases.append(
        (canvas, center_x, canvas.w_cm, h_cm, p_cfg.get("name", ""))
    )
    xs.extend([center_x - canvas.w_cm / 2, center_x + canvas.w_cm / 2])
    ys.extend([-h_cm / 2, h_cm / 2])

  for px in plus_xs + product_plus_xs:
    xs.append(px)

  for a_label, pts in audit_arrows:
    for pt in pts:
      xs.append(float(pt[0]))
      ys.append(float(pt[1]))

  x_min = min(xs) - 0.8
  x_max = max(xs) + 0.8
  y_min = min(ys) - 1.2
  y_max = max(ys) + (1.2 if title else 0.8)

  w_in = max((x_max - x_min) / 2.54, 5.0)
  h_in = max((y_max - y_min) / 2.54, 3.0)

  fig, ax = plt.subplots(figsize=(w_in, h_in), dpi=300)
  ax.set_xlim(x_min, x_max)
  ax.set_ylim(y_min, y_max)
  ax.set_aspect("equal")
  ax.axis("off")

  # Render each molecule image
  for canvas, center_x, w_cm, h_cm, c_name in all_canvases:
    if os.path.exists(canvas.png_path):
      img = plt.imread(canvas.png_path)
      ax.imshow(
          img,
          extent=[
              center_x - w_cm / 2,
              center_x + w_cm / 2,
              -h_cm / 2,
              h_cm / 2,
          ],
          aspect="equal",
          zorder=2,
      )
    if c_name:
      clean_name = (
          re.sub(r"\\[a-zA-Z]+", "", str(c_name))
          .replace("{", "")
          .replace("}", "")
          .strip()
      )
      ax.text(
          center_x,
          -h_cm / 2 - 0.35,
          clean_name,
          ha="center",
          va="top",
          fontsize=9,
          fontweight="bold",
          color="#0f172a",
          zorder=3,
      )

  # Render plus signs
  for px in plus_xs + product_plus_xs:
    ax.text(
        px,
        0,
        "+",
        ha="center",
        va="center",
        fontsize=20,
        fontweight="bold",
        color="#334155",
        zorder=3,
    )

  # Render main reaction arrow with sleek stealth tip
  _draw_matplotlib_reaction_arrow(ax, rxn_arrow_start_x, rxn_arrow_end_x, y=0.0)

  # Render reaction conditions
  top_c = cond.get("top", "")
  bot_c = cond.get("bottom", "")
  if top_c:
    clean_top = (
        re.sub(r"\\[a-zA-Z]+", "", str(top_c))
        .replace("{", "")
        .replace("}", "")
        .replace("$", "")
        .strip()
    )
    ax.text(
        (rxn_arrow_start_x + rxn_arrow_end_x) / 2,
        0.28,
        clean_top,
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="#1e293b",
        zorder=4,
    )
  if bot_c:
    clean_bot = (
        re.sub(r"\\[a-zA-Z]+", "", str(bot_c))
        .replace("{", "")
        .replace("}", "")
        .replace("$", "")
        .strip()
    )
    ax.text(
        (rxn_arrow_start_x + rxn_arrow_end_x) / 2,
        -0.28,
        clean_bot,
        ha="center",
        va="top",
        fontsize=8.5,
        color="#1e293b",
        zorder=4,
    )

  # Render electron-pushing curved arrows with sleek stealth arrowheads
  for a_label, (p0, p1, p2, p3) in audit_arrows:
    lbl_low = str(a_label).lower()
    electrons_cnt = (
        1
        if "fishhook" in lbl_low or "single" in lbl_low or "1e" in lbl_low
        else 2
    )
    _draw_matplotlib_stealth_arrow(
        ax,
        p0,
        p1,
        p2,
        p3,
        color="#E11D48",
        base_lw=1.05,
        electrons=electrons_cnt,
    )

  if title:
    clean_title = (
        re.sub(r"\\[a-zA-Z]+", "", str(title))
        .replace("{", "")
        .replace("}", "")
        .replace("$", "")
        .strip()
    )
    fig.suptitle(
        clean_title, fontsize=11, fontweight="bold", y=0.96, color="#0f172a"
    )

  os.makedirs(os.path.dirname(png_final), exist_ok=True)
  plt.savefig(png_final, dpi=300, bbox_inches="tight", facecolor="white")
  plt.close(fig)
  return True


def _resolve_mol_and_idx(
    default_key: str, tok: Any, mol_map: dict[str, Any]
) -> tuple[str, int]:
  """Resolves an atom token across `mol_map`, honoring `[mol, tok]` or `:N`."""
  m_key = default_key
  if (
      isinstance(tok, (tuple, list))
      and len(tok) == 2
      and isinstance(tok[0], str)
      and tok[0] in mol_map
  ):
    m_key, tok = tok[0], tok[1]
  if m_key not in mol_map:
    raise KeyError(m_key)
  if isinstance(tok, str) and tok.startswith(":") and tok[1:].isdigit():
    want_map = int(tok[1:])
    if want_map not in mol_map[m_key].map_to_idx:
      alt = [
          k2
          for k2, c2 in mol_map.items()
          if not k2.startswith("p") and want_map in c2.map_to_idx
      ]
      if len(alt) == 1:
        m_key = alt[0]
  return m_key, resolve_atom_token(mol_map[m_key], tok)


def _extract_nascent_1e_candidates(
    a_spec: dict[str, Any], mol_map: dict[str, Any]
) -> list[tuple[frozenset[tuple[str, int]], str, Any, Any, tuple[str, int]]]:
  """Extracts candidate nascent-bond sinks for an unpaired 1e fishhook arrow."""
  if int(a_spec.get("electrons", 2)) != 1 or a_spec.get("side") is not None:
    return []
  atype = a_spec.get("type")
  try:
    if atype == "attack":
      src = a_spec.get("src")
      tgt = a_spec.get("tgt")
      if not (src and tgt) or isinstance(tgt[1], (list, tuple)):
        return []
      s_k, s_tok = src
      t_k, t_tok = tgt
      if s_k not in mol_map or t_k not in mol_map:
        return []
      s_mol, t_mol = mol_map[s_k], mol_map[t_k]
      t_idx = resolve_atom_token(t_mol, t_tok)
      if not isinstance(s_tok, (list, tuple)):
        s_idx = resolve_atom_token(s_mol, s_tok)
        if (s_k, s_idx) == (t_k, t_idx) or (
            s_k == t_k and s_mol.find_bond(s_idx, t_idx) is not None
        ):
          return []
        ep_s = (s_k, s_idx)
        return [(frozenset((ep_s, (t_k, t_idx))), "atom", s_mol, s_idx, ep_s)]
      if len(s_tok) != 2:
        return []
      u_i = resolve_atom_token(s_mol, s_tok[0])
      v_i = resolve_atom_token(s_mol, s_tok[1])
      if s_mol.find_bond(u_i, v_i) is None:
        return []
      c_cands = (
          [resolve_atom_token(s_mol, a_spec["conn"])]
          if a_spec.get("conn") is not None
          else [u_i, v_i]
      )
      out = []
      for c_i in c_cands:
        if (s_k, c_i) != (t_k, t_idx) and not (
            s_k == t_k and s_mol.find_bond(c_i, t_idx) is not None
        ):
          ep_s = (s_k, c_i)
          out.append((
              frozenset((ep_s, (t_k, t_idx))),
              "bond",
              s_mol,
              (u_i, v_i),
              ep_s,
          ))
      return out

    if atype == "bond_to_bond":
      s_k = a_spec.get("src_mol") or a_spec.get("mol", "r1")
      t_k = a_spec.get("tgt_mol") or a_spec.get("mol", s_k)
      if s_k not in mol_map or t_k not in mol_map:
        return []
      m = mol_map[s_k]
      b1_raw = a_spec.get("b1") or a_spec.get("src_bond")
      b2_raw = a_spec.get("b2") or a_spec.get("tgt_bond")
      if not (b1_raw and b2_raw):
        return []
      b1 = (
          resolve_atom_token(m, b1_raw[0]),
          resolve_atom_token(m, b1_raw[1]),
      )
      if m.find_bond(*b1) is None:
        return []
      ep_b2_0 = _resolve_mol_and_idx(t_k, b2_raw[0], mol_map)
      ep_b2_1 = _resolve_mol_and_idx(t_k, b2_raw[1], mol_map)
      if ep_b2_0 == ep_b2_1 or (
          ep_b2_0[0] == ep_b2_1[0]
          and mol_map[ep_b2_0[0]].find_bond(ep_b2_0[1], ep_b2_1[1]) is not None
      ):
        return []
      shared_eps = {(s_k, b1[0]), (s_k, b1[1])} & {ep_b2_0, ep_b2_1}
      if len(shared_eps) == 1:
        ep_s = next(iter(shared_eps))
      elif s_k == ep_b2_0[0] == ep_b2_1[0]:
        b1_mid = m.find_bond(*b1)["mid"]
        ep_s = min(
            (ep_b2_0, ep_b2_1),
            key=lambda ep: float(
                np.linalg.norm(m.atoms[ep[1]]["pos"] - b1_mid)
            ),
        )
      else:
        return []
      ep_other = ep_b2_1 if ep_b2_0 == ep_s else ep_b2_0
      return [(frozenset((ep_s, ep_other)), "bond", m, b1, ep_s)]

    if atype == "lone_pair_to_bond":
      m_k = a_spec.get("mol", "r1")
      bond_raw = a_spec.get("bond")
      if m_k not in mol_map or not bond_raw or "atom" not in a_spec:
        return []
      ep_a = _resolve_mol_and_idx(m_k, a_spec["atom"], mol_map)
      ep_b0 = _resolve_mol_and_idx(m_k, bond_raw[0], mol_map)
      ep_b1 = _resolve_mol_and_idx(m_k, bond_raw[1], mol_map)
      if (
          ep_b0 == ep_b1
          or ep_a not in (ep_b0, ep_b1)
          or (
              ep_b0[0] == ep_b1[0]
              and mol_map[ep_b0[0]].find_bond(ep_b0[1], ep_b1[1]) is not None
          )
      ):
        return []
      ep_other = ep_b1 if ep_b0 == ep_a else ep_b0
      return [
          (frozenset((ep_a, ep_other)), "atom", mol_map[ep_a[0]], ep_a[1], ep_a)
      ]
  except (KeyError, ValueError, TypeError, IndexError):
    return []
  return []


def match_nascent_1e_pairs(
    arrow_specs: list[dict[str, Any]], mol_map: dict[str, Any]
) -> dict[int, tuple[int, tuple[Any, ...], tuple[Any, ...]]]:
  """Matches pairs of 1e arrows contributing to the same nascent bond."""
  cands_by_idx = {
      idx: _extract_nascent_1e_candidates(a_s, mol_map)
      for idx, a_s in enumerate(arrow_specs)
  }
  paired: dict[int, tuple[int, tuple[Any, ...], tuple[Any, ...]]] = {}
  used: set[int] = set()
  for i, cands_i in cands_by_idx.items():
    if i in used or not cands_i:
      continue
    found = False
    for ci in cands_i:
      for j, cands_j in cands_by_idx.items():
        if j <= i or j in used:
          continue
        for cj in cands_j:
          if ci[0] == cj[0] and ci[4] != cj[4]:
            paired[i] = (j, ci, cj)
            paired[j] = (i, cj, ci)
            used.update((i, j))
            found = True
            break
        if found:
          break
      if found:
        break
  return paired


def _bezier_has_inflection(
    bez: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
) -> bool:
  """Returns True if cubic Bézier `bez` loops backward or changes turn sign."""
  c0, c1, c2, c3 = bez
  if float(np.dot(c2 - c1, c3 - c0)) < 0.0:
    return True
  pts = sample_bezier_pts(c0, c1, c2, c3, n=20, t_min=0.0, t_max=1.0)
  dp = pts[1:] - pts[:-1]
  tr = dp[:-1, 0] * dp[1:, 1] - dp[:-1, 1] * dp[1:, 0]
  return bool(float(np.min(tr)) < -1e-5 and float(np.max(tr)) > 1e-5)


def _tail_clearance_below(
    poly: np.ndarray,
    segs: list[tuple[np.ndarray, np.ndarray]],
    min_cm: float = 0.060,
) -> bool:
  """Returns True if departure tail of `poly` is within `min_cm` of `segs`."""
  for pt_t in (poly[0], poly[1]):
    for bp1, bp2 in segs:
      if dist_pt_to_seg(pt_t, bp1, bp2) < min_cm:
        return True
  return False


def _compute_paired_fishhook_beziers(
    desc_u: tuple[Any, ...],
    desc_v: tuple[Any, ...],
    mol_map: dict[str, Any],
    scene_canvases: list[Any],
    plus_xs: list[float],
    prior_polys: list[np.ndarray],
) -> (
    tuple[
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    ]
    | None
):
  """Computes symmetric midpoint-meeting 1e fishhooks for a nascent bond."""
  _, kind_u, mol_u, src_u, ep_u = desc_u
  _, kind_v, mol_v, src_v, ep_v = desc_v
  p_u = mol_map[ep_u[0]].atoms[ep_u[1]]["pos"]
  p_v = mol_map[ep_v[0]].atoms[ep_v[1]]["pos"]
  v_uv = p_v - p_u
  l_uv = float(np.linalg.norm(v_uv))
  if l_uv < 0.45:
    return None
  u_uv = v_uv / l_uv
  n_base = np.array([-u_uv[1], u_uv[0]])

  def _filtered_bond_segs(kind_s, mol_s, src_id, *, exclude_incident: bool):
    inc_ids = {src_id} if kind_s == "atom" else set(src_id)
    segs = []
    for c_obj in scene_canvases:
      for b_d in c_obj.bonds:
        if c_obj is mol_s:
          if exclude_incident and (
              b_d["a1"] in inc_ids or b_d["a2"] in inc_ids
          ):
            continue
          if (
              not exclude_incident
              and kind_s == "bond"
              and {b_d["a1"], b_d["a2"]} == inc_ids
          ):
            continue
        segs.extend(b_d["stroke_segments"])
    return segs

  non_inc_u = _filtered_bond_segs(kind_u, mol_u, src_u, exclude_incident=True)
  non_inc_v = _filtered_bond_segs(kind_v, mol_v, src_v, exclude_incident=True)
  non_src_u = _filtered_bond_segs(kind_u, mol_u, src_u, exclude_incident=False)
  non_src_v = _filtered_bond_segs(kind_v, mol_v, src_v, exclude_incident=False)

  def _departures(kind_s, mol_s, src_id, p_aim, arch_n):
    out = []
    if kind_s == "atom":
      p_a = mol_s.atoms[src_id]["pos"]
      r_a = mol_s.atoms[src_id]["glyph_r"]
      nbr_angs = [
          math.atan2(
              mol_s.atoms[nb.GetIdx()]["pos"][1] - p_a[1],
              mol_s.atoms[nb.GetIdx()]["pos"][0] - p_a[0],
          )
          for nb in mol_s.mol.GetAtomWithIdx(src_id).GetNeighbors()
      ]
      v_aim = p_aim - p_a
      u_aim = v_aim / max(float(np.linalg.norm(v_aim)), 1e-6)
      raw_u = []
      for w_n in (0.42, 0.68, 0.18):
        u_c = (1.0 - w_n) * u_aim + w_n * arch_n
        n_c = float(np.linalg.norm(u_c))
        if n_c > 1e-6:
          raw_u.append(u_c / n_c)
      for b_deg in _open_sector_bisector_degs(nbr_angs):
        b_rad = math.radians(b_deg)
        u_b = np.array([math.cos(b_rad), math.sin(b_rad)])
        if float(np.dot(u_b, u_aim)) > -0.10:
          raw_u.append(u_b)
      pref_dir = 0.65 * u_aim + 0.35 * arch_n
      pref_dir = pref_dir / max(float(np.linalg.norm(pref_dir)), 1e-6)
      for u_c in raw_u:
        c_ang = math.atan2(u_c[1], u_c[0])
        if nbr_angs and min(
            abs((c_ang - na + math.pi) % (2.0 * math.pi) - math.pi)
            for na in nbr_angs
        ) < math.radians(20.0):
          continue
        p0_c = p_a + (r_a + 0.055) * u_c
        pen_c = 4.0 * (1.0 - float(np.dot(u_c, pref_dir)))
        out.append((p0_c, u_c, pen_c))
      out.sort(key=lambda item: item[2])
      return out[:3]

    b_obj = mol_s.find_bond(*src_id)
    if b_obj is None:
      return []
    n_b = b_obj["norm"]
    n_ext = _ring_exterior_normal(mol_s, src_id[0], src_id[1], n_b)
    if n_ext is not None:
      bn_cands = [(n_ext, 0.0)]
    else:
      pref_bn = n_b if float(np.dot(p_aim - b_obj["mid"], n_b)) >= 0.0 else -n_b
      bn_cands = [(pref_bn, 0.0), (-pref_bn, 12.0)]
    for dep_n, bn_pen in bn_cands:
      outer_off = bond_outer_offset_cm(b_obj, dep_n)
      p0_c = b_obj["mid"] + (outer_off + 0.085) * dep_n
      v_aim = p_aim - p0_c
      u_aim = v_aim / max(float(np.linalg.norm(v_aim)), 1e-6)
      t_dir = 0.70 * u_aim + 0.30 * arch_n
      t_dir = t_dir / max(float(np.linalg.norm(t_dir)), 1e-6)
      for w_bn in (0.65, 0.35):
        u_c = w_bn * dep_n + (1.0 - w_bn) * t_dir
        n_c = float(np.linalg.norm(u_c))
        if n_c < 1e-6:
          continue
        u_c = u_c / n_c
        if abs(float(np.dot(u_c, b_obj["tangent"]))) > 0.94:
          continue
        if n_ext is not None and float(np.dot(u_c, n_ext)) < -0.25:
          continue
        out.append((p0_c, u_c, bn_pen + 2.0 * abs(w_bn - 0.65)))
    out.sort(key=lambda item: item[2])
    return out[:3]

  h_meets = (
      (0.14, 0.26, 0.05, 0.38, 0.52, 0.66)
      if plus_xs
      else (0.14, 0.26, 0.05, 0.38)
  )
  survivors = []
  for sign in (1.0, -1.0):
    arch_n = sign * n_base
    arch_bonus = 2.5 * float(arch_n[1])
    for h_idx, h_meet in enumerate(h_meets):
      for s_shift in (0.0, -0.07, 0.07):
        p_meet = 0.5 * (p_u + p_v) + s_shift * l_uv * u_uv + h_meet * arch_n
        gap = min(0.125, max(0.095, 0.13 * l_uv))
        p3_u = p_meet - gap * u_uv
        p3_v = p_meet + gap * u_uv
        deps_u = _departures(
            kind_u, mol_u, src_u, p3_u + 0.18 * l_uv * arch_n, arch_n
        )
        deps_v = _departures(
            kind_v, mol_v, src_v, p3_v + 0.18 * l_uv * arch_n, arch_n
        )
        if not deps_u or not deps_v:
          continue
        for theta_deg in (48.0, 60.0):
          th_r = math.radians(theta_deg)
          u_app_u = -math.cos(th_r) * u_uv + math.sin(th_r) * arch_n
          u_app_v = math.cos(th_r) * u_uv + math.sin(th_r) * arch_n
          for p0_u, u_dep_u, pen_u in deps_u:
            d_u = float(np.linalg.norm(p3_u - p0_u))
            if d_u < 0.28:
              continue
            for p0_v, u_dep_v, pen_v in deps_v:
              d_v = float(np.linalg.norm(p3_v - p0_v))
              if d_v < 0.28:
                continue
              for lift_s in (0.38, 0.54):
                l1_u = max(0.16, min(2.4, d_u * lift_s))
                l2_u = max(0.16, min(1.8, d_u * (lift_s * 0.85)))
                l1_v = max(0.16, min(2.4, d_v * lift_s))
                l2_v = max(0.16, min(1.8, d_v * (lift_s * 0.85)))
                bez_u = (
                    p0_u,
                    p0_u + l1_u * u_dep_u,
                    p3_u + l2_u * u_app_u,
                    p3_u,
                )
                bez_v = (
                    p0_v,
                    p0_v + l1_v * u_dep_v,
                    p3_v + l2_v * u_app_v,
                    p3_v,
                )
                if _bezier_has_inflection(bez_u) or _bezier_has_inflection(
                    bez_v
                ):
                  continue
                wings_u = arrowhead_wing_segments(*bez_u, electrons=1)
                wings_v = arrowhead_wing_segments(*bez_v, electrons=1)
                if _wing_bond_penalty(wings_u, wings_v, min_dist_cm=0.060) > 0:
                  continue
                if (
                    _wing_bond_penalty(wings_u, non_src_u, min_dist_cm=0.032)
                    > 0
                    or _wing_bond_penalty(wings_v, non_src_v, min_dist_cm=0.032)
                    > 0
                ):
                  continue
                poly_u_f = sample_bezier_pts(*bez_u, n=32, t_min=0.0, t_max=1.0)
                poly_v_f = sample_bezier_pts(*bez_v, n=32, t_min=0.0, t_max=1.0)
                pair_cm = _polyline_min_dist_cm(poly_u_f, poly_v_f)
                if pair_cm < 0.085:
                  continue
                if pair_cm < 0.25 and (
                    _wings_hit_polyline(wings_v, poly_u_f)
                    or _wings_hit_polyline(wings_u, poly_v_f)
                ):
                  continue
                if prior_polys:
                  cr_p1, sf_p1 = _prior_arrow_penalty(
                      poly_u_f, wings_u, prior_polys
                  )
                  cr_p2, sf_p2 = _prior_arrow_penalty(
                      poly_v_f, wings_v, prior_polys
                  )
                  if cr_p1 > 0 or cr_p2 > 0:
                    continue
                  prior_pen = sf_p1 + sf_p2
                else:
                  prior_pen = 0.0
                if hits_plus_sign(poly_u_f, plus_xs, wings_u) or hits_plus_sign(
                    poly_v_f, plus_xs, wings_v
                ):
                  continue
                if _tail_clearance_below(
                    poly_u_f, non_inc_u
                ) or _tail_clearance_below(poly_v_f, non_inc_v):
                  continue
                geom_score = (
                    arch_bonus
                    - pen_u
                    - pen_v
                    - 1.8 * h_idx
                    - 8.0 * abs(s_shift)
                    - 3.0 * abs(lift_s - 0.38)
                    - 0.02 * prior_pen
                )
                survivors.append((geom_score, bez_u, bez_v))

  if not survivors:
    return None
  survivors.sort(key=lambda item: item[0], reverse=True)
  best_pair = None
  best_total = -1e9
  for geom_score, bez_u, bez_v in survivors[:18]:
    aud = audit_arrow_clearances(
        [("u", bez_u), ("v", bez_v)], scene_canvases, plus_xs=plus_xs
    )
    if not aud["all_clear"]:
      continue
    tot = (
        geom_score
        + 10.0 * min(aud["min_bond_clearance_mm"], 1.2)
        + 6.0 * min(aud["min_glyph_clearance_mm"], 0.8)
        + 4.0 * min(aud["closest_arrow_pair_mm"], 2.0)
    )
    if tot > best_total:
      best_total = tot
      best_pair = (bez_u, bez_v)
  return best_pair


def compile_clean_reaction(
    rxn_def,
    out_dir: str | None = None,
    spec=None,
    auto_rotate: bool = False,
    strict_audit: bool = False,
):
  """Renders one reaction mechanism to a 300 DPI PNG.

  Args:
      rxn_def: Reaction definition with `id`, reactants/products configs
        (specified either as `reactants`/`products` lists/dicts or as `r1`,
        `r2`, ... / `p1`, `p2`, ...), and optional `cond` labels.
      out_dir: Directory to write `<id>.png` into.
      spec: Arrow specification. Defaults to the built-in catalog entry for
        `rxn_def["id"]`; pass explicitly to render a custom mechanism that is
        not in the catalog.
      auto_rotate: When True (or when `rxn_def["auto_rotate"]` /
        `spec["auto_rotate"]` is True), automatically solves ChemDraw/IUPAC
        reactant orientations and MCS product alignments for any molecule not
        explicitly listed in `spec["rot"]`.
      strict_audit: When True, raises `ValueError` if any curved arrow fails
        bond, atom-label, or arrow-to-arrow clearance thresholds.

  Returns:
      A result tuple on success, or a "[SKIP]" string if no spec is available.

  Raises:
      ValueError: If mechanism verification, molecule extraction, or strict
        clearance auditing fails.
  """
  ensure_dirs(out_dir=out_dir)
  rid = rxn_def["id"]
  if spec is None:
    spec = reaction_specs_100.SPECS.get(rid)
  if not spec:
    return f"[SKIP] No spec for {rid}"

  # Formally verify the elementary step with ElectronGraph before rendering
  audit_info = verify_mechanism_with_electron_graph(rxn_def, spec)
  r_kekule_map = audit_info.get("r_kekule_map", {}) if audit_info else {}
  p_kekule_map = audit_info.get("p_kekule_map", {}) if audit_info else {}
  p_ref_map = audit_info.get("p_ref_map", {}) if audit_info else {}

  explicit_rot_keys = set((spec.get("rot") or {}).keys())
  rot_map = dict(spec.get("rot", {}))
  if auto_rotate or rxn_def.get("auto_rotate") or spec.get("auto_rotate"):
    solved_rots = solve_reaction_rotations(rxn_def, spec)
    for mol_k, solved_deg in solved_rots.items():
      if mol_k not in rot_map:
        rot_map[mol_k] = solved_deg

  add_h_map = spec.get("add_h", {})
  coord_maps = spec.get("coord_map", {})
  cond = rxn_def.get("cond", {})

  reactant_items = _extract_molecules(rxn_def, "r")
  if not reactant_items:
    raise ValueError(f"{rid}: No reactants found in rxn_def")

  product_items = _extract_molecules(rxn_def, "p")
  if not product_items:
    raise ValueError(f"{rid}: No products found in rxn_def")

  for k, cfg in reactant_items:
    if "rot" in cfg:
      explicit_rot_keys.add(k)
    cfg["rot"] = rot_map.get(k, cfg.get("rot", 0.0))
    if k in add_h_map:
      cfg["add_h"] = add_h_map[k]

  for k, cfg in product_items:
    if "rot" in cfg:
      explicit_rot_keys.add(k)
    cfg["rot"] = rot_map.get(k, cfg.get("rot", 0.0))
    if k in add_h_map:
      cfg["add_h"] = add_h_map[k]

  # Layout Reactants
  plus_xs: list[float] = []
  reactant_canvases: list[tuple[str, dict[str, Any], MoleculeCanvas, float]] = (
      []
  )
  curr_right = 0.0

  def _build_reactant_canvas(
      mol_k: str,
      r_cfg: dict[str, Any],
      w_cm: float,
      cx: float,
      *,
      flip_y: bool = False,
  ) -> MoleculeCanvas:
    intramol_for_k = []
    btb_for_k = []
    for a_spec in spec.get("arrows", []):
      if a_spec.get("type") == "attack":
        src_spec = a_spec.get("src")
        tgt_spec = a_spec.get("tgt")
        if (
            isinstance(src_spec, (list, tuple))
            and isinstance(tgt_spec, (list, tuple))
            and len(src_spec) == 2
            and len(tgt_spec) == 2
            and src_spec[0] == mol_k
            and tgt_spec[0] == mol_k
        ):
          intramol_for_k.append((src_spec[1], tgt_spec[1], a_spec.get("conn")))
      elif a_spec.get("type") == "bond_to_bond":
        s_m = a_spec.get("src_mol") or a_spec.get("mol", "r1")
        t_m = a_spec.get("tgt_mol") or a_spec.get("mol", s_m)
        if s_m == mol_k and t_m == mol_k:
          b1_s = a_spec.get("b1") or a_spec.get("src_bond")
          b2_s = a_spec.get("b2") or a_spec.get("tgt_bond")
          if b1_s and b2_s:
            btb_for_k.append((b1_s, b2_s))
    return MoleculeCanvas(
        r_cfg["smi"],
        f"{rid}_{mol_k}",
        rot_deg=r_cfg["rot"],
        add_h=r_cfg.get("add_h"),
        coord_map=None if flip_y else coord_maps.get(mol_k),
        w_cm=w_cm,
        canvas_pos=(cx, 0.0),
        atom_labels=r_cfg.get("atom_labels"),
        intramol_targets=intramol_for_k or None,
        kekule_bonds=r_kekule_map.get(mol_k),
        bond_shift_pairs=btb_for_k or None,
        flip_y=flip_y,
    )

  for i, (k, r_cfg) in enumerate(reactant_items):
    default_w = 4.0 if i == 0 else 3.0
    r_w = r_cfg.get("w", default_w)
    if i == 0:
      center_x = 0.0
      curr_right = center_x + r_w / 2.0
    else:
      plus_x = curr_right + 0.6
      plus_xs.append(plus_x)
      center_x = plus_x + 0.6 + r_w / 2.0
      curr_right = center_x + r_w / 2.0

    canvas = _build_reactant_canvas(k, r_cfg, r_w, center_x)
    reactant_canvases.append((k, r_cfg, canvas, center_x))

  if len(reactant_canvases) == 2 and (
      auto_rotate or rxn_def.get("auto_rotate") or spec.get("auto_rotate")
  ):
    k1, _, c1, _ = reactant_canvases[0]
    k2, _, c2, _ = reactant_canvases[1]
    inter_pts = []
    for a_s in spec.get("arrows", []):
      at_c = a_s.get("type")
      if at_c == "attack":
        s_c, t_c = a_s.get("src"), a_s.get("tgt")
        if s_c and t_c and {s_c[0], t_c[0]} == {k1, k2}:
          tok1, tok2 = (s_c[1], t_c[1]) if s_c[0] == k1 else (t_c[1], s_c[1])
          inter_pts.append((
              resolve_arrow_endpoint(c1, tok1)[0],
              resolve_arrow_endpoint(c2, tok2)[0],
          ))
      elif at_c == "bond_to_atom" and {
          a_s.get("src_mol"),
          a_s.get("tgt_mol"),
      } == {k1, k2}:
        b_c, ta_c = a_s.get("bond"), a_s.get("tgt_atom")
        tok1, tok2 = (b_c, ta_c) if a_s.get("src_mol") == k1 else (ta_c, b_c)
        inter_pts.append((
            resolve_arrow_endpoint(c1, tok1)[0],
            resolve_arrow_endpoint(c2, tok2)[0],
        ))
      elif at_c == "bond_to_bond":
        s_mol = a_s.get("src_mol") or a_s.get("mol", "r1")
        t_mol = a_s.get("tgt_mol") or a_s.get("mol", s_mol)
        b1_c = a_s.get("b1") or a_s.get("src_bond")
        b2_c = a_s.get("b2") or a_s.get("tgt_bond")
        if {s_mol, t_mol} == {k1, k2} and b1_c and b2_c:
          s_cv, t_cv = (c1, c2) if s_mol == k1 else (c2, c1)
          try:
            pt_s = resolve_arrow_endpoint(s_cv, b1_c)[0]
            try:
              pt_t = resolve_arrow_endpoint(t_cv, b2_c)[0]
            except (ValueError, KeyError, TypeError):
              t_tok = b2_c[0] if b2_c[1] in b1_c else b2_c[1]
              pt_t = resolve_arrow_endpoint(t_cv, t_tok)[0]
            inter_pts.append((pt_s, pt_t) if s_mol == k1 else (pt_t, pt_s))
          except (ValueError, KeyError, TypeError):
            pass
    if len(inter_pts) >= 2:
      d1 = inter_pts[0][0] - inter_pts[1][0]
      d2 = inter_pts[0][1] - inter_pts[1][1]
      if (
          abs(float(d1[1])) > max(0.12, 1.5 * abs(float(d1[0])))
          and abs(float(d2[1])) > max(0.12, 1.5 * abs(float(d2[0])))
          and float(d1[1]) * float(d2[1]) < -0.02
      ):
        for flip_i in (1, 0):
          kf, rf_cfg, cf, cxf = reactant_canvases[flip_i]
          if (
              kf not in explicit_rot_keys
              and coord_maps.get(kf) is None
              and _is_planar_achiral_mol(cf.mol)
          ):
            w_f = rf_cfg.get("w", 4.0 if flip_i == 0 else 3.0)
            reactant_canvases[flip_i] = (
                kf,
                rf_cfg,
                _build_reactant_canvas(kf, rf_cfg, w_f, cxf, flip_y=True),
                cxf,
            )
            break

  if len(reactant_canvases) == 1:
    rxn_arrow_start_x = curr_right + 1.2
  else:
    rxn_arrow_start_x = curr_right + 0.8
  cond_text_len = max(len(cond.get("top", "")), len(cond.get("bottom", "")))
  arrow_len = max(2.8, min(6.8, cond_text_len * 0.145))
  rxn_arrow_end_x = rxn_arrow_start_x + arrow_len

  # Layout Products
  product_plus_xs = []
  product_canvases = []
  curr_prod_right = rxn_arrow_end_x + 0.8

  for j, (k, p_cfg) in enumerate(product_items):
    default_w = 4.0 if j == 0 else 2.2
    p_w = p_cfg.get("w", default_w)
    if j == 0:
      center_x = curr_prod_right + p_w / 2.0
      curr_prod_right = center_x + p_w / 2.0
    else:
      prod_plus_x = curr_prod_right + 0.5
      product_plus_xs.append(prod_plus_x)
      center_x = prod_plus_x + 0.5 + p_w / 2.0
      curr_prod_right = center_x + p_w / 2.0

    r_partner_key = p_ref_map.get(k)
    ref_c = None
    if r_partner_key and k not in explicit_rot_keys:
      ref_c = next(
          (rc for rk, _, rc, _ in reactant_canvases if rk == r_partner_key),
          None,
      )

    canvas = MoleculeCanvas(
        p_cfg["smi"],
        f"{rid}_{k}",
        rot_deg=p_cfg["rot"],
        add_h=p_cfg.get("add_h"),
        coord_map=coord_maps.get(k),
        w_cm=p_w,
        canvas_pos=(center_x, 0.0),
        atom_labels=p_cfg.get("atom_labels"),
        kekule_bonds=p_kekule_map.get(k),
        ref_canvas=ref_c,
        is_product=True,
    )
    product_canvases.append((k, p_cfg, canvas, center_x))

  # Build Lookup Maps
  mol_map = {}
  for k, _, canvas, _ in reactant_canvases:
    mol_map[k] = canvas
  for k, _, canvas, _ in product_canvases:
    mol_map[k] = canvas

  all_bonds = []
  for _, _, canvas, _ in reactant_canvases:
    for b in canvas.bonds:
      b_copy = dict(b)
      b_copy["mol"] = canvas
      all_bonds.append(b_copy)

  # Compute Explicit Bézier Curved Arrows
  tikz_arrows = []
  audit_arrows = []

  def _style_for_curve(pts_tuple, electrons_cnt):
    poly_c = _sample_bezier_pts_2d(*pts_tuple, n=40)
    arc_l = float(
        sum(
            np.linalg.norm(poly_c[idx_c + 1] - poly_c[idx_c])
            for idx_c in range(len(poly_c) - 1)
        )
    )
    if electrons_cnt == 1:
      ccw_turn = bezier_turns_ccw(*pts_tuple)
      if arc_l < 0.45:
        return "fishhookarrowshortswap" if ccw_turn else "fishhookarrowshort"
      return "fishhookarrowswap" if ccw_turn else "fishhookarrow"
    return "chemarrowshort" if arc_l < 0.48 else "chemarrow"

  def _arrow_label(num_idx: int, spec_item: dict[str, Any]) -> str:
    el_cnt = int(spec_item.get("electrons", 2))
    kind = spec_item["type"]
    return (
        f"arrow {num_idx} ({kind}, 1e)"
        if el_cnt == 1
        else f"arrow {num_idx} ({kind})"
    )

  scene_canvases = [c for c in mol_map.values() if c is not None]
  arrow_specs = list(spec["arrows"])
  for a_num, a_spec in enumerate(arrow_specs, 1):
    atype = a_spec["type"]
    if atype not in _ARROW_TYPES:
      raise ValueError(
          f"{rid}: {_arrow_label(a_num, a_spec)} has unknown type {atype!r}."
          f" Supported types are {', '.join(_ARROW_TYPES)}."
      )

  # Route tight local shifts first, then attacks ordered top-to-bottom, and
  # cleavages last so every arrow avoids previously routed trajectories while
  # emitting TikZ and audit entries in the spec's original 1..N order.
  route_order = sorted(
      range(len(arrow_specs)),
      key=lambda idx_s: _arrow_route_sort_key(arrow_specs[idx_s], mol_map),
  )
  routed_results: dict[
      int,
      tuple[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray], str],
  ] = {}
  prior_polys: list[np.ndarray] = []
  paired_1e = match_nascent_1e_pairs(arrow_specs, mol_map)

  def _record_routed(idx_r: int, bez_r, el_r: int) -> None:
    p0_r, p1_r, p2_r, p3_r = bez_r
    style_r = _style_for_curve(bez_r, el_r)
    cmd_r = (
        f"\\draw[{style_r}] ({p0_r[0]:.3f}, {p0_r[1]:.3f}) .. controls"
        f" ({p1_r[0]:.3f}, {p1_r[1]:.3f}) and ({p2_r[0]:.3f},"
        f" {p2_r[1]:.3f}) .. ({p3_r[0]:.3f}, {p3_r[1]:.3f});"
    )
    routed_results[idx_r] = (
        _arrow_label(idx_r + 1, arrow_specs[idx_r]),
        bez_r,
        cmd_r,
    )
    prior_polys.append(sample_bezier_pts(*bez_r, n=60, t_min=0.0, t_max=1.0))

  for spec_idx in route_order:
    if spec_idx in routed_results:
      continue
    if spec_idx in paired_1e:
      partner_idx, desc_self, desc_part = paired_1e[spec_idx]
      if partner_idx not in routed_results:
        pair_res = _compute_paired_fishhook_beziers(
            desc_self,
            desc_part,
            mol_map,
            scene_canvases,
            plus_xs,
            prior_polys,
        )
        if pair_res is not None:
          _record_routed(spec_idx, pair_res[0], 1)
          _record_routed(partner_idx, pair_res[1], 1)
          continue
    a_spec = arrow_specs[spec_idx]
    atype = a_spec["type"]
    electrons = int(a_spec.get("electrons", 2))
    res = None

    if atype == "attack":
      src_mol = mol_map[a_spec["src"][0]]
      tgt_mol = mol_map[a_spec["tgt"][0]]
      res = compute_attack_bezier(
          src_mol,
          a_spec["src"][1],
          tgt_mol,
          a_spec["tgt"][1],
          all_bonds=all_bonds,
          plus_x=plus_xs if plus_xs else None,
          side=a_spec.get("side"),
          lift=a_spec.get("lift"),
          canvases=scene_canvases,
          prior_arrows=prior_polys,
          electrons=electrons,
      )

    elif atype == "bond_cleave":
      mol_key = a_spec.get("mol", "r1")
      m = mol_map[mol_key]
      a1_idx = resolve_atom_token(m, a_spec["bond"][0])
      a2_idx = resolve_atom_token(m, a_spec["bond"][1])
      tgt_idx = resolve_atom_token(m, a_spec.get("tgt", a_spec["bond"][1]))
      c_idx_cleave = a1_idx if a1_idx != tgt_idx else a2_idx
      p_c_cleave = m.atoms[c_idx_cleave]["pos"]
      att_pos = None
      att_dir = None
      for other_idx, other_a in enumerate(arrow_specs):
        if other_idx == spec_idx:
          continue
        o_type = other_a.get("type")
        hits_cleave = False
        if o_type == "attack":
          tgt = other_a.get("tgt")
          if (
              isinstance(tgt, (list, tuple))
              and len(tgt) == 2
              and tgt[0] == mol_key
          ):
            tgt_tok = tgt[1]
            tgt_set = (
                {resolve_atom_token(m, sub_t) for sub_t in tgt_tok}
                if isinstance(tgt_tok, (list, tuple))
                else {resolve_atom_token(m, tgt_tok)}
            )
            if tgt_set & {a1_idx, a2_idx}:
              hits_cleave = True
              src_m = mol_map[other_a["src"][0]]
              att_pos = resolve_arrow_endpoint(src_m, other_a["src"][1])[0]
        elif o_type == "lone_pair_to_bond":
          o_mk = other_a.get("mol", "r1")
          try:
            ep0 = _resolve_mol_and_idx(o_mk, other_a["bond"][0], mol_map)
            ep1 = _resolve_mol_and_idx(o_mk, other_a["bond"][1], mol_map)
            if (mol_key, c_idx_cleave) in (ep0, ep1):
              hits_cleave = True
              ep_att = _resolve_mol_and_idx(o_mk, other_a["atom"], mol_map)
              att_pos = mol_map[ep_att[0]].atoms[ep_att[1]]["pos"]
          except (KeyError, ValueError, TypeError, IndexError):
            pass
        elif o_type == "bond_to_bond":
          b2_tok = other_a.get("b2") or other_a.get("tgt_bond")
          o_sm = other_a.get("src_mol") or other_a.get("mol", "r1")
          o_tm = other_a.get("tgt_mol") or other_a.get("mol", o_sm)
          if b2_tok:
            try:
              ep0 = _resolve_mol_and_idx(o_tm, b2_tok[0], mol_map)
              ep1 = _resolve_mol_and_idx(o_tm, b2_tok[1], mol_map)
              if (mol_key, c_idx_cleave) in (ep0, ep1):
                hits_cleave = True
                b1_tok = other_a.get("b1") or other_a.get("src_bond")
                if b1_tok and o_sm in mol_map:
                  att_pos = resolve_arrow_endpoint(mol_map[o_sm], b1_tok)[0]
            except (KeyError, ValueError, TypeError, IndexError):
              pass
        elif o_type == "bond_to_atom" and other_a.get("tgt_mol") == mol_key:
          if resolve_atom_token(m, other_a["tgt_atom"]) in {a1_idx, a2_idx}:
            hits_cleave = True
        if hits_cleave:
          if other_idx in routed_results:
            att_poly = sample_bezier_pts(*routed_results[other_idx][1], n=30)
            att_dir = att_poly[-6] - p_c_cleave
          break
      res = compute_bond_cleave_bezier(
          m,
          a1_idx,
          a2_idx,
          tgt_idx,
          incoming_from_right=True,
          attack_pos=att_pos,
          side=a_spec.get("side"),
          attack_dir=att_dir,
          canvases=scene_canvases,
          prior_arrows=prior_polys,
          electrons=electrons,
      )

    elif atype == "bond_to_bond":
      src_mol_k = a_spec.get("src_mol") or a_spec.get("mol", "r1")
      tgt_mol_k = a_spec.get("tgt_mol") or a_spec.get("mol", src_mol_k)
      m = mol_map[src_mol_k]
      tgt_m = mol_map[tgt_mol_k]
      b1_tok = a_spec.get("b1") or a_spec["src_bond"]
      b2_tok = a_spec.get("b2") or a_spec["tgt_bond"]
      res = compute_bond_to_bond_bezier(
          m,
          b1_tok,
          b2_tok,
          canvases=scene_canvases,
          allow_incipient=bool(a_spec.get("incipient", True)),
          tgt_mol=tgt_m,
          all_bonds=all_bonds,
          plus_x=plus_xs if plus_xs else None,
          side=a_spec.get("side"),
          lift=a_spec.get("lift"),
          prior_arrows=prior_polys,
          electrons=electrons,
      )

    elif atype == "bond_to_atom":
      res = compute_bond_to_atom_bezier(
          mol_map[a_spec.get("src_mol", "r1")],
          a_spec["bond"],
          mol_map[a_spec.get("tgt_mol", "r2")],
          a_spec["tgt_atom"],
          lift=a_spec.get("lift"),
          canvases=scene_canvases,
          prior_arrows=prior_polys,
          electrons=electrons,
          plus_x=plus_xs if plus_xs else None,
      )

    elif atype == "lone_pair_to_bond":
      m = mol_map[a_spec.get("mol", "r1")]
      res = compute_lone_pair_to_bond_bezier(
          m,
          a_spec["atom"],
          a_spec["bond"],
          side=a_spec.get("side"),
          canvases=scene_canvases,
          allow_incipient=bool(a_spec.get("incipient", electrons == 1)),
          prior_arrows=prior_polys,
          electrons=electrons,
          all_bonds=all_bonds,
          plus_x=plus_xs if plus_xs else None,
      )

    if res is not None:
      _record_routed(spec_idx, res, electrons)

  for spec_idx in range(len(arrow_specs)):
    if spec_idx in routed_results:
      a_label, res, tikz_cmd = routed_results[spec_idx]
      audit_arrows.append((a_label, res))
      tikz_arrows.append(tikz_cmd)

  audit = audit_arrow_clearances(
      audit_arrows, scene_canvases, plus_xs=[*plus_xs, *product_plus_xs]
  )
  report_arrow_audit(rid, audit, len(spec["arrows"]), strict=strict_audit)

  # Build Standalone TikZ
  tex_lines = [
      r"\documentclass[border=15pt]{standalone}",
      r"\usepackage{amsmath}",
      r"\usepackage{amssymb}",
      r"\usepackage{graphicx}",
      r"\usepackage{tikz}",
      r"\usetikzlibrary{arrows.meta,bending,calc,positioning}",
      r"\begin{document}",
      r"\begin{tikzpicture}[",
      r"    font=\sffamily,",
      r"    line cap=round,",
      r"    line join=round,",
      r"    chemarrow/.style={",
      r"      -{Stealth[length=2.0mm, width=1.30mm, flex=0.75]},",
      r"      line width=1.0pt,",
      r"      color=red!85!black",
      r"    },",
      r"    chemarrowshort/.style={",
      r"      -{Stealth[length=1.35mm, width=0.98mm, flex=0.75]},",
      r"      line width=0.92pt,",
      r"      color=red!85!black",
      r"    },",
      r"    fishhookarrow/.style={",
      r"      -{Stealth[harpoon, length=2.45mm, width=2.25mm, flex=0.75]},",
      r"      line width=1.0pt,",
      r"      color=red!85!black",
      r"    },",
      r"    fishhookarrowswap/.style={",
      (
          r"      -{Stealth[harpoon, swap, length=2.45mm, width=2.25mm,"
          r" flex=0.75]},"
      ),
      r"      line width=1.0pt,",
      r"      color=red!85!black",
      r"    },",
      r"    fishhookarrowshort/.style={",
      r"      -{Stealth[harpoon, length=1.85mm, width=1.75mm, flex=0.75]},",
      r"      line width=0.95pt,",
      r"      color=red!85!black",
      r"    },",
      r"    fishhookarrowshortswap/.style={",
      (
          r"      -{Stealth[harpoon, swap, length=1.85mm, width=1.75mm,"
          r" flex=0.75]},"
      ),
      r"      line width=0.95pt,",
      r"      color=red!85!black",
      r"    }",
      r"  ]",
      r"",
  ]

  def _caption_below_str(canvas_c, cx_c: float, raw_label: str) -> str:
    y_south = float(canvas_c.canvas_pos[1]) - 0.5 * float(canvas_c.h_cm)
    half_w = (
        0.5
        * max(
            float(canvas_c.w_cm),
            min(float(canvas_c.w_cm) + 2.4, len(str(raw_label)) * 0.19),
        )
        + 0.15
    )
    min_ay = y_south
    for poly_p in prior_polys:
      in_col = np.abs(poly_p[:, 0] - cx_c) <= half_w
      if np.any(in_col):
        min_ay = min(min_ay, float(np.min(poly_p[in_col, 1])))
    below_pt = 8.0 + max(0.0, y_south - min_ay) * (72.27 / 2.54)
    return "8pt" if below_pt <= 8.05 else f"{below_pt:.1f}pt"

  for i, (k, r_cfg, canvas, center_x) in enumerate(reactant_canvases):
    if i > 0:
      px = plus_xs[i - 1]
      tex_lines.append(
          f"  \\node[fill=white, inner sep=2pt] at ({px:.2f}, 0) {{\\huge"
          " $+$};"
      )
    tex_lines.append(
        f"  \\node[inner sep=0pt] ({k}) at ({center_x:.2f}, 0)"
        f" {{\\includegraphics[width={canvas.w_cm:.2f}cm]{{{canvas.png_path}}}}};"
    )
    if r_cfg.get("name"):
      r_name = latex_escape_label(r_cfg["name"])
      b_opt = _caption_below_str(canvas, center_x, r_cfg["name"])
      tex_lines.append(
          f"  \\node[below={b_opt} of {k}, align=center,"
          f" font=\\small\\bfseries] {{{r_name}}};"
      )
    tex_lines.append(r"")

  top_c = latex_escape_label(cond.get("top", ""))
  bot_c = latex_escape_label(cond.get("bottom", ""))
  tex_lines.extend([
      (
          f"  \\draw[->, line width=1.3pt, >=Stealth] ({rxn_arrow_start_x:.2f},"
          f" 0) -- ({rxn_arrow_end_x:.2f}, 0)"
      ),
      "    node[midway, above=2pt, align=center, font=\\footnotesize]"
      f" {{{top_c}}}"
      if top_c
      else "",
      "    node[midway, below=2pt, align=center, font=\\footnotesize]"
      f" {{{bot_c}}};"
      if bot_c
      else "    ;",
      r"",
  ])

  for j, (k, p_cfg, canvas, center_x) in enumerate(product_canvases):
    if j > 0:
      ppx = product_plus_xs[j - 1]
      tex_lines.append(
          f"  \\node[fill=white, inner sep=2pt] at ({ppx:.2f}, 0) {{\\huge"
          r" $+$};"
      )
    tex_lines.append(
        f"  \\node[inner sep=0pt] ({k}) at ({center_x:.2f}, 0)"
        f" {{\\includegraphics[width={canvas.w_cm:.2f}cm]{{{canvas.png_path}}}}};"
    )
    if p_cfg.get("name"):
      p_name = latex_escape_label(p_cfg["name"])
      b_opt_p = _caption_below_str(canvas, center_x, p_cfg["name"])
      tex_lines.append(
          f"  \\node[below={b_opt_p} of {k}, align=center,"
          f" font=\\small\\bfseries] {{{p_name}}};"
      )
    tex_lines.append(r"")

  tex_lines.append("  % Curved Arrows")
  for arr in tikz_arrows:
    tex_lines.append(f"  {arr}")

  tex_lines.extend([r"\end{tikzpicture}", r"\end{document}", r""])

  tex_code = "\n".join(tex_lines)
  tex_file = os.path.join(TEX_DIR, f"{rid}.tex")
  png_final = os.path.join(IMG_DIR, f"{rid}.png")

  _atomic_write_bytes(tex_file, tex_code.encode("utf-8"))

  rendered_png = False
  has_pdflatex = shutil.which("pdflatex") is not None
  if has_pdflatex:
    with tempfile.TemporaryDirectory(
        dir=TEX_DIR, prefix=f"job_{rid}_"
    ) as job_dir:
      job_tex = os.path.join(job_dir, f"{rid}.tex")
      job_pdf = os.path.join(job_dir, f"{rid}.pdf")
      job_png_prefix = os.path.join(job_dir, rid)
      job_png = f"{job_png_prefix}.png"
      with open(job_tex, "w", encoding="utf-8") as jf:
        jf.write(tex_code)
      subprocess.run(
          ["pdflatex", "-interaction=nonstopmode", f"{rid}.tex"],
          cwd=job_dir,
          stdout=subprocess.DEVNULL,
          stderr=subprocess.DEVNULL,
          check=False,
      )
      if os.path.exists(job_pdf):
        subprocess.run(
            [
                "pdftoppm",
                "-png",
                "-r",
                "300",
                "-singlefile",
                job_pdf,
                job_png_prefix,
            ],
            cwd=job_dir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if os.path.exists(job_png):
          with open(job_png, "rb") as png_in:
            _atomic_write_bytes(png_final, png_in.read())
          rendered_png = True

  if not rendered_png:
    # Fallback to high-resolution Matplotlib curved-arrow renderer
    _render_mechanism_matplotlib(
        rid=rid,
        title=spec.get("title", "") if spec else "",
        reactant_canvases=reactant_canvases,
        product_canvases=product_canvases,
        plus_xs=plus_xs,
        product_plus_xs=product_plus_xs,
        rxn_arrow_start_x=rxn_arrow_start_x,
        rxn_arrow_end_x=rxn_arrow_end_x,
        cond=cond,
        audit_arrows=audit_arrows,
        png_final=png_final,
    )

  if os.path.exists(png_final):
    return (rid, tikz_arrows)
  else:
    return f"[ERROR] Rendering failed for {rid}"


def run_all_parallel(cores: int | None = None, out_dir: str | None = None):
  t0 = time.time()
  ensure_dirs(out_dir=out_dir)

  # Remove old cached mol pngs to force fresh rotated rendering
  for p in glob.glob(os.path.join(CACHE_DIR, "*.png")):
    os.remove(p)

  num_cores = cores or mp.cpu_count()
  print(
      f"Compiling all {len(reactions_data.REACTIONS)} reactions using"
      f" {num_cores} CPU cores..."
  )

  with mp.Pool(num_cores) as pool:
    results = pool.map(compile_clean_reaction, reactions_data.REACTIONS)

  dt = time.time() - t0
  success = [r for r in results if isinstance(r, tuple)]
  errors = [
      r for r in results if isinstance(r, str) and r.startswith("[ERROR]")
  ]
  print(
      f"Done in {dt:.2f}s! Successfully compiled"
      f" {len(success)}/{len(reactions_data.REACTIONS)} reactions."
  )
  if errors:
    print("Errors encountered:", file=sys.stderr)
    for e in errors:
      print(" ", e)

  # Update reactions_data.py with generated arrows
  arrow_dict = dict(success)
  for r in reactions_data.REACTIONS:
    rid = r["id"]
    if rid in arrow_dict:
      r["arrows_tikz"] = arrow_dict[rid]

  # Write back reactions_data.py
  with open(os.path.join(BASE_DIR, "reactions_data.py"), "w") as f:
    f.write(
        "# Auto-generated reactions database with verified cubic Bezier"
        " ChemDraw arrows\n\n"
    )
    f.write(f"REACTIONS = {repr(reactions_data.REACTIONS)}\n")
  print("Updated reactions_data.py successfully!")


def compile_multistep_cascade(
    cascade_def: dict,
    out_dir: str | None = None,
    auto_rotate: bool = False,
    strict_audit: bool = False,
):
  """Verifies and renders a multi-step reaction cascade to a composite PNG.

  Verifies each elementary step with ElectronGraph, enforces topological
  intermediate continuity between consecutive steps (at least one product of
  step k must appear in step k+1 reactants unless verify_continuity is False),
  renders each step, and composes a vertical multi-step scheme diagram.

  Args:
      cascade_def: Multi-step cascade dictionary containing `'id'` and
        `'steps'`.
      out_dir: Optional directory for rendered PNG output.
      auto_rotate: If True, optimizes 2D reactant orientations per step.
      strict_audit: If True, raises ValueError on any arrow clearance violation.

  Returns:
      Tuple `(cid, step_results, final_png)` (or `(cid, step_results, None)` if
      `pdflatex` is unavailable).

  Raises:
      ValueError: If `'steps'` is empty, intermediate continuity fails, or any
        step fails verification or strict clearance auditing.
  """
  ensure_dirs(out_dir=out_dir)
  cid = cascade_def["id"]
  steps = cascade_def.get("steps", [])
  if not isinstance(steps, list) or not steps:
    raise ValueError(f"{cid}: 'steps' must be a non-empty list of step dicts.")

  verify_continuity = cascade_def.get("verify_continuity", True)
  cascade_auto_rot = bool(auto_rotate or cascade_def.get("auto_rotate"))

  step_results = []
  step_panels = []
  prev_prod_canons: set[str] = set()

  for idx, step in enumerate(steps, 1):
    step_id = step.get("id", f"{cid}_step_{idx}")
    step_rxn = dict(step)
    step_rxn["id"] = step_id
    step_spec = step_rxn.get("spec")
    if step_spec is None:
      step_spec = {
          k: step_rxn[k]
          for k in ("arrows", "rot", "add_h", "coord_map", "auto_rotate")
          if k in step_rxn
      }

    if idx > 1 and verify_continuity and prev_prod_canons:
      curr_r_items = _extract_molecules(step_rxn, "r")
      curr_r_canons = {
          _canon_smi_no_iso(r["smi"]) for _, r in curr_r_items if r.get("smi")
      }
      if not (curr_r_canons & prev_prod_canons):
        raise ValueError(
            f"[CASCADE ERROR] Disconnected cascade between step {idx - 1} and"
            f" step {idx} in '{cid}': no intermediate from step {idx - 1}"
            f" products ({sorted(prev_prod_canons)}) appears in step {idx}"
            f" reactants ({sorted(curr_r_canons)})."
        )

    res = compile_clean_reaction(
        step_rxn,
        out_dir=out_dir,
        spec=step_spec,
        auto_rotate=cascade_auto_rot,
        strict_audit=strict_audit,
    )
    if not isinstance(res, tuple):
      raise RuntimeError(f"Failed to render step {idx} ({step_id}): {res}")

    step_results.append(res)
    prod_items = _extract_molecules(step_rxn, "p")
    prev_prod_canons = {
        _canon_smi_no_iso(p["smi"]) for _, p in prod_items if p.get("smi")
    }

    step_png = os.path.join(IMG_DIR, f"{step_id}.png")
    step_title = step.get("title") or step.get("name") or f"Step {idx}"
    step_panels.append((step_title, step_png))

  loaded_imgs = []
  for step_title, step_png in step_panels:
    im = Image.open(step_png).convert("RGB")
    loaded_imgs.append((step_title, im))

  pad_x = 48
  pad_y = 36
  header_h = 54
  max_w = max(im.width for _, im in loaded_imgs) + 2 * pad_x
  total_h = sum(im.height + header_h + pad_y for _, im in loaded_imgs) + pad_y

  canvas_im = Image.new("RGB", (max_w, total_h), (255, 255, 255))
  draw = ImageDraw.Draw(canvas_im)
  try:
    font = ImageFont.truetype(
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 26
    )
  except Exception:
    font = ImageFont.load_default()

  curr_y = pad_y
  for idx, (step_title, im) in enumerate(loaded_imgs, 1):
    banner_label = (
        step_title
        if step_title.lower().startswith("step")
        else f"Step {idx}: {step_title}"
    )
    draw.text((pad_x, curr_y), banner_label, fill=(30, 41, 59), font=font)
    curr_y += header_h
    x_off = (max_w - im.width) // 2
    canvas_im.paste(im, (x_off, curr_y))
    curr_y += im.height + pad_y
    if idx < len(loaded_imgs):
      draw.line(
          [(pad_x, curr_y - pad_y // 2), (max_w - pad_x, curr_y - pad_y // 2)],
          fill=(226, 232, 240),
          width=2,
      )

  final_png = os.path.join(IMG_DIR, f"{cid}.png")
  canvas_im.save(final_png, dpi=(300, 300))
  return (cid, step_results, final_png)


def compile_catalytic_cycle(
    cycle_def: dict,
    out_dir: str | None = None,
    auto_rotate: bool = False,
):
  """Renders a circular catalytic cycle diagram to a 300 DPI PNG via TikZ.

  Supports N >= 2 catalytic intermediates placed around a circle with curved
  main-cycle arcs, transition labels, substrate entry arrows, byproduct exit
  arrows, and optional ElectronGraph verification of attached elementary steps.
  """
  ensure_dirs(out_dir=out_dir)
  cid = cycle_def["id"]
  cycle_cfg = (
      cycle_def["catalytic_cycle"]
      if isinstance(cycle_def.get("catalytic_cycle"), dict)
      else cycle_def
  )
  intermediates = cycle_cfg.get("intermediates", [])
  if not isinstance(intermediates, list) or len(intermediates) < 2:
    raise ValueError(
        f"{cid}: catalytic cycle requires at least 2 entries in"
        " 'intermediates'."
    )

  transitions = cycle_cfg.get("transitions", cycle_cfg.get("steps", []))
  for tr_idx, tr in enumerate(transitions):
    rxn_step = tr.get("rxn_step")
    if isinstance(rxn_step, dict) and (
        rxn_step.get("arrows") or (rxn_step.get("spec") or {}).get("arrows")
    ):
      tr_rxn = dict(rxn_step)
      tr_rxn.setdefault("id", f"{cid}_tr_{tr_idx + 1}")
      tr_spec = tr_rxn.get("spec") or {
          k: tr_rxn[k]
          for k in ("arrows", "rot", "add_h", "coord_map")
          if k in tr_rxn
      }
      verify_mechanism_with_electron_graph(tr_rxn, tr_spec)

  n_inter = len(intermediates)
  radius_cm = max(float(cycle_cfg.get("radius_cm", 0.0)), 1.25 * n_inter, 5.2)

  inter_nodes = []
  id_to_idx = {}
  for i, inter in enumerate(intermediates):
    iid = inter.get("id", f"I{i + 1}")
    id_to_idx[iid] = i
    theta_deg = 90.0 - i * (360.0 / n_inter)
    theta_rad = math.radians(theta_deg)
    cx = radius_cm * math.cos(theta_rad)
    cy = radius_cm * math.sin(theta_rad)
    w_cm = float(inter.get("w", 3.0))
    canvas = MoleculeCanvas(
        inter["smi"],
        f"{cid}_{iid}",
        rot_deg=float(inter.get("rot", 0.0)),
        add_h=inter.get("add_h"),
        w_cm=w_cm,
        canvas_pos=(cx, cy),
        atom_labels=inter.get("atom_labels"),
    )
    inter_nodes.append((iid, inter, canvas, cx, cy, theta_deg))

  tex_lines = [
      r"\documentclass[border=18pt]{standalone}",
      r"\usepackage{amsmath}",
      r"\usepackage{amssymb}",
      r"\usepackage{graphicx}",
      r"\usepackage{tikz}",
      r"\usetikzlibrary{arrows.meta,calc,positioning,bending}",
      r"\begin{document}",
      r"\begin{tikzpicture}[",
      r"    font=\sffamily,",
      r"    line cap=round,",
      r"    line join=round,",
      r"    cyclearrow/.style={",
      r"      -{Stealth[length=2.6mm, width=1.8mm, bend]},",
      r"      line width=1.3pt,",
      r"      color=black!85",
      r"    },",
      r"    sidearrow/.style={",
      r"      -{Stealth[length=2.0mm, width=1.4mm, bend]},",
      r"      line width=0.95pt,",
      r"      color=blue!72!black",
      r"    }",
      r"  ]",
      r"",
  ]

  title = latex_escape_label(
      cycle_cfg.get("title")
      or cycle_def.get("title")
      or cycle_def.get("name")
      or ""
  )
  if title:
    tex_lines.append(
        "  \\node[align=center, text width=6.8cm, font=\\large\\bfseries] at"
        f" (0, 0) {{{title}}};"
    )

  step_deg = 360.0 / n_inter
  pad_deg = min(0.27 * step_deg, max(12.0, math.degrees(1.55 / radius_cm)))
  for i in range(n_inter):
    j = (i + 1) % n_inter
    iid_from = inter_nodes[i][0]
    iid_to = inter_nodes[j][0]
    tr = next(
        (
            t
            for t in transitions
            if (t.get("from") == iid_from and t.get("to") == iid_to)
        ),
        transitions[i] if i < len(transitions) else {},
    )

    theta_start = 90.0 - i * step_deg - pad_deg
    theta_end = 90.0 - (i + 1) * step_deg + pad_deg
    theta_mid = 0.5 * (theta_start + theta_end)

    rad_s = math.radians(theta_start)
    rad_m = math.radians(theta_mid)

    sx = radius_cm * math.cos(rad_s)
    sy = radius_cm * math.sin(rad_s)
    tex_lines.append(
        f"  \\draw[cyclearrow] ({sx:.2f}, {sy:.2f}) arc"
        f" ({theta_start:.1f}:{theta_end:.1f}:{radius_cm:.2f}cm);"
    )

    tr_label = latex_escape_label(tr.get("label") or tr.get("cond") or "")
    in_lbl = latex_escape_label(
        tr.get("in_label")
        or tr.get("substrate_in")
        or tr.get("substrate")
        or ""
    )
    out_lbl = latex_escape_label(
        tr.get("out_label")
        or tr.get("byproduct_out")
        or tr.get("byproduct")
        or ""
    )

    # Polar unit vectors at transition midpoint: u_r (radial outward), u_t (clockwise tangent)
    ur_x, ur_y = math.cos(rad_m), math.sin(rad_m)
    ut_x, ut_y = math.sin(rad_m), -math.cos(rad_m)

    if tr_label:
      lbl_r = max(2.2, radius_cm - 2.15)
      lx = lbl_r * ur_x
      ly = lbl_r * ur_y
      tex_lines.append(
          "  \\node[align=center, text width=2.8cm,"
          " font=\\footnotesize\\bfseries, fill=white, inner sep=1.5pt] at"
          f" ({lx:.2f}, {ly:.2f}) {{{tr_label}}};"
      )

    if in_lbl:
      t_kiss = -0.16 if out_lbl else 0.0
      t_tail = -1.15 if out_lbl else -0.85
      r_tail = 1.45
      ix_start = (radius_cm + r_tail) * ur_x + t_tail * ut_x
      iy_start = (radius_cm + r_tail) * ur_y + t_tail * ut_y
      ix_end = (radius_cm + 0.08) * ur_x + t_kiss * ut_x
      iy_end = (radius_cm + 0.08) * ur_y + t_kiss * ut_y
      ic1_x = ix_start + 0.50 * ut_x - 0.45 * ur_x
      ic1_y = iy_start + 0.50 * ut_y - 0.45 * ur_y
      ic2_x = ix_end - 0.55 * ut_x + 0.12 * ur_x
      ic2_y = iy_end - 0.55 * ut_y + 0.12 * ur_y
      tex_lines.append(
          f"  \\draw[sidearrow] ({ix_start:.2f}, {iy_start:.2f}) .. controls"
          f" ({ic1_x:.2f}, {ic1_y:.2f}) and ({ic2_x:.2f}, {ic2_y:.2f}) .."
          f" ({ix_end:.2f}, {iy_end:.2f});"
      )
      lx_in = (radius_cm + 2.15) * ur_x + (t_tail - 0.20) * ut_x
      ly_in = (radius_cm + 2.15) * ur_y + (t_tail - 0.20) * ut_y
      tex_lines.append(
          "  \\node[align=center, text width=3.2cm, font=\\scriptsize,"
          f" text=blue!75!black, fill=white, inner sep=1pt] at ({lx_in:.2f},"
          f" {ly_in:.2f}) {{{in_lbl}}};"
      )

    if out_lbl:
      t_kiss_o = 0.16 if in_lbl else 0.0
      t_tip = 1.15 if in_lbl else 0.85
      r_tip = 1.45
      ox_start = (radius_cm + 0.08) * ur_x + t_kiss_o * ut_x
      oy_start = (radius_cm + 0.08) * ur_y + t_kiss_o * ut_y
      ox_end = (radius_cm + r_tip) * ur_x + t_tip * ut_x
      oy_end = (radius_cm + r_tip) * ur_y + t_tip * ut_y
      oc1_x = ox_start + 0.55 * ut_x + 0.12 * ur_x
      oc1_y = oy_start + 0.55 * ut_y + 0.12 * ur_y
      oc2_x = ox_end - 0.50 * ut_x - 0.45 * ur_x
      oc2_y = oy_end - 0.50 * ut_y - 0.45 * ur_y
      tex_lines.append(
          f"  \\draw[sidearrow] ({ox_start:.2f}, {oy_start:.2f}) .. controls"
          f" ({oc1_x:.2f}, {oc1_y:.2f}) and ({oc2_x:.2f}, {oc2_y:.2f}) .."
          f" ({ox_end:.2f}, {oy_end:.2f});"
      )
      lx_out = (radius_cm + 2.15) * ur_x + (t_tip + 0.20) * ut_x
      ly_out = (radius_cm + 2.15) * ur_y + (t_tip + 0.20) * ut_y
      tex_lines.append(
          "  \\node[align=center, text width=3.2cm, font=\\scriptsize,"
          f" text=blue!75!black, fill=white, inner sep=1pt] at ({lx_out:.2f},"
          f" {ly_out:.2f}) {{{out_lbl}}};"
      )

  for iid, inter, canvas, cx, cy, theta_deg in inter_nodes:
    lbl = latex_escape_label(inter.get("label") or inter.get("name") or "")
    if lbl:
      tex_lines.append(
          "  \\node[inner sep=2.5pt, fill=white, rounded corners=4pt,"
          f" align=center, font=\\footnotesize\\bfseries] ({iid}) at"
          f" ({cx:.2f}, {cy:.2f})"
          f" {{\\includegraphics[width={canvas.w_cm:.2f}cm]{{{canvas.png_path}}}\\\\[2pt]"
          f"{lbl}}};"
      )
    else:
      tex_lines.append(
          f"  \\node[inner sep=2.5pt, fill=white, rounded corners=4pt] ({iid})"
          f" at ({cx:.2f}, {cy:.2f})"
          f" {{\\includegraphics[width={canvas.w_cm:.2f}cm]{{{canvas.png_path}}}}};"
      )

  tex_lines.extend([r"\end{tikzpicture}", r"\end{document}", r""])
  tex_code = "\n".join(tex_lines)
  tex_file = os.path.join(TEX_DIR, f"{cid}.tex")
  pdf_file = os.path.join(TEX_DIR, f"{cid}.pdf")
  png_final = os.path.join(IMG_DIR, f"{cid}.png")

  if os.path.exists(pdf_file):
    os.remove(pdf_file)
  with open(tex_file, "w") as f:
    f.write(tex_code)

  proc = subprocess.run(
      ["pdflatex", "-interaction=nonstopmode", f"{cid}.tex"],
      cwd=TEX_DIR,
      stdout=subprocess.DEVNULL,
      stderr=subprocess.DEVNULL,
  )
  if proc.returncode != 0 or not os.path.exists(pdf_file):
    return f"[ERROR] Failed pdflatex for catalytic cycle {cid}"

  subprocess.run(
      [
          "pdftoppm",
          "-png",
          "-r",
          "300",
          "-singlefile",
          pdf_file,
          os.path.join(IMG_DIR, cid),
      ],
      stdout=subprocess.DEVNULL,
      stderr=subprocess.DEVNULL,
  )
  return (cid, png_final)


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Render chemical reaction mechanisms with curved electron-pushing"
          " arrows."
      )
  )
  parser.add_argument(
      "--id", type=str, help="Reaction ID to render (e.g. rxn_001)"
  )
  parser.add_argument(
      "--all", action="store_true", help="Render all 100 reactions in parallel"
  )
  parser.add_argument(
      "--out-dir",
      "--out_dir",
      type=str,
      default=None,
      help="Output directory for generated reaction PNG images",
  )
  parser.add_argument(
      "--cores",
      type=int,
      default=None,
      help="Number of CPU cores for parallel rendering",
  )
  parser.add_argument(
      "--base-dir",
      "--base_dir",
      type=str,
      default=None,
      help="Base working directory",
  )
  parser.add_argument(
      "--rxn-json",
      "--rxn_json",
      dest="rxn_json",
      type=str,
      default=None,
      help="Path to a custom mechanism specification in JSON.",
  )
  parser.add_argument(
      "--rxn-yaml",
      "--rxn_yaml",
      dest="rxn_yaml",
      type=str,
      default=None,
      help="Path to a custom mechanism specification in YAML.",
  )
  parser.add_argument(
      "--auto-rotate",
      "--auto_rotate",
      dest="auto_rotate",
      action="store_true",
      help=(
          "Automatically compute ChemDraw/IUPAC reactant orientations and MCS"
          " product alignments via rotation_solver."
      ),
  )
  args = parser.parse_args()

  if not args.id and not args.all and not args.rxn_json and not args.rxn_yaml:
    parser.error(
        "Must specify one of --id <reaction_id>, --all, --rxn-json or"
        " --rxn-yaml."
    )

  ensure_dirs(base_dir=args.base_dir, out_dir=args.out_dir)

  def _run_cli_step(fn, *f_args, **f_kwargs):
    try:
      return fn(*f_args, **f_kwargs)
    except (ValueError, KeyError, TypeError, RuntimeError, OSError) as exc:
      print(f"[ERROR] {type(exc).__name__}: {exc}", file=sys.stderr)
      sys.exit(1)

  if args.rxn_json or args.rxn_yaml:
    path = args.rxn_json or args.rxn_yaml
    if not os.path.exists(path):
      print(f"[ERROR] Specification file not found: {path}", file=sys.stderr)
      sys.exit(1)
    with open(path, encoding="utf-8") as fh:
      if args.rxn_yaml:
        rxn = yaml.safe_load(fh)
      else:
        rxn = json.load(fh)
    if not isinstance(rxn, dict) or "id" not in rxn:
      print(
          f"[ERROR] {path} must be a mapping with at least an 'id' key.",
          file=sys.stderr,
      )
      sys.exit(1)

    # Dispatch catalytic cycle vs multi-step cascade vs single elementary step
    if (
        rxn.get("mode") == "catalytic_cycle"
        or "catalytic_cycle" in rxn
        or "intermediates" in rxn
    ):
      res = _run_cli_step(
          compile_catalytic_cycle,
          rxn,
          out_dir=args.out_dir,
          auto_rotate=args.auto_rotate,
      )
      if isinstance(res, tuple):
        print(f"[SUCCESS] Rendered catalytic cycle {rxn['id']} -> {res[1]}")
      else:
        print(f"[ERROR] {res}", file=sys.stderr)
        sys.exit(1)
      return

    if isinstance(rxn.get("steps"), list):
      res = _run_cli_step(
          compile_multistep_cascade,
          rxn,
          out_dir=args.out_dir,
          auto_rotate=args.auto_rotate,
          strict_audit=True,
      )
      if isinstance(res, tuple):
        print(f"[SUCCESS] Rendered multi-step cascade {rxn['id']} -> {res[2]}")
      else:
        print(f"[ERROR] {res}", file=sys.stderr)
        sys.exit(1)
      return

    # Single elementary step: the arrow spec may be nested under "spec" or
    # written inline at the top level alongside the structures.
    spec = rxn.get("spec")
    if spec is None:
      spec = {
          k: rxn[k]
          for k in ("arrows", "rot", "add_h", "coord_map", "auto_rotate")
          if k in rxn
      }
    if not spec.get("arrows"):
      print(
          f"[ERROR] {path} defines no 'arrows'. A mechanism needs at "
          "least one curved arrow.",
          file=sys.stderr,
      )
      sys.exit(1)
    res = _run_cli_step(
        compile_clean_reaction,
        rxn,
        out_dir=args.out_dir,
        spec=spec,
        auto_rotate=args.auto_rotate,
        strict_audit=True,
    )
    if isinstance(res, tuple):
      print(
          f"[SUCCESS] Rendered {rxn['id']} ->"
          f" {os.path.join(args.out_dir or IMG_DIR, rxn['id'] + '.png')}"
      )
    else:
      print(f"[ERROR] {res}", file=sys.stderr)
      sys.exit(1)
  elif args.id:
    rxn = next(
        (r for r in reactions_data.REACTIONS if r["id"] == args.id), None
    )
    if not rxn:
      print(
          f"[ERROR] Reaction {args.id} not found in catalog.", file=sys.stderr
      )
      sys.exit(1)
    res = _run_cli_step(
        compile_clean_reaction,
        rxn,
        out_dir=args.out_dir,
        auto_rotate=args.auto_rotate,
        strict_audit=True,
    )
    if isinstance(res, tuple):
      print(
          f"[SUCCESS] Rendered {args.id} ->"
          f" {os.path.join(IMG_DIR, args.id + '.png')}"
      )
    else:
      print(f"[ERROR] {res}", file=sys.stderr)
      sys.exit(1)
  elif args.all:
    run_all_parallel(cores=args.cores, out_dir=args.out_dir)


if __name__ == "__main__":
  main()
