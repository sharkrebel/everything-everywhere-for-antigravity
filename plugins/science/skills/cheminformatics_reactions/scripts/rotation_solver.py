# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy",
#     "pandas",
#     "rdkit",
# ]
# ///
"""Automatic Reaction Rotation Solver for RDKit Reaction Mechanisms.

Autonomously determines optimal 2D ChemDraw/IUPAC orientations for reactants
(R1, R2) and products (P1, P2) in reaction mechanism diagrams:
  1. R1 reacting site naturally faces rightward (+x) towards R2.
  2. R2 reacting site naturally faces leftward (-x) towards R1.
  3. Intermolecular attack chord distance across '+' sign is minimized.
  4. Standard IUPAC/ChemDraw orientation priors (cardinal 0/90/180/270 deg,
     hexagonal 60/120/240/300 deg) prevent awkward arbitrary tilts.
  5. Product P1 (and P2) core scaffold is aligned to R1 (and R2) via
     Maximum Common Substructure (MCS) minimum-RMSD alignment, preserving
     chemical readability across the reaction arrow.
"""

import argparse
import json
import os
import sys
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
  sys.path.insert(0, SCRIPT_DIR)

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
import reaction_specs_100
import reactions_data

STANDARD_ANGLES = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]


def _resolve_mol_token(mol, token: int | str) -> int:
  """Resolves an atom index or `:N` atom-map token on `mol`."""
  if isinstance(token, str) and token.startswith(":"):
    if not token[1:].isdigit():
      raise ValueError(f"Invalid atom map label {token!r}")
    want = int(token[1:])
    if want <= 0:
      raise ValueError(f"Invalid atom map label {token!r}: must be >= 1")
    for a in mol.GetAtoms():
      if a.GetAtomMapNum() == want:
        return a.GetIdx()
    raise ValueError(f"Atom map label {token!r} matches no atom on molecule")
  return int(token)


def _has_map_token(mol, token: int | str) -> bool:
  """Returns True if `mol` owns the `:N` atom-map label `token`."""
  if not (
      isinstance(token, str) and token.startswith(":") and token[1:].isdigit()
  ):
    return False
  want = int(token[1:])
  return want > 0 and any(a.GetAtomMapNum() == want for a in mol.GetAtoms())


def _resolve_owned_tokens(mol_self, mol_other, tokens):
  """Resolves tokens owned by `mol_self`, skipping `:N` on partner reactants."""
  others = (
      (mol_other,)
      if isinstance(mol_other, Chem.Mol)
      else tuple(m for m in mol_other or () if m is not None)
  )
  out = []
  for tok in tokens:
    try:
      out.append(_resolve_mol_token(mol_self, tok))
    except ValueError:
      if (
          isinstance(tok, str)
          and tok.startswith(":")
          and any(_has_map_token(m, tok) for m in others)
      ):
        continue
      raise
  return out


def _parse_reactant_mol(smi: str):
  """Parses a SMILES string while preserving atom-mapped explicit hydrogens."""
  if "[H:" in smi or "[2H:" in smi or "[3H:" in smi:
    ps = Chem.SmilesParserParams()
    ps.removeHs = False
    mol = Chem.MolFromSmiles(smi, ps)
    if mol is not None:
      rm_params = Chem.RemoveHsParameters()
      rm_params.removeDegreeZero = False
      rm_params.removeOnlyHNeighbors = False
      rm_params.removeIsotopes = False
      rm_params.removeMapped = False
      rm_params.removeDefiningBondStereo = False
      rm_params.showWarnings = False
      mol = Chem.RemoveHs(mol, rm_params)
  else:
    mol = Chem.MolFromSmiles(smi)
  if mol is None:
    raise ValueError(f"Invalid SMILES: {smi}")
  return mol


def _extract_side_cfgs(rxn_def, prefix: str) -> dict[str, dict[str, Any]]:
  """Extracts keyed molecule config dicts (`r1..` or `p1..`) from `rxn_def`."""
  plural_key = "reactants" if prefix == "r" else "products"
  if plural_key in rxn_def:
    raw = rxn_def[plural_key]
    if isinstance(raw, list):
      return {
          item.get("key", f"{prefix}{idx + 1}"): dict(item)
          for idx, item in enumerate(raw)
          if isinstance(item, dict)
      }
    if isinstance(raw, dict):
      return {
          key: dict(item) for key, item in raw.items() if isinstance(item, dict)
      }
    raise ValueError(
        f"Invalid {plural_key} format: expected list or dict, got {type(raw)}"
    )
  keys = sorted(
      (
          key
          for key, item in rxn_def.items()
          if key.startswith(prefix)
          and key[len(prefix) :].isdigit()
          and isinstance(item, dict)
      ),
      key=lambda key: int(key[len(prefix) :]),
  )
  return {key: dict(rxn_def[key]) for key in keys}


def get_mol_atom_coords(
    smi,
    rot_deg=0.0,
    add_h=None,
    coord_map=None,
    w_cm=3.6,
    canvas_pos=(0.0, 0.0),
    w_px=350,
    h_px=300,
):
  """Computes 2D atom draw coordinates in cm in memory without disk I/O."""
  mol = _parse_reactant_mol(smi)

  if add_h is not None:
    rw = Chem.RWMol(mol)
    add_h_list = [add_h] if isinstance(add_h, (int, str)) else list(add_h)
    for h_tok in add_h_list:
      h_target = _resolve_mol_token(mol, h_tok)
      target_atom = rw.GetAtomWithIdx(h_target)
      num_h = target_atom.GetTotalNumHs()
      if num_h > 0:
        target_atom.SetNumExplicitHs(max(0, num_h - 1))
        target_atom.SetNoImplicit(True)
      h_idx = rw.AddAtom(Chem.Atom(1))
      rw.AddBond(h_target, h_idx, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)

  if coord_map is not None:
    pt_map = {}
    for k, v in coord_map.items():
      k_idx = _resolve_mol_token(mol, k)
      if isinstance(v, (tuple, list)):
        pt_map[k_idx] = Point2D(float(v[0]), float(v[1]))
      else:
        pt_map[k_idx] = v
    rdDepictor.Compute2DCoords(mol, coordMap=pt_map)
  else:
    rdDepictor.Compute2DCoords(mol)
    regularize_nonconvex_small_rings(mol)

  d2d = rdMolDraw2D.MolDraw2DCairo(w_px, h_px)
  opts = d2d.drawOptions()
  opts.clearBackground = True
  opts.bondLineWidth = 2.0
  if float(rot_deg) != 0.0:
    opts.rotate = float(rot_deg)
  d2d.DrawMolecule(mol)
  d2d.FinishDrawing()

  scale = float(w_cm) / float(w_px)
  canvas_x, canvas_y = float(canvas_pos[0]), float(canvas_pos[1])

  atoms = {}
  for i in range(mol.GetNumAtoms()):
    pt = d2d.GetDrawCoords(i)
    cx = canvas_x + (pt.x - w_px / 2.0) * scale
    cy = canvas_y - (pt.y - h_px / 2.0) * scale
    atoms[i] = np.array([cx, cy])

  return mol, atoms


def extract_reacting_atoms(spec):
  """Extracts r1/r2 reacting atom tokens and any intermolecular attack."""
  by_mol = {"r1": set(), "r2": set()}
  inter_att = None

  for arr in spec.get("arrows", []):
    atype = arr.get("type")
    if atype == "attack":
      src, tgt = arr.get("src"), arr.get("tgt")
      if src and tgt:
        s_mol, s_idx = src
        t_mol, t_idx = tgt
        s_toks = (s_idx,) if isinstance(s_idx, (int, str)) else tuple(s_idx)
        t_toks = (t_idx,) if isinstance(t_idx, (int, str)) else tuple(t_idx)
        if s_mol in by_mol:
          by_mol[s_mol].update(s_toks)
        if t_mol in by_mol:
          by_mol[t_mol].update(t_toks)
        if {s_mol, t_mol} == {"r1", "r2"} and inter_att is None:
          inter_att = arr
    elif atype == "bond_to_bond":
      s_mol = arr.get("src_mol") or arr.get("mol", "r1")
      t_mol = arr.get("tgt_mol") or arr.get("mol", s_mol)
      b1 = tuple(arr.get("b1") or arr.get("src_bond") or ())
      b2 = tuple(arr.get("b2") or arr.get("tgt_bond") or ())
      if s_mol in by_mol:
        by_mol[s_mol].update(b1)
      if t_mol in by_mol:
        by_mol[t_mol].update(b2)
    elif atype == "bond_cleave":
      mol_name = arr.get("mol", "r1")
      if mol_name in by_mol:
        by_mol[mol_name].update(arr.get("bond") or ())
        if arr.get("tgt") is not None:
          by_mol[mol_name].add(arr["tgt"])
    elif atype == "bond_to_atom":
      s_mol = arr.get("src_mol", "r1")
      t_mol = arr.get("tgt_mol", "r2")
      if s_mol in by_mol:
        by_mol[s_mol].update(arr.get("bond") or ())
      if t_mol in by_mol and arr.get("tgt_atom") is not None:
        by_mol[t_mol].add(arr["tgt_atom"])
    elif atype == "lone_pair_to_bond":
      mol_name = arr.get("mol", "r1")
      if mol_name in by_mol:
        if arr.get("atom") is not None:
          by_mol[mol_name].add(arr["atom"])
        by_mol[mol_name].update(arr.get("bond") or ())

  return by_mol["r1"], by_mol["r2"], inter_att


def _extract_intermolecular_endpoints(spec):
  """Extracts `(s_mol, s_toks, t_mol, t_toks)` tuples between `r1` and `r2`."""
  pairs = []
  for arr in spec.get("arrows", []):
    atype = arr.get("type")
    if atype == "attack":
      src, tgt = arr.get("src"), arr.get("tgt")
      if src and tgt and {src[0], tgt[0]} == {"r1", "r2"}:
        s_toks = (src[1],) if isinstance(src[1], (int, str)) else tuple(src[1])
        t_toks = (tgt[1],) if isinstance(tgt[1], (int, str)) else tuple(tgt[1])
        pairs.append((src[0], s_toks, tgt[0], t_toks))
    elif (
        atype == "bond_to_atom"
        and {arr.get("src_mol", "r1"), arr.get("tgt_mol", "r2")} == {"r1", "r2"}
        and arr.get("bond")
        and arr.get("tgt_atom") is not None
    ):
      pairs.append((
          arr.get("src_mol", "r1"),
          tuple(arr["bond"]),
          arr.get("tgt_mol", "r2"),
          (arr["tgt_atom"],),
      ))
    elif atype == "bond_to_bond":
      s_mol = arr.get("src_mol") or arr.get("mol", "r1")
      t_mol = arr.get("tgt_mol") or arr.get("mol", s_mol)
      b1 = tuple(arr.get("b1") or arr.get("src_bond") or ())
      b2 = tuple(arr.get("b2") or arr.get("tgt_bond") or ())
      if {s_mol, t_mol} == {"r1", "r2"} and b1 and b2:
        pairs.append((s_mol, b1, t_mol, b2))
  return pairs


def solve_reaction_rotations(rxn_def, spec, angles=None):
  """Solves optimal rotation angles (in degrees) for all reaction participants.

  Returns a dict, e.g. {'r1': 270, 'r2': 0, 'p1': 270, 'p2': 0}.
  """
  if angles is None:
    angles = STANDARD_ANGLES

  r_cfgs = _extract_side_cfgs(rxn_def, "r")
  p_cfgs = _extract_side_cfgs(rxn_def, "p")
  r1_cfg = r_cfgs["r1"]
  r2_cfg = r_cfgs.get("r2")
  p1_cfg = p_cfgs["p1"]
  p2_cfg = p_cfgs.get("p2")

  add_h_map = spec.get("add_h", {})
  coord_maps = spec.get("coord_map", {})
  r1_add_h = add_h_map.get("r1")
  r2_add_h = add_h_map.get("r2")

  r1_w = float(r1_cfg.get("w", 4.0))
  r2_w = float(r2_cfg.get("w", 3.0)) if r2_cfg else 0.0
  r2_x = (r1_w / 2.0) + 1.2 + (r2_w / 2.0)

  r1_rxn_atoms, r2_rxn_atoms, _ = extract_reacting_atoms(spec)
  extra_mols = [
      _parse_reactant_mol(cfg["smi"])
      for key, cfg in r_cfgs.items()
      if key not in ("r1", "r2") and cfg.get("smi")
  ]

  # 1. Precompute R1 coordinate sets
  r1_coords = {}
  mol_r1 = None
  for th in angles:
    mol, at_dict = get_mol_atom_coords(
        r1_cfg["smi"],
        rot_deg=th,
        add_h=r1_add_h,
        coord_map=coord_maps.get("r1"),
        w_cm=r1_w,
        canvas_pos=(0.0, 0.0),
    )
    r1_coords[th] = at_dict
    if mol_r1 is None:
      mol_r1 = mol

  if r2_cfg:
    r2_coords = {}
    mol_r2 = None
    for th in angles:
      mol, at_dict = get_mol_atom_coords(
          r2_cfg["smi"],
          rot_deg=th,
          add_h=r2_add_h,
          coord_map=coord_maps.get("r2"),
          w_cm=r2_w,
          canvas_pos=(r2_x, 0.0),
      )
      r2_coords[th] = at_dict
      if mol_r2 is None:
        mol_r2 = mol

    if mol_r1 is not None:
      r1_rxn_atoms = set(
          _resolve_owned_tokens(mol_r1, (mol_r2, *extra_mols), r1_rxn_atoms)
      )
    if mol_r2 is not None:
      r2_rxn_atoms = set(
          _resolve_owned_tokens(mol_r2, (mol_r1, *extra_mols), r2_rxn_atoms)
      )

    inter_pairs = []
    if mol_r1 is not None and mol_r2 is not None:
      for s_mol, s_toks, _, t_toks in _extract_intermolecular_endpoints(spec):
        r1_raw, r2_raw = (s_toks, t_toks) if s_mol == "r1" else (t_toks, s_toks)
        r1_ids = tuple(
            _resolve_owned_tokens(mol_r1, (mol_r2, *extra_mols), r1_raw)
        )
        r2_ids = tuple(
            _resolve_owned_tokens(mol_r2, (mol_r1, *extra_mols), r2_raw)
        )
        if r1_ids and r2_ids:
          inter_pairs.append((r1_ids, r2_ids))
      mols = {"r1": mol_r1, "r2": mol_r2}
      for arr in spec.get("arrows", []):
        own_k = arr.get("mol", "r1")
        if (
            arr.get("type") != "lone_pair_to_bond"
            or not arr.get("bond")
            or own_k not in mols
        ):
          continue
        oth_k = "r2" if own_k == "r1" else "r1"
        own_ids = _resolve_owned_tokens(
            mols[own_k], (mols[oth_k], *extra_mols), arr["bond"]
        )
        map_toks = [
            tok
            for tok in arr["bond"]
            if isinstance(tok, str) and tok.startswith(":")
        ]
        oth_ids = _resolve_owned_tokens(
            mols[oth_k], (mols[own_k], *extra_mols), map_toks
        )
        if own_ids and oth_ids:
          r1_ids, r2_ids = (
              (own_ids, oth_ids) if own_k == "r1" else (oth_ids, own_ids)
          )
          r1_rxn_atoms.update(r1_ids)
          r2_rxn_atoms.update(r2_ids)
          inter_pairs.append((tuple(r1_ids), tuple(r2_ids)))

    is_cycloaddition = False
    if len(inter_pairs) >= 2 and angles:

      def _pair_loci(at_ref, side_idx):
        return [
            np.mean([at_ref[i] for i in pair[side_idx]], axis=0)
            for pair in inter_pairs
            if all(i in at_ref for i in pair[side_idx])
        ]

      r1_loci = _pair_loci(r1_coords[angles[0]], 0)
      r2_loci = _pair_loci(r2_coords[angles[0]], 1)
      if len(r1_loci) >= 2 and len(r2_loci) >= 2:
        is_cycloaddition = (
            float(np.linalg.norm(r1_loci[0] - r1_loci[1])) > 0.15
            and float(np.linalg.norm(r2_loci[0] - r2_loci[1])) > 0.15
        )

    best_score = -1e9
    best_th1 = 0
    best_th2 = 0

    for th1 in angles:
      at1 = r1_coords[th1]
      all_pts1 = np.array(list(at1.values()))
      c1 = np.mean(all_pts1, axis=0)

      if r1_rxn_atoms:
        pts_rxn1 = np.array([at1[idx] for idx in r1_rxn_atoms if idx in at1])
        v1 = (
            np.mean(pts_rxn1, axis=0) - c1
            if len(pts_rxn1) > 0
            else np.array([1.0, 0.0])
        )
      else:
        v1 = np.array([1.0, 0.0])

      score_r1 = (
          v1[0] * 4.5 - abs(v1[1]) * 2.0
          if is_cycloaddition
          else v1[0] * 4.0 + max(0.0, v1[1]) * 1.5
      )
      prior_r1 = 0.0 if th1 % 90 == 0 else (-0.3 if th1 % 60 == 0 else -0.8)

      for th2 in angles:
        at2 = r2_coords[th2]
        all_pts2 = np.array(list(at2.values()))
        c2 = np.mean(all_pts2, axis=0)

        if r2_rxn_atoms:
          pts_rxn2 = np.array([at2[idx] for idx in r2_rxn_atoms if idx in at2])
          v2 = (
              np.mean(pts_rxn2, axis=0) - c2
              if len(pts_rxn2) > 0
              else np.array([-1.0, 0.0])
          )
        else:
          v2 = np.array([-1.0, 0.0])

        score_r2 = (
            -v2[0] * 4.5 - abs(v2[1]) * 2.0
            if is_cycloaddition
            else -v2[0] * 4.0 + max(0.0, v2[1]) * 1.5
        )
        prior_r2 = 0.0 if th2 % 90 == 0 else (-0.3 if th2 % 60 == 0 else -0.8)

        # Attack chord distance across the '+' sign
        pair_slice = inter_pairs if is_cycloaddition else inter_pairs[:1]
        chords = [
            float(
                np.linalg.norm(
                    np.mean([at2[j] for j in r2_ids], axis=0)
                    - np.mean([at1[i] for i in r1_ids], axis=0)
                )
            )
            for r1_ids, r2_ids in pair_slice
            if all(i in at1 for i in r1_ids) and all(j in at2 for j in r2_ids)
        ]
        dist_score = -2.5 * float(np.mean(chords)) if chords else 0.0

        total = score_r1 + score_r2 + dist_score + prior_r1 + prior_r2
        if total > best_score:
          best_score = total
          best_th1 = th1
          best_th2 = th2
  else:
    # Unimolecular reaction
    if mol_r1 is not None:
      r1_rxn_atoms = set(
          _resolve_owned_tokens(mol_r1, extra_mols, r1_rxn_atoms)
      )
    best_th2 = 0
    best_score = -1e9
    best_th1 = 0
    for th1 in angles:
      at1 = r1_coords[th1]
      all_pts1 = np.array(list(at1.values()))
      c1 = np.mean(all_pts1, axis=0)
      if r1_rxn_atoms:
        pts_rxn1 = np.array([at1[idx] for idx in r1_rxn_atoms if idx in at1])
        v1 = (
            np.mean(pts_rxn1, axis=0) - c1
            if len(pts_rxn1) > 0
            else np.array([0.0, 1.0])
        )
      else:
        v1 = np.array([0.0, 1.0])
      score = v1[1] * 3.0 + (
          0.0 if th1 % 90 == 0 else (-0.3 if th1 % 60 == 0 else -0.8)
      )
      if score > best_score:
        best_score = score
        best_th1 = th1

  # 2. Align Product P1 to R1 via Maximum Common Substructure (MCS)
  best_th_p1 = best_th1
  try:
    m1_rdk, _ = get_mol_atom_coords(
        r1_cfg["smi"],
        rot_deg=0,
        add_h=add_h_map.get("r1"),
        coord_map=coord_maps.get("r1"),
    )
    mp1_rdk, _ = get_mol_atom_coords(p1_cfg["smi"], rot_deg=0)
    mcs = rdFMCS.FindMCS([m1_rdk, mp1_rdk], timeout=2)
    mcs_any = rdFMCS.FindMCS(
        [m1_rdk, mp1_rdk],
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=True,
        timeout=2,
    )
    if mcs_any.numAtoms > mcs.numAtoms:
      mcs = mcs_any
    if mcs.numAtoms >= 3:
      q = Chem.MolFromSmarts(mcs.smartsString)
      match1 = m1_rdk.GetSubstructMatch(q)
      matchp = mp1_rdk.GetSubstructMatch(q)

      at1_opt = r1_coords[best_th1]
      pts_r1 = np.array([at1_opt[idx] for idx in match1 if idx in at1_opt])
      if len(pts_r1) == len(match1):
        c_r1 = np.mean(pts_r1, axis=0)
        pts_r1_c = pts_r1 - c_r1
        norm_r1 = max(1e-6, float(np.linalg.norm(pts_r1_c)))

        best_rmsd = 1e9
        p1_w = float(p1_cfg.get("w", 4.0))
        for th_p in angles:
          _, at_p = get_mol_atom_coords(p1_cfg["smi"], rot_deg=th_p, w_cm=p1_w)
          pts_p = np.array([at_p[idx] for idx in matchp if idx in at_p])
          if len(pts_p) == len(matchp):
            c_p = np.mean(pts_p, axis=0)
            pts_p_c = pts_p - c_p
            norm_p = max(1e-6, float(np.linalg.norm(pts_p_c)))
            pts_p_scaled = pts_p_c * (norm_r1 / norm_p)
            rmsd = np.sqrt(
                np.mean(np.sum((pts_p_scaled - pts_r1_c) ** 2, axis=1))
            )
            if rmsd < best_rmsd:
              best_rmsd = rmsd
              best_th_p1 = th_p
  except (ValueError, RuntimeError, KeyError, ZeroDivisionError) as exc:
    sys.stderr.write(
        "[rotation_solver] Warning: P1-to-R1 MCS rotation alignment failed"
        f" (keeping default rotation {best_th_p1}°): {exc}\n"
    )

  # 3. Align Product P2 to R2 via MCS if present
  best_th_p2 = 0
  if p2_cfg and r2_cfg:
    best_th_p2 = best_th2
    try:
      m2_rdk, _ = get_mol_atom_coords(
          r2_cfg["smi"],
          rot_deg=0,
          add_h=add_h_map.get("r2"),
          coord_map=coord_maps.get("r2"),
      )
      mp2_rdk, _ = get_mol_atom_coords(p2_cfg["smi"], rot_deg=0)
      mcs2 = rdFMCS.FindMCS([m2_rdk, mp2_rdk], timeout=2)
      mcs2_any = rdFMCS.FindMCS(
          [m2_rdk, mp2_rdk],
          bondCompare=rdFMCS.BondCompare.CompareAny,
          ringMatchesRingOnly=True,
          timeout=2,
      )
      if mcs2_any.numAtoms > mcs2.numAtoms:
        mcs2 = mcs2_any
      if mcs2.numAtoms >= 2:
        q2 = Chem.MolFromSmarts(mcs2.smartsString)
        match2 = m2_rdk.GetSubstructMatch(q2)
        matchp2 = mp2_rdk.GetSubstructMatch(q2)

        at2_opt = r2_coords[best_th2]
        pts_r2 = np.array([at2_opt[idx] for idx in match2 if idx in at2_opt])
        if len(pts_r2) == len(match2):
          c_r2 = np.mean(pts_r2, axis=0)
          pts_r2_c = pts_r2 - c_r2
          norm_r2 = max(1e-6, float(np.linalg.norm(pts_r2_c)))

          best_rmsd2 = 1e9
          p2_w = float(p2_cfg.get("w", 2.2))
          for th_p2 in angles:
            _, at_p2 = get_mol_atom_coords(
                p2_cfg["smi"], rot_deg=th_p2, w_cm=p2_w
            )
            pts_p2 = np.array([at_p2[idx] for idx in matchp2 if idx in at_p2])
            if len(pts_p2) == len(matchp2):
              c_p2 = np.mean(pts_p2, axis=0)
              pts_p2_c = pts_p2 - c_p2
              norm_p2 = max(1e-6, float(np.linalg.norm(pts_p2_c)))
              pts_p2_scaled = pts_p2_c * (norm_r2 / norm_p2)
              rmsd2 = np.sqrt(
                  np.mean(np.sum((pts_p2_scaled - pts_r2_c) ** 2, axis=1))
              )
              if rmsd2 < best_rmsd2:
                best_rmsd2 = rmsd2
                best_th_p2 = th_p2
    except (ValueError, RuntimeError, KeyError, ZeroDivisionError) as exc:
      sys.stderr.write(
          "[rotation_solver] Warning: P2-to-R2 MCS rotation alignment failed"
          f" (keeping default rotation {best_th_p2}°): {exc}\n"
      )

  solved_rot = {"r1": best_th1, "p1": best_th_p1}
  if r2_cfg:
    solved_rot["r2"] = best_th2
  if p2_cfg:
    solved_rot["p2"] = best_th_p2

  return solved_rot


def main() -> None:
  """CLI entry point.

  This module is primarily a library used by `render_mechanism.py`. It
  previously had no argparse at all, so `--help` was ignored and the script
  silently ran a hardcoded 10-reaction demo — indistinguishable, to a caller,
  from a tool that had processed their input.
  """
  parser = argparse.ArgumentParser(
      description=(
          "Solve 2D depiction rotations for reaction components. Primarily a "
          "library for render_mechanism.py; this CLI is for inspection."
      )
  )
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument(
      "--reaction-id",
      action="append",
      default=[],
      help="Catalog reaction id to solve, e.g. rxn_001 (repeatable).",
  )
  group.add_argument(
      "--all", action="store_true", help="Solve every reaction in the catalog."
  )
  parser.add_argument(
      "--limit", type=int, default=None, help="Stop after N reactions."
  )
  parser.add_argument(
      "--output-json", default=None, help="Write the solved rotations to JSON."
  )
  args = parser.parse_args()

  by_id = {r["id"]: r for r in reactions_data.REACTIONS}
  if args.all:
    selected = list(reactions_data.REACTIONS)
  else:
    selected = []
    for rid in args.reaction_id:
      if rid not in by_id:
        print(
            f"Error: unknown reaction id '{rid}'. The catalog holds "
            f"{len(by_id)} reactions, rxn_001..rxn_{len(by_id):03d}.",
            file=sys.stderr,
        )
        sys.exit(1)
      selected.append(by_id[rid])
  if args.limit is not None:
    selected = selected[: args.limit]

  solved_all = {}
  for rxn in selected:
    rid = rxn["id"]
    spec = reaction_specs_100.SPECS.get(rid, {})
    solved = solve_reaction_rotations(rxn, spec)
    solved_all[rid] = solved
    print(
        f"{rid} ({rxn.get('title', '')}): solved={solved} | "
        f"spec={spec.get('rot', {})}"
    )

  if args.output_json:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_json)) or ".", exist_ok=True
    )
    with open(args.output_json, "w") as fh:
      json.dump(solved_all, fh, indent=2)
      fh.write("\n")
    print(f"Wrote solved rotations to: {args.output_json}")


if __name__ == "__main__":
  main()
