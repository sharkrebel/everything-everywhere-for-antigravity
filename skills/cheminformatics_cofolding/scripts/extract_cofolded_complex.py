# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "gemmi",
#     "rdkit",
#     "numpy",
#     "pandas",
# ]
# ///

"""AlphaFold 3 mmCIF Pose Extractor and Protein-Ligand Interaction Analyzer.

Parses AlphaFold 3 output mmCIF structure files, isolates predicted
small-molecule ligand poses into RDKit SDF files using template-based bond
assignment, and computes pocket contacts, Protein-Ligand Interaction
Fingerprints (PLIF: hydrogen bonds, salt bridges, pi-stacking, hydrophobic
contacts), clashes, and pLDDT/ipTM confidence metrics.
"""

import argparse
import json
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compute_plif  # noqa: E402

_STANDARD_AMINO_ACIDS = frozenset({
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
})


def parse_mmcif_atom_sites(cif_path: str) -> List[Dict[str, Any]]:
  """Lightweight, robust parser for mmCIF _atom_site loop tables."""
  atoms = []
  in_atom_site = False
  columns = []

  with open(cif_path, "r", encoding="utf-8") as f:
    for line in f:
      line_str = line.strip()
      if not line_str or line_str.startswith("#"):
        continue

      if line_str.startswith("loop_"):
        in_atom_site = False
        columns = []
        continue

      if line_str.startswith("_atom_site."):
        in_atom_site = True
        col_name = line_str.split(".")[1]
        columns.append(col_name)
        continue

      if in_atom_site:
        if line_str.startswith("_") or line_str.startswith("loop_"):
          in_atom_site = False
          continue

        parts = line_str.split()
        if len(parts) >= len(columns):
          rec = {columns[i]: parts[i] for i in range(len(columns))}
          atoms.append(rec)
  return atoms


def _perceive_connectivity(rw_mol):
  """Infers bonds from 3D coordinates, in place, preserving atom order.

  Args:
      rw_mol: RWMol carrying a single conformer and no bonds.

  Returns:
      A Mol with geometry-perceived single bonds.
  """
  m = Chem.Mol(rw_mol)
  rdDetermineBonds.DetermineConnectivity(m)
  return m


def _finalize_valences(m):
  """Clears the placeholder radicals left by connectivity perception.

  DetermineConnectivity marks every unfilled valence as a radical and
  sets noImplicit, which suppresses implicit hydrogens and makes every
  such atom serialize as a bracket atom. Resetting these lets RDKit
  re-derive the correct implicit H counts.

  Args:
      m: Mol to fix up in place.

  Returns:
      The sanitized molecule.
  """
  for i in range(m.GetNumAtoms()):
    a = m.GetAtomWithIdx(i)
    a.SetNoImplicit(False)
    a.SetNumExplicitHs(0)
    a.SetNumRadicalElectrons(0)
  m.UpdatePropertyCache(strict=False)
  Chem.SanitizeMol(m)
  return m


def _cif_value(rec: Dict[str, Any], *keys: str, default: str = "") -> str:
  """First meaningful value among `keys`, treating mmCIF '?' and '.' as unset."""
  for k in keys:
    v = rec.get(k)
    if v not in (None, "", "?", "."):
      return str(v)
  return default


def _record_pos(rec: Dict[str, Any]) -> np.ndarray:
  return np.array([
      float(rec.get("Cartn_x", 0.0)),
      float(rec.get("Cartn_y", 0.0)),
      float(rec.get("Cartn_z", 0.0)),
  ])


def compute_plif_interactions(
    ligand_atoms: List[Dict[str, Any]],
    protein_atoms: List[Dict[str, Any]],
    extracted_mol: Optional[Chem.Mol] = None,
    ligand_coords: Optional[np.ndarray] = None,
    hbond_cutoff: float = 3.5,
    salt_bridge_cutoff: float = 4.0,
    pi_stack_cutoff: float = 5.5,
    hydrophobic_cutoff: float = 4.0,
) -> List[Dict[str, Any]]:
  """Compute a PLIF from raw mmCIF atom records.

  This is a thin adapter: it reshapes the mmCIF records into the plain
  dictionaries that :func:`compute_plif.detect_interactions` expects and then
  defers entirely to that single implementation of the geometric rules.

  Ligand element identity, aromatic ring perception and formal charges are
  taken from `extracted_mol` (which carries template-assigned bond orders),
  not guessed from the mmCIF, so pi-stacking and salt bridges are evaluated
  against real chemistry. When `ligand_coords` is supplied it overrides the
  mmCIF positions, keeping the fingerprint consistent with the pre-relaxed
  pose that is actually written to the output SDF.

  Args:
    ligand_atoms: mmCIF `_atom_site` records for the ligand chain, in the same
      order as the atoms of `extracted_mol`.
    protein_atoms: mmCIF `_atom_site` records for the protein partner.
    extracted_mol: the reconstructed ligand with correct bond orders.
    ligand_coords: optional (N, 3) array of pre-relaxed ligand coordinates.

  Returns:
    A deterministically sorted list of interaction records.
  """
  # --- Protein side -------------------------------------------------------
  prot_records: List[Dict[str, Any]] = []
  residue_atom_index: Dict[str, List[Dict[str, Any]]] = {}
  residue_atom_coords: Dict[str, Dict[str, np.ndarray]] = {}
  residue_names: Dict[str, str] = {}

  for pa in protein_atoms:
    chain = _cif_value(pa, "auth_asym_id", "label_asym_id")
    res_name = _cif_value(pa, "label_comp_id", "auth_comp_id")
    res_num = _cif_value(pa, "auth_seq_id", "label_seq_id")
    atom_name = _cif_value(pa, "label_atom_id", "auth_atom_id").strip('"')
    element = _cif_value(pa, "type_symbol", default="C").upper()
    pos = _record_pos(pa)
    res_id = f"{chain}:{res_name}{res_num}"

    residue_names[res_id] = res_name
    # Hydrogens are excluded from contact counting but are required for
    # D-H...A directionality, so they stay in the per-residue index.
    residue_atom_index.setdefault(res_id, []).append(
        {"pos": pos, "element": element.capitalize(), "name": atom_name}
    )
    residue_atom_coords.setdefault(res_id, {})[atom_name] = pos

    if element == "H":
      continue
    prot_records.append({
        "chain": chain,
        "res_id": res_id,
        "res_num": res_num,
        "res_name": res_name,
        "atom_name": atom_name,
        "element": element.capitalize(),
        "pos": pos,
    })

  protein_rings = {
      res_id: compute_plif.ring_centroids_from_atoms(
          res_name, residue_atom_coords.get(res_id, {})
      )
      for res_id, res_name in residue_names.items()
      if res_name in compute_plif.AROMATIC_RESIDUES
  }
  protein_rings = {k: v for k, v in protein_rings.items() if v}

  # --- Ligand side --------------------------------------------------------
  coords = ligand_coords
  if (
      coords is None
      and extracted_mol is not None
      and extracted_mol.GetNumConformers()
  ):
    conf = extracted_mol.GetConformer()
    coords = np.array([
        list(conf.GetAtomPosition(i))
        for i in range(extracted_mol.GetNumAtoms())
    ])
  if coords is None:
    coords = np.array([_record_pos(la) for la in ligand_atoms])

  lig_records: List[Dict[str, Any]] = []
  n_lig = min(len(ligand_atoms), len(coords))
  for i in range(n_lig):
    la = ligand_atoms[i]
    name = _cif_value(
        la, "label_atom_id", "auth_atom_id", default=f"A{i}"
    ).strip('"')
    element = _cif_value(la, "type_symbol", default="C").capitalize()
    charge = 0
    is_don = element in ("N", "O")
    is_acc = element in compute_plif.HBOND_DONORS_ACCEPTORS
    if extracted_mol is not None and i < extracted_mol.GetNumAtoms():
      atom = extracted_mol.GetAtomWithIdx(i)
      element = atom.GetSymbol()
      charge = atom.GetFormalCharge()
      is_don, is_acc = compute_plif.classify_ligand_hbond_atom(atom)
    lig_records.append({
        "name": name,
        "element": element,
        "pos": np.asarray(coords[i], dtype=float),
        "formal_charge": charge,
        "is_donor": is_don,
        "is_acceptor": is_acc,
    })

  lig_rings: List[Tuple[np.ndarray, np.ndarray]] = []
  if extracted_mol is not None:
    for ring in extracted_mol.GetRingInfo().AtomRings():
      if not all(extracted_mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
        continue
      if any(i >= len(coords) for i in ring):
        continue
      pts = np.array([coords[i] for i in ring], dtype=float)
      center = pts.mean(axis=0)
      normal = np.cross(pts[1] - pts[0], pts[2] - pts[0])
      normal = normal / (np.linalg.norm(normal) + 1e-8)
      lig_rings.append((center, normal))

  interactions, _ = compute_plif.detect_interactions(
      protein_atoms=prot_records,
      ligand_atoms=lig_records,
      residue_atom_index=residue_atom_index,
      protein_rings=protein_rings,
      lig_rings=lig_rings,
      hbond_cutoff=hbond_cutoff,
      hydrophobic_cutoff=hydrophobic_cutoff,
      salt_bridge_cutoff=salt_bridge_cutoff,
      pi_stack_cutoff=pi_stack_cutoff,
  )

  interactions.sort(
      key=lambda x: (
          x["chain"],
          int(x["residue_number"]) if str(x["residue_number"]).isdigit() else 0,
          x["interaction_type"],
          x["distance_angstrom"],
      )
  )
  return interactions


def prerelax_ligand_coordinates(
    lig_coords: np.ndarray,
    prot_coords: np.ndarray,
    clash_cutoff: float = 2.0,
    max_iters: int = 50,
    step_size: float = 0.05,
    k_restraint: float = 1.0,
) -> np.ndarray:
  """Pre-relaxes ligand coordinates in static protein pocket using restrained gradient descent

  to relieve micro-clashes (<2.0 A) before contact and clash evaluation.
  """
  if len(lig_coords) == 0 or len(prot_coords) == 0:
    return lig_coords
  coords = np.array(lig_coords, dtype=float).copy()
  init_coords = coords.copy()

  for _ in range(max_iters):
    forces = np.zeros_like(coords)
    any_clash = False
    for i, lc in enumerate(coords):
      dists = np.linalg.norm(prot_coords - lc, axis=1)
      clash_mask = (dists < clash_cutoff) & (dists > 1e-4)
      if np.any(clash_mask):
        any_clash = True
        diffs = lc - prot_coords[clash_mask]
        d_vals = dists[clash_mask, np.newaxis]
        repulsion = (diffs / d_vals) * (clash_cutoff - d_vals) * 4.0
        forces[i] += np.sum(repulsion, axis=0)

      # Restraint to initial pose (harmonic k=1.0)
      forces[i] -= k_restraint * (lc - init_coords[i])

    if not any_clash:
      break

    coords += step_size * np.clip(forces, -1.0, 1.0)

  return coords


def extract_ligand_and_pocket(
    cif_path: str,
    ligand_chain: Optional[str] = None,
    protein_chain: Optional[str] = None,
    ligand_smiles: Optional[str] = None,
    confidences_json: Optional[str] = None,
    output_sdf: Optional[str] = None,
    output_csv: Optional[str] = None,
    pocket_radius_angstrom: float = 4.5,
) -> Dict[str, Any]:
  """Isolates ligand from AF3 mmCIF, calculates pocket contacts, PLIF, and confidences."""
  atoms = parse_mmcif_atom_sites(cif_path)
  if not atoms:
    raise ValueError(f"No atom records found in mmCIF file: {cif_path}")

  solvent_comps = {
      "HOH",
      "WAT",
      "CL",
      "NA",
      "SO4",
      "GOL",
      "EDO",
      "DMS",
      "ACT",
      "PO4",
  }

  # Determine ligand chain if not specified
  all_chains = sorted(
      list(
          set(a.get("auth_asym_id", a.get("label_asym_id", "")) for a in atoms)
      )
  )
  target_chain = ligand_chain
  if not target_chain:
    candidate_chains = []
    for a in atoms:
      comp = a.get("label_comp_id", a.get("auth_comp_id", ""))
      chain = a.get("label_asym_id", a.get("auth_asym_id", ""))
      if (
          comp not in _STANDARD_AMINO_ACIDS
          and comp not in solvent_comps
          and chain not in candidate_chains
      ):
        candidate_chains.append(chain)
    target_chain = candidate_chains[0] if candidate_chains else all_chains[-1]

  def _is_target_lig(a: Dict[str, Any]) -> bool:
    if a.get("auth_asym_id", "") == target_chain:
      if protein_chain and target_chain == protein_chain:
        comp = a.get("label_comp_id", a.get("auth_comp_id", ""))
        return comp not in _STANDARD_AMINO_ACIDS and comp not in solvent_comps
      return True
    return (
        a.get("label_asym_id", "") == target_chain
        or a.get("label_comp_id", a.get("auth_comp_id", "")) == target_chain
    )

  ligand_atoms = [a for a in atoms if _is_target_lig(a)]
  if protein_chain:
    lig_auth_chains = [a.get("auth_asym_id", "") for a in ligand_atoms]
    if len(set(lig_auth_chains)) > 1:
      chosen_chain = (
          protein_chain
          if protein_chain in lig_auth_chains
          else lig_auth_chains[0]
      )
      ligand_atoms = [
          a for a in ligand_atoms if a.get("auth_asym_id", "") == chosen_chain
      ]
    protein_atoms = [
        a
        for a in atoms
        if a.get("auth_asym_id", a.get("label_asym_id", "")) == protein_chain
        and not _is_target_lig(a)
        and a.get("label_comp_id", a.get("auth_comp_id", ""))
        not in solvent_comps
    ]
  else:
    protein_atoms = [
        a
        for a in atoms
        if not _is_target_lig(a)
        and a.get("label_comp_id", a.get("auth_comp_id", ""))
        not in solvent_comps
    ]

  if not ligand_atoms:
    raise ValueError(
        f"No ligand atoms found for chain '{target_chain}'. Available chains:"
        f" {all_chains}"
    )

  if not ligand_chain:
    heavy_count = sum(
        1 for a in ligand_atoms if a.get("type_symbol", "C").capitalize() != "H"
    )
    if heavy_count > 150:
      sys.stderr.write(
          f"WARNING: Auto-detected ligand chain '{target_chain}' contains"
          f" {heavy_count} heavy atoms (>150). This strongly suggests a full"
          " macromolecule or protein chain was auto-selected as the ligand."
          " Please specify --ligand-chain explicitly to select the intended"
          " small molecule ligand.\n"
      )

  # Ligand coordinates & elements
  lig_coords = []
  lig_elements = []
  lig_plddts = []
  for a in ligand_atoms:
    x = float(a.get("Cartn_x", 0.0))
    y = float(a.get("Cartn_y", 0.0))
    z = float(a.get("Cartn_z", 0.0))
    elem = a.get("type_symbol", "C").capitalize()
    plddt = float(a.get("B_iso_or_equiv", 0.0))
    lig_coords.append([x, y, z])
    lig_elements.append(elem)
    lig_plddts.append(plddt)

  lig_coords_arr = np.array(lig_coords)
  mean_ligand_plddt = float(np.mean(lig_plddts)) if lig_plddts else 0.0

  # Build RDKit Mol for ligand.
  #
  # The coordinates arrive as a bare list of elements and positions with no
  # connectivity. Bond perception therefore has to happen in two stages:
  # first infer which atoms are bonded from the 3D geometry, then upgrade
  # those single bonds to the correct orders/aromaticity using the input
  # SMILES as a template. Skipping the first stage hands
  # AssignBondOrdersFromTemplate a bond-less molecule, which can never match
  # and yields an SDF with zero bonds.

  extracted_mol = None
  raw_mol = Chem.RWMol()
  for elem in lig_elements:
    raw_mol.AddAtom(Chem.Atom(elem))
  conf = Chem.Conformer(len(lig_elements))
  for i, (x, y, z) in enumerate(lig_coords):
    conf.SetAtomPosition(i, (x, y, z))
  raw_mol.AddConformer(conf)

  if ligand_smiles:
    ref_mol = Chem.MolFromSmiles(ligand_smiles)
    if ref_mol is None:
      raise ValueError(f"Invalid --ligand-smiles template: {ligand_smiles!r}")
    ref_no_h = Chem.RemoveHs(ref_mol)
    if raw_mol.GetNumAtoms() != ref_no_h.GetNumAtoms():
      raise ValueError(
          f"Ligand atom count mismatch: structure has {raw_mol.GetNumAtoms()}"
          " heavy atoms but --ligand-smiles template has"
          f" {ref_no_h.GetNumAtoms()}. Bond orders cannot be assigned; check"
          " that --ligand-chain selects only the ligand."
      )
    try:
      connected = _perceive_connectivity(raw_mol)
      extracted_mol = _finalize_valences(
          AllChem.AssignBondOrdersFromTemplate(ref_no_h, connected)
      )
    except Exception as exc:
      # Do not silently emit a bond-less pose: downstream PLIF and
      # pharmacophore scoring would see no rings, no charges and no
      # aromaticity, and would report confidently wrong numbers.
      print(
          "WARNING: could not assign bond orders from --ligand-smiles "
          f"({type(exc).__name__}: {exc}). Falling back to "
          "geometry-perceived single bonds; aromaticity and formal "
          "charges will be missing, so PLIF pi-stacking and salt-bridge "
          "counts from this pose are NOT reliable.",
          file=sys.stderr,
      )
      try:
        extracted_mol = _finalize_valences(_perceive_connectivity(raw_mol))
      except Exception:
        extracted_mol = raw_mol.GetMol()
  else:
    print(
        "WARNING: no --ligand-smiles given; bond orders will be inferred "
        "from geometry alone. Aromaticity and formal charges may be wrong. "
        "Pass --ligand-smiles for a chemically exact pose.",
        file=sys.stderr,
    )
    try:
      extracted_mol = _finalize_valences(_perceive_connectivity(raw_mol))
    except Exception:
      extracted_mol = raw_mol.GetMol()

  # Pre-relax ligand coordinates in static pocket to relieve micro-clashes (< 2.0 A)
  prot_heavy_coords = np.array([
      [
          float(pa.get("Cartn_x", 0.0)),
          float(pa.get("Cartn_y", 0.0)),
          float(pa.get("Cartn_z", 0.0)),
      ]
      for pa in protein_atoms
      if pa.get("type_symbol", "").upper() != "H"
  ])
  if len(prot_heavy_coords) > 0 and len(lig_coords_arr) > 0:
    lig_coords_arr = prerelax_ligand_coordinates(
        lig_coords_arr, prot_heavy_coords
    )
    if extracted_mol and extracted_mol.GetNumConformers() > 0:
      for i, (x, y, z) in enumerate(lig_coords_arr):
        extracted_mol.GetConformer().SetAtomPosition(
            i, (float(x), float(y), float(z))
        )

  # Identify binding pocket residues (< pocket_radius_angstrom)
  pocket_residues = set()
  hbond_count = 0
  clash_count = 0
  donor_acceptor_elements = {"N", "O", "S"}

  for pa in protein_atoms:
    px = float(pa.get("Cartn_x", 0.0))
    py = float(pa.get("Cartn_y", 0.0))
    pz = float(pa.get("Cartn_z", 0.0))
    pelem = pa.get("type_symbol", "").capitalize()
    p_coord = np.array([px, py, pz])

    dists = np.linalg.norm(lig_coords_arr - p_coord, axis=1)
    min_dist = float(np.min(dists))

    if min_dist <= pocket_radius_angstrom:
      res_name = pa.get("label_comp_id", pa.get("auth_comp_id", ""))
      res_num = pa.get("auth_seq_id", pa.get("label_seq_id", ""))
      chain = pa.get("auth_asym_id", pa.get("label_asym_id", ""))
      pocket_residues.add(f"{chain}:{res_name}{res_num}")

      if min_dist <= 3.5 and pelem in donor_acceptor_elements:
        closest_lig_idx = int(np.argmin(dists))
        if lig_elements[closest_lig_idx] in donor_acceptor_elements:
          hbond_count += 1

      if min_dist < 2.0 and pelem != "H":
        clash_count += 1

  # Compute geometric PLIF interactions. The rules themselves live in
  # compute_plif.py; passing the pre-relaxed coordinates keeps the
  # fingerprint consistent with the pose written to the output SDF.
  interactions = compute_plif_interactions(
      ligand_atoms=ligand_atoms,
      protein_atoms=protein_atoms,
      extracted_mol=extracted_mol,
      ligand_coords=lig_coords_arr,
      hbond_cutoff=3.5,
      salt_bridge_cutoff=4.0,
      pi_stack_cutoff=5.5,
      hydrophobic_cutoff=4.0,
  )

  plif_summary = {
      "hydrogen_bonds": sum(
          1 for i in interactions if i["interaction_type"] == "hydrogen_bond"
      ),
      "salt_bridges": sum(
          1 for i in interactions if i["interaction_type"] == "salt_bridge"
      ),
      "pi_stacking": sum(
          1 for i in interactions if i["interaction_type"] == "pi_stacking"
      ),
      "hydrophobic": sum(
          1 for i in interactions if i["interaction_type"] == "hydrophobic"
      ),
  }

  # Write PLIF CSV if requested
  if output_csv:
    compute_plif.write_plif_csv(interactions, output_csv)

  # Parse confidence summary if provided or auto-detected alongside CIF
  resolved_conf_path = confidences_json
  if not resolved_conf_path:
    cif_dir = os.path.dirname(os.path.abspath(cif_path))
    cif_base = os.path.splitext(os.path.basename(cif_path))[0]
    job_prefix = cif_base[:-6] if cif_base.endswith("_model") else cif_base
    auto_candidates = [
        os.path.join(cif_dir, f"{job_prefix}_confidences.json"),
        os.path.join(cif_dir, "metrics.json"),
        os.path.join(cif_dir, "run.json"),
        os.path.join(os.path.dirname(cif_dir), "run.json"),
        os.path.join(os.path.dirname(os.path.dirname(cif_dir)), "run.json"),
    ]
    for cand in auto_candidates:
      if os.path.isfile(cand):
        resolved_conf_path = cand
        break

  confidence_data: Dict[str, Any] = {}
  if resolved_conf_path and os.path.exists(resolved_conf_path):
    try:
      with open(resolved_conf_path, "r", encoding="utf-8") as f:
        c_data = json.load(f)
      if isinstance(c_data.get("output"), dict):
        out_block = c_data["output"]
        metrics_src = (out_block.get("best_sample") or {}).get("metrics") or {}
        binding_src = out_block.get("binding_metrics") or {}
      else:
        metrics_src = c_data
        binding_src = c_data.get("binding_metrics") or {}

      confidence_data = {
          "iptm": metrics_src.get("iptm"),
          "ptm": metrics_src.get("ptm"),
          "ranking_score": c_data.get(
              "ranking_score", metrics_src.get("structure_confidence")
          ),
          "has_clash": c_data.get("has_clash"),
      }
      for extra_key in (
          "ligand_iptm",
          "protein_iptm",
          "structure_confidence",
          "complex_plddt",
          "complex_iplddt",
          "complex_pde",
          "complex_ipde",
      ):
        if metrics_src.get(extra_key) is not None:
          confidence_data[extra_key] = metrics_src.get(extra_key)

      b_conf = c_data.get(
          "binding_confidence", binding_src.get("binding_confidence")
      )
      opt_score = c_data.get(
          "optimization_score", binding_src.get("optimization_score")
      )
      if b_conf is not None:
        confidence_data["binding_confidence"] = b_conf
      if opt_score is not None:
        confidence_data["optimization_score"] = opt_score
    except Exception:
      pass

  # Save to SDF if requested
  if output_sdf and extracted_mol:
    os.makedirs(os.path.dirname(os.path.abspath(output_sdf)), exist_ok=True)
    with Chem.SDWriter(output_sdf) as writer:
      extracted_mol.SetProp("AF3_Ligand_Chain", target_chain)
      extracted_mol.SetProp("AF3_Mean_pLDDT", f"{mean_ligand_plddt:.2f}")
      if confidence_data.get("iptm") is not None:
        extracted_mol.SetProp("AF3_ipTM", f"{confidence_data['iptm']:.3f}")
      if confidence_data.get("ptm") is not None:
        extracted_mol.SetProp("AF3_pTM", f"{confidence_data['ptm']:.3f}")
      if confidence_data.get("structure_confidence") is not None:
        extracted_mol.SetProp(
            "Boltz_Structure_Confidence",
            f"{confidence_data['structure_confidence']:.3f}",
        )
      if confidence_data.get("binding_confidence") is not None:
        extracted_mol.SetProp(
            "Boltz_Binding_Confidence",
            f"{confidence_data['binding_confidence']:.3f}",
        )
      if confidence_data.get("optimization_score") is not None:
        extracted_mol.SetProp(
            "Boltz_Optimization_Score",
            f"{confidence_data['optimization_score']:.3f}",
        )
      writer.write(extracted_mol)

  return {
      "cif_file": cif_path,
      "ligand_chain": target_chain,
      "ligand_atom_count": len(ligand_atoms),
      "mean_ligand_plddt": round(mean_ligand_plddt, 2),
      "pocket_residue_count": len(pocket_residues),
      "pocket_residues": sorted(list(pocket_residues)),
      "potential_hbonds": hbond_count,
      "severe_clashes": clash_count,
      "confidence_summary": confidence_data,
      "output_sdf": output_sdf,
      "output_csv": output_csv,
      "plif_interaction_count": len(interactions),
      "plif_summary": plif_summary,
  }


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Extract predicted ligand pose, pocket metrics, and PLIF from AF3"
          " mmCIF."
      )
  )
  parser.add_argument(
      "--cif",
      required=True,
      help="Input AlphaFold 3 output mmCIF file (*_model.cif).",
  )
  parser.add_argument(
      "--confidences", default=None, help="Path to summary_confidences_*.json."
  )
  parser.add_argument(
      "--ligand-smiles",
      default=None,
      help="Reference ligand SMILES for bond order reconstruction.",
  )
  parser.add_argument(
      "--ligand-chain",
      default=None,
      help="Chain ID of ligand (auto-detected if omitted).",
  )
  parser.add_argument(
      "--protein-chain",
      default=None,
      help="Chain ID of protein partner (defaults to all non-ligand chains).",
  )
  parser.add_argument(
      "--output-sdf",
      default=None,
      help="Output SDF file path for isolated ligand pose.",
  )
  parser.add_argument(
      "--output-csv",
      default=None,
      help=(
          "Output CSV file path for Protein-Ligand Interaction Fingerprint"
          " (PLIF)."
      ),
  )
  parser.add_argument(
      "--output",
      default=None,
      help=(
          "Unified output file path (SDF if ending in .sdf, CSV if ending in"
          " .csv)."
      ),
  )
  parser.add_argument(
      "--pocket-radius",
      type=float,
      default=4.5,
      help="Radius (Angstrom) for pocket residue contact definition.",
  )

  args = parser.parse_args()

  # Route unified --output flag if provided
  output_sdf = args.output_sdf
  output_csv = args.output_csv
  if args.output:
    if args.output.endswith(".csv") and not output_csv:
      output_csv = args.output
    elif args.output.endswith(".sdf") and not output_sdf:
      output_sdf = args.output

  try:
    report = extract_ligand_and_pocket(
        cif_path=args.cif,
        ligand_chain=args.ligand_chain,
        protein_chain=args.protein_chain,
        ligand_smiles=args.ligand_smiles,
        confidences_json=args.confidences,
        output_sdf=output_sdf,
        output_csv=output_csv,
        pocket_radius_angstrom=args.pocket_radius,
    )
    print(json.dumps(report, indent=2))
  except Exception as e:
    print(f"Error extracting AF3 complex: {e}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
  main()
