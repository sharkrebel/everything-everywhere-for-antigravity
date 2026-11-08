# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "gemmi",
# ]
# ///

"""Compute Protein-Ligand Interaction Fingerprints (PLIF) from mmCIF/PDB complexes."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

import gemmi
import numpy as np
from rdkit import Chem

INTERACTION_TYPES = [
    "hbond_p_don",  # Protein donor -> Ligand acceptor
    "hbond_l_don",  # Ligand donor -> Protein acceptor
    "salt_p_cat",  # Protein cation -> Ligand anion
    "salt_l_cat",  # Ligand cation -> Protein anion
    "pi_stack",  # Aromatic pi-pi stacking
    "hydrophobic",  # Non-polar hydrophobic contact
]

HBOND_DONORS_ACCEPTORS = {"N", "O", "F", "ZN", "MG", "MN", "FE", "CA"}
BASIC_RESIDUES = {"LYS", "ARG", "HIS"}
ACIDIC_RESIDUES = {"ASP", "GLU"}
AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS"}

# Only the atoms that actually carry the formal charge at physiological pH.
# Using whole-residue membership (any N of an Arg, any O of an Asp) counts
# backbone and side-chain amide atoms as ionic partners, which massively
# over-reports salt bridges.
CATIONIC_PROTEIN_ATOMS = {
    "ARG": {"NE", "NH1", "NH2"},
    "LYS": {"NZ"},
    "HIS": {
        "ND1",
        "NE2",
    },  # only when protonated; reported as a weak salt bridge
}
ANIONIC_PROTEIN_ATOMS = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
}

# Geometric thresholds (Angstroms / degrees).
HBOND_MIN_DH_A_ANGLE = 120.0  # D-H...A when explicit hydrogens are present
HBOND_MIN_ANTECEDENT_ANGLE = 90.0  # R-D...A fallback for hydrogen-free models
PI_PARALLEL_MAX_ANGLE = 30.0
PI_TSHAPED_MIN_ANGLE = 60.0
PI_PARALLEL_MAX_OFFSET = 2.0


def _angle_deg(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
  """Returns the angle a-b-c in degrees."""
  v1 = a - b
  v2 = c - b
  n1 = np.linalg.norm(v1)
  n2 = np.linalg.norm(v2)
  if n1 < 1e-8 or n2 < 1e-8:
    return 0.0
  cos = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
  return float(np.degrees(np.arccos(cos)))


def _interplanar_angle_deg(n1: np.ndarray, n2: np.ndarray) -> float:
  """Returns the acute angle between two ring planes, from their normals."""
  cos = abs(float(np.clip(np.dot(n1, n2), -1.0, 1.0)))
  return float(np.degrees(np.arccos(cos)))


def _ring_offset(c1: np.ndarray, n1: np.ndarray, c2: np.ndarray) -> float:
  """Lateral centroid offset of ring 2 projected onto the plane of ring 1."""
  d = c2 - c1
  perp = float(np.dot(d, n1))
  return float(np.sqrt(max(np.dot(d, d) - perp * perp, 0.0)))


def _hbond_geometry(
    donor_pos: np.ndarray,
    acceptor_pos: np.ndarray,
    residue_atoms: List[Dict[str, Any]],
) -> Tuple[float, str]:
  """Returns the best hydrogen bond angle and which criterion produced it.

  Predicted complexes (AlphaFold 3, Boltz-2) are usually delivered without
  explicit hydrogens, so two criteria are supported:

  - ``D-H...A``: used when a hydrogen is bonded to the donor. This is the
    real directionality test.
  - ``R-D...A``: fallback using the donor's antecedent heavy atom. It cannot
    confirm a hydrogen bond, but it does reject contacts where the acceptor
    sits behind the donor relative to its own substituent.

  Args:
    donor_pos: Coordinates of the protein donor/acceptor heavy atom.
    acceptor_pos: Coordinates of the ligand partner atom.
    residue_atoms: All atoms of the donor's residue, each with "pos" and
      "element" keys.

  Returns:
    Tuple of (angle in degrees, criterion label).
  """
  hydrogens = [
      a
      for a in residue_atoms
      if a["element"] == "H" and np.linalg.norm(a["pos"] - donor_pos) <= 1.3
  ]
  if hydrogens:
    return (
        max(_angle_deg(donor_pos, h["pos"], acceptor_pos) for h in hydrogens),
        "D-H...A",
    )

  antecedents = [
      a
      for a in residue_atoms
      if a["element"] not in ("H", "D")
      and 1e-3 < np.linalg.norm(a["pos"] - donor_pos) <= 1.8
  ]
  if antecedents:
    return (
        max(_angle_deg(a["pos"], donor_pos, acceptor_pos) for a in antecedents),
        "R-D...A",
    )

  # Isolated atom (ion, water oxygen): no directionality available.
  return 180.0, "distance_only"


def parse_structure(filepath: str) -> Any:
  """Loads a structure from mmCIF or PDB."""
  if filepath.endswith(".cif") or filepath.endswith(".mmcif"):
    doc = gemmi.cif.read(filepath)
    return gemmi.make_structure_from_block(doc.sole_block())
  return gemmi.read_structure(filepath)


def get_ring_centroids(res) -> List[Tuple[np.ndarray, np.ndarray]]:
  """Extract aromatic ring centroids and normal vectors for standard amino acids."""
  atoms = {a.name.strip(): np.array([a.pos.x, a.pos.y, a.pos.z]) for a in res}
  return ring_centroids_from_atoms(res.name, atoms)


def ring_centroids_from_atoms(
    rname: str, atoms: Dict[str, np.ndarray]
) -> List[Tuple[np.ndarray, np.ndarray]]:
  """Aromatic ring centroids/normals from a plain {atom_name: xyz} mapping.

  Kept separate from :func:`get_ring_centroids` so that callers which hold
  raw parsed coordinates rather than gemmi residues (e.g.
  ``extract_cofolded_complex.py``) can reuse the identical geometry instead
  of maintaining a second copy of it.
  """
  rings = []

  if rname in ("PHE", "TYR"):
    req = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
    if all(k in atoms for k in req):
      pts = np.array([atoms[k] for k in req])
      center = np.mean(pts, axis=0)
      v1 = pts[1] - pts[0]
      v2 = pts[2] - pts[0]
      normal = np.cross(v1, v2)
      normal /= np.linalg.norm(normal) + 1e-8
      rings.append((center, normal))
  elif rname == "TRP":
    req6 = ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"]
    if all(k in atoms for k in req6):
      pts = np.array([atoms[k] for k in req6])
      center = np.mean(pts, axis=0)
      v1 = pts[1] - pts[0]
      v2 = pts[2] - pts[0]
      normal = np.cross(v1, v2)
      normal /= np.linalg.norm(normal) + 1e-8
      rings.append((center, normal))
    req5 = ["CG", "CD1", "NE1", "CE2", "CD2"]
    if all(k in atoms for k in req5):
      pts = np.array([atoms[k] for k in req5])
      center = np.mean(pts, axis=0)
      v1 = pts[1] - pts[0]
      v2 = pts[2] - pts[0]
      normal = np.cross(v1, v2)
      normal /= np.linalg.norm(normal) + 1e-8
      rings.append((center, normal))
  elif rname == "HIS":
    req = ["CG", "ND1", "CD2", "CE1", "NE2"]
    if all(k in atoms for k in req):
      pts = np.array([atoms[k] for k in req])
      center = np.mean(pts, axis=0)
      v1 = pts[1] - pts[0]
      v2 = pts[2] - pts[0]
      normal = np.cross(v1, v2)
      normal /= np.linalg.norm(normal) + 1e-8
      rings.append((center, normal))

  return rings


def _protein_hbond_roles(
    res_name: str, atom_name: str, elem: str
) -> Tuple[bool, bool]:
  """Returns (is_donor, is_acceptor) for a protein heavy atom at physiological pH."""
  if elem not in ("N", "O") and not (res_name == "CYS" and atom_name == "SG"):
    return False, False
  if atom_name in ("O", "OXT"):
    return False, True
  if atom_name == "N":
    return (res_name != "PRO"), False
  if res_name == "ASP" and atom_name in ("OD1", "OD2"):
    return False, True
  if res_name == "GLU" and atom_name in ("OE1", "OE2"):
    return False, True
  if res_name in ("ASN", "GLN"):
    if atom_name in ("OD1", "OE1"):
      return False, True
    if atom_name in ("ND2", "NE2"):
      return True, False
  if res_name in ("ARG", "LYS", "TRP") and atom_name in (
      "NE",
      "NH1",
      "NH2",
      "NZ",
      "NE1",
  ):
    return True, False
  if res_name in ("SER", "THR", "TYR", "HIS", "CYS") and atom_name in (
      "OG",
      "OG1",
      "OH",
      "ND1",
      "NE2",
      "SG",
  ):
    return True, True
  return (elem == "N"), True


def classify_ligand_hbond_atom(atom: Chem.Atom) -> Tuple[bool, bool]:
  """Returns (is_donor, is_acceptor) for an RDKit ligand atom at physiological pH."""
  sym = atom.GetSymbol()
  if sym not in ("N", "O", "F"):
    return False, False
  is_amide_n = any(
      nbr.GetSymbol() in ("C", "S", "P")
      and any(
          b.GetBondType() == Chem.BondType.DOUBLE
          and b.GetOtherAtom(nbr).GetSymbol() in ("O", "S")
          for b in nbr.GetBonds()
      )
      for nbr in atom.GetNeighbors()
  )
  basic_sp3_n = (
      sym == "N"
      and not atom.GetIsAromatic()
      and atom.GetHybridization() == Chem.HybridizationType.SP3
      and not is_amide_n
  )
  is_donor = (
      (atom.GetTotalNumHs() > 0) or (atom.GetFormalCharge() > 0) or basic_sp3_n
  )
  is_acceptor = atom.GetFormalCharge() <= 0
  return is_donor, is_acceptor


def detect_interactions(
    protein_atoms: List[Dict[str, Any]],
    ligand_atoms: List[Dict[str, Any]],
    residue_atom_index: Dict[str, List[Dict[str, Any]]],
    protein_rings: Dict[str, List[Tuple[np.ndarray, np.ndarray]]],
    lig_rings: List[Tuple[np.ndarray, np.ndarray]],
    hbond_cutoff: float = 3.5,
    hydrophobic_cutoff: float = 4.0,
    salt_bridge_cutoff: float = 4.0,
    pi_stack_cutoff: float = 5.5,
) -> Tuple[List[Dict[str, Any]], Dict[str, Dict[str, int]]]:
  """Detect protein-ligand contacts from pre-parsed atom records.

  This is the single implementation of the interaction rules. It works on
  plain dictionaries rather than a particular structure-file reader so that
  every caller -- :func:`compute_plif` and ``extract_cofolded_complex.py``
  alike -- applies identical geometric criteria. Hydrogen bonds and
  pi-stacking are angle-gated, and salt bridges require a real formal charge
  on both partners rather than merely an acidic/basic residue name.

  Args:
    protein_atoms: dicts with keys ``chain``, ``res_id``, ``res_num``,
      ``res_name``, ``atom_name``, ``element``, ``pos`` (heavy atoms only).
    ligand_atoms: dicts with keys ``name``, ``element``, ``pos``,
      ``formal_charge``.
    residue_atom_index: ``res_id`` -> atom dicts *including* hydrogens, used for
      H-bond directionality.
    protein_rings: ``res_id`` -> list of ``(centroid, normal)`` for aromatic
      side chains.
    lig_rings: list of ``(centroid, normal)`` for the ligand's aromatic rings.

  Returns:
    ``(detailed_interactions, residue_fingerprints)``.
  """
  lig_has_carbon = any(la["element"] == "C" for la in ligand_atoms)

  detailed_interactions: List[Dict[str, Any]] = []
  residue_fingerprints: Dict[str, Dict[str, int]] = {}

  atoms_by_res: Dict[str, List[Dict[str, Any]]] = {}
  for a in protein_atoms:
    atoms_by_res.setdefault(a["res_id"], []).append(a)

  for p_res_id in sorted(atoms_by_res):
    res_contacts = {t: 0 for t in INTERACTION_TYPES}
    p_atoms = atoms_by_res[p_res_id]
    p_res_name = p_atoms[0]["res_name"]

    for pa in p_atoms:
      p_pos = pa["pos"]
      p_elem = pa["element"]
      p_rname = pa["res_name"]

      for la in ligand_atoms:
        l_pos = la["pos"]
        l_elem = la["element"]
        d = float(np.linalg.norm(p_pos - l_pos))

        # Hydrogen bonds: distance, donor-acceptor complementarity, and
        # directionality.
        if d <= hbond_cutoff and l_elem in HBOND_DONORS_ACCEPTORS:
          p_don, p_acc = _protein_hbond_roles(p_rname, pa["atom_name"], p_elem)
          l_don = la.get("is_donor", l_elem in ("N", "O"))
          l_acc = la.get("is_acceptor", l_elem in HBOND_DONORS_ACCEPTORS)
          can_p_donate = p_don and l_acc
          can_l_donate = l_don and p_acc
          if can_p_donate or can_l_donate:
            angle, criterion = _hbond_geometry(
                p_pos, l_pos, residue_atom_index.get(pa["res_id"], [])
            )
            ok = (
                angle >= HBOND_MIN_DH_A_ANGLE
                if criterion == "D-H...A"
                else angle >= HBOND_MIN_ANTECEDENT_ANGLE
            )
            if ok:
              if can_p_donate:
                res_contacts["hbond_p_don"] += 1
              else:
                res_contacts["hbond_l_don"] += 1
              detailed_interactions.append({
                  "chain": pa["chain"],
                  "residue_number": pa["res_num"],
                  "residue_name": p_rname,
                  "interaction_type": "hydrogen_bond",
                  "distance_angstrom": round(d, 2),
                  "angle_degrees": round(angle, 1),
                  "angle_criterion": criterion,
                  "ligand_atom": la["name"],
                  "protein_atom": pa["atom_name"],
              })

        # Salt bridges: both partners must actually carry a formal
        # charge. Matching on residue identity alone counts backbone
        # and amide atoms as ionic partners.
        if d <= salt_bridge_cutoff:
          p_atom_name = pa["atom_name"]
          p_cationic = p_atom_name in CATIONIC_PROTEIN_ATOMS.get(p_rname, ())
          p_anionic = p_atom_name in ANIONIC_PROTEIN_ATOMS.get(p_rname, ())
          if p_cationic and la.get("formal_charge", 0) < 0:
            res_contacts["salt_p_cat"] += 1
            detailed_interactions.append({
                "chain": pa["chain"],
                "residue_number": pa["res_num"],
                "residue_name": p_rname,
                "interaction_type": "salt_bridge",
                "distance_angstrom": round(d, 2),
                "ligand_atom": la["name"],
                "protein_atom": p_atom_name,
            })
          elif p_anionic and la.get("formal_charge", 0) > 0:
            res_contacts["salt_l_cat"] += 1
            detailed_interactions.append({
                "chain": pa["chain"],
                "residue_number": pa["res_num"],
                "residue_name": p_rname,
                "interaction_type": "salt_bridge",
                "distance_angstrom": round(d, 2),
                "ligand_atom": la["name"],
                "protein_atom": p_atom_name,
            })

        # Hydrophobic contacts
        if d <= hydrophobic_cutoff and p_elem == "C":
          if (lig_has_carbon and l_elem == "C") or (not lig_has_carbon):
            res_contacts["hydrophobic"] += 1
            detailed_interactions.append({
                "chain": pa["chain"],
                "residue_number": pa["res_num"],
                "residue_name": p_rname,
                "interaction_type": "hydrophobic",
                "distance_angstrom": round(d, 2),
                "ligand_atom": la["name"],
                "protein_atom": pa["atom_name"],
            })

    # Pi-Stacking: centroid proximity alone treats two rings meeting edge-on
    # at an unproductive angle as stacked. Require a genuine parallel
    # (sandwich / parallel-displaced) or T-shaped (edge-to-face) geometry.
    if lig_rings and protein_rings.get(p_res_id):
      p_rings = protein_rings[p_res_id]
      for p_cent, p_norm in p_rings:
        for l_cent, l_norm in lig_rings:
          c_dist = float(np.linalg.norm(p_cent - l_cent))
          if c_dist > pi_stack_cutoff:
            continue
          theta = _interplanar_angle_deg(p_norm, l_norm)
          offset = _ring_offset(p_cent, p_norm, l_cent)
          if (
              theta <= PI_PARALLEL_MAX_ANGLE
              and offset <= PI_PARALLEL_MAX_OFFSET
          ):
            geometry = "parallel"
          elif theta >= PI_TSHAPED_MIN_ANGLE:
            geometry = "t_shaped"
          else:
            continue
          res_contacts["pi_stack"] += 1
          detailed_interactions.append({
              "chain": p_res_id.split(":")[0],
              "residue_number": "".join(filter(str.isdigit, p_res_id)),
              "residue_name": p_res_name,
              "interaction_type": "pi_stacking",
              "distance_angstrom": round(c_dist, 2),
              "angle_degrees": round(theta, 1),
              "offset_angstrom": round(offset, 2),
              "geometry": geometry,
              "ligand_atom": "RING",
              "protein_atom": "RING",
          })

    if any(v > 0 for v in res_contacts.values()):
      residue_fingerprints[p_res_id] = res_contacts

  return detailed_interactions, residue_fingerprints


def compute_plif(
    complex_path: str,
    ligand_chain: str = "B",
    ligand_sdf_path: Optional[str] = None,
    hbond_cutoff: float = 3.5,
    hydrophobic_cutoff: float = 4.0,
    salt_bridge_cutoff: float = 4.0,
    pi_stack_cutoff: float = 5.5,
) -> Dict[str, Any]:
  """Computes comprehensive PLIF from structure."""
  st = parse_structure(complex_path)
  model = st[0]

  protein_atoms = []
  ligand_atoms = []
  protein_residues = {}
  # Per-residue atom lists (including hydrogens) used for H-bond geometry.
  residue_atom_index: Dict[str, List[Dict[str, Any]]] = {}

  solvent_names = {
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

  for chain in model:
    for res in chain:
      res_id = f"{chain.name}:{res.name}{res.seqid.num}"
      is_lig = (
          chain.name == ligand_chain
          or getattr(res, "subchain", "") == ligand_chain
          or res.name == ligand_chain
      )
      if is_lig:
        for atom in res:
          if atom.element.name == "H":
            continue
          pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
          ligand_atoms.append({
              "name": atom.name.strip(),
              "element": atom.element.name,
              "pos": pos,
              "res_name": res.name,
              "formal_charge": 0,
          })
      else:
        if res.name in solvent_names:
          continue
        if res_id not in protein_residues:
          protein_residues[res_id] = res
        for atom in res:
          pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
          residue_atom_index.setdefault(res_id, []).append({
              "pos": pos,
              "element": atom.element.name,
              "name": atom.name.strip(),
          })
          if atom.element.name == "H":
            continue
          protein_atoms.append({
              "chain": chain.name,
              "res_id": res_id,
              "res_num": str(res.seqid.num),
              "res_name": res.name,
              "atom_name": atom.name.strip(),
              "element": atom.element.name,
              "pos": pos,
          })

  if not ligand_atoms:
    raise ValueError(
        f"No ligand atoms found on chain '{ligand_chain}' in {complex_path}"
    )

  # Extract ligand aromatic rings and H-bond roles from SDF if provided
  lig_rings = []
  lig_atom_roles: List[Tuple[int, bool, bool, np.ndarray]] = []
  if ligand_sdf_path:
    suppl = Chem.SDMolSupplier(ligand_sdf_path)
    if suppl and len(suppl) > 0 and suppl[0] is not None:
      l_mol = suppl[0]
      conf = l_mol.GetConformer()
      for atom in l_mol.GetAtoms():
        fc = atom.GetFormalCharge()
        is_don, is_acc = classify_ligand_hbond_atom(atom)
        q = conf.GetAtomPosition(atom.GetIdx())
        lig_atom_roles.append((fc, is_don, is_acc, np.array([q.x, q.y, q.z])))
      ring_info = l_mol.GetRingInfo()
      for r_indices in ring_info.AtomRings():
        if all(l_mol.GetAtomWithIdx(i).GetIsAromatic() for i in r_indices):
          pts = np.array([conf.GetAtomPosition(i) for i in r_indices])
          center = np.mean(pts, axis=0)
          v1 = pts[1] - pts[0]
          v2 = pts[2] - pts[0]
          normal = np.cross(v1, v2)
          normal /= np.linalg.norm(normal) + 1e-8
          lig_rings.append((center, normal))

  if ligand_sdf_path and lig_atom_roles:
    for charge, is_don, is_acc, cpos in lig_atom_roles:
      nearest = min(
          ligand_atoms, key=lambda a: float(np.linalg.norm(a["pos"] - cpos))
      )
      if float(np.linalg.norm(nearest["pos"] - cpos)) <= 0.6:
        if charge:
          nearest["formal_charge"] = charge
        nearest["is_donor"] = is_don
        nearest["is_acceptor"] = is_acc

  if ligand_sdf_path and not any(a["formal_charge"] for a in ligand_atoms):
    print(
        "Note: ligand carries no formal charges; salt bridges will only be "
        "reported if the SDF encodes ionised groups.",
        file=sys.stderr,
    )

  protein_rings = {
      r_id: get_ring_centroids(r)
      for r_id, r in protein_residues.items()
      if r.name in AROMATIC_RESIDUES
  }

  detailed_interactions, residue_fingerprints = detect_interactions(
      protein_atoms=protein_atoms,
      ligand_atoms=ligand_atoms,
      residue_atom_index=residue_atom_index,
      protein_rings=protein_rings,
      lig_rings=lig_rings,
      hbond_cutoff=hbond_cutoff,
      hydrophobic_cutoff=hydrophobic_cutoff,
      salt_bridge_cutoff=salt_bridge_cutoff,
      pi_stack_cutoff=pi_stack_cutoff,
  )

  bit_vector = []
  bit_labels = []
  for r_id in sorted(residue_fingerprints.keys()):
    for itype in INTERACTION_TYPES:
      bit = 1 if residue_fingerprints[r_id][itype] > 0 else 0
      bit_vector.append(bit)
      bit_labels.append(f"{r_id}_{itype}")

  return {
      "complex": complex_path,
      "ligand_chain": ligand_chain,
      "total_contacts": len(detailed_interactions),
      "pocket_residue_count": len(residue_fingerprints),
      "pocket_residues": sorted(list(residue_fingerprints.keys())),
      "fingerprint_bits": bit_vector,
      "fingerprint_labels": bit_labels,
      "per_residue_plif": residue_fingerprints,
      "interactions": detailed_interactions,
  }


def write_plif_csv(
    interactions: List[Dict[str, Any]], output_path: str
) -> None:
  """Exports structured PLIF interactions table to CSV."""
  os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
  fieldnames = [
      "chain",
      "residue_number",
      "residue_name",
      "interaction_type",
      "distance_angstrom",
      "ligand_atom",
      "protein_atom",
  ]
  with open(output_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for row in interactions:
      writer.writerow({
          "chain": row.get("chain", ""),
          "residue_number": row.get("residue_number", ""),
          "residue_name": row.get("residue_name", ""),
          "interaction_type": row.get("interaction_type", ""),
          "distance_angstrom": (
              f"{float(row.get('distance_angstrom', 0.0)):.2f}"
          ),
          "ligand_atom": row.get("ligand_atom", ""),
          "protein_atom": row.get("protein_atom", ""),
      })


def compute_tanimoto(fp1: List[int], fp2: List[int]) -> float:
  """Computes binary Tanimoto similarity between two PLIF bit vectors."""
  a = sum(fp1)
  b = sum(fp2)
  c = sum(x and y for x, y in zip(fp1, fp2))
  denom = a + b - c
  return float(c / denom) if denom > 0 else 0.0


def main():
  parser = argparse.ArgumentParser(
      description="Compute Protein-Ligand Interaction Fingerprints (PLIF)."
  )
  parser.add_argument(
      "--complex",
      "--cif",
      dest="complex",
      required=True,
      help="Path to complex (.cif or .pdb).",
  )
  parser.add_argument(
      "--ligand-chain", default="B", help="Chain ID of ligand (default: 'B')."
  )
  parser.add_argument(
      "--ligand-sdf", default=None, help="Path to extracted ligand SDF."
  )
  parser.add_argument(
      "--output-json", default=None, help="Path to output JSON summary."
  )
  parser.add_argument(
      "--output-csv",
      default=None,
      help="Path to output CSV interactions table.",
  )
  parser.add_argument(
      "--reference-complex",
      default=None,
      help="Reference crystal complex for Tanimoto.",
  )
  parser.add_argument(
      "--reference-ligand-sdf",
      default=None,
      help=(
          "Ligand SDF for the reference complex. Defaults to --ligand-sdf when"
          " --reference-complex is the same file. Without it the reference"
          " carries no formal charges or aromatic rings, so its salt-bridge and"
          " pi-stacking bits are always zero and the Tanimoto is only a lower"
          " bound."
      ),
  )
  parser.add_argument(
      "--reference-ligand-chain",
      default=None,
      help=(
          "Ligand chain ID in the reference complex (defaults to"
          " --ligand-chain)."
      ),
  )
  args = parser.parse_args()

  res = compute_plif(
      complex_path=args.complex,
      ligand_chain=args.ligand_chain,
      ligand_sdf_path=args.ligand_sdf,
  )

  if args.reference_complex:
    ref_sdf = args.reference_ligand_sdf
    if ref_sdf is None and os.path.abspath(
        args.reference_complex
    ) == os.path.abspath(args.complex):
      # Same structure: reuse the ligand so a complex scores 1.000 against
      # itself instead of being penalised for bits the reference was never
      # given the chemistry to express.
      ref_sdf = args.ligand_sdf
    if ref_sdf is None:
      print(
          "WARNING: no --reference-ligand-sdf given; the reference will be "
          "scored without ligand formal charges or aromatic rings, so its "
          "salt-bridge and pi-stacking bits are forced to zero and "
          "tanimoto_to_reference is a LOWER BOUND, not a like-for-like "
          "comparison.",
          file=sys.stderr,
      )
    ref_res = compute_plif(
        complex_path=args.reference_complex,
        ligand_chain=args.reference_ligand_chain or args.ligand_chain,
        ligand_sdf_path=ref_sdf,
    )
    all_labels = sorted(
        list(set(res["pocket_residues"] + ref_res["pocket_residues"]))
    )
    v1, v2 = [], []
    for r in all_labels:
      for itype in INTERACTION_TYPES:
        v1.append(
            1 if res["per_residue_plif"].get(r, {}).get(itype, 0) > 0 else 0
        )
        v2.append(
            1 if ref_res["per_residue_plif"].get(r, {}).get(itype, 0) > 0 else 0
        )
    res["tanimoto_to_reference"] = round(compute_tanimoto(v1, v2), 3)
    res["reference_complex"] = args.reference_complex

  if args.output_csv:
    write_plif_csv(res["interactions"], args.output_csv)
    print(f"PLIF CSV saved to {args.output_csv}")

  if args.output_json:
    with open(args.output_json, "w", encoding="utf-8") as f:
      json.dump(res, f, indent=2)
    print(f"PLIF results saved to {args.output_json}")

  print(
      json.dumps(
          {
              "complex": res["complex"],
              "pocket_residue_count": res["pocket_residue_count"],
              "total_contacts": res["total_contacts"],
              "pocket_residues": res["pocket_residues"],
              "per_residue_plif": res["per_residue_plif"],
              "tanimoto_to_reference": res.get("tanimoto_to_reference"),
              "output_csv": args.output_csv,
          },
          indent=2,
      )
  )


if __name__ == "__main__":
  main()
