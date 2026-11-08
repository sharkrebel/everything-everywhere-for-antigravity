# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "gemmi",
#     "numpy",
#     "openmm",
#     "pdbfixer",
#     "polite-http",
# ]
# ///

"""Standalone pH-Dependent Protein/Receptor Protonation & PDBQT Preparation CLI.

Prepares experimental or predicted protein structures (PDB / mmCIF / RCSB ID)
at a target solution pH independently of molecular docking or MD workflows:

1. Strips crystallographic waters and buffer cosolvents (configurable via
   --keep-waters and --keep-resnames) and resolves alternate conformations
   (altLoc).
2. Repairs missing terminal or sidechain heavy atoms via PDBFixer (without
   fabricating missing loops) and assigns pH-dependent titratable residue
   states and explicit 3D hydrogen coordinates via
   `openmm.app.Modeller.addHydrogens(forcefield, pH=ph)` (`amber14-all.xml`).
3. Exports a protonated receptor `.pdb`, an optional AutoDock4/Vina rigid
   receptor `.pdbqt` (with verified polar `HD` hydrogen-bond donors), and a JSON
   summary of assigned titratable residue states and binding-pocket center.
"""

import argparse
import io
import json
import os
import re
import sys
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
from polite_http import http_client

try:
  import gemmi
except ImportError:
  gemmi = None

try:
  from openmm import app
  import pdbfixer
except ImportError:
  app = None
  pdbfixer = None

_RCSB_CLIENT = http_client.HttpClient("https://files.rcsb.org/", qps=5.0)

_DEFAULT_SOLVENT_RES = frozenset({
    "HOH",
    "WAT",
    "DOD",
    "DMS",
    "SO4",
    "GOL",
    "EDO",
    "ACT",
    "PO4",
    "CL",
    "NA",
    "K",
    "IOD",
    "PEG",
    "MPD",
    "TRS",
    "MES",
    "EPE",
})

_AROMATIC_RES_ATOMS = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "CD2", "CE1"},
    "HID": {"CG", "CD2", "CE1"},
    "HIE": {"CG", "CD2", "CE1"},
    "HIP": {"CG", "CD2", "CE1"},
}


def download_pdb_structure(pdb_id: str, output_path: str) -> str:
  """Downloads an experimental PDB structure from RCSB PDB via polite_http."""
  if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
    return output_path
  url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
  os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
  data = _RCSB_CLIENT.fetch_bytes(url, timeout=30.0)
  with open(output_path, "wb") as f:
    f.write(data)
  return output_path


def _clean_pdb_lines(
    pdb_path: str,
    ligand_resname: Optional[str] = None,
    keep_waters: bool = False,
    keep_hetatm: bool = False,
    keep_resnames: Optional[Set[str]] = None,
) -> Tuple[List[str], np.ndarray]:
  """Filters receptor ATOM/HETATM lines, resolves altLocs, and locates pocket center."""
  keep_set = {r.upper() for r in keep_resnames or set()}
  protein_lines: List[str] = []
  ligand_coords: List[List[float]] = []
  prot_coords: List[List[float]] = []
  seen_atom_keys: Set[Tuple[str, str, str, str]] = set()

  if pdb_path.lower().endswith((".cif", ".mmcif")):
    if gemmi is None:
      raise ImportError(
          "gemmi is required to parse .cif/.mmcif receptor structures in"
          " protonate_receptor.py."
      )
    st = gemmi.read_structure(pdb_path)
    pdb_text = st.make_minimal_pdb()
    raw_lines = pdb_text.splitlines(keepends=True)
  else:
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
      raw_lines = f.readlines()

  for line in raw_lines:
    if line.startswith("TER"):
      protein_lines.append(line if line.endswith("\n") else line + "\n")
      continue
    if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
      continue

    alt_loc = line[16:17]
    if alt_loc not in (" ", "A", "1"):
      continue

    atom_name = line[12:16].strip()
    res_name = line[17:20].strip().upper()
    chain_id = line[21:22]
    res_seq = line[22:27]
    atom_key = (chain_id, res_seq, res_name, atom_name)

    try:
      x = float(line[30:38])
      y = float(line[38:46])
      z = float(line[46:54])
    except ValueError:
      continue

    is_water = res_name in ("HOH", "WAT", "DOD")
    is_solvent = res_name in _DEFAULT_SOLVENT_RES

    if line.startswith("HETATM"):
      if ligand_resname and res_name == ligand_resname.upper():
        ligand_coords.append([x, y, z])
      elif not ligand_resname and not is_solvent:
        ligand_coords.append([x, y, z])

      should_keep = (
          (is_water and keep_waters)
          or (res_name in keep_set)
          or (keep_hetatm and not is_solvent)
      )
      if not should_keep:
        continue
    else:
      if is_water and not keep_waters:
        continue
      if (
          is_solvent
          and res_name not in keep_set
          and not (is_water and keep_waters)
      ):
        continue

    if atom_key in seen_atom_keys:
      continue
    seen_atom_keys.add(atom_key)

    clean_line = line[:16] + " " + line[17:]
    if not clean_line.endswith("\n"):
      clean_line += "\n"
    protein_lines.append(clean_line)
    prot_coords.append([x, y, z])

  if ligand_coords:
    center = np.mean(ligand_coords, axis=0)
  elif prot_coords:
    center = np.mean(prot_coords, axis=0)
  else:
    center = np.array([0.0, 0.0, 0.0])

  return protein_lines, center


def _unit_vec(v: np.ndarray) -> np.ndarray:
  norm = float(np.linalg.norm(v))
  if norm < 1e-6:
    return np.array([1.0, 0.0, 0.0], dtype=float)
  return v / norm


def _format_pdb_atom_line(
    serial: int,
    atom_name: str,
    res_name: str,
    chain_id: str,
    res_seq: str,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
  """Formats a standard 80-column PDB ATOM record."""
  name_field = (
      f" {atom_name:<3s}" if len(atom_name) < 4 else f"{atom_name[:4]:4s}"
  )
  chain_char = (chain_id or "A")[0]
  res_seq_field = f"{res_seq:>4s}"[:4]
  elem_field = f"{element.strip().upper():>2s}"[:2]
  return (
      f"ATOM  {serial:5d} {name_field:4s} {res_name:>3s} {chain_char:1s}"
      f"{res_seq_field:4s}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          "
      f"{elem_field:>2s}\n"
  )


def _geometric_protonate_pdb_lines(
    protein_lines: List[str],
    ph: float = 7.4,
) -> Tuple[List[str], Dict[str, int]]:
  """Fallback geometric polar-hydrogen placement for unparameterizable structures."""
  residues: List[Dict[str, Any]] = []
  res_map: Dict[Tuple[str, str], Dict[str, Any]] = {}

  for line in protein_lines:
    if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
      continue
    atom_name = line[12:16].strip()
    res_name = line[17:20].strip().upper()
    chain_id = line[21:22]
    res_seq = line[22:26].strip()
    try:
      coord = np.array(
          [float(line[30:38]), float(line[38:46]), float(line[46:54])],
          dtype=float,
      )
    except ValueError:
      continue
    elem = line[76:78].strip().upper() if len(line) >= 78 else ""
    if not elem:
      elem = re.sub(r"[^A-Za-z]", "", atom_name)[:1].upper()

    key = (chain_id, res_seq)
    if key not in res_map:
      rec = {
          "chain_id": chain_id,
          "res_seq": res_seq,
          "res_name": res_name,
          "atoms": {},
          "lines": [],
      }
      res_map[key] = rec
      residues.append(rec)
    res_map[key]["atoms"][atom_name] = coord
    res_map[key]["lines"].append((atom_name, elem, coord, line))

  state_counts: Dict[str, int] = {}
  out_lines: List[str] = []
  serial = 0

  for i, res in enumerate(residues):
    rname = res["res_name"]
    atoms = res["atoms"]
    chain_id = res["chain_id"]
    res_seq = res["res_seq"]

    assigned_name = rname
    if rname in ("HIS", "HID", "HIE", "HIP"):
      assigned_name = "HIP" if ph < 6.2 else "HIE"
      state_counts[assigned_name] = state_counts.get(assigned_name, 0) + 1
    elif rname in ("ASP", "ASH"):
      assigned_name = "ASH" if ph < 4.0 else "ASP"
      state_counts[assigned_name] = state_counts.get(assigned_name, 0) + 1
    elif rname in ("GLU", "GLH"):
      assigned_name = "GLH" if ph < 4.4 else "GLU"
      state_counts[assigned_name] = state_counts.get(assigned_name, 0) + 1
    elif rname in ("LYS", "LYN"):
      assigned_name = "LYN" if ph > 10.5 else "LYS"
      state_counts[assigned_name] = state_counts.get(assigned_name, 0) + 1

    has_any_h = any(el == "H" for _, el, _, _ in res["lines"])

    for _, _, _, orig_line in res["lines"]:
      serial += 1
      updated_line = (
          f"{orig_line[:6]}{serial:5d}{orig_line[11:17]}"
          f"{assigned_name:>3s}{orig_line[20:]}"
      )
      if not updated_line.endswith("\n"):
        updated_line += "\n"
      out_lines.append(updated_line)

    if has_any_h:
      continue

    if (
        rname != "PRO"
        and "N" in atoms
        and "CA" in atoms
        and i > 0
        and residues[i - 1]["chain_id"] == chain_id
        and "C" in residues[i - 1]["atoms"]
    ):
      n_pos = atoms["N"]
      ca_pos = atoms["CA"]
      prev_c = residues[i - 1]["atoms"]["C"]
      v_nh = _unit_vec(_unit_vec(n_pos - ca_pos) + _unit_vec(n_pos - prev_c))
      h_pos = n_pos + 1.01 * v_nh
      serial += 1
      out_lines.append(
          _format_pdb_atom_line(
              serial,
              "H",
              assigned_name,
              chain_id,
              res_seq,
              h_pos[0],
              h_pos[1],
              h_pos[2],
              "H",
          )
      )
    elif "N" in atoms and "CA" in atoms:
      n_pos = atoms["N"]
      ca_pos = atoms["CA"]
      h_pos = n_pos + 1.01 * _unit_vec(n_pos - ca_pos)
      serial += 1
      out_lines.append(
          _format_pdb_atom_line(
              serial,
              "H1",
              assigned_name,
              chain_id,
              res_seq,
              h_pos[0],
              h_pos[1],
              h_pos[2],
              "H",
          )
      )

    if rname == "SER" and "OG" in atoms and "CB" in atoms:
      h_pos = atoms["OG"] + 0.96 * _unit_vec(atoms["OG"] - atoms["CB"])
      serial += 1
      out_lines.append(
          _format_pdb_atom_line(
              serial,
              "HG",
              assigned_name,
              chain_id,
              res_seq,
              h_pos[0],
              h_pos[1],
              h_pos[2],
              "H",
          )
      )
    elif rname == "THR" and "OG1" in atoms and "CB" in atoms:
      h_pos = atoms["OG1"] + 0.96 * _unit_vec(atoms["OG1"] - atoms["CB"])
      serial += 1
      out_lines.append(
          _format_pdb_atom_line(
              serial,
              "HG1",
              assigned_name,
              chain_id,
              res_seq,
              h_pos[0],
              h_pos[1],
              h_pos[2],
              "H",
          )
      )
    elif rname == "TYR" and "OH" in atoms and "CZ" in atoms and ph < 10.0:
      h_pos = atoms["OH"] + 0.96 * _unit_vec(atoms["OH"] - atoms["CZ"])
      serial += 1
      out_lines.append(
          _format_pdb_atom_line(
              serial,
              "HH",
              assigned_name,
              chain_id,
              res_seq,
              h_pos[0],
              h_pos[1],
              h_pos[2],
              "H",
          )
      )

  out_lines.append("END\n")
  return out_lines, state_counts


def _build_openmm_modeller_from_pdb_lines(protein_lines: List[str]):
  """Builds an OpenMM Modeller, repairing missing terminal/sidechain heavy atoms via PDBFixer."""
  if app is None or pdbfixer is None:
    raise ImportError(
        "openmm and pdbfixer are required for OpenMM Amber14 receptor"
        " protonation."
    )
  pdb_text = "".join(protein_lines)
  fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(pdb_text))
  fixer.findMissingResidues()
  fixer.missingResidues = {}  # Do not fabricate unmodeled loop regions
  fixer.findNonstandardResidues()
  fixer.replaceNonstandardResidues()
  fixer.findMissingAtoms()
  fixer.addMissingAtoms()
  return app.Modeller(fixer.topology, fixer.positions)


def _openmm_protonate_pdb_lines(
    protein_lines: List[str],
    ph: float = 7.4,
) -> Tuple[List[str], Dict[str, int]]:
  """Protonates receptor PDB lines using OpenMM Modeller.addHydrogens(pH=ph)."""
  modeller = _build_openmm_modeller_from_pdb_lines(protein_lines)
  forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
  variants = modeller.addHydrogens(forcefield, pH=ph)

  state_counts: Dict[str, int] = {}
  for res, var in zip(modeller.topology.residues(), variants):
    label = var if var else res.name
    if res.name in (
        "HIS",
        "HID",
        "HIE",
        "HIP",
        "ASP",
        "ASH",
        "GLU",
        "GLH",
        "LYS",
        "LYN",
        "CYS",
        "CYM",
    ):
      state_counts[label] = state_counts.get(label, 0) + 1

  pdb_out = io.StringIO()
  app.PDBFile.writeFile(
      modeller.topology, modeller.positions, pdb_out, keepIds=True
  )
  out_lines = pdb_out.getvalue().splitlines(keepends=True)
  return out_lines, state_counts


def pdb_to_receptor_pdbqt(
    receptor_pdb: str,
    output_pdbqt: Optional[str] = None,
) -> Tuple[str, int]:
  """Converts a protonated receptor PDB into an AutoDock4/Vina rigid receptor PDBQT."""
  if receptor_pdb.endswith(".pdbqt") and os.path.exists(receptor_pdb):
    return receptor_pdb, 0
  if output_pdbqt is None:
    output_pdbqt = os.path.splitext(receptor_pdb)[0] + ".pdbqt"

  heavy_coords: List[Tuple[np.ndarray, str]] = []
  raw_records: List[Tuple[str, str, str, float, float, float, str]] = []

  with open(receptor_pdb, "r", encoding="utf-8", errors="replace") as f:
    for line in f:
      if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
        continue
      atom_name = line[12:16].strip()
      res_name = line[17:20].strip().upper()
      try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
      except ValueError:
        continue
      elem = line[76:78].strip().upper() if len(line) >= 78 else ""
      if not elem:
        elem = re.sub(r"[^A-Za-z]", "", atom_name)[:1].upper()
      raw_records.append((line, atom_name, res_name, x, y, z, elem))
      if elem != "H":
        heavy_coords.append((np.array([x, y, z], dtype=float), elem))

  out_lines: List[str] = []
  serial = 0
  polar_hd_count = 0

  for line, atom_name, res_name, x, y, z, elem in raw_records:
    if elem == "H":
      h_pt = np.array([x, y, z], dtype=float)
      min_d = 1.35
      nearest_elem = ""
      for hc, helem in heavy_coords:
        if (
            abs(hc[0] - x) > min_d
            or abs(hc[1] - y) > min_d
            or abs(hc[2] - z) > min_d
        ):
          continue
        d = float(np.linalg.norm(hc - h_pt))
        if d < min_d:
          min_d = d
          nearest_elem = helem
      if nearest_elem in ("N", "O", "S") or (
          not nearest_elem
          and not atom_name.startswith(
              ("HA", "HB", "HG2", "HG3", "HD2", "HD3", "HE2", "HE3")
          )
      ):
        ad4_type = "HD"
        polar_hd_count += 1
      else:
        continue
    elif elem == "C":
      if atom_name in _AROMATIC_RES_ATOMS.get(res_name, set()):
        ad4_type = "A"
      else:
        ad4_type = "C"
    elif elem == "O":
      ad4_type = "OA"
    elif elem == "N":
      if res_name in ("HIS", "HID", "HIE") and atom_name in ("ND1", "NE2"):
        ad4_type = "NA"
      else:
        ad4_type = "N"
    elif elem == "S":
      ad4_type = "SA"
    elif elem == "P":
      ad4_type = "P"
    elif elem in ("F", "CL", "BR", "I", "MG", "ZN", "FE", "CA", "MN"):
      ad4_type = elem.capitalize()
    else:
      ad4_type = "C"

    serial += 1
    prefix = f"{line[:6]}{serial:5d}{line[11:54]}"
    out_lines.append(f"{prefix}  1.00  0.00     0.000 {ad4_type:<2s}\n")

  os.makedirs(os.path.dirname(os.path.abspath(output_pdbqt)), exist_ok=True)
  with open(output_pdbqt, "w", encoding="utf-8") as f:
    f.writelines(out_lines)
  return output_pdbqt, polar_hd_count


def protonate_receptor_pdb(
    pdb_path: str,
    output_pdb: Optional[str] = None,
    output_pdbqt: Optional[str] = None,
    ph: float = 7.4,
    ligand_resname: Optional[str] = None,
    keep_waters: bool = False,
    keep_hetatm: bool = False,
    keep_resnames: Optional[Set[str]] = None,
    protonate: bool = True,
) -> Dict[str, Any]:
  """Cleans and protonates a receptor PDB at target pH using OpenMM + PDBFixer.

  Args:
    pdb_path: Input PDB/mmCIF file path or 4-letter RCSB PDB ID.
    output_pdb: Output protonated PDB path.
    output_pdbqt: Optional output AutoDock4/Vina rigid receptor PDBQT path.
    ph: Target pH for titratable residues (default: 7.4).
    ligand_resname: Optional 3-letter ligand residue ID for pocket centering.
    keep_waters: Retain crystallographic water molecules.
    keep_hetatm: Retain non-solvent HETATM cofactors/ions.
    keep_resnames: Explicit set of 3-letter residue names to retain.
    protonate: If False, only cleans solvent/altLocs without adding hydrogens.

  Returns:
    Dictionary containing `output_pdb`, `output_pdbqt`, `pocket_center`,
    `ph`, `backend`, `titratable_states`, `total_hydrogens`, and
    `polar_hd_donors`.

  Raises:
    FileNotFoundError: If `pdb_path` does not exist and is not a valid PDB ID.
    ValueError: If no valid protein `ATOM`/`HETATM` records are found.
  """
  resolved_in = pdb_path
  if len(pdb_path) == 4 and pdb_path.isalnum() and not os.path.exists(pdb_path):
    base_dir = os.path.dirname(os.path.abspath(output_pdb or "receptor.pdb"))
    resolved_in = os.path.join(base_dir, f"{pdb_path.lower()}_raw.pdb")
    download_pdb_structure(pdb_path, resolved_in)

  if not os.path.exists(resolved_in):
    raise FileNotFoundError(f"Receptor file not found: {pdb_path}")

  if output_pdb is None:
    stem = os.path.splitext(resolved_in)[0]
    output_pdb = f"{stem}_pH{ph:.1f}.pdb" if protonate else f"{stem}_clean.pdb"

  clean_lines, center = _clean_pdb_lines(
      resolved_in,
      ligand_resname=ligand_resname,
      keep_waters=keep_waters,
      keep_hetatm=keep_hetatm,
      keep_resnames=keep_resnames,
  )
  if not clean_lines:
    raise ValueError(
        f"No valid protein ATOM/HETATM records found in receptor: {pdb_path}"
    )

  backend = "unprotonated_clean_only"
  state_counts: Dict[str, int] = {}
  final_lines = clean_lines

  if protonate and clean_lines:
    try:
      final_lines, state_counts = _openmm_protonate_pdb_lines(
          clean_lines, ph=ph
      )
      backend = "openmm_amber14"
    except Exception:
      final_lines, state_counts = _geometric_protonate_pdb_lines(
          clean_lines, ph=ph
      )
      backend = "geometric_polar_h"

  os.makedirs(
      os.path.dirname(os.path.abspath(output_pdb)) or ".", exist_ok=True
  )
  with open(output_pdb, "w", encoding="utf-8") as f:
    f.writelines(final_lines)

  total_h = 0
  for line in final_lines:
    if line.startswith(("ATOM  ", "HETATM")):
      elem = line[76:78].strip().upper() if len(line) >= 78 else ""
      if not elem:
        elem = re.sub(r"[^A-Za-z]", "", line[12:16].strip())[:1].upper()
      if elem == "H":
        total_h += 1

  pdbqt_out, polar_hd = (None, 0)
  if output_pdbqt:
    pdbqt_out, polar_hd = pdb_to_receptor_pdbqt(output_pdb, output_pdbqt)
  else:
    tmp_pdbqt = output_pdb + ".tmp.pdbqt"
    try:
      _, polar_hd = pdb_to_receptor_pdbqt(output_pdb, tmp_pdbqt)
    finally:
      if os.path.exists(tmp_pdbqt):
        os.remove(tmp_pdbqt)

  return {
      "input_receptor": pdb_path,
      "output_pdb": os.path.abspath(output_pdb),
      "output_pdbqt": os.path.abspath(pdbqt_out) if pdbqt_out else None,
      "ph": ph,
      "protonated": protonate,
      "backend": backend,
      "pocket_center": [round(float(c), 3) for c in center],
      "titratable_states": state_counts,
      "total_hydrogens": total_h,
      "polar_hd_donors": polar_hd,
  }


def main() -> None:
  parser = argparse.ArgumentParser(
      description=(
          "Standalone pH-dependent protein/receptor protonation via OpenMM"
          " Modeller + PDBFixer."
      )
  )
  parser.add_argument(
      "--receptor",
      "--input",
      "--pdb",
      dest="receptor",
      type=str,
      required=True,
      help="Input receptor .pdb/.cif file path or 4-letter RCSB PDB ID.",
  )
  parser.add_argument(
      "--ph",
      type=float,
      default=7.4,
      help="Target solution pH for titratable residues (default: 7.4).",
  )
  parser.add_argument(
      "--ligand-resname",
      "--ligand_resname",
      dest="ligand_resname",
      type=str,
      default=None,
      help=(
          "Optional 3-letter co-crystal ligand residue name for pocket"
          " centering."
      ),
  )
  parser.add_argument(
      "--keep-waters",
      "--keep_waters",
      dest="keep_waters",
      action="store_true",
      help="Retain crystallographic water molecules (HOH/WAT).",
  )
  parser.add_argument(
      "--keep-hetatm",
      "--keep_hetatm",
      dest="keep_hetatm",
      action="store_true",
      help="Retain non-solvent HETATM cofactors/ions.",
  )
  parser.add_argument(
      "--keep-resnames",
      "--keep_resnames",
      dest="keep_resnames",
      type=str,
      default=None,
      help=(
          "Comma-separated 3-letter HETATM residue names to retain (e.g."
          " 'HEM,ZN,MG')."
      ),
  )
  parser.add_argument(
      "--output-pdb",
      "--output_pdb",
      "--output",
      dest="output_pdb",
      type=str,
      default=None,
      help="Output protonated PDB file path.",
  )
  parser.add_argument(
      "--output-pdbqt",
      "--output_pdbqt",
      dest="output_pdbqt",
      type=str,
      default=None,
      help="Optional output AutoDock4/Vina rigid receptor PDBQT file path.",
  )
  parser.add_argument(
      "--summary-json",
      "--summary_json",
      "--output-json",
      "--output_json",
      dest="summary_json",
      type=str,
      default=None,
      help="Optional output JSON summary report path.",
  )
  args = parser.parse_args()

  keep_res = (
      {r.strip().upper() for r in args.keep_resnames.split(",") if r.strip()}
      if args.keep_resnames
      else None
  )

  try:
    summary = protonate_receptor_pdb(
        pdb_path=args.receptor,
        output_pdb=args.output_pdb,
        output_pdbqt=args.output_pdbqt,
        ph=args.ph,
        ligand_resname=args.ligand_resname,
        keep_waters=args.keep_waters,
        keep_hetatm=args.keep_hetatm,
        keep_resnames=keep_res,
        protonate=True,
    )
  except (FileNotFoundError, ValueError, RuntimeError) as exc:
    print(f"Error: {exc}", file=sys.stderr)
    sys.exit(1)

  if args.summary_json:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.summary_json)) or ".",
        exist_ok=True,
    )
    with open(args.summary_json, "w", encoding="utf-8") as f:
      json.dump(summary, f, indent=2)

  print(f"Prepared receptor at pH {args.ph:.1f} ({summary['backend']}):")
  print(f"  - Protonated PDB:     {summary['output_pdb']}")
  if summary["output_pdbqt"]:
    print(f"  - Receptor PDBQT:     {summary['output_pdbqt']}")
  print(f"  - Total Hydrogens:    {summary['total_hydrogens']}")
  print(f"  - Polar HD Donors:    {summary['polar_hd_donors']}")
  print(f"  - Titratable States:  {summary['titratable_states']}")
  print(f"  - Pocket Center (A):  {summary['pocket_center']}")


if __name__ == "__main__":
  main()
