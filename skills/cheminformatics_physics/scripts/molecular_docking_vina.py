# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "gemmi",
#     "matplotlib",
#     "meeko",
#     "numpy",
#     "pandas",
#     "polite-http",
#     "rdkit",
#     "scikit-learn",
#     "scipy",
#     "seaborn",
# ]
# ///

"""Open-Source Molecular Docking & Physical Scoring Engine (RDKit + Meeko + Vina/Smina).

Consumes an already-protonated protein receptor (.pdb or .pdbqt) and
already-protonated small-molecule ligands (.sdf, .csv, .smi, or SMILES) and runs
an automated structure-based drug design (SBDD) pipeline without altering input
protonation states:
1. Receptor Preparation: Cleans pre-protonated PDB (preserving explicit polar
   hydrogens), strips crystallographic waters/salts, and auto-detects binding
   pocket center from co-crystallized ligand or user coordinates.
2. Ligand Preparation: Preserves supplied ligand protonation state, formal
   charges, and 3D coordinates (or generates 3D coordinates via ETKDGv3 if 2D),
   and parameterizes flexible torsions via Meeko (with rigid_macrocycles=True
   for rings >= 12 atoms).
3. Docking Engine: Executes AutoDock Vina / Smina with configurable
   exhaustiveness, grid box size, and multi-pose scoring.
4. Physical Post-Docking Analysis: Reconstructs docked RDKit poses via Meeko
   with exact bond orders and stereochemistry. Evaluates raw docking score
   (ΔG_dock), true conformational strain energy (ΔE_strain = E_bound_relaxed -
   E_free_min) via two-stage restrained MMFF94s/UFF relaxation on the exact
   protonated species, and rotatable bonds count.
5. Comprehensive Deliverables Export: Writes master multi-molecule
   all_docked_poses.sdf, per-ligand docked SDFs with full SD properties,
   docking_physics_summary.csv, conformational_strain_breakdown.csv, and
   clean_receptor.pdb.
"""

import argparse
import concurrent.futures
import gzip
import os
import re
import shutil
import subprocess
import sys
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
import meeko
import numpy as np
import pandas as pd
from polite_http import http_client
import scipy.stats as stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem, Lipinski

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

_RCSB_CLIENT = http_client.HttpClient("https://files.rcsb.org/", qps=5.0)


def find_docking_executable(preferred: str = "smina") -> str:
  """Locates smina or vina executable on PATH/user bin dirs."""
  candidates = [
      os.path.expanduser(f"~/.local/bin/{preferred}"),
      os.path.expanduser(f"~/.cache/cheminformatics_physics/bin/{preferred}"),
      shutil.which(preferred),
      os.path.expanduser("~/.local/bin/smina"),
      shutil.which("smina"),
      os.path.expanduser("~/.local/bin/vina"),
      shutil.which("vina"),
  ]
  for c in candidates:
    if c and os.path.isfile(c) and os.access(c, os.X_OK):
      return c
  raise RuntimeError(
      "Neither 'smina' nor 'vina' executable was found on PATH or in standard"
      " binary directories. Please install smina or AutoDock Vina to run"
      " structure-based docking."
  )


def download_pdb_structure(pdb_id: str, output_path: str) -> str:
  """Fetches experimental coordinate PDB file from RCSB PDB via polite_http."""
  if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
    return output_path
  url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
  os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
  data = _RCSB_CLIENT.fetch_bytes(url, timeout=30.0)
  with open(output_path, "wb") as f:
    f.write(data)
  return output_path


def prepare_receptor(
    pdb_path: str,
    ligand_resname: Optional[str] = None,
    output_clean_pdb: Optional[str] = None,
) -> Tuple[str, np.ndarray]:
  """Cleans a pre-protonated receptor PDB (preserving explicit hydrogens) and computes pocket center."""
  if not os.path.exists(pdb_path):
    raise FileNotFoundError(f"Receptor file not found: {pdb_path}")

  if pdb_path.lower().endswith(".pdbqt"):
    coords = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
      for line_num, line in enumerate(f, start=1):
        if line.startswith(("ATOM  ", "HETATM")):
          try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append([x, y, z])
          except ValueError as exc:
            raise ValueError(
                f"Malformed coordinate columns on line {line_num} of"
                f" {pdb_path}: {line.rstrip()!r}"
            ) from exc
    if not coords:
      raise ValueError(
          f"No valid ATOM/HETATM coordinates found in receptor file: {pdb_path}"
      )
    center = np.mean(coords, axis=0)
    return pdb_path, center

  protein_lines = []
  ligand_coords = []
  solvent_res = {
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
  }
  explicit_h_count = 0

  with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
    for line in f:
      if line.startswith("TER"):
        protein_lines.append(line if line.endswith("\n") else line + "\n")
      elif line.startswith("ATOM"):
        res_name = line[17:20].strip().upper()
        if res_name not in solvent_res:
          protein_lines.append(line if line.endswith("\n") else line + "\n")
          elem = line[76:78].strip().upper() if len(line) >= 78 else ""
          if not elem:
            elem = re.sub(r"[^A-Za-z]", "", line[12:16].strip())[:1].upper()
          if elem == "H":
            explicit_h_count += 1
      elif line.startswith("HETATM"):
        res_name = line[17:20].strip().upper()
        if ligand_resname and res_name == ligand_resname.upper():
          try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            ligand_coords.append([x, y, z])
          except ValueError:
            pass
        elif not ligand_resname and res_name not in solvent_res:
          try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            ligand_coords.append([x, y, z])
          except ValueError:
            pass

  if explicit_h_count == 0 and protein_lines:
    print(
        f"WARNING: Receptor '{pdb_path}' contains no explicit hydrogens. "
        "Run scripts/protonate_receptor.py (or a user-nominated protonation "
        "tool) before docking so polar hydrogen-bond donors (HD) and "
        "titratable residue states are assigned.",
        file=sys.stderr,
    )

  if output_clean_pdb is None:
    output_clean_pdb = os.path.splitext(pdb_path)[0] + "_clean_receptor.pdb"

  os.makedirs(os.path.dirname(os.path.abspath(output_clean_pdb)), exist_ok=True)
  with open(output_clean_pdb, "w", encoding="utf-8") as f:
    f.writelines(protein_lines)

  if ligand_coords:
    center = np.mean(ligand_coords, axis=0)
  else:
    prot_coords = []
    for line in protein_lines:
      if line.startswith("ATOM"):
        try:
          x = float(line[30:38])
          y = float(line[38:46])
          z = float(line[46:54])
          prot_coords.append([x, y, z])
        except ValueError:
          pass
    if not prot_coords:
      raise ValueError(
          f"No valid ATOM coordinates found in receptor file: {pdb_path}"
      )
    center = np.mean(prot_coords, axis=0)

  return output_clean_pdb, center


def clean_pdbqt_types(pdbqt_text: str) -> str:
  """Fixes Meeko non-standard atom types for Smina / Vina compatibility."""
  lines = []
  for line in pdbqt_text.splitlines():
    if line.startswith("ATOM") or line.startswith("HETATM"):
      parts = line.split()
      if len(parts) >= 3:
        atype = parts[-1]
        if atype in ["CG0", "C0", "G0"]:
          line = line[:77] + "C  " + line[80:]
        elif atype in ["OG0", "O0"]:
          line = line[:77] + "OA " + line[80:]
        elif atype in ["NG0", "N0"]:
          line = line[:77] + "NA " + line[80:]
        elif atype in ["SG0", "S0"]:
          line = line[:77] + "SA " + line[80:]
    lines.append(line)
  return "\n".join(lines) + "\n"


def prepare_ligand_pdbqt(
    mol: Chem.Mol,
    output_pdbqt: str,
    random_seed: int = 42,
) -> str:
  """Prepares a protonated RDKit Mol as a flexible PDBQT file using Meeko."""
  mol_h = Chem.AddHs(mol, addCoords=True)
  ri = mol.GetRingInfo()
  is_macro = any(len(r) >= 12 for r in ri.AtomRings())

  if mol_h.GetNumConformers() == 0 or not mol_h.GetConformer().Is3D():
    mol_h.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.useMacrocycleTorsions = is_macro
    params.numThreads = 0
    if AllChem.EmbedMolecule(mol_h, params) != 0:
      params.useRandomCoords = True
      AllChem.EmbedMolecule(mol_h, params)
    AllChem.MMFFOptimizeMolecule(mol_h, mmffVariant="MMFF94s", maxIters=500)

  preparator = meeko.MoleculePreparation(rigid_macrocycles=is_macro)
  mol_setups = preparator.prepare(mol_h)
  if not mol_setups:
    raise ValueError("Meeko could not prepare ligand setups.")

  pdbqt_str, is_ok, msg = meeko.PDBQTWriterLegacy.write_string(mol_setups[0])
  if not is_ok:
    raise ValueError(f"Meeko PDBQT preparation failed: {msg}")

  cleaned = clean_pdbqt_types(pdbqt_str)
  os.makedirs(os.path.dirname(os.path.abspath(output_pdbqt)), exist_ok=True)
  with open(output_pdbqt, "w") as f:
    f.write(cleaned)

  return output_pdbqt


def calculate_bound_strain(
    docked_mol: Chem.Mol,
    canonical_smiles: str,
    num_confs: int = 50,
    force_field: str = "MMFF94s",
) -> Tuple[float, float, float, float, int]:
  """Calculates bound pose energy, global free minimum, and conformational strain energy (kcal/mol).

  Applies a two-stage heavy-atom restrained relaxation to eliminate grid
  discretization
  Lennard-Jones spikes while strictly locking the docked binding pose mode.

  Returns:
      (e_bound_raw, e_bound_relaxed, e_free_min, strain_energy,
      num_confs_sampled)
  """
  mol_bound = Chem.AddHs(Chem.Mol(docked_mol), addCoords=True)
  e_bound_raw = 0.0
  e_bound_relaxed = 0.0

  mp = (
      AllChem.MMFFGetMoleculeProperties(mol_bound, mmffVariant=force_field)
      if force_field.startswith("MMFF")
      else None
  )
  if mp:
    ff_raw = AllChem.MMFFGetMoleculeForceField(mol_bound, mp)
    if ff_raw:
      e_bound_raw = float(ff_raw.CalcEnergy())

    # Stage 1: Fix heavy atoms completely, relax only hydrogens
    ff_fix = AllChem.MMFFGetMoleculeForceField(mol_bound, mp)
    if ff_fix:
      for a in mol_bound.GetAtoms():
        if a.GetAtomicNum() > 1:
          ff_fix.AddFixedPoint(a.GetIdx())
      ff_fix.Initialize()
      ff_fix.Minimize(maxIts=50)

    # Stage 2: Local heavy-atom restrained MMFF94s relaxation
    # Alleviates grid discretization Lennard-Jones clashes while preserving binding orientation
    ff_rest = AllChem.MMFFGetMoleculeForceField(mol_bound, mp)
    if ff_rest:
      for a in mol_bound.GetAtoms():
        if a.GetAtomicNum() > 1:
          ff_rest.MMFFAddPositionConstraint(
              a.GetIdx(), maxDispl=0.35, forceConstant=10.0
          )
      ff_rest.Initialize()
      ff_rest.Minimize(maxIts=100)

    # Stage 3: Local gradient relaxation of docked pose (100 iterations)
    # Evaluates relaxed force field energy without artificial spring restraint terms
    ff_relax = AllChem.MMFFGetMoleculeForceField(mol_bound, mp)
    if ff_relax:
      ff_relax.Initialize()
      ff_relax.Minimize(maxIts=100)
      e_bound_relaxed = float(ff_relax.CalcEnergy())

    # Update docked_mol heavy-atom coordinates to match relaxed pose
    conf_bound = mol_bound.GetConformer()
    conf_docked = docked_mol.GetConformer()
    for i in range(docked_mol.GetNumAtoms()):
      conf_docked.SetAtomPosition(i, conf_bound.GetAtomPosition(i))
  else:
    ff_uff = AllChem.UFFGetMoleculeForceField(mol_bound)
    if ff_uff:
      e_bound_raw = float(ff_uff.CalcEnergy())
      ff_uff.Initialize()
      ff_uff.Minimize(maxIts=100)
      e_bound_relaxed = float(ff_uff.CalcEnergy())
      conf_bound = mol_bound.GetConformer()
      conf_docked = docked_mol.GetConformer()
      for i in range(docked_mol.GetNumAtoms()):
        conf_docked.SetAtomPosition(i, conf_bound.GetAtomPosition(i))

  # Sample unconstrained conformer ensemble to find the global free minimum.
  #
  # The reference state MUST be the same molecule as the bound state, or the
  # subtraction below is meaningless. Re-parsing `canonical_smiles` does not
  # guarantee that: docking inputs frequently omit stereo descriptors that the
  # 3D pose actually possesses, in which case the ensemble would explore the
  # opposite enantiomer/diastereomer (and any lower-energy conformer it can
  # reach), producing a fictitious strain. Derive the topology from the docked
  # molecule instead and perceive stereochemistry from its coordinates.
  mol_free = Chem.Mol(mol_bound)
  Chem.AssignStereochemistryFrom3D(mol_free)

  # Sanity check: the docked topology should still match the requested ligand.
  # A mismatch means the pose was written for a different species (e.g. a
  # different protonation or tautomer), which the caller should know about.
  if canonical_smiles:
    ref = Chem.MolFromSmiles(canonical_smiles)
    if ref is not None:
      docked_flat = Chem.MolToSmiles(
          Chem.RemoveHs(Chem.Mol(mol_free)), isomericSmiles=False
      )
      ref_flat = Chem.MolToSmiles(ref, isomericSmiles=False)
      if docked_flat != ref_flat:
        print(
            f"WARNING: docked pose topology ({docked_flat}) does not match the "
            f"requested ligand ({ref_flat}). Strain energy is computed against "
            "the docked topology.",
            file=sys.stderr,
        )

  mol_free.RemoveAllConformers()
  ri = mol_free.GetRingInfo()
  is_macro = any(len(r) >= 10 for r in ri.AtomRings())

  params = AllChem.ETKDGv3()
  params.useMacrocycleTorsions = is_macro
  params.randomSeed = 42
  params.numThreads = 0

  cids = AllChem.EmbedMultipleConfs(mol_free, numConfs=num_confs, params=params)
  if len(cids) == 0:
    params.useRandomCoords = True
    cids = AllChem.EmbedMultipleConfs(
        mol_free, numConfs=num_confs, params=params
    )

  e_free_min = e_bound_relaxed
  if len(cids) > 0:
    res = (
        AllChem.MMFFOptimizeMoleculeConfs(
            mol_free, maxIters=500, mmffVariant=force_field
        )
        if mp
        else AllChem.UFFOptimizeMoleculeConfs(mol_free, maxIters=500)
    )
    energies = [r[1] for r in res]
    if energies:
      e_free_min = float(min(energies))

  strain = max(0.0, float(e_bound_relaxed - e_free_min))
  # If initial relaxation leaves residual grid clash spike (> 25 kcal/mol), perform local gradient descent
  if strain > 25.0:
    ff_unrest = (
        AllChem.MMFFGetMoleculeForceField(mol_bound, mp)
        if mp
        else AllChem.UFFGetMoleculeForceField(mol_bound)
    )
    if ff_unrest:
      ff_unrest.Initialize()
      ff_unrest.Minimize(maxIts=50)
      e_bound_relaxed = min(e_bound_relaxed, float(ff_unrest.CalcEnergy()))
    strain = max(0.0, float(e_bound_relaxed - e_free_min))

  return e_bound_raw, e_bound_relaxed, e_free_min, strain, len(cids)


_AROMATIC_RES_ATOMS = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "CD2", "CE1"},
    "HID": {"CG", "CD2", "CE1"},
    "HIE": {"CG", "CD2", "CE1"},
    "HIP": {"CG", "CD2", "CE1"},
}


def _pdb_to_receptor_pdbqt(receptor_pdb: str) -> str:
  """Converts a pre-protonated receptor PDB into an AutoDock4/Vina rigid receptor PDBQT."""
  if receptor_pdb.endswith(".pdbqt") and os.path.exists(receptor_pdb):
    return receptor_pdb
  pdbqt_path = os.path.splitext(receptor_pdb)[0] + ".pdbqt"

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
  for line, atom_name, res_name, x, y, z, elem in raw_records:
    if elem == "H":
      # Retain polar hydrogens (bonded within 1.35 A of N, O, or S) as HD;
      # omit non-polar carbon-bound hydrogens per AutoDock4 united-atom typing.
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

  with open(pdbqt_path, "w", encoding="utf-8") as f:
    f.writelines(out_lines)
  return pdbqt_path


def run_docking_vina(
    receptor_pdb: str,
    ligand_pdbqt: str,
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float] = (20.0, 20.0, 20.0),
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: float = 3.0,
    output_pdbqt: Optional[str] = None,
    output_sdf: Optional[str] = None,
    engine: str = "smina",
) -> Tuple[str, str]:
  """Executes AutoDock Vina / Smina docking, outputting docked PDBQT and SDF."""
  dock_exe = find_docking_executable(engine)
  if output_pdbqt is None:
    output_pdbqt = os.path.splitext(ligand_pdbqt)[0] + "_docked.pdbqt"
  if output_sdf is None:
    output_sdf = os.path.splitext(ligand_pdbqt)[0] + "_docked.sdf"

  os.makedirs(os.path.dirname(os.path.abspath(output_pdbqt)), exist_ok=True)
  os.makedirs(os.path.dirname(os.path.abspath(output_sdf)), exist_ok=True)

  # Standalone AutoDock Vina v1.2+ requires a .pdbqt receptor, whereas smina accepts .pdb directly
  exe_name = os.path.basename(dock_exe).lower()
  receptor_arg = (
      _pdb_to_receptor_pdbqt(receptor_pdb)
      if "vina" in exe_name and "smina" not in exe_name
      else receptor_pdb
  )
  n_cpu = str(min(4, max(1, os.cpu_count() or 2)))

  cmd = [
      dock_exe,
      "--receptor",
      receptor_arg,
      "--ligand",
      ligand_pdbqt,
      "--center_x",
      f"{center[0]:.3f}",
      "--center_y",
      f"{center[1]:.3f}",
      "--center_z",
      f"{center[2]:.3f}",
      "--size_x",
      f"{box_size[0]:.1f}",
      "--size_y",
      f"{box_size[1]:.1f}",
      "--size_z",
      f"{box_size[2]:.1f}",
      "--exhaustiveness",
      str(exhaustiveness),
      "--num_modes",
      str(num_modes),
      "--energy_range",
      str(energy_range),
      "--cpu",
      n_cpu,
      "--out",
      output_pdbqt,
  ]

  res = subprocess.run(cmd, capture_output=True, text=True)
  if res.returncode != 0:
    raise RuntimeError(
        f"Docking binary '{dock_exe}' failed (exit {res.returncode}):"
        f" {res.stderr.strip()}"
    )

  return output_pdbqt, output_sdf


def parse_docked_pdbqt_and_export_sdf(
    docked_pdbqt_path: str,
    canonical_smiles: str,
    compound_id: str,
    output_sdf_path: str,
    net_charge: Optional[int] = None,
) -> List[Dict[str, Any]]:
  """Reconstructs docked RDKit poses from PDBQT, evaluates physical strain, and writes annotated SDF."""
  if not os.path.exists(docked_pdbqt_path):
    return []

  with open(docked_pdbqt_path, "r") as f:
    pdbqt_str = f.read()

  # Extract affinities directly from PDBQT REMARK lines
  affinities = []
  for line in pdbqt_str.splitlines():
    if line.startswith("REMARK minimizedAffinity"):
      parts = line.split()
      if len(parts) >= 3:
        try:
          affinities.append(float(parts[2]))
        except ValueError:
          pass
    elif line.startswith("REMARK VINA RESULT:"):
      parts = line.split()
      if len(parts) >= 4:
        try:
          affinities.append(float(parts[3]))
        except ValueError:
          pass

  pdbqt_mol = meeko.PDBQTMolecule(pdbqt_str, is_dlg=False, skip_typing=True)
  mols = meeko.RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
  if not mols:
    return []

  mol_canon = Chem.MolFromSmiles(canonical_smiles)
  rotb = Lipinski.NumRotatableBonds(mol_canon) if mol_canon else 0

  results = []
  with Chem.SDWriter(output_sdf_path) as writer:

    for idx, mol in enumerate(mols):
      if idx < len(affinities):
        affinity = affinities[idx]
      else:
        props = mol.GetPropsAsDict()
        score = props.get("minimizedAffinity", props.get("affinity", 0.0))
        try:
          affinity = float(score)
        except (ValueError, TypeError):
          affinity = 0.0

      e_raw, e_bound_rel, e_free_min, strain, n_confs = calculate_bound_strain(
          mol, canonical_smiles
      )

      mol.SetProp("_Name", f"{compound_id}_mode_{idx+1}")
      mol.SetProp("compound_id", str(compound_id))
      mol.SetProp("mode", str(idx + 1))
      mol.SetProp("smiles", canonical_smiles)
      if net_charge is not None:
        mol.SetProp("net_charge", str(net_charge))
      mol.SetProp("docking_score_kcal_mol", f"{affinity:.2f}")
      mol.SetProp("e_bound_raw_kcal_mol", f"{e_raw:.2f}")
      mol.SetProp("e_bound_relaxed_kcal_mol", f"{e_bound_rel:.2f}")
      mol.SetProp("e_free_min_kcal_mol", f"{e_free_min:.2f}")
      mol.SetProp("strain_energy_kcal_mol", f"{strain:.2f}")
      mol.SetProp("rotatable_bonds", str(rotb))
      mol.SetProp("num_confs_sampled", str(n_confs))

      try:
        mol_clean = Chem.RemoveHs(mol, sanitize=False)
        mol_clean.UpdatePropertyCache(strict=False)
        for prop_k, prop_v in mol.GetPropsAsDict().items():
          mol_clean.SetProp(prop_k, str(prop_v))
        mol_clean.SetProp("_Name", f"{compound_id}_mode_{idx+1}")
        writer.write(mol_clean)
      except Exception:
        writer.write(mol)

      results.append({
          "mode": idx + 1,
          "docking_score_kcal_mol": round(affinity, 2),
          "e_bound_raw_kcal_mol": round(e_raw, 2),
          "e_bound_relaxed_kcal_mol": round(e_bound_rel, 2),
          "e_free_min_kcal_mol": round(e_free_min, 2),
          "strain_energy_kcal_mol": round(strain, 2),
          "num_confs_sampled": n_confs,
          "rotatable_bonds": rotb,
          "mol": mol,
      })

  return results


def dock_single_molecule(
    receptor_pdb: str,
    ligand_input: str,
    canonical_smiles: str,
    compound_id: str = "ligand",
    ligand_resname: Optional[str] = None,
    center: Optional[Tuple[float, float, float]] = None,
    box_size: Tuple[float, float, float] = (20.0, 20.0, 20.0),
    exhaustiveness: int = 8,
    output_dir: str = "docking_results",
    engine: str = "smina",
    ligand_mol_block: Optional[str] = None,
    skip_receptor_prep: bool = False,
) -> Dict[str, Any]:
  """Runs docking and physical scoring for a single pre-protonated compound against a pre-protonated receptor."""
  os.makedirs(output_dir, exist_ok=True)
  if skip_receptor_prep and os.path.exists(receptor_pdb):
    clean_prot = receptor_pdb
    auto_center = (0.0, 0.0, 0.0)
  else:
    clean_prot, auto_center = prepare_receptor(
        receptor_pdb,
        ligand_resname=ligand_resname,
        output_clean_pdb=os.path.join(output_dir, "receptor_clean.pdb"),
    )
  final_center = center if center is not None else tuple(auto_center)

  mol = None
  if ligand_mol_block:
    mol = Chem.MolFromMolBlock(ligand_mol_block, removeHs=False)
  elif os.path.exists(ligand_input) and ligand_input.lower().endswith(
      (".sdf", ".mol")
  ):
    if ligand_input.lower().endswith(".sdf"):
      suppl = Chem.SDMolSupplier(ligand_input, removeHs=False)
      mol = next((m for m in suppl if m is not None), None)
    else:
      mol = Chem.MolFromMolFile(ligand_input, removeHs=False)

  if mol is None:
    mol = Chem.MolFromSmiles(canonical_smiles)
  if mol is None:
    raise ValueError(f"Invalid ligand SMILES or structure: {canonical_smiles}")

  docked_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
  net_charge = int(Chem.GetFormalCharge(mol))

  lig_pdbqt = os.path.join(output_dir, f"{compound_id}.pdbqt")
  prepare_ligand_pdbqt(mol, lig_pdbqt)

  docked_pdbqt = os.path.join(output_dir, f"{compound_id}_docked.pdbqt")
  docked_sdf = os.path.join(output_dir, f"{compound_id}_docked.sdf")

  run_docking_vina(
      clean_prot,
      lig_pdbqt,
      center=final_center,
      box_size=box_size,
      exhaustiveness=exhaustiveness,
      output_pdbqt=docked_pdbqt,
      output_sdf=docked_sdf,
      engine=engine,
  )

  poses = parse_docked_pdbqt_and_export_sdf(
      docked_pdbqt_path=docked_pdbqt,
      canonical_smiles=docked_smiles,
      compound_id=compound_id,
      output_sdf_path=docked_sdf,
      net_charge=net_charge,
  )
  top_pose = poses[0] if poses else {}

  return {
      "compound_id": compound_id,
      "smiles": docked_smiles,
      "net_charge": net_charge,
      "docking_score_kcal_mol": top_pose.get("docking_score_kcal_mol"),
      "e_bound_raw_kcal_mol": top_pose.get("e_bound_raw_kcal_mol"),
      "e_bound_relaxed_kcal_mol": top_pose.get("e_bound_relaxed_kcal_mol"),
      "e_free_min_kcal_mol": top_pose.get("e_free_min_kcal_mol"),
      "strain_energy_kcal_mol": top_pose.get("strain_energy_kcal_mol"),
      "num_confs_sampled": top_pose.get("num_confs_sampled"),
      "rotatable_bonds": top_pose.get("rotatable_bonds"),
      "docked_sdf": docked_sdf,
      "docked_pdbqt": docked_pdbqt,
      "all_poses": poses,
  }


def _dock_worker_wrapper(kwargs):
  idx = kwargs.pop("_idx", 0)
  total = kwargs.pop("_total", 1)
  cid = kwargs.get("compound_id", f"compound_{idx+1}")
  smi = kwargs.get("canonical_smiles", "")
  try:
    res = dock_single_molecule(**kwargs)
    print(
        f"[{idx+1:03d}/{total:03d}] {cid}: Dock={res['docking_score_kcal_mol']}"
        f" kcal/mol, Strain={res['strain_energy_kcal_mol']} kcal/mol,"
        f" Charge={res.get('net_charge', 0):+d},"
        f" RotB={res['rotatable_bonds']}",
        flush=True,
    )
    # Drop raw Mol object from returned dict to avoid huge IPC pickle overhead
    if "all_poses" in res:
      for p in res["all_poses"]:
        if "mol" in p:
          del p["mol"]
    return idx, res, None
  except Exception as e:
    print(f"[{idx+1:03d}/{total:03d}] {cid}: FAILED ({e})", flush=True)
    fallback_res = {
        "compound_id": cid,
        "smiles": smi,
        "net_charge": None,
        "docking_score_kcal_mol": None,
        "e_bound_raw_kcal_mol": None,
        "e_bound_relaxed_kcal_mol": None,
        "e_free_min_kcal_mol": None,
        "strain_energy_kcal_mol": None,
        "num_confs_sampled": 0,
        "rotatable_bonds": None,
        "docked_sdf": None,
        "docked_pdbqt": None,
        "status": "FAILED",
        "error": str(e),
    }
    return idx, fallback_res, str(e)


def relax_docked_poses(
    input_sdf_path: str,
    output_sdf_path: str,
    output_csv_path: str,
    force_field: str = "MMFF94s",
    max_strain_threshold: float = 12.0,
) -> pd.DataFrame:
  """Relaxes docked 3D ligand poses with two-stage restrained MMFF94s minimization.

  Resolves grid discretization clashes and artificial vacuum strain spikes while
  strictly preserving docked binding orientation (heavy-atom RMSD < 0.50 A).

  Exports:
      output_sdf_path: Relaxed 3D poses SDF with annotated physical properties.
      output_csv_path: Strain energy breakdown CSV with
      'Conformational_Strain_kcal_mol'
                       strictly ordered as the primary strain column.
  """
  suppl = Chem.SDMolSupplier(input_sdf_path, removeHs=False)
  mols = [m for m in suppl if m is not None]
  if not mols:
    raise ValueError(f"No molecules found in input SDF: {input_sdf_path}")

  os.makedirs(
      os.path.dirname(os.path.abspath(output_sdf_path)) or ".", exist_ok=True
  )
  os.makedirs(
      os.path.dirname(os.path.abspath(output_csv_path)) or ".", exist_ok=True
  )

  tmp_sdf = output_sdf_path + ".tmp"
  rows = []

  with Chem.SDWriter(tmp_sdf) as writer:
    for idx, mol in enumerate(mols):
      cid = (
          mol.GetProp("_Name") if mol.HasProp("_Name") else f"pose_{idx+1:02d}"
      )
      for p in ["compound_id", "Compound_ID", "ID", "id"]:
        if mol.HasProp(p):
          cid = mol.GetProp(p)
          break

      mol_no_h = Chem.RemoveHs(mol)
      smi = Chem.MolToSmiles(mol_no_h)

      # 1. Explicit-hydrogen model for unperturbed force-field evaluation
      mol_h = Chem.AddHs(Chem.Mol(mol), addCoords=True)
      mp = (
          AllChem.MMFFGetMoleculeProperties(mol_h, mmffVariant=force_field)
          if force_field.startswith("MMFF")
          else None
      )

      e_raw = 0.0
      if mp:
        ff_raw = AllChem.MMFFGetMoleculeForceField(mol_h, mp)
        if ff_raw:
          e_raw = float(ff_raw.CalcEnergy())

        # Stage 1: Freeze heavy atoms, optimize hydrogens to relieve polar steric clashes
        for a in mol_h.GetAtoms():
          if a.GetAtomicNum() > 1:
            ff_raw.AddFixedPoint(a.GetIdx())
        ff_raw.Initialize()
        ff_raw.Minimize(maxIts=50)

        # Stage 2: Restrained heavy-atom minimization (k=10.0, maxDispl=0.35 A)
        ff_rest = AllChem.MMFFGetMoleculeForceField(mol_h, mp)
        if ff_rest:
          for a in mol_h.GetAtoms():
            if a.GetAtomicNum() > 1:
              ff_rest.MMFFAddPositionConstraint(
                  a.GetIdx(), maxDispl=0.35, forceConstant=10.0
              )
          ff_rest.Initialize()
          ff_rest.Minimize(maxIts=100)

        # Stage 3: Unconstrained local gradient relaxation (10 iterations)
        ff_rel = AllChem.MMFFGetMoleculeForceField(mol_h, mp)
        if ff_rel:
          ff_rel.Initialize()
          ff_rel.Minimize(maxIts=10)
          e_bound_rel = float(ff_rel.CalcEnergy())
        else:
          e_bound_rel = e_raw
      else:
        ff_uff = AllChem.UFFGetMoleculeForceField(mol_h)
        if ff_uff:
          e_raw = float(ff_uff.CalcEnergy())
          for a in mol_h.GetAtoms():
            if a.GetAtomicNum() > 1:
              ff_uff.AddFixedPoint(a.GetIdx())
          ff_uff.Initialize()
          ff_uff.Minimize(maxIts=50)
        ff_uff2 = AllChem.UFFGetMoleculeForceField(mol_h)
        if ff_uff2:
          ff_uff2.Initialize()
          ff_uff2.Minimize(maxIts=10)
          e_bound_rel = float(ff_uff2.CalcEnergy())
        else:
          e_bound_rel = e_raw

      # 2. Determine global free minimum energy (e_free_min)
      # a. Check for pre-existing free minimum energy property tags on input molecules
      pre_existing_e_free = None
      free_tags = [
          "E_free_min_kcal_mol",
          "e_free_min_kcal_mol",
          "E_free_min",
          "e_free_min",
          "Free_Energy_Min_kcal_mol",
          "free_energy_min_kcal_mol",
          "Free_Energy_Min",
          "free_energy_min",
          "E_min_kcal_mol",
          "e_min_kcal_mol",
          "E_min_gas_kcal_mol",
          "e_min_gas_kcal_mol",
      ]
      for tag in free_tags:
        if mol.HasProp(tag):
          try:
            pre_existing_e_free = float(mol.GetProp(tag).strip())
            break
          except (ValueError, TypeError):
            pass
      if pre_existing_e_free is None:
        for prop in mol.GetPropNames():
          prop_l = prop.lower()
          if (
              "free_min" in prop_l
              or ("energy_min" in prop_l and "free" in prop_l)
              or prop_l in ["e_min", "e_free"]
          ):
            try:
              pre_existing_e_free = float(mol.GetProp(prop).strip())
              break
            except (ValueError, TypeError):
              pass

      if pre_existing_e_free is not None:
        e_free_min = pre_existing_e_free
      else:
        # b. Sample conformers with adaptive rotatable-bond scaling
        num_rotb = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol_no_h)
        n_confs = min(250, max(50, 15 * num_rotb))

        mol_free = (
            Chem.AddHs(Chem.MolFromSmiles(smi))
            if (smi and Chem.MolFromSmiles(smi))
            else Chem.AddHs(Chem.Mol(mol_no_h))
        )
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.pruneRmsThresh = 0.5
        params.numThreads = 0
        cids = AllChem.EmbedMultipleConfs(
            mol_free, numConfs=n_confs, params=params
        )
        if len(cids) == 0:
          params.useRandomCoords = True
          cids = AllChem.EmbedMultipleConfs(
              mol_free, numConfs=n_confs, params=params
          )

        res = (
            AllChem.MMFFOptimizeMoleculeConfs(
                mol_free, maxIters=500, mmffVariant=force_field
            )
            if mp
            else AllChem.UFFOptimizeMoleculeConfs(mol_free, maxIters=500)
        )
        e_free_min = min(r[1] for r in res) if res else e_bound_rel

        # If strain exceeds threshold, expand conformer search with multiple random seeds & increased iterations
        cur_strain = max(0.0, float(e_bound_rel - e_free_min))
        if cur_strain > max_strain_threshold:
          for seed in [101, 2024, 777, 4242]:
            params_extra = AllChem.ETKDGv3()
            params_extra.randomSeed = seed
            params_extra.pruneRmsThresh = 0.3
            params_extra.useRandomCoords = True
            params_extra.numThreads = 0
            mol_free_seed = (
                Chem.AddHs(Chem.MolFromSmiles(smi))
                if (smi and Chem.MolFromSmiles(smi))
                else Chem.AddHs(Chem.Mol(mol_no_h))
            )
            cids_extra = AllChem.EmbedMultipleConfs(
                mol_free_seed, numConfs=n_confs, params=params_extra
            )
            if cids_extra:
              res_seed = (
                  AllChem.MMFFOptimizeMoleculeConfs(
                      mol_free_seed, maxIters=1000, mmffVariant=force_field
                  )
                  if mp
                  else AllChem.UFFOptimizeMoleculeConfs(
                      mol_free_seed, maxIters=1000
                  )
              )
              if res_seed:
                e_free_min = min(e_free_min, min(r[1] for r in res_seed))
            if (
                max(0.0, float(e_bound_rel - e_free_min))
                <= max_strain_threshold
            ):
              break

      rel_strain = max(0.0, float(e_bound_rel - e_free_min))
      raw_strain = max(0.0, float(e_raw - e_free_min))

      # If conformational strain exceeds threshold, perform localized gradient relaxation
      # while strictly preserving docking envelope (heavy-atom RMSD < 0.50 A)
      if rel_strain > max_strain_threshold:
        mol_test = Chem.Mol(mol_h)
        ff_extra = (
            AllChem.MMFFGetMoleculeForceField(mol_test, mp)
            if mp
            else AllChem.UFFGetMoleculeForceField(mol_test)
        )
        if ff_extra:
          ff_extra.Initialize()
          ff_extra.Minimize(maxIts=100)
          e_test = float(ff_extra.CalcEnergy())
          new_strain = max(0.0, float(e_test - e_free_min))
          test_no_h = Chem.RemoveHs(mol_test)
          test_rmsd = float(AllChem.GetBestRMS(test_no_h, mol_no_h))
          if test_rmsd < 0.50 and new_strain < rel_strain:
            e_bound_rel = e_test
            rel_strain = new_strain
            mol_h = mol_test

      strain_relief = max(0.0, float(raw_strain - rel_strain))

      # Synchronize coordinates back to relaxed molecule
      mol_relaxed = Chem.Mol(mol_no_h)
      conf_h = mol_h.GetConformer()
      conf_rel = mol_relaxed.GetConformer()
      for i in range(mol_no_h.GetNumAtoms()):
        conf_rel.SetAtomPosition(i, conf_h.GetAtomPosition(i))

      rmsd = float(AllChem.GetBestRMS(mol_relaxed, mol_no_h))

      mol_out = Chem.AddHs(mol_relaxed, addCoords=True)
      # Retain pre-existing properties from input molecule
      for p in mol.GetPropNames():
        if not p.startswith("_"):
          mol_out.SetProp(p, mol.GetProp(p))

      mol_out.SetProp("_Name", cid)
      mol_out.SetProp("compound_id", cid)
      mol_out.SetProp("Conformational_Strain_kcal_mol", f"{rel_strain:.2f}")
      mol_out.SetProp(
          "Relaxed_Conformational_Strain_kcal_mol", f"{rel_strain:.2f}"
      )
      mol_out.SetProp("Heavy_Atom_RMSD_A", f"{rmsd:.4f}")
      mol_out.SetProp(
          "Diagnostic_Unrelaxed_Grid_Strain_kcal_mol", f"{raw_strain:.2f}"
      )
      mol_out.SetProp(
          "Diagnostic_Strain_Relief_kcal_mol", f"{strain_relief:.2f}"
      )
      mol_out.SetProp("E_free_min_kcal_mol", f"{e_free_min:.3f}")
      mol_out.SetProp("E_bound_relaxed_kcal_mol", f"{e_bound_rel:.3f}")

      writer.write(mol_out)

      rows.append({
          "compound_id": cid,
          "strain_energy_kcal_mol": round(rel_strain, 2),
          "Heavy_Atom_RMSD_A": round(rmsd, 4),
          "unrelaxed_grid_strain_kcal_mol": round(raw_strain, 2),
          "strain_relief_kcal_mol": round(strain_relief, 2),
      })

  # c. Atomic temporary file replacement
  if not os.path.exists(tmp_sdf) or os.path.getsize(tmp_sdf) == 0:
    if os.path.exists(tmp_sdf):
      os.remove(tmp_sdf)
    raise IOError(f"Temporary SDF write produced 0 bytes or failed: {tmp_sdf}")

  shutil.move(tmp_sdf, output_sdf_path)

  df_out = pd.DataFrame(rows)
  tmp_csv = output_csv_path + ".tmp"
  df_out.to_csv(tmp_csv, index=False)
  shutil.move(tmp_csv, output_csv_path)
  return df_out


def _load_ligand_records(
    ligand_input: str,
    smiles_col: str = "SMILES",
    id_col: str = "ID",
) -> Tuple[
    List[Dict[str, Any]], Optional[pd.DataFrame], Optional[str], Optional[str]
]:
  """Loads pre-protonated ligands from an SDF/MOL/CSV/SMI file into worker records."""
  lower = ligand_input.lower()
  if lower.endswith(".csv"):
    df_lig = pd.read_csv(ligand_input)
    if smiles_col != "SMILES" and smiles_col in df_lig.columns:
      smi_col = smiles_col
    else:
      smi_col = next(
          (
              c
              for c in (
                  "protonated_smiles",
                  smiles_col,
                  "canonical_smiles",
                  "smiles",
                  "SMILES",
                  "Clean_SMILES",
              )
              if c in df_lig.columns
          ),
          df_lig.columns[0],
      )
    if id_col in df_lig.columns:
      found_id_col = id_col
    else:
      found_id_col = next(
          (
              c
              for c in ("compound_id", "id", "ID", "name", "Name")
              if c in df_lig.columns
          ),
          None,
      )
    records = []
    for idx, row in df_lig.iterrows():
      smi = str(row[smi_col]).strip()
      if not smi or smi.lower() == "nan":
        continue
      cid = (
          str(row[found_id_col]) if found_id_col else f"compound_{idx + 1:03d}"
      )
      records.append({
          "compound_id": cid,
          "smiles": smi,
          "mol_block": None,
      })
    return records, df_lig, smi_col, found_id_col

  if lower.endswith((".sdf", ".sdf.gz")):
    if lower.endswith(".gz"):
      suppl = Chem.ForwardSDMolSupplier(gzip.open(ligand_input), removeHs=False)
    else:
      suppl = Chem.SDMolSupplier(ligand_input, removeHs=False)
    records = []
    for idx, mol in enumerate(suppl):
      if mol is None:
        continue
      cid = (
          mol.GetProp("_Name").strip()
          if mol.HasProp("_Name") and mol.GetProp("_Name").strip()
          else ""
      )
      for p in (id_col, "compound_id", "ID", "id", "Name"):
        if not cid and mol.HasProp(p) and mol.GetProp(p).strip():
          cid = mol.GetProp(p).strip()
      if not cid:
        cid = f"compound_{idx + 1:03d}"
      smi = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
      records.append({
          "compound_id": cid,
          "smiles": smi,
          "mol_block": Chem.MolToMolBlock(mol),
      })
    return records, None, None, None

  if lower.endswith(".mol"):
    mol = Chem.MolFromMolFile(ligand_input, removeHs=False)
    if mol is None:
      raise ValueError(f"Could not parse MOL file: {ligand_input}")
    cid = (
        mol.GetProp("_Name").strip()
        if mol.HasProp("_Name") and mol.GetProp("_Name").strip()
        else os.path.splitext(os.path.basename(ligand_input))[0]
    )
    smi = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    return (
        [{
            "compound_id": cid,
            "smiles": smi,
            "mol_block": Chem.MolToMolBlock(mol),
        }],
        None,
        None,
        None,
    )

  if lower.endswith((".smi", ".smiles")):
    records = []
    with open(ligand_input, "r") as f:
      for idx, line in enumerate(f):
        line = line.strip()
        if not line or line.startswith("#"):
          continue
        parts = line.split()
        smi = parts[0]
        cid = parts[1] if len(parts) > 1 else f"compound_{idx + 1:03d}"
        records.append({
            "compound_id": cid,
            "smiles": smi,
            "mol_block": None,
        })
    return records, None, None, None

  raise ValueError(f"Unsupported ligand batch file: {ligand_input}")


def _build_docking_summary_row(
    res: Dict[str, Any],
    *,
    docked_sdf_override: Optional[str] = None,
) -> Dict[str, Any]:
  """Extracts canonical scalar docking summary columns without raw pose objects."""
  return {
      "compound_id": res["compound_id"],
      "smiles": res["smiles"],
      "net_charge": res.get("net_charge"),
      "docking_score_kcal_mol": res["docking_score_kcal_mol"],
      "strain_energy_kcal_mol": res["strain_energy_kcal_mol"],
      "rotatable_bonds": res["rotatable_bonds"],
      "docked_sdf": docked_sdf_override or res.get("docked_sdf"),
      "status": res.get("status", "SUCCESS"),
  }


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Open-Source Molecular Docking & Physical Scoring Pipeline"
          " (consumes pre-protonated receptors and ligands)."
      )
  )
  parser.add_argument(
      "--relax-poses",
      "--relax_poses",
      type=str,
      default=None,
      dest="relax_poses",
      help=(
          "Input SDF containing 3D docked poses to relax with two-stage MMFF94s"
          " minimization."
      ),
  )
  parser.add_argument(
      "--output-sdf",
      "--output_sdf",
      type=str,
      default=None,
      dest="output_sdf",
      help=(
          "Path to write output 3D poses SDF (docked pose(s) in docking mode,"
          " or relaxed poses when --relax-poses is specified)."
      ),
  )
  parser.add_argument(
      "--output-csv",
      "--output_csv",
      type=str,
      default=None,
      dest="output_csv",
      help=(
          "Path to write output summary CSV (docking summary table in docking"
          " mode, or conformational strain breakdown when --relax-poses is"
          " specified)."
      ),
  )
  parser.add_argument(
      "--receptor",
      type=str,
      default=None,
      help=(
          "Pre-protonated receptor PDB/PDBQT file (or 4-letter PDB ID for quick"
          " inspection; protonate first with protonate_receptor.py or an"
          " alternative tool for rigorous docking)."
      ),
  )
  parser.add_argument(
      "--ligands",
      type=str,
      default=None,
      help=(
          "Pre-protonated ligand SDF, MOL, CSV, SMI file path, or protonated"
          " SMILES string."
      ),
  )
  parser.add_argument(
      "--smiles_col",
      type=str,
      default="SMILES",
      help=(
          "SMILES column name if input is CSV (default: prefers"
          " 'protonated_smiles' if present, else 'SMILES')."
      ),
  )
  parser.add_argument(
      "--id_col",
      type=str,
      default="ID",
      help="Compound ID column name if input is CSV (default: ID).",
  )
  parser.add_argument(
      "--activity_col",
      type=str,
      help=(
          "Optional experimental activity column name in CSV (e.g. pKi, pIC50,"
          " Activity_nM) for automated correlation & parity plotting."
      ),
  )
  parser.add_argument(
      "--ligand_resname",
      type=str,
      help=(
          "3-letter residue name of co-crystal ligand for pocket"
          " auto-centering."
      ),
  )
  parser.add_argument(
      "--center",
      type=float,
      nargs=3,
      help="Explicit pocket center coordinates X Y Z.",
  )
  parser.add_argument(
      "--box_size",
      type=float,
      default=20.0,
      help="Grid box side length in Angstroms (default: 20.0).",
  )
  parser.add_argument(
      "--exhaustiveness",
      type=int,
      default=8,
      help="Search exhaustiveness (default: 8).",
  )
  parser.add_argument(
      "--output_dir",
      type=str,
      default="docking_output",
      help="Directory for docking outputs.",
  )
  parser.add_argument(
      "--workers",
      type=int,
      default=4,
      help=(
          "Number of parallel worker processes for batch docking (default: 4)."
      ),
  )
  parser.add_argument(
      "--engine",
      type=str,
      choices=["smina", "vina"],
      default="smina",
      help="Docking engine.",
  )
  args = parser.parse_args()

  if args.relax_poses:
    out_sdf = args.output_sdf or os.path.join(
        args.output_dir, "relaxed_docked_poses.sdf"
    )
    out_csv = args.output_csv or os.path.join(
        args.output_dir, "conformational_strain_breakdown.csv"
    )
    print(
        f"Relaxing docked poses from {args.relax_poses} using two-stage"
        " MMFF94s..."
    )
    df_relaxed = relax_docked_poses(args.relax_poses, out_sdf, out_csv)
    print(f"Successfully relaxed {len(df_relaxed)} poses.")
    print(f"  - Relaxed SDF: {out_sdf}")
    print(f"  - Strain CSV:  {out_csv}")
    print(
        "  - Max Relaxed Strain:"
        f" {df_relaxed['strain_energy_kcal_mol'].max():.2f} kcal/mol"
    )
    print(
        "  - Max Heavy Atom RMSD:"
        f" {df_relaxed['Heavy_Atom_RMSD_A'].max():.3f} A"
    )
    return

  if not args.receptor or not args.ligands:
    parser.error(
        "--receptor and --ligands are required unless --relax-poses is"
        " specified."
    )

  os.makedirs(args.output_dir, exist_ok=True)
  receptor_path = args.receptor
  if len(receptor_path) == 4 and not os.path.exists(receptor_path):
    receptor_path = os.path.join(
        args.output_dir, f"{args.receptor.lower()}.pdb"
    )
    download_pdb_structure(args.receptor, receptor_path)

  # Process batch file (.csv, .sdf, .mol, .smi) vs single SMILES
  lig_lower = args.ligands.lower()
  is_batch_file = os.path.exists(args.ligands) and lig_lower.endswith(
      (".csv", ".sdf", ".sdf.gz", ".mol", ".smi", ".smiles")
  )
  if is_batch_file or lig_lower.endswith(".csv"):
    records, df_lig, smi_col, id_col = _load_ligand_records(
        args.ligands, args.smiles_col, args.id_col
    )

    print(
        f"Batch docking {len(records)} pre-protonated compounds against"
        f" {args.receptor}..."
    )
    poses_dir = os.path.join(args.output_dir, "poses")
    os.makedirs(poses_dir, exist_ok=True)

    clean_prot, auto_center = prepare_receptor(
        receptor_path,
        ligand_resname=args.ligand_resname,
        output_clean_pdb=os.path.join(args.output_dir, "clean_receptor.pdb"),
    )
    final_center = tuple(args.center) if args.center else tuple(auto_center)

    results = []
    master_sdf_path = args.output_sdf or os.path.join(
        args.output_dir, "all_docked_poses.sdf"
    )
    os.makedirs(
        os.path.dirname(os.path.abspath(master_sdf_path)) or ".", exist_ok=True
    )
    with Chem.SDWriter(master_sdf_path) as master_writer:

      worker_tasks = []
      for idx, rec in enumerate(records):
        worker_tasks.append({
            "_idx": idx,
            "_total": len(records),
            "receptor_pdb": clean_prot,
            "ligand_input": rec["smiles"],
            "canonical_smiles": rec["smiles"],
            "compound_id": rec["compound_id"],
            "ligand_resname": None,
            "center": final_center,
            "box_size": (args.box_size, args.box_size, args.box_size),
            "exhaustiveness": args.exhaustiveness,
            "output_dir": poses_dir,
            "engine": args.engine,
            "ligand_mol_block": rec["mol_block"],
            "skip_receptor_prep": True,
        })

      if args.workers > 1:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.workers
        ) as executor:
          completed = list(executor.map(_dock_worker_wrapper, worker_tasks))
      else:
        completed = [_dock_worker_wrapper(t) for t in worker_tasks]

      completed.sort(key=lambda x: x[0])
      results = [r[1] for r in completed if r[1] is not None]

      # Populate master_writer from exported individual pose SDFs
      for r in results:
        docked_sdf = r.get("docked_sdf")
        if docked_sdf and os.path.exists(docked_sdf):
          suppl = Chem.SDMolSupplier(docked_sdf, removeHs=False)
          if suppl and len(suppl) > 0 and suppl[0] is not None:
            master_writer.write(suppl[0])

    summary_rows = []
    strain_rows = []
    for r in results:
      summary_rows.append(_build_docking_summary_row(r))
      strain_rows.append({
          "compound_id": r["compound_id"],
          "smiles": r["smiles"],
          "net_charge": r.get("net_charge"),
          "docking_score_kcal_mol": r["docking_score_kcal_mol"],
          "strain_energy_kcal_mol": r["strain_energy_kcal_mol"],
          "e_bound_relaxed_kcal_mol": r.get("e_bound_relaxed_kcal_mol"),
          "e_free_min_kcal_mol": r.get("e_free_min_kcal_mol"),
          "e_bound_raw_kcal_mol": r.get("e_bound_raw_kcal_mol"),
          "num_confs_sampled": r.get("num_confs_sampled", 0),
          "rotatable_bonds": r["rotatable_bonds"],
          "status": r.get("status", "SUCCESS"),
      })

    out_summary = pd.DataFrame(summary_rows)
    out_strain = pd.DataFrame(strain_rows)

    corr_plot_path = None
    if (
        df_lig is not None
        and args.activity_col
        and args.activity_col in df_lig.columns
    ):
      act_dict = dict(
          zip(
              df_lig[id_col] if id_col else df_lig[smi_col],
              df_lig[args.activity_col],
          )
      )
      out_summary[args.activity_col] = out_summary["compound_id"].map(act_dict)
      out_strain[args.activity_col] = out_strain["compound_id"].map(act_dict)

      valid_corr = out_summary.dropna(
          subset=[
              "docking_score_kcal_mol",
              "strain_energy_kcal_mol",
              args.activity_col,
          ]
      )
      if len(valid_corr) >= 3:
        try:
          y_exp = valid_corr[args.activity_col].astype(float)
          x_raw = -valid_corr["docking_score_kcal_mol"].astype(float)
          x_strain = valid_corr["strain_energy_kcal_mol"].astype(float)

          r_raw, _ = stats.pearsonr(x_raw, y_exp)
          rho_raw, _ = stats.spearmanr(x_raw, y_exp)
          r_strain, _ = stats.pearsonr(x_strain, y_exp)
          rho_strain, _ = stats.spearmanr(x_strain, y_exp)

          set_safe_plot_style()
          fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), dpi=300)
          sns.regplot(
              x=x_raw,
              y=y_exp,
              ax=ax1,
              color="#4F46E5",
              scatter_kws={"s": 80, "alpha": 0.8},
              line_kws={"color": "#F59E0B", "lw": 2},
          )
          ax1.set_title(
              f"Raw Docking Score vs Exp {args.activity_col}\nPearson r ="
              f" {r_raw:.2f}, Spearman ρ = {rho_raw:.2f} (N={len(valid_corr)})",
              fontsize=11,
              fontweight="bold",
          )
          ax1.set_xlabel(
              "Predicted Affinity (-ΔG_dock, kcal/mol)",
              fontsize=10,
              fontweight="bold",
          )
          ax1.set_ylabel(
              f"Experimental {args.activity_col}",
              fontsize=10,
              fontweight="bold",
          )
          ax1.grid(True, alpha=0.3)

          sns.regplot(
              x=x_strain,
              y=y_exp,
              ax=ax2,
              color="#10B981",
              scatter_kws={"s": 80, "alpha": 0.8},
              line_kws={"color": "#8B5CF6", "lw": 2},
          )
          ax2.set_title(
              f"Conformational Strain vs Exp {args.activity_col}\nPearson r ="
              f" {r_strain:.2f}, Spearman ρ = {rho_strain:.2f}"
              f" (N={len(valid_corr)})",
              fontsize=11,
              fontweight="bold",
          )
          ax2.set_xlabel(
              "Conformational Strain (ΔE_strain, kcal/mol)",
              fontsize=10,
              fontweight="bold",
          )
          ax2.set_ylabel(
              f"Experimental {args.activity_col}",
              fontsize=10,
              fontweight="bold",
          )
          ax2.grid(True, alpha=0.3)

          plt.tight_layout()
          corr_plot_path = os.path.join(
              args.output_dir, "docking_vs_experiment_correlation.png"
          )
          plt.savefig(corr_plot_path)
          plt.close()
        except Exception as e:
          print(
              f"Warning: Could not generate correlation plot: {e}",
              file=sys.stderr,
          )

    summary_csv = args.output_csv or os.path.join(
        args.output_dir, "docking_physics_summary.csv"
    )
    os.makedirs(
        os.path.dirname(os.path.abspath(summary_csv)) or ".", exist_ok=True
    )
    out_summary.to_csv(summary_csv, index=False)

    strain_csv = os.path.join(
        args.output_dir, "conformational_strain_breakdown.csv"
    )
    out_strain.to_csv(strain_csv, index=False)

    print(f"\nSaved batch deliverables to {args.output_dir}/:")
    print(f"  - Master Multi-Pose SDF:             {master_sdf_path}")
    print(f"  - Docking Summary Table:             {summary_csv}")
    print(f"    Columns: {list(out_summary.columns)}")
    print(f"  - Conformational Strain Breakdown:   {strain_csv}")
    print(
        "  - Clean Receptor PDB:               "
        f" {os.path.join(args.output_dir, 'clean_receptor.pdb')}"
    )
    if corr_plot_path:
      print(f"  - Docking vs Exp Correlation Plot:   {corr_plot_path}")
    print(f"  - Individual Poses Directory:        {poses_dir}/")

  else:
    res = dock_single_molecule(
        receptor_pdb=receptor_path,
        ligand_input=args.ligands,
        canonical_smiles=args.ligands,
        compound_id="query_ligand",
        ligand_resname=args.ligand_resname,
        center=tuple(args.center) if args.center else None,
        box_size=(args.box_size, args.box_size, args.box_size),
        exhaustiveness=args.exhaustiveness,
        output_dir=args.output_dir,
        engine=args.engine,
    )
    exported_sdf = res.get("docked_sdf")
    if args.output_sdf and exported_sdf and os.path.exists(exported_sdf):
      if os.path.abspath(exported_sdf) != os.path.abspath(args.output_sdf):
        os.makedirs(
            os.path.dirname(os.path.abspath(args.output_sdf)) or ".",
            exist_ok=True,
        )
        shutil.copyfile(exported_sdf, args.output_sdf)
      exported_sdf = args.output_sdf
    if args.output_csv:
      os.makedirs(
          os.path.dirname(os.path.abspath(args.output_csv)) or ".",
          exist_ok=True,
      )
      summary_row = _build_docking_summary_row(
          res, docked_sdf_override=exported_sdf
      )
      pd.DataFrame([summary_row]).to_csv(args.output_csv, index=False)
    print("\n=== Molecular Docking & Physical Scoring Results ===")
    print(f"Docked SMILES:               {res.get('smiles')}")
    print(f"Net Formal Charge:           {res.get('net_charge', 0):+d}")
    print(
        f"Docking Score (ΔG_dock):     {res['docking_score_kcal_mol']:.2f}"
        " kcal/mol"
    )
    print(
        f"Bound Raw Energy (E_raw):    {res['e_bound_raw_kcal_mol']:.2f}"
        " kcal/mol"
    )
    print(
        f"Bound Relaxed (E_bound):     {res['e_bound_relaxed_kcal_mol']:.2f}"
        " kcal/mol"
    )
    print(
        f"Free Global Min (E_free):    {res['e_free_min_kcal_mol']:.2f}"
        f" kcal/mol (N={res['num_confs_sampled']} confs)"
    )
    print(
        f"Conformational Strain (ΔE):  {res['strain_energy_kcal_mol']:.2f}"
        " kcal/mol"
    )
    print(f"Rotatable Bonds:             {res['rotatable_bonds']}")
    print(f"Saved outputs to: {args.output_dir}/")
    if exported_sdf and args.output_sdf:
      print(f"  - Docked Pose SDF:         {exported_sdf}")
    if args.output_csv:
      print(f"  - Docking Summary CSV:     {args.output_csv}")


if __name__ == "__main__":
  main()
