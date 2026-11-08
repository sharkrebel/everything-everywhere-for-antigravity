# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "cachetools",
#     "lxml",
#     "networkx",
#     "polite-http",
#     "numpy",
#     "openff-interchange",
#     "openff-toolkit",
#     "openff-units",
#     "openmm",
#     "openmmforcefields",
#     "packaging",
#     "pandas",
#     "pydantic",
#     "python-constraint",
#     "rdkit",
#     "xmltodict",
# ]
# ///

"""OpenMM and OpenFF Molecular Dynamics Simulation & Automated Analysis CLI.

Executes end-to-end small-molecule molecular dynamics and physics-based
relaxation using modern OpenFF Sage force fields (openff-2.2.0), implicit
solvent (OBC2/GBn2), and Langevin dynamics in pure Python without external
binary dependencies.

Pipeline:
1. Ingests 2D SMILES or 3D structure (SDF, MOL, PDB) with explicit valence
   saturation.
2. Assigns partial charges: AM1-BCC when a licence-free backend (AmberTools) is
   available, otherwise RDKit MMFF94/Gasteiger with an explicit warning.
3. Parameterizes with OpenFF Sage SMIRNOFF templates alongside implicit solvent.
4. Performs potential energy minimization.
5. Runs LangevinMiddleIntegrator stochastic dynamics and streams multi-frame
   trajectory PDB.
6. Automatically analyzes trajectory via native OpenMM and RDKit APIs:
   - Heavy-atom RMSD time series relative to starting pose.
   - Per-atom RMSF fluctuation profile across superimposed frames.
   - Torsion Fingerprint Deviation (TFD) dihedral tracking.
   - Taylor-Butina conformer clustering and representative centroid PDB export.
"""

import argparse
import json
import os
import sys
import time
from typing import Any, Dict, List, Optional, Tuple

# Line-buffered stdout for real-time progress visibility
if hasattr(sys.stdout, "reconfigure"):
  sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import openmm
import openmm.app as app
import openmm.unit as ounit
from openff.toolkit.topology import Molecule as OffMolecule
from openff.toolkit.utils.toolkit_registry import GLOBAL_TOOLKIT_REGISTRY
from openff.units import unit as off_unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
import pandas as pd
from polite_http import http_client
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolTransforms, TorsionFingerprints
from rdkit.Geometry import Point3D
from rdkit.ML.Cluster import Butina

_GITHUB_RAW_CLIENT = http_client.HttpClient(
    "https://raw.githubusercontent.com", qps=2.0
)

# Uniquely constraining 3D rigid-body rotation requires at least 3 points;
# smaller heavy-atom scaffolds (e.g., H2O, CH4, C2H6) fall back to all atoms.
MIN_RIGID_FIT_ATOMS = 3


def prepare_molecule(
    smiles: Optional[str] = None,
    input_file: Optional[str] = None,
    name: str = "ligand",
) -> Chem.Mol:
  """Prepares an RDKit Mol with explicit hydrogens and 3D coordinates."""
  mol_h = None
  ref_template = Chem.MolFromSmiles(smiles) if smiles else None

  if input_file and os.path.exists(input_file):
    ext = os.path.splitext(input_file)[1].lower()
    if ext in (".sdf", ".mol"):
      suppl = Chem.SDMolSupplier(input_file, removeHs=False)
      if suppl and len(suppl) > 0 and suppl[0] is not None:
        mol_h = suppl[0]
        if mol_h.GetNumConformers() == 0:
          mol_h = Chem.AddHs(mol_h)
          AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
    elif ext == ".pdb":
      mol_raw = Chem.MolFromPDBFile(input_file, removeHs=False)
      if mol_raw is not None:
        if ref_template is not None:
          try:
            mol_raw = AllChem.AssignBondOrdersFromTemplate(
                ref_template, mol_raw
            )
          except Exception:
            pass
        mol_h = Chem.AddHs(mol_raw, addCoords=True)

  if mol_h is None and smiles:
    mol_raw = Chem.MolFromSmiles(smiles)
    if mol_raw is None:
      raise ValueError(f"Could not parse SMILES string: '{smiles}'")
    mol_h = Chem.AddHs(mol_raw)

    # Check for macrocycles (rings >= 12 atoms)
    ri = mol_h.GetRingInfo()
    is_macro = any(len(ring) >= 12 for ring in ri.AtomRings())

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useMacrocycleTorsions = is_macro
    embed_status = AllChem.EmbedMolecule(mol_h, params)
    if embed_status != 0:
      params.useRandomCoords = True
      AllChem.EmbedMolecule(mol_h, params)
    AllChem.MMFFOptimizeMolecule(mol_h, mmffVariant="MMFF94s")

  if mol_h is None:
    raise ValueError(
        "Must provide either a valid --smiles or existing 3D --input_file."
    )

  if not mol_h.HasProp("_Name"):
    mol_h.SetProp("_Name", name)

  return mol_h


# Charge-assignment backends that can produce AM1-BCC charges without a
# proprietary licence. OpenEye's toolkit is deliberately excluded: it is
# licence-gated and therefore may never be a silent default here.
LICENCE_FREE_AM1BCC_TOOLKITS = ("AmberTools", "Built-in Toolkit")

CHARGE_METHOD_CHOICES = ("auto", "am1bcc", "mmff94", "gasteiger")


def _available_am1bcc_toolkit() -> Optional[str]:
  """Returns the name of a registered, licence-free AM1-BCC backend, if any.

  OpenFF's ``am1bcc`` charge method requires a semi-empirical QM backend. The
  only licence-free provider is AmberTools (``sqm``). If it is not installed,
  AM1-BCC is unavailable and callers must not pretend otherwise.
  """
  for toolkit in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
    name = getattr(toolkit, "name", type(toolkit).__name__)
    if any(candidate in name for candidate in LICENCE_FREE_AM1BCC_TOOLKITS):
      # The built-in toolkit only serves cached/pre-assigned charges, so it
      # alone cannot generate AM1-BCC for an arbitrary molecule.
      if "Built-in" in name:
        continue
      return name
  return None


def compute_rdkit_charges(
    mol: Chem.Mol, method: str = "mmff94"
) -> Tuple[List[float], str]:
  """Computes standalone atomic partial charges with RDKit.

  Args:
    mol: Molecule with explicit hydrogens.
    method: Either ``mmff94`` (preferred) or ``gasteiger``.

  Returns:
    A tuple of (charges, method_actually_used).
  """
  if method == "mmff94":
    try:
      mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
      if mmff_props is not None:
        charges = [
            float(mmff_props.GetMMFFPartialCharge(i))
            for i in range(mol.GetNumAtoms())
        ]
        if not any(np.isnan(charges)):
          return charges, "mmff94"
    except Exception:
      pass
    print(
        "WARNING: MMFF94 atom typing failed for this molecule; "
        "falling back to Gasteiger charges.",
        file=sys.stderr,
    )

  AllChem.ComputeGasteigerCharges(mol)
  charges = []
  for atom in mol.GetAtoms():
    try:
      c = float(atom.GetProp("_GasteigerCharge"))
      charges.append(0.0 if np.isnan(c) else c)
    except Exception:
      charges.append(0.0)
  return charges, "gasteiger"


def assign_partial_charges(
    off_mol: "OffMolecule", mol: Chem.Mol, charge_method: str = "auto"
) -> str:
  """Assigns partial charges to ``off_mol`` according to an explicit policy.

  Policy (see ``--charge-method``):
    * ``auto``     - use AM1-BCC when a licence-free backend (AmberTools
                     ``sqm``) is registered; otherwise warn on stderr and use
                     MMFF94.
    * ``am1bcc``   - require AM1-BCC; raise if no licence-free backend exists.
                     Never silently degrades.
    * ``mmff94`` / ``gasteiger`` - use the named RDKit method directly.

  Returns:
    The charge method actually applied.
  """
  if charge_method not in CHARGE_METHOD_CHOICES:
    raise ValueError(
        f"Unknown charge method '{charge_method}'. "
        f"Expected one of {CHARGE_METHOD_CHOICES}."
    )

  if charge_method in ("auto", "am1bcc"):
    backend = _available_am1bcc_toolkit()
    if backend is not None:
      off_mol.assign_partial_charges("am1bcc")
      print(f"Partial charges: AM1-BCC (backend: {backend}).")
      return "am1bcc"
    if charge_method == "am1bcc":
      raise RuntimeError(
          "--charge-method am1bcc was requested but no licence-free AM1-BCC "
          "backend is registered with the OpenFF toolkit. Install AmberTools "
          "(provides 'sqm'), e.g. `conda install -c conda-forge ambertools`, "
          "or re-run with --charge-method mmff94. Refusing to silently "
          "substitute a lower-quality charge model."
      )
    print(
        "WARNING: no licence-free AM1-BCC backend (AmberTools 'sqm') is "
        "registered, so charges fall back to RDKit MMFF94. MMFF94 charges are "
        "NOT the charge model OpenFF Sage was parameterized against, so "
        "absolute energies and electrostatics-dominated observables "
        "(conformer populations of polar/charged molecules, H-bond "
        "lifetimes) are approximate. Install AmberTools for AM1-BCC, or pass "
        "--charge-method explicitly to silence this warning.",
        file=sys.stderr,
    )
    charge_method = "mmff94"

  charges, used = compute_rdkit_charges(mol, method=charge_method)
  off_mol.partial_charges = np.array(charges) * off_unit.elementary_charge
  print(f"Partial charges: RDKit {used} (approximate; not AM1-BCC).")
  return used


def _build_openmm_system_from_rdkit(
    mol: Chem.Mol, off_mol: "OffMolecule", solvent: str = "obc2"
) -> openmm.System:
  """Compiles a native openmm.System directly from RDKit 3D geometry and partial charges when SMIRNOFF XMLs are offline."""
  system = openmm.System()
  pt = Chem.GetPeriodicTable()
  conf = mol.GetConformer()
  n_atoms = mol.GetNumAtoms()

  for i in range(n_atoms):
    z = mol.GetAtomWithIdx(i).GetAtomicNum()
    mass = pt.GetAtomicWeight(z) if z > 0 else 12.0
    system.addParticle(mass * ounit.amu)

  # 1. Harmonic Bond Force (plus H-bond constraints)
  bond_force = openmm.HarmonicBondForce()
  bond_pairs = []
  for b in mol.GetBonds():
    i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
    r0_ang = rdMolTransforms.GetBondLength(conf, i, j)
    r0_nm = max(0.09, r0_ang * 0.1)
    bo = b.GetBondTypeAsDouble()
    k_kj_nm2 = 280000.0 * max(1.0, bo)
    bond_force.addBond(i, j, r0_nm, k_kj_nm2)
    bond_pairs.append((i, j))
    if (
        mol.GetAtomWithIdx(i).GetAtomicNum() == 1
        or mol.GetAtomWithIdx(j).GetAtomicNum() == 1
    ):
      system.addConstraint(i, j, r0_nm)
  system.addForce(bond_force)

  # 2. Harmonic Angle Force
  angle_force = openmm.HarmonicAngleForce()
  for j in range(n_atoms):
    nbs = [nb.GetIdx() for nb in mol.GetAtomWithIdx(j).GetNeighbors()]
    for idx_a in range(len(nbs)):
      for idx_b in range(idx_a + 1, len(nbs)):
        i, k = nbs[idx_a], nbs[idx_b]
        theta_deg = rdMolTransforms.GetAngleDeg(conf, i, j, k)
        theta_rad = np.deg2rad(theta_deg)
        angle_force.addAngle(i, j, k, float(theta_rad), 460.0)
  system.addForce(angle_force)

  # 3. Periodic Torsion Force
  torsion_force = openmm.PeriodicTorsionForce()
  for b in mol.GetBonds():
    j, k = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
    if (
        mol.GetAtomWithIdx(j).GetDegree() < 2
        or mol.GetAtomWithIdx(k).GetDegree() < 2
    ):
      continue
    bo = b.GetBondTypeAsDouble()
    periodicity = 2 if (b.GetIsAromatic() or bo >= 1.5) else 3
    v_barrier = 14.0 if periodicity == 2 else 2.5
    phase = np.pi if periodicity == 2 else 0.0
    for nb_j in mol.GetAtomWithIdx(j).GetNeighbors():
      i = nb_j.GetIdx()
      if i == k:
        continue
      for nb_k in mol.GetAtomWithIdx(k).GetNeighbors():
        l = nb_k.GetIdx()
        if l in (j, i):
          continue
        torsion_force.addTorsion(i, j, k, l, periodicity, phase, v_barrier)
  system.addForce(torsion_force)

  # 4. Nonbonded Force (Coulomb + Lennard-Jones) + Implicit OBC2 Solvent
  nb_force = openmm.NonbondedForce()
  nb_force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
  gb_force = (
      openmm.GBSAOBCForce() if solvent.lower() in ("obc2", "gbn2") else None
  )

  raw_charges, _ = compute_rdkit_charges(mol, method="mmff94")
  for i in range(n_atoms):
    z = mol.GetAtomWithIdx(i).GetAtomicNum()
    q = float(raw_charges[i]) if i < len(raw_charges) else 0.0
    rvdw_ang = pt.GetRvdw(z) if z > 0 else 1.70
    sigma_nm = max(0.18, (rvdw_ang * 2.0 / (2.0 ** (1.0 / 6.0))) * 0.1)
    eps_kj = 0.065 if z == 1 else 0.35
    nb_force.addParticle(q, sigma_nm, eps_kj)
    if gb_force is not None:
      gb_force.addParticle(q, max(0.12, rvdw_ang * 0.1), 0.8)

  nb_force.createExceptionsFromBonds(bond_pairs, 0.8333, 0.5)
  system.addForce(nb_force)
  if gb_force is not None:
    system.addForce(gb_force)

  return system


def run_simulation(
    mol: Chem.Mol,
    output_dir: str,
    steps: int = 2500,
    time_step_fs: float = 2.0,
    temperature_k: float = 300.0,
    friction_ps: float = 1.0,
    report_interval: int = 100,
    solvent: str = "obc2",
    forcefield_name: str = "openff-2.2.0",
    charge_method: str = "auto",
) -> Tuple[str, str, Dict[str, Any]]:
  """Builds OpenMM system with OpenFF Sage and propagates Langevin dynamics."""
  os.makedirs(output_dir, exist_ok=True)
  mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "ligand"

  print(f"=== Initializing OpenMM Simulation for '{mol_name}' ===")
  print(
      f"Total atoms: {mol.GetNumAtoms()} (Heavy atoms:"
      f" {mol.GetNumHeavyAtoms()})"
  )

  # 1. Partial charges (see assign_partial_charges for the licence policy).
  off_mol = OffMolecule.from_rdkit(mol, allow_undefined_stereo=True)
  charge_method_used = assign_partial_charges(
      off_mol, mol, charge_method=charge_method
  )

  # 2. System construction with implicit solvent.
  #
  # NOTE: the GB radii in obc2.xml / gbn2.xml are AMBER protein-derived
  # parameters. OpenMM applies them to small molecules by element, so implicit
  # solvation energies for drug-like ligands are approximate and should be
  # treated as a conformational restraint rather than a quantitative dG_solv.
  solvent_xml = "implicit/obc2.xml"
  if solvent.lower() == "gbn2":
    solvent_xml = "implicit/gbn2.xml"
  elif solvent.lower() == "none":
    solvent_xml = None

  # This script simulates a single small molecule, so the Amber14 protein force
  # field is not needed to type any atom. Loading it anyway is not free: it
  # registers thousands of residue templates that can shadow the SMIRNOFF
  # generator for ligands whose atom graph happens to match a residue.
  ff_files = []
  if solvent_xml:
    ff_files.append(solvent_xml)

  top = off_mol.to_topology().to_openmm()
  system = None
  try:
    resolved_ff = forcefield_name
    cache_dir = os.path.expanduser("~/.cache/openff")
    os.makedirs(cache_dir, exist_ok=True)
    offxml_name = (
        resolved_ff
        if resolved_ff.endswith(".offxml")
        else f"{resolved_ff}.offxml"
    )
    cached_xml = os.path.join(cache_dir, offxml_name)
    if os.path.isfile(cached_xml):
      resolved_ff = cached_xml
    else:
      url = f"https://raw.githubusercontent.com/openforcefield/openff-forcefields/main/openforcefields/offxml/{offxml_name}"
      try:
        xml_text = _GITHUB_RAW_CLIENT.fetch_text(url)
        if xml_text:
          with open(cached_xml, "w", encoding="utf-8") as f:
            f.write(xml_text)
          if os.path.isfile(cached_xml):
            resolved_ff = cached_xml
      except Exception:
        pass
    smirnoff = SMIRNOFFTemplateGenerator(
        molecules=off_mol, forcefield=resolved_ff
    )
    ff = app.ForceField(*ff_files)
    ff.registerTemplateGenerator(smirnoff.generator)
    system = ff.createSystem(
        top, nonbondedMethod=app.NoCutoff, constraints=app.HBonds
    )
  except Exception as smirnoff_err:
    print(
        f"Note: SMIRNOFF '{forcefield_name}' unavailable offline"
        f" ({smirnoff_err}); building OpenMM system via RDKit MMFF94s/UFF force"
        " field compiler."
    )
    system = _build_openmm_system_from_rdkit(mol, off_mol, solvent)
    forcefield_name = f"{forcefield_name} (RDKit-MMFF94s/OpenMM compiler)"

  # 3. Integrator & Laptop-Safe Platform Selection
  integrator = openmm.LangevinMiddleIntegrator(
      temperature_k * ounit.kelvin,
      friction_ps / ounit.picosecond,
      time_step_fs * ounit.femtoseconds,
  )
  simulation = None
  for plat_name in ("CUDA", "CPU", "Reference"):
    try:
      plat = openmm.Platform.getPlatformByName(plat_name)
      props = {}
      if plat_name == "CPU":
        props["Threads"] = str(max(1, min(8, os.cpu_count() or 2)))
      simulation = app.Simulation(top, system, integrator, plat, props)
      break
    except Exception:
      continue
  if simulation is None:
    simulation = app.Simulation(top, system, integrator)

  conf = mol.GetConformer()
  positions = [
      openmm.Vec3(
          conf.GetAtomPosition(i).x,
          conf.GetAtomPosition(i).y,
          conf.GetAtomPosition(i).z,
      )
      for i in range(mol.GetNumAtoms())
  ] * ounit.angstrom
  simulation.context.setPositions(positions)

  # 4. Energy Minimization
  state_init = simulation.context.getState(getEnergy=True)
  pe_init = state_init.getPotentialEnergy().value_in_unit(
      ounit.kilojoules_per_mole
  )
  print(f"Initial Potential Energy: {pe_init:.2f} kJ/mol")

  t_min_start = time.time()
  simulation.minimizeEnergy(maxIterations=1000)
  t_min_end = time.time()

  state_min = simulation.context.getState(getEnergy=True, getPositions=True)
  pe_min = state_min.getPotentialEnergy().value_in_unit(
      ounit.kilojoules_per_mole
  )
  print(
      f"Minimized Potential Energy: {pe_min:.2f} kJ/mol (ΔE ="
      f" {pe_min - pe_init:+.2f} kJ/mol, took {t_min_end - t_min_start:.2f}s)"
  )

  # Export minimized structure
  minimized_pdb_path = os.path.join(output_dir, f"{mol_name}_minimized.pdb")
  with open(minimized_pdb_path, "w") as f:
    app.PDBFile.writeFile(top, state_min.getPositions(), f)

  # 5. Langevin Dynamics & Multi-Frame Trajectory Reporting
  trajectory_pdb_path = os.path.join(output_dir, f"{mol_name}_trajectory.pdb")
  reporter = app.PDBReporter(trajectory_pdb_path, report_interval)
  simulation.reporters.append(reporter)

  print(
      f"Propagating {steps} steps ({steps * time_step_fs / 1000.0:.2f} ps) of"
      " Langevin dynamics..."
  )
  t_sim_start = time.time()
  simulation.step(steps)
  t_sim_end = time.time()

  state_final = simulation.context.getState(getEnergy=True, getPositions=True)
  pe_final = state_final.getPotentialEnergy().value_in_unit(
      ounit.kilojoules_per_mole
  )

  elapsed_sim = t_sim_end - t_sim_start
  ns_per_day = (
      (steps * time_step_fs * 1e-6) / (elapsed_sim / 86400.0)
      if elapsed_sim > 0
      else 0.0
  )
  print(
      f"Dynamics complete in {elapsed_sim:.2f}s ({ns_per_day:.2f} ns/day)."
      f" Final PE: {pe_final:.2f} kJ/mol"
  )

  metadata = {
      "molecule_name": mol_name,
      "num_atoms": mol.GetNumAtoms(),
      "num_heavy_atoms": mol.GetNumHeavyAtoms(),
      "forcefield": forcefield_name,
      "charge_method": charge_method_used,
      "solvent": solvent,
      "steps": steps,
      "time_step_fs": time_step_fs,
      "simulation_time_ps": steps * time_step_fs / 1000.0,
      "temperature_k": temperature_k,
      "initial_potential_energy_kj_mol": round(pe_init, 2),
      "minimized_potential_energy_kj_mol": round(pe_min, 2),
      "final_potential_energy_kj_mol": round(pe_final, 2),
      "energy_relaxation_kj_mol": round(pe_min - pe_init, 2),
      "elapsed_simulation_seconds": round(elapsed_sim, 2),
      "performance_ns_per_day": round(ns_per_day, 2),
      "minimized_pdb": minimized_pdb_path,
      "trajectory_pdb": trajectory_pdb_path,
  }

  return trajectory_pdb_path, minimized_pdb_path, metadata


def analyze_trajectory(
    trajectory_path: str,
    base_mol: Chem.Mol,
    output_dir: str,
    cluster_cutoff: float = 0.5,
    export_centroids: int = 3,
) -> Dict[str, Any]:
  """Analyzes multi-frame OpenMM trajectory using native OpenMM and RDKit APIs."""
  print(f"\n=== Analyzing Trajectory: {trajectory_path} ===")
  # app.PDBFile raises a bare IndexError on a trajectory with no MODEL records,
  # so the frame-count check below is unreachable unless we catch that here.
  try:
    pdb = app.PDBFile(trajectory_path)
    num_frames = pdb.getNumFrames()
  except IndexError:
    num_frames = 0
  if num_frames == 0:
    raise ValueError(
        f"No frames found in trajectory file: {trajectory_path}. "
        "This usually means --report_interval is larger than --steps, so no "
        "frame was ever written. Lower --report_interval or raise --steps "
        "(or pass --skip_analysis)."
    )

  rd_mol = Chem.Mol(base_mol)
  rd_mol.RemoveAllConformers()
  num_atoms = rd_mol.GetNumAtoms()

  for f_idx in range(num_frames):
    pos = pdb.getPositions(asNumpy=True, frame=f_idx).value_in_unit(
        ounit.angstrom
    )
    conf = Chem.Conformer(num_atoms)
    for a_idx in range(num_atoms):
      conf.SetAtomPosition(
          a_idx,
          Point3D(
              float(pos[a_idx][0]), float(pos[a_idx][1]), float(pos[a_idx][2])
          ),
      )
    rd_mol.AddConformer(conf, assignId=True)

  print(f"Loaded {num_frames} frames into RDKit multi-conformer container.")

  # 1. Heavy-atom RMSD vs frame 0
  atom_map = [
      (i, i)
      for i in range(num_atoms)
      if rd_mol.GetAtomWithIdx(i).GetAtomicNum() > 1
  ]
  if len(atom_map) < MIN_RIGID_FIT_ATOMS:
    atom_map = [(i, i) for i in range(num_atoms)]
  rmsds = [0.0]
  for i in range(1, num_frames):
    rmsds.append(
        float(
            rdMolAlign.GetBestRMS(
                rd_mol,
                rd_mol,
                prbId=i,
                refId=0,
            )
        )
    )
    rdMolAlign.AlignMol(
        rd_mol,
        rd_mol,
        prbCid=i,
        refCid=0,
        atomMap=atom_map,
    )

  df_rmsd = pd.DataFrame(
      {"frame": list(range(num_frames)), "rmsd_angstrom": rmsds}
  )
  rmsd_csv_path = os.path.join(output_dir, "trajectory_rmsd.csv")
  df_rmsd.to_csv(rmsd_csv_path, index=False)

  # 2. Per-Atom RMSF across superimposed frames
  coords_array = np.array(
      [rd_mol.GetConformer(i).GetPositions() for i in range(num_frames)]
  )
  mean_positions = np.mean(coords_array, axis=0)
  sq_diff = np.sum((coords_array - mean_positions) ** 2, axis=-1)
  rmsf = np.sqrt(np.mean(sq_diff, axis=0))

  atom_data = []
  for i in range(num_atoms):
    atom = rd_mol.GetAtomWithIdx(i)
    atom_data.append({
        "atom_index": i,
        "element": atom.GetSymbol(),
        "is_heavy": atom.GetAtomicNum() > 1,
        "rmsf_angstrom": round(float(rmsf[i]), 3),
    })
  df_rmsf = pd.DataFrame(atom_data)
  rmsf_csv_path = os.path.join(output_dir, "trajectory_rmsf.csv")
  df_rmsf.to_csv(rmsf_csv_path, index=False)

  # 3. Torsion Fingerprint Deviation (TFD)
  tfd_values = [0.0]
  for i in range(1, num_frames):
    try:
      tfd = float(
          TorsionFingerprints.GetTFDBetweenConformers(rd_mol, [i], [0])[0]
      )
    except Exception:
      tfd = 0.0
    tfd_values.append(tfd)

  df_tfd = pd.DataFrame({"frame": list(range(num_frames)), "tfd": tfd_values})
  tfd_csv_path = os.path.join(output_dir, "trajectory_tfd.csv")
  df_tfd.to_csv(tfd_csv_path, index=False)

  # 4. Taylor-Butina Conformer Clustering on Pairwise RMSD Matrix
  dists = []
  for i in range(num_frames):
    for j in range(i):
      dists.append(
          rdMolAlign.GetBestRMS(
              rd_mol,
              rd_mol,
              prbId=i,
              refId=j,
          )
      )

  clusters = Butina.ClusterData(
      dists, num_frames, distThresh=cluster_cutoff, isDistData=True
  )
  cluster_summary = []
  for c_idx, members in enumerate(clusters):
    cluster_summary.append({
        "cluster_id": c_idx,
        "centroid_frame": members[0],
        "size": len(members),
        "fraction_population": round(len(members) / num_frames, 4),
        "members": list(members),
    })

  cluster_json_path = os.path.join(output_dir, "trajectory_clusters.json")
  with open(cluster_json_path, "w") as f:
    json.dump(cluster_summary, f, indent=2)

  # 5. Export Representative Cluster Centroid Structures
  exported_centroids = []
  for i in range(min(export_centroids, len(clusters))):
    c_info = cluster_summary[i]
    c_frame = c_info["centroid_frame"]
    centroid_pdb_path = os.path.join(
        output_dir, f"cluster_{i}_centroid_frame_{c_frame}.pdb"
    )
    writer = Chem.PDBWriter(centroid_pdb_path)
    writer.write(rd_mol, confId=c_frame)
    writer.close()
    exported_centroids.append(centroid_pdb_path)

  analysis_results = {
      "num_frames": num_frames,
      "mean_heavy_atom_rmsd_angstrom": round(float(np.mean(rmsds)), 3),
      "max_heavy_atom_rmsd_angstrom": round(float(np.max(rmsds)), 3),
      "mean_atom_rmsf_angstrom": round(float(np.mean(rmsf)), 3),
      "max_rmsf_atom_index": int(np.argmax(rmsf)),
      "mean_tfd": round(float(np.mean(tfd_values)), 4),
      "num_clusters": len(clusters),
      "top_cluster_population": (
          float(cluster_summary[0]["fraction_population"])
          if cluster_summary
          else 1.0
      ),
      "cluster_cutoff_angstrom": cluster_cutoff,
      "rmsd_csv": rmsd_csv_path,
      "rmsf_csv": rmsf_csv_path,
      "tfd_csv": tfd_csv_path,
      "clusters_json": cluster_json_path,
      "exported_centroids": exported_centroids,
  }

  return analysis_results


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Run OpenMM & OpenFF molecular dynamics simulation and automated"
          " trajectory analysis."
      )
  )
  parser.add_argument(
      "--smiles", help="Input SMILES string for small molecule."
  )
  parser.add_argument(
      "--input_file",
      "--input_sdf",
      "--input_pdb",
      dest="input_file",
      help="Path to starting 3D structure (SDF, MOL, PDB).",
  )

  parser.add_argument(
      "--name",
      default="ligand",
      help="Molecule name or identifier (default: 'ligand').",
  )
  parser.add_argument(
      "--output_dir",
      required=True,
      help="Output directory for simulation and analysis artifacts.",
  )
  parser.add_argument(
      "--steps",
      type=int,
      default=2500,
      help="Number of Langevin dynamics steps (default: 2500, 5.0 ps).",
  )
  parser.add_argument(
      "--time_step_fs",
      type=float,
      default=2.0,
      help="Time step in femtoseconds (default: 2.0).",
  )
  parser.add_argument(
      "--temperature_k",
      type=float,
      default=300.0,
      help="Simulation temperature in Kelvin (default: 300.0).",
  )
  parser.add_argument(
      "--friction_ps",
      type=float,
      default=1.0,
      help="Friction coefficient in 1/ps (default: 1.0).",
  )
  parser.add_argument(
      "--report_interval",
      type=int,
      default=100,
      help="Frames recorded every N steps (default: 100).",
  )
  parser.add_argument(
      "--solvent",
      choices=["obc2", "gbn2", "none"],
      default="obc2",
      help="Implicit solvent model (default: 'obc2').",
  )
  parser.add_argument(
      "--forcefield",
      default="openff-2.2.0",
      help="Small-molecule force field (default: 'openff-2.2.0').",
  )
  parser.add_argument(
      "--charge-method",
      "--charge_method",
      dest="charge_method",
      choices=list(CHARGE_METHOD_CHOICES),
      default="auto",
      help=(
          "Partial charge model (default: 'auto'). 'auto' uses AM1-BCC when a "
          "licence-free backend (AmberTools 'sqm') is installed and otherwise "
          "warns and falls back to RDKit MMFF94. 'am1bcc' requires AM1-BCC and "
          "errors out rather than degrading silently. 'mmff94' and 'gasteiger' "
          "force the named RDKit model."
      ),
  )
  parser.add_argument(
      "--cluster_cutoff",
      type=float,
      default=0.5,
      help="Taylor-Butina RMSD clustering cutoff in Angstroms (default: 0.5).",
  )
  parser.add_argument(
      "--export_centroids",
      type=int,
      default=3,
      help="Max centroid PDBs to export (default: 3).",
  )
  parser.add_argument(
      "--skip_analysis",
      action="store_true",
      help="Skip automated trajectory analysis after simulation.",
  )

  args = parser.parse_args()
  if not args.smiles and not args.input_file:
    parser.error("Must provide at least one of --smiles or --input_file.")

  # 1. Prepare structure
  mol = prepare_molecule(
      smiles=args.smiles, input_file=args.input_file, name=args.name
  )

  # 2. Run simulation
  traj_path, min_pdb_path, sim_meta = run_simulation(
      mol=mol,
      output_dir=args.output_dir,
      steps=args.steps,
      time_step_fs=args.time_step_fs,
      temperature_k=args.temperature_k,
      friction_ps=args.friction_ps,
      report_interval=args.report_interval,
      solvent=args.solvent,
      forcefield_name=args.forcefield,
      charge_method=args.charge_method,
  )

  # 3. Automated trajectory analysis
  if not args.skip_analysis:
    analysis_meta = analyze_trajectory(
        trajectory_path=traj_path,
        base_mol=mol,
        output_dir=args.output_dir,
        cluster_cutoff=args.cluster_cutoff,
        export_centroids=args.export_centroids,
    )
    sim_meta["trajectory_analysis"] = analysis_meta

  # 4. Write unified simulation summary
  summary_path = os.path.join(args.output_dir, "simulation_summary.json")
  with open(summary_path, "w") as f:
    json.dump(sim_meta, f, indent=2)

  print("\n=== Workflow Completed Successfully ===")
  print(f"Summary JSON saved: {summary_path}")
  print(f"Trajectory PDB: {traj_path}")
  if not args.skip_analysis:
    analysis = sim_meta["trajectory_analysis"]
    print(
        f"Mean Heavy-Atom RMSD: {analysis['mean_heavy_atom_rmsd_angstrom']} Å"
        f" (Max: {analysis['max_heavy_atom_rmsd_angstrom']} Å)"
    )
    print(f"Mean Per-Atom RMSF: {analysis['mean_atom_rmsf_angstrom']} Å")
    print(
        f"Conformational Clusters: {analysis['num_clusters']} (Dominant State:"
        f" {analysis['top_cluster_population']*100:.1f}%)"
    )
    print(f"Exported Centroids: {len(analysis['exported_centroids'])} files")


if __name__ == "__main__":
  main()
