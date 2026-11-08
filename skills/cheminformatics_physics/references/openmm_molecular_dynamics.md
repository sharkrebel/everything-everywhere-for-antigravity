# OpenMM and OpenFF Molecular Dynamics for Small Molecules

A comprehensive reference guide for parameterizing small molecules using OpenFF
Sage force fields (`openff-2.2.0`), constructing OpenMM simulation systems with
implicit solvent models (OBC2/GBn2), and running energy minimization, Langevin
dynamics, and native trajectory analysis in pure Python.

--------------------------------------------------------------------------------

## 1. Overview & Architecture

Molecular dynamics (MD) and physics-based conformational relaxation provide
rigorous thermodynamic and conformational sampling beyond static empirical force
fields. To simulate arbitrary drug-like ligands and small-molecule complexes
without external binary compilation:

1.  **RDKit 3D Embedding & Charge Assignment**: Prepares 3D coordinates with
    explicit hydrogens, then assigns partial charges under an explicit policy
    (see §2): AM1-BCC via a licence-free backend when one is installed,
    otherwise RDKit MMFF94/Gasteiger with a warning.
2.  **OpenFF Molecule & SMIRNOFF Parameterization**: Converts RDKit molecules to
    OpenFF topologies and uses
    `openmmforcefields.generators.SMIRNOFFTemplateGenerator` to apply modern
    Sage force fields (`openff-2.2.0`).
3.  **OpenMM System & Implicit Solvent**: Applies implicit solvent (OBC2 / GBn2)
    for rapid relaxation without expensive periodic water boxes. Add a
    macromolecular force field (e.g. Amber14) only when the system actually
    contains a protein; a ligand-only system does not need it, and loading it
    registers residue templates that can shadow the SMIRNOFF generator.
4.  **Langevin Dynamics & Structure Output**: Propagates Newton's equations of
    motion with stochastic friction and thermal fluctuation via
    `LangevinMiddleIntegrator`, minimizing potential energy and exporting
    trajectory frames to PDB format.
5.  **Native OpenMM & RDKit Trajectory Analysis**: Reads multi-frame PDB
    trajectories directly via OpenMM and maps frames into an RDKit
    multi-conformer container to compute RMSD, per-atom RMSF, TFD dihedral
    changes, and Butina cluster centroids without external trajectory packages.

--------------------------------------------------------------------------------

## 2. Invariants & Execution Rules

> [!IMPORTANT] - **Charge Model Is a Declared Choice, Never a Silent One**:
> OpenFF Sage was parameterized against AM1-BCC charges. Prefer AM1-BCC
> (`off_mol.assign_partial_charges("am1bcc")`) whenever a **licence-free**
> backend is registered — that means AmberTools (`sqm`). Never make a
> licence-gated toolkit (e.g. OpenEye) a default. When no such backend exists,
> compute charges via RDKit (`MMFFGetMoleculeProperties` or
> `ComputeGasteigerCharges`), assign them with `numpy_array *
> openff.units.unit.elementary_charge`, and **emit a warning stating that the
> charge model does not match the force field's parameterization**. If the
> caller explicitly demanded AM1-BCC, raise instead of downgrading. - **Implicit
> Solvent Radii Are Approximate for Ligands**: The GB radii in `obc2.xml` /
> `gbn2.xml` are AMBER protein-derived and are applied to small molecules by
> element. Treat implicit solvation as a conformational restraint, not a
> quantitative ΔG_solv. - **Explicit Hydrogens Required**: Classical force
> fields require all explicit valence hydrogens. Always ensure `Chem.AddHs(mol)`
> has been called and all hydrogens have 3D coordinates before constructing
> OpenFF molecules or OpenMM topologies. - **Energy Minimization Before
> Dynamics**: Always perform energy minimization
> (`simulation.minimizeEnergy(maxIterations=1000)`) before propagating
> stochastic dynamics to resolve clashes and prevent numerical coordinates from
> exploding (NaNs). - **Native Trajectory Analysis Pipeline**: Analyze
> multi-frame simulation trajectories using OpenMM's native `PDBFile` reader and
> RDKit's alignment and clustering modules, avoiding external trajectory
> libraries.

--------------------------------------------------------------------------------

## 3. Modular Python 3 Recipes

### 3.1 3D Conformer Generation and Standalone Charge Calculation

```python
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from openff.units import unit as off_unit
from openff.toolkit.topology import Molecule as OffMolecule

def prepare_parameterized_openff_molecule(
    mol: Chem.Mol,
    charge_method: str = "mmff",
    random_seed: int = 42,
) -> tuple[OffMolecule, Chem.Mol]:
    """Prepares an OpenFF Molecule with explicit 3D coordinates and assigned partial charges.

    Args:
        mol: Input RDKit molecule (with or without hydrogens/3D conformers).
        charge_method: 'mmff' (MMFF94 charges with Gasteiger fallback) or 'gasteiger'.
        random_seed: Random seed for 3D conformer embedding.

    Returns:
        Tuple of (openff.toolkit.topology.Molecule, rdkit.Chem.Mol with hydrogens and 3D coords).
    """
    # 1. Ensure explicit hydrogens
    mol_h = Chem.AddHs(mol)

    # 2. Embed 3D conformer if not present
    if mol_h.GetNumConformers() == 0:
        params = AllChem.ETKDGv3()
        params.randomSeed = random_seed
        status = AllChem.EmbedMolecule(mol_h, params)
        if status != 0:
            AllChem.EmbedMolecule(mol_h, useRandomCoords=True, randomSeed=random_seed)
        AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)

    # 3. Standalone charge calculation (pure Python, zero binary dependencies)
    charges = None
    if charge_method.lower() == "mmff":
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol_h)
        if mmff_props is not None:
            charges = [mmff_props.GetMMFFPartialCharge(i) for i in range(mol_h.GetNumAtoms())]

    if charges is None:
        # Fallback to Gasteiger charges
        AllChem.ComputeGasteigerCharges(mol_h)
        charges = [float(atom.GetProp("_GasteigerCharge")) for atom in mol_h.GetAtoms()]
        charges = [0.0 if np.isnan(c) else c for c in charges]

    # 4. Convert to OpenFF Molecule
    off_mol = OffMolecule.from_rdkit(mol_h, allow_undefined_stereo=True)
    off_mol.partial_charges = np.array(charges) * off_unit.elementary_charge

    return off_mol, mol_h
```

### 3.2 OpenMM System Construction with OpenFF Sage & Implicit Solvent

```python
import openmm
import openmm.app as app
from openff.toolkit.topology import Molecule as OffMolecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

def create_openmm_system(
    off_mol: OffMolecule,
    forcefield_name: str = "openff-2.2.0",
    solvent_model: str = "obc2",
) -> tuple[app.Topology, openmm.System]:
    """Constructs an OpenMM System for a small molecule using SMIRNOFF and implicit solvent.

    Args:
        off_mol: Parameterized OpenFF Molecule.
        forcefield_name: SMIRNOFF force field name (e.g. 'openff-2.2.0').
        solvent_model: Implicit solvent XML identifier ('obc2' -> 'implicit/obc2.xml').

    Returns:
        Tuple of (OpenMM Topology, OpenMM System).
    """
    # 1. Register OpenFF SMIRNOFF template generator
    smirnoff_gen = SMIRNOFFTemplateGenerator(molecules=off_mol, forcefield=forcefield_name)

    # 2. Implicit solvent only - a ligand-only system needs no protein force field.
    solvent_xml = f"implicit/{solvent_model.lower()}.xml"
    ff = app.ForceField(solvent_xml)
    ff.registerTemplateGenerator(smirnoff_gen.generator)

    # 3. Create OpenMM Topology and System
    topology = off_mol.to_topology().to_openmm()
    system = ff.createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
    )

    return topology, system
```

### 3.3 Energy Minimization and Langevin Dynamics

```python
from rdkit import Chem
import openmm
import openmm.app as app
import openmm.unit as omm_unit

def run_minimization_and_dynamics(
    mol_h: Chem.Mol,
    topology: app.Topology,
    system: openmm.System,
    temperature_kelvin: float = 300.0,
    friction_per_ps: float = 1.0,
    time_step_fs: float = 2.0,
    max_minimization_iters: int = 1000,
    simulation_steps: int = 5000,
    output_pdb_path: str = None,
) -> dict:
    """Performs energy minimization and Langevin dynamics on an OpenMM system.

    Args:
        mol_h: RDKit molecule with 3D conformer matching topology.
        topology: OpenMM Topology.
        system: OpenMM System.
        temperature_kelvin: Simulation temperature in Kelvin.
        friction_per_ps: Langevin friction coefficient in 1/ps.
        time_step_fs: Integration time step in femtoseconds.
        max_minimization_iters: Max iterations for potential energy minimization.
        simulation_steps: Number of integration steps for production dynamics.
        output_pdb_path: Optional file path to export final relaxed structure as PDB.

    Returns:
        Dict containing initial and minimized potential energies (kJ/mol) and simulation stats.
    """
    # 1. Setup LangevinMiddleIntegrator
    temperature = temperature_kelvin * omm_unit.kelvin
    friction = friction_per_ps / omm_unit.picosecond
    time_step = time_step_fs * omm_unit.femtoseconds
    integrator = openmm.LangevinMiddleIntegrator(temperature, friction, time_step)

    # 2. Initialize Simulation and set positions from RDKit 3D coordinates
    simulation = app.Simulation(topology, system, integrator)
    conf = mol_h.GetConformer()
    positions = [
        openmm.Vec3(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z)
        for i in range(mol_h.GetNumAtoms())
    ] * omm_unit.angstrom
    simulation.context.setPositions(positions)

    # 3. Measure initial energy
    state_initial = simulation.context.getState(getEnergy=True)
    initial_pe = state_initial.getPotentialEnergy().value_in_unit(omm_unit.kilojoules_per_mole)

    # 4. Perform Energy Minimization
    simulation.minimizeEnergy(maxIterations=max_minimization_iters)
    state_minimized = simulation.context.getState(getEnergy=True, getPositions=True)
    minimized_pe = state_minimized.getPotentialEnergy().value_in_unit(omm_unit.kilojoules_per_mole)

    # 5. Run Langevin Dynamics
    if simulation_steps > 0:
        simulation.step(simulation_steps)
        state_final = simulation.context.getState(getEnergy=True, getPositions=True)
        final_pe = state_final.getPotentialEnergy().value_in_unit(omm_unit.kilojoules_per_mole)
        final_positions = state_final.getPositions()
    else:
        final_pe = minimized_pe
        final_positions = state_minimized.getPositions()

    # 6. Export relaxed structure to PDB if requested
    if output_pdb_path:
        with open(output_pdb_path, "w") as f:
            app.PDBFile.writeFile(topology, final_positions, f)

    return {
        "initial_potential_energy_kj_mol": round(initial_pe, 2),
        "minimized_potential_energy_kj_mol": round(minimized_pe, 2),
        "final_potential_energy_kj_mol": round(final_pe, 2),
        "energy_delta_kj_mol": round(minimized_pe - initial_pe, 2),
        "steps_completed": simulation_steps,
        "output_pdb": output_pdb_path,
    }
```

### 3.4 Native Trajectory Recording & RDKit Analysis (Zero External Trajectory Dependencies)

OpenMM natively writes multi-frame trajectories using `openmm.app.PDBReporter`,
and OpenMM's native `openmm.app.PDBFile` reads them back in without needing
external, unmaintained trajectory libraries . Converting trajectory frames
directly into conformers on a single RDKit `Mol` object enables heavy-atom RMSD
tracking, per-atom RMSF fluctuation profiles, Torsion Fingerprint Deviation
(TFD), and Butina conformer clustering.

#### Core Trajectory Metrics & Principles

1.  **Multi-Conformer Ingestion**: Load coordinates using OpenMM's native parser
    (`pdb = openmm.app.PDBFile(traj_path)`). Iterate across all frames
    (`pdb.getNumFrames()`) and append each frame's coordinates to an RDKit `Mol`
    as a distinct conformer (`Chem.Conformer`).
2.  **Heavy-Atom Root-Mean-Square Deviation (RMSD)**: Measures overall
    structural drift and equilibration across simulation time $t$ relative to a
    reference geometry (frame 0): $$\text{RMSD}(t) = \sqrt{\frac{1}{N}
    \sum_{i=1}^N \|\mathbf{r}_i(t) - \mathbf{r}_i(0)\|^2}$$ Computed with
    rigid-body superposition using `rdkit.Chem.rdMolAlign.GetBestRMS`.
3.  **Per-Atom Root-Mean-Square Fluctuation (RMSF)**: Quantifies the local
    positional flexibility of each atom $a$ across $T$ frames after rigid-body
    superposition: $$\text{RMSF}_a = \sqrt{\frac{1}{T}\sum_{t=1}^T
    \|\mathbf{r}_{a,t} - \bar{\mathbf{r}}_a\|^2}$$ where $\bar{\mathbf{r}}_a$ is
    the average coordinate of atom $a$. Distinguishes rigid aromatic or
    macrocyclic cores from flexible aliphatic chains and rotatable substituents.
4.  **Torsion Fingerprint Deviation (TFD)**: Measures dihedral rotamer
    transitions using `rdkit.Chem.TorsionFingerprints.GetTFDBetweenConformers`.
    Unlike Cartesian RMSD, TFD isolates internal conformational transitions
    independently of rigid-body translation or molecular size.
5.  **Taylor-Butina Sphere-Exclusion Clustering**: Constructs a condensed
    pairwise heavy-atom RMSD distance matrix and applies the Taylor-Butina
    sphere exclusion algorithm (`rdkit.ML.Cluster.Butina.ClusterData`) with a
    distance threshold (e.g. 0.5 Å). Automatically partitions continuous MD
    trajectories into discrete conformational clusters and identifies
    representative centroid frames for downstream docking or energy analysis.

```python
import numpy as np
import openmm.app as app
import openmm.unit as omm_unit
from rdkit import Chem
from rdkit.Chem import rdMolAlign, TorsionFingerprints
from rdkit.Geometry import Point3D
from rdkit.ML.Cluster import Butina

def analyze_openmm_pdb_trajectory(
    pdb_trajectory_path: str,
    base_mol: Chem.Mol,
    cluster_cutoff: float = 0.5,
) -> dict:
    """Loads a multi-frame OpenMM PDB trajectory into RDKit and performs conformational analysis."""
    pdb_traj = app.PDBFile(pdb_trajectory_path)
    num_frames = pdb_traj.getNumFrames()

    rd_mol = Chem.Mol(base_mol)
    rd_mol.RemoveAllConformers()
    num_atoms = rd_mol.GetNumAtoms()
    coords_list = []

    for f_idx in range(num_frames):
        pos_angstrom = pdb_traj.getPositions(asNumpy=True, frame=f_idx).value_in_unit(omm_unit.angstrom)
        coords_list.append(pos_angstrom)
        conf = Chem.Conformer(num_atoms)
        for a_idx in range(num_atoms):
            p = pos_angstrom[a_idx]
            conf.SetAtomPosition(a_idx, Point3D(float(p[0]), float(p[1]), float(p[2])))
        rd_mol.AddConformer(conf, assignId=True)

    # 1. Heavy-Atom RMSD vs Frame 0 & Superposition for RMSF
    atom_map = [(i, i) for i in range(num_atoms) if rd_mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
    # Uniquely constraining 3D rotation requires >= 3 points; fall back to all atoms for H2O/CH4/C2H6
    if len(atom_map) < 3:
        atom_map = [(i, i) for i in range(num_atoms)]
    rmsds = [0.0]
    for i in range(1, num_frames):
        rmsds.append(float(rdMolAlign.GetBestRMS(rd_mol, rd_mol, prbId=i, refId=0)))
        rdMolAlign.AlignMol(rd_mol, rd_mol, prbCid=i, refCid=0, atomMap=atom_map)

    # 2. Per-Atom RMSF (Root Mean Square Fluctuation) across superimposed frames
    coords_arr = np.array([rd_mol.GetConformer(i).GetPositions() for i in range(num_frames)])
    mean_coords = np.mean(coords_arr, axis=0)
    rmsf = np.sqrt(np.mean(np.sum((coords_arr - mean_coords)**2, axis=-1), axis=0))

    # 3. Torsion Fingerprint Deviation (TFD) vs Frame 0
    tfd = [TorsionFingerprints.GetTFDBetweenConformers(rd_mol, [i], [0])[0] for i in range(num_frames)]

    # 4. Butina Conformer Clustering on Pairwise RMSD Matrix
    dists = []
    for i in range(num_frames):
        for j in range(i):
            dists.append(rdMolAlign.GetBestRMS(rd_mol, rd_mol, prbId=i, refId=j))
    clusters = Butina.ClusterData(dists, num_frames, distThresh=cluster_cutoff, isDistData=True)

    return {
        "num_frames": num_frames,
        "mean_rmsd_angstrom": float(np.mean(rmsds)),
        "max_rmsd_angstrom": float(np.max(rmsds)),
        "rmsd_time_series": rmsds,
        "rmsf_per_atom": rmsf.tolist(),
        "tfd_time_series": tfd,
        "num_clusters": len(clusters),
        "centroid_frames": [c[0] for c in clusters],
    }
```

### 3.5 End-to-End Simulation & Trajectory Generation via `run_openmm_simulation.py`

To run a complete molecular dynamics simulation and automatically generate the
multi-frame trajectory along with RMSD, RMSF, TFD, and clustering metrics in a
single automated step:

```bash
# Run 2,000 steps of Langevin dynamics from SMILES with OBC2 implicit solvent
uv run scripts/run_openmm_simulation.py \
    --smiles "Cc1nn(C)c2c1-c1cnc(N)c(c1)O[C@H](C)c1cc(F)ccc1C(=O)N(C)C2" \
    --steps 2000 \
    --report_interval 100 \
    --solvent obc2 \
    --output_dir path/to/output_dir

# Or run from an existing 3D structure file (SDF, MOL, or PDB):
uv run scripts/run_openmm_simulation.py \
    --input_file ligand.sdf \
    --steps 5000 \
    --temperature_k 300.0 \
    --cluster_cutoff 0.5 \
    --output_dir path/to/output_dir
```

This pipeline automatically performs:

1.  Valence completion (`Chem.AddHs`) and 3D embedding (ETKDGv3 if needed).
2.  Charge assignment per `--charge-method` (default `auto`: AM1-BCC if a
    licence-free backend is present, else MMFF94 with a warning).
3.  OpenFF Sage parameterization (`openff-2.2.0`) + implicit solvent (OBC2 or
    GBn2).
4.  Potential energy minimization (up to 1,000 iterations).
5.  Langevin dynamics propagation (2.0 fs time step, 300 K, 1.0 ps⁻¹ friction).
6.  Real-time trajectory streaming to `trajectory.pdb`.
7.  Automatic conformational trajectory analysis producing all 6 analytical
    outputs (`trajectory_rmsd.csv`, `trajectory_rmsf.csv`, `trajectory_tfd.csv`,
    `trajectory_clusters.json`, `simulation_summary.json`, and centroid PDB
    conformers).

### 3.6 Standalone Offline Trajectory Analysis via `openmm_trajectory_analyzer.py`

If you already have a pre-existing multi-frame simulation PDB trajectory file,
analyze it directly:

```bash
uv run scripts/openmm_trajectory_analyzer.py \
    --trajectory path/to/trajectory.pdb \
    --output_dir path/to/output_dir \
    --cluster_cutoff 0.5 \
    --export_centroids 3
```

This produces:

-   `trajectory_rmsd.csv`: Frame-by-frame heavy-atom RMSD vs starting pose.
-   `trajectory_rmsf.csv`: Per-atom fluctuation profile (rigid vs flexible
    atoms).
-   `trajectory_tfd.csv`: Dihedral rotamer transitions across simulation time.
-   `trajectory_clusters.json`: Butina cluster memberships, populations, and
    centroid frame IDs.
-   `cluster_*_centroid_frame_*.pdb`: Exported representative conformer
    structures.

--------------------------------------------------------------------------------

## 4. Diagnostic & Troubleshooting Guide

-   **`IncompatibleTypeError: Expected openff.units.unit.Quantity`**
    -   **Root Cause**: Passing OpenMM `openmm.unit.Quantity` directly into
        `off_mol.partial_charges`
    -   **Remediation**: Use `openff.units.unit` (e.g. `import openff.units.unit
        as off_unit; off_mol.partial_charges = np.array(charges) *
        off_unit.elementary_charge`).
-   **`NoTemplateFoundError: Could not find template for residue`**
    -   **Root Cause**: OpenMM cannot recognize non-standard small-molecule
        residue
    -   **Remediation**: Ensure `SMIRNOFFTemplateGenerator(molecules=off_mol,
        forcefield='openff-2.2.0')` is registered with the `ForceField` via
        `ff.registerTemplateGenerator(generator.generator)`.
-   **`ValueError: Molecule must contain explicit hydrogens`**
    -   **Root Cause**: Attempting to parameterize a SMILES or mol without
        explicit hydrogens
    -   **Remediation**: Call `Chem.AddHs(mol)` and embed a 3D conformer before
        constructing `openff.toolkit.topology.Molecule`.
-   **`subprocess.CalledProcessError: antechamber command not found` /
    `ChargeMethodUnavailableError`**
    -   **Root Cause**: AM1-BCC needs a semi-empirical QM backend, and the only
        licence-free provider (AmberTools, which supplies `sqm`/`antechamber`)
        is not installed.
    -   **Remediation**: Install AmberTools (`conda install -c conda-forge
        ambertools`) to get the charge model Sage was parameterized against. If
        that is not possible, run with `--charge-method mmff94` and record in
        your report that charges are RDKit-derived and therefore approximate. Do
        **not** substitute a licence-gated toolkit as a default.
-   **`OpenMMException: Particle coordinate is NaN`**
    -   **Root Cause**: Extreme clashes in unminimized starting conformers
        causing numerical blowup
    -   **Remediation**: Always execute
        `simulation.minimizeEnergy(maxIterations=1000)` before initiating
        `simulation.step()`.
