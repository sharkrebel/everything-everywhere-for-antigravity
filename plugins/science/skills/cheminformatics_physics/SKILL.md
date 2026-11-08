---
name: cheminformatics-physics
description: >-
  Simulates 3D molecular physics, conformer ensembles, and CADD workflows:
  pH-dependent ligand and receptor protonation, ETKDGv3 3D conformer generation,
  dihedral/relaxed PES torsion scans, AutoDock Vina/smina docking,
  conformational ligand strain energy, PROTAC/macrocycle 3D geometry, and OpenMM
  molecular dynamics.
metadata:
  publisher: google
---

# RDKit Molecular Physics & 3D CADD Hub (`cheminformatics-physics`)

## Prerequisites

1.  **`rdkit-general`**: Read the `rdkit-general` skill before running any
    workflow here.
2.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH. Always execute skill scripts as `uv run
    /path/to/script.py ...` (or `python3 /path/to/script.py ...`); **never** run
    `uv run python3 /path/to/script.py`, as placing `python3` after `uv run`
    bypasses the script's PEP 723 `# /// script` dependency block.
3.  **User Notification**: If `.licenses/cheminformatics_physics_LICENSE.txt`
    does not already exist in the workspace root then (1) prominently notify the
    user of the terms of use for the tools and data sources this skill relies
    on, then (2) create the file recording the notification text and timestamp:
    -   **AutoDock Vina**: Apache 2.0 (Trott & Olson 2010; Eberhardt et
        al. 2021)
    -   **Open Babel**: GPLv2 (https://openbabel.org/; O'Boyle et al. 2011)
    -   **OpenMM & PDBFixer**: MIT and LGPL (Eastman et al. 2017)
    -   **OpenFF Sage**: Open Force Field Initiative (MIT)
    -   **RCSB PDB** (https://www.rcsb.org/pages/policies), used when
        `--receptor` is a 4-letter PDB ID

## Core Rules

-   **Shared foundation rules**: Follow all core rules in
    `../rdkit_general/SKILL.md`.
-   **Decoupled Protonation & User-Specified Alternatives**: Protonate receptors
    and small-molecule ligands **before** passing them to downstream docking or
    simulation scripts. By default, use `scripts/protonate_ligands.py` (Open
    Babel `OBPhModel::CorrectForPH`) for small molecules and
    `scripts/protonate_receptor.py` (OpenMM `Modeller` + `PDBFixer`) for protein
    receptors. However, the user may ask to use an alternative protonation
    method or external tool (e.g., Dimorphite-DL, Epik, Marvin `cxcalc`, OpenEye
    `fixpka`, PROPKA, PDB2PQR, `reduce`, or custom rules); when they do, run the
    user's requested protonation method to produce the protonated `.sdf`/`.csv`
    or `.pdb`/`.pdbqt` files first. Downstream scripts that consume protonated
    systems (`scripts/molecular_docking_vina.py`,
    `scripts/run_openmm_simulation.py`, etc.) accept an already-protonated
    system as input rather than internally calling protonation code.
-   **100% Series Coverage & Batch Verification**: Never arbitrarily downsample
    a chemical series when running 3D conformer generation, Vina docking, or
    strain calculation. Wait for all batch jobs to finish completely and verify
    that output CSV/SDF files contain all expected records before running
    downstream ML comparisons or writing reports.

--------------------------------------------------------------------------------

## Task-to-Domain Routing Directory

Match the user's objective to the domain below, **read the corresponding
reference guide**, and execute from `scripts/`:

1.  **Protonation, Conformers, PES Torsion Drives, Force Fields &
    Structure-Based Docking**

    -   **Reference Guides**:
        -   `references/3d_geometry_and_forcefields.md` — read when generating
            ETKDGv3 conformer ensembles, minimizing with MMFF94s/UFF, computing
            RMSD / Torsion Fingerprint Deviation (TFD) clusters, or profiling 3D
            shape descriptors (PMI, radius of gyration, 3D PSA).
        -   `references/relaxed_pes_strain_and_docking.md` — read when
            protonating small-molecule ligands or protein receptors at target
            pH, running wavefront relaxed PES torsion drives, computing bound
            conformational strain ($\Delta E_{\text{strain}}$), or docking with
            AutoDock Vina/smina.
    -   **Standard Executables**:
        -   `uv run scripts/protonate_ligands.py` — standalone pH-dependent
            small-molecule protonation (`--ph 7.4` via Open Babel
            `OBPhModel::CorrectForPH`), pH-window microstate enumeration
            (`--enumerate-microstates`), 3D ETKDGv3 embedding (`--embed-3d`),
            and `.sdf`/`.csv`/`.pdbqt` export. (If the user asks to use an
            alternative ligand protonation method, run that tool instead to
            produce the protonated `.sdf`/`.csv`.)
        -   `uv run scripts/protonate_receptor.py` — standalone pH-dependent
            protein/receptor protonation (`--ph 7.4` via OpenMM `Modeller` +
            `PDBFixer`), solvent/altLoc cleanup, pocket centroid calculation
            (`--ligand-resname`), and `.pdb`/`.pdbqt` export. (If the user asks
            to use an alternative receptor protonation method, run that tool
            instead to produce the protonated `.pdb`/`.pdbqt`.)
        -   `uv run scripts/molecular_docking_vina.py` — consumes an
            already-protonated receptor (`.pdb` or `.pdbqt`) and
            already-protonated ligands (`.sdf`, `.mol`, `.csv`, `.smi`, or
            protonated SMILES), runs Vina/smina docking, and computes
            matched-state bound conformational strain (`--workers 4` for batch
            series docking).
        -   `uv run scripts/generate_3d_conformers.py` — ETKDGv3 conformer
            ensembles with MMFF94s/UFF minimization and RMSD pruning.
        -   `uv run scripts/torsion_drive_scanner.py` — relaxed PES torsion
            drive scans by wavefront propagation.
        -   `uv run scripts/torsion_fingerprint_deviation.py` — TFD matrix and
            conformer clustering by internal dihedrals.
        -   `uv run scripts/compute_3d_descriptors.py` — 3D conformer geometry
            and shape descriptors (PMI, radius of gyration, spherocity,
            asphericity, inertial shape factor, and macrocycle conformer
            sampling).

2.  **OpenMM Molecular Dynamics & Trajectory Analysis**

    -   **Reference Guide**: `references/openmm_molecular_dynamics.md` — read
        when parameterizing ligands with OpenFF Sage (`openff-2.2.0`), selecting
        `--charge-method auto` vs. `am1bcc`, or running OpenMM MD and trajectory
        analysis.
    -   **Standard Executables**:
        -   `uv run scripts/run_openmm_simulation.py` — end-to-end
            small-molecule MD and physics relaxation via OpenMM and OpenFF Sage.
        -   `uv run scripts/openmm_trajectory_analyzer.py` — analyzes
            multi-frame trajectories (RMSD, RMSF, TFD, Butina clustering, and
            representative centroid export).

--------------------------------------------------------------------------------

## Progressive Execution Protocol

1.  **Inspect Domain Reference**: Call `view_file` on the target reference guide
    above to load CLI syntax, force-field flags, and charge-method rules.
2.  **Protonation, Structure & Charge Preparation**: Before docking or
    simulation, assign physiological protonation states at target pH (`--ph
    7.4`) using `protonate_ligands.py` (Open Babel) and `protonate_receptor.py`
    (OpenMM + PDBFixer), or using an alternative protonation method if requested
    by the user, then pass the resulting protonated `.sdf`/`.csv` and
    `.pdb`/`.pdbqt` files into the downstream consumer script. Ensure explicit
    hydrogens (`Chem.AddHs`) and 3D ETKDGv3 coordinates before force-field
    minimization or OpenMM parameterization (`--charge-method auto` uses
    AmberTools `sqm` AM1-BCC when installed and falls back to RDKit MMFF94 with
    an explicit warning; pass `--charge-method am1bcc` to require AM1-BCC).
3.  **Execute Simulation / Docking & Verify Outputs**: Run the designated script
    from `scripts/` and verify all output SDF/CSV/PDB files are populated on
    disk.
4.  **Summarize Findings in Markdown**: Summarize docking affinities ($\Delta
    G$), strain energies ($\Delta E_{\text{strain}}$), dihedral barriers, or MD
    trajectory metrics in the output directory using GFM tables and generated
    figures
    (`../rdkit_general/references/reporting_and_nomenclature_invariants.md`).
