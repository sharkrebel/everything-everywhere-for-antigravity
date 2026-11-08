# Relaxed Dihedral PES Scanning, Conformational Strain & Docking

This guide covers physical modeling of rotational barriers, bound conformational
strain ($\Delta E_{\text{strain}}$), and protein-ligand docking using the
production engines in `scripts/`.

--------------------------------------------------------------------------------

## 1. Relaxed Dihedral PES Scanning vs. Rigid 1D Sweeps

### The Physics of Rotamer Jamming

A rigid 1D scan rotates a dihedral angle $\phi$ while holding all other
molecular coordinates frozen. In flexible systems (biaryls, amides, substituted
rings), this rigid rotation drives adjacent substituents into severe steric
overlap, generating **unphysical energy spikes of $+30\text{ to }+100\text{
kcal/mol}$**.

*   **Rigid 1D Sweep (Unphysical)**: Holds orthogonal coordinates frozen while
    stepping the target dihedral angle $\phi$, causing artificial steric clash
    spikes ($+30\text{ to }+100\text{ kcal/mol}$).
*   **Relaxed Dihedral Scan (Thermodynamic)**: Minimizes all orthogonal degrees
    of freedom under a harmonic dihedral restraint at each step, propagating
    wavefront conformers from adjacent relaxed minima to recover the true
    rotational barrier ($\sim 2\text{--}5\text{ kcal/mol}$).

### Production CLI (`torsion_drive_scanner.py`):

Select the dihedral either by four explicit zero-based atom indices
(`--dihedral`) or by a SMARTS pattern plus the indices *within* that match
(`--smarts` + `--smarts_indices`). Use `--list_rotors` first if you do not know
the indices.

```bash
uv run scripts/torsion_drive_scanner.py \
    --smiles "c1ccccc1-c2ccccc2" \
    --list_rotors
```

```bash
uv run scripts/torsion_drive_scanner.py \
    --smiles "c1ccccc1-c2ccccc2" \
    --dihedral 4 5 6 7 \
    --algorithm bidirectional \
    --step_deg 10.0 \
    --force_field MMFF94s \
    --output_csv scratch/biaryl_pes_scan.csv \
    --plot_png scratch/biaryl_energy_profile.png
```

*   **Bidirectional Wavefront Minimization** (`--algorithm bidirectional`):
    Propagates coordinates in both clockwise ($0^\circ \to 360^\circ$) and
    counter-clockwise directions to quantify hysteresis and eliminate artificial
    barriers.
*   **Outputs**: Tabulated angles, relative energies ($\Delta E$ in
    $\text{kcal/mol}$), rotational barrier height, and an energy profile plot
    (`--plot_png`). Optional `--output_json` and a multi-conformer trajectory
    `--output_sdf`.

--------------------------------------------------------------------------------

## 2. Bound Conformational Strain ($\Delta E_{\text{strain}}$) & Preorganization

Rigidification and preorganization improve binding affinity by reducing both
bound conformational strain and the entropic penalty of binding: $$\Delta
G_{\text{eff}} = \Delta G_{\text{dock}} + \Delta E_{\text{strain}} - T\Delta
S_{\text{rot}}$$ where:

*   **Conformational Strain**: $\Delta E_{\text{strain}} = E_{\text{bound}} -
    E_{\text{free\_min}}$ (bound conformer local-minimum energy minus the
    lowest-energy free conformer). Strains $< 3.0\text{ kcal/mol}$ are readily
    tolerated; strains $> 5.0\text{ kcal/mol}$ severely penalize binding.
*   **Rotational Entropy Penalty**: $-T\Delta S_{\text{rot}} \approx +0.65\text{
    kcal/mol} \times \text{RotB}$ (loss of conformational freedom upon binding).

### Bound strain from a 3D pose (`molecular_docking_vina.py --relax-poses`)

$\Delta E_{\text{strain}}$ requires a 3D bound pose and is computed by the
docking engine's two-stage MMFF94s pose-relaxation workflow:

```bash
uv run scripts/molecular_docking_vina.py \
    --relax-poses scratch/docking_results/poses.sdf \
    --output-sdf scratch/poses_relaxed.sdf \
    --output-csv scratch/strain_breakdown.csv
```

*   Performs constrained local MMFF94s minimization of each bound pose and
    compares its energy against the global minimum of a free ETKDGv3 conformer
    ensemble of the same stereoisomer.

--------------------------------------------------------------------------------

## 3. Standalone pH-Dependent Protonation (`protonate_ligands.py` & `protonate_receptor.py`)

Crystal structures from the PDB typically lack explicit hydrogens, and raw
ligand SMILES strings are usually written in neutral unionized forms. Docking or
running force-field simulations without pH-dependent protonation omits explicit
polar hydrogen-bond donors (`HD`) on the protein receptor and misses critical
salt-bridge and charge-assisted hydrogen-bond interactions on acidic/basic
ligands.

Always protonate small molecules and protein receptors **before** invoking
downstream docking (`molecular_docking_vina.py`) or simulation scripts. Down-
stream consumers take an already-protonated system as input rather than calling
protonation code internally.

> **User-Specified Alternative Protonation Methods**: By default, this skill
> uses Open Babel (`protonate_ligands.py`) for ligands and OpenMM + PDBFixer
> (`protonate_receptor.py`) for receptors. However, the user may ask to use an
> alternative protonation method or external tool (such as Dimorphite-DL,
> Schrödinger Epik, ChemAxon Marvin `cxcalc`, OpenEye `fixpka` / `quacpac`,
> MolGpka, PROPKA, PDB2PQR, `reduce`, or AmberTools `tleap`). When the user
> requests an alternative method, run that tool directly to produce the
> protonated `.sdf`/`.csv` and `.pdb`/`.pdbqt` files, then pass those files to
> `molecular_docking_vina.py` or other downstream scripts.

### Small-Molecule Ligand Protonation (`protonate_ligands.py`)

Assigns physiological ionization states at any target pH using **Open Babel
`OBPhModel::CorrectForPH`** (`phmodel.txt` empirical SMARTS transformation
rules; O'Boyle et al. 2011) while mapping formal charges and explicit hydrogens
back onto the RDKit molecule to preserve input 3D coordinates, stereochemistry,
and Kekulization:

```bash
uv run scripts/protonate_ligands.py \
    --input ligands.csv \
    --smiles-col SMILES \
    --id-col compound_id \
    --ph 7.4 \
    --embed-3d \
    --output-sdf scratch/ligands_ph74.sdf \
    --output-csv scratch/ligands_ph74.csv
```

*   `--input` accepts a single SMILES string, `.csv`, `.sdf`/`.mol`, or `.smi`
    file.
*   `--enumerate-microstates` (with `--ph-window 1.5` and `--max-microstates 8`)
    enumerates near-$\text{p}K_a$ ionization microstates by sweeping the target
    pH across $[\text{pH} - \Delta\text{pH},\; \text{pH} + \Delta\text{pH}]$.
*   `--output-pdbqt` exports ligand `.pdbqt` file(s) with explicit polar `HD`
    atoms via `meeko`.
*   `--keep-input-charges` preserves pre-existing formal charges instead of
    neutralizing before re-ionizing; `--no-desalt` retains counterions.

### Protein Receptor Protonation (`protonate_receptor.py`)

Strips waters (`--keep-waters` to retain), non-protein HETATMs (`--keep-hetatm`
to retain), and alternate conformations (`altLoc`), repairs missing terminal or
heavy sidechain atoms with `PDBFixer` (without fabricating missing loops), and
assigns pH-dependent residue protonation states (His `HID`/`HIE`/`HIP`, Asp,
Glu, Lys, Cys, Tyr) via **OpenMM `Modeller.addHydrogens(pH=ph)`**
(`amber14-all.xml` by default, with a deterministic geometric polar-`HD`
fallback):

```bash
uv run scripts/protonate_receptor.py \
    --receptor 1iep \
    --ligand-resname STI \
    --ph 7.4 \
    --output-pdb scratch/receptor_ph74.pdb \
    --output-pdbqt scratch/receptor_ph74.pdbqt \
    --output-json scratch/receptor_ph74_summary.json
```

*   `--receptor` accepts a local `.pdb`/`.cif` file or a 4-letter PDB ID.
*   `--ligand-resname` computes and reports the co-crystal ligand pocket
    centroid (`pocket_center_xyz` in `--output-json`) before stripping HETATMs.
*   `--output-pdbqt` writes an AutoDock4/Vina receptor `.pdbqt` containing
    explicit polar hydrogen-bond donors (`HD`) and Gasteiger-style partial
    charges.

--------------------------------------------------------------------------------

## 4. Protein-Ligand Docking (`molecular_docking_vina.py`)

`molecular_docking_vina.py` consumes an **already-protonated receptor** (`.pdb`
or `.pdbqt`) and **already-protonated ligands** (`.sdf`, `.mol`, `.csv`, `.smi`,
or a protonated SMILES string) and evaluates AutoDock Vina / smina binding
affinities alongside matched-state bound conformational strain ($\Delta
E_{\text{strain}}$):

```bash
# Step 1: Protonate receptor and ligands (or use a user-requested alternative tool)
uv run scripts/protonate_receptor.py \
    --receptor 1iep \
    --ligand-resname STI \
    --ph 7.4 \
    --output-pdb scratch/receptor_ph74.pdb \
    --output-json scratch/receptor_ph74_summary.json

uv run scripts/protonate_ligands.py \
    --input "Cc1ccc(Nc2nccc(-c3cccnc3)n2)cc1NC(=O)c4ccc(CN5CCN(C)CC5)cc4" \
    --ph 7.4 \
    --embed-3d \
    --output-sdf scratch/ligands_ph74.sdf

# Step 2: Dock the pre-protonated system
uv run scripts/molecular_docking_vina.py \
    --receptor scratch/receptor_ph74.pdb \
    --ligands scratch/ligands_ph74.sdf \
    --center 15.2 28.4 -5.1 \
    --box_size 20.0 \
    --exhaustiveness 8 \
    --engine smina \
    --output_dir scratch/docking_results/
```

*   `--receptor` accepts a pre-protonated `.pdb` or `.pdbqt` path (or a 4-letter
    PDB ID for quick inspection; `prepare_receptor` warns on `stderr` if the
    receptor contains zero explicit hydrogens).
*   `--ligands` accepts a pre-protonated `.sdf`/`.mol` path (preserving exact
    formal charges, protonation states, and 3D conformers), a `.csv` path
    (preferring `protonated_smiles` if present, or `--smiles_col`), a `.smi`
    path, or a protonated SMILES string. Batch docking parallelizes over
    `--workers`.
*   Pocket centering: either pass `--center X Y Z` explicitly (e.g., from
    `protonate_receptor.py --output-json`'s `pocket_center_xyz`), or pass
    `--ligand_resname LIG` if the receptor PDB still contains the co-crystal
    ligand.
*   `--box_size` is a single side length in Angstroms (cubic box).
*   Converts receptor and ligand into PDBQT files (preserving explicit polar
    `HD` atoms and net formal charges), evaluates AutoDock Vina / smina binding
    affinities ($\text{kcal/mol}$), records `net_charge` alongside `smiles`, and
    computes matched-state bound ligand strain ($\Delta E_{\text{strain}}$).
