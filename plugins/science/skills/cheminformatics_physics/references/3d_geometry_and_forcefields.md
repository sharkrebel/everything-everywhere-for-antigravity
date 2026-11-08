# 3D Conformer Generation, Alignment, & Force Field Minimization

This reference covers 3D conformer ensemble generation (ETKDGv3), force field
energy minimization (MMFF94s/UFF), geometric restraints, and 3D alignment using
modern RDKit.

--------------------------------------------------------------------------------

## 1. Production CLI: Batch Conformer Generation (`generate_3d_conformers.py`)

When ionizable acids or bases are present, assign physiological protonation
states first (`uv run scripts/protonate_ligands.py --input ... --ph 7.4
--output-csv ...`) so electrostatic and hydrogen-bonding terms in MMFF94s
reflect the target pH.

Generate low-energy 3D conformers with RMSD pruning and force-field
minimization:

```bash
uv run scripts/generate_3d_conformers.py \
    --input_smiles "CC(C)Cc1ccc(cc1)C(C)C(=O)O" \
    --num_confs 50 \
    --prune_rmsd 0.5 \
    --force_field MMFF94s \
    --output_sdf scratch/conformers.sdf
```

*   Uses **ETKDGv3** with small-ring and macrocycle torsion rules.
*   Minimizes conformers using **MMFF94s** (with UFF fallback for non-MMFF
    elements like B, Si, metals).
*   Prunes redundant conformers within $0.5\text{ Å}$ RMSD threshold.

--------------------------------------------------------------------------------

## 2. ETKDG Parameters & Python Usage

```python
from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.AddHs(Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O"))

params = AllChem.ETKDGv3()
params.randomSeed = 42
params.pruneRmsThresh = 0.5            # Prune redundant conformers (< 0.5 Å RMSD)
params.useSmallRingTorsions = True     # Specialized potentials for 3-5 membered rings
params.useMacrocycleTorsions = True    # Specialized potentials for rings >= 10 atoms
params.numThreads = 0                  # Utilize all available CPU cores

cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=50, params=params))
```

--------------------------------------------------------------------------------

## 3. Force Field Minimization (MMFF94s vs. UFF)

*   **MMFF94s**: Recommended static parameterization for drug-like organic
    molecules (planar anilines, delocalized nitrogens).
*   **UFF**: Universal fallback across the periodic table when MMFF lacks
    parameters (e.g. Boron, metals).

```python
if AllChem.MMFFHasAllMoleculeParams(mol):
    results = AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant='MMFF94s', numThreads=0)
else:
    results = AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=0)
# results contains [(status, energy_kcal_mol), ...]
```

--------------------------------------------------------------------------------

## 4. Geometric Restraints & Constraints

Restraints are attached to the `ForceField` instance:

```python
import math
from rdkit.Chem import AllChem

ff = AllChem.UFFGetMoleculeForceField(mol, confId=0)

# 1. Freeze Atom Position (Fixed point)
ff.AddFixedPoint(0)

# 2. Distance Restraint (min_dist, max_dist, force_constant)
ff.AddDistanceConstraint(2, 5, 2.80, 3.20, 100.0)

# 3. Angle Restraint: MUST PASS RADIANS!
min_rad = math.radians(105.0)
max_rad = math.radians(115.0)
ff.AddAngleConstraint(2, 5, 6, min_rad, max_rad, 50.0)

# 4. Torsion / Dihedral Restraint: DEGREES, low force constant (0.01 - 0.1)
ff.UFFAddTorsionConstraint(2, 5, 6, 7, False, -5.0, 5.0, 0.05)

ff.Initialize()
ff.Minimize(maxIts=500)
```

--------------------------------------------------------------------------------

## 5. 3D Alignment & RMSD Calculation

```python
from rdkit.Chem import rdMolAlign

# Best RMSD with topological symmetry perception
best_rms = rdMolAlign.GetBestRMS(probe_mol, ref_mol, prbId=0, refId=0)

# Align all conformers to conformer 0
for cid in cids[1:]:
    rdMolAlign.AlignMol(mol, mol, prbCid=cid, refCid=0)
```

--------------------------------------------------------------------------------

## 6. Torsion Fingerprint Deviation (TFD) & Conformer Clustering (`torsion_fingerprint_deviation.py`)

TFD compares internal rotatable dihedral angles without being confounded by
global rigid-body translation/rotation or peripheral side-chain orientation:

```bash
uv run scripts/torsion_fingerprint_deviation.py \
    --smiles "CCCC(=O)OCC" \
    --num_confs 50 \
    --tfd_thresh 0.05 \
    --output_csv scratch/tfd_clusters.csv
```

--------------------------------------------------------------------------------

## 7. 3D Shape & Conformational Descriptors (`compute_3d_descriptors.py`)

Generates a free-conformer ensemble with ETKDGv3 (`useMacrocycleTorsions=True`
for rings $\ge 12$ atoms) and computes 3D shape descriptors
(`pmi1`/`pmi2`/`pmi3`, `radius_of_gyration`, `spherocity_index`, `asphericity`,
`inertial_shape_factor`, and rigid-rotor entropy penalty):

```bash
uv run scripts/compute_3d_descriptors.py \
    --smiles "CC(C)Cc1ccc(cc1)C(C)C(=O)O" \
    --num_confs 50 \
    --output_csv scratch/shape_descriptors.csv
```
