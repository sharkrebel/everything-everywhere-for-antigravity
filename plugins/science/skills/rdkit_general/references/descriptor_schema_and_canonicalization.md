# Molecular Descriptor Schema, API Directory, & Column Canonicalization

This reference defines the RDKit descriptor API directory (and common
`AttributeError` traps), the 13-property canonical column naming schema, the
`DescriptorCanonicalizer` engine, and `descriptors.py` recipes.

--------------------------------------------------------------------------------

## 1. Molecular Descriptor API Directory & Hallucination Traps

Conflating property abbreviations (`HBD`, `HBA`, `MW`, `TPSA`) with RDKit module
names (`Lipinski` vs `Descriptors` vs `rdMolDescriptors`) raises
`AttributeError`. Use these exact calls (or call
`compute_standard_descriptors(mol)` from `descriptors.py` to compute all 13 at
once):

-   **H-Bond Donors (`HBD`)**: `Lipinski.NumHDonors(mol)` or
    `rdMolDescriptors.CalcNumHBD(mol)` — *NOT* `Lipinski.NumHBD(mol)` or
    `Descriptors.NumHBD(mol)`.
-   **H-Bond Acceptors (`HBA`)**: `Lipinski.NumHAcceptors(mol)` or
    `rdMolDescriptors.CalcNumHBA(mol)` — *NOT* `Lipinski.NumHBA(mol)` or
    `Descriptors.NumHBA(mol)`.
-   **Rotatable Bonds (`RotB`)**: `Lipinski.NumRotatableBonds(mol)` or
    `rdMolDescriptors.CalcNumRotatableBonds(mol)` — *NOT*
    `mol.GetNumRotatableBonds()` or `Descriptors.RotatableBonds(mol)`.
-   **Molecular Weight (`MW`)**: `Descriptors.MolWt(mol)` — *NOT*
    `Descriptors.MW(mol)` or `mol.GetMolWt()`.
-   **Lipophilicity (`cLogP`)**: `Crippen.MolLogP(mol)` or
    `Descriptors.MolLogP(mol)` — *NOT* `Descriptors.cLogP(mol)` or
    `mol.GetLogP()`.
-   **Polar Surface Area (`TPSA`)**: `rdMolDescriptors.CalcTPSA(mol)` — *NOT*
    `Descriptors.TPSA(mol)`.
-   **Fraction Csp3 (`Fsp3`)**: `rdMolDescriptors.CalcFractionCSP3(mol)` — *NOT*
    `Lipinski.FractionCSP3(mol)`.
-   **Heavy Atoms (`HeavyAtoms`)**: `mol.GetNumHeavyAtoms()` — *NOT*
    `Descriptors.NumHeavyAtoms(mol)`.
-   **Ring Count (`Rings`)**: `Lipinski.RingCount(mol)` or
    `rdMolDescriptors.CalcNumRings(mol)` — *NOT* `mol.GetRingCount()`.
-   **Aromatic Rings (`AromaticRings`)**:
    `rdMolDescriptors.CalcNumAromaticRings(mol)` — *NOT*
    `Lipinski.AromaticRingCount(mol)`.
-   **Aliphatic Rings (`AliphaticRings`)**:
    `rdMolDescriptors.CalcNumAliphaticRings(mol)` — *NOT*
    `Lipinski.AliphaticRingCount(mol)`.
-   **Drug-Likeness (`QED`)**: `QED.qed(mol)` — *NOT* `Descriptors.QED(mol)`.

--------------------------------------------------------------------------------

## 2. Canonical 13-Property Column Naming Schema

To eliminate `KeyError` exceptions when joining tables or feeding featurized
DataFrames into QSAR models, all workflows use one of two canonical column
styles:

-   **Canonical Short Style (Default for Tables & Plots)**: `MW`, `cLogP`,
    `TPSA`, `HBD`, `HBA`, `RotB`, `Fsp3`, `HeavyAtoms`, `Rings`,
    `AromaticRings`, `AliphaticRings`, `QED`, `SA_Score`
-   **Canonical Snake Style (Alternative for Python Dicts & Backend Schemas)**:
    `molecular_weight`, `clogp`, `tpsa`, `hbd`, `hba`, `rotatable_bonds`,
    `fraction_csp3`, `heavy_atom_count`, `ring_count`, `aromatic_ring_count`,
    `aliphatic_ring_count`, `qed`, `sa_score`

### Multi-Tier Column Normalization Rules

When harmonizing external compound tables, normalize column headers via six
rules (implemented by `DescriptorCanonicalizer` in the `cheminformatics-core`
skill):

1.  **Bracket & Unit Stripping**: Removes units (`MW (Da)` $\to$ `MW`, `TPSA
    (Å²)` $\to$ `TPSA`, `Ki (nM)` $\to$ `Ki`).
2.  **Acronym-Preserving CamelCase Splitting**: Splits `PolarSurfaceArea` while
    preserving `MW`, `TPSA`, `QED`, `cLogP`, `RotB`.
3.  **Tool/Vendor Prefix Stripping**: Strips `rdkit_`, `calc_`, `calculated_`,
    `chembl_`, `pfizer_`, `pred_`, `prop_`, `lipinski_`.
4.  **Filler & Suffix Stripping**: Normalizes `number_of_rotatable_bonds`,
    `total_rings`, `hba_count`, `qed_score`.
5.  **Abbreviation Expansion**: Maps `wt` $\to$ `weight`, `arom` $\to$
    `aromatic`, `aliph` $\to$ `aliphatic`, `rot` $\to$ `rotatable`, `acc` $\to$
    `acceptor`, `don` $\to$ `donor`.
6.  **Target-Qualifying Non-Descriptor Invariant**: Metadata (`Compound_ID`,
    `Yield_%`, `target_chembl_id`, `assay_type`) and target-qualified
    bioactivity columns (`IC50_ALK`, `Ki_Carbonic_anhydrase_4`) are never
    collapsed to bare `IC50` or `Ki`.

--------------------------------------------------------------------------------

## 3. Direct RDKit Python Recipes & Scope Boundary

> [!IMPORTANT] **`compat.py` Scope Boundary**: `scripts/compat.py` is strictly
> an RDKit C++ version shim and plot theming module (`set_safe_plot_style`). It
> does **not** export `calculate_cns_mpo` or `compute_standard_descriptors`. For
> batch CSV featurization, CNS MPO, and SA score calculation, use the
> **`cheminformatics-core`** skill.

```python
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors

mol = Chem.MolFromSmiles("CC1=C(C=C(C=C1)C(=O)N2CCC(CC2)NC3=NC=NC4=C3C=CN4)C")

props_short = {
    "MW": round(Descriptors.MolWt(mol), 2),
    "cLogP": round(Crippen.MolLogP(mol), 2),
    "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
    "HBD": Lipinski.NumHDonors(mol),
    "HBA": Lipinski.NumHAcceptors(mol),
    "RotB": Lipinski.NumRotatableBonds(mol),
    "Fsp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
    "HeavyAtoms": mol.GetNumHeavyAtoms(),
    "Rings": Lipinski.RingCount(mol),
    "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
    "AliphaticRings": rdMolDescriptors.CalcNumAliphaticRings(mol),
    "QED": round(QED.qed(mol), 3),
}
```
