# RDKit API Migration, Sanitization, & Molecular Graph Editing

This reference covers RDKit runtime API version migration (2022+ / 2024+),
`scripts/compat.py` wrappers, molecular sanitization recovery, and safe
`Chem.RWMol` editing.

--------------------------------------------------------------------------------

## 1. Cross-Version Discrepancies & Modern API Migration

Always prefer `scripts/compat.py` (`get_morgan_fingerprint`,
`standardize_molecule`, `get_etkdg_params`, `safe_substruct_match`,
`set_safe_plot_style`) instead of branching on the RDKit version in caller code.

-   **Morgan Fingerprints (RDKit 2022.09+)**:

    -   **Legacy (Deprecated)**: `AllChem.GetMorganFingerprintAsBitVect(mol, 2,
        nBits=2048)`
    -   **Modern Generator API**:

        ```python
        from rdkit.Chem import rdFingerprintGenerator
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        fp = gen.GetFingerprint(mol)       # or gen.GetCountFingerprint(mol)
        ```
-   **Molecular Standardization (RDKit 2021.03+)**:

    -   **Modern API**: `from rdkit.Chem.MolStandardize import rdMolStandardize`
        (`rdMolStandardize.SuperParent(mol)`, `Normalizer()`,
        `LargestFragmentChooser()`, `TautomerEnumerator()`).
-   **Substructure Match Parameters (RDKit 2020.03+)**:

    -   **Modern API**: `params = Chem.SubstructMatchParameters()`;
        `mol.GetSubstructMatches(query, params)`.
-   **Pandas Mol Representation (RDKit 2023.03+)**:

    -   **Modern API**: `PandasTools.molRepresentation = 'svg'`.
-   **CoordGen 2D Layout (RDKit 2020.09+)**:

    -   **Modern API**: `rdDepictor.SetPreferCoordGen(True)`.

--------------------------------------------------------------------------------

## 2. Molecular Sanitization, Valence Traps, & Ring Inspection

### Formal Charges on Nitro & N-Oxide Functional Groups

RDKit enforces strict octet rules during sanitization:

-   Nitro groups must be drawn with formal charges (`[N+](=O)[O-]`), not
    pentavalent `N(=O)=O`.
-   Pyridine N-oxides must be drawn `[n+]1ccccc1[O-]`, not `n1ccccc1=O`.

### Safe Property Cache Before Ring Inspection (`UpdatePropertyCache`)

Calling `Chem.GetSymmSSSR(mol)` or `mol.GetRingInfo()` on an unsanitized
molecule (`sanitize=False`) before computing implicit valences throws a C++
`Pre-condition Violation` exception. Always update the property cache first:

```python
mol.UpdatePropertyCache(strict=False)
Chem.GetSymmSSSR(mol)
```

### Partial Sanitization Recovery

When parsing molecules with non-standard aromaticity or nitrogen valences:

```python
from rdkit import Chem
flags = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
Chem.SanitizeMol(mol, sanitizeOps=flags)
```

--------------------------------------------------------------------------------

## 3. Safe Molecule Editing with `Chem.RWMol`

### Modifying Atoms or Bonds During Traversal (Index Shifting Bug)

Calling `rwmol.RemoveAtom(atom.GetIdx())` inside a forward loop invalidates
subsequent atom indices. Always collect target indices first and delete in
**strictly descending order**:

```python
rwmol = Chem.RWMol(mol)
to_remove = [atom.GetIdx() for atom in rwmol.GetAtoms() if atom.GetSymbol() == "Cl"]
for idx in sorted(to_remove, reverse=True):
    rwmol.RemoveAtom(idx)
Chem.SanitizeMol(rwmol)
```

--------------------------------------------------------------------------------

## 4. Chemical Reaction & Bulk I/O Mechanics

-   **SMARTS vs. SMILES**: `Chem.MolFromSmarts` creates query patterns without
    perceived valence. Never run 3D conformer generation or force fields on a
    molecule constructed via `MolFromSmarts`.
-   **Nested Tuples in `rxn.RunReactants`**: `rxn.RunReactants` returns a tuple
    of product sets (`((prodA_1, prodB_1), (prodA_2, prodB_2))`), and reaction
    products are unsanitized by default; call `Chem.SanitizeMol(p)` on every
    product molecule.
-   **Corrupted SDF Records**: `Chem.SDMolSupplier` and
    `Chem.ForwardSDMolSupplier` yield `None` on malformed records rather than
    raising Python exceptions; always guard loops with `if mol is None:
    continue`.
