# Standard AGY Science Skills Integration for Machine Learning & QSAR

This guide details interoperability between `cheminformatics_ml` and canonical
AGY science skills for predictive modeling workflows.

--------------------------------------------------------------------------------

## 1. Science Skills Integration Map for ML

-   ****`chembl-database`****
    -   **Role in Predictive Modeling & QSAR**: Bioactivity extraction &
        temporal tracking
    -   **Standard Command / Tool**: `chembl_api.py activity --filter
        target_chembl_id=...`
    -   **Output Used in ML**: Standard values (`pIC50`, `pKi`), `document_year`
        (for Sheridan temporal split)
-   ****`pubchem-database`****
    -   **Role in Predictive Modeling & QSAR**: Bioassay concentration &
        structure mapping
    -   **Standard Command / Tool**: `pubchem_api.py` / PUG-REST bioassays
    -   **Output Used in ML**: $pK_i$ / $pIC_{50}$ values with verified
        $\mu\text{M}$ unit conversion
-   ****`cheminformatics_core`****
    -   **Role in Predictive Modeling & QSAR**: 2D standardization & property
        profiling
    -   **Standard Command / Tool**: `standardize_molecule.py`,
        `chemical_vendor_resolver.py`
    -   **Output Used in ML**: Cleaned parent structures (`SuperParent`),
        verified CAS, MW, cLogP, TPSA
-   ****`cheminformatics_physics`****
    -   **Role in Predictive Modeling & QSAR**: 3D physics & thermodynamic
        descriptors
    -   **Standard Command / Tool**: `torsion_drive_scanner.py`,
        `molecular_docking_vina.py`
    -   **Output Used in ML**: Strain energy ($\Delta E_{\text{strain}}$),
        docking scores ($\Delta G_{\text{dock}}$), effective free energy
        ($\Delta G_{\text{eff}}$)
-   ****`pdb-database`****
    -   **Role in Predictive Modeling & QSAR**: Experimental structural
        templates
    -   **Standard Command / Tool**: `search_pdb.py`
    -   **Output Used in ML**: Co-crystal structures for docking grid boundaries
        and bound conformation strain

--------------------------------------------------------------------------------

## 2. Standard Predictive Modeling Pipeline

```
1. Extract Bioassay Data (chembl-database / pubchem-database)
   │
   ▼
2. Molecular Standardization (cheminformatics_core)
   │  - SuperParent: Salt stripping, charge neutralization, tautomer canonicalization
   ▼
3. Optional 3D Physical Descriptors (cheminformatics_physics)
   │  - Conformational strain (MMFF94s), docking affinities, dihedral rotational barriers
   ▼
4. Principled Dataset Partitioning & ML Benchmarking (cheminformatics_ml)
   │  - congeneric_series_ml.py
   │  - Sheridan prospective temporal split, nested cross-validation, or scaffold split
   ▼
5. Multi-Target Parity & Reporting (cheminformatics_ml)
   │  - plot_archetypes.py (multi-target parity grids, coefficient ranking)
   │  - report_builder.py (Markdown deliverable with verbatim metric concordance)
```

--------------------------------------------------------------------------------

## 3. Mandatory Invariants for Cross-Skill ML

*   **No Raw Crawlers**: Always use canonical skill scripts (`chembl_api.py`)
    rather than raw `urllib`/`curl` scripts.
*   **Unit Verification**: Verify whether assay endpoints are in $\text{nM}$
    ($pK_i = 9 - \log_{10}X$) or $\mu\text{M}$ ($pK_i = 6 - \log_{10}X$).
*   **Zero Leakage**: Feature computation, scaling, and hyperparameter tuning
    must occur strictly inside cross-validation loops or on the training split.
