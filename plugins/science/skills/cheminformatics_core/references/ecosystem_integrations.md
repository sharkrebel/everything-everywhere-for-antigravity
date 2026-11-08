# RDKit Ecosystem Integrations: Datamol, Interactive Grids, & Open Babel

This guide details the integration of RDKit with ecosystem packages
(**Datamol**, **Interactive Grid Viewer**, and **Open Babel**), including format
conversions and accelerated multiprocessing.

--------------------------------------------------------------------------------

## 1. Ecosystem Feature & Capability Overview

-   **`datamol` (`dm`)**: Modern Pythonic RDKit API and multiprocessing
    (`dm.parallelized(fn, mols)`).
    -   *Fallback when not installed*: Native `concurrent.futures` or
        `rdMolStandardize.SuperParent`.
-   **`interactive_grid_viewer.py`**: Standalone 2D Chemical Grid + SAR Scatter
    Plot (`uv run ../rdkit_general/scripts/interactive_grid_viewer.py`).
    -   *Fallback when not installed*: Zero-dependency, pure inline SVG + HTML
        (CSP-compliant).
-   **`openbabel` (`pybel`)**: 100+ Chemical Format Interconversion
    (`format_converter_openbabel.py`).
    -   *Fallback when not installed*: Native RDKit I/O (SMILES, SDF, MOL, PDB).

--------------------------------------------------------------------------------

## 2. Interactive 2D Structure Grid Viewer (`interactive_grid_viewer.py`)

Generate self-contained, CSP-compliant 2D HTML chemical structure grids with
embedded SAR scatter plots (MW vs. cLogP), Lipinski/PAINS property badges, modal
zoom, and client-side filtering:

```bash
uv run ../rdkit_general/scripts/interactive_grid_viewer.py \
    --input scratch/candidates.csv \
    --smiles_col SMILES \
    --align_core "c1ccccc1" \
    --output "/absolute/path/to/artifacts/candidates_grid.html"
```

*   Generates zero-CDN, standalone HTML viewable directly in browsers or IDE
    previewers.
*   Mandatory deliverable whenever profiling compound sets ($\ge 5$ molecules).

--------------------------------------------------------------------------------

## 3. Open Babel Multi-Format Interconversion (`format_converter_openbabel.py`)

Convert between chemical formats (MOL2, XYZ, GRO, CJSON, PDBQT, SDF, SMILES):

```bash
uv run scripts/format_converter_openbabel.py \
    --input docking_pose.mol2 \
    --output docking_pose.pdbqt
```
