---
name: rdkit-general
description: >-
  Generates 2D chemical structure depictions (PNG/SVG), interactive HTML
  compound grids, standalone chemical applications/visualizers, and ChemDraw
  (.cdxml/.cdx) documents using ACS 1996 publication styling; validates
  IUPAC/stereochemical nomenclature and provides shared RDKit version
  compatibility helpers.
metadata:
  publisher: google
---

# RDKit General Foundation (`rdkit-general`)

## Overview

`rdkit-general` provides 2D molecular depictions (`scripts/draw_molecules.py`),
sortable interactive HTML structure grids
(`scripts/interactive_grid_viewer.py`), ChemDraw `.cdxml`/`.cdx` document
conversion (`scripts/chemdraw_cdxml.py`), programmatic IUPAC and stereochemical
nomenclature verification (`scripts/verify_chemical_nomenclature.py`), shared
RDKit version compatibility and plotting helpers (`scripts/compat.py`), and
architectural patterns for building standalone chemical applications,
dashboards, and visualizers on top of the `cheminformatics-*` skill suite
(`references/rdkit_applications_and_visualizers.md`).

## Prerequisites

1.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH. Always execute skill scripts as `uv run
    /path/to/script.py ...` (or `python3 /path/to/script.py ...`); **never** run
    `uv run python3 /path/to/script.py`, as placing `python3` after `uv run`
    bypasses the script's PEP 723 `# /// script` dependency block.
2.  **User Notification**: If `.licenses/rdkit_general_LICENSE.txt` does not
    already exist in the workspace root then (1) prominently notify the user of
    the terms of use for the tools and data sources this skill relies on, then
    (2) create the file recording the notification text and timestamp:
    -   **OPSIN** (Apache 2.0, https://opsin.ch.cam.ac.uk) and **PubChem PUG
        REST** (https://pubchem.ncbi.nlm.nih.gov/docs/citation-guidelines), used
        by `scripts/verify_chemical_nomenclature.py`
    -   **RDKit**: BSD 3-Clause License (Landrum et al.)

## Core Rules

-   **Run via `uv run <script.py>`**: Execute every script directly as `uv run
    <path/to/script.py>` (without `python3` after `uv run`) so inline PEP 723
    metadata resolves dependencies (`rdkit`, `polite-http`, `pandas`, etc.).
-   **Write results to a file**: Pass `--output` (or `--output_csv` /
    `--output_json` / `--output_sdf`) so large result sets stay out of the
    context window, then read back only what you need.
-   **Never invent chemistry**: Do not report a descriptor, price, CAS number,
    binding affinity, or vendor that a script did not actually return. Cite the
    data source in the deliverable.
-   **Do not modify skill scripts in place, and never write scratch scripts for
    series harvesting or curation**: Run existing scripts in `scripts/`
    (`--help`) directly. For harvesting ChEMBL/PubMed literature, merging patent
    tables, aligning/coalescing assay columns, splitting cohorts, computing
    physicochemical descriptors, and plotting chemical space distributions, read
    the **`cheminformatics-core`** skill and use its
    `curate_congeneric_corpus.py` and `descriptors.py` workflows rather than
    writing custom Python scripts in `/tmp` or `scratch/`.
-   **No ASCII Structures, Monospace Tables, or Mermaid Chemical/SAR Diagrams**:
    Never render chemical structures, scaffolds, lead-optimization progressions
    (`Lead A -> Precursor B -> Drug C`), reaction schemes, or property tables
    using ASCII art, monospace code blocks, or `mermaid` flowcharts (`flowchart
    LR`/`TD` boxes). Always render 2D structure depictions with
    `scripts/draw_molecules.py` (or `render_mechanism.py` for reactions) and
    format all tables and compound comparisons using GitHub Flavored Markdown
    (GFM). Read `references/reporting_and_nomenclature_invariants.md` when
    formatting technical reports, GFM tables, `file://` links, or verifying
    IUPAC/stereochemical nomenclature.
-   **Publication Visual Styling (`scripts/compat.py`)**: Prefer built-in
    `--plot` flags over writing freeform plotting code. Whenever you write
    custom Matplotlib/Seaborn code, import the styling helpers from
    `scripts/compat.py`, call `set_safe_plot_style()` before creating figures
    (`layout="constrained"`), use `PUBLICATION_PALETTE` colors, and call
    `apply_publication_spines(ax)` on every subplot axis:

    ```python
    from compat import PUBLICATION_PALETTE, apply_publication_spines, set_safe_plot_style
    set_safe_plot_style()  # Call before plt.subplots(..., layout="constrained"); call apply_publication_spines(ax) on every axis
    ```

    Keep default opaque white backgrounds (`clearBackground = True`) or wrap 2D
    depictions in a white `.chem-svg-well` (`background: #ffffff; border-radius:
    8px; padding: 8px;`) so ACS 1996 black (`#000000`) bonds remain legible in
    both light and dark UIs; only pair `--transparent` with `--dark-mode` on
    `scripts/draw_molecules.py` when explicitly rendering light strokes onto a
    dark canvas. Use `scripts/interactive_grid_viewer.py` for interactive HTML
    compound grids ($\ge 5$ molecules).

--------------------------------------------------------------------------------

## 1. Empirical Consistency & Full-Dataset Integrity

-   **Full Dataset Utilization (No Arbitrary Downsampling or Limits)**: Never
    downsample, slice, or impose artificial row limits on a chemical series.
    Featurize, evaluate, and model 100% of the available dataset.
-   **Verify Before Writing**: Complete all batch computations and verify that
    output CSV/SDF/JSON files exist and match reported metrics ($r, \rho, R^2,
    \text{RMSE}, N$) verbatim before drafting report text.
-   **Primary-Key Joins & Partition Isolation**: Preserve primary keys
    (`compound_id`, `inchikey`, `canonical_smiles`) when merging tables, and
    keep `train_cv` and `held_out_test` strictly isolated.
-   **Scoped Missing-Data Handling**: Never call global `df.dropna()` across an
    entire dataset. Scope missing-value filtering strictly to the target column
    currently being modeled (`df.dropna(subset=[target])`), route unlabeled
    compounds to prospective scoring, and wrap per-molecule loops in
    `try/except`.
-   **Programmatic Nomenclature & Topological Verification
    (`scripts/verify_chemical_nomenclature.py`)**: Before citing systematic
    IUPAC names, macrocyclic ring sizes, or $(R)/(S)/(E)/(Z)$ designations,
    verify them programmatically via `uv run
    scripts/verify_chemical_nomenclature.py` (see
    `references/reporting_and_nomenclature_invariants.md`). If prochirality
    (`pro-R`/`pro-S`, `Re`/`Si`) or stereotopic/NMR equivalence (`homotopic` vs.
    `enantiotopic` vs. `diastereotopic`, which `Chem.CanonicalRankAtoms` does
    not distinguish even after `AssignStereochemistryFrom3D`) is important, use
    `uv run ../cheminformatics_reactions/scripts/prochirality.py --smiles
    "<SMILES>"` (`achiral_classes_by_element`).

--------------------------------------------------------------------------------

## 2. Modern RDKit API Invariants & `compat.py` Scope

Read `references/rdkit_api_and_sanitization.md` when using RDKit 2024.03+
fingerprint generators (`rdFingerprintGenerator`), `rdMolStandardize`,
`SanitizeMol` valence recovery, `UpdatePropertyCache(strict=False)` before ring
perception on unsanitized molecules, or descending index deletion in
`Chem.RWMol`.

-   **`scripts/compat.py` Scope Boundary**: `scripts/compat.py` is an RDKit
    runtime C++ version compatibility shim (`get_morgan_fingerprint`,
    `standardize_molecule`, `get_etkdg_params`, `safe_substruct_match`) and
    publication plot theming helper (`set_safe_plot_style`). **`compat.py` does
    NOT export domain descriptors or `calculate_cns_mpo`.**

--------------------------------------------------------------------------------

## 3. Canonical Descriptor Schema & Normalization

Read `references/descriptor_schema_and_canonicalization.md` when resolving RDKit
descriptor module locations (`Lipinski` vs `rdMolDescriptors` vs `Descriptors`),
normalizing dataframe column names, or writing custom descriptor calculation
code.

-   **Canonical Short Style (Default for Tables & Plots)**: `MW`, `cLogP`,
    `TPSA`, `HBD`, `HBA`, `RotB`, `Fsp3`, `HeavyAtoms`, `Rings`,
    `AromaticRings`, `AliphaticRings`, `QED`, `SA_Score`
-   **Canonical Snake Style (Alternative for Python Dicts)**:
    `molecular_weight`, `clogp`, `tpsa`, `hbd`, `hba`, `rotatable_bonds`,
    `fraction_csp3`, `heavy_atom_count`, `ring_count`, `aromatic_ring_count`,
    `aliphatic_ring_count`, `qed`, `sa_score`

--------------------------------------------------------------------------------

## 4. Module Imports & Visualization Executables

-   **Module Imports (`scripts/skill_paths.py`)**: When writing custom
    post-processing scripts, prepend the colon-separated directories printed by
    `scripts/skill_paths.py` to `PYTHONPATH` in the **same command** (each tool
    call runs in a fresh subshell, and `scripts/skill_paths.py` is relative to
    the `rdkit_general` skill directory—pass the full path to `skill_paths.py`
    when running from another working directory) and import shared helpers
    (`from compat import set_safe_plot_style`) rather than hardcoding absolute
    installation paths:

    ```bash
    PYTHONPATH="$(uv run scripts/skill_paths.py)${PYTHONPATH:+:$PYTHONPATH}" uv run scratch/apply_selection.py
    ```
-   **2D Structure Depictions (`scripts/draw_molecules.py`)**: Render 2D
    compound and grid depictions (`--input_csv`, `--ids`, `--legend_cols`; keep
    default opaque white backgrounds for dark-mode legibility, or pass
    `--dark-mode` with `--transparent` for light strokes on dark surfaces).
-   **Interactive 2D HTML Grid Viewer (`scripts/interactive_grid_viewer.py`)**:
    Generate sortable interactive HTML structure grids and property scatter
    plots for compound sets ($\ge 5$ molecules). **Critical UI Rendering Rule**:

    1.  First call `write_to_file` to create a placeholder
        `<artifact_dir>/interactive_compound_browser.html` with
        `ArtifactMetadata={"Summary": "...", "UserFacing": True,
        "RequestFeedback": False}` so the agent UI registers the artifact card
        in the side panel, then run `uv run scripts/interactive_grid_viewer.py
        --output <artifact_dir>/interactive_compound_browser.html` to populate
        it.
    2.  In all Markdown reports and chat responses, **always** link the HTML
        viewer via its `<artifact_dir>/` path
        (`[interactive_compound_browser.html](file://<artifact_dir>/interactive_compound_browser.html)`),
        **never** a repository source workspace path (source repository `.html`
        links open raw HTML source code instead of the interactive preview).
-   **ChemDraw `.cdxml`/`.cdx` Interoperability (`scripts/chemdraw_cdxml.py`)**:
    Convert SMILES/Molfiles to publication-styled `.cdxml` (`--style
    {acs1996,rsc,wiley}`), parse `.cdxml`/`.cdx` (`--format
    {smiles,sdf,json,rxn}`), or audit `.cdxml` files (`--inspect`). Read
    `references/chemdraw_cdxml_and_acs1996.md` when generating or parsing
    ChemDraw files.

--------------------------------------------------------------------------------

## 5. Building Standalone Chemical Applications & Visualizers

When asked to build a standalone chemical application, interactive molecule or
SAR visualizer, or web dashboard, **always build on top of `rdkit-general` and
the companion `cheminformatics-*` skills** rather than writing raw RDKit
boilerplate from scratch:

-   **Read `references/rdkit_applications_and_visualizers.md`**: Consult this
    reference for existing visualization tools across the skills, turnkey HTML
    viewers, `PYTHONPATH` module setup, responsive SVG generation, and
    Antigravity (AGY) host theme styling.
-   **Turnkey Zero-CDN HTML Visualizers**: Sandboxed artifact webviews block
    external JS CDNs (`3Dmol.js`, `SmilesDrawer`, `RDKit.js`, `Plotly`,
    `Chart.js`). Use `scripts/interactive_grid_viewer.py` for 2D compound
    browsers and SAR scatter dashboards.
-   **Import Verified Skill Modules via `skill_paths.py`**: In custom Python
    backends or multi-tab app generators, prepend the skill script directories
    printed by `scripts/skill_paths.py` (relative to the `rdkit_general` skill
    directory) to `PYTHONPATH` in the same command (`PYTHONPATH="$(uv run
    scripts/skill_paths.py)${PYTHONPATH:+:$PYTHONPATH}" uv run app.py`) and
    reuse `standardize_molecule`, `compute_standard_descriptors`,
    `calculate_cns_mpo`, and `screen_structural_alerts` instead of re-coding
    descriptor or alert logic.
-   **AGY Theme Variables, Structure Legibility & Plot Sizing**: Include
    `<script src="/preload.js"></script>` in `<head>` and style custom HTML UIs
    with host-injected AGY CSS variables (`var(--background)`, `var(--card)`,
    `var(--border)`, `var(--foreground)`, `var(--muted-foreground)`,
    `var(--primary)`, `var(--accent)`) without overriding `:root`. Always place
    standard ACS 1996 2D depictions inside a high-contrast white
    `.chem-svg-well` (`background: #ffffff`) or enable
    `rdMolDraw2D.SetDarkMode(opts)` (`--dark-mode`), and constrain all embedded
    viewer plots with `max-width: 100%; max-height: min(360px, 42vh);
    object-fit: contain;` so 300-DPI figures never overflow the viewport.
