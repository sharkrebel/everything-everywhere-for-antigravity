---
name: cheminformatics-core
description: >-
  Computes RDKit physicochemical descriptors (MW, cLogP, TPSA, CNS MPO, SA
  score), PAINS/Brenk structural alerts, pH 7.4 ionization, congeneric series
  curation, activity cliffs, combinatorial library enumeration, vendor pricing,
  and organometallic %V_bur sterics.
metadata:
  publisher: google
---

# RDKit Cheminformatics Hub (`cheminformatics-core`)

## Prerequisites

1.  **`rdkit-general`**: Read the `rdkit-general` skill before running any
    workflow here.
2.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH. Always execute skill scripts as `uv run
    /path/to/script.py ...` (or `python3 /path/to/script.py ...`); **never** run
    `uv run python3 /path/to/script.py`, as placing `python3` after `uv run`
    bypasses the script's PEP 723 `# /// script` dependency block.
3.  **User Notification**: If `.licenses/cheminformatics_core_LICENSE.txt` does
    not already exist in the workspace root then (1) prominently notify the user
    of the terms of use for the tools and data sources this skill relies on,
    then (2) create the file recording the notification text and timestamp:
    -   **RDKit**: BSD 3-Clause License (https://www.rdkit.org/)
    -   **PubChem PUG REST**
        (https://pubchem.ncbi.nlm.nih.gov/docs/citation-guidelines), used by
        `bioassay_data_curator.py`, `find_coupling_partners.py`,
        `enumerate_matrix_library.py`, `chemical_vendor_resolver.py`, and
        `validate_chemical_table.py`
    -   **ChEMBL / EMBL-EBI** (CC BY-SA 3.0,
        https://chembl.gitbook.io/chembl-interface-documentation), used by
        `curate_congeneric_corpus.py` and `bioassay_data_curator.py`
    -   **NIST Chemistry WebBook** (https://webbook.nist.gov/chemistry/), used
        by `nist_webbook_fetcher.py`
    -   **Enamine** (https://enamine.net), **MolPort**
        (https://www.molport.com), and **Chemspace** (https://chem-space.com)
        vendor catalogues, used by `screen_compound_pricer.py` and
        `chemical_vendor_resolver.py`
4.  **Optional Live Vendor Pricing Credentials (`credentials`)**: If querying
    live MolPort or Chemspace commercial pricing via `screen_compound_pricer.py`
    (rather than offline catalog CSVs), read the `credentials` skill and check
    for `MOLPORT_API_KEY` / `CHEMSPACE_API_KEY` in `~/.env` without printing
    secret values (`grep -sq '^MOLPORT_API_KEY=.\+' ~/.env` and `grep -sq
    '^CHEMSPACE_API_KEY=.\+' ~/.env`).

## Core Rules

-   **Shared foundation rules**: Follow all core rules in
    `../rdkit_general/SKILL.md`.
-   **No ASCII Structures or Monospace Tables**: Format all tables using GitHub
    Flavored Markdown (GFM) and render 2D depictions with
    `../rdkit_general/scripts/draw_molecules.py` (see
    `../rdkit_general/references/reporting_and_nomenclature_invariants.md`).
-   **Large Datasets**: For out-of-core SDF streaming, two-tier substructure
    prescreening, or bulk similarity, follow
    `references/performance_and_optimizations.md`.

--------------------------------------------------------------------------------

## Task-to-Domain Routing Directory

Match the user's objective to the domain below, **read the corresponding
reference guide**, and execute the designated script from `scripts/`:

1.  **Molecular Standardization, Descriptors & Corpus Curation**

    -   **Reference Guides**:
        -   `references/corpus_curation_and_scaffold_triage.md` — read when
            curating chemical series (`curate_congeneric_corpus.py`), triaging
            patent harvests with `--scaffold_survey` / `--query_scaffold`,
            configuring `--panels` / `--plot-only`, or reconciling PubChem
            $\mu\text{M}$ vs $\text{nM}$ units.
        -   `references/api_reference.md` — read when computing 2D RDKit
            descriptors, Pfizer CNS MPO (`calculate_cns_mpo`), or molecular
            standardization and pH 7.4 ionization states.
    -   **End-to-End Series Pipeline**:
        -   When curating a patent or catalog harvest, **always triage first**
            with `uv run scripts/library_search.py --scaffold_survey`, then
            filter to the target series cores with `--query_scaffold
            '<scaffold>'` (raising `--max_mw 900 --max_logp 8.0 --max_tpsa 250
            --max_rotb 20` so macrocycles and larger analogues are retained;
            never paste a Murcko SMILES into `--query_smarts`).
        -   Run `uv run scripts/curate_congeneric_corpus.py` (`--input-csv`,
            `--chembl-docs` / `--pmids`, `--col-map` / `--rename-cols`,
            `--group-by`, `--split-by`, `--panels`, `--highlight`,
            `--plot-only`) to merge literature and triaged patent tables,
            coalesce heterogeneous endpoint columns via `--col-map` (e.g.
            `--col-map "Unified_Ki_nM=colA|colB"`), deduplicate by
            `canonical_smiles` while merging provenance (`all_sources`,
            `patent_example_id`), compute the 13 canonical descriptors + CNS MPO
            +   SA_Score, and emit master + per-cohort (`--split-by`) CSVs/SDFs,
                a GFM summary table, and a multi-panel chemical space plot.
    -   **Single-Purpose Executables**:
        -   `uv run scripts/descriptors.py` — canonical 13-property featurizer,
            column normalizer (`DescriptorCanonicalizer`), and Pfizer CNS MPO
            scorer (`calculate_cns_mpo` lives here, never in `compat.py`).
        -   `uv run scripts/filter_catalogs.py` — PAINS/Brenk/NIH/ZINC
            structural alerts and CNS MPO screening.
        -   `uv run scripts/standardize_molecule.py` — salt stripping, charge
            neutralization, and tautomer canonicalization.
        -   `uv run scripts/standardize_physiological_ionization.py` — applies
            predominant pH 7.4 protonation microstates across a CSV.
        -   `uv run scripts/peptide_helm_converter.py` — converts peptide
            sequences, FASTA, and HELM notation to molecules.
        -   `uv run scripts/synthetic_accessibility_scorer.py` — Ertl &
            Schuffenhauer SA_Score (1 easy to 10 hard).
        -   `uv run scripts/compute_descriptors.py` — single-SMILES 2D
            properties, complexity indices, and alert flags.
        -   `uv run scripts/bioassay_data_curator.py` — fetches, normalizes, and
            caches ChEMBL/PubChem bioassay SAR datasets and activity cliffs
            (`--pmids`, `--chembl-docs`, `--pubchem_aid`, `--input_csv`).
        -   `uv run scripts/nist_webbook_fetcher.py` — NIST WebBook experimental
            thermochemistry ($\Delta_f H^\circ$, $T_{\text{boil}}$,
            $T_{\text{fus}}$, Antoine constants, 2D MOL).
        -   `uv run scripts/pareto_frontier_optimizer.py` — non-dominated
            multi-objective lead selection.
        -   `uv run ../rdkit_general/scripts/draw_molecules.py` — 2D depictions
            (PNG/SVG) for single compounds and grids (`--input_csv`, `--ids`,
            `--legend_cols`).

2.  **DataFrames, Commercial Pricing, Ecosystem Integration & Performance**

    -   **Reference Guides**:
        -   `references/chemical_pricing_and_procurement.md` — read when
            querying commercial vendor pricing, pack tiers, or normalized price
            per milligram (`$/mg`).
        -   `references/chemical_supplier_apis.md` — read when resolving CAS
            registry numbers or looking up commercial building block suppliers
            (`resolve-cas`, `lookup-suppliers`).
        -   `references/pandas_and_dataframe_integration.md` — read when
            validating chemical DataFrames, loading/writing SDF files, or
            rendering interactive HTML grids.
        -   `references/ecosystem_integrations.md` — read when converting
            chemical file formats with Open Babel or integrating external
            cheminformatics toolkits.
        -   `references/performance_and_optimizations.md` — read when
            parallelizing descriptor computation or screening large
            multi-thousand-compound libraries.
    -   **Standard Executables**:
        -   `uv run scripts/screen_compound_pricer.py` — supplier pricing and
            pack tiers across MolPort, Chemspace, Mcule, Enamine, and
            Sigma-Aldrich, normalized to price per milligram.
        -   `uv run scripts/chemical_vendor_resolver.py` — CAS registry
            resolution and supplier lookup (`resolve-cas`, `lookup-suppliers`).
        -   `uv run ../rdkit_general/scripts/interactive_grid_viewer.py` —
            interactive HTML compound grid with sortable properties.
        -   `uv run scripts/validate_chemical_table.py` — parallel integrity
            audit of SMILES, CAS, property, and vendor columns.
        -   `uv run scripts/format_converter_openbabel.py` — format conversion
            via Open Babel with native RDKit fallback.

3.  **Chemical Reactions, SMARTS, Library Space & Coupling Partner Discovery**

    -   **Reference Guides**:
        -   `references/chemical_reactions_and_smarts.md` — read when writing
            Reaction SMARTS, performing substructure or MCS searches, or running
            one-step retrosynthetic disconnections.
        -   `references/library_and_reaction_space_search.md` — read when
            discovering coupling partners, running $M \times N$ matrix library
            enumerations, or virtual screening.
    -   **Standard Executables**:
        -   `uv run scripts/chemical_reaction.py` — builds, validates, and
            applies Reaction SMARTS.
        -   `uv run scripts/substructure_search.py` — SMARTS matching and MCS
            alignment with 2D highlighting.
        -   `uv run scripts/find_coupling_partners.py` — streams commercial
            coupling partners from PubChem, runs the in-silico reaction, and
            profiles the products.
        -   `uv run scripts/library_search.py` — scaffold survey
            (`--scaffold_survey`), scaffold filtering (`--query_scaffold`), and
            virtual screening.
        -   `uv run scripts/enumerate_matrix_library.py` — $M \times N$
            combinatorial enumeration with property profiling.
        -   `uv run scripts/retrosynthesis_decomposer.py` — one-step
            disconnections, RECAP tree, and BRICS synthon sourcing.

4.  **Organometallics, Ligand Sterics, Coordination Stereochemistry & Catalytic
    Cycles**

    -   **Reference Guides**:
        -   `references/coordination_chemistry_and_dative_bonds.md` — read when
            standardizing dative bonds (`->`), counting CBC MLX valence
            electrons, or rendering 2D `horizon_perspective` complexes.
        -   `references/ligand_descriptors_and_buried_volume.md` — read when
            computing SambVca buried volume ($\%V_{\text{bur}}$), Tolman cone
            angles, or 2D steric topographic maps.
    -   **Standard Executables**:
        -   `uv run scripts/organometallic_standardizer.py` — dative bond
            normalization, oxidation state, and 18e counting (`--render
            --render_style horizon_perspective`).
        -   `uv run scripts/ligand_steric_descriptors.py` — buried volume
            $\%V_{\text{bur}}$ and Tolman cone angle.
        -   `uv run scripts/catalytic_reaction_engine.py` — curated Reaction
            SMARTS for cross-couplings, C-H activation, and metathesis.
        -   `uv run scripts/coordination_complex_builder.py` — assembles and
            minimizes 3D transition metal coordination polyhedra.
        -   `uv run scripts/render_steric_map.py` — 2D topographic steric
            contour maps around the metal coordination sphere.

5.  **Patent Literature Retrieval & Document Series Curation
    (`cheminformatics-doc2mol`)**

    -   When extracting chemical series or assay tables from patents or papers,
        read `../cheminformatics_doc2mol/SKILL.md` via `view_file` and execute
        `uv run ../cheminformatics_doc2mol/scripts/patent_fetcher.py` and `uv
        run ../cheminformatics_doc2mol/scripts/document_mol_curator.py` (never
        write custom `requests`/`BeautifulSoup` Google Patents scrapers).

--------------------------------------------------------------------------------

## Progressive Execution Protocol

1.  **Align on Scope & Procurement Quantity**: When pricing compounds with
    `screen_compound_pricer.py`, if the user has not specified a purchase
    quantity, ask how much material they need per compound (e.g. 1–2 mg vs.
    20–50 mg) before passing `--target-amount-mg <amount>`.
2.  **Inspect Domain Reference**: Call `view_file` on the target reference guide
    above to load exact CLI flags and edge cases.
3.  **Execute or Adapt via Scratch**: Run `uv run scripts/...`. If custom logic
    is needed, copy the template script to `scratch/` first (`cp scripts/...
    scratch/...`).
4.  **Summarize Findings & Structure Depictions**: When producing a summary
    report, format tables in GitHub Flavored Markdown (GFM) with 2D structure
    depictions (`../rdkit_general/scripts/draw_molecules.py`). For compound sets
    ($\ge 5$ molecules), register a placeholder
    `<artifact_dir>/interactive_compound_browser.html` via `write_to_file`
    (`UserFacing: True`), populate it with
    `../rdkit_general/scripts/interactive_grid_viewer.py --output
    <artifact_dir>/interactive_compound_browser.html`, and always link
    `[interactive_compound_browser.html](file://<artifact_dir>/interactive_compound_browser.html)`
    from the `<artifact_dir>/` artifact directory (never link a source code
    workspace path for `.html` files).
