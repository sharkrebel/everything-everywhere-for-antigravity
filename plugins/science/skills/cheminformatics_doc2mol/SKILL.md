---
name: cheminformatics-doc2mol
description: >-
  Extracts and curates validated chemical structures from patents, literature,
  synthesis procedures, and Markush core/R-group tables: resolves IUPAC names
  via OPSIN and PubChem, retrieves Google Patents full text and assay tables,
  and exports curated SDF/CSV.
metadata:
  publisher: google
---

# RDKit Chemical Document Intelligence Hub (`cheminformatics-doc2mol`)

## Prerequisites

1.  **`rdkit-general`**: Read the `rdkit-general` skill before running any
    workflow here.
2.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH. Always execute skill scripts as `uv run
    /path/to/script.py ...` (or `python3 /path/to/script.py ...`); **never** run
    `uv run python3 /path/to/script.py`, as placing `python3` after `uv run`
    bypasses the script's PEP 723 `# /// script` dependency block
    (`polite-http`, `rdkit`, etc.).
3.  **User Notification**: If `.licenses/cheminformatics_doc2mol_LICENSE.txt`
    does not already exist in the workspace root then (1) prominently notify the
    user of the terms of use for the tools and data sources this skill relies
    on, then (2) create the file recording the notification text and timestamp:
    -   **RDKit**: BSD 3-Clause License (https://www.rdkit.org/)
    -   **OPSIN** (Apache 2.0, https://opsin.ch.cam.ac.uk) and **NCI CACTUS**
        (https://cactus.nci.nih.gov/chemical/structure), used by
        `iupac_resolver.py` and `patent_fetcher.py`
    -   **PubChem PUG REST** and **NCBI E-utilities**
        (https://www.ncbi.nlm.nih.gov/home/about/policies/), used by
        `iupac_resolver.py`, `patent_fetcher.py`, `paper_fetcher.py`, and
        `precedent_lookup.py`
    -   **Google Patents** (https://patents.google.com) and **ChEMBL /
        EMBL-EBI** (CC BY-SA 3.0,
        https://chembl.gitbook.io/chembl-interface-documentation), used by
        `patent_fetcher.py`, `paper_fetcher.py`, and `precedent_lookup.py`
    -   **Open Reaction Database (ORD)** (CC BY-SA 4.0,
        https://open-reaction-database.org), used by `precedent_lookup.py`

## Core Rules

-   **Shared foundation rules**: Follow all core rules in
    `../rdkit_general/SKILL.md`.
-   **Exhaustive Production Patent Retrieval**: In production curation
    workflows, leave `--limit` unset (`None`, the default) when running
    `scripts/patent_fetcher.py` (the `--limit` flag is intended only for quick
    connectivity or smoke testing). Patent databases index synthetic
    intermediates before exemplified active leads; letting multi-page API
    pagination and OPSIN resolution run to completion ensures no active examples
    are missed.
-   **Predecessor & Lead Composition-of-Matter Patent Families**: When curating
    a medicinal chemistry series to study structural transitions (e.g. scaffold
    hops, ring closures, or bioisosteres), ingest the primary
    composition-of-matter patent families (`--include-family`, $150 \le
    \text{MW} \le 900\text{ Da}$) for **both** the structural predecessor series
    and the optimized lead series alongside target-wide ChEMBL bioactivity
    records (`target_chembl_id`). Read
    `references/synthesis_and_patent_retrieval.md` for full patent family, HTML
    table, and grouped chiral peak rules.
-   **Full-Text Google Patents HTML Tables & Grouped Chiral Separations**:
    Always run `scripts/patent_fetcher.py` (never write custom
    `requests`/`BeautifulSoup` Google Patents scrapers). It automatically parses
    full-text Google Patents HTML `<table>` potency rows (`Example N`),
    propagates shared parent synthesis headings across grouped chiral SFC peaks
    (`Peak 1` / `Peak 2`), repairs OCR locant typos
    (`repair_chemical_ocr_text()`), and resolves IUPAC names via OPSIN.
-   **Comprehensive Multi-Parameter Bioactivity & ADMET Extraction**: Capture
    all available secondary bioassays, mutant panels, cellular assays, and
    ADMET/PK endpoints ($\text{CL}_{\text{int}}$, $P_{\text{app}}$, efflux
    ratios, $\log D$, off-target selectivity) in addition to primary $K_i$ /
    $\text{IC}_{50}$.
-   **Scaffold Triage & Two-Tier Patent Role Retention**: Because patent
    harvests are associative (including reagents and prior art), always triage
    `patent_fetcher.py` outputs with
    `../cheminformatics_core/scripts/library_search.py` (`--scaffold_survey`,
    then `--query_scaffold` with `--max_mw 900 --max_logp 8.0 --max_tpsa 250
    --max_rotb 20`) before passing to `curate_congeneric_corpus.py`. Retain both
    quantitatively labeled compounds (`dataset_role='labeled_train_eval'`) and
    unlabeled exemplified patent leads (`dataset_role='prospective_patent'`).

--------------------------------------------------------------------------------

## Task-to-Domain Routing Directory

Match the user's objective to the domain below, **read the corresponding
reference guide**, and execute from `scripts/`:

1.  **Chemical Nomenclature & Text/OCR Sanitization**

    -   **Reference Guide**: `references/nomenclature_and_text.md` — read when
        sanitizing OCR/PDF chemical text, repairing broken SMILES, or resolving
        systematic IUPAC names via OPSIN/PubChem/CACTUS.
    -   **Standard Executables**:
        -   `uv run scripts/chemical_text_sanitizer.py` — repairs typography,
            OCR, and line-break/zero-to-oxygen SMILES corruptions.
        -   `uv run scripts/iupac_resolver.py` — resolves IUPAC names, trade
            names, development codes, and CAS numbers via OPSIN, PubChem, and
            CACTUS in cascade.

2.  **Patent Tables & Markush Enumeration**

    -   **Reference Guide**: `references/markush_and_patent_tables.md` — read
        when assembling Markush core scaffolds (`[*:1]`, `[*:2]`) with R-group
        substituent tables.
    -   **Standard Executables**:
        -   `uv run scripts/markush_assembler.py` — combines a generic core with
            an R-group substituent table into enumerated, validated structures.

3.  **Synthesis Procedures, Academic Paper Harvesting & Public Patent
    Retrieval**

    -   **Reference Guide**: `references/synthesis_and_patent_retrieval.md` —
        read when harvesting patents (`patent_fetcher.py`), harvesting academic
        papers (`paper_fetcher.py`), looking up ORD/USPTO or patent reaction
        precedents, or extracting experimental procedures.
    -   **Standard Executables**:
        -   `uv run scripts/patent_fetcher.py` — retrieves patent structures
            from ChEMBL/PubChem, expands patent families (`--include-family`),
            and parses Google Patents full-text HTML biological assay tables
            (`Example N`) via OPSIN into one unified CSV.
        -   `uv run scripts/paper_fetcher.py` — retrieves academic literature
            structures, `paper_compound_id` labels (`1a`, `2b`), and pivoted
            WT/mutant/ADMET bioassay columns from ChEMBL and PubMed/PubChem
            (`--pmids`, `--chembl-docs`, `--dois`, `--output-csv`).
        -   `uv run scripts/precedent_lookup.py` — searches ORD/USPTO reaction
            corpora via the pre-indexed ORD v2 API + vectorized DRFP similarity
            and patent literature (`--api all`) via batched PubChem + Google
            Patents XHR/HTML, clustering condition/yield distributions across
            single reactions or multi-step routes (`--batch-json`).
        -   `uv run scripts/procedure_extractor.py` — parses experimental
            synthesis paragraphs into products, salts, yields, and LC-MS parent
            ions.

4.  **SAR Table Extraction & Series Curation**

    -   **Reference Guide**: `references/sar_tables_and_curation.md` — read when
        normalizing bioassay SAR tables (`<`, `>`, `=`), converting molar
        potency to pActivity, or auditing series completeness.
    -   **Standard Executables**:
        -   `uv run scripts/sar_table_extractor.py` — parses assay tables,
            normalizing inequality operators (`<`, `>`, `=`) and converting
            molarity to pActivity.
        -   `uv run scripts/document_mol_curator.py` — canonicalizes structures,
            deduplicates via InChIKey, and exports CSV/SDF/HTML.

--------------------------------------------------------------------------------

## Progressive Execution Protocol

1.  **Inspect Domain Reference**: Call `view_file` on the target reference guide
    in `references/` for detailed flags and edge cases.
2.  **Execute Modular Script (Do Not Rush or Cap)**: Run the `scripts/` utility
    with `--output-csv` and `--output-sdf` without `--limit` caps, allowing
    multi-page API/OPSIN resolution to complete.
3.  **Triage & Curate**: Triage patent harvests with `library_search.py`
    (`--scaffold_survey` $\to$ `--query_scaffold`), then run
    `document_mol_curator.py` or `curate_congeneric_corpus.py`.
4.  **Summarize Findings & Depictions**: Write extraction yields, provenance
    breakdowns, and GFM tables to a markdown summary in the output directory,
    generate 2D structure depictions
    (`../rdkit_general/scripts/draw_molecules.py`), and for compound sets ($\ge
    5$ molecules) generate an interactive HTML viewer via
    `../rdkit_general/scripts/interactive_grid_viewer.py` (or pass
    `--output-html` and `--output-structures-png` to
    `../cheminformatics_core/scripts/curate_congeneric_corpus.py`).
