# Congeneric Corpus Curation, Scaffold Triage, & Bioassay Unit Standards

This reference details end-to-end chemical series curation with
`scripts/curate_congeneric_corpus.py`, two-step patent and catalog scaffold
triage with `scripts/library_search.py`, and bioassay unit conversion rules.

--------------------------------------------------------------------------------

## 1. Two-Step Scaffold Triage for Patent & Vendor Harvests (`library_search.py`)

Patent and vendor harvests (e.g. from `patent_fetcher.py`) are **associative** —
they return prior art, comparators, synthetic intermediates, and reagents
alongside the target series. Always triage a patent or catalog harvest with
`library_search.py` before passing it to `curate_congeneric_corpus.py`:

### Step 1: Discover Scaffolds (`--scaffold_survey`)

Run `--scaffold_survey` first on any harvest to inspect the Murcko scaffolds
actually present, their member counts, and macrocyclic flags:

```bash
uv run scripts/library_search.py \
    --input patent_harvest.csv \
    --scaffold_survey
```

It prints each scaffold as a ready-to-paste `--query_scaffold` argument and
reveals when a patent covers multiple distinct series (e.g. an acyclic
predecessor series and its macrocyclized successor).

### Step 2: Filter to Series Cores (`--query_scaffold` / `--query_smarts`)

Pass one or more `--query_scaffold` flags (repeatable, OR-combined) to filter
the harvest to the series while preserving all original assay and provenance
columns (`Matched_Core` records which core each hit matched):

```bash
uv run scripts/library_search.py \
    --input patent_harvest.csv \
    --query_scaffold "<ACYCLIC_SCAFFOLD_SMILES>" \
    --query_scaffold "<MACROCYCLIC_SCAFFOLD_SMILES>" \
    --max_mw 900 --max_logp 8.0 --max_tpsa 250 --max_rotb 20 \
    --output_csv triaged_series.csv
```

### Critical Matching & Threshold Rules:

-   **Never paste a Murcko scaffold SMILES into `--query_smarts`**: Stripping an
    N-substituent leaves an `[nH]` whose explicit hydrogen count the substituted
    parent molecule fails, causing `--query_smarts` to match almost nothing.
    `--query_scaffold` compares Murcko scaffold to Murcko scaffold so the hit
    count matches `--scaffold_survey` exactly.
-   **When to use `--query_smarts`**: Use `--query_smarts '<core>'` (also
    repeatable and OR-combined) for a hand-written SMARTS core when the series
    varies at the ring framework itself (e.g. adding a fused ring).
-   **Why to avoid `--query_smiles` for series triage**: Tanimoto similarity
    anchored to a single lead excludes distant members of the same series (a
    macrocyclized analogue typically scores $< 0.2$ Tanimoto against its acyclic
    parent). Fall back to `--query_smiles <lead> --threshold 0.5` only when no
    set of cores describes the series.
-   **Raise Property Cutoffs for Macrocycles / Beyond-Ro5 Series**: Default
    cutoffs (`--max_mw 550`, `--max_logp 5.0`, `--max_tpsa 160`, `--max_rotb
    12`) reject macrocycles and larger analogues. Always raise them (`--max_mw
    900 --max_logp 8.0 --max_tpsa 250 --max_rotb 20`) when triaging a series
    rather than virtual-screening drug-like fragments.

--------------------------------------------------------------------------------

## 2. End-to-End Series Assembly (`curate_congeneric_corpus.py`)

`scripts/curate_congeneric_corpus.py` ingests ChEMBL documents (`--chembl-docs`
/ `--pmids`, or pre-harvested paper CSVs from
`../cheminformatics_doc2mol/scripts/paper_fetcher.py`) and local patent/vendor
CSVs (`--input-csv`) in a single call. Before deduplicating by
`canonical_smiles`, it pivots all distinct `(standard_type, target_name,
mutant)` bioassay endpoints into wide columns (`Ki_<target>`,
`IC50_<target>_<mutant>`, `Ratio_MDCK`, `LogD_...`, `CL_...`, plus
`all_endpoints` JSON), coalesces non-null metadata across colliding rows into
`all_sources`, computes the canonical 13-property descriptor set, Pfizer CNS
MPO, SA_Score, and InChIKeys, and outputs CSV, SDF, a GFM summary table, a
multi-panel distribution figure (`--output-plot`), an interactive HTML compound
browser (`--output-html`), and a 2D representative structures grid
(`--output-structures-png`).

```bash
uv run scripts/curate_congeneric_corpus.py \
    --chembl-docs CHEMBL<DOC_ID> \
    --input-csv triaged_series.csv \
    --group-by series \
    --output-csv curated_corpus.csv \
    --output-sdf curated_corpus.sdf \
    --output-html interactive_compound_grid.html \
    --output-structures-png representative_structures.png \
    --panels "scatter:MW:cLogP" "kde:TPSA" "box:QED" "hist:potency" \
    --output-plot chemical_space.png
```

> **Note on `--panels` and `--output-plot`**: `--panels` has no hardcoded
> default. If you pass `--output-plot` (or `--plot-only`) without `--panels`,
> `curate_congeneric_corpus.py` still writes `curated_corpus.csv`,
> `curated_corpus.sdf`, `--output-html`, and `--output-structures-png` and
> prints the cohort summary table first, then outputs a helper error listing all
> available physicochemical descriptors and discovered quantitative/discrete
> bioassay columns in the dataset so you can re-run with `--plot-only --panels
> ...`.

### Column Coalescing (`--col-map` / `--rename-cols`), Cohort Splitting (`--split-by`), & Configurable Panels (`--panels` / `--plot-only`)

-   **Aligning & Coalescing Heterogeneous Assay Columns (`--col-map`)**: When
    merging ChEMBL documents (`Ki_<target>`) with patent tables
    (`patent_assay_<header>_nM` or `<PATENT_ID>_col_1_nM`), use `--col-map` (or
    `--rename-cols`) to coalesce source-specific columns into a single unified
    endpoint column (`NewCol=old_col1|old_col2` or `file.csv:old=new`). You can
    pass `--col-map` upfront during curation or re-run with `--plot-only` on an
    existing CSV to coalesce columns and export `--split-by` CSV/SDF subsets in
    place without re-fetching ChEMBL:

    ```bash
    uv run scripts/curate_congeneric_corpus.py \
        --plot-only --input-csv curated_corpus.csv --output-csv curated_corpus.csv \
        --output-sdf curated_corpus.sdf \
        --col-map "Primary_Ki_nM=Ki_Target_Name|patent_assay_target_ki_nM" \
        --group-by series --split-by series \
        --panels "scatter:MW:TPSA" "box:SA_Score" "kde:cLogP" "hist:potency" \
        --output-plot custom_panels.png
    ```

-   **Keep Biochemical/Enzymatic and Cellular Endpoints Separate in
    `--col-map`**: Never coalesce cell-free biochemical/enzymatic
    $\text{IC}_{50}$ columns (e.g. `patent_assay_*_enzymatic_assay_*_nM` or
    biochemical ChEMBL `IC50_<Target_Protein>`) into a cellular
    proliferation/ELISA column (e.g. `Cellular_<Target>_IC50_nM`), or label an
    enzymatic patent $\text{IC}_{50}$ column as a cellular endpoint. Map
    enzymatic $K_i$, enzymatic $\text{IC}_{50}$, and cellular $\text{IC}_{50}$
    into distinct columns so intrinsic target affinity is never confounded with
    cellular permeability or efflux.

-   **Include Heteroaryl Bioisosteres in `--series-smarts`**: When partitioning
    sub-series via `--series-smarts` (`Label=SMARTS`), use `[c,n]` rather than
    carbon-only `c` at peripheral aryl positions where a phenyl ring may be
    replaced by a pyridine or pyrazine bioisostere, so heteroaromatic analogues
    from the same generation do not fall through into an earlier fallback
    cohort.

-   **Never Write Custom Post-Processing Scripts**: Do not write custom scratch
    scripts to rename/coalesce endpoint columns, split CSV/SDF cohorts, or
    redraw property panels — `curate_congeneric_corpus.py` handles all three
    natively via `--col-map`, `--split-by`, and `--panels`.

--------------------------------------------------------------------------------

## 3. Bioassay Units & PubChem Ingestion Traps

-   **Unit Assumption Trap**: Never assume raw PubChem bioassay numbers are
    nanomolar ($\text{nM}$) and blindly compute $pK_i = 9.0 - \log_{10}(X)$.
-   **Standard PubChem Bioassay Units**: In PubChem bioassay CSVs (e.g. AID
    1153093 / 1153094), activity values are typically reported in **Micromolar
    ($\mu\text{M}$)** ($10^{-6}\text{ M}$): $$pK_i = 6.0 -
    \log_{10}(\text{Value}_{\mu\text{M}}), \qquad pK_i = 9.0 -
    \log_{10}(\text{Value}_{\text{nM}})$$ Always inspect the `standard_units`
    column (`nM` vs `uM` / `µM`) before converting to logarithmic potency
    (`pActivity` / `pChEMBL`).
