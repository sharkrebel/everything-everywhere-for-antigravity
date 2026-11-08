# 03. Synthesis Procedures, Patent Family Retrieval, & Full-Text Curation Guide

## Overview

Extraction of chemical synthesis procedures from experimental paragraphs,
characterization data verification, zero-auth patent family & full-text HTML
table retrieval (`patent_fetcher.py`), and multi-parameter bioactivity curation.

--------------------------------------------------------------------------------

## 1. Patent Family Scope & Exhaustive Retrieval (`patent_fetcher.py`)

-   **No-Limit Exhaustive Retrieval**: Never pass `--limit` values to
    `scripts/patent_fetcher.py`. Patent databases index synthetic intermediates
    before exemplified active products, and large pharmaceutical patents
    frequently contain thousands of structures. Leave `--limit` unset
    (defaulting to `None`, no limit) and allow multi-page API pagination and
    nomenclature resolution passes to finish so late-stage exemplified leads are
    never truncated.
-   **Predecessor & Lead Composition-of-Matter Patent Families
    (`--include-family`, enabled by default)**: When curating a medicinal
    chemistry series to analyze structural transitions (such as scaffold hops,
    ring closures/expansions, or bioisosteric replacements), ingest primary
    composition-of-matter patents for **both** the structural predecessor series
    and the optimized lead series (`150 <= MW <= 900 Da`), alongside target-wide
    ChEMBL bioactivity records (`target_chembl_id`). By default
    (`--include-family`), `patent_fetcher.py` parses Google Patents metadata for
    every seed patent to discover and harvest compounds from:
    1.  **Continuations & Divisionals**
        (`itemprop="continuation|childApplications"`).
    2.  **DOCDB & PCT (`WO`) Family Counterparts** (`itemprop="docdbFamily"`).
    3.  **Foundational Parent Patents** (`backwardReferences` by the same
        assignee lineage).
    4.  **Follow-On & Bridge Series Patents** (`forwardReferences` by the same
        assignee lineage).

--------------------------------------------------------------------------------

## 2. Full-Text Google Patents HTML Tables & Grouped Chiral Separations

-   **Automated HTML Biological Table & Example Parsing**: Public databases
    (PubChem/SureChEMBL) index patent structures but frequently omit tabular
    potency data (`TABLE 1`, `Example N` $K_i$/$\text{IC}_{50}$ tables) printed
    in the specification. `patent_fetcher.py` automatically downloads the
    full-text Google Patents HTML
    (`https://patents.google.com/patent/<PATENT_ID>/en`), scans every `<table>`
    for example rows (`12`, `12a`, `Example 12`), resolves `Example N` synthesis
    IUPAC headings via parallel OPSIN queries, verifies `LCMS [M+H]+` parent
    ions within $1.5\text{ Da}$, and attaches tabular potencies
    (`patent_table_col_1`, `pchembl_value`) directly to the output CSV. Never
    write ad-hoc `requests`/`BeautifulSoup` Google Patents scrapers.
-   **Grouped Chiral Separations (`Peak 1` / `Peak 2`) & OCR Repair**: Patent
    synthesis sections often group chiral SFC separations, enantiomers, or
    diastereomers under a shared parent heading (e.g. `Example A and Example B`
    followed by `Peak 1` and `Peak 2`), placing the shared IUPAC name several
    paragraphs above the individual `Example N` tags. `patent_fetcher.py`
    propagates shared parent synthesis headings across paired sub-examples and
    sanitizes numeric locant OCR typos (such as letter `O` for digit `0`) via
    `repair_chemical_ocr_text()` before calling OPSIN.

--------------------------------------------------------------------------------

## 3. Multi-Parameter Bioactivity & Two-Tier Patent Role Curation

-   **Comprehensive Multi-Parameter Bioactivity & ADMET Extraction**: When
    extracting data from ChEMBL documents or literature tables, do not filter
    solely for primary target potency ($K_i$ or $\text{IC}_{50}$). Capture all
    available secondary bioassays, mutant panels, cellular assays, and ADMET/PK
    endpoints (intrinsic microsomal clearance $\text{CL}_{\text{int}}$, passive
    permeability $P_{\text{app}}$, efflux ratios, $\log D$, and off-target
    selectivity) so multi-parameter profiles and selectivity margins can be
    evaluated completely.
-   **Two-Tier Patent Role Separation & Scaffold Triage**:
    1.  First triage associative patent harvests with
        `../cheminformatics_core/scripts/library_search.py`
        (`--scaffold_survey`, then `--query_scaffold` with `--max_mw 900
        --max_logp 8.0 --max_tpsa 250 --max_rotb 20`) to strip reagents,
        intermediates, and unrelated prior art.
    2.  Retain **both** quantitatively labeled structures
        (`dataset_role='labeled_train_eval'`, with measured
        $K_i$/$\text{IC}_{50}$ for supervised QSAR training) and unlabeled
        exemplified patent leads (`dataset_role='prospective_patent'`) so
        unlabeled patent leads can be scored prospectively across downstream
        physics and QSAR models.

--------------------------------------------------------------------------------

## 4. Reaction Precedent Lookup (`precedent_lookup.py`)

-   **Concurrent Pre-Indexed ORD v2 + Patent Precedent Lookup**: When
    `--api all` (or `--api ord --include-patents`) is specified,
    `precedent_lookup.py` queries the live pre-indexed ORD v2 search API
    (`POST /api/search` + batch `GET /api/reactions`) concurrently alongside
    batched PubChem InChIKey resolution, Google Patents full-text synthesis XHR
    search, and selective Robertson-Sparck Jones empirical IDF `PatentID`
    co-occurrence on the highest-Bertz-complexity organic structures. Passing
    `--include-patents` without `--api` (or `--api patents`) queries only
    PubChem and Google Patents.
-   **Vectorized Batch DRFP, Target-Product Synthesis Lookup (`--product`), &
    Multi-Step `--batch-json`**: Candidate reactions are normalized to unmapped
    canonical reaction SMILES, scored in a single vectorized
    `DrfpEncoder.encode` pass (or product Morgan Tanimoto similarity when
    querying syntheses of a target molecule via `--product "<TARGET_SMILES>"` or
    `--reaction-smiles ">><TARGET_SMILES>"` without reactants), and parsed
    concurrently via `ProcedureExtractor`. For multi-step synthesis routes, pass
    `--batch-json <steps.json>` to evaluate all steps concurrently in one
    process, amortizing Python/RDKit startup and sharing in-memory PubChem
    compound/patent cross-reference and Google Patents HTML caches across steps.
    The JSON file may be either a non-empty list `[...]` of step objects or
    `{"steps": [...]}`; each step object requires non-empty `reaction_smiles`,
    `reactants`, or `product`, preserves route metadata (`proposal_id`,
    `step_id`, `step_number`, `title`), and accepts optional per-step overrides
    (`api`, `include_patents`, `keywords`, `min_similarity`, `limit`,
    `max_patents`, `precedent_db`; `null` values inherit the batch defaults).

--------------------------------------------------------------------------------

## 5. CLI Examples

```bash
# Retrieve seed patent structures + automatically discover and harvest linked family, parent, and follow-on patents
uv run scripts/patent_fetcher.py \
  --patent-id "SEED_PATENT_1,SEED_PATENT_2" \
  --include-family \
  --max-related-patents 3 \
  --output-csv /tmp/patent_mols.csv \
  --output-sdf /tmp/patent_mols.sdf

# Extract synthesis procedure, salts, and MS characterization from experimental text
uv run scripts/procedure_extractor.py \
  --input-text "Example 12: Preparation of 4-(4-chlorophenyl)thiazol-2-amine... MS (ESI+): m/z 211.0 [M+H]+." \
  --output-json /tmp/extracted.json --output-csv /tmp/extracted.csv

# Look up existing syntheses of a target molecule (--product)
uv run scripts/precedent_lookup.py \
  --product "COc1ccc(-c2ccccc2)cc1" \
  --api all \
  --output-json /tmp/target_syntheses.json --output-markdown /tmp/target_syntheses.md

# Look up reaction precedents across ORD v2 and Google Patents (--api all)
uv run scripts/precedent_lookup.py \
  --reaction-smiles "CC(=O)O.CCN>>CC(=O)NCC" \
  --api all \
  --keywords "amidation" \
  --output-json /tmp/precedents.json --output-markdown /tmp/precedents.md

# Evaluate a multi-step synthesis route in a single process (--batch-json)
uv run scripts/precedent_lookup.py \
  --batch-json /tmp/route_steps.json \
  --api all \
  --limit 25 \
  --output-json /tmp/route_precedents.json --output-markdown /tmp/route_precedents.md
```
