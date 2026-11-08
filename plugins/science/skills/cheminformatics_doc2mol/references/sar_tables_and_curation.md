# 04. SAR Bioactivity Tables & Curation Guide

## Overview

Extraction and standardization of SAR bioactivity tables from literature/patents
and comprehensive chemical library curation.

## Key Capabilities & CLI Recipes

-   **SAR Bioactivity Extraction**: Parses inequality operators (`<`, `>`, `<=`,
    `>=`, `=`), separates values from units (`pM`, `nM`, `uM`, `mM`), normalizes
    concentrations to molarity (M), and calculates logarithmic potencies
    (`pIC50`, `pKi`, `pActivity`).
-   **Library Curation & Deduplication**: Canonicalizes SMILES, computes
    standard InChIKeys, aggregates multi-source provenance, computes Lipinski
    Ro5 and ADMET descriptors, and exports to SDF, CSV, and interactive HTML
    card views with SVG depictions.

## CLI Examples

```bash
# Standardize raw SAR table with mixed concentration units and inequalities
uv run scripts/sar_table_extractor.py \
  --input-csv /tmp/sar_raw.csv --activity-cols IC50_raw Ki_raw --output-csv /tmp/sar_normalized.csv

# Consolidate and curate multi-source chemical library into interactive HTML report
uv run scripts/document_mol_curator.py \
  --input /tmp/patent_mols.csv --output-csv /tmp/curated.csv --output-html /tmp/curated_report.html
```
