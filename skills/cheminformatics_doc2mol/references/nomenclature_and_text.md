# 01. Nomenclature & Text Sanitization Guide

## Overview

Automated chemical text normalization and multi-tier name-to-structure
resolution using public zero-auth services and local parsing algorithms.

## Key Capabilities & CLI Recipes

-   **OCR Stem Repair**: Corrects optical character recognition confusion (e.g.
    `1 -> l` in `-y1`, `ch1oro`, `pheny1`, `am1ne`, `piperaz1ne`; `0 -> O` in
    `0xazole`, `di0xane`, `benz0ic`; prime quotes on locants like `2'-deoxy`)
    without corrupting numeric positional locants (e.g. `1-methyl`).
-   **IUPAC / Chemical Name Resolution**: Multi-tier resolution cascade:
    Local/REST OPSIN -> PubChem PUG REST -> NIH Chemical Identifier Resolver
    (CIR).
-   **PDF SMILES Typographic Repair**: Repairs broken SMILES from PDF text
    streams (line breaks, soft hyphens `-\n`, intra-SMILES whitespace, unicode
    minus `\u2212`, uncharged nitro `N(=O)=O`, 5-membered aromatic carbocycle
    dropouts).

## CLI Examples

```bash
# Single name resolution
uv run scripts/iupac_resolver.py --name "4-(4-chlorobenzoyl)pyridine" --output-csv /tmp/resolved.csv

# Batch name resolution from CSV (auto-detects chemical name column)
uv run scripts/iupac_resolver.py --input-csv /tmp/names.csv --output-csv /tmp/resolved.csv --output-sdf /tmp/resolved.sdf

# Repair PDF-corrupted SMILES strings
uv run scripts/chemical_text_sanitizer.py --repair-smiles-file /tmp/corrupted.txt --output-csv /tmp/repaired.csv
```
