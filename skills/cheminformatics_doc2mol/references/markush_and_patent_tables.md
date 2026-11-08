# 02. Patent Tables & Markush Assembly Guide

## Overview

Robust graph-based assembly of core scaffolds and substituent R-groups using
RDKit's `Chem.RWMol`.

## Critical Graph Assembly Rules

-   **No SMARTS Reaction Replacement**: Never use `AllChem.ReactionFromSmarts`
    for attachment point substitution because SMARTS enforces strict 1:1
    atom-mapping and fails on complex substituents. Always use
    `MarkushAssembler` (atomic graph recombination via `Chem.RWMol` and
    `Chem.CombineMols`).
-   **Homologous Series Range Expansion**: Condenses aliphatic ranges (e.g.
    `C1-C3 alkyl`, `C1-C2 alkoxy`) into discrete structural isomers (Me, Et,
    nPr, iPr) before SMILES parsing to prevent unclosed ring errors (`[*:2]C1`).
-   **Medicinal Chemistry Shorthand**: Automatically expands abbreviations
    (`iPr`, `tBu`, `cPr`, `Bn`, `Ph`, `Ac`, `NHAc`, `OAc`, `CO2Et`, `CO2Me`,
    `SO2Me`, `Boc`, `Cbz`, `Fmoc`) into valid covalent fragments.
-   **On-the-Fly Property Filtering (Ro5)**: Supports `--max-mw`, `--max-logp`,
    `--max-rotb`, `--max-hbd`, `--max-hba` to prune non-drug-like combinatorial
    explosions during generation.

## CLI Examples

```bash
# Expand Markush table with on-the-fly Ro5 filtering
uv run scripts/markush_assembler.py \
  --core "c1cc([*:1])cc([*:2])c1" \
  --input-csv /tmp/r_groups.csv \
  --max-mw 550 --max-logp 5.0 \
  --output-csv /tmp/expanded.csv --output-sdf /tmp/expanded.sdf
```
