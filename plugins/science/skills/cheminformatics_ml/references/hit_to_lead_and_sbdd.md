# Hit-to-Lead SAR, Matched Molecular Pairs & Series Partitioning

This guide covers medicinal chemistry SAR analysis, matched molecular pair
analysis (MMPA), activity cliffs (SALI), diversity selection, positional
analogue scanning, and chemical clustering using the production tools in
`scripts/`.

--------------------------------------------------------------------------------

## 1. Matched Molecular Pair Analysis (MMPA) & Activity Cliffs (SALI)

Identifies structural transformations between pairs of compounds and quantifies
their impact on potency:

$$\text{SALI}_{ij} = \frac{|\text{Activity}_i - \text{Activity}_j|}{1 - \text{Tanimoto}(i, j)}$$

### Production CLI:

```bash
uv run scripts/mmpa_sar_analysis.py \
    --input_csv scratch/series.csv \
    --smiles_col SMILES \
    --activity_col pIC50 \
    --min_similarity 0.6 \
    --output_csv scratch/mmpa_activity_cliffs.csv
```

*   Identifies single-cut chemical transformations (e.g. `[H] >> F`, `Me >>
    OMe`).
*   Computes SALI scores to prioritize sharp activity cliffs over smooth SAR.

--------------------------------------------------------------------------------

## 2. Positional Analogue Scanning (PAS)

Systematically walks functional group substituents (Fluorine, Methyl, Chlorine,
Hydroxyl, Methoxy, CF3, Nitrogen-walk) across all aromatic C–H vectors:

### Production CLI:

```bash
uv run scripts/positional_analogue_scanner.py \
    --smiles "c1ccc2c(c1)nccn2" \
    --scan_type fluorine \
    --output_csv scratch/fluorine_walk.csv
```

*   Supported scans: `fluorine`, `methyl`, `chlorine`, `hydroxyl`, `methoxy`,
    `trifluoromethyl`, `nitrogen_walk`.
*   Calculates exact property deltas ($\Delta\text{MW}$, $\Delta\text{cLogP}$,
    $\Delta\text{TPSA}$) for each analogue.

--------------------------------------------------------------------------------

## 3. Bemis-Murcko Scaffold & Chronological Temporal Partitioning

Partitions compound libraries by their 2D molecular frameworks or strict
chronological publication/patent cutoffs to eliminate train/test leakage:

### Production CLI (Scaffold Split):

```bash
uv run scripts/scaffold_splitter.py \
    --input_csv scratch/library.csv \
    --smiles_col SMILES \
    --train_ratio 0.8 \
    --val_ratio 0.0 \
    --test_ratio 0.2 \
    --output_csv scratch/library_with_splits.csv
```

### Production CLI (Strict Cutoff Temporal Split):

```bash
uv run scripts/scaffold_splitter.py \
    --input_csv scratch/library.csv \
    --split-method time \
    --time-col document_year \
    --time-cutoff 2013 \
    --output_csv scratch/library_temporal_splits.csv
```

*   Writes a single CSV: the input rows plus a `Split` column
    (`train`/`val`/`test`) with zero Bemis-Murcko scaffold overlap (scaffold
    mode) or zero timestamp overlap $\max(t_{\text{train}}) <
    \min(t_{\text{test}})$ (temporal mode: `train <= cutoff`, `test > cutoff`).
*   Automatically rejects boolean `False`/`True` placeholders, backfills missing
    dates from companion document/patent metadata (`publication_date`,
    `filing_date`, `priority_date`, sibling `chembl_id`/`document_id`, or
    WIPO/US-A1 patent publication years), and prevents tied years from splitting
    across train and test sets.

--------------------------------------------------------------------------------

## 4. Chemical Clustering & Centroid Selection (Butina)

Clusters compound collections using Taylor-Butina clustering on Morgan circular
fingerprints (ECFP4) and identifies cluster centroids:

### Production CLI:

```bash
uv run scripts/cluster_compounds.py \
    --input_csv scratch/hits.csv \
    --smiles_col SMILES \
    --cutoff_distance 0.4 \
    --output_csv scratch/clustered_hits.csv
```

*   Groups molecules by chemical similarity (Tanimoto distance threshold $1 -
    \text{Sim}$).
*   Assigns `Cluster_ID` and marks `Is_Centroid=True` for representative cluster
    leads.

--------------------------------------------------------------------------------

## 5. MaxMin Diversity Selection

Picks a maximally diverse subset from a large compound library using C++
accelerated `LazyBitVectorPick`:

### Production CLI:

```bash
uv run scripts/diversity_picker.py \
    --input_csv scratch/deck.csv \
    --num_picks 50 \
    --smiles_col SMILES \
    --output_csv scratch/diverse_subset.csv
```

*   Guaranteed maximal Tanimoto diversity across selected compounds.
*   Scales to $100,000+$ candidate molecules without storing dense distance
    matrices.

--------------------------------------------------------------------------------

## 6. Ring System & Stereoisomer Analysis

*   **Ring System Enumeration**: Use `scripts/ring_system_analyzer.py` to
    dissect ring frameworks and identify novel ring topologies.
*   **Chiral Stereoisomer Enumeration**: Use
    `scripts/stereoisomer_enumerator.py` to enumerate all unassigned
    stereocenters and double bond geometries.

--------------------------------------------------------------------------------

## 7. Structure-Based & 3D Modeling Cross-Reference

For 3D physical modeling, receptor-ligand docking, conformational strain,
relaxed dihedral PES scans, and macrocycle/PROTAC preorganization, refer to the
**`cheminformatics_physics`** companion skill (`references/`).
