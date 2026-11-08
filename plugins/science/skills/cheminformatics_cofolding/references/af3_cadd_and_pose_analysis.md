# AlphaFold 3 Complex Analysis & CADD Pose Scoring Reference

This guide covers post-processing AlphaFold 3 predicted mmCIF structures,
extracting small-molecule poses to RDKit SDF format, and quantifying binding
pocket interactions.

## 1. Extracting Predicted Ligand Poses

AlphaFold 3 writes complete multimolecular complexes into `.cif` (mmCIF) files.
To isolate the small molecule and reconstruct chemical bond orders from the
input 2D SMILES:

```bash
uv run scripts/extract_cofolded_complex.py \
  --cif predicted_model.cif \
  --confidences summary_confidences_0.json \
  --ligand-smiles "Cc1ccc(NC(=O)c2...)cc1" \
  --output-sdf predicted_ligand_pose.sdf
```

## 2. Confidence Metrics & Interpretation

AlphaFold 3 outputs multi-dimensional confidence scores in
`summary_confidences_*.json`:

-   ****ipTM (Interface pTM)****
    -   **Target Range**: $\ge 0.80$
    -   **Biological / CADD Interpretation**: High confidence in relative
        binding mode and interface packing.
-   ****pLDDT (Pocket Residues)****
    -   **Target Range**: $\ge 70$
    -   **Biological / CADD Interpretation**: High local backbone and side-chain
        coordinate accuracy.
-   ****PAE (Ligand-to-Pocket)****
    -   **Target Range**: $< 5$ Å
    -   **Biological / CADD Interpretation**: Low positional uncertainty between
        ligand atoms and binding pocket.
-   ****has_clash****
    -   **Target Range**: `false`
    -   **Biological / CADD Interpretation**: Verifies no severe steric overlap
        at the protein-ligand interface.

## 3. Pocket Interaction Metrics

`extract_cofolded_complex.py` and `compute_plif.py` automatically measure:

-   **Pocket Residues**: All protein residues with heavy atoms within 4.5 Å of
    any ligand atom.
-   **Potential H-Bonds**: Polar contacts between N/O/S atoms with distance $\le
    3.5$ Å.
-   **Severe Clashes**: Interatomic distances $< 2.0$ Å.
-   **Confidence Summary**: `iptm`, `ptm`, `plddt`, `ranking_score`, and
    `has_clash` are stored in `summary_confidences_*.json` files produced by the
    co-folding engine.
