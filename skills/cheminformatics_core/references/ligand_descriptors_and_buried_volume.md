# Ligand Steric & Electronic Descriptors: SambVca %V_bur, Tolman & Sterimol

## 1. SambVca Buried Volume (%V_bur)

The Percent Buried Volume ($\%V_{\text{bur}}$) measures the fraction of a sphere
of radius $R$ centered on the metal atom occupied by the ligand: \[
\%V_{\text{bur}} = \frac{V_{\text{ligand\_in\_sphere}}}{V_{\text{sphere}}}
\times 100\% \]

*   **Standard Cavallo Parameters**: Sphere radius $R = 3.5\text{ \AA}$,
    $d_{\text{M-L}} = 2.28\text{ \AA}$, Bondi VdW radii scaled by $1.17$.
*   **Benchmark Values**:
    *   $\text{PPh}_3$: $31.5\%$
    *   $\text{PCy}_3$: $36.2\%$
    *   $\text{P}(t\text{-Bu})_3$: $44.1\%$
    *   $\text{IPr}$ (NHC): $51.2\%$

## 2. Tolman Cone Angle (θ) & Sterimol (L, B1, B5)

*   **Tolman Cone Angle**: Apex angle at the metal vertex enclosing the van der
    Waals boundary of the ligand.
*   **Sterimol $L$**: Total ligand length along the coordination bond axis.
*   **Sterimol $B_1, B_5$**: Minimum and maximum perpendicular widths.

## 3. Standard CLI Executables

### Quantitative Ligand Sterics Engine (`ligand_steric_descriptors.py`)

Computes sphere-grid buried volume ($\%V_{\text{bur}}$), 4-quadrant breakdown,
Tolman cone angle $\theta$, and Sterimol $(L, B_1, B_5)$ descriptors, saving
JSON output to a file:

```bash
uv run scripts/ligand_steric_descriptors.py \
  --smiles "P(C1CCCCC1)(C2CCCCC2)C3CCCCC3" \
  --metal_dist 2.28 \
  --output pcy3_sterics.json
```

### 2D Topographic Steric Map Renderer (`render_steric_map.py`)

Renders a publication-quality 2D topographic steric contour map (altitude map of
ligand steric bulk around the metal coordination sphere) to PNG/SVG:

```bash
uv run scripts/render_steric_map.py \
  --smiles "P(C1CCCCC1)(C2CCCCC2)C3CCCCC3" \
  --output pcy3_steric_map.png
```
