# Coordination Chemistry, Dative Bonds, & 2D Perspective Depiction in RDKit

## 1. The Dative Bond Challenge in Cheminformatics

Classical cheminformatics representations (SMILES, SDF V2000) rely on
fixed-valence Lewis structures. When transition metals ($d$-block) are
represented with standard covalent bonds, RDKit's strict sanitization routine
fails with `KekulizeException` or `Explicit valence exceeds normal`.

To resolve this, RDKit provides `Chem.BondType.DATIVE`:

*   **Directional Coordinate Bond**: Formed from neutral electron-pair donor
    atoms ($\text{PR}_3$, $\text{NR}_3$, $\text{CO}$, ethers) to empty metal
    $d$-orbitals: `[Donor] -> [Metal]`.
*   **Zero Formal Charge Buildup**: Dative bonds allow neutral ligands to
    coordinate without creating unphysical $+1/-1$ formal charge pairs on the
    metal.

## 2. Electron Counting & The MLX Method

Using Green's Covalent Bond Classification (CBC / MLX):

*   **L-type Ligands**: Neutral 2-electron donors ($\text{PR}_3$, amines,
    $\text{CO}$, alkenes).
*   **X-type Ligands**: Anionic 1-electron donors (halides $\text{Cl}^-$, alkyls
    $\text{Me}^-$, hydrides $\text{H}^-$, alkoxides $\text{RO}^-$).
*   **Z-type Ligands**: Neutral 0-electron acceptors (Lewis acids, e.g.
    $\text{BR}_3$).

### Equations:

1.  **Oxidation State (OS)**: \[ OS = N_X + q_{\text{complex}} \]
2.  **d-Electron Count ($d^n$)**: \[ d^n = \text{Group Number} - OS \]
3.  **Valence Electron Count (VEC)**: \[ VEC = d^n + 2 N_L + 2 N_X \]

## 3. Standard CLI Executables

### Organometallic Standardizer & CBC Electron Counter (`organometallic_standardizer.py`)

Converts transition metal complexes to standardized dative bond notation,
computes CBC MLX classification, oxidation state, and 18-electron rule valence
counts, and persists JSON reports to disk:

```bash
uv run scripts/organometallic_standardizer.py \
  --smiles "Cl[Pd](Cl)(P(c1ccccc1)(c2ccccc2)c3ccccc3)P(c4ccccc4)(c5ccccc5)c6ccccc6" \
  --electron_count \
  --output pd_complex_analysis.json
```

### 3D Coordination Complex Builder (`coordination_complex_builder.py`)

Assembles idealized 3D coordination geometries (octahedral, square planar,
tetrahedral, trigonal bipyramidal) around transition metal centers from ligand
SMILES:

```bash
uv run scripts/coordination_complex_builder.py \
  --metal Pd \
  --geometry square_planar \
  --ligands "Cl" "Cl" "P(C1CCCCC1)(C2CCCCC2)C3CCCCC3" "P(C1CCCCC1)(C2CCCCC2)C3CCCCC3" \
  --output pd_pc3_complex.sdf
```

--------------------------------------------------------------------------------

## 4. Organometallic & Coordination Depiction Gotchas

### Conflating Dative Bond Type with Stereochemical Depth (`horizon_perspective`)

*   ❌ **The Bug**: Dashing bonds based on chemical donor type (e.g. rendering
    both $\text{L} \rightarrow \text{M}$ amine arrows as dashed lines while
    leaving covalent $\text{M}-\text{X}$ bonds as plain solid lines).
*   ⚠️ **Why It Fails**: In IUPAC chemical depiction, dashed bonds strictly
    represent depth receding behind the viewing plane ($-z$). Dashing both
    cis-amines while leaving both cis-chlorides in-plane transforms a planar
    square-planar complex into an unphysical, distorted **umbrella / puckered /
    pseudo-tetrahedral geometry** ($158.2^\circ$ trans angle, $0.29\text{ \AA}$
    metal out-of-plane displacement).
*   ✅ **The Fix**: In 3D perspective projection (`--render_style
    horizon_perspective`), the two rear receding bonds are a heterogeneous pair
    (one covalent, one dative), preserving the exact $180.0^\circ$ trans
    linearity and coplanarity.

### SVG Regex Word Boundary Collisions

*   ❌ **The Bug**: `pattern = rf"bond-{b_idx}"` in SVG string replacement.
*   ⚠️ **Why It Fails**: Substring matching without word boundaries causes
    `bond-1` to match `bond-10`, `bond-11`, `bond-12`, etc., unintentionally
    dashing unrelated organic rings (e.g. cyclobutane rings in Carboplatin).
*   ✅ **The Fix**: Always anchor with word boundaries: `pattern =
    rf"\bbond-{b_idx}\b"`.

### OpenSMILES Square-Planar Neighbor Ordering (`@SP1` / `@SP2`)

*   ❌ **The Bug**: Expecting `@SP1` vs `@SP2` to be an absolute configuration
    regardless of branch order in SMILES.
*   ⚠️ **Why It Fails**: OpenSMILES square-planar chiral classes index neighbors
    in the exact order they appear in the SMILES string. Swapping branch orders
    can inadvertently invert cis and trans isomers.
*   ✅ **The Fix**: Call `Chem.AssignStereochemistry(mol, force=True,
    cleanIt=True)` and verify cis ($90^\circ$) vs. trans ($180^\circ$)
    ligand-metal-ligand angles after 2D coordinate generation
    (`rdDepictor.Compute2DCoords`).
