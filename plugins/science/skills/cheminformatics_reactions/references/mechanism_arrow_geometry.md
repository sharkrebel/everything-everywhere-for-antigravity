# Mathematical Principles of Collision-Free Reaction Mechanism Geometry

This document outlines the geometric, stereochemical, and algorithmic standards
implemented in `cheminformatics_reactions` (`scripts/render_mechanism.py` and
`scripts/rotation_solver.py`) for rendering publication-quality curved
electron-pushing arrows.

--------------------------------------------------------------------------------

## 1. Coordinate Frame & Molecule Canvas Normalization

-   **Metric Units**: All molecular depictions, atom coordinates, glyph
    exclusion radii, and Bézier control points are evaluated in physical
    centimeters (cm) mapped to a standard TikZ coordinate system (`[x=1cm,
    y=1cm]`).
-   **Scale Factor**: A molecule canvas of pixel width $W_{\text{px}}$ and
    metric width $W_{\text{cm}}$ defines the uniform conversion:

    $$\text{scale} = \frac{W_{\text{cm}}}{W_{\text{px}}}$$

-   **Atom Position Extraction**: Atom center coordinates $(c_x,\,c_y)$ are
    extracted directly from RDKit's 2D depictor via `d2d.GetDrawCoords(i)`:

    $$c_x = x_{\text{canvas}} + (p_x - W_{\text{px}}/2) \times \text{scale}$$

    $$c_y = y_{\text{canvas}} - (p_y - H_{\text{px}}/2) \times \text{scale}$$

-   **Glyph Exclusion Disks**: Every non-hydrogen or charged atom possesses a
    radial exclusion buffer $r_{\text{glyph}}$:
    -   Heteroatoms / charged centers: $r_{\text{glyph}}=0.20$–$0.26$ cm
    -   Skeletal carbon vertices: $r_{\text{glyph}}=0.12$ cm. Arrows must never
        originate or terminate within $r_{\text{glyph}}$ to prevent stroke
        collision with atom labels.

--------------------------------------------------------------------------------

## 2. Inter-Molecular Attack Curves (`compute_attack_bezier`)

Used when a nucleophilic lone pair or bond attacks an electrophilic atom across
reagents (e.g. over a $+$ sign or reaction coordinate gap).

### Geometric Invariants:

1.  **Departure Angle**: Departs along the nucleophile's open back-vector
    $\mathbf{u}_{\text{back}}$ pointing away from internal substituent bonds:

    $$\begin{aligned}
      \mathbf{u}_{\text{back}} = -\frac{
        \sum_{j \in \text{nbrs}} (\mathbf{p}_j - \mathbf{p}_{\text{src}})
      }{
        \left\|
          \sum_{j \in \text{nbrs}} (\mathbf{p}_j - \mathbf{p}_{\text{src}})
        \right\|
      }
    \end{aligned}$$

    $$\begin{aligned}
      p_0 &= \mathbf{p}_{\text{src}} \\
      &\quad + (r_{\text{glyph}} + 0.04\,\text{cm})\,\mathbf{u}_{\text{dep}}
    \end{aligned}$$

2.  **Arrival Angle**: Enters the electrophilic center along its open quadrant
    $\mathbf{u}_{\text{open}}$, landing on the glyph perimeter:

    $$\begin{aligned}
      p_3 &= \mathbf{p}_{\text{tgt}} \\
      &\quad + (r_{\text{glyph}} + 0.04\,\text{cm})\,\mathbf{u}_{\text{arr}}
    \end{aligned}$$

3.  **Obstacle Clearance (Plus Signs)**:
    -   For attack crossing a $+$ sign located at $(x_{+},\,0)$, the trajectory
        must pass either above ($y>+0.30$ cm) or below ($y<-0.30$ cm).
    -   The candidate path generator evaluates both "over" and "under" arches,
        penalizing any line segment intersecting the bounding box
        $[x_{+}\pm{0.25},\,\pm{0.30}]$.

--------------------------------------------------------------------------------

## 3. Double-Bond & Polar Bond Cleavage (`compute_bond_cleave_bezier`)

Used for cleavage of $\pi$-electrons or polar $\sigma$-bonds onto an
electronegative atom (e.g., $\mathrm{C{=}O}\to\text{C–O}^-$,
$\text{C–X}\to\text{C}^++\text{X}^-$).

### The Double-Bond Overlap Bug & The Normal Offset Invariant:

-   In RDKit 2D depictions (`multipleBondOffset` of `0.15`), centered double
    bonds (e.g., `Fe=O`, ketones) draw two parallel strokes at
    $\pm{0.075}\,L_{\text{bond}}$ ($\approx\pm{0.06}$–$0.08$ cm), ring and
    asymmetric double bonds draw strokes at $0$ and
    $\pm{0.15}\,L_{\text{bond}}$, and triple bonds draw strokes at $0$ and
    $\pm{0.15}\,L_{\text{bond}}$.
-   ChemDraw/TikZ arrow strokes possess a finite line width ($1.0$ pt,
    $\approx{0.035}$ cm).
-   If the arrow origin $p_0$ is measured only from the bond centerline, any
    normal offset $\le{0.08}$ cm physically grazes or overlaps the outer
    parallel stroke of a double or triple bond.
-   **Invariant**: `MoleculeCanvas` computes the exact `stroke_segments` for
    each single, double, and triple bond, and all bond-originating or
    bond-targeting arrow builders (`compute_bond_cleave_bezier`,
    `compute_attack_bezier`, `compute_bond_to_bond_bezier`,
    `compute_lone_pair_to_bond_bezier`) add the outer parallel-stroke offset
    `bond_outer_offset_cm(b, n)` to the base clearance ($0.105$–$0.14$ cm),
    guaranteeing a crisp gap of pure whitespace outside all drawn bond strokes.

--------------------------------------------------------------------------------

## 4. Bond-to-Bond Transitions (`compute_bond_to_bond_bezier`)

Used for concerted electron shifts: $\alpha$-deprotonations
($\text{C–H}\to\mathrm{C{=}C}$), eliminations, allylic shifts, and pericyclic
cascades.

-   **Ring Inward Curvature**: Inward curvature toward `ring_center` is ONLY
    valid when **both** bonds $b_1$ and $b_2$ belong to the same ring (e.g.
    cyclopentadiene in Diels-Alder reactions).
-   **Exocyclic Transitions**: If either bond is exocyclic (e.g., $\text{C–H}$
    deprotonation forming an enolate or alkene), the arrow must **never** loop
    into or through the ring cavity. It is strictly restricted to the open acute
    exterior sector between the two bonds.

--------------------------------------------------------------------------------

## 5. Reaction Orientation Solver (`scripts/rotation_solver.py` & `--auto-rotate`)

-   By default, `render_mechanism.py` draws each molecule at the rotation in the
    spec's `rot` map, else the molecule entry's own `rot`, else $0^\circ$.
-   Pass `--auto-rotate` on the CLI (or `"auto_rotate": true` in JSON/YAML) to
    automatically invoke `scripts/rotation_solver.py`
    (`solve_reaction_rotations`) for any molecule without an explicit `rot`
    entry.
-   It scores candidates across the 12-angle ChemDraw/IUPAC candidate space
    $\Theta=\{0^\circ,\,30^\circ,\,\dots,\,330^\circ\}$ so $R_1$'s reacting site
    faces rightward ($+x$) toward $R_2$, $R_2$'s reacting site faces leftward
    ($-x$) toward $R_1$ (aligning multi-locus cycloaddition partners such as
    $[4+2]$ Diels–Alder and $[3+2]$ dipolar cycloadditions horizontally across
    all intermolecular chords), and product scaffolds ($P_1$, $P_2$) are aligned
    to $R_1$, $R_2$ via Maximum Common Substructure (MCS) minimum-RMSD
    alignment.

--------------------------------------------------------------------------------

## 6. Explicit Hydrogens (`add_h`), Dummy Labels (`atom_labels`), & Amino Acid Side-Chain Identity Gate

-   **Explicit Hydrogens (`add_h`)**: For proton transfers, deprotonations, and
    hydride additions, pass `add_h: {"r1": parent_atom_idx}` so an explicit
    hydrogen atom is appended at index `N_atoms`. `MoleculeCanvas` unfolds
    crowded intramolecular pockets ($d<2.35$ Å) and orients the added hydrogen
    directly toward the attacking atom. Omit `"side"` on short intramolecular
    proton-transfer attack arrows so `compute_attack_bezier` selects the direct
    interior arc.
-   **Macromolecular Attachments (`atom_labels`)**: When using dummy atoms (`*`)
    for macromolecular attachments (`tRNA`, `Peptidyl`, `Enz`), provide
    `atom_labels: {atom_idx: "Label"}` on every molecule/intermediate.
-   **Amino Acid Side-Chain Identity Gate (`_validate_amino_acid_sidechain`)**:
    `verify_mechanism_with_electron_graph` validates that any molecule whose
    `name` or `atom_labels` names a canonical amino acid residue uses its
    authentic chemical functional group (rejecting mismatched surrogates like
    `*[NH3+]` for `His` or `*[O-]` for `Asp`/`Glu`/`Lys`):
    -   **His**: `*c1c[nH+]c[nH]1` (`His-H+` acid) / `*c1cnc[nH]1` (`His` base)
    -   **Lys**: `*CC[NH3+]` or `*[NH3+]` (`Lys-H+`) / `*CCN` or `*N` (`Lys`)
    -   **Arg**: `*NC(=[NH2+])N` (`Arg-H+`) / `*NC(=N)N` (`Arg`)
    -   **Asp / Glu**: `*CC(=O)[O-]` or `*C(=O)[O-]` (base) / `*CC(=O)O` or
        `*C(=O)O` (acid)
    -   **Cys**: `*C[S-]` or `*[S-]` / `*CS` or `*S`
    -   **Ser / Thr**: `*CO` or `*O` (`*C(C)O` for `Thr`) / `*C[O-]` or `*[O-]`
    -   **Tyr**: `*c1ccc(O)cc1` (`Tyr-OH`) / `*c1ccc([O-])cc1` (`Tyr-O-`)
    -   **Asn / Gln / Trp / Phe / Met / Pro**: `*CC(=O)N`, `*c1c[nH]c2ccccc12`,
        `*Cc1ccccc1`, `*CCSC`, `*C1CCCN1`

--------------------------------------------------------------------------------

## 7. Intramolecular Macrocyclization Folding & Multi-Arrow Coordination

-   **Stereoretentive $180^\circ$ Rotamer Folding**
    (`_fold_macrocycle_rotamers`): RDKit's default 2D depictor lays out long
    acyclic chains (`len(path) >= 5`, such as lorlatinib-class macrocycle
    precursors) in extended zig-zag conformations where the nucleophile and
    electrophile sit far apart across intervening aryl or amide bonds.
    `MoleculeCanvas` automatically searches exact $180^\circ$ reflections across
    acyclic single bonds along the shortest intramolecular path (preserving all
    bond lengths, internal bond angles, and restricted atropisomer bonds while
    clearing cached wedges so RDKit recomputes CIP-faithful wedges/dashes) to
    form an unobstructed C-shaped reactive bay.
-   **Open-Sector Bisectors & Stealth Wing Clearance**: `compute_attack_bezier`
    supplements horizontal `over`/`under` and chord `interior` arcs with
    open-sector bisector candidates (`_open_sector_bisector_degs`) so vertical
    or tilted intramolecular bays never hit a $\pm{180}^\circ$ angular dead
    zone, and penalizes Stealth arrowhead wing segments
    (`_arrowhead_wing_segments`) against non-endpoint bond strokes.
-   **Priority Multi-Arrow Routing & Quantitative Pairwise Audit**:
    `compile_clean_reaction` routes arrows in topological priority order
    (`lone_pair_to_bond` / intramolecular `bond_to_bond` $\to$ intramolecular
    `attack` $\to$ top-to-bottom intermolecular `attack` / `bond_to_atom` $\to$
    `bond_cleave`) while emitting TikZ and audit records in specification order,
    passing previously routed polylines (`prior_arrows`) and incoming driver
    vectors (`attack_dir`) so concerted cascades and tetrahedral collapses stay
    cleanly separated. `audit_arrow_clearances` reports pairwise
    `arrow_clearance_mm` / `closest_arrow_pair_mm` and gates `all_clear` on
    $\ge{0.50}$ mm arrow-to-arrow separation; CLI invocations (`--rxn-json`,
    `--rxn-yaml`, `--id`) run with `strict_audit=True` and exit non-zero
    whenever `not audit["all_clear"]`.
-   **Paired $1\text{e}^-$ Nascent-Bond Fishhooks**: When two `"electrons": 1`
    arrows (`attack`, `bond_to_bond`, or incipient `lone_pair_to_bond`)
    contribute single electrons to the same nascent bond (e.g., Baldwin
    5-exo-trig radical cyclizations or H-atom transfer steps),
    `_compute_paired_fishhook_beziers` routes both harpoon arrows to converge
    symmetrically at the nascent-bond midpoint with non-crossing convex barbs.
-   **Planar Cycloaddition Face Uncrossing & Dynamic Caption Clearance**: When
    two vertically aligned intermolecular chords in a bimolecular reaction cross
    in vertical endpoint order (`dy1 * dy2 < 0`), `compile_clean_reaction`
    reflects an achiral planar reactant canvas vertically (`flip_y=True`) so the
    cycloaddition arcs frame the reaction bay symmetrically above and below
    without crossing, and shifts any molecule caption dynamically below
    under-arcing trajectories (`_caption_below_str`).

--------------------------------------------------------------------------------

## 8. Antigravity Visual Styling & Crimson Arrow Aesthetics

Curved electron-pushing arrows are rendered in Antigravity Crimson (`#E11D48`,
`\definecolor{antigravitycrimson}{HTML}{E11D48}`) with refined
`-Stealth[length=2.0mm, width=1.30mm, flex=0.75]` arrowheads (`chemarrow` for
2-electron pairs) and prominent single-barbed `-Stealth[harpoon, length=2.45mm,
width=2.25mm, flex=0.75]` arrowheads (`fishhookarrow` / `fishhookarrowswap` for
`"electrons": 1` radical steps, automatically placing the single lateral barb on
the convex exterior side of the curve).
