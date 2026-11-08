# ChemDraw CDXML/CDX Interoperability & ACS Document 1996 Styling

This reference covers reading, writing, auditing, and automating ChemDraw files
(`.cdxml`, `.cdx`) using `scripts/chemdraw_cdxml.py` and
`../cheminformatics_reactions/scripts/reaction_to_cdxml.py`.

--------------------------------------------------------------------------------

## 1. Publication Styling Presets (`ACS Document 1996`, `RSC`, `Wiley`)

All spatial measurements in CDXML are in points ($1\text{ pt} =
\frac{1}{72}\text{ inch} = 0.3528\text{ mm}$).

-   **`acs1996` (ACS Document 1996 — Default)**:
    -   **Chain Angle**: $120^\circ$
    -   **Bond Length**: `14.4 pt` (5.08 mm)
    -   **Bond Spacing**: `18%` of length (2.59 pt)
    -   **Line Width**: `0.6 pt` (0.212 mm)
    -   **Bold Width**: `2.0 pt` (0.706 mm)
    -   **Hash Spacing**: `2.5 pt` (0.88 mm)
    -   **Margin Width**: `1.6 pt` (0.56 mm)
    -   **Atom Label Font**: Arial (`font="3"`), `10 pt`
    -   **Caption Font**: Times New Roman (`font="4"`), `10 pt`
-   **`rsc` (Royal Society of Chemistry)**:
    -   **Chain Angle**: $120^\circ$
    -   **Bond Length**: `15.0 pt` (5.30 mm)
    -   **Bond Spacing**: `18%` (2.70 pt)
    -   **Line Width / Bold / Hash / Margin**: `0.6 pt` / `2.0 pt` / `2.5 pt` /
        `1.6 pt`
    -   **Atom Label & Caption Fonts**: Arial `9 pt` / Times New Roman `9 pt`
-   **`wiley` (Wiley Journals)**:
    -   **Chain Angle**: $120^\circ$
    -   **Bond Length**: `17.0 pt` (6.00 mm)
    -   **Bond Spacing**: `18%` (3.06 pt)
    -   **Line Width / Bold / Hash / Margin**: `0.6 pt` / `2.0 pt` / `2.5 pt` /
        `1.6 pt`
    -   **Atom Label & Caption Fonts**: Arial `10 pt` / Times New Roman `10 pt`

--------------------------------------------------------------------------------

## 2. Mandatory CDXML Geometry & Valence Rules

1.  **Vertical Coordinate Inversion**: In CDXML, the origin $(0, 0)$ is at the
    top-left and $y$ increases **downwards**: $$x_{\text{cdxml}} =
    x_{\text{cartesian}} \times \text{scale} + x_{\text{offset}}, \qquad
    y_{\text{cdxml}} = -y_{\text{cartesian}} \times \text{scale} +
    y_{\text{offset}}$$
2.  **Heteroatom `NumHydrogens` (No Red Valence Warning Boxes)**: When exporting
    molecules to CDXML, ChemDraw evaluates atom valences from explicit `<b>`
    bonds plus the `NumHydrogens` attribute on `<n>`. Every non-carbon atom
    (`Element != 6`, e.g. `N`, `O`, `S`, `P`) with implicit hydrogens MUST
    declare `NumHydrogens="{atom.GetTotalNumHs()}"`. Without it, ChemDraw treats
    $-\text{NH}_2$, $-\text{OH}$, $-\text{NH}-$, and $-\text{SH}$ as open
    radicals and paints bright red warning boxes.
3.  **Explicit Kekulization (`Order="2"`)**: ChemDraw Desktop does not render
    alternating double bonds for `Order="1.5"` aromatic query bonds. Always call
    `Chem.Kekulize(mol, clearAromaticFlags=True)` before emitting `<b>` elements
    so double bonds carry `Order="2"`.
4.  **Dual Reaction Arrow Tags (Desktop + Web Compatibility)**: ChemDraw Desktop
    renders visible arrows from `<graphic GraphicType="Line"
    ArrowType="FullHead">`, whereas ChemDraw JS and
    `rdChemDraw.ReactionsFromChemDrawFile` parse `<step>` and `<arrow>`. Always
    emit both:
    -   In `<graphic GraphicType="Line" ArrowType="FullHead">`, `<point>` #1 is
        the **arrow head** (`arrow_head_x`) and `<point>` #2 is the **arrow
        tail** (`arrow_tail_x`), so the arrowhead points rightward toward the
        product.
    -   In `<step ReactionType="Normal">`, include `<arrow ArrowheadHead="Full"
        Head3D="..." Tail3D="..."/>` and `<objects><object id="..."/></objects>`
        referencing all reactant and product `<fragment>` IDs.
5.  **Collision-Free Reaction Scheme Clearances (`reaction_to_cdxml.py`)**:
    -   **Arrow length**: $\ge 115.0\text{ pt}$, dynamically expanded to
        $\max(115.0, \text{text\_width} + 30.0\text{ pt})$.
    -   **Reactant/Product clearance**: $\ge 40.0\text{ pt}$ from arrow tail and
        arrow head; $25.0\text{ pt}$ around `+` signs.
    -   **Multi-line reagent wrapping**: Wrap reagent/condition strings to $\le
        35$ characters per line using native CDXML newlines (`&#10;`) inside
        `<s>` with `LineHeight="auto"` and `Justification="Center"`.
    -   **Dynamic compound labels**: Place bottom compound labels at $Y =
        Y_{\text{bbox\_max}} + 18.0\text{ pt}$.
    -   **Multi-row separation**: Separate rows in multi-step schemes by $\ge
        280.0\text{ pt}$ vertically.

--------------------------------------------------------------------------------

## 3. CLI & Python Recipes

### Single Molecule SMILES/Molfile $\to$ Styled CDXML

```bash
uv run scripts/chemdraw_cdxml.py \
    --smiles "CC(=O)Nc1ccc(O)cc1" \
    --style acs1996 \
    --output scratch/paracetamol.cdxml
```

### Parse `.cdxml` / `.cdx` to SMILES, SDF, JSON, or Reaction SMILES

```bash
uv run scripts/chemdraw_cdxml.py --input scheme.cdxml --format smiles --output scratch/mols.smi
uv run scripts/chemdraw_cdxml.py --input scheme.cdxml --format rxn --output scratch/rxns.smi
```

### Audit a `.cdxml` File for Valence Warnings & Dual Arrows

```bash
uv run scripts/chemdraw_cdxml.py --inspect --input scratch/scheme.cdxml
```

### Multi-Step Reaction Scheme $\to$ Styled CDXML (`cheminformatics_reactions`)

```bash
uv run ../cheminformatics_reactions/scripts/reaction_to_cdxml.py \
    --json scratch/synthesis_steps.json \
    --output scratch/synthesis_scheme.cdxml
```

--------------------------------------------------------------------------------

## 4. Desktop Automation (Windows COM & macOS AppleScript)

When running locally on a Windows or macOS workstation with ChemDraw Desktop
installed:

-   **Windows (`pywin32`)**:

    ```python
    import win32com.client
    cd_app = win32com.client.Dispatch("ChemDraw.Application")
    doc = cd_app.Documents.Open(r"C:\path\to\scheme.cdxml")
    doc.Objects.Select()
    doc.Objects.CleanUpStructure()
    doc.SaveAs(r"C:\path\to\scheme.pdf")
    doc.Close()
    ```
-   **macOS (`osascript`)**:

    ```bash
    osascript -e 'tell application "ChemDraw" to open POSIX file "/path/to/scheme.cdxml"'
    ```
