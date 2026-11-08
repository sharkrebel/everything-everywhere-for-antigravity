# Building Standalone Applications, Dashboards, and Visualizers with the Cheminformatics Skill Suite

When a user asks to generate a standalone chemical application, interactive
molecule visualizer, SAR explorer, reaction dashboard, or web app, **do not
build the chemistry, descriptor, or molecular rendering stack from scratch**.
Always compose the application on top of `rdkit-general` and the companion
`cheminformatics-*` skills.

--------------------------------------------------------------------------------

## 1. Existing Visualization Tools Across the Skill Suite

Before writing custom UI or rendering code, check whether these existing
visualization scripts already cover the view you need:

-   **`rdkit-general`** — Interactive HTML 2D compound grids/SAR scatter plots
    and ACS-styled 2D structure depictions
    (`scripts/interactive_grid_viewer.py`, `scripts/draw_molecules.py`).
-   **`cheminformatics-reactions`** — Curved-arrow reaction mechanism diagrams
    and HTML reaction galleries
    (`../cheminformatics_reactions/scripts/render_mechanism.py`,
    `../cheminformatics_reactions/scripts/build_gallery.py`).
-   **`cheminformatics-ml`** — Atom-level similarity attribution maps and
    multi-parameter radar plots
    (`../cheminformatics_ml/scripts/similarity_map_generator.py`,
    `../cheminformatics_ml/scripts/mpo_pareto_radar.py`).

--------------------------------------------------------------------------------

## 2. Turnkey Zero-CDN Standalone HTML Visualizers

Sandboxed artifact webviews enforce strict Content Security Policy (CSP) rules
that block external JavaScript CDNs (`unpkg.com`, `cdn.jsdelivr.net`,
`3dmol.org`, `cdnjs.cloudflare.com`). Relying on external JS libraries for
SMILES drawing or 3D molecular rendering results in blank panels.

### 2D Compound Explorer & SAR Scatter Dashboard

Use `scripts/interactive_grid_viewer.py` to generate a self-contained, zero-CDN,
CSP-compliant HTML application featuring:

-   Dual-view toggle between a responsive 2D structure card grid and an
    interactive SAR scatter plot (`MW` vs. `cLogP` with Lipinski Rule-of-5
    sweet-spot overlay).
-   Automatic calculation of `MW`, `cLogP`, `TPSA`, `HBD`, `HBA`, `RotB`, `QED`,
    Pfizer CNS MPO (`calculate_cns_mpo`), and PAINS/Brenk structural alert
    badges.
-   Common-core scaffold alignment (`--align_core`), smart search/filtering
    (multi-term `AND`, `OR`, `-` exclusion, quoted phrases, `field:value`,
    numeric comparisons like `mw<450`, `pki>9`, `ic50<=10`, and ranges like
    `mw:300..450`), modal zoom with one-click SMILES copy, and client-side CSV
    export.

```bash
# 1. First create the placeholder artifact file via write_to_file (UserFacing: True)
# 2. Populate the standalone HTML app directly in the artifact directory:
uv run scripts/interactive_grid_viewer.py \
  --input curated_series.csv \
  --smiles_col canonical_smiles \
  --id_col compound_id \
  --activity_col pIC50 \
  --align_core "c1ccc2ncnc(N)c2c1" \
  --title "Kinase Series SAR Explorer" \
  --output <artifact_dir>/kinase_sar_explorer.html
```

--------------------------------------------------------------------------------

## 3. Integrating the Skill Suite into Custom Python Applications & Backends

When building a custom multi-tab web application, REST backend (FastAPI /
Flask), interactive dashboard (Streamlit / Gradio / Dash), or custom single-file
HTML generator that requires bespoke UI workflows beyond the turnkey viewers,
**always import the skill modules as your cheminformatics backend engine**
rather than writing raw RDKit boilerplate from scratch.

### Step 1: Add the Skill Suite to `PYTHONPATH`

Every script in the suite has a unique module basename and can be imported
directly once the colon-separated path list printed by `scripts/skill_paths.py`
is prepended to `PYTHONPATH`. Set `PYTHONPATH` in the **same command** that
launches your script or server (agent harnesses execute each shell call in a
fresh subshell, so a standalone `export` will not persist), and note that
`scripts/skill_paths.py` is relative to the `rdkit_general` skill directory
(pass the full path to `skill_paths.py` when running from another working
directory):

```bash
PYTHONPATH="$(uv run scripts/skill_paths.py)${PYTHONPATH:+:$PYTHONPATH}" uv run app.py
```

Inside your application (`app.py`), import the verified helpers directly at the
top of the file:

```python
from compat import get_morgan_fingerprint, safe_substruct_match, standardize_molecule
from descriptors import calculate_cns_mpo, compute_standard_descriptors
from filter_catalogs import screen_structural_alerts
from mmpa_sar_analysis import find_matched_molecular_pairs
```

### Step 2: Use Verified Pipeline Functions Instead of Ad-Hoc RDKit Calls

-   **Molecule Ingestion & Standardization**: Always pass user-supplied SMILES
    through `standardize_molecule` (`scripts/compat.py`) or
    `../cheminformatics_core/scripts/standardize_molecule.py` to strip
    counterions/salts, normalize functional groups, and recover sanitization
    errors before computing properties or rendering depictions.
-   **Descriptor & Multiparameter Profiling**: Call
    `compute_standard_descriptors` and `calculate_cns_mpo` from
    `../cheminformatics_core/scripts/descriptors.py` so column names (`MW`,
    `cLogP`, `TPSA`, `HBD`, `HBA`, `RotB`, `Fsp3`, `QED`, `SA_Score`) conform to
    the canonical schema across all tables, plots, and JSON APIs.
-   **MedChem Liability Screening**: Call `screen_structural_alerts` from
    `../cheminformatics_core/scripts/filter_catalogs.py` to attach PAINS (A/B/C)
    and Brenk alert annotations and Rule-of-5 / Rule-of-3 / Veber compliance
    flags to every compound card.
-   **Substructure & Similarity Filtering**: Use `safe_substruct_match` and
    `get_morgan_fingerprint` from `scripts/compat.py` (which uses RDKit 2024.03+
    `rdFingerprintGenerator.GetMorganGenerator` instead of deprecated `AllChem`
    calls).

--------------------------------------------------------------------------------

## 4. Generating Responsive, Theme-Adaptive 2D SVGs & Sizing Plots for Web UIs

When embedding 2D chemical structures or analytical figures into custom HTML/CSS
cards, interactive tables, or dashboards:

1.  **Use CoordGen & Scaffold Alignment**: Enable CoordGen
    (`rdDepictor.SetPreferCoordGen(True)`) for clean macrocycle and polycyclic
    layouts, and align congeneric series to their shared maximum common
    substructure (MCS) or user-specified SMARTS core using
    `rdDepictor.GenerateDepictionMatching2DStructure` so substituents appear in
    consistent orientations across the UI.
2.  **Dark-Mode Structure Legibility Invariant (No Black Bonds on Dark
    Surfaces)**: Standard ACS 1996 RDKit depictions draw carbon bonds and atom
    labels in black (`#000000`). **Never** place a transparent black-stroke SVG
    or PNG (`clearBackground = False` / `--transparent` alone) directly onto a
    dark UI surface (`var(--card)` / `#1e293b` or `#0f172a`). Follow one of
    these two safe patterns:
    -   **Default (Recommended — White Structure Well)**: Keep
        `clearBackground = True` (or wrap the SVG in `.chem-svg-well {
        background: #ffffff; border-radius: 8px; padding: 8px; border: 1px solid
        var(--border); }`) so black bonds and CPK heteroatom colors retain crisp
        publication contrast in both light and dark UI modes.
    -   **Explicit Dark Canvas (Opt-In)**: If rendering directly onto a dark
        background without a white well, enable RDKit's dark-mode palette via
        `rdMolDraw2D.SetDarkMode(opts)` (or pass `--dark-mode --transparent` to
        `scripts/draw_molecules.py`).
3.  **Responsive SVG Sizing**: Strip fixed pixel `width` and `height` attributes
    from the root `<svg>` element while preserving `viewBox="0 0 W H"` so the
    SVG scales fluidly inside CSS Grid or Flexbox containers:

```python
import re
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D


def render_responsive_svg(
    mol: Chem.Mol,
    width: int = 280,
    height: int = 200,
    highlight_atoms: list[int] | None = None,
    dark_mode: bool = False,
) -> str:
  """Renders a responsive ACS-styled SVG string for web UIs.

  When dark_mode=False (default), keeps an opaque white background so black
  (#000000) bonds remain legible in both light and dark themes (wrap in a
  .chem-svg-well container). When dark_mode=True, switches to light strokes via
  rdMolDraw2D.SetDarkMode on a transparent canvas.
  """
  rdDepictor.SetPreferCoordGen(True)
  if not mol.GetNumConformers():
    rdDepictor.Compute2DCoords(mol)
  drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
  opts = drawer.drawOptions()
  if dark_mode:
    rdMolDraw2D.SetDarkMode(opts)
  opts.clearBackground = not dark_mode
  opts.bondLineWidth = 1.8
  opts.multipleBondOffset = 0.18
  drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms or [])
  drawer.FinishDrawing()
  svg = drawer.GetDrawingText()
  svg = re.sub(r"<\?xml[^>]*\?>\s*", "", svg)
  svg = re.sub(r"<!DOCTYPE[^>]*>\s*", "", svg)
  svg = re.sub(
      r'(<svg\b[^>]*?)\s+width=[\'"][^\'"]*[\'"]\s+height=[\'"][^\'"]*[\'"]',
      r'\1 width="100%" height="100%"',
      svg,
      count=1,
  )
  return svg
```

4.  **Plot & Chart Sizing Invariants for HTML Viewers**: Many skill scripts
    export 300-DPI publication figures (e.g., `4800 x 3600` pixels). Embedding
    those images or wide SVGs into an HTML viewer without strict CSS bounds
    renders oversized plots that overflow the viewport. Always apply these CSS
    sizing rules in custom HTML visualizers:
    -   Bound every embedded plot (`<img>`, `<svg>`, `<canvas>`) with
        `max-width: 100%; max-height: min(360px, 42vh); width: auto; height:
        auto; object-fit: contain; display: block; margin: 0 auto;`.
    -   Arrange multiple plots in a responsive CSS grid (`grid-template-columns:
        repeat(auto-fit, minmax(320px, 1fr)); gap: 16px;`) with click-to-zoom
        lightbox support rather than stacking full-width figures vertically.
    -   When generating Matplotlib figures specifically for inline HTML
        dashboards, prefer compact web dimensions (`figsize=(6, 4)`, `dpi=144`),
        reserving `dpi=300` for standalone print/report artifacts.

```css
.chem-svg-well {
  background: #ffffff;
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 8px;
  display: flex;
  align-items: center;
  justify-content: center;
}
.plot-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
  gap: 16px;
}
.plot-card {
  background: var(--card);
  border: 1px solid var(--border);
  border-radius: 10px;
  padding: 12px;
  overflow: hidden;
}
.plot-card img,
.plot-card svg,
.plot-card canvas {
  max-width: 100%;
  max-height: min(360px, 42vh);
  width: auto;
  height: auto;
  object-fit: contain;
  display: block;
  margin: 0 auto;
}
```

--------------------------------------------------------------------------------

## 5. Antigravity (AGY) Host Theme Variables & Card Layout Rules

When building a custom HTML application beyond the turnkey viewers, bind all
surfaces, text, borders, and interactive controls to the **host-injected
Antigravity (AGY) CSS variables** rather than hardcoding hex colors or declaring
high-specificity local `:root` color overrides:

1.  **Zero-Specificity Theme Fallbacks (`:where(:root)`) & No `/preload.js`**:
    -   **Never include `<script src="/preload.js"></script>` in HTML
        artifacts**: `/preload.js` is only for UI plugin sidecars. In HTML
        artifacts rendered inside sandboxed `SafeContentFrame`
        (`*.scf.usercontent.goog`) iframes, the injected CSP is `script-src
        'unsafe-inline'` (plus optional `<script
        src="https://www.gstatic.com/antigravity/web/dev/tailwindcss.min.js"></script>`),
        so `/preload.js` is blocked by CSP and 404s.
    -   The AGY host automatically injects its active `:root` theme variables in
        `<head>` and sets `class="dark"` or `class="light"` on `<html>`. Wrap
        any standalone fallback tokens in `:where(:root)` and
        `:where(html.dark)` (zero specificity) so they never clobber the
        host-injected `:root` variables:
    -   **Surfaces & Borders**: `var(--background)` (page background),
        `var(--content)` (main content area), `var(--card)` (elevated cards,
        modals, panels), `var(--sidebar)` (navigation/filter rails),
        `var(--border)` (dividers, card borders, input outlines).
    -   **Typography**: `var(--foreground)` (headings, primary values),
        `var(--muted-foreground)` (secondary labels, captions, axis ticks),
        `var(--placeholder)` (input placeholders).
    -   **Interactive Accents**: `var(--primary)` with
        `var(--primary-foreground)` (primary buttons, active tab pills),
        `var(--secondary)` with `var(--secondary-foreground)` (secondary chips
        and controls), `var(--accent)` (hover highlights, focus rings, selected
        compound cards).
2.  **Molecule Card Property & Metadata Formatting (No Right-Edge Smooshing)**:
    -   **Stacked Property Cells**: Inside 240–280px molecule cards, never use
        `display: flex; justify-content: space-between;` inside multi-column
        property cells or inline `<b>MW:</b> 416.151` text, which pushes numbers
        against the right border of each cell. Instead, stack `.prop-label`
        (uppercase `9.5px`, `var(--muted-foreground)`) above `.prop-val`
        (`11.5px` bold, `font-variant-numeric: tabular-nums`) and round
        `MW`/`TPSA` to 1 decimal place and `cLogP`/`QED`/`CNS_MPO` to 2 decimal
        places.
    -   **Humanized Metadata Keys**: Convert `snake_case` column names (e.g.,
        `binding_affinity_nM` to `binding affinity nM`) while preserving the
        caller's casing and keeping raw column names on hover tooltips, and
        align key-value rows with `justify-content: flex-start; gap: 6px;` and
        `max-width: 68%; text-overflow: ellipsis;` on `.meta-key` so long column
        names never push numeric values against the right edge of the card.
    -   **Atomic File Writes**: When writing `.html` and `.metadata.json` files
        from Python scripts, always write to a temporary file in the same
        directory and `os.replace()` atomically so the host file watcher never
        mounts a truncated or half-written HTML payload.
3.  **Data Visualization & Status Palette (`PUBLICATION_PALETTE`)**: For inline
    SVG/`<canvas>` scatter plots, histograms, radar charts, and medicinal
    chemistry status badges (Ro5, PAINS/Brenk alerts, CNS MPO), use the shared
    palette from `scripts/compat.py`:
    -   **Primary / Series 1**: `#4F46E5` (Indigo)
    -   **Favorable / Clean PAINS / Ro5 Pass**: `#10B981` (Emerald)
    -   **Intermediate / 1 Ro5 Violation**: `#F59E0B` (Amber)
    -   **Alert / Activity Cliff / PAINS Liability**: `#F43F5E` (Rose)
    -   **Additional Series**: `#8B5CF6` (Violet), `#06B6D4` (Cyan)

--------------------------------------------------------------------------------

## 6. Precomputing Rich Payloads for Multi-View Standalone Apps

When a user requests a custom multi-tab standalone HTML application (for
example, combining a 2D SAR compound browser and an activity-cliff network in
one deliverable):

1.  **Run the Skill CLI Scripts First to Generate Verified Data**:
    -   Curate and featurize the dataset with
        `../cheminformatics_core/scripts/descriptors.py` and
        `../cheminformatics_core/scripts/filter_catalogs.py`.
    -   Compute matched molecular pairs and SALI activity cliffs with
        `../cheminformatics_ml/scripts/mmpa_sar_analysis.py`.
    -   Generate 3D conformers or docked poses with
        `../cheminformatics_physics/scripts/generate_3d_conformers.py`.
2.  **Embed Precomputed JSON + Inline SVGs (Zero External CDNs)**:
    -   Serialize the curated compound records, inline responsive SVGs, and MMPA
        pair tables directly into the standalone HTML bundle (or generate a
        separate linked artifact viewer via
        `scripts/interactive_grid_viewer.py`).
    -   Register every user-facing `.html` application in the root of
        `<artifact_dir>/` (`UserFacing: True`) and link it using
        `[app_name.html](file://<artifact_dir>/app_name.html)`.

--------------------------------------------------------------------------------

## 7. Anti-Patterns to Avoid

-   **Never build a chemistry app without consulting the `cheminformatics-*`
    skills**: Do not write ad-hoc descriptor calculators, naive regex SMILES
    parsers, or uncalibrated MPO formulas when `cheminformatics-core`,
    `cheminformatics-physics`, `cheminformatics-ml`,
    `cheminformatics-reactions`, `cheminformatics-doc2mol`, and
    `cheminformatics-cofolding` provide tested implementations.
-   **Never load external JS CDNs in artifact HTML files**: External
    `<script src="https://...">` tags for `3Dmol.js`, `SmilesDrawer`,
    `RDKit.js`, `Plotly`, or `Chart.js` are blocked by CSP in sandboxed
    webviews. Always pre-render 2D SVGs in Python via RDKit and use the
    self-contained SVG architecture in `scripts/interactive_grid_viewer.py`.
-   **Never substitute ASCII art or Mermaid diagrams for chemical structures**:
    Every molecule or reaction displayed in an application UI must be rendered
    as an authentic 2D chemical graph SVG/PNG.
