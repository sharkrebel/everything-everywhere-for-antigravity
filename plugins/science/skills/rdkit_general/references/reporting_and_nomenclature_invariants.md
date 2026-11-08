# Scientific Reporting, Nomenclature, and Table Invariants

This reference document defines mandatory standards for technical reports,
chemical nomenclature verification, tabular formatting, and scientific tone
across the RDKit skill suite.

--------------------------------------------------------------------------------

## 1. Mandatory GitHub Flavored Markdown (GFM) Tables

-   **Strict Prohibition on ASCII Code Blocks**: Never render compound property
    tables, SAR summaries, candidate comparisons, or synthetic routes as raw
    monospace text inside fenced code blocks.
-   **Authentic GFM Table Formatting**: All tabular deliverables must use
    standard GitHub Flavored Markdown table syntax:

    ```markdown
    | Compound ID | Canonical SMILES | MW (Da) | cLogP | Predicted pKi | Macrocycle Ring Size | Chiral Centers |
    | :---------- | :--------------- | ------: | ----: | ------------: | :------------------: | :------------: |
    | DESIGN_01   | Cc1cc(O)cc(C)c1  |   406.4 |  2.85 |          8.42 |     14-membered      |      (R)       |
    | DESIGN_02   | Cc1cc(F)cc(C)c1  |   408.4 |  3.10 |          8.79 |     14-membered      |      (S)       |
    ```

-   **Consistent Column Alignment**: Left-align text and identifiers (`:---`),
    right-align numeric metrics and descriptors (`---:`), and center categorical
    flags or stereocenters (`:---:`).

--------------------------------------------------------------------------------

## 2. Authoritative Scientific Voice & Zero Automation Artifacts

-   **Peer-Reviewed Scientific Tone**: Reports must read as authoritative
    medicinal chemistry dossiers, publication manuscripts, or executive project
    reviews.
-   **Zero Tooling / Script Artifacts in Narrative**:
    -   Never mention internal python automation scripts (e.g.
        `document_mol_curator.py`, `congeneric_series_ml.py`,
        `markush_assembler.py`) in user-facing reports.
    -   Never mention test runners, exit codes, or automated assertion messages
        (e.g. `status: PASSED`, `score: 6/6`, `All 5 assertions succeeded`).
    -   Instead, report the scientific outcome directly:
    -   *Incorrect*: "We ran `document_mol_curator.py` and it passed all 6
        assertions. Then we ran `congeneric_series_ml.py`..."
    -   *Correct*: "Congeneric series curation yielded 42 stereochemically
        defined analogs with a 4.2 log-unit dynamic range in target binding
        affinity. Multi-model QSAR benchmarking identified ElasticNetCV with 2D
        physicochemical descriptors as the optimal predictive architecture
        ($R^2_{\text{CV}} = 0.74$, $\text{RMSE} = 0.41$)..."

--------------------------------------------------------------------------------

## 3. Programmatic Verification of All Naming & Topological Descriptors

-   **Zero Unverified Assertions**: Every systematic chemical name, IUPAC name,
    ring size, stereocenter designation, bridgehead description, and rotatable
    bond count must be programmatically verified against the RDKit molecular
    graph using `scripts/verify_chemical_nomenclature.py`.
-   **Ring Perception & Macrocycle Sizing**:
    -   Never state a molecule is an "$N$-membered ring" or "$N$-membered
        macrocycle" based on visual estimation or analogy.
    -   Run `analyze_molecular_topology(mol)` or `verify_ring_size(mol, N)` to
        extract the exact atom count of each ring via
        `mol.GetRingInfo().AtomRings()`.
    -   Every ring or macrocycle size cited in prose must match the exact atom
        count returned by `mol.GetRingInfo().AtomRings()`.
-   **Stereochemical Verification**:
    -   Every $(R)/(S)$ stereocenter and $(E)/(Z)$ olefin descriptor must be
        verified using `Chem.FindMolChiralCenters(mol, includeUnassigned=True)`
        and double-bond stereo flags.
    -   Undefined or racemic stereocenters must be explicitly designated as such
        rather than assigned arbitrary $(R)/(S)$ prefixes.
-   **IUPAC / Systematic Name Round-Trip Resolution**:
    -   For any proposed systematic name, verify round-trip consistency using
        `verify_iupac_name(mol, claimed_name)` which queries OPSIN / PubChem PUG
        REST.
    -   Ensure the standard InChIKey matches: $\text{InChIKey}(\text{mol}) ==
        \text{InChIKey}(\text{resolved name})$.

--------------------------------------------------------------------------------

## 4. Visual Embed Protocol

-   **Feature Attribution Plots**: Always embed feature importance / attribution
    as high-resolution plots (`![Feature Attribution](...)`), never as raw ASCII
    text tables.
-   **Reaction Mechanism Diagrams**: Whenever proposing a chemical synthesis,
    cyclization, or key bond-forming step, render an authentic 2D curved
    electron-pushing arrow diagram using `cheminformatics-reactions`
    (`../cheminformatics_reactions/scripts/render_mechanism.py`) and embed it in
    the report:

    ```markdown
    ![Key Bond-Forming Macrocyclization Mechanism](images/macrocyclization_mechanism.png)
    ```
-   **Validated Literature Precedents**: Every proposed reaction step, reagent,
    catalyst, and set of conditions must cite authenticated peer-reviewed
    literature or patent precedents (author, journal, year, volume, page/DOI or
    patent number).

### Linking to generated files and avoiding Mermaid chemical diagrams

A report that lists its own outputs must make each one clickable. Three distinct
forms apply, and they are not interchangeable:

-   **Interactive HTML viewers (`interactive_grid_viewer.py` outputs)** must be
    registered in the root of the conversation's artifact directory
    (`<artifact_dir>/<viewer>.html`) by first calling `write_to_file`
    (`UserFacing: True`) to create the placeholder artifact, then running
    `interactive_grid_viewer.py --output <artifact_dir>/<viewer>.html` (which
    writes the HTML and camelCase `.html.metadata.json` with `"userFacing":
    true`), and linked strictly to that root artifact path:

    ```markdown
    [interactive_compound_explorer.html](file://<artifact_dir>/interactive_compound_explorer.html)
    ```

    **Never** link an `.html` viewer to a source code repository path or a
    nested subdirectory; clicking a non-artifact `.html` link opens raw HTML/SVG
    source code in the code editor instead of rendering the interactive webview.

-   **Data deliverables the reader opens** (CSV, SDF, CDXML, PDFs) use a plain
    link with an absolute `file://` URL:

    ```markdown
    [curated_corpus.csv](file:///abs/path/to/curated_corpus.csv)
    ```

-   **Images the reader sees inline** use the embed form with an absolute path
    inside the conversation's artifact directory and no scheme:

    ```markdown
    ![Property distribution](/abs/path/to/artifact_dir/property_distributions.png)
    ```

-   **No Mermaid Diagrams for Chemistry, SAR Progressions, or Compound
    Summaries**: **Never** use `mermaid` (`flowchart LR` / `flowchart TD`, etc.)
    to depict molecules, lead-optimization progressions (`Compound A ->
    Precursor B -> Clinical Candidate C`), reaction schemes, or multi-parameter
    compound cards. Mermaid draws text boxes instead of chemical graphs, and raw
    HTML tags (`<b>`, `<br/>`, `<sub>`, `<sup>`) inside Mermaid node/edge labels
    are unsupported by the artifact renderer and produce distorted boxes with
    raw markup. Always use RDKit 2D structure depictions
    (`scripts/draw_molecules.py` or `render_mechanism.py`) and GFM Markdown
    tables instead.

Never wrap links or image embeds in backticks. A backtick code span suppresses
all inline markup inside it, so `` `[name](file:///path)` `` renders as literal
text and the link is lost. Relative paths and bare filenames do not resolve for
the reader; always use the absolute path.
