---
name: biomolecular-structure-viewer
description: >
  Interactive 3D biomolecular structure viewer via Mol*. Renders 3D molecular structures
  from local workspace files (.cif, .mmcif, .pdb) or conversation artifacts,
  and manages user 3D annotations and highlighted residues.
  Use when the user asks to visualise a protein or nucleic acid, display a 3D structure,
  or asks about 3D biomolecular structure annotations and highlighted residues.
metadata:
  publisher: google
---

# Biomolecular Structure Viewer

The Biomolecular Structure Viewer sidecar renders interactive 3D molecular
structures in the Antigravity side pane using Mol*.

All interaction with the viewer is through the
**`science_biomolecular_structure_viewer_mcp`** MCP tool catalog. Never
open or edit the viewer's own code, config or saved state; the MCP tools are the
only way to operate it. They validate inputs, apply atomic state patches, and
prevent race conditions.

> [!IMPORTANT] **Local files only**: The viewer only loads structure files saved
> locally. Always download the structure into the current conversation's
> artifact directory first (using `pdb-database` or
> `alphafold-database-fetch-and-analyze`, see
> [Where files live](#where-files-live)), then pass the file path to
> `load_structure`. Database IDs (e.g. `1CRN`, `P01308`, `AF-P01308-F1`) and
> URLs are rejected; download URLs into the artifact directory first with
> `curl`, not a web-page fetch tool such as `read_url_content`, which can
> silently truncate large files.

--------------------------------------------------------------------------------

## Where files live

The viewer is scoped to the **current conversation**: it only lists and loads
structure files from that conversation's artifact directory.

Viewer state (including user notes) is managed entirely by the MCP tools; never
read or write it directly.

--------------------------------------------------------------------------------

## Standard Biomolecular Structure Visualization Workflow

When the user asks to **visualise**, **show**, or **display** a biomolecular
structure:

1.  **Obtain Local Structure File (Required)**:
    -   If the structure file (`.cif`, `.mmcif`, or `.pdb`) is already in the
        conversation's artifact directory (including `.user_uploaded/`), use its
        path relative to that directory.
    -   If the file is elsewhere on disk (e.g. the user's workspace), copy it
        into the conversation's artifact directory first (e.g. `cp`): the viewer
        does not load files from outside it.
    -   If given a protein name, gene symbol, UniProt ID, or PDB ID, use
        appropriate search/fetch skills or tools (e.g. `pdb-database`,
        `alphafold-database-fetch-and-analyze`) to download the 3D coordinate
        file into the conversation's artifact directory. Decompress `.gz` files
        (e.g. `gunzip`); the viewer does not read gzipped files.
    -   If the download fails, tell the user why and stop. Do not pass the
        database ID to `load_structure`: it is rejected.
2.  **Load into Viewer**: Call the **`load_structure`** MCP tool with the local
    file path relative to the conversation directory (e.g.
    `P00533/AF-P00533-F1-model_v6.cif` or `4HHB/pdb_00004hhb.cif`).
3.  **Surface Sidecar Link**: Always include this link in your chat response so
    the user can open the viewer with one click:

    ```markdown
    [Biomolecular Structure Viewer](sidecar://science/biomolecular_structure_viewer/ui.html)
    ```

--------------------------------------------------------------------------------

## MCP Tool Catalog (`science_biomolecular_structure_viewer_mcp`)



| Tool                 | Parameters                | Description               |
| -------------------- | ------------------------- | ------------------------- |
| `get_state`          | *(none)*                  | Returns the sidecar's     |
:                      :                           : current state (active     :
:                      :                           : structure, `selections`,  :
:                      :                           : the user's                :
:                      :                           : `activeSelection`) and    :
:                      :                           : state version.            :
| `load_structure`     | `id` (required), `source` | Loads a local structure   |
:                      : (optional\: `'workspace' | : file into the viewer:     :
:                      : | 'artifact'`)            : a workspace/artifact file :
:                      :                           : path or artifact          :
:                      :                           : filename. Database IDs    :
:                      :                           : and URLs are rejected     :
:                      :                           : (download first).         :
:                      :                           : Each structure keeps its  :
:                      :                           : own selections, style and :
:                      :                           : overlays (see             :
:                      :                           : **Per-structure view      :
:                      :                           : state**): switching back  :
:                      :                           : restores them, and        :
:                      :                           : loading the structure     :
:                      :                           : already shown changes     :
:                      :                           : nothing. `source` is      :
:                      :                           : auto-detected when        :
:                      :                           : omitted.                  :
| `set_representation` | `style` (required:        | Changes the 3D rendering  |
:                      : `'cartoon' | 'stick' |    : representation and        :
:                      : 'sphere' | 'surface' |    : optional molecular        :
:                      : 'line'`), `components`    : component visibility      :
:                      : (optional\: `{waters,     : flags.                    :
:                      : sugars, ligands,          :                           :
:                      : sidechains}`)             :                           :
| `set_overlays`       | `annotations`,            | Shows or hides the        |
:                      : `confidence` (optional    : annotation pins / notes   :
:                      : booleans; pass at least   : overlay and the pLDDT     :
:                      : one)                      : confidence coloring and   :
:                      :                           : legend. See **Overlays**  :
:                      :                           : below.                    :
| `highlight_region`   | `residues` (see below) or | Appends ONE unlabelled    |
:                      : `start_auth_seq_id` /     : color entry (Mol\*        :
:                      : `end_auth_seq_id`;        : overpaint) covering every :
:                      : `auth_asym_id` (default   : listed residue to         :
:                      : `'A'`), `color` (named    : `selections`, without     :
:                      : color or `#rrggbb`)       : creating an annotation    :
:                      :                           : note. Returns its `id`.   :
:                      :                           : It is drawn on top of     :
:                      :                           : older entries, notes      :
:                      :                           : included, which keep      :
:                      :                           : their residues and        :
:                      :                           : color.                    :
| `annotate_region`    | `text` (required);        | Appends ONE labelled      |
:                      : `residues` (see below) or : entry (annotation note)   :
:                      : `start_auth_seq_id` /     : covering every listed     :
:                      : `end_auth_seq_id`;        : residue, possibly on      :
:                      : `auth_asym_id` (default   : several chains, to        :
:                      : `'A'`); `color`           : `selections` and turns    :
:                      : (optional; defaults to    : the `annotations` overlay :
:                      : the residue's color,      : on. Returns its `id` and  :
:                      : else `#ff4fa8`)           : the resolved `color`. If  :
:                      :                           : a note already covers     :
:                      :                           : exactly these residues,   :
:                      :                           : it is updated in place    :
:                      :                           : instead (`updated` is     :
:                      :                           : true; same `id`).         :
| `remove_annotation`  | `id` (required, e.g.      | Removes the annotation    |
:                      : `'s2'`)                   : note with that `id`.      :
:                      :                           : Errors (listing valid     :
:                      :                           : ids) if it is not a       :
:                      :                           : labelled entry.           :
| `clear_annotations`  | *(none)*                  | Removes all labelled      |
:                      :                           : entries (notes). Color    :
:                      :                           : highlights stay.          :
| `clear_highlights`   | *(none)*                  | Removes all unlabelled    |
:                      :                           : entries (color            :
:                      :                           : highlights). Notes stay.  :



**`residues`** (both `highlight_region` and `annotate_region`) is a list of
residue numbers and/or `{"auth_asym_id", "start_auth_seq_id",
"end_auth_seq_id"}` objects, in any mix. A number is one residue on
`auth_asym_id`: `[6, 10]` means residues 6 and 10, not the span 6–10. Use an
object only for a span (`{"start_auth_seq_id": 57, "end_auth_seq_id": 60}`) or a
residue on another chain; `end_auth_seq_id` defaults to `start_auth_seq_id` and
a missing `auth_asym_id` defaults to the top-level `auth_asym_id`.
`start_auth_seq_id`/`end_auth_seq_id` is shorthand for a single span. Entries of
`activeSelection` (`{"startResidue", "endResidue"}`, see
[Current User Selection](#current-user-selection-stateactiveselection)) are
accepted too, as they are. Only this form carries residue names: when the user
gives them, or you read them from the structure file
(`_atom_site.auth_comp_id`), pass them as `auth_comp_id`, e.g. `{"startResidue":
{"auth_seq_id": 57, "auth_comp_id": "HIS"}, "endResidue": {"auth_seq_id": 60,
"auth_comp_id": "GLY"}}`. They are stored as given (not checked against the
structure).

**`chains`** (both tools) lists **whole chains** by `auth_asym_id`, e.g.
`chains=["B"]`. Use it whenever the user means an entire chain ("color chain B
red", "note on chain C") rather than listing the chain's residues; it can be
combined with `residues` on other chains. The viewer fills in the chain's
residues from the loaded structure and shows it as `Chain B`.

When the same note applies to several residues, make **one** `annotate_region`
call listing all of them — do not create one note per residue. Likewise, one
`highlight_region` call colors all listed residues.

**Choose by whether the request includes text, not by the verb.** Users may say
"highlight", "color", "annotate", "note", "label" or anything else:

*   The request attaches **text** to residues (a label, note, name or comment),
    even when it also asks for a color ("highlight the catalytic triad in blue
    and label it"): make **one** `annotate_region` call, passing the color as
    `color`. An annotation colors its residues itself.
*   The request asks for **color only**, even when it says what the residues are
    ("color the binding-site residues green"): call `highlight_region`. Don't
    invent a label the user did not ask for.
*   Never call both on the same residues: that adds a second, unlabelled entry,
    which the viewer lists as a separate "No annotation" card.

Calling `annotate_region` on exactly the residues of an existing note (e.g. one
the user added in the viewer) **overwrites** that note's text, keeping its `id`,
position and color unless you pass `color`. So when the user asks to add
information to a note they already made, call `annotate_region` on the same
residues with the full new text (include anything worth keeping from the old
note, returned as `previousNote`). A note on different or only partly
overlapping residues is added alongside.

### Examples

**Load a downloaded AlphaFold model:**
`load_structure(id="P00533/AF-P00533-F1-model_v6.cif")`

**Load a downloaded PDB entry:** `pdb-database`'s `download_coordinate_files.py
--ids 4HHB --format mmcif --output_dir <conv artifact dir>/4HHB` saves
`4HHB/pdb_00004hhb.cif.gz`. The viewer does not read gzipped files, so `gunzip`
it first, then `load_structure(id="4HHB/pdb_00004hhb.cif")`.

**Load from a conversation artifact:** `load_structure(id="prediction.cif",
source="artifact")`

**Rejected:** `load_structure(id="1CRN")`,
`load_structure(id="https://files.rcsb.org/download/1CRN.cif")`. Download the
file first.

**Change representation and hide waters:** `set_representation(style="stick",
components={"waters": false})`

**Show annotation pins and hide pLDDT coloring:**
`set_overlays(annotations=true, confidence=false)`

**Highlight an active site in blue:** `highlight_region(start_auth_seq_id=50,
end_auth_seq_id=60, auth_asym_id="A", color="#0000ff")`

**One note on two separate residues:** `annotate_region(text="hi 2",
auth_asym_id="A", residues=[6, 10])`

**Annotate a binding pocket:** `annotate_region(start_auth_seq_id=100,
end_auth_seq_id=115, text="ATP binding pocket", auth_asym_id="A")`

**Annotate a site spanning two chains:**

```
annotate_region(residues=[{"auth_asym_id": "A", "start_auth_seq_id": 57,
                           "end_auth_seq_id": 60},
                          {"auth_asym_id": "B", "start_auth_seq_id": 102}],
                text="Catalytic triad")
```

**Remove that note:** `remove_annotation(id="s3")` (use the `id` returned by
`annotate_region` or listed in `get_state`)

**Update or recolor a note:** call `annotate_region` again with the same
`residues`, the new `text`, and optionally the new `color`; the existing note is
overwritten in place. (The user can also rename a note from the viewer's side
panel, or recolor it by saving it again on the same residues.)

--------------------------------------------------------------------------------

## Rendering Styles

The `style` parameter on `load_structure` and `set_representation` controls the
3D rendering representation:

Value     | Description
--------- | -----------------------------------------
`cartoon` | Ribbon / helix / sheet backbone (default)
`stick`   | Bond sticks with small atom spheres
`sphere`  | Van der Waals spacefill spheres
`surface` | Gaussian molecular surface
`line`    | Thin lines connecting bonded atoms

## Molecular Components

The optional `components` parameter on `set_representation` controls visibility
of specific molecular component categories. All sub-properties are boolean; a
structure shown for the first time has only `ligands` on:

Sub-property | Description
------------ | -------------------------------------
`waters`     | Solvent water molecules
`sugars`     | Carbohydrate / glycan chains
`ligands`    | Non-polymer small-molecule ligands
`sidechains` | Amino acid sidechains (Cβ and beyond)

## Overlays (`annotations` and `confidence`)

Use **`set_overlays`** to switch the two viewport overlays on or off. Pass only
the overlays you want to change:

Overlay       | Default | Description
------------- | ------- | -----------------------------------------
`annotations` | `false` | 3D annotation pins and ribbon highlights.
`confidence`  | `true`  | Per-residue pLDDT confidence 3D coloring.

```
set_overlays(annotations=true)
set_overlays(confidence=false)
```

Both flags are kept per structure (see
[Per-structure view state](#per-structure-view-state-statestructures)): a
structure shown for the first time starts from the defaults, and switching back
to one restores its flags. `annotate_region` turns `annotations` on by itself.
The viewer does not turn it on just because a structure has saved notes: call
`set_overlays(annotations=true)` when the user wants to see existing notes. The
user can also show or hide each overlay with the eye button in the header of the
**Selections** and **Confidence** sections of the viewer's side panel (opened
with the button at the top right of the viewport); don't call `set_overlays`
just to preserve their choice.

Predicted models that carry per-residue confidence scores (AlphaFold DB,
ColabFold, ESMFold, and similar) are colored by pLDDT; the side panel's
**Confidence** section shows the legend and the share of residues in each band.
Detection is based on the contents of the file, not its name, so `confidence` is
a no-op on an experimental structure — it keeps its default coloring and no
Confidence section appears.

Call `set_overlays(confidence=false)` when the user wants to see a predicted
model in the standard color scheme, for example to compare it side by side with
an experimental structure.

--------------------------------------------------------------------------------

## Per-structure view state (`state.structures`)

Each structure the viewer has shown keeps its own view state in
**`structures`**, keyed by the structure's workspace path (or artifact
filename): its `style`, `components`, `selections`, `annotations` and
`confidence`. The structure shown is the one named by the top-level `source` and
`path`/`filename`; `nextSelectionId` and the user's `activeSelection` stay at
the top level:

```json
{"source": "workspace", "path": "4HHB/4hhb.cif",
 "nextSelectionId": 4, "activeSelection": null,
 "structures": {
   "4HHB/4hhb.cif": {"style": "stick", "annotations": true,
                     "selections": {"s3": {"label": "Heme pocket", ...}}},
   "1crn.cif": {"confidence": false,
                "selections": {"s1": {"color": "#ff1f3d", ...}}}}}
```

*   Switching structures (`load_structure`, the viewer's file list or an upload)
    only changes which structure is shown; switching back restores its
    selections, style, components and overlays. A structure with no entry is
    shown with the defaults: `cartoon`, only ligands visible, `annotations` off,
    `confidence` on.
*   `set_representation`, `set_overlays` and the selection tools act on the
    structure shown only. To change another structure, load it first.
*   An entry holding nothing but defaults is removed.
*   Those tools need a structure to be shown: until one is, they return an
    error asking you to call `load_structure` first.

## Selections (`state.structures[...].selections`)

Color highlights and annotation notes share one id-keyed map, **`selections`**,
in the shown structure's entry of `structures`. Each entry is stored under its
selection id (`s<n>`):

```json
{"s3": {"label": "Catalytic triad", "color": "#06b6d4",
        "residues": [
          {"startResidue": {"auth_seq_id": 57, "auth_comp_id": "HIS",
                            "auth_asym_id": "A"},
           "endResidue": {"auth_seq_id": 60, "auth_comp_id": "GLY",
                          "auth_asym_id": "A"}},
          {"startResidue": {"auth_seq_id": 102, "auth_comp_id": "ASP",
                            "auth_asym_id": "B"},
           "endResidue": {"auth_seq_id": 102, "auth_comp_id": "ASP",
                          "auth_asym_id": "B"}}]}}
```

`residues` uses the same residue spans as `activeSelection` (see below): runs of
residues on one chain, from `startResidue` to `endResidue`. `auth_comp_id` is
the name passed to the tools (stored as given), replaced by the loaded
structure's name whenever the viewer saves the entry; residues named by neither
have none.

*   **Unlabelled** entry (`color`, no `label`): a color highlight.
*   **Labelled** entry: an annotation note over every range in `residues`. Its
    pin sits on the first range; ribbons and click-to-zoom cover all of them. A
    labelled entry always has a `color`: its ribbon, 3D pin, and notes-panel
    marker are drawn in it, and the user can change it in the viewer. When a
    note is created without one, it takes the color of the last color highlight
    covering its first residue (the first residue of its first range, after
    sorting), or `#ff4fa8` (hot pink) if that residue is not colored. Entries
    read without a `color` get one by the same rule.
*   Entries are painted in map order: where a note and a color highlight
    overlap, whichever comes later in `selections` is on top.
*   The key (`s<n>`) is the selection's stable `id` and is what
    `remove_annotation` takes. Spans are sorted by chain and start, and
    overlapping or contiguous spans are merged.
*   Ids are unique across structures and never reused, even after an entry is
    removed. The next free number is kept in the top-level `nextSelectionId`;
    don't edit it.
*   An entry the user saved with whole chains selected (every residue of the
    chain, however it was picked), or that `highlight_region`/ `annotate_region`
    wrote with `chains`, also has `"chains": ["A", ...]`: the chains selected
    whole, which the viewer lists as `Chain A` instead of their residue ranges
    (residues selected on other chains are listed one by one, in the same
    entry). `residues` lists every span of those chains too, so read the
    residues from `residues`. Right after a tool call a chain may have no spans
    in `residues` yet: the viewer fills in the chain's polymer residues from the
    loaded structure, or all of its residues for a chain with no polymer (e.g.
    only ligands), and drops chains the structure does not have. The viewer
    drops a chain from `chains` once its spans change.

## Current User Selection (`state.activeSelection`)

The residues the user has selected right now (clicked in 3D or the sequence
panel, not yet saved) are reported by the viewer as **`activeSelection`**, a
list of runs of residues on one chain, each from `startResidue` to `endResidue`
(the same residue for a single one). A residue is `{"auth_seq_id",
"auth_comp_id", "auth_asym_id"}` (the mmCIF author residue number, residue name
and chain), with `auth_comp_id` read from the loaded structure (left out when
unknown); `[]` when nothing is selected:

```json
{"activeSelection": [
  {"startResidue": {"auth_seq_id": 1, "auth_comp_id": "ADP",
                    "auth_asym_id": "A"},
   "endResidue": {"auth_seq_id": 1, "auth_comp_id": "ADP",
                  "auth_asym_id": "A"}},
  {"startResidue": {"auth_seq_id": 308, "auth_comp_id": "THR",
                    "auth_asym_id": "A"},
   "endResidue": {"auth_seq_id": 310, "auth_comp_id": "LYS",
                  "auth_asym_id": "A"}}]}
```

When the user refers to "this residue", "the selected residues", or "what I
clicked", call `get_state` and use `activeSelection`, e.g. pass its entries as
they are as `residues` to `annotate_region` or `highlight_region`. If it is
empty, ask the user which residues they mean. It is read-only and transient: the
viewer overwrites it on every selection change, so don't write it, and don't
treat it as a saved highlight or note.

> [!IMPORTANT] **Author numbering**: all ids are author ids (`auth_asym_id`,
> `auth_seq_id`, `auth_comp_id`), never `label_*`. In 2GXA, author chain A
> residue 1 is ADP; label chain A residue 1 is a DT of author chain M.

--------------------------------------------------------------------------------

## Coloring Specific Residues

Use **`highlight_region`** to color specific residues, functional sites, active
pockets, or mutation sites **without creating annotation labels**:

```
highlight_region(residues=[50, 51, 52, 53], auth_asym_id="A",
                   color="#0000ff")
highlight_region(start_auth_seq_id=100, end_auth_seq_id=110,
                   auth_asym_id="B", color="#ff0000")
highlight_region(chains=["B"], color="#ff0000")
```

Use highlights when the user only wants color, with no text. These do not create
annotation labels or notes; when the request includes text, call only
`annotate_region` (with `color`), which colors the residues too. `auth_asym_id`
is optional and defaults to `A`. To color an **entire chain**, pass `chains`
instead of guessing its residue range: the viewer fills in the chain's residues
from the loaded structure and lists it as `Chain B`.

`#ff4fa8` is reserved for annotation ribbons. Avoid it, or your highlights will
be indistinguishable from annotation markers.

`highlight_region` stacks: the new color entry is drawn on top of older entries
on the same residues, notes included, which keep their residues and color, so
deleting it (the user can, in the Selections list) shows them again. Each call
adds an entry, so don't recolor residues just to restate a color they already
have.

Users color residues themselves from the viewer's side panel. A click in 3D or
in the sequence panel selects a residue (or a whole chain, with **Select:
Chain**), Shift+click adds the range from the last clicked residue to the
clicked one and zooms to everything selected (in ball-and-stick), Ctrl /
Cmd+click adds the clicked residue or removes it and zooms to everything
selected too, and a click on empty space or Escape clears it. The **Select &
customize** section lists what is selected as chips (a residue each, or `Chain
X` for a chain whose residues are all selected, whichever mode picked them; each
with a button that removes it from the selection), plus **Clear**. Below the
chips the user picks a color (swatches of the color names above and the note
pink, or a custom color from a picker or `#rrggbb`; the first swatch is
preselected), may type an annotation, and presses the single **Save**. Save
without an annotation paints the selection; with one it saves the note on those
residues, drawn in the chosen color (see **User 3D Annotations & Context
Notes**). These colors follow the same rules as `highlight_region` and are saved
to `selections` as unlabelled entries, so the agent and the user share them:
`get_state` returns what the user painted, and `clear_highlights` removes the
user's colors too.

To clear all highlights, call **`clear_highlights`**.

Highlights are scoped to the structure they were applied to: each structure
keeps its own in its `structures` entry. Switching structures shows only the new
structure's highlights, so a range like "chain A, 50-60" never colors unrelated
residues of another structure, and switching back brings them back.

--------------------------------------------------------------------------------

## User 3D Annotations & Context Notes

Annotation notes are the labelled entries of **`selections`**, shared by the
agent and the user. Users write notes in the side panel's **Select & customize**
section: select residues (in 3D or the sequence, across ranges and chains), pick
a color, type in the **Annotation** field and press Enter or Save. A note saved
on exactly the residues of an existing note replaces it, as with
`annotate_region`. Every entry of `selections`, notes and color highlights, is
listed in the **Selections** section below; a color highlight's card reads "No
annotation". Each card has a pencil to write the note in place (Enter or the
check saves, Escape cancels; giving a color highlight a note makes it a note), a
trash icon to delete the entry, and a chip per residue range; clicking a card, a
chip or a pin zooms to it. **Delete All** removes every entry. Notes the user
writes and notes you add with `annotate_region` end up in the same list, and
`get_state` always returns the current list. No annotation file is written to
disk.



| Task                   | Tool                                   |
| ---------------------- | -------------------------------------- |
| Add a note             | `annotate_region(text,                 |
:                        : start_auth_seq_id, end_auth_seq_id,    :
:                        : auth_asym_id)` or                      :
:                        : `annotate_region(text, residues)`,     :
:                        : which also turns the `annotations`     :
:                        : overlay on                             :
| Read the notes         | `get_state`, then use the `selections` |
:                        : entries that have a `label`            :
| Update a note          | `annotate_region` on the same residues |
:                        : with the new `text` (and `color`)      :
| Remove one note        | `remove_annotation(id)`                |
| Remove every note      | `clear_annotations()`                  |
| Show or hide the notes | `set_overlays(annotations=true/false)` |



Notes belong to the loaded structure, and each structure keeps its own list.
`load_structure` and switching structures in the viewer show the new structure's
notes (none the first time it is shown); switching back restores the old
structure's notes.

A note added with `annotate_region` and no `color` takes the color its residues
are painted in (pink otherwise); a note saved in the viewer takes the color
chosen next to it. Saving a color without a note from the viewer's **Select &
customize** section (see **Coloring Specific Residues**) leaves notes alone: the
color is drawn on top of any note on those residues. Users change a note's color
by saving the note again on the same residues.

### Agent Instructions for Handling Annotations

When the user asks about 3D biomolecular structure annotations, highlighted
residues, structure notes, or says *"tell me about the regions I've annotated"*:

1.  **READ ACTIVE STATE**: Call `get_state`. Check which structure is displayed
    and read the `selections` of its `structures` entry (labelled entries are
    notes, unlabelled ones are color highlights).
2.  **STRICT FIDELITY RULE**: Comment ONLY on the notes and highlights returned
    by `get_state`. Do NOT guess, hallucinate, or recall residues from previous
    turns or unselected regions.
3.  **SINGLE SOURCE OF TRUTH**: Before each answer, call `get_state` again,
    because the user can add or delete notes in the viewer at any time.
4.  **Synthesize Structural Analysis**: Provide detailed per-residue structural
    analysis (secondary structure location, hydrophobic packing, hydrogen
    bonding, disulfide bonds, domain boundaries, and mutation impacts) strictly
    for active annotations.

--------------------------------------------------------------------------------

## Key Principles

1.  **MCP tools are the only interface**: Always use the tools in the MCP Tool
    Catalog above (`load_structure`, `set_representation`, `set_overlays`,
    `highlight_region`, `clear_highlights`, `annotate_region`,
    `remove_annotation`, `clear_annotations`, `get_state`). Never open or edit
    the viewer's own code, config or saved state: nothing in them helps operate
    the viewer.
2.  **Scope files to the current conversation**: Store downloaded structure
    files as conversation artifacts.
3.  **Check active state**: Call `get_state` to verify what the viewer is
    actively displaying before answering questions about the visualization.
4.  **Delegate data retrieval**: Use `pdb-database` or
    `alphafold-database-fetch-and-analyze` when files need to be searched,
    downloaded, or analyzed on disk.
5.  **Keep visualizations in the sidecar**: Avoid dumping heavy raw coordinate
    text or ASCII drawings into chat.
6.  **Fix failed loads**: When the viewer can't show the state, you receive the
    notification `BiomolecularStructureViewer: <message>`. Check the file and
    that `id` is its path, fix the arguments or the file (download a bad file
    again with `curl`), and call the tool again. Save a fixed file under a new
    name: the same `id` does not reload. If both are fine, tell the user the
    viewer couldn't load it and stop.
