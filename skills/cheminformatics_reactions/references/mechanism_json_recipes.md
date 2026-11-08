# Mechanism Specification Format and Worked Recipes

Read this when a mechanism needs more than a single nucleophilic attack: several
arrows in one step, an intermolecular pair of reactants, a multi-step sequence,
or manual control over layout. `SKILL.md` carries the minimal single-arrow case;
everything below is the full surface of the same format.

Specs are read by `scripts/render_mechanism.py` from either `--rxn-json` or
`--rxn-yaml`. JSON and YAML use identical keys.

--------------------------------------------------------------------------------

## 1. Top-Level Keys

`id` (required)
:   Identifier for the step. The rendered PNG is written to
    `<out-dir>/<id>.png`, so it doubles as the output filename stem.

`title`
:   Caption printed above the scheme.

`r1`, `r2`, `r3`, ... (or `reactants`)
:   Reactant entries (`r1`, `r2`, `r3`, ... or a `reactants` list/dict). `r1` is
    required; additional reactant keys draw a `+` between consecutive reactants.

`p1`, `p2`, `p3`, ... (or `products`)
:   Product entries (`p1`, `p2`, `p3`, ... or a `products` list/dict). `p1` is
    required; additional product keys draw discrete co-products or leaving
    groups separated by `+`.

`cond`
:   Reaction-arrow annotations, `{"top": ..., "bottom": ...}`. Both accept raw
    LaTeX, e.g. `r"$\text{THF}, -78^\circ\text{C}$"`.

`arrows`
:   List of curved-arrow definitions. See section 2.

`add_h`
:   Maps a molecule key to an atom index (or `":N"` map label, or a list of
    indices/labels) whose implicit hydrogens should be drawn explicitly, e.g.
    `{"r1": 7}` or `{"r1": [":1", ":3"]}`. Needed when an arrow starts or ends
    on a C-H bond.

`rot`
:   Per-molecule rotation in degrees, e.g. `{"r1": 90}`. Any molecule omitted
    from `rot` is drawn at $0^\circ$ unless `--auto-rotate` or `"auto_rotate":
    true` is enabled, in which case `scripts/rotation_solver.py` automatically
    solves ChemDraw/IUPAC reactant orientations and MCS product alignments.

`auto_rotate`
:   Boolean (`true`/`false`). When `true` (or when `--auto-rotate` is passed on
    the CLI), automatically invokes `rotation_solver.solve_reaction_rotations`
    for any molecule not explicitly present in `rot`.

`steps`
:   List of elementary reaction step objects for rendering a multi-step cascade
    scheme (`compile_multistep_cascade`). Each step is verified with
    `ElectronGraph`, intermediate continuity between Step $k$ products and Step
    $k+1$ reactants is checked automatically, and a composite vertical scheme
    PNG is written to `<out-dir>/<id>.png`.

`catalytic_cycle` (or `"mode": "catalytic_cycle"`)
:   Renders a circular organometallic/enzymatic catalytic cycle
    (`compile_catalytic_cycle`) from `"intermediates": [{"id": "I1", "smi": ...,
    "label": ...}, ...]` and `"transitions": [{"from": "I1", "to": "I2",
    "label": ..., "in_label": ..., "out_label": ...}, ...]`.

`coord_map`
:   Per-molecule template coordinates for cases where the depictor's layout is
    still wrong. A last resort.

### Molecule Entries

Each of `r1`, `r2`, ..., `p1`, `p2`, ... is an object:

`smi` (required)
:   SMILES for the species. Atoms that arrows refer to are easiest to write with
    map labels, `[O:1]`, `[C:2]`. Map numbers are used to resolve arrow
    endpoints and are automatically stripped from the rendered 2D depiction so
    they do not clutter atom labels.

`name`
:   Label drawn beneath the structure.

`w`
:   Drawing width in cm. Defaults are 4.0 for `r1`, 3.0 for `r2`+, 4.0 for `p1`
    and 2.2 for `p2`+. Widen a crowded species rather than shrinking its
    neighbours.

`atom_labels`
:   Optional mapping from atom index or `":N"` map label to custom display text
    (e.g. `{"0": "R1", ":2": "Ar"}`), passed to RDKit's `opts.atomLabels` for
    generic R-group or protecting-group depictions.

--------------------------------------------------------------------------------

## 2. Arrow Types & Formal Pre-Render Verification

Every atom reference is either a zero-based atom index or an atom-map label such
as `":1"`, matching an atom written `[X:1]` in that molecule's SMILES. Map
labels survive edits to the SMILES string; indices do not, so prefer labels.

By default, every arrow moves an electron pair (`"electrons": 2`, rendered with
a full double-barbed `chemarrow` arrowhead). Setting `"electrons": 1` on any
arrow models a single-electron radical shift and renders a single-barbed harpoon
fishhook arrow (`fishhookarrow`).

Before rendering, `render_mechanism.py` automatically runs
`verify_mechanism_with_electron_graph(rxn_def, spec)` to enforce valence
electron conservation, formal charge balance, octet/non-hypervalency rules,
spectator locality (atoms not touched by arrows cannot silently gain/lose
hydrogens or charge), stereochemical Tiers 1–3, and product matching.

`attack`
:   Lone pair, radical, or bond attacking an atom. Takes `src` and `tgt`, each a
    `[molecule_key, atom]` pair (or `[molecule_key, [atom1, atom2]]` when a bond
    is the electron source). The two may name different molecules
    (intermolecular) or the same one (cyclization).

`bond_cleave`
:   Heterolysis (or homolysis with `"electrons": 1`) of a bond, electrons
    collapsing onto one atom. Takes `mol` and `bond` as `[atom, atom]`, plus an
    optional `tgt` naming the atom that keeps the electrons. `tgt` defaults to
    the second atom of `bond`.

`bond_to_bond`
:   A pi or sigma bond shifting to become another bond, as in enolization or a
    pericyclic step. Takes `mol`, `b1` (source bond) and `b2` (destination
    bond), each `[atom, atom]`.

`bond_to_atom`
:   A bond's electrons moving onto an atom in a different species. Takes
    `src_mol`, `bond`, `tgt_mol`, `tgt_atom`, and an optional `lift` (default
    0.8) controlling how far the arc bows away from the structures.

`lone_pair_to_bond`
:   A lone pair (or radical with `"electrons": 1`) forming a pi bond, as in
    carbonyl or oxocarbenium formation. Takes `mol`, `atom` and `bond`.

--------------------------------------------------------------------------------

## 3. Intermolecular Step

Two reactants (`r1` + `r2`), concerted nucleophilic attack and leaving-group
expulsion, drawn with a `+` between the reactants.

```json
{
  "id": "sn2_attack",
  "title": "Nucleophilic Substitution",
  "r1": {"smi": "C[CH2:1][Br:2]"},
  "r2": {"smi": "[I-:3]"},
  "p1": {"smi": "CCI"},
  "cond": {"top": "acetone", "bottom": "25 °C, - Br-"},
  "arrows": [
    {"type": "attack", "src": ["r2", ":3"], "tgt": ["r1", ":1"]},
    {"type": "bond_cleave", "mol": "r1", "bond": [":1", ":2"], "tgt": ":2"}
  ]
}
```

--------------------------------------------------------------------------------

## 4. Multi-Step Sequence

Each step is its own spec file and its own PNG. The product of step $k$ is the
reactant of step $k+1$: the second step below opens with the chlorosulfite
intermediate / acyl chloride, not with the carboxylate it came from. A sequence
that re-attacks the original starting material is depicting two competing
reactions, not one route.

```json
// scratch/step1.json — activation of a carboxylate with SOCl2
{
  "id": "step1_activation",
  "title": "Step 1: Carboxylate Attack on Thionyl Chloride",
  "r1": {"smi": "CC(=O)[O-:1]"},
  "r2": {"smi": "O=[S:2]([Cl:3])Cl"},
  "p1": {"smi": "CC(=O)OS(=O)Cl"},
  "cond": {"top": "CH2Cl2", "bottom": "- Cl-"},
  "arrows": [
    {"type": "attack", "src": ["r1", ":1"], "tgt": ["r2", ":2"]},
    {"type": "bond_cleave", "mol": "r2", "bond": [":2", ":3"], "tgt": ":3"}
  ]
}
```

```json
// scratch/step2.json — nucleophilic acyl substitution on the acyl chloride
{
  "id": "step2_substitution",
  "title": "Step 2: Alkoxide Acyl Substitution",
  "r1": {"smi": "C[C:1](=[O:2])[Cl:3]"},
  "r2": {"smi": "[O-:4]CC"},
  "p1": {"smi": "CC(=O)OCC"},
  "cond": {"top": "THF", "bottom": "- Cl-"},
  "arrows": [
    {"type": "attack", "src": ["r2", ":4"], "tgt": ["r1", ":1"]},
    {"type": "bond_cleave", "mol": "r1", "bond": [":1", ":3"], "tgt": ":3"}
  ]
}
```

Substitution and elimination steps carry two arrows, as step 2 does: the attack
on the reactive centre and the cleavage of the bond to the leaving group.
Drawing only the attack leaves the diagram with a pentavalent carbon.

Compile each step to a 300 DPI PNG:

```bash
uv run scripts/render_mechanism.py --rxn-json scratch/step1.json --out-dir ./images
uv run scripts/render_mechanism.py --rxn-json scratch/step2.json --out-dir ./images
```

Atom conservation is checked across the whole transformation rather than per
step, so run `scripts/verify_atom_conservation.py` once on the net
precursor-product pair, not on each intermediate.

--------------------------------------------------------------------------------

## 5. Programmatic API

The same spec can be assembled in Python and compiled without touching disk.
Splitting `rxn_def` (species and conditions) from `spec` (arrows and layout)
lets one set of structures be rendered with several arrow hypotheses.

```python
import sys
sys.path.insert(0, "scripts")  # run from the skill directory
from render_mechanism import compile_clean_reaction

rxn_def = {
    "id": "kinetic_enolization",
    "r1": {"smi": "CC1CCCC(=O)C1", "name": "3-methylcyclohexan-1-one", "w": 4.5},
    "r2": {"smi": "CC(C)[N-]C(C)C", "name": "LDA", "w": 3.8},
    "p1": {"smi": "CC1CCCC([O-])=C1", "name": "Lithium Enolate", "w": 4.5},
    "cond": {
        "top": r"$\text{THF}, -78^\circ\text{C}$",
        "bottom": r"$-\text{HN}(i\text{-Pr})_2$",
    },
}

spec = {
    "add_h": {"r1": 7},
    "arrows": [
        {"type": "attack", "src": ("r2", 3), "tgt": ("r1", 8)},
        {"type": "bond_to_bond", "mol": "r1", "b1": (7, 8), "b2": (5, 7)},
        {"type": "bond_cleave", "mol": "r1", "bond": (5, 6), "tgt": 6},
    ],
}

# Renders TikZ and writes a 300 DPI PNG to out_dir.
# Returns a result tuple on success, or a "[SKIP] ..." string.
result = compile_clean_reaction(rxn_def, out_dir="output", spec=spec)
```

Each call prints `[AUDIT]` lines to stdout: one per arrow that passes closer
than 0.5 mm to a bond or crosses an atom label, one if an arrow produced no
geometry and is therefore missing from the diagram, and otherwise a single line
giving the closest approach in the scene. An arrow whose `type` is not one of
the five above raises `ValueError` rather than being dropped.

This example uses bare indices to show that they work; `add_h` forces the C-H at
index 7 to be drawn so the `bond_to_bond` arrow has a bond to start from.
