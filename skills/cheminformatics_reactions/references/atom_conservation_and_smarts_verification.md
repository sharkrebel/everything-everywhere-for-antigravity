# Atom Conservation, Intramolecular SMARTS, & Formal Mechanism Verification

This reference details reaction stoichiometry validation
(`scripts/verify_atom_conservation.py`), intramolecular Reaction SMARTS syntax
and 6-tier validation (`scripts/reaction_smarts.py`), and formal graph-theoretic
electron-pushing verification (`scripts/electron_graph.py`).

--------------------------------------------------------------------------------

## 1. Strict Atom Conservation (`scripts/verify_atom_conservation.py`)

In organic transformations, cyclizations, and macrocyclizations, the product
cannot gain or lose heavy atoms ($\Delta C = 0$ and heavy-atom balance after
accounting for explicit byproducts). Always validate every precursor-product
pair via CLI:

```bash
uv run scripts/verify_atom_conservation.py \
    --precursor "<SMILES>" \
    --product "<MACROCYCLE_SMILES>" \
    --reaction-type amide \
    --macrocycle-size 12
```

### Multi-Component & Byproduct Rules:

-   **Multi-component reactions**: Pass each reactant separately, e.g.
    `--precursor "CC(=O)O" "CCO" --product "CC(=O)OCC" --byproducts H2O`.
-   **Byproduct stoichiometry**: Prefix with `<N>*` when more than one
    equivalent is released, e.g. `--byproducts "2*H2O"` for a bis-amidation.
-   **Formula vs. SMILES disambiguation**: Byproducts may be molecular formulas
    (`B(OH)3`, `Al(OH)3`) or SMILES. Tokens valid as both (e.g. `CCO` =
    $\text{C}_2\text{O}$ vs. ethanol) are flagged as ambiguous; prefix with
    `smiles:CCO` or `formula:CCO` to disambiguate.
-   **Hydrogen balance (`--strict-hydrogen`)**: By default, $\Delta C = 0$ and
    heavy-atom balance are strictly enforced while non-zero $\Delta H$ logs a
    `[NOTE]` to accommodate ionic states. Pass `--strict-hydrogen` to require
    exact $\Delta H = 0$.

--------------------------------------------------------------------------------

## 2. Intramolecular Reaction SMARTS Parentheses Rule (`scripts/reaction_smarts.py`)

In RDKit Reaction SMARTS, top-level dots (`.`) define separate reactant
molecules (`GetNumReactantTemplates() > 1`). For an **intramolecular ring
closure on a single bifunctional precursor**, disconnected reacting functional
groups **must be enclosed in outer parentheses**:

-   ❌ **Wrong (Intermolecular Dimerization Trap)**:
    `[C:1](=[O:2])[OH].[N;H1:3]>>[C:1](=[O:2])[N:3]`
-   ✅ **Correct (Intramolecular Ring Closure)**:
    `([C:1](=[O:2])[OH].[N;H1:3])>>[C:1](=[O:2])[N:3]`

```bash
# Validate an intramolecular Reaction SMARTS
uv run scripts/reaction_smarts.py --smarts "([C:1](=[O:2])[OH].[N;H1:3])>>[C:1](=[O:2])[N:3]" --intramolecular --output /tmp/smarts_report.json

# Export ReactionSMARTSConverter electron-pushing graph metadata (broken/formed bonds & charge deltas)
uv run scripts/reaction_smarts.py --smarts "([C:1](=[O:2])[OH].[N;H1:3])>>[C:1](=[O:2])[N:3]" --to-graph-json /tmp/electron_graph.json

# Automatically wrap an intermolecular SMARTS into an intramolecular ring-closure SMARTS
uv run scripts/reaction_smarts.py --make-intramolecular "[C:1](=[O:2])[OH].[N;H1:3]>>[C:1](=[O:2])[N:3]"

# List built-in reaction templates
uv run scripts/reaction_smarts.py --list-templates
```

--------------------------------------------------------------------------------

## 3. Formal Mechanism Graph Verification (`ElectronGraph` Python API)

`scripts/electron_graph.py` enforces valence electron conservation, formal
charge balance, octet/non-hypervalency capacity, spectator locality, Tier 1–4
stereochemical and ring-strain rules, and product matching:

```python
from cheminformatics_reactions.scripts.electron_graph import (
    CurvedArrow,
    MolecularGraph,
    ReactionVerifierEngine,
    normalize_dipolar_rdkit_mol,
)

# 1. Build molecular state (or ingest via MolecularGraph.from_smiles("CCO", strict_lewis=True))
g = MolecularGraph("SN2 Reactants")
g.add_atom("O_nuc", "O", lone_pairs=3, formal_charge=-1)
g.add_atom("H_nuc", "H")
g.add_bond("O_nuc", "H_nuc", order=1)

g.add_atom("C", "C")
g.add_atom("Br", "Br", lone_pairs=3)
g.add_bond("C", "Br", order=1)
for i in range(3):
    g.add_atom(f"H{i}", "H")
    g.add_bond("C", f"H{i}", order=1)

# 2. Define electron-pushing curved arrows (2e- pairs or 1e- fishhooks)
arrows = [
    CurvedArrow(
        source_type="lone_pair",
        source_id="O_nuc",
        sink_type="bond",
        sink_id=("O_nuc", "C"),
        label="Nucleophilic attack",
    ),
    CurvedArrow(
        source_type="bond",
        source_id=("C", "Br"),
        sink_type="atom",
        sink_id="Br",
        label="Leaving group expulsion",
    ),
]

# 3. Apply graph rewriting and verify physical & stereochemical invariants (Tiers 1-4)
valid, prod, audit = ReactionVerifierEngine.apply_and_verify_step(
    g, arrows, step_name="SN2 Displacement"
)
print(f"Verification Status: {audit['status']}")

# 4. Validate against expected product SMILES
is_match, note = ReactionVerifierEngine.verify_product_match(prod, "CO.[Br-]")
```

### A. Automated Curved-Arrow Synthesis from Atom-Mapped SMILES (`solve_arrows_from_mapped_smiles`)

When given atom-mapped reactant and product SMILES (or a single
`"reactants>>products"` reaction SMILES, or `MolecularGraph` instances via
`solve_arrows_from_mapped_graphs`),
`ReactionVerifierEngine.solve_arrows_from_mapped_smiles(reactant_smiles,
product_smiles=None, allow_multiproton=False)` solves the Dugundji–Ugi reaction
difference matrix $\Delta R = R_{\text{prod}} - R_{\text{react}}$ via directed
vertex-flow conservation and returns `(reactant_graph, product_graph, arrows)`
(automatically detecting $2e^-$ polar / pericyclic / vinylogous steps vs. $1e^-$
radical fishhook steps from the reactant/product radical states):

```python
# Polar / concerted 2e- arrow synthesis (SN2, E2, Diels-Alder, 1,5-H shifts, etc.)
rg, pg, arrows = ReactionVerifierEngine.solve_arrows_from_mapped_smiles(
    "[OH-:1].[CH3:2][Br:3]",
    "[OH:1][CH3:2].[Br-:3]",
)

# Radical 1e- fishhook arrow synthesis (homolysis, H-abstraction, radical addition)
rg_rad, pg_rad, rad_arrows = ReactionVerifierEngine.solve_arrows_from_mapped_smiles(
    "[Cl:1][Cl:2]>>[Cl:1].[Cl:2]"
)
```

### B. Strict Lewis Electron-Pair Parity & Dipolar Normalization (`strict_lewis=True`)

Pass `strict_lewis=True` to `MolecularGraph.from_smiles(smi, strict_lewis=True)`
(or `MolecularGraph.from_rdkit(mol, strict_lewis=True)`) to:

1.  Normalize pentavalent 2nd-row nitrogen representations in place via
    `normalize_dipolar_rdkit_mol(mol)` (converting `-N(=O)=O`, `-N=[N+]=N`,
    `C=[N+]=N`, $N$-oxides, and nitrile oxides into octet-compliant
    charge-separated Lewis forms without altering atom indices or atom-map
    numbers).
2.  Enforce strict Lewis electron-pair parity (`rem_e >= 0` and `rem_e % 2 ==
    0`), raising a `ValueError` on inconsistent formal charge / radical states
    rather than silently truncating via integer division.

### C. Tier 4 Ring Strain & Bredt's Rule Bridgehead Verification

During `apply_and_verify_step`, Tier 4 (`audit["errors"]` and
`audit["invariants"]["stereochemical_integrity"]`) automatically checks:

1.  **Small-ring sp strain**: Rejects triple bonds or cumulated allenes inside
    3- to 7-membered rings.
2.  **Trans-cycloalkene strain**: Rejects `trans` (`/` / `\`) double bonds
    inside 3- to 7-membered rings.
3.  **Bredt's Rule Bridgehead Double-Bond Strain**: Computes the bridged
    bicyclic path lengths $(x, y, z)$ ($x \ge y \ge z \ge 1$) via
    `ReactionVerifierEngine._check_bredt_bridgehead_strain` for every newly
    formed ring double bond at a bridgehead atom and flags anti-Bredt bridgehead
    alkenes with Fawcett $S = x + y + z < 7$ (e.g., $\Delta^{1,2}$-norbornene,
    bicyclo[2.2.1]hept-1-ene, $S = 4$), while permitting isolable $S \ge 7$
    systems (e.g., bicyclo[3.3.1]non-1-ene, $S = 7$) and non-bridged fused rings
    ($z = 0$, such as decalin).

--------------------------------------------------------------------------------

## 4. Automated Reaction SMARTS Generation & 6-Tier Validation (`Python API`)

```python
from cheminformatics_reactions.scripts.reaction_smarts import (
    ReactionSMARTSGenerator,
    ReactionSMARTSValidator,
)

# 1. Derive Reaction SMARTS directly from an ElectronGraph step
res = ReactionSMARTSGenerator.from_electron_graph(
    g, arrows, radius=1, collapse_spectator_hydrogens=True
)
print("Generated Reaction SMARTS:", res["reaction_smarts"])

# 2. Perform 6-tier deep validation with live dry-run on sample reactants
report = ReactionSMARTSValidator.validate(
    res["reaction_smarts"], sample_reactants=["[OH-]", "CCBr"]
)
print("Valid:", report["valid"], "Products:", report["tier6_live_execution"]["products"])
```

--------------------------------------------------------------------------------

## 5. Prochirality, Topicity, & Tier 3D Stereogenesis (`scripts/prochirality.py`)

RDKit's default `Chem.CanonicalRankAtoms` only ranks 2D topology and existing
chiral centers, leaving prochiral (`pro-R`/`pro-S`, `Re`/`Si`) and geminal
diastereotopic (`-CH2-`) groups unranked. `scripts/prochirality.py` and
`ElectronGraph` Tier 3D enforce:

1.  **Tetrahedral ($sp^3$) Topicity (`pro-R` / `pro-S`)**: Applies Hanson's Rule
    (IUPAC P-92.5.1) via valid heavier NIST isotopes (`rdCIPLabeler`) to
    classify paired ligands as `homotopic`, `enantiotopic`, or `diastereotopic`,
    and computes `achiral_classes_by_element` vs. `chiral_classes_by_element`.
2.  **Trigonal ($sp^2$) Facial Topicity (`Re` / `Si`)**: Ranks $sp^2$
    substituents ($1 > 2 > 3$) with a `Z=0` phantom probe (IUPAC P-92.5.2) and
    resolves `CurvedArrow(..., stereo_policy="Re" | "Si" | "R" | "S")` into the
    exact product tetrahedral configuration, while rejecting stereoselective
    policies on homotopic centers.
3.  **Reaction Stereogenesis & False Chirality Audit**: Rejects `@`/`@@` tags
    placed on homotopic atoms and flags achiral-to-single-enantiomer reactions
    that omit a chiral catalyst/auxiliary.
