# First-Principles Reaction SMARTS Construction & Building Block Sourcing

This guide provides a first-principles manual for **constructing Reaction SMARTS
for any arbitrary chemical reaction**, validating them programmatically in
RDKit, extracting unmapped reactant search templates, and querying commercial
suppliers for real building blocks.

## Standard CLI Executables (`scripts/`)

-   `chemical_reaction.py`: Starting core reactive handle analysis, reaction
    SMARTS execution, and regioselective product enumeration.
-   `substructure_search.py`: SMARTS substructure filtering, core scaffold
    matching, and regioisomer disambiguation.
-   `enumerate_matrix_library.py`: Combinatorial matrix library enumeration (`R1
    x R2`) with Lipinski/Ro5 property gating.
-   `find_coupling_partners.py`: High-throughput commercial coupling partner
    discovery across PubChem and vendor catalogs.

--------------------------------------------------------------------------------

## 1. The First-Principles Algorithm for Constructing Arbitrary Reaction SMARTS

A Reaction SMARTS represents a chemical transformation as:
$$\text{Reactant}_1\text{ . }\text{Reactant}_2\text{ }>>\text{
}\text{Product}_1\text{ . }\text{Product}_2$$

To construct a valid Reaction SMARTS for **any arbitrary reaction**:

### Step 1: Identify Reacting Centers (Bonds Broken & Bonds Formed)

1.  Determine which specific atoms participate in the bond breakage or bond
    formation.
2.  Only include the reacting functional group atoms and their immediate
    connectivity in the reactant templates. Do **not** hardcode spectator alkyl
    chains or full scaffold structures.

### Step 2: Assign Atom Mapping Numbers (`:n`)

*   Every atom that persists from reactants into the product **must** have an
    identical integer map number (`:1`, `:2`, `:3`, ...) inside square brackets.
*   Example (Amide bond formation):
    *   Carboxylic carbon: `[C:1]`
    *   Carbonyl oxygen: `[O:2]`
    *   Amine nitrogen: `[N:3]`
    *   Product amide: `[C:1](=[O:2])[N:3]`

### Step 3: Handle Leaving Groups & Unmapped Atoms

*   **Leaving Groups**: Atoms in the reactant template that are discarded during
    the reaction (e.g. `[OH]` on a carboxylic acid, `[Cl,Br,I]` on an aryl
    halide, `[B](OH)2` on a boronic acid, or `[H]` on a nucleophile) must
    **NOT** have atom map numbers. RDKit automatically eliminates unmapped
    reactant atoms.
*   **New Atoms in Products**: Any unmapped atom in the product template is
    synthesized with default implicit valences.

### Step 4: Specify Atom Environment & Specificity

Use standard SMARTS atom primitives to avoid over-matching or under-matching:

*   **Aromatic vs Aliphatic**: `[c:1]` (aromatic carbon) vs `[C:1]` (aliphatic
    carbon).
*   **Hydrogen Counts**: `[N;H2,H1:3]` (matches only primary and secondary
    amines, excluding tertiary).
*   **Heteroatom Disjunctions**: `[O,S:2]` (matches ether or thioether),
    `[Cl,Br,I:3]` (matches reactive halogens).
*   **Ring / Charge Constraints**: `[r5,r6:1]` (ring size 5 or 6), `[N+:2]`
    (quaternary/charged nitrogen).

--------------------------------------------------------------------------------

## 2. Programmatic Construction & Validation in Python

Use `AllChem.ReactionFromSmarts` to construct and validate arbitrary reactions:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Programmatic Assembly:
reactant1 = "[C:1](=[O:2])[OH]"     # Carboxylic Acid (OH is unmapped leaving group)
reactant2 = "[N;H2,H1:3]"           # 1° or 2° Amine
product   = "[C:1](=[O:2])[N:3]"     # Amide

# Construct full SMARTS string:
rxn_smarts = f"{reactant1}.{reactant2}>>{product}"

# Parse and Validate:
rxn = AllChem.ReactionFromSmarts(rxn_smarts)
assert rxn is not None, "Syntax error in constructed SMARTS"

n_warn, n_err = rxn.Validate()
if n_err > 0:
    raise ValueError(f"Reaction SMARTS validation failed with {n_err} error(s)")

print(f"Validated Reaction SMARTS: {rxn_smarts}")
print(f"  Reactant Templates: {rxn.GetNumReactantTemplates()}")
print(f"  Product Templates:  {rxn.GetNumProductTemplates()}")
```

--------------------------------------------------------------------------------

## 3. Extracting Search SMARTS & Sourcing Commercial Building Blocks

Once a reaction SMARTS is constructed, extract the **unmapped search query** for
each reactant template to search commercial vendor catalogs (PubChem, Enamine,
Chemspace, MolPort) and execute transformations at high throughput:

### Production CLI Recipe: High-Throughput Bulk Reaction Discovery

Use `find_coupling_partners.py` to query PubChem's fastsubstructure database,
download properties via vectorized POST, and stream RDKit in-silico reactions
directly to disk:

```bash
uv run scripts/find_coupling_partners.py \
    --core "<CORE_SMILES>" \
    --reaction_smarts "<REACTION_SMARTS>" \
    --partner_smarts "<OPTIONAL_PARTNER_SMARTS>" \
    --output_csv scratch/commercial_derivatives.csv \
    --output_html scratch/commercial_derivatives.html
```

-   Automatically retrieves all matching commercial coupling partners from
    PubChem without artificial caps.
-   Executes in-silico reactions and profiles MW, cLogP, TPSA, HBD, HBA, RotB,
    QED, and SA score.
-   Provides direct commercial vendor purchase linkouts (PubChem chemical vendor
    directory) for each reactant.
-   Generates interactive 2D structure grid (`.html`) with SAR scatter plots.

--------------------------------------------------------------------------------

## 4. R-Group Decomposition (`rdRGroupDecomposition`)

Decomposes a series of related analogs into a common core scaffold and
substituent R-groups:

```python
from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition

core = Chem.MolFromSmarts("c1cc([*:1])cc([*:2])c1")
molecules = [
    Chem.MolFromSmiles("c1cc(Cl)cc(O)c1"),
    Chem.MolFromSmiles("c1cc(F)cc(N)c1"),
]

rgd = rdRGroupDecomposition.RGroupDecomposition(core)
for mol in molecules:
    rgd.Add(mol)
rgd.Process()

rgroups = rgd.GetRGroupsAsColumns()
```

--------------------------------------------------------------------------------

## 5. Fragment-Based Decomposition (Murcko & BRICS)

### Murcko Scaffold Extraction

```python
from rdkit.Chem.Scaffolds import MurckoScaffold

mol = Chem.MolFromSmiles("Cc1ccc(Nc2nccc(-c3cccnc3)n2)cc1")
scaffold = MurckoScaffold.GetScaffoldForMol(mol)
generic_scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)
```

### BRICS Synthon Decomposition

```python
from rdkit.Chem import BRICS

mol = Chem.MolFromSmiles("CC(=O)Nc1ccc(OC(=O)c2ccccc2)cc1")
fragments = BRICS.BRICSDecompose(mol)
```

--------------------------------------------------------------------------------

## 6. Strategic Retrosynthesis & Building Block Sourcing

Use `scripts/retrosynthesis_decomposer.py` to perform strategic 1-step
disconnections, RECAP tree decompositions, and commercial precursor availability
checks:

```bash
# 1. Full Retrosynthetic Analysis & Commercial Precursor Sourcing:
uv run scripts/retrosynthesis_decomposer.py \
    --smiles "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C" \
    --recap \
    --brics \
    --output_csv retrosynthesis_plan.csv
```

### In-Memory Retrosynthetic Disconnections in Python

```python
from rdkit import Chem
from rdkit.Chem import AllChem, Recap, BRICS
import pandas as pd

# Strategic 1-Step Disconnection Rules (with exocyclic !@ constraints)
DISCONNECTIONS = [
    ("Amide Coupling", "[C:1](=[O:2])!@[N;!$(N=O):3]>>[C:1](=[O:2])O.[N:3]"),
    ("Suzuki-Miyaura Biaryl", "[c:1]!@[c:2]>>[c:1]B(O)O.[c:2]Br"),
    ("Buchwald-Hartwig C-N", "[c:1]!@[N;!$(NC=O);!$(N=O):2]>>[c:1]Br.[N:2]"),
    ("Sulfonamide Formation", "[S:1](=[O:2])(=[O:3])!@[N:4]>>[S:1](=[O:2])(=[O:3])Cl.[N:4]"),
    ("Esterification", "[C:1](=[O:2])!@[O:3][C:4]>>[C:1](=[O:2])O.[O:3][C:4]"),
    ("Ether / SNAr Alkylation", "[c:1]!@[O:2][C:3]>>[c:1][O:2].[C:3]Br"),
    ("Reductive Amination", "[c:1][C;H2:2]!@[N;!$(NC=O);!$(N=O):3]>>[c:1][C:2]=O.[N:3]"),
]

def disconnect_target(target_smiles: str):
    mol = Chem.MolFromSmiles(target_smiles)
    pathways = []
    for rxn_name, smarts in DISCONNECTIONS:
        rxn = AllChem.ReactionFromSmarts(smarts)
        for outcome in rxn.RunReactants((mol,)):
            p1, p2 = outcome[0], outcome[1]
            Chem.SanitizeMol(p1)
            Chem.SanitizeMol(p2)
            pathways.append({
                "Reaction": rxn_name,
                "Precursor_A": Chem.MolToSmiles(p1),
                "Precursor_B": Chem.MolToSmiles(p2),
            })
    return pd.DataFrame(pathways)
```
