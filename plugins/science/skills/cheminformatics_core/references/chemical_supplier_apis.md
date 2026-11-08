# Commercial Chemical Sourcing, CAS Registry, & Building Block APIs

This guide details commercial chemical vendor sourcing (**Chemspace**,
**Enamine**, **MolPort**, **PubChem**), substructure-based building block
discovery, and verified CAS Registry Number identification using the production
tools in `scripts/`.

--------------------------------------------------------------------------------

## 1. Summary of Chemical Sourcing & Building Block Platforms

-   **PubChem PUG-REST** (>110M compounds):
    -   *Access Protocol*: Open REST (Zero Auth).
    -   *Key Capabilities*: Fast substructure SMARTS search, vendor registries,
        verified CAS numbers.
    -   *Best Used For*: Zero-cost vendor availability checking & unconstrained
        catalog search.
-   **Chemspace** (>100B molecules):
    -   *Access Protocol*: REST API v2 (`api.chem-space.com`).
    -   *Key Capabilities*: In-stock Building Blocks (`CSBB`), screening
        compounds, Freedom spaces.
    -   *Best Used For*: Multi-vendor aggregated purchasing & building block
        sourcing.
-   **EnamineStore** (~3.8M stock + >38B REAL):
    -   *Access Protocol*: REST API (`enaminestore.com/api/v1/`).
    -   *Key Capabilities*: In-stock Building Blocks (BB), screening decks,
        custom synthesis.
    -   *Best Used For*: Primary building block sourcing & REAL library orders.
-   **MolPort** (~7M stock + >100M make):
    -   *Access Protocol*: REST API (`api.molport.com/v1/`).
    -   *Key Capabilities*: Multi-supplier stock levels (>70 vendors), tier
        pricing, delivery times.
    -   *Best Used For*: Comprehensive vendor aggregation & tier pricing.

--------------------------------------------------------------------------------

## 2. Robust CAS Registry Number Identification Protocol

A CAS Registry Number has the format `[2-7 digits]-[2 digits]-[1 digit]`.

### The Mod-10 Checksum Algorithm

For a CAS number $d_k d_{k-1} \dots d_2 d_1 - R$: $$\text{Checksum} = \left(
\sum_{i=1}^k i \cdot d_i \right) \bmod 10 == R$$ *Example*: For 50-00-0
(Formaldehyde), $5 \times 1 = 5 \equiv 0 \pmod{10}$ (Valid). Never invent CAS
numbers or accept unverified strings.

### Production CLI (`chemical_vendor_resolver.py`):

This tool is **single-query** and exposes two subcommands. To enrich a whole
table, loop over its rows in a shell or Python driver and collect the per-query
JSON.

Verify a CAS Registry Number (optionally round-tripping it against an expected
structure):

```bash
uv run scripts/chemical_vendor_resolver.py resolve-cas \
    --cas 50-00-0 \
    --smiles "C=O" \
    --output scratch/formaldehyde_cas.json
```

Look up commercial vendors for a SMILES, chemical name, or PubChem CID:

```bash
uv run scripts/chemical_vendor_resolver.py lookup-suppliers \
    --query "Cc1ccc(Nc2nccc(-c3cccnc3)n2)cc1" \
    --output scratch/vendors.json
```

*   Queries PubChem PUG-REST by InChIKey for official synonyms and verified CAS
    Registry Numbers.
*   Verifies Mod-10 checksum on every extracted CAS string.
*   Generates live vendor section linkouts
    (`https://pubchem.ncbi.nlm.nih.gov/compound/<CID>#section=Chemical-Vendors`).

--------------------------------------------------------------------------------

## 3. High-Throughput Commercial Reactant Discovery (`find_coupling_partners.py`)

When searching for real commercial building blocks to couple with a core
starting material:

```bash
uv run scripts/find_coupling_partners.py \
    --core "B(O)Oc1ccccc1" \
    --reaction_smarts "[c:1][B]([OH])[OH].[c:2][Cl,Br,I]>>[c:1][c:2]" \
    --all \
    --output_csv scratch/commercial_aryl_halides.csv \
    --output_html scratch/commercial_aryl_halides.html
```

> [!IMPORTANT] **NO SILENT THROTTLING**: Never insert artificial
> `max_records=50` limits into building block search loops. Use `--all` for
> unconstrained streaming, or align with the user via `ask_question` on the
> target hit depth.
