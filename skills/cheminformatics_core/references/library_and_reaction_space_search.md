# High-Throughput Virtual Screening, Library Search, & Reaction Spaces

This guide details the execution of **actual high-throughput virtual screens**
across compound libraries (SDF, SDF.GZ, CSV, SMILES) and live chemical databases
(PubChem, ChEMBL), as well as on-demand synthon-space traversal
($10^{10}-10^{36}$ compounds).

--------------------------------------------------------------------------------

## 1. Virtual Screening Invariants & Mandatory User Scope Alignment

When prompted to "run a screen", "find analogs", or "filter a library":

1.  **Mandatory Screening Strategy Alignment**:

    *   You MUST clarify before executing any search:
        -   **Screening Archetype**: Broad commercial building blocks, targeted
            core heterocycles (e.g. pyrazoles, indazoles), polar/solubilizing
            modifications (e.g. morpholines, piperazines), or lead-like property
            space.
        -   **Screening Source**: PubChem commercial vendors (Enamine,
            Chemspace, MolPort), ChEMBL bioactives, local screening decks, or
            on-demand synthon spaces.
        -   **Target Depth**: Exhaustive unconstrained streaming natively vs.
            Top $N$ diverse commercial compounds.
    *   *Never guess the screening strategy or silently insert arbitrary
        category caps without aligning naturally based on intent first.*

2.  **Never Fabricate Candidates**: You must stream or iterate across a real
    compound file (`.sdf`, `.csv`, `.smi`) or live API stream
    (`pubchem-database` / PubChem fastsubstructure). Never hardcode a static
    array of molecules in Python (`candidates = [...]`).

3.  **Multi-Tier Filtering Pipeline**:

    *   **Tier 1 (Substructure Filter)**: Fast bitwise Pattern Fingerprint
        inclusion (`AllProbeBitsMatch`) $\to$ Graph isomorphism.
    *   **Tier 2 (Similarity Filter)**: Vectorized Morgan fingerprint Tanimoto /
        Tversky distance.
    *   **Tier 3 (Physicochemical Filters)**: Molecular weight, cLogP, TPSA,
        Rotatable Bonds, HBD, HBA, and PAINS filters.

4.  **Screening Attrition Reporting**: Always compute and report formal
    screening attrition statistics:

    *   Total molecules evaluated ($N$)
    *   Substructure / Similarity filter passes ($k_1$)
    *   Property filter passes ($k_2$)
    *   Overall Hit Rate (%) and timing.

--------------------------------------------------------------------------------

## 2. High-Throughput Screening CLI Usage

Use `scripts/library_search.py` to execute production screens:

```bash
# 1. Screen an SDF or compressed SDF library with Substructure and Property Filters:
uv run scripts/library_search.py \
    --input screening_deck.sdf.gz \
    --query_smarts "c1ncc(Br)cn1" \
    --max_mw 500.0 \
    --max_logp 4.5 \
    --output_csv screen_hits.csv

# 2. Screen Live PubChem Database for Purchasable Building Blocks:
uv run scripts/library_search.py \
    --source pubchem \
    --query_smarts "c1cnn(C2CCNCC2)c1" \
    --max_hits 50 \
    --output_csv pubchem_hits.csv
```

--------------------------------------------------------------------------------

## 3. High-Throughput Streaming in Python

For screening multi-gigabyte chemical files without memory exhaustion:

```python
import gzip
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski, rdFingerprintGenerator

def stream_and_screen(file_path, query_smarts, max_mw=500.0, max_logp=4.5):
    q_mol = Chem.MolFromSmarts(query_smarts)
    q_pat_fp = Chem.PatternFingerprint(q_mol)

    total_screened = 0
    hits = []

    # Open plain or gzip compressed stream
    supplier = Chem.ForwardSDMolSupplier(gzip.open(file_path, 'rb') if file_path.endswith('.gz') else open(file_path, 'rb'))

    for mol in supplier:
        if mol is None:
            continue
        total_screened += 1

        # Tier 1: C++ Fast Bitwise Prescreen
        target_pat_fp = Chem.PatternFingerprint(mol)
        if not DataStructs.AllProbeBitsMatch(q_pat_fp, target_pat_fp):
            continue

        # Tier 2: Substructure Isomorphism Match
        if not mol.HasSubstructMatch(q_mol):
            continue

        # Tier 3: Physicochemical Property Filters
        mw = Descriptors.ExactMolWt(mol)
        logp = Crippen.MolLogP(mol)
        if mw > max_mw or logp > max_logp:
            continue

        hits.append(mol)

    print(f"Screened {total_screened} molecules -> Found {len(hits)} hits ({len(hits)/total_screened*100:.2f}% hit rate)")
    return hits
```
