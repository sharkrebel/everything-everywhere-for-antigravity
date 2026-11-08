# RDKit 2D Cheminformatics API Reference

Verified module paths, function signatures, argument types, and return objects
for modern RDKit (2022+ / 2024+ / 2026+).

## Standard CLI Executables (`scripts/`)

-   `descriptors.py`: Canonical 13-property descriptor engine, column name
    normalization (DescriptorCanonicalizer), and Pfizer CNS MPO scorer.
-   `filter_catalogs.py`: PAINS (A/B/C), Brenk, NIH, and ZINC structural alert
    screening, Lipinski/Veber rules, and CNS MPO reporting.
-   `standardize_molecule.py`: Salt stripping, charge neutralization, and
    tautomer canonicalization.
-   `compute_descriptors.py`: Lipinski Ro5, Veber, QED, TPSA, and
    physicochemical descriptors.
-   `bioassay_data_curator.py`: Automated ChEMBL/PubChem target bioactivity
    fetching and pChEMBL normalization.
-   `nist_webbook_fetcher.py`: Official NIST Chemistry WebBook experimental
    thermochemical & physical property client ($\Delta_f H^\circ$,
    $T_{\text{boil}}$, $T_{\text{fus}}$, Antoine vapor pressure constants,
    ionization energy, and official 2D `.mol` structures).
-   `pareto_frontier_optimizer.py`: Multi-objective Pareto frontier ranking
    across potency, LipE, and ADMET properties.
-   `../rdkit_general/scripts/draw_molecules.py`: Publication-quality 2D
    structure and reaction scheme rendering.

--------------------------------------------------------------------------------

## 1. Core Molecular I/O & Manipulation (`rdkit.Chem`)

```python
from rdkit import Chem

# Creation from SMILES / SMARTS / InChI
mol = Chem.MolFromSmiles(smiles_str, sanitize=True)  # Returns Mol or None
query = Chem.MolFromSmarts(smarts_str)               # Query Mol for substructure
mol = Chem.MolFromInchi(inchi_str, sanitize=True)
inchi_str = Chem.MolToInchi(mol)
inchikey = Chem.MolToInchiKey(mol)

# Bulk Suppliers (SDF / SMILES)
suppl = Chem.SDMolSupplier(sdf_path, sanitize=True, removeHs=True)
# Streaming for large files:
f_suppl = Chem.ForwardSDMolSupplier(file_obj, sanitize=True)

# Export
smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
Chem.MolToMolFile(mol, "output.mol")

# Sanitization & Stereochemistry Perception
Chem.SanitizeMol(mol)
Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)  # [(idx, 'R'|'S')]
```

--------------------------------------------------------------------------------

## 2. Molecular Standardization (`rdkit.Chem.MolStandardize`)

```python
from rdkit.Chem.MolStandardize import rdMolStandardize

# All-in-One Parent Standardizer (Removes salts, neutralizes, tautomer canonicalization)
standard_parent = rdMolStandardize.SuperParent(mol)

# Stepwise Components
clean_mol = rdMolStandardize.Cleanup(mol)
parent_mol = rdMolStandardize.FragmentParent(clean_mol)  # Largest organic fragment
uncharged_mol = rdMolStandardize.Uncharger().uncharge(parent_mol)
tautomer = rdMolStandardize.TautomerEnumerator().Canonicalize(uncharged_mol)
```

--------------------------------------------------------------------------------

## 3. Descriptors & Physicochemical Properties

```python
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski, Crippen, QED

# Mass & Lipophilicity
mw = Descriptors.MolWt(mol)
logp = Crippen.MolLogP(mol)               # or Descriptors.MolLogP(mol)
mr = Descriptors.MolMR(mol)
tpsa = rdMolDescriptors.CalcTPSA(mol)

# Hydrogen Bonding & Flexibility
hbd = Lipinski.NumHDonors(mol)            # or rdMolDescriptors.CalcNumHBD(mol)
hba = Lipinski.NumHAcceptors(mol)         # or rdMolDescriptors.CalcNumHBA(mol)
rot_bonds = Lipinski.NumRotatableBonds(mol)
fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)

# Drug-Likeness & Structural Alerts (PAINS)
qed_score = QED.qed(mol)
from rdkit.Chem import FilterCatalog
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
catalog = FilterCatalog.FilterCatalog(params)
has_pains = catalog.HasMatch(mol)
# Note: Structural alert catalogs (PAINS, Brenk) assess drug-likeness of final products or lead candidates.
# NEVER apply them to building blocks/reactants, whose reactive synthetic handles will be falsely flagged.
```

### Canonical Descriptors & Pfizer CNS MPO Engine (`descriptors.py`)

`scripts/descriptors.py` is the single authoritative source of truth for Pfizer
CNS MPO scoring, canonical 13-property featurization, and descriptor column
normalization. **Do NOT import `calculate_cns_mpo` from `compat.py`**
(`compat.py` is strictly an RDKit runtime C++ version compatibility shim).

```python
from descriptors import (
    calculate_cns_mpo,
    compute_standard_descriptors,
    normalize_descriptor_columns,
    DescriptorCanonicalizer,
)
from rdkit import Chem

mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")

# Pfizer CNS MPO Score [0.0 - 6.0] (Wager et al., ACS Chem. Neurosci. 2010/2016)
# Evaluates 6 desirability components: cLogP, cLogD, MW, TPSA, HBD, pKa.
# Scores >= 4.0 indicate high probability of blood-brain barrier penetrance.
cns_mpo = calculate_cns_mpo(mol)  # returns 5.83 for aspirin

# Canonical 13-property featurization
props = compute_standard_descriptors(mol, style="short")  # or style="snake"

# Column name normalization for DataFrames
df_norm = normalize_descriptor_columns(df, target_style="short")
```

--------------------------------------------------------------------------------

## 4. Modern Fingerprints & Similarity Metrics

```python
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator, MACCSkeys

# Explicit Generator Pattern (RDKit 2022+)
morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fp1 = morgan_gen.GetFingerprint(mol1)                    # ExplicitBitVect
fp_numpy = morgan_gen.GetFingerprintAsNumPy(mol1)        # numpy uint8 array

# Similarity
tanimoto = DataStructs.TanimotoSimilarity(fp1, fp2)
bulk_sims = DataStructs.BulkTanimotoSimilarity(query_fp, fp_list)  # Fast C++ vector
```

--------------------------------------------------------------------------------

## 5. Substructure Searching & Chemical Reactions

```python
from rdkit.Chem import AllChem, rdFMCS

# Substructure Match
query = Chem.MolFromSmarts("c1ccccc1")
has_match = mol.HasSubstructMatch(query)
matches = mol.GetSubstructMatches(query, uniquify=True)

# Chemical Reactions from SMARTS
rxn_smarts = "[C:1](=[O:2])[OH].[N;H2,H1:3]>>[C:1](=[O:2])[N:3]"
rxn = AllChem.ReactionFromSmarts(rxn_smarts)
products = rxn.RunReactants((acid_mol, amine_mol))
# Always sanitize reaction products:
for outcome in products:
    for prod in outcome:
        Chem.SanitizeMol(prod)
```

> **ChEMBL Substructure API Gotcha**: ChEMBL's REST `/substructure` endpoint
> runs `Chem.SanitizeMol()` server-side and crashes with an HTTP 500 error
> (`Boost.Python.ArgumentError: SanitizeMol(NoneType)`) if given a SMARTS query
> containing atom/bond query primitives (`[#6]`, `~`, `:`). Always pass a valid
> kekulized SMILES string to ChEMBL `/substructure`, or filter downloaded ChEMBL
> target bioactivities locally using RDKit
> `mol.HasSubstructMatch(Chem.MolFromSmarts(...))`.

--------------------------------------------------------------------------------

## 6. 2D Structure Rendering (`rdkit.Chem.Draw`)

```python
from rdkit.Chem.Draw import rdMolDraw2D, rdDepictor

rdDepictor.Compute2DCoords(mol)
prep_mol = rdMolDraw2D.PrepareMolForDrawing(mol)

# High-Resolution Headless PNG
drawer_png = rdMolDraw2D.MolDraw2DCairo(400, 400)
drawer_png.DrawMolecule(prep_mol)
drawer_png.FinishDrawing()
png_bytes = drawer_png.GetDrawingText()
```

--------------------------------------------------------------------------------

## 7. NIST Chemistry WebBook Thermochemical & Physical Property Client (`nist_webbook_fetcher.py`)

Queries the official NIST Chemistry WebBook
(`https://webbook.nist.gov/chemistry/`) by CAS Registry Number, InChIKey,
SMILES, or name. Extracts evaluated heats of formation ($\Delta_f
H^\circ_{\text{gas/liquid/solid}}$), boiling/melting points ($T_{\text{boil}},
T_{\text{fus}}$), Antoine vapor pressure constants ($A, B, C$), ionization
energy ($\text{IE}$), and official NIST 2D `.mol` structures.

```bash
# Query by CAS Registry Number
uv run scripts/nist_webbook_fetcher.py \
  --cas "64-17-5" \
  --output-json ethanol_nist.json \
  --output-mol ethanol_nist.mol

# Query by SMILES or Name
uv run scripts/nist_webbook_fetcher.py \
  --smiles "c1ccccc1" \
  --output-json benzene_nist.json
```

--------------------------------------------------------------------------------

## 8. 3D Conformer & Force Field Cross-Reference

For 3D conformer embedding (ETKDGv3), force field minimization (MMFF94s/UFF),
and protein docking, refer to the **`cheminformatics_physics`** skill
(`references/`).
