# RDKit and Pandas Integration Guide (`PandasTools`)

This guide covers modern, error-free integration of RDKit with Pandas DataFrames
(`rdkit.Chem.PandasTools`), covering DataFrame manipulation, high-performance
SDF I/O, substructure filtering, Excel export with embedded structures, and
serialization best practices.

## Standard CLI Executables (`scripts/`)

-   `validate_chemical_table.py`: Schema validation, primary-key deduplication,
    and SMILES/InChIKey sanity checks on chemical CSV/TSV tables.
-   `../rdkit_general/scripts/interactive_grid_viewer.py`: Interactive 2D
    chemical structure grid and SAR scatter explorer HTML generator. This lives
    in `rdkit_general`; there is no copy in this skill.

--------------------------------------------------------------------------------

## 1. Core `PandasTools` Methods & Signatures

```python
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools

# 1. Add Molecule Column from SMILES
# Modifies df in-place, creating 'ROMol' objects
PandasTools.AddMoleculeColumnToFrame(df, smilesCol='SMILES', molCol='ROMol', includeFingerprints=False)

# 2. Add Bemis-Murcko Scaffolds to Frame
PandasTools.AddMurckoToFrame(df, molCol='ROMol', MurckoCol='Murcko_SMILES', Generic=False)

# 3. Strip Salts / Solvents In-Place
PandasTools.RemoveSaltsFromFrame(df, molCol='ROMol')

# 4. Align 2D Coordinates to Core Scaffold
PandasTools.AlignToScaffold(df, molCol='ROMol', scaffoldCol='Murcko_SMILES')
```

--------------------------------------------------------------------------------

## 2. High-Performance SDF I/O

`PandasTools.LoadSDF` reads multi-property SDF files directly into Pandas
DataFrames, populating all SD tags as typed DataFrame columns:

```python
# Load SDF directly to DataFrame
df = PandasTools.LoadSDF(
    "dataset.sdf",
    idName="Compound_ID",     # Tag to use as compound ID
    molColName="ROMol",        # Column to store Mol objects
    includeFingerprints=False,
    isomericSmiles=True
)

# Write DataFrame to SDF preserving specified properties
PandasTools.WriteSDF(
    df,
    "output.sdf",
    molColName="ROMol",
    idName="Compound_ID",
    properties=["IC50_uM", "cLogP", "TPSA", "Murcko_SMILES"],
    allNumeric=False
)
```

--------------------------------------------------------------------------------

## 3. Substructure Filtering on DataFrames

`PandasTools` overloads the `>=` operator on DataFrame columns for concise
substructure search:

```python
# 1. Overloaded Operator Syntax (Fast)
query = Chem.MolFromSmarts("c1ccncc1")
pyridine_hits_df = df[df["ROMol"] >= query]

# 2. Safe Lambda Filtering (Guards against None values)
def has_substruct_safe(mol, query_mol):
    return mol.HasSubstructMatch(query_mol) if mol is not None else False

filtered_df = df[df["ROMol"].apply(lambda m: has_substruct_safe(m, query))]
```

--------------------------------------------------------------------------------

## 4. Exporting with Embedded Molecular Structures

### Excel (`.xlsx`) with Embedded Image Thumbnails

`PandasTools.SaveXlsxFromFrame` embeds actual high-resolution 2D chemical
structure images directly into Excel spreadsheet cells:

```python
# Exports Excel sheet with 150x150 px embedded PNGs in the 'ROMol' column
PandasTools.SaveXlsxFromFrame(df, "compounds.xlsx", molCol="ROMol", size=(150, 150))
```

### HTML Reports with Crisp Vector SVGs

```python
# Set global representation to vector SVG
PandasTools.molRepresentation = "svg"

# CRITICAL: escape=False is mandatory to prevent escaping raw <svg> tags to &lt;svg
html_content = df[["Compound_ID", "SMILES", "MW", "cLogP", "ROMol"]].to_html(escape=False)

with open("report.html", "w") as f:
    f.write(html_content)
```

### Grid Images from DataFrame

```python
from rdkit.Chem import Draw

mols = [m for m in df["ROMol"] if m is not None]
legends = list(df["Compound_ID"])

# Render SVG grid directly
svg_grid = Draw.MolsToGridImage(mols, legends=legends, molsPerRow=4, useSVG=True)
```

--------------------------------------------------------------------------------

## 5. Serialization & Storage Best Practices

`rdkit.Chem.rdchem.Mol` objects cannot be stored natively in columnar formats
like Parquet or Feather.

### Recommended Storage Protocols:

1.  **Canonical SMILES Storage (Default)**:

    ```python
    df_to_save = df.drop(columns=["ROMol"])
    df_to_save.to_parquet("dataset.parquet")
    ```
2.  **Binary Pickle / Byte String Storage (Preserves 3D coordinates &
    conformers)**:

    ```python
    df["Mol_Bytes"] = df["ROMol"].apply(lambda m: m.ToBinary() if m else None)
    df_to_save = df.drop(columns=["ROMol"])
    df_to_save.to_parquet("conformers.parquet")

    # Reloading from binary:
    df["ROMol"] = df["Mol_Bytes"].apply(lambda b: Chem.Mol(b) if b else None)
    ```

--------------------------------------------------------------------------------

## 6. Presentation & Reporting Protocol: Compact Summaries, Visuals & Artifact Delegation

When communicating cheminformatics and drug design results:

1.  **Keep inline tables concise**: Limit summary markdown tables to 3–5
    representative lead molecules.
2.  **Save Full Datasets & Interactive Grids to Disk**: Save complete CSV tables
    (`<name>.csv`) and interactive HTML grids (`<name>.html`) into the output
    directory and reference them with markdown links
    (`[filename.csv](file:///...)`).
3.  **Focus on Visuals & 2D Depictions**: Use `uv run
    ../rdkit_general/scripts/draw_molecules.py` to generate high-resolution PNGs
    of key leads or reaction schemes (`![caption](/absolute/path/to/mol.png)`).
4.  **Executive Summaries First**: Highlight high-level metrics (library
    attrition %, hit counts, $\text{MW}/\text{cLogP}/\text{TPSA}$ ranges, SAR
    activity cliff insights).
