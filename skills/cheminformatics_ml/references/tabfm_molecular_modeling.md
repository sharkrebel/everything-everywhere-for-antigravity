# TabFM Zero-Shot Tabular Foundation Modeling for Molecular Data

A comprehensive domain guide for applying Google Research's pretrained Tabular
Foundation Models (`tabfm.TabFMRegressor` and `tabfm.TabFMClassifier`) to
chemical datasets using RDKit and pandas.

--------------------------------------------------------------------------------

## 1. Overview & Conceptual Architecture

TabFM is a foundation model architecture pretrained on diverse tabular datasets,
capable of performing **zero-shot in-context inference** directly on tabular
data without task-specific neural network training, gradient descent, or
hyperparameter optimization.

Chemical datasets are inherently multimodal:

1.  **Continuous Numerical Features**: 2D physical and topological properties
    (e.g. molecular weight, topological polar surface area, calculated logP,
    fractional CSP3, hydrogen bond counts) and Morgan count fingerprints.
2.  **Categorical Structural Features**: Core chemotype topology and ring
    systems captured via Bemis-Murcko scaffold SMILES strings.
3.  **TabFM Synergy**: TabFM natively supports mixed continuous and
    string-categorical columns without requiring manual one-hot encoding, target
    encoding, or category mapping.

--------------------------------------------------------------------------------

## 2. Invariants & Data Hygiene

> [!IMPORTANT] - **C++ Object Sanitization Invariant**: In-memory RDKit
> molecules (`rdkit.Chem.rdchem.Mol` / `ROMol`) are C++ wrapper objects. They
> cannot be serialized by tabular models and will raise exceptions if passed in
> feature matrices. You MUST purge `ROMol` columns and raw SMILES before passing
> DataFrames to TabFM estimators. - **Strict Partition Separation**: When
> evaluating TabFM zero-shot models, fit exclusively on defined training
> partitions ($X_{\text{train}}, y_{\text{train}}$) and evaluate on untouched
> prospective test partitions. Never mix test rows into the in-context training
> context.

--------------------------------------------------------------------------------

## 3. End-to-End Five-Stage Pipeline

```
   ┌──────────────────────────────────────────────────────────┐
   │ 1. DataFrame Ingestion via PandasTools                   │
   │    rdkit.Chem.PandasTools.AddMoleculeColumnToFrame       │
   │    Filter invalid/unparsable ROMol entries               │
   └─────────────────────────────┬────────────────────────────┘
                                 │
   ┌─────────────────────────────▼────────────────────────────┐
   │ 2. Mixed Feature Engineering                             │
   │    - Continuous 2D descriptors (CalcMolDescriptors)      │
   │    - Categorical Bemis-Murcko scaffolds (scaffold SMILES)│
   │    - Optional count fingerprints (Morgan counts)         │
   └─────────────────────────────┬────────────────────────────┘
                                 │
   ┌─────────────────────────────▼────────────────────────────┐
   │ 3. Feature Space Sanitization                            │
   │    Purge C++ ROMol objects, raw SMILES, IDs, metadata    │
   │    Retain only numeric scalars and categorical strings   │
   └─────────────────────────────┬────────────────────────────┘
                                 │
   ┌─────────────────────────────▼────────────────────────────┐
   │ 4. Partition Isolation                                   │
   │    Isolate train partition from prospective test set     │
   │    Zero data leakage or pre-fit feature contamination    │
   └─────────────────────────────┬────────────────────────────┘
                                 │
   ┌─────────────────────────────▼────────────────────────────┐
   │ 5. Zero-Shot TabFM Estimation                            │
   │    Fit TabFMRegressor / TabFMClassifier on (X_train, y)  │
   │    Predict on X_test; compute generalization metrics     │
   └──────────────────────────────────────────────────────────┘
```

--------------------------------------------------------------------------------

## 4. Standard Python Recipes

### 4.1 Molecular Ingestion & Validation

```python
import pandas as pd
from rdkit.Chem import PandasTools

def load_and_validate_molecules(
    csv_path: str,
    smiles_column: str = "smiles",
    mol_column: str = "ROMol",
) -> pd.DataFrame:
    """Loads CSV dataset, converts SMILES to RDKit ROMol, and drops invalid structures."""
    df = pd.read_csv(csv_path)
    if smiles_column not in df.columns:
        raise ValueError(f"SMILES column '{smiles_column}' not found in {csv_path}")

    # Add RDKit Mol column without automatic fingerprint generation
    PandasTools.AddMoleculeColumnToFrame(
        df,
        smilesCol=smiles_column,
        molCol=mol_column,
        includeFingerprints=False,
    )

    # Filter out unparsable SMILES (where ROMol is None / NaN)
    n_initial = len(df)
    df = df[df[mol_column].notna()].reset_index(drop=True)
    n_valid = len(df)
    if n_valid < n_initial:
        print(f"Filtered {n_initial - n_valid} unparsable molecules ({n_valid}/{n_initial} retained).")

    return df
```

### 4.2 Mixed Feature Space Construction

```python
import numpy as np
import pandas as pd
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def extract_mixed_molecular_features(
    df: pd.DataFrame,
    mol_column: str = "ROMol",
    include_descriptors: bool = True,
    include_scaffolds: bool = True,
) -> pd.DataFrame:
    """Extracts continuous 2D physical descriptors and categorical Murcko scaffold strings."""
    df_features = df.copy()

    # 1. Continuous 2D physical descriptors
    if include_descriptors:
        desc_list = [
            (name, func) for name, func in Descriptors.descList
            if not name.startswith("MaxPartialCharge") and not name.startswith("MinPartialCharge")
        ]
        for name, func in desc_list:
            vals = []
            for mol in df_features[mol_column]:
                try:
                    val = float(func(mol))
                    vals.append(val if np.isfinite(val) else np.nan)
                except Exception:
                    vals.append(np.nan)
            df_features[f"desc_{name}"] = vals

    # 2. Categorical Bemis-Murcko scaffold identities
    if include_scaffolds:
        scaffolds = []
        for mol in df_features[mol_column]:
            try:
                scaff = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
                scaffolds.append(str(scaff) if scaff else "acyclic")
            except Exception:
                scaffolds.append("unknown")
        df_features["scaffold_smiles"] = pd.Series(scaffolds, dtype="category")

    return df_features
```

### 4.3 Sanitization & Partition Isolation

```python
from typing import List, Tuple
import pandas as pd

def prepare_feature_partitions(
    df_engineered: pd.DataFrame,
    target_column: str,
    split_column: str,
    train_label: str = "train",
    test_label: str = "test",
    mol_column: str = "ROMol",
    smiles_column: str = "smiles",
    id_columns: List[str] = None,
) -> Tuple[pd.DataFrame, pd.Series, pd.DataFrame, pd.Series]:
    """Sanitizes feature matrix by purging C++ Mol objects and non-feature columns, then splits."""
    if id_columns is None:
        id_columns = ["id", "compound_id", "sample_id", "name"]

    # Identify metadata and non-feature columns to drop
    drop_cols = {mol_column, smiles_column, split_column, target_column}
    for col in df_engineered.columns:
        if col.lower() in [i.lower() for i in id_columns]:
            drop_cols.add(col)

    feature_cols = [c for c in df_engineered.columns if c not in drop_cols]

    # Partition masks
    train_mask = df_engineered[split_column] == train_label
    test_mask = df_engineered[split_column] == test_label

    X_train = df_engineered.loc[train_mask, feature_cols].copy()
    y_train = df_engineered.loc[train_mask, target_column].copy()
    X_test = df_engineered.loc[test_mask, feature_cols].copy()
    y_test = df_engineered.loc[test_mask, target_column].copy()

    # Final verification: ensure NO C++ objects remain in feature columns
    for col in feature_cols:
        sample_val = X_train[col].dropna().iloc[0] if len(X_train[col].dropna()) > 0 else None
        if hasattr(sample_val, "__class__") and "rdkit" in str(sample_val.__class__).lower():
            raise TypeError(f"Column '{col}' contains unpurged RDKit C++ objects: {sample_val.__class__}")

    return X_train, y_train, X_test, y_test
```

### 4.4 Model Fitting & Inference

#### Continuous Targets (Regression)

```python
import tabfm
from tabfm import TabFMRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from scipy.stats import spearmanr

# 1. Load the pre-trained TabFM regression model
tabfm_base_model = tabfm.tabfm_v1_0_0.load(model_type="regression")

# 2. Instantiate TabFMRegressor (optionally with ensemble preset)
regressor = TabFMRegressor(model=tabfm_base_model, random_state=42)
# Alternatively for heavier ensembling:
# regressor = TabFMRegressor.ensemble(model=tabfm_base_model, random_state=42)

# 3. Fit zero-shot model on training partition
regressor.fit(X_train, y_train)

# 4. Predict on held-out test partition
y_pred = regressor.predict(X_test)

# 5. Evaluate generalization metrics
r2 = r2_score(y_test, y_pred)
rho, _ = spearmanr(y_test, y_pred)
mae = mean_absolute_error(y_test, y_pred)

print(f"Held-Out Test Results: R²={r2:.3f}, Spearman ρ={rho:.3f}, MAE={mae:.3f}")
```

#### Discrete / Categorical Targets (Classification)

```python
import tabfm
from tabfm import TabFMClassifier
from sklearn.metrics import roc_auc_score, average_precision_score, balanced_accuracy_score

# 1. Load the pre-trained TabFM classification model
tabfm_base_model = tabfm.tabfm_v1_0_0.load(model_type="classification")

# 2. Instantiate TabFMClassifier
classifier = TabFMClassifier(model=tabfm_base_model, random_state=42)

# 3. Fit zero-shot classifier on training partition
classifier.fit(X_train, y_train)

# 4. Predict probabilities and discrete classes on test partition
y_prob = classifier.predict_proba(X_test)[:, 1]
y_pred = classifier.predict(X_test)

# 5. Evaluate classification performance
roc_auc = roc_auc_score(y_test, y_prob)
pr_auc = average_precision_score(y_test, y_prob)
bacc = balanced_accuracy_score(y_test, y_pred)

print(f"Held-Out Test Results: ROC-AUC={roc_auc:.3f}, PR-AUC={pr_auc:.3f}, BalAcc={bacc:.3f}")
```

### 4.5 Automated Multi-Model Benchmarking via `congeneric_series_ml.py`

To evaluate TabFM alongside scikit-learn baselines (RidgeCV, ElasticNetCV, SVR,
Random Forest, Gradient Boosting) across multiple feature spaces without manual
scripting, use `congeneric_series_ml.py` with `--include_tabfm`:

```bash
uv run scripts/congeneric_series_ml.py \
    --input_csv congeneric_series.csv \
    --target_col pIC50 \
    --include_tabfm \
    --tabfm_estimators 2 \
    --output_csv ml_benchmark_summary.csv \
    --output_predictions ml_predictions.csv \
    --plot_parity ml_parity.png
```

Key arguments:

-   `--include_tabfm`: Activates TabFM zero-shot regression within candidate
    model sweeps. Pre-loads model weights once to eliminate per-fold load
    overhead.
-   `--tabfm_estimators <N>`: Sets number of ensemble members for in-context
    predictions (default: 4; use 2 for fast CPU evaluation).

--------------------------------------------------------------------------------

## 5. Diagnostic & Troubleshooting Guide

-   **`TypeError: Object of type Mol is not JSON serializable`**
    -   **Root Cause**: C++ `ROMol` column left in DataFrame
    -   **Remediation**: Drop the `ROMol` column (`df.drop(columns=['ROMol'])`)
        prior to calling `.fit()` or `.predict()`.
-   **`KeyError: 'model'` or `TypeError: missing 1 required positional argument:
    'model'`**
    -   **Root Cause**: Attempting to instantiate `TabFMRegressor()` without
        passing pre-trained model weights
    -   **Remediation**: Call `tabfm.tabfm_v1_0_0.load(model_type='regression')`
        and pass the returned model object into `TabFMRegressor(model=...)`.
-   **NaN predictions or feature shape mismatch**
    -   **Root Cause**: Unfiltered unparsable SMILES resulting in NaN feature
        rows
    -   **Remediation**: Filter molecules with `df =
        df[df['ROMol'].notna()].reset_index(drop=True)` immediately after
        `AddMoleculeColumnToFrame`.
