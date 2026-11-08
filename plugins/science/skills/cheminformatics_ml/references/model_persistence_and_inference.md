# Model Persistence, Feature Attribution Visuals, and Forward Inference Guide

This reference guide documents the model serialization standards, automated
feature attribution plotting, and forward potency inference protocols for
`cheminformatics-ml`.

--------------------------------------------------------------------------------

## 1. Model Persistence & Production "Max-Data" Retraining

When benchmarking congeneric series via `congeneric_series_ml.py`:

-   **Winning Architecture Selection**: The engine evaluates candidate model
    families (RidgeCV, ElasticNetCV, SVR, Random Forest, Gradient Boosting,
    TabFM) across feature sets (2D physicochemical, 3D shape, Morgan circular
    fingerprints) using leakage-free CV or nested holdouts.
-   **Max-Data Retraining Invariant**: In discovery projects, model evaluation
    determines the optimal inductive bias and regularization strength. Once the
    winning architecture is identified, the engine fits a fresh preprocessor and
    retrains that exact architecture on **100% of the available labeled data**
    (train + holdout combined). Both the benchmark and the production model come
    from the same architecture factory, so they can never diverge.
-   **Self-Contained Model Bundles (`.joblib`)**: Persistence is **opt-in**:
    pass `--model_dir <dir>` to write bundles. Without it, no model files are
    created and the benchmark leaves the working directory clean.
    -   `<target>_maxdata_model.joblib`: retrained on 100% of labeled data, for
        production prospective screening.
    -   `<target>_best_cv_model.joblib`: the winning benchmark model, fitted on
        the training split only, preserving cross-validation statistics.
-   **Bundle Architecture**:

    ```python
    bundle = {
        "model": trained_estimator,          # Fitted model
        "preprocessor": fitted_preprocessor, # Impute + VarianceThreshold + StandardScaler
        "feature_set": best_feature_set,     # e.g. "2D Descriptors", "Morgan Circular FP"
        "feature_cols": list_of_raw_columns, # Input feature columns required
        "kept_feature_names": kept_names,    # Features surviving variance filtering
        "target": target_name,               # Endpoint name (e.g. "pIC50", "pKi")
        "model_name": algorithm_name,        # e.g. "ElasticNetCV", "Random Forest"
        "is_fingerprint": bool,
        "is_p_scale": bool,                  # True only if the endpoint is a negative
                                             # log molar affinity AND the observed labels
                                             # are consistent with one. Inference clamps
                                             # p-scale predictions to [0, 15]; recording
                                             # this here stops the consumer guessing from
                                             # the target's name alone and clamping, say,
                                             # a nanomolar IC50 into oblivion.
        "training_range": [y_min, y_max],    # Observed label range, or None. Drives the
                                             # Extrapolated_Beyond_Train output column.
        "training_stats": {
            "n_samples": n_obs,
            "cv_r2": cv_r2,
            "test_r2": test_r2,
        },
        "training_smiles": list_of_training_smiles, # For Applicability Domain checks
    }
    ```
-   **Portability**: bundles are written so that a plain `joblib.load` works in
    any process. When the training script runs as `__main__`, its custom classes
    would otherwise pickle as `__main__.BoundedRegressor` and fail to load
    elsewhere; the writer re-points them at the real module name first.

### CLI Execution

```bash
uv run scripts/congeneric_series_ml.py \
  --input_csv series.csv \
  --smiles_col SMILES \
  --target_col pKi \
  --model_dir models \
  --plot_feature_importance
```

--------------------------------------------------------------------------------

## 2. Visual Feature Attribution Protocol (No ASCII Tables)

-   **Mandatory High-Resolution Plots**: Feature importance / attribution must
    always be presented as an embedded visual plot (`![Feature
    Attribution](...)`), never as raw ASCII text tables inside code blocks.
-   **Universal Attribution Engine**: the method is chosen from what the fitted
    model actually exposes:
    -   **Linear Models (RidgeCV, ElasticNetCV)**: standardized regression
        coefficients ($\beta$), colored by directional effect (positive in blue,
        negative in red).
    -   **Tree Models (Random Forest, Gradient Boosting)**: mean decrease in
        impurity (MDI / Gini importance).
    -   **Kernel & Non-Linear Estimators (SVR RBF, TabFM)**: permutation feature
        importance, measuring $R^2$ degradation on the training data when each
        feature is shuffled.
-   **Opt-in Execution**: pass `--plot_feature_importance` with no value to
    write `feature_importance_<target>.png` in the current directory, or supply
    an explicit path. With multiple `--target_col` endpoints an explicit path is
    suffixed per target so plots do not overwrite each other. The chosen
    attribution method is printed in the axis label and plot title, and the
    markdown embed line is emitted on stdout. If no attribution can be computed
    for the winning model, a warning is printed to stderr and no plot is written
    — an empty plot is never produced.

--------------------------------------------------------------------------------

## 3. Forward Potency Scoring (`predict_potency.py`)

To score novel designed analogs, virtual screening hits, or prospective series
without refitting models:

### CLI Execution

```bash
uv run scripts/predict_potency.py \
  --model_path models/pKi_maxdata_model.joblib \
  --input_csv candidate_designs.csv \
  --smiles_col SMILES \
  --output_csv scored_candidates.csv \
  --min_tanimoto 0.40
```

### Protocol & Safeguards:

1.  **Authoritative Featurization**: Extracts features using the exact pipeline
    (`compute_molecular_features`) used during model training.
2.  **Pipeline Preprocessing**: Applies the fitted preprocessor from the bundle
    without recomputing training-fold statistics.
3.  **Physical Bounds Safeguards**: When the bundle records `is_p_scale: true`,
    enforces the thermodynamic / aqueous diffusion window $[0.0, 15.0]$ and
    clamps sub-molecular fragments ($\leq 2$ heavy atoms) to inactive baselines.
    Non-p-scale endpoints are left completely unconstrained, so a nanomolar or
    percent-inhibition endpoint is never silently crushed into the p-scale
    window.
4.  **Applicability Domain (AD) Verification**:
    -   For every query molecule, calculates maximum Tanimoto similarity
        ($S_{\max}$) against all molecules in `training_smiles` using Morgan
        circular fingerprints (radius=2, 2048 bits).
    -   If $S_{\max} < 0.40$ (or `--min_tanimoto`), flags the compound as
        out-of-domain in the output CSV and prints an explicit warning.
5.  **Label Extrapolation Flagging**: when the bundle carries `training_range`,
    adds `Extrapolated_Beyond_Train` and `Training_Range` columns marking
    predictions that fall outside the observed label span. This is orthogonal to
    the AD check: a compound can be structurally in-domain yet predicted at a
    potency never observed in training.
6.  **Output Deliverables**: Generates `scored_candidates.csv` containing
    predicted potencies, maximum Tanimoto similarities, applicability domain
    boolean flags, and — when available — the extrapolation flags above.
