# Prospective Candidate Scoring & Empirical Physics-ML Benchmarking Guide

## 1. Unified Physics vs. ML Empirical Benchmarking

When a study uses both physical modeling (docking, conformational strain) and
empirical machine learning:

-   **Batch Physical Scoring**: Never restrict physics to a few exemplar
    molecules. Run `molecular_docking_vina.py` across the test split or series
    to compute empirical docking free energy $\Delta G_{\text{dock}}$ and bound
    conformational strain $\Delta E_{\text{strain}}$.
-   **Independent Quantities & No Naive Addition**: Never compute uncalibrated
    addition $\Delta G_{\text{eff}} = \Delta G_{\text{dock}} + \Delta
    E_{\text{strain}}$. Report docking and force field strain independently.
-   **Strict Split Partition Isolation**: When benchmarking against machine
    learning models, always supply `--split-col <col> --eval-split test` (or
    rely on default holdout detection) to prevent evaluating across pooled
    train+test data.
-   **Head-to-Head Side-by-Side Comparison**: Run `physics_ml_comparator.py` to
    evaluate $R^2$, Pearson $r$, Spearman $\rho$, and RMSE for physical
    estimators alongside ML models strictly on the evaluation split.
-   **Multi-Variable Delta-Learning ($\Delta$-ML)**: Train and evaluate hybrid
    models blending physical priors with ML predictions ($y \approx
    w_{\text{ml}} \cdot y_{\text{ml}} + w_{\text{dock}} \cdot (-\Delta
    G_{\text{dock}}) + w_{\text{strain}} \cdot \Delta E_{\text{strain}} + b$)
    via multi-variable regression to test for incremental generalization gains.

### CLI Example

```bash
uv run scripts/physics_ml_comparator.py \
  --csv series_predictions.csv \
  --target-col pIC50 \
  --docking-col dGdock \
  --strain-col dEstrain \
  --ml-cols pred_RidgeCV pred_RandomForest \
  --split-col split \
  --eval-split test \
  --output-csv /tmp/comparator_metrics.csv \
  --output-md /tmp/comparator_report.md
```

## 2. Prospective Candidate Multi-Endpoint Scoring

When proposing novel chemical analogs:

-   **Comprehensive Multi-Endpoint Inference**: Forward inference must be
    executed across every trained model endpoint, not just a single favored
    assay.
-   **Parent-Relative Delta Shifts**: Calculate $\Delta \hat{y} =
    \hat{y}_{\text{candidate}} - \hat{y}_{\text{parent}}$ to measure predicted
    potency gain.
-   **Selectivity & Resistance Margins**: Calculate $\Delta\Delta =
    \Delta\text{Mutant} - \Delta\text{WT}$ to assess resistance evasion or
    selectivity.
-   **Consensus Pareto Prioritization**: Derive objective priority tiers (Tier 1
    High, Tier 2 Moderate, Tier 3 Low) integrating predicted potency, physical
    interactions, and drug-likeness.

### CLI Example

```bash
uv run scripts/predict_potency.py \
  --model-path trained_models.joblib \
  --input-csv proposals.csv \
  --smiles-col SMILES \
  --output-csv /tmp/scored_proposals.csv
```
