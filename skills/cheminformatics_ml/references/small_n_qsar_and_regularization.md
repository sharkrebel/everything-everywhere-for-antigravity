# Small-$N$ Machine Learning & Regularized Modeling on Congeneric Series

This guide provides principles and an authoritative Python reference template
for training, evaluating, and attributing machine learning models on small
congeneric series ($N < 100$, $P \gg N$).

--------------------------------------------------------------------------------

## 1. Fundamental Principles for Small-$N$ Cheminformatics ML

When modeling small congeneric lead series ($N = 20\text{--}80$):

1.  **High-Capacity Tree Overfitting**:
    -   Random Forest and Gradient Boosting models have high capacity. When
        trained on high-dimensional sparse representations (e.g.
        128/1024/2048-bit Morgan Fingerprints), unregularized tree estimators
        memorize training bit noise, resulting in poor or negative
        cross-validated $R^2$.
2.  **Adaptive Regularization Across Dimensions**:
    -   A static regularization penalty (e.g. fixed $\alpha = 1.0$) behaves
        completely differently on a 7-dimensional physical descriptor set than
        on a 147-dimensional hybrid representation.
    -   Always use cross-validated adaptive regularization (`RidgeCV`,
        `ElasticNetCV`, `LassoCV`) with log-spaced $\alpha$ grids ($\alpha \in
        [10^{-3}, 10^3]$).
3.  **Transparent Evaluation Across Model Families**:
    -   Compare linear regularized baselines (`RidgeCV`), kernel methods
        (`SVR-RBF`), and tree ensembles (`RandomForest`, `GradientBoosting`).
    -   Report the full benchmark table transparently rather than selectively
        presenting only the top-performing rows.
4.  **Disciplined Feature Attribution**:
    -   Never extract feature importances (e.g. Gini impurity or tree
        permutation) from models that have poor or unvalidated cross-validated
        performance (e.g. LOOCV $R^2 < 0.35$ or below the linear baseline).
    -   For regularized linear models, extract standardized regression
        coefficients ($\beta_j$) to evaluate property contributions.

--------------------------------------------------------------------------------

## 2. Production Execution: Congeneric Series ML Benchmark

To prevent data leakage, evaluate independent holdout test sets, and compare
models across representations, always run the production engine
`scripts/congeneric_series_ml.py`:

```bash
uv run scripts/congeneric_series_ml.py \
    --input_csv scratch/curated_series.csv \
    --smiles_col smiles \
    --target_col "pKi_WT,pKi_Mutant,pIC50" \
    --extra_cols "Docking & Strain:docking_score,strain_energy,delta_g_eff" \
    --extra_cols "Pharmacophore:is_macrocycle,c3_nitrile" \
    --output_csv scratch/ml_benchmark_summary.csv \
    --output_predictions scratch/predictions_cv_and_test.csv \
    --plot_parity scratch/multitarget_parity_grid.png \
    --plot_coefficients scratch/feature_weights.png
```

### Key CLI Capabilities & Parameters:

*   **Adaptive Splitting Protocol (Mandatory Default)**:
    -   **Sheridan Temporal Split (`--time-cutoff`)**: Activated when valid
        non-boolean chronological metadata (`document_year`, `publication_date`,
        `filing_date`, or `priority_date`) is available across $\ge 2$ distinct
        dates. Boolean values (`False` / `True`) are rejected (`False < True` is
        never a valid temporal split), missing dates are backfilled from true
        document/patent metadata before partitioning, and the split is always
        cut strictly at `--time-cutoff <year_or_YYYY-MM-DD>` (`train <= cutoff`,
        `test > cutoff`) so tied dates never straddle partitions and
        $\max(t_{\text{train}}) < \min(t_{\text{test}})$ (zero time overlap).
    -   **Nested Cross-Validation ($N < 30$)**: Outer 5-fold CV splits data into
        untouched test folds, while inner LOOCV tunes hyperparameters on
        training folds (Cawley & Talbot 2010 protocol). Delivers 100%
        out-of-sample test coverage ($N_{\text{test}} = N$) with zero selection
        bias.
    -   **Bemis-Murcko Scaffold Split ($N \ge 30$)**: Sizes the test set
        adaptively ($N_{\text{test}} = \max(6, \sqrt{N})$) to evaluate scaffold
        hopping with zero core overlap.
*   **`--extra_cols "<Category Name>:col1,col2,..."` (Repeatable)**: Adds one or
    more custom feature categories (e.g., `Pharmacophore`, `Macrocyclization`,
    `Docking Physics`). A category name is **required** for every custom block
    (via the `<Category Name>:` prefix or `--extra_name`). Pass `--extra_cols`
    multiple times (or separate blocks with `;`) when custom columns come from
    different categories; the script automatically benchmarks each category
    individually (`<Name> Features`, `Hybrid (2D/3D + <Name>)`) as well as their
    joint combination (`Combined (<Name1> + <Name2>)`,
    `Hybrid (2D/3D + <Name1> + <Name2>)`, and `Full Hybrid`).
*   **`--plot_parity <path.png>`**: Generates a publication-grade multi-target
    parity grid with separate markers for `Train (CV)` and `Holdout Test`
    points.
*   **`--plot_coefficients <path.png>`**: Extracts and ranks standardized
    regression coefficients ($\beta_j$) for the best linear regularized model.

> [!CAUTION] **ABSOLUTE BAN ON PURE-CV SHORTCUTS**: Never write custom scripts
> that execute pure Leave-One-Out CV on an entire dataset without an independent
> holdout test set or nested outer test partition. Reporting only
> cross-validation $Q^2$ without a test set is scientifically invalid and
> prohibited.

--------------------------------------------------------------------------------

## 3. Root-Cause Diagnosis & Iterative Feature Engineering Protocol ($Q^2 < 0.35$)

When standard 2D topological descriptors or Morgan fingerprints yield poor
predictive performance ($Q^2 \le 0.10$ or negative), agents must NOT passively
accept the poor score. Execute a structured 3-stage troubleshooting protocol:

### Stage 1: Structural Root-Cause Diagnosis

Evaluate why 2D / topological features fail on the series:

1.  **Stereoisomeric Potency Cliffs**: Chiral inversions (e.g. eutomer vs
    distomer) share identical 2D connectivity graphs and $\ge 99\%$ identical
    ECFP4 bit representations despite substantial potency differences.
2.  **Macrocyclic Preorganization & Rigidification**: 2D descriptors treat
    flexible acyclics and rigid macrocycles similarly, ignoring frozen rotor
    entropy gains ($-T\Delta S_{\text{rot}}$) and bound conformational strain
    ($\Delta E_{\text{strain}}$).
3.  **Sparse Bit Noise Overfitting**: On small $N$ ($N \le 30$), 512/1024-bit
    fingerprint vectors have $P \gg N$, causing non-linear models (SVR, RF,
    GBDT) to memorize sparse bit combinations.

### Stage 2: Hypothesis-Driven Feature Engineering

1.  **Incorporate 3D Physical & Thermodynamic Descriptors**:
    -   MMFF94s bound conformational strain ($\Delta E_{\text{strain}}$)
    -   Rotational entropy penalty ($-T\Delta S_{\text{rot}} = -0.65 \times
        \text{RotB}$)
    -   Receptor contact docking scores ($\Delta G_{\text{dock}}$)
    -   Composite effective free energy ($\Delta G_{\text{eff}} = \Delta
        G_{\text{dock}} + \Delta E_{\text{strain}} - T\Delta S_{\text{rot}}$)
2.  **Stereochemical & Congeneric Indicators**:
    -   Chiral center flags, dihedral alignment match metrics, and matched
        molecular pair $\Delta$-substituent properties.
3.  **Dimensionality Filtering**:
    -   Use `VarianceThreshold(0.01)` and regularized linear models (`RidgeCV`,
        `ElasticNetCV`) to prune uninformative features inside cross-validation
        folds.

### Stage 3: Iteration Progression Reporting

Document the optimization progression using the exact computed metrics from
execution:

-   ****Attempt 1 (Baseline)****
    -   **Feature Representation & Model Strategy**: 2D Descriptors (RidgeCV)
    -   **LOOCV $Q^2$**: *[Computed $Q^2$]*
    -   **Pearson $r$**: *[Computed $r$]*
    -   **Status / Diagnosis**: Assess topological baseline and sensitivity to
        2D properties
-   ****Attempt 1 (Baseline)****
    -   **Feature Representation & Model Strategy**: Morgan Fingerprints ECFP4
        (RidgeCV)
    -   **LOOCV $Q^2$**: *[Computed $Q^2$]*
    -   **Pearson $r$**: *[Computed $r$]*
    -   **Status / Diagnosis**: Assess bit-level substructure patterns and
        sparsity issues
-   ****Attempt 2 (Iteration 1)****
    -   **Feature Representation & Model Strategy**: Pruned 2D + 3D Shape
        Descriptors
    -   **LOOCV $Q^2$**: *[Computed $Q^2$]*
    -   **Pearson $r$**: *[Computed $r$]*
    -   **Status / Diagnosis**: Test whether overall molecular shape/dimensions
        improve accuracy
-   ****Attempt 3 (Iteration 2)****
    -   **Feature Representation & Model Strategy**: Physics-Based ($\Delta
        G_{\text{eff}}, \Delta E_{\text{strain}}, \text{RotB}$)
    -   **LOOCV $Q^2$**: *[Computed $Q^2$]*
    -   **Pearson $r$**: *[Computed $r$]*
    -   **Status / Diagnosis**: Test physical binding hypothesis (docking +
        preorganization)

--------------------------------------------------------------------------------

## 4. Validation vs. Holdout Testing & Reporting Integrity

Understanding the distinction between internal cross-validation and external
holdout testing is essential for rigorous cheminformatics research:

### 1. Internal Cross-Validation ($Q^2_{\text{LOO}}$ / K-Fold CV)

*   **Purpose**: Evaluates algorithmic stability, tunes hyperparameters
    ($\alpha$), and guides representation selection on the training set.
*   **Terminology**: Cross-validation folds are **validation folds**, not an
    independent test set. Never report LOOCV $Q^2$ as "holdout test set
    performance".
*   **Selection Bias**: When benchmarking 30+ combinations of representations
    and algorithms on the same dataset, the maximum observed $Q^2$ is subject to
    multiple hypothesis testing (optimistic bias). Always report the full
    benchmark table rather than selectively quoting only the highest number.

### 2. Independent Holdout Test Sets ($R^2_{\text{test}}$)

*   **Purpose**: Evaluates true prospective generalization to unseen chemical
    series using Bemis-Murcko scaffold splitting or temporal splitting.
*   **The Small-$N$ Holdout Reality**: In small series ($N < 30$), holding out
    20% leaves only 5 to 6 compounds in the test set. If the test set compounds
    have a narrow potency range (e.g. all within $0.5 \log$ units) or contain
    unmatched structural modifications, $R^2_{\text{test}}$ will frequently be
    zero or negative ($< 0$).
*   **Absolute Prohibition on Suppressing Negative Test Results**:
    -   If an explicit train/test split is evaluated, **you must report the
        exact test set metrics ($R^2_{\text{test}}$,
        $\text{RMSE}_{\text{test}}$, Pearson $r$) faithfully**.
    -   Never hide negative test set performance, omit test columns from summary
        tables, or pretend that cross-validation folds substitute for an
        independent test set.
    -   Transparently discuss the scientific reasons for test set failure (e.g.
        narrow target variance in the test set, out-of-domain stereocenters, or
        lack of training examples for key substituents).

### 3. Leakage-Free Preprocessing Rules

*   **Imputation**: Calculate column means strictly on the training fold; never
    use the full dataset mean to impute missing values.
*   **Variance Filtering**: Fit `VarianceThreshold` strictly on the training
    fold; apply `.transform()` to validation or test samples.
*   **Standardization**: Fit `StandardScaler` strictly on the training fold;
    apply `.transform()` to validation or test samples.
*   **Zero Contamination**: Test samples must never be seen during any stage of
    feature selection, scaling, or model fitting.

### 4. Partitioning Strategies: Random, Scaffold, and Temporal (Time) Splits

*   **Random Split (PROHIBITED - DATA LEAKAGE RISK)**: Random train/test
    splitting on congeneric small-molecule series leaks identical Murcko
    scaffolds and congeneric analogues between train and test sets, artificially
    inflating test metrics ($R^2$, RMSE) by evaluating nearest-neighbor
    interpolation rather than prospective generalization. Use `--split_method
    auto` (scaffold, temporal, or nested CV).
*   **Bemis-Murcko Scaffold Split**: Clusters molecules by Bemis-Murcko core
    frameworks. Evaluates scaffold-hopping and structural extrapolation with
    zero core scaffold overlap.
*   **Temporal / Time Split**:
    -   **The Sheridan Prospective Paradigm**: Real-world drug discovery is
        inherently temporal—models trained on historical lead compounds must
        predict future designs. Evaluating on a chronological holdout test set
        tests true prospective utility.
    -   **ChEMBL & Patent Time Information**: ChEMBL and Google Patents provide
        `document_year`, `publication_date`, `filing_date`, and `priority_date`
        records. When source provenance is available, missing or placeholder
        dates are backfilled with true chronological dates before splitting, and
        partitioning is performed strictly with `--time-cutoff` (e.g. train $\le
        2013$, test $> 2013$ via `--time-cutoff 2013`, resolving plain-year
        cutoffs to Dec 31 when sub-year ISO dates are present) so there is zero
        time overlap. Boolean flags (`False` / `True`) and out-of-range
        sentinels (`0`, `-1`, `9999`) are rejected (`False < True` is never a
        valid temporal split). Any residual rows that remain undated (`NaN`)
        after backfilling trigger a warning and are placed in `train` only
        (never in `test`); automatic `--split-method auto` requires $\ge 80\%$
        dated coverage before choosing a temporal split.
    -   **Single-Document Datasets**: Within a single publication, all compounds
        share the same `document_year`. In this context, Bemis-Murcko scaffold
        splitting or nested cross-validation should be used rather than forcing
        a degenerate temporal split.

### 5. Principled Adaptive Partitioning & Nested Cross-Validation (Default Engine)

*   **Why 80/20 is Unprincipled**: Arbitrary 80/20 splits are scale-blind. On
    $N = 16$, 20% leaves only 3 test compounds ($df = 1$, wildly unstable
    metrics).
*   **The Adaptive Engine**:
    1.  **Sheridan Temporal Split (`--time-cutoff`)**: Triggered when valid
        non-boolean chronological metadata (`document_year`, `publication_date`,
        etc.) is present across $\ge 2$ distinct dates, snapping to a strict
        `--time-cutoff` so $\max(t_{\text{train}}) < \min(t_{\text{test}})$.
    2.  **Nested Cross-Validation ($N < 30$)**: Based on Cawley & Talbot (2010,
        *JMLR*), an outer $K$-fold loop partitions data into untouched test
        folds, while an inner LOOCV loop tunes hyperparameters strictly on
        training folds. Aggregating outer test predictions provides **100%
        out-of-sample coverage ($N_{\text{test}} = N$)** with zero selection
        bias and complete degrees of freedom.
    3.  **Adaptive Scaffold Sizing ($N \ge 30$)**: Sizes the test set adaptively
        via $N_{\text{test}} = \max(6, \sqrt{N})$ to guarantee statistical test
        power while preserving training diversity.

### 6. Production CLI Recipe: Automated Congeneric Series Benchmarking

Use `congeneric_series_ml.py` to benchmark regularized models across
representations and generate parity plots:

```bash
uv run scripts/congeneric_series_ml.py \
    --input_csv scratch/curated_series.csv \
    --smiles_col smiles \
    --target_col "pKi,pIC50" \
    --extra_cols "Physics:docking_score_kcal_mol,strain_energy_kcal_mol,effective_delta_g_kcal_mol" \
    --output_csv scratch/ml_benchmark.csv \
    --output_predictions scratch/predictions_cv_and_test.csv \
    --plot_parity scratch/multitarget_ml_parity_grid.png \
    --plot_coefficients scratch/feature_weights.png
```

-   Executes adaptive splitting: Sheridan temporal split if `document_year`
    present, nested 5-fold CV for $N < 30$, or adaptive Bemis-Murcko scaffold
    split for $N \ge 30$.
-   Evaluates across Morgan fingerprints, 2D physicochemical descriptors, 3D
    physics features, and hybrid sets with zero data leakage.

### 7. Matched Molecular Pair Analysis (MMPA) & Macrocyclization Matching

Run `mmpa_sar_analysis.py` to identify single-transformation activity cliffs and
matched pairs:

```bash
uv run scripts/mmpa_sar_analysis.py \
    --input_csv scratch/curated_series.csv \
    --smiles_col smiles \
    --activity_col pIC50 \
    --output_mmpa_csv scratch/mmpa_pairs.csv \
    --output_sali_csv scratch/sali_cliffs.csv
```

*   **Strict Size-Bounded Invariants**:
    -   Core scaffold must have $\ge 12$ heavy atoms and $\ge 50\%$ of the
        parent molecule's heavy atom count.
    -   Substituent transformation must have $\le 10$ heavy atoms and $\le 35\%$
        of the parent molecule.
    -   Directional single-cut assignment: larger fragment is strictly the core;
        smaller is strictly the substituent. Rejects trivial 1-atom cores
        (`[*:1]C`, `[*:1]N`, `[*:1]Cl`).
*   **Native Macrocyclization Matching**:
    -   Use `--match_macrocycles` to pair open-chain acyclic precursors with
        macrocyclic leads (ring size $\ge 12$) sharing a conserved
        core/pharmacophore.
    -   Computes conformational rigidification ($\Delta \text{RotB}$), brain
        penetration optimization ($\Delta \text{CNS MPO}$), and fold potency
        gain.
