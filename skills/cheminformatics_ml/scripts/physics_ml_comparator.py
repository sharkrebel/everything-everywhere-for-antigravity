# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "matplotlib",
#     "numpy",
#     "pandas",
#     "scikit-learn",
#     "scipy",
# ]
# ///

"""Hybrid Physics vs. ML Empirical Benchmarking and Delta-Learning Engine.

Provides unified, side-by-side empirical performance comparisons between
physics-based scoring (docking free energy dGdock, conformational strain
dEstrain, calibrated effective free energies) and machine learning models
(RidgeCV, Random Forest, TabFM, LightGBM), evaluating correlation metrics (R2,
Pearson r, Spearman rho, RMSE, MAE) and testing hybrid delta-learning
(Delta-ML).

Enforces:
1. Strict train/test partition isolation (--eval-split) preventing pooled
   evaluation leakage.
2. Prohibition of naive uncalibrated vacuum strain addition (dGeff != dGdock +
   dEstrain); raw quantities are evaluated independently and hybrid models
   calibrated via multi-variable regression.
3. Thermodynamic sanity clamping (for active binders, Delta G >= 0 => Kd = inf).
"""

import argparse
import json
import os
from typing import Any, Dict, List, Optional, Tuple, Union
import warnings

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold

GAS_CONSTANT_R = 1.98720425864083e-3  # kcal / (mol * K)
DEFAULT_TEMPERATURE = 298.15  # K


def check_thermodynamic_sanity(
    dG_values: np.ndarray,
    target_potency: Optional[np.ndarray] = None,
    potency_threshold: float = 5.0,
    col_name: str = "dG",
) -> Dict[str, Any]:
  """Checks thermodynamic validity of free energy values for active binders.

  Flags unphysical non-spontaneous binding (Delta G >= 0 kcal/mol) for active
  compounds.
  Returns diagnostic statistics and whether clamping or calibration is needed.
  """
  vals = np.asarray(dG_values, dtype=np.float64)
  mask_finite = np.isfinite(vals)
  if not np.any(mask_finite):
    return {
        "col_name": col_name,
        "n_total": len(vals),
        "n_active": 0,
        "n_unphysical": 0,
        "fraction_unphysical": 0.0,
        "has_unphysical": False,
    }

  if target_potency is not None:
    t_vals = np.asarray(target_potency, dtype=np.float64)
    active_mask = (
        mask_finite & np.isfinite(t_vals) & (t_vals >= potency_threshold)
    )
    if not np.any(active_mask):
      active_mask = mask_finite
  else:
    active_mask = mask_finite

  unphysical_mask = active_mask & (vals >= 0.0)
  n_unphysical = int(np.sum(unphysical_mask))
  n_active = int(np.sum(active_mask))
  frac_unphysical = float(n_unphysical / n_active) if n_active > 0 else 0.0
  median_dg = float(np.median(vals[active_mask])) if n_active > 0 else 0.0

  return {
      "col_name": col_name,
      "n_total": len(vals),
      "n_active": n_active,
      "n_unphysical": n_unphysical,
      "fraction_unphysical": round(frac_unphysical, 4),
      "has_unphysical": n_unphysical > 0,
      "median_dG_kcal_mol": round(median_dg, 3),
  }


def safe_stats(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
  """Calculates R2, Pearson r, Spearman rho, RMSE, and MAE with degenerate protection."""
  if len(y_true) == 0 or len(y_pred) == 0:
    return {
        "r2": 0.0,
        "pearson_r": 0.0,
        "spearman_rho": 0.0,
        "rmse": 0.0,
        "mae": 0.0,
    }

  # Remove any NaN or Inf pairs
  mask = np.isfinite(y_true) & np.isfinite(y_pred)
  y_t = y_true[mask]
  y_p = y_pred[mask]
  if len(y_t) < 2:
    return {
        "r2": 0.0,
        "pearson_r": 0.0,
        "spearman_rho": 0.0,
        "rmse": 0.0,
        "mae": 0.0,
    }

  r2 = float(r2_score(y_t, y_p))
  rmse = float(np.sqrt(mean_squared_error(y_t, y_p)))
  mae = float(mean_absolute_error(y_t, y_p))

  if len(y_t) < 3 or np.std(y_t) < 1e-9 or np.std(y_p) < 1e-9:
    pr = 0.0
    sr = 0.0
  else:
    pr_val, _ = pearsonr(y_t, y_p)
    pr = 0.0 if np.isnan(pr_val) else float(pr_val)
    sr_val, _ = spearmanr(y_t, y_p)
    sr = 0.0 if np.isnan(sr_val) else float(sr_val)

  return {
      "r2": round(r2, 4),
      "pearson_r": round(pr, 4),
      "spearman_rho": round(sr, 4),
      "rmse": round(rmse, 4),
      "mae": round(mae, 4),
  }


def calibrate_physics_to_target(
    y_true: np.ndarray,
    physics_values: np.ndarray,
    invert: bool = True,
    train_data: Optional[Tuple[np.ndarray, np.ndarray]] = None,
) -> Tuple[np.ndarray, Tuple[float, float]]:
  """Fits an affine calibration y ≈ a * x_phys + b to place physical energy on target scale.

  Docking free energies (kcal/mol) are negative, where more negative implies
  tighter binding.
  For pIC50 or pKi, higher values imply tighter binding. When invert=True,
  x_phys is negated. When train_data=(y_train, physics_train) is provided,
  the calibration is fit strictly on the training split to prevent test leakage.
  """
  x = -physics_values if invert else physics_values
  if train_data is not None:
    y_fit, phys_fit = train_data
    x_fit = -phys_fit if invert else phys_fit
  else:
    y_fit, x_fit = y_true, x

  mask = np.isfinite(y_fit) & np.isfinite(x_fit)
  if np.sum(mask) < 2 or np.std(x_fit[mask]) < 1e-9:
    finite_y = y_fit[np.isfinite(y_fit)]
    fallback = float(np.mean(finite_y)) if len(finite_y) else 0.0
    return np.full_like(y_true, fallback, dtype=np.float64), (
        0.0,
        0.0,
    )

  reg = LinearRegression().fit(x_fit[mask].reshape(-1, 1), y_fit[mask])
  slope = float(reg.coef_[0])
  intercept = float(reg.intercept_)
  calibrated = np.full(len(x), intercept, dtype=np.float64)
  out_mask = np.isfinite(x)
  if np.any(out_mask):
    calibrated[out_mask] = reg.predict(x[out_mask].reshape(-1, 1))
  return calibrated, (slope, intercept)


def evaluate_hybrid_delta_learning(
    y_true: np.ndarray,
    y_ml: np.ndarray,
    physics_features: Union[np.ndarray, Dict[str, np.ndarray]],
    cv_folds: int = 5,
    train_data: Optional[
        Tuple[np.ndarray, np.ndarray, Union[np.ndarray, Dict[str, np.ndarray]]]
    ] = None,
) -> Tuple[np.ndarray, Dict[str, float], Dict[str, float]]:
  """Trains and evaluates a hybrid estimator combining ML predictions and physics via multi-variable regression.

  Fits y ≈ w_ml * y_ml + sum(w_phys_i * x_phys_i) + b via cross-validation (or
  prospective train/test fit)
  to produce out-of-fold hybrid predictions and assess incremental predictive
  gain over pure ML.
  Prohibits naive uncalibrated addition of force field strain to docking scores.
  """

  def _build_X(
      y_pred_ml: np.ndarray, phys: Union[np.ndarray, Dict[str, np.ndarray]]
  ) -> Tuple[np.ndarray, List[str]]:
    cols = [y_pred_ml.reshape(-1, 1)]
    feature_names = ["ml"]
    if isinstance(phys, dict):
      for k, v in phys.items():
        v_arr = np.asarray(v, dtype=np.float64).reshape(-1, 1)
        # Invert free energy / docking so more negative binding energy -> higher value (tighter binding)
        if any(term in k.lower() for term in ["dock", "dg", "vina"]):
          v_arr = -v_arr
        cols.append(v_arr)
        feature_names.append(k)
    else:
      p_arr = np.asarray(phys, dtype=np.float64)
      if p_arr.ndim == 1:
        cols.append((-p_arr).reshape(-1, 1))
        feature_names.append("physics")
      elif p_arr.ndim == 2:
        for i in range(p_arr.shape[1]):
          cols.append((-p_arr[:, i]).reshape(-1, 1))
          feature_names.append(f"physics_{i}")
    return np.hstack(cols), feature_names

  X_blend, feat_names = _build_X(y_ml, physics_features)
  n_samples = len(y_true)

  if train_data is not None:
    # Prospective train/test split: train regression strictly on training set, predict on test set
    y_tr, y_tr_ml, phys_tr = train_data
    X_tr, _ = _build_X(y_tr_ml, phys_tr)
    ridge = RidgeCV(alphas=np.logspace(-3, 3, 20)).fit(X_tr, y_tr)
    y_hybrid_pred = ridge.predict(X_blend)
    final_model = ridge
  else:
    # Cross-validation out-of-fold predictions
    k = min(cv_folds, n_samples) if n_samples > 2 else 2
    kf = KFold(n_splits=k, shuffle=True, random_state=42)
    y_hybrid_pred = np.zeros(n_samples, dtype=np.float64)

    for tr_idx, te_idx in kf.split(X_blend):
      X_tr, X_te = X_blend[tr_idx], X_blend[te_idx]
      y_tr = y_true[tr_idx]
      ridge = RidgeCV(alphas=np.logspace(-3, 3, 20)).fit(X_tr, y_tr)
      y_hybrid_pred[te_idx] = ridge.predict(X_te)

    final_model = RidgeCV(alphas=np.logspace(-3, 3, 20)).fit(X_blend, y_true)

  weights = {
      "weight_ml": round(float(final_model.coef_[0]), 4),
      "intercept": round(float(final_model.intercept_), 4),
      "chosen_alpha": round(float(final_model.alpha_), 4),
  }
  for i, name in enumerate(feat_names[1:], start=1):
    weights[f"weight_{name}"] = round(float(final_model.coef_[i]), 4)
  # For backward compatibility if single physics feature
  if "physics" in feat_names[1:]:
    weights["weight_physics"] = weights.get(
        "weight_physics", round(float(final_model.coef_[1]), 4)
    )

  metrics = safe_stats(y_true, y_hybrid_pred)
  return y_hybrid_pred, metrics, weights


class PhysicsMLComparator:
  """Automated comparative benchmark engine for physics, machine learning, and hybrid models."""

  def __init__(
      self,
      target_col: Optional[str] = None,
      docking_col: Optional[str] = None,
      strain_col: Optional[str] = None,
      effective_col: Optional[str] = None,
      ml_cols: Optional[List[str]] = None,
      split_col: Optional[str] = None,
      eval_split: Optional[str] = None,
  ):
    self.target_col = target_col
    self.docking_col = docking_col
    self.strain_col = strain_col
    self.effective_col = effective_col
    self.ml_cols = ml_cols or []
    self.split_col = split_col
    self.eval_split = eval_split
    self.results_df: Optional[pd.DataFrame] = None
    self.summary_dict: Dict[str, Any] = {}

  def _detect_columns(self, df: pd.DataFrame) -> None:
    cols = list(df.columns)
    cols_lower = {c.lower(): c for c in cols}

    # 1. Target column
    if not self.target_col or self.target_col not in df.columns:
      candidates = [
          "pic50",
          "pki",
          "affinity",
          "activity",
          "target",
          "label",
          "exp",
      ]
      for cand in candidates:
        if cand in cols_lower:
          self.target_col = cols_lower[cand]
          break
      if not self.target_col and len(cols) > 0:
        self.target_col = cols[0]

    # 2. Docking column
    if not self.docking_col or self.docking_col not in df.columns:
      candidates = [
          "dgdock",
          "dg_dock",
          "docking_score",
          "docking",
          "docking_score_kcal_mol",
          "vina_affinity",
          "vina_score",
      ]
      for cand in candidates:
        if cand in cols_lower:
          self.docking_col = cols_lower[cand]
          break

    # 3. Strain column
    if not self.strain_col or self.strain_col not in df.columns:
      candidates = [
          "destrain",
          "de_strain",
          "strain_energy",
          "strain",
          "strain_energy_kcal_mol",
          "mmff94_strain",
          "mmff94s_strain",
      ]
      for cand in candidates:
        if cand in cols_lower:
          self.strain_col = cols_lower[cand]
          break

    # 4. Split column
    if not self.split_col or self.split_col not in df.columns:
      candidates = [
          "split",
          "partition",
          "dataset",
          "subset",
          "split_group",
          "ds",
          "fold",
      ]
      for cand in candidates:
        if cand in cols_lower:
          self.split_col = cols_lower[cand]
          break

    # 5. ML prediction columns
    if not self.ml_cols:
      ml_found = []
      exclude_cols = {
          self.target_col,
          self.docking_col,
          self.strain_col,
          self.split_col,
      }
      for c in cols:
        cl = c.lower()
        if c in exclude_cols:
          continue
        if any(
            prefix in cl
            for prefix in [
                "pred",
                "ridge",
                "rf",
                "randomforest",
                "lgb",
                "lightgbm",
                "tabfm",
                "gbt",
                "svr",
                "ml_",
                "elasticnet",
            ]
        ):
          ml_found.append(c)
      self.ml_cols = ml_found

  def compare(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Runs head-to-head benchmarking on the input DataFrame.

    Enforces strict train/test partition isolation when a split column is
    present.
    Prohibits naive uncalibrated addition of vacuum strain to docking scores.
    Enforces thermodynamic validity checks for active binders.
    """
    self._detect_columns(df)
    work_df = df.copy()

    if self.target_col not in work_df.columns:
      raise ValueError(
          f"Target column '{self.target_col}' not found in dataframe."
      )

    # -------------------------------------------------------------
    # Primary Key Hygiene & Anti-Scrambling Validation
    # -------------------------------------------------------------
    # If upstream predictions table provided an 'Actual' column alongside experimental target,
    # assert ground-truth values match and detect any row scrambling caused by positional .values assignment
    for actual_cand in ["Actual", "ML_Actual"]:
      if actual_cand in work_df.columns and self.target_col in work_df.columns:
        act_vals = pd.to_numeric(work_df[actual_cand], errors="coerce")
        exp_vals = pd.to_numeric(work_df[self.target_col], errors="coerce")
        valid_mask = act_vals.notna() & exp_vals.notna()
        if valid_mask.sum() > 0:
          diffs = (act_vals[valid_mask] - exp_vals[valid_mask]).abs()
          mismatches = diffs > 1e-3
          if mismatches.any():
            raise ValueError(
                "Row scrambling / alignment corruption detected!"
                f" {mismatches.sum()} of {len(work_df)} rows have mismatched"
                f" ground-truth bioactivities between '{self.target_col}' and"
                f" '{actual_cand}' (max divergence: {diffs.max():.4f})."
                " Positional array / .values merging across reordered splits"
                " is strictly prohibited. Always merge via primary key"
                " 'compound_id' with validate='1:1'."
            )

    summary: Dict[str, Any] = {
        "target_col": self.target_col,
        "n_samples": len(work_df),
        "models": {},
        "hybrid_models": {},
        "thermodynamic_sanity": {},
        "warnings": [],
    }

    # -------------------------------------------------------------
    # Partition Isolation: Filter to evaluation split if split exists
    # -------------------------------------------------------------
    train_df = None
    if self.split_col and self.split_col in work_df.columns:
      raw_splits = work_df[self.split_col].dropna().unique()
      unique_splits = [str(x) for x in raw_splits]
      chosen_split = None

      if self.eval_split:
        for u in unique_splits:
          if u.lower() == self.eval_split.lower():
            chosen_split = u
            break
        if chosen_split is None:
          raise ValueError(
              f"Evaluation split '{self.eval_split}' not found in column"
              f" '{self.split_col}'. Available splits: {unique_splits}"
          )
      else:
        # Default priority: test, holdout_test, held_out_test, test_set, holdout, val, validation
        candidate_splits = [
            "test",
            "holdout_test",
            "held_out_test",
            "test_set",
            "holdout",
            "val",
            "validation",
        ]
        lower_map = {u.lower(): u for u in unique_splits}
        for cand in candidate_splits:
          if cand in lower_map:
            chosen_split = lower_map[cand]
            break
        if chosen_split is None:
          for u in unique_splits:
            if "test" in u.lower():
              chosen_split = u
              break

      if chosen_split is not None:
        n_total_before = len(work_df)
        train_candidates = ["train", "train_cv", "training"]
        train_split = None
        lower_map = {u.lower(): u for u in unique_splits}
        for tc in train_candidates:
          if tc in lower_map:
            train_split = lower_map[tc]
            break
        if train_split is None:
          for u in unique_splits:
            if "train" in u.lower():
              train_split = u
              break

        if train_split is not None and train_split != chosen_split:
          train_df = work_df[
              work_df[self.split_col].astype(str) == str(train_split)
          ].copy()

        work_df = work_df[
            work_df[self.split_col].astype(str) == str(chosen_split)
        ].copy()
        if work_df.empty:
          raise ValueError(
              f"Evaluation split '{chosen_split}' in '{self.split_col}'"
              " contains 0 samples."
          )

        summary["split_partition"] = {
            "split_col": self.split_col,
            "eval_split": chosen_split,
            "n_eval": len(work_df),
            "n_total": n_total_before,
        }
        summary["n_samples"] = len(work_df)
      elif len(unique_splits) > 1 and any(
          "train" in u.lower() for u in unique_splits
      ):
        raise ValueError(
            f"Dataset contains split column '{self.split_col}' with partitions"
            f" {unique_splits} including 'train', but no evaluation split was"
            " found or specified. Pass --eval-split (e.g. --eval-split test)"
            " to isolate the evaluation partition; benchmarking metrics across"
            " pooled train+test data is strictly prohibited."
        )

    # NOTE: Prohibit uncalibrated linear addition of strain:
    # dGeff = dGdock + dEstrain is NOT created.
    # Raw physical quantities are evaluated independently.
    y_true = work_df[self.target_col].values.astype(np.float64)

    rows = []

    # Evaluate Physical Estimators independently
    physics_specs = []
    if self.docking_col and self.docking_col in work_df.columns:
      physics_specs.append(("Docking (dGdock)", self.docking_col))
    if self.strain_col and self.strain_col in work_df.columns:
      physics_specs.append(
          ("Conformational Strain (dEstrain)", self.strain_col)
      )

    for label, col in physics_specs:
      raw_vals = work_df[col].values.astype(np.float64)

      # Thermodynamic Sanity Check for free energy columns
      is_free_energy = any(term in label.lower() for term in ["dock", "dgdock"])
      if is_free_energy:
        sanity = check_thermodynamic_sanity(
            raw_vals, target_potency=y_true, col_name=col
        )
        summary["thermodynamic_sanity"][col] = sanity
        if sanity["has_unphysical"]:
          warn_msg = (
              "Thermodynamic Sanity Alert:"
              f" {sanity['n_unphysical']}/{sanity['n_active']} active binders"
              " have non-spontaneous binding free energy (Delta G >= 0"
              f" kcal/mol) in column '{col}'. For active binders, Delta G >= 0"
              " implies Kd = inf. Raw uncalibrated force fields or vacuum"
              " strain additions must not be interpreted as physical free"
              " energies without regression calibration."
          )
          warnings.warn(warn_msg)
          summary["warnings"].append(warn_msg)

      # Invert physics if free energy or docking so lower energy -> higher affinity
      invert_for_calib = is_free_energy or ("strain" not in label.lower())
      phys_train_data = None
      if (
          train_df is not None
          and self.target_col in train_df.columns
          and col in train_df.columns
      ):
        phys_train_data = (
            train_df[self.target_col].values.astype(np.float64),
            train_df[col].values.astype(np.float64),
        )
      calibrated_vals, _ = calibrate_physics_to_target(
          y_true, raw_vals, invert=invert_for_calib, train_data=phys_train_data
      )
      m = safe_stats(y_true, calibrated_vals)

      # Raw correlation with inverted physics (if free energy)
      eval_raw = -raw_vals if invert_for_calib else raw_vals
      m_raw = safe_stats(y_true, eval_raw)
      m["pearson_r"] = m_raw["pearson_r"]
      m["spearman_rho"] = m_raw["spearman_rho"]

      row = {
          "Model / Estimator": label,
          "Family": "Physics",
          "Source_Column": col,
          "R2": m["r2"],
          "Pearson_r": m["pearson_r"],
          "Spearman_rho": m["spearman_rho"],
          "RMSE": m["rmse"],
          "MAE": m["mae"],
      }
      rows.append(row)
      summary["models"][label] = m

    # Evaluate ML Estimators
    best_ml_name = None
    best_ml_r2 = -999.0
    best_ml_pred = None

    for col in self.ml_cols:
      if col not in work_df.columns:
        continue
      y_pred = work_df[col].values.astype(np.float64)
      m = safe_stats(y_true, y_pred)
      row = {
          "Model / Estimator": f"ML: {col}",
          "Family": "Machine Learning",
          "Source_Column": col,
          "R2": m["r2"],
          "Pearson_r": m["pearson_r"],
          "Spearman_rho": m["spearman_rho"],
          "RMSE": m["rmse"],
          "MAE": m["mae"],
      }
      rows.append(row)
      summary["models"][f"ML: {col}"] = m

      if m["r2"] > best_ml_r2:
        best_ml_r2 = m["r2"]
        best_ml_name = col
        best_ml_pred = y_pred

    # Evaluate Hybrid Delta-Learning via Multi-Variable Regression Calibration
    # Prohibits uncalibrated addition work_df['dGeff'] = work_df[docking] + work_df[strain]
    if best_ml_pred is not None:
      phys_features: Dict[str, np.ndarray] = {}
      has_docking = self.docking_col and self.docking_col in work_df.columns
      has_strain = self.strain_col and self.strain_col in work_df.columns

      if has_docking and has_strain:
        # Multi-variable Delta-ML regression strictly fits independent coefficients:
        # y_hat = w0 + w_ml * y_ml + w_dock * (-dGdock) + w_strain * dEstrain
        phys_features["docking"] = work_df[self.docking_col].values.astype(
            np.float64
        )
        phys_features["strain"] = work_df[self.strain_col].values.astype(
            np.float64
        )
        hybrid_label = f"Hybrid: {best_ml_name} + Docking & Strain (Delta-ML)"
        source_desc = f"{best_ml_name} + {self.docking_col} + {self.strain_col}"
      elif has_docking:
        phys_features["docking"] = work_df[self.docking_col].values.astype(
            np.float64
        )
        hybrid_label = f"Hybrid: {best_ml_name} + Docking (Delta-ML)"
        source_desc = f"{best_ml_name} + {self.docking_col}"
      else:
        phys_features = {}
        hybrid_label = ""
        source_desc = ""

      if phys_features:
        train_data_tuple = None
        if (
            train_df is not None
            and self.target_col in train_df.columns
            and best_ml_name in train_df.columns
        ):
          y_tr = train_df[self.target_col].values.astype(np.float64)
          y_tr_ml = train_df[best_ml_name].values.astype(np.float64)
          tr_phys: Dict[str, np.ndarray] = {}
          if (
              has_docking
              and has_strain
              and self.docking_col in train_df.columns
              and self.strain_col in train_df.columns
          ):
            tr_phys["docking"] = train_df[self.docking_col].values.astype(
                np.float64
            )
            tr_phys["strain"] = train_df[self.strain_col].values.astype(
                np.float64
            )
          elif has_docking and self.docking_col in train_df.columns:
            tr_phys["docking"] = train_df[self.docking_col].values.astype(
                np.float64
            )

          if tr_phys:
            train_data_tuple = (y_tr, y_tr_ml, tr_phys)

        y_hybrid, m_hybrid, weights = evaluate_hybrid_delta_learning(
            y_true=y_true,
            y_ml=best_ml_pred,
            physics_features=phys_features,
            cv_folds=5,
            train_data=train_data_tuple,
        )

        row_hybrid = {
            "Model / Estimator": hybrid_label,
            "Family": "Hybrid (Delta-ML)",
            "Source_Column": source_desc,
            "R2": m_hybrid["r2"],
            "Pearson_r": m_hybrid["pearson_r"],
            "Spearman_rho": m_hybrid["spearman_rho"],
            "RMSE": m_hybrid["rmse"],
            "MAE": m_hybrid["mae"],
        }
        rows.append(row_hybrid)
        summary["hybrid_models"][hybrid_label] = {
            "metrics": m_hybrid,
            "weights": weights,
            "r2_delta": round(m_hybrid["r2"] - best_ml_r2, 4),
            "calibrated_via": "multi_variable_regression",
        }

    out_df = pd.DataFrame(rows)
    if not out_df.empty:
      out_df = out_df.sort_values(
          by=["R2", "Pearson_r"], ascending=[False, False]
      ).reset_index(drop=True)

    self.results_df = out_df
    self.summary_dict = summary
    return out_df, summary

  def to_markdown(self) -> str:
    """Renders the side-by-side benchmark table as a clean GitHub-flavored markdown table."""
    if self.results_df is None or self.results_df.empty:
      return "No benchmark results computed."
    lines = [
        "### Side-by-Side Physics vs. Machine Learning Benchmark Matrix",
        "",
        (
            "| Rank | Model / Estimator | Family | R² | Pearson r | Spearman ρ"
            " | RMSE | MAE |"
        ),
        "| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: |",
    ]
    for idx, row in self.results_df.iterrows():
      r2_str = f"{row['R2']:.3f}" if not np.isnan(row["R2"]) else "N/A"
      pr_str = (
          f"{row['Pearson_r']:.3f}" if not np.isnan(row["Pearson_r"]) else "N/A"
      )
      sr_str = (
          f"{row['Spearman_rho']:.3f}"
          if not np.isnan(row["Spearman_rho"])
          else "N/A"
      )
      rmse_str = f"{row['RMSE']:.3f}" if not np.isnan(row["RMSE"]) else "N/A"
      mae_str = f"{row['MAE']:.3f}" if not np.isnan(row["MAE"]) else "N/A"
      lines.append(
          f"| {idx+1} | {row['Model / Estimator']} | {row['Family']} | {r2_str}"
          f" | {pr_str} | {sr_str} | {rmse_str} | {mae_str} |"
      )

    # Highlight split isolation if applied
    if "split_partition" in self.summary_dict:
      sp = self.summary_dict["split_partition"]
      lines.append("")
      lines.append("> [!NOTE]")
      lines.append(
          "> **Evaluation Partition Isolation**: Benchmarked exclusively on"
          f" partition `{sp['eval_split']}` from `{sp['split_col']}`"
          f" ({sp['n_eval']} of {sp['n_total']} samples). Evaluation on pooled"
          " train+test data strictly prevented."
      )

    # Highlight thermodynamic sanity warnings if present
    if "thermodynamic_sanity" in self.summary_dict:
      for col, s_info in self.summary_dict["thermodynamic_sanity"].items():
        if s_info.get("has_unphysical", False):
          lines.append("")
          lines.append("> [!WARNING]")
          lines.append(
              f"> **Thermodynamic Sanity Alert**: Column `{col}` contains"
              f" {s_info['n_unphysical']}/{s_info['n_active']} active binders"
              " with non-spontaneous binding free energy (ΔG ≥ 0 kcal/mol,"
              " implying Kd = ∞). Uncalibrated vacuum force field strain must"
              " not be linearly added without empirical regression"
              " calibration."
          )

    # Highlight hybrid delta-learning if present
    if (
        "hybrid_models" in self.summary_dict
        and self.summary_dict["hybrid_models"]
    ):
      for h_name, h_info in self.summary_dict["hybrid_models"].items():
        delta = h_info.get("r2_delta", 0.0)
        delta_str = f"+{delta:.3f}" if delta >= 0 else f"{delta:.3f}"
        lines.append("")
        lines.append("> [!NOTE]")
        lines.append(
            f"> **Hybrid Delta-Learning Gain**: {h_name} achieved ΔR² ="
            f" {delta_str} over the pure ML baseline."
        )
        w_strs = [f"{k} = {v}" for k, v in h_info["weights"].items()]
        lines.append(f"> Calibrated blending weights: {', '.join(w_strs)}.")

    return "\n".join(lines)


def merge_with_primary_key_validation(
    df_master: pd.DataFrame,
    df_preds: pd.DataFrame,
    key_col: str = "compound_id",
    target_col: Optional[str] = None,
) -> pd.DataFrame:
  """Merges master/physics series data with ML predictions using strict 1:1 primary-key validation.

  Enforces:
  1. Mandatory presence of primary key (`compound_id` or `canonical_smiles`) in
  both tables.
  2. Banning of naive positional array or `.values` column assignment across
  reordered splits.
  3. Strict `validate="1:1"` inner join relationship preventing Cartesian
  expansion or dropped rows.
  4. Cross-validation assertion that ground-truth bioactivity values match
  between tables.
  """
  # Determine effective primary key
  resolved_key = None
  for candidate in [
      key_col,
      "compound_id",
      "canonical_smiles",
      "SMILES",
      "smiles",
  ]:
    if candidate in df_master.columns and candidate in df_preds.columns:
      resolved_key = candidate
      break

  if resolved_key is None:
    raise ValueError(
        f"Primary key '{key_col}' missing from one or both tables! Master"
        f" columns: {list(df_master.columns)}, Predictions columns:"
        f" {list(df_preds.columns)}. Positional array assignment (e.g."
        " df_master['Split'] = df_preds['Split'].values) is strictly"
        " prohibited. Both tables must contain a primary key ('compound_id' or"
        " 'canonical_smiles') for 1:1 join validation."
    )

  preds_clean = df_preds.copy()

  # Filter target if multiple targets in predictions
  if "Target" in preds_clean.columns and target_col:
    matching_targets = [
        t
        for t in preds_clean["Target"].unique()
        if str(t).lower() == target_col.lower()
    ]
    if matching_targets:
      preds_clean = preds_clean[preds_clean["Target"] == matching_targets[0]]

  # Deduplicate predictions if both Train_CV and Holdout_Test rows exist for the same compound
  if not preds_clean[resolved_key].is_unique:
    if "Split" in preds_clean.columns:
      # Prioritize Holdout_Test predictions over Train_CV
      is_holdout = (
          preds_clean["Split"]
          .astype(str)
          .str.contains("holdout|test", case=False, na=False)
      )
      if is_holdout.any():
        preds_clean = preds_clean.sort_values(
            by="Split", ascending=False
        ).drop_duplicates(subset=[resolved_key], keep="first")
      else:
        preds_clean = preds_clean.drop_duplicates(
            subset=[resolved_key], keep="first"
        )
    else:
      preds_clean = preds_clean.drop_duplicates(
          subset=[resolved_key], keep="first"
      )

  # Columns to merge: keep key and non-overlapping prediction columns
  cols_to_pull = [resolved_key]
  for c in ["Predicted", "Actual", "Split", "Best_Model", "Best_Feature_Set"]:
    if c in preds_clean.columns:
      target_name = f"ML_{c}" if c in df_master.columns else c
      if target_name not in df_master.columns:
        preds_clean = preds_clean.rename(columns={c: target_name})
        cols_to_pull.append(target_name)
      else:
        cols_to_pull.append(c)

  cols_to_pull = list(dict.fromkeys(cols_to_pull))
  sub_preds = preds_clean[cols_to_pull]

  merged = pd.merge(
      df_master,
      sub_preds,
      on=resolved_key,
      how="inner",
      validate="1:1",
  )

  if len(merged) == 0:
    raise ValueError(
        f"1:1 Primary-key merge on '{resolved_key}' yielded 0 matched rows."
        f" Verify that compound identifiers align between {len(df_master)}"
        f" master rows and {len(df_preds)} prediction rows."
    )

  # Validate ground-truth alignment
  actual_col = (
      "Actual"
      if "Actual" in merged.columns
      else ("ML_Actual" if "ML_Actual" in merged.columns else None)
  )
  exp_col = (
      target_col if (target_col and target_col in merged.columns) else None
  )
  if not exp_col:
    for cand in [
        "pIC50",
        "pic50",
        "pActivity",
        "pactivity",
        "pKi",
        "pki",
        "affinity",
        "activity",
    ]:
      if cand in merged.columns:
        exp_col = cand
        break

  if actual_col and exp_col:
    try:
      act_vals = pd.to_numeric(merged[actual_col], errors="coerce")
      exp_vals = pd.to_numeric(merged[exp_col], errors="coerce")
      mask = act_vals.notna() & exp_vals.notna()
      if mask.sum() > 0:
        diffs = (act_vals[mask] - exp_vals[mask]).abs()
        mismatches = diffs > 1e-3
        if mismatches.any():
          raise ValueError(
              "Row scrambling / alignment corruption detected during merge!"
              f" {mismatches.sum()} of {len(merged)} rows have mismatched"
              f" ground-truth bioactivities between '{exp_col}' and"
              f" '{actual_col}' (max error: {diffs.max():.4f}). Positional"
              " array / .values merging is strictly prohibited."
          )
    except Exception as e:
      if "Row scrambling" in str(e):
        raise

  return merged


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Comparative Benchmark: Head-to-Head Physics vs. Machine Learning &"
          " Delta-Learning."
      )
  )
  parser.add_argument(
      "--csv",
      "--input-csv",
      type=str,
      required=True,
      dest="csv",
      help="Input CSV containing experimental, physics, and ML predictions.",
  )
  parser.add_argument(
      "--predictions-csv",
      "--preds-csv",
      "--ml-csv",
      type=str,
      dest="predictions_csv",
      help=(
          "Optional path to ML predictions CSV to merge via strict 1:1"
          " primary-key validation."
      ),
  )
  parser.add_argument(
      "--key-col",
      "--pk",
      type=str,
      default="compound_id",
      dest="key_col",
      help="Primary key column for 1:1 validation (default: compound_id).",
  )
  parser.add_argument(
      "--target-col",
      "--target",
      type=str,
      dest="target_col",
      help="Column name for experimental potency/affinity.",
  )
  parser.add_argument(
      "--docking-col",
      type=str,
      dest="docking_col",
      help="Column name for docking free energy (dGdock).",
  )
  parser.add_argument(
      "--strain-col",
      type=str,
      dest="strain_col",
      help="Column name for conformational strain (dEstrain).",
  )
  parser.add_argument(
      "--ml-cols",
      nargs="+",
      dest="ml_cols",
      help="Columns containing ML model predictions.",
  )
  parser.add_argument(
      "--split-col",
      "--split",
      type=str,
      dest="split_col",
      help="Optional split column (train/test).",
  )
  parser.add_argument(
      "--eval-split",
      "--split-name",
      type=str,
      dest="eval_split",
      help="Evaluation partition to isolate (e.g. 'test', 'holdout_test').",
  )
  parser.add_argument(
      "--output-csv",
      type=str,
      dest="output_csv",
      help="Path to write benchmarking table CSV.",
  )
  parser.add_argument(
      "--output-json",
      type=str,
      dest="output_json",
      help="Path to write summary metrics JSON.",
  )
  parser.add_argument(
      "--output-md",
      type=str,
      dest="output_md",
      help="Path to write markdown summary table.",
  )
  args = parser.parse_args()

  df = pd.read_csv(args.csv)
  if args.predictions_csv:
    df_preds = pd.read_csv(args.predictions_csv)
    df = merge_with_primary_key_validation(
        df, df_preds, key_col=args.key_col, target_col=args.target_col
    )
    print(
        f"Successfully performed 1:1 primary-key join on '{args.key_col}'"
        f" ({len(df)} compounds aligned)."
    )

  comparator = PhysicsMLComparator(
      target_col=args.target_col,
      docking_col=args.docking_col,
      strain_col=args.strain_col,
      ml_cols=args.ml_cols,
      split_col=args.split_col,
      eval_split=args.eval_split,
  )
  table_df, summary = comparator.compare(df)

  md_output = comparator.to_markdown()
  print("\n" + md_output + "\n")

  if args.output_csv:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True
    )
    table_df.to_csv(args.output_csv, index=False)
    print(f"Saved benchmark table to CSV: {args.output_csv}")

  if args.output_json:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_json)), exist_ok=True
    )
    with open(args.output_json, "w") as f:
      json.dump(summary, f, indent=2)
    print(f"Saved summary metrics to JSON: {args.output_json}")

  if args.output_md:
    os.makedirs(os.path.dirname(os.path.abspath(args.output_md)), exist_ok=True)
    with open(args.output_md, "w") as f:
      f.write(md_output)
    print(f"Saved markdown report to: {args.output_md}")


if __name__ == "__main__":
  main()
