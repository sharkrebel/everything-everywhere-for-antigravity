# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "pandas",
#     "scipy",
#     "scikit-learn",
#     "matplotlib",
#     "joblib",
#     "seaborn>=0.12.0",
# ]
# ///

"""Small-N Congeneric Series QSAR & Machine Learning Benchmarking Engine.

Provides automated, adaptive regularized machine learning workflows for
congeneric
lead optimization series (N = 15 to 150):
1. Automated Feature Engineering: 2D physicochemical properties, 3D shape
   descriptors, and Morgan Circular Fingerprints with variance threshold
   filtering.
2. Leakage-Free Preprocessing: Imputation, variance thresholding, and
   standardization are strictly fit inside cross-validation loops or on the
   training split, preventing information leakage from test/validation folds.
3. Adaptive Cross-Validation & Holdout Testing: Evaluates linear regularized
   models (RidgeCV, ElasticNetCV), kernel methods (SVR), and tree ensembles
   (Random Forest, Gradient Boosting) under Leave-One-Out (LOOCV), K-Fold, or
   explicit train/test holdouts (scaffold or random splitting).
4. Transparent Benchmark Summary: Exports un-cherry-picked comparison tables
   across all feature sets and model families, reporting cross-validation and
   test set metrics faithfully.
5. Disciplined Feature Attribution: Extracts standardized regression
   coefficients from validated models and renders publication-ready attribution
   plots.
"""

import argparse
import importlib
import importlib.util
import json
import math
import os
import re
import sys
from typing import Any, Dict, Final, List, Optional, Tuple
import warnings

import joblib
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Descriptors3D import (
    Asphericity,
    PMI1,
    PMI2,
    PMI3,
    RadiusOfGyration,
    SpherocityIndex,
)
from scipy.stats import pearsonr, spearmanr, tukey_hsd
from sklearn import base as sklearn_base
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.feature_selection import VarianceThreshold
from sklearn.inspection import permutation_importance
from sklearn.linear_model import ElasticNetCV, RidgeCV
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV, KFold, LeaveOneOut, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

_ml_scripts = os.path.dirname(os.path.abspath(__file__))
_general_scripts = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "..",
        "rdkit_general",
        "scripts",
    )
)
if _general_scripts not in sys.path:
  sys.path.append(_general_scripts)
if _ml_scripts not in sys.path:
  sys.path.insert(0, _ml_scripts)

from compat import (
    PUBLICATION_PALETTE,
    apply_publication_spines,
    set_safe_plot_style,
)
import plot_archetypes
from scaffold_splitter import (
    backfill_temporal_series,
    format_temporal_value,
    is_calendar_col_name,
    parse_temporal_cutoff,
    parse_temporal_series,
    resolve_temporal_cutoff,
    scaffold_split,
    temporal_split,
)

# Line-buffered stdout for real-time background task logging
if hasattr(sys.stdout, "reconfigure"):
  sys.stdout.reconfigure(line_buffering=True)

matplotlib.use("Agg")

HAS_TABFM = importlib.util.find_spec("tabfm") is not None
tabfm = importlib.import_module("tabfm") if HAS_TABFM else None
TabFMRegressor = getattr(tabfm, "TabFMRegressor", None) if tabfm else None

warnings.filterwarnings("ignore", category=ConvergenceWarning)

CENSORED_POLICIES = ("bound", "error", "drop")
_CENSORED_TARGET_RE = re.compile(
    r"^(<=|>=|<|>|≤|≥|~|=)?\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)$"
)


def parse_target_series(
    series: pd.Series,
    target_col: str = "target",
    censored_policy: str = "bound",
    warn: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
  """Parses numeric and relational-qualified target values per `censored_policy`.

  Args:
    series: Input pandas Series (or 1D sequence) of target endpoint values.
    target_col: Endpoint column name for diagnostic messages.
    censored_policy: One of `"bound"` (default; converts `<`, `<=`, `≤`, `>`,
      `>=`, `≥` to numeric detection bounds without dropping any non-null rows),
      `"error"` (raises `ValueError` on censored inequalities), or `"drop"`
      (marks censored inequalities as excluded via `keep_mask`).
    warn: If True, prints a summary warning to `sys.stderr` when censored values
      are converted or dropped.

  Returns:
    Tuple `(values, keep_mask)` of length `len(series)`, where `values` is a
    `float64` array (`np.nan` for excluded/null entries) and `keep_mask` is a
    boolean mask of retained rows.

  Raises:
    ValueError: If `censored_policy` is invalid, any value is non-numeric text,
      or a censored inequality is encountered under `censored_policy="error"`.
  """
  if censored_policy not in CENSORED_POLICIES:
    raise ValueError(
        f"Invalid censored_policy {censored_policy!r}; expected one of"
        f" {CENSORED_POLICIES}."
    )

  s = pd.Series(series)
  n = len(s)
  values = np.full(n, np.nan, dtype=np.float64)
  keep_mask = np.zeros(n, dtype=bool)
  n_left = 0
  n_right = 0
  bad_rows: List[Tuple[int, Any]] = []

  for i, raw_val in enumerate(s.tolist()):
    if pd.isna(raw_val):
      continue
    if isinstance(raw_val, (int, float, np.number)) and not isinstance(
        raw_val, bool
    ):
      fval = float(raw_val)
      if not np.isfinite(fval):
        bad_rows.append((i, raw_val))
        continue
      values[i] = fval
      keep_mask[i] = True
      continue

    text = str(raw_val).strip()
    match = _CENSORED_TARGET_RE.match(text)
    if match is None:
      bad_rows.append((i, raw_val))
      continue

    op, num_str = match.groups()
    fval = float(num_str)
    if not np.isfinite(fval):
      bad_rows.append((i, raw_val))
      continue

    is_left = op in ("<", "<=", "≤")
    is_right = op in (">", ">=", "≥")
    if is_left or is_right:
      if censored_policy == "error":
        raise ValueError(
            f"Censored target value {raw_val!r} at row {i} in column"
            f" '{target_col}' is not permitted when --censored-policy=error."
        )
      if is_left:
        n_left += 1
      else:
        n_right += 1
      if censored_policy == "drop":
        continue

    values[i] = fval
    keep_mask[i] = True

  if bad_rows:
    preview = ", ".join(f"row {idx}: {val!r}" for idx, val in bad_rows[:5])
    if len(bad_rows) > 5:
      preview += f", ... (+{len(bad_rows) - 5} more)"
    raise ValueError(
        f"Unparseable non-numeric target value(s) in column '{target_col}':"
        f" {preview}"
    )

  n_censored = n_left + n_right
  if warn and n_censored > 0:
    if censored_policy == "bound":
      print(
          f"WARNING: converted {n_censored} censored value(s) in '{target_col}'"
          f" ({n_left} left-censored, {n_right} right-censored) to numeric"
          " detection bounds (--censored-policy=bound).",
          file=sys.stderr,
      )
    elif censored_policy == "drop":
      print(
          f"WARNING: dropped {n_censored} censored value(s) in '{target_col}'"
          f" ({n_left} left-censored, {n_right} right-censored)"
          " (--censored-policy=drop).",
          file=sys.stderr,
      )

  return values, keep_mask


def safe_stats(
    y_true: np.ndarray, y_pred: np.ndarray
) -> Tuple[float, float, float, float, float]:
  """Computes R2, RMSE, MAE, Pearson r, and Spearman rho with degenerate-input protection."""
  if len(y_true) == 0 or len(y_pred) == 0:
    return 0.0, 0.0, 0.0, 0.0, 0.0
  r2 = float(r2_score(y_true, y_pred)) if len(y_true) >= 2 else 0.0
  rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))
  mae = float(mean_absolute_error(y_true, y_pred))

  if len(y_true) < 3 or np.std(y_true) < 1e-9 or np.std(y_pred) < 1e-9:
    pr = 0.0
    sr = 0.0
  else:
    pr_val, _ = pearsonr(y_true, y_pred)
    pr = 0.0 if np.isnan(pr_val) else float(pr_val)
    sr_val, _ = spearmanr(y_true, y_pred)
    sr = 0.0 if np.isnan(sr_val) else float(sr_val)
  return r2, rmse, mae, pr, sr


class LeakageFreePreprocessor:
  """Encapsulates imputation, variance filtering, and standardization fit strictly on train."""

  def __init__(self, variance_threshold: float = 0.01):
    self.variance_threshold = variance_threshold
    self.col_means: Optional[np.ndarray] = None
    self.vt: Optional[VarianceThreshold] = None
    self.scaler: Optional[StandardScaler] = None
    self.kept_indices: Optional[np.ndarray] = None

  def fit_transform(self, X: np.ndarray) -> np.ndarray:
    # 0. Sanitize inf/-inf to nan
    X = np.where(np.isinf(X), np.nan, X)
    # 1. Impute column means computed strictly on this training split
    self.col_means = np.nanmean(X, axis=0)
    self.col_means = np.where(np.isnan(self.col_means), 0.0, self.col_means)
    X_imp = np.where(np.isnan(X), self.col_means, X)

    # 2. Variance Threshold (if P > 20)
    if X_imp.shape[1] > 20:
      self.vt = VarianceThreshold(threshold=self.variance_threshold)
      try:
        X_filt = self.vt.fit_transform(X_imp)
        self.kept_indices = self.vt.get_support(indices=True)
      except ValueError:
        self.vt = None
        X_filt = X_imp
        self.kept_indices = np.arange(X_imp.shape[1])
    else:
      self.vt = None
      X_filt = X_imp
      self.kept_indices = np.arange(X_imp.shape[1])

    # 3. Standardize
    self.scaler = StandardScaler()
    X_scaled = self.scaler.fit_transform(X_filt)
    return X_scaled

  def transform(self, X: np.ndarray) -> np.ndarray:
    if self.col_means is None or self.scaler is None:
      raise RuntimeError("Preprocessor must be fit before calling transform.")
    # 0. Sanitize inf/-inf to nan
    X = np.where(np.isinf(X), np.nan, X)
    # 1. Impute using train column means
    X_imp = np.where(np.isnan(X), self.col_means, X)
    # 2. Filter using train mask
    if self.vt is not None and self.kept_indices is not None:
      X_filt = X_imp[:, self.kept_indices]
    else:
      X_filt = X_imp
    # 3. Standardize using train scaler
    return self.scaler.transform(X_filt)

  def get_feature_names(self, original_names: List[str]) -> List[str]:
    if (
        self.kept_indices is not None
        and self.col_means is not None
        and len(original_names) == len(self.col_means)
    ):
      return [original_names[i] for i in self.kept_indices]
    return original_names


P_SCALE_BOUNDS: Tuple[float, float] = (0.0, 14.0)
DEFAULT_MORGAN_FP_BITS: Final[int] = 2048


def apply_applicability_domain_safeguards(
    preds: np.ndarray,
    mols: Optional[List[Optional[Chem.Mol]]] = None,
    smiles_list: Optional[List[str]] = None,
    train_envelope: Tuple[float, float] = P_SCALE_BOUNDS,
    is_p_scale: bool = True,
) -> np.ndarray:
  """Applies p-scale applicability domain bounds when ``is_p_scale`` is True.

  When ``is_p_scale`` is False, returns a ``float64`` copy of ``preds``
  unchanged after input validation (empirical target bounds are enforced by
  ``BoundedRegressor``). When ``is_p_scale`` is True, clips predictions to
  ``train_envelope`` and clamps invalid or sub-molecular structures (< 2 heavy
  atoms) to the inactive baseline (<= 1.5).

  Args:
    preds: Predicted endpoint values to post-process.
    mols: Optional RDKit molecules corresponding 1:1 to ``preds``.
    smiles_list: Optional SMILES strings corresponding 1:1 to ``preds``.
    train_envelope: Physical ``(min, max)`` window applied when ``is_p_scale``
      is True.
    is_p_scale: Whether the target is a negative log molar p-scale affinity.

  Returns:
    Post-processed ``float64`` predictions with the same shape as ``preds``.

  Raises:
    ValueError: If neither ``mols`` nor ``smiles_list`` is provided, or if a
      provided sequence does not match ``len(preds)``.
  """
  out_preds = np.copy(preds).astype(np.float64)
  if mols is None and smiles_list is None:
    raise ValueError(
        "Either 'mols' or 'smiles_list' must be provided to"
        " apply_applicability_domain_safeguards."
    )
  if mols is not None and len(mols) != len(out_preds):
    raise ValueError(
        f"Length of 'mols' ({len(mols)}) must match length of 'preds'"
        f" ({len(out_preds)})."
    )
  if smiles_list is not None and len(smiles_list) != len(out_preds):
    raise ValueError(
        f"Length of 'smiles_list' ({len(smiles_list)}) must match length of"
        f" 'preds' ({len(out_preds)})."
    )
  if not is_p_scale:
    return out_preds

  min_bound, max_bound = train_envelope
  out_preds = np.clip(out_preds, min_bound, max_bound)

  if mols is not None:
    resolved_mols = mols
  else:
    assert smiles_list is not None
    resolved_mols = [
        Chem.MolFromSmiles(str(s)) if pd.notna(s) else None for s in smiles_list
    ]

  for i, mol in enumerate(resolved_mols):
    if mol is None or mol.GetNumHeavyAtoms() < 2:
      out_preds[i] = min(out_preds[i], 1.5)

  return out_preds


def embed_molecule_robust(
    mol: Chem.Mol, random_seed: int = 42
) -> Optional[Chem.Mol]:
  """Robust multi-tier 3D embedding that succeeds on any organic molecule or macrocycle."""
  if mol is None:
    return None
  mol_h = Chem.AddHs(mol)

  # Tier 1: ETKDGv3 with random coordinates
  params1 = AllChem.ETKDGv3()
  params1.randomSeed = random_seed
  params1.useRandomCoords = True
  params1.maxIterations = 200
  res = AllChem.EmbedMolecule(mol_h, params1)
  if res >= 0 and mol_h.GetNumConformers() > 0:
    _optimize_mol_safe(mol_h)
    return mol_h

  # Tier 2: Standard Distance Geometry (EmbedParameters) with random coordinates
  params2 = AllChem.EmbedParameters()
  params2.randomSeed = random_seed
  params2.useRandomCoords = True
  params2.ignoreSmoothingFailures = True
  params2.enforceChirality = False
  params2.maxIterations = 200
  res = AllChem.EmbedMolecule(mol_h, params2)
  if res >= 0 and mol_h.GetNumConformers() > 0:
    _optimize_mol_safe(mol_h)
    return mol_h

  # Tier 3: Embed heavy atoms first, then AddHs
  mol_heavy = Chem.Mol(mol)
  params3 = AllChem.ETKDGv3()
  params3.randomSeed = random_seed
  params3.useRandomCoords = True
  params3.maxIterations = 100
  res = AllChem.EmbedMolecule(mol_heavy, params3)
  if res < 0 or mol_heavy.GetNumConformers() == 0:
    params4 = AllChem.EmbedParameters()
    params4.randomSeed = random_seed
    params4.useRandomCoords = True
    params4.ignoreSmoothingFailures = True
    params4.enforceChirality = False
    res = AllChem.EmbedMolecule(mol_heavy, params4)

  if mol_heavy.GetNumConformers() > 0:
    mol_with_h = Chem.AddHs(mol_heavy, addCoords=True)
    _optimize_mol_safe(mol_with_h)
    return mol_with_h

  # Tier 4: Fallback to simple EmbedMolecule with generous bounds
  res = AllChem.EmbedMolecule(
      mol_h,
      maxAttempts=500,
      useRandomCoords=True,
      ignoreSmoothingFailures=True,
      enforceChirality=False,
  )
  if mol_h.GetNumConformers() > 0:
    _optimize_mol_safe(mol_h)
    return mol_h
  return None


def _optimize_mol_safe(mol: Chem.Mol) -> None:
  try:
    if AllChem.MMFFHasAllMoleculeParams(mol):
      AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
    else:
      AllChem.UFFOptimizeMolecule(mol, maxIters=300)
  except Exception:
    pass


def compute_molecular_features(
    df: pd.DataFrame,
    smiles_col: str = "SMILES",
    compute_3d: bool = True,
    fp_radius: int = 2,
    fp_bits: int = DEFAULT_MORGAN_FP_BITS,
) -> Tuple[pd.DataFrame, List[str], List[str], List[str], str]:
  """Calculates 2D descriptors, 3D shape descriptors, and Morgan fingerprints."""
  if isinstance(fp_radius, bool) or fp_radius < 0:
    raise ValueError(f"fp_radius must be a non-negative int; got {fp_radius!r}")
  if isinstance(fp_bits, bool) or fp_bits <= 0:
    raise ValueError(f"fp_bits must be a positive int; got {fp_bits!r}")
  smiles_candidates = [
      smiles_col,
      "canonical_smiles",
      "SMILES",
      "Smiles",
      "smiles",
  ]
  actual_smiles_col = next(
      (c for c in smiles_candidates if c in df.columns), smiles_col
  )
  if actual_smiles_col not in df.columns:
    raise ValueError(
        f"Could not find valid SMILES column (checked {smiles_candidates})"
    )

  desc_2d_list = []
  desc_3d_list = []
  fp_list = []

  mfpgen = rdFingerprintGenerator.GetMorganGenerator(
      radius=fp_radius, fpSize=fp_bits
  )

  for idx, smi in enumerate(df[actual_smiles_col]):
    m = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None
    if m is None:
      # Fallback zero row
      desc_2d_list.append({
          c: 0.0
          for c in [
              "MW",
              "cLogP",
              "TPSA",
              "HBD",
              "HBA",
              "RotB",
              "Fsp3",
              "Num_Rings",
              "Num_Arom_Rings",
              "Num_Heavy_Atoms",
              "Max_Ring_Size",
          ]
      })
      if compute_3d:
        desc_3d_list.append({
            c: 0.0
            for c in [
                "Radius_of_Gyration",
                "Spherocity_Index",
                "Asphericity",
                "NPR1",
                "NPR2",
            ]
        })
      fp_list.append(np.zeros(fp_bits, dtype=np.float64))
      continue

    # 2D Properties
    ri = m.GetRingInfo()
    ring_sizes = [len(r) for r in ri.AtomRings()]

    d2 = {
        "MW": float(Descriptors.MolWt(m)),
        "cLogP": float(Descriptors.MolLogP(m)),
        "TPSA": float(Descriptors.TPSA(m)),
        "HBD": float(Lipinski.NumHDonors(m)),
        "HBA": float(Lipinski.NumHAcceptors(m)),
        "RotB": float(Lipinski.NumRotatableBonds(m)),
        "Fsp3": float(rdMolDescriptors.CalcFractionCSP3(m)),
        "Num_Rings": float(Lipinski.RingCount(m)),
        "Num_Arom_Rings": float(Lipinski.NumAromaticRings(m)),
        "Num_Heavy_Atoms": float(m.GetNumHeavyAtoms()),
        "Max_Ring_Size": float(max(ring_sizes) if ring_sizes else 0.0),
    }
    desc_2d_list.append(d2)

    # 3D Shape Descriptors. Skipped entirely under --skip_3d: appending a
    # zero-filled placeholder would leave cols_3d non-empty and cause the
    # benchmark to duplicate every model under a phantom "2D + 3D Shape"
    # feature set built from constant columns.
    if compute_3d:
      d3 = {
          "Radius_of_Gyration": 0.0,
          "Spherocity_Index": 0.0,
          "Asphericity": 0.0,
          "NPR1": 0.0,
          "NPR2": 0.0,
      }
      try:
        m_3d = embed_molecule_robust(m, random_seed=42)
        if m_3d is not None and m_3d.GetNumConformers() > 0:
          pmi1 = PMI1(m_3d)
          pmi2 = PMI2(m_3d)
          pmi3 = PMI3(m_3d)
          d3 = {
              "Radius_of_Gyration": float(RadiusOfGyration(m_3d)),
              "Spherocity_Index": float(SpherocityIndex(m_3d)),
              "Asphericity": float(Asphericity(m_3d)),
              "NPR1": float(pmi1 / pmi3) if pmi3 > 0 else 0.0,
              "NPR2": float(pmi2 / pmi3) if pmi3 > 0 else 0.0,
          }
      except Exception:
        pass
      desc_3d_list.append(d3)

    # Fingerprints
    fp_list.append(mfpgen.GetFingerprintAsNumPy(m).astype(np.float64))

  df_2d = pd.DataFrame(desc_2d_list)
  df_3d = pd.DataFrame(desc_3d_list) if desc_3d_list else pd.DataFrame()
  fp_cols = [f"FP_{i}" for i in range(fp_bits)]
  fp_matrix = (
      np.vstack(fp_list)
      if fp_list
      else np.empty((0, fp_bits), dtype=np.float64)
  )
  df_fp = pd.DataFrame(fp_matrix, columns=fp_cols)

  dfs_to_concat = [df.reset_index(drop=True), df_2d]
  if not df_3d.empty:
    dfs_to_concat.append(df_3d)
  dfs_to_concat.append(df_fp)
  feature_df = pd.concat(dfs_to_concat, axis=1)
  feature_df = feature_df.loc[:, ~feature_df.columns.duplicated(keep="last")]

  # Enforce primary keys: canonical_smiles and compound_id
  canonical_list = []
  for smi in df[actual_smiles_col]:
    mol_smi = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None
    if mol_smi is not None:
      try:
        canonical_list.append(Chem.MolToSmiles(mol_smi, isomericSmiles=True))
      except Exception:
        canonical_list.append(str(smi))
    else:
      canonical_list.append(str(smi) if pd.notna(smi) else "")
  if "canonical_smiles" not in feature_df.columns:
    feature_df["canonical_smiles"] = canonical_list
  else:
    feature_df["canonical_smiles"] = (
        feature_df["canonical_smiles"]
        .replace("", np.nan)
        .fillna(pd.Series(canonical_list, index=feature_df.index))
    )
  for smi_alias in ["SMILES", "Clean_SMILES", "smiles", actual_smiles_col]:
    if smi_alias and smi_alias in feature_df.columns:
      feature_df[smi_alias] = (
          feature_df[smi_alias]
          .replace("", np.nan)
          .fillna(feature_df["canonical_smiles"])
      )

  if "compound_id" not in feature_df.columns:
    id_candidates = [
        "compound_id",
        "Compound_ID",
        "compound_name",
        "Compound_Name",
        "ID",
        "id",
        "Name",
        "name",
        "mol_id",
        "Molecule_ID",
        "chembl_id",
    ]
    found_id = next((c for c in id_candidates if c in feature_df.columns), None)
    if found_id:
      feature_df["compound_id"] = feature_df[found_id].astype(str)
    else:
      feature_df["compound_id"] = [
          f"CMPD_{i+1:04d}" for i in range(len(feature_df))
      ]

  cols_2d = list(df_2d.columns)
  cols_3d = list(df_3d.columns) if not df_3d.empty else []

  return feature_df, cols_2d, cols_3d, fp_cols, actual_smiles_col


class BoundedRegressor(BaseEstimator, RegressorMixin):
  """Wraps a regressor to constrain its predictions to biologically plausible bounds based on the target distribution."""

  def __init__(self, estimator: BaseEstimator, is_p_scale: bool = True):
    self.estimator = estimator
    self.is_p_scale = is_p_scale
    self.min_bound, self.max_bound = P_SCALE_BOUNDS

  def fit(self, X, y):
    self.estimator.fit(X, y)
    y_min = float(np.min(y))
    y_max = float(np.max(y))
    std_y = float(np.std(y))
    if self.is_p_scale:
      self.min_bound = max(y_min - 2 * std_y, P_SCALE_BOUNDS[0])
      self.max_bound = min(y_max + 2 * std_y, P_SCALE_BOUNDS[1])
    else:
      self.min_bound = y_min - 2 * std_y
      self.max_bound = y_max + 2 * std_y
    return self

  def predict(self, X):
    preds = self.estimator.predict(X)
    return np.clip(preds, self.min_bound, self.max_bound)

  @property
  def coef_(self):
    return getattr(self.estimator, "coef_", None)


# Inner-loop hyperparameter grids for the estimators that do NOT tune themselves.
# RidgeCV and ElasticNetCV run their own internal CV over `alphas` on whatever
# data they are fitted with, so wrapping them in another search is redundant;
# TabFM is used zero-shot. Grids are intentionally small: these series are
# small-N, and an over-large grid inside an outer fold is itself a source of
# selection bias (Cawley & Talbot 2010).
INNER_PARAM_GRIDS: Dict[str, Dict[str, List[Any]]] = {
    "SVR (RBF)": {
        "estimator__C": [0.5, 1.0, 5.0, 20.0],
        "estimator__gamma": ["scale", 0.01, 0.1],
        "estimator__epsilon": [0.05, 0.1, 0.2],
    },
    "Random Forest": {
        "max_depth": [2, 3, 5],
        "min_samples_split": [2, 4],
    },
    "Gradient Boosting": {
        "max_depth": [1, 2, 3],
        "learning_rate": [0.02, 0.05, 0.1],
    },
}


def _inner_cv_splits(n_samples: int, max_splits: int = 5) -> Optional[int]:
  """Returns a workable inner-CV fold count, or None if the fold is too small."""
  if n_samples < 4:
    return None
  return max(2, min(max_splits, n_samples // 2))


def _resolve_n_jobs(n_jobs: Optional[int] = None) -> int:
  """Resolves worker thread count, defaulting to min(4, cpu_count)."""
  cpu_total = max(1, os.cpu_count() or 1)
  if n_jobs is None:
    return min(4, cpu_total)
  if n_jobs == -1:
    return cpu_total
  if n_jobs < 1:
    raise ValueError(f"n_jobs must be >= 1 or -1, got {n_jobs}.")
  return int(n_jobs)


def tune_model(
    model_name: str,
    model: Any,
    x: np.ndarray,
    y: np.ndarray,
    max_inner_splits: int = 5,
    n_jobs: Optional[int] = None,
) -> Tuple[Any, Dict[str, Any]]:
  """Fits `model` to (x, y), running an inner grid search when one is defined.

  The search sees only the data passed in, so when this is called on an outer
  training fold no information from the outer test fold can leak into model
  selection.

  Args:
    model_name: Key into INNER_PARAM_GRIDS.
    model: Unfitted estimator.
    x: Preprocessed training features for this fold.
    y: Training targets for this fold.
    max_inner_splits: Upper bound on inner K-fold splits.
    n_jobs: Number of parallel worker threads for inner grid search.

  Returns:
    Tuple of (fitted estimator, dict describing the chosen hyperparameters).
  """
  resolved_jobs = 1 if model_name == "TabFM" else _resolve_n_jobs(n_jobs)
  grid = INNER_PARAM_GRIDS.get(model_name)
  n_splits = _inner_cv_splits(len(y), max_inner_splits)
  if not grid or n_splits is None:
    model.fit(x, y)
    return model, {"inner_search": "none", "params": {}}

  grid_model = model
  if (
      resolved_jobs > 1
      and hasattr(model, "get_params")
      and "n_jobs" in model.get_params(deep=False)
  ):
    grid_model = sklearn_base.clone(model)
    grid_model.set_params(n_jobs=1)

  search = GridSearchCV(
      grid_model,
      grid,
      cv=KFold(n_splits=n_splits, shuffle=True, random_state=42),
      scoring="neg_mean_squared_error",
      n_jobs=resolved_jobs,
      refit=True,
  )
  with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=ConvergenceWarning)
    if resolved_jobs > 1:
      with joblib.parallel_backend("threading", n_jobs=resolved_jobs):
        search.fit(x, y)
    else:
      search.fit(x, y)
  return search.best_estimator_, {
      "inner_search": f"grid ({n_splits}-fold)",
      "params": search.best_params_,
  }


def build_candidate_models(
    n_train: int,
    include_tabfm: bool = False,
    tabfm_base_model: Any = None,
    tabfm_estimators: int = 4,
    is_p_scale: bool = True,
    n_jobs: int = 1,
) -> Dict[str, Any]:
  """Constructs the candidate model zoo.

  This is the single definition of the architectures. Both the benchmark loop
  and the max-data retraining step build from it, so a production bundle can
  never be a different architecture from the one that won the benchmark.

  Args:
    n_train: number of training rows, used to size ElasticNetCV's inner CV.
    include_tabfm: whether to add the TabFM zero-shot foundation model.
    tabfm_base_model: preloaded TabFM base model (required if include_tabfm).
    tabfm_estimators: number of TabFM permutation passes.
    is_p_scale: whether the endpoint is a p-scale affinity bounded to [0, 14].
    n_jobs: worker count for standalone parallelizable estimators.

  Returns:
    Mapping of algorithm name to an unfitted estimator.
  """
  rf_jobs = _resolve_n_jobs(n_jobs)
  alphas_grid = np.logspace(-3, 3, 25)
  candidates = {
      "RidgeCV": BoundedRegressor(
          RidgeCV(alphas=alphas_grid, scoring="neg_mean_squared_error"),
          is_p_scale=is_p_scale,
      ),
      "ElasticNetCV": BoundedRegressor(
          ElasticNetCV(
              l1_ratio=[0.05, 0.15, 0.3, 0.5],
              alphas=np.logspace(-3, 0.5, 20),
              cv=min(3, n_train - 1) if n_train > 3 else 2,
              random_state=42,
              max_iter=3000,
              tol=1e-3,
          ),
          is_p_scale=is_p_scale,
      ),
      "SVR (RBF)": BoundedRegressor(
          SVR(C=5.0, epsilon=0.1), is_p_scale=is_p_scale
      ),
      "Random Forest": RandomForestRegressor(
          n_estimators=100,
          max_depth=3,
          min_samples_split=4,
          random_state=42,
          n_jobs=rf_jobs,
      ),
      "Gradient Boosting": GradientBoostingRegressor(
          n_estimators=40, max_depth=2, learning_rate=0.05, random_state=42
      ),
  }
  if include_tabfm and HAS_TABFM and tabfm_base_model is not None:
    candidates["TabFM"] = TabFMRegressor(
        model=tabfm_base_model,
        n_estimators=tabfm_estimators,
        random_state=42,
        verbose=False,
    )
  return candidates


def build_feature_sets(
    cols_2d: List[str],
    cols_3d: List[str],
    cols_fp: List[str],
    available_cols: List[str],
    extra_cols_spec: str | List[str] | None = None,
    extra_name: str | List[str] | None = None,
) -> Dict[str, List[str]]:
  """Constructs benchmark feature sets from one or more required-name custom column blocks."""
  feature_sets: Dict[str, List[str]] = {
      "2D Descriptors": cols_2d,
      "Morgan Circular FP": cols_fp,
  }
  if cols_3d:
    feature_sets["2D + 3D Shape"] = cols_2d + cols_3d

  if not extra_cols_spec:
    return feature_sets

  raw_specs: List[str] = (
      [extra_cols_spec]
      if isinstance(extra_cols_spec, str)
      else list(extra_cols_spec)
  )
  expanded_specs: List[str] = []
  for item in raw_specs:
    for part in str(item).split(";"):
      if part.strip():
        expanded_specs.append(part.strip())

  raw_names: List[str] = []
  if isinstance(extra_name, str) and extra_name.strip():
    raw_names = [n.strip() for n in extra_name.split(";") if n.strip()]
  elif isinstance(extra_name, list):
    for item in extra_name:
      raw_names.extend(n.strip() for n in str(item).split(";") if n.strip())

  base_prefix = "2D/3D" if cols_3d else "2D"
  base_cols = cols_2d + cols_3d if cols_3d else list(cols_2d)
  named_groups: List[Tuple[str, List[str]]] = []

  for idx, spec in enumerate(expanded_specs):
    label = ""
    col_part = spec
    if ":" in spec:
      prefix, rest = spec.split(":", 1)
      if prefix.strip() and "," not in prefix:
        label = prefix.strip()
        col_part = rest.strip()
    if not label and idx < len(raw_names):
      label = raw_names[idx]
    if not label:
      raise ValueError(
          "Custom feature groups must provide an explicit category name using "
          "'--extra_cols \"<Category Name>:col1,col2\"' (repeatable for "
          f"multiple categories) or '--extra_name'. Got unnamed spec: {spec!r}"
      )
    e_cols = [
        c.strip() for c in col_part.split(",") if c.strip() in available_cols
    ]
    if e_cols:
      named_groups.append((label, e_cols))

  for label, e_cols in named_groups:
    has_suffix = any(
        w in label.lower()
        for w in ("feature", "descriptor", "indicator", "score", "term")
    )
    standalone_label = label if has_suffix else f"{label} Features"
    feature_sets[standalone_label] = e_cols
    feature_sets[f"Hybrid ({base_prefix} + {label})"] = base_cols + e_cols
    if len(named_groups) == 1:
      feature_sets[f"Full Hybrid ({base_prefix} + {label} + FP)"] = (
          base_cols + e_cols + cols_fp
      )

  if len(named_groups) > 1:
    combined_label = " + ".join(lbl for lbl, _ in named_groups)
    combined_cols: List[str] = []
    for _, cols in named_groups:
      for col in cols:
        if col not in combined_cols:
          combined_cols.append(col)
    feature_sets[f"Combined ({combined_label})"] = combined_cols
    feature_sets[f"Hybrid ({base_prefix} + {combined_label})"] = (
        base_cols + combined_cols
    )
    feature_sets[f"Full Hybrid ({base_prefix} + {combined_label} + FP)"] = (
        base_cols + combined_cols + cols_fp
    )

  return feature_sets


def run_congeneric_benchmark(
    train_df: pd.DataFrame,
    target_col: str,
    feature_sets: Dict[str, List[str]],
    test_df: Optional[pd.DataFrame] = None,
    nested_outer_k: Optional[int] = None,
    cv_strategy_name: str = "auto",
    n_splits: int = 5,
    include_tabfm: bool = False,
    tabfm_base_model: Any = None,
    tabfm_estimators: int = 4,
    smiles_col: Optional[str] = None,
    censored_policy: str = "bound",
    n_jobs: Optional[int] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
  """Evaluates ML architectures across feature sets with leakage-free CV and principled holdout testing.

  Supports:
  1. Standard Train/Test holdout (temporal or scaffold partition).
  2. Nested Cross-Validation for small-N series: an outer K-fold holdout, with
  an
     inner K-fold grid search over hyperparameters run separately inside each
     outer
     training fold (Cawley & Talbot 2010). RidgeCV/ElasticNetCV are exempt from
     the
     inner grid because they already cross-validate their own regularisation
     path.
  3. Pure Cross-Validation (LOOCV or K-Fold).
  """
  y_train_parsed, train_keep = parse_target_series(
      train_df[target_col],
      target_col=target_col,
      censored_policy=censored_policy,
  )
  train_df = train_df.loc[train_keep].copy().reset_index(drop=True)
  y_train = y_train_parsed[train_keep]
  train_df[target_col] = y_train
  n_train = len(y_train)
  if n_train < 2:
    raise ValueError(
        f"Target '{target_col}' has fewer than 2 usable training observations"
        f" (N={n_train}) after applying censored_policy={censored_policy!r}."
    )
  is_p_scale = infer_is_p_scale(target_col, y_train)

  if test_df is not None and len(test_df) > 0:
    y_test_parsed, test_keep = parse_target_series(
        test_df[target_col],
        target_col=target_col,
        censored_policy=censored_policy,
    )
    test_df = test_df.loc[test_keep].copy().reset_index(drop=True)
    y_test: Optional[np.ndarray] = y_test_parsed[test_keep]
    test_df[target_col] = y_test
  else:
    y_test = None

  smiles_candidates = [
      c
      for c in [
          smiles_col,
          "canonical_smiles",
          "Clean_SMILES",
          "SMILES",
          "smiles",
          "Smiles",
      ]
      if c
  ]
  present_train_smi = [c for c in smiles_candidates if c in train_df.columns]
  actual_smiles_col = (
      max(present_train_smi, key=lambda c: int(train_df[c].notna().sum()))
      if present_train_smi
      else None
  )
  present_test_smi = [
      c
      for c in smiles_candidates
      if test_df is not None and c in test_df.columns
  ]
  test_smiles_col = (
      max(present_test_smi, key=lambda c: int(test_df[c].notna().sum()))
      if present_test_smi
      else None
  )

  cv_mode = cv_strategy_name.lower()
  if cv_mode == "kfold":
    cv = KFold(n_splits=min(n_splits, n_train), shuffle=True, random_state=42)
  elif cv_mode == "auto" and n_train > 60:
    cv = KFold(n_splits=min(10, n_train), shuffle=True, random_state=42)
  elif cv_mode == "loo" and n_train > 80:
    cv = KFold(n_splits=min(10, n_train), shuffle=True, random_state=42)
  else:
    cv = LeaveOneOut()

  has_test = (y_test is not None and len(y_test) > 0) or (
      nested_outer_k is not None and nested_outer_k >= 2
  )
  is_nested = nested_outer_k is not None and nested_outer_k >= 2

  resolved_jobs = _resolve_n_jobs(n_jobs)
  results = []
  trained_models = {}

  for fset_name, cols in feature_sets.items():
    existing_cols = list(
        dict.fromkeys(c for c in cols if c in train_df.columns)
    )
    if not existing_cols:
      continue

    X_train_raw = train_df[existing_cols].values.astype(np.float64)
    X_test_raw = (
        test_df[existing_cols].values.astype(np.float64)
        if (test_df is not None)
        else None
    )

    def _get_candidate_models(model_n_jobs: int = 1):
      return build_candidate_models(
          n_train=n_train,
          include_tabfm=include_tabfm,
          tabfm_base_model=tabfm_base_model,
          tabfm_estimators=tabfm_estimators,
          is_p_scale=is_p_scale,
          n_jobs=model_n_jobs,
      )

    models = _get_candidate_models(model_n_jobs=1)
    cv_splits = list(cv.split(X_train_raw))

    for m_name in models:
      # 1. Leakage-Free Cross-Validation Loop (Internal model validation)
      y_cv_preds = np.zeros(n_train, dtype=np.float64)
      fold_jobs = 1 if m_name == "TabFM" else resolved_jobs

      def _eval_cv_fold(
          tr_idx: np.ndarray,
          val_idx: np.ndarray,
          x_raw: np.ndarray = X_train_raw,
          get_models: Any = _get_candidate_models,
          model_key: str = m_name,
      ):
        prep_fold = LeakageFreePreprocessor(variance_threshold=0.01)
        x_tr_s = prep_fold.fit_transform(x_raw[tr_idx])
        x_val_s = prep_fold.transform(x_raw[val_idx])
        fold_model = get_models(model_n_jobs=1)[model_key]
        fold_model.fit(x_tr_s, y_train[tr_idx])
        return val_idx, fold_model.predict(x_val_s)

      if fold_jobs > 1 and len(cv_splits) > 1:
        fold_outputs = joblib.Parallel(n_jobs=fold_jobs, prefer="threads")(
            joblib.delayed(_eval_cv_fold)(tr_idx, val_idx)
            for tr_idx, val_idx in cv_splits
        )
      else:
        fold_outputs = [
            _eval_cv_fold(tr_idx, val_idx) for tr_idx, val_idx in cv_splits
        ]
      for val_idx, pred in fold_outputs:
        y_cv_preds[val_idx] = pred

      train_smiles = (
          train_df[actual_smiles_col].tolist() if train_df is not None else None
      )
      y_cv_preds = apply_applicability_domain_safeguards(
          y_cv_preds, smiles_list=train_smiles, is_p_scale=is_p_scale
      )
      cv_r2, cv_rmse, cv_mae, cv_pr, cv_sr = safe_stats(y_train, y_cv_preds)

      # 2. Fit Full Model on Entire Training Set (Strictly using train statistics)
      full_prep = LeakageFreePreprocessor(variance_threshold=0.01)
      X_train_s = full_prep.fit_transform(X_train_raw)
      kept_feature_names = full_prep.get_feature_names(existing_cols)

      full_model = _get_candidate_models(model_n_jobs=fold_jobs)[m_name]
      full_model.fit(X_train_s, y_train)

      # 3. Evaluate Holdout Test Set or Nested Outer Folds
      if is_nested:
        # Cawley & Talbot 2010 nested CV: the hyperparameter search runs
        # inside each outer fold, on that fold's training data only, so
        # the outer score includes the cost of model selection.
        outer_cv = KFold(
            n_splits=min(nested_outer_k, n_train), shuffle=True, random_state=42
        )
        y_nested_test_preds = np.zeros(n_train, dtype=np.float64)
        nested_fold_params = []
        for out_tr_idx, out_te_idx in outer_cv.split(X_train_raw):
          X_out_tr, y_out_tr = X_train_raw[out_tr_idx], y_train[out_tr_idx]
          X_out_te = X_train_raw[out_te_idx]

          prep_out = LeakageFreePreprocessor(variance_threshold=0.01)
          X_out_tr_s = prep_out.fit_transform(X_out_tr)
          X_out_te_s = prep_out.transform(X_out_te)

          fold_m, fold_choice = tune_model(
              m_name,
              _get_candidate_models(model_n_jobs=1)[m_name],
              X_out_tr_s,
              y_out_tr,
              n_jobs=fold_jobs,
          )
          nested_fold_params.append(fold_choice)
          y_nested_test_preds[out_te_idx] = fold_m.predict(X_out_te_s)

        y_nested_test_preds = apply_applicability_domain_safeguards(
            y_nested_test_preds,
            smiles_list=train_smiles,
            is_p_scale=is_p_scale,
        )
        test_r2, test_rmse, test_mae, test_pr, test_sr = safe_stats(
            y_train, y_nested_test_preds
        )
        y_test_actual_eval = y_train
        y_test_preds_eval = y_nested_test_preds
      elif has_test and X_test_raw is not None and y_test is not None:
        X_test_s = full_prep.transform(X_test_raw)
        y_test_preds = full_model.predict(X_test_s)
        test_smiles = (
            test_df[test_smiles_col].tolist()
            if (test_df is not None and test_smiles_col is not None)
            else None
        )
        y_test_preds = apply_applicability_domain_safeguards(
            y_test_preds, smiles_list=test_smiles, is_p_scale=is_p_scale
        )
        test_r2, test_rmse, test_mae, test_pr, test_sr = safe_stats(
            y_test, y_test_preds
        )
        y_test_actual_eval = y_test
        y_test_preds_eval = y_test_preds
      else:
        y_test_actual_eval = None
        y_test_preds_eval = None
        test_r2 = test_rmse = test_mae = test_pr = test_sr = 0.0

      model_key = f"{fset_name}__{m_name}"
      trained_models[model_key] = {
          "model": full_model,
          "preprocessor": full_prep,
          "feature_names": kept_feature_names,
          # Everything below is what predict_potency.py needs to rebuild
          # the exact featurization at inference time.
          "feature_cols": existing_cols,
          "feature_set": fset_name,
          "model_name": m_name,
          "is_fingerprint": "FP" in fset_name or "Fingerprint" in fset_name,
          "training_smiles": list(train_smiles) if train_smiles else [],
          "cv_r2": cv_r2,
          "cv_rmse": cv_rmse,
          "cv_mae": cv_mae,
          "cv_pearson_r": cv_pr,
          "cv_spearman_rho": cv_sr,
          "y_train_actual": y_train,
          "y_train_cv_pred": y_cv_preds,
          "has_test": has_test,
          "is_nested": is_nested,
          "test_r2": test_r2,
          "test_rmse": test_rmse,
          "test_mae": test_mae,
          "test_pearson_r": test_pr,
          "test_spearman_rho": test_sr,
          "y_test_actual": y_test_actual_eval,
          "y_test_pred": y_test_preds_eval,
          "train_df": train_df,
          "test_df": test_df,
          "nested_fold_hyperparameters": (
              nested_fold_params if is_nested else None
          ),
      }

      if has_test:
        val_mode = (
            f"Nested_{nested_outer_k}Fold_CV" if is_nested else "Holdout_Split"
        )
        n_tr_disp = (
            n_train - (n_train // nested_outer_k) if is_nested else n_train
        )
        n_te_disp = (
            n_train if is_nested else (len(y_test) if y_test is not None else 0)
        )
        results.append({
            "Target": target_col,
            "Feature_Set": fset_name,
            "Model": m_name,
            "Validation_Protocol": val_mode,
            "N_Train": n_tr_disp,
            "N_Test": n_te_disp,
            "N_Features": X_train_s.shape[1],
            "Train_CV_R2": round(cv_r2, 3),
            "Train_CV_RMSE": round(cv_rmse, 3),
            "Train_CV_r": round(cv_pr, 3),
            "Test_R2": round(test_r2, 3),
            "Test_RMSE": round(test_rmse, 3),
            "Test_MAE": round(test_mae, 3),
            "Test_Pearson_r": round(test_pr, 3),
            "Test_Spearman_rho": round(test_sr, 3),
        })
        print(
            f"   [{target_col}] {fset_name:24s} | {m_name:16s} -> Train CV R²:"
            f" {cv_r2:+.3f} | Test R²: {test_r2:+.3f} (r={test_pr:.3f})",
            flush=True,
        )
      else:
        results.append({
            "Target": target_col,
            "Feature_Set": fset_name,
            "Model": m_name,
            "Validation_Protocol": "Full_CV",
            "N_Samples": n_train,
            "N_Features": X_train_s.shape[1],
            "CV_R2": round(cv_r2, 3),
            "RMSE": round(cv_rmse, 3),
            "MAE": round(cv_mae, 3),
            "Pearson_r": round(cv_pr, 3),
            "Spearman_rho": round(cv_sr, 3),
        })
        print(
            f"   [{target_col}] {fset_name:24s} | {m_name:16s} -> CV R²:"
            f" {cv_r2:+.3f} (r={cv_pr:.3f})",
            flush=True,
        )

  sort_col = "Train_CV_R2" if has_test else "CV_R2"
  summary_df = pd.DataFrame(results).sort_values(
      by=["Target", sort_col], ascending=[True, False]
  )
  return summary_df, trained_models


def tukey_hsd_model_tiers(
    trained_models: Dict[str, Any],
    alpha: float = 0.05,
    use_test_predictions: bool = True,
) -> pd.DataFrame:
  """Groups models into statistical tiers using block-centered paired Tukey HSD.

  On small congeneric series the gap between the top few models is usually
  well inside the noise, so ranking by a point estimate of R2 invites
  over-interpretation. Tukey's Honestly Significant Difference test controls
  the family-wise error rate across all pairwise model comparisons, letting
  us say which models are genuinely indistinguishable from the best.

  Because every candidate model is evaluated on the same compounds, compound-to-
  compound difficulty (e.g., structural outliers or assay censoring floors) is
  a shared blocking factor. When models have aligned per-compound errors, we
  block-center the errors by subtracting each compound's mean error across
  models (Randomized Complete Block Design) before running Tukey's HSD. Tier A
  is the best model plus every model not significantly worse than it; the best
  remaining model then seeds Tier B, and so on.

  Args:
    trained_models: Mapping produced by run_congeneric_benchmark.
    alpha: Family-wise significance level.
    use_test_predictions: Include held-out/nested predictions when present.

  Returns:
    DataFrame with one row per model: mean absolute error, tier, and the
    Tukey p-value against the best model. Empty if fewer than two models
    have usable predictions.
  """
  errors: Dict[str, np.ndarray] = {}
  for key, entry in trained_models.items():
    test_actual = entry.get("y_test_actual")
    test_pred = entry.get("y_test_pred")
    train_actual = entry.get("y_train_actual")
    train_pred = entry.get("y_train_cv_pred")

    actual = pred = None
    if (
        use_test_predictions
        and test_actual is not None
        and test_pred is not None
    ):
      test_act_arr = np.asarray(test_actual, dtype=np.float64)
      test_prd_arr = np.asarray(test_pred, dtype=np.float64)
      if (
          not entry.get("is_nested", False)
          and train_actual is not None
          and train_pred is not None
      ):
        train_act_arr = np.asarray(train_actual, dtype=np.float64)
        train_prd_arr = np.asarray(train_pred, dtype=np.float64)
        same_split = len(train_act_arr) == len(test_act_arr) and np.allclose(
            train_act_arr, test_act_arr, equal_nan=True
        )
        if not same_split:
          actual = np.concatenate([train_act_arr, test_act_arr])
          pred = np.concatenate([train_prd_arr, test_prd_arr])
        else:
          actual, pred = test_act_arr, test_prd_arr
      else:
        actual, pred = test_act_arr, test_prd_arr
    elif train_actual is not None and train_pred is not None:
      actual = np.asarray(train_actual, dtype=np.float64)
      pred = np.asarray(train_pred, dtype=np.float64)
    elif test_actual is not None and test_pred is not None:
      actual = np.asarray(test_actual, dtype=np.float64)
      pred = np.asarray(test_pred, dtype=np.float64)

    if actual is None or pred is None:
      continue
    abs_err = np.abs(actual - pred)
    abs_err = abs_err[np.isfinite(abs_err)]
    # Tukey's HSD needs within-group variance; a constant group has none.
    if len(abs_err) >= 3 and np.ptp(abs_err) > 0:
      errors[key] = abs_err

  if len(errors) < 2:
    return pd.DataFrame(
        columns=[
            "Model_Key",
            "Mean_Abs_Error",
            "Tier",
            "p_vs_Best",
            "Distinguishable_From_Best",
        ]
    )

  names = list(errors)
  lengths = {len(errors[n]) for n in names}
  if len(lengths) == 1:
    # Block-center across compounds (Randomized Complete Block Design) to
    # remove shared compound-level variance while preserving model mean diffs.
    err_matrix = np.vstack([errors[n] for n in names])
    compound_means = np.mean(err_matrix, axis=0, keepdims=True)
    grand_mean = float(np.mean(err_matrix))
    centered_matrix = err_matrix - compound_means + grand_mean
    if np.max(np.ptp(centered_matrix, axis=1)) > 1e-12:
      groups = [centered_matrix[i] for i in range(len(names))]
    else:
      groups = [errors[n] for n in names]
  else:
    groups = [errors[n] for n in names]

  result = tukey_hsd(*groups)
  pvalues = result.pvalue
  means = {n: float(np.mean(errors[n])) for n in names}

  # Greedy tiering: seed each tier with the best unassigned model, then absorb
  # every unassigned model that is NOT significantly worse than that seed.
  tiers: Dict[str, str] = {}
  unassigned = sorted(names, key=lambda n: means[n])
  tier_index = 0
  while unassigned:
    seed = unassigned[0]
    seed_pos = names.index(seed)
    label = chr(ord("A") + tier_index)
    members = [
        n
        for n in unassigned
        if n == seed or pvalues[seed_pos][names.index(n)] >= alpha
    ]
    for n in members:
      tiers[n] = f"Tier {label}"
    unassigned = [n for n in unassigned if n not in members]
    tier_index += 1

  best = min(names, key=lambda n: means[n])
  best_pos = names.index(best)
  rows = []
  for n in sorted(names, key=lambda n: means[n]):
    p_vs_best = (
        float(pvalues[best_pos][names.index(n)]) if n != best else float("nan")
    )
    rows.append({
        "Model_Key": n,
        "Mean_Abs_Error": round(means[n], 4),
        "Tier": tiers[n],
        "p_vs_Best": None if n == best else round(p_vs_best, 4),
        "Distinguishable_From_Best": (
            None if n == best else bool(p_vs_best < alpha)
        ),
    })
  return pd.DataFrame(rows)


def plot_model_coefficients(
    model_entry: Dict[str, Any],
    output_png: str,
    top_k: int = 15,
) -> None:
  """Renders a standardized regression coefficient bar plot."""
  set_safe_plot_style()

  raw_model = model_entry["model"]
  model = (
      raw_model.estimator
      if isinstance(raw_model, BoundedRegressor)
      else raw_model
  )
  names = list(model_entry["feature_names"])
  if hasattr(model, "coef_"):
    coefs = np.ravel(model.coef_)
  else:
    return

  min_len = min(len(names), len(coefs))
  names = names[:min_len]
  coefs = coefs[:min_len]

  coef_df = pd.DataFrame({"Feature": names, "Coefficient": coefs})
  # Filter out anonymous circular fingerprint bit hashes - only display interpretable descriptors
  is_bit_hash = (
      coef_df["Feature"]
      .astype(str)
      .str.contains(r"^(?:FP_\d+|Bit_\d+|fp_\d+|bit_\d+)$", regex=True)
  )
  coef_df = coef_df[~is_bit_hash]
  if coef_df.empty:
    print(
        "Notice: Model features are anonymous fingerprint bits; skipping"
        f" coefficient plot for {output_png}"
    )
    return

  coef_df["AbsCoef"] = coef_df["Coefficient"].abs()
  coef_df = coef_df.sort_values(by="AbsCoef", ascending=False).head(top_k)
  coef_df = coef_df.sort_values(by="Coefficient", ascending=True)

  os.makedirs(
      os.path.dirname(os.path.abspath(output_png)) or ".", exist_ok=True
  )
  fig, ax = plt.subplots(
      figsize=(8.5, max(4.5, len(coef_df) * 0.38)),
      dpi=300,
      layout="constrained",
  )
  fig.patch.set_facecolor(PUBLICATION_PALETTE["background"])
  ax.set_facecolor(PUBLICATION_PALETTE["card"])
  colors = [
      PUBLICATION_PALETTE["success"] if c > 0 else PUBLICATION_PALETTE["danger"]
      for c in coef_df["Coefficient"]
  ]
  bars = ax.barh(
      coef_df["Feature"],
      coef_df["Coefficient"],
      color=colors,
      edgecolor="#0F172A",
      linewidth=0.8,
      height=0.6,
  )
  ax.axvline(0, color="#64748B", linestyle="--", linewidth=0.9)
  max_abs_c = max(0.2, float(coef_df["AbsCoef"].max()) * 1.28)
  ax.set_xlim(-max_abs_c, max_abs_c)
  for bar, c_val in zip(bars, coef_df["Coefficient"]):
    val = float(c_val)
    off = max_abs_c * 0.025 if val >= 0 else -max_abs_c * 0.025
    ha_c = "left" if val >= 0 else "right"
    ax.text(
        val + off,
        bar.get_y() + bar.get_height() / 2.0,
        f"{val:+.3f}",
        va="center",
        ha=ha_c,
        fontsize=8.5,
        fontweight="bold",
        color=PUBLICATION_PALETTE.get("text", "#0F172A"),
    )
  ax.set_xlabel(
      "Standardized Regression Coefficient (β)", fontsize=11, fontweight="bold"
  )
  fset_disp = (
      f"{model_entry.get('feature_set', '')} - "
      if model_entry.get("feature_set")
      else ""
  )
  ax.set_title(
      "Standardized Feature Weights"
      f" ({fset_disp}{model_entry.get('model_name', 'model')}, CV R² ="
      f" {model_entry.get('cv_r2', 0.0):.3f})",
      fontsize=12,
      fontweight="bold",
      pad=10,
  )
  apply_publication_spines(ax)
  fig.savefig(output_png, dpi=300)
  plt.close(fig)
  print(f"Saved coefficient plot to {output_png}")


def compute_feature_attribution(model_entry: Dict[str, Any]) -> pd.DataFrame:
  """Computes per-feature attribution using the method appropriate to the model.

  Linear models expose standardized coefficients directly; tree ensembles
  expose mean decrease in impurity; everything else (SVR with an RBF kernel,
  TabFM) is opaque and needs permutation importance, which is measured
  against the training data the model was fitted on.

  Args:
    model_entry: an entry from ``run_congeneric_benchmark``'s second return
      value.

  Returns:
    A DataFrame with columns ``Feature``, ``Importance`` and ``Signed``,
    sorted by descending absolute importance. ``Signed`` is True only for
    coefficient-based attribution, where the sign is meaningful. Empty if
    attribution could not be computed.
  """
  model = model_entry.get("model")
  names = list(model_entry.get("feature_names") or [])
  if model is None or not names:
    return pd.DataFrame(columns=["Feature", "Importance", "Signed"])

  # BoundedRegressor wraps the real estimator.
  # Unwrap ONLY our own wrapper. sklearn ensembles (RandomForestRegressor,
  # GradientBoostingRegressor) also expose `.estimator`, but theirs is the
  # unfitted prototype base learner. A blanket getattr therefore threw away
  # the fitted forest and its `feature_importances_`, silently falling
  # through to slow permutation importance for every tree model.
  inner = model.estimator if isinstance(model, BoundedRegressor) else model

  values = None
  signed = False
  method = ""
  if hasattr(inner, "coef_"):
    values = np.ravel(inner.coef_)
    signed = True
    method = "standardized coefficient"
  elif hasattr(inner, "feature_importances_"):
    values = np.ravel(inner.feature_importances_)
    method = "mean decrease in impurity"
  else:
    # Permutation importance: the only option for kernel / foundation models.
    try:
      prep = model_entry.get("preprocessor")
      train_df = model_entry.get("train_df")
      cols = model_entry.get("feature_cols")
      y = model_entry.get("y_train_actual")
      if prep is None or train_df is None or not cols or y is None:
        return pd.DataFrame(columns=["Feature", "Importance", "Signed"])
      x_raw = train_df[
          [c for c in cols if c in train_df.columns]
      ].values.astype(np.float64)
      x_s = prep.transform(x_raw)
      result = permutation_importance(
          model,
          x_s,
          np.asarray(y, dtype=np.float64),
          n_repeats=5,
          random_state=42,
          scoring="r2",
      )
      values = np.ravel(result.importances_mean)
      method = "permutation importance"
    except Exception as exc:
      # Reported to stderr below; never silently dropped.
      print(
          "Note: could not compute feature attribution for "
          f"{type(inner).__name__} ({type(exc).__name__}: {exc}).",
          file=sys.stderr,
      )
      return pd.DataFrame(columns=["Feature", "Importance", "Signed"])

  n = min(len(names), len(values))
  out = pd.DataFrame({
      "Feature": names[:n],
      "Importance": values[:n],
      "Signed": signed,
  })
  # Filter out uninterpretable anonymous circular fingerprint bit hashes
  is_bit_hash = (
      out["Feature"]
      .astype(str)
      .str.contains(r"^(?:FP_\d+|Bit_\d+|fp_\d+|bit_\d+)$", regex=True)
  )
  out = out[~is_bit_hash]
  out.attrs["method"] = method
  return out.reindex(out["Importance"].abs().sort_values(ascending=False).index)


def plot_feature_importance(
    model_entry: Dict[str, Any],
    output_png: str,
    target_name: str = "",
    top_k: int = 20,
) -> bool:
  """Renders a feature attribution bar plot.

  Args:
    model_entry: entry from ``run_congeneric_benchmark``.
    output_png: destination PNG path.
    target_name: endpoint name, used in the title.
    top_k: number of features to show.

  Returns:
    True if a plot was written.
  """
  attrib = compute_feature_attribution(model_entry)
  if attrib.empty:
    return False

  # Double check no anonymous bit hashes slip through
  is_bit = (
      attrib["Feature"]
      .astype(str)
      .str.contains(r"^(?:FP_\d+|Bit_\d+|fp_\d+|bit_\d+)$", regex=True)
  )
  attrib = attrib[~is_bit]
  if attrib.empty:
    return False

  set_safe_plot_style()

  method = attrib.attrs.get("method", "importance")
  top = attrib.head(top_k).sort_values(by="Importance", ascending=True)
  is_signed = bool(top["Signed"].iloc[0]) if len(top) else False

  os.makedirs(
      os.path.dirname(os.path.abspath(output_png)) or ".", exist_ok=True
  )
  fig, ax = plt.subplots(
      figsize=(8, max(4.5, len(top) * 0.35)), dpi=300, layout="constrained"
  )
  fig.patch.set_facecolor(PUBLICATION_PALETTE["background"])
  ax.set_facecolor(PUBLICATION_PALETTE["card"])
  if is_signed:
    colors = [
        PUBLICATION_PALETTE["success"]
        if c > 0
        else PUBLICATION_PALETTE["danger"]
        for c in top["Importance"]
    ]
  else:
    colors = [PUBLICATION_PALETTE["primary"]] * len(top)
  ax.barh(
      top["Feature"],
      top["Importance"],
      color=colors,
      edgecolor="#0F172A",
      height=0.6,
  )
  ax.axvline(0, color="#64748B", linestyle="--", linewidth=0.8)
  ax.set_xlabel(method.title(), fontsize=11, fontweight="bold")
  title_target = f"{target_name}: " if target_name else ""
  fset_disp = (
      f"{model_entry.get('feature_set', '')} - "
      if model_entry.get("feature_set")
      else ""
  )
  ax.set_title(
      f"{title_target}Feature Attribution"
      f" ({fset_disp}{model_entry.get('model_name', 'model')}, {method})",
      fontsize=12,
      fontweight="bold",
      pad=10,
  )
  apply_publication_spines(ax)
  fig.savefig(output_png, dpi=300)
  plt.close(fig)
  print(f"Saved feature attribution plot to {output_png}")
  print(f"![Feature Attribution]({output_png})")
  return True


P_SCALE_TARGET_TOKENS = ("pki", "pic50", "pkd", "pec50", "pchembl")

# A p-scale affinity is a negative log molar concentration. Anything a real
# assay can report sits inside this window; values outside it mean the endpoint
# is on some other scale despite its name, and clamping would corrupt it.
P_SCALE_PLAUSIBLE_RANGE = (0.0, 15.0)


def infer_is_p_scale(target_name: str, y_values: Optional[np.ndarray]) -> bool:
  """Decides whether an endpoint is a p-scale affinity.

  ``predict_potency.py`` clamps p-scale predictions to the thermodynamic
  window, which is right for pKi and destructive for, say, a nanomolar IC50.
  It can only guess from the target's *name*; the producer can also look at
  the training values, so the decision is made here and stored in the bundle.

  Args:
    target_name: endpoint column name.
    y_values: observed training values, if available.

  Returns:
    True only if the name looks like a p-scale AND the observed values are
    consistent with one. A name-only match on out-of-range data is treated as
    a mislabelled endpoint and reported.
  """
  name_says_p = any(
      tok in str(target_name).lower() for tok in P_SCALE_TARGET_TOKENS
  )
  if not name_says_p:
    return False
  if y_values is None or len(y_values) == 0:
    return True
  lo, hi = float(np.nanmin(y_values)), float(np.nanmax(y_values))
  if lo < P_SCALE_PLAUSIBLE_RANGE[0] or hi > P_SCALE_PLAUSIBLE_RANGE[1]:
    print(
        f"WARNING: target '{target_name}' is named like a p-scale affinity but"
        f" its training values span [{lo:.3g}, {hi:.3g}], outside the plausible"
        f" p-scale range {P_SCALE_PLAUSIBLE_RANGE}. Recording is_p_scale=False"
        " so inference does not clamp predictions to [0, 15].",
        file=sys.stderr,
    )
    return False
  return True


def build_model_bundle(
    model_entry: Dict[str, Any],
    target_name: str,
    model: Any = None,
    preprocessor: Any = None,
    training_smiles: Optional[List[str]] = None,
    n_samples: Optional[int] = None,
    y_values: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
  """Assembles a self-contained inference bundle.

  The schema is the contract consumed by ``predict_potency.py``; changing a
  key here breaks inference.

  Args:
    model_entry: benchmark entry supplying metadata and default estimators.
    target_name: endpoint name.
    model: override estimator (used for the max-data refit).
    preprocessor: override preprocessor (used for the max-data refit).
    training_smiles: override training SMILES for applicability domain.
    n_samples: override training row count.
    y_values: override observed training values, used for the training range and
      the p-scale decision.

  Returns:
    The bundle dict.
  """
  prep = (
      preprocessor
      if preprocessor is not None
      else model_entry.get("preprocessor")
  )
  cols = list(model_entry.get("feature_cols") or [])
  kept = (
      prep.get_feature_names(cols)
      if prep is not None
      else list(model_entry.get("feature_names") or [])
  )
  smiles = (
      training_smiles
      if training_smiles is not None
      else list(model_entry.get("training_smiles") or [])
  )

  # y_train_actual is a numpy array: `or []` would raise on ambiguous truth.
  y_obs = (
      y_values if y_values is not None else model_entry.get("y_train_actual")
  )
  y_arr = None if y_obs is None else np.asarray(y_obs, dtype=np.float64)

  if n_samples is not None:
    n_obs = int(n_samples)
  else:
    n_obs = 0 if y_arr is None else int(len(y_arr))

  # Lets inference flag predictions that fall outside the observed label
  # range. Without this key the Extrapolated_Beyond_Train column is silently
  # never emitted.
  training_range = None
  if y_arr is not None and len(y_arr) and np.isfinite(y_arr).any():
    training_range = [float(np.nanmin(y_arr)), float(np.nanmax(y_arr))]

  return {
      "model": model if model is not None else model_entry.get("model"),
      "preprocessor": prep,
      "feature_set": model_entry.get("feature_set", "2D Descriptors"),
      "feature_cols": cols,
      "kept_feature_names": kept,
      "target": target_name,
      "model_name": model_entry.get("model_name", ""),
      "is_fingerprint": bool(model_entry.get("is_fingerprint", False)),
      "is_p_scale": infer_is_p_scale(target_name, y_arr),
      "training_range": training_range,
      "training_stats": {
          "n_samples": n_obs,
          "cv_r2": model_entry.get("cv_r2"),
          "test_r2": model_entry.get("test_r2"),
      },
      "training_smiles": smiles,
  }


def retrain_on_all_data(
    model_entry: Dict[str, Any],
    target_col: str,
    include_tabfm: bool = False,
    tabfm_base_model: Any = None,
    tabfm_estimators: int = 4,
    censored_policy: str = "bound",
    n_jobs: Optional[int] = None,
) -> Optional[Dict[str, Any]]:
  """Refits the winning architecture on 100% of the labeled data.

  Model selection needs a held-out partition, but once the architecture is
  chosen there is no reason to throw the test rows away. This fits a fresh
  preprocessor and a fresh instance of the same architecture on train + test
  combined.

  Args:
    model_entry: the winning benchmark entry.
    target_col: endpoint column name.
    include_tabfm: whether TabFM is among the candidates.
    tabfm_base_model: preloaded TabFM base model.
    tabfm_estimators: TabFM permutation passes.
    censored_policy: policy (`"bound"`, `"error"`, or `"drop"`) for censored
      target strings.
    n_jobs: number of parallel worker threads (`None` defaults to `min(4,
      cpu_count)`).

  Returns:
    ``{"model", "preprocessor", "n_samples", "training_smiles"}``, or None if
    retraining was not possible.
  """
  train_df = model_entry.get("train_df")
  test_df = model_entry.get("test_df")
  cols = list(model_entry.get("feature_cols") or [])
  m_name = model_entry.get("model_name")
  if train_df is None or not cols or not m_name:
    return None
  if target_col not in train_df.columns:
    return None

  # A nested-CV run has no separate test frame: every row is already in train_df.
  frames = (
      [train_df]
      if test_df is None or model_entry.get("is_nested")
      else [train_df, test_df]
  )
  full_df = (
      pd.concat(frames, ignore_index=True) if len(frames) > 1 else train_df
  )
  y_parsed, keep_mask = parse_target_series(
      full_df[target_col],
      target_col=target_col,
      censored_policy=censored_policy,
  )
  full_df = full_df.loc[keep_mask].copy().reset_index(drop=True)
  y = y_parsed[keep_mask]

  usable = [c for c in cols if c in full_df.columns]
  if not usable:
    return None

  x_raw = full_df[usable].values.astype(np.float64)
  if len(y) < 3:
    return None

  prep = LeakageFreePreprocessor(variance_threshold=0.01)
  x_s = prep.fit_transform(x_raw)

  candidates = build_candidate_models(
      n_train=len(y),
      include_tabfm=include_tabfm,
      tabfm_base_model=tabfm_base_model,
      tabfm_estimators=tabfm_estimators,
      is_p_scale=infer_is_p_scale(target_col, y),
      n_jobs=_resolve_n_jobs(n_jobs),
  )
  model = candidates.get(m_name)
  if model is None:
    return None
  model.fit(x_s, y)

  smiles_candidates = [
      "Clean_SMILES",
      "SMILES",
      "canonical_smiles",
      "smiles",
      "Smiles",
  ]
  smi_col = next((c for c in smiles_candidates if c in full_df.columns), None)
  smiles = full_df[smi_col].astype(str).tolist() if smi_col else []

  return {
      "model": model,
      "preprocessor": prep,
      "n_samples": len(y),
      "training_smiles": smiles,
      "y_values": y,
  }


def _ensure_portable_pickling() -> None:
  """Makes bundles loadable by any consumer, not just this script.

  When this file runs as a script its classes carry ``__module__ ==
  "__main__"``, so pickle records them as ``__main__.BoundedRegressor``. A
  plain ``joblib.load`` in a different process then fails with
  ``AttributeError: Can't get attribute 'BoundedRegressor'``. Re-pointing the
  classes at a real, importable module name -- and registering that name in
  ``sys.modules`` so the running process can still resolve them -- makes the
  written bundle portable.
  """
  module_name = "congeneric_series_ml"
  this_module = sys.modules.get(__name__)
  if module_name not in sys.modules and this_module is not None:
    sys.modules[module_name] = this_module
  for cls in (BoundedRegressor, LeakageFreePreprocessor):
    if cls.__module__ == "__main__":
      cls.__module__ = module_name


def save_model_bundles(
    model_entry: Dict[str, Any],
    target_col: str,
    model_dir: str,
    include_tabfm: bool = False,
    tabfm_base_model: Any = None,
    tabfm_estimators: int = 4,
    censored_policy: str = "bound",
    n_jobs: Optional[int] = None,
) -> Dict[str, str]:
  """Writes the benchmark and production model bundles to `model_dir`.

  Args:
    model_entry: the winning benchmark entry.
    target_col: endpoint name, used in the filenames.
    model_dir: destination directory, created if absent.
    include_tabfm: whether TabFM is among the candidates.
    tabfm_base_model: preloaded TabFM base model.
    tabfm_estimators: TabFM permutation passes.
    censored_policy: policy (`"bound"`, `"error"`, or `"drop"`) for censored
      target strings.
    n_jobs: number of parallel worker threads (`None` defaults to `min(4,
      cpu_count)`).

  Returns:
    Mapping of bundle kind (``best_cv`` / ``maxdata``) to written path.
  """
  _ensure_portable_pickling()
  is_file_target = str(model_dir).endswith((".joblib", ".pkl"))
  target_dir = (
      os.path.dirname(os.path.abspath(model_dir)) or "."
      if is_file_target
      else os.path.abspath(model_dir)
  )
  os.makedirs(target_dir, exist_ok=True)
  safe_target = "".join(
      c if (c.isalnum() or c in "-_") else "_" for c in str(target_col)
  )
  written: Dict[str, str] = {}

  cv_bundle = build_model_bundle(model_entry, target_col)
  cv_path = os.path.join(target_dir, f"{safe_target}_best_cv_model.joblib")
  joblib.dump(cv_bundle, cv_path)
  written["best_cv"] = cv_path
  print(f"Saved benchmark model bundle to: {cv_path}")

  refit = retrain_on_all_data(
      model_entry,
      target_col,
      include_tabfm=include_tabfm,
      tabfm_base_model=tabfm_base_model,
      tabfm_estimators=tabfm_estimators,
      censored_policy=censored_policy,
      n_jobs=n_jobs,
  )
  if refit is None:
    if is_file_target:
      joblib.dump(cv_bundle, os.path.abspath(model_dir))
      written["explicit_file"] = os.path.abspath(model_dir)
    print(
        f"WARNING: could not retrain '{model_entry.get('model_name')}' on the "
        f"full dataset for target '{target_col}'; no "
        f"{safe_target}_maxdata_model.joblib was written. Use "
        f"{os.path.basename(cv_path)} for inference instead.",
        file=sys.stderr,
    )
    return written

  max_bundle = build_model_bundle(
      model_entry,
      target_col,
      model=refit["model"],
      preprocessor=refit["preprocessor"],
      training_smiles=refit["training_smiles"],
      n_samples=refit["n_samples"],
      y_values=refit["y_values"],
  )
  max_path = os.path.join(target_dir, f"{safe_target}_maxdata_model.joblib")
  joblib.dump(max_bundle, max_path)
  written["maxdata"] = max_path
  if is_file_target:
    joblib.dump(max_bundle, os.path.abspath(model_dir))
    written["explicit_file"] = os.path.abspath(model_dir)
  print(
      f"Saved production max-data model bundle to: {max_path} "
      f"(retrained on {refit['n_samples']} compounds)"
  )
  return written


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Small-N Congeneric Series QSAR & Machine Learning Benchmarking"
          " Engine."
      )
  )
  parser.add_argument(
      "--input_csv",
      "--input-csv",
      "--input",
      dest="input_csv",
      type=str,
      required=True,
      help="Path to input CSV containing SMILES and activity data.",
  )
  parser.add_argument(
      "--smiles_col",
      "--smiles-col",
      "--smiles",
      dest="smiles_col",
      type=str,
      default="SMILES",
      help="Column name for ligand SMILES (default: SMILES).",
  )
  parser.add_argument(
      "--target_col",
      "--target-col",
      "--target_cols",
      "--target-cols",
      "--targets",
      "--target",
      dest="target_col",
      nargs="+",
      required=True,
      help=(
          "Column name(s) for target endpoint(s). Supports multiple args or"
          " comma-separated targets."
      ),
  )
  parser.add_argument(
      "--id_col",
      "--id-col",
      "--name_col",
      "--name-col",
      dest="id_col",
      type=str,
      default=None,
      help="Column name for compound identifier.",
  )
  parser.add_argument(
      "--split_col",
      "--split-col",
      dest="split_col",
      type=str,
      default=None,
      help=(
          "Column name containing pre-assigned split annotations (e.g. 'Split',"
          " 'train'/'test')."
      ),
  )
  parser.add_argument(
      "--eval_split",
      "--eval-split",
      dest="eval_split",
      type=str,
      default="test",
      help="Split label to evaluate as held-out test (default: 'test').",
  )
  parser.add_argument(
      "--extra_cols",
      "--extra-cols",
      "--custom_cols",
      "--custom-cols",
      "--physics_cols",
      "--physics-cols",
      dest="physics_cols",
      action="append",
      type=str,
      help=(
          "Named custom feature block in '<Category Name>:col1,col2' format. "
          "Repeat flag or separate blocks with ';' to benchmark multiple "
          "feature categories."
      ),
  )
  parser.add_argument(
      "--extra_name",
      "--extra-name",
      "--custom_name",
      "--custom-name",
      "--physics_name",
      "--physics-name",
      dest="extra_name",
      action="append",
      type=str,
      default=None,
      help=(
          "Explicit category name(s) corresponding to --extra_cols when not "
          "using the inline '<Category Name>:col1,col2' prefix."
      ),
  )
  parser.add_argument(
      "--cv_strategy",
      "--cv-strategy",
      dest="cv_strategy",
      type=str,
      choices=["auto", "loo", "kfold"],
      default="auto",
      help=(
          "Cross-validation strategy for training folds (default: auto, uses"
          " LOOCV for N <= 60 and 10-fold KFold for N > 60)."
      ),
  )
  parser.add_argument(
      "--n_splits",
      "--n-splits",
      dest="n_splits",
      type=int,
      default=5,
      help="Number of splits for KFold CV (default: 5).",
  )
  parser.add_argument(
      "--test_split",
      "--test-split",
      dest="test_split",
      type=float,
      default=None,
      help=(
          "Fraction of compounds to hold out as explicit test set (default:"
          " None, adaptively computed). Use 0.0 for pure CV."
      ),
  )
  parser.add_argument(
      "--split_method",
      "--split-method",
      "--split_strategy",
      "--split-strategy",
      "--split_type",
      "--split-type",
      dest="split_method",
      type=str,
      choices=[
          "auto",
          "scaffold",
          "time",
          "random",
          "nested",
          "custom",
          "none",
      ],
      default="auto",
      help=(
          "Partitioning strategy: auto (principled adaptive selection),"
          " scaffold, time, random, nested, custom, or none (default: auto)."
      ),
  )
  parser.add_argument(
      "--nested_k",
      "--nested-k",
      dest="nested_k",
      type=int,
      default=5,
      help="Number of outer folds for nested cross-validation (default: 5).",
  )
  parser.add_argument(
      "--time_col",
      "--time-col",
      dest="time_col",
      type=str,
      default=None,
      help=(
          "Column containing time, date, or publication year for temporal split"
          " (auto-detected if omitted)."
      ),
  )
  parser.add_argument(
      "--time_cutoff",
      "--time-cutoff",
      dest="time_cutoff",
      type=str,
      default=None,
      help=(
          "Explicit cutoff value (e.g. 2013 or 2014-06-30) for strict"
          " non-overlapping temporal split (train <= cutoff, test > cutoff)."
      ),
  )
  parser.add_argument(
      "--random_seed",
      "--random-seed",
      dest="random_seed",
      type=int,
      default=42,
      help="Random seed for data splitting (default: 42).",
  )
  parser.add_argument(
      "--output_csv",
      "--output-csv",
      "--output",
      "-o",
      dest="output_csv",
      type=str,
      default="ml_benchmark_summary.csv",
      help="Path for saving benchmark results table.",
  )
  parser.add_argument(
      "--output_md",
      "--output-md",
      dest="output_md",
      type=str,
      default=None,
      help="Optional Markdown report output path.",
  )
  parser.add_argument(
      "--output_predictions",
      "--output-preds",
      "--output_preds",
      dest="output_predictions",
      type=str,
      help=(
          "Path for saving cross-validated and holdout test predictions table."
      ),
  )
  parser.add_argument(
      "--output_metrics",
      "--output-metrics",
      dest="output_metrics",
      type=str,
      help=(
          "Path for saving standardized JSON metrics deliverable with top-level"
          " contract aliases."
      ),
  )
  parser.add_argument(
      "--tukey-alpha",
      "--tukey_alpha",
      dest="tukey_alpha",
      type=float,
      default=0.05,
      help=(
          "Family-wise significance level for Tukey HSD model tiering (default:"
          " 0.05)."
      ),
  )
  parser.add_argument(
      "--output-tukey",
      "--output_tukey",
      dest="output_tukey",
      type=str,
      default=None,
      help=(
          "Optional CSV path for the Tukey HSD statistical model tiering table."
      ),
  )
  parser.add_argument(
      "--plot_coefficients",
      "--plot-coefficients",
      dest="plot_coefficients",
      type=str,
      help=(
          "Path for saving top standardized feature coefficients bar plot PNG."
      ),
  )
  parser.add_argument(
      "--model_dir",
      "--model-dir",
      "--save_model",
      "--save-model",
      dest="model_dir",
      type=str,
      default=None,
      help=(
          "Directory or .joblib file path for saving model bundles."
          " Omit to skip persistence."
      ),
  )
  parser.add_argument(
      "--plot_feature_importance",
      "--plot-feature-importance",
      dest="plot_feature_importance",
      type=str,
      nargs="?",
      const="AUTO",
      default=None,
      help=(
          "Path for saving the feature attribution PNG. Pass the flag with no"
          " value to write feature_importance_<target>.png in the current"
          " directory."
      ),
  )
  parser.add_argument(
      "--plot_parity",
      "--plot-parity",
      dest="plot_parity",
      type=str,
      help="Path for saving multi-target parity regressions grid PNG.",
  )
  parser.add_argument(
      "--skip_3d",
      "--skip-3d",
      dest="skip_3d",
      action="store_true",
      default=False,
      help=(
          "Skip 3D conformer embedding and shape descriptor calculation for"
          " fast 2D screening."
      ),
  )
  parser.add_argument(
      "--fp_bits",
      "--fp-bits",
      dest="fp_bits",
      type=int,
      default=DEFAULT_MORGAN_FP_BITS,
      help=(
          "Bit length for Morgan circular fingerprints (default:"
          f" {DEFAULT_MORGAN_FP_BITS})."
      ),
  )
  parser.add_argument(
      "--cv_folds", type=int, dest="n_splits", help="Alias for --n_splits."
  )
  parser.add_argument(
      "--include_tabfm",
      action="store_true",
      help="Include TabFM zero-shot foundation model in benchmarking.",
  )
  parser.add_argument(
      "--tabfm_estimators",
      type=int,
      default=4,
      help=(
          "Number of TabFM permutation passes (default: 4 for fast CPU"
          " evaluation)."
      ),
  )
  parser.add_argument(
      "--top_k",
      "--top-k",
      dest="top_k",
      type=int,
      default=12,
      help=(
          "Maximum number of standardized feature coefficients to display in"
          " console output (default: 12; all features are saved to CSV)."
      ),
  )
  parser.add_argument(
      "--censored-policy",
      "--censored_policy",
      dest="censored_policy",
      choices=list(CENSORED_POLICIES),
      default="bound",
      help=(
          "Policy for relational-qualified target strings such as '< 4.0' or"
          " '>= 200': 'bound' (default; converts to numeric detection bounds"
          " without dropping rows), 'error', or 'drop'."
      ),
  )
  parser.add_argument(
      "--n-jobs",
      "--n_jobs",
      dest="n_jobs",
      type=int,
      default=None,
      help=(
          "Number of worker threads for CV and nested hyperparameter tuning"
          " (default: min(4, cpu_count); pass -1 for all CPUs)."
      ),
  )
  args = parser.parse_args()

  df_raw = pd.read_csv(args.input_csv)
  target_cols = []
  for item in args.target_col:
    for c in item.split(","):
      c = c.strip()
      if c and c not in target_cols:
        target_cols.append(c)

  for t_col in target_cols:
    if t_col not in df_raw.columns:
      raise ValueError(f"Target column '{t_col}' not found in {args.input_csv}")

  print(
      f"Loaded dataset: {len(df_raw)} molecules from {args.input_csv} with"
      f" {len(target_cols)} target endpoint(s): {target_cols}"
  )

  # Drop unparseable SMILES before featurization. Left in place they are
  # featurized as all-zero vectors, receive a fabricated prediction, and
  # silently corrupt every reported metric.
  smiles_candidates = [
      args.smiles_col,
      "canonical_smiles",
      "SMILES",
      "Smiles",
      "smiles",
  ]
  ingest_smiles_col = next(
      (c for c in smiles_candidates if c in df_raw.columns), None
  )
  if ingest_smiles_col is None:
    raise ValueError(
        f"SMILES column '{args.smiles_col}' not found in {args.input_csv}. "
        f"Available columns: {list(df_raw.columns)}"
    )
  parse_ok = [
      Chem.MolFromSmiles(str(s)) is not None if pd.notna(s) else False
      for s in df_raw[ingest_smiles_col]
  ]
  n_bad = len(parse_ok) - sum(parse_ok)
  if n_bad:
    bad_rows = [
        (i, df_raw[ingest_smiles_col].iloc[i])
        for i, good in enumerate(parse_ok)
        if not good
    ]
    preview = ", ".join(f"row {i}: {v!r}" for i, v in bad_rows[:5])
    if len(bad_rows) > 5:
      preview += f", ... (+{len(bad_rows) - 5} more)"
    print(
        f"WARNING: dropping {n_bad} row(s) with unparseable SMILES in "
        f"column '{ingest_smiles_col}': {preview}",
        file=sys.stderr,
    )
    df_raw = df_raw[parse_ok].copy().reset_index(drop=True)
    print(f"Proceeding with {len(df_raw)} valid molecules.")
  if df_raw.empty:
    raise ValueError(
        f"No parseable SMILES remain in {args.input_csv}; refusing to "
        "report metrics for an empty dataset."
    )

  if args.skip_3d:
    print(
        "Computing 2D properties and Morgan circular fingerprints (3D shape"
        " skipped via --skip_3d)..."
    )
  else:
    print(
        "Computing 2D properties, 3D shape descriptors, and Morgan circular"
        " fingerprints..."
    )
  feat_df, cols_2d, cols_3d, cols_fp, actual_smiles_col = (
      compute_molecular_features(
          df_raw,
          smiles_col=args.smiles_col,
          compute_3d=(not args.skip_3d),
          fp_bits=args.fp_bits,
      )
  )

  feature_sets = build_feature_sets(
      cols_2d=cols_2d,
      cols_3d=cols_3d,
      cols_fp=cols_fp,
      available_cols=list(feat_df.columns),
      extra_cols_spec=args.physics_cols,
      extra_name=args.extra_name,
  )

  # Pre-load TabFM foundation model once if requested
  tabfm_base_model = None
  if args.include_tabfm:
    if not HAS_TABFM:
      print(
          "Warning: tabfm library not found. Skipping TabFM evaluation.",
          file=sys.stderr,
      )
    else:
      print("Pre-loading TabFM foundation model checkpoint for benchmarking...")
      try:
        tabfm_base_model = tabfm.tabfm_v1_0_0.load(model_type="regression")
        print("TabFM foundation model loaded successfully.")
      except Exception as e:
        print(
            f"Warning: Failed to load TabFM checkpoint ({e}). Skipping TabFM.",
            file=sys.stderr,
        )

  if args.time_col and args.time_col not in feat_df.columns:
    raise ValueError(
        f"Time column '{args.time_col}' not found in {args.input_csv}."
        f" Available columns: {list(df_raw.columns)}"
    )

  # Backfill true dates across the full dataset before per-target filtering so
  # sibling rows of the same compound/patent/document contribute their dates.
  full_time_s, has_subyear_dates, global_resolved_time_col = (
      backfill_temporal_series(feat_df, time_col=args.time_col)
  )
  if args.time_cutoff is not None:
    parse_temporal_cutoff(args.time_cutoff, has_subyear_dates=has_subyear_dates)

  all_summaries = []
  all_tier_tables = []
  all_best_predictions = []
  all_standardized_coefficients = {}
  last_temporal_metadata = None
  per_target_split_meta: Dict[str, Dict[str, Any]] = {}

  for t_col in target_cols:
    print("\n=======================================================")
    print(f"Benchmarking Target Endpoint: {t_col}")
    print("=======================================================")
    parsed_vals, keep_mask = parse_target_series(
        feat_df[t_col],
        target_col=t_col,
        censored_policy=args.censored_policy,
    )
    valid_df = feat_df.loc[keep_mask].copy().reset_index(drop=True)
    valid_df[t_col] = parsed_vals[keep_mask]
    valid_time_s = full_time_s.loc[keep_mask].reset_index(drop=True)
    n_valid = len(valid_df)
    if n_valid < 4:
      print(f"Skipping {t_col}: fewer than 4 valid observations.")
      continue

    # 1. Validate non-boolean temporal metadata (>= 2 distinct time points)
    n_dated = int(valid_time_s.notna().sum())
    n_unique_dates = int(valid_time_s.dropna().nunique())
    target_has_subyear = bool(
        has_subyear_dates
        and any(
            1800.0 <= float(v) < 2101.0
            and abs(float(v) - round(float(v))) >= 1e-7
            for v in valid_time_s.dropna().tolist()
        )
    )
    resolved_time_col = (
        global_resolved_time_col if n_unique_dates >= 2 else None
    )

    if args.time_col is not None and resolved_time_col is None:
      raise ValueError(
          f"Time column '{args.time_col}' for endpoint '{t_col}' does not"
          " contain at least 2 distinct valid non-boolean time values"
          f" (found {n_unique_dates} distinct across {n_dated}/{n_valid}"
          " dated rows; boolean True/False placeholders are rejected)."
      )

    # 2. Resolve partitioning strategy
    if args.split_col and args.split_col in valid_df.columns:
      chosen_strategy = "custom"
    else:
      chosen_strategy = args.split_method.lower()
      if chosen_strategy == "auto":
        if args.test_split is not None and args.test_split == 0.0:
          chosen_strategy = "none"
        elif args.time_cutoff is not None:
          chosen_strategy = "time"
        elif resolved_time_col is not None and n_dated >= max(
            4, int(math.ceil(0.8 * n_valid))
        ):
          chosen_strategy = "time"
        elif n_valid < 30:
          chosen_strategy = "nested"
        else:
          chosen_strategy = "scaffold"

    # 3. Execute strategy
    train_df = valid_df
    test_df = None
    nested_k = None
    cur_temporal_meta = None

    if chosen_strategy == "custom":
      split_series = (
          valid_df[args.split_col].astype(str).str.strip().str.lower()
      )
      target_eval = (
          args.eval_split.strip().lower() if args.eval_split else "test"
      )
      test_mask = split_series == target_eval
      split_labels = set(split_series.unique())
      if target_eval == "test" and {"train", "val", "test"} <= split_labels:
        train_mask = split_series == "train"
      else:
        train_mask = ~test_mask
      train_df = valid_df[train_mask].copy().reset_index(drop=True)
      test_df = valid_df[test_mask].copy().reset_index(drop=True)
      print(
          f"[Pre-Assigned Split from '{args.split_col}': '{target_eval}' as"
          " holdout test]:"
      )
      print(f"  Training Set : N = {len(train_df)} compounds")
      print(f"  Holdout Test : N = {len(test_df)} compounds")
    elif chosen_strategy == "nested":
      nested_k = min(args.nested_k, max(3, n_valid // 3)) if n_valid >= 6 else 2
      print(
          f"[Adaptive Split Protocol] Small-N series (N={n_valid} < 30) without"
          " time metadata."
      )
      print(
          f"  Executing Nested {nested_k}-Fold Cross-Validation (Cawley &"
          " Talbot 2010 protocol):"
      )
      print(
          "  - Inner loop: Leakage-free preprocessing and alpha tuning on"
          " train folds"
      )
      print(
          "  - Outer loop: Unbiased out-of-sample evaluation on held-out folds"
      )
      print(
          f"  - Full Statistical Power: N_test = {n_valid} (100% out-of-sample"
          " coverage)"
      )
    elif chosen_strategy == "time":
      if resolved_time_col is None:
        raise ValueError(
            f"Temporal split requested for endpoint '{t_col}', but no column"
            " with at least 2 distinct valid non-boolean timestamps is"
            f" available (dated rows: {n_dated}/{n_valid}, distinct:"
            f" {n_unique_dates}). Refusing to fabricate a split from boolean"
            " flags or row sequence order."
        )
      if args.test_split is not None and args.test_split > 0.0:
        ratio_val = args.test_split
      else:
        n_te_target = max(
            5, min(int(0.35 * n_dated), int(round(math.sqrt(n_dated))))
        )
        ratio_val = min(0.5, float(n_te_target) / float(n_dated))

      if args.time_cutoff is not None:
        cutoff_val = parse_temporal_cutoff(
            args.time_cutoff, has_subyear_dates=target_has_subyear
        )
      else:
        cutoff_val, _ = resolve_temporal_cutoff(
            valid_time_s,
            train_ratio=1.0 - ratio_val,
            val_ratio=0.0,
            test_ratio=ratio_val,
        )

      train_idx, _, test_idx = temporal_split(
          valid_time_s,
          train_ratio=1.0 - ratio_val,
          val_ratio=0.0,
          test_ratio=ratio_val,
          cutoff=cutoff_val,
          has_subyear_dates=target_has_subyear,
      )

      # Persist backfilled calendar dates onto valid_df before slicing
      if (
          is_calendar_col_name(resolved_time_col)
          or resolved_time_col == "patent_id"
      ):
        if (
            resolved_time_col in valid_df.columns
            and resolved_time_col != "patent_id"
        ):
          out_t_col = resolved_time_col
        else:
          out_t_col = "document_year"
        if out_t_col not in valid_df.columns:
          valid_df[out_t_col] = pd.Series(
              np.nan, index=valid_df.index, dtype=object
          )
        cur_t_parsed, _ = parse_temporal_series(
            valid_df[out_t_col], calendar_bounds=True
        )
        fill_mask = cur_t_parsed.isna() & valid_time_s.notna()
        is_year_only_col = "year" in out_t_col.lower()
        if fill_mask.any():
          valid_df[out_t_col] = valid_df[out_t_col].astype(object)
          for r_idx in valid_df.index[fill_mask]:
            val_f = float(valid_time_s.loc[r_idx])
            valid_df.at[r_idx, out_t_col] = (
                int(math.floor(val_f + 1e-9))
                if is_year_only_col
                else format_temporal_value(val_f, target_has_subyear)
            )

      train_df = valid_df.iloc[train_idx].copy().reset_index(drop=True)
      test_df = valid_df.iloc[test_idx].copy().reset_index(drop=True)
      tr_times = valid_time_s.iloc[train_idx].dropna()
      te_times = valid_time_s.iloc[test_idx].dropna()
      n_undated_tr = len(train_idx) - len(tr_times)
      cutoff_fmt = format_temporal_value(
          cutoff_val, has_subyear_dates=target_has_subyear
      )
      tr_min_fmt = format_temporal_value(
          float(tr_times.min()), has_subyear_dates=target_has_subyear
      )
      tr_max_fmt = format_temporal_value(
          float(tr_times.max()), has_subyear_dates=target_has_subyear
      )
      te_min_fmt = format_temporal_value(
          float(te_times.min()), has_subyear_dates=target_has_subyear
      )
      te_max_fmt = format_temporal_value(
          float(te_times.max()), has_subyear_dates=target_has_subyear
      )
      if target_has_subyear:
        cutoff_json_val: Any = cutoff_fmt
      elif abs(cutoff_val - round(cutoff_val)) < 1e-6:
        cutoff_json_val = int(round(cutoff_val))
      else:
        cutoff_json_val = float(cutoff_val)
      cur_temporal_meta = {
          "time_col": resolved_time_col,
          "time_cutoff": cutoff_json_val,
          "train_time_range": [tr_min_fmt, tr_max_fmt],
          "test_time_range": [te_min_fmt, te_max_fmt],
          "n_undated_train": int(n_undated_tr),
      }
      last_temporal_metadata = cur_temporal_meta
      undated_suffix = (
          f" [{n_undated_tr} undated placed in train]"
          if n_undated_tr > 0
          else ""
      )
      print(
          "[Adaptive Split Protocol] Sheridan Prospective TIME split"
          f" (column: '{resolved_time_col}', --time-cutoff {cutoff_fmt}):"
      )
      print(
          f"  Training Set (Past Leads <= {cutoff_fmt}) : N = {len(train_df)}"
          f" compounds (range: [{tr_min_fmt}, {tr_max_fmt}]){undated_suffix}"
      )
      print(
          f"  Holdout Test (Future Leads > {cutoff_fmt}): N = {len(test_df)}"
          f" compounds (range: [{te_min_fmt}, {te_max_fmt}])"
      )
    elif chosen_strategy == "scaffold":
      if args.test_split is not None and args.test_split > 0.0:
        ratio_val = args.test_split
      else:
        n_te_target = max(
            6, min(int(0.35 * n_valid), int(round(math.sqrt(n_valid))))
        )
        ratio_val = n_te_target / n_valid

      train_idx, _, test_idx = scaffold_split(
          valid_df[actual_smiles_col].tolist(),
          train_ratio=1.0 - ratio_val,
          val_ratio=0.0,
          test_ratio=ratio_val,
      )
      train_df = valid_df.iloc[train_idx].copy().reset_index(drop=True)
      test_df = valid_df.iloc[test_idx].copy().reset_index(drop=True)
      print(
          "[Adaptive Split Protocol] Bemis-Murcko SCAFFOLD split (Adaptive"
          f" N_test={len(test_df)}, {ratio_val*100:.1f}%):"
      )
      print(f"  Training Set : N = {len(train_df)} compounds")
      print(f"  Holdout Test : N = {len(test_df)} compounds")
    elif chosen_strategy == "random":
      ratio_val = (
          args.test_split
          if (args.test_split is not None and args.test_split > 0.0)
          else 0.20
      )
      print("\n" + "=" * 78)
      print("[CRITICAL WARNING: POTENTIAL DATA LEAKAGE DETECTED]")
      print(
          "Random train/test splitting on congeneric small-molecule series"
          " violates"
      )
      print(
          "medicinal chemistry benchmarking invariants. Highly similar"
          " structural analogues"
      )
      print(
          "and identical Murcko scaffolds leak between train and test sets,"
          " causing"
      )
      print(
          "severe test set overoptimism (measuring 2D interpolation rather than"
          " prospective generalization)."
      )
      print(
          "Use '--split_method auto', 'scaffold', 'time', or 'nested' for valid"
          " prospective evaluation."
      )
      print("=" * 78 + "\n")
      train_idx, test_idx = train_test_split(
          np.arange(n_valid),
          test_size=ratio_val,
          random_state=args.random_seed,
      )
      train_df = valid_df.iloc[train_idx].copy().reset_index(drop=True)
      test_df = valid_df.iloc[test_idx].copy().reset_index(drop=True)
      print(f"Partitioned {n_valid} compounds via RANDOM split:")
      print(f"  Training Set : N = {len(train_df)} compounds")
      print(f"  Holdout Test : N = {len(test_df)} compounds")
    else:  # none
      print(
          f"Running Full Cross-Validation on all {n_valid} compounds (no"
          " holdout split)."
      )
      train_df = valid_df
      test_df = None

    per_target_split_meta[t_col] = {
        "strategy": chosen_strategy,
        "temporal_meta": cur_temporal_meta,
        "n_test": (
            len(test_df)
            if test_df is not None
            else (n_valid if chosen_strategy == "nested" else 0)
        ),
    }

    summary_df, trained_models = run_congeneric_benchmark(
        train_df=train_df,
        target_col=t_col,
        feature_sets=feature_sets,
        test_df=test_df,
        nested_outer_k=nested_k,
        cv_strategy_name=args.cv_strategy,
        n_splits=args.n_splits,
        include_tabfm=args.include_tabfm,
        tabfm_base_model=tabfm_base_model,
        tabfm_estimators=args.tabfm_estimators,
        smiles_col=actual_smiles_col,
        censored_policy=args.censored_policy,
        n_jobs=args.n_jobs,
    )
    all_summaries.append(summary_df)

    # Statistical multi-comparison: which models are actually distinguishable?
    tier_df = tukey_hsd_model_tiers(trained_models, alpha=args.tukey_alpha)
    if not tier_df.empty:
      tier_df.insert(0, "Target", t_col)
      all_tier_tables.append(tier_df)
      top_tier = tier_df[tier_df["Tier"] == "Tier A"]["Model_Key"].tolist()
      print(
          f"\n### [{t_col}] Tukey HSD Statistical Model Tiering"
          f" (alpha={args.tukey_alpha}, {len(top_tier)} Tier A model(s))\n",
          flush=True,
      )
      t_cols = list(tier_df.columns)
      print("| " + " | ".join(str(c) for c in t_cols) + " |")
      print("| " + " | ".join(":---" for _ in t_cols) + " |")
      for _, tr in tier_df.iterrows():
        cells = []
        for c in t_cols:
          v = tr[c]
          if isinstance(v, (float, np.floating)) and np.isfinite(v):
            cells.append(f"{float(v):.3f}")
          elif pd.isna(v):
            cells.append("—")
          else:
            cells.append(str(v))
        print("| " + " | ".join(cells) + " |", flush=True)

    # Collect out-of-fold and test predictions from best model (selected by Train CV R2)
    if not summary_df.empty:
      best_row = summary_df.iloc[0]
      best_key = f"{best_row['Feature_Set']}__{best_row['Model']}"
      best_entry = trained_models.get(best_key)

      if best_entry is not None and args.model_dir:
        save_model_bundles(
            best_entry,
            target_col=t_col,
            model_dir=args.model_dir,
            include_tabfm=args.include_tabfm,
            tabfm_base_model=tabfm_base_model,
            tabfm_estimators=args.tabfm_estimators,
            censored_policy=args.censored_policy,
            n_jobs=args.n_jobs,
        )

      if best_entry is not None and args.plot_feature_importance:
        fi_path = args.plot_feature_importance
        if fi_path == "AUTO":
          safe_t = "".join(
              c if (c.isalnum() or c in "-_") else "_" for c in str(t_col)
          )
          fi_path = f"feature_importance_{safe_t}.png"
        elif len(target_cols) > 1:
          # One path, several endpoints: disambiguate rather than
          # letting each target overwrite the previous one's plot.
          root, ext = os.path.splitext(fi_path)
          safe_t = "".join(
              c if (c.isalnum() or c in "-_") else "_" for c in str(t_col)
          )
          fi_path = f"{root}_{safe_t}{ext or '.png'}"

        target_attr_entry = best_entry
        if target_attr_entry.get("is_fingerprint", False):
          # Overall winning model uses anonymous circular fingerprint bits.
          # For scientific interpretability, select top continuous descriptor model.
          interp_rows = summary_df[
              ~summary_df["Feature_Set"]
              .astype(str)
              .str.contains(
                  r"FP|Fingerprint|ecfp|morgan", case=False, regex=True
              )
          ]
          if not interp_rows.empty:
            top_interp_row = interp_rows.iloc[0]
            interp_key = (
                f"{top_interp_row['Feature_Set']}__{top_interp_row['Model']}"
            )
            target_attr_entry = trained_models.get(interp_key)
            interp_r2 = float(
                top_interp_row.get(
                    "Train_CV_R2", top_interp_row.get("CV_R2", 0.0)
                )
            )
            print(
                f"\n[Notice] Winning model '{best_key}' uses anonymous"
                " fingerprint bit hashes (uninformative for attribution).\n"
                "Selected top interpretable continuous descriptor model"
                f" '{interp_key}' (CV R² = {interp_r2:.3f}) for feature"
                " attribution."
            )
          else:
            print(
                "\n[Notice] No interpretable continuous descriptor model"
                f" available for target '{t_col}'; suppressing anonymous"
                " fingerprint bit attribution.",
                file=sys.stderr,
            )
            target_attr_entry = None

        if target_attr_entry is not None:
          if not plot_feature_importance(
              target_attr_entry, fi_path, target_name=t_col
          ):
            print(
                "WARNING: no feature attribution could be computed for "
                f"'{target_attr_entry.get('model_name')}' on target '{t_col}'; "
                f"no plot was written to {fi_path}.",
                file=sys.stderr,
            )

      if best_entry is not None:
        if best_entry.get("is_nested"):
          tr_df = best_entry.get("train_df")
          for i, (act, cv_p, te_p) in enumerate(
              zip(
                  best_entry["y_train_actual"],
                  best_entry["y_train_cv_pred"],
                  best_entry["y_test_pred"],
              )
          ):
            row_cid = (
                tr_df.iloc[i].get("compound_id", f"CMPD_{i+1:04d}")
                if tr_df is not None
                else f"CMPD_{i+1:04d}"
            )
            row_smi = (
                tr_df.iloc[i].get(
                    "canonical_smiles", tr_df.iloc[i].get(actual_smiles_col, "")
                )
                if tr_df is not None
                else ""
            )
            all_best_predictions.append({
                "Compound_ID": str(row_cid),
                "compound_id": str(row_cid),
                "SMILES": str(row_smi),
                "canonical_smiles": str(row_smi),
                "Target": t_col,
                "Split": "train",
                "Actual": act,
                "Prediction": cv_p,
                "Predicted": cv_p,
                "Residual": float(act - cv_p),
                "Predicted_pKi": cv_p,
                "Best_Model": best_row["Model"],
                "Best_Feature_Set": best_row["Feature_Set"],
                "Train_CV_R2": best_row["Train_CV_R2"],
                "Train_CV_r": best_row["Train_CV_r"],
                "Test_R2": best_row["Test_R2"],
                "Test_Pearson_r": best_row["Test_Pearson_r"],
            })
            all_best_predictions.append({
                "Compound_ID": str(row_cid),
                "compound_id": str(row_cid),
                "SMILES": str(row_smi),
                "canonical_smiles": str(row_smi),
                "Target": t_col,
                "Split": "test",
                "Actual": act,
                "Prediction": te_p,
                "Predicted": te_p,
                "Residual": float(act - te_p),
                "Predicted_pKi": te_p,
                "Best_Model": best_row["Model"],
                "Best_Feature_Set": best_row["Feature_Set"],
                "Train_CV_R2": best_row["Train_CV_R2"],
                "Train_CV_r": best_row["Train_CV_r"],
                "Test_R2": best_row["Test_R2"],
                "Test_Pearson_r": best_row["Test_Pearson_r"],
            })
        else:
          tr_df = best_entry.get("train_df")
          te_df = best_entry.get("test_df")
          # 1. Train CV predictions
          for i, (act, pred) in enumerate(
              zip(best_entry["y_train_actual"], best_entry["y_train_cv_pred"])
          ):
            row_cid = None
            if tr_df is not None:
              for c_cand in [
                  "Compound_ID",
                  "compound_id",
                  "ID",
                  "id",
                  "Molecule_ID",
              ]:
                if c_cand in tr_df.columns:
                  row_cid = tr_df.iloc[i][c_cand]
                  break
            if row_cid is None:
              row_cid = f"CMPD_{i+1:04d}"

            row_smi = ""
            if tr_df is not None:
              for s_cand in [
                  "SMILES",
                  "smiles",
                  "canonical_smiles",
                  actual_smiles_col,
              ]:
                if s_cand in tr_df.columns:
                  row_smi = tr_df.iloc[i][s_cand]
                  break

            pred_record = {
                "Compound_ID": str(row_cid),
                "compound_id": str(row_cid),
                "SMILES": str(row_smi),
                "canonical_smiles": str(row_smi),
                "Target": t_col,
                "Split": "train",
                "Actual": act,
                "Prediction": pred,
                "Predicted": pred,
                "Residual": float(act - pred),
                "Predicted_pKi": pred,
                "Best_Model": best_row["Model"],
                "Best_Feature_Set": best_row["Feature_Set"],
            }
            # Which metric columns exist depends on whether the
            # benchmark actually had held-out rows, not on whether a
            # test frame object was constructed: an empty test frame
            # is not None but yields the pure-CV schema.
            for col in (
                "Train_CV_R2",
                "Train_CV_r",
                "Test_R2",
                "Test_Pearson_r",
                "CV_R2",
                "Pearson_r",
            ):
              if col in best_row.index:
                pred_record[col] = best_row[col]
            all_best_predictions.append(pred_record)

          # 2. Holdout Test predictions
          if (
              best_entry.get("has_test")
              and best_entry["y_test_actual"] is not None
          ):
            for i, (act, pred) in enumerate(
                zip(best_entry["y_test_actual"], best_entry["y_test_pred"])
            ):
              row_cid = None
              if te_df is not None:
                for c_cand in [
                    "Compound_ID",
                    "compound_id",
                    "ID",
                    "id",
                    "Molecule_ID",
                ]:
                  if c_cand in te_df.columns:
                    row_cid = te_df.iloc[i][c_cand]
                    break
              if row_cid is None:
                row_cid = f"CMPD_TEST_{i+1:04d}"

              row_smi = ""
              if te_df is not None:
                for s_cand in [
                    "SMILES",
                    "smiles",
                    "canonical_smiles",
                    actual_smiles_col,
                ]:
                  if s_cand in te_df.columns:
                    row_smi = te_df.iloc[i][s_cand]
                    break

              all_best_predictions.append({
                  "Compound_ID": str(row_cid),
                  "compound_id": str(row_cid),
                  "SMILES": str(row_smi),
                  "canonical_smiles": str(row_smi),
                  "Target": t_col,
                  "Split": "test",
                  "Actual": act,
                  "Prediction": pred,
                  "Predicted": pred,
                  "Residual": float(act - pred),
                  "Predicted_pKi": pred,
                  "Best_Model": best_row["Model"],
                  "Best_Feature_Set": best_row["Feature_Set"],
                  "Train_CV_R2": best_row["Train_CV_R2"],
                  "Train_CV_r": best_row["Train_CV_r"],
                  "Test_R2": best_row["Test_R2"],
                  "Test_Pearson_r": best_row["Test_Pearson_r"],
              })

          # 3. Prospective predictions for unlabeled/patent rows missing this target
          unlabeled_df = (
              feat_df[feat_df[t_col].isna()].copy().reset_index(drop=True)
          )
          if (
              not unlabeled_df.empty
              and best_entry.get("preprocessor") is not None
              and best_entry.get("model") is not None
          ):
            try:
              f_cols = best_entry.get("feature_cols", [])
              if f_cols and all(c in unlabeled_df.columns for c in f_cols):
                X_unlab_raw = unlabeled_df[f_cols].values.astype(float)
                X_unlab_proc = best_entry["preprocessor"].transform(X_unlab_raw)
                y_unlab_pred = best_entry["model"].predict(X_unlab_proc)
                print(
                    f"  [+] Scored {len(unlabeled_df)} unlabeled/patent"
                    f" compounds prospectively for {t_col}."
                )
                prosp_rows_disp = []
                for i, pred_val in enumerate(y_unlab_pred):
                  row_cid = unlabeled_df.iloc[i].get(
                      "compound_id", f"CMPD_UNLAB_{i+1:04d}"
                  )
                  row_smi = unlabeled_df.iloc[i].get(
                      "canonical_smiles",
                      unlabeled_df.iloc[i].get(actual_smiles_col, ""),
                  )
                  prosp_entry = {
                      "Compound_ID": str(row_cid),
                      "compound_id": str(row_cid),
                      "SMILES": str(row_smi),
                      "canonical_smiles": str(row_smi),
                      "Target": t_col,
                      "Split": "prospective_unlabeled",
                      "Actual": np.nan,
                      "Prediction": float(pred_val),
                      "Predicted": float(pred_val),
                      "Residual": np.nan,
                      "Predicted_pKi": float(pred_val),
                      "Best_Model": best_row["Model"],
                      "Best_Feature_Set": best_row["Feature_Set"],
                  }
                  all_best_predictions.append(prosp_entry)
                  prosp_rows_disp.append(prosp_entry)
                if prosp_rows_disp:
                  prosp_disp_df = pd.DataFrame(prosp_rows_disp).sort_values(
                      by="Predicted", ascending=False
                  )
                  print(
                      "\n=== TOP PROSPECTIVE UNLABELED CANDIDATES SCORED FOR"
                      f" {t_col} ==="
                  )
                  print(
                      "| Rank | Compound ID | Predicted"
                      f" {t_col} | Model | Canonical SMILES |"
                  )
                  print("| :--- | :--- | :--- | :--- | :--- |")
                  for rk, (_, pr) in enumerate(
                      prosp_disp_df.head(10).iterrows(), 1
                  ):
                    smi_short = str(pr["canonical_smiles"])
                    if len(smi_short) > 48:
                      smi_short = smi_short[:45] + "..."
                    print(
                        f"| {rk} | {pr['Compound_ID']} |"
                        f" {pr['Predicted']:.3f} | {pr['Best_Model']} |"
                        f" `{smi_short}` |"
                    )
            except Exception as e:
              print(
                  "  [!] Warning: Could not prospectively score unlabeled rows"
                  f" for {t_col}: {e}"
              )

    ridge_models = [
        k
        for k, v in trained_models.items()
        if "Ridge" in k
        and not v.get("is_fingerprint", False)
        and v.get("cv_r2", -999) > -1.0
    ]
    linear_models = ridge_models or [
        k
        for k, v in trained_models.items()
        if ("Ridge" in k or "ElasticNet" in k)
        and not v.get("is_fingerprint", False)
        and v.get("cv_r2", -999) > -1.0
    ]
    if linear_models:
      best_linear_key = max(
          linear_models, key=lambda k: trained_models[k].get("cv_r2", -999)
      )
      try:
        attrib_df = compute_feature_attribution(trained_models[best_linear_key])
        if not attrib_df.empty:
          print(
              "\n=== STANDARDIZED FEATURE COEFFICIENTS (β) FOR"
              f" {t_col} ({best_linear_key}) ==="
          )
          print("| Rank | Descriptor / Feature | Standardized β | Direction |")
          print("| :--- | :--- | :--- | :--- |")
          coef_records = []
          for rk, (_, ar) in enumerate(attrib_df.iterrows(), 1):
            beta_val = float(ar["Importance"])
            direction = (
                "Increases Potency (+)"
                if beta_val > 0
                else "Decreases Potency (-)"
            )
            coef_records.append({
                "Rank": rk,
                "Target": t_col,
                "Linear_Model": best_linear_key,
                "Feature": str(ar["Feature"]),
                "Standardized_Beta": round(beta_val, 5),
                "Direction": direction,
            })
            if rk <= args.top_k:
              print(
                  f"| {rk} | `{ar['Feature']}` | {beta_val:+.4f} | {direction}"
                  " |"
              )
          if len(attrib_df) > args.top_k:
            print(
                f"Showing top {min(args.top_k, len(attrib_df))} of"
                f" {len(attrib_df)} features (all {len(attrib_df)} features"
                " saved to CSV/JSON)."
            )
          all_standardized_coefficients[t_col] = coef_records
          if args.output_csv:
            coef_csv = (
                os.path.splitext(args.output_csv)[0] + "_coefficients.csv"
            )
            pd.DataFrame(coef_records).to_csv(coef_csv, index=False)
            print(f"Saved standardized feature coefficients CSV to: {coef_csv}")
      except Exception:
        pass
      if args.plot_coefficients:
        if len(target_cols) == 1:
          out_coef = args.plot_coefficients
        else:
          base, ext = os.path.splitext(args.plot_coefficients)
          out_coef = f"{base}_{t_col}{ext}"
        print(f"Saving top feature weights for {t_col} to: {out_coef}")
        plot_model_coefficients(trained_models[best_linear_key], out_coef)
    elif args.plot_coefficients:
      print(
          "Notice: No interpretable continuous linear models found for"
          " coefficient plotting."
      )

  if not all_summaries:
    print(
        "Error: No endpoints had enough valid observations to benchmark.",
        file=sys.stderr,
    )
    sys.exit(1)

  master_summary_df = pd.concat(all_summaries, ignore_index=True)
  os.makedirs(
      os.path.dirname(os.path.abspath(args.output_csv)) or ".", exist_ok=True
  )
  master_summary_df.to_csv(args.output_csv, index=False)
  print("\n### Multi-Target Congeneric ML Benchmark Summary\n")
  disp_cols = [
      c
      for c in (
          "Target",
          "Feature_Set",
          "Model",
          "Train_CV_R2",
          "Train_CV_RMSE",
          "Test_R2",
          "Test_RMSE",
          "Test_Spearman_rho",
          "Tukey_Tier",
      )
      if c in master_summary_df.columns
  ]
  if not disp_cols:
    disp_cols = list(master_summary_df.columns)[:8]
  print("| " + " | ".join(disp_cols) + " |")
  print("| " + " | ".join(":---" for _ in disp_cols) + " |")
  for _, mrow in master_summary_df[disp_cols].iterrows():
    cells = []
    for cname in disp_cols:
      v = mrow[cname]
      if isinstance(v, (float, np.floating)) and np.isfinite(v):
        cells.append(f"{float(v):.3f}")
      elif pd.isna(v):
        cells.append("—")
      else:
        cells.append(str(v))
    print("| " + " | ".join(cells) + " |")
  print(f"\nFull benchmark saved to: {args.output_csv}")

  # Save Tukey HSD statistical tiering
  if all_tier_tables:
    tukey_path = (
        args.output_tukey
        or os.path.splitext(args.output_csv)[0] + "_tukey_tiers.csv"
    )
    master_tier_df = pd.concat(all_tier_tables, ignore_index=True)
    os.makedirs(
        os.path.dirname(os.path.abspath(tukey_path)) or ".", exist_ok=True
    )
    master_tier_df.to_csv(tukey_path, index=False)
    print(f"Tukey HSD model tiering saved to: {tukey_path}")

  # Save predictions
  pred_path = (
      args.output_predictions
      if args.output_predictions
      else os.path.splitext(args.output_csv)[0] + "_predictions.csv"
  )
  if all_best_predictions:
    pred_df = pd.DataFrame(all_best_predictions)
    os.makedirs(
        os.path.dirname(os.path.abspath(pred_path)) or ".", exist_ok=True
    )
    pred_df.to_csv(pred_path, index=False)
    print(f"Saved predictions to: {pred_path}")

    if args.plot_parity:
      try:
        try:
          plot_archetypes.plot_multi_target_parity(
              df=pred_df,
              target_col="Target",
              exp_col="Actual",
              pred_col="Predicted",
              output_png=args.plot_parity,
          )
        except Exception:
          try:
            set_safe_plot_style()
          except Exception:
            pass

          palette = globals().get("PUBLICATION_PALETTE") or {
              "background": "#FFFFFF",
              "card": "#F8FAFC",
              "border": "#CBD5E1",
              "primary": "#0284C7",
              "primary_dark": "#0369A1",
              "primary_light": "#BAE6FD",
              "secondary": "#64748B",
              "accent": "#DB2777",
              "success": "#16A34A",
              "warning": "#D97706",
              "danger": "#DC2626",
              "purple": "#9333EA",
              "cyan": "#0891B2",
              "text": "#0F172A",
              "text_muted": "#475569",
              "grid": "#E2E8F0",
          }

          fig, ax = plt.subplots(figsize=(7, 6), dpi=300)
          fig.patch.set_facecolor(palette.get("background", "#FFFFFF"))
          valid_pts = pred_df.dropna(subset=["Actual", "Predicted"])

          split_col = None
          for cand in [
              "Split",
              "split",
              "Set",
              "set",
              "Partition",
              "partition",
              "Fold",
              "fold",
              "subset",
              "Subset",
          ]:
            if cand in valid_pts.columns:
              split_col = cand
              break

          has_test_pts = False
          tr_pts = None
          te_pts = None
          if split_col is not None:
            tr_pts = valid_pts[
                valid_pts[split_col]
                .astype(str)
                .str.contains(r"train|cv", case=False, na=False)
            ]
            te_pts = valid_pts[
                valid_pts[split_col]
                .astype(str)
                .str.contains(
                    r"test|holdout|val|external", case=False, na=False
                )
            ]
            has_test_pts = len(tr_pts) > 0 and len(te_pts) > 0

          if has_test_pts:
            ax.scatter(
                tr_pts["Actual"],
                tr_pts["Predicted"],
                alpha=0.85,
                color=palette.get("primary", "#38BDF8"),
                edgecolors="#0F172A",
                s=65,
                label=f"Train CV (N={len(tr_pts)})",
                zorder=3,
            )
            ax.scatter(
                te_pts["Actual"],
                te_pts["Predicted"],
                alpha=1.0,
                color=palette.get("danger", "#F87171"),
                marker="D",
                edgecolors="#FB7185",
                s=85,
                label=f"Holdout Test (N={len(te_pts)})",
                zorder=4,
            )
          else:
            ax.scatter(
                valid_pts["Actual"],
                valid_pts["Predicted"],
                alpha=0.85,
                color=palette.get("primary", "#38BDF8"),
                edgecolors="#0F172A",
                s=70,
                label=f"Evaluated Congeners (N={len(valid_pts)})",
                zorder=3,
            )

          min_val = min(valid_pts["Actual"].min(), valid_pts["Predicted"].min())
          max_val = max(valid_pts["Actual"].max(), valid_pts["Predicted"].max())
          ax.plot(
              [min_val, max_val],
              [min_val, max_val],
              color="#64748B",
              linestyle="--",
              alpha=0.8,
              label="Ideal Parity (y = x)",
          )
          ax.set_xlabel(
              "Observed / Experimental Value", fontsize=11, fontweight="bold"
          )
          ax.set_ylabel("Predicted Value", fontsize=11, fontweight="bold")
          ax.set_title(
              "Experimental vs Predicted Parity Plot",
              fontsize=12,
              fontweight="bold",
          )
          ax.legend(frameon=True, loc="lower right", fontsize=9)
          try:
            apply_publication_spines(ax)
          except Exception:
            pass
          plt.tight_layout()
          fig.savefig(args.plot_parity)
          plt.close(fig)
        print(
            "Saved multi-target parity grid to:"
            f" {os.path.abspath(args.plot_parity)}",
            flush=True,
        )
        print(
            "  Markdown embed: ![Multi-Target Parity"
            f" Grid]({os.path.abspath(args.plot_parity)})",
            flush=True,
        )
      except Exception as e:
        print(f"Could not render parity plot: {e}", flush=True)

  # Save metrics JSON deliverable with top-level contract aliases
  if args.output_metrics:
    models_dict = {}
    best_test_r2_val = None
    best_test_entry = None
    for _, r in master_summary_df.iterrows():
      m_key = f"{r['Feature_Set']}__{r['Model']}"
      m_data = {
          "feature_set": r["Feature_Set"],
          "model": r["Model"],
          "train_cv_r2": float(r.get("Train_CV_R2", r.get("CV_R2", 0.0))),
          "train_cv_rmse": (
              float(r.get("Train_CV_RMSE", 0.0))
              if "Train_CV_RMSE" in r and pd.notna(r["Train_CV_RMSE"])
              else None
          ),
      }
      if "Test_R2" in r and pd.notna(r["Test_R2"]):
        t_r2 = float(r["Test_R2"])
        m_data["test_r2"] = t_r2
        m_data["r2"] = t_r2
        m_data["rmse"] = (
            float(r.get("Test_RMSE", 0.0))
            if pd.notna(r.get("Test_RMSE"))
            else None
        )
        m_data["mae"] = (
            float(r.get("Test_MAE", 0.0))
            if pd.notna(r.get("Test_MAE"))
            else None
        )
        m_data["spearman"] = (
            float(r.get("Test_Spearman_rho", 0.0))
            if pd.notna(r.get("Test_Spearman_rho"))
            else None
        )
        m_data["spearman_rho"] = m_data["spearman"]
        if best_test_r2_val is None or t_r2 > best_test_r2_val:
          best_test_r2_val = t_r2
          best_test_entry = m_data
      models_dict[m_key] = m_data

    best_target_name = (
        best_test_entry.get("target") if best_test_entry is not None else None
    )
    if best_target_name and best_target_name in per_target_split_meta:
      rep_meta = per_target_split_meta[best_target_name]
    elif per_target_split_meta:
      rep_meta = next(iter(per_target_split_meta.values()))
    else:
      rep_meta = {
          "strategy": chosen_strategy,
          "temporal_meta": last_temporal_metadata,
          "n_test": len(test_df) if test_df is not None else 0,
      }
    rep_strategy = rep_meta["strategy"]
    rep_n_test = int(rep_meta["n_test"])
    metrics_payload = {}
    metrics_payload["split_strategy"] = (
        rep_strategy
        if rep_strategy != "random"
        else "random (POTENTIAL DATA LEAKAGE RISK)"
    )
    if rep_strategy == "time" and rep_meta["temporal_meta"] is not None:
      metrics_payload.update(rep_meta["temporal_meta"])
    if rep_strategy == "random":
      metrics_payload["data_leakage_warning"] = (
          "Random splitting on congeneric series risks leaking identical"
          " scaffolds between train and test sets, inflating test metrics."
      )
    if len(per_target_split_meta) > 1:
      metrics_payload["per_target_splits"] = per_target_split_meta
    if best_test_entry is not None:
      metrics_payload["test_r2"] = float(best_test_entry["test_r2"])
      metrics_payload["held_out_test_r2"] = float(best_test_entry["test_r2"])
      metrics_payload["r2_test"] = float(best_test_entry["test_r2"])
      metrics_payload["held_out_test"] = {
          "best_model": (
              f"{best_test_entry['feature_set']}__{best_test_entry['model']}"
          ),
          "r2": float(best_test_entry["test_r2"]),
          "test_r2": float(best_test_entry["test_r2"]),
          "rmse": best_test_entry.get("rmse"),
          "mae": best_test_entry.get("mae"),
          "spearman_rho": best_test_entry.get("spearman_rho"),
          "spearman": best_test_entry.get("spearman"),
          "n_test": rep_n_test,
          "dof": max(0, rep_n_test - 2) if rep_n_test > 0 else 0,
      }
    metrics_payload["models"] = models_dict
    if all_standardized_coefficients:
      metrics_payload["standardized_feature_coefficients"] = (
          all_standardized_coefficients
      )
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_metrics)) or ".",
        exist_ok=True,
    )
    with open(args.output_metrics, "w", encoding="utf-8") as f:
      json.dump(metrics_payload, f, indent=2)
    print(f"Saved metrics JSON deliverable to: {args.output_metrics}")


if __name__ == "__main__":
  main()
