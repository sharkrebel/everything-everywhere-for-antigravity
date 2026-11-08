# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "pandas",
#     "scipy",
#     "scikit-learn",
#     "joblib",
# ]
# ///
"""Inference engine for small-N congeneric series QSAR and ML models.

Loads serialized model bundles (.joblib) produced by congeneric_series_ml.py,
extracts matching molecular features using the authoritative feature computation
pipeline, applies the fitted preprocessor, predicts target endpoint values, and
enforces rigorous Applicability Domain (AD) safeguards using maximum Tanimoto
similarity to the training distribution.
"""

import argparse
from collections.abc import Sequence
import os
import sys
import types
from typing import Any, Final, Optional

import __main__
import joblib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

# Add script directory to sys.path so congeneric_series_ml can be imported.
_script_dir = os.path.dirname(os.path.abspath(__file__))
if _script_dir not in sys.path:
  sys.path.insert(0, _script_dir)

from congeneric_series_ml import (
    BoundedRegressor,
    DEFAULT_MORGAN_FP_BITS,
    LeakageFreePreprocessor,
    compute_molecular_features,
)

# Ensure unpickling compatibility for bundles pickled when congeneric_series_ml
# ran as __main__.
__main__.BoundedRegressor = BoundedRegressor
__main__.LeakageFreePreprocessor = LeakageFreePreprocessor

if "congeneric_series_ml" not in sys.modules:
  _c_mod = types.ModuleType("congeneric_series_ml")
  _c_mod.BoundedRegressor = BoundedRegressor
  _c_mod.LeakageFreePreprocessor = LeakageFreePreprocessor
  sys.modules["congeneric_series_ml"] = _c_mod
else:
  sys.modules["congeneric_series_ml"].BoundedRegressor = BoundedRegressor
  sys.modules["congeneric_series_ml"].LeakageFreePreprocessor = (
      LeakageFreePreprocessor
  )


try:
  from rdkit.Chem import rdFingerprintGenerator

  _mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

  def _mol_to_fp(m: Chem.Mol):
    return _mfpgen.GetFingerprint(m)

except (ImportError, AttributeError):

  def _mol_to_fp(m: Chem.Mol):
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)


_SMILES_COL_CANDIDATES: Final[tuple[str, ...]] = (
    "Clean_SMILES",
    "SMILES",
    "canonical_smiles",
    "smiles",
    "Smiles",
)
_ID_COL_CANDIDATES: Final[tuple[str, ...]] = (
    "Compound_ID",
    "compound_id",
    "ID",
    "id",
    "Name",
    "name",
    "Molecule_ID",
)
_P_SCALE_KEYWORDS: Final[tuple[str, ...]] = (
    "pki",
    "pic50",
    "pkd",
    "pec50",
    "pchembl",
)
# Physical thermodynamic & aqueous diffusion limits for p-scale affinities.
_P_SCALE_MIN: Final[float] = 0.0
_P_SCALE_MAX: Final[float] = 15.0


def compute_max_tanimoto(
    query_mol: Optional[Chem.Mol], train_fps: Sequence[Any]
) -> float:
  """Computes maximum Tanimoto similarity of query molecule to training set.

  Args:
    query_mol: Parsed RDKit molecule, or None if the candidate SMILES failed
      parsing.
    train_fps: Precomputed Morgan fingerprint bitvectors for the training set.

  Returns:
    Maximum Tanimoto similarity in [0.0, 1.0] for a valid `query_mol` (`0.0` if
    `train_fps` is empty or shares zero fingerprint bits), or `NaN` when
    `query_mol` is None.
  """
  if query_mol is None:
    return float("nan")
  if not train_fps:
    return 0.0
  q_fp = _mol_to_fp(query_mol)
  sims = DataStructs.BulkTanimotoSimilarity(q_fp, list(train_fps))
  return float(max(sims)) if sims else 0.0


def _parse_candidate_smiles(value: Any) -> Optional[Chem.Mol]:
  """Parses a candidate SMILES value, rejecting blank or invalid strings."""
  if pd.isna(value):
    return None
  text = str(value).strip()
  if not text:
    return None
  mol = Chem.MolFromSmiles(text)
  if mol is None or mol.GetNumAtoms() == 0:
    return None
  return mol


def _resolve_smiles_col(df: pd.DataFrame, smiles_col: str) -> str:
  """Resolves and validates the SMILES column name against `df.columns`."""
  if smiles_col in df.columns:
    return smiles_col
  if not smiles_col or smiles_col in _SMILES_COL_CANDIDATES:
    resolved = next(
        (c for c in _SMILES_COL_CANDIDATES if c in df.columns), None
    )
    if resolved is not None:
      return resolved
    raise ValueError(
        f"SMILES column {smiles_col!r} not found in DataFrame columns"
        f" {list(df.columns)} (also checked {list(_SMILES_COL_CANDIDATES)})."
    )
  raise ValueError(
      f"SMILES column {smiles_col!r} not found in DataFrame columns"
      f" {list(df.columns)}."
  )


def _resolve_id_col(df: pd.DataFrame, id_col: Optional[str]) -> Optional[str]:
  """Resolves an identifier column in `df.columns`, or None if absent."""
  if id_col:
    if id_col in df.columns:
      return id_col
    raise ValueError(
        f"ID column {id_col!r} not found in DataFrame columns"
        f" {list(df.columns)}."
    )
  return next((c for c in _ID_COL_CANDIDATES if c in df.columns), None)


def score_candidate_dataframe(
    df: pd.DataFrame,
    bundle: dict[str, Any],
    smiles_col: str = "SMILES",
    id_col: Optional[str] = None,
    min_tanimoto: float = 0.40,
) -> pd.DataFrame:
  """Scores candidate molecules with applicability-domain safeguards.

  Args:
    df: Candidate molecule table.
    bundle: Serialized QSAR model dictionary produced by
      `congeneric_series_ml.py`.
    smiles_col: Requested SMILES column name. When set to a standard alias (such
      as default `"SMILES"`) and absent from `df`, falls back across
      `Clean_SMILES`, `SMILES`, `canonical_smiles`, `smiles`, and `Smiles`.
    id_col: Optional compound identifier column name. When None, preserves the
      first matching standard identifier column (`Compound_ID`, `compound_id`,
      `ID`, `id`, `Name`, `name`, `Molecule_ID`) if present in `df`, or
      synthesizes `CMPD_0001...` identifiers otherwise.
    min_tanimoto: Minimum nearest-neighbor Morgan Tanimoto similarity required
      for `In_Applicability_Domain=True`.

  Returns:
    Scored DataFrame aligned 1:1 with `df` rows, retaining unparseable SMILES
    rows with `NaN` predictions and `NaN` `Max_Tanimoto_Similarity` under
    `Confidence_Tier="INVALID_SMILES"`.

  Raises:
    ValueError: If no valid SMILES column can be resolved in `df.columns`, or if
      an explicitly specified `id_col` is absent from `df.columns`.
  """
  df = df.reset_index(drop=True)
  actual_smiles_col = _resolve_smiles_col(df, smiles_col)
  resolved_id_col = _resolve_id_col(df, id_col)

  model = bundle["model"]
  prep = bundle.get("preprocessor")
  fset_name = bundle.get("feature_set", "2D Descriptors")
  target_name = bundle.get("target", "Target")
  required_cols = bundle.get("feature_cols", [])
  model_name = bundle.get("model_name", type(model).__name__)
  training_smiles = bundle.get("training_smiles", [])
  train_range = bundle.get("training_range")

  if df.empty:
    cols = [
        resolved_id_col or "Compound_ID",
        "SMILES",
        f"Predicted_{target_name}",
        "Max_Tanimoto_Similarity",
        "In_Applicability_Domain",
        "Confidence_Tier",
    ]
    if train_range:
      cols.extend(["Extrapolated_Beyond_Train", "Training_Range"])
    cols.extend(["Model_Architecture", "Feature_Set"])
    return pd.DataFrame(columns=cols)

  needs_3d = "3D" in fset_name or "Shape" in fset_name
  fp_bit_indices = [
      int(c[3:])
      for c in required_cols
      if c.startswith("FP_") and c[3:].isdigit()
  ]
  fp_bits = (
      max(fp_bit_indices) + 1 if fp_bit_indices else DEFAULT_MORGAN_FP_BITS
  )
  print(
      f"Extracting features for {len(df)} candidate molecules"
      f" (compute_3d={needs_3d})..."
  )
  feat_df, _, _, _, actual_smiles_col = compute_molecular_features(
      df,
      smiles_col=actual_smiles_col,
      compute_3d=needs_3d,
      fp_bits=fp_bits,
  )

  # Populate external columns (e.g. physics_cols) or missing features.
  for c in required_cols:
    if c not in feat_df.columns:
      feat_df[c] = df[c] if c in df.columns else 0.0

  x_raw = feat_df[required_cols].values.astype(np.float64)
  x_raw = np.nan_to_num(x_raw, nan=0.0, posinf=0.0, neginf=0.0)

  x_proc = prep.transform(x_raw) if prep is not None else x_raw
  preds = np.array(model.predict(x_proc), dtype=np.float64, copy=True)

  # Physical safeguards (only for p-scale affinity, diffusion limits).
  # Does NOT clamp to empirical training set maximum!
  is_p_scale = bundle.get("is_p_scale")
  if is_p_scale is None:
    t_low = target_name.lower()
    is_p_scale = any(p in t_low for p in _P_SCALE_KEYWORDS)

  mols = [_parse_candidate_smiles(s) for s in df[actual_smiles_col]]
  for i, m in enumerate(mols):
    if m is None:
      preds[i] = np.nan
    elif is_p_scale:
      # Enforce physical thermodynamic & aqueous diffusion limits [0.0, 15.0].
      preds[i] = float(np.clip(preds[i], _P_SCALE_MIN, _P_SCALE_MAX))
    # Non-p-scale endpoints remain completely unconstrained.

  # Precompute training set bitvectors for Applicability Domain check.
  train_bitvectors = [
      _mol_to_fp(m)
      for m in (_parse_candidate_smiles(s) for s in training_smiles)
      if m is not None
  ]

  max_sims = [compute_max_tanimoto(m, train_bitvectors) for m in mols]
  in_ad = [
      bool(
          m is not None
          and m.GetNumHeavyAtoms() > 2
          and np.isfinite(s)
          and s >= min_tanimoto
      )
      for m, s in zip(mols, max_sims, strict=True)
  ]
  confidence_tiers = [
      "INVALID_SMILES"
      if m is None
      else ("IN_DOMAIN" if ad else "OUT_OF_DOMAIN")
      for m, ad in zip(mols, in_ad, strict=True)
  ]

  out_df = pd.DataFrame()
  if resolved_id_col is not None:
    out_df[resolved_id_col] = df[resolved_id_col]
  else:
    out_df["Compound_ID"] = [f"CMPD_{i+1:04d}" for i in range(len(df))]

  out_df["SMILES"] = df[actual_smiles_col]
  out_df[f"Predicted_{target_name}"] = np.round(preds, 3)
  out_df["Max_Tanimoto_Similarity"] = np.round(max_sims, 3)
  out_df["In_Applicability_Domain"] = in_ad
  out_df["Confidence_Tier"] = confidence_tiers
  if train_range:
    tr_min, tr_max = float(train_range[0]), float(train_range[1])
    out_df["Extrapolated_Beyond_Train"] = np.isfinite(preds) & (
        (preds < tr_min) | (preds > tr_max)
    )
    out_df["Training_Range"] = f"[{tr_min:.2f}, {tr_max:.2f}]"
  out_df["Model_Architecture"] = model_name
  out_df["Feature_Set"] = fset_name
  return out_df


def main() -> None:
  """CLI entry point for scoring candidate SMILES with a QSAR model bundle."""
  parser = argparse.ArgumentParser(
      description=(
          "Predict potency / activity for candidate molecules using trained"
          " QSAR model bundles."
      )
  )
  parser.add_argument(
      "--model_path",
      "--model-path",
      "--model",
      "--bundle",
      dest="model_path",
      type=str,
      required=True,
      help="Path to serialized .joblib model bundle.",
  )
  parser.add_argument(
      "--input_csv",
      "--input-csv",
      "--input",
      dest="input_csv",
      type=str,
      required=True,
      help="Path to CSV containing candidate molecule SMILES.",
  )
  parser.add_argument(
      "--smiles_col",
      "--smiles-col",
      dest="smiles_col",
      type=str,
      default="SMILES",
      help="Column name for molecule SMILES (default: SMILES).",
  )
  parser.add_argument(
      "--id_col",
      "--id-col",
      dest="id_col",
      type=str,
      default=None,
      help="Column name for compound identifier.",
  )
  parser.add_argument(
      "--output_csv",
      "--output-csv",
      "--output",
      dest="output_csv",
      type=str,
      default="predicted_potency.csv",
      help="Path for saving scored predictions CSV.",
  )
  parser.add_argument(
      "--min_tanimoto",
      "--min-tanimoto",
      dest="min_tanimoto",
      type=float,
      default=0.40,
      help=(
          "Minimum maximum-Tanimoto similarity threshold for Applicability"
          " Domain (default: 0.40)."
      ),
  )
  args = parser.parse_args()

  if not os.path.exists(args.model_path):
    print(
        f"Error: Model bundle not found at: {args.model_path}", file=sys.stderr
    )
    sys.exit(1)

  if not os.path.exists(args.input_csv):
    print(f"Error: Input CSV not found at: {args.input_csv}", file=sys.stderr)
    sys.exit(1)

  print(f"Loading model bundle from: {args.model_path}")
  bundle = joblib.load(args.model_path)

  model = bundle["model"]
  fset_name = bundle.get("feature_set", "2D Descriptors")
  target_name = bundle.get("target", "Target")
  required_cols = bundle.get("feature_cols", [])
  model_name = bundle.get("model_name", type(model).__name__)

  print(
      f"Loaded {model_name} model for target '{target_name}' using feature set"
      f" '{fset_name}'."
  )
  print(f"Number of required model features: {len(required_cols)}")

  df = pd.read_csv(args.input_csv)
  try:
    out_df = score_candidate_dataframe(
        df,
        bundle,
        smiles_col=args.smiles_col,
        id_col=args.id_col,
        min_tanimoto=args.min_tanimoto,
    )
  except ValueError as exc:
    print(f"Error: {exc}", file=sys.stderr)
    sys.exit(1)

  os.makedirs(
      os.path.dirname(os.path.abspath(args.output_csv)) or ".", exist_ok=True
  )
  out_df.to_csv(args.output_csv, index=False)
  print(
      f"\nScored {len(out_df)} compounds. Results saved to: {args.output_csv}"
  )

  print("\n--- Predicted Potency Summary ---")
  print(out_df.to_string(index=False))

  id_col_name = out_df.columns[0]
  tiers = out_df["Confidence_Tier"].tolist()
  invalid_indices = [
      i for i, tier in enumerate(tiers) if tier == "INVALID_SMILES"
  ]
  if invalid_indices:
    print(
        f"\n[!] WARNING: {len(invalid_indices)} compound(s) have invalid or"
        " unparseable SMILES (predictions set to NaN):"
    )
    for idx in invalid_indices:
      cid = out_df.iloc[idx][id_col_name]
      raw_smi = out_df.iloc[idx]["SMILES"]
      print(f"    - {cid}: SMILES = {raw_smi!r}")

  out_of_domain = [i for i, tier in enumerate(tiers) if tier == "OUT_OF_DOMAIN"]
  if out_of_domain:
    print(
        f"\n[!] WARNING: {len(out_of_domain)} compound(s) fall outside the"
        f" Applicability Domain (Max Tanimoto < {args.min_tanimoto:.2f}):"
    )
    for idx in out_of_domain:
      cid = out_df.iloc[idx][id_col_name]
      sim = out_df.iloc[idx]["Max_Tanimoto_Similarity"]
      print(f"    - {cid}: Max Tanimoto = {sim:.3f}")


if __name__ == "__main__":
  main()
