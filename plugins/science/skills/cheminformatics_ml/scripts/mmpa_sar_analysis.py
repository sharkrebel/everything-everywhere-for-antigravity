# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "pandas",
#     "scipy",
#     "scikit-learn",
#     "matplotlib",
#     "seaborn>=0.12.0",
# ]
# ///

"""Matched Molecular Pair Analysis (MMPA), Topological Pair Matching, and SALI Cliffs.

Performs:
1. Hussain-Rea matched molecular pair (MMP) fragmentation (1-cut R-group
   substitutions and 2-cut internal linker/core replacements) with medicinal
   chemistry invariants:
   - Minimum conserved scaffold size (>= 12 heavy atoms, >= 50% parent fraction)
   - Maximum variable fragment size (<= 10 heavy atoms, <= 35% parent fraction)
   - Rejection of trivial 1-2 atom cores (N[*:1], C[*:1], Cl[*:1], F[*:1])
2. General topological & ring-transformation pair matching:
   - Because Hussain-Rea (rdMMPA.FragmentMol) only cuts acyclic bonds, any
     transformation that forms or breaks a ring (Delta_Rings != 0, e.g. ring
     closure/cyclization of any ring size, ring expansion/contraction, or
     scaffold hop) or compares compounds across user-defined cohorts
     (--group_col) is captured via topological path similarity
     (RDKit Topological Fingerprints).
   - Computes Delta_Activity, Fold_Change, Delta_RotB, Delta_Rings,
     Delta_MaxRingSize, Delta_CNS_MPO, Delta_cLogP, Delta_TPSA, and Delta_MW.
3. Structure-Activity Landscape Index (SALI) activity cliff detection:
   SALI(i, j) = |Act_i - Act_j| / (1 - Tanimoto(i, j))
"""

import argparse
import collections
import os
import re
import sys
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator, rdMMPA, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

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

from compat import get_morgan_fingerprint
from skill_paths import add_skill_scripts

add_skill_scripts("cheminformatics_core")
add_skill_scripts("cheminformatics_ml")

from descriptors import calculate_cns_mpo
import plot_archetypes

_calculate_cns_mpo = calculate_cns_mpo


def _max_ring_size(mol: Chem.Mol) -> int:
  """Returns the size of the largest SSSR ring in `mol` (0 if acyclic)."""
  rings = mol.GetRingInfo().AtomRings()
  return max((len(r) for r in rings), default=0)


def fragment_molecule_mmpa(
    mol_or_smiles: Union[str, Chem.Mol],
    max_cuts: int = 2,
    max_cut_bonds: int = 20,
    min_core_heavies: int = 12,
    min_core_ratio: float = 0.50,
    max_substituent_heavies: int = 10,
    max_substituent_ratio: float = 0.35,
) -> List[Tuple[str, str]]:
  """Generates Hussain-Rea (conserved_scaffold, variable_fragment) pairs.

  Supports both:
  1. Single-cut and multi-cut peripheral R-group substitutions (where `core` is
     the conserved scaffold and `chains` is the variable R-group).
  2. Two-cut internal linker/core replacements (where `chains` holds the two
     conserved flanking pharmacophores >= min_core_heavies and `core` is the
     variable internal linker/ring <= max_substituent_heavies).

  Strictly prohibits 1-2 atom trivial cores (e.g. N[*:1], C[*:1], Cl[*:1]) from
  being indexed as scaffolds.
  """
  if isinstance(mol_or_smiles, Chem.Mol):
    mol = mol_or_smiles
  elif isinstance(mol_or_smiles, str):
    mol = Chem.MolFromSmiles(mol_or_smiles)
    if mol is None:
      return []
  else:
    return []

  parent_ha = mol.GetNumHeavyAtoms()
  if parent_ha < (min_core_heavies + 1):
    return []

  try:
    frags = rdMMPA.FragmentMol(mol, maxCuts=max_cuts, maxCutBonds=max_cut_bonds)
  except Exception:
    return []

  results = []
  seen = set()
  for core, chains in frags:
    if core is not None and chains is not None:
      c_ha = core.GetNumHeavyAtoms()
      ch_ha = chains.GetNumHeavyAtoms()
      # Orientation A: `core` is the large conserved scaffold, `chains` is the small substituent
      if (
          c_ha >= min_core_heavies
          and (c_ha / parent_ha) >= min_core_ratio
          and ch_ha <= max_substituent_heavies
          and (ch_ha / parent_ha) <= max_substituent_ratio
          and c_ha > ch_ha
      ):
        pair = (Chem.MolToSmiles(core), Chem.MolToSmiles(chains))
        if pair not in seen:
          seen.add(pair)
          results.append(pair)
      # Orientation B (2-cut linker/ring replacement): `chains` holds the large conserved flanks,
      # and `core` is the small variable internal linker or ring
      elif (
          ch_ha >= min_core_heavies
          and (ch_ha / parent_ha) >= min_core_ratio
          and c_ha <= max_substituent_heavies
          and (c_ha / parent_ha) <= max_substituent_ratio
          and ch_ha > c_ha
      ):
        pair = (Chem.MolToSmiles(chains), Chem.MolToSmiles(core))
        if pair not in seen:
          seen.add(pair)
          results.append(pair)
    elif core is None and chains is not None:
      ch_smi = Chem.MolToSmiles(chains)
      parts = ch_smi.split(".")
      if len(parts) == 2:
        m0 = Chem.MolFromSmiles(parts[0])
        m1 = Chem.MolFromSmiles(parts[1])
        if not m0 or not m1:
          continue
        ha0 = m0.GetNumHeavyAtoms()
        ha1 = m1.GetNumHeavyAtoms()

        if ha0 >= ha1:
          core_smi, sub_smi = parts[0], parts[1]
          c_ha, s_ha = ha0, ha1
        else:
          core_smi, sub_smi = parts[1], parts[0]
          c_ha, s_ha = ha1, ha0

        if (
            c_ha >= min_core_heavies
            and (c_ha / parent_ha) >= min_core_ratio
            and s_ha <= max_substituent_heavies
            and (s_ha / parent_ha) <= max_substituent_ratio
            and c_ha > s_ha
        ):
          pair = (core_smi, sub_smi)
          if pair not in seen:
            seen.add(pair)
            results.append(pair)

  return results


def find_matched_molecular_pairs(
    df: pd.DataFrame,
    smiles_col: str = "SMILES",
    activity_col: Optional[str] = None,
    id_col: Optional[str] = None,
    max_cuts: int = 2,
    min_core_heavies: int = 12,
    min_core_ratio: float = 0.50,
    max_substituent_heavies: int = 10,
    max_substituent_ratio: float = 0.35,
) -> pd.DataFrame:
  """Finds all matched molecular pairs in a dataframe that share an identical conserved scaffold.

  Computes activity delta (Act_2 - Act_1), fold change, and property shifts.
  """
  core_index = collections.defaultdict(list)

  for idx, row in df.iterrows():
    smi = str(row[smiles_col]) if pd.notna(row.get(smiles_col)) else ""
    if not smi:
      continue
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
      continue
    cid = (
        str(row[id_col])
        if (id_col and id_col in df.columns and pd.notna(row.get(id_col)))
        else f"CMPD_{idx+1:04d}"
    )
    act = (
        float(row[activity_col])
        if (
            activity_col
            and activity_col in df.columns
            and pd.notna(row.get(activity_col))
        )
        else None
    )

    rotb = float(Descriptors.NumRotatableBonds(mol))
    mw = float(Descriptors.MolWt(mol))
    clogp = float(Descriptors.MolLogP(mol))
    tpsa = float(Descriptors.TPSA(mol))
    n_rings = int(rdMolDescriptors.CalcNumRings(mol))
    max_ring = _max_ring_size(mol)
    cns_mpo = _calculate_cns_mpo(mol)

    frags = fragment_molecule_mmpa(
        mol,
        max_cuts=max_cuts,
        min_core_heavies=min_core_heavies,
        min_core_ratio=min_core_ratio,
        max_substituent_heavies=max_substituent_heavies,
        max_substituent_ratio=max_substituent_ratio,
    )
    for core, substituent in frags:
      core_index[core].append({
          "df_idx": idx,
          "compound_id": cid,
          "smiles": smi,
          "substituent": substituent,
          "activity": act,
          "rotb": rotb,
          "mw": mw,
          "clogp": clogp,
          "tpsa": tpsa,
          "rings": n_rings,
          "max_ring": max_ring,
          "cns_mpo": cns_mpo,
      })

  best_by_compound_pair: Dict[
      Tuple[Any, Any], Tuple[Tuple[int, int, int], Dict[str, Any]]
  ] = {}

  for core, entries in core_index.items():
    if len(entries) < 2:
      continue
    core_mol = Chem.MolFromSmiles(core)
    core_ha = core_mol.GetNumHeavyAtoms() if core_mol is not None else 0
    for i in range(len(entries)):
      for j in range(i + 1, len(entries)):
        e1, e2 = entries[i], entries[j]
        if e1["df_idx"] == e2["df_idx"]:
          continue
        r1, r2 = e1["substituent"], e2["substituent"]
        if r1 == r2:
          continue
        # Reject spectator 2-cuts where one of the dot-separated fragments did not change
        if "." in r1 or "." in r2:
          s1 = {re.sub(r"\[\*:\d+\]", "[*]", frag) for frag in r1.split(".")}
          s2 = {re.sub(r"\[\*:\d+\]", "[*]", frag) for frag in r2.split(".")}
          if s1 & s2:
            continue

        m_r1 = Chem.MolFromSmiles(r1)
        m_r2 = Chem.MolFromSmiles(r2)
        ha1 = m_r1.GetNumHeavyAtoms() if m_r1 is not None else 99
        ha2 = m_r2.GetNumHeavyAtoms() if m_r2 is not None else 99
        # Canonicalize transformation direction so (R1, R2) is consistently ordered across pairs
        if (ha1, r1) > (ha2, r2):
          e1, e2 = e2, e1
          r1, r2 = r2, r1
          ha1, ha2 = ha2, ha1

        n_cuts = (
            2
            if ("." in r1 or "." in r2 or "[*:2]" in r1 or "[*:2]" in r2)
            else 1
        )
        cut_score = (n_cuts, ha1 + ha2, -core_ha)
        comp_pair_key = tuple(sorted([e1["df_idx"], e2["df_idx"]]))
        if (
            comp_pair_key in best_by_compound_pair
            and cut_score >= best_by_compound_pair[comp_pair_key][0]
        ):
          continue

        p_data = {
            "Compound_1": e1["compound_id"],
            "Compound_2": e2["compound_id"],
            "SMILES_1": e1["smiles"],
            "SMILES_2": e2["smiles"],
            "Core": core,
            "R1": r1,
            "R2": r2,
            "Transformation": f"{r1} >> {r2}",
            "Delta_RotB": round(e2["rotb"] - e1["rotb"], 1),
            "Delta_Rings": int(e2["rings"] - e1["rings"]),
            "Delta_MaxRingSize": int(e2["max_ring"] - e1["max_ring"]),
            "Delta_CNS_MPO": round(e2["cns_mpo"] - e1["cns_mpo"], 2),
            "Delta_MW": round(e2["mw"] - e1["mw"], 2),
            "Delta_cLogP": round(e2["clogp"] - e1["clogp"], 2),
            "Delta_TPSA": round(e2["tpsa"] - e1["tpsa"], 2),
        }

        if e1["activity"] is not None and e2["activity"] is not None:
          delta_act = e2["activity"] - e1["activity"]
          p_data["Activity_1"] = e1["activity"]
          p_data["Activity_2"] = e2["activity"]
          p_data["Delta_Activity"] = round(delta_act, 3)
          p_data["Abs_Delta_Activity"] = round(abs(delta_act), 3)
          p_data["Fold_Change"] = round(float(10.0**delta_act), 2)

        best_by_compound_pair[comp_pair_key] = (cut_score, p_data)

  pairs = [val[1] for val in best_by_compound_pair.values()]
  pair_df = pd.DataFrame(pairs)
  if not pair_df.empty and "Abs_Delta_Activity" in pair_df.columns:
    pair_df = pair_df.sort_values(by="Abs_Delta_Activity", ascending=False)
  return pair_df


def find_topological_matched_pairs(
    df: pd.DataFrame,
    smiles_col: str = "SMILES",
    activity_col: Optional[str] = None,
    id_col: Optional[str] = None,
    group_col: Optional[str] = None,
    core_smarts: Optional[str] = None,
    min_similarity: float = 0.60,
    max_ha_difference: int = 16,
    max_pairs_per_compound: Optional[int] = None,
) -> pd.DataFrame:
  """Identifies general topological matched pairs across ring transformations, scaffold hops, or cohorts.

  Because Hussain-Rea (`rdMMPA.FragmentMol`) only cuts acyclic bonds, any
  transformation that forms or opens a ring (`Delta_Rings != 0` or
  `Delta_MaxRingSize != 0`, e.g., cyclization of any ring size, ring
  expansion/contraction, or fused-ring formation) cannot be found by cutting
  acyclic bonds alone. This function uses RDKit Topological Path Fingerprints
  to match structurally conserved pairs across:
  - Ring closures / cyclizations (`Rings_2 > Rings_1`) of any ring size
  - Ring expansions / contractions (`MaxRingSize_2 != MaxRingSize_1`)
  - User-specified cohort columns (`--group-col`)
  """
  patt = Chem.MolFromSmarts(core_smarts) if core_smarts else None
  rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(
      minPath=1, maxPath=7, fpSize=2048
  )

  def _extract_records(require_pattern: bool) -> List[Dict[str, Any]]:
    recs = []
    for idx, row in df.iterrows():
      smi = str(row[smiles_col]) if pd.notna(row.get(smiles_col)) else ""
      if not smi or smi == "nan":
        continue
      mol = Chem.MolFromSmiles(smi)
      if mol is None:
        continue
      if (
          require_pattern
          and patt is not None
          and not mol.HasSubstructMatch(patt)
      ):
        continue

      cid = (
          str(row[id_col])
          if (id_col and id_col in df.columns and pd.notna(row.get(id_col)))
          else f"CMPD_{idx+1:04d}"
      )
      act = (
          float(row[activity_col])
          if (
              activity_col
              and activity_col in df.columns
              and pd.notna(row.get(activity_col))
          )
          else None
      )
      grp_val = (
          str(row[group_col]).strip()
          if (
              group_col
              and group_col in df.columns
              and pd.notna(row.get(group_col))
          )
          else None
      )
      n_rings = int(rdMolDescriptors.CalcNumRings(mol))
      max_ring = _max_ring_size(mol)
      try:
        murcko = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
      except Exception:
        murcko = ""

      recs.append({
          "df_idx": idx,
          "compound_id": cid,
          "smiles": smi,
          "mol": mol,
          "fp": rdkit_gen.GetFingerprint(mol),
          "ha": mol.GetNumHeavyAtoms(),
          "activity": act,
          "group": grp_val,
          "rings": n_rings,
          "max_ring": max_ring,
          "murcko": murcko,
          "rotb": float(Descriptors.NumRotatableBonds(mol)),
          "cns_mpo": _calculate_cns_mpo(mol),
          "clogp": float(Descriptors.MolLogP(mol)),
          "tpsa": float(Descriptors.TPSA(mol)),
          "mw": float(Descriptors.MolWt(mol)),
      })
    return recs

  records = _extract_records(require_pattern=True)
  if len(records) < 2 and patt is not None:
    records = _extract_records(require_pattern=False)
  if len(records) < 2:
    return pd.DataFrame()

  eff_threshold = float(min_similarity)
  eff_ha_diff = min(max_ha_difference, 8)

  # Index 1-cut acyclic R-group pairs so peripheral ring substituents are not
  # misclassified as topological ring closures when group_col is omitted.
  one_cut_pairs: set[Tuple[Any, Any]] = set()
  if not group_col:
    one_cut_cores = collections.defaultdict(list)
    for r in records:
      for core_smi, sub_smi in fragment_molecule_mmpa(r["mol"], max_cuts=1):
        one_cut_cores[core_smi].append((r["df_idx"], sub_smi))
    for entries in one_cut_cores.values():
      if len(entries) < 2:
        continue
      for idx_a in range(len(entries)):
        for idx_b in range(idx_a + 1, len(entries)):
          df_a, sub_a = entries[idx_a]
          df_b, sub_b = entries[idx_b]
          if df_a != df_b and sub_a != sub_b:
            one_cut_pairs.add((min(df_a, df_b), max(df_a, df_b)))

  pairs: List[Dict[str, Any]] = []
  n_rec = len(records)
  for i in range(n_rec):
    for j in range(i + 1, n_rec):
      r_a, r_b = records[i], records[j]
      if abs(r_a["ha"] - r_b["ha"]) > eff_ha_diff:
        continue

      # Determine whether (r_a, r_b) represents a topological / cohort transformation
      if group_col and r_a["group"] is not None and r_b["group"] is not None:
        if r_a["group"] == r_b["group"]:
          continue
      else:
        # General ring-topology transformation: different ring count or different max ring size
        if r_a["rings"] == r_b["rings"] and r_a["max_ring"] == r_b["max_ring"]:
          continue
        # Exclude 1-cut peripheral R-group ring additions (e.g. N-Me -> N-cPr)
        # that leave the largest ring size unchanged; those belong in mmpa_df.
        pair_key = (
            min(r_a["df_idx"], r_b["df_idx"]),
            max(r_a["df_idx"], r_b["df_idx"]),
        )
        if r_a["max_ring"] == r_b["max_ring"] and pair_key in one_cut_pairs:
          continue

      sim = DataStructs.TanimotoSimilarity(r_a["fp"], r_b["fp"])
      if sim < eff_threshold:
        continue

      # Order pair canonically so Compound_1 is the lower-ring-count / smaller-ring / smaller-MW precursor
      # and Compound_2 is the cyclized / ring-expanded derivative
      key_a = (r_a["max_ring"], r_a["rings"], r_a["ha"], r_a["compound_id"])
      key_b = (r_b["max_ring"], r_b["rings"], r_b["ha"], r_b["compound_id"])
      if key_a <= key_b:
        r1, r2 = r_a, r_b
      else:
        r1, r2 = r_b, r_a

      d_rings = r2["rings"] - r1["rings"]
      d_max_ring = r2["max_ring"] - r1["max_ring"]
      if d_rings > 0 and d_max_ring != 0:
        trans_type = (
            f"Ring Closure (+{d_rings} ring, max ring"
            f" {r1['max_ring']}->{r2['max_ring']})"
        )
      elif d_max_ring != 0:
        trans_type = f"Ring Size Shift ({r1['max_ring']}->{r2['max_ring']})"
      elif d_rings > 0:
        trans_type = (
            f"Internal Ring Addition (+{d_rings} ring, max ring"
            f" {r2['max_ring']})"
        )
      elif r1["murcko"] != r2["murcko"]:
        trans_type = "Scaffold Hop"
      else:
        trans_type = f"Cohort Shift ({r1['group']}->{r2['group']})"

      pair_entry: Dict[str, Any] = {
          "Compound_1": r1["compound_id"],
          "Compound_2": r2["compound_id"],
          "SMILES_1": r1["smiles"],
          "SMILES_2": r2["smiles"],
          "Transformation_Type": trans_type,
          "Topological_Similarity": round(sim, 3),
          "Rings_1": r1["rings"],
          "Rings_2": r2["rings"],
          "Delta_Rings": int(d_rings),
          "MaxRingSize_1": r1["max_ring"],
          "MaxRingSize_2": r2["max_ring"],
          "Delta_MaxRingSize": int(d_max_ring),
          "RotB_1": r1["rotb"],
          "RotB_2": r2["rotb"],
          "Delta_RotB": round(r2["rotb"] - r1["rotb"], 1),
          "CNS_MPO_1": r1["cns_mpo"],
          "CNS_MPO_2": r2["cns_mpo"],
          "Delta_CNS_MPO": round(r2["cns_mpo"] - r1["cns_mpo"], 2),
          "Delta_cLogP": round(r2["clogp"] - r1["clogp"], 2),
          "Delta_TPSA": round(r2["tpsa"] - r1["tpsa"], 2),
          "Delta_MW": round(r2["mw"] - r1["mw"], 2),
      }

      has_both_act = r1["activity"] is not None and r2["activity"] is not None
      if has_both_act:
        d_act = r2["activity"] - r1["activity"]
        pair_entry["Activity_1"] = r1["activity"]
        pair_entry["Activity_2"] = r2["activity"]
        pair_entry["Delta_Activity"] = round(d_act, 3)
        pair_entry["Delta_pActivity"] = round(d_act, 3)
        pair_entry["Fold_Change"] = round(float(10.0**d_act), 1)

      pairs.append(pair_entry)

  df_pairs = pd.DataFrame(pairs)
  if not df_pairs.empty:
    if "Delta_Activity" in df_pairs.columns:
      assayed_df = df_pairs[df_pairs["Delta_Activity"].notna()].sort_values(
          by=["Topological_Similarity", "Delta_Activity"],
          ascending=[False, False],
      )
      unassayed_df = df_pairs[df_pairs["Delta_Activity"].isna()].sort_values(
          by="Topological_Similarity", ascending=False
      )
      df_pairs = pd.concat([assayed_df, unassayed_df], ignore_index=True)
    else:
      df_pairs = df_pairs.sort_values(
          by="Topological_Similarity", ascending=False
      ).reset_index(drop=True)
  return df_pairs


def calculate_sali_matrix(
    smiles_list: List[str],
    activities: List[float],
    ids: Optional[List[str]] = None,
    min_similarity: float = 0.60,
) -> List[Dict[str, Any]]:
  """Calculates Structure-Activity Landscape Index (SALI) for activity cliffs:

  SALI = |Act_i - Act_j| / (1 - Tanimoto(i, j))
  """
  mols = []
  valid_map = {}
  valid_acts = []
  valid_ids = []

  for idx, s in enumerate(smiles_list):
    if idx >= len(activities) or pd.isna(activities[idx]):
      continue
    try:
      act_val = float(activities[idx])
    except (TypeError, ValueError):
      continue
    if not np.isfinite(act_val):
      continue
    m = Chem.MolFromSmiles(str(s)) if pd.notna(s) else None
    if m is not None:
      valid_map[len(mols)] = idx
      mols.append(m)
      valid_acts.append(act_val)
      cid = (
          str(ids[idx])
          if (ids and idx < len(ids) and pd.notna(ids[idx]))
          else f"CMPD_{idx+1:04d}"
      )
      valid_ids.append(cid)

  if len(mols) < 2:
    return []

  fps = [
      get_morgan_fingerprint(m, radius=2, n_bits=2048, use_chirality=True)
      for m in mols
  ]

  n = len(mols)
  cliffs = []
  for i in range(n):
    for j in range(i + 1, n):
      sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
      if min_similarity <= sim < 0.9999:
        delta_act = abs(valid_acts[i] - valid_acts[j])
        sali = delta_act / max(1e-4, (1.0 - sim))
        cliffs.append({
            "Compound_1": valid_ids[i],
            "Compound_2": valid_ids[j],
            "idx_i": valid_map[i],
            "idx_j": valid_map[j],
            "smiles_i": smiles_list[valid_map[i]],
            "smiles_j": smiles_list[valid_map[j]],
            "Activity_1": valid_acts[i],
            "Activity_2": valid_acts[j],
            "delta_activity": round(delta_act, 2),
            "Delta_Activity": round(delta_act, 2),
            "similarity": round(sim, 3),
            "Similarity": round(sim, 3),
            "sali": round(sali, 2),
            "SALI": round(sali, 2),
        })

  cliffs.sort(key=lambda x: x["SALI"], reverse=True)
  return cliffs


def calculate_sali_dataframe(
    df: pd.DataFrame,
    smiles_col: str = "SMILES",
    activity_col: str = "pIC50",
    id_col: Optional[str] = None,
    min_similarity: float = 0.60,
) -> pd.DataFrame:
  """Computes SALI activity cliffs for a pandas DataFrame."""
  smiles = df[smiles_col].tolist()
  activities = df[activity_col].tolist()
  ids = df[id_col].tolist() if (id_col and id_col in df.columns) else None
  cliffs = calculate_sali_matrix(
      smiles, activities, ids=ids, min_similarity=min_similarity
  )
  return pd.DataFrame(cliffs)


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Matched Molecular Pair Analysis (MMPA), Topological Pair Matching,"
          " and SALI Activity Cliff Detector."
      )
  )
  parser.add_argument(
      "--input_csv",
      "--input-csv",
      "--input",
      dest="input_csv",
      type=str,
      help="Input CSV file containing SMILES and activity columns.",
  )
  parser.add_argument(
      "--smiles_col",
      "--smiles-col",
      dest="smiles_col",
      type=str,
      default="SMILES",
      help="Column name for SMILES strings.",
  )
  parser.add_argument(
      "--activity_col",
      "--activity-col",
      dest="activity_col",
      type=str,
      default=None,
      help="Column name for activity values (e.g. pIC50, pchembl_value).",
  )
  parser.add_argument(
      "--id_col",
      "--id-col",
      dest="id_col",
      type=str,
      default=None,
      help="Column name for compound identifiers.",
  )
  parser.add_argument(
      "--group_col",
      "--group-col",
      dest="group_col",
      type=str,
      default=None,
      help=(
          "Optional grouping/cohort column for cross-cohort topological pair"
          " matching."
      ),
  )
  parser.add_argument(
      "--smiles",
      nargs="+",
      help="Input SMILES list (for direct command-line execution).",
  )
  parser.add_argument(
      "--activities",
      nargs="+",
      type=float,
      help="Corresponding pIC50 or activity values.",
  )
  parser.add_argument(
      "--min_similarity",
      "--min-sim",
      dest="min_similarity",
      type=float,
      default=0.60,
      help=(
          "Minimum Tanimoto similarity for SALI cliffs and topological pairs"
          " (default: 0.60)."
      ),
  )
  parser.add_argument(
      "--max_cuts",
      dest="max_cuts",
      type=int,
      default=2,
      help="Max cuts for Hussain-Rea MMPA fragmentation (default: 2).",
  )
  parser.add_argument(
      "--min_core_heavies",
      "--min-core-heavies",
      dest="min_core_heavies",
      type=int,
      default=12,
      help=(
          "Minimum heavy atoms required for a valid core scaffold (default:"
          " 12)."
      ),
  )
  parser.add_argument(
      "--max_substituent_heavies",
      "--max-substituent-heavies",
      dest="max_substituent_heavies",
      type=int,
      default=10,
      help=(
          "Maximum heavy atoms allowed for a variable R-group or linker"
          " (default: 10)."
      ),
  )
  parser.add_argument(
      "--match_topological_pairs",
      "--match-topological-pairs",
      dest="match_topological_pairs",
      action="store_true",
      help="Run topological / ring-transformation matched pair analysis.",
  )
  parser.add_argument(
      "--core_smarts",
      "--core-smarts",
      dest="core_smarts",
      type=str,
      default=None,
      help=(
          "Optional core SMARTS pattern to restrict topological pair matching."
      ),
  )
  parser.add_argument(
      "--output_mmpa_csv",
      "--output-mmpa-csv",
      "--output-mmpa",
      dest="output_mmpa_csv",
      type=str,
      default=None,
      help="Output CSV path for matched molecular pairs.",
  )
  parser.add_argument(
      "--output_sali_csv",
      "--output-sali-csv",
      "--output-sali",
      dest="output_sali_csv",
      type=str,
      default=None,
      help="Output CSV path for SALI activity cliffs.",
  )
  parser.add_argument(
      "--output_paired_csv",
      "--output-paired-csv",
      "--output_topological_csv",
      "--output-topological-csv",
      dest="output_paired_csv",
      type=str,
      default=None,
      help=(
          "Output CSV path for topological / ring-transformation matched pairs."
      ),
  )
  parser.add_argument(
      "--output_csv",
      "--output-csv",
      "--output",
      "-o",
      dest="output_csv",
      type=str,
      default=None,
      help="Output CSV path for primary results (SALI cliffs).",
  )
  parser.add_argument(
      "--plot",
      "--output_png",
      "--output-png",
      "--plot_png",
      "--plot-png",
      "--output_plot_png",
      "--output-plot-png",
      dest="plot_png",
      type=str,
      default=None,
      help=(
          "Output PNG path to render publication-ready MMPA transformation"
          " shifts and SALI landscape plot."
      ),
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
      "--top_n",
      "--top-n",
      dest="top_n",
      type=int,
      default=20,
      help="Number of top cliffs/pairs to display (default: 20).",
  )
  parser.add_argument(
      "--max_pairs_per_compound",
      "--max-pairs-per-compound",
      dest="max_pairs_per_compound",
      type=int,
      default=None,
      help="Optional cap on topological matched pairs retained per compound.",
  )
  parser.add_argument(
      "--paired_x_col",
      "--paired-x-col",
      "--x_col",
      "--x-col",
      dest="paired_x_col",
      type=str,
      default=None,
      help=(
          "Property shift column for Panel C X-axis (e.g. Delta_CNS_MPO,"
          " Delta_cLogP, Delta_TPSA, Delta_MW, Delta_RotB)."
      ),
  )
  parser.add_argument(
      "--paired_y_col",
      "--paired-y-col",
      "--y_col",
      "--y-col",
      dest="paired_y_col",
      type=str,
      default=None,
      help=(
          "Activity/potency shift column for Panel C Y-axis and Panel D"
          " ranking (default: Delta_pActivity / Delta_Activity)."
      ),
  )
  args = parser.parse_args()

  def _df_to_gfm(sub_df: pd.DataFrame) -> str:
    if sub_df is None or sub_df.empty:
      return "_No pairs found._"
    cols = list(sub_df.columns)
    lines = [
        "| " + " | ".join(str(c) for c in cols) + " |",
        "| " + " | ".join(":---" for _ in cols) + " |",
    ]
    for _, r in sub_df.iterrows():
      row_cells = []
      for c in cols:
        v = r[c]
        if isinstance(v, float) and pd.notna(v):
          row_cells.append(f"{float(v):.3f}")
        elif pd.isna(v):
          row_cells.append("—")
        else:
          row_cells.append(str(v).replace("|", "/"))
      lines.append("| " + " | ".join(row_cells) + " |")
    return "\n".join(lines)

  md_sections = [
      "# Matched Molecular Pair (MMPA) & SALI Activity Cliff Analysis\n"
  ]

  if args.input_csv:
    df = pd.read_csv(args.input_csv)
    smiles_col = args.smiles_col
    if smiles_col not in df.columns:
      for c in [
          "smiles",
          "SMILES",
          "canonical_smiles",
          "Structure",
          "Clean_SMILES",
      ]:
        if c in df.columns:
          smiles_col = c
          break

    act_col = args.activity_col
    if act_col and act_col not in df.columns:
      # Case-insensitive match first
      ci_map = {c.lower(): c for c in df.columns}
      if act_col.lower() in ci_map:
        act_col = ci_map[act_col.lower()]
      else:
        act_col = None
    if not act_col or act_col not in df.columns:
      for c in [
          "pActivity_primary",
          "pActivity",
          "pActivity_composite",
          "pchembl_value",
          "pIC50",
          "pic50",
          "pKi",
          "pki",
          "potency",
          "Activity",
          "Standard_Value",
      ]:
        if c in df.columns and df[c].notna().any():
          act_col = c
          break
      if not act_col:
        for c in df.columns:
          cl = c.lower()
          if (
              any(k in cl for k in ("pactivity", "pic50", "pki", "pchembl"))
              and pd.api.types.is_numeric_dtype(df[c])
              and df[c].notna().any()
          ):
            act_col = c
            break

    def _resolve_short_compound_ids(
        frame: pd.DataFrame, primary_col: Optional[str]
    ) -> pd.Series:
      """Resolves concise human-readable compound IDs (<=22 chars) per row."""
      raw_Labels: list[str] = []
      for idx_r, row in frame.iterrows():
        # 1. Check primary column first if it is already concise and non-IUPAC
        prim_val = (
            str(row.get(primary_col, "")).strip()
            if primary_col
            and primary_col in frame.columns
            and pd.notna(row.get(primary_col))
            else ""
        )
        if prim_val.lower() in ("nan", "none", "null"):
          prim_val = ""
        # Extract trivial drug/lead name before parenthetical code, e.g. '<LeadName> (<CodeID>)' -> '<LeadName>'
        m_triv = re.match(
            r"^([A-Za-z][A-Za-z0-9_.\s-]{2,16})\s*\(([^)]+)\)$", prim_val
        )
        if m_triv and not re.search(
            r"\b(rac|rel|cis|trans)\b", m_triv.group(1), re.I
        ):
          raw_Labels.append(m_triv.group(1).strip())
          continue
        is_iupac_like = (
            len(prim_val) > 22
            or prim_val.count("(") >= 2
            or prim_val.count("-") >= 3
            or bool(
                re.search(
                    r"yl\)|methanone|carboxamide|sulfonamide|phenyl|pyridin|pyrimidin|quinazolin|pyrazol|indol|oxy",
                    prim_val,
                    re.I,
                )
            )
        )
        if prim_val and not is_iupac_like:
          raw_Labels.append(prim_val)
          continue

        # 2. Otherwise search short identifier columns on the same row
        chosen = ""
        paper_id = (
            str(row.get("paper_compound_id", "")).strip()
            if "paper_compound_id" in frame.columns
            and pd.notna(row.get("paper_compound_id"))
            else ""
        )
        if (
            paper_id
            and paper_id.lower() not in ("nan", "none")
            and len(paper_id) <= 16
        ):
          is_code = bool(
              re.match(
                  r"^(?:cmpd|compound|example|ex\b|chembl|[a-z]{2,}[-_]?\d)",
                  paper_id,
                  re.IGNORECASE,
              )
          )
          chosen = paper_id if is_code else f"Cmpd {paper_id}"
        if not chosen:
          for cand_c in (
              "patent_example_id",
              "compound_id",
              "Compound_ID",
              "ID",
              "molecule_chembl_id",
              "chembl_id",
          ):
            if cand_c in frame.columns and pd.notna(row.get(cand_c)):
              cv = str(row.get(cand_c, "")).strip()
              if cv and cv.lower() not in ("nan", "none") and len(cv) <= 22:
                chosen = cv
                break
        if (
            not chosen
            and "pubchem_cid" in frame.columns
            and pd.notna(row.get("pubchem_cid"))
        ):
          try:
            chosen = f"CID_{int(float(row.get('pubchem_cid')))}"
          except Exception:
            pass
        if not chosen:
          chosen = (
              (prim_val[:17].rstrip(" -_,([") + "…")
              if len(prim_val) > 18
              else (prim_val or f"CMPD_{len(raw_Labels)+1:04d}")
          )
        raw_Labels.append(chosen)

      # Disambiguate any duplicate short labels (e.g. two papers both having 'Cmpd 6a') using ChEMBL/compound ID
      counts = pd.Series(raw_Labels).value_counts()
      dups = set(counts[counts > 1].index)
      if dups:
        disambiguated: list[str] = []
        seen_counts: dict[str, int] = {}
        for i_r, lbl in enumerate(raw_Labels):
          if lbl in dups:
            row = frame.iloc[i_r]
            alt_id = ""
            for alt_c in (
                "chembl_id",
                "molecule_chembl_id",
                "compound_id",
                "patent_example_id",
            ):
              if alt_c in frame.columns and pd.notna(row.get(alt_c)):
                av = str(row.get(alt_c, "")).strip()
                if (
                    av
                    and av.lower() not in ("nan", "none")
                    and av != lbl
                    and len(av) <= 18
                ):
                  alt_id = av
                  break
            if alt_id:
              lbl = f"{lbl} ({alt_id})"
            else:
              seen_counts[lbl] = seen_counts.get(lbl, 0) + 1
              lbl = f"{lbl}#{seen_counts[lbl]}"
          disambiguated.append(lbl)
        raw_Labels = disambiguated

      return pd.Series(raw_Labels, index=frame.index, dtype=str)

    df["_resolved_id"] = _resolve_short_compound_ids(df, args.id_col)
    resolved_id_col = "_resolved_id"

    print(f"Loaded {len(df)} compounds from {args.input_csv}")
    print(f"Using SMILES column: '{smiles_col}', Activity column: '{act_col}'")

    # 1. Standard Matched Molecular Pair Analysis (Hussain-Rea 1-cut & 2-cut)
    mmpa_df = find_matched_molecular_pairs(
        df,
        smiles_col=smiles_col,
        activity_col=act_col,
        id_col=resolved_id_col,
        max_cuts=args.max_cuts,
        min_core_heavies=args.min_core_heavies,
        max_substituent_heavies=args.max_substituent_heavies,
    )
    print(
        f"\n### 1. Matched Molecular Pairs (N={len(mmpa_df)} sharing core"
        f" >= {args.min_core_heavies} heavy atoms)\n"
    )
    if not mmpa_df.empty:
      disp_cols = ["Compound_1", "Compound_2", "Transformation"]
      if "Delta_Activity" in mmpa_df.columns:
        disp_cols.extend(["Activity_1", "Activity_2", "Delta_Activity"])
      tbl_mmpa = _df_to_gfm(mmpa_df[disp_cols].head(args.top_n))
      print(tbl_mmpa)
      md_sections.append(
          "## 1. Top Matched Molecular Pair Transformations"
          f" (N={len(mmpa_df)})\n\n"
          + tbl_mmpa
          + "\n"
      )

      if args.output_mmpa_csv:
        os.makedirs(
            os.path.dirname(os.path.abspath(args.output_mmpa_csv)) or ".",
            exist_ok=True,
        )
        mmpa_df.to_csv(args.output_mmpa_csv, index=False)
        print(f"\nSaved MMPA transformations to: {args.output_mmpa_csv}")

    # 2. General Topological & Ring-Transformation Pair Analysis
    topo_pairs = find_topological_matched_pairs(
        df,
        smiles_col=smiles_col,
        activity_col=act_col,
        id_col=resolved_id_col,
        group_col=args.group_col,
        core_smarts=args.core_smarts,
        min_similarity=args.min_similarity,
        max_pairs_per_compound=args.max_pairs_per_compound,
    )
    if (
        not topo_pairs.empty
        or args.match_topological_pairs
        or args.output_paired_csv
    ):
      print(
          "\n### 2. Topological & Ring-Transformation Matched Pairs"
          f" (N={len(topo_pairs)})\n"
      )
      if not topo_pairs.empty:
        t_disp = [
            "Compound_1",
            "Compound_2",
            "Transformation_Type",
            "Topological_Similarity",
            "Delta_Rings",
            "Delta_RotB",
            "Delta_CNS_MPO",
        ]
        if "Delta_Activity" in topo_pairs.columns:
          t_disp.extend([
              "Activity_1",
              "Activity_2",
              "Delta_Activity",
              "Fold_Change",
          ])
        t_disp = [c for c in t_disp if c in topo_pairs.columns]
        tbl_topo = _df_to_gfm(topo_pairs[t_disp].head(args.top_n))
        print(tbl_topo)
        md_sections.append(
            "## 2. Topological & Ring-Transformation Matched Pairs"
            f" (N={len(topo_pairs)})\n\n"
            + tbl_topo
            + "\n"
        )

        if args.output_paired_csv:
          os.makedirs(
              os.path.dirname(os.path.abspath(args.output_paired_csv)) or ".",
              exist_ok=True,
          )
          topo_pairs.to_csv(args.output_paired_csv, index=False)
          print(
              f"\nSaved topological matched pairs to: {args.output_paired_csv}"
          )

    # 3. SALI Activity Cliff Detection
    if act_col and act_col in df.columns:
      cliffs_df = calculate_sali_dataframe(
          df,
          smiles_col=smiles_col,
          activity_col=act_col,
          id_col=resolved_id_col,
          min_similarity=args.min_similarity,
      )
      out_sali = args.output_sali_csv or args.output_csv
      if out_sali:
        os.makedirs(
            os.path.dirname(os.path.abspath(out_sali)) or ".", exist_ok=True
        )
        cliffs_df.to_csv(out_sali, index=False)
        print(f"\nSaved {len(cliffs_df)} SALI activity cliffs to {out_sali}")

      print(
          "\n### 3. Top SALI Activity Cliffs (min"
          f" similarity={args.min_similarity})\n"
      )
      if not cliffs_df.empty:
        sali_disp = [
            c
            for c in (
                "Compound_1",
                "Compound_2",
                "Similarity",
                "Delta_Activity",
                "SALI",
            )
            if c in cliffs_df.columns
        ]
        tbl_sali = _df_to_gfm(cliffs_df[sali_disp].head(args.top_n))
        print(tbl_sali)
        md_sections.append(
            f"## 3. SALI Activity Cliffs (N={len(cliffs_df)})\n\n"
            + tbl_sali
            + "\n"
        )
  elif args.smiles and args.activities:
    cliffs = calculate_sali_matrix(
        args.smiles, args.activities, min_similarity=args.min_similarity
    )
    cliffs_df = pd.DataFrame(cliffs)
    out_sali = args.output_sali_csv or args.output_csv
    if out_sali:
      os.makedirs(
          os.path.dirname(os.path.abspath(out_sali)) or ".", exist_ok=True
      )
      cliffs_df.to_csv(out_sali, index=False)
      print(f"Saved {len(cliffs_df)} SALI activity cliffs to {out_sali}")

    print(f"Top SALI Activity Cliffs (min similarity={args.min_similarity}):")
    for c in cliffs[:5]:
      print(
          f"  SALI: {c['sali']:6.2f} | ΔAct: {c['delta_activity']:4.2f} | Sim:"
          f" {c['similarity']:.2f}"
      )
      print(f"    Mol A (idx {c['idx_i']}): {c['smiles_i']}")
      print(f"    Mol B (idx {c['idx_j']}): {c['smiles_j']}")
  else:
    parser.error(
        "Either --input_csv or both --smiles and --activities must be provided."
    )

  if args.plot_png:
    plot_file = os.path.abspath(args.plot_png)
    os.makedirs(os.path.dirname(plot_file) or ".", exist_ok=True)
    try:
      target_lbl = act_col if "act_col" in locals() and act_col else "pIC50"
      mmpa_to_plot = (
          mmpa_df if "mmpa_df" in locals() and not mmpa_df.empty else None
      )
      sali_to_plot = (
          cliffs_df if "cliffs_df" in locals() and not cliffs_df.empty else None
      )
      paired_to_plot = (
          topo_pairs
          if "topo_pairs" in locals() and not topo_pairs.empty
          else None
      )

      plot_archetypes.plot_mmpa_sar_landscape(
          mmpa_df=mmpa_to_plot,
          sali_df=sali_to_plot,
          output_png=plot_file,
          target_label=target_lbl,
          paired_shifts_df=paired_to_plot,
          paired_x_col=args.paired_x_col,
          paired_y_col=args.paired_y_col,
      )
      print(
          "\n[Plotting] Saved publication-ready MMPA & SALI SAR landscape to:"
          f" {plot_file}"
      )
      print(f"  Markdown embed: ![MMPA SAR Landscape]({plot_file})")
      md_sections.append(
          "## 4. MMPA & SALI SAR Landscape Visualization\n\n![MMPA SAR"
          f" Landscape]({plot_file})\n"
      )
    except Exception as e:
      print(
          f"Warning: Could not render MMPA landscape plot: {e}",
          file=sys.stderr,
      )

  if args.output_md:
    md_path = os.path.abspath(args.output_md)
    os.makedirs(os.path.dirname(md_path) or ".", exist_ok=True)
    with open(md_path, "w", encoding="utf-8") as fh:
      fh.write("\n".join(md_sections))
    print(f"Saved Markdown report to: {md_path}")


if __name__ == "__main__":
  main()
