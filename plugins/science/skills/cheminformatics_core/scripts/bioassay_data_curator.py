# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "polite-http",
#     "rdkit",
#     "numpy",
#     "pandas",
#     "scipy",
#     "scikit-learn",
#     "matplotlib",
# ]
# ///

"""Bioassay Data Curator & Batch Ingestion CLI.

Fetches, normalizes, and locally caches bioassay and structure-activity
relationship (SAR) datasets from ChEMBL, PubChem, and PubMed with automated unit
conversion (uM/nM -> pKi/pIC50), multi-endpoint consolidation, RDKit descriptor
calculation, and persistent disk caching.
"""

import argparse
import collections
import io
import json
import math
import os
import re
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from polite_http import http_client
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, MACCSkeys, rdFMCS, rdFingerprintGenerator, rdMolDescriptors

_CHEMBL_CLIENT = http_client.HttpClient(
    'https://www.ebi.ac.uk/chembl/api/data', qps=5.0
)
_PUBCHEM_CLIENT = http_client.HttpClient(
    'https://pubchem.ncbi.nlm.nih.gov/rest/pug', qps=5.0
)

# Matches single-letter amino acid point mutations (e.g. L858R, T790M, G12C) in assay descriptions.
_AMINO_ACID_POINT_MUTATION_RE = re.compile(
    r'\b[acdefghiklmnpqrstvwy]\d{2,4}[acdefghiklmnpqrstvwy]\b'
)

# Canonical CNS MPO implementation (single source of truth) and companion doc2mol fetcher.
_script_dir = os.path.dirname(os.path.abspath(__file__))
if _script_dir not in sys.path:
  sys.path.insert(0, _script_dir)

_doc2mol_dir = os.path.abspath(
    os.path.join(_script_dir, '..', '..', 'cheminformatics_doc2mol', 'scripts')
)
if os.path.isdir(_doc2mol_dir) and _doc2mol_dir not in sys.path:
  sys.path.append(_doc2mol_dir)

import curate_congeneric_corpus
from descriptors import calculate_cns_mpo
import paper_fetcher


def fetch_chembl_document_by_pmid(
    pmid: str, timeout: int = 10
) -> Optional[str]:
  """Resolves a PubMed ID to a ChEMBL document ID."""
  data = _CHEMBL_CLIENT.fetch_json(
      f'document.json?pubmed_id={pmid}', timeout=timeout
  )
  docs = data.get('documents', [])
  if docs:
    return docs[0].get('document_chembl_id')
  return None


def _clean_doc_year(raw_val: Any) -> Optional[int]:
  """Extracts a valid calendar year in [1800, 2100], rejecting booleans."""
  if raw_val is None or isinstance(raw_val, (bool, np.bool_)):
    return None
  try:
    if pd.isna(raw_val):
      return None
  except (TypeError, ValueError):
    pass
  try:
    fval = float(raw_val)
    if not math.isfinite(fval):
      return None
    yr = int(round(fval))
    return yr if 1800 <= yr <= 2100 else None
  except (TypeError, ValueError):
    m_yr = re.match(r'^\s*(18\d{2}|19\d{2}|20\d{2})\b', str(raw_val))
    return int(m_yr.group(1)) if m_yr else None


def fetch_chembl_activities(
    doc_chembl_id: str, timeout: int = 10, limit: int = 1000
) -> List[Dict[str, Any]]:
  """Fetches all bioactivities associated with a ChEMBL document ID."""
  activities = []
  offset = 0
  while True:
    path = f'activity.json?document_chembl_id={doc_chembl_id}&limit={limit}&offset={offset}'
    data = _CHEMBL_CLIENT.fetch_json(path, timeout=timeout)
    batch = data.get('activities', [])
    activities.extend(batch)
    total = data.get('page_meta', {}).get('total_count', len(activities))
    if len(activities) >= total or not batch:
      break
    offset += len(batch)
  if activities and any(
      _clean_doc_year(a.get('document_year')) is None for a in activities
  ):
    try:
      doc_meta = _CHEMBL_CLIENT.fetch_json(
          f'document/{doc_chembl_id}.json', timeout=timeout
      )
      doc_yr = _clean_doc_year((doc_meta or {}).get('year'))
      if doc_yr is not None:
        for a in activities:
          if _clean_doc_year(a.get('document_year')) is None:
            a['document_year'] = doc_yr
    except (
        http_client.HttpError,
        ValueError,
        KeyError,
        TypeError,
        OSError,
    ):
      pass
  return activities


def fetch_pubchem_assay_csv(aid: str, timeout: int = 10) -> pd.DataFrame:
  """Fetches and parses a PubChem BioAssay CSV by AID."""
  raw = _PUBCHEM_CLIENT.fetch_bytes(f'assay/aid/{aid}/CSV', timeout=timeout)
  df = pd.read_csv(io.StringIO(raw.decode('utf-8')))
  return df


def compute_molecule_descriptors(smiles: str) -> Dict[str, Any]:
  """Computes 2D physicochemical descriptors, CNS MPO, and macrocycle indicators."""
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    return {}

  ring_info = mol.GetRingInfo()
  ring_sizes = [len(r) for r in ring_info.AtomRings()]
  macro_sizes = [s for s in ring_sizes if s >= 12]
  is_macro = 1 if len(macro_sizes) > 0 else 0
  macro_size = max(macro_sizes) if macro_sizes else 0

  chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
  chiral_tag = '; '.join([f'{idx}:{tag}' for idx, tag in chiral_centers])

  mw = Descriptors.MolWt(mol)
  logp = Crippen.MolLogP(mol)
  tpsa = rdMolDescriptors.CalcTPSA(mol)
  hbd = Lipinski.NumHDonors(mol)
  hba = Lipinski.NumHAcceptors(mol)
  rotb = Lipinski.NumRotatableBonds(mol)
  fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
  cns_mpo = calculate_cns_mpo(mol)

  return {
      'MW': round(mw, 2),
      'LogP': round(logp, 2),
      'TPSA': round(tpsa, 2),
      'HBD': hbd,
      'HBA': hba,
      'RotB': rotb,
      'Fsp3': round(fsp3, 3),
      'CNS_MPO': round(cns_mpo, 2),
      'Is_Macrocycle': is_macro,
      'Macrocycle_Size': macro_size,
      'Chiral_Centers': chiral_tag,
  }


def curate_chembl_activities(activities: List[Dict[str, Any]]) -> pd.DataFrame:
  """Normalizes and consolidates ChEMBL activity records into a single wide DataFrame."""
  mol_dict = {}

  for a in activities:
    smi = a.get('canonical_smiles')
    mid = a.get('molecule_chembl_id')
    if not smi or not mid:
      continue

    cand_yr = _clean_doc_year(a.get('document_year'))
    if mid not in mol_dict:
      mol_dict[mid] = {
          'compound_id': mid,
          'canonical_smiles': smi,
          'pref_name': a.get('molecule_pref_name') or mid,
          'document_chembl_id': a.get('document_chembl_id'),
          'document_year': cand_yr,
      }
    elif cand_yr is not None and (
        mol_dict[mid].get('document_year') is None
        or cand_yr < mol_dict[mid]['document_year']
    ):
      mol_dict[mid]['document_year'] = cand_yr
      mol_dict[mid]['document_chembl_id'] = a.get('document_chembl_id')

    st = a.get('standard_type', '')
    val_str = a.get('standard_value')
    units = (a.get('standard_units') or '').lower()
    desc = (a.get('assay_description') or '').lower()

    if val_str is None:
      continue

    try:
      val = float(val_str)
    except ValueError:
      continue

    # In ChEMBL, standard_units are typically nM
    # Calculate pActivity (-log10 M)
    if units.lower() in ('nm', 'nmol/l'):
      p_val = 9.0 - math.log10(max(val, 1e-12))
    elif units.lower() in ('um', 'umol/l', 'micromolar'):
      p_val = 6.0 - math.log10(max(val, 1e-12))
    elif units.lower() in ('m', 'mol/l'):
      p_val = -math.log10(max(val, 1e-12))
    elif units.lower() in ('pm', 'pmol/l'):
      p_val = 12.0 - math.log10(max(val, 1e-12))
    else:
      # Default ChEMBL standard_value assumption is nM
      p_val = 9.0 - math.log10(max(val, 1e-12))

    # Classify assay target endpoint
    if st in ('Ki', 'IC50', 'EC50', 'Kd'):
      prefix = st
      if (
          'mutant' in desc
          or 'mutation' in desc
          or 'variant' in desc
          or bool(_AMINO_ACID_POINT_MUTATION_RE.search(desc))
      ):
        col_raw = f'{prefix}_mutant_nM'
        col_p = f'p{prefix}_mutant'
      elif any(
          t in desc
          for t in (
              'cell',
              'fusion',
              'viability',
              'proliferation',
              'cytotox',
              'growth',
          )
      ):
        col_raw = f'{prefix}_cell_nM'
        col_p = f'p{prefix}_cell'
      else:
        col_raw = f'{prefix}_wt_nM'
        col_p = f'p{prefix}_wt'

      mol_dict[mid][col_raw] = val
      mol_dict[mid][col_p] = round(p_val, 3)

    elif st == 'LogD':
      mol_dict[mid]['LogD'] = val
    elif 'ratio' in st.lower() or 'efflux' in desc:
      mol_dict[mid]['Efflux_Ratio'] = val

  df = pd.DataFrame(list(mol_dict.values()))
  if df.empty:
    return df

  # Compute RDKit descriptors for all molecules
  desc_rows = [compute_molecule_descriptors(s) for s in df['canonical_smiles']]
  desc_df = pd.DataFrame(desc_rows)
  curated_df = pd.concat(
      [df.reset_index(drop=True), desc_df.reset_index(drop=True)], axis=1
  )

  return curated_df


def detect_activity_cliffs(
    df: pd.DataFrame,
    smi_col: str = 'canonical_smiles',
    act_col: str = 'pActivity',
    id_col: str = 'compound_id',
    sim_threshold: float = 0.8,
    cliff_threshold: float = 2.0,
) -> Tuple[pd.DataFrame, List[Dict[str, Any]]]:
  """Detects activity cliffs between pairs of compounds with high structural similarity."""
  mols = []
  valid_indices = []
  for idx, smi in enumerate(df[smi_col]):
    m = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None
    if m is not None:
      mols.append(m)
      valid_indices.append(idx)

  if len(mols) < 2:
    return df, []

  fcfp_gen = rdFingerprintGenerator.GetMorganGenerator(
      radius=2,
      fpSize=2048,
      atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen(),
  )
  ecfp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

  fps_fcfp = [fcfp_gen.GetFingerprint(m) for m in mols]
  fps_ecfp = [ecfp_gen.GetFingerprint(m) for m in mols]
  fps_maccs = [MACCSkeys.GenMACCSKeys(m) for m in mols]
  fps_rdk = [Chem.RDKFingerprint(m) for m in mols]

  n = len(mols)
  cliffs = []
  cliff_flags = [False] * len(df)
  cliff_partners = collections.defaultdict(list)
  cliff_max_delta = collections.defaultdict(float)
  cliff_max_sim = collections.defaultdict(float)

  for i in range(n):
    idx_i = valid_indices[i]
    id_i = str(df.loc[idx_i, id_col])
    act_i = df.loc[idx_i, act_col]
    if pd.isna(act_i):
      continue
    act_i = float(act_i)

    for j in range(i + 1, n):
      idx_j = valid_indices[j]
      id_j = str(df.loc[idx_j, id_col])
      act_j = df.loc[idx_j, act_col]
      if pd.isna(act_j):
        continue
      act_j = float(act_j)

      delta_act = abs(act_i - act_j)

      sim_fcfp = DataStructs.TanimotoSimilarity(fps_fcfp[i], fps_fcfp[j])
      sim_ecfp = DataStructs.TanimotoSimilarity(fps_ecfp[i], fps_ecfp[j])
      sim_maccs = DataStructs.TanimotoSimilarity(fps_maccs[i], fps_maccs[j])
      sim_rdk = DataStructs.TanimotoSimilarity(fps_rdk[i], fps_rdk[j])

      # MCS atom ratio
      res = rdFMCS.FindMCS([mols[i], mols[j]])
      num_mcs = res.numAtoms
      sim_mcs = num_mcs / max(
          1, (mols[i].GetNumHeavyAtoms() + mols[j].GetNumHeavyAtoms() - num_mcs)
      )

      max_sim = max(sim_fcfp, sim_maccs, sim_rdk, sim_mcs)

      if (max_sim >= sim_threshold or sim_ecfp >= sim_threshold) and (
          delta_act >= cliff_threshold
      ):
        cliff_flags[idx_i] = True
        cliff_flags[idx_j] = True
        cliff_partners[idx_i].append(
            f'{id_j} (Δ={delta_act:.2f}, Sim={max_sim:.2f})'
        )
        cliff_partners[idx_j].append(
            f'{id_i} (Δ={delta_act:.2f}, Sim={max_sim:.2f})'
        )
        cliff_max_delta[idx_i] = max(cliff_max_delta[idx_i], delta_act)
        cliff_max_delta[idx_j] = max(cliff_max_delta[idx_j], delta_act)
        cliff_max_sim[idx_i] = max(cliff_max_sim[idx_i], max_sim)
        cliff_max_sim[idx_j] = max(cliff_max_sim[idx_j], max_sim)

        cliffs.append({
            'Compound_1': id_i,
            'Compound_2': id_j,
            'SMILES_1': df.loc[idx_i, smi_col],
            'SMILES_2': df.loc[idx_j, smi_col],
            'pActivity_1': act_i,
            'pActivity_2': act_j,
            'Delta_Activity': round(delta_act, 3),
            'Similarity_FCFP4': round(sim_fcfp, 3),
            'Similarity_MACCS': round(sim_maccs, 3),
            'Similarity_RDKit': round(sim_rdk, 3),
            'Similarity_MCS': round(sim_mcs, 3),
            'Similarity_ECFP4': round(sim_ecfp, 3),
            'Max_Similarity': round(max_sim, 3),
        })

  df['is_activity_cliff'] = cliff_flags
  df['activity_cliff_partner'] = [
      '; '.join(cliff_partners[i]) if cliff_partners[i] else 'None'
      for i in range(len(df))
  ]
  df['cliff_max_delta'] = [
      round(cliff_max_delta[i], 3) if i in cliff_max_delta else 0.0
      for i in range(len(df))
  ]
  df['cliff_max_similarity'] = [
      round(cliff_max_sim[i], 3) if i in cliff_max_sim else 0.0
      for i in range(len(df))
  ]

  return df, cliffs


def main():
  parser = argparse.ArgumentParser(
      description='Bioassay Data Curator & Batch Ingestion CLI'
  )
  parser.add_argument(
      '--input',
      '--input_csv',
      type=str,
      dest='input_csv',
      help='Path to local raw bioassay CSV dataset',
  )
  parser.add_argument(
      '--chembl_doc',
      '--chembl-doc',
      '--chembl_docs',
      '--chembl-docs',
      dest='chembl_doc',
      nargs='*',
      help='ChEMBL Document ID(s) (e.g. CHEMBL<DOC_ID>)',
  )
  parser.add_argument(
      '--pmid',
      '--pmids',
      dest='pmid',
      nargs='*',
      help='PubMed ID(s) (e.g. <PMID_1> <PMID_2>)',
  )
  parser.add_argument(
      '--doi',
      '--dois',
      dest='dois',
      nargs='*',
      help='DOI(s) (e.g. 10.1021/<DOI_SUFFIX>)',
  )
  parser.add_argument(
      '--pubchem_aid',
      '--pubchem-aid',
      dest='pubchem_aid',
      type=str,
      help='PubChem Assay ID (e.g. 1153093)',
  )
  parser.add_argument(
      '--output',
      '--output_csv',
      '--output-csv',
      type=str,
      dest='output_csv',
      required=True,
      help='Path to output curated CSV dataset',
  )
  parser.add_argument(
      '--output_json', type=str, help='Optional path to cache raw JSON response'
  )
  parser.add_argument(
      '--cliff_threshold',
      type=float,
      default=2.0,
      help='Potency delta threshold for activity cliffs (default: 2.0)',
  )
  parser.add_argument(
      '--sim_threshold',
      type=float,
      default=0.8,
      help='Structural similarity threshold for activity cliffs (default: 0.8)',
  )
  parser.add_argument(
      '--timeout', type=int, default=10, help='HTTP request timeout in seconds'
  )

  args = parser.parse_args()

  if args.input_csv and os.path.exists(args.input_csv):
    print(f'Curating local bioassay dataset: {args.input_csv}...')
    df_raw = pd.read_csv(args.input_csv)
    smi_col = None
    for c in ['canonical_smiles', 'smiles', 'SMILES', 'Smiles', 'structure']:
      if c in df_raw.columns:
        smi_col = c
        break
    if not smi_col:
      raise ValueError(f'Could not find SMILES column in {args.input_csv}')

    id_col = None
    for c in [
        'molecule_id',
        'compound_id',
        'molecule_chembl_id',
        'ID',
        'id',
        'Name',
        'name',
    ]:
      if c in df_raw.columns:
        id_col = c
        break
    if not id_col:
      id_col = 'compound_id'
      df_raw[id_col] = [f'CMPD_{i+1}' for i in range(len(df_raw))]

    # Find activity column
    act_col = None
    for c in [
        'raw_value',
        'value',
        'pIC50',
        'pActivity',
        'pChEMBL',
        'standard_value',
        'IC50',
        'activity',
        'Ki',
        'EC50',
    ]:
      if c in df_raw.columns:
        act_col = c
        break

    units_col = None
    for c in ['units', 'unit', 'standard_units', 'Units']:
      if c in df_raw.columns:
        units_col = c
        break

    df_curated = df_raw.copy()
    if 'compound_id' not in df_curated.columns:
      df_curated['compound_id'] = df_curated[id_col]
    if 'canonical_smiles' not in df_curated.columns:
      df_curated['canonical_smiles'] = [
          Chem.MolToSmiles(
              Chem.MolFromSmiles(str(s)), canonical=True, isomericSmiles=True
          )
          if Chem.MolFromSmiles(str(s)) is not None
          else str(s)
          for s in df_curated[smi_col]
      ]
    df_curated = curate_congeneric_corpus.backfill_document_dates(df_curated)
    if 'document_year' not in df_curated.columns:
      df_curated['document_year'] = pd.NA

    def calc_pactivity(row):
      if act_col is None:
        return np.nan
      val = row[act_col]
      if pd.isna(val):
        return np.nan
      try:
        v = float(val)
      except (ValueError, TypeError):
        return np.nan
      if 'p' in act_col.lower() and act_col.lower() in (
          'pic50',
          'pactivity',
          'pchembl',
          'pki',
          'pec50',
      ):
        return round(v, 3)
      u = (
          str(row[units_col]).strip().lower()
          if units_col and pd.notna(row.get(units_col))
          else ''
      )
      if u in ('nm', 'nmol/l', 'nanomolar'):
        molar = v * 1e-9
      elif u in ('um', 'umol/l', 'micromolar', 'µm', 'umolar'):
        molar = v * 1e-6
      elif u in ('mm', 'mmol/l', 'millimolar'):
        molar = v * 1e-3
      elif u in ('pm', 'pmol/l', 'picomolar'):
        molar = v * 1e-12
      elif u in ('m', 'mol/l', 'molar'):
        molar = v
      else:
        if v > 1.0:
          molar = v * 1e-9
        elif v > 1e-3:
          molar = v * 1e-6
        else:
          molar = v
      return round(float(-np.log10(max(molar, 1e-15))), 3)

    df_curated['pActivity'] = [
        calc_pactivity(r) for _, r in df_curated.iterrows()
    ]
    if 'activity_type' in df_curated.columns:
      for idx, r in df_curated.iterrows():
        atype = str(r['activity_type']).upper()
        if atype == 'IC50':
          df_curated.loc[idx, 'pIC50'] = df_curated.loc[idx, 'pActivity']
          df_curated.loc[idx, 'pIC50_wt'] = df_curated.loc[idx, 'pActivity']
        elif atype == 'KI':
          df_curated.loc[idx, 'pKi'] = df_curated.loc[idx, 'pActivity']
          df_curated.loc[idx, 'pKi_wt'] = df_curated.loc[idx, 'pActivity']

    desc_rows = [
        compute_molecule_descriptors(str(s))
        for s in df_curated['canonical_smiles']
    ]
    desc_df = pd.DataFrame(desc_rows)
    for col in desc_df.columns:
      if col not in df_curated.columns:
        df_curated[col] = desc_df[col]

    df_curated, cliffs = detect_activity_cliffs(
        df_curated,
        smi_col='canonical_smiles',
        act_col='pActivity',
        id_col='compound_id',
        sim_threshold=args.sim_threshold,
        cliff_threshold=args.cliff_threshold,
    )

    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True
    )
    df_curated.to_csv(args.output_csv, index=False)
    print(f'Curated {len(df_curated)} compounds to {args.output_csv}')
    print(
        f'Activity cliffs detected: {len(cliffs)} matched pair(s) with Sim >='
        f' {args.sim_threshold} and Δ >= {args.cliff_threshold} log units.'
    )
    if cliffs:
      print(
          pd.DataFrame(cliffs)[[
              'Compound_1',
              'Compound_2',
              'pActivity_1',
              'pActivity_2',
              'Delta_Activity',
              'Similarity_FCFP4',
              'Similarity_MACCS',
              'Similarity_RDKit',
              'Similarity_MCS',
          ]].to_string(index=False)
      )
    return

  doc_tokens: List[str] = []
  for grp in (args.chembl_doc, args.pmid, getattr(args, 'dois', None)):
    if grp:
      if isinstance(grp, str):
        grp = [grp]
      for item in grp:
        for part in str(item).split(','):
          if part.strip():
            doc_tokens.append(part.strip())

  if doc_tokens:
    all_activities: List[Dict[str, Any]] = []
    for tok in doc_tokens:
      doc_id = tok
      if not doc_id.upper().startswith('CHEMBL') and '.' not in doc_id:
        clean_pmid = re.sub(r'^(?i:pmid\s*:?\s*)', '', doc_id)
        print(f'Resolving PMID {clean_pmid} to ChEMBL document ID...')
        resolved = fetch_chembl_document_by_pmid(
            clean_pmid, timeout=args.timeout
        )
        if resolved:
          print(f'Found ChEMBL Document ID: {resolved}')
          doc_id = resolved
        else:
          print(
              f'Could not resolve PMID {clean_pmid} in ChEMBL document index.'
          )
          doc_id = ''
      if doc_id and doc_id.upper().startswith('CHEMBL'):
        print(f'Fetching bioactivities for ChEMBL document {doc_id}...')
        acts = fetch_chembl_activities(doc_id, timeout=args.timeout)
        print(f'Fetched {len(acts)} activity records for {doc_id}.')
        all_activities.extend(acts)

    if args.output_json and all_activities:
      os.makedirs(
          os.path.dirname(os.path.abspath(args.output_json)), exist_ok=True
      )
      with open(args.output_json, 'w', encoding='utf-8') as f:
        json.dump(all_activities, f, indent=2)
      print(f'Cached raw JSON to {args.output_json}')

    df_curated = (
        curate_chembl_activities(all_activities)
        if all_activities
        else pd.DataFrame()
    )
    try:
      df_wide = paper_fetcher.fetch_chembl_document_activities(doc_tokens)
      if not df_wide.empty:
        if df_curated.empty:
          df_curated = df_wide
        elif (
            'compound_id' in df_curated.columns
            and 'compound_id' in df_wide.columns
        ):
          wide_dedup = df_wide.groupby(
              'compound_id', sort=False, as_index=False
          ).first()
          if (
              'document_year' in wide_dedup.columns
              and 'document_year' in df_curated.columns
          ):
            yr_map = {
                cid: yr
                for cid, raw_y in zip(
                    wide_dedup['compound_id'], wide_dedup['document_year']
                )
                if (yr := _clean_doc_year(raw_y)) is not None
            }
            if yr_map:
              df_curated['document_year'] = [
                  cur_y
                  if _clean_doc_year(cur_y) is not None
                  else yr_map.get(cid, cur_y)
                  for cid, cur_y in zip(
                      df_curated['compound_id'], df_curated['document_year']
                  )
              ]
          extra_cols = [
              c
              for c in wide_dedup.columns
              if c not in df_curated.columns or c == 'compound_id'
          ]
          if len(extra_cols) > 1:
            df_curated = df_curated.merge(
                wide_dedup[extra_cols], on='compound_id', how='left'
            )
    except Exception as exc:
      print(f'Note: paper_fetcher enrichment skipped ({exc})', file=sys.stderr)

    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True
    )
    df_curated.to_csv(args.output_csv, index=False)
    print(f'Curated {len(df_curated)} unique compounds to {args.output_csv}')
    print(f'Columns ({len(df_curated.columns)}): {list(df_curated.columns)}')

  elif args.pubchem_aid:
    print(f'Fetching PubChem BioAssay AID {args.pubchem_aid}...')
    df_aid = fetch_pubchem_assay_csv(args.pubchem_aid, timeout=args.timeout)
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True
    )
    df_aid.to_csv(args.output_csv, index=False)
    print(f'Saved PubChem AID dataset to {args.output_csv}')

  else:
    print(
        'Error: Must specify --chembl_doc, --pmid, --pubchem_aid, or --input.',
        file=sys.stderr,
    )
    sys.exit(1)


if __name__ == '__main__':
  main()
