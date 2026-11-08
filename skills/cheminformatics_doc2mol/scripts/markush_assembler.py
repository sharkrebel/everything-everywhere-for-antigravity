# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pandas",
#     "rdkit",
# ]
# ///

"""Patent Markush Combinatorial Structure Assembler.

Assembles generic core scaffolds with variable R-group substituent tables
from patent exemplified compounds and literature series into validated RDKit
Mols.
"""

import argparse
import itertools
import os
import re
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

COMMON_SUBSTITUENTS = {
    "h": "[H]",
    "c": "C[*:X]",
    "me": "C[*:X]",
    "ch3": "C[*:X]",
    "et": "CC[*:X]",
    "ch2ch3": "CC[*:X]",
    "ipr": "[*:X]C(C)C",
    "i-pr": "[*:X]C(C)C",
    "ch(ch3)2": "[*:X]C(C)C",
    "npr": "[*:X]CCC",
    "n-pr": "[*:X]CCC",
    "propyl": "[*:X]CCC",
    "tbu": "[*:X]C(C)(C)C",
    "t-bu": "[*:X]C(C)(C)C",
    "c(ch3)3": "[*:X]C(C)(C)C",
    "tert-butyl": "[*:X]C(C)(C)C",
    "nbu": "[*:X]CCCC",
    "n-bu": "[*:X]CCCC",
    "butyl": "[*:X]CCCC",
    "ibu": "[*:X]CC(C)C",
    "i-bu": "[*:X]CC(C)C",
    "isobutyl": "[*:X]CC(C)C",
    "sbu": "[*:X]C(C)CC",
    "s-bu": "[*:X]C(C)CC",
    "cpr": "[*:X]C1CC1",
    "c-pr": "[*:X]C1CC1",
    "cyclopropyl": "[*:X]C1CC1",
    "cbu": "[*:X]C1CCC1",
    "c-bu": "[*:X]C1CCC1",
    "cyclobutyl": "[*:X]C1CCC1",
    "cpent": "[*:X]C1CCCC1",
    "c-pent": "[*:X]C1CCCC1",
    "cyclopentyl": "[*:X]C1CCCC1",
    "chex": "[*:X]C1CCCCC1",
    "c-hex": "[*:X]C1CCCCC1",
    "cy": "[*:X]C1CCCCC1",
    "cyclohexyl": "[*:X]C1CCCCC1",
    "cf3": "C(F)(F)(F)[*:X]",
    "ch2cf3": "C(F)(F)(F)C[*:X]",
    # Alkoxy groups attach through OXYGEN. Putting the attachment point on the
    # carbon instead yields an alcohol with the same molecular formula (e.g.
    # a benzylic alcohol rather than an anisole), which MW checks cannot catch.
    "ome": "[*:X]OC",
    "och3": "[*:X]OC",
    "oet": "[*:X]OCC",
    "onpr": "[*:X]OCCC",
    "o-npr": "[*:X]OCCC",
    "n-propoxy": "[*:X]OCCC",
    "propoxy": "[*:X]OCCC",
    "oipr": "[*:X]OC(C)C",
    "otbu": "[*:X]OC(C)(C)C",
    "ocf3": "[*:X]OC(F)(F)F",
    "f": "F[*:X]",
    "cl": "Cl[*:X]",
    "br": "Br[*:X]",
    "i": "I[*:X]",
    "cn": "C(#N)[*:X]",
    "c#n": "C(#N)[*:X]",
    # Nitro attaches through NITROGEN. The previous form attached via the
    # anionic oxygen, which is both an unparseable valence and the wrong
    # connectivity (a nitrite ester).
    "no2": "[*:X][N+](=O)[O-]",
    "oh": "O[*:X]",
    "nh2": "N[*:X]",
    "nhme": "[*:X]NC",
    "nhet": "[*:X]NCC",
    "nme2": "[*:X]N(C)C",
    "nhet2": "[*:X]N(CC)CC",
    "ph": "[*:X]c1ccccc1",
    "phenyl": "[*:X]c1ccccc1",
    "bn": "[*:X]Cc1ccccc1",
    "benzyl": "[*:X]Cc1ccccc1",
    "tol": "[*:X]c1ccc(C)cc1",
    "tolyl": "[*:X]c1ccc(C)cc1",
    "p-tolyl": "[*:X]c1ccc(C)cc1",
    "4-tolyl": "[*:X]c1ccc(C)cc1",
    "4-f-ph": "c1cc(F)ccc1[*:X]",
    "4-fluorophenyl": "c1cc(F)ccc1[*:X]",
    "4-cl-ph": "c1cc(Cl)ccc1[*:X]",
    "4-chlorophenyl": "c1cc(Cl)ccc1[*:X]",
    "4-pyridyl": "c1ccncc1[*:X]",
    "3-pyridyl": "c1cnccc1[*:X]",
    "2-pyridyl": "c1ncccc1[*:X]",
    "ac": "[*:X]C(=O)C",
    "acetyl": "[*:X]C(=O)C",
    "nhac": "[*:X]NC(=O)C",
    "acetamido": "[*:X]NC(=O)C",
    "oac": "[*:X]OC(=O)C",
    "acetoxy": "[*:X]OC(=O)C",
    "co2h": "[*:X]C(=O)O",
    "cooh": "[*:X]C(=O)O",
    "co2me": "[*:X]C(=O)OC",
    "coome": "[*:X]C(=O)OC",
    "co2et": "[*:X]C(=O)OCC",
    "cooet": "[*:X]C(=O)OCC",
    "co2tbu": "[*:X]C(=O)OC(C)(C)C",
    "cootbu": "[*:X]C(=O)OC(C)(C)C",
    "conh2": "[*:X]C(=O)N",
    "conhme": "[*:X]C(=O)NC",
    "conme2": "[*:X]C(=O)N(C)C",
    "so2me": "[*:X]S(=O)(=O)C",
    "ms": "[*:X]S(=O)(=O)C",
    "mesyl": "[*:X]S(=O)(=O)C",
    "so2ph": "[*:X]S(=O)(=O)c1ccccc1",
    "ts": "[*:X]S(=O)(=O)c1ccc(C)cc1",
    "tosyl": "[*:X]S(=O)(=O)c1ccc(C)cc1",
    "tf": "[*:X]S(=O)(=O)C(F)(F)F",
    "triflyl": "[*:X]S(=O)(=O)C(F)(F)F",
    "so2cf3": "[*:X]S(=O)(=O)C(F)(F)F",
    "so2nh2": "[*:X]S(=O)(=O)N",
    "so2nme2": "[*:X]S(=O)(=O)N(C)C",
    "boc": "[*:X]C(=O)OC(C)(C)C",
    "t-boc": "[*:X]C(=O)OC(C)(C)C",
    "nhboc": "[*:X]NC(=O)OC(C)(C)C",
    "cbz": "[*:X]C(=O)OCc1ccccc1",
    "nhcbz": "[*:X]NC(=O)OCc1ccccc1",
    "fmoc": "[*:X]C(=O)OCC1c2ccccc2-c3ccccc13",
    "nhfmoc": "[*:X]NC(=O)OCC1c2ccccc2-c3ccccc13",
    "alloc": "[*:X]C(=O)OCC=C",
    "thp": "[*:X]C1CCCCO1",
    "tms": "[*:X][Si](C)(C)C",
    "tbs": "[*:X][Si](C)(C)C(C)(C)C",
    "tbdms": "[*:X][Si](C)(C)C(C)(C)C",
    "tips": "[*:X][Si](C(C)C)(C(C)C)C(C)C",
    "morpholino": "O1CCN([*:X])CC1",
    "piperidino": "C1CCN([*:X])CC1",
    "pyrrolidino": "C1CCN([*:X])C1",
    "4-methylpiperazin-1-yl": "CN1CCN([*:X])CC1",
    "piperazin-1-yl": "N1CCN([*:X])CC1",
}


def normalize_rgroup_string(raw_val: str, r_idx: int) -> str:
  """Converts raw substituent string into an attachment fragment with [*:r_idx]."""
  if not raw_val or pd.isna(raw_val):
    return ""

  val = str(raw_val).strip()
  if not val:
    return ""
  target_tag = f"[*:{r_idx}]"

  # Strip leading/trailing wildcards, dashes, and whitespace to extract core abbreviation
  core_token = re.sub(
      r"^\[\*:[0-9]+\]\s*|^\[\*\]\s*|^[\*\-]\s*", "", val
  ).strip()
  core_token = re.sub(
      r"\s*\[\*:[0-9]+\]$|\s*\[\*\]$|\s*[\*\-]$", "", core_token
  ).strip()
  if core_token.startswith("(") and core_token.endswith(")"):
    core_token = core_token[1:-1].strip()
  core_lower = core_token.lower()
  if not core_lower:
    return "[H]"

  # 1. Check common substituent dictionary against stripped core token
  if core_lower in COMMON_SUBSTITUENTS:
    frag_smi = COMMON_SUBSTITUENTS[core_lower]
    if frag_smi == "[H]":
      return "[H]"
    return frag_smi.replace("[*:X]", target_tag)

  val_clean = val.lstrip("-").strip()
  val_lower = val_clean.lower()
  if val_lower in COMMON_SUBSTITUENTS:
    frag_smi = COMMON_SUBSTITUENTS[val_lower]
    if frag_smi == "[H]":
      return "[H]"
    return frag_smi.replace("[*:X]", target_tag)

  # 2. Check if already has attachment point [*:X] or [*]
  if "[*" in val:
    return re.sub(r"\[\*:[0-9]+\]|\[\*\]", target_tag, val)

  # 3. Check if starts or ends with attachment dash or star
  if val.startswith("*") or val.startswith("-"):
    smi = val.lstrip("*-").strip()
    smi = re.sub(r"(?<!\[)CF3(?![a-zA-Z0-9\]])", "C(F)(F)F", smi)
    return f"{target_tag}{smi}"
  elif val.endswith("*") or val.endswith("-"):
    smi = val.rstrip("*-").strip()
    smi = re.sub(r"(?<!\[)CF3(?![a-zA-Z0-9\]])", "C(F)(F)F", smi)
    return f"{smi}{target_tag}"

  # 4. Fallback: treat as fragment with attachment at first atom
  smi = re.sub(r"(?<!\[)CF3(?![a-zA-Z0-9\]])", "C(F)(F)F", val)
  return f"{target_tag}{smi}"


def _attach_fragment(
    core: Chem.Mol, frag: Chem.Mol, r_map_num: int
) -> Chem.Mol:
  """Attaches a single fragment to the core at r_map_num."""
  core_dummy = None
  for atom in core.GetAtoms():
    if atom.GetAtomicNum() == 0 and (
        atom.GetAtomMapNum() == r_map_num or atom.GetIsotope() == r_map_num
    ):
      core_dummy = atom
      break
  if not core_dummy:
    return core

  core_nbr = core_dummy.GetNeighbors()[0].GetIdx()

  frag_dummy = None
  for atom in frag.GetAtoms():
    if atom.GetAtomicNum() == 0 and (
        atom.GetAtomMapNum() == r_map_num or atom.GetIsotope() == r_map_num
    ):
      frag_dummy = atom
      break
  if not frag_dummy:
    return core

  frag_nbr = frag_dummy.GetNeighbors()[0].GetIdx()

  combined = Chem.CombineMols(core, frag)
  rw = Chem.RWMol(combined)
  shift = core.GetNumAtoms()
  rw.AddBond(core_nbr, frag_nbr + shift, Chem.BondType.SINGLE)

  dummy_indices = []
  for atom in rw.GetAtoms():
    if atom.GetAtomicNum() == 0 and (
        atom.GetAtomMapNum() == r_map_num or atom.GetIsotope() == r_map_num
    ):
      dummy_indices.append(atom.GetIdx())
  for d_idx in sorted(dummy_indices, reverse=True):
    rw.RemoveAtom(d_idx)

  mol = rw.GetMol()
  Chem.SanitizeMol(mol)
  return mol


def _remove_dummy_as_h(core: Chem.Mol, r_map_num: int) -> Chem.Mol:
  """Replaces a dummy atom at r_map_num with hydrogen."""
  rw = Chem.RWMol(core)
  core_dummy_idx = None
  for atom in rw.GetAtoms():
    if atom.GetAtomicNum() == 0 and (
        atom.GetAtomMapNum() == r_map_num or atom.GetIsotope() == r_map_num
    ):
      core_dummy_idx = atom.GetIdx()
      break
  if core_dummy_idx is not None:
    rw.RemoveAtom(core_dummy_idx)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return mol
  return core


def assemble_markush_molecule(
    scaffold_mol: Chem.Mol,
    r_groups: Dict[int, str],
    compound_id: str = "CMPD",
) -> Optional[Chem.Mol]:
  """Recombines a generic core scaffold with numbered R-group attachment points."""
  try:
    mol = Chem.Mol(scaffold_mol)

    for r_idx, r_smi in r_groups.items():
      if not r_smi or r_smi == "[H]":
        mol = _remove_dummy_as_h(mol, r_idx)
        continue

      frag_mol = Chem.MolFromSmiles(r_smi)
      if frag_mol is None:
        continue

      mol = _attach_fragment(mol, frag_mol, r_idx)

    # Clean up any residual unmapped dummy atoms
    dummy_indices = [
        a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0
    ]
    if dummy_indices:
      rw = Chem.RWMol(mol)
      for d_idx in sorted(dummy_indices, reverse=True):
        rw.RemoveAtom(d_idx)
      mol = rw.GetMol()
      Chem.SanitizeMol(mol)

    mol.SetProp("_Name", compound_id)
    mol.SetProp("Compound_ID", compound_id)
    mol.SetProp("CanonicalSMILES", Chem.MolToSmiles(mol))
    return mol
  except Exception:
    return None


CARBON_RANGE_REGEX = re.compile(
    r"^\s*(?:C|c)\s*([0-9]+)\s*(?:-|to)\s*(?:C|c)?\s*([0-9]+)\s*(alkyl|alkoxy|alkenyl|haloalkyl)?\s*$",
    re.IGNORECASE,
)

HOMOLOGOUS_ALKYL = {
    1: ["Me"],
    2: ["Et"],
    3: ["nPr", "iPr"],
    4: ["nBu", "iBu", "sBu", "tBu"],
    5: ["CCCCC", "CCCC(C)", "CC(C)(C)C"],
    6: ["CCCCCC", "CCCC(C)C", "CCC(C)CC", "CC(C)CCC"],
}

HOMOLOGOUS_ALKOXY = {
    1: ["OMe"],
    2: ["OEt"],
    3: ["OnPr", "OiPr"],
    4: ["OnBu", "OiBu", "OsBu", "OtBu"],
    5: ["OCCCCC", "OCCCC(C)", "OCC(C)(C)C"],
    6: ["OCCCCCC", "OCCCC(C)C", "OCCC(C)CC", "OCC(C)CCC"],
}

HOMOLOGOUS_HALOGEN = ["F", "Cl", "Br", "I"]


def expand_homologous_token(token: Any) -> List[str]:
  """Expands condensed carbon-range homologous tokens (e.g.

  'C1-C3 alkyl') into discrete substituent lists.
  """
  if not token or pd.isna(token):
    return [str(token)] if token is not None else [""]
  val_str = str(token).strip()

  # Check halogen range
  if re.match(r"^\s*(?:halogen|halo)\s*$", val_str, re.IGNORECASE):
    return list(HOMOLOGOUS_HALOGEN)

  # Check carbon range e.g. C1-C3 alkyl, C1-C2 alkoxy
  m = CARBON_RANGE_REGEX.match(val_str)
  if not m:
    return [val_str]

  n_min, n_max = int(m.group(1)), int(m.group(2))
  cat = (m.group(3) or "alkyl").lower()
  table = HOMOLOGOUS_ALKOXY if "alkoxy" in cat else HOMOLOGOUS_ALKYL

  expanded = []
  for n in range(n_min, n_max + 1):
    expanded.extend(table.get(n, []))

  return expanded if expanded else [val_str]


def expand_table_homologues(
    df: pd.DataFrame, id_col: str, r_cols: Dict[int, str]
) -> pd.DataFrame:
  """Expands table rows containing homologous substituent ranges into discrete Cartesian product rows."""
  expanded_rows = []
  r_col_names = [r_cols[idx] for idx in sorted(r_cols.keys())]

  for _, row in df.iterrows():
    token_lists = [expand_homologous_token(row[c]) for c in r_col_names]
    if all(len(tl) == 1 for tl in token_lists):
      expanded_rows.append(dict(row))
    else:
      base_id = str(row[id_col])
      prod = list(itertools.product(*token_lists))
      for i, combo in enumerate(prod, 1):
        r = dict(row)
        r[id_col] = f"{base_id}_{i}"
        for c_name, val in zip(r_col_names, combo):
          r[c_name] = val
        expanded_rows.append(r)

  return pd.DataFrame(expanded_rows)


def assemble_patent_table(
    scaffold_smiles: str,
    table_csv: str,
    id_col: Optional[str] = None,
    output_csv: Optional[str] = None,
    output_sdf: Optional[str] = None,
    max_mw: Optional[float] = None,
    max_logp: Optional[float] = None,
    max_rotb: Optional[int] = None,
    max_hbd: Optional[int] = None,
    max_hba: Optional[int] = None,
) -> Tuple[pd.DataFrame, List[Chem.Mol]]:
  """Assembles an entire table of exemplified compounds from a patent CSV/TSV with optional property filters."""
  df = pd.read_csv(table_csv)

  cols = list(df.columns)
  if not id_col:
    cand_id = [
        c
        for c in cols
        if any(
            k in c.lower() for k in ["id", "compound", "example", "no", "name"]
        )
    ]
    id_col = cand_id[0] if cand_id else cols[0]

  r_cols = {}
  for c in cols:
    if c == id_col:
      continue
    m = re.search(r"r\s*_?\[?([0-9]+)\]?|\[\*:([0-9]+)\]", c.lower())
    if m:
      idx = int(m.group(1) or m.group(2))
      r_cols[idx] = c

  if not r_cols:
    idx = 1
    for c in cols:
      if c != id_col:
        r_cols[idx] = c
        idx += 1

  # Expand homologous series ranges (e.g. C1-C3 alkyl, C1-C2 alkoxy) into discrete rows
  df = expand_table_homologues(df, id_col, r_cols)

  scaff_mol = Chem.MolFromSmiles(scaffold_smiles)
  if scaff_mol is None:
    raise ValueError(f"Invalid scaffold SMILES: {scaffold_smiles}")

  for atom in scaff_mol.GetAtoms():
    if (
        atom.GetAtomicNum() == 0
        and atom.GetAtomMapNum() == 0
        and atom.GetIsotope() > 0
    ):
      atom.SetAtomMapNum(atom.GetIsotope())

  has_filtering = any(
      f is not None for f in [max_mw, max_logp, max_rotb, max_hbd, max_hba]
  )

  results = []
  mols = []

  for _, row in df.iterrows():
    cmpd_id = str(row[id_col])
    r_map = {}
    for r_idx, c_name in r_cols.items():
      raw_val = row[c_name]
      frag = normalize_rgroup_string(raw_val, r_idx)
      r_map[r_idx] = frag

    assembled_mol = assemble_markush_molecule(scaff_mol, r_map, cmpd_id)
    if assembled_mol is not None:
      smi = Chem.MolToSmiles(assembled_mol)
      res_row = dict(row)
      res_row["assembled_smiles"] = smi
      res_row["assembly_status"] = "SUCCESS"

      # Calculate physicochemical descriptors
      mw = round(Descriptors.MolWt(assembled_mol), 2)
      logp = round(Descriptors.MolLogP(assembled_mol), 2)
      rotb = rdMolDescriptors.CalcNumRotatableBonds(assembled_mol)
      hbd = rdMolDescriptors.CalcNumHBD(assembled_mol)
      hba = rdMolDescriptors.CalcNumHBA(assembled_mol)

      res_row["mw"] = mw
      res_row["logp"] = logp
      res_row["rotatable_bonds"] = rotb
      res_row["hbd"] = hbd
      res_row["hba"] = hba

      pass_filter = True
      if max_mw is not None and mw > max_mw:
        pass_filter = False
      if max_logp is not None and logp > max_logp:
        pass_filter = False
      if max_rotb is not None and rotb > max_rotb:
        pass_filter = False
      if max_hbd is not None and hbd > max_hbd:
        pass_filter = False
      if max_hba is not None and hba > max_hba:
        pass_filter = False

      res_row["filter_passed"] = pass_filter

      if pass_filter:
        mols.append(assembled_mol)
        results.append(res_row)
      elif not has_filtering:
        mols.append(assembled_mol)
        results.append(res_row)
    else:
      res_row = dict(row)
      res_row["assembled_smiles"] = ""
      res_row["assembly_status"] = "FAILED"
      if not has_filtering:
        results.append(res_row)

  out_df = pd.DataFrame(results)

  if output_csv:
    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    out_df.to_csv(output_csv, index=False)
    print(f"Saved {len(mols)} assembled molecules to CSV: {output_csv}")

  if output_sdf and mols:
    os.makedirs(os.path.dirname(os.path.abspath(output_sdf)), exist_ok=True)
    with Chem.SDWriter(output_sdf) as w:
      for m in mols:
        w.write(m)
    print(f"Saved {len(mols)} assembled molecules to SDF: {output_sdf}")

  return out_df, mols


def main():
  parser = argparse.ArgumentParser(
      description="Patent Markush Combinatorial Structure Assembler."
  )
  parser.add_argument(
      "--scaffold_smiles",
      "--scaffold-smiles",
      "--scaffold",
      "--core",
      type=str,
      dest="scaffold_smiles",
      required=True,
      help="Scaffold SMILES with [*:1], [*:2] wildcards.",
  )
  parser.add_argument(
      "--table_csv",
      "--table-csv",
      "--table",
      "--input-csv",
      "--input_csv",
      "--input",
      type=str,
      dest="table_csv",
      required=True,
      help="Path to patent table CSV/TSV.",
  )
  parser.add_argument(
      "--id_col",
      "--id-col",
      type=str,
      dest="id_col",
      help="Column name for compound identifier.",
  )
  parser.add_argument(
      "--output_csv",
      "--output-csv",
      type=str,
      dest="output_csv",
      help="Path to write assembled CSV.",
  )
  parser.add_argument(
      "--output_sdf",
      "--output-sdf",
      type=str,
      dest="output_sdf",
      help="Path to write assembled SDF.",
  )
  parser.add_argument(
      "--max-mw",
      "--max_mw",
      type=float,
      default=None,
      dest="max_mw",
      help="Maximum molecular weight filter.",
  )
  parser.add_argument(
      "--max-logp",
      "--max_logp",
      type=float,
      default=None,
      dest="max_logp",
      help="Maximum cLogP filter.",
  )
  parser.add_argument(
      "--max-rotb",
      "--max_rotb",
      type=int,
      default=None,
      dest="max_rotb",
      help="Maximum rotatable bonds filter.",
  )
  parser.add_argument(
      "--max-hbd",
      "--max_hbd",
      type=int,
      default=None,
      dest="max_hbd",
      help="Maximum H-bond donors filter.",
  )
  parser.add_argument(
      "--max-hba",
      "--max_hba",
      type=int,
      default=None,
      dest="max_hba",
      help="Maximum H-bond acceptors filter.",
  )
  args = parser.parse_args()

  out_df, mols = assemble_patent_table(
      scaffold_smiles=args.scaffold_smiles,
      table_csv=args.table_csv,
      id_col=args.id_col,
      output_csv=args.output_csv,
      output_sdf=args.output_sdf,
      max_mw=args.max_mw,
      max_logp=args.max_logp,
      max_rotb=args.max_rotb,
      max_hbd=args.max_hbd,
      max_hba=args.max_hba,
  )
  print(f"Assembly completed: {len(mols)}/{len(out_df)} succeeded.")


if __name__ == "__main__":
  main()
