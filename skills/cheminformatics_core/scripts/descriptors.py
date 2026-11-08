# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "pandas",
#     "numpy",
# ]
# ///

"""Standardized Molecular Descriptor Engine & Column Normalization Hub.

Eliminates RDKit API attribute traps (Lipinski vs Descriptors vs
rdMolDescriptors) and unifies column naming across chemical datasets and machine
learning workflows.
"""

import argparse
import collections
import json
import os
import re
import sys
from typing import Dict, Optional, Sequence, Union

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors

_script_dir = os.path.dirname(os.path.abspath(__file__))
if _script_dir not in sys.path:
  sys.path.insert(0, _script_dir)

from synthetic_accessibility_scorer import calculate_sa_score

# Ensure Lipinski module has NumHBD and NumHBA aliases to prevent AttributeError
if not hasattr(Lipinski, "NumHBD"):
  Lipinski.NumHBD = Lipinski.NumHDonors
if not hasattr(Lipinski, "NumHBA"):
  Lipinski.NumHBA = Lipinski.NumHAcceptors

# -----------------------------------------------------------------------------
# Canonical Column Schemas & Intelligent Normalization Engine
# -----------------------------------------------------------------------------
DESCRIPTOR_CANONICAL_SHORT = [
    "MW",
    "cLogP",
    "TPSA",
    "HBD",
    "HBA",
    "RotB",
    "Fsp3",
    "HeavyAtoms",
    "Rings",
    "AromaticRings",
    "AliphaticRings",
    "MaxRingSize",
    "is_macrocycle",
    "QED",
    "SA_Score",
    "CNS_MPO",
]

DESCRIPTOR_CANONICAL_SNAKE = [
    "molecular_weight",
    "clogp",
    "tpsa",
    "hbd",
    "hba",
    "rotatable_bonds",
    "fraction_csp3",
    "heavy_atom_count",
    "ring_count",
    "aromatic_ring_count",
    "aliphatic_ring_count",
    "max_ring_size",
    "is_macrocycle",
    "qed",
    "sa_score",
    "cns_mpo",
]


class DescriptorCanonicalizer:
  """Expansive descriptor column normalizer.

  Uses multi-tier tokenization, regex heuristics, unit-stripping, CamelCase
  splitting, and prefix/suffix cleansing to robustly map arbitrary real-world
  column variations (e.g. 'MW (Da)', 'PolarSurfaceArea', 'calc_tpsa',
  'number_of_rotatable_bonds', 'chembl_logp', 'ExactMolWt [g/mol]') into
  canonical short or snake-case identifiers. Preserves non-descriptor metadata
  (like compound IDs, yields, assay metadata) completely intact.
  """

  KNOWN_SHORT = {
      "rotb",
      "clogp",
      "tpsa",
      "psa",
      "hbd",
      "hba",
      "fsp3",
      "fcsp3",
      "qed",
      "mw",
      "molwt",
      "exactmw",
      "pic50",
      "ic50",
      "ki",
      "kd",
      "rdkit",
      "chembl",
      "pfizer",
      "lipinski",
      "smiles",
      "smi",
  }

  CANONICAL_TARGETS = {
      "MW": {
          "short": "MW",
          "snake": "molecular_weight",
          "keywords": ["mw", "molwt", "exactmw", "exactmolwt", "molweight"],
          "phrases": [
              ("molecular", "weight"),
              ("mol", "weight"),
              ("exact", "weight"),
          ],
          "regex": [
              r"\bmw\b",
              r"\bmol_?wt\b",
              r"\bmolecular_?weight\b",
              r"\bexact_?mw\b",
              r"\bexact_?mol_?wt\b",
          ],
      },
      "cLogP": {
          "short": "cLogP",
          "snake": "clogp",
          "keywords": [
              "clogp",
              "logp",
              "slogp",
              "alogp",
              "mollogp",
              "lipophilicity",
          ],
          "phrases": [("log", "p"), ("c", "log", "p")],
          "regex": [
              r"\bc?logp\b",
              r"\bslogp\b",
              r"\balogp\b",
              r"\bmol_?logp\b",
              r"\blipophilicity\b",
          ],
      },
      "TPSA": {
          "short": "TPSA",
          "snake": "tpsa",
          "keywords": ["tpsa", "psa"],
          "phrases": [("polar", "surface", "area"), ("polar", "surface")],
          "regex": [r"\btpsa\b", r"\bpsa\b", r"\bpolar_?surface_?area\b"],
      },
      "HBD": {
          "short": "HBD",
          "snake": "hbd",
          "keywords": ["hbd", "numhbd", "hdonors", "numhdonors", "hdonor"],
          "phrases": [
              ("h", "donor"),
              ("hydrogen", "donor"),
              ("hbond", "donor"),
              ("h", "donors"),
              ("hydrogen", "donors"),
          ],
          "regex": [
              r"\bhbd\b",
              r"\bnum_?h_?donors?\b",
              r"\bh_?donors?\b",
              r"\bhydrogen_?donors?\b",
          ],
      },
      "HBA": {
          "short": "HBA",
          "snake": "hba",
          "keywords": [
              "hba",
              "numhba",
              "hacceptors",
              "numhacceptors",
              "hacceptor",
          ],
          "phrases": [
              ("h", "acceptor"),
              ("hydrogen", "acceptor"),
              ("hbond", "acceptor"),
              ("h", "acceptors"),
              ("hydrogen", "acceptors"),
          ],
          "regex": [
              r"\bhba\b",
              r"\bnum_?h_?acceptors?\b",
              r"\bh_?acceptors?\b",
              r"\bhydrogen_?acceptors?\b",
          ],
      },
      "RotB": {
          "short": "RotB",
          "snake": "rotatable_bonds",
          "keywords": [
              "rotb",
              "nrotb",
              "numrotb",
              "rotatablebonds",
              "rotbonds",
          ],
          "phrases": [
              ("rotatable", "bond"),
              ("rotatable", "bonds"),
              ("rot", "bonds"),
          ],
          "regex": [
              r"\brotb\b",
              r"\bn_?rotb\b",
              r"\brotatable_?bonds?\b",
              r"\bnum_?rotatable_?bonds?\b",
          ],
      },
      "Fsp3": {
          "short": "Fsp3",
          "snake": "fraction_csp3",
          "keywords": [
              "fsp3",
              "fcsp3",
              "fractioncsp3",
              "fractions_sp3",
              "csp3",
          ],
          "phrases": [
              ("fraction", "csp3"),
              ("fraction", "sp3"),
              ("fraction", "carbon", "sp3"),
          ],
          "regex": [
              r"\bf_?c?sp3\b",
              r"\bfraction_?c?sp3\b",
              r"\bsp3_?fraction\b",
          ],
      },
      "HeavyAtoms": {
          "short": "HeavyAtoms",
          "snake": "heavy_atom_count",
          "keywords": [
              "heavyatoms",
              "heavyatomcount",
              "nheavy",
              "numheavyatoms",
              "heavyatom",
          ],
          "phrases": [
              ("heavy", "atom"),
              ("heavy", "atoms"),
              ("num", "heavy"),
              ("n", "heavy"),
          ],
          "regex": [
              r"\bheavy_?atoms?\b",
              r"\bheavy_?atom_?count\b",
              r"\bn_?heavy\b",
              r"\bnum_?heavy\b",
          ],
      },
      "AromaticRings": {
          "short": "AromaticRings",
          "snake": "aromatic_ring_count",
          "keywords": ["aromaticrings", "numaromaticrings"],
          "phrases": [("aromatic", "ring"), ("aromatic", "rings")],
          "regex": [
              r"\baromatic_?rings?\b",
              r"\barom_?rings?\b",
              r"\bnum_?aromatic_?rings?\b",
          ],
      },
      "AliphaticRings": {
          "short": "AliphaticRings",
          "snake": "aliphatic_ring_count",
          "keywords": ["aliphaticrings", "numaliphaticrings"],
          "phrases": [("aliphatic", "ring"), ("aliphatic", "rings")],
          "regex": [
              r"\baliphatic_?rings?\b",
              r"\baliph_?rings?\b",
              r"\bnum_?aliphatic_?rings?\b",
          ],
      },
      "Rings": {
          "short": "Rings",
          "snake": "ring_count",
          "keywords": ["rings", "ringcount", "numrings"],
          "phrases": [("ring", "count"), ("total", "rings")],
          "regex": [r"\brings?\b", r"\bring_?count\b", r"\bnum_?rings?\b"],
      },
      "QED": {
          "short": "QED",
          "snake": "qed",
          "keywords": ["qed", "druglikeness"],
          "phrases": [("drug", "likeness"), ("qed", "score")],
          "regex": [r"\bqed\b", r"\bdrug_?likeness\b"],
      },
      "SA_Score": {
          "short": "SA_Score",
          "snake": "sa_score",
          "keywords": ["sascore", "sa_score", "synthaccess", "sa"],
          "phrases": [("synthetic", "accessibility"), ("synth", "access")],
          "regex": [
              r"\bsa_?score\b",
              r"\bsa\b",
              r"\bsynthetic_?accessibility\b",
          ],
      },
      "pIC50": {
          "short": "pIC50",
          "snake": "pic50",
          "keywords": ["pic50", "p_ic50"],
          "phrases": [("p", "ic50"), ("neg", "log", "ic50")],
          "regex": [r"\bp_?ic50\b", r"\bpic50\b"],
      },
      "IC50": {
          "short": "IC50",
          "snake": "ic50",
          "keywords": ["ic50"],
          "phrases": [("ic50", "value"), ("ic50", "nm"), ("ic50", "um")],
          "regex": [r"\bic50\b"],
      },
      "Ki": {
          "short": "Ki",
          "snake": "ki",
          "keywords": ["ki"],
          "phrases": [("ki", "value"), ("ki", "nm"), ("ki", "um")],
          "regex": [r"\bki\b"],
      },
      "Kd": {
          "short": "Kd",
          "snake": "kd",
          "keywords": ["kd"],
          "phrases": [("kd", "value"), ("kd", "nm"), ("kd", "um")],
          "regex": [r"\bkd\b"],
      },
  }

  ABBREVIATIONS = {
      "wt": "weight",
      "arom": "aromatic",
      "aliph": "aliphatic",
      "rot": "rotatable",
      "acc": "acceptor",
      "don": "donor",
  }

  # Tokens that carry no meaning of their own. A column is only renamed when
  # every token outside the match is one of these, so `IC50 (nM)` resolves to
  # IC50 while `IC50_ALK` keeps the target that distinguishes it.
  FILLER_TOKENS = frozenset({
      "value",
      "val",
      "score",
      "count",
      "number",
      "num",
      "total",
      "mean",
      "median",
      "std",
      "nm",
      "um",
      "mm",
      "da",
      "g",
      "mol",
      "unit",
      "units",
      "calc",
      "calculated",
      "pred",
      "predicted",
      "prop",
      "property",
      "desc",
      "descriptor",
      "rdkit",
      "chembl",
      "pfizer",
      "lipinski",
      "of",
  })

  @classmethod
  def clean_name(cls, name: str) -> str:
    """Standardizes arbitrary raw column strings into space-separated normalized tokens."""
    s = str(name).strip()
    # 1. Strip brackets, parens, braces and their contents (units, qualifiers)
    s = re.sub(r"\s*[\(\[\{][^\)\]\}]*[\)\]\}]", "", s)
    # 2. Replace non-alphanumerics with space
    s = re.sub(r"[^a-zA-Z0-9]+", " ", s)

    # 3. Process each space-separated word, preserving acronyms while splitting CamelCase
    words = []
    for word in s.split():
      w_lower = word.lower()
      if w_lower in cls.KNOWN_SHORT:
        words.append(w_lower)
      else:
        split_w = re.sub(r"([a-z0-9])([A-Z])", r"\1 \2", word)
        split_w = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1 \2", split_w)
        for sub_w in split_w.split():
          words.append(sub_w.lower())

    # 4. Strip common tool/library prefixes
    while words and words[0] in (
        "rdkit",
        "calc",
        "calculated",
        "chembl",
        "pfizer",
        "pred",
        "predicted",
        "prop",
        "property",
        "desc",
        "descriptor",
        "lipinski",
    ):
      words.pop(0)

    # 5. Expand abbreviations
    words = [cls.ABBREVIATIONS.get(w, w) for w in words]

    # 6. Filter out common suffixes
    if len(words) > 1 and words[-1] in (
        "score",
        "value",
        "val",
        "da",
        "g_mol",
        "nm",
        "um",
        "count",
        "v1",
        "v2",
        "v3",
        "mean",
        "median",
        "std",
        "druglike",
        "str",
    ):
      words.pop(-1)

    # 7. Strip filler phrases
    res = " ".join(words)
    res = re.sub(
        r"\b(number\s+of|num\s+of|count\s+of|total\s+number\s+of)\b", "", res
    ).strip()
    return res

  @classmethod
  def _regex_covered_indices(cls, tokens, pattern):
    """Indices of the tokens a regex match spans, or None if it does not match.

    Args:
      tokens: Cleaned tokens of the column name, in order.
      pattern: One canonical-target regex, written against the underscore-
        joined form of the name.

    Returns:
      Set of token indices overlapped by the match, or None.
    """
    match = re.search(pattern, "_".join(tokens))
    if not match:
      return None
    start, end = match.span()
    covered = set()
    offset = 0
    for index, token in enumerate(tokens):
      if offset < end and offset + len(token) > start:
        covered.add(index)
      offset += len(token) + 1
    return covered

  @classmethod
  def _canonical_if_fully_explained(cls, info, tokens, covered, target_style):
    """Returns the canonical name only if every unmatched token is filler.

    A match on part of a name does not make the name that quantity.
    `IC50 (nM)` is an IC50; `IC50_ALK` is one assay among many and keeps its
    own name, because collapsing it would both lose the target and collide
    with every other IC50 column in the table.

    Args:
      info: The CANONICAL_TARGETS entry that matched.
      tokens: Cleaned tokens of the column name, in order.
      covered: Indices of the tokens the match accounts for.
      target_style: "short" or "snake".

    Returns:
      The canonical column name, or None to leave the column alone.
    """
    residual = {t for i, t in enumerate(tokens) if i not in covered}
    if residual - cls.FILLER_TOKENS:
      return None
    return info["snake"] if target_style == "snake" else info["short"]

  @classmethod
  def resolve(cls, col_name: str, target_style: str = "short") -> str:
    """Resolves an arbitrary column name to its canonical identifier.

    If the column does not match any known descriptor/property, or matches
    only part of itself, returns the original column name.
    """
    target_style = target_style.lower()
    cleaned = cls.clean_name(col_name)
    tokens = cleaned.split()
    if not tokens:
      return col_name
    token_set = set(tokens)
    cleaned_no_space = cleaned.replace(" ", "")

    # Priority check: AromaticRings vs AliphaticRings vs general Rings
    for key, markers in (
        ("AromaticRings", {"aromatic", "arom"}),
        ("AliphaticRings", {"aliphatic", "aliph"}),
    ):
      if token_set & markers and token_set & {"ring", "rings"}:
        covered = {
            i for i, t in enumerate(tokens) if t in markers | {"ring", "rings"}
        }
        canonical = cls._canonical_if_fully_explained(
            cls.CANONICAL_TARGETS[key], tokens, covered, target_style
        )
        if canonical:
          return canonical

    # Priority check: pIC50 vs IC50
    if (
        "pic50" in cleaned_no_space
        or {"p", "ic50"} <= token_set
        or {"neg", "ic50"} <= token_set
    ):
      covered = {
          i
          for i, t in enumerate(tokens)
          if t in {"pic50", "p", "neg", "log", "ic50"}
      }
      canonical = cls._canonical_if_fully_explained(
          cls.CANONICAL_TARGETS["pIC50"], tokens, covered, target_style
      )
      if canonical:
        return canonical

    # Target rule checks
    for key, info in cls.CANONICAL_TARGETS.items():
      if key in ("AromaticRings", "AliphaticRings", "pIC50"):
        continue
      if cleaned_no_space in info["keywords"]:
        return info["snake"] if target_style == "snake" else info["short"]
      keyword_hits = {i for i, t in enumerate(tokens) if t in info["keywords"]}
      if keyword_hits:
        canonical = cls._canonical_if_fully_explained(
            info, tokens, keyword_hits, target_style
        )
        if canonical:
          return canonical
      for phrase in info["phrases"]:
        if all(p in token_set for p in phrase):
          covered = {i for i, t in enumerate(tokens) if t in phrase}
          canonical = cls._canonical_if_fully_explained(
              info, tokens, covered, target_style
          )
          if canonical:
            return canonical
      for pattern in info["regex"]:
        covered = cls._regex_covered_indices(tokens, pattern)
        if covered:
          canonical = cls._canonical_if_fully_explained(
              info, tokens, covered, target_style
          )
          if canonical:
            return canonical

    return col_name

  @classmethod
  def get_column_mapping(
      cls, columns: Sequence[str], target_style: str = "short"
  ) -> Dict[str, str]:
    """Returns a renaming dictionary mapping raw column names to canonical names."""
    mapping = {}
    for col in columns:
      canon = cls.resolve(col, target_style=target_style)
      if canon != col:
        mapping[col] = canon
    return mapping

  @classmethod
  def normalize_dataframe(
      cls, df: pd.DataFrame, target_style: str = "short"
  ) -> pd.DataFrame:
    """Renames recognized descriptor columns, never creating duplicate labels.

    Two columns resolving to one canonical name would leave the frame with a
    repeated label, which breaks `pd.concat` and silently makes label-based
    indexing ambiguous. When that would happen, every column in the clashing
    group keeps its original name.

    Args:
      df: Input DataFrame.
      target_style: "short" or "snake".

    Returns:
      A DataFrame with unique column labels.
    """
    mapping = cls.get_column_mapping(df.columns, target_style=target_style)
    if not mapping:
      return df.copy()
    claimed = collections.Counter(mapping.values())
    untouched = {c for c in df.columns if c not in mapping}
    safe = {
        source: canonical
        for source, canonical in mapping.items()
        if claimed[canonical] == 1 and canonical not in untouched
    }
    if not safe:
      return df.copy()
    return df.rename(columns=safe)


# Backward compatibility dictionary
DESCRIPTOR_ALIAS_MAP = {
    k: (v["short"], v["snake"])
    for k, v in DescriptorCanonicalizer.CANONICAL_TARGETS.items()
}


def calculate_cns_mpo(
    mol: "Chem.Mol",
    clogp: Optional[float] = None,
    mw: Optional[float] = None,
    tpsa: Optional[float] = None,
    hbd: Optional[int] = None,
    pka: Optional[float] = None,
) -> float:
  """Canonical Pfizer CNS Multiparameter Optimization (CNS MPO) score [0.0 - 6.0].

  This is the single authoritative CNS MPO implementation for the RDKit skill
  suite. Ref: Wager et al., ACS Chem. Neurosci. 2010, 1, 435-449; 2016, 7,
  767-775.

  Six desirability components, each scored 0.0-1.0 and summed:
    1. cLogP  (<=3.0 -> 1.0, >=5.0 -> 0.0)
    2. cLogD  (<=2.0 -> 1.0, >=4.0 -> 0.0)
    3. MW     (<=360 -> 1.0, >=500 -> 0.0)
    4. TPSA   (40-90 -> 1.0; 20-40 and 90-120 ramps; outside -> 0.0)
    5. HBD    (<=0.5 -> 1.0, >3.5 -> 0.0)
    6. pKa    (<=8.0 -> 1.0, >=10.0 -> 0.0)

  Scores >= 4.0 indicate favorable alignment for blood-brain barrier
  penetration.

  Args:
    mol: Input molecule.
    clogp: Optional precomputed cLogP; derived via Crippen.MolLogP if omitted.
    mw: Optional precomputed exact molecular weight.
    tpsa: Optional precomputed topological polar surface area.
    hbd: Optional precomputed hydrogen bond donor count.
    pka: Optional measured/predicted most-basic pKa. When omitted, estimated as
      9.5 if a basic amine is present and 5.0 otherwise.

  Returns:
    CNS MPO score in [0.0, 6.0], rounded to 2 decimal places.
  """
  if mol is None:
    return 0.0

  if mw is None:
    mw = float(Descriptors.MolWt(mol))
  if clogp is None:
    clogp = float(Crippen.MolLogP(mol))
  if tpsa is None:
    tpsa = float(rdMolDescriptors.CalcTPSA(mol))
  if hbd is None:
    hbd = int(Lipinski.NumHDonors(mol))

  # Heuristic most-basic pKa when no measured value is supplied.
  if pka is None:
    basic_amine = Chem.MolFromSmarts(
        "[NX3;H2,H1,H0;!$(NC=[O,S]);!$(NS=O);!$(Nc);!$(N#*)]"
    )
    has_basic_amine = (
        mol.HasSubstructMatch(basic_amine) if basic_amine else False
    )
    pka = 9.5 if has_basic_amine else 5.0

  # cLogD at physiological pH from the Henderson-Hasselbalch shift for bases.
  clogd = clogp - (pka - 7.4) if pka > 7.4 else clogp

  def _rescale(x: float, min_x: float, max_x: float) -> float:
    return min(1.0, max(0.0, (max_x - x) / (max_x - min_x)))

  s_clogp = _rescale(clogp, 3.0, 5.0)
  s_clogd = _rescale(clogd, 2.0, 4.0)
  s_mw = _rescale(mw, 360.0, 500.0)

  if 40.0 <= tpsa <= 90.0:
    s_tpsa = 1.0
  elif 20.0 < tpsa < 40.0:
    s_tpsa = (tpsa - 20.0) / 20.0
  elif 90.0 < tpsa < 120.0:
    s_tpsa = _rescale(tpsa, 90.0, 120.0)
  else:
    s_tpsa = 0.0

  s_hbd = _rescale(float(hbd), 0.5, 3.5)
  s_pka = _rescale(pka, 8.0, 10.0)

  return round(float(s_clogp + s_clogd + s_mw + s_tpsa + s_hbd + s_pka), 2)


def get_macrocycle_info(
    mol: "Chem.Mol", min_ring_size: int = 12
) -> tuple[bool, int]:
  """Detects whether a molecule is macrocyclic and returns (is_macrocycle, max_ring_size)."""
  if mol is None:
    return False, 0
  ring_info = mol.GetRingInfo()
  atom_rings = ring_info.AtomRings()
  if not atom_rings:
    return False, 0
  max_size = max(len(r) for r in atom_rings)
  return max_size >= min_ring_size, int(max_size)


def compute_standard_descriptors(
    mol: Chem.Mol,
    style: str = "short",
    include_sa_score: bool = True,
) -> Dict[str, Union[float, int, bool]]:
  """Computes physicochemical properties using verified, canonical RDKit APIs.

  Guarantees correct function calls:
  - HBD: Lipinski.NumHDonors(mol) (NEVER Lipinski.NumHBD)
  - HBA: Lipinski.NumHAcceptors(mol) (NEVER Lipinski.NumHBA)
  - RotB: Lipinski.NumRotatableBonds(mol)
  - MW: Descriptors.MolWt(mol)
  - cLogP: Crippen.MolLogP(mol)
  - TPSA: rdMolDescriptors.CalcTPSA(mol)
  - Fsp3: rdMolDescriptors.CalcFractionCSP3(mol)

  Args:
    mol: Valid RDKit Mol object.
    style: "short" (e.g. MW, cLogP, RotB) or "snake" (e.g. molecular_weight,
      clogp, rotatable_bonds).
    include_sa_score: If True, computes synthetic accessibility score.

  Returns:
    Dictionary of property name to numeric/boolean value.
  """
  if mol is None:
    raise ValueError("Cannot compute descriptors on None molecule.")

  mw = float(Descriptors.MolWt(mol))
  logp = float(Crippen.MolLogP(mol))
  tpsa = float(rdMolDescriptors.CalcTPSA(mol))
  hbd = int(Lipinski.NumHDonors(mol))
  hba = int(Lipinski.NumHAcceptors(mol))
  rotb = int(Lipinski.NumRotatableBonds(mol))
  fsp3 = float(rdMolDescriptors.CalcFractionCSP3(mol))
  heavy_atoms = int(mol.GetNumHeavyAtoms())
  rings = int(Lipinski.RingCount(mol))
  aromatic_rings = int(rdMolDescriptors.CalcNumAromaticRings(mol))
  aliphatic_rings = int(rdMolDescriptors.CalcNumAliphaticRings(mol))
  is_macro, max_ring_size = get_macrocycle_info(mol)
  qed_score = float(QED.qed(mol))
  cns_mpo = calculate_cns_mpo(mol, clogp=logp, mw=mw, tpsa=tpsa, hbd=hbd)

  sa_val: Optional[float] = None
  if include_sa_score:
    try:
      sa_val = float(calculate_sa_score(mol))
    except Exception as exc:  # pragma: no cover - defensive
      print(
          f"Warning: SA_Score could not be computed ({type(exc).__name__}:"
          f" {exc}).",
          file=sys.stderr,
      )
      sa_val = None

  if style.lower() == "snake":
    d: Dict[str, Union[float, int, bool]] = {
        "molecular_weight": round(mw, 3),
        "clogp": round(logp, 2),
        "tpsa": round(tpsa, 1),
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rotb,
        "fraction_csp3": round(fsp3, 3),
        "heavy_atom_count": heavy_atoms,
        "ring_count": rings,
        "aromatic_ring_count": aromatic_rings,
        "aliphatic_ring_count": aliphatic_rings,
        "max_ring_size": max_ring_size,
        "is_macrocycle": is_macro,
        "qed": round(qed_score, 4),
    }
    if sa_val is not None:
      d["sa_score"] = round(sa_val, 2)
    d["cns_mpo"] = cns_mpo
  else:
    d = {
        "MW": round(mw, 3),
        "cLogP": round(logp, 2),
        "TPSA": round(tpsa, 1),
        "HBD": hbd,
        "HBA": hba,
        "RotB": rotb,
        "Fsp3": round(fsp3, 3),
        "HeavyAtoms": heavy_atoms,
        "Rings": rings,
        "AromaticRings": aromatic_rings,
        "AliphaticRings": aliphatic_rings,
        "MaxRingSize": max_ring_size,
        "is_macrocycle": is_macro,
        "QED": round(qed_score, 4),
    }
    if sa_val is not None:
      d["SA_Score"] = round(sa_val, 2)
    d["CNS_MPO"] = cns_mpo

  return d


def normalize_descriptor_columns(
    df: pd.DataFrame,
    target_style: str = "short",
) -> pd.DataFrame:
  """Normalizes heterogeneous descriptor column names in a DataFrame to canonical schema.

  Resolves naming inconsistencies across arbitrary real-world formats using
  DescriptorCanonicalizer:
  - RotB, num_rotatable_bonds, rotatable_bonds, n_rotb -> RotB
  - MW (Da), molecular_weight [g/mol], ExactMolWt -> MW
  - cLogP, logp, chembl_logp, lipophilicity -> cLogP
  - TPSA (Å²), PolarSurfaceArea, calc_tpsa -> TPSA

  Args:
    df: Input pandas DataFrame.
    target_style: "short" (MW, cLogP, RotB) or "snake" (molecular_weight, clogp,
      rotatable_bonds).

  Returns:
    DataFrame with renamed canonical columns.
  """
  return DescriptorCanonicalizer.normalize_dataframe(
      df, target_style=target_style
  )


def featurize_smiles_list(
    smiles_list: Sequence[str],
    style: str = "short",
) -> pd.DataFrame:
  """Featurizes a list of SMILES strings into a standardized descriptor DataFrame."""
  records = []
  for idx, smi in enumerate(smiles_list):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
      records.append({"smiles": smi, "valid": False})
      continue
    props = compute_standard_descriptors(mol, style=style)
    props["smiles"] = smi
    props["valid"] = True
    records.append(props)
  return pd.DataFrame(records)


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Standardized Molecular Descriptor Engine & Column Normalization Hub."
      )
  )
  parser.add_argument(
      "--smiles",
      nargs="*",
      default=None,
      help="One or more SMILES strings to calculate.",
  )
  parser.add_argument(
      "--input",
      "--input-csv",
      dest="input_csv",
      default=None,
      help="Input CSV to featurize or normalize.",
  )
  parser.add_argument(
      "--output",
      "--output-csv",
      dest="output_csv",
      default=None,
      help="Output CSV path.",
  )
  parser.add_argument(
      "--smiles-col",
      default="smiles",
      help="Column name containing SMILES in input CSV.",
  )
  parser.add_argument(
      "--normalize-only",
      action="store_true",
      default=False,
      help=(
          "Only rename heterogeneous descriptor columns without recalculating."
      ),
  )
  parser.add_argument(
      "--style",
      choices=["short", "snake"],
      default="short",
      help=(
          "Target column naming style: 'short' (RotB, MW) or 'snake'"
          " (rotatable_bonds, molecular_weight)."
      ),
  )
  parser.add_argument(
      "--json",
      action="store_true",
      default=False,
      help="Output JSON for single SMILES calculation.",
  )

  args = parser.parse_args()

  if args.smiles:
    if len(args.smiles) == 1 and args.json:
      m = Chem.MolFromSmiles(args.smiles[0])
      if m is None:
        print(json.dumps({"error": f"Invalid SMILES: {args.smiles[0]}"}))
        sys.exit(1)
      print(
          json.dumps(
              compute_standard_descriptors(m, style=args.style), indent=2
          )
      )
      sys.exit(0)
    else:
      df = featurize_smiles_list(args.smiles, style=args.style)
      if args.output_csv:
        df.to_csv(args.output_csv, index=False)
        print(f"Saved featurized molecules to {args.output_csv}")
      else:
        print(df.to_string(index=False))
      sys.exit(0)

  if args.input_csv:
    df = pd.read_csv(args.input_csv)
    if args.normalize_only:
      df_norm = normalize_descriptor_columns(df, target_style=args.style)
      out_path = args.output_csv or args.input_csv
      df_norm.to_csv(out_path, index=False)
      print(f"Normalized descriptor columns in {out_path}")
    else:
      # Normalize existing columns first
      df_norm = normalize_descriptor_columns(df, target_style=args.style)
      # Find smiles column
      smi_col = None
      for c in df_norm.columns:
        if c.lower() in ("smiles", "smi", "canonical_smiles"):
          smi_col = c
          break
      if not smi_col:
        smi_col = args.smiles_col

      if smi_col not in df_norm.columns:
        print(
            f"Error: Could not find SMILES column '{smi_col}' in"
            f" {args.input_csv}",
            file=sys.stderr,
        )
        sys.exit(1)

      mols = [
          Chem.MolFromSmiles(s) if isinstance(s, str) else None
          for s in df_norm[smi_col]
      ]
      prop_records = []
      for m in mols:
        if m is not None:
          prop_records.append(compute_standard_descriptors(m, style=args.style))
        else:
          prop_records.append({})
      prop_df = pd.DataFrame(prop_records)

      for col in prop_df.columns:
        df_norm[col] = prop_df[col]

      out_path = args.output_csv or args.input_csv
      df_norm.to_csv(out_path, index=False)
      print(f"Featurized and normalized {len(df_norm)} molecules to {out_path}")
    sys.exit(0)

  parser.error("Must provide --smiles, --input_csv, or --print_schema.")


if __name__ == "__main__":
  main()
