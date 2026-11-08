# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pandas",
#     "polite-http",
#     "rdkit",
# ]
# ///

"""Synthetic Experimental Procedure Chemical Extractor.

Parses organic chemistry synthesis paragraphs from papers and patent
specifications, extracts target products, starting materials, reagents,
solvents, and reaction yields, and reconstructs validated RDKit Mol objects.
"""

import argparse
import html
import json
import os
import re
import sys
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

ADDUCT_MASS_SHIFTS = {
    "[m+h]+": 1.007276,
    "[m-h]-": -1.007276,
    "[m+na]+": 22.989218,
    "[m+k]+": 38.963158,
    "[m+nh4]+": 18.033823,
    "[m+h]": 1.007276,
    "[m-h]": -1.007276,
    "m+h": 1.007276,
    "m-h": -1.007276,
    "m+na": 22.989218,
}

# Import IUPACResolver
_script_dir = os.path.dirname(os.path.abspath(__file__))
if _script_dir not in sys.path:
  sys.path.insert(0, _script_dir)
from iupac_resolver import IUPACResolver
from chemical_text_sanitizer import clean_chemical_text

SALT_PREFIX_MULTIPLIERS = [
    (re.compile(r"^(?:hemi|semi)[-\s]?", re.IGNORECASE), 0.5),
    (re.compile(r"^sesqui[-\s]?", re.IGNORECASE), 1.5),
    (re.compile(r"^(?:di|bis)[-\s]?", re.IGNORECASE), 2.0),
    (re.compile(r"^(?:tri|tris)[-\s]?", re.IGNORECASE), 3.0),
    (re.compile(r"^(?:tetra|tetrakis)[-\s]?", re.IGNORECASE), 4.0),
    (re.compile(r"^mono[-\s]?", re.IGNORECASE), 1.0),
]

KNOWN_SALTS = [
    "trifluoroacetate",
    "trifluoroacetic acid",
    "tfa",
    "hydrochloride",
    "chloride",
    "dihydrochloride",
    "trihydrochloride",
    "hydrobromide",
    "bromide",
    "sulfate",
    "sulphate",
    "hydrogensulfate",
    "bisulfate",
    "hemisulfate",
    "mesylate",
    "methanesulfonate",
    "dimesylate",
    "tosylate",
    "p-toluenesulfonate",
    "acetate",
    "diacetate",
    "citrate",
    "tartrate",
    "maleate",
    "fumarate",
    "oxalate",
    "succinate",
    "benzoate",
    "phosphate",
    "nitrate",
    "formate",
]

KNOWN_HYDRATES = [
    "hydrate",
    "hemihydrate",
    "sesquihydrate",
    "dihydrate",
    "trihydrate",
    "tetrahydrate",
    "solvate",
    "ethanolate",
    "methanolate",
]

KNOWN_SOLVATES = [
    "ethanol solvate",
    "methanol solvate",
    "isopropanol solvate",
    "acetone solvate",
    "ether solvate",
    "ethyl acetate solvate",
    "solvate",
    "ethanolate",
    "methanolate",
    "isopropanolate",
]


_SUBSCRIPT_DIGITS = str.maketrans(
    "₀₁₂₃₄₅₆₇₈₉⁰¹²³⁴⁵⁶⁷⁸⁹",
    "01234567890123456789",
)

_INLINE_CHEM_TAG_RE = re.compile(
    r"</?(?:i|b|sub|sup|em|strong|span)\b[^>]*>",
    re.IGNORECASE,
)

# Strips NMR characterization headers/parentheticals (e.g. "1H NMR (400 MHz,
# DMSO-d6)", "13C{1H} NMR (101 MHz, CDCl3)") and standalone deuterated NMR
# solvent tokens (e.g. "DMSO-d6", "[D6]DMSO", "THF-d8") before solvent lookup.
_NMR_PAREN_BLOCK_RE = re.compile(
    r"\b(?:[0-9]+\s*[A-Za-z]+(?:\{[0-9A-Za-z]+\})?\s*-?\s*)?NMR\b\s*\([^)]*\)",
    re.IGNORECASE,
)
_DEUTERATED_SOLVENT_TOKEN_RE = re.compile(
    r"(?:\b[dD][0-9]+-|\[\s*[dD][0-9]+\s*\]-?)[A-Za-z0-9,-]+\b"
    r"|\b[A-Za-z0-9,]+-[dD][0-9]+\b"
)

_FORMULA_SOLVENT_DISPLAY_BY_TOKEN = {
    "DMF": "DMF",
    "DMSO": "DMSO",
    "THF": "THF",
    "DCM": "DCM",
    "CH2Cl2": "CH2Cl2",
    "CHCl3": "CHCl3",
    "CCl4": "CCl4",
    "EtOAc": "EtOAc",
    "ETOAC": "EtOAc",
    "EtAc": "ETAC",
    "ETAC": "ETAC",
    "EA": "EA",
    "MeOH": "MeOH",
    "MEOH": "MeOH",
    "EtOH": "EtOH",
    "ETOH": "EtOH",
    "i-PrOH": "i-PrOH",
    "iPrOH": "i-PrOH",
    "IPA": "IPA",
    "t-BuOH": "t-BuOH",
    "tBuOH": "t-BuOH",
    "H2O": "H2O",
    "PhMe": "PhMe",
    "CH3CN": "CH3CN",
    "MeCN": "MeCN",
    "MECN": "MeCN",
    "ACN": "ACN",
    "Et2O": "Et2O",
    "MTBE": "MTBE",
    "2-MeTHF": "2-MeTHF",
    "1,2-DCE": "1,2-DCE",
    "DCE": "DCE",
    "DMAC": "DMAC",
    "DMAc": "DMAC",
    "NMP": "NMP",
}

_NAME_SOLVENT_DISPLAY_BY_TOKEN = {
    "dichloromethane": "Dichloromethane",
    "chloroform": "Chloroform",
    "ethyl acetate": "Ethyl Acetate",
    "methanol": "Methanol",
    "ethanol": "Ethanol",
    "isopropanol": "Isopropanol",
    "water": "Water",
    "toluene": "Toluene",
    "acetonitrile": "Acetonitrile",
    "acetone": "Acetone",
    "diethyl ether": "Diethyl Ether",
    "ether": "Ether",
    "hexane": "Hexane",
    "hexanes": "Hexanes",
    "heptane": "Heptane",
    "1,4-dioxane": "1,4-Dioxane",
    "dioxane": "Dioxane",
}

_SOLVENT_DISPLAY_BY_TOKEN = {
    **_FORMULA_SOLVENT_DISPLAY_BY_TOKEN,
    **_NAME_SOLVENT_DISPLAY_BY_TOKEN,
}

SOLVENTS = tuple(_SOLVENT_DISPLAY_BY_TOKEN)

_FORMULA_SOLVENT_RE = re.compile(
    r"\b(?:"
    + "|".join(
        re.escape(tok)
        for tok in sorted(
            _FORMULA_SOLVENT_DISPLAY_BY_TOKEN, key=len, reverse=True
        )
    )
    + r")\b"
)

_NAME_SOLVENT_PATTERNS = tuple(
    (re.compile(r"\b" + re.escape(tok) + r"\b", re.IGNORECASE), display)
    for tok, display in _NAME_SOLVENT_DISPLAY_BY_TOKEN.items()
)


class ProcedureExtractor:
  """Rule-based and heuristic chemical entity extractor for synthetic procedures."""

  def __init__(self):
    self.resolver = IUPACResolver()

    # Known common solvents
    self.solvents = list(SOLVENTS)

    # Known common reagents/catalysts
    self.reagents = [
        "edci",
        "edc",
        "hobt",
        "hatu",
        "hbtu",
        "dcc",
        "cdi",
        "pd(pph3)4",
        "pd(dppf)cl2",
        "pd/c",
        "pd(oac)2",
        "k2co3",
        "na2co3",
        "cs2co3",
        "nahco3",
        "naoh",
        "koh",
        "et3n",
        "tea",
        "dipea",
        "hunig's base",
        "pyridine",
        "dmap",
        "tfa",
        "hcl",
        "nacl",
        "mgso4",
        "na2so4",
        "cucl",
        "cui",
        "boc2o",
    ]

  def deformulate_salt_and_hydrate(
      self, name: str, full_text: str = ""
  ) -> Tuple[
      str, Optional[str], Optional[float], Optional[str], Optional[float]
  ]:
    """Separates salt counterions and hydrate solvates from product name and extracts stoichiometric multipliers."""
    detected_salt = None
    salt_mult = None
    detected_hydrate = None
    hydrate_mult = None
    cleaned_name = name

    # 0. Look for solvates in name and strip them
    for solv in sorted(KNOWN_SOLVATES, key=len, reverse=True):
      pattern = (
          r"[\s,-]+((?:hemi|semi|sesqui|di|bis|tri|tris|tetra|tetrakis|mono)?[-"
          rf" ]?{solv})\b"
      )
      m = re.search(pattern, cleaned_name, re.IGNORECASE)
      if m:
        cleaned_name = re.sub(
            pattern, "", cleaned_name, count=1, flags=re.IGNORECASE
        ).strip(" -.,:")

    # 1. Look for hydrate token in name
    for hyd in sorted(KNOWN_HYDRATES, key=len, reverse=True):
      pattern = (
          r"[\s,-]+((?:hemi|semi|sesqui|di|bis|tri|tris|tetra|tetrakis|mono)?[-"
          rf" ]?{hyd})\b"
      )
      m = re.search(pattern, cleaned_name, re.IGNORECASE)
      if m:
        full_hyd = m.group(1).strip("- ")
        detected_hydrate = full_hyd
        hydrate_mult = 1.0
        for pref_re, mult in SALT_PREFIX_MULTIPLIERS:
          if pref_re.search(full_hyd):
            hydrate_mult = mult
            break
        cleaned_name = re.sub(
            pattern, "", cleaned_name, count=1, flags=re.IGNORECASE
        ).strip(" -.,:")
        break

    # 2. Look for salt tokens in name (handles multi-salt mixtures e.g. dimaleate hemisulfate)
    detected_salts = []
    for s in sorted(KNOWN_SALTS, key=len, reverse=True):
      pattern = (
          r"[\s,-]+((?:hemi|semi|sesqui|di|bis|tri|tris|tetra|tetrakis|mono)?[-"
          rf" ]?{s}(?:\s+salt)?)\b"
      )
      for match in re.finditer(pattern, cleaned_name, re.IGNORECASE):
        full_salt = match.group(1).strip("- ")
        curr_mult = 1.0
        for pref_re, mult in SALT_PREFIX_MULTIPLIERS:
          if pref_re.search(full_salt):
            curr_mult = mult
            break
        detected_salts.append(
            (match.start(), match.end(), full_salt, curr_mult, pattern)
        )

    if detected_salts:
      # Sort by start position in the name so first mentioned salt is primary
      detected_salts.sort(key=lambda x: x[0])
      primary_salt = detected_salts[0]
      detected_salt = primary_salt[2]
      salt_mult = primary_salt[3]

      # Strip all matched salt patterns from cleaned_name
      for _, _, _, _, pat in detected_salts:
        cleaned_name = re.sub(pat, "", cleaned_name, flags=re.IGNORECASE).strip(
            " -.,:"
        )

    # 3. Fallback: inspect full_text if salt/hydrate was not present in the name
    if not detected_salt and full_text:
      m_eq = re.search(
          r"treated\s+with\s+([0-9\.]+)\s+equivalents?\s+of\s+([A-Za-z0-9\s]+(?:acid|salt))",
          full_text,
          re.IGNORECASE,
      )
      if m_eq:
        salt_mult = float(m_eq.group(1))
        acid_phrase = m_eq.group(2).strip()
        if "sulfur" in acid_phrase.lower():
          detected_salt = "hemisulfate" if salt_mult == 0.5 else "sulfate"
        elif "trifluoro" in acid_phrase.lower():
          detected_salt = (
              "bis-trifluoroacetate" if salt_mult == 2.0 else "trifluoroacetate"
          )
        elif (
            "methanesulfon" in acid_phrase.lower()
            or "mesyl" in acid_phrase.lower()
        ):
          detected_salt = "dimesylate" if salt_mult == 2.0 else "mesylate"
        else:
          detected_salt = acid_phrase

    return (
        cleaned_name,
        detected_salt,
        salt_mult,
        detected_hydrate,
        hydrate_mult,
    )

  def isolate_free_base(self, mol: Chem.Mol) -> Tuple[Optional[Chem.Mol], str]:
    """Extracts the largest organic parent API fragment from a multi-component salt/solvate molecule."""
    if mol is None:
      return None, ""
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) == 1:
      return mol, Chem.MolToSmiles(mol)

    frags_sorted = sorted(
        frags, key=lambda f: f.GetNumHeavyAtoms(), reverse=True
    )
    best_frag = frags_sorted[0]
    return best_frag, Chem.MolToSmiles(best_frag)

  def extract_target_product_name(self, text: str) -> Optional[str]:
    """Extracts target product chemical name from title or yield sentence."""
    # Pattern 1: Title pattern: "Preparation of <Compound Name>:" or "Synthesis and chiral resolution of <Compound Name>:"
    m = re.search(
        r"(?:preparation|synthesis)(?:\s+and\s+[A-Za-z\s]+?)?\s+of\s+((?:[A-Za-z0-9\-\,\(\)\[\]\s]|\.(?!\s|$))+?)(?::|\.(?:\s|$)|\n|To"
        r" a\b)",
        text,
        re.IGNORECASE,
    )
    if m:
      cand = m.group(1).strip().rstrip(".:")
      if (
          len(cand) > 3
          and not cand.lower().startswith("the ")
          and not re.search(r"\bet\s+al\b", cand, re.IGNORECASE)
      ):
        return cand

    # Pattern 2: Yield pattern: "to afford/yield/give <Compound Name> (<quantity>, <yield>%)"
    m_yield = re.search(
        r"(?:afford|yield|give|provide)\s+(?:the\s+title\s+compound,?\s*)?((?:[A-Za-z0-9\-\,\(\)\[\]\s]|\.(?!\s|$))+?)\s*\([0-9\.]+\s*(?:mg|g|mmol)",
        text,
        re.IGNORECASE,
    )
    if m_yield:
      cand = m_yield.group(1).strip().rstrip(".:")
      if len(cand) > 3 and not re.search(r"\bet\s+al\b", cand, re.IGNORECASE):
        return cand

    return None

  def extract_yield(self, text: str) -> Dict[str, str]:
    """Extracts yield percentage and mass/mole quantities."""
    yield_info = {"yield_percent": "", "yield_mass": "", "yield_moles": ""}

    # Percent: prioritize explicit yield statements over reagent concentrations (e.g. 15% NaOH)
    m_pct = re.findall(
        r"([0-9]+(?:\.[0-9]+)?)\s*%\s*(?:yield|over\b)", text, re.IGNORECASE
    )
    if m_pct:
      yield_info["yield_percent"] = m_pct[-1]
    else:
      m_paren = re.findall(
          r"\([^)]*?([0-9]+(?:\.[0-9]+)?)\s*%\s*(?:yield)?\s*\)",
          text,
          re.IGNORECASE,
      )
      if m_paren:
        yield_info["yield_percent"] = m_paren[-1]
      else:
        m_word = re.findall(
            r"yield(?:\s*of|\s*:|\s*is)?\s*([0-9]+(?:\.[0-9]+)?)\s*%",
            text,
            re.IGNORECASE,
        )
        if m_word:
          yield_info["yield_percent"] = m_word[-1]
        else:
          m_all = re.findall(
              r"([0-9]+(?:\.[0-9]+)?)\s*%(?!\s*(?:naoh|koh|hcl|h2so4|aq|aqueous|wt|w/v|v/v|solution|soln|purity|ee|de)\b)",
              text,
              re.IGNORECASE,
          )
          if m_all:
            yield_info["yield_percent"] = m_all[-1]

    # Mass: check concluding parenthetical mass (isolated product) or last mass
    m_pmass = re.findall(
        r"\(\s*([0-9]+(?:\.[0-9]+)?)\s*(mg|g|kg)\b", text, re.IGNORECASE
    )
    if m_pmass:
      yield_info["yield_mass"] = f"{m_pmass[-1][0]} {m_pmass[-1][1]}"
    else:
      m_mass = re.findall(
          r"\b([0-9]+(?:\.[0-9]+)?)\s*(mg|g|kg)\b", text, re.IGNORECASE
      )
      if m_mass:
        yield_info["yield_mass"] = f"{m_mass[-1][0]} {m_mass[-1][1]}"

    # Moles
    m_mol = re.findall(
        r"\b([0-9]+(?:\.[0-9]+)?)\s*(mmol|mol|umol)\b", text, re.IGNORECASE
    )
    if m_mol:
      yield_info["yield_moles"] = f"{m_mol[-1][0]} {m_mol[-1][1]}"

    return yield_info

  def extract_solvents(self, text: str) -> List[str]:
    """Identifies reaction and extraction solvents while ignoring NMR blocks."""
    if not text:
      return []

    normalized = _INLINE_CHEM_TAG_RE.sub("", html.unescape(text))
    normalized = normalized.translate(_SUBSCRIPT_DIGITS)
    normalized = _NMR_PAREN_BLOCK_RE.sub(" ", normalized)
    normalized = _DEUTERATED_SOLVENT_TOKEN_RE.sub(" ", normalized)

    found = set()
    if tuple(self.solvents) == SOLVENTS:
      for match in _FORMULA_SOLVENT_RE.finditer(normalized):
        found.add(_FORMULA_SOLVENT_DISPLAY_BY_TOKEN[match.group(0)])
      for pattern, display in _NAME_SOLVENT_PATTERNS:
        if pattern.search(normalized):
          found.add(display)
      return sorted(found)

    text_lower = normalized.lower()
    for s in self.solvents:
      tok = s.translate(_SUBSCRIPT_DIGITS)
      if any(c.isupper() for c in tok):
        if re.search(r"\b" + re.escape(tok) + r"\b", normalized):
          found.add(_SOLVENT_DISPLAY_BY_TOKEN.get(tok, tok))
      else:
        key = tok.lower()
        if re.search(r"\b" + re.escape(key) + r"\b", text_lower):
          found.add(
              _SOLVENT_DISPLAY_BY_TOKEN.get(
                  key, tok.upper() if len(tok) <= 4 else tok.title()
              )
          )
    return sorted(found)

  def extract_reagents(self, text: str) -> List[str]:
    """Identifies coupling reagents, catalysts, and bases."""
    found = set()
    text_lower = text.lower()
    for r in self.reagents:
      pattern = r"\b" + re.escape(r) + r"\b"
      if re.search(pattern, text_lower):
        found.add(r.upper() if len(r) <= 6 else r.title())
    return sorted(list(found))

  def extract_mass_spec(
      self, text: str, mol: Optional[Chem.Mol] = None
  ) -> Dict[str, Any]:
    """Extracts mass spectrometry characterization data and cross-validates calculated m/z."""
    res = {
        "observed_mz": None,
        "adduct": "",
        "calculated_mz": None,
        "delta_mz": None,
        "mass_spec_match": None,
    }
    m_ms = re.search(
        r"m/z\s*(?:=|:)?\s*([0-9]+(?:\.[0-9]+)?)\s*\[?\s*(M\s*[+-]\s*[A-Za-z0-9]+)\s*\]?([+-])?",
        text,
        re.IGNORECASE,
    )
    if not m_ms:
      m_acs = re.search(
          r"m/z\s*(?:=|:)?\s*\[?\s*(M\s*[+-]\s*[A-Za-z0-9]+)\s*\]?([+-])?.*?\b(?:found|obsd|observed)\b\s*(?:=|:|,)?\s*([0-9]{2,4}\.[0-9]+)",
          text,
          re.IGNORECASE,
      ) or re.search(
          r"m/z\s*(?:=|:)?\s*\[?\s*(M\s*[+-]\s*[A-Za-z0-9]+)\s*\]?([+-])?\s*(?:=|:|,)?\s*([0-9]{2,4}\.[0-9]+)",
          text,
          re.IGNORECASE,
      )
      if m_acs:
        adduct_core = m_acs.group(1).replace(" ", "")
        sign = m_acs.group(2) or (
            "+" if "+" in adduct_core else ("-" if "-" in adduct_core else "+")
        )
        adduct = (
            f"[{adduct_core}]{sign}"
            if not adduct_core.endswith(("+", "-"))
            else f"[{adduct_core}]"
        )
        res["observed_mz"] = float(m_acs.group(3))
        res["adduct"] = adduct
      else:
        m_simple = re.search(
            r"m/z\s*(?:=|:)?\s*([0-9]+(?:\.[0-9]+)?)\b", text, re.IGNORECASE
        )
        if m_simple:
          res["observed_mz"] = float(m_simple.group(1))
          res["adduct"] = "[M+H]+"
        else:
          return res
    else:
      obs_mz = float(m_ms.group(1))
      adduct_core = m_ms.group(2).replace(" ", "")
      sign = m_ms.group(3) or (
          "+" if "+" in adduct_core else ("-" if "-" in adduct_core else "+")
      )
      adduct = (
          f"[{adduct_core}]{sign}"
          if not adduct_core.endswith(("+", "-"))
          else f"[{adduct_core}]"
      )
      res["observed_mz"] = obs_mz
      res["adduct"] = adduct

    obs_mz = res["observed_mz"]
    adduct = res["adduct"]

    if mol is not None and obs_mz is not None:
      exact_mw = Descriptors.ExactMolWt(mol)
      shift = ADDUCT_MASS_SHIFTS.get(adduct.lower(), 1.007276)
      calc_mz = round(exact_mw + shift, 4)
      delta_mz = round(abs(obs_mz - calc_mz), 4)
      res["calculated_mz"] = calc_mz
      res["delta_mz"] = delta_mz
      res["mass_spec_match"] = bool(delta_mz <= 0.5)

    return res

  def extract_chiral_metadata(self, text: str) -> Dict[str, Any]:
    """Extracts chiral chromatography fraction, optical rotation ([alpha]D), and enantiomeric excess (% ee)."""
    chiral_peak = None
    m_peak = re.search(r"\b(?:Peak|Fraction)\s*([0-9]+)\b", text, re.IGNORECASE)
    if m_peak:
      chiral_peak = int(m_peak.group(1))
    else:
      if re.search(r"\bfirst[ -]eluting\b", text, re.IGNORECASE):
        chiral_peak = 1
      elif re.search(r"\bsecond[ -]eluting\b", text, re.IGNORECASE):
        chiral_peak = 2
      elif re.search(r"\bthird[ -]eluting\b", text, re.IGNORECASE):
        chiral_peak = 3

    opt_rot = None
    m_rot = re.search(
        r"\[\s*(?:alpha|a)\s*\]\s*D[0-9]*\s*=\s*([+-]?[0-9]+(?:\.[0-9]+)?)",
        text,
        re.IGNORECASE,
    )
    if m_rot:
      opt_rot = float(m_rot.group(1))

    ee_val = None
    m_ee = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*%\s*ee", text, re.IGNORECASE)
    if m_ee:
      ee_val = float(m_ee.group(1))

    return {
        "chiral_peak": chiral_peak,
        "optical_rotation": opt_rot,
        "ee_percent": ee_val,
    }

  def parse_procedure(
      self, text: str, procedure_id: str = "Procedure_01"
  ) -> Dict[str, Any]:
    """Parses full procedure paragraph, extracting components and building target RDKit Mol."""
    cleaned = clean_chemical_text(text)
    prod_name = self.extract_target_product_name(cleaned)
    yields = self.extract_yield(cleaned)
    solvents = self.extract_solvents(cleaned)
    reagents = self.extract_reagents(cleaned)

    free_base_name = prod_name or ""
    detected_salt = None
    salt_stoichiometry = None
    detected_hydrate = None
    hydrate_stoichiometry = None

    if prod_name:
      (
          free_base_name,
          detected_salt,
          salt_stoichiometry,
          detected_hydrate,
          hydrate_stoichiometry,
      ) = self.deformulate_salt_and_hydrate(prod_name, cleaned)

    product_mol = None
    product_smiles = ""
    resolution_engine = "NONE"

    # Build candidate names to resolve
    resolve_candidates = []
    for cand in [free_base_name, prod_name]:
      if not cand:
        continue
      if cand not in resolve_candidates:
        resolve_candidates.append(cand)
      # Clean 'derivative', 'analog', etc.
      c_clean = re.sub(
          r"\s+(?:derivative|analog|salt)\b", "", cand, flags=re.IGNORECASE
      ).strip()
      if c_clean and c_clean not in resolve_candidates:
        resolve_candidates.append(c_clean)
      # Strip leading stereodescriptor prefix e.g. (S)-, (R)-, (2R)-, (±)-
      c_nostereo = re.sub(r"^\s*\([0-9A-Za-z,+-±]+\)-", "", c_clean).strip()
      if c_nostereo and c_nostereo not in resolve_candidates:
        resolve_candidates.append(c_nostereo)

    for cand in resolve_candidates:
      mol = self.resolver.resolve(cand)
      if mol is not None:
        free_mol, free_smi = self.isolate_free_base(mol)
        if free_mol is not None:
          product_mol = free_mol
          product_smiles = free_smi
          resolution_engine = (
              mol.GetProp("ResolutionEngine")
              if mol.HasProp("ResolutionEngine")
              else "RESOLVED"
          )
          break

    # Extract analytical mass spec and chiral metadata
    ms_data = self.extract_mass_spec(cleaned, product_mol)
    chiral_data = self.extract_chiral_metadata(cleaned)

    if product_mol is not None:
      # Attach procedure metadata directly to the RDKit Mol
      product_mol.SetProp("_Name", procedure_id)
      product_mol.SetProp("TargetProductName", prod_name or free_base_name)
      if free_base_name:
        product_mol.SetProp("TargetFreeBaseName", free_base_name)
      if detected_salt:
        product_mol.SetProp("DetectedSalt", detected_salt)
      if salt_stoichiometry is not None:
        product_mol.SetProp("SaltStoichiometry", str(salt_stoichiometry))
      if detected_hydrate:
        product_mol.SetProp("DetectedHydrate", detected_hydrate)
      if hydrate_stoichiometry is not None:
        product_mol.SetProp("HydrateStoichiometry", str(hydrate_stoichiometry))
      if yields.get("yield_percent"):
        product_mol.SetProp("YieldPercent", yields["yield_percent"])
      if yields.get("yield_mass"):
        product_mol.SetProp("YieldMass", yields["yield_mass"])
      if solvents:
        product_mol.SetProp("Solvents", ", ".join(solvents))
      if reagents:
        product_mol.SetProp("Reagents", ", ".join(reagents))

      # Mass spec props
      if ms_data.get("observed_mz") is not None:
        product_mol.SetProp("Observed_m_z", str(ms_data["observed_mz"]))
        product_mol.SetProp("Adduct", ms_data["adduct"])
      if ms_data.get("calculated_mz") is not None:
        product_mol.SetProp("Calculated_m_z", str(ms_data["calculated_mz"]))
        product_mol.SetProp("Delta_m_z", str(ms_data["delta_mz"]))
      if ms_data.get("mass_spec_match") is not None:
        product_mol.SetProp("MassSpecMatch", str(ms_data["mass_spec_match"]))

      # Chiral metadata props
      if chiral_data.get("chiral_peak") is not None:
        product_mol.SetProp("ChiralPeak", str(chiral_data["chiral_peak"]))
      if chiral_data.get("optical_rotation") is not None:
        product_mol.SetProp(
            "OpticalRotation", str(chiral_data["optical_rotation"])
        )
      if chiral_data.get("ee_percent") is not None:
        product_mol.SetProp(
            "EnantiomericExcess", str(chiral_data["ee_percent"])
        )

    return {
        "procedure_id": procedure_id,
        "raw_text": cleaned,
        "target_product_name": prod_name or "",
        "product_name": prod_name or "",
        "target_free_base_name": free_base_name,
        "free_base_name": free_base_name,
        "target_product_smiles": product_smiles,
        "product_smiles": product_smiles,
        "free_base_smiles": product_smiles,
        "Clean_SMILES": product_smiles,
        "canonical_smiles": product_smiles,
        "detected_salt": detected_salt or "",
        "salt_stoichiometry": salt_stoichiometry,
        "stoichiometry_multiplier": (
            float(salt_stoichiometry) if salt_stoichiometry is not None else 1.0
        ),
        "detected_hydrate": detected_hydrate or "",
        "hydrate_stoichiometry": hydrate_stoichiometry,
        "resolution_engine": resolution_engine,
        "yield_percent": yields.get("yield_percent", ""),
        "yield_mass": yields.get("yield_mass", ""),
        "yield_moles": yields.get("yield_moles", ""),
        "solvents": solvents,
        "reagents": reagents,
        "observed_mz": ms_data.get("observed_mz"),
        "adduct": ms_data.get("adduct", ""),
        "calculated_mz": ms_data.get("calculated_mz"),
        "delta_mz": ms_data.get("delta_mz"),
        "mass_spec_match": ms_data.get("mass_spec_match"),
        "chiral_peak": chiral_data.get("chiral_peak"),
        "optical_rotation": chiral_data.get("optical_rotation"),
        "ee_percent": chiral_data.get("ee_percent"),
        "product_mol": product_mol,
    }


def main():
  parser = argparse.ArgumentParser(
      description="Synthetic Experimental Procedure Chemical Extractor."
  )
  parser.add_argument(
      "--input_text",
      "--input-text",
      "--text",
      type=str,
      dest="input_text",
      help="Raw procedure text paragraph.",
  )
  parser.add_argument(
      "--input_file",
      "--input-file",
      "--file",
      type=str,
      dest="input_file",
      help="File containing procedure paragraphs.",
  )
  parser.add_argument(
      "--output_json",
      "--output-json",
      type=str,
      dest="output_json",
      help="Path to write JSON extracted results.",
  )
  parser.add_argument(
      "--output_sdf",
      "--output-sdf",
      type=str,
      dest="output_sdf",
      help="Path to write resolved product molecules.",
  )
  parser.add_argument(
      "--output_csv",
      "--output-csv",
      type=str,
      dest="output_csv",
      help="Path to write CSV extracted results.",
  )
  args = parser.parse_args()

  extractor = ProcedureExtractor()

  if args.input_file:
    with open(args.input_file, "r", encoding="utf-8") as f:
      text_to_parse = f.read()
  elif args.input_text:
    text_to_parse = args.input_text
  else:
    parser.error("Either --input_text or --input_file is required.")

  res = extractor.parse_procedure(text_to_parse)

  print("=" * 60)
  print("Procedure Extraction Summary:")
  print(f"Target Product Name:   {res['target_product_name']}")
  print(f"Target Free Base Name: {res['target_free_base_name']}")
  print(f"Target Product SMILES: {res['target_product_smiles']}")
  print(
      f"Detected Salt:         {res['detected_salt']} (multiplier:"
      f" {res['salt_stoichiometry']})"
  )
  print(
      f"Detected Hydrate:      {res['detected_hydrate']} (multiplier:"
      f" {res['hydrate_stoichiometry']})"
  )
  print(f"Resolution Engine:     {res['resolution_engine']}")
  print(f"Yield:                 {res['yield_percent']}% ({res['yield_mass']})")
  print(f"Solvents:              {', '.join(res['solvents'])}")
  print(f"Reagents:              {', '.join(res['reagents'])}")
  print("=" * 60)

  if args.output_json:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_json)), exist_ok=True
    )
    out_dict = dict(res)
    out_dict.pop("product_mol", None)
    with open(args.output_json, "w") as f:
      json.dump(out_dict, f, indent=2)
    print(f"Saved JSON extraction to: {args.output_json}")

  if args.output_csv:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True
    )
    csv_dict = {
        "target_product_name": [res["target_product_name"]],
        "target_free_base_name": [res["target_free_base_name"]],
        "target_product_smiles": [res["target_product_smiles"]],
        "Clean_SMILES": [res["target_product_smiles"]],
        "canonical_smiles": [res["target_product_smiles"]],
        "detected_salt": [res["detected_salt"]],
        "salt_stoichiometry": [res["salt_stoichiometry"]],
        "stoichiometry_multiplier": [res["stoichiometry_multiplier"]],
        "detected_hydrate": [res["detected_hydrate"]],
        "hydrate_stoichiometry": [res["hydrate_stoichiometry"]],
        "yield_percent": [res["yield_percent"]],
        "yield_mass": [res["yield_mass"]],
        "solvents": ["; ".join(res["solvents"])],
        "reagents": ["; ".join(res["reagents"])],
        "resolution_engine": [res["resolution_engine"]],
        "observed_mz": [res["observed_mz"]],
        "adduct": [res["adduct"]],
        "calculated_mz": [res["calculated_mz"]],
        "delta_mz": [res["delta_mz"]],
        "mass_spec_match": [res["mass_spec_match"]],
        "chiral_peak": [res["chiral_peak"]],
        "optical_rotation": [res["optical_rotation"]],
        "ee_percent": [res["ee_percent"]],
    }
    pd.DataFrame(csv_dict).to_csv(args.output_csv, index=False)
    print(f"Saved CSV extraction to: {args.output_csv}")

  if args.output_sdf and res.get("product_mol") is not None:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_sdf)), exist_ok=True
    )
    with Chem.SDWriter(args.output_sdf) as w:
      w.write(res["product_mol"])
    print(f"Saved product molecule to SDF: {args.output_sdf}")


if __name__ == "__main__":
  main()
