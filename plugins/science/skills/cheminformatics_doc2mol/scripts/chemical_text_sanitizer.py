# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pandas",
#     "rdkit",
# ]
# ///

"""Chemical Text & SMILES Normalization / Sanitization Utility.

Repairs common typography, OCR, and encoding artifacts found in scientific
papers, patent disclosures, and chemical text documents.
"""

import argparse
import html
import os
import re
from typing import Optional

import pandas as pd
from rdkit import Chem

_LOCANT_ZERO_TYPO_RE = re.compile(
    r"(?<![A-Za-z0-9)\]])([0-9]+)(O+)(?=[\],-]|$|\b|[0-9]|[RSEZrs])"
)


def repair_chemical_ocr_text(text: str) -> str:
  """Repairs OCR character defects in chemical names using general IUPAC structural rules.

  Principles:
  1. Prime quotes normalization: normalizes prime symbols (′, ″, ´, `) attached
  to locants.
  2. Contextual Numeric Locant Principle: in numeric locants, stereodescriptors,
     and ring brackets, uppercase 'O' adjacent to digits is digit '0' (1O ->
     10),
     and standalone 'l'/'I' bounded by locant punctuation or stereodescriptors
     is digit '1' ((lR,2S) -> (1R,2S), [l,2-a] -> [1,2-a], lH- -> 1H-).
  3. Contextual Alphabetic Morpheme Principle: inside lowercase alphabetic
  chemical
     stems ([a-z]+), digits '0' and '1' never occur legitimately:
     - Embedded/terminating '0' in lowercase stems is 'o' (fluor0 -> fluoro,
     meth0x -> methox).
     - Stem-ending '-y1' is '-yl' (methy1 -> methyl, indo1-5-y1 -> indol-5-yl).
     - Embedded '1' before n/d in nitrogenous suffixes is 'i' (pyrid1n ->
     pyridin, am1ne -> amine).
     - Remaining embedded/terminating '1' in lowercase stems is 'l' (ch1or ->
     chlor, cyc1o -> cyclo).
  """
  if not text:
    return ""

  s = text

  # 1. Prime quotes normalization
  s = (
      s.replace("\u2032", "'")
      .replace("\u2033", "''")
      .replace("\u00b4", "'")
      .replace("`", "'")
  )
  s = re.sub(r'([0-9A-Za-z])"', r"\1''", s)

  # 2. Contextual Numeric Locant Principle (O -> 0 and l/I -> 1 in locant positions)
  # Compound locants like lO or IO -> 10
  s = re.sub(r"(?:^|(?<=[,(\[-]))[lI]O(?=[\],-]|$|\b)", "10", s)
  # Letter O adjacent to digits in numeric locants (e.g. 1O -> 10,
  # 2,1O- -> 2,10-), excluding formula tokens where the digit run follows
  # letters or closing group delimiters (e.g. Et2O, H2O, Boc2O, C2H6O,
  # (i-Pr)2O).
  while _LOCANT_ZERO_TYPO_RE.search(s):
    s = _LOCANT_ZERO_TYPO_RE.sub(
        lambda m: m.group(1) + "0" * len(m.group(2)), s
    )
  # Standalone l or I bounded by locant punctuation/stereodescriptors (e.g. (lR,2S), [l,2-a], lH-, [2.2.l])
  s = re.sub(r"(?:^|(?<=[,(\[-]))([lI])(?=[,)\]-]|[RSEZrsH]\b)", "1", s)
  s = re.sub(r"(?<=\.)([lI])(?=\]|\.)", "1", s)
  s = re.sub(r"\b([lI])(?=,[0-9]|-[0-9a-zA-Z])", "1", s)

  # 3. Contextual Alphabetic Morpheme Principle (0 -> o and 1 -> l/i in stems)
  # Leading 0 before x in oxa/oxo/oxid stems (e.g. 15-0xo -> 15-oxo)
  s = re.sub(r"(?i)(?:^|(?<=[\s(\[-]))0(?=x[a-z]+)", "o", s)
  # Digit 0 inside or terminating a lowercase stem -> o (e.g. fluor0, pyrazol0)
  while re.search(r"(?<=[a-z])0(?=[a-z]|[-()\[\],\s]|$)", s):
    s = re.sub(r"(?<=[a-z])0(?=[a-z]|[-()\[\],\s]|$)", "o", s)
  s = re.sub(r"(?<=-)0(?=ne\b|l\b)", "o", s)
  # Suffix -y1 at the end of a lowercase stem -> -yl
  s = re.sub(r"(?<=[a-z])y1(?=[^0-9a-zA-Z]|$)", "yl", s)
  # Digit 1 before n/d in amine/amide/azine/idine/aniline stems -> i
  s = re.sub(r"(?<=[a-z]{2})1(?=[nd](?:[a-z]|[-()\[\],\s]|$))", "i", s)
  # Remaining digit 1 embedded inside or terminating lowercase stems -> l (e.g. ch1or, f1uor, cyc1o, indo1)
  while re.search(r"(?<=[a-z])1(?=[a-z]|[-()\[\],\s]|$)", s):
    s = re.sub(r"(?<=[a-z])1(?=[a-z]|[-()\[\],\s]|$)", "l", s)

  # 4. Chemical formula shorthand OCR fixes
  for k, v in [
      ("C02Et", "CO2Et"),
      ("C02Me", "CO2Me"),
      ("C02H", "CO2H"),
      ("C00H", "COOH"),
  ]:
    s = s.replace(k, v)
  s = re.sub(r"\bl-(Pr|Bu)\b", r"i-\1", s)
  s = re.sub(r"\bl(Pr|Bu)\b", r"i\1", s)

  return s


def clean_chemical_text(text: str) -> str:
  """Decodes HTML/XML entities, strips markup, normalizes Greek characters, and repairs OCR stems."""
  if not text:
    return ""

  # Decode HTML entities (&alpha;, &minus;, &prime;, etc.)
  cleaned = html.unescape(text.strip())

  # Strip inline HTML/XML formatting tags (<i>, <b>, <sub>, <sup>, <span>)
  cleaned = re.sub(
      r"</?(?:i|b|sub|sup|em|strong|span)\b[^>]*>",
      "",
      cleaned,
      flags=re.IGNORECASE,
  )

  # Remove soft hyphens and zero-width spaces
  cleaned = (
      cleaned.replace("\xad", "").replace("\u200b", "").replace("\ufeff", "")
  )

  # Replace typographic dashes and linebreaks inside words
  cleaned = re.sub(r"(\w+)-\s*\n\s*(\w+)", r"\1-\2", cleaned)
  cleaned = cleaned.replace("—", "-").replace("–", "-").replace("−", "-")

  # Collapse whitespace around hyphens inside chemical locants (e.g. "2 - methyl" -> "2-methyl", "N - (2-" -> "N-(2-")
  cleaned = re.sub(
      r"\b([0-9]+|[NOPS])\s+-\s+([a-zA-Z0-9(\[])", r"\1-\2", cleaned
  )

  # Replace typographic quotes
  cleaned = (
      cleaned.replace("“", '"')
      .replace("”", '"')
      .replace("‘", "'")
      .replace("’", "'")
  )

  # Transliterate common Greek symbols to English words for IUPAC parsing
  greek_map = {
      "α": "alpha",
      "β": "beta",
      "γ": "gamma",
      "δ": "delta",
      "ε": "epsilon",
      "ζ": "zeta",
      "η": "eta",
      "θ": "theta",
      "κ": "kappa",
      "λ": "lambda",
      "μ": "mu",
      "ν": "nu",
      "ξ": "xi",
      "π": "pi",
      "ρ": "rho",
      "σ": "sigma",
      "τ": "tau",
      "υ": "upsilon",
      "φ": "phi",
      "χ": "chi",
      "ψ": "psi",
      "ω": "omega",
      "±": "+/-",
  }
  for g_char, g_word in greek_map.items():
    cleaned = cleaned.replace(g_char, g_word)

  # Contextual chemical OCR stem repair
  cleaned = repair_chemical_ocr_text(cleaned)

  return cleaned


def normalize_smiles_string(smiles: str) -> str:
  """Repairs common OCR, shorthand substituents, and notation defects in SMILES strings.

  Key principles:
  1. Substituent Shorthand Expansion: expands parenthesized shorthands ((CF3),
  (Boc), (Ac), etc.)
     and unbracketed perfluoroalkyls (CF3, OCF3, SCF3).
  2. Hypervalent Uncharged Functional Groups: normalizes N(=O)=O -> [N+](=O)[O-]
     and N=N=N -> N=[N+]=[N-].
  3. Ring-Closure Parity Principle for Digit '0': in SMILES, ring closure digits
  must
     appear an even number of times. '=0' is always carbonyl/sulfonyl oxygen
     '=O', and
     if digit '0' has an odd count, it cannot be a ring closure and is repaired
     to 'O'.
  """
  if not smiles:
    return ""

  smi = smiles.strip()

  # 1. Expand common substituent shorthands when enclosed in parentheses
  shorthand_map = {
      "(CF3)": "(C(F)(F)F)",
      "(OCF3)": "(OC(F)(F)F)",
      "(SCF3)": "(SC(F)(F)F)",
      "(C2F5)": "(C(F)(F)C(F)(F)F)",
      "(NO2)": "([N+](=O)[O-])",
      "(N3)": "(N=[N+]=[N-])",
      "(CN)": "(C#N)",
      "(CHO)": "(C=O)",
      "(COOH)": "(C(=O)O)",
      "(CO2H)": "(C(=O)O)",
      "(Me)": "(C)",
      "(Et)": "(CC)",
      "(nPr)": "(CCC)",
      "(iPr)": "(C(C)C)",
      "(nBu)": "(CCCC)",
      "(iBu)": "(CC(C)C)",
      "(tBu)": "(C(C)(C)C)",
      "(OMe)": "(OC)",
      "(OEt)": "(OCC)",
      "(OiPr)": "(OC(C)C)",
      "(OtBu)": "(OC(C)(C)C)",
      "(Ac)": "(C(=O)C)",
      "(OAc)": "(OC(=O)C)",
      "(NHAc)": "(NC(=O)C)",
      "(Boc)": "(C(=O)OC(C)(C)C)",
      "(Ph)": "(c1ccccc1)",
      "(Bn)": "(Cc1ccccc1)",
  }
  for k, v in shorthand_map.items():
    smi = smi.replace(k, v)

  # Expand unbracketed CF3, OCF3, SCF3 anywhere outside brackets
  smi = re.sub(r"(?<!\[)OCF3(?![a-zA-Z0-9\]])", "OC(F)(F)F", smi)
  smi = re.sub(r"(?<!\[)SCF3(?![a-zA-Z0-9\]])", "SC(F)(F)F", smi)
  smi = re.sub(r"(?<!\[)CF3(?![a-zA-Z0-9\]])", "C(F)(F)F", smi)

  # 2. Repair pentavalent uncharged nitro and uncharged azide
  smi = re.sub(r"N\(=O\)=O", "[N+](=O)[O-]", smi)
  smi = re.sub(r"N=N=N", "N=[N+]=[N-]", smi)

  # 3. Combinatorial Ring-Closure vs. Oxygen Disambiguation for Digit '0'
  # Unambiguous syntactic positions: '=0', '(0', '^0', '[0' are always oxygen 'O'
  smi = smi.replace("=0", "=O").replace("[0-]", "[O-]").replace("[0H]", "[OH]")
  smi = re.sub(r"(^|[(.\-])0", r"\1O", smi)

  # If digit '0' is present, disambiguate which '0's are ring closures vs. oxygen 'O'
  # by scoring valid RDKit parses based on ring strain (5/6/7-membered rings score +10; 3/4-membered score -10)
  zero_indices = [m.start() for m in re.finditer(r"0", smi)]
  if zero_indices:
    best_candidate = smi
    best_score = -999
    candidates_to_test = [smi]
    for idx in zero_indices:
      candidates_to_test.append(smi[:idx] + "O" + smi[idx + 1 :])
    candidates_to_test.append(smi.replace("0", "O"))

    for cand in candidates_to_test:
      mol_cand = Chem.MolFromSmiles(cand)
      if mol_cand is not None:
        ring_sizes = [len(r) for r in mol_cand.GetRingInfo().AtomRings()]
        score = sum(
            10 if 5 <= r <= 7 else (-10 if r < 5 else 2) for r in ring_sizes
        )
        if score > best_score:
          best_score = score
          best_candidate = cand
    smi = best_candidate

  # 4. Strip whitespace inside SMILES
  smi = re.sub(r"\s+", "", smi)

  return smi


def repair_and_sanitize_smiles(raw_smiles: str) -> Optional[Chem.Mol]:
  """Parses and sanitizes a potentially defective SMILES string using graph-based valence healing."""
  if not raw_smiles:
    return None

  # Pass 1: Normalize shorthands and disambiguate ring digits vs OCR oxygens
  normalized = normalize_smiles_string(raw_smiles)
  mol = Chem.MolFromSmiles(normalized)
  if mol is not None:
    return mol

  # Pass 3: Non-sanitized parse to diagnose and heal valencies on the molecular graph
  try:
    mol = Chem.MolFromSmiles(normalized, sanitize=False)
    if mol is None:
      return None

    mol.UpdatePropertyCache(strict=False)

    # Heal uncharged pentavalent nitro groups and hypervalent N-oxides
    for atom in mol.GetAtoms():
      if atom.GetSymbol() == "N" and atom.GetFormalCharge() == 0:
        valence = int(sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()))
        if valence == 5:
          oxygens = [
              nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "O"
          ]
          if len(oxygens) >= 2:
            atom.SetFormalCharge(1)
            for o in oxygens:
              bond = mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx())
              if bond.GetBondType() == Chem.BondType.DOUBLE:
                bond.SetBondType(Chem.BondType.SINGLE)
                o.SetFormalCharge(-1)
                break
          elif len(oxygens) == 1:
            # Tertiary or aromatic N-oxide written with double bond (R3N=O -> R3[N+][O-])
            o = oxygens[0]
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx())
            if bond.GetBondType() == Chem.BondType.DOUBLE:
              atom.SetFormalCharge(1)
              bond.SetBondType(Chem.BondType.SINGLE)
              o.SetFormalCharge(-1)

    # Now attempt standard RDKit sanitization
    Chem.SanitizeMol(mol)
    return mol
  except Exception:
    return None


def sanitize_pdf_smiles(raw_smiles: str) -> str:
  """Repairs PDF/OCR layout-corrupted SMILES strings.

  Handles:
  1. Soft hyphens and line-break splits ('-\\n', '-\\r\\n', '\\xad').
  2. Internal whitespace within SMILES lines.
  3. Unicode minuses, en-dashes, em-dashes ('\\u2212', '\\u2013', '\\u2014').
  4. Pentavalent uncharged nitro syntax (N(=O)=O -> [N+](=O)[O-]).
  5. 5-carbon aromatic carbocycle dropouts (e.g. c1cc(...)cc1 -> c1ccc(...)cc1).
  """
  if not raw_smiles:
    return ""

  s = str(raw_smiles).strip()
  # Remove soft hyphens and hidden zero-width characters
  s = s.replace("\xad", "").replace("\u200b", "").replace("\ufeff", "")
  # Remove hyphen linebreaks
  s = re.sub(r"-\s*[\r\n]+\s*", "", s)
  # Remove all internal whitespace and regular newlines
  s = re.sub(r"\s+", "", s)
  # Normalize typographic dashes and minus signs
  s = s.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-")
  # Expand shorthand unbracketed CF3 and uncharged nitro
  s = re.sub(r"(?<!\[)CF3(?![a-zA-Z0-9\]])", "C(F)(F)F", s)
  s = re.sub(r"(?<!\[)N\(=O\)=O", "[N+](=O)[O-]", s)

  # Check if directly valid
  m = Chem.MolFromSmiles(s)
  if m is not None:
    return Chem.MolToSmiles(m)

  # Check for 5-membered aromatic carbocycle dropout: c1cc(...)cc1
  if "c1cc(" in s:
    candidate = re.sub(r"c1cc\(", "c1ccc(", s)
    m = Chem.MolFromSmiles(candidate)
    if m is not None:
      return Chem.MolToSmiles(m)

  # Attempt fallback through repair_and_sanitize_smiles
  m_fallback = repair_and_sanitize_smiles(s)
  if m_fallback is not None:
    return Chem.MolToSmiles(m_fallback)

  return s


def main():
  parser = argparse.ArgumentParser(
      description="Chemical Text & SMILES Normalization / Sanitization Utility."
  )
  parser.add_argument(
      "--repair-smiles-file",
      "--repair_smiles_file",
      type=str,
      help="Path to text file with corrupted SMILES separated by '---'.",
  )
  parser.add_argument(
      "--output-csv",
      "--output_csv",
      type=str,
      help="Output CSV path for repaired SMILES.",
  )
  parser.add_argument(
      "--smiles", type=str, help="Single SMILES string to repair."
  )
  args = parser.parse_args()

  if args.repair_smiles_file:
    with open(args.repair_smiles_file, "r", encoding="utf-8") as f:
      content = f.read()

    # Split by '---' separators
    chunks = [
        c.strip()
        for c in re.split(r"(?:^|\n)\s*---\s*(?:\n|$)", content)
        if c.strip()
    ]

    repaired_list = []
    for c in chunks:
      rep = sanitize_pdf_smiles(c)
      m = Chem.MolFromSmiles(rep)
      repaired_list.append({
          "original_smiles": c,
          "repaired_smiles": rep,
          "is_valid": m is not None,
      })

    df = pd.DataFrame(repaired_list)
    if args.output_csv:
      os.makedirs(
          os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True
      )
      df.to_csv(args.output_csv, index=False)
      print(f"Repaired {len(df)} SMILES records saved to {args.output_csv}")
    else:
      for r in repaired_list:
        print(f"Original: {r['original_smiles']}")
        print(f"Repaired: {r['repaired_smiles']} (Valid: {r['is_valid']})\n")

  elif args.smiles:
    rep = sanitize_pdf_smiles(args.smiles)
    m = Chem.MolFromSmiles(rep)
    print(f"Original: {args.smiles}")
    print(f"Repaired: {rep} (Valid: {m is not None})")

  else:
    parser.error("Provide --smiles, --text, --input_file, or --input_csv.")


if __name__ == "__main__":
  main()
