# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
# ]
# ///

"""Lightweight Reaction Atom Conservation & Stoichiometry Validator.

Verifies fundamental physical and chemical conservation laws between synthetic
precursors, byproducts, and products (especially for macrocyclizations and
prospective analog design):
1. Strict Carbon Conservation: ΔC = 0 (prohibits magical creation or loss of
   carbon atoms).
2. Heavy Atom Stoichiometry: Balances every element (N, O, F, Cl, Br, S, P, B,
   ...) against expected reaction byproducts, with multiplicity.
3. Cyclization Ring Count: Verifies ring count increases by +1 and detects
   macrocyclic ring size (>= 9).
4. Formula & Mass Balance: Reports exact elemental delta.

Multi-component reactions are supported: pass several precursors/products.
Byproducts accept an integer coefficient (``2*H2O``) and may be given as a
name, a formula (parentheses and nesting supported, e.g. ``Al(OH)3``), or a
SMILES. A token that is valid as both a formula and a SMILES (e.g. ``CCO``) is
reported as ambiguous; prefix it with ``formula:`` or ``smiles:`` to be
explicit.
"""

import argparse
from collections import Counter
import json
import os
import re
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from rdkit import Chem
from rdkit.Chem import Lipinski, rdMolDescriptors

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
  sys.path.insert(0, _SCRIPT_DIR)
import prochirality

COMMON_BYPRODUCT_FORMULAS: Dict[str, Dict[str, int]] = {
    "h2o": {"H": 2, "O": 1},
    "water": {"H": 2, "O": 1},
    "hf": {"H": 1, "F": 1},
    "hcl": {"H": 1, "Cl": 1},
    "hbr": {"H": 1, "Br": 1},
    "hi": {"H": 1, "I": 1},
    "nh3": {"N": 1, "H": 3},
    "b(oh)3": {"B": 1, "H": 3, "O": 3},
    "boric_acid": {"B": 1, "H": 3, "O": 3},
    "pinacol": {"C": 6, "H": 14, "O": 2},
    "pinacol_borate": {"C": 6, "H": 13, "B": 1, "O": 3},
    "acetic_acid": {"C": 2, "H": 4, "O": 2},
    "co2": {"C": 1, "O": 2},
    "n2": {"N": 2},
}

PRESET_REACTION_BYPRODUCTS: Dict[str, List[Dict[str, int]]] = {
    "amide": [COMMON_BYPRODUCT_FORMULAS["h2o"]],
    "lactam": [COMMON_BYPRODUCT_FORMULAS["h2o"]],
    "lactamization": [COMMON_BYPRODUCT_FORMULAS["h2o"]],
    "ester": [COMMON_BYPRODUCT_FORMULAS["h2o"]],
    "lactone": [COMMON_BYPRODUCT_FORMULAS["h2o"]],
    "mitsunobu": [COMMON_BYPRODUCT_FORMULAS["h2o"]],
    "snar_fluoro": [COMMON_BYPRODUCT_FORMULAS["hf"]],
    "snar_chloro": [COMMON_BYPRODUCT_FORMULAS["hcl"]],
    "snar": [COMMON_BYPRODUCT_FORMULAS["hf"]],
    "suzuki": [
        COMMON_BYPRODUCT_FORMULAS["hbr"],
        COMMON_BYPRODUCT_FORMULAS["b(oh)3"],
    ],
    "suzuki_pinacol": [
        COMMON_BYPRODUCT_FORMULAS["hbr"],
        COMMON_BYPRODUCT_FORMULAS["pinacol_borate"],
    ],
    "buchwald": [COMMON_BYPRODUCT_FORMULAS["hbr"]],
    "click": [],
    "azide_alkyne": [],
    "ring_closing_metathesis": [{"C": 2, "H": 4}],  # ethylene
    "rcm": [{"C": 2, "H": 4}],
    "macrocyclization": [],  # auto-detected
}


_PERIODIC_TABLE = Chem.GetPeriodicTable()

# Tokens: an element symbol with optional count, "(" , ")" with optional count.
_FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\d*)|(\()|(\))(\d*)")


def _is_element(symbol: str) -> bool:
  """Returns True if `symbol` is a real element symbol."""
  try:
    return _PERIODIC_TABLE.GetAtomicNumber(symbol) > 0
  except Exception:
    return False


def parse_formula(formula_str: str) -> Optional[Counter]:
  """Parses a molecular formula into element counts.

  Supports a leading stoichiometric coefficient and nested parenthesised
  groups, e.g. ``"2 H2O"`` -> ``{H: 4, O: 2}`` and ``"Al(OH)3"`` ->
  ``{Al: 1, O: 3, H: 3}``.

  Every symbol is validated against the periodic table, so a string that is not
  a well-formed formula (including most SMILES) yields ``None`` rather than a
  plausible-looking but wrong count.

  Args:
    formula_str: Formula text, optionally prefixed by an integer coefficient.

  Returns:
    A Counter of element -> count, or None if the string is not a valid formula.
  """
  text = formula_str.strip()
  if not text:
    return None

  # Leading stoichiometric coefficient, e.g. "2 H2O", "2H2O", "2*H2O".
  coefficient = 1
  lead = re.match(r"^(\d+)\s*[*x]?\s*", text)
  if lead:
    coefficient = int(lead.group(1))
    if coefficient == 0:
      return None
    text = text[lead.end() :].strip()
    if not text:
      return None

  counts: Counter = Counter()
  stack: List[Counter] = [counts]
  pos = 0
  while pos < len(text):
    if text[pos].isspace():
      pos += 1
      continue
    match = _FORMULA_TOKEN_RE.match(text, pos)
    if not match or match.end() == pos:
      return (
          None  # Unparseable character (bond symbol, lowercase aromatic, ...).
      )
    symbol, num, open_paren, close_paren, group_num = match.groups()
    if symbol:
      if not _is_element(symbol):
        return None
      stack[-1][symbol] += int(num) if num else 1
    elif open_paren:
      stack.append(Counter())
    elif close_paren:
      if len(stack) == 1:
        return None  # Unbalanced ")".
      group = stack.pop()
      multiplier = int(group_num) if group_num else 1
      for elem, n in group.items():
        stack[-1][elem] += n * multiplier
    pos = match.end()

  if len(stack) != 1:
    return None  # Unbalanced "(".
  if not counts:
    return None

  if coefficient != 1:
    for elem in list(counts):
      counts[elem] *= coefficient
  return counts


def resolve_species(token: str) -> Tuple[Optional[Counter], str, List[str]]:
  """Resolves a byproduct token to element counts.

  A bare token such as ``"CCO"`` is both a valid formula (C2O) and a valid
  SMILES (ethanol, C2H6O). Guessing silently is how hydrogens go missing, so
  ambiguous tokens are reported. Prefix a token with ``smiles:`` or
  ``formula:`` to state the intent explicitly.

  Args:
    token: Byproduct name, formula, or SMILES, with optional coefficient and
      optional ``smiles:``/``formula:`` prefix.

  Returns:
    Tuple of (counts or None, human-readable description, warnings).
  """
  warnings: List[str] = []
  text = token.strip()

  # Split off a leading coefficient so it can be applied to either branch.
  coefficient = 1
  lead = re.match(r"^(\d+)\s*[*x]\s*", text)
  if lead:
    coefficient = int(lead.group(1))
    text = text[lead.end() :].strip()

  def scaled(counts: Counter) -> Counter:
    if coefficient == 1:
      return counts
    return Counter({e: n * coefficient for e, n in counts.items()})

  prefix_label = f"{coefficient}x " if coefficient != 1 else ""

  if text.lower().startswith("smiles:"):
    smi = text[len("smiles:") :].strip()
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
      return None, text, [f"Could not parse SMILES byproduct '{smi}'."]
    return (
        scaled(get_mol_atom_counts(mol, include_hydrogens=True)),
        f"{prefix_label}{smi} (SMILES)",
        warnings,
    )

  if text.lower().startswith("formula:"):
    formula = text[len("formula:") :].strip()
    parsed = parse_formula(formula)
    if parsed is None:
      return None, text, [f"Could not parse formula byproduct '{formula}'."]
    return scaled(parsed), f"{prefix_label}{formula} (formula)", warnings

  if text.lower() in COMMON_BYPRODUCT_FORMULAS:
    known = Counter(COMMON_BYPRODUCT_FORMULAS[text.lower()])
    return scaled(known), f"{prefix_label}{text}", warnings

  as_formula = parse_formula(text)
  mol = Chem.MolFromSmiles(text)
  as_smiles = (
      get_mol_atom_counts(mol, include_hydrogens=True)
      if mol is not None
      else None
  )

  if as_formula is None and as_smiles is None:
    return (
        None,
        text,
        [f"Could not parse byproduct '{text}' as a formula or SMILES."],
    )
  if as_formula is None:
    return scaled(as_smiles), f"{prefix_label}{text} (SMILES)", warnings
  if as_smiles is None:
    return scaled(as_formula), f"{prefix_label}{text} (formula)", warnings

  if as_formula != as_smiles:
    warnings.append(
        f"Byproduct '{text}' is ambiguous: as a formula it is"
        f" {format_counts(as_formula)}, as SMILES it is"
        f" {format_counts(as_smiles)}. Interpreting it as a formula. Use"
        f" 'smiles:{text}' or 'formula:{text}' to state which you mean."
    )
  return scaled(as_formula), f"{prefix_label}{text} (formula)", warnings


def format_counts(counts: Dict[str, int]) -> str:
  """Renders element counts in Hill order for messages."""
  ordered: List[str] = []
  for elem in ("C", "H"):
    if counts.get(elem):
      ordered.append(f"{elem}{counts[elem] if counts[elem] > 1 else ''}")
  for elem in sorted(e for e in counts if e not in ("C", "H")):
    if counts[elem]:
      ordered.append(f"{elem}{counts[elem] if counts[elem] > 1 else ''}")
  return "".join(ordered) if ordered else "(empty)"


def get_mol_atom_counts(
    mol: Chem.Mol, include_hydrogens: bool = True
) -> Counter:
  """Returns elemental counts of a molecule including explicit and implicit hydrogens."""
  counts: Counter = Counter()
  if mol is None:
    return counts

  # Add explicit hydrogens to guarantee exact H perception
  m_with_h = Chem.AddHs(mol) if include_hydrogens else mol
  for atom in m_with_h.GetAtoms():
    sym = atom.GetSymbol()
    counts[sym] += 1
  return counts


def get_macrocyclic_ring_sizes(
    mol: Chem.Mol, min_ring_size: int = 12
) -> List[int]:
  """Returns list of ring sizes for macrocyclic rings (size >= 12)."""
  if mol is None:
    return []
  ring_info = mol.GetRingInfo()
  atom_rings = ring_info.AtomRings()
  macro_sizes = [len(r) for r in atom_rings if len(r) >= min_ring_size]
  return sorted(macro_sizes)


def verify_reaction_stoichiometry(
    precursors: Union[str, Chem.Mol, Sequence[Union[str, Chem.Mol]]],
    products: Union[str, Chem.Mol, Sequence[Union[str, Chem.Mol]]],
    reaction_type: Optional[str] = None,
    byproducts: Optional[Sequence[Union[str, Dict[str, int]]]] = None,
    is_cyclization: Optional[bool] = None,
    expected_macrocycle_size: Optional[int] = None,
    strict_hydrogen: bool = False,
    chiral_catalyst: bool = False,
) -> Dict[str, Any]:
  """Verifies conservation of atoms and mass between precursors and products.

  Args:
    precursors: Single SMILES/Mol or list of precursors.
    products: Single SMILES/Mol or list of products.
    reaction_type: Name of preset reaction (e.g. 'amide', 'mitsunobu', 'snar',
      'suzuki').
    byproducts: Explicit list of byproduct formulas (e.g. ['H2O', 'HBr',
      'B(OH)3']) or dicts.
    is_cyclization: True if validating an intramolecular cyclization /
      macrocyclization.
    expected_macrocycle_size: Expected ring size of the newly formed macrocycle
      (e.g. 12, 14, 16).

  Returns:
    Diagnostic dictionary with verification status, delta atom counts, and
    errors.
  """
  errors: List[str] = []
  warnings: List[str] = []

  # 1. Parse Precursors
  prec_list = (
      [precursors]
      if isinstance(precursors, (str, Chem.Mol))
      else list(precursors)
  )
  prec_mols: List[Chem.Mol] = []
  prec_atoms: Counter = Counter()

  for idx, p in enumerate(prec_list):
    if isinstance(p, str):
      m = Chem.MolFromSmiles(p)
      if m is None:
        errors.append(f"Invalid precursor #{idx+1} SMILES: '{p}'")
        continue
    else:
      m = p
    prec_mols.append(m)
    prec_atoms.update(get_mol_atom_counts(m, include_hydrogens=True))

  # 2. Parse Products
  prod_list = (
      [products] if isinstance(products, (str, Chem.Mol)) else list(products)
  )
  prod_mols: List[Chem.Mol] = []
  prod_atoms: Counter = Counter()

  for idx, p in enumerate(prod_list):
    if isinstance(p, str):
      m = Chem.MolFromSmiles(p)
      if m is None:
        errors.append(f"Invalid product #{idx+1} SMILES: '{p}'")
        continue
    else:
      m = p
    prod_mols.append(m)
    prod_atoms.update(get_mol_atom_counts(m, include_hydrogens=True))

  if errors:
    return {
        "valid": False,
        "status": "PARSING_ERROR",
        "errors": errors,
        "warnings": warnings,
    }

  # 3. Resolve Byproducts
  byprod_atoms: Counter = Counter()
  applied_byproducts: List[str] = []

  if byproducts:
    for bp in byproducts:
      if isinstance(bp, str):
        counts, label, bp_warnings = resolve_species(bp)
        warnings.extend(bp_warnings)
        if counts is not None:
          byprod_atoms.update(counts)
          applied_byproducts.append(label)
      elif isinstance(bp, dict):
        byprod_atoms.update(bp)
        applied_byproducts.append(format_counts(bp))
  elif reaction_type:
    rxn_key = reaction_type.lower().strip()
    if rxn_key in PRESET_REACTION_BYPRODUCTS:
      for bp_dict in PRESET_REACTION_BYPRODUCTS[rxn_key]:
        byprod_atoms.update(bp_dict)
      applied_byproducts.append(rxn_key)
    else:
      warnings.append(
          f"Unknown preset reaction_type '{reaction_type}'. Specify byproducts"
          " explicitly."
      )
  elif is_cyclization:
    # Auto-infer standard condensation byproduct for explicit intramolecular cyclizations
    delta_heavy = sum(
        prec_atoms[elem] for elem in prec_atoms if elem != "H"
    ) - sum(prod_atoms[elem] for elem in prod_atoms if elem != "H")
    if delta_heavy == 1 and (prec_atoms["O"] - prod_atoms["O"] == 1):
      byprod_atoms.update(COMMON_BYPRODUCT_FORMULAS["h2o"])
      applied_byproducts.append("H2O (inferred condensation)")
    elif delta_heavy == 1 and (prec_atoms["F"] - prod_atoms["F"] == 1):
      byprod_atoms.update(COMMON_BYPRODUCT_FORMULAS["hf"])
      applied_byproducts.append("HF (inferred SNAr)")
    elif delta_heavy == 1 and (prec_atoms["Cl"] - prod_atoms["Cl"] == 1):
      byprod_atoms.update(COMMON_BYPRODUCT_FORMULAS["hcl"])
      applied_byproducts.append("HCl (inferred condensation)")
    elif delta_heavy == 1 and (prec_atoms["Br"] - prod_atoms["Br"] == 1):
      byprod_atoms.update(COMMON_BYPRODUCT_FORMULAS["hbr"])
      applied_byproducts.append("HBr (inferred condensation)")
  elif len(prec_mols) == 1 and len(prod_mols) == 1:
    delta_heavy = sum(
        prec_atoms[elem] for elem in prec_atoms if elem != "H"
    ) - sum(prod_atoms[elem] for elem in prod_atoms if elem != "H")
    for elem, bp_name in (
        ("O", "H2O"),
        ("F", "HF"),
        ("Cl", "HCl"),
        ("Br", "HBr"),
    ):
      if delta_heavy == 1 and (prec_atoms[elem] - prod_atoms[elem] == 1):
        warnings.append(
            f"Hint: Product is missing {bp_name} relative to precursor; if this"
            " is a condensation or cyclization reaction, pass --byproducts"
            f" {bp_name} or --cyclization."
        )

  # 4. Compute Elemental Balances
  all_elements = (
      set(prec_atoms.keys()) | set(prod_atoms.keys()) | set(byprod_atoms.keys())
  )
  delta_atoms: Dict[str, int] = {}
  carbon_conserved = True
  heavy_atoms_conserved = True

  for elem in sorted(all_elements):
    n_prec = prec_atoms[elem]
    n_prod_total = prod_atoms[elem] + byprod_atoms[elem]
    diff = (
        n_prod_total - n_prec
    )  # positive: gained atoms in product, negative: lost atoms
    delta_atoms[elem] = diff

    if elem == "C" and diff != 0:
      carbon_conserved = False
      verb = "gained" if diff > 0 else "lost"
      errors.append(
          f"Carbon count mismatch! Precursor has {n_prec} C, Product has"
          f" {prod_atoms['C']} C (net {verb} {abs(diff)} carbon"
          f" atom{'s' if abs(diff) > 1 else ''}). Stoichiometric violation of"
          " carbon conservation."
      )

    if elem != "H" and diff != 0:
      heavy_atoms_conserved = False
      verb = "gained" if diff > 0 else "lost"
      errors.append(
          f"Heavy atom mismatch on '{elem}': Precursors have {n_prec}, Products"
          f" + Byproducts have {n_prod_total} (net {verb} {abs(diff)}"
          f" '{elem}')."
      )

  # Check Hydrogen balance (warn if small difference due to formal protonation state)
  delta_h = delta_atoms.get("H", 0)
  notes: List[str] = []
  if delta_h != 0:
    if strict_hydrogen:
      verb = "gained" if delta_h > 0 else "lost"
      errors.append(
          "Hydrogen count mismatch under --strict-hydrogen! Precursor has"
          f" {prec_atoms.get('H', 0)} H, Product + Byproducts has"
          f" {prod_atoms.get('H', 0) + byprod_atoms.get('H', 0)} H (net {verb}"
          f" {abs(delta_h)} H). Stoichiometric violation of hydrogen"
          " conservation."
      )
    else:
      notes.append(
          f"hydrogen imbalance: {delta_h:+d} (expected for ionic steps)"
      )
      warnings.append(
          f"Hydrogen balance offset (ΔH = {delta_h:+d}). Ensure formal"
          " ionization states and proton transfers match the balanced"
          " equation."
      )

  # 5. Cyclization & Macrocyclic Ring Verification
  ring_diagnostics: Dict[str, Any] = {}
  if is_cyclization or (
      len(prec_mols) == 1
      and len(prod_mols) == 1
      and expected_macrocycle_size is not None
  ):
    prec_mol = prec_mols[0]
    prod_mol = prod_mols[0]

    prec_rings = Lipinski.RingCount(prec_mol)
    prod_rings = Lipinski.RingCount(prod_mol)
    delta_rings = prod_rings - prec_rings
    macro_sizes = get_macrocyclic_ring_sizes(prod_mol)

    ring_diagnostics = {
        "precursor_ring_count": prec_rings,
        "product_ring_count": prod_rings,
        "delta_rings": delta_rings,
        "product_macrocycle_ring_sizes": macro_sizes,
    }

    if delta_rings < 1:
      errors.append(
          "Intramolecular cyclization failed to form a new ring! Precursor"
          f" ring count: {prec_rings}, Product ring count: {prod_rings} (Δ ="
          f" {delta_rings})."
      )

    if expected_macrocycle_size is not None:
      if expected_macrocycle_size not in macro_sizes:
        errors.append(
            f"Expected {expected_macrocycle_size}-membered macrocyclic ring, "
            f"but detected macrocycle ring sizes in product are: {macro_sizes}."
        )

  # 6. Prochirality & Stereogenesis Audit
  prec_smis = [
      p if isinstance(p, str) else Chem.MolToSmiles(p, isomericSmiles=True)
      for p in prec_list
  ]
  prod_smis = [
      p if isinstance(p, str) else Chem.MolToSmiles(p, isomericSmiles=True)
      for p in prod_list
  ]
  stereogenesis_diagnostics: Dict[str, Any] = (
      prochirality.verify_reaction_stereogenesis(
          prec_smis, prod_smis, chiral_catalyst=chiral_catalyst
      )
  )
  errors.extend(stereogenesis_diagnostics.get("errors", []))
  warnings.extend(stereogenesis_diagnostics.get("warnings", []))
  notes.extend(stereogenesis_diagnostics.get("notes", []))

  hydrogen_conserved = delta_h == 0
  is_valid = (
      len(errors) == 0
      and carbon_conserved
      and heavy_atoms_conserved
      and (not strict_hydrogen or hydrogen_conserved)
  )

  prec_formula = (
      rdMolDescriptors.CalcMolFormula(prec_mols[0])
      if len(prec_mols) == 1
      else " + ".join([rdMolDescriptors.CalcMolFormula(m) for m in prec_mols])
  )
  prod_formula = (
      rdMolDescriptors.CalcMolFormula(prod_mols[0])
      if len(prod_mols) == 1
      else " + ".join([rdMolDescriptors.CalcMolFormula(m) for m in prod_mols])
  )

  return {
      "valid": is_valid,
      "status": "PASSED" if is_valid else "FAILED_STOICHIOMETRY",
      "carbon_conserved": carbon_conserved,
      "heavy_atoms_conserved": heavy_atoms_conserved,
      "hydrogen_conserved": hydrogen_conserved,
      "precursor_formula": prec_formula,
      "product_formula": prod_formula,
      "applied_byproducts": applied_byproducts,
      "delta_atoms": {k: v for k, v in delta_atoms.items() if v != 0},
      "ring_diagnostics": ring_diagnostics,
      "stereogenesis_diagnostics": stereogenesis_diagnostics,
      "errors": errors,
      "warnings": warnings,
      "notes": notes,
  }


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Lightweight Reaction Atom Conservation & Stoichiometry Validator."
      )
  )
  parser.add_argument(
      "--precursor",
      "--prec",
      required=True,
      nargs="+",
      help=(
          "Precursor SMILES. Repeat for multi-component reactions (e.g."
          " --precursor A B)."
      ),
  )
  parser.add_argument(
      "--product",
      "--prod",
      required=True,
      nargs="+",
      help=(
          "Product SMILES. Repeat if the reaction yields more than one product."
      ),
  )
  parser.add_argument(
      "--reaction-type",
      "--rxn",
      default=None,
      help="Preset reaction type (e.g. amide, mitsunobu, snar, suzuki).",
  )
  parser.add_argument(
      "--byproducts",
      "--byprod",
      nargs="*",
      default=None,
      help=(
          "Explicit byproducts as formulas, names, or SMILES (e.g. H2O HCl"
          " 'B(OH)3'). Prefix an integer coefficient with '*' or 'x' to repeat"
          " (e.g. '2*H2O'). Disambiguate a token that is both a formula and a"
          " SMILES with 'formula:' or 'smiles:' (e.g. 'smiles:CCO' for ethanol,"
          " not C2O)."
      ),
  )
  parser.add_argument(
      "--cyclization",
      action="store_true",
      default=False,
      help="Assert intramolecular cyclization (ring count +1).",
  )
  parser.add_argument(
      "--macrocycle-size",
      type=int,
      default=None,
      help="Assert specific macrocyclic ring size (e.g. 12, 14).",
  )
  parser.add_argument(
      "--strict-hydrogen",
      action="store_true",
      default=False,
      help=(
          "Assert strict hydrogen balance (ΔH = 0). By default, hydrogen"
          " imbalance is reported as a note rather than a failure to"
          " accommodate ionic steps."
      ),
  )
  parser.add_argument(
      "--chiral-catalyst",
      action="store_true",
      default=False,
      help=(
          "Indicate that a chiral catalyst or auxiliary is present for"
          " enantioselective stereogenesis."
      ),
  )
  parser.add_argument(
      "--json",
      action="store_true",
      default=False,
      help="Output results as JSON.",
  )
  parser.add_argument(
      "--output",
      "-o",
      type=str,
      default=None,
      help="Output file path to save stoichiometry report as JSON.",
  )

  args = parser.parse_args()

  res = verify_reaction_stoichiometry(
      precursors=args.precursor,
      products=args.product,
      reaction_type=args.reaction_type,
      byproducts=args.byproducts,
      is_cyclization=args.cyclization or (args.macrocycle_size is not None),
      expected_macrocycle_size=args.macrocycle_size,
      strict_hydrogen=args.strict_hydrogen,
      chiral_catalyst=args.chiral_catalyst,
  )

  if args.output:
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
      json.dump(res, f, indent=2)

  if args.json:
    print(json.dumps(res, indent=2))
  else:
    print("=" * 70)
    print("REACTION ATOM CONSERVATION & STOICHIOMETRY REPORT")
    print("=" * 70)
    print(f"Status:             {'[PASSED]' if res['valid'] else '[FAILED]'}")
    print(f"Precursor Formula:  {res.get('precursor_formula')}")
    print(f"Product Formula:    {res.get('product_formula')}")
    if res.get("applied_byproducts"):
      print(f"Applied Byproducts: {', '.join(res['applied_byproducts'])}")
    print(
        "Carbon Conserved:  "
        f" {'YES' if res['carbon_conserved'] else 'NO (STOLEN / CREATED CARBONS!)'}"
    )
    print(
        f"Heavy Conserved:    {'YES' if res['heavy_atoms_conserved'] else 'NO'}"
    )
    print(
        "Hydrogen Conserved:"
        f" {'YES' if res.get('hydrogen_conserved') else 'NO'}"
    )

    if res.get("notes"):
      for n in res["notes"]:
        print(f"[NOTE] {n}")

    if res.get("delta_atoms"):
      print(f"Net Atom Imbalance: {res['delta_atoms']}")

    ring = res.get("ring_diagnostics", {})
    if ring:
      print(
          f"Rings:              {ring.get('precursor_ring_count')} ->"
          f" {ring.get('product_ring_count')} (Δ ="
          f" {ring.get('delta_rings'):+d})"
      )
      if ring.get("product_macrocycle_ring_sizes"):
        print(
            "Macrocycle Sizes:  "
            f" {ring.get('product_macrocycle_ring_sizes')}-membered rings"
        )

    if res.get("errors"):
      print("\n[ERRORS]:")
      for e in res["errors"]:
        print(f"  ❌ {e}")

    if res.get("warnings"):
      print("\n[WARNINGS]:")
      for w in res["warnings"]:
        print(f"  ⚠️  {w}")
    print("=" * 70)

  sys.exit(0 if res["valid"] else 1)


if __name__ == "__main__":
  main()
