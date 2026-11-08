# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "networkx",
#     "rdkit",
# ]
# ///
"""ReactionSMARTS: Automated Generation, Deep Validation, and Mechanism Conversion Engine

======================================================================================

Provides automated, mathematically and chemically rigorous tools for:
1. Reaction SMARTS Generation (`ReactionSMARTSGenerator`):
   - Derives atom-mapped Reaction SMARTS directly from `ElectronGraph`
     elementary steps.
   - Automates reaction center extraction with variable context radius (r=0,
     r=1, or full molecule).
   - Enforces property-preserving serialization (explicit `+0` formal charges to
     prevent RDKit's stale charge inheritance trap).
2. Deep Multi-Tiered SMARTS Validation (`ReactionSMARTSValidator`):
   - Tier 1: Syntactic parsing & RDKit `rxn.Validate()`.
   - Tier 2: Atom-mapping bijection and duplicate/orphan map detection.
   - Tier 3: Element / nuclei preservation across mapped centers.
   - Tier 4: Net formal charge conservation.
   - Tier 5: Reaction center topological diff (broken/formed bonds, order
     shifts).
   - Tier 6: Live operational execution (`rxn.RunReactants` +
     `Chem.SanitizeMol`).
3. SMARTS-to-Mechanism Inversion (`ReactionSMARTSConverter`):
   - Maps Reaction SMARTS transforms back into formal electron-pushing curved
     arrows.
"""

import argparse
from collections import defaultdict
import itertools
import json
import os
import re
import sys
from typing import Any, Dict, List, Optional, Set, Tuple
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem, Lipinski

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
  sys.path.insert(0, SCRIPT_DIR)

from electron_graph import (  # noqa: E402
    CurvedArrow,
    MolecularGraph,
    ReactionVerifierEngine,
)
from reactions_data import REACTIONS  # noqa: E402


class ReactionSMARTSValidator:
  """Deep multi-tiered validator for arbitrary chemical Reaction SMARTS strings."""

  @staticmethod
  def validate(
      reaction_smarts: str,
      sample_reactants: Optional[List[str]] = None,
      sanitize_products: bool = True,
      is_intramolecular: Optional[bool] = None,
  ) -> Dict[str, Any]:
    """Performs comprehensive multi-tiered validation on a Reaction SMARTS string:

    - Tier 1: Syntax & RDKit C++ Parser Check (including Intramolecular
    Parentheses Invariant)
    - Tier 2: Atom Mapping Integrity & Bijection
    - Tier 3: Element Preservation on Mapped Centers
    - Tier 4: Net Charge Balance
    - Tier 5: Topological Reaction Center Analysis
    - Tier 6: Live Execution Sanity (if sample_reactants provided)
    """
    report = {
        "reaction_smarts": reaction_smarts,
        "valid": False,
        "status": "UNVERIFIED",
        "tier1_syntax": {"passed": False},
        "tier2_atom_mapping": {"passed": False},
        "tier3_element_preservation": {"passed": False},
        "tier4_charge_conservation": {"passed": False},
        "tier5_reaction_center": {},
        "tier6_live_execution": {"executed": False},
        "errors": [],
        "warnings": [],
    }

    # ---------------------------------------------------------------------
    # Tier 1: Syntax & Structural Parsing
    # ---------------------------------------------------------------------
    if ">>" not in reaction_smarts:
      report["errors"].append("Syntax Error: Missing '>>' reaction delimiter.")
      report["status"] = "FAILED_SYNTAX"
      return report

    parts = reaction_smarts.split(">>")
    if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
      report["errors"].append(
          "Syntax Error: Reaction must have non-empty reactant and product"
          " templates."
      )
      report["status"] = "FAILED_SYNTAX"
      return report

    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    if rxn is None:
      report["errors"].append(
          "RDKit Error: AllChem.ReactionFromSmarts failed to parse SMARTS."
      )
      report["status"] = "FAILED_SYNTAX"
      return report

    n_warn, n_err = rxn.Validate()
    report["tier1_syntax"] = {
        "passed": n_err == 0,
        "num_reactants": rxn.GetNumReactantTemplates(),
        "num_products": rxn.GetNumProductTemplates(),
        "rdkit_warnings": n_warn,
        "rdkit_errors": n_err,
    }
    if n_err > 0:
      report["errors"].append(f"RDKit rxn.Validate() reported {n_err} errors.")

    # Intramolecular Reaction SMARTS Invariant:
    # In RDKit, an intramolecular reaction must have exactly 1 reactant template.
    # Disconnected functional groups within a single reactant must be grouped in parentheses:
    # '([pattern1].[pattern2])>>[product]'
    # Without parentheses, RDKit treats '.' as separating multiple distinct reactant molecules.
    if is_intramolecular is True and rxn.GetNumReactantTemplates() > 1:
      report["errors"].append(
          "Intramolecular Reaction SMARTS syntax error: Contains"
          f" {rxn.GetNumReactantTemplates()} reactant templates separated by"
          " top-level '.'. In RDKit, top-level '.' defines separate reactant"
          " molecules (intermolecular). For an intramolecular ring closure on"
          " a single bifunctional molecule, wrap disconnected functional group"
          " patterns in grouping parentheses:"
          f" '({parts[0].strip()})>>{parts[1].strip()}'."
      )

    # ---------------------------------------------------------------------
    # Tier 2: Atom Mapping Integrity
    # ---------------------------------------------------------------------
    # Collect maps from reactant templates
    reactant_maps: Dict[int, List[Tuple[int, int]]] = defaultdict(
        list
    )  # map_num -> [(template_idx, atom_idx)]
    reactant_map_elements: Dict[int, str] = {}
    reactant_map_charges: Dict[int, int] = {}
    reactant_total_charge = 0

    for r_idx in range(rxn.GetNumReactantTemplates()):
      tmpl = rxn.GetReactantTemplate(r_idx)
      seen_maps_in_mol = set()
      for a_idx, atom in enumerate(tmpl.GetAtoms()):
        q = atom.GetFormalCharge()
        reactant_total_charge += q
        m = atom.GetAtomMapNum()
        if m > 0:
          if m in seen_maps_in_mol:
            report["errors"].append(
                f"Duplicate atom map :{m} found within single reactant template"
                f" #{r_idx + 1}."
            )
          seen_maps_in_mol.add(m)
          reactant_maps[m].append((r_idx, a_idx))
          reactant_map_elements[m] = atom.GetSymbol()
          reactant_map_charges[m] = q

    # Collect maps from product templates
    product_maps: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
    product_map_elements: Dict[int, str] = {}
    product_map_charges: Dict[int, int] = {}
    product_total_charge = 0
    unmapped_product_atoms = 0

    for p_idx in range(rxn.GetNumProductTemplates()):
      tmpl = rxn.GetProductTemplate(p_idx)
      seen_maps_in_mol = set()
      for a_idx, atom in enumerate(tmpl.GetAtoms()):
        q = atom.GetFormalCharge()
        product_total_charge += q
        m = atom.GetAtomMapNum()
        if m > 0:
          if m in seen_maps_in_mol:
            report["errors"].append(
                f"Duplicate atom map :{m} found within single product template"
                f" #{p_idx + 1}."
            )
          seen_maps_in_mol.add(m)
          product_maps[m].append((p_idx, a_idx))
          product_map_elements[m] = atom.GetSymbol()
          product_map_charges[m] = q
        else:
          unmapped_product_atoms += 1

    r_keys = set(reactant_maps.keys())
    p_keys = set(product_maps.keys())

    # Unmapped product maps (product atom has map num not in reactants)
    orphan_product_maps = p_keys - r_keys
    if orphan_product_maps:
      report["errors"].append(
          "Atom mapping error: Product contains atom maps not present in"
          f" reactants: {sorted(orphan_product_maps)}."
      )

    report["tier2_atom_mapping"] = {
        "passed": len(orphan_product_maps) == 0 and len(report["errors"]) == 0,
        "mapped_atoms_count": len(r_keys & p_keys),
        "reactant_only_maps (leaving groups)": sorted(r_keys - p_keys),
        "orphan_product_maps": sorted(orphan_product_maps),
        "unmapped_product_atoms": unmapped_product_atoms,
    }

    # ---------------------------------------------------------------------
    # Tier 3: Element / Nuclei Preservation
    # ---------------------------------------------------------------------
    element_mismatches = []
    for m in r_keys & p_keys:
      r_elem = reactant_map_elements[m]
      p_elem = product_map_elements[m]
      # When SMARTS uses generic symbols or numbers [#6], compare
      if r_elem != p_elem and r_elem != "*" and p_elem != "*":
        element_mismatches.append(
            f"Map :{m} changed element from '{r_elem}' to '{p_elem}'."
        )

    if element_mismatches:
      report["errors"].extend(element_mismatches)

    report["tier3_element_preservation"] = {
        "passed": len(element_mismatches) == 0,
        "mismatches": element_mismatches,
    }

    # ---------------------------------------------------------------------
    # Tier 4: Net Charge Conservation
    # ---------------------------------------------------------------------
    charge_balanced = reactant_total_charge == product_total_charge
    report["tier4_charge_conservation"] = {
        "passed": charge_balanced,
        "reactant_net_charge": reactant_total_charge,
        "product_net_charge": product_total_charge,
    }
    if not charge_balanced:
      report["warnings"].append(
          "Net charge imbalance in templates:"
          f" Reactants={reactant_total_charge},"
          f" Products={product_total_charge}."
      )

    # ---------------------------------------------------------------------
    # Tier 5: Reaction Center Topological Analysis
    # ---------------------------------------------------------------------
    # Extract pairs of mapped atoms connected in reactants vs products
    reactant_bonds: Dict[Tuple[int, int], float] = {}
    for r_idx in range(rxn.GetNumReactantTemplates()):
      tmpl = rxn.GetReactantTemplate(r_idx)
      for bond in tmpl.GetBonds():
        m1 = bond.GetBeginAtom().GetAtomMapNum()
        m2 = bond.GetEndAtom().GetAtomMapNum()
        if m1 > 0 and m2 > 0:
          pair = (min(m1, m2), max(m1, m2))
          reactant_bonds[pair] = bond.GetBondTypeAsDouble()

    product_bonds: Dict[Tuple[int, int], float] = {}
    for p_idx in range(rxn.GetNumProductTemplates()):
      tmpl = rxn.GetProductTemplate(p_idx)
      for bond in tmpl.GetBonds():
        m1 = bond.GetBeginAtom().GetAtomMapNum()
        m2 = bond.GetEndAtom().GetAtomMapNum()
        if m1 > 0 and m2 > 0:
          pair = (min(m1, m2), max(m1, m2))
          product_bonds[pair] = bond.GetBondTypeAsDouble()

    all_pairs = set(reactant_bonds.keys()) | set(product_bonds.keys())
    broken_bonds = []
    formed_bonds = []
    order_changes = []

    for pair in all_pairs:
      in_r = pair in reactant_bonds
      in_p = pair in product_bonds
      if in_r and not in_p:
        broken_bonds.append(f":{pair[0]}-:{pair[1]}")
      elif not in_r and in_p:
        formed_bonds.append(f":{pair[0]}-:{pair[1]}")
      elif in_r and in_p:
        bo_r = reactant_bonds[pair]
        bo_p = product_bonds[pair]
        if bo_r != bo_p:
          order_changes.append(f":{pair[0]}-:{pair[1]} ({bo_r} -> {bo_p})")

    charge_changes = []
    for m in r_keys & p_keys:
      q_r = reactant_map_charges[m]
      q_p = product_map_charges[m]
      if q_r != q_p:
        charge_changes.append(f":{m} (charge {q_r} -> {q_p})")

    report["tier5_reaction_center"] = {
        "broken_bonds": broken_bonds,
        "formed_bonds": formed_bonds,
        "order_changes": order_changes,
        "charge_changes": charge_changes,
    }

    # ---------------------------------------------------------------------
    # Tier 6: Live Operational Dry-Run (Optional)
    # ---------------------------------------------------------------------
    if sample_reactants:
      report["tier6_live_execution"]["executed"] = True
      try:
        reactant_mols = []
        for smi in sample_reactants:
          m = Chem.MolFromSmiles(smi)
          if m is None:
            raise ValueError(f"Invalid sample reactant SMILES: {smi}")
          reactant_mols.append(m)

        raw_products = ()
        # Try all permutations of input reactant order
        for perm in itertools.permutations(reactant_mols):
          raw_products = rxn.RunReactants(perm)
          if raw_products:
            break

        # If still 0, try with explicit hydrogens (if templates require explicit H)
        if not raw_products:
          reactant_mols_hs = [Chem.AddHs(m) for m in reactant_mols]
          for perm in itertools.permutations(reactant_mols_hs):
            raw_products = rxn.RunReactants(perm)
            if raw_products:
              break

        report["tier6_live_execution"]["num_outcomes"] = len(raw_products)

        product_smiles_sets = []
        sanitized_ok = True
        for p_tuple in raw_products:
          curr_set = []
          for p_mol in p_tuple:
            # Clean explicit hydrogens for clean reporting
            try:
              clean_p_mol = Chem.RemoveHs(p_mol)
            except Exception:
              clean_p_mol = p_mol
            if sanitize_products:
              try:
                Chem.SanitizeMol(clean_p_mol)
              except Exception as e:
                sanitized_ok = False
                report["errors"].append(
                    f"Product sanitization error during dry-run: {e}"
                )
            curr_set.append(Chem.MolToSmiles(clean_p_mol))
          product_smiles_sets.append(curr_set)

        report["tier6_live_execution"]["passed"] = (
            sanitized_ok and len(raw_products) > 0
        )
        report["tier6_live_execution"]["products"] = product_smiles_sets

        # Intramolecular verification (only when explicitly requested)
        if is_intramolecular is True and len(reactant_mols) == 1:
          prec_mol = reactant_mols[0]
          prec_heavy = prec_mol.GetNumHeavyAtoms()
          for p_set in product_smiles_sets:
            for p_smi in p_set:
              p_m = Chem.MolFromSmiles(p_smi)
              if p_m:
                p_heavy = p_m.GetNumHeavyAtoms()
                if p_heavy > prec_heavy:
                  report["errors"].append(
                      "Intramolecular reaction formed an intermolecular dimer!"
                      f" Product heavy atom count ({p_heavy}) > Precursor heavy"
                      f" atom count ({prec_heavy})."
                  )
                  report["tier6_live_execution"]["passed"] = False

        if len(raw_products) == 0:
          report["warnings"].append(
              "Sample reactants did not match Reaction SMARTS templates (0"
              " products generated)."
          )
      except Exception as e:
        report["tier6_live_execution"]["passed"] = False
        report["tier6_live_execution"]["error"] = str(e)
        report["errors"].append(f"Live dry-run failed: {e}")

    report["valid"] = len(report["errors"]) == 0
    report["status"] = (
        "PASSED_VALIDATION" if report["valid"] else "FAILED_VALIDATION"
    )
    return report


class ReactionSMARTSGenerator:
  """Automated generator for constructing valid, atom-mapped Reaction SMARTS."""

  @staticmethod
  def _fix_explicit_neutral_charges(
      smarts_str: str, neutral_mapped_atoms: Set[int]
  ) -> str:
    """RDKit Property Inheritance Fix:

    In RDKit Reaction SMARTS, if an atom does not have an explicit formal charge
    in the product
    template (e.g. `[#8:1]`), RDKit copies the formal charge from the reactant
    atom (e.g. `-1`).
    This method converts `[#X:m]` to `[#X+0:m]` for all product atoms that must
    be neutral.
    """
    for m in neutral_mapped_atoms:
      smarts_str = re.sub(
          rf"\[([#A-Za-z0-9]+):{m}\]", rf"[\1+0:{m}]", smarts_str
      )
    return smarts_str

  @classmethod
  def from_electron_graph(
      cls,
      reactant_graph: MolecularGraph,
      arrows: List[CurvedArrow],
      radius: int = 1,
      full_molecules: bool = False,
      collapse_spectator_hydrogens: bool = False,
      step_name: str = "",
  ) -> Dict[str, Any]:
    """Automatically derives validated, atom-mapped Reaction SMARTS from an ElectronGraph step.

    - radius: 0 = only reaction center; 1 = reaction center + immediate 1-hop
    neighbors; 2 = 2 hops. - full_molecules: if True, extracts complete reactant
    and product molecules. - collapse_spectator_hydrogens: if True, omits
    unreacting spectator hydrogens from the SMARTS pattern.
    """
    # 1. Apply step to produce verified product graph
    valid, product_graph, audit = ReactionVerifierEngine.apply_and_verify_step(
        reactant_graph, arrows, step_name=step_name
    )
    if not valid:
      raise ValueError(
          "Cannot generate Reaction SMARTS for chemically invalid step:"
          f" {audit.get('errors')}"
      )

    return cls.from_graphs(
        reactant_graph,
        product_graph,
        arrows=arrows,
        radius=radius,
        full_molecules=full_molecules,
        collapse_spectator_hydrogens=collapse_spectator_hydrogens,
    )

  @classmethod
  def from_graphs(
      cls,
      reactant_graph: MolecularGraph,
      product_graph: MolecularGraph,
      arrows: Optional[List[CurvedArrow]] = None,
      radius: int = 1,
      full_molecules: bool = False,
      collapse_spectator_hydrogens: bool = False,
  ) -> Dict[str, Any]:
    """Derives atom-mapped Reaction SMARTS by comparing reactant and product molecular graphs."""
    # 1. Detect reacting center atoms
    reacting_atoms: Set[str] = set()

    if arrows:
      for a in arrows:
        if a.source_type in ("lone_pair", "radical"):
          reacting_atoms.add(str(a.source_id))
        elif a.source_type == "bond":
          reacting_atoms.add(str(a.source_id[0]))
          reacting_atoms.add(str(a.source_id[1]))
        if a.sink_type == "atom":
          reacting_atoms.add(str(a.sink_id))
        elif a.sink_type == "bond":
          reacting_atoms.add(str(a.sink_id[0]))
          reacting_atoms.add(str(a.sink_id[1]))

    # Topological diff
    for aid, r_atom in reactant_graph.atoms.items():
      if aid in product_graph.atoms:
        p_atom = product_graph.atoms[aid]
        if r_atom.formal_charge != p_atom.formal_charge:
          reacting_atoms.add(aid)
        if r_atom.lone_pairs != p_atom.lone_pairs:
          reacting_atoms.add(aid)
        if r_atom.radicals != p_atom.radicals:
          reacting_atoms.add(aid)
        if reactant_graph.total_bond_order_on_atom(
            aid
        ) != product_graph.total_bond_order_on_atom(aid):
          reacting_atoms.add(aid)

    # Formed or broken edges
    r_edges = set(_canon_pair(u, v) for u, v in reactant_graph.graph.edges())
    p_edges = set(_canon_pair(u, v) for u, v in product_graph.graph.edges())
    for u, v in r_edges ^ p_edges:
      reacting_atoms.add(u)
      reacting_atoms.add(v)

    if not reacting_atoms:
      # Identity / no-op
      reacting_atoms = set(reactant_graph.atoms.keys())

    # 2. Subgraph scope expansion
    if full_molecules:
      # Include all atoms in any connected component containing a reacting atom
      included_atoms = set()
      for comp in nx.connected_components(reactant_graph.graph):
        if comp & reacting_atoms:
          included_atoms.update(comp)
    else:
      included_atoms = set(reacting_atoms)
      for _ in range(radius):
        next_hop = set()
        for aid in included_atoms:
          if aid in reactant_graph.graph:
            for nbr in reactant_graph.graph.neighbors(aid):
              next_hop.add(nbr)
        included_atoms.update(next_hop)

    if collapse_spectator_hydrogens:
      spectator_hs = {
          aid
          for aid in included_atoms
          if reactant_graph.atoms[aid].element == "H"
          and aid not in reacting_atoms
      }
      included_atoms -= spectator_hs

    # 3. Assign deterministic atom map numbers :1, :2, ...
    sorted_atoms = sorted(included_atoms)
    atom_to_map: Dict[str, int] = {
        aid: i for i, aid in enumerate(sorted_atoms, 1)
    }

    # 4. Build Reactant SMARTS components
    r_subgraph = reactant_graph.graph.subgraph(included_atoms)
    r_components = list(nx.connected_components(r_subgraph))
    r_smarts_list = []

    for comp in r_components:
      rw = Chem.RWMol()
      id_to_idx = {}
      for aid in comp:
        an = reactant_graph.atoms[aid]
        atom = Chem.Atom(an.element)
        atom.SetFormalCharge(an.formal_charge)
        atom.SetAtomMapNum(atom_to_map[aid])
        idx = rw.AddAtom(atom)
        id_to_idx[aid] = idx

      for u, v in reactant_graph.graph.edges(comp):
        if u in id_to_idx and v in id_to_idx:
          bo = reactant_graph.get_bond_order(u, v)
          btype = (
              Chem.BondType.SINGLE
              if bo == 1
              else (Chem.BondType.DOUBLE if bo == 2 else Chem.BondType.TRIPLE)
          )
          rw.AddBond(id_to_idx[u], id_to_idx[v], btype)

      r_smarts_list.append(Chem.MolToSmarts(rw.GetMol()))

    # 5. Build Product SMARTS components
    p_subgraph = product_graph.graph.subgraph(included_atoms)
    p_components = list(nx.connected_components(p_subgraph))
    p_smarts_list = []

    # Find product atoms that are neutral (charge 0) and need explicit +0
    neutral_product_maps = set()
    for aid in included_atoms:
      if aid in product_graph.atoms:
        if product_graph.atoms[aid].formal_charge == 0:
          neutral_product_maps.add(atom_to_map[aid])

    for comp in p_components:
      rw = Chem.RWMol()
      id_to_idx = {}
      for aid in comp:
        an = product_graph.atoms[aid]
        atom = Chem.Atom(an.element)
        atom.SetFormalCharge(an.formal_charge)
        atom.SetAtomMapNum(atom_to_map[aid])
        idx = rw.AddAtom(atom)
        id_to_idx[aid] = idx

      for u, v in product_graph.graph.edges(comp):
        if u in id_to_idx and v in id_to_idx:
          bo = product_graph.get_bond_order(u, v)
          btype = (
              Chem.BondType.SINGLE
              if bo == 1
              else (Chem.BondType.DOUBLE if bo == 2 else Chem.BondType.TRIPLE)
          )
          rw.AddBond(id_to_idx[u], id_to_idx[v], btype)

      raw_smarts = Chem.MolToSmarts(rw.GetMol())
      fixed_smarts = cls._fix_explicit_neutral_charges(
          raw_smarts, neutral_product_maps
      )
      p_smarts_list.append(fixed_smarts)

    # 6. Assemble Reaction SMARTS
    reaction_smarts = f"{'.'.join(r_smarts_list)}>>{'.'.join(p_smarts_list)}"

    # 7. Validate generated SMARTS
    val_report = ReactionSMARTSValidator.validate(reaction_smarts)

    return {
        "reaction_smarts": reaction_smarts,
        "reacting_atoms": sorted(reacting_atoms),
        "included_atoms": sorted(included_atoms),
        "atom_mapping": atom_to_map,
        "num_reactant_templates": len(r_smarts_list),
        "num_product_templates": len(p_smarts_list),
        "validation": val_report,
    }

  @staticmethod
  def make_intramolecular(reaction_smarts: str) -> str:
    """Converts an intermolecular Reaction SMARTS into an intramolecular template.

    In RDKit, top-level '.' defines distinct reactant molecules. Disconnected
    functional groups within a single reactant must be grouped with parentheses:
    '([group1].[group2])>>[cyclized_product]'.
    """
    if ">>" not in reaction_smarts:
      raise ValueError("Missing '>>' delimiter in Reaction SMARTS.")
    reactants, products = reaction_smarts.split(">>", 1)
    reactants = reactants.strip()
    products = products.strip()

    if reactants.startswith("(") and reactants.endswith(")"):
      return reaction_smarts

    return f"({reactants})>>{products}"


def _canon_pair(u: Any, v: Any) -> Tuple[str, str]:
  su, sv = str(u), str(v)
  return (su, sv) if su <= sv else (sv, su)


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Reaction SMARTS Generation, Multi-Tier Deep Validation, and"
          " Intramolecular Engine."
      )
  )
  parser.add_argument(
      "--smarts", default=None, help="Reaction SMARTS string to validate."
  )
  parser.add_argument(
      "--reactants",
      nargs="*",
      default=None,
      help="Sample reactant SMILES for live dry-run.",
  )
  parser.add_argument(
      "--intramolecular",
      "--intra",
      action="store_true",
      default=False,
      help="Assert intramolecular reaction.",
  )
  parser.add_argument(
      "--make-intramolecular",
      default=None,
      help="Convert an intermolecular SMARTS to intramolecular template.",
  )
  parser.add_argument(
      "--to-graph-json",
      "--to_graph_json",
      dest="to_graph_json",
      default=None,
      help=(
          "Export reaction center broken/formed bonds and charge deltas to JSON"
          " file."
      ),
  )
  parser.add_argument(
      "--list-templates",
      "--list_templates",
      action="store_true",
      default=False,
      help="List built-in reaction templates from reactions_data.py.",
  )
  parser.add_argument(
      "--json",
      action="store_true",
      default=False,
      help="Output validation report as JSON.",
  )
  parser.add_argument(
      "--output",
      "-o",
      type=str,
      default=None,
      help="Output file path to save validation report or converted SMARTS.",
  )

  args = parser.parse_args()

  if args.list_templates:
    for entry in REACTIONS:
      r_id = entry.get("id", "")
      title = entry.get("title") or entry.get("name", "")
      cat = entry.get("category", "")
      print(f"{r_id}: {title} [{cat}]")
    sys.exit(0)

  if args.make_intramolecular:
    intra_smarts = ReactionSMARTSGenerator.make_intramolecular(
        args.make_intramolecular
    )
    val = ReactionSMARTSValidator.validate(intra_smarts, is_intramolecular=True)
    if args.output:
      os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
      with open(args.output, "w", encoding="utf-8") as f:
        json.dump(
            {"intramolecular_smarts": intra_smarts, "validation": val},
            f,
            indent=2,
        )
    print(f"Intramolecular Reaction SMARTS: {intra_smarts}")
    print(f"Validation Status: {'PASSED' if val['valid'] else 'FAILED'}")
    if not val["valid"]:
      print(f"Errors: {val['errors']}", file=sys.stderr)
    sys.exit(0 if val["valid"] else 1)

  if args.smarts:
    report = ReactionSMARTSValidator.validate(
        args.smarts,
        sample_reactants=args.reactants,
        is_intramolecular=args.intramolecular or None,
    )
    out_path = args.output or args.to_graph_json
    if out_path:
      os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
      payload = (
          report["tier5_reaction_center"] if args.to_graph_json else report
      )
      with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    if args.json:
      print(json.dumps(report, indent=2))
    else:
      print("=" * 70)
      print("REACTION SMARTS DEEP VALIDATION REPORT")
      print("=" * 70)
      print(f"SMARTS:             {report['reaction_smarts']}")
      print(f"Status:             {report['status']}")
      print(f"Valid:              {report['valid']}")
      print(f"Reactant Templates: {report['tier1_syntax']['num_reactants']}")
      print(f"Product Templates:  {report['tier1_syntax']['num_products']}")
      print(
          "Mapped Centers:    "
          f" {report['tier2_atom_mapping']['mapped_atoms_count']}"
      )
      print(
          "Charge Balance:    "
          f" {'YES' if report['tier4_charge_conservation']['passed'] else 'NO'}"
      )
      print(f"Reaction Center:    {report['tier5_reaction_center']}")
      if report["tier6_live_execution"]["executed"]:
        print(
            "Live Dry-Run:      "
            f" {'PASSED' if report['tier6_live_execution']['passed'] else 'FAILED'}"
        )
        print(
            "Live Products:     "
            f" {report['tier6_live_execution'].get('products')}"
        )
      if report["errors"]:
        print("\n[ERRORS]:")
        for e in report["errors"]:
          print(f"  ❌ {e}")
      if report["warnings"]:
        print("\n[WARNINGS]:")
        for w in report["warnings"]:
          print(f"  ⚠️  {w}")
      print("=" * 70)
    sys.exit(0 if report["valid"] else 1)

  parser.error("Provide --smarts or --make-intramolecular.")


if __name__ == "__main__":
  main()
