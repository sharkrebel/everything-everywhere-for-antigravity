# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "networkx",
# ]
# ///
"""ElectronGraph: Formal Graph-Theoretic Reaction Verification & Electron-Pushing Engine

=====================================================================================

A chemically rigorous graph engine for electron-pushing reaction mechanisms.
Formal Definition:
  A chemical state is a multigraph G = (V, E, L, R, Q) where:
    - V is the set of atomic nuclei with element Z(v), formal charge Q(v), lone
      pairs L(v), radical electrons R(v), and implicit hydrogens H_impl(v).
    - E is the set of bonds with bond order BO(u, v) in {1, 2, 3}.
    - Valence electrons N_val(v) = 2 * (sum_{u ~ v} BO(u, v) + H_impl(v)) + 2 *
      L(v) + R(v).
    - Valence capacity N_cap(v): H/He -> 2, 2nd row (B, C, N, O, F) -> 8 (strict
      octet), 3rd row and heavier (Si, P, S, Cl, Br, I) -> 12 (expanded octet),
      metals (Li, Na, K, Mg, Zn, etc.) -> 8 or 18 (covalent/coordination
      models).

An Elementary Step is an ordered set of curved electron-pushing arrows:
    A = {a_1, a_2, ..., a_k} where each arrow a_i represents a transfer of n
    electrons (n=2 for double-barbed curved arrows,
n=1 for single-barbed radical fishhooks) from a Source to a Sink:
    Source in { LonePair(u), Radical(u), Bond(u, v) } Sink in { Atom(w), Bond(w,
    x) }

Chemical Invariants & Physical Principles Enforced:
  1. Topological Locality Invariant:
     - Bond(u, v) -> Atom(w): w MUST be in {u, v} (bond cleaves onto one of its
       own atoms).
     - LonePair(u) / Radical(u) -> Bond(w, x): u MUST be in {w, x} (bond forms
       with source atom).
     - Bond(u, v) -> Bond(w, x): |{u, v} ∩ {w, x}| >= 1 (bond shift must share a
       pivot atom).
     - LonePair(u) -> Atom(w): w == u (electrons cannot teleport between
       unconnected atoms).
  2. Deterministic Differential Formal Charge Conservation:
     Delta Q(u) = - (2 * Delta L(u) + Delta R(u) + Delta BO(u)) Guarantees exact
     formal charge updates on reacting atoms while leaving unreacting spectator
     atoms invariant and conserving global charge.
  3. Total Electron Conservation:
     Initial valence electrons == Final valence electrons.
  4. Global Net Charge Conservation:
     Initial net charge == Final net charge.
  5. Nuclei / Mass Conservation:
     Multiset of atomic nuclei (elements and counts) is strictly invariant
     across the step.
  6. Lewis Octet Rule & Valence Capacity:
     Strict non-hypervalency on 2nd row atoms (C, N, O, F <= 8 e-).
  7. Spectator Locality Invariant:
     All atoms not explicitly involved as arrow sources or sinks strictly
     maintain their incident edge list (neighbors and bond orders), formal
     charges, lone pairs, and radicals.
  8. Non-negative Lone Pairs, Radicals, and Positive Bond Orders.
  9. Bidirectional RDKit Interoperability:
     - from_smiles(smiles, make_hs_explicit=False) -> MolecularGraph
     - to_rdkit(sanitize=True) -> Chem.Mol
     - to_smiles() -> str
     - verify_product_match(product_graph, target_smiles_or_mol) -> (bool, str)
     - verify_step_against_target(reactant_graph, arrows, target_product) ->
       (bool, MolecularGraph, dict)
"""

import argparse
from collections import Counter
import copy
from dataclasses import dataclass, field
import itertools
import json
import math
import sys
from typing import Any, Dict, FrozenSet, List, Optional, Sequence, Tuple, Union

import networkx as nx
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

# Standard valence electrons in neutral free atoms (Group valence)
VALENCE_ELECTRONS: Dict[str, int] = {
    "H": 1,
    "He": 2,
    "Li": 1,
    "Be": 2,
    "B": 3,
    "C": 4,
    "N": 5,
    "O": 6,
    "F": 7,
    "Ne": 8,
    "Na": 1,
    "Mg": 2,
    "Al": 3,
    "Si": 4,
    "P": 5,
    "S": 6,
    "Cl": 7,
    "Ar": 8,
    "K": 1,
    "Ca": 2,
    "Sc": 3,
    "Ti": 4,
    "V": 5,
    "Cr": 6,
    "Mn": 7,
    "Fe": 8,
    "Co": 9,
    "Ni": 10,
    "Cu": 11,
    "Zn": 2,
    "Ga": 3,
    "Ge": 4,
    "As": 5,
    "Se": 6,
    "Br": 7,
    "Kr": 8,
    "Rb": 1,
    "Sr": 2,
    "Y": 3,
    "Zr": 4,
    "Nb": 5,
    "Mo": 6,
    "Tc": 7,
    "Ru": 8,
    "Rh": 9,
    "Pd": 10,
    "Ag": 11,
    "Cd": 2,
    "In": 3,
    "Sn": 4,
    "Sb": 5,
    "Te": 6,
    "I": 7,
    "Xe": 8,
    "Cs": 1,
    "Ba": 2,
    "La": 3,
    "Ce": 3,
    "Hf": 4,
    "Ta": 5,
    "W": 6,
    "Re": 7,
    "Os": 8,
    "Ir": 9,
    "Pt": 10,
    "Au": 11,
    "Hg": 4,
    "Pb": 4,
    "Bi": 5,
}

# Pauling electronegativities for canonical dipole-directed arrow synthesis
PAULING_ELECTRONEGATIVITY: Dict[str, float] = {
    "H": 2.20,
    "Li": 0.98,
    "B": 2.04,
    "C": 2.55,
    "N": 3.04,
    "O": 3.44,
    "F": 3.98,
    "Na": 0.93,
    "Mg": 1.31,
    "Si": 1.90,
    "P": 2.19,
    "S": 2.58,
    "Cl": 3.16,
    "K": 0.82,
    "Cu": 1.90,
    "Zn": 1.65,
    "Br": 2.96,
    "Pd": 2.20,
    "I": 2.66,
    "Hg": 2.00,
}

# Maximum valence capacity (maximum electrons in valence shell)
# 1st row: H (2 duet), He (2)
# 2nd row: B (8), C (8), N (8), O (8), F (8) -> strict octet
# 3rd row and higher: Si (12), P (12), S (12), Cl (12), Br (12), I (12) -> hypervalent allowed
# Metals in organometallic / coordination models: Li (8), Na (8), K (8), Mg (8), Zn (18)
MAX_VALENCE_CAPACITY: Dict[str, int] = {
    "H": 2,
    "He": 2,
    "Li": 8,
    "Be": 8,
    "B": 8,
    "C": 8,
    "N": 8,
    "O": 8,
    "F": 8,
    "Ne": 8,
    "Na": 8,
    "Mg": 8,
    "Al": 8,
    "Si": 12,
    "P": 12,
    "S": 12,
    "Cl": 12,
    "Ar": 8,
    "K": 8,
    "Ca": 8,
    "Sc": 18,
    "Ti": 18,
    "V": 18,
    "Cr": 18,
    "Mn": 18,
    "Fe": 18,
    "Co": 18,
    "Ni": 18,
    "Cu": 18,
    "Zn": 18,
    "Ga": 18,
    "Ge": 18,
    "As": 12,
    "Se": 12,
    "Br": 12,
    "Kr": 8,
    "Rb": 8,
    "Sr": 8,
    "Y": 18,
    "Zr": 18,
    "Nb": 18,
    "Mo": 18,
    "Tc": 18,
    "Ru": 18,
    "Rh": 18,
    "Pd": 18,
    "Ag": 18,
    "Cd": 18,
    "Sn": 18,
    "Sb": 12,
    "Te": 12,
    "I": 12,
    "Xe": 12,
    "La": 18,
    "Ce": 18,
    "Hf": 18,
    "Ta": 18,
    "W": 18,
    "Re": 18,
    "Os": 18,
    "Ir": 18,
    "Pt": 18,
    "Au": 18,
    "Hg": 18,
}

# Open-shell d/f-block metals supporting unpaired d/f-electrons and sigma bonds
D_F_BLOCK_METALS: FrozenSet[str] = frozenset({
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "La",
    "Ce",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
})


def _canon_pair(u: Any, v: Any) -> Tuple[str, str]:
  """Returns a canonical string tuple representing an undirected edge."""
  su, sv = str(u), str(v)
  return (su, sv) if su <= sv else (sv, su)


@dataclass
class AtomNode:
  """Represents a single atomic nucleus and its immediate electronic state."""

  id: str
  element: str
  lone_pairs: int = 0
  radicals: int = 0
  unpaired_d: int = 0
  formal_charge: int = 0
  implicit_h: int = 0
  chiral_tag: int = (
      0  # 0=CHI_UNSPECIFIED, 1=CHI_TETRAHEDRAL_CW, 2=CHI_TETRAHEDRAL_CCW
  )
  chiral_neighbors: List[str] = field(default_factory=list)
  metadata: Dict[str, Any] = field(default_factory=dict)

  def total_valence_electrons(self, total_bond_order: float) -> int:
    """Total electrons currently in the valence shell of this atom."""
    return int(
        2 * (total_bond_order + self.implicit_h)
        + 2 * self.lone_pairs
        + self.radicals
        + self.unpaired_d
    )

  def calculate_formal_charge(self, total_bond_order: float) -> int:
    """Calculates formal charge deterministically as a pure state function:

    FC = V_free - 2*L - radicals - unpaired_d - (total_bond_order + implicit_h)
    """
    v_free = VALENCE_ELECTRONS.get(self.element, 4)
    v_bonded = int(
        total_bond_order
        + self.implicit_h
        + 2 * self.lone_pairs
        + self.radicals
        + self.unpaired_d
    )
    return v_free - v_bonded


@dataclass
class CurvedArrow:
  """Represents an electron-pushing arrow in an elementary reaction step.

  - electrons=2: standard curved double-barbed arrow (electron pair) -
  electrons=1: single-barbed radical fishhook arrow (single electron)
  """

  source_type: str  # 'lone_pair', 'radical', or 'bond'
  source_id: Any  # atom_id (str) if lone_pair/radical, or (u, v) tuple if bond
  sink_type: str  # 'atom' or 'bond'
  sink_id: Any  # atom_id (str) if atom, or (w, x) tuple if bond
  electrons: int = 2  # 2 = pair, 1 = single electron
  label: str = ""
  stereo_policy: Optional[str] = (
      "auto"  # 'auto', 'invert', 'retain', 'racemize'
  )
  is_dative: Optional[bool] = None

  def __post_init__(self):
    if self.source_type == "bond" and isinstance(self.source_id, (list, tuple)):
      self.source_id = (str(self.source_id[0]), str(self.source_id[1]))
    elif isinstance(self.source_id, (int, str)):
      self.source_id = str(self.source_id)

    if self.sink_type == "bond" and isinstance(self.sink_id, (list, tuple)):
      self.sink_id = (str(self.sink_id[0]), str(self.sink_id[1]))
    elif isinstance(self.sink_id, (int, str)):
      self.sink_id = str(self.sink_id)

  def __repr__(self):
    arrow_sym = "↷ (2e⁻)" if self.electrons == 2 else "⇀ (1e⁻)"
    if self.source_type == "lone_pair":
      src = f"LP({self.source_id})"
    elif self.source_type == "radical":
      src = f"Rad({self.source_id})"
    else:
      src = f"Bond({self.source_id[0]}-{self.source_id[1]})"

    if self.sink_type == "atom":
      snk = f"Atom({self.sink_id})"
    else:
      snk = f"Bond({self.sink_id[0]}-{self.sink_id[1]})"

    lbl = f" '{self.label}'" if self.label else ""
    return f"CurvedArrow({src} {arrow_sym} {snk}{lbl})"


def normalize_dipolar_rdkit_mol(mol: Chem.Mol) -> Chem.Mol:
  """Normalizes pentavalent 2nd-row nitrogen dipolar groups (nitro, N-oxide, azide, diazo, isonitrile) in place.

  Preserves all fragments, explicit hydrogens, and atom indices (0..N-1) while
  converting un-normalized pentavalent nitrogen representations into strict
  Lewis octet-compliant zwitterionic states.
  """
  if mol is None or mol.GetNumAtoms() == 0:
    return mol

  mol_work = Chem.Mol(mol)
  mol_work.UpdatePropertyCache(strict=False)

  def _bo_sum(atom: Chem.Atom) -> float:
    total = 0.0
    for b in atom.GetBonds():
      bt = b.GetBondType()
      if bt == Chem.BondType.SINGLE:
        total += 1.0
      elif bt == Chem.BondType.DOUBLE:
        total += 2.0
      elif bt == Chem.BondType.TRIPLE:
        total += 3.0
      elif bt == Chem.BondType.AROMATIC:
        total += 1.5
      else:
        total += b.GetBondTypeAsDouble()
    return total

  needs_norm = any(
      a.GetSymbol() == "N"
      and a.GetFormalCharge() == 0
      and round(_bo_sum(a) + a.GetNumExplicitHs() + a.GetNumImplicitHs(), 2)
      >= 4.5
      for a in mol_work.GetAtoms()
  )
  if not needs_norm:
    return mol

  rw = Chem.RWMol(mol_work)
  for idx in range(rw.GetNumAtoms()):
    atom = rw.GetAtomWithIdx(idx)
    if atom.GetSymbol() != "N" or atom.GetFormalCharge() != 0:
      continue
    total_val = round(
        _bo_sum(atom) + atom.GetNumExplicitHs() + atom.GetNumImplicitHs(), 2
    )
    if total_val < 4.5:
      continue

    double_targets: List[Tuple[int, int, Chem.Atom]] = []
    triple_targets: List[Tuple[int, Chem.Atom]] = []
    for b in atom.GetBonds():
      nb = b.GetOtherAtom(atom)
      nb_idx = nb.GetIdx()
      if b.GetBondType() == Chem.BondType.DOUBLE and nb.GetFormalCharge() == 0:
        prio = (
            0
            if (nb.GetSymbol() == "O" and nb.GetDegree() == 1)
            else (1 if nb.GetSymbol() in ("O", "S") else 2)
        )
        double_targets.append((prio, nb_idx, nb))
      elif (
          b.GetBondType() == Chem.BondType.TRIPLE and nb.GetFormalCharge() == 0
      ):
        triple_targets.append((nb_idx, nb))

    if double_targets and not triple_targets:
      double_targets.sort(key=lambda t: (t[0], t[1]))
      _, nb_idx, nb = double_targets[0]
      bond = rw.GetBondBetweenAtoms(idx, nb_idx)
      bond.SetBondType(Chem.BondType.SINGLE)
      atom.SetFormalCharge(1)
      nb.SetFormalCharge(-1)
      nb.SetNoImplicit(True)
      nb.SetNumExplicitHs(0)
    elif triple_targets:
      nb_idx, nb = triple_targets[0]
      bond = rw.GetBondBetweenAtoms(idx, nb_idx)
      if nb.GetSymbol() == "N" and nb.GetDegree() == 1:
        bond.SetBondType(Chem.BondType.DOUBLE)
        atom.SetFormalCharge(1)
        nb.SetFormalCharge(-1)
        nb.SetNoImplicit(True)
        nb.SetNumExplicitHs(0)
      elif nb.GetSymbol() == "C" and nb.GetDegree() == 1:
        atom.SetFormalCharge(1)
        nb.SetFormalCharge(-1)
        nb.SetNoImplicit(True)
        nb.SetNumExplicitHs(0)

  norm_mol = rw.GetMol()
  try:
    Chem.SanitizeMol(norm_mol)
  except Exception:
    norm_mol.UpdatePropertyCache(strict=False)
  return norm_mol


class MolecularGraph:
  """Rigorous chemical molecular graph data structure supporting electron-pushing rewrites,

  chemical invariant audits, and seamless bidirectional RDKit conversion.
  """

  def __init__(self, name: str = ""):
    self.name = name
    self.graph = nx.Graph()
    self.atoms: Dict[str, AtomNode] = {}

  def add_atom(
      self,
      atom_id: str,
      element: str,
      lone_pairs: int = 0,
      radicals: int = 0,
      unpaired_d: int = 0,
      formal_charge: Optional[int] = None,
      implicit_h: int = 0,
      chiral_tag: int = 0,
      chiral_neighbors: Optional[List[str]] = None,
      **meta,
  ) -> "MolecularGraph":
    aid = str(atom_id)
    node = AtomNode(
        id=aid,
        element=element,
        lone_pairs=int(lone_pairs),
        radicals=int(radicals),
        unpaired_d=int(unpaired_d),
        formal_charge=int(formal_charge) if formal_charge is not None else 0,
        implicit_h=int(implicit_h),
        chiral_tag=int(chiral_tag),
        chiral_neighbors=(
            list(chiral_neighbors) if chiral_neighbors is not None else []
        ),
        metadata=dict(meta),
    )
    self.atoms[aid] = node
    self.graph.add_node(aid, data=node)
    return self

  def add_bond(
      self,
      u: str,
      v: str,
      order: int = 1,
      stereo: str = "NONE",
      stereo_atoms: Optional[Sequence[str]] = None,
      is_dative: bool = False,
      dative_begin: Optional[str] = None,
      dative_end: Optional[str] = None,
  ) -> "MolecularGraph":
    su, sv = str(u), str(v)
    if su not in self.atoms or sv not in self.atoms:
      raise KeyError(
          f"Both atoms '{su}' and '{sv}' must exist in graph before creating"
          " bond."
      )
    sa = tuple(str(x) for x in stereo_atoms) if stereo_atoms else ()
    self.graph.add_edge(
        su,
        sv,
        order=int(order),
        stereo=str(stereo),
        stereo_atoms=sa,
        is_dative=bool(is_dative),
        dative_begin=str(dative_begin) if dative_begin is not None else None,
        dative_end=str(dative_end) if dative_end is not None else None,
    )
    return self

  def remove_bond(self, u: str, v: str) -> "MolecularGraph":
    su, sv = str(u), str(v)
    if self.graph.has_edge(su, sv):
      self.graph.remove_edge(su, sv)
    return self

  def get_bond_order(self, u: str, v: str) -> int:
    su, sv = str(u), str(v)
    if self.graph.has_edge(su, sv):
      return int(self.graph[su][sv].get("order", 1))
    return 0

  def total_bond_order_on_atom(self, u: str) -> int:
    su = str(u)
    if su not in self.graph:
      return 0
    return sum(
        int(self.graph[su][v].get("order", 1)) for v in self.graph.neighbors(su)
    )

  def effective_bond_valence_on_atom(self, u: str) -> int:
    """Bonding valence for formal-charge accounting (dative donor +2, acceptor 0)."""
    su = str(u)
    if su not in self.graph:
      return 0
    total = 0
    for v in self.graph.neighbors(su):
      edata = self.graph[su][v]
      order = int(edata.get("order", 1))
      if edata.get("is_dative"):
        if edata.get("dative_begin") == su:
          total += 2 + max(0, order - 1)
        else:
          total += max(0, order - 1)
      else:
        total += order
    return total

  def recompute_formal_charges(self) -> "MolecularGraph":
    """Deterministically recomputes formal charges for all atoms via Lewis state function."""
    for aid, atom in self.atoms.items():
      total_bo = self.effective_bond_valence_on_atom(aid)
      atom.formal_charge = atom.calculate_formal_charge(total_bo)
    return self

  def net_charge(self) -> int:
    """Computes net formal charge across all atoms."""
    return sum(a.formal_charge for a in self.atoms.values())

  def total_electrons(self) -> int:
    """Computes total valence electrons in the entire molecular system:

    E_total = 2 * (Sum of explicit bond orders + Sum of implicit H bonds)
              + 2 * (Sum of lone pairs)
              + Sum of radical electrons
              + Sum of unpaired d/f-shell electrons on metals
    """
    bond_e = 2 * sum(
        int(data.get("order", 1)) for _, _, data in self.graph.edges(data=True)
    )
    impl_h_e = 2 * sum(an.implicit_h for an in self.atoms.values())
    lp_e = 2 * sum(an.lone_pairs for an in self.atoms.values())
    rad_e = sum(an.radicals + an.unpaired_d for an in self.atoms.values())
    return int(bond_e + impl_h_e + lp_e + rad_e)

  def nuclei_multiset(self) -> Counter:
    """Returns element stoichiometry (including implicit hydrogens) for mass conservation."""
    counts: Counter = Counter()
    for an in self.atoms.values():
      counts[an.element] += 1
      if an.implicit_h > 0:
        counts["H"] += an.implicit_h
    return counts

  def copy(self) -> "MolecularGraph":
    """Creates a deep copy of this MolecularGraph."""
    return copy.deepcopy(self)

  def to_summary_dict(self) -> Dict[str, Any]:
    """Serializes current state for audit logs."""
    return {
        "name": self.name,
        "total_valence_electrons": self.total_electrons(),
        "net_charge": self.net_charge(),
        "nuclei": dict(self.nuclei_multiset()),
        "atoms": {
            aid: {
                "element": an.element,
                "lone_pairs": an.lone_pairs,
                "radicals": an.radicals,
                "unpaired_d": an.unpaired_d,
                "formal_charge": an.formal_charge,
                "implicit_h": an.implicit_h,
                "total_bond_order": self.total_bond_order_on_atom(aid),
                "valence_shell_e": an.total_valence_electrons(
                    self.total_bond_order_on_atom(aid)
                ),
            }
            for aid, an in self.atoms.items()
        },
        "bonds": [
            {
                "u": u,
                "v": v,
                "order": int(d.get("order", 1)),
                "stereo": d.get("stereo", "NONE"),
                "is_dative": bool(d.get("is_dative", False)),
            }
            for u, v, d in self.graph.edges(data=True)
        ],
    }

  # -------------------------------------------------------------------------
  # Bidirectional RDKit Interoperability
  # -------------------------------------------------------------------------
  @classmethod
  def from_smiles(
      cls,
      smiles: str,
      name: str = "",
      make_hs_explicit: bool = False,
      strict_lewis: bool = False,
  ) -> "MolecularGraph":
    """Constructs a MolecularGraph directly from a SMILES string via RDKit.

    Automatically Kekulizes aromatic systems to explicit single/double bonds,
    normalizes pentavalent 2nd-row nitrogen dipolar representations, and infers
    lone pairs from standard group valence and formal charges.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
      mol_raw = Chem.MolFromSmiles(smiles, sanitize=False)
      if mol_raw is not None and mol_raw.GetNumAtoms() > 0:
        try:
          mol_raw.UpdatePropertyCache(strict=False)
          mol_norm = normalize_dipolar_rdkit_mol(mol_raw)
          Chem.SanitizeMol(mol_norm)
          mol = mol_norm
        except Exception:
          mol = None
    if mol is None:
      raise ValueError(f"RDKit failed to parse SMILES: '{smiles}'")
    return cls.from_rdkit(
        mol,
        name=name or smiles,
        make_hs_explicit=make_hs_explicit,
        strict_lewis=strict_lewis,
    )

  @classmethod
  def from_rdkit(
      cls,
      mol: Chem.Mol,
      name: str = "",
      make_hs_explicit: bool = False,
      strict_lewis: bool = False,
  ) -> "MolecularGraph":
    """Constructs a MolecularGraph from an RDKit Mol object."""
    mol = normalize_dipolar_rdkit_mol(Chem.Mol(mol))
    if make_hs_explicit:
      mol = Chem.AddHs(mol)
    try:
      Chem.AssignStereochemistry(mol, cleanIt=False, force=True)
    except Exception:
      pass
    try:
      Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
      pass

    mg = cls(name=name)
    for atom in mol.GetAtoms():
      idx_str = str(atom.GetIdx())
      elem = atom.GetSymbol()
      q = atom.GetFormalCharge()
      rad = 0 if elem in D_F_BLOCK_METALS else atom.GetNumRadicalElectrons()
      impl_h = 0 if make_hs_explicit else atom.GetTotalNumHs()
      chiral_tag = int(atom.GetChiralTag())
      chiral_neighbors = [
          str(b.GetOtherAtomIdx(atom.GetIdx())) for b in atom.GetBonds()
      ]

      # Calculate lone pairs (accounting for dative donor/acceptor bonds)
      v_free = VALENCE_ELECTRONS.get(elem, 4) if elem != "*" else 0
      explicit_val = 0
      for b in atom.GetBonds():
        if b.GetBondType() == Chem.BondType.DATIVE:
          if b.GetBeginAtomIdx() == atom.GetIdx():
            explicit_val += 2
        else:
          explicit_val += int(round(b.GetBondTypeAsDouble()))
      rem_e = v_free - q - rad - (explicit_val + impl_h) if elem != "*" else 0
      if strict_lewis and elem != "*":
        if rem_e < 0 or (rem_e % 2 != 0 and elem not in D_F_BLOCK_METALS):
          raise ValueError(
              f"Strict Lewis electron-pair parity violation on atom {idx_str}"
              f" ({elem}, q={q}, rad={rad}, bonds={explicit_val + impl_h}):"
              f" remaining non-bonding electrons rem_e={rem_e} (must be"
              " non-negative and even)."
          )
      lp = max(0, rem_e // 2)
      unpaired_d = (
          (rem_e % 2) if (elem in D_F_BLOCK_METALS and rem_e > 0) else 0
      )

      mg.add_atom(
          atom_id=idx_str,
          element=elem,
          lone_pairs=lp,
          radicals=rad,
          unpaired_d=unpaired_d,
          formal_charge=q,
          implicit_h=impl_h,
          chiral_tag=chiral_tag,
          chiral_neighbors=chiral_neighbors,
          rdkit_idx=atom.GetIdx(),
          is_aromatic=atom.GetIsAromatic(),
      )

    for bond in mol.GetBonds():
      u = str(bond.GetBeginAtomIdx())
      v = str(bond.GetEndAtomIdx())
      is_dat = bond.GetBondType() == Chem.BondType.DATIVE
      bo = 1 if is_dat else int(round(bond.GetBondTypeAsDouble()))
      if bo < 1:
        bo = 1
      stereo_str = "NONE"
      stereo_atoms_pair: Tuple[str, ...] = ()
      if bo == 2 and bond.GetStereo() in (
          Chem.BondStereo.STEREOE,
          Chem.BondStereo.STEREOZ,
          Chem.BondStereo.STEREOTRANS,
          Chem.BondStereo.STEREOCIS,
      ):
        sa = [str(idx) for idx in bond.GetStereoAtoms()]
        if len(sa) == 2:
          stereo_str = (
              "E"
              if bond.GetStereo()
              in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS)
              else "Z"
          )
          stereo_atoms_pair = (sa[0], sa[1])
      mg.add_bond(
          u,
          v,
          order=bo,
          stereo=stereo_str,
          stereo_atoms=stereo_atoms_pair,
          is_dative=is_dat,
          dative_begin=u if is_dat else None,
          dative_end=v if is_dat else None,
      )

    return mg

  def to_rdkit(self, sanitize: bool = True) -> Chem.Mol:
    """Converts this MolecularGraph into an RDKit Mol object and optionally validates

    with RDKit's C++ chemical sanitization engine.
    """
    rw = Chem.RWMol()
    node_to_idx: Dict[str, int] = {}

    for aid, an in self.atoms.items():
      rd_atom = Chem.Atom(an.element)
      rd_atom.SetFormalCharge(an.formal_charge)
      rd_atom.SetNumRadicalElectrons(an.radicals)
      if an.chiral_tag:
        rd_atom.SetChiralTag(Chem.ChiralType(an.chiral_tag))
      if an.element != "*":
        has_dummy_nb = any(
            self.atoms[nb].element == "*" for nb in self.graph.neighbors(aid)
        )
        if an.implicit_h > 0 or has_dummy_nb:
          rd_atom.SetNoImplicit(True)
          rd_atom.SetNumExplicitHs(int(max(0, an.implicit_h)))
      idx = rw.AddAtom(rd_atom)
      node_to_idx[aid] = idx

    added_edges = set()
    # First pass: for chiral atoms, add bonds in chiral_neighbors order to preserve CW/CCW frame
    for aid, an in self.atoms.items():
      if an.chiral_tag and an.chiral_neighbors:
        for nb in an.chiral_neighbors:
          if self.graph.has_edge(aid, nb):
            edge = _canon_pair(aid, nb)
            if edge not in added_edges:
              edata = self.graph[aid][nb]
              order = int(edata.get("order", 1))
              if (
                  edata.get("is_dative")
                  and edata.get("dative_begin")
                  and edata.get("dative_end")
              ):
                rw.AddBond(
                    node_to_idx[edata["dative_begin"]],
                    node_to_idx[edata["dative_end"]],
                    Chem.BondType.DATIVE,
                )
              else:
                bt = Chem.BondType.SINGLE
                if order == 2:
                  bt = Chem.BondType.DOUBLE
                elif order == 3:
                  bt = Chem.BondType.TRIPLE
                elif order == 4:
                  bt = Chem.BondType.QUADRUPLE
                rw.AddBond(node_to_idx[edge[0]], node_to_idx[edge[1]], bt)
              added_edges.add(edge)

    # Second pass: add remaining bonds
    for u, v, data in self.graph.edges(data=True):
      edge = _canon_pair(u, v)
      if edge not in added_edges:
        order = int(data.get("order", 1))
        if (
            data.get("is_dative")
            and data.get("dative_begin")
            and data.get("dative_end")
        ):
          rw.AddBond(
              node_to_idx[data["dative_begin"]],
              node_to_idx[data["dative_end"]],
              Chem.BondType.DATIVE,
          )
        else:
          idx_u = node_to_idx[u]
          idx_v = node_to_idx[v]
          if order == 1:
            rw.AddBond(idx_u, idx_v, Chem.BondType.SINGLE)
          elif order == 2:
            rw.AddBond(idx_u, idx_v, Chem.BondType.DOUBLE)
          elif order == 3:
            rw.AddBond(idx_u, idx_v, Chem.BondType.TRIPLE)
          elif order == 4:
            rw.AddBond(idx_u, idx_v, Chem.BondType.QUADRUPLE)
          else:
            rw.AddBond(idx_u, idx_v, Chem.BondType.SINGLE)
        added_edges.add(edge)

    # Restore Tier 4 E/Z double-bond stereochemistry
    for u, v, data in self.graph.edges(data=True):
      if int(data.get("order", 1)) == 2 and data.get("stereo") in ("E", "Z"):
        sa = data.get("stereo_atoms", ())
        if len(sa) == 2 and sa[0] in node_to_idx and sa[1] in node_to_idx:
          rd_bond = rw.GetBondBetweenAtoms(node_to_idx[u], node_to_idx[v])
          if rd_bond is not None:
            b_begin = rd_bond.GetBeginAtomIdx()
            b_end = rd_bond.GetEndAtomIdx()
            s0, s1 = node_to_idx[sa[0]], node_to_idx[sa[1]]
            if (
                rw.GetBondBetweenAtoms(b_begin, s0) is not None
                and rw.GetBondBetweenAtoms(b_end, s1) is not None
            ):
              rd_bond.SetStereoAtoms(s0, s1)
              rd_bond.SetStereo(
                  Chem.BondStereo.STEREOE
                  if data["stereo"] == "E"
                  else Chem.BondStereo.STEREOZ
              )
            elif (
                rw.GetBondBetweenAtoms(b_begin, s1) is not None
                and rw.GetBondBetweenAtoms(b_end, s0) is not None
            ):
              rd_bond.SetStereoAtoms(s1, s0)
              rd_bond.SetStereo(
                  Chem.BondStereo.STEREOE
                  if data["stereo"] == "E"
                  else Chem.BondStereo.STEREOZ
              )

    mol = rw.GetMol()
    try:
      Chem.SetDoubleBondNeighborDirections(mol)
    except Exception:
      pass
    if sanitize:
      Chem.SanitizeMol(mol)
    return mol

  def to_smiles(self, canonical: bool = True, isomeric: bool = True) -> str:
    """Exports the graph to a SMILES string via RDKit."""
    try:
      mol = self.to_rdkit(sanitize=False)
      return Chem.MolToSmiles(mol, canonical=canonical, isomericSmiles=isomeric)
    except Exception:
      mol = self.to_rdkit(sanitize=False)
      return Chem.MolToSmiles(mol, isomericSmiles=isomeric)


class ReactionVerifierEngine:
  """Formal verifier and graph-rewriting engine for electron pushing mechanisms."""

  @staticmethod
  def apply_and_verify_step(
      reactant_graph: MolecularGraph,
      arrows: List[CurvedArrow],
      step_name: str = "",
      sanitize_with_rdkit: bool = False,
      allow_multiproton: bool = False,
  ) -> Tuple[bool, MolecularGraph, Dict]:
    """Applies an elementary set of curved arrows to a molecular graph,

    enforces chemical locality and physical conservation laws,
    and returns (is_valid, product_graph, audit_log).
    """
    audit = {
        "step_name": step_name,
        "initial_state": reactant_graph.to_summary_dict(),
        "arrows_applied": len(arrows),
        "arrow_details": [str(a) for a in arrows],
        "invariants": {},
        "errors": [],
        "warnings": [],
        "status": "UNVERIFIED",
    }

    # ---------------------------------------------------------------------
    # 1. Pre-Flight & Topological Locality Checks
    # ---------------------------------------------------------------------
    for i, arrow in enumerate(arrows, 1):
      e_count = arrow.electrons

      # A. Source Verification
      if arrow.source_type == "lone_pair":
        aid = str(arrow.source_id)
        if aid not in reactant_graph.atoms:
          audit["errors"].append(
              f"Arrow {i}: Source atom '{aid}' does not exist in graph."
          )
        else:
          src_atom = reactant_graph.atoms[aid]
          avail_nonbond = 2 * src_atom.lone_pairs + (
              src_atom.unpaired_d
              if src_atom.element in D_F_BLOCK_METALS and e_count == 1
              else 0
          )
          if avail_nonbond < e_count:
            audit["errors"].append(
                f"Arrow {i}: Source atom '{aid}' has"
                f" {src_atom.lone_pairs} lone pairs; cannot push"
                f" {e_count}e⁻."
            )
      elif arrow.source_type == "radical":
        aid = str(arrow.source_id)
        if aid not in reactant_graph.atoms:
          audit["errors"].append(
              f"Arrow {i}: Source atom '{aid}' does not exist in graph."
          )
        else:
          src_atom = reactant_graph.atoms[aid]
          avail_rad = src_atom.radicals + (
              src_atom.unpaired_d if src_atom.element in D_F_BLOCK_METALS else 0
          )
          if avail_rad < 1:
            audit["errors"].append(
                f"Arrow {i}: Source atom '{aid}' has"
                f" {src_atom.radicals} radical electrons; cannot"
                f" push {e_count}e⁻."
            )
      elif arrow.source_type == "bond":
        u, v = str(arrow.source_id[0]), str(arrow.source_id[1])
        if not reactant_graph.graph.has_edge(u, v):
          audit["errors"].append(
              f"Arrow {i}: Source bond '{u}-{v}' does not exist in graph."
          )
        elif reactant_graph.get_bond_order(u, v) < 1:
          audit["errors"].append(
              f"Arrow {i}: Source bond '{u}-{v}' has order 0."
          )
      else:
        audit["errors"].append(
            f"Arrow {i}: Unknown source_type '{arrow.source_type}'."
        )

      # B. Sink Verification
      if arrow.sink_type == "atom":
        aid = str(arrow.sink_id)
        if aid not in reactant_graph.atoms:
          audit["errors"].append(
              f"Arrow {i}: Sink atom '{aid}' does not exist in graph."
          )
      elif arrow.sink_type == "bond":
        w, x = str(arrow.sink_id[0]), str(arrow.sink_id[1])
        if w not in reactant_graph.atoms or x not in reactant_graph.atoms:
          audit["errors"].append(
              f"Arrow {i}: Sink bond atoms '{w}' and '{x}' must exist in graph."
          )
      else:
        audit["errors"].append(
            f"Arrow {i}: Unknown sink_type '{arrow.sink_type}'."
        )

      # C. Chemical Locality Invariant (Physical non-teleportation)
      if arrow.source_type == "bond" and arrow.sink_type == "atom":
        u, v = str(arrow.source_id[0]), str(arrow.source_id[1])
        w = str(arrow.sink_id)
        if w not in (u, v):
          audit["errors"].append(
              f"Arrow {i}: Locality violation: Source bond '{u}-{v}' cannot"
              f" cleave onto unrelated atom '{w}'. Sink must be in source bond:"
              f" '{u}' or '{v}'."
          )

      elif (
          arrow.source_type in ("lone_pair", "radical")
          and arrow.sink_type == "bond"
      ):
        u = str(arrow.source_id)
        w, x = str(arrow.sink_id[0]), str(arrow.sink_id[1])
        if u not in (w, x):
          audit["errors"].append(
              f"Arrow {i}: Locality violation: Source atom '{u}' cannot form"
              f" disjoint bond '{w}-{x}' without connecting to source atom"
              f" '{u}'."
          )

      elif arrow.source_type == "bond" and arrow.sink_type == "bond":
        u, v = str(arrow.source_id[0]), str(arrow.source_id[1])
        w, x = str(arrow.sink_id[0]), str(arrow.sink_id[1])
        shared = set((u, v)) & set((w, x))
        if len(shared) == 0:
          audit["errors"].append(
              f"Arrow {i}: Locality violation: Source bond '{u}-{v}' cannot"
              f" shift to disjoint bond '{w}-{x}' without a shared pivot atom."
          )

      elif (
          arrow.source_type in ("lone_pair", "radical")
          and arrow.sink_type == "atom"
      ):
        u = str(arrow.source_id)
        w = str(arrow.sink_id)
        if u != w:
          audit["errors"].append(
              f"Arrow {i}: Locality violation: Electron on atom '{u}' cannot"
              f" teleport to unrelated atom '{w}' without forming a bond."
          )

    if audit["errors"]:
      audit["status"] = "REJECTED_LOCALITY_OR_SOURCE_VIOLATION"
      return False, reactant_graph, audit

    # ---------------------------------------------------------------------
    # 2. Graph Rewriting Transformation
    # ---------------------------------------------------------------------
    product_graph = reactant_graph.copy()
    product_graph.name = (
        f"Product of {step_name}"
        if step_name
        else f"Product({reactant_graph.name})"
    )

    # Track net bond order changes across all arrows
    net_bo_delta: Dict[Tuple[str, str], float] = Counter()

    for i, arrow in enumerate(arrows, 1):
      e_count = arrow.electrons

      # A. Source Depletion
      if arrow.source_type == "lone_pair":
        u = str(arrow.source_id)
        u_atom = product_graph.atoms[u]
        if e_count == 2:
          u_atom.lone_pairs -= 1
        elif e_count == 1:
          if u_atom.element in D_F_BLOCK_METALS:
            if u_atom.unpaired_d >= 1:
              u_atom.unpaired_d -= 1
            else:
              u_atom.lone_pairs -= 1
              u_atom.unpaired_d += 1
          else:
            u_atom.lone_pairs -= 1
            u_atom.radicals += 1
      elif arrow.source_type == "radical":
        u = str(arrow.source_id)
        u_atom = product_graph.atoms[u]
        if (
            u_atom.element in D_F_BLOCK_METALS
            and u_atom.radicals == 0
            and u_atom.unpaired_d >= 1
        ):
          u_atom.unpaired_d -= 1
        else:
          u_atom.radicals -= 1
      elif arrow.source_type == "bond":
        edge_pair = _canon_pair(arrow.source_id[0], arrow.source_id[1])
        if e_count == 2:
          net_bo_delta[edge_pair] -= 1.0
        elif e_count == 1:
          net_bo_delta[edge_pair] -= 0.5

      # B. Sink Acquisition
      if arrow.sink_type == "atom":
        w = str(arrow.sink_id)
        w_atom = product_graph.atoms[w]
        if e_count == 2:
          w_atom.lone_pairs += 1
        elif e_count == 1:
          if w_atom.element in D_F_BLOCK_METALS and w_atom.radicals == 0:
            if w_atom.unpaired_d >= 1:
              w_atom.unpaired_d -= 1
              w_atom.lone_pairs += 1
            else:
              w_atom.unpaired_d += 1
          else:
            w_atom.radicals += 1
      elif arrow.sink_type == "bond":
        edge_pair = _canon_pair(arrow.sink_id[0], arrow.sink_id[1])
        if e_count == 2:
          net_bo_delta[edge_pair] += 1.0
        elif e_count == 1:
          net_bo_delta[edge_pair] += 0.5

    # Apply net bond order changes
    for (u, v), delta in net_bo_delta.items():
      if not math.isclose(delta, round(delta), abs_tol=1e-5):
        audit["errors"].append(
            f"Unbalanced radical electron flow on bond '{u}-{v}': fractional"
            f" net order change {delta}."
        )
        audit["status"] = "FAILED_UNBALANCED_RADICALS"
        return False, product_graph, audit

      int_delta = int(round(delta))
      current_bo = product_graph.get_bond_order(u, v)
      new_bo = current_bo + int_delta
      if new_bo < 0:
        audit["errors"].append(
            f"Bond order on '{u}-{v}' became negative ({new_bo})."
        )
      elif new_bo == 0:
        product_graph.remove_bond(u, v)
      else:
        if product_graph.graph.has_edge(u, v):
          product_graph.graph[u][v]["order"] = new_bo
        else:
          is_dat = False
          dat_b: Optional[str] = None
          dat_e: Optional[str] = None
          if new_bo == 1:
            edge_pair = _canon_pair(u, v)
            edge_arrows = [
                arr
                for arr in arrows
                if arr.sink_type == "bond"
                and _canon_pair(arr.sink_id[0], arr.sink_id[1]) == edge_pair
            ]
            if (
                len(edge_arrows) == 1
                and edge_arrows[0].electrons == 2
                and edge_arrows[0].source_type == "lone_pair"
            ):
              arr0 = edge_arrows[0]
              donor = str(arr0.source_id)
              acceptor = v if u == donor else u
              if arr0.is_dative:
                is_dat = True
                dat_b, dat_e = donor, acceptor
              elif arr0.is_dative is None:
                donor_init = reactant_graph.atoms[donor]
                acceptor_init = reactant_graph.atoms[acceptor]
                other_donor_bo_change = any(
                    donor in ep and not math.isclose(dval, 0.0, abs_tol=1e-5)
                    for ep, dval in net_bo_delta.items()
                    if ep != edge_pair
                )
                if (
                    acceptor_init.element in D_F_BLOCK_METALS
                    and donor_init.element not in D_F_BLOCK_METALS
                    and donor_init.formal_charge == 0
                    and not other_donor_bo_change
                    and product_graph.atoms[donor].lone_pairs
                    == donor_init.lone_pairs - 1
                ):
                  is_dat = True
                  dat_b, dat_e = donor, acceptor
          product_graph.add_bond(
              u,
              v,
              order=new_bo,
              is_dative=is_dat,
              dative_begin=dat_b,
              dative_end=dat_e,
          )

    # ---------------------------------------------------------------------
    # 3. Deterministic Formal Charge Update via Differential Conservation
    # ---------------------------------------------------------------------
    # Delta Q(u) = - (2 * Delta L(u) + Delta R(u) + Delta U_d(u) + Delta BO(u))
    # This derives directly from Q = V_free - 2*L - R - U_d - BO - H_impl.
    for aid, prod_atom in product_graph.atoms.items():
      init_atom = reactant_graph.atoms[aid]
      delta_bo = product_graph.effective_bond_valence_on_atom(
          aid
      ) - reactant_graph.effective_bond_valence_on_atom(aid)
      delta_lp = prod_atom.lone_pairs - init_atom.lone_pairs
      delta_rad = prod_atom.radicals - init_atom.radicals
      delta_ud = prod_atom.unpaired_d - init_atom.unpaired_d
      delta_q = -int(round(2 * delta_lp + delta_rad + delta_ud + delta_bo))
      prod_atom.formal_charge = init_atom.formal_charge + delta_q

    # ---------------------------------------------------------------------
    # 3B. Stereochemical Invariant & State Update (Tiers 1–3 & Tier 3D)
    # ---------------------------------------------------------------------
    def _arrow_touches_atom(arr: CurvedArrow, atom_id: str) -> bool:
      return (
          arr.sink_id == atom_id
          or arr.source_id == atom_id
          or (isinstance(arr.sink_id, (tuple, list)) and atom_id in arr.sink_id)
          or (
              isinstance(arr.source_id, (tuple, list))
              and atom_id in arr.source_id
          )
      )

    tier3d_arrows = {
        idx: arr
        for idx, arr in enumerate(arrows)
        if arr.stereo_policy in ("R", "S", "Re", "Si", "re", "si")
    }
    satisfied_tier3d_arrows: set[int] = set()

    for aid, prod_atom in product_graph.atoms.items():
      init_atom = reactant_graph.atoms[aid]
      curr_neighbors = list(product_graph.graph.neighbors(aid))
      total_coordination = len(curr_neighbors) + prod_atom.implicit_h
      has_multiple_bond = any(
          product_graph.get_bond_order(aid, nb) > 1 for nb in curr_neighbors
      )
      old_nbs = set(reactant_graph.graph.neighbors(aid))
      new_nbs = set(curr_neighbors)
      broken_nbs = list(old_nbs - new_nbs)
      formed_nbs = list(new_nbs - old_nbs)

      if not init_atom.chiral_tag:
        # Tier 3D: Chirality Creation from Prochiral Center (Re/Si or R/S)
        if (
            not has_multiple_bond
            and total_coordination == 4
            and len(formed_nbs) >= 1
        ):
          incident_t3d = [
              (idx, arr)
              for idx, arr in tier3d_arrows.items()
              if _arrow_touches_atom(arr, aid)
          ]
          if incident_t3d:
            policy = incident_t3d[0][1].stereo_policy
            cw_tag = int(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            ccw_tag = int(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
            pol_up = policy.upper()
            ordered_old = [nb for nb in curr_neighbors if nb in old_nbs]
            prod_atom.chiral_neighbors = ordered_old + formed_nbs

            # Map Re/Si on sp2 -> sp3 addition to CW/CCW via Z=0 phantom probe
            target_tag = 0
            if (
                pol_up in ("RE", "SI")
                and len(ordered_old) + init_atom.implicit_h == 3
                and len(formed_nbs) == 1
            ):
              r_probe = copy.deepcopy(reactant_graph)
              r_probe.atoms[aid].chiral_neighbors = ordered_old
              for cand_tag in (cw_tag, ccw_tag):
                r_probe.atoms[aid].chiral_tag = cand_tag
                rw_p = Chem.RWMol(r_probe.to_rdkit(sanitize=False))
                c_i = list(r_probe.atoms.keys()).index(aid)
                d_i = rw_p.AddAtom(Chem.Atom(0))
                rw_p.AddBond(c_i, d_i, Chem.BondType.SINGLE)
                rw_p.GetAtomWithIdx(c_i).SetChiralTag(Chem.ChiralType(cand_tag))
                rw_p.UpdatePropertyCache(strict=False)
                rw_p_h = Chem.RWMol(Chem.AddHs(rw_p))
                for s in Chem.FindPotentialStereo(rw_p_h):
                  if (
                      s.type == Chem.StereoType.Atom_Tetrahedral
                      and s.centeredOn != c_i
                  ):
                    a_s = rw_p_h.GetAtomWithIdx(s.centeredOn)
                    if a_s.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                      a_s.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                Chem.AssignStereochemistry(rw_p_h, cleanIt=True, force=True)
                rdCIPLabeler.AssignCIPLabels(rw_p_h)
                cip = (
                    rw_p_h.GetAtomWithIdx(c_i).GetPropsAsDict().get("_CIPCode")
                )
                if (pol_up == "RE" and cip in ("S", "s")) or (
                    pol_up == "SI" and cip in ("R", "r")
                ):
                  target_tag = cand_tag
                  break
            elif pol_up in ("R", "S"):
              c_i = list(product_graph.atoms.keys()).index(aid)
              for cand_tag in (cw_tag, ccw_tag):
                prod_atom.chiral_tag = cand_tag
                m_p = product_graph.to_rdkit(sanitize=False)
                m_p.UpdatePropertyCache(strict=False)
                m_p_h = Chem.RWMol(Chem.AddHs(m_p))
                Chem.AssignStereochemistry(m_p_h, cleanIt=True, force=True)
                rdCIPLabeler.AssignCIPLabels(m_p_h)
                cip = m_p_h.GetAtomWithIdx(c_i).GetPropsAsDict().get("_CIPCode")
                if cip and cip.upper() == pol_up:
                  target_tag = cand_tag
                  break

            # Verify resulting product center is genuinely stereogenic
            if target_tag:
              prod_atom.chiral_tag = target_tag
              m_chk = product_graph.to_rdkit(sanitize=False)
              m_chk.UpdatePropertyCache(strict=False)
              c_i = list(product_graph.atoms.keys()).index(aid)
              pot_chk = {
                  s.centeredOn
                  for s in Chem.FindPotentialStereo(m_chk)
                  if s.type == Chem.StereoType.Atom_Tetrahedral
              }
              Chem.AssignStereochemistry(m_chk, cleanIt=True, force=True)
              if (
                  m_chk.GetAtomWithIdx(c_i).GetChiralTag()
                  == Chem.ChiralType.CHI_UNSPECIFIED
                  and c_i not in pot_chk
              ):
                target_tag = 0
                prod_atom.chiral_tag = 0
                prod_atom.chiral_neighbors = []

            if target_tag:
              for idx, _ in incident_t3d:
                satisfied_tier3d_arrows.add(idx)
            else:
              prod_atom.chiral_tag = 0
              prod_atom.chiral_neighbors = []
        continue

      # Tier 2: Planarization / Chiral Destruction
      # If atom forms multiple bonds (order >= 2) or coordination drops below 4
      if has_multiple_bond or total_coordination < 4:
        prod_atom.chiral_tag = 0
        prod_atom.chiral_neighbors = []
        continue

      # Tier 1: Spectator center (no incident bonds broken or formed)
      if not broken_nbs and not formed_nbs:
        prod_atom.chiral_tag = init_atom.chiral_tag
        prod_atom.chiral_neighbors = list(init_atom.chiral_neighbors)
        continue

      # Tier 3: Substitution at tetrahedral center (1 bond broken, 1 bond formed)
      if len(broken_nbs) == 1 and len(formed_nbs) == 1:
        lg = broken_nbs[0]
        nuc = formed_nbs[0]

        # Determine stereo_policy from arrows incident to aid
        policy = "auto"
        incident_t3_idxs = []
        for idx, arr in enumerate(arrows):
          if (
              arr.stereo_policy
              and arr.stereo_policy != "auto"
              and _arrow_touches_atom(arr, aid)
          ):
            if policy == "auto":
              policy = arr.stereo_policy
            if idx in tier3d_arrows:
              incident_t3_idxs.append(idx)

        new_chiral_nbs = list(init_atom.chiral_neighbors)
        if lg in new_chiral_nbs:
          k = new_chiral_nbs.index(lg)
          new_chiral_nbs[k] = nuc
        else:
          new_chiral_nbs.append(nuc)

        if policy in ("auto", "invert"):
          # Slot replacement: substitute nuc into lg's slot in chiral_neighbors
          prod_atom.chiral_neighbors = new_chiral_nbs

          # Physical Walden inversion: flip CW <-> CCW
          cw_tag = int(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
          ccw_tag = int(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
          if init_atom.chiral_tag == cw_tag:
            prod_atom.chiral_tag = ccw_tag
          elif init_atom.chiral_tag == ccw_tag:
            prod_atom.chiral_tag = cw_tag
          else:
            prod_atom.chiral_tag = 0

        elif policy == "retain":
          # Slot replacement with retention of physical configuration
          prod_atom.chiral_neighbors = new_chiral_nbs
          prod_atom.chiral_tag = init_atom.chiral_tag

        elif policy == "racemize":
          prod_atom.chiral_tag = 0
          prod_atom.chiral_neighbors = []

        elif policy.upper() in ("R", "S"):
          prod_atom.chiral_neighbors = new_chiral_nbs
          cw_tag = int(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
          ccw_tag = int(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
          pol_up = policy.upper()
          c_i = list(product_graph.atoms.keys()).index(aid)
          target_tag = 0
          for cand_tag in (cw_tag, ccw_tag):
            prod_atom.chiral_tag = cand_tag
            m_p = product_graph.to_rdkit(sanitize=False)
            m_p.UpdatePropertyCache(strict=False)
            m_p_h = Chem.RWMol(Chem.AddHs(m_p))
            Chem.AssignStereochemistry(m_p_h, cleanIt=True, force=True)
            rdCIPLabeler.AssignCIPLabels(m_p_h)
            cip = m_p_h.GetAtomWithIdx(c_i).GetPropsAsDict().get("_CIPCode")
            if cip and cip.upper() == pol_up:
              target_tag = cand_tag
              break
          if target_tag:
            prod_atom.chiral_tag = target_tag
            for idx in incident_t3_idxs:
              satisfied_tier3d_arrows.add(idx)
          else:
            prod_atom.chiral_tag = 0
            prod_atom.chiral_neighbors = []
            for idx in incident_t3_idxs:
              satisfied_tier3d_arrows.add(idx)
            audit["errors"].append(
                f"Tier 3 Stereochemistry Violation: Atom '{aid}' cannot be"
                f" assigned stereo_policy='{policy}' because the substituted"
                " center is achiral in the product."
            )
        else:
          prod_atom.chiral_tag = 0
          prod_atom.chiral_neighbors = []
      else:
        # Elimination, addition, or multiple bond changes at the center
        # (len(broken_nbs) != 1 or len(formed_nbs) != 1): clear tetrahedral
        # stereochemistry.
        prod_atom.chiral_tag = 0
        prod_atom.chiral_neighbors = []

    for idx, arr in tier3d_arrows.items():
      if idx not in satisfied_tier3d_arrows:
        audit["errors"].append(
            f"Tier 3D Prochirality/Stereogenesis Violation: Arrow {arr}"
            f" specifies stereo_policy='{arr.stereo_policy}', but no incident"
            " reaction center creates a valid stereocenter (ligands/faces are"
            " homotopic)."
        )

    # ---------------------------------------------------------------------
    # 3C. Tier 4 Alkene E/Z Stereochemistry & Ring-Strain Invariants
    # ---------------------------------------------------------------------
    def _best_alkene_ref_neighbor(
        g_obj: MolecularGraph, center: str, other: str
    ) -> Optional[str]:
      candidates = [nb for nb in g_obj.graph.neighbors(center) if nb != other]
      if not candidates:
        return None
      # Prefer heavy atoms over explicit hydrogens, then by group valence / id
      candidates.sort(
          key=lambda nb: (
              0 if g_obj.atoms[nb].element == "H" else 1,
              VALENCE_ELECTRONS.get(g_obj.atoms[nb].element, 4),
              nb,
          ),
          reverse=True,
      )
      return candidates[0]

    for u, v, p_data in list(product_graph.graph.edges(data=True)):
      p_order = int(p_data.get("order", 1))
      r_order = reactant_graph.get_bond_order(u, v)

      # Tier 4B: Pi-bond cleavage or addition destroys E/Z stereochemistry
      if p_order != 2:
        p_data["stereo"] = "NONE"
        p_data["stereo_atoms"] = ()
        continue

      # Tier 4A: Spectator double bond (order 2 -> 2) preserves E/Z stereochemistry
      if r_order == 2 and reactant_graph.graph.has_edge(u, v):
        r_data = reactant_graph.graph[u][v]
        r_stereo = r_data.get("stereo", "NONE")
        r_sa = r_data.get("stereo_atoms", ())
        if r_stereo in ("E", "Z") and len(r_sa) == 2:
          ref_u, ref_v = r_sa
          if (
              product_graph.graph.has_edge(u, ref_u)
              and product_graph.graph.has_edge(v, ref_v)
          ) or (
              product_graph.graph.has_edge(u, ref_v)
              and product_graph.graph.has_edge(v, ref_u)
          ):
            p_data["stereo"] = r_stereo
            p_data["stereo_atoms"] = (ref_u, ref_v)
          else:
            p_data["stereo"] = "NONE"
            p_data["stereo_atoms"] = ()
        continue

      # Tier 4C: Newly formed double bond (order < 2 -> 2)
      if r_order < 2 and p_order == 2:
        # Check for arrow-level stereo_policy targeting this bond or global step
        ez_policy = None
        for arr in arrows:
          if arr.stereo_policy in ("E", "Z", "trans", "cis", "anti", "syn"):
            if (
                arr.sink_type == "bond"
                and isinstance(arr.sink_id, (tuple, list))
                and _canon_pair(arr.sink_id[0], arr.sink_id[1])
                == _canon_pair(u, v)
            ):
              ez_policy = arr.stereo_policy
              break
            if ez_policy is None:
              ez_policy = arr.stereo_policy

        # Check shortest cycle containing (u, v)
        temp_g = product_graph.graph.copy()
        temp_g.remove_edge(u, v)
        cycle_path = None
        if nx.has_path(temp_g, u, v):
          cycle_path = nx.shortest_path(temp_g, u, v)

        if cycle_path is not None and 3 <= len(cycle_path) <= 7:
          ring_size = len(cycle_path)
          bredt_err = ReactionVerifierEngine._check_bredt_bridgehead_strain(
              product_graph.graph, u, v
          )
          if bredt_err is not None:
            audit["errors"].append(bredt_err)
          elif ez_policy in ("E", "trans"):
            audit["errors"].append(
                "Tier 4 Stereochemistry Violation: Newly formed trans (E)"
                f" double bond ({u}={v}) in a {ring_size}-membered ring"
                " violates cycloalkene ring-strain constraints (rings <= 7"
                " must be Z/cis)."
            )
          else:
            p_data["stereo"] = "Z"
            p_data["stereo_atoms"] = (cycle_path[1], cycle_path[-2])
        elif ez_policy in ("E", "Z", "trans", "cis", "anti", "syn"):
          mapped_ez = "E" if ez_policy in ("E", "trans", "anti") else "Z"
          sub_u = _best_alkene_ref_neighbor(product_graph, u, v)
          sub_v = _best_alkene_ref_neighbor(product_graph, v, u)
          if sub_u is not None and sub_v is not None:
            p_data["stereo"] = mapped_ez
            p_data["stereo_atoms"] = (sub_u, sub_v)

    # ---------------------------------------------------------------------
    # 4. Rigorous Chemical Invariant Audits
    # ---------------------------------------------------------------------
    # Invariant I: Total Valence Electron Conservation
    init_e = reactant_graph.total_electrons()
    final_e = product_graph.total_electrons()
    e_conserved = init_e == final_e
    audit["invariants"]["electron_pair_conservation"] = {
        "passed": e_conserved,
        "initial_valence_electrons": init_e,
        "final_valence_electrons": final_e,
    }
    if not e_conserved:
      audit["errors"].append(
          f"Violation: Valence electrons changed from {init_e} to {final_e}."
      )

    # Invariant II: Global Net Charge Conservation
    init_q = reactant_graph.net_charge()
    final_q = product_graph.net_charge()
    q_conserved = init_q == final_q
    audit["invariants"]["charge_conservation"] = {
        "passed": q_conserved,
        "initial_net_charge": init_q,
        "final_net_charge": final_q,
    }
    if not q_conserved:
      audit["errors"].append(
          f"Violation: Net formal charge changed from {init_q} to {final_q}."
      )

    # Invariant III: Nuclei / Mass Conservation
    init_nuclei = reactant_graph.nuclei_multiset()
    final_nuclei = product_graph.nuclei_multiset()
    nuclei_conserved = init_nuclei == final_nuclei
    audit["invariants"]["nuclei_conservation"] = {
        "passed": nuclei_conserved,
        "initial_nuclei": dict(init_nuclei),
        "final_nuclei": dict(final_nuclei),
    }
    if not nuclei_conserved:
      audit["errors"].append(
          f"Violation: Nuclei/mass changed. Initial: {dict(init_nuclei)},"
          f" Final: {dict(final_nuclei)}."
      )

    # Invariant IV: Lewis Octet Rule / Capacity Check (Non-Hypervalency)
    octet_violations = []
    for aid, atom in product_graph.atoms.items():
      total_bo = product_graph.total_bond_order_on_atom(aid)
      tot_val = atom.total_valence_electrons(total_bo)
      cap = MAX_VALENCE_CAPACITY.get(atom.element, 8)
      if tot_val > cap:
        octet_violations.append(
            f"Atom '{aid}' ({atom.element}) has {tot_val} valence electrons"
            f" (exceeds max capacity {cap})"
        )

    audit["invariants"]["octet_rule_invariance"] = {
        "passed": len(octet_violations) == 0,
        "violations": octet_violations,
    }
    if octet_violations:
      audit["errors"].extend(octet_violations)

    # Invariant V: Non-negative Lone Pairs and Radicals
    negative_lp = [
        aid for aid, an in product_graph.atoms.items() if an.lone_pairs < 0
    ]
    if negative_lp:
      audit["errors"].append(
          f"Violation: Negative lone pairs on atoms: {negative_lp}"
      )

    negative_rad = [
        aid
        for aid, an in product_graph.atoms.items()
        if an.radicals < 0 or an.unpaired_d < 0
    ]
    if negative_rad:
      audit["errors"].append(
          f"Violation: Negative radicals on atoms: {negative_rad}"
      )

    # Invariant VI: Tier 5 Elementary Step (Single-Proton Transfer Rule)
    transferred_hydrogens = []
    for aid, an in reactant_graph.atoms.items():
      if an.element == "H" and set(reactant_graph.graph.neighbors(aid)) != set(
          product_graph.graph.neighbors(aid)
      ):
        is_hydride_shift = any(
            arr.source_type == "bond"
            and arr.sink_type == "bond"
            and isinstance(arr.source_id, (tuple, list))
            and isinstance(arr.sink_id, (tuple, list))
            and (aid in arr.source_id)
            and (aid in arr.sink_id)
            for arr in arrows
        )
        if not is_hydride_shift:
          transferred_hydrogens.append(aid)
    single_proton_ok = len(transferred_hydrogens) <= 1 or bool(
        allow_multiproton
    )
    audit["invariants"]["elementary_step_single_proton_transfer"] = {
        "passed": single_proton_ok,
        "transferred_hydrogens": transferred_hydrogens,
    }
    if not single_proton_ok:
      audit["errors"].append(
          "Tier 5 Elementary Step Violation: Step transfers"
          f" {len(transferred_hydrogens)} distinct protons/hydrogens"
          f" ({transferred_hydrogens}) simultaneously. A single elementary step"
          " may transfer at most 1 proton (split multi-proton composite steps"
          " such as enolization + tautomerization into discrete sequential"
          " elementary steps, or pass allow_multiproton=True for a synchronous"
          " concerted shuttle)."
      )

    # Invariant VI: Optional RDKit Chemical Sanitization Validator
    if sanitize_with_rdkit:
      try:
        product_graph.to_rdkit(sanitize=True)
        audit["invariants"]["rdkit_sanitization"] = {"passed": True}
      except Exception as e:
        audit["invariants"]["rdkit_sanitization"] = {
            "passed": False,
            "error": str(e),
        }
        audit["warnings"].append(f"RDKit sanitization warning: {e}")

    # Invariant VII: Spectator Locality Invariance
    passed_locality, locality_violations = (
        ReactionVerifierEngine.verify_spectator_locality(
            reactant_graph, product_graph, arrows
        )
    )
    audit["invariants"]["spectator_locality"] = {
        "passed": passed_locality,
        "violations": locality_violations,
    }
    if not passed_locality:
      audit["errors"].extend(locality_violations)

    audit["final_state"] = product_graph.to_summary_dict()

    if not audit["errors"]:
      audit["status"] = "PASSED_RIGOROUS_VERIFICATION"
      return True, product_graph, audit
    else:
      audit["status"] = "FAILED_VERIFICATION"
      return False, product_graph, audit

  # -------------------------------------------------------------------------
  # Spectator Locality Verification
  # -------------------------------------------------------------------------
  @staticmethod
  def verify_spectator_locality(
      reactant_graph: MolecularGraph,
      product_graph: MolecularGraph,
      arrows: Sequence[CurvedArrow],
  ) -> Tuple[bool, List[str]]:
    """Invariant VII: Spectator Locality Invariance.

    Asserts that all atoms not explicitly involved as arrow sources or sinks
    strictly maintain their incident edge list (neighbors and bond orders),
    formal charges, lone pairs, and radical electrons.
    """
    active_atoms = set()
    for arr in arrows:
      for target in (arr.source_id, arr.sink_id):
        if isinstance(target, (str, int)):
          active_atoms.add(str(target))
        elif isinstance(target, (tuple, list)):
          active_atoms.update(str(x) for x in target)

    violations = []

    for aid, r_atom in reactant_graph.atoms.items():
      if aid in active_atoms:
        continue

      p_atom = product_graph.atoms.get(aid)
      if p_atom is None:
        violations.append(
            f"Spectator atom '{aid}' was removed from the product graph."
        )
        continue

      # A. Invariant local electronic state
      if p_atom.formal_charge != r_atom.formal_charge:
        violations.append(
            f"Spectator atom '{aid}' formal charge mutated from"
            f" {r_atom.formal_charge} to {p_atom.formal_charge}."
        )
      if p_atom.lone_pairs != r_atom.lone_pairs:
        violations.append(
            f"Spectator atom '{aid}' lone pairs mutated from"
            f" {r_atom.lone_pairs} to {p_atom.lone_pairs}."
        )
      if p_atom.radicals != r_atom.radicals:
        violations.append(
            f"Spectator atom '{aid}' radicals mutated from"
            f" {r_atom.radicals} to {p_atom.radicals}."
        )
      if p_atom.unpaired_d != r_atom.unpaired_d:
        violations.append(
            f"Spectator atom '{aid}' unpaired_d mutated from"
            f" {r_atom.unpaired_d} to {p_atom.unpaired_d}."
        )

      # B. Invariant incident edge list: {neighbor: bond_order}
      r_edges = {
          nb: reactant_graph.get_bond_order(aid, nb)
          for nb in reactant_graph.graph.neighbors(aid)
      }
      p_edges = {
          nb: product_graph.get_bond_order(aid, nb)
          for nb in product_graph.graph.neighbors(aid)
      }
      if r_edges != p_edges:
        violations.append(
            f"Spectator atom '{aid}' incident edge list mutated.\n"
            f"  Reactant incident edges: {r_edges}\n"
            f"  Product incident edges:  {p_edges}"
        )

    return len(violations) == 0, violations

  # -------------------------------------------------------------------------
  # Target Product SMILES / Graph Isomorphism Verification
  # -------------------------------------------------------------------------
  @staticmethod
  def verify_product_match(
      product_graph: MolecularGraph,
      target_product: Union[str, Chem.Mol, MolecularGraph],
  ) -> Tuple[bool, str]:
    """Validates whether product_graph matches an expected target product

    via canonical SMILES equality (with implicit/explicit H normalization).
    """
    try:
      prod_rdkit = product_graph.to_rdkit(sanitize=False)
      # Remove explicit hydrogens that were only added for mechanism arrow tracking
      try:
        prod_mol_clean = Chem.RemoveHs(prod_rdkit)
      except Exception:
        prod_mol_clean = None
      if prod_mol_clean is None:
        prod_mol_clean = prod_rdkit

      prod_smi = Chem.MolToSmiles(
          prod_mol_clean, canonical=True, isomericSmiles=True
      )

      if isinstance(target_product, str):
        raw_t = Chem.MolFromSmiles(target_product, sanitize=False)
        params_t = Chem.SmilesParserParams()
        params_t.removeHs = False
        target_mol = Chem.MolFromSmiles(target_product, params_t)
        if target_mol is None or raw_t is None:
          return False, f"Invalid target SMILES: '{target_product}'"
        Chem.AssignStereochemistry(target_mol, cleanIt=True, force=True)
        pot_t = {
            s.centeredOn
            for s in Chem.FindPotentialStereo(target_mol)
            if s.type == Chem.StereoType.Atom_Tetrahedral
        }
        false_chiral = [
            a.GetIdx()
            for a in raw_t.GetAtoms()
            if a.GetChiralTag()
            in (
                Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
            )
            and target_mol.GetAtomWithIdx(a.GetIdx()).GetChiralTag()
            == Chem.ChiralType.CHI_UNSPECIFIED
            and a.GetIdx() not in pot_t
        ]
        if false_chiral:
          return (
              False,
              (
                  f"Target SMILES '{target_product}' specifies @/@@ chirality"
                  f" on homotopic/achiral atom(s) {false_chiral}."
              ),
          )
        try:
          target_mol_clean = Chem.RemoveHs(target_mol)
        except Exception:
          target_mol_clean = None
        if target_mol_clean is None:
          target_mol_clean = target_mol
        target_smi = Chem.MolToSmiles(
            target_mol_clean, canonical=True, isomericSmiles=True
        )
      elif isinstance(target_product, Chem.Mol):
        try:
          target_mol_clean = Chem.RemoveHs(target_product)
        except Exception:
          target_mol_clean = None
        if target_mol_clean is None:
          target_mol_clean = target_product
        target_smi = Chem.MolToSmiles(
            target_mol_clean, canonical=True, isomericSmiles=True
        )
      elif isinstance(target_product, MolecularGraph):
        target_smi = target_product.to_smiles(canonical=True, isomeric=True)
      else:
        return False, f"Unsupported target_product type: {type(target_product)}"

      # Canonical SMILES comparison
      if Chem.CanonSmiles(prod_smi) == Chem.CanonSmiles(target_smi):
        return (
            True,
            (
                f"Product canonical SMILES '{prod_smi}' matches target"
                f" '{target_smi}'."
            ),
        )

      return (
          False,
          (
              f"Product SMILES '{prod_smi}' does not match expected target"
              f" '{target_smi}'."
          ),
      )

    except Exception as e:
      return False, f"Error comparing products: {e}"

  @staticmethod
  def verify_step_against_target(
      reactant_graph: MolecularGraph,
      arrows: Sequence[CurvedArrow],
      target_product: Union[str, Chem.Mol, MolecularGraph],
      step_name: str = "",
  ) -> Tuple[bool, MolecularGraph, Dict[str, Any]]:
    """Applies arrows, verifies all invariants (including spectator locality),

    and strictly asserts that the resulting product matches the target product.
    """
    valid, prod_graph, audit = ReactionVerifierEngine.apply_and_verify_step(
        reactant_graph, arrows, step_name=step_name
    )
    if not valid:
      return False, prod_graph, audit

    match, msg = ReactionVerifierEngine.verify_product_match(
        prod_graph, target_product
    )
    audit["invariants"]["target_product_match"] = {
        "passed": match,
        "message": msg,
    }
    if not match:
      audit["errors"].append(f"Target product mismatch: {msg}")
      audit["status"] = "FAILED_TARGET_MISMATCH"
      return False, prod_graph, audit

    return True, prod_graph, audit

  @staticmethod
  def _check_bredt_bridgehead_strain(
      graph: nx.Graph, u: str, v: str
  ) -> Optional[str]:
    """Checks whether newly formed double bond (u=v) violates Bredt's rule.

    A double bond at a bridgehead vertex b1 of a bridged bicyclic system
    bicyclo[x.y.z] (with three internally vertex-disjoint bridges of interior
    atom counts x >= y >= z >= 1 between bridgeheads b1 and b2, one of which
    contains edge (u, v)) violates Bredt's rule when Fawcett's S = x + y + z <
    7.
    Fused bicyclic systems (z == 0) are not flagged.
    """
    for b1 in (u, v):
      if graph.degree(b1) < 3:
        continue
      other_end = v if b1 == u else u
      for b2 in sorted(graph.nodes()):
        if b2 == b1 or graph.degree(b2) < 3:
          continue
        # Collect simple paths between b1 and b2 of edge length 2..6 (1..5 interior atoms)
        paths_with_uv: List[Tuple[int, FrozenSet[str]]] = []
        paths_other: List[Tuple[int, FrozenSet[str]]] = []
        for path in nx.all_simple_paths(graph, b1, b2, cutoff=6):
          if len(path) < 3:
            continue
          interior = frozenset(path[1:-1])
          if path[1] == other_end:
            paths_with_uv.append((len(interior), interior))
          else:
            paths_other.append((len(interior), interior))
        if not paths_with_uv or len(paths_other) < 2:
          continue
        paths_with_uv.sort(key=lambda t: t[0])
        paths_other.sort(key=lambda t: t[0])
        for len1, int1 in paths_with_uv:
          for i, (len2, int2) in enumerate(paths_other):
            if int1 & int2:
              continue
            for len3, int3 in paths_other[i + 1 :]:
              if (int1 & int3) or (int2 & int3):
                continue
              x, y, z = sorted((len1, len2, len3), reverse=True)
              s_val = x + y + z
              if z >= 1 and s_val < 7:
                return (
                    "Tier 4 Stereochemistry Violation (Bredt's Rule): Newly"
                    f" formed bridgehead double bond ({u}={v}) at bridgehead"
                    f" '{b1}' in bridged bicyclic system"
                    f" bicyclo[{x}.{y}.{z}] has Fawcett S = {s_val} < 7."
                )
    return None

  @staticmethod
  def solve_arrows_from_mapped_graphs(
      reactant_graph: MolecularGraph,
      product_graph: MolecularGraph,
      atom_map: Optional[Dict[str, str]] = None,
      allow_multiproton: bool = False,
  ) -> Optional[List[CurvedArrow]]:
    """Synthesizes a verified CurvedArrow list connecting atom-mapped reactant and product graphs.

    Enforces the Dugundji-Ugi vertex-flow conservation law at every atom v:
      Delta B_eff(v) + Delta q(v) + Delta R(v) + Delta U_d(v) + 2 * Delta L(v)
      == 0
    Supports both 2-electron polar/pericyclic/vinylogous steps (arbitrary arrow
    count m >= 1) and 1-electron radical/open-shell metal fishhook steps
    (electrons=1).
    """
    if atom_map is not None:
      p_aligned = MolecularGraph(name=f"{product_graph.name}_aligned")
      inv_map = {p_id: r_id for r_id, p_id in atom_map.items()}
      for r_id, p_id in atom_map.items():
        pa = product_graph.atoms[p_id]
        p_aligned.add_atom(
            r_id,
            pa.element,
            lone_pairs=pa.lone_pairs,
            radicals=pa.radicals,
            unpaired_d=pa.unpaired_d,
            formal_charge=pa.formal_charge,
            implicit_h=pa.implicit_h,
        )
      for pu, pv, d in product_graph.graph.edges(data=True):
        if pu in inv_map and pv in inv_map:
          is_dat = bool(d.get("is_dative", False))
          dat_b = (
              inv_map[d["dative_begin"]]
              if is_dat and d.get("dative_begin") in inv_map
              else None
          )
          dat_e = (
              inv_map[d["dative_end"]]
              if is_dat and d.get("dative_end") in inv_map
              else None
          )
          p_aligned.add_bond(
              inv_map[pu],
              inv_map[pv],
              order=int(d.get("order", 1)),
              is_dative=is_dat,
              dative_begin=dat_b,
              dative_end=dat_e,
          )
      p_graph = p_aligned
    else:
      p_graph = product_graph

    atoms = list(reactant_graph.atoms.keys())
    if set(atoms) != set(p_graph.atoms.keys()):
      return None

    def _edge_orders(g: MolecularGraph) -> Dict[Tuple[str, str], int]:
      return {
          _canon_pair(u, v): int(d.get("order", 1))
          for u, v, d in g.graph.edges(data=True)
      }

    r_edges = _edge_orders(reactant_graph)
    p_edges = _edge_orders(p_graph)
    dec_edges: List[Tuple[str, str]] = []
    inc_edges: List[Tuple[str, str]] = []

    for u, v in sorted(set(r_edges.keys()) | set(p_edges.keys())):
      b_r = r_edges.get((u, v), 0)
      b_p = p_edges.get((u, v), 0)
      if b_p < b_r:
        for _ in range(b_r - b_p):
          dec_edges.append((u, v))
      elif b_p > b_r:
        for _ in range(b_p - b_r):
          inc_edges.append((u, v))

    if not dec_edges and not inc_edges:
      return None

    delta_l: Dict[str, int] = {}
    delta_r: Dict[str, int] = {}
    delta_u: Dict[str, int] = {}
    delta_nonbond: Dict[str, int] = {}
    for a in atoms:
      ra = reactant_graph.atoms[a]
      pa = p_graph.atoms[a]
      if ra.element != pa.element or ra.implicit_h != pa.implicit_h:
        return None
      db = p_graph.effective_bond_valence_on_atom(
          a
      ) - reactant_graph.effective_bond_valence_on_atom(a)
      dq = int(pa.formal_charge) - int(ra.formal_charge)
      dr = int(pa.radicals) - int(ra.radicals)
      du = int(pa.unpaired_d) - int(ra.unpaired_d)
      dl = int(pa.lone_pairs) - int(ra.lone_pairs)
      if db + dq + dr + du + 2 * dl != 0:
        return None
      if int(ra.lone_pairs) + dl < 0:
        return None
      delta_l[a] = dl
      delta_r[a] = dr
      delta_u[a] = du
      delta_nonbond[a] = 2 * dl + dr + du

    # 1. 1-Electron Radical Fishhook Mode (electrons == 1)
    active_edge_atoms = {x for e in dec_edges + inc_edges for x in e}
    has_1e_change = any(
        delta_nonbond[a] % 2 != 0 or (delta_r[a] + delta_u[a]) != 0
        for a in atoms
    )
    has_active_radical = any(
        reactant_graph.atoms[a].radicals > 0 or p_graph.atoms[a].radicals > 0
        for a in active_edge_atoms
    )
    can_try_1e = all(
        delta_l[a] == 0
        or (
            has_1e_change
            and reactant_graph.atoms[a].element in D_F_BLOCK_METALS
            and abs(delta_nonbond[a]) <= 2
        )
        or (delta_l[a] < 0 and (delta_r[a] + delta_u[a]) > 0)
        for a in atoms
    )
    if (has_1e_change or has_active_radical) and can_try_1e:
      in_1e: Dict[str, List[Tuple[str, Any]]] = {a: [] for a in atoms}
      out_1e: Dict[str, List[Tuple[str, Any]]] = {a: [] for a in atoms}
      for a, dn in delta_nonbond.items():
        if dn < 0:
          ra = reactant_graph.atoms[a]
          avail_rad = int(ra.radicals) + (
              int(ra.unpaired_d) if ra.element in D_F_BLOCK_METALS else 0
          )
          n_rad_src = min(-dn, avail_rad)
          n_lp_src = (-dn) - n_rad_src
          for _ in range(n_rad_src):
            in_1e[a].append(("radical", a))
          for _ in range(n_lp_src):
            in_1e[a].append(("lone_pair", a))
        elif dn > 0:
          for _ in range(dn):
            out_1e[a].append(("atom", a))
      for e in dec_edges:
        in_1e[e[0]].append(("bond", e))
        in_1e[e[1]].append(("bond", e))
      for e in inc_edges:
        out_1e[e[0]].append(("bond", e))
        out_1e[e[1]].append(("bond", e))

      if all(len(in_1e[a]) == len(out_1e[a]) for a in atoms):
        active_1e = [a for a in atoms if in_1e[a]]

        def _pair_1e(
            v_idx: int, acc: List[CurvedArrow]
        ) -> Optional[List[CurvedArrow]]:
          if v_idx == len(active_1e):
            ok, prod, _ = ReactionVerifierEngine.apply_and_verify_step(
                reactant_graph, acc, allow_multiproton=allow_multiproton
            )
            if ok and _edge_orders(prod) == p_edges:
              if all(
                  prod.atoms[a].formal_charge == p_graph.atoms[a].formal_charge
                  and (prod.atoms[a].radicals + prod.atoms[a].unpaired_d)
                  == (p_graph.atoms[a].radicals + p_graph.atoms[a].unpaired_d)
                  and prod.atoms[a].lone_pairs == p_graph.atoms[a].lone_pairs
                  for a in atoms
              ):
                return list(acc)
            return None
          v = active_1e[v_idx]
          ins, outs = in_1e[v], out_1e[v]
          for perm in set(itertools.permutations(outs)):
            if any(
                s[0] == "bond" and t[0] == "bond" and set(s[1]) == set(t[1])
                for s, t in zip(ins, perm)
            ):
              continue
            step_a = [
                CurvedArrow(s[0], s[1], t[0], t[1], electrons=1)
                for s, t in zip(ins, perm)
            ]
            res = _pair_1e(v_idx + 1, acc + step_a)
            if res is not None:
              return res
          return None

        rad_arrows = _pair_1e(0, [])
        if rad_arrows is not None:
          return rad_arrows

    if has_1e_change:
      return None

    # 2. 2-Electron Polar / Pericyclic / Vinylogous Mode (electrons == 2)
    items = [("dec", e) for e in dec_edges] + [("inc", e) for e in inc_edges]
    rem_deg = {a: 0 for a in atoms}
    for _, (u, v) in items:
      rem_deg[u] += 1
      rem_deg[v] += 1

    cur_diff = {a: 0 for a in atoms}
    q_choices: List[str] = [""] * len(dec_edges)
    p_choices: List[str] = [""] * len(inc_edges)
    solutions: List[Tuple[List[str], List[str]]] = []

    def _backtrack(idx: int) -> None:
      if len(solutions) >= 8:
        return
      if idx == len(items):
        if all(cur_diff[a] == delta_l[a] for a in atoms):
          solutions.append((list(q_choices), list(p_choices)))
        return
      kind, (u, v) = items[idx]
      rem_deg[u] -= 1
      rem_deg[v] -= 1
      delta = 1 if kind == "dec" else -1
      en_u = PAULING_ELECTRONEGATIVITY.get(
          reactant_graph.atoms[u].element, 2.55
      ) + 0.5 * int(reactant_graph.atoms[u].formal_charge)
      en_v = PAULING_ELECTRONEGATIVITY.get(
          reactant_graph.atoms[v].element, 2.55
      ) + 0.5 * int(reactant_graph.atoms[v].formal_charge)
      choices = (
          (u, v)
          if (kind == "dec" and en_u >= en_v)
          or (kind == "inc" and en_u <= en_v)
          else (v, u)
      )
      for chosen in choices:
        other = v if chosen == u else u
        cur_diff[chosen] += delta
        if (
            abs(delta_l[chosen] - cur_diff[chosen]) <= rem_deg[chosen]
            and abs(delta_l[other] - cur_diff[other]) <= rem_deg[other]
        ):
          if kind == "dec":
            q_choices[idx] = chosen
          else:
            p_choices[idx - len(dec_edges)] = chosen
          _backtrack(idx + 1)
        cur_diff[chosen] -= delta
      rem_deg[u] += 1
      rem_deg[v] += 1

    _backtrack(0)

    for q_c, p_c in solutions:
      in_s: Dict[str, List[Tuple[str, Any]]] = {a: [] for a in atoms}
      out_s: Dict[str, List[Tuple[str, Any]]] = {a: [] for a in atoms}
      for a, dl in delta_l.items():
        if dl < 0:
          for _ in range(-dl):
            in_s[a].append(("lone_pair", a))
        elif dl > 0:
          for _ in range(dl):
            out_s[a].append(("atom", a))
      for e, v in zip(dec_edges, q_c):
        in_s[v].append(("bond", e))
      for e, v in zip(inc_edges, p_c):
        out_s[v].append(("bond", e))

      active = [a for a in atoms if in_s[a]]

      def _pair_v(
          v_idx: int, acc: List[CurvedArrow]
      ) -> Optional[List[CurvedArrow]]:
        if v_idx == len(active):
          ok, prod, _ = ReactionVerifierEngine.apply_and_verify_step(
              reactant_graph, acc, allow_multiproton=allow_multiproton
          )
          if ok and _edge_orders(prod) == p_edges:
            if all(
                prod.atoms[a].formal_charge == p_graph.atoms[a].formal_charge
                and prod.atoms[a].lone_pairs == p_graph.atoms[a].lone_pairs
                and (prod.atoms[a].radicals + prod.atoms[a].unpaired_d)
                == (p_graph.atoms[a].radicals + p_graph.atoms[a].unpaired_d)
                for a in atoms
            ):
              return list(acc)
          return None
        v = active[v_idx]
        ins, outs = in_s[v], out_s[v]
        for perm in set(itertools.permutations(outs)):
          if any(
              s[0] == "bond" and t[0] == "bond" and set(s[1]) == set(t[1])
              for s, t in zip(ins, perm)
          ):
            continue
          step_a: List[CurvedArrow] = []
          for s, t in zip(ins, perm):
            dat_flag: Optional[bool] = None
            if (
                s[0] == "lone_pair"
                and t[0] == "bond"
                and not reactant_graph.graph.has_edge(t[1][0], t[1][1])
                and p_graph.graph.has_edge(t[1][0], t[1][1])
            ):
              dat_flag = bool(
                  p_graph.graph[t[1][0]][t[1][1]].get("is_dative", False)
              )
            step_a.append(
                CurvedArrow(
                    s[0], s[1], t[0], t[1], electrons=2, is_dative=dat_flag
                )
            )
          res = _pair_v(v_idx + 1, acc + step_a)
          if res is not None:
            return res
        return None

      ans = _pair_v(0, [])
      if ans is not None:
        return ans

    return None

  @staticmethod
  def solve_arrows_from_mapped_smiles(
      reactant_smiles: str,
      product_smiles: Optional[str] = None,
      allow_multiproton: bool = False,
  ) -> Optional[Tuple[MolecularGraph, MolecularGraph, List[CurvedArrow]]]:
    """Synthesizes a verified CurvedArrow mechanism from atom-mapped SMILES.

    Accepts either `(reactant_smiles, product_smiles)` or a single reaction
    SMILES `"reactants>>products"` in `reactant_smiles`. Automatically promotes
    any implicit hydrogen transferred between mapped heavy atoms to an explicit
    hydrogen node prior to solving.
    """
    if product_smiles is None:
      if ">>" not in reactant_smiles:
        raise ValueError(
            "Provide product_smiles or a 'reactants>>products' reaction SMILES."
        )
      parts = reactant_smiles.split(">>")
      reactant_smiles, product_smiles = parts[0].strip(), parts[-1].strip()

    def _build_mapped_graph(smi: str, name: str) -> MolecularGraph:
      mol = Chem.MolFromSmiles(smi)
      if mol is None:
        mol_raw = Chem.MolFromSmiles(smi, sanitize=False)
        if mol_raw is None:
          raise ValueError(f"RDKit failed to parse SMILES: '{smi}'")
        mol_raw.UpdatePropertyCache(strict=False)
        mol = normalize_dipolar_rdkit_mol(mol_raw)
        Chem.SanitizeMol(mol)
      else:
        mol = normalize_dipolar_rdkit_mol(mol)
      try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
      except Exception:
        pass
      mg = MolecularGraph(name=name)
      idx_to_aid: Dict[int, str] = {}
      for atom in mol.GetAtoms():
        amap = atom.GetAtomMapNum()
        aid = f":{amap}" if amap > 0 else str(atom.GetIdx())
        idx_to_aid[atom.GetIdx()] = aid
        elem = atom.GetSymbol()
        q = atom.GetFormalCharge()
        rad = 0 if elem in D_F_BLOCK_METALS else atom.GetNumRadicalElectrons()
        impl_h = atom.GetTotalNumHs()
        v_free = VALENCE_ELECTRONS.get(elem, 4) if elem != "*" else 0
        explicit_val = 0
        for b in atom.GetBonds():
          if b.GetBondType() == Chem.BondType.DATIVE:
            if b.GetBeginAtomIdx() == atom.GetIdx():
              explicit_val += 2
          else:
            explicit_val += int(round(b.GetBondTypeAsDouble()))
        rem_e = v_free - q - rad - (explicit_val + impl_h) if elem != "*" else 0
        lp = max(0, rem_e // 2)
        unpaired_d = (
            (rem_e % 2) if (elem in D_F_BLOCK_METALS and rem_e > 0) else 0
        )
        mg.add_atom(
            atom_id=aid,
            element=elem,
            lone_pairs=lp,
            radicals=rad,
            unpaired_d=unpaired_d,
            formal_charge=q,
            implicit_h=impl_h,
        )
      for bond in mol.GetBonds():
        u = idx_to_aid[bond.GetBeginAtomIdx()]
        v = idx_to_aid[bond.GetEndAtomIdx()]
        is_dat = bond.GetBondType() == Chem.BondType.DATIVE
        bo = 1 if is_dat else max(1, int(round(bond.GetBondTypeAsDouble())))
        mg.add_bond(
            u,
            v,
            order=bo,
            is_dative=is_dat,
            dative_begin=u if is_dat else None,
            dative_end=v if is_dat else None,
        )
      return mg

    rg = _build_mapped_graph(reactant_smiles, "reactants")
    pg = _build_mapped_graph(product_smiles, "products")
    if set(rg.atoms.keys()) != set(pg.atoms.keys()):
      return None

    # Promote transferred implicit hydrogens between mapped heavy atoms
    donors: List[str] = []
    acceptors: List[str] = []
    for aid in sorted(rg.atoms.keys()):
      dh = int(pg.atoms[aid].implicit_h) - int(rg.atoms[aid].implicit_h)
      if dh < 0:
        donors.extend([aid] * (-dh))
      elif dh > 0:
        acceptors.extend([aid] * dh)
    if len(donors) != len(acceptors):
      return None

    for idx_h, (d_aid, a_aid) in enumerate(zip(donors, acceptors), start=1):
      h_aid = f"H_xfer_{idx_h}"
      rg.atoms[d_aid].implicit_h -= 1
      rg.add_atom(h_aid, "H", lone_pairs=0, radicals=0, formal_charge=0)
      rg.add_bond(d_aid, h_aid, order=1)
      pg.atoms[a_aid].implicit_h -= 1
      pg.add_atom(h_aid, "H", lone_pairs=0, radicals=0, formal_charge=0)
      pg.add_bond(a_aid, h_aid, order=1)

    arrows = ReactionVerifierEngine.solve_arrows_from_mapped_graphs(
        rg, pg, allow_multiproton=allow_multiproton
    )
    if arrows is None:
      return None
    return rg, pg, arrows


def main() -> None:
  parser = argparse.ArgumentParser(
      description=(
          "ElectronGraph: Formal Bond-Electron Matrix & Curved Arrow Mechanism"
          " Verifier."
      )
  )
  parser.add_argument(
      "--smiles",
      type=str,
      default=None,
      help="Construct and inspect a MolecularGraph from a SMILES string.",
  )
  parser.add_argument(
      "--json",
      type=str,
      default=None,
      help="Path to a mechanism step JSON specification to verify.",
  )
  parser.add_argument(
      "--output",
      "-o",
      type=str,
      default=None,
      help="Output JSON file path for the verification report.",
  )
  args = parser.parse_args()

  if args.smiles:
    mg = MolecularGraph.from_smiles(args.smiles)
    payload = mg.to_summary_dict()
    if args.output:
      with open(args.output, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(json.dumps(payload, indent=2))
    return

  if args.json:
    with open(args.json, "r", encoding="utf-8") as f:
      spec = json.load(f)
    mg = MolecularGraph.from_smiles(
        spec["reactant_smiles"],
        make_hs_explicit=bool(spec.get("make_hs_explicit", False)),
    )
    arrows = [
        CurvedArrow(
            source_type=a["source_type"],
            source_id=(
                tuple(a["source_id"])
                if isinstance(a["source_id"], list)
                else a["source_id"]
            ),
            sink_type=a["sink_type"],
            sink_id=(
                tuple(a["sink_id"])
                if isinstance(a["sink_id"], list)
                else a["sink_id"]
            ),
            electrons=int(a.get("electrons", 2)),
            label=str(a.get("label", "")),
        )
        for a in spec.get("arrows", [])
    ]
    valid, prod, audit = ReactionVerifierEngine.apply_and_verify_step(
        mg, arrows, step_name=str(spec.get("step_name", "Mechanism Step"))
    )
    audit["product_smiles"] = prod.to_smiles() if valid else None
    if args.output:
      with open(args.output, "w", encoding="utf-8") as f:
        json.dump(audit, f, indent=2)
    print(json.dumps(audit, indent=2))
    sys.exit(0 if valid else 1)

  parser.error("Provide --smiles or --json to construct or verify a mechanism.")


if __name__ == "__main__":
  main()
