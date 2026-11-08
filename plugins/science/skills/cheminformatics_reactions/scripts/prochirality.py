# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
# ]
# ///

"""Prochirality, Topicity (pro-R/pro-S, Re/Si), & Reaction Stereogenesis."""

import argparse
import collections
import json
import sys
from typing import Any, Dict, List, Sequence, Set, Tuple, Union

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdCIPLabeler

_PT = Chem.GetPeriodicTable()
_CW = Chem.ChiralType.CHI_TETRAHEDRAL_CW
_CCW = Chem.ChiralType.CHI_TETRAHEDRAL_CCW


def _heavier_isotope(z: int, iso: int = 0) -> int:
  """Returns the next valid heavier NIST isotope for CIP Rule 1b perturbation."""
  base = iso if iso > 0 else round(_PT.GetMostCommonIsotopeMass(z))
  base_mass = (
      _PT.GetMassForIsotope(z, iso)
      if iso > 0
      else _PT.GetMostCommonIsotopeMass(z)
  )
  for cand in range(base + 1, base + 12):
    if _PT.GetMassForIsotope(z, cand) > base_mass:
      return cand
  return base + 1


def _invert_tetrahedral(mol: Chem.Mol) -> Chem.Mol:
  """Inverts all specified tetrahedral stereocenters to form the mirror image."""
  m = Chem.Mol(mol)
  for a in m.GetAtoms():
    t = a.GetChiralTag()
    if t == _CW:
      a.SetChiralTag(_CCW)
    elif t == _CCW:
      a.SetChiralTag(_CW)
  return m


def _perturb_atom(ref_mol: Chem.Mol, idx: int) -> Tuple[Chem.Mol, str, str]:
  """Perturbs `idx` with a heavier isotope and returns (mol, smi, mirror_smi)."""
  rw = Chem.Mol(ref_mol)
  a = rw.GetAtomWithIdx(idx)
  a.SetIsotope(_heavier_isotope(a.GetAtomicNum(), a.GetIsotope()))
  if rw.GetNumConformers() > 0:
    Chem.AssignStereochemistryFrom3D(rw)
  Chem.AssignStereochemistry(rw, cleanIt=True, force=True)
  rdCIPLabeler.AssignCIPLabels(rw)
  s_orig = Chem.MolToSmiles(rw, isomericSmiles=True, canonical=True)
  s_mir = Chem.MolToSmiles(
      _invert_tetrahedral(rw), isomericSmiles=True, canonical=True
  )
  return rw, s_orig, s_mir


def detect_false_chirality(smiles: str) -> List[int]:
  """Returns atom indices carrying @/@@ in SMILES whose ligands are homotopic."""
  raw = Chem.MolFromSmiles(smiles, sanitize=False)
  params = Chem.SmilesParserParams()
  params.removeHs = False
  clean = Chem.MolFromSmiles(smiles, params)
  if raw is None or clean is None:
    return []
  tagged = [
      a.GetIdx() for a in raw.GetAtoms() if a.GetChiralTag() in (_CW, _CCW)
  ]
  if not tagged:
    return []
  Chem.AssignStereochemistry(clean, cleanIt=True, force=True)
  pot = {
      s.centeredOn
      for s in Chem.FindPotentialStereo(clean)
      if s.type == Chem.StereoType.Atom_Tetrahedral
  }
  return [
      i
      for i in tagged
      if clean.GetAtomWithIdx(i).GetChiralTag() not in (_CW, _CCW)
      and i not in pot
  ]


def analyze_prochirality(mol_or_smiles: Union[str, Chem.Mol]) -> Dict[str, Any]:
  """Perceives sp3 (pro-R/pro-S) and sp2 (Re/Si) prochirality and topicity."""
  mol = (
      Chem.MolFromSmiles(mol_or_smiles)
      if isinstance(mol_or_smiles, str)
      else Chem.Mol(mol_or_smiles)
  )
  if mol is None:
    raise ValueError(f"Invalid SMILES: {mol_or_smiles!r}")
  mol_h = Chem.AddHs(mol)
  Chem.AssignStereochemistry(mol_h, cleanIt=True, force=True)
  const_ranks = list(
      Chem.CanonicalRankAtoms(mol_h, breakTies=False, includeChirality=False)
  )

  # Embed a 3D reference conformer so polycyclic/bridged stereocenters have
  # physically consistent 3D relative configurations before perturbation.
  ref_mol = Chem.Mol(mol_h)
  params = AllChem.ETKDGv3()
  params.randomSeed = 42
  if AllChem.EmbedMolecule(ref_mol, params) == 0:
    try:
      AllChem.MMFFOptimizeMolecule(ref_mol)
    except Exception:
      pass
    Chem.AssignStereochemistryFrom3D(ref_mol)
  else:
    for a in ref_mol.GetAtoms():
      if (
          a.GetDegree() == 4
          and a.GetTotalValence() == 4
          and a.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
      ):
        a.SetChiralTag(_CW)

  # 1. Whole-molecule stereotopic perturbation (homotopic/enantiotopic)
  perturbed = [_perturb_atom(ref_mol, i) for i in range(mol_h.GetNumAtoms())]
  achiral_classes: Dict[str, Set[Tuple[str, ...]]] = collections.defaultdict(
      set
  )
  chiral_classes: Dict[str, Set[str]] = collections.defaultdict(set)
  for i, (_, s_orig, s_mir) in enumerate(perturbed):
    sym = mol_h.GetAtomWithIdx(i).GetSymbol()
    achiral_classes[sym].add(tuple(sorted([s_orig, s_mir])))
    chiral_classes[sym].add(s_orig)

  # 2. Tetrahedral (sp3) prochiral centers (Hanson's Rule, IUPAC P-92.5.1)
  sp3_centers = []
  for a in mol_h.GetAtoms():
    if a.GetDegree() != 4 or a.GetTotalValence() != 4:
      continue
    c_idx = a.GetIdx()
    by_rank: Dict[int, List[int]] = collections.defaultdict(list)
    for nb in a.GetNeighbors():
      by_rank[const_ranks[nb.GetIdx()]].append(nb.GetIdx())
    for pair in (idxs for idxs in by_rank.values() if len(idxs) == 2):
      ia, ib = pair
      rw_a, smi_a, mir_a = perturbed[ia]
      rw_b, smi_b, _ = perturbed[ib]
      if smi_a == smi_b:
        top = "homotopic"
      elif mir_a == smi_b:
        top = "enantiotopic"
      else:
        top = "diastereotopic"
      cip_a = rw_a.GetAtomWithIdx(c_idx).GetPropsAsDict().get("_CIPCode")
      cip_b = rw_b.GetAtomWithIdx(c_idx).GetPropsAsDict().get("_CIPCode")
      sp3_centers.append({
          "center_idx": c_idx,
          "topicity": top,
          "is_prochiral": top != "homotopic",
          "ligands": [
              {"atom_idx": ia, "descriptor": f"pro-{cip_a}" if cip_a else top},
              {"atom_idx": ib, "descriptor": f"pro-{cip_b}" if cip_b else top},
          ],
      })

  # 3. Trigonal planar (sp2) prochiral centers (Re/Si faces, IUPAC P-92.5.2)
  sp2_centers = []
  for a in mol_h.GetAtoms():
    c_idx = a.GetIdx()
    if a.GetDegree() != 3 or a.GetIsAromatic():
      continue
    if not any(b.GetBondTypeAsDouble() >= 2.0 for b in a.GetBonds()) and (
        a.GetFormalCharge() <= 0
    ):
      continue
    probes = {}
    for name, tag in (("CW", _CW), ("CCW", _CCW)):
      rw = Chem.RWMol(ref_mol)
      d_idx = rw.AddAtom(Chem.Atom(0))
      rw.AddBond(c_idx, d_idx, Chem.BondType.SINGLE)
      rw.GetAtomWithIdx(c_idx).SetChiralTag(tag)
      rw.UpdatePropertyCache(strict=False)
      if rw.GetNumConformers() > 0:
        Chem.AssignStereochemistryFrom3D(rw, -1, False)
      for s in Chem.FindPotentialStereo(rw):
        if s.type == Chem.StereoType.Atom_Tetrahedral and s.centeredOn != c_idx:
          a_s = rw.GetAtomWithIdx(s.centeredOn)
          if a_s.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            a_s.SetChiralTag(_CW)
      Chem.AssignStereochemistry(rw, cleanIt=True, force=True)
      rdCIPLabeler.AssignCIPLabels(rw)
      ca = rw.GetAtomWithIdx(c_idx)
      cip = ca.GetPropsAsDict().get("_CIPCode")
      s_orig = Chem.MolToSmiles(rw, isomericSmiles=True, canonical=True)
      s_mir = Chem.MolToSmiles(
          _invert_tetrahedral(rw), isomericSmiles=True, canonical=True
      )
      orig_nbs = [nb for nb in ca.GetNeighbors() if nb.GetIdx() != d_idx]
      nbs = sorted(
          orig_nbs,
          key=lambda x: (
              int(x.GetPropsAsDict().get("_CIPRank", -1)),
              x.GetAtomicNum(),
              x.GetIsotope(),
              const_ranks[x.GetIdx()],
          ),
          reverse=True,
      )
      prio_idxs = [nb.GetIdx() for nb in nbs]
      if cip in ("R", "S", "r", "s") and len(prio_idxs) == 3:
        orig_idxs = [nb.GetIdx() for nb in orig_nbs]
        inv = 0
        for i in range(3):
          for j in range(i + 1, 3):
            if orig_idxs.index(prio_idxs[i]) > orig_idxs.index(prio_idxs[j]):
              inv += 1
        expected_mod = 0 if cip in ("R", "r") else 1
        if inv % 2 != expected_mod:
          prio_idxs[1], prio_idxs[2] = prio_idxs[2], prio_idxs[1]
      probes[name] = (
          cip,
          s_orig,
          s_mir,
          prio_idxs,
      )
    (cw_cip, cw_s, cw_m, prio), (ccw_cip, ccw_s, _, _) = (
        probes["CW"],
        probes["CCW"],
    )
    if (cw_s == ccw_s and cw_cip == ccw_cip) or (
        cw_cip is None and ccw_cip is None
    ):
      face_top, faces = "homotopic", {"Re": None, "Si": None}
    else:
      face_top = (
          "enantiotopic"
          if (cw_m == ccw_s and cw_cip in ("R", "S"))
          else "diastereotopic"
      )
      faces = (
          {"Re": "CW", "Si": "CCW"}
          if cw_cip in ("S", "s")
          else {"Re": "CCW", "Si": "CW"}
      )
    sp2_centers.append({
        "center_idx": c_idx,
        "facial_topicity": face_top,
        "is_prochiral": face_top != "homotopic",
        "priority_1_to_3_neighbors": prio,
        "slot4_face_tags": faces,
    })

  return {
      "smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
      "tetrahedral_centers": sp3_centers,
      "trigonal_centers": sp2_centers,
      "achiral_classes_by_element": {
          k: len(v) for k, v in achiral_classes.items()
      },
      "chiral_classes_by_element": {
          k: len(v) for k, v in chiral_classes.items()
      },
  }


def verify_reaction_stereogenesis(
    precursors: Union[str, Sequence[str]],
    products: Union[str, Sequence[str]],
    chiral_catalyst: bool = False,
) -> Dict[str, Any]:
  """Validates reaction stereogenesis, false chirality, and asymmetric induction."""
  precs = [precursors] if isinstance(precursors, str) else list(precursors)
  prods = [products] if isinstance(products, str) else list(products)
  errors, warnings, notes = [], [], []

  def _inspect(smi_list: List[str], role: str):
    out = []
    for smi in smi_list:
      false_idxs = detect_false_chirality(smi)
      if false_idxs:
        errors.append(
            f"{role} '{smi}' has false @/@@ chirality on homotopic atom(s)"
            f" {false_idxs}."
        )
      m = Chem.MolFromSmiles(smi)
      if m is None:
        errors.append(f"Invalid {role} SMILES: '{smi}'")
        continue
      Chem.AssignStereochemistry(m, cleanIt=True, force=True)
      pot = [
          s
          for s in Chem.FindPotentialStereo(m)
          if s.type == Chem.StereoType.Atom_Tetrahedral
      ]
      spec = [s for s in pot if s.specified == Chem.StereoSpecified.Specified]
      s_can = Chem.MolToSmiles(m, isomericSmiles=True, canonical=True)
      s_mir = Chem.MolToSmiles(
          _invert_tetrahedral(m), isomericSmiles=True, canonical=True
      )
      out.append({
          "smiles": smi,
          "potential_stereocenters": len(pot),
          "specified_stereocenters": len(spec),
          "optically_active": s_can != s_mir,
          "meso": len(spec) >= 2 and s_can == s_mir,
      })
    return out

  r_info = _inspect(precs, "Precursor")
  p_info = _inspect(prods, "Product")
  if errors:
    return {"valid": False, "errors": errors, "warnings": warnings}

  delta_pot = sum(x["potential_stereocenters"] for x in p_info) - sum(
      x["potential_stereocenters"] for x in r_info
  )
  r_opt = any(x["optically_active"] for x in r_info)
  p_opt = any(x["optically_active"] for x in p_info)
  r_meso = any(x["meso"] for x in r_info)

  mode = "none"
  if delta_pot > 0 or (not r_opt and p_opt):
    if r_opt:
      mode = "diastereogenic"
      notes.append("Substrate-controlled diastereogenic stereocenter creation.")
    elif r_meso and p_opt:
      mode = "meso_desymmetrization"
    elif p_opt:
      mode = "enantioselective_de_novo"
      if not chiral_catalyst:
        warnings.append(
            "Optically inactive precursors cannot yield a single enantiomer"
            " product without a chiral catalyst/auxiliary (achiral conditions"
            " yield a racemate)."
        )
    else:
      mode = "racemic_de_novo"
  elif delta_pot < 0:
    mode = "chirality_destruction"

  return {
      "valid": True,
      "stereogenesis_mode": mode,
      "delta_potential_stereocenters": delta_pot,
      "precursors_optically_active": r_opt,
      "products_optically_active": p_opt,
      "errors": errors,
      "warnings": warnings,
      "notes": notes,
  }


def main() -> None:
  ap = argparse.ArgumentParser(description="Prochirality & Stereogenesis CLI.")
  ap.add_argument("--smiles", "-s", help="Molecule SMILES to analyze.")
  ap.add_argument("--precursor", "--prec", nargs="+", help="Precursor SMILES.")
  ap.add_argument("--product", "--prod", nargs="+", help="Product SMILES.")
  ap.add_argument("--chiral-catalyst", action="store_true")
  args = ap.parse_args()
  if args.precursor and args.product:
    res = verify_reaction_stereogenesis(
        args.precursor, args.product, args.chiral_catalyst
    )
  elif args.smiles:
    res = analyze_prochirality(args.smiles)
  else:
    ap.error("Provide --smiles or (--precursor and --product).")
    return
  print(json.dumps(res, indent=2))
  sys.exit(0 if res.get("valid", True) else 1)


if __name__ == "__main__":
  main()
