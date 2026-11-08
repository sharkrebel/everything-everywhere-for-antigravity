# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "pandas",
#     "scipy",
#     "scikit-learn",
#     "matplotlib",
# ]
# ///

"""Organometallic Structure Standardization, Dative Bonds & Electron Counting Engine.

Handles transition metal complexes, coordination chemistry representations,
dative bond normalization (BondType.DATIVE), oxidation state assignment, and
18-electron rule / d^n valence counting in RDKit.
"""

import argparse
import json
import os
import re
import subprocess
import sys
from typing import Any, Dict, List, Optional, Set, Tuple

import rdkit
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

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

from draw_molecules import regularize_nonconvex_small_rings

# Periodic Table & Transition Metal Classification
TRANSITION_METALS: Set[int] = {
    # 3d metals
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,  # Sc to Zn
    # 4d metals
    39,
    40,
    41,
    42,
    43,
    44,
    45,
    46,
    47,
    48,  # Y to Cd
    # 5d metals & Lanthanides
    57,
    58,
    59,
    60,
    61,
    62,
    63,
    64,
    65,
    66,
    67,
    68,
    69,
    70,
    71,  # La-Lu
    72,
    73,
    74,
    75,
    76,
    77,
    78,
    79,
    80,  # Hf to Hg
    # Actinides & 6d metals
    89,
    90,
    91,
    92,
    93,
    94,
    95,
    96,
    97,
    98,
    99,
    100,
    101,
    102,
    103,  # Ac-Lr
    104,
    105,
    106,
    107,
    108,
    109,
    110,
    111,
    112,  # Rf to Cn
}

MAIN_GROUP_METALS: Set[int] = {
    3,
    11,
    12,
    13,
    19,
    20,
    31,
    37,
    38,
    49,
    50,
    55,
    56,
    81,
    82,
    83,  # Li, Na, Mg, Al, K, Ca, Ga, Rb, Sr, In, Sn, Cs, Ba, Tl, Pb, Bi
}

ALL_METALS: Set[int] = TRANSITION_METALS | MAIN_GROUP_METALS

# Valence / Group Numbers for Transition Metals (Column in Periodic Table)
METAL_GROUP_NUMBERS: Dict[str, int] = {
    "Sc": 3,
    "Y": 3,
    "La": 3,
    "Lu": 3,
    "Ti": 4,
    "Zr": 4,
    "Hf": 4,
    "V": 5,
    "Nb": 5,
    "Ta": 5,
    "Cr": 6,
    "Mo": 6,
    "W": 6,
    "Mn": 7,
    "Tc": 7,
    "Re": 7,
    "Fe": 8,
    "Ru": 8,
    "Os": 8,
    "Co": 9,
    "Rh": 9,
    "Ir": 9,
    "Ni": 10,
    "Pd": 10,
    "Pt": 10,
    "Cu": 11,
    "Ag": 11,
    "Au": 11,
    "Zn": 12,
    "Cd": 12,
    "Hg": 12,
}

# Neutral / L-type donor atom patterns (e.g. phosphines, amines, carbonyls, ethers, thioethers, carbenes)
L_TYPE_DONORS: Set[int] = {7, 8, 15, 16, 33, 34}  # N, O, P, S, As, Se


def is_metal_atom(atom: Chem.Atom) -> bool:
  """Returns True if the atom is a transition metal or main group metal."""
  return atom.GetAtomicNum() in ALL_METALS


def find_metal_centers(mol: Chem.Mol) -> List[int]:
  """Returns the atom indices of all metal atoms in the molecule."""
  return [atom.GetIdx() for atom in mol.GetAtoms() if is_metal_atom(atom)]


def is_l_type_carbon(atom: Chem.Atom) -> bool:
  """Returns True if the carbon is an L-type donor (CO or carbene), False for alkyl/aryl."""
  if atom.GetAtomicNum() != 6:
    return False
  # Check if CO: carbon bonded to oxygen
  for nbr in atom.GetNeighbors():
    if nbr.GetAtomicNum() == 8:
      return True
  # Check if carbene (e.g. NHC): carbon bonded to two nitrogen atoms
  n_count = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 7)
  if n_count >= 2:
    return True
  return False


def standardize_dative_bonds(
    mol: Chem.Mol,
    donor_atomic_nums: Optional[Set[int]] = None,
) -> Chem.Mol:
  """Converts covalent bonds between neutral donor atoms and metal centers into directional DATIVE bonds.

  Direction: Donor Atom -> Metal Center.
  Preserves the exact neighbor bond ordering around each metal center so that
  square-planar (@SP1/@SP2/@SP3), trigonal-bipyramidal (@TB), octahedral (@OH),
  and tetrahedral (@/@@) stereochemical permutations are never inverted when
  converting single bonds to directed dative arrows.
  """
  if donor_atomic_nums is None:
    donor_atomic_nums = L_TYPE_DONORS

  rwmol = Chem.RWMol(mol)
  metal_indices = find_metal_centers(rwmol)

  for metal_idx in metal_indices:
    metal_atom = rwmol.GetAtomWithIdx(metal_idx)
    ordered_metal_bonds: List[Tuple[int, int, Chem.BondType]] = []
    any_converted = False

    for bond in metal_atom.GetBonds():
      donor_atom = bond.GetOtherAtom(metal_atom)
      donor_idx = donor_atom.GetIdx()
      btype = bond.GetBondType()

      if is_metal_atom(donor_atom):
        ordered_metal_bonds.append(
            (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), btype)
        )
        continue

      convert_to_dative = False
      if donor_atom.GetAtomicNum() == 6:
        if is_l_type_carbon(donor_atom) and btype in (
            Chem.BondType.SINGLE,
            Chem.BondType.DATIVE,
            Chem.BondType.UNSPECIFIED,
        ):
          convert_to_dative = True
      elif (
          donor_atom.GetAtomicNum() in donor_atomic_nums
          and donor_atom.GetFormalCharge() >= 0
          and btype
          in (
              Chem.BondType.SINGLE,
              Chem.BondType.DATIVE,
              Chem.BondType.UNSPECIFIED,
          )
      ):
        convert_to_dative = True

      if convert_to_dative:
        any_converted = True
        ordered_metal_bonds.append((donor_idx, metal_idx, Chem.BondType.DATIVE))
      else:
        ordered_metal_bonds.append((donor_idx, metal_idx, btype))

    if any_converted:
      # Remove and re-add all bonds on this metal in their exact original sequence
      # so that @SP1/@SP2/@OH neighbor permutations remain invariant.
      for b_begin, b_end, _ in ordered_metal_bonds:
        rwmol.RemoveBond(b_begin, b_end)
      for b_begin, b_end, new_btype in ordered_metal_bonds:
        rwmol.AddBond(b_begin, b_end, new_btype)

  try:
    Chem.SanitizeMol(rwmol)
  except Exception:
    rwmol.UpdatePropertyCache(strict=False)

  out_mol = rwmol.GetMol()
  try:
    Chem.AssignStereochemistry(out_mol, force=True, cleanIt=True)
  except Exception:
    pass
  return out_mol


def disconnect_metal_complex(
    mol: Chem.Mol,
) -> Tuple[Chem.Mol, List[Dict[str, Any]]]:
  """Disconnects all metal-ligand coordination bonds to yield isolated metal centers and free ligands.

  Useful for separate ligand-based virtual screening, descriptor calculation,
  and salt stripping.
  """
  rwmol = Chem.RWMol(mol)
  metal_indices = find_metal_centers(mol)
  bonds_to_remove = []

  for midx in metal_indices:
    matom = rwmol.GetAtomWithIdx(midx)
    for nbr in matom.GetNeighbors():
      bonds_to_remove.append((midx, nbr.GetIdx()))

  for a1, a2 in bonds_to_remove:
    rwmol.RemoveBond(a1, a2)

  try:
    Chem.SanitizeMol(rwmol)
  except Exception:
    rwmol.UpdatePropertyCache(strict=False)

  disc_mol = rwmol.GetMol()
  frags = Chem.GetMolFrags(disc_mol, asMols=True, sanitizeFrags=False)

  frag_metadata = []
  for f in frags:
    f_metals = find_metal_centers(f)
    is_metal_frag = len(f_metals) > 0
    smi = Chem.MolToSmiles(f)
    frag_metadata.append({
        "smiles": smi,
        "is_metal_core": is_metal_frag,
        "num_atoms": f.GetNumAtoms(),
        "formula": rdkit.Chem.rdMolDescriptors.CalcMolFormula(f),
    })

  return disc_mol, frag_metadata


def count_electrons_and_oxidation_state(mol: Chem.Mol) -> List[Dict[str, Any]]:
  """Calculates coordination number (CN), oxidation state (OS), d^n count, and valence electron count (VEC).

  Uses the Green MLX covalent-bond classification method:
  - X-ligands (1e- anionic donors): Halides (Cl, Br, I, F), Alkyls (Me, Et, Ar),
  Hydrides (H), Alkoxides (OR), Amides (NR2).
  - L-ligands (2e- neutral donors): Phosphines (PR3), Amines (NR3), Carbonyls
  (CO), Ethers (OR2), Carbenes (NHC), Alkenes.
  - Z-ligands (0e- neutral Lewis acid acceptors): Boranes (BR3), Alanes (AlR3).
  """
  results = []
  metal_indices = find_metal_centers(mol)

  for midx in metal_indices:
    matom = mol.GetAtomWithIdx(midx)
    symbol = matom.GetSymbol()
    group_num = METAL_GROUP_NUMBERS.get(symbol, 0)

    x_ligands = 0
    l_ligands = 0
    coordination_number = 0
    ligand_details = []

    for bond in matom.GetBonds():
      coordination_number += 1
      nbr = bond.GetOtherAtom(matom)
      btype = bond.GetBondType()

      # Classify by bond type and neighbor characteristics
      if btype == Chem.BondType.DATIVE:
        l_ligands += 1
        ligand_type = "L (2e- Neutral Dative)"
      elif nbr.GetAtomicNum() in (9, 17, 35, 53):  # Halides (F, Cl, Br, I)
        x_ligands += 1
        ligand_type = "X (1e- Anionic Halide)"
      elif nbr.GetAtomicNum() == 1:  # Hydride
        x_ligands += 1
        ligand_type = "X (1e- Anionic Hydride)"
      elif (
          btype == Chem.BondType.SINGLE and nbr.GetAtomicNum() == 6
      ):  # Alkyl / Aryl
        x_ligands += 1
        ligand_type = "X (1e- Anionic Alkyl/Aryl)"
      elif btype == Chem.BondType.SINGLE and nbr.GetAtomicNum() in (
          7,
          8,
          15,
          16,
      ):
        # Anionic or neutral single bond
        if nbr.GetFormalCharge() < 0:
          x_ligands += 1
          ligand_type = "X (1e- Anionic Heteroatom)"
        else:
          l_ligands += 1
          ligand_type = "L (2e- Neutral Heteroatom)"
      else:
        # Default fallback
        l_ligands += 1
        ligand_type = "L (2e- Neutral)"

      ligand_details.append({
          "neighbor_symbol": nbr.GetSymbol(),
          "neighbor_idx": nbr.GetIdx(),
          "bond_type": str(btype),
          "classification": ligand_type,
      })

    # Formal Oxidation State = Total X-type ligands + net formal charge on metal
    formal_charge = matom.GetFormalCharge()
    oxidation_state = x_ligands + formal_charge

    # d-electron count: d^n = Group Number - Oxidation State
    d_electrons = max(0, group_num - oxidation_state) if group_num > 0 else None

    # Total Valence Electron Count (VEC) = d^n + 2*L + 2*X
    if d_electrons is not None:
      valence_electron_count = d_electrons + (2 * l_ligands) + (2 * x_ligands)
    else:
      valence_electron_count = None

    # Stability rules: 18-electron (octahedral / tetrahedral), 16-electron (square planar d8)
    if valence_electron_count == 18:
      electron_rule_status = "Satisfies 18-Electron Rule (Closed Shell)"
    elif valence_electron_count == 16 and d_electrons == 8:
      electron_rule_status = (
          "Satisfies 16-Electron Rule (Stable d8 Square Planar)"
      )
    elif valence_electron_count is not None and valence_electron_count < 18:
      electron_rule_status = (
          f"Coordinatively Unsaturated ({valence_electron_count}e- < 18e-)"
      )
    else:
      electron_rule_status = f"{valence_electron_count}e- Complex"

    results.append({
        "metal_symbol": symbol,
        "metal_idx": midx,
        "group_number": group_num,
        "formal_charge": formal_charge,
        "coordination_number": coordination_number,
        "x_type_ligands": x_ligands,
        "l_type_ligands": l_ligands,
        "oxidation_state": oxidation_state,
        "d_electron_count": (
            f"d^{d_electrons}" if d_electrons is not None else "N/A"
        ),
        "valence_electron_count": valence_electron_count,
        "electron_rule_status": electron_rule_status,
        "ligands": ligand_details,
    })

  return results


def render_organometallic_depiction(
    mol: Chem.Mol,
    output_path: str,
    render_style: str = "horizon_perspective",
    legend: str = "",
    width: int = 450,
    height: int = 400,
) -> str:
  """Renders a 2D depiction of an organometallic complex with optional 3D horizon-perspective dashed rear bonds."""
  m = Chem.Mol(mol)
  try:
    Chem.AssignStereochemistry(m, force=True, cleanIt=True)
  except Exception:
    pass
  if m.GetNumConformers() == 0:
    rdDepictor.Compute2DCoords(m)
    regularize_nonconvex_small_rings(m)

  dashed_bonds: List[int] = []
  if render_style == "horizon_perspective" and m.GetNumConformers() > 0:
    conf = m.GetConformer()
    for midx in find_metal_centers(m):
      matom = m.GetAtomWithIdx(midx)
      if len(matom.GetBonds()) < 4:
        continue
      m_pos = conf.GetAtomPosition(midx)
      candidates = []
      for b in matom.GetBonds():
        other_idx = b.GetOtherAtomIdx(midx)
        opos = conf.GetAtomPosition(other_idx)
        dy = opos.y - m_pos.y
        if dy > 0.1:
          candidates.append((b.GetIdx(), dy))
      candidates.sort(key=lambda item: item[1], reverse=True)
      for b_idx, _ in candidates[:2]:
        dashed_bonds.append(b_idx)

  d2d = rdMolDraw2D.MolDraw2DSVG(width, height)
  d2d.drawOptions().bondLineWidth = 2.0
  d2d.DrawMolecule(m, legend=legend)
  d2d.FinishDrawing()
  svg = d2d.GetDrawingText()

  for b_idx in dashed_bonds:
    pattern = (
        rf"(<path[^>]*class=[\x27\x22][^\x27\x22]*\bbond-{b_idx}\b"
        rf"[^\x27\x22]*[\x27\x22][^>]*style=[\x27\x22][^\x27\x22]*fill:none"
        rf"[^\x27\x22]*)([\x27\x22])"
    )

    def _repl(match):
      orig = match.group(1)
      quote = match.group(2)
      if "stroke-dasharray" not in orig:
        return orig + ";stroke-dasharray:6,4" + quote
      return orig + quote

    svg = re.sub(pattern, _repl, svg)

  abs_out = os.path.abspath(output_path)
  os.makedirs(os.path.dirname(abs_out) or ".", exist_ok=True)
  if abs_out.lower().endswith(".svg"):
    with open(abs_out, "w", encoding="utf-8") as f:
      f.write(svg)
  else:
    svg_tmp = abs_out + ".tmp.svg"
    with open(svg_tmp, "w", encoding="utf-8") as f:
      f.write(svg)
    try:
      subprocess.run(
          ["convert", "-density", "300", svg_tmp, abs_out],
          check=True,
          capture_output=True,
      )
    except Exception:
      d_png = rdMolDraw2D.MolDraw2DCairo(width, height)
      d_png.drawOptions().bondLineWidth = 2.0
      d_png.DrawMolecule(m, legend=legend)
      d_png.FinishDrawing()
      with open(abs_out, "wb") as f:
        f.write(d_png.GetDrawingText())
    finally:
      if os.path.exists(svg_tmp):
        os.remove(svg_tmp)
  return abs_out


def main():
  parser = argparse.ArgumentParser(
      description="Organometallic Standardization & Electron Counting Engine."
  )
  parser.add_argument(
      "--smiles", type=str, required=True, help="Organometallic complex SMILES."
  )
  parser.add_argument(
      "--dative",
      action="store_true",
      help="Convert to standardized dative bonds.",
  )
  parser.add_argument(
      "--disconnect",
      action="store_true",
      help="Disconnect metal and return free ligands.",
  )
  parser.add_argument(
      "--electron_count",
      action="store_true",
      help="Print d^n and 18-electron rule analysis.",
  )
  parser.add_argument(
      "--json", action="store_true", help="Output summary as JSON."
  )
  parser.add_argument(
      "--output",
      "--output-json",
      dest="output_path",
      type=str,
      default=None,
      help="Save JSON analysis output to a file.",
  )
  parser.add_argument(
      "--render",
      type=str,
      nargs="?",
      const="organometallic_depiction.svg",
      default=None,
      help="Render 2D depiction of standardized complex to SVG/PNG path.",
  )
  parser.add_argument(
      "--render_style",
      "--render-style",
      dest="render_style",
      choices=["standard", "horizon_perspective"],
      default="horizon_perspective",
      help="2D coordination depiction style (default: horizon_perspective).",
  )
  args = parser.parse_args()

  mol = Chem.MolFromSmiles(args.smiles)
  if mol is None:
    # Try relaxed parsing
    mol = Chem.MolFromSmiles(args.smiles, sanitize=False)
    if mol is not None:
      mol.UpdatePropertyCache(strict=False)

  if mol is None:
    print(f"Error: Could not parse SMILES '{args.smiles}'", file=sys.stderr)
    sys.exit(1)

  std_mol = standardize_dative_bonds(mol)
  elec_info = count_electrons_and_oxidation_state(std_mol)

  if args.render:
    rendered_path = render_organometallic_depiction(
        std_mol,
        output_path=args.render,
        render_style=args.render_style,
    )
    print(
        f"[+] Rendered organometallic depiction ({args.render_style}) to:"
        f" {rendered_path}"
    )

  if args.disconnect:
    disc_mol, frags = disconnect_metal_complex(std_mol)
    print(f"Disconnected SMILES: {Chem.MolToSmiles(disc_mol)}")
    for i, frag in enumerate(frags):
      print(
          f"  Frag {i+1}: {frag['smiles']} ({frag['formula']},"
          f" Metal={frag['is_metal_core']})"
      )
    return

  payload = {
      "input_smiles": args.smiles,
      "standardized_dative_smiles": Chem.MolToSmiles(std_mol),
      "metals": elec_info,
  }

  if args.output_path:
    with open(args.output_path, "w") as f:
      json.dump(payload, f, indent=2)
    print(f"[+] Saved organometallic analysis JSON to: {args.output_path}")

  if args.json:
    print(json.dumps(payload, indent=2))
    return

  print("\n=== Organometallic Complex Analysis ===")
  print(f"Input SMILES:        {args.smiles}")
  print(f"Dative SMILES:       {Chem.MolToSmiles(std_mol)}")
  for m in elec_info:
    print(
        f"\nMetal Center:        {m['metal_symbol']} (Idx {m['metal_idx']},"
        f" Group {m['group_number']})"
    )
    print(f"Coordination Number: {m['coordination_number']}")
    print(f"Oxidation State:     {m['metal_symbol']}({m['oxidation_state']})")
    print(f"d-Electron Count:    {m['d_electron_count']}")
    print(
        f"Valence Electrons:   {m['valence_electron_count']}e-"
        f" ({m['electron_rule_status']})"
    )
    print(
        f"Ligand Breakdown:    {m['x_type_ligands']} X-type (1e-),"
        f" {m['l_type_ligands']} L-type (2e-)"
    )
    for l in m["ligands"]:
      print(
          f"  • {l['neighbor_symbol']} (Idx {l['neighbor_idx']}):"
          f" {l['classification']}"
      )


if __name__ == "__main__":
  main()
