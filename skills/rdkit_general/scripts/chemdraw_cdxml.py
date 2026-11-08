# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy",
#     "pandas",
#     "rdkit",
# ]
# ///
"""ChemDraw CDXML/CDX Converter, Parser, and Publication Style Auditor.

Provides headless interoperability with ChemDraw (.cdxml, .cdx) files:
1. Export (`--smiles` / `--molfile` -> `.cdxml`):
   - Enforces ACS Document 1996 (`acs1996`, 14.4 pt bond length, 0.6 pt line
     width, 2.0 pt bold width, Arial 10 pt atom labels, Times New Roman 10 pt
     captions), RSC (`rsc`, 15.0 pt), or Wiley (`wiley`, 17.0 pt) styling.
   - Inverts Cartesian y-coordinates (`y_cdxml = -y_cartesian * scale`).
   - Kekulizes aromatic rings (`Order="2"`) so desktop ChemDraw renders
     alternating double bonds rather than blank `Order="1.5"` query rings.
   - Injects explicit `NumHydrogens` on all non-carbon atoms with implicit
     hydrogens (`-NH2`, `-OH`, `-NH-`, `-SH`) to eliminate bright red valence
     warning boxes in ChemDraw Desktop and ChemDraw JS.
2. Parse (`--parse` / `--format {smiles,sdf,json,rxn}`):
   - Extracts RDKit molecules (`rdChemDraw.MolsFromChemDrawFile`), reactions
     (`rdChemDraw.ReactionsFromChemDrawFile`), and XML metadata from `.cdxml`
     and `.cdx` files, with a pure-XML fallback parser.
3. Inspect (`--inspect`):
   - Audits `.cdxml` documents for heteroatom valence satisfaction
     (`NumHydrogens` vs bond orders), unkekulized `Order="1.5"` bonds, and dual
     reaction arrow tags (`<graphic GraphicType="Line" ArrowType="FullHead">` +
     `<step>/<arrow>`).
"""

import argparse
from io import StringIO
import json
import math
import os
import sys
from typing import Any, Dict, List, Sequence, Tuple, Union
import xml.etree.ElementTree as ET

from rdkit import Chem
from rdkit.Chem import AllChem

try:
  from rdkit.Chem import rdChemDraw

  _HAS_RDCHEMDRAW = True
except ImportError:
  rdChemDraw = None
  _HAS_RDCHEMDRAW = False

_script_dir = os.path.dirname(os.path.abspath(__file__))
if _script_dir not in sys.path:
  sys.path.insert(0, _script_dir)

from draw_molecules import regularize_nonconvex_small_rings


STYLE_PRESETS: Dict[str, Dict[str, Any]] = {
    "acs1996": {
        "name": "ACS Document 1996",
        "bond_length": 14.4,  # pt (0.2 in / 5.08 mm)
        "line_width": 0.6,
        "bold_width": 2.0,
        "hash_spacing": 2.5,
        "margin_width": 1.6,
        "bond_spacing": 18,
        "chain_angle": 120,
        "label_font": "3",  # Arial
        "label_size": "10",
        "caption_font": "4",  # Times New Roman
        "caption_size": "10",
    },
    "rsc": {
        "name": "RSC Template",
        "bond_length": 15.0,
        "line_width": 0.6,
        "bold_width": 2.0,
        "hash_spacing": 2.5,
        "margin_width": 1.6,
        "bond_spacing": 18,
        "chain_angle": 120,
        "label_font": "3",
        "label_size": "9",
        "caption_font": "4",
        "caption_size": "9",
    },
    "wiley": {
        "name": "Wiley Template",
        "bond_length": 17.0,
        "line_width": 0.6,
        "bold_width": 2.0,
        "hash_spacing": 2.5,
        "margin_width": 1.6,
        "bond_spacing": 18,
        "chain_angle": 120,
        "label_font": "3",
        "label_size": "10",
        "caption_font": "4",
        "caption_size": "10",
    },
}

PERIODIC_TABLE: Dict[int, str] = {
    1: "H",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    11: "Na",
    12: "Mg",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    19: "K",
    20: "Ca",
    26: "Fe",
    29: "Cu",
    30: "Zn",
    34: "Se",
    35: "Br",
    46: "Pd",
    53: "I",
    78: "Pt",
}


def mol_to_scaled_fragment(
    mol: Chem.Mol,
    offset_x: float = 0.0,
    offset_y: float = 0.0,
    start_node_id: int = 1,
    start_bond_id: int = 1,
    target_bond_length: float = 14.4,
) -> Tuple[List[str], List[str], int, int, Tuple[float, float, float, float]]:
  """Converts an RDKit Mol into CDXML `<n>` and `<b>` lines scaled to `target_bond_length`.

  Inverts Cartesian y-coordinates, Kekulizes aromatic bonds to `Order="2"`, and
  injects explicit `NumHydrogens` on heteroatoms to prevent ChemDraw red valence
  warning boxes.
  """
  mol_copy = Chem.Mol(mol)
  if mol_copy.GetNumConformers() == 0:
    AllChem.Compute2DCoords(mol_copy)
    regularize_nonconvex_small_rings(mol_copy)
  try:
    Chem.Kekulize(mol_copy, clearAromaticFlags=True)
  except (
      Chem.rdchem.KekulizeException,
      Chem.rdchem.MolSanitizeException,
      ValueError,
      RuntimeError,
  ):
    pass
  conf = mol_copy.GetConformer()

  bond_lengths: List[float] = []
  for b in mol_copy.GetBonds():
    p1 = conf.GetAtomPosition(b.GetBeginAtomIdx())
    p2 = conf.GetAtomPosition(b.GetEndAtomIdx())
    d = math.hypot(p2.x - p1.x, p2.y - p1.y)
    if d > 0.1:
      bond_lengths.append(d)

  if bond_lengths:
    avg_len = sum(bond_lengths) / len(bond_lengths)
    scale = target_bond_length / avg_len
  else:
    scale = 9.6

  num_atoms = mol_copy.GetNumAtoms()
  coords: List[Tuple[float, float]] = []
  for i in range(num_atoms):
    pos = conf.GetAtomPosition(i)
    coords.append((pos.x * scale, -pos.y * scale))

  if coords:
    raw_xs = [c[0] for c in coords]
    raw_ys = [c[1] for c in coords]
    center_x = (min(raw_xs) + max(raw_xs)) / 2.0
    center_y = (min(raw_ys) + max(raw_ys)) / 2.0
  else:
    center_x, center_y = 0.0, 0.0

  dx = offset_x - center_x
  dy = offset_y - center_y
  placed_coords = [(c[0] + dx, c[1] + dy) for c in coords]

  if placed_coords:
    xs = [c[0] for c in placed_coords]
    ys = [c[1] for c in placed_coords]
    bbox = (min(xs), max(xs), min(ys), max(ys))
  else:
    bbox = (offset_x, offset_x, offset_y, offset_y)

  node_elements: List[str] = []
  atom_id_map: Dict[int, int] = {}
  current_node_id = start_node_id

  for atom in mol_copy.GetAtoms():
    a_idx = atom.GetIdx()
    x, y = placed_coords[a_idx]
    elem = atom.GetAtomicNum()
    charge = atom.GetFormalCharge()
    tot_h = atom.GetTotalNumHs()

    node_id = current_node_id
    atom_id_map[a_idx] = node_id
    current_node_id += 1

    attrs = [f'id="{node_id}"', f'p="{x:.2f} {y:.2f}"']
    if elem != 6:
      attrs.append(f'Element="{elem}"')
      if tot_h > 0:
        attrs.append(f'NumHydrogens="{tot_h}"')
    if charge != 0:
      attrs.append(f'Charge="{charge}"')

    node_elements.append(f'  <n {" ".join(attrs)}/>')

  bond_elements: List[str] = []
  current_bond_id = start_bond_id
  for bond in mol_copy.GetBonds():
    b_node = atom_id_map[bond.GetBeginAtomIdx()]
    e_node = atom_id_map[bond.GetEndAtomIdx()]
    b_type = bond.GetBondType()

    order_attr = ""
    if b_type == Chem.BondType.DOUBLE:
      order_attr = ' Order="2"'
    elif b_type == Chem.BondType.TRIPLE:
      order_attr = ' Order="3"'

    display = ""
    b_dir = bond.GetBondDir()
    if b_dir == Chem.BondDir.BEGINWEDGE:
      display = ' Display="WedgeBegin"'
    elif b_dir == Chem.BondDir.BEGINDASH:
      display = ' Display="WedgedHashBegin"'

    bond_id = current_bond_id
    current_bond_id += 1
    bond_elements.append(
        f'  <b id="{bond_id}" B="{b_node}" E="{e_node}"{order_attr}{display}/>'
    )

  return node_elements, bond_elements, current_node_id, current_bond_id, bbox


def generate_styled_cdxml(
    mol: Union[Chem.Mol, Sequence[Chem.Mol]],
    style_name: str = "ACS Document 1996",
) -> str:
  """Generates publication-quality ACS 1996 / RSC / Wiley compliant CDXML from an RDKit Mol or sequence of Mols."""
  style_key = "acs1996"
  for key, val in STYLE_PRESETS.items():
    if key == style_name.lower() or val["name"].lower() == style_name.lower():
      style_key = key
      break
  cfg = STYLE_PRESETS[style_key]
  target_bond_length = float(cfg["bond_length"])
  line_width_val = float(cfg["line_width"])

  raw_mols = [mol] if isinstance(mol, Chem.Mol) else list(mol)
  frags: List[Chem.Mol] = []
  for m in raw_mols:
    if m is None:
      continue
    sub_frags = list(Chem.GetMolFrags(m, asMols=True))
    if sub_frags:
      frags.extend(sub_frags)
    else:
      frags.append(m)

  total_atoms_all = sum(f.GetNumAtoms() for f in frags)
  curr_node_id = len(frags) + 2
  curr_bond_id = curr_node_id + total_atoms_all + 2
  curr_offset_x = 150.0

  frag_blocks: List[str] = []
  all_min_x: List[float] = []
  all_max_x: List[float] = []
  all_min_y: List[float] = []
  all_max_y: List[float] = []

  for frag_idx, frag_mol in enumerate(frags, start=1):
    n_els, b_els, curr_node_id, curr_bond_id, bbox = mol_to_scaled_fragment(
        frag_mol,
        offset_x=curr_offset_x,
        offset_y=150.0,
        start_node_id=curr_node_id,
        start_bond_id=curr_bond_id,
        target_bond_length=target_bond_length,
    )
    f_min_x, f_max_x, f_min_y, f_max_y = bbox
    all_min_x.append(f_min_x)
    all_max_x.append(f_max_x)
    all_min_y.append(f_min_y)
    all_max_y.append(f_max_y)
    curr_offset_x = f_max_x + 90.0

    frag_lines = [
        f' <fragment id="{frag_idx}"'
        f' BoundingBox="{f_min_x:.2f} {f_min_y:.2f}'
        f' {f_max_x:.2f} {f_max_y:.2f}">'
    ]
    frag_lines.extend(n_els)
    frag_lines.extend(b_els)
    frag_lines.append(" </fragment>")
    frag_blocks.extend(frag_lines)

  min_x = (min(all_min_x) - 20.0) if all_min_x else 100.0
  max_x = (max(all_max_x) + 20.0) if all_max_x else 300.0
  min_y = (min(all_min_y) - 20.0) if all_min_y else 100.0
  max_y = (max(all_max_y) + 20.0) if all_max_y else 300.0

  cdxml_lines = [
      '<?xml version="1.0" encoding="UTF-8" ?>',
      '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >',
      "<CDXML",
      ' CreationProgram="ChemDraw 23.1"',
      f' Name="{cfg["name"]}"',
      f' BoundingBox="{min_x:.2f} {min_y:.2f} {max_x:.2f} {max_y:.2f}"',
      ' WindowPosition="0 0"',
      ' WindowSize="0 0"',
      ' FractionalWidths="yes"',
      ' InterpretChemically="yes"',
      ' ShowAtomQuery="yes"',
      ' ShowAtomStereo="no"',
      ' ShowAtomEnhancedStereo="yes"',
      ' ShowAtomNumber="no"',
      ' ShowResidueID="no"',
      ' ShowBondQuery="yes"',
      ' ShowBondRxn="yes"',
      ' ShowBondStereo="no"',
      ' ShowTerminalCarbonLabels="no"',
      ' ShowNonTerminalCarbonLabels="no"',
      ' HideImplicitHydrogens="no"',
      f' LabelFont="{cfg["label_font"]}"',
      f' LabelSize="{cfg["label_size"]}"',
      ' LabelFace="0"',
      f' CaptionFont="{cfg["caption_font"]}"',
      f' CaptionSize="{cfg["caption_size"]}"',
      f' HashSpacing="{cfg["hash_spacing"]}"',
      f' MarginWidth="{cfg["margin_width"]}"',
      f' LineWidth="{line_width_val:.2f}"',
      f' BoldWidth="{cfg["bold_width"]}"',
      f' BondLength="{target_bond_length:.2f}"',
      f' BondSpacing="{cfg["bond_spacing"]}"',
      f' ChainAngle="{cfg["chain_angle"]}"',
      ' LabelJustification="Auto"',
      ' CaptionJustification="Left"',
      ' color="0"',
      ' bgcolor="1"',
      ">",
      "<colortable>",
      ' <color r="1" g="1" b="1"/>',
      ' <color r="0" g="0" b="0"/>',
      "</colortable>",
      "<fonttable>",
      ' <font id="3" charset="x-mac-roman" name="Arial"/>',
      ' <font id="4" charset="x-mac-roman" name="Times New Roman"/>',
      "</fonttable>",
      (
          f'<page id="0" BoundingBox="{min_x:.2f} {min_y:.2f} {max_x:.2f}'
          f' {max_y:.2f}">'
      ),
  ]
  cdxml_lines.extend(frag_blocks)
  cdxml_lines.extend(["</page>", "</CDXML>"])
  return "\n".join(cdxml_lines)


def parse_cdxml_fallback(cdxml_path: str) -> Dict[str, Any]:
  """Pure-Python XML parser for CDXML files when rdChemDraw is unavailable."""
  tree = ET.parse(cdxml_path)
  root = tree.getroot()
  data: Dict[str, Any] = {
      "program": root.attrib.get("CreationProgram", "Unknown"),
      "style": root.attrib.get("Name", "Default"),
      "fragments": [],
      "text_boxes": [],
      "reaction_steps": [],
  }
  for t_el in root.iter("t"):
    spans = [s.text for s in t_el.findall("s") if s.text]
    text = "".join(spans).strip()
    if text:
      data["text_boxes"].append({"text": text, "pos": t_el.attrib.get("p", "")})

  for frag_idx, frag_el in enumerate(root.iter("fragment")):
    nodes = []
    bonds = []
    for n_el in frag_el.findall("n"):
      nodes.append({
          "id": n_el.attrib.get("id"),
          "element": int(n_el.attrib.get("Element", 6)),
          "charge": int(n_el.attrib.get("Charge", 0)),
          "num_hydrogens": int(n_el.attrib.get("NumHydrogens", 0)),
          "pos": n_el.attrib.get("p", "0 0"),
      })
    for b_el in frag_el.findall("b"):
      bonds.append({
          "id": b_el.attrib.get("id"),
          "begin": b_el.attrib.get("B"),
          "end": b_el.attrib.get("E"),
          "order": b_el.attrib.get("Order", "1"),
          "display": b_el.attrib.get("Display", "Solid"),
      })
    data["fragments"].append({
        "index": frag_idx,
        "node_count": len(nodes),
        "bond_count": len(bonds),
        "nodes": nodes,
        "bonds": bonds,
    })
  return data


def inspect_cdxml(filepath: str, verbose: bool = True) -> Dict[str, Any]:
  """Audits a CDXML document for ACS styling, heteroatom valences, Kekulization, and dual arrows."""
  tree = ET.parse(filepath)
  root = tree.getroot()

  fragments = list(root.iter("fragment"))
  steps = list(root.iter("step"))
  arrows = list(root.iter("arrow"))
  graphics_arrows = [
      g
      for g in root.iter("graphic")
      if g.attrib.get("ArrowType") in ("FullHead", "HalfHead", "Head")
      or g.attrib.get("GraphicType") == "Line"
  ]

  total_atoms = 0
  total_bonds = 0
  valence_warnings = 0
  unkekulized_bonds = 0
  fragments_detail: List[Dict[str, Any]] = []

  for frag_idx, frag in enumerate(fragments, start=1):
    nodes = frag.findall("n")
    bonds = frag.findall("b")
    total_atoms += len(nodes)
    total_bonds += len(bonds)
    fragments_detail.append({
        "fragment_id": frag.attrib.get("id", str(frag_idx)),
        "atoms": len(nodes),
        "bonds": len(bonds),
    })

    node_bonds: Dict[str, float] = {}
    for b in bonds:
      b_id = b.attrib.get("B", "")
      e_id = b.attrib.get("E", "")
      order_str = b.attrib.get("Order", "1")
      if order_str == "1.5":
        unkekulized_bonds += 1
      order = float(order_str) if order_str else 1.0
      node_bonds[b_id] = node_bonds.get(b_id, 0.0) + order
      node_bonds[e_id] = node_bonds.get(e_id, 0.0) + order

    for n in nodes:
      n_id = n.attrib.get("id", "")
      elem_num = int(n.attrib.get("Element", 6))
      charge = int(n.attrib.get("Charge", 0))
      num_h = int(n.attrib.get("NumHydrogens", 0))
      tot_order = node_bonds.get(n_id, 0.0) + num_h

      if elem_num == 7:
        exp_val = 3 + charge
        if abs(tot_order - exp_val) > 0.01 and abs(tot_order - 5) > 0.01:
          valence_warnings += 1
      elif elem_num == 8:
        exp_val = 2 + charge
        if abs(tot_order - exp_val) > 0.01:
          valence_warnings += 1
      elif elem_num == 16:
        if tot_order < 2.0 + charge:
          valence_warnings += 1

  bond_length_pt = float(root.attrib.get("BondLength", 14.4))
  line_width_pt = float(root.attrib.get("LineWidth", 0.6))
  chain_angle = int(float(root.attrib.get("ChainAngle", 120)))
  dual_arrows_ok = (len(arrows) == 0) or (len(graphics_arrows) >= len(arrows))
  acs_1996_compliant = (
      abs(bond_length_pt - 14.4) < 0.05
      and abs(line_width_pt - 0.6) < 0.05
      and chain_angle == 120
      and valence_warnings == 0
      and unkekulized_bonds == 0
  )
  report = {
      "file": os.path.basename(filepath),
      "style": root.attrib.get("Name", "Default"),
      "bond_length_pt": bond_length_pt,
      "line_width_pt": line_width_pt,
      "chain_angle_deg": chain_angle,
      "acs_1996_compliant": acs_1996_compliant,
      "fragment_count": len(fragments),
      "fragments": fragments_detail,
      "total_atoms": total_atoms,
      "total_bonds": total_bonds,
      "step_count": len(steps),
      "arrow_count": len(arrows),
      "graphic_arrow_count": len(graphics_arrows),
      "valence_warnings": valence_warnings,
      "unkekulized_bonds": unkekulized_bonds,
      "dual_arrows_ok": dual_arrows_ok,
      "valid": (
          valence_warnings == 0 and unkekulized_bonds == 0 and dual_arrows_ok
      ),
  }

  if verbose:
    print("=" * 60)
    print(f" ChemDraw Document Inspection: {report['file']}")
    print("=" * 60)
    print(f"  Style Template     : {report['style']}")
    print(f"  Bond Length        : {report['bond_length_pt']:.2f} pt")
    print(f"  Line Width         : {report['line_width_pt']:.2f} pt")
    print(f"  Chain Angle        : {report['chain_angle_deg']} deg")
    print(f"  ACS 1996 Compliant : {report['acs_1996_compliant']}")
    print(
        f"  Fragments          : {report['fragment_count']} ({total_atoms}"
        f" atoms, {total_bonds} bonds)"
    )
    for fd in fragments_detail:
      print(
          f"    - Fragment {fd['fragment_id']}: {fd['atoms']} atoms,"
          f" {fd['bonds']} bonds"
      )
    print(
        f"  Reaction Steps     : {len(steps)} steps, {len(arrows)} <arrow>,"
        f" {len(graphics_arrows)} <graphic> lines"
    )
    if valence_warnings == 0:
      print(
          "  [OK] Valence Attributes: All heteroatoms have satisfied valences"
          " (no red warning boxes)."
      )
    else:
      print(
          f"  [!] Valence Warning: {valence_warnings} heteroatoms missing"
          " explicit NumHydrogens."
      )
    if unkekulized_bonds > 0:
      print(
          f"  [!] Kekulization Warning: {unkekulized_bonds} bonds use"
          " Order='1.5' instead of Kekulized Order='2'."
      )
    if len(arrows) > 0:
      if dual_arrows_ok:
        print(
            "  [OK] Dual Arrows: Document contains both visual <graphic> lines"
            " and semantic <step>/<arrow> tags."
        )
      else:
        print(
            "  [!] Arrow Warning: Missing <graphic GraphicType='Line'"
            " ArrowType='FullHead'> for desktop rendering."
        )
  return report


def main() -> None:
  parser = argparse.ArgumentParser(
      description="ChemDraw CDXML/CDX Converter, Parser, and Quality Inspector."
  )
  parser.add_argument(
      "command",
      nargs="?",
      choices=["write", "inspect", "parse"],
      help="Optional subcommand (write, inspect, or parse)",
  )
  parser.add_argument(
      "-s",
      "--smiles",
      nargs="+",
      help="Input SMILES string(s) to export to CDXML",
  )
  parser.add_argument(
      "-m", "--molfile", help="Input Mol/SDF file path to export to CDXML"
  )
  parser.add_argument(
      "-i", "--input", help="Input .cdxml or .cdx file to parse or inspect"
  )
  parser.add_argument(
      "-o",
      "--output",
      help="Output file path (.cdxml, .cdx, .smi, .sdf, .json)",
  )
  parser.add_argument(
      "--style",
      choices=["acs1996", "rsc", "wiley"],
      default="acs1996",
      help="Journal style preset for CDXML export (default: acs1996)",
  )
  parser.add_argument(
      "--inspect",
      action="store_true",
      help=(
          "Audit input .cdxml file for ACS geometry, valence warnings, and dual"
          " arrow tags"
      ),
  )
  parser.add_argument(
      "-f",
      "--format",
      choices=["smiles", "sdf", "json", "rxn"],
      default="smiles",
      help=(
          "Output format when parsing an input .cdxml/.cdx file (default:"
          " smiles)"
      ),
  )
  args = parser.parse_args()

  if args.inspect or args.command == "inspect":
    target = args.input or args.output
    if not target:
      sys.stderr.write("Error: --inspect requires --input <file.cdxml>\n")
      sys.exit(1)
    inspect_cdxml(target, verbose=True)
    return

  molfile_input = args.molfile or (
      args.input if args.command == "write" else None
  )
  if args.smiles or molfile_input:
    mols: List[Chem.Mol] = []
    if args.smiles:
      for smi_arg in args.smiles:
        for token in [t.strip() for t in smi_arg.split(",") if t.strip()]:
          m = Chem.MolFromSmiles(token)
          if m is None:
            sys.stderr.write(f"Error: Could not parse input SMILES: {token}\n")
            sys.exit(1)
          mols.extend(Chem.GetMolFrags(m, asMols=True) or [m])
    elif molfile_input:
      if molfile_input.lower().endswith(".sdf"):
        for m in Chem.SDMolSupplier(molfile_input):
          if m is not None:
            mols.append(m)
      elif molfile_input.lower().endswith((".smi", ".smiles")):
        for m in Chem.SmilesMolSupplier(molfile_input, titleLine=False):
          if m is not None:
            mols.append(m)
      else:
        m = Chem.MolFromMolFile(molfile_input)
        if m is not None:
          mols.append(m)
    if not mols:
      sys.stderr.write("Error: Could not parse input molecule.\n")
      sys.exit(1)
    for m in mols:
      AllChem.Compute2DCoords(m)
      regularize_nonconvex_small_rings(m)
    style_name = STYLE_PRESETS[args.style]["name"]
    cdxml_str = generate_styled_cdxml(mols, style_name=style_name)
    if args.output:
      _write_cdxml_and_optional_cdx(cdxml_str, mols, args.output)
      print(
          f"Successfully generated styled CDXML ({style_name}): {args.output}"
      )
    else:
      print(cdxml_str)
    return

  if args.command == "write":
    parser.error(
        "write requires --smiles, --molfile, or --input <file.sdf|mol|smi>."
    )

  if args.input:
    output_str = _format_parsed_cdxml_input(args.input, args.format)
    if args.output:
      os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
      with open(args.output, "w", encoding="utf-8") as f:
        f.write(output_str + "\n")
      print(f"Successfully wrote parsed output to {args.output}")
    else:
      print(output_str)
    return

  parser.error("Must provide --smiles, --molfile, --input, or --inspect.")


def _write_cdx_binary(mols: Sequence[Chem.Mol], cdx_path: str) -> None:
  """Combines molecules and writes binary ChemDraw (.cdx) via rdChemDraw."""
  combined = mols[0]
  for extra_m in mols[1:]:
    combined = Chem.CombineMols(combined, extra_m)
  cdx_bytes = rdChemDraw.MolToChemDrawBlock(
      combined, format=rdChemDraw.CDXFormat.CDX
  )
  with open(cdx_path, "wb") as bf:
    bf.write(
        cdx_bytes
        if isinstance(cdx_bytes, bytes)
        else cdx_bytes.encode("latin1")
    )


def _write_cdxml_and_optional_cdx(
    cdxml_str: str, mols: Sequence[Chem.Mol], output_arg: str
) -> None:
  """Writes styled CDXML output and optional companion binary CDX file."""
  out_path = os.path.abspath(output_arg)
  os.makedirs(os.path.dirname(out_path), exist_ok=True)
  has_cdx_writer = _HAS_RDCHEMDRAW and hasattr(rdChemDraw, "MolToChemDrawBlock")

  if out_path.lower().endswith(".cdx"):
    # Always write a companion .cdxml alongside .cdx so XML tools can inspect it.
    cdxml_path = out_path[:-4] + ".cdxml"
    with open(cdxml_path, "w", encoding="utf-8") as f:
      f.write(cdxml_str + "\n")
    if has_cdx_writer:
      _write_cdx_binary(mols, out_path)
    else:
      with open(out_path, "w", encoding="utf-8") as f:
        f.write(cdxml_str + "\n")
    return

  with open(out_path, "w", encoding="utf-8") as f:
    f.write(cdxml_str + "\n")
  if has_cdx_writer and out_path.lower().endswith(".cdxml"):
    cdx_path = out_path[:-6] + ".cdx"
    try:
      _write_cdx_binary(mols, cdx_path)
    except (OSError, ValueError, RuntimeError) as exc:
      print(
          f"Warning: Failed to write companion binary CDX file {cdx_path}:"
          f" {exc}",
          file=sys.stderr,
      )


def _format_parsed_cdxml_input(input_path: str, output_format: str) -> str:
  """Parses a CDXML/CDX file into the requested output format (rxn/json/sdf/smi)."""
  if (
      output_format == "rxn"
      and _HAS_RDCHEMDRAW
      and hasattr(rdChemDraw, "ReactionsFromChemDrawFile")
  ):
    rxns = rdChemDraw.ReactionsFromChemDrawFile(
        input_path, sanitize=True, removeHs=True
    )
    return "\n".join(Chem.ReactionToSmiles(r) for r in rxns)

  if output_format == "json":
    meta = parse_cdxml_fallback(input_path)
    # rdChemDraw.MolsFromChemDrawFile(filename, sanitize=True, removeHs=True)
    # was added in RDKit 2024.09; sanitize=True and removeHs=True are defaults.
    if _HAS_RDCHEMDRAW and hasattr(rdChemDraw, "MolsFromChemDrawFile"):
      try:
        mols = rdChemDraw.MolsFromChemDrawFile(input_path)
        meta["rdkit_smiles"] = [
            Chem.MolToSmiles(m) for m in mols if m is not None
        ]
      except (OSError, ValueError, RuntimeError) as e:
        meta["rdkit_error"] = str(e)
    return json.dumps(meta, indent=2)

  mols: List[Chem.Mol] = []
  if _HAS_RDCHEMDRAW and hasattr(rdChemDraw, "MolsFromChemDrawFile"):
    mols = [
        m for m in rdChemDraw.MolsFromChemDrawFile(input_path) if m is not None
    ]
  if output_format == "sdf":
    sio = StringIO()
    with Chem.SDWriter(sio) as w:
      for m in mols:
        w.write(m)
    return sio.getvalue()
  return "\n".join(Chem.MolToSmiles(m) for m in mols)


if __name__ == "__main__":
  main()
