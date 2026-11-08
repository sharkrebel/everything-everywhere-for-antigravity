# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy",
#     "pandas",
#     "rdkit",
# ]
# ///
"""Reaction Scheme to Styled ChemDraw CDXML Converter.

Generates publication-quality single- and multi-step reaction schemes in
ChemDraw
CDXML format conforming to ACS Document 1996 standards:
1. ACS Document 1996 geometry (14.4 pt bond length, 0.6 pt line width, 2.0 pt
   bold width, Arial 10 pt atom labels, Times New Roman 10 pt captions).
2. Dual arrow tags for 100% Desktop and Web compatibility:
   - `<graphic GraphicType="Line" ArrowType="FullHead" HeadSize="1000">` with
     point 1 at `arrow_head_x` (right) and point 2 at `arrow_tail_x` (left) so
     ChemDraw Desktop renders arrowheads pointing forward toward products.
   - `<step ReactionType="Normal"><arrow ArrowheadHead="Full" .../><objects>...`
     for ChemDraw JS and `rdChemDraw.ReactionsFromChemDrawFile`.
3. Explicit `NumHydrogens` on all heteroatoms with implicit hydrogens and
   explicit Kekulization (`Order="2"`) to prevent ChemDraw red valence warning
   boxes.
4. Collision-free layout:
   - Dynamic arrow length `max(115.0 pt, text_width + 30.0 pt)` with >= 40.0 pt
     clearance from reactant and product bounding boxes.
   - Reagent/condition strings wrapped to <= 35 characters per line using native
     CDXML newlines (`&#10;`) with `LineHeight="auto"` and
     `Justification="Center"`.
   - Dynamic compound bottom labels at `Y = Y_bbox_max + 18.0 pt`.
   - Multi-row scheme wrapping separated by `>= 280.0 pt`
     (`ROW_Y_SEPARATION_PT`).
"""

import argparse
import json
import math
import os
import sys
import textwrap
from typing import Any, Dict, List, Tuple
from xml.sax.saxutils import escape

from rdkit import Chem
from rdkit.Chem import AllChem

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

BOND_LENGTH_PT = 14.4
LINE_WIDTH_PT = 0.6
BOLD_WIDTH_PT = 2.0
HASH_SPACING_PT = 2.5
MARGIN_WIDTH_PT = 1.6
CHAIN_ANGLE_DEG = 120
BOND_SPACING_PCT = 18

ARROW_LENGTH_PT = 115.0
ARROW_HEAD_SIZE = 1000
MOL_CLEARANCE_PT = 40.0
PLUS_CLEARANCE_PT = 25.0
LABEL_Y_OFFSET_PT = 18.0
REAGENT_WRAP_WIDTH = 35
ROW_Y_SEPARATION_PT = 280.0


def wrap_text_lines(text: str, width: int = REAGENT_WRAP_WIDTH) -> List[str]:
  """Wraps text cleanly into lines with at most `width` characters."""
  if not text:
    return []
  lines: List[str] = []
  for rline in str(text).split("\n"):
    rline = rline.strip()
    if not rline:
      continue
    lines.extend(textwrap.wrap(rline, width=width))
  return lines


def mol_to_scaled_fragment(
    mol: Chem.Mol,
    offset_x: float = 0.0,
    offset_y: float = 0.0,
    start_node_id: int = 1,
    start_bond_id: int = 1,
) -> Tuple[List[str], List[str], int, int, Tuple[float, float, float, float]]:
  """Converts an RDKit Mol into CDXML nodes and bonds scaled to 14.4 pt bond length.

  Places the fragment so its left bounding edge (`min_x`) starts at `offset_x`
  and its vertical center aligns with `offset_y`.
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
    scale = BOND_LENGTH_PT / avg_len
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
    min_raw_x = min(raw_xs)
    center_y = (min(raw_ys) + max(raw_ys)) / 2.0
  else:
    min_raw_x, center_y = 0.0, 0.0

  # Align left edge to offset_x so clearance from previous element is exact (>= 40 pt)
  dx = offset_x - min_raw_x
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

    node_elements.append(f'    <n {" ".join(attrs)}/>')

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
        f'    <b id="{bond_id}" B="{b_node}"'
        f' E="{e_node}"{order_attr}{display}/>'
    )

  return node_elements, bond_elements, current_node_id, current_bond_id, bbox


def normalize_steps_input(raw_data: Any) -> List[Dict[str, Any]]:
  """Normalizes either a list of synthesis step dicts or a render_mechanism.py spec dict."""
  if isinstance(raw_data, list):
    return raw_data
  if isinstance(raw_data, dict):
    if "steps" in raw_data and isinstance(raw_data["steps"], list):
      normalized: List[Dict[str, Any]] = []
      for s in raw_data["steps"]:
        reacts = s.get("reactants") or s.get("species") or []
        prods = s.get("products") or []
        r_smi = [m["smiles"] if isinstance(m, dict) else str(m) for m in reacts]
        p_smi = [m["smiles"] if isinstance(m, dict) else str(m) for m in prods]
        r_lbl = [
            m.get("name", "") if isinstance(m, dict) else "" for m in reacts
        ]
        p_lbl = [
            m.get("name", "") if isinstance(m, dict) else "" for m in prods
        ]
        normalized.append({
            "reactants": r_smi,
            "products": p_smi,
            "reagents": (
                s.get("reagents")
                or s.get("conditions")
                or s.get("step_title")
                or ""
            ),
            "reactant_labels": s.get("reactant_labels") or r_lbl,
            "product_labels": s.get("product_labels") or p_lbl,
        })
      return normalized
    if "reactants" in raw_data or "species" in raw_data:
      reacts = raw_data.get("reactants") or raw_data.get("species") or []
      prods = raw_data.get("products") or []
      r_smi = [m["smiles"] if isinstance(m, dict) else str(m) for m in reacts]
      p_smi = [m["smiles"] if isinstance(m, dict) else str(m) for m in prods]
      r_lbl = [m.get("name", "") if isinstance(m, dict) else "" for m in reacts]
      p_lbl = [m.get("name", "") if isinstance(m, dict) else "" for m in prods]
      return [{
          "reactants": r_smi,
          "products": p_smi,
          "reagents": (
              raw_data.get("reagents")
              or raw_data.get("conditions")
              or raw_data.get("title")
              or ""
          ),
          "reactant_labels": raw_data.get("reactant_labels") or r_lbl,
          "product_labels": raw_data.get("product_labels") or p_lbl,
      }]
  raise ValueError(
      "Unsupported reaction steps structure; expected list of step dicts or"
      " mechanism spec dict."
  )


def generate_reaction_cdxml(
    steps_data: Any,
    title: str = "Synthesis Scheme",
    steps_per_row: int = 2,
) -> str:
  """Generates publication-quality ACS 1996 CDXML for a single- or multi-step reaction scheme."""
  steps_list = normalize_steps_input(steps_data)
  id_counter = 1
  page_elements: List[str] = []
  step_elements: List[str] = []

  left_margin_x = 72.0
  current_x = left_margin_x
  current_y = 120.0
  max_seen_x = 720.0

  for step_idx, step in enumerate(steps_list):
    if step_idx > 0 and steps_per_row > 0 and (step_idx % steps_per_row == 0):
      current_x = left_margin_x
      current_y += ROW_Y_SEPARATION_PT

    reactants_smiles = step.get("reactants", [])
    products_smiles = step.get("products", [])
    reagents_text = step.get("reagents", "")
    reactant_labels = step.get("reactant_labels", [])
    product_labels = step.get("product_labels", [])

    reactant_frag_ids: List[int] = []
    product_frag_ids: List[int] = []

    # 1. Place Reactants
    for r_idx, smi in enumerate(reactants_smiles):
      mol = Chem.MolFromSmiles(smi)
      if mol is None:
        continue
      if r_idx > 0:
        plus_x = current_x + PLUS_CLEARANCE_PT
        page_elements.append(
            f'  <t id="{id_counter}" p="{plus_x:.2f} {current_y:.2f}"'
            ' CaptionJustification="Center"><s font="3" size="12"'
            ' face="Bold">+</s></t>'
        )
        id_counter += 1
        current_x += PLUS_CLEARANCE_PT * 2.0

      frag_id = id_counter
      id_counter += 1
      reactant_frag_ids.append(frag_id)

      n_els, b_els, id_counter, b_id_cnt, bbox = mol_to_scaled_fragment(
          mol,
          offset_x=current_x,
          offset_y=current_y,
          start_node_id=id_counter,
          start_bond_id=id_counter + mol.GetNumAtoms() + 10,
      )
      id_counter = b_id_cnt

      page_elements.append(f'  <fragment id="{frag_id}">')
      page_elements.extend(n_els)
      page_elements.extend(b_els)
      page_elements.append("  </fragment>")

      if r_idx < len(reactant_labels) and reactant_labels[r_idx]:
        lbl = escape(str(reactant_labels[r_idx]))
        lbl_x = (bbox[0] + bbox[1]) / 2.0
        lbl_y = bbox[3] + LABEL_Y_OFFSET_PT
        page_elements.append(
            f'  <t id="{id_counter}" p="{lbl_x:.2f} {lbl_y:.2f}"'
            ' CaptionJustification="Center"><s font="4" size="10"'
            f' face="Bold">{lbl}</s></t>'
        )
        id_counter += 1

      current_x = bbox[1] + MOL_CLEARANCE_PT

    # 2. Place Arrow and Conditions
    lines = wrap_text_lines(reagents_text, width=REAGENT_WRAP_WIDTH)
    if lines:
      approx_text_width = max(len(l) for l in lines) * 5.6
      active_arrow_len = max(ARROW_LENGTH_PT, approx_text_width + 30.0)
    else:
      active_arrow_len = ARROW_LENGTH_PT

    arrow_tail_x = current_x
    arrow_head_x = arrow_tail_x + active_arrow_len
    arrow_y = current_y

    arrow_graphic_id = id_counter
    id_counter += 1
    p1_id = id_counter
    id_counter += 1
    p2_id = id_counter
    id_counter += 1
    arrow_step_id = id_counter
    id_counter += 1

    # DUAL ARROW TAG 1: <graphic GraphicType="Line" ArrowType="FullHead">
    # Point 1 is the arrow HEAD (right) and Point 2 is the TAIL (left) so ChemDraw
    # Desktop renders the arrowhead pointing forward toward the product.
    page_elements.append(
        f'  <graphic id="{arrow_graphic_id}" SupersededBy="{arrow_step_id}"'
        f' GraphicType="Line" ArrowType="FullHead" HeadSize="{ARROW_HEAD_SIZE}"'
        f' BoundingBox="{arrow_head_x:.2f} {arrow_y:.2f}'
        f' {arrow_tail_x:.2f} {arrow_y:.2f}">\n    <point id="{p1_id}"'
        f' p="{arrow_head_x:.2f} {arrow_y:.2f}"/>\n    <point id="{p2_id}"'
        f' p="{arrow_tail_x:.2f} {arrow_y:.2f}"/>\n  </graphic>'
    )

    if lines:
      text_center_x = (arrow_tail_x + arrow_head_x) / 2.0
      text_start_y = arrow_y - 12.0 - ((len(lines) - 1) * 11.0)
      text_content = "&#10;".join(escape(l) for l in lines)
      page_elements.append(
          f'  <t id="{id_counter}" p="{text_center_x:.2f} {text_start_y:.2f}"'
          ' CaptionJustification="Center" Justification="Center"'
          f' LineHeight="auto"><s font="3" size="9.5">{text_content}</s></t>'
      )
      id_counter += 1

    current_x = arrow_head_x + MOL_CLEARANCE_PT

    # 3. Place Products
    for p_idx, smi in enumerate(products_smiles):
      mol = Chem.MolFromSmiles(smi)
      if mol is None:
        continue
      if p_idx > 0:
        plus_x = current_x + PLUS_CLEARANCE_PT
        page_elements.append(
            f'  <t id="{id_counter}" p="{plus_x:.2f} {current_y:.2f}"'
            ' CaptionJustification="Center"><s font="3" size="12"'
            ' face="Bold">+</s></t>'
        )
        id_counter += 1
        current_x += PLUS_CLEARANCE_PT * 2.0

      frag_id = id_counter
      id_counter += 1
      product_frag_ids.append(frag_id)

      n_els, b_els, id_counter, b_id_cnt, bbox = mol_to_scaled_fragment(
          mol,
          offset_x=current_x,
          offset_y=current_y,
          start_node_id=id_counter,
          start_bond_id=id_counter + mol.GetNumAtoms() + 10,
      )
      id_counter = b_id_cnt

      page_elements.append(f'  <fragment id="{frag_id}">')
      page_elements.extend(n_els)
      page_elements.extend(b_els)
      page_elements.append("  </fragment>")

      if p_idx < len(product_labels) and product_labels[p_idx]:
        lbl = escape(str(product_labels[p_idx]))
        lbl_x = (bbox[0] + bbox[1]) / 2.0
        lbl_y = bbox[3] + LABEL_Y_OFFSET_PT
        page_elements.append(
            f'  <t id="{id_counter}" p="{lbl_x:.2f} {lbl_y:.2f}"'
            ' CaptionJustification="Center"><s font="4" size="10"'
            f' face="Bold">{lbl}</s></t>'
        )
        id_counter += 1

      current_x = bbox[1] + MOL_CLEARANCE_PT

    max_seen_x = max(max_seen_x, current_x + 72.0)

    # DUAL ARROW TAG 2: <step> and <arrow> reaction metadata anchor
    objects_xml = ""
    for r_fid in reactant_frag_ids:
      objects_xml += f'\n      <object id="{r_fid}"/>'
    for p_fid in product_frag_ids:
      objects_xml += f'\n      <object id="{p_fid}"/>'

    step_elements.append(
        f'  <step id="{id_counter}" ReactionType="Normal">\n    <arrow'
        f' id="{arrow_step_id}" ArrowheadHead="Full" ArrowheadType="Solid"'
        f' HeadSize="{ARROW_HEAD_SIZE}" Head3D="{arrow_head_x:.2f}'
        f' {arrow_y:.2f} 0" Tail3D="{arrow_tail_x:.2f} {arrow_y:.2f} 0"'
        f' Center3D="{(arrow_tail_x + arrow_head_x)/2.0:.2f} {arrow_y:.2f}'
        f' 0"/>\n    <objects>{objects_xml}\n    </objects>\n  </step>'
    )
    id_counter += 1

  doc_max_x = max(max_seen_x, 720.0)
  doc_max_y = max(current_y + 180.0, 540.0)

  cdxml_lines = [
      '<?xml version="1.0" encoding="UTF-8" ?>',
      '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >',
      "<CDXML",
      ' CreationProgram="ChemDraw 23.1"',
      ' Name="ACS Document 1996"',
      f' BoundingBox="0 0 {doc_max_x:.2f} {doc_max_y:.2f}"',
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
      ' LabelFont="3"',
      ' LabelSize="10"',
      ' LabelFace="0"',
      ' CaptionFont="4"',
      ' CaptionSize="10"',
      f' HashSpacing="{HASH_SPACING_PT:.1f}"',
      f' MarginWidth="{MARGIN_WIDTH_PT:.1f}"',
      f' LineWidth="{LINE_WIDTH_PT:.1f}"',
      f' BoldWidth="{BOLD_WIDTH_PT:.1f}"',
      f' BondLength="{BOND_LENGTH_PT:.1f}"',
      f' BondSpacing="{BOND_SPACING_PCT}"',
      f' ChainAngle="{CHAIN_ANGLE_DEG}"',
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
      f'<page id="0" BoundingBox="0 0 {doc_max_x:.2f} {doc_max_y:.2f}">',
  ]
  cdxml_lines.extend(page_elements)
  cdxml_lines.extend(step_elements)
  cdxml_lines.extend(["</page>", "</CDXML>"])
  return "\n".join(cdxml_lines)


def main() -> None:
  parser = argparse.ArgumentParser(
      description=(
          "Generate publication-grade ACS 1996 ChemDraw (.cdxml) reaction"
          " schemes with collision avoidance and dual desktop/web arrows."
      )
  )
  parser.add_argument(
      "-j",
      "--json",
      help=(
          "Input JSON file defining reaction steps (or render_mechanism.py JSON"
          " spec)"
      ),
  )
  parser.add_argument(
      "--rxn-smiles",
      help=(
          "Single reaction SMILES string (e.g. 'react1.react2>reagents>prod1')"
      ),
  )
  parser.add_argument(
      "-o", "--output", required=True, help="Output .cdxml file path"
  )
  parser.add_argument(
      "--title", default="Synthesis Scheme", help="Scheme title"
  )
  parser.add_argument(
      "--steps-per-row",
      type=int,
      default=2,
      help=(
          "Maximum steps per horizontal row before wrapping by 280 pt"
          " (default: 2)"
      ),
  )
  args = parser.parse_args()

  if args.json:
    with open(args.json, "r", encoding="utf-8") as f:
      steps_data = json.load(f)
  elif args.rxn_smiles:
    parts = args.rxn_smiles.split(">")
    if len(parts) == 3:
      r_part, a_part, p_part = parts
    elif len(parts) == 2:
      r_part, a_part, p_part = parts[0], "", parts[1]
    else:
      sys.stderr.write(f"Error: Invalid reaction SMILES '{args.rxn_smiles}'\n")
      sys.exit(1)
    steps_data = [{
        "reactants": [s for s in r_part.split(".") if s],
        "products": [s for s in p_part.split(".") if s],
        "reagents": a_part,
    }]
  else:
    parser.error("Either --rxn or --json is required.")

  cdxml_content = generate_reaction_cdxml(
      steps_data, title=args.title, steps_per_row=args.steps_per_row
  )
  os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
  with open(args.output, "w", encoding="utf-8") as f:
    f.write(cdxml_content + "\n")

  print(f"Successfully generated reaction scheme CDXML: {args.output}")


if __name__ == "__main__":
  main()
