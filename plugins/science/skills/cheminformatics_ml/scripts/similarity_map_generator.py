# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "pandas",
#     "matplotlib",
# ]
# ///

"""Generates 2D atom attribution similarity maps and property heatmaps using modern RDKit rdMolDraw2D.

Renders per-atom contour attribution maps highlighting:
1. Fingerprint similarity to a reference compound (Morgan fingerprints).
2. Physicochemical property contributions (Crippen atomic LogP or Gasteiger
charges).
3. External ML model attribution weights (e.g. Random Forest, GNN, or
delta-potency weights).
"""

import argparse
import json
import os
import sys
from typing import List, Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Draw import SimilarityMaps, rdMolDraw2D

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

from draw_molecules import compute_clean_2d_coords


def create_draw2d(width: int = 400, height: int = 400, is_svg: bool = False):
  """Creates appropriate MolDraw2D canvas."""
  if is_svg:
    return rdMolDraw2D.MolDraw2DSVG(width, height)
  else:
    return rdMolDraw2D.MolDraw2DCairo(width, height)


def save_drawing(d2d, output_path: str, is_svg: bool = False):
  """Finishes drawing and writes to output file."""
  d2d.FinishDrawing()
  text = d2d.GetDrawingText()
  mode = "w" if is_svg else "wb"
  with open(output_path, mode) as f:
    f.write(text)


def generate_fingerprint_similarity_map(
    reference_smiles: str,
    target_smiles: str,
    output_path: str = "similarity_map.png",
    width: int = 400,
    height: int = 400,
) -> float:
  """Renders atomic contribution map showing which atoms drive fingerprint similarity to reference."""
  ref_mol = Chem.MolFromSmiles(reference_smiles)
  target_mol = Chem.MolFromSmiles(target_smiles)
  if ref_mol is None or target_mol is None:
    raise ValueError("Invalid reference or target SMILES.")
  compute_clean_2d_coords(ref_mol)
  compute_clean_2d_coords(target_mol, core_mol=ref_mol)

  is_svg = output_path.endswith(".svg")
  d2d = create_draw2d(width, height, is_svg)

  RDLogger.DisableLog("rdApp.warning")
  weights, max_weight = SimilarityMaps.GetSimilarityMapForFingerprint(
      ref_mol,
      target_mol,
      SimilarityMaps.GetMorganFingerprint,
      draw2d=d2d,
  )
  RDLogger.EnableLog("rdApp.warning")
  save_drawing(d2d, output_path, is_svg)
  return float(max_weight)


def generate_property_attribution_map(
    target_smiles: str,
    property_type: str = "logp",
    custom_weights: Optional[List[float]] = None,
    output_path: str = "attribution_map.png",
    width: int = 400,
    height: int = 400,
) -> None:
  """Renders atomic contribution heatmap for LogP, Gasteiger charges, or custom ML weights."""
  target_mol = Chem.MolFromSmiles(target_smiles)
  if target_mol is None:
    raise ValueError("Invalid target SMILES.")
  compute_clean_2d_coords(target_mol)

  if custom_weights is not None:
    if len(custom_weights) != target_mol.GetNumAtoms():
      raise ValueError(
          f"Custom weights count ({len(custom_weights)}) does not match atom"
          f" count ({target_mol.GetNumAtoms()})."
      )
    weights = custom_weights
  elif property_type.lower() == "logp":
    contribs = rdMolDescriptors._CalcCrippenContribs(target_mol)
    weights = [float(c[0]) for c in contribs]
  elif property_type.lower() == "charge":
    AllChem.ComputeGasteigerCharges(target_mol)
    weights = [
        float(target_mol.GetAtomWithIdx(i).GetDoubleProp("_GasteigerCharge"))
        for i in range(target_mol.GetNumAtoms())
    ]
  else:
    raise ValueError(f"Unknown property type: {property_type}")

  is_svg = output_path.endswith(".svg")
  d2d = create_draw2d(width, height, is_svg)

  SimilarityMaps.GetSimilarityMapFromWeights(target_mol, weights, draw2d=d2d)
  save_drawing(d2d, output_path, is_svg)


def main():
  parser = argparse.ArgumentParser(
      description="Generate 2D atom attribution heatmaps and similarity maps."
  )
  parser.add_argument(
      "--target-smiles",
      "--target_smiles",
      "--smiles",
      dest="target_smiles",
      required=True,
      help="Target molecule SMILES.",
  )
  parser.add_argument(
      "--ref-smiles",
      "--ref_smiles",
      dest="ref_smiles",
      default=None,
      help="Reference molecule SMILES for fingerprint comparison.",
  )
  parser.add_argument(
      "--property",
      default="logp",
      choices=["logp", "charge", "custom"],
      help="Property type for atomic heatmap.",
  )
  parser.add_argument(
      "--weights",
      default=None,
      help=(
          "Comma-separated or JSON list of custom atom weights (for --property"
          " custom)."
      ),
  )
  parser.add_argument(
      "--output",
      "--output-png",
      "--output_png",
      dest="output",
      default="attribution_map.png",
      help="Output PNG or SVG file path.",
  )
  parser.add_argument(
      "--width", type=int, default=400, help="Image width (default: 400)."
  )
  parser.add_argument(
      "--height", type=int, default=400, help="Image height (default: 400)."
  )
  args = parser.parse_args()

  if args.ref_smiles:
    print(
        "Generating fingerprint similarity map to reference"
        f" {args.ref_smiles}..."
    )
    max_wt = generate_fingerprint_similarity_map(
        args.ref_smiles,
        args.target_smiles,
        output_path=args.output,
        width=args.width,
        height=args.height,
    )
    print(f"Similarity map saved to {args.output} (Max weight: {max_wt:.3f})")
  else:
    custom_w = None
    if args.weights:
      try:
        custom_w = json.loads(args.weights)
      except Exception:
        custom_w = [float(x.strip()) for x in args.weights.split(",")]
    print(
        f"Generating property heatmap ({args.property}) for"
        f" {args.target_smiles}..."
    )
    generate_property_attribution_map(
        args.target_smiles,
        property_type=args.property,
        custom_weights=custom_w,
        output_path=args.output,
        width=args.width,
        height=args.height,
    )
    print(f"Attribution map saved to {args.output}")


if __name__ == "__main__":
  main()
