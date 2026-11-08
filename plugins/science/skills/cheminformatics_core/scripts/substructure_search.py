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

"""Performs SMARTS substructure search and MCS alignment with 2D highlighting."""

import argparse
import os
import sys
from typing import List, Tuple
from rdkit import Chem
from rdkit.Chem import rdFMCS
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

from draw_molecules import compute_clean_2d_coords


def find_substructure(
    mol: Chem.Mol,
    smarts_query: str,
    max_matches: int = 1000,
    uniquify: bool = True,
) -> List[Tuple[int, ...]]:
  """Performs substructure match of a SMARTS query against a Mol object.

  Args:
      mol: RDKit Mol object.
      smarts_query: SMARTS pattern string.
      max_matches: Maximum number of matches to return.
      uniquify: Return only unique matches.

  Returns:
      List of tuples of matched atom indices.
  """
  query = Chem.MolFromSmarts(smarts_query)
  if query is None:
    raise ValueError(f"Invalid SMARTS pattern: {smarts_query}")

  return list(
      mol.GetSubstructMatches(query, uniquify=uniquify, maxMatches=max_matches)
  )


def find_maximum_common_substructure(
    mols: List[Chem.Mol],
    atom_compare: str = "elements",
    bond_compare: str = "order",
    ring_matches_ring_only: bool = True,
    complete_rings_only: bool = False,
    timeout: int = 60,
) -> Tuple[str, int, int]:
  """Finds the Maximum Common Substructure (MCS) across a list of molecules.

  Args:
      mols: List of RDKit Mol objects.
      atom_compare: 'elements', 'any', or 'isotopes'.
      bond_compare: 'order' or 'any'.
      ring_matches_ring_only: Ring atoms match only ring atoms.
      complete_rings_only: Only match complete rings.
      timeout: Timeout in seconds.

  Returns:
      Tuple of (smarts_string, num_atoms, num_bonds).
  """
  atom_cmp = {
      "elements": rdFMCS.AtomCompare.CompareElements,
      "any": rdFMCS.AtomCompare.CompareAny,
      "isotopes": rdFMCS.AtomCompare.CompareIsotopes,
  }.get(atom_compare.lower(), rdFMCS.AtomCompare.CompareElements)

  bond_cmp = {
      "order": rdFMCS.BondCompare.CompareOrder,
      "any": rdFMCS.BondCompare.CompareAny,
  }.get(bond_compare.lower(), rdFMCS.BondCompare.CompareOrder)

  res = rdFMCS.FindMCS(
      mols,
      atomCompare=atom_cmp,
      bondCompare=bond_cmp,
      ringMatchesRingOnly=ring_matches_ring_only,
      completeRingsOnly=complete_rings_only,
      timeout=timeout,
  )
  return res.smartsString, res.numAtoms, res.numBonds


def draw_molecule_with_highlights(
    mol: Chem.Mol,
    highlight_atoms: List[int],
    output_path: str,
    img_size: Tuple[int, int] = (400, 400),
) -> None:
  """Renders 2D molecule with highlighted atom and bond indices to SVG or PNG."""
  # Find matching bonds between highlighted atoms
  highlight_bonds = []
  for bond in mol.GetBonds():
    if (
        bond.GetBeginAtomIdx() in highlight_atoms
        and bond.GetEndAtomIdx() in highlight_atoms
    ):
      highlight_bonds.append(bond.GetIdx())

  if mol.GetNumConformers() == 0:
    mol = Chem.Mol(mol)
    compute_clean_2d_coords(mol)
  prep_mol = rdMolDraw2D.PrepareMolForDrawing(mol)
  if output_path.lower().endswith(".svg"):
    drawer = rdMolDraw2D.MolDraw2DSVG(img_size[0], img_size[1])
    drawer.DrawMolecule(
        prep_mol, highlightAtoms=highlight_atoms, highlightBonds=highlight_bonds
    )
    drawer.FinishDrawing()
    with open(output_path, "w") as f:
      f.write(drawer.GetDrawingText())
  else:
    drawer = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
    drawer.DrawMolecule(
        prep_mol, highlightAtoms=highlight_atoms, highlightBonds=highlight_bonds
    )
    drawer.FinishDrawing()
    with open(output_path, "wb") as f:
      f.write(drawer.GetDrawingText())


def main():
  parser = argparse.ArgumentParser(
      description="Substructure search and MCS finder."
  )
  parser.add_argument(
      "--smiles", type=str, required=True, help="Target molecule SMILES."
  )
  parser.add_argument("--smarts", type=str, help="SMARTS query pattern.")
  parser.add_argument(
      "--mcs_with", type=str, help="Second SMILES for MCS comparison."
  )
  parser.add_argument(
      "--output_img", type=str, help="Optional image path (PNG or SVG)."
  )
  args = parser.parse_args()

  mol = Chem.MolFromSmiles(args.smiles)
  if mol is None:
    sys.stderr.write(f"Invalid SMILES: {args.smiles}\n")
    sys.exit(1)

  if args.smarts:
    matches = find_substructure(mol, args.smarts)
    print(f"Found {len(matches)} match(es) for SMARTS '{args.smarts}':")
    for i, match in enumerate(matches):
      print(f"  Match {i+1}: Atom indices {match}")

    if args.output_img and matches:
      draw_molecule_with_highlights(mol, list(matches[0]), args.output_img)
      print(f"Saved highlighted depiction to {args.output_img}")

  if args.mcs_with:
    mol2 = Chem.MolFromSmiles(args.mcs_with)
    if mol2 is None:
      sys.stderr.write(f"Invalid second SMILES: {args.mcs_with}\n")
      sys.exit(1)
    smarts, n_atoms, n_bonds = find_maximum_common_substructure([mol, mol2])
    print("MCS Results:")
    print(f"  SMARTS:    {smarts}")
    print(f"  Num Atoms: {n_atoms}")
    print(f"  Num Bonds: {n_bonds}")


if __name__ == "__main__":
  main()
