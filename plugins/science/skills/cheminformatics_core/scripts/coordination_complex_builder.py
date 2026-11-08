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

"""3D Coordination Complex Assembler & Force Field Minimization Engine.

Assembles 3D transition metal coordination complexes from central metal symbols
and ligand SMILES using ideal coordination polyhedra templates:
1. Linear (CN=2) - e.g. Au(I), Cu(I)
2. Square Planar (CN=4) - e.g. Pd(II), Pt(II), Rh(I), Ir(I)
3. Tetrahedral (CN=4) - e.g. Ni(0), Zn(II), Cu(I)
4. Trigonal Bipyramidal (CN=5) - e.g. Fe(0), Ru(II)
5. Octahedral (CN=6) - e.g. Ru(II), Ir(III), Rh(III), Fe(II/III)

Applies Universal Force Field (UFF) energy minimization with harmonic
coordination core restraints to relax ligand substituents without distorting the
metal geometry.
"""

import argparse
import math
import os
from typing import Any, Dict, List, Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Ideal Coordination Vectors (Unit vectors for coordination polyhedra slots)
COORDINATION_TEMPLATES: Dict[str, Dict[str, Any]] = {
    "linear": {
        "coordination_number": 2,
        "description": "Linear (180°)",
        "vectors": [
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ],
    },
    "square_planar": {
        "coordination_number": 4,
        "description": "Square Planar (90° in XY plane)",
        "vectors": [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ],
    },
    "tetrahedral": {
        "coordination_number": 4,
        "description": "Tetrahedral (109.5°)",
        "vectors": [
            [1.0 / math.sqrt(3), 1.0 / math.sqrt(3), 1.0 / math.sqrt(3)],
            [1.0 / math.sqrt(3), -1.0 / math.sqrt(3), -1.0 / math.sqrt(3)],
            [-1.0 / math.sqrt(3), 1.0 / math.sqrt(3), -1.0 / math.sqrt(3)],
            [-1.0 / math.sqrt(3), -1.0 / math.sqrt(3), 1.0 / math.sqrt(3)],
        ],
    },
    "trigonal_bipyramidal": {
        "coordination_number": 5,
        "description": (
            "Trigonal Bipyramidal (3 equatorial at 120°, 2 axial at 180°)"
        ),
        "vectors": [
            [1.0, 0.0, 0.0],
            [-0.5, math.sqrt(3) / 2.0, 0.0],
            [-0.5, -math.sqrt(3) / 2.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ],
    },
    "octahedral": {
        "coordination_number": 6,
        "description": "Octahedral (90° / 180°)",
        "vectors": [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ],
    },
}

# Default typical Metal-Donor bond lengths in Angstroms
DEFAULT_BOND_LENGTHS: Dict[str, float] = {
    "P": 2.30,
    "N": 2.05,
    "C": 2.00,
    "O": 2.05,
    "S": 2.35,
    "Cl": 2.35,
    "Br": 2.45,
    "I": 2.60,
    "F": 1.95,
    "H": 1.60,
}


def build_coordination_complex(
    metal_symbol: str,
    geometry: str,
    ligand_smiles_list: List[str],
    output_sdf: Optional[str] = None,
    minimize: bool = True,
) -> Chem.Mol:
  """Assembles a 3D coordination complex from metal symbol and ligand SMILES."""
  geom_key = geometry.lower().replace("-", "_").replace(" ", "_")
  if geom_key not in COORDINATION_TEMPLATES:
    raise ValueError(
        f"Unknown geometry '{geometry}'. Supported:"
        f" {list(COORDINATION_TEMPLATES.keys())}"
    )

  template = COORDINATION_TEMPLATES[geom_key]
  cn = template["coordination_number"]
  slot_vectors = [np.array(v) for v in template["vectors"]]

  if len(ligand_smiles_list) != cn:
    raise ValueError(
        f"Geometry '{geometry}' requires exactly {cn} ligands, but"
        f" {len(ligand_smiles_list)} were provided."
    )

  # 1. Create master RWMol and add central metal atom at origin
  master_mol = Chem.RWMol()
  metal_atom = Chem.Atom(metal_symbol)
  metal_idx = master_mol.AddAtom(metal_atom)

  conformer = Chem.Conformer(1)
  conformer.SetAtomPosition(metal_idx, (0.0, 0.0, 0.0))

  donor_atom_indices = []

  # 2. Process each ligand slot
  for slot_idx, (lig_smi, slot_vec) in enumerate(
      zip(ligand_smiles_list, slot_vectors)
  ):
    # Handle monoatomic anions/ligands (Cl, Br, I, F, H) directly
    if lig_smi in (
        "Cl",
        "[Cl-]",
        "Br",
        "[Br-]",
        "I",
        "[I-]",
        "F",
        "[F-]",
        "[H-]",
        "H",
    ):
      sym = lig_smi.strip("[]-+")
      if not sym:
        sym = "Cl"
      lig_atom = Chem.Atom(sym)
      atom_idx = master_mol.AddAtom(lig_atom)
      bond_len = DEFAULT_BOND_LENGTHS.get(sym, 2.35)
      pos = slot_vec * bond_len
      conformer.SetAtomPosition(atom_idx, tuple(pos))
      master_mol.AddBond(atom_idx, metal_idx, Chem.BondType.SINGLE)
      donor_atom_indices.append(atom_idx)
      continue

    # Multi-atom ligand: generate 3D conformer
    lig_mol = Chem.MolFromSmiles(lig_smi)
    if lig_mol is None:
      raise ValueError(f"Could not parse ligand SMILES '{lig_smi}'")
    lig_mol = Chem.AddHs(lig_mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42 + slot_idx
    params.numThreads = 0
    AllChem.EmbedMolecule(lig_mol, params)
    AllChem.MMFFOptimizeMolecule(lig_mol)

    # Identify donor atom (P > N > S > C > O)
    donor_local_idx = 0
    for atomic_num in (15, 7, 16, 6, 8):
      found = False
      for a in lig_mol.GetAtoms():
        if a.GetAtomicNum() == atomic_num:
          donor_local_idx = a.GetIdx()
          found = True
          break
      if found:
        break

    donor_atom_obj = lig_mol.GetAtomWithIdx(donor_local_idx)
    donor_sym = donor_atom_obj.GetSymbol()
    bond_len = DEFAULT_BOND_LENGTHS.get(donor_sym, 2.20)

    # Align ligand coordinate frame: Donor placed along slot_vec at distance bond_len
    lig_conf = lig_mol.GetConformer()
    d_pos = np.array(lig_conf.GetAtomPosition(donor_local_idx))

    # Vector from donor to ligand center of mass
    nbrs = [
        np.array(lig_conf.GetAtomPosition(n.GetIdx()))
        for n in donor_atom_obj.GetNeighbors()
    ]
    if nbrs:
      sub_vec = np.mean(nbrs, axis=0) - d_pos
      norm = np.linalg.norm(sub_vec)
      sub_vec = sub_vec / norm if norm > 1e-4 else np.array([0.0, 0.0, 1.0])
    else:
      sub_vec = np.array([0.0, 0.0, 1.0])

    # Rotation matrix to align sub_vec with slot_vec (pointing outward from metal)
    target_dir = slot_vec / np.linalg.norm(slot_vec)
    v = np.cross(sub_vec, target_dir)
    c = np.dot(sub_vec, target_dir)
    s = np.linalg.norm(v)

    if s > 1e-5:
      vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
      rot_matrix = np.eye(3) + vx + np.dot(vx, vx) * ((1 - c) / (s**2))
    elif c < -0.99:
      rot_matrix = -np.eye(3)
    else:
      rot_matrix = np.eye(3)

    target_donor_pos = target_dir * bond_len

    # Add ligand atoms to master molecule
    atom_map = {}
    for l_atom in lig_mol.GetAtoms():
      old_idx = l_atom.GetIdx()
      new_idx = master_mol.AddAtom(Chem.Atom(l_atom.GetSymbol()))
      atom_map[old_idx] = new_idx

      # Transform coordinate
      orig_pos = np.array(lig_conf.GetAtomPosition(old_idx)) - d_pos
      rotated_pos = np.dot(rot_matrix, orig_pos) + target_donor_pos
      conformer.SetAtomPosition(new_idx, tuple(rotated_pos))

    # Add ligand internal bonds
    for b in lig_mol.GetBonds():
      master_mol.AddBond(
          atom_map[b.GetBeginAtomIdx()],
          atom_map[b.GetEndAtomIdx()],
          b.GetBondType(),
      )

    # Add coordination dative bond from donor to metal
    master_donor_idx = atom_map[donor_local_idx]
    master_mol.AddBond(master_donor_idx, metal_idx, Chem.BondType.DATIVE)
    donor_atom_indices.append(master_donor_idx)

  final_mol = master_mol.GetMol()
  final_mol.AddConformer(conformer, assignId=True)
  try:
    final_mol.UpdatePropertyCache(strict=False)
  except Exception:
    pass

  # 3. UFF Energy Minimization with Restrained Core
  if minimize:
    try:
      uff = AllChem.UFFGetMoleculeForceField(final_mol)
      if uff is not None:
        # Fix central metal position
        uff.AddFixedPoint(metal_idx)
        # Restrain donor atom distances
        for d_idx in donor_atom_indices:
          uff.AddDistanceConstraint(metal_idx, d_idx, 1.8, 2.6, 1000.0)
        uff.Minimize(maxIts=500)
    except Exception:
      pass

  if output_sdf:
    os.makedirs(os.path.dirname(os.path.abspath(output_sdf)), exist_ok=True)
    with Chem.SDWriter(output_sdf) as writer:
      writer.write(final_mol)

  return final_mol


def main():
  parser = argparse.ArgumentParser(
      description="3D Coordination Complex Assembler & UFF Minimizer."
  )
  parser.add_argument(
      "--metal",
      type=str,
      required=True,
      help="Central metal symbol (e.g. Pd, Pt, Ru, Ir, Fe, Ni).",
  )
  parser.add_argument(
      "--geometry",
      type=str,
      required=True,
      choices=[
          "linear",
          "square_planar",
          "tetrahedral",
          "trigonal_bipyramidal",
          "octahedral",
      ],
      help="Coordination geometry.",
  )
  parser.add_argument(
      "--ligands",
      type=str,
      nargs="+",
      required=True,
      help="List of ligand SMILES strings.",
  )
  parser.add_argument(
      "--output_sdf",
      "--output",
      dest="output_sdf",
      type=str,
      default="complex.sdf",
      help="Output 3D SDF filepath.",
  )
  parser.add_argument(
      "--no_minimize", action="store_true", help="Skip UFF minimization."
  )
  args = parser.parse_args()

  mol = build_coordination_complex(
      metal_symbol=args.metal,
      geometry=args.geometry,
      ligand_smiles_list=args.ligands,
      output_sdf=args.output_sdf,
      minimize=not args.no_minimize,
  )

  print("\n=== Coordination Complex Assembled ===")
  print(f"Formula:         {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
  print(f"Total Atoms:     {mol.GetNumAtoms()}")
  print(
      f"Geometry:        {COORDINATION_TEMPLATES[args.geometry]['description']}"
  )
  print(f"Saved 3D SDF to: {args.output_sdf}")


if __name__ == "__main__":
  main()
