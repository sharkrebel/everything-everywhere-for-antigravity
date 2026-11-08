# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "pandas",
#     "scipy",
#     "scikit-learn",
#     "matplotlib",
#     "seaborn",
# ]
# ///

"""2D Topographic Steric Contour Map Renderer.

Generates 2D topographic contour elevation maps (isocontour plots of ligand
steric hindrance around the metal coordination sphere) for organometallic
complexes.
"""

import argparse
import math
import os
import sys
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
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

from compat import set_safe_plot_style

# Bondi Van der Waals Radii
BONDI_VDW: Dict[int, float] = {
    1: 1.20,
    6: 1.70,
    7: 1.55,
    8: 1.52,
    9: 1.47,
    14: 2.10,
    15: 1.80,
    16: 1.80,
    17: 1.75,
    35: 1.85,
    53: 1.98,
}


def generate_steric_map_grid(
    mol_3d: Chem.Mol,
    donor_idx: int,
    metal_dist: float = 2.28,
    sphere_radius: float = 3.5,
    vdw_scale: float = 1.17,
    grid_size: int = 100,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
  """Computes a 2D topographic elevation grid Z(x, y) across the coordination sphere."""
  conf = mol_3d.GetConformer()
  donor_pos = np.array(conf.GetAtomPosition(donor_idx))

  nbrs = [
      nbr.GetIdx() for nbr in mol_3d.GetAtomWithIdx(donor_idx).GetNeighbors()
  ]
  if nbrs:
    nbr_coords = [np.array(conf.GetAtomPosition(i)) for i in nbrs]
    centroid = np.mean(nbr_coords, axis=0)
    z_axis = centroid - donor_pos
    norm = np.linalg.norm(z_axis)
    z_axis = z_axis / norm if norm > 1e-4 else np.array([0.0, 0.0, 1.0])
  else:
    z_axis = np.array([0.0, 0.0, 1.0])

  metal_pos = donor_pos - z_axis * metal_dist

  # Orthogonal frame
  if abs(z_axis[0]) < 0.9:
    x_axis = np.cross(z_axis, [1.0, 0.0, 0.0])
  else:
    x_axis = np.cross(z_axis, [0.0, 1.0, 0.0])
  x_axis /= np.linalg.norm(x_axis)
  y_axis = np.cross(z_axis, x_axis)

  # 2D XY Grid spanning [-sphere_radius, sphere_radius]
  lin = np.linspace(-sphere_radius, sphere_radius, grid_size)
  grid_x, grid_y = np.meshgrid(lin, lin)
  z_elevation = np.zeros_like(grid_x)

  atoms = []
  for atom in mol_3d.GetAtoms():
    idx = atom.GetIdx()
    atomic_num = atom.GetAtomicNum()
    radius = BONDI_VDW.get(atomic_num, 1.70) * vdw_scale
    pos = np.array(conf.GetAtomPosition(idx)) - metal_pos
    # Project into local frame
    lx = np.dot(pos, x_axis)
    ly = np.dot(pos, y_axis)
    lz = np.dot(pos, z_axis)
    atoms.append((lx, ly, lz, radius))

  # For each (x, y) point on the grid, find the highest Z coordinate of the ligand VdW surface
  for i in range(grid_size):
    for j in range(grid_size):
      gx = grid_x[i, j]
      gy = grid_y[i, j]
      if gx**2 + gy**2 > sphere_radius**2:
        continue

      max_z = 0.0
      for lx, ly, lz, r in atoms:
        d_xy2 = (gx - lx) ** 2 + (gy - ly) ** 2
        if d_xy2 <= r**2:
          dz = math.sqrt(r**2 - d_xy2)
          top_z = lz + dz
          if top_z > max_z:
            max_z = top_z
      z_elevation[i, j] = min(sphere_radius, max_z)

  return grid_x, grid_y, z_elevation


def render_2d_steric_contour_map(
    mol_3d: Chem.Mol,
    donor_idx: int,
    output_png: str = "steric_map.png",
    title: str = "Topographic Steric Map",
) -> str:
  """Renders and saves a 2D topographic steric contour map."""
  set_safe_plot_style()
  grid_x, grid_y, z_elev = generate_steric_map_grid(mol_3d, donor_idx)

  fig, ax = plt.subplots(figsize=(6.5, 6), dpi=300)
  fig.patch.set_facecolor("#FFFFFF")
  ax.set_facecolor("#FCFDFE")

  # Polar circle perimeter
  circle = plt.Circle(
      (0, 0), 3.5, color="#1E293B", fill=False, linestyle="--", linewidth=1.5
  )
  ax.add_patch(circle)

  # Quadrant dividing lines
  ax.axhline(0, color="#94A3B8", linestyle=":", linewidth=1.0)
  ax.axvline(0, color="#94A3B8", linestyle=":", linewidth=1.0)

  # Mask outside sphere
  r_grid = np.sqrt(grid_x**2 + grid_y**2)
  z_masked = np.ma.masked_where(r_grid > 3.5, z_elev)

  # Filled contour plot
  levels = np.linspace(0, 3.5, 15)
  cs = ax.contourf(
      grid_x, grid_y, z_masked, levels=levels, cmap="coolwarm", extend="both"
  )
  cbar = fig.colorbar(cs, ax=ax, shrink=0.85, pad=0.04)
  cbar.set_label(
      "Topographic Elevation (Å)",
      fontsize=11,
      fontweight="bold",
      color="#1E293B",
  )

  # Contours lines
  ax.contour(
      grid_x,
      grid_y,
      z_masked,
      levels=levels[::2],
      colors="#334155",
      linewidths=0.5,
      alpha=0.6,
  )

  # Metal center marker
  ax.plot(
      0,
      0,
      marker="o",
      markersize=8,
      color="#0F172A",
      markeredgecolor="#FFFFFF",
      markeredgewidth=1.5,
      label="Metal Vertex",
  )

  # Quadrant annotations
  ax.text(
      2.6,
      2.6,
      "Q1",
      fontsize=11,
      fontweight="bold",
      color="#475569",
      ha="center",
  )
  ax.text(
      -2.6,
      2.6,
      "Q2",
      fontsize=11,
      fontweight="bold",
      color="#475569",
      ha="center",
  )
  ax.text(
      -2.6,
      -2.6,
      "Q3",
      fontsize=11,
      fontweight="bold",
      color="#475569",
      ha="center",
  )
  ax.text(
      2.6,
      -2.6,
      "Q4",
      fontsize=11,
      fontweight="bold",
      color="#475569",
      ha="center",
  )

  ax.set_xlim(-3.8, 3.8)
  ax.set_ylim(-3.8, 3.8)
  ax.set_aspect("equal")
  ax.set_xlabel("X (Å)", fontsize=11, fontweight="bold")
  ax.set_ylabel("Y (Å)", fontsize=11, fontweight="bold")
  ax.set_title(title, fontsize=13, fontweight="bold", pad=12, color="#0F172A")
  ax.legend(loc="upper right", framealpha=0.9)

  plt.tight_layout()
  os.makedirs(os.path.dirname(os.path.abspath(output_png)), exist_ok=True)
  plt.savefig(output_png, dpi=300)
  plt.close()

  return output_png


def main():
  parser = argparse.ArgumentParser(
      description="2D Topographic Steric Map Renderer."
  )
  parser.add_argument(
      "--smiles", type=str, required=True, help="Ligand SMILES string."
  )
  parser.add_argument(
      "--donor_idx",
      type=int,
      default=0,
      help="Primary donor atom index (default: 0).",
  )
  parser.add_argument(
      "--output_png",
      "--output",
      dest="output_png",
      type=str,
      default="steric_map.png",
      help="Output PNG path.",
  )
  parser.add_argument(
      "--title", type=str, default="Topographic Steric Map", help="Plot title."
  )
  args = parser.parse_args()

  mol = Chem.MolFromSmiles(args.smiles)
  mol = Chem.AddHs(mol)
  AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
  AllChem.MMFFOptimizeMolecule(mol)

  out = render_2d_steric_contour_map(
      mol,
      donor_idx=args.donor_idx,
      output_png=args.output_png,
      title=args.title,
  )
  print(f"\n[SAVED] 2D Topographic Steric Map rendered to: {out}")


if __name__ == "__main__":
  main()
