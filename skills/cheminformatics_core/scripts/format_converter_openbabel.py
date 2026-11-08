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
# Note: `openbabel-wheel` is an optional runtime accelerator (`uv run --with
# openbabel-wheel ...`) omitted from mandatory PEP 723 dependencies so `uv run`
# does not fail package resolution on Python versions/platforms lacking a
# prebuilt wheel before reaching the native RDKit fallback below.

"""Chemical Format Converter supporting Open Babel / Pybel with Graceful Native RDKit Fallback."""

import argparse
import os
import sys
from typing import Optional
from rdkit import Chem
from rdkit.Chem import AllChem

try:
  from openbabel import pybel

  HAS_OPENBABEL = True
except ImportError:
  pybel = None
  HAS_OPENBABEL = False


def convert_format_openbabel(
    input_path: str,
    output_path: str,
    in_fmt: Optional[str] = None,
    out_fmt: Optional[str] = None,
    add_hs: bool = False,
    gen_3d: bool = False,
    partial_charges: Optional[str] = None,
) -> bool:
  """Converts chemical formats using Open Babel (supports 100+ formats)."""
  if not HAS_OPENBABEL or pybel is None:
    return False
  in_format = in_fmt or os.path.splitext(input_path)[1].lstrip(".")
  out_format = out_fmt or os.path.splitext(output_path)[1].lstrip(".")

  mols = list(pybel.readfile(in_format, input_path))
  if not mols:
    print(
        f"Error: No molecules found in {input_path} as format '{in_format}'",
        file=sys.stderr,
    )
    return False

  out_file = pybel.Outputfile(out_format, output_path, overwrite=True)
  count = 0
  for mol in mols:
    if add_hs:
      mol.addh()
    if gen_3d:
      mol.make3D()
    if partial_charges:
      mol.OBMol.SetPartialChargesPerceived()
    out_file.write(mol)
    count += 1
  out_file.close()

  print(
      f"[Open Babel Engine] Converted {count} molecule(s) from '{in_format}' to"
      f" '{out_format}' -> {output_path}"
  )
  return True


def convert_format_rdkit_fallback(
    input_path: str,
    output_path: str,
    in_fmt: Optional[str] = None,
    out_fmt: Optional[str] = None,
    add_hs: bool = False,
    gen_3d: bool = False,
) -> bool:
  """Converts standard chemical formats (SMILES, SDF, MOL, PDB) using pure RDKit."""
  in_format = (in_fmt or os.path.splitext(input_path)[1].lstrip(".")).lower()
  out_format = (out_fmt or os.path.splitext(output_path)[1].lstrip(".")).lower()

  # Supported native RDKit formats
  native_formats = {"smi", "smiles", "sdf", "mol", "pdb"}
  if in_format not in native_formats or out_format not in native_formats:
    print(
        f"Error: Format '{in_format}' -> '{out_format}' requires Open Babel.",
        file=sys.stderr,
    )
    return False

  # Read molecules
  mols = []
  if in_format in ("smi", "smiles"):
    suppl = Chem.SmilesMolSupplier(input_path, titleLine=False)
    mols = [m for m in suppl if m is not None]
  elif in_format == "sdf":
    suppl = Chem.SDMolSupplier(input_path)
    mols = [m for m in suppl if m is not None]
  elif in_format == "mol":
    m = Chem.MolFromMolFile(input_path)
    if m:
      mols.append(m)
  elif in_format == "pdb":
    m = Chem.MolFromPDBFile(input_path)
    if m:
      mols.append(m)

  if not mols:
    print(f"Error: No valid molecules read from {input_path}", file=sys.stderr)
    return False

  # Optional 3D / Hydrogens
  processed_mols = []
  for m in mols:
    if add_hs:
      m = Chem.AddHs(m)
    if gen_3d and m.GetNumConformers() == 0:
      m_h = Chem.AddHs(m)
      AllChem.EmbedMolecule(m_h, AllChem.ETKDGv3())
      m = m_h
    processed_mols.append(m)

  # Write output
  os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
  if out_format in ("smi", "smiles"):
    with open(output_path, "w") as f:
      for m in processed_mols:
        f.write(f"{Chem.MolToSmiles(m)}\n")
  elif out_format == "sdf":
    with Chem.SDWriter(output_path) as writer:
      for m in processed_mols:
        writer.write(m)
  elif out_format == "pdb":
    Chem.MolToPDBFile(processed_mols[0], output_path)
  elif out_format == "mol":
    Chem.MolToMolFile(processed_mols[0], output_path)

  print(
      f"[RDKit Fallback Engine] Converted {len(processed_mols)} molecule(s)"
      f" from '{in_format}' to '{out_format}' -> {output_path}"
  )
  return True


def convert_file(
    input_path: str,
    output_path: str,
    in_fmt: Optional[str] = None,
    out_fmt: Optional[str] = None,
    add_hs: bool = False,
    gen_3d: bool = False,
) -> bool:
  """Converts chemical file formats via Open Babel with RDKit fallback."""
  if HAS_OPENBABEL and pybel is not None:
    if convert_format_openbabel(
        input_path, output_path, in_fmt, out_fmt, add_hs, gen_3d
    ):
      return True
  return convert_format_rdkit_fallback(
      input_path, output_path, in_fmt, out_fmt, add_hs, gen_3d
  )


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Convert chemical structure file formats using Open Babel or RDKit."
      )
  )
  parser.add_argument(
      "--input", type=str, required=True, help="Input chemical file path."
  )
  parser.add_argument(
      "--output", type=str, required=True, help="Output chemical file path."
  )
  parser.add_argument(
      "--in_fmt",
      type=str,
      default=None,
      help="Explicit input format (e.g. smi, sdf, mol2, pdb).",
  )
  parser.add_argument(
      "--out_fmt",
      type=str,
      default=None,
      help="Explicit output format (e.g. sdf, mol2, pdb, xyz).",
  )
  parser.add_argument(
      "--add_hs", action="store_true", help="Add explicit hydrogen atoms."
  )
  parser.add_argument(
      "--gen_3d", action="store_true", help="Generate 3D coordinates."
  )
  args = parser.parse_args()

  success = convert_file(
      args.input,
      args.output,
      in_fmt=args.in_fmt,
      out_fmt=args.out_fmt,
      add_hs=args.add_hs,
      gen_3d=args.gen_3d,
  )
  if not success:
    sys.exit(1)


if __name__ == "__main__":
  main()
