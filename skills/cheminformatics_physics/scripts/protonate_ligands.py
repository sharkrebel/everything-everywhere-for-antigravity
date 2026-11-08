# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "gemmi",
#     "meeko",
#     "openbabel-wheel",
#     "pandas",
#     "rdkit",
#     "scipy",
# ]
# ///

"""Standalone pH-Dependent Small-Molecule Protonation & Microstate Engine.

Assigns predominant pH-dependent ionization states and explicit 3D hydrogens to
small molecules using Open Babel's `OBPhModel::CorrectForPH` (`phmodel.txt`),
producing protonated `.sdf`, `.csv`, `.pdbqt`, or `.json` files that can be
consumed directly by downstream docking, conformer, or MD scripts.

Capabilities:
1. Open Babel `OBPhModel::CorrectForPH(ph)` ionization across arbitrary target
   pH values with 1-to-1 heavy-atom formal-charge mapping back onto RDKit
   molecules (preserving stereochemistry, Kekulization, and 3D coordinates).
2. Optional pH-window microstate enumeration (`--enumerate-microstates
   --ph-window 1.5`) for ionizable sites whose charge state transitions near the
   target pH.
3. 3D coordinate preservation for input `.sdf`/`.mol` files or fresh ETKDGv3 +
   MMFF94s 3D conformer generation (`--embed-3d`) with explicit hydrogens.
4. Direct export to `.sdf`, `.csv`, `.pdbqt` (via Meeko), and `.json` reports.
"""

import argparse
import gzip
import itertools
import json
import os
import sys
from typing import Any, Dict, List, Set, Tuple

from openbabel import pybel
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Geometry import Point3D

try:
  import meeko
except ImportError:
  meeko = None

_HYDROGEN_ATOMIC_NUM = 1
_NITROGEN_ATOMIC_NUM = 7

_FILE_EXTENSIONS = (
    ".csv",
    ".sdf",
    ".sdf.gz",
    ".mol",
    ".smi",
    ".smiles",
    ".txt",
)


def _desalt_molecule(mol: Chem.Mol) -> Chem.Mol:
  """Extracts the largest covalent fragment to strip counterions/salts."""
  try:
    lfc = rdMolStandardize.LargestFragmentChooser()
    clean = lfc.choose(mol)
    if clean is not None and clean.GetNumAtoms() > 0:
      return clean
  except Exception:
    pass
  return mol


def _uncharge_preserving_quaternary(mol: Chem.Mol) -> Chem.Mol:
  """Neutralizes reversible formal charges while preserving permanent cations."""
  mol_no_h = Chem.RemoveHs(mol)
  try:
    uncharger = rdMolStandardize.Uncharger()
    neutral = uncharger.uncharge(mol_no_h)
    Chem.SanitizeMol(neutral)
    return neutral
  except Exception as exc:
    print(
        f"WARNING: Uncharger sanitization failed ({exc}); retaining input"
        " formal charges.",
        file=sys.stderr,
    )
    return mol_no_h


def _apply_atom_charge_map_to_rdkit(
    base_mol: Chem.Mol,
    rd_idx_to_new_charge: Dict[int, int],
) -> Tuple[Chem.Mol, List[Dict[str, Any]]]:
  """Applies atom-level formal charge assignments to an RDKit molecule."""
  rw = Chem.RWMol(base_mol)
  atom_rings = base_mol.GetRingInfo().AtomRings()
  ionized_sites: List[Dict[str, Any]] = []

  for rd_idx, new_q in sorted(rd_idx_to_new_charge.items()):
    orig_atom = base_mol.GetAtomWithIdx(rd_idx)
    old_q = orig_atom.GetFormalCharge()
    if new_q == old_q:
      continue

    target_idx = rd_idx
    # When Open Babel deprotonates a 5-membered aromatic ring (e.g. tetrazole),
    # ensure the negative charge is placed on the [nH]-bearing ring nitrogen
    # so RDKit Kekulization remains valid.
    if (
        new_q < old_q
        and orig_atom.GetIsAromatic()
        and orig_atom.GetTotalNumHs() == 0
    ):
      for ring in atom_rings:
        if rd_idx in ring:
          for r_idx in ring:
            r_orig = base_mol.GetAtomWithIdx(r_idx)
            r_cur = rw.GetAtomWithIdx(r_idx)
            if (
                r_orig.GetIsAromatic()
                and r_orig.GetAtomicNum() == _NITROGEN_ATOMIC_NUM
                and r_orig.GetTotalNumHs() > 0
                and r_cur.GetFormalCharge() == 0
            ):
              target_idx = r_idx
              break

    tgt = rw.GetAtomWithIdx(target_idx)
    tgt.SetFormalCharge(new_q)
    if tgt.GetIsAromatic() and tgt.GetAtomicNum() == _NITROGEN_ATOMIC_NUM:
      tgt.SetNumExplicitHs(1 if new_q > 0 else 0)
      tgt.SetNoImplicit(True)
    else:
      tgt.SetNumExplicitHs(0)
      tgt.SetNoImplicit(False)

    ionized_sites.append({
        "atom_idx": target_idx,
        "element": tgt.GetSymbol(),
        "old_charge": old_q,
        "new_charge": new_q,
    })

  Chem.SanitizeMol(rw)
  return rw.GetMol(), ionized_sites


def _openbabel_charge_map_at_ph(
    base_mol: Chem.Mol,
    ph: float,
) -> Dict[int, int]:
  """Computes RDKit atom-index -> formal-charge map using Open Babel OBPhModel."""
  smi = Chem.MolToSmiles(base_mol, isomericSmiles=True)
  props = base_mol.GetPropsAsDict(includePrivate=True, includeComputed=True)
  atom_order = list(
      props.get("_smilesAtomOutputOrder", range(base_mol.GetNumAtoms()))
  )

  ob_mol = pybel.readstring("smi", smi).OBMol
  ob_mol.CorrectForPH(float(ph))

  charge_map: Dict[int, int] = {}
  for smi_pos, rd_idx in enumerate(atom_order):
    ob_atom = ob_mol.GetAtom(smi_pos + 1)
    if ob_atom is not None:
      charge_map[rd_idx] = int(ob_atom.GetFormalCharge())
  return charge_map


def _protonate_with_openbabel(
    base_mol: Chem.Mol,
    ph: float,
    enumerate_microstates: bool = False,
    ph_window: float = 1.5,
    max_microstates: int = 8,
) -> Tuple[Chem.Mol, List[Dict[str, Any]], List[Dict[str, Any]]]:
  """Protonates an RDKit Mol using Open Babel's OBPhModel::CorrectForPH."""
  pred_map = _openbabel_charge_map_at_ph(base_mol, ph)
  mol_ion, ionized_sites = _apply_atom_charge_map_to_rdkit(base_mol, pred_map)
  pred_smi = Chem.MolToSmiles(mol_ion, isomericSmiles=True)
  pred_q = sum(a.GetFormalCharge() for a in mol_ion.GetAtoms())

  microstates: List[Dict[str, Any]] = [{
      "protonated_smiles": pred_smi,
      "net_charge": pred_q,
      "is_predominant": True,
  }]

  if enumerate_microstates and ph_window > 0:
    ph_low = max(0.0, ph - ph_window)
    ph_high = min(14.0, ph + ph_window)
    n_steps = max(3, int(round((ph_high - ph_low) / 0.25)) + 1)
    step_size = (ph_high - ph_low) / max(1, n_steps - 1)
    test_phs = [ph_low + k * step_size for k in range(n_steps)]

    atom_charge_choices: Dict[int, Set[int]] = {
        i: {pred_map.get(i, base_mol.GetAtomWithIdx(i).GetFormalCharge())}
        for i in range(base_mol.GetNumAtoms())
    }
    for t_ph in test_phs:
      cmap = _openbabel_charge_map_at_ph(base_mol, t_ph)
      for rd_idx, q_val in cmap.items():
        atom_charge_choices[rd_idx].add(q_val)

    variable_atoms = [
        idx
        for idx, q_set in sorted(atom_charge_choices.items())
        if len(q_set) > 1
    ]
    if variable_atoms:
      seen_smiles = {pred_smi}
      choice_lists = [
          sorted(atom_charge_choices[idx]) for idx in variable_atoms
      ]
      for combo in itertools.product(*choice_lists):
        if len(microstates) >= max_microstates:
          break
        cand_map = dict(pred_map)
        for idx, q_val in zip(variable_atoms, combo):
          cand_map[idx] = q_val
        try:
          cand_mol, _ = _apply_atom_charge_map_to_rdkit(base_mol, cand_map)
          cand_smi = Chem.MolToSmiles(cand_mol, isomericSmiles=True)
        except Exception:
          continue
        if cand_smi in seen_smiles:
          continue
        seen_smiles.add(cand_smi)
        cand_q = sum(a.GetFormalCharge() for a in cand_mol.GetAtoms())
        microstates.append({
            "protonated_smiles": cand_smi,
            "net_charge": cand_q,
            "is_predominant": False,
        })

  return mol_ion, ionized_sites, microstates


def _embed_or_preserve_3d(
    source_mol: Chem.Mol,
    ionized_mol: Chem.Mol,
    embed_3d: bool = False,
    random_seed: int = 42,
) -> Chem.Mol:
  """Preserves input 3D coordinates or generates an ETKDGv3 3D conformer."""
  has_source_3d = (
      source_mol.GetNumConformers() > 0 and source_mol.GetConformer().Is3D()
  )
  if not has_source_3d and not embed_3d:
    return ionized_mol

  if has_source_3d:
    src_no_h = Chem.RemoveHs(source_mol)
    if src_no_h.GetNumAtoms() == ionized_mol.GetNumAtoms():
      src_conf = src_no_h.GetConformer()
      new_conf = Chem.Conformer(ionized_mol.GetNumAtoms())
      new_conf.Set3D(True)
      for i in range(ionized_mol.GetNumAtoms()):
        pt = src_conf.GetAtomPosition(i)
        new_conf.SetAtomPosition(i, Point3D(pt.x, pt.y, pt.z))
      mol_3d = Chem.Mol(ionized_mol)
      mol_3d.RemoveAllConformers()
      mol_3d.AddConformer(new_conf, assignId=True)
      mol_h = Chem.AddHs(mol_3d, addCoords=True)
      try:
        mp = AllChem.MMFFGetMoleculeProperties(mol_h, mmffVariant="MMFF94s")
        if mp:
          ff = AllChem.MMFFGetMoleculeForceField(mol_h, mp)
          if ff:
            for a in mol_h.GetAtoms():
              if a.GetAtomicNum() > _HYDROGEN_ATOMIC_NUM:
                ff.AddFixedPoint(a.GetIdx())
            ff.Initialize()
            ff.Minimize(maxIts=200)
      except Exception:
        pass
      return mol_h

  mol_h = Chem.AddHs(ionized_mol)
  ri = ionized_mol.GetRingInfo()
  is_macro = any(len(r) >= 12 for r in ri.AtomRings())
  params = (
      AllChem.ETKDGv3() if hasattr(AllChem, "ETKDGv3") else AllChem.ETKDGv2()
  )
  params.randomSeed = random_seed
  if hasattr(params, "useSmallRingTorsions"):
    params.useSmallRingTorsions = True
  if is_macro and hasattr(params, "useMacrocycleTorsions"):
    params.useMacrocycleTorsions = True
  if AllChem.EmbedMolecule(mol_h, params) != 0:
    params.useRandomCoords = True
    AllChem.EmbedMolecule(mol_h, params)
  try:
    if AllChem.MMFFHasAllMoleculeParams(mol_h):
      AllChem.MMFFOptimizeMolecule(mol_h, mmffVariant="MMFF94s", maxIters=500)
    else:
      AllChem.UFFOptimizeMolecule(mol_h, maxIters=500)
  except Exception:
    pass
  return mol_h


def protonate_ligand_mol(
    mol: Chem.Mol,
    ph: float = 7.4,
    enumerate_microstates: bool = False,
    ph_window: float = 1.5,
    max_microstates: int = 8,
    embed_3d: bool = False,
    desalt: bool = True,
    uncharge_first: bool = True,
    random_seed: int = 42,
) -> Dict[str, Any]:
  """Protonates an RDKit Mol to its predominant microstate(s) at target pH.

  Args:
    mol: Input RDKit Mol.
    ph: Target solution pH (default: 7.4).
    enumerate_microstates: If True, enumerates alternative microstates across
      [ph - ph_window, ph + ph_window].
    ph_window: pH window around target pH for microstate enumeration.
    max_microstates: Maximum number of microstates to enumerate.
    embed_3d: If True, ensures the returned molecule has 3D coordinates and
      explicit hydrogens.
    desalt: If True, strips salt counterions before protonation.
    uncharge_first: If True, neutralizes reversible charges before applying
      pH-dependent site ionization.
    random_seed: Deterministic PRNG seed for RDKit ETKDGv3 3D conformer
      embedding so generated coordinates are reproducible across runs.

  Returns:
    Dictionary containing the predominant RDKit `mol`, `input_smiles`,
    `protonated_smiles`, `net_charge`, `positive_charges`, `negative_charges`,
    `is_zwitterion`, `species_class`, `ionized_sites`, and `microstates`.
  """
  if mol is None:
    raise ValueError("Input RDKit Mol is None.")

  working_mol = _desalt_molecule(mol) if desalt else Chem.Mol(mol)
  input_smiles = Chem.MolToSmiles(
      Chem.RemoveHs(working_mol), isomericSmiles=True
  )
  if uncharge_first:
    base_mol = _uncharge_preserving_quaternary(working_mol)
  else:
    base_mol = Chem.RemoveHs(working_mol)
    Chem.SanitizeMol(base_mol)
  clogp_neutral = round(float(Descriptors.MolLogP(base_mol)), 3)

  mol_ion, ionized_sites, microstates = _protonate_with_openbabel(
      base_mol,
      ph=ph,
      enumerate_microstates=enumerate_microstates,
      ph_window=ph_window,
      max_microstates=max_microstates,
  )

  protonated_smiles = Chem.MolToSmiles(mol_ion, isomericSmiles=True)
  pos_count = sum(1 for a in mol_ion.GetAtoms() if a.GetFormalCharge() > 0)
  neg_count = sum(1 for a in mol_ion.GetAtoms() if a.GetFormalCharge() < 0)
  net_charge = sum(a.GetFormalCharge() for a in mol_ion.GetAtoms())
  is_zwitterion = bool(pos_count > 0 and neg_count > 0)
  if is_zwitterion:
    species_class = "zwitterion"
  elif pos_count > 0:
    species_class = "cation"
  elif neg_count > 0:
    species_class = "anion"
  else:
    species_class = "neutral"

  mol_final = _embed_or_preserve_3d(
      working_mol, mol_ion, embed_3d=embed_3d, random_seed=random_seed
  )
  if mol.HasProp("_Name"):
    mol_final.SetProp("_Name", mol.GetProp("_Name"))

  return {
      "mol": mol_final,
      "input_smiles": input_smiles,
      "protonated_smiles": protonated_smiles,
      "ph": ph,
      "net_charge": net_charge,
      "positive_charges": pos_count,
      "negative_charges": neg_count,
      "is_zwitterion": is_zwitterion,
      "species_class": species_class,
      "cLogP_neutral": clogp_neutral,
      "ionized_sites": ionized_sites,
      "microstates": microstates,
  }


def protonate_smiles(
    smiles: str,
    ph: float = 7.4,
    enumerate_microstates: bool = False,
    ph_window: float = 1.5,
    max_microstates: int = 8,
    embed_3d: bool = False,
    desalt: bool = True,
    uncharge_first: bool = True,
    random_seed: int = 42,
) -> Dict[str, Any]:
  """Protonates a SMILES string at target pH using Open Babel OBPhModel.

  Args:
    smiles: Input SMILES string.
    ph: Target solution pH (default: 7.4).
    enumerate_microstates: If True, enumerates alternative microstates across
      [ph - ph_window, ph + ph_window].
    ph_window: pH window around target pH for microstate enumeration.
    max_microstates: Maximum number of microstates to enumerate.
    embed_3d: If True, generates 3D coordinates and explicit hydrogens.
    desalt: If True, strips salt counterions before protonation.
    uncharge_first: If True, neutralizes reversible charges before applying
      pH-dependent site ionization.
    random_seed: Deterministic PRNG seed for RDKit ETKDGv3 3D conformer
      embedding so generated coordinates are reproducible across runs.

  Returns:
    Dictionary of protonated RDKit Mol and ionization metadata from
    `protonate_ligand_mol`.
  """
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    raise ValueError(f"Invalid SMILES string: {smiles!r}")
  return protonate_ligand_mol(
      mol,
      ph=ph,
      enumerate_microstates=enumerate_microstates,
      ph_window=ph_window,
      max_microstates=max_microstates,
      embed_3d=embed_3d,
      desalt=desalt,
      uncharge_first=uncharge_first,
      random_seed=random_seed,
  )


def load_input_molecules(
    input_path: str,
    smiles_col: str = "smiles",
    id_col: str = "compound_id",
) -> List[Tuple[str, Chem.Mol, str]]:
  """Loads (compound_id, mol, smiles) records from CSV, SDF, MOL, or SMI files."""
  if not os.path.exists(input_path):
    raise FileNotFoundError(f"Input file not found: {input_path}")

  lower_path = input_path.lower()
  out: List[Tuple[str, Chem.Mol, str]] = []

  if lower_path.endswith(".csv"):
    df = pd.read_csv(input_path)
    if smiles_col not in df.columns:
      raise ValueError(
          f"Could not find SMILES column {smiles_col!r} in {input_path}."
          f" Available columns: {list(df.columns)}"
      )
    has_id_col = id_col in df.columns
    for idx, row in df.iterrows():
      raw_val = row[smiles_col]
      if pd.isna(raw_val):
        continue
      smi = str(raw_val).strip()
      if not smi:
        continue
      mol = Chem.MolFromSmiles(smi)
      if mol is None:
        continue
      cid = (
          str(row[id_col]).strip()
          if has_id_col and not pd.isna(row[id_col])
          else f"compound_{idx+1:03d}"
      )
      out.append((cid, mol, Chem.MolToSmiles(mol, isomericSmiles=True)))
    return out

  if lower_path.endswith((".sdf", ".sdf.gz", ".mol")):
    if lower_path.endswith(".sdf.gz"):
      with gzip.open(input_path, "rb") as gz:
        suppl = list(Chem.ForwardSDMolSupplier(gz, removeHs=False))
    elif lower_path.endswith(".sdf"):
      suppl = list(Chem.SDMolSupplier(input_path, removeHs=False))
    else:
      suppl = [Chem.MolFromMolFile(input_path, removeHs=False)]

    for idx, mol in enumerate(suppl):
      if mol is None:
        continue
      cid = (
          mol.GetProp("_Name").strip()
          if mol.HasProp("_Name") and mol.GetProp("_Name").strip()
          else f"compound_{idx+1:03d}"
      )
      for p in ("compound_id", "Compound_ID", "ID", "id"):
        if mol.HasProp(p) and mol.GetProp(p).strip():
          cid = mol.GetProp(p).strip()
          break
      smi = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
      out.append((cid, mol, smi))
    return out

  if lower_path.endswith((".smi", ".smiles", ".txt")):
    with open(input_path, "r", encoding="utf-8") as f:
      for idx, line in enumerate(f):
        parts = line.strip().split()
        if not parts or parts[0].startswith("#"):
          continue
        smi = parts[0]
        cid = parts[1] if len(parts) > 1 else f"compound_{idx+1:03d}"
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
          continue
        out.append((cid, mol, Chem.MolToSmiles(mol, isomericSmiles=True)))
    return out

  raise ValueError(f"Unsupported molecule input format: {input_path}")


def write_ligand_pdbqt(mol: Chem.Mol, output_pdbqt: str) -> str:
  """Writes a protonated RDKit Mol with 3D coordinates to an AutoDock PDBQT."""
  if meeko is None:
    raise ImportError(
        "meeko is required for --output-pdbqt export in protonate_ligands.py."
    )
  has_h = any(a.GetAtomicNum() == _HYDROGEN_ATOMIC_NUM for a in mol.GetAtoms())
  mol_h = mol if has_h else Chem.AddHs(mol)
  if mol_h.GetNumConformers() == 0:
    mol_h = _embed_or_preserve_3d(mol_h, Chem.RemoveHs(mol_h), embed_3d=True)

  ri = mol_h.GetRingInfo()
  is_macro = any(len(r) >= 12 for r in ri.AtomRings())
  preparator = meeko.MoleculePreparation(rigid_macrocycles=is_macro)
  mol_setups = preparator.prepare(mol_h)
  if not mol_setups:
    raise ValueError("Meeko could not prepare ligand setups.")

  pdbqt_str, is_ok, msg = meeko.PDBQTWriterLegacy.write_string(mol_setups[0])
  if not is_ok:
    raise ValueError(f"Meeko PDBQT preparation failed: {msg}")

  lines = []
  for line in pdbqt_str.splitlines():
    if line.startswith(("ATOM", "HETATM")):
      parts = line.split()
      if len(parts) >= 3:
        atype = parts[-1]
        if atype in ("CG0", "C0", "G0"):
          line = line[:77] + "C  " + line[80:]
        elif atype in ("OG0", "O0"):
          line = line[:77] + "OA " + line[80:]
        elif atype in ("NG0", "N0"):
          line = line[:77] + "NA " + line[80:]
        elif atype in ("SG0", "S0"):
          line = line[:77] + "SA " + line[80:]
    lines.append(line)

  os.makedirs(
      os.path.dirname(os.path.abspath(output_pdbqt)) or ".", exist_ok=True
  )
  with open(output_pdbqt, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
  return output_pdbqt


def main() -> None:
  parser = argparse.ArgumentParser(
      description=(
          "Standalone pH-dependent small-molecule protonation via Open Babel"
          " OBPhModel::CorrectForPH."
      )
  )
  parser.add_argument(
      "--smiles",
      type=str,
      default=None,
      help="Single input SMILES string to protonate.",
  )
  parser.add_argument(
      "--input",
      "--input-csv",
      "--input_csv",
      "--input-sdf",
      "--input_sdf",
      dest="input_path",
      type=str,
      default=None,
      help="Input SMILES string or .csv, .sdf, .smi, or .mol file.",
  )
  parser.add_argument(
      "--smiles-col",
      "--smiles_col",
      dest="smiles_col",
      type=str,
      default="smiles",
      help="SMILES column name when --input is a CSV (default: smiles).",
  )
  parser.add_argument(
      "--id-col",
      "--id_col",
      dest="id_col",
      type=str,
      default="compound_id",
      help=(
          "Compound ID column name when --input is a CSV (default:"
          " compound_id)."
      ),
  )
  parser.add_argument(
      "--ph",
      type=float,
      default=7.4,
      help="Target solution pH for ionization (default: 7.4).",
  )
  parser.add_argument(
      "--enumerate-microstates",
      "--enumerate_microstates",
      dest="enumerate_microstates",
      action="store_true",
      help=(
          "Enumerate alternative protonation microstates for sites"
          " transitioning within --ph-window of --ph."
      ),
  )
  parser.add_argument(
      "--ph-window",
      "--ph_window",
      dest="ph_window",
      type=float,
      default=1.5,
      help="pH window around --ph for microstate enumeration (default: 1.5).",
  )
  parser.add_argument(
      "--max-microstates",
      "--max_microstates",
      dest="max_microstates",
      type=int,
      default=8,
      help="Maximum microstates per molecule when enumerating (default: 8).",
  )
  parser.add_argument(
      "--embed-3d",
      "--embed_3d",
      dest="embed_3d",
      action="store_true",
      help="Generate or preserve 3D coordinates with explicit hydrogens.",
  )
  parser.add_argument(
      "--no-desalt",
      "--no_desalt",
      dest="no_desalt",
      action="store_true",
      help="Retain counterions/salts rather than extracting largest fragment.",
  )
  parser.add_argument(
      "--keep-input-charges",
      "--keep_input_charges",
      dest="keep_input_charges",
      action="store_true",
      help=(
          "Preserve pre-existing formal charges instead of neutralizing first."
      ),
  )
  parser.add_argument(
      "--output",
      type=str,
      default=None,
      help="Primary output file path (.sdf, .csv, .pdbqt, or .json).",
  )
  parser.add_argument(
      "--output-sdf",
      "--output_sdf",
      dest="output_sdf",
      type=str,
      default=None,
      help="Output SDF file path for protonated molecules.",
  )
  parser.add_argument(
      "--output-csv",
      "--output_csv",
      dest="output_csv",
      type=str,
      default=None,
      help="Output CSV file path with protonated SMILES and charge properties.",
  )
  parser.add_argument(
      "--output-pdbqt",
      "--output_pdbqt",
      dest="output_pdbqt",
      type=str,
      default=None,
      help="Output AutoDock PDBQT file path (or directory for batch input).",
  )
  parser.add_argument(
      "--summary-json",
      "--summary_json",
      "--output-json",
      "--output_json",
      dest="summary_json",
      type=str,
      default=None,
      help="Output JSON summary report path.",
  )
  args = parser.parse_args()

  if not args.smiles and not args.input_path:
    parser.error("Provide either --smiles or --input (.csv/.sdf/.smi/.mol).")

  if args.output:
    out_lower = args.output.lower()
    if out_lower.endswith(".sdf") and not args.output_sdf:
      args.output_sdf = args.output
    elif out_lower.endswith(".csv") and not args.output_csv:
      args.output_csv = args.output
    elif out_lower.endswith(".pdbqt") and not args.output_pdbqt:
      args.output_pdbqt = args.output
    elif out_lower.endswith(".json") and not args.summary_json:
      args.summary_json = args.output
    elif not args.output_csv:
      args.output_csv = args.output

  need_3d = bool(args.embed_3d or args.output_sdf or args.output_pdbqt)
  desalt = not args.no_desalt
  uncharge_first = not args.keep_input_charges

  records: List[Dict[str, Any]] = []
  mols_out: List[Chem.Mol] = []

  single_smiles = args.smiles
  if not single_smiles and args.input_path:
    looks_like_file = (
        os.path.exists(args.input_path)
        or "/" in args.input_path
        or "\\" in args.input_path
        or args.input_path.lower().endswith(_FILE_EXTENSIONS)
    )
    if not looks_like_file:
      single_smiles = args.input_path

  if single_smiles:
    try:
      res = protonate_smiles(
          single_smiles,
          ph=args.ph,
          enumerate_microstates=args.enumerate_microstates,
          ph_window=args.ph_window,
          max_microstates=args.max_microstates,
          embed_3d=need_3d,
          desalt=desalt,
          uncharge_first=uncharge_first,
      )
    except (ValueError, RuntimeError) as exc:
      print(f"Error: {exc}", file=sys.stderr)
      sys.exit(1)

    mol_out = res.pop("mol")
    mol_out.SetProp("_Name", "query_ligand")
    mol_out.SetProp("compound_id", "query_ligand")
    mol_out.SetProp("input_smiles", res["input_smiles"])
    mol_out.SetProp("protonated_smiles", res["protonated_smiles"])
    mol_out.SetProp("ph", f"{args.ph:.2f}")
    mol_out.SetProp("net_charge", str(res["net_charge"]))
    res["compound_id"] = "query_ligand"
    records.append(res)
    mols_out.append(mol_out)
  else:
    try:
      loaded = load_input_molecules(
          args.input_path, smiles_col=args.smiles_col, id_col=args.id_col
      )
    except (FileNotFoundError, ValueError) as exc:
      print(f"Error: {exc}", file=sys.stderr)
      sys.exit(1)

    if not loaded:
      print(
          f"Error: No valid molecules found in {args.input_path}",
          file=sys.stderr,
      )
      sys.exit(1)

    for cid, mol_in, _ in loaded:
      try:
        res = protonate_ligand_mol(
            mol_in,
            ph=args.ph,
            enumerate_microstates=args.enumerate_microstates,
            ph_window=args.ph_window,
            max_microstates=args.max_microstates,
            embed_3d=need_3d,
            desalt=desalt,
            uncharge_first=uncharge_first,
        )
      except (ValueError, RuntimeError) as exc:
        print(f"Error on {cid}: {exc}", file=sys.stderr)
        sys.exit(1)

      mol_out = res.pop("mol")
      mol_out.SetProp("_Name", cid)
      mol_out.SetProp("compound_id", cid)
      mol_out.SetProp("input_smiles", res["input_smiles"])
      mol_out.SetProp("protonated_smiles", res["protonated_smiles"])
      mol_out.SetProp("ph", f"{args.ph:.2f}")
      mol_out.SetProp("net_charge", str(res["net_charge"]))

      rec = {
          "compound_id": cid,
          "input_smiles": res["input_smiles"],
          "protonated_smiles": res["protonated_smiles"],
          "SMILES_pH": res["protonated_smiles"],
          "ph": args.ph,
          "net_charge": res["net_charge"],
          "positive_charges": res["positive_charges"],
          "negative_charges": res["negative_charges"],
          "is_zwitterion": res["is_zwitterion"],
          "species_class": res["species_class"],
          "cLogP_neutral": res["cLogP_neutral"],
          "ionized_sites": json.dumps(res["ionized_sites"]),
      }
      records.append(rec)
      mols_out.append(mol_out)

  if args.output_csv:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_csv)) or ".", exist_ok=True
    )
    pd.DataFrame(records).to_csv(args.output_csv, index=False)
    print(
        f"Wrote protonated CSV ({len(records)} records) to: {args.output_csv}"
    )

  if args.output_sdf:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_sdf)) or ".", exist_ok=True
    )
    with Chem.SDWriter(args.output_sdf) as writer:
      for m in mols_out:
        writer.write(m)
    print(
        f"Wrote protonated 3D SDF ({len(mols_out)} molecules) to:"
        f" {args.output_sdf}"
    )

  if args.output_pdbqt:
    if not mols_out:
      print(
          "Error: No valid molecule available to write PDBQT.", file=sys.stderr
      )
      sys.exit(1)
    if len(mols_out) == 1 and args.output_pdbqt.lower().endswith(".pdbqt"):
      write_ligand_pdbqt(mols_out[0], args.output_pdbqt)
      print(f"Wrote protonated PDBQT to: {args.output_pdbqt}")
    else:
      out_dir = (
          os.path.dirname(args.output_pdbqt)
          if args.output_pdbqt.lower().endswith(".pdbqt")
          else args.output_pdbqt
      ) or "."
      os.makedirs(out_dir, exist_ok=True)
      for m in mols_out:
        cid = m.GetProp("compound_id") if m.HasProp("compound_id") else "ligand"
        write_ligand_pdbqt(m, os.path.join(out_dir, f"{cid}.pdbqt"))
      print(f"Wrote {len(mols_out)} protonated PDBQT files to: {out_dir}/")

  if args.summary_json:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.summary_json)) or ".",
        exist_ok=True,
    )
    with open(args.summary_json, "w", encoding="utf-8") as f:
      json.dump(records, f, indent=2)
    print(f"Wrote protonation summary JSON to: {args.summary_json}")

  if not (
      args.output_csv
      or args.output_sdf
      or args.output_pdbqt
      or args.summary_json
  ):
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
  main()
