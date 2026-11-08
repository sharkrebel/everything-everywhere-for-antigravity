# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
# ]
# ///

"""Multi-Engine Biomolecular Co-Folding Input Generator & Pre-Flight Ligand Validator.

Prepares validated co-folding input specifications for:
- AlphaFold 3 (JSON format, default engine)
- Boltz-2 (YAML format, open-source MIT biomolecular diffusion backup)
- Boltz API (`boltz-2.1` JSON format, hosted structure & binding affinity API)
- ESMFold + RDKit/Vina (FASTA + SDF format, single-sequence backup)

Runs pre-flight RDKit 3D conformer checks to prevent pipeline aborts during
featurization.
"""

import argparse
import json
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem


def check_rdkit_ligand_embed(smiles: str) -> Tuple[bool, str]:
  """Pre-flight check: tests whether RDKit can generate a valid 3D conformer."""
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    return (
        False,
        "INVALID_SMILES: RDKit cannot parse the provided SMILES string.",
    )

  mol_h = Chem.AddHs(mol)
  params = AllChem.ETKDGv3()
  params.randomSeed = 42

  res = AllChem.EmbedMolecule(mol_h, params)
  if res >= 0:
    return True, "PASSED: Standard ETKDGv3 embedding successful."

  params.useRandomCoords = True
  params.maxAttempts = 1000
  params.enforceChirality = False
  res2 = AllChem.EmbedMolecule(mol_h, params)
  if res2 >= 0:
    return (
        True,
        (
            "PASSED_WITH_WARNING: Required random coordinates. Suggest setting"
            " --conformer_max_iterations=1000 in AF3."
        ),
    )

  return (
      False,
      (
          "FAILED: RDKit could not embed a 3D conformer. Provide reference"
          " coordinates via userCCDPath or supply a pre-built CCD mmCIF."
      ),
  )


def build_boltz_api_payload(
    proteins: List[Tuple[str, str]],
    ligands_smiles: List[Tuple[str, str]],
    ligands_ccd: List[Tuple[str, List[str]]],
    dna_chains: Optional[List[Tuple[str, str]]] = None,
    rna_chains: Optional[List[Tuple[str, str]]] = None,
    bonded_atom_pairs: Optional[List[Any]] = None,
    msa_free: bool = False,
    binder_chain_id: Optional[str] = None,
    enable_binding: bool = True,
    pocket_residues: Optional[Dict[str, List[int]]] = None,
    pocket_max_distance: float = 6.0,
    pocket_force: bool = False,
    num_samples: int = 1,
) -> Dict[str, Any]:
  """Builds a Boltz API (boltz-2.1) structure-and-binding input payload."""
  entities: List[Dict[str, Any]] = []
  polymer_chain_ids = set()
  smiles_chain_ids = set()
  ligand_chain_ids: List[str] = []

  for cid, seq in proteins:
    polymer_chain_ids.add(cid)
    ent: Dict[str, Any] = {
        "type": "protein",
        "value": seq,
        "chain_ids": [cid],
    }
    if msa_free:
      ent["msa"] = {"type": "empty"}
    entities.append(ent)

  for cid, smi in ligands_smiles:
    smiles_chain_ids.add(cid)
    ligand_chain_ids.append(cid)
    entities.append({
        "type": "ligand_smiles",
        "value": smi,
        "chain_ids": [cid],
    })

  for cid, ccds in ligands_ccd:
    ligand_chain_ids.append(cid)
    entities.append({
        "type": "ligand_ccd",
        "value": ccds[0],
        "chain_ids": [cid],
    })

  if dna_chains:
    for cid, seq in dna_chains:
      polymer_chain_ids.add(cid)
      entities.append({
          "type": "dna",
          "value": seq,
          "chain_ids": [cid],
      })

  if rna_chains:
    for cid, seq in rna_chains:
      polymer_chain_ids.add(cid)
      entities.append({
          "type": "rna",
          "value": seq,
          "chain_ids": [cid],
      })

  payload: Dict[str, Any] = {"entities": entities}

  target_binder = binder_chain_id or (
      ligand_chain_ids[0] if ligand_chain_ids else None
  )
  has_nucleic = bool(dna_chains or rna_chains)
  if enable_binding:
    if proteins and ligand_chain_ids and not has_nucleic and target_binder:
      payload["binding"] = {
          "type": "ligand_protein_binding",
          "binder_chain_id": target_binder,
      }
    elif (
        len(proteins) >= 2
        and not ligand_chain_ids
        and not has_nucleic
        and binder_chain_id
    ):
      payload["binding"] = {
          "type": "protein_protein_binding",
          "binder_chain_ids": [binder_chain_id],
      }

  if pocket_residues and target_binder:
    contact_0idx: Dict[str, List[int]] = {}
    for p_cid, res_list in pocket_residues.items():
      contact_0idx[p_cid] = [max(0, int(r) - 1) for r in res_list]
    payload["constraints"] = [{
        "type": "pocket",
        "binder_chain_id": target_binder,
        "contact_residues": contact_0idx,
        "max_distance_angstrom": float(pocket_max_distance),
        "force": bool(pocket_force),
    }]

  if bonded_atom_pairs:
    bonds = []
    for pair in bonded_atom_pairs:
      if isinstance(pair, (list, tuple)) and len(pair) == 2:
        endpoints = []
        for ep in pair:
          c_id, r_idx, atom_id = str(ep[0]), int(ep[1]), str(ep[2])
          if c_id in polymer_chain_ids:
            endpoints.append({
                "type": "polymer_atom",
                "chain_id": c_id,
                "residue_index": max(0, r_idx - 1),
                "atom_name": atom_id,
            })
          elif c_id in smiles_chain_ids and atom_id.isdigit():
            endpoints.append({
                "type": "smiles_atom",
                "chain_id": c_id,
                "atom_map": int(atom_id),
            })
          else:
            endpoints.append({
                "type": "ligand_atom",
                "chain_id": c_id,
                "atom_name": atom_id,
            })
        bonds.append({"atom1": endpoints[0], "atom2": endpoints[1]})
    if bonds:
      payload["bonds"] = bonds

  if num_samples and int(num_samples) > 1:
    payload["num_samples"] = int(num_samples)

  return payload


def create_cofolding_inputs(
    *,
    name: str,
    proteins: List[Tuple[str, str]],
    ligands_smiles: List[Tuple[str, str]],
    ligands_ccd: List[Tuple[str, List[str]]],
    dna_chains: Optional[List[Tuple[str, str]]] = None,
    rna_chains: Optional[List[Tuple[str, str]]] = None,
    bonded_atom_pairs: Optional[List[Any]] = None,
    seeds: Optional[List[int]] = None,
    user_ccd_path: Optional[str] = None,
    output_json: Optional[str] = None,
    output_path: Optional[str] = None,
    output_boltz_yaml: Optional[str] = None,
    output_boltz_api_json: Optional[str] = None,
    msa_free: bool = False,
    skip_preflight: bool = False,
    binder_chain_id: Optional[str] = None,
    enable_binding: bool = True,
    pocket_residues: Optional[Dict[str, List[int]]] = None,
    pocket_max_distance: float = 6.0,
    pocket_force: bool = False,
    num_samples: int = 1,
) -> Dict[str, Any]:
  """Builds AlphaFold 3 JSON, optional Boltz-2 YAML, and Boltz API JSON co-folding inputs."""
  output_json = output_json or output_path
  preflight_reports = {}
  if not skip_preflight:
    for chain_id, smi in ligands_smiles:
      ok, msg = check_rdkit_ligand_embed(smi)
      preflight_reports[chain_id] = {
          "smiles": smi,
          "passed": ok,
          "message": msg,
      }
      if not ok:
        raise ValueError(
            f"Ligand pre-flight check failed for chain {chain_id}: {msg}"
        )

  sequences = []
  boltz_sequences = []

  for cid, seq in proteins:
    prot_entry: Dict[str, Any] = {"id": cid, "sequence": seq}
    if msa_free:
      prot_entry["unpairedMsa"] = f">{cid}\n{seq}\n"
      prot_entry["pairedMsa"] = ""
      prot_entry["templates"] = []
    sequences.append({"protein": prot_entry})
    boltz_sequences.append({"protein": {"id": cid, "sequence": seq}})

  for cid, smi in ligands_smiles:
    sequences.append({"ligand": {"id": cid, "smiles": smi}})
    boltz_sequences.append({"ligand": {"id": cid, "smiles": smi}})

  for cid, ccds in ligands_ccd:
    sequences.append({"ligand": {"id": cid, "ccdCodes": ccds}})
    boltz_sequences.append({"ligand": {"id": cid, "ccd": ccds[0]}})

  # Nucleic acids must reach BOTH payloads. Emitting them only into the AF3
  # JSON would silently hand Boltz-2 a protein-only system that folds fine
  # and answers a different question than the one asked.
  if dna_chains:
    for cid, seq in dna_chains:
      sequences.append({"dna": {"id": cid, "sequence": seq}})
      boltz_sequences.append({"dna": {"id": cid, "sequence": seq}})
  if rna_chains:
    for cid, seq in rna_chains:
      sequences.append({"rna": {"id": cid, "sequence": seq}})
      boltz_sequences.append({"rna": {"id": cid, "sequence": seq}})

  af3_payload: Dict[str, Any] = {
      "name": name,
      "modelSeeds": seeds or [42],
      "sequences": sequences,
      "dialect": "alphafold3",
      "version": 1,
  }
  if bonded_atom_pairs:
    af3_payload["bondedAtomPairs"] = bonded_atom_pairs
  if user_ccd_path:
    af3_payload["userCCDPath"] = user_ccd_path

  boltz_api_payload = build_boltz_api_payload(
      proteins=proteins,
      ligands_smiles=ligands_smiles,
      ligands_ccd=ligands_ccd,
      dna_chains=dna_chains,
      rna_chains=rna_chains,
      bonded_atom_pairs=bonded_atom_pairs,
      msa_free=msa_free,
      binder_chain_id=binder_chain_id,
      enable_binding=enable_binding,
      pocket_residues=pocket_residues,
      pocket_max_distance=pocket_max_distance,
      pocket_force=pocket_force,
      num_samples=num_samples,
  )

  if output_json:
    os.makedirs(
        os.path.dirname(os.path.abspath(output_json)) or ".", exist_ok=True
    )
    with open(output_json, "w") as f:
      json.dump(af3_payload, f, indent=2)
      f.write("\n")

  if output_boltz_api_json:
    os.makedirs(
        os.path.dirname(os.path.abspath(output_boltz_api_json)) or ".",
        exist_ok=True,
    )
    with open(output_boltz_api_json, "w") as f:
      json.dump(boltz_api_payload, f, indent=2)
      f.write("\n")

  if output_boltz_yaml:
    os.makedirs(
        os.path.dirname(os.path.abspath(output_boltz_yaml)) or ".",
        exist_ok=True,
    )
    yaml_lines = ["version: 1", "sequences:"]
    for item in boltz_sequences:
      if "protein" in item:
        p = item["protein"]
        yaml_lines.append("  - protein:")
        yaml_lines.append(f"      id: {p['id']}")
        yaml_lines.append(f"      sequence: {p['sequence']}")
        if msa_free:
          yaml_lines.append("      msa: empty")
      elif "ligand" in item:
        lig = item["ligand"]
        yaml_lines.append("  - ligand:")
        yaml_lines.append(f"      id: {lig['id']}")
        if "smiles" in lig:
          yaml_lines.append(f"      smiles: '{lig['smiles']}'")
        elif "ccd" in lig:
          yaml_lines.append(f"      ccd: {lig['ccd']}")
      else:
        for kind in ("dna", "rna"):
          if kind in item:
            na = item[kind]
            yaml_lines.append(f"  - {kind}:")
            yaml_lines.append(f"      id: {na['id']}")
            yaml_lines.append(f"      sequence: {na['sequence']}")
    with open(output_boltz_yaml, "w") as f:
      f.write("\n".join(yaml_lines) + "\n")

  return {
      "output_path": output_json,
      "af3_json_path": output_json,
      "boltz_yaml_path": output_boltz_yaml,
      "boltz_api_json_path": output_boltz_api_json,
      "preflight_reports": preflight_reports,
      "payload": af3_payload,
      "boltz_api_payload": boltz_api_payload,
  }


create_af3_input_json = create_cofolding_inputs


# AF3 accepts the unambiguous nucleotide alphabets plus N. Anything else is
# almost always a protein sequence pasted into the wrong flag, or RNA passed
# as DNA -- both of which fold happily and answer the wrong question.
DNA_ALPHABET = set("ACGTN")
RNA_ALPHABET = set("ACGUN")


def parse_chain_specs(
    raw_items: List[List[str]], default_chain_id: str
) -> List[Tuple[str, str]]:
  """Parses repeated ``--flag ID:VALUE`` / ``--flag ID VALUE`` arguments.

  Args:
    raw_items: the argparse ``action="append", nargs="+"`` accumulator.
    default_chain_id: chain ID to assume when only a bare value is given.

  Returns:
    List of ``(chain_id, value)`` pairs.
  """
  out: List[Tuple[str, str]] = []
  for item in raw_items:
    if len(item) == 1 and ":" in item[0]:
      cid, val = item[0].split(":", 1)
      out.append((cid, val))
    elif len(item) == 1:
      out.append((default_chain_id, item[0]))
    elif len(item) >= 2:
      out.append((item[0], item[1]))
  return out


def validate_nucleic_acid_chains(
    chains: List[Tuple[str, str]], alphabet: set, kind: str
) -> List[Tuple[str, str]]:
  """Rejects sequences that are not in the expected nucleotide alphabet.

  Args:
    chains: ``(chain_id, sequence)`` pairs.
    alphabet: permitted single-letter codes.
    kind: ``"DNA"`` or ``"RNA"``, used in the error message.

  Returns:
    The chains, uppercased.

  Raises:
    SystemExit: on any invalid character, naming the offending letters. A
      co-folding run is expensive; failing here is far cheaper than folding
      a mistyped sequence and trusting the result.
  """
  validated: List[Tuple[str, str]] = []
  for cid, seq in chains:
    upper = seq.strip().upper()
    bad = sorted(set(upper) - alphabet)
    if not upper:
      print(f"Error: {kind} chain '{cid}' is empty.", file=sys.stderr)
      sys.exit(1)
    if bad:
      hint = ""
      if kind == "DNA" and "U" in bad:
        hint = " Did you mean --rna?"
      elif kind == "RNA" and "T" in bad:
        hint = " Did you mean --dna?"
      print(
          f"Error: {kind} chain '{cid}' contains characters outside the {kind}"
          f" alphabet {''.join(sorted(alphabet))}: {', '.join(bad)}.{hint}",
          file=sys.stderr,
      )
      sys.exit(1)
    validated.append((cid, upper))
  return validated


def main() -> None:
  parser = argparse.ArgumentParser(
      description="Multi-Engine Biomolecular Co-Folding Input Generator"
  )
  parser.add_argument(
      "--name",
      "--job-name",
      dest="name",
      required=True,
      help="Job name identifier",
  )
  parser.add_argument(
      "--protein",
      action="append",
      default=[],
      nargs="+",
      help="Protein chain ID:SEQUENCE",
  )
  parser.add_argument(
      "--protein-fasta", default=None, help="Path to protein FASTA file"
  )
  parser.add_argument(
      "--ligand-smiles",
      action="append",
      default=[],
      nargs="+",
      help="Ligand chain ID:SMILES",
  )
  parser.add_argument(
      "--ligand-ccd",
      action="append",
      default=[],
      nargs="+",
      help="Ligand CCD code ID:CCD",
  )
  parser.add_argument(
      "--dna",
      action="append",
      default=[],
      nargs="+",
      help=(
          "DNA chain ID:SEQUENCE (repeatable; supply both strands separately"
          " for a duplex)"
      ),
  )
  parser.add_argument(
      "--rna",
      action="append",
      default=[],
      nargs="+",
      help="RNA chain ID:SEQUENCE (repeatable)",
  )
  parser.add_argument(
      "--covalent-bond",
      action="append",
      default=[],
      help="Covalent bond spec C1,R1,A1,C2,R2,A2",
  )
  parser.add_argument(
      "--output-json", required=True, help="Output AlphaFold 3 JSON file path"
  )
  parser.add_argument(
      "--output-boltz-yaml",
      default=None,
      help="Optional output Boltz-2 YAML file path",
  )
  parser.add_argument(
      "--output-boltz-api-json",
      default=None,
      help="Optional output Boltz API (boltz-2.1) JSON file path",
  )
  parser.add_argument(
      "--binder-chain-id",
      default=None,
      help="Binder chain ID for Boltz API binding affinity evaluation",
  )
  parser.add_argument(
      "--disable-binding",
      action="store_true",
      help="Disable automatic Boltz API binding affinity calculation",
  )
  parser.add_argument(
      "--pocket-residues",
      action="append",
      default=[],
      help=(
          "Pocket constraint specification CHAIN:RES1,RES2,... (1-based"
          " residue numbers; converted to 0-based for Boltz API)"
      ),
  )
  parser.add_argument(
      "--num-samples",
      type=int,
      default=1,
      help="Number of diffusion structure samples for Boltz API (default: 1)",
  )
  parser.add_argument(
      "--msa-free",
      action="store_true",
      help="Configure single-sequence MSA-free mode (no 1.5 TB DB needed)",
  )
  parser.add_argument("--seed", type=int, default=42, help="Random seed")
  parser.add_argument(
      "--skip-preflight",
      action="store_true",
      help="Skip RDKit 3D ligand pre-flight check",
  )
  args = parser.parse_args()

  proteins: List[Tuple[str, str]] = []
  if args.protein_fasta:
    if os.path.exists(args.protein_fasta):
      with open(args.protein_fasta, "r") as f:
        seq = "".join(
            l.strip() for l in f if l.strip() and not l.strip().startswith(">")
        )
        proteins.append(("A", seq))
    else:
      seq = "".join(
          l.strip()
          for l in args.protein_fasta.splitlines()
          if l.strip() and not l.strip().startswith(">")
      )
      proteins.append(("A", seq))

  proteins.extend(parse_chain_specs(args.protein, "A"))

  dna_chains = validate_nucleic_acid_chains(
      parse_chain_specs(args.dna, "C"), DNA_ALPHABET, "DNA"
  )
  rna_chains = validate_nucleic_acid_chains(
      parse_chain_specs(args.rna, "D"), RNA_ALPHABET, "RNA"
  )

  ligands_smiles = parse_chain_specs(args.ligand_smiles, "B")

  ligands_ccd: List[Tuple[str, List[str]]] = [
      (cid, [ccd.upper()])
      for cid, ccd in parse_chain_specs(args.ligand_ccd, "B")
  ]

  bonded_atom_pairs = []
  for bond_str in args.covalent_bond:
    delim = "," if "," in bond_str else ":"
    parts = [p.strip() for p in bond_str.split(delim)]
    if len(parts) == 6:
      c1, r1, a1, c2, r2, a2 = parts
      bonded_atom_pairs.append([[c1, int(r1), a1], [c2, int(r2), a2]])

  pocket_residues: Dict[str, List[int]] = {}
  for p_spec in args.pocket_residues:
    if ":" in p_spec:
      p_cid, r_csv = p_spec.split(":", 1)
    else:
      p_cid, r_csv = "A", p_spec
    res_nums = [int(x.strip()) for x in r_csv.split(",") if x.strip()]
    if res_nums:
      pocket_residues.setdefault(p_cid.strip(), []).extend(res_nums)

  res = create_cofolding_inputs(
      name=args.name,
      proteins=proteins,
      ligands_smiles=ligands_smiles,
      ligands_ccd=ligands_ccd,
      dna_chains=dna_chains or None,
      rna_chains=rna_chains or None,
      bonded_atom_pairs=bonded_atom_pairs if bonded_atom_pairs else None,
      seeds=[args.seed],
      output_json=args.output_json,
      output_boltz_yaml=args.output_boltz_yaml,
      output_boltz_api_json=args.output_boltz_api_json,
      msa_free=args.msa_free,
      skip_preflight=args.skip_preflight,
      binder_chain_id=args.binder_chain_id,
      enable_binding=not args.disable_binding,
      pocket_residues=pocket_residues or None,
      num_samples=args.num_samples,
  )
  print(f"Successfully generated co-folding input JSON: {res['af3_json_path']}")
  if res.get("boltz_yaml_path"):
    print(f"Successfully generated Boltz-2 YAML: {res['boltz_yaml_path']}")
  if res.get("boltz_api_json_path"):
    print(
        f"Successfully generated Boltz API JSON: {res['boltz_api_json_path']}"
    )


if __name__ == "__main__":
  main()
