# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "boltz-api",
#     "gemmi",
#     "numpy",
#     "openprotein-python",
#     "polite-http",
#     "python-dotenv",
#     "rdkit",
# ]
# ///

"""Multi-Engine Biomolecular Co-Folding Runner (AlphaFold 3, Boltz-2, Boltz API, ESMFold2).

Generates and executes co-folding inference commands across local Python
environments, Docker/Podman containers, hosted Boltz API (boltz-2.1), and
OpenProtein APIs.
"""

import argparse
import glob
import importlib
import importlib.util
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from typing import Any, Dict, List, Optional

import dotenv
import numpy as np
from polite_http import http_client

try:
  import openprotein
  from openprotein import molecules as op_molecules
except ImportError:
  openprotein = None
  op_molecules = None

dotenv.load_dotenv(os.path.expanduser("~/.env"))

_DEFAULT_AF3_MODEL_DIRS = (
    os.path.expanduser("~/.cache/alphafold3"),
    os.path.expanduser("~/.local/share/alphafold3"),
    "./models",
)


def prerelax_ligand_coordinates(
    lig_coords: np.ndarray,
    prot_coords: np.ndarray,
    clash_cutoff: float = 2.0,
    max_iters: int = 50,
    step_size: float = 0.05,
    k_restraint: float = 1.0,
) -> np.ndarray:
  """Pre-relaxes ligand coordinates in static protein pocket to relieve initial clashes (<2.0 A)."""
  if len(lig_coords) == 0 or len(prot_coords) == 0:
    return lig_coords
  coords = np.array(lig_coords, dtype=float).copy()
  init_coords = coords.copy()

  for _ in range(max_iters):
    forces = np.zeros_like(coords)
    any_clash = False
    for i, lc in enumerate(coords):
      dists = np.linalg.norm(prot_coords - lc, axis=1)
      clash_mask = (dists < clash_cutoff) & (dists > 1e-4)
      if np.any(clash_mask):
        any_clash = True
        diffs = lc - prot_coords[clash_mask]
        d_vals = dists[clash_mask, np.newaxis]
        repulsion = (diffs / d_vals) * (clash_cutoff - d_vals) * 4.0
        forces[i] += np.sum(repulsion, axis=0)

      # Restraint to initial pose (harmonic k=1.0)
      forces[i] -= k_restraint * (lc - init_coords[i])

    if not any_clash:
      break

    coords += step_size * np.clip(forces, -1.0, 1.0)

  return coords


def validate_input_json(json_path: str) -> Dict[str, Any]:
  """Validates AF3 or Boltz API input JSON syntax and counts molecular entities."""
  if not os.path.exists(json_path):
    raise FileNotFoundError(f"Input JSON not found: {json_path}")

  with open(json_path, "r", encoding="utf-8") as f:
    data = json.load(f)

  stem_name = os.path.splitext(os.path.basename(json_path))[0]
  name = (
      data.get("name")
      or data.get("job_id")
      or data.get("job_name")
      or (stem_name if "entities" in data else "unnamed_job")
  )

  if "entities" in data or (
      isinstance(data.get("input"), dict) and "entities" in data["input"]
  ):
    ent_list = (
        data["entities"]
        if "entities" in data
        else data["input"].get("entities", [])
    )
    proteins = [e for e in ent_list if e.get("type") == "protein"]
    ligands = [
        e
        for e in ent_list
        if e.get("type") in ("ligand_smiles", "ligand_ccd", "glycan")
    ]
    dna = [e for e in ent_list if e.get("type") == "dna"]
    rna = [e for e in ent_list if e.get("type") == "rna"]
    return {
        "job_name": name,
        "dialect": "boltz-api",
        "protein_count": len(proteins),
        "ligand_count": len(ligands),
        "dna_count": len(dna),
        "rna_count": len(rna),
        "total_entities": len(ent_list),
    }

  dialect = data.get("dialect", "unknown")
  sequences = data.get("sequences", [])

  ligands = [s for s in sequences if "ligand" in s]
  proteins = [s for s in sequences if "protein" in s]
  dna = [s for s in sequences if "dna" in s]
  rna = [s for s in sequences if "rna" in s]

  if not proteins:
    for k in ["protein", "partner_protein", "proteins"]:
      if k in data:
        val = data[k]
        if isinstance(val, list):
          proteins.extend([{"protein": p} for p in val])
        elif isinstance(val, dict):
          proteins.append({"protein": val})

  if not ligands:
    for k in ["ligand", "ligands"]:
      if k in data:
        val = data[k]
        if isinstance(val, list):
          ligands.extend([{"ligand": l} for l in val])
        elif isinstance(val, dict):
          ligands.append({"ligand": val})

  total_ent = (
      len(sequences)
      if sequences
      else (len(proteins) + len(ligands) + len(dna) + len(rna))
  )

  return {
      "job_name": name,
      "dialect": dialect,
      "protein_count": len(proteins),
      "ligand_count": len(ligands),
      "dna_count": len(dna),
      "rna_count": len(rna),
      "total_entities": total_ent,
  }


def build_local_python_command(
    json_path: str,
    output_dir: str,
    af3_dir: Optional[str] = None,
    model_dir: Optional[str] = None,
    db_dir: Optional[str] = None,
    run_data_pipeline: bool = True,
    run_inference: bool = True,
) -> List[str]:
  """Builds a local Python command to run run_alphafold.py."""
  python_bin = sys.executable
  script_path = "run_alphafold.py"
  if af3_dir:
    script_path = os.path.join(af3_dir, "run_alphafold.py")

  cmd = [
      python_bin,
      script_path,
      f"--json_path={os.path.abspath(json_path)}",
      f"--output_dir={os.path.abspath(output_dir)}",
      f"--run_data_pipeline={'true' if run_data_pipeline else 'false'}",
      f"--run_inference={'true' if run_inference else 'false'}",
  ]
  if model_dir:
    cmd.append(f"--model_dir={os.path.abspath(model_dir)}")
  if db_dir:
    cmd.append(f"--db_dir={os.path.abspath(db_dir)}")

  return cmd


def build_container_command(
    container_engine: str,
    json_path: str,
    output_dir: str,
    model_dir: str,
    db_dir: Optional[str] = None,
    image_name: str = "alphafold3",
    use_gpu: bool = True,
    run_data_pipeline: bool = True,
    run_inference: bool = True,
) -> List[str]:
  """Builds a containerized (Docker or Podman) execution command for local laptops/workstations."""
  abs_json = os.path.abspath(json_path)
  abs_input_dir = os.path.dirname(abs_json)
  abs_output_dir = os.path.abspath(output_dir)
  os.makedirs(abs_output_dir, exist_ok=True)
  abs_model_dir = os.path.abspath(model_dir)
  json_filename = os.path.basename(json_path)

  has_nvidia_gpu = shutil.which("nvidia-smi") is not None and use_gpu
  effective_data_pipeline = bool(run_data_pipeline and db_dir)

  cmd = [container_engine, "run", "--rm"]
  if has_nvidia_gpu and container_engine == "docker":
    cmd.extend(["--gpus", "all"])
  elif has_nvidia_gpu and container_engine == "podman":
    cmd.append("--device=nvidia.com/gpu=all")
  else:
    cmd.extend([
        "-e",
        "XLA_PYTHON_CLIENT_PREALLOCATE=false",
        "-e",
        "XLA_FLAGS=--xla_cpu_multi_thread_eigen=true",
    ])

  cmd.extend([
      "--volume",
      f"{abs_input_dir}:/root/af_input",
      "--volume",
      f"{abs_output_dir}:/root/af_output",
      "--volume",
      f"{abs_model_dir}:/root/models",
  ])

  if db_dir:
    abs_db_dir = os.path.abspath(db_dir)
    cmd.extend(["--volume", f"{abs_db_dir}:/root/public_databases"])

  if not has_nvidia_gpu:
    wrapper_src = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "resources",
        "af3_cpu_wrapper.py",
    )
    wrapper_path = os.path.join(abs_input_dir, "_af3_cpu_wrapper.py")
    shutil.copyfile(wrapper_src, wrapper_path)
    af3_args = [
        f"--json_path=/root/af_input/{json_filename}",
        "--output_dir=/root/af_output",
        "--model_dir=/root/models",
        f"--run_data_pipeline={'true' if effective_data_pipeline else 'false'}",
        f"--run_inference={'true' if run_inference else 'false'}",
        "--num_diffusion_samples=1",
    ]
    if db_dir:
      af3_args.append("--db_dir=/root/public_databases")
    cmd.extend([
        image_name,
        "python",
        "/root/af_input/_af3_cpu_wrapper.py",
        *af3_args,
    ])
  else:
    cmd.extend([
        image_name,
        "python",
        "run_alphafold.py",
        f"--json_path=/root/af_input/{json_filename}",
        "--output_dir=/root/af_output",
        "--model_dir=/root/models",
        f"--run_data_pipeline={'true' if effective_data_pipeline else 'false'}",
        f"--run_inference={'true' if run_inference else 'false'}",
    ])
    if db_dir:
      cmd.append("--db_dir=/root/public_databases")

  return cmd


def resolve_af3_model_dir(
    custom_model_dir: Optional[str] = None,
) -> Optional[str]:
  """Resolves AlphaFold 3 model weights directory from explicit flag, env vars, or standard XDG paths."""
  candidates: List[str] = []
  if custom_model_dir:
    candidates.append(custom_model_dir)
  if os.environ.get("AF3_MODEL_DIR"):
    candidates.append(os.environ["AF3_MODEL_DIR"])
  candidates.extend(_DEFAULT_AF3_MODEL_DIRS)
  for d in candidates:
    if os.path.isdir(d) and (
        os.path.isfile(os.path.join(d, "af3.bin.zst"))
        or os.path.isfile(os.path.join(d, "af3.bin"))
        or any("af3" in f for f in os.listdir(d))
    ):
      return os.path.abspath(d)

  sys.stderr.write(
      "WARNING: AlphaFold 3 model weights (af3.bin.zst / af3.bin) were not"
      f" found in searched paths: {candidates}.\n"
  )
  return None


def execute_boltz2_cofolding(json_path: str, output_dir: str) -> Dict[str, Any]:
  """Executes Boltz-2 MIT open-source joint protein-ligand diffusion co-folding."""
  os.makedirs(output_dir, exist_ok=True)
  with open(json_path, "r", encoding="utf-8") as f:
    data = json.load(f)

  job_name = data.get("name", "boltz2_job")
  yaml_path = os.path.join(output_dir, f"{job_name}_boltz2.yaml")

  lines = ["version: 1", "sequences:"]
  for seq_entry in data.get("sequences", []):
    if "protein" in seq_entry:
      p = seq_entry["protein"]
      lines.append("  - protein:")
      lines.append(f"      id: {p.get('id', 'A')}")
      lines.append(f"      sequence: {p.get('sequence', '')}")
      if p.get("unpairedMsa") == "":
        lines.append("      msa: empty")
    elif "ligand" in seq_entry:
      l = seq_entry["ligand"]
      lines.append("  - ligand:")
      lines.append(f"      id: {l.get('id', 'B')}")
      if "smiles" in l:
        lines.append(f"      smiles: '{l['smiles']}'")
      elif "ccdCodes" in l and l["ccdCodes"]:
        lines.append(f"      ccd: {l['ccdCodes'][0]}")
  with open(yaml_path, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")

  local_boltz = os.path.expanduser("~/.local/bin/boltz")
  boltz_bin = shutil.which("boltz") or (
      local_boltz if os.path.isfile(local_boltz) else None
  )
  uvx_bin = shutil.which("uvx")
  if boltz_bin:
    cmd = [boltz_bin, "predict", yaml_path, "--out_dir", output_dir]
  elif uvx_bin:
    cmd = [
        uvx_bin,
        "--python",
        "3.11",
        "--from",
        "boltz",
        "boltz",
        "predict",
        yaml_path,
        "--out_dir",
        output_dir,
    ]
  else:
    raise RuntimeError(
        "Boltz-2 binary ('boltz') and 'uvx' were not found on PATH. Install"
        " the `boltz` package or `uv` (`uvx --python 3.11 --from boltz`)."
        " CRITICAL: Do NOT fake co-folding or substitute template alignment."
    )

  env = os.environ.copy()
  is_apple_silicon = platform.system() == "Darwin" and platform.machine() in (
      "arm64",
      "aarch64",
  )
  if shutil.which("nvidia-smi") is not None:
    accelerator = "gpu"
  elif is_apple_silicon:
    accelerator = "mps"
    env["PYTORCH_ENABLE_MPS_FALLBACK"] = "1"
    cmd.extend([
        "--accelerator",
        "mps",
        "--recycling_steps",
        "1",
        "--sampling_steps",
        "50",
        "--diffusion_samples",
        "1",
    ])
  else:
    accelerator = "cpu"
    n_threads = str(min(8, max(1, os.cpu_count() or 4)))
    env["OMP_NUM_THREADS"] = n_threads
    env["MKL_NUM_THREADS"] = n_threads
    cmd.extend([
        "--accelerator",
        "cpu",
        "--no_kernels",
        "--recycling_steps",
        "1",
        "--sampling_steps",
        "50",
        "--diffusion_samples",
        "1",
    ])
  print(
      "Launching Boltz-2 (MIT Open-Source Diffusion,"
      f" accelerator={accelerator}): {' '.join(cmd)}"
  )
  subprocess.run(cmd, check=True, env=env)

  metadata_path = os.path.join(output_dir, "inference_metadata.json")
  meta = {
      "job_name": job_name,
      "engine_used": "boltz2",
      "accelerator": accelerator,
      "model_type": "BOLTZ2_DIFFUSION",
      "boltz_yaml_path": yaml_path,
  }
  with open(metadata_path, "w", encoding="utf-8") as f:
    json.dump(meta, f, indent=2)
  return meta


def convert_af3_or_boltz_to_boltz_api_input(
    data: Dict[str, Any],
    msa_free: bool = False,
    num_samples: int = 1,
) -> Dict[str, Any]:
  """Converts an AlphaFold 3 JSON dict or native Boltz API dict into a boltz-2.1 prediction_input dict."""
  if "entities" in data:
    payload = dict(data)
    if num_samples > 1 and "num_samples" not in payload:
      payload["num_samples"] = int(num_samples)
    return payload
  if isinstance(data.get("input"), dict) and "entities" in data["input"]:
    payload = dict(data["input"])
    if num_samples > 1 and "num_samples" not in payload:
      payload["num_samples"] = int(num_samples)
    return payload

  entities: List[Dict[str, Any]] = []
  polymer_chain_ids = set()
  smiles_chain_ids = set()
  ligand_chain_ids: List[str] = []
  protein_chain_ids: List[str] = []
  has_nucleic = False

  for seq_entry in data.get("sequences", []):
    if "protein" in seq_entry:
      p = seq_entry["protein"]
      cid_raw = p.get("id", "A")
      cids = cid_raw if isinstance(cid_raw, list) else [str(cid_raw)]
      for c in cids:
        polymer_chain_ids.add(c)
        protein_chain_ids.append(c)
      ent: Dict[str, Any] = {
          "type": "protein",
          "value": p.get("sequence", ""),
          "chain_ids": cids,
      }
      if msa_free or ("unpairedMsa" in p and p.get("pairedMsa") == ""):
        ent["msa"] = {"type": "empty"}
      entities.append(ent)
    elif "ligand" in seq_entry:
      lig = seq_entry["ligand"]
      cid_raw = lig.get("id", "B")
      cids = cid_raw if isinstance(cid_raw, list) else [str(cid_raw)]
      for c in cids:
        ligand_chain_ids.append(c)
      if "smiles" in lig:
        for c in cids:
          smiles_chain_ids.add(c)
        entities.append({
            "type": "ligand_smiles",
            "value": lig["smiles"],
            "chain_ids": cids,
        })
      elif "ccdCodes" in lig and lig["ccdCodes"]:
        entities.append({
            "type": "ligand_ccd",
            "value": lig["ccdCodes"][0],
            "chain_ids": cids,
        })
    elif "dna" in seq_entry:
      d = seq_entry["dna"]
      cid_raw = d.get("id", "C")
      cids = cid_raw if isinstance(cid_raw, list) else [str(cid_raw)]
      for c in cids:
        polymer_chain_ids.add(c)
      has_nucleic = True
      entities.append({
          "type": "dna",
          "value": d.get("sequence", ""),
          "chain_ids": cids,
      })
    elif "rna" in seq_entry:
      r = seq_entry["rna"]
      cid_raw = r.get("id", "D")
      cids = cid_raw if isinstance(cid_raw, list) else [str(cid_raw)]
      for c in cids:
        polymer_chain_ids.add(c)
      has_nucleic = True
      entities.append({
          "type": "rna",
          "value": r.get("sequence", ""),
          "chain_ids": cids,
      })

  payload: Dict[str, Any] = {"entities": entities}
  if protein_chain_ids and ligand_chain_ids and not has_nucleic:
    payload["binding"] = {
        "type": "ligand_protein_binding",
        "binder_chain_id": ligand_chain_ids[0],
    }

  bonded_pairs = data.get("bondedAtomPairs", [])
  if bonded_pairs:
    bonds = []
    for pair in bonded_pairs:
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


def build_boltz_api_command(
    payload_json_path: str,
    job_name: str,
    model: str = "boltz-2.1",
) -> List[str]:
  """Builds the boltz-api CLI command for structure-and-binding predictions."""
  return [
      "boltz-api",
      "predictions:structure-and-binding",
      "run",
      "--model",
      model,
      "--name",
      job_name,
      "--input",
      f"@json://{os.path.abspath(payload_json_path)}",
  ]


def _normalize_boltz_api_outputs(
    output_dir: str,
    job_name: str,
    model: str,
    run_data: Dict[str, Any],
    source_cif_path: Optional[str] = None,
    cif_text: Optional[str] = None,
) -> Dict[str, Any]:
  """Normalizes Boltz API outputs into <job_name>_model.cif, <job_name>_confidences.json, and inference_metadata.json."""
  dest_cif_path = os.path.join(output_dir, f"{job_name}_model.cif")
  if source_cif_path and os.path.isfile(source_cif_path):
    if os.path.abspath(source_cif_path) != os.path.abspath(dest_cif_path):
      shutil.copyfile(source_cif_path, dest_cif_path)
  elif cif_text:
    with open(dest_cif_path, "w", encoding="utf-8") as f:
      f.write(cif_text)

  pred_id = run_data.get("id")
  out_block = run_data.get("output") or {}
  best_sample = out_block.get("best_sample") or {}
  metrics = best_sample.get("metrics") or {}
  binding_metrics = out_block.get("binding_metrics") or {}

  confidences: Dict[str, Any] = {
      "model": run_data.get("model", model),
      "prediction_id": pred_id,
      "livemode": run_data.get("livemode"),
      "ptm": metrics.get("ptm"),
      "iptm": metrics.get("iptm"),
      "ligand_iptm": metrics.get("ligand_iptm"),
      "protein_iptm": metrics.get("protein_iptm"),
      "structure_confidence": metrics.get("structure_confidence"),
      "complex_plddt": metrics.get("complex_plddt"),
      "complex_iplddt": metrics.get("complex_iplddt"),
      "complex_pde": metrics.get("complex_pde"),
      "complex_ipde": metrics.get("complex_ipde"),
      "ranking_score": (
          metrics.get("structure_confidence")
          if metrics.get("structure_confidence") is not None
          else metrics.get("iptm")
      ),
      "binding_metrics": binding_metrics,
      "binding_confidence": binding_metrics.get("binding_confidence"),
      "optimization_score": binding_metrics.get("optimization_score"),
  }
  confidences_path = os.path.join(output_dir, f"{job_name}_confidences.json")
  with open(confidences_path, "w", encoding="utf-8") as f:
    json.dump(confidences, f, indent=2)

  metadata_path = os.path.join(output_dir, "inference_metadata.json")
  meta = {
      "job_name": job_name,
      "engine_used": "boltz-api",
      "model_type": model,
      "prediction_id": pred_id,
      "cif_path": dest_cif_path if os.path.isfile(dest_cif_path) else None,
      "confidences_path": confidences_path,
      "metrics": metrics,
      "binding_metrics": binding_metrics,
  }
  with open(metadata_path, "w", encoding="utf-8") as f:
    json.dump(meta, f, indent=2)
  return meta


def execute_boltz_api_cofolding(
    json_path: str,
    output_dir: str,
    model: str = "boltz-2.1",
    msa_free: bool = False,
    num_samples: int = 1,
    dry_run: bool = False,
) -> Dict[str, Any]:
  """Executes structure and binding predictions via the hosted Boltz API (boltz-2.1)."""
  os.makedirs(output_dir, exist_ok=True)
  with open(json_path, "r", encoding="utf-8") as f:
    raw_data = json.load(f)

  stem_name = os.path.splitext(os.path.basename(json_path))[0]
  job_name = (
      raw_data.get("name")
      or raw_data.get("job_id")
      or raw_data.get("job_name")
      or stem_name
  )
  prediction_input = convert_af3_or_boltz_to_boltz_api_input(
      raw_data, msa_free=msa_free, num_samples=num_samples
  )

  payload_path = os.path.join(output_dir, f"{job_name}_boltz_api_input.json")
  with open(payload_path, "w", encoding="utf-8") as f:
    json.dump(prediction_input, f, indent=2)

  cli_cmd = build_boltz_api_command(payload_path, job_name, model=model)
  print(f"\nExecution Command (boltz-api):\n{' '.join(cli_cmd)}\n")
  if dry_run:
    print("Dry-run mode: Exiting without launching Boltz API job.")
    return {
        "job_name": job_name,
        "engine_used": "boltz-api",
        "model_type": model,
        "boltz_api_input_path": payload_path,
        "dry_run": True,
    }

  api_key = os.environ.get("BOLTZ_API_KEY", "").strip()
  base_url = os.environ.get("BOLTZ_BASE_URL", "https://api.boltz.bio").strip()
  cli_bin = shutil.which("boltz-api")
  sdk_spec = importlib.util.find_spec("boltz_api")

  if not api_key and not cli_bin:
    raise RuntimeError(
        "BOLTZ_API_KEY is not set in the environment or ~/.env. Configure"
        " BOLTZ_API_KEY (and optional BOLTZ_BASE_URL, default"
        " https://api.boltz.bio) in ~/.env via the credentials skill before"
        " running `--engine boltz-api`."
    )

  # Path 1: Official boltz_api Python SDK
  if sdk_spec is not None and api_key:
    boltz_mod = importlib.import_module("boltz_api")
    client = boltz_mod.Boltz(base_url=base_url, api_key=api_key)
    print(f"Submitting {model} prediction '{job_name}' via boltz_api SDK...")
    prediction = client.predictions.structure_and_binding.start(
        model=model,
        input=prediction_input,
    )
    print(
        f"Prediction started (id={prediction.id}). Polling and downloading"
        f" results to {output_dir}..."
    )
    run_dir = str(
        client.experiments.download_results(
            id=prediction.id,
            name=job_name,
            root_dir=output_dir,
        )
    )
    run_json_path = os.path.join(run_dir, "run.json")
    run_data: Dict[str, Any] = {"id": prediction.id, "model": model}
    if os.path.isfile(run_json_path):
      with open(run_json_path, "r", encoding="utf-8") as f:
        run_data = json.load(f)
    cif_candidates = sorted(
        glob.glob(os.path.join(run_dir, "**", "*.cif"), recursive=True)
    )
    return _normalize_boltz_api_outputs(
        output_dir=output_dir,
        job_name=job_name,
        model=model,
        run_data=run_data,
        source_cif_path=cif_candidates[0] if cif_candidates else None,
    )

  # Path 2: Official boltz-api CLI (supports API key or OAuth session)
  if cli_bin:
    env = os.environ.copy()
    env["BOLTZ_BASE_URL"] = base_url
    start_cmd = [
        cli_bin,
        "--format",
        "raw",
        "predictions:structure-and-binding",
        "start",
        "--model",
        model,
        "--input",
        f"@json://{os.path.abspath(payload_path)}",
    ]
    start_res = subprocess.run(
        start_cmd, capture_output=True, text=True, check=True, env=env
    )
    start_obj = json.loads(start_res.stdout)
    pred_id = start_obj["id"]
    print(
        f"Prediction started via boltz-api CLI (id={pred_id}). Downloading"
        " results..."
    )
    dl_cmd = [
        cli_bin,
        "download-results",
        "--id",
        pred_id,
        "--name",
        job_name,
    ]
    subprocess.run(dl_cmd, cwd=output_dir, check=True, env=env)
    run_dir = os.path.join(output_dir, "boltz-experiments", job_name)
    run_json_path = os.path.join(run_dir, "run.json")
    run_data = start_obj
    if os.path.isfile(run_json_path):
      with open(run_json_path, "r", encoding="utf-8") as f:
        run_data = json.load(f)
    cif_candidates = sorted(
        glob.glob(os.path.join(run_dir, "**", "*.cif"), recursive=True)
    )
    return _normalize_boltz_api_outputs(
        output_dir=output_dir,
        job_name=job_name,
        model=model,
        run_data=run_data,
        source_cif_path=cif_candidates[0] if cif_candidates else None,
    )

  # Path 3: Direct REST via polite_http.http_client.HttpClient
  api_client = http_client.HttpClient(
      base_url,
      qps=2.0,
      max_retries=3,
      timeout=30.0,
  )
  headers = {
      "x-api-key": api_key,
      "Content-Type": "application/json",
  }
  start_url = (
      f"{base_url.rstrip('/')}/compute/v1/predictions/structure-and-binding"
  )
  print(f"Submitting {model} prediction '{job_name}' to {start_url}...")
  start_resp = api_client.fetch_json(
      start_url,
      method="POST",
      headers=headers,
      json_body={"model": model, "input": prediction_input},
  )
  pred_id = start_resp["id"]
  status_url = f"{start_url}/{pred_id}"
  print(f"Prediction started (id={pred_id}). Polling {status_url}...")

  run_data = start_resp
  for _ in range(180):
    status = run_data.get("status")
    if status == "succeeded":
      break
    if status == "failed":
      err_info = run_data.get("error") or {}
      raise RuntimeError(
          f"Boltz API prediction {pred_id} failed:"
          f" {err_info.get('code', 'ERROR')} -"
          f" {err_info.get('message', run_data)}"
      )
    time.sleep(5.0)
    run_data = api_client.fetch_json(status_url, method="GET", headers=headers)
  else:
    raise RuntimeError(
        f"Boltz API prediction {pred_id} timed out waiting for completion."
    )

  run_json_path = os.path.join(output_dir, "run.json")
  with open(run_json_path, "w", encoding="utf-8") as f:
    json.dump(run_data, f, indent=2)

  cif_text = None
  out_block = run_data.get("output") or {}
  best_sample = out_block.get("best_sample") or {}
  struct_info = best_sample.get("structure") or {}
  struct_url = struct_info.get("url")
  if struct_url:
    dl_client = http_client.HttpClient(struct_url, qps=5.0, timeout=60.0)
    cif_text = dl_client.fetch_text(struct_url)

  return _normalize_boltz_api_outputs(
      output_dir=output_dir,
      job_name=job_name,
      model=model,
      run_data=run_data,
      cif_text=cif_text,
  )


def execute_esmfold2_cofolding(
    json_path: str, output_dir: str, fast_mode: bool = False
) -> Dict[str, Any]:
  """Executes ESMFold2 all-atom protein-ligand complex co-folding via OpenProtein API."""
  if openprotein is None or op_molecules is None:
    raise ImportError(
        "openprotein-python is required for --engine esmfold2; install"
        " openprotein-python or use --engine af3, boltz2, or boltz-api."
    )

  os.makedirs(output_dir, exist_ok=True)
  with open(json_path, "r", encoding="utf-8") as f:
    data = json.load(f)

  job_name = data.get("name", "esmfold2_job")
  try:
    session = openprotein.connect()
    chains: Dict[str, Any] = {}
    for seq_entry in data.get("sequences", []):
      if "protein" in seq_entry:
        p = seq_entry["protein"]
        prot = op_molecules.Protein(sequence=p.get("sequence", ""))
        if fast_mode:
          prot = prot.single_sequence_mode()
        chains[p.get("id", "A")] = prot
      elif "ligand" in seq_entry:
        l = seq_entry["ligand"]
        if "smiles" in l:
          chains[l.get("id", "L")] = op_molecules.Ligand(smiles=l["smiles"])
        elif "ccdCodes" in l and l["ccdCodes"]:
          chains[l.get("id", "L")] = op_molecules.Ligand(ccd=l["ccdCodes"][0])
      elif "dna" in seq_entry:
        d = seq_entry["dna"]
        chains[d.get("id", "D")] = op_molecules.DNA(d.get("sequence", ""))
      elif "rna" in seq_entry:
        r = seq_entry["rna"]
        chains[r.get("id", "R")] = op_molecules.RNA(r.get("sequence", ""))

    complex_obj = op_molecules.Complex(chains=chains)
    model = session.fold.esmfold2_fast if fast_mode else session.fold.esmfold2
    print(
        f"Launching OpenProtein {'ESMFold2-Fast' if fast_mode else 'ESMFold2'}"
        " complex co-folding..."
    )
    future = model.fold([complex_obj], num_steps=50, seed=42)
    future.wait_until_done(verbose=True, timeout=900)
    results = future.get()
    pred = results[0]

    cif_path = os.path.join(output_dir, f"{job_name}_esmfold2.cif")
    with open(cif_path, "w", encoding="utf-8") as f:
      f.write(pred.to_cif())

    meta = {
        "job_name": job_name,
        "engine_used": "esmfold2_fast" if fast_mode else "esmfold2",
        "model_type": "ESMFOLD2_OPENPROTEIN",
        "cif_path": cif_path,
    }
    metadata_path = os.path.join(output_dir, "inference_metadata.json")
    with open(metadata_path, "w", encoding="utf-8") as f:
      json.dump(meta, f, indent=2)
    return meta
  except Exception as e:
    raise RuntimeError(
        f"ESMFold2 execution failed ({e}). Ensure valid OpenProtein API"
        " credentials (OPENPROTEIN_USERNAME / OPENPROTEIN_API_KEY in ~/.env)."
    ) from e


def main() -> None:
  parser = argparse.ArgumentParser(
      description=(
          "Multi-Engine Biomolecular Co-Folding Runner (AlphaFold 3 default,"
          " Boltz-2, Boltz API, ESMFold2)."
      )
  )
  parser.add_argument(
      "--json-path",
      "--input",
      dest="json_path",
      required=True,
      help="Path to input co-folding JSON (AlphaFold 3 or Boltz API format).",
  )
  parser.add_argument(
      "--output-dir", required=True, help="Directory to save predicted outputs."
  )
  parser.add_argument(
      "--engine",
      choices=["af3", "boltz2", "boltz-api", "esmfold2", "auto"],
      default="af3",
      help="Co-folding model engine to execute (default: af3).",
  )
  parser.add_argument(
      "--boltz-model",
      default="boltz-2.1",
      help="Boltz API model version (default: boltz-2.1).",
  )
  parser.add_argument(
      "--num-samples",
      type=int,
      default=1,
      help="Number of structure samples for Boltz API (default: 1).",
  )
  parser.add_argument(
      "--runner",
      choices=["auto", "python", "docker", "podman"],
      default="auto",
      help=(
          "Container/runtime execution mode for AF3 (default: auto-detect"
          " podman/docker)."
      ),
  )
  parser.add_argument(
      "--model-dir",
      default=None,
      help="Directory containing AF3 weights (af3.bin.zst).",
  )
  parser.add_argument(
      "--db-dir",
      default=None,
      help=(
          "Directory containing reference genetic databases (optional for"
          " MSA-free)."
      ),
  )
  parser.add_argument(
      "--af3-dir",
      default=None,
      help="Path to cloned alphafold3 repository if using --runner=python.",
  )
  parser.add_argument(
      "--image",
      default="alphafold3",
      help="Docker/Podman image name (default: 'alphafold3').",
  )
  parser.add_argument(
      "--no-gpu",
      action="store_true",
      help="Disable GPU acceleration (CPU-only mode).",
  )
  parser.add_argument(
      "--skip-data-pipeline",
      action="store_true",
      help="Skip data pipeline (MSA search / single-sequence mode).",
  )
  parser.add_argument(
      "--skip-inference",
      action="store_true",
      help="Skip model inference (data pipeline featurization only).",
  )
  parser.add_argument(
      "--dry-run", action="store_true", help="Print command without executing."
  )
  parser.add_argument(
      "--accept-af3-terms",
      action="store_true",
      help=(
          "Acknowledge Google DeepMind AlphaFold 3 Model Parameters Terms of"
          " Use (Non-Commercial Use Only)."
      ),
  )

  args = parser.parse_args()

  try:
    val = validate_input_json(args.json_path)
    print(
        f"Validated co-folding input: '{val['job_name']}' (Proteins:"
        f" {val['protein_count']}, Ligands: {val['ligand_count']})"
    )

    if args.engine == "boltz-api":
      execute_boltz_api_cofolding(
          json_path=args.json_path,
          output_dir=args.output_dir,
          model=args.boltz_model,
          msa_free=args.skip_data_pipeline,
          num_samples=args.num_samples,
          dry_run=args.dry_run,
      )
      return

    resolved_model_dir = resolve_af3_model_dir(args.model_dir)
    active_runner = args.runner
    if active_runner == "auto":
      for cand in ("podman", "docker"):
        if shutil.which(cand):
          active_runner = cand
          break
      else:
        active_runner = "docker"

    if args.engine == "boltz2":
      execute_boltz2_cofolding(args.json_path, args.output_dir)
      return

    if args.engine == "esmfold2":
      execute_esmfold2_cofolding(
          args.json_path, args.output_dir, fast_mode=args.skip_data_pipeline
      )
      return

    local_boltz = os.path.expanduser("~/.local/bin/boltz")
    boltz_bin = shutil.which("boltz") or (
        local_boltz if os.path.isfile(local_boltz) else None
    )
    if args.engine == "auto" and not resolved_model_dir and not args.dry_run:
      if boltz_bin or shutil.which("uvx"):
        print(
            "NOTICE: AlphaFold 3 weights not found in default paths;"
            " --engine=auto selected Boltz-2 (MIT open-source diffusion)..."
        )
        execute_boltz2_cofolding(args.json_path, args.output_dir)
        return
      if os.environ.get("BOLTZ_API_KEY"):
        print(
            "NOTICE: AlphaFold 3 weights not found in default paths;"
            " --engine=auto selected hosted Boltz API (boltz-2.1)..."
        )
        execute_boltz_api_cofolding(
            json_path=args.json_path,
            output_dir=args.output_dir,
            model=args.boltz_model,
            msa_free=args.skip_data_pipeline,
            num_samples=args.num_samples,
            dry_run=args.dry_run,
        )
        return

    if not resolved_model_dir:
      err_msg = (
          "AlphaFold 3 model weights (af3.bin.zst) were not found in default"
          " paths.\nPlease prompt the user to choose one of the following"
          " options:\n  1. Provide a local path to DeepMind AlphaFold 3 weights"
          " via `--model-dir /path/to/af3_weights` (or `export"
          " AF3_MODEL_DIR=...`) under the Google DeepMind AlphaFold 3 Model"
          " Parameters Terms of Use (Non-Commercial Use Only).\n  2. Run"
          " MIT-licensed open-source co-folding instead by passing `--engine"
          " boltz2` (weights auto-cached in ~/.boltz/).\n  3. Run hosted"
          " Boltz-2.1 structure and binding API by passing `--engine"
          " boltz-api` (requires BOLTZ_API_KEY in ~/.env).\n  4. Run"
          " OpenProtein ESMFold2 by passing `--engine esmfold2`."
      )
      raise RuntimeError(err_msg)

    license_file = os.path.join(
        ".licenses", "cheminformatics_cofolding_LICENSE.txt"
    )
    terms_accepted = args.accept_af3_terms or os.path.isfile(license_file)
    if not terms_accepted and not args.dry_run:
      print(
          "\n[NOTICE] Using Google DeepMind AlphaFold 3 Model Parameters"
          f" ({resolved_model_dir}).\nAlphaFold 3 model parameters are subject"
          " to the AlphaFold 3 Model Parameters Terms of Use (Non-Commercial"
          " Research Only; https://github.com/google-deepmind/alphafold3).\nFor"
          " commercial or unrestricted open-source use, re-run with `--engine"
          " boltz2` (MIT License) or `--engine boltz-api`.\n"
      )
      if sys.stdin.isatty():
        resp = (
            input(
                "Confirm use of DeepMind AlphaFold 3 weights under"
                " Non-Commercial Terms of Use? [y/N]: "
            )
            .strip()
            .lower()
        )
        if resp not in ("y", "yes"):
          print(
              "Aborted AF3 run. Re-run with `--engine boltz2` or `--engine"
              " boltz-api`."
          )
          sys.exit(1)
      else:
        raise RuntimeError(
            "Using DeepMind AlphaFold 3 weights requires user confirmation of"
            " the AlphaFold 3 Model Parameters Terms of Use (Non-Commercial Use"
            " Only). Prompt the user to confirm AF3 non-commercial use (then"
            " pass `--accept-af3-terms` or create"
            " `.licenses/cheminformatics_cofolding_LICENSE.txt`), or run"
            " `--engine boltz2` (MIT License) / `--engine boltz-api` instead."
        )

    if active_runner == "python":
      cmd = build_local_python_command(
          json_path=args.json_path,
          output_dir=args.output_dir,
          af3_dir=args.af3_dir,
          model_dir=resolved_model_dir,
          db_dir=args.db_dir,
          run_data_pipeline=not args.skip_data_pipeline,
          run_inference=not args.skip_inference,
      )
    else:
      cmd = build_container_command(
          container_engine=active_runner,
          json_path=args.json_path,
          output_dir=args.output_dir,
          model_dir=resolved_model_dir,
          db_dir=args.db_dir,
          image_name=args.image,
          use_gpu=not args.no_gpu,
          run_data_pipeline=not args.skip_data_pipeline,
          run_inference=not args.skip_inference,
      )

    cmd_str = " ".join(cmd)
    print(f"\nExecution Command ({active_runner}):\n{cmd_str}\n")

    if args.dry_run:
      print("Dry-run mode: Exiting without launching command.")
      return

    print(f"Launching AlphaFold 3 via {active_runner}...")
    subprocess.run(cmd, check=True)

  except Exception as e:
    print(f"Error executing co-folding runner: {e}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
  main()
