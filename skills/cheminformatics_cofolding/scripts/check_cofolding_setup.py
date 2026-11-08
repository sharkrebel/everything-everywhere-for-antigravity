# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "boltz-api",
#     "openprotein-python",
#     "polite-http",
#     "python-dotenv",
#     "rdkit",
# ]
# ///

"""Multi-Engine Biomolecular Co-Folding Environment Diagnostic Utility.

Checks system readiness across all supported co-folding engines:
1. AlphaFold 3 (af3 - default): Container runtime (podman/docker), image, model
   parameters via --model-dir / $AF3_MODEL_DIR, and GPU acceleration.
2. Boltz-2 (boltz2): Open-source MIT biomolecular diffusion engine (CLI/Python).
3. Boltz API (boltz-api): Hosted Boltz-2.1 cloud API for structure & binding
   affinity predictions (binding_confidence and optimization_score).
4. ESMFold2 (esmfold2): All-atom complex co-folding via OpenProtein API.
"""

import argparse
import importlib.util
import json
import os
import platform
import shutil
import subprocess
from typing import Any, Dict, List, Optional

import dotenv

dotenv.load_dotenv(os.path.expanduser("~/.env"))

_DEFAULT_AF3_MODEL_DIRS = (
    os.path.expanduser("~/.cache/alphafold3"),
    os.path.expanduser("~/.local/share/alphafold3"),
    "./models",
)


def check_container_engine() -> Dict[str, Any]:
  """Detects Podman or Docker container runtimes."""
  for engine in ("podman", "docker"):
    path = shutil.which(engine)
    if path:
      try:
        res = subprocess.run(
            [engine, "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if res.returncode == 0:
          return {
              "available": True,
              "engine": engine,
              "version": res.stdout.strip(),
          }
      except Exception:
        pass
  return {"available": False, "engine": None, "version": None}


def check_af3_weights(custom_model_dir: Optional[str] = None) -> Dict[str, Any]:
  """Checks for AlphaFold 3 model weights via explicit flags, env vars, or standard XDG paths."""
  candidates: List[str] = []
  if custom_model_dir:
    candidates.append(custom_model_dir)
  if os.environ.get("AF3_MODEL_DIR"):
    candidates.append(os.environ["AF3_MODEL_DIR"])
  candidates.extend(_DEFAULT_AF3_MODEL_DIRS)

  for d in candidates:
    for fname in ("af3.bin.zst", "af3.bin"):
      p = os.path.join(d, fname)
      if os.path.isfile(p):
        return {
            "found": True,
            "path": os.path.abspath(p),
            "directory": os.path.abspath(d),
        }

  return {
      "found": False,
      "searched_paths": candidates,
      "remediation": (
          "Pass --model-dir /path/to/af3_weights or set export"
          " AF3_MODEL_DIR=/path/to/af3_weights. Weights can be requested from"
          " https://github.com/google-deepmind/alphafold3"
      ),
  }


def _get_system_ram_gb() -> float:
  """Returns total physical RAM in GiB (or 16.0 if unavailable)."""
  try:
    pages = os.sysconf("SC_PHYS_PAGES")
    page_size = os.sysconf("SC_PAGE_SIZE")
    if pages > 0 and page_size > 0:
      return round((pages * page_size) / (1024.0**3), 1)
  except (AttributeError, ValueError, OSError):
    pass
  return 16.0


def check_gpu() -> Dict[str, Any]:
  """Checks for NVIDIA GPU acceleration (CUDA), Apple Silicon Metal (MPS), or CPU-only laptop mode."""
  ram_gb = _get_system_ram_gb()
  if shutil.which("nvidia-smi"):
    try:
      res = subprocess.run(
          [
              "nvidia-smi",
              "--query-gpu=name,memory.total",
              "--format=csv,noheader",
          ],
          capture_output=True,
          text=True,
          timeout=5,
      )
      if res.returncode == 0 and res.stdout.strip():
        return {
            "available": True,
            "accelerator": "gpu",
            "device": res.stdout.strip().splitlines()[0],
            "ram_gb": ram_gb,
        }
    except Exception:
      pass

  if platform.system() == "Darwin" and platform.machine() in (
      "arm64",
      "aarch64",
  ):
    return {
        "available": True,
        "accelerator": "mps",
        "device": (
            f"Apple Silicon Metal (MPS, {ram_gb} GB unified RAM; laptop flags"
            " enabled)"
        ),
        "ram_gb": ram_gb,
    }

  return {
      "available": False,
      "accelerator": "cpu",
      "device": (
          f"CPU-only laptop mode ({ram_gb} GB RAM; JAX/Boltz CPU flags enabled"
          " automatically)"
      ),
      "ram_gb": ram_gb,
  }


def check_boltz2() -> Dict[str, Any]:
  """Checks whether Boltz-2 (MIT open-source biomolecular diffusion) is installed."""
  local_boltz = os.path.expanduser("~/.local/bin/boltz")
  cli_path = shutil.which("boltz") or (
      local_boltz if os.path.isfile(local_boltz) else None
  )
  spec = importlib.util.find_spec("boltz")
  if cli_path or spec is not None:
    return {
        "available": True,
        "cli_path": cli_path,
        "license": "MIT (100% open-source weights auto-cached in ~/.boltz/)",
    }
  return {
      "available": False,
      "remediation": (
          "Install the `boltz` package (or ensure `uvx` is on PATH for"
          " on-demand `uvx --python 3.11 --from boltz` execution) for"
          " MIT-licensed joint protein-ligand diffusion co-folding."
      ),
  }


def check_boltz_api() -> Dict[str, Any]:
  """Checks whether hosted Boltz API (boltz-2.1) credentials and client are configured."""
  api_key = os.environ.get("BOLTZ_API_KEY", "").strip()
  base_url = os.environ.get("BOLTZ_BASE_URL", "https://api.boltz.bio").strip()
  cli_path = shutil.which("boltz-api")
  sdk_installed = importlib.util.find_spec("boltz_api") is not None
  http_fallback = importlib.util.find_spec("polite_http") is not None

  key_mode = None
  if api_key:
    if "_test_" in api_key:
      key_mode = "test"
    elif "_live_" in api_key:
      key_mode = "live"
    else:
      key_mode = "configured"

  if api_key and (sdk_installed or cli_path or http_fallback):
    backend_parts = []
    if sdk_installed:
      backend_parts.append("boltz_api Python SDK")
    if cli_path:
      backend_parts.append(f"CLI ({cli_path})")
    if http_fallback and not (sdk_installed or cli_path):
      backend_parts.append("polite_http REST client")
    return {
        "available": True,
        "base_url": base_url,
        "key_mode": key_mode,
        "cli_path": cli_path,
        "sdk_installed": sdk_installed,
        "backend": ", ".join(backend_parts),
        "model": "boltz-2.1",
        "notes": (
            "Hosted Boltz-2.1 structure & binding affinity predictions"
            " (binding_confidence and optimization_score)."
        ),
    }

  return {
      "available": False,
      "base_url": base_url,
      "key_mode": key_mode,
      "cli_path": cli_path,
      "sdk_installed": sdk_installed,
      "remediation": (
          "Configure BOLTZ_API_KEY (and optional BOLTZ_BASE_URL, default"
          " https://api.boltz.bio) in ~/.env via the credentials skill to"
          " enable hosted Boltz-2.1 structure and binding predictions."
      ),
  }


def check_esmfold2() -> Dict[str, Any]:
  """Checks whether OpenProtein Python SDK (esmfold2 / esmfold2_fast) is configured for all-atom complex co-folding."""
  has_creds = bool(
      os.environ.get("OPENPROTEIN_USERNAME")
      or os.environ.get("OPENPROTEIN_API_KEY")
  )
  spec = importlib.util.find_spec("openprotein")
  if spec is not None and has_creds:
    return {
        "available": True,
        "backend": "openprotein (session.fold.esmfold2 / esmfold2_fast)",
        "credentials_detected": has_creds,
        "notes": (
            "All-atom complex co-folding of Protein + Ligand(smiles=...) +"
            " DNA/RNA (supports MSA and single-sequence fast mode)."
        ),
    }
  if spec is None:
    return {
        "available": False,
        "sdk_installed": False,
        "credentials_detected": has_creds,
        "remediation": (
            "Install `openprotein-python` (declared in PEP 723 dependencies;"
            " run via `uv run scripts/check_cofolding_setup.py`) to enable"
            " ESMFold2 all-atom co-folding."
        ),
    }
  return {
      "available": False,
      "sdk_installed": True,
      "credentials_detected": False,
      "remediation": (
          "Configure OPENPROTEIN_USERNAME / OPENPROTEIN_API_KEY in ~/.env"
          " (via the credentials skill) to enable ESMFold2 all-atom co-folding."
      ),
  }


def main() -> None:
  parser = argparse.ArgumentParser(
      description="Multi-Engine Biomolecular Co-Folding Setup Diagnostic"
  )
  parser.add_argument(
      "--model-dir",
      default=None,
      help="Explicit path to AlphaFold 3 model weights directory",
  )
  parser.add_argument(
      "--json", action="store_true", help="Output diagnostic report as JSON"
  )
  parser.add_argument(
      "-o",
      "--output",
      "--output-json",
      dest="output_json",
      default=None,
      help="Write diagnostic JSON report to file path",
  )
  args = parser.parse_args()

  container = check_container_engine()
  af3_weights = check_af3_weights(args.model_dir)
  gpu = check_gpu()
  boltz2 = check_boltz2()
  boltz_api = check_boltz_api()
  esmfold2 = check_esmfold2()

  af3_ready = bool(container["available"] and af3_weights["found"])
  any_ready = (
      af3_ready
      or boltz2["available"]
      or boltz_api["available"]
      or esmfold2["available"]
  )

  report = {
      "default_engine": "af3",
      "any_engine_ready": any_ready,
      "engines": {
          "af3": {
              "ready": af3_ready,
              "container": container,
              "weights": af3_weights,
              "gpu": gpu,
          },
          "boltz2": boltz2,
          "boltz_api": boltz_api,
          "esmfold2": esmfold2,
      },
  }

  if args.output_json:
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_json)), exist_ok=True
    )
    with open(args.output_json, "w", encoding="utf-8") as f:
      json.dump(report, f, indent=2)

  if args.json:
    print(json.dumps(report, indent=2))
    return

  print(
      "================ Biomolecular Co-Folding Environment Diagnostic"
      " ================"
  )
  print(
      f"[{'✓' if container['available'] else '✗'}] Container Runtime:  "
      f" {container['version'] or 'Not found (podman/docker)'}"
  )
  print(
      f"[{'✓' if gpu['available'] else '!'}] Hardware Accel:     "
      f" {gpu['device']}"
  )
  print(
      "--------------------------------------------------------------------------------"
  )
  print("Engine 1 (Default): AlphaFold 3 (af3)")
  if af3_weights["found"]:
    print(f"  [✓] Model Weights:     Found at {af3_weights['path']}")
  else:
    print(
        "  [✗] Model Weights:     Not found in $AF3_MODEL_DIR or standard XDG"
        " paths"
    )
    print(f"      --> Remediation:   {af3_weights['remediation']}")
  print(
      "--------------------------------------------------------------------------------"
  )
  print("Engine 2 (Backup #1 - Open-Source Local Diffusion): Boltz-2 (boltz2)")
  if boltz2["available"]:
    print(
        "  [✓] Status:            Available"
        f" ({boltz2['cli_path'] or 'Python module'}) - {boltz2['license']}"
    )
  else:
    print(f"  [!] Status:            Not installed ({boltz2['remediation']})")
  print(
      "--------------------------------------------------------------------------------"
  )
  print(
      "Engine 3 (Cloud Structure & Binding API): Boltz API (boltz-api /"
      " boltz-2.1)"
  )
  if boltz_api["available"]:
    print(
        f"  [✓] Status:            Available ({boltz_api['backend']};"
        f" key_mode={boltz_api['key_mode']}, base_url={boltz_api['base_url']})"
    )
  else:
    print(
        f"  [!] Status:            Not configured ({boltz_api['remediation']})"
    )
  print(
      "--------------------------------------------------------------------------------"
  )
  print(
      "Engine 4 (Backup #2 - All-Atom Complex Co-Folding via OpenProtein):"
      " ESMFold2 (esmfold2)"
  )
  if esmfold2["available"]:
    print(f"  [✓] Status:            Available ({esmfold2['backend']})")
  else:
    print(
        f"  [!] Status:            Not configured ({esmfold2['remediation']})"
    )
  print(
      "================================================================================"
  )
  if not any_ready:
    print("Status: No co-folding engine is currently configured on this host.")
    print(
        "Configure AlphaFold 3 ($AF3_MODEL_DIR), Boltz-2, Boltz API"
        " (BOLTZ_API_KEY), or ESMFold2 (OPENPROTEIN_API_KEY) before running"
        " inference."
    )


if __name__ == "__main__":
  main()
