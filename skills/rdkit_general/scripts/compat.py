# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "matplotlib",
#     "seaborn",
# ]
# ///

"""RDKit Runtime Version Detection, Ecosystem Integrations, and Compatibility Shims.

Provides seamless backward and forward compatibility across RDKit versions, plus
graceful fallbacks, environment inspection, and capability detection for
optional companion packages and binaries (AutoDock Vina, Smina, Meeko, Gemmi,
Seaborn, Datamol, Mols2Grid, Open Babel).
"""

import contextlib
import importlib.util
import os
import shutil
import sys
from typing import Any, Dict, Iterator, List, Optional, Tuple

# NOTE: Do NOT import `matplotlib`, `cycler`, or `seaborn` at module load.
# Pure chemistry scripts (`descriptors.py`, `standardize_molecule.py`, etc.)
# run in isolated PEP 723 environments that declare only `rdkit`/`numpy` without
# `matplotlib` or `seaborn` installed, and importing `matplotlib` eagerly adds
# ~400 ms backend/font startup overhead to non-plotting chemistry calls.
# Keep plotting imports deferred inside `_load_matplotlib()` and
# `_try_apply_seaborn_theme()`.
import numpy as np
import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Lipinski
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import rdSynthonSpaceSearch
from rdkit.Chem.MolStandardize import rdMolStandardize

# Ensure Lipinski module has NumHBD and NumHBA aliases across RDKit versions
if not hasattr(Lipinski, "NumHBD"):
  Lipinski.NumHBD = Lipinski.NumHDonors
if not hasattr(Lipinski, "NumHBA"):
  Lipinski.NumHBA = Lipinski.NumHAcceptors

# -------------------------------------------------------------------------
# 1. Parse Version Tuple (Major, Minor, [Patch])
# -------------------------------------------------------------------------
try:
  _raw_parts = rdkit.__version__.split(".")
  RDKIT_VERSION: Tuple[int, int] = (int(_raw_parts[0]), int(_raw_parts[1]))
  RDKIT_VERSION_STR: str = rdkit.__version__
except Exception:
  RDKIT_VERSION = (2020, 0)
  RDKIT_VERSION_STR = "unknown"

# -------------------------------------------------------------------------
# 2. Feature & Module Capability Flags with Direct Import Checks
# -------------------------------------------------------------------------
HAS_FP_GENERATORS: bool = rdFingerprintGenerator is not None
HAS_SUPER_PARENT: bool = hasattr(rdMolStandardize, "SuperParent")

try:
  HAS_MACROCYCLE_ETKDG: bool = hasattr(
      AllChem.ETKDGv3(), "useMacrocycleTorsions"
  )
except Exception:
  HAS_MACROCYCLE_ETKDG = False

HAS_SRETKDG: bool = hasattr(AllChem, "srETKDGv3")
HAS_SYNTHON_SPACE: bool = hasattr(rdSynthonSpaceSearch, "SynthonSpace")
HAS_SUBSTRUCT_MATCH_PARAMS: bool = hasattr(Chem, "SubstructMatchParameters")


def _has_mol_representation_svg() -> bool:
  """Detects RDKit PandasTools SVG support without eagerly importing pandas."""
  try:
    return (
        RDKIT_VERSION >= (2019, 9)
        and importlib.util.find_spec("pandas") is not None
        and importlib.util.find_spec("rdkit.Chem.PandasTools") is not None
    )
  except (ModuleNotFoundError, ValueError):
    return False


HAS_MOL_REPRESENTATION_SVG: bool = _has_mol_representation_svg()

# -------------------------------------------------------------------------
# 3. Optional Ecosystem Companion Packages & Binaries
# -------------------------------------------------------------------------
HAS_SEABORN: bool = importlib.util.find_spec("seaborn") is not None
HAS_MEEKO: bool = importlib.util.find_spec("meeko") is not None
HAS_GEMMI: bool = importlib.util.find_spec("gemmi") is not None
HAS_DATAMOL: bool = importlib.util.find_spec("datamol") is not None
HAS_MOLS2GRID: bool = importlib.util.find_spec("mols2grid") is not None
HAS_OPENBABEL: bool = importlib.util.find_spec("openbabel") is not None
HAS_OPENMM: bool = importlib.util.find_spec("openmm") is not None
HAS_PY3DMOL: bool = importlib.util.find_spec("py3Dmol") is not None
HAS_XGBOOST: bool = importlib.util.find_spec("xgboost") is not None
HAS_LIGHTGBM: bool = importlib.util.find_spec("lightgbm") is not None
HAS_STATSMODELS: bool = importlib.util.find_spec("statsmodels") is not None


def _find_binary(name: str) -> Optional[str]:
  """Checks for binary executable on PATH or in ~/.local/bin."""
  home_local = os.path.expanduser(f"~/.local/bin/{name}")
  candidates = [
      home_local,
      shutil.which(name),
  ]
  for c in candidates:
    if c and os.path.isfile(c) and os.access(c, os.X_OK):
      return c
  return None


HAS_VINA: bool = _find_binary("vina") is not None
HAS_SMINA: bool = _find_binary("smina") is not None

# -------------------------------------------------------------------------
# 4. Dependency Metadata Catalog
# -------------------------------------------------------------------------
DEPENDENCY_CATALOG: Dict[str, Dict[str, Any]] = {
    "vina": {
        "type": "binary",
        "description": (
            "AutoDock Vina 1.2.5 Static Linux Binary (Molecular Docking Engine)"
        ),
        "source_url": (
            "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64"
        ),
        "target_path": os.path.expanduser("~/.local/bin/vina"),
        "check_fn": lambda: _find_binary("vina") is not None,
    },
    "smina": {
        "type": "binary",
        "description": (
            "Smina Static Linux Binary (AutoDock Vina Fork with Custom Scoring)"
        ),
        "source_url": (
            "https://sourceforge.net/projects/smina/files/smina.static/download"
        ),
        "target_path": os.path.expanduser("~/.local/bin/smina"),
        "check_fn": lambda: _find_binary("smina") is not None,
    },
    "meeko": {
        "type": "python_package",
        "package_name": "meeko",
        "description": (
            "Meeko (Preparation of Flexible Ligand PDBQT for AutoDock Vina)"
        ),
        "pip_args": ["meeko"],
        "check_fn": lambda: importlib.util.find_spec("meeko") is not None,
    },
    "seaborn": {
        "type": "python_package",
        "package_name": "seaborn",
        "description": "Seaborn (Statistical Data & Publication Plot Styling)",
        "pip_args": ["seaborn"],
        "check_fn": lambda: importlib.util.find_spec("seaborn") is not None,
    },
    "gemmi": {
        "type": "python_package",
        "package_name": "gemmi",
        "description": (
            "Gemmi (Macromolecular Crystallography & PDB/mmCIF Processing)"
        ),
        "pip_args": ["gemmi"],
        "check_fn": lambda: importlib.util.find_spec("gemmi") is not None,
    },
    "datamol": {
        "type": "python_package",
        "package_name": "datamol",
        "description": (
            "Datamol (High-level Molecular Manipulation & Parallel Processing)"
        ),
        "pip_args": ["datamol"],
        "check_fn": lambda: importlib.util.find_spec("datamol") is not None,
    },
    "mols2grid": {
        "type": "python_package",
        "package_name": "mols2grid",
        "description": (
            "Mols2Grid (Interactive 2D Molecular Grid HTML Generation)"
        ),
        "pip_args": ["mols2grid"],
        "check_fn": lambda: importlib.util.find_spec("mols2grid") is not None,
    },
    "openbabel": {
        "type": "python_package",
        "package_name": "openbabel",
        "description": "Open Babel (Chemical File Format Interconversions)",
        "pip_args": ["openbabel-wheel"],
        "check_fn": lambda: importlib.util.find_spec("openbabel") is not None,
    },
}


# -------------------------------------------------------------------------
# 5. Permission-Gated Dependency Checking & Installation
# -------------------------------------------------------------------------
def check_dependency(name: str) -> bool:
  """Checks if a dependency (binary or python package) is currently installed."""
  key = name.lower()
  if key in DEPENDENCY_CATALOG:
    return DEPENDENCY_CATALOG[key]["check_fn"]()
  # Generic module check
  return importlib.util.find_spec(name) is not None


def get_missing_dependencies(names: List[str]) -> List[str]:
  """Returns a list of dependencies from 'names' that are not yet installed."""
  return [name for name in names if not check_dependency(name)]


def get_dependency_install_command(name: str) -> Optional[str]:
  """Returns the recommended dependency specification string for an agent or user."""
  key = name.lower()
  if key not in DEPENDENCY_CATALOG:
    return f"Declare '{name}' in PEP 723 script dependencies"
  meta = DEPENDENCY_CATALOG[key]
  dep_type = meta.get("type")
  if dep_type == "python_package":
    pip_args = meta.get("pip_args", [key])
    return f"Declare {' '.join(pip_args)} in PEP 723 script dependencies"
  elif dep_type == "binary":
    source_url = meta.get("source_url", "")
    return f"Download {key} binary from {source_url} and place on PATH"
  return None


# -------------------------------------------------------------------------
# 6. Core Info & Plot Styling
# -------------------------------------------------------------------------
def get_version_info() -> Dict[str, Any]:
  """Returns a dictionary summary of RDKit version and active capability flags."""
  return {
      "version_str": RDKIT_VERSION_STR,
      "version_tuple": RDKIT_VERSION,
      "core_capabilities": {
          "has_fp_generators": HAS_FP_GENERATORS,
          "has_super_parent": HAS_SUPER_PARENT,
          "has_macrocycle_etkdg": HAS_MACROCYCLE_ETKDG,
          "has_sretkdg": HAS_SRETKDG,
          "has_synthon_space": HAS_SYNTHON_SPACE,
          "has_substruct_match_params": HAS_SUBSTRUCT_MATCH_PARAMS,
          "has_mol_representation_svg": HAS_MOL_REPRESENTATION_SVG,
      },
      "optional_ecosystem": {
          "has_vina": _find_binary("vina") is not None,
          "has_smina": _find_binary("smina") is not None,
          "has_meeko": importlib.util.find_spec("meeko") is not None,
          "has_gemmi": importlib.util.find_spec("gemmi") is not None,
          "has_seaborn": importlib.util.find_spec("seaborn") is not None,
          "has_datamol": importlib.util.find_spec("datamol") is not None,
          "has_mols2grid": importlib.util.find_spec("mols2grid") is not None,
          "has_openbabel": importlib.util.find_spec("openbabel") is not None,
      },
  }


# -------------------------------------------------------------------------
# Publication Plot Styling & Palette
# -------------------------------------------------------------------------
PUBLICATION_PALETTE: Dict[str, str] = {
    "primary": "#4F46E5",  # Indigo 600
    "primary_light": "#6366F1",  # Indigo 500
    "primary_dark": "#4338CA",  # Indigo 700
    "success": "#10B981",  # Emerald 500
    "success_light": "#34D399",  # Emerald 400
    "warning": "#F59E0B",  # Amber 500
    "danger": "#F43F5E",  # Rose 500
    "purple": "#8B5CF6",  # Violet 500
    "cyan": "#06B6D4",  # Cyan 500
    "muted": "#64748B",  # Slate 500
    "background": "#FFFFFF",
    "card": "#F8FAFC",  # Slate 50
    "border": "#E2E8F0",  # Slate 200
    "foreground": "#0F172A",  # Slate 900
    "muted_foreground": "#64748B",
    "text": "#0F172A",
    "text_muted": "#475569",
    "secondary": "#64748B",
    "accent": "#DB2777",
    "grid": "#E2E8F0",
}

PUBLICATION_CYCLE: List[str] = [
    "#4F46E5",  # Indigo
    "#10B981",  # Emerald
    "#F59E0B",  # Amber
    "#F43F5E",  # Rose
    "#8B5CF6",  # Violet
    "#06B6D4",  # Cyan
    "#64748B",  # Slate
    "#EC4899",  # Pink
]


def apply_publication_spines(ax: Any, grid_axis: Optional[str] = None) -> None:
  """Formats axis spines and optionally configures dashed gridlines.

  Args:
    ax: Matplotlib Axes object to style.
    grid_axis: Axis on which to enable dashed gridlines ('both', 'x', or 'y').
      If None, existing axis gridline state is left unchanged.

  Raises:
    ValueError: If `grid_axis` is not None, 'both', 'x', or 'y'.
  """
  if grid_axis is not None and (
      not isinstance(grid_axis, str) or grid_axis not in ("both", "x", "y")
  ):
    raise ValueError(
        f"grid_axis must be None, 'both', 'x', or 'y'; got {grid_axis!r}."
    )
  if hasattr(ax, "spines"):
    if "top" in ax.spines:
      ax.spines["top"].set_visible(False)
    if "right" in ax.spines:
      ax.spines["right"].set_visible(False)
    if "left" in ax.spines:
      ax.spines["left"].set_color("#CBD5E1")
      ax.spines["left"].set_linewidth(0.9)
    if "bottom" in ax.spines:
      ax.spines["bottom"].set_color("#CBD5E1")
      ax.spines["bottom"].set_linewidth(0.9)
  if grid_axis is not None and hasattr(ax, "grid"):
    ax.grid(False)
    ax.grid(True, axis=grid_axis, linestyle="--", alpha=0.5, color="#E2E8F0")


def _load_matplotlib() -> Tuple[Any, Any, Any]:
  """Lazily imports Matplotlib components when plotting helpers are called."""
  try:
    from cycler import cycler
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    return plt, MaxNLocator, cycler
  except ImportError:
    return None, None, None


def _extract_numeric(data_obj: Any) -> Optional[np.ndarray]:
  """Extracts finite numeric values from an array, Series, or sequence."""
  if data_obj is None:
    return None
  try:
    if (
        isinstance(data_obj, (list, tuple))
        and data_obj
        and hasattr(data_obj[0], "__iter__")
    ):
      arrays = []
      for item in data_obj:
        arr = np.asarray(item, dtype=float)
        arr = arr[~np.isnan(arr)]
        if arr.size > 0:
          arrays.append(arr)
      return np.concatenate(arrays) if arrays else None
    arr = np.asarray(data_obj, dtype=float)
    arr = arr[~np.isnan(arr)]
    return arr if arr.size > 0 else None
  except (TypeError, ValueError):
    return None


def autodetect_integer_axes(
    ax: Any,
    x_data: Any = None,
    y_data: Any = None,
    is_histogram: bool = False,
) -> None:
  """Automatically enforces integer tick locators and non-negative bounds when plotted data are integers.

  Args:
    ax: Matplotlib Axes object.
    x_data: Array, Series, or list of arrays/Series plotted on the x-axis.
    y_data: Array, Series, or list of arrays/Series plotted on the y-axis.
    is_histogram: If True, forces integer ticks and bottom=0 on the y-axis
      (count frequency), and if x_data is non-negative integers, sets x-limits
      to [-0.5, max + 0.5] so bars are centered on integer ticks without
      negative tick labels.
  """
  _, max_n_locator, _ = _load_matplotlib()
  if max_n_locator is None or ax is None:
    return

  x_arr = _extract_numeric(x_data)
  if x_arr is not None and x_arr.size > 0:
    if np.all(np.equal(np.mod(x_arr, 1), 0)):
      ax.xaxis.set_major_locator(max_n_locator(integer=True))
      if np.min(x_arr) >= 0:
        if is_histogram:
          ax.set_xlim(left=-0.5, right=float(np.max(x_arr)) + 0.5)
        else:
          cur_xlim = ax.get_xlim()
          if cur_xlim[0] < 0:
            ax.set_xlim(left=-0.5 if np.min(x_arr) == 0 else 0.0)

  if is_histogram:
    ax.yaxis.set_major_locator(max_n_locator(integer=True))
    ax.set_ylim(bottom=0)
  else:
    y_arr = _extract_numeric(y_data)
    if y_arr is not None and y_arr.size > 0:
      if np.all(np.equal(np.mod(y_arr, 1), 0)):
        ax.yaxis.set_major_locator(max_n_locator(integer=True))
        if np.min(y_arr) >= 0:
          cur_ylim = ax.get_ylim()
          if cur_ylim[0] < 0:
            ax.set_ylim(bottom=-0.5 if np.min(y_arr) == 0 else 0.0)


def enforce_discrete_axis_ticks(
    ax: Any, values: Any = None, axis: str = "x"
) -> None:
  """Forces discrete integer tick positions on the requested axis."""
  _, max_n_locator, _ = _load_matplotlib()
  if max_n_locator is None or ax is None:
    return

  if axis in ("x", "both"):
    ax.xaxis.set_major_locator(max_n_locator(integer=True))
    if values is not None:
      autodetect_integer_axes(ax, x_data=values)
  if axis in ("y", "both"):
    ax.yaxis.set_major_locator(max_n_locator(integer=True))
    if values is not None:
      autodetect_integer_axes(ax, y_data=values)


def _try_apply_seaborn_theme(style: str, font_scale: float) -> None:
  """Applies optional Seaborn theme defaults when `seaborn` is installed."""
  try:
    import seaborn as sns
  except ImportError:
    return
  try:
    sns.set_theme(style=style, font_scale=font_scale)
    sns.set_palette(PUBLICATION_CYCLE)
  except (OSError, ValueError):
    pass


def _build_publication_rc_params(
    cycler_cls: Any,
    font_scale: float = 1.0,
    despine: bool = True,
) -> Dict[str, Any]:
  """Builds the dictionary of Matplotlib `rcParams` for publication figures."""
  rc_updates: Dict[str, Any] = {
      "figure.facecolor": "#FFFFFF",
      "axes.facecolor": "#F8FAFC",
      "axes.edgecolor": "#CBD5E1",
      "axes.linewidth": 0.9,
      "axes.grid": True,
      "axes.axisbelow": True,
      "axes.prop_cycle": cycler_cls(color=PUBLICATION_CYCLE),
      "grid.color": "#E2E8F0",
      "grid.linestyle": "--",
      "grid.linewidth": 0.7,
      "grid.alpha": 0.7,
      "font.family": "sans-serif",
      "font.sans-serif": [
          "Google Sans",
          "Inter",
          "system-ui",
          "-apple-system",
          "BlinkMacSystemFont",
          "Segoe UI",
          "Roboto",
          "DejaVu Sans",
          "Helvetica",
          "Arial",
          "sans-serif",
      ],
      "axes.labelcolor": "#334155",
      "axes.labelsize": 10.0 * font_scale,
      "axes.labelweight": "bold",
      "axes.titlesize": 12.0 * font_scale,
      "axes.titleweight": "bold",
      "axes.titlecolor": "#0F172A",
      "xtick.color": "#475569",
      "ytick.color": "#475569",
      "xtick.labelsize": 9.0 * font_scale,
      "ytick.labelsize": 9.0 * font_scale,
      "text.color": "#0F172A",
      "legend.frameon": True,
      "legend.framealpha": 0.95,
      "legend.facecolor": "#FFFFFF",
      "legend.edgecolor": "#CBD5E1",
      "legend.fontsize": 9.0 * font_scale,
      "legend.title_fontsize": 9.5 * font_scale,
      "savefig.bbox": "tight",
      "savefig.dpi": 300,
  }
  if despine:
    rc_updates["axes.spines.top"] = False
    rc_updates["axes.spines.right"] = False
  return rc_updates


@contextlib.contextmanager
def publication_style_context(
    font_scale: float = 1.0,
    despine: bool = True,
) -> Iterator[None]:
  """Applies publication `rcParams` locally and restores prior state on exit."""
  plt, _, cycler_cls = _load_matplotlib()
  if plt is None:
    yield
    return
  rc_updates = _build_publication_rc_params(
      cycler_cls, font_scale=font_scale, despine=despine
  )
  with plt.rc_context(rc=rc_updates):
    yield


def set_safe_plot_style(
    style: str = "whitegrid",
    font_scale: float = 1.0,
    theme: str = "publication",
    despine: bool = True,
) -> None:
  """Applies publication plot styling to process-global Matplotlib state.

  Mutates process-wide `matplotlib.rcParams` (and Seaborn theme state when
  installed). Intended for standalone CLI scripts initializing process styling;
  library callers that must preserve caller `rcParams` across chained calls
  should use `with publication_style_context():` instead.
  """
  del theme
  if not isinstance(style, str):
    raise TypeError(
        f"style must be a str (pass font_scale=... by keyword); got {style!r}."
    )
  plt, _, cycler_cls = _load_matplotlib()
  if plt is None:
    return
  _try_apply_seaborn_theme(style=style, font_scale=font_scale)
  plt.rcParams.update(
      _build_publication_rc_params(
          cycler_cls, font_scale=font_scale, despine=despine
      )
  )


def get_morgan_fingerprint(
    mol: Chem.Mol,
    radius: int = 2,
    n_bits: int = 2048,
    use_features: bool = False,
    use_chirality: bool = False,
    use_counts: bool = False,
    as_numpy: bool = False,
) -> Any:
  """Unified Morgan / Circular Fingerprint generator with automatic fallback."""
  if HAS_FP_GENERATORS:
    if use_features:
      gen = rdFingerprintGenerator.GetMorganGenerator(
          radius=radius,
          fpSize=n_bits,
          includeChirality=use_chirality,
          atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen(),
      )
    else:
      gen = rdFingerprintGenerator.GetMorganGenerator(
          radius=radius,
          fpSize=n_bits,
          includeChirality=use_chirality,
      )
    if as_numpy:
      if use_counts:
        return gen.GetCountFingerprintAsNumPy(mol)
      return gen.GetFingerprintAsNumPy(mol)
    if use_counts:
      return gen.GetCountFingerprint(mol)
    return gen.GetFingerprint(mol)
  else:
    if use_counts:
      fp = AllChem.GetHashedMorganFingerprint(
          mol,
          radius=radius,
          nBits=n_bits,
          useFeatures=use_features,
          useChirality=use_chirality,
      )
    else:
      fp = AllChem.GetMorganFingerprintAsBitVect(
          mol,
          radius=radius,
          nBits=n_bits,
          useFeatures=use_features,
          useChirality=use_chirality,
      )
    if as_numpy:
      arr = np.zeros((n_bits,), dtype=np.int32 if use_counts else np.int8)
      DataStructs.ConvertToNumpyArray(fp, arr)
      return arr
    return fp


def get_etkdg_params(
    version: int = 3,
    use_macrocycles: bool = False,
    strained_rings: bool = False,
    random_seed: int = 42,
    num_threads: int = 0,
    prune_rms_thresh: float = -1.0,
    use_small_ring_torsions: bool = True,
    enforce_chirality: bool = True,
) -> AllChem.EmbedParameters:
  """Unified ETKDG 3D Conformer Generation Parameters."""
  if strained_rings and HAS_SRETKDG:
    params = AllChem.srETKDGv3()
  elif version >= 3 and hasattr(AllChem, "ETKDGv3"):
    params = AllChem.ETKDGv3()
    if use_macrocycles and HAS_MACROCYCLE_ETKDG:
      params.useMacrocycleTorsions = True
  elif version >= 2 and hasattr(AllChem, "ETKDGv2"):
    params = AllChem.ETKDGv2()
  else:
    params = AllChem.ETKDG()

  params.randomSeed = random_seed
  if hasattr(params, "numThreads"):
    params.numThreads = num_threads
  if prune_rms_thresh > 0:
    params.pruneRmsThresh = prune_rms_thresh
  if hasattr(params, "useSmallRingTorsions"):
    params.useSmallRingTorsions = use_small_ring_torsions
  if hasattr(params, "enforceChirality"):
    params.enforceChirality = enforce_chirality
  return params


def standardize_molecule(
    mol: Chem.Mol,
    use_super_parent: bool = False,
    largest_fragment: bool = True,
    uncharge: bool = True,
    canonical_tautomer: bool = True,
    disconnect_metals: bool = True,
) -> Optional[Chem.Mol]:
  """Applies modern RDKit standardization pipeline to a molecule."""
  if mol is None:
    return None
  if not isinstance(mol, Chem.Mol):
    raise TypeError(
        "standardize_molecule expects an rdkit.Chem.Mol instance or None, got"
        f" {type(mol).__name__}."
    )

  try:
    if disconnect_metals:
      mol = rdMolStandardize.DisconnectOrganometallics(mol)

    if use_super_parent and HAS_SUPER_PARENT:
      return rdMolStandardize.SuperParent(mol)

    # 1. Basic cleanup (standardize functional groups, aromaticity)
    mol = rdMolStandardize.Cleanup(mol)

    # 2. Extract largest fragment (strips salts, solvents)
    if largest_fragment:
      mol = rdMolStandardize.FragmentParent(mol)

    # 3. Uncharge / neutralize
    if uncharge:
      uncharger = rdMolStandardize.Uncharger()
      mol = uncharger.uncharge(mol)

    # 4. Canonical tautomer
    if canonical_tautomer:
      te = rdMolStandardize.TautomerEnumerator()
      mol = te.Canonicalize(mol)

    Chem.SanitizeMol(mol)
    return mol
  except Exception as e:
    sys.stderr.write(f"Standardization error: {e}\n")
    return None


def safe_substruct_match(
    mol: Chem.Mol,
    query: Chem.Mol,
    use_chirality: bool = True,
    max_matches: int = 1000,
) -> Tuple[Tuple[int, ...], ...]:
  """Unified substructure matching supporting SubstructMatchParameters or legacy args.

  Args:
    mol: molecule to search.
    query: query as an RDKit Mol. Pass ``Chem.MolFromSmarts(...)`` output, not a
      raw SMARTS string.
    use_chirality: whether matching respects stereochemistry.
    max_matches: cap on the number of matches returned.

  Returns:
    Tuple of atom-index tuples, one per match.

  Raises:
    ValueError: if `query` is None or is not a Mol. A malformed query is not
      silently reported as "no match", because that would turn a typo into a
      confident negative result.
  """
  if query is None:
    raise ValueError(
        "safe_substruct_match received query=None. Chem.MolFromSmarts / "
        "MolFromSmiles returns None on a malformed pattern; check the "
        "pattern rather than treating this as 'no match'."
    )
  if not isinstance(query, Chem.Mol):
    raise ValueError(
        "safe_substruct_match expects an RDKit Mol query, got "
        f"{type(query).__name__}. Convert it first, e.g. "
        "Chem.MolFromSmarts(pattern)."
    )
  if HAS_SUBSTRUCT_MATCH_PARAMS:
    params = Chem.SubstructMatchParameters()
    params.useChirality = use_chirality
    params.maxMatches = max_matches
    return mol.GetSubstructMatches(query, params)
  else:
    return mol.GetSubstructMatches(
        query, useChirality=use_chirality, maxMatches=max_matches
    )


def save_publication_figure(fig, output_path: str, dpi: int = 300) -> str:
  """Saves a Matplotlib figure at high DPI with an opaque white background and tight bounding box."""
  os.makedirs(
      os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True
  )
  fig.savefig(
      output_path,
      dpi=dpi,
      bbox_inches="tight",
      facecolor="white",
      transparent=False,
  )
  return os.path.abspath(output_path)


def print_environment_report():
  """Prints a diagnostic inspection of the active RDKit environment and companion packages."""
  info = get_version_info()
  print(f"RDKit Version: {info['version_str']}")
  print("\nCore RDKit Capabilities:")
  for feat, status in info["core_capabilities"].items():
    status_str = "SUPPORTED (Modern)" if status else "FALLBACK (Legacy/Missing)"
    print(f"  {feat:<28}: {status_str}")

  print("\nOptional Companion Packages & Binaries:")
  for pkg, status in info["optional_ecosystem"].items():
    cmd = get_dependency_install_command(pkg)
    install_hint = f" (Install: `{cmd}`)" if cmd and not status else ""
    status_str = (
        "INSTALLED (Ready)" if status else f"NOT INSTALLED{install_hint}"
    )
    print(f"  {pkg:<28}: {status_str}")


if __name__ == "__main__":
  print_environment_report()
