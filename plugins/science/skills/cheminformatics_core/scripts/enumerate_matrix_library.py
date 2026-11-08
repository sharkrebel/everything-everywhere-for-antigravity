# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "pandas",
#     "scipy",
#     "scikit-learn",
#     "matplotlib",
#     "polite-http",
# ]
# ///

"""Combinatorial Matrix Library Enumerator (M x N Reactant Matrix, Parallel Property Profiling & Multi-Format Export).

Generates complete combinatorial virtual libraries by reacting Pool A (M
reactants) with Pool B (N reactants) across named reaction templates (Amide,
Suzuki, Reductive Amination, SNAr, Sulfonamide, Urea, Esterification) or custom
reaction SMARTS. Supports direct PubChem commercial building block discovery,
CAS verification, vendor catalog matching, C++ FilterCatalog (PAINS/Brenk)
screening, Lipinski Rule of 5, Pfizer CNS MPO, and exports to CSV, SDF, and
interactive HTML grids.
"""

import argparse
import concurrent.futures
import importlib.util
import itertools
import os
import sys
from typing import Any, Dict, List, Optional
import urllib.parse

import pandas as pd
from polite_http import http_client
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

_PUBCHEM_CLIENT = http_client.HttpClient(
    "https://pubchem.ncbi.nlm.nih.gov/rest/", qps=5.0
)

# Canonical CNS MPO implementation (single source of truth).
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

from skill_paths import add_skill_scripts, skill_dir

add_skill_scripts("cheminformatics_core")

from descriptors import calculate_cns_mpo
from draw_molecules import compute_clean_2d_coords
from synthetic_accessibility_scorer import calculate_sa_score

COMMERCIAL_VENDORS = [
    "Enamine",
    "Chemspace",
    "MolPort",
    "Sigma-Aldrich",
    "WuXi AppTec",
    "Selleck Chemicals",
    "MedChemExpress (MCE)",
    "Tocris Bioscience",
    "ChemDiv",
    "Life Chemicals",
    "Specs",
    "TargetMol",
    "BOC Sciences",
    "Combi-Blocks",
    "Fluorochem",
    "Oakwood Chemical",
    "Thermo Fisher Scientific",
    "Alfa Aesar",
    "Ambeed",
    "BLDpharm",
]

NAMED_REACTIONS = {
    "amide": {
        "name": "Amide Coupling (Acid + Amine)",
        "smarts": "[C:1](=O)[OH].[N;!H0:2]>>[C:1](=O)[N:2]",
        "desc": "Carboxylic acid + Primary/Secondary amine -> Amide",
        "default_query_a": "C(=O)O",
        "default_query_b": "[N;H2,H1;!$(NC=O);!$(NS=O)]",
    },
    "suzuki": {
        "name": "Suzuki-Miyaura Coupling (Halide + Boronic Acid/Ester)",
        "smarts": "[c:1][Br,I,Cl].[c,a:2][B,b]([OH,O])[OH,O]>>[c:1]-[c,a:2]",
        "desc": "Aryl/heteroaryl halide + Boronic acid/pinacol ester -> Biaryl",
        "default_query_a": "[c,a][Br,I]",
        "default_query_b": "B(O)O",
    },
    "reductive_amination": {
        "name": "Reductive Amination (Carbonyl + Amine)",
        "smarts": "[C:1]=O.[N;!H0;!$(NC=O):2]>>[C:1][N:2]",
        "desc": "Aldehyde/Ketone + Primary/Secondary amine -> Alkylamine",
        "default_query_a": "[CX3H1,CX3H2]=O",
        "default_query_b": "[N;H2,H1;!$(NC=O)]",
    },
    "snar": {
        "name": "SNAr Nucleophilic Aromatic Substitution",
        "smarts": "[c:1][F,Cl,Br,I].[N,O,S;!H0:2]>>[c:1][N,O,S:2]",
        "desc": "Aryl halide + Amine/Alcohol/Thiol nucleophile",
        "default_query_a": "c[F,Cl,Br]",
        "default_query_b": "[N,O;H1,H2]",
    },
    "sulfonamide": {
        "name": "Sulfonamide Formation (Sulfonyl Halide + Amine)",
        "smarts": "[S:1](=O)(=O)[Cl,OH].[N;!H0;!$(NC=O):2]>>[S:1](=O)(=O)[N:2]",
        "desc": "Sulfonyl chloride/acid + Amine -> Sulfonamide",
        "default_query_a": "S(=O)(=O)Cl",
        "default_query_b": "[N;H2,H1;!$(NC=O)]",
    },
    "urea": {
        "name": "Urea / Carbamate Formation",
        "smarts": "[C:1](=O)[Cl,Oc1ccccc1].[N;!H0:2]>>[C:1](=O)[N:2]",
        "desc": "Chloroformate/Phenyl carbamate + Amine",
        "default_query_a": "C(=O)Cl",
        "default_query_b": "[N;H2,H1]",
    },
    "esterification": {
        "name": "Esterification (Acid + Alcohol)",
        "smarts": "[C:1](=O)[OH].[O;!H0:2]>>[C:1](=O)[O:2]",
        "desc": "Carboxylic acid + Alcohol -> Ester",
        "default_query_a": "C(=O)O",
        "default_query_b": "[OX2H1]",
    },
}


def _resolve_grid_viewer():
  """Locates the authoritative interactive_grid_viewer.py in rdkit_general.

  Returns:
    Absolute path to the viewer script, or None if it cannot be located.
  """
  candidate = os.path.join(
      skill_dir("rdkit_general", "scripts"), "interactive_grid_viewer.py"
  )
  return candidate if os.path.exists(candidate) else None


def is_valid_cas(cas_str: str) -> bool:
  """Validates format and Mod-10 checksum of a CAS Registry Number."""
  if not isinstance(cas_str, str):
    return False
  parts = cas_str.strip().split("-")
  if len(parts) != 3:
    return False
  p1, p2, p3 = parts
  if not (
      p1.isdigit()
      and p2.isdigit()
      and p3.isdigit()
      and 2 <= len(p1) <= 7
      and len(p2) == 2
      and len(p3) == 1
  ):
    return False
  digits = p1 + p2
  return sum(
      int(d) * (len(digits) - i) for i, d in enumerate(digits)
  ) % 10 == int(p3)


def get_pains_filter_catalog() -> FilterCatalog:
  """Initializes RDKit C++ FilterCatalog with PAINS (A/B/C) and Brenk alerts."""
  params = FilterCatalogParams()
  params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
  params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
  params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
  params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
  return FilterCatalog(params)


def query_pubchem_reactants(
    search_smarts: str,
    max_records: int = 500,
) -> List[Dict[str, Any]]:
  """Directly queries PubChem database for commercial compounds with vendor and CAS resolution."""
  enc = urllib.parse.quote(search_smarts, safe="")
  cids = []
  try:
    data = _PUBCHEM_CLIENT.fetch_json(
        f"pug/compound/fastsubstructure/smarts/{enc}/cids/JSON?MaxRecords={max_records}",
        timeout=16.0,
    )
    cids = data.get("IdentifierList", {}).get("CID", [])
  except Exception:
    pass

  if not cids:
    return []

  cids = list(dict.fromkeys(cids))[:max_records]
  chunks = [cids[i : i + 100] for i in range(0, len(cids), 100)]

  def process_chunk(chunk):
    cid_str = ",".join(str(c) for c in chunk)
    chunk_data = {
        c: {"cid": c, "smiles": None, "name": None, "cas": None, "vendors": []}
        for c in chunk
    }

    # 1. Properties
    try:
      d = _PUBCHEM_CLIENT.fetch_json(
          f"pug/compound/cid/{cid_str}/property/ConnectivitySMILES,CanonicalSMILES,IUPACName/JSON",
          timeout=16.0,
      )
      for p in d.get("PropertyTable", {}).get("Properties", []):
        c = p.get("CID")
        if c in chunk_data:
          chunk_data[c]["smiles"] = p.get("ConnectivitySMILES") or p.get(
              "CanonicalSMILES"
          )
          chunk_data[c]["name"] = p.get("IUPACName") or f"CID_{c}"
    except Exception:
      pass

    # 2. CAS Synonyms
    try:
      d = _PUBCHEM_CLIENT.fetch_json(
          f"pug/compound/cid/{cid_str}/synonyms/JSON",
          timeout=16.0,
      )
      for item in d.get("InformationList", {}).get("Information", []):
        c = item.get("CID")
        syns = item.get("Synonym", [])
        cas_list = [s for s in syns if is_valid_cas(s)]
        if c in chunk_data and cas_list:
          chunk_data[c]["cas"] = cas_list[0]
    except Exception:
      pass

    # 3. Vendors
    try:
      d = _PUBCHEM_CLIENT.fetch_json(
          f"pug/compound/cid/{cid_str}/xrefs/SourceName/JSON",
          timeout=16.0,
      )
      for item in d.get("InformationList", {}).get("Information", []):
        c = item.get("CID")
        matched = set()
        for s in item.get("SourceName", []):
          for vendor in COMMERCIAL_VENDORS:
            if vendor.lower() in s.lower():
              matched.add(vendor)
        if c in chunk_data:
          chunk_data[c]["vendors"] = sorted(list(matched))
    except Exception:
      pass

    return [v for v in chunk_data.values() if v["smiles"]]

  with concurrent.futures.ThreadPoolExecutor(
      max_workers=min(4, len(chunks) or 1)
  ) as ex:
    results = list(ex.map(process_chunk, chunks))

  return [item for sublist in results for item in sublist]


def parse_reactant_pool(
    source: Optional[str],
    query_smarts: Optional[str] = None,
    max_count: int = 50,
    pool_label: str = "Pool",
) -> List[Dict[str, Any]]:
  """Loads reactant pool from comma-separated SMILES, CSV file, or queries PubChem."""
  reactants = []

  if source:
    # Check if source is a file path
    if os.path.exists(source):
      if source.endswith(".csv"):
        df = pd.read_csv(source)
        smi_candidates = [
            c
            for c in df.columns
            if "smi" in c.lower() or "structure" in c.lower()
        ]
        if not smi_candidates:
          raise ValueError(
              f"Input CSV {source} for {pool_label} has no SMILES/structure "
              f"column (found columns: {list(df.columns)})."
          )
        smi_col = smi_candidates[0]
        name_col = [
            c for c in df.columns if "name" in c.lower() or "id" in c.lower()
        ]
        name_col = name_col[0] if name_col else smi_col
        for idx, row in df.iterrows():
          smi = str(row[smi_col]).strip()
          nm = str(row[name_col]).strip()
          cas = str(row.get("CAS_Number", "N/A"))
          vendors = str(row.get("Vendors", ""))
          reactants.append({
              "name": nm,
              "smiles": smi,
              "cas": cas,
              "vendors": (
                  [v.strip() for v in vendors.split(",") if v.strip()]
                  if vendors
                  else []
              ),
          })
      elif source.endswith(".sdf"):
        suppl = Chem.SDMolSupplier(source)
        for idx, mol in enumerate(suppl):
          if mol is not None:
            nm = (
                mol.GetProp("_Name")
                if mol.HasProp("_Name")
                else f"{pool_label}_{idx+1}"
            )
            reactants.append({
                "name": nm,
                "smiles": Chem.MolToSmiles(mol),
                "cas": mol.GetProp("CAS") if mol.HasProp("CAS") else "N/A",
                "vendors": [],
            })
    else:
      # Assume comma-separated list of SMILES or names
      items = [s.strip() for s in source.split(",") if s.strip()]
      for idx, item in enumerate(items):
        if "=" in item and Chem.MolFromSmiles(item) is None:
          # Format: Name=SMILES
          nm, smi = item.split("=", 1)
        else:
          nm = f"{pool_label}_{idx+1}"
          smi = item
        reactants.append({
            "name": nm.strip(),
            "smiles": smi.strip(),
            "cas": "N/A",
            "vendors": [],
        })

  elif query_smarts:
    print(
        f"   Querying PubChem for {pool_label} matching '{query_smarts}' (Max:"
        f" {max_count})..."
    )
    reactants = query_pubchem_reactants(query_smarts, max_records=max_count)
    print(
        f"   Retrieved {len(reactants)} commercial reactants for {pool_label}."
    )

  return reactants


def enumerate_matrix_library(
    pool_a: List[Dict[str, Any]],
    pool_b: List[Dict[str, Any]],
    reaction_smarts: str,
    filter_pains: bool = True,
    filter_ro5: bool = False,
    max_products: Optional[int] = None,
) -> pd.DataFrame:
  """Performs combinatorial Cartesian product (M x N) reaction enumeration with complete C++ profiling."""
  rxn = AllChem.ReactionFromSmarts(reaction_smarts)
  if rxn is None:
    raise ValueError(f"Invalid Reaction SMARTS: {reaction_smarts}")
  if rxn.GetNumReactantTemplates() != 2:
    raise ValueError(
        "Reaction SMARTS must have exactly 2 reactant templates, found"
        f" {rxn.GetNumReactantTemplates()}"
    )

  pains_catalog = get_pains_filter_catalog() if filter_pains else None

  # Pre-parse RDKit mols
  parsed_a = []
  for item in pool_a:
    m = Chem.MolFromSmiles(item["smiles"])
    if m is not None:
      parsed_a.append({**item, "mol": m})

  parsed_b = []
  for item in pool_b:
    m = Chem.MolFromSmiles(item["smiles"])
    if m is not None:
      parsed_b.append({**item, "mol": m})

  if not parsed_a or not parsed_b:
    print(
        "Warning: One or both reactant pools have 0 valid chemical structures.",
        file=sys.stderr,
    )
    return pd.DataFrame()

  total_combinations = len(parsed_a) * len(parsed_b)
  print(
      f"Enumerating {len(parsed_a)} (Pool A) x {len(parsed_b)} (Pool B) ="
      f" {total_combinations} theoretical combinations..."
  )

  records = []
  seen_products = set()

  for item_a, item_b in itertools.product(parsed_a, parsed_b):
    if max_products and len(records) >= max_products:
      break

    mol_a, mol_b = item_a["mol"], item_b["mol"]

    # Attempt reaction (A + B)
    outcomes = rxn.RunReactants((mol_a, mol_b))
    # If no match and reaction might be symmetric or inverted order, try (B + A)
    if not outcomes or not outcomes[0]:
      outcomes = rxn.RunReactants((mol_b, mol_a))

    if not outcomes or not outcomes[0]:
      continue

    prod_mol = outcomes[0][0]
    try:
      Chem.SanitizeMol(prod_mol)
    except Exception:
      continue

    prod_smi = Chem.MolToSmiles(prod_mol)
    if prod_smi in seen_products:
      continue
    seen_products.add(prod_smi)

    # PAINS / Alert screening
    pains_alert = "Clean"
    pains_matches = []
    if pains_catalog:
      entries = pains_catalog.GetMatches(prod_mol)
      if entries:
        pains_matches = [e.GetDescription() for e in entries]
        pains_alert = "; ".join(pains_matches)

    if filter_pains and pains_alert != "Clean":
      continue

    # Rule of 5 evaluation
    mw = round(Descriptors.ExactMolWt(prod_mol), 1)
    clogp = round(Crippen.MolLogP(prod_mol), 2)
    hbd = Lipinski.NumHDonors(prod_mol)
    hba = Lipinski.NumHAcceptors(prod_mol)
    rotb = Lipinski.NumRotatableBonds(prod_mol)
    tpsa = round(rdMolDescriptors.CalcTPSA(prod_mol), 1)
    fsp3 = round(rdMolDescriptors.CalcFractionCSP3(prod_mol), 2)
    qed_val = round(QED.qed(prod_mol), 2)
    sa_val = round(calculate_sa_score(prod_mol), 2)
    cns_mpo = calculate_cns_mpo(prod_mol)

    ro5_violations = 0
    if mw > 500.0:
      ro5_violations += 1
    if clogp > 5.0:
      ro5_violations += 1
    if hbd > 5:
      ro5_violations += 1
    if hba > 10:
      ro5_violations += 1

    if filter_ro5 and ro5_violations > 0:
      continue

    prod_id = f"PROD_{len(records)+1:04d}_{item_a['name']}_{item_b['name']}"

    # Combine vendors and CAS
    combined_vendors = sorted(
        list(set(item_a.get("vendors", []) + item_b.get("vendors", [])))
    )

    records.append({
        "Product_ID": prod_id,
        "Product_SMILES": prod_smi,
        "Reactant_A_Name": item_a["name"],
        "Reactant_A_SMILES": item_a["smiles"],
        "Reactant_A_CAS": item_a.get("cas", "N/A"),
        "Reactant_B_Name": item_b["name"],
        "Reactant_B_SMILES": item_b["smiles"],
        "Reactant_B_CAS": item_b.get("cas", "N/A"),
        "Vendors": ", ".join(combined_vendors) if combined_vendors else "N/A",
        "MW": mw,
        "cLogP": clogp,
        "TPSA": tpsa,
        "HBD": hbd,
        "HBA": hba,
        "RotB": rotb,
        "Fsp3": fsp3,
        "QED": qed_val,
        "SA_Score": sa_val,
        "CNS_MPO": cns_mpo,
        "PAINS_Status": pains_alert,
        "Ro5_Violations": ro5_violations,
    })

  return pd.DataFrame(records)


def export_sdf(df: pd.DataFrame, smiles_col: str, output_path: str):
  """Exports DataFrame to an SDFile with full embedded properties."""
  with Chem.SDWriter(output_path) as writer:
    for _, row in df.iterrows():
      smi = str(row[smiles_col])
      mol = Chem.MolFromSmiles(smi)
      if mol is not None:
        compute_clean_2d_coords(mol)
        for col in df.columns:
          if col != smiles_col:
            mol.SetProp(str(col), str(row[col]))
        writer.write(mol)


def dataframe_to_markdown_table(df: pd.DataFrame) -> str:
  """Renders a pandas DataFrame as a Markdown table without external tabulate dependency."""
  if df.empty:
    return ""
  cols = list(df.columns)
  lines = [
      "| " + " | ".join(str(c) for c in cols) + " |",
      "| " + " | ".join("---" for _ in cols) + " |",
  ]
  for _, row in df.iterrows():
    lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
  return "\n".join(lines)


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Combinatorial Matrix Library Enumerator (M x N Reactants, Reaction"
          " SMARTS, C++ Profiling)."
      )
  )
  parser.add_argument(
      "--reaction",
      type=str,
      default="amide",
      choices=list(NAMED_REACTIONS.keys()) + ["custom"],
      help="Named reaction preset or 'custom'.",
  )
  parser.add_argument(
      "--reaction_smarts",
      type=str,
      default=None,
      help="Custom reaction SMARTS string (required if --reaction custom).",
  )
  parser.add_argument(
      "--pool_a",
      type=str,
      default=None,
      help=(
          "Pool A reactants: comma-separated SMILES, 'Name=SMILES', or path to"
          " CSV/SDF."
      ),
  )
  parser.add_argument(
      "--pool_a_query",
      type=str,
      default=None,
      help=(
          "SMARTS pattern to discover Pool A commercial reactants via PubChem."
      ),
  )
  parser.add_argument(
      "--pool_b",
      type=str,
      default=None,
      help=(
          "Pool B reactants: comma-separated SMILES, 'Name=SMILES', or path to"
          " CSV/SDF."
      ),
  )
  parser.add_argument(
      "--pool_b_query",
      type=str,
      default=None,
      help=(
          "SMARTS pattern to discover Pool B commercial reactants via PubChem."
      ),
  )
  parser.add_argument(
      "--max_pool_a",
      type=int,
      default=20,
      help="Max building blocks to query for Pool A (default: 20).",
  )
  parser.add_argument(
      "--max_pool_b",
      type=int,
      default=50,
      help="Max building blocks to query for Pool B (default: 50).",
  )
  parser.add_argument(
      "--max_products",
      type=int,
      default=None,
      help="Maximum total products to generate (default: all).",
  )
  parser.add_argument(
      "--filter_pains",
      action="store_true",
      default=True,
      help="Filter out PAINS / Brenk structural alerts (default: True).",
  )
  parser.add_argument(
      "--no_filter_pains",
      action="store_false",
      dest="filter_pains",
      help="Disable PAINS / Brenk filtering.",
  )
  parser.add_argument(
      "--filter_ro5",
      action="store_true",
      help="Filter out compounds violating Lipinski's Rule of 5.",
  )
  parser.add_argument(
      "--output_csv",
      type=str,
      default=None,
      help="Path to save enumerated library as CSV.",
  )
  parser.add_argument(
      "--output_sdf",
      type=str,
      default=None,
      help="Path to save enumerated library as SDF.",
  )
  parser.add_argument(
      "--output_html",
      type=str,
      default=None,
      help="Path to generate interactive molecular grid HTML viewer.",
  )
  args = parser.parse_args()

  # Determine Reaction SMARTS
  if args.reaction != "custom" and not args.reaction_smarts:
    rxn_info = NAMED_REACTIONS[args.reaction]
    rxn_smarts = rxn_info["smarts"]
    print(f"Reaction Mode: {rxn_info['name']}")
    print(f"Reaction SMARTS: {rxn_smarts}")
    query_a = args.pool_a_query or (
        rxn_info["default_query_a"] if not args.pool_a else None
    )
    query_b = args.pool_b_query or (
        rxn_info["default_query_b"] if not args.pool_b else None
    )
  else:
    if not args.reaction_smarts:
      print(
          "Error: --reaction_smarts is required when using --reaction custom",
          file=sys.stderr,
      )
      sys.exit(1)
    rxn_smarts = args.reaction_smarts
    print(f"Custom Reaction SMARTS: {rxn_smarts}")
    query_a = args.pool_a_query
    query_b = args.pool_b_query

  if not (args.pool_a or query_a) or not (args.pool_b or query_b):
    parser.error(
        "Provide --pool_a (or --pool_a_query / named --reaction) and "
        "--pool_b (or --pool_b_query / named --reaction)."
    )

  # Parse Pool A
  print("\n1. Loading Reactant Pool A...")
  pool_a = parse_reactant_pool(
      args.pool_a,
      query_smarts=query_a,
      max_count=args.max_pool_a,
      pool_label="Pool_A",
  )
  if not pool_a:
    parser.error("Reactant Pool A (--pool_a) contains no valid reactants.")
  print(f"   Pool A contains {len(pool_a)} starting materials.")

  # Parse Pool B
  print("\n2. Loading Reactant Pool B...")
  pool_b = parse_reactant_pool(
      args.pool_b,
      query_smarts=query_b,
      max_count=args.max_pool_b,
      pool_label="Pool_B",
  )
  if not pool_b:
    parser.error("Reactant Pool B (--pool_b) contains no valid reactants.")
  print(f"   Pool B contains {len(pool_b)} building blocks.")

  # Enumerate Cartesian Product
  print("\n3. Performing Combinatorial Reaction Enumeration...")
  df = enumerate_matrix_library(
      pool_a=pool_a,
      pool_b=pool_b,
      reaction_smarts=rxn_smarts,
      filter_pains=args.filter_pains,
      filter_ro5=args.filter_ro5,
      max_products=args.max_products,
  )

  print("\n" + "=" * 85)
  print(f"  ENUMERATED MATRIX LIBRARY: {len(df)} VALIDATED PRODUCTS")
  print("=" * 85)
  if df.empty:
    print("No products formed or all products were eliminated by filters.")
    return

  display_cols = [
      "Product_ID",
      "Reactant_A_Name",
      "Reactant_B_Name",
      "MW",
      "cLogP",
      "TPSA",
      "CNS_MPO",
      "QED",
  ]
  print(dataframe_to_markdown_table(df[display_cols].head(12)))
  if len(df) > 12:
    print(f"\n... and {len(df) - 12} additional products.")

  # Export outputs
  if args.output_csv:
    df.to_csv(args.output_csv, index=False)
    print(
        f"\nSaved CSV ({len(df)} records) to:"
        f" {os.path.abspath(args.output_csv)}"
    )

  if args.output_sdf:
    export_sdf(df, smiles_col="Product_SMILES", output_path=args.output_sdf)
    print(
        f"Saved SDF ({len(df)} records) to: {os.path.abspath(args.output_sdf)}"
    )

  if args.output_html:
    # Import interactive_grid_viewer dynamically
    viewer_path = _resolve_grid_viewer()
    spec = (
        importlib.util.spec_from_file_location(
            "interactive_grid_viewer", viewer_path
        )
        if viewer_path
        else None
    )
    if spec and spec.loader:
      viewer_mod = importlib.util.module_from_spec(spec)
      spec.loader.exec_module(viewer_mod)
      viewer_mod.create_molecular_grid(
          df=df,
          smiles_col="Product_SMILES",
          id_col="Product_ID",
          output_html=args.output_html,
          title=f"Combinatorial Matrix Library ({len(df)} Products)",
      )
    else:
      print(
          f"Warning: Could not import grid viewer from {viewer_path}",
          file=sys.stderr,
      )


if __name__ == "__main__":
  main()
