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

"""Production Virtual Screening Engine: High-Throughput Substructure, Similarity, and Physicochemical Filtering."""

import argparse
import gzip
import json
import os
import sys
import time
from typing import Any, Dict, Generator, Iterable, List, Optional, Sequence, Tuple, Union
import urllib.parse
import pandas as pd
from polite_http import http_client
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.Scaffolds import MurckoScaffold

_PUBCHEM_CLIENT = http_client.HttpClient(
    "https://pubchem.ncbi.nlm.nih.gov/rest/", qps=5.0
)


def stream_molecules(
    file_path: str,
) -> Generator[Tuple[str, Chem.Mol, Optional[Dict[str, Any]]], None, None]:
  """Streams molecules with constant memory usage from SDF, SDF.GZ, CSV, or SMILES files.

  Yields:
    A (name, mol, record) triple. `record` carries the molecule's source row so
    that screening can return hits with their input data intact; it is None for
    formats where the fields live on the molecule itself (SDF) or where there
    are none (SMILES), in which case the caller reads them from the molecule.
  """
  ext = os.path.splitext(file_path)[1].lower()
  if file_path.endswith(".sdf.gz"):
    with gzip.open(file_path, "rb") as gz:
      suppl = Chem.ForwardSDMolSupplier(gz)
      for idx, mol in enumerate(suppl):
        if mol is not None:
          name = (
              mol.GetProp("_Name") if mol.HasProp("_Name") else f"Mol_{idx+1}"
          )
          yield name, mol, None
  elif ext == ".sdf":
    suppl = Chem.SDMolSupplier(file_path)
    for idx, mol in enumerate(suppl):
      if mol is not None:
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Mol_{idx+1}"
        yield name, mol, None
  elif ext in (".smi", ".smiles", ".txt"):
    with open(file_path, "r", encoding="utf-8") as f:
      for idx, line in enumerate(f):
        line = line.strip()
        if not line or line.startswith("#"):
          continue
        parts = line.split()
        smi = parts[0]
        name = parts[1] if len(parts) > 1 else f"Mol_{idx+1}"
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
          yield name, mol, None
  elif ext == ".csv":
    df_chunks = pd.read_csv(file_path, chunksize=1000, low_memory=False)
    mol_idx = 0
    skipped_empty = 0
    skipped_unparseable = 0
    for chunk in df_chunks:
      # Prefer the canonical column names. Falling straight to "first column
      # whose name contains 'smi'" misfires on wide assay tables, where
      # columns such as Inhibition_Malate_dehydrogenase also match.
      smi_col = next(
          (
              c
              for c in ("canonical_smiles", "smiles", "SMILES", "Smiles")
              if c in chunk.columns
          ),
          None,
      )
      if smi_col is None:
        smi_col = next(
            (c for c in chunk.columns if "smi" in c.lower()), chunk.columns[0]
        )
      id_col = next(
          (
              c
              for c in chunk.columns
              if "id" in c.lower() or "name" in c.lower()
          ),
          None,
      )
      for _, row in chunk.iterrows():
        mol_idx += 1
        raw = row[smi_col]
        name = str(row[id_col]) if id_col else f"Mol_{mol_idx}"
        if pd.isna(raw) or not str(raw).strip():
          # Structureless rows are legitimate (an assay-only patent example,
          # or a Markush record), but they cannot be screened.
          skipped_empty += 1
          continue
        mol = Chem.MolFromSmiles(str(raw))
        if mol is None:
          skipped_unparseable += 1
          continue
        yield name, mol, row.to_dict()
    if skipped_empty or skipped_unparseable:
      print(
          f"[library_search] Read {mol_idx} rows from {smi_col}; skipped"
          f" {skipped_empty} with no structure and {skipped_unparseable}"
          " whose SMILES did not parse.",
          file=sys.stderr,
      )


def search_pubchem_live(
    smarts_query: str, max_records: int = 50
) -> List[Tuple[str, Chem.Mol, Optional[Dict[str, Any]]]]:
  """Executes live substructure search on PubChem returning real commercial candidate molecules."""
  enc_smarts = urllib.parse.quote(smarts_query, safe="")
  mols = []
  try:
    data = _PUBCHEM_CLIENT.fetch_json(
        f"pug/compound/fastsubstructure/smarts/{enc_smarts}/property/CanonicalSMILES,IUPACName,MolecularWeight/JSON?MaxRecords={max_records}",
        timeout=18.0,
    )
    props = data.get("PropertyTable", {}).get("Properties", [])
    for p in props:
      cid = str(p.get("CID", ""))
      smi = p.get("ConnectivitySMILES") or p.get("CanonicalSMILES")
      name = p.get("IUPACName") or f"CID_{cid}"
      if smi:
        mol = Chem.MolFromSmiles(smi)
        if mol:
          mol.SetProp("CID", cid)
          mol.SetProp("IUPACName", name)
          mols.append((f"CID_{cid}", mol, dict(p)))
  except Exception as e:
    sys.stderr.write(f"Live PubChem query error: {e}\n")
  return mols


def run_virtual_screen(
    mol_stream: Iterable[Tuple[Any, ...]],
    query_smarts: Optional[Union[str, Sequence[str]]] = None,
    query_scaffolds: Optional[Union[str, Sequence[str]]] = None,
    query_smiles: Optional[str] = None,
    sim_threshold: float = 0.6,
    max_mw: float = 600.0,
    max_logp: float = 5.5,
    max_tpsa: float = 160.0,
    max_rotb: int = 12,
    filter_pains: bool = False,
    max_hits: int = 0,
) -> Dict[str, Any]:
  """Executes multi-tier virtual screening pipeline with formal attrition statistics.

  Args:
    mol_stream: Tuples of (name, mol) or (name, mol, source_record).
    query_smarts: One SMARTS pattern, or several. Several are combined with OR,
      because a single patent commonly claims more than one chemical series and
      no one core selects them all.
    query_scaffolds: One Murcko scaffold SMILES, or several, as reported by
      survey_scaffolds. A molecule matches when its own Murcko scaffold is
      identical, so the hit count reproduces the survey's member count exactly.
      Combined with query_smarts by OR. Use this for a scaffold taken from the
      survey and query_smarts for a hand-generalized core.
    query_smiles: Lead structure for Tanimoto similarity screening.
    sim_threshold: Tanimoto cutoff applied when query_smiles is given.
    max_mw: Molecular weight ceiling.
    max_logp: cLogP ceiling.
    max_tpsa: TPSA ceiling.
    max_rotb: Rotatable bond ceiling.
    filter_pains: Whether to reject PAINS matches.
    max_hits: Stop after this many hits; 0 or negative means exhaustive.

  Returns:
    A dict with "metrics" and "hits". Each hit carries a Matched_Core column
    naming the pattern(s) it matched, which is what tells a multi-series result
    apart from a single-series one.
  """
  is_unlimited = max_hits is None or max_hits <= 0
  total_screened = 0
  passed_substruct = 0
  passed_similarity = 0
  passed_properties = 0
  passed_pains = 0
  hits = []

  # Compile Substructure Queries. A bare string stays valid so existing callers
  # and the single-core case are unchanged.
  if query_smarts is None:
    smarts_list: List[str] = []
  elif isinstance(query_smarts, str):
    smarts_list = [query_smarts]
  else:
    smarts_list = [s for s in query_smarts if s]

  compiled_queries: List[Tuple[str, Any, Any]] = []
  for smarts in smarts_list:
    q = Chem.MolFromSmarts(smarts)
    if q is None:
      # Silently ignoring an unparseable core would report zero hits and look
      # like a chemistry result rather than a typo.
      raise ValueError(
          f"--query_smarts pattern did not parse as SMARTS: {smarts}"
      )
    compiled_queries.append((smarts, q, Chem.PatternFingerprint(q)))

  # Scaffold queries are matched scaffold-to-scaffold rather than as
  # substructures, so a scaffold copied from the survey selects exactly the
  # molecules the survey counted under it.
  if query_scaffolds is None:
    scaffold_list: List[str] = []
  elif isinstance(query_scaffolds, str):
    scaffold_list = [query_scaffolds]
  else:
    scaffold_list = [s for s in query_scaffolds if s]
  scaffold_lookup: Dict[str, str] = {}
  for scaf in scaffold_list:
    canon = Chem.CanonSmiles(scaf)
    if canon is None:
      raise ValueError(f"--query_scaffold did not parse as SMILES: {scaf}")
    scaffold_lookup[canon] = scaf

  core_hit_counts: Dict[str, int] = {s: 0 for s in smarts_list + scaffold_list}

  # Compile Similarity Query (ECFP4)
  gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
  sim_q_mol = Chem.MolFromSmiles(query_smiles) if query_smiles else None
  sim_q_fp = gen.GetFingerprint(sim_q_mol) if sim_q_mol else None

  # Compile PAINS Catalog
  pains_catalog = None
  if filter_pains:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    pains_catalog = FilterCatalog(params)

  start_time = time.time()

  for item in mol_stream:
    name, mol = item[0], item[1]
    record = item[2] if len(item) > 2 else None
    total_screened += 1

    # 1. Core Filter. SMARTS cores are substructure matches; scaffold cores are
    # compared against this molecule's own Murcko scaffold. Any match passes.
    matched_cores: List[str] = []
    if compiled_queries or scaffold_lookup:
      if compiled_queries:
        mol_pat_fp = Chem.PatternFingerprint(mol)
        for smarts, q_mol, q_pat_fp in compiled_queries:
          if not DataStructs.AllProbeBitsMatch(q_pat_fp, mol_pat_fp):
            continue
          if mol.HasSubstructMatch(q_mol):
            matched_cores.append(smarts)
      if scaffold_lookup:
        try:
          mol_scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        except Exception:
          mol_scaffold = ""
        if mol_scaffold in scaffold_lookup:
          matched_cores.append(scaffold_lookup[mol_scaffold])
      if not matched_cores:
        continue
      for core in matched_cores:
        core_hit_counts[core] += 1
    passed_substruct += 1

    # 2. Similarity Filter
    sim_score = 1.0
    if sim_q_fp is not None:
      mol_fp = gen.GetFingerprint(mol)
      sim_score = DataStructs.TanimotoSimilarity(sim_q_fp, mol_fp)
      if sim_score < sim_threshold:
        continue
    passed_similarity += 1

    # 3. Physicochemical Property Filters
    mw = Descriptors.ExactMolWt(mol)
    clogp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rotb = Lipinski.NumRotatableBonds(mol)

    if mw > max_mw or clogp > max_logp or tpsa > max_tpsa or rotb > max_rotb:
      continue
    passed_properties += 1

    # 4. PAINS Filter
    if pains_catalog is not None:
      if pains_catalog.HasMatch(mol):
        continue
    passed_pains += 1

    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)

    # Carry the input record through. A filter must not destroy the data it
    # filters: assay values, provenance and patent example ids all live in the
    # source row, and dropping them here silently strips the corpus down to
    # structures alone. Screen-computed fields below take precedence on any
    # name collision.
    if record is None:
      record = mol.GetPropsAsDict(includePrivate=False, includeComputed=False)
    hits.append({
        **record,
        "ID": name,
        "SMILES": Chem.MolToSmiles(mol),
        "Matched_Core": ";".join(matched_cores) if matched_cores else "",
        "Similarity_Score": round(sim_score, 3) if sim_q_fp else 1.0,
        "MW": round(mw, 1),
        "cLogP": round(clogp, 2),
        "TPSA": round(tpsa, 1),
        "RotB": rotb,
        "HBD": hbd,
        "HBA": hba,
        "Fsp3": round(fsp3, 2),
    })

    if not is_unlimited and len(hits) >= max_hits:
      break

  elapsed = round(time.time() - start_time, 2)
  hit_rate = round(
      (len(hits) / total_screened * 100) if total_screened > 0 else 0, 2
  )

  return {
      "metrics": {
          "total_screened": total_screened,
          "substructure_passes": passed_substruct,
          "substructure_passes_by_core": core_hit_counts,
          "similarity_passes": passed_similarity,
          "property_passes": passed_properties,
          "pains_passes": passed_pains if filter_pains else "N/A",
          "final_hit_count": len(hits),
          "hit_rate_pct": hit_rate,
          "screen_duration_sec": elapsed,
      },
      "hits": hits,
  }


def survey_scaffolds(
    mol_stream: Iterable[Tuple[Any, ...]],
    min_heavy_atoms: int = 20,
    top_n: int = 25,
) -> Dict[str, Any]:
  """Reports the Murcko scaffolds present in a library, most populated first.

  A patent harvest is a mixed bag: the claimed series, its synthetic
  intermediates, and often a second series entirely. Choosing a core SMARTS
  without looking first means guessing. This reports what is actually there so
  the cores can be read off the data.

  Args:
    mol_stream: Tuples of (name, mol) or (name, mol, source_record).
    min_heavy_atoms: Molecules below this are skipped as reagents and fragments,
      which otherwise dominate the counts with trivial scaffolds.
    top_n: How many scaffolds to report.

  Returns:
    A dict with "total", "skipped_small", "acyclic", "distinct_scaffolds" and
    "scaffolds", a list of per-scaffold summaries.
  """
  counts: Dict[str, int] = {}
  macro: Dict[str, int] = {}
  examples: Dict[str, str] = {}
  total = 0
  skipped_small = 0
  acyclic = 0

  for item in mol_stream:
    name, mol = item[0], item[1]
    total += 1
    if mol.GetNumHeavyAtoms() < min_heavy_atoms:
      skipped_small += 1
      continue
    try:
      scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
      continue
    if not scaffold:
      # Fully acyclic molecules have no Murcko scaffold and cannot anchor a
      # core query at all.
      acyclic += 1
      continue
    counts[scaffold] = counts.get(scaffold, 0) + 1
    ring_sizes = [len(r) for r in mol.GetRingInfo().AtomRings()]
    if ring_sizes and max(ring_sizes) >= 12:
      macro[scaffold] = macro.get(scaffold, 0) + 1
    examples.setdefault(scaffold, str(name))

  ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
  scaffolds = [
      {
          "scaffold_smiles": s,
          "members": n,
          "macrocyclic_members": macro.get(s, 0),
          "example_id": examples.get(s, ""),
      }
      for s, n in ranked[:top_n]
  ]
  return {
      "total": total,
      "skipped_small": skipped_small,
      "acyclic": acyclic,
      "distinct_scaffolds": len(counts),
      "scaffolds": scaffolds,
  }


def main():
  parser = argparse.ArgumentParser(
      description="High-Throughput Virtual Screening Engine."
  )
  parser.add_argument(
      "--input",
      type=str,
      default=None,
      help="Input chemical library (.sdf, .sdf.gz, .csv, .smi).",
  )
  parser.add_argument(
      "--source",
      type=str,
      choices=["file", "pubchem"],
      default="file",
      help="Screening source.",
  )
  parser.add_argument(
      "--query_smarts",
      type=str,
      action="append",
      default=None,
      help=(
          "Substructure SMARTS filter (e.g. 'c[Br,I]'). Repeat the flag to"
          " accept any of several cores, which is what a patent covering more"
          " than one chemical series needs. Hits record which core matched in"
          " a Matched_Core column. Use --scaffold_survey first to find out"
          " which cores are present."
      ),
  )
  parser.add_argument(
      "--query_scaffold",
      type=str,
      action="append",
      default=None,
      help=(
          "Murcko scaffold SMILES, as printed by --scaffold_survey. Matched"
          " scaffold-to-scaffold, so the hit count reproduces the survey's"
          " member count. Repeat for several series. Combined with"
          " --query_smarts by OR."
      ),
  )
  parser.add_argument(
      "--scaffold_survey",
      action="store_true",
      help=(
          "Report the Murcko scaffolds present in the input and exit without"
          " screening. Use this before choosing --query_smarts."
      ),
  )
  parser.add_argument(
      "--survey_min_heavy_atoms",
      type=int,
      default=20,
      help=(
          "Ignore molecules below this size during --scaffold_survey"
          " (default: 20), so reagents and fragments do not dominate."
      ),
  )
  parser.add_argument(
      "--survey_top",
      type=int,
      default=25,
      help="How many scaffolds --scaffold_survey reports (default: 25).",
  )
  parser.add_argument(
      "--query_smiles",
      type=str,
      default=None,
      help="Lead SMILES for similarity screening.",
  )
  parser.add_argument(
      "--threshold",
      type=float,
      default=0.5,
      help="Tanimoto similarity cutoff [0.0-1.0].",
  )
  parser.add_argument(
      "--max_mw",
      type=float,
      default=550.0,
      help="Maximum Molecular Weight cutoff.",
  )
  parser.add_argument(
      "--max_logp", type=float, default=5.0, help="Maximum cLogP cutoff."
  )
  parser.add_argument(
      "--max_tpsa", type=float, default=160.0, help="Maximum TPSA cutoff."
  )
  parser.add_argument(
      "--max_rotb",
      type=int,
      default=12,
      help="Maximum rotatable-bond count.",
  )
  parser.add_argument(
      "--filter_pains",
      action="store_true",
      help="Filter out PAINS assay-interference structures.",
  )
  parser.add_argument(
      "--max_hits",
      type=int,
      default=0,
      help="Maximum hits to retain (0 or negative for exhaustive/all).",
  )
  parser.add_argument(
      "--all",
      action="store_true",
      help="Retain ALL screening hits without limit.",
  )
  parser.add_argument(
      "--output_csv",
      "--output",
      "-o",
      dest="output_csv",
      type=str,
      default=None,
      help="Path to save hit table CSV.",
  )
  parser.add_argument(
      "--json", action="store_true", help="Output raw JSON metrics and hits."
  )
  args = parser.parse_args()

  limit = 0 if args.all else args.max_hits
  is_unlimited = limit is None or limit <= 0

  # 1. Setup Molecule Stream
  if args.source == "pubchem" or (args.input is None and args.query_smarts):
    # PubChem's substructure endpoint takes exactly one pattern, so a
    # multi-core query cannot be pushed down to it.
    if args.query_smarts and len(args.query_smarts) > 1:
      print(
          "Note: PubChem live search accepts a single SMARTS. Using the first"
          f" of {len(args.query_smarts)}; the rest are applied locally to the"
          " fetched records.",
          file=sys.stderr,
      )
    query = args.query_smarts[0] if args.query_smarts else "c[Br,I]"
    target_fetch = 1000 if is_unlimited else max(100, limit * 2)
    print(
        f"Executing Live PubChem Substructure Stream for SMARTS: '{query}'"
        f" (Target fetch: {target_fetch})..."
    )
    mol_list = search_pubchem_live(query, max_records=target_fetch)
    mol_stream = (item for item in mol_list)
  elif args.input and os.path.exists(args.input):
    mol_stream = stream_molecules(args.input)
  else:
    print(
        "Error: Please provide an --input file (SDF, CSV, SMI) or specify"
        " --source pubchem with --query_smarts.",
        file=sys.stderr,
    )
    sys.exit(1)

  # 1b. Scaffold survey. This answers "what is in here" and is deliberately a
  # terminal mode: the point is to choose --query_smarts before filtering, not
  # to filter and inspect at once.
  if args.scaffold_survey:
    survey = survey_scaffolds(
        mol_stream,
        min_heavy_atoms=args.survey_min_heavy_atoms,
        top_n=args.survey_top,
    )
    if args.json:
      print(json.dumps(survey, indent=2))
      return
    print("\n" + "=" * 55)
    print("             MURCKO SCAFFOLD SURVEY")
    print("=" * 55)
    print(f"  Molecules read:            {survey['total']}")
    print(
        f"  Below {args.survey_min_heavy_atoms} heavy atoms (skipped):"
        f" {survey['skipped_small']}"
    )
    print(f"  Acyclic (no scaffold):     {survey['acyclic']}")
    print(f"  Distinct scaffolds:        {survey['distinct_scaffolds']}")
    print("=" * 55)
    if not survey["scaffolds"]:
      print(
          "  No scaffolds found. Try lowering --survey_min_heavy_atoms.",
          file=sys.stderr,
      )
      return
    for s in survey["scaffolds"]:
      macro_note = (
          f", {s['macrocyclic_members']} macrocyclic"
          if s["macrocyclic_members"]
          else ""
      )
      print(f"  {s['members']:5d} members{macro_note}")
      print(f"        --query_scaffold '{s['scaffold_smiles']}'")
    print(
        "\nPassing these back as --query_scaffold matches scaffold to scaffold,"
        " so the hit count reproduces the member count above. That is an exact"
        " framework match: a molecule carrying an extra fused ring has a"
        " different scaffold and will not match. Use --query_smarts with a"
        " hand-written core when the series varies at the framework itself."
    )
    return

  # 2. Run Screen
  results = run_virtual_screen(
      mol_stream,
      query_smarts=args.query_smarts,
      query_scaffolds=args.query_scaffold,
      query_smiles=args.query_smiles,
      sim_threshold=args.threshold,
      max_mw=args.max_mw,
      max_logp=args.max_logp,
      max_tpsa=args.max_tpsa,
      max_rotb=args.max_rotb,
      filter_pains=args.filter_pains,
      max_hits=limit,
  )

  # 3. Output Reporting
  if args.json:
    print(json.dumps(results, indent=2))
    return

  m = results["metrics"]
  print("\n" + "=" * 55)
  print("           VIRTUAL SCREENING ATTRITION REPORT")
  print("=" * 55)
  print(f"  Total Molecules Screened:  {m['total_screened']}")
  print(f"  Passed Substructure:       {m['substructure_passes']}")
  by_core = m.get("substructure_passes_by_core") or {}
  if len(by_core) > 1:
    for core, n in sorted(by_core.items(), key=lambda kv: -kv[1]):
      print(f"        {n:6d}  {core}")
  for core, n in by_core.items():
    if n == 0:
      print(f"  [!] Core matched nothing: {core}")
  print(f"  Passed Similarity Filter:  {m['similarity_passes']}")
  print(f"  Passed Property Filters:   {m['property_passes']}")
  if args.filter_pains:
    print(f"  Passed PAINS Filters:      {m['pains_passes']}")
  print(f"  Final Retained Hits:       {m['final_hit_count']}")
  print(f"  Hit Rate:                  {m['hit_rate_pct']}%")
  print(f"  Screening Time:            {m['screen_duration_sec']}s")
  print("=" * 55)

  if results["hits"]:
    df_hits = pd.DataFrame(results["hits"])
    print("\nTop Retained Hits:")
    print(
        df_hits[["ID", "SMILES", "MW", "cLogP", "TPSA", "RotB"]]
        .head(10)
        .to_string(index=False)
    )

    if args.output_csv:
      df_hits.to_csv(args.output_csv, index=False)
      print(
          f"\nSaved {len(df_hits)} hits to CSV:"
          f" {os.path.abspath(args.output_csv)}"
      )


if __name__ == "__main__":
  main()
