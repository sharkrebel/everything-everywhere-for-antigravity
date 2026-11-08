# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "rdkit>=2024.3.1",
#   "networkx>=3.0",
# ]
# ///
"""Electronic Lab Notebook (ELN) Reaction Table & ElectronGraph Sanity Checker.

Generates a structured lab-notebook stoichiometry table (MW, equivalents, mmol,
mass in mg/g, optional liquid dispense volume from density or molarity, solvent
volume, theoretical/actual yield, Atom Economy, and PMI) while formally sanity-
checking the reaction stoichiometry, Lewis octets, and formal charges using
ElectronGraph (MolecularGraph & ReactionVerifierEngine).
"""

import argparse
import json
import os
import sys
from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
  sys.path.insert(0, SCRIPT_DIR)

import electron_graph  # noqa: E402
import render_mechanism  # noqa: E402

COMMON_CONDENSATES: list[tuple[str, str, dict[str, int]]] = [
    ("H2O", "O", {"H": 2, "O": 1}),
    ("HCl", "Cl", {"H": 1, "Cl": 1}),
    ("HBr", "Br", {"H": 1, "Br": 1}),
    ("HI", "I", {"H": 1, "I": 1}),
    ("HF", "F", {"H": 1, "F": 1}),
    ("NH3", "N", {"N": 1, "H": 3}),
    ("N2", "N#N", {"N": 2}),
    ("CO2", "O=C=O", {"C": 1, "O": 2}),
    ("CO", "[C-]#[O+]", {"C": 1, "O": 1}),
    ("MeOH", "CO", {"C": 1, "H": 4, "O": 1}),
    ("EtOH", "CCO", {"C": 2, "H": 6, "O": 1}),
    ("AcOH", "CC(=O)O", {"C": 2, "H": 4, "O": 2}),
]

SOLVENT_DENSITIES_G_ML: dict[str, float] = {
    "thf": 0.889,
    "tetrahydrofuran": 0.889,
    "2-methf": 0.855,
    "toluene": 0.867,
    "phme": 0.867,
    " xylene": 0.864,
    "dcm": 1.326,
    "ch2cl2": 1.326,
    "dichloromethane": 1.326,
    "chloroform": 1.489,
    "chcl3": 1.489,
    "dce": 1.253,
    "1,2-dichloroethane": 1.253,
    "dmf": 0.944,
    "dmac": 0.937,
    "dma": 0.937,
    "nmp": 1.028,
    "dmso": 1.100,
    "water": 1.000,
    "h2o": 1.000,
    "meoh": 0.792,
    "methanol": 0.792,
    "etoh": 0.789,
    "ethanol": 0.789,
    "ipa": 0.786,
    "isopropanol": 0.786,
    "2-propanol": 0.786,
    "t-buoh": 0.781,
    "mecn": 0.786,
    "ch3cn": 0.786,
    "acetonitrile": 0.786,
    "etoac": 0.902,
    "ethyl acetate": 0.902,
    "dioxane": 1.033,
    "1,4-dioxane": 1.033,
    "ether": 0.713,
    "et2o": 0.713,
    "diethyl ether": 0.713,
    "mtbe": 0.740,
    "acetone": 0.784,
    "hexane": 0.655,
    "hexanes": 0.655,
    "heptane": 0.684,
    "pentane": 0.626,
    "acoh": 1.049,
    "acetic acid": 1.049,
}


def sanity_check_reaction_with_electron_graph(
    reactant_smiles_list: list[str],
    product_smiles_list: list[str],
    mechanism_spec: dict[str, Any] | None = None,
) -> dict[str, Any]:
  """Audits reactant/product graphs, octets, and atom conservation via ElectronGraph."""
  r_joined = ".".join(s for s in reactant_smiles_list if s)
  p_joined = ".".join(s for s in product_smiles_list if s)

  r_mg = electron_graph.MolecularGraph.from_smiles(r_joined)
  p_mg = electron_graph.MolecularGraph.from_smiles(p_joined)

  r_ok, _, r_audit = (
      electron_graph.ReactionVerifierEngine.apply_and_verify_step(
          r_mg, [], step_name="reactants"
      )
  )
  p_ok, _, p_audit = (
      electron_graph.ReactionVerifierEngine.apply_and_verify_step(
          p_mg, [], step_name="products"
      )
  )
  r_errors = r_audit.get("errors", [])
  p_errors = p_audit.get("errors", [])

  r_counts = r_mg.nuclei_multiset()
  p_counts = p_mg.nuclei_multiset()

  # Determine atom differences (missing leaving group/condensate vs excess atoms)
  missing_in_products = {
      k: r_counts[k] - p_counts.get(k, 0)
      for k in r_counts
      if r_counts[k] > p_counts.get(k, 0)
  }
  excess_in_products = {
      k: p_counts[k] - r_counts.get(k, 0)
      for k in p_counts
      if p_counts[k] > r_counts.get(k, 0)
  }

  inferred_byproduct = None
  if missing_in_products and not excess_in_products:
    for name, smi, comp in COMMON_CONDENSATES:
      if comp == missing_in_products:
        inferred_byproduct = {"name": name, "smiles": smi, "atoms": comp}
        break

  mech_verified = None
  mech_errors: list[str] = []
  if mechanism_spec and mechanism_spec.get("arrows"):
    rxn_def = dict(mechanism_spec)
    rxn_def.setdefault("id", "eln_rxn")
    mg, arrows = render_mechanism.build_graph_and_arrows_from_spec(
        rxn_def, mechanism_spec
    )
    ok, _, audit = electron_graph.ReactionVerifierEngine.apply_and_verify_step(
        mg,
        arrows,
        step_name=rxn_def["id"],
        allow_multiproton=bool(mechanism_spec.get("allow_multiproton")),
    )
    mech_verified = ok
    mech_errors = audit.get("errors", [])

  is_valid = (
      r_ok and p_ok and not excess_in_products and (mech_verified is not False)
  )

  return {
      "electron_graph_valid": is_valid,
      "reactant_lewis_octet_valid": r_ok,
      "product_lewis_octet_valid": p_ok,
      "reactant_audit_errors": r_errors,
      "product_audit_errors": p_errors,
      "reactant_net_charge": r_mg.net_charge(),
      "product_net_charge": p_mg.net_charge(),
      "reactant_nuclei": dict(r_counts),
      "product_nuclei": dict(p_counts),
      "exact_atom_balance": r_counts == p_counts,
      "unlisted_byproduct_atoms": missing_in_products,
      "inferred_condensate_byproduct": inferred_byproduct,
      "unbalanced_excess_product_atoms": excess_in_products,
      "mechanism_arrow_verified": mech_verified,
      "mechanism_arrow_errors": mech_errors,
  }


def build_eln_table(
    rxn_smiles: str,
    limiting_mass_mg: float | None = None,
    limiting_mmol: float | None = None,
    reagents_spec: list[dict[str, Any]] | None = None,
    solvent: str = "THF",
    solvent_density_g_ml: float | None = None,
    concentration_molar: float = 0.20,
    actual_yield_mg: float | None = None,
    procedure_text: str | None = None,
    mechanism_spec: dict[str, Any] | None = None,
) -> dict[str, Any]:
  """Builds a complete ELN stoichiometry table and ElectronGraph sanity audit."""
  if ">>" not in rxn_smiles:
    raise ValueError(
        f"Expected reaction SMILES with '>>' separator, got: {rxn_smiles!r}"
    )
  lhs, _, rhs = rxn_smiles.split(">")
  reactant_smis = [s.strip() for s in lhs.split(".") if s.strip()]
  product_smis = [s.strip() for s in rhs.split(".") if s.strip()]
  if not reactant_smis or not product_smis:
    raise ValueError(
        "Reaction SMILES must contain at least one reactant and one product."
    )

  # 1. Run ElectronGraph formal sanity check
  sanity = sanity_check_reaction_with_electron_graph(
      reactant_smis, product_smis, mechanism_spec=mechanism_spec
  )

  # 2. Limiting reactant (first reactant by default)
  lim_mol = Chem.MolFromSmiles(reactant_smis[0])
  if lim_mol is None:
    raise ValueError(f"Invalid limiting reactant SMILES: {reactant_smis[0]!r}")
  lim_mw = float(Descriptors.MolWt(lim_mol))
  if limiting_mmol is not None:
    base_mmol = float(limiting_mmol)
    base_mass_mg = base_mmol * lim_mw
  elif limiting_mass_mg is not None:
    base_mass_mg = float(limiting_mass_mg)
    base_mmol = base_mass_mg / lim_mw
  else:
    base_mmol = 1.0
    base_mass_mg = base_mmol * lim_mw

  # 3. Assemble reagent rows
  rows: list[dict[str, Any]] = []
  if reagents_spec:
    raw_entries = reagents_spec
  else:
    raw_entries = []
    for idx, smi in enumerate(reactant_smis):
      raw_entries.append({
          "name": "Limiting Reactant" if idx == 0 else f"Reactant {idx + 1}",
          "role": "Limiting Reactant" if idx == 0 else "Reactant",
          "smiles": smi,
          "equiv": 1.0,
      })

  total_input_mass_mg = 0.0
  sum_stoich_reactant_mw = 0.0

  for idx, entry in enumerate(raw_entries):
    smi = entry.get("smiles", "")
    mol = Chem.MolFromSmiles(smi) if smi else None
    if mol is not None:
      mw = float(Descriptors.MolWt(mol))
      exact_mw = float(rdMolDescriptors.CalcExactMolWt(mol))
      formula = rdMolDescriptors.CalcMolFormula(mol)
    else:
      mw = float(entry["mw"])
      exact_mw = mw
      formula = entry.get("formula", "N/A")

    eq = float(entry.get("equiv", entry.get("eq", 1.0)))
    role = entry.get(
        "role", "Limiting Reactant" if idx == 0 else "Reactant/Reagent"
    )
    raw_wt = entry.get("wt_pct", entry.get("purity", 100.0))
    wt_pct = (
        float(raw_wt) * 100.0 if 0.0 < float(raw_wt) <= 1.0 else float(raw_wt)
    )
    density = entry.get("density_g_ml", entry.get("density"))
    molarity = entry.get("molarity_m", entry.get("molarity"))

    mmol = base_mmol * eq
    active_mass_mg = mmol * mw
    dispense_mass_mg = active_mass_mg / (wt_pct / 100.0)
    total_input_mass_mg += dispense_mass_mg

    if "Reactant" in role:
      sum_stoich_reactant_mw += mw

    vol_ul = None
    if molarity is not None and float(molarity) > 0:
      vol_ul = (mmol / float(molarity)) * 1000.0
    elif density is not None and float(density) > 0:
      vol_ul = dispense_mass_mg / float(density)

    rows.append({
        "compound": entry.get("name", f"Reagent {idx + 1}"),
        "role": role,
        "smiles": smi,
        "formula": formula,
        "mw_g_mol": round(mw, 2),
        "exact_mw": round(exact_mw, 4),
        "equiv": round(eq, 3),
        "mmol": round(mmol, 4),
        "active_mass_mg": round(active_mass_mg, 2),
        "dispense_mass_mg": round(dispense_mass_mg, 2),
        "wt_pct": wt_pct,
        "volume_ul": round(float(vol_ul), 1) if vol_ul is not None else None,
    })

  # 4. Desired product & theoretical yield
  prod_mol = Chem.MolFromSmiles(product_smis[0])
  if prod_mol is None:
    raise ValueError(f"Invalid target product SMILES: {product_smis[0]!r}")
  prod_mw = float(Descriptors.MolWt(prod_mol))
  prod_exact_mw = float(rdMolDescriptors.CalcExactMolWt(prod_mol))
  prod_formula = rdMolDescriptors.CalcMolFormula(prod_mol)
  theo_yield_mg = base_mmol * prod_mw

  solvent_vol_ml = (
      (base_mmol / concentration_molar) if concentration_molar > 0 else 0.0
  )
  resolved_solvent_density = (
      float(solvent_density_g_ml)
      if solvent_density_g_ml is not None
      else SOLVENT_DENSITIES_G_ML.get(solvent.strip().lower(), 0.889)
  )
  solvent_mass_mg = solvent_vol_ml * resolved_solvent_density * 1000.0
  total_process_mass_mg = total_input_mass_mg + solvent_mass_mg

  if sum_stoich_reactant_mw <= 0:
    sum_stoich_reactant_mw = sum(
        float(Descriptors.MolWt(Chem.MolFromSmiles(s))) for s in reactant_smis
    )
  atom_economy_pct = (prod_mw / sum_stoich_reactant_mw) * 100.0

  actual_yield_pct = None
  pmi = None
  if actual_yield_mg is not None and theo_yield_mg > 0:
    actual_yield_pct = (float(actual_yield_mg) / theo_yield_mg) * 100.0
    if float(actual_yield_mg) > 0:
      pmi = total_process_mass_mg / float(actual_yield_mg)
  elif theo_yield_mg > 0:
    pmi = total_process_mass_mg / theo_yield_mg

  # 5. Format clean ELN Markdown table
  md_lines = [
      f"## Electronic Lab Notebook (ELN) Entry — `{rxn_smiles}`",
      "",
      "### 1. ElectronGraph Reaction Sanity Check",
      (
          "- **ElectronGraph Status**:"
          f" `{'PASS' if sanity['electron_graph_valid'] else 'WARNING/FAIL'}`"
      ),
      (
          "- **Lewis Octet & Valence Valid**:"
          f" Reactants=`{sanity['reactant_lewis_octet_valid']}`,"
          f" Products=`{sanity['product_lewis_octet_valid']}`"
      ),
      (
          "- **Atom Balance**:"
          f" `{'Exact 100% match' if sanity['exact_atom_balance'] else 'Condensate/Leaving Group Omitted'}`"
      ),
  ]
  if sanity["inferred_condensate_byproduct"]:
    byp = sanity["inferred_condensate_byproduct"]
    md_lines.append(
        f"- **Inferred Condensate / Byproduct**: `{byp['name']}`"
        f" (`{byp['smiles']}`)"
    )
  elif sanity["unlisted_byproduct_atoms"]:
    md_lines.append(
        "- **Unlisted Leaving Group Atoms**:"
        f" `{sanity['unlisted_byproduct_atoms']}`"
    )
  if sanity["unbalanced_excess_product_atoms"]:
    md_lines.append(
        "- **ERROR — Excess Atoms in Product**:"
        f" `{sanity['unbalanced_excess_product_atoms']}`"
    )

  md_lines.extend([
      "",
      "### 2. Stoichiometry & Dispense Table",
      (
          "| Compound | Role | Formula | MW (g/mol) | Equiv | mmol | Mass (mg)"
          " | Volume (µL) |"
      ),
      "| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: |",
  ])
  for r in rows:
    vol_str = f"{r['volume_ul']:.1f}" if r["volume_ul"] is not None else "—"
    md_lines.append(
        f"| {r['compound']} | {r['role']} | `{r['formula']}` |"
        f" {r['mw_g_mol']:.2f} | {r['equiv']:.2f} | {r['mmol']:.4f} |"
        f" **{r['dispense_mass_mg']:.1f}** | {vol_str} |"
    )

  md_lines.extend([
      "",
      "### 3. Reaction Scale, Solvent & Yield Summary",
      (
          f"- **Solvent**: `{solvent}` (`{solvent_vol_ml:.2f} mL` for"
          f" `{concentration_molar:.2f} M` reaction concentration)"
      ),
      (
          f"- **Target Product**: `{prod_formula}` (`MW = {prod_mw:.2f} g/mol`,"
          f" `Exact Mass = {prod_exact_mw:.4f} Da`)"
      ),
      (
          f"- **Theoretical Yield**: **`{theo_yield_mg:.1f} mg`**"
          f" (`{base_mmol:.4f} mmol`)"
      ),
      f"- **Atom Economy**: `{atom_economy_pct:.1f}%`",
      f"- **Estimated Process Mass Intensity (PMI)**: `{pmi:.1f}`"
      if pmi
      else "",
  ])
  if actual_yield_pct is not None:
    md_lines.append(
        f"- **Isolated Yield**: **`{actual_yield_mg:.1f} mg`"
        f" (`{actual_yield_pct:.1f}%`)**"
    )
  if procedure_text:
    md_lines.extend(["", "### 4. Procedure Notes", procedure_text])

  return {
      "reaction_smiles": rxn_smiles,
      "electron_graph_sanity_check": sanity,
      "stoichiometry_table": rows,
      "reaction_summary": {
          "solvent": solvent,
          "solvent_density_g_ml": round(resolved_solvent_density, 4),
          "concentration_molar": concentration_molar,
          "solvent_volume_ml": round(solvent_vol_ml, 3),
          "product_formula": prod_formula,
          "product_mw_g_mol": round(prod_mw, 2),
          "product_exact_mass": round(prod_exact_mw, 4),
          "limiting_mmol": round(base_mmol, 4),
          "theoretical_yield_mg": round(theo_yield_mg, 2),
          "actual_yield_mg": actual_yield_mg,
          "actual_yield_pct": (
              round(actual_yield_pct, 2)
              if actual_yield_pct is not None
              else None
          ),
          "atom_economy_pct": round(atom_economy_pct, 2),
          "estimated_pmi": round(pmi, 2) if pmi is not None else None,
      },
      "eln_markdown": "\n".join(line for line in md_lines if line is not None),
  }


def main() -> None:
  parser = argparse.ArgumentParser(
      description=(
          "Generate an Electronic Lab Notebook (ELN) stoichiometry table and"
          " sanity-check reaction balance/octets via ElectronGraph."
      )
  )
  parser.add_argument(
      "--rxn",
      required=True,
      help="Reaction SMILES in 'R1.R2>>P1' format.",
  )
  parser.add_argument(
      "--limiting-mass-mg",
      type=float,
      default=None,
      help=(
          "Mass of limiting reactant in mg (default: 1.0 mmol scale if"
          " omitted)."
      ),
  )
  parser.add_argument(
      "--limiting-mmol",
      type=float,
      default=None,
      help="Molar scale of limiting reactant in mmol.",
  )
  parser.add_argument(
      "--reagents-json",
      default=None,
      help=(
          "Optional JSON string or file path listing reagent dicts with keys:"
          " name, role, smiles, equiv, density_g_ml, molarity_m, wt_pct."
      ),
  )
  parser.add_argument(
      "--solvent",
      default="THF",
      help="Reaction solvent name (default: THF).",
  )
  parser.add_argument(
      "--solvent-density",
      type=float,
      default=None,
      help=(
          "Optional solvent density in g/mL (defaults to physical lookup for"
          " --solvent)."
      ),
  )
  parser.add_argument(
      "--concentration-molar",
      type=float,
      default=0.20,
      help="Target molarity of limiting reactant (default: 0.20 M).",
  )
  parser.add_argument(
      "--actual-yield-mg",
      type=float,
      default=None,
      help="Optional isolated product mass in mg to calculate percent yield.",
  )
  parser.add_argument(
      "--procedure-text",
      default=None,
      help="Optional plain-text procedure notes to include in the ELN report.",
  )
  parser.add_argument(
      "--output",
      "-o",
      required=True,
      help="Output JSON file path.",
  )
  parser.add_argument(
      "--output-md",
      default=None,
      help=(
          "Optional output Markdown (.md) file path for the formatted ELN"
          " table."
      ),
  )
  args = parser.parse_args()

  reagents_spec = None
  if args.reagents_json:
    if os.path.exists(args.reagents_json):
      with open(args.reagents_json, "r", encoding="utf-8") as f:
        reagents_spec = json.load(f)
    else:
      reagents_spec = json.loads(args.reagents_json)

  try:
    result = build_eln_table(
        rxn_smiles=args.rxn,
        limiting_mass_mg=args.limiting_mass_mg,
        limiting_mmol=args.limiting_mmol,
        reagents_spec=reagents_spec,
        solvent=args.solvent,
        solvent_density_g_ml=args.solvent_density,
        concentration_molar=args.concentration_molar,
        actual_yield_mg=args.actual_yield_mg,
        procedure_text=args.procedure_text,
    )
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
      json.dump(result, f, indent=2)
    if args.output_md:
      os.makedirs(
          os.path.dirname(os.path.abspath(args.output_md)), exist_ok=True
      )
      with open(args.output_md, "w", encoding="utf-8") as f:
        f.write(result["eln_markdown"] + "\n")
    print(
        "[SUCCESS] ELN table & ElectronGraph audit complete"
        f" (valid={result['electron_graph_sanity_check']['electron_graph_valid']},"
        f" theoretical_yield={result['reaction_summary']['theoretical_yield_mg']:.1f}"
        f" mg) -> {args.output}"
    )
  except Exception as exc:
    print(f"[ERROR] {exc}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
  main()
