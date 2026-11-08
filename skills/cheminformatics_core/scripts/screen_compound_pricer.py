# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "numpy",
#     "pandas",
#     "polite-http",
#     "python-dotenv",
# ]
# ///

"""Unified Commercial Chemical Pricing and Screening Re-Supply Optimizer.

Queries commercial supplier pricing and available package tiers (MolPort,
Chemspace, Mcule, Enamine, Sigma-Aldrich), normalizes price-per-milligram,
optimizes purchase packaging for target assay amounts, and computes campaign
procurement budgets for virtual screening and med chem hit lists.
"""

import abc
import argparse
import json
import os
from typing import Any, Dict, List, Optional
import urllib.parse

import dotenv
import pandas as pd
from polite_http import http_client
from rdkit import Chem

dotenv.load_dotenv(os.path.expanduser("~/.env"))


class PriceTier:
  """Represents a discrete commercial packaging tier."""

  def __init__(
      self,
      vendor: str,
      catalog_id: str,
      package_mg: float,
      price_usd: float,
      in_stock: bool = True,
      lead_time_days: Optional[int] = None,
      purity_percent: Optional[float] = None,
  ):
    self.vendor = vendor
    self.catalog_id = catalog_id
    self.package_mg = float(package_mg)
    self.price_usd = float(price_usd)
    self.price_per_mg = (
        self.price_usd / self.package_mg if self.package_mg > 0 else 0.0
    )
    self.in_stock = bool(in_stock)
    self.lead_time_days = (
        int(lead_time_days) if lead_time_days is not None else None
    )
    self.purity_percent = (
        float(purity_percent) if purity_percent is not None else None
    )

  def to_dict(self) -> Dict[str, Any]:
    return {
        "vendor": self.vendor,
        "catalog_id": self.catalog_id,
        "package_mg": round(self.package_mg, 2),
        "price_usd": round(self.price_usd, 2),
        "price_per_mg": round(self.price_per_mg, 3),
        "in_stock": self.in_stock,
        "lead_time_days": self.lead_time_days,
        "purity_percent": self.purity_percent,
    }


class BasePricingProvider(abc.ABC):
  """Abstract base class for chemical pricing providers."""

  @abc.abstractmethod
  def get_provider_name(self) -> str:
    pass

  @abc.abstractmethod
  def query_compound_prices(
      self, mol: Chem.Mol, smiles: str, inchikey: str
  ) -> List[PriceTier]:
    pass


class LocalCatalogCsvProvider(BasePricingProvider):
  """Ingests verified vendor quotes or institutional catalog files from a local CSV."""

  def __init__(self, catalog_csv: Optional[str] = None):
    self.catalog_csv = catalog_csv
    self.df: Optional[pd.DataFrame] = None
    if catalog_csv and os.path.exists(catalog_csv):
      self.df = pd.read_csv(catalog_csv)

  def get_provider_name(self) -> str:
    return "LocalCatalogCsv"

  def query_compound_prices(
      self, mol: Chem.Mol, smiles: str, inchikey: str
  ) -> List[PriceTier]:
    if self.df is None or self.df.empty:
      return []

    tiers: List[PriceTier] = []
    # Match by inchikey or smiles
    matches = pd.DataFrame()
    if "inchikey" in self.df.columns and inchikey:
      matches = self.df[self.df["inchikey"] == inchikey]
    elif "smiles" in self.df.columns and smiles:
      matches = self.df[self.df["smiles"] == smiles]

    for _, row in matches.iterrows():
      if pd.isna(row.get("package_mg")) or pd.isna(row.get("price_usd")):
        continue
      lead_val = row.get("lead_time_days")
      purity_val = row.get("purity_percent")
      tiers.append(
          PriceTier(
              vendor=str(row.get("vendor", "CatalogCSV")),
              catalog_id=str(row.get("catalog_id", "UNKNOWN")),
              package_mg=float(row["package_mg"]),
              price_usd=float(row["price_usd"]),
              in_stock=bool(row.get("in_stock", True)),
              lead_time_days=int(lead_val) if pd.notna(lead_val) else None,
              purity_percent=float(purity_val)
              if pd.notna(purity_val)
              else None,
          )
      )
    return tiers


class LiveApiPricingProvider(BasePricingProvider):
  """Integrates live vendor REST endpoints when API keys are configured."""

  def __init__(self):
    self.molport_api_key = os.environ.get("MOLPORT_API_KEY")
    self.chemspace_api_key = os.environ.get("CHEMSPACE_API_KEY")
    self._molport_client = http_client.HttpClient(
        "https://api.molport.com/", qps=5.0
    )
    self._chemspace_client = http_client.HttpClient(
        "https://api.chem-space.com/", qps=5.0
    )

  def get_provider_name(self) -> str:
    return "LiveVendorApi"

  def query_compound_prices(
      self, mol: Chem.Mol, smiles: str, inchikey: str
  ) -> List[PriceTier]:
    tiers = []
    # MolPort REST query if key present
    if self.molport_api_key and inchikey:
      try:
        data = self._molport_client.fetch_json(
            f"v1/compounds/inchikey/{inchikey}",
            headers={"Authorization": f"Bearer {self.molport_api_key}"},
            timeout=5.0,
        )
        for offer in data.get("offers", []):
          if "amount_mg" not in offer or "price_usd" not in offer:
            continue
          lead_val = offer.get("lead_time_days")
          purity_val = offer.get("purity")
          tiers.append(
              PriceTier(
                  vendor="MolPort",
                  catalog_id=offer.get("catalog_number", "MP-UNKNOWN"),
                  package_mg=float(offer["amount_mg"]),
                  price_usd=float(offer["price_usd"]),
                  in_stock=offer.get("stock", True),
                  lead_time_days=int(lead_val)
                  if lead_val is not None
                  else None,
                  purity_percent=float(purity_val)
                  if purity_val is not None
                  else None,
              )
          )
      except Exception:
        pass

    # Chemspace REST query if key present
    if self.chemspace_api_key and smiles:
      try:
        enc_smi = urllib.parse.quote(smiles)
        data = self._chemspace_client.fetch_json(
            f"v2/search/exact?smiles={enc_smi}",
            headers={"Authorization": f"Bearer {self.chemspace_api_key}"},
            timeout=5.0,
        )
        for item in data.get("items", []):
          for pack in item.get("pack_sizes", []):
            if "mg" not in pack or "price_usd" not in pack:
              continue
            lead_val = pack.get("lead_time_days")
            purity_val = item.get("purity")
            tiers.append(
                PriceTier(
                    vendor="Chemspace",
                    catalog_id=item.get("catalog_number", "CS-UNKNOWN"),
                    package_mg=float(pack["mg"]),
                    price_usd=float(pack["price_usd"]),
                    in_stock=item.get("in_stock", True),
                    lead_time_days=int(lead_val)
                    if lead_val is not None
                    else None,
                    purity_percent=float(purity_val)
                    if purity_val is not None
                    else None,
                )
            )
      except Exception:
        pass

    return tiers


class UnifiedScreeningPricer:
  """Coordinates multi-provider commercial queries, package optimization, and screening budgets."""

  def __init__(
      self, use_live_api: bool = True, catalog_csv: Optional[str] = None
  ):
    self.providers: List[BasePricingProvider] = []
    if catalog_csv:
      self.providers.append(LocalCatalogCsvProvider(catalog_csv))
    if use_live_api:
      self.providers.append(LiveApiPricingProvider())

  def get_compound_pricing(self, smiles: str) -> Dict[str, Any]:
    """Resolves all available pricing tiers for a single molecule."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
      return {"smiles": smiles, "error": "INVALID_SMILES", "tiers": []}

    canonical_smiles = Chem.MolToSmiles(mol)
    inchikey = Chem.MolToInchiKey(mol)

    all_tiers: List[PriceTier] = []
    for provider in self.providers:
      res = provider.query_compound_prices(mol, canonical_smiles, inchikey)
      if res:
        all_tiers.extend(res)

    all_tiers.sort(key=lambda t: (t.price_per_mg, t.price_usd))

    return {
        "canonical_smiles": canonical_smiles,
        "inchikey": inchikey,
        "total_options_found": len(all_tiers),
        "tiers": [t.to_dict() for t in all_tiers],
    }

  def optimize_procurement_for_amount(
      self, tiers: List[Dict[str, Any]], target_amount_mg: float
  ) -> Dict[str, Any]:
    """Selects the single lowest-cost commercial package that fulfills target_amount_mg."""
    if not tiers:
      return {
          "fulfilled": False,
          "status": "unverified_missing_api_key_or_catalog",
          "vendor": None,
          "catalog_id": None,
          "package_mg": None,
          "price_usd": None,
          "price_per_mg": None,
          "in_stock": False,
          "lead_time_days": None,
          "note": (
              "No verified commercial pricing found. Configure MOLPORT_API_KEY"
              " or CHEMSPACE_API_KEY, or pass --catalog-csv with vendor quotes."
          ),
      }

    # Filter packs that individually cover the target amount
    qualifying = [t for t in tiers if t["package_mg"] >= target_amount_mg]
    if qualifying:
      best = min(qualifying, key=lambda t: t["price_usd"])
      return {
          "fulfilled": True,
          "status": "verified",
          "vendor": best["vendor"],
          "catalog_id": best["catalog_id"],
          "package_mg": best["package_mg"],
          "price_usd": best["price_usd"],
          "price_per_mg": best["price_per_mg"],
          "in_stock": best["in_stock"],
          "lead_time_days": best["lead_time_days"],
      }

    # If no single package covers target amount, take largest available
    largest = max(tiers, key=lambda t: t["package_mg"])
    return {
        "fulfilled": False,
        "status": "insufficient_package_size",
        "note": (
            f"Largest pack ({largest['package_mg']} mg) is smaller than target"
            f" ({target_amount_mg} mg)"
        ),
        "vendor": largest["vendor"],
        "catalog_id": largest["catalog_id"],
        "package_mg": largest["package_mg"],
        "price_usd": largest["price_usd"],
        "price_per_mg": largest["price_per_mg"],
        "in_stock": largest["in_stock"],
        "lead_time_days": largest["lead_time_days"],
    }

  def price_screening_table(
      self,
      input_csv: str,
      output_csv: Optional[str] = None,
      target_amount_mg: float = 10.0,
      smiles_col: Optional[str] = None,
      id_col: Optional[str] = None,
  ) -> Dict[str, Any]:
    """Processes an entire screening table, optimizing pricing and aggregating budget."""
    df = pd.read_csv(input_csv)

    # Detect SMILES column
    if smiles_col and smiles_col in df.columns:
      s_col = smiles_col
    else:
      candidates = [
          c
          for c in df.columns
          if c.lower() in ("smiles", "structure", "canonical_smiles", "smi")
      ]
      s_col = candidates[0] if candidates else df.columns[1]

    # Detect ID column
    if id_col and id_col in df.columns:
      _ = id_col
    else:
      candidates = [
          c
          for c in df.columns
          if c.lower()
          in ("compound_id", "id", "name", "molecule_id", "chembl_id")
      ]
      _ = candidates[0] if candidates else df.columns[0]

    results = []
    for _, row in df.iterrows():
      smi = str(row[s_col])
      price_info = self.get_compound_pricing(smi)
      best_opt = self.optimize_procurement_for_amount(
          price_info.get("tiers", []), target_amount_mg
      )

      rec = dict(row)
      rec["canonical_smiles"] = price_info.get("canonical_smiles", smi)
      rec["inchikey"] = price_info.get("inchikey", "")
      rec["target_procure_mg"] = target_amount_mg
      rec["pricing_status"] = best_opt.get("status", "unverified")
      rec["optimal_vendor"] = best_opt["vendor"]
      rec["optimal_catalog_id"] = best_opt["catalog_id"]
      rec["optimal_pack_mg"] = best_opt["package_mg"]
      rec["optimal_price_usd"] = best_opt["price_usd"]
      rec["optimal_price_per_mg"] = best_opt["price_per_mg"]
      rec["optimal_in_stock"] = best_opt["in_stock"]
      rec["optimal_lead_time_days"] = best_opt["lead_time_days"]
      results.append(rec)

    res_df = pd.DataFrame(results)

    if output_csv:
      os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
      res_df.to_csv(output_csv, index=False)

    total_compounds = len(res_df)
    verified_mask = res_df["optimal_price_usd"].notna()
    verified_count = int(verified_mask.sum())
    total_budget_usd = (
        float(res_df.loc[verified_mask, "optimal_price_usd"].sum())
        if verified_count > 0
        else 0.0
    )
    avg_price_per_mg = (
        float(res_df.loc[verified_mask, "optimal_price_per_mg"].mean())
        if verified_count > 0
        else 0.0
    )
    in_stock_pct = (
        (res_df["optimal_in_stock"].sum() / total_compounds * 100.0)
        if total_compounds > 0
        else 0.0
    )

    summary = {
        "total_screen_compounds": total_compounds,
        "verified_priced_compounds": verified_count,
        "unverified_compounds": total_compounds - verified_count,
        "target_amount_per_compound_mg": target_amount_mg,
        "total_procurement_budget_usd": round(total_budget_usd, 2),
        "average_price_per_mg_usd": round(avg_price_per_mg, 3),
        "in_stock_availability_percent": round(float(in_stock_pct), 1),
        "output_table": output_csv if output_csv else "NOT_SAVED",
        "vendor_breakdown": (
            res_df["optimal_vendor"].dropna().value_counts().to_dict()
        ),
    }
    return summary


def main():
  parser = argparse.ArgumentParser(
      description=(
          "Unified Commercial Screening Compound Pricer & Budget Optimizer."
      )
  )
  parser.add_argument(
      "--catalog-csv",
      default=None,
      help="Optional path to a local vendor quote/catalog CSV file.",
  )
  parser.add_argument(
      "--smiles", default=None, help="Compound SMILES string (single lookup)."
  )
  parser.add_argument(
      "--input",
      "--input-csv",
      "--input_csv",
      dest="input",
      default=None,
      help="Input CSV/TSV table containing screening hits.",
  )
  parser.add_argument(
      "--output",
      "--output-csv",
      "--output_csv",
      dest="output",
      default=None,
      help="Output enriched CSV with pricing and vendor recommendations.",
  )
  parser.add_argument(
      "--target-amount-mg",
      "--target_amount_mg",
      dest="target_amount_mg",
      type=float,
      default=10.0,
      help="Target assay quantity in mg.",
  )
  parser.add_argument(
      "--smiles-col",
      "--smiles_col",
      dest="smiles_col",
      default=None,
      help="Name of SMILES column (auto-detected if omitted).",
  )
  parser.add_argument(
      "--id-col",
      "--id_col",
      dest="id_col",
      default=None,
      help="Name of Compound ID column (auto-detected if omitted).",
  )
  subparsers = parser.add_subparsers(dest="subcommand", required=False)

  # Subcommand: lookup-single
  single_p = subparsers.add_parser(
      "lookup-single",
      help="Query commercial pricing tiers for a single compound.",
  )
  single_p.add_argument(
      "--smiles", required=True, help="Compound SMILES string."
  )
  single_p.add_argument(
      "--target-amount-mg",
      "--target_amount_mg",
      dest="target_amount_mg",
      type=float,
      default=10.0,
      help="Target assay quantity in mg.",
  )

  # Subcommand: price-screen
  screen_p = subparsers.add_parser(
      "price-screen",
      help="Price and optimize procurement for a screening table.",
  )
  screen_p.add_argument(
      "--input",
      "--input-csv",
      "--input_csv",
      dest="input",
      required=True,
      help="Input CSV/TSV table containing screening hits.",
  )
  screen_p.add_argument(
      "--output",
      "--output-csv",
      "--output_csv",
      dest="output",
      required=True,
      help="Output enriched CSV with pricing and vendor recommendations.",
  )
  screen_p.add_argument(
      "--target-amount-mg",
      "--target_amount_mg",
      dest="target_amount_mg",
      type=float,
      default=10.0,
      help="Target assay quantity in mg per compound.",
  )
  screen_p.add_argument(
      "--smiles-col",
      "--smiles_col",
      dest="smiles_col",
      default=None,
      help="Name of SMILES column (auto-detected if omitted).",
  )
  screen_p.add_argument(
      "--id-col",
      "--id_col",
      dest="id_col",
      default=None,
      help="Name of Compound ID column (auto-detected if omitted).",
  )

  args = parser.parse_args()
  pricer = UnifiedScreeningPricer(catalog_csv=args.catalog_csv)

  mode = args.subcommand
  if not mode:
    if args.input:
      mode = "price-screen"
    elif args.smiles:
      mode = "lookup-single"
    else:
      parser.error(
          "Specify a subcommand (lookup-single / price-screen) or pass --smiles"
          " or --input."
      )

  if mode == "lookup-single":
    pricing = pricer.get_compound_pricing(args.smiles)
    best = pricer.optimize_procurement_for_amount(
        pricing.get("tiers", []), args.target_amount_mg
    )
    output = {
        "query_smiles": args.smiles,
        "inchikey": pricing.get("inchikey"),
        "target_amount_mg": args.target_amount_mg,
        "optimal_procurement": best,
        "all_available_tiers": pricing.get("tiers", []),
    }
    if args.output:
      with open(args.output, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)
    print(json.dumps(output, indent=2))

  elif mode == "price-screen":
    out_csv = args.output or "priced_screening_hits.csv"
    summary = pricer.price_screening_table(
        input_csv=args.input,
        output_csv=out_csv,
        target_amount_mg=args.target_amount_mg,
        smiles_col=args.smiles_col,
        id_col=args.id_col,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
  main()
