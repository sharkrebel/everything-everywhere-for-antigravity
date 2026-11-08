# Commercial Chemical Pricing & Screening Procurement Reference

This reference documents the multi-provider commercial pricing architecture, API
authorization mechanisms, package optimization algorithms, and standard
screening library catalog fallbacks.

## 1. Provider Ecosystem & Architecture

Chemical synthesis and compound procurement workflows require navigating
multiple distinct commercial databases:

-   **MolPort**
    -   **Library Coverage**: Global aggregated screening compounds (>7M
        stock, >20M MTO)
    -   **Primary Endpoints**: REST API `/v1/compounds/inchikey`
    -   **Auth Mechanism**: Bearer Token (`MOLPORT_API_KEY`)
-   **Chemspace**
    -   **Library Coverage**: Enamine REAL, ChemDiv, UORSY (>2B
        make-on-demand, >15M stock)
    -   **Primary Endpoints**: REST API `/v2/search/exact`
    -   **Auth Mechanism**: Bearer Token (`CHEMSPACE_API_KEY`)
-   **Enamine Store**
    -   **Library Coverage**: Building blocks & curated screening libraries
    -   **Primary Endpoints**: REST API
    -   **Auth Mechanism**: Customer API Key
-   **Mcule**
    -   **Library Coverage**: Curated stock & diversity libraries
    -   **Primary Endpoints**: REST API `/api/v1/`
    -   **Auth Mechanism**: API Token (`MCULE_API_KEY`)

## 2. Empirical Pricing Verification & Local Catalog Ingestion

Commercial chemical prices and stock availability change continuously and must
never be algorithmically estimated or simulated.

-   **Live Vendor APIs**: Configure `MOLPORT_API_KEY` or `CHEMSPACE_API_KEY`
    environment variables to query real-time commercial stock and pricing tiers.
-   **Offline / Local Vendor Quotes (`--catalog-csv`)**: When working offline or
    with pre-negotiated institutional vendor quotes, pass a CSV file via
    `--catalog-csv` containing columns `smiles` (or `inchikey`), `vendor`,
    `catalog_id`, `package_mg`, and `price_usd`.
-   **Missing API Keys / Unlisted Compounds**: If neither live API credentials
    nor a matching `--catalog-csv` entry are available,
    `screen_compound_pricer.py` explicitly reports `fulfilled: false`,
    `price_usd: null`, and `status: "unverified_missing_api_key_or_catalog"`.

## 3. Package Optimization Strategy

Given a target assay requirement (e.g., `--target-amount-mg 10.0`), the
optimizer selects the single package satisfying: $$\text{argmin}_{t \in
\text{tiers}} \{ t.\text{price\_usd} \mid t.\text{package\_mg} \ge
\text{target\_mg} \}$$

This guarantees researchers do not inadvertently purchase sub-optimal multi-vial
combinations when a single bulk tier offers lower net expenditure.

## 4. CLI Usage Examples

### Single Molecule Price Query

```bash
uv run scripts/screen_compound_pricer.py lookup-single \
  --smiles "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1" \
  --target-amount-mg 10.0
```

### Batch Virtual Screen / Med Chem Table Pricing

```bash
uv run scripts/screen_compound_pricer.py price-screen \
  --input screen_hits.csv \
  --output screen_hits_priced.csv \
  --target-amount-mg 10.0
```
