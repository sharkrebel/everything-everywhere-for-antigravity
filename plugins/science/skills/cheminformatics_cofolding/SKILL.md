---
name: cheminformatics-cofolding
description: >-
  Co-folds protein-ligand, protein-DNA/RNA, and covalent inhibitor complexes:
  prepares AlphaFold 3 JSON, Boltz-2 YAML, and Boltz API (boltz-2.1) inputs,
  extracts bound ligand poses from mmCIF/PDB, and computes 3D protein-ligand
  interaction fingerprints (PLIF) and pharmacophore overlaps.
metadata:
  publisher: google
---

# RDKit Biomolecular Co-Folding Hub (`cheminformatics-cofolding`)

## Prerequisites

1.  **`rdkit-general`**: Read the `rdkit-general` skill before running any
    workflow here.
2.  **`uv`**: Read the `uv` skill and follow its Setup instructions to ensure
    `uv` is installed and on PATH. Always execute skill scripts as `uv run
    /path/to/script.py ...` (or `python3 /path/to/script.py ...`); **never** run
    `uv run python3 /path/to/script.py`, as placing `python3` after `uv run`
    bypasses the script's PEP 723 `# /// script` dependency block.
3.  **Engine Readiness & License Prompt (`scripts/check_cofolding_setup.py`)**:
    Before running co-folding inference, run `uv run
    scripts/check_cofolding_setup.py` and prompt the user (unless the user
    already specified the engine and terms):
    -   **If AlphaFold 3 weights (`af3.bin.zst`) are detected** (via
        `--model-dir`, `AF3_MODEL_DIR`, or `~/.cache/alphafold3`): Ask the user
        to confirm running **AlphaFold 3** under the **Google DeepMind AlphaFold
        3 Model Parameters Terms of Use (Non-Commercial Use Only)**
        (`--accept-af3-terms`), running **Boltz-2 (`--engine boltz2`, MIT
        License, commercial-safe)**, or running the hosted **Boltz API
        (`--engine boltz-api`, `boltz-2.1`)**.
    -   **If AlphaFold 3 weights are not found**: Ask whether the user wants to
        (a) supply a local `--model-dir` (or set `AF3_MODEL_DIR`) containing
        `af3.bin.zst`, (b) run **Boltz-2 (`--engine boltz2`, MIT open-source
        weights)**, or (c) run the hosted **Boltz API (`--engine boltz-api`,
        `boltz-2.1`)**. Never silently substitute Boltz-2 when `--engine af3`
        was requested without asking the user.
    -   **Boltz API Credentials (`--engine boltz-api`, Optional)**: When using
        the hosted Boltz API (`https://api.boltz.bio`, model `boltz-2.1`), check
        `grep -sq "^BOLTZ_API_KEY=" ~/.env` via the `credentials` skill
        (optional `BOLTZ_BASE_URL` defaults to `https://api.boltz.bio`). Live
        keys (`sk_bc_*_live_*`) run GPU inference and return real structures
        plus `binding_confidence` and `optimization_score`, whereas test keys
        (`sk_bc_*_test_*`) return zero-cost synthetic placeholders for pipeline
        testing. Never print or inspect the contents of `~/.env`;
        `scripts/check_cofolding_setup.py` and
        `scripts/run_cofolding_inference.py` load `~/.env` automatically via
        `python-dotenv`.
    -   **OpenProtein Credentials (`--engine esmfold2`, Optional)**: When using
        the OpenProtein ESMFold2 engine, check `grep -sq
        "^OPENPROTEIN_USERNAME=" ~/.env && grep -sq "^OPENPROTEIN_API_KEY="
        ~/.env` via the `credentials` skill. Never print or inspect the contents
        of `~/.env`; `scripts/check_cofolding_setup.py` and
        `scripts/run_cofolding_inference.py` load `~/.env` automatically via
        `python-dotenv`.
    -   Record the accepted license notice and timestamp in
        `.licenses/cheminformatics_cofolding_LICENSE.txt`.

## Core Rules

-   **Shared foundation rules**: Follow all core rules in
    `../rdkit_general/SKILL.md`.
-   **CPU / Apple Silicon Execution (`--skip-data-pipeline`)**: Do not skip
    co-folding on CPU-only or Apple Silicon hosts; pass `--skip-data-pipeline`
    (`--msa-free` / `msa: empty`) when running without an NVIDIA GPU so
    `run_cofolding_inference.py` applies CPU/MPS flags, or use `--engine
    boltz-api` to offload structure and binding inference to the Boltz cloud
    API.
-   **100% Dataset Processing**: Process all quantitatively labeled compounds
    and score all prospective patent structures in a series without arbitrary
    subsampling.
-   **Strict Engine & License Disclosure**: Inspect `inference_metadata.json` /
    `confidence_summary.json` and explicitly report which co-folding engine
    (**AlphaFold 3**, **Boltz-2**, **Boltz API (`boltz-2.1`)**, or **ESMFold2**)
    and license generated the 3D complex.
-   **Explicit Ligand Chain Selection (`--ligand-chain`)**: When extracting
    poses or contacts with `scripts/extract_cofolded_complex.py`, always pass
    `--ligand-chain` and `--ligand-smiles` explicitly whenever known so a
    protein chain is never misidentified as a small-molecule ligand.

--------------------------------------------------------------------------------

## Task-to-Domain Routing Directory

Match the user's objective to the domain below, **read the corresponding
reference guide**, and execute from `scripts/`:

1.  **Environment Check, Input Preparation & Multi-Engine Inference**

    -   **Reference Guide**: `references/af3_input_schema_and_ligands.md` — read
        when preparing AF3 JSON, Boltz-2 YAML, or Boltz API (`boltz-2.1`) JSON
        inputs, covalent `bondedAtomPairs`, pocket constraints
        (`--pocket-residues`), DNA/RNA duplex chains, RDKit 3D pre-flight
        checks, or resolving multi-engine weights/caches for `af3`, `boltz2`,
        `boltz-api`, and `esmfold2`.
    -   **Standard Executables**:
        -   `uv run scripts/check_cofolding_setup.py` — diagnoses container
            runtime, GPU/MPS availability, and installed engines/weights/API
            credentials.
        -   `uv run scripts/prepare_cofolding_input.py` — writes AF3 JSON
            (`--output-json`), Boltz-2 YAML (`--output-boltz-yaml`), and Boltz
            API JSON (`--output-boltz-api-json`) for protein, SMILES/CCD ligand,
            DNA/RNA, and covalent complexes, with an RDKit 3D conformer
            pre-flight check.
        -   `uv run scripts/run_cofolding_inference.py` — executes co-folding
            inference (`--engine {af3,boltz2,boltz-api,esmfold2,auto}`,
            `--json-path`, `--output-dir`, `--boltz-model`, `--num-samples`).

2.  **Pose Extraction, PLIF Contacts & Pharmacophore Overlap**

    -   **Reference Guide**: `references/af3_cadd_and_pose_analysis.md` — read
        when extracting predicted ligand poses from mmCIF/PDB, computing 3D
        protein-ligand interaction fingerprints (PLIF), or measuring 3D
        pharmacophore feature map overlap.
    -   **Standard Executables**:
        -   `uv run scripts/extract_cofolded_complex.py` — extracts the
            predicted ligand from an output mmCIF/PDB with exact bond orders and
            confidence/binding metrics (`--ligand-chain`, `--ligand-smiles`,
            `--output-sdf`, `--output-csv`).
        -   `uv run scripts/compute_plif.py` — computes protein-ligand
            interaction fingerprints (H-bonds, $\pi$-stacking, salt bridges,
            hydrophobic contacts) from an mmCIF/PDB complex.
        -   `uv run scripts/feature_map_3d_scorer.py` — computes 3D RDKit
            pharmacophore feature map overlap between aligned poses.

--------------------------------------------------------------------------------

## Progressive Execution Protocol

1.  **Check Setup & Confirm Engine**: Run `uv run
    scripts/check_cofolding_setup.py` and confirm engine selection (`af3` vs.
    `boltz2` vs. `boltz-api` vs. `esmfold2`; see
    `references/af3_input_schema_and_ligands.md`).
2.  **Prepare Inputs & Run Inference**: Run `prepare_cofolding_input.py`
    (`--output-json`, `--output-boltz-yaml`, `--output-boltz-api-json`,
    `--msa-free`) and execute `run_cofolding_inference.py`.
3.  **Extract Pose & Profile Interactions**: Run `extract_cofolded_complex.py`
    (with explicit `--ligand-chain` and `--ligand-smiles`) and
    `compute_plif.py`.
4.  **Summarize Findings in Markdown**: Report the active co-folding engine,
    `iptm`, `plddt`, Boltz `binding_confidence` / `optimization_score` (when run
    via `boltz-api`), binding contacts, and GFM tables in a markdown summary in
    the output directory.
