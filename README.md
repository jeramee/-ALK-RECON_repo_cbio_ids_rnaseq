# ALK-RECON

**ALK-RECON** is a research / education pipeline that converts ALK-related
ctDNA / NGS findings into an **auditable resistance-mechanism dossier**.

The core idea is simple but strict:

> Separate **what was observed**, **what was inferred**, and **what is speculative** —
> and preserve the evidence trail for each.

⚠️ Research and education only — **not medical advice**.

---

## What this repository does

ALK-RECON:
- Ingests mutation tables (synthetic or real-world)
- Optionally ingests RNA-seq expression data
- Groups findings into deterministic **case snapshots**
- Routes each case through explicit resistance-mechanism logic
- Produces human-readable + machine-readable dossiers

---

## How to read this repository

| If you want to… | Start here |
|---|---|
| Run the pipeline quickly | [Quick start](docs/README_quickstart.md) |
| Understand the test data | [Synthetic test data](docs/README_test_data.md) |
| Understand ingest internals | [Ingest internals](docs/README_ingest_internals.md) |
| Trust behavior via tests | [Pytest regression](docs/README_pytest_regression.md) |
| Swap in real data | [cBioPortal swap](docs/README_cbioportal_swap.md) |

---

## License

MIT
