# Quick Start — ALK-RECON

## Minimal run
```bash
python -m alk_recon.cli run   --variants data/test/variants_synthetic_alk.tsv   --out out
```

## With RNA-seq
```bash
python -m alk_recon.cli run   --variants data/test/variants_synthetic_alk.tsv   --rnaseq-counts data/test/rnaseq_expression.tsv   --rnaseq-meta data/test/sample_map.tsv   --out out
```
