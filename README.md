ALK-RECON (cbio_ids + rnaseq)
ALK-RECON is a research/education pipeline starter kit that turns ALK-related ctDNA/NGS findings into an auditable resistance-mechanism report.
This repo flavor supports cBioPortal study + sample IDs and optional RNA-seq counts + metadata (“cbio_ids rnaseq”).
Not medical advice.

Project title
ALK Lockpick — Stop resistance from changing the locks.

What this project is
ALK-RECON organizes one core idea: don’t treat ALK resistance as “one mutation = one answer.”
Classify why escape happened (mutation vs bypass vs persistence) and generate a transparent evidence-led dossier + research routing buckets.

Use cases (non-clinical):

Scientist: mechanism framing + validation steps + failure modes
Clinician: plain-English map + data-driven logic (not a treatment recommendation)
Investor: why the problem is hard + why the approach is coherent + what “success” means
Scope boundary: hypotheses + audit trails for research/education. No diagnosis. No patient-specific recommendations.

Why ALK resistance is hard
ALK-targeted therapy can work extremely well—until it doesn’t. Failures often mix three modes:

On-target ALK changes: binding/conformation shifts (e.g., solvent-front G1202R, gatekeeper L1196M, compound ALK mutations)
Off-target bypass: alternate drivers restore survival signaling (common examples: MET, EGFR, MAPK)
Tolerance / persistence: drug-tolerant survivors without stable resistance mutations (often described with HDAC/DNMT/BET-type state logic)
ALK-RECON exists to separate these modes explicitly and make reasoning reusable.

KNOWN / INFERRED / SPECULATIVE rule
Every dossier separates claims into:

KNOWN: directly observed in inputs (variants, CNAs, RNA counts, metadata)
INFERRED: derived features (e.g., “solvent-front category” from G1202R)
SPECULATIVE: clearly labeled hypotheses (e.g., “epigenetic adjunct might suppress persisters”)
What the pipeline does (MVP)
Inputs (“front doors”)
Variant table (CSV/TSV) — recommended MVP input
Optional cBioPortal columns: study_id, sample_id
Optional RNA-seq: count matrix (genes×samples) + sample metadata
Outputs (per case / sample / timepoint)
Mechanism calls (ranked): on_target_alk, bypass, persistence
Evidence ledger: which inputs/features triggered each call + provenance hooks
Strategy routing buckets (research framing; not advice):
A mutation-aware ATP-site inhibitor (next-gen)
B allosteric / conformation-lock thesis (αC-helix / A-loop concepts)
C targeted degradation (PROTAC / molecular glue / degrader)
D persistence suppression adjuncts (HDAC/DNMT/BET-type logic)
E sequencing / monitoring logic (portfolio rules)
Exports: Markdown dossier + JSON dossier
What “auditable” means here
If a claim can’t point to an evidence ID, it shouldn’t be in the report.
Every dossier aims to answer: What did we see (KNOWN)? What did we infer? What’s speculative? What would falsify it?

Repository anatomy
schema/ — CaseSnapshot object + JSON schema export (case_snapshot.py)
ingest/ — variant table import; optional RNA-seq ingest
features/ — engineered flags (ALK categories, bypass indicators, persistence scores)
mechanism_engine/ — rule-based scoring + explanation assembly (auditable)
reports/ — dossier generator (Markdown + JSON)
llm_layer/ — “safe narrator” (structured-in → structured-out; no invention)
tests/ — golden cases (G1202R, L1196M, compound, bypass, persistence)
Optional docs you can add under docs/: 01_GLOSSARY.md, 02_MECHANISM_MAP.md, 03_STRATEGY_OPTIONS.md, 04_DECISION_TREE.md, 05_VALIDATION_PLAN.md, 07_EXEC_SUMMARY_1PAGE.md

Quickstart
Assumes a CLI entry at alk_recon.cli. If yours differs, adjust.

Install
git clone https://github.com/jeramee/ALK-RECON.git
cd ALK-RECON
python -m venv .venv
# Windows: .venv\Scripts\activate
# Mac/Linux: source .venv/bin/activate
pip install -e .
Run (variants only)
python -m alk_recon.cli --input examples/variants_sample.tsv --out out
Run (variants + RNA-seq)
python -m alk_recon.cli --input examples/variants_sample.tsv --rnaseq-counts examples/rnaseq_counts_example.tsv --rnaseq-meta examples/rnaseq_meta_example.tsv --out out
Outputs
out/dossiers/<case_id>.md, out/dossiers/<case_id>.json, out/index.json

Input formats (minimal)
Variant table (CSV/TSV)
Required: case_id, gene
Recommended: study_id, sample_id, timepoint_id, protein_change, variant_type, vaf, copy_number, notes

Example (TSV):

study_id	case_id	sample_id	timepoint_id	gene	protein_change	variant_type	vaf
luad_tcga_pan_can_atlas_2018	CASE_0001	TCGA-XX-YYYY-01	baseline	ALK	G1202R	SNV	0.22
luad_tcga_pan_can_atlas_2018	CASE_0001	TCGA-XX-YYYY-01	baseline	MET		AMP	
RNA-seq counts (optional)
genes × samples; first column gene id; remaining columns sample IDs matching metadata sample_id.

RNA-seq metadata (optional)
Required: sample_id
Recommended: case_id, study_id, timepoint_id, condition, batch

Output dossier layout (what you should expect)
A typical out/dossiers/<case_id>.md is structured like:

Header: case/sample/timepoint identifiers
One-paragraph summary (plain English)
Mechanism ranking table (on_target_alk / bypass / persistence)
“Top evidence” bullets (3–10 items)
Evidence ledger table (feature → source datum → evidence_id)
Routing buckets (A–E) with short “why” lines
Uncertainty tags (KNOWN/INFERRED/SPECULATIVE) applied inline
“What to test next” (research validation suggestions)
“How this fails” (what would falsify the leading explanation)
The JSON output mirrors the same structure, but keeps everything machine-readable.

Evidence ledger (minimal schema idea)
The evidence ledger is the heart of “auditable.” Each entry is meant to be:

small
specific
attributable to an input row / gene / sample
Typical fields:

evidence_id (stable within a dossier)
case_id, sample_id, timepoint_id
layer (KNOWN / INFERRED / SPECULATIVE)
feature_name (e.g., has_G1202R)
feature_value (boolean / numeric / string)
source_type (variants / cnv / rnaseq / metadata / note)
source_ref (file + row index OR equivalent pointer)
gene (if relevant)
note (short, optional)
This is intentionally simple so you can export, audit, and explain it.

CaseSnapshot (unit of analysis)
If you can’t attach a datum to a CaseSnapshot, it’s not pipeline-ready.
A snapshot bundles identity (case/sample/timepoint), alterations, optional expression summary, an evidence ledger, mechanism calls, and routing.

Mechanism engine (rule-first)
Rule-first is intentional (labels are scarce/noisy; interpretability matters). Typical feature families:

On-target ALK: G1202R, L1196M, compound patterns, mutation counts, fusion presence (if provided)
Bypass: MET/EGFR/MAPK node events, “other driver emerged”
Persistence (optional): signature score(s), marker hits, longitudinal patterns (if present)
Contradictions are recorded; routing is non-prescriptive.

Configuration (optional, recommended)
If you want this to be reusable across cohorts, keep config in a YAML file (e.g., config/alk_recon.yaml) and keep the CLI clean.

Example YAML sketch:

genes:
  core: [ALK, MET, EGFR, KRAS, BRAF, MAP2K1]
thresholds:
  vaf_min: 0.01
  cn_amp_min: 6
persistence:
  enabled: true
  signatures:
    persister_v1:
      genes_up: [AXL, VIM]
      genes_dn: [EPCAM]
routing:
  prefer_allosteric_if: [has_G1202R, has_compound_alk_mutations]
LLM layer (narrator only)
LLMs can improve readability, but must cite evidence IDs, never invent facts, and label uncertainty.
The structured pipeline output remains the source of truth.

Golden cases (what tests should cover)
At minimum, tests should include:

ALK solvent-front dominant (G1202R)
ALK gatekeeper dominant (L1196M)
ALK compound mutation (two or more ALK resistance mutations)
bypass-dominant (strong MET/EGFR/MAPK evidence without strong ALK mutation)
mixed mode (ALK mutation + bypass evidence together)
persistence-leaning (weak on-target/bypass, RNA signature elevated)
The mystery “M-” drug (uncertainty handling)
If a remembered “M-” adjunct exists but is unknown, treat it as a hypothesis tag until identified.
Fast recall candidates: mocetinostat (HDAC), mirdametinib (MEK), momelotinib (JAK).
Disambiguation cues: suffix -stat / -tinib / -mab, and whether the discussion was chromatin/HDAC/BET vs MAPK rebound vs inflammation/JAK.

Validation (program reality check)
This repo is “AI makes hypotheses auditable,” not “AI cures cancer.”
Meaningful evidence usually requires: engagement → signaling shutdown → phenotype, mutant coverage (single+compound), and durability assays (days–weeks), not just short-term viability.

Non-goals (deliberate)
Not a clinical decision support tool
Not a treatment recommender
Not a black-box ML model that hides its reasoning
Not a replacement for wet-lab validation
Troubleshooting
VAF units: don’t mix 0–1 and 0–100 without normalization.
Protein change formatting: keep consistent (e.g., G1202R vs p.G1202R).
RNA-seq: counts columns must match metadata sample_id.
Timepoints: keep timepoint_id consistent (baseline/progression/etc.).
Safety / disclaimer
No medical advice. Research/education only. If you use real patient data, you are responsible for compliance and de-identification.

Roadmap (practical)
Near-term: tighter validators, more golden cases, simple report viewer.
Medium: optional cBioPortal fetch mode, configurable signatures (YAML), cohort summaries.
Later: optional structure sidecar for pocket adjacency / allosteric hypothesis mapping (not required for MVP).

Contributing (lightweight)
Keep new outputs auditable (evidence IDs + provenance).
Don’t add clinical language that reads like “recommendation.”
Prefer small, testable rules over large opaque models early on.
Add/extend a golden test case whenever a new rule is introduced.
License
MIT

GitHub description (short tagline)
ALK-RECON is a research/education pipeline starter kit that turns ALK-related ctDNA/NGS findings into an auditable mechanism report. This project is not medical advice.