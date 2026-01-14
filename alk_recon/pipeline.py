from __future__ import annotations

from pathlib import Path
import json

from schema.case_snapshot import save_json_schema

from ingest.variant_table_import import (
    load_variant_table,
    dataframe_to_case_snapshots,
)
from ingest.rnaseq_import import build_expression_index
from ingest.cna_import import build_cna_index
from features.apply_features import apply_all_features
from mechanism_engine.rule_engine import score_mechanisms_and_route
from reports.dossier import write_dossier_bundle


def _candidate_keys(
    study_id: str | None,
    case_id: str | None,
    sample_id: str | None,
    timepoint_id: str | None,
) -> list[tuple[str | None, str | None, str, str | None]]:
    """
    Generate fallback keys for attaching aux data (RNA/CNA).

    TCGA files often disagree on identifiers:
      - TCGA-XX-YYYY-01 (aliquot/sample)
      - TCGA-XX-YYYY     (patient/case; first 12 chars)
      - TCGA.XX.YYYY     (dot form)

    We try a few deterministic fallbacks.
    """
    if not sample_id:
        return []

    case_id = case_id or ""

    def norm12(x: str) -> str:
        return str(x).replace(".", "-")[:12]

    st = study_id or "unspecified"
    tp = timepoint_id

    out: list[tuple[str | None, str | None, str, str | None]] = [(st, case_id, sample_id, tp)]

    s12 = norm12(sample_id)
    c12 = norm12(case_id)

    # patient-id fallbacks
    if s12 != sample_id:
        out.append((st, case_id, s12, tp))
    if c12 != case_id:
        out.append((st, c12, s12, tp))

    # timepoint fallbacks
    for tpf in ("", "unspecified", None):
        out.append((st, case_id, sample_id, tpf))
        out.append((st, c12, s12, tpf))

    # study fallbacks
    for stf in ("unspecified", None):
        out.append((stf, case_id, sample_id, tp))
        out.append((stf, c12, s12, tp))

    # de-dup preserve order
    seen: set[tuple[str | None, str | None, str, str | None]] = set()
    uniq: list[tuple[str | None, str | None, str, str | None]] = []
    for k in out:
        if k not in seen:
            uniq.append(k)
            seen.add(k)
    return uniq


def run_pipeline(
    input_path: Path | str,
    out_dir: Path | str | None = None,
    # ---- back-compat alias (older tests/callers) ----
    outdir: Path | str | None = None,
    delimiter: str | None = None,
    case_col: str | None = None,
    sample_col: str | None = None,
    timepoint_col: str | None = None,
    min_vaf: float | None = None,
    rnaseq_counts_path: Path | str | None = None,
    rnaseq_metadata_path: Path | str | None = None,
    # ---- back-compat alias ----
    rnaseq_meta_path: Path | str | None = None,
    rnaseq_gene_col: str | None = None,
    rnaseq_sample_id_col: str | None = None,
    rnaseq_case_id_col: str | None = None,
    rnaseq_study_id_col: str | None = None,
    rnaseq_timepoint_id_col: str | None = None,
    rnaseq_signature_path: Path | str | None = None,
    cna_thresholded_path: Path | str | None = None,
    cna_linear_path: Path | str | None = None,
    cna_metadata_path: Path | str | None = None,
    cna_sample_id_col: str | None = None,
    cna_case_id_col: str | None = None,
    cna_study_id_col: str | None = None,
    cna_timepoint_id_col: str | None = None,
    cna_genes: str | None = None,
    # ---- back-compat alias ----
    cna_gene_whitelist: str | None = None,
    # optional tuning (tests pass this in some branches)
    persister_threshold: float | None = None,
    write_json_schema: bool = False,
) -> list[dict]:
    """End-to-end run: import -> CaseSnapshot -> features -> mechanism -> reports."""
    # resolve back-compat aliases
    if out_dir is None:
        out_dir = outdir
    if rnaseq_metadata_path is None:
        rnaseq_metadata_path = rnaseq_meta_path
    if cna_genes is None:
        cna_genes = cna_gene_whitelist

    if out_dir is None:
        raise ValueError("out_dir is required (or provide outdir=...).")

    input_path = Path(input_path)
    out_dir = Path(out_dir)

    # ---- Optional RNA-seq index ----
    expr_index = None
    if rnaseq_counts_path and rnaseq_metadata_path:
        expr_index = build_expression_index(
            counts_path=Path(rnaseq_counts_path),
            metadata_path=Path(rnaseq_metadata_path),
            gene_col=rnaseq_gene_col,
            sample_id_col=rnaseq_sample_id_col,
            case_id_col=rnaseq_case_id_col,
            study_id_col=rnaseq_study_id_col,
            timepoint_id_col=rnaseq_timepoint_id_col,
            signature_path=Path(rnaseq_signature_path) if rnaseq_signature_path else None,
        )

    # ---- Optional CNA index ----
    cna_index = None
    if cna_thresholded_path or cna_linear_path:
        meta_path = cna_metadata_path or rnaseq_metadata_path
        if not meta_path:
            raise ValueError(
                "CNA ingest requested but no metadata path provided. "
                "Provide --cna-meta or reuse --rnaseq-meta."
            )

        genes = None
        if cna_genes:
            genes = {g.strip() for g in cna_genes.split(",") if g.strip()}

        cna_index = build_cna_index(
            metadata_path=Path(meta_path),
            thresholded_path=Path(cna_thresholded_path) if cna_thresholded_path else None,
            linear_path=Path(cna_linear_path) if cna_linear_path else None,
            sample_id_col=cna_sample_id_col,
            case_id_col=cna_case_id_col,
            study_id_col=cna_study_id_col,
            timepoint_id_col=cna_timepoint_id_col,
            genes_of_interest=genes,
        )

    # ---- Load variants and build snapshots ----
    df = load_variant_table(
        input_path,
        delimiter=delimiter,
        case_col=case_col,
        sample_col=sample_col,
        timepoint_col=timepoint_col,
        min_vaf=min_vaf,
    )
    snapshots = dataframe_to_case_snapshots(df, source_id=str(input_path))

    outputs: list[dict] = []
    for snap in snapshots:
        # attach expression
        if expr_index is not None:
            snap.expression = None
            for key in _candidate_keys(snap.study_id, snap.case_id, snap.sample_id, snap.timepoint_id):
                if key in expr_index:
                    snap.expression = expr_index[key]
                    break

        # attach CNA
        if cna_index is not None and getattr(snap, "genomic", None) is not None:
            existing = list(getattr(snap.genomic, "copy_number_events", []) or [])
            added = None
            for key in _candidate_keys(snap.study_id, snap.case_id, snap.sample_id, snap.timepoint_id):
                if key in cna_index:
                    added = cna_index[key]
                    break
            if added:
                snap.genomic.copy_number_events = existing + list(added)

        # features + routing
        # (persister_threshold is optional and only used if your feature fn supports it)
        try:
            apply_all_features(snap, persister_threshold=persister_threshold)
        except TypeError:
            apply_all_features(snap)

        score_mechanisms_and_route(snap)
        bundle = write_dossier_bundle(snap, out_dir)
        outputs.append(bundle)

    if write_json_schema:
        save_json_schema(str(out_dir / "case_snapshot.schema.json"))

    # index.json for convenience
    idx_path = out_dir / "index.json"
    with idx_path.open("w", encoding="utf-8") as f:
        json.dump(outputs, f, indent=2, ensure_ascii=False)

    return outputs
