from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import json

from schema.case_snapshot import save_json_schema

from ingest.variant_table_import import (
    load_variant_table,
    dataframe_to_case_snapshots,
)
from ingest.rnaseq_import import build_expression_index
from features.apply_features import apply_all_features
from mechanism_engine.rule_engine import score_mechanisms_and_route
from reports.dossier import write_dossier_bundle


def run_pipeline(
    input_path: Path | str,
    out_dir: Path | str,
    delimiter: str | None = None,
    case_col: str | None = None,
    sample_col: str | None = None,
    timepoint_col: str | None = None,
    min_vaf: float | None = None,
    rnaseq_counts_path: Path | str | None = None,
    rnaseq_metadata_path: Path | str | None = None,
    rnaseq_gene_col: str | None = None,
    rnaseq_sample_id_col: str | None = None,
    rnaseq_case_id_col: str | None = None,
    rnaseq_study_id_col: str | None = None,
    rnaseq_timepoint_id_col: str | None = None,
    rnaseq_signature_path: Path | str | None = None,
    write_json_schema: bool = False,
) -> list[dict]:
    """End-to-end run: import -> CaseSnapshot -> features -> mechanism -> reports.

    Returns list of output dicts (one per snapshot).
    """
    input_path = Path(input_path)
    out_dir = Path(out_dir)

    # Optional RNA-seq module: create an index of ExpressionSummary objects keyed by
    # (study_id, case_id, sample_id, timepoint_id). We attach these to snapshots later.
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
        if expr_index is not None:
            key = (
                (snap.study_id or "unspecified"),
                (snap.case_id or ""),
                (snap.sample_id or ""),
                (snap.timepoint_id or "unspecified"),
            )
            # Attach expression if available. If not found, we leave it as None.
            snap.expression = expr_index.get(key)

        apply_all_features(snap)
        score_mechanisms_and_route(snap)
        bundle = write_dossier_bundle(snap, out_dir)
        outputs.append(bundle)

    if write_json_schema:
        save_json_schema(str(out_dir / "case_snapshot.schema.json"))

    # Also write an index file for convenience
    idx_path = out_dir / "index.json"
    with idx_path.open("w", encoding="utf-8") as f:
        json.dump(outputs, f, indent=2, ensure_ascii=False)

    return outputs

