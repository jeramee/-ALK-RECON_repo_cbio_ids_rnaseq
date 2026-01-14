from __future__ import annotations

import argparse
from pathlib import Path

from .pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="alk-recon",
        description="Generate ALK mechanism dossiers from a ctDNA/NGS variant table (research/education only).",
    )
    p.add_argument(
        "--input",
        required=True,
        help="Path to ctDNA/NGS variant table (CSV or TSV).",
    )
    p.add_argument(
        "--outdir",
        required=True,
        help="Output directory. Will be created if missing.",
    )
    p.add_argument(
        "--delimiter",
        default=None,
        help="Optional delimiter override. Default: auto-detect from file extension.",
    )
    p.add_argument(
        "--case-col",
        default=None,
        help="Optional: name of the case/patient column if auto-detection fails.",
    )
    p.add_argument(
        "--sample-col",
        default=None,
        help="Optional: name of the sample column if auto-detection fails.",
    )
    p.add_argument(
        "--timepoint-col",
        default=None,
        help="Optional: name of the timepoint column if present.",
    )
    p.add_argument(
        "--min-vaf",
        type=float,
        default=None,
        help="Optional: minimum VAF threshold (0..1) applied when a VAF column exists.",
    )
    p.add_argument(
        "--write-json-schema",
        action="store_true",
        help="Also write schema/case_snapshot.schema.json into the output directory.",
    )

    # Optional: RNA-seq counts + metadata (persistence signatures)
    p.add_argument(
        "--rnaseq-counts",
        default=None,
        help="Optional: RNA-seq count matrix (CSV/TSV). Rows=genes, columns=samples (or vice versa; see docs).",
    )
    p.add_argument(
        "--rnaseq-meta",
        default=None,
        help="Optional: RNA-seq sample metadata (CSV/TSV) with sample_id and optional case_id/study_id/timepoint_id.",
    )
    p.add_argument(
        "--rnaseq-gene-col",
        default=None,
        help="Optional: gene column name for the counts matrix (default: auto-detect first column / 'gene').",
    )
    p.add_argument(
        "--rnaseq-sample-id-col",
        default=None,
        help="Optional: sample_id column in the RNA-seq metadata (default: auto-detect).",
    )
    p.add_argument(
        "--rnaseq-case-id-col",
        default=None,
        help="Optional: case_id column in the RNA-seq metadata (default: auto-detect).",
    )
    p.add_argument(
        "--rnaseq-study-id-col",
        default=None,
        help="Optional: study_id column in the RNA-seq metadata (default: auto-detect).",
    )
    p.add_argument(
        "--rnaseq-timepoint-id-col",
        default=None,
        help="Optional: timepoint_id column in the RNA-seq metadata (default: auto-detect).",
    )
    p.add_argument(
        "--rnaseq-signature",
        default=None,
        help="Optional: path to a newline-delimited gene list used to compute a simple 'persister_score'.",
    )

    # Optional: Copy-number matrices (cBioPortal bundles)
    p.add_argument(
        "--cna-thresholded",
        default=None,
        help="Optional: thresholded CNA matrix (e.g., cBioPortal data_cna.txt).",
    )
    p.add_argument(
        "--cna-linear",
        default=None,
        help="Optional: linear CNA matrix (e.g., cBioPortal data_linear_cna.txt).",
    )
    p.add_argument(
        "--cna-meta",
        default=None,
        help=(
            "Optional: CNA sample metadata mapping (sample_id, optional case_id/study_id/timepoint_id). "
            "If omitted, --rnaseq-meta will be reused."
        ),
    )
    p.add_argument(
        "--cna-study-id-col",
        default=None,
        help="Optional: study_id column in CNA metadata (default: auto-detect).",
    )
    p.add_argument(
        "--cna-timepoint-id-col",
        default=None,
        help="Optional: timepoint_id column in CNA metadata (default: auto-detect).",
    )
    p.add_argument(
        "--cna-genes",
        default=None,
        help="Optional: comma-separated genes to extract from CNA matrices (default: ALK,MET,EGFR,ERBB2,KRAS,BRAF,RET,ROS1).",
    )
    return p


def main() -> None:
    args = build_parser().parse_args()

    in_path = Path(args.input).expanduser().resolve()
    out_dir = Path(args.outdir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    run_pipeline(
        input_path=in_path,
        out_dir=out_dir,
        delimiter=args.delimiter,
        case_col=args.case_col,
        sample_col=args.sample_col,
        timepoint_col=args.timepoint_col,
        min_vaf=args.min_vaf,
        rnaseq_counts_path=args.rnaseq_counts,
        rnaseq_metadata_path=args.rnaseq_meta,
        rnaseq_gene_col=args.rnaseq_gene_col,
        rnaseq_sample_id_col=args.rnaseq_sample_id_col,
        rnaseq_case_id_col=args.rnaseq_case_id_col,
        rnaseq_study_id_col=args.rnaseq_study_id_col,
        rnaseq_timepoint_id_col=args.rnaseq_timepoint_id_col,
        rnaseq_signature_path=args.rnaseq_signature,
        cna_thresholded_path=args.cna_thresholded,
        cna_linear_path=args.cna_linear,
        cna_metadata_path=args.cna_meta,
        cna_study_id_col=args.cna_study_id_col,
        cna_timepoint_id_col=args.cna_timepoint_id_col,
        cna_genes=args.cna_genes,
        write_json_schema=args.write_json_schema,
    )


if __name__ == "__main__":
    main()
