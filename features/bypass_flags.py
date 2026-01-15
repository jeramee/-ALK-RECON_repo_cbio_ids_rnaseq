from __future__ import annotations

from schema.case_snapshot import CaseSnapshot, EvidenceItem, EvidenceLevel


# A small starter list of bypass-related genes/labels. Tune as you learn.
BYPASS_GENES = {
    "MET",
    "EGFR",
    "KRAS",
    "NRAS",
    "BRAF",
    "MAP2K1",  # MEK1
    "MAP2K2",  # MEK2
    "ERBB2",
    "PIK3CA",
    "PTEN",
}


def apply_bypass_flags(cs: CaseSnapshot) -> CaseSnapshot:
    """Derive simple bypass flags from copy-number and bypass_events tables."""
    flags = cs.genomic.flags

    # Normalize events to uppercase gene keys where possible
    genes_present: set[str] = set()
    for ev in cs.genomic.bypass_events or []:
        gene = str(ev.get("gene") or ev.get("GENE") or ev.get("Hugo_Symbol") or "").upper()
        if gene:
            genes_present.add(gene)
    for cn in cs.genomic.copy_number_events or []:
        gene = str(cn.get("gene") or cn.get("GENE") or cn.get("Hugo_Symbol") or "").upper()
        if gene:
            genes_present.add(gene)

    flags["has_MET_amp"] = False
    if cs.genomic and cs.genomic.bypass_events:
        for ev in cs.genomic.bypass_events:
            gene = (ev.get("gene") or "").upper()
            et = (ev.get("event_type") or ev.get("type") or "").upper()
            if gene == "MET" and ("AMP" in et or "AMPL" in et):
                flags["has_MET_amp"] = True
                break


    # Simple starter heuristics
    flags["has_bypass_gene_event"] = any(g in BYPASS_GENES for g in genes_present)
    flags["has_MET_event"] = "MET" in genes_present
    flags["has_EGFR_event"] = "EGFR" in genes_present
    flags["has_MAPK_pathway_event"] = any(g in {"KRAS", "NRAS", "BRAF", "MAP2K1", "MAP2K2"} for g in genes_present)

    # Add evidence for derived flags
    cs.evidence.append(EvidenceItem(
        id=f"F_bypass_any",
        level=EvidenceLevel.INFERRED,
        kind="feature_flag",
        label="Bypass-related genomic event present (starter heuristic)",
        value={k: flags.get(k) for k in ["has_bypass_gene_event", "has_MET_event", "has_EGFR_event", "has_MAPK_pathway_event"]},
    ))

    return cs

