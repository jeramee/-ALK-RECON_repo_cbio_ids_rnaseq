from __future__ import annotations

from schema.case_snapshot import CaseSnapshot, EvidenceItem, EvidenceLevel


DEFAULT_PERSISTER_SCORE_KEY = "persister_score"
DEFAULT_PERSISTER_THRESHOLD = 1.0


def apply_persistence_flags(
    cs: CaseSnapshot,
    score_key: str = DEFAULT_PERSISTER_SCORE_KEY,
    threshold: float = DEFAULT_PERSISTER_THRESHOLD,
) -> None:
    """Set boolean flags derived from expression signature scores.

    If no expression summary exists, this does nothing.
    """
    
    flags = cs.genomic.flags
    flags.setdefault("has_persistence_evidence", False)
    flags.setdefault("persister_signature_score_high", False)
    flags["has_persister_score"] = bool(
    cs.expression and cs.expression.signature_scores and
    ("persister_score" in cs.expression.signature_scores or "PERSISTER_SCORE" in cs.expression.signature_scores)
)


    if cs.expression is None:
        return

    if not cs.expression.signature_scores:
        return

    score = cs.expression.signature_scores.get(score_key)
    if score is None:
        return

    flags["has_persistence_evidence"] = True
    flags["persister_signature_score_high"] = bool(score >= threshold)

    cs.evidence.append(
        EvidenceItem(
            id=f"E_PERSIST_{len(cs.evidence)+1}",
            level=EvidenceLevel.INFERRED,
            kind="expression_signature",
            label=f"Persistence signature score {score_key}={score} (threshold {threshold})",
            value={"score_key": score_key, "score": float(score), "threshold": float(threshold)},
        )
    )
