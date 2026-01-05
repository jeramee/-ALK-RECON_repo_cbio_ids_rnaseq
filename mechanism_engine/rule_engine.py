from __future__ import annotations

from dataclasses import replace

from schema.case_snapshot import (
    CaseSnapshot,
    MechanismCall,
    MechanismType,
    StrategyBucket,
    StrategyRouting,
)


def score_mechanisms_and_route(cs: CaseSnapshot) -> CaseSnapshot:
    """Score mechanisms (rule-based) and compute strategy routing.

    This is intentionally interpretable and auditable for early MVP.
    You can later replace scoring with an ML model, but keep the same
    inputs/outputs.
    """

    flags = cs.genomic.flags or {}

    # --- Mechanism scoring (0..1, heuristic) ---
    on_target = 0.0
    bypass = 0.0
    persistence = 0.0

    # On-target ALK evidence
    if flags.get("has_any_alk_mutation"):
        on_target += 0.40
    if flags.get("has_G1202R"):
        on_target += 0.35
    if flags.get("has_L1196M"):
        on_target += 0.25
    if flags.get("has_compound_alk_mutations"):
        on_target += 0.35

    # Bypass evidence
    if flags.get("has_any_bypass_event"):
        bypass += 0.35
    if flags.get("has_MET_event"):
        bypass += 0.30
    if flags.get("has_EGFR_event"):
        bypass += 0.25
    if flags.get("has_MAPK_event"):
        bypass += 0.20

    # Persistence evidence
    if flags.get("has_expression_data"):
        persistence += 0.15
    if flags.get("has_persister_score_high"):
        persistence += 0.60

    # Normalize into 0..1 by simple cap
    on_target = min(1.0, on_target)
    bypass = min(1.0, bypass)
    persistence = min(1.0, persistence)

    # Build calls (ranked)
    calls = []
    calls.append(
        MechanismCall(
            mechanism=MechanismType.ON_TARGET_ALK,
            score=on_target,
            rationale=_rationale_on_target(cs),
            supporting_evidence_ids=_supporting_evidence_ids(cs, "on_target_alk"),
            contradicting_evidence_ids=[],
        )
    )
    calls.append(
        MechanismCall(
            mechanism=MechanismType.BYPASS,
            score=bypass,
            rationale=_rationale_bypass(cs),
            supporting_evidence_ids=_supporting_evidence_ids(cs, "bypass"),
            contradicting_evidence_ids=[],
        )
    )
    calls.append(
        MechanismCall(
            mechanism=MechanismType.PERSISTENCE,
            score=persistence,
            rationale=_rationale_persistence(cs),
            supporting_evidence_ids=_supporting_evidence_ids(cs, "persistence"),
            contradicting_evidence_ids=[],
        )
    )

    calls_sorted = sorted(calls, key=lambda c: c.score, reverse=True)
    cs.mechanism_calls = calls_sorted

    cs.routing = _route(cs)
    return cs


def _route(cs: CaseSnapshot) -> StrategyRouting:
    flags = cs.genomic.flags or {}
    top = cs.mechanism_calls[0] if cs.mechanism_calls else None

    ranked: list[StrategyBucket] = []
    what_to_test: list[str] = []
    what_to_avoid: list[str] = []

    # Some stable guardrails (keeps reports honest)
    what_to_avoid.append(
        "Do not treat this as a clinical recommendation; this is a research summary of evidence." 
    )

    # If we have on-target mutations, prioritize mutation-aware + degradation
    if flags.get("has_any_alk_mutation"):
        ranked.append(StrategyBucket.A_ATP_MUTATION_AWARE)
        ranked.append(StrategyBucket.C_DEGRADATION)
        what_to_test.append("Confirm effect across a WT + mutant panel (include solvent-front, gatekeeper, compound if possible).")

        if flags.get("has_compound_alk_mutations"):
            what_to_test.append("If compound mutations exist, explicitly test compound constructs; single-mutant results can mislead.")

    # If bypass features, add bypass-aware logic
    if flags.get("has_any_bypass_event"):
        ranked.append(StrategyBucket.E_SEQUENCING_LOGIC)
        what_to_test.append("Validate bypass signal (e.g., MET/EGFR/MAPK) with orthogonal evidence when available.")

    # If persistence evidence, include persistence adjunct bucket
    if flags.get("has_persister_score_high"):
        ranked.append(StrategyBucket.D_PERSISTENCE_ADJUNCT)
        what_to_test.append("Run a durability assay (days–weeks) to quantify a persister fraction and rebound kinetics.")

    # If weak evidence overall, emphasize sequencing/monitoring logic
    if top is not None and top.score < 0.35:
        ranked = [StrategyBucket.E_SEQUENCING_LOGIC]
        what_to_test = [
            "Acquire higher-yield evidence: ctDNA/tissue NGS for ALK mutations and key bypass events; optional expression signatures if available."
        ]

    # Deduplicate while preserving order
    ranked_unique: list[StrategyBucket] = []
    for b in ranked:
        if b not in ranked_unique:
            ranked_unique.append(b)

    return StrategyRouting(
        ranked_buckets=ranked_unique,
        what_to_test_next=what_to_test,
        what_to_avoid=what_to_avoid,
    )


def _rationale_on_target(cs: CaseSnapshot) -> str:
    f = cs.genomic.flags or {}
    if not f.get("has_any_alk_mutation"):
        return "No ALK kinase-domain variant was detected in the provided variant table."
    bits = []
    if f.get("has_G1202R"):
        bits.append("solvent-front pattern (G1202R)")
    if f.get("has_L1196M"):
        bits.append("gatekeeper pattern (L1196M)")
    if f.get("has_compound_alk_mutations"):
        bits.append("multiple ALK variants (possible compound context)")
    if not bits:
        bits.append("ALK variant(s) detected")
    return "On-target ALK resistance is supported by: " + ", ".join(bits) + "."


def _rationale_bypass(cs: CaseSnapshot) -> str:
    f = cs.genomic.flags or {}
    if not f.get("has_any_bypass_event"):
        return "No obvious bypass-related alterations were detected in the provided variant table."
    bits = []
    if f.get("has_MET_event"):
        bits.append("MET-related alteration")
    if f.get("has_EGFR_event"):
        bits.append("EGFR/ERBB-family alteration")
    if f.get("has_MAPK_event"):
        bits.append("MAPK-pathway alteration")
    if not bits:
        bits.append("bypass-related alteration")
    return "Bypass signaling is supported by: " + ", ".join(bits) + "."


def _rationale_persistence(cs: CaseSnapshot) -> str:
    f = cs.genomic.flags or {}
    if not f.get("has_expression_data"):
        return "No expression-level persistence evidence was provided (expression block is missing)."
    if f.get("has_persister_score_high"):
        return "Expression signature score suggests a possible drug-tolerant persister program (threshold exceeded)."
    return "Expression data was provided, but no persistence signature crossed the current threshold."


def _supporting_evidence_ids(cs: CaseSnapshot, mechanism: str) -> list[str]:
    """Heuristic mapping from flags to evidence items.

    For the MVP we keep this simple: any evidence label containing a key term
    is considered supporting. Later you can make this explicit and exact.
    """
    mech_terms = {
        "on_target_alk": ["ALK", "G1202R", "L1196M"],
        "bypass": ["MET", "EGFR", "ERBB", "KRAS", "NRAS", "BRAF", "MAP2K", "MAPK"],
        "persistence": ["persister", "signature", "expression"],
    }
    terms = mech_terms.get(mechanism, [])
    out: list[str] = []
    for ev in cs.evidence:
        text = f"{ev.label} {ev.value}".upper()
        if any(t.upper() in text for t in terms):
            out.append(ev.id)
    return out
