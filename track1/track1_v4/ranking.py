"""Track 1 V4 ranking corrections over an isolated copy of the V3 engine.

No patient files are read here. Callers supply already locally annotated variants
and phenotype profiles. Numeric utility weights retain their V3 meaning; the
new eligibility fields describe an automated research shortlist, not diagnosis.

Primary references: ClinVar review precedence,
https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/ ; GRCh38 PAR coordinates,
https://www.ncbi.nlm.nih.gov/grc/human .
"""
from __future__ import annotations

from collections import Counter
import importlib.util
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))
_spec = importlib.util.spec_from_file_location("_track1_v4_isolated_v3_ranker", ROOT / "scripts/track1_v3.py")
_v3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_v3)
_parse_condition_assertions = _v3.rcv_evidence
_base_evidence = _v3.evidence

# One-based inclusive GRCh38 coordinates, as used by VCF.
GRCH38_X_PAR = ((10001, 2781479), (155701383, 156030895))


def _classification(assertion):
    if assertion.get("pathogenic_assertion"):
        return "pathogenic"
    if assertion.get("benign"):
        return "benign"
    if "uncertain" in str(assertion.get("classification", "")).lower():
        return "uncertain"
    return "other"


def _review_class(assertion):
    """Expert/guideline precedence; ordinary submissions retain their dissent."""
    stars = assertion["review_stars"]
    return stars if stars >= 3 else 0


def rcv_evidence(variant, disease):
    """Respect higher review classes without deleting lower-class assertions.

    V3 supplies strict condition matching and classification parsing, including
    noncontributing and non-germline exclusions. RCV and SCV are not counted as
    independent evidence; support is still the strongest active assertion.
    """
    parsed = _parse_condition_assertions(variant, disease)
    assertions = parsed["assertions"]
    if not assertions:
        return {**parsed, "active_review_class": None, "lower_review_dissent": False}
    highest = max(map(_review_class, assertions))
    active = [a for a in assertions if _review_class(a) == highest]
    classes = {_classification(a) for a in active}
    conflict = any(a["conflicting"] for a in active)
    conflict |= len(classes & {"pathogenic", "benign", "uncertain"}) > 1
    lower = [a for a in assertions if _review_class(a) < highest]
    dissent = any(a["conflicting"] or _classification(a) not in classes for a in lower)
    supported = [a for a in active if a["pathogenic_assertion"] and not conflict]
    return {
        **parsed,
        "support": max((a["support"] for a in supported), default=0),
        "benign": any(a["benign"] for a in active) and not conflict and not supported,
        "conflict": bool(conflict),
        "assertions": [{**a, "active_review_class_assertion": _review_class(a) == highest}
                       for a in assertions],
        "active_review_class": {0: "ordinary", 3: "expert_panel", 4: "practice_guideline"}[highest],
        "active_assertion_count": len(active),
        "lower_review_dissent": bool(dissent),
        "review_policy": "practice_guideline > expert_panel > ordinary; retain all eligible same-condition assertions",
    }


def _evidence(variant, annotation, disease, config, mode="unknown", **options):
    evidence = _base_evidence(variant, annotation, disease, config, mode, **options)
    if not options.get("no_clinvar", False):
        clinical = rcv_evidence(variant, disease)
        evidence["clinvar_active_review_class"] = clinical["active_review_class"]
        evidence["clinvar_lower_review_dissent"] = clinical["lower_review_dissent"]
        if clinical["lower_review_dissent"]:
            evidence["flags"] = sorted(set(evidence["flags"]) | {"lower_review_ClinVar_dissent_retained"})
    return evidence


# These objects belong only to the separately loaded module. An independently
# imported track1_v3 module is never modified.
_v3.rcv_evidence = rcv_evidence
_v3.evidence = _evidence


def _exact_rules(row, config, key="mechanism_rules"):
    return [rule for rule in config.get(key, [])
            if rule.get("source") and rule.get("gene") == row["gene"]
            and _v3.normalize_condition_id(rule.get("disease_id", "")) == _v3.normalize_condition_id(row["disease_id"])
            and rule.get("inheritance") == row["inheritance"]]


def _in_par(variant, config):
    if _v3.canonical_chrom(variant["chrom"]) != "X":
        return False
    if config.get("assembly", "GRCh38") != "GRCh38":
        return None
    start = int(variant["pos"])
    end = start + len(variant.get("ref", "N")) - 1
    return any(start <= right and end >= left for left, right in GRCH38_X_PAR)


def _inheritance_status(row, members, config):
    mode = row["inheritance"]
    if mode not in {"XL", "XLR", "XLD"}:
        return "not_X_linked", [], []
    contexts = _exact_rules(row, config, "inheritance_context")
    flags, deferred = [], []
    for variant in members:
        par = _in_par(variant, config)
        if par is None:
            deferred.append("X_pseudoautosomal_status_unknown_for_assembly")
        elif par and not any(c.get("par_inheritance_documented") is True for c in contexts):
            deferred.append("pseudoautosomal_variant_requires_documented_inheritance_context")
        gt = variant.get("gt") or []
        dosage = gt.count(variant.get("alt_index", 1))
        observed_haploid_alt = len(gt) == 1 and dosage == 1
        observed_hom_alt = len(gt) == 2 and dosage == 2
        observed_diploid_het = len(gt) == 2 and dosage == 1 and all(a is not None for a in gt)
        if not (observed_haploid_alt or observed_hom_alt or observed_diploid_het):
            deferred.append("X_allele_ploidy_or_dosage_unresolved")
        if mode == "XL":
            deferred.append("X_linked_inheritance_subtype_unresolved")
        if mode == "XLR" and observed_diploid_het:
            documented = any(c.get("manifesting_heterozygote_compatible") is True for c in contexts)
            flags.append("diploid_XLR_heterozygote_may_manifest_context_required")
            if not documented:
                deferred.append("diploid_XLR_heterozygote_expression_context_unresolved")
    flags.append("observed_genotype_does_not_establish_sex_or_parental_origin")
    return ("compatibility_unresolved" if deferred else "observed_genotype_compatible"), flags, deferred


def _mechanism_status(row, members, config):
    rules = [r for r in _exact_rules(row, config) if r.get("mechanism")]
    if not rules:
        return "no_exact_curated_rule", [], []
    flags, deferred, status = [], [], []
    for variant, evidence in zip(members, row["allele_evidence"]):
        annotations = [a for a in variant["annotations"] if a["gene"] == row["gene"]
                       and a.get("transcript") == evidence["selected_transcript"]]
        compatible = []
        for ann in annotations:
            for rule in rules:
                allowed = set(rule.get("supported_consequences", [])) | set(rule.get("plausible_consequences", []))
                # An unspecified list provides no negative compatibility evidence.
                if not allowed:
                    continue
                transcripts = {str(t).split(".")[0] for t in rule.get("transcripts", [])}
                transcript_ok = not transcripts or str(ann.get("transcript", "")).split(".")[0] in transcripts
                compatible.append(bool(set(ann.get("consequences", [])) & allowed) and transcript_ok)
        if compatible and not any(compatible):
            flags.append("consequence_compatibility_with_exact_curated_mechanism_unknown")
            if evidence["clinical_support"] > 0:
                status.append("condition_supported_despite_mechanism_review")
                flags.append("condition_specific_assertion_retained_despite_mechanism_mismatch")
            else:
                status.append("compatibility_unresolved")
                deferred.append("curated_mechanism_compatibility_unresolved_without_condition_support")
        else:
            status.append("consequence_class_plausible_not_functionally_confirmed" if compatible else "curated_consequence_scope_unspecified")
    return ";".join(sorted(set(status))), flags, deferred


def _pair_transcript_status(row, members, config):
    if len(members) < 2:
        return {"status": "single_allele", "common_relevant_transcripts": []}, [], []
    threshold = config.get("minimum_allele_utility", .3)
    rules = _exact_rules(row, config)
    relevant_sets = []
    for variant, evidence in zip(members, row["allele_evidence"]):
        transcripts = set()
        for ann in variant["annotations"]:
            if ann["gene"] != row["gene"] or not ann.get("transcript"):
                continue
            transcript = str(ann["transcript"]).split(".")[0]
            if rules and any(r.get("transcripts") for r in rules):
                permitted = set().union(*({str(t).split(".")[0] for t in r.get("transcripts", [])} for r in rules))
                if transcript not in permitted:
                    continue
            plausible_rule = any(set(ann.get("consequences", [])) &
                                 (set(r.get("supported_consequences", [])) | set(r.get("plausible_consequences", [])))
                                 for r in rules)
            splice_review = evidence.get("splice_prediction_rescue") and ann.get("transcript") == evidence.get("selected_transcript")
            if float(ann.get("effect_score", 0)) >= threshold or plausible_rule or splice_review:
                transcripts.add(transcript)
        relevant_sets.append(transcripts)
    common = set.intersection(*relevant_sets) if all(relevant_sets) else set()
    evidence = {"common_relevant_transcripts": sorted(common),
                "relevant_transcripts_by_member": [sorted(s) for s in relevant_sets],
                "same_transcript_does_not_establish_trans_phase": True}
    if not all(relevant_sets):
        return {**evidence, "status": "unknown"}, [], ["pair_relevant_transcript_context_unresolved"]
    if not common:
        return {**evidence, "status": "disjoint_available_relevant_transcripts"}, [], ["pair_available_relevant_transcripts_disjoint"]
    selected = {str(t).split(".")[0] for t in row["selected_transcripts"]}
    flags = [] if len(selected) == 1 else ["selected_transcripts_differ_but_common_relevant_annotation_exists"]
    return {**evidence, "status": "common_relevant_annotation_available_not_functional_confirmation"}, flags, []


def rank(variants, profiles, config, **options):
    """Return automatic hypotheses, deferred hypotheses, and operational counts.

    Input profiles must preserve XLR/XLD/XL rather than collapsing subtypes.
    Sourced optional ``inheritance_context`` rules use exact gene, disease_id and
    inheritance, plus ``manifesting_heterozygote_compatible`` and/or
    ``par_inheritance_documented``; neither value is inferred from the genome.
    """
    variants = list(variants)
    lookup = {variant["id"]: variant for variant in variants}
    if len(lookup) != len(variants):
        raise ValueError("Ranking requires normalized, unique variant IDs")
    for diseases in profiles.values():
        for disease in diseases:
            score = disease["score"]
            if not isinstance(score, (int, float)) or not math.isfinite(score) or not 0 <= score <= 1:
                raise ValueError("Phenotype scores must be finite and within zero to one")
    rows, base_counts = _v3.model(variants, profiles, config, **options)
    counts = Counter(base_counts)
    automatic, review = [], []
    for row in rows:
        members = [lookup[identifier] for identifier in row["members"]]
        flags, deferred = set(row["flags"]), []
        inheritance, extra, reasons = _inheritance_status(row, members, config)
        flags.update(extra); deferred.extend(reasons)
        mechanism, extra, reasons = _mechanism_status(row, members, config)
        flags.update(extra); deferred.extend(reasons)
        transcript, extra, reasons = _pair_transcript_status(row, members, config)
        flags.update(extra); deferred.extend(reasons)
        if not row["eligible_for_manual_selection"]:
            deferred.append("incomplete_or_unknown_inheritance_hypothesis")
        eligible = not deferred
        row.update({"combined_review_rank": row["rank"], "flags": sorted(flags | set(deferred)),
                    "inheritance_compatibility": inheritance, "mechanism_compatibility": mechanism,
                    "pair_transcript_evidence": transcript,
                    "automatic_shortlist_eligible": eligible,
                    "automatic_shortlist_deferral_reasons": sorted(set(deferred)),
                    "clinical_confirmation": False,
                    "review_status": "automated_research_candidate" if eligible else "deferred_for_context_review"})
        # Automated submission eligibility is distinct from clinical confirmation
        # and from the V3 optional human-review decision. Nothing attests review.
        row["submission_eligible"] = False
        target = automatic if eligible else review
        row["rank"] = len(target) + 1
        target.append(row)
        for reason in set(deferred):
            counts["v4_deferred_" + reason] += 1
    counts["automatic_shortlist_hypotheses"] = len(automatic)
    counts["deferred_review_hypotheses"] = len(review)
    return automatic, review, dict(counts)


def distinct_rows(rows, limit=None):
    """Preserve ordering while selecting each complete allele set at most once."""
    if limit is not None and (not isinstance(limit, int) or limit < 0):
        raise ValueError("Limit must be a nonnegative integer or None")
    if limit == 0:
        return []
    selected, seen = [], set()
    for row in rows:
        key = tuple(sorted(row["members"]))
        if key in seen:
            continue
        seen.add(key)
        selected.append(row)
        if limit is not None and len(selected) >= limit:
            break
    return selected if limit is None else selected[:limit]
