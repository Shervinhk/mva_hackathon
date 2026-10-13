"""Pure local rendering of saved Track 1 evidence; no file or network access.

Returned text is patient-derived when real evidence is supplied. Save it locally;
do not print it into an assistant session. This renderer neither reranks nor
estimates causality. Missing fields remain explicitly unavailable.
"""
from __future__ import annotations

import html
import math


def _text(value):
    text = html.escape(" ".join(str(value).split()), quote=False)
    for token in ("\\", "`", "*", "_", "[", "]", "|"):
        text = text.replace(token, "\\" + token)
    return text


def _number(value, digits=3):
    if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value):
        return f"{value:.{digits}g}"
    return "unavailable"


def _key(row):
    return (row.get("gene"), row.get("disease_id"), row.get("inheritance"),
            tuple(sorted(row.get("members", []))))


def _phenotype(row, profiles):
    matching = [p for p in profiles.get(row.get("gene"), [])
                if p.get("disease_id") == row.get("disease_id")
                and row.get("inheritance") in p.get("modes", [row.get("inheritance")])]
    score = row.get("components", {}).get("phenotype")
    lines = [f"Phenotype similarity is {_number(score)} on the model's 0–1 scale; it is not a causal probability."]
    if len(matching) != 1:
        lines.append("An unambiguous term-level disease profile was not supplied to this report.")
        return " ".join(lines)
    profile = matching[0]
    terms = profile.get("terms", [])
    scored = [t for t in terms if isinstance(t, dict) and t.get("scored") is True]
    exact = [t for t in scored if t.get("term") == t.get("matched_ancestor") and t.get("matched_information", 0) > 0]
    related = [t for t in scored if t.get("term") != t.get("matched_ancestor") and t.get("matched_information", 0) > 0]
    if scored:
        lines.append(f"Among {len(scored)} scored positive observations, {len(exact)} match at the queried HPO term "
                     f"and {len(related)} have support only through a shared ancestor term in this disease profile.")
        if exact:
            lines.append("Query-term match examples: " + ", ".join(_text(t["term"]) for t in exact[:3]) + ". "
                         "These matches can also arise from a more specific annotated disease term.")
        elif related:
            examples = [f"{_text(t['term'])} → {_text(t['matched_ancestor'])}" for t in related[:3]]
            lines.append("Shared-ancestor examples: " + "; ".join(examples) + ".")
        if len(scored) > len(exact) + len(related):
            lines.append("The remaining scored observations contribute no positive matching information to this profile; "
                         "this does not establish that the disease is excluded.")
    else:
        lines.append("No scored term-level matches were supplied to this report.")
    unscored = profile.get("unscored_terms", [])
    if unscored:
        lines.append(f"{len(unscored)} query term(s) were excluded from the similarity calculation as unscored.")
    return " ".join(lines)


def _allele(identifier, variant, evidence):
    pieces = [f"{_text(identifier)}"]
    if not variant:
        return pieces[0] + ": variant-level details were not supplied."
    genotype = variant.get("gt")
    if isinstance(genotype, (list, tuple)) and genotype:
        # Slash separates observed allele indexes without inventing phase.
        pieces.append("observed allele indexes " + "/".join("." if a is None else str(a) for a in genotype))
    if variant.get("quality_pass") is True:
        pieces.append("passed the recorded call-quality filters")
    elif variant.get("quality_pass") is False:
        pieces.append("did not pass the recorded call-quality filters")
    else:
        pieces.append("call-quality status unavailable")
    transcript = evidence.get("selected_transcript")
    annotations = [a for a in variant.get("annotations", []) if a.get("transcript") == transcript]
    consequences = sorted({c for a in annotations for c in a.get("consequences", [])})
    if transcript and consequences:
        pieces.append("annotated " + "/".join(_text(c) for c in consequences) + " on " + _text(transcript))
    elif transcript:
        pieces.append("selected transcript " + _text(transcript) + "; consequence details unavailable")
    population = variant.get("population", {})
    if population.get("known") is True:
        pieces.append("recorded maximum usable population allele frequency " + _number(population.get("af"), 4))
        if population.get("all_measurement_quality_verified") is not True:
            pieces.append("per-allele population filter verification is incomplete")
    else:
        pieces.append("population frequency unknown or unusable, not established as absent")
    support = evidence.get("clinical_support")
    if isinstance(support, (int, float)) and support > 0:
        pieces.append("matching condition-specific pathogenic/likely pathogenic assertion support is recorded")
        review_class = evidence.get("clinvar_active_review_class")
        if review_class:
            pieces.append("active ClinVar review class " + _text(review_class).replace("\\_", " "))
    elif support == 0:
        pieces.append("no positive condition-specific ClinVar assertion support contributes to this score")
    else:
        pieces.append("condition-specific ClinVar support status unavailable")
    if evidence.get("clinical_conflict"):
        pieces.append("active-level ClinVar disagreement remains unresolved")
    if evidence.get("clinvar_lower_review_dissent"):
        pieces.append("lower-review ClinVar dissent is retained")
    return "; ".join(pieces) + "."


def _uncertainty(row):
    phase = row.get("phase")
    sentences = []
    if phase == "unknown":
        sentences.append("Phase is unresolved: the two alleles are not established to be on opposite chromosome copies, "
                         "so confirmed compound heterozygosity is not claimed.")
    elif phase == "trans":
        sentences.append("The saved source-record phase evidence supports trans; this report does not add parental or experimental validation.")
    elif phase == "homozygous":
        sentences.append("The saved hypothesis is homozygous based on the called diploid genotype.")
    elif phase == "not_applicable":
        sentences.append("This is a single-allele hypothesis; de novo status and parental origin are not established by this ranking.")
    else:
        sentences.append("Saved phase category: " + _text(phase if phase is not None else "unavailable") + ".")
    if row.get("phase_evidence"):
        sentences.append("Recorded phase basis: " + _text(row["phase_evidence"]) + ".")
    pair = row.get("pair_transcript_evidence", {})
    common = pair.get("common_relevant_transcripts", [])
    if len(row.get("members", [])) == 2 and common:
        sentences.append("A shared plausibly affected transcript annotation is available (" +
                         ", ".join(_text(t) for t in common[:3]) + "); this does not establish phase or protein function.")
    if "selected_transcripts_differ_but_common_relevant_annotation_exists" in row.get("flags", []):
        sentences.append("The individually selected highest-scoring transcripts differ, so the shared-transcript interpretation remains a review point.")
    mechanisms = []
    for evidence in row.get("allele_evidence", []):
        detail = evidence.get("mechanism_evidence", {})
        if detail.get("mechanism"):
            item = _text(detail["mechanism"]) + " (recorded status: " + _text(detail.get("status", "unavailable")) + ")"
            if item not in mechanisms:
                mechanisms.append(item)
    if mechanisms:
        sentences.append("Curated mechanism context: " + "; ".join(mechanisms) + ".")
    else:
        sentences.append("An exact curated disease-mechanism description was not supplied in the selected allele evidence.")
    sentences.append("Sequence consequences and database assertions do not establish the allele's functional effect or causality in this child.")
    return " ".join(sentences)


def render_case_rationale(selected_rows, evidence_by_id, disease_profiles, sensitivity=None, limit=3):
    """Render leading selected rows from exact saved V4 evidence, without IO.

    ``selected_rows`` contains hypothesis dictionaries, not (row, choice) tuples.
    ``evidence_by_id`` maps normalized allele ID to its V3/V4 evidence dictionary.
    ``disease_profiles`` is gene -> list of disease phenotype-score dictionaries.
    Optional ``sensitivity`` is scenario name -> compact ranked hypothesis list;
    reported ranks match gene, disease, inheritance, and the complete allele set.
    """
    if not isinstance(limit, int) or isinstance(limit, bool) or limit < 0:
        raise ValueError("Report limit must be a nonnegative integer")
    lines = ["## Evidence behind the leading hypotheses", "",
             "These paragraphs describe the saved evidence for the selected predictions. Entry order and similarity "
             "scores are research prioritization outputs, not calibrated probabilities or a confirmed diagnosis.", ""]
    modes = {"AR": "autosomal recessive", "AD": "autosomal dominant", "XLR": "X-linked recessive",
             "XLD": "X-linked dominant", "XL": "X-linked, subtype unresolved"}
    for number, row in enumerate(list(selected_rows)[:limit], 1):
        members = row.get("members", [])
        mode = row.get("inheritance", "unknown")
        lines.extend([f"### Entry {number}: {_text(row.get('gene', 'unavailable'))} / {_text(row.get('disease_id', 'unavailable'))}", "",
                      "The ranked explanation uses " + _text(modes.get(mode, mode)) + " inheritance and " +
                      f"{len(members)} allele(s): " + " + ".join(_text(m) for m in members) + ". " +
                      _phenotype(row, disease_profiles), ""])
        allele_evidence = row.get("allele_evidence", [])
        for index, identifier in enumerate(members):
            evidence = allele_evidence[index] if index < len(allele_evidence) else {}
            lines.extend([_allele(identifier, evidence_by_id.get(identifier), evidence), ""])
        lines.extend([_uncertainty(row), ""])
        if sensitivity:
            results = []
            for name, scenario_rows in sensitivity.items():
                matched = [r for r in scenario_rows if _key(r) == _key(row)]
                value = "rank " + _text(matched[0].get("rank", "unavailable")) if len(matched) == 1 else "not uniquely listed in the supplied scenario"
                results.append(_text(name) + ": " + value)
            lines.extend(["Sensitivity of this exact gene/disease/inheritance/allele hypothesis: " + "; ".join(results) +
                          ". These are hypothesis ranks, not necessarily distinct submission-entry ranks; stability does not establish accuracy.", ""])
    return "\n".join(lines)
