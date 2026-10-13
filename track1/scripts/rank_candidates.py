"""Untrained, candidate-level research rankers. No network or patient data bundled."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config/ranking.json"


def unit(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number in [0,1]")
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{name} must be a finite number in [0,1]")
    return float(value)


def validate(data, config):
    if data.get("assembly") != "GRCh38":
        raise ValueError("Only independently verified GRCh38 inputs are supported")
    if not isinstance(data.get("synthetic"), bool):
        raise ValueError("synthetic must explicitly be true or false")
    weights = config["evidence_weights"]
    if set(weights) != {"rarity", "effect", "clinical"}:
        raise ValueError("Expected rarity, effect and clinical weights")
    for name, value in weights.items():
        unit(value, name)
    if not math.isclose(sum(weights.values()), 1.0):
        raise ValueError("Evidence weights must sum to one")
    unit(config["phenotype_weight"], "phenotype_weight")
    unit(config["unknown_phase_multiplier"], "unknown_phase_multiplier")
    variants = data.get("variants", [])
    if not variants:
        raise ValueError("No candidate variants supplied")
    ids, alleles = set(), set()
    for v in variants:
        if not isinstance(v.get("id"), str) or not v["id"] or v["id"] in ids:
            raise ValueError("Each candidate needs a unique, nonempty id")
        ids.add(v["id"])
        if not isinstance(v.get("gene"), str) or not v["gene"]:
            raise ValueError(f"{v['id']}: gene is required")
        if v.get("genotype") not in {"het", "hom_alt", "hemizygous"}:
            raise ValueError(f"{v['id']}: unsupported/missing genotype")
        if v.get("inheritance") not in {"AR", "AD", "X", "unknown"}:
            raise ValueError(f"{v['id']}: inheritance must be AR, AD, X or unknown")
        if v.get("qc") not in {"pass", "review", "fail"}:
            raise ValueError(f"{v['id']}: qc must be pass, review or fail")
        if v.get("chrom") not in {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}:
            raise ValueError(f"{v['id']}: use canonical chr-prefixed contigs")
        if isinstance(v.get("pos"), bool) or not isinstance(v.get("pos"), int) or v["pos"] <= 0:
            raise ValueError(f"{v['id']}: pos must be a positive integer")
        for field in ("ref", "alt"):
            if not v.get(field) or not set(v[field]) <= set("ACGT"):
                raise ValueError(f"{v['id']}: only normalized simple A/C/G/T alleles supported")
        if v["ref"] == v["alt"]:
            raise ValueError(f"{v['id']}: REF and ALT must differ")
        allele = tuple(v[k] for k in ("chrom", "pos", "ref", "alt"))
        if allele in alleles:
            raise ValueError("Duplicate allele; collapse transcript rows before ranking")
        alleles.add(allele)
        for field in (*weights, "phenotype"):
            if field not in v:
                raise ValueError(f"{v['id']}: {field} required; use null when unknown")
            if v[field] is not None:
                unit(v[field], f"{v['id']}.{field}")
        if not isinstance(v.get("provenance"), str) or not v["provenance"].strip():
            raise ValueError(f"{v['id']}: provenance is required")
    by_id = {v["id"]: v for v in variants}
    seen_pairs = set()
    for pair in data.get("pairs", []):
        members = pair.get("members", [])
        if len(members) != 2 or len(set(members)) != 2 or not set(members) <= ids:
            raise ValueError("Pair needs two distinct known candidate IDs")
        key = tuple(sorted(members))
        if key in seen_pairs:
            raise ValueError("Duplicate pair")
        seen_pairs.add(key)
        a, b = [by_id[x] for x in members]
        if a["gene"] != b["gene"] or a["chrom"] != b["chrom"]:
            raise ValueError("A recessive pair must be in the same gene and chromosome")
        if any(v["genotype"] != "het" or v["inheritance"] != "AR" for v in (a, b)):
            raise ValueError("Pair model currently supports AR heterozygotes only")
        if pair.get("phase") not in {"trans", "cis", "unknown"}:
            raise ValueError("Pair phase must be trans, cis or unknown")
        if not isinstance(pair.get("phase_evidence"), str) or not pair["phase_evidence"].strip():
            raise ValueError("Provide phase_evidence; state why phase is unknown if unresolved")


def evidence_score(v, config):
    """Missingness has an explicit lower bound; never renormalize missing features."""
    weights = config["evidence_weights"]
    observed = {k: w * v[k] for k, w in weights.items() if v[k] is not None}
    score = sum(observed.values())
    missing_weight = sum(w for k, w in weights.items() if v[k] is None)
    return {"score": score, "missingness_upper_bound": score + missing_weight,
            "contributions": observed, "missing": [k for k in weights if v[k] is None]}


def phenotype_score(v, config):
    e = evidence_score(v, config)
    w = config["phenotype_weight"]
    present = v["phenotype"] is not None
    return {"score": (1-w)*e["score"] + w*(v["phenotype"] if present else 0),
            "missingness_upper_bound": (1-w)*e["missingness_upper_bound"] + w*(v["phenotype"] if present else 1),
            "contributions": {**{k: (1-w)*value for k, value in e["contributions"].items()},
                              **({"phenotype": w*v["phenotype"]} if present else {})},
            "missing": e["missing"] + ([] if present else ["phenotype"])}


def ordered(rows):
    rows = sorted(rows, key=lambda r: (-r["score"], tuple(r["members"])))
    return [{"rank": i+1, **row} for i, row in enumerate(rows)]


def rank(data, config):
    validate(data, config)
    variants = {v["id"]: v for v in data["variants"]}
    evidence, phenotype, hypotheses, excluded = [], [], [], []
    for v in variants.values():
        if v["qc"] != "pass":
            excluded.append({"members": [v["id"]], "reason": f"QC {v['qc']}; review separately"})
            continue
        base = {"members": [v["id"]], "gene": v["gene"]}
        evidence.append({**base, **evidence_score(v, config)})
        p = phenotype_score(v, config)
        phenotype.append({**base, **p})
        complete = (v["inheritance"] == "AR" and v["genotype"] == "hom_alt") or v["inheritance"] == "AD"
        # X-linked context and unknown inheritance need more information than this prototype has.
        if complete:
            hypotheses.append({**base, **p, "kind": "single", "phase": "not_applicable",
                               "caveat": "Inheritance-compatible candidate only; causality unconfirmed"})
        else:
            excluded.append({**base, "reason": "Single-allele hypothesis incomplete or requires inheritance review"})
    for pair in data.get("pairs", []):
        members = sorted(pair["members"])
        a, b = [variants[x] for x in members]
        if pair["phase"] == "cis" or any(v["qc"] != "pass" for v in (a,b)):
            excluded.append({"members": members, "reason": "Cis pair or candidate QC not passed"})
            continue
        pa, pb = phenotype_score(a, config), phenotype_score(b, config)
        factor = 1.0 if pair["phase"] == "trans" else config["unknown_phase_multiplier"]
        hypotheses.append({"members": members, "gene": a["gene"], "kind": "pair",
            "score": factor * min(pa["score"], pb["score"]),
            "missingness_upper_bound": factor * min(pa["missingness_upper_bound"], pb["missingness_upper_bound"]),
            "phase": pair["phase"], "phase_evidence": pair["phase_evidence"],
            "phase_multiplier": factor, "member_scores": [pa, pb],
            "caveat": "Candidate pair; trans NOT established" if pair["phase"] == "unknown" else "Trans evidence supplied by user; not verified by code"})
    return {"synthetic": data["synthetic"], "assembly": data["assembly"],
            "status": "Untrained heuristic research rankings; NOT calibrated probabilities or challenge submissions",
            "evidence_only": ordered(evidence), "phenotype_aware": ordered(phenotype),
            "inheritance_aware": ordered(hypotheses), "excluded_or_review": excluded}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.resolve() in {args.input.resolve(), args.config.resolve()}:
        parser.error("Output cannot overwrite an input or configuration")
    try:
        source, config_bytes = args.input.read_bytes(), args.config.read_bytes()
        result = rank(json.loads(source), json.loads(config_bytes))
    except (ValueError, KeyError, TypeError) as exc:
        parser.error(str(exc))
    result["provenance"] = {
        "input_sha256": hashlib.sha256(source).hexdigest(),
        "config_sha256": hashlib.sha256(config_bytes).hexdigest(),
        "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"Wrote three {'SYNTHETIC' if result['synthetic'] else 'UNVALIDATED research'} rankings to {args.output}")


if __name__ == "__main__":
    main()
