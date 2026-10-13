"""Fabricated inputs only; no private callsets, outputs, or clinical findings."""
import copy
import importlib.util
from pathlib import Path
import unittest

from track1_v4 import ranking as m


def config():
    return {"assembly": "GRCh38", "phenotype_weight": .4, "mechanism_weight": 0,
            "predictor_weight": 0, "unknown_phase_factor": 1,
            "maximum_population_af_ar": .01, "maximum_population_af_ad": .001,
            "minimum_allele_utility": .3, "mechanism_rules": []}


def annotation(transcript="ENST_TOY.1", consequence="missense", score=.35):
    return {"gene": "TOY", "transcript": transcript, "consequences": [consequence],
            "effect_score": score, "amino_acid_change": "2E>2K"}


def variant(pos=4000000, chrom="chr1", genotype="het", anns=None, assertions=None):
    gt = {"het": [0, 1], "hom_alt": [1, 1], "hemizygous": [1]}[genotype]
    v = {"id": f"{chrom}:{pos}:G:A", "chrom": chrom, "pos": pos, "ref": "G", "alt": "A",
         "genotype": genotype, "gt": gt, "alt_index": 1, "quality_pass": True,
         "record_id": f"{chrom}:{pos}:G:A", "phase_set": None, "haplotype": None,
         "annotations": anns or [annotation()], "clinvar_records": [],
         "public_annotation": {"rcv_assertions": assertions or []}}
    v["population"] = m._v3.population_evidence(v)
    return v


def disease(mode="AD", score=.8):
    return {"disease_id": "OMIM:999999", "modes": [mode], "score": score}


def assertion(label="Pathogenic", stars=1, **extra):
    return {"accession": "SYNTHETIC_ASSERTION", "assertion_type": "SCV",
            "disease_ids": ["OMIM:999999"], "classification": label,
            "review_stars": stars, "review_status": "criteria provided, single submitter", **extra}


def run(variants, mode="AD", cfg=None, **options):
    return m.rank(variants, {"TOY": [disease(mode)]}, cfg or config(), **options)


class ClinVarHierarchyTests(unittest.TestCase):
    def test_lower_uncertainty_does_not_erase_expert_pathogenic_support(self):
        v = variant(assertions=[assertion(stars=3), assertion("Uncertain significance")])
        e = m.rcv_evidence(v, disease())
        self.assertEqual(e["support"], .95)
        self.assertFalse(e["conflict"])
        self.assertTrue(e["lower_review_dissent"])
        self.assertEqual(len(e["assertions"]), 2)
        self.assertEqual([a["active_review_class_assertion"] for a in e["assertions"]], [True, False])
        ranked, _, _ = run([v])
        self.assertIn("lower_review_ClinVar_dissent_retained", ranked[0]["flags"])

    def test_conflicting_equal_experts_remain_conflicting(self):
        e = m.rcv_evidence(variant(assertions=[assertion(stars=3), assertion("Benign", stars=3)]), disease())
        self.assertTrue(e["conflict"])
        self.assertEqual(e["support"], 0)
        self.assertFalse(e["benign"])

    def test_guideline_precedence_and_ordinary_conflict_are_distinct(self):
        e = m.rcv_evidence(variant(assertions=[assertion("Benign", stars=4), assertion(stars=3)]), disease())
        self.assertTrue(e["benign"])
        self.assertTrue(e["lower_review_dissent"])
        e = m.rcv_evidence(variant(assertions=[assertion(stars=2), assertion("Uncertain significance", stars=1)]), disease())
        self.assertTrue(e["conflict"])

    def test_ignored_context_and_other_condition_never_gain_support(self):
        items = [assertion(stars=4, contributes_to_aggregate=False),
                 assertion(stars=3, classification_type="somatic"),
                 assertion(stars=3, disease_ids=["OMIM:999998"])]
        e = m.rcv_evidence(variant(assertions=items), disease())
        self.assertEqual(e["support"], 0)
        self.assertEqual(e["assertions"], [])

    def test_import_does_not_monkeypatch_v3_module(self):
        spec = importlib.util.spec_from_file_location("v3_isolation_regression", Path(m.__file__).parents[1] / "scripts/track1_v3.py")
        original = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(original)
        v = variant(assertions=[assertion(stars=3), assertion("Uncertain significance")])
        self.assertEqual(original.rcv_evidence(v, disease())["support"], 0)
        self.assertEqual(m.rcv_evidence(v, disease())["support"], .95)


class EligibilityTests(unittest.TestCase):
    def test_xlr_heterozygote_is_deferred_not_declared_impossible(self):
        ranked, review, _ = run([variant(chrom="chrX")], "XLR")
        self.assertFalse(ranked)
        self.assertEqual(review[0]["inheritance"], "XLR")
        self.assertEqual(review[0]["inheritance_compatibility"], "compatibility_unresolved")
        self.assertIn("diploid_XLR_heterozygote_may_manifest_context_required", review[0]["flags"])
        self.assertFalse(review[0]["clinical_confirmation"])

    def test_documented_manifesting_heterozygote_can_remain_candidate(self):
        cfg = config()
        cfg["inheritance_context"] = [{"gene": "TOY", "disease_id": "OMIM:999999", "inheritance": "XLR",
                                        "source": "https://example.org/fabricated", "manifesting_heterozygote_compatible": True}]
        ranked, review, _ = run([variant(chrom="chrX")], "XLR", cfg)
        self.assertTrue(ranked[0]["automatic_shortlist_eligible"])
        self.assertFalse(review)
        self.assertFalse(ranked[0]["clinical_confirmation"])

    def test_xlr_haploid_and_xld_diploid_do_not_infer_sex(self):
        for mode, genotype in [("XLR", "hemizygous"), ("XLD", "het")]:
            ranked, review, _ = run([variant(chrom="chrX", genotype=genotype)], mode)
            self.assertFalse(review)
            self.assertEqual(ranked[0]["inheritance"], mode)
            self.assertIn("observed_genotype_does_not_establish_sex_or_parental_origin", ranked[0]["flags"])

    def test_generic_xl_and_par_need_context(self):
        for mode, position in [("XL", 4000000), ("XLD", 10001), ("XLR", 155701383)]:
            ranked, review, _ = run([variant(pos=position, chrom="chrX", genotype="hemizygous")], mode)
            self.assertFalse(ranked)
            self.assertTrue(review)

    def test_curated_missense_mechanism_defers_unsupported_truncation(self):
        cfg = config()
        cfg["mechanism_rules"] = [{"gene": "TOY", "disease_id": "OMIM:999999", "inheritance": "AD",
                                  "mechanism": "specific altered protein function", "source": "https://example.org/fabricated",
                                  "supported_consequences": [], "plausible_consequences": ["missense"]}]
        v = variant(anns=[annotation(consequence="stop_gained", score=.8)])
        ranked, review, _ = run([v], cfg=cfg)
        self.assertFalse(ranked)
        self.assertEqual(review[0]["mechanism_compatibility"], "compatibility_unresolved")
        self.assertEqual(review[0]["score"], .6 * .8 + .4 * .8)
        # Condition support can warrant review despite the generic mechanism warning.
        v["public_annotation"]["rcv_assertions"] = [assertion(stars=3)]
        ranked, review, _ = run([v], cfg=cfg)
        self.assertTrue(ranked)
        self.assertFalse(review)

    def test_unsupported_selected_effect_is_not_rescued_by_other_transcript(self):
        cfg = config()
        cfg["mechanism_rules"] = [{"gene": "TOY", "disease_id": "OMIM:999999", "inheritance": "AD",
                                  "mechanism": "specific altered protein function", "source": "https://example.org/fabricated",
                                  "plausible_consequences": ["missense"]}]
        v = variant(anns=[annotation("ENST_TRUNC.1", "stop_gained", .8), annotation("ENST_MISSENSE.1")])
        ranked, review, _ = run([v], cfg=cfg)
        self.assertFalse(ranked)
        self.assertTrue(review)

    def test_pair_disjoint_functional_transcripts_are_deferred(self):
        a = variant(anns=[annotation("ENST_A.1"), annotation("ENST_B.1", "synonymous", .05)])
        b = variant(pos=4000001, anns=[annotation("ENST_B.1"), annotation("ENST_A.1", "synonymous", .05)])
        ranked, review, _ = run([a, b], "AR")
        self.assertFalse(ranked)
        self.assertEqual(review[0]["pair_transcript_evidence"]["status"], "disjoint_available_relevant_transcripts")
        self.assertEqual(review[0]["phase"], "unknown")

    def test_pair_annotation_intersection_not_only_selected_transcripts(self):
        a = variant(anns=[annotation("ENST_A.1", "stop_gained", .8), annotation("ENST_COMMON.1")])
        b = variant(pos=4000001, anns=[annotation("ENST_B.1", "stop_gained", .8), annotation("ENST_COMMON.2")])
        ranked, review, _ = run([a, b], "AR")
        self.assertFalse(review)
        row = ranked[0]
        self.assertEqual(row["pair_transcript_evidence"]["common_relevant_transcripts"], ["ENST_COMMON"])
        self.assertEqual(row["phase"], "unknown")
        self.assertFalse(row["clinical_confirmation"])

    def test_splice_rescue_relevant_transcript_is_retained_with_warning(self):
        a = variant(anns=[annotation(consequence="synonymous", score=.05)])
        a["public_annotation"]["other_predictors"] = {"SPLICE_AI": {"score": .8}}
        b = variant(pos=4000001)
        ranked, review, _ = run([a, b], "AR")
        self.assertFalse(review)
        self.assertEqual(ranked[0]["pair_transcript_evidence"]["common_relevant_transcripts"], ["ENST_TOY"])
        self.assertIn("splice_prediction_rescue_requires_transcript_and_RNA_review", ranked[0]["flags"])

    def test_no_rule_keeps_default_utility_and_inputs_unchanged(self):
        v = variant()
        before = copy.deepcopy(v)
        ranked, review, _ = run([v])
        self.assertFalse(review)
        self.assertAlmostEqual(ranked[0]["score"], .6 * .35 + .4 * .8)
        self.assertEqual(v, before)

    def test_invalid_phenotype_score_rejected(self):
        for invalid in [2, -1, float("nan"), float("inf")]:
            with self.assertRaises(ValueError):
                m.rank([variant()], {"TOY": [disease(score=invalid)]}, config())

    def test_distinct_rows_preserve_complete_pair_and_order(self):
        rows = [{"members": ["a", "b"]}, {"members": ["b", "a"]}, {"members": ["a"]}]
        self.assertEqual(m.distinct_rows(rows), [rows[0], rows[2]])
        self.assertEqual(m.distinct_rows(rows, limit=0), [])
        self.assertEqual(m.distinct_rows(rows, limit=1), [rows[0]])


if __name__ == "__main__":
    unittest.main()
