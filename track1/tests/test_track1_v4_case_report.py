"""Case-rationale rendering tests use only fabricated dictionaries."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from track1_v4_case_report import render_case_rationale


def fabricated():
    row = {"gene": "TOY", "disease_id": "OMIM:999999", "inheritance": "AR",
           "members": ["chr1:100:G:A", "chr1:200:T:C"], "phase": "unknown", "rank": 1,
           "phase_evidence": "No shared informative synthetic phase block",
           "components": {"phenotype": .75}, "flags": [],
           "pair_transcript_evidence": {"common_relevant_transcripts": ["ENST_TOY"]},
           "allele_evidence": [{"selected_transcript": "ENST_TOY.1", "clinical_support": .95,
                                "clinvar_active_review_class": "expert_panel", "clinvar_lower_review_dissent": True},
                               {"selected_transcript": "ENST_TOY.1", "clinical_support": 0}]}
    evidence = {m: {"gt": [0, 1], "quality_pass": True,
                   "annotations": [{"transcript": "ENST_TOY.1", "consequences": ["missense"]}],
                   "population": {"known": False}} for m in row["members"]}
    evidence[row["members"][0]]["population"] = {"known": True, "af": .00002, "all_measurement_quality_verified": True}
    profiles = {"TOY": [{"disease_id": "OMIM:999999", "modes": ["AR"], "score": .75,
                         "terms": [{"term": "HP:9000001", "matched_ancestor": "HP:9000001", "matched_information": 2, "scored": True},
                                   {"term": "HP:9000002", "matched_ancestor": "HP:9000003", "matched_information": 1, "scored": True},
                                   {"term": "HP:9000004", "matched_ancestor": None, "matched_information": 0, "scored": False}],
                         "unscored_terms": ["HP:9000004"]}]}
    return row, evidence, profiles


class CaseReportTests(unittest.TestCase):
    def test_grounded_evidence_phase_population_and_no_input_mutation(self):
        row, evidence, profiles = fabricated()
        before = copy.deepcopy((row, evidence, profiles))
        text = render_case_rationale([row], evidence, profiles)
        self.assertIn("Phenotype similarity is 0.75", text)
        self.assertIn("1 match at the queried HPO term and 1 have support only", text)
        self.assertIn("HP:9000001", text)
        self.assertIn("2e-05", text)
        self.assertIn("population frequency unknown or unusable, not established as absent", text)
        self.assertIn("Phase is unresolved", text)
        self.assertIn("confirmed compound heterozygosity is not claimed", text)
        self.assertIn("active ClinVar review class expert panel", text)
        self.assertIn("lower-review ClinVar dissent is retained", text)
        self.assertEqual((row, evidence, profiles), before)

    def test_sensitivity_requires_complete_alleles_disease_and_inheritance(self):
        row, evidence, profiles = fabricated()
        scenarios = {"match": [{**row, "rank": 7}], "different disease": [{**row, "disease_id": "OMIM:999998", "rank": 1}],
                     "only one allele": [{**row, "members": row["members"][:1], "rank": 1}]}
        text = render_case_rationale([row], evidence, profiles, scenarios)
        self.assertIn("match: rank 7", text)
        self.assertIn("different disease: not uniquely listed", text)
        self.assertIn("only one allele: not uniquely listed", text)
        self.assertIn("stability does not establish accuracy", text)

    def test_single_allele_has_no_invented_de_novo_claim(self):
        row, evidence, profiles = fabricated()
        row.update(inheritance="AD", members=row["members"][:1], phase="not_applicable")
        text = render_case_rationale([row], evidence, profiles)
        self.assertIn("single-allele hypothesis", text)
        self.assertIn("de novo status and parental origin are not established", text)
        self.assertIn("An unambiguous term-level disease profile was not supplied", text)

    def test_limit_missing_evidence_and_markup_are_explicit(self):
        row, _, _ = fabricated()
        row["gene"] = "<script>\nTOY*"
        text = render_case_rationale([row] * 5, {}, {}, limit=3)
        self.assertEqual(text.count("### Entry"), 3)
        self.assertNotIn("<script>", text)
        self.assertIn("variant-level details were not supplied", text)
        self.assertNotIn("passed the recorded call-quality filters", text)
        self.assertIn("Phenotype similarity is 0.75", text)

    def test_trans_evidence_does_not_invent_parental_confirmation(self):
        row, evidence, profiles = fabricated()
        row["phase"] = "trans"
        text = render_case_rationale([row], evidence, profiles)
        self.assertIn("saved source-record phase evidence supports trans", text)
        self.assertIn("does not add parental or experimental validation", text)


if __name__ == "__main__":
    unittest.main()
