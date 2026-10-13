import copy
import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rank_candidates", ROOT / "scripts/rank_candidates.py")
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)


class RankingTests(unittest.TestCase):
    def setUp(self):
        self.data = json.loads((ROOT / "examples/synthetic_candidates.json").read_text())
        self.config = json.loads((ROOT / "config/ranking.json").read_text())

    def test_phenotype_changes_leading_candidate(self):
        result = model.rank(self.data, self.config)
        self.assertEqual(result["evidence_only"][0]["members"], ["toy_C"])
        self.assertEqual(result["phenotype_aware"][0]["members"], ["toy_A"])

    def test_missing_evidence_is_explicit_and_not_renormalized(self):
        v = self.data["variants"][1]
        score = model.evidence_score(v, self.config)
        self.assertEqual(score["missing"], ["clinical"])
        self.assertAlmostEqual(score["missingness_upper_bound"] - score["score"], 0.35)
        v["clinical"] = 1
        self.assertGreater(model.evidence_score(v, self.config)["score"], score["score"])

    def test_cis_pair_excluded(self):
        self.data["pairs"][0]["phase"] = "cis"
        rows = model.rank(self.data, self.config)["inheritance_aware"]
        self.assertFalse(any(len(r["members"]) == 2 for r in rows))

    def test_unknown_phase_is_flagged_and_penalized(self):
        before = next(r for r in model.rank(self.data, self.config)["inheritance_aware"] if r["kind"] == "pair")
        self.assertIn("NOT established", before["caveat"])
        self.data["pairs"][0].update(phase="trans", phase_evidence="Synthetic parental phase")
        after = next(r for r in model.rank(self.data, self.config)["inheritance_aware"] if r["kind"] == "pair")
        self.assertGreater(after["score"], before["score"])

    def test_qc_review_does_not_enter_rankings(self):
        for key in ["evidence_only", "phenotype_aware", "inheritance_aware"]:
            self.assertFalse(any("toy_E" in r["members"] for r in model.rank(self.data, self.config)[key]))

    def test_single_ar_heterozygote_not_complete_hypothesis(self):
        rows = model.rank(self.data, self.config)["inheritance_aware"]
        self.assertFalse(any(r["members"] == ["toy_A"] for r in rows))

    def test_rejects_nan_duplicate_wrong_build_cross_gene_pairs(self):
        invalid = []
        d = copy.deepcopy(self.data); d["variants"][0]["effect"] = float("nan"); invalid.append(d)
        d = copy.deepcopy(self.data); d["variants"].append(copy.deepcopy(d["variants"][0])); invalid.append(d)
        d = copy.deepcopy(self.data); d["assembly"] = "GRCh37"; invalid.append(d)
        d = copy.deepcopy(self.data); d["pairs"][0]["members"] = ["toy_A","toy_C"]; invalid.append(d)
        for d in invalid:
            with self.subTest(data=d), self.assertRaises(ValueError):
                model.rank(d, self.config)

    def test_deterministic_order_and_bounds(self):
        result = model.rank(self.data, self.config)
        self.data["variants"].reverse()
        for key in ["evidence_only", "phenotype_aware", "inheritance_aware"]:
            self.assertEqual(result[key], model.rank(self.data, self.config)[key])
            for row in result[key]:
                self.assertTrue(0 <= row["score"] <= row["missingness_upper_bound"] <= 1)


if __name__ == "__main__":
    unittest.main()
