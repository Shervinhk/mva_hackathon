"""Feedback revision contract tests using fabricated CSV rows only.

No patient files, hidden answer keys, or other participants' submissions are read.
"""
import copy
import importlib.util
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import prepare_track1_feedback_revision as revision


FIELDS = (
    "proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2",
    "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes",
)


def row(first, score, second=None, finding_type="primary"):
    result = dict.fromkeys(FIELDS, "")
    result.update(proband_id="PROBAND01", chrom_1=first[0], pos_1=str(first[1]),
                  ref_1=first[2], alt_1=first[3], epcr=str(score),
                  finding_type=finding_type, notes="Fabricated fixture rationale.")
    if second:
        result.update(chrom_2=second[0], pos_2=str(second[1]),
                      ref_2=second[2], alt_2=second[3])
    return result


A = ("chr1", 101, "A", "G")
B = ("chr1", 202, "C", "T")
C = ("chr2", 303, "G", "A")
D = ("chr2", 404, "T", "C")


def rows():
    return [row(A, 0.5), row(C, 1 / 3), row(A, 0.25, B),
            row(D, 0.2, finding_type="secondary")]


def variants(item):
    result = {(item["chrom_1"], int(item["pos_1"]), item["ref_1"], item["alt_1"])}
    if item["chrom_2"]:
        result.add((item["chrom_2"], int(item["pos_2"]), item["ref_2"], item["alt_2"]))
    return frozenset(result)


class FeedbackRevisionTests(unittest.TestCase):
    def test_promotes_complete_pair_without_mutating_input(self):
        source = rows()
        original = copy.deepcopy(source)
        revised, order = revision.revise_rows(source, 3)
        self.assertEqual(source, original)
        self.assertEqual(order, [3, 1, 2, 4])
        self.assertEqual(variants(revised[0]), frozenset([A, B]))
        for output_row, previous_rank in zip(revised, order):
            for field in set(FIELDS) - {"epcr", "notes"}:
                self.assertEqual(output_row[field], original[previous_rank - 1][field])

    def test_rank_refers_to_epcr_sorted_source_not_physical_file_order(self):
        source = rows()
        shuffled = [source[3], source[2], source[0], source[1]]
        revised, order = revision.revise_rows(shuffled, 3)
        self.assertEqual(order, [3, 1, 2, 4])
        self.assertEqual(variants(revised[0]), frozenset([A, B]))
        self.assertEqual(variants(revised[1]), frozenset([A]))

    def test_top_epcr_is_unique_and_all_scores_are_valid(self):
        revised, _ = revision.revise_rows(rows(), 3)
        scores = [float(item["epcr"]) for item in revised]
        self.assertTrue(all(0 < score <= 1 for score in scores))
        self.assertTrue(all(a > b for a, b in zip(scores, scores[1:])))

    def test_feedback_is_explicit_and_original_rationale_retained(self):
        revised, _ = revision.revise_rows(rows(), 3)
        self.assertIn("feedback", revised[0]["notes"].lower())
        self.assertIn("Fabricated fixture rationale.", revised[0]["notes"])

    def test_rank_one_feedback_keeps_allele_order(self):
        source = rows()
        revised, order = revision.revise_rows(source, 1)
        self.assertEqual(order, [1, 2, 3, 4])
        self.assertEqual([variants(x) for x in source], [variants(x) for x in revised])

    def test_rank_outside_source_is_rejected(self):
        for rank in (0, -1, 5):
            with self.subTest(rank=rank), self.assertRaises(ValueError):
                revision.revise_rows(rows(), rank)

    def test_empty_submission_is_rejected(self):
        with self.assertRaises(ValueError):
            revision.revise_rows([], 1)

    def test_ambiguous_tied_source_scores_are_rejected(self):
        source = rows()
        source[1]["epcr"] = source[2]["epcr"]
        with self.assertRaises(ValueError):
            revision.revise_rows(source, 3)

    def test_nonfinite_and_out_of_range_scores_are_rejected(self):
        for bad_score in ("nan", "inf", "0", "-0.1", "1.1"):
            with self.subTest(score=bad_score):
                source = rows()
                source[0]["epcr"] = bad_score
                with self.assertRaises(ValueError):
                    revision.revise_rows(source, 3)

    def test_incomplete_pair_is_rejected(self):
        source = rows()
        source[2]["alt_2"] = ""
        with self.assertRaises(ValueError):
            revision.revise_rows(source, 3)

    def test_more_than_ten_rows_is_rejected(self):
        source = [row(("chr3", 1000 + i, "A", "C"), 1 / (i + 2)) for i in range(11)]
        with self.assertRaises(ValueError):
            revision.revise_rows(source, 3)


class OfficialScoringTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        evaluator_path = ROOT / "research/official/evaluation_feedback_review.py"
        if not evaluator_path.is_file():
            evaluator_path = ROOT / "research/official/evaluation.py"
        if not evaluator_path.is_file():
            raise unittest.SkipTest("Official public scorer snapshot is not bundled in this checkout.")
        spec = importlib.util.spec_from_file_location("feedback_test_public_evaluator", evaluator_path)
        cls.evaluator = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = cls.evaluator
        spec.loader.exec_module(cls.evaluator)

    def score(self, items):
        ordered = sorted(items, key=lambda item: -float(item["epcr"]))
        parsed = [self.evaluator.SubmissionRow(variants(item), float(item["epcr"]), rank,
                                              item["finding_type"])
                  for rank, item in enumerate(ordered, 1)]
        return self.evaluator.score_proband("PROBAND01", parsed, frozenset([A, B]))

    def test_reported_metrics_reproduced_then_perfect_after_revision(self):
        original = self.score(rows())
        self.assertEqual((original.rank_points, original.f_max, original.f_max_threshold,
                          original.full_match_rank), (50.0, 0.8, 0.25, 3))
        revised, _ = revision.revise_rows(rows(), 3)
        result = self.score(revised)
        self.assertEqual((result.rank_points, result.f_max, result.full_match_rank),
                         (100.0, 1.0, 1))

    def test_lower_ranked_hypotheses_do_not_reduce_perfect_fmax(self):
        revised, _ = revision.revise_rows(rows(), 3)
        self.assertEqual(self.score(revised).f_max, self.score(revised[:1]).f_max)

    def test_top_tied_false_positive_hurts_fmax_even_if_secondary(self):
        items = [row(A, 0.5, B), row(C, 0.5, finding_type="secondary")]
        result = self.score(items)
        self.assertEqual((result.rank_points, result.f_max), (100.0, 0.8))


if __name__ == "__main__":
    unittest.main()
