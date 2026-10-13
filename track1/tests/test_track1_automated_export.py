"""Fabricated hypotheses only; no patient files are opened."""
import copy
import importlib.util
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from track1_automated_export import automated_choices


def row(rank, members=None, **updates):
    result = {'rank': rank, 'score': 1 - rank / 100, 'inheritance': 'AR',
              'phase': 'unknown', 'members': members or [f'chr1:{rank}:A:C'],
              'eligible_for_manual_selection': True}
    result.update(updates)
    return result


class AutomatedExportTests(unittest.TestCase):
    def test_dedup_complete_pairs_and_preserve_order_without_review_attestation(self):
        pair = ['chr1:10:A:C', 'chr1:20:G:T']
        rows = [row(1, pair), row(2, list(reversed(pair))),
                row(3, phase='second_allele_unresolved'), row(4, inheritance='unknown'),
                row(5, eligible_for_manual_selection=False), row(6)]
        before = copy.deepcopy(rows)
        selected, manifest = automated_choices(rows, 'fabricated-hash', lambda r: str(r['rank']))
        self.assertEqual([r['rank'] for r, _ in selected], [1, 6])
        self.assertEqual(selected[0][0]['members'], pair)
        self.assertEqual(selected[0][0]['phase'], 'unknown')
        self.assertEqual([c['epcr'] for _, c in selected], [.5, 1/3])
        self.assertFalse(manifest['reviewed_locally'])
        self.assertTrue(all(not c['evidence_reviewed'] for c in manifest['choices']))
        self.assertIn('not estimated or calibrated probabilities', manifest['confidence_policy'])
        self.assertEqual(rows, before)

    def test_top_ten_and_strict_relative_order(self):
        selected, _ = automated_choices([row(n) for n in range(1, 21)], 'toy', lambda r: str(r['rank']))
        self.assertEqual(len(selected), 10)
        weights = [c['epcr'] for _, c in selected]
        self.assertTrue(all(0 < value <= 1 for value in weights))
        self.assertTrue(all(a > b for a, b in zip(weights, weights[1:])))

    def test_reject_malformed_or_changed_ranking(self):
        cases = [[row(1, score=float('nan'))], [row(1), row(2, score=1.5)],
                 [row(2), row(1)], [row(1, ['chr1:1:A:C'] * 2)],
                 [row(1, ['chr1:1:A:C', 'chr1:2:A:C', 'chr1:3:A:C'])]]
        for rows in cases:
            with self.subTest(case=cases.index(rows)), self.assertRaises(ValueError):
                automated_choices(rows, 'toy', lambda r: str(r['rank']))

    def test_empty_eligible_set_fails_closed(self):
        with self.assertRaises(ValueError):
            automated_choices([row(1, eligible_for_manual_selection=False)], 'toy', lambda r: str(r['rank']))

    def test_official_metrics_depend_on_order_not_weight_magnitudes(self):
        path = Path(__file__).resolve().parents[1] / 'research/official/evaluation.py'
        if not path.exists():
            self.skipTest('Public official scorer snapshot is not included in the source archive.')
        spec = importlib.util.spec_from_file_location('public_track1_scorer_for_test', path)
        scorer = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = scorer
        spec.loader.exec_module(scorer)
        variants = [frozenset({('chr1', n, 'A', 'C')}) for n in range(1, 11)]
        ordinal = [scorer.SubmissionRow(v, 1/(rank+1), rank) for rank, v in enumerate(variants, 1)]
        alternative = [scorer.SubmissionRow(v, (11-rank)/11, rank) for rank, v in enumerate(variants, 1)]
        for truth in variants:
            first = scorer.score_proband('SYNTHETIC', ordinal, truth)
            second = scorer.score_proband('SYNTHETIC', alternative, truth)
            self.assertEqual(first.rank_points, second.rank_points)
            self.assertEqual(first.f_max, second.f_max)
            self.assertEqual(first.full_match_rank, second.full_match_rank)


if __name__ == '__main__':
    unittest.main()
