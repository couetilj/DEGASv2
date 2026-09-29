import unittest
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from DEGAS_python.validation import (holdout_split, validation_folds, select_one_se,
                                     average_rank_scores, classification_metrics)


def metadata():
    return pd.DataFrame(dict(patient_id=[f'p{i}' for i in range(80)],
                             study=np.repeat(['A', 'B', 'C', 'D'], 20),
                             label=np.tile([0, 1], 40)))


def predictions():
    rng = np.random.default_rng(7)
    m = metadata()
    tables = []
    for size in (10, 20, 40):
        t = m.copy()
        t['size'] = size
        t['score'] = np.clip(rng.random(len(m)) * .8 + m.label * .2, 0, 1)
        tables.append(t)
    return pd.concat(tables, ignore_index=True)


class ValidationTest(unittest.TestCase):
    def test_patient_split_keeps_repeated_rows_together(self):
        m = pd.concat([metadata(), metadata()], ignore_index=True)
        tr, te = holdout_split(m, mode='patient')
        self.assertEqual(m.iloc[te].patient_id.nunique(), 8)
        self.assertFalse(set(m.iloc[tr].patient_id) & set(m.iloc[te].patient_id))
        dev = m.iloc[tr].reset_index(drop=True)
        folds = validation_folds(dev, mode='patient', n_splits=4)
        seen = []
        for a, b in folds:
            self.assertFalse(set(dev.iloc[a].patient_id) & set(dev.iloc[b].patient_id))
            self.assertEqual(set(dev.iloc[b].label), {0, 1})
            seen.extend(b)
        self.assertEqual(sorted(seen), list(range(len(dev))))
        np.testing.assert_array_equal(te, holdout_split(m, mode='patient')[1])

    def test_study_split_and_inner_disjointness(self):
        m = metadata()
        tr, te = holdout_split(m)
        self.assertEqual(m.iloc[te].study.nunique(), 1)
        self.assertFalse(set(m.iloc[tr].study) & set(m.iloc[te].study))
        dev = m.iloc[tr].reset_index(drop=True)
        folds = validation_folds(dev)
        self.assertEqual(len(folds), 3)
        for a, b in folds:
            self.assertFalse(set(dev.iloc[a].study) & set(dev.iloc[b].study))
        # Changing untouched holdout outcomes cannot change inner fold membership.
        m.loc[te, 'label'] = 1 - m.loc[te, 'label']
        for old, new in zip(folds, validation_folds(m.iloc[tr].reset_index(drop=True))):
            np.testing.assert_array_equal(old[1], new[1])

    def test_invalid_designs(self):
        m = metadata()
        for mode in ('oops',):
            with self.assertRaises(ValueError):
                holdout_split(m, mode=mode)
        with self.assertRaises(ValueError):
            holdout_split(m.assign(study='one'))
        with self.assertRaises(ValueError):
            validation_folds(m.assign(label=0), mode='patient')
        with self.assertRaises(ValueError):
            validation_folds(pd.concat([m, m.iloc[:1].assign(study='another')]))
        with self.assertRaises(ValueError):
            holdout_split(m.assign(patient_id=None))
        with self.assertRaises(ValueError):
            validation_folds(m.assign(label=(m.study == 'A').astype(int)))

    def test_metrics_known_values_and_invalid_probabilities(self):
        r = classification_metrics([0, 1, 0, 1], [.1, .9, .2, .8])
        self.assertEqual(r['AUROC'], 1)
        self.assertEqual(r['average_precision'], 1)
        self.assertAlmostEqual(r['brier'], .025)
        self.assertEqual(r['balanced_accuracy'], 1)
        with self.assertRaises(ValueError):
            classification_metrics([0, 1], [-.1, .8])

    def test_selection_matches_equal_study_means_and_bootstrap(self):
        p = predictions()
        metrics, selection = select_one_se(p, n_bootstrap=400, seed=19)
        pd.testing.assert_frame_equal(selection, select_one_se(p, n_bootstrap=400, seed=19)[1])
        np.testing.assert_allclose(selection.mean_AUROC, metrics.groupby('size').AUROC.mean())
        best = selection.iloc[selection.mean_AUROC.argmax()]
        np.testing.assert_allclose(selection.cutoff, best.mean_AUROC - best.SE)
        np.testing.assert_array_equal(selection.selected, selection.mean_AUROC >= best.mean_AUROC - best.SE)
        # Independent ordinary resampling estimate for the first size.
        rng = np.random.default_rng(19)
        estimates = []
        for _ in range(1000):
            aucs = []
            for _, part in p[p['size'] == 10].groupby('study'):
                sample = pd.concat([g.iloc[rng.integers(len(g), size=len(g))] for _, g in part.groupby('label')])
                aucs.append(roc_auc_score(sample.label, sample.score))
            estimates.append(np.mean(aucs))
        self.assertAlmostEqual(selection.SE.iloc[0], np.std(estimates, ddof=1), delta=.008)

    def test_equal_study_weighting_with_unequal_cohorts(self):
        p = predictions()
        p = p[~p.patient_id.isin([f'p{i}' for i in range(4, 20)])].copy()
        p.loc[p.study == 'A', 'score'] = p.loc[p.study == 'A', 'label']
        metrics, selection = select_one_se(p, n_bootstrap=20)
        np.testing.assert_allclose(selection.mean_AUROC, metrics.groupby('size').AUROC.mean())
        first = metrics[metrics['size'] == 10]
        self.assertNotAlmostEqual(selection.mean_AUROC.iloc[0],
                                  np.average(first.AUROC, weights=first.n_patients))

    def test_selection_rejects_missing_duplicate_and_tiny_strata(self):
        p = predictions()
        for bad in (p.iloc[1:], pd.concat([p, p.iloc[:1]]), p.assign(score=np.nan), p.assign(score=2)):
            with self.assertRaises(ValueError):
                select_one_se(bad, n_bootstrap=5)
        with self.assertRaises(ValueError):
            select_one_se(p[p.patient_id.isin(['p0', 'p1'])], n_bootstrap=5)

    def test_perfect_tied_sizes_all_retained(self):
        p = predictions()
        p['score'] = p.label
        _, selection = select_one_se(p, n_bootstrap=10)
        self.assertTrue(selection.selected.all())
        np.testing.assert_allclose(selection.SE, 0, atol=1e-15)
        np.testing.assert_allclose(selection.cutoff, 1)

    def test_average_ranks_ties_and_order(self):
        x = pd.DataFrame({10: [0, 0, 1], 20: [1, .2, 0]}, index=['a', 'b', 'c'])
        result = average_rank_scores(x, [10, 20])
        np.testing.assert_allclose(result, [.75, 7/12, 2/3])
        self.assertEqual(list(result.index), ['a', 'b', 'c'])
        self.assertFalse(np.allclose(result, x.mean(axis=1).rank()/3))
        with self.assertRaises(ValueError):
            average_rank_scores(x, [])
        with self.assertRaises(ValueError):
            average_rank_scores(x.rename(index={'b': 'a'}), [10])

    def test_plots_and_exports(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from DEGAS_python.evaluation_plots import plot_validation
        metrics, selection = select_one_se(predictions(), n_bootstrap=20)
        with tempfile.TemporaryDirectory() as tmp:
            figures = plot_validation(metrics, selection, output_dir=tmp)
            for file in ('auroc_selection.png', 'evaluation_metrics.pdf', 'feature_size_selection.csv'):
                self.assertGreater((Path(tmp)/file).stat().st_size, 100)
            for f in figures:
                plt.close(f)
        with self.assertRaises(ValueError):
            plot_validation(metrics, selection.assign(cutoff=0))


if __name__ == '__main__':
    unittest.main()
