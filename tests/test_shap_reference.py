"""Reference choice changes the scientific contrast, independently of cell coverage."""
import unittest
import numpy as np
from DEGAS_python.shap_workflow import choose_reference


class ReferenceTest(unittest.TestCase):
    def test_low_quartile_and_all_observations(self):
        result = choose_reference(np.arange(100), [str(i) for i in range(100)])
        self.assertEqual(result['background_indices'], list(range(25)))
        self.assertEqual(result['observation_indices'], list(range(100)))
        self.assertEqual(result['cutoff'], 24.75)

    def test_boundary_ties_and_constant_scores(self):
        result = choose_reference([0, 1, 1, 1, 2, 3, 4, 5], list('abcdefgh'), background_size=0)
        self.assertEqual(result['background_ids'], list('abcd'))
        constant = choose_reference([1]*8, list('abcdefgh'), background_size=0)
        self.assertEqual(constant['background_ids'], list('abcdefgh'))

    def test_sampling_is_reproducible_and_separate(self):
        ids = [str(i) for i in range(1000)]
        a = choose_reference(np.arange(1000), ids)
        b = choose_reference(np.arange(1000), ids, max_cells=30)
        self.assertEqual(a['background_ids'], b['background_ids'])
        self.assertEqual(len(a['background_ids']), 128)
        self.assertEqual(len(b['observation_ids']), 30)
        self.assertTrue(set(a['background_ids']) <= set(a['eligible_ids']))
        self.assertEqual(b, choose_reference(np.arange(1000), ids, max_cells=30))

    def test_population_and_custom(self):
        result = choose_reference([3, 2, 1], list('abc'), mode='population', background_size=0)
        self.assertEqual(result['background_ids'], list('abc'))
        result = choose_reference([3, 2, 1], list('abc'), mode='custom', custom_ids=['c', 'a'])
        self.assertEqual(result['background_ids'], ['a', 'c'])
        for bad in ([], ['missing'], ['a', 'a']):
            with self.assertRaises(ValueError):
                choose_reference([3, 2, 1], list('abc'), mode='custom', custom_ids=bad)
