"""Actual DEGAS classifier/Cox inference with a singleton final batch."""
import tempfile
import unittest
import numpy as np
import torch
from DEGAS_python.models import load_models
from DEGAS_python.datasets import load_datasets
from DEGAS_python.options import BlankClass_opt, BlankCox_opt
from DEGAS_python.direct_shap import DiseaseScore


class NativePredictionTest(unittest.TestCase):
    def test_singleton_classifier_batch_matches_direct_wrapper(self):
        with tempfile.TemporaryDirectory() as tmp:
            opt = dict(BlankClass_opt, save_dir=tmp, input_shape=3, feature_dim=4, batch_size=2)
            model = load_models(opt)
            x = np.array([[.1, .2, .3], [.7, .3, .2], [.4, .8, .1]], dtype='float32')
            _, loader = load_datasets('eval', opt, x, np.array([0, 1, 1]), x)
            native, _ = model.linear_eval(loader)
            wrapper = DiseaseScore(model)
            with torch.no_grad():
                direct = wrapper(torch.as_tensor(x, device=model.device)).cpu().numpy().ravel()
            np.testing.assert_allclose(native.hazard, direct, atol=1e-6)
            self.assertEqual(len(native), 3)

    def test_cox_singleton_batch_retains_scalar_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            opt = dict(BlankCox_opt, save_dir=tmp, input_shape=3, feature_dim=4, batch_size=2)
            model = load_models(opt)
            x = np.ones((3, 3), dtype='float32')
            _, loader = load_datasets('eval', opt, x, np.array([[3, 1], [4, 0], [6, 1]]), x)
            native, _ = model.linear_eval(loader)
            self.assertEqual(len(native), 3)
            self.assertTrue(np.isfinite(native.hazard).all())
