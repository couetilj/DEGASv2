import unittest
import numpy as np
import torch
from DEGAS_python.direct_shap import explain, DiseaseScore

class DirectShapTest(unittest.TestCase):
    def test_linear_attribution(self):
        model=torch.nn.Linear(3,1)
        with torch.no_grad():
            model.weight.copy_(torch.tensor([[2.,-1.,0.]]));model.bias.fill_(.3)
        x=np.array([[1.,2.,3.],[2.,1.,4.]],dtype='float32')
        result=explain(model,np.zeros((4,3),dtype='float32'),x,nsamples=32)
        np.testing.assert_allclose(result['values'],x*np.array([2.,-1.,0.]),atol=1e-6)
        np.testing.assert_allclose(result['residual'],0,atol=1e-6)

    def test_wrapper(self):
        class Extract(torch.nn.Module):
            def forward(self,x): return x,x*2
        class Mock: pass
        model=Mock();model.feature_extractor_layer=Extract();model.low_reso_pred_layer=torch.nn.Linear(3,2)
        wrapper=DiseaseScore(model);x=torch.randn(8,3,requires_grad=True)
        torch.testing.assert_close(wrapper(x)[:,0],torch.softmax(model.low_reso_pred_layer(x*2),1)[:,1])
        wrapper(x).sum().backward();self.assertTrue(torch.isfinite(x.grad).all())

class InvalidShapTest(unittest.TestCase):
    def test_empty_and_nonfinite_inputs(self):
        model = torch.nn.Linear(3, 1)
        for background, observations in [
            (np.empty((0, 3)), np.ones((1, 3))),
            (np.ones((1, 3)), np.empty((0, 3))),
            (np.ones((2, 3)), np.full((1, 3), np.nan)),
            (np.ones((2, 2)), np.ones((1, 3))),
        ]:
            with self.assertRaises(ValueError):
                explain(model, background, observations)

    def test_rejects_cox_and_invalid_class(self):
        class Extract(torch.nn.Module):
            def forward(self, x):
                return x, x
        class Mock:
            pass
        model = Mock()
        model.feature_extractor_layer = Extract()
        model.low_reso_pred_layer = torch.nn.Linear(3, 1)
        with self.assertRaises(ValueError):
            DiseaseScore(model)(torch.ones((2, 3)))
        model.low_reso_pred_layer = torch.nn.Linear(3, 2)
        with self.assertRaises(ValueError):
            DiseaseScore(model, class_index=2)(torch.ones((2, 3)))
        with self.assertRaises(ValueError):
            DiseaseScore(model, class_index=-1)

if __name__=='__main__':unittest.main()
