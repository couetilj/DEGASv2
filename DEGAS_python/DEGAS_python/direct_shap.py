"""Direct expected-gradient explanations of DEGAS classification scores.

Optional dependencies: torch, numpy, shap. Inputs must already have exactly
training-compatible preprocessing and feature order. Does not explain the
raw-count preprocessing operation. No surrogate is fitted.
"""
import numpy as np
import torch
from torch import nn


class DiseaseScore(nn.Module):
    """Differentiable equivalent of softmax(low_reso_pred_layer(embedding))."""
    def __init__(self, degas, class_index=1):
        super().__init__()
        self.extractor = degas.feature_extractor_layer
        self.head = degas.low_reso_pred_layer
        if not isinstance(class_index, int) or class_index < 0:
            raise ValueError('class_index must be a nonnegative integer')
        self.class_index = class_index
        self.eval()

    def forward(self, x):
        _, embedding = self.extractor(x)
        logits = self.head(embedding)
        if logits.ndim != 2 or logits.shape[1] < 2 or self.class_index >= logits.shape[1]:
            raise ValueError('DiseaseScore requires a classification head and valid class_index')
        return torch.softmax(logits, dim=1)[:, self.class_index:self.class_index+1]


def explain(model, background, observations, nsamples=512, seed=42, batch_size=128):
    """Return approximate attributions, prediction, base value and residual.

    Residual = prediction - base - sum(attributions); callers must inspect
    Monte Carlo accuracy before interpretation. Model is put in eval mode.
    Background is a declared scientific reference, not inferred here.
    """
    if not isinstance(nsamples, int) or nsamples < 1 or not isinstance(batch_size, int) or batch_size < 1:
        raise ValueError('nsamples and batch_size must be positive integers')
    import shap
    model.eval()
    device = next(model.parameters()).device
    background = torch.as_tensor(background, dtype=torch.float32, device=device)
    observations = torch.as_tensor(observations, dtype=torch.float32, device=device)
    if background.ndim != 2 or observations.ndim != 2 or background.shape[1] != observations.shape[1]:
        raise ValueError('Expected aligned samples-by-features matrices')
    if not len(background) or not len(observations) or not background.shape[1]:
        raise ValueError('Background and observations must be nonempty')
    if not torch.isfinite(background).all() or not torch.isfinite(observations).all():
        raise ValueError('Nonfinite input')
    with torch.no_grad():
        base = sum(model(background[i:i+batch_size]).sum().item()
                   for i in range(0, len(background), batch_size)) / len(background)
        prediction = model(observations).flatten().cpu().numpy()
    explainer = shap.GradientExplainer(model, background, batch_size=batch_size)
    with torch.enable_grad():
        values = explainer.shap_values(observations, nsamples=nsamples, rseed=seed)
    if isinstance(values, list):
        if len(values) != 1: raise ValueError('Expected one model output')
        values = values[0]
    values = np.asarray(values)
    if values.ndim == 3 and values.shape[-1] == 1: values = values[..., 0]
    if values.shape != tuple(observations.shape) or not np.isfinite(values).all():
        raise ValueError('Invalid attribution shape/values')
    return dict(values=values, prediction=prediction, base_value=base,
                residual=prediction-base-values.sum(axis=1), nsamples=nsamples)
