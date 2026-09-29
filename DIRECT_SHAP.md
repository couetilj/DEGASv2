# Direct SHAP for DEGAS classifiers

`DEGAS_python.direct_shap.DiseaseScore` wraps the trained feature extractor and
low-resolution classification head in a differentiable scalar class score.
`explain` uses SHAP `GradientExplainer` (expected gradients). It explains the
actual DEGAS classifier, without fitting a surrogate. These are approximate,
background-dependent attributions, not exact combinatorial Shapley values.

Install from this checkout with `pip install -e './DEGAS_python[validation,explain]'`.
The [complete tutorial](examples/validation/README.md) trains DEGAS and writes
SHAP outputs with `--shap`. With an already fitted classifier:

```python
import numpy as np
from DEGAS_python.direct_shap import DiseaseScore, explain

# model is the fitted DEGAS model from run_model(..., return_model=True).
# Both matrices are samples x selected genes, already preprocessed with the
# training transform and ordered exactly as the saved training gene list.
score = DiseaseScore(model, class_index=1)
result = explain(score, background, observations, nsamples=512, seed=42)
np.savez_compressed('direct_shap.npz', **result)
print('Mean absolute residual:', np.abs(result['residual']).mean())
print('95th-percentile absolute residual:', np.quantile(np.abs(result['residual']), .95))
```

`values` has shape observations × genes; `base_value` is the mean background
prediction; `prediction` is the actual model output; `residual` is
`prediction - base_value - values.sum(axis=1)`. Verify prediction parity with
saved DEGAS probabilities. Inspect residuals, increase `nsamples` if needed,
repeat the estimator with a second seed, and check an alternate scientifically
reasonable background. Choose numerical tolerances appropriate to the output
scale before interpreting explanations; Monte Carlo error is not biological
uncertainty. The helper rejects empty/nonfinite/misaligned matrices and Cox
heads. Class 1 must actually be the outcome you intend to explain.

For seed ensembles, use identical background and observation rows, then average
attributions, base values and raw predictions across seeds. Linearity applies
to the mean raw probability. Keep per-seed results to assess variability. Use a
common declared background to compare groups with the same model; background
composition changes the question the explanation answers.

```python
import shap
import matplotlib.pyplot as plt

explanation = shap.Explanation(
    values=result['values'],
    base_values=np.repeat(result['base_value'], len(observations)),
    data=observations,
    feature_names=gene_names,
)
shap.plots.beeswarm(explanation, show=False)
plt.savefig('shap_beeswarm.png', dpi=200, bbox_inches='tight')
plt.close()
```

These explanations concern **normalized model inputs and raw probabilities**.
They do not explain the raw-count preprocessing operation or the cohort-relative
average-rank ensemble. Rank transformation is nonlinear, so averaging raw-score
SHAP values does not explain percentile scores. Keep explanations for different
gene lists separate; do not concatenate them or silently zero-fill missing
features. Correlated genes and the reference distribution affect attribution.
Cell composition is not itself attributed unless it is an input feature.
Attribution is not independent evidence of a disease-associated gene.

Validation: `PYTHONPATH=DEGAS_python python -m unittest discover -s tests -p test_direct_shap.py`.
