# Direct SHAP for DEGAS classifiers

`DEGAS_python.direct_shap.DiseaseScore` exposes the feature extractor and
low-resolution classification head as a differentiable scalar class score.
`explain` uses SHAP GradientExplainer (expected gradients); no surrogate model
is trained. Install SHAP separately; core training dependencies are unchanged.

Supply already-preprocessed matrices in the trained gene order. For five-seed
ensembles, explain each seed using identical background rows and average the
attributions and base values. Check original prediction parity, Monte Carlo
sum residuals, background sensitivity and attribution stability before use.
The output explains normalized model inputs, not raw-count interventions.
Correlated genes and the chosen reference affect attribution. Cell composition
is not attributed unless composition is actually an input to the model.

Use a common declared background when comparing spatial and randomized pools
with one model; do not concatenate attributions from different gene spaces.
PCA should be fitted on a declared reference and then project comparison data.
Model attribution is not independent evidence of disease association.

Validation: `PYTHONPATH=DEGAS_python python -m unittest discover -s tests -p test_direct_shap.py`.
