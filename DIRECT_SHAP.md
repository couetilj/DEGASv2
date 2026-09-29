# Direct SHAP for DEGAS classifiers

`DEGAS_python.direct_shap.DiseaseScore` wraps the trained feature extractor and
low-resolution classification head in a differentiable scalar class score.
`explain` uses SHAP `GradientExplainer` (expected gradients). It explains the
actual DEGAS classifier, without fitting a surrogate. These are approximate,
background-dependent attributions, not exact combinatorial Shapley values.

Install from this checkout with `pip install -e './DEGAS_python[validation,explain]'`.
The [complete tutorial](examples/validation/README.md) trains DEGAS and writes
SHAP outputs with `--shap`. The workflow explains **all supplied cells**, in resumable batches of 256.
DEGAS prediction also covers all cells. The background is a separate comparison
population; it is not a limit on which cells receive explanations.

The default `--shap-background low-risk-quartile` uses the final ensemble's
average-rank cell scores to define eligibility at or below the empirical 25th
percentile. Boundary ties are included (constant scores make every cell eligible).
It then samples up to 128 eligible cells without replacement, reproducibly using
`--seed`. The same cell IDs are used across every retained size and seed, with
each model's own gene list and preprocessing. This reference is frozen once,
after final pooling, in `reports/final/shap_reference.json`; it records the cutoff,
all eligible IDs, actual background IDs and explained IDs. Holdout patient labels
do not define this reference.

This answers: “Which model inputs account for this cell's predicted probability
relative to cells the ensemble scores as low risk?” **Low predicted risk does not
mean healthy controls.** A pooled quartile can differ in cell type, donor or
study composition. A matched, externally defined reference may be preferable;
compare plausible backgrounds before interpreting biology. Each model's base
value is its mean probability over the common background, so base values can
differ between models.

Choose these options when preparing or running a workflow:

```bash
# Explain all cells relative to the entire low-risk quartile:
python -m DEGAS_python.sweep run --output runs/shap_quartile --shap \
  --shap-background-size 0

# Explain all cells relative to a sample of the full cell population:
python -m DEGAS_python.sweep run --output runs/shap_population --shap \
  --shap-background population --shap-background-size 128

# Prespecified controls or matched reference: one input cell ID per line:
python -m DEGAS_python.sweep run --input data --output runs/shap_custom --shap \
  --shap-background custom --shap-background-ids reference_cells.txt \
  --shap-background-size 0
```

`--shap-max-cells N` explicitly opts into a reproducible sample of cells to
explain; the default 0 means all. `--shap-chunk-size` controls memory and restart
granularity. `--shap-samples` controls expected-gradient Monte Carlo samples per
explained cell, independently of background size. Increasing these budgets can
increase runtime and memory; all-cell attribution output scales as cells × genes
× retained models. The defaults 128 and 512 are computation budgets, not validated
biological or accuracy thresholds. Inspect the residual diagnostics.

Local `run --shap` runs all phases. With SLURM, submit `--phase explain` after
final pooling, using the same resource flags as training. One job restores each
retained model and writes immutable chunks; retrying reuses completed chunks.
Pooling requires all tasks and averages seeds chunkwise into `.npy` arrays:

```python
import json
import numpy as np
from pathlib import Path
folder = Path('runs/shap_quartile/reports/explain/size_20')  # retained size
values = np.load(folder / 'values.npy', mmap_mode='r')  # cells × genes
inputs = np.load(folder / 'inputs.npy', mmap_mode='r')
metadata = json.loads((folder / 'metadata.json').read_text())
# metadata: gene order, observation IDs, background IDs and mean base value
# prediction.npy and residual.npy accompany values; diagnostics.csv is per seed/chunk.
```

The selected ensemble defines the reference only. Attributions explain each
model's raw class-1 probability and its seed mean, **not the nonlinear rank
ensemble**. Different feature sizes remain separate. Expected gradients uses
the supplied background distribution; see the official
[SHAP GradientExplainer documentation](https://shap.readthedocs.io/en/latest/generated/shap.GradientExplainer.html).

With an already fitted classifier (call in batches for large datasets):

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
