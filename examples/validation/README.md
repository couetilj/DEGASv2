# Patient validation, one-SE ensembles and direct SHAP

This binary-classification example takes raw bulk and single-cell counts all the
way through DEGAS training, held-out predictions, gene-set-size selection,
publication-style evaluation plots, and optional direct SHAP. It uses the
existing `BlankClass` backend; it does not change the default training API.
Class 1 is the disease/outcome of interest. Cox/survival models need different
metrics and are outside this example.

## Run the synthetic example

From the repository root, in a Python environment with PyTorch support:

```bash
python -m pip install -e './DEGAS_python[validation,explain]'
python examples/validation/run.py --output runs/study --split study --shap
python examples/validation/run.py --output runs/patient --split patient --shap
```

Each command requires a new output directory. Omit `--shap` if not needed;
`[validation]` alone supplies the training and plotting dependencies. For a
quick execution check, append `--iters 2 --seeds 1 --bootstrap 100` (these are
smoke settings, not a trained analysis). The default uses five seeds, 300
iterations and 20,000 bootstrap replicates. Inspect saved training losses and
choose training settings on development data for your own application.

## Your inputs

Use `--input path/to/csvs --sizes 10 20 40 60 120`. The directory contains:

| File | Contents |
|---|---|
| `bulk_counts.csv` | Rows = patients, columns = unique gene IDs; first column is patient ID. |
| `patients.csv` | `patient_id,study,label`; one row per patient, binary label 0/1. |
| `cell_counts.csv` | Rows = reference cells, columns = unique gene IDs; first column is cell ID. |
| `cells.csv` | `cell_id,patient_id`; known donor ID for every reference cell. |

Use globally unique patient IDs across studies and assays. Matrix and metadata
rows are joined by ID. The tutorial requires an independent cell-reference
cohort and rejects overlapping bulk/cell donors. If cohorts overlap, implement
fold-specific reference-cell exclusions too: held-out patients must not enter
through the cell arm. Unknown donor identity cannot establish independence.

The example intersects measured bulk/cell gene IDs and saves each selected gene
list in order. It ranks training-bulk genes by variance of `log1p(counts)` as a
simple illustrative selector, separately inside every fold. This is not a
replacement for the R native marker selector or a recommended selector for
every assay. Replace the selector in `fit()` to use your study design; fit
supervised selection, batch correction and other learned transformations on
training data only. Do not zero-fill unmeasured genes at prediction time.

After selection, `preprocess_counts` applies the existing DEGAS transform to
**each row across the selected genes**:
`normalize_scale(1.5 ** log2(counts + 1))`. Do not preprocess the complete gene
universe and then subset, or normalize each new assay differently. Inputs in
R genes-by-samples orientation must be transposed. For actual log2(count+1)
inputs the preprocessing API supports `already_log=True`; the CSV example
expects raw counts. Already-normalized or negative matrices require a deliberate
compatible preprocessing path, not another count transform.

## Choose the split design

- `--split study`: reserve `ceil(0.1 * number_of_studies)` complete studies
  (at least one); leave each remaining study out once for validation. This
  example therefore needs at least three studies, each with both classes.
  With unequal study sizes, the final holdout is not necessarily 10% of patients.
- `--split patient`: reserve approximately 10% of unique patients, stratified by
  outcome, then use stratified patient K-fold validation (`--folds 5`). Studies
  can appear on both sides; this estimates within-cohort patient generalization,
  not transfer to unseen studies. It is useful with few or very unequal studies.

Change the reserved fraction with `--holdout-fraction` and the prespecified split
seed with `--seed`. Tiny fractions/cohorts may not supply both outcome classes;
the workflow fails explicitly. Review membership and balance rather than trying
seeds until the performance looks favorable. Native helpers also keep repeated
rows of a patient together, but this example and metric table require one score
per patient. Aggregate technical replicates with a prespecified rule first.

Reusable API (indices refer to the supplied metadata row order):

```python
from DEGAS_python.validation import holdout_split, validation_folds

development, final_test = holdout_split(metadata, mode='patient', fraction=0.1)
dev_metadata = metadata.iloc[development].reset_index(drop=True)
for train_local, valid_local in validation_folds(dev_metadata, mode='patient'):
    train = development[train_local]
    valid = development[valid_local]
    # Select genes on train; preprocess; fit DEGAS on train; predict valid.
```

The existing `tot_folds` option concerns cell subsampling and is **not** this
patient validation design. The example sets `tot_folds=1, fold=-1` and supplies
explicit bulk training/evaluation matrices to `run_model`. Its new optional
`return_model=True` returns `(directory, fitted_model)` for prediction/SHAP;
old callers still receive only a directory.

## One-SE selection and score construction

1. Within each size and validation fold, average raw class-1 probabilities over
   seeds. Combine these out-of-fold predictions to obtain exactly one prediction
   per development patient and size.
2. Compute AUROC per study and take the **equal-study mean**, including when
   folds were split at the patient level. Require both classes and at least two
   patients per study/class for the bootstrap. Do not silently drop a study.
3. Bootstrap patients within study × diagnosis, using the same resampled
   patients for every size. The default is 20,000 replicates; SE is the sample
   SD of bootstrap equal-study means. This conditions on fitted models and the
   observed studies; it does not capture all training or study-level uncertainty.
4. Let `best` have the largest mean AUROC (smallest size breaks an exact tie).
   Retain **all** sizes with `mean_AUROC >= mean_AUROC[best] - SE[best]`.
   This is an ensemble rule, not the conventional smallest-model one-SE rule.
5. Refit each retained size using development patients only. For the chosen
   reference population, average seeds within size, compute average-tie rank/N
   within each size, then average those ranks equally across retained sizes.
   Do not rank the average raw score or re-rank the final mean.
6. Evaluate the frozen ensemble on the untouched final holdout once. Holdout
   labels never select sizes, features, seeds, epochs, backgrounds or cutoffs.

The cell reference is **all supplied reference cells across donors/types**, not
one donor or cell type at a time. Filtering after ranking preserves those
reference ranks. The example ranks final test patients separately over their
entire holdout cohort. Both scores depend on their declared reference cohort;
adding observations can change ranks. A percentile ensemble is not a calibrated
probability: the example reports its AUROC and average precision, and reports
log loss/Brier/balanced accuracy only for raw per-size probability predictions.
A fixed deployment rank reference must be specified separately.

If there is no independent final test set, call `validation_folds(metadata)`
directly and use its out-of-fold predictions for tuning; call the resulting
performance post-selection/exploratory, not an unbiased final-test estimate.

## Native evaluation figures

```python
from DEGAS_python.validation import select_one_se
from DEGAS_python.evaluation_plots import plot_validation

# One seed-averaged OOF prediction per patient and size:
# patient_id, study, label, size, score
metrics, selection = select_one_se(predictions, n_bootstrap=20000, seed=42)
auroc_figure, metric_figure = plot_validation(
    metrics, selection, title='My cohort: bulk patient validation',
    output_dir='figures/my_cohort',
)
```

The figures follow the T2D `SC_plotting_starter.R` evaluation design: faint study
curves, a black equal-study mean, a gold ring at the best AUROC, a dashed one-SE
cutoff, and a vertical marker from the best mean down by exactly one SE. Green
squares show retained sizes. SE is **not** a 95% confidence interval. The other
figure shows AUROC, average precision, log loss, Brier score, and balanced
accuracy (threshold 0.5), with PNG/PDF files and the supporting CSV tables.

Run separately for broad and panel-intersection gene universes. Never reuse one
universe's cutoff for another. These are held-out **bulk patient** metrics,
not validation of cellular diagnosis or spatial biology. SC/Xenium pseudobulk
transfer comparisons in the original plotting script remain separate,
descriptive analyses; cells/neighborhoods are not independent patients.

## Direct SHAP

See [the direct SHAP guide](../../DIRECT_SHAP.md). `--shap` explains a small,
fixed set of 16 reference cells for every retained final model, with 32 common
background cells (or all available when fewer). Increase/stratify these budgets
for your question. These synthetic defaults are not a biologically selected
reference. Every seed for a given size uses the same background and observations.

The tutorial checks predictions against the original DEGAS evaluation output,
writes per-seed SHAP values and seed-mean explanations, and reports mean/p95
absolute additivity residuals. SHAP uses the trained feature extractor and
classifier directly, with no fitted surrogate. It explains the **raw class-1
probability as a function of normalized input**, not the rank ensemble or raw
count interventions. Attributions from different gene sets remain separate.

## Outputs and verification

`fold_membership.csv`, `validation_predictions.csv`, `feature_size_selection.csv`
and `metrics_by_study.csv` record selection. `auroc_selection.{png,pdf}` and
`evaluation_metrics.{png,pdf}` are the figures. `final/size_*/` contains ordered
gene lists, per-seed settings/checkpoints/losses, and optional SHAP artifacts.
The top-level raw-score, rank-score and `holdout_metrics.json` files retain final
evaluation separately. `run.json` records options and dependency versions.
Preserve the input files and their provenance with the run.

Run tests from the repository root:

```bash
python -m pip install pytest
PYTHONPATH=DEGAS_python python -m pytest DEGAS_python/tests tests -q
```

### Figure previews (synthetic data only)

These previews come from the synthetic example with `--sizes 10 20 40 --seeds 2
--iters 20 --bootstrap 1000`. They illustrate the layout, not biological performance.

![Synthetic AUROC selection plot](figures/auroc_selection.png)
![Synthetic five-metric plot](figures/evaluation_metrics.png)
