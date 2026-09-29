# Verification record

Local CPU checks on September 29, 2026 (macOS arm64, Python 3.13):

- 38 tests passed, covering patient/study separation, repeat-patient grouping,
  invalid designs, equal-study weighting with unequal cohort sizes, paired
  bootstrap SE, complete size coverage, tie ranks, plot exports, direct SHAP,
  preprocessing, existing OOF calibration, and singleton classifier/Cox batches.
- All-cell SHAP checks cover the empirical low-risk quartile (including ties),
  population/custom backgrounds, reproducible independent sampling caps, all 100
  synthetic cells across two seeds, identical background IDs, checkpoint prediction
  parity, chunk reuse on restart and seed-mean attribution equality.
- New workflow tests cover automatic grids through 16,561 genes, exact-list
  order/availability, training-only feature selection, all-seed pooling,
  incomplete/corrupt/stale result rejection, duplicate IDs, changed inputs,
  idempotent retries and generated SLURM commands/dependencies/resources.
- Real independent Python worker processes trained an exact list concurrently;
  their results pooled through final holdout scoring and SHAP. A separate real
  training task used 3,000 input genes. Generated bash scripts passed syntax checks.
- The current local CLI completed validation, final fitting and all-cell SHAP for
  two candidate sizes, with 40-cell chunks and 32 integration samples.
- Before the all-cell default change, a full synthetic empirical run evaluated 12 sizes through 3,000 genes across
  three patient folds (36 fits), pooled them, trained all five retained sizes,
  and completed final holdout scoring and SHAP. Both wide-grid figures were
  visually inspected. This two-iteration smoke run checks execution only.
- Earlier study and patient end-to-end synthetic workflows completed with direct
  SHAP. Each checks native DEGAS versus wrapper prediction parity.
- A single-study CSV input run completed using patient splits, including ID-based
  input alignment, independent-reference metadata and the final holdout.
- The earlier 3-size / 2-seed / 20-iteration baseline workflow retained all three sizes.
  Saved rank scores matched a separate rank/N calculation; saved seed-mean SHAP
  matched averaging individual explanations. Per-seed mean absolute additivity
  residuals ranged from 0.00245 to 0.00487; p95 ranged from 0.00552 to 0.01151
  with 512 integration samples. These are synthetic numerical checks only.
- AUROC and five-metric PNGs were visually inspected; PNG/PDF exports were tested.

Environment: numpy 2.5.3, pandas 3.0.6, scipy 1.18.1, scikit-learn 1.9.1,
matplotlib 3.11.2, torch 2.14.0, shap 0.52.0. Existing upstream dependency
warnings remain. SLURM submission was tested with dry-runs/mocked sbatch; no
actual cluster jobs were submitted. GPU execution and real-cohort performance
were not assessed.

The Python validation GitHub workflow reruns tests and both split-mode smoke
checks. This local record does not imply that remote CI has passed.
