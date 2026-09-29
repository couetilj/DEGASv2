# Verification record

Local CPU checks on September 29, 2026 (macOS arm64, Python 3.13):

- 25 tests passed, covering patient/study separation, repeat-patient grouping,
  invalid designs, equal-study weighting with unequal cohort sizes, paired
  bootstrap SE, complete size coverage, tie ranks, plot exports, direct SHAP,
  preprocessing, existing OOF calibration, and singleton classifier/Cox batches.
- Both study and patient end-to-end synthetic workflows completed with direct
  SHAP. Each checks native DEGAS versus wrapper prediction parity.
- A single-study CSV input run completed using patient splits, including ID-based
  input alignment, independent-reference metadata and the final holdout.
- A 3-size / 2-seed / 20-iteration study workflow retained all three sizes.
  Saved rank scores matched a separate rank/N calculation; saved seed-mean SHAP
  matched averaging individual explanations. Per-seed mean absolute additivity
  residuals ranged from 0.00245 to 0.00487; p95 ranged from 0.00552 to 0.01151
  with 512 integration samples. These are synthetic numerical checks only.
- AUROC and five-metric PNGs were visually inspected; PNG/PDF exports were tested.

Environment: numpy 2.5.3, pandas 3.0.6, scipy 1.18.1, scikit-learn 1.9.1,
matplotlib 3.11.2, torch 2.14.0, shap 0.52.0. Existing upstream dependency
warnings remain. GPU execution and real-cohort performance were not assessed.

The Python validation GitHub workflow reruns tests and both split-mode smoke
checks. This local record does not imply that remote CI has passed.
