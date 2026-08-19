"""Tests for native patient-level OOF calibration outputs."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from DEGAS_python.oof import (
    apply_calibrator,
    binary_metrics,
    cross_fit_calibration,
    fit_ridge_logistic,
    patient_group_folds,
    write_calibration_outputs,
)


def synthetic_ensemble():
    groups = np.repeat(["study_a", "study_b", "study_c", "study_d"], 6)
    labels = np.tile([0, 0, 0, 1, 1, 1], 4)
    raw = np.array([0.20, 0.25, 0.35, 0.65, 0.75, 0.80] * 4)
    raw += np.repeat([-0.05, 0.02, 0.06, -0.02], 6)
    return pd.DataFrame(
        {
            "patient_id": ["p{:02d}".format(i) for i in range(len(labels))],
            "label": labels,
            "heldout_group": groups,
            "raw_probability": np.clip(raw, 0.01, 0.99),
            "raw_probability_sd": 0.03,
            "n_submodels": 5,
        }
    )


def test_patient_group_folds_exclude_complete_studies():
    labels = np.tile([0, 1], 8)
    groups = np.repeat(["a", "b", "c", "d"], 4)
    folds = patient_group_folds(labels, groups)
    assert len(folds) == 4
    held = []
    for fold in folds:
        assert not np.intersect1d(fold["train_indices"], fold["test_indices"]).size
        assert set(groups[fold["test_indices"]]) == {fold["heldout_group"]}
        held.extend(fold["test_indices"].tolist())
    assert sorted(held) == list(range(len(labels)))


def test_cross_fit_calibration_never_uses_heldout_group():
    frame = synthetic_ensemble()
    calibrated, metrics, reliability, parameters, deployment = cross_fit_calibration(frame)
    assert calibrated.calibrated_probability.between(0, 1).all()
    for _, row in calibrated.iterrows():
        assert row.heldout_group not in row.calibration_train_groups.split(";")
    assert set(parameters.evaluation_group) == set(frame.heldout_group)
    assert set(metrics.probability) == {"raw_probability", "calibrated_probability"}
    assert set(reliability.probability) == {"raw_probability", "calibrated_probability"}
    assert deployment["application_order"] == "average_matching_ensemble_then_calibrate"


def test_calibration_and_metrics_are_finite():
    frame = synthetic_ensemble()
    parameters = fit_ridge_logistic(frame.raw_probability, frame.label)
    values = apply_calibrator(frame.raw_probability, parameters)
    metrics = binary_metrics(frame.label, values)
    assert np.isfinite(list(parameters.values())).all()
    assert np.isfinite(list(metrics.values())).all()
    assert 0 <= metrics["roc_auc"] <= 1
    assert 0 <= metrics["brier_score"] <= 1


def test_write_calibration_outputs(tmp_path):
    frame = synthetic_ensemble()
    write_calibration_outputs(tmp_path, frame)
    expected = {
        "patient_oof_calibrated.csv",
        "patient_oof_metrics.csv",
        "patient_oof_reliability.csv",
        "patient_oof_calibration_folds.csv",
        "patient_calibrator.json",
    }
    assert expected == {path.name for path in tmp_path.iterdir()}
    payload = json.loads((tmp_path / "patient_calibrator.json").read_text())
    assert payload["schema_version"] == "degas_patient_calibrator/v1"
