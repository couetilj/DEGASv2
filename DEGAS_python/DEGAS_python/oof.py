"""Patient-level out-of-fold aggregation and calibration for DEGAS.

The functions in this module are independent of torch so that fold integrity
and calibration can be tested without training a neural network.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logit
from scipy.stats import rankdata


def patient_group_folds(labels, groups):
    """Return deterministic leave-one-group-out patient indices."""
    labels = np.asarray(labels).reshape(-1)
    groups = np.asarray(groups).reshape(-1)
    if len(labels) != len(groups) or not len(labels):
        raise ValueError("patient labels and groups must have the same nonzero length")
    if pd.isna(groups).any():
        raise ValueError("patient groups cannot be missing")
    folds = []
    for fold, heldout in enumerate(sorted(np.unique(groups).tolist(), key=str)):
        test = np.flatnonzero(groups == heldout)
        train = np.flatnonzero(groups != heldout)
        if not len(test) or not len(train):
            raise ValueError("every patient fold requires nonempty train and test sets")
        if len(np.unique(labels[train])) < 2:
            raise ValueError("every patient training fold requires both outcome classes")
        if np.intersect1d(train, test).size:
            raise AssertionError("a patient occurs in both train and test")
        folds.append(
            {
                "fold": fold,
                "heldout_group": heldout,
                "train_indices": train,
                "test_indices": test,
            }
        )
    if np.unique(np.concatenate([x["test_indices"] for x in folds])).size != len(labels):
        raise AssertionError("patient folds do not cover every patient exactly once")
    return folds


def _calibration_design(probability):
    probability = np.asarray(probability, dtype=float)
    if not np.isfinite(probability).all() or np.any((probability < 0) | (probability > 1)):
        raise ValueError("probabilities must be finite and lie in [0,1]")
    clipped = np.clip(probability, 1e-6, 1 - 1e-6)
    return np.column_stack([np.ones(len(clipped)), logit(clipped)])


def fit_ridge_logistic(probability, label, penalty=1.0):
    """Fit Platt-style intercept/slope calibration with a ridge slope penalty."""
    design = _calibration_design(probability)
    label = np.asarray(label, dtype=float).reshape(-1)
    if len(label) != len(design) or not set(np.unique(label)).issubset({0.0, 1.0}):
        raise ValueError("labels must be binary and aligned to probabilities")
    if len(np.unique(label)) != 2:
        raise ValueError("calibration requires both outcome classes")
    if not np.isfinite(penalty) or penalty < 0:
        raise ValueError("ridge penalty must be finite and nonnegative")

    def objective(parameters):
        eta = design @ parameters
        likelihood = np.logaddexp(0, eta).sum() - np.dot(label, eta)
        return likelihood + 0.5 * penalty * parameters[1] ** 2

    result = minimize(objective, np.array([0.0, 1.0]), method="BFGS")
    if not result.success or not np.isfinite(result.x).all():
        raise RuntimeError("ridge-logistic calibration failed to converge")
    return {"intercept": float(result.x[0]), "slope": float(result.x[1]), "penalty": float(penalty)}


def apply_calibrator(probability, parameters):
    design = _calibration_design(probability)
    values = expit(design @ np.array([parameters["intercept"], parameters["slope"]]))
    return np.clip(values, 0.0, 1.0)


def binary_metrics(label, probability):
    label = np.asarray(label, dtype=int)
    probability = np.asarray(probability, dtype=float)
    if len(label) != len(probability) or len(np.unique(label)) != 2:
        raise ValueError("binary metrics require aligned predictions and both classes")
    predicted = probability >= 0.5
    sensitivity = np.mean(predicted[label == 1])
    specificity = np.mean(~predicted[label == 0])
    positive = probability[label == 1]
    negative = probability[label == 0]
    ranks = rankdata(np.concatenate([positive, negative]), method="average")
    n_positive = len(positive)
    auc = (ranks[:n_positive].sum() - n_positive * (n_positive + 1) / 2) / (
        n_positive * len(negative)
    )
    return {
        "n_patients": int(len(label)),
        "balanced_accuracy": float((sensitivity + specificity) / 2),
        "roc_auc": float(auc),
        "brier_score": float(np.mean((probability - label) ** 2)),
    }


def reliability_table(frame, probability_column, bins=5):
    if bins < 2:
        raise ValueError("reliability table requires at least two bins")
    work = frame[["label", probability_column]].copy()
    work["reliability_bin"] = pd.qcut(
        work[probability_column], q=min(bins, len(work)), duplicates="drop"
    )
    result = (
        work.groupby("reliability_bin", observed=True)
        .agg(
            mean_prediction=(probability_column, "mean"),
            observed_fraction=("label", "mean"),
            n_patients=("label", "size"),
        )
        .reset_index(drop=True)
    )
    result.insert(0, "probability", probability_column)
    result.insert(1, "bin", np.arange(1, len(result) + 1))
    return result


def cross_fit_calibration(ensemble, group_column="heldout_group", penalty=1.0, bins=5):
    """Calibrate each group using only labels from the other OOF groups."""
    required = {"patient_id", "label", "raw_probability", group_column}
    missing = required - set(ensemble)
    if missing:
        raise ValueError("OOF ensemble omits {}".format(sorted(missing)))
    if ensemble.patient_id.astype(str).duplicated().any():
        raise ValueError("OOF ensemble must contain one row per patient")
    result = ensemble.copy()
    result["calibrated_probability"] = np.nan
    result["calibration_train_groups"] = ""
    parameter_rows = []
    groups = sorted(result[group_column].unique().tolist(), key=str)
    for heldout in groups:
        test = result[group_column] == heldout
        train = ~test
        train_groups = sorted(result.loc[train, group_column].astype(str).unique())
        parameters = fit_ridge_logistic(
            result.loc[train, "raw_probability"], result.loc[train, "label"], penalty
        )
        result.loc[test, "calibrated_probability"] = apply_calibrator(
            result.loc[test, "raw_probability"], parameters
        )
        result.loc[test, "calibration_train_groups"] = ";".join(train_groups)
        parameter_rows.append(
            {"evaluation_group": heldout, "training_groups": ";".join(train_groups), **parameters}
        )
    if result.calibrated_probability.isna().any():
        raise AssertionError("cross-fitted calibration omitted patients")

    metrics = []
    for name in ("raw_probability", "calibrated_probability"):
        row = {"probability": name, "evaluation_group": "ALL"}
        row.update(binary_metrics(result.label, result[name]))
        diagnostic = fit_ridge_logistic(result[name], result.label, penalty=1e-8)
        row["calibration_intercept"] = diagnostic["intercept"]
        row["calibration_slope"] = diagnostic["slope"]
        metrics.append(row)
        for group, part in result.groupby(group_column, sort=True):
            if len(np.unique(part.label)) < 2:
                continue
            subgroup = {"probability": name, "evaluation_group": group}
            subgroup.update(binary_metrics(part.label, part[name]))
            subgroup["calibration_intercept"] = np.nan
            subgroup["calibration_slope"] = np.nan
            metrics.append(subgroup)

    deployment = fit_ridge_logistic(result.raw_probability, result.label, penalty)
    deployment.update(
        {
            "schema_version": "degas_patient_calibrator/v1",
            "input": "bagged_patient_oof_raw_probability",
            "application_order": "average_matching_ensemble_then_calibrate",
            "n_patients": int(len(result)),
            "groups": [str(x) for x in groups],
        }
    )
    reliability = pd.concat(
        [reliability_table(result, name, bins) for name in ("raw_probability", "calibrated_probability")],
        ignore_index=True,
    )
    return result, pd.DataFrame(metrics), reliability, pd.DataFrame(parameter_rows), deployment


def write_calibration_outputs(directory, ensemble, penalty=1.0, bins=5):
    directory = Path(directory)
    calibrated, metrics, reliability, fold_parameters, deployment = cross_fit_calibration(
        ensemble, penalty=penalty, bins=bins
    )
    calibrated.to_csv(directory / "patient_oof_calibrated.csv", index=False)
    metrics.to_csv(directory / "patient_oof_metrics.csv", index=False)
    reliability.to_csv(directory / "patient_oof_reliability.csv", index=False)
    fold_parameters.to_csv(directory / "patient_oof_calibration_folds.csv", index=False)
    (directory / "patient_calibrator.json").write_text(
        json.dumps(deployment, indent=2, sort_keys=True) + "\n"
    )
    return calibrated
