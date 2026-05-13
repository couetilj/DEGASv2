"""Tests for patient-grouped CV split in STSCDataset.

The default behavior (no groups) is preserved: each fold gets a random
1/N slice of the data with no group awareness. When groups are provided,
sklearn.GroupKFold ensures no group is split across folds — the property
that matters for SC/ST data where many cells/spots come from the same
patient.
"""
from __future__ import annotations

import numpy as np
import pytest

from DEGAS_python.datasets.dataset import STSCDataset


def test_no_groups_preserves_default_random_shuffle() -> None:
    """When groups is None, fold splitting falls back to the original random shuffle."""
    rng = np.random.default_rng(0)
    n_samples = 90
    X = rng.standard_normal((n_samples, 20)).astype(np.float32)
    y = (rng.random(n_samples) > 0.5).astype(np.int8)

    fold_indices_by_fold = []
    for fold in range(3):
        ds = STSCDataset(X, y, fold=fold, tot_folds=3, tot_iters=10,
                          batch_size=8, phase='train')
        fold_indices_by_fold.append(np.sort(ds.fold_indices))

    # Three disjoint slices that together cover all samples (modulo fold_size truncation)
    all_kept = np.concatenate(fold_indices_by_fold)
    assert len(np.unique(all_kept)) == len(all_kept), 'folds overlap'


def test_grouped_split_no_group_spans_two_folds() -> None:
    """Core contract: with groups supplied, no group's samples cross folds."""
    rng = np.random.default_rng(0)
    n_groups = 30
    samples_per_group = 3
    n_samples = n_groups * samples_per_group
    groups = np.repeat(np.arange(n_groups), samples_per_group)
    X = rng.standard_normal((n_samples, 20)).astype(np.float32)
    y = (rng.random(n_samples) > 0.5).astype(np.int8)

    # Track which fold each group landed in
    group_to_fold: dict[int, int] = {}
    all_fold_indices = []
    for fold in range(3):
        ds = STSCDataset(X, y, fold=fold, tot_folds=3, tot_iters=10,
                          batch_size=8, phase='train', groups=groups)
        all_fold_indices.append(ds.fold_indices)
        for idx in ds.fold_indices:
            g = int(groups[idx])
            if g in group_to_fold and group_to_fold[g] != fold:
                pytest.fail(f'group {g} appears in both fold {group_to_fold[g]} and {fold}')
            group_to_fold[g] = fold

    # All groups should have been assigned to a fold
    assert len(group_to_fold) == n_groups, f'only {len(group_to_fold)} groups covered'

    # Folds are disjoint at the sample level too
    pooled = np.concatenate(all_fold_indices)
    assert len(np.unique(pooled)) == len(pooled), 'sample-level overlap'


def test_grouped_split_uneven_group_sizes() -> None:
    """Groups with very different sizes should still split cleanly."""
    rng = np.random.default_rng(0)
    # group 0 has 50 samples, groups 1..9 have 5 each → 100 total
    sizes = [50] + [5] * 9
    groups = np.concatenate([np.full(s, i) for i, s in enumerate(sizes)])
    n_samples = len(groups)
    X = rng.standard_normal((n_samples, 20)).astype(np.float32)
    y = (rng.random(n_samples) > 0.5).astype(np.int8)

    group_to_fold: dict[int, int] = {}
    for fold in range(3):
        ds = STSCDataset(X, y, fold=fold, tot_folds=3, tot_iters=10,
                          batch_size=8, phase='train', groups=groups)
        for idx in ds.fold_indices:
            g = int(groups[idx])
            assert g not in group_to_fold or group_to_fold[g] == fold,                 f'group {g} split across folds'
            group_to_fold[g] = fold
    assert len(group_to_fold) == 10


def test_too_few_groups_raises() -> None:
    """Fewer unique groups than folds is not splittable — should raise."""
    groups = np.array([0, 0, 1, 1, 1])  # 2 groups
    X = np.zeros((5, 10), dtype=np.float32)
    y = np.zeros(5, dtype=np.int8)
    with pytest.raises(ValueError, match='at least tot_folds'):
        STSCDataset(X, y, fold=0, tot_folds=3, tot_iters=5, batch_size=2,
                     phase='train', groups=groups)


def test_groups_length_mismatch_raises() -> None:
    X = np.zeros((10, 5), dtype=np.float32)
    y = np.zeros(10, dtype=np.int8)
    with pytest.raises(ValueError, match='groups length'):
        STSCDataset(X, y, fold=0, tot_folds=2, tot_iters=5, batch_size=2,
                     phase='train', groups=np.arange(7))
