"""Tests for patient-grouped fold subsampling on the BULK side (PatDataset).

DEGAS's SC/ST side has its own random-shuffle fold subsampling (STSCDataset),
which we intentionally do NOT patient-group. Per-cancer scientific motivation
for SC/ST is to bag across distinct cell subsets, not to hold patients out.

The BULK side, by contrast, can have matched T+N aliquots from the same
patient (e.g. TCGA-COAD  + ). When that is the case, fold
subsampling must respect patient boundaries, or matched pairs get split
across folds and the model effectively trains-and-evaluates on the same
patient.

This is OFF by default (groups=None preserves the original behavior of
every fold seeing the full bulk set). When groups is supplied and
tot_folds > 1, each fold is restricted to one GroupKFold patient slice.
"""
from __future__ import annotations

import numpy as np
import pytest

from DEGAS_python.datasets.dataset import PatDataset


def test_no_groups_uses_full_bulk_every_fold() -> None:
    """groups=None preserves original DEGAS behavior: full bulk per fold."""
    rng = np.random.default_rng(0)
    n_samples = 60
    X = rng.standard_normal((n_samples, 20)).astype(np.float32)
    y = np.zeros((n_samples, 2), dtype=np.float32)
    y[:30, 0] = 1.0
    y[30:, 1] = 1.0
    ds = PatDataset(X, y, fold=0, tot_folds=3, tot_iters=10,
                     batch_size=8, phase='train')
    assert ds.fold_indices is None, 'no-groups path should not set fold_indices'


def test_grouped_split_no_group_spans_two_folds() -> None:
    """Core contract: with groups, no patient's aliquots cross folds."""
    rng = np.random.default_rng(0)
    n_pat = 30
    aliquots_per_pat = 2   # simulates matched T+N
    n_samples = n_pat * aliquots_per_pat
    groups = np.repeat(np.arange(n_pat), aliquots_per_pat)
    X = rng.standard_normal((n_samples, 20)).astype(np.float32)
    y = np.zeros((n_samples, 2), dtype=np.float32)
    y[:n_samples // 2, 0] = 1.0
    y[n_samples // 2:, 1] = 1.0

    group_to_fold: dict[int, int] = {}
    all_fold_indices = []
    for fold in range(3):
        ds = PatDataset(X, y, fold=fold, tot_folds=3, tot_iters=10,
                         batch_size=8, phase='train', groups=groups)
        assert ds.fold_indices is not None
        all_fold_indices.append(ds.fold_indices)
        for idx in ds.fold_indices:
            g = int(groups[idx])
            if g in group_to_fold and group_to_fold[g] != fold:
                pytest.fail(f'patient {g} appears in both fold {group_to_fold[g]} and {fold}')
            group_to_fold[g] = fold

    assert len(group_to_fold) == n_pat
    pooled = np.concatenate(all_fold_indices)
    assert len(np.unique(pooled)) == len(pooled), 'sample-level overlap across folds'


def test_grouped_split_matched_TN_pairs_stay_together() -> None:
    """Concrete TCGA scenario: each patient contributes one tumor + one normal."""
    rng = np.random.default_rng(0)
    n_pat = 24
    # patient k contributes samples 2k (tumor) and 2k+1 (normal)
    n_samples = 2 * n_pat
    groups = np.repeat(np.arange(n_pat), 2)
    X = rng.standard_normal((n_samples, 50)).astype(np.float32)
    # All aliquots from a patient share the patient's label
    pat_labels = (rng.random(n_pat) > 0.5).astype(int)
    y_int = np.repeat(pat_labels, 2)
    y = np.zeros((n_samples, 2), dtype=np.float32)
    y[y_int == 0, 0] = 1.0
    y[y_int == 1, 1] = 1.0

    for fold in range(4):
        ds = PatDataset(X, y, fold=fold, tot_folds=4, tot_iters=10,
                         batch_size=8, phase='train', groups=groups)
        kept_patients = set(int(groups[i]) for i in ds.fold_indices)
        for p in kept_patients:
            tumor_idx = 2 * p
            normal_idx = 2 * p + 1
            assert tumor_idx in ds.fold_indices and normal_idx in ds.fold_indices, \
                f'patient {p} fold {fold}: T+N pair was split'


def test_cox_model_path_with_groups() -> None:
    """Cox path uses self.status for the y arg — verify it still splits cleanly."""
    rng = np.random.default_rng(0)
    n_pat = 20
    n_samples = n_pat * 2
    groups = np.repeat(np.arange(n_pat), 2)
    X = rng.standard_normal((n_samples, 30)).astype(np.float32)
    # [time, event] schema
    times = rng.integers(100, 2000, size=n_samples)
    events = rng.integers(0, 2, size=n_samples)
    y_cox = np.column_stack([times, events]).astype(np.float32)

    for fold in range(3):
        ds = PatDataset(X, y_cox, fold=fold, tot_folds=3, tot_iters=10,
                         batch_size=8, phase='train',
                         model_type='BlankCox', groups=groups)
        assert ds.fold_indices is not None
        kept_groups = set(int(groups[i]) for i in ds.fold_indices)
        # No partial patient
        for g in kept_groups:
            members = np.where(groups == g)[0]
            for m in members:
                assert m in ds.fold_indices, f'Cox fold {fold} split patient {g}'


def test_too_few_groups_raises() -> None:
    groups = np.array([0, 0, 1, 1, 1])
    X = np.zeros((5, 10), dtype=np.float32)
    y = np.zeros((5, 2), dtype=np.float32)
    y[:, 0] = 1.0
    with pytest.raises(ValueError, match='at least tot_folds'):
        PatDataset(X, y, fold=0, tot_folds=3, tot_iters=5, batch_size=2,
                    phase='train', groups=groups)


def test_groups_length_mismatch_raises() -> None:
    X = np.zeros((10, 5), dtype=np.float32)
    y = np.zeros((10, 2), dtype=np.float32)
    y[:, 0] = 1.0
    with pytest.raises(ValueError, match='groups length'):
        PatDataset(X, y, fold=0, tot_folds=2, tot_iters=5, batch_size=2,
                    phase='train', groups=np.arange(7))
