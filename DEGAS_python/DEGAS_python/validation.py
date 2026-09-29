"""Grouped binary-classifier validation, one-SE selection and rank ensembles.

Independent of torch. Indices always refer to the supplied metadata row order.
"""
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import (roc_auc_score, average_precision_score, log_loss,
                             brier_score_loss, balanced_accuracy_score)
from sklearn.model_selection import StratifiedKFold, train_test_split

METRICS = ('AUROC', 'average_precision', 'log_loss', 'brier', 'balanced_accuracy')


def _patients(metadata):
    required = ['patient_id', 'study', 'label']
    if not set(required).issubset(metadata.columns) or metadata[required].isna().any().any():
        raise ValueError('Metadata requires nonmissing patient_id, study, label')
    if not len(metadata) or not set(metadata.label.unique()).issubset({0, 1}):
        raise ValueError('Use binary labels 0 and 1')
    if (metadata.groupby('patient_id')[['study', 'label']].nunique() > 1).any().any():
        raise ValueError('Each globally unique patient_id must have one study and label')
    return metadata[required].drop_duplicates('patient_id').reset_index(drop=True)


def _split(metadata, patients, heldout_ids):
    mask = metadata.patient_id.isin(heldout_ids).to_numpy()
    train, test = np.flatnonzero(~mask), np.flatnonzero(mask)
    for ids in (metadata.iloc[train].patient_id, metadata.iloc[test].patient_id):
        part = patients[patients.patient_id.isin(ids)]
        if set(part.label) != {0, 1}:
            raise ValueError('Each training and validation partition needs both classes; '
                             'change the design, not the seed after viewing model performance')
    return train, test


def holdout_split(metadata, *, mode='study', fraction=0.1, seed=42):
    """Reserve patients (stratified) or whole studies (random, at least one).

    Study mode holds ceil(fraction * number_of_studies) studies, not 10% of
    rows. Repeated rows of a patient always stay together. No outcome-based
    search for a favorable study split is performed.
    """
    patients = _patients(metadata)
    if not 0 < fraction < 1:
        raise ValueError('fraction must lie strictly between 0 and 1')
    if mode == 'patient':
        _, ids = train_test_split(patients.patient_id, test_size=fraction,
                                  stratify=patients.label, random_state=seed)
    elif mode == 'study':
        studies = sorted(patients.study.unique(), key=str)
        n = int(np.ceil(fraction * len(studies)))
        if n >= len(studies):
            raise ValueError('Study holdout requires at least two studies')
        chosen = np.random.default_rng(seed).choice(len(studies), n, replace=False)
        ids = patients.loc[patients.study.isin([studies[i] for i in chosen]), 'patient_id']
    else:
        raise ValueError("mode must be 'study' or 'patient'")
    return _split(metadata, patients, ids)


def validation_folds(metadata, *, mode='study', n_splits=5, seed=42):
    """Leave-one-study-out or stratified patient K-fold splits.

    Call on development metadata AFTER reserving the final holdout. Study mode
    uses every remaining study once. Patient mode can mix studies across folds.
    """
    patients = _patients(metadata)
    if mode == 'study':
        heldout = [patients.loc[patients.study == s, 'patient_id']
                   for s in sorted(patients.study.unique(), key=str)]
    elif mode == 'patient':
        if n_splits < 2 or patients.label.value_counts().min() < n_splits:
            raise ValueError('Need at least n_splits patients in each class')
        splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        heldout = [patients.iloc[test].patient_id for _, test in
                   splitter.split(patients, patients.label)]
    else:
        raise ValueError("mode must be 'study' or 'patient'")
    return [_split(metadata, patients, ids) for ids in heldout]


def classification_metrics(label, score):
    """Metrics on raw class-1 probabilities; balanced accuracy uses 0.5."""
    y, p = np.asarray(label), np.asarray(score, dtype=float)
    if y.ndim != 1 or p.shape != y.shape or set(y) != {0, 1}:
        raise ValueError('Metrics require aligned 1D binary labels with both classes')
    if not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError('Expected finite probabilities in [0, 1]')
    return dict(AUROC=roc_auc_score(y, p), average_precision=average_precision_score(y, p),
                log_loss=log_loss(y, p, labels=[0, 1]), brier=brier_score_loss(y, p),
                balanced_accuracy=balanced_accuracy_score(y, p >= 0.5))


def select_one_se(predictions, *, n_bootstrap=20000, seed=42):
    """Return study metrics and one-SE size selection from OOF predictions.

    Input: one seed-averaged probability per patient and size, columns
    patient_id, study, label, size, score. All sizes must cover the same patients.
    Equal-study mean AUROC; paired patient bootstrap within study x diagnosis.
    SE is conditional on the fitted OOF models and observed studies, not a CI
    or an estimate of training/study-population uncertainty. Seeds are not
    independent patient replicates. Requires >=2 patients per study/class.
    """
    _patients(predictions)
    if not {'size', 'score'}.issubset(predictions.columns):
        raise ValueError('Predictions require size and score')
    if predictions[['size', 'score']].isna().any().any() or predictions.duplicated(['patient_id', 'size']).any():
        raise ValueError('Supply one finite seed-mean prediction per patient and size')
    sizes = np.sort(predictions['size'].unique())
    if np.any(sizes < 2) or np.any(sizes != sizes.astype(int)):
        raise ValueError('Gene-set sizes must be integers >=2')
    if not isinstance(n_bootstrap, int) or n_bootstrap < 2:
        raise ValueError('n_bootstrap must be an integer >=2')
    rng = np.random.default_rng(seed)
    rows, means = [], []
    bootstrap = np.zeros((n_bootstrap, len(sizes)))
    for study, part in predictions.groupby('study', sort=True):
        table = part.pivot(index=['patient_id', 'label'], columns='size', values='score').reindex(columns=sizes)
        if table.isna().any().any():
            raise ValueError('All sizes must cover exactly the same patients')
        y = table.index.get_level_values('label').to_numpy()
        scores = table.to_numpy()
        for j, size in enumerate(sizes):
            rows.append(dict(study=study, size=int(size), n_patients=len(y),
                             **classification_metrics(y, scores[:, j])))
        positive, negative = scores[y == 1], scores[y == 0]
        if min(len(positive), len(negative)) < 2:
            raise ValueError('Bootstrap SE needs >=2 patients per study and class')
        # Multinomial bootstrap weights avoid repeatedly computing sorted AUCs.
        wins = (positive[:, None, :] > negative[None, :, :]).astype(float)
        wins += 0.5 * (positive[:, None, :] == negative[None, :, :])
        means.append(wins.mean(axis=(0, 1)))
        for start in range(0, n_bootstrap, 256):
            count = min(256, n_bootstrap - start)
            a = rng.multinomial(len(positive), np.full(len(positive), 1 / len(positive)), count) / len(positive)
            b = rng.multinomial(len(negative), np.full(len(negative), 1 / len(negative)), count) / len(negative)
            bootstrap[start:start+count] += np.einsum('bi,ijs,bj->bs', a, wins, b, optimize=True)
    average = np.mean(means, axis=0)
    se = (bootstrap / len(means)).std(axis=0, ddof=1)
    best = int(np.argmax(average))  # deterministic: smallest size on exact ties
    cutoff = average[best] - se[best]
    selection = pd.DataFrame(dict(size=sizes, mean_AUROC=average, SE=se,
                                   cutoff=cutoff, selected=average >= cutoff))
    return pd.DataFrame(rows), selection


def average_rank_scores(scores_by_size, selected_sizes):
    """Mean of per-size average-tie ranks/N over one declared reference cohort.

    Rows are the SAME observations in the SAME order; columns are integer sizes.
    Average seeds before calling. No re-ranking after averaging across sizes.
    This cohort-relative score is not a calibrated probability or a pointwise
    function suitable for the direct probability SHAP wrapper.
    """
    selected = list(selected_sizes)
    if not selected or len(set(selected)) != len(selected):
        raise ValueError('selected_sizes must be nonempty and unique')
    if not scores_by_size.index.is_unique or not scores_by_size.columns.is_unique:
        raise ValueError('Observation IDs and size columns must be unique')
    x = scores_by_size.loc[:, selected].to_numpy(dtype=float)
    if not len(x) or not np.isfinite(x).all():
        raise ValueError('Scores must be finite and nonempty')
    ranks = rankdata(x, method='average', axis=0) / len(x)
    return pd.Series(ranks.mean(axis=1), index=scores_by_size.index, name='average_rank_score')
