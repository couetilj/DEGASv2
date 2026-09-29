"""Runnable binary-classification tutorial. See README.md in this directory."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import pandas as pd
import torch
import DEGAS_python as degas
from DEGAS_python.direct_shap import DiseaseScore, explain
from DEGAS_python.validation import (holdout_split, validation_folds, select_one_se,
                                     classification_metrics, average_rank_scores)
from DEGAS_python.evaluation_plots import plot_validation


def demo():
    rng = np.random.default_rng(42)
    y = np.tile([0, 1], 60)
    bulk = rng.poisson(10, (120, 40)).astype(float)
    bulk[:, :8] += y[:, None] * rng.poisson(8, (120, 8))
    cells = rng.poisson(10, (100, 40)).astype(float)
    cells[:50, :8] += rng.poisson(8, (50, 8))
    genes = [f'gene_{i}' for i in range(40)]
    meta = pd.DataFrame(dict(patient_id=[f'patient_{i}' for i in range(120)],
                             study=np.repeat(['A', 'B', 'C', 'D'], 30), label=y))
    return (pd.DataFrame(bulk, index=meta.patient_id, columns=genes), meta,
            pd.DataFrame(cells, index=[f'cell_{i}' for i in range(100)], columns=genes))


def read_inputs(folder):
    bulk = pd.read_csv(folder / 'bulk_counts.csv', index_col=0)
    cells = pd.read_csv(folder / 'cell_counts.csv', index_col=0)
    meta = pd.read_csv(folder / 'patients.csv', dtype={'patient_id': str, 'study': str})
    bulk.index = bulk.index.astype(str)
    if not bulk.index.is_unique or not cells.index.is_unique or meta.patient_id.duplicated().any():
        raise ValueError('Example requires one row per patient and unique cell IDs')
    if set(bulk.index) != set(meta.patient_id):
        raise ValueError('Bulk matrix IDs must exactly match patient metadata')
    # A separate reference cohort prevents patient leakage via the target arm.
    cell_meta = pd.read_csv(folder / 'cells.csv', dtype={'cell_id': str, 'patient_id': str})
    cells.index = cells.index.astype(str)
    if (cell_meta[['cell_id', 'patient_id']].isna().any().any()
            or cell_meta.cell_id.duplicated().any() or set(cell_meta.cell_id) != set(cells.index)):
        raise ValueError('cells.csv must identify the known donor of every reference cell')
    if set(cell_meta.patient_id) & set(meta.patient_id):
        raise ValueError('Use an independent cell-reference cohort for this tutorial; '
                         'overlapping donors require fold-specific cell exclusion')
    return bulk.loc[meta.patient_id], meta, cells


def fit(bulk, meta, cells, train, evaluate, sizes, args, folder, explain_final=False):
    """Recompute the example selector on training bulk only, then preprocess."""
    from DEGAS_python.run_module import run_model
    genes = bulk.columns
    # Illustrative training-bulk variance selector. Not DEGAS's R marker selector.
    rank = np.argsort(-np.log1p(bulk.iloc[train].to_numpy()).var(axis=0), kind='stable')
    patient_scores, cell_scores = {}, {}
    for size in sizes:
        selected = genes[rank[:size]].tolist()
        root = folder / f'size_{size}'
        root.mkdir(parents=True)
        (root / 'genes.json').write_text(json.dumps(selected, indent=2))
        x = degas.preprocess_counts(bulk[selected].to_numpy()).astype('float32')
        z = degas.preprocess_counts(cells[selected].to_numpy()).astype('float32')
        ps, cs, explanations = [], [], []
        rng = np.random.default_rng(args.seed)
        bg_ids = rng.choice(len(z), min(32, len(z)), replace=False)
        obs_ids = np.arange(min(16, len(z)))
        for seed in range(args.seeds):
            opt = degas.BlankClass_opt.copy()
            opt.update(save_dir=str(root), seed=seed, fold=-1, tot_folds=1,
                       tot_iters=args.iters, save_freq=args.iters, is_save=True,
                       feature_dim=16, batch_size=32, pat_batch_size=32)
            path, model = run_model(opt, x[train], meta.label.to_numpy()[train], z,
                                    pat_eval_expr_mat=x[evaluate],
                                    pat_eval_lab_mat=meta.label.to_numpy()[evaluate],
                                    return_model=True)
            score = DiseaseScore(model)
            device = next(score.parameters()).device
            with torch.no_grad():
                p = score(torch.as_tensor(x[evaluate], device=device)).cpu().numpy().ravel()
                c = np.concatenate([score(torch.as_tensor(chunk, device=device)).cpu().numpy().ravel()
                                    for chunk in np.array_split(z, max(1, int(np.ceil(len(z)/1024))))])
            native = pd.read_csv(Path(path) / f'low_reso_results_epoch_{args.iters}.csv').sort_values('pid')
            np.testing.assert_allclose(p, native.hazard, atol=1e-6, rtol=1e-5)
            ps.append(p)
            cs.append(c)
            if explain_final and args.shap:
                result = explain(score, z[bg_ids], z[obs_ids], nsamples=args.shap_samples, seed=args.seed)
                np.testing.assert_allclose(result['prediction'], c[obs_ids], atol=1e-6)
                np.savez_compressed(root / f'shap_seed_{seed}.npz', **result,
                                    genes=np.asarray(selected), background=z[bg_ids], inputs=z[obs_ids],
                                    background_ids=cells.index.to_numpy(dtype=str)[bg_ids],
                                    observation_ids=cells.index.to_numpy(dtype=str)[obs_ids])
                explanations.append(result)
        patient_scores[size] = np.mean(ps, axis=0)
        cell_scores[size] = np.mean(cs, axis=0)
        if explanations:
            values = np.mean([e['values'] for e in explanations], axis=0)
            base = np.mean([e['base_value'] for e in explanations])
            prediction = cell_scores[size][obs_ids]
            residual = prediction - base - values.sum(axis=1)
            np.savez_compressed(root / 'shap_seed_mean.npz', values=values, base_value=base,
                                prediction=prediction, residual=residual, genes=np.asarray(selected),
                                inputs=z[obs_ids], observation_ids=cells.index.to_numpy(dtype=str)[obs_ids])
            diagnostics = [{'seed': i, 'mean_abs_residual': float(np.abs(e['residual']).mean()),
                            'p95_abs_residual': float(np.quantile(np.abs(e['residual']), .95))}
                           for i, e in enumerate(explanations)]
            pd.DataFrame(diagnostics).to_csv(root / 'shap_diagnostics.csv', index=False)
    return (pd.DataFrame(patient_scores, index=meta.iloc[evaluate].patient_id),
            pd.DataFrame(cell_scores, index=cells.index))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, help='CSV directory; omit for synthetic demonstration')
    parser.add_argument('--output', type=Path, required=True, help='New directory; existing paths refused')
    parser.add_argument('--split', choices=['study', 'patient'], default='study')
    parser.add_argument('--holdout-fraction', type=float, default=.1)
    parser.add_argument('--folds', type=int, default=5, help='Patient-mode inner folds')
    parser.add_argument('--sizes', type=int, nargs='+', default=[10, 20, 40])
    parser.add_argument('--seeds', type=int, default=5)
    parser.add_argument('--seed', type=int, default=42, help='Split/bootstrap/reference seed')
    parser.add_argument('--iters', type=int, default=300)
    parser.add_argument('--bootstrap', type=int, default=20000)
    parser.add_argument('--shap', action='store_true')
    parser.add_argument('--shap-samples', type=int, default=512)
    args = parser.parse_args()
    if min(args.seeds, args.iters, args.shap_samples) < 1:
        parser.error('seeds, iters, and shap-samples must be positive')
    torch.set_num_threads(1)
    bulk, meta, cells = read_inputs(args.input) if args.input else demo()
    if not bulk.columns.is_unique or not cells.columns.is_unique:
        raise ValueError('Gene IDs must be unique')
    common = bulk.columns.intersection(cells.columns, sort=False)
    bulk, cells = bulk[common], cells[common]
    if any(size < 2 or size > len(common) for size in args.sizes) or len(set(args.sizes)) != len(args.sizes):
        raise ValueError('Each size must be unique and between 2 and the shared-gene count')
    for matrix in (bulk, cells):
        if not len(matrix) or not np.isfinite(matrix.to_numpy()).all() or (matrix.to_numpy() < 0).any():
            raise ValueError('Supply finite, nonnegative raw counts')
    dev, test = holdout_split(meta, mode=args.split, fraction=args.holdout_fraction, seed=args.seed)
    development = meta.iloc[dev].reset_index(drop=True)
    folds = validation_folds(development, mode=args.split, n_splits=args.folds, seed=args.seed)
    args.output.mkdir(parents=True, exist_ok=False)
    versions = {p: importlib.metadata.version(p) for p in ['numpy', 'pandas', 'scipy', 'scikit-learn', 'torch', 'matplotlib']}
    if args.shap:
        versions['shap'] = importlib.metadata.version('shap')
    (args.output / 'run.json').write_text(json.dumps(dict(options={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}, versions=versions), indent=2))
    membership = meta.copy()
    membership['partition'] = 'development'
    membership.loc[test, 'partition'] = 'final_holdout'
    membership['validation_fold'] = -1
    tables = []
    for fold, (tr, va) in enumerate(folds):
        train, valid = dev[tr], dev[va]
        membership.loc[valid, 'validation_fold'] = fold
        p, _ = fit(bulk, meta, cells, train, valid, args.sizes, args, args.output / f'fold_{fold}')
        for size in args.sizes:
            part = meta.iloc[valid].copy()
            part['size'], part['score'], part['fold'] = size, p[size].to_numpy(), fold
            tables.append(part)
    membership.to_csv(args.output / 'fold_membership.csv', index=False)
    predictions = pd.concat(tables, ignore_index=True)
    predictions.to_csv(args.output / 'validation_predictions.csv', index=False)
    metrics, selection = select_one_se(predictions, n_bootstrap=args.bootstrap, seed=args.seed)
    figures = plot_validation(metrics, selection, output_dir=args.output,
                              title=f'DEGAS {args.split}-level cross-validation')
    import matplotlib.pyplot as plt
    for figure in figures:
        plt.close(figure)
    selected = selection.loc[selection.selected, 'size'].astype(int).tolist()
    # Final models only use development patients; final holdout labels never tune sizes.
    p, c = fit(bulk, meta, cells, dev, test, selected, args, args.output / 'final', explain_final=True)
    p.to_csv(args.output / 'holdout_raw_scores_by_size.csv')
    c.to_csv(args.output / 'cell_raw_scores_by_size.csv')
    average_rank_scores(c, selected).to_csv(args.output / 'cell_average_rank_scores.csv')
    ranked = average_rank_scores(p, selected)
    ranked.to_csv(args.output / 'holdout_average_rank_scores.csv')
    # A percentile ensemble is not a probability: only discrimination is reported for it.
    from sklearn.metrics import roc_auc_score, average_precision_score
    y = meta.iloc[test].label
    report = {'rank_ensemble': {'AUROC': roc_auc_score(y, ranked),
                                'average_precision': average_precision_score(y, ranked)},
              'raw_probability_by_size': {size: classification_metrics(y, p[size]) for size in selected}}
    (args.output / 'holdout_metrics.json').write_text(json.dumps(report, indent=2))
    print(f'Completed: {args.output}; eligible sizes: {selected}')


if __name__ == '__main__':
    main()
