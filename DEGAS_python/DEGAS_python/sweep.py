"""Run, distribute and pool DEGAS validation: python -m DEGAS_python.sweep --help."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import uuid

import numpy as np
import pandas as pd
from .validation import (holdout_split, validation_folds, select_one_se,
                         classification_metrics, average_rank_scores)
from .workflow_inputs import demo, read_inputs

GRID = (10, 20, 40, 60, 120, 250, 500, 750, 1000, 1500, 2500, 5000, 7500, 10000)


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def source_hashes():
    root = Path(__file__).parent
    return {str(p.relative_to(root)): digest(p) for p in sorted(root.rglob('*.py'))}


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name('.' + path.name + '.' + uuid.uuid4().hex + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    os.replace(temporary, path)


def empirical_sizes(n_genes, max_genes=None):
    """Broad default grid plus the available/configured maximum; no small-gene cap."""
    limit = n_genes if max_genes is None else min(n_genes, max_genes)
    if limit < 2:
        raise ValueError('Need at least two shared genes')
    return sorted({n for n in GRID if n <= limit} | {limit})


def feature_candidates(common_genes, *, mode='user', sizes=None, genes=None, max_genes=None):
    """Validate exact gene lists or HVG counts; never silently truncate a request."""
    if mode not in ('user', 'empirical'):
        raise ValueError('feature_mode must be user or empirical')
    if genes is not None:
        if mode != 'user' or sizes is not None or max_genes is not None:
            raise ValueError('--genes is an exact list: use user mode without sizes/max-genes')
        names = [line.strip() for line in Path(genes).read_text().splitlines() if line.strip()]
        if len(names) < 2 or len(set(names)) != len(names):
            raise ValueError('Gene list must contain at least two unique IDs, one per line')
        missing = set(names) - set(common_genes)
        if missing:
            raise ValueError(f'Requested genes absent from bulk or cells: {sorted(missing)[:20]}')
        return [len(names)], names
    if max_genes is not None and mode != 'empirical':
        raise ValueError('--max-genes applies only to empirical mode')
    if sizes is None:
        sizes = empirical_sizes(len(common_genes), max_genes) if mode == 'empirical' else [10, 20, 40]
    if (not sizes or len(set(sizes)) != len(sizes)
            or any(not isinstance(n, int) or n < 2 or n > len(common_genes) for n in sizes)):
        raise ValueError('Every requested size must be a unique integer from 2 to shared-gene count')
    if max_genes is not None and max(sizes) > max_genes:
        raise ValueError('Explicit sizes exceed --max-genes')
    return sorted(sizes), None


def _tasks(folds, sizes, seeds):
    return [dict(task_id=i, fold=f, size=n, seed=s)
            for i, (f, n, s) in enumerate((f, n, s) for f in folds for n in sizes for s in range(seeds))]


def prepare(args):
    """Freeze inputs, splits, training-only feature lists and task manifests once."""
    if min(args.seeds, args.iters, args.shap_samples, args.batch_size, args.feature_dim) < 1:
        raise ValueError('Training counts and dimensions must be positive')
    if args.batch_size < 2 or args.feature_dim < 2 or args.bootstrap < 2:
        raise ValueError('Batch size, feature dimension and bootstrap count must be >=2')
    if args.shap_background_size < 0 or args.shap_max_cells < 0 or args.shap_chunk_size < 1:
        raise ValueError('SHAP sizes must be nonnegative; chunk size must be positive')
    if (args.shap_background == 'custom') != (args.shap_background_ids is not None):
        raise ValueError('Custom background requires --shap-background custom and --shap-background-ids')
    bulk, meta, cells = read_inputs(args.input) if args.input else demo(args.demo_genes)
    if not bulk.columns.is_unique or not cells.columns.is_unique:
        raise ValueError('Gene IDs must be unique')
    common = bulk.columns.intersection(cells.columns, sort=False)
    bulk, cells = bulk[common], cells[common]
    for matrix in (bulk, cells):
        if not len(matrix) or not np.isfinite(matrix.to_numpy()).all() or (matrix.to_numpy() < 0).any():
            raise ValueError('Supply finite, nonnegative raw counts')
    custom_ids = None
    if args.shap_background_ids is not None:
        custom_ids = [v.strip() for v in args.shap_background_ids.read_text().splitlines() if v.strip()]
        if not custom_ids or len(set(custom_ids)) != len(custom_ids) or set(custom_ids) - set(cells.index.astype(str)):
            raise ValueError('Background IDs must be nonempty, unique and present in cells')
    sizes, exact = feature_candidates(common, mode=args.feature_mode, sizes=args.sizes,
                                      genes=args.genes, max_genes=args.max_genes)
    dev, test = holdout_split(meta, mode=args.split, fraction=args.holdout_fraction, seed=args.seed)
    development = meta.iloc[dev].reset_index(drop=True)
    if (pd.crosstab(development.study, development.label).reindex(columns=[0, 1], fill_value=0) < 2).any().any():
        raise ValueError('Development data need >=2 patients per study/class for bootstrap SE')
    folds = validation_folds(development, mode=args.split, n_splits=args.folds, seed=args.seed)
    run = args.output.resolve()
    run.mkdir(parents=True, exist_ok=False)
    (run/'inputs').mkdir()
    (run/'features').mkdir()
    np.save(run/'inputs/bulk.npy', bulk.to_numpy(dtype='float32'))
    np.save(run/'inputs/cells.npy', cells.to_numpy(dtype='float32'))
    meta.to_csv(run/'inputs/patients.csv', index=False)
    write_json(run/'inputs/genes.json', list(common))
    write_json(run/'inputs/cell_ids.json', list(cells.index.astype(str)))
    if custom_ids is not None:
        write_json(run/'inputs/shap_background_ids.json', custom_ids)
    split_map = {str(i): dict(train=dev[tr].tolist(), evaluate=dev[va].tolist())
                 for i, (tr, va) in enumerate(folds)}
    split_map['final'] = dict(train=dev.tolist(), evaluate=test.tolist())
    membership = meta.copy()
    membership['partition'] = 'development'
    membership.loc[test, 'partition'] = 'final_holdout'
    membership['validation_fold'] = -1
    for i, (_, va) in enumerate(folds):
        membership.loc[dev[va], 'validation_fold'] = i
    membership.to_csv(run/'fold_membership.csv', index=False)
    feature_hashes = {}
    for fold, split in split_map.items():
        if exact is None:
            variance = np.log1p(bulk.iloc[split['train']].to_numpy()).var(axis=0)
            ranked = common[np.argsort(-variance, kind='stable')].tolist()
        else:
            ranked = exact
        for size in sizes:
            relative = f'features/{fold}_{size}.json'
            write_json(run/relative, ranked[:size])
            feature_hashes[relative] = digest(run/relative)
    options = {k: str(v.resolve()) if isinstance(v, Path) else v for k, v in vars(args).items()}
    versions = {p: importlib.metadata.version(p) for p in ['numpy', 'pandas', 'scipy', 'scikit-learn', 'torch']}
    config = dict(schema=1, run_id=uuid.uuid4().hex, options=options, sizes=sizes,
                  selector='exact-list' if exact is not None else 'training-log1p-variance', source_hashes=source_hashes(),
                  splits=split_map, versions=versions, feature_hashes=feature_hashes,
                  input_hashes={str(p.relative_to(run)): digest(p) for p in (run/'inputs').iterdir()})
    write_json(run/'run.json', config)
    write_json(run/'validation_tasks.json', _tasks(range(len(folds)), sizes, args.seeds))
    print(f'Prepared {len(folds)*len(sizes)*args.seeds} validation tasks; sizes={sizes}; run={run}')
    return run


def _config(run):
    config = json.loads((run/'run.json').read_text())
    if config.get('source_hashes') != source_hashes():
        raise ValueError('DEGAS source changed since preparation; restore the recorded source or prepare a new run')
    return config


def _task_list(run, phase):
    config = _config(run)
    if phase == 'explain':
        _report(run, 'final')
        if not config['options']['shap']:
            raise ValueError('Prepare with --shap to enable explanations')
        expected = _task_list(run, 'final')
    elif phase == 'validation':
        expected = _tasks(range(len(config['splits'])-1), config['sizes'], config['options']['seeds'])
    else:
        report = _report(run, 'validation')
        selection = pd.read_csv(report/'feature_size_selection.csv')
        sizes = selection.loc[selection.selected, 'size'].astype(int).tolist()
        expected = _tasks(['final'], sizes, config['options']['seeds'])
    saved = json.loads((run/f'{phase}_tasks.json').read_text())
    if saved != expected:
        raise ValueError(f'{phase} task manifest does not match frozen configuration/selection')
    return saved


def _fingerprint(run, phase):
    value = digest(run/'run.json') + ':' + digest(run/f'{phase}_tasks.json')
    if phase == 'explain':
        value += ':' + digest(run/'final_pooled.json')
    return value


def _completion(run, phase, task):
    directory = run/phase/f"task_{task['task_id']:06d}"
    marker = directory/'complete.json'
    if not marker.exists():
        return None
    info = json.loads(marker.read_text())
    if info['task'] != task or info['fingerprint'] != _fingerprint(run, phase):
        raise ValueError(f'Stale or mismatched task: {marker}')
    attempt = directory/info['attempt']
    for file, sha in info['hashes'].items():
        if digest(attempt/file) != sha:
            raise ValueError(f'Corrupt task artifact: {attempt/file}')
    return attempt


def _publish(marker, info):
    """First successful attempt wins; a partial write never marks a task complete."""
    temp = marker.with_name(uuid.uuid4().hex + '.json')
    write_json(temp, info)
    try:
        os.link(temp, marker)
    except FileExistsError:
        pass
    finally:
        temp.unlink()


def _preprocess(raw, columns, path, batch=1024):
    from .preprocess import preprocess_counts
    output = np.lib.format.open_memmap(path, mode='w+', dtype='float32', shape=(len(raw), len(columns)))
    for start in range(0, len(raw), batch):
        output[start:start+batch] = preprocess_counts(raw[start:start+batch][:, columns]).astype('float32')
    output.flush()
    return output


def worker(run, phase, task_id):
    """One independent fold/size/seed fit, restartable without overwriting attempts."""
    run = Path(run).resolve()
    tasks = _task_list(run, phase)
    if task_id < 0 or task_id >= len(tasks):
        raise ValueError('Task ID is outside the manifest')
    task = tasks[task_id]
    if _completion(run, phase, task) is not None:
        print(f'{phase} task {task_id} already complete')
        return
    if phase == 'explain':
        from .shap_workflow import explain_worker
        return explain_worker(run, task)
    config = _config(run)
    settings = config['options']
    feature_file = f"features/{task['fold']}_{task['size']}.json"
    if digest(run/feature_file) != config['feature_hashes'][feature_file]:
        raise ValueError('Frozen feature list changed')
    selected = json.loads((run/feature_file).read_text())
    genes = json.loads((run/'inputs/genes.json').read_text())
    gene_index = {gene: i for i, gene in enumerate(genes)}
    columns = [gene_index[gene] for gene in selected]
    split = config['splits'][str(task['fold'])]
    train, evaluate = np.asarray(split['train']), np.asarray(split['evaluate'])
    meta = pd.read_csv(run/'inputs/patients.csv', dtype={'patient_id': str, 'study': str})
    parent = run/phase/f'task_{task_id:06d}'
    attempt = parent/('attempt_' + uuid.uuid4().hex)
    attempt.mkdir(parents=True)
    import torch
    from .options import BlankClass_opt
    from .run_module import run_model
    from .direct_shap import DiseaseScore, explain
    torch.set_num_threads(int(os.environ.get('SLURM_CPUS_PER_TASK', '1')))
    bulk = np.load(run/'inputs/bulk.npy', mmap_mode='r')
    cells = np.load(run/'inputs/cells.npy', mmap_mode='r')
    x = _preprocess(bulk, columns, attempt/'bulk_preprocessed.npy')
    z = _preprocess(cells, columns, attempt/'cells_preprocessed.npy')
    opt = dict(BlankClass_opt, save_dir=str(attempt/'model'), seed=task['seed'],
               fold=-1, tot_folds=1, tot_iters=settings['iters'], save_freq=settings['iters'],
               is_save=True, feature_dim=settings['feature_dim'], batch_size=settings['batch_size'],
               pat_batch_size=settings['batch_size'])
    path, model = run_model(opt, x[train], meta.label.to_numpy()[train], z,
                            pat_eval_expr_mat=x[evaluate], pat_eval_lab_mat=meta.label.to_numpy()[evaluate],
                            sc_eval_expr_mat=z[:min(32, len(z))], return_model=True)
    score = DiseaseScore(model)
    device = next(score.parameters()).device
    with torch.no_grad():
        p = score(torch.as_tensor(x[evaluate].copy(), device=device)).cpu().numpy().ravel()
    native = pd.read_csv(Path(path)/f"low_reso_results_epoch_{settings['iters']}.csv").sort_values('pid')
    np.testing.assert_allclose(p, native.hazard, atol=1e-6, rtol=1e-5)
    result = meta.iloc[evaluate].copy()
    result['score'] = p
    result.to_csv(attempt/'predictions.csv', index=False)
    artifacts = ['predictions.csv', 'genes.json']
    write_json(attempt/'genes.json', selected)
    if phase == 'final':
        c = np.empty(len(z), dtype='float32')
        with torch.no_grad():
            for start in range(0, len(z), 1024):
                c[start:start+1024] = score(torch.as_tensor(np.array(z[start:start+1024]), device=device)).cpu().numpy().ravel()
        np.save(attempt/'cell_scores.npy', c)
        artifacts.append('cell_scores.npy')
        artifacts.extend(str(p.relative_to(attempt)) for p in (attempt/'model').rglob('*')
                         if p.is_file() and (p.suffix == '.pth' or p.name == 'configs.json'))
    # Large normalized scratch matrices can be reconstructed from frozen raw inputs.
    del x, z
    (attempt/'bulk_preprocessed.npy').unlink()
    (attempt/'cells_preprocessed.npy').unlink()
    _publish(parent/'complete.json', dict(task=task, fingerprint=_fingerprint(run, phase),
                                         attempt=attempt.name, hashes={f: digest(attempt/f) for f in artifacts},
                                         slurm_job_id=os.environ.get('SLURM_JOB_ID'), device=str(device),
                                         versions={p: importlib.metadata.version(p) for p in config['versions']}))
    print(f'Completed {phase} task {task_id}: size={task["size"]}, seed={task["seed"]}')


def _report(run, phase):
    marker = run/f'{phase}_pooled.json'
    info = json.loads(marker.read_text())
    if info['fingerprint'] != _fingerprint(run, phase):
        raise ValueError('Pooled report does not match current run/task manifest')
    report = run/info['directory']
    for name, sha in info['hashes'].items():
        if digest(report/name) != sha:
            raise ValueError('Pooled report was modified')
    link = run/'reports'/phase
    try:
        link.symlink_to(report.name, target_is_directory=True)
    except FileExistsError:
        if link.resolve() != report.resolve():
            raise ValueError('Report link points to another result')
    return report


def collect(run, phase):
    """Pool complete seed fits only; validate exact coverage before any selection."""
    run = Path(run).resolve()
    if (run/f'{phase}_pooled.json').exists():
        return _report(run, phase)
    config = _config(run)
    tasks = _task_list(run, phase)
    missing = [t['task_id'] for t in tasks if _completion(run, phase, t) is None]
    if missing:
        raise ValueError(f'Missing {phase} tasks: {missing}; resubmit these IDs; nothing was pooled')
    for file, sha in {**config['input_hashes'], **config['feature_hashes']}.items():
        if digest(run/file) != sha:
            raise ValueError(f'Frozen input changed: {file}')
    if phase == 'explain':
        from .shap_workflow import collect_explanations
        return collect_explanations(run, tasks)
    metadata = pd.read_csv(run/'inputs/patients.csv', dtype={'patient_id': str, 'study': str})
    tables, cell_sums = [], {}
    for task in tasks:
        attempt = _completion(run, phase, task)
        part = pd.read_csv(attempt/'predictions.csv', dtype={'patient_id': str, 'study': str})
        expected = metadata.iloc[config['splits'][str(task['fold'])]['evaluate']].reset_index(drop=True)
        if list(part.patient_id) != list(expected.patient_id) or not part[['patient_id', 'study', 'label']].equals(expected[['patient_id', 'study', 'label']]):
            raise ValueError(f'Task {task["task_id"]} predictions have incorrect IDs/labels/order')
        classification_metrics(part.label, part.score)
        if json.loads((attempt/'genes.json').read_text()) != json.loads((run/f"features/{task['fold']}_{task['size']}.json").read_text()):
            raise ValueError('Task used incorrect features')
        part['size'], part['seed'], part['fold'] = task['size'], task['seed'], task['fold']
        tables.append(part)
        if phase == 'final':
            c = np.load(attempt/'cell_scores.npy')
            n_cells = len(json.loads((run/'inputs/cell_ids.json').read_text()))
            if c.shape != (n_cells,) or not np.isfinite(c).all() or ((c < 0) | (c > 1)).any():
                raise ValueError('Invalid cell prediction coverage/values')
            cell_sums.setdefault(task['size'], np.zeros(n_cells))
            cell_sums[task['size']] += c / config['options']['seeds']
    submodels = pd.concat(tables, ignore_index=True)
    if submodels.duplicated(['patient_id', 'size', 'seed']).any():
        raise ValueError('Duplicate patient/size/seed predictions')
    counts = submodels.groupby(['patient_id', 'size']).size()
    if not (counts == config['options']['seeds']).all():
        raise ValueError('Incomplete seed coverage')
    predictions = submodels.groupby(['patient_id', 'study', 'label', 'size'], as_index=False).score.mean()
    report = run/'reports'/f'{phase}_{uuid.uuid4().hex}'
    report.mkdir(parents=True)
    submodels.to_csv(report/'submodel_predictions.csv', index=False)
    predictions.to_csv(report/'validation_predictions.csv' if phase == 'validation' else report/'holdout_predictions.csv', index=False)
    if phase == 'validation':
        from .evaluation_plots import plot_validation
        import matplotlib.pyplot as plt
        metrics, selection = select_one_se(predictions, n_bootstrap=config['options']['bootstrap'], seed=config['options']['seed'])
        figures = plot_validation(metrics, selection, title=f"DEGAS {config['options']['split']}-level validation", output_dir=report)
        for figure in figures:
            plt.close(figure)
        sizes = selection.loc[selection.selected, 'size'].astype(int).tolist()
        # Exact gene lists have a single candidate, whose size is retained trivially.
        write_json(run/'final_tasks.json', _tasks(['final'], sizes, config['options']['seeds']))
    else:
        _final_report(run, report, predictions, cell_sums)
    _publish(run/f'{phase}_pooled.json', dict(directory=str(report.relative_to(run)),
                                            fingerprint=_fingerprint(run, phase),
                                            hashes={str(p.relative_to(report)): digest(p) for p in report.rglob('*') if p.is_file()}))
    report = _report(run, phase)
    print(f'Pooled {phase} results: {run/"reports"/phase}')
    return report


def _final_report(run, report, predictions, cell_sums):
    from sklearn.metrics import roc_auc_score, average_precision_score
    p = predictions.pivot(index=['patient_id', 'label'], columns='size', values='score')
    y = p.index.get_level_values('label').to_numpy()
    p.index = p.index.get_level_values('patient_id')
    c = pd.DataFrame(cell_sums, index=json.loads((run/'inputs/cell_ids.json').read_text()))
    p.to_csv(report/'holdout_raw_scores_by_size.csv')
    c.to_csv(report/'cell_raw_scores_by_size.csv')
    cell_ranks = average_rank_scores(c, list(c.columns))
    cell_ranks.to_csv(report/'cell_average_rank_scores.csv')
    if _config(run)['options']['shap']:
        from .shap_workflow import reference_manifest
        reference = reference_manifest(run, cell_ranks.to_numpy())
        write_json(report/'shap_reference.json', reference)
        write_json(run/'explain_tasks.json', _task_list(run, 'final'))
    ranked = average_rank_scores(p, list(p.columns))
    ranked.to_csv(report/'holdout_average_rank_scores.csv')
    write_json(report/'holdout_metrics.json', dict(
        aggregation='pooled final holdout patients', n_patients=len(y),
        rank_ensemble=dict(AUROC=roc_auc_score(y, ranked), average_precision=average_precision_score(y, ranked)),
        raw_probability_by_size={int(size): classification_metrics(y, p[size]) for size in p.columns}))


def submit(args):
    """Submit an array and an afterok CPU pooling job; dry-run prints exact argv."""
    run = args.run.resolve()
    tasks = _task_list(run, args.phase)
    pending = [t['task_id'] for t in tasks if _completion(run, args.phase, t) is None]
    if not pending:
        print('All tasks complete; run collect for this phase if needed.')
        return
    if args.parallel < 1 or args.cpus < 1 or args.gpus < 0:
        raise ValueError('parallel/cpus must be positive; gpus cannot be negative')
    scripts = run/'slurm'
    scripts.mkdir(exist_ok=True)
    (run/'logs').mkdir(exist_ok=True)
    worker_script = scripts/'worker.sh'
    pool_script = scripts/'pool.sh'
    worker_script.write_text('#!/usr/bin/env bash\nset -euo pipefail\n'
                             'export MPLBACKEND=Agg\n'
                             'export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"\n'
                             'export MKL_NUM_THREADS="$OMP_NUM_THREADS" OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"\n'
                             'exec "$1" -m DEGAS_python.sweep worker --run "$2" --phase "$3" --task-id "${SLURM_ARRAY_TASK_ID:?}"\n')
    pool_script.write_text('#!/usr/bin/env bash\nset -euo pipefail\n'
                           'export MPLBACKEND=Agg\n'
                             'export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"\n'
                             'export MKL_NUM_THREADS="$OMP_NUM_THREADS" OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"\n'
                           'exec "$1" -m DEGAS_python.sweep collect --run "$2" --phase "$3"\n')
    base = ['sbatch', '--parsable', '--chdir', str(run)]
    if args.account:
        base += ['--account', args.account]
    train = base + ['--job-name', f'degas-{args.phase}', '--array',
                    ','.join(map(str, pending)) + f'%{args.parallel}',
                    '--cpus-per-task', str(args.cpus), '--mem', args.mem, '--time', args.time,
                    '--output', str(run/'logs'/f'{args.phase}_%A_%a.log')]
    if args.partition:
        train += ['--partition', args.partition]
    if args.gpus:
        train += ['--gpus', str(args.gpus)]
    train += [str(worker_script), str(Path(args.python).resolve()), str(run), args.phase]
    pool = base + ['--job-name', f'degas-pool-{args.phase}', '--kill-on-invalid-dep=yes', '--cpus-per-task', '1',
                   '--mem', args.pool_mem, '--time', args.pool_time,
                   '--output', str(run/'logs'/f'{args.phase}_pool_%j.log')]
    if args.pool_partition:
        pool += ['--partition', args.pool_partition]
    tail = [str(pool_script), str(Path(args.python).resolve()), str(run), args.phase]
    if args.dry_run:
        print(shlex.join(train))
        print(shlex.join(pool + ['--dependency', 'afterok:ARRAY_JOB_ID'] + tail))
        return
    # Record a successful array submission even if submission of its pool job fails.
    job = subprocess.run(train, check=True, text=True, capture_output=True).stdout.strip().split(';')[0]
    if not job.isdigit():
        raise RuntimeError(f'Unexpected sbatch job ID: {job!r}')
    record = scripts/f'submission_{args.phase}_{uuid.uuid4().hex}.json'
    commands = dict(array_job_id=job, array_command=train, pending_task_ids=pending)
    write_json(record, commands)
    pool_command = pool + ['--dependency', f'afterok:{job}'] + tail
    pooled_job = subprocess.run(pool_command, check=True, text=True, capture_output=True).stdout.strip().split(';')[0]
    commands.update(pool_job_id=pooled_job, pool_command=pool_command)
    write_json(record, commands)
    print(f'Submitted {args.phase} array {job} and pooling job {pooled_job}; record={record}')


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for action in ('run', 'prepare'):
        p = commands.add_parser(action, help='Prepare frozen inputs' + (' and execute locally' if action == 'run' else ' for parallel jobs'))
        p.add_argument('--input', type=Path, help='CSV directory; omit for synthetic data')
        p.add_argument('--output', type=Path, required=True)
        p.add_argument('--feature-mode', choices=['user', 'empirical'], default='user')
        feature = p.add_mutually_exclusive_group()
        feature.add_argument('--sizes', type=int, nargs='+', help='Requested HVG counts or override empirical grid')
        feature.add_argument('--genes', type=Path, help='Exact gene list, one gene per line without header')
        p.add_argument('--max-genes', type=int, help='Optional empirical-search budget cap')
        p.add_argument('--split', choices=['study', 'patient'], default='study')
        p.add_argument('--holdout-fraction', type=float, default=.1)
        p.add_argument('--folds', type=int, default=5)
        p.add_argument('--seeds', type=int, default=5)
        p.add_argument('--seed', type=int, default=42)
        p.add_argument('--iters', type=int, default=300)
        p.add_argument('--batch-size', type=int, default=32)
        p.add_argument('--feature-dim', type=int, default=16)
        p.add_argument('--bootstrap', type=int, default=20000)
        p.add_argument('--shap', action='store_true')
        p.add_argument('--shap-samples', type=int, default=512)
        p.add_argument('--shap-background', choices=['low-risk-quartile', 'population', 'custom'], default='low-risk-quartile')
        p.add_argument('--shap-background-ids', type=Path, help='Custom background: one cell ID per line')
        p.add_argument('--shap-background-size', type=int, default=0, help='Optional background sample cap; default 0 uses all eligible cells')
        p.add_argument('--shap-max-cells', type=int, default=0, help='Optional explained-cell sample cap; default 0 explains all')
        p.add_argument('--shap-chunk-size', type=int, default=256, help='Explained cells per resumable chunk')
        p.add_argument('--demo-genes', type=int, default=40)
    for action in ('worker', 'collect', 'status', 'submit'):
        p = commands.add_parser(action)
        p.add_argument('--run', type=Path, required=True)
        p.add_argument('--phase', choices=['validation', 'final', 'explain'], default='validation')
        if action == 'worker':
            p.add_argument('--task-id', type=int, required=True)
        if action == 'submit':
            p.add_argument('--account')
            p.add_argument('--partition', help='Training partition')
            p.add_argument('--pool-partition', help='CPU pooling partition (site default if omitted)')
            p.add_argument('--parallel', type=int, default=8)
            p.add_argument('--cpus', type=int, default=1)
            p.add_argument('--gpus', type=int, default=0)
            p.add_argument('--mem', default='16G', help='Memory per training task; size for your largest model')
            p.add_argument('--time', default='04:00:00')
            p.add_argument('--pool-mem', default='8G')
            p.add_argument('--pool-time', default='01:00:00')
            p.add_argument('--python', default=sys.executable, help='Shared compute-node Python environment')
            p.add_argument('--dry-run', action='store_true')
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command in ('run', 'prepare'):
        run = prepare(args)
        if args.command == 'run':
            for phase in (('validation', 'final', 'explain') if args.shap else ('validation', 'final')):
                for task in _task_list(run, phase):
                    worker(run, phase, task['task_id'])
                collect(run, phase)
    elif args.command == 'worker':
        worker(args.run, args.phase, args.task_id)
    elif args.command == 'collect':
        print(collect(args.run.resolve(), args.phase))
    elif args.command == 'submit':
        submit(args)
    else:
        tasks = _task_list(args.run, args.phase)
        missing = [t['task_id'] for t in tasks if _completion(args.run, args.phase, t) is None]
        print(json.dumps(dict(phase=args.phase, total=len(tasks), complete=len(tasks)-len(missing), missing=missing)))


if __name__ == '__main__':
    main()
