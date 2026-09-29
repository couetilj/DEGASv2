"""Frozen cohort references and resumable, all-cell expected-gradient explanations."""
import json
import os
from pathlib import Path
import uuid

import numpy as np
import pandas as pd


def choose_reference(scores, cell_ids, *, mode='low-risk-quartile', custom_ids=None,
                     background_size=128, max_cells=0, seed=42):
    """Use <= the empirical 25th percentile, retaining all boundary ties."""
    scores = np.asarray(scores, dtype=float)
    ids = list(cell_ids)
    if scores.shape != (len(ids),) or not len(ids) or not np.isfinite(scores).all():
        raise ValueError('Need one finite ensemble score per cell')
    if len(set(ids)) != len(ids) or background_size < 0 or max_cells < 0:
        raise ValueError('IDs must be unique and sampling caps nonnegative')
    cutoff = None
    if mode == 'low-risk-quartile':
        cutoff = float(np.quantile(scores, .25))
        eligible = np.flatnonzero(scores <= cutoff)
    elif mode == 'population':
        eligible = np.arange(len(ids))
    elif mode == 'custom':
        if not custom_ids or len(set(custom_ids)) != len(custom_ids) or set(custom_ids) - set(ids):
            raise ValueError('Invalid custom background IDs')
        wanted = set(custom_ids)
        eligible = np.array([i for i, v in enumerate(ids) if v in wanted])
    else:
        raise ValueError('Unknown background mode')
    rng = np.random.default_rng(seed)
    background = eligible
    if background_size and len(background) > background_size:
        background = np.sort(rng.choice(background, background_size, replace=False))
    observations = np.arange(len(ids))
    if max_cells and len(ids) > max_cells:
        observations = np.sort(np.random.default_rng(seed).choice(observations, max_cells, replace=False))
    return dict(mode=mode, cutoff=cutoff, tie_policy='include scores <= cutoff',
                score_definition='final average of within-size ranks after averaging seed probabilities',
                seed=seed, background_size_cap=background_size, max_cells=max_cells,
                eligible_indices=eligible.tolist(), background_indices=background.tolist(),
                observation_indices=observations.tolist(),
                eligible_ids=[ids[i] for i in eligible], background_ids=[ids[i] for i in background],
                observation_ids=[ids[i] for i in observations])


def reference_manifest(run, scores):
    from .sweep import _config
    settings = _config(run)['options']
    custom = run/'inputs/shap_background_ids.json'
    return choose_reference(scores, json.loads((run/'inputs/cell_ids.json').read_text()),
                            mode=settings['shap_background'],
                            custom_ids=json.loads(custom.read_text()) if custom.exists() else None,
                            background_size=settings['shap_background_size'],
                            max_cells=settings['shap_max_cells'], seed=settings['seed'])


def explain_worker(run, task):
    """One saved final model per job; each cell chunk is independently restartable."""
    import torch
    from .sweep import _config, _report, _completion, _fingerprint, _publish, digest, write_json
    from .models import load_models
    from .direct_shap import DiseaseScore, explain
    from .preprocess import preprocess_counts
    config = _config(run)
    settings = config['options']
    for file, sha in {**config['input_hashes'], **config['feature_hashes']}.items():
        if digest(run/file) != sha:
            raise ValueError(f'Frozen input changed: {file}')
    reference = json.loads((_report(run, 'final')/'shap_reference.json').read_text())
    trained = _completion(run, 'final', task)
    if trained is None:
        raise ValueError('Final model is incomplete')
    parent = run/'explain'/f"task_{task['task_id']:06d}"
    # Stable chunks are retained on retry; final completion is published only after all chunks.
    attempt = parent/'chunks'
    attempt.mkdir(parents=True, exist_ok=True)
    checkpoint = next((trained/'model').rglob('configs.json')).parent
    options = json.loads((checkpoint/'configs.json').read_text())
    options['save_dir'] = str(attempt/('restore_' + uuid.uuid4().hex))
    torch.set_num_threads(int(os.environ.get('SLURM_CPUS_PER_TASK', '1')))
    model = load_models(options)
    for name in model.net_name_list:
        getattr(model, name).load_state_dict(torch.load(
            checkpoint/f"{settings['iters']}_net_{name}.pth", map_location=model.device, weights_only=True))
    score = DiseaseScore(model)
    genes = json.loads((trained/'genes.json').read_text())
    all_genes = json.loads((run/'inputs/genes.json').read_text())
    lookup = {g: i for i, g in enumerate(all_genes)}
    columns = [lookup[g] for g in genes]
    cells = np.load(run/'inputs/cells.npy', mmap_mode='r')
    bg_indices = reference['background_indices']
    background = preprocess_counts(cells[bg_indices][:, columns]).astype('float32')
    observations = reference['observation_indices']
    predictions = np.load(trained/'cell_scores.npy', mmap_mode='r')
    fingerprint = _fingerprint(run, 'explain')
    artifacts = []
    for start in range(0, len(observations), settings['shap_chunk_size']):
        filename = f'chunk_{start:09d}.npz'
        marker = attempt/(filename + '.json')
        if marker.exists():
            info = json.loads(marker.read_text())
            if info['fingerprint'] != fingerprint or digest(attempt/info['file']) != info['sha256']:
                raise ValueError('Corrupt or stale SHAP chunk')
        else:
            indices = observations[start:start+settings['shap_chunk_size']]
            inputs = preprocess_counts(cells[indices][:, columns]).astype('float32')
            result = explain(score, background, inputs, nsamples=settings['shap_samples'], seed=settings['seed']+start)
            np.testing.assert_allclose(result['prediction'], predictions[indices], atol=1e-6, rtol=1e-5)
            # A unique immutable artifact plus first-wins marker handles duplicate workers.
            actual = f'{filename}.{uuid.uuid4().hex}.npz'
            np.savez_compressed(attempt/actual, **result, inputs=inputs, genes=np.asarray(genes),
                                observation_ids=np.asarray(reference['observation_ids'][start:start+len(indices)]),
                                background_ids=np.asarray(reference['background_ids']))
            _publish(marker, dict(fingerprint=fingerprint, file=actual, sha256=digest(attempt/actual)))
            info = json.loads(marker.read_text())
        artifacts.extend([marker.name, info['file']])
    write_json(attempt/'genes.json', genes)
    artifacts.append('genes.json')
    _publish(parent/'complete.json', dict(task=task, fingerprint=fingerprint, attempt=attempt.name,
                                         hashes={f: digest(attempt/f) for f in artifacts}))
    print(f"Explained {len(observations)} cells: size={task['size']}, seed={task['seed']}")


def collect_explanations(run, tasks):
    """Average seeds chunkwise into memory-mapped arrays, keeping gene sizes separate."""
    from .sweep import _config, _report, _completion, _fingerprint, _publish, digest, write_json
    settings = _config(run)['options']
    reference = json.loads((_report(run, 'final')/'shap_reference.json').read_text())
    report = run/'reports'/('explain_' + uuid.uuid4().hex)
    report.mkdir(parents=True)
    write_json(report/'shap_reference.json', reference)
    attempts = {t['task_id']: _completion(run, 'explain', t) for t in tasks}
    n = len(reference['observation_indices'])
    for size in sorted({t['size'] for t in tasks}):
        group = [t for t in tasks if t['size'] == size]
        folder = report/f'size_{size}'
        folder.mkdir()
        genes = json.loads((run/f'features/final_{size}.json').read_text())
        values = np.lib.format.open_memmap(folder/'values.npy', mode='w+', dtype='float64', shape=(n, size))
        inputs = np.lib.format.open_memmap(folder/'inputs.npy', mode='w+', dtype='float32', shape=(n, size))
        prediction, residual = np.zeros(n), np.zeros(n)
        bases, diagnostics = {}, []
        for start in range(0, n, settings['shap_chunk_size']):
            stop = min(n, start+settings['shap_chunk_size'])
            values[start:stop] = 0
            for task in group:
                attempt = attempts[task['task_id']]
                info = json.loads((attempt/f'chunk_{start:09d}.npz.json').read_text())
                with np.load(attempt/info['file']) as data:
                    if (data['observation_ids'].tolist() != reference['observation_ids'][start:stop]
                            or data['background_ids'].tolist() != reference['background_ids']
                            or data['genes'].tolist() != genes):
                        raise ValueError('SHAP chunk IDs/features do not match frozen reference')
                    if data['values'].shape != (stop-start, size) or not all(np.isfinite(data[k]).all() for k in ('values', 'inputs', 'prediction', 'base_value', 'residual')):
                        raise ValueError('Invalid SHAP chunk shape/values')
                    base = float(data['base_value'])
                    if task['seed'] in bases and not np.isclose(bases[task['seed']], base, atol=1e-7):
                        raise ValueError('Background base changed between chunks')
                    bases[task['seed']] = base
                    if data['inputs'].shape != (stop-start, size) or data['prediction'].shape != (stop-start,) or data['residual'].shape != (stop-start,):
                        raise ValueError('Invalid SHAP chunk input/prediction/residual shape')
                    np.testing.assert_allclose(data['residual'], data['prediction']-base-data['values'].sum(axis=1))
                    values[start:stop] += data['values']/len(group)
                    inputs[start:stop] = data['inputs']
                    prediction[start:stop] += data['prediction']/len(group)
                    diagnostics.append(dict(seed=task['seed'], start=start, stop=stop,
                                            mean_abs_residual=float(np.abs(data['residual']).mean()),
                                            p95_abs_residual=float(np.quantile(np.abs(data['residual']), .95))))
            residual[start:stop] = prediction[start:stop]-np.mean(list(bases.values()))-values[start:stop].sum(axis=1)
        values.flush()
        inputs.flush()
        np.save(folder/'prediction.npy', prediction)
        np.save(folder/'residual.npy', residual)
        write_json(folder/'metadata.json', dict(genes=genes, observation_ids=reference['observation_ids'],
                                               background_ids=reference['background_ids'],
                                               base_value=float(np.mean(list(bases.values()))),
                                               explained_output='mean raw class-1 probability across seeds'))
        pd.DataFrame(diagnostics).to_csv(folder/'diagnostics.csv', index=False)
        del values, inputs
    _publish(run/'explain_pooled.json', dict(directory=str(report.relative_to(run)),
                                            fingerprint=_fingerprint(run, 'explain'),
                                            hashes={str(p.relative_to(report)): digest(p) for p in report.rglob('*') if p.is_file()}))
    return _report(run, 'explain')
