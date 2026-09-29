"""Features, frozen split integrity, complete pooling, restart and SLURM commands."""
from concurrent.futures import ThreadPoolExecutor
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from DEGAS_python import sweep


class SweepTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def prepare(self, name='run', extras=()):
        args = sweep.build_parser().parse_args(['prepare', '--output', str(self.root/name),
            '--sizes', '10', '20', '--seeds', '2', '--iters', '2', '--bootstrap', '30', *extras])
        return sweep.prepare(args)

    def fake_results(self, run):
        meta = pd.read_csv(run/'inputs/patients.csv', dtype={'patient_id': str, 'study': str})
        cfg = sweep._config(run)
        for task in sweep._task_list(run, 'validation'):
            parent = run/'validation'/f'task_{task["task_id"]:06d}'
            attempt = parent/'attempt_test'
            attempt.mkdir(parents=True)
            part = meta.iloc[cfg['splits'][str(task['fold'])]['evaluate']].copy()
            part['score'] = .2 + .6*part.label + task['seed']*.01
            part.to_csv(attempt/'predictions.csv', index=False)
            (attempt/'genes.json').write_bytes((run/f'features/{task["fold"]}_{task["size"]}.json').read_bytes())
            sweep.write_json(parent/'complete.json', dict(task=task, fingerprint=sweep._fingerprint(run, 'validation'),
                attempt=attempt.name, hashes={f: sweep.digest(attempt/f) for f in ['predictions.csv', 'genes.json']}))

    def test_broad_automatic_grid_and_strict_user_sizes(self):
        self.assertEqual(sweep.empirical_sizes(16561)[-4:], [5000, 7500, 10000, 16561])
        self.assertEqual(sweep.empirical_sizes(337)[-2:], [250, 337])
        self.assertEqual(sweep.empirical_sizes(16561, 3000)[-2:], [2500, 3000])
        self.assertEqual(sweep.empirical_sizes(3), [3])
        genes = ['a', 'b', 'c']
        for sizes in ([4], [1], [2, 2]):
            with self.assertRaises(ValueError):
                sweep.feature_candidates(genes, sizes=sizes)
        with self.assertRaises(ValueError):
            sweep.feature_candidates(genes, mode='empirical', max_genes=1)

    def test_exact_gene_order_and_missing_genes(self):
        file = self.root/'genes.txt'
        file.write_text('gene_9\ngene_1\ngene_4\n')
        sizes, genes = sweep.feature_candidates([f'gene_{i}' for i in range(10)], genes=file)
        self.assertEqual((sizes, genes), ([3], ['gene_9', 'gene_1', 'gene_4']))
        (self.root/'bulk_counts.csv').write_text('patient_id,gene_1,gene_1\np0,1,2\n')
        with self.assertRaisesRegex(ValueError, 'unique, nonempty gene IDs'):
            sweep.read_inputs(self.root)
        (self.root/'bulk_counts.csv').write_text('patient_id,gene_1,gene_2\n0001,1,2\n0002,3,4\n')
        (self.root/'patients.csv').write_text('patient_id,study,label\n0001,A,0\n0002,A,1\n')
        (self.root/'cell_counts.csv').write_text(',gene_1,gene_2\n0003,1,2\n')
        (self.root/'cells.csv').write_text('cell_id,patient_id\n0003,reference_donor\n')
        bulk, _, cells = sweep.read_inputs(self.root)
        self.assertEqual(bulk.index.tolist(), ['0001', '0002'])
        self.assertEqual(cells.index.tolist(), ['0003'])
        for content in ('absent\ngene_1\n', 'gene_1\ngene_1\n'):
            file.write_text(content)
            with self.assertRaises(ValueError):
                sweep.feature_candidates([f'gene_{i}' for i in range(10)], genes=file)
        file.write_text('gene_1\ngene_4\n')
        with self.assertRaises(ValueError):
            sweep.feature_candidates([f'gene_{i}' for i in range(10)], genes=file, sizes=[2])

    def test_features_and_inner_folds_never_use_final_test(self):
        run = self.prepare()
        cfg = sweep._config(run)
        final_test = set(cfg['splits']['final']['evaluate'])
        bulk = np.load(run/'inputs/bulk.npy')
        genes = json.loads((run/'inputs/genes.json').read_text())
        for fold, split in cfg['splits'].items():
            self.assertFalse(set(split['train']) & final_test)
            if fold != 'final':
                self.assertFalse(set(split['evaluate']) & final_test)
            # Deliberately destroy all non-training expression: selection is unchanged.
            changed = bulk.copy()
            changed[list(set(range(len(bulk)))-set(split['train']))] = 1e20
            ranks = np.argsort(-np.log1p(changed[split['train']].astype(float)).var(axis=0), kind='stable')
            expected = [genes[i] for i in ranks[:10]]
            self.assertEqual(json.loads((run/f'features/{fold}_10.json').read_text()), expected)

    def test_collection_requires_every_seed_and_is_idempotent(self):
        run = self.prepare()
        with self.assertRaisesRegex(ValueError, 'Missing validation tasks'):
            sweep.collect(run, 'validation')
        self.assertFalse((run/'validation_pooled.json').exists())
        self.fake_results(run)
        report = sweep.collect(run, 'validation')
        part = pd.read_csv(report/'validation_predictions.csv')
        np.testing.assert_allclose(part.score, .205 + .6*part.label)
        self.assertEqual(len(part), 90*2)
        self.assertEqual(len(sweep._task_list(run, 'final')), 4)
        self.assertEqual(sweep.collect(run, 'validation'), report)
        self.assertTrue((report/'auroc_selection.pdf').exists())
        self.assertTrue((report/'evaluation_metrics.png').exists())

    def test_corrupt_and_stale_results_are_rejected(self):
        run = self.prepare()
        self.fake_results(run)
        tasks = sweep._task_list(run, 'validation')
        attempt = sweep._completion(run, 'validation', tasks[0])
        (attempt/'predictions.csv').write_text('wrong')
        with self.assertRaisesRegex(ValueError, 'Corrupt task artifact'):
            sweep.collect(run, 'validation')
        # A changed configuration cannot silently reuse previously completed tasks.
        cfg = sweep._config(run)
        cfg['options']['iters'] += 1
        sweep.write_json(run/'run.json', cfg)
        with self.assertRaisesRegex(ValueError, 'Stale or mismatched'):
            sweep._completion(run, 'validation', tasks[1])
        cfg['source_hashes']['sweep.py'] = 'changed'
        sweep.write_json(run/'run.json', cfg)
        with self.assertRaisesRegex(ValueError, 'source changed'):
            sweep._task_list(run, 'validation')

    def test_duplicate_predictions_and_modified_inputs_are_rejected(self):
        run = self.prepare()
        self.fake_results(run)
        task = sweep._task_list(run, 'validation')[0]
        attempt = sweep._completion(run, 'validation', task)
        part = pd.read_csv(attempt/'predictions.csv')
        part.loc[1, 'patient_id'] = part.loc[0, 'patient_id']
        part.to_csv(attempt/'predictions.csv', index=False)
        marker = attempt.parent/'complete.json'
        info = json.loads(marker.read_text())
        info['hashes']['predictions.csv'] = sweep.digest(attempt/'predictions.csv')
        sweep.write_json(marker, info)
        with self.assertRaisesRegex(ValueError, 'incorrect IDs'):
            sweep.collect(run, 'validation')
        with (run/'inputs/bulk.npy').open('ab') as f:
            f.write(b'changed')
        with self.assertRaisesRegex(ValueError, 'Frozen input changed'):
            sweep.collect(run, 'validation')

    def test_slurm_dependency_resources_and_dry_run(self):
        run = self.prepare('with spaces')
        args = sweep.build_parser().parse_args(['submit', '--run', str(run), '--parallel', '3',
             '--gpus', '1', '--partition', 'gpu', '--pool-partition', 'cpu', '--account', 'my-account'])
        with patch('DEGAS_python.sweep.subprocess.run') as command:
            command.side_effect = [subprocess.CompletedProcess([], 0, '123;cluster\n'),
                                   subprocess.CompletedProcess([], 0, '124\n')]
            sweep.submit(args)
            train, pool = [c.args[0] for c in command.call_args_list]
        self.assertTrue(train[train.index('--array')+1].endswith('%3'))
        self.assertIn('afterok:123', pool)
        self.assertIn('--kill-on-invalid-dep=yes', pool)
        self.assertIn('--gpus', train)
        self.assertNotIn('--gpus', pool)
        self.assertIn(str(run), train)
        self.assertIn(str(run), pool)
        args.dry_run = True
        with patch('DEGAS_python.sweep.subprocess.run') as command, contextlib.redirect_stdout(io.StringIO()) as output:
            sweep.submit(args)
            command.assert_not_called()
            self.assertIn('afterok:ARRAY_JOB_ID', output.getvalue())
        self.fake_results(run)
        (run/'validation/task_000007/complete.json').unlink()
        with patch('DEGAS_python.sweep.subprocess.run') as command, contextlib.redirect_stdout(io.StringIO()) as output:
            sweep.submit(args)
            command.assert_not_called()
            self.assertIn('--array 7%3', output.getvalue())
        args.dry_run = False
        with patch('DEGAS_python.sweep.subprocess.run') as command:
            command.side_effect = [subprocess.CompletedProcess([], 0, '987\n'),
                                   subprocess.CalledProcessError(1, 'sbatch')]
            with self.assertRaises(subprocess.CalledProcessError):
                sweep.submit(args)
        records = [json.loads(p.read_text()) for p in (run/'slurm').glob('submission_*.json')]
        self.assertTrue(any(r['array_job_id'] == '987' and r['pending_task_ids'] == [7] for r in records))
        for script in (run/'slurm').glob('*.sh'):
            subprocess.run(['bash', '-n', str(script)], check=True)

    def test_real_exact_list_parallel_workers_and_restart(self):
        gene_file = self.root/'genes.txt'
        gene_file.write_text('gene_9\ngene_1\ngene_4\n')
        args = sweep.build_parser().parse_args(['prepare', '--output', str(self.root/'exact'),
                    '--genes', str(gene_file), '--seeds', '1', '--iters', '2', '--bootstrap', '30',
                    '--shap', '--shap-samples', '32'])
        run = sweep.prepare(args)
        def call(task_id):
            result = subprocess.run([sys.executable, '-m', 'DEGAS_python.sweep', 'worker', '--run', str(run),
                                     '--phase', 'validation', '--task-id', str(task_id)], capture_output=True,
                                     text=True, env=dict(os.environ, MPLBACKEND='Agg', OMP_NUM_THREADS='1'))
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        with ThreadPoolExecutor(max_workers=2) as pool:
            list(pool.map(call, range(3)))
        before = list((run/'validation/task_000000').glob('attempt_*'))
        sweep.worker(run, 'validation', 0)
        self.assertEqual(before, list((run/'validation/task_000000').glob('attempt_*')))
        # One failed attempt can remain without being mixed into a completed result.
        (run/'validation/task_000001/attempt_failed').mkdir()
        report = sweep.collect(run, 'validation')
        self.assertEqual(pd.read_csv(report/'feature_size_selection.csv')['size'].tolist(), [3])
        sweep.worker(run, 'final', 0)
        final = sweep.collect(run, 'final')
        self.assertTrue((final/'holdout_metrics.json').exists())
        self.assertTrue((final/'shap_seed_mean_3.npz').exists())
        for path in (run/'features').glob('*.json'):
            self.assertEqual(json.loads(path.read_text()), ['gene_9', 'gene_1', 'gene_4'])

    def test_empirical_task_really_trains_thousands_of_genes(self):
        args = sweep.build_parser().parse_args(['prepare', '--output', str(self.root/'large'),
               '--feature-mode', 'empirical', '--demo-genes', '3000', '--seeds', '1', '--iters', '2', '--bootstrap', '30'])
        run = sweep.prepare(args)
        task = next(t for t in sweep._task_list(run, 'validation') if t['size'] == 3000)
        sweep.worker(run, 'validation', task['task_id'])
        attempt = sweep._completion(run, 'validation', task)
        self.assertEqual(len(json.loads((attempt/'genes.json').read_text())), 3000)
        options = json.loads(next((attempt/'model').glob('*/configs.json')).read_text())
        self.assertEqual(options['input_shape'], 3000)
        self.assertFalse((attempt/'cells_preprocessed.npy').exists())


if __name__ == '__main__':
    unittest.main()
