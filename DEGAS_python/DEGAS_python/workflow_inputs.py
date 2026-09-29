"""CSV inputs and deterministic synthetic data for the validation workflow."""
import csv
import numpy as np
import pandas as pd


def demo(n_genes=40):
    rng = np.random.default_rng(42)
    y = np.tile([0, 1], 60)
    bulk = rng.poisson(10, (120, n_genes)).astype(float)
    bulk[:, :min(8, n_genes)] += y[:, None] * rng.poisson(8, (120, min(8, n_genes)))
    cells = rng.poisson(10, (100, n_genes)).astype(float)
    cells[:50, :min(8, n_genes)] += rng.poisson(8, (50, min(8, n_genes)))
    genes = [f'gene_{i}' for i in range(n_genes)]
    meta = pd.DataFrame(dict(patient_id=[f'patient_{i}' for i in range(120)],
                             study=np.repeat(['A', 'B', 'C', 'D'], 30), label=y))
    return (pd.DataFrame(bulk, index=meta.patient_id, columns=genes), meta,
            pd.DataFrame(cells, index=[f'cell_{i}' for i in range(100)], columns=genes))


def _read_counts(path):
    # pandas otherwise silently renames duplicate CSV column headers (gene, gene.1).
    with path.open(newline='') as stream:
        header = next(csv.reader(stream), [])
    genes = header[1:]
    if not genes or len(set(genes)) != len(genes) or any(not gene.strip() for gene in genes):
        raise ValueError(f'{path.name} requires unique, nonempty gene IDs in its header')
    return pd.read_csv(path, index_col=0, converters={0: str})


def read_inputs(folder):
    bulk = _read_counts(folder / 'bulk_counts.csv')
    cells = _read_counts(folder / 'cell_counts.csv')
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
