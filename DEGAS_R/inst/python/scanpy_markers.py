"""Scanpy Wilcoxon with Seurat-compatible output/filtering for DEGAS.

Input is Seurat's log-normalized data, not counts or scaled PCA input.
Seurat provides fold changes/detection fractions to preserve installed semantics.
Scanpy asymptotic p-values can differ from Seurat (continuity/exact tests).
"""
import argparse
import json
import time
from pathlib import Path
import importlib.metadata
from contextlib import contextmanager
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
import scanpy as sc


@contextmanager
def cached_tie_correction(enabled):
    """Cache one immutable rank block in Scanpy 1.10.4's single-process call.

    That version calls _tiecorrect on the SAME ranks object for every cluster.
    Keep a strong reference (not an id) to avoid object-id reuse across chunks.
    Restore the original function even when rank_genes_groups raises.
    This CLI is single-threaded; this private-API optimization is version-pinned.
    """
    if not enabled:
        yield
        return
    if importlib.metadata.version('scanpy') != '1.10.4':
        raise RuntimeError('Tie cache requires verified Scanpy 1.10.4; disable it for other versions')
    from scanpy.tools import _rank_genes_groups as module
    original = module._tiecorrect
    previous = [None, None]
    def once(ranks):
        if ranks is not previous[0]:
            previous[:] = [ranks, original(ranks)]
        return previous[1]
    module._tiecorrect = once
    try:
        yield
    finally:
        module._tiecorrect = original


def find_markers(x, cells, genes, fc, min_pct=.25, logfc_threshold=.25,
                 return_thresh=.01, only_pos=False, cache_ties=True):
    if x.shape != (len(cells), len(genes)):
        raise ValueError('Matrix/identifier shape mismatch')
    if cells.cell_id.duplicated().any() or len(set(genes)) != len(genes):
        raise ValueError('Duplicate identifiers')
    if not np.isfinite(x.data).all() or (x.data < 0).any():
        raise ValueError('Expected finite nonnegative log-normalized expression')
    categories = list(dict.fromkeys(cells.cluster.astype(str)))
    cells = cells.copy()
    cells['cluster'] = pd.Categorical(cells.cluster.astype(str), categories=categories)
    sizes = cells.cluster.value_counts()
    if len(categories) < 2 or min(sizes.min(), len(cells)-sizes.max()) < 3:
        raise ValueError('Each cluster and its rest need at least three cells')
    data = ad.AnnData(x, obs=cells.set_index('cell_id'),
                     var=pd.DataFrame(index=pd.Index(genes, name='gene')))
    # Both signs and ALL genes: never pre-truncate before the DEGAS global rank.
    with cached_tie_correction(cache_ties):
        sc.tl.rank_genes_groups(data, 'cluster', method='wilcoxon', reference='rest',
                               use_raw=False, tie_correct=True, n_genes=len(genes),
                               rankby_abs=False, corr_method='bonferroni', pts=False)
    tables = []
    for cluster in categories:
        result = sc.get.rank_genes_groups_df(data, group=cluster)
        result = result.rename(columns={'names':'gene', 'pvals':'p_val'})
        stats = fc[fc.cluster.astype(str).eq(cluster)].copy()
        result = stats.merge(result[['gene','p_val']], on='gene', validate='one_to_one')
        # Seurat corrects over all assay genes, irrespective of prefiltering.
        result['p_val_adj'] = np.minimum(result.p_val * len(genes), 1.)
        keep = result[['pct.1','pct.2']].max(axis=1).ge(min_pct)
        keep &= (result.avg_log2FC if only_pos else result.avg_log2FC.abs()).ge(logfc_threshold)
        keep &= result.p_val.lt(return_thresh)
        result = result.loc[keep].dropna()
        tables.append(result.sort_values(['p_val','avg_log2FC'], ascending=[True,False], kind='stable'))
    return pd.concat(tables, ignore_index=True)[['p_val','avg_log2FC','pct.1','pct.2','p_val_adj','cluster','gene']]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    args = parser.parse_args(); d = args.directory; start = time.monotonic()
    cfg = json.loads((d/'parameters.json').read_text())
    genes = (d/'genes.txt').read_text().splitlines()
    cells = pd.read_csv(d/'cells.csv', dtype=str)
    x = sparse.csr_matrix((np.fromfile(d/'x.f64',dtype='<f8'),
                           np.fromfile(d/'i.i32',dtype='<i4'),
                           np.fromfile(d/'p.i32',dtype='<i4')),
                          shape=(len(cells),len(genes)))
    fc = pd.read_csv(d/'fold_changes.csv.gz', dtype={'cluster':str,'gene':str})
    result = find_markers(x,cells,genes,fc,**cfg)
    result.to_csv(d/'markers.csv',index=False)
    versions = {p:importlib.metadata.version(p) for p in ['scanpy','anndata','numpy','scipy','pandas']}
    (d/'scanpy_metrics.json').write_text(json.dumps(dict(seconds=time.monotonic()-start,
        cells=len(cells),genes=len(genes),clusters=cells.cluster.nunique(),marker_rows=len(result),
        unique_genes=result.gene.nunique(),versions=versions,tie_cache=True),indent=2))

if __name__ == '__main__':
    main()
