import sys,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scanpy_markers import find_markers,cached_tie_correction
from scanpy.tools import _rank_genes_groups as module
p=Path(sys.argv[1]);genes=(p/'genes.txt').read_text().splitlines();cells=pd.read_csv(p/'cells.csv',dtype=str)
x=sparse.csr_matrix((np.fromfile(p/'x.f64',dtype='<f8'),np.fromfile(p/'i.i32',dtype='<i4'),np.fromfile(p/'p.i32',dtype='<i4')),shape=(len(cells),len(genes)))
fc=pd.read_csv(p/'fold_changes.csv.gz',dtype={'cluster':str,'gene':str})
# Force multiple rank chunks; catch accidental reuse of a prior chunk's answer.
# The 1.10.4 _ranks uses a local constant, so supply its mathematically identical
# rank generator with 13-gene chunks for this bounded test only.
original_ranks=module._ranks
original_tie=module._tiecorrect
calls=[0]
def chunked(X,*args):
    assert not args
    for lo in range(0,X.shape[1],13):
        hi=min(lo+13,X.shape[1])
        yield pd.DataFrame(X[:,lo:hi].toarray()).rank(),lo,hi

def counting(ranks):
    calls[0]+=1
    return original_tie(ranks)
module._ranks=chunked;module._tiecorrect=counting
try:
    plain=find_markers(x,cells,genes,fc,cache_ties=False);plain_calls=calls[0]
    calls[0]=0
    cached=find_markers(x,cells,genes,fc,cache_ties=True);cached_calls=calls[0]
    pd.testing.assert_frame_equal(plain,cached,check_exact=True)
    assert plain_calls==21 and cached_calls==7
    assert module._tiecorrect is counting
    try:
        with cached_tie_correction(True):raise ValueError('test restoration')
    except ValueError:pass
    assert module._tiecorrect is counting
finally:
    module._ranks=original_ranks;module._tiecorrect=original_tie
print(json.dumps(dict(status='passed',exact_frame_equality=True,rows=len(cached),plain_tie_calls=plain_calls,cached_tie_calls=cached_calls,exception_restoration=True)))
