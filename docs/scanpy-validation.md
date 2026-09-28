# Scanpy integration validation

The Python marker implementation is unchanged from the T2D run validated on
2026-09-28: Seurat comparison retained the same 92 marker rows and top-20 genes,
with maximum absolute p-value difference 4.0844712213600545e-5. Cached versus
uncached Scanpy output matched exactly on a multi-block fixture (21 versus 7 tie
correction calls), including restoration after an exception. Full-data markers
and 400 downstream patient-validation fits completed in that analysis.

For the fork integration, R source and Rd parsing, Python compilation, and
`git diff --check` passed. `R CMD INSTALL` passed on Quartz, including loading
from both temporary and final installation locations. The installed package's
adapter and bundled script were found and invoked by the integration test.

The fresh numerical test rerun did not finish: the cluster Python process stalled
in dependency imports before marker calculation. A diagnostic traceback located
blocking reads in h5py and later pandas on shared storage. Copying h5py to local
storage allowed that import to complete, but other shared dependency reads
remained slow. The incomplete retry was stopped; it is not claimed as a new
numerical test pass. The previous validation above applies to the unchanged
marker implementation. No completed analysis outputs were altered.

Quartz test logs and installed package:
`/N/project/degas_st/MarkOne/T2D_runs/20260928_degas_fork_scanpy_test_v1`.
The reproducible tests and commands are included in the repository README.
