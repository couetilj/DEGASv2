# Minimal DEGAS BlankCox Example

A minimal, end-to-end harness for training a DEGAS BlankCox survival model and
extracting per-spot risk scores — with no `DEGAS_preprocessing()` and no
`run_DEGAS_SCST()`. It calls `DEGAS_python.bagging_all_results()` directly on
pre-normalized expression matrices and a two-column survival label array. Two
execution paths (pure Python and R via reticulate) are provided so their
outputs can be compared as a parity check.

---

## Setup

Tested on macOS arm64, R 4.5.3, Python 3.13.11 (Homebrew).

### Python dependencies

```bash
pip install --upgrade setuptools
pip install numpy pandas scipy scikit-learn tqdm torch
pip install -e ../DEGAS_python          # editable install from repo
```

> **Note:** NumPy 2.x works fine with reticulate ≥ 1.40. The older `numpy<2`
> constraint only applies if you are running reticulate < 1.40.

### R dependency

```r
install.packages("reticulate")   # tested with reticulate 1.46
```

### RETICULATE_PYTHON

Tell reticulate which Python interpreter to use (the same one you used for
`pip install -e DEGAS_python`):

```bash
export RETICULATE_PYTHON=/opt/homebrew/bin/python3.13
```

Set this in your shell profile or prefix it to every `Rscript` call (see "How
to run" below).

---

## How to run

Run all three commands **from inside `minimal_example/`** (see Gotchas):

```bash
cd minimal_example/

python simulate_data.py

RETICULATE_PYTHON=/opt/homebrew/bin/python3.13 Rscript minimal_test.R

python minimal_test.py
```

The final line of `minimal_test.py` prints the parity verdict:

```
PARITY CHECK: max |R - py| = 1.11e-16  pearson r = 1.000000 — PASS
```

---

## Inputs

`simulate_data.py` writes `sim_data.npz` (inside `minimal_example/`) with
three arrays:

| Array | Shape | dtype | Values |
|-------|-------|-------|--------|
| `pat_expr` | (50, 100) | float64 | Uniform [0, 1] — bulk patient gene expression, z-scored then min-max scaled |
| `pat_lab` | (50, 2) | int64 | Column 0 = survival time (1–120 months, integer); column 1 = event status (0 = censored, 1 = event) |
| `sc_expr` | (300, 100) | float64 | Uniform [0, 1] — spatial/single-cell gene expression, same normalization as `pat_expr` |

Both expression matrices share the same 100 gene columns in the same order.
The same `sim_data.npz` is loaded by both test scripts so the parity check
compares truly identical inputs.

---

## Outputs

### `risk_R.csv` / `risk_py.csv`

Two-column CSV written by each test script:

| Column | Type | Meaning |
|--------|------|---------|
| `index` | int | Spot index (0-based row of `sc_expr`) |
| `hazard` | float | Predicted hazard score (ensemble mean across seeds) |

### `out_R/` and `out_py/`

One subdirectory per random seed:

```
out_R/
  fold_-1_random_seed_0/
    configs.json                        model hyperparameters
    losses.csv                          per-iteration loss (low_reso_loss, transfer_loss)
    high_reso_results_epoch_{N}.csv     per-spot hazard at save checkpoints
    low_reso_results_epoch_{N}.csv      per-patient predictions at save checkpoints
  fold_-1_random_seed_1/
    ...
  summary.csv                           all seeds concatenated
  summary_mean.csv                      per-spot mean hazard across seeds (columns: index, hazard, fold, seed)
```

`summary_mean.csv` is what gets returned by `bagging_all_results()`; the
`hazard` column there is the final ensemble score.

---

## Adapting to real data

```python
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import zscore
import DEGAS_python

# 1. Load your expression matrices (samples × genes)
pat_expr_raw = pd.read_csv("bulk_expr.csv", index_col=0)   # patients × all genes
sc_expr_raw  = pd.read_csv("st_expr.csv",  index_col=0)    # spots   × all genes

# 2. Z-score then min-max scale each separately so values are in [0, 1]
pat_z  = zscore(pat_expr_raw.values, axis=0)
sc_z   = zscore(sc_expr_raw.values,  axis=0)
scaler = MinMaxScaler()
pat_expr = scaler.fit_transform(pat_z)
sc_expr  = MinMaxScaler().fit_transform(sc_z)

# 3. Subset to the top 100 most variable bulk genes; apply to both matrices
gene_var   = np.var(pat_expr_raw.values, axis=0)
top_genes  = pat_expr_raw.columns[np.argsort(gene_var)[-100:]]
pat_expr   = MinMaxScaler().fit_transform(zscore(pat_expr_raw[top_genes].values, axis=0))
sc_expr    = MinMaxScaler().fit_transform(zscore(sc_expr_raw[top_genes].values,  axis=0))
#   ↑ both matrices must have the same gene columns in the same order

# 4. Survival labels: 2-column int matrix [time, status]
pat_lab = np.column_stack([
    survival_df["time"].astype(int),
    survival_df["status"].astype(int),   # 0 = censored, 1 = event
]).astype(np.int64)

# 5. Configure and call bagging_all_results
opt = DEGAS_python.BlankCox_opt.copy()
opt["save_dir"]       = "out_mydata"
opt["data_name"]      = "mydata"
opt["tot_seeds"]      = 5
opt["tot_iters"]      = 500
opt["pat_batch_size"] = 20   # must be ≤ n_patients
opt["batch_size"]     = 256
opt["save_freq"]      = 100

results = DEGAS_python.bagging_all_results(opt, pat_expr, pat_lab, sc_expr, sc_lab_mat=None)

# 6. Read per-spot hazard
hazard = results[["index", "hazard"]].sort_values("index").reset_index(drop=True)

# 7. Center by median for visualization
hazard["hazard_centered"] = hazard["hazard"] - hazard["hazard"].median()

# 8. Merge back to spatial metadata
st_meta = pd.read_csv("st_coordinates.csv")   # must have a row for each spot in sc_expr order
st_meta["hazard"] = hazard["hazard_centered"].values
```

---

## Gotchas

* **Run from `minimal_example/`, not the repo root.** When Python is launched
  from the repo root, `DEGASv2/DEGAS_python/` lands on `sys.path` as an empty
  namespace package, shadowing the pip-installed package, and `import
  DEGAS_python` resolves to the wrong (empty) module.

* **`pat_batch_size` must be ≤ `n_patients`.** Cox sampling uses
  `np.random.choice(..., replace=False)`; a batch larger than the population
  raises an error.

* **`pat_expr.shape[1] == sc_expr.shape[1]` is required.** Both matrices must
  have the same number of gene columns in the same order.

* **The published R wrapper `run_DEGAS_SCST()` references an undefined
  `n_st_classes` at `DEGAS_R/R/run_DEGAS.R:110`.** This harness sidesteps the
  issue by calling Python directly via reticulate.

* **reticulate ≥ 1.40 auto-converts pandas DataFrames to R `data.frame`.**
  The R script uses `results$index` and `results$hazard` directly. If you use
  reticulate < 1.40 (which does *not* auto-convert), replace those lines with
  `results$__getitem__("index")$to_numpy()` etc. and pass `convert = FALSE`
  when importing numpy.

* **numpy arrays must be C-contiguous and correctly typed before handing them
  to DEGAS.** The R script calls `np$ascontiguousarray(...)` on each array and
  forces `dtype = np$int64` on `pat_lab`. Without this, reticulate's
  round-trip to R matrices produces F-contiguous float64 buffers and silently
  promotes `pat_lab` from int64 to float64, desyncing from the pure-Python
  path.

* **Set `RETICULATE_PYTHON` to the interpreter used for `pip install -e
  DEGAS_python`.** Without it, reticulate may pick a system Python that
  doesn't have the package installed.
