"""
Generate a fixed, deterministic synthetic dataset for BlankCox training.

Writes `sim_data.npz` at the repo root with four arrays:
    pat_expr : (n_pat,  n_gene) float64, uniform [0, 1]
    pat_lab  : (n_pat, 2)       int64,   col 0 = survival time, col 1 = status
    sc_expr  : (n_spot, n_gene) float64, uniform [0, 1]

The same `sim_data.npz` is loaded by both minimal_test.R and minimal_test.py
so the R/Python parity check compares training on bit-identical inputs.
"""

import os

import numpy as np

N_PAT = 50
N_SPOT = 300
N_GENE = 100
SEED = 1234

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sim_data.npz")


def main() -> None:
    rng = np.random.default_rng(SEED)

    pat_expr = rng.uniform(0.0, 1.0, size=(N_PAT, N_GENE))
    sc_expr = rng.uniform(0.0, 1.0, size=(N_SPOT, N_GENE))

    time = rng.integers(1, 121, size=N_PAT)
    status = rng.binomial(1, 0.6, size=N_PAT)
    pat_lab = np.column_stack([time, status]).astype(np.int64)

    assert pat_expr.min() >= 0.0 and pat_expr.max() <= 1.0
    assert sc_expr.min() >= 0.0 and sc_expr.max() <= 1.0
    assert pat_expr.shape[1] == sc_expr.shape[1]
    assert pat_lab.shape == (N_PAT, 2)

    np.savez(OUT, pat_expr=pat_expr, pat_lab=pat_lab, sc_expr=sc_expr)
    print(
        f"wrote {OUT}: pat_expr={pat_expr.shape}, "
        f"pat_lab={pat_lab.shape}, sc_expr={sc_expr.shape}"
    )


if __name__ == "__main__":
    main()
