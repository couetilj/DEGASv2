"""
Minimal BlankCox training via the pure-Python interface.

Mirror of minimal_test.R: same inputs (sim_data.npz), same opt overrides,
same DEGAS_python.bagging_all_results() entry point — but with no
reticulate involved. Also runs the R<->Python parity check at the end.

Inputs:  sim_data.npz
Outputs: risk_py.csv (columns: index, hazard)
         parity result printed to stdout (vs risk_R.csv if present)
"""

import os

# Pin threading *before* torch is imported (DEGAS_python pulls it in).
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import torch

torch.set_num_threads(1)

import DEGAS_python

HERE = os.path.dirname(os.path.abspath(__file__))


def run_blankcox() -> pd.DataFrame:
    # Match the R path: reticulate hands DEGAS fresh, writable numpy arrays,
    # whereas np.load() returns read-only views. Torch issues warnings and
    # some in-place ops behave differently on non-writable tensors, so copy
    # here to keep parity honest.
    npz = np.load(os.path.join(HERE, "sim_data.npz"))
    pat_expr = np.array(npz["pat_expr"])
    pat_lab = np.array(npz["pat_lab"])
    sc_expr = np.array(npz["sc_expr"])
    print(
        f"pat_expr={pat_expr.shape}  pat_lab={pat_lab.shape}  sc_expr={sc_expr.shape}"
    )

    opt = DEGAS_python.BlankCox_opt.copy()
    opt["save_dir"] = os.path.join(HERE, "out_py")
    opt["data_name"] = "minimal_test_py"
    opt["tot_seeds"] = 2
    opt["tot_iters"] = 100
    opt["pat_batch_size"] = 20
    opt["batch_size"] = 100
    opt["save_freq"] = 50
    os.makedirs(opt["save_dir"], exist_ok=True)

    results = DEGAS_python.bagging_all_results(
        opt, pat_expr, pat_lab, sc_expr, sc_lab_mat=None
    )
    risk = results[["index", "hazard"]].sort_values("index").reset_index(drop=True)

    out_csv = os.path.join(HERE, "risk_py.csv")
    risk.to_csv(out_csv, index=False)
    print(f"wrote {out_csv} ({len(risk)} rows)")

    print("head(risk):")
    print(risk.head())
    centered = risk["hazard"] - risk["hazard"].median()
    print("summary(hazard - median(hazard)):")
    print(centered.describe())
    return risk


def parity_check(risk_py: pd.DataFrame) -> None:
    r_csv = os.path.join(HERE, "risk_R.csv")
    if not os.path.exists(r_csv):
        print(f"PARITY CHECK: skipped ({r_csv} not found — run minimal_test.R first)")
        return

    risk_r = pd.read_csv(r_csv).sort_values("index").reset_index(drop=True)
    if len(risk_r) != len(risk_py):
        print(
            f"PARITY CHECK: FAIL — row count mismatch "
            f"(R={len(risk_r)}, py={len(risk_py)})"
        )
        return
    if not np.array_equal(risk_r["index"].to_numpy(), risk_py["index"].to_numpy()):
        print("PARITY CHECK: FAIL — index columns differ")
        return

    # Tolerance note: even with OMP_NUM_THREADS=1 and deterministic torch,
    # running the same code under Rscript+reticulate vs standalone Python
    # produces ~1e-4 level drift on loss_1 from step 0 onward (see
    # losses.csv). The root cause is CPU-level FP non-determinism when
    # Python is embedded in another process (ASLR-dependent BLAS paths,
    # different OpenMP init state). That drift compounds across ~100
    # training steps. 1e-2 tolerance still catches wrong-model / wrong-
    # loss / wrong-inputs bugs while allowing the embedding noise.
    hazard_r = risk_r["hazard"].to_numpy()
    hazard_py = risk_py["hazard"].to_numpy()
    max_abs = float(np.max(np.abs(hazard_r - hazard_py)))
    pearson = float(np.corrcoef(hazard_r, hazard_py)[0, 1])
    ok = np.allclose(hazard_r, hazard_py, atol=1e-2, rtol=1e-2)
    verdict = "PASS" if ok else "FAIL"
    print(
        f"PARITY CHECK: max |R - py| = {max_abs:.6e}  "
        f"pearson r = {pearson:.6f} — {verdict}"
    )


if __name__ == "__main__":
    risk_py = run_blankcox()
    parity_check(risk_py)
