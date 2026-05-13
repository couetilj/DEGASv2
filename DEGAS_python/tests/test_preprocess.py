"""Tests for DEGAS_python.preprocess against the R reference."""

from __future__ import annotations

import pathlib

import numpy as np
import pytest

from DEGAS_python.preprocess import normalize_scale, preprocess_counts


HERE = pathlib.Path(__file__).parent
REF = HERE / "r_reference"


def test_already_log_equivalence() -> None:
    """preprocess_counts(raw) must equal preprocess_counts(log2(raw+1), already_log=True).

    This contract is the entire justification for the already_log branch.
    Without it, UCSC Xena bulk gets double-log-transformed.
    """
    rng = np.random.default_rng(0)
    raw = rng.integers(0, 10_000, size=(20, 50)).astype(np.float64)
    a = preprocess_counts(raw)
    b = preprocess_counts(np.log2(raw + 1.0), already_log=True)
    np.testing.assert_allclose(a, b, rtol=1e-9, atol=1e-12)


def test_output_in_unit_interval() -> None:
    """Row mins are exactly 0; row maxes are <= 1 (slightly below 1 with +1e-3)."""
    rng = np.random.default_rng(0)
    raw = rng.integers(0, 10_000, size=(20, 50)).astype(np.float64)
    out = preprocess_counts(raw)
    assert out.shape == raw.shape
    assert (out >= 0).all()
    assert (out <= 1).all()
    np.testing.assert_allclose(out.min(axis=1), 0.0, atol=1e-12)


def test_normalize_scale_is_per_row() -> None:
    """Each row's output depends only on that row, not on the rest of the matrix."""
    rng = np.random.default_rng(0)
    X = rng.standard_normal((5, 8))
    out_full = normalize_scale(X)
    # Normalizing just row 0 alone (as a 1-row matrix) must give the same output for row 0.
    out_alone = normalize_scale(X[:1])
    np.testing.assert_allclose(out_full[0:1], out_alone, rtol=1e-12, atol=1e-15)


def test_rejects_non_2d_input() -> None:
    with pytest.raises(ValueError):
        normalize_scale(np.zeros(5))


@pytest.mark.skipif(
    not (REF / "fixture_input.csv").exists(),
    reason="R reference fixture not generated — run r_reference/generate_reference.R",
)
def test_matches_r_degas_reference() -> None:
    """Python preprocess_counts ~ R preprocessCounts on a deterministic fixture.

    R reads (genes x samples) and returns (samples x genes). Load the R-shape
    input, transpose to Python convention, and assert agreement.
    """
    raw_R = np.loadtxt(REF / "fixture_input.csv", delimiter=",")
    expected = np.loadtxt(REF / "fixture_expected_output.csv", delimiter=",")
    assert raw_R.shape == (50, 20), f"R input fixture shape {raw_R.shape}"
    assert expected.shape == (20, 50), f"R output fixture shape {expected.shape}"
    out = preprocess_counts(raw_R.T)
    np.testing.assert_allclose(out, expected, rtol=1e-6, atol=1e-9)
