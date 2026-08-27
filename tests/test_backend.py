"""Dual-backend validation: NumpyLaurent vs PolyToolsLaurent.

The numpy backend must be numerically interchangeable with the compiled
poly_tools backend on every surface CharPoly consumes: minor degrees,
value/partial/second-partial evaluation, and 1-D root solving (including
the degree-deficiency padding channel).  These tests are the acceptance
gate recorded in TODO/charpoly-pluggable-backend.md.
"""

import numpy as np
import pytest
from cmath import exp

from pygbz2d import backend as bz_backend
from pygbz2d.backend import (
    NumpyLaurent,
    PolyToolsLaurent,
    make_laurent,
)
from pygbz2d.core import CharPoly, TWO_PI

try:
    import poly_tools  # noqa: F401
    HAS_POLY_TOOLS = True
except Exception:
    HAS_POLY_TOOLS = False

needs_poly_tools = pytest.mark.skipif(
    not HAS_POLY_TOOLS, reason="poly_tools not importable")


# ---------------------------------------------------------------------------
# Shared models
# ---------------------------------------------------------------------------

def hn2d(basis):
    from conftest import build_HN2D_polynomial
    return build_HN2D_polynomial(
        1.0, 0.7, 0.35, -0.2, 0.1, -0.05, basis=basis)


def random_laurent(n_terms=24, dmax=3, seed=42):
    rng = np.random.default_rng(seed)
    degs = rng.integers(-dmax, dmax + 1, size=(n_terms, 3))
    degs = np.unique(degs, axis=0)          # unique degree tuples
    coeffs = rng.normal(size=len(degs)) + 1j * rng.normal(size=len(degs))
    return coeffs.astype(complex), degs.astype(int)


MODELS = {
    "hn2d_10": hn2d("10"),
    "hn2d_11": hn2d("11"),
    "random24": random_laurent(),
}
MODEL_IDS = list(MODELS)


# ---------------------------------------------------------------------------
# Numerical equivalence
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("model", MODEL_IDS)
@needs_poly_tools
def test_minor_degrees_match(model):
    coeffs, degs = MODELS[model]
    cp_pt = CharPoly(coeffs, degs, backend="poly_tools")
    cp_np = CharPoly(coeffs, degs, backend="numpy")
    assert (cp_pt.M, cp_pt.N) == (cp_np.M, cp_np.N)
    assert cp_pt.dim == cp_np.dim == 3


@pytest.mark.parametrize("model", MODEL_IDS)
@needs_poly_tools
def test_eval_and_partials_match(model):
    coeffs, degs = MODELS[model]
    cp_pt = CharPoly(coeffs, degs, backend="poly_tools")
    cp_np = CharPoly(coeffs, degs, backend="numpy")
    rng = np.random.default_rng(7)

    for _ in range(50):
        var = (rng.normal() * 2 + 1j * rng.normal(),
               exp(rng.normal() * 0.5 + 1j * rng.uniform(0, TWO_PI)),
               exp(rng.normal() * 0.5 + 1j * rng.uniform(0, TWO_PI)))
        f_pt, f_np = cp_pt.eval_val(var), cp_np.eval_val(var)
        assert abs(f_pt - f_np) <= 1e-12 * max(1.0, abs(f_pt))

        p_pt = cp_pt.eval_partials(var)
        p_np = cp_np.eval_partials(var)
        for a, b in zip(p_pt, p_np):
            assert abs(a - b) <= 1e-12 * max(1.0, abs(a))

        s_pt = cp_pt.eval_partials_2(var, 1, 2)
        s_np = cp_np.eval_partials_2(var, 1, 2)
        assert abs(s_pt - s_np) <= 1e-12 * max(1.0, abs(s_pt))


@pytest.mark.parametrize("model", MODEL_IDS)
@needs_poly_tools
def test_solve_roots_1d_match(model):
    coeffs, degs = MODELS[model]
    cp_pt = CharPoly(coeffs, degs, backend="poly_tools")
    cp_np = CharPoly(coeffs, degs, backend="numpy")
    rng = np.random.default_rng(11)

    for _ in range(25):
        beta1 = exp(rng.normal() * 0.5 + 1j * rng.uniform(0, TWO_PI))
        r_pt = cp_pt.solve_roots_1d((0, 1), (1.1 + 0.2j, beta1), (2,))
        r_np = cp_np.solve_roots_1d((0, 1), (1.1 + 0.2j, beta1), (2,))
        assert len(r_pt) == len(r_np)
        if len(r_pt):
            a = np.sort_complex(r_pt)
            b = np.sort_complex(r_np)
            assert np.max(np.abs(a - b)) < 1e-9


# ---------------------------------------------------------------------------
# NumpyLaurent semantics (the three details that must never regress)
# ---------------------------------------------------------------------------

def test_same_degree_terms_merge():
    # Two beta2^0 terms (E and beta1^{-1}) must merge into one numerator
    # term after partial evaluation — the coefficient assembly writes one
    # coefficient per degree.
    la = NumpyLaurent(
        np.array([2.0 + 0j, 3.0 + 0j, 1.0 + 0j]),
        np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]], dtype=int),
    )
    coeffs, degs, denom = la.partial_terms_1d((1.0 + 0j, 1.0 + 0j), (0, 1), (2,))
    # With E=1, beta1=1: terms 2 + 3 merge at degree 0; term 1 stays at 1.
    assert degs == [0, 1]
    assert denom == 0
    assert set(coeffs) == {5.0 + 0j, 1.0 + 0j}


def test_exact_zero_coeff_dropped():
    la = NumpyLaurent(
        np.array([0.0 + 0j, 1.0 + 0j]),
        np.array([[0, 0, -1], [0, 0, 1]], dtype=int),
    )
    # Zero coefficient dropped at construction.
    assert len(la.coeffs) == 1

    # Two same-degree terms cancelling at partial evaluation are dropped
    # entirely (the merged-to-zero channel).
    la2 = NumpyLaurent(
        np.array([1.0 + 0j, -1.0 + 0j]),
        np.array([[0, 0, -1], [0, 0, -1]], dtype=int),
    )
    coeffs, degs, denom = la2.partial_terms_1d((1.0,), (0,), (2,))
    assert coeffs == [] and degs == [] and denom == 0


def test_denom_orders_clear_negative_powers_only():
    la = NumpyLaurent(
        np.array([1.0 + 0j, 2.0 + 0j, 3.0 + 0j]),
        np.array([[0, -2, -1], [0, 0, 1], [0, 1, 2]], dtype=int),
    )
    assert list(la.denom_orders) == [0, 2, 1]   # max(0, -min_deg) per var
    assert la.num_max_degrees() == (0, 3, 3)    # cleared numerator maxima


@needs_poly_tools
def test_empty_partial_polynomial_both_backends():
    # Every term carries E^1, so E=0 annihilates the polynomial: both
    # backends must report the empty-term channel identically.  With
    # M=1, N=1 the padding rule then materialises one 0-root (from the
    # cleared β₂ denominator) and one ∞-root (degree deficiency).
    coeffs = np.array([1.0 + 0j, -1.0 + 0j])
    degs = np.array([[1, 0, -1], [1, 1, 1]], dtype=int)
    cp_pt = CharPoly(coeffs, degs, backend="poly_tools")
    cp_np = CharPoly(coeffs, degs, backend="numpy")
    assert (cp_pt.M, cp_pt.N) == (cp_np.M, cp_np.N) == (1, 1)
    args = ((0, 1), (0.0 + 0j, 1.0 + 0j), (2,))
    r_pt = cp_pt.solve_roots_1d(*args)
    r_np = cp_np.solve_roots_1d(*args)
    assert len(r_pt) == len(r_np) == 2
    assert np.array_equal(np.isinf(r_pt), np.isinf(r_np))
    assert np.count_nonzero(r_pt == 0) == np.count_nonzero(r_np == 0) == 1
    assert np.count_nonzero(np.isinf(r_pt)) == 1


# ---------------------------------------------------------------------------
# Backend selection logic
# ---------------------------------------------------------------------------

def test_explicit_numpy_backend():
    coeffs, degs = MODELS["hn2d_11"]
    la = make_laurent(coeffs, degs, "numpy")
    assert isinstance(la, NumpyLaurent)


@needs_poly_tools
def test_explicit_poly_tools_backend():
    coeffs, degs = MODELS["hn2d_11"]
    la = make_laurent(coeffs, degs, "poly_tools")
    assert isinstance(la, PolyToolsLaurent)


@needs_poly_tools
def test_auto_prefers_poly_tools(monkeypatch):
    monkeypatch.delenv("POLY_BACKEND", raising=False)
    monkeypatch.setattr(bz_backend, "_POLY_TOOLS_AVAILABLE", True)
    coeffs, degs = MODELS["hn2d_11"]
    assert isinstance(make_laurent(coeffs, degs), PolyToolsLaurent)


@needs_poly_tools
def test_auto_falls_back_to_numpy_with_one_warning(monkeypatch):
    monkeypatch.delenv("POLY_BACKEND", raising=False)
    monkeypatch.setattr(bz_backend, "_POLY_TOOLS_AVAILABLE", False)
    monkeypatch.setattr(bz_backend, "_numpy_fallback_warned", False)
    coeffs, degs = MODELS["hn2d_11"]

    with pytest.warns(UserWarning, match="pure-numpy Laurent backend"):
        la = make_laurent(coeffs, degs)
    assert isinstance(la, NumpyLaurent)

    # One-time: the second construction stays silent.
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        make_laurent(coeffs, degs)


@needs_poly_tools
def test_missing_poly_tools_explicit_backend_raises(monkeypatch):
    monkeypatch.setattr(bz_backend, "_POLY_TOOLS_AVAILABLE", False)
    coeffs, degs = MODELS["hn2d_11"]
    with pytest.raises(ImportError, match="poly_tools"):
        make_laurent(coeffs, degs, "poly_tools")


def test_env_var_backend(monkeypatch):
    coeffs, degs = MODELS["hn2d_11"]
    monkeypatch.setenv("POLY_BACKEND", "numpy")
    assert isinstance(make_laurent(coeffs, degs), NumpyLaurent)
    # Explicit argument beats the environment variable.
    la = make_laurent(coeffs, degs, "numpy")
    assert isinstance(la, NumpyLaurent)


def test_unknown_backend_name_raises():
    coeffs, degs = MODELS["hn2d_11"]
    with pytest.raises(ValueError, match="unknown backend"):
        make_laurent(coeffs, degs, "sympy")


def test_custom_backend_class():
    coeffs, degs = MODELS["hn2d_11"]
    cp = CharPoly(coeffs, degs, backend=NumpyLaurent)
    assert cp.backend_name == "NumpyLaurent"
    beta1 = exp(0.2 + 0.7j)
    roots = cp.solve_roots_1d((0, 1), (1.1 + 0.2j, beta1), (2,))
    assert len(roots) == cp.M + cp.N
