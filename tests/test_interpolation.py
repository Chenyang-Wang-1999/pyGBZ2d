"""Tests for continuation.interpolation — unified Hermite polynomial builders.

All returned arrays follow the numpy ``poly`` convention (highest power
first), so the tests verify them through ``np.polyval`` / ``np.polyder``.
"""

import numpy as np
import pytest

from continuation.interpolation import cubic_hermite_poly, hermite_interp_poly


def test_cubic_hermite_poly_reproduces_cubic():
    """For a genuine cubic, the Hermite poly equals the cubic itself."""
    h = 0.7
    # f(x) = x^3 - 2 x^2 + x + 1
    v0 = 1.0
    dv0 = 1.0
    v1 = h ** 3 - 2 * h ** 2 + h + 1.0
    dv1 = 3 * h ** 2 - 4 * h + 1.0

    poly = cubic_hermite_poly(h, v0, dv0, v1, dv1)
    assert poly.shape == (4,)
    assert np.allclose(poly, [1.0, -2.0, 1.0, 1.0])


def test_cubic_hermite_poly_matches_endpoint_data():
    h = 0.4
    v0, dv0, v1, dv1 = 2.0, -0.5, -1.3, 3.25
    poly = cubic_hermite_poly(h, v0, dv0, v1, dv1)
    deriv = np.polyder(poly)

    assert np.polyval(poly, 0.0) == pytest.approx(v0)
    assert np.polyval(deriv, 0.0) == pytest.approx(dv0)
    assert np.polyval(poly, h) == pytest.approx(v1)
    assert np.polyval(deriv, h) == pytest.approx(dv1)


def test_hermite_interp_poly_finite_derivs_equals_cubic():
    h = 0.3
    args = (h, 1.0, 2.0, -0.5, 4.0)
    assert len(hermite_interp_poly(*args)) == 4
    assert np.allclose(hermite_interp_poly(*args), cubic_hermite_poly(*args))


@pytest.mark.parametrize("bad", [np.inf, -np.inf, np.nan])
def test_hermite_interp_poly_divergent_dv0_falls_back_to_linear(bad):
    h = 0.2
    v0, v1 = 1.0, -0.5
    poly = hermite_interp_poly(h, v0, bad, v1, 2.0)
    assert len(poly) == 2
    assert np.allclose(poly, [(v1 - v0) / h, v0])


@pytest.mark.parametrize("bad", [np.inf, -np.inf, np.nan])
def test_hermite_interp_poly_divergent_dv1_falls_back_to_linear(bad):
    h = 0.2
    v0, v1 = 1.0, -0.5
    poly = hermite_interp_poly(h, v0, 2.0, v1, bad)
    assert len(poly) == 2
    assert np.allclose(poly, [(v1 - v0) / h, v0])


def test_linear_fallback_matches_endpoint_values():
    h = 0.6
    v0, v1 = 0.3, 0.9
    poly = hermite_interp_poly(h, v0, np.inf, v1, 1.0)
    assert np.polyval(poly, 0.0) == pytest.approx(v0)
    assert np.polyval(poly, h) == pytest.approx(v1)


def test_complex_values_supported():
    """Root-curve derivatives are complex (dβ₂/dθ = V·β₂)."""
    h = 0.25
    v0 = 1.0 + 1.0j
    v1 = -0.5 + 2.0j
    m0 = 0.3 - 0.4j
    m1 = 1.2 + 0.1j
    poly = hermite_interp_poly(h, v0, m0, v1, m1)
    assert len(poly) == 4
    assert np.polyval(poly, 0.0) == pytest.approx(v0)
    assert np.polyval(poly, h) == pytest.approx(v1)


def test_polysub_roots_give_curve_intersection():
    """The intended intersection workflow: roots(polysub(p1, p2))."""
    h = 1.0
    # p1: line y = x (slope 1 at both endpoints → degenerate cubic).
    p1 = hermite_interp_poly(h, 0.0, 1.0, 1.0, 1.0)
    # p2: constant y = 0.5 (zero slope at both endpoints).
    p2 = hermite_interp_poly(h, 0.5, 0.0, 0.5, 0.0)
    roots = np.roots(np.polysub(p1, p2))
    # The only intersection of these two lines lies at x = 0.5.
    assert any(abs(r - 0.5) < 1e-12 for r in roots)
