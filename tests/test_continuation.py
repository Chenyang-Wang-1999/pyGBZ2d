"""Tests for continuation — pseudo-arclength root tracking.

Tests low-level functions (compute_tangent, predict_roots, estimate_error,
arclength_step), multiple-root detection (detect_cluster,
solve_multiple_roots_in_interval), and segment integration (integrate_segment).

Uses synthetic polynomials A–F and the 2D Hatano-Nelson model.
"""

import numpy as np
import pytest
from math import pi
from cmath import exp
from collections import defaultdict

from gbz_types import (
    TWO_PI,
    CharPoly, hungarian_match_indices, sort_by_root_abs,
    to_sphere_r3, cost_from_sphere_r3,
)
from continuation.arclength import (
    compute_tangent,
    predict_roots,
    predict_roots_hermite,
    estimate_error,
    arclength_step,
    StepResult,
    StepControl,
    ZERO_THRESHOLD,
    INF_THRESHOLD,
)
from continuation.multiple_roots import (
    multiple_root_point_trigger,
    MultipleRootIntervalTrigger,
    _closest_pair_deriv,
    detect_cluster,
    solve_multiple_roots_in_interval,
)
from continuation.zero_manager import (
    integrate_segment,
    SegmentResult,
    StopReason,
)


# ===========================================================================
# Synthetic polynomial builders
# ===========================================================================

def _expand_monomial_product(roots_spec):
    poly = {(0, 0): 1.0 + 0j}
    for c, k in roots_spec:
        new_poly = defaultdict(complex)
        for (e2, e1), coeff in poly.items():
            new_poly[(e2 + 1, e1)] += coeff
            new_poly[(e2, e1 + k)] -= c * coeff
        poly = new_poly
    terms = sorted(poly.items())
    coeffs = np.array([v for _, v in terms], dtype=complex)
    degs = np.array([[0, e1, e2] for (e2, e1), _ in terms], dtype=int)
    return coeffs, degs


def make_synthetic_poly(roots_spec):
    coeffs, degs = _expand_monomial_product(roots_spec)
    return CharPoly(coeffs, degs)


def _make_poly_A():
    c = 2.0
    coeffs = np.array([1, -1, -1, -c, c, c], dtype=complex)
    degs = np.array([
        [0, 0, 2], [0, 1, 1], [0, -1, 1],
        [0, 0, 1], [0, 1, 0], [0, -1, 0],
    ], dtype=int)
    return CharPoly(coeffs, degs)


POLY_B_ROOTS_SPEC = [(1.0, 1), (1.0, 0)]
POLY_C_ROOTS_SPEC = [(1.0, 1), (1.0, -1), (1.0, 0)]
POLY_D_ROOTS_SPEC = [(0.5, 1), (1.5, 1), (2.0, 1)]


def _make_poly_E(eps=1e-8):
    return make_synthetic_poly([(1.0, 1), (1.0 + eps, 1)])


def _make_poly_F():
    coeffs = np.array([1, -2, -1, 2], dtype=complex)
    degs = np.array([
        [0, 0, 2], [0, 0, 1], [0, 1, 0], [0, 0, 0],
    ], dtype=int)
    return CharPoly(coeffs, degs)


# ===========================================================================
# HN model
# ===========================================================================

from conftest import build_HN2D_polynomial as _build_hn2d_raw


def build_HN2D_polynomial(J1, J2, gamma_1, gamma_2, delta_1, delta_2):
    """CharPoly wrapper over the shared (coeffs, degs) builder in conftest."""
    coeffs, degs = _build_hn2d_raw(J1, J2, gamma_1, gamma_2, delta_1, delta_2)
    return CharPoly(coeffs, degs)


# ===========================================================================
# Fixtures
# ===========================================================================

@pytest.fixture
def poly_A():
    return _make_poly_A()

@pytest.fixture
def poly_B():
    return make_synthetic_poly(POLY_B_ROOTS_SPEC)

@pytest.fixture
def poly_C():
    return make_synthetic_poly(POLY_C_ROOTS_SPEC)

@pytest.fixture
def poly_D():
    return make_synthetic_poly(POLY_D_ROOTS_SPEC)

@pytest.fixture
def poly_E():
    return _make_poly_E(eps=1e-8)

@pytest.fixture
def poly_F():
    return _make_poly_F()

@pytest.fixture
def hn_poly():
    return build_HN2D_polynomial(
        J1=1.0, J2=1.0, gamma_1=0.2, gamma_2=0.3,
        delta_1=0.0, delta_2=0.0,
    )

@pytest.fixture
def hn_params():
    return {"E_ref": 1.0 + 0j, "mu1": 0.2}


# ===========================================================================
# Synthetic polynomial construction
# ===========================================================================

class TestSyntheticPolynomials:
    def test_poly_A_roots_at_theta1_0(self, poly_A):
        beta1 = exp(0.0 + 1j * 0.0)
        roots = np.asarray(poly_A.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        assert len(roots) == 2
        for r in roots:
            assert abs(r - 2.0) < 1e-6

    def test_poly_A_root_passes_zero(self, poly_A):
        beta1 = exp(0.0 + 1j * pi / 2)
        roots = np.asarray(poly_A.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        roots_abs = np.sort(np.abs(roots))
        assert roots_abs[0] < 1e-10

    def test_poly_B_double_root_at_origin(self, poly_B):
        beta1 = exp(0.0 + 1j * 0.0)
        roots = np.asarray(poly_B.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        for r in roots:
            assert abs(r - 1.0) < 1e-6

    def test_poly_C_triple_root(self, poly_C):
        beta1 = exp(0.0 + 1j * 0.0)
        roots = np.asarray(poly_C.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        assert len(roots) == 3
        for r in roots:
            assert abs(r - 1.0) < 2e-5

    def test_poly_C_double_root_at_pi(self, poly_C):
        beta1 = exp(0.0 + 1j * pi)
        roots = np.asarray(poly_C.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        n_near_m1 = sum(1 for r in roots if abs(r + 1.0) < 1e-6)
        n_near_p1 = sum(1 for r in roots if abs(r - 1.0) < 1e-6)
        assert n_near_m1 == 2
        assert n_near_p1 == 1

    def test_poly_D_roots_ordering(self, poly_D):
        for theta1 in [0.0, 1.0, 2.0, 3.0]:
            beta1 = exp(0.0 + 1j * theta1)
            roots = np.asarray(poly_D.solve_roots_1d((0, 1), (0j, beta1), (2,)))
            abs_sorted = np.sort(np.abs(roots))
            assert len(abs_sorted) == 3
            np.testing.assert_allclose(
                abs_sorted / abs_sorted[0],
                [1.0, 3.0, 4.0],
                rtol=1e-10,
            )

    def test_poly_F_double_root_at_theta1_0(self, poly_F):
        beta1 = exp(0.0 + 1j * 0.0)
        roots = np.asarray(poly_F.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        assert len(roots) == 2
        for r in roots:
            assert abs(r - 1.0) < 1e-6

    def test_poly_F_roots_split_away_from_origin(self, poly_F):
        for theta1 in [0.1, 0.5, 1.0]:
            beta1 = exp(0.0 + 1j * theta1)
            roots = np.asarray(poly_F.solve_roots_1d((0, 1), (0j, beta1), (2,)))
            assert len(roots) == 2
            dist = abs(roots[0] - roots[1])
            assert dist > 1e-8


# ===========================================================================
# compute_tangent
# ===========================================================================

class TestComputeTangentSynthetic:
    def test_poly_A_tangent_diverges_near_pi_half(self, poly_A):
        mu1 = 0.0
        for theta1, min_norm in [(1.56, 10.0), (1.57, 100.0), (1.5707, 1000.0)]:
            beta1 = exp(mu1 + 1j * theta1)
            roots = np.asarray(poly_A.solve_roots_1d((0, 1), (0j, beta1), (2,)))
            V, norm_V = compute_tangent(poly_A, 0j, beta1, roots)
            assert norm_V > min_norm

    def test_poly_B_tangent_is_constant(self, poly_B):
        mu1 = 0.0
        for theta1 in [0.5, 1.0, 2.0, 3.0]:
            beta1 = exp(mu1 + 1j * theta1)
            roots = np.asarray(poly_B.solve_roots_1d((0, 1), (0j, beta1), (2,)))
            V, norm_V = compute_tangent(poly_B, 0j, beta1, roots)
            V_abs_sorted = np.sort(np.abs(V))
            np.testing.assert_allclose(V_abs_sorted, [0.0, 1.0], atol=1e-10)
            np.testing.assert_allclose(norm_V, np.sqrt(2.0), rtol=1e-10)

    def test_poly_D_all_tangents_equal(self, poly_D):
        mu1 = 0.0
        beta1 = exp(mu1 + 1j * 1.0)
        roots = np.asarray(poly_D.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        V, norm_V = compute_tangent(poly_D, 0j, beta1, roots)
        for v in V:
            np.testing.assert_allclose(v, 1j, atol=1e-10)
        np.testing.assert_allclose(norm_V, np.sqrt(1 + 3), rtol=1e-10)

    def test_finite_difference_poly_A(self, poly_A):
        mu1 = 0.0
        theta1 = 1.0
        dtheta = 1e-6
        beta1 = exp(mu1 + 1j * theta1)
        roots = np.asarray(poly_A.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        V, _ = compute_tangent(poly_A, 0j, beta1, roots)
        beta1_fwd = exp(mu1 + 1j * (theta1 + dtheta))
        roots_fwd = np.asarray(poly_A.solve_roots_1d((0, 1), (0j, beta1_fwd), (2,)))
        matches = hungarian_match_indices(roots, roots_fwd)
        for from_idx, to_idx in enumerate(matches):
            b2 = roots[from_idx]
            b2_fwd = roots_fwd[to_idx]
            if abs(b2) < 1e-14:
                continue
            ln_b2 = np.log(b2)
            ln_b2_fwd = np.log(b2_fwd)
            diff = ln_b2_fwd - ln_b2
            diff = diff - 2j * pi * round(diff.imag / (TWO_PI))
            fd = diff / dtheta
            np.testing.assert_allclose(V[from_idx], fd, rtol=1e-3, atol=1e-6)

    def test_poly_F_tangent_diverges_near_double_root(self, poly_F):
        mu1 = 0.0
        norms = []
        for theta1 in [0.1, 0.05, 0.01, 0.005]:
            beta1 = exp(mu1 + 1j * theta1)
            roots = np.asarray(poly_F.solve_roots_1d((0, 1), (0j, beta1), (2,)))
            _, norm_V = compute_tangent(poly_F, 0j, beta1, roots)
            norms.append(norm_V)
        for i in range(len(norms) - 1):
            assert norms[i] < norms[i + 1]

    def test_padding_roots_get_nan_and_norm_ignores_nan(self, poly_A):
        """0/∞ padding roots have UNDEFINED tangents (nan, not 0); norm_V
        stays finite by ignoring the nan components only."""
        mu1 = 0.0
        beta1 = exp(mu1 + 1j * 0.3)
        roots = np.asarray(poly_A.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        padded = np.append(roots, [0.0, np.inf])

        V, norm_V = compute_tangent(poly_A, 0j, beta1, padded)
        assert np.isnan(V[-2])
        assert np.isnan(V[-1])

        # norm over the padded array equals norm over the finite roots.
        _, norm_finite = compute_tangent(poly_A, 0j, beta1, roots)
        assert norm_V == pytest.approx(norm_finite)


# ===========================================================================
# predict_roots
# ===========================================================================

class TestPredictRootsSynthetic:
    def test_poly_B_prediction_is_exact(self, poly_B):
        mu1 = 0.0
        for theta1 in [0.5, 1.0, 2.0]:
            dtheta = 0.3
            beta1 = exp(mu1 + 1j * theta1)
            roots = np.asarray(poly_B.solve_roots_1d((0, 1), (0j, beta1), (2,)))
            V, _ = compute_tangent(poly_B, 0j, beta1, roots)
            predicted = predict_roots(roots, V, dtheta)
            beta1_new = exp(mu1 + 1j * (theta1 + dtheta))
            roots_actual = np.asarray(poly_B.solve_roots_1d((0, 1), (0j, beta1_new), (2,)))
            matches = hungarian_match_indices(predicted, roots_actual)
            for fi, ti in enumerate(matches):
                np.testing.assert_allclose(predicted[fi], roots_actual[ti],
                                           rtol=1e-12, atol=1e-12)

    def test_poly_A_prediction_has_error(self, poly_A):
        mu1 = 0.0
        theta1 = 1.0
        dtheta = 0.1
        beta1 = exp(mu1 + 1j * theta1)
        roots = np.asarray(poly_A.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        V, _ = compute_tangent(poly_A, 0j, beta1, roots)
        predicted = predict_roots(roots, V, dtheta)
        beta1_new = exp(mu1 + 1j * (theta1 + dtheta))
        roots_actual = np.asarray(poly_A.solve_roots_1d((0, 1), (0j, beta1_new), (2,)))
        err = estimate_error(predicted, roots_actual)
        assert err > 1e-12


# predict_roots_hermite
# ===========================================================================

class TestPredictRootsHermiteSynthetic:
    """Two-endpoint cubic Hermite prediction: exact for cubic-in-θ tracks,
    with lerp / hold-fixed fallbacks for singular tracks."""

    def test_cubic_track_is_exact(self):
        """A track β₂(θ) = a + b·θ + c·θ² + d·θ³ is reproduced exactly by a
        two-endpoint Hermite given (value, derivative) at both ends."""

        def beta(theta):
            return 1.0 + 0.5 * theta - 0.2 * theta ** 2 + 0.1 * theta ** 3

        def dbeta(theta):
            return 0.5 - 0.4 * theta + 0.3 * theta ** 2

        # V = dβ/dθ / β = (dβ/dθ) / β
        theta0, theta1 = 0.0, 1.0
        theta_target = 0.37
        roots0 = np.array([beta(theta0)])
        roots1 = np.array([beta(theta1)])
        V0 = np.array([dbeta(theta0) / beta(theta0)])
        V1 = np.array([dbeta(theta1) / beta(theta1)])

        predicted = predict_roots_hermite(
            theta_target, theta0, roots0, V0, theta1, roots1, V1,
        )
        np.testing.assert_allclose(predicted[0], beta(theta_target),
                                    atol=1e-12)

    def test_singtrack_singular_root_held_fixed(self):
        """A 0/∞ padding root (|β| < ZERO_THRESHOLD) is held fixed regardless
        of the supplied tangents — same contract as predict_roots."""
        eps = 1e-12  # below ZERO_THRESHOLD (1e-6)
        roots0 = np.array([eps + 0j])
        roots1 = np.array([2.0 + 0j])
        V0 = np.array([1.0 + 0j])
        V1 = np.array([1.0 + 0j])
        predicted = predict_roots_hermite(0.5, 0.0, roots0, V0, 1.0, roots1, V1)
        np.testing.assert_allclose(predicted[0], roots0[0], atol=0.0)

    def test_single_endpoint_matches_predict_roots(self):
        """One-endpoint call is the same tangent extrapolation as predict_roots."""
        roots = np.array([2.0 + 1j, 0.5 - 0.5j])
        V = np.array([0.3 + 0.1j, -0.2 + 0.4j])
        dtheta = 0.27
        a = predict_roots(roots, V, dtheta)
        b = predict_roots_hermite(roots[0].real * 0 + dtheta,  # theta_target = dtheta
                                  0.0, roots, V)
        np.testing.assert_allclose(a, b, atol=0.0)


# ===========================================================================
# arclength_step
# ===========================================================================

class TestArclengthStepSynthetic:
    def test_step_shrinks_near_singularity(self, poly_A):
        mu1 = 0.0
        theta1_safe = 1.0
        beta1 = exp(mu1 + 1j * theta1_safe)
        roots = np.asarray(poly_A.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        result_safe = arclength_step(poly_A, 0j, mu1, theta1_safe, roots, h=0.1)
        dtheta_safe = abs(result_safe.theta1_new - theta1_safe)
        theta1_near = 1.56
        beta1 = exp(mu1 + 1j * theta1_near)
        roots = np.asarray(poly_A.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        result_near = arclength_step(poly_A, 0j, mu1, theta1_near, roots, h=0.1)
        dtheta_near = abs(result_near.theta1_new - theta1_near)
        assert dtheta_near < dtheta_safe

    def test_step_constant_for_monomial_roots(self, poly_B):
        mu1 = 0.0
        h0 = 0.1
        theta1 = 0.5
        beta1 = exp(mu1 + 1j * theta1)
        roots = np.asarray(poly_B.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        result = arclength_step(poly_B, 0j, mu1, theta1, roots, h=h0)
        assert result.accepted
        assert result.h_new >= h0


# ===========================================================================
# detect_cluster
# ===========================================================================

class TestDetectCluster:
    def test_poly_A_double_root_detected(self, poly_A):
        """detect_cluster returns list[tuple[int,...]]; empty list if none."""
        mu1 = 0.0; theta1 = 0.0
        beta1 = exp(mu1 + 1j * theta1)
        roots = np.asarray(poly_A.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        clusters = detect_cluster(roots, cluster_tol=1e-6)
        assert len(clusters) >= 1
        assert len(clusters[0]) >= 2

    def test_poly_B_double_root_detected(self, poly_B):
        mu1 = 0.0; theta1 = 0.0
        beta1 = exp(mu1 + 1j * theta1)
        roots = np.asarray(poly_B.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        clusters = detect_cluster(roots, cluster_tol=1e-6)
        assert len(clusters) >= 1

    def test_poly_C_triple_root_detected(self, poly_C):
        mu1 = 0.0; theta1 = 0.0
        beta1 = exp(mu1 + 1j * theta1)
        roots = np.asarray(poly_C.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        clusters = detect_cluster(roots, cluster_tol=1e-4)
        assert len(clusters) >= 1
        assert len(clusters[0]) >= 3

    def test_poly_A_no_false_positive_away(self, poly_A):
        mu1 = 0.0; theta1 = 1.0
        beta1 = exp(mu1 + 1j * theta1)
        roots = np.asarray(poly_A.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        clusters = detect_cluster(roots, cluster_tol=1e-6)
        assert len(clusters) == 0

    def test_poly_D_no_false_positive(self, poly_D):
        for theta1 in [0.0, 0.5, 1.0, 2.0]:
            beta1 = exp(0.0 + 1j * theta1)
            roots = np.asarray(poly_D.solve_roots_1d((0, 1), (0j, beta1), (2,)))
            clusters = detect_cluster(roots, cluster_tol=1e-6)
            assert len(clusters) == 0

    def test_poly_E_near_degenerate_detected(self, poly_E):
        mu1 = 0.0; theta1 = 0.0
        beta1 = exp(mu1 + 1j * theta1)
        roots = np.asarray(poly_E.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        clusters = detect_cluster(roots, cluster_tol=1e-6)
        assert len(clusters) >= 1

    def test_poly_F_double_root_detected(self, poly_F):
        mu1 = 0.0; theta1 = 0.0
        beta1 = exp(mu1 + 1j * theta1)
        roots = np.asarray(poly_F.solve_roots_1d((0, 1), (0j, beta1), (2,)))
        clusters = detect_cluster(roots, cluster_tol=1e-6)
        assert len(clusters) >= 1
        assert len(clusters[0]) >= 2


# ===========================================================================
# solve_multiple_roots_in_interval
# ===========================================================================

class TestSolveMultipleRootsInInterval:
    def test_refine_poly_F(self, poly_F):
        """Brent solver should find θ₁ ≈ 0 for Poly F generic double root."""
        theta1_mr = solve_multiple_roots_in_interval(
            poly_F, 0j, 0.0, -0.1, 0.1, np.array([1.0 + 0j, 1.0 + 0j]),
        )
        assert min(abs(theta1_mr), abs(theta1_mr - TWO_PI)) < 0.01

    def test_refine_poly_B(self, poly_B):
        """Non-generic double root at θ₁=0: solve in a small bracket."""
        # Use a bracket that contains 0
        beta1_left = exp(1j * (-0.05))
        roots_left = np.asarray(poly_B.solve_roots_1d((0, 1), (0j, beta1_left), (2,)))
        theta1_mr = solve_multiple_roots_in_interval(
            poly_B, 0j, 0.0, -0.05, 0.05, roots_left,
        )
        assert min(abs(theta1_mr), abs(theta1_mr - TWO_PI)) < 0.01


# ===========================================================================
# integrate_segment
# ===========================================================================

class TestIntegrateSegment:
    def test_poly_D_no_mr(self, poly_D):
        """Poly D has well-separated roots — should complete without MR."""
        roots_0 = np.asarray(poly_D.solve_roots_1d((0, 1), (0j, exp(0j)), (2,)))
        seg = integrate_segment(poly_D, 0j, 0.0, 0.0, roots_0, TWO_PI,
                                h0=0.1)
        assert seg.stop_reason == StopReason.completed
        assert seg.tracked_roots.shape[1] == 3
        assert len(seg.theta1_arr) == seg.tracked_roots.shape[0]

    def test_poly_F_stops_at_mr(self, poly_F):
        """Integrating from 0.2 to 2π should complete (MR at 0=2π is boundary)."""
        roots_0 = np.asarray(poly_F.solve_roots_1d((0, 1), (0j, exp(0.001j)), (2,)))
        seg = integrate_segment(poly_F, 0j, 0.0, 0.001, roots_0, TWO_PI,
                                h0=0.1, min_dtheta=1e-6)
        assert seg.stop_reason == StopReason.completed

    def test_poly_B_segment(self, poly_B):
        """Poly B: integrate a segment."""
        roots_0 = np.asarray(poly_B.solve_roots_1d((0, 1), (0j, exp(0.5j)), (2,)))
        seg = integrate_segment(poly_B, 0j, 0.0, 0.5, roots_0, TWO_PI,
                                h0=0.1)
        assert seg.stop_reason == StopReason.completed
        assert seg.tracked_roots.shape[1] == 2

    def test_hn_model_segment(self, hn_poly, hn_params):
        """HN model: integrate a segment."""
        E_ref = hn_params["E_ref"]
        mu1 = hn_params["mu1"]
        roots_0 = np.asarray(hn_poly.solve_roots_1d((0, 1), (E_ref, exp(mu1)), (2,)))
        seg = integrate_segment(hn_poly, E_ref, mu1, 0.0, roots_0, TWO_PI,
                                h0=0.1)
        # Interval trigger may fire on genuine close approaches.
        assert seg.stop_reason in (StopReason.completed, StopReason.multiple_root_in_interval)
        assert seg.tracked_roots.shape[1] == hn_poly.M + hn_poly.N

    def test_tracks_continuous(self, poly_D):
        """Tracked roots within a segment should not have large jumps."""
        roots_0 = np.asarray(poly_D.solve_roots_1d((0, 1), (0j, exp(0j)), (2,)))
        seg = integrate_segment(poly_D, 0j, 0.0, 0.0, roots_0, TWO_PI,
                                h0=0.1)
        tracked = seg.tracked_roots
        for i in range(len(tracked) - 1):
            for j in range(tracked.shape[1]):
                b2_c = tracked[i, j]
                b2_n = tracked[i + 1, j]
                if abs(b2_c) < 1e-14 or not np.isfinite(b2_c):
                    continue
                if abs(b2_n) < 1e-14 or not np.isfinite(b2_n):
                    continue
                pc = to_sphere_r3(np.array([b2_c]))
                pn = to_sphere_r3(np.array([b2_n]))
                dist = np.linalg.norm(pc[0] - pn[0])
                assert dist < 1.0


# ===========================================================================
# HN model low-level tests
# ===========================================================================

class TestComputeTangentHN:
    def test_shape_and_types(self, hn_poly, hn_params):
        E_ref = hn_params["E_ref"]; mu1 = hn_params["mu1"]
        beta1 = exp(mu1 + 1j * 0.0)
        roots = np.asarray(hn_poly.solve_roots_1d((0, 1), (E_ref, beta1), (2,)))
        V, norm_V = compute_tangent(hn_poly, E_ref, beta1, roots)
        assert V.shape == roots.shape
        assert isinstance(norm_V, float)
        assert norm_V >= 1.0

    def test_finite_difference_agreement(self, hn_poly, hn_params):
        E_ref = hn_params["E_ref"]; mu1 = hn_params["mu1"]
        theta1 = 0.5; dtheta = 1e-6
        beta1 = exp(mu1 + 1j * theta1)
        roots = np.asarray(hn_poly.solve_roots_1d((0, 1), (E_ref, beta1), (2,)))
        V, _ = compute_tangent(hn_poly, E_ref, beta1, roots)
        beta1_fwd = exp(mu1 + 1j * (theta1 + dtheta))
        roots_fwd = np.asarray(hn_poly.solve_roots_1d((0, 1), (E_ref, beta1_fwd), (2,)))
        matches = hungarian_match_indices(roots, roots_fwd)
        for from_idx, to_idx in enumerate(matches):
            b2 = roots[from_idx]; b2_fwd = roots_fwd[to_idx]
            if abs(b2) < ZERO_THRESHOLD or abs(b2) > INF_THRESHOLD:
                # Undefined tangent for a 0/∞ padding root: nan, not 0.
                assert np.isnan(V[from_idx])
                continue
            ln_b2 = np.log(b2); ln_b2_fwd = np.log(b2_fwd)
            diff = ln_b2_fwd - ln_b2
            diff = diff - 2j * pi * round(diff.imag / (TWO_PI))
            np.testing.assert_allclose(V[from_idx], diff / dtheta, rtol=1e-3, atol=1e-6)


class TestPredictRootsHN:
    def test_zero_displacement(self, hn_poly, hn_params):
        E_ref = hn_params["E_ref"]; mu1 = hn_params["mu1"]
        beta1 = exp(mu1 + 1j * 0.0)
        roots = np.asarray(hn_poly.solve_roots_1d((0, 1), (E_ref, beta1), (2,)))
        V, _ = compute_tangent(hn_poly, E_ref, beta1, roots)
        predicted = predict_roots(roots, V, 0.0)
        np.testing.assert_allclose(predicted, roots)

    def test_small_step_approximation(self, hn_poly, hn_params):
        E_ref = hn_params["E_ref"]; mu1 = hn_params["mu1"]
        theta1 = 0.3; dtheta = 0.01
        beta1 = exp(mu1 + 1j * theta1)
        roots = np.asarray(hn_poly.solve_roots_1d((0, 1), (E_ref, beta1), (2,)))
        V, _ = compute_tangent(hn_poly, E_ref, beta1, roots)
        predicted = predict_roots(roots, V, dtheta)
        beta1_new = exp(mu1 + 1j * (theta1 + dtheta))
        roots_actual = np.asarray(hn_poly.solve_roots_1d((0, 1), (E_ref, beta1_new), (2,)))
        matches = hungarian_match_indices(predicted, roots_actual)
        max_rel_error = 0.0
        for from_idx, to_idx in enumerate(matches):
            p = predicted[from_idx]; a = roots_actual[to_idx]
            rel_err = abs(p - a) / max(abs(a), 1e-12)
            max_rel_error = max(max_rel_error, rel_err)
        assert max_rel_error < 0.1


class TestEstimateError:
    def test_identical_roots(self):
        roots = np.array([1.0 + 0j, 2.0 + 0j, -1.0 + 0j])
        err = estimate_error(roots, roots)
        assert err < 1e-10

    def test_different_roots(self):
        r1 = np.array([1.0 + 0j, 2.0 + 0j])
        r2 = np.array([100.0 + 0j, 200.0 + 0j])
        err = estimate_error(r1, r2)
        assert err > 1.0


class TestArclengthStepHN:
    def test_single_step_accepted(self, hn_poly, hn_params):
        E_ref = hn_params["E_ref"]; mu1 = hn_params["mu1"]
        beta1 = exp(mu1 + 1j * 0.0)
        roots = np.asarray(hn_poly.solve_roots_1d((0, 1), (E_ref, beta1), (2,)))
        result = arclength_step(hn_poly, E_ref, mu1, 0.0, roots, h=0.1, ctrl=StepControl(min_step=1e-14))
        assert isinstance(result, StepResult)
        assert result.accepted
        assert result.theta1_new > 0
        assert result.h_new > 0


class TestDetectClusterHN:
    def test_no_false_positive(self, hn_poly, hn_params):
        E_ref = hn_params["E_ref"]; mu1 = hn_params["mu1"]
        beta1 = exp(mu1 + 1j * 0.0)
        roots = np.asarray(hn_poly.solve_roots_1d((0, 1), (E_ref, beta1), (2,)))
        clusters = detect_cluster(roots)
        assert len(clusters) == 0


class TestEdgeCases:
    def test_hermitian_limit(self):
        coeffs = np.array([1, -1, -1, -1, -1], dtype=complex)
        degs = np.array([
            [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1],
        ], dtype=int)
        poly = CharPoly(coeffs, degs)
        roots_0 = np.asarray(poly.solve_roots_1d((0, 1), (0j, 1.0+0j), (2,)))
        seg = integrate_segment(poly, 0j, 0.0, 0.0, roots_0, TWO_PI, h0=0.1)
        assert seg.stop_reason == StopReason.completed

    def test_empty_polynomial_handling(self):
        coeffs = np.array([1.0 + 0j])
        degs = np.array([[1, 0, 0]], dtype=int)
        poly = CharPoly(coeffs, degs)
        if poly.M == 0 and poly.N == 0:
            roots = poly.solve_roots_1d((0, 1), (1.0 + 0j, 1.0 + 0j), (2,))
            assert len(roots) == 0


# ===========================================================================
# MultipleRootIntervalTrigger — unit tests
# ===========================================================================


class TestMultipleRootIntervalTrigger:
    """Unit tests for the interval trigger.

    Uses synthetic root/V arrays to verify trigger logic without
    requiring full CharPoly integration.
    """

    def test_no_trigger_when_far_apart(self):
        """Roots far apart → state reset, no trigger."""
        trigger = MultipleRootIntervalTrigger(min_dist_threshold=0.1)
        roots = np.array([1.0 + 0j, 2.0 + 0j])  # distance 1.0 > 0.1
        V = np.array([-1.0 + 0j, 1.0 + 0j])
        ok, interval = trigger(roots, V, 0.0)
        assert not ok
        assert interval is None

    def test_trigger_on_approach_then_separate(self):
        """Same pair approaches (deriv<0) then separates (deriv>0) → trigger."""
        trigger = MultipleRootIntervalTrigger(min_dist_threshold=0.1)

        # Step 1: roots at distance 0.05 (< 0.1), approaching.
        roots1 = np.array([0.0 + 0j, 0.05 + 0j])
        V1 = np.array([1.0 + 0j, -1.0 + 0j])
        ok, interval = trigger(roots1, V1, 0.0)
        assert not ok  # first call, prev_sign=0

        # Step 2: roots at distance 0.02, separating.
        roots2 = np.array([0.015 + 0j, 0.035 + 0j])
        V2 = np.array([-1.0 + 0j, 1.0 + 0j])
        ok, interval = trigger(roots2, V2, 0.1)
        assert ok
        assert interval is not None
        assert interval[0] == 0.0
        assert interval[1] == 0.1

    def test_no_trigger_on_pair_change(self):
        """Closest pair changes identity between steps → no trigger."""
        trigger = MultipleRootIntervalTrigger(min_dist_threshold=0.1)

        # Step 1: pair (0,1) closest, approaching.
        roots1 = np.array([0.0 + 0j, 0.05 + 0j, 1.0 + 0j])
        V1 = np.array([0.0 + 0j, -1.0 + 0j, 0.0 + 0j])
        ok, _ = trigger(roots1, V1, 0.0)
        assert not ok

        # Step 2: pair (1,2) now closest (different pair), separating.
        roots2 = np.array([0.0 + 0j, 0.03 + 0j, 0.04 + 0j])
        V2 = np.array([0.0 + 0j, 0.0 + 0j, 1.0 + 0j])
        ok, _ = trigger(roots2, V2, 0.1)
        assert not ok  # pair changed → no trigger

    def test_state_reset_above_threshold(self):
        """Distance exceeds threshold → state reset, no false trigger."""
        trigger = MultipleRootIntervalTrigger(min_dist_threshold=0.1)

        # Step 1: roots close, approaching.
        roots1 = np.array([0.0 + 0j, 0.05 + 0j])
        V1 = np.array([1.0 + 0j, -1.0 + 0j])
        ok, _ = trigger(roots1, V1, 0.0)
        assert not ok

        # Step 2: roots far apart (> threshold) → state reset.
        roots2 = np.array([0.0 + 0j, 1.0 + 0j])
        V2 = np.array([0.0 + 0j, 0.0 + 0j])
        ok, _ = trigger(roots2, V2, 0.1)
        assert not ok

        # Step 3: roots close again, separating → no trigger (prev_sign=0).
        roots3 = np.array([0.04 + 0j, 0.05 + 0j])
        V3 = np.array([1.0 + 0j, 1.0 + 0j])
        ok, _ = trigger(roots3, V3, 0.2)
        assert not ok

    def test_closest_pair_deriv_returns_pair(self):
        """_closest_pair_deriv returns the correct closest-pair indices."""
        roots = np.array([1.0 + 0j, 2.0 + 0j, 1.05 + 0j])  # pair (0,2) closest
        V = np.array([0.0 + 0j, 0.0 + 0j, 0.0 + 0j])
        min_dist, sign, pair = _closest_pair_deriv(roots, V)
        assert pair == (0, 2) or pair == (2, 0)
        assert min_dist == pytest.approx(0.05)
        assert sign == 0  # all V=0 → no movement

    def test_point_trigger_basic(self):
        """Point trigger fires when dtheta < min_dtheta."""
        assert multiple_root_point_trigger(1e-12, min_dtheta=1e-10)
        assert not multiple_root_point_trigger(1e-8, min_dtheta=1e-10)
