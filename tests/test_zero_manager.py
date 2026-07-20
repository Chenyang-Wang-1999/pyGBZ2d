"""Tests for continuation.zero_manager — ZeroManager class."""

import numpy as np
import pytest
from math import pi
from cmath import exp
from collections import defaultdict

from gbz_types import CharPoly, hungarian_match_indices, to_sphere_r3, cost_from_sphere_r3
from continuation.zero_manager import ZeroManager, SegmentData


# ===========================================================================
# Synthetic polynomial builders
# ===========================================================================

def make_synthetic_poly(roots_spec):
    poly_dict = {(0, 0): 1.0 + 0j}
    for c, k in roots_spec:
        new_poly = defaultdict(complex)
        for (e2, e1), coeff in poly_dict.items():
            new_poly[(e2 + 1, e1)] += coeff
            new_poly[(e2, e1 + k)] -= c * coeff
        poly_dict = new_poly
    terms = sorted(poly_dict.items())
    coeffs = np.array([v for _, v in terms], dtype=complex)
    degs = np.array([[0, e1, e2] for (e2, e1), _ in terms], dtype=int)
    return CharPoly(coeffs, degs)


def _make_poly_A():
    c = 2.0
    coeffs = np.array([1, -1, -1, -c, c, c], dtype=complex)
    degs = np.array([
        [0, 0, 2], [0, 1, 1], [0, -1, 1],
        [0, 0, 1], [0, 1, 0], [0, -1, 0],
    ], dtype=int)
    return CharPoly(coeffs, degs)


def _make_poly_F():
    coeffs = np.array([1, -2, -1, 2], dtype=complex)
    degs = np.array([
        [0, 0, 2], [0, 0, 1], [0, 1, 0], [0, 0, 0],
    ], dtype=int)
    return CharPoly(coeffs, degs)


def build_HN2D_polynomial(J1, J2, gamma_1, gamma_2, delta_1, delta_2):
    J11 = exp(gamma_1 + 1j * delta_1) * J1
    J12 = exp(-gamma_1 + 1j * delta_1) * np.conj(J1)
    J21 = exp(gamma_2 + 1j * delta_2) * J2
    J22 = exp(-gamma_2 + 1j * delta_2) * np.conj(J2)
    coeffs = np.array([1, -J11, -J12, -J21, -J22], dtype=complex)
    degs = np.array([
        [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1],
    ], dtype=int)
    return CharPoly(coeffs, degs)


@pytest.fixture
def poly_A():
    return _make_poly_A()

@pytest.fixture
def poly_B():
    return make_synthetic_poly([(1.0, 1), (1.0, 0)])

@pytest.fixture
def poly_C():
    return make_synthetic_poly([(1.0, 1), (1.0, -1), (1.0, 0)])

@pytest.fixture
def poly_D():
    return make_synthetic_poly([(0.5, 1), (1.5, 1), (2.0, 1)])

@pytest.fixture
def poly_F():
    return _make_poly_F()

@pytest.fixture
def hn_poly():
    return build_HN2D_polynomial(
        J1=1.0, J2=1.0, gamma_1=0.2, gamma_2=0.3,
        delta_1=0.0, delta_2=0.0,
    )


# ===========================================================================
# Basic structure tests
# ===========================================================================

class TestZeroManagerBasic:
    def test_no_crash_poly_D(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run(h0=0.1, max_step=0.5)
        assert zm.K == 3
        assert zm.n_multiple_roots == 0
        assert zm.n_segments == 1
        seg = zm.segments[0]
        assert seg.left_mr == -1
        assert seg.right_mr == -1
        assert seg.tracked_roots.shape[1] == 3

    def test_data_consistency(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        seg = zm.segments[0]
        n = len(seg.theta1_arr)
        assert seg.tracked_roots.shape == (n, zm.K)
        assert seg.abs_argsort.shape == (n, zm.K)
        for i in range(n):
            assert set(seg.abs_argsort[i]) == set(range(zm.K))


# ===========================================================================
# Multiple root detection
# ===========================================================================

class TestMultipleRootDetection:
    def test_poly_B_mr_at_start(self, poly_B):
        """Poly B has a non-generic double root at θ₁=0.  The initial
        detect_cluster check catches it."""
        zm = ZeroManager(poly_B, 0j, 0.0)
        zm.run(h0=0.1, max_step=0.5, cluster_tol=1e-4)
        assert zm.n_multiple_roots == 1
        assert zm.mr_thetas[0] == pytest.approx(0.0, abs=1e-6)
        assert zm.mr_cluster_masks[0].sum() == 2

    def test_poly_F_generic_double_root(self, poly_F):
        """Poly F has a generic double root at θ₁=0.  The initial check
        catches it, and the integrator detects it as well."""
        zm = ZeroManager(poly_F, 0j, 0.0)
        zm.run(h0=0.1, max_step=0.5, cluster_tol=1e-4, min_dtheta=1e-6)
        assert zm.n_multiple_roots >= 1
        assert zm.mr_thetas[0] == pytest.approx(0.0, abs=1e-6)

    def test_poly_A_no_generic_mr(self, poly_A):
        """Poly A at μ₁=0 has a non-generic double root at 0.  The initial
        detect_cluster at θ₁=0 catches it."""
        zm = ZeroManager(poly_A, 0j, 0.0)
        zm.run(h0=0.1, max_step=0.5, cluster_tol=1e-4)
        # Non-generic but still detectable by explicit cluster check.
        assert zm.n_multiple_roots == 1

    def test_poly_C_no_generic_mr(self, poly_C):
        """Poly C has non-generic triple/double roots — detect_cluster
        at the start catches the triple root at θ₁=0."""
        zm = ZeroManager(poly_C, 0j, 0.0)
        zm.run(h0=0.1, max_step=0.5, cluster_tol=1e-3)
        assert zm.n_multiple_roots >= 1

    def test_hn_no_false_positive(self, hn_poly):
        zm = ZeroManager(hn_poly, 1.0 + 0j, 0.2)
        zm.run(h0=0.1, max_step=0.5)
        assert zm.n_multiple_roots == 0


# ===========================================================================
# Segment tests
# ===========================================================================

class TestSegments:
    def test_one_segment_no_mr(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        assert zm.n_segments == 1

    def test_segment_mr_linking(self, poly_F):
        zm = ZeroManager(poly_F, 0j, 0.0)
        zm.run(h0=0.1, max_step=0.5, cluster_tol=1e-4, min_dtheta=1e-6)
        n_mr = zm.n_multiple_roots
        for seg in zm.segments:
            if seg.left_mr >= 0:
                assert seg.left_mr < n_mr
            if seg.right_mr >= 0:
                assert seg.right_mr < n_mr

    def test_theta1_monotonic(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        for seg in zm.segments:
            diffs = np.diff(seg.theta1_arr)
            assert np.all(diffs >= -1e-12)


# ===========================================================================
# Insert tests
# ===========================================================================

class TestInsert:
    def test_insert_returns_correct_shape(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        seg_idx, local_idx, roots = zm.insert(1.0)
        assert 0 <= seg_idx < zm.n_segments
        assert roots.shape == (zm.K,)

    def test_insert_modifies_segment(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        seg = zm.segments[0]
        n_before = len(seg.theta1_arr)
        zm.insert(1.5)
        assert len(seg.theta1_arr) == n_before + 1

    def test_insert_track_continuity(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        for theta1 in [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]:
            seg_idx, local_idx, roots = zm.insert(theta1)
            seg = zm.segments[seg_idx]
            if local_idx > 0:
                neighbor = seg.tracked_roots[local_idx - 1]
                matches = hungarian_match_indices(neighbor, roots)
                r3_n = to_sphere_r3(neighbor)
                r3_r = to_sphere_r3(roots)
                for fi, ti in matches:
                    dist = np.linalg.norm(r3_n[fi] - r3_r[ti])
                    assert dist < 1.0

    def test_insert_near_mr(self, poly_F):
        zm = ZeroManager(poly_F, 0j, 0.0)
        zm.run(h0=0.1, max_step=0.5, cluster_tol=1e-4, min_dtheta=1e-6)
        seg_idx, local_idx, roots = zm.insert(0.01)
        assert roots.shape == (zm.K,)
        if zm.n_multiple_roots > 0:
            mr_roots = zm.mr_roots[0]
            mr_theta = zm.mr_thetas[0]
            if abs(mr_theta) < 0.1:
                matches = hungarian_match_indices(mr_roots, roots)
                for fi, ti in matches:
                    dist = abs(mr_roots[fi] - roots[ti])
                    assert dist < 0.5


# ===========================================================================
# Edge case tests
# ===========================================================================

class TestEdgeCases:
    def test_hermitian_limit(self):
        coeffs = np.array([1, -1, -1, -1, -1], dtype=complex)
        degs = np.array([
            [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1],
        ], dtype=int)
        poly = CharPoly(coeffs, degs)
        zm = ZeroManager(poly, 0j, 0.0)
        zm.run(h0=0.1, max_step=0.5)
        assert zm.n_segments >= 1

    def test_single_segment_after_insert(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        n_before = len(zm.segments[0].theta1_arr)
        for _ in range(5):
            zm.insert(np.random.uniform(0, 2 * pi))
        assert len(zm.segments[0].theta1_arr) == n_before + 5

    def test_insert_all_around_circle(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        for theta1 in np.linspace(0, 2 * pi, 20, endpoint=False):
            seg_idx, _, roots = zm.insert(theta1)
            assert roots.shape == (zm.K,)

    def test_empty_mr_list(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        assert zm.mr_thetas == []
        assert zm.n_multiple_roots == 0

    def test_mr_theta1_normalized(self, poly_F):
        zm = ZeroManager(poly_F, 0j, 0.0)
        zm.run(h0=0.1, max_step=0.5, cluster_tol=1e-4, min_dtheta=1e-6)
        for t in zm.mr_thetas:
            assert 0 <= t < 2 * pi
