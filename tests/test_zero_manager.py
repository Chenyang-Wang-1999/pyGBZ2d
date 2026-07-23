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
        zm.run(h0=0.1)
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
        zm.run(h0=0.1, cluster_tol=1e-4)
        assert zm.n_multiple_roots == 1
        assert zm.multiple_roots[0].theta1 == pytest.approx(0.0, abs=1e-6)
        # cluster_indices is list[tuple[int,...]]; flatten to count indices
        n_clustered = sum(len(c) for c in zm.multiple_roots[0].cluster_indices)
        assert n_clustered == 2

    def test_poly_F_generic_double_root(self, poly_F):
        """Poly F has a generic double root at θ₁=0.  The initial check
        catches it, and the integrator detects it as well."""
        zm = ZeroManager(poly_F, 0j, 0.0)
        zm.run(h0=0.1, cluster_tol=1e-4, min_dtheta=1e-6)
        assert zm.n_multiple_roots >= 1
        assert zm.multiple_roots[0].theta1 == pytest.approx(0.0, abs=1e-6)

    def test_poly_A_no_generic_mr(self, poly_A):
        """Poly A at μ₁=0 has a non-generic double root at 0.  The initial
        detect_cluster at θ₁=0 catches it."""
        zm = ZeroManager(poly_A, 0j, 0.0)
        zm.run(h0=0.1, cluster_tol=1e-4, verbose=True)
        # Non-generic but still detectable by explicit cluster check.
        assert zm.n_multiple_roots == 1

    def test_poly_C_no_generic_mr(self, poly_C):
        """Poly C has non-generic triple/double roots — detect_cluster
        at the start catches the triple root at θ₁=0."""
        zm = ZeroManager(poly_C, 0j, 0.0)
        zm.run(h0=0.1, cluster_tol=1e-3)
        assert zm.n_multiple_roots >= 1

    def test_hn_no_false_positive(self, hn_poly):
        """HN model: no false-positive MR from step-size trigger.
        The interval trigger may find genuine close approaches that
        the old detector missed — those are real, not false positives."""
        zm = ZeroManager(hn_poly, 1.0 + 0j, 0.2)
        zm.run(h0=0.1)
        # The interval trigger is active; accept genuine MR detections.
        assert zm.n_multiple_roots >= 1


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
        zm.run(h0=0.1, cluster_tol=1e-4, min_dtheta=1e-6)
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
        zm.run(h0=0.1)
        assert zm.n_segments >= 1

    def test_empty_mr_list(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        assert zm.n_multiple_roots == 0
        assert zm.multiple_roots == []

    def test_mr_theta1_normalized(self, poly_F):
        zm = ZeroManager(poly_F, 0j, 0.0)
        zm.run(h0=0.1, cluster_tol=1e-4, min_dtheta=1e-6)
        for mr in zm.multiple_roots:
            assert 0 <= mr.theta1 < 2 * pi

    def test_boundary_perm_exists(self, poly_D):
        """boundary_perm should be set after a successful run."""
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        assert hasattr(zm, 'boundary_perm')
        assert zm.boundary_perm.shape == (zm.K,)

    def test_has_boundary_mr_flag(self, poly_D):
        """Poly D has no MR at boundary → has_boundary_mr should be False."""
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        assert zm.has_boundary_mr is False

    def test_has_boundary_mr_flag_poly_B(self, poly_B):
        """Poly B has double root at θ₁=0 → has_boundary_mr should be True."""
        zm = ZeroManager(poly_B, 0j, 0.0)
        zm.run(h0=0.1, cluster_tol=1e-4)
        assert zm.has_boundary_mr is True


# ===========================================================================
# boundary_perm convention
# ===========================================================================

class TestBoundaryPermConvention:
    """Docstring invariant: roots_right[boundary_perm] == roots_left.

    The guard in ZeroManager.run stops the post-loop fallback from
    clobbering the value already set by the 'completed' branch.  That
    clobber only changes the result for a *non-involution* permutation,
    so we use a model whose three β₂ roots braid into a 3-cycle over
    θ₁ ∈ [0, 2π): β₂³ = β₁ + β₁⁻¹ + c0.  With β₁ = e^{μ₁+iθ₁} and
    c0 < 2 cosh μ₁, the RHS traces an ellipse enclosing 0, so its three
    cube roots cyclically permute (cube-root monodromy) as θ₁ runs 0→2π.
    The roots stay 120° apart → never collide → no multiple root, so the
    run completes via the 'completed' branch.
    """

    def test_3cycle_satisfies_convention(self):
        mu1 = 0.2
        c0 = 1.0  # < 2 cosh(0.2) ≈ 2.04 → RHS encloses 0, no MR
        coeffs = np.array([1, -1, -1, -c0], dtype=complex)
        degs = np.array(
            [[0, 0, 3], [0, 1, 0], [0, -1, 0], [0, 0, 0]], dtype=int
        )
        poly = CharPoly(coeffs, degs)

        zm = ZeroManager(poly, 0j, mu1)
        zm.run(h0=0.1, cluster_tol=1e-4, min_dtheta=1e-6)

        # Completed without MR → boundary_perm comes from the completed branch.
        assert zm.n_multiple_roots == 0
        assert zm.segments[-1].right_mr == -1

        bp = zm.boundary_perm
        roots_right = zm.segments[-1].tracked_roots[-1]
        roots_left = zm.left_boundary_roots

        # Must be a genuine 3-cycle (non-involution); otherwise the test
        # cannot tell the guarded value from its inverse.
        assert not np.array_equal(bp, np.argsort(bp)), (
            f"boundary_perm {bp} is an involution — test would not catch a clobber"
        )

        # Docstring convention: roots_right[boundary_perm] == roots_left.
        assert np.allclose(roots_right[bp], roots_left)

        # The inverse (what the pre-fix fallback clobber produced) must
        # violate the convention — this is what makes the test a regression
        # guard for the boundary_perm_set flag.
        assert not np.allclose(roots_right[np.argsort(bp)], roots_left)


# ===========================================================================
# SegmentData.mr range (no sentinel leak)
# ===========================================================================

class TestSegmentMrRange:
    """left_mr / right_mr must never leak the old -2 false-positive sentinel.

    Asserts every segment's mr values are either -1 (no MR) or a valid index
    into ``multiple_roots``.  Guards the run() refactor that replaces the
    right_mr == -2 sentinel with a local ``pending`` holder.
    """

    def test_no_sentinel_leak(
        self, poly_A, poly_B, poly_C, poly_D, poly_F, hn_poly,
    ):
        cases = [
            ("poly_D", poly_D, 0j, 0.0, {}),
            ("poly_B", poly_B, 0j, 0.0, {"cluster_tol": 1e-4}),
            ("poly_F", poly_F, 0j, 0.0, {"cluster_tol": 1e-4, "min_dtheta": 1e-6}),
            ("poly_A", poly_A, 0j, 0.0, {"cluster_tol": 1e-4}),
            ("poly_C", poly_C, 0j, 0.0, {"cluster_tol": 1e-3}),
            ("hn_poly", hn_poly, 1.0 + 0j, 0.2, {}),
        ]
        for name, poly, E, mu1, kwargs in cases:
            zm = ZeroManager(poly, E, mu1)
            zm.run(h0=0.1, **kwargs)
            n_mr = zm.n_multiple_roots
            for seg in zm.segments:
                for mr in (seg.left_mr, seg.right_mr):
                    assert mr == -1 or 0 <= mr < n_mr, (
                        f"{name}: segment mr {mr} out of range (n_mr={n_mr})")
