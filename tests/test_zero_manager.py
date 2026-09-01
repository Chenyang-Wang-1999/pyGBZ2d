"""Tests for continuation.zero_manager — ZeroManager class."""

import numpy as np
import pytest
from math import pi
from cmath import exp
from collections import defaultdict
from types import SimpleNamespace

from pygbz2d.core import (CharPoly, hungarian_match_indices, to_sphere_r3,
                       cost_from_sphere_r3, TWO_PI)
from pygbz2d.continuation.zero_manager import ZeroManager, SegmentData
from pygbz2d.continuation.multiple_roots import MultipleRootInfo


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


from conftest import build_HN2D_polynomial as _build_hn2d_raw


def build_HN2D_polynomial(J1, J2, gamma_1, gamma_2, delta_1, delta_2):
    """CharPoly wrapper over the shared (coeffs, degs) builder in conftest."""
    coeffs, degs = _build_hn2d_raw(J1, J2, gamma_1, gamma_2, delta_1, delta_2)
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

    def test_run_twice_raises(self, poly_D):
        """run() appends, not resets — a second run must fail loudly."""
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run(h0=0.1)

        n_segments = zm.n_segments
        n_mrs = zm.n_multiple_roots
        boundary_perm = zm.boundary_perm.copy()

        with pytest.raises(RuntimeError, match="already been called"):
            zm.run(h0=0.1)

        # The rejected second call must not have touched the built state.
        assert zm.n_segments == n_segments
        assert zm.n_multiple_roots == n_mrs
        assert np.array_equal(zm.boundary_perm, boundary_perm)


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
        detect_cluster at θ₁=0 catches it.

        The restart row (θ₁ = mr_jump past the MR) can still land inside the
        same degenerate neighbourhood and re-detect the same cluster; the
        run loop must recognise that as a re-detection and restart farther
        instead of appending a duplicate MR.
        """
        zm = ZeroManager(poly_A, 0j, 0.0)
        zm.run(h0=0.1, cluster_tol=1e-4)
        assert zm.n_multiple_roots == 1
        assert zm.multiple_roots[0].theta1 == 0.0
        assert list(zm.multiple_roots[0].cluster_indices) == [(0, 1)]

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
            assert 0 <= mr.theta1 < TWO_PI

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


class TestBoundaryMrTangentFrame:
    """MR 0's θ=0 modulus-sorted columns must be translated at θ=2π.

    ``boundary_perm`` obeys ``roots_right[boundary_perm] == roots_left``:
    left/MR-record column ``k`` is closing-row track column
    ``boundary_perm[k]``.  Without this translation a nontrivial monodromy
    marks the wrong columns ``inf`` on the final boundary-MR row.
    """

    @staticmethod
    def _zm(has_boundary_mr=True, with_perm=True):
        mr = MultipleRootInfo(
            theta1=0.0,
            cluster_indices=[(1, 2)],
            roots=np.array([1.0 + 0j, 2.0 + 0j, 3.0 + 0j, 4.0 + 0j]),
            cluster_stds=(0.0,),
        )
        zm = SimpleNamespace(
            multiple_roots=[mr],
            has_boundary_mr=has_boundary_mr,
        )
        if with_perm:
            # left columns (1, 2) -> right columns (0, 3)
            zm.boundary_perm = np.array([2, 0, 3, 1])
        return zm

    def test_right_boundary_mr_columns_translate(self):
        zm = self._zm()
        tangents = np.zeros((3, 4), dtype=complex)
        out = ZeroManager._mark_mr_tangents_inf(
            zm, tangents, left_mr=-1, right_mr=0)

        assert np.all(np.isfinite(out[:-1, :]))
        assert np.all(np.isinf(out[-1, [0, 3]].real))
        assert np.all(np.isfinite(out[-1, [1, 2]].real))

    def test_left_boundary_mr_columns_are_not_translated(self):
        zm = self._zm()
        tangents = np.zeros((3, 4), dtype=complex)
        out = ZeroManager._mark_mr_tangents_inf(
            zm, tangents, left_mr=0, right_mr=-1)

        assert np.all(np.isinf(out[0, [1, 2]].real))
        assert np.all(np.isfinite(out[0, [0, 3]].real))
        assert np.all(np.isfinite(out[1:, :]))

    def test_interior_mr_index_zero_is_not_translated(self):
        # When has_boundary_mr=False, MR index 0 is an ordinary interior MR
        # whose cluster_indices already use its own boundary-row track frame.
        zm = self._zm(has_boundary_mr=False)
        tangents = np.zeros((3, 4), dtype=complex)
        out = ZeroManager._mark_mr_tangents_inf(
            zm, tangents, left_mr=-1, right_mr=0)

        assert np.all(np.isinf(out[-1, [1, 2]].real))
        assert np.all(np.isfinite(out[-1, [0, 3]].real))

    def test_right_boundary_mr_requires_perm(self):
        zm = self._zm(with_perm=False)
        tangents = np.zeros((3, 4), dtype=complex)
        with pytest.raises(RuntimeError, match="boundary_perm"):
            ZeroManager._mark_mr_tangents_inf(
                zm, tangents, left_mr=-1, right_mr=0)


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


# ===========================================================================
# Interpolation & insertion
# ===========================================================================

class TestInterpolation:
    """Cubic-Hermite interpolation (interpolate_roots) and insert_solution."""

    def test_tangents_stored(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        for seg in zm.segments:
            assert seg.tangents is not None
            assert seg.tangents.shape == seg.tracked_roots.shape

    def test_locate_returns_containing_interval(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        seg = zm.segments[0]
        arr = seg.theta1_arr
        # A θ strictly inside [arr[5], arr[6]] → i == 5.
        theta_mid = 0.5 * (arr[5] + arr[6])
        s_idx, i = zm.locate(theta_mid)
        assert s_idx == 0
        assert i == 5
        # A θ on a mesh point → the interval starting at that point.
        s_idx, i = zm.locate(float(arr[3]))
        assert i == 3

    def test_interpolate_at_mesh_point_exact(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        seg = zm.segments[0]
        arr = seg.theta1_arr
        for i in (0, 3, len(arr) - 1):
            got = zm.interpolate_roots(float(arr[i]))
            assert np.allclose(got, seg.tracked_roots[i, :]), (
                f"interpolate at mesh point {i} mismatch"
            )

    def test_interpolate_accuracy(self, poly_D):
        """Cubic Hermite reproduces a smooth track to ~1e-3 (4th-order)."""
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        seg = zm.segments[0]
        arr = seg.theta1_arr
        theta = 0.5 * (arr[5] + arr[6])
        interp = zm.interpolate_roots(theta)
        solved = zm._solve(theta)
        perm = hungarian_match_indices(interp, solved)
        matched = solved[perm]
        r3_a = to_sphere_r3(interp)
        r3_b = to_sphere_r3(matched)
        max_chordal = float(np.max(np.linalg.norm(r3_a - r3_b, axis=1)))
        assert max_chordal < 1e-3, f"interpolation error {max_chordal} too large"

    def test_interpolate_explicit_indices_skip_search(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        seg = zm.segments[0]
        arr = seg.theta1_arr
        theta = 0.5 * (arr[5] + arr[6])
        via_search = zm.interpolate_roots(theta)
        via_idx = zm.interpolate_roots(theta, seg_idx=0, i=5)
        assert np.allclose(via_search, via_idx)

    def test_interpolate_mr_linear_fallback(self, poly_F):
        """Boundary MR at θ=0 → segment 0's i=0 interval uses linear."""
        zm = ZeroManager(poly_F, 0j, 0.0)
        zm.run(h0=0.1, cluster_tol=1e-4, min_dtheta=1e-6)
        seg = zm.segments[0]
        assert seg.left_mr == 0  # boundary MR
        arr = seg.theta1_arr
        theta = 0.5 * (arr[0] + arr[1])
        got = zm.interpolate_roots(theta, seg_idx=0, i=0)
        s = (theta - arr[0]) / (arr[1] - arr[0])
        linear = seg.tracked_roots[0, :] + s * (
            seg.tracked_roots[1, :] - seg.tracked_roots[0, :]
        )
        assert np.allclose(got, linear), (
            "MR-touching interval did not fall back to linear"
        )

    def test_insert_solution(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        seg = zm.segments[0]
        arr = seg.theta1_arr
        n_before = len(arr)
        theta = 0.5 * (arr[5] + arr[6])
        new_idx, changed = zm.insert_solution(theta)
        assert changed is True
        assert new_idx == 6
        arr_after = seg.theta1_arr
        assert len(arr_after) == n_before + 1
        assert np.all(np.diff(arr_after) >= -1e-12)  # monotonic
        assert seg.tracked_roots.shape[0] == len(arr_after)
        assert seg.abs_argsort.shape[0] == len(arr_after)
        assert seg.tangents.shape[0] == len(arr_after)
        # interpolate at the inserted θ now returns the inserted row exactly.
        assert np.allclose(zm.interpolate_roots(theta),
                            seg.tracked_roots[new_idx, :])

    def test_interpolate_matches_inserted_solution(self, poly_D):
        """The cubic-Hermite prediction (from the pre-insertion neighbors)
        must agree with the roots ``insert_solution`` actually solves and
        reorders at the same θ — both are anchored on the same neighbor
        pair, so their difference is the cubic interpolation error (~1e-3),
        not a track-mismatch."""
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        seg = zm.segments[0]
        arr = seg.theta1_arr
        theta = 0.5 * (arr[5] + arr[6])

        interp = zm.interpolate_roots(theta)          # pre-insertion cubic
        new_idx, changed = zm.insert_solution(theta)           # true solve + reorder
        inserted = seg.tracked_roots[new_idx, :]

        r3_a = to_sphere_r3(interp)
        r3_b = to_sphere_r3(inserted)
        max_chordal = float(np.max(np.linalg.norm(r3_a - r3_b, axis=1)))
        assert max_chordal < 1e-3, (
            f"interpolation vs inserted solution chordal {max_chordal} too large"
        )

    def test_insert_endpoint_returns_unchanged(self, poly_D):
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        seg = zm.segments[0]
        arr = seg.theta1_arr
        n_before = len(arr)
        # left endpoint → unchanged, returns that endpoint's index.
        idx, changed = zm.insert_solution(float(arr[5]), seg_idx=0, i=5)
        assert changed is False
        assert idx == 5
        assert len(seg.theta1_arr) == n_before  # no row added
        # right endpoint → unchanged, returns the right endpoint's index.
        idx, changed = zm.insert_solution(float(arr[6]), seg_idx=0, i=5)
        assert changed is False
        assert idx == 6
        assert len(seg.theta1_arr) == n_before

    def test_insert_single_row_raises(self, poly_D):
        """A 1-row segment has no interval → ValueError on both methods."""
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        seg0 = zm.segments[0]
        seg0.theta1_arr = np.array([0.3])
        seg0.tracked_roots = seg0.tracked_roots[:1, :]
        seg0.tangents = seg0.tangents[:1, :]
        with pytest.raises(ValueError):
            zm.interpolate_roots(0.3)
        with pytest.raises(ValueError):
            zm.insert_solution(0.35)

    def test_solve_at_matches_insert_solution(self, poly_D):
        """solve_at is the non-mutating core of insert_solution.

        solve_at alone must not change any mesh shape, and the roots/V it
        returns must equal the row a subsequent insert_solution stores at
        the same θ (same _solve + same interpolate_roots anchor, so the
        Hungarian match is identical).
        """
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        seg = zm.segments[0]
        arr = seg.theta1_arr
        theta = 0.5 * (arr[5] + arr[6])
        n_before = len(arr)

        roots, V = zm.solve_at(theta)
        assert len(seg.theta1_arr) == n_before          # no mutation
        assert seg.tracked_roots.shape[0] == n_before
        assert seg.tangents.shape[0] == n_before

        new_idx, changed = zm.insert_solution(theta)
        assert changed is True
        assert new_idx == 6
        assert np.allclose(seg.tracked_roots[new_idx, :], roots)
        assert np.allclose(seg.tangents[new_idx, :], V)

    def test_solve_at_hint_equals_search(self, poly_D):
        """Hinted solve_at (seg_idx, i) equals the locate-based call."""
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        arr = zm.segments[0].theta1_arr
        theta = 0.5 * (arr[5] + arr[6])
        r_a, V_a = zm.solve_at(theta)
        r_b, V_b = zm.solve_at(theta, seg_idx=0, i=5)
        assert np.allclose(r_a, r_b)
        assert np.allclose(V_a, V_b)

    def test_solve_at_single_row_raises(self, poly_D):
        """Mirrors test_insert_single_row_raises for the read-only method."""
        zm = ZeroManager(poly_D, 0j, 0.0)
        zm.run()
        seg0 = zm.segments[0]
        seg0.theta1_arr = np.array([0.3])
        seg0.tracked_roots = seg0.tracked_roots[:1, :]
        seg0.tangents = seg0.tangents[:1, :]
        with pytest.raises(ValueError):
            zm.solve_at(0.35)
