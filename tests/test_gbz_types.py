"""Direct tests for the shared primitives in bfgbz2d/core.py.

The old test_bfgbz2d/core.py was deleted in commit a118614 without a
replacement; these tools (cyclic intervals, chordal costs, Hungarian
matching, probe ladder, zero-plateau helpers) are the seam-boundary
semantic core shared by both GBZ modules and deserve direct coverage.
"""

import math
import numpy as np
import pytest

from bfgbz2d.core import (
    TWO_PI,
    PointSubset, LineSubset, GBZResult,
    sort_by_root_abs, to_sphere_r3, chordal_cost_matrix,
    hungarian_match_indices, find_cyclic_true_intervals,
    generate_probe_steps, circ_dist,
    check_points_clustered_on_torus, probe_zero_plateau,
    JoinableLinePiece, is_mr_cluster_endpoint,
    get_minor_degrees,
)
from bfgbz2d.core import CharPoly


# ---------------------------------------------------------------------------
# PointSubset / LineSubset / GBZResult
# ---------------------------------------------------------------------------

class TestDataclasses:
    def test_point_subset_theta1_wraps_negative(self):
        p = PointSubset(E=0j, beta1=complex(math.exp(0.3) * math.cos(-0.5),
                                            math.exp(0.3) * math.sin(-0.5)),
                        beta2=1 + 0j)
        assert p.mu1 == pytest.approx(0.3)
        assert p.theta1 == pytest.approx(TWO_PI - 0.5)

    def test_line_subset_width_wraps(self):
        # theta from 3π/2 to π/2 (crossing the seam): width = π, not -π
        ls = LineSubset(E=0j, mu1=0.0,
                        theta1_arr=np.array([3 * math.pi / 2, TWO_PI - 0.1,
                                             0.0, math.pi / 2]),
                        beta2_arr=np.ones(4))
        assert ls.theta1_width == pytest.approx(math.pi)

    def test_line_subset_validation(self):
        with pytest.raises(ValueError):
            LineSubset(E=0j, mu1=0.0,
                       theta1_arr=np.zeros((2, 2)), beta2_arr=np.zeros(4))
        with pytest.raises(ValueError):
            LineSubset(E=0j, mu1=0.0,
                       theta1_arr=np.zeros(3), beta2_arr=np.zeros(4))

    def test_gbz_result_is_gbz_semantics(self):
        empty = GBZResult(E_ref=0j)
        assert empty.subsets == []
        assert not empty.is_gbz
        # is_continuum alone means in-spectrum even with no subsets
        cont = GBZResult(E_ref=0j, is_continuum=True)
        assert cont.is_gbz
        failed = GBZResult(E_ref=0j, success=False, index=(2, 0))
        assert not failed.is_gbz

    def test_joinable_line_piece_carries_span(self):
        p = JoinableLinePiece(0j, 0.0, np.linspace(0, 1, 3), np.ones(3),
                              ml=0, mr=2)
        assert isinstance(p, LineSubset)
        assert (p.ml, p.mr) == (0, 2)


# ---------------------------------------------------------------------------
# Root-array primitives
# ---------------------------------------------------------------------------

class TestRootPrimitives:
    def test_sort_by_root_abs(self):
        r = np.array([3 + 0j, 1 + 0j, 2 + 0j])
        assert list(np.abs(sort_by_root_abs(r))) == [1, 2, 3]

    def test_to_sphere_r3_north_pole_for_inf(self):
        r3 = to_sphere_r3(np.array([0 + 0j, np.inf + 0j, 1 + 0j]))
        assert tuple(r3[1]) == (0.0, 0.0, 1.0)          # ∞ → north pole
        assert r3[0, 2] == pytest.approx(-1.0)          # 0 → south pole
        assert np.allclose(np.linalg.norm(r3, axis=1), 1.0)

    def test_chordal_cost_matrix_symmetric(self):
        a = np.array([1 + 1j, 0.1 + 0j])
        b = np.array([2 + 0j, 0.1 + 0.01j])
        c = chordal_cost_matrix(a, b)
        assert c.shape == (2, 2)
        assert c[0, 0] >= 0
        # identical arrays → zero diagonal
        c2 = chordal_cost_matrix(a, a)
        assert np.allclose(np.diag(c2), 0.0)

    def test_hungarian_match_indices_perm_semantics(self):
        # roots_to was permuted from roots_from by perm p: matching must
        # recover p exactly.
        rng = np.random.default_rng(42)
        base = rng.normal(size=6) + 1j * rng.normal(size=6)
        p = rng.permutation(6)
        perm = hungarian_match_indices(base, base[p])
        assert np.array_equal(base[p][perm], base)

    def test_hungarian_match_rejects_nan(self):
        with pytest.raises(ValueError):
            hungarian_match_indices(np.array([np.nan + 0j, 1 + 0j]),
                                    np.array([1 + 0j, 2 + 0j]))
        with pytest.raises(ValueError):
            hungarian_match_indices(np.array([1 + 0j, 2 + 0j]),
                                    np.array([np.nan + 0j, 1 + 0j]))


# ---------------------------------------------------------------------------
# Cyclic interval semantics (the seam core)
# ---------------------------------------------------------------------------

class TestCyclicIntervals:
    def test_basic_interior_interval(self):
        m = np.array([False, True, True, False])
        assert find_cyclic_true_intervals(m) == [(1, 2)]

    def test_wraparound_interval(self):
        m = np.array([True, False, False, True, True])
        assert find_cyclic_true_intervals(m) == [(3, 0)]

    def test_all_true_and_all_false(self):
        assert find_cyclic_true_intervals(np.ones(5, bool)) == [(0, 4)]
        assert find_cyclic_true_intervals(np.zeros(5, bool)) == []
        assert find_cyclic_true_intervals(np.array([], dtype=bool)) == []

    def test_circ_dist(self):
        assert circ_dist(0.0, TWO_PI - 1e-9) == pytest.approx(1e-9)
        assert circ_dist(0.1, TWO_PI + 0.1) == pytest.approx(0.0)
        assert circ_dist(0.0, math.pi) == pytest.approx(math.pi)


# ---------------------------------------------------------------------------
# Probe ladder
# ---------------------------------------------------------------------------

class TestProbeLadder:
    def test_steps_sorted_positive_bounded(self):
        steps = generate_probe_steps(bracket_width=1e-3, probe_radius=0.5,
                                     zero_tol=1e-10)
        assert steps == sorted(steps)
        assert all(s > 0 for s in steps)
        assert steps[0] >= 1e-12
        assert steps[-1] <= 0.5 * (1 + 1e-9)
        # the bracket-derived steps must be present
        assert any(abs(s - 1e-3) < 1e-12 for s in steps)

    def test_probe_zero_plateau_found_stops_ladder(self):
        calls = []

        def evaluator(mu1_probe):
            calls.append(mu1_probe)
            # plateau on the RIGHT side of the second step only
            found = mu1_probe > 0.0 and len(calls) > 2
            return {"success": True, "is_continuum": False,
                    "is_plateau": found}

        res = probe_zero_plateau(0.0, (0.0, 1.0), zero_tol=1e-10,
                                 probe_radius=1.0, bracket_width=1.0,
                                 evaluator=evaluator)
        assert res["found"] is True
        assert res["status"] == "found"
        # ladder stopped early: not all probe points evaluated
        assert len(calls) < 2 * len(res["steps"])

    def test_probe_zero_plateau_not_found_and_inconclusive(self):
        def ev_ok(mu1):
            return {"success": True, "is_continuum": False,
                    "is_plateau": False}

        res = probe_zero_plateau(0.0, (0.0, 1.0), zero_tol=1e-10,
                                 probe_radius=1e-3, bracket_width=0.0,
                                 evaluator=ev_ok)
        assert res["status"] == "not_found"
        assert res["found"] is False

        def ev_left_only(mu1):
            # right-side probes FAIL (bisection error on the out-of-spectrum
            # side) — exactly the M2 scenario; the ladder must not raise
            return ({"success": False, "is_plateau": False}
                    if mu1 > 0 else
                    {"success": True, "is_continuum": False,
                     "is_plateau": False})

        res2 = probe_zero_plateau(0.0, None, zero_tol=1e-10,
                                  probe_radius=1e-3, bracket_width=0.0,
                                  evaluator=ev_left_only)
        assert res2["status"] == "inconclusive"


class TestClusteredOnTorus:
    def test_paired_points_cluster(self):
        pts = [(0.0, 0.0), (1e-4, 0.0), (3.0, 3.0), (3.0, 3.0 + 1e-4)]
        assert check_points_clustered_on_torus(pts, 1e-3) is True

    def test_isolated_point_fails(self):
        pts = [(0.0, 0.0), (3.0, 3.0)]
        assert check_points_clustered_on_torus(pts, 1e-3) is False

    def test_wraparound_neighbour_counts(self):
        # θ₁ = 0 and θ₁ = 2π − ε are torus-neighbours
        pts = [(0.0, 0.0), (TWO_PI - 1e-5, 0.0)]
        assert check_points_clustered_on_torus(pts, 1e-3) is True


class TestIsMrClusterEndpoint:
    def _mk(self, cluster=(1, 2), *, has_boundary_mr=True,
            with_perm=True, left_mr=0, right_mr=-1):
        from types import SimpleNamespace
        from bfgbz2d.continuation.multiple_roots import MultipleRootInfo
        mr = MultipleRootInfo(
            theta1=1.0, cluster_indices=[tuple(cluster)],
            roots=np.array([1.0 + 0j, 2.0 + 0j, 3.0 + 0j, 4.0 + 0j]),
        )
        seg = SimpleNamespace(left_mr=left_mr, right_mr=right_mr)
        zm = SimpleNamespace(
            multiple_roots=[mr], K=4, has_boundary_mr=has_boundary_mr)
        if with_perm:
            zm.boundary_perm = np.array([2, 0, 3, 1])
        return zm, seg

    def test_left_cluster_member_matches_by_column(self):
        zm, seg = self._mk((1, 2))
        assert is_mr_cluster_endpoint(zm, seg, 'left', 1) is True
        assert is_mr_cluster_endpoint(zm, seg, 'left', 2) is True

    def test_left_regular_column_passes(self):
        zm, seg = self._mk((1, 2))
        assert is_mr_cluster_endpoint(zm, seg, 'left', 0) is False
        assert is_mr_cluster_endpoint(zm, seg, 'left', 3) is False

    def test_no_mr_side_is_false(self):
        zm, seg = self._mk((1, 2), left_mr=-1)
        assert is_mr_cluster_endpoint(zm, seg, 'left', 1) is False

    def test_empty_cluster_is_false(self):
        zm, seg = self._mk(())
        assert is_mr_cluster_endpoint(zm, seg, 'left', 1) is False

    def test_boundary_mr_right_columns_translate(self):
        # Convention: roots_right[boundary_perm] == roots_left.  Left-frame
        # cluster columns (1, 2) are right-frame columns (0, 3).
        zm, seg = self._mk((1, 2), right_mr=0)
        assert is_mr_cluster_endpoint(zm, seg, 'right', 0) is True
        assert is_mr_cluster_endpoint(zm, seg, 'right', 3) is True
        assert is_mr_cluster_endpoint(zm, seg, 'right', 1) is False
        assert is_mr_cluster_endpoint(zm, seg, 'right', 2) is False

    def test_interior_mr_index_zero_right_is_not_translated(self):
        zm, seg = self._mk(
            (1, 2), has_boundary_mr=False, right_mr=0)
        assert is_mr_cluster_endpoint(zm, seg, 'right', 1) is True
        assert is_mr_cluster_endpoint(zm, seg, 'right', 2) is True
        assert is_mr_cluster_endpoint(zm, seg, 'right', 0) is False

    def test_boundary_mr_right_requires_perm(self):
        zm, seg = self._mk((1, 2), right_mr=0, with_perm=False)
        with pytest.raises(RuntimeError, match="boundary_perm"):
            is_mr_cluster_endpoint(zm, seg, 'right', 0)


class TestCharPolyMinorDegrees:
    def test_square_lattice_nn(self):
        coeffs = np.array([1, -1, -1, -1, -1], dtype=complex)
        degs = np.array([
            [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1],
        ], dtype=int)
        poly = CharPoly(coeffs, degs)
        assert poly.get_minor_degrees(2) == (1, 1)
        assert poly.get_minor_degrees(1) == (1, 1)
        assert get_minor_degrees(poly) == (1, 1)


class TestCharPolyRootPadding:
    """Degree-deficient partial polynomials must still return M+N roots."""

    @staticmethod
    def _roots(coeffs, degs, beta1=1.0 + 0j):
        poly = CharPoly(np.asarray(coeffs, dtype=complex),
                        np.asarray(degs, dtype=int))
        return poly, np.asarray(poly.solve_roots_1d(
            (0, 1), (0j, beta1), (2,)))

    def test_high_degree_deficiency_pads_all_infinite_roots(self):
        # β₁ β₂² + β₂ − 3 β₂⁻¹: M=1, N=2.  At β₁=0 both high-degree terms
        # implied by N=2 vanish except β₂; the finite polynomial degenerates
        # from degree 3 to degree 1, so two roots are at infinity.
        coeffs = [1, 1, -3]
        degs = [[0, 1, 2], [0, 0, 0], [0, 0, -1]]
        poly, roots = self._roots(coeffs, degs, beta1=0j)
        assert (poly.M, poly.N) == (1, 2)
        assert len(roots) == poly.M + poly.N
        assert np.count_nonzero(np.isinf(roots)) == 2
        assert np.count_nonzero(np.isfinite(roots)) == 1

    def test_low_degree_deficiency_pads_all_zero_roots(self):
        # β₁ β₂⁻² + β₂ − 3: M=2, N=1.  At β₁=0 the two denominator roots
        # vanish with the β₂⁻² term, leaving one finite root and two β₂=0
        # padding roots.
        coeffs = [1, 1, -3]
        degs = [[0, 1, -2], [0, 0, 0], [0, 0, 1]]
        poly, roots = self._roots(coeffs, degs, beta1=0j)
        assert (poly.M, poly.N) == (2, 1)
        assert len(roots) == poly.M + poly.N
        assert np.count_nonzero(roots == 0) == 2
        assert np.count_nonzero(np.isinf(roots)) == 0
        assert np.count_nonzero((roots != 0) & np.isfinite(roots)) == 1

    def test_both_sides_deficient_pads_multiple_zero_and_infinite(self):
        # β₁ β₂³ + 1 + β₁ β₂⁻³: at β₁=0 both the three denominator roots and
        # three numerator roots disappear.  The old code appended only one 0
        # and one ∞, returning 2 roots instead of K=6.
        coeffs = [1, 1, 1]
        degs = [[0, 1, 3], [0, 0, 0], [0, 1, -3]]
        poly, roots = self._roots(coeffs, degs, beta1=0j)
        assert (poly.M, poly.N) == (3, 3)
        assert len(roots) == poly.M + poly.N
        assert np.count_nonzero(roots == 0) == 3
        assert np.count_nonzero(np.isinf(roots)) == 3

    def test_completely_vanished_partial_polynomial(self):
        # β₁(β₂³ + β₂⁻³): at β₁=0 the partial polynomial has no terms at
        # all.  It must still return the fixed K=6 track array, not crash on
        # max([]), and the padding is M zeros + N infinities.
        coeffs = [1, 1]
        degs = [[0, 1, 3], [0, 1, -3]]
        poly, roots = self._roots(coeffs, degs, beta1=0j)
        assert (poly.M, poly.N) == (3, 3)
        assert len(roots) == poly.M + poly.N
        assert np.count_nonzero(roots == 0) == 3
        assert np.count_nonzero(np.isinf(roots)) == 3

    def test_undepleted_polynomial_is_not_padded(self):
        coeffs = [1, 1, 1]
        degs = [[0, 1, 3], [0, 0, 0], [0, 1, -3]]
        poly, roots = self._roots(coeffs, degs, beta1=1.0 + 0j)
        assert len(roots) == poly.M + poly.N
        assert np.all(np.isfinite(roots))
