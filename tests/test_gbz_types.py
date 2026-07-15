"""Unit tests for gbz_types data classes and utility functions."""

import math
import numpy as np
import pytest
from cmath import exp

from gbz_types import (
    PointSubset, LineSubset, GBZResult, ConnectedSubset,
    sort_by_root_abs, to_sphere_r3, chordal_cost_matrix,
    hungarian_match_indices, find_cyclic_true_intervals,
    generate_probe_steps,
)


class TestPointSubset:
    def test_creation(self):
        ps = PointSubset(E=1.0 + 0j, beta1=exp(0.2 + 1j * 0.5), beta2=exp(0.3 + 1j * 1.0))
        assert ps.E == 1.0 + 0j
        assert abs(ps.beta1) == pytest.approx(exp(0.2))
        assert abs(ps.beta2) == pytest.approx(exp(0.3))

    def test_mu1_property(self):
        ps = PointSubset(E=0j, beta1=exp(0.5 + 1j * 1.0), beta2=1.0 + 0j)
        assert ps.mu1 == pytest.approx(0.5)

    def test_theta1_property(self):
        ps = PointSubset(E=0j, beta1=exp(1j * math.pi / 2), beta2=1.0 + 0j)
        assert ps.theta1 == pytest.approx(math.pi / 2)

    def test_frozen(self):
        ps = PointSubset(E=0j, beta1=1j, beta2=1j)
        with pytest.raises(Exception):
            ps.beta1 = 2j  # type: ignore

    def test_as_triplet(self):
        ps = PointSubset(E=1.0 + 0.5j, beta1=exp(0.2 + 1j * 1.0), beta2=exp(0.3 + 1j * 2.0))
        E, k1, k2 = ps.as_triplet()
        assert E == 1.0 + 0.5j
        assert abs(k1.imag + 0.2) < 1e-10  # Im(k1) = -mu1
        assert abs(k2.imag + 0.3) < 1e-10  # Im(k2) = -mu2


class TestLineSubset:
    def test_creation(self):
        ls = LineSubset(E=1.0 + 0j, mu1=0.2, theta1_start=0.0, theta1_end=math.pi)
        assert ls.E == 1.0 + 0j
        assert ls.mu1 == 0.2
        assert ls.beta2_arr is None
        assert not ls.is_loaded()

    def test_theta1_width_normal(self):
        ls = LineSubset(E=0j, mu1=0.0, theta1_start=0.5, theta1_end=2.0)
        assert ls.theta1_width == pytest.approx(1.5)

    def test_theta1_width_wraparound(self):
        ls = LineSubset(E=0j, mu1=0.0, theta1_start=6.0, theta1_end=0.5)
        assert ls.theta1_width == pytest.approx(2 * math.pi - 5.5)

    def test_endpoints_require_fill(self):
        ls = LineSubset(E=0j, mu1=0.0, theta1_start=0.0, theta1_end=1.0)
        with pytest.raises(RuntimeError, match="fill_beta2"):
            _ = ls.left_endpoint


class TestGBZResult:
    def test_empty(self):
        gbz = GBZResult(E_ref=1.0 + 0j, subsets=[], index=(0, 0))
        assert gbz.is_empty
        assert not gbz.is_gbz
        assert gbz.E_ref == 1.0 + 0j

    def test_nonempty_point(self):
        ps = PointSubset(E=1.0 + 0j, beta1=1j, beta2=2j)
        gbz = GBZResult(E_ref=1.0 + 0j, subsets=[ps], index=(1, 0))
        assert not gbz.is_empty
        assert gbz.is_gbz
        assert gbz.index == (1, 0)

    def test_nonempty_line(self):
        ls = LineSubset(E=1.0 + 0j, mu1=0.2, theta1_start=0.0, theta1_end=1.0)
        gbz = GBZResult(E_ref=1.0 + 0j, subsets=[ls], index=(0, 1))
        assert not gbz.is_empty
        assert gbz.is_gbz
        assert gbz.index == (0, 1)

    def test_mixed(self):
        ps = PointSubset(E=0j, beta1=1j, beta2=1j)
        ls = LineSubset(E=0j, mu1=0.0, theta1_start=0.0, theta1_end=1.0)
        gbz = GBZResult(E_ref=0j, subsets=[ps, ls], index=(1, 1))
        assert gbz.index == (1, 1)
        assert len(gbz.subsets) == 2


class TestSortByRootAbs:
    def test_sorting(self):
        roots = np.array([2 + 0j, 1 + 0j, 3 + 0j, 1 + 1j])
        sorted_roots = sort_by_root_abs(roots)
        assert abs(sorted_roots[0]) == pytest.approx(1.0)
        assert abs(sorted_roots[-1]) == pytest.approx(3.0)

    def test_infinity_last(self):
        roots = np.array([2 + 0j, np.inf, 1 + 0j])
        sorted_roots = sort_by_root_abs(roots)
        assert np.isinf(abs(sorted_roots[-1]))


class TestToSphereR3:
    def test_finite_point(self):
        r3 = to_sphere_r3(np.array([1 + 0j]))
        assert r3.shape == (1, 3)
        # z=1 maps to (0, 0, 0)? No: stereographic from south pole
        # For z=1: x=1,y=0 → r3 = (2*1/2, 0, (1-1)/2) = (1, 0, 0)
        assert r3[0, 0] == pytest.approx(1.0)
        assert r3[0, 1] == pytest.approx(0.0)
        assert r3[0, 2] == pytest.approx(0.0)

    def test_origin(self):
        r3 = to_sphere_r3(np.array([0 + 0j]))
        assert r3[0, 0] == pytest.approx(0.0)
        assert r3[0, 2] == pytest.approx(-1.0)

    def test_infinity(self):
        r3 = to_sphere_r3(np.array([np.inf]))
        assert r3[0, 2] == pytest.approx(1.0)


class TestChordalCostMatrix:
    def test_identical_arrays(self):
        arr = np.array([1 + 0j, 2 + 0j])
        cost = chordal_cost_matrix(arr, arr)
        assert cost[0, 0] == pytest.approx(0.0, abs=1e-10)
        assert cost[1, 1] == pytest.approx(0.0, abs=1e-10)

    def test_shape(self):
        a = np.array([1 + 0j, 2 + 0j, 3 + 0j])
        b = np.array([1j, 2j])
        cost = chordal_cost_matrix(a, b)
        assert cost.shape == (3, 2)


class TestHungarianMatchIndices:
    def test_simple_match(self):
        a = np.array([1 + 0j, 2 + 0j])
        b = np.array([2 + 0j, 1 + 0j])
        matches = hungarian_match_indices(a, b)
        assert len(matches) == 2
        # a[0]=1 should match b[1]=1, a[1]=2 should match b[0]=2
        match_dict = dict(matches)
        assert match_dict[0] == 1
        assert match_dict[1] == 0

    def test_nan_raises(self):
        a = np.array([np.nan, 1 + 0j])
        b = np.array([1 + 0j, 2 + 0j])
        with pytest.raises(ValueError, match="NaN"):
            hungarian_match_indices(a, b)


class TestFindCyclicTrueIntervals:
    def test_all_true(self):
        mask = np.array([True, True, True])
        intervals = find_cyclic_true_intervals(mask)
        assert intervals == [(0, 2)]

    def test_all_false(self):
        mask = np.array([False, False, False])
        intervals = find_cyclic_true_intervals(mask)
        assert intervals == []

    def test_single_interval(self):
        mask = np.array([False, True, True, False])
        intervals = find_cyclic_true_intervals(mask)
        assert intervals == [(1, 2)]

    def test_cross_boundary(self):
        mask = np.array([True, True, False, False, True])
        intervals = find_cyclic_true_intervals(mask)
        # One interval wrapping around: indices 4, 0, 1 are True
        assert len(intervals) == 1
        assert intervals[0] == (4, 1)

    def test_empty(self):
        assert find_cyclic_true_intervals(np.array([], dtype=bool)) == []


class TestGenerateProbeSteps:
    def test_sorted_positive(self):
        steps = generate_probe_steps(1.0, 0.5, 1e-10)
        assert steps == sorted(steps)
        assert all(s > 0 for s in steps)

    def test_contains_expected(self):
        steps = generate_probe_steps(1.0, 0.5, 1e-10)
        assert 0.25 in steps or any(abs(s - 0.25) < 1e-15 for s in steps)
        assert 0.5 in steps or any(abs(s - 0.5) < 1e-15 for s in steps)
        assert 1.0 in steps or any(abs(s - 1.0) < 1e-15 for s in steps)
