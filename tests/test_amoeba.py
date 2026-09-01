"""Tests for pygbz2d.amoeba module using 2D HN model analytic solution."""

import numpy as np
import pytest
from cmath import exp

import pygbz2d.amoeba as bfa
from pygbz2d.amoeba.bisect import (
    bisect_amoeba_ronkin_min,
    _find_mu2_for_w2_zero,
    _try_fast_mu2,
)
from pygbz2d.amoeba.zm_extract import (
    AmoebaZeroManager,
    calculate_a2_average_winding,
    detect_continuum,
    find_crossings,
)
from pygbz2d.core import PointSubset, LineSubset, GBZResult, CharPoly, TWO_PI

from conftest import build_HN2D_polynomial


@pytest.fixture
def params_A():
    return {"J1": 1.0, "J2": 1.0, "gamma_1": 0.2, "gamma_2": 0.3,
            "delta_1": 0.0, "delta_2": 0.0}


@pytest.fixture
def poly_A(params_A):
    return build_HN2D_polynomial(**params_A)


class TestAmoeba:
    def test_returns_gbzresult(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfa.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0)
        assert isinstance(gbz, GBZResult)

    def test_inside_spectrum(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfa.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0)
        assert gbz.is_gbz
        assert not gbz.is_empty
        assert gbz.index != (0, 0)

    def test_outside_spectrum(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfa.collect_GBZ_subsets(coeffs, degs, 5.0 + 0j, 0.0)
        assert not gbz.is_gbz or gbz.is_empty
        if not gbz.is_gbz:
            assert gbz.index == (0, 0)

    def test_beta_magnitudes_analytic(self, poly_A, params_A):
        """At E=0 the HN amoeba GBZ is a full-circle continuum."""
        coeffs, degs = poly_A
        gbz = bfa.collect_GBZ_subsets(coeffs, degs, 0.0 + 0j, 0.0)
        assert gbz.is_gbz, (
            f"E=0 full-circle continuum should be GBZ, got index={gbz.index}"
        )
        assert gbz.index[0] == 0, (
            f"E=0 should have no discrete points, got {gbz.index[0]}"
        )
        assert gbz.index[1] > 0, (
            f"E=0 should have continuum LineSubsets, got {gbz.index[1]}"
        )
        gamma_1, gamma_2 = params_A["gamma_1"], params_A["gamma_2"]
        for s in gbz.subsets:
            if isinstance(s, PointSubset):
                assert abs(s.beta1) == pytest.approx(exp(gamma_1), rel=1e-3)
                assert abs(s.beta2) == pytest.approx(exp(gamma_2), rel=1e-3)
            elif isinstance(s, LineSubset):
                assert s.mu1 == pytest.approx(gamma_1, rel=1e-3)

    def test_subsets_valid_types(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfa.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0)
        for s in gbz.subsets:
            assert isinstance(s, (PointSubset, LineSubset))

    def test_index_consistent(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfa.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0)
        n_0d = sum(1 for s in gbz.subsets if isinstance(s, PointSubset))
        n_1d = sum(1 for s in gbz.subsets if isinstance(s, LineSubset))
        assert gbz.index == (n_0d, n_1d)


def test_all_exports():
    expected = [
        "bisect_amoeba_ronkin_min",
        "AmoebaZeroManager",
        "calculate_a2_average_winding",
        "collect_GBZ_subsets",
        "PointSubset", "LineSubset", "GBZResult", "ConnectedSubset",
    ]
    for name in expected:
        assert hasattr(bfa, name), f"Missing export: {name}"


# ---- crossing / winding primitives ----

class TestCrossingPrimitives:
    @pytest.fixture
    def zm_hn(self, poly_A):
        coeffs, degs = poly_A
        poly = CharPoly(coeffs, degs)
        zm = AmoebaZeroManager(poly, 1.0 + 0j, 0.0)
        zm.run()
        return zm, poly

    def test_find_crossings_returns_beta_pairs(self, zm_hn):
        zm, poly = zm_hn
        crossings = find_crossings(zm, 0.0, 0.3, return_refined=True)
        assert isinstance(crossings, list)
        assert len(crossings) == 2
        for b1, b2 in crossings:
            assert abs(b1) == pytest.approx(1.0, rel=1e-9)
            assert abs(b2) == pytest.approx(exp(0.3), rel=1e-9)

    def test_find_crossings_avoided_segments(self, zm_hn):
        zm, poly = zm_hn
        crossings = find_crossings(
            zm, 0.0, 0.3,
            avoided_segments=[(0, 0), (0, 1)],
            return_refined=True,
        )
        assert crossings == []

    def test_calculate_a2_average_winding(self, zm_hn):
        zm, poly = zm_hn
        w2 = calculate_a2_average_winding(zm, 0.0, 0.3)
        assert abs(w2) < 1e-9

    def test_try_fast_mu2_has_expected_keys(self, zm_hn):
        zm, poly = zm_hn
        fast = _try_fast_mu2(zm, poly)
        assert set(["ok", "A", "B"]) <= set(fast.keys())


class TestBisection:
    def test_known_continuum_boundary(self, poly_A):
        coeffs, degs = poly_A
        poly = CharPoly(coeffs, degs)
        res = bisect_amoeba_ronkin_min(poly, 1.0 + 0j)
        assert res["is_continuum"] is True
        assert res["mu1"] == pytest.approx(0.2, abs=1e-3)
        assert res["mu2"] == pytest.approx(0.3, abs=1e-3)
        assert "_continuum_members" in res

    def test_inner_fast_path_discrete(self, poly_A):
        coeffs, degs = poly_A
        poly = CharPoly(coeffs, degs)
        zm = AmoebaZeroManager(poly, 1.0 + 0j, -0.5)
        zm.run()
        inner = _find_mu2_for_w2_zero(poly, 1.0 + 0j, -0.5, -1.0, 1.0, _zm=zm)
        assert inner["is_continuum"] is False
        assert inner["mu2"] == pytest.approx(0.3, abs=1e-6)


# ---- plateau edge detection tests ----
# Reference case from diagnostics/plateau-check-commit-report.md:
# the next-nearest-coupling model at E = 1.648 + 0.0294j exposes a
# zero-plateau edge where the bisection finds no discrete crossings.


@pytest.fixture
def nnc_char_poly():
    """Build the next-nearest-coupling model."""
    tb = pytest.importorskip("BerryPy").TightBinding

    gamma = 0.2
    u1, v1, w = 1.0, 0.8, 0.5
    u2, v2 = 0.1, 0.1
    u3, v3 = 0.05, 0.05
    inter_cell = [
        [0, 0, u1 * np.exp(gamma), (1, 0)],
        [0, 0, u1 * np.exp(-gamma), (-1, 0)],
        [0, 0, u2, (2, 0)],
        [0, 0, u2, (-2, 0)],
        [0, 0, u3, (3, 0)],
        [0, 0, u3, (-3, 0)],
        [0, 0, v1, (0, 1)],
        [0, 0, v1, (0, -1)],
        [0, 0, v2, (0, 2)],
        [0, 0, v2, (0, -2)],
        [0, 0, v3, (0, 3)],
        [0, 0, v3, (0, -3)],
        [0, 0, w, (1, 1)],
        [0, 0, w, (1, -1)],
        [0, 0, w, (-1, 1)],
        [0, 0, w, (-1, -1)],
    ]
    model = tb.TightBindingModel(2, 1, np.eye(2), [], inter_cell, [(0, 0)])
    coeffs, degs = model.get_characteristic_polynomial_data()
    return CharPoly(coeffs, degs), coeffs, degs


class TestPlateauEdge:
    """Verify the plateau pre-check correctly handles a known zero-plateau edge."""

    E_PLATEAU_EDGE = 1.648 + 0.0294j

    def test_plateau_edge_outside_spectrum(self, nnc_char_poly):
        char_poly, coeffs, degs = nnc_char_poly
        gbz = bfa.collect_GBZ_subsets(
            coeffs, degs, self.E_PLATEAU_EDGE, 0.0, plateau_check=True,
        )
        assert not gbz.is_gbz, (
            f"Plateau edge should be outside spectrum, got is_gbz={gbz.is_gbz}"
        )
        assert gbz.index == (0, 0)

    def test_plateau_edge_no_discrete_crossings(self, nnc_char_poly):
        char_poly, coeffs, degs = nnc_char_poly
        res = bisect_amoeba_ronkin_min(char_poly, self.E_PLATEAU_EDGE)
        assert res["is_continuum"] is False
        assert len(res["zeros"]) == 0
