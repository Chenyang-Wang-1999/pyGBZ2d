"""Tests for brute_force_amoeba module using 2D HN model analytic solution."""

import numpy as np
import pytest
from cmath import exp
import poly_tools as pt

import brute_force_amoeba as bfa
from brute_force_amoeba.bisect import _resolve_continuum
from brute_force_amoeba.tracks import _compute_root_tracks
from gbz_types import PointSubset, LineSubset, GBZResult


def build_HN2D_polynomial(J1, J2, gamma_1, gamma_2, delta_1, delta_2, basis="10"):
    J11 = exp(gamma_1 + 1j * delta_1) * J1
    J12 = exp(-gamma_1 + 1j * delta_1) * np.conj(J1)
    J21 = exp(gamma_2 + 1j * delta_2) * J2
    J22 = exp(-gamma_2 + 1j * delta_2) * np.conj(J2)

    coeffs = np.array([1, -J11, -J12, -J21, -J22], dtype=complex)
    degs = np.array([
        [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1],
    ], dtype=int)
    return coeffs, degs


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
        # Amoeba check may return non-empty but is_amoeba check should fail
        # Result: empty GBZResult
        assert not gbz.is_gbz or gbz.is_empty
        # At minimum, outside E should have index (0,0)
        if not gbz.is_gbz:
            assert gbz.index == (0, 0)

    def test_beta_magnitudes_analytic(self, poly_A, params_A):
        coeffs, degs = poly_A
        gbz = bfa.collect_GBZ_subsets(coeffs, degs, 0.0 + 0j, 0.0)
        if not gbz.is_gbz:
            pytest.skip("Amoeba returned empty for E=0")
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
        "get_hungarian_sorted_roots",
        "get_a2_average_winding", "get_a1_average_winding",
        "bisect_amoeba_ronkin_min",
        "collect_GBZ_subsets",
        "PointSubset", "LineSubset", "GBZResult", "ConnectedSubset",
    ]
    for name in expected:
        assert hasattr(bfa, name), f"Missing export: {name}"


# ---- _resolve_continuum tests ----


def _build_char_poly(coeffs, degs):
    """Build a CLaurent polynomial from coefficient / degree arrays."""
    cp = pt.CLaurent(3)
    coeffs_ct = pt.CScalarVec(coeffs)
    degs_ct = pt.CLaurentIndexVec(degs.flatten())
    cp.set_Laurent_by_terms(coeffs_ct, degs_ct)
    return cp


class TestResolveContinuum:
    """Tests for _resolve_continuum using 2D HN model.

    HN model: J1=J2=1, gamma_1=0.2, gamma_2=0.3.
    Analytic SGBZ / amoeba GBZ at mu1=0.2, mu2=0.3.
    """

    @pytest.fixture
    def hn_char_poly(self, poly_A):
        coeffs, degs = poly_A
        return _build_char_poly(coeffs, degs)

    def test_at_known_boundary(self, hn_char_poly):
        """At (mu1=0.2, mu2=0.3) with E=0, the amoeba boundary is a
        continuum.  Both w2 and w1 limits should straddle zero."""
        E_ref = 0.0 + 0j
        mu1 = 0.2
        mu2 = 0.3

        tracks = _compute_root_tracks(hn_char_poly, E_ref, mu1, N_points=301)
        result = _resolve_continuum(
            hn_char_poly, E_ref, mu1, mu2, tracks,
            continuum_perturb=1e-4,
        )

        # w2 limits should exist and straddle zero
        assert result["w2_left"] is not None
        assert result["w2_right"] is not None
        assert result["w2_opposite"] is True, (
            f"w2_left={result['w2_left']}, w2_right={result['w2_right']}"
        )

        # w1 should have been resolved and also straddle zero
        assert result["w1_resolved"] is True
        assert result["w1_left"] is not None
        assert result["w1_right"] is not None
        assert result["w1_opposite"] is True, (
            f"w1_left={result['w1_left']}, w1_right={result['w1_right']}"
        )

        assert result["is_boundary"] is True

    def test_outside_spectrum_no_crash(self, hn_char_poly):
        """E=5 is outside the HN spectrum.  The function should not crash
        and should return a well-formed dict regardless of the physics."""
        E_ref = 5.0 + 0j
        mu1 = 0.2
        mu2 = 0.3

        tracks = _compute_root_tracks(hn_char_poly, E_ref, mu1, N_points=301)
        result = _resolve_continuum(
            hn_char_poly, E_ref, mu1, mu2, tracks,
            continuum_perturb=1e-4,
        )

        # Verify well-formed return dict
        for key in ("w2_left", "w2_right", "w1_left", "w1_right"):
            assert key in result
        assert isinstance(result["w2_opposite"], bool)
        assert isinstance(result["is_boundary"], bool)
        # w1_opposite and w1_resolved may be None when w2 not opposite
        for key in ("w1_opposite", "w1_resolved"):
            assert key in result

        # w2 limits should be valid (not None) — even without continuum,
        # _compute_winding_from_tracks returns normal winding values
        assert result["w2_left"] is not None
        assert result["w2_right"] is not None

    def test_w2_not_opposite_no_w1(self, hn_char_poly):
        """When w2 limits don't straddle 0, w1 should not be resolved.

        At mu2=0 (far below the Ronkin minimum at mu2=0.3), the w2
        winding is on the same side of zero at mu2 ± ε, so the limits
        don't straddle 0 and w1 resolution is skipped."""
        E_ref = 0.0 + 0j
        mu1 = 0.2
        mu2 = 0.0

        tracks = _compute_root_tracks(hn_char_poly, E_ref, mu1, N_points=301)
        result = _resolve_continuum(
            hn_char_poly, E_ref, mu1, mu2, tracks,
            continuum_perturb=1e-4,
        )

        assert result["w2_opposite"] is False
        assert result["w1_resolved"] is False
        assert result["w1_left"] is None
        assert result["w1_right"] is None
        assert result["w1_opposite"] is None
        assert result["is_boundary"] is False


# ---- plateau edge detection tests ----
# Reference case from diagnostics/plateau-check-commit-report.md:
# the next-nearest-coupling model at E = 1.648 + 0.0294j exposes a
# zero-plateau edge where the bisection finds 4 clustered zeros with
# canceling jumps (net_zero_count = 0).  The plateau pre-check must
# correctly trigger the probe and classify the point as outside the
# amoeba spectrum.


@pytest.fixture
def nnc_char_poly():
    """Build the next-nearest-coupling model from demos/imaginary-degeneracy-splitting.py."""
    from BerryPy import TightBinding as tb

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
    return _build_char_poly(coeffs, degs), coeffs, degs


class TestPlateauEdge:
    """Verify the plateau pre-check correctly handles a known zero-plateau edge."""

    E_PLATEAU_EDGE = 1.648 + 0.0294j

    def test_plateau_edge_outside_spectrum(self, nnc_char_poly):
        """At the plateau-edge E, should return is_gbz=False with empty subsets."""
        char_poly, coeffs, degs = nnc_char_poly
        gbz = bfa.collect_GBZ_subsets(
            coeffs, degs, self.E_PLATEAU_EDGE, 0.0, plateau_check=True,
        )
        assert not gbz.is_gbz, (
            f"Plateau edge should be outside spectrum, got is_gbz={gbz.is_gbz}"
        )
        assert gbz.index == (0, 0)

    def test_plateau_pre_check_triggers_probe(self, nnc_char_poly):
        """The plateau pre-check should trigger the probe: zeros are clustered
        with tiny w1_area at the plateau edge."""
        char_poly, coeffs, degs = nnc_char_poly
        from brute_force_amoeba.bisect import bisect_amoeba_ronkin_min
        from brute_force_amoeba.ronkin_winding import _get_average_winding_from_zeros
        from brute_force_amoeba.amoeba import _check_zeros_are_clustered

        res = bisect_amoeba_ronkin_min(char_poly, self.E_PLATEAU_EDGE, N_points=301)
        zeros = res["zeros"]

        # Condition (a): non-zero winding area must be tiny
        w1_area = res["_w1_area"]
        _, w2_area = _get_average_winding_from_zeros(
            char_poly, self.E_PLATEAU_EDGE, res["mu1"], res["mu2"],
            zeros, direction=2,
        )
        plateau_area_threshold = 1e-2
        assert w1_area < plateau_area_threshold, f"w1_area={w1_area} should be tiny at plateau edge"
        assert w2_area < plateau_area_threshold, f"w2_area={w2_area} should be tiny at plateau edge"

        # Condition (b): zeros must be clustered
        assert len(zeros) > 0, "Should find zeros at plateau edge"
        assert _check_zeros_are_clustered(zeros, plateau_area_threshold), (
            "Zeros should be clustered at plateau edge"
        )

    def test_plateau_edge_net_zero_count(self, nnc_char_poly):
        """Zeros at the plateau edge should have canceling jump directions
        (net_zero_count = 0), distinguishing them from genuine GBZ zeros."""
        char_poly, coeffs, degs = nnc_char_poly
        from brute_force_amoeba.bisect import bisect_amoeba_ronkin_min

        res = bisect_amoeba_ronkin_min(char_poly, self.E_PLATEAU_EDGE, N_points=301)
        zeros = res["zeros"]

        # Sum of jump directions should be 0 (zeros cancel in pairs)
        jump_sum = sum(z[2] for z in zeros)
        assert jump_sum == 0, (
            f"Plateau-edge zeros should have canceling jumps, got sum={jump_sum}"
        )
