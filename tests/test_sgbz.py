"""Tests for brute_force_SGBZ module using 2D HN model analytic solution."""

import numpy as np
import pytest
from cmath import exp, log
import poly_tools as pt

import brute_force_SGBZ as bfs
from gbz_types import PointSubset, LineSubset, GBZResult


# ---- shared helpers ----

def build_HN2D_polynomial(J1, J2, gamma_1, gamma_2, delta_1, delta_2, basis="10"):
    """Build 2D HN model characteristic polynomial."""
    J11 = exp(gamma_1 + 1j * delta_1) * J1
    J12 = exp(-gamma_1 + 1j * delta_1) * np.conj(J1)
    J21 = exp(gamma_2 + 1j * delta_2) * J2
    J22 = exp(-gamma_2 + 1j * delta_2) * np.conj(J2)

    coeffs = np.array([1, -J11, -J12, -J21, -J22], dtype=complex)

    if basis == "10":
        degs = np.array([
            [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1],
        ], dtype=int)
    elif basis == "11":
        degs = np.array([
            [1, 0, 0], [0, -1, 1], [0, 1, -1], [0, 0, -1], [0, 0, 1],
        ], dtype=int)
    else:
        raise ValueError(f"Unknown basis: {basis}")

    return coeffs, degs


# ---- fixtures ----

@pytest.fixture
def params_A():
    return {"J1": 1.0, "J2": 1.0, "gamma_1": 0.2, "gamma_2": 0.3,
            "delta_1": 0.0, "delta_2": 0.0}


@pytest.fixture
def poly_A(params_A):
    return build_HN2D_polynomial(**params_A, basis="10")


@pytest.fixture
def poly_A_11(params_A):
    return build_HN2D_polynomial(**params_A, basis="11")


# ---- [10]-SGBZ tests ----

class TestSGBZ10:
    """Tests for [10]-SGBZ of 2D HN model."""

    def test_returns_gbzresult(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        assert isinstance(gbz, GBZResult)
        assert gbz.E_ref == 1.0 + 0j

    def test_inside_spectrum(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        assert gbz.is_gbz
        assert not gbz.is_empty
        assert gbz.index != (0, 0)

    def test_outside_spectrum(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfs.check_SGBZ(coeffs, degs, 5.0 + 0j, 0.0, N_points=101)
        assert not gbz.is_gbz
        assert gbz.is_empty
        assert gbz.index == (0, 0)

    def test_mu1_matches_analytic(self, poly_A, params_A):
        coeffs, degs = poly_A
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        assert gbz.is_gbz
        for s in gbz.subsets:
            if isinstance(s, PointSubset):
                assert s.mu1 == pytest.approx(params_A["gamma_1"], rel=1e-5)
            elif isinstance(s, LineSubset):
                assert s.mu1 == pytest.approx(params_A["gamma_1"], rel=1e-5)

    def test_beta_magnitudes_analytic(self, poly_A, params_A):
        coeffs, degs = poly_A
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        gamma_1, gamma_2 = params_A["gamma_1"], params_A["gamma_2"]
        for s in gbz.subsets:
            if isinstance(s, PointSubset):
                assert abs(s.beta1) == pytest.approx(exp(gamma_1), rel=1e-4)
                assert abs(s.beta2) == pytest.approx(exp(gamma_2), rel=1e-4)
            elif isinstance(s, LineSubset):
                b1 = exp(s.mu1 + 1j * s.theta1_start)
                assert abs(b1) == pytest.approx(exp(gamma_1), rel=1e-4)

    def test_subsets_valid_types(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        for s in gbz.subsets:
            assert isinstance(s, (PointSubset, LineSubset))

    def test_index_consistent(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        n_0d = sum(1 for s in gbz.subsets if isinstance(s, PointSubset))
        n_1d = sum(1 for s in gbz.subsets if isinstance(s, LineSubset))
        assert gbz.index == (n_0d, n_1d)


# ---- [11]-SGBZ tests ----

class TestSGBZ11:
    """Tests for [11]-SGBZ of 2D HN model."""

    def test_returns_gbzresult(self, poly_A_11):
        coeffs, degs = poly_A_11
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        assert isinstance(gbz, GBZResult)

    def test_inside_spectrum(self, poly_A_11):
        coeffs, degs = poly_A_11
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        assert gbz.is_gbz
        assert not gbz.is_empty

    def test_outside_spectrum(self, poly_A_11):
        coeffs, degs = poly_A_11
        gbz = bfs.check_SGBZ(coeffs, degs, 5.0 + 0j, 0.0, N_points=101)
        assert not gbz.is_gbz

    def test_mu1_matches_analytic(self, poly_A_11, params_A):
        coeffs, degs = poly_A_11
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        assert gbz.is_gbz
        expected_mu = params_A["gamma_1"] + params_A["gamma_2"]  # = 0.5
        for s in gbz.subsets:
            if isinstance(s, PointSubset):
                assert s.mu1 == pytest.approx(expected_mu, rel=1e-5)
            elif isinstance(s, LineSubset):
                assert s.mu1 == pytest.approx(expected_mu, rel=1e-5)

    def test_beta1_magnitude_analytic(self, poly_A_11, params_A):
        coeffs, degs = poly_A_11
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        expected_r = exp(params_A["gamma_1"] + params_A["gamma_2"])
        for s in gbz.subsets:
            if isinstance(s, PointSubset):
                assert abs(s.beta1) == pytest.approx(expected_r, rel=1e-4)


# ---- triplet conversion tests ----

class TestTripletConversion:
    def test_convert_returns_list(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        triplets = bfs.convert_gbz_to_triplets(gbz)
        assert isinstance(triplets, list)
        assert len(triplets) > 0

    def test_triplet_structure(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        triplets = bfs.convert_gbz_to_triplets(gbz)
        for t in triplets:
            assert len(t) == 3
            E, k1, k2 = t
            assert isinstance(E, complex)
            assert isinstance(k1, complex)
            assert isinstance(k2, complex)
            # Im(k1) = -mu1 (-gamma_1 for [10]-SGBZ)
            assert k1.imag == pytest.approx(-0.2, rel=1e-4)

    def test_convert_empty(self):
        gbz = GBZResult(E_ref=0j, subsets=[], index=(0, 0))
        triplets = bfs.convert_gbz_to_triplets(gbz)
        assert triplets == []

    def test_convert_list(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        triplets = bfs.convert_gbz_list_to_triplets([gbz])
        assert len(triplets) > 0


# ---- Hermitian limit tests (γ = 0 → GBZ = BZ) ----

class TestHermitianLimit:
    """Hermitian system: GBZ = BZ for all SGBZ directions and amoeba."""

    @pytest.fixture
    def poly_hermitian(self):
        return build_HN2D_polynomial(
            J1=1.0, J2=1.0,
            gamma_1=0.0, gamma_2=0.0,  # Hermitian!
            delta_1=0.0, delta_2=0.0,
            basis="10",
        )

    @pytest.fixture
    def poly_hermitian_11(self):
        return build_HN2D_polynomial(
            J1=1.0, J2=1.0,
            gamma_1=0.0, gamma_2=0.0,
            delta_1=0.0, delta_2=0.0,
            basis="11",
        )

    def test_sgbz10_mu1_zero(self, poly_hermitian):
        coeffs, degs = poly_hermitian
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        assert gbz.is_gbz
        for s in gbz.subsets:
            mu = s.mu1 if isinstance(s, PointSubset) else s.mu1
            assert mu == pytest.approx(0.0, abs=1e-5)

    def test_sgbz10_beta_on_unit_circle(self, poly_hermitian):
        coeffs, degs = poly_hermitian
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        for s in gbz.subsets:
            if isinstance(s, PointSubset):
                assert abs(s.beta1) == pytest.approx(1.0, rel=1e-4)
                assert abs(s.beta2) == pytest.approx(1.0, rel=1e-4)

    def test_sgbz11_mu1_zero(self, poly_hermitian_11):
        coeffs, degs = poly_hermitian_11
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        assert gbz.is_gbz
        for s in gbz.subsets:
            mu = s.mu1 if isinstance(s, PointSubset) else s.mu1
            assert mu == pytest.approx(0.0, abs=1e-5)

    def test_sgbz11_beta_on_unit_circle(self, poly_hermitian_11):
        coeffs, degs = poly_hermitian_11
        gbz = bfs.check_SGBZ(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        for s in gbz.subsets:
            if isinstance(s, PointSubset):
                assert abs(s.beta1) == pytest.approx(1.0, rel=1e-4)

    def test_hermitian_all_equal(self, poly_hermitian, poly_hermitian_11):
        """Hermitian limit: [10]-SGBZ, [11]-SGBZ, Amoeba all give same GBZ (=BZ)."""
        coeffs10, degs10 = poly_hermitian
        coeffs11, degs11 = poly_hermitian_11

        import brute_force_amoeba as bfa

        gbz10 = bfs.check_SGBZ(coeffs10, degs10, 1.0 + 0j, 0.0, N_points=201)
        gbz11 = bfs.check_SGBZ(coeffs11, degs11, 1.0 + 0j, 0.0, N_points=201)
        gbz_amoeba = bfa.collect_GBZ_subsets(coeffs10, degs10, 1.0 + 0j, 0.0)

        assert gbz10.is_gbz
        assert gbz11.is_gbz
        assert gbz_amoeba.is_gbz

        # All should have β on the unit circle
        for gbz in [gbz10, gbz11, gbz_amoeba]:
            for s in gbz.subsets:
                if isinstance(s, PointSubset):
                    assert abs(s.beta1) == pytest.approx(1.0, rel=1e-4)


# ---- export verification ----

def test_all_exports():
    """Verify all expected symbols are exported."""
    expected = [
        "ComplexEqConverter", "complex_root", "poly_to_np_coefficients",
        "calculate_point_roots",
        "get_minor_degrees", "get_roots_and_PMGBZ", "get_loop_winding",
        "get_strip_winding",
        "PolyDiffContext", "WindingFun", "MatWindingFun", "get_winding_number",
        "SGBZSolver", "SGBZChecker", "check_SGBZ",
        "convert_gbz_to_triplets", "convert_gbz_list_to_triplets",
        "PointSubset", "LineSubset", "GBZResult", "ConnectedSubset",
    ]
    for name in expected:
        assert hasattr(bfs, name), f"Missing export: {name}"
