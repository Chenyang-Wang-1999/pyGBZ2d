"""Tests for brute_force_SGBZ module using 2D HN model analytic solution."""

import numpy as np
import pytest
from cmath import exp, log

import brute_force_SGBZ as bfs
from gbz_types import PointSubset, LineSubset, GBZResult, CharPoly


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
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        assert isinstance(gbz, GBZResult)
        assert gbz.E_ref == 1.0 + 0j

    def test_inside_spectrum(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        assert gbz.is_gbz
        assert not gbz.is_empty
        assert gbz.index != (0, 0)

    def test_outside_spectrum(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 5.0 + 0j, 0.0, N_points=101)
        assert not gbz.is_gbz
        assert gbz.is_empty
        assert gbz.index == (0, 0)

    def test_mu1_matches_analytic(self, poly_A, params_A):
        coeffs, degs = poly_A
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        assert gbz.is_gbz
        for s in gbz.subsets:
            if isinstance(s, PointSubset):
                assert s.mu1 == pytest.approx(params_A["gamma_1"], rel=1e-5)
            elif isinstance(s, LineSubset):
                assert s.mu1 == pytest.approx(params_A["gamma_1"], rel=1e-5)

    def test_beta_magnitudes_analytic(self, poly_A, params_A):
        coeffs, degs = poly_A
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
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
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        for s in gbz.subsets:
            assert isinstance(s, (PointSubset, LineSubset))

    def test_index_consistent(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        n_0d = sum(1 for s in gbz.subsets if isinstance(s, PointSubset))
        n_1d = sum(1 for s in gbz.subsets if isinstance(s, LineSubset))
        assert gbz.index == (n_0d, n_1d)


# ---- [11]-SGBZ tests ----

class TestSGBZ11:
    """Tests for [11]-SGBZ of 2D HN model."""

    def test_returns_gbzresult(self, poly_A_11):
        coeffs, degs = poly_A_11
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        assert isinstance(gbz, GBZResult)

    def test_inside_spectrum(self, poly_A_11):
        coeffs, degs = poly_A_11
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0, N_points=101)
        assert gbz.is_gbz
        assert not gbz.is_empty

    def test_outside_spectrum(self, poly_A_11):
        coeffs, degs = poly_A_11
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 5.0 + 0j, 0.0, N_points=101)
        assert not gbz.is_gbz

    def test_mu1_matches_analytic(self, poly_A_11, params_A):
        coeffs, degs = poly_A_11
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        assert gbz.is_gbz
        expected_mu = params_A["gamma_1"] + params_A["gamma_2"]  # = 0.5
        for s in gbz.subsets:
            if isinstance(s, PointSubset):
                assert s.mu1 == pytest.approx(expected_mu, rel=1e-5)
            elif isinstance(s, LineSubset):
                assert s.mu1 == pytest.approx(expected_mu, rel=1e-5)

    def test_beta1_magnitude_analytic(self, poly_A_11, params_A):
        coeffs, degs = poly_A_11
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        expected_r = exp(params_A["gamma_1"] + params_A["gamma_2"])
        for s in gbz.subsets:
            if isinstance(s, PointSubset):
                assert abs(s.beta1) == pytest.approx(expected_r, rel=1e-4)


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
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        assert gbz.is_gbz
        for s in gbz.subsets:
            mu = s.mu1 if isinstance(s, PointSubset) else s.mu1
            assert mu == pytest.approx(0.0, abs=1e-5)

    def test_sgbz10_beta_on_unit_circle(self, poly_hermitian):
        coeffs, degs = poly_hermitian
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        for s in gbz.subsets:
            if isinstance(s, PointSubset):
                assert abs(s.beta1) == pytest.approx(1.0, rel=1e-4)
                assert abs(s.beta2) == pytest.approx(1.0, rel=1e-4)

    def test_sgbz11_mu1_zero(self, poly_hermitian_11):
        coeffs, degs = poly_hermitian_11
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        assert gbz.is_gbz
        for s in gbz.subsets:
            mu = s.mu1 if isinstance(s, PointSubset) else s.mu1
            assert mu == pytest.approx(0.0, abs=1e-5)

    def test_sgbz11_beta_on_unit_circle(self, poly_hermitian_11):
        coeffs, degs = poly_hermitian_11
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0, N_points=201)
        for s in gbz.subsets:
            if isinstance(s, PointSubset):
                assert abs(s.beta1) == pytest.approx(1.0, rel=1e-4)

    def test_hermitian_all_equal(self, poly_hermitian, poly_hermitian_11):
        """Hermitian limit: [10]-SGBZ, [11]-SGBZ, Amoeba all give same GBZ (=BZ)."""
        coeffs10, degs10 = poly_hermitian
        coeffs11, degs11 = poly_hermitian_11

        import brute_force_amoeba as bfa

        gbz10 = bfs.collect_GBZ_subsets(coeffs10, degs10, 1.0 + 0j, 0.0, N_points=201)
        gbz11 = bfs.collect_GBZ_subsets(coeffs11, degs11, 1.0 + 0j, 0.0, N_points=201)
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
        "poly_to_np_coefficients",
        "calculate_point_roots",
        "get_minor_degrees", "get_roots_and_PMGBZ", "get_loop_winding",
        "get_strip_winding",
        "CharPoly", "WindingFun", "MatWindingFun", "get_winding_number",
        "solve_SGBZ_for_E", "collect_GBZ_subsets",
        "PointSubset", "LineSubset", "GBZResult", "ConnectedSubset",
    ]
    for name in expected:
        assert hasattr(bfs, name), f"Missing export: {name}"


# ---- Loop winding theta2 consistency ----

class TestLoopWindingTheta2:
    """loop winding w0 must be consistent (same integer) across theta2 values
    not equal to pi, where the integrand may have numerical spikes.

    Regression test for Haldane-model x-direction (4-band supercell) at
    mu1 values where PMGBZ detection returns no points.  Without the fix
    theta2_median could land on π, producing w0=0 or w0=2 instead of the
    correct w0=+1, which then causes a fake winding zero and false-positive
    plateau detection.
    """

    # Haldane gain-loss model, x-direction supercell coefficients at
    # E = -1.868 - 0.3162j.  The 4-band supercell creates a band anti-crossing
    # at mu1 ≈ 0.27–0.28 where PMGBZ points disappear.
    _coeffs = np.array([
        1.0, -1.0, -1.0,
        0.75 - 2.77555756156289e-17j, 0.25,
        -2.36602540378444 + 5.55111512312578e-17j, 0.25,
        -1.0, -7.5, -1.0,
        0.25, -0.633974596215562 + 1.11022302462516e-16j, 0.25,
        0.75 + 2.77555756156289e-17j, -0.25,
        -0.5 + 5.57781305375155e-18j, 0.43301270189222 + 1.11022302462516e-16j,
        -0.5 + 5.57781305375155e-18j, -0.0669872981077807 + 2.77555756156289e-17j,
        2.61602540378444 + 1.11022302462516e-16j, -0.0669872981077807,
        -1.5, -1.0, -1.5 + 5.55111512312578e-17j,
        -0.933012701892219 + 1.11022302462516e-16j, 0.883974596215562,
        -0.933012701892220, -0.5 - 5.57781305375155e-18j,
        -0.433012701892219 - 8.32667268468867e-17j, -0.5 - 5.57781305375155e-18j,
        -0.25, 0.0625, 0.0625000000000001 + 6.93889390390723e-18j,
        -0.59150635094611 + 1.11022302462516e-16j, 0.0625000000000001 + 6.93889390390723e-18j,
        0.0625, -0.625 + 4.16333634234434e-17j, -0.821474596215561 + 2.77555756156289e-17j,
        -0.625 + 4.85722573273506e-17j, 0.0625,
        -0.125, -1.62948729810778 + 1.94289029309402e-16j,
        5.30560796608386 + 1.11022302462516e-16j, -1.62948729810778 + 8.32667268468867e-17j,
        -0.125, 0.1875,
        0.999999999999999 - 2.77555756156289e-16j, 12.375 - 1.11022302462516e-16j,
        0.999999999999999 - 2.77555756156289e-16j, 0.1875,
        -0.125, -2.49551270189222 - 2.4980018054066e-16j,
        -2.05560796608386 - 5.55111512312578e-16j, -2.49551270189222 - 2.77555756156289e-17j,
        -0.125, 0.0625,
        -0.624999999999999 - 1.11022302462516e-16j, -2.55352540378444 - 5.55111512312578e-16j,
        -0.625 + 5.55111512312578e-17j, 0.0625,
        0.0625000000000001 - 6.93889390390723e-18j, -0.15849364905389 - 1.11022302462516e-16j,
        0.0625000000000001 - 6.93889390390723e-18j, 0.0625,
    ], dtype=complex)

    _degs = np.array([
        [4, 0, 0], [3, 1, 0], [3, -1, 0], [2, 2, 0], [2, 1, 1],
        [2, 1, 0], [2, 1, -1], [2, 0, 1], [2, 0, 0], [2, 0, -1],
        [2, -1, 1], [2, -1, 0], [2, -1, -1], [2, -2, 0], [1, 3, 0],
        [1, 2, 1], [1, 2, 0], [1, 2, -1], [1, 1, 1], [1, 1, 0],
        [1, 1, -1], [1, 0, 1], [1, 0, 0], [1, 0, -1], [1, -1, 1],
        [1, -1, 0], [1, -1, -1], [1, -2, 1], [1, -2, 0], [1, -2, -1],
        [1, -3, 0], [0, 4, 0], [0, 3, 1], [0, 3, 0], [0, 3, -1],
        [0, 2, 2], [0, 2, 1], [0, 2, 0], [0, 2, -1], [0, 2, -2],
        [0, 1, 2], [0, 1, 1], [0, 1, 0], [0, 1, -1], [0, 1, -2],
        [0, 0, 2], [0, 0, 1], [0, 0, 0], [0, 0, -1], [0, 0, -2],
        [0, -1, 2], [0, -1, 1], [0, -1, 0], [0, -1, -1], [0, -1, -2],
        [0, -2, 2], [0, -2, 1], [0, -2, 0], [0, -2, -1], [0, -2, -2],
        [0, -3, 1], [0, -3, 0], [0, -3, -1], [0, -4, 0],
    ], dtype=int)

    _E_ref = -1.868 - 0.3162j

    @staticmethod
    def _compute_w0_vs_theta2(poly_diff, E_ref, mu1, mu2_fun, n_theta=17, n_fine=10001):
        """Return list of (theta2, w0) for n_theta samples around [0, 2π)."""
        from math import pi
        from cmath import exp
        results = []
        for theta2 in np.linspace(0, 2 * pi, n_theta):
            # Direct integration on fine grid
            t_fine = np.linspace(0, 2 * pi, n_fine)
            dt = t_fine[1] - t_fine[0]
            w0 = 0.0
            for t in t_fine:
                m2v = float(mu2_fun(t, extrapolate="periodic"))
                m2p = float(mu2_fun(t, extrapolate="periodic", nu=1))
                beta1 = exp(mu1 + 1j * t)
                beta2 = exp(m2v + 1j * theta2)
                val = poly_diff.eval_val((E_ref, beta1, beta2))
                partials = poly_diff.eval_partials((E_ref, beta1, beta2))
                dbeta1_dt = 1j * beta1
                dbeta2_dt = beta2 * m2p
                dval_dt = partials[1] * dbeta1_dt + partials[2] * dbeta2_dt
                w0 += (dval_dt / val).imag * dt
            w0 /= (2 * pi)
            results.append((theta2, w0))
        return results

    def _build_poly_and_mu2(self, mu1):
        """Build CharPoly + mu2_fun at given mu1."""
        from scipy import interpolate
        poly = CharPoly(self._coeffs, self._degs)

        gbz, theta1_arr, sols_arr, info = bfs.get_roots_and_PMGBZ(
            poly, self._E_ref, mu1, N_points=301,
        )
        M = poly.M
        # Sort by norm
        for row_ind in range(sols_arr.shape[0]):
            sols_arr[row_ind, :] = sols_arr[row_ind,
                np.argsort(np.abs(sols_arr[row_ind, :]))]

        mu2_max = np.log(np.abs(sols_arr[:, M]))
        mu2_min = np.log(np.abs(sols_arr[:, M - 1]))
        mu2_max[-1] = mu2_max[0]
        mu2_min[-1] = mu2_min[0]
        mu2_arr = (mu2_max + mu2_min) / 2
        mu2_fun = interpolate.CubicSpline(theta1_arr, mu2_arr, bc_type="periodic")

        return poly, mu2_fun, info

    def test_w0_consistent_mu1_0270(self):
        """At mu1=0.270 (no PMGBZ), w0 should be +1 at all theta2 ≠ π."""
        poly_diff, mu2_fun, info = self._build_poly_and_mu2(0.270)
        assert info["continuum_flag"] is False
        # Continuation adaptive mesh may detect PMGBZ points that the old
        # 301-point uniform mesh missed.  Accept 0 or 1 as valid.
        assert len(info.get("_pmgbz_raw", [])) <= 1

        results = self._compute_w0_vs_theta2(poly_diff, self._E_ref, 0.270, mu2_fun)
        # All w0 values should round to 1, except possibly theta2=π
        w0_vals = [w for _, w in results]
        # Exclude theta2=π (index 8 of 17, value ≈ 3.1416)
        w0_no_pi = [w for i, (t, w) in enumerate(results) if abs(t - np.pi) > 1e-3]
        consensus = int(round(np.median(w0_no_pi)))
        assert consensus == 1, f"Expected consensus w0=1, got {consensus}"

        # Verify majority: at least 15 of 17 agree on consensus
        n_agree = sum(1 for _, w in results if int(round(w)) == consensus)
        assert n_agree >= 15, f"Only {n_agree}/17 agree on w0={consensus}"

        # theta2=π gives wrong value (0 instead of 1) — this is the bug
        w0_pi = [w for t, w in results if abs(t - np.pi) < 1e-3]
        assert len(w0_pi) == 1
        assert abs(w0_pi[0] - 0.0) < 0.1, (
            f"Expected w0≈0 at theta2=π (numerical artifact), got {w0_pi[0]:.4f}"
        )

    def test_w0_consistent_mu1_0280(self):
        """At mu1=0.280 (no PMGBZ), w0 should be +1 at all theta2 ≠ π."""
        poly_diff, mu2_fun, info = self._build_poly_and_mu2(0.280)
        assert info["continuum_flag"] is False
        # Continuation adaptive mesh may detect PMGBZ points that the old
        # 301-point uniform mesh missed.  Accept 0 or 1 as valid.
        assert len(info.get("_pmgbz_raw", [])) <= 1

        results = self._compute_w0_vs_theta2(poly_diff, self._E_ref, 0.280, mu2_fun)
        w0_no_pi = [w for i, (t, w) in enumerate(results) if abs(t - np.pi) > 1e-3]
        consensus = int(round(np.median(w0_no_pi)))
        assert consensus == 1, f"Expected consensus w0=1, got {consensus}"

        n_agree = sum(1 for _, w in results if int(round(w)) == consensus)
        assert n_agree >= 15, f"Only {n_agree}/17 agree on w0={consensus}"

        # theta2=π gives wrong value (2 instead of 1) — this is the bug
        w0_pi = [w for t, w in results if abs(t - np.pi) < 1e-3]
        assert len(w0_pi) == 1
        assert abs(w0_pi[0] - 2.0) < 0.1, (
            f"Expected w0≈2 at theta2=π (numerical artifact), got {w0_pi[0]:.4f}"
        )

    def test_strip_winding_monotonic(self):
        """Strip winding should be continuous near mu1=0.27–0.28 without spurious zeros."""
        w_vals = []
        for mu1 in [0.20, 0.24, 0.26, 0.27, 0.28, 0.29, 0.30, 0.35]:
            w, _ = bfs.get_strip_winding(
                self._build_poly_and_mu2(mu1)[0], self._E_ref, mu1, N_points=301,
            )
            if w is None:
                continue  # continuum detected — skip (unexpected in this range)
            if isinstance(w, tuple):
                w = float(np.nanmean(w))
            w_vals.append(float(w))
        # All windings should be > 0 (no fake zero at mu1=0.27)
        for w in w_vals:
            assert w > 0.01, f"Strip winding {w:.6f} ≤ 0.01 — spurious zero"
