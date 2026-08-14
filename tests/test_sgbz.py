"""Tests for the ZeroManager-based brute_force_SGBZ module.

The 2D HN model's SGBZ boundary at ``mu1 = gamma_1`` is a *continuum*
(a 1D LineSubset of equal-modulus degeneracy), not isolated points — so
the inside-spectrum tests assert ``is_continuum`` and the materialized
``LineSubset``s.  Spectrum membership (in vs out) is decided by the
winding-zero / left-right-limit bisection and is fully exercised here.
"""

import numpy as np
import pytest
from cmath import exp

import brute_force_SGBZ as bfs
from gbz_types import PointSubset, LineSubset, GBZResult, CharPoly
from continuation import ZeroManager


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


# ---- [10]-SGBZ spectrum-membership tests ----

class TestSGBZ10:
    """Tests for [10]-SGBZ of 2D HN model.

    The SGBZ boundary at mu1 = gamma_1 is a continuum (1D LineSubset);
    the materialized LineSubsets are asserted on (not just the flag).
    """

    def test_returns_gbzresult(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0)
        assert isinstance(gbz, GBZResult)
        assert gbz.E_ref == 1.0 + 0j

    def test_inside_spectrum_is_continuum(self, poly_A, params_A):
        """Continuum boundary materializes as 2 LineSubsets at |β₂|=exp(γ₂).

        Analytic (2D HN, J1=J2=1, δ=0): the M-1/M boundary pair is the
        ±β₂ pair on the circle |β₂| = exp(γ₂), valid over a sub-arc of
        θ₁.  Two LineSubsets (one per ±β₂), mu1 = γ₁, |β₂| = exp(γ₂).
        """
        coeffs, degs = poly_A
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0)
        assert gbz.is_gbz            # in spectrum
        assert gbz.is_continuum
        assert gbz.success

        lines = [s for s in gbz.subsets if isinstance(s, LineSubset)]
        assert len(lines) == 2
        assert gbz.index == (0, 2)
        # mu1 = gamma_1 (the bisection lands the continuum boundary).
        for s in lines:
            assert s.mu1 == pytest.approx(params_A["gamma_1"], abs=2e-3)
            # |β₂| = exp(γ₂) along the whole line (boundary-pair modulus).
            abs_b2 = np.abs(s.beta2_arr)
            assert np.allclose(abs_b2, np.exp(params_A["gamma_2"]), atol=1e-5)
        # The two lines are the ±β₂ pair: same |β₂|, β₂ endpoints conjugate
        # at the shared seam (θ₁ = 2π).
        assert lines[0].theta1_start == pytest.approx(lines[1].theta1_start)
        assert lines[0].theta1_end == pytest.approx(lines[1].theta1_end)

    def test_outside_spectrum(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 5.0 + 0j, 0.0)
        assert not gbz.is_gbz
        assert not gbz.is_continuum
        assert gbz.is_empty
        assert gbz.index == (0, 0)

    def test_solver_mu1_matches_analytic(self, poly_A, params_A):
        """The bisection locates mu1 = gamma_1 (the continuum boundary).

        The continuum is detected via the left/right-winding-limit straddle,
        which lands at mu1 = γ₁ − ~7e-9 (not exactly γ₁) because the
        ZeroManager root solver is slightly unstable near the degenerate
        continuum — see log/2026-08-13-sgbz-mu2mid-refactor.md "Blocked".
        The 2e-3 tolerance covers that offset and the continuum_perturb
        geometry (1e-2 × scales 1..8); it does NOT accept the non-physical
        eps~1e-8 regime flagged in that log, which remains an open issue.
        """
        coeffs, degs = poly_A
        poly = CharPoly(coeffs, degs)
        res = bfs.solve_SGBZ_for_E(poly, 1.0 + 0j)
        assert res["is_continuum"]
        assert res["mu1"] == pytest.approx(params_A["gamma_1"], abs=2e-3)

    def test_winding_signs_straddle_boundary(self, poly_A):
        """W is negative below gamma_1 and positive above it (plateau ±1)."""
        coeffs, degs = poly_A
        poly = CharPoly(coeffs, degs)
        E_ref = 1.0 + 0j

        zm_lo = ZeroManager(poly, E_ref, -0.5)
        zm_lo.run()
        assert not bfs.detect_continuum_simple(zm_lo, poly)
        _, w_lo = bfs.detect_crossings_and_winding(zm_lo, poly)
        assert w_lo < 0

        zm_hi = ZeroManager(poly, E_ref, 0.5)
        zm_hi.run()
        assert not bfs.detect_continuum_simple(zm_hi, poly)
        _, w_hi = bfs.detect_crossings_and_winding(zm_hi, poly)
        assert w_hi > 0


# ---- [11]-SGBZ tests ----

class TestSGBZ11:
    """Tests for [11]-SGBZ of 2D HN model."""

    def test_returns_gbzresult(self, poly_A_11):
        coeffs, degs = poly_A_11
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0)
        assert isinstance(gbz, GBZResult)

    def test_inside_spectrum_is_continuum(self, poly_A_11, params_A):
        """[11]-SGBZ: continuum at mu1 = γ₁+γ₂, materialized as 2 LineSubsets."""
        coeffs, degs = poly_A_11
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0)
        assert gbz.is_gbz
        assert gbz.is_continuum
        assert gbz.success
        lines = [s for s in gbz.subsets if isinstance(s, LineSubset)]
        assert len(lines) == 2
        assert gbz.index == (0, 2)
        for s in lines:
            assert s.mu1 == pytest.approx(
                params_A["gamma_1"] + params_A["gamma_2"], abs=2e-3)

    def test_outside_spectrum(self, poly_A_11):
        coeffs, degs = poly_A_11
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 5.0 + 0j, 0.0)
        assert not gbz.is_gbz

    def test_solver_mu1_matches_analytic(self, poly_A_11, params_A):
        """mu1 = gamma_1 + gamma_2 for the [11] direction."""
        coeffs, degs = poly_A_11
        poly = CharPoly(coeffs, degs)
        res = bfs.solve_SGBZ_for_E(poly, 1.0 + 0j)
        expected = params_A["gamma_1"] + params_A["gamma_2"]
        assert res["is_continuum"]
        assert res["mu1"] == pytest.approx(expected, abs=2e-3)


# ---- Hermitian limit (γ = 0 → GBZ = BZ, a continuum on the unit circle) ----

class TestHermitianLimit:
    """Hermitian system: GBZ = BZ, realized as a continuum (|β|=1 everywhere)."""

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

    def test_sgbz10_in_spectrum(self, poly_hermitian):
        coeffs, degs = poly_hermitian
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0)
        assert gbz.is_gbz
        assert gbz.is_continuum

    def test_sgbz11_in_spectrum(self, poly_hermitian_11):
        coeffs, degs = poly_hermitian_11
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0)
        assert gbz.is_gbz
        assert gbz.is_continuum

    def test_hermitian_all_equal(self, poly_hermitian, poly_hermitian_11):
        """Hermitian limit: [10]-SGBZ, [11]-SGBZ, Amoeba all in spectrum."""
        coeffs10, degs10 = poly_hermitian
        coeffs11, degs11 = poly_hermitian_11

        import brute_force_amoeba as bfa

        gbz10 = bfs.collect_GBZ_subsets(coeffs10, degs10, 1.0 + 0j, 0.0)
        gbz11 = bfs.collect_GBZ_subsets(coeffs11, degs11, 1.0 + 0j, 0.0)
        gbz_amoeba = bfa.collect_GBZ_subsets(coeffs10, degs10, 1.0 + 0j, 0.0)

        assert gbz10.is_gbz
        assert gbz11.is_gbz
        assert gbz_amoeba.is_gbz


# ---- 0D crossing detection (off the continuum boundary) ----

class TestCrossingDetection:
    """detect_crossings_simple on a non-continuum mu1: 0D PMGBZ points.

    At mu1 off the boundary the roots are non-degenerate; the M-1/M sorted
    positions cross transversally, yielding 2 PointSubsets per crossing
    with opposite charges.  This exercises the §2 cubic+Newton path that
    the inside-spectrum (continuum) tests do not reach.
    """

    def test_two_points_one_crossing(self, poly_A):
        """At mu1 = 0.1 (below gamma_1=0.2) one crossing → 2 PointSubsets.

        The two boundary tracks cross μ₂_mid at the same θ₁* but distinct
        θ₂ (two β₂) — two independent zeros, each with its own charge
        (§3.1: charge = sign(g'), one per zero-curve).  So one crossing of
        the boundary pair yields 2 PointSubsets AND 2 charge dicts with
        opposite ±1 charges.
        """
        coeffs, degs = poly_A
        poly = CharPoly(coeffs, degs)
        zm = ZeroManager(poly, 1.0 + 0j, 0.1)
        zm.run()
        assert not bfs.detect_continuum_simple(zm, poly)

        subsets, charges = bfs.detect_crossings_simple(zm, poly)
        assert len(subsets) == 2
        assert all(isinstance(s, PointSubset) for s in subsets)
        # Both points share the same beta1 (same theta1 crossing).
        assert subsets[0].beta1 == pytest.approx(subsets[1].beta1, rel=1e-9)
        # One charge dict per crossing track (§3.1): the boundary pair has
        # two tracks, so 2 ordinary charges with opposite sign.
        assert len(charges) == 2
        assert all(c["kind"] == "ordinary" for c in charges)
        assert {c["charge"] for c in charges} == {+1, -1}

    def test_charges_propagate_winding(self, poly_A):
        """compute_average_winding returns a finite, sign-consistent value."""
        coeffs, degs = poly_A
        poly = CharPoly(coeffs, degs)
        zm = ZeroManager(poly, 1.0 + 0j, 0.1)
        zm.run()
        subsets, W = bfs.detect_crossings_and_winding(zm, poly)
        assert np.isfinite(W)
        assert W < 0  # below gamma_1 → negative winding

    def test_detection_leaves_mesh_unchanged(self, poly_A):
        """The crossing phase must not mutate the mesh (2026-08-14 contract).

        A previous version inserted one row per bracket probe
        (``insert_solution``), shifting indices mid-sweep and letting the
        touch loop pair a β₂ with a stale θ₁.  Detection now solves probes
        transiently (``ZeroManager.solve_at``) — only ``build_mu2_mid``
        mutates, so the mesh row counts must be identical before and after
        ``detect_crossings_simple``.
        """
        coeffs, degs = poly_A
        poly = CharPoly(coeffs, degs)
        zm = bfs.Mu2MidZM(poly, 1.0 + 0j, 0.1)
        zm.run()
        zm.build_mu2_mid()
        assert not zm.has_continuum
        n_before = [len(seg.theta1_arr) for seg in zm.segments]

        subsets, charges = bfs.detect_crossings_simple(zm, poly)
        assert len(subsets) == 2  # sanity: the known crossing is still found

        n_after = [len(seg.theta1_arr) for seg in zm.segments]
        assert n_after == n_before


# ---- continuum LineSubset materialization (direct) ----

class TestContinuumMaterialization:
    """Direct exercise of extract_continuum_linesubsets on the 2D HN continuum.

    Localized signal for the §1/§2 build → boundary-run → join-across-MR
    pipeline, independent of the slow end-to-end collect_GBZ_subsets path.
    """

    def test_extract_two_lines_at_exp_gamma2(self, poly_A, params_A):
        """Build Mu2MidZM on the continuum and extract directly.

        The boundary-pair modulus on the continuum is exp(γ₂); two LineSubsets
        (the ±β₂ pair) result, with mu1 = γ₁.
        """
        coeffs, degs = poly_A
        poly = CharPoly(coeffs, degs)
        zm = bfs.Mu2MidZM(poly, 1.0 + 0j, params_A["gamma_1"])
        zm.run()
        zm.build_mu2_mid()
        assert zm.has_continuum

        lines = bfs.extract_continuum_linesubsets(zm, poly)
        assert len(lines) == 2
        for s in lines:
            assert s.mu1 == pytest.approx(params_A["gamma_1"], abs=1e-9)
            assert np.allclose(np.abs(s.beta2_arr),
                               np.exp(params_A["gamma_2"]), atol=1e-5)
        # both lines span the same θ₁ sub-arc
        assert lines[0].theta1_start == pytest.approx(lines[1].theta1_start)
        assert lines[0].theta1_end == pytest.approx(lines[1].theta1_end)


# ---- 0/∞ root truncation ----

def test_logabs_clamped_keeps_mu2_mid_finite():
    """A 0/∞ root in the boundary pair must not make μ₂_mid ±∞.

    ``solve_roots_1d`` pads 0/∞ roots into a degree-deficient 1-D
    polynomial; if such a padded root lands at sorted position M-1 or M the
    naive ``np.log|β|`` gives ±∞ and μ₂_mid = mean → ±∞, overflowing the
    winding loop.  ``logabs_clamped`` clamps to ±14 (aligned with
    arclength.INF_THRESHOLD=1e6) so every consumer sees a finite curve.

    Directly unit-tested because the HN models used elsewhere never produce a
    0/∞ boundary root; a synthetic root array is the honest fixture.
    """
    from brute_force_SGBZ.mu2mid import logabs_clamped, _LOGABS_CLAMP_L

    # one finite root, one 0, one ∞ — the 0/∞ would be ±∞ unclamped.
    roots = np.array([1.5 + 0.3j, 0.0 + 0.0j, np.inf + 0j])
    la = logabs_clamped(roots)
    assert np.all(np.isfinite(la))
    assert la[1] == -_LOGABS_CLAMP_L          # 0 root → −L
    assert la[2] == _LOGABS_CLAMP_L           # ∞ root → +L
    assert la[0] == pytest.approx(np.log(abs(1.5 + 0.3j)))  # finite root untouched

    # μ₂_mid = mean of the M-1/M boundary pair stays finite even when one is 0/∞
    mu2_mid = (la[1] + la[0]) / 2.0
    assert np.isfinite(mu2_mid)


# ---- export verification ----

def test_all_exports():
    """Verify all expected symbols are exported."""
    expected = [
        "CharPoly", "get_minor_degrees",
        "detect_continuum_simple",
        "extract_continuum_linesubsets",
        "detect_crossings_simple", "detect_crossings_and_winding",
        "compute_average_winding",
        "solve_SGBZ_for_E", "collect_GBZ_subsets",
        "PointSubset", "LineSubset", "GBZResult", "ConnectedSubset",
    ]
    for name in expected:
        assert hasattr(bfs, name), f"Missing export: {name}"
