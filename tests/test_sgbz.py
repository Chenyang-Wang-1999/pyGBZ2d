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
from brute_force_SGBZ import sgbz_solver, winding as sgbz_winding
from brute_force_SGBZ import pairwise as sgbz_pairwise
from gbz_types import PointSubset, LineSubset, GBZResult, CharPoly, TWO_PI
from continuation import ZeroManager

from conftest import build_HN2D_polynomial


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

    def test_extract_requires_continuum_detected(self, poly_A):
        """Materialization without a detected continuum must raise, not return []."""
        from brute_force_SGBZ import continuum_lines

        coeffs, degs = poly_A
        poly = CharPoly(coeffs, degs)
        zm = bfs.Mu2MidZM(poly, 1.0 + 0j, 0.1)
        zm.run()
        zm.build_mu2_mid()
        assert not zm.has_continuum

        with pytest.raises(RuntimeError, match="has_continuum=True"):
            continuum_lines.extract_continuum_linesubsets(zm, poly)

    def test_extract_raises_when_no_boundary_runs(
        self, poly_A, params_A, monkeypatch,
    ):
        """has_continuum=True with zero boundary runs is an invariant violation."""
        from brute_force_SGBZ import continuum_lines

        coeffs, degs = poly_A
        poly = CharPoly(coeffs, degs)
        zm = bfs.Mu2MidZM(poly, 1.0 + 0j, params_A["gamma_1"])
        zm.run()
        zm.build_mu2_mid()
        assert zm.has_continuum

        monkeypatch.setattr(continuum_lines, "_find_boundary_runs", lambda m: [])
        with pytest.raises(RuntimeError, match="invariant broken"):
            continuum_lines.extract_continuum_linesubsets(zm, poly)


# ---- pairwise EventGroup analysis (2026-08-16 refactor) ----

class TestPairwiseAnalysis:
    """Pairwise ItemView intersections, EventGroups and side-change charges."""

    def test_mu01_has_single_seam_event_group(self, poly_A):
        poly = CharPoly(*poly_A)
        zm = bfs.Mu2MidZM(poly, 1.0 + 0j, 0.1)
        zm.run()
        zm.analyze()
        assert len(zm._pair_events) == 1
        ev = zm._pair_events[0]
        assert ev.kind == 'cross'
        assert ev.converged
        assert ev.pair_kind == 'M-1_M'
        assert len(zm._event_groups) == 1
        g = zm._event_groups[0]
        assert set(g.point_columns) == {0, 1}
        assert g.column_q[0] == 1
        assert g.column_q[1] == -1

    def test_minimum_direction_deriv_protection(self):
        assert sgbz_pairwise._protected_direction(0.0, 1e-12) is None
        assert sgbz_pairwise._protected_direction(float('nan'), 1e-12) is None
        assert sgbz_pairwise._protected_direction(1e-11, 1e-12) == 1
        assert sgbz_pairwise._protected_direction(-1e-11, 1e-12) == -1

    def test_close_events_merge_without_dedup(self):
        class FakeSeg:
            theta1_arr = np.array([0.0, 0.1, 0.2])

        class FakeZM:
            segments = [FakeSeg]

        ev1 = sgbz_pairwise.PairEvent(
            seg_idx=0, i=0, ia=0, ib=1, cols_a=(0,), cols_b=(1,),
            rep_a=0, rep_b=1, kind='cross', pair_kind='M-1_M',
            theta_star=0.1, direction=1)
        ev2 = sgbz_pairwise.PairEvent(
            seg_idx=0, i=1, ia=1, ib=2, cols_a=(1,), cols_b=(2,),
            rep_a=1, rep_b=2, kind='cross', pair_kind='M-1_M',
            theta_star=0.1 + 0.5e-10, direction=-1)
        groups = sgbz_pairwise.group_events(FakeZM, [ev1, ev2], merge_tol=1e-10)
        assert len(groups) == 1
        assert len(groups[0].events) == 2
        # the merged group snapped to the straddled mesh data point 0.1
        assert groups[0].theta == pytest.approx(0.1)
        assert groups[0].column_components == ((0, 1, 2),)

    def test_mr_boundary_touch_is_collected_not_skipped(self):
        """A modulus coincidence on an MR boundary row is an event (2026-08-18).

        ``collect_pair_events`` historically skipped touch events landing on
        a segment's MR boundary rows (``left_mr/right_mr >= 0``), deferring
        to the MR materialization channel — but that channel only
        materialized clusters straddling M-1/M, so any OTHER modulus
        coincidence at the boundary row was dropped by both channels (the
        seam missed-detection behind the E=1.212 failure).  Now the touch
        flows through the same EventGroup machinery, marked ``is_mr``, with
        ``direction=None`` from the divergent cluster tangents (→ hard,
        charge None downstream).

        poly_F ``f = β₂² − 2β₂ − β₁ + 2`` at μ₁=0: the roots
        ``β₂ = 1 ± √(β₁ − 1)`` coincide exactly at θ₁=0 (β₁=1), so run()
        records a boundary MR and segment 0's row 0 carries the snapped
        equal-modulus pair → one touch at i=0 with left_mr=0.
        """
        coeffs = np.array([1, -2, -1, 2], dtype=complex)
        degs = np.array([
            [0, 0, 2], [0, 0, 1], [0, 1, 0], [0, 0, 0],
        ], dtype=int)
        poly = CharPoly(coeffs, degs)
        zm = bfs.Mu2MidZM(poly, 0j, 0.0)
        zm.run(h0=0.1, cluster_tol=1e-4, min_dtheta=1e-6)
        zm.analyze()

        assert zm.segments[0].left_mr == 0  # row 0 IS the boundary MR row
        mr_events = [e for e in zm._pair_events if e.is_mr]
        assert len(mr_events) == 1
        ev = mr_events[0]
        assert ev.kind == 'touch'
        assert ev.seg_idx == 0 and ev.i == 0
        assert ev.theta_star == pytest.approx(0.0, abs=1e-9)
        # Cluster tangents at the MR row are inf → direction is None
        # (hard boundary, charge None downstream).
        assert ev.direction is None

        mr_groups = [g for g in zm._event_groups if g.is_mr]
        assert len(mr_groups) == 1
        assert mr_groups[0].theta == pytest.approx(0.0, abs=1e-9)

    def test_seam_merge_relabels_2pi_side_events_to_theta0_frame(self):
        """2π-side events merged into the θ=0 group must swap column frames.

        ``group_events``' seam merge anchors the merged group at θ=0, but
        events detected near 2π carry segment track-frame labels, which
        differ from the θ=0 frame by ``boundary_perm`` (the closing row
        stores ``left_boundary_roots`` permuted into the track frame).
        Without relabelling, finalize's frame-0 interpretation of the
        2π-side labels hit the wrong tracks and the boundary coverage test
        silently dropped the event (the missed seam crossing behind the
        E=1.212 y-SGBZ failure).

        Synthetic two-column setup: moduli coincide exactly at both seam
        rows (θ=0 and θ=2π — physically the same point), splitting in
        between, with a NON-trivial boundary_perm [1, 0] so the two frames
        genuinely differ.  The 2π-side event (columns 0,1 in track frame)
        must be relabelled to (1,0) in the θ=0-anchored merged group.
        """
        from brute_force_SGBZ.mu2mid import ItemView
        from dataclasses import replace as _replace

        th = np.array([0.0, 0.1, 0.2, 6.183185307179586,
                       6.283185307179586])  # 2π - 0.1, then 2π
        twopi = 6.283185307179586

        # Modulus curves: equal at row 0 (θ=0) and row 4 (θ=2π), split in
        # between — crossing in interval [3,4] (2π side) only.
        la = np.array([
            [0.0, 0.0],       # θ=0: equal (seam degeneracy)
            [-0.05, 0.05],
            [-0.08, 0.08],
            [-0.04, 0.04],
            [0.0, 0.0],       # θ=2π: equal again
        ])

        class FakeSeg:
            theta1_arr = th
            left_mr = -1
            right_mr = -1

        class FakeZM:
            segments = [FakeSeg()]
            M = 1
            K = 2
            boundary_perm = np.array([1, 0])   # NON-trivial monodromy
            _continuum_clusters = [[]]

            def __init__(self):
                self._item_views = [ItemView(
                    rep_cols=np.array([0, 1]),
                    mults=np.array([1, 1]),
                    item_logabs=la,
                    item_tang_re=np.zeros_like(la),
                    sort_to_item=np.argsort(la, axis=1),
                    j_lo=np.array([0, 0, 0, 0, 0]),
                    j_hi=np.array([1, 1, 1, 1, 1]),
                )]

        zm = FakeZM()
        # Event detected on the 2π side (interval [3,4]) with track-frame
        # labels (0, 1); θ* refines to exactly 2π.
        ev_2pi = sgbz_pairwise.PairEvent(
            seg_idx=0, i=3, ia=0, ib=1, cols_a=(0,), cols_b=(1,),
            rep_a=0, rep_b=1, kind='cross', pair_kind='M-1_M',
            theta_star=twopi, direction=1)

        groups = sgbz_pairwise.group_events(zm, [ev_2pi], merge_tol=1e-10)
        # Single event, no seam counterpart on the 0 side: no merge — and
        # crucially NO relabelling either (nothing moves to θ=0).
        assert len(groups) == 1
        assert groups[0].events[0].cols_a == (0,)

        # Now the full scenario: a matching event on the 0 side (touch at
        # row 0, frame-0 labels) so the seam merge fires.
        ev_0 = sgbz_pairwise.PairEvent(
            seg_idx=0, i=0, ia=0, ib=1, cols_a=(0,), cols_b=(1,),
            rep_a=0, rep_b=1, kind='touch', pair_kind='M-1_M',
            theta_star=0.0, direction=None)

        groups = sgbz_pairwise.group_events(zm, [ev_0, ev_2pi], merge_tol=1e-10)
        merged = groups[0]
        assert merged.theta == 0.0
        # The 2π-side event comes FIRST in the tuple (last.events +
        # first.events) and must now carry θ=0-frame labels: inv_perm maps
        # 0→1, 1→0 for boundary_perm=[1,0].
        e2, e0 = merged.events[0], merged.events[1]
        assert e2.theta_star == pytest.approx(twopi)   # provenance kept
        assert e2.cols_a == (1,) and e2.cols_b == (0,)
        assert e2.rep_a == 1 and e2.rep_b == 0
        # The 0-side event is already frame-0: untouched.
        assert e0.cols_a == (0,) and e0.cols_b == (1,)
        assert e0.theta_star == pytest.approx(0.0)

    def test_transient_inf_row_does_not_discard_pair_events(self):
        """A transient 0/∞ root must not kill the pair's other crossings.

        ``collect_pair_events`` historically skipped a pair whenever ANY
        row of ``ln|β₂_a| − ln|β₂_b|`` was non-finite, so one transient
        padding root (a degree-deficient solve at a single θ, later
        re-matched to a real root) silently discarded the pair's genuine
        crossings elsewhere on the segment.  The scan must instead drop
        only the mesh intervals that touch the bad row.  A persistent
        padding column (every row non-finite) still yields no events.

        Synthetic fixture (HN models never produce transient 0/∞ roots):
        item a is ∞ at row 0 only; the pair has one transversal crossing
        at θ* = 1/3 inside interval [0.3, 0.4].
        """
        from brute_force_SGBZ.mu2mid import ItemView

        th = np.array([0.0, 0.1, 0.2, 0.3, 0.4])

        def roots_at(t):
            # item a: e^{1-3t} (inf at t=0 via the padded row below);
            # item b: e^{-1+3t} i — log-moduli cross exactly at t = 1/3.
            return np.array([np.exp(1.0 - 3.0 * t),
                             np.exp(-1.0 + 3.0 * t) * 1j])

        def build_zm(transient: bool):
            class FakeSeg:
                theta1_arr = th
                tracked_roots = np.array([
                    (np.array([np.inf + 0j, np.exp(-1.0) * 1j])
                     if transient and i == 0 else roots_at(t))
                    for i, t in enumerate(th)
                ])
                tangents = np.array([
                    [np.nan, 3.0] if (transient and i == 0) else [-3.0, 3.0]
                    for i in range(len(th))
                ], dtype=complex)
                left_mr = -1
                right_mr = -1

            la = np.log(np.abs(FakeSeg.tracked_roots))
            tang_re = FakeSeg.tangents.real
            sort_to_item = np.argsort(la, axis=1)  # mults all 1

            class FakeZM:
                segments = [FakeSeg()]
                M = 1
                _continuum_clusters = [[]]

                def __init__(self):
                    self._item_views = [ItemView(
                        rep_cols=np.array([0, 1]),
                        mults=np.array([1, 1]),
                        item_logabs=la,
                        item_tang_re=tang_re,
                        sort_to_item=sort_to_item,
                        j_lo=sort_to_item[:, 0],
                        j_hi=sort_to_item[:, 1],
                    )]

                def solve_at(self, t, seg_idx=0, i=0, interp=None):
                    return roots_at(float(t)), np.array([-3.0, 3.0],
                                                         dtype=complex)

            return FakeZM()

        events = sgbz_pairwise.collect_pair_events(build_zm(True))
        assert len(events) == 1
        ev = events[0]
        assert ev.kind == 'cross'
        assert ev.converged
        assert ev.theta_star == pytest.approx(1.0 / 3.0, abs=1e-9)
        assert ev.direction == -1
        assert ev.pair_kind == 'M-1_M'

        # Persistent padding (item a non-finite on EVERY row): still no
        # events — the historical behaviour is the degenerate case.
        zm = build_zm(False)
        zm._item_views[0].item_logabs[:, 0] = np.inf
        assert sgbz_pairwise.collect_pair_events(zm) == []

    def test_isolated_tie_row_does_not_set_jlo_equal_jhi(self, poly_A):
        """ItemView clusters are whole-segment continua, not isolated ties."""
        poly = CharPoly(*poly_A)
        zm = bfs.Mu2MidZM(poly, 1.0 + 0j, 0.1)
        zm.run()
        zm.analyze()
        assert not zm.has_continuum
        # the seam event row is a tie, yet it is NOT injected into ItemView:
        # every row still has two distinct items at M-1 and M.
        for view in zm._item_views:
            assert not np.any(view.j_lo == view.j_hi)

    def test_mu2_mid_is_bounded_and_piecewise(self, poly_A):
        poly = CharPoly(*poly_A)
        zm = bfs.Mu2MidZM(poly, 1.0 + 0j, 0.1)
        zm.run()
        zm.analyze()
        path = zm.mu2_mid
        assert isinstance(path, bfs.Mu2Mid)
        assert len(path.pieces) > 0
        vals = np.array([path.value(t) for t in np.linspace(0, TWO_PI, 401)])
        assert vals.min() >= -14.0
        assert vals.max() <= 14.0


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


# ---- 2026-08-15 review fixes (B.1–B.5 + charge=None) ----

class TestReviewFixes:
    """Regression tests for the 2026-08-15 solver/winding review fixes."""

    # ---- B.1: the old cubic-Hermite bracket iteration is gone; analyze() replaces it ----

    def test_analyze_replaces_old_cubic_hermite_iteration(self, poly_A):
        poly = CharPoly(*poly_A)
        zm = bfs.Mu2MidZM(poly, 1.0 + 0j, 0.1)
        zm.run()
        zm.analyze()
        assert not hasattr(zm, '_cubic_hermite_iterate')
        # the known mu1=0.1 boundary crossing is still found through
        # pairwise EventGroups, with the seam crossing materialized at θ=0.
        assert len(zm._event_groups) == 1
        assert set(zm._event_groups[0].point_columns) == {0, 1}

    # ---- B.3: bracket expansion caps (left side) ----

    def test_left_bracket_expansion_raises_after_cap(self, poly_A, monkeypatch):
        poly = CharPoly(*poly_A)
        calls = {"n": 0}

        def fake_eval(poly_, E_ref, mu1, zm_run_kwargs, *,
                      continuum_tol, crossing_tol):
            calls["n"] += 1
            return 1.0, [], None   # W never crosses zero

        monkeypatch.setattr(sgbz_solver, "_evaluate_winding", fake_eval)
        monkeypatch.setattr(sgbz_solver, "_MAX_BRACKET_EXPANSIONS", 3)

        with pytest.raises(RuntimeError, match="left bracket expansion"):
            sgbz_solver.solve_SGBZ_for_E(poly, 1.0 + 0j, mu1_guess=(0.0, 1.0))
        assert calls["n"] == 3

    # ---- B.2: right bracket keeps expanding after a same-sign continuum proxy ----

    def test_right_bracket_keeps_expanding_after_same_sign_continuum_proxy(
            self, poly_A, monkeypatch):
        poly = CharPoly(*poly_A)
        seen = []

        def fake_eval(poly_, E_ref, mu1, zm_run_kwargs, *,
                      continuum_tol, crossing_tol):
            seen.append(float(mu1))
            if mu1 < 0.0:
                return -1.0, [], None
            return None, None, None   # continuum on the right side

        def fake_resolve(*args, **kwargs):
            # same-sign non-boundary: proxy = (mu1 + 0.1, -1.0)
            return -1.0, -1.0, 0.1

        monkeypatch.setattr(sgbz_solver, "_evaluate_winding", fake_eval)
        monkeypatch.setattr(
            sgbz_solver, "_resolve_continuum_winding", fake_resolve)
        monkeypatch.setattr(sgbz_solver, "_MAX_BRACKET_EXPANSIONS", 3)

        with pytest.raises(RuntimeError, match="right bracket expansion"):
            sgbz_solver.solve_SGBZ_for_E(poly, 1.0 + 0j, mu1_guess=(-1.0, 0.0))
        # The corrected band-edge proxy is re-evaluated (not skipped past by
        # the generic +1 step): -1 establishes the left bracket, then the
        # right loop visits 0.0, 0.1, 0.2 before the cap fires.
        assert seen == [-1.0, 0.0, 0.1, 0.2]

    # ---- B.2: a zero right endpoint found via a continuum proxy exits cleanly ----

    def test_right_endpoint_zero_from_continuum_proxy(self, poly_A, monkeypatch):
        poly = CharPoly(*poly_A)

        def fake_eval(poly_, E_ref, mu1, zm_run_kwargs, *,
                      continuum_tol, crossing_tol):
            if mu1 == -1.0:
                return -1.0, [], None
            if mu1 == 0.0:
                return None, None, None   # continuum at the right guess
            if mu1 == 0.1:
                return 0.0, [], None      # proxy correction lands on W = 0
            raise AssertionError(f"unexpected mu1={mu1!r}")

        def fake_resolve(*args, **kwargs):
            # non-boundary: w_l = -1, w_r = 0 at mu1 + eps = 0.1
            return -1.0, 0.0, 0.1

        monkeypatch.setattr(sgbz_solver, "_evaluate_winding", fake_eval)
        monkeypatch.setattr(
            sgbz_solver, "_resolve_continuum_winding", fake_resolve)

        res = sgbz_solver.solve_SGBZ_for_E(
            poly, 1.0 + 0j, mu1_guess=(-1.0, 0.0))
        assert res["_exit_reason"] == "right_endpoint_zero"
        assert res["mu1"] == pytest.approx(0.1)
        assert res["subsets"] == []

    # ---- B.4: ensure_mu2mid analyzes with CONTINUUM_TOL ----

    def test_ensure_mu2mid_uses_continuum_tol(self, poly_A, monkeypatch):
        poly = CharPoly(*poly_A)
        zm = ZeroManager(poly, 1.0 + 0j, 0.1)
        zm.run()

        seen = {}
        orig_analyze = bfs.Mu2MidZM.analyze

        def spy_analyze(self_, *args, **kwargs):
            seen["tie_tol"] = kwargs.get("tie_tol")
            return orig_analyze(self_, *args, **kwargs)

        monkeypatch.setattr(bfs.Mu2MidZM, "analyze", spy_analyze)
        subsets, charges = bfs.detect_crossings_simple(zm, poly)

        assert seen["tie_tol"] == bfs.CONTINUUM_TOL
        assert len(subsets) == 2   # the known mu1=0.1 crossing is still found

    # ---- hard-boundary charge is None (unknown), not a numeric placeholder ----
    #
    # The legacy _classify_charge direct test was removed with crossings.py.
    # Both semantics it guarded are covered on the real path:
    #   * ordinary charges are ±1 — asserted in TestSgbzSolver10 above;
    #   * the None sentinel never enters arithmetic — that is exactly what
    #     test_compute_average_winding_skips_none_hard_charge (B.5) checks.

    # ---- B.5: soft-only charge conservation ----

    def test_compute_average_winding_rejects_nonconserved_charges(self, poly_A):
        poly = CharPoly(*poly_A)
        zm = bfs.Mu2MidZM(poly, 1.0 + 0j, 0.1)
        zm.run()
        zm.build_mu2_mid()

        bad = [
            dict(theta1=0.0, theta2=0.5, charge=1, kind="ordinary"),
            dict(theta1=1.0, theta2=1.5, charge=1, kind="ordinary"),
        ]
        with pytest.raises(RuntimeError, match="do not sum to zero"):
            bfs.compute_average_winding(zm, poly, bad)

    def test_compute_average_winding_skips_none_hard_charge(
            self, poly_A, monkeypatch):
        """A hard boundary (charge=None) disables the conservation check and
        must never enter arithmetic: the winding path stays TypeError-free."""
        poly = CharPoly(*poly_A)
        zm = bfs.Mu2MidZM(poly, 1.0 + 0j, 0.1)
        zm.run()
        zm.build_mu2_mid()

        quad_calls = []

        def fake_quad(m, poly_, E_ref, mu1, theta2):
            quad_calls.append(theta2)
            return 1.0

        monkeypatch.setattr(sgbz_winding, "_loop_winding_quad", fake_quad)
        charges = [
            dict(theta1=0.0, theta2=0.5, charge=None, kind="mr"),
            dict(theta1=1.0, theta2=2.5, charge=1, kind="ordinary"),
        ]
        W = bfs.compute_average_winding(zm, poly, charges)
        assert np.isfinite(W)
        assert len(quad_calls) >= 1
