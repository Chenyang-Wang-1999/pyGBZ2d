"""Tests for debug_tool — fixed-(E_ref, mu1) GBZ subsets + winding loops.

The HN-2D chain pair (conftest's shared model) is fast (<1 s per method) and
has analytic answers: mu1 = gamma1 = 0.2 is the true GBZ radius (continuum
case), while any other mu1 (e.g. 0.25) gives discrete crossings whose
charges must sum to zero and reproduce the loop-winding profile jumps.

The Haldane gain-loss supercell scenario (the E=1.212 debug point that
motivated the tool) is marked slow.
"""

import math

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from conftest import build_HN2D_polynomial

from pygbz2d.core import CharPoly, PointSubset, LineSubset, TWO_PI

from debug_tool import (
    MethodDebug, GBZDebugReport, LoopWinding, GapCheck, MethodLoops,
    collect_debug_subsets, auto_theta2_grid, compute_loop_windings,
    plot_winding_debug,
)


GAMMA1, GAMMA2 = 0.2, 0.3


@pytest.fixture(scope="module")
def poly_hn2d():
    coeffs, degs = build_HN2D_polynomial(1.0, 1.0, GAMMA1, GAMMA2, 0.0, 0.0)
    return CharPoly(coeffs, degs)


# ---------------------------------------------------------------------------
# Function 1: collect_debug_subsets
# ---------------------------------------------------------------------------

class TestCollectDebugSubsets:
    def test_report_structure(self, poly_hn2d):
        report = collect_debug_subsets(poly_hn2d, 1.0 + 0j, 0.25)
        assert isinstance(report, GBZDebugReport)
        assert set(report.methods) == {"sgbz", "amoeba"}
        assert report.sgbz is report["sgbz"]
        assert report.amoeba is report["amoeba"]
        for md in report.methods.values():
            assert md.success, md.error
            assert md.mu1 == 0.25
            # subsets are the unified ConnectedSubset types
            for s in md.subsets:
                assert isinstance(s, (PointSubset, LineSubset))

    def test_sgbz_charges_aligned_with_points(self, poly_hn2d):
        """SGBZ charges: one per PointSubset, same order, ordinary ±1."""
        md = collect_debug_subsets(poly_hn2d, 1.0 + 0j, 0.25)["sgbz"]
        assert not md.is_continuum
        n_pts = len(md.point_subsets)
        assert n_pts == len(md.charges) >= 2
        for sub, ch in zip(md.subsets, md.charges):
            assert isinstance(sub, PointSubset)
            assert ch["kind"] == "ordinary"
            assert ch["charge"] in (+1, -1)
            assert ch["theta1"] == pytest.approx(sub.theta1, abs=1e-6)
            assert ch["theta2"] == pytest.approx(
                np.angle(sub.beta2) % TWO_PI, abs=1e-6)
        # one full theta2 circle returns the winding to itself
        assert sum(c["charge"] for c in md.charges) == 0
        # W(E, mu1) != 0 here: 0.25 is off the GBZ radius (gamma1 = 0.2)
        assert md.W == pytest.approx(0.667126, abs=1e-3)

    def test_amoeba_slice_has_mu2_and_w2_zero(self, poly_hn2d):
        md = collect_debug_subsets(poly_hn2d, 1.0 + 0j, 0.25)["amoeba"]
        assert md.mu2 is not None
        assert md.charges == []          # amoeba has no charge concept
        assert md.w_secondary == pytest.approx(0.0, abs=1e-6)
        # w1 (the outer-loop winding) agrees with the SGBZ W on this model
        assert md.W == pytest.approx(0.667126, abs=1e-3)

    def test_continuum_case_at_true_gbz_radius(self, poly_hn2d):
        """mu1 = gamma1 is the analytic GBZ radius: SGBZ reports a continuum."""
        report = collect_debug_subsets(poly_hn2d, 1.0 + 0j, GAMMA1)
        md = report["sgbz"]
        assert md.is_continuum
        assert md.charges == []
        assert md.W is None
        assert md.line_subsets, "continuum must materialize LineSubsets"
        for line in md.line_subsets:
            assert line.mu1 == pytest.approx(GAMMA1)
        # the amoeba lands on the analytic mu2 = gamma2 band
        assert report["amoeba"].mu2 == pytest.approx(GAMMA2, abs=1e-4)

    def test_method_failure_is_captured_not_raised(self, poly_hn2d):
        """A broken method reports error=... instead of killing the report."""
        import debug_tool.gbz_debug as gd
        original = gd._collect_amoeba
        try:
            gd._collect_amoeba = lambda *a, **k: (_ for _ in ()).throw(
                RuntimeError("boom"))
            report = collect_debug_subsets(poly_hn2d, 1.0 + 0j, 0.25)
        finally:
            gd._collect_amoeba = original
        assert report["amoeba"].error and "boom" in report["amoeba"].error
        assert not report["amoeba"].success
        assert report["sgbz"].success, "the other method's evidence survives"
        with pytest.raises(ValueError):
            compute_loop_windings(report["amoeba"], [0.5])


# ---------------------------------------------------------------------------
# Function 2: loop windings + gap checks
# ---------------------------------------------------------------------------

class TestComputeLoopWindings:
    def test_sgbz_loop_windings_and_charge_consistency(self, poly_hn2d):
        """The profile jump between consecutive loops equals the enclosed
        charges — the core invariant the tool exists to check."""
        md = collect_debug_subsets(poly_hn2d, 1.0 + 0j, 0.25)["sgbz"]
        ml = compute_loop_windings(md, [0.5, 2.5, 4.5])
        assert isinstance(ml, MethodLoops)
        assert [lw.theta2 for lw in ml.loops] == [0.5, 2.5, 4.5]
        for lw in ml.loops:
            assert isinstance(lw, LoopWinding)
            assert lw.reliable
            assert lw.winding == round(lw.winding_raw)
            assert lw.min_abs_f > 0
        # charges: -1 at theta2≈2.096, +1 at theta2≈4.187.  u(theta2):
        #   gap (0.5→2.5) encloses the -1 at 2.096  → jump -1
        #   gap (2.5→4.5) encloses the +1 at 4.187  → jump +1
        #   gap (4.5→0.5+2π) encloses no crossing    → jump 0
        w = {lw.theta2: lw.winding for lw in ml.loops}
        assert w[2.5] - w[0.5] == -1
        assert w[4.5] - w[2.5] == +1
        assert w[0.5] - w[4.5] == 0
        assert all(g.status == "ok" for g in ml.gap_checks)
        for g in ml.gap_checks:
            assert isinstance(g, GapCheck)
            assert g.expected == g.actual
        # the charge-propagated profile reconstruction reproduces the
        # solver's W exactly, even though this coarse grid's gaps each
        # enclose a full cancelling +1/-1 charge pair
        assert ml.W_profile == pytest.approx(md.W, abs=1e-9)

    def test_amoeba_root_count_matches_quad(self, poly_hn2d):
        md = collect_debug_subsets(poly_hn2d, 1.0 + 0j, 0.25)["amoeba"]
        ml = compute_loop_windings(md, [0.5, 2.5, 4.5])
        for lw in ml.loops:
            assert lw.winding_quad == pytest.approx(lw.winding, abs=0.1)
            assert lw.margin > 1e-6
        # discrete points: the SIGNED point-jump sum across a gap must equal
        # the winding jump (cancelling +1/−1 pairs inside one gap are legal)
        assert all(g.status == "ok" for g in ml.gap_checks)
        # amoeba's own w1 is the same arc-weighted mean
        assert ml.W_profile == pytest.approx(md.W, abs=1e-9)

    def test_loop_through_crossing_is_flagged(self, poly_hn2d):
        """A loop landing exactly on a crossing theta2 must be flagged
        unreliable rather than silently returning a winding."""
        md = collect_debug_subsets(poly_hn2d, 1.0 + 0j, 0.25)["sgbz"]
        t2_cross = md.charges[0]["theta2"]
        ml = compute_loop_windings(md, [t2_cross])
        assert not ml.loops[0].reliable
        assert ml.loops[0].min_abs_f < 1e-8

    def test_continuum_gap_checks_skipped(self, poly_hn2d):
        md = collect_debug_subsets(poly_hn2d, 1.0 + 0j, GAMMA1)["sgbz"]
        ml = compute_loop_windings(md, [0.5, 3.0])
        assert all(lw.note for lw in ml.loops)      # continuum caveat shown
        assert ml.gap_checks[0].status == "skipped"


# ---------------------------------------------------------------------------
# Auto theta2 grid
# ---------------------------------------------------------------------------

class TestAutoTheta2Grid:
    def test_midpoints_avoid_crossings(self, poly_hn2d):
        md = collect_debug_subsets(poly_hn2d, 1.0 + 0j, 0.25)["sgbz"]
        grid = auto_theta2_grid(md)
        assert 1 <= len(grid) <= 8
        occupied = [c["theta2"] for c in md.charges]
        for t2 in grid:
            assert min(abs(t2 - o) for o in occupied) > 0.1

    def test_no_subsets_falls_back_to_quadrants(self, poly_hn2d):
        md = MethodDebug(method="sgbz", E_ref=1j, mu1=0.0, poly=poly_hn2d)
        assert auto_theta2_grid(md) == [
            pytest.approx(a * math.pi / 4) for a in (1, 3, 5, 7)]


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

class TestPlotWindingDebug:
    def test_report_plot_produces_artists(self, poly_hn2d):
        report = collect_debug_subsets(poly_hn2d, 1.0 + 0j, 0.25)
        fig, ax, results = plot_winding_debug(report, [0.5, 2.5, 4.5])
        try:
            assert set(results) == {"sgbz", "amoeba"}
            # subsets + loops + MR guides are all on the axes
            assert len(ax.lines) > 4
            texts = " ".join(t.get_text() for t in ax.texts)
            assert "+1" in texts and "-1" in texts      # charge annotations
            assert "w=" in texts                        # loop winding labels
            assert ax.get_xlabel() == r"$\theta_1$"
            # torus frame with a small visibility margin beyond [0, 2π]
            xl = ax.get_xlim()
            assert -0.5 < xl[0] <= 0 and TWO_PI <= xl[1] < TWO_PI + 0.5
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)

    def test_single_method_and_auto_grid(self, poly_hn2d):
        md = collect_debug_subsets(poly_hn2d, 1.0 + 0j, 0.25, methods=("sgbz",))
        fig, ax, results = plot_winding_debug(md)   # auto theta2
        try:
            assert set(results) == {"sgbz"}
            assert results["sgbz"].loops
        finally:
            import matplotlib.pyplot as plt
            plt.close(fig)


# ---------------------------------------------------------------------------
# Slow: the Haldane gain-loss debug point that motivated the tool
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestHaldaneDebugPoint:
    """E_ref = 1.212, mu1 = 0.135328598265 (log/2026-08-16, simple-plot.py).

    14 SGBZ crossings with ±1 ordinary charges; every charge-consistency
    gap check must pass — a regression here is exactly the class of bug the
    tool is meant to catch.
    """

    def test_full_session(self):
        import importlib.util
        from pathlib import Path
        root = Path(__file__).resolve().parents[1]
        spec = importlib.util.spec_from_file_location(
            "haldane", root / "playground" / "Haldane-model-gainloss.py")
        h = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(h)
        model = h.Haldane_non_Hermitian_phase(*h.ALL_PARAMS)
        model = model.get_supercell(
            [(0, 0), (1, 0)], np.array([[1, 1], [-1, 1]], dtype=int))
        coeffs, degs = model.get_characteristic_polynomial_data()
        poly = CharPoly(coeffs, degs)

        report = collect_debug_subsets(poly, 1.2120000000000002 + 0j,
                                       0.135328598265)
        assert report["sgbz"].success, report["sgbz"].error
        assert len(report["sgbz"].point_subsets) == 14
        assert all(c["kind"] == "ordinary" for c in report["sgbz"].charges)
        assert sum(c["charge"] for c in report["sgbz"].charges) == 0

        grid = auto_theta2_grid(report["sgbz"])
        ml = compute_loop_windings(report["sgbz"], grid)
        assert all(g.status == "ok" for g in ml.gap_checks)
        assert ml.W_profile == pytest.approx(report["sgbz"].W, abs=1e-9)
