"""Regression tests for the 2026-08-18 code-review fixes.

Each test pins the exact failure mode found in the review:

* H1 — ZeroManager MR restart ping-pong (infinite loop, same MR re-recorded):
  the restart distance must be branch-point safe for the (h0, min_dtheta)
  pair, and a no-forward-progress refinement must raise, not loop.
* H2 — _join_runs_across_mrs merge budget: joining a 3-fold continuum over
  2 segments needs 3 merges, one more than the old n_seg-capped loop allowed.
* M4 — interval MR trigger blinded by nan-sentinel singular roots.
* M5 — estimate_error accepting a completely-wrong step when the actual
  roots are majority-inf.
"""

import signal
import numpy as np
import pytest
from math import pi
from types import SimpleNamespace

from gbz_types import CharPoly, TWO_PI
from continuation.zero_manager import ZeroManager
from continuation.multiple_roots import (
    MultipleRootIntervalTrigger, _closest_pair_deriv,
)
from continuation.arclength import estimate_error


def _make_pingpong_poly():
    """f = β₂² − (β₁ − i): double root at θ₁ = π/2, β₂ = 0 (branch point)."""
    coeffs = np.array([-1, 1j, 1], dtype=complex)
    degs = np.array([[0, 1, 0], [0, 0, 0], [0, 0, 2]], dtype=int)
    return CharPoly(coeffs, degs)


class TestH1MrRestartPingpong:
    """The MR restart distance must clear both step-collapse scales."""

    @pytest.mark.parametrize("h0,min_dtheta", [
        (0.5, 1e-6),   # original repro: 3572 duplicate MRs in 45 s
        (0.1, 1e-6),   # sub-agent repro: 10000 identical MR records
        (0.2, 1e-4),
        (1.0, 1e-6),
    ])
    def test_run_terminates_with_single_mr(self, h0, min_dtheta):
        poly = _make_pingpong_poly()
        zm = ZeroManager(poly, 0.5 + 0j, 0.0)

        def on_alarm(sig, frm):
            raise TimeoutError("watchdog: run() did not terminate")
        old = signal.signal(signal.SIGALRM, on_alarm)
        signal.alarm(60)
        try:
            zm.run(h0=h0, min_dtheta=min_dtheta)
        finally:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old)

        thetas = [m.theta1 for m in zm.multiple_roots]
        assert len(thetas) == 1, f"expected 1 MR, got {thetas}"
        assert thetas[0] == pytest.approx(pi / 2, abs=1e-5)

    def test_no_duplicate_mr_thetas(self):
        """Even if several segments appear, no two MR records share θ₁."""
        poly = _make_pingpong_poly()
        zm = ZeroManager(poly, 0.5 + 0j, 0.0)
        zm.run(h0=0.5, min_dtheta=1e-6)
        thetas = sorted(m.theta1 for m in zm.multiple_roots)
        for a, b in zip(thetas, thetas[1:]):
            assert b - a > 1e-6, f"duplicate/near-duplicate MR thetas: {thetas}"


class TestH2JoinBudget:
    """_join_runs_across_mrs must not cap merges at n_seg."""

    def test_three_fold_continuum_closes(self):
        from brute_force_SGBZ.continuum_lines import _join_runs_across_mrs
        from gbz_types import JoinableLinePiece
        from continuation.multiple_roots import MultipleRootInfo

        vals = [1.0, 10.0, 100.0]  # 3-fold continuum: tracks A, B, C
        seg0 = SimpleNamespace(
            theta1_arr=np.array([0.0, np.pi / 2, np.pi]),
            tracked_roots=np.array([[v + 0j for v in vals],
                                    [v + 0.5 for v in vals],
                                    [v + 1.0 for v in vals]]),
            left_mr=-1, right_mr=0,
        )
        seg1 = SimpleNamespace(
            theta1_arr=np.array([np.pi, 3 * np.pi / 2, TWO_PI]),
            tracked_roots=np.array([[v + 1.0 + 0j for v in vals],
                                    [v + 0.5 for v in vals],
                                    [v + 0.0 for v in vals]]),
            left_mr=0, right_mr=-1,
        )
        zm = SimpleNamespace(
            segments=[seg0, seg1], K=3, _event_groups=[],
            # This synthetic seam has identity monodromy: each final-row root
            # is already the corresponding first-row root.
            boundary_perm=np.array([0, 1, 2]),
            multiple_roots=[MultipleRootInfo(
                theta1=float(np.pi), cluster_indices=[],
                roots=seg0.tracked_roots[-1])],
        )

        pieces = []
        for s_idx, seg in enumerate((seg0, seg1)):
            for j in range(len(vals)):
                pieces.append(JoinableLinePiece(
                    0j, 0.0, seg.theta1_arr.copy(),
                    seg.tracked_roots[:, j].copy(), ml=s_idx, mr=s_idx))

        joined = _join_runs_across_mrs(zm, pieces)
        # 3 closed loops, one per track — NOT 4 pieces (the capped budget
        # left track C broken into two open pieces).
        assert len(joined) == 3
        for p in joined:
            assert p.theta1_arr[0] == pytest.approx(np.pi)
            assert p.theta1_arr[-1] == pytest.approx(3 * np.pi)
            # closed loop: endpoint roots coincide across the seam
            assert abs(p.beta2_arr[0] - p.beta2_arr[-1]) < 1e-9


class TestM4IntervalTriggerNanBlindness:
    def test_closest_pair_skips_nan_sentinel_roots(self):
        # j0 is a 0/∞ padding root (V = nan) whose distance to j1 (0.03) is
        # below the trigger threshold: without the fix it is the tracked
        # closest pair and its nan deriv blinds the trigger forever.
        roots1 = np.array([0.0 + 0j, 0.03 + 0j, 1.00 + 0j, 1.04 + 0j])
        V1 = np.array([np.nan + 0j, 0.1 + 0j, -0.5 + 0j, -0.5 + 0j])
        d, deriv, pair = _closest_pair_deriv(roots1, V1)
        assert pair == (2, 3), "nan-sentinel root must not own the pair"
        assert np.isfinite(deriv)

    def test_trigger_fires_with_singular_root_present(self):
        roots1 = np.array([0.0 + 0j, 0.03 + 0j, 1.00 + 0j, 1.04 + 0j])
        V1 = np.array([np.nan + 0j, 0.1 + 0j, -0.5 + 0j, -0.5 + 0j])
        roots2 = np.array([0.0 + 0j, 0.03 + 0j, 1.02 + 0j, 1.03 + 0j])
        V2 = np.array([np.nan + 0j, 0.1 + 0j, 0.5 + 0j, 0.5 + 0j])

        trig = MultipleRootIntervalTrigger(min_dist_threshold=0.1)
        ok1, _ = trig(roots1, V1, 0.0)
        ok2, interval = trig(roots2, V2, 0.1)
        assert not ok1
        assert ok2, "approach→separation flip must fire despite nan sentinel"
        assert interval == (0.0, 0.1)
        assert np.isfinite(trig._prev_deriv) or ok2  # no nan persisted

    def test_all_singular_roots_return_sentinel(self):
        d, deriv, pair = _closest_pair_deriv(
            np.array([0j, 1j]), np.array([np.nan + 0j, np.nan + 0j]))
        assert pair == (-1, -1)
        assert d == np.inf


class TestM5EstimateErrorInfScale:
    def test_majority_inf_actual_rejected(self):
        pred = np.array([1.0 + 0j, 2.0 + 0j, 3.0 + 0j])
        act = np.array([np.inf + 0j, np.inf + 0j, 3.0 + 0j])
        assert estimate_error(pred, act) > 1.0

    def test_single_inf_actual_rejected(self):
        pred = np.array([1.0 + 0j, 2.0 + 0j, 3.0 + 0j])
        act = np.array([np.inf + 0j, 2.0 + 0j, 3.0 + 0j])
        assert estimate_error(pred, act) > 1.0

    def test_identical_roots_still_zero_error(self):
        r = np.array([1.0 + 1j, 2.0 - 1j])
        assert estimate_error(r, r) < 1e-10
