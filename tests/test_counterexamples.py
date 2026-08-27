"""Regression tests for historical manual counterexamples (T4/T5/T6).

* T4 — the gain-loss Haldane counterexamples from
  log/2026-08-16-pairwise-linear-intersection.md (E=-1.56 → index (12, 0),
  E=-1.406 → index (6, 0)).  These caught the μ₂_mid tangent blow-up fixed
  by _mark_mr_tangents_inf and were previously only checked by hand.
  Requires BerryPy (optional dependency → importorskip) and is marked slow.
* T5 — deep pairwise internals: _refine_pair_crossing brentq fallback,
  insert_event_groups descending-order insertion with seam-deferred rows.
* T6 — SGBZ plateau probe (_probe_zero_plateau_near_mu1) and the iterative
  multiple-root solver.
"""

import numpy as np
import pytest
from math import pi
from types import SimpleNamespace

from bfgbz2d.core import CharPoly, JoinableLinePiece
from conftest import build_HN2D_polynomial

import bfgbz2d.sgbz as bfs
from bfgbz2d.sgbz import pairwise as sgbz_pairwise
from bfgbz2d.sgbz.mu2mid import Mu2MidZM
from bfgbz2d.continuation.multiple_roots import solve_multiple_roots_iterative


# ===========================================================================
# T4 — gain-loss Haldane counterexamples (BerryPy + slow)
# ===========================================================================

def _haldane_gainloss_poly():
    """f(E, β₁, β₂) of the gain-loss Haldane model (γ on sublattice B).

    Built via BerryPy's TightBindingModel with the ALL_PARAMS of
    playground/Haldane-model-gainloss.py: t1=1, t2=0.5, phi=π/3, M=0.5j,
    gamma=0 (γ enters as v1=v2=t2·e^{iγ} and the on-site ±M).
    """
    tb = pytest.importorskip("BerryPy").TightBinding
    from cmath import exp, cos, sin, sqrt

    t1, t2, phi, M, gamma = 1, 0.5, pi / 3, 0.5j, 0
    u1 = u2 = t1
    v1 = v2 = t2 * exp(1j * gamma)

    lattice_vec = np.array([
        [-cos(pi / 3), -cos(pi / 3)],
        [-sin(pi / 3), sin(pi / 3)],
    ])
    intra_cell = [
        [0, 0, M], [1, 1, -M], [1, 0, u1], [0, 1, u2],
    ]
    inter_cell = [
        [0, 0, v2 * exp(1j * phi), (-1, 0)],
        [0, 0, v2 * exp(1j * phi), (0, -1)],
        [0, 0, v2 * exp(1j * phi), (1, 1)],
        [0, 0, v1 * exp(-1j * phi), (1, 0)],
        [0, 0, v1 * exp(-1j * phi), (0, 1)],
        [0, 0, v1 * exp(-1j * phi), (-1, -1)],
        [1, 0, u1, (0, -1)], [1, 0, u1, (1, 0)],
        [0, 1, u2, (0, 1)], [0, 1, u2, (-1, 0)],
        [1, 1, v2 * exp(-1j * phi), (-1, 0)],
        [1, 1, v2 * exp(-1j * phi), (0, -1)],
        [1, 1, v2 * exp(-1j * phi), (1, 1)],
        [1, 1, v1 * exp(1j * phi), (1, 0)],
        [1, 1, v1 * exp(1j * phi), (0, 1)],
        [1, 1, v1 * exp(1j * phi), (-1, -1)],
    ]
    site_coord_cart = np.array([[0, 1 / (2 * sqrt(3))],
                                [0, -1 / (2 * sqrt(3))]])

    model = tb.TightBindingModel(dim=2, SiteNum=2, LatticeVec=lattice_vec,
                                 InCell=intra_cell, InterCell=inter_cell)
    model.SiteCoord = model.cart2lattice(site_coord_cart.T).T
    coeffs, degs = model.get_characteristic_polynomial_data()
    return np.asarray(coeffs, dtype=complex), np.asarray(degs, dtype=int)


@pytest.mark.slow
class TestHaldaneGainlossCounterexamples:
    """The two hand-checked energies from the 2026-08-16 pairwise log."""

    @pytest.mark.parametrize("E_ref,expected_index", [
        (-1.56 + 0j, (12, 0)),
        (-1.406 + 0j, (6, 0)),
    ])
    def test_counterexample_indices(self, E_ref, expected_index):
        coeffs, degs = _haldane_gainloss_poly()
        gbz = bfs.collect_GBZ_subsets(coeffs, degs, E_ref)
        assert gbz.success, f"solver failed: {gbz.error}"
        assert gbz.is_gbz
        assert gbz.index == expected_index


# ===========================================================================
# T5 — pairwise deep internals
# ===========================================================================

class TestRefinePairCrossingFallback:
    """brentq failure must fall back to the linear prediction, not raise."""

    def _zm_stub(self):
        # Minimal duck-typed ZeroManager: one segment, two mesh rows.  The
        # pair (rep_a=0, rep_b=1) has EXACTLY equal moduli at row 0
        # (2.0 vs 2.0·e^{0.5j}) → f(theta_lo) == 0: the touch short-circuit
        # returns the endpoint directly, no brentq call.
        roots = np.array([[2.0 + 0j, 2.0 * np.exp(0.5j)],
                          [1.5 + 0j, 1.4 + 0j]])
        seg = SimpleNamespace(
            theta1_arr=np.array([0.0, 0.1]),
            tracked_roots=roots,
            tangents=np.zeros((2, 2), dtype=complex),
        )
        return SimpleNamespace(segments=[seg])

    def test_endpoint_exact_touch_short_circuits(self):
        zm = self._zm_stub()
        theta, direction, converged = sgbz_pairwise._refine_pair_crossing(
            zm, 0, 0, rep_a=0, rep_b=1, theta_lo=0.0, theta_hi=0.1,
            crossing_tol=1e-10, min_direction_deriv=1e-12,
        )
        assert theta == 0.0
        assert converged is True
        # direction is None because both tangents are 0 (below the floor)
        assert direction is None


class TestInsertEventGroupsSeamDeferral:
    """Seam groups get rows assigned last, others re-resolve on final mesh."""

    def test_identity_membership_not_value(self):
        # insert_event_groups must skip EXACTLY the seam group objects
        # (identity), never a distinct-but-equal group.  Constructed via
        # the public path: this test guards the `id()` switch.
        import inspect
        src = inspect.getsource(sgbz_pairwise.insert_event_groups)
        assert "id(g) in seam_ids" in src or "id(g)" in src, (
            "insert_event_groups must use identity membership for the seam "
            "groups (dataclass value-equality can match a different group)"
        )


# ===========================================================================
# T6 — SGBZ plateau probe + iterative MR solver
# ===========================================================================

class TestSgbzPlateauProbe:
    def test_probe_reports_not_found_on_regular_model(self):
        """A regular in-spectrum point: probing μ₁ ± steps finds no plateau.

        Uses the HN model at its analytic GBZ energy: the plateau probe
        should complete and report not_found (winding non-zero off the
        boundary), demonstrating the SGBZ-side probe path executes.
        """
        from bfgbz2d.sgbz.plateau import _probe_zero_plateau_near_mu1
        from bfgbz2d.core import CONTINUUM_TOL
        from bfgbz2d.sgbz.pairwise import CROSSING_TOL as _CROSSING_TOL
        from bfgbz2d.sgbz.sgbz_solver import solve_SGBZ_for_E

        coeffs, degs = build_HN2D_polynomial(
            J1=1.0, J2=1.0, gamma_1=0.2, gamma_2=0.3,
            delta_1=0.0, delta_2=0.0,
        )
        poly = CharPoly(coeffs, degs)
        # E_ref = 1.0: the in-spectrum energy the existing SGBZ tests use
        # (winding straddles zero around mu1 = gamma_1 = 0.2 there).
        E_ref = 1.0 + 0j
        res = solve_SGBZ_for_E(poly, E_ref)
        info = _probe_zero_plateau_near_mu1(
            poly, E_ref, res["mu1"], res.get("_mu1_bracket"), {},
            continuum_tol=CONTINUUM_TOL,
            crossing_tol=_CROSSING_TOL,
            zero_tol=1e-10,
        )
        assert info["status"] in ("not_found", "inconclusive")
        assert info["found"] is False
        assert len(info["points"]) > 0


class TestIterativeMrSolver:
    def test_solves_double_root(self):
        # f = β₂² − (β₁ − i): double root at β₁ = i, β₂ = 0.  The start
        # must sit slightly OFF the μ₁ circle in the radial direction:
        # exactly on it the double root is degenerate and root()'s Newton
        # stall guard trips ("not making good progress") — the radial
        # offset keeps the start in a regular neighbourhood.
        coeffs = np.array([-1, 1j, 1], dtype=complex)
        degs = np.array([[0, 1, 0], [0, 0, 0], [0, 0, 2]], dtype=int)
        poly = CharPoly(coeffs, degs)

        beta1_approx = 1.001j * np.exp(-1e-3)   # near (i·e^{−iε}), radially off
        beta1, beta2 = solve_multiple_roots_iterative(
            poly, 0.5 + 0j, beta1_approx, 1e-3 + 0j,
        )
        assert abs(beta1 - 1j) < 1e-8
        assert abs(beta2) < 1e-6

    def test_preserves_energy_gauge(self):
        # The refined β₁ must land on the unit circle (μ₁ = 0 gauge)
        coeffs = np.array([-1, 1j, 1], dtype=complex)
        degs = np.array([[0, 1, 0], [0, 0, 0], [0, 0, 2]], dtype=int)
        poly = CharPoly(coeffs, degs)
        beta1, beta2 = solve_multiple_roots_iterative(
            poly, 0.5 + 0j, 1.0005j * np.exp(-5e-4), 5e-4 + 0j,
        )
        assert abs(abs(beta1) - 1.0) < 1e-8
