'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-10
Copyright © Department of Physics, Tsinghua University. All rights reserved

Topological charge + average major-axis winding (§3, §6.4).

Computes the average major-axis winding number ``W(E_ref, mu1)`` — the
quantity whose zero in ``mu1`` defines the SGBZ — from a
``Mu2MidZM`` (a ``ZeroManager`` with the piecewise-smooth μ₂_mid built) at
fixed ``(E_ref, mu1)`` together with the charge list from
:mod:`brute_force_SGBZ.crossings`.

Design (``log/2026-08-13-SGBZ算法梳理.md`` §3/§6.4):

  * The loop is ``β₂ = exp(μ₂_mid(θ₁) + iθ₂)`` with ``θ₁ ∈ [0, 2π)`` and
    fixed ``θ₂``; ``μ₂_mid`` is the piecewise-smooth boundary-pair mean — a
    first-class object with analytic values and derivatives and left/right
    derivatives at its breakpoints (§2.0).  The loop winding is the integral
    of ``Im[f'/f]`` over ``θ₁``, split at the μ₂_mid breakpoints so each quad
    segment is smooth (§6.4).  This replaces the old mesh-row argument-change
    summation: the path is the *exact* piecewise-smooth μ₂_mid, not a polyline
    approximation, and the integration honours the analytic derivative.
  * Crossings come in two kinds.  **Ordinary** (charge ±1, SOFT): the winding
    across it is fixed by its charge, so it only partitions ``θ₂`` into
    intervals *within* a region.  **MR / tangent / unknown** (charge 0, HARD):
    charge unknown, so it DELIMITES regions; across a hard boundary the winding
    is recomputed independently.
  * Therefore: partition the circle into regions delimited by hard
    boundaries; within each region, ordinary boundaries split it into
    intervals.  Pick ONE seed interval per region (the safest — farthest from
    all roots), compute ``w₀`` there, and propagate across the region's soft
    boundaries via charges.  Each region contributes exactly one loop-winding
    evaluation; ±1 numerical noise is confined to the per-region seed.
'''

from __future__ import annotations

import math

import numpy as np
from cmath import exp
from scipy import integrate

from gbz_types import CharPoly, PointSubset
from continuation import ZeroManager

from .crossings import detect_crossings_simple, _ensure_mu2mid, _MAX_NEWTON_ITER
from .mu2mid import Mu2MidZM, _Mu2MidPath


# ---------------------------------------------------------------------------
# WindingFun + get_winding_number (quad of Im[f'/f])
# ---------------------------------------------------------------------------

class WindingFun:
    """Winding integrand of the characteristic polynomial along a loop.

    Calling the instance at parameter t returns ``Im[f'(t)/f(t)]``, i.e.
    ``d/dt Im log f(param(t))``; its integral over the loop range divided by
    ``2π`` is the winding number of f around zero (see
    :func:`get_winding_number`).
    """

    def __init__(self, char_poly: CharPoly, loop_fun, loop_range: tuple[float, float]):
        """
        Parameters:
            char_poly: characteristic Laurent polynomial.
            loop_fun: ``t -> (param, dparam/dt)`` where ``param`` is the
                tuple of polynomial variables at t.
            loop_range: ``(t_start, t_end)`` parameter range of the closed loop.
        """
        self.char_poly = char_poly
        self.loop_fun = loop_fun
        self.loop_range = loop_range

    def __call__(self, t: float) -> float:
        param, dparam_dt = self.loop_fun(t)
        val = self.char_poly.eval_val(param)
        partials = self.char_poly.eval_partials(param)
        dval_dt = sum(partials[k] * dparam_dt[k] for k in range(len(partials)))
        return (dval_dt / val).imag


def get_winding_number(
    winding_fun: WindingFun,
    seg_bounds: list[float] | None = None,
    *,
    n_seg: int = 1,
) -> float:
    """Integrate ``winding_fun`` over its loop range, divide by ``2π``.

    If *seg_bounds* is given (a sorted list of θ₁ breakpoints, endpoints
    excluded), each interval between consecutive breakpoints — plus the
    ``[start, first_bp]`` and ``[last_bp, end]`` ends — is handed to
    ``scipy.integrate.quad`` separately.  Splitting at the μ₂_mid breakpoints
    keeps every quad segment on a single smooth piece (§6.4).  Otherwise the
    range is split into *n_seg* equal pieces.

    Returns the real-valued winding number (unrounded; callers round).
    """
    a, b = winding_fun.loop_range
    if seg_bounds:
        bounds = [a] + sorted(set(seg_bounds)) + [b]
    else:
        bounds = list(np.linspace(a, b, n_seg + 1))
    total = 0.0
    # epsabs/epsrel = 1e-3 looks loose, but every caller rounds the result
    # to an integer: the total absolute error is ≲ (#intervals)·1e-3, i.e.
    # ~1e-2 per loop → ~1.6e-3 in winding units, two orders below the 0.5
    # rounding margin.  Tighter tolerances only buy speed loss.
    for i in range(len(bounds) - 1):
        total += integrate.quad(
            winding_fun, bounds[i], bounds[i + 1],
            epsabs=1e-3, epsrel=1e-3, limit=200,
        )[0]
    return total / (2 * math.pi)


# ---------------------------------------------------------------------------
# Piecewise-smooth μ₂_mid loop path
# ---------------------------------------------------------------------------
# _Mu2MidPath now lives in mu2mid.py (shared with crossings.py, which needs it
# to evaluate the built μ₂_mid curve at arbitrary probe θ without mutating the
# mesh — the 2026-08-14 zero-mutation crossing refinement).


def _loop_winding_quad(
    zm: Mu2MidZM,
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    theta2: float,
) -> float:
    """Loop winding number via quad of ``Im[f'/f]`` over the μ₂_mid loop.

    The loop is ``β₂ = exp(μ₂_mid(θ₁) + iθ₂)``, ``θ₁ ∈ [0, 2π)``, with the
    analytic ``dβ₂/dθ₁ = β₂ · μ₂_mid'(θ₁)`` from the piecewise-smooth path
    (§6.4).  The integration is split at the μ₂_mid breakpoints so every quad
    segment lies on one smooth piece.  Consecutive segments share their
    MR-boundary row exactly, so the cross-segment seam contributes nothing
    extra (the breakpoint split already isolates it).
    """
    path = _Mu2MidPath(zm)
    twopi = 2 * math.pi
    bp_thetas = [float(b.theta1) for b in zm.mu2_mid_breakpoints
                 if 0.0 < b.theta1 < twopi]

    def loop_fun(t: float):
        mu2v, dmu2 = path.value_deriv(t)
        beta1 = exp(mu1 + 1j * t)
        beta2 = exp(mu2v + 1j * theta2)
        dbeta1 = 1j * beta1
        dbeta2 = beta2 * dmu2
        return (E_ref, beta1, beta2), (0.0, dbeta1, dbeta2)

    wf = WindingFun(poly, loop_fun, (0.0, twopi))
    return get_winding_number(wf, seg_bounds=bp_thetas)


# ---------------------------------------------------------------------------
# Seed selection (safety metric: min |f| along the loop)
# ---------------------------------------------------------------------------

def _loop_min_f(
    theta2: float,
    zm: Mu2MidZM,
    poly: CharPoly,
) -> float:
    """Minimum ``|f(E, β₁, β₂)|`` along the μ₂_mid loop at fixed *theta2*.

    The TRUE safety metric for a loop-winding seed (§3.2 "离零点最远" =
    farthest from the zeros of f, not from the β₂ roots).  The loop is
    ``β₁=exp(μ₁+iθ₁), β₂=exp(μ₂_mid(θ₁)+iθ₂)``; the quad of ``Im[f'/f]``
    is reliable only while the loop stays clear of char-poly zeros, and a
    zero of f is reached when ``β₂_loop`` meets a tracked root — but
    ``|f|`` also depends on β₁ and on ``|∂f/∂β₂|``, so the β₂-distance
    proxy underestimates the danger (a loop 0.5 from a β₂ root can still
    have ``|f|≈0`` where ``|∂f/∂β₂|`` is large).  Maximising ``min |f|``
    picks a θ₂ where the whole loop is genuinely far from any zero.
    """
    worst = math.inf
    for s_idx, seg in enumerate(zm.segments):
        th = seg.theta1_arr
        mu2 = zm.seg_mu2_values[s_idx]
        beta1 = np.exp(zm.mu1 + 1j * th)
        beta2 = np.exp(mu2 + 1j * theta2)
        for i in range(len(th)):
            v = abs(poly.eval_val((zm.E_ref, beta1[i], beta2[i])))
            if v < worst:
                worst = v
                if worst == 0.0:
                    return 0.0
    return float(worst) if math.isfinite(worst) else 0.0


def _pick_seed_theta2(
    intervals: list[tuple[float, float]],
    zm: Mu2MidZM,
    poly: CharPoly,
    *,
    n_per_interval: int = 4,
) -> tuple[float, int]:
    """Pick the θ₂ maximizing ``min |f|`` along the loop over *intervals*.

    The loop-winding seed must lie inside one of the region's intervals and
    stay clear of char-poly zeros for the quad to be reliable, so the safest
    θ₂ maximises ``min |f(E, β₁(θ₁), β₂_loop(θ₁))|`` over θ₁ (§3.2 — the
    true safety metric, accounting for β₁ and |∂f/∂β₂| that a β₂-distance
    proxy ignores).  Samples are strictly interior (excluding the endpoints,
    which are crossing root phases).

    Returns ``(theta2, interval_index)`` — the interval index lets the caller
    place the seed in the region's cyclic order without a fragile
    circular-containment test.
    """
    twopi = 2 * math.pi
    fracs = np.arange(1, n_per_interval + 1) / (n_per_interval + 1)

    best_t2 = float(intervals[0][0])
    best_f = -1.0
    best_i = 0
    for i, (a, b) in enumerate(intervals):
        if b < a:
            b = b + twopi
        for fr in fracs:
            t2 = (a + fr * (b - a)) % twopi
            fmin = _loop_min_f(t2, zm, poly)
            if fmin > best_f:
                best_f = fmin
                best_t2 = t2
                best_i = i
    return best_t2, best_i


def compute_average_winding(
    zm: ZeroManager,
    poly: CharPoly,
    charges: list[dict],
) -> float:
    """Compute the average major-axis winding number ``W(E_ref, mu1)``.

    Topology (§3): crossings come in two kinds:
      * ordinary (charge ±1, SOFT): the winding across it is fixed by its
        charge, so it only partitions θ₂ into intervals *within* a region.
      * mr / tangent / unknown (charge 0, HARD): charge unknown, so it
        DELIMITES regions.  Across a hard boundary the winding is recomputed
        independently.

    Therefore: partition the circle into REGIONS delimited by hard
    boundaries; within each region, ordinary boundaries split it into
    intervals.  Pick ONE seed interval per region (the safest — farthest
    from all boundaries), compute w₀ there via :func:`_loop_winding_quad`,
    and propagate across the region's soft boundaries via charges.

    No boundary grouping: coincident boundaries (an MR's θ₂_a == θ₂_b, or
    two ordinary crossings sharing a θ₂) create zero-width intervals, which
    contribute 0 to the arc-weighted mean; sequential charge propagation
    handles a +1/−1 pair at one θ₂ correctly.
    """
    m = _ensure_mu2mid(zm, poly)
    E_ref = m.E_ref
    mu1 = m.mu1
    twopi = 2 * math.pi

    # Boundary list (θ₂, is_hard, dc): the winding change as θ₂ increases
    # PAST this boundary.  Each zero of f is a SINGLE boundary at its θ₂ with
    # its own charge (§3.1): dc = charge = sign(g') = the zero-curve's
    # direction through μ₂_mid, which equals the loop-winding jump at that θ₂
    # (verified by direct evaluation on both sides).  No pairing — the two
    # boundary tracks at the same θ₁ are two independent zeros (two θ₂).
    boundaries: list[tuple[float, bool, int]] = []
    for ch in charges:
        hard = ch['kind'] in ('mr', 'tangent', 'unknown')
        boundaries.append((ch['theta2'] % twopi, hard, ch['charge']))

    # ≤ 1 boundary: the full θ₂ circle is one region with constant winding.
    if len(boundaries) <= 1:
        t2, _ = _pick_seed_theta2([(0.0, twopi)], m, poly)
        w0 = _loop_winding_quad(m, poly, E_ref, mu1, t2)
        return float(round(w0))

    boundaries.sort(key=lambda x: x[0])
    N = len(boundaries)
    bdry_t2 = [b[0] for b in boundaries]
    bdry_hard = [b[1] for b in boundaries]
    bdry_dc = [b[2] for b in boundaries]

    def _interval(k: int) -> tuple[float, float]:
        a = bdry_t2[k]
        b = bdry_t2[(k + 1) % N]
        if b < a:
            b += twopi
        return a, b

    def _width(k: int) -> float:
        a, b = _interval(k)
        return b - a

    # regions: split at hard boundaries
    hard_idx = [k for k in range(N) if bdry_hard[k]]
    if not hard_idx:
        regions: list[list[int]] = [list(range(N))]
    else:
        regions = []
        for r in range(len(hard_idx)):
            start = hard_idx[r]
            end = hard_idx[(r + 1) % len(hard_idx)]
            if end > start:
                regions.append(list(range(start, end)))
            else:
                regions.append(list(range(start, N)) + list(range(0, end)))

    windings: list[int] = [0] * N
    for region in regions:
        valid = [(p, idx) for p, idx in enumerate(region) if _width(idx) > 0]
        if not valid:
            continue
        valid_intervals = [_interval(idx) for _, idx in valid]
        t2_seed, best_valid = _pick_seed_theta2(
            valid_intervals, m, poly,
        )
        best_p = valid[best_valid][0]
        seed_idx = region[best_p]
        w0 = int(round(_loop_winding_quad(m, poly, E_ref, mu1, t2_seed)))
        windings[seed_idx] = w0

        # forward propagation across soft boundaries
        w = w0
        for p in range(best_p + 1, len(region)):
            idx = region[p]
            w = w + bdry_dc[idx]
            windings[idx] = w
        # backward propagation
        w = w0
        for p in range(best_p - 1, -1, -1):
            idx = region[p]
            w = w - bdry_dc[region[p + 1]]
            windings[idx] = w

    total = 0.0
    for k in range(N):
        total += windings[k] * _width(k)
    return total / twopi


# ---------------------------------------------------------------------------
# Convenience: detection + winding in one call
# ---------------------------------------------------------------------------

def detect_crossings_and_winding(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    crossing_tol: float = 1e-10,
    max_newton: int = _MAX_NEWTON_ITER,
    zm_run_kwargs: dict | None = None,
) -> tuple[list[PointSubset], float]:
    """Crossing detection + average major-axis winding.

    Builds μ₂_mid once (in :func:`detect_crossings_simple` via
    :func:`_ensure_mu2mid`); the winding reuses that built ``Mu2MidZM`` — the
    shared instance is read from the crossing detector's bootstrap.  When the
    caller passes a plain ``ZeroManager`` both stages rebuild a fresh
    ``Mu2MidZM`` (the bisection path avoids this by constructing one itself).

    Returns ``(subsets, W_avg)``.
    """
    m = _ensure_mu2mid(zm, poly, **(zm_run_kwargs or {}))
    subsets, charges = detect_crossings_simple(
        m, poly,
        crossing_tol=crossing_tol,
        max_newton=max_newton,
    )
    W_avg = compute_average_winding(m, poly, charges)
    return subsets, W_avg
