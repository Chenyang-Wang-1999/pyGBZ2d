'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-10
Copyright © Department of Physics, Tsinghua University. All rights reserved

Topological charge + average major-axis winding, plus the 0D
PointSubset materialization from analyzed EventGroups.

Computes the average major-axis winding number ``W(E_ref, mu1)`` — the
quantity whose zero in ``mu1`` defines the SGBZ — from a
``Mu2MidZM`` (a ``ZeroManager`` with the piecewise-smooth μ₂_mid built) at
fixed ``(E_ref, mu1)`` together with the charge list from
:func:`detect_crossings_simple` below.

Design rationale:

  * The loop is ``β₂ = exp(μ₂_mid(θ₁) + iθ₂)`` with ``θ₁ ∈ [0, 2π)`` and
    fixed ``θ₂``; ``μ₂_mid`` is the piecewise-smooth boundary-pair mean — a
    first-class object with analytic values and derivatives and left/right
    derivatives at its breakpoints. The loop winding is the integral
    of ``Im[f'/f]`` over ``θ₁``, split at the μ₂_mid breakpoints so each quad
    segment is smooth. A global interpolant across a derivative jump would
    distort the loop. The path instead uses local Hermite interpolation
    (linear fallback at singular tangents), and the integrand uses the
    derivative of that interpolated path.
  * Crossings come in two kinds. **Ordinary** (charge -1/0/+1, SOFT): the winding
    across it is fixed by its charge, so it only partitions ``θ₂`` into
    intervals *within* a region.  **MR / tangent / unknown** (HARD): the
    charge is UNKNOWN (the charge dict stores ``None``, never a numeric
    value), so it DELIMITES regions; across a hard boundary the winding is
    recomputed independently, and any accidental arithmetic on its charge
    fails with ``TypeError``.
  * Therefore: partition the circle into regions delimited by hard
    boundaries; within each region, ordinary boundaries split it into
    intervals. Pick one seed per region with positive angular width by
    maximizing the sampled minimum of ``|f|``, compute ``w₀`` there, and
    propagate across soft boundaries via charges. This reduces the number
    of quadratures near zeros, where ``f'/f`` is ill-conditioned.
'''

from __future__ import annotations
from typing import Optional

import math

import numpy as np
from cmath import exp
from scipy import integrate

from ..core import CharPoly, PointSubset, TWO_PI
from ..continuation import ZeroManager

from ..core import live_defaults

# Loop-winding integral settings: every caller rounds the result to an
# integer, so these stay loose on purpose.
WINDING_QUAD_EPSABS: float = 1e-3
WINDING_QUAD_EPSREL: float = 1e-3
WINDING_QUAD_LIMIT: int = 200
#: Samples per interval when picking the θ₂ loop-winding seed.
SEED_N_PER_INTERVAL: int = 4
from .mu2mid import Mu2MidZM, ensure_mu2mid


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
    keeps each quad segment away from derivative jumps. Otherwise the
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
            epsabs=WINDING_QUAD_EPSABS, epsrel=WINDING_QUAD_EPSREL,
            limit=WINDING_QUAD_LIMIT,
        )[0]
    return total / (TWO_PI)


# ---------------------------------------------------------------------------
# Piecewise-smooth μ₂_mid loop path
# ---------------------------------------------------------------------------
# The independent Mu2Mid object lives in mu2mid.py and is built by
# Mu2MidZM.analyze(); winding evaluates it directly.


def _loop_winding_quad(
    zm: Mu2MidZM,
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    theta2: float,
) -> float:
    """Loop winding number via quad of ``Im[f'/f]`` over the μ₂_mid loop.

    The loop is ``β₂ = exp(μ₂_mid(θ₁) + iθ₂)``, ``θ₁ ∈ [0, 2π)``, with the
    analytic ``dβ₂/dθ₁ = β₂ · μ₂_mid'(θ₁)`` of the interpolated path.
    The integration is split at value/derivative discontinuities; ordinary
    C1 knots need no split. Consecutive root-tracking segments share their
    MR-boundary row exactly, so the cross-segment seam contributes nothing
    extra (the breakpoint split already isolates it).
    """
    if zm.mu2_mid is None:
        raise RuntimeError("compute_average_winding requires Mu2MidZM.analyze() "
                           "to have built the mu2_mid path")
    path = zm.mu2_mid
    twopi = TWO_PI
    bp_thetas = [float(t) for t in path.breakpoints
                 if 0.0 < t < twopi]

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
    """Minimum sampled ``|f|`` at fixed *theta2* on the μ₂_mid mesh.

    Uses ``seg_mu2_values`` at the stored θ₁ rows. A small denominator
    makes the winding integrand ``Im[f'/f]`` sensitive to numerical error.
    Distance to a β₂ root alone omits the polynomial's local scale: near
    a simple root, ``|f|`` also scales with ``|∂f/∂β₂|``. Maximizing this
    sampled minimum is a seed-selection heuristic, not a certified lower
    bound on ``|f|`` between mesh rows or a geometric distance to its zeros.
    """
    if zm.mu2_mid is None:
        raise RuntimeError("_loop_min_f requires Mu2MidZM.analyze()")
    worst = math.inf
    for s_idx, seg in enumerate(zm.segments):
        th = seg.theta1_arr
        if len(th) == 0:
            continue
        mu2 = zm.seg_mu2_values[s_idx]
        beta1 = np.exp(zm.mu1 + 1j * th)
        beta2 = np.exp(mu2 + 1j * theta2)
        for i in range(len(th)):
            v = abs(poly.eval_val((zm.E_ref, complex(beta1[i]), complex(beta2[i]))))
            if v < worst:
                worst = v
                if worst == 0.0:
                    return 0.0
    return float(worst) if math.isfinite(worst) else 0.0


@live_defaults(n_per_interval="sgbz.winding:SEED_N_PER_INTERVAL")
def _pick_seed_theta2(
    intervals: list[tuple[float, float]],
    zm: Mu2MidZM,
    poly: CharPoly,
    *,
    n_per_interval: Optional[int] = None,
) -> tuple[float, int]:
    """Pick the candidate θ₂ maximizing the mesh-sampled minimum of ``|f|``.

    The loop-winding seed must lie inside one of the region's intervals and
    stay clear of char-poly zeros for the quad to be reliable, so the safest
    θ₂ maximises the sampled ``min |f(E, β₁(θ₁), β₂_loop(θ₁))|``.
    Candidate θ₂ samples are strictly interior (excluding the endpoints,
    which are crossing root phases).

    Returns ``(theta2, interval_index)`` — the interval index lets the caller
    place the seed in the region's cyclic order without a fragile
    circular-containment test.
    """
    twopi = TWO_PI
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

    Crossings come in two kinds:
      * ordinary (charge -1/0/+1, SOFT): the winding across it is fixed by its
        charge, so it only partitions θ₂ into intervals *within* a region.
      * mr / tangent / unknown (HARD): the charge is UNKNOWN — the charge
        dict stores ``None`` — so it DELIMITES regions.  Across a hard
        boundary the winding is recomputed independently; any accidental
        arithmetic on its charge raises ``TypeError``.

    Charge conservation: when EVERY boundary is soft, the charges must sum
    to zero (one full θ₂ circle returns the winding to itself); a non-zero
    sum means the crossing detector missed or duplicated a zero and raises
    ``RuntimeError``.  Any hard boundary disables this check because its
    charge is unknown.

    Therefore: partition the circle into REGIONS delimited by hard
    boundaries; within each region, ordinary boundaries split it into
    intervals. Pick one seed per region with positive width using the
    sampled minimum of ``|f|``, compute w₀ via :func:`_loop_winding_quad`,
    and propagate across the region's soft boundaries via charges.

    No boundary grouping: coincident boundaries (an MR's θ₂_a == θ₂_b, or
    two ordinary crossings sharing a θ₂) create zero-width intervals, which
    contribute 0 to the arc-weighted mean; sequential charge propagation
    handles a +1/−1 pair at one θ₂ correctly.
    """
    m = ensure_mu2mid(zm)
    E_ref = m.E_ref
    mu1 = m.mu1
    twopi = TWO_PI

    # Boundary list (θ₂, is_hard, dc).  For an ORDINARY boundary
    # dc = charge = sign(g') is the loop-winding jump as θ₂ increases PAST
    # it (verified by direct evaluation on both sides).  For a HARD boundary
    # (mr/tangent/unknown) the physical charge is unknown: the dict stores
    # None — never a numeric placeholder — so dc is None and any accidental
    # arithmetic on it raises TypeError.  Hard boundaries delimit regions and
    # are never used for propagation.  No pairing — the two boundary tracks
    # at the same θ₁ are two independent zeros (two θ₂).
    boundaries: list[tuple[float, bool, int | None]] = []
    for ch in charges:
        hard = ch['kind'] in ('mr', 'tangent', 'unknown')
        boundaries.append((ch['theta2'] % twopi, hard, ch['charge']))

    # Charge conservation applies ONLY when every boundary is SOFT
    # (ordinary, charge ±1): one full θ₂ circle must return the winding to
    # itself, so the soft charges must sum to zero.  A HARD boundary has
    # UNKNOWN charge (stored as None), so its presence disables the check.
    if not any(b[1] for b in boundaries):
        charge_sum = sum(b[2] for b in boundaries)
        if charge_sum != 0:
            raise RuntimeError(
                f"crossing charges do not sum to zero (sum={charge_sum}, "
                f"n_boundaries={len(boundaries)}) at E_ref={E_ref}, "
                f"mu1={mu1}: crossing detection is incomplete"
            )

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
# 0D PointSubset materialization from analyzed EventGroups
# ---------------------------------------------------------------------------
#
# The crossing *detection* lives in pairwise.py (EventGroups, built by
# Mu2MidZM.analyze); MR boundary rows flow through the same EventGroup
# machinery (marked ``is_mr``), so this section ONLY turns the analyzed
# groups into PointSubsets + charge dicts.  There is deliberately no
# separate MR materialization channel and no MR-echo drop: both were
# removed because the MR channel only materialized clusters straddling
# M-1/M, silently dropping any other modulus coincidence at a boundary
# row (the seam missed-detection at E=1.212), and the echo-drop
# proximity rule discarded legitimate dense events.


@live_defaults(crossing_tol="sgbz.pairwise:CROSSING_TOL")
def detect_crossings_simple(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    crossing_tol: Optional[float] = None,
) -> tuple[list[PointSubset], list[dict]]:
    """Materialize 0D PMGBZ PointSubsets from analyzed EventGroups.

    An EventGroup whose induced item/column component covers BOTH sorted
    positions M-1 and M is an SGBZ boundary point; every REAL column of that
    component gets its own PointSubset (representative items are expanded
    through their full clusters).  MR-row events (``is_mr``) and tangent
    events materialize with ``charge=None`` -- a hard boundary for the
    winding propagation.  Assumes no boundary continuum (the solver gates
    with ``has_continuum``).

    *crossing_tol* is accepted for backward compatibility and has no
    effect: the refinement tolerance is fixed at ``Mu2MidZM.analyze`` time
    (by the caller that built the instance), not here.

    Returns
    -------
    subsets : list[PointSubset]
        One PointSubset per REAL column of every M-1/M boundary event
        component.
    charges : list[dict]
        One dict per subset: ``charge`` is the side-change charge
        (+1 / -1 / 0) for ordinary events and ``None`` for mr/tangent.
    """
    M = poly.M
    K = poly.M + poly.N
    if M >= K:
        raise ValueError(f"M={M} >= K={K}: no PMGBZ boundary to check")
    if M <= 0:
        raise ValueError(f"M={M} <= 0: invalid boundary index")

    m = ensure_mu2mid(zm)

    E_ref = m.E_ref
    mu1 = m.mu1

    subsets: list[PointSubset] = []
    charges: list[dict] = []

    for g in m._event_groups:
        if not g.point_columns:
            continue
        if g.row < 0:
            continue
        seg = m.segments[g.seg_idx]
        theta1 = float(g.theta) % (TWO_PI)

        for col in g.point_columns:
            beta2 = complex(seg.tracked_roots[g.row, col])
            theta2 = float(np.angle(beta2)) % (TWO_PI)
            q = g.column_q.get(col)
            try:
                kind = g.column_kind[col]
            except KeyError as exc:
                raise RuntimeError(
                    f"EventGroup at θ₁={theta1} has no kind for column "
                    f"{col}; was finalize_event_groups skipped?"
                ) from exc
            subsets.append(PointSubset(
                E=E_ref, beta1=exp(mu1 + 1j * theta1), beta2=beta2))
            charges.append(dict(
                theta1=theta1, theta2=theta2, charge=q, kind=kind))
    return subsets, charges


# ---------------------------------------------------------------------------
# Convenience: detection + winding in one call
# ---------------------------------------------------------------------------

@live_defaults(crossing_tol="sgbz.pairwise:CROSSING_TOL")
def detect_crossings_and_winding(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    crossing_tol: Optional[float] = None,
) -> tuple[list[PointSubset], float]:
    """Crossing detection + average major-axis winding.

    Calls :func:`ensure_mu2mid` once, then passes the same analyzed
    ``Mu2MidZM`` to both point materialization and winding. A plain
    ``ZeroManager`` requires one fresh build; an already analyzed
    ``Mu2MidZM`` is reused directly.

    Returns ``(subsets, W_avg)``.
    """
    m = ensure_mu2mid(zm)
    subsets, charges = detect_crossings_simple(
        m, poly, crossing_tol=crossing_tol,
    )
    W_avg = compute_average_winding(m, poly, charges)
    return subsets, W_avg
