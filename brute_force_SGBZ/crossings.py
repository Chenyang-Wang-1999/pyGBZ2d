'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-13
Copyright © Department of Physics, Tsinghua University. All rights reserved

§2 交点探测与求解 — per-column vs μ₂_mid (unified with amoeba).

Refactored 2026-08-13 (see ``log/2026-08-13-SGBZ算法梳理.md`` §2).  The
pairwise column-pair sweep and Newton refinement are gone.  Now:

  * μ₂_mid is a first-class piecewise-smooth object
    (``brute_force_SGBZ.mu2mid.Mu2MidZM``).  ``build_mu2_mid`` refines every
    sort-change crossing into the mesh, so the post-build mesh is clean.
  * Crossing detection = each track column vs μ₂_mid, mirroring amoeba's
    ``extract_amoeba_subsets`` structure — the only difference is μ₂ is the
    curve ``μ₂_mid(θ₁)`` instead of a constant (§2.1).

Why per-column has no false positives (§2.2): μ₂_mid sits in the gap between
the boundary-pair moduli, where no other track lives.  A track crossing
μ₂_mid must become a boundary column at the crossing instant ⟺ the boundary
gap closes ⟺ an SGBZ point.  So "track j crosses μ₂_mid" is *equivalent* to
the SGBZ point set — including multi-index jumps (a track leaping M+3→M-2
still crosses μ₂_mid once and is caught), which the old sort-change bracket
missed.

The topological charge is read inline from the crossing direction
(§3.1): ``charge = sign(g')`` where ``g = ln|β_j| − μ₂_mid``; at the
crossing ``j`` is a boundary column so ``sign(g') = sign(½ gap')``.

Mesh-mutation contract (2026-08-14): the mesh is mutated ONLY by
``Mu2MidZM.build_mu2_mid``.  The crossing refinement probes θ by solving
roots there transiently (:func:`_eval_g` via :meth:`ZeroManager.solve_at`)
and never inserts a row — so every row index the detection sweep records
stays valid for the whole sweep, and the touch/charge reads can never pair
a β₂ with a stale θ₁ (an earlier version inserted one row per bracket
probe, which shifted indices mid-sweep).

Assumes NO continuum is present — the caller must gate with
``brute_force_SGBZ.continuum_lines.detect_continuum_simple`` (or
``Mu2MidZM.has_continuum``) first.
'''

from __future__ import annotations

import math

import numpy as np
from cmath import exp

from gbz_types import CharPoly, PointSubset, circ_dist
from continuation import ZeroManager

from .mu2mid import (
    Mu2MidZM, _Mu2MidPath, _cubic_hermite_coeffs, _cubic_roots_in_interval,
)


# ---------------------------------------------------------------------------
# Tunables
# ---------------------------------------------------------------------------

# Convergence tolerance for the cubic-Hermite bracketing of a crossing.
_CROSSING_TOL: float = 1e-10

# Max cubic-Hermite bracketing iterations per crossing.
_MAX_BRACKET_ITER: int = 100
# Backward-compat alias (sgbz_solver / plateau import the old Newton name).
_MAX_NEWTON_ITER: int = _MAX_BRACKET_ITER

# θ₁ within this (circular) distance of a boundary MR → the crossing is the
# MR's own echo and is dropped in favour of the exact MR record.
_MR_PROXIMITY_TOL: float = 1e-4


# ---------------------------------------------------------------------------
# Mu2MidZM bootstrapping
# ---------------------------------------------------------------------------

def _ensure_mu2mid(zm: ZeroManager, poly: CharPoly, **run_kwargs) -> Mu2MidZM:
    """Return a Mu2MidZM with μ₂_mid built, reusing *zm* if possible.

    A plain ``ZeroManager`` cannot be mutated in place (it lacks
    ``build_mu2_mid``), so a fresh ``Mu2MidZM`` is constructed from *zm*'s
    ``(poly, E_ref, mu1)`` and built.  When *zm* is already a built
    ``Mu2MidZM`` it is returned as-is (the bisection path builds once and
    reuses, avoiding a rebuild).
    """
    if isinstance(zm, Mu2MidZM) and getattr(zm, '_mu2_mid_built', False):
        return zm
    m = Mu2MidZM(zm.poly, zm.E_ref, zm.mu1)
    m.run(**run_kwargs)
    m.build_mu2_mid()
    return m


# ---------------------------------------------------------------------------
# Per-column cubic-Hermite bracketing on g = ln|β_j| − μ₂_mid
# ---------------------------------------------------------------------------

def _track_g_and_gp(zm: Mu2MidZM, s_idx: int, row: int, j: int) -> tuple[float, float]:
    """``g = ln|β_j| − μ₂_mid`` and ``g' = Re(V_j) − dμ₂_mid/dθ₁`` at a row,
    using the BUILT μ₂_mid curve (``seg_mu2_values`` / ``seg_mu2_derivs``).

    This is the track-vs-μ₂_mid comparison the design mandates (§2): the
    zero-curve ``ln|β_j(θ₁)|`` against the piecewise-smooth μ₂_mid curve.
    It is numerically *not* the same as a track-vs-track (gap) comparison —
    μ₂_mid is an average of two boundary tracks, and using the pre-built
    curve keeps it consistent across the seam (where the boundary pair swaps)
    rather than re-deriving a per-interval mean that drifts between the two
    seam copies.  NaN derivative when tangents are unavailable.
    """
    seg = zm.segments[s_idx]
    # RAW track ln|β₂| (see mu2mid._LOGABS_CLAMP_L): only μ₂_mid itself is
    # clamped.  A 0/∞ padding track gives g = ±∞ — it never changes sign,
    # so it contributes no crossing (and M = 0 / N = 0 models are
    # unsolvable by design — see doc/SGBZ.md).
    la = float(np.log(np.abs(seg.tracked_roots[row, j])))
    g = la - float(zm.seg_mu2_values[s_idx][row])
    if seg.tangents is None:
        return g, float('nan')
    gp = float(seg.tangents[row, j].real) - float(zm.seg_mu2_derivs[s_idx][row])
    return g, gp


def _charge_at_row(
    zm: Mu2MidZM, s_idx: int, row: int, j: int, theta1: float,
) -> tuple[complex, float, int, str]:
    """Atomically read ``(beta2, g', charge, kind)`` at a crossing row.

    The crossing's β₂, its track-vs-μ₂_mid derivative ``g'``, and the
    resulting topological charge are all taken from the **same mesh row at
    the same instant**.  The mesh is static during detection (only
    ``build_mu2_mid`` mutates — 2026-08-14 contract), so no later bracket
    can shift row indices; the charge is still classified *here* rather
    than by a post-hoc ``m.locate(theta1)`` in the materialisation loop so
    that ``g'`` stays bound to the β₂ of this row.
    """
    seg = zm.segments[s_idx]
    beta2 = complex(seg.tracked_roots[row, j])
    _, gp = _track_g_and_gp(zm, s_idx, row, j)
    near_mr = _is_near_mr(theta1, zm)
    ch = _classify_charge(theta1, beta2, gp, near_mr)
    return beta2, gp, ch['charge'], ch['kind']


def _eval_g(
    zm: Mu2MidZM,
    path: _Mu2MidPath,
    s_idx: int,
    i: int,
    j: int,
    theta: float,
    cache: dict,
) -> tuple[float, float, complex]:
    """``(g, g', β₂_j)`` at θ inside mesh interval ``[th[i], th[i+1]]`` — no mesh mutation.

    The zero-mutation evaluation the crossing refinement uses instead of the
    old ``insert_solution`` probe: an exact mesh row is read from the stored
    track/μ₂_mid arrays (:func:`_track_g_and_gp`); any other θ is SOLVED
    transiently via :meth:`ZeroManager.solve_at` (solve + Hungarian onto the
    track frame + tangent) and the row discarded.  Results are cached per θ
    — a predicted point becomes one of the next iteration's bracket
    endpoints, whose true (g, g') must be reused, not re-derived.  The cache
    is keyed by the θ float as assigned (``theta_pred`` is passed around
    unchanged, so dict equality is exact); a ulp mismatch would at worst
    cause one redundant deterministic solve, never a wrong value.

    The mesh-row hit test is EXACT float equality (not a tolerance): inside
    the 2π seam's ``_BOUNDARY_THETA_TOL`` wrap band ``solve_at`` raises
    ``ValueError`` (the wrapped θ leaves the interval), which is what
    reproduces the old near-seam abandonment — a tolerance-based hit would
    instead treat such a θ_pred as a mesh row and emit a duplicate of the
    interval-0 crossing.
    """
    if theta in cache:
        return cache[theta]
    seg = zm.segments[s_idx]
    th = seg.theta1_arr
    row = None
    if theta == th[i]:
        row = i
    elif theta == th[i + 1]:
        row = i + 1
    if row is not None:
        g, gp = _track_g_and_gp(zm, s_idx, row, j)
        result = (g, gp, complex(seg.tracked_roots[row, j]))
    else:
        roots, V = zm.solve_at(theta, seg_idx=s_idx, i=i)
        la = float(np.log(np.abs(roots[j])))  # RAW track ln|β₂|
        mu2, dmu2 = path.value_deriv(theta)
        g = la - mu2
        gp = float(V[j].real) - dmu2
        result = (g, gp, complex(roots[j]))
    cache[theta] = result
    return result


def _bracket_crossing(
    zm: Mu2MidZM,
    j: int,
    s_idx: int,
    theta_lo: float,
    theta_hi: float,
    *,
    xtol: float,
    max_iter: int,
    path: _Mu2MidPath | None = None,
) -> tuple[float, complex, float, int, str, bool] | None:
    """Refine ``g(θ) = ln|β_j| − μ₂_mid(θ) = 0`` in ``[theta_lo, theta_hi]``.

    ZERO MESH MUTATION (2026-08-14): the bracket probes θ by solving roots
    there transiently (:func:`_eval_g`) and discards the row — only
    ``build_mu2_mid`` ever mutates the mesh.  The bracket therefore lives
    entirely inside ONE original adjacent mesh interval ``[th[i], th[i+1]]``
    (tightening only shrinks it), located once at entry; the old version
    inserted each probe as a row and had to re-resolve both endpoints
    against the mutated mesh every iteration.

    μ₂_mid is the BUILT piecewise-smooth curve, evaluated through
    :class:`_Mu2MidPath` (values + analytic derivative, correct one-sided
    derivatives at breakpoints).  So::

        g      = ln|β_j| − μ₂_mid(θ)
        g'     = Re(V_j) − dμ₂_mid/dθ₁     (track vs curve)

    Cubic-Hermite bracketing (§2.3): build the cubic from the endpoints'
    (value, g') → ``np.roots`` predicts θ_pred → ``_eval_g`` takes the true
    g_pred and true direction → derivative sign + g_pred sign fixes which
    side of the root θ_pred is on → tighten bracket.  Converges on bracket
    width < xtol.  Robust where Newton diverges (branch points) and where a
    sign-change check deadlocks (g_pred agreeing with both ends).

    Returns ``(theta1, beta2, g_prime, charge, kind, converged)``; ``None``
    if the interval has no transversal root (cubic false positive) or the
    direction is undefined (tangency / branch point).  On width-convergence
    the root is reported at the LAST predicted point — its β₂ is the true
    solved root there, more accurate than the mesh-row fallback the old
    version returned.  The charge is classified here from the SAME
    (θ, β₂, g') triple the bracket converged on, so g' stays bound to the
    β₂ it describes.
    """
    if path is None:
        path = _Mu2MidPath(zm)

    # locate the original containing interval ONCE — the mesh is static
    # after build_mu2_mid and every θ_pred stays inside [th[i], th[i+1]].
    si, i = zm.locate(theta_lo)
    # the bracket stays inside the originating segment (μ₂_mid is smooth
    # there; crossing a segment boundary = an MR breakpoint, handled by
    # _mr_boundary_entries, not here).
    if si != s_idx:
        return None
    seg = zm.segments[si]
    th = seg.theta1_arr
    if i + 1 >= len(th) or abs(th[i] - theta_lo) > 1e-15 \
       or abs(th[i + 1] - theta_hi) > 1e-15:
        return None

    cache: dict[float, tuple[float, float, complex]] = {}
    last_pred: tuple[float, complex, float] | None = None

    def _finish(theta: float, b2: complex, gp: float, converged: bool):
        ch = _classify_charge(theta, b2, gp, _is_near_mr(theta, zm))
        return theta, b2, gp, ch['charge'], ch['kind'], converged

    for _ in range(max_iter):
        v0, dv0, b2_0 = _eval_g(zm, path, si, i, j, theta_lo, cache)
        v1, dv1, b2_1 = _eval_g(zm, path, si, i, j, theta_hi, cache)

        # endpoint touches: the root IS this θ (a mesh row or a previously
        # predicted point whose true g vanished) — no further iteration.
        if abs(v0) < xtol:
            return _finish(theta_lo, b2_0, dv0, True)
        if abs(v1) < xtol:
            return _finish(theta_hi, b2_1, dv1, True)

        # bracket converged by width → root at the best-known point: the
        # last prediction (true solved β₂), else the left row.
        if theta_hi - theta_lo < xtol:
            if last_pred is not None:
                theta, b2, gp = last_pred
            else:
                theta, b2, gp = theta_lo, b2_0, dv0
            return _finish(theta, b2, gp, True)

        # a transversal crossing must have opposite signs; same sign ⇒ no
        # root (abandon — the caller only passes sign-change intervals, but
        # cubic refinement can land on a same-sign sub-range).
        if v0 * v1 > 0:
            return None

        h = theta_hi - theta_lo
        # divergent / unavailable tangent → linear fallback (bounded, §2.3).
        if not (np.isfinite(dv0) and np.isfinite(dv1)):
            dv0 = (v1 - v0) / h
            dv1 = dv0

        coeffs = _cubic_hermite_coeffs(h, v0, dv0, v1, dv1)
        s_roots = _cubic_roots_in_interval(coeffs, h)
        if not s_roots:
            return None

        # pick the transversal (cubic g' ≠ 0) root nearest the midpoint;
        # a cubic g'≈0 root is a tangent touch, not a transversal crossing.
        s_pred = None
        for s in s_roots:
            dg_cubic = (3.0 * coeffs[0] * s * s
                        + 2.0 * coeffs[1] * s + coeffs[2])
            if abs(dg_cubic) < 1e-15:
                continue
            if s_pred is None or abs(s - h / 2) < abs(s_pred - h / 2):
                s_pred = s
        if s_pred is None:
            return None

        theta_pred = theta_lo + s_pred
        if abs(theta_pred - theta_lo) < 1e-15 \
           or abs(theta_pred - theta_hi) < 1e-15:
            return None

        try:
            # true g_pred and true direction sign from a TRANSIENT solve at
            # θ_pred — the row is discarded, the mesh is never mutated.
            # ValueError reproduces the old insert_solution guard (e.g. a
            # θ_pred inside the 2π seam's wrap band) → abandon the crossing.
            g_pred, gp_pred, b2_pred = _eval_g(
                zm, path, si, i, j, theta_pred, cache)
        except (ValueError, RuntimeError):
            return None

        # compute_tangent returns V_j = inf (not NaN) at a multiple root,
        # so the finiteness check must stay alongside the ~0 check.
        if not np.isfinite(gp_pred) or abs(gp_pred) < 1e-15:
            return None  # direction undefined — abandon

        last_pred = (theta_pred, b2_pred, gp_pred)

        # tighten: direction sign + g_pred sign fixes the side.
        increasing = gp_pred > 0
        if (g_pred > 0) == increasing:
            theta_hi = theta_pred
        else:
            theta_lo = theta_pred

    return None



# ---------------------------------------------------------------------------
# Dedup & charge helpers (carried over from the previous version)
# ---------------------------------------------------------------------------


def _mr_boundary_entries(zm: Mu2MidZM, M: int) -> list[dict]:
    """Multiple roots sitting exactly on the PMGBZ boundary.

    An MR whose cluster covers sorted positions M-1 AND M is a subset point:
    those positions hold the snapped degenerate root, so ``|β_M| == |β_{M+1}|``
    holds exactly there.  ``MultipleRootInfo.roots`` is the modulus-sorted
    (padded, snapped) root array and ``cluster_indices`` indexes into it, so
    the boundary test is a direct position lookup — no mesh scan, no
    refinement.  This matters because the mesh scan has a systematic blind
    spot at such rows: the snapped cluster makes μ₂_mid exactly equal to both
    boundary columns' logabs in floating point.
    """
    entries: list[dict] = []
    for mr in zm.multiple_roots:
        for cluster in mr.cluster_indices:
            if (M - 1) in cluster and M in cluster:
                beta2 = complex(mr.roots[cluster[0]])
                entries.append({
                    'theta1': float(mr.theta1),
                    'beta2': beta2,
                    'theta2': float(np.angle(beta2)) % (2 * math.pi),
                })
                break  # a sorted position belongs to at most one cluster
    return entries


def _is_near_mr(theta1: float, zm: Mu2MidZM) -> bool:
    """Whether *theta1* is within _MR_PROXIMITY_TOL of any known MR."""
    for mr in zm.multiple_roots:
        if circ_dist(theta1, mr.theta1) < _MR_PROXIMITY_TOL:
            return True
    return False


def _classify_charge(
    theta1: float,
    beta2: complex,
    g_prime: float,
    near_mr: bool,
) -> dict:
    """Topological charge of ONE zero of f at ``(θ₁, θ₂)`` (§3.1).

    Each zero-curve (track j) crossing μ₂_mid is a single zero of the char
    poly at ``(θ₁*, θ₂_j(θ₁*))`` — NOT a pair.  Its charge is the direction
    of the zero-curve through μ₂_mid, ``charge = sign(g')`` where
    ``g' = d(ln|β_j| − μ₂_mid)/dθ₁`` (the track-vs-curve derivative).  That
    equals the loop-winding jump as θ₂ crosses ``θ₂_j(θ₁*)`` (verified by
    direct evaluation on both sides).

    A sign-change crossing is transversal by construction, so its charge is
    ``sign(g')`` regardless of ``|g'|``.  Near-MR branch points (``g'``
    divergent) carry charge 0 as a hard region boundary.

    Returns a dict with ``theta1, theta2, charge, kind``.  *kind* is
    'ordinary' (charge from sign(g')), 'mr' (near a multiple root), or
    'tangent' (g' not finite).  The latter two carry charge 0 and act as hard
    region boundaries — they delimit regions but don't propagate winding.
    """
    theta2 = float(np.angle(beta2)) % (2 * math.pi)
    base = dict(theta1=theta1, theta2=theta2)

    if near_mr:
        return {**base, 'charge': 0, 'kind': 'mr'}
    if not math.isfinite(g_prime):
        return {**base, 'charge': 0, 'kind': 'tangent'}
    return {**base, 'charge': 1 if g_prime > 0 else -1, 'kind': 'ordinary'}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def detect_crossings_simple(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    crossing_tol: float = _CROSSING_TOL,
    max_newton: int = _MAX_BRACKET_ITER,
    zm_run_kwargs: dict | None = None,
) -> tuple[list[PointSubset], list[dict]]:
    """0D PMGBZ-boundary crossing detection + charge classification.

    Per-column vs μ₂_mid (§2): builds μ₂_mid (the piecewise-smooth
    boundary-pair mean, refined at every sort-change), then — AFTER the
    build — detects each zero-curve (track) crossing the μ₂_mid curve.
    Assumes no continuum — gate with
    :func:`brute_force_SGBZ.continuum_lines.detect_continuum_simple` first.

    Algorithm
    ---------
    1. **Build μ₂_mid** (:meth:`Mu2MidZM.build_mu2_mid`): walls + sort-change
       refinement via cubic-Hermite bracketing.  The post-build mesh is clean
       and the μ₂_mid curve is available as the first-class arrays
       ``seg_mu2_values`` / ``seg_mu2_derivs`` (values + analytic derivative).
    2. **Per-column detection (post-build)**: for each track j compute
       ``g_j = ln|β_j| − μ₂_mid_values`` against the BUILT curve (track-vs-curve,
       numerically distinct from a track-vs-track gap comparison, §2.1).  A
       genuine **sign change** of g_j → bracket & refine (cubic Hermite).  No
       near-zero / near-miss shortcut — that used to fire spuriously where g_j
       hovers at machine-ε (the degenerate seam) without an actual crossing.
    3. **No dedup**: each detected zero-curve crossing (track j at θ₁*) is
       kept as-is — distinct tracks at the same θ₁ are distinct zeros (the
       two boundary tracks at a PMGBZ point are two independent zeros, two
       θ₂), and any seam/echo double-counting is a detection issue to fix
       at the source, not papered over with a merge.
    4. **MR echo drop**: crossings within _MR_PROXIMITY_TOL of a boundary MR
       are replaced by the exact ZeroManager MR record.
    5. **Charge + PointSubset construction**: charge = sign(g') where
       ``g' = Re(V_j) − μ₂_mid_derivs`` (the track-vs-curve direction, §3.1);
       one PointSubset per detected zero-curve (a PMGBZ point crossed by
       both boundary tracks yields two) + one per MR.

    .. note::
       The detection phase does NOT mutate the mesh (2026-08-14 contract):
       only :meth:`Mu2MidZM.build_mu2_mid` — run inside :func:`_ensure_mu2mid`
       — inserts rows; the crossing brackets solve transiently and discard.
       Consequently:

       - run continuum detection BEFORE this function;
       - the scan is re-callable on the same built *zm*.

    Returns
    -------
    subsets : list[PointSubset]
        One PointSubset per detected zero-curve crossing (track j crossing
        μ₂_mid at θ₁*) plus one per boundary MR, in detection order.  The
        two boundary tracks of a PMGBZ point are two independent zeros, so
        one PMGBZ point typically yields two PointSubsets.
    charges : list[dict]
        One charge dict per surviving crossing (``kind='ordinary'``,
        charge ±1) plus one per boundary MR (``kind='mr'``, charge 0).
    """
    M = poly.M
    K = poly.M + poly.N
    if M >= K:
        raise ValueError(f"M={M} >= K={K}: no PMGBZ boundary to check")
    if M <= 0:
        raise ValueError(f"M={M} <= 0: invalid boundary index")

    m = _ensure_mu2mid(zm, poly, **(zm_run_kwargs or {}))

    E_ref = m.E_ref
    mu1 = m.mu1
    max_iter = max_newton
    # One path evaluator for the whole sweep: the mesh is static during
    # detection, so the μ₂_mid curve does not change between brackets.
    path = _Mu2MidPath(m)

    # ---- collect per-column crossing candidates ----
    # Each candidate is ONE zero-curve (track j) crossing the BUILT μ₂_mid
    # curve at (θ₁*, θ₂_j(θ₁*)); the charge is the crossing direction
    # sign(g') (§3.1).  The two boundary tracks at the same θ₁* are TWO
    # independent zeros (two θ₂) — they are not paired.  Only a genuine sign
    # change of g_j counts (no near-zero shortcut, which fired spuriously at
    # the degenerate seam).
    raw: list[tuple[float, int, float]] = []

    for s_idx, seg in enumerate(m.segments):
        th = seg.theta1_arr
        N = len(th)
        if N < 2:
            continue
        la = np.log(np.abs(seg.tracked_roots))           # (N, K) RAW track ln|β₂|
        mu2 = m.seg_mu2_values[s_idx]                    # (N,) built μ₂_mid curve
        for j in range(K):
            g = la[:, j] - mu2                            # (N,) track j vs μ₂_mid
            # Sign changes (transversal crossings) and exact touches (g == 0
            # at a mesh row — a tangency or a crossing landing on a sample)
            # are disjoint: an exact 0 makes g[i]·g[i+1] == 0, never < 0.
            #
            # Touch detection runs on g[:-1] (left-closed / right-open): the
            # right endpoint row is shared with the next segment (an MR) or is
            # the θ₁=2π seam ≡ θ₁=0 of the first segment, so counting it here
            # would double the point.  An exact touch IS the root — no bracket
            # iteration needed, just read it off the row.
            sc = np.flatnonzero(g[:-1] * g[1:] < 0)       # sign-change intervals
            touch = np.flatnonzero(g[:-1] == 0)           # exact-touch rows, excl. last
            for i in sc:
                res = _bracket_crossing(
                    m, j, s_idx, float(th[i]), float(th[i + 1]),
                    xtol=crossing_tol, max_iter=max_iter, path=path,
                )
                if res is None:
                    continue
                t1, beta2, _gp, charge, kind, _ = res
                raw.append((t1, beta2, charge, kind))
            for i in touch:
                beta2, _gp, charge, kind = _charge_at_row(
                    m, s_idx, i, j, float(th[i]))
                raw.append((float(th[i]), beta2, charge, kind))

    # No dedup: each detected zero-curve crossing (track j at θ₁*) is kept
    # as-is — distinct tracks at the same θ₁ are distinct zeros, and any
    # seam/echo double-counting is a detection issue to fix at the source,
    # not papered over with a merge.
    mr_entries = _mr_boundary_entries(m, M)
    if mr_entries:
        raw = [
            r for r in raw
            if not any(circ_dist(r[0], e['theta1']) < _MR_PROXIMITY_TOL
                       for e in mr_entries)
        ]

    # ---- one zero → one PointSubset + one single-θ₂ charge ----
    # charge (β₂, sign(g')) was fixed inside _bracket_crossing at the root
    # row; here we only materialise it.  Re-locating would still pair the
    # charge with the wrong β₂ — the crossing θ is a solved point, not a
    # mesh row.
    subsets: list[PointSubset] = []
    charges: list[dict] = []
    for t1, beta2, charge, kind in raw:
        theta2 = float(np.angle(beta2)) % (2 * math.pi)
        subsets.append(PointSubset(E=E_ref, beta1=exp(mu1 + 1j * t1), beta2=beta2))
        charges.append(dict(theta1=t1, theta2=theta2,
                            charge=charge, kind=kind))

    for e in mr_entries:
        subsets.append(PointSubset(
            E=E_ref, beta1=exp(mu1 + 1j * e['theta1']), beta2=e['beta2']))
        charges.append(dict(theta1=e['theta1'], theta2=e['theta2'],
                            charge=0, kind='mr'))
    return subsets, charges
