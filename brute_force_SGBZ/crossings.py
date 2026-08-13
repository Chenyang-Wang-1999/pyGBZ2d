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
    Mu2MidZM, _cubic_hermite_coeffs, _cubic_roots_in_interval, _dv_column,
)


# ---------------------------------------------------------------------------
# Tunables
# ---------------------------------------------------------------------------

# Convergence tolerance for the cubic-Hermite bracketing of a crossing.
_CROSSING_TOL: float = 1e-10

# Generous threshold for suspicious-interval near-miss detection (tangency
# candidates).  Intervals where |ln|β₂_j| − μ₂_mid| < this at both endpoints
# are flagged; the cubic refinement filters false positives.
_DETECT_THRESHOLD: float = 1e-2

# Max cubic-Hermite bracketing iterations per crossing.
_MAX_BRACKET_ITER: int = 100
# Backward-compat alias (sgbz_solver / plateau import the old Newton name).
_MAX_NEWTON_ITER: int = _MAX_BRACKET_ITER

# L²-distance threshold in (β₁, β₂) space for duplicate-crossing dedup.
_DEDUP_TOL: float = 1e-6

# |f| at the bracket's best iterate below this (without width-converging)
# → tangent touch: f reaches ~0 but the direction sign is unreliable.
# Genuine cubic false positives (tracks passing at distance) stay well above.
_TANGENT_F_TOL: float = 1e-6

# θ₁ within this (circular) distance of a boundary MR → the crossing is the
# MR's own echo and is dropped in favour of the exact MR record.
_MR_PROXIMITY_TOL: float = 1e-4

# |g'| below this at a crossing → tangency (charge 0, hard region boundary).
_TANGENCY_THRESHOLD: float = 1e-3


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
    the same instant** — the caller (``_bracket_crossing``) holds the row
    right after inserting/touching the root, before any later bracket on
    another track can shift row indices.  This is why charge is classified
    *here* rather than by a post-hoc ``m.locate(theta1)`` in the
    materialisation loop: that loop would re-locate against a mesh mutated
    by subsequent crossings and bind ``g'`` to the wrong β₂.
    """
    seg = zm.segments[s_idx]
    beta2 = complex(seg.tracked_roots[row, j])
    _, gp = _track_g_and_gp(zm, s_idx, row, j)
    near_mr = _is_near_mr(theta1, zm)
    ch = _classify_charge(theta1, beta2, gp, near_mr)
    return beta2, gp, ch['charge'], ch['kind']


def _bracket_crossing(
    zm: Mu2MidZM,
    j: int,
    s_idx: int,
    theta_lo: float,
    theta_hi: float,
    *,
    tol: float,
    xtol: float,
    max_iter: int,
) -> tuple[float, complex, float, int, str, bool] | None:
    """Refine ``g(θ) = ln|β_j| − μ₂_mid(θ) = 0`` in ``[theta_lo, theta_hi]``.

    μ₂_mid is the BUILT piecewise-smooth curve (values ``seg_mu2_values``,
    derivatives ``seg_mu2_derivs``) — a first-class object, not re-derived per
    interval.  So::

        g      = ln|β_j| − μ₂_mid_values[row]
        g'     = Re(V_j) − μ₂_mid_derivs[row]     (track vs curve)

    Cubic-Hermite bracketing (§2.3): build the cubic from the endpoints'
    (value, g') → ``np.roots`` predicts θ_pred → ``insert_solution`` takes
    the true g_pred and true direction → derivative sign + g_pred sign fixes
    which side of the root θ_pred is on → tighten bracket.  Converges on
    bracket width < xtol.  Robust where Newton diverges (branch points) and
    where a sign-change check deadlocks (g_pred agreeing with both ends).

    Returns ``(theta1, beta2, g_prime, charge, kind, converged)``; ``None``
    if the interval has no transversal root (cubic false positive) or the
    direction is undefined (tangency / branch point).  *beta2*, *g_prime*
    and *charge* are read atomically from the (inserted or touched) root
    row **at the instant it is identified** — before any later bracket on
    another track can insert rows and shift indices.  The charge is
    classified here (not in the materialisation loop) precisely so that
    *g'* and the β₂ it describes stay bound to the same row; a post-hoc
    ``locate(theta1)`` would re-resolve against a mutated mesh and pair
    *g'* with the wrong β₂.
    """
    for _ in range(max_iter):
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

        v0, _ = _track_g_and_gp(zm, si, i, j)
        v1, _ = _track_g_and_gp(zm, si, i + 1, j)

        # endpoint touches: the root IS a mesh row (e.g. a crossing landing
        # on the θ₁ = 2π seam).  Read its true direction and return.
        if abs(v0) < xtol:
            b2, gp, q, k = _charge_at_row(zm, si, i, j, theta_lo)
            return theta_lo, b2, gp, q, k, True
        if abs(v1) < xtol:
            b2, gp, q, k = _charge_at_row(zm, si, i + 1, j, theta_hi)
            return theta_hi, b2, gp, q, k, True

        # bracket converged by width → root at the left row.
        if theta_hi - theta_lo < xtol:
            b2, gp, q, k = _charge_at_row(zm, si, i, j, theta_lo)
            return theta_lo, b2, gp, q, k, True

        # a transversal crossing must have opposite signs; same sign ⇒ no
        # root (abandon — the caller only passes sign-change intervals, but
        # cubic refinement can land on a same-sign sub-range).
        if v0 * v1 > 0:
            return None

        h = theta_hi - theta_lo
        N = len(th)
        touches_mr = (
            (i == 0 and seg.left_mr >= 0)
            or (i == N - 2 and seg.right_mr >= 0)
        )
        _, dv0 = _track_g_and_gp(zm, si, i, j)
        _, dv1 = _track_g_and_gp(zm, si, i + 1, j)
        # MR / divergent tangent → linear fallback (bounded, §2.3).
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
            insert_at, _ = zm.insert_solution(theta_pred, si, i)
        except (ValueError, RuntimeError):
            return None

        seg = zm.segments[si]
        la = np.log(np.abs(seg.tracked_roots))
        g_pred = float(la[insert_at, j]) - float(zm.seg_mu2_values[si][insert_at])
        if seg.tangents is not None:
            g_prime = float(seg.tangents[insert_at, j].real) \
                      - float(zm.seg_mu2_derivs[si][insert_at])
        else:
            g_prime = (3.0 * coeffs[0] * s_pred * s_pred
                       + 2.0 * coeffs[1] * s_pred + coeffs[2])

        if not np.isfinite(g_prime) or abs(g_prime) < 1e-15:
            return None  # direction undefined — abandon

        # tighten: direction sign + g_pred sign fixes the side.
        increasing = g_prime > 0
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
    detect_threshold: float = _DETECT_THRESHOLD,
    max_newton: int = _MAX_BRACKET_ITER,
    dedup_tol: float = _DEDUP_TOL,
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
    3. **Dedup by θ*** (circular distance): multiple tracks detect the same
       PMGBZ point (§2.2 — every boundary-crossing track hits μ₂_mid there).
    4. **MR echo drop**: crossings within _MR_PROXIMITY_TOL of a boundary MR
       are replaced by the exact ZeroManager MR record.
    5. **Charge + PointSubset construction**: charge = sign(g') where
       ``g' = Re(V_j) − μ₂_mid_derivs`` (the track-vs-curve direction, §3.1);
       two PointSubsets per crossing (the boundary β₂ pair) + one per MR.

    .. note::
       This function MUTATES *zm* via :meth:`Mu2MidZM.build_mu2_mid` and the
       bracket's :meth:`insert_solution`.  Consequently:

       - run continuum detection BEFORE this function;
       - do not call this function twice on the same *zm*.

    Returns
    -------
    subsets : list[PointSubset]
        Two PointSubset per crossing (β₂_a and β₂_b) plus one per boundary
        MR, in detection order.
    charges : list[dict]
        One charge dict per surviving crossing plus one per boundary MR
        (``kind='mr'``, charge 0).
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
        la = np.log(np.abs(seg.tracked_roots))           # (N, K)
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
                    tol=crossing_tol, xtol=crossing_tol, max_iter=max_iter,
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
    # row; here we only materialise it.  Re-locating against the (now
    # mutation-rich) mesh would pair the charge with the wrong β₂.
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
