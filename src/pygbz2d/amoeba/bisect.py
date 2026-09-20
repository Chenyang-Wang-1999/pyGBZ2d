"""
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-14
Copyright © Department of Physics, Tsinghua University. All rights reserved

Bisection algorithms for amoeba Ronkin function critical-point search.

The μ₂ bisection is now continuum-first:

  1. ``_try_fast_mu2``      — cheap gap test at θ₁=0 (no continuum involved).
  2. ``detect_continuum``   — std-based flat-track detection, run up front.
  3. continuum probes       — for each merged μ₂_c, perturb by ±ε and test
     whether w2 straddles zero.  If yes the inner solve returns the continuum
     boundary; if no, the signed probe points tighten the μ₂ bracket.
  4. ``_bisect_mu2_discrete`` — plain discrete μ₂ bisection on the tightened
     bracket, using ``calculate_a2_average_winding`` (which is
     continuum-unaware).

The outer μ₁ bisection receives ``is_continuum=True`` results and resolves
them with μ₁ ± ε perturbations: opposite w1 signs end the search; equal signs
update the μ₁ bracket.
"""

from typing import Optional
import warnings

import numpy as np
from cmath import exp

from ..core import CharPoly
from ..continuation.interpolation import hermite_interp_poly

from .. import core
from ..core import live_defaults

from .ronkin_winding import _get_average_winding_from_zeros
from .zm_extract import (
    AmoebaZeroManager,
    _calculate_a2_winding_and_zeros,
    calculate_a2_average_winding,
    detect_continuum,
    find_crossings,
)

# μ₂ bisection budget (the winding is monotonic in μ₂, so range
# expansion guarantees a sign bracket eventually).
BISECT_MAX_ITER: int = 60
#: Coarse-stage tolerance for the first, unrefined μ₂ bisection.
BISECT_COARSE_XTOL: float = 1e-3
MAX_RANGE_EXPANSIONS: int = 10
RANGE_EXPAND_FACTOR: float = 2.0
#: Relative tolerance for mesh-insertion dedup in extremum refinement.
EXTREMUM_INSERT_REL_TOL: float = 1e-12


# ---------------------------------------------------------------------------
# Extremum refinement for _try_fast_mu2 (Hermite prediction + mesh insertion)
# ---------------------------------------------------------------------------

def _screen_extremum_intervals(
    zm: AmoebaZeroManager,
    cols: np.ndarray,
    mode: str,
) -> list[tuple[int, int, int]]:
    """Vectorized first screen for intervals that may contain an extremum.

    Returns ``[(seg_idx, i, col), ...]`` where ``i`` is the left endpoint of a
    suspicious interval ``[i, i+1]`` for track ``col``.  Only these intervals
    are later processed with cubic-Hermite root solving (which is not
    vectorizable and therefore must be kept small).

    Suspicious intervals for ``mode='max'``:
      - derivative sign changes from + to − across the interval;
      - either endpoint derivative is exactly 0;
      - either endpoint derivative is non-finite (divergent at an MR row).
    ``mode='min'`` uses − to +.  Intervals adjacent to a non-finite derivative
    row are always marked (the later refinement inserts a midpoint separator
    there to resolve the sharp feature with an exact solve).
    """
    out: list[tuple[int, int, int]] = []
    cols_arr = np.asarray(cols, dtype=int)
    for s, seg in enumerate(zm.segments):
        n = len(seg.theta1_arr)
        if n < 2:
            continue
        la = zm.seg_logabs[s][:, cols_arr]
        d = seg.tangents.real[:, cols_arr]          # (N, C)

        finite = np.isfinite(d)
        bad = ~finite                                # (N, C)

        if mode == 'max':
            sign_change = (d[:-1] > 0) & (d[1:] < 0)
        else:
            sign_change = (d[:-1] < 0) & (d[1:] > 0)

        zero_deriv = (d[:-1] == 0) | (d[1:] == 0)
        # Any interval adjacent to a non-finite derivative row is suspicious:
        # the Hermite fallback inserts a midpoint separator there so the sharp
        # feature gets resolved by an exact solve.  This single vectorized OR
        # covers all bad-row cases (interior, first, last) — no per-row loop
        # is needed.
        bad_interval = bad[:-1] | bad[1:]
        susp = sign_change | zero_deriv | bad_interval   # (N-1, C)

        rows, cidx = np.where(susp)
        for i, c in zip(rows, cidx):
            out.append((s, int(i), int(cols_arr[int(c)])))
    return out


def _refine_track_extrema(
    zm: AmoebaZeroManager,
    cols: np.ndarray,
    mode: str,
) -> None:
    """Insert mesh rows at predicted extrema of selected tracks.

    ``mode='max'`` inserts predicted maxima; ``mode='min'`` inserts predicted
    minima.  The heavy Hermite-polynomial work is only done on intervals that
    survived :func:`_screen_extremum_intervals`.
    """
    cols_arr = np.asarray(cols, dtype=int)
    suspicious = _screen_extremum_intervals(zm, cols_arr, mode)
    if not suspicious:
        return

    # Current global extreme over all rows (skipping the duplicated θ=2π seam
    # row of the last segment) — used as a cheap filter: candidates that cannot
    # improve the current extreme are not inserted.
    current = _extreme_over_segments(zm, cols_arr, mode)
    if not np.isfinite(current):
        current = -np.inf if mode == 'max' else np.inf

    candidates: list[tuple[float, int, float]] = []  # (theta, seg_idx, pred)
    for s, i, j in suspicious:
        seg = zm.segments[s]
        th = seg.theta1_arr
        la = zm.seg_logabs[s]
        h = float(th[i + 1] - th[i])
        if not np.isfinite(h) or h <= 0.0:
            continue
        v0 = float(la[i, j])
        v1 = float(la[i + 1, j])
        d0 = float(seg.tangents.real[i, j])
        d1 = float(seg.tangents.real[i + 1, j])

        finite_deriv = bool(np.isfinite(d0) and np.isfinite(d1))
        poly = hermite_interp_poly(h, v0, d0, v1, d1)

        if not finite_deriv:
            # Divergent tangent: the Hermite falls back to a linear piece with
            # no interior extremum.  Insert the interval midpoint as a mesh
            # separator so the sharp feature is resolved by an exact solve.
            theta = float(th[i] + 0.5 * h)
            pred = float(np.polyval(poly, 0.5 * h))
            candidates.append((theta, s, pred))
            continue

        deriv = np.polyder(poly)
        for root in np.roots(deriv):
            if abs(float(root.imag)) > 1e-12 * max(1.0, h):
                continue
            x = float(root.real)
            if not (0.0 < x < h):
                continue
            pred = float(np.polyval(poly, x))
            if mode == 'max' and pred <= current:
                continue
            if mode == 'min' and pred >= current:
                continue
            candidates.append((float(th[i] + x), s, pred))

    if not candidates:
        return

    # Insert descending per segment so earlier insertions (larger θ) do not
    # invalidate the located intervals of later ones.  Dedup against existing
    # mesh rows with a per-interval relative tolerance, mirroring the SGBZ
    # refinement-grid insertion.
    candidates.sort(key=lambda x: (x[1], -x[0]))
    for s, grp in _group_by_segment(candidates):
        seg = zm.segments[s]
        th = np.asarray(seg.theta1_arr, dtype=float)
        inserted: list[float] = []
        for theta, _, _ in sorted(grp, key=lambda x: -x[0]):
            tol = EXTREMUM_INSERT_REL_TOL * max(1.0, float(theta))
            if np.any(np.abs(th - theta) <= tol):
                continue
            if inserted and min(abs(theta - t) for t in inserted) <= tol:
                continue
            try:
                _, changed = zm.insert_solution(theta, seg_idx=s, interp='hermite')
            except Exception as exc:
                warnings.warn(
                    f"extremum refinement insert failed at theta="
                    f"{theta:.6e} in segment {s}: {exc}"
                )
                continue
            if changed:
                inserted.append(theta)
    zm.refresh_logabs()


def _group_by_segment(candidates):
    """Group (theta, seg_idx, pred) candidates by segment index."""
    by_seg: dict[int, list] = {}
    for theta, s, pred in candidates:
        by_seg.setdefault(s, []).append((theta, s, pred))
    return list(by_seg.items())


def _extreme_over_segments(
    zm: AmoebaZeroManager,
    cols: np.ndarray,
    mode: str,
) -> float:
    """Current global extreme of selected track columns over all segments.

    The last segment's final row is the θ=2π copy of the θ=0 row, but its
    columns are permuted by ``boundary_perm``; it is skipped here because the
    θ=0 row already covers that physical point in the correct frame.
    """
    parts: list[np.ndarray] = []
    n_seg = len(zm.segments)
    for s, seg in enumerate(zm.segments):
        la = zm.seg_logabs[s]
        end = len(seg.theta1_arr)
        if s == n_seg - 1 and end > 1:
            end -= 1
        if end <= 0:
            continue
        parts.append(la[:end, cols])
    if not parts:
        return float('nan')
    all_vals = np.concatenate(parts)
    return float(np.max(all_vals) if mode == 'max' else np.min(all_vals))


# ---------------------------------------------------------------------------
# Step 2: cheap gap test at θ₁ = 0
# ---------------------------------------------------------------------------

def _try_fast_mu2(
    zm: AmoebaZeroManager,
    char_poly: CharPoly,
) -> dict:
    """Try the simple μ₂ gap criterion.

    At θ₁ = 0, sort the roots by |β₂|.  Let ``lo_cols`` be the first M
    columns and ``hi_cols`` the remaining N columns.  If the largest ln|β₂|
    over all ``lo_cols`` tracks is strictly below the smallest ln|β₂| over all
    ``hi_cols`` tracks, then any μ₂ between the two bounds has exactly M roots
    below it everywhere, so w2 = 0 without bisection.

    Returns ``{ok, A, B, lo_cols, hi_cols}``.  On failure, A and B are still
    returned so the caller can use them to tighten the discrete μ₂ bracket.
    """
    M, N = char_poly.get_minor_degrees(2)
    K = zm.K
    if not (0 < M < K):
        return {"ok": False, "A": float('nan'), "B": float('nan'),
                "lo_cols": np.array([], dtype=int),
                "hi_cols": np.array([], dtype=int)}

    if len(zm.segments) == 0 or len(zm.segments[0].theta1_arr) == 0:
        return {"ok": False, "A": float('nan'), "B": float('nan'),
                "lo_cols": np.array([], dtype=int),
                "hi_cols": np.array([], dtype=int)}

    roots0 = zm.segments[0].tracked_roots[0, :]
    order = np.argsort(np.abs(roots0))
    lo_cols = order[:M].astype(int)
    hi_cols = order[M:].astype(int)

    # Refine the mesh near the extrema that determine the gap.
    _refine_track_extrema(zm, lo_cols, 'max')
    _refine_track_extrema(zm, hi_cols, 'min')

    A = _extreme_over_segments(zm, lo_cols, 'max')
    B = _extreme_over_segments(zm, hi_cols, 'min')

    ok = bool(np.isfinite(A) and np.isfinite(B) and A < B)
    return {"ok": ok, "A": float(A), "B": float(B),
            "lo_cols": lo_cols, "hi_cols": hi_cols}


# ---------------------------------------------------------------------------
# Pure discrete μ₂ bisection
# ---------------------------------------------------------------------------

def _bisect_mu2_discrete(
    char_poly: CharPoly,
    E_ref: complex,
    mu1: float,
    mu2_low: float,
    mu2_high: float,
    max_iter: int,
    wtol: float,
    coarse_xtol: float,
    max_range_expansions: int,
    range_expand_factor: float,
    _zm: AmoebaZeroManager,
) -> dict:
    """Two-stage discrete μ₂ bisection on a continuum-free bracket.

    The caller has already run :func:`detect_continuum` and tightened the
    bracket away from flat tracks, so this function is continuum-unaware.

    Stage 1: coarse bisection with ``return_refined=False`` and tolerance
    ``coarse_xtol`` — no crossing refinement. Stage 2: fine bisection on the
    coarse bracket with ``return_refined=True`` and tolerance ``wtol``. Root
    rows inserted during refinement remain in the ZM for later μ₂ values.
    """
    zm = _zm
    low, high = float(mu2_low), float(mu2_high)

    def _winding_at(mu2_val, return_refined=False):
        return calculate_a2_average_winding(
            zm, mu1, mu2_val, return_refined=return_refined,
        )

    # Adaptive range expansion (coarse, no refine)
    for _ in range(max_range_expansions):
        w_low = _winding_at(low)
        w_high = _winding_at(high)

        if w_low * w_high <= 0:
            break

        width = high - low
        low = low - range_expand_factor * width
        high = high + range_expand_factor * width
    else:
        raise ValueError(f"Maximal range expansion reached! E_ref:{E_ref}, mu1:{mu1}, low:{low}, high:{high}")

    # Stage 1: coarse bisection, no refine.
    for _ in range(max_iter):
        mu2_mid = 0.5 * (low + high)
        w_mid = _winding_at(mu2_mid)

        if (high - low) < coarse_xtol:
            break

        if w_low * w_mid <= 0:
            high = mu2_mid
            w_high = w_mid
        else:
            low = mu2_mid
            w_low = w_mid
    else:
        mu2_final = 0.5 * (low + high)
        raise RuntimeError(
            f"mu2 coarse bisection failed to converge after {max_iter} "
            f"iterations: E_ref={E_ref}, mu1={mu1}, mu2={mu2_final:.12g}, "
            f"bracket=({low:.12g}, {high:.12g})"
        )

    # Stage 2: fine bisection with refined crossings.  Refinement may move
    # the winding-zero slightly, so the coarse bracket's endpoints are not
    # guaranteed to straddle after re-evaluating with ``return_refined=True``.
    # Expand the fine bracket outward until the refined windings straddle
    # (w2 is monotone in μ₂, so the required direction is known).
    low_fine, high_fine = low, high
    for _ in range(max_range_expansions):
        w_low = _winding_at(low_fine, return_refined=True)
        w_high = _winding_at(high_fine, return_refined=True)
        if w_low * w_high <= 0:
            break
        if w_low > 0:
            low_fine -= coarse_xtol
        else:
            high_fine += coarse_xtol
    else:
        raise RuntimeError(
            f"fine μ₂ bracket expansion failed: E_ref={E_ref}, mu1={mu1}, "
            f"coarse bracket=({low:.12g}, {high:.12g})"
        )

    for _ in range(max_iter):
        mu2_mid = 0.5 * (low_fine + high_fine)
        w_mid, zeros = _calculate_a2_winding_and_zeros(
            zm, mu1, mu2_mid, return_refined=True)

        # SGBZ-parity exit rule (sgbz_solver.py:387): drop the bracket-width
        # exit, keep only the winding-zero exit.  Both halves of the old
        # condition were needed to *terminate*, but only `abs(w_mid) < wtol`
        if abs(w_mid) < wtol:
            return {
                "mu2": mu2_mid,
                "zeros": zeros,
                "is_continuum": False,
                "winding": w_mid,
                "_zm": zm,
            }

        if w_low * w_mid <= 0:
            high_fine = mu2_mid
            w_high = w_mid
        else:
            low_fine = mu2_mid
            w_low = w_mid

    mu2_final = 0.5 * (low_fine + high_fine)
    raise RuntimeError(
        f"mu2 fine bisection failed to converge after {max_iter} iterations: "
        f"E_ref={E_ref}, mu1={mu1}, mu2={mu2_final:.12g}, winding={w_mid}"
    )


# ---------------------------------------------------------------------------
# Inner μ₂ solve: continuum-first dispatcher
# ---------------------------------------------------------------------------

@live_defaults(continuum_tol="core:CONTINUUM_TOL", continuum_perturb="core:CONTINUUM_PERTURB",
               max_iter="amoeba.bisect:BISECT_MAX_ITER", wtol="core:WINDING_ZERO_TOL",
               coarse_xtol="amoeba.bisect:BISECT_COARSE_XTOL",
               max_range_expansions="amoeba.bisect:MAX_RANGE_EXPANSIONS",
               range_expand_factor="amoeba.bisect:RANGE_EXPAND_FACTOR")
def _find_mu2_for_w2_zero(
    char_poly: CharPoly,
    E_ref: complex,
    mu1: float,
    mu2_low: float,
    mu2_high: float,
    continuum_tol: Optional[float] = None,
    continuum_perturb: Optional[float] = None,
    max_iter: Optional[int] = None,
    wtol: Optional[float] = None,
    coarse_xtol: Optional[float] = None,
    max_range_expansions: Optional[int] = None,
    range_expand_factor: Optional[float] = None,
    _zm: Optional[AmoebaZeroManager] = None,
    frac: Optional[float] = None,
) -> dict:
    """Find μ₂ where w2 = 0, continuum-first.

    Flow:
      1. build / reuse the μ₂-independent ``AmoebaZeroManager``;
      2. ``_try_fast_mu2`` — a wide μ₂ gap makes the answer trivial;
      3. ``detect_continuum`` — std-based flat-track groups;
      4. for each merged μ₂_c, probe w2 at μ₂_c ± ε:
           opposite signs → return the continuum boundary;
           equal signs     → tighten the μ₂ bracket with the signed probes;
      5. run ``_bisect_mu2_discrete`` on the tightened bracket.
    """
    if _zm is None:
        zm = AmoebaZeroManager(char_poly, E_ref, mu1)
        zm.run()
    else:
        zm = _zm

    # ---- step 2: cheap gap test ----
    fast = _try_fast_mu2(zm, char_poly)
    if fast["ok"]:
        return {
            "mu2": 0.5 * (fast["A"] + fast["B"]),
            "zeros": [],
            "is_continuum": False,
            "winding": 0.0,
            "_zm": zm,
            "_fast": True,
            "_fast_A": fast["A"],
            "_fast_B": fast["B"],
        }

    # ---- bracket seed from the fast test's A / B ----
    low, high = float(mu2_low), float(mu2_high)
    A, B = fast["A"], fast["B"]
    if np.isfinite(A) and np.isfinite(B) and A >= B:
        # fast path failed ⇒ A >= B; [min(A,B), max(A,B)] = [B, A] brackets
        # the w2 sign change.
        if B < A and (A - B) > wtol:
            low, high = B, A

    # ---- continuum detection (std, up front) ----
    groups = detect_continuum(zm, continuum_tol)

    if groups:
        for mu2_c, members in groups:
            w_left = w_right = None
            for scale in core.ESCAPE_LADDER:
                eps = continuum_perturb * scale
                w_left = calculate_a2_average_winding(
                    zm, mu1, mu2_c - eps,
                )
                w_right = calculate_a2_average_winding(
                    zm, mu1, mu2_c + eps,
                )
                if np.isfinite(w_left) and np.isfinite(w_right):
                    break
            else:
                raise RuntimeError(
                    f"continuum at mu2_c={mu2_c:.8g} (E_ref={E_ref}, "
                    f"mu1={mu1}) could not be resolved: all perturbation "
                    f"scales still degenerate."
                )

            if w_left * w_right < 0:
                # Continuum boundary: w2 changes sign across the flat track.
                return {
                    "mu2": mu2_c,
                    "zeros": [],
                    "is_continuum": True,
                    "winding": (w_left, w_right),
                    "_continuum_members": members,
                    "_zm": zm,
                }

            # w2 not opposite: signed probe.  w2 is monotone increasing in μ₂,
            # so both positive ⇒ the zero lies below the flat band; both
            # negative ⇒ it lies above.
            if w_left > 0:
                high = min(high, mu2_c - eps)
            else:
                low = max(low, mu2_c + eps)

    # ---- bracket sanity after continuum tightening ----
    if not (np.isfinite(low) and np.isfinite(high)) or low >= high:
        low, high = float(mu2_low), float(mu2_high)

    return _bisect_mu2_discrete(
        char_poly, E_ref, mu1, low, high,
        max_iter, wtol, coarse_xtol,
        max_range_expansions, range_expand_factor,
        _zm=zm,
    )


# ---------------------------------------------------------------------------
# Continuum handling for the outer μ₁ bisection
# ---------------------------------------------------------------------------

def _handle_continuum(
    char_poly: CharPoly,
    E_ref: complex,
    mu1_val: float,
    mu2_c: float,
    *,
    continuum_perturb: float,
) -> dict:
    """Resolve a continuum inner result by perturbing μ₁.

    Keeps ``μ₂`` fixed at the continuum level.  At each of ``mu1_val ± eps``
    a fresh ZM is built and ``find_crossings`` locates the zeros of
    ``ln|β₂| = μ₂``; those zeros are then fed to the a1 average-winding
    integral.  Returns:

      - ``is_boundary``: w1_left and w1_right have opposite signs — this μ₁ is
        a Ronkin-minimum continuum boundary.
      - ``w1_left`` / ``w1_right``: the two limits.
    """
    eps = continuum_perturb

    def _w1_at(mu1_probe: float) -> float:
        zm = AmoebaZeroManager(char_poly, E_ref, mu1_probe)
        zm.run()
        crossings = find_crossings(zm, mu1_probe, mu2_c)
        zeros = [
            (float(np.angle(b1) % (2.0 * np.pi)),
             float(np.angle(b2) % (2.0 * np.pi)))
            for b1, b2 in crossings
        ]
        w1, _ = _get_average_winding_from_zeros(
            char_poly, E_ref, mu1_probe, mu2_c, zeros, direction=1,
        )
        return w1

    w1_left = _w1_at(mu1_val - eps)
    w1_right = _w1_at(mu1_val + eps)

    return {
        "is_boundary": bool(w1_left * w1_right < 0),
        "w1_left": w1_left,
        "w1_right": w1_right,
        "mu2_c": mu2_c,
    }


# ---------------------------------------------------------------------------
# Outer μ₁ bisection for the Ronkin minimum
# ---------------------------------------------------------------------------

@live_defaults(continuum_tol="core:CONTINUUM_TOL", continuum_perturb="core:CONTINUUM_PERTURB",
               max_iter="amoeba.bisect:BISECT_MAX_ITER", wtol="core:WINDING_ZERO_TOL",
               max_range_expansions="amoeba.bisect:MAX_RANGE_EXPANSIONS",
               range_expand_factor="amoeba.bisect:RANGE_EXPAND_FACTOR")
def bisect_amoeba_ronkin_min(
    char_poly: CharPoly,
    E_ref: complex,
    mu1_low: float = -1,
    mu1_high: float = 1,
    mu2_low: float = -1,
    mu2_high: float = 1,
    continuum_tol: Optional[float] = None,
    continuum_perturb: Optional[float] = None,
    max_iter: Optional[int] = None,
    wtol: Optional[float] = None,
    max_range_expansions: Optional[int] = None,
    range_expand_factor: Optional[float] = None,
    frac: Optional[float] = None,
) -> dict:
    """Find the Ronkin function minimum by bisecting μ₁ and μ₂.

    Outer loop: bisect μ₁.
    Inner loop: ``_find_mu2_for_w2_zero`` (continuum-first).

    When the inner solve returns a continuum boundary, the outer loop resolves
    it with μ₁ ± ε perturbations: opposite w1 signs end the search; equal
    signs update the μ₁ bracket.
    """
    # Evaluate w1 at the mu1 endpoints, with adaptive range expansion
    low, high = float(mu1_low), float(mu1_high)
    w1_low = w1_high = 0.0
    inner_low = inner_high = None

    for _ in range(max_range_expansions):
        inner_low = _find_mu2_for_w2_zero(
            char_poly, E_ref, low, mu2_low, mu2_high,
            continuum_tol=continuum_tol,
            continuum_perturb=continuum_perturb, max_iter=max_iter, wtol=wtol,
            max_range_expansions=max_range_expansions,
            range_expand_factor=range_expand_factor,
        )
        inner_high = _find_mu2_for_w2_zero(
            char_poly, E_ref, high, mu2_low, mu2_high,
            continuum_tol=continuum_tol,
            continuum_perturb=continuum_perturb, max_iter=max_iter, wtol=wtol,
            max_range_expansions=max_range_expansions,
            range_expand_factor=range_expand_factor,
        )

        # Every inner return may be a continuum boundary.  For an endpoint we
        # only need a scalar w1 sign; if the endpoint itself is a w1 boundary,
        # it is already the solution.
        if inner_low["is_continuum"]:
            hc_low = _handle_continuum(
                char_poly, E_ref, low, inner_low["mu2"],
                continuum_perturb=continuum_perturb,
            )
            if hc_low["is_boundary"]:
                return {
                    "mu1": low, "mu2": hc_low["mu2_c"], "zeros": [],
                    "is_continuum": True,
                    "_mu1_bracket": (low, high),
                    "_w1_bracket": (hc_low["w1_left"], hc_low["w1_right"]),
                    "_exit_reason": "continuum_boundary",
                    "_continuum_members": inner_low.get("_continuum_members"),
                    "_zm": inner_low.get("_zm"),
                }
            w1_low = hc_low["w1_right"]  # interior-side limit
        else:
            w1_low, _ = _get_average_winding_from_zeros(
                char_poly, E_ref, low, inner_low["mu2"],
                inner_low["zeros"], direction=1,
            )

        if inner_high["is_continuum"]:
            hc_high = _handle_continuum(
                char_poly, E_ref, high, inner_high["mu2"],
                continuum_perturb=continuum_perturb,
            )
            if hc_high["is_boundary"]:
                return {
                    "mu1": high, "mu2": hc_high["mu2_c"], "zeros": [],
                    "is_continuum": True,
                    "_mu1_bracket": (low, high),
                    "_w1_bracket": (hc_high["w1_left"], hc_high["w1_right"]),
                    "_exit_reason": "continuum_boundary",
                    "_continuum_members": inner_high.get("_continuum_members"),
                    "_zm": inner_high.get("_zm"),
                }
            w1_high = hc_high["w1_left"]  # interior-side limit
        else:
            w1_high, _ = _get_average_winding_from_zeros(
                char_poly, E_ref, high, inner_high["mu2"],
                inner_high["zeros"], direction=1,
            )

        if w1_low * w1_high <= 0:
            break

        width = high - low
        low = low - range_expand_factor * width
        high = high + range_expand_factor * width
    else:
        raise ValueError(f"bisect_amoeba_ronkin_min: reach maximum expansion. E_ref:{E_ref}, low:{low}, high:{high}")

    mu1_low, mu1_high = low, high
    inner_mid = None
    w1_mid = None

    for _ in range(max_iter):
        mu1_mid = 0.5 * (mu1_low + mu1_high)

        # Build the ZeroManager once for this mu1_mid — reused across all
        # μ₂ evaluations in the inner bisection.
        zm = AmoebaZeroManager(char_poly, E_ref, mu1_mid)
        zm.run()

        inner_mid = _find_mu2_for_w2_zero(
            char_poly, E_ref, mu1_mid, mu2_low, mu2_high,
            continuum_tol=continuum_tol,
            continuum_perturb=continuum_perturb, max_iter=max_iter, wtol=wtol,
            max_range_expansions=max_range_expansions,
            range_expand_factor=range_expand_factor, _zm=zm,
        )

        mu2_mid = inner_mid["mu2"]

        # ---- continuum path ----
        if inner_mid["is_continuum"]:
            hc = _handle_continuum(
                char_poly, E_ref, mu1_mid, mu2_mid,
                continuum_perturb=continuum_perturb,
            )

            if hc["is_boundary"]:
                # Both axes straddle zero → Ronkin minimum at the continuum
                # boundary.  Outer bisection ends here.
                return {
                    "mu1": mu1_mid, "mu2": hc["mu2_c"], "zeros": [],
                    "is_continuum": True,
                    "_mu1_bracket": (mu1_low, mu1_high),
                    "_w1_bracket": (w1_low, w1_high),
                    "_exit_reason": "continuum_boundary",
                    "_continuum_members": inner_mid.get("_continuum_members"),
                    "_zm": zm,
                }

            # w1 not opposite: update the μ₁ bracket by the sign of w1.
            if hc["w1_left"] >= 0:
                mu1_high = mu1_mid - continuum_perturb
                w1_high = hc["w1_left"]
            else:
                mu1_low = mu1_mid + continuum_perturb
                w1_low = hc["w1_right"]
            continue

        # ---- normal (discrete) path ----
        zeros_mid = inner_mid["zeros"]

        # Compute w1 at (mu1_mid, mu2_mid), reusing zeros from the
        # inner bisection — no extra root tracking needed.
        w1_mid, w1_area = _get_average_winding_from_zeros(
            char_poly, E_ref, mu1_mid, mu2_mid,
            zeros_mid, direction=1,
        )

        # SGBZ-parity exit rule (sgbz_solver.py:387, 2026-08-15): the bracket
        # width is NOT an exit reason. Bisection now ends
        # only on a genuine winding zero, matching SGBZ.
        if abs(w1_mid) < wtol:
            return {
                "mu1": mu1_mid,
                "mu2": mu2_mid,
                "zeros": zeros_mid,
                "is_continuum": inner_mid["is_continuum"],
                "_mu1_bracket": (mu1_low, mu1_high),
                "_w1_bracket": (w1_low, w1_high),
                "_exit_reason": "w1_zero",
                "_w1_area": w1_area,
                "_zm": zm,
            }

        if w1_low * w1_mid < 0:
            mu1_high = mu1_mid
            w1_high = w1_mid
        else:
            mu1_low = mu1_mid
            w1_low = w1_mid

    # Max iterations exhausted — the outer μ₁ bisection failed to converge.
    mu1_final = 0.5 * (mu1_low + mu1_high)
    zeros_final = (inner_mid.get("zeros") or []) if inner_mid is not None else []
    raise RuntimeError(
        f"outer mu1 bisection failed to converge after {max_iter} iterations: "
        f"E_ref={E_ref}, mu1={(mu1_high, mu1_low)}, mu2={mu2_mid}, len(subsets)={len(zeros_final)}, "
        f"is_continuum={inner_mid['is_continuum'] if inner_mid is not None else None}, "
        f"winding={w1_high}, {w1_low}"
    )
