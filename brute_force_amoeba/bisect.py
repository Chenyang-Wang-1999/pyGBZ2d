"""
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-14
Copyright © Department of Physics, Tsinghua University. All rights reserved

Bisection algorithms for amoeba Ronkin function critical-point search.

The winding source is now ``continuation.ZeroManager`` (via
``AmoebaZeroManager`` + ``zm_extract.amoeba_windings``) rather than the old
fixed-grid Hungarian tracks.  Each ``(E, μ₁)`` builds one
``AmoebaZeroManager`` and reuses it across all μ₂ evaluations (~30×).

Provides μ₂ bisection (w2 winding crossing) and the outer μ₁/μ₂ bisection
for the Ronkin minimum.
"""

from typing import Optional
from gbz_types import CharPoly

from .ronkin_winding import _get_average_winding_from_zeros
from .zm_extract import AmoebaZeroManager, amoeba_windings, CONTINUUM_TOL, CONTINUUM_FRAC


def _refine_and_correct(
    char_poly, E_ref, mu1, mu2_0,
    zm, low, high, low_init, high_init,
    continuum_tol, min_continuum_pts_unused, xtol,
    frac,
):
    """Refine crossings at mu2_0, then apply Newton correction if needed.

    Uses the analytical dW/dmu2 computed from zero derivatives to take one
    Newton step, avoiding re-bisection after refinement.
    """
    w_ref, zeros_ref, has_cont, dW_dmu2 = amoeba_windings(
        zm, char_poly, E_ref, mu1, mu2_0,
        tol=continuum_tol, frac=frac, refine=True,
    )

    # Defensive: refined detection may find continuum that the unrefined
    # bisection missed (e.g. at band edges).  Return empty zeros and let
    # the caller handle the continuum via perturbation.
    if has_cont:
        return {
            "mu2": mu2_0, "zeros": [], "is_continuum": True,
            "winding": 0.0, "_zm": zm,
        }

    if abs(w_ref) < xtol:
        return {
            "mu2": mu2_0, "zeros": zeros_ref, "is_continuum": False,
            "winding": w_ref, "_zm": zm,
        }

    # Newton correction using analytical derivative
    if abs(dW_dmu2) > 1e-15:
        delta = -w_ref / dW_dmu2
        # Clamp to initial bracket with margin
        bracket_width = high_init - low_init
        max_step = 0.5 * bracket_width
        delta = max(-max_step, min(max_step, delta))
        mu2_new = mu2_0 + delta
        # Ensure within initial bracket
        mu2_new = max(low_init, min(high_init, mu2_new))
    else:
        # dW/dmu2 near zero — use bisection bracket midpoint
        mu2_new = 0.5 * (low + high)

    # Avoid re-refining at essentially the same point
    if abs(mu2_new - mu2_0) < xtol:
        return {
            "mu2": mu2_0, "zeros": zeros_ref, "is_continuum": False,
            "winding": w_ref, "_zm": zm,
        }

    # Refine at the corrected mu2
    w_ref2, zeros_ref2, _, _ = amoeba_windings(
        zm, char_poly, E_ref, mu1, mu2_new,
        tol=continuum_tol, frac=frac, refine=True,
    )

    if abs(w_ref2) < xtol:
        return {
            "mu2": mu2_new, "zeros": zeros_ref2, "is_continuum": False,
            "winding": w_ref2, "_zm": zm,
        }

    # Fall back: return the better of the two refined points
    if abs(w_ref2) < abs(w_ref):
        return {
            "mu2": mu2_new, "zeros": zeros_ref2, "is_continuum": False,
            "winding": w_ref2, "_zm": zm,
        }
    return {
        "mu2": mu2_0, "zeros": zeros_ref, "is_continuum": False,
        "winding": w_ref, "_zm": zm,
    }


def _find_mu2_for_w2_zero(
    char_poly: CharPoly,
    E_ref: complex,
    mu1: float,
    mu2_low: float,
    mu2_high: float,
    N_points: int = 301,
    continuum_tol: float = 1e-6,
    min_continuum_pts: int = 3,
    continuum_perturb: float = 1e-4,
    max_iter: int = 60,
    xtol: float = 1e-10,
    max_range_expansions: int = 10,
    range_expand_factor: float = 2.0,
    _zm: Optional[AmoebaZeroManager] = None,
    frac: float = CONTINUUM_FRAC,
) -> dict:
    """Find mu2 where w2 winding crosses 0, with adaptive range.

    The w2 winding is monotonic in mu2, so expanding the search range
    guarantees a sign change will eventually be found.

    Uses unrefined (cheap) winding during bisection, then refines crossings
    at the final mu2 and applies a Newton correction using the analytical
    derivative dW/dmu2.

    _zm: optionally pre-built AmoebaZeroManager at (E, mu1).  When provided,
    avoids rebuilding root tracks (the ZM is mu2-independent).  When None,
    a fresh one is built here.
    """
    low, high = float(mu2_low), float(mu2_high)

    # Build ZeroManager once (independent of mu2).  Reused across ~30 mu2
    # evaluations — the caching role the old `tracks` dict played.
    if _zm is None:
        zm = AmoebaZeroManager(char_poly, E_ref, mu1)
        zm.run()
    else:
        zm = _zm

    def _winding_at(mu2_val, refine=False):
        """Evaluate winding at mu2_val.  Strips dW_dmu2 for the bisection loop."""
        w, z, c, _ = amoeba_windings(
            zm, char_poly, E_ref, mu1, mu2_val,
            tol=continuum_tol, frac=frac, refine=refine,
        )
        return w, z, c

    # Adaptive range expansion (unrefined)
    for _ in range(max_range_expansions):
        w_low, _, has_cont_low = _winding_at(low)
        w_high, _, has_cont_high = _winding_at(high)

        # If either endpoint encounters a continuum subset the winding
        # value is None — expand the range and retry.  The tolerance is
        # O(continuum_tol), so a single expansion typically escapes it.
        if has_cont_low or has_cont_high:
            width = high - low
            low = low - range_expand_factor * width
            high = high + range_expand_factor * width
            continue

        if w_low * w_high <= 0:
            low_init, high_init = low, high
            break

        # Expand outward
        width = high - low
        low = low - range_expand_factor * width
        high = high + range_expand_factor * width
    else:
        raise ValueError(f"Maximal range expansion reached! E_ref:{E_ref}, mu1:{mu1}, low:{low}, high:{high}")

    # Bisection with unrefined winding (fast — no fsolve)
    for _ in range(max_iter):
        mu2_mid = 0.5 * (low + high)

        w_mid, zeros, has_continuum = _winding_at(mu2_mid)

        if has_continuum:
            w_left, _, _ = _winding_at(mu2_mid - continuum_perturb)
            w_right, _, _ = _winding_at(mu2_mid + continuum_perturb)

            if w_left * w_right < 0:

                # At the continuum mu2 there are no discrete zeros
                # (every point in the band satisfies |β2| ≈ exp(mu2)).
                # Return zeros=[] so the outer bisection knows to compute
                # w1 limits via mu1 perturbation rather than from zeros.
                return {
                    "mu2": mu2_mid, "zeros": [], "is_continuum": True,
                    "winding": (w_left, w_right), "_zm": zm,
                }
            else:
                if w_low * w_left > 0:
                    low = mu2_mid
                    w_low = w_left
                else:
                    high = mu2_mid
                    w_high = w_left
                continue

        if abs(w_mid) < xtol or (high - low) < xtol:
            # --- Post-refinement + Newton correction ---
            return _refine_and_correct(
                char_poly, E_ref, mu1, mu2_mid,
                zm, low, high, low_init, high_init,
                continuum_tol, min_continuum_pts, xtol, frac,
            )

        if w_low * w_mid < 0:
            high = mu2_mid
            w_high = w_mid
        else:
            low = mu2_mid
            w_low = w_mid

    # Max iterations reached
    raise ValueError("Max iterations reached in mu2 bisection.")


def _resolve_continuum(
    char_poly: CharPoly,
    E_ref: complex,
    mu1: float,
    mu2: float,
    zm: AmoebaZeroManager,
    mu2_low: float = -1.0,
    mu2_high: float = 1.0,
    continuum_perturb: float = 1e-4,
    N_points: int = 301,
    continuum_tol: float = 1e-6,
    min_continuum_pts: int = 3,
    max_iter: int = 60,
    xtol: float = 1e-10,
    max_range_expansions: int = 10,
    range_expand_factor: float = 2.0,
    frac: float = CONTINUUM_FRAC,
) -> dict:
    """Resolve a continuum point by computing winding left/right limits.

    Step 1 — a2 axis: perturb mu2 ± ε, compute w2 winding at each.
    If the two values have opposite signs, the mu2 jump crosses 0 →
    this mu2 is the w2=0 boundary.

    Step 2 — a1 axis (only when Step 1 succeeds): perturb mu1 ± ε,
    re-run the full inner mu2 bisection at each perturbed mu1, then
    compute w1 winding from the resulting (mu2, zeros).  If the two
    w1 values straddle zero, this (mu1, mu2) is the Ronkin minimum.

    Returns a dict with keys:
        w2_left:      float | None
        w2_right:     float | None
        w2_opposite:  bool
        w1_left:      float | None
        w1_right:     float | None
        w1_opposite:  bool | None
        is_boundary:  bool
        w1_resolved:  bool
    """
    # ---- Step 1: w2 limits via mu2 perturbation ----
    # The continuum band is O(continuum_tol), so continuum_perturb
    # (default 1e-4) should escape it.  If the perturbed mu2 still
    # falls in the band, expand the perturbation and retry.
    w2_left = None
    w2_right = None
    w2_opposite = False

    for scale in (1.0, 2.0, 4.0, 8.0):
        eps = continuum_perturb * scale

        w_left, _, has_cont_left, _ = amoeba_windings(
            zm, char_poly, E_ref, mu1, mu2 - eps,
            tol=continuum_tol, frac=frac, refine=False,
        )
        w_right, _, has_cont_right, _ = amoeba_windings(
            zm, char_poly, E_ref, mu1, mu2 + eps,
            tol=continuum_tol, frac=frac, refine=False,
        )

        if (not has_cont_left) and (not has_cont_right):
            w2_left, w2_right = w_left, w_right
            w2_opposite = bool(
                w2_left * w2_right < 0
            )
            break
    else:
        # All perturbation scales fell inside the continuum band.
        return {
            "w2_left": None, "w2_right": None, "w2_opposite": False,
            "w1_left": None, "w1_right": None, "w1_opposite": None,
            "is_boundary": False, "w1_resolved": False,
        }

    # ---- Step 2: w1 limits via mu1 perturbation ----
    if not w2_opposite:
        return {
            "w2_left": w2_left, "w2_right": w2_right, "w2_opposite": False,
            "w1_left": None, "w1_right": None, "w1_opposite": None,
            "is_boundary": False, "w1_resolved": False,
        }

    # Re-run inner mu2 bisection at mu1 ± ε.  Each perturbed mu1 builds
    # its own ZeroManager.  If the perturbed mu1 also yields continuum
    # (is_continuum=True), the zeros list is empty —
    # _get_average_winding_from_zeros falls back to single-point sampling,
    # which is correct: with no zeros the winding is constant on the full
    # theta2 circle.
    inner_left = _find_mu2_for_w2_zero(
        char_poly, E_ref, mu1 - continuum_perturb, mu2_low, mu2_high,
        N_points=N_points,
        continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
        continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
        max_range_expansions=max_range_expansions,
        range_expand_factor=range_expand_factor, frac=frac,
    )
    inner_right = _find_mu2_for_w2_zero(
        char_poly, E_ref, mu1 + continuum_perturb, mu2_low, mu2_high,
        N_points=N_points,
        continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
        continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
        max_range_expansions=max_range_expansions,
        range_expand_factor=range_expand_factor, frac=frac,
    )

    zeros_left = inner_left.get("zeros") or []
    zeros_right = inner_right.get("zeros") or []

    w1_left, _ = _get_average_winding_from_zeros(
        char_poly, E_ref, mu1 - continuum_perturb, inner_left["mu2"],
        zeros_left, direction=1,
    )
    w1_right, _ = _get_average_winding_from_zeros(
        char_poly, E_ref, mu1 + continuum_perturb, inner_right["mu2"],
        zeros_right, direction=1,
    )

    w1_opposite = bool(w1_left * w1_right < 0)

    return {
        "w2_left": w2_left, "w2_right": w2_right, "w2_opposite": w2_opposite,
        "w1_left": w1_left, "w1_right": w1_right, "w1_opposite": w1_opposite,
        "is_boundary": w2_opposite and w1_opposite,
        "w1_resolved": True,
    }


def bisect_amoeba_ronkin_min(
    char_poly: CharPoly,
    E_ref: complex,
    mu1_low: float = -1,
    mu1_high: float = 1,
    mu2_low: float = -1,
    mu2_high: float = 1,
    N_points: int = 301,
    continuum_tol: float = 1e-6,
    min_continuum_pts: int = 3,
    continuum_perturb: float = 1e-4,
    max_iter: int = 60,
    xtol: float = 1e-10,
    max_range_expansions: int = 10,
    range_expand_factor: float = 2.0,
    frac: float = CONTINUUM_FRAC,
    zm_run_kwargs: Optional[dict] = None,
) -> dict:
    """
    Find the Ronkin function minimum by bisecting mu1 and mu2.

    Outer loop: bisect mu1.
    Inner loop: for each mu1, build an AmoebaZeroManager (μ₂-independent,
    reused across all μ₂ evaluations), then find mu2 where a2 average
    winding = 0, then evaluate w1 average winding at (mu1, mu2).

    The Ronkin minimum satisfies w1 = w2 = 0 simultaneously.

    Continuum handling:
    - When w1 is degenerate at (mu1_mid, mu2_0), perturb mu1 ± epsilon.
      For each perturbed mu1, re-run the inner mu2 bisection to find w2=0,
      then compute w1.
      * Opposite signs → this is the boundary, stop.
      * Same sign → use the sign to continue the outer bisection.

    Returns a dict with keys:
        mu1: critical mu1 value
        mu2: critical mu2 value
        zeros: list of (theta1, theta2, jump) crossing pairs at the critical point
        is_continuum: whether the result is a continuum point
        _zm: the AmoebaZeroManager built at the solved mu1 (reused by
             collect_GBZ_subsets for subset extraction, avoiding a 2nd run)
    """
    # Evaluate w1 at the mu1 endpoints, with adaptive range expansion
    low, high = float(mu1_low), float(mu1_high)
    w1_low = w1_high = 0.0
    inner_low = inner_high = None

    for _ in range(max_range_expansions):
        inner_low = _find_mu2_for_w2_zero(
            char_poly, E_ref, low, mu2_low, mu2_high,
            N_points=N_points,
            continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
            continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
            max_range_expansions=max_range_expansions,
            range_expand_factor=range_expand_factor, frac=frac,
        )
        inner_high = _find_mu2_for_w2_zero(
            char_poly, E_ref, high, mu2_low, mu2_high,
            N_points=N_points,
            continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
            continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
            max_range_expansions=max_range_expansions,
            range_expand_factor=range_expand_factor, frac=frac,
        )

        w1_low, _ = _get_average_winding_from_zeros(
            char_poly, E_ref, low, inner_low["mu2"],
            inner_low["zeros"], direction=1,
        )
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

    for _ in range(max_iter):
        mu1_mid = 0.5 * (mu1_low + mu1_high)

        # Build the ZeroManager once for this mu1_mid — reused across all
        # mu2 evaluations in the inner bisection.
        zm = AmoebaZeroManager(char_poly, E_ref, mu1_mid)
        zm.run(**(zm_run_kwargs or {}))

        inner_mid = _find_mu2_for_w2_zero(
            char_poly, E_ref, mu1_mid, mu2_low, mu2_high,
            N_points=N_points,
            continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
            continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
            max_range_expansions=max_range_expansions,
            range_expand_factor=range_expand_factor, frac=frac, _zm=zm,
        )

        mu2_mid = inner_mid["mu2"]

        # ---- continuum path ----
        # When the inner bisection hits a continuum, zeros are [] and
        # w1 cannot be computed from them.  Compute w1 left/right limits
        # via mu1 perturbation instead.
        w1_area = 0.0  # plateau pre-check area (only meaningful in non-continuum path)
        if inner_mid["is_continuum"]:
            resolved = _resolve_continuum(
                char_poly, E_ref, mu1_mid, mu2_mid, zm,
                mu2_low=mu2_low, mu2_high=mu2_high,
                continuum_perturb=continuum_perturb,
                N_points=N_points, continuum_tol=continuum_tol,
                min_continuum_pts=min_continuum_pts,
                max_iter=max_iter, xtol=xtol,
                max_range_expansions=max_range_expansions,
                range_expand_factor=range_expand_factor, frac=frac,
            )

            if resolved["is_boundary"]:
                # Both w2 and w1 limits straddle zero → Ronkin minimum.
                return {
                    "mu1": mu1_mid, "mu2": mu2_mid, "zeros": [],
                    "is_continuum": True,
                    "_mu1_bracket": (mu1_low, mu1_high),
                    "_w1_bracket": (w1_low, w1_high),
                    "_exit_reason": "continuum_boundary",
                    "_zm": zm,
                }

            if resolved["w1_resolved"]:
                # w2 opposite but w1 not — use w1 sign for bracket update.
                # w1_left is at mu1_mid − ε, between mu1_low and mu1_mid.
                # Compare with w1_low to determine the zero's location.
                w1_proxy = resolved["w1_left"]
                if w1_low * w1_proxy < 0:
                    mu1_high = mu1_mid
                    w1_high = w1_proxy
                else:
                    mu1_low = mu1_mid
                    w1_low = w1_proxy
                continue

            # w2 limits not opposite or w1 not resolved — should not
            # happen here (inner bisection already confirmed a2 opposite
            # via w_left * w_right < 0).  Fall through to normal path
            # as a safety measure.
            zeros_mid = []
            w1_mid = resolved.get("w1_left", 0.0)
        else:
            # ---- normal (non-continuum) path ----
            zeros_mid = inner_mid["zeros"]

            # Compute w1 at (mu1_mid, mu2_mid), reusing zeros from the
            # inner bisection — no extra root tracking needed.
            w1_mid, w1_area = _get_average_winding_from_zeros(
                char_poly, E_ref, mu1_mid, mu2_mid,
                zeros_mid, direction=1,
            )

        if abs(w1_mid) < xtol or (mu1_high - mu1_low) < xtol:
            exit_reason = "w1_zero" if abs(w1_mid) < xtol else "bracket_xtol"
            return {
                "mu1": mu1_mid,
                "mu2": mu2_mid,
                "zeros": zeros_mid,
                "is_continuum": inner_mid["is_continuum"],
                "_mu1_bracket": (mu1_low, mu1_high),
                "_w1_bracket": (w1_low, w1_high),
                "_exit_reason": exit_reason,
                "_w1_area": w1_area,
                "_zm": zm,
            }

        if w1_low * w1_mid < 0:
            mu1_high = mu1_mid
            w1_high = w1_mid
        else:
            mu1_low = mu1_mid
            w1_low = w1_mid

    mu1_mid = 0.5 * (mu1_low + mu1_high)
    raise ValueError(f"bisect_amoeba_ronkin_min: Bisection failed. E_ref{E_ref}")
