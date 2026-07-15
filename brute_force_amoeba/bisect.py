"""
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-14
Copyright © Department of Physics, Tsinghua University. All rights reserved

Bisection algorithms for amoeba Ronkin function critical-point search.

Provides mu2 bisection (a2 winding crossing) and the outer mu1/mu2
bisection for the Ronkin minimum.
"""

from typing import Optional
import numpy as np
import poly_tools as pt


from .ronkin_winding import (
    _compute_winding_from_tracks, _get_average_winding_from_zeros,
)
from .tracks import _compute_root_tracks

def _refine_and_correct(
    char_poly, E_ref, mu1, mu2_0, target,
    tracks, low, high, low_init, high_init, f_low_unref, f_high_unref,
    continuum_tol, min_continuum_pts, xtol,
):
    """Refine crossings at mu2_0, then apply Newton correction if needed.

    Uses the analytical dW/dmu2 computed from zero derivatives to take one
    Newton step, avoiding re-bisection after refinement.
    """
    w_ref, zeros_ref, _, dW_dmu2 = _compute_winding_from_tracks(
        char_poly, E_ref, mu1, mu2_0, tracks,
        continuum_tol, min_continuum_pts, refine_crossings=True,
    )
    f_ref = w_ref - target

    if abs(f_ref) < xtol:
        return {
            "mu2": mu2_0, "zeros": zeros_ref, "is_continuum": False,
            "winding": w_ref, "_tracks": tracks,
        }

    # Newton correction using analytical derivative
    if abs(dW_dmu2) > 1e-15:
        delta = -f_ref / dW_dmu2
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
            "winding": w_ref, "_tracks": tracks,
        }

    # Refine at the corrected mu2
    w_ref2, zeros_ref2, _, _ = _compute_winding_from_tracks(
        char_poly, E_ref, mu1, mu2_new, tracks,
        continuum_tol, min_continuum_pts, refine_crossings=True,
    )
    f_ref2 = w_ref2 - target

    if abs(f_ref2) < xtol:
        return {
            "mu2": mu2_new, "zeros": zeros_ref2, "is_continuum": False,
            "winding": w_ref2, "_tracks": tracks,
        }

    # Fall back: return the better of the two refined points
    if abs(f_ref2) < abs(f_ref):
        return {
            "mu2": mu2_new, "zeros": zeros_ref2, "is_continuum": False,
            "winding": w_ref2, "_tracks": tracks,
        }
    return {
        "mu2": mu2_0, "zeros": zeros_ref, "is_continuum": False,
        "winding": w_ref, "_tracks": tracks,
    }


def _find_mu2_for_a2_zero(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1: float,
    mu2_low: float,
    mu2_high: float,
    target_winding: float = 0.0,
    N_points: int = 301,
    continuum_tol: float = 1e-8,
    min_continuum_pts: int = 3,
    continuum_perturb: float = 1e-4,
    max_iter: int = 60,
    xtol: float = 1e-10,
    max_range_expansions: int = 10,
    range_expand_factor: float = 2.0,
    _root_tracks: Optional[dict] = None,
) -> dict:
    """Find mu2 where a2 winding crosses target_winding, with adaptive range.

    The a2 winding is monotonic in mu2, so expanding the search range
    guarantees a sign change will eventually be found.

    Uses unrefined (cheap) winding during bisection, then refines crossings
    at the final mu2 and applies a Newton correction using the analytical
    derivative dW/dmu2.

    _root_tracks: optionally pre-computed root tracks from _compute_root_tracks.
    When provided, avoids redundant Hungarian matching.
    """
    low, high = float(mu2_low), float(mu2_high)

    # Pre-compute root tracks once (Independent of mu2)
    if _root_tracks is not None:
        tracks = _root_tracks
    else:
        tracks = _compute_root_tracks(char_poly, E_ref, mu1, N_points)

    def _winding_at(mu2_val, refine=False):
        """Evaluate winding at mu2_val.  Strips dW_dmu2 for the bisection loop."""
        w, z, c, _ = _compute_winding_from_tracks(
            char_poly, E_ref, mu1, mu2_val, tracks,
            continuum_tol, min_continuum_pts,
            refine_crossings=refine,
        )
        return w, z, c

    # Adaptive range expansion (unrefined)
    for _ in range(max_range_expansions):
        w_low, _, _ = _winding_at(low)
        w_high, _, _ = _winding_at(high)
        f_low = w_low - target_winding
        f_high = w_high - target_winding

        if f_low * f_high <= 0:
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
            f_left = w_left - target_winding
            f_right = w_right - target_winding

            if f_left * f_right < 0:
                # Refine zeros at the continuum boundary
                w_ref, zeros_ref, _, _ = _compute_winding_from_tracks(
                    char_poly, E_ref, mu1, mu2_mid, tracks,
                    continuum_tol, min_continuum_pts, refine_crossings=True,
                )
                return {
                    "mu2": mu2_mid, "zeros": zeros_ref, "is_continuum": True,
                    "winding": (w_left, w_right), "_tracks": tracks,
                }
            else:
                if f_low * f_left > 0:
                    low = mu2_mid
                    f_low = f_left
                else:
                    high = mu2_mid
                    f_high = f_left
                continue

        f_mid = w_mid - target_winding

        if abs(f_mid) < xtol or (high - low) < xtol:
            # --- Post-refinement + Newton correction ---
            return _refine_and_correct(
                char_poly, E_ref, mu1, mu2_mid, target_winding,
                tracks, low, high, low_init, high_init, f_low, f_high,
                continuum_tol, min_continuum_pts, xtol,
            )

        if f_low * f_mid < 0:
            high = mu2_mid
            f_high = f_mid
        else:
            low = mu2_mid
            f_low = f_mid

    # Max iterations reached — refine at the final midpoint
    mu2_mid = 0.5 * (low + high)
    w_ref, zeros_ref, _, _ = _compute_winding_from_tracks(
        char_poly, E_ref, mu1, mu2_mid, tracks,
        continuum_tol, min_continuum_pts, refine_crossings=True,
    )
    return {
        "mu2": mu2_mid, "zeros": zeros_ref, "is_continuum": False,
        "winding": w_ref, "_tracks": tracks,
    }


def _make_ronkin_result(
    mu1: float,
    mu2: float,
    zeros: list[tuple[float, float]],
    is_continuum: bool,
    mu1_low: float,
    mu1_high: float,
    a1_low: float,
    a1_high: float,
    exit_reason: str,
) -> dict:
    return {
        "mu1": mu1,
        "mu2": mu2,
        "zeros": zeros,
        "is_continuum": is_continuum,
        "_mu1_bracket": (mu1_low, mu1_high),
        "_a1_bracket": (a1_low, a1_high),
        "_exit_reason": exit_reason,
    }


def bisect_amoeba_ronkin_min(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1_low: float = -1,
    mu1_high: float = 1,
    mu2_low: float = -1,
    mu2_high: float = 1,
    N_points: int = 301,
    continuum_tol: float = 1e-8,
    min_continuum_pts: int = 3,
    continuum_perturb: float = 1e-4,
    max_iter: int = 60,
    xtol: float = 1e-10,
    max_range_expansions: int = 10,
    range_expand_factor: float = 2.0,
) -> dict:
    """
    Find the Ronkin function minimum by bisecting mu1 and mu2.

    Outer loop: bisect mu1.
    Inner loop: for each mu1, find mu2 where a2 average winding = 0,
    then evaluate a1 average winding at (mu1, mu2).

    The Ronkin minimum satisfies a1 = a2 = 0 simultaneously.

    Continuum handling:
    - When a1 is degenerate at (mu1_mid, mu2_0), perturb mu1 ± epsilon.
      For each perturbed mu1, re-run the inner mu2 bisection to find a2=0,
      then compute a1.
      * Opposite signs → this is the boundary, stop.
      * Same sign → use the sign to continue the outer bisection.

    Returns a dict with keys:
        mu1: critical mu1 value
        mu2: critical mu2 value
        zeros: list of (theta1, theta2) crossing pairs at the critical point
        is_continuum: whether the result is a continuum point
    """
    # Evaluate a1 at the mu1 endpoints, with adaptive range expansion
    low, high = float(mu1_low), float(mu1_high)
    a1_low = a1_high = 0.0
    inner_low = inner_high = None

    for _ in range(max_range_expansions):
        inner_low = _find_mu2_for_a2_zero(
            char_poly, E_ref, low, mu2_low, mu2_high,
            target_winding=0.0, N_points=N_points,
            continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
            continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
            max_range_expansions=max_range_expansions,
            range_expand_factor=range_expand_factor,
        )
        inner_high = _find_mu2_for_a2_zero(
            char_poly, E_ref, high, mu2_low, mu2_high,
            target_winding=0.0, N_points=N_points,
            continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
            continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
            max_range_expansions=max_range_expansions,
            range_expand_factor=range_expand_factor,
        )

        a1_low = _get_average_winding_from_zeros(
            char_poly, E_ref, low, inner_low["mu2"],
            inner_low["zeros"], direction=1,
        )
        a1_high = _get_average_winding_from_zeros(
            char_poly, E_ref, high, inner_high["mu2"],
            inner_high["zeros"], direction=1,
        )

        if a1_low * a1_high <= 0:
            break

        width = high - low
        low = low - range_expand_factor * width
        high = high + range_expand_factor * width
    else:
        raise ValueError(f"bisect_amoeba_ronkin_min: reach maximum expansion. E_ref:{E_ref}, low:{low}, high:{high}")

    mu1_low, mu1_high = low, high

    for _ in range(max_iter):
        mu1_mid = 0.5 * (mu1_low + mu1_high)

        inner_mid = _find_mu2_for_a2_zero(
            char_poly, E_ref, mu1_mid, mu2_low, mu2_high,
            target_winding=0.0, N_points=N_points,
            continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
            continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
            max_range_expansions=max_range_expansions,
            range_expand_factor=range_expand_factor,
        )

        mu2_mid = inner_mid["mu2"]
        zeros_mid = inner_mid["zeros"]

        # Normal case: compute a1 at (mu1_mid, mu2_mid)
        # Reuse zeros from inner bisection — no extra Hungarian matching
        a1_mid = _get_average_winding_from_zeros(
            char_poly, E_ref, mu1_mid, mu2_mid,
            zeros_mid, direction=1,
        )

        if abs(a1_mid) < xtol or (mu1_high - mu1_low) < xtol:
            exit_reason = "a1_zero" if abs(a1_mid) < xtol else "bracket_xtol"
            return _make_ronkin_result(
                mu1_mid, mu2_mid, zeros_mid, inner_mid["is_continuum"],
                mu1_low, mu1_high, a1_low, a1_high,
                exit_reason,
            )

        if a1_low * a1_mid < 0:
            mu1_high = mu1_mid
            a1_high = a1_mid
        else:
            mu1_low = mu1_mid
            a1_low = a1_mid

    mu1_mid = 0.5 * (mu1_low + mu1_high)
    inner_mid = _find_mu2_for_a2_zero(
        char_poly, E_ref, mu1_mid, mu2_low, mu2_high,
        target_winding=0.0, N_points=N_points,
        continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
        continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
        max_range_expansions=max_range_expansions,
        range_expand_factor=range_expand_factor,
    )
    raise ValueError(f"bisect_amoeba_ronkin_min: Bisection failed. E_ref{E_ref}")


