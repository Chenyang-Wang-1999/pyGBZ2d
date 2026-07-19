'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-05-19
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

from typing import Optional
import numpy as np
from math import pi, sqrt
from cmath import exp

from gbz_types import (
    PointSubset, LineSubset, GBZResult, CharPoly,
    generate_probe_steps, find_cyclic_true_intervals,
)

from .bisect import (
    bisect_amoeba_ronkin_min, _find_mu2_for_w2_zero,
)
from .ronkin_winding import (
    _get_average_winding_from_zeros,
    _find_exact_crossing,
)
from .tracks import _compute_root_tracks


# ---- plateau pre-check helpers ----

def _check_zeros_are_clustered(
    zeros: list[tuple[float, float, int]],
    tol_normalized: float,
) -> bool:
    """Check whether every zero has another zero within tol_normalized distance.

    Each zero is a 3-tuple (θ₁, θ₂, jump);
    only the first two components are used.

    Distance is the Euclidean metric on the (θ₁, θ₂)-torus [0, 2π)²,
    normalized by 2π so that the full torus diagonal is √2.

    At a genuine GBZ point zeros are well-separated (they partition the
    circle into meaningful segments).  At a zero-plateau boundary the
    winding changes sign over a vanishingly narrow angular region, so
    the zeros cluster into nearly degenerate pairs — each zero sits
    within tol_normalized of a neighbour.

    Returns False when any zero is isolated (no neighbour within the
    threshold), which rules out a plateau and allows the caller to skip
    expensive probing.
    """
    if len(zeros) < 2:
        return False
    # Convert normalized tolerance back to radians
    tol_rad = tol_normalized * (2 * pi)
    for i in range(len(zeros)):
        t1_i, t2_i = zeros[i][0], zeros[i][1]
        has_neighbor = False
        for j in range(len(zeros)):
            if i == j:
                continue
            t1_j, t2_j = zeros[j][0], zeros[j][1]
            # Torus distance in each angular dimension
            d1 = abs(t1_i - t1_j)
            d1 = min(d1, 2 * pi - d1)
            d2 = abs(t2_i - t2_j)
            d2 = min(d2, 2 * pi - d2)
            if sqrt(d1 * d1 + d2 * d2) < tol_rad:
                has_neighbor = True
                break
        if not has_neighbor:
            return False
    return True


# ---- plateau detection ----

def _is_zero_plateau_probe(point: dict, winding_tol: float) -> bool:
    return (
        point["success"]
        and (not point["is_continuum"])
        and point["zero_count"] == 0
        and abs(point["w1"]) <= winding_tol
    )


def _probe_zero_plateau_near_mu1(
    char_poly: CharPoly,
    E_ref: complex,
    mu1: float,
    mu1_bracket: Optional[tuple[float, float]],
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
    winding_tol: Optional[float] = None,
    probe_radius: Optional[float] = None,
) -> dict:
    """Check whether a nonempty-zero candidate sits next to a zero plateau.

    The decisive plateau signature is strict: after solving a2=0 at a nearby
    mu1, w1 is zero and there are no a2-crossing zeros.
    """
    if winding_tol is None:
        winding_tol = max(xtol, 1e-10)
    if probe_radius is None:
        probe_radius = continuum_perturb

    bracket_width = 0.0
    if mu1_bracket is not None:
        bracket_width = abs(float(mu1_bracket[1]) - float(mu1_bracket[0]))

    steps = generate_probe_steps(bracket_width, probe_radius, winding_tol)

    probe_points = []
    found_plateau = False
    saw_left_nonplateau = False
    saw_right_nonplateau = False

    for step in steps:
        for side in (-1, 1):
            mu1_probe = mu1 + side * step
            point = {
                "mu1": mu1_probe,
                "side": side,
                "step": step,
                "success": False,
            }
            inner = _find_mu2_for_w2_zero(
                char_poly, E_ref, mu1_probe, mu2_low, mu2_high, N_points=N_points,
                continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
                continuum_perturb=continuum_perturb, max_iter=max_iter,
                xtol=xtol, max_range_expansions=max_range_expansions,
                range_expand_factor=range_expand_factor,
            )
            zeros_probe = inner.get("zeros") or []
            w1, _ = _get_average_winding_from_zeros(
                char_poly, E_ref, mu1_probe, inner["mu2"],
                zeros_probe, direction=1,
            )
            point.update({
                "success": True,
                "mu2": inner["mu2"],
                "w1": float(w1),
                "zero_count": len(zeros_probe),
                "is_continuum": bool(inner["is_continuum"]),
            })
            if _is_zero_plateau_probe(point, winding_tol):
                found_plateau = True
            elif point["success"] and (not point["is_continuum"]):
                if side < 0:
                    saw_left_nonplateau = True
                else:
                    saw_right_nonplateau = True
            probe_points.append(point)
            if found_plateau:
                return {
                    "status": "found",
                    "found": True,
                    "winding_tol": winding_tol,
                    "probe_radius": probe_radius,
                    "bracket_width": bracket_width,
                    "steps": steps,
                    "points": probe_points,
                }

    if found_plateau:
        status = "found"
    elif saw_left_nonplateau and saw_right_nonplateau:
        status = "not_found"
    else:
        status = "inconclusive"

    return {
        "status": status,
        "found": found_plateau,
        "winding_tol": winding_tol,
        "probe_radius": probe_radius,
        "bracket_width": bracket_width,
        "steps": steps,
        "points": probe_points,
    }


# ---- continuum extraction ----

def _extract_continuum_intervals(
    tracks: dict,
    mu2: float,
    continuum_tol: float = 1e-8,
    min_continuum_pts: int = 3,
) -> list[tuple[float, float]]:
    """Extract continuum intervals where |ln|beta2| - mu2| < continuum_tol.

    Reuses the continuum-detection logic from _compute_winding_from_tracks.
    """
    theta1_arr = tracks["theta1_arr"]
    tracked = tracks["tracked"]
    n_roots = tracked.shape[1]
    n_pts = len(theta1_arr)

    near_boundary = np.abs(np.log(np.abs(tracked)) - mu2) < continuum_tol

    intervals = []
    for j in range(n_roots):
        for start_idx, end_idx in find_cyclic_true_intervals(near_boundary[:n_pts, j]):
            width = (end_idx - start_idx) % n_pts + 1
            if width >= min_continuum_pts:
                t_start = float(theta1_arr[start_idx])
                # end_idx is inclusive, convert to exclusive for LineSubset
                t_end = float(theta1_arr[(end_idx + 1) % n_pts])
                if t_end <= t_start:
                    t_end += 2 * pi
                intervals.append((t_start, t_end % (2 * pi)))
    return intervals


def _detect_crossings_outside_continuum(
    tracks: dict,
    char_poly: CharPoly,
    E: complex,
    mu1: float,
    mu2: float,
    continuum_tol: float = 1e-8,
    min_continuum_pts: int = 3,
) -> list[tuple[float, float, int]]:
    """Detect refined (theta1, theta2, jump) crossings in non-continuum regions.

    When a continuum is present, not all root tracks stay within the
    continuum band — individual tracks may cross mu2 at isolated theta1
    points in the gaps between continuum intervals.  This function detects
    and refines those discrete crossings, using the same vectorised
    crossover logic as ``_compute_winding_from_tracks`` but skipping
    crossings that fall inside any continuum interval.

    Returns:
        List of (theta1, theta2, jump) tuples for discrete crossings.
    """
    theta1_arr = tracks["theta1_arr"]
    tracked = tracks["tracked"]
    theta1_ext = tracks["theta1_ext"]
    tracked_ext = tracks["tracked_ext"]
    poly = tracks["char_poly"]

    n_roots = tracked.shape[1]
    n_pts = len(theta1_arr)

    # Use pre-computed ln|beta2| when available
    log_abs_all = tracks.get("log_abs_tracked_ext")
    if log_abs_all is None:
        log_abs_all = np.log(np.abs(tracked_ext))

    near_boundary = np.abs(log_abs_all - mu2) < continuum_tol  # (n_pts+1, n_roots)

    refined_crossings = []

    for j in range(n_roots):
        # ---- find continuum intervals for this track ----
        continuum_intervals_j = []
        for start_idx, end_idx in find_cyclic_true_intervals(
            near_boundary[:n_pts, j]
        ):
            width = (end_idx - start_idx) % n_pts + 1
            if width >= min_continuum_pts:
                t_start = float(theta1_arr[start_idx])
                t_end = float(theta1_arr[(end_idx + 1) % n_pts])
                if t_end <= t_start:
                    t_end += 2 * pi
                t_end = t_end % (2 * pi)
                continuum_intervals_j.append((t_start, t_end))

        # Helper: is a theta1 value inside any continuum interval of this track?
        def _in_continuum(t1: float) -> bool:
            for t_s, t_e in continuum_intervals_j:
                # Normalise t1 to [0, 2π)
                t1_norm = t1 % (2 * pi)
                if t_e > t_s:
                    if t_s <= t1_norm <= t_e:
                        return True
                else:
                    # Wrap-around interval
                    if t1_norm >= t_s or t1_norm <= t_e:
                        return True
            return False

        # ---- detect crossings for this track in non-continuum gaps ----
        d = log_abs_all[:, j] - mu2  # (n_pts+1,)
        sign_change = (d[:-1] * d[1:]) < 0
        nb = near_boundary[:, j]
        noise = nb[:-1] & nb[1:]  # both ends in continuum band
        valid = sign_change & ~noise

        for i in np.where(valid)[0]:
            # Linear interpolation
            frac = (mu2 - log_abs_all[i, j]) / (
                log_abs_all[i + 1, j] - log_abs_all[i, j]
            )
            theta1_approx = theta1_ext[i] + frac * (
                theta1_ext[i + 1] - theta1_ext[i]
            )
            beta2_approx = tracked_ext[i, j] + frac * (
                tracked_ext[i + 1, j] - tracked_ext[i, j]
            )

            # Skip crossings inside continuum intervals
            if _in_continuum(theta1_approx):
                continue

            jump = 1 if log_abs_all[i, j] < mu2 else -1
            theta2_guess = float(np.angle(beta2_approx))

            # Refine
            result = _find_exact_crossing(
                poly, E, mu1, mu2,
                theta1_approx % (2 * pi), theta2_guess,
            )
            if result is not None:
                refined_crossings.append((result[0], result[1], jump))
            else:
                refined_crossings.append(
                    (theta1_approx % (2 * pi), theta2_guess % (2 * pi), jump)
                )

    return refined_crossings


# ---- main entry point ----

def collect_GBZ_subsets(
    coeffs: np.ndarray,
    degs: np.ndarray,
    E_ref: complex,
    perc: float,
    debug_mode: bool = False,
    **options,
) -> GBZResult:
    """Check amoeba condition and return GBZ points for a reference energy.

    Returns:
        GBZResult with connected subsets.  ``gbz.is_empty`` means E_ref is
        outside the amoeba GBZ spectrum.
    """
    print("%.2f" % (perc * 100) + r"%")
    char_poly = CharPoly(coeffs, degs)

    solver_options = dict(options)
    plateau_check = solver_options.pop("plateau_check", True)
    plateau_winding_tol = solver_options.pop("plateau_winding_tol", None)
    plateau_probe_radius = solver_options.pop("plateau_probe_radius", None)
    plateau_area_threshold = solver_options.pop("plateau_area_threshold", 1e-2)

    try:
        amoeba_res = bisect_amoeba_ronkin_min(
            char_poly, E_ref, **solver_options
        )

        # Build subsets from amoeba result
        subsets = []
        if amoeba_res["is_continuum"]:
            # Extract continuum intervals from root tracks.
            # Reuse _tracks from the bisection when available (avoids a 3rd
            # Hungarian matching for the same (E, mu1)), otherwise compute fresh.
            tracks = amoeba_res.get("_tracks")
            if tracks is None:
                tracks = _compute_root_tracks(char_poly, E_ref, amoeba_res["mu1"])
            intervals = _extract_continuum_intervals(tracks, amoeba_res["mu2"])
            poly = tracks.get("char_poly", char_poly)

            # 1. Collect continuum intervals → LineSubset
            for t_start, t_end in intervals:
                subsets.append(LineSubset(
                    E=E_ref, mu1=amoeba_res["mu1"],
                    theta1_start=t_start, theta1_end=t_end,
                ))

            # 2. Collect discrete crossings in non-continuum gaps → PointSubset
            crossings = _detect_crossings_outside_continuum(
                tracks, char_poly, E_ref, amoeba_res["mu1"], amoeba_res["mu2"],
            )
            for theta1, theta2, jump_direction in crossings:
                subsets.append(PointSubset(
                    E=E_ref,
                    beta1=exp(amoeba_res["mu1"] + 1j * theta1),
                    beta2=exp(amoeba_res["mu2"] + 1j * theta2),
                ))
        else:
            for theta1, theta2, jump_direction in (amoeba_res.get("zeros") or []):
                subsets.append(PointSubset(
                    E=E_ref,
                    beta1=exp(amoeba_res["mu1"] + 1j * theta1),
                    beta2=exp(amoeba_res["mu2"] + 1j * theta2),
                ))

        # Determine if this is actually a zero plateau (no GBZ)
        is_amoeba = True
        if (not amoeba_res["is_continuum"]) and (len(amoeba_res["zeros"]) == 0):
            is_amoeba = False
        elif plateau_check and (not amoeba_res["is_continuum"]):
            # Lightweight plateau pre-check: at a zero-plateau boundary
            # three necessary conditions hold simultaneously:
            #
            #   (a) w1 & w2 non-zero winding areas are tiny — most of
            #       the angular circle has u = 0.
            #   (b) zeros cluster into nearly degenerate pairs — at a
            #       genuine GBZ point zeros are well-separated, while a
            #       plateau boundary concentrates all sign changes in a
            #       vanishingly narrow angular region.
            #
            # Skip the expensive plateau probe when any condition fails.
            _should_probe = False
            w1_area = amoeba_res.get("_w1_area")
            if w1_area is None:
                _, w1_area = _get_average_winding_from_zeros(
                    char_poly, E_ref, amoeba_res["mu1"], amoeba_res["mu2"],
                    amoeba_res["zeros"], direction=1,
                )
            if w1_area < plateau_area_threshold:
                _, w2_area = _get_average_winding_from_zeros(
                    char_poly, E_ref, amoeba_res["mu1"], amoeba_res["mu2"],
                    amoeba_res["zeros"], direction=2,
                )
                if w2_area < plateau_area_threshold:
                    _should_probe = _check_zeros_are_clustered(
                        amoeba_res["zeros"], plateau_area_threshold,
                    )

            if _should_probe:
                plateau_info = _probe_zero_plateau_near_mu1(
                    char_poly, E_ref, amoeba_res["mu1"],
                    amoeba_res.get("_mu1_bracket"),
                    mu2_low=solver_options.get("mu2_low", -1),
                    mu2_high=solver_options.get("mu2_high", 1),
                    N_points=solver_options.get("N_points", 301),
                    continuum_tol=solver_options.get("continuum_tol", 1e-8),
                    min_continuum_pts=solver_options.get("min_continuum_pts", 3),
                    continuum_perturb=solver_options.get("continuum_perturb", 1e-4),
                    max_iter=solver_options.get("max_iter", 60),
                    xtol=solver_options.get("xtol", 1e-10),
                    max_range_expansions=solver_options.get("max_range_expansions", 10),
                    range_expand_factor=solver_options.get("range_expand_factor", 2.0),
                    winding_tol=plateau_winding_tol,
                    probe_radius=plateau_probe_radius,
                )
                if plateau_info["found"]:
                    is_amoeba = False

        if not is_amoeba:
            subsets = []

        n_0d = sum(1 for s in subsets if isinstance(s, PointSubset))
        n_1d = sum(1 for s in subsets if isinstance(s, LineSubset))
        return GBZResult(E_ref=E_ref, subsets=subsets, index=(n_0d, n_1d))

    except Exception as e:
        if debug_mode:
            raise e
        return GBZResult(E_ref=E_ref, success=False, error=str(e))
