'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-05-19
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

from typing import Optional
import numpy as np
import poly_tools as pt
from math import pi
from cmath import exp

from brute_force_SGBZ.winding import PolyDiffContext

from gbz_types import (
    PointSubset, LineSubset, GBZResult, ConnectedSubset,
    generate_probe_steps, find_cyclic_true_intervals,
)

from .bisect import (
    bisect_amoeba_ronkin_min, _find_mu2_for_a2_zero,
)
from .ronkin_winding import (
    _get_average_winding_from_zeros,
)
from .tracks import _compute_root_tracks


# ---- plateau detection ----

def _is_zero_plateau_probe(point: dict, winding_tol: float) -> bool:
    return (
        point["success"]
        and (not point["is_continuum"])
        and point["zero_count"] == 0
        and abs(point["a1"]) <= winding_tol
    )


def _probe_zero_plateau_near_mu1(
    char_poly: pt.CLaurent,
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
    mu1, a1 is zero and there are no a2-crossing zeros.
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
            inner = _find_mu2_for_a2_zero(
                char_poly, E_ref, mu1_probe, mu2_low, mu2_high,
                target_winding=0.0, N_points=N_points,
                continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
                continuum_perturb=continuum_perturb, max_iter=max_iter,
                xtol=xtol, max_range_expansions=max_range_expansions,
                range_expand_factor=range_expand_factor,
            )
            a1 = _get_average_winding_from_zeros(
                char_poly, E_ref, mu1_probe, inner["mu2"],
                inner["zeros"], direction=1,
            )
            point.update({
                "success": True,
                "mu2": inner["mu2"],
                "a1": float(a1),
                "zero_count": len(inner["zeros"]),
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


# ---- main entry point ----

def check_amoeba(
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
    char_poly = pt.CLaurent(3)
    coeffs_ct = pt.CScalarVec(coeffs)
    degs_ct = pt.CLaurentIndexVec(degs.flatten())
    char_poly.set_Laurent_by_terms(coeffs_ct, degs_ct)

    solver_options = dict(options)
    plateau_check = solver_options.pop("plateau_check", True)
    plateau_winding_tol = solver_options.pop("plateau_winding_tol", None)
    plateau_probe_radius = solver_options.pop("plateau_probe_radius", None)

    try:
        amoeba_res = bisect_amoeba_ronkin_min(
            char_poly, E_ref, **solver_options
        )

        # Build subsets from amoeba result
        subsets = []
        if amoeba_res["is_continuum"]:
            # Extract continuum intervals from root tracks
            tracks = _compute_root_tracks(char_poly, E_ref, amoeba_res["mu1"])
            intervals = _extract_continuum_intervals(tracks, amoeba_res["mu2"])
            for t_start, t_end in intervals:
                subsets.append(LineSubset(
                    E=E_ref, mu1=amoeba_res["mu1"],
                    theta1_start=t_start, theta1_end=t_end,
                    _M=tracks["M"], _N=tracks["N"],
                    _poly_diff=PolyDiffContext(char_poly),
                ))
        else:
            for theta1, theta2, jump_direction in amoeba_res["zeros"]:
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
