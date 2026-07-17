'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2025-11-23
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

import numpy as np
import poly_tools as pt
from math import pi
from cmath import exp, log
from typing import Optional

from scipy import optimize

from .strip_winding_number import get_strip_winding
from .winding import PolyDiffContext
from .pmgbz_detector import get_roots_and_PMGBZ

from gbz_types import (
    PointSubset, LineSubset, GBZResult, ConnectedSubset,
    get_minor_degrees, generate_probe_steps,
)


# ---- internal helpers ----

def _scalar_winding(winding) -> float:
    if isinstance(winding, tuple):
        return float(np.nanmean(winding))
    return float(winding)


def _is_zero_plateau_probe(point: dict, zero_tol: float) -> bool:
    return (
        point["success"]
        and (not point["is_continuum"])
        and point["gbz_count"] == 0
        and abs(point["winding"]) <= zero_tol
    )


def _probe_zero_plateau_near_mu1(
    poly_diff: PolyDiffContext,
    E_ref: complex,
    mu1: float,
    mu1_bracket: Optional[tuple[float, float]],
    N_points: int = 301,
    zero_tol: float = 1e-10,
    continuum_perturb: float = 1e-2,
    probe_radius: Optional[float] = None,
) -> dict:
    """Check whether a nonempty-PMGBZ candidate sits next to a zero plateau."""
    if probe_radius is None:
        probe_radius = continuum_perturb

    bracket_width = 0.0
    if mu1_bracket is not None:
        bracket_width = abs(float(mu1_bracket[1]) - float(mu1_bracket[0]))

    steps = generate_probe_steps(bracket_width, probe_radius, zero_tol)

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
            winding, gbz = get_strip_winding(
                poly_diff, E_ref, mu1_probe, N_points,
                continuum_perturb=continuum_perturb,
                zero_tol=zero_tol,
            )
            is_continuum = isinstance(winding, tuple)
            point.update({
                "success": True,
                "winding": _scalar_winding(winding),
                "raw_winding": winding,
                "gbz_count": len(gbz.subsets),
                "is_continuum": is_continuum,
            })
            if _is_zero_plateau_probe(point, zero_tol):
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
                    "zero_tol": zero_tol,
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
        "zero_tol": zero_tol,
        "probe_radius": probe_radius,
        "bracket_width": bracket_width,
        "steps": steps,
        "points": probe_points,
    }


# ---- solver classes ----

class SGBZSolver:
    char_poly: pt.CLaurent
    poly_diff: PolyDiffContext

    def __init__(self, char_poly: pt.CLaurent):
        self.char_poly = char_poly
        self.poly_diff = PolyDiffContext(char_poly)

    def solve_for_E(
        self,
        E_ref: complex,
        mu1_guess: tuple[float, float] = (-1, 1),
        zero_tol: float = 1e-10,
        N_points: int = 101
    ):
        result = self.solve_for_E_info(E_ref, mu1_guess, zero_tol, N_points)
        return result["mu1"], result["gbz"]

    def solve_for_E_info(
        self,
        E_ref: complex,
        mu1_guess: tuple[float, float] = (-1, 1),
        zero_tol: float = 1e-10,
        N_points: int = 101,
    ) -> dict:
        # 1. Determine mu1_left and mu1_right
        mu1_left = mu1_guess[0]
        mu1_right = None

        while True:
            left_winding, left_gbz = get_strip_winding(self.poly_diff, E_ref, mu1_left, N_points)
            left_winding_scalar = _scalar_winding(left_winding)
            if left_winding_scalar < zero_tol:
                break
            else:
                mu1_right = mu1_left
                mu1_left -= 1

        if left_winding_scalar > -zero_tol:
            return {
                "mu1": mu1_left,
                "gbz": left_gbz,
                "winding": left_winding,
                "_mu1_bracket": (mu1_left, mu1_left),
                "_winding_bracket": (left_winding_scalar, left_winding_scalar),
                "_exit_reason": "left_endpoint_zero",
            }

        # 1.2 mu1_right
        if mu1_right is None:
            mu1_right = mu1_guess[1]
            while True:
                right_winding, right_gbz = get_strip_winding(self.poly_diff, E_ref, mu1_right, N_points)
                right_winding_scalar = _scalar_winding(right_winding)
                if right_winding_scalar > -zero_tol:
                    break
                else:
                    mu1_right += 1

            if right_winding_scalar < zero_tol:
                return {
                    "mu1": mu1_right,
                    "gbz": right_gbz,
                    "winding": right_winding,
                    "_mu1_bracket": (mu1_right, mu1_right),
                    "_winding_bracket": (right_winding_scalar, right_winding_scalar),
                    "_exit_reason": "right_endpoint_zero",
                }
        else:
            right_winding, _ = get_strip_winding(self.poly_diff, E_ref, mu1_right, N_points)
            right_winding_scalar = _scalar_winding(right_winding)

        # 2. Find zero of the strip winding number
        def strip_winding_fun(mu1: float):
            w = get_strip_winding(self.poly_diff, E_ref, mu1, N_points)[0]
            return _scalar_winding(w)

        mu1_0 = optimize.brentq(strip_winding_fun, mu1_left, mu1_right)
        winding_0, gbz_0 = get_strip_winding(self.poly_diff, E_ref, mu1_0, N_points)

        return {
            "mu1": mu1_0,
            "gbz": gbz_0,
            "winding": winding_0,
            "_mu1_bracket": (mu1_left, mu1_right),
            "_winding_bracket": (left_winding_scalar, right_winding_scalar),
            "_exit_reason": "brentq_zero",
        }


# ---- main entry point ----

def check_SGBZ(
    coeffs: np.ndarray,
    degs: np.ndarray,
    E_ref: complex,
    perc: float,
    debug_mode: bool = False,
    **options,
) -> GBZResult:
    """Check SGBZ condition and return GBZ points for a reference energy.

    Returns:
        GBZResult with connected subsets.  ``gbz.is_empty`` / ``gbz.index == (0,0)``
        means E_ref is outside the SGBZ spectrum.
    """
    print("%.2f" % (perc * 100) + r"%")
    coeffs_ct = pt.CScalarVec(coeffs)
    degs_ct = pt.CLaurentIndexVec(degs.flatten())
    char_poly = pt.CLaurent(3)
    char_poly.set_Laurent_by_terms(coeffs_ct, degs_ct)
    solver = SGBZSolver(char_poly)

    solver_options = dict(options)
    plateau_check = solver_options.pop("plateau_check", True)
    plateau_probe_radius = solver_options.pop("plateau_probe_radius", None)
    mu1_guess = solver_options.pop("mu1_guess", (-1, 1))
    zero_tol = solver_options.pop("zero_tol", 1e-10)
    N_points = solver_options.pop("N_points", 301)
    continuum_perturb = solver_options.pop("continuum_perturb", 1e-2)

    try:
        sgbz_res = solver.solve_for_E_info(
            E_ref, mu1_guess=mu1_guess, zero_tol=zero_tol, N_points=N_points,
        )
        gbz: GBZResult = sgbz_res["gbz"]
        mu1 = sgbz_res["mu1"]
    except Exception as e:
        if debug_mode:
            raise e
        print("Error: %s" % str(e))
        return GBZResult(E_ref=E_ref, success=False, error=str(e))

    if not gbz.is_empty:
        # Run plateau check
        classification_reason = "nonempty_gbz"
        if plateau_check:
            plateau_info = _probe_zero_plateau_near_mu1(
                solver.poly_diff, E_ref, mu1, sgbz_res.get("_mu1_bracket"),
                N_points=N_points, zero_tol=zero_tol,
                continuum_perturb=continuum_perturb,
                probe_radius=plateau_probe_radius,
            )
            if plateau_info["found"]:
                gbz = GBZResult(E_ref=E_ref, subsets=[], index=(0, 0))
                classification_reason = "nearby_zero_plateau"
            elif plateau_info["status"] == "not_found":
                classification_reason = "nonempty_gbz_no_plateau"
            else:
                classification_reason = "gbz_plateau_check_inconclusive"
    else:
        classification_reason = "empty_gbz_zero_plateau"

    return gbz


# ---- triplet conversion ----

def convert_gbz_to_triplets(
    gbz: GBZResult,
) -> list[tuple[complex, complex, complex]]:
    """Convert a single GBZResult to (E, k1, k2) triplets.

    PointSubset: directly from beta1, beta2.
    LineSubset: calls fill_beta2() (lazy load, cached on first call).
    """
    triplets: list[tuple[complex, complex, complex]] = []
    for subset in gbz.subsets:
        if isinstance(subset, PointSubset):
            k1 = -1j * log(subset.beta1)
            k2 = -1j * log(subset.beta2)
            triplets.append((subset.E, k1, k2))
        elif isinstance(subset, LineSubset):
            subset.fill_beta2()
            if subset.beta2_arr is None:
                continue
            n_pts = len(subset.beta2_arr)
            theta_arr = np.linspace(
                subset.theta1_start, subset.theta1_end, n_pts, endpoint=False,
            )
            for i in range(n_pts):
                theta1 = theta_arr[i]
                k1 = theta1 - 1j * subset.mu1
                for j in range(subset.beta2_arr.shape[1]):
                    k2 = -1j * log(subset.beta2_arr[i, j])
                    triplets.append((subset.E, k1, k2))
    return triplets


def convert_gbz_list_to_triplets(
    gbz_list: list[GBZResult],
) -> list[tuple[complex, complex, complex]]:
    """Convert a list of GBZResult to (E, k1, k2) triplets."""
    all_triplets: list[tuple[complex, complex, complex]] = []
    for gbz in gbz_list:
        all_triplets.extend(convert_gbz_to_triplets(gbz))
    return all_triplets
