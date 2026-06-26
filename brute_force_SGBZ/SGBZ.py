'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2025-11-23
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

import numpy as np
import poly_tools as pt
from typing import Optional
from scipy import optimize
from .strip_winding_number import get_strip_winding, get_minor_degrees, get_roots_and_PMGBZ
from .winding import PolyDiffContext


def _scalar_winding(winding) -> float:
    if isinstance(winding, tuple):
        return float(np.nanmean(winding))
    return float(winding)


def _is_zero_plateau_probe(point: dict, zero_tol: float) -> bool:
    return (
        point["success"]
        and (not point["is_continuum"])
        and point["PMGBZ_count"] == 0
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

    eps_min = max(10.0 * zero_tol, 1e-12)
    eps_max = max(float(probe_radius), 4.0 * bracket_width, 100.0 * zero_tol, eps_min)

    steps = {eps_min, eps_max}
    if bracket_width > 0:
        steps.update([
            0.25 * bracket_width,
            0.5 * bracket_width,
            bracket_width,
            2.0 * bracket_width,
            4.0 * bracket_width,
        ])
    if probe_radius > 0:
        steps.update([
            0.25 * probe_radius,
            0.5 * probe_radius,
            float(probe_radius),
        ])

    step = eps_min
    while step < eps_max:
        steps.add(step)
        step *= 2.0
    steps = sorted(step for step in steps if step > 0 and step <= eps_max * (1 + 1e-12))

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
            try:
                winding, PMGBZ_points = get_strip_winding(
                    poly_diff, E_ref, mu1_probe, N_points,
                    continuum_perturb=continuum_perturb,
                    zero_tol=zero_tol,
                )
                is_continuum = isinstance(winding, tuple)
                point.update({
                    "success": True,
                    "winding": _scalar_winding(winding),
                    "raw_winding": winding,
                    "PMGBZ_count": len(PMGBZ_points),
                    "is_continuum": is_continuum,
                })
                if _is_zero_plateau_probe(point, zero_tol):
                    found_plateau = True
                elif point["success"] and (not point["is_continuum"]):
                    if side < 0:
                        saw_left_nonplateau = True
                    else:
                        saw_right_nonplateau = True
            except Exception as exc:
                point["error"] = str(exc)
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
        return result["mu1"], result["PMGBZ_points"]

    def solve_for_E_info(
        self,
        E_ref: complex,
        mu1_guess: tuple[float, float] = (-1, 1),
        zero_tol: float = 1e-10,
        N_points: int = 101,
    ) -> dict:
        # 1. Determine mu1_left and mu1_right
        # 1.1 mu1_left
        mu1_left = mu1_guess[0]
        mu1_right = None

        while True:
            left_winding, PMGBZ_points = get_strip_winding(self.poly_diff, E_ref, mu1_left, N_points)
            left_winding_scalar = _scalar_winding(left_winding)
            if left_winding_scalar < zero_tol:
                break
            else:
                mu1_right = mu1_left
                mu1_left -= 1

        # Check zero
        if left_winding_scalar > - zero_tol:
            return {
                "mu1": mu1_left,
                "PMGBZ_points": PMGBZ_points,
                "winding": left_winding,
                "_mu1_bracket": (mu1_left, mu1_left),
                "_winding_bracket": (left_winding_scalar, left_winding_scalar),
                "_exit_reason": "left_endpoint_zero",
            }

        # 1.2 mu1_right
        if mu1_right is None:
            mu1_right = mu1_guess[1]
            while True:
                right_winding, PMGBZ_points = get_strip_winding(self.poly_diff, E_ref, mu1_right, N_points)
                right_winding_scalar = _scalar_winding(right_winding)
                if right_winding_scalar > -zero_tol:
                    break
                else:
                    mu1_right += 1

            if right_winding_scalar < zero_tol:
                return {
                    "mu1": mu1_right,
                    "PMGBZ_points": PMGBZ_points,
                    "winding": right_winding,
                    "_mu1_bracket": (mu1_right, mu1_right),
                    "_winding_bracket": (right_winding_scalar, right_winding_scalar),
                    "_exit_reason": "right_endpoint_zero",
                }
        else:
            right_winding, _ = get_strip_winding(self.poly_diff, E_ref, mu1_right, N_points)
            right_winding_scalar = _scalar_winding(right_winding)

        # 2. Find zero of the strip winding number between mu1_left and mu1_right
        def strip_winding_fun(mu1: float):
            w = get_strip_winding(self.poly_diff, E_ref, mu1, N_points)[0]
            return _scalar_winding(w)
        
        mu1_0 = optimize.brentq(strip_winding_fun, mu1_left, mu1_right)
        winding_0, PMGBZ_points = get_strip_winding(self.poly_diff, E_ref, mu1_0, N_points)

        return {
            "mu1": mu1_0,
            "PMGBZ_points": PMGBZ_points,
            "winding": winding_0,
            "_mu1_bracket": (mu1_left, mu1_right),
            "_winding_bracket": (left_winding_scalar, right_winding_scalar),
            "_exit_reason": "brentq_zero",
        }

class SGBZChecker:
    char_poly: pt.CLaurent
    poly_diff: PolyDiffContext

    def __init__(self, char_poly: pt.CLaurent):
        self.char_poly = char_poly
        self.poly_diff = PolyDiffContext(char_poly)
    
    def check_for_E_and_mu1(self, E: complex, mu1: float, N_points: int = 101, zero_tol: float = 1e-10):
        winding, PMGBZ_points = get_strip_winding(self.poly_diff, E, mu1, N_points)

        if isinstance(winding, tuple):
            if winding[0] < 0 and winding[1] > 0:
                return True
            else:
                return False
        else:
            if abs(winding) < zero_tol and len(PMGBZ_points) > 0:
                return True
            else:
                return False


def check_SGBZ(
    coeffs: np.ndarray,
    degs: np.ndarray,
    E_ref: complex,
    perc: float,
    debug_mode: bool = False,
    **options,
):
    print("%.2f" % (perc * 100) + r"%")
    coeffs = pt.CScalarVec(coeffs)
    degs = pt.CLaurentIndexVec(degs.flatten())
    char_poly = pt.CLaurent(3)
    char_poly.set_Laurent_by_terms(coeffs, degs)
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
        success = True
    except Exception as e:
        if debug_mode:
            raise e
        else:
            print("Error: %s" % str(e))
        success = False
        return {"success": success, "error": str(e)}

    mu1 = sgbz_res["mu1"]
    SGBZ_points = sgbz_res["PMGBZ_points"]
    if SGBZ_points:
        is_SGBZ = True
        plateau_info = {
            "status": "skipped",
            "found": False,
            "reason": "disabled",
            "points": [],
        }
        classification_reason = "nonempty_SGBZ_points"
        if plateau_check:
            plateau_info = _probe_zero_plateau_near_mu1(
                solver.poly_diff, E_ref, mu1, sgbz_res.get("_mu1_bracket"),
                N_points=N_points, zero_tol=zero_tol,
                continuum_perturb=continuum_perturb,
                probe_radius=plateau_probe_radius,
            )
            if plateau_info["found"]:
                is_SGBZ = False
                classification_reason = "nearby_zero_plateau"
            elif plateau_info["status"] == "not_found":
                classification_reason = "nonempty_SGBZ_points_no_plateau"
            else:
                classification_reason = "SGBZ_plateau_check_inconclusive"
        return {
            "success": success,
            "is_SGBZ": is_SGBZ,
            "mu1": mu1,
            "SGBZ_points": SGBZ_points,
            "_plateau_check": plateau_info["status"],
            "_plateau_check_found": plateau_info["found"],
            "_plateau_probe_points": plateau_info["points"],
            "_classification_reason": classification_reason,
            "_mu1_bracket": sgbz_res.get("_mu1_bracket"),
            "_winding_bracket": sgbz_res.get("_winding_bracket"),
            "_exit_reason": sgbz_res.get("_exit_reason"),
            "winding": sgbz_res.get("winding"),
        }
    else:
        return {
            "success": success,
            "is_SGBZ": False,
            "mu1": mu1,
            "_plateau_check": "found",
            "_plateau_check_found": True,
            "_plateau_probe_points": [],
            "_classification_reason": "empty_SGBZ_zero_plateau",
            "_mu1_bracket": sgbz_res.get("_mu1_bracket"),
            "_winding_bracket": sgbz_res.get("_winding_bracket"),
            "_exit_reason": sgbz_res.get("_exit_reason"),
            "winding": sgbz_res.get("winding"),
        }


def convert_results_to_triplet(
    results: list[dict], # results given by check_SGBZ
    E_list: np.ndarray, # E_list used in check_SGBZ
    coeffs, degs
) -> list[tuple[complex, complex, complex]]:
    """
    Convert the results of check_SGBZ to a list of triplets (E, k1, k2) for each SGBZ point.
    """

    char_poly = pt.CLaurent(3)
    char_poly.set_Laurent_by_terms(
        pt.CScalarVec(coeffs),
        pt.CLaurentIndexVec(degs.flatten())
    )
    char_poly_context = PolyDiffContext(char_poly)
    M, N = get_minor_degrees(char_poly_context)

    new_data_list = []
    for ind in range(len(results)):
        curr_res = results[ind]
        if curr_res["success"] and curr_res["is_SGBZ"]:
            curr_E = E_list[ind]
            curr_mu1 = curr_res["mu1"]
            sols_arr = None
            for item in curr_res["SGBZ_points"]:
                if item["is_continuum"]:
                    if sols_arr is None:
                        # calculate zeros
                        _, theta1_arr, sols_arr, info = get_roots_and_PMGBZ(
                            char_poly_context,
                            curr_E,
                            curr_mu1,
                        )
                    exp_theta1_start = np.exp(1j * item["theta1_start"])
                    theta1_diff = np.angle(np.exp(1j * item["theta1_end"]) / exp_theta1_start) % (2 * np.pi)
                    theta1_mask = np.angle(np.exp(1j * theta1_arr) / exp_theta1_start) % (2 * np.pi) <= theta1_diff
                    theta1_selected = theta1_arr[theta1_mask]
                    sols_selected = sols_arr[theta1_mask, :]
                    for theta1_ind in range(len(theta1_selected)):
                        new_data_list.append((
                            curr_E, 
                            theta1_selected[theta1_ind] - 1j * curr_mu1,
                            np.log(sols_selected[theta1_ind, M-1]) / 1j))
                        new_data_list.append((
                            curr_E, 
                            theta1_selected[theta1_ind] - 1j * curr_mu1,
                            np.log(sols_selected[theta1_ind, M]) / 1j))
                else:
                    curr_theta1 = (item["theta1"] + np.pi) % (2 * np.pi) - np.pi
                    beta2_pos, beta2_neg, beta2_zero = item["beta2_sols"]
                    for beta2 in np.concatenate([beta2_pos, beta2_neg, beta2_zero]):
                        new_data_list.append((
                            curr_E, 
                            curr_theta1 - 1j * curr_mu1,
                            np.log(beta2) / 1j))
    return new_data_list
