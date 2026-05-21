'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2025-11-23 00:00:00
Copyright © YourCompanyName All rights reserved
'''

import numpy as np
import poly_tools as pt
from scipy import optimize
from .strip_winding_number import get_strip_winding, get_minor_degrees, get_roots_and_PMGBZ
from .winding import PolyDiffContext


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
        # 1. Determine mu1_left and mu1_right
        # 1.1 mu1_left
        mu1_left = mu1_guess[0]
        mu1_right = None

        while True:
            left_winding, PMGBZ_points = get_strip_winding(self.poly_diff, E_ref, mu1_left, N_points)
            if left_winding < zero_tol:
                break
            else:
                mu1_right = mu1_left
                mu1_left -= 1

        # Check zero
        if left_winding > - zero_tol:
            return mu1_left, PMGBZ_points

        # 1.2 mu1_right
        if mu1_right is None:
            mu1_right = mu1_guess[1]
            while True:
                right_winding, PMGBZ_points = get_strip_winding(self.poly_diff, E_ref, mu1_right, N_points)
                if right_winding > -zero_tol:
                    break
                else:
                    mu1_right += 1

            if right_winding < zero_tol:
                return mu1_right, PMGBZ_points

        # 2. Find zero of the strip winding number between mu1_left and mu1_right
        def strip_winding_fun(mu1: float):
            w = get_strip_winding(self.poly_diff, E_ref, mu1, N_points)[0]
            if isinstance(w, float):
                return w
            else:
                return np.mean(w)
        
        mu1_0 = optimize.brentq(strip_winding_fun, mu1_left, mu1_right)
        _, PMGBZ_points = get_strip_winding(self.poly_diff, E_ref, mu1_0, N_points)

        return mu1_0, PMGBZ_points

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
    debug_mode: bool = False
):
    print("%.2f" % (perc * 100) + r"%")
    coeffs = pt.CScalarVec(coeffs)
    degs = pt.CLaurentIndexVec(degs.flatten())
    char_poly = pt.CLaurent(3)
    char_poly.set_Laurent_by_terms(coeffs, degs)
    solver = SGBZSolver(char_poly)
    try:
        mu1, PMGBZ_points = solver.solve_for_E(E_ref, N_points=301)
        success = True
    except Exception as e:
        if debug_mode:
            raise e
        else:
            print("Error: %s" % str(e))
        success = False
        return {"success": success, "error": str(e)}

    if PMGBZ_points:
        return {"success": success, "is_PMGBZ": True, "mu1": mu1, "PMGBZ_points": PMGBZ_points}
    else:
        return {"success": success, "is_PMGBZ": False, "mu1": mu1}


def convert_results_to_triplet(
    results: list[dict], # results given by check_SGBZ
    E_list: np.ndarray, # E_list used in check_SGBZ
    coeffs, degs
) -> list[tuple[complex, complex, complex]]:
    """
    Convert the results of check_SGBZ to a list of triplets (E, k1, k2) for each PMGBZ point.
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
        if curr_res["success"] and curr_res["is_PMGBZ"]:
            curr_E = E_list[ind]
            curr_mu1 = curr_res["mu1"]
            sols_arr = None
            for item in curr_res["PMGBZ_points"]:
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
