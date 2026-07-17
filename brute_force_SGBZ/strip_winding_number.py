'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2025-11-18
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

''' Calculations of the strip winding number '''

import numpy as np
import poly_tools as pt
from math import pi
from cmath import exp
from scipy import interpolate
from scipy import optimize

from .winding import WindingFun, PolyDiffContext, get_winding_number
from .pmgbz_detector import get_roots_and_PMGBZ
from gbz_types import GBZResult


def get_loop_winding(
    poly_diff: PolyDiffContext,
    E_ref: complex,
    mu1: float,
    mu2_fun: callable,
    theta2: float,
    N_seg: int = 5
) -> float:
    def param_fun(theta1: float):
        mu2_val = mu2_fun(theta1, extrapolate="periodic")
        mu2_prime = mu2_fun(theta1, extrapolate="periodic", nu=1)
        param = np.array([
            E_ref,
            exp(mu1 + 1j * theta1),
            exp(mu2_val + 1j * theta2)
        ])
        dparam = np.array([
            0,
            1j * exp(mu1 + 1j * theta1),
            exp(mu2_val + 1j * theta2) * mu2_prime
        ])
        return param, dparam

    winding_fun = WindingFun(poly_diff.char_poly, param_fun, (0, 2 * pi))
    return get_winding_number(winding_fun, N_seg=N_seg)


def _strip_winding_from_result(
    poly_diff: PolyDiffContext,
    E_ref: complex,
    mu1: float,
    gbz_result: GBZResult,
    theta1_arr: np.ndarray,
    sols_arr: np.ndarray,
    info: dict,
    N_points: int = 301,
    with_gap_info: bool = False,
    continuum_perturb: float = 1e-2,
    zero_tol: float = 1e-10,
    GBZ_check_tol: float = 1e-6,
):
    """Compute strip winding from pre-computed root results.

    Accepts pre-computed (gbz_result, theta1_arr, sols_arr, info) so
    that continuum recursion can pass through already-computed results,
    avoiding redundant root solving (optimization C).
    """
    # Sort sols_arr by norm at each row
    for row_ind in range(sols_arr.shape[0]):
        sols_arr[row_ind, :] = sols_arr[row_ind, np.argsort(np.abs(sols_arr[row_ind, :]))]

    M = info["M"]
    PMGBZ_raw = info.get("_pmgbz_raw", [])

    if info["continuum_flag"]:
        # --- Continuum case: perturb mu1 to find winding limits ---
        gbz_left, t1_l, sols_l, info_l = get_roots_and_PMGBZ(
            poly_diff, E_ref, mu1 - continuum_perturb, N_points, zero_tol, GBZ_check_tol,
        )
        if info_l["continuum_flag"]:
            W_left = np.nan
        else:
            W_left = _strip_winding_from_result(
                poly_diff, E_ref, mu1 - continuum_perturb,
                gbz_left, t1_l, sols_l, info_l,
                N_points, with_gap_info, continuum_perturb, zero_tol, GBZ_check_tol,
            )[0]

        gbz_right, t1_r, sols_r, info_r = get_roots_and_PMGBZ(
            poly_diff, E_ref, mu1 + continuum_perturb, N_points, zero_tol, GBZ_check_tol,
        )
        if info_r["continuum_flag"]:
            W_right = np.nan
        else:
            W_right = _strip_winding_from_result(
                poly_diff, E_ref, mu1 + continuum_perturb,
                gbz_right, t1_r, sols_r, info_r,
                N_points, with_gap_info, continuum_perturb, zero_tol, GBZ_check_tol,
            )[0]

        if with_gap_info:
            return (W_left, W_right), gbz_result, None
        else:
            return (W_left, W_right), gbz_result

    # --- Normal winding computation ---
    mu2_max_arr = np.log(np.abs(sols_arr[:, M]))
    mu2_min_arr = np.log(np.abs(sols_arr[:, M - 1]))
    mu2_max_arr[-1] = mu2_max_arr[0]
    mu2_min_arr[-1] = mu2_min_arr[0]
    mu2_arr = (mu2_max_arr + mu2_min_arr) / 2
    mu2_fun = interpolate.CubicSpline(theta1_arr, mu2_arr, bc_type="periodic")

    def calc_loop_root_distance(theta2_target: float) -> float:
        """Min Euclidean distance from beta2_loop to ALL roots across theta1.

        Larger = safer theta2.  Replaces the old calculate_min_mu2_diff
        which only checked boundary-root theta2 crossings and returned inf
        for avoided-crossing regions where the loop passes between closely-
        spaced roots without crossing them.
        """
        n_pts = len(theta1_arr) - 1
        beta2_loop = np.exp(mu2_arr[:n_pts] + 1j * theta2_target)
        dist = np.abs(beta2_loop[:, None] - sols_arr[:n_pts, :])
        return float(np.min(dist))

    # Select best theta2
    if not PMGBZ_raw:
        N_samples = 4
        offset = np.random.uniform(-pi, pi)
        theta2_samples = (np.linspace(0, 2 * pi, N_samples, endpoint=False) + offset) % (2 * pi)
        mu2_diff_samples = [calc_loop_root_distance(theta2) for theta2 in theta2_samples]
        theta2_median = theta2_samples[np.argmax(mu2_diff_samples)]
        W_strip = np.round(get_loop_winding(poly_diff, E_ref, mu1, mu2_fun, theta2_median))
    else:
        # Collect beta2 values from raw PMGBZ classification (pos / neg / zero)
        PMGBZ_beta2_pos = []
        PMGBZ_beta2_neg = []
        PMGBZ_beta2_zero = []
        for v in PMGBZ_raw:
            pos_sols, neg_sols, zero_sols = v["beta2_sols"]
            PMGBZ_beta2_pos.extend(list(pos_sols))
            PMGBZ_beta2_neg.extend(list(neg_sols))
            PMGBZ_beta2_zero.extend(list(zero_sols))

        PMGBZ_beta2_all = PMGBZ_beta2_pos + PMGBZ_beta2_neg + PMGBZ_beta2_zero
        PMGBZ_beta2 = np.asarray(PMGBZ_beta2_all, dtype=complex)
        PMGBZ_beta2_pos = np.asarray(PMGBZ_beta2_pos, dtype=complex)
        PMGBZ_beta2_neg = np.asarray(PMGBZ_beta2_neg, dtype=complex)

        theta2_arr = np.angle(PMGBZ_beta2)
        theta2_arr = (theta2_arr + 2 * pi) % (2 * pi)
        theta2_arr.sort()
        if theta2_arr.size == 1:
            theta2_samples = np.array([(theta2_arr[0] + pi) % (2 * pi)])
        else:
            theta2_samples = (theta2_arr[1:] + theta2_arr[:-1]) / 2
            theta2_samples = np.hstack(
                (theta2_samples, [(theta2_arr[0] + 2 * pi + theta2_arr[-1]) / 2])
            )
        mu2_diff_samples = [calc_loop_root_distance(theta2) for theta2 in theta2_samples]
        theta2_median = theta2_samples[np.argmax(mu2_diff_samples)]

        w0 = np.round(get_loop_winding(poly_diff, E_ref, mu1, mu2_fun, theta2_median))
        theta2_pos = np.angle(PMGBZ_beta2_pos / exp(1j * theta2_median)) % (2 * pi)
        theta2_neg = np.angle(PMGBZ_beta2_neg / exp(1j * theta2_median)) % (2 * pi)
        W_strip = w0 + (np.sum(theta2_neg) - np.sum(theta2_pos)) / (2 * pi)

    if with_gap_info:
        if PMGBZ_raw:
            gap_info = None
        else:
            mu2_max_ind = np.argmin(mu2_max_arr)
            mu2_min_ind = np.argmax(mu2_min_arr)

            def beta2_solving_fun(beta2, theta1):
                var_ctype = pt.CScalarVec([E_ref, exp(mu1 + 1j * theta1), beta2])
                return poly_diff.char_poly.eval(var_ctype), poly_diff.dchar_poly[2].eval(var_ctype)

            def derivative_fun(theta1, beta2_initial):
                res = optimize.root_scalar(
                    beta2_solving_fun,
                    args=(theta1,),
                    x0=beta2_initial,
                    fprime=True
                )
                beta2_val = res.root
                return poly_diff.eval_dmu2([E_ref, exp(mu1 + 1j * theta1), beta2_val])[1]

            theta1_extended = np.hstack(([theta1_arr[-2] - 2 * pi], theta1_arr, [theta1_arr[1] + 2 * pi]))
            theta1_max = optimize.brentq(derivative_fun, theta1_extended[mu2_max_ind], theta1_extended[mu2_max_ind + 2],
                    args=(sols_arr[mu2_max_ind, M],))
            theta1_min = optimize.brentq(derivative_fun, theta1_extended[mu2_min_ind], theta1_extended[mu2_min_ind + 2],
                    args=(sols_arr[mu2_min_ind, M - 1],))

            beta2_max = optimize.root_scalar(
                beta2_solving_fun, args=(theta1_max,),
                x0=sols_arr[mu2_max_ind, M], fprime=True
            ).root
            beta2_min = optimize.root_scalar(
                beta2_solving_fun, args=(theta1_min,),
                x0=sols_arr[mu2_min_ind, M - 1], fprime=True
            ).root

            gap_info = {
                "point_max": (E_ref, exp(mu1 + 1j * theta1_max), beta2_max),
                "point_min": (E_ref, exp(mu1 + 1j * theta1_min), beta2_min),
                "mu2_diff": mu2_max_arr[mu2_max_ind] - mu2_min_arr[mu2_min_ind],
                "dmu2_diff": poly_diff.eval_dmu2([E_ref, exp(mu1 + 1j * theta1_max), beta2_max])[0]
                           - poly_diff.eval_dmu2([E_ref, exp(mu1 + 1j * theta1_min), beta2_min])[0]
            }

        return W_strip, gbz_result, gap_info
    else:
        return W_strip, gbz_result


def get_strip_winding(
    poly_diff: PolyDiffContext,
    E_ref: complex,
    mu1: float,
    N_points: int = 301,
    with_gap_info: bool = False,
    continuum_perturb: float = 1e-2,
    zero_tol: float = 1e-10,
    GBZ_check_tol: float = 1e-6,
):
    """Compute the strip winding number and GBZ points at (E_ref, mu1).

    Returns:
        (winding, GBZResult) — winding is float for normal points or
        (W_left, W_right) tuple for continuum points.
        If with_gap_info, a third gap_info element is included.
    """
    gbz_result, theta1_arr, sols_arr, info = get_roots_and_PMGBZ(
        poly_diff, E_ref, mu1, N_points, zero_tol, GBZ_check_tol,
    )
    return _strip_winding_from_result(
        poly_diff, E_ref, mu1,
        gbz_result, theta1_arr, sols_arr, info,
        N_points, with_gap_info, continuum_perturb, zero_tol, GBZ_check_tol,
    )
