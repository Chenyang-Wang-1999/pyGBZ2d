'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2025-11-18
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

''' Strip winding number calculation.

The strip winding number at (E_ref, mu1) is the winding of the
characteristic polynomial f(E_ref, beta1, beta2) along the loop
beta1 = exp(mu1 + i*theta1), beta2 = exp(mu2(theta1) + i*theta2), where
mu2(theta1) is a periodic spline through the middle of the gap between
the M-th and (M+1)-th smallest |beta2| roots.  theta2 is chosen to keep
the loop as far from all roots as possible; PMGBZ points (accidental
equal-modulus crossings) contribute fractional angle corrections.

The zero of the strip winding number as a function of mu1 defines the
SGBZ (see solve_SGBZ_for_E).  Continuum-degenerate mu1 values have no
well-defined winding; the winding value is returned as None.  The
caller is responsible for resolving the continuum via left/right limits
at mu1 ± epsilon.
'''

import numpy as np
from math import pi
from cmath import exp
from typing import Optional
from scipy import interpolate

from .winding import WindingFun, PolyDiffContext, get_winding_number
from .pmgbz_detector import get_roots_and_PMGBZ
from gbz_types import GBZResult

# Winding value: float at normal points, None at continuum points.
StripWinding = Optional[float]
StripWindingResult = tuple[StripWinding, GBZResult]


def get_loop_winding(
    poly_diff: PolyDiffContext,
    E_ref: complex,
    mu1: float,
    mu2_fun: callable,
    theta2: float,
    N_seg: int = 5
) -> float:
    """Winding number of the characteristic polynomial along a mid-gap loop.

    The loop is parameterized by theta1 in [0, 2*pi):
    beta1 = exp(mu1 + i*theta1), beta2 = exp(mu2(theta1) + i*theta2).

    Parameters:
        poly_diff: polynomial evaluation context.
        E_ref: reference energy.
        mu1: log|beta1| of the strip.
        mu2_fun: periodic spline theta1 -> mu2 (mid-gap log-radius of beta2).
        theta2: fixed phase of beta2 along the loop.
        N_seg: number of integration segments passed to get_winding_number.

    Returns:
        Real-valued winding number of f around zero along the loop
        (unrounded; the caller rounds to an integer).
    """
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
) -> StripWindingResult:
    """Compute strip winding from pre-computed root results.

    At a continuum-degenerate mu1 the winding is undefined; the function
    returns None for the winding value.  The caller (solve_SGBZ_for_E) is
    responsible for resolving the continuum via left/right limits.
    """
    # Sort sols_arr by norm at each row
    for row_ind in range(sols_arr.shape[0]):
        sols_arr[row_ind, :] = sols_arr[row_ind, np.argsort(np.abs(sols_arr[row_ind, :]))]

    M = info["M"]
    PMGBZ_raw = info.get("_pmgbz_raw", [])

    if info["continuum_flag"]:
        # Continuum-degenerate mu1 — winding undefined, caller resolves via limits.
        return None, gbz_result

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

    return W_strip, gbz_result


def get_strip_winding(
    poly_diff: PolyDiffContext,
    E_ref: complex,
    mu1: float,
    N_points: int = 301,
    continuum_perturb: float = 1e-2,
    zero_tol: float = 1e-10,
    GBZ_check_tol: float = 1e-6,
    refine_continuum: bool = True,
) -> StripWindingResult:
    """Compute the strip winding number and GBZ points at (E_ref, mu1).

    Convenience wrapper: runs get_roots_and_PMGBZ to solve the beta2 roots
    and detect PMGBZ points, then delegates to _strip_winding_from_result.

    Parameters:
        poly_diff: polynomial evaluation context.
        E_ref: reference energy.
        mu1: log|beta1| of the strip.
        N_points: number of theta1 mesh points on [0, 2*pi).
        continuum_perturb: mu1 offset for continuum left/right limits.
        zero_tol: PMGBZ gap zero-threshold.
        GBZ_check_tol: equal-modulus cluster detection tolerance.
        refine_continuum: forwarded to get_roots_and_PMGBZ.  When False,
            continuum detection skips boundary refinement and stage-3
            accidental-point detection.

    Returns:
        (winding, GBZResult) — winding is float for normal points or
        None for continuum points.
    """
    gbz_result, theta1_arr, sols_arr, info = get_roots_and_PMGBZ(
        poly_diff, E_ref, mu1, N_points, zero_tol, GBZ_check_tol,
        refine_continuum=refine_continuum,
    )
    return _strip_winding_from_result(
        poly_diff, E_ref, mu1,
        gbz_result, theta1_arr, sols_arr, info,
    )
