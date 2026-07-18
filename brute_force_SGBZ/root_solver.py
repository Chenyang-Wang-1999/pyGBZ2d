'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2025-11-18
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

'''
Calculations related to root solving of polynomials
'''

from typing import Optional

import numpy as np
import poly_tools as pt

def poly_to_np_coefficients(coeffs: list[complex], degs: list[int]) -> np.ndarray:
    '''
    Convert a polynomial to a numpy array of coefficients.

    Parameters:
    coeffs (list[complex]): The coefficients of the polynomial.
    degs (list[int]): The degrees of the coefficients.

    Returns:
    np.ndarray: The numpy array of coefficients.
    '''
    max_deg = max(degs)
    np_coeffs = np.zeros(max_deg + 1, dtype=complex)
    for coeff, deg in zip(coeffs, degs):
        np_coeffs[max_deg - deg] = coeff
    return np_coeffs


def calculate_point_roots(
    char_poly: pt.CLaurent,
    param_ind_ctype: pt.CIndexVec,
    param_val: np.ndarray,
    var_ind_ctype: pt.CIndexVec,
    M_max: int,
    N_max: int
):
    poly_1d = char_poly.partial_eval(
        pt.CScalarVec(param_val),
        param_ind_ctype,
        var_ind_ctype
    )
    coeffs = pt.CScalarVec([])
    degs = pt.CIndexVec([])
    poly_1d.num.batch_get_data(coeffs, degs)
    deg_M = poly_1d.denom_orders[0]
    np_coeffs = poly_to_np_coefficients(coeffs, degs)
    curr_roots = list(np.roots(np_coeffs))
    new_root_len = len(curr_roots)
    if new_root_len < M_max + N_max:
        if deg_M < M_max:
            curr_roots.append(0)
        if new_root_len - deg_M < N_max:
            curr_roots.append(np.inf)
    return curr_roots


# ---- mesh root solving (shared by pmgbz_detector and fill_beta2) ----

def _solve_sorted_roots_at(
    poly_diff,
    E_ref: complex,
    mu1: float,
    theta1: float,
    M: int,
    N: int,
    param_ind,
    var_ind,
    cache: Optional[dict] = None,
) -> np.ndarray:
    """Solve and sort beta2 roots at a single theta1 point.

    Parameters:
        poly_diff: PolyDiffContext (provides char_poly).
        E_ref, mu1: fixed energy and log|beta1|.
        theta1: phase angle.
        M, N: minor degrees for root padding (see calculate_point_roots).
        param_ind, var_ind: CIndexVec for partial evaluation.
        cache: optional dict keyed by round(theta1, 14).  When provided,
            hits are returned directly and misses are stored after solving.

    Returns:
        Sorted complex array of beta2 roots.
    """
    from cmath import exp
    from gbz_types import sort_by_root_abs

    theta_norm = float(theta1 % (2 * __import__('math').pi))
    cache_key = round(theta_norm, 14)
    if cache is not None and cache_key in cache:
        return cache[cache_key]

    beta1 = exp(mu1 + 1j * theta_norm)
    roots = calculate_point_roots(
        poly_diff.char_poly, param_ind, (E_ref, beta1), var_ind, M, N,
    )
    result = sort_by_root_abs(roots)
    if cache is not None:
        cache[cache_key] = result
    return result


def solve_roots_on_mesh(
    poly_diff,
    E_ref: complex,
    mu1: float,
    N_points: int = 301,
    extra_thetas: tuple[float, ...] = (),
    cache: Optional[dict] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve beta2 roots on a uniform theta1 mesh + extra points.

    Minimal version of get_roots_and_PMGBZ with NO PMGBZ detection:
    uniform mesh ∪ extra_thetas, calculate_point_roots + sort_by_root_abs
    at each point, periodic closure.

    Parameters:
        poly_diff: polynomial evaluation context (PolyDiffContext).
        E_ref: reference energy.
        mu1: log|beta1|.
        N_points: number of uniform theta1 mesh points on [0, 2π).
        extra_thetas: additional theta1 values to include in the mesh.
        cache: optional dict for memoisation (shared across calls).

    Returns:
        theta1_arr: (K+1,) sorted mesh, closed (last = first + 2π).
        sols_arr: (K+1, M+N) roots sorted by modulus.
    """
    from gbz_types import get_minor_degrees
    from math import pi

    M, N = get_minor_degrees(poly_diff)
    theta_base = np.linspace(0.0, 2 * pi, N_points, endpoint=False)

    param_ind = pt.CIndexVec((0, 1))
    var_ind = pt.CIndexVec([2])

    extra_norm = [float(t) % (2 * pi) for t in extra_thetas]
    all_thetas: list[float] = list(theta_base)
    all_thetas.extend(extra_norm)
    all_unique = np.unique(all_thetas)
    all_unique.sort()

    sols = np.vstack([
        _solve_sorted_roots_at(
            poly_diff, E_ref, mu1, float(t), M, N,
            param_ind, var_ind, cache,
        )
        for t in all_unique
    ])
    theta1_arr = np.hstack((all_unique, [all_unique[0] + 2 * pi]))
    sols_arr = np.vstack((sols, sols[0]))
    return theta1_arr, sols_arr
