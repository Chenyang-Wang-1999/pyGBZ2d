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
from gbz_types import CharPoly, sort_by_root_abs


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
    char_poly: CharPoly,
    param_indices,
    param_val: np.ndarray,
    var_indices,
    M_max: int,
    N_max: int
):
    """Solve roots at a point.  Thin wrapper delegating to CharPoly.

    Kept for backward compatibility.  Prefer ``char_poly.solve_roots_1d()``
    in new code.
    """
    return char_poly.solve_roots_1d(param_indices, param_val, var_indices, M_max, N_max)


# ---- mesh root solving (shared by pmgbz_detector and fill_beta2) ----

def _solve_sorted_roots_at(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    theta1: float,
    M: Optional[int] = None,
    N: Optional[int] = None,
    cache: Optional[dict] = None,
) -> np.ndarray:
    """Solve and sort beta2 roots at a single theta1 point.

    Parameters:
        poly: CharPoly instance.
        E_ref, mu1: fixed energy and log|beta1|.
        theta1: phase angle.
        M, N: minor degrees for root padding.  Defaults to ``poly.M, poly.N``.
        cache: optional dict keyed by round(theta1, 14).  When provided,
            hits are returned directly and misses are stored after solving.

    Returns:
        Sorted complex array of beta2 roots.
    """
    from cmath import exp
    from math import pi as _pi

    theta_norm = float(theta1 % (2 * _pi))
    cache_key = round(theta_norm, 14)
    if cache is not None and cache_key in cache:
        return cache[cache_key]

    beta1 = exp(mu1 + 1j * theta_norm)
    roots = poly.solve_roots_1d((0, 1), (E_ref, beta1), (2,), M, N)
    result = sort_by_root_abs(roots)
    if cache is not None:
        cache[cache_key] = result
    return result


def solve_roots_on_mesh(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    N_points: int = 301,
    extra_thetas: tuple[float, ...] = (),
    cache: Optional[dict] = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve beta2 roots on a uniform theta1 mesh + extra points.

    Minimal version of get_roots_and_PMGBZ with NO PMGBZ detection:
    uniform mesh ∪ extra_thetas, solve_roots_1d + sort_by_root_abs
    at each point, periodic closure.

    Parameters:
        poly: CharPoly instance.
        E_ref: reference energy.
        mu1: log|beta1|.
        N_points: number of uniform theta1 mesh points on [0, 2π).
        extra_thetas: additional theta1 values to include in the mesh.
        cache: optional dict for memoisation (shared across calls).

    Returns:
        theta1_arr: (K+1,) sorted mesh, closed (last = first + 2π).
        sols_arr: (K+1, M+N) roots sorted by modulus.
    """
    from math import pi

    theta_base = np.linspace(0.0, 2 * pi, N_points, endpoint=False)

    extra_norm = [float(t) % (2 * pi) for t in extra_thetas]
    all_thetas: list[float] = list(theta_base)
    all_thetas.extend(extra_norm)
    all_unique = np.unique(all_thetas)
    all_unique.sort()

    sols = np.vstack([
        _solve_sorted_roots_at(poly, E_ref, mu1, float(t), cache=cache)
        for t in all_unique
    ])
    theta1_arr = np.hstack((all_unique, [all_unique[0] + 2 * pi]))
    sols_arr = np.vstack((sols, sols[0]))
    return theta1_arr, sols_arr
