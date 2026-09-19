"""
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-14
Copyright © Department of Physics, Tsinghua University. All rights reserved

Ronkin-function winding helpers for amoeba GBZ.

The mesh-based winding source (``_compute_winding_from_tracks`` and friends)
was removed when the zero-solving layer moved to ``continuation.ZeroManager``;
the winding is now computed from ZM zeros via
``zm_extract.calculate_a2_average_winding``.

Retained here are the zero- and polynomial-level primitives that the new
pipeline still consumes:
  - ``_find_exact_crossing``         fsolve refinement of a (θ₁, θ₂) crossing.
  - ``_get_average_winding_from_zeros``  average winding from a zero partition.
"""

from typing import Optional
import warnings
import numpy as np
from cmath import exp
from scipy.optimize import fsolve

# fsolve refinement settings for amoeba (θ₁, θ₂) crossings.
FSOLVE_XTOL: float = 1e-12
FSOLVE_MAXFEV: int = 500
#: Residual gate that accepts a refined crossing.
CROSSING_RESIDUAL_TOL: float = 1e-10
from pygbz2d.core import CharPoly, TWO_PI


def _find_exact_crossing(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    mu2: float,
    theta1_guess: float,
    theta2_guess: float,
) -> Optional[tuple[float, float]]:
    """
    Refine a crossing point using ``scipy.optimize.fsolve`` with the exact
    Jacobian.

    Solves f(E, exp(mu1 + i*t1), exp(mu2 + i*t2)) = 0 for (t1, t2).

    Returns (theta1, theta2) if converged, None if the root finder failed.
    """
    def func(theta):
        t1, t2 = theta
        beta1 = exp(mu1 + 1j * t1)
        beta2 = exp(mu2 + 1j * t2)
        val = poly.eval_val((E_ref, beta1, beta2))
        return [val.real, val.imag]

    def jac(theta):
        t1, t2 = theta
        beta1 = exp(mu1 + 1j * t1)
        beta2 = exp(mu2 + 1j * t2)
        partials = poly.eval_partials((E_ref, beta1, beta2))
        df_dt1 = 1j * beta1 * partials[1]
        df_dt2 = 1j * beta2 * partials[2]
        return [[df_dt1.real, df_dt2.real],
                [df_dt1.imag, df_dt2.imag]]

    # fsolve may emit RuntimeWarning when the initial guess falls in a
    # flat or ill-conditioned region (e.g. near band edges explored by
    # the adaptive bisection).  The residual check below is the actual
    # quality gate; suppress the scipy warning since non-convergence is
    # handled correctly by returning None and falling back to the
    # unrefined linear-interpolation estimate.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol = fsolve(func, [theta1_guess, theta2_guess], fprime=jac,
                     xtol=FSOLVE_XTOL, maxfev=FSOLVE_MAXFEV)
    t1, t2 = sol
    beta1 = exp(mu1 + 1j * t1)
    beta2 = exp(mu2 + 1j * t2)
    residual = abs(poly.eval_val((E_ref, beta1, beta2)))
    if residual < CROSSING_RESIDUAL_TOL:
        return (float(t1 % (TWO_PI)), float(t2 % (TWO_PI)))
    # return None
    raise RuntimeError(f"fsolve failed to converge after {FSOLVE_MAXFEV} iterations: ")


def _get_average_winding_from_zeros(
    char_poly: CharPoly,
    E_ref: complex,
    mu1: float,
    mu2: float,
    zeros: list[tuple[float, float]],  # (theta1, theta2)
    direction: int,  # 1 = w1 (solve beta1), 2 = w2 (solve beta2)
) -> tuple[float, float]:
    """Compute average winding and normalized non-zero interval area.

    Returns (avg_winding, non_zero_area) where both are divided by (2π):

      avg_winding   = Σ(u · width) / (2π)
      non_zero_area = Σ(|u| · width) / (2π)

    The zeros partition the angular circle in the `direction` variable.
    u_d is constant on each segment between consecutive zero crossings.

    direction=2 (w2): partition theta1, solve beta2 at each segment midpoint.
    direction=1 (w1): partition theta2, solve beta1 at each segment midpoint.

    non_zero_area serves as a lightweight plateau indicator: at a
    zero-plateau boundary, both w1 and w2 non-zero areas are tiny
    (most of the circle has u = 0).  Callers can skip expensive plateau
    probing when either area exceeds a safe threshold (e.g. 1e-2).
    """
    M, N = char_poly.get_minor_degrees(direction)

    if direction == 2:
        partition_thetas = np.unique([z[0] for z in zeros])
        mu_solve = mu1
        mu_count = mu2
        param_inds, var_inds = (0, 1), (2,)
    else:
        partition_thetas = np.unique([z[1] for z in zeros])
        mu_solve = mu2
        mu_count = mu1
        param_inds, var_inds = (0, 2), (1,)

    partition_thetas.sort()
    n_seg = len(partition_thetas)

    if n_seg == 0:
        if direction == 2:
            beta_param = exp(mu1 + 1j * 0.0)
        else:
            beta_param = exp(mu2 + 1j * 0.0)
        roots = char_poly.solve_roots_1d(
            param_inds, (E_ref, beta_param), var_inds, M, N,
        )
        count_below = np.sum(np.log(np.abs(np.asarray(roots, dtype=complex))) < mu_count)
        u_const = count_below - M
        # Single segment spans the full circle: width = 2π, so non_zero_area / (2π) = |u|
        return u_const, abs(u_const)

    total = 0.0
    non_zero_area = 0.0
    for i in range(n_seg):
        left = partition_thetas[i]
        right = partition_thetas[(i + 1) % n_seg]
        if i == n_seg - 1:
            right += TWO_PI
        width = right - left
        mid = (0.5 * (left + right)) % (TWO_PI)

        beta_param = exp(mu_solve + 1j * mid)
        roots = char_poly.solve_roots_1d(
            param_inds, (E_ref, beta_param), var_inds, M, N,
        )
        count_below = np.sum(np.log(np.abs(np.asarray(roots, dtype=complex))) < mu_count)
        u = count_below - M
        total += u * width
        non_zero_area += abs(u) * width

    return total / (TWO_PI), non_zero_area / (TWO_PI)
