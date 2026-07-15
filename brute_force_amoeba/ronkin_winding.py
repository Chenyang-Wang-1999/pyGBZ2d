"""
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-14
Copyright © Department of Physics, Tsinghua University. All rights reserved

Ronkin-function winding number computation for amoeba GBZ.

Provides a1 and a2 average winding numbers, zero-crossing detection, and
root-track-based winding computation.
"""

from typing import Optional
import warnings
import numpy as np
import poly_tools as pt
from math import pi
from cmath import exp
from scipy.optimize import fsolve

from brute_force_SGBZ.root_solver import calculate_point_roots
from brute_force_SGBZ.winding import PolyDiffContext
from gbz_types import (
    get_minor_degrees, find_cyclic_true_intervals,
    hungarian_match_indices,
)

from .tracks import get_hungarian_sorted_roots, _compute_root_tracks

def _find_exact_crossing(
    poly_diff: PolyDiffContext,
    E_ref: complex,
    mu1: float,
    mu2: float,
    theta1_guess: float,
    theta2_guess: float,
) -> Optional[tuple[float, float]]:
    """
    Refine a crossing point using scipy.optimize.root with the exact Jacobian.

    Solves f(E, exp(mu1 + i*t1), exp(mu2 + i*t2)) = 0 for (t1, t2).

    Returns (theta1, theta2) if converged, None if the root finder failed.
    """
    def func(theta):
        t1, t2 = theta
        beta1 = exp(mu1 + 1j * t1)
        beta2 = exp(mu2 + 1j * t2)
        val = poly_diff.eval_val((E_ref, beta1, beta2))
        return [val.real, val.imag]

    def jac(theta):
        t1, t2 = theta
        beta1 = exp(mu1 + 1j * t1)
        beta2 = exp(mu2 + 1j * t2)
        partials = poly_diff.eval_partials((E_ref, beta1, beta2))
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
                     xtol=1e-12, maxfev=500)
    t1, t2 = sol
    beta1 = exp(mu1 + 1j * t1)
    beta2 = exp(mu2 + 1j * t2)
    residual = abs(poly_diff.eval_val((E_ref, beta1, beta2)))
    if residual < 1e-10:
        return (float(t1 % (2 * pi)), float(t2 % (2 * pi)))
    return None


def _get_average_winding_from_zeros(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1: float,
    mu2: float,
    zeros: list[tuple[float, float]],  # (theta1, theta2)
    direction: int,  # 1 = a1 (solve beta1), 2 = a2 (solve beta2)
) -> float:
    """Compute average winding from pre-computed (theta1, theta2) zero points.

    The zeros partition the angular circle in the `direction` variable.
    u_d is constant on each segment between consecutive zero crossings.

    direction=2 (a2): partition theta1, solve beta2 at each segment midpoint.
    direction=1 (a1): partition theta2, solve beta1 at each segment midpoint.
    """
    M, N = get_minor_degrees(PolyDiffContext(char_poly), direction)

    if direction == 2:
        partition_thetas = np.unique([z[0] for z in zeros])
        mu_solve = mu1
        mu_count = mu2
        param_ind = pt.CIndexVec((0, 1))
        var_ind = pt.CIndexVec([2])
    else:
        partition_thetas = np.unique([z[1] for z in zeros])
        mu_solve = mu2
        mu_count = mu1
        param_ind = pt.CIndexVec((0, 2))
        var_ind = pt.CIndexVec([1])

    partition_thetas.sort()
    n_seg = len(partition_thetas)

    if n_seg == 0:
        if direction == 2:
            beta_param = exp(mu1 + 1j * 0.0)
        else:
            beta_param = exp(mu2 + 1j * 0.0)
        roots = calculate_point_roots(
            char_poly, param_ind, (E_ref, beta_param), var_ind, M, N,
        )
        count_below = np.sum(np.log(np.abs(np.asarray(roots, dtype=complex))) < mu_count)
        return float(count_below - M)

    total = 0.0
    for i in range(n_seg):
        left = partition_thetas[i]
        right = partition_thetas[(i + 1) % n_seg]
        if i == n_seg - 1:
            right += 2 * pi
        width = right - left
        mid = (0.5 * (left + right)) % (2 * pi)

        beta_param = exp(mu_solve + 1j * mid)
        roots = calculate_point_roots(
            char_poly, param_ind, (E_ref, beta_param), var_ind, M, N,
        )
        count_below = np.sum(np.log(np.abs(np.asarray(roots, dtype=complex))) < mu_count)
        u = count_below - M
        total += u * width

    return total / (2 * pi)


def _compute_root_tracks(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1: float,
    N_points: int = 301,
) -> dict:
    """Pre-compute Hungarian-matched beta2 root tracks across theta1.

    This is the expensive part of the pipeline — depends only on (E, mu1),
    NOT on mu2.  Cache this when scanning over mu2 at fixed (E, mu1).

    Returns a dict with keys:
        theta1_arr, tracked, M, N, theta1_ext, tracked_ext, poly_diff
    """
    theta1_arr, tracked, M, N = get_hungarian_sorted_roots(
        char_poly, E_ref, mu1, N_points,
    )
    theta1_ext = np.hstack([theta1_arr, [theta1_arr[0] + 2 * pi]])

    # Ensure the periodic extension respects physical root identity, not
    # just array index.  Hungarian matching tracks roots by chordal
    # distance on the Riemann sphere; when |beta2| ordering crosses
    # (a genuine PMGBZ-like degeneracy), tracked[-1, j] and tracked[0, j]
    # may belong to *different* physical roots.  Naively stacking
    # tracked[0] at the end would connect unrelated roots across the
    # wrap-around boundary, producing spurious ln|beta2| crossings and
    # corrupting both the coarse crossing detection and the fsolve
    # refinement step.
    matches_wrap = hungarian_match_indices(tracked[-1], tracked[0])
    first_periodic = np.zeros(tracked.shape[1], dtype=complex)
    for from_idx, to_idx in matches_wrap:
        first_periodic[from_idx] = tracked[0, to_idx]
    tracked_ext = np.vstack([tracked, first_periodic[np.newaxis, :]])
    return {
        "theta1_arr": theta1_arr,
        "tracked": tracked,
        "M": M,
        "N": N,
        "theta1_ext": theta1_ext,
        "tracked_ext": tracked_ext,
        "poly_diff": PolyDiffContext(char_poly),
    }


def _compute_zero_dtheta1_dmu2(
    poly_diff: PolyDiffContext,
    E_ref: complex,
    mu1: float,
    mu2: float,
    theta1: float,
    theta2: float,
) -> float:
    """Compute d(theta1)/d(mu2) at a refined zero (theta1, theta2).

    From the implicit equation f(E, beta1, beta2) = 0 with beta_j = exp(mu_j + i*theta_j),
    treating theta1, theta2 as functions of mu2:

        i*a*theta1_dot + i*b*theta2_dot = -b

    where a = df/dbeta1 * beta1, b = df/dbeta2 * beta2.

    Solves the 2x2 real linear system for theta1_dot.
    Uses eval_partials data that is already computed during fsolve refinement.
    """
    beta1 = exp(mu1 + 1j * theta1)
    beta2 = exp(mu2 + 1j * theta2)
    partials = poly_diff.eval_partials((E_ref, beta1, beta2))
    a = partials[1] * beta1  # df/dbeta1 * beta1
    b = partials[2] * beta2  # df/dbeta2 * beta2

    # Solve [[Re(a), Re(b)], [Im(a), Im(b)]] * [theta1_dot, theta2_dot] = [-Im(b), Re(b)]
    det = a.real * b.imag - a.imag * b.real
    if abs(det) < 1e-30:
        return 0.0
    # Cramer's rule for theta1_dot
    return float(-(b.real * b.real + b.imag * b.imag) / det)


def _compute_winding_from_tracks(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1: float,
    mu2: float,
    tracks: dict,
    continuum_tol: float = 1e-8,
    min_continuum_pts: int = 3,
    refine_crossings: bool = True,
):
    """Detect a2 crossings and compute a2 average winding from pre-computed tracks.

    Uses pre-computed root tracks (from _compute_root_tracks) — does NOT redo
    Hungarian matching.  Safe to call many times at different mu2 for the same
    (E, mu1).

    When refine_crossings=False, skips fsolve refinement and uses linear
    interpolation directly — much faster, suitable for bisection iterations
    where only the winding value (not exact zero positions) is needed.

    Returns:
        winding: float, a2 average winding number
        crossing_pairs: list of (theta1, theta2) crossing points
        has_continuum: bool, whether extended continuum regions exist
        dW_dmu2: float, analytical derivative d(winding)/d(mu2) (only when
                 refine_crossings=True and crossings exist, else 0.0)
    """
    theta1_arr = tracks["theta1_arr"]
    tracked = tracks["tracked"]
    M = tracks["M"]
    N = tracks["N"]
    theta1_ext = tracks["theta1_ext"]
    tracked_ext = tracks["tracked_ext"]
    poly_diff = tracks["poly_diff"]

    n_roots = M + N
    n_pts = len(theta1_arr)

    # Per-root continuum mask: True where |ln|beta2_j| - mu2| < continuum_tol
    near_boundary = np.abs(np.log(np.abs(tracked_ext)) - mu2) < continuum_tol

    # Check for extended continuum regions
    has_continuum = False
    for j in range(n_roots):
        intervals = find_cyclic_true_intervals(near_boundary[:n_pts, j])
        for start, end in intervals:
            width = (end - start) % n_pts + 1
            if width >= min_continuum_pts:
                has_continuum = True
                break
        if has_continuum:
            break

    # Find all approximate crossings, filtering out continuum noise.
    # Tag each crossing with its winding jump direction:
    #   jump = +1: root crosses mu2 going UP   (log|beta2| increasing, winding decreases)
    #   jump = -1: root crosses mu2 going DOWN (log|beta2| decreasing, winding increases)
    # The jump is used to compute dW/dmu2 = (1/2pi) * sum(jump_i * dtheta1_i/dmu2).
    crossing_approx = []  # (theta1_approx, beta2_approx, jump)
    for j in range(n_roots):
        log_abs = np.log(np.abs(tracked_ext[:, j]))
        near_j = near_boundary[:, j]
        for i in range(n_pts):
            d0 = log_abs[i] - mu2
            d1 = log_abs[i + 1] - mu2
            if d0 * d1 < 0:
                if near_j[i] and near_j[i + 1]:
                    continue  # noise crossing inside continuum band
                frac = (mu2 - log_abs[i]) / (log_abs[i + 1] - log_abs[i])
                theta1_approx = theta1_ext[i] + frac * (theta1_ext[i + 1] - theta1_ext[i])
                beta2_approx = tracked_ext[i, j] + frac * (tracked_ext[i + 1, j] - tracked_ext[i, j])
                jump = 1 if d0 < 0 else -1
                crossing_approx.append((theta1_approx, beta2_approx, jump))

    if not crossing_approx:
        winding = _get_average_winding_from_zeros(
            char_poly, E_ref, mu1, mu2, [], direction=2,
        )
        return winding, [], has_continuum, 0.0

    # Refine or skip
    if refine_crossings:
        refined_full = []  # (theta1, theta2, jump)
        for theta1_approx, beta2_approx, jump in crossing_approx:
            theta2_guess = float(np.angle(beta2_approx))
            result = _find_exact_crossing(
                poly_diff, E_ref, mu1, mu2,
                theta1_approx % (2 * pi), theta2_guess,
            )
            if result is not None:
                refined_full.append((result[0], result[1], jump))
            else:
                refined_full.append((theta1_approx % (2 * pi), theta2_guess, jump))
    else:
        refined_full = [
            (theta1_approx % (2 * pi), float(np.angle(beta2_approx)) % (2 * pi), jump)
            for theta1_approx, beta2_approx, jump in crossing_approx
        ]

    # Deduplicate crossing pairs (by theta1)
    dW_dmu2 = 0.0
    if refine_crossings:
        for t1, t2, jump in sorted(refined_full, key=lambda x: x[0]):
            theta1_dot = _compute_zero_dtheta1_dmu2(
                poly_diff, E_ref, mu1, mu2, t1, t2,
            )
            dW_dmu2 += jump * theta1_dot

        dW_dmu2 /= (2 * pi)

    # Compute a2 average winding from the zero points
    winding = _get_average_winding_from_zeros(
        char_poly, E_ref, mu1, mu2, refined_full, direction=2,
    )

    return winding, refined_full, has_continuum, dW_dmu2


def _compute_crossings_and_winding(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1: float,
    mu2: float,
    N_points: int = 301,
    continuum_tol: float = 1e-8,
    min_continuum_pts: int = 3,
):
    """Detect a2 crossings and compute a2 average winding.

    Thin wrapper: computes root tracks then detects crossings.
    For repeated calls at different mu2 (same E, mu1), use
    _compute_root_tracks + _compute_winding_from_tracks directly.
    """
    tracks = _compute_root_tracks(char_poly, E_ref, mu1, N_points)
    return _compute_winding_from_tracks(
        char_poly, E_ref, mu1, mu2, tracks, continuum_tol, min_continuum_pts,
    )


def get_a2_average_winding(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1: float,
    mu2: float,
    N_points: int = 301,
    continuum_tol: float = 1e-8,
    min_continuum_pts: int = 3,
) -> float:
    """
    Compute the a2-direction average winding number at given (E, mu1, mu2).

    The a2 winding number u2(theta1) = #{beta2 roots with |beta2| < exp(mu2)} - M.
    Crossings of ln|beta2| = mu2 on root tracks partition theta1 into segments
    where u2 is constant. The average = sum(segment_u2 * segment_width) / (2*pi).

    Continuum handling: when a root track stays at |beta2| ≈ exp(mu2) over an
    extended theta1 range, spurious crossings from numerical noise are filtered
    out.  The winding at exactly a continuum mu2 corresponds to one limit value;
    callers should perturb mu2 to resolve the ambiguity.

    This corresponds to dR/d(mu2), the mu2-derivative of the Ronkin function.
    """
    winding, _, _, _ = _compute_crossings_and_winding(
        char_poly, E_ref, mu1, mu2, N_points, continuum_tol, min_continuum_pts,
    )
    return winding


def get_a1_average_winding(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1: float,
    mu2: float,
    N_points: int = 301,
    continuum_tol: float = 1e-8,
    min_continuum_pts: int = 3,
) -> float:
    """
    Compute the a1-direction average winding number at given (E, mu1, mu2).

    Reuses the (theta1, theta2) zero points computed by the a2 crossing
    detection.  These zeros' theta2 values partition theta2; u1 is constant
    on each segment.  No separate Hungarian matching is needed.

    This corresponds to dR/d(mu1), the mu1-derivative of the Ronkin function.
    """
    _, zeros, _, _ = _compute_crossings_and_winding(
        char_poly, E_ref, mu1, mu2, N_points, continuum_tol, min_continuum_pts,
    )
    return _get_average_winding_from_zeros(
        char_poly, E_ref, mu1, mu2, zeros, direction=1,
    )

