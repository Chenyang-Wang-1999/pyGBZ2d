'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-05-19
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

from typing import Optional
import warnings

import numpy as np
import poly_tools as pt
from math import pi
from cmath import exp
from scipy.optimize import linear_sum_assignment, fsolve

from brute_force_solver.root_solver import calculate_point_roots
from brute_force_solver.winding import PolyDiffContext


# ---- inlined utilities to avoid pulling in BerryPy via brute_force_solver.winding ----

def _get_minor_degrees_for_direction(char_poly: pt.CLaurent, direction: int):
    """Get (M, N) for a given direction (1 = beta1, 2 = beta2).

    Reads the numerator polynomial data from char_poly.num, where triplets
    are [E_exp, beta1_exp, beta2_exp].
    """
    coeffs = pt.CScalarVec([])
    degs = pt.CIndexVec([])
    char_poly.num.batch_get_data(coeffs, degs)
    M_plus_N = max(degs[direction::3])
    M = char_poly.denom_orders[direction]
    N = M_plus_N - M
    return M, N


def _get_minor_degrees(char_poly: pt.CLaurent):
    return _get_minor_degrees_for_direction(char_poly, 2)


def _sort_by_root_abs(roots) -> np.ndarray:
    return np.array(sorted(roots, key=lambda x: abs(x)), dtype=complex)


def _chordal_cost_matrix(roots_from: np.ndarray, roots_to: np.ndarray) -> np.ndarray:
    def to_sphere_r3(roots: np.ndarray) -> np.ndarray:
        x = roots.real
        y = roots.imag
        is_inf = np.isinf(x) | np.isinf(y)
        abs_sq = x * x + y * y

        r3 = np.zeros((roots.size, 3), dtype=float)
        finite = ~is_inf
        denom = 1.0 + abs_sq[finite]
        r3[finite, 0] = 2.0 * x[finite] / denom
        r3[finite, 1] = 2.0 * y[finite] / denom
        r3[finite, 2] = (abs_sq[finite] - 1.0) / denom
        r3[is_inf, 2] = 1.0
        return r3

    p_from = to_sphere_r3(roots_from)
    p_to = to_sphere_r3(roots_to)
    diff = p_from[:, None, :] - p_to[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=2), dtype=float)


def _hungarian_match_indices(roots_from: np.ndarray, roots_to: np.ndarray):
    if np.any(np.isnan(roots_from.real)) or np.any(np.isnan(roots_from.imag)):
        raise ValueError(f"NaN root encountered in roots_from: {roots_from}")
    if np.any(np.isnan(roots_to.real)) or np.any(np.isnan(roots_to.imag)):
        raise ValueError(f"NaN root encountered in roots_to: {roots_to}")

    cost = _chordal_cost_matrix(roots_from, roots_to)
    row_ind, col_ind = linear_sum_assignment(cost)
    return [(int(i), int(j)) for i, j in zip(row_ind, col_ind)]


# ---- continuum detection ----

def _find_cyclic_true_intervals(mask: np.ndarray) -> list[tuple[int, int]]:
    """Find contiguous True intervals in a cyclic boolean array.

    Returns list of (start_idx, end_idx_inclusive). Empty list if no True values,
    [(0, n-1)] if all True.
    """
    n = len(mask)
    if n == 0:
        return []
    if np.all(mask):
        return [(0, n - 1)]
    if not np.any(mask):
        return []

    starts = np.where(mask & (~np.roll(mask, 1)))[0]
    intervals = []
    for start in starts:
        end = start
        while mask[end]:
            end = (end + 1) % n
            if end == start:
                break
        intervals.append((int(start), int((end - 1) % n)))
    return intervals


# ---- main API ----

def get_hungarian_sorted_roots(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1: float,
    N_points: int = 301,
):
    """
    Solve beta2 roots on a theta1 mesh and sort using Hungarian matching,
    establishing continuous root tracks across theta1.

    Unlike SGBZ (which sorts roots by |beta2|), this uses Hungarian matching
    to track individual roots continuously as theta1 varies.

    Returns:
        theta1_arr: (N_points,) theta1 mesh
        tracked_roots: (N_points, M+N) root tracks, each column is continuous
        M: denominator order in beta2
        N: numerator max degree minus M
    """
    M, N = _get_minor_degrees(char_poly)
    n_roots = M + N
    theta1_arr = np.linspace(0.0, 2 * pi, N_points, endpoint=False)

    param_ind = pt.CIndexVec((0, 1))
    var_ind = pt.CIndexVec([2])

    all_roots = []
    for theta1 in theta1_arr:
        beta1 = exp(mu1 + 1j * theta1)
        roots = calculate_point_roots(
            char_poly, param_ind, (E_ref, beta1), var_ind, M, N,
        )
        all_roots.append(np.asarray(roots, dtype=complex))

    tracked = np.empty((N_points, n_roots), dtype=complex)
    tracked[0] = _sort_by_root_abs(all_roots[0])

    for i in range(N_points - 1):
        roots_next = all_roots[i + 1]
        matches = _hungarian_match_indices(tracked[i], roots_next)
        reordered = np.zeros(n_roots, dtype=complex)
        for from_idx, to_idx in matches:
            reordered[from_idx] = roots_next[to_idx]
        tracked[i + 1] = reordered

    return theta1_arr, tracked, M, N


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
    M, N = _get_minor_degrees_for_direction(char_poly, direction)

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
    matches_wrap = _hungarian_match_indices(tracked[-1], tracked[0])
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
        intervals = _find_cyclic_true_intervals(near_boundary[:n_pts, j])
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


def _refine_and_correct(
    char_poly, E_ref, mu1, mu2_0, target,
    tracks, low, high, low_init, high_init, f_low_unref, f_high_unref,
    continuum_tol, min_continuum_pts, xtol,
):
    """Refine crossings at mu2_0, then apply Newton correction if needed.

    Uses the analytical dW/dmu2 computed from zero derivatives to take one
    Newton step, avoiding re-bisection after refinement.
    """
    w_ref, zeros_ref, _, dW_dmu2 = _compute_winding_from_tracks(
        char_poly, E_ref, mu1, mu2_0, tracks,
        continuum_tol, min_continuum_pts, refine_crossings=True,
    )
    f_ref = w_ref - target

    if abs(f_ref) < xtol:
        return {
            "mu2": mu2_0, "zeros": zeros_ref, "is_continuum": False,
            "winding": w_ref, "_tracks": tracks,
        }

    # Newton correction using analytical derivative
    if abs(dW_dmu2) > 1e-15:
        delta = -f_ref / dW_dmu2
        # Clamp to initial bracket with margin
        bracket_width = high_init - low_init
        max_step = 0.5 * bracket_width
        delta = max(-max_step, min(max_step, delta))
        mu2_new = mu2_0 + delta
        # Ensure within initial bracket
        mu2_new = max(low_init, min(high_init, mu2_new))
    else:
        # dW/dmu2 near zero — use bisection bracket midpoint
        mu2_new = 0.5 * (low + high)

    # Avoid re-refining at essentially the same point
    if abs(mu2_new - mu2_0) < xtol:
        return {
            "mu2": mu2_0, "zeros": zeros_ref, "is_continuum": False,
            "winding": w_ref, "_tracks": tracks,
        }

    # Refine at the corrected mu2
    w_ref2, zeros_ref2, _, _ = _compute_winding_from_tracks(
        char_poly, E_ref, mu1, mu2_new, tracks,
        continuum_tol, min_continuum_pts, refine_crossings=True,
    )
    f_ref2 = w_ref2 - target

    if abs(f_ref2) < xtol:
        return {
            "mu2": mu2_new, "zeros": zeros_ref2, "is_continuum": False,
            "winding": w_ref2, "_tracks": tracks,
        }

    # Fall back: return the better of the two refined points
    if abs(f_ref2) < abs(f_ref):
        return {
            "mu2": mu2_new, "zeros": zeros_ref2, "is_continuum": False,
            "winding": w_ref2, "_tracks": tracks,
        }
    return {
        "mu2": mu2_0, "zeros": zeros_ref, "is_continuum": False,
        "winding": w_ref, "_tracks": tracks,
    }


def _find_mu2_for_a2_zero(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1: float,
    mu2_low: float,
    mu2_high: float,
    target_winding: float = 0.0,
    N_points: int = 301,
    continuum_tol: float = 1e-8,
    min_continuum_pts: int = 3,
    continuum_perturb: float = 1e-4,
    max_iter: int = 60,
    xtol: float = 1e-10,
    max_range_expansions: int = 10,
    range_expand_factor: float = 2.0,
    _root_tracks: Optional[dict] = None,
) -> dict:
    """Find mu2 where a2 winding crosses target_winding, with adaptive range.

    The a2 winding is monotonic in mu2, so expanding the search range
    guarantees a sign change will eventually be found.

    Uses unrefined (cheap) winding during bisection, then refines crossings
    at the final mu2 and applies a Newton correction using the analytical
    derivative dW/dmu2.

    _root_tracks: optionally pre-computed root tracks from _compute_root_tracks.
    When provided, avoids redundant Hungarian matching.
    """
    low, high = float(mu2_low), float(mu2_high)

    # Pre-compute root tracks once (Independent of mu2)
    if _root_tracks is not None:
        tracks = _root_tracks
    else:
        tracks = _compute_root_tracks(char_poly, E_ref, mu1, N_points)

    def _winding_at(mu2_val, refine=False):
        """Evaluate winding at mu2_val.  Strips dW_dmu2 for the bisection loop."""
        w, z, c, _ = _compute_winding_from_tracks(
            char_poly, E_ref, mu1, mu2_val, tracks,
            continuum_tol, min_continuum_pts,
            refine_crossings=refine,
        )
        return w, z, c

    # Adaptive range expansion (unrefined)
    for _ in range(max_range_expansions):
        w_low, _, _ = _winding_at(low)
        w_high, _, _ = _winding_at(high)
        f_low = w_low - target_winding
        f_high = w_high - target_winding

        if f_low * f_high <= 0:
            low_init, high_init = low, high
            break

        # Expand outward
        width = high - low
        low = low - range_expand_factor * width
        high = high + range_expand_factor * width
    else:
        raise ValueError(f"Maximal range expansion reached! E_ref:{E_ref}, mu1:{mu1}, low:{low}, high:{high}")

    # Bisection with unrefined winding (fast — no fsolve)
    for _ in range(max_iter):
        mu2_mid = 0.5 * (low + high)

        w_mid, zeros, has_continuum = _winding_at(mu2_mid)

        if has_continuum:
            w_left, _, _ = _winding_at(mu2_mid - continuum_perturb)
            w_right, _, _ = _winding_at(mu2_mid + continuum_perturb)
            f_left = w_left - target_winding
            f_right = w_right - target_winding

            if f_left * f_right < 0:
                # Refine zeros at the continuum boundary
                w_ref, zeros_ref, _, _ = _compute_winding_from_tracks(
                    char_poly, E_ref, mu1, mu2_mid, tracks,
                    continuum_tol, min_continuum_pts, refine_crossings=True,
                )
                return {
                    "mu2": mu2_mid, "zeros": zeros_ref, "is_continuum": True,
                    "winding": (w_left, w_right), "_tracks": tracks,
                }
            else:
                if f_low * f_left > 0:
                    low = mu2_mid
                    f_low = f_left
                else:
                    high = mu2_mid
                    f_high = f_left
                continue

        f_mid = w_mid - target_winding

        if abs(f_mid) < xtol or (high - low) < xtol:
            # --- Post-refinement + Newton correction ---
            return _refine_and_correct(
                char_poly, E_ref, mu1, mu2_mid, target_winding,
                tracks, low, high, low_init, high_init, f_low, f_high,
                continuum_tol, min_continuum_pts, xtol,
            )

        if f_low * f_mid < 0:
            high = mu2_mid
            f_high = f_mid
        else:
            low = mu2_mid
            f_low = f_mid

    # Max iterations reached — refine at the final midpoint
    mu2_mid = 0.5 * (low + high)
    w_ref, zeros_ref, _, _ = _compute_winding_from_tracks(
        char_poly, E_ref, mu1, mu2_mid, tracks,
        continuum_tol, min_continuum_pts, refine_crossings=True,
    )
    return {
        "mu2": mu2_mid, "zeros": zeros_ref, "is_continuum": False,
        "winding": w_ref, "_tracks": tracks,
    }


def bisect_a2_winding(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1: float,
    mu2_low: float,
    mu2_high: float,
    target_winding: float = 0.0,
    N_points: int = 301,
    continuum_tol: float = 1e-8,
    min_continuum_pts: int = 3,
    continuum_perturb: float = 1e-4,
    max_iter: int = 60,
    xtol: float = 1e-10,
) -> dict:
    """
    Bisect mu2 to find where the average a2 winding crosses target_winding.

    Thin wrapper around _find_mu2_for_a2_zero without adaptive range
    expansion (the caller provides the bracket).

    Returns a dict with keys:
        mu2: critical mu2 value
        zeros: list of (theta1, theta2) crossing pairs at mu2
        is_continuum: whether the result is a continuum point
        winding: average winding at the returned mu2.
                 For continuum points this is a tuple (w_left, w_right).
    """
    return _find_mu2_for_a2_zero(
        char_poly, E_ref, mu1, mu2_low, mu2_high,
        target_winding, N_points, continuum_tol, min_continuum_pts,
        continuum_perturb, max_iter, xtol,
        max_range_expansions=0,  # no adaptive expansion for the thin wrapper
    )


def _a1_is_degenerate(zeros: list[tuple[float, float]]) -> bool:
    """Heuristic: a1 winding is degenerate if zeros can't partition theta2.

    True when all zero points share essentially the same theta2 value,
    meaning the a2 crossings don't provide enough angular diversity to
    segment theta2 for a meaningful a1 winding computation.
    """
    if len(zeros) <= 1:
        return True
    theta2_vals = np.unique([z[1] for z in zeros])
    return len(theta2_vals) <= 1


def bisect_amoeba_ronkin_min(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1_low: float = -1,
    mu1_high: float = 1,
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
) -> dict:
    """
    Find the Ronkin function minimum by bisecting mu1 and mu2.

    Outer loop: bisect mu1.
    Inner loop: for each mu1, find mu2 where a2 average winding = 0,
    then evaluate a1 average winding at (mu1, mu2).

    The Ronkin minimum satisfies a1 = a2 = 0 simultaneously.

    Continuum handling:
    - When a1 is degenerate at (mu1_mid, mu2_0), perturb mu1 ± epsilon.
      For each perturbed mu1, re-run the inner mu2 bisection to find a2=0,
      then compute a1.
      * Opposite signs → this is the boundary, stop.
      * Same sign → use the sign to continue the outer bisection.

    Returns a dict with keys:
        mu1: critical mu1 value
        mu2: critical mu2 value
        zeros: list of (theta1, theta2) crossing pairs at the critical point
        is_continuum: whether the result is a continuum point
    """
    # Evaluate a1 at the mu1 endpoints, with adaptive range expansion
    low, high = float(mu1_low), float(mu1_high)
    a1_low = a1_high = 0.0
    inner_low = inner_high = None

    for _ in range(max_range_expansions):
        inner_low = _find_mu2_for_a2_zero(
            char_poly, E_ref, low, mu2_low, mu2_high,
            target_winding=0.0, N_points=N_points,
            continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
            continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
            max_range_expansions=max_range_expansions,
            range_expand_factor=range_expand_factor,
        )
        inner_high = _find_mu2_for_a2_zero(
            char_poly, E_ref, high, mu2_low, mu2_high,
            target_winding=0.0, N_points=N_points,
            continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
            continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
            max_range_expansions=max_range_expansions,
            range_expand_factor=range_expand_factor,
        )

        a1_low = _get_average_winding_from_zeros(
            char_poly, E_ref, low, inner_low["mu2"],
            inner_low["zeros"], direction=1,
        )
        a1_high = _get_average_winding_from_zeros(
            char_poly, E_ref, high, inner_high["mu2"],
            inner_high["zeros"], direction=1,
        )

        if a1_low * a1_high <= 0:
            break

        width = high - low
        low = low - range_expand_factor * width
        high = high + range_expand_factor * width
    else:
        raise ValueError(f"bisect_amoeba_ronkin_min: reach maximum expansion. E_ref:{E_ref}, low:{low}, high:{high}")

    mu1_low, mu1_high = low, high

    for _ in range(max_iter):
        mu1_mid = 0.5 * (mu1_low + mu1_high)

        inner_mid = _find_mu2_for_a2_zero(
            char_poly, E_ref, mu1_mid, mu2_low, mu2_high,
            target_winding=0.0, N_points=N_points,
            continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
            continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
            max_range_expansions=max_range_expansions,
            range_expand_factor=range_expand_factor,
        )

        mu2_mid = inner_mid["mu2"]
        zeros_mid = inner_mid["zeros"]

        # Check for a1 degeneracy
        if _a1_is_degenerate(zeros_mid):
            # a1 may already be at zero (plateau) without crossing sign.
            # Catch this before perturbing to avoid stalling the bisection.
            a1_mid = _get_average_winding_from_zeros(
                char_poly, E_ref, mu1_mid, mu2_mid,
                zeros_mid, direction=1,
            )
            if abs(a1_mid) < xtol:
                return {
                    "mu1": mu1_mid, "mu2": mu2_mid,
                    "zeros": zeros_mid,
                    "is_continuum": inner_mid["is_continuum"],
                }

            inner_left = _find_mu2_for_a2_zero(
                char_poly, E_ref, mu1_mid - continuum_perturb,
                mu2_low, mu2_high,
                target_winding=0.0, N_points=N_points,
                continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
                continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
                max_range_expansions=max_range_expansions,
                range_expand_factor=range_expand_factor,
            )
            inner_right = _find_mu2_for_a2_zero(
                char_poly, E_ref, mu1_mid + continuum_perturb,
                mu2_low, mu2_high,
                target_winding=0.0, N_points=N_points,
                continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
                continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
                max_range_expansions=max_range_expansions,
                range_expand_factor=range_expand_factor,
            )

            a1_left = _get_average_winding_from_zeros(
                char_poly, E_ref, mu1_mid - continuum_perturb,
                inner_left["mu2"], inner_left["zeros"], direction=1,
            )
            a1_right = _get_average_winding_from_zeros(
                char_poly, E_ref, mu1_mid + continuum_perturb,
                inner_right["mu2"], inner_right["zeros"], direction=1,
            )

            if a1_left * a1_right < 0:
                return {
                    "mu1": mu1_mid, "mu2": mu2_mid,
                    "zeros": zeros_mid,
                    "is_continuum": True,
                }
            elif a1_left < 0 and a1_right < 0:
                mu1_low = mu1_mid + continuum_perturb
                a1_low = a1_right
            else:
                mu1_high = mu1_mid - continuum_perturb
                a1_high = a1_left
            continue

        # Normal case: compute a1 at (mu1_mid, mu2_mid)
        # Reuse zeros from inner bisection — no extra Hungarian matching
        a1_mid = _get_average_winding_from_zeros(
            char_poly, E_ref, mu1_mid, mu2_mid,
            zeros_mid, direction=1,
        )

        if abs(a1_mid) < xtol or (mu1_high - mu1_low) < xtol:
            return {
                "mu1": mu1_mid, "mu2": mu2_mid,
                "zeros": zeros_mid,
                "is_continuum": inner_mid["is_continuum"],
            }

        if a1_low * a1_mid < 0:
            mu1_high = mu1_mid
            a1_high = a1_mid
        else:
            mu1_low = mu1_mid
            a1_low = a1_mid

    mu1_mid = 0.5 * (mu1_low + mu1_high)
    inner_mid = _find_mu2_for_a2_zero(
        char_poly, E_ref, mu1_mid, mu2_low, mu2_high,
        target_winding=0.0, N_points=N_points,
        continuum_tol=continuum_tol, min_continuum_pts=min_continuum_pts,
        continuum_perturb=continuum_perturb, max_iter=max_iter, xtol=xtol,
        max_range_expansions=max_range_expansions,
        range_expand_factor=range_expand_factor,
    )
    raise ValueError(f"bisect_amoeba_ronkin_min: Bisection failed. E_ref{E_ref}")


def check_amoeba(
    coeffs: np.ndarray,
    degs: np.ndarray,
    E_ref: complex,
    perc: float,
    debug_mode: bool = False,
    **options,
):
    ''' Check amoeba condition for a given reference enenrgy. '''
    print("%.2f" % (perc * 100) + r"%")
    char_poly = pt.CLaurent(3)
    coeffs = pt.CScalarVec(coeffs)
    degs = pt.CLaurentIndexVec(degs.flatten())
    char_poly.set_Laurent_by_terms(coeffs, degs)

    try:
        amoeba_res = bisect_amoeba_ronkin_min(
            char_poly, E_ref, **options
        )
        if (not amoeba_res["is_continuum"]) and (len(amoeba_res["zeros"]) == 0):
            is_amoeba = False
        else:
            is_amoeba = True

        res = {
            "success": True,
            "is_amoeba": is_amoeba,
        }
        res.update(amoeba_res)
        return res

    except Exception as e:
        if debug_mode:
            raise e
        return {"success": False, "error": str(e)}
