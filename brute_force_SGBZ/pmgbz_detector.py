'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2025-11-18 00:00:00
Copyright © YourCompanyName All rights reserved
'''

''' PMGBZ detection pipeline — root solving, continuum detection, accidental point refinement.

For fixed (E_ref, mu1), solve the beta2 roots of the characteristic
polynomial on a theta1 mesh around the loop |beta1| = exp(mu1), sort them
by modulus, and find every theta1 where the PMGBZ condition
|beta2_M| == |beta2_{M+1}| holds (M = denominator order in beta2):

- continuum intervals — whole theta1 ranges of degeneracy, boundaries
  refined by bisection on the degeneracy predicate → LineSubset;
- accidental points — isolated crossings where root moduli exchange order,
  located by Hungarian-matching bisection between mesh points, then
  validated (fake-crossing filtering via analytic dmu2/dtheta1 or a
  finite-difference fallback) → PointSubset;
- boundary double roots — near-degenerate root pairs that defeat confident
  matching → PointSubset with a "zero" classification.

Entry point: get_roots_and_PMGBZ.  Everything else is a helper.
'''

import numpy as np
import poly_tools as pt
from math import pi
from cmath import exp
from itertools import chain
from scipy.optimize import linear_sum_assignment

from .root_solver import _solve_sorted_roots_at, solve_roots_on_mesh
from .winding import PolyDiffContext

from gbz_types import (
    PointSubset, LineSubset, GBZResult,
    get_minor_degrees,
    chordal_cost_matrix, hungarian_match_indices,
    find_cyclic_true_intervals,
)

def _pmgbz_gap_and_scale(sorted_roots: np.ndarray, M: int) -> tuple[float, float]:
    abs_roots = np.abs(sorted_roots)
    lower = abs_roots[M - 1]
    upper = abs_roots[M]
    scale = max(1.0, 0.5 * (lower + upper))
    return upper - lower, scale

def _is_pmgbz_degenerate(sorted_roots: np.ndarray, M: int, zero_tol: float) -> bool:
    gap, scale = _pmgbz_gap_and_scale(sorted_roots, M)
    return gap <= zero_tol * scale

def _cross_boundary_matches(roots_left: np.ndarray, roots_right: np.ndarray, M: int) -> list[tuple[int, int]]:
    matches = hungarian_match_indices(roots_left, roots_right)
    return [
        (i, j)
        for i, j in matches
        if (i < M and j >= M) or (i >= M and j < M)
    ]

def _interval_width(theta_left: float, theta_right: float) -> float:
    width = float(theta_right) - float(theta_left)
    if width <= 0:
        width += 2 * pi
    return width

def _unwrap_theta(theta: float, theta_left: float) -> float:
    theta = float(theta)
    while theta <= theta_left:
        theta += 2 * pi
    return theta

def _get_boundary_cluster_bounds(
    roots_point: np.ndarray,
    M: int,
    GBZ_check_tol: float,
) -> tuple[int, int]:
    abs_roots = np.abs(roots_point)
    lower = abs_roots[M - 1]
    upper = abs_roots[M]
    scale = max(1.0, 0.5 * (lower + upper))
    eq_tol = GBZ_check_tol * scale

    candidate_inds = [
        ind for ind, val in enumerate(abs_roots)
        if (abs(val - lower) <= eq_tol) or (abs(val - upper) <= eq_tol)
    ]
    if not candidate_inds:
        return M - 1, M + 1

    start = min(candidate_inds)
    end = max(candidate_inds) + 1
    return start, end

def _analyze_boundary_matching(
    roots_left: np.ndarray,
    roots_right: np.ndarray,
    M: int,
    GBZ_check_tol: float,
    match_confidence_tol: float,
    double_root_tol: float,
) -> dict:
    def detect_boundary_double_root(roots_point: np.ndarray):
        start, end = _get_boundary_cluster_bounds(roots_point, M, GBZ_check_tol)
        if end - start < 2:
            return None

        cluster_roots = roots_point[start:end]
        pair_cost = chordal_cost_matrix(cluster_roots, cluster_roots)
        np.fill_diagonal(pair_cost, np.inf)
        min_ind = int(np.argmin(pair_cost))
        min_dist = float(pair_cost.flat[min_ind])
        if (not np.isfinite(min_dist)) or (min_dist > double_root_tol):
            return None

        i_local, j_local = np.unravel_index(min_ind, pair_cost.shape)
        if i_local > j_local:
            i_local, j_local = j_local, i_local
        beta_a = cluster_roots[i_local]
        beta_b = cluster_roots[j_local]
        if np.isfinite(beta_a.real) and np.isfinite(beta_a.imag) and np.isfinite(beta_b.real) and np.isfinite(beta_b.imag):
            beta_zero = 0.5 * (beta_a + beta_b)
        elif np.isfinite(beta_a.real) and np.isfinite(beta_a.imag):
            beta_zero = beta_a
        else:
            beta_zero = beta_b

        return {
            "inds": (start + int(i_local), start + int(j_local)),
            "distance": min_dist,
            "beta2_zero": beta_zero,
        }

    matches = hungarian_match_indices(roots_left, roots_right)
    cross_matches = [
        (i, j)
        for i, j in matches
        if (i < M and j >= M) or (i >= M and j < M)
    ]

    left_start, left_end = _get_boundary_cluster_bounds(roots_left, M, GBZ_check_tol)
    right_start, right_end = _get_boundary_cluster_bounds(roots_right, M, GBZ_check_tol)
    block_start = min(left_start, right_start)
    block_end = max(left_end, right_end)
    local_left = roots_left[block_start:block_end]
    local_right = roots_right[block_start:block_end]

    local_cost = chordal_cost_matrix(local_left, local_right)
    row_ind, col_ind = linear_sum_assignment(local_cost)
    local_match = {
        block_start + int(i): block_start + int(j)
        for i, j in zip(row_ind, col_ind)
    }
    matched_costs = local_cost[row_ind, col_ind]
    cost_scale = max(1e-12, float(np.mean(matched_costs)))

    exchange_margins = []
    left_rows = range(block_start, min(M, block_end))
    right_rows = range(max(M, block_start), block_end)
    for i in left_rows:
        if i not in local_match:
            continue
        j = local_match[i]
        for k in right_rows:
            if k not in local_match:
                continue
            l = local_match[k]

            curr_cross = int(j >= M) + int(l < M)
            swap_cross = int(l >= M) + int(j < M)
            if curr_cross == swap_cross:
                continue

            curr_cost = local_cost[i - block_start, j - block_start] + local_cost[k - block_start, l - block_start]
            swap_cost = local_cost[i - block_start, l - block_start] + local_cost[k - block_start, j - block_start]
            exchange_margins.append(float(swap_cost - curr_cost))

    if exchange_margins:
        min_exchange_margin = min(exchange_margins)
        confidence = min_exchange_margin / cost_scale
        is_confident = confidence > match_confidence_tol
    else:
        min_exchange_margin = np.inf
        confidence = np.inf
        is_confident = True

    return {
        "left_double_root": detect_boundary_double_root(roots_left),
        "right_double_root": detect_boundary_double_root(roots_right),
        "matches": matches,
        "cross_matches": cross_matches,
        "is_confident": is_confident,
        "confidence": confidence,
        "exchange_margin": min_exchange_margin,
        "cluster_block": (block_start, block_end),
    }

def _get_pmgbz_cluster_data(
    roots_point: np.ndarray,
    M: int,
    GBZ_check_tol: float,
) -> tuple[list[int], int, int]:
    start, end = _get_boundary_cluster_bounds(roots_point, M, GBZ_check_tol)
    candidate_inds = list(range(start, end))
    if not candidate_inds:
        raise ValueError("No PMGBZ-equivalent roots found at a degenerate point.")

    i_min = min(candidate_inds)
    if candidate_inds != list(range(i_min, i_min + len(candidate_inds))):
        raise ValueError(f"PMGBZ equal-modulus cluster is not contiguous: {candidate_inds}")

    n_pos = M - i_min
    if n_pos <= 0 or n_pos > len(candidate_inds):
        raise ValueError(
            f"Invalid PMGBZ classification window: M={M}, i_min={i_min}, "
            f"n_pos={n_pos}, n_candidates={len(candidate_inds)}."
        )
    return candidate_inds, i_min, n_pos

def _build_limit_order(
    roots_point: np.ndarray,
    cluster_inds: list[int],
    ordered_cluster_inds: list[int],
) -> np.ndarray:
    ordered_roots = np.array(roots_point, copy=True)
    start = cluster_inds[0]
    ordered_roots[start:start + len(cluster_inds)] = np.array(
        [roots_point[ind] for ind in ordered_cluster_inds],
        dtype=complex
    )
    return ordered_roots

def _validate_pmgbz_point(
    roots_point: np.ndarray,
    M: int,
    poly_diff: PolyDiffContext,
    E_ref: complex,
    beta1: complex,
    theta_point: float,
    theta_left: float,
    theta_right: float,
    solve_sorted_roots: callable,
    GBZ_check_tol: float,
):
    cluster_inds, _, n_pos = _get_pmgbz_cluster_data(roots_point, M, GBZ_check_tol)

    scored = []
    for ind in cluster_inds:
        dmu2_dtheta1 = poly_diff.eval_dmu2((E_ref, beta1, roots_point[ind]))[1]
        scored.append((ind, dmu2_dtheta1))

    scored_desc = sorted(scored, key=lambda item: item[1], reverse=True)
    derivative_scale = max(1.0, max(abs(val) for _, val in scored_desc))
    derivative_tol = 1e-10 * derivative_scale
    analytic_distinguishable = all(
        abs(scored_desc[k][1] - scored_desc[k + 1][1]) > derivative_tol
        for k in range(len(scored_desc) - 1)
    )

    if analytic_distinguishable:
        left_cluster_inds = [ind for ind, _ in scored_desc]
        right_cluster_inds = list(reversed(left_cluster_inds))
        is_fake = False
    else:
        delta = 0.1 * min(_interval_width(theta_left, theta_point), _interval_width(theta_point, theta_right))
        if delta <= 1e-12:
            raise ValueError(
                f"PMGBZ validation interval too small for fallback finite difference: "
                f"left={theta_left}, point={theta_point}, right={theta_right}"
            )

        roots_minus = solve_sorted_roots(theta_point - delta)
        roots_plus = solve_sorted_roots(theta_point + delta)
        point_to_minus = dict(hungarian_match_indices(roots_point, roots_minus))
        point_to_plus = dict(hungarian_match_indices(roots_point, roots_plus))
        left_cluster_inds = [
            ind for ind, _ in sorted(
                ((ind, point_to_minus[ind]) for ind in cluster_inds),
                key=lambda item: item[1]
            )
        ]
        right_cluster_inds = [
            ind for ind, _ in sorted(
                ((ind, point_to_plus[ind]) for ind in cluster_inds),
                key=lambda item: item[1]
            )
        ]
        is_fake = (left_cluster_inds == right_cluster_inds)

    left_limit_roots = _build_limit_order(roots_point, cluster_inds, left_cluster_inds)
    right_limit_roots = _build_limit_order(roots_point, cluster_inds, right_cluster_inds)
    pos_sols = tuple(roots_point[ind] for ind in left_cluster_inds[:n_pos])
    neg_sols = tuple(roots_point[ind] for ind in left_cluster_inds[n_pos:])

    return {
        "is_fake": is_fake,
        "left_limit_roots": left_limit_roots,
        "right_limit_roots": right_limit_roots,
        "beta2_sols": [pos_sols, neg_sols, tuple()],
    }


def get_roots_and_PMGBZ(
    poly_diff: PolyDiffContext,  # f(E, beta1, beta2)
    E_ref: complex,  # reference energy
    mu1: float,  # mu1 = log(|beta1|)
    N_points: int = 301,  # number of points around the loop |beta1| = exp(mu1)
    zero_tol: float = 1e-10,  # tolerance for the zero of the characteristic polynomial
    GBZ_check_tol: float = 1e-6,
    refine_continuum: bool = True,
) -> tuple[GBZResult, np.ndarray, np.ndarray, dict]:
    """Solve beta2 roots on a theta1 mesh and detect all PMGBZ points at (E_ref, mu1).

    Parameters:
        poly_diff: evaluation context of the characteristic polynomial
            f(E, beta1, beta2).
        E_ref: reference energy.
        mu1: log|beta1| defining the loop beta1 = exp(mu1 + i*theta1).
        N_points: number of uniform theta1 mesh points on [0, 2*pi).
        zero_tol: relative tolerance on the gap |beta2_{M+1}| - |beta2_M|
            for declaring a point PMGBZ-degenerate.
        GBZ_check_tol: relative tolerance for grouping roots into the
            equal-modulus boundary cluster during matching/validation.

    Returns:
        (gbz_result, theta1_arr, sols_arr, info):
        gbz_result: GBZResult whose subsets are PointSubset (accidental
            points / double roots) and LineSubset (continuum intervals),
            with index = (n_0D, n_1D).
        theta1_arr: (K+1,) sorted theta1 mesh including refined special
            points, closed periodically (last = first + 2*pi).
        sols_arr: (K+1, M+N) beta2 roots at each theta1, sorted by modulus.
        info: {"M", "N", "continuum_flag", "_pmgbz_raw"} — minor degrees,
            whether any continuum interval was found, and the raw PMGBZ
            point dicts consumed by the strip winding formula.
    """
    M, N = get_minor_degrees(poly_diff)
    theta_base = np.linspace(0.0, 2 * pi, N_points, endpoint=False)
    if theta_base.size == 0:
        raise ValueError("N_points must be a positive integer.")

    param_ind = pt.CIndexVec((0, 1))
    var_ind = pt.CIndexVec([2])

    def normalize_theta(theta: float) -> float:
        return float(theta % (2 * pi))

    # Shared cache across the entire PMGBZ pipeline: mesh scan, bisection
    # refinement, and final assembly all reuse the same solved roots.
    cache: dict = {}

    def solve_sorted_roots(theta1: float) -> np.ndarray:
        return _solve_sorted_roots_at(
            poly_diff, E_ref, mu1, theta1, M, N,
            param_ind, var_ind, cache,
        )

    def is_degenerate(theta1: float) -> bool:
        return _is_pmgbz_degenerate(solve_sorted_roots(theta1), M, zero_tol)

    def refine_continuum_boundary(
        theta_nondeg: float,
        theta_deg: float,
        deg_on_right: bool,
        max_iter: int = 60
    ) -> float:
        left = float(theta_nondeg)
        right = float(theta_deg)
        if right <= left:
            right += 2 * pi

        for _ in range(max_iter):
            mid = 0.5 * (left + right)
            mid_is_deg = is_degenerate(mid)
            if deg_on_right:
                if mid_is_deg:
                    right = mid
                else:
                    left = mid
            else:
                if mid_is_deg:
                    left = mid
                else:
                    right = mid
        return normalize_theta(0.5 * (left + right))

    def refine_accidental_point(
        theta_left: float,
        roots_left: np.ndarray,
        theta_right: float,
        max_iter: int = 60
    ) -> tuple[float, np.ndarray]:
        left = float(theta_left)
        right = float(theta_right)
        if right <= left:
            right += 2 * pi

        left_roots = roots_left

        for _ in range(max_iter):
            mid = 0.5 * (left + right)
            mid_roots = solve_sorted_roots(mid)

            if _is_pmgbz_degenerate(mid_roots, M, zero_tol):
                return normalize_theta(mid), mid_roots

            cross_left_mid = _cross_boundary_matches(left_roots, mid_roots, M)
            if cross_left_mid:
                right = mid
            else:
                left = mid
                left_roots = mid_roots

        theta_refined = normalize_theta(0.5 * (left + right))
        return theta_refined, solve_sorted_roots(theta_refined)

    # 1) Solve on a uniform theta1 mesh (fills cache) and sort by |beta2|.
    _, sols_mesh = solve_roots_on_mesh(
        poly_diff, E_ref, mu1, N_points, cache=cache,
    )
    roots_base = sols_mesh[:-1]

    # 2) Detect continuum-degenerate intervals from |beta_M| == |beta_{M+1}|.
    degenerate_mask = np.array(
        [_is_pmgbz_degenerate(roots_base[i], M, zero_tol) for i in range(theta_base.size)],
        dtype=bool
    )
    continuum_intervals = find_cyclic_true_intervals(degenerate_mask)

    PMGBZ_points = []
    special_thetas = []
    continuum_flag = False
    theta_step = 2 * pi / theta_base.size
    point_like_theta_tol = max(1e-8, 0.25 * theta_step)
    recursive_theta_tol = max(1e-12, 1e-4 * theta_step)
    match_confidence_tol = 1e-3
    double_root_tol = GBZ_check_tol

    def should_recurse_interval(
        theta_a: float,
        theta_b: float,
        parent_width: float,
    ) -> bool:
        child_width = _interval_width(theta_a, theta_b)
        progress_tol = max(recursive_theta_tol, 1e-8 * parent_width)
        return (child_width > progress_tol) and (child_width < parent_width - progress_tol)

    def process_interval(
        theta_left: float,
        roots_left: np.ndarray,
        theta_right: float,
        roots_right: np.ndarray,
        depth: int = 0,
    ) -> list[dict]:
        if depth > 600:
            raise RecursionError(
                f"Exceeded PMGBZ recursion depth on interval ({theta_left}, {theta_right})."
            )
        analysis = _analyze_boundary_matching(
            roots_left,
            roots_right,
            M,
            GBZ_check_tol,
            match_confidence_tol,
            double_root_tol,
        )

        double_root_points = []
        if not analysis["is_confident"]:
            for theta_point, root_info in (
                (theta_left, analysis["left_double_root"]),
                (theta_right, analysis["right_double_root"]),
            ):
                if root_info is None:
                    continue
                theta_norm = normalize_theta(theta_point)
                if any(
                    abs(((p["theta1"] - theta_norm + pi) % (2 * pi)) - pi) < 1e-8
                    for p in double_root_points
                ):
                    continue
                double_root_points.append(
                    {
                        "theta1": theta_norm,
                        "beta2_sols": [tuple(), tuple(), (root_info["beta2_zero"],)],
                        "is_continuum": False,
                    }
                )
        if double_root_points:
            return double_root_points

        if _interval_width(theta_left, theta_right) <= recursive_theta_tol:
            if analysis["is_confident"] and analysis["cross_matches"]:
                raise ValueError(
                    f"Unresolved PMGBZ crossings remain in a tiny interval: "
                    f"left={theta_left}, right={theta_right}"
                )
            return []

        if not analysis["is_confident"]:
            theta_mid = 0.5 * (theta_left + theta_right)
            roots_mid = solve_sorted_roots(theta_mid)
            left_points = process_interval(
                theta_left,
                roots_left,
                theta_mid,
                roots_mid,
                depth + 1,
            )
            right_points = process_interval(
                theta_mid,
                roots_mid,
                theta_right,
                roots_right,
                depth + 1,
            )
            return left_points + right_points

        if not analysis["cross_matches"]:
            return []

        parent_width = _interval_width(theta_left, theta_right)
        theta_cross, roots_cross = refine_accidental_point(theta_left, roots_left, theta_right)
        theta_cross = _unwrap_theta(theta_cross, theta_left)
        beta1_cross = exp(mu1 + 1j * normalize_theta(theta_cross))
        validation = _validate_pmgbz_point(
            roots_cross,
            M,
            poly_diff,
            E_ref,
            beta1_cross,
            theta_cross,
            theta_left,
            theta_right,
            solve_sorted_roots,
            GBZ_check_tol,
        )

        left_points = []
        if should_recurse_interval(theta_left, theta_cross, parent_width):
            left_points = process_interval(
                theta_left,
                roots_left,
                theta_cross,
                validation["left_limit_roots"],
                depth + 1,
            )

        right_points = []
        if should_recurse_interval(theta_cross, theta_right, parent_width):
            right_points = process_interval(
                theta_cross,
                validation["right_limit_roots"],
                theta_right,
                roots_right,
                depth + 1,
            )

        curr_points = []
        if not validation["is_fake"]:
            curr_points.append(
                {
                    "theta1": normalize_theta(theta_cross),
                    "beta2_sols": validation["beta2_sols"],
                    "is_continuum": False,
                }
            )

        return left_points + curr_points + right_points

    # --- sweep-mode pre-pass: skip refinement for runs of >=2 degenerate mesh points ---
    run_lengths = [
        (end_idx - start_idx) % theta_base.size + 1
        for start_idx, end_idx in continuum_intervals
    ]
    if not refine_continuum and any(L >= 2 for L in run_lengths):
        continuum_flag = True
        for start_idx, end_idx in continuum_intervals:
            L = (end_idx - start_idx) % theta_base.size + 1
            if L >= 2:
                start_theta = normalize_theta(theta_base[start_idx])
                end_theta = normalize_theta(theta_base[end_idx])
                PMGBZ_points.append({
                    "theta1_start": start_theta,
                    "theta1_end": end_theta,
                    "is_continuum": True,
                })
        # Assemble provisional gbz (non-empty, replaced later by precise solve).
        subsets_early = []
        for p in PMGBZ_points:
            if p.get("is_continuum", False):
                subsets_early.append(LineSubset(
                    E=E_ref, mu1=mu1,
                    theta1_start=p["theta1_start"],
                    theta1_end=p["theta1_end"],
                    _M=M, _N=N, _poly_diff=poly_diff,
                ))
            else:
                beta1_pt = exp(mu1 + 1j * p["theta1"])
                for beta2_val in chain(*p["beta2_sols"]):
                    subsets_early.append(PointSubset(E=E_ref, beta1=beta1_pt, beta2=beta2_val))
        n_0d_early = sum(1 for s in subsets_early if isinstance(s, PointSubset))
        n_1d_early = sum(1 for s in subsets_early if isinstance(s, LineSubset))
        gbz_result = GBZResult(
            E_ref=E_ref, subsets=subsets_early, index=(n_0d_early, n_1d_early),
        )
        theta1_arr = np.hstack((theta_base, [theta_base[0] + 2 * pi]))
        sols_arr = np.vstack((roots_base, roots_base[0]))
        return gbz_result, theta1_arr, sols_arr, {
            "M": M, "N": N, "continuum_flag": True,
            "_pmgbz_raw": PMGBZ_points,
        }

    for start_idx, end_idx in continuum_intervals:
        if np.all(degenerate_mask):
            start_theta = 0.0
            end_theta = 2 * pi
        else:
            prev_idx = (start_idx - 1) % theta_base.size
            next_idx = (end_idx + 1) % theta_base.size
            start_theta = refine_continuum_boundary(
                theta_base[prev_idx],
                theta_base[start_idx],
                deg_on_right=True
            )
            end_theta = refine_continuum_boundary(
                theta_base[end_idx],
                theta_base[next_idx],
                deg_on_right=False
            )

        if np.all(degenerate_mask):
            interval_width = 2 * pi
        else:
            interval_width = (end_theta - start_theta) % (2 * pi)

        # A "continuum interval" with vanishing width is just one accidental point
        # (including the periodic case start~2pi and end~0 for the same theta).
        if interval_width <= point_like_theta_tol:
            theta_point = normalize_theta(start_theta)
            roots_point = solve_sorted_roots(theta_point)
            beta1_point = exp(mu1 + 1j * theta_point)
            validation = _validate_pmgbz_point(
                roots_point,
                M,
                poly_diff,
                E_ref,
                beta1_point,
                theta_point,
                theta_point - theta_step,
                theta_point + theta_step,
                solve_sorted_roots,
                GBZ_check_tol,
            )
            if not validation["is_fake"]:
                PMGBZ_points.append(
                    {
                        "theta1": theta_point,
                        "beta2_sols": validation["beta2_sols"],
                        "is_continuum": False
                    }
                )
                special_thetas.append(theta_point)
            continue

        continuum_flag = True
        PMGBZ_points.append(
            {
                "theta1_start": start_theta,
                "theta1_end": end_theta,
                "is_continuum": True
            }
        )
        special_thetas.extend([start_theta, end_theta])

    # --- sweep-mode second early return: all run==1 intervals processed, some turned out continuum ---
    if not refine_continuum and continuum_flag:
        subsets_early2 = []
        for p in PMGBZ_points:
            if p.get("is_continuum", False):
                subsets_early2.append(LineSubset(
                    E=E_ref, mu1=mu1,
                    theta1_start=p["theta1_start"],
                    theta1_end=p["theta1_end"],
                    _M=M, _N=N, _poly_diff=poly_diff,
                ))
            else:
                beta1_pt = exp(mu1 + 1j * p["theta1"])
                for beta2_val in chain(*p["beta2_sols"]):
                    subsets_early2.append(PointSubset(E=E_ref, beta1=beta1_pt, beta2=beta2_val))
        n_0d_e2 = sum(1 for s in subsets_early2 if isinstance(s, PointSubset))
        n_1d_e2 = sum(1 for s in subsets_early2 if isinstance(s, LineSubset))
        gbz_result = GBZResult(
            E_ref=E_ref, subsets=subsets_early2, index=(n_0d_e2, n_1d_e2),
        )
        theta1_arr = np.hstack((theta_base, [theta_base[0] + 2 * pi]))
        sols_arr = np.vstack((roots_base, roots_base[0]))
        return gbz_result, theta1_arr, sols_arr, {
            "M": M, "N": N, "continuum_flag": True,
            "_pmgbz_raw": PMGBZ_points,
        }

    # 3) In non-continuum segments, use Hungarian matching to find accidental PMGBZ points.
    for i in range(theta_base.size):
        j = (i + 1) % theta_base.size
        if degenerate_mask[i] or degenerate_mask[j]:
            continue

        roots_i = roots_base[i]
        roots_j = roots_base[j]

        theta_left = float(theta_base[i])
        theta_right = float(theta_base[j])
        if theta_right <= theta_left:
            theta_right += 2 * pi
        analysis = _analyze_boundary_matching(
            roots_i,
            roots_j,
            M,
            GBZ_check_tol,
            match_confidence_tol,
            double_root_tol,
        )
        if analysis["is_confident"] and not analysis["cross_matches"]:
            continue

        interval_points = process_interval(
            theta_left, roots_i, theta_right, roots_j
        )
        for point in interval_points:
            theta_curr = point["theta1"]
            already_added = any(
                (not p.get("is_continuum", False))
                and abs(((p["theta1"] - theta_curr + pi) % (2 * pi)) - pi) < 1e-8
                for p in PMGBZ_points
            )
            if already_added:
                continue
            PMGBZ_points.append(point)
            special_thetas.append(theta_curr)

    # 4) Build output arrays with added special points and periodic closure.
    if special_thetas:
        theta_core = np.unique(np.hstack((theta_base, np.array(special_thetas))))
    else:
        theta_core = theta_base
    theta_core.sort()

    sols_core = np.vstack([solve_sorted_roots(theta) for theta in theta_core])
    theta1_arr = np.hstack((theta_core, [theta_core[0] + 2 * pi]))
    sols_arr = np.vstack((sols_core, sols_core[0]))

    # Build GBZResult with ConnectedSubsets
    subsets = []
    for p in PMGBZ_points:
        if p.get("is_continuum", False):
            subsets.append(LineSubset(
                E=E_ref, mu1=mu1,
                theta1_start=p["theta1_start"],
                theta1_end=p["theta1_end"],
                _M=M, _N=N, _poly_diff=poly_diff,
            ))
        else:
            beta1_pt = exp(mu1 + 1j * p["theta1"])
            pos_sols, neg_sols, zero_sols = p["beta2_sols"]
            for beta2_val in chain(pos_sols, neg_sols, zero_sols):
                subsets.append(PointSubset(E=E_ref, beta1=beta1_pt, beta2=beta2_val))

    n_0d = sum(1 for s in subsets if isinstance(s, PointSubset))
    n_1d = sum(1 for s in subsets if isinstance(s, LineSubset))
    gbz_result = GBZResult(E_ref=E_ref, subsets=subsets, index=(n_0d, n_1d))

    return gbz_result, theta1_arr, sols_arr, {
        "M": M, "N": N, "continuum_flag": continuum_flag,
        "_pmgbz_raw": PMGBZ_points,  # kept for internal winding formula
    }

