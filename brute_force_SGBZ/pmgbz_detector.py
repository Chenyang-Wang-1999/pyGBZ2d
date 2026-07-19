'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2025-11-18
Copyright © Department of Physics, Tsinghua University. All rights reserved
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
from math import pi
from cmath import exp
from itertools import chain
from scipy.optimize import linear_sum_assignment

from .root_solver import _solve_sorted_roots_at, solve_roots_on_mesh

from gbz_types import (
    PointSubset, LineSubset, GBZResult, CharPoly,
    to_sphere_r3, cost_from_sphere_r3, hungarian_match_indices,
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

def _normalize_theta(theta: float) -> float:
    return float(theta % (2 * pi))

def _theta_close(theta_a: float, theta_b: float, tol: float = 1e-8) -> bool:
    """Cyclic closeness of two theta1 values (2*pi wrap-aware)."""
    return abs(((theta_a - theta_b + pi) % (2 * pi)) - pi) < tol

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
    poly: CharPoly,
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
        dmu2_dtheta1 = poly.eval_dmu2((E_ref, beta1, roots_point[ind]))[1]
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


# ---------------------------------------------------------------------------
# Structure: one entry point orchestrating _PmgbzScan
# ---------------------------------------------------------------------------
# LineSubset (continuum) detection is cheap (contiguous True runs of the
# degenerate mask + boundary bisection); PointSubset (accidental) detection
# is expensive (Hungarian matching + recursive splitting).  Two ingredients
# keep the single pipeline fast enough for the mu1-bisection hot path:
#
#   - refine_continuum=False (sweep mode) short-circuits as soon as a run of
#     >= 2 degenerate mesh points proves a continuum exists, skipping
#     boundary refinement and Stage 3 entirely — the mu1 bisection only
#     needs the continuum_flag, and the precise solve happens later via
#     refine_continuum=True in handle_continuum;
#   - the wide-gap pre-filter skips Hungarian matching on mesh pairs where
#     the |beta2| boundary gap stays large across consecutive mesh points:
#     a crossing needs the gap to reach zero, so it cannot hide inside such
#     a pair without the gap dipping visibly at a neighbouring sample.
# ---------------------------------------------------------------------------

class _PmgbzScan:
    """Shared per-call machinery for the PMGBZ pipeline.

    Holds the state every stage needs (polynomial, energy, mu1, root cache,
    mesh, derived tolerances) and provides the building blocks the entry
    point orchestrates: mesh solving, boundary-cluster matching analysis,
    continuum-boundary bisection, accidental-point bisection, the unified
    interval recursion, and output assembly.
    """

    def __init__(
        self,
        poly: CharPoly,
        E_ref: complex,
        mu1: float,
        N_points: int,
        zero_tol: float,
        GBZ_check_tol: float,
    ):
        theta_base = np.linspace(0.0, 2 * pi, N_points, endpoint=False)
        if theta_base.size == 0:
            raise ValueError("N_points must be a positive integer.")

        self.poly = poly
        self.E_ref = E_ref
        self.mu1 = mu1
        self.M = poly.M
        self.zero_tol = zero_tol
        self.GBZ_check_tol = GBZ_check_tol

        self.theta_base = theta_base
        self.theta_step = 2 * pi / theta_base.size
        self.point_like_theta_tol = max(1e-8, 0.25 * self.theta_step)
        self.recursive_theta_tol = max(1e-12, 1e-4 * self.theta_step)
        self.match_confidence_tol = 1e-3
        self.double_root_tol = GBZ_check_tol
        # Wide-gap pre-filter threshold (relative to the boundary modulus
        # scale).  Deliberately large: a mesh pair is skipped only when the
        # gap exceeds 50% of the root scale at the pair AND its two outer
        # neighbours, i.e. where a crossing (gap -> 0) cannot plausibly hide
        # between samples.  Replaces the Hungarian-only sweep pre-filter.
        self.gap_skip_tol = 0.5

        # Shared cache across the entire PMGBZ pipeline: mesh scan, bisection
        # refinement, and final assembly all reuse the same solved roots.
        self.cache: dict = {}

        # Filled by solve_mesh().
        self.roots_base: np.ndarray = None
        self.degenerate_mask: np.ndarray = None
        self.wide_gap_mask: np.ndarray = None
        self.continuum_intervals: list[tuple[int, int]] = None

    # ---- root solving ----

    def solve(self, theta1: float) -> np.ndarray:
        return _solve_sorted_roots_at(
            self.poly, self.E_ref, self.mu1, theta1, cache=self.cache,
        )

    def is_degenerate(self, theta1: float) -> bool:
        return _is_pmgbz_degenerate(self.solve(theta1), self.M, self.zero_tol)

    def solve_mesh(self) -> None:
        """Solve on the uniform theta1 mesh (fills cache), sort by |beta2|,
        and derive the degeneracy / wide-gap masks from the boundary gap
        |beta_M| - |beta_{M-1}| at every mesh point."""
        _, sols_mesh = solve_roots_on_mesh(
            self.poly, self.E_ref, self.mu1, self.theta_base.size, cache=self.cache,
        )
        self.roots_base = sols_mesh[:-1]

        abs_roots = np.abs(self.roots_base)
        lower = abs_roots[:, self.M - 1]
        upper = abs_roots[:, self.M]
        gaps = upper - lower
        scales = np.maximum(1.0, 0.5 * (lower + upper))
        self.degenerate_mask = gaps <= self.zero_tol * scales
        self.wide_gap_mask = gaps > self.gap_skip_tol * scales
        self.continuum_intervals = find_cyclic_true_intervals(self.degenerate_mask)

    # ---- boundary-cluster matching analysis ----

    def analyze_boundary_matching(
        self,
        roots_left: np.ndarray,
        roots_right: np.ndarray,
    ) -> dict:
        """Hungarian matching + confidence + double-root detection for one interval.

        The sphere projections and the full cost matrix are computed once and
        reused everywhere: the global Hungarian match, the local boundary
        block (a plain slice of the full matrix), and the double-root
        self-cost matrices (built from slices of the cached projections).
        """
        M = self.M
        r3_left = to_sphere_r3(roots_left)
        r3_right = to_sphere_r3(roots_right)
        cost = cost_from_sphere_r3(r3_left, r3_right)
        if np.any(np.isnan(cost)):
            raise ValueError(
                f"NaN cost in boundary matching: roots_left={roots_left}, "
                f"roots_right={roots_right}"
            )

        row_ind, col_ind = linear_sum_assignment(cost)
        cross_matches = [
            (int(i), int(j))
            for i, j in zip(row_ind, col_ind)
            if (i < M and j >= M) or (i >= M and j < M)
        ]

        left_bounds = _get_boundary_cluster_bounds(roots_left, M, self.GBZ_check_tol)
        right_bounds = _get_boundary_cluster_bounds(roots_right, M, self.GBZ_check_tol)
        block_start = min(left_bounds[0], right_bounds[0])
        block_end = max(left_bounds[1], right_bounds[1])

        local_cost = cost[block_start:block_end, block_start:block_end]
        local_row, local_col = linear_sum_assignment(local_cost)
        local_match = {
            block_start + int(i): block_start + int(j)
            for i, j in zip(local_row, local_col)
        }
        matched_costs = local_cost[local_row, local_col]
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
            confidence = min(exchange_margins) / cost_scale
            is_confident = confidence > self.match_confidence_tol
        else:
            is_confident = True

        def detect_boundary_double_root(roots_point, r3_point, bounds):
            start, end = bounds
            if end - start < 2:
                return None

            pair_cost = cost_from_sphere_r3(r3_point[start:end], r3_point[start:end])
            np.fill_diagonal(pair_cost, np.inf)
            min_ind = int(np.argmin(pair_cost))
            min_dist = float(pair_cost.flat[min_ind])
            if (not np.isfinite(min_dist)) or (min_dist > self.double_root_tol):
                return None

            i_local, j_local = np.unravel_index(min_ind, pair_cost.shape)
            if i_local > j_local:
                i_local, j_local = j_local, i_local
            beta_a = roots_point[start + i_local]
            beta_b = roots_point[start + j_local]
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

        return {
            "left_double_root": detect_boundary_double_root(roots_left, r3_left, left_bounds),
            "right_double_root": detect_boundary_double_root(roots_right, r3_right, right_bounds),
            "cross_matches": cross_matches,
            "is_confident": is_confident,
        }

    # ---- bisection refinement ----

    def refine_continuum_boundary(
        self,
        theta_nondeg: float,
        theta_deg: float,
        deg_on_right: bool,
        max_iter: int = 60,
    ) -> float:
        left = float(theta_nondeg)
        right = float(theta_deg)
        if right <= left:
            right += 2 * pi

        for _ in range(max_iter):
            mid = 0.5 * (left + right)
            mid_is_deg = self.is_degenerate(mid)
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
        return _normalize_theta(0.5 * (left + right))

    def refine_accidental_point(
        self,
        theta_left: float,
        roots_left: np.ndarray,
        theta_right: float,
        max_iter: int = 60,
    ) -> tuple[float, np.ndarray]:
        left = float(theta_left)
        right = float(theta_right)
        if right <= left:
            right += 2 * pi

        left_roots = roots_left

        for _ in range(max_iter):
            mid = 0.5 * (left + right)
            mid_roots = self.solve(mid)

            if _is_pmgbz_degenerate(mid_roots, self.M, self.zero_tol):
                return _normalize_theta(mid), mid_roots

            if _cross_boundary_matches(left_roots, mid_roots, self.M):
                right = mid
            else:
                left = mid
                left_roots = mid_roots

        theta_refined = _normalize_theta(0.5 * (left + right))
        return theta_refined, self.solve(theta_refined)

    # ---- accidental-point recursion ----

    def validate_point(
        self,
        roots_point: np.ndarray,
        theta_point: float,
        theta_left: float,
        theta_right: float,
    ) -> dict:
        beta1_point = exp(self.mu1 + 1j * _normalize_theta(theta_point))
        return _validate_pmgbz_point(
            roots_point,
            self.M,
            self.poly,
            self.E_ref,
            beta1_point,
            theta_point,
            theta_left,
            theta_right,
            self.solve,
            self.GBZ_check_tol,
        )

    def should_recurse_interval(
        self,
        theta_a: float,
        theta_b: float,
        parent_width: float,
    ) -> bool:
        child_width = _interval_width(theta_a, theta_b)
        progress_tol = max(self.recursive_theta_tol, 1e-8 * parent_width)
        return (child_width > progress_tol) and (child_width < parent_width - progress_tol)

    def process_interval(
        self,
        theta_left: float,
        roots_left: np.ndarray,
        theta_right: float,
        roots_right: np.ndarray,
        depth: int = 0,
    ) -> list[dict]:
        """Unified interval recursion (see log-old.md):

        - matching not confident → emit boundary double roots if present,
          otherwise subdivide;
        - confident, no cross-boundary match → interval done;
        - confident with crossings → one PMGBZ search, then recurse on both
          sides using the limit-ordered roots at the crossing.
        """
        if depth > 600:
            raise RecursionError(
                f"Exceeded PMGBZ recursion depth on interval ({theta_left}, {theta_right})."
            )
        analysis = self.analyze_boundary_matching(roots_left, roots_right)

        double_root_points = []
        if not analysis["is_confident"]:
            for theta_point, root_info in (
                (theta_left, analysis["left_double_root"]),
                (theta_right, analysis["right_double_root"]),
            ):
                if root_info is None:
                    continue
                theta_norm = _normalize_theta(theta_point)
                if any(_theta_close(p["theta1"], theta_norm) for p in double_root_points):
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

        if _interval_width(theta_left, theta_right) <= self.recursive_theta_tol:
            if analysis["is_confident"] and analysis["cross_matches"]:
                raise ValueError(
                    f"Unresolved PMGBZ crossings remain in a tiny interval: "
                    f"left={theta_left}, right={theta_right}"
                )
            return []

        if not analysis["is_confident"]:
            theta_mid = 0.5 * (theta_left + theta_right)
            roots_mid = self.solve(theta_mid)
            return (
                self.process_interval(theta_left, roots_left, theta_mid, roots_mid, depth + 1)
                + self.process_interval(theta_mid, roots_mid, theta_right, roots_right, depth + 1)
            )

        if not analysis["cross_matches"]:
            return []

        parent_width = _interval_width(theta_left, theta_right)
        theta_cross, roots_cross = self.refine_accidental_point(theta_left, roots_left, theta_right)
        theta_cross = _unwrap_theta(theta_cross, theta_left)
        validation = self.validate_point(roots_cross, theta_cross, theta_left, theta_right)

        left_points = []
        if self.should_recurse_interval(theta_left, theta_cross, parent_width):
            left_points = self.process_interval(
                theta_left,
                roots_left,
                theta_cross,
                validation["left_limit_roots"],
                depth + 1,
            )

        right_points = []
        if self.should_recurse_interval(theta_cross, theta_right, parent_width):
            right_points = self.process_interval(
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
                    "theta1": _normalize_theta(theta_cross),
                    "beta2_sols": validation["beta2_sols"],
                    "is_continuum": False,
                }
            )

        return left_points + curr_points + right_points

    # ---- per-stage drivers ----

    def handle_degenerate_interval(
        self,
        start_idx: int,
        end_idx: int,
        PMGBZ_points: list,
        special_thetas: list,
    ) -> None:
        """Refine one degenerate mesh run into a continuum interval or, when
        its refined width vanishes, a single accidental point."""
        theta_base = self.theta_base
        if np.all(self.degenerate_mask):
            # Whole loop degenerate.  (Only reachable with
            # refine_continuum=True: a full-circle run has length >= 2 and
            # triggers the sweep-mode early return.)
            start_theta = 0.0
            end_theta = 2 * pi
            interval_width = 2 * pi
        else:
            prev_idx = (start_idx - 1) % theta_base.size
            next_idx = (end_idx + 1) % theta_base.size
            start_theta = self.refine_continuum_boundary(
                theta_base[prev_idx],
                theta_base[start_idx],
                deg_on_right=True,
            )
            end_theta = self.refine_continuum_boundary(
                theta_base[end_idx],
                theta_base[next_idx],
                deg_on_right=False,
            )
            interval_width = (end_theta - start_theta) % (2 * pi)

        # A "continuum interval" with vanishing width is just one accidental point
        # (including the periodic case start~2pi and end~0 for the same theta).
        if interval_width <= self.point_like_theta_tol:
            theta_point = _normalize_theta(start_theta)
            roots_point = self.solve(theta_point)
            validation = self.validate_point(
                roots_point,
                theta_point,
                theta_point - self.theta_step,
                theta_point + self.theta_step,
            )
            if not validation["is_fake"]:
                PMGBZ_points.append(
                    {
                        "theta1": theta_point,
                        "beta2_sols": validation["beta2_sols"],
                        "is_continuum": False,
                    }
                )
                special_thetas.append(theta_point)
            return

        PMGBZ_points.append(
            {
                "theta1_start": start_theta,
                "theta1_end": end_theta,
                "is_continuum": True,
            }
        )
        special_thetas.extend([start_theta, end_theta])

    def accidental_scan(
        self,
        PMGBZ_points: list,
        special_thetas: list,
    ) -> None:
        """Stage 3: accidental PMGBZ points between non-degenerate mesh pairs.

        Wide-gap pre-filter: a pair is skipped when the boundary gap exceeds
        gap_skip_tol at the pair AND at both outer neighbours (4 consecutive
        mesh points) — a crossing needs |beta_M| == |beta_{M-1}| somewhere
        inside, which cannot happen without the gap collapsing at a nearby
        sample.  All other pairs go through the full interval recursion.
        """
        theta_base = self.theta_base
        n = theta_base.size
        wide = self.wide_gap_mask
        for i in range(n):
            j = (i + 1) % n
            if self.degenerate_mask[i] or self.degenerate_mask[j]:
                continue
            if wide[(i - 1) % n] and wide[i] and wide[j] and wide[(j + 1) % n]:
                continue

            roots_i = self.roots_base[i]
            roots_j = self.roots_base[j]

            theta_left = float(theta_base[i])
            theta_right = float(theta_base[j])
            if theta_right <= theta_left:
                theta_right += 2 * pi

            interval_points = self.process_interval(theta_left, roots_i, theta_right, roots_j)
            for point in interval_points:
                theta_curr = point["theta1"]
                already_added = any(
                    (not p.get("is_continuum", False)) and _theta_close(p["theta1"], theta_curr)
                    for p in PMGBZ_points
                )
                if already_added:
                    continue
                PMGBZ_points.append(point)
                special_thetas.append(theta_curr)

    def assemble(
        self,
        PMGBZ_points: list,
        special_thetas: list,
    ) -> tuple[GBZResult, np.ndarray, np.ndarray, dict]:
        """Stage 4: output arrays with added special points and periodic
        closure, plus the GBZResult of ConnectedSubsets."""
        if special_thetas:
            theta_core = np.unique(np.hstack((self.theta_base, np.array(special_thetas))))
        else:
            theta_core = self.theta_base
        theta_core.sort()

        sols_core = np.vstack([self.solve(theta) for theta in theta_core])
        theta1_arr = np.hstack((theta_core, [theta_core[0] + 2 * pi]))
        sols_arr = np.vstack((sols_core, sols_core[0]))

        subsets = []
        for p in PMGBZ_points:
            if p.get("is_continuum", False):
                subsets.append(LineSubset(
                    E=self.E_ref, mu1=self.mu1,
                    theta1_start=p["theta1_start"],
                    theta1_end=p["theta1_end"],
                ))
            else:
                beta1_pt = exp(self.mu1 + 1j * p["theta1"])
                pos_sols, neg_sols, zero_sols = p["beta2_sols"]
                for beta2_val in chain(pos_sols, neg_sols, zero_sols):
                    subsets.append(PointSubset(E=self.E_ref, beta1=beta1_pt, beta2=beta2_val))

        n_0d = sum(1 for s in subsets if isinstance(s, PointSubset))
        n_1d = sum(1 for s in subsets if isinstance(s, LineSubset))
        gbz_result = GBZResult(E_ref=self.E_ref, subsets=subsets, index=(n_0d, n_1d))

        return gbz_result, theta1_arr, sols_arr, {
            "continuum_flag": any(p.get("is_continuum", False) for p in PMGBZ_points),
            "_pmgbz_raw": PMGBZ_points,  # kept for internal winding formula
        }


def get_roots_and_PMGBZ(
    poly: CharPoly,  # f(E, beta1, beta2)
    E_ref: complex,  # reference energy
    mu1: float,  # mu1 = log(|beta1|)
    N_points: int = 301,  # number of points around the loop |beta1| = exp(mu1)
    zero_tol: float = 1e-10,  # tolerance for the zero of the characteristic polynomial
    GBZ_check_tol: float = 1e-6,
    refine_continuum: bool = True,
) -> tuple[GBZResult, np.ndarray, np.ndarray, dict]:
    """Find all PMGBZ points (LineSubset + PointSubset) at (E_ref, mu1).

    Pipeline: mesh root-solve → continuum-degenerate runs (boundaries
    refined by bisection; vanishing-width runs degrade to accidental
    points) → Stage 3 accidental-point recursion on non-degenerate mesh
    pairs, gated by the wide-gap pre-filter (see accidental_scan).

    refine_continuum=False (sweep mode, used by the mu1 bisection in
    solve_SGBZ_for_E): as soon as a run of >= 2 consecutive degenerate mesh
    points proves a continuum exists, return immediately.  The returned
    GBZResult is then a bare carrier for the continuum_flag=True info dict —
    handle_continuum computes the precise subsets via refine_continuum=True,
    so the placeholder is never exposed to users.

    Parameters / returns: see module docstring.
    """
    scan = _PmgbzScan(poly, E_ref, mu1, N_points, zero_tol, GBZ_check_tol)
    scan.solve_mesh()

    if not refine_continuum:
        for start_idx, end_idx in scan.continuum_intervals:
            run_length = (end_idx - start_idx) % scan.theta_base.size + 1
            if run_length >= 2:
                theta1_arr = np.hstack((scan.theta_base, [scan.theta_base[0] + 2 * pi]))
                sols_arr = np.vstack((scan.roots_base, scan.roots_base[0]))
                return GBZResult(E_ref=E_ref), theta1_arr, sols_arr, {
                    "continuum_flag": True,
                    "_pmgbz_raw": [{
                        "theta1_start": _normalize_theta(scan.theta_base[start_idx]),
                        "theta1_end": _normalize_theta(scan.theta_base[end_idx]),
                        "is_continuum": True,
                    }],
                }

    PMGBZ_points: list[dict] = []
    special_thetas: list[float] = []

    for start_idx, end_idx in scan.continuum_intervals:
        scan.handle_degenerate_interval(start_idx, end_idx, PMGBZ_points, special_thetas)

    scan.accidental_scan(PMGBZ_points, special_thetas)

    return scan.assemble(PMGBZ_points, special_thetas)
