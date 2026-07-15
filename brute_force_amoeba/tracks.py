'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-14
Copyright © Department of Physics, Tsinghua University. All rights reserved

Root tracking via Hungarian matching for amoeba GBZ computation.
'''

import numpy as np
import poly_tools as pt
from math import pi
from cmath import exp

from brute_force_SGBZ.root_solver import calculate_point_roots
from brute_force_SGBZ.winding import PolyDiffContext
from gbz_types import (
    get_minor_degrees, sort_by_root_abs,
    chordal_cost_matrix, hungarian_match_indices,
)


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
    M, N = get_minor_degrees(PolyDiffContext(char_poly))
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
    tracked[0] = sort_by_root_abs(all_roots[0])

    for i in range(N_points - 1):
        roots_next = all_roots[i + 1]
        matches = hungarian_match_indices(tracked[i], roots_next)
        reordered = np.zeros(n_roots, dtype=complex)
        for from_idx, to_idx in matches:
            reordered[from_idx] = roots_next[to_idx]
        tracked[i + 1] = reordered

    return theta1_arr, tracked, M, N


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

    # Periodic extension respecting physical root identity
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
