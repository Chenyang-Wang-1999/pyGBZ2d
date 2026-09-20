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
  - ``_find_exact_crossing``         bracketed refinement with persistent ZM rows.
  - ``_get_average_winding_from_zeros``  average winding from a zero partition.
"""

from __future__ import annotations

from typing import TYPE_CHECKING
import numpy as np
from cmath import exp
from scipy.optimize import brentq

# Angular accuracy of the bracketed solve, independent of polynomial scaling.
CROSSING_XTOL: float = 1e-12
CROSSING_MAXITER: int = 500
from ..core import CharPoly, TWO_PI, live_defaults

if TYPE_CHECKING:
    from ..continuation.zero_manager import ZeroManager

@live_defaults(xtol="amoeba.ronkin_winding:CROSSING_XTOL", max_iter="amoeba.ronkin_winding:CROSSING_MAXITER")
def _find_exact_crossing(
    zm: ZeroManager,
    mu2: float,
    seg_idx: int,
    col_idx: int,
    theta_lo: float,
    theta_hi: float,
    *,
    xtol: float = None,
    max_iter: int = None,
) -> tuple[float, float]:
    """Solve ln|β₂_j(θ₁)| = μ₂ and retain every new root row in ``zm``.

    Interpolation only identifies the track when solving the polynomial;
    all function values come from actual roots. Endpoints are read from
    the mesh, including θ₁=2π, which ``solve_at`` cannot accept as an
    interior point. The returned θ₁ stays unwrapped for mesh bookkeeping.
    A collapsed bracket returns an existing exact mesh crossing. Brent's
    angular convergence is the stopping criterion; a separate residual gate
    would impose a different accuracy requirement on steep root tracks.
    """
    context = (f"E_ref={zm.E_ref}, mu1={zm.mu1}, mu2={mu2}, "
               f"segment={seg_idx}, track={col_idx}, "
               f"bracket=({theta_lo:.17g}, {theta_hi:.17g})")

    def root_at(theta):
        seg = zm.segments[seg_idx]
        row = int(np.searchsorted(seg.theta1_arr, theta))
        if row == len(seg.theta1_arr) or seg.theta1_arr[row] != theta:
            # Re-locate on every call: earlier evaluations have inserted rows.
            row, _ = zm.insert_solution(theta, seg_idx=seg_idx, i=row - 1)
        root = complex(seg.tracked_roots[row, col_idx])
        if not np.isfinite(root) or abs(root) == 0:
            raise ValueError(f"non-finite log-modulus at theta1={theta:.17g}")
        return root

    def func(theta):
        root = root_at(float(theta))
        return float(np.log(np.abs(root)) - mu2)

    try:
        if theta_lo == theta_hi:
            t1 = float(theta_lo)
        else:
            t1 = float(brentq(func, theta_lo, theta_hi,
                             xtol=xtol,
                             rtol=4.0 * np.finfo(float).eps,
                             maxiter=max_iter))
        t2 = float(np.angle(root_at(t1)) % TWO_PI)
        return t1, t2
    except (ValueError, RuntimeError, FloatingPointError) as exc:
        raise RuntimeError(f"crossing refinement failed: {context}; {exc}") from exc


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
