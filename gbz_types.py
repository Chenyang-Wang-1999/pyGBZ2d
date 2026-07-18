'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-14
Copyright © Department of Physics, Tsinghua University. All rights reserved

Unified data types and shared utility functions for GBZ computation.

Provides:
  - ConnectedSubset abstraction: PointSubset (0D) and LineSubset (1D)
  - GBZResult container for per-E_ref computation results
  - Shared tool functions: chordal distance, Hungarian matching,
    cyclic interval detection, minor degree extraction, probe step generation
'''

from __future__ import annotations

import math
from cmath import exp, log
from dataclasses import dataclass
from itertools import chain
from typing import Any, Optional, TYPE_CHECKING, Union

import numpy as np
import poly_tools as pt
from scipy.optimize import linear_sum_assignment

if TYPE_CHECKING:
    from brute_force_SGBZ.winding import PolyDiffContext


# ---- data classes ----

@dataclass(frozen=True)
class PointSubset:
    """0D connected subset: an isolated GBZ point.

    Represented by the (E, beta1, beta2) triplet.  Frozen (immutable)
    because it is pure data with no lazy-loaded fields.

    Attributes:
        E: Reference energy.
        beta1: Bloch factor in the principal direction.
        beta2: Representative Bloch factor in the secondary direction.
    """
    E: complex
    beta1: complex
    beta2: complex

    @property
    def mu1(self) -> float:
        """|beta1| = exp(mu1)."""
        return float(math.log(abs(self.beta1)))

    @property
    def theta1(self) -> float:
        """Phase angle of beta1 in [0, 2π)."""
        return float(math.atan2(self.beta1.imag, self.beta1.real) % (2 * math.pi))

    def as_triplet(self) -> tuple[complex, complex, complex]:
        """Return (E, k1, k2) where k_j = -i * log(beta_j)."""
        k1 = -1j * log(self.beta1)
        k2 = -1j * log(self.beta2)
        return (self.E, k1, k2)


@dataclass
class LineSubset:
    """1D connected subset: a continuous degenerate line segment.

    beta2_arr is lazy-loaded via :meth:`fill_beta2`.  Call it before
    accessing :attr:`left_endpoint` or :attr:`right_endpoint`.

    Attributes:
        E: Reference energy.
        mu1: Fixed |beta1| radius (= ln|beta1|) across the segment.
        theta1_start: Start phase angle in [0, 2π).
        theta1_end: End phase angle in [0, 2π).
        beta2_arr: (N, 2) array of boundary beta2 roots at each theta1
                   sample point.  ``None`` until :meth:`fill_beta2` is called.
    """
    E: complex
    mu1: float
    theta1_start: float
    theta1_end: float
    beta2_arr: Optional[np.ndarray] = None  # shape (N, 2), lazy

    # Internal storage for lazy fill — not part of the public API.
    _M: int = 0
    _N: int = 0
    _poly_diff: Any = None  # PolyDiffContext

    @property
    def theta1_width(self) -> float:
        """Angular width of the interval (radians), handling 2π wrap."""
        w = self.theta1_end - self.theta1_start
        return float(w if w > 0 else w + 2 * math.pi)

    @property
    def left_endpoint(self) -> tuple[complex, complex]:
        """(beta1_left, beta2_left).  Requires beta2_arr to be filled."""
        if self.beta2_arr is None:
            raise RuntimeError("beta2_arr not filled; call fill_beta2() first")
        b1 = exp(self.mu1 + 1j * self.theta1_start)
        return (b1, self.beta2_arr[0, 0])

    @property
    def right_endpoint(self) -> tuple[complex, complex]:
        """(beta1_right, beta2_right).  Requires beta2_arr to be filled."""
        if self.beta2_arr is None:
            raise RuntimeError("beta2_arr not filled; call fill_beta2() first")
        b1 = exp(self.mu1 + 1j * self.theta1_end)
        return (b1, self.beta2_arr[-1, 1])

    def is_loaded(self) -> bool:
        """True if beta2_arr has been filled."""
        return self.beta2_arr is not None

    def fill_beta2(self, N_points: int = 301) -> None:
        """Lazy-load beta2_arr by solving roots on a uniform theta1 mesh.

        Uses a lightweight roots-only solver (no PMGBZ detection) that
        explicitly includes the interval endpoints, producing a result
        identical to the old get_roots_and_PMGBZ-based path.

        Uses a lazy import to avoid a module-level circular dependency
        between gbz_types and brute_force_SGBZ.pmgbz_detector.
        """
        if self.beta2_arr is not None:
            return
        # Lazy import — breaks circular dependency at module level.
        from brute_force_SGBZ.root_solver import solve_roots_on_mesh  # noqa: E402

        theta1_arr, sols_arr = solve_roots_on_mesh(
            self._poly_diff, self.E, self.mu1, N_points,
            extra_thetas=(self.theta1_start, self.theta1_end),
        )

        # Build a boolean mask for theta1 values inside the interval.
        t_start = self.theta1_start % (2 * math.pi)
        t_end = self.theta1_end % (2 * math.pi)
        if t_end > t_start:
            mask = (theta1_arr[:-1] >= t_start) & (theta1_arr[:-1] <= t_end)
        else:
            mask = (theta1_arr[:-1] >= t_start) | (theta1_arr[:-1] <= t_end)

        # Store the two boundary roots (M-1 and M) for each theta1.
        M = self._M
        self.beta2_arr = np.column_stack([
            sols_arr[:-1, M - 1][mask],
            sols_arr[:-1, M][mask],
        ])


@dataclass
class GBZResult:
    """Per-E_ref GBZ computation result.

    Attributes:
        E_ref: Reference energy.
        success: Whether the computation completed without error.
        error: Error message if ``success`` is False.
        subsets: List of connected subsets (PointSubset / LineSubset).
        index: (n_0D, n_1D) counts.  ``(0, 0)`` means the energy is
               outside the GBZ.
    """
    E_ref: complex
    success: bool = True
    error: Optional[str] = None
    subsets: list[ConnectedSubset] = None  # type: ignore[assignment]
    index: tuple[int, int] = (0, 0)

    def __post_init__(self):
        if self.subsets is None:
            object.__setattr__(self, 'subsets', [])

    @property
    def is_empty(self) -> bool:
        return len(self.subsets) == 0

    @property
    def is_gbz(self) -> bool:
        """True if this energy point lies on the GBZ."""
        return self.success and self.index != (0, 0)


ConnectedSubset = Union[PointSubset, LineSubset]


# ---- shared utility functions ----

def sort_by_root_abs(roots: np.ndarray) -> np.ndarray:
    """Sort complex roots by absolute value (modulus), ascending."""
    return np.array(sorted(roots, key=lambda x: abs(x)), dtype=complex)


def to_sphere_r3(roots: np.ndarray) -> np.ndarray:
    """Map complex roots to Riemann sphere R^3 coordinates.

    Finite roots map to the unit sphere via stereographic projection.
    Roots at infinity map to the north pole (0, 0, 1).

    Parameters:
        roots: 1-D array of complex numbers.

    Returns:
        Array of shape (len(roots), 3) with columns (x, y, z).
    """
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


def chordal_cost_matrix(arr1: np.ndarray, arr2: np.ndarray) -> np.ndarray:
    """Chordal-distance cost matrix between two root arrays.

    The chordal distance is the Euclidean distance on the Riemann sphere,
    computed via :func:`to_sphere_r3`.

    Parameters:
        arr1, arr2: 1-D arrays of complex roots.

    Returns:
        2-D array ``cost[i, j]`` = chordal distance between arr1[i] and arr2[j].
    """
    p1 = to_sphere_r3(arr1)
    p2 = to_sphere_r3(arr2)
    diff = p1[:, None, :] - p2[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=2), dtype=float)


def hungarian_match_indices(
    roots_from: np.ndarray,
    roots_to: np.ndarray,
) -> list[tuple[int, int]]:
    """Hungarian (linear-sum-assignment) matching between two root sets.

    Uses chordal distance as the cost metric.  Raises ``ValueError`` if
    either array contains NaN entries.

    Returns:
        List of ``(from_idx, to_idx)`` pairs.
    """
    if np.any(np.isnan(roots_from.real)) or np.any(np.isnan(roots_from.imag)):
        raise ValueError(f"NaN root encountered in roots_from: {roots_from}")
    if np.any(np.isnan(roots_to.real)) or np.any(np.isnan(roots_to.imag)):
        raise ValueError(f"NaN root encountered in roots_to: {roots_to}")

    cost = chordal_cost_matrix(roots_from, roots_to)
    row_ind, col_ind = linear_sum_assignment(cost)
    return [(int(i), int(j)) for i, j in zip(row_ind, col_ind)]


def find_cyclic_true_intervals(mask: np.ndarray) -> list[tuple[int, int]]:
    """Find contiguous ``True`` intervals in a cyclic boolean array.

    Returns:
        List of ``(start_idx, end_idx_inclusive)``.
        Empty list if no ``True`` values.
        ``[(0, n-1)]`` if all ``True``.
    """
    n = len(mask)
    if n == 0:
        return []
    if np.all(mask):
        return [(0, n - 1)]
    if not np.any(mask):
        return []

    starts = np.where(mask & (~np.roll(mask, 1)))[0]
    intervals: list[tuple[int, int]] = []
    for start in starts:
        end = start
        while mask[end]:
            end = (end + 1) % n
            if end == start:
                break
        intervals.append((int(start), int((end - 1) % n)))
    return intervals


def get_minor_degrees(
    poly_diff: PolyDiffContext,
    direction: int = 2,
) -> tuple[int, int]:
    """Extract (M, N) denominator / numerator minor degrees.

    For a CLaurent polynomial in 3 variables (E, beta1, beta2):
      - M = denominator order in the ``direction`` variable.
      - N = (max numerator degree) - M.

    Parameters:
        poly_diff: Polynomial context wrapping the characteristic polynomial.
        direction: Variable index (1 = beta1, 2 = beta2).  Default is 2.

    Returns:
        (M, N) tuple.
    """
    coeffs = pt.CScalarVec([])
    degs = pt.CIndexVec([])
    poly_diff.char_poly.num.batch_get_data(coeffs, degs)
    M_plus_N = max(degs[direction::3])
    M = poly_diff.char_poly.denom_orders[direction]
    N = M_plus_N - M
    return M, N


def generate_probe_steps(
    bracket_width: float,
    probe_radius: float,
    zero_tol: float,
    xtol: float = 1e-10,
) -> list[float]:
    """Generate probe step distances for zero-plateau detection.

    Shared by both SGBZ and Amoeba :func:`_probe_zero_plateau_near_mu1`
    implementations.  Returns a sorted list of positive step sizes
    ranging from ~zero_tol up to max(probe_radius, 4*bracket_width).

    Parameters:
        bracket_width: Width of the current bisection bracket in mu1.
        probe_radius: User-specified maximum probe distance.
        zero_tol: Tolerance for zero detection.
        xtol: Bisection tolerance (used to floor the minimum step).

    Returns:
        Sorted list of positive step sizes.
    """
    eps_min = max(10.0 * zero_tol, 1e-12)
    eps_max = max(float(probe_radius), 4.0 * bracket_width, 100.0 * zero_tol, eps_min)

    steps: set[float] = {eps_min, eps_max}
    if bracket_width > 0:
        steps.update([
            0.25 * bracket_width,
            0.5 * bracket_width,
            bracket_width,
            2.0 * bracket_width,
            4.0 * bracket_width,
        ])
    if probe_radius > 0:
        steps.update([
            0.25 * probe_radius,
            0.5 * probe_radius,
            float(probe_radius),
        ])

    step = eps_min
    while step < eps_max:
        steps.add(step)
        step *= 2.0

    return sorted(s for s in steps if s > 0 and s <= eps_max * (1 + 1e-12))
