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

# ---- characteristic polynomial wrapper ----

class CharPoly:
    """Characteristic Laurent polynomial f(E, beta1, beta2).

    The single entry point for polynomial construction, evaluation, and
    root-solving.  Wraps poly_tools.CLaurent internally — no other file
    in the project needs to import poly_tools directly.

    Parameters:
        coeffs: 1-D complex ndarray of polynomial coefficients.
        degs: (n_terms, 3) integer ndarray of (E, beta1, beta2) exponents.
    """

    def __init__(self, coeffs: np.ndarray, degs: np.ndarray):
        self._coeffs = np.asarray(coeffs, dtype=complex)
        self._degs = np.asarray(degs, dtype=int)
        self._claurent = pt.CLaurent(3)
        self._claurent.set_Laurent_by_terms(
            pt.CScalarVec(self._coeffs),
            pt.CLaurentIndexVec(self._degs.flatten()),
        )
        # Pre-compute partial derivatives (absorbs PolyDiffContext logic).
        self._dclaurent = [
            self._claurent.derivative(i) for i in range(self._claurent.dim)
        ]
        self._d2claurent = [
            [self._dclaurent[i].derivative(j) for j in range(self._claurent.dim)] 
            for i in range(self._claurent.dim)
        ]
        # Pre-compute minor degrees for both directions.
        coeffs_ct = pt.CScalarVec([])
        degs_ct = pt.CIndexVec([])
        self._claurent.num.batch_get_data(coeffs_ct, degs_ct)
        for d in (1, 2):
            M_plus_N = max(degs_ct[d::3])
            M_val = self._claurent.denom_orders[d]
            N_val = M_plus_N - M_val
            if d == 2:
                self._M, self._N = M_val, N_val
            else:
                self._M1, self._N1 = M_val, N_val

    # -- Properties --

    @property
    def dim(self) -> int:
        """Number of variables (always 3)."""
        return self._claurent.dim

    @property
    def M(self) -> int:
        """Denominator order in beta2 (direction=2)."""
        return self._M

    @property
    def N(self) -> int:
        """Numerator max degree minus M in beta2 (direction=2)."""
        return self._N

    # -- Evaluation --

    def eval_val(self, var: tuple) -> complex:
        """Evaluate f(E, beta1, beta2) at the given variable tuple."""
        return self._claurent.eval(pt.CScalarVec(var))

    def eval_partials(self, var: tuple) -> list:
        """Evaluate all first partial derivatives at the given variable tuple."""
        var_ctype = pt.CScalarVec(var)
        return [dcl.eval(var_ctype) for dcl in self._dclaurent]
    
    def eval_partials_2(self, var: tuple, var_id1: int, var_id2: int):
        '''Evaluate second order partials'''
        var_ctype = pt.CScalarVec(var)
        return self._d2claurent[var_id1][var_id2].eval(var_ctype)

    def eval_dmu2(self, var: tuple) -> tuple:
        """Return (dmu2/dmu1, dmu2/dtheta1) from partial derivative ratios."""
        partials = self.eval_partials(var)
        complex_diff = var[1] * partials[1] / (var[2] * partials[2])
        return (-complex_diff.real, complex_diff.imag)

    # -- Degree info --

    def get_minor_degrees(self, direction: int = 2) -> tuple:
        """Extract (M, N) denominator / numerator minor degrees.

        Parameters:
            direction: Variable index (1 = beta1, 2 = beta2).  Default is 2.

        Returns:
            (M, N) tuple.
        """
        if direction == 2:
            return self._M, self._N
        return self._M1, self._N1

    # -- Root solving --

    def solve_roots_1d(self, param_indices, param_vals, var_indices, M=None, N=None) -> np.ndarray:
        """Partial-evaluate fixing param variables, solve 1D polynomial.

        Replaces the ``calculate_point_roots`` / ``partial_eval`` /
        ``batch_get_data`` pattern with a single method.

        Parameters:
            param_indices: e.g. ``(0, 1)`` to fix E and beta1.
            param_vals: values for the fixed variables, e.g. ``(E_ref, beta1)``.
            var_indices: e.g. ``(2,)`` to solve for beta2.
            M, N: minor degrees for root padding.  Defaults to ``self.M, self.N``
                (direction=2).

        Returns:
            Array of complex roots (unsorted).
        """
        if M is None:
            M = self._M
        if N is None:
            N = self._N
        poly_1d = self._claurent.partial_eval(
            pt.CScalarVec(param_vals),
            pt.CIndexVec(param_indices),
            pt.CIndexVec(var_indices),
        )
        coeffs_ct = pt.CScalarVec([])
        degs_ct = pt.CIndexVec([])
        poly_1d.num.batch_get_data(coeffs_ct, degs_ct)
        deg_M = poly_1d.denom_orders[0]

        # Convert poly_tools containers → numpy for np.roots.
        coeffs_list = list(coeffs_ct)
        degs_list = list(degs_ct)
        max_deg = max(degs_list)
        np_coeffs = np.zeros(max_deg + 1, dtype=complex)
        for c, d in zip(coeffs_list, degs_list):
            np_coeffs[max_deg - d] = c
        curr_roots = np.roots(np_coeffs)

        # Pad with 0 / inf for deficient root count.
        if len(curr_roots) < M + N:
            if deg_M < M:
                curr_roots = np.append(curr_roots, 0)
            if len(curr_roots) - deg_M < N:
                curr_roots = np.append(curr_roots, np.inf)
        return curr_roots


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
    """1D connected subset: a continuous GBZ curve segment.

    Stores the curve data eagerly (no lazy fill_beta2).  One instance
    represents a single β₂ curve; n-fold degeneracies produce n LineSubsets.

    Attributes:
        E: Reference energy.
        mu1: Fixed |beta1| radius (= ln|beta1|) across the segment.
        theta1_arr: (N,) θ₁ sampling points (monotonic).
        beta2_arr: (N,) β₂ values along this single curve.
    """
    E: complex
    mu1: float
    theta1_arr: np.ndarray   # (N,)
    beta2_arr: np.ndarray    # (N,)

    def __post_init__(self):
        if self.theta1_arr.ndim != 1:
            raise ValueError("theta1_arr must be 1-D")
        if self.beta2_arr.ndim != 1:
            raise ValueError("beta2_arr must be 1-D")
        if len(self.theta1_arr) != len(self.beta2_arr):
            raise ValueError("theta1_arr and beta2_arr must have same length")

    @property
    def theta1_start(self) -> float:
        return float(self.theta1_arr[0])

    @property
    def theta1_end(self) -> float:
        return float(self.theta1_arr[-1])

    @property
    def theta1_width(self) -> float:
        w = self.theta1_end - self.theta1_start
        return float(w if w > 0 else w + 2 * math.pi)

    @property
    def left_endpoint(self) -> tuple[complex, complex]:
        b1 = exp(self.mu1 + 1j * self.theta1_start)
        return (b1, self.beta2_arr[0])

    @property
    def right_endpoint(self) -> tuple[complex, complex]:
        b1 = exp(self.mu1 + 1j * self.theta1_end)
        return (b1, self.beta2_arr[-1])


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
    return roots[np.argsort(np.abs(roots))]


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


def cost_from_sphere_r3(p1: np.ndarray, p2: np.ndarray) -> np.ndarray:
    """Pairwise Euclidean distances between two sets of R^3 sphere points.

    Projection-free core of :func:`chordal_cost_matrix` — use it directly
    when the :func:`to_sphere_r3` projections are shared across several
    cost matrices (see ``_PmgbzScan.analyze_boundary_matching``).
    """
    diff = p1[:, None, :] - p2[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=2), dtype=float)


def chordal_cost_matrix(arr1: np.ndarray, arr2: np.ndarray) -> np.ndarray:
    """Chordal-distance cost matrix between two root arrays.

    The chordal distance is the Euclidean distance on the Riemann sphere,
    computed via :func:`to_sphere_r3`.

    Parameters:
        arr1, arr2: 1-D arrays of complex roots.

    Returns:
        2-D array ``cost[i, j]`` = chordal distance between arr1[i] and arr2[j].
    """
    return cost_from_sphere_r3(to_sphere_r3(arr1), to_sphere_r3(arr2))


def hungarian_match_indices(
    roots_from: np.ndarray,
    roots_to: np.ndarray,
) -> np.ndarray:
    """Hungarian (linear-sum-assignment) matching between two root sets.

    Uses chordal distance as the cost metric.  Raises ``ValueError`` if
    either array contains NaN entries.

    Returns:
        1-D int array ``perm`` such that ``perm[from_idx] = to_idx``, i.e.
        ``roots_from[i]`` is matched to ``roots_to[perm[i]]``.  Reorder the
        destination onto the source ordering with ``roots_to[perm]``.
    """
    if np.any(np.isnan(roots_from.real)) or np.any(np.isnan(roots_from.imag)):
        raise ValueError(f"NaN root encountered in roots_from: {roots_from}")
    if np.any(np.isnan(roots_to.real)) or np.any(np.isnan(roots_to.imag)):
        raise ValueError(f"NaN root encountered in roots_to: {roots_to}")

    cost = chordal_cost_matrix(roots_from, roots_to)
    row_ind, col_ind = linear_sum_assignment(cost)
    return col_ind[np.argsort(row_ind)]


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
    poly: CharPoly,
    direction: int = 2,
) -> tuple[int, int]:
    """Extract (M, N) denominator / numerator minor degrees.

    Thin wrapper around :meth:`CharPoly.get_minor_degrees`.  Kept for
    backward compatibility; prefer ``poly.M, poly.N`` in new code.

    Parameters:
        poly: CharPoly instance.
        direction: Variable index (1 = beta1, 2 = beta2).  Default is 2.

    Returns:
        (M, N) tuple.
    """
    return poly.get_minor_degrees(direction)


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
