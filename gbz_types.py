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
from typing import Optional, Union

import numpy as np
import poly_tools as pt
from scipy.optimize import linear_sum_assignment

# The single 2π constant for the whole project.  Every module imports it
# from here instead of spelling `2 * pi` locally, so all seam comparisons
# (θ % 2π, circ_dist, closing-row θ=2π, boundary-MR tolerance) operate on
# one bit-identical float.  Equal to math.tau; asserted in tests.
TWO_PI: float = 2.0 * math.pi

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

        Replaces the old ``partial_eval`` / ``batch_get_data`` pattern with a
        single method.

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

        # Convert poly_tools containers → numpy for np.roots.  A completely
        # vanished partial polynomial (all β₂ terms carried a parameter that
        # became exactly 0) has no terms; its roots are entirely the padding
        # roots below, so max() must not be attempted.  poly_tools keeps the
        # global denominator order on that empty object, but for padding the
        # *actual* denominator degree is zero.
        coeffs_list = list(coeffs_ct)
        degs_list = list(degs_ct)
        if degs_list:
            deg_M = int(poly_1d.denom_orders[0])
            max_deg = max(degs_list)
            np_coeffs = np.zeros(max_deg + 1, dtype=complex)
            for c, d in zip(coeffs_list, degs_list):
                np_coeffs[max_deg - d] = c
            curr_roots = np.roots(np_coeffs)
        else:
            deg_M = 0
            max_deg = 0
            curr_roots = np.array([], dtype=complex)

        # A parameter value can cancel more than one leading/trailing term.
        # Count, rather than merely detect, the deficiency:
        #   * cancelled low-side terms → roots at β₂ = 0;
        #   * cancelled high-side terms → roots at β₂ = ∞.
        # After clearing the β₂ denominator, the finite polynomial degree is
        # max_deg, while deg_M + N is the degree of the uncancellated
        # Laurent numerator.  Their difference is the ∞-root count, and
        # M - deg_M is the 0-root count.
        target_count = M + N
        n_zero = max(0, int(M) - deg_M)
        n_inf = max(0, deg_M + int(N) - int(max_deg))
        if len(curr_roots) + n_zero + n_inf != target_count:
            raise RuntimeError(
                f"root padding failed: polynomial degree {max_deg}, "
                f"denominator degree {deg_M}, M={M}, N={N} give "
                f"{len(curr_roots)} finite + {n_zero} zero + {n_inf} "
                f"infinite roots (expected {target_count})"
            )

        if n_zero:
            curr_roots = np.append(curr_roots, np.zeros(n_zero, dtype=complex))
        if n_inf:
            curr_roots = np.append(
                curr_roots, np.full(n_inf, np.inf, dtype=complex))
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
        return float(math.atan2(self.beta1.imag, self.beta1.real) % (TWO_PI))

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
        return float(w if w > 0 else w + TWO_PI)

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
               outside the GBZ — unless ``is_continuum`` is set (see below).
        is_continuum: True when the SGBZ subset at this energy is a continuum
               (1D).  SGBZ materializes the LineSubsets into ``subsets``
               (``index == (0, n_1d)``); the flag additionally marks the
               result as *in spectrum* (``is_gbz`` is True) even if
               ``subsets`` happens to be empty.  Amoeba (which has its own
               LineSubset extractor) never sets this flag.
    """
    E_ref: complex
    success: bool = True
    error: Optional[str] = None
    subsets: list[ConnectedSubset] = None  # type: ignore[assignment]
    index: tuple[int, int] = (0, 0)
    is_continuum: bool = False

    def __post_init__(self):
        if self.subsets is None:
            object.__setattr__(self, 'subsets', [])

    @property
    def is_empty(self) -> bool:
        return len(self.subsets) == 0

    @property
    def is_gbz(self) -> bool:
        """True if this energy point lies on the GBZ."""
        return self.success and (self.index != (0, 0) or self.is_continuum)


ConnectedSubset = Union[PointSubset, LineSubset]


# ---------------------------------------------------------------------------
# Cross-module LineSubset joining helpers
# ---------------------------------------------------------------------------
#
# Both the amoeba extractor (brute_force_amoeba.zm_extract) and the SGBZ
# continuum extractor (brute_force_SGBZ.continuum_lines) join per-segment
# continuum LineSubsets across MR boundaries into closed curves.  The join
# unit and the MR-cluster endpoint test are module-agnostic (they only need
# the ZeroManager protocol: segments with left_mr/right_mr, multiple_roots),
# so they live here — no extractor should import the other's private pieces.

class JoinableLinePiece(LineSubset):
    """A per-segment continuum LineSubset being joined across MR boundaries.

    Carries ``ml``/``mr`` — the leftmost/rightmost original segment indices
    spanned so far — so merges can be chained and the join at the cyclic seam
    (segment 0 ↔ last segment) detected.  Behaviourally a ``LineSubset`` once
    joining is done.
    """

    def __init__(self, E, mu1, theta1_arr, beta2_arr, ml: int, mr: int):
        super().__init__(E=E, mu1=mu1,
                         theta1_arr=theta1_arr, beta2_arr=beta2_arr)
        self.ml = ml
        self.mr = mr


def is_mr_cluster_endpoint(zm, seg, side: str, col: int) -> bool:
    """Whether endpoint track *col* belongs to the boundary MR cluster.

    A segment boundary is an MR, but only the tracks listed in
    ``multiple_roots[mr].cluster_indices`` are genuinely multiple there; every
    other track is regular and passes straight through.  ``mr`` is the
    segment's ``left_mr`` / ``right_mr``; ``mr < 0`` is the only "no MR" case
    — it marks the θ₁=0/2π circle seam (segment 0's left / last segment's
    right, set to -1 by ``ZeroManager.run``).

    Membership is decided by COLUMN IDENTITY, not by a nearest-value match.
    Interior MR records share the track frame of both adjacent boundary rows.
    The one exception is the θ=0 boundary MR reused as the final segment's
    right boundary at θ=2π: its cluster indices use the θ=0 modulus-sorted
    frame, while ``col`` indexes the final segment's track frame.  The class
    convention ``roots_right[boundary_perm] == roots_left`` translates left
    column ``k`` to right column ``boundary_perm[k]`` there.
    """
    mr = seg.left_mr if side == 'left' else seg.right_mr
    if mr < 0:
        return False
    cluster = zm.multiple_roots[mr].cluster_indices
    if not cluster:
        return False

    out_col = int(col)
    if (side == 'right' and mr == 0 and zm.has_boundary_mr):
        if not hasattr(zm, "boundary_perm"):
            raise RuntimeError(
                "cannot test the θ=2π boundary MR cluster before "
                "boundary_perm is set"
            )
        inv_perm = np.empty(zm.K, dtype=int)
        inv_perm[zm.boundary_perm] = np.arange(zm.K)
        out_col = int(inv_perm[out_col])

    return any(out_col in c for c in cluster)


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


def circ_dist(a: float, b: float) -> float:
    """Shortest arc distance on any 2π-periodic circle (θ₁ or θ₂).

    Needed because a boundary MR sits at θ₁ = 0 while its neighbourhood
    extends to θ₁ ≈ 2π; likewise for θ₂ boundaries near the seam.  Robust
    to inputs outside [0, 2π) via the leading ``% (2π)``.
    """
    d = abs(a - b) % (TWO_PI)
    return min(d, TWO_PI - d)


# ---------------------------------------------------------------------------
# Zero-plateau detection (shared by SGBZ and Amoeba)
# ---------------------------------------------------------------------------
#
# A zero-plateau boundary sits at the edge of the spectrum: the winding
# changes sign over a vanishingly narrow angular region, so the GBZ points
# cluster into nearly degenerate pairs rather than partitioning the circle.
# Both modules detect this the same way — a torus-distance clustering
# pre-check, then a geometric probe ladder in mu1 ± step — differing only in
# how a single probe is *evaluated* (amoeba: mu2 bisection + Ronkin winding;
# SGBZ: ZeroManager build + crossing winding).  The loop and the
# classification are module-agnostic, so they live here.

def check_points_clustered_on_torus(
    points: list[tuple[float, float]],
    tol_normalized: float,
) -> bool:
    """Whether every (θ₁, θ₂) point has a neighbour within *tol_normalized*.

    Distance is the Euclidean metric on the (θ₁, θ₂)-torus ``[0, 2π)²``,
    normalized by ``2π`` so the full torus diagonal is ``√2``.  At a genuine
    GBZ point the points are well-separated (they partition the circle into
    meaningful segments); at a zero-plateau boundary they cluster into nearly
    degenerate pairs, each within *tol_normalized* of a neighbour.

    Returns ``False`` when any point is isolated (no neighbour within the
    threshold), which rules out a plateau and lets the caller skip the
    expensive probe.  ``tol_normalized`` is in units of the torus period
    (i.e. ``1.0`` = one full ``2π`` circle).

    Parameters:
        points: list of ``(theta1, theta2)`` in radians.
        tol_normalized: neighbour threshold, in units of ``2π``.
    """
    if len(points) < 2:
        return False
    tol_rad = tol_normalized * (TWO_PI)
    for i, (t1_i, t2_i) in enumerate(points):
        has_neighbor = False
        for j, (t1_j, t2_j) in enumerate(points):
            if i == j:
                continue
            d1 = circ_dist(t1_i, t1_j)
            d2 = circ_dist(t2_i, t2_j)
            if math.sqrt(d1 * d1 + d2 * d2) < tol_rad:
                has_neighbor = True
                break
        if not has_neighbor:
            return False
    return True


def probe_zero_plateau(
    mu1: float,
    mu1_bracket: Optional[tuple[float, float]],
    *,
    zero_tol: float,
    probe_radius: Optional[float],
    bracket_width: Optional[float],
    evaluator,
) -> dict:
    """Geometric probe ladder checking ``mu1 ± step`` for a zero plateau.

    Walks outward in ``±step`` over the ladder from
    :func:`generate_probe_steps`.  A plateau is declared as soon as one probe
    is a plateau point.  *evaluator* is a callable ``mu1_probe -> dict`` whose
    return carries the module-specific probe; the shared loop reads only:

      * ``success`` (bool) — the probe completed.
      * ``is_continuum`` (bool) — the probe landed on a continuum (undefined
        winding, neither plateau nor ordinary GBZ).
      * ``is_plateau`` (bool) — the module-specific plateau criterion already
        collapsed to one boolean (empty GBZ + zero winding, however the
        module defines "winding" / "GBZ count").

    Both sides of a step are probed before the plateau check, so the returned
    ``points`` are symmetric and the left/right non-plateau bookkeeping is
    correct.

    Returns a dict with ``found`` (bool), ``status``
    (``"found"``/``"not_found"``/``"inconclusive"``), ``steps``, and the
    per-probe ``points`` list.
    """
    steps = generate_probe_steps(
        float(bracket_width or 0.0),
        float(probe_radius or 0.0),
        zero_tol,
    )
    probe_points: list[dict] = []
    found_plateau = False
    saw_left_nonplateau = False
    saw_right_nonplateau = False

    for step in steps:
        for side in (-1, 1):
            mu1_probe = mu1 + side * step
            point = {"mu1": mu1_probe, "side": side, "step": step}
            point.update(evaluator(mu1_probe))
            probe_points.append(point)
            if point.get("is_plateau"):
                found_plateau = True
            elif point.get("success") and not point.get("is_continuum"):
                if side < 0:
                    saw_left_nonplateau = True
                else:
                    saw_right_nonplateau = True
        if found_plateau:
            return {
                "status": "found",
                "found": True,
                "zero_tol": zero_tol,
                "probe_radius": probe_radius,
                "bracket_width": float(bracket_width or 0.0),
                "steps": steps,
                "points": probe_points,
            }

    if found_plateau:
        status = "found"
    elif saw_left_nonplateau and saw_right_nonplateau:
        status = "not_found"
    else:
        status = "inconclusive"
    return {
        "status": status,
        "found": found_plateau,
        "zero_tol": zero_tol,
        "probe_radius": probe_radius,
        "bracket_width": float(bracket_width or 0.0),
        "steps": steps,
        "points": probe_points,
    }
