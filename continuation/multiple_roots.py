'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-21
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

''' Solvers for multiple roots '''

import numpy as np
from typing import Optional, NamedTuple
from cmath import exp, pi
from gbz_types import (
    CharPoly,
    hungarian_match_indices
)
from gbz_types import (
    CharPoly,
    to_sphere_r3,
    cost_from_sphere_r3,
)
from .arclength import (
    compute_tangent
)
from scipy import optimize

class MultipleRootInfo(NamedTuple):
    """Information about a detected multiple root."""
    theta1: float
    cluster_indices: list[tuple[int, ...]]  # indices into the modulus-sorted root array
    roots: np.ndarray            # modulus-sorted β₂ roots at this θ₁


# ---------------------------------------------------------------------------
# Multiple-root detection & refinement
# ---------------------------------------------------------------------------
#
# Two complementary triggers detect multiple roots during integration:
#
#   point_trigger  — step size collapses → we are AT the MR.  The
#                    implicit-function derivative diverges (∂f/∂β₂ → 0),
#                    so any finite arclength step in θ₁ space becomes tiny.
#
#   interval_trigger — closest-pair distance derivative flips sign
#                    (approaching → separating) → an MR lies BETWEEN
#                    two successive steps.  This catches "accidental"
#                    multiple roots where the tangent stays well-conditioned
#                    and the step size never collapses.
#
# Both are called inside ``integrate_segment`` after each accepted step.
# The distinction is surfaced through ``SegmentResult.stop_reason`` so
# that ZeroManager can handle each case appropriately (direct refine vs.
# bisection + refine).


def multiple_root_point_trigger(
    dtheta: float,
    *,
    min_dtheta: float = 1e-10,
) -> bool:
    """Trigger when the arclength step size has collapsed below *min_dtheta*.

    A collapsed step indicates that the implicit-function derivative
    diverged (∂f/∂β₂ → 0), meaning we are effectively **at** a multiple
    root.  ZeroManager should refine directly at *mr_approx_theta*.
    """
    return dtheta < min_dtheta


class MultipleRootIntervalTrigger:
    """Detect multiple roots by tracking the closest-pair distance derivative.

    Called after each accepted integration step.  Internally computes the
    θ₁-derivative sign of min |β_i − β_j|² via :func:`_min_pairwise_deriv`
    and tracks it across steps.  Triggers when the sign flips from negative
    (approaching) to positive (separating), indicating a local minimum —
    a multiple root — in the interval between the previous and current θ₁.

    Tracking is only active when the minimum Euclidean distance among roots
    is below *min_dist_threshold*; above it the internal state is reset.
    This avoids triggering on trajectory wiggles among well-separated roots.

    Usage::

        trigger = MultipleRootIntervalTrigger(min_dist_threshold=0.1)
        for step in integration:
            ...
            ok, interval = trigger(roots, V, theta1)
            if ok:
                bisect(*interval)  # MR candidate in [interval[0], interval[1]]
    """

    def __init__(self, min_dist_threshold: float = 0.1) -> None:
        self._min_dist_threshold = min_dist_threshold
        self._prev_deriv: int = 0          # −1=approaching, 0=unset, +1=separating
        self._prev_pair: tuple[int, int] | None = None
        self._prev_theta: float | None = None

    def __call__(
        self, roots: np.ndarray, V: np.ndarray, theta1: float,
    ) -> tuple[bool, tuple[float, float]]:
        """Check for sign flip.

        Returns
        -------
        (triggered, interval)
            *triggered* is True when the SAME pair that was closest at the
            previous step has its distance derivative flip from negative
            to positive, AND the current minimum distance is below the
            threshold.  *interval* is ``(start, end)`` — the θ₁ range
            containing the local minimum — or None if not triggered.
        """
        min_dist, deriv, pair = _min_pairwise_deriv(roots, V)

        # Only track when roots are close enough for a meaningful signal.
        if min_dist >= self._min_dist_threshold:
            self._prev_deriv = 0
            self._prev_pair = None
            self._prev_theta = None
            return False, None

        # Same-pair guard: avoid spurious sign flips when the closest pair
        # changes identity between steps.
        triggered = (
            self._prev_deriv < 0 and deriv > 0
            and pair == self._prev_pair
        )
        interval = (self._prev_theta, theta1) if triggered else None

        self._prev_deriv = deriv
        self._prev_pair = pair
        self._prev_theta = theta1

        return triggered, interval


def detect_cluster(
    roots: np.ndarray,
    *,
    cluster_tol: float = 1e-6,
) -> list[tuple[int, ...]]:
    """Find all root clusters as connected components of the proximity graph.

    Two roots are connected when their chordal distance is below
    *cluster_tol*.  Clusters are the connected components of this graph,
    computed via :func:`scipy.sparse.csgraph.connected_components`.

    Parameters
    ----------
    roots : np.ndarray, shape (n_roots,)
        β₂ roots (any ordering).
    cluster_tol : float
        Chordal-distance threshold below which two roots are considered
        "touching".

    Returns
    -------
    list[tuple[int, ...]]
        Each tuple holds the sorted indices of one cluster.  Empty list
        when all roots are well-separated.
    """
    from scipy.sparse.csgraph import connected_components

    r3 = to_sphere_r3(roots)
    pw = cost_from_sphere_r3(r3, r3)

    adj = pw < cluster_tol
    np.fill_diagonal(adj, False)

    n_components, labels = connected_components(adj, directed=False)

    clusters: list[tuple[int, ...]] = []
    for label in range(n_components):
        indices = np.where(labels == label)[0]
        if len(indices) >= 2:
            clusters.append(tuple(indices))

    return clusters


def _pairwise_deriv(
    roots: np.ndarray,
    V: np.ndarray,
    pair: tuple[int, int],
):
    ''' 
        Compute the pairwise derivative of the pair distance against theta1 

        d/dθ₁ |β_i − β_j|² = 2 Re[(β̇_i − β̇_j) · (β_i − β_j)†]
        where β̇_k = V_k · β_k (from the tangent vector).

    '''
    min_i, min_j = pair

    beta_dot_i = V[min_i] * roots[min_i]
    beta_dot_j = V[min_j] * roots[min_j]

    return 2.0 * np.real(
        (beta_dot_i - beta_dot_j) * np.conj(roots[min_i] - roots[min_j])
    )


def _min_pairwise_deriv(
    roots: np.ndarray,
    V: np.ndarray,
) -> tuple[float, float, tuple[int, int]]:
    """Compute min Euclidean distance among roots, its θ₁-derivative,
    and the pair achieving the minimum.

    Returns
    -------
    min_dist : float
        Minimum |β_i − β_j| over all root pairs.
    deriv_sign : int
        −1 (approaching), 0 (stationary), or +1 (separating).
    pair : tuple[int, int]
        Track indices (i, j) of the closest pair.
    """
    n = len(roots)
    min_dist = np.inf
    min_i, min_j = -1, -1

    for i in range(n):
        for j in range(i + 1, n):
            d = np.abs(roots[i] - roots[j])
            if d < min_dist:
                min_dist = d
                min_i, min_j = i, j

    deriv = _pairwise_deriv(roots, V, (min_i, min_j))

    return float(min_dist), deriv, (min_i, min_j)


def solve_multiple_roots_in_interval(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    theta1_left: float,
    theta1_right: float,
    roots_ref: np.ndarray,
    min_pair: tuple[int, int] = None,
    cluster_tol: float = 1e-6,
):
    """Locate a multiple root inside a θ₁ interval via Brent's method.

    The θ₁-derivative of the minimum pairwise distance among β₂ roots
    crosses zero at a multiple root (roots stop approaching and start
    separating).  This function brackets that zero-crossing with Brent's
    method and verifies the result with :func:`detect_cluster`.

    ** Note ** we assume that the multiple roots are well-separated in theta1 axis.
    If any exceptions are raised in this function, it probably means that several multiple roots are close to each other in theta1 axis.
    You can use try-except to avoid these exceptions if the bunching of multiple roots do not influence your calculations.

    Parameters
    ----------
    poly : CharPoly
        Characteristic polynomial f(E, β₁, β₂).
    E_ref : complex
        Fixed reference energy.
    mu1 : float
        Log-modulus of β₁ (β₁ = exp(μ₁ + iθ₁)).
    theta1_left, theta1_right : float
        Bracketing interval.  The derivative must have opposite signs at
        the two endpoints: ``f'(left) < 0`` (approaching) and
        ``f'(right) > 0`` (separating), or vice versa.
    min_pair : tuple[int, int] or None
        If given, the *only* root pair whose derivative is tracked.
        A :class:`ValueError` is raised if the closest pair changes
        identity inside the interval.  When None (default), any pair
        is accepted — useful when the identity of the merging pair is
        unknown ahead of time.
    cluster_tol : float
        Chordal-distance threshold passed to :func:`detect_cluster`.

    Returns
    -------
    theta1_mr : float
        The θ₁ value where the derivative crosses zero, i.e. the
        multiple-root location (not wrapped to [0, 2π)).
    clusters : list[tuple[int, ...]]
        Clusters detected at *theta1_mr*, as returned by
        :func:`detect_cluster`.  Empty if the bracket contained a
        derivative zero-crossing without roots actually touching
        (e.g. a near-miss).
    """
    def _compute_deriv(theta1: float):
        beta1 = exp(mu1 + 1j * theta1)
        roots = poly.solve_roots_1d((0, 1), (E_ref, beta1), (2,))
        inds = hungarian_match_indices(roots_ref, roots)
        roots = roots[inds]
        V_list, _ = compute_tangent(poly, E_ref, beta1, roots)
        _, deriv, new_pair = _min_pairwise_deriv(roots, V_list)
        if min_pair is not None:
            if new_pair != min_pair:
                raise ValueError(f"Pair {new_pair} is not the minimum pair {min_pair}")
        return deriv

    if theta1_right < theta1_left:
        theta1_right += 2 * pi

    theta1_mr = optimize.brentq(_compute_deriv, theta1_left, theta1_right)

    return theta1_mr


def solve_multiple_roots_iterative(
    poly: CharPoly,
    E_ref: complex,
    beta1_approx: complex,
    beta2_approx: complex,
):
    ''' solve multiple roots by iterative solver '''
    def _mr_fun(x: np.ndarray):
        beta1 = complex(x[0], x[1])
        beta2 = complex(x[2], x[3])
        curr_pt = (E_ref, beta1, beta2)
        f_val = poly.eval_val(curr_pt)
        df = poly.eval_partials(curr_pt)
        df_12 = poly.eval_partials_2(curr_pt, 1, 2)
        df_22 = poly.eval_partials_2(curr_pt, 2, 2)

        eq_val = np.array(
            [f_val.real, f_val.imag, df[2].real, df[2].imag]
        )

        eq_jac_cc = np.array([
            [df[1], df[2]],
            [df_12, df_22]
        ])

        eq_jac_cr = eq_jac_cc @ np.array([
            [1, 1j, 0, 0],
            [0, 0, 1, 1j]
        ])
        eq_jac = np.vstack([
            eq_jac_cr[0, :].real,
            eq_jac_cr[0, :].imag,
            eq_jac_cr[1, :].real,
            eq_jac_cr[1, :].imag,
        ])
        return eq_val, eq_jac
    
    res = optimize.root(
        _mr_fun, 
        np.array([ 
                beta1_approx.real, 
                beta1_approx.imag, 
                beta2_approx.real, 
                beta2_approx.imag]), 
        jac=True
    )

    if not res.success:
        raise ValueError(
            f"Root solver failed. Message={res.message}"
        )

    x = res.x
    return (
        complex(x[0], x[1]),
        complex(x[2], x[3]),
    )

