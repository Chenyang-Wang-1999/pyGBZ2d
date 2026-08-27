'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-21
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

''' Solvers for multiple roots '''

import numpy as np
from typing import Optional, NamedTuple
from cmath import exp
# Chordal-distance threshold for detect_cluster — the single cluster
# predicate of the MR machinery (unified 2026-08 to the ZeroManager.run
# side value; the old direct-call default 1e-6 is retired).
CLUSTER_TOL: float = 1e-4
#: Closest-pair distance below which the MR interval trigger arms.
MIN_DIST_THRESHOLD: float = 0.1

from pygbz2d.core import live_defaults
from pygbz2d.core import (
    TWO_PI,
    CharPoly,
    hungarian_match_indices,
    to_sphere_r3,
    cost_from_sphere_r3,
)
from .arclength import (
    compute_tangent,
    predict_roots_hermite,
)
from scipy import optimize
from scipy.sparse.csgraph import connected_components

class MultipleRootInfo(NamedTuple):
    """Information about a detected multiple root."""
    theta1: float
    cluster_indices: list[tuple[int, ...]]  # indices into the modulus-sorted root array
    roots: np.ndarray            # modulus-sorted β₂ roots at this θ₁
    # Per-cluster spread of the *raw* roots before snapping to the mean
    # (``sqrt(mean(|β − mean|²))``), same order as ``cluster_indices``.
    # Empty when no cluster was detected.  A small value means the numerical
    # roots were already nearly coincident; a large one flags a loose cluster.
    # tuple (not list): a NamedTuple literal default must be immutable —
    # a shared mutable list default would alias across records.  Consumers
    # only read it.
    cluster_stds: tuple[float, ...] = ()


def snap_clusters_to_mean(
    roots: np.ndarray,
    cluster_indices: list[tuple[int, ...]],
) -> tuple[np.ndarray, list[float]]:
    """Snap every cluster's roots to their mean and report each cluster's spread.

    A detected multiple root is *exactly* degenerate in theory, but the
    numerical roots that land in a cluster are only approximately equal
    (finite solver tolerance, θ₁ bracketing error, etc.).  This enforces
    exact degeneracy by replacing every root in a cluster with the cluster's
    complex mean.  The standard deviation of the *original* roots in each
    cluster — ``sqrt(mean(|β − mean|²))``, a real RMS distance from the
    mean — is returned so the caller can record how loose the cluster was.

    Roots not in any cluster are left untouched.

    Parameters
    ----------
    roots : np.ndarray, shape (n_roots,)
        β₂ roots in whatever ordering the caller works in (modulus-sorted
        at the boundary MR, track-ordered at an in-loop MR).
    cluster_indices : list[tuple[int, ...]]
        Output of :func:`detect_cluster` — indices into *roots*.

    Returns
    -------
    snapped_roots : np.ndarray
        Copy of *roots* with each cluster's entries replaced by their mean.
    cluster_stds : list[float]
        Per-cluster standard deviation, same order as *cluster_indices*.
        Empty when *cluster_indices* is empty.
    """
    snapped = roots.copy()
    stds: list[float] = []
    for indices in cluster_indices:
        idx = np.array(indices, dtype=int)
        cluster_roots = roots[idx]
        mean = cluster_roots.mean()
        stds.append(float(np.std(cluster_roots)))
        snapped[idx] = mean
    return snapped, stds


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


@live_defaults(min_dtheta="continuation.zero_manager:MIN_DTHETA")
def multiple_root_point_trigger(
    dtheta: float,
    *,
    min_dtheta: Optional[float] = None,
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
    θ₁-derivative sign of min |β_i − β_j|² via :func:`_closest_pair_deriv`
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

    def __init__(self, min_dist_threshold: Optional[float] = None) -> None:
        if min_dist_threshold is None:
            min_dist_threshold = MIN_DIST_THRESHOLD
        self._min_dist_threshold = min_dist_threshold
        self._prev_deriv: int = 0          # −1=approaching, 0=unset, +1=separating
        self._prev_pair: tuple[int, int] | None = None
        # State of the *previous* step (the interval's left endpoint when a
        # trigger fires).  roots/V are kept alongside theta so the caller
        # can build a two-endpoint Hermite anchor without re-solving or
        # re-deriving the tangent at the interval start.
        self._prev_theta: float | None = None
        self._prev_roots: np.ndarray | None = None
        self._prev_V: np.ndarray | None = None

    def __call__(
        self, roots: np.ndarray, V: np.ndarray, theta1: float,
    ) -> tuple[bool, Optional[tuple[float, float]]]:
        """Check for sign flip.

        Returns
        -------
        (triggered, interval)
            *triggered* is True when the SAME pair that was closest at the
            previous step has its distance derivative flip from negative
            to positive, AND the current minimum distance is below the
            threshold.  *interval* is ``(start, end)`` — the θ₁ range
            containing the local minimum — or None if not triggered.

            On a trigger, ``self._prev_roots`` / ``self._prev_V`` /
            ``self._prev_theta`` hold the interval's *left* endpoint state
            (they are not overwritten before the trigger returns).
        """
        min_dist, deriv, pair = _closest_pair_deriv(roots, V)

        # Only track when roots are close enough for a meaningful signal.
        if min_dist >= self._min_dist_threshold:
            self._reset()
            return False, None

        # Same-pair guard: avoid spurious sign flips when the closest pair
        # changes identity between steps.
        triggered = (
            self._prev_deriv < 0 and deriv > 0
            and pair == self._prev_pair
        )
        if triggered:
            # Early-return WITHOUT overwriting prev_*: the interval's left
            # endpoint is the previous step's state, which the caller reads.
            return True, (self._prev_theta, theta1)

        self._prev_deriv = deriv
        self._prev_pair = pair
        self._prev_theta = theta1
        self._prev_roots = roots
        self._prev_V = V

        return False, None

    def _reset(self) -> None:
        self._prev_deriv = 0
        self._prev_pair = None
        self._prev_theta = None
        self._prev_roots = None
        self._prev_V = None


@live_defaults(cluster_tol="continuation.multiple_roots:CLUSTER_TOL")
def detect_cluster(
    roots: np.ndarray,
    *,
    cluster_tol: Optional[float] = None,
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


def _pair_distance_deriv(
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


def _closest_pair_deriv(
    roots: np.ndarray,
    V: np.ndarray,
) -> tuple[float, float, tuple[int, int]]:
    """Compute min Euclidean distance among roots, its θ₁-derivative,
    and the pair achieving the minimum.

    Roots whose tangent is nan (the 0/∞ padding roots held fixed by the
    integrator) or inf (divergent, at a branch point) are excluded from
    the pair search: a nan/inf derivative would poison the interval
    trigger's ``_prev_deriv`` and blind it permanently (``nan < 0`` is
    False forever after).  The point trigger owns the branch-point regime.

    Returns
    -------
    min_dist : float
        Minimum |β_i − β_j| over finite-tangent root pairs (``inf`` when
        fewer than two such roots exist).
    deriv : float
        Raw d|β_i − β_j|²/dθ₁ of that pair; only its sign is consumed.
    pair : tuple[int, int]
        Track indices (i, j) of the closest pair, or (-1, -1) when none.
    """
    n = len(roots)
    finite = [
        j for j in range(n)
        if np.isfinite(V[j].real) and np.isfinite(V[j].imag)
    ]

    min_dist = np.inf
    min_i, min_j = -1, -1

    for a in range(len(finite)):
        for b in range(a + 1, len(finite)):
            i, j = finite[a], finite[b]
            d = np.abs(roots[i] - roots[j])
            if d < min_dist:
                min_dist = d
                min_i, min_j = i, j

    if min_i < 0:
        return float(min_dist), 0.0, (-1, -1)

    deriv = _pair_distance_deriv(roots, V, (min_i, min_j))

    return float(min_dist), deriv, (min_i, min_j)


def solve_multiple_roots_in_interval(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    theta1_left: float,
    theta1_right: float,
    roots_ref: np.ndarray,
    min_pair: tuple[int, int] = None,
):
    """Locate a multiple root inside a θ₁ interval via Brent's method.

    The θ₁-derivative of the minimum pairwise distance among β₂ roots
    crosses zero at a multiple root (roots stop approaching and start
    separating).  This function locates that zero-crossing with Brent's
    method.  Cluster verification at the located θ₁ is left to the
    caller (see :func:`detect_cluster`).

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

    Returns
    -------
    theta1_mr : float
        The θ₁ value where the derivative crosses zero, i.e. the
        multiple-root location (not wrapped to [0, 2π)).
    """
    def _compute_deriv(theta1: float):
        beta1 = exp(mu1 + 1j * theta1)
        roots = poly.solve_roots_1d((0, 1), (E_ref, beta1), (2,))
        # Match predicted → solved, not ref → solved (a bare match that
        # swaps near-degenerate tracks).  roots_ref sits at theta1_right;
        # extrapolate its tangent to the trial θ₁.  Δt is small (Brent's
        # bracket is two consecutive integration steps), so first-order
        # tangent extrapolation is adequate.
        predicted = predict_roots_hermite(theta1, theta1_right, roots_ref, V_ref)
        inds = hungarian_match_indices(predicted, roots)
        roots = roots[inds]
        V_list, _ = compute_tangent(poly, E_ref, beta1, roots)
        _, deriv, new_pair = _closest_pair_deriv(roots, V_list)
        if min_pair is not None:
            if new_pair != min_pair:
                raise ValueError(f"Pair {new_pair} is not the minimum pair {min_pair}")
        return deriv

    if theta1_right < theta1_left:
        theta1_right += TWO_PI

    # Tangent at the reference (right) endpoint — reused across Brent trials.
    V_ref, _ = compute_tangent(
        poly, E_ref, exp(mu1 + 1j * theta1_right), roots_ref,
    )

    theta1_mr = optimize.brentq(_compute_deriv, theta1_left, theta1_right)

    return theta1_mr


def solve_multiple_roots_iterative(
    poly: CharPoly,
    E_ref: complex,
    beta1_approx: complex,
    beta2_approx: complex,
):
    """Refine a multiple root by solving f = 0 and ∂f/∂β₂ = 0 simultaneously.

    A multiple root in β₂ is a point where f(E, β₁, β₂) = 0 *and* the
    β₂-derivative vanishes (∂f/∂β₂ = 0) — i.e. f has a repeated β₂-root.
    With E fixed, (β₁, β₂) are the 2 complex unknowns (4 real), matched by
    the 4 real equations Re/Im of {f, ∂f/∂β₂}.  Solved as a real nonlinear
    system via :func:`scipy.optimize.root` with an analytic Jacobian.
    """
    def _mr_fun(x: np.ndarray):
        beta1 = complex(x[0], x[1])
        beta2 = complex(x[2], x[3])
        curr_pt = (E_ref, beta1, beta2)
        f_val = poly.eval_val(curr_pt)
        df = poly.eval_partials(curr_pt)            # [∂f/∂E, ∂f/∂β₁, ∂f/∂β₂]
        df_12 = poly.eval_partials_2(curr_pt, 1, 2)  # ∂²f/∂β₁∂β₂
        df_22 = poly.eval_partials_2(curr_pt, 2, 2)  # ∂²f/∂β₂²

        # Residuals: [Re f, Im f, Re ∂f/∂β₂, Im ∂f/∂β₂].
        eq_val = np.array(
            [f_val.real, f_val.imag, df[2].real, df[2].imag]
        )

        # Complex Jacobian of (f, ∂f/∂β₂) w.r.t. (β₁, β₂):
        #   row 0 = [∂f/∂β₁,    ∂f/∂β₂]
        #   row 1 = [∂²f/∂β₁∂β₂, ∂²f/∂β₂²]
        eq_jac_cc = np.array([
            [df[1], df[2]],
            [df_12, df_22]
        ])

        # Convert to derivatives w.r.t. the real coordinates
        # (Re β₁, Im β₁, Re β₂, Im β₂): dβ = dReβ + i·dImβ, so the map
        # (Re β₁, Im β₁, Re β₂, Im β₂) → (β₁, β₂) is the matrix below.
        eq_jac_cr = eq_jac_cc @ np.array([
            [1, 1j, 0, 0],
            [0, 0, 1, 1j]
        ])
        # Split each complex row into [Re; Im] → 4×4 real Jacobian.
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

