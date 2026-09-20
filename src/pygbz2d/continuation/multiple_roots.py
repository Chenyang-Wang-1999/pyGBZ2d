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
    _is_singular_root,
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
    geometric mean: branch-unwrapped ``ln β₂`` is averaged, then ``exp``
    maps the mean back to β₂.  The standard deviation of the *original*
    branch-unwrapped ``ln β₂`` values is returned so the caller can record
    how loose the cluster was.

    Roots not in any cluster are left untouched.

    **Guard for 0/∞ padding roots.**  A cluster containing a singular root
    (exact ``0`` / ``∞`` padding, or a root classified by
    :func:`_is_singular_root`) is left untouched and its ``cluster_stds``
    entry is ``nan``.  These are degree-deficiency roots at the boundary of
    the β₂-sphere, not finite branch points; snapping them is meaningless
    (``np.mean`` over ``∞`` yields ``∞+nanj``) and, worse, a mixed cluster
    such as ``[0, 1e-5]`` or ``[∞, 2e4]`` would drag a genuine finite root
    to a fake mean.  ZeroManager only needs finite MRs to prevent track
    swapping, so these clusters are deliberately not snapped.

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
        Copy of *roots* with each finite cluster's entries replaced by their
        geometric mean ``exp(mean(ln β₂))``; clusters containing a 0/∞
        padding root are left unchanged.
    cluster_stds : list[float]
        Per-cluster standard deviation of the branch-unwrapped ``ln β₂``
        values, same order as *cluster_indices*.  ``nan`` for clusters
        skipped by the 0/∞ guard.  Empty when *cluster_indices* is empty.
    """
    snapped = roots.copy()
    stds: list[float] = []
    for indices in cluster_indices:
        idx = np.array(indices, dtype=int)
        cluster_roots = roots[idx]

        # Guard: do not snap clusters involving 0/∞ padding roots.  See the
        # docstring for the mixed-cluster corruption this prevents.
        if any(_is_singular_root(r) for r in cluster_roots):
            stds.append(float('nan'))
            continue

        # Geometric mean via branch-unwrapped log.  Plain np.log uses the
        # principal branch, so a cluster straddling the negative real axis
        # would have log imag parts split across ±π; unwrapping relative to
        # the first member puts all log values on the same branch before
        # averaging.
        logs = np.log(cluster_roots)
        ref_imag = float(logs[0].imag)
        shifts = 1j * TWO_PI * np.round((ref_imag - logs.imag) / (TWO_PI))
        unwrapped = logs + shifts
        mean = np.exp(unwrapped.mean())
        stds.append(float(np.std(unwrapped)))
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


class MRTriggerRecord(NamedTuple):
    """One root pair whose pairwise distance derivative flipped sign.

    Emitted by :class:`MultipleRootIntervalTrigger` when a pair went from
    approaching (deriv < 0) at *theta_lo* to separating (deriv > 0) at
    *theta_hi*: a local minimum of that pair's distance — a multiple root
    of those two tracks — lies inside ``(*theta_lo*, *theta_hi*)``.
    """

    pair: tuple[int, int]        # track indices (i < j)
    theta_lo: float              # θ₁ of the previous accepted row
    theta_hi: float              # θ₁ of the current accepted row
    dist_hi: float               # |β_i − β_j| at theta_hi (diagnostics)
    deriv_lo: float              # d|β_i−β_j|²/dθ₁ at theta_lo
    deriv_hi: float              # d|β_i−β_j|²/dθ₁ at theta_hi


class MultipleRootIntervalTrigger:
    """Detect multiple roots by tracking EVERY pair's distance derivative.

    Called after each accepted integration step.  Computes the θ₁-derivative
    of ``|β_i − β_j|²`` for all root pairs at once (vectorized
    :func:`_pairwise_dist_deriv`) and tracks each pair's sign across steps.
    A pair triggers when ITS OWN derivative flips from negative
    (approaching) to positive (separating) — a local minimum of that pair's
    distance, i.e. a multiple root of those two tracks, lies between the
    previous and current θ₁.

    Per-pair tracking replaces the old single closest-pair state with its
    same-pair guard: when two pairs degenerate simultaneously (a symmetric
    double MR) their distances are tied, the *argmin* pair identity
    flickers row by row, and the old guard blocked the flip detection
    exactly at the sign change.  Each pair carrying its own state makes
    the argmin identity irrelevant.

    A pair's tracking is armed only while its own distance is below
    *min_dist_threshold*; above it that pair's state resets, avoiding
    triggers on trajectory wiggles among well-separated roots.

    Usage::

        trigger = MultipleRootIntervalTrigger(min_dist_threshold=0.1)
        for step in integration:
            ...
            records = trigger(roots, V, theta1)
            if records:
                # one or more MRs inside (records[0].theta_lo,
                # records[0].theta_hi); all records share that interval
    """

    def __init__(self, min_dist_threshold: Optional[float] = None) -> None:
        if min_dist_threshold is None:
            min_dist_threshold = MIN_DIST_THRESHOLD
        self._min_dist_threshold = min_dist_threshold
        # Distance derivative per pair, symmetric (n, n); 0 = unarmed.
        # Floats (only their signs are consumed), matching the raw
        # comparisons of the old single-pair state.
        self._prev_deriv: Optional[np.ndarray] = None
        # State of the *previous* step (the interval's left endpoint when a
        # trigger fires).  roots/V are kept alongside theta so the caller
        # can build a two-endpoint Hermite anchor without re-solving or
        # re-deriving the tangent at the interval start.
        self._prev_theta: Optional[float] = None
        self._prev_roots: Optional[np.ndarray] = None
        self._prev_V: Optional[np.ndarray] = None

    def __call__(
        self, roots: np.ndarray, V: np.ndarray, theta1: float,
    ) -> list[MRTriggerRecord]:
        """Check every armed pair for a −→+ derivative flip.

        Returns
        -------
        records : list[MRTriggerRecord]
            One record per pair whose own derivative flipped from negative
            to positive between the previous call and this one, with its
            distance below the threshold.  Empty when nothing flipped.
            All records of one call share the interval
            ``(self._prev_theta, theta1)`` — the stop happens at the first
            step where ANY pair flips, so every simultaneous flip is caught
            in the same record set.

            On a non-empty return the ``self._prev_*`` snapshot still holds
            the interval's *left* endpoint state (not overwritten).
        """
        n = len(roots)
        pairs, dist, deriv = _pairwise_dist_deriv(roots, V)

        prev = self._prev_deriv
        if prev is not None and prev.shape != (n, n):
            # Root count changed between calls (defensive): stale state.
            prev = None

        new_prev = np.zeros((n, n))
        records: list[MRTriggerRecord] = []
        for k in range(len(pairs)):
            i, j = int(pairs[k, 0]), int(pairs[k, 1])
            d = float(dist[k])
            if d >= self._min_dist_threshold:
                continue  # unarmed this row → per-pair reset (stays 0)
            dv = float(deriv[k])
            new_prev[i, j] = new_prev[j, i] = dv
            if prev is not None and prev[i, j] < 0.0 and dv > 0.0:
                records.append(MRTriggerRecord(
                    pair=(i, j),
                    theta_lo=float(self._prev_theta),
                    theta_hi=theta1,
                    dist_hi=d,
                    deriv_lo=float(prev[i, j]),
                    deriv_hi=dv,
                ))

        if records:
            # Early-return WITHOUT overwriting the snapshot: the interval's
            # left endpoint is the previous step's state, which the caller
            # reads.
            self._prev_deriv = new_prev
            return records

        if not np.any(new_prev):
            # Nothing armed → full reset (the old threshold reset).
            self._reset()
            return []

        self._prev_deriv = new_prev
        self._prev_theta = theta1
        self._prev_roots = roots
        self._prev_V = V
        return []

    def _reset(self) -> None:
        self._prev_deriv = None
        self._prev_theta = None
        self._prev_roots = None
        self._prev_V = None


@live_defaults(cluster_tol="continuation.multiple_roots:CLUSTER_TOL")
def detect_cluster(
    roots: np.ndarray,
    *,
    cluster_tol: Optional[float] = None,
) -> list[tuple[int, ...]]:
    """Find all finite-root clusters as connected components of the proximity graph.

    Two roots are connected when their chordal distance is below
    *cluster_tol*.  Clusters are the connected components of this graph,
    computed via :func:`scipy.sparse.csgraph.connected_components`.

    Roots classified as singular by :func:`_is_singular_root` — exact
    ``0`` / ``∞`` padding roots and roots next to them — are excluded from
    the proximity graph.  The MR solver only models finite branch points;
    a 0/∞ padding root would otherwise cluster with a nearby finite root
    (e.g. ``[0, 1e-5]`` or ``[∞, 2e4]``) and drag it into a fake MR whose
    snapping step corrupts the finite root.  Excluding them keeps 0/∞-side
    degeneracies out of the MR flow entirely.

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
        Each tuple holds the sorted ORIGINAL indices of one finite-root
        cluster.  Empty list when all roots are well-separated, or when
        fewer than two finite (non-singular) roots exist.
    """
    roots = np.asarray(roots)
    finite_positions = np.flatnonzero(
        [not _is_singular_root(r) for r in roots]
    )
    if finite_positions.size < 2:
        return []

    finite_roots = roots[finite_positions]
    r3 = to_sphere_r3(finite_roots)
    pw = cost_from_sphere_r3(r3, r3)

    adj = pw < cluster_tol
    np.fill_diagonal(adj, False)

    n_components, labels = connected_components(adj, directed=False)

    clusters: list[tuple[int, ...]] = []
    for label in range(n_components):
        local_indices = np.where(labels == label)[0]
        if len(local_indices) >= 2:
            clusters.append(tuple(finite_positions[local_indices]))

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


def _pairwise_dist_deriv(
    roots: np.ndarray,
    V: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """All pairwise distances and their θ₁-derivatives, vectorized.

    Computes for every ``i < j`` pair of tracks
    ``dist = |β_i − β_j|`` and
    ``deriv = d|β_i − β_j|²/dθ₁ = 2 Re[(β̇_i − β̇_j)·(β_i − β_j)†]``
    with ``β̇_k = V_k·β_k`` — one broadcasting block instead of a Python
    double loop, returning ALL pairs so callers can track them per pair
    (not just the closest one).

    Roots whose tangent is nan (the 0/∞ padding roots held fixed by the
    integrator) or inf (divergent, at a branch point) are excluded from
    the pair list: a nan/inf derivative would poison per-pair sign state
    (``nan < 0`` is False forever after).  The point trigger owns the
    branch-point regime.

    Returns
    -------
    pairs : np.ndarray (m, 2) int
        Track index pairs ``(i, j)`` with ``i < j``, row-major order
        (the same order — and tie-breaking — as the retired scalar loop).
    dist : np.ndarray (m,)
        Euclidean distances ``|β_i − β_j|``.
    deriv : np.ndarray (m,)
        Raw ``d|β_i − β_j|²/dθ₁`` of each pair.
    """
    roots = np.asarray(roots)
    n = len(roots)
    if n < 2:
        return np.empty((0, 2), dtype=int), np.empty(0), np.empty(0)

    finite = np.isfinite(V.real) & np.isfinite(V.imag)
    # nan/inf tangents make their matrix rows nan on purpose (masked out
    # below); silence the expected invalid-multiply warnings.
    with np.errstate(invalid="ignore"):
        diff = roots[:, None] - roots[None, :]
        dist = np.abs(diff)
        beta_dot = V * roots
        beta_dot_diff = beta_dot[:, None] - beta_dot[None, :]
        deriv = 2.0 * np.real(beta_dot_diff * np.conj(diff))

    iu, ju = np.triu_indices(n, k=1)
    ok = finite[iu] & finite[ju]
    pairs = np.stack([iu[ok], ju[ok]], axis=1)
    return pairs, dist[iu, ju][ok], deriv[iu, ju][ok]


def _closest_pair_deriv(
    roots: np.ndarray,
    V: np.ndarray,
) -> tuple[float, float, tuple[int, int]]:
    """Compute min Euclidean distance among roots, its θ₁-derivative,
    and the pair achieving the minimum.

    Thin wrapper over :func:`_pairwise_dist_deriv`, kept for callers that
    need only the closest pair.  Tie-breaking matches the retired scalar
    loop: the row-major first minimum wins.

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
    pairs, dist, deriv = _pairwise_dist_deriv(roots, V)
    if len(dist) == 0:
        return float(np.inf), 0.0, (-1, -1)
    k = int(np.argmin(dist))
    return (float(dist[k]), float(deriv[k]),
            (int(pairs[k, 0]), int(pairs[k, 1])))


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

    If Brent fails, solve ``f = ∂f/∂β₂ = 0`` from the midpoint of its
    last enclosing bracket.  A failed initial sign check leaves the
    original interval intact.  Cluster verification remains with the
    caller; failure of the point solver propagates to the caller too.

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
        If given, the root pair whose own distance derivative drives the
        root search (per-pair g(θ), continuous across closest-pair
        identity changes — with several pairs degenerating together the
        argmin flickers between them and an argmin-based g is
        discontinuous exactly at the flip).  When None (default), the
        closest pair's derivative is used — useful when the identity of
        the merging pair is unknown ahead of time.

    Returns
    -------
    theta1_mr : float
        The θ₁ value where the derivative crosses zero, i.e. the
        multiple-root location (not wrapped to [0, 2π)).
    """
    if theta1_right < theta1_left:
        theta1_right += TWO_PI

    # scipy's Brent interface does not expose its bracket on failure.
    # Retain the enclosing sign bracket, not merely its last two trials,
    # which may both lie on the same side of the candidate.
    bracket = [theta1_left, theta1_right]
    bracket_values = [None, None]

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
        if min_pair is not None:
            deriv = _pair_distance_deriv(roots, V_list, min_pair)
        else:
            _, deriv, _ = _closest_pair_deriv(roots, V_list)
        if np.isfinite(deriv):
            if theta1 == bracket[0]:
                bracket_values[0] = deriv
            elif theta1 == bracket[1]:
                bracket_values[1] = deriv
            elif (bracket[0] < theta1 < bracket[1]
                  and all(v is not None for v in bracket_values)
                  and np.signbit(bracket_values[0]) != np.signbit(bracket_values[1])):
                side = 0 if np.signbit(deriv) == np.signbit(bracket_values[0]) else 1
                bracket[side] = theta1
                bracket_values[side] = deriv
        return deriv

    # Tangent at the reference (right) endpoint — reused across Brent trials.
    V_ref, _ = compute_tangent(
        poly, E_ref, exp(mu1 + 1j * theta1_right), roots_ref,
    )

    try:
        theta1_mr = optimize.brentq(_compute_deriv, theta1_left, theta1_right)
    except (ValueError, RuntimeError):
        theta_mid = 0.5 * (bracket[0] + bracket[1])
        beta1_mid = exp(mu1 + 1j * theta_mid)
        roots_mid = poly.solve_roots_1d((0, 1), (E_ref, beta1_mid), (2,))
        # The triggering pair may contain a spectator of the actual MR.
        # Seed from the polynomial derivative, as the point-trigger path
        # does, without carrying that pair's ambiguous branch labels.
        finite_roots = [r for r in roots_mid if not _is_singular_root(r)]
        if not finite_roots:
            raise
        beta2_mid = min(
            finite_roots,
            key=lambda r: abs(poly.eval_partials((E_ref, beta1_mid, r))[2]),
        )
        beta1_mr, _ = solve_multiple_roots_iterative(
            poly, E_ref, beta1_mid, beta2_mid,
        )
        theta1_mr = float(np.angle(beta1_mr / exp(1j * theta_mid)) + theta_mid)
        # An unconstrained point solve can reach a different MR.  Such a
        # point cannot split the segment belonging to this interval.
        if not theta1_left <= theta1_mr <= theta1_right:
            raise ValueError(
                f"MR point fallback left the trigger interval: "
                f"theta1={theta1_mr}, interval=({theta1_left}, {theta1_right})"
            )

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

