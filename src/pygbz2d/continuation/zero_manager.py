"""
ZeroManager — polynomial β₂-root topology manager over θ₁ ∈ [0, 2π).

Orchestrates ``integrate_segment`` and ``refine_multiple_root_theta`` to
build a picture of the root manifold at fixed (E, μ₁): a list of refined
multiple roots and a list of per-segment curve data.

Each call to ``integrate_segment`` produces one SegmentData; the loop
stops at multiple roots, refines them, records them, and restarts.
"""

from __future__ import annotations

import numpy as np
from enum import Enum
from math import pi
from cmath import exp
from dataclasses import dataclass
from typing import NamedTuple, Optional
import warnings


# Run-level knobs of the ZeroManager loop.  The MR restart trio below is
# ONE formula's parameters: effective restart distance past an MR =
# max(MR_JUMP, MR_RESTART_FACTOR_H0*min_dtheta/h0,
#     MR_RESTART_FACTOR_ABS*min_dtheta) — keep the three together.
H0: float = 0.1
MIN_DTHETA: float = 1e-10
MR_JUMP: float = 1e-6
MR_RESTART_FACTOR_H0: float = 10.0
MR_RESTART_FACTOR_ABS: float = 100.0
#: Same-MR re-detection retry: each retry multiplies the restart jump by
#: this factor (bounded by MR_REDETECT_MAX_RETRIES).
MR_REDETECT_RETRY_FACTOR: float = 2.0
MR_REDETECT_MAX_RETRIES: int = 4
#: Refined MR θ₁ closer than this to the previous one = no forward progress.
MR_STUCK_TOL: float = 1e-12
#: |θ − 2π| below which a boundary MR is pinned to exactly 2π.
BOUNDARY_THETA_TOL: float = 1e-6
#: Warn when the iterative MR solver's θ₁ drifts more than this.
MR_GAUGE_TOL: float = 0.1
#: Hard cap on segment count (runaway protection).
MAX_SEGMENTS: int = 10000
#: Dense sampling between MR events found inside ONE trigger bracket:
#: maximum θ₁ spacing of the regular rows sampled between two
#: consecutive events.
MR_DENSE_MAX_STEP: float = 1e-4
#: Dense sampling between bracket events: minimum number of interior
#: sample rows between two consecutive events.
MR_DENSE_MIN_SAMPLES: int = 8
from ..core import (
    TWO_PI,
    CharPoly,
    hungarian_match_indices,
    sort_by_root_abs,
)
from .arclength import (
    arclength_step,
    compute_tangent,
    predict_roots_hermite,
    StepControl,
)

from . import multiple_roots
from .multiple_roots import (
    MultipleRootInfo,
    MRTriggerRecord,
    multiple_root_point_trigger,
    MultipleRootIntervalTrigger,
    detect_cluster,
    snap_clusters_to_mean,
    solve_multiple_roots_in_interval,
    solve_multiple_roots_iterative,
)



# ---------------------------------------------------------------------------
# Segment integration
# ---------------------------------------------------------------------------


class StopReason(Enum):
    completed = "completed"
    multiple_root_encountered = "multiple_root_encountered"
    multiple_root_in_interval = "multiple_root_in_interval"


class _MREndpoint(NamedTuple):
    """One endpoint of the bracket straddling a multiple root.

    Carries the (θ₁, track-ordered roots, tangent V) at that endpoint so the
    MR refinement can build a derivative-based prediction anchor without
    re-solving roots or re-deriving the tangent there.
    """
    theta: float
    roots: np.ndarray
    V: np.ndarray


class SegmentResult(NamedTuple):
    """Result of integrating one curve segment until a stop condition.

    ``stop_reason`` is one of:

    * ``'completed'`` — reached *theta_end* without hitting a multiple root.
    * ``'multiple_root_encountered'`` — step size collapsed or integrator
      rejected; a multiple root is at or very near *mr_approx_theta*.
      ZeroManager should refine directly.
    * ``'multiple_root_in_interval'`` — one or more root pairs had their
      closest-pair distance derivative flip sign, so each flipped pair
      has a local minimum (multiple root) in
      [*mr_triggers[k].theta_lo*, *mr_triggers[k].theta_hi*].
      ZeroManager should bisect each pair's interval to locate its MR,
      then refine.

    For both MR stop reasons, *mr_ref* is the MR-adjacent regular endpoint
    (the segment's last accepted row) — used as the single-endpoint
    prediction anchor for ``multiple_root_encountered``.  For
    ``multiple_root_in_interval``, *mr_ref2* is the trigger interval's
    other endpoint so the caller can build a two-endpoint cubic Hermite
    anchor straddling the bracketed MRs.
    """
    theta1_arr: np.ndarray
    tracked_roots: np.ndarray
    abs_argsort: np.ndarray
    stop_reason: StopReason
    mr_approx_theta: float
    n_steps: int
    n_rejected: int
    # MR-adjacent regular endpoint (always set for MR stop reasons).
    mr_ref: Optional[_MREndpoint] = None
    # Other endpoint of the trigger interval (only for in_interval).
    mr_ref2: Optional[_MREndpoint] = None
    # Valid only for 'multiple_root_in_interval': one record per root pair
    # whose own distance derivative flipped −→+ at this stop.  All records
    # share the bracket (theta_lo, theta_hi) = (mr_ref2.theta, mr_ref.theta).
    mr_triggers: tuple[MRTriggerRecord, ...] = ()


def _abs_argsort(roots_2d: np.ndarray) -> np.ndarray:
    """Per-row argsort by |β₂| (column j → |β₂ⱼ| ascending)."""
    return np.argsort(np.abs(roots_2d), axis=1)


def integrate_segment(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    theta_start: float,
    roots_start: np.ndarray,
    theta_end: float,
    *,
    h0: Optional[float] = None,
    ctrl: Optional[StepControl] = None,
    min_dtheta: Optional[float] = None,
    min_dist_threshold: Optional[float] = None,
) -> SegmentResult:
    """Integrate β₂ roots from *theta_start* toward *theta_end*.

    Stops when *theta_end* is reached (``stop_reason=StopReason.completed``) or when
    a multiple root is detected.  Two trigger mechanisms run after each
    accepted step:

    * ``'multiple_root_encountered'`` — step-size collapse (point trigger).
    * ``'multiple_root_in_interval'`` — closest-pair distance derivative
      sign flip (interval trigger).

    *roots_start* must be track-ordered (column j = physical root j) and
    is included as row 0 of the returned ``tracked_roots``.  Each
    subsequent row is Hungarian-matched to the previous one to maintain
    track continuity within the segment.
    """
    if ctrl is None:
        ctrl = StepControl()
    if h0 is None:
        h0 = H0
    if min_dtheta is None:
        min_dtheta = MIN_DTHETA
    if min_dist_threshold is None:
        min_dist_threshold = multiple_roots.MIN_DIST_THRESHOLD
    h = h0
    theta1 = theta_start
    roots = roots_start

    n_roots = len(roots_start)

    theta1_list: list[float] = [theta1]
    tracked_list: list[np.ndarray] = [roots.copy()]

    n_accepted = 0
    n_rejected = 0

    # Interval trigger: tracks closest-pair distance derivative sign
    # across steps.  Encapsulates all internal state and computation.
    interval_trigger = MultipleRootIntervalTrigger(min_dist_threshold)

    def _finish(stop_reason: StopReason, mr_approx_theta: float,
                mr_triggers: tuple[MRTriggerRecord, ...] = (),
                mr_ref: Optional[_MREndpoint] = None,
                mr_ref2: Optional[_MREndpoint] = None,
                ) -> SegmentResult:
        # Assemble the SegmentResult from the accumulated lists.  Shared by
        # every exit path so the (theta1_arr, tracked_roots, abs_argsort)
        # assembly lives in one place instead of being copied four times.
        tracked_arr = np.array(tracked_list)
        return SegmentResult(
            theta1_arr=np.array(theta1_list),
            tracked_roots=tracked_arr,
            abs_argsort=(_abs_argsort(tracked_arr) if len(tracked_list)
                         else np.empty((0, n_roots), dtype=int)),
            stop_reason=stop_reason,
            mr_approx_theta=mr_approx_theta,
            n_steps=n_accepted,
            n_rejected=n_rejected,
            mr_ref=mr_ref,
            mr_ref2=mr_ref2,
            mr_triggers=tuple(mr_triggers),
        )

    while theta1 < theta_end:
        result = arclength_step(
            poly, E_ref, mu1, theta1, roots, h,
            ctrl=ctrl,
        )

        if result.accepted:
            dtheta = abs(result.theta1_new - theta1)

            # ---- Point trigger: step-size collapse ----
            # We are AT a multiple root (tangent diverged → step collapsed).
            if multiple_root_point_trigger(dtheta, min_dtheta=min_dtheta):
                # Single-endpoint anchor: the segment's last accepted row
                # (theta1, roots, V=result.V) sits next to the collapsed MR.
                ref = _MREndpoint(theta1, roots.copy(), result.V.copy())
                return _finish(StopReason.multiple_root_encountered, theta1,
                               mr_ref=ref)

            # ---- Interval trigger: derivative sign flip ----
            # One or more pairs went from shrinking to growing distance —
            # each flipped pair has a local minimum (multiple root) inside
            # [prev_theta, theta1].  Catches "accidental" multiple roots
            # where the tangent stays well-conditioned and step size never
            # collapses; per-pair tracking also catches simultaneous
            # degeneracies whose closest-pair identity flickers row by row.
            records = interval_trigger(roots, result.V, theta1)
            if records:
                # Two-endpoint anchor: the trigger's prev state (interval
                # start) and the current accepted row (interval end) bracket
                # the MRs — both endpoints share the segment's track frame,
                # so a cubic Hermite between them is well-posed.
                ref2 = _MREndpoint(interval_trigger._prev_theta,
                                   interval_trigger._prev_roots.copy(),
                                   interval_trigger._prev_V.copy())
                ref = _MREndpoint(theta1, roots.copy(), result.V.copy())
                return _finish(StopReason.multiple_root_in_interval, theta1,
                               mr_triggers=tuple(records),
                               mr_ref=ref, mr_ref2=ref2)

            theta1_new = result.theta1_new
            h = result.h_new
            # arclength_step already returned roots_new in track order,
            # anchored on the tangent prediction (see arclength_step).  No
            # second matching here — re-matching on the raw old→new chordal
            # distance swaps two near-degenerate tracks at a closest
            # approach and fabricates a spurious |b2|=1 crossing.
            roots_new = result.roots_new

            theta1_list.append(theta1_new)
            tracked_list.append(roots_new)

            theta1 = theta1_new
            roots = roots_new
            n_accepted += 1

        else:
            # max_iter exhausted without acceptable step → near singularity.
            n_rejected += 1
            ref = _MREndpoint(theta1, roots.copy(), result.V.copy())
            return _finish(StopReason.multiple_root_encountered, theta1,
                           mr_ref=ref)

    # Completed without hitting a multiple root.
    return _finish(StopReason.completed, float('nan'))



def _normalize_theta(theta: float) -> float:
    """Normalize to [0, 2π)."""
    return float(theta % (TWO_PI))


def _same_mr_cluster(
    cluster_a: list[tuple[int, ...]],
    cluster_b: list[tuple[int, ...]],
) -> bool:
    """Whether two MR cluster-index lists describe the same cluster set.

    The comparison is order-independent (both the list of clusters and the
    column tuples inside each cluster), because the track order on a
    re-detection may be permuted by the restart matching.
    """
    norm_a = frozenset(frozenset(c) for c in cluster_a)
    norm_b = frozenset(frozenset(c) for c in cluster_b)
    return norm_a == norm_b


def _chordal_dist(a: complex, b: complex) -> float:
    """Chordal distance on the Riemann sphere — the same metric as
    detect_cluster's proximity graph, used by the zero-coincidence test
    that groups trigger candidates into events."""
    return abs(a - b) / np.sqrt((1.0 + abs(a) ** 2) * (1.0 + abs(b) ** 2))


# ---------------------------------------------------------------------------
# SegmentData
# ---------------------------------------------------------------------------

@dataclass
class SegmentData:
    """Curve segment between two consecutive multiple roots (or circle ends).

    Attributes
    ----------
    theta1_arr : np.ndarray (float, shape (N,))
        Local θ₁ mesh, monotonic.  Includes the MR at the right boundary.
    tracked_roots : np.ndarray (complex, shape (N, K))
        Track-ordered β₂ roots.  Column j is consistent within this segment.
    abs_argsort : np.ndarray (int, shape (N, K))
        Per-row argsort by |β₂|.
    left_mr : int
        MR index at the left boundary.  Indexing convention (see
        ``ZeroManager.run``): ``mr == -1`` is the ONLY "no MR" case — it marks
        the θ₁ = 0 / 2π circle seam (segment 0's left, last segment's right).
        Any ``mr >= 0`` is a genuine MR listed in ``multiple_roots[mr]``.
        When ``has_boundary_mr`` is True the boundary MR at θ₁ = 0 occupies
        index 0, so interior MRs start at index 1; **when ``has_boundary_mr``
        is False there is no boundary MR, and interior MRs start at index 0**
        — i.e. ``multiple_roots[0]`` is then the first interior MR, NOT the
        seam.  Do not treat ``mr == 0`` as the seam without checking
        ``has_boundary_mr``; the seam is identified solely by ``mr < 0``,
        with no exceptions.
    right_mr : int
        MR index at the right boundary.  Same convention as ``left_mr``:
        ``-1`` is the θ₁ = 2π / 0 circle seam; any ``>= 0`` is a real MR.
    tangents : np.ndarray (complex, shape (N, K))
        Per-mesh-row analytic tangent V_j = d(ln β₂ⱼ)/dθ₁, always built by
        ``ZeroManager.run``.  Sentinel semantics follow :func:`compute_tangent`:
        ``nan`` = undefined tangent of a 0/∞ padding root, ``inf`` = divergent
        tangent at a multiple root.  MR cluster tangents are set to ``inf``
        by hand during segment finalization (the snapped degenerate roots
        make the analytic value unreliable); only a finite value is a real
        tangent (0 included — the root is stationary).
    """
    theta1_arr: np.ndarray
    tracked_roots: np.ndarray
    abs_argsort: np.ndarray
    left_mr: int
    right_mr: int
    tangents: np.ndarray

    def __post_init__(self):
        n = len(self.theta1_arr)
        if self.tracked_roots.shape[0] != n:
            raise ValueError(
                f"theta1_arr length {n} != tracked_roots rows "
                f"{self.tracked_roots.shape[0]}"
            )
        if self.tangents.shape != self.tracked_roots.shape:
            raise ValueError(
                f"tangents shape {self.tangents.shape} != tracked_roots "
                f"shape {self.tracked_roots.shape}"
            )


@dataclass
class _PendingSeg:
    """A segment held back because its right-end MR was a false positive.

    Its data is merged into the next segment instead of being appended with
    a sentinel ``right_mr``.  Carries ``left_mr`` so the merge preserves it.
    Replaces the old ``SegmentData(right_mr=-2)`` + post-loop patch.
    """
    theta1_arr: np.ndarray
    tracked_roots: np.ndarray
    left_mr: int


class _MREvent(NamedTuple):
    """One materialized multiple-root event inside a trigger bracket.

    ``roots`` keeps the raw (unsnapped) track-frame roots at *theta1* so
    later candidates can run the zero-coincidence test against them —
    the snapped row would trivially pass it for its own clusters.
    """

    theta1: float
    roots: np.ndarray            # unsnapped track-frame roots at theta1
    roots_snapped: np.ndarray    # cluster-snapped closing row
    cluster: list[tuple[int, ...]]
    cluster_stds: list[float]


class _MRRefine(NamedTuple):
    """Outcome of one ``_refine_mr`` pass.

    ``events`` is θ-ascending and empty when every triggered "MR" was a
    false positive (caller holds the segment back as pending).  ``dense``
    holds, for each pair of consecutive events, the ``(theta1_arr,
    tracked_roots)`` rows sampled densely between them.
    """
    events: list[_MREvent]
    dense: list[tuple[np.ndarray, np.ndarray]]
    theta_temp: float
    roots_temp: np.ndarray


# ---------------------------------------------------------------------------
# ZeroManager
# ---------------------------------------------------------------------------

class ZeroManager:
    """β₂-root topology over θ₁ ∈ [0, 2π) at fixed (E, μ₁).

    After ``.run()``:
      - ``multiple_roots`` : list[MultipleRootInfo] — each MR's θ₁, cluster
        indices, and modulus-sorted β₂ roots.
      - ``segments`` : list[SegmentData] — curve segments, each with
        ``left_mr`` / ``right_mr`` indices into ``multiple_roots``
        (-1 = the θ₁=0/2π circle seam, the ONLY non-MR boundary).
      - ``boundary_perm``: np.ndarray (int, K) permutation from right boundary
        to left boundary: roots_right[boundary_perm] == roots_left.
      - ``has_boundary_mr``: bool — whether there is a multiple root at θ₁=0.
        This is also the indexing base shift: when True, the boundary MR is
        ``multiple_roots[0]`` and interior MRs start at index 1; when False,
        there is no boundary MR and interior MRs start at index 0.  Either
        way, a segment's ``left_mr``/``right_mr == -1`` always means the seam,
        never ``mr == 0`` — see ``SegmentData`` for the full convention.

    ``run()`` may only be called once per instance — it appends, not resets.
    Create a fresh ``ZeroManager`` for a second run.
    """
    has_boundary_mr: bool
    left_boundary_roots: np.ndarray
    boundary_perm: np.ndarray

    multiple_roots: list[MultipleRootInfo]
    segments: list[SegmentData]
    poly: CharPoly
    E_ref: complex
    mu1: float
    M: int
    K: int
    _cluster_tol: float

    def __init__(self, poly: CharPoly, E_ref: complex, mu1: float):
        self.poly = poly
        self.E_ref = E_ref
        self.mu1 = mu1
        self.M = poly.M
        self.K = poly.M + poly.N

        self.multiple_roots: list[MultipleRootInfo] = []
        self.segments: list[SegmentData] = []
        self._cluster_tol: float = multiple_roots.CLUSTER_TOL
        # run() is single-shot: it appends to ``multiple_roots``/``segments``
        # instead of resetting them, so a second call would mix two topologies
        # into one half-built state.  Guard against that instead of silently
        # returning garbage.
        self._has_run = False

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(
        self,
        *,
        h0: Optional[float] = None,
        ctrl: Optional[StepControl] = None,
        min_dtheta: Optional[float] = None,
        cluster_tol: Optional[float] = None,
        mr_jump: Optional[float] = None,
        verbose: bool = False,
    ) -> None:
        """Execute the full pipeline.

        Parameters
        ----------
        h0 :
            Initial arclength step passed to ``integrate_segment``.
        ctrl :
            Step-size controller (tolerances and RK45 factors) passed to
            ``integrate_segment``.
        min_dtheta :
            Point-trigger threshold for multiple-root detection.
        cluster_tol :
            Chordal-distance threshold for multiple-root detection.
        mr_jump :
            θ₁ step used to jump past a refined multiple root.

        Raises
        ------
        RuntimeError
            If ``run`` was already called on this instance.  Each run needs a
            fresh ``ZeroManager``: a second run would append a new topology on
            top of the first one instead of replacing it.
        """
        if ctrl is None:
            ctrl = StepControl()
        if h0 is None:
            h0 = H0
        if min_dtheta is None:
            min_dtheta = MIN_DTHETA
        if cluster_tol is None:
            cluster_tol = multiple_roots.CLUSTER_TOL
        if mr_jump is None:
            mr_jump = MR_JUMP
        if self._has_run:
            raise RuntimeError(
                "ZeroManager.run() has already been called on this instance; "
                "create a fresh ZeroManager for another run."
            )
        self._has_run = True
        self._cluster_tol = cluster_tol

        # Effective restart distance past an MR: never smaller than the
        # branch-point-safe floors (see MR_RESTART_FACTOR_H0 / _ABS), or the
        # steps right after the restart collapse below min_dtheta and the
        # point trigger re-detects the SAME MR in an infinite ping-pong.
        mr_jump_eff = max(
            mr_jump,
            MR_RESTART_FACTOR_H0 * min_dtheta / h0,
            MR_RESTART_FACTOR_ABS * min_dtheta,
        )

        # ---- Init at θ₁ = 0 ----
        theta = 0.0
        roots = self._solve(theta)  # modulus-sorted
        roots = sort_by_root_abs(roots)

        if verbose:
            print("Initial roots: ", roots)

        cluster = detect_cluster(roots, cluster_tol=cluster_tol)
        self.has_boundary_mr = len(cluster) > 0

        # Enforce exact degeneracy at a boundary MR: snap each cluster's roots
        # to their mean *before* recording left_boundary_roots, so the snapped
        # value propagates consistently into the boundary MR record, the
        # segment's closing row, and the boundary_perm matching below.
        cluster_stds: list[float] = []
        if self.has_boundary_mr:
            roots, cluster_stds = snap_clusters_to_mean(roots, cluster)
        self.left_boundary_roots = roots

        # Segment-joining state across iterations:
        #   left_mr  — index of the MR on the current segment's LEFT boundary
        #              (-1 = none).  Stable across a false-positive MR.
        #   pending  — a segment held back because its right-end MR was a false
        #              positive; merged into the next segment instead of being
        #              appended with a sentinel.  None when nothing is pending.
        #              (Replaces the old right_mr == -2 + post-loop patch.)
        pending: _PendingSeg | None = None

        if self.has_boundary_mr:
            left_mr = 0
            # Add a multiple root at theta1 = 0.
            self.multiple_roots.append(
                MultipleRootInfo(
                    theta1=0.0,
                    cluster_indices=cluster,
                    roots=roots,
                    cluster_stds=tuple(cluster_stds),
                )
            )

            # Reinitialize the solver to avoid MR (mr_jump_eff, not the bare
            # mr_jump: the boundary MR at θ=0 is a branch point like any
            # other — a too-close restart re-triggers the point trigger).
            theta += mr_jump_eff
            # The boundary MR at θ = 0 is a degenerate cluster (snapped above),
            # so the tangent there is undefined — compute_tangent zeros the
            # divergent tracks and the prediction degrades to holding the
            # cluster fixed.  The cluster→split matching is then genuinely
            # ambiguous; this is the one site where the prediction anchor
            # cannot help, but it does not regress the bare match either.
            roots = self._to_track_order(self._solve(theta), theta, 0.0, roots)
        else:
            left_mr = -1

        # Set True by the 'completed' or boundary-MR branches when the right
        # boundary (θ₁ = 2π) is reached.  integrate_segment itself runs to
        # θ_end = 2π or stops at an MR, so the loop below only terminates via
        # those two branches — no θ < 2π guard here, which would duplicate the
        # integrator's own and need a post-loop fallback.
        boundary_perm_set = False

        # Same-MR re-detection retry counter (see the MR branch below).
        redetect_retries = 0

        # ---- Integrate segments ----
        # `for` with a hard cap instead of `while theta < 2π`: the integrator
        # already advances θ to 2π (completed) or an MR; the cap only catches
        # a runaway refine loop that never converges to 2π.
        for _ in range(MAX_SEGMENTS):
            seg = integrate_segment(
                self.poly, self.E_ref, self.mu1,
                theta, roots, TWO_PI,
                h0=h0, ctrl=ctrl, min_dtheta=min_dtheta,
            )
            if verbose:
                print("====== Segment solved =======")
                print("theta1 range: ", (seg.theta1_arr[0], seg.theta1_arr[-1]))
                print("stop_reason: ", seg.stop_reason.value)

            # Common processes: 
            # Attach the left boundary: resume a pending (false-positive)
            # segment, else prepend the left MR, else use seg as-is.
            if pending is not None:
                new_seg_theta1 = np.concatenate([
                    pending.theta1_arr, seg.theta1_arr
                ])
                new_seg_tracked_roots = np.vstack([
                    pending.tracked_roots, seg.tracked_roots
                ])
                pending = None
            elif left_mr >= 0:
                new_seg_theta1 = np.concatenate([
                    [self.multiple_roots[left_mr].theta1],
                    seg.theta1_arr
                ])
                new_seg_tracked_roots = np.vstack([
                    self.multiple_roots[left_mr].roots,
                    seg.tracked_roots
                ])
            else:
                new_seg_theta1 = seg.theta1_arr
                new_seg_tracked_roots = seg.tracked_roots


            if seg.stop_reason == StopReason.completed:
                # Predict roots at θ₁ = 2π from the segment rows bracketing
                # 2π, and match predicted → left_boundary_roots (≡ solve(0)).
                # Anchoring on a smooth continuation of each track (rather
                # than a bare match between two solved sets) avoids swapping
                # near-degenerate tracks at a closest approach.
                predicted_2pi = self._predict_roots_at_2pi(
                    new_seg_theta1, new_seg_tracked_roots,
                )

                # Convention (see ZeroManager docstring): roots_right[boundary_perm] == roots_left.
                self.boundary_perm, perm = self._boundary_perm_from_right(predicted_2pi)
                boundary_perm_set = True

                new_seg_tracked_roots = np.vstack([
                    new_seg_tracked_roots[:-1, :],
                    self.left_boundary_roots[perm]
                ])
                self._append_segment(
                    np.concatenate([new_seg_theta1[:-1], [TWO_PI]]),
                    new_seg_tracked_roots,
                    left_mr,
                    0 if self.has_boundary_mr else -1,
                )
                break

            elif (seg.stop_reason == StopReason.multiple_root_in_interval
                  or seg.stop_reason == StopReason.multiple_root_encountered
            ):
                refined = self._refine_mr(
                    seg, new_seg_theta1, new_seg_tracked_roots,
                    mr_jump, mr_jump_eff, verbose,
                )
                if refined.events:
                    ev0 = refined.events[0]
                    if verbose:
                        for ev in refined.events:
                            print("Find multiple roots. θ₁ = ", ev.theta1,
                                  " Cluster = ", ev.cluster)
                    # Same-MR re-detection guard: when the restart row still
                    # lands inside the previous MR's degenerate neighbourhood
                    # and the same cluster is detected again, do not append a
                    # duplicate MR.  Restart farther away and retry.
                    # (Checked on the FIRST event only — later events of the
                    # same bracket lie strictly past it, hence past the
                    # guard window too.)
                    if left_mr >= 0:
                        prev_mr = self.multiple_roots[left_mr]
                        same_cluster = _same_mr_cluster(
                            ev0.cluster, prev_mr.cluster_indices)
                        if (same_cluster
                                and ev0.theta1 <= prev_mr.theta1 + mr_jump_eff):
                            redetect_retries += 1
                            if redetect_retries > MR_REDETECT_MAX_RETRIES:
                                raise RuntimeError(
                                    f"Same MR re-detected {redetect_retries} "
                                    f"times at θ≈{ev0.theta1!r} (prev "
                                    f"θ={prev_mr.theta1!r}); the restart "
                                    f"distance may still be too small for the "
                                    f"current (h0, min_dtheta) pair."
                                )
                            jump = mr_jump_eff * (
                                MR_REDETECT_RETRY_FACTOR ** redetect_retries)
                            theta = min(prev_mr.theta1 + jump, TWO_PI)
                            roots = self._to_track_order(
                                self._solve(theta), theta,
                                prev_mr.theta1, prev_mr.roots,
                            )
                            continue
                    # Forward-progress guard: the refined MR must lie strictly
                    # past the MR that closed the previous segment.  A θ₁ at or
                    # behind it means the restart landed before the MR again
                    # (ping-pong) — fail loudly instead of looping forever.
                    if (left_mr >= 0 and ev0.theta1
                            <= self.multiple_roots[left_mr].theta1 + MR_STUCK_TOL):
                        raise RuntimeError(
                            f"Multiple-root refinement made no forward "
                            f"progress: refined θ₁={ev0.theta1!r} is not past "
                            f"the previous MR (θ₁="
                            f"{self.multiple_roots[left_mr].theta1!r}). The "
                            f"MR restart distance may be too small for the "
                            f"current (h0, min_dtheta) pair."
                        )
                    # Trim the segment up to the FIRST event and close it
                    # there; each later event closes the densely sampled
                    # inter-event rows instead (see _dense_rows_between).
                    seg_used = new_seg_theta1 < ev0.theta1
                    seg_theta1 = np.concatenate([
                        new_seg_theta1[seg_used],
                        [ev0.theta1]
                    ])
                    seg_tracked_roots = np.vstack([
                        new_seg_tracked_roots[seg_used, :],
                        ev0.roots_snapped
                    ])

                    boundary_reached = False
                    for k, ev in enumerate(refined.events):
                        if k > 0:
                            seg_theta1, seg_tracked_roots = refined.dense[k - 1]
                        if abs(ev.theta1 - TWO_PI) < BOUNDARY_THETA_TOL:
                            if verbose:
                                print("Boundary multiple root detected. right_mr = 0")
                            # Pin the closing row's θ to exactly 2π (the boundary
                            # MR sits within BOUNDARY_THETA_TOL = 1e-6 of it):
                            # leaving the refined θ_mr in the array opens a sliver
                            # gap (θ_mr, 2π) that locate() would refuse to cover.
                            seg_theta1 = np.concatenate([
                                seg_theta1[:-1], [TWO_PI]
                            ])
                            # Same convention as 'completed': predict at 2π and
                            # match onto left_boundary_roots BEFORE appending the
                            # segment.  MR 0 was recorded in the θ=0 modulus-sorted
                            # frame, while this closing row is in the final
                            # segment's track frame; _append_segment needs
                            # boundary_perm to translate MR-cluster columns before
                            # marking their tangents inf.
                            predicted_2pi = self._predict_roots_at_2pi(
                                seg_theta1, seg_tracked_roots,
                            )
                            self.boundary_perm, _ = self._boundary_perm_from_right(
                                predicted_2pi
                            )
                            boundary_perm_set = True
                            self._append_segment(
                                seg_theta1, seg_tracked_roots, left_mr, 0
                            )
                            # The boundary MR at θ₁_mr ≈ 2π is degenerate (the
                            # closing row is the snapped cluster), so the Hermite
                            # prediction degrades toward lerp / hold-fixed — no
                            # regression versus a bare match, and the
                            # cluster↔cluster correspondence stays ambiguous.
                            boundary_reached = True
                            break
                        right_mr = len(self.multiple_roots)
                        # Append the MR record BEFORE _append_segment so the
                        # tangent finalizer can look up its cluster indices
                        # and set the MR row's cluster tangents to inf.
                        self.multiple_roots.append(
                            MultipleRootInfo(
                                theta1=ev.theta1,
                                cluster_indices=ev.cluster,
                                roots=ev.roots_snapped,
                                cluster_stds=tuple(ev.cluster_stds),
                            )
                        )
                        self._append_segment(
                            seg_theta1, seg_tracked_roots, left_mr, right_mr
                        )
                        left_mr = right_mr
                    if boundary_reached:
                        break
                else:
                    if verbose:
                        print(f"No multiple roots found with cluster_tol = "
                              f"{self._cluster_tol: .2e}.")
                    # False positive: hold this segment back as pending instead
                    # of appending it with a sentinel right_mr.
                    pending = _PendingSeg(
                        new_seg_theta1[:-1],
                        new_seg_tracked_roots[:-1, :],
                        left_mr,
                    )

                theta = refined.theta_temp
                roots = refined.roots_temp

            else:
                raise RuntimeError(
                    f"Unknown stop_reason: {seg.stop_reason.value}"
                )
        else:
            # for-loop exhausted without reaching θ₁ = 2π (completed) or a
            # boundary MR — a runaway refine loop.  Surface it loudly rather
            # than silently producing a half-built topology.
            raise RuntimeError(
                f"ZeroManager did not reach θ₁ = 2π within "
                f"{MAX_SEGMENTS} segments; stuck near θ₁ = {theta}"
            )

        # A pending false-positive segment can survive only if the loop body
        # broke on a completed/boundary branch in the same iteration that
        # resolved the pending — but those branches append before breaking, so
        # there is nothing to flush.  Guard defensively anyway: if a pending
        # segment somehow lingers, it has no right MR (the false positive was
        # never resolved into a real one), so close it with right_mr = -1.
        if pending is not None:
            self._append_segment(
                pending.theta1_arr, pending.tracked_roots, pending.left_mr, -1
            )

        # Sanity: completed or boundary-MR must have set boundary_perm.
        # A raise (not assert) so it survives `python -O` — a half-built
        # topology here would otherwise AttributeError downstream on
        # self.boundary_perm.
        if not boundary_perm_set:
            raise RuntimeError(
                "ZeroManager loop exited without setting boundary_perm "
                "(neither completed nor boundary-MR was reached)"
            )

    @property
    def n_multiple_roots(self) -> int:
        return len(self.multiple_roots)

    @property
    def n_segments(self) -> int:
        return len(self.segments)

    # ------------------------------------------------------------------
    # Interpolation & insertion
    # ------------------------------------------------------------------

    def locate(self, theta1: float) -> tuple[int, int]:
        """Find the segment and mesh interval containing *theta1*.

        Returns ``(seg_idx, i)`` such that
        ``segments[seg_idx].theta1_arr[i] <= theta1 < theta1_arr[i+1]``
        (left-closed / right-open; the final mesh point is clamped to the
        last interval).

        θ₁ is normalized to ``[0, 2π)`` first (``_normalize_theta``), so a
        request at exactly 2π wraps to 0 → segment 0, ``i = 0``
        (≡ ``left_boundary_roots``, consistent with the closing-row
        convention).  Raises ``ValueError`` if no segment covers *theta1*
        — surfaces a half-built topology loudly rather than silently.

        Containment is tested with exact ``lo <= θ <= hi``: adjacent segments
        share their endpoint θ *exactly* (the MR θ is stored once in
        ``multiple_roots`` and prepended to the next segment by ``run()``),
        so no floating-point cushion is needed — a boundary θ is caught by
        the first segment whose range contains it.

        *i* is clamped to ``[0, N-2]``; a single-row segment (``N == 1``)
        yields ``i = 0`` and callers that need a real interval
        (``interpolate_roots``, ``insert_solution``) raise on it themselves.
        """
        theta1_wrapped = _normalize_theta(theta1)
        for s_idx, seg in enumerate(self.segments):
            arr = seg.theta1_arr
            lo, hi = arr[0], arr[-1]
            if lo <= theta1_wrapped < hi:
                n = len(arr)
                if n < 2:
                    return s_idx, 0
                i = int(np.searchsorted(arr, theta1_wrapped, side='right') - 1)
                i = max(0, min(i, n - 2))
                return s_idx, i
        raise ValueError(
            f"theta1={theta1} (normalized {theta1_wrapped}) is not covered by "
            f"any of the {len(self.segments)} segments"
        )

    def interpolate_roots(
        self,
        theta1: float,
        seg_idx: Optional[int] = None,
        i: Optional[int] = None,
    ) -> np.ndarray:
        """Cubic-Hermite interpolate the track-ordered β₂ roots at *theta1*.

        Uses the value *and* analytic tangent (``SegmentData.tangents``) at
        both ends of the mesh sub-interval ``[θ_i, θ_{i+1}]`` of
        ``segments[seg_idx]`` — per-track cubic Hermite, the same kernel as
        :func:`predict_roots_hermite`.  Exact for tracks that are cubic in θ₁.

        *seg_idx* / *i* may be passed to skip the θ₁→index search; either
        left ``None`` is resolved via :meth:`locate`.

        **Linear fallback (whole interval):** when the sub-interval touches a
        segment MR boundary — ``i == 0`` with ``left_mr >= 0`` or
        ``i == N-2`` with ``right_mr >= 0`` — ``dβ₂/dθ₁`` diverges on the
        cluster tracks as ``(θ − θ_MR)^{-1/2}`` (square-root branch), so the
        cubic Hermite tangent is unusable.  Linear interpolation
        ``r_a + s·(r_b − r_a)`` is used for *all* tracks instead of a
        per-track fallback: per-track fallback would silently hold the
        cluster fixed while Hermite-interpolating the smooth tracks, hiding
        the kink.  Linear is honest and bounded.

        Raises ``ValueError`` for a single-row segment (no interval to
        interpolate).
        """
        if seg_idx is None or i is None:
            seg_idx, i = self.locate(theta1)
        seg = self.segments[seg_idx]
        arr = seg.theta1_arr
        n = len(arr)
        if n < 2:
            raise ValueError(
                f"segment {seg_idx} has {n} row(s); cannot interpolate"
            )
        theta_a, theta_b = arr[i], arr[i + 1]
        r_a = seg.tracked_roots[i, :]
        r_b = seg.tracked_roots[i + 1, :]
        theta1_wrapped = _normalize_theta(theta1)

        # An interval touching an MR boundary has a divergent dβ₂/dθ₁ on the
        # cluster tracks → linear for everyone.
        touches_mr = (
            (i == 0 and seg.left_mr >= 0)
            or (i == n - 2 and seg.right_mr >= 0)
        )
        if touches_mr:
            s = (theta1_wrapped - theta_a) / (theta_b - theta_a)
            return r_a + s * (r_b - r_a)

        return predict_roots_hermite(
            theta1_wrapped, theta_a, r_a, seg.tangents[i, :],
            theta_b, r_b, seg.tangents[i + 1, :],
        )

    def _interpolate_roots_linear(
        self,
        theta1: float,
        seg_idx: int,
        i: int,
    ) -> np.ndarray:
        """Two-point linear interpolation of the track-ordered β₂ roots.

        The matching anchor used by the SGBZ pairwise crossing refinement:
        simple, bounded, and independent of cubic Hermite.  Kept as a
        separate method so :meth:`solve_at`'s prediction anchor is an
        explicit caller choice.
        """
        seg = self.segments[seg_idx]
        arr = seg.theta1_arr
        n = len(arr)
        if n < 2:
            raise ValueError(
                f"segment {seg_idx} has {n} row(s); cannot interpolate"
            )
        theta_a, theta_b = arr[i], arr[i + 1]
        r_a = seg.tracked_roots[i, :]
        r_b = seg.tracked_roots[i + 1, :]
        theta1_wrapped = _normalize_theta(theta1)
        s = (theta1_wrapped - theta_a) / (theta_b - theta_a)
        return r_a + s * (r_b - r_a)

    def solve_at(
        self,
        theta1: float,
        seg_idx: Optional[int] = None,
        i: Optional[int] = None,
        *,
        interp: str = 'hermite',
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solve + track-match β₂ roots at *theta1* — WITHOUT inserting a row.

        The non-mutating core of :meth:`insert_solution` (which is this
        method plus the mesh insert): solves ``self._solve(theta1)``, reorders
        the ``np.roots``-order result onto the segment's track frame via
        Hungarian matching against ``interp``-chosen prediction anchor and
        computes the analytic tangent ``V``.  *interp* selects the matching
        anchor: ``'hermite'`` (default, :meth:`interpolate_roots`) or
        ``'linear'`` (:meth:`_interpolate_roots_linear`).

        Unlike :meth:`insert_solution` there is no coincidence early-return:
        ``theta1 == theta_a`` solves fine (equivalent to reading row *i* up to
        Hungarian tie-breaking), while ``theta1 == theta_b`` raises the same
        containment ``ValueError`` as any out-of-interval θ — the coincidence
        semantics are handled by ``insert_solution`` *before* calling this.

        Returns ``(roots, V)``; *roots* is track-ordered and *V* is the
        ``compute_tangent`` tangent (always computed as a pure function of
        the solved roots).  Nothing is written into ``SegmentData``.
        """
        if seg_idx is None or i is None:
            seg_idx, i = self.locate(theta1)
        seg = self.segments[seg_idx]
        arr = seg.theta1_arr
        n = len(arr)
        if n < 2:
            raise ValueError(
                f"segment {seg_idx} has {n} row(s); cannot solve"
            )

        theta1_wrapped = _normalize_theta(theta1)
        theta_a, theta_b = arr[i], arr[i + 1]
        if not (theta_a <= theta1_wrapped < theta_b):
            raise ValueError(
                f"theta1={theta1} (normalized {theta1_wrapped}) outside interval "
                f"[{theta_a}, {theta_b}] of segment {seg_idx}, i={i}"
            )

        roots = self._solve(theta1_wrapped)
        # Reorder np.roots-order onto the segment track frame via Hungarian
        # matching against ``interpolate_roots`` as the prediction anchor.
        # Using interpolate_roots rather than the raw _to_track_order Hermite
        # prediction is deliberate: when the bracketing interval is near a
        # multiple root, interpolate_roots falls back to whole-interval
        # linear interpolation (see its docstring), yielding a bounded,
        # honest prediction for ALL tracks including the cluster ones — the
        # per-track Hermite→lerp fallback in _to_track_order would silently
        # hold cluster tracks at one endpoint while Hermite-interpolating the
        # smooth ones, hiding the near-MR divergence.
        if interp == 'linear':
            predicted = self._interpolate_roots_linear(
                theta1_wrapped, seg_idx=seg_idx, i=i,
            )
        elif interp == 'hermite':
            predicted = self.interpolate_roots(
                theta1_wrapped, seg_idx=seg_idx, i=i,
            )
        else:
            raise ValueError(f"unknown interp={interp!r}; use 'hermite' or 'linear'")
        perm = hungarian_match_indices(predicted, roots)
        roots = roots[perm]

        beta1 = exp(self.mu1 + 1j * theta1_wrapped)
        V, _ = compute_tangent(self.poly, self.E_ref, beta1, roots)
        return roots, V

    def insert_solution(
        self,
        theta1: float,
        seg_idx: Optional[int] = None,
        i: Optional[int] = None,
        *,
        interp: str = 'hermite',
    ) -> tuple[int, bool]:
        """Solve β₂ roots at an interior *theta1* and insert the new row.

        :meth:`solve_at` (solve + track-frame matching + tangent) followed by
        the mesh insert: the new row is inserted at mesh position ``i + 1``
        so ``theta1_arr`` stays monotonic; ``abs_argsort`` is refreshed and
        the tangent spliced into ``SegmentData.tangents``.

        Returns ``(insert_at, changed)`` where *insert_at* is the mesh index
        of the (existing or inserted) row and *changed* is ``True`` when a
        new row was actually inserted.  When *theta1* coincides with an
        existing mesh point (``theta_a`` or ``theta_b``), *insert_at* is the
        index of that point, *changed* is ``False``, and no mutation occurs.

        *seg_idx* / *i* follow :meth:`interpolate_roots` (either may be
        ``None`` → resolved via :meth:`locate`).  *interp* selects the
        matching anchor passed to :meth:`solve_at` (``'hermite'`` default,
        or ``'linear'``).

        Raises ``ValueError`` when the segment has fewer than 2 rows or
        *theta1* is outside the specified interval.
        """
        if seg_idx is None or i is None:
            seg_idx, i = self.locate(theta1)
        seg = self.segments[seg_idx]
        arr = seg.theta1_arr
        n = len(arr)
        if n < 2:
            raise ValueError(
                f"segment {seg_idx} has {n} row(s); cannot insert"
            )

        theta1_wrapped = _normalize_theta(theta1)
        theta_a, theta_b = arr[i], arr[i + 1]
        # Coincidence with an existing mesh point: return its index, no mutation.
        if theta1_wrapped == theta_a:
            return i, False
        if theta1_wrapped == theta_b:
            return i + 1, False
        if not (theta_a <= theta1_wrapped < theta_b):
            raise ValueError(
                f"theta1={theta1} (normalized {theta1_wrapped}) outside interval "
                f"[{theta_a}, {theta_b}] of segment {seg_idx}, i={i}"
            )

        # Solve + reorder onto the track frame (the non-mutating core).
        roots, V_new = self.solve_at(
            theta1_wrapped, seg_idx=seg_idx, i=i, interp=interp,
        )

        insert_at = i + 1
        seg.theta1_arr = np.insert(arr, insert_at, theta1_wrapped)
        seg.tracked_roots = np.insert(
            seg.tracked_roots, insert_at, roots, axis=0,
        )
        seg.abs_argsort = _abs_argsort(seg.tracked_roots)
        seg.tangents = np.insert(seg.tangents, insert_at, V_new, axis=0)

        return insert_at, True

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _solve(self, theta1: float) -> np.ndarray:
        """Solve β₂ roots at *theta1*."""
        beta1 = exp(self.mu1 + 1j * _normalize_theta(theta1))
        return self.poly.solve_roots_1d((0, 1), (self.E_ref, beta1), (2,))

    def _to_track_order(
        self,
        roots: np.ndarray,
        theta_target: float,
        ref_theta: float, ref_roots: np.ndarray,
        ref_theta2: float | None = None, ref_roots2: np.ndarray | None = None,
        ref_V: np.ndarray | None = None, ref_V2: np.ndarray | None = None,
    ) -> np.ndarray:
        """Reorder *roots* (solved at *theta_target*, np.roots order) onto the
        track frame of the reference endpoint(s).

        The anchor for the Hungarian match is a derivative-based prediction,
        not a bare match between two solved sets — the latter swaps two
        near-degenerate tracks at a closest approach, exactly the failure
        ``arclength_step`` avoids by matching *predicted → solved*.  Here we
        reuse the same principle via :func:`predict_roots_hermite`.

        Two-endpoint (``ref_theta2``/``ref_roots2`` given, endpoints share a
        track frame, target inside the interval) → cubic Hermite.
        One-endpoint → tangent extrapolation.  Singular tracks fall back to
        lerp / hold-fixed inside :func:`predict_roots_hermite`.

        ``ref_V`` / ``ref_V2`` let the caller pass pre-computed tangents
        (e.g. from a :class:`_MREndpoint`); when omitted the tangent is
        recomputed here.
        """
        V = (ref_V if ref_V is not None
             else compute_tangent(self.poly, self.E_ref,
                                  exp(self.mu1 + 1j * ref_theta), ref_roots)[0])
        if ref_theta2 is not None and ref_roots2 is not None:
            V2 = (ref_V2 if ref_V2 is not None
                  else compute_tangent(self.poly, self.E_ref,
                                       exp(self.mu1 + 1j * ref_theta2),
                                       ref_roots2)[0])
            predicted = predict_roots_hermite(
                theta_target, ref_theta, ref_roots, V,
                ref_theta2, ref_roots2, V2,
            )
        else:
            predicted = predict_roots_hermite(theta_target, ref_theta, ref_roots, V)
        perm = hungarian_match_indices(predicted, roots)
        return roots[perm]

    def _boundary_perm_from_right(
        self, roots_right: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Permutation mapping the segment's right boundary (track frame,
        e.g. a Hermite prediction at θ₁ = 2π) onto ``left_boundary_roots``.

        Returns ``(boundary_perm, perm)`` where:

        * ``perm = hungarian_match_indices(roots_right, left_boundary_roots)``
          — reorder the left boundary into the segment frame with
          ``left_boundary_roots[perm]`` (the closing row).
        * ``boundary_perm`` — the inverse of ``perm``, satisfying the class
          convention ``roots_right[boundary_perm] == roots_left``.
        """
        perm = hungarian_match_indices(roots_right, self.left_boundary_roots)
        boundary_perm = np.zeros(self.K, dtype=int)
        boundary_perm[perm] = np.arange(self.K)
        return boundary_perm, perm

    def _predict_roots_at_2pi(
        self,
        theta1_arr: np.ndarray,
        tracked_roots: np.ndarray,
    ) -> np.ndarray:
        """Predict track-ordered β₂ roots at θ₁ = 2π from segment rows.

        Locates the last row with θ₁ < 2π (``idx``).  If the next row exists
        it lies past 2π, so the pair brackets 2π and a two-endpoint cubic
        Hermite (value + tangent at both rows, same track frame) is exact for
        cubic-in-θ tracks.  Otherwise (segment ends at or before 2π, or a
        single-row segment) falls back to single-end tangent extrapolation
        from the nearest row — the same prediction the arclength integrator
        uses.  Singular tracks degrade to lerp / hold-fixed inside
        :func:`predict_roots_hermite`.
        """
        below = theta1_arr < TWO_PI
        idx = int(np.flatnonzero(below)[-1]) if below.any() else -1

        if idx >= 0 and idx < len(theta1_arr) - 1:
            # Bracket: rows idx (< 2π) and idx+1 (≥ 2π) straddle 2π.
            theta_a = theta1_arr[idx]
            theta_b = theta1_arr[idx + 1]
            roots_a = tracked_roots[idx, :]
            roots_b = tracked_roots[idx + 1, :]
            Va = compute_tangent(
                self.poly, self.E_ref,
                exp(self.mu1 + 1j * theta_a), roots_a,
            )[0]
            Vb = compute_tangent(
                self.poly, self.E_ref,
                exp(self.mu1 + 1j * theta_b), roots_b,
            )[0]
            return predict_roots_hermite(
                TWO_PI, theta_a, roots_a, Va, theta_b, roots_b, Vb,
            )

        # No bracketing pair: single-end extrapolate from the row nearest 2π.
        src = idx if idx >= 0 else 0
        theta_s = theta1_arr[src]
        roots_s = tracked_roots[src, :]
        Vs = compute_tangent(
            self.poly, self.E_ref, exp(self.mu1 + 1j * theta_s), roots_s,
        )[0]
        return predict_roots_hermite(TWO_PI, theta_s, roots_s, Vs)

    # ------------------------------------------------------------------
    # run() helpers
    # ------------------------------------------------------------------

    def _compute_tangents(
        self, theta1_arr: np.ndarray, tracked_roots: np.ndarray,
    ) -> np.ndarray:
        """Per-mesh-row analytic tangent V for every row of a finalized segment.

        ``tangents[i] = compute_tangent(...)[0]`` evaluated at
        ``theta1_arr[i]`` on the track-ordered roots.  Recomputed over the
        finalized arrays rather than collected during ``integrate_segment``
        because ``run()`` transforms the mesh/roots after integration
        (pending-segment merges, boundary-MR row prepend, the
        ``completed``-branch closing-row replacement); mirroring those
        transforms on a collected-V array is error-prone, while a one-pass
        recompute over the final arrays is simple and correct.  Cost is
        ``O(N·K)`` ``eval_partials`` calls per segment — negligible vs the
        integration itself.

        Sentinel semantics: ``nan`` marks the undefined tangent of a 0/∞
        padding root (from :func:`compute_tangent`).  ``inf`` marks a
        divergent multiple-root tangent — ``compute_tangent`` may produce it,
        and ``_mark_mr_tangents_inf`` additionally sets it MANUALLY on MR
        cluster columns after ``snap_clusters_to_mean``.  Interpolation sites
        use these sentinels to fall back to linear / hold fixed.
        """
        n = len(theta1_arr)
        tangents = np.zeros((n, self.K), dtype=complex)
        for i in range(n):
            theta = theta1_arr[i]
            beta1 = exp(self.mu1 + 1j * _normalize_theta(theta))
            V, _ = compute_tangent(
                self.poly, self.E_ref, beta1, tracked_roots[i, :]
            )
            tangents[i, :] = V
        return tangents

    def _mark_mr_tangents_inf(
        self,
        tangents: np.ndarray,
        left_mr: int,
        right_mr: int,
    ) -> np.ndarray:
        """Set MR cluster tangents to ``inf`` on both boundary rows.

        ``snap_clusters_to_mean`` enforces exact root degeneracy, so the
        analytic derivative on cluster tracks is no longer trustworthy (at a
        true multiple root it diverges).  The ``inf`` here is therefore
        MANUALLY SET during segment finalization, not returned by
        ``compute_tangent`` — do not remove it under the assumption that it
        was computed analytically.

        Frame rule: every interior MR record shares the track frame of its
        boundary row.  The one exception is the θ=0 boundary MR reused as the
        last segment's RIGHT boundary at θ=2π: its ``cluster_indices`` index
        the θ=0 modulus-sorted ``left_boundary_roots``, while the closing row
        is in that segment's track frame.  The class convention
        ``roots_right[boundary_perm] == roots_left`` means left column ``k``
        is right column ``boundary_perm[k]`` there, so the columns must be
        translated before being marked.
        """
        n = tangents.shape[0]
        for row, mr_idx in ((0, left_mr), (n - 1, right_mr)):
            if mr_idx < 0:
                continue
            seam_translate = (
                row == n - 1
                and mr_idx == 0
                and self.has_boundary_mr
            )
            if seam_translate and not hasattr(self, "boundary_perm"):
                raise RuntimeError(
                    "cannot mark the θ=2π boundary MR tangents before "
                    "boundary_perm is set; compute the right→left root "
                    "monodromy before appending the final segment"
                )

            for cluster in self.multiple_roots[mr_idx].cluster_indices:
                for col in cluster:
                    out_col = (int(self.boundary_perm[int(col)])
                               if seam_translate else int(col))
                    tangents[row, out_col] = np.inf + 0j
        return tangents

    def _append_segment(
        self,
        theta1_arr: np.ndarray,
        tracked_roots: np.ndarray,
        left_mr: int,
        right_mr: int,
    ) -> None:
        """Append a SegmentData with abs_argsort and tangents precomputed."""
        tangents = self._compute_tangents(theta1_arr, tracked_roots)
        tangents = self._mark_mr_tangents_inf(tangents, left_mr, right_mr)
        self.segments.append(
            SegmentData(
                theta1_arr=theta1_arr,
                tracked_roots=tracked_roots,
                abs_argsort=_abs_argsort(tracked_roots),
                left_mr=left_mr,
                right_mr=right_mr,
                tangents=tangents,
            )
        )

    def _refine_mr(
        self,
        seg: SegmentResult,
        new_seg_theta1: np.ndarray,
        new_seg_tracked_roots: np.ndarray,
        mr_jump: float,
        mr_jump_eff: float,
        verbose: bool,
    ) -> _MRRefine:
        """Locate the multiple root(s) near *seg*'s stop point and the
        restart point just past them.

        Returns an :class:`_MRRefine`: the materialized *events*
        (θ-ascending; empty when every triggered "MR" was a false
        positive — the caller then holds the segment back as pending),
        the densely sampled *dense* rows bridging consecutive events of
        one bracket, and the restart point ``(theta_temp, roots_temp)``.

        For each event, ``roots_snapped`` is snapped to each cluster's
        mean (see :func:`snap_clusters_to_mean`) so the exact-degeneracy
        correction flows into the segment's closing row, the boundary_perm
        matching, and the MultipleRootInfo record alike.
        ``cluster_stds`` carries the per-cluster spread of the raw roots.
        """
        # MR-adjacent regular endpoint — the anchor for prediction-based
        # matching.  Comes straight from integrate_segment (no dependence on
        # the segment's last row, which may be a prepended boundary MR row).
        ref = seg.mr_ref
        ref2 = seg.mr_ref2
        roots_ref = ref.roots
        theta_ref = ref.theta
        if verbose:
            print("roots_ref = ", roots_ref)
        if seg.stop_reason == StopReason.multiple_root_in_interval:
            # Restart from the interval's right endpoint (MR-adjacent
            # regular row, past every event inside the bracket).  No
            # re-matching needed — ref.roots is already track-ordered.
            theta_temp = ref.theta
            roots_temp = ref.roots
            events = self._events_from_triggers(seg, ref, ref2, verbose)
            dense = (self._dense_rows_between(events)
                     if len(events) > 1 else [])
            return _MRRefine(events, dense, theta_temp, roots_temp)

        # multiple_root_encountered — refine with the iterative solver.
        beta1_approx = exp(self.mu1 + 1j * seg.mr_approx_theta)
        df2 = [self.poly.eval_partials((self.E_ref, beta1_approx, beta2))[2]
               for beta2 in roots_ref]
        min_idx = np.argmin(np.abs(df2))
        beta1, _ = solve_multiple_roots_iterative(
            self.poly,
            self.E_ref,
            beta1_approx,
            roots_ref[min_idx],
        )

        # Fix the 2π gauge to seg.mr_approx_theta.
        theta1_mr = np.angle(beta1 / exp(1j * seg.mr_approx_theta)) + seg.mr_approx_theta
        if abs(theta1_mr - seg.mr_approx_theta) > MR_GAUGE_TOL:
            warnings.warn(
                f"Segment ends at {seg.mr_approx_theta}, "
                f"but the iteration solver gives {theta1_mr}"
            )

        # Mirror the MR-adjacent point across the MR to land on the far
        # side, then single-end tangent-extrapolate the anchor to it
        # (consistent with the arclength integrator).
        #
        # The bare mirror 2θ_mr − θ_ref is only valid while θ_ref sits
        # BEFORE the MR; once refinement places θ_mr just behind the last
        # accepted row the mirror lands back BEFORE the MR and the restart
        # re-detects the same MR forever.  The lower clamp θ_mr +
        # mr_jump_eff (≥ the branch-point-safe distance, see run()) keeps
        # the restart strictly past the MR; the upper `min` against
        # θ_mr + mr_jump preserves the original behaviour of not jumping
        # farther than mr_jump when the mirror is distant.
        theta_temp = max(
            min(2 * theta1_mr - ref.theta, theta1_mr + mr_jump),
            theta1_mr + mr_jump_eff,
        )
        roots_temp = self._to_track_order(
            self._solve(theta_temp), theta_temp, theta_ref, roots_ref,
            ref_V=ref.V,
        )

        # θ₁_mr sits just past the last accepted row; only the MR-adjacent
        # regular point is nearby: single-end tangent extrapolation,
        # consistent with arclength.
        roots_mr = self._to_track_order(
            self._solve(theta1_mr), theta1_mr,
            theta_ref, roots_ref, ref_V=ref.V,
        )
        cluster = detect_cluster(roots_mr, cluster_tol=self._cluster_tol)
        cluster_stds: list[float] = []
        if cluster:
            roots_mr, cluster_stds = snap_clusters_to_mean(roots_mr, cluster)
            events = [_MREvent(theta1_mr, roots_mr, roots_mr, cluster,
                               list(cluster_stds))]
        else:
            events = []
        return _MRRefine(events, [], theta_temp, roots_temp)

    def _events_from_triggers(
        self,
        seg: SegmentResult,
        ref: _MREndpoint,
        ref2: _MREndpoint,
        verbose: bool,
    ) -> list[_MREvent]:
        """Materialize the MR events of one interval-trigger stop.

        Each triggered pair gets its own Brent refinement of ITS distance
        derivative inside its bracket, then the candidate θ₁s are grouped
        into events by ZERO COINCIDENCE, not by a θ tolerance: a candidate
        belongs to an existing event when its triggered pair is already
        within the cluster tolerance at that event's θ₁ — the same
        physical degeneracy reached through different track columns.  This
        merges simultaneous degeneracies of different zero pairs into ONE
        event with several clusters (detect_cluster enumerates all of them
        globally at that θ₁) while genuinely distinct nearby MRs stay
        separate events.  Candidates whose θ₁ shows no cluster at all are
        false positives and are dropped.

        All events of one stop share the trigger bracket, hence one
        two-endpoint Hermite anchor (ref2, ref) for the track-frame
        mapping.
        """
        # Per-pair Brent over the shared bracket: several pairs may flip
        # inside the same step (simultaneous degeneracies) — each is the
        # root of its own continuous per-pair g(θ).
        candidates: list[tuple[float, tuple[int, int]]] = []
        for rec in seg.mr_triggers:
            theta_c = solve_multiple_roots_in_interval(
                self.poly,
                self.E_ref,
                self.mu1,
                rec.theta_lo,
                rec.theta_hi,
                ref.roots,
                min_pair=rec.pair,
            )
            candidates.append((theta_c, rec.pair))

        events: list[_MREvent] = []
        for theta_c, pair in sorted(candidates, key=lambda c: c[0]):
            # Zero-coincidence test against earlier events: evaluate THIS
            # candidate's pair at the EVENT's θ₁, on the unsnapped roots.
            # The implicit θ resolution of the test, ~(CLUSTER_TOL/C)² for
            # a √-type approach with constant C, sits between Brent noise
            # and any physically distinct MR separation — no extra
            # tolerance constant needed.
            if any(_chordal_dist(ev.roots[pair[0]], ev.roots[pair[1]])
                   < self._cluster_tol
                   for ev in events):
                if verbose:
                    print(f"Trigger pair {pair}: coincides with an earlier "
                          f"event — duplicate θ₁ dropped")
                continue
            roots_ev = self._to_track_order(
                self._solve(theta_c), theta_c,
                ref2.theta, ref2.roots,
                ref_theta2=ref.theta, ref_roots2=ref.roots,
                ref_V=ref2.V, ref_V2=ref.V,
            )
            cluster = detect_cluster(roots_ev, cluster_tol=self._cluster_tol)
            if not cluster:
                if verbose:
                    print(f"Trigger pair {pair} at θ₁={theta_c!r}: "
                          "no cluster — false positive, dropped")
                continue
            roots_snapped, stds = snap_clusters_to_mean(roots_ev, cluster)
            events.append(_MREvent(theta_c, roots_ev, roots_snapped,
                                   cluster, list(stds)))
        return events

    def _dense_rows_between(
        self,
        events: list[_MREvent],
    ) -> list[tuple[np.ndarray, np.ndarray]]:
        """Regular rows sampled densely between consecutive bracket events.

        Two events inside one bracket can sit far below the θ resolution
        an adaptive restart could survive, so instead of restarting the
        integrator between them we sample roots directly on a uniform
        grid (spacing bounded by ``MR_DENSE_MAX_STEP``, at least
        ``MR_DENSE_MIN_SAMPLES`` interior rows) and flank the grid with
        the two MR closing rows — the same [left MR row | regular rows |
        right MR row] shape as any integrated segment.

        The first interior row is matched onto the left event's snapped
        (degenerate) closing row: like the θ=0 boundary-MR restart, the
        cluster columns' prediction degrades to hold-fixed there; every
        later row anchors on a regular row.
        """
        dense: list[tuple[np.ndarray, np.ndarray]] = []
        for ev_a, ev_b in zip(events[:-1], events[1:]):
            delta = ev_b.theta1 - ev_a.theta1
            n_interior = max(
                MR_DENSE_MIN_SAMPLES,
                int(np.ceil(delta / MR_DENSE_MAX_STEP)) - 1,
            )
            grid = np.linspace(ev_a.theta1, ev_b.theta1,
                               n_interior + 2)[1:-1]
            rows_theta = [ev_a.theta1]
            rows_roots = [ev_a.roots_snapped]
            anchor_theta, anchor_roots = ev_a.theta1, ev_a.roots_snapped
            for theta_s in grid:
                roots_s = self._to_track_order(
                    self._solve(theta_s), theta_s,
                    anchor_theta, anchor_roots,
                )
                rows_theta.append(theta_s)
                rows_roots.append(roots_s)
                anchor_theta, anchor_roots = theta_s, roots_s
            rows_theta.append(ev_b.theta1)
            rows_roots.append(ev_b.roots_snapped)
            dense.append((np.array(rows_theta), np.vstack(rows_roots)))
        return dense
