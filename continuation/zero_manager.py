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
from typing import NamedTuple
import warnings


from gbz_types import (
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

from .multiple_roots import (
    MultipleRootInfo,
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
    * ``'multiple_root_in_interval'`` — closest-pair distance derivative
      flipped sign, so a local minimum (multiple root) lies in
      [*mr_interval_start*, *mr_interval_end*].  ZeroManager should bisect
      the interval to locate it, then refine.

    For both MR stop reasons, *mr_ref* is the MR-adjacent regular endpoint
    (the segment's last accepted row) — used as the single-endpoint
    prediction anchor for ``multiple_root_encountered``.  For
    ``multiple_root_in_interval``, *mr_ref2* is the interval's other
    endpoint so the caller can build a two-endpoint cubic Hermite anchor
    straddling the bracketed MR.
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
    # Valid only for 'multiple_root_in_interval':
    mr_interval_start: float = float('nan')
    mr_interval_end: float = float('nan')


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
    h0: float = 0.1,
    ctrl: StepControl = StepControl(),
    min_dtheta: float = 1e-10,
    min_dist_threshold: float = 0.1,
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
                interval: tuple[float, float] = (float('nan'), float('nan')),
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
            mr_interval_start=interval[0],
            mr_interval_end=interval[1],
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
            # Closest-pair distance went from shrinking to growing — a
            # local minimum (multiple root) lies in [interval[0], interval[1]].
            # Catches "accidental" multiple roots where the tangent stays
            # well-conditioned and step size never collapses.
            ok, interval = interval_trigger(roots, result.V, theta1)
            if ok:
                # Two-endpoint anchor: the trigger's prev state (interval
                # start) and the current accepted row (interval end) bracket
                # the MR — both endpoints share the segment's track frame,
                # so a cubic Hermite between them is well-posed.
                ref2 = _MREndpoint(interval_trigger._prev_theta,
                                   interval_trigger._prev_roots.copy(),
                                   interval_trigger._prev_V.copy())
                ref = _MREndpoint(theta1, roots.copy(), result.V.copy())
                return _finish(StopReason.multiple_root_in_interval, theta1,
                               interval, mr_ref=ref, mr_ref2=ref2)

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
    """Normalize to [0, 2π).  Values within _BOUNDARY_THETA_TOL of 2π snap to 0."""
    t = float(theta % (2 * pi))
    if t > 2 * pi - _BOUNDARY_THETA_TOL:
        t = 0.0
    return t


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
        MR index at the left boundary (-1 = start of circle).
    right_mr : int
        MR index at the right boundary (-1 = end of circle).
    """
    theta1_arr: np.ndarray
    tracked_roots: np.ndarray
    abs_argsort: np.ndarray
    left_mr: int
    right_mr: int

    def __post_init__(self):
        n = len(self.theta1_arr)
        if self.tracked_roots.shape[0] != n:
            raise ValueError(
                f"theta1_arr length {n} != tracked_roots rows "
                f"{self.tracked_roots.shape[0]}"
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


# ---------------------------------------------------------------------------
# ZeroManager
# ---------------------------------------------------------------------------

# Small θ₁ step used to jump past a multiple root after refinement.
_MR_JUMP = 1e-6
# θ₁ within this of 2π (≡ 0) is treated as the θ₁ = 0 boundary.
_BOUNDARY_THETA_TOL = 1e-6
# Warn if the iterative MR solver's θ₁ drifts more than this from the trigger.
_MR_GAUGE_TOL = 0.1
# Hard cap on segment count — guards against a runaway MR-refine loop
# that never reaches θ₁ = 2π.  Generous: K roots admit at most O(K) MRs.
_MAX_SEGMENTS = 10000


class ZeroManager:
    """β₂-root topology over θ₁ ∈ [0, 2π) at fixed (E, μ₁).

    After ``.run()``:
      - ``multiple_roots`` : list[MultipleRootInfo] — each MR's θ₁, cluster
        indices, and modulus-sorted β₂ roots.
      - ``segments`` : list[SegmentData] — curve segments, each with
        ``left_mr`` / ``right_mr`` indices into ``multiple_roots``
        (-1 = no MR at that boundary).
      - ``boundary_perm``: np.ndarray (int, K) permutation from right boundary
        to left boundary: roots_right[boundary_perm] == roots_left.
      - ``has_boundary_mr``: bool — whether there is a multiple root at θ₁=0.
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
        self._cluster_tol: float = 1e-6

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(
        self,
        *,
        h0: float = 0.1,
        ctrl: StepControl = StepControl(),
        min_dtheta: float = 1e-10,
        cluster_tol: float = 1e-4,
        mr_jump: float = _MR_JUMP,
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
        """
        self._cluster_tol = cluster_tol

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
                    cluster_stds=cluster_stds,
                )
            )

            # Reinitialize the solver to avoid MR
            theta += mr_jump
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

        # ---- Integrate segments ----
        # `for` with a hard cap instead of `while theta < 2π`: the integrator
        # already advances θ to 2π (completed) or an MR; the cap only catches
        # a runaway refine loop that never converges to 2π.
        for _ in range(_MAX_SEGMENTS):
            seg = integrate_segment(
                self.poly, self.E_ref, self.mu1,
                theta, roots, 2 * pi,
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
                    np.concatenate([new_seg_theta1[:-1], [2 * pi]]),
                    new_seg_tracked_roots,
                    left_mr,
                    0 if self.has_boundary_mr else -1,
                )
                break

            elif (seg.stop_reason == StopReason.multiple_root_in_interval
                  or seg.stop_reason == StopReason.multiple_root_encountered
            ):
                theta1_mr, roots_mr, cluster, cluster_stds, theta_temp, roots_temp = self._refine_mr(
                    seg, new_seg_theta1, new_seg_tracked_roots, mr_jump, verbose
                )
                if cluster:
                    if verbose:
                        print("Find multiple roots. Cluster = ", cluster)
                    # Trim the segment up to the MR and close it there.
                    seg_used = new_seg_theta1 < theta1_mr
                    new_seg_theta1 = np.concatenate([
                        new_seg_theta1[seg_used],
                        [theta1_mr]
                    ])
                    new_seg_tracked_roots = np.vstack([
                        new_seg_tracked_roots[seg_used, :],
                        roots_mr
                    ])

                    if abs(theta1_mr - 2 * pi) < _BOUNDARY_THETA_TOL:
                        if verbose:
                            print("Boundary multiple root detected. right_mr = 0")
                        self._append_segment(
                            new_seg_theta1, new_seg_tracked_roots, left_mr, 0
                        )
                        # Same convention as 'completed': predict at 2π and
                        # match onto left_boundary_roots.  The boundary MR at
                        # θ₁_mr ≈ 2π is degenerate (closing row is the snapped
                        # cluster), so the bracket endpoint carrying it has a
                        # zeroed tangent and the Hermite degrades toward lerp
                        # / hold-fixed — no regression versus a bare match, and
                        # the cluster↔cluster correspondence stays ambiguous.
                        predicted_2pi = self._predict_roots_at_2pi(
                            new_seg_theta1, new_seg_tracked_roots,
                        )
                        self.boundary_perm, _ = self._boundary_perm_from_right(
                            predicted_2pi
                        )
                        boundary_perm_set = True
                        break
                    else:
                        right_mr = len(self.multiple_roots)
                        self._append_segment(
                            new_seg_theta1, new_seg_tracked_roots, left_mr, right_mr
                        )
                        left_mr = right_mr
                        self.multiple_roots.append(
                            MultipleRootInfo(
                                theta1=theta1_mr,
                                cluster_indices=cluster,
                                roots=roots_mr,
                                cluster_stds=cluster_stds,
                            )
                        )
                else:
                    if verbose:
                        print(f"No multiple roots found with cluster_tol = {self._cluster_tol: .2e}. \n Roots:", roots_mr)
                    # False positive: hold this segment back as pending instead
                    # of appending it with a sentinel right_mr.
                    pending = _PendingSeg(
                        new_seg_theta1[:-1],
                        new_seg_tracked_roots[:-1, :],
                        left_mr,
                    )

                theta = theta_temp
                roots = roots_temp

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
                f"{_MAX_SEGMENTS} segments; stuck near θ₁ = {theta}"
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
        assert boundary_perm_set, (
            "ZeroManager loop exited without setting boundary_perm"
        )

    @property
    def n_multiple_roots(self) -> int:
        return len(self.multiple_roots)

    @property
    def n_segments(self) -> int:
        return len(self.segments)

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
        below = theta1_arr < 2 * pi
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
                2 * pi, theta_a, roots_a, Va, theta_b, roots_b, Vb,
            )

        # No bracketing pair: single-end extrapolate from the row nearest 2π.
        src = idx if idx >= 0 else 0
        theta_s = theta1_arr[src]
        roots_s = tracked_roots[src, :]
        Vs = compute_tangent(
            self.poly, self.E_ref, exp(self.mu1 + 1j * theta_s), roots_s,
        )[0]
        return predict_roots_hermite(2 * pi, theta_s, roots_s, Vs)

    # ------------------------------------------------------------------
    # run() helpers
    # ------------------------------------------------------------------

    def _append_segment(
        self,
        theta1_arr: np.ndarray,
        tracked_roots: np.ndarray,
        left_mr: int,
        right_mr: int,
    ) -> None:
        """Append a SegmentData with abs_argsort computed from tracked_roots."""
        self.segments.append(
            SegmentData(
                theta1_arr=theta1_arr,
                tracked_roots=tracked_roots,
                abs_argsort=_abs_argsort(tracked_roots),
                left_mr=left_mr,
                right_mr=right_mr,
            )
        )

    def _refine_mr(
        self,
        seg: SegmentResult,
        new_seg_theta1: np.ndarray,
        new_seg_tracked_roots: np.ndarray,
        mr_jump: float,
        verbose: bool,
    ) -> tuple[float, np.ndarray, list, list[float], float, np.ndarray]:
        """Locate the multiple root near *seg*'s stop point and the restart
        point just past it.

        Returns ``(theta1_mr, roots_mr, cluster, cluster_stds, theta_temp,
        roots_temp)``.  *cluster* is empty when the triggered "MR" was a
        false positive (no roots actually touch at *theta1_mr*); the caller
        then holds the segment back as pending instead of closing it.

        When *cluster* is non-empty, *roots_mr* is snapped to each cluster's
        mean (see :func:`snap_clusters_to_mean`) so the exact-degeneracy
        correction flows into the segment's closing row, the boundary_perm
        matching, and the MultipleRootInfo record alike.  *cluster_stds*
        carries the per-cluster spread of the raw roots.
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
            # Restart from the interval's right endpoint (MR-adjacent regular
            # row).  No re-matching needed — ref.roots is already track-ordered.
            theta_temp = ref.theta
            roots_temp = ref.roots
            theta1_mr = solve_multiple_roots_in_interval(
                self.poly,
                self.E_ref,
                self.mu1,
                seg.mr_interval_start,
                seg.mr_interval_end,
                roots_ref,
            )
        else:
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
            if abs(theta1_mr - seg.mr_approx_theta) > _MR_GAUGE_TOL:
                warnings.warn(
                    f"Segment ends at {seg.mr_approx_theta}, "
                    f"but the iteration solver gives {theta1_mr}"
                )

            # Mirror the MR-adjacent point across the MR to land on the far
            # side, then single-end tangent-extrapolate the anchor to it
            # (consistent with the arclength integrator).
            theta_temp = min(2 * theta1_mr - ref.theta, theta1_mr + mr_jump)
            roots_temp = self._to_track_order(
                self._solve(theta_temp), theta_temp, theta_ref, roots_ref,
                ref_V=ref.V,
            )

        # Locate the MR itself.
        #   multiple_root_in_interval → θ₁_mr is interior to the trigger
        #     bracket [ref2.theta, ref.theta] whose endpoints both come from
        #     integrate_segment and share the segment's track frame:
        #     two-endpoint cubic Hermite.
        #   multiple_root_encountered → θ₁_mr sits just past the last accepted
        #     row; only the MR-adjacent regular point is nearby: single-end
        #     tangent extrapolation, consistent with arclength.
        if seg.stop_reason == StopReason.multiple_root_in_interval:
            roots_mr = self._to_track_order(
                self._solve(theta1_mr), theta1_mr,
                ref2.theta, ref2.roots,
                ref_theta2=ref.theta, ref_roots2=ref.roots,
                ref_V=ref2.V, ref_V2=ref.V,
            )
        else:
            roots_mr = self._to_track_order(
                self._solve(theta1_mr), theta1_mr,
                theta_ref, roots_ref, ref_V=ref.V,
            )
        cluster = detect_cluster(roots_mr, cluster_tol=self._cluster_tol)
        cluster_stds: list[float] = []
        if cluster:
            roots_mr, cluster_stds = snap_clusters_to_mean(roots_mr, cluster)
        return theta1_mr, roots_mr, cluster, cluster_stds, theta_temp, roots_temp
