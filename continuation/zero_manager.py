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
    StepControl,
)

from .multiple_roots import (
    MultipleRootInfo,
    multiple_root_point_trigger,
    MultipleRootIntervalTrigger,
    detect_cluster,
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
    """
    theta1_arr: np.ndarray
    tracked_roots: np.ndarray
    abs_argsort: np.ndarray
    stop_reason: StopReason
    mr_approx_theta: float
    n_steps: int
    n_rejected: int
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
                return _finish(StopReason.multiple_root_encountered, theta1)

            # ---- Interval trigger: derivative sign flip ----
            # Closest-pair distance went from shrinking to growing — a
            # local minimum (multiple root) lies in [interval[0], interval[1]].
            # Catches "accidental" multiple roots where the tangent stays
            # well-conditioned and step size never collapses.
            ok, interval = interval_trigger(roots, result.V, theta1)
            if ok:
                return _finish(StopReason.multiple_root_in_interval, theta1, interval)

            theta1_new = result.theta1_new
            h = result.h_new
            roots_new = result.roots_new

            # Hungarian-match to maintain track continuity.
            matches = hungarian_match_indices(roots, roots_new)
            reordered = roots_new[matches]

            theta1_list.append(theta1_new)
            tracked_list.append(reordered)

            theta1 = theta1_new
            roots = reordered
            n_accepted += 1

        else:
            # max_iter exhausted without acceptable step → near singularity.
            n_rejected += 1
            return _finish(StopReason.multiple_root_encountered, theta1)

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
        self.left_boundary_roots = roots

        if verbose:
            print("Initial roots: ", roots)

        cluster = detect_cluster(roots, cluster_tol=cluster_tol)
        self.has_boundary_mr = len(cluster) > 0

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
                )
            )

            # Reinitialize the solver to avoid MR
            theta += mr_jump
            roots = self._to_track_order(self._solve(theta), roots)  # modulus-sorted
        else:
            left_mr = -1

        # Set True by the 'completed' or boundary-MR branches; the post-loop
        # fallback runs only when neither was hit.
        boundary_perm_set = False

        # ---- Integrate segments ----
        while theta < 2 * pi:
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
                # The left boundary point should be attached to the right boundary
                boundary_perm_inv = hungarian_match_indices(
                    new_seg_tracked_roots[-2, :],
                    self.left_boundary_roots,
                )

                # Convention (see ZeroManager docstring): roots_right[boundary_perm] == roots_left.
                self.boundary_perm = np.zeros(self.K, dtype=int)
                self.boundary_perm[boundary_perm_inv] = np.arange(self.K)
                boundary_perm_set = True

                new_seg_tracked_roots = np.vstack([
                    new_seg_tracked_roots[:-1, :],
                    self.left_boundary_roots[boundary_perm_inv]
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
                theta1_mr, roots_mr, cluster, theta_temp, roots_temp = self._refine_mr(
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
                        # Same convention as 'completed': left→right map.
                        self.boundary_perm = hungarian_match_indices(
                            self.left_boundary_roots, roots_mr
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
                                roots=roots_mr
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

        # Flush a pending false-positive segment if the loop ended before it
        # resolved (e.g. theta overshot 2π right after a false positive).
        if pending is not None:
            self._append_segment(
                pending.theta1_arr, pending.tracked_roots, pending.left_mr, -1
            )

        # Fall back: if neither 'completed' nor a boundary MR set boundary_perm
        # (the loop exited via theta >= 2π), match right to left by hand.
        if not boundary_perm_set:
            self.boundary_perm = hungarian_match_indices(
                self.left_boundary_roots, self.segments[-1].tracked_roots[-1]
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
        self, roots: np.ndarray, ref_left: np.ndarray, ref_right: np.ndarray = None
    ) -> np.ndarray:
        """Convert *roots* to track order.

        Hungarian-matches to the reference row *ref* (track-ordered).
        """
        perm_left = hungarian_match_indices(ref_left, roots)
        if ref_right is not None:
            perm_right = hungarian_match_indices(ref_right, roots)
            if np.any(perm_left != perm_right):
                raise ValueError(f"Roots do not match reference boundary. \nLeft: {perm_left}. Right: {perm_right}")
        return roots[perm_left]

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
    ) -> tuple[float, np.ndarray, list, float, np.ndarray]:
        """Locate the multiple root near *seg*'s stop point and the restart
        point just past it.

        Returns ``(theta1_mr, roots_mr, cluster, theta_temp, roots_temp)``.
        *cluster* is empty when the triggered "MR" was a false positive (no
        roots actually touch at *theta1_mr*); the caller then holds the
        segment back as pending instead of closing it.
        """
        roots_ref = new_seg_tracked_roots[-1, :]
        if verbose:
            print("roots_ref = ", roots_ref)
        if seg.stop_reason == StopReason.multiple_root_in_interval:
            # Restart from the segment's right edge.
            theta_temp = new_seg_theta1[-1]
            roots_temp = new_seg_tracked_roots[-1, :]
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

            # Mirror the previous sample across the MR to land on the far side.
            theta_temp = min(2 * theta1_mr - new_seg_theta1[-2], theta1_mr + mr_jump)
            roots_temp = self._to_track_order(self._solve(theta_temp), roots_ref)

        roots_mr = self._to_track_order(self._solve(theta1_mr), roots_ref)
        cluster = detect_cluster(roots_mr, cluster_tol=self._cluster_tol)
        return theta1_mr, roots_mr, cluster, theta_temp, roots_temp
