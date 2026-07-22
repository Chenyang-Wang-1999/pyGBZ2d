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
    SAFETY,
    MIN_FACTOR,
    MAX_FACTOR,
    ERROR_EXPONENT,
    arclength_step
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





def integrate_segment(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    theta_start: float,
    roots_start: np.ndarray,
    theta_end: float,
    *,
    h0: float = 0.1,
    max_step: float = 0.5,
    min_step: float = 1e-12,
    min_dtheta: float = 1e-10,
    min_dist_threshold: float = 0.1,
    atol: float = 1e-12,
    rtol: float = 1e-3,
    safety: float = SAFETY,
    min_factor: float = MIN_FACTOR,
    max_factor: float = MAX_FACTOR,
    error_exponent: float = ERROR_EXPONENT,
    max_step_iter: int = 20,
) -> SegmentResult:
    """Integrate β₂ roots from *theta_start* toward *theta_end*.

    Stops when *theta_end* is reached (``stop_reason=StopReason.completed``) or when
    a multiple root is detected.  Two trigger mechanisms run after each
    accepted step:

    * ``'multiple_root_encountered'`` — step-size collapse (point trigger).
    * ``'multiple_root_in_interval'`` — closest-pair distance derivative
      sign flip (interval trigger).

    *roots_start* must be track-ordered (column j = physical root j).
    The returned ``tracked_roots`` does **not** include the start point —
    the caller already has it.  Each subsequent row is Hungarian-matched
    to maintain track continuity within the segment.
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

    while theta1 < theta_end:
        result = arclength_step(
            poly, E_ref, mu1, theta1, roots, h,
            direction=1.0,
            max_step=max_step,
            min_step=min_step,
            atol=atol,
            rtol=rtol,
            safety=safety,
            min_factor=min_factor,
            max_factor=max_factor,
            error_exponent=error_exponent,
            max_iter=max_step_iter,
        )

        if result.accepted:
            dtheta = abs(result.theta1_new - theta1)

            # ---- Point trigger: step-size collapse ----
            # We are AT a multiple root (tangent diverged → step collapsed).
            if multiple_root_point_trigger(dtheta, min_dtheta=min_dtheta):
                theta1_arr = np.array(theta1_list)
                tracked_roots_arr = np.array(tracked_list)
                abs_argsort = np.argsort(np.abs(tracked_roots_arr), axis=1) if len(tracked_list) > 0 else np.empty((0, n_roots), dtype=int)
                return SegmentResult(
                    theta1_arr=theta1_arr,
                    tracked_roots=tracked_roots_arr,
                    abs_argsort=abs_argsort,
                    stop_reason=StopReason.multiple_root_encountered,
                    mr_approx_theta=theta1,
                    n_steps=n_accepted,
                    n_rejected=n_rejected,
                )

            # ---- Interval trigger: derivative sign flip ----
            # Closest-pair distance went from shrinking to growing — a
            # local minimum (multiple root) lies in [interval[0], interval[1]].
            # Catches "accidental" multiple roots where the tangent stays
            # well-conditioned and step size never collapses.
            ok, interval = interval_trigger(roots, result.V, theta1)
            if ok:
                theta1_arr = np.array(theta1_list)
                tracked_roots_arr = np.array(tracked_list)
                abs_argsort = np.argsort(np.abs(tracked_roots_arr), axis=1) if len(tracked_list) > 0 else np.empty((0, n_roots), dtype=int)
                return SegmentResult(
                    theta1_arr=theta1_arr,
                    tracked_roots=tracked_roots_arr,
                    abs_argsort=abs_argsort,
                    stop_reason=StopReason.multiple_root_in_interval,
                    mr_approx_theta=theta1,
                    mr_interval_start=interval[0],
                    mr_interval_end=interval[1],
                    n_steps=n_accepted,
                    n_rejected=n_rejected,
                )

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
            theta1_arr = np.array(theta1_list)
            tracked_roots_arr = np.array(tracked_list)
            abs_argsort = np.argsort(np.abs(tracked_roots_arr), axis=1) if len(tracked_list) > 0 else np.empty((0, n_roots), dtype=int)
            return SegmentResult(
                theta1_arr=theta1_arr,
                tracked_roots=tracked_roots_arr,
                abs_argsort=abs_argsort,
                stop_reason=StopReason.multiple_root_encountered,
                mr_approx_theta=theta1,
                n_steps=n_accepted,
                n_rejected=n_rejected,
            )

    # Completed without hitting a multiple root.
    theta1_arr = np.array(theta1_list)
    tracked_roots_arr = np.array(tracked_list)
    abs_argsort = np.argsort(np.abs(tracked_roots_arr), axis=1) if len(tracked_list) > 0 else np.empty((0, n_roots), dtype=int)
    return SegmentResult(
        theta1_arr=theta1_arr,
        tracked_roots=tracked_roots_arr,
        abs_argsort=abs_argsort,
        stop_reason=StopReason.completed,
        mr_approx_theta=float('nan'),
        n_steps=n_accepted,
        n_rejected=n_rejected,
    )



def _normalize_theta(theta: float) -> float:
    """Normalize to [0, 2π).  Values within 1e-6 of 2π snap to 0."""
    t = float(theta % (2 * pi))
    if t > 2 * pi - 1e-6:
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


# ---------------------------------------------------------------------------
# ZeroManager
# ---------------------------------------------------------------------------

# Small θ₁ step used to jump past a multiple root after refinement.
_MR_JUMP = 1e-6


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
        max_step: float = 0.5,
        min_step: float = 1e-12,
        min_dtheta: float = 1e-10,
        cluster_tol: float = 1e-4,
        atol: float = 1e-12,
        rtol: float = 1e-3,
        mr_jump: float = _MR_JUMP,
        verbose: bool = False,
    ) -> None:
        """Execute the full pipeline.

        Parameters
        ----------
        h0, max_step, min_step, min_dtheta :
            Passed to ``integrate_segment``.
        cluster_tol :
            Chordal-distance threshold for multiple-root detection.
        atol, rtol :
            Tolerances for error estimation.
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

        if self.has_boundary_mr:
            prev_mr = 0
            # Add an multiple root at theta1 = 0
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
            prev_mr = -1


        # ---- Integrate segments ----
        while theta < 2 * pi:
            seg = integrate_segment(
                self.poly, self.E_ref, self.mu1,
                theta, roots, 2 * pi,
                h0=h0, max_step=max_step,
                min_step=min_step, min_dtheta=min_dtheta,
                atol=atol, rtol=rtol,
            )
            if verbose:
                print("====== Segment solved =======")
                print("theta1 range: ", (seg.theta1_arr[0], seg.theta1_arr[-1]))
                print("stop_reason: ", seg.stop_reason.value)

            # Common processes: 
            #  add left boundary points to this segment
            if prev_mr >= 0:
                # Has multiple root on its left endpoint
                new_seg_theta1 = np.concatenate([
                    [self.multiple_roots[prev_mr].theta1],
                    seg.theta1_arr
                ])
                new_seg_tracked_roots = np.vstack([
                    self.multiple_roots[prev_mr].roots,
                    seg.tracked_roots
                ])
            elif prev_mr == -2:
                # False positive trigger
                prev_seg = self.segments.pop(-1)

                new_seg_theta1 = np.concatenate([
                    prev_seg.theta1_arr,
                    seg.theta1_arr
                ])
                new_seg_tracked_roots = np.vstack([
                    prev_seg.tracked_roots,
                    seg.tracked_roots
                ])
                prev_mr = prev_seg.left_mr
            else:
                new_seg_theta1 = seg.theta1_arr
                new_seg_tracked_roots = seg.tracked_roots


            if seg.stop_reason == StopReason.completed:
                # The left boundary point should be attached to the right boundary

                boundary_perm_inv = hungarian_match_indices(
                    new_seg_tracked_roots[-2, :],
                    self.left_boundary_roots,
                )

                # Here, we adopt the convension that left_boundary_roots = right_boundary_roots[boundary_perm]
                self.boundary_perm = np.zeros(self.K, dtype=int)
                self.boundary_perm[boundary_perm_inv] = np.arange(self.K)

                new_seg_tracked_roots = np.vstack([
                    new_seg_tracked_roots[:-1, :],
                    self.left_boundary_roots[boundary_perm_inv]
                ])
                self.segments.append(
                    SegmentData(
                        theta1_arr=np.concatenate([
                            new_seg_theta1[:-1],
                            [2 * pi]
                        ]),
                        tracked_roots=new_seg_tracked_roots,
                        abs_argsort=np.argsort(np.abs(new_seg_tracked_roots), axis=1),
                        left_mr=prev_mr,
                        right_mr=0 if self.has_boundary_mr else -1
                    )
                )

                break

            elif (seg.stop_reason == StopReason.multiple_root_in_interval
                  or seg.stop_reason == StopReason.multiple_root_encountered
            ):
                roots_ref = new_seg_tracked_roots[-1, :]
                if verbose:
                    print("roots_ref = ", roots_ref)
                if seg.stop_reason == StopReason.multiple_root_in_interval:
                    # Save right boundary point
                    theta_temp = new_seg_theta1[-1]
                    roots_temp = new_seg_tracked_roots[-1, :]
                    theta1_mr = solve_multiple_roots_in_interval(
                        self.poly,
                        self.E_ref,
                        self.mu1,
                        seg.mr_interval_start,
                        seg.mr_interval_end,
                        roots_ref,
                        cluster_tol=self._cluster_tol,
                    )
                else:
                    # Solve multiple roots using iterative solver
                    beta1_approx = exp(self.mu1 + 1j * seg.mr_approx_theta)
                    df2 = [self.poly.eval_partials((self.E_ref, beta1_approx, beta2))[2] for beta2 in roots_ref]
                    min_idx = np.argmin(np.abs(df2))
                    beta1, _ = solve_multiple_roots_iterative(
                        self.poly,
                        self.E_ref,
                        beta1_approx,
                        roots_ref[min_idx]
                    )

                    # Fix 2 pi gauge according to seg.mr_approx_theta
                    theta1_mr = np.angle(beta1 / exp(1j * seg.mr_approx_theta)) + seg.mr_approx_theta
                    if abs(theta1_mr - seg.mr_approx_theta) > 0.1:
                        warnings.warn(f"Segment ends at {seg.mr_approx_theta}, but the iteration solver gives {theta1_mr}")

                    theta_temp = min(2 * theta1_mr - new_seg_theta1[-2], theta1_mr + mr_jump)
                    roots_temp = self._to_track_order(self._solve(theta_temp), roots_ref)
 
                roots_mr = self._to_track_order(self._solve(theta1_mr), roots_ref)
                cluster = detect_cluster(roots_mr, cluster_tol=self._cluster_tol)
                if len(cluster) > 0:
                    if verbose:
                        print("Find multiple roots. Cluster = ", cluster)
                    # Process segment data
                    seg_used = new_seg_theta1 < theta1_mr
                    new_seg_theta1 = np.concatenate([
                        new_seg_theta1[seg_used],
                        [theta1_mr]
                    ])
                    new_seg_tracked_roots = np.vstack([
                        new_seg_tracked_roots[seg_used, :],
                        roots_mr
                    ])

                    if abs(theta1_mr - 2 * pi) < 1e-6:
                        if verbose:
                            print("Boundary multiple root detected. right_mr = 0")
                        # Right mr = initial mr
                        self.segments.append(
                            SegmentData(
                                theta1_arr=new_seg_theta1,
                                tracked_roots=new_seg_tracked_roots,
                                abs_argsort=np.argsort(np.abs(new_seg_tracked_roots), axis=1),
                                left_mr=prev_mr,
                                right_mr=0
                            )
                        )
                        self.boundary_perm = hungarian_match_indices(roots_mr, self.left_boundary_roots)
                        break

                    else:
                        self.segments.append(
                            SegmentData(
                                theta1_arr=new_seg_theta1,
                                tracked_roots=new_seg_tracked_roots,
                                abs_argsort=np.argsort(np.abs(new_seg_tracked_roots), axis=1),
                                left_mr=prev_mr,
                                right_mr=len(self.multiple_roots)
                            )
                        )

                        # Add multiple root to the segment
                        prev_mr = len(self.multiple_roots)

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
                    # restart from the right edge
                    self.segments.append(
                        SegmentData(
                            theta1_arr=new_seg_theta1[:-1],
                            tracked_roots=new_seg_tracked_roots[:-1, :],
                            abs_argsort=np.argsort(np.abs(new_seg_tracked_roots[:-1, :]), axis=1),
                            left_mr=prev_mr,
                            right_mr=-2
                        )
                    )
                    prev_mr = -2

                theta = theta_temp
                roots = roots_temp

            else:
                raise RuntimeError(
                    f"Unknown stop_reason: {seg.stop_reason.value}"
                )

        # Fall back: if StopReason.completed is not encoutered, match right to left by hand
        self.boundary_perm = hungarian_match_indices(self.segments[-1].tracked_roots[-1], self.left_boundary_roots)
        if self.segments[-1].right_mr == -2:
            # False positive trigger of multiple roots. No need to do anything except for attaching it to the left boundary
            seg = self.segments[-1]
            self.segments[-1] = SegmentData(
                theta1_arr=seg.theta1_arr,
                tracked_roots=seg.tracked_roots,
                abs_argsort=seg.abs_argsort,
                left_mr=seg.left_mr,
                right_mr=-1,
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
