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
from math import pi
from cmath import exp
from dataclasses import dataclass

from typing import Optional

from gbz_types import (
    CharPoly,
    hungarian_match_indices,
    sort_by_root_abs,
    to_sphere_r3,
    cost_from_sphere_r3,
)
from continuation.arclength import (
    integrate_segment,
    detect_cluster,
    refine_multiple_root_theta,
    SegmentResult,
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
      - ``mr_thetas`` : list[float] — refined MR θ₁ values.
      - ``mr_cluster_masks`` : list[bool[K]] — which tracks are in the cluster.
      - ``mr_roots`` : list[complex[K]] — track-ordered roots at each MR.
      - ``segments`` : list[SegmentData] — curve segments, each with
        ``left_mr`` / ``right_mr`` indices into the three lists above.

    ``.insert(theta1)`` appends a point to the appropriate segment.
    """

    def __init__(self, poly: CharPoly, E_ref: complex, mu1: float):
        self.poly = poly
        self.E_ref = E_ref
        self.mu1 = mu1
        self.M = poly.M
        self.K = poly.M + poly.N

        self.mr_thetas: list[float] = []
        self.mr_cluster_masks: list[np.ndarray] = []
        self.mr_roots: list[np.ndarray] = []
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
        cluster_tol: float = 1e-6,
        atol: float = 1e-12,
        rtol: float = 1e-3,
        mr_jump: float = _MR_JUMP,
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

        # ---- Step 1: initialise at θ₁ = 0 ----
        theta = 0.0
        roots = self._solve(theta)  # modulus-sorted

        # Check for multiple root at the starting point.
        initial_cluster = detect_cluster(roots, cluster_tol=cluster_tol)
        has_initial_mr = initial_cluster is not None
        if has_initial_mr:
            theta_mr = _normalize_theta(
                refine_multiple_root_theta(
                    self.poly, self.E_ref, self.mu1, theta,
                    cluster_tol=cluster_tol,
                ))
            self._record_mr(theta_mr, roots)
            theta = _normalize_theta(theta_mr + mr_jump)
            roots = self._solve(theta)

        prev_mr = 0 if has_initial_mr else -1

        # ---- Step 2: integrate segments ----
        while theta < 2 * pi - 1e-12:
            seg = integrate_segment(
                self.poly, self.E_ref, self.mu1,
                theta, roots, 2 * pi,
                h0=h0, max_step=max_step,
                min_step=min_step, min_dtheta=min_dtheta,
                cluster_tol=cluster_tol,
                atol=atol, rtol=rtol,
            )

            if seg.stop_reason == 'completed':
                # Build final segment.  If there is an initial MR, the
                # segment wraps around to it (θ=2π ≡ θ=0).
                right_mr = 0 if has_initial_mr else -1
                append = None
                if has_initial_mr:
                    ref = (seg.tracked_roots[-1] if len(seg.tracked_roots) > 0
                           else roots)
                    append = self._to_track_order(self.mr_roots[0], ref)
                seg_data = self._make_segment(
                    seg, theta, roots, prev_mr, right_mr, append_mr=append,
                )
                self.segments.append(seg_data)
                break

            # ---- Multiple root detected ----
            # Refine the MR θ₁.
            theta_mr = _normalize_theta(
                refine_multiple_root_theta(
                    self.poly, self.E_ref, self.mu1,
                    seg.mr_approx_theta,
                    cluster_tol=cluster_tol,
                ))

            # Solve roots at the refined MR.
            roots_mr = self._solve(theta_mr)

            # Record the MR.
            self._record_mr(theta_mr, roots_mr)
            new_mr = len(self.mr_thetas) - 1

            # Build segment: start point + integration data + MR end point.
            ref = (seg.tracked_roots[-1] if len(seg.tracked_roots) > 0
                   else roots)
            tracked_mr = self._to_track_order(roots_mr, ref)
            seg_data = self._make_segment(
                seg, theta, roots, prev_mr, new_mr, append_mr=tracked_mr,
            )
            self.segments.append(seg_data)

            # Jump past the MR for the next segment.
            prev_mr = new_mr
            theta = _normalize_theta(theta_mr + mr_jump)
            roots = self._solve(theta)

    def insert(self, theta1: float) -> tuple[int, int, np.ndarray]:
        """Solve roots at *theta1*, insert into the appropriate segment.

        Returns
        -------
        (seg_idx, local_idx, track_ordered_roots)
        """
        theta_norm = _normalize_theta(theta1)
        roots_new = self._solve(theta_norm)
        seg_idx = self._find_segment(theta_norm)
        seg = self.segments[seg_idx]
        tracked = self._to_track_order(roots_new, seg.tracked_roots[-1])

        local_idx = int(np.searchsorted(seg.theta1_arr, theta_norm))
        if local_idx >= len(seg.theta1_arr):
            local_idx = len(seg.theta1_arr)

        seg.theta1_arr = np.insert(seg.theta1_arr, local_idx, theta_norm)
        seg.tracked_roots = np.insert(
            seg.tracked_roots, local_idx,
            tracked[np.newaxis, :], axis=0,
        )
        seg.abs_argsort = np.argsort(np.abs(seg.tracked_roots), axis=1)

        return seg_idx, local_idx, tracked

    @property
    def n_multiple_roots(self) -> int:
        return len(self.mr_thetas)

    @property
    def n_segments(self) -> int:
        return len(self.segments)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _solve(self, theta1: float) -> np.ndarray:
        """Solve β₂ roots at *theta1*, return modulus-sorted."""
        beta1 = exp(self.mu1 + 1j * _normalize_theta(theta1))
        roots = self.poly.solve_roots_1d((0, 1), (self.E_ref, beta1), (2,))
        return sort_by_root_abs(np.asarray(roots, dtype=complex))

    def _to_track_order(
        self, roots: np.ndarray, ref: np.ndarray,
    ) -> np.ndarray:
        """Convert modulus-sorted *roots* to track order.

        Hungarian-matches to the reference row *ref* (track-ordered).
        """
        matches = hungarian_match_indices(ref, roots)
        tracked = np.zeros(self.K, dtype=complex)
        for from_idx, to_idx in matches:
            tracked[from_idx] = roots[to_idx]
        return tracked

    def _record_mr(
        self, theta_mr: float, roots_mod: np.ndarray,
    ) -> None:
        """Record a multiple root: θ₁, cluster mask (mod-sorted), roots."""
        cluster = detect_cluster(roots_mod, cluster_tol=self._cluster_tol)
        mask = np.zeros(self.K, dtype=bool)
        if cluster is not None:
            for ci in cluster:
                mask[int(ci)] = True
        self.mr_thetas.append(theta_mr)
        self.mr_cluster_masks.append(mask)
        # Store modulus-sorted roots; caller can convert to track order.
        self.mr_roots.append(roots_mod.copy())

    def _make_segment(
        self,
        seg: SegmentResult,
        theta_start: float,
        roots_start: np.ndarray,
        left_mr: int,
        right_mr: int,
        append_mr: Optional[np.ndarray],
    ) -> SegmentData:
        """Build a SegmentData from a SegmentResult + start point + optional MR end.

        *seg* does not include the start point; we prepend it here.
        """
        # Prepend the start point (track order defined by roots_start).
        theta1_arr = np.concatenate([
            [theta_start], seg.theta1_arr,
        ])
        tracked_roots = np.vstack([
            roots_start[np.newaxis, :], seg.tracked_roots,
        ]) if len(seg.tracked_roots) > 0 else roots_start[np.newaxis, :]

        # Optionally append the MR end point.
        if append_mr is not None:
            mr_theta = self.mr_thetas[right_mr] if right_mr >= 0 else 0.0
            if mr_theta <= theta1_arr[-1]:
                mr_theta += 2 * pi
            theta1_arr = np.append(theta1_arr, mr_theta)
            tracked_roots = np.vstack([
                tracked_roots, append_mr[np.newaxis, :],
            ])

        abs_argsort = np.argsort(np.abs(tracked_roots), axis=1)

        return SegmentData(
            theta1_arr=theta1_arr,
            tracked_roots=tracked_roots,
            abs_argsort=abs_argsort,
            left_mr=left_mr,
            right_mr=right_mr,
        )

    def _find_segment(self, theta1: float) -> int:
        """Find which segment contains *theta1*."""
        theta_norm = _normalize_theta(theta1)
        for i, seg in enumerate(self.segments):
            t0 = _normalize_theta(seg.theta1_arr[0])
            t1 = _normalize_theta(seg.theta1_arr[-1])
            if t1 >= t0:
                if t0 <= theta_norm <= t1:
                    return i
            else:
                if theta_norm >= t0 or theta_norm <= t1:
                    return i
        best = 0
        best_dist = float("inf")
        for i, seg in enumerate(self.segments):
            mid = _normalize_theta(
                0.5 * (seg.theta1_arr[0] + seg.theta1_arr[-1])
            )
            d = abs(theta_norm - mid)
            d = min(d, 2 * pi - d)
            if d < best_dist:
                best_dist = d
                best = i
        return best
