'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-16
Copyright © Department of Physics, Tsinghua University. All rights reserved

SGBZ continuum: detection + 1D LineSubset materialization.

With exact pairwise event rows available, a LineSubset is no longer read off
raw ``j_lo == j_hi`` row runs: events that change the boundary pair
(``EventGroup.changes_boundary``) are exact terminators.  An event that does
NOT change ``(j_lo, j_hi)`` (e.g. an M+1/M+2 exchange) does not terminate the
line and the run simply continues through it.
'''

from __future__ import annotations

import numpy as np
from typing import Optional

from ..core import CharPoly, LineSubset
from ..continuation import ZeroManager

from ..core import live_defaults
from .mu2mid import ensure_mu2mid
from .mu2mid import Mu2MidZM
from ..core import JoinableLinePiece as _LinePiece, is_mr_cluster_endpoint as _is_cluster_endpoint


# ---------------------------------------------------------------------------
# Continuum presence -- the bisection gate for plain ZeroManager callers
# ---------------------------------------------------------------------------

def _check_boundary_indices(poly: CharPoly) -> None:
    M = poly.M
    K = poly.M + poly.N
    if M >= K:
        raise ValueError(f"M={M} >= K={K}: no PMGBZ boundary to check")
    if M <= 0:
        raise ValueError(f"M={M} <= 0: invalid boundary index")


@live_defaults(continuum_tol="core:CONTINUUM_TOL")
def detect_continuum_simple(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    continuum_tol: Optional[float] = None,
) -> bool:
    """Presence-only continuum detection.

    Builds a fresh ``Mu2MidZM``, runs ``analyze``, and returns
    ``has_continuum``.
    """
    _check_boundary_indices(poly)
    m = Mu2MidZM(zm.poly, zm.E_ref, zm.mu1)
    m.run()
    m.analyze(tie_tol=continuum_tol)
    return m.has_continuum


# ---------------------------------------------------------------------------
# LineSubset materialization
# ---------------------------------------------------------------------------

def _item_columns(
    zm: Mu2MidZM, s_idx: int, view, item: int,
) -> np.ndarray:
    """Real tracked_roots columns belonging to ItemView item *item*."""
    rep = int(view.rep_cols[item])
    seg_clusters = (zm._continuum_clusters[s_idx]
                    if s_idx < len(zm._continuum_clusters) else [])
    for cl in seg_clusters:
        if rep in cl:
            return np.array(list(cl), dtype=int)
    return np.array([rep], dtype=int)


def _find_boundary_runs(
    zm: Mu2MidZM,
) -> list[tuple[int, int, np.ndarray, int, int]]:
    """Contiguous LineSubset runs delimited by MRs, seams and events.

    Event rows with ``changes_boundary == True`` terminate a run.  Event rows
    that leave ``(j_lo, j_hi)`` unchanged do NOT split the run.
    """
    runs: list[tuple[int, int, np.ndarray, int, int]] = []
    for s_idx, seg in enumerate(zm.segments):
        n = len(seg.theta1_arr)
        if n == 0:
            continue
        cuts = {0, n - 1}
        for g in zm._event_groups:
            if g.seg_idx == s_idx and g.changes_boundary and 0 <= g.row < n:
                cuts.add(g.row)
        ordered = sorted(cuts)
        view = zm._item_views[s_idx]
        for a, b in zip(ordered[:-1], ordered[1:]):
            if b - a < 1:
                continue  # no regular row between two terminators → no run
            rep = a + 1 if a + 1 < b else a
            j_lo = int(view.j_lo[rep])
            j_hi = int(view.j_hi[rep])
            if j_lo != j_hi:
                continue
            cols = _item_columns(zm, s_idx, view, j_lo)
            runs.append((s_idx, j_lo, cols, int(a), int(b)))
    return runs


def _runs_to_pieces(
    zm: Mu2MidZM,
    runs: list[tuple[int, int, np.ndarray, int, int]],
) -> list[_LinePiece]:
    """One ``_LinePiece`` per run per continuum track."""
    pieces: list[_LinePiece] = []
    for s_idx, _, cols, a, b in runs:
        seg = zm.segments[s_idx]
        th = seg.theta1_arr[a:b + 1]
        tr = seg.tracked_roots[a:b + 1]
        for j in cols:
            pieces.append(_LinePiece(
                E=zm.E_ref, mu1=zm.mu1,
                theta1_arr=th.copy(), beta2_arr=tr[:, int(j)].copy(),
                ml=s_idx, mr=s_idx,
            ))
    return pieces


def _join_runs_across_mrs(
    zm: Mu2MidZM,
    pieces: list[_LinePiece],
) -> list[_LinePiece]:
    """Join pieces whose endpoints touch a segment boundary (MR / seam).

    Interior event endpoints cannot join across segments. Boundary-changing
    events on the seam are also excluded explicitly. At an MR a track ends if
    its endpoint track column is in the MR's cluster; otherwise it continues.
    Interior MR boundaries share one track frame, while the first/last-segment
    seam is translated through ``boundary_perm`` — no nearest-value root
    matching is used for cluster membership or track continuation.
    """
    n_seg = len(zm.segments)
    if n_seg <= 1:
        return pieces

    # A boundary-changing event is a hard LineSubset terminator, including
    # when it sits on the θ=0≡2π seam row.
    event_rows_changed = {
        (g.seg_idx, g.row)
        for g in getattr(zm, '_event_groups', [])
        if g.changes_boundary and g.row >= 0
    }

    def _mod(k: int) -> int:
        return k % n_seg

    def find_by_left(seg_s: int, root: complex) -> int | None:
        for idx, p in enumerate(pieces):
            if p.ml == seg_s and np.abs(p.beta2_arr[0] - root) < 1e-9:
                return idx
        return None

    def find_by_right(seg_s: int, root: complex) -> int | None:
        for idx, p in enumerate(pieces):
            if p.mr == seg_s and np.abs(p.beta2_arr[-1] - root) < 1e-9:
                return idx
        return None

    def _match_column_right(cur_s: int, prev_s: int, col: int) -> int:
        # Interior MR rows share one track frame.  The only frame change is
        # the cyclic first/last-segment seam, where the documented monodromy
        # is roots_right[boundary_perm] == roots_left.
        if cur_s == 0 and prev_s == n_seg - 1:
            return int(zm.boundary_perm[int(col)])
        return int(col)

    # Restart the scan after every merge: a merge invalidates the piece
    # indices (the list is rebuilt in place), so a plain nested loop would
    # read stale positions.  `while changed` (NOT a fixed n_seg budget —
    # the required number of merges scales with continuum multiplicity ×
    # segment boundaries, not with n_seg; a capped budget left high-fold
    # continua as broken open pieces) is safe: each merge strictly
    # decreases the piece count, so the loop terminates.
    while True:
        changed = False
        for s in range(n_seg):
            seg = zm.segments[s]
            prev = _mod(s - 1 + n_seg)
            prev_seg = zm.segments[prev]
            left_mr = seg.left_mr
            right_mr = prev_seg.right_mr
            is_circle_seam = (left_mr < 0 and right_mr < 0)
            if not is_circle_seam and left_mr != right_mr:
                continue
            if (s, 0) in event_rows_changed:
                continue
            if (prev, len(prev_seg.theta1_arr) - 1) in event_rows_changed:
                continue

            left_b = seg.tracked_roots[0, :]

            for j_l in range(zm.K):
                root = complex(left_b[j_l])
                if _is_cluster_endpoint(zm, seg, 'left', int(j_l)):
                    continue
                li_idx = find_by_left(s, root)
                if li_idx is None:
                    continue
                j_prev = _match_column_right(s, prev, int(j_l))
                root_prev = complex(prev_seg.tracked_roots[-1, j_prev])
                if _is_cluster_endpoint(zm, prev_seg, 'right', j_prev):
                    raise ValueError(
                        f"Continuum track {j_l} of segment {s} ends at the MR "
                        f"at θ₁={seg.theta1_arr[0]:.4f} as a non-cluster root, "
                        f"but the matched root (track {j_prev}) of segment "
                        f"{prev} is a cluster root there."
                    )
                pi_idx = find_by_right(prev, root_prev)
                if pi_idx is None:
                    continue
                if pi_idx == li_idx:
                    continue
                # cyclic=True wraps the segment-0 side past 2π so the merged
                # θ array stays monotonic across the seam.
                _merge_two(pieces, pi_idx, li_idx, cyclic=(s == 0))
                changed = True
                break
            if changed:
                break
        if not changed:
            break
    return pieces


def _merge_two(
    line_pieces: list[_LinePiece],
    prev_idx: int, cur_idx: int, *, cyclic: bool,
) -> None:
    """Merge ``line_pieces[prev_idx]`` (right side) with ``[cur_idx]`` (left)."""
    rp = line_pieces[prev_idx]
    cp = line_pieces[cur_idx]
    th = np.concatenate([rp.theta1_arr, cp.theta1_arr[1:]])
    b2 = np.concatenate([rp.beta2_arr, cp.beta2_arr[1:]])
    if cyclic:
        twopi = 2.0 * float(np.pi)
        n_rp = len(rp.theta1_arr)
        th[n_rp:] += twopi
    merged = _LinePiece(
        E=cp.E, mu1=cp.mu1, theta1_arr=th, beta2_arr=b2,
        ml=rp.ml, mr=cp.mr,
    )
    keep = [i for i in range(len(line_pieces)) if i not in (prev_idx, cur_idx)]
    line_pieces[:] = [line_pieces[i] for i in keep] + [merged]


def extract_continuum_linesubsets(
    zm: ZeroManager, poly: CharPoly,
) -> list[LineSubset]:
    """Materialise the 1D continuum LineSubsets of *zm*.

    Precondition: continuum detection has already run and returned
    ``has_continuum == True``.
    """
    m = ensure_mu2mid(zm)
    if not m.has_continuum:
        raise RuntimeError(
            "extract_continuum_linesubsets requires has_continuum=True; "
            "run continuum detection and gate on it before materializing "
            "LineSubsets."
        )

    runs = _find_boundary_runs(m)
    if not runs:
        raise RuntimeError(
            "has_continuum=True but _find_boundary_runs found no boundary "
            "continuum runs: continuum detection invariant broken."
        )

    pieces = _runs_to_pieces(m, runs)
    pieces = _join_runs_across_mrs(m, pieces)
    return [
        LineSubset(E=p.E, mu1=p.mu1,
                   theta1_arr=p.theta1_arr, beta2_arr=p.beta2_arr)
        for p in pieces
    ]
