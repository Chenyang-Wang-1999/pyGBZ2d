'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-13
Copyright © Department of Physics, Tsinghua University. All rights reserved

SGBZ continuum: detection + 1D LineSubset materialization.

A continuum arises when ``|β_M| = |β_{M+1}|`` holds identically over a θ₁
range (a 1D subset).  Continuum detection is folded into the μ₂_mid build
itself (:meth:`Mu2MidZM.build_mu2_mid` sets :attr:`has_continuum` — §1/§6.3):
the whole-segment same-modulus ItemView criterion has no false negatives, so
no separate two-point gate is needed.  This module provides:

  * :func:`detect_continuum_simple` — a presence-only flag for callers that
    hold a plain ``ZeroManager`` and only need "is there a continuum?" (the
    plateau probe).  The μ₁ bisection itself does NOT call this — it builds
    its own ``Mu2MidZM`` and reads ``has_continuum`` inline.
  * :func:`extract_continuum_linesubsets` — materialise the 1D continuum tracks
    as ``LineSubset``s.

SGBZ LineSubset semantics differ from amoeba's.  An SGBZ LineSubset is the
stretch where a specific **boundary pair** (sorted positions M-1/M) is the
degenerate continuum item — it terminates at EITHER:

  * a **modulus-sort change** (another root overtakes the boundary pair, so
    they swap out of M-1/M even though still same-modulus) — a hard terminator
    *inside* a segment, OR
  * a **multiple root** (MR) at a segment boundary — same as amoeba; the track
    ends there if it is in the MR's cluster, otherwise continues into the
    adjacent segment.

So a piece is a **contiguous run** of ``j_lo == j_hi == item`` rows (the
boundary pair IS this continuum item), NOT the whole segment.  Amoeba's
whole-segment join logic cannot be reused: it has no notion of a sort-change
terminator and would take the entire segment's θ₁ range.  Only the MR-side
cluster test (``_is_cluster_endpoint``) and the piece container
(``_LinePiece``) are reused from amoeba — both read only the shared
``ZeroManager`` fields present on ``Mu2MidZM``.
'''

from __future__ import annotations

import numpy as np

from gbz_types import CharPoly, LineSubset
from continuation import ZeroManager

from .crossings import _ensure_mu2mid
from .mu2mid import Mu2MidZM, CONTINUUM_TOL
from brute_force_amoeba.zm_extract import _LinePiece, _is_cluster_endpoint


# ---------------------------------------------------------------------------
# Continuum presence (§1) — the bisection gate for plain ZeroManager callers
# ---------------------------------------------------------------------------

def _check_boundary_indices(poly: CharPoly) -> None:
    """Validate the M-1/M boundary indices against the polynomial degree."""
    M = poly.M
    K = poly.M + poly.N
    if M >= K:
        raise ValueError(f"M={M} >= K={K}: no PMGBZ boundary to check")
    if M <= 0:
        raise ValueError(f"M={M} <= 0: invalid boundary index")


def detect_continuum_simple(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    continuum_tol: float = CONTINUUM_TOL,
    zm_run_kwargs: dict | None = None,
) -> bool:
    """Simplified continuous-modulus-equality detection — presence only.

    Builds a fresh ``Mu2MidZM`` from *zm*'s ``(poly, E_ref, mu1)`` (a plain
    ``ZeroManager`` cannot be mutated in place) and returns
    :attr:`Mu2MidZM.has_continuum`.  The whole-segment same-modulus ItemView
    criterion (cluster columns share one modulus across the whole segment)
    replaces the old two-point gate — no false negatives at sub-grid
    continua (§1).  The build cost is paid at most once per probe; the
    bisection path avoids it by constructing its own ``Mu2MidZM`` and
    reusing the build for crossing detection and winding.
    """
    _check_boundary_indices(poly)
    m = Mu2MidZM(zm.poly, zm.E_ref, zm.mu1)
    m.run(**(zm_run_kwargs or {}))
    m.build_mu2_mid(tie_tol=continuum_tol)
    return m.has_continuum


def _item_columns(
    zm: Mu2MidZM, s_idx: int, view, item: int,
) -> np.ndarray:
    """The tracked_roots columns belonging to item *item* of *view*.

    The cluster item's representative is ``view.rep_cols[item]``; the full
    cluster membership is recovered from ``zm._continuum_clusters[s_idx]``
    — the per-segment list of column-tuples built by ``build_mu2_mid``.
    Scoped to *this* segment: the same representative column could belong to
    different clusters in different segments (a continuum that re-clusters
    across an MR), and a global first-match search would return the wrong
    segment's cluster.
    """
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
    """Contiguous runs where a continuum item occupies the M-1/M boundary.

    Returns ``(seg, item, cols, row_a, row_b)`` per run — track ``j`` is on
    the SGBZ boundary (the degenerate item ``item`` holds sorted positions
    M-1 and M) for rows ``[row_a, row_b]`` inclusive.  A run ends when
    ``j_lo``/``j_hi`` change (a sort change swaps the boundary pair) or at
    the segment edge.  Multiple runs per segment are possible: the same
    continuum pair can drop out of the boundary and re-enter later.
    """
    runs: list[tuple[int, int, np.ndarray, int, int]] = []
    for s_idx, view in enumerate(zm._item_views):
        seg = zm.segments[s_idx]
        N = len(seg.theta1_arr)
        if N == 0 or not np.any(view.j_lo == view.j_hi):
            continue
        for item in range(len(view.mults)):
            if view.mults[item] < 2:
                continue
            in_boundary = (view.j_lo == view.j_hi) & (view.j_lo == item)
            rows = np.flatnonzero(in_boundary)
            if len(rows) == 0:
                continue
            cols = _item_columns(zm, s_idx, view, item)
            # split into maximal contiguous runs
            start = rows[0]
            prev = rows[0]
            for r in rows[1:]:
                if r != prev + 1:
                    runs.append((s_idx, item, cols, int(start), int(prev)))
                    start = r
                prev = r
            runs.append((s_idx, item, cols, int(start), int(prev)))
    return runs


def _runs_to_pieces(
    zm: Mu2MidZM,
    runs: list[tuple[int, int, np.ndarray, int, int]],
) -> list[_LinePiece]:
    """Build one ``_LinePiece`` per run per continuum track.

    Each run yields one piece per column in the continuum item's cluster
    (every column is a distinct β₂ curve at the same |β₂|).  The piece's
    θ₁/β₂ arrays span exactly ``[row_a, row_b]`` — the run, not the whole
    segment.
    """
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

    A piece's endpoint is joinable only if it sits on a segment edge — i.e.
    the run started at row 0 (left edge) or ended at the last row (right
    edge).  An endpoint strictly inside a segment is a sort-change
    terminator: the continuum pair dropped out of the boundary there, so the
    LineSubset ends — no join.

    At a segment edge that is an MR: if the endpoint root is in the MR's
    cluster the track terminates (a genuine LineSubset end); otherwise it
    continues into the adjacent segment and is matched by root value.  At the
    θ₁=0≡2π seam (both sides ``mr < 0``) the last segment's right end matches
    segment 0's left end via ``boundary_perm``.

    Iterated to a fixpoint so chains and the cyclic seam converge.
    """
    n_seg = len(zm.segments)
    if n_seg <= 1:
        return pieces

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

    def _match_root_right(prev_seg, root: complex) -> tuple[int, complex] | None:
        """Track index + root on prev_seg's right boundary matching *root* by value."""
        right_b = prev_seg.tracked_roots[-1, :]
        j_prev = int(np.argmin(np.abs(right_b - root)))
        return j_prev, complex(right_b[j_prev])

    for _ in range(n_seg):
        changed = False
        for s in range(n_seg):
            seg = zm.segments[s]
            prev = _mod(s - 1 + n_seg)
            prev_seg = zm.segments[prev]
            left_mr = seg.left_mr
            right_mr = prev_seg.right_mr
            is_circle_seam = (left_mr < 0 and right_mr < 0)
            if not is_circle_seam and left_mr != right_mr:
                continue  # not a shared boundary

            left_b = seg.tracked_roots[0, :]

            for j_l in range(zm.K):
                root = complex(left_b[j_l])
                if _is_cluster_endpoint(zm, seg, 'left', root):
                    continue  # terminates at the MR
                li_idx = find_by_left(s, root)
                if li_idx is None:
                    continue  # not the leftmost piece here, or no run touches
                j_prev, root_prev = _match_root_right(prev_seg, root)
                if _is_cluster_endpoint(zm, prev_seg, 'right', root_prev):
                    # ends on one side, cluster on the other: inconsistent
                    raise ValueError(
                        f"Continuum track {j_l} of segment {s} ends at the MR "
                        f"at θ₁={seg.theta1_arr[0]:.4f} as a non-cluster root, "
                        f"but the matched root (track {j_prev}) of segment "
                        f"{prev} is a cluster root there."
                    )
                pi_idx = find_by_right(prev, root_prev)
                if pi_idx is None:
                    # The matched track has no run touching prev's right edge:
                    # it ends at a sort-change inside prev, so this piece truly
                    # terminates at the MR — not an error, just no join.
                    continue
                if pi_idx == li_idx:
                    continue  # already joined (e.g. the cyclic seam)
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
    """Merge ``line_pieces[prev_idx]`` (right side) with ``[cur_idx]`` (left).

    Non-cyclic (interior MR): ``prev`` (segment s-1) is to the LEFT of ``cur``
    (segment s) in θ₁, so the array is ``[prev, cur[1:]]`` — θ₁ stays
    monotonic, the shared MR row (cur's first row) is dropped.

    Cyclic seam (θ₁=0 ≡ 2π): ``prev`` is the LAST segment (right end at 2π),
    ``cur`` is segment 0 (left end at 0).  Align ``rp[-1]`` (θ=2π) with
    ``cp[0]`` (θ=0) by putting rp first: ``[rp, cp[1:]]``.

    ``ml``/``mr`` track the leftmost/rightmost original segment spanned: the
    merged piece's left end is ``rp``'s left end (``rp.ml``), right end is
    ``cp``'s right end (``cp.mr``) — same for both branches.
    """
    rp = line_pieces[prev_idx]
    cp = line_pieces[cur_idx]
    th = np.concatenate([rp.theta1_arr, cp.theta1_arr[1:]])
    b2 = np.concatenate([rp.beta2_arr, cp.beta2_arr[1:]])
    if cyclic:
        # rp ends at θ=2π, cp continues from θ=0.  Unwrap the cp portion by
        # +2π so the merged θ₁ stays monotonically increasing (LineSubset
        # requires a monotonic theta1_arr).  The seam point (θ=2π ≡ 0) lands
        # inside the array where the β₂ curve is continuous.
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
    zm: ZeroManager, poly: CharPoly, *, zm_run_kwargs: dict | None = None,
) -> list[LineSubset]:
    """Materialise the 1D continuum LineSubsets of *zm*.

    Returns ``[]`` when *zm* is not a continuum (``has_continuum`` False).
    Otherwise: one ``LineSubset`` per continuum track per boundary run, joined
    across MR boundaries where the track passes through as a non-cluster root.
    Each LineSubset's θ₁ range is exactly where its track held the M-1/M
    boundary — sort-change terminators inside a segment end the piece, unlike
    the old whole-segment materialization.
    """
    m = _ensure_mu2mid(zm, poly, **(zm_run_kwargs or {}))
    if not m.has_continuum:
        return []

    runs = _find_boundary_runs(m)
    if not runs:
        return []

    pieces = _runs_to_pieces(m, runs)
    pieces = _join_runs_across_mrs(m, pieces)
    return [
        LineSubset(E=p.E, mu1=p.mu1,
                   theta1_arr=p.theta1_arr, beta2_arr=p.beta2_arr)
        for p in pieces
    ]
