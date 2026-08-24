'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-16
Copyright © Department of Physics, Tsinghua University. All rights reserved

Pairwise intersection of ItemView representative modulus curves.

The post-``ZeroManager.run`` crossing detector.  On every segment's ItemView
it compares representative item ``ln|β₂|`` curves pairwise:

  * ``d[:-1] == 0``             -- exact touch on the left mesh point
    (left-closed / right-open, so the right endpoint belongs to the next
    interval and the θ=0/2π seam is counted once);
  * ``d[:-1] * d[1:] < 0``      -- transversal sign change.

Only sign-change events are refined: a linear root prediction seeds
``scipy.optimize.brentq`` on the true ``ln|β₂_a| − ln|β₂_b|`` function.
There is NO θ deduplication -- every sign-change interval produces exactly
one refined point.  Close refined θ's (gap < ``crossing_tol``) are merged
into one *EventGroup* mesh row; the group keeps every member event and the
item/column connectivity they induce, but no numerical modulus check is ever
used to decide cluster membership.

MR boundary rows are NOT special-cased: a modulus coincidence there is
detected by the same touch/cross scan as anywhere else (``is_mr=True`` on
the event marks it, and its topological charge becomes ``None`` -- the
divergent cluster tangents make any direction unreliable).  There is no
separate MR→PointSubset materialization channel and no echo-drop rule.

Topological charge is NOT obtained by summing pair directions.  It is the
side change of the curve relative to the M-1/M cut:

    side(θ) = +1  if the item sits at sorted positions >= M
            = -1  if the item sits at sorted positions <  M

    charge = (side_right − side_left) / 2

so a single event is the same rule specialised to one bracket.
'''

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import brentq

from gbz_types import TWO_PI
from continuation import ZeroManager

if TYPE_CHECKING:
    from .mu2mid import Mu2MidZM

# Minimum |Re(V_a) - Re(V_b)| below which a pair crossing is treated as a
# tangent touch (hard boundary): the crossing direction is not trustworthy.
MIN_DIRECTION_DERIV: float = 1e-12

# brentq extra iteration budget (the bracket is guaranteed by the sign scan).
_BRENTQ_MAXITER: int = 100

# Real-root filter for nothing here (brentq operates on the true function);
# this tolerance only guards exact-endpoint float comparison.
_THETA_EQ_TOL: float = 1e-15


# ---------------------------------------------------------------------------
# Event data
# ---------------------------------------------------------------------------

@dataclass
class PairEvent:
    """One pairwise intersection of two ItemView representative items.

    *cols_a* / *cols_b* are the REAL tracked_roots columns represented by the
    two items -- the representative columns are only the refinement handle.
    PointSubset materialisation must expand the whole cluster, never just the
    representative.
    """
    seg_idx: int
    i: int                       # left mesh interval index [θ_i, θ_{i+1})
    ia: int                      # ItemView item index, left operand
    ib: int                      # ItemView item index, right operand
    cols_a: tuple[int, ...]      # real columns behind item ia
    cols_b: tuple[int, ...]      # real columns behind item ib
    rep_a: int                   # representative column used for refinement
    rep_b: int
    kind: str                    # 'touch' | 'cross'
    pair_kind: str               # 'M-2_M-1' | 'M-1_M' | 'M_M+1' | 'multi'
    theta_star: float
    direction: int | None        # sign(Re(V_a) - Re(V_b)), None = tangent
    is_mr: bool = False          # event sits on a segment MR boundary row
    converged: bool = True


@dataclass
class EventGroup:
    """One event mesh row -- the ONLY unit consumed by later extraction.

    A group may contain a single event (``len(events) == 1``) or several
    close events merged into one θ; the downstream logic is identical for
    both cases.
    """
    seg_idx: int
    theta: float
    events: tuple[PairEvent, ...]
    # Connectivity induced by the member events (nodes = ITEMS).
    item_components: tuple[tuple[int, ...], ...] = ()
    # Same components expanded to real tracked_roots columns.
    column_components: tuple[tuple[int, ...], ...] = ()
    # Any member event sits on a segment MR boundary row.
    is_mr: bool = False

    # Filled after mesh insertion and ItemView rebuild.
    row: int = -1
    boundary_pair_left: tuple[int, int] | None = None   # (j_lo, j_hi) items
    boundary_pair_right: tuple[int, int] | None = None
    changes_boundary: bool = False
    point_columns: tuple[int, ...] = ()
    column_q: dict[int, int | None] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# ItemView helpers (duck-typed on Mu2MidZM)
# ---------------------------------------------------------------------------

def _item_columns(zm, s_idx: int, item: int) -> tuple[int, ...]:
    """Real tracked_roots columns of an ItemView item on one segment."""
    view = zm._item_views[s_idx]
    rep = int(view.rep_cols[item])
    clusters = zm._continuum_clusters[s_idx]
    for cl in clusters:
        if rep in cl:
            return tuple(int(c) for c in cl)
    return (rep,)


def _item_positions(view, row: int, item: int) -> np.ndarray:
    """Sorted positions occupied by *item* in one ItemView row."""
    return np.flatnonzero(view.sort_to_item[row] == int(item))


def _pair_kind_from_row(view, row: int, ia: int, ib: int, M: int) -> str:
    """M-2_M-1 / M-1_M / M_M+1 classification from a REGULAR row.

    Never call this on a tie row: the classification reads the boundary
    between the two items' occupied sort positions.
    """
    pos_a = _item_positions(view, row, ia)
    pos_b = _item_positions(view, row, ib)
    if len(pos_a) == 0 or len(pos_b) == 0:
        return 'multi'
    boundary = min(int(pos_a[-1]), int(pos_b[-1]))
    if boundary == M - 2:
        return 'M-2_M-1'
    if boundary == M - 1:
        return 'M-1_M'
    if boundary == M:
        return 'M_M+1'
    return 'multi'


def _touch_pair_kind(view, seg, i: int, ia: int, ib: int, M: int) -> str:
    """Classification of a touch event from the nearest non-tie side row."""
    n = len(seg.theta1_arr)
    d = (view.item_logabs[:, ia] - view.item_logabs[:, ib])
    for row in (i - 1, i + 1):
        if 0 <= row < n and d[row] != 0.0:
            return _pair_kind_from_row(view, row, ia, ib, M)
    # Degenerate: the two items are equal on every accessible row; the event
    # is a continuum-adjacent touch and must not be classified as a point.
    return 'multi'


def _direction_from_tangents(
    seg, row: int, rep_a: int, rep_b: int, min_direction_deriv: float,
) -> int | None:
    """Pair derivative difference with minimum-derivative protection."""
    gp = float(seg.tangents[row, rep_a].real - seg.tangents[row, rep_b].real)
    if not np.isfinite(gp) or abs(gp) < min_direction_deriv:
        return None
    return 1 if gp > 0 else -1


def _normalize_theta(theta: float) -> float:
    """Normalize to [0, 2π)."""
    return float(theta % (TWO_PI))


# ---------------------------------------------------------------------------
# Event collection
# ---------------------------------------------------------------------------

def collect_pair_events(
    zm: Mu2MidZM,
    *,
    crossing_tol: float = 1e-10,
    min_direction_deriv: float = MIN_DIRECTION_DERIV,
) -> list[PairEvent]:
    """Scan every representative item pair on every segment.

    Zero deduplication: one sign-change interval → one event; one exact touch
    row → one event.  ``touch`` events are exact and are never refined.
    """
    events: list[PairEvent] = []
    for s_idx, seg in enumerate(zm.segments):
        th = seg.theta1_arr
        n = len(th)
        if n < 2:
            raise ValueError(f"Segment {s_idx} has less than 2 items. E_ref={zm.E_ref}, mu1={zm.mu1}")
        view = zm._item_views[s_idx]
        n_items = len(view.rep_cols)
        for ia in range(n_items):
            cols_a = _item_columns(zm, s_idx, ia)
            rep_a = int(view.rep_cols[ia])
            for ib in range(ia + 1, n_items):
                cols_b = _item_columns(zm, s_idx, ib)
                rep_b = int(view.rep_cols[ib])
                d = view.item_logabs[:, ia] - view.item_logabs[:, ib]
                # A non-finite d row (0/∞ padding root → log|β₂| = ±∞)
                # poisons only the mesh intervals it touches, not the whole
                # pair: a transient singular root at one θ (a
                # degree-deficient solve later re-matched to a real root)
                # must not silently discard this pair's genuine crossings
                # elsewhere on the segment.  Persistent padding columns
                # still yield zero valid intervals, so the historical
                # whole-pair skip is the degenerate case of this scan.
                # touch needs no finite mask: d[i] == 0.0 already implies
                # both moduli are finite there (±∞ − ±∞ is nan, never 0).
                finite = np.isfinite(d)
                touch = np.flatnonzero(d[:-1] == 0.0)
                cross = np.flatnonzero(
                    finite[:-1] & finite[1:] & (d[:-1] * d[1:] < 0.0))

                for i in touch:
                    # A touch landing on a segment's MR boundary row is an
                    # MR event, not a skipped one: the MR channel used to
                    # own it, but that channel only materializes clusters
                    # straddling M-1/M, so any other modulus coincidence at
                    # the boundary row was silently dropped (the seam
                    # missed-detection at E=1.212).  It flows through the
                    # same EventGroup machinery; only its charge becomes
                    # None (is_mr=True).
                    is_mr = ((i == 0 and seg.left_mr >= 0)
                             or (i == n - 1 and seg.right_mr >= 0))
                    events.append(PairEvent(
                        seg_idx=s_idx, i=int(i), ia=ia, ib=ib,
                        cols_a=cols_a, cols_b=cols_b,
                        rep_a=rep_a, rep_b=rep_b,
                        kind='touch',
                        pair_kind=_touch_pair_kind(
                            view, seg, int(i), ia, ib, zm.M),
                        theta_star=float(th[i]),
                        direction=_direction_from_tangents(
                            seg, int(i), rep_a, rep_b, min_direction_deriv),
                        is_mr=is_mr,
                        converged=True,
                    ))

                for i in cross:
                    theta_star, direction, converged = _refine_pair_crossing(
                        zm, s_idx, int(i), rep_a, rep_b,
                        float(th[i]), float(th[i + 1]),
                        crossing_tol=crossing_tol,
                        min_direction_deriv=min_direction_deriv,
                    )
                    events.append(PairEvent(
                        seg_idx=s_idx, i=int(i), ia=ia, ib=ib,
                        cols_a=cols_a, cols_b=cols_b,
                        rep_a=rep_a, rep_b=rep_b,
                        kind='cross',
                        pair_kind=_pair_kind_from_row(
                            view, int(i), ia, ib, zm.M),
                        theta_star=theta_star,
                        direction=direction,
                        converged=converged,
                    ))
    return events


def _refine_pair_crossing(
    zm: ZeroManager,
    s_idx: int,
    i: int,
    rep_a: int,
    rep_b: int,
    theta_lo: float,
    theta_hi: float,
    *,
    crossing_tol: float,
    min_direction_deriv: float,
) -> tuple[float, int | None, bool]:
    """brentq refinement of ``ln|β₂_rep_a| − ln|β₂_rep_b| = 0``.

    Endpoints are read straight off the mesh rows; interior points are solved
    transiently with ``ZeroManager.solve_at(..., interp='linear')`` and
    cached.  No mesh mutation happens here.
    """
    cache: dict[float, tuple[float, float]] = {}

    def eval_point(t: float) -> tuple[float, float]:
        t = float(t)
        if t in cache:
            return cache[t]
        seg = zm.segments[s_idx]
        th = seg.theta1_arr
        if abs(t - float(th[i])) < _THETA_EQ_TOL:
            roots = seg.tracked_roots[i, :]
            V = seg.tangents[i, :]
        elif abs(t - float(th[i + 1])) < _THETA_EQ_TOL:
            roots = seg.tracked_roots[i + 1, :]
            V = seg.tangents[i + 1, :]
        else:
            roots, V = zm.solve_at(t, seg_idx=s_idx, i=i, interp='linear')
        la = np.log(np.abs(roots))
        gp = float(V[rep_a].real - V[rep_b].real)
        out = (float(la[rep_a] - la[rep_b]), gp)
        cache[t] = out
        return out

    f0, gp0 = eval_point(theta_lo)
    f1, gp1 = eval_point(theta_hi)
    # The caller guarantees a sign change; an exact endpoint root is a touch.
    if f0 == 0.0:
        return theta_lo, _protected_direction(gp0, min_direction_deriv), True
    if f1 == 0.0:
        return theta_hi, _protected_direction(gp1, min_direction_deriv), True

    theta_lin = theta_lo + (0.0 - f0) / (f1 - f0) * (theta_hi - theta_lo)
    try:
        root = brentq(
            lambda t: eval_point(t)[0],
            theta_lo, theta_hi,
            xtol=crossing_tol, rtol=4.0 * np.finfo(float).eps,
            maxiter=_BRENTQ_MAXITER,
        )
        _, gp = eval_point(root)
        return float(root), _protected_direction(gp, min_direction_deriv), True
    except Exception:
        warnings.warn(
            f"pairwise refinement failed in segment {s_idx} interval "
            f"[{theta_lo:.6e}, {theta_hi:.6e}] pair ({rep_a},{rep_b}); "
            f"falling back to the linear prediction {theta_lin:.6e}"
        )
        secant = (f1 - f0) / (theta_hi - theta_lo)
        return theta_lin, _protected_direction(secant, min_direction_deriv), False


def _protected_direction(gp: float, min_direction_deriv: float) -> int | None:
    if not np.isfinite(gp) or abs(gp) < min_direction_deriv:
        return None
    return 1 if gp > 0 else -1


# ---------------------------------------------------------------------------
# Event grouping (close θ* → one EventGroup)
# ---------------------------------------------------------------------------

def group_events(
    zm: ZeroManager,
    events: list[PairEvent],
    *,
    merge_tol: float,
) -> list[EventGroup]:
    """Merge events whose θ* are closer than *merge_tol*.

    A merged group lands on an EXISTING mesh data point whenever its member
    θ's straddle one; a seam merge lands on 0 (≡ 2π).  Member events are
    always kept inside the group.
    """
    groups: list[EventGroup] = []
    by_seg: dict[int, list[PairEvent]] = {}
    for ev in events:
        by_seg.setdefault(ev.seg_idx, []).append(ev)

    for s_idx, evs in sorted(by_seg.items()):
        th = np.asarray(zm.segments[s_idx].theta1_arr, dtype=float)
        evs = sorted(evs, key=lambda e: e.theta_star)
        seg_groups: list[EventGroup] = []
        for ev in evs:
            if (seg_groups
                    and ev.theta_star - seg_groups[-1].events[-1].theta_star
                    < merge_tol):
                seg_groups[-1].events += (ev,)
            else:
                seg_groups.append(EventGroup(
                    seg_idx=s_idx, theta=ev.theta_star, events=(ev,),
                ))

        for g in seg_groups:
            _fill_group_components(g)
            if len(g.events) == 1:
                g.theta = g.events[0].theta_star
                continue
            thetas = [e.theta_star for e in g.events]
            lo, hi = min(thetas), max(thetas)
            # A multi-event merge lands on an existing mesh data point when
            # the events straddle one; otherwise it keeps the midpoint.
            inside = th[(th >= lo) & (th <= hi)]
            if inside.size:
                mid = float((lo + hi) / 2.0)
                g.theta = float(inside[int(np.argmin(np.abs(inside - mid)))])
            else:
                g.theta = float((lo + hi) / 2.0)

        # θ=0 and θ=2π are the same physical point: merge through the seam
        # and place the group on 0.
        if len(seg_groups) >= 2:
            first, last = seg_groups[0], seg_groups[-1]
            if _circ_gap(first.theta, last.theta) < merge_tol:
                # Events detected on the 2π side carry segment track-frame
                # column labels, which differ from the θ=0 frame by
                # boundary_perm (the closing row stores left_boundary_roots
                # permuted into the track frame).  The merged group is
                # anchored at θ=0, where finalize reads every member's
                # labels as frame-0 — so relabel ONLY events genuinely
                # detected on the 2π side (last.theta > π); events detected
                # near 0⁺ are already frame-0 and must pass through
                # untouched.
                if last.theta > math.pi:
                    seam_events = tuple(
                        _relabel_event_to_theta0_frame(zm, s_idx, e)
                        for e in last.events
                    )
                else:
                    seam_events = last.events
                merged = EventGroup(
                    seg_idx=s_idx, theta=0.0,
                    events=seam_events + first.events,
                )
                _fill_group_components(merged)
                seg_groups = [merged] + seg_groups[1:-1]

        groups.extend(seg_groups)
    return groups


def _relabel_event_to_theta0_frame(
    zm: ZeroManager, s_idx: int, ev: PairEvent,
) -> PairEvent:
    """Translate one 2π-side event's identity labels to the θ=0 frame.

    The closing row of a full-circle segment stores ``left_boundary_roots``
    permuted into the track frame (run(): ``left_boundary_roots[perm]``), so
    track column j at the 2π end carries the θ=0 root ``inv_perm[j]`` — the
    same physical root that track ``inv_perm[j]`` owns at row 0.  Without
    this translation a seam-merged event names its pair by 2π-side labels
    while the group lives in the θ=0 frame, so finalize's side/charge reads
    hit the wrong tracks and the boundary coverage test silently drops the
    event.  ``theta_star`` is kept as detected (provenance); only the
    identity labels move.
    """
    inv_perm = np.empty(zm.K, dtype=int)
    inv_perm[zm.boundary_perm] = np.arange(zm.K)

    view = zm._item_views[s_idx]
    rep_a = int(inv_perm[ev.rep_a])
    rep_b = int(inv_perm[ev.rep_b])
    ia = int(np.flatnonzero(view.rep_cols == rep_a)[0])
    ib = int(np.flatnonzero(view.rep_cols == rep_b)[0])
    return replace(
        ev,
        ia=ia, ib=ib,
        cols_a=tuple(int(inv_perm[c]) for c in ev.cols_a),
        cols_b=tuple(int(inv_perm[c]) for c in ev.cols_b),
        rep_a=rep_a, rep_b=rep_b,
    )


def _circ_gap(a: float, b: float) -> float:
    twopi = TWO_PI
    d = abs((a - b) % twopi)
    return min(d, twopi - d)


def _fill_group_components(g: EventGroup) -> None:
    """Item/column connectivity induced by member events (no modulus check)."""
    g.is_mr = any(e.is_mr for e in g.events)
    item_edges = [(e.ia, e.ib) for e in g.events]
    item_components = _connected_components(item_edges)
    g.item_components = tuple(tuple(c) for c in item_components)

    col_components: list[tuple[int, ...]] = []
    # The item graph components already cover every event column; expand each
    # item component through the events' stored real-column sets.  This keeps
    # the "no extra modulus check" rule explicit.
    for comp in item_components:
        cols: list[int] = []
        for e in g.events:
            if e.ia in comp:
                for c in e.cols_a:
                    if c not in cols:
                        cols.append(c)
            if e.ib in comp:
                for c in e.cols_b:
                    if c not in cols:
                        cols.append(c)
        col_components.append(tuple(sorted(cols)))
    g.column_components = tuple(col_components)


def _connected_components(edges: list[tuple[int, int]]) -> list[list[int]]:
    parent: dict[int, int] = {}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for a, b in edges:
        if a not in parent:
            parent[a] = a
        if b not in parent:
            parent[b] = b
        union(a, b)

    comps: dict[int, list[int]] = {}
    for x in parent:
        comps.setdefault(find(x), []).append(x)
    return [sorted(v) for v in comps.values()]


def insert_event_groups(
    zm: ZeroManager,
    groups: list[EventGroup],
) -> None:
    """Insert one mesh row per EventGroup (descending θ per segment)."""
    by_seg: dict[int, list[EventGroup]] = {}
    for g in groups:
        by_seg.setdefault(g.seg_idx, []).append(g)

    for s_idx, gs in sorted(by_seg.items()):
        seam_groups: list[EventGroup] = []
        for g in sorted(gs, key=lambda x: x.theta, reverse=True):
            theta_norm = _normalize_theta(g.theta)
            if theta_norm == 0.0 and g.theta > math.pi:
                # θ* collapsed onto the 2π seam: defer row assignment until
                # every interior insertion has happened (indices shift).
                seam_groups.append(g)
                continue
            zm.insert_solution(g.theta, seg_idx=s_idx, interp='linear')
        for g in seam_groups:
            last = len(zm.segments) - 1
            g.seg_idx = last
            g.row = len(zm.segments[last].theta1_arr) - 1
        # Row indices shift while later (smaller-θ) insertions happen, so
        # resolve every group's row against the FINAL mesh.  Identity (not
        # `in`) membership: EventGroup is a dataclass whose value-equality
        # could match a DIFFERENT group that happens to carry equal fields
        # (symmetric models); the seam groups were already assigned rows
        # above and must be skipped, not re-resolved.
        seam_ids = {id(g) for g in seam_groups}
        for g in gs:
            if id(g) in seam_ids:
                continue
            g.row = _find_mesh_row(zm, g.seg_idx, g.theta)


def _find_mesh_row(zm: ZeroManager, seg_idx: int, theta: float) -> int:
    """Exact-float row index of an inserted/refined θ in the final mesh."""
    th = zm.segments[seg_idx].theta1_arr
    theta_norm = _normalize_theta(theta)
    idx = int(np.argmin(np.abs(th - theta_norm)))
    if abs(float(th[idx]) - theta_norm) > 1e-14:
        raise RuntimeError(
            f"event θ={theta} (normalized {theta_norm}) not found in final "
            f"mesh of segment {seg_idx}"
        )
    return idx


def finalize_event_groups(zm: Mu2MidZM, groups: list[EventGroup]) -> None:
    """Fill boundary-pair change, point columns and side-change charges.

    Must be called AFTER the event rows are inserted and ItemView rebuilt.
    Sort information comes only from the nearest REGULAR rows; event rows are
    never used for sorting.  The θ=0/2π seam uses the documented convention
    ``roots_right[boundary_perm] == roots_left``.
    """
    event_rows: set[tuple[int, int]] = {
        (g.seg_idx, g.row) for g in groups if g.row >= 0
    }
    M = zm.M
    inv_perm = np.empty(zm.K, dtype=int)
    inv_perm[zm.boundary_perm] = np.arange(zm.K)

    for g in groups:
        s_idx, row = g.seg_idx, g.row
        if row < 0:
            continue
        n = len(zm.segments[s_idx].theta1_arr)

        left = _regular_side(zm, event_rows, s_idx, row, n, -1)
        right = _regular_side(zm, event_rows, s_idx, row, n, +1)
        g.boundary_pair_left = _item_pair(zm, s_idx, row, left, inv_perm)
        g.boundary_pair_right = _item_pair(zm, s_idx, row, right, inv_perm)
        g.changes_boundary = (
            g.boundary_pair_left is not None
            and g.boundary_pair_right is not None
            and g.boundary_pair_left != g.boundary_pair_right
        )

        for comp_idx, col_comp in enumerate(g.column_components):
            positions_left = _component_positions(
                zm, s_idx, row, col_comp, left, inv_perm)
            positions_right = _component_positions(
                zm, s_idx, row, col_comp, right, inv_perm)
            if not ({M - 1, M} <= positions_left
                    or {M - 1, M} <= positions_right):
                continue

            hard = any(
                e.direction is None
                for e in g.events
                if (comp_idx < len(g.item_components)
                    and (e.ia in g.item_components[comp_idx]
                         or e.ib in g.item_components[comp_idx]))
            )
            for c in sorted(col_comp):
                if hard:
                    q = None
                else:
                    # μ₂_mid passes between sort positions M-1 and M:
                    # side = -1 below it (j ≤ M-1), +1 above it (j ≥ M).
                    side_left = _side_of_column(
                        zm, s_idx, row, c, left, inv_perm, M)
                    side_right = _side_of_column(
                        zm, s_idx, row, c, right, inv_perm, M)
                    q = int((side_right - side_left) / 2)
                g.column_q.setdefault(c, q)
                if c not in g.point_columns:
                    g.point_columns += (c,)
        g.point_columns = tuple(sorted(g.point_columns))


def _regular_side(
    zm: Mu2MidZM,
    event_rows: set[tuple[int, int]],
    s_idx: int,
    row: int,
    n: int,
    step: int,
) -> tuple[int, int, int | None] | None:
    """Nearest regular row on one side.

    Returns ``(seg_idx, row, frame)``; *frame* is 0 on the θ=0 side, 1 on
    the θ=2π side, and None for an ordinary same-frame row.
    """
    r = _nearest_regular_row(event_rows, s_idx, row, n, step)
    if r is not None:
        return s_idx, r, None

    n_seg = len(zm.segments)
    if step < 0 and s_idx == 0 and row == 0:
        other = n_seg - 1
        oth_n = len(zm.segments[other].theta1_arr)
        r = _nearest_regular_row(event_rows, other, oth_n - 1, oth_n, -1)
        return (other, r, 1) if r is not None else None
    if step > 0 and s_idx == n_seg - 1 and row == n - 1:
        r = _nearest_regular_row(
            event_rows, 0, 0, len(zm.segments[0].theta1_arr), +1)
        return (0, r, 0) if r is not None else None
    return None


def _group_frame(zm: Mu2MidZM, s_idx: int, row: int) -> int | None:
    """0 if this row is the θ=0 seam side, 1 if the θ=2π seam side."""
    n_seg = len(zm.segments)
    if s_idx == 0 and row == 0:
        return 0
    if s_idx == n_seg - 1 and row == len(zm.segments[s_idx].theta1_arr) - 1:
        return 1
    return None


def _map_col_across_seam(
    zm: Mu2MidZM,
    col: int,
    src_frame: int | None,
    dst_frame: int | None,
    inv_perm: np.ndarray,
) -> int:
    """Map a column index between the θ=0 and θ=2π track frames.

    Convention (ZeroManager docstring): ``roots_right[boundary_perm] ==
    roots_left``, i.e. left column j ↔ right column ``boundary_perm[j]``.
    """
    if src_frame is None or dst_frame is None or src_frame == dst_frame:
        return int(col)
    if src_frame == 0 and dst_frame == 1:
        return int(zm.boundary_perm[col])
    if src_frame == 1 and dst_frame == 0:
        return int(inv_perm[col])
    raise ValueError(f"unexpected frame transition {src_frame} → {dst_frame}")


def _item_pair(
    zm: Mu2MidZM,
    group_seg: int,
    group_row: int,
    side: tuple[int, int, int | None] | None,
    inv_perm: np.ndarray,
) -> tuple[int, int] | None:
    """``(j_lo, j_hi)`` items at *side*, expressed in the group's frame."""
    if side is None:
        return None
    src_seg, src_row, src_frame = side
    view_src = zm._item_views[src_seg]
    if src_seg == group_seg and src_frame is None:
        return int(view_src.j_lo[src_row]), int(view_src.j_hi[src_row])

    dst_frame = _group_frame(zm, group_seg, group_row)
    cols = [
        int(view_src.rep_cols[view_src.j_lo[src_row]]),
        int(view_src.rep_cols[view_src.j_hi[src_row]]),
    ]
    mapped = [
        _map_col_across_seam(zm, c, src_frame, dst_frame, inv_perm)
        for c in cols
    ]
    return tuple(_item_of_column(zm, group_seg, c) for c in mapped)


def _component_positions(
    zm: Mu2MidZM,
    group_seg: int,
    group_row: int,
    col_comp: tuple[int, ...],
    side: tuple[int, int, int | None] | None,
    inv_perm: np.ndarray,
) -> set[int]:
    """Sort positions occupied by the component on one side."""
    if side is None:
        return set()
    src_seg, src_row, src_frame = side
    dst_frame = _group_frame(zm, group_seg, group_row)
    items = {
        _item_of_column(
            zm, src_seg,
            _map_col_across_seam(zm, c, dst_frame, src_frame, inv_perm),
        )
        for c in col_comp
    }
    return _positions_of_items(zm._item_views[src_seg], src_row, items)


def _side_of_column(
    zm: Mu2MidZM,
    group_seg: int,
    group_row: int,
    col: int,
    side: tuple[int, int, int | None] | None,
    inv_perm: np.ndarray,
    M: int,
) -> int:
    """-1 below μ₂_mid (j ≤ M-1), +1 above it (j ≥ M), 0 if straddling."""
    if side is None:
        return 0
    src_seg, src_row, src_frame = side
    dst_frame = _group_frame(zm, group_seg, group_row)
    mapped = _map_col_across_seam(zm, col, dst_frame, src_frame, inv_perm)
    item = _item_of_column(zm, src_seg, mapped)
    return _item_side(zm._item_views[src_seg], src_row, item, M)


def _nearest_regular_row(
    event_rows: set[tuple[int, int]],
    s_idx: int,
    row: int,
    n: int,
    step: int,
) -> int | None:
    """Nearest row in direction *step* that is not an event row."""
    r = row + step
    while 0 <= r < n:
        if (s_idx, r) not in event_rows:
            return r
        r += step
    return None


def _item_of_column(zm: Mu2MidZM, s_idx: int, col: int) -> int:
    """Rebuilt-ItemView item index owning a real column on this segment."""
    for item in range(len(zm._item_views[s_idx].rep_cols)):
        if col in _item_columns(zm, s_idx, item):
            return item
    raise ValueError(
        f"column {col} not found in rebuilt ItemView of segment {s_idx}"
    )


def _positions_of_items(view, row: int | None, items: set[int]) -> set[int]:
    if row is None:
        return set()
    out: set[int] = set()
    for item in items:
        out.update(int(p) for p in _item_positions(view, row, item))
    return out


def _item_side(view, row: int | None, item: int, M: int) -> int:
    """side relative to μ₂_mid: -1 below (j ≤ M-1), +1 above (j ≥ M)."""
    if row is None:
        return 0
    pos = _item_positions(view, row, item)
    if len(pos) == 0:
        return 0
    if int(pos[-1]) <= M - 1:
        return -1
    if int(pos[0]) >= M:
        return 1
    return 0
