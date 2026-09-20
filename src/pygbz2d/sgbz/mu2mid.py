'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-16
Copyright © Department of Physics, Tsinghua University. All rights reserved

SGBZ ItemView + simplified μ₂_mid path.

``ItemView`` keeps its original semantics: a per-segment representative-item
view whose clusters are exactly "columns that stay equal-modulus over the
whole segment".  ``j_lo == j_hi`` therefore means one continuum item occupies
both boundary positions M-1 and M -- a 1D LineSubset.  Isolated tie rows
(single-θ crossings) are NOT injected into ItemView; they remain different
items and only swap sort positions.

``Mu2Mid`` is the loop-winding path, deliberately independent of
``ZeroManager``.  It stores piecewise-smooth intervals with endpoint values
and derivatives; each interval is interpolated with
``continuation.interpolation.hermite_interp_poly`` so the derivative is
smooth inside an interval.  Values are bounded to ±14 AT BUILD TIME: every
interval's Hermite polynomial is split at its ±14 crossings, and subintervals
outside the band become constant ``v=±14, dv=0`` pieces.  No post-hoc clip
exists anywhere downstream.

``Mu2MidZM.analyze()`` runs the pairwise intersection pipeline
(:mod:`pygbz2d.sgbz.pairwise`) and then builds this path.
'''

from __future__ import annotations

import warnings

from dataclasses import dataclass
from typing import NamedTuple, Optional

import numpy as np

from ..continuation import ZeroManager
from ..continuation.interpolation import hermite_interp_poly

from .. import core
from ..core import live_defaults

#: Clamp band for ln|β₂| when building the μ₂_mid path (|β₂| = e^±14).
LOGABS_CLAMP: float = 14.0
from .pairwise import (
    EventGroup,
    collect_pair_events,
    finalize_event_groups,
    group_events,
    insert_event_groups,
    refine_mesh_for_multiple_crossings,
)


def logabs_clamped(roots: np.ndarray) -> np.ndarray:
    """``np.log(np.abs(roots))`` with ±∞ clamped to ``±LOGABS_CLAMP``."""
    return np.clip(np.log(np.abs(roots)), -LOGABS_CLAMP, LOGABS_CLAMP)


def ensure_mu2mid(zm: ZeroManager, **run_kwargs) -> "Mu2MidZM":
    """Return an analyzed Mu2MidZM, reusing *zm* if it is already analyzed.

    A plain ``ZeroManager`` cannot be analyzed in place, so a fresh
    ``Mu2MidZM`` is constructed from its ``(poly, E_ref, mu1)``, run and
    analyzed with ``tie_tol=core.CONTINUUM_TOL``.  This is the single bootstrap
    every consumer (winding, continuum extraction) goes through, kept next
    to the class it builds.
    """
    if isinstance(zm, Mu2MidZM) and getattr(zm, '_analyzed', False):
        if run_kwargs:
            warnings.warn(
                f"zm is already analyzed; run kwargs {run_kwargs} are "
                f"ignored (analyze() only runs once per Mu2MidZM)"
            )
        return zm
    m = Mu2MidZM(zm.poly, zm.E_ref, zm.mu1)
    m.run(**run_kwargs)
    m.analyze(tie_tol=core.CONTINUUM_TOL)
    return m


# Back-compat alias: the bootstrap historically lived in crossings.py.
_ensure_mu2mid = ensure_mu2mid


# ---------------------------------------------------------------------------
# ItemView
# ---------------------------------------------------------------------------

class ItemView(NamedTuple):
    """Per-segment representative-item view of the K β₂ tracks.

    Cluster semantics: a cluster is a set of columns whose moduli stay equal
    over the WHOLE segment (continuum).  Isolated crossings are never
    clustered here.
    """
    rep_cols: np.ndarray
    mults: np.ndarray
    item_logabs: np.ndarray
    item_tang_re: np.ndarray
    sort_to_item: np.ndarray
    j_lo: np.ndarray
    j_hi: np.ndarray


# ---------------------------------------------------------------------------
# Simplified μ₂_mid path
# ---------------------------------------------------------------------------

@dataclass
class Mu2MidPiece:
    """One smooth interval of the μ₂_mid path.

    *poly* is in numpy poly order and evaluated at ``x = theta1 - theta0``;
    *dpoly* is its cached derivative.
    """
    theta0: float
    theta1: float
    v0: float
    dv0: float
    v1: float
    dv1: float
    poly: np.ndarray
    dpoly: np.ndarray | None = None


def _make_piece(t0, t1, v0, dv0, v1, dv1) -> Mu2MidPiece:
    h = t1 - t0
    poly = hermite_interp_poly(h, v0, dv0, v1, dv1)
    return Mu2MidPiece(
        theta0=t0, theta1=t1, v0=v0, dv0=dv0, v1=v1, dv1=dv1,
        poly=poly, dpoly=np.polyder(poly),
    )


def _make_constant_piece(t0, t1, value) -> Mu2MidPiece:
    return Mu2MidPiece(
        theta0=t0, theta1=t1, v0=value, dv0=0.0, v1=value, dv1=0.0,
        poly=np.array([0.0, value]), dpoly=np.array([0.0]),
    )


class Mu2Mid:
    """Independent piecewise-smooth μ₂_mid path for loop-winding.

    NOT a ZeroManager subclass and holds no reference to one.
    """

    def __init__(self, pieces: list[Mu2MidPiece]):
        self.pieces = pieces
        self._ends = np.array([p.theta1 for p in pieces], dtype=float)

    @property
    def theta1(self) -> np.ndarray:
        return np.array(
            [self.pieces[0].theta0] + [p.theta1 for p in self.pieces],
            dtype=float,
        )

    @property
    def values(self) -> np.ndarray:
        return np.array(
            [self.pieces[0].v0] + [p.v1 for p in self.pieces],
            dtype=float,
        )

    @property
    def derivs(self) -> np.ndarray:
        return np.array(
            [self.pieces[0].dv0] + [p.dv1 for p in self.pieces],
            dtype=float,
        )

    @property
    def breakpoints(self) -> np.ndarray:
        """True derivative-discontinuity points (quad split points).

        Ordinary C1 knots are NOT returned: a quad segment may span several
        smooth pieces as long as value and derivative join continuously.
        """
        bps: list[float] = []
        for p, q in zip(self.pieces[:-1], self.pieces[1:]):
            same_value = np.isclose(p.v1, q.v0, rtol=1e-12, atol=1e-12)
            same_deriv = np.isclose(p.dv1, q.dv0, rtol=1e-12, atol=1e-12)
            if not (same_value and same_deriv):
                bps.append(float(q.theta0))
        return np.array(bps, dtype=float)

    def value_deriv(self, theta1: float) -> tuple[float, float]:
        if not self.pieces:
            return 0.0, 0.0
        t = float(theta1)
        if t <= self.pieces[0].theta0:
            p = self.pieces[0]
            return float(p.v0), 0.0
        if t >= self.pieces[-1].theta1:
            p = self.pieces[-1]
            return float(p.v1), 0.0
        idx = int(np.searchsorted(self._ends, t, side='right'))
        idx = max(0, min(idx, len(self.pieces) - 1))
        p = self.pieces[idx]
        x = t - p.theta0
        dpoly = p.dpoly if p.dpoly is not None else np.polyder(p.poly)
        return float(np.polyval(p.poly, x)), float(np.polyval(dpoly, x))

    def value(self, theta1: float) -> float:
        return self.value_deriv(theta1)[0]

    def values_at(self, theta1: np.ndarray) -> np.ndarray:
        """Vectorized μ₂_mid evaluation on a 1-D θ array."""
        t = np.asarray(theta1, dtype=float)
        if not self.pieces:
            return np.zeros_like(t)
        idx = np.searchsorted(self._ends, t, side='right')
        idx = np.clip(idx, 0, len(self.pieces) - 1)
        out = np.empty_like(t)
        for j, p in enumerate(self.pieces):
            mask = idx == j
            if not np.any(mask):
                continue
            x = t[mask] - p.theta0
            out[mask] = np.polyval(p.poly, x)
        return out


# ---------------------------------------------------------------------------
# Mu2MidZM
# ---------------------------------------------------------------------------

class Mu2MidZM(ZeroManager):
    """``ZeroManager`` + ItemView analysis + pairwise crossings + μ₂_mid path."""

    def __init__(self, poly, E_ref, mu1):
        super().__init__(poly, E_ref, mu1)
        self._analyzed = False
        self.has_continuum = False
        self._continuum_clusters: list = []
        self._item_views: list[ItemView] = []
        self._pair_events: list = []
        self._event_groups: list[EventGroup] = []
        self.mu2_mid: Mu2Mid | None = None

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    @live_defaults(tie_tol="core:CONTINUUM_TOL", crossing_tol="sgbz.pairwise:CROSSING_TOL",
                     min_direction_deriv="sgbz.pairwise:MIN_DIRECTION_DERIV",
                     refine_max_rounds="sgbz.pairwise:REFINE_MAX_ROUNDS",
                     refine_safety_factor="sgbz.pairwise:REFINE_SAFETY_FACTOR",
                     refine_max_subintervals="sgbz.pairwise:REFINE_MAX_SUBINTERVALS",
                     refine_max_total_inserts="sgbz.pairwise:REFINE_MAX_TOTAL_INSERTS")
    def analyze(
        self,
        continuum_clusters: list | None = None,
        *,
        tie_tol: Optional[float] = None,
        crossing_tol: Optional[float] = None,
        min_direction_deriv: Optional[float] = None,
        verbose: bool = False,
        refine_multi_crossings: bool = True,
        refine_max_rounds: Optional[int] = None,
        refine_safety_factor: Optional[float] = None,
        refine_max_subintervals: Optional[int] = None,
        refine_max_total_inserts: Optional[int] = None,
    ) -> None:
        """Run pairwise crossing analysis and build the μ₂_mid path.

        ``continuum_clusters``: per-segment column-tuples (the ItemView
        clustering source).  ``None`` runs the internal whole-segment vote.

        ``refine_multi_crossings`` enables the pre-crossing mesh refinement
        pass (:func:`pygbz2d.sgbz.pairwise.refine_mesh_for_multiple_crossings`)
        that isolates intervals with two or more close crossings before the
        sign-change scan runs.
        """
        if self._analyzed:
            return
        self._analyzed = True

        if continuum_clusters is None:
            continuum_clusters = self._detect_continuum_clusters_internal(tie_tol)
        self._continuum_clusters = continuum_clusters
        self._item_views = self._build_item_views(continuum_clusters)

        # 0. pre-crossing mesh refinement: isolate multiple crossings per
        # interval so the sign-change scan below sees each one.
        if refine_multi_crossings:
            n_refined = refine_mesh_for_multiple_crossings(
                self,
                tie_tol=tie_tol,
                crossing_tol=crossing_tol,
                max_rounds=refine_max_rounds,
                safety_factor=refine_safety_factor,
                max_subintervals=refine_max_subintervals,
                max_total_inserts=refine_max_total_inserts,
            )
            if verbose:
                print(f"[analyze] multi-crossing refinement inserted "
                      f"{n_refined} rows")
            self._continuum_clusters = self._detect_continuum_clusters_internal(
                tie_tol)
            self._item_views = self._build_item_views(self._continuum_clusters)

        # 1. representative-item pairwise intersections (zero dedup).
        events = collect_pair_events(
            self,
            crossing_tol=crossing_tol,
            min_direction_deriv=min_direction_deriv,
        )
        if verbose:
            print(f"[analyze] pairwise events: {len(events)}")

        # 2. merge close θ* (< crossing_tol) into one EventGroup per mesh row.
        groups = group_events(self, events, merge_tol=crossing_tol)
        insert_event_groups(self, groups)
        if verbose:
            print(f"[analyze] event groups: {len(groups)}")

        # 3. rebuild ItemView on the refined mesh, then finalize groups.
        self._continuum_clusters = self._detect_continuum_clusters_internal(tie_tol)
        self._item_views = self._build_item_views(self._continuum_clusters)
        finalize_event_groups(self, groups)

        self._pair_events = events
        self._event_groups = groups

        # 4. inline continuum detection (j_lo == j_hi ⟺ LineSubset).
        self.has_continuum = any(
            bool(np.any(view.j_lo == view.j_hi)) for view in self._item_views
        )

        # 5. loop-winding path.
        self.mu2_mid = build_mu2_mid(self, groups)
        self._sync_compat_arrays()

    @live_defaults(tie_tol="core:CONTINUUM_TOL", crossing_tol="sgbz.pairwise:CROSSING_TOL",
                     min_direction_deriv="sgbz.pairwise:MIN_DIRECTION_DERIV",
                     refine_max_rounds="sgbz.pairwise:REFINE_MAX_ROUNDS",
                     refine_safety_factor="sgbz.pairwise:REFINE_SAFETY_FACTOR",
                     refine_max_subintervals="sgbz.pairwise:REFINE_MAX_SUBINTERVALS",
                     refine_max_total_inserts="sgbz.pairwise:REFINE_MAX_TOTAL_INSERTS")
    def build_mu2_mid(
        self,
        continuum_clusters: list | None = None,
        *,
        tie_tol: Optional[float] = None,
        crossing_tol: Optional[float] = None,
        min_direction_deriv: Optional[float] = None,
        verbose: bool = False,
        refine_multi_crossings: bool = True,
        refine_max_rounds: Optional[int] = None,
        refine_safety_factor: Optional[float] = None,
        refine_max_subintervals: Optional[int] = None,
        refine_max_total_inserts: Optional[int] = None,
    ) -> None:
        """Legacy-compatible alias for :meth:`analyze`."""
        self.analyze(
            continuum_clusters=continuum_clusters,
            tie_tol=tie_tol,
            crossing_tol=crossing_tol,
            min_direction_deriv=min_direction_deriv,
            verbose=verbose,
            refine_multi_crossings=refine_multi_crossings,
            refine_max_rounds=refine_max_rounds,
            refine_safety_factor=refine_safety_factor,
            refine_max_subintervals=refine_max_subintervals,
            refine_max_total_inserts=refine_max_total_inserts,
        )

    # ------------------------------------------------------------------
    # ItemView construction (unchanged semantics)
    # ------------------------------------------------------------------

    def _detect_continuum_clusters_for_segment(
        self, seg, tie_tol: float,
    ) -> list:
        """Same-modulus column clusters of ONE segment.

        Per-segment kernel of :meth:`_detect_continuum_clusters_internal`;
        extracted so the multi-crossing refinement loop can rebuild a single
        segment's clusters after inserting rows there, without paying for
        the untouched segments again.
        """
        K = self.K
        N = len(seg.theta1_arr)
        if N == 0:
            return []
        logabs = np.log(np.abs(seg.tracked_roots))
        same = np.zeros((K, K), dtype=bool)
        for j in range(K):
            for k in range(j + 1, K):
                frac_in_band = np.mean(
                    np.abs(logabs[:, j] - logabs[:, k]) < tie_tol)
                if frac_in_band > core.CONTINUUM_FRAC:
                    same[j, k] = same[k, j] = True
        visited = [False] * K
        clusters: list[tuple] = []
        for j in range(K):
            if visited[j]:
                continue
            cluster: list[int] = []
            stack = [j]
            visited[j] = True
            while stack:
                c = stack.pop()
                cluster.append(c)
                for k in range(K):
                    if not visited[k] and same[c, k]:
                        visited[k] = True
                        stack.append(k)
            if len(cluster) > 1:
                clusters.append(tuple(cluster))
        return clusters

    def _detect_continuum_clusters_internal(self, tie_tol: float) -> list:
        """Per-segment same-modulus clusters of zero-curve COLUMNS."""
        return [
            self._detect_continuum_clusters_for_segment(seg, tie_tol)
            for seg in self.segments
        ]

    def _build_item_view_for_segment(
        self, s_idx: int, seg, clusters: list,
    ) -> ItemView:
        """ItemView of ONE segment from its (optional) continuum clusters.

        Per-segment kernel of :meth:`_build_item_views`; *clusters* must
        already be resolved for this segment (no index-based fallback here).
        """
        M = self.M
        K = self.K
        N = len(seg.theta1_arr)
        if N == 0:
            return ItemView(
                np.array([], dtype=int), np.array([], dtype=int),
                np.empty((0, 0)), np.empty((0, 0)),
                np.empty((0, K), dtype=int),
                np.array([], dtype=int), np.array([], dtype=int),
            )
        logabs = np.log(np.abs(seg.tracked_roots))
        tang_re = seg.tangents.real

        cluster_cols: set[int] = set()
        for c in clusters:
            cluster_cols.update(c)

        rep_cols: list[int] = []
        mults: list[int] = []
        for c in clusters:
            rep_cols.append(int(c[0]))
            mults.append(len(c))
        for j in range(K):
            if j not in cluster_cols:
                rep_cols.append(j)
                mults.append(1)
        rep_cols_arr = np.array(rep_cols, dtype=int)
        mults_arr = np.array(mults, dtype=int)

        item_logabs = logabs[:, rep_cols_arr]
        item_tang_re = tang_re[:, rep_cols_arr]

        sort_to_item = np.empty((N, K), dtype=int)
        order = np.argsort(item_logabs, axis=1)
        for i in range(N):
            sort_to_item[i] = np.repeat(order[i], mults_arr[order[i]])

        j_lo = sort_to_item[:, M - 1]
        j_hi = sort_to_item[:, M]

        # Padding roots legitimately live at the outer sort positions
        # (-∞ for β₂=0, +∞ for β₂=∞).  They may not, however, become the
        # finite μ₂_mid boundary pair: that would make the SGBZ modulus
        # itself 0 or ∞, for which the ±14 clamp is only an unphysical
        # finite stand-in.  Fail explicitly instead of building such a
        # path.
        rows = np.arange(N)
        lo_vals = item_logabs[rows, j_lo]
        hi_vals = item_logabs[rows, j_hi]
        bad_rows = np.flatnonzero(
            ~np.isfinite(lo_vals) | ~np.isfinite(hi_vals))
        if bad_rows.size:
            details = []
            for i in bad_rows[:3]:
                li, hi = int(j_lo[i]), int(j_hi[i])
                details.append(
                    f"row={int(i)} items=({li},{hi}) "
                    f"logabs=({float(lo_vals[i]):.6g},"
                    f"{float(hi_vals[i]):.6g})")
            more = f" (+{bad_rows.size - len(details)} more)" \
                if bad_rows.size > len(details) else ""
            raise ValueError(
                f"padding β₂ root occupies the M-1/M boundary pair in "
                f"segment {s_idx}: {'; '.join(details)}{more}; the "
                f"SGBZ boundary is degenerate at |β₂|=0 or |β₂|=∞ and "
                f"cannot be represented by the finite μ₂_mid path"
            )

        return ItemView(
            rep_cols_arr, mults_arr, item_logabs, item_tang_re,
            sort_to_item, j_lo, j_hi,
        )

    def _build_item_views(
        self, continuum_clusters: list | None,
    ) -> list[ItemView]:
        """Per-segment ItemView from the (optional) continuum clusters.

        Dispatches through the CLASS so unbound calls with a duck-typed
        ``zm`` (SimpleNamespace in tests) keep working — only ``self.M`` /
        ``self.K`` are read, never instance methods.
        """
        views: list[ItemView] = []
        for s_idx, seg in enumerate(self.segments):
            clusters = (continuum_clusters[s_idx]
                        if continuum_clusters and s_idx < len(continuum_clusters)
                        else [])
            views.append(Mu2MidZM._build_item_view_for_segment(
                self, s_idx, seg, clusters))
        return views

    # ------------------------------------------------------------------
    # Compatibility accessors (demos / older callers)
    # ------------------------------------------------------------------

    def _sync_compat_arrays(self) -> None:
        path = self.mu2_mid
        if path is None:
            return
        self.mu2_mid_theta1 = path.theta1
        self.mu2_mid_values = path.values
        self.mu2_mid_derivs = path.derivs

        jls, jhs = [], []
        seg_vals, seg_ders = [], []
        for s_idx, view in enumerate(self._item_views):
            seg = self.segments[s_idx]
            N = len(seg.theta1_arr)
            if N == 0:
                seg_vals.append(np.array([]))
                seg_ders.append(np.array([]))
                continue
            rows = np.arange(N)
            mu = (
                np.clip(view.item_logabs[rows, view.j_lo],
                        -LOGABS_CLAMP, LOGABS_CLAMP)
                + np.clip(view.item_logabs[rows, view.j_hi],
                          -LOGABS_CLAMP, LOGABS_CLAMP)
            ) / 2.0
            dm = (view.item_tang_re[rows, view.j_lo]
                  + view.item_tang_re[rows, view.j_hi]) / 2.0
            seg_vals.append(mu)
            seg_ders.append(dm)
            jls.append(view.j_lo)
            jhs.append(view.j_hi)
        self.mu2_mid_jlo = (np.concatenate(jls) if jls
                            else np.array([], dtype=int))
        self.mu2_mid_jhi = (np.concatenate(jhs) if jhs
                            else np.array([], dtype=int))
        self.seg_mu2_values = seg_vals
        self.seg_mu2_derivs = seg_ders


# ---------------------------------------------------------------------------
# μ₂_mid path construction
# ---------------------------------------------------------------------------

def _boundary_mean(
    view: ItemView,
    row: int,
    item_lo: int,
    item_hi: int,
) -> float:
    a = np.clip(float(view.item_logabs[row, item_lo]),
                -LOGABS_CLAMP, LOGABS_CLAMP)
    b = np.clip(float(view.item_logabs[row, item_hi]),
                -LOGABS_CLAMP, LOGABS_CLAMP)
    return float((a + b) / 2.0)


def _boundary_deriv(
    view: ItemView,
    row: int,
    item_lo: int,
    item_hi: int,
) -> float:
    def contrib(item: int) -> float:
        raw = float(view.item_logabs[row, item])
        if raw <= -LOGABS_CLAMP or raw >= LOGABS_CLAMP:
            return 0.0  # clipped value is saturated → derivative zero
        d = float(view.item_tang_re[row, item])
        return d if np.isfinite(d) else float('inf')

    return float((contrib(item_lo) + contrib(item_hi)) / 2.0)


def _event_rows_by_segment(groups: list[EventGroup]) -> dict:
    out: dict[int, dict[int, EventGroup]] = {}
    for g in groups:
        if g.row >= 0:
            out.setdefault(g.seg_idx, {})[g.row] = g
    return out


def _knots_for_segment(
    zm: Mu2MidZM,
    s_idx: int,
    row_groups: dict[int, EventGroup],
) -> list[tuple[float, float, float, float]]:
    """(theta, value, deriv_left, deriv_right) for every mesh row."""
    seg = zm.segments[s_idx]
    view = zm._item_views[s_idx]
    th = seg.theta1_arr
    n = len(th)
    knots: list[tuple[float, float, float, float]] = []

    for i in range(n):
        if i in row_groups:
            g = row_groups[i]
            lo = g.boundary_pair_left
            hi = g.boundary_pair_right
            if lo is None and hi is None:
                # fall back to the ItemView row (should not happen for
                # interior events; only defensively for a seam event).
                lo = hi = (int(view.j_lo[i]), int(view.j_hi[i]))
            if lo is not None:
                v_left = _boundary_mean(view, i, lo[0], lo[1])
                d_left = _boundary_deriv(view, i, lo[0], lo[1])
            else:
                v_left = 0.0
                d_left = float('inf')
            if hi is not None:
                v_right = _boundary_mean(view, i, hi[0], hi[1])
                d_right = _boundary_deriv(view, i, hi[0], hi[1])
            else:
                v_right = 0.0
                d_right = float('inf')
            if lo is None:
                v_left = v_right
                d_left = d_right
            if hi is None:
                v_right = v_left
                d_right = d_left
            value = float((v_left + v_right) / 2.0)
        else:
            lo = int(view.j_lo[i])
            hi = int(view.j_hi[i])
            value = _boundary_mean(view, i, lo, hi)
            d_left = d_right = _boundary_deriv(view, i, lo, hi)

        knots.append((float(th[i]), value, d_left, d_right))
    return knots


def _split_piece_at_clamp(
    t0: float, t1: float, v0: float, dv0: float, v1: float, dv1: float,
) -> list[Mu2MidPiece]:
    """Split one Hermite interval at its ±14 crossings."""
    h = t1 - t0
    poly = hermite_interp_poly(h, v0, dv0, v1, dv1)
    deriv = np.polyder(poly)

    # Amplitude guard before the two root solves: the Bernstein control
    # points bound the piece's whole value range (convex hull property),
    # so a piece whose hull stays inside the clamp band can never cross
    # ±14 and both np.roots calls would provably find nothing.  The band
    # only matters next to MRs where a boundary root's ln|β₂| diverges
    # (|β₂| = e^14 ≈ 1.2e6); ordinary pieces sit at O(1) and skip both
    # solves.  Divergent endpoint slopes degrade the piece to linear,
    # whose hull is just its two endpoint values.  Closed inequalities
    # are safe: a tangential interior touch of the bound is an output
    # no-op — splitting there and refitting from the split endpoint
    # values/derivatives reproduces the same cubic (cubic-Hermite
    # uniqueness), so missing such a root changes nothing.
    if np.isfinite(dv0) and np.isfinite(dv1):
        hull = (v0, v1, v0 + h * dv0 / 3.0, v1 - h * dv1 / 3.0)
    else:
        hull = (v0, v1)
    if max(hull) <= LOGABS_CLAMP and min(hull) >= -LOGABS_CLAMP:
        # Bit-identical to the unguarded no-crossing path below (same
        # _make_piece refit over the full interval).
        return [_make_piece(
            t0, t1,
            float(np.polyval(poly, 0.0)), float(np.polyval(deriv, 0.0)),
            float(np.polyval(poly, h)), float(np.polyval(deriv, h)),
        )]

    xs = [0.0, h]
    for bound in (-LOGABS_CLAMP, LOGABS_CLAMP):
        q = np.array(poly, dtype=complex)
        q[-1] -= bound
        for r in np.roots(q):
            if abs(r.imag) > 1e-12 * max(1.0, abs(r.real)):
                continue
            x = float(r.real)
            if 0.0 < x < h and min(abs(x - y) for y in xs) > 1e-14:
                xs.append(x)
    xs.sort()

    pieces: list[Mu2MidPiece] = []
    for a, b in zip(xs[:-1], xs[1:]):
        if b - a <= 0.0:
            continue
        mid = (a + b) / 2.0
        val_mid = float(np.polyval(poly, mid))
        if -LOGABS_CLAMP <= val_mid <= LOGABS_CLAMP:
            va = float(np.polyval(poly, a))
            vb = float(np.polyval(poly, b))
            da = float(np.polyval(deriv, a))
            db = float(np.polyval(deriv, b))
            pieces.append(_make_piece(t0 + a, t0 + b, va, da, vb, db))
        else:
            bound = (LOGABS_CLAMP if val_mid > LOGABS_CLAMP
                     else -LOGABS_CLAMP)
            pieces.append(_make_constant_piece(t0 + a, t0 + b, bound))
    return pieces


def build_mu2_mid(
    zm: Mu2MidZM,
    groups: list[EventGroup] | None = None,
) -> Mu2Mid:
    """Build the independent, piecewise-smooth μ₂_mid path.

    Requires ``zm._item_views`` to be rebuilt on the final (event-refined)
    mesh.  *groups* supplies event-row boundary-pair metadata; pass the same
    groups used by :meth:`Mu2MidZM.analyze`.
    """
    groups = groups or []
    row_groups = _event_rows_by_segment(groups)

    raw: list[tuple[float, float, float, float]] = []
    for s_idx in range(len(zm.segments)):
        raw.extend(_knots_for_segment(zm, s_idx, row_groups.get(s_idx, {})))

    if not raw:
        return Mu2Mid([])

    # Sort by θ and merge exact duplicates (shared MR rows, seam copies).
    raw.sort(key=lambda x: x[0])
    knots: list[tuple[float, float, float, float]] = []
    for t, v, dl, dr in raw:
        if knots and abs(t - knots[-1][0]) < 1e-14:
            t0, v0, dl0, dr0 = knots[-1]
            knots[-1] = (t0, (v0 + v) / 2.0, dl0, dr)
        else:
            knots.append((t, v, dl, dr))

    pieces: list[Mu2MidPiece] = []
    for k in range(len(knots) - 1):
        t0, v0, _, d0 = knots[k]
        t1, v1, d1, _ = knots[k + 1]
        pieces.extend(_split_piece_at_clamp(t0, t1, v0, d0, v1, d1))

    return Mu2Mid(pieces)
