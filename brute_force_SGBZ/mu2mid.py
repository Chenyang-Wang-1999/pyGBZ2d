'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-13
Copyright © Department of Physics, Tsinghua University. All rights reserved

μ₂_mid — the piecewise-smooth mid-boundary modulus as a first-class object.

``μ₂_mid(θ₁) = (ln|β_{j_lo}| + ln|β_{j_hi}|) / 2`` where ``j_lo``/``j_hi`` are
the item indices occupying sorted positions ``M-1``/``M`` (0-based).  It is
piecewise smooth; the breakpoints are the θ₁ where ``j_lo``/``j_hi`` change,
plus the multiple-root (MR) rows where the sort order is undefined.

This module is the production port of the prototype in
``demos/demo_zm_gbz.py`` (reviewed 2026-08-13, see
``log/2026-08-13-mu2mid-build.md``).  It implements the unified detector of
``log/2026-08-13-SGBZ算法梳理.md`` §1/§2:

  * **ItemView** — a per-segment representative-item view that collapses
    continuum (whole-segment same-modulus) column clusters into one item with
    a multiplicity, so ``sort_to_item`` is stable where ``abs_argsort`` is not.
    The no-continuum case is the special case ``n_items = K, mult = 1``.
  * **continuum inline** — after building, ``has_continuum`` is read directly
    off the ItemView (``j_lo == j_hi`` over a fraction of the segment), so the
    bisection gate is the build itself (§1/§6.3), with no separate two-point
    detector and its false-negative risk.
  * **cubic-Hermite bracketing** — sort-change crossings are refined by a
    bracketing iteration (not Newton), which stays finite and monotone near
    tangencies / MR branch points (§2.3).

Usage::

    zm = Mu2MidZM(poly, E_ref, mu1)
    zm.run()
    zm.build_mu2_mid()
    if zm.has_continuum: ...   # 1D subset (LineSubset)
    else: ...                  # 0D subsets via crossings.py
'''

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from continuation import ZeroManager


# ---------------------------------------------------------------------------
# Tunables
# ---------------------------------------------------------------------------

# Wall distance as a fraction of the local grid spacing.
_WALL_FRAC = 0.1
# Near-tie threshold (rows this close in modulus get walled so sort-change
# detection only runs on clean finite-tangent endpoints).
_TIE_TOL_DEFAULT = 1e-8
# cubic-Hermite bracketing convergence tolerance.
_MU2MID_XTOL_DEFAULT = 1e-12
_MU2MID_MAX_ITER = 50

# A genuine continuum is constant-modulus to ~1e-9; a transversal crossing
# leaves the level after one sample.  1e-6 separates them.
CONTINUUM_TOL = 1e-6
# Fraction of a segment's rows that must be in-band for a track to count as a
# continuum (real-analyticity ⇒ either degenerate everywhere or at isolated
# points; 0.9 separates them, with a warning on a near-threshold pass).
CONTINUUM_FRAC = 0.9

# Clamp band for ln|β₂|.  Roots at 0 or ∞ (e.g. the 0/∞ padding roots
# ``CharPoly.solve_roots_1d`` appends to a degree-deficient 1-D polynomial,
# or an interior root hitting 0/∞) would make ``np.log|β|`` = ±∞, and the
# boundary-pair mean μ₂_mid = (ln|β_{j_lo}| + ln|β_{j_hi}|)/2 would go ±∞,
# overflowing the winding loop (β₂ = exp(±∞+iθ₂)) and corrupting the
# track-vs-μ₂_mid comparison in crossing detection.  Clamp to ±14 — |β| ∈
# [e⁻¹⁴, e¹⁴] ≈ [1.2e-6, 1.2e6], aligned with arclength.ZERO_THRESHOLD /
# INF_THRESHOLD (1e-6 / 1e6).  A 0/∞ boundary root's ln|β| becomes the band
# edge, so μ₂_mid is the mean of the band edge and the other (finite) root's
# ln|β| — a finite, honest "the boundary ran to the band edge" value.
_LOGABS_CLAMP_L = 14.0


def logabs_clamped(roots: np.ndarray) -> np.ndarray:
    """``np.log(np.abs(roots))`` with ±∞ clamped to ``±_LOGABS_CLAMP_L``.

    All SGBZ μ₂_mid / crossing consumers must read ln|β₂| through this so a
    0/∞ root never propagates ±∞ into μ₂_mid (and from there into the winding
    loop or the track-vs-μ₂_mid sign checks).  See ``_LOGABS_CLAMP_L``.
    """
    return np.clip(np.log(np.abs(roots)), -_LOGABS_CLAMP_L, _LOGABS_CLAMP_L)


# ---------------------------------------------------------------------------
# Cubic Hermite helpers (shared with crossings.py)
# ---------------------------------------------------------------------------

def _cubic_hermite_coeffs(h, v0, dv0, v1, dv1):
    """cubic Hermite: f(0)=v0, f'(0)=dv0, f(h)=v1, f'(h)=dv1.

    Returns ``[a, b, c, d]`` for ``a·t³ + b·t² + c·t + d`` (descending
    powers, ready for ``np.roots``).
    """
    h2 = h * h
    h3 = h2 * h
    a = (2.0 * (v0 - v1)) / h3 + (dv0 + dv1) / h2
    b = (3.0 * (v1 - v0)) / h2 - (2.0 * dv0 + dv1) / h
    return np.array([a, b, dv0, v0])


def _cubic_roots_in_interval(coeffs, h):
    """ALL real roots of the cubic *coeffs* in ``[0, h)``.

    Left-closed, right-open: the right endpoint is owned by the adjacent
    interval, avoiding duplicate detections across sub-interval and segment
    boundaries.
    """
    def _valid(t):
        return 0.0 <= t < h

    scale = max(1.0, float(np.max(np.abs(coeffs))))
    a, b, c, d = coeffs

    # degenerate to quadratic / linear
    if abs(a) < 1e-15 * scale:
        if abs(b) < 1e-15 * scale:
            if abs(c) < 1e-15 * scale:
                return []
            t = -d / c
            return [float(t)] if _valid(t) else []
        disc = c * c - 4.0 * b * d
        if disc < 0:
            return []
        sqrt_disc = np.sqrt(disc)
        roots = []
        for t in ((-c + sqrt_disc) / (2 * b), (-c - sqrt_disc) / (2 * b)):
            if _valid(t):
                roots.append(float(t))
        return sorted(roots)

    # cubic: np.roots, keep all real roots in [0, h)
    raw = np.roots(coeffs)
    result = []
    for r in raw:
        if abs(r.imag) > 1e-10 * max(1.0, abs(r.real)):
            continue
        t = float(r.real)
        if _valid(t):
            result.append(t)
    return sorted(result)


def _dv_column(tang_row, col, va, vb, h, use_linear):
    """d(ln|β₂_col|)/dθ₁ at a mesh row, or linear fallback.

    The analytic tangent ``Re(V_col) = d(ln|β₂_col|)/dθ₁``; when the row
    touches an MR (tangent diverges) or tangents are absent, fall back to the
    secant slope ``(vb - va) / h`` — honest and bounded (§2.3).
    """
    if use_linear or tang_row is None:
        return (vb - va) / h
    d = float(tang_row[col].real)
    return (vb - va) / h if not np.isfinite(d) else d


# ---------------------------------------------------------------------------
# Data objects
# ---------------------------------------------------------------------------

class Mu2MidBreakpoint(NamedTuple):
    """One breakpoint of the piecewise-smooth μ₂_mid curve.

    μ₂_mid is continuous at a breakpoint (``value`` shared) but its derivative
    may jump (``deriv_left ≠ deriv_right``).  ``pair_kind`` names the sort pair
    that triggers the break; ``is_pmgbz`` is True iff ``pair_kind == 'M-1_M'``
    (the PMGBZ boundary crossing ``|β_M| = |β_{M+1}|``).  ``gap`` is the
    refined residual ``|ln|β_a| − ln|β_b||`` at the breakpoint (≈0 for a true
    crossing).
    """
    theta1: float
    value: float
    deriv_left: float
    deriv_right: float
    pair_kind: str            # 'M-2_M-1' | 'M-1_M' | 'M_M+1' | 'multi'
    columns: tuple            # (col_a, col_b) of the crossing pair
    is_pmgbz: bool
    gap: float                # refined residual (≈0 at a true crossing)


class ItemView(NamedTuple):
    """Per-segment representative-item view of the K β₂ tracks.

    Continuum clusters (whole-segment same-modulus columns) collapse to one
    item with multiplicity = cluster size; non-cluster columns are singleton
    items (mult 1).  ``sort_to_item[i, p]`` maps sorted position ``p`` to the
    item index occupying it (items sorted by modulus, repeated ``mult`` times),
    so it is stable where ``abs_argsort`` is not (cluster snapping makes
    ``abs_argsort`` jump arbitrarily between equal-modulus columns).

    ``j_lo = sort_to_item[:, M-1]``, ``j_hi = sort_to_item[:, M]``; ``j_lo ==
    j_hi`` ⟺ the M-1/M boundary pair is one continuum item (1D subset).  The
    no-continuum case is the special case ``n_items = K, mult = 1,
    sort_to_item == abs_argsort``.
    """
    rep_cols: np.ndarray      # (n_items,) representative column per item
    mults: np.ndarray         # (n_items,) multiplicity
    item_logabs: np.ndarray   # (N, n_items) representative ln|β₂|
    item_tang_re: np.ndarray  # (N, n_items) representative Re(V)
    sort_to_item: np.ndarray  # (N, K) sort position → item index
    j_lo: np.ndarray          # (N,) item index at sort pos M-1
    j_hi: np.ndarray          # (N,) item index at sort pos M


# ---------------------------------------------------------------------------
# Mu2MidZM
# ---------------------------------------------------------------------------

class Mu2MidZM(ZeroManager):
    """``ZeroManager`` + piecewise-smooth μ₂_mid construction.

    Usage: ``zm = Mu2MidZM(poly, E_ref, mu1); zm.run(); zm.build_mu2_mid()``.
    After ``build_mu2_mid`` the ``mu2_mid_*`` arrays, ``mu2_mid_breakpoints``,
    ``has_continuum`` and ``_item_views`` are available for crossing detection
    (§2) and winding (§6.4).  The build inserts walls and refined sort-change
    rows into the mesh (the intended side effect — breakpoint rows get exact
    roots + finite tangents).
    """

    def __init__(self, poly, E_ref, mu1):
        super().__init__(poly, E_ref, mu1)
        self._mu2_mid_built = False
        self.has_continuum = False

    # ------------------------------------------------------------------
    # Public build entry point
    # ------------------------------------------------------------------

    def build_mu2_mid(
        self,
        continuum_clusters: list | None = None,
        *,
        tie_tol: float = _TIE_TOL_DEFAULT,
        wall_frac: float = _WALL_FRAC,
        xtol: float = _MU2MID_XTOL_DEFAULT,
        max_iter: int = _MU2MID_MAX_ITER,
        verbose: bool = False,
    ) -> None:
        """Build the piecewise-smooth μ₂_mid representation.

        ``continuum_clusters``: per-segment list of column-tuples (the output
        of the continuum clusters).  ``None``
        (default) runs the cheap internal whole-segment same-modulus detection,
        which covers all three sort-adjacent pairs (M-2,M-1)/(M-1,M)/(M,M+1) —
        not just the boundary pair — because their continua also destabilise
        ``abs_argsort`` (§1).  Either source builds the ItemView with
        representatives; the simplified (``None``) and full versions are unified.

        Stages: 0 walls → 1 sort-change detection (outside walls) → 2
        cubic-Hermite bracketing refinement of each sort-change → 3 assemble
        flat arrays + breakpoints.  Inline continuum detection sets
        ``has_continuum``.
        """
        if continuum_clusters is None:
            continuum_clusters = self._detect_continuum_clusters_internal(tie_tol)
        self._continuum_clusters = continuum_clusters

        # Per-segment ItemView (no-continuum → column view; continuum → rep view)
        self._item_views = self._build_item_views(continuum_clusters)

        # ---- inline continuum detection (§1/§6.3) ----
        # j_lo == j_hi over a fraction of a segment ⟹ the M-1/M boundary pair
        # is one continuum item (mult ≥ 2) ⟹ 1D subset (|β_M|=|β_{M+1}| holds
        # identically).  Whole-segment frac criterion ⇒ no false negatives.
        self.has_continuum = self._detect_continuum_inline()

        # ---- Stage 0: walls ----
        walled_rows = self._find_near_tie_rows(tie_tol)
        if verbose:
            print(f"[mu2_mid] near-tie/MR rows: {len(walled_rows)}")
        self._build_walls(walled_rows, wall_frac)

        # walls inserted → rebuild views + re-identify wall rows
        self._item_views = self._build_item_views(continuum_clusters)
        walled_rows = self._find_near_tie_rows(tie_tol)

        # ---- Stage 1: sort-change detection (outside walls) ----
        events = self._detect_sort_changes(walled_rows)
        if verbose:
            print(f"[mu2_mid] sort-change events: {len(events)}")

        # ---- Stage 2: cubic-Hermite bracketing refinement ----
        refined = []  # (theta1, a, b, pair_kind, converged)
        for t_mid, a, b, t_lo, t_hi, pk, _s_idx, _i_idx in events:
            res = self._cubic_hermite_iterate(
                t_lo, t_hi, a, b, xtol, max_iter,
            )
            if res is None:
                if verbose:
                    print(f"[mu2_mid] no root in [{t_lo:.4e}, {t_hi:.4e}] "
                          f"pair {pk} cols ({a},{b})")
                continue
            theta_star, converged = res
            if not converged and verbose:
                print(f"[mu2_mid] bracket not converged at θ≈{theta_star:.6e} "
                      f"pair {pk}")
            refined.append((theta_star, a, b, pk, converged))

        self._refined_points = refined

        # ---- Stage 3: assemble (wall + refinement inserts are now in mesh) ----
        self._item_views = self._build_item_views(continuum_clusters)
        self._assemble()

        self._mu2_mid_built = True

    # ------------------------------------------------------------------
    # Mesh mutation hook — keep μ₂_mid in sync after a post-build insert
    # ------------------------------------------------------------------

    def insert_solution(self, theta1, seg_idx=None, i=None):
        """Insert a root row AND refresh the μ₂_mid representation.

        ``ZeroManager.insert_solution`` only updates ``SegmentData``
        (``theta1_arr`` / ``tracked_roots`` / ``tangents`` / ``abs_argsort``).
        The μ₂_mid representation — ``_item_views``, the per-segment
        ``seg_mu2_values`` / ``seg_mu2_derivs``, the flat ``mu2_mid_*``
        arrays, and ``mu2_mid_breakpoints`` — is *derived* in ``_assemble``
        and would otherwise go stale: a later read like
        ``crossings._bracket_crossing``'s ``seg_mu2_values[si][insert_at]``
        would index a row that no longer lines up with the grown mesh, and
        the winding's ``_Mu2MidPath`` would run on a mesh whose μ₂_mid curve
        ignores every inserted row.

        So after a **post-build** insertion (``_mu2_mid_built`` True — e.g.
        ``_bracket_crossing`` refining a multi-index crossing the build's
        pairwise bracket missed), rebuild the ItemView from the frozen
        ``_continuum_clusters`` (a transversal crossing does not change
        which columns are same-modulus, so the clusters are invariant) and
        re-``_assemble``.  Every later read then sees the inserted row —
        both the bracket's own ``g_pred`` / ``g'`` and the winding path.

        During ``build_mu2_mid`` itself (``_mu2_mid_built`` still False —
        walls in Stage 0, sort-change refinement in Stage 2) the refresh is
        skipped: those inserts mutate the mesh *before* Stage 3 ``_assemble``
        rebuilds μ₂_mid from the final mesh anyway, and syncing here would
        rebuild arrays that do not yet exist.  The cost of a post-build
        refresh is one ``_assemble`` (pure numpy bookkeeping on the ~10²–10³
        row mesh, µs) — negligible next to the ``_solve`` polynomial
        root-find that ``insert_solution`` already pays.
        """
        insert_at, changed = super().insert_solution(theta1, seg_idx, i)
        if changed and self._mu2_mid_built:
            self._item_views = self._build_item_views(self._continuum_clusters)
            self._assemble()
        return insert_at, changed

    # ------------------------------------------------------------------
    # ItemView construction
    # ------------------------------------------------------------------

    def _detect_continuum_clusters_internal(self, tie_tol: float) -> list:
        """Per-segment whole-segment same-modulus clusters (cheap, no voting).

        Covers **all** same-modulus column pairs (including the three
        sort-adjacent pairs (M-2,M-1)/(M-1,M)/(M,M+1)), not just the boundary
        pair — their continua also make ``abs_argsort`` jump.  Criterion: a
        pair (j,k) with ``max over rows |ln|β_j|−ln|β_k|| < tie_tol`` is
        same-modulus; transitive closure (BFS) merges them.  Whole-segment
        same-modulus = continuum (real-analyticity); an accidental isolated
        touch fails the whole-segment max and is auto-excluded.
        """
        clusters_per_seg: list[list] = []
        for seg in self.segments:
            K = self.K
            N = len(seg.theta1_arr)
            if N == 0:
                clusters_per_seg.append([])
                continue
            logabs = logabs_clamped(seg.tracked_roots)  # (N, K)
            same = np.zeros((K, K), dtype=bool)
            for j in range(K):
                for k in range(j + 1, K):
                    if np.max(np.abs(logabs[:, j] - logabs[:, k])) < tie_tol:
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
            clusters_per_seg.append(clusters)
        return clusters_per_seg

    def _build_item_views(
        self, continuum_clusters: list | None,
    ) -> list[ItemView]:
        """Per-segment ItemView from the (optional) continuum clusters."""
        M = self.M
        K = self.K
        views: list[ItemView] = []
        for s_idx, seg in enumerate(self.segments):
            N = len(seg.theta1_arr)
            if N == 0:
                views.append(ItemView(
                    np.array([], dtype=int), np.array([], dtype=int),
                    np.empty((0, 0)), np.empty((0, 0)),
                    np.empty((0, K), dtype=int),
                    np.array([], dtype=int), np.array([], dtype=int),
                ))
                continue
            logabs = logabs_clamped(seg.tracked_roots)  # (N, K)
            if seg.tangents is not None:
                tang_re = seg.tangents.real  # (N, K)
            else:
                tang_re = np.full((N, K), np.nan)

            clusters = (continuum_clusters[s_idx]
                        if continuum_clusters and s_idx < len(continuum_clusters)
                        else [])
            cluster_cols: set[int] = set()
            for c in clusters:
                cluster_cols.update(c)

            # items: cluster representative + non-cluster singletons
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

            item_logabs = logabs[:, rep_cols_arr]       # (N, n_items)
            item_tang_re = tang_re[:, rep_cols_arr]     # (N, n_items)

            # per row: items by ascending modulus, each repeated mult times
            sort_to_item = np.empty((N, K), dtype=int)
            order = np.argsort(item_logabs, axis=1)     # (N, n_items)
            for i in range(N):
                sort_to_item[i] = np.repeat(order[i], mults_arr[order[i]])

            j_lo = sort_to_item[:, M - 1]
            j_hi = sort_to_item[:, M]

            views.append(ItemView(
                rep_cols_arr, mults_arr, item_logabs, item_tang_re,
                sort_to_item, j_lo, j_hi,
            ))
        return views

    def _detect_continuum_inline(self) -> bool:
        """True iff any segment has a 1D continuum (j_lo == j_hi, mult ≥ 2).

        A boundary pair collapsing to one item over a fraction of the segment
        means ``|β_M| = |β_{M+1}|`` holds identically there — a 1D subset, where
        the average winding is undefined (§1/§6.3).

        Whole-segment same-modulus clustering (``_detect_continuum_clusters_internal``
        uses ``max over rows < tie_tol``) guarantees that a clustered pair is
        genuinely same-modulus everywhere, so a ``j_lo == j_hi`` (mult ≥ 2) row is
        a real continuum region, never a transversal 0D PMGBZ point (those are
        same-modulus at one θ only and fail the whole-segment max).  The frac
        threshold is therefore only a guard against a single spurious row; it is
        kept low (0.1) because a missed continuum would make the bisection
        mis-apply the 0D winding and run off the true boundary (§1: no false
        negatives).
        """
        for view in self._item_views:
            if len(view.j_lo) == 0:
                continue
            one_item = view.j_lo == view.j_hi
            if not np.any(one_item):
                continue
            one_item_idx = view.j_lo[one_item]
            if np.any(view.mults[one_item_idx] >= 2):
                if float(np.mean(one_item)) > 0.1:
                    return True
        return False

    # ------------------------------------------------------------------
    # Stage 0: walls
    # ------------------------------------------------------------------

    def _find_near_tie_rows(self, tie_tol: float) -> set:
        """Rows that need walls: MR rows + rows near-tie in any of the three
        sort-adjacent pairs (M-2,M-1)/(M-1,M)/(M,M+1).

        Same-item pairs (continuum, ``j_lo == j_hi``) are skipped — that is
        continuous-modulus equality, not a sort-change to detect.
        """
        walled: set[tuple[int, int]] = set()
        M = self.M
        K = self.K
        for s_idx, view in enumerate(self._item_views):
            seg = self.segments[s_idx]
            N = len(seg.theta1_arr)
            if N == 0:
                continue
            for i in range(N):
                is_mr_row = (
                    (i == 0 and seg.left_mr >= 0)
                    or (i == N - 1 and seg.right_mr >= 0)
                )
                near_tie = False
                for p in (M - 2, M - 1, M):
                    if p < 0 or p + 1 >= K:
                        continue
                    ia = int(view.sort_to_item[i, p])
                    ib = int(view.sort_to_item[i, p + 1])
                    if ia == ib:
                        continue  # continuum, skip
                    if abs(float(view.item_logabs[i, ia])
                           - float(view.item_logabs[i, ib])) < tie_tol:
                        near_tie = True
                        break
                if is_mr_row or near_tie:
                    walled.add((s_idx, i))
        return walled

    def _build_walls(self, walled_rows: set, wall_frac: float) -> None:
        """Insert a point at ``wall_frac·grid`` on each side of every wall row.

        Snapshot all wall rows' θ (mesh unchanged), then insert in descending θ
        order — descending keeps already-processed large θ free of later small-θ
        insertions, and ``insert_solution`` re-locates internally so row-index
        drift is harmless.  Boundary rows (segment endpoints) get a single side.
        """
        if not walled_rows:
            return
        snaps: list[tuple[float, float | None, float | None]] = []
        for s_idx, i in sorted(walled_rows):
            seg = self.segments[s_idx]
            th = seg.theta1_arr
            theta_row = float(th[i])
            tl = float(th[i - 1]) if i > 0 else None
            tr = float(th[i + 1]) if i < len(th) - 1 else None
            snaps.append((theta_row, tl, tr))
        wall_thetas: list[float] = []
        for theta_row, tl, tr in snaps:
            if tl is not None:
                wall_thetas.append(theta_row - wall_frac * (theta_row - tl))
            if tr is not None:
                wall_thetas.append(theta_row + wall_frac * (tr - theta_row))
        wall_thetas.sort(reverse=True)
        for wt in wall_thetas:
            try:
                self.insert_solution(wt)
            except (ValueError, RuntimeError):
                # out-of-range or degenerate interval — skip this wall point
                pass

    # ------------------------------------------------------------------
    # Stage 1: sort-change detection (outside walls)
    # ------------------------------------------------------------------

    def _detect_sort_changes(self, walled_rows: set) -> list:
        """Scan adjacent rows for j_lo/j_hi (item index) changes → events.

        Each event = ``(t_mid, a, b, t_lo, t_hi, pair_kind, seg_idx, row_idx)``
        where ``a, b`` are representative columns (``tracked_roots`` indices)
        for the cubic-Hermite step.  ``pair_kind`` from the changing item's last
        sort position in the left row (counting multiplicity): ``= M-2/M-1/M``
        ⟹ ``M-2_M-1 / M-1_M / M_M+1``; ``mult=1`` degenerates to ``min(inv_l)``.
        Dedup by ``(a, b, round(t_mid, 8))``.
        """
        M = self.M
        events: dict[tuple, tuple] = {}
        for s_idx, view in enumerate(self._item_views):
            seg = self.segments[s_idx]
            th = seg.theta1_arr
            N = len(th)
            if N < 2:
                continue
            sort_to_item = view.sort_to_item  # (N, K)
            rep_cols = view.rep_cols

            for i in range(N - 1):
                if (s_idx, i) in walled_rows or (s_idx, i + 1) in walled_rows:
                    continue
                changed = np.where(sort_to_item[i] != sort_to_item[i + 1])[0]
                for p in changed:
                    if p not in (M - 2, M - 1, M):
                        continue
                    ia = int(sort_to_item[i, p])
                    ib = int(sort_to_item[i + 1, p])
                    if ia == ib:
                        continue
                    # boundary = min(last sort pos of ia, last sort pos of ib)
                    pos_a = np.where(sort_to_item[i] == ia)[0]
                    pos_b = np.where(sort_to_item[i] == ib)[0]
                    if len(pos_a) == 0 or len(pos_b) == 0:
                        continue
                    boundary = min(int(pos_a[-1]), int(pos_b[-1]))
                    if boundary == M - 2:
                        pk = 'M-2_M-1'
                    elif boundary == M - 1:
                        pk = 'M-1_M'
                    elif boundary == M:
                        pk = 'M_M+1'
                    else:
                        continue
                    ra, rb = int(rep_cols[ia]), int(rep_cols[ib])
                    a, b = (ra, rb) if ra < rb else (rb, ra)
                    t_mid = float((th[i] + th[i + 1]) / 2)
                    key = (a, b, round(t_mid, 8))
                    if key not in events:
                        events[key] = (t_mid, a, b, float(th[i]),
                                       float(th[i + 1]), pk, s_idx, i)
        return list(events.values())

    # ------------------------------------------------------------------
    # Stage 2: cubic-Hermite bracketing iteration
    # ------------------------------------------------------------------

    def _cubic_hermite_iterate(
        self,
        theta_lo: float,
        theta_hi: float,
        a: int,
        b: int,
        xtol: float,
        max_iter: int,
    ) -> tuple[float, bool] | None:
        """Refine ``f = ln|β_a| − ln|β_b| = 0`` in bracket ``[theta_lo, theta_hi]``.

        Cubic-Hermite bracketing (not Newton): build the cubic from the two
        endpoints' (value, ``f' = Re(V_a) − Re(V_b)``) → ``np.roots`` predicts
        ``θ_pred`` → ``insert_solution`` takes the true ``f_pred`` and true
        derivative → **derivative sign (direction) + ``f_pred`` sign** fixes
        which side of the root ``θ_pred`` is on → tighten bracket → rebuild
        cubic.  Converges when the bracket width ``< xtol``.

        The derivative-sign criterion is what makes this robust where Newton and
        sign-change checks fail (§2.3): Newton's ``f/f'`` diverges near a
        branch point; sign-change fails when ``f_pred`` agrees with both ends
        (cubic error).  The derivative sign depends only on the direction, so
        the bracket still contains the root after every step.

        Returns ``(theta_star, converged)``; ``None`` if the interval has no
        transversal root (cubic false positive or tangent ``f'≈0``).

        Invariant: the bracket endpoints are always adjacent mesh rows.
        """
        for _ in range(max_iter):
            if theta_hi - theta_lo < xtol:
                return (theta_lo + theta_hi) / 2, True

            si, i = self.locate(theta_lo)
            seg = self.segments[si]
            th = seg.theta1_arr
            if i + 1 >= len(th) or abs(th[i] - theta_lo) > 1e-15 \
               or abs(th[i + 1] - theta_hi) > 1e-15:
                return None

            h = theta_hi - theta_lo
            logabs = logabs_clamped(seg.tracked_roots)
            v0 = float(logabs[i, a] - logabs[i, b])
            v1 = float(logabs[i + 1, a] - logabs[i + 1, b])

            # initial bracket must contain a root: transversal sort-change ⇒
            # opposite signs.  same sign ⇒ not a real crossing → abandon.
            if v0 * v1 > 0:
                return None

            if abs(v0) < xtol:
                return theta_lo, True
            if abs(v1) < xtol:
                return theta_hi, True

            N = len(th)
            touches_mr = (
                (i == 0 and seg.left_mr >= 0)
                or (i == N - 2 and seg.right_mr >= 0)
            )
            use_linear = touches_mr or seg.tangents is None
            tang_i = seg.tangents[i] if seg.tangents is not None else None
            tang_ip1 = seg.tangents[i + 1] if seg.tangents is not None else None
            va0, vb0 = float(logabs[i, a]), float(logabs[i, b])
            va1, vb1 = float(logabs[i + 1, a]), float(logabs[i + 1, b])
            dv0 = (_dv_column(tang_i, a, va0, va1, h, use_linear)
                   - _dv_column(tang_i, b, vb0, vb1, h, use_linear))
            dv1 = (_dv_column(tang_ip1, a, va0, va1, h, use_linear)
                   - _dv_column(tang_ip1, b, vb0, vb1, h, use_linear))

            coeffs = _cubic_hermite_coeffs(h, v0, dv0, v1, dv1)
            s_roots = _cubic_roots_in_interval(coeffs, h)
            if not s_roots:
                return None

            # pick the transversal (cubic f' ≠ 0) root nearest the midpoint;
            # a cubic f'≈0 root is a tangent touch, which a sort-change
            # (transversal crossing) should not produce — skip it.
            s_pred = None
            for s in s_roots:
                dv_cubic = (3.0 * coeffs[0] * s * s
                            + 2.0 * coeffs[1] * s + coeffs[2])
                if abs(dv_cubic) < 1e-15:
                    continue
                if s_pred is None or abs(s - h / 2) < abs(s_pred - h / 2):
                    s_pred = s
            if s_pred is None:
                return None

            theta_pred = theta_lo + s_pred
            if abs(theta_pred - theta_lo) < 1e-15 \
               or abs(theta_pred - theta_hi) < 1e-15:
                return None

            try:
                insert_at, _ = self.insert_solution(theta_pred, si, i)
            except (ValueError, RuntimeError):
                return None

            # true f_pred and true direction sign (from the inserted row's
            # tangent — more accurate than the cubic approximation).
            seg = self.segments[si]
            logabs = logabs_clamped(seg.tracked_roots)
            f_pred = float(logabs[insert_at, a] - logabs[insert_at, b])
            if seg.tangents is not None:
                f_prime = float(seg.tangents[insert_at, a].real
                                - seg.tangents[insert_at, b].real)
            else:
                f_prime = (3.0 * coeffs[0] * s_pred * s_pred
                           + 2.0 * coeffs[1] * s_pred + coeffs[2])

            if not np.isfinite(f_prime) or abs(f_prime) < 1e-15:
                return None  # direction undefined (branch point) — abandon

            # tighten bracket: derivative sign + f_pred sign fixes the side.
            increasing = f_prime > 0
            if (f_pred > 0) == increasing:
                theta_hi = theta_pred
            else:
                theta_lo = theta_pred

        return (theta_lo + theta_hi) / 2, True

    # ------------------------------------------------------------------
    # Stage 3: assemble
    # ------------------------------------------------------------------

    def _assemble(self) -> None:
        """Flatten the mesh into μ₂_mid arrays (values, derivs, j_lo, j_hi).

        Value/deriv use the ItemView's representative ``item_logabs`` /
        ``item_tang_re`` at the ``j_lo``/``j_hi`` item indices.  The per-segment
        arrays are also stored (``seg_mu2_values`` / ``seg_mu2_derivs``) so the
        post-build crossing detector treats μ₂_mid as a first-class curve —
        comparing each zero-curve (track) against *this* built curve rather than
        re-deriving a boundary-pair mean per interval (which is numerically a
        track-vs-track comparison, not track-vs-μ₂_mid, and goes inconsistent
        at the seam where the boundary pair swaps).
        """
        ths, vals, drvs, jls, jhs = [], [], [], [], []
        seg_row: list[tuple[int, int]] = []
        self.seg_mu2_values: list[np.ndarray] = []
        self.seg_mu2_derivs: list[np.ndarray] = []
        for s_idx, view in enumerate(self._item_views):
            seg = self.segments[s_idx]
            th = seg.theta1_arr
            N = len(th)
            if N == 0:
                self.seg_mu2_values.append(np.array([]))
                self.seg_mu2_derivs.append(np.array([]))
                continue
            rows = np.arange(N)
            j_lo = view.j_lo
            j_hi = view.j_hi
            mu = (view.item_logabs[rows, j_lo]
                  + view.item_logabs[rows, j_hi]) / 2.0
            dm = (view.item_tang_re[rows, j_lo]
                  + view.item_tang_re[rows, j_hi]) / 2.0
            ths.append(th)
            vals.append(mu)
            drvs.append(dm)
            jls.append(j_lo)
            jhs.append(j_hi)
            seg_row.extend([(s_idx, int(r)) for r in range(N)])
            self.seg_mu2_values.append(mu)
            self.seg_mu2_derivs.append(dm)

        self.mu2_mid_theta1 = np.concatenate(ths) if ths else np.array([])
        self.mu2_mid_values = np.concatenate(vals) if vals else np.array([])
        self.mu2_mid_derivs = np.concatenate(drvs) if drvs else np.array([])
        self.mu2_mid_jlo = np.concatenate(jls) if jls else np.array([], dtype=int)
        self.mu2_mid_jhi = np.concatenate(jhs) if jhs else np.array([], dtype=int)
        self._seg_row = seg_row

        self._find_breakpoints()

    def _is_mr_row(self, s_idx: int, row: int) -> bool:
        """Whether this row is a segment's MR boundary (cluster snapping makes
        ``abs_argsort`` unstable there — sort-change detection is skipped and
        MR breakpoints are taken from ``multiple_roots`` directly)."""
        seg = self.segments[s_idx]
        N = len(seg.theta1_arr)
        if row == 0 and seg.left_mr >= 0:
            return True
        if row == N - 1 and seg.right_mr >= 0:
            return True
        return False

    def _find_breakpoints(self) -> None:
        """Sort-change breakpoints = adjacent rows where (j_lo, j_hi) change.

        ``abs_argsort`` only **identifies** (which item changed, the crossing
        pair, ``pair_kind``); the **precise** value/gap/derivatives come from
        the column indices + the θ* row's tangent (§2.4):

          * ``theta1``/``value``/``gap`` taken at θ* (the smaller-gap row, i.e.
            the refinement point);
          * ``deriv_left`` = ``(V_{j_lo_l} + V_{j_hi_l})/2`` at θ* with the
            LEFT segment's item indices; ``deriv_right`` with the RIGHT's.
            A swap (same {a,b} set) ⇒ strictly continuous; an internal break
            (set changes) ⇒ real jump.

        MR rows are skipped here (handled by :meth:`_find_mr_breakpoints`).
        """
        M = self.M
        th = self.mu2_mid_theta1
        jl = self.mu2_mid_jlo
        jh = self.mu2_mid_jhi
        bps: list[Mu2MidBreakpoint] = []
        for k in range(len(th) - 1):
            if jl[k] == jl[k + 1] and jh[k] == jh[k + 1]:
                continue
            s_idx_l, r_l = self._seg_row[k]
            s_idx_r, r_r = self._seg_row[k + 1]
            if self._is_mr_row(s_idx_l, r_l) or self._is_mr_row(s_idx_r, r_r):
                continue
            view_l = self._item_views[s_idx_l]
            view_r = self._item_views[s_idx_r]
            set_l = {int(jl[k]), int(jh[k])}
            set_r = {int(jl[k + 1]), int(jh[k + 1])}
            symdiff = set_l ^ set_r
            if not symdiff:
                ca, cb = int(jl[k]), int(jh[k])
            elif len(symdiff) == 2:
                ca, cb = sorted(symdiff)
            else:
                ca, cb = -1, -1
            pk = 'multi'
            if ca >= 0:
                s2i_l = view_l.sort_to_item[r_l]
                pos_a = np.where(s2i_l == ca)[0]
                pos_b = np.where(s2i_l == cb)[0]
                if len(pos_a) > 0 and len(pos_b) > 0:
                    boundary = min(int(pos_a[-1]), int(pos_b[-1]))
                    if boundary == M - 2:
                        pk = 'M-2_M-1'
                    elif boundary == M - 1:
                        pk = 'M-1_M'
                    elif boundary == M:
                        pk = 'M_M+1'

            j_lo_l, j_hi_l = int(jl[k]), int(jh[k])
            j_lo_r, j_hi_r = int(jl[k + 1]), int(jh[k + 1])

            # θ* row = the smaller-gap (more refined) side
            if ca >= 0:
                la_l = float(view_l.item_logabs[r_l, ca]
                             - view_l.item_logabs[r_l, cb])
                la_r = float(view_r.item_logabs[r_r, ca]
                             - view_r.item_logabs[r_r, cb])
                if abs(la_l) <= abs(la_r):
                    s_star, r_star, view_star = s_idx_l, r_l, view_l
                else:
                    s_star, r_star, view_star = s_idx_r, r_r, view_r
            else:
                s_star, r_star, view_star = s_idx_l, r_l, view_l
            seg_star = self.segments[s_star]
            theta_bp = float(seg_star.theta1_arr[r_star])
            la_star = view_star.item_logabs[r_star]
            tang_star = view_star.item_tang_re[r_star]

            value_bp = float((la_star[j_lo_l] + la_star[j_hi_l]) / 2.0)
            gap = abs(float(la_star[ca] - la_star[cb])) if ca >= 0 else float('nan')

            def _deriv(j1: int, j2: int) -> float:
                v = (float(tang_star[j1]) + float(tang_star[j2])) / 2.0
                return v if np.isfinite(v) else float('inf')

            deriv_left = _deriv(j_lo_l, j_hi_l)
            deriv_right = _deriv(j_lo_r, j_hi_r)

            rep_cols = view_star.rep_cols
            cols_bp = (int(rep_cols[ca]), int(rep_cols[cb])) if ca >= 0 else (-1, -1)

            bps.append(Mu2MidBreakpoint(
                theta1=theta_bp,
                value=value_bp,
                deriv_left=deriv_left,
                deriv_right=deriv_right,
                pair_kind=pk,
                columns=cols_bp,
                is_pmgbz=(pk == 'M-1_M'),
                gap=gap,
            ))
        self.mu2_mid_breakpoints = bps

        self.mu2_mid_breakpoints.extend(self._find_mr_breakpoints())
        self.mu2_mid_breakpoints.sort(key=lambda b: b.theta1)

    def _find_mr_breakpoints(self) -> list:
        """MR breakpoints from ``multiple_roots`` directly.

        MR rows have cluster snapping (exact equal modulus) so ``abs_argsort``
        is fully unstable; sort-change detection is unreliable there.  An MR
        whose cluster covers sorted position M-1 or M touches the μ₂_mid
        boundary pair → a breakpoint.  The derivative is singular at the MR
        branch point (``dβ/dθ`` diverges ⇒ ``V → ∞``), so
        ``deriv_left = deriv_right = inf``; downstream interpolation
        special-cases MR.
        """
        M = self.M
        th = self.mu2_mid_theta1
        bps: list[Mu2MidBreakpoint] = []
        for mr in self.multiple_roots:
            idx = int(np.argmin(np.abs(th - mr.theta1)))
            if not np.isclose(th[idx], mr.theta1, atol=1e-6):
                continue
            value = float(self.mu2_mid_values[idx])
            deriv = float(self.mu2_mid_derivs[idx])
            if not np.isfinite(deriv):
                deriv = float('inf')

            for cluster in mr.cluster_indices:
                covers = set(int(c) for c in cluster)
                if (M - 1) not in covers and M not in covers:
                    continue
                if {M - 1, M} <= covers:
                    pk = 'M-1_M'; is_pmgbz = True
                elif {M - 2, M - 1} <= covers:
                    pk = 'M-2_M-1'; is_pmgbz = False
                elif {M, M + 1} <= covers:
                    pk = 'M_M+1'; is_pmgbz = False
                else:
                    pk = 'multi'; is_pmgbz = False
                bps.append(Mu2MidBreakpoint(
                    theta1=float(mr.theta1),
                    value=value,
                    deriv_left=deriv,
                    deriv_right=deriv,
                    pair_kind=pk,
                    columns=(-1, -1),
                    is_pmgbz=is_pmgbz,
                    gap=0.0,
                ))
        return bps
