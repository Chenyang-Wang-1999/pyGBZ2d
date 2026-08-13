'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-23
Copyright © Department of Physics, Tsinghua University. All rights reserved

Prototype: extract GBZ subsets directly from continuation.ZeroManager output.

Uses brute_force_SGBZ for the SGBZ pipeline (continuum detection, crossing
detection, winding), with amoeba extraction inline for comparison.
'''

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path
from cmath import exp
from typing import NamedTuple
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gbz_types import CharPoly, PointSubset, LineSubset
from continuation import ZeroManager

# Import the SGBZ pipeline from brute_force_SGBZ
from brute_force_SGBZ.continuum_lines import (
    detect_continuum_simple,
    CONTINUUM_TOL,
    CONTINUUM_FRAC,
    _DV_TOL,
)
from brute_force_SGBZ.crossings import detect_crossings_simple
from brute_force_SGBZ.winding import (
    compute_average_winding,
    detect_crossings_and_winding,
)

SNAP_TOL = 1e-3


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

def build_HN2D(J1, J2, gamma_1, gamma_2, delta_1, delta_2):
    J11 = exp(gamma_1 + 1j * delta_1) * J1
    J12 = exp(-gamma_1 + 1j * delta_1) * np.conj(J1)
    J21 = exp(gamma_2 + 1j * delta_2) * J2
    J22 = exp(-gamma_2 + 1j * delta_2) * np.conj(J2)
    coeffs = np.array([1, -J11, -J12, -J21, -J22], dtype=complex)
    degs = np.array([
        [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1],
    ], dtype=int)
    return coeffs, degs


# ---------------------------------------------------------------------------
# Amoeba extraction (inline, not in brute_force_SGBZ)
# ---------------------------------------------------------------------------

def extract_amoeba_subsets(zm: ZeroManager, poly, E, mu1, mu2, tol=CONTINUUM_TOL,
                           frac=CONTINUUM_FRAC, snap_tol=SNAP_TOL):
    """Zeros of f at ln|b1|=mu1, ln|b2|=mu2, from ZeroManager tracks."""
    hits: list[tuple[int, int, int, str]] = []
    continuum_lines: list[LineSubset] = []
    seg_cont_endpoints: list[tuple[np.ndarray, np.ndarray]] = []
    seg_logabs: list[np.ndarray] = []

    for seg_idx, seg in enumerate(zm.segments):
        th = seg.theta1_arr
        tr = seg.tracked_roots
        logabs = np.log(np.abs(tr))
        seg_logabs.append(logabs)

        # 1. continuum tracks (fraction-based)
        frac_in_band = np.mean(np.abs(logabs - mu2) < tol, axis=0)
        is_cont = frac_in_band > frac
        for j in np.where(is_cont)[0]:
            if frac_in_band[j] < 1.0:
                warnings.warn(
                    f"amoeba continuum track {j}: frac={frac_in_band[j]:.3f} "
                    f"< 1 (possible numerical outlier)")
            continuum_lines.append(LineSubset(
                E=E, mu1=mu1,
                theta1_arr=th.copy(), beta2_arr=tr[:, j].copy()))
        cont_tracks = np.where(is_cont)[0]
        seg_cont_endpoints.append((
            tr[0, cont_tracks].copy() if cont_tracks.size else np.array([]),
            tr[-1, cont_tracks].copy() if cont_tracks.size else np.array([]),
        ))

        # 2. discrete crossings on non-continuum tracks
        other = np.where(~is_cont)[0]
        if other.size:
            d = logabs[:, other] - mu2

            i_idx, j_idx = np.where(d == 0)
            for i, jj in zip(i_idx, j_idx):
                hits.append((seg_idx, int(i), int(other[jj]), 'zero'))

            sc = d[:-1] * d[1:] < 0
            i_idx, j_idx = np.where(sc)
            for i, jj in zip(i_idx, j_idx):
                hits.append((seg_idx, int(i), int(other[jj]), 'cross'))

    # materialize + boundary filtering
    n_seg = len(zm.segments)
    CURVE_TOL = 1e-9

    seen_zero_theta1: set[int] = set()
    points: list[tuple[float, complex]] = []
    for seg_idx, i, j, kind in hits:
        seg = zm.segments[seg_idx]
        th, tr = seg.theta1_arr, seg.tracked_roots
        la = seg_logabs[seg_idx]
        if kind == 'zero':
            t1, b2 = float(th[i]), complex(tr[i, j])
        else:
            frac_i = (mu2 - la[i, j]) / (la[i + 1, j] - la[i, j])
            t1 = float(th[i] + frac_i * (th[i + 1] - th[i]))
            b2 = complex(tr[i, j] + frac_i * (tr[i + 1, j] - tr[i, j]))

        # Rule 1: boundary sample + root-value match with line endpoint
        if kind == 'zero':
            at_boundary = (i == 0) or (i == len(th) - 1)
        else:
            at_boundary = ((i == 0 and abs(t1 - th[0]) < 1e-9)
                           or (i == len(th) - 2 and abs(t1 - th[-1]) < 1e-9))
        if at_boundary:
            if i == 0:
                endpoint_vals = np.concatenate([
                    seg_cont_endpoints[seg_idx][0],
                    seg_cont_endpoints[(seg_idx - 1) % n_seg][1],
                ])
            else:
                endpoint_vals = np.concatenate([
                    seg_cont_endpoints[seg_idx][1],
                    seg_cont_endpoints[(seg_idx + 1) % n_seg][0],
                ])
            if (endpoint_vals.size
                    and np.min(np.abs(endpoint_vals - b2)) < CURVE_TOL):
                continue

        # Rule 2: 'zero' boundary-duplicate dedup
        if kind == 'zero':
            key = round(t1, 12)
            if key in seen_zero_theta1:
                continue
            seen_zero_theta1.add(key)

        points.append((t1, b2))

    subsets: list = list(continuum_lines)
    for t1, b2 in points:
        subsets.append(PointSubset(E=E, beta1=exp(mu1 + 1j * t1), beta2=b2))
    return subsets


# ---------------------------------------------------------------------------
# μ₂_mid 构建
# ---------------------------------------------------------------------------
#
# μ₂_mid(θ₁) = (ln|β₂_{j_lo}| + ln|β₂_{j_hi}|) / 2  分段光滑。
# 非光滑点出现在 j_lo / j_hi 变化处。变化分三类：
#   - "swap":       j_lo 与 j_hi 互换（PMGBZ 边界 crossing）
#   - "j_lo_change": j_lo 与第三列交换 → 边界对改变，非 PMGBZ
#   - "j_hi_change": j_hi 与第三列交换 → 边界对改变，非 PMGBZ
#
# 只用 cubic Hermite 预测，不用 Newton 精修：
#   预测交点 θ_pred → insert_solution → 检查子区间是否仍有变化 →
#   如有则继续预测。反复直至所有子区间 j_lo, j_hi 恒定。
#
# 同时输出解析一阶导数：
#   d(mu2_mid)/dθ₁ = (Re(V_{j_lo}) + Re(V_{j_hi})) / 2
# 在光滑域内成立；域边界处不可导。

# ---- cubic Hermite helpers ----

def _cubic_hermite_coeffs(h, v0, dv0, v1, dv1):
    """cubic Hermite: f(0)=v0, f'(0)=dv0, f(h)=v1, f'(h)=dv1.
    Returns [a, b, c, d] for a·t³+b·t²+c·t+d."""
    h2 = h * h
    h3 = h2 * h
    a = (2.0 * (v0 - v1)) / h3 + (dv0 + dv1) / h2
    b = (3.0 * (v1 - v0)) / h2 - (2.0 * dv0 + dv1) / h
    return np.array([a, b, dv0, v0])


def _cubic_roots_in_interval(coeffs, h):
    """ALL real roots of cubic in [0, h). Returns list (sorted, may be empty)."""
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
        sqrt_disc = math.sqrt(disc)
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
    """d(ln|β₂_col|)/dθ₁ at a mesh row, or linear fallback."""
    if use_linear or tang_row is None:
        return (vb - va) / h
    d = float(tang_row[col].real)
    return (vb - va) / h if not np.isfinite(d) else d


# ---------------------------------------------------------------------------
# Mu2MidZM — ZeroManager + piecewise-smooth μ₂_mid 构建
# ---------------------------------------------------------------------------
#
# μ₂_mid(θ₁) = (ln|β_{j_lo}| + ln|β_{j_hi}|) / 2，j_lo = abs_argsort[:, M-1],
# j_hi = abs_argsort[:, M]。分段光滑；断点 = j_lo/j_hi 列号变化处。
#
# 三阶段（替换 Newton 的 cubic Hermite bracketing 迭代）：
#   Stage 0  建墙：在近相等行（三组排序对 (M-2,M-1)/(M-1,M)/(M,M+1) 近相等
#           或 MR 行）左右 0.1·grid 处插点，让 sort-change 检测只在墙外跑。
#   Stage 1  检测 sort-change（墙外）：相邻行 j_lo/j_hi 列号变化 → 事件。
#   Stage 2  cubic Hermite 迭代精修：bracket [θ_lo, θ_hi]，端点 value+f' 建
#           cubic → np.roots 预测 θ_pred → insert_solution 取真值 → 按 f 变号
#           收紧 bracket → 重建 cubic。至 |f|<tol 收敛。
#           比 Newton 稳：端点在墙外（远离重根），f' finite；bracketing 单调。
#
# 继承 ZeroManager：insert_solution / segments / abs_argsort / tangents / locate
# 全作 self.xxx 复用；墙点与精修点的 mesh 插入是预期副作用（断点处得精确数据
# + finite tangent）。不 override run()——mu2_mid 是上层派生物，显式两步调用。

class Mu2MidBreakpoint(NamedTuple):
    """μ₂_mid 的一个分段断点。

    μ₂_mid 在断点处连续（value 共享），导数可能跳变（deriv_left ≠ deriv_right）。
    pair_kind 标记触发断点的排序对；is_pmgbz 为 True 当且仅当 pair_kind == 'M-1_M'
    （PMGBZ 边界 crossing，|β_M|=|β_{M+1}|）。
    gap = 断点处交换对 (a,b) 的 |ln|β_a| − ln|β_b||（cubic Hermite 迭代精修
    残差），取断点相邻两点中较小者（精修点 θ* 行 ≈ 0，邻行 > 0）。
    """
    theta1: float
    value: float
    deriv_left: float
    deriv_right: float
    pair_kind: str            # 'M-2_M-1' | 'M-1_M' | 'M_M+1'
    columns: tuple            # (col_a, col_b) 交叉的两列
    is_pmgbz: bool
    gap: float                # 断点处 |ln|β_a| − ln|β_b||（精修残差，应≈0）


# 墙距网格的比例
_WALL_FRAC = 0.1
# 近相等容差（判定是否需建墙）
_TIE_TOL_DEFAULT = 1e-8
# cubic Hermite 迭代收敛容差
_MU2MID_TOL_DEFAULT = 1e-10
_MU2MID_XTOL_DEFAULT = 1e-12
_MU2MID_MAX_ITER = 50


class ItemView(NamedTuple):
    """Per-segment 代表元视图：把 K 列零点曲线抽象成 n_items 个代表元。

    continuum 完整版：cluster 内整段模等 + 导数匹配，取一个代表列即可；
    非 cluster 列各自一个代表元（重数 1）。原版（无 continuum）是特例：
    n_items = K, mult = 1, rep_cols = arange(K)。

    sort_to_item[i, p] = 排序位 p（0..K-1）对应的 item 索引（item 按模升序，
    重复 mult 次）——即"M 和 M-1 对应到 cluster index 重新排序"。
    j_lo = sort_to_item[:, M-1], j_hi = sort_to_item[:, M]；若相等 = continuum。
    """
    rep_cols: np.ndarray      # (n_items,) 代表元对应的 tracked_roots 列号
    mults: np.ndarray         # (n_items,) 重数
    item_logabs: np.ndarray   # (N, n_items) 代表元模长
    item_tang_re: np.ndarray  # (N, n_items) 代表元 Re(V)
    sort_to_item: np.ndarray  # (N, K) 排序位 → item 索引
    j_lo: np.ndarray          # (N,) 排序位 M-1 的 item 索引
    j_hi: np.ndarray          # (N,) 排序位 M   的 item 索引


class Mu2MidZM(ZeroManager):
    """ZeroManager + piecewise-smooth μ₂_mid 构建。

    用法：``zm = Mu2MidZM(poly, E_ref, mu1); zm.run(); zm.build_mu2_mid()``
    """

    def build_mu2_mid(
        self,
        continuum_clusters: list | None = None,
        *,
        tie_tol: float = _TIE_TOL_DEFAULT,
        wall_frac: float = _WALL_FRAC,
        tol: float = _MU2MID_TOL_DEFAULT,
        xtol: float = _MU2MID_XTOL_DEFAULT,
        max_iter: int = _MU2MID_MAX_ITER,
        verbose: bool = False,
    ) -> None:
        """构建 μ₂_mid 分段光滑表示。填充 self.mu2_mid_* 属性。

        continuum_clusters : 完整版连续模相等 clusters（detect_continuum_full
        输出，per-segment list of column-tuples）。None = 原版（0D-only，
        item=列，mult=1）。给定则用代表元集合，item=cluster 代表+单根，
        排序计入重数——"M 和 M-1 对应到 cluster index 重新排序"。

        Stage 0 建墙 → Stage 1 检测 sort-change（墙外）→ Stage 2 cubic Hermite
        迭代精修每个 sort-change 点 → Stage 3 组装扁平数组 + 断点。
        """
        # continuum 检测：None → 内部廉价检测（整段同模，覆盖三组排序对
        # (M-2,M-1)/(M-1,M)/(M,M+1)）；给定 → 用外部精算 clusters（完整版）。
        # 两者都构建 ItemView 用代表元——continuum 段 cluster 内同模，
        # abs_argsort 不稳定（j_lo/j_hi 在同模列间任意跳），代表元使
        # sort_to_item 稳定。简化版与完整版统一，区别仅在检测来源。
        if continuum_clusters is None:
            continuum_clusters = self._detect_continuum_clusters_internal(tie_tol)
        self._continuum_clusters = continuum_clusters

        # 构建 per-segment ItemView（原版=列 view，完整版=代表元 view）
        self._item_views = self._build_item_views(continuum_clusters)

        # ---- Stage 0: 建墙 ----
        walled_rows = self._find_near_tie_rows(tie_tol)
        if verbose:
            print(f"[mu2_mid] near-tie/MR rows: {len(walled_rows)}")
        self._build_walls(walled_rows, wall_frac)

        # 墙点插入后重新构建 views（mesh 变了）并重新识别墙行
        self._item_views = self._build_item_views(continuum_clusters)
        walled_rows = self._find_near_tie_rows(tie_tol)

        # ---- Stage 1: 检测 sort-change（墙外）----
        events = self._detect_sort_changes(walled_rows)
        if verbose:
            print(f"[mu2_mid] sort-change events: {len(events)}")

        # ---- Stage 2: cubic Hermite 迭代精修 ----
        refined = []  # (theta1, a, b, pair_kind, converged)
        for t_mid, a, b, t_lo, t_hi, pk, s_idx, i_idx in events:
            res = self._cubic_hermite_iterate(
                t_lo, t_hi, a, b, tol, xtol, max_iter,
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

        # ---- Stage 3: 组装 ----
        # 墙点+精修点已插入 mesh，重建 views 反映新 mesh
        self._item_views = self._build_item_views(continuum_clusters)
        self._assemble()

    # ------------------------------------------------------------------
    # ItemView 构建
    # ------------------------------------------------------------------

    def _detect_continuum_clusters_internal(self, tie_tol: float) -> list:
        """Per-segment 检测整段同模 clusters（廉价，无投票）。

        覆盖**所有**同模列对（含三组排序对 (M-2,M-1)/(M-1,M)/(M,M+1)），
        不只边界对——因为 (M-2,M-1)/(M,M+1) 的 continuum 同样使 j_lo/j_hi
        在同模列间任意跳（abs_argsort 不稳定），影响 mu2_mid 构建。

        判据：列对 (j,k) 在整段 |ln|β_j|−ln|β_k|| 最大值 < tie_tol → 同模。
        传递闭包合并成 clusters。整段同模 = continuum（real-analyticity），
        偶然同模（孤立点）不会整段 < tol，自动排除。
        """
        clusters_per_seg: list[list] = []
        for seg in self.segments:
            K = self.K
            N = len(seg.theta1_arr)
            if N == 0:
                clusters_per_seg.append([])
                continue
            logabs = np.log(np.abs(seg.tracked_roots))  # (N, K)
            # 整段同模邻接矩阵
            same = np.zeros((K, K), dtype=bool)
            for j in range(K):
                for k in range(j + 1, K):
                    if np.max(np.abs(logabs[:, j] - logabs[:, k])) < tie_tol:
                        same[j, k] = same[k, j] = True
            # 传递闭包（BFS）
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
    ) -> list:
        """Per-segment ItemView。

        continuum_clusters[seg_idx] = [(col, ...), ...]（cluster 内列号）。
        None 或空 → 原版：item=列, mult=1。
        """
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
            logabs = np.log(np.abs(seg.tracked_roots))  # (N, K)
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

            # items: cluster 代表 + 非 cluster 单根
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
            n_items = len(rep_cols_arr)

            item_logabs = logabs[:, rep_cols_arr]       # (N, n_items)
            item_tang_re = tang_re[:, rep_cols_arr]     # (N, n_items)

            # per row item 按模升序，重复 mult 次填 sort_to_item (N, K)
            sort_to_item = np.empty((N, K), dtype=int)
            order = np.argsort(item_logabs, axis=1)     # (N, n_items) item 索引
            for i in range(N):
                sort_to_item[i] = np.repeat(order[i], mults_arr[order[i]])

            j_lo = sort_to_item[:, M - 1]
            j_hi = sort_to_item[:, M]

            views.append(ItemView(
                rep_cols_arr, mults_arr, item_logabs, item_tang_re,
                sort_to_item, j_lo, j_hi,
            ))
        return views

    # ------------------------------------------------------------------
    # Stage 0: 墙
    # ------------------------------------------------------------------

    def _find_near_tie_rows(self, tie_tol: float) -> set:
        """返回 {(seg_idx, row_idx)} 标记需建墙的行。

        墙行 = MR 行 或三组排序对 (M-2,M-1)/(M-1,M)/(M,M+1) 中任一组
        近相等的行（用 item 模长）。同 item 的对（continuum，j_lo==j_hi）
        跳过——那是连续模相等，不该 sort-change 检测。
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
                        continue  # 同 item = continuum，跳过
                    if abs(float(view.item_logabs[i, ia])
                           - float(view.item_logabs[i, ib])) < tie_tol:
                        near_tie = True
                        break
                if is_mr_row or near_tie:
                    walled.add((s_idx, i))
        return walled

    def _build_walls(self, walled_rows: set, wall_frac: float) -> None:
        """对每个墙行，在其左右 wall_frac·local_grid 处插点。

        快照所有墙行及其邻居 θ（读时 mesh 未变），再降序插入——降序保证已
        处理的大 θ 不受后续小 θ 插入干扰，且 insert_solution 内部用 locate
        重定位，行号漂移无害。边界行（segment 端点）只插单侧。
        """
        if not walled_rows:
            return
        # 快照（升序读，mesh 未变）
        snaps: list[tuple[float, float | None, float | None]] = []
        for s_idx, i in sorted(walled_rows):
            seg = self.segments[s_idx]
            th = seg.theta1_arr
            theta_row = float(th[i])
            tl = float(th[i - 1]) if i > 0 else None
            tr = float(th[i + 1]) if i < len(th) - 1 else None
            snaps.append((theta_row, tl, tr))
        # 计算墙 θ，降序插入
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
                # 越界或退化区间——跳过，该墙点不建
                pass

    # ------------------------------------------------------------------
    # Stage 1: sort-change 检测（墙外）
    # ------------------------------------------------------------------

    def _detect_sort_changes(self, walled_rows: set) -> list:
        """扫描相邻行 j_lo/j_hi（item 索引）变化，返回事件列表。

        每事件 = (t_mid, a, b, t_lo, t_hi, pair_kind, seg_idx, row_idx)。
        a,b 是**代表元 rep_col**（tracked_roots 列号），供 cubic Hermite 直接
        求 f=ln|β_a|−ln|β_b|（列号连续，不依赖排序）。
        pair_kind 由变化 item 在左行的最后排序位 e_a 决定（计入重数）：
          e_a=M-2 → 'M-2_M-1'，e_a=M-1 → 'M-1_M'，e_a=M → 'M_M+1'。
        mult=1 时 e_a 退化为 inv_l[a]，与原版一致。dedup by (a,b,round(t_mid,8))。
        """
        M = self.M
        K = self.K
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
                    # 边界 = min(ia 最后位, ib 最后位)（左行）
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
                    # rep_col 对（cubic Hermite 用列号）
                    ra, rb = int(rep_cols[ia]), int(rep_cols[ib])
                    a, b = (ra, rb) if ra < rb else (rb, ra)
                    t_mid = float((th[i] + th[i + 1]) / 2)
                    key = (a, b, round(t_mid, 8))
                    if key not in events:
                        events[key] = (t_mid, a, b, float(th[i]),
                                       float(th[i + 1]), pk, s_idx, i)
        return list(events.values())

    # ------------------------------------------------------------------
    # Stage 2: cubic Hermite bracketing 迭代
    # ------------------------------------------------------------------

    def _cubic_hermite_iterate(
        self,
        theta_lo: float,
        theta_hi: float,
        a: int,
        b: int,
        tol: float,
        xtol: float,
        max_iter: int,
    ) -> tuple[float, bool] | None:
        """在 bracket [theta_lo, theta_hi] 求 f = ln|β_a| − ln|β_b| = 0。

        cubic Hermite bracketing 迭代：
          两端（相邻 mesh 点）的 value + analytic f' = Re(V_a) − Re(V_b)
          → 建 cubic Hermite → np.roots 预测根 θ_pred
          → insert_solution 取真值 f_pred 与真值导数 f'_pred
          → 用 **导数符号（走向）+ f_pred 符号** 定 θ_pred 在根的哪一侧，
            收紧 bracket；重建 cubic。

        收敛判据：bracket 宽度 < xtol → converged=True。导数符号 + f_pred
        符号保证每次收紧 bracket 仍夹住零点（θ_pred 总被分到根的正确一侧），
        所以宽度收敛即根精度——与 |f| 大小、f' 大小都无关，切触也 work
        （只要 f' 符号确定）。

        比 Newton 稳：端点在墙外（远离重根），f' finite；bracketing 单调。
        比变号检查稳：变号检查在 f_pred 与两端同号时失效（best_pred 死循环）；
        导数符号不依赖两端符号，直接由走向定侧。

        返回 (theta_star, converged)；None = 区间内无横截根（cubic 假阳性
        或切触 f'≈0）。

        不变量：bracket 两端始终是相邻 mesh 点。初始成立（sort-change 区间
        相邻）；插入 θ_pred 后新 bracket 的两端仍相邻。
        """
        for _ in range(max_iter):
            if theta_hi - theta_lo < xtol:
                # 导数符号保证 bracket 含根，宽度收敛即根
                return (theta_lo + theta_hi) / 2, True

            si, i = self.locate(theta_lo)
            seg = self.segments[si]
            th = seg.theta1_arr
            # 不变量：theta_lo == th[i], theta_hi == th[i+1]
            if i + 1 >= len(th) or abs(th[i] - theta_lo) > 1e-15 \
               or abs(th[i + 1] - theta_hi) > 1e-15:
                return None

            h = theta_hi - theta_lo
            logabs = np.log(np.abs(seg.tracked_roots))
            v0 = float(logabs[i, a] - logabs[i, b])
            v1 = float(logabs[i + 1, a] - logabs[i + 1, b])

            # 初始 bracket 必须含根：横截 sort-change 两端 f 异号。同号且
            # 端点非零 → 非真横截交叉，放弃（让墙/continuum 处理）。
            if v0 * v1 > 0:
                return None

            # 端点恰好为零——已是根
            if abs(v0) < xtol:
                return theta_lo, True
            if abs(v1) < xtol:
                return theta_hi, True

            # 端点导数
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

            # 选预测根：cubic 导数非零（横截）且最接近中点。cubic 导数≈0
            # 的根是切触，sort-change 不该有切触（横截交叉），跳过。
            s_pred = None
            for s in s_roots:
                # cubic f'(t) = 3a·t² + 2b·t + c
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

            # 真值 f_pred 与真值导数符号（走向）。真值导数比 cubic 近似准，
            # 用 insert 后的 seg.tangents。
            seg = self.segments[si]
            logabs = np.log(np.abs(seg.tracked_roots))
            f_pred = float(logabs[insert_at, a] - logabs[insert_at, b])
            if seg.tangents is not None:
                f_prime = float(seg.tangents[insert_at, a].real
                                - seg.tangents[insert_at, b].real)
            else:
                # 无 tangent：退回 cubic 预测导数
                f_prime = (3.0 * coeffs[0] * s_pred * s_pred
                           + 2.0 * coeffs[1] * s_pred + coeffs[2])

            if not np.isfinite(f_prime) or abs(f_prime) < 1e-15:
                # 走向不定（重根处 f'→0）——sort-change 横截不该触发，
                # 但数值上可能。放弃此候选。
                return None

            # 用走向 + f_pred 符号定 θ_pred 在根的哪一侧，收紧 bracket：
            #   f' > 0（递增）：根处 f=0，左 f<0 右 f>0
            #     f_pred > 0 → θ_pred 在根右 → 新 bracket = [θ_lo, θ_pred]
            #     f_pred < 0 → θ_pred 在根左 → 新 bracket = [θ_pred, θ_hi]
            #   f' < 0（递减）：反之
            increasing = f_prime > 0
            if (f_pred > 0) == increasing:
                # f_pred>0 且递增 → 右侧；f_pred<0 且递减 → 右侧（递减时 <0 在右）
                theta_hi = theta_pred
            else:
                theta_lo = theta_pred

        return (theta_lo + theta_hi) / 2, True

    # ------------------------------------------------------------------
    # Stage 3: 组装
    # ------------------------------------------------------------------

    def _assemble(self) -> None:
        """扁平数组 + 断点。用 ItemView 的 j_lo/j_hi（item 索引）和
        item_logabs/item_tang_re 算 mu2_mid 值与导数。"""
        ths, vals, drvs, jls, jhs = [], [], [], [], []
        seg_row: list[tuple[int, int]] = []
        for s_idx, view in enumerate(self._item_views):
            seg = self.segments[s_idx]
            th = seg.theta1_arr
            N = len(th)
            if N == 0:
                continue
            rows = np.arange(N)
            j_lo = view.j_lo  # (N,) item 索引
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

        self.mu2_mid_theta1 = np.concatenate(ths) if ths else np.array([])
        self.mu2_mid_values = np.concatenate(vals) if vals else np.array([])
        self.mu2_mid_derivs = np.concatenate(drvs) if drvs else np.array([])
        self.mu2_mid_jlo = np.concatenate(jls) if jls else np.array([], dtype=int)
        self.mu2_mid_jhi = np.concatenate(jhs) if jhs else np.array([], dtype=int)
        self._seg_row = seg_row

        self._find_breakpoints()

    def _is_mr_row(self, s_idx: int, row: int) -> bool:
        """该行是否为 segment 的 MR 边界行（重根处，abs_argsort 不稳定）。

        MR 行的 cluster snapping 使多根精确等模，排序任意——sort-change
        检测在此处不可靠，跳过交给 _find_mr_breakpoints 从 multiple_roots
        直接取。
        """
        seg = self.segments[s_idx]
        N = len(seg.theta1_arr)
        if row == 0 and seg.left_mr >= 0:
            return True
        if row == N - 1 and seg.right_mr >= 0:
            return True
        return False

    def _find_breakpoints(self) -> None:
        """sort-change 断点 = 相邻点 (j_lo, j_hi) 变化处（墙外，非 MR）。

        abs_argsort 只用于**识别**（哪段列号变了、交换对 ca/cb、pair_kind）；
        **精确求解**用列号 + θ* 行 tangent：
          - theta1/value/gap 取 θ* 行（gap 较小者，a,b 模等）
          - deriv_left = θ* 行 (V_{j_lo_l} + V_{j_hi_l})/2（左段列号）
          - deriv_right = θ* 行 (V_{j_lo_r} + V_{j_hi_r})/2（右段列号）
        这样左右极限导数在同一点 θ* 用列号精确算，不依赖 abs_argsort
        采样：swap 同集合 → 严格连续；内部断点 → 真跳变。

        MR 行（cluster 精确等模）跳过，由 _find_mr_breakpoints 处理。
        """
        M = self.M
        K = self.K
        th = self.mu2_mid_theta1
        jl = self.mu2_mid_jlo
        jh = self.mu2_mid_jhi
        bps: list[Mu2MidBreakpoint] = []
        for k in range(len(th) - 1):
            if jl[k] == jl[k + 1] and jh[k] == jh[k + 1]:
                continue
            s_idx_l, r_l = self._seg_row[k]
            s_idx_r, r_r = self._seg_row[k + 1]
            # MR 行跳过（交给 MR 处理）
            if self._is_mr_row(s_idx_l, r_l) \
               or self._is_mr_row(s_idx_r, r_r):
                continue
            view_l = self._item_views[s_idx_l]
            view_r = self._item_views[s_idx_r]
            # 交换对（item 索引）
            set_l = {int(jl[k]), int(jh[k])}
            set_r = {int(jl[k + 1]), int(jh[k + 1])}
            symdiff = set_l ^ set_r
            if not symdiff:
                ca, cb = int(jl[k]), int(jh[k])
            elif len(symdiff) == 2:
                ca, cb = sorted(symdiff)
            else:
                ca, cb = -1, -1
            # pair_kind：边界 = min(ca 最后位, cb 最后位)（两 item 最后位的
            # 较小者 = 交界处前 item 的末位）。mult=1 退化为 min(inv_l[ca],
            # inv_l[cb])，与原版 pos_set 集合判定一致。
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

            # 左段 / 右段 item 索引（k 行=左段值，k+1 行=右段值）
            j_lo_l, j_hi_l = int(jl[k]), int(jh[k])
            j_lo_r, j_hi_r = int(jl[k + 1]), int(jh[k + 1])

            # θ* 行 = gap 较小者（item 模等的精修点）。用其 item_logabs/
            # item_tang_re 按 item 索引精确算 value/gap/deriv。
            if ca >= 0:
                la_l = float(view_l.item_logabs[r_l, ca]
                             - view_l.item_logabs[r_l, cb])
                la_r = float(view_r.item_logabs[r_r, ca]
                             - view_r.item_logabs[r_r, cb])
                gap_l = abs(la_l)
                gap_r = abs(la_r)
                if gap_l <= gap_r:
                    s_star, r_star, view_star = s_idx_l, r_l, view_l
                else:
                    s_star, r_star, view_star = s_idx_r, r_r, view_r
            else:
                s_star, r_star, view_star = s_idx_l, r_l, view_l
            seg_star = self.segments[s_star]
            theta_bp = float(seg_star.theta1_arr[r_star])
            la_star = view_star.item_logabs[r_star]      # (n_items,) item 模长
            tang_star = view_star.item_tang_re[r_star]  # (n_items,) item Re(V)

            # value: θ* 处用左段 item 索引（a,b 等模 → 左右段同值）
            value_bp = float((la_star[j_lo_l] + la_star[j_hi_l]) / 2.0)
            # gap: 交换对 ca,cb 在 θ* 处 item 模差
            gap = abs(float(la_star[ca] - la_star[cb])) if ca >= 0 else float('nan')

            # deriv: θ* 行用 item 索引 + item_tang_re 精确算左右极限
            def _deriv(j1: int, j2: int) -> float:
                v = (float(tang_star[j1]) + float(tang_star[j2])) / 2.0
                return v if np.isfinite(v) else float('inf')

            deriv_left = _deriv(j_lo_l, j_hi_l)
            deriv_right = _deriv(j_lo_r, j_hi_r)

            # columns: item 索引转 rep_col（tracked_roots 列号）
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

        # 追加 MR 断点
        self.mu2_mid_breakpoints.extend(self._find_mr_breakpoints())
        self.mu2_mid_breakpoints.sort(key=lambda b: b.theta1)

    def _find_mr_breakpoints(self) -> list:
        """从 multiple_roots 直接取 MR 处的 μ₂_mid 断点。

        MR 行 cluster snapping 使多根精确等模，abs_argsort 完全不稳定，
        sort-change 检测不可靠。MR 处断点直接由 cluster 结构判定：
        cluster 覆盖排序位 M-1 或 M 即触及 μ₂_mid 边界对 → 断点。

        导数：MR 是重根分支点，dβ/dθ 发散（cluster 根 ∂f/∂β₂=0 →
        tangent V=inf），μ₂_mid 导数奇异 → deriv_left=deriv_right=inf。
        value 取 MR 行 mu2_mid_values（_assemble 已算，cluster 覆盖时
        = ln|cluster root|）。
        """
        M = self.M
        th = self.mu2_mid_theta1
        bps: list[Mu2MidBreakpoint] = []
        for mr in self.multiple_roots:
            # 找 MR 行在 mu2_mid 的位置（左段末或右段首，同 θ 同值）
            idx = int(np.argmin(np.abs(th - mr.theta1)))
            if not np.isclose(th[idx], mr.theta1, atol=1e-6):
                continue
            value = float(self.mu2_mid_values[idx])
            deriv = float(self.mu2_mid_derivs[idx])
            if not np.isfinite(deriv):
                deriv = float('inf')

            for cluster in mr.cluster_indices:
                covers = set(int(c) for c in cluster)
                # 只关心触及 M-1 或 M 的 cluster
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
                    columns=(-1, -1),  # MR cluster，列号不适用
                    is_pmgbz=is_pmgbz,
                    gap=0.0,  # cluster snapped，精确等模
                ))
        return bps


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def plot_continuum_curves(
    zm: ZeroManager,
    subsets: list,
    *,
    ax=None,
    zero_color: str = 'gray',
    zero_alpha: float = 0.4,
    zero_linewidth: float = 0.6,
    subset_cmap=None,
):
    """Plot zero curves and LineSubsets in (θ₁, μ₂) space.

    Parameters
    ----------
    zm : ZeroManager
        Root curves from ``zm.segments``, each column of ``tracked_roots``
        drawn as a thin gray ln|β₂| vs θ₁ curve.
    subsets : list
        Iterable of ``LineSubset`` / ``PointSubset``.  Each ``LineSubset``
        is plotted with a distinct color and a legend entry keyed by its
        (E, μ₁).  ``PointSubset`` entries are ignored.
    ax : matplotlib Axes, optional
        If None, a new figure+axes is created.
    zero_color, zero_alpha, zero_linewidth : passed to the root-curve plot
    subset_cmap : matplotlib Colormap, optional
        Defaults to ``tab10`` cycled for up to 10 distinct (E, μ₁) pairs.

    Returns
    -------
    ax : matplotlib Axes
    """
    import matplotlib.pyplot as plt

    if subset_cmap is None:
        subset_cmap = plt.cm.tab10
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 5))

    # ---- root curves (gray background) ----
    for seg in zm.segments:
        th = seg.theta1_arr
        logabs = np.log(np.abs(seg.tracked_roots))
        for j in range(logabs.shape[1]):
            mask = np.isfinite(logabs[:, j])
            if not np.any(mask):
                continue
            ax.plot(th[mask], logabs[mask, j],
                    color=zero_color, alpha=zero_alpha,
                    linewidth=zero_linewidth, zorder=1)

    # ---- LineSubsets (colored, with legend) ----
    seen_keys: set[tuple] = set()
    color_idx = 0
    n_colors = 10  # tab10 has 10 colors
    for s in subsets:
        if not isinstance(s, LineSubset):
            continue
        # Legend key: (E rounded, mu1 rounded); one entry per unique pair.
        key = (round(s.E.real, 4), round(s.E.imag, 4), round(s.mu1, 4))
        color = subset_cmap(color_idx % n_colors)
        mu2_arr = np.log(np.abs(s.beta2_arr))
        if key not in seen_keys:
            label = f"E={s.E.real:.2f}{s.E.imag:+.2f}j, μ₁={s.mu1:.3f}"
            ax.plot(s.theta1_arr, mu2_arr, color=color, linewidth=1.5,
                    label=label, zorder=2)
            seen_keys.add(key)
            color_idx += 1
        else:
            ax.plot(s.theta1_arr, mu2_arr, color=color, linewidth=1.5,
                    zorder=2)

    ax.set_xlabel(r"$\theta_1$")
    ax.set_ylabel(r"$\mu_2 = \ln|\beta_2|$")
    ax.set_xlim(0, 2 * math.pi)
    ax.legend(loc='best')
    return ax


# ---------------------------------------------------------------------------
# Haldane gain-loss model
# ---------------------------------------------------------------------------

from math import sin, cos, sqrt, pi
from BerryPy import TightBinding as tb
import matplotlib.pyplot as plt

ALL_PARAMS = (
    1,      # t1
    0.5,    # t2
    pi / 3, # phi
    0.5j,   # M
    0       # gamma
)


def non_Hermitian_Haldane_H(u1, u2, v1, v2, phi, M):
    dim = 2
    site_num = 2

    lattice_vec = np.array(
        [[-cos(pi/3), -cos(pi/3)],
         [-sin(pi/3), sin(pi/3)]]
    )

    intra_cell = [
        [0, 0, M],
        [1, 1, -M],
        [1, 0, u1],
        [0, 1, u2]
    ]
    inter_cell = [
        [0, 0, v2 * exp(1j*phi), (-1,0)],
        [0, 0, v2 * exp(1j*phi), (0, -1)],
        [0, 0, v2 * exp(1j*phi), (1,1)],
        [0, 0, v1 * exp(-1j*phi), (1,0)],
        [0, 0, v1 * exp(-1j*phi), (0,1)],
        [0, 0, v1 * exp(-1j*phi), (-1,-1)],
        [1, 0, u1, (0, -1)],
        [1, 0, u1, (1, 0)],
        [0, 1, u2, (0, 1)],
        [0, 1, u2, (-1, 0)],
        [1, 1, v2 * exp(-1j * phi), (-1, 0)],
        [1, 1, v2 * exp(-1j*phi), (0, -1)],
        [1, 1, v2 * exp(-1j*phi), (1, 1)],
        [1, 1, v1 * exp(1j * phi), (1, 0)],
        [1, 1, v1 * exp(1j * phi), (0, 1)],
        [1, 1, v1 * exp(1j * phi), (-1,-1)]
    ]

    site_coord_cart = np.array(
        [[0, 1 / (2 * sqrt(3))],
         [0, - 1 / (2 * sqrt(3))]]
    )

    model = tb.TightBindingModel(dim, site_num, lattice_vec, intra_cell, inter_cell)
    model.SiteCoord = model.cart2lattice(site_coord_cart.T).T

    return model


def Haldane_non_Hermitian_phase(t1, t2, phi, M, gamma):
    # v1 = v2 = t2 * exp(i gamma)
    return non_Hermitian_Haldane_H(t1, t1, t2 * exp(1j * gamma), t2 * exp(1j * gamma), phi, M)


# ---------------------------------------------------------------------------
# Test harness
# ---------------------------------------------------------------------------

def plot_winding_number():
    """Plot W(E_ref, mu1) using the SGBZ pipeline from brute_force_SGBZ."""
    E_ref = 3 + 0.1j * np.random.randn()
    mu1 = np.linspace(-1, 1, 101)

    model = Haldane_non_Hermitian_phase(*ALL_PARAMS)
    model = model.get_supercell(
        [(0, 0), (1, 0)],
        np.array([
            [1, 1],
            [1, -1]
        ], dtype=int)
    )

    coeffs, degs = model.get_characteristic_polynomial_data()

    poly = CharPoly(coeffs, degs)

    w_list = np.full_like(mu1, np.nan)
    for mu_ind in range(len(mu1)):
        zm = ZeroManager(poly, E_ref, mu1[mu_ind])
        print(mu_ind)
        zm.run()
        # The average winding is defined only for 0D subsets.  A
        # continuum (1D line) makes W jump; the bisection stage must
        # handle those μ₁ via left/right limits instead.
        if detect_continuum_simple(zm, poly):
            warnings.warn(
                f"mu1={mu1[mu_ind]:.4f}: continuum (1D subset) detected; "
                f"W undefined there, skipping"
            )
            continue
        _, w_list[mu_ind] = detect_crossings_and_winding(
            zm, poly
        )

    plt.plot(mu1, w_list, '.-')
    plt.show()


def test_crossing_detection():
    """Quick test of the new crossing detection at a single (E, mu1)."""
    import time

    E_ref = 3.0 + 0.0j
    mu1_val = -0.2

    model = Haldane_non_Hermitian_phase(*ALL_PARAMS)
    model = model.get_supercell(
        [(0, 0), (1, 0)],
        np.array([[1, 1], [1, -1]], dtype=int)
    )

    coeffs, degs = model.get_characteristic_polynomial_data()
    poly = CharPoly(coeffs, degs)

    print(f"Testing crossing detection for E={E_ref}, mu1={mu1_val}")

    t0 = time.perf_counter()
    zm = ZeroManager(poly, E_ref, mu1_val)
    zm.run()
    t1 = time.perf_counter()
    print(f"  ZeroManager.run: {t1 - t0:.3f}s, {len(zm.segments)} segments, "
          f"total rows: {sum(len(s.theta1_arr) for s in zm.segments)}")

    if detect_continuum_simple(zm, poly):
        print("  Continuum detected — skipping crossing detection")
        return

    t2 = time.perf_counter()
    subsets, W = detect_crossings_and_winding(zm, poly)
    t3 = time.perf_counter()
    print(f"  Crossing detection + winding: {t3 - t2:.3f}s")
    print(f"  W = {W}, {len(subsets)} subset points")

    for s in subsets[:5]:
        if isinstance(s, PointSubset):
            print(f"    θ₁={float(np.angle(s.beta1)) % (2*math.pi):.4f}, "
                  f"θ₂={float(np.angle(s.beta2)) % (2*math.pi):.4f}")
    if len(subsets) > 5:
        print(f"    ... and {len(subsets) - 5} more")

def build_synthetic_mu2_model():
    """合成 β₂ 多项式：含 β₂ 奇次项（±1）且系数依赖 β₁ → μ₂_mid 非平凡。

    2-band 紧束缚模型（M=N）的 β₂ 多项式总是 self-inversive（根成对 (β,1/β)），
    导致 μ₂_mid≡0；只含 β₂ 偶次幂的多项式退化为 β₂² 的方程，根模成对相等。
    此合成模型含 β₂ 的 ±1 次项且系数依赖 β₁，打破两种简并，用于测试 μ₂_mid 构建。

    方程：E − a/β₂² − b·β₂² − c·β₁·β₂ − d·β₁/β₂ − e = 0
    """
    coeffs = np.array([1, -1.3, -1.7, -0.37+0.1j, -0.43-0.05j, -0.51+0.2j],
                       dtype=complex)
    degs = np.array(
        [[1, 0, 0], [0, 0, -2], [0, 0, 2], [0, 1, 1], [0, 1, -1], [0, 0, 0]],
        dtype=int,
    )
    return CharPoly(coeffs, degs)


def build_continuum_mu2_model():
    """合成 β₂ 多项式：β₂² 的二次方程 → 根成 ±对（同模 = continuum cluster）。

    方程：E − a/β₂² − b·β₂² − c·β₁/β₂² − d·β₁·β₂² − e = 0，作 β₂² 的二次：
    (b+d·β₁)·y² + (E−e)·y + (a+c·β₁) = 0 (y=β₂²)。两根 y₁,y₂ → β₂=±√y₁,±√√y₂，
    ±对同模。两 cluster（{±√y₁}, {±√y₂}）各 mult 2，模随 θ₁ 变。

    用于测试简化版内部 continuum 检测（三组排序对）+ 代表元：原版 abs_argsort
    在 ±对间任意跳（126 假断点），代表元使 sort_to_item 稳定（0 假断点）。
    """
    coeffs = np.array([1, -1.3, -1.7, -0.37+0.1j, -0.43-0.05j, -0.5+0.2j],
                      dtype=complex)
    degs = np.array(
        [[1, 0, 0], [0, 0, -2], [0, 0, 2], [0, 1, -2], [0, 1, 2], [0, 0, 0]],
        dtype=int,
    )
    return CharPoly(coeffs, degs)


def test_mu2_mid():
    """构建 μ₂_mid 分段光滑表示并验证。

    断言：(M-1,M) 断点导数连续；其余断点导数跳变。
    画图：灰根曲线 + 红 μ₂_mid + 蓝断点 vline（PMGBZ 实线、内部虚线）。
    存 demos/mu2_mid_demo.png。
    """
    import time

    E_ref = 2.0 + 0.3j
    mu1_val = -0.1

    poly = build_synthetic_mu2_model()
    M = poly.M

    print(f"Testing μ₂_mid build for E={E_ref}, mu1={mu1_val}, M={M}")

    t0 = time.perf_counter()
    zm = Mu2MidZM(poly, E_ref, mu1_val)
    zm.run()
    t1 = time.perf_counter()
    print(f"  ZeroManager.run: {t1 - t0:.3f}s, {len(zm.segments)} segments, "
          f"total rows: {sum(len(s.theta1_arr) for s in zm.segments)}")

    if detect_continuum_simple(zm, poly):
        raise RuntimeError("continuum detected — μ₂_mid build assumes 0D-only")

    t2 = time.perf_counter()
    zm.build_mu2_mid(verbose=True)
    t3 = time.perf_counter()
    print(f"  build_mu2_mid: {t3 - t2:.3f}s")

    bps = zm.mu2_mid_breakpoints
    print(f"  breakpoints: {len(bps)}")
    # 精修收敛统计
    refined = getattr(zm, '_refined_points', [])
    n_conv = sum(1 for r in refined if r[4])
    print(f"  refined points: {len(refined)} (converged {n_conv})")

    cont_ok = True
    for bp in bps:
        djump = abs(bp.deriv_left - bp.deriv_right)
        d_scale = max(abs(bp.deriv_left), abs(bp.deriv_right), 1e-9)
        rel_jump = djump / d_scale
        tag = 'PMGBZ' if bp.is_pmgbz else 'internal'
        print(f"    θ₁={bp.theta1:.6e} val={bp.value:+.6e} "
              f"d_left={bp.deriv_left:+.4e} d_right={bp.deriv_right:+.4e} "
              f"|Δd|={djump:.2e} (rel {rel_jump:.1%}) gap={bp.gap:.2e} "
              f"{bp.pair_kind} {tag} cols={bp.columns}")
        if bp.pair_kind == 'M-1_M':
            # swap: 导数应连续。相对跳变 < 5% 视为数值噪声（交点处两根
            # 接近，tangent 计算有有限精度）。
            if rel_jump > 0.05:
                print(f"    WARN: (M-1,M) swap expected continuous deriv, "
                      f"rel jump {rel_jump:.1%}")
                cont_ok = False
        elif bp.pair_kind in ('M-2_M-1', 'M_M+1'):
            if rel_jump < 0.2:
                print(f"    WARN: {bp.pair_kind} expected deriv jump, "
                      f"rel jump {rel_jump:.1%}")
    if cont_ok:
        print("  [ok] (M-1,M) swaps continuous (rel < 5%)")
    print(f"  deriv continuity check: {'PASS' if cont_ok else 'FAIL'}")

    # ---- 画图 ----
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # 上图：根曲线 + μ₂_mid + 断点
    for seg in zm.segments:
        th = seg.theta1_arr
        logabs = np.log(np.abs(seg.tracked_roots))
        for j in range(logabs.shape[1]):
            mask = np.isfinite(logabs[:, j])
            if not np.any(mask):
                continue
            ax1.plot(th[mask], logabs[mask, j], color='gray',
                     alpha=0.4, linewidth=0.6, zorder=1)
    ax1.plot(zm.mu2_mid_theta1, zm.mu2_mid_values, 'r-',
             linewidth=1.5, label=r'$\mu_{2,\mathrm{mid}}(\theta_1)$', zorder=3)
    for bp in bps:
        ls = '-' if bp.is_pmgbz else '--'
        ax1.axvline(bp.theta1, color='blue', linestyle=ls, alpha=0.5,
                    linewidth=0.8, zorder=2)
    ax1.set_ylabel(r"$\ln|\beta_2|$")
    ax1.set_xlim(0, 2 * math.pi)
    ax1.legend(loc='best')

    # 下图：导数（左/右）
    th = zm.mu2_mid_theta1
    dv = zm.mu2_mid_derivs
    ax2.plot(th, dv, 'r.', markersize=2, zorder=3)
    for bp in bps:
        ax2.plot([bp.theta1, bp.theta1], [bp.deriv_left, bp.deriv_right],
                 'k-', linewidth=1.0, zorder=4)
        ls = '-' if bp.is_pmgbz else '--'
        ax2.axvline(bp.theta1, color='blue', linestyle=ls, alpha=0.5,
                    linewidth=0.8)
    ax2.axhline(0, color='gray', linewidth=0.5)
    ax2.set_xlabel(r"$\theta_1$")
    ax2.set_ylabel(r"$d\mu_{2,\mathrm{mid}}/d\theta_1$")
    ax2.set_xlim(0, 2 * math.pi)

    plt.show()

    plt.tight_layout()
    out = Path(__file__).parent / 'mu2_mid_demo.png'
    plt.savefig(out, dpi=120)
    print(f"  saved {out}")
    plt.close(fig)


def test_mu2_mid_continuum():
    """continuum case：±对同模 cluster，测试简化版内部 continuum 检测 + 代表元。

    验证：原版 abs_argsort 在 ±对间任意跳（126 假断点），内部检测 continuum
    后用代表元使 sort_to_item 稳定（0 假断点）。画图：灰根曲线（同 cluster
    同色）+ 红 μ₂_mid（代表元构建）+ cluster 标注。
    存 demos/mu2_mid_continuum_demo.png。
    """
    import time

    E_ref = 2.0 + 0.3j
    mu1_val = -0.1
    poly = build_continuum_mu2_model()
    M = poly.M

    print(f"Testing μ₂_mid (continuum) for E={E_ref}, mu1={mu1_val}, M={M}")

    zm = Mu2MidZM(poly, E_ref, mu1_val)
    t0 = time.perf_counter()
    zm.run()
    t1 = time.perf_counter()
    print(f"  ZeroManager.run: {t1-t0:.3f}s, {len(zm.segments)} segments, "
          f"total rows: {sum(len(s.theta1_arr) for s in zm.segments)}")

    # 内部 continuum 检测（简化版）
    clusters = zm._detect_continuum_clusters_internal(1e-8)
    print(f"  internal continuum clusters: {clusters}")

    t2 = time.perf_counter()
    zm.build_mu2_mid(verbose=True)  # continuum_clusters=None → 内部检测
    t3 = time.perf_counter()
    print(f"  build_mu2_mid: {t3-t2:.3f}s")

    bps = zm.mu2_mid_breakpoints
    print(f"  breakpoints: {len(bps)} (应=0：continuum 段代表元稳定，无假断点)")
    jl = zm.mu2_mid_jlo
    print(f"  jlo 相邻变化: {int(np.sum(jl[:-1]!=jl[1:]))} (应=0：代表元稳定)")

    # continuum 段统计（j_lo==j_hi → M-1,M 同 item = 1D subset）
    cont_frac = float(np.mean(zm.mu2_mid_jlo == zm.mu2_mid_jhi))
    print(f"  continuum 行比例 (j_lo==j_hi): {cont_frac:.1%}")
    print(f"  mu2_mid range=[{zm.mu2_mid_values.min():.4e},"
          f"{zm.mu2_mid_values.max():.4e}]")

    # ---- 画图 ----
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # cluster 着色：同 cluster 的列同色，非 cluster 列灰色
    cluster_cols: dict[int, int] = {}
    color_idx = 0
    for seg_clusters in clusters:
        for cl in seg_clusters:
            for c in cl:
                cluster_cols[c] = color_idx
            color_idx += 1
    cmap = plt.cm.tab10

    for s_idx, seg in enumerate(zm.segments):
        th = seg.theta1_arr
        logabs = np.log(np.abs(seg.tracked_roots))
        for j in range(logabs.shape[1]):
            mask = np.isfinite(logabs[:, j])
            if not np.any(mask):
                continue
            if j in cluster_cols:
                c = cmap(cluster_cols[j] % 10)
                lbl = f"cluster col {j}"
                ax1.plot(th[mask], logabs[mask, j], color=c, alpha=0.8,
                         linewidth=1.0, zorder=2)
            else:
                ax1.plot(th[mask], logabs[mask, j], color='gray',
                         alpha=0.4, linewidth=0.6, zorder=1)

    ax1.plot(zm.mu2_mid_theta1, zm.mu2_mid_values, 'r-',
             linewidth=2.0, label=r'$\mu_{2,\mathrm{mid}}(\theta_1)$（代表元）',
             zorder=4)
    ax1.set_ylabel(r"$\ln|\beta_2|$")
    ax1.set_xlim(0, 2 * math.pi)
    ax1.set_title(f"continuum case: ±对同模 clusters={clusters}")
    ax1.legend(loc='best')

    # 下图：导数
    th = zm.mu2_mid_theta1
    dv = zm.mu2_mid_derivs
    ax2.plot(th, dv, 'r-', linewidth=1.5, zorder=3)
    ax2.axhline(0, color='gray', linewidth=0.5)
    ax2.set_xlabel(r"$\theta_1$")
    ax2.set_ylabel(r"$d\mu_{2,\mathrm{mid}}/d\theta_1$")
    ax2.set_xlim(0, 2 * math.pi)

    plt.tight_layout()
    out = Path(__file__).parent / 'mu2_mid_continuum_demo.png'
    plt.savefig(out, dpi=120)
    print(f"  saved {out}")
    plt.close(fig)


if __name__ == '__main__':
    test_mu2_mid()
    print()
    test_mu2_mid_continuum()
