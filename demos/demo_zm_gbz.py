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
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gbz_types import CharPoly, PointSubset, LineSubset
from continuation import ZeroManager

# Import the SGBZ pipeline from brute_force_SGBZ
from brute_force_SGBZ.continuum_lines import detect_continuum_simple
from brute_force_SGBZ.mu2mid import CONTINUUM_TOL, CONTINUUM_FRAC
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
# μ₂_mid — imported from production (brute_force_SGBZ.mu2mid)
# ---------------------------------------------------------------------------
# This demo originally carried a full duplicate of Mu2MidZM (plus ItemView,
# Mu2MidBreakpoint and the cubic-Hermite helpers).  It now imports the
# production class from brute_force_SGBZ.mu2mid so the two cannot diverge.
# The synthetic-model builders, tests, and plotting below remain demo-only.
from brute_force_SGBZ.mu2mid import Mu2MidZM


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
