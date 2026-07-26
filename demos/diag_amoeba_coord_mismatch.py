'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-24
Copyright © Department of Physics, Tsinghua University. All rights reserved

Minimal reproduction + diagnosis of the amoeba-spectrum mismatch between the
gain-loss Haldane model in two coordinate systems:

  - BASE : the primitive 2-site cell  -> data/Haldane-gain-loss-amoeba.pkl
  - XY   : a 2-site supercell with basis [[1,1],[1,-1]] -> ...-xy.pkl

Both describe the SAME physical model, so the amoeba spectrum (whether E is in
the GBZ) must agree.  87 grid points disagree on the PointSubset *count*; of
those, 5 are genuine spectrum mismatches (BASE says outside, XY says inside).

This script reproduces one spectrum-mismatch point (E = 0.365) and tests the
hypothesis:

  At the solved (mu1, mu2) = (0, 0) the BASE root tracks approach ln|b2| = 0
  from opposite sides but, on the adaptive ZM grid, never cross (gap ~0.002).
  The XY tracks (same physics, reparameterised theta1) do cross.  If the gap
  is a sampling artifact, then refining the theta1 mesh around the closest
  approach should make the BASE track cross mu2 = 0 (or touch it), recovering
  the missing two points.

The test:
  1. Run the bisection + ZM for BASE and XY at E=0.365, print the per-track
     ln|beta2| ranges (the tell-tale "touch but don't cross" signature).
  2. For BASE, find the theta1 of the closest approach of the two near-zero
     tracks, then re-solve beta2 roots on a fine uniform theta1 grid straddling
     that point, and check whether ln|beta2| - 0 changes sign (cross) or
     reaches 0 (touch).
'''

import sys
from pathlib import Path
import numpy as np
from math import sin, cos, sqrt, pi
from cmath import exp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from BerryPy import TightBinding as tb
from gbz_types import CharPoly
from brute_force_amoeba.bisect import bisect_amoeba_ronkin_min
from brute_force_amoeba.zm_extract import AmoebaZeroManager
from brute_force_amoeba.zm_extract import amoeba_windings
import brute_force_amoeba as bfa


# Haldane gain-loss model, same params as demos/Haldane-model-gainloss.py
ALL_PARAMS = (1, 0.5, pi / 3, 0.5j, 0)


def non_Hermitian_Haldane_H(u1, u2, v1, v2, phi, M):
    lattice_vec = np.array([[-cos(pi / 3), -cos(pi / 3)],
                           [-sin(pi / 3), sin(pi / 3)]])
    intra = [[0, 0, M], [1, 1, -M], [1, 0, u1], [0, 1, u2]]
    inter = [
        [0, 0, v2 * exp(1j * phi), (-1, 0)],
        [0, 0, v2 * exp(1j * phi), (0, -1)],
        [0, 0, v2 * exp(1j * phi), (1, 1)],
        [0, 0, v1 * exp(-1j * phi), (1, 0)],
        [0, 0, v1 * exp(-1j * phi), (0, 1)],
        [0, 0, v1 * exp(-1j * phi), (-1, -1)],
        [1, 0, u1, (0, -1)],
        [1, 0, u1, (1, 0)],
        [0, 1, u2, (0, 1)],
        [0, 1, u2, (-1, 0)],
        [1, 1, v2 * exp(-1j * phi), (-1, 0)],
        [1, 1, v2 * exp(-1j * phi), (0, -1)],
        [1, 1, v2 * exp(-1j * phi), (1, 1)],
        [1, 1, v1 * exp(1j * phi), (1, 0)],
        [1, 1, v1 * exp(1j * phi), (0, 1)],
        [1, 1, v1 * exp(1j * phi), (-1, -1)],
    ]
    return tb.TightBindingModel(2, 2, lattice_vec, intra, inter)


def Haldane_non_Hermitian_phase(t1, t2, phi, M, gamma):
    return non_Hermitian_Haldane_H(t1, t1, t2 * exp(1j * gamma),
                                   t2 * exp(1j * gamma), phi, M)


def model_base():
    return Haldane_non_Hermitian_phase(*ALL_PARAMS)


def model_xy():
    m = Haldane_non_Hermitian_phase(*ALL_PARAMS)
    return m.get_supercell([(0, 0), (1, 0)],
                           np.array([[1, 1], [1, -1]], dtype=int))


def per_track_ranges(zm, mu2):
    """Per-track ln|beta2| range, sign-change count, exact-zero count."""
    la = zm.seg_logabs[0]
    out = []
    d = la - mu2
    sc = (d[:-1] * d[1:]) < 0
    for j in range(la.shape[1]):
        col = la[:, j]
        out.append({
            'min': float(col.min()), 'max': float(col.max()),
            'n_cross': int(sc[:, j].sum()),
            'n_exact0': int((d[:, j] == 0).sum()),
        })
    return out


def solve_beta2_logabs_on_grid(poly, E, mu1, theta1_grid, mu2):
    """Solve beta2 roots on a uniform theta1 grid (sorted by |beta2|) and
    return ln|beta2| arrays.  Each column is a track sorted by modulus at
    that sample — NOT Hungarian-matched across samples.  Good enough to see
    whether ANY root touches/crosses mu2 on the refined grid."""
    n = len(theta1_grid)
    all_logabs = []
    for t1 in theta1_grid:
        beta1 = exp(mu1 + 1j * t1)
        roots = np.asarray(poly.solve_roots_1d((0, 1), (E, beta1), (2,)),
                           dtype=complex)
        all_logabs.append(np.sort(np.log(np.abs(roots))))
    return np.array(all_logabs)  # (n, K), columns sorted by modulus per row


def main():
    E = 0.365 + 0j
    print("=" * 70)
    print(f"Reproduction: E = {E}")
    print("=" * 70)

    results = {}
    for name, model_fn in [('BASE', model_base), ('XY', model_xy)]:
        coeffs, degs = model_fn().get_characteristic_polynomial_data()
        poly = CharPoly(coeffs, degs)
        res = bisect_amoeba_ronkin_min(poly, E, N_points=301)
        zm = res['_zm']
        mu1, mu2 = res['mu1'], res['mu2']
        w2, zs, hc, _ = amoeba_windings(zm, poly, E, mu1, mu2, refine=False)
        gbz = bfa.collect_GBZ_subsets(coeffs, degs, E, 0.0, True)
        results[name] = dict(poly=poly, zm=zm, mu1=mu1, mu2=mu2,
                             w2=w2, n_cross=len(zs), gbz=gbz)
        print(f"\n[{name}] mu1={mu1:.6f} mu2={mu2:.6f} K={zm.K} "
              f"w2={w2} n_cross={len(zs)} gbz_idx={gbz.index} is_gbz={gbz.is_gbz}")
        for j, r in enumerate(per_track_ranges(zm, mu2)):
            print(f"    track {j}: ln|b2| in [{r['min']:+.4f}, {r['max']:+.4f}]  "
                  f"cross={r['n_cross']} exact0={r['n_exact0']}")

    # ---- refined-grid test on BASE ----
    # The two near-zero tracks (1 and 2) touch ln|b2|=0 around their closest
    # approach.  Find the sample where |ln|b2|| is smallest for the track that
    # approaches 0 from below, refine theta1 around it, and check for a cross.
    print("\n" + "=" * 70)
    print("Refined-grid test on BASE: does the near-touch actually cross?")
    print("=" * 70)

    zm = results['BASE']['zm']
    poly = results['BASE']['poly']
    mu1, mu2 = results['BASE']['mu1'], results['BASE']['mu2']
    th = zm.segments[0].theta1_arr
    la = zm.seg_logabs[0]

    # closest approach of the two tracks that bracket mu2 (the touching pair).
    # Identify the track with max < mu2 (approaches from below) and its min |max-mu2|.
    below = [j for j in range(la.shape[1]) if la[:, j].max() < mu2]
    above = [j for j in range(la.shape[1]) if la[:, j].min() > mu2]
    print(f"  tracks below mu2 (max<mu2): {below}")
    print(f"  tracks above mu2 (min>mu2): {above}")

    # gap = min over below-track of (mu2 - max) + min over above of (min - mu2)
    # The touching pair: below track whose max is closest to mu2,
    # above track whose min is closest to mu2.
    jb = max(below, key=lambda j: la[:, j].max()) if below else None
    ja = min(above, key=lambda j: la[:, j].min()) if above else None
    if jb is None or ja is None:
        print("  No below/above pair found; cannot localise touch.")
        return
    ib = int(np.argmax(la[:, jb]))   # sample where below track is highest
    ia = int(np.argmin(la[:, ja]))  # sample where above track is lowest
    print(f"  touching pair: below track {jb} (max={la[ib, jb]:+.6f} @th1={th[ib]:.4f}), "
          f"above track {ja} (min={la[ia, ja]:+.6f} @th1={th[ia]:.4f})")
    print(f"  gap below->mu2 = {mu2 - la[ib, jb]:.6e}, gap above->mu2 = {la[ia, ja] - mu2:.6e}")

    # Refine theta1 on a fine grid straddling BOTH closest-approach samples.
    center = 0.5 * (th[ib] + th[ia]) if ib != ia else th[ib]
    half = 0.5 * abs(th[1] - th[0]) * 3   # a few ZM-step widths each side
    th_fine = np.linspace(center - half, center + half, 2001)
    la_fine = solve_beta2_logabs_on_grid(poly, E, mu1, th_fine, mu2)
    # per-row sorted; check whether any row has a root with ln|b2| crossing 0
    # i.e. min ln|b2| <= mu2 <= max ln|b2| on the fine grid
    d_fine = la_fine - mu2
    below_any = (d_fine < 0).any(axis=1)
    above_any = (d_fine > 0).any(axis=1)
    both = below_any & above_any          # a sample where roots straddle mu2
    n_touch = int(np.sum(np.abs(d_fine).min(axis=1) < 1e-12))

    print(f"  refined grid: {len(th_fine)} pts in [{th_fine[0]:.4f}, {th_fine[-1]:.4f}]")
    print(f"  min |ln|b2| - mu2| over refined grid = {np.abs(d_fine).min():.6e}")
    print(f"  samples with a root touching mu2 (|d|<1e-12): {n_touch}")
    print(f"  samples where roots straddle mu2 (both signs present): {int(both.sum())}")

    # sign-change detection per column on the refined grid (sorted columns are
    # not tracks, so look at the envelope: does the multiset of ln|b2| cross?)
    # Robust check: at each refined sample, count roots below mu2.  A crossing
    # shows up as the count changing between adjacent samples.
    n_below = (d_fine < 0).sum(axis=1)
    n_equal = (d_fine == 0).sum(axis=1)
    count_changes = int(np.sum(np.abs(np.diff(n_below)) > 0))
    print(f"  refined: #{'roots below mu2'} changes between adjacent samples: {count_changes}")
    print(f"  refined: #samples with an exact-zero root: {int(n_equal.sum())}")

    print()
    print("Interpretation:")
    print("  - If min|ln|b2|-mu2| ~ 0 AND n_equal>0 on the fine grid: it is a")
    print("    true TOUCH (degenerate at |b2|=exp(mu2)).  BASE's coarse ZM grid")
    print("    missed it -> sampling artifact, confirms the hypothesis.")
    print("  - If count_changes>0 on the fine grid: it is a true CROSS that the")
    print("    ZM grid stepped over -> sampling artifact, confirms the hypothesis.")
    print("  - If the gap stays open (~1e-3) on the fine grid: it is a genuine")
    print("    non-crossing; XY's crossing is the artifact (or a different point).")


if __name__ == "__main__":
    main()
