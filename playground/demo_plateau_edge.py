'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-24
Copyright © Department of Physics, Tsinghua University. All rights reserved

Debug demo: the zero-plateau edge of a next-nearest-coupling HN model.

Reproduces the ``TestPlateauEdge`` case from ``tests/test_amoeba.py``:
the next-nearest-coupling model at E = 1.648 + 0.0294j sits on the edge of a
zero-w2 plateau (μ₂ = 0).  With the ZeroManager-based amoeba solver this is
resolved directly:

  - two β₂ tracks merely *touch* |β₂| = exp(μ₂) = 1 from opposite sides
    (one from below, one from above) without crossing;
  - w2 is therefore uniformly 0 and there are no discrete crossings;
  - the point is correctly classified as non-GBZ.

The old fixed-grid solver instead produced 4 canceling spurious zeros at the
edge; the functional outcome (edge → non-GBZ) is unchanged.

Run this script to inspect the bisection result, the ZM tracks, the winding,
and the final GBZ classification step by step.
'''

import sys
from pathlib import Path
import numpy as np
from cmath import exp, pi

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from BerryPy import TightBinding as tb

from bfgbz2d.core import CharPoly, PointSubset, LineSubset
import bfgbz2d.amoeba as bfa
from bfgbz2d.amoeba.bisect import bisect_amoeba_ronkin_min
from bfgbz2d.amoeba.zm_extract import AmoebaZeroManager
from bfgbz2d.amoeba.zm_extract import amoeba_windings
from bfgbz2d.amoeba.ronkin_winding import _get_average_winding_from_zeros
from bfgbz2d.amoeba.amoeba import _check_zeros_are_clustered


E_PLATEAU_EDGE = 1.648 + 0.0294j
PLATEAU_AREA_THRESHOLD = 1e-2


def build_nnc_char_poly():
    """Next-nearest-coupling model (mirrors tests/test_amoeba.py::nnc_char_poly)."""
    gamma = 0.2
    u1, v1, w = 1.0, 0.8, 0.5
    u2, v2 = 0.1, 0.1
    u3, v3 = 0.05, 0.05
    inter_cell = [
        [0, 0, u1 * np.exp(gamma), (1, 0)],
        [0, 0, u1 * np.exp(-gamma), (-1, 0)],
        [0, 0, u2, (2, 0)],
        [0, 0, u2, (-2, 0)],
        [0, 0, u3, (3, 0)],
        [0, 0, u3, (-3, 0)],
        [0, 0, v1, (0, 1)],
        [0, 0, v1, (0, -1)],
        [0, 0, v2, (0, 2)],
        [0, 0, v2, (0, -2)],
        [0, 0, v3, (0, 3)],
        [0, 0, v3, (0, -3)],
        [0, 0, w, (1, 1)],
        [0, 0, w, (1, -1)],
        [0, 0, w, (-1, 1)],
        [0, 0, w, (-1, -1)],
    ]
    model = tb.TightBindingModel(2, 1, np.eye(2), [], inter_cell, [(0, 0)])
    coeffs, degs = model.get_characteristic_polynomial_data()
    return CharPoly(coeffs, degs), coeffs, degs


def main():
    char_poly, coeffs, degs = build_nnc_char_poly()
    E = E_PLATEAU_EDGE
    print("=" * 70)
    print(f"Plateau-edge debug: NNC model, E = {E}")
    print("=" * 70)

    # ---- 1. bisection finds (mu1, mu2) ----
    res = bisect_amoeba_ronkin_min(char_poly, E)
    mu1, mu2 = res["mu1"], res["mu2"]
    print(f"\n[1] bisection:")
    print(f"    mu1 = {mu1}")
    print(f"    mu2 = {mu2}")
    print(f"    is_continuum = {res['is_continuum']}")
    print(f"    exit_reason  = {res.get('_exit_reason')}")
    print(f"    n_zeros      = {len(res['zeros'])}")

    # ---- 2. reuse the ZM built inside the bisection ----
    zm = res["_zm"]
    print(f"\n[2] ZeroManager @ solved mu1:")
    print(f"    n_segments = {zm.n_segments}")
    print(f"    n_multiple_roots = {zm.n_multiple_roots}")
    for s in range(len(zm.segments)):
        seg = zm.segments[s]
        la = zm.seg_logabs[s]
        print(f"    seg{s}: {la.shape[0]} samples, {la.shape[1]} tracks, "
              f"th1=[{seg.theta1_arr[0]:.4f}, {seg.theta1_arr[-1]:.4f}]")

    # ---- 3. per-track ln|beta2| range at the solved mu2 ----
    #   A *touch* shows up as one track's max ln|beta2| == mu2 and another's
    #   min ln|beta2| == mu2, both approached from one side only.
    print(f"\n[3] per-track ln|beta2| range vs mu2 = {mu2:.6f}:")
    la = zm.seg_logabs[0]
    for j in range(la.shape[1]):
        col = la[:, j]
        d = col - mu2
        n_up = int(np.sum(d > 0))
        n_dn = int(np.sum(d < 0))
        side = ("all above" if n_dn == 0
                else "all below" if n_up == 0
                else "crosses")
        print(f"    track {j}: range [{col.min():+.4f}, {col.max():+.4f}]  "
              f"|d|_min={np.min(np.abs(d)):.3e}  {side}")

    # ---- 4. winding from ZM tracks ----
    w2, zeros, has_cont, dW = amoeba_windings(
        zm, char_poly, E, mu1, mu2, refine=False,
    )
    print(f"\n[4] amoeba_windings @ (mu1, mu2):")
    print(f"    w2           = {w2}")
    print(f"    n_zeros      = {len(zeros)}")
    print(f"    has_continuum= {has_cont}")

    # ---- 5. plateau pre-check conditions ----
    zeros_ref = res["zeros"]
    w1_area = res.get("_w1_area")
    if w1_area is None:
        _, w1_area = _get_average_winding_from_zeros(
            char_poly, E, mu1, mu2, zeros_ref, direction=1,
        )
    _, w2_area = _get_average_winding_from_zeros(
        char_poly, E, mu1, mu2, zeros_ref, direction=2,
    )
    clustered = (_check_zeros_are_clustered(zeros_ref, PLATEAU_AREA_THRESHOLD)
                 if len(zeros_ref) >= 2 else False)
    print(f"\n[5] plateau pre-check:")
    print(f"    w1_area      = {w1_area:.3e}  (< {PLATEAU_AREA_THRESHOLD}: {w1_area < PLATEAU_AREA_THRESHOLD})")
    print(f"    w2_area      = {w2_area:.3e}  (< {PLATEAU_AREA_THRESHOLD}: {w2_area < PLATEAU_AREA_THRESHOLD})")
    print(f"    zeros_clustered = {clustered}  (n_zeros = {len(zeros_ref)})")

    # ---- 6. final GBZ classification ----
    gbz = bfa.collect_GBZ_subsets(coeffs, degs, E, 0.0, plateau_check=True)
    print(f"\n[6] collect_GBZ_subsets:")
    print(f"    index        = {gbz.index}")
    print(f"    is_gbz       = {gbz.is_gbz}")
    print(f"    success      = {gbz.success}")
    n_pts = sum(1 for s in gbz.subsets if isinstance(s, PointSubset))
    n_lines = sum(1 for s in gbz.subsets if isinstance(s, LineSubset))
    print(f"    subsets      = {n_pts} PointSubset, {n_lines} LineSubset")
    print()
    print("Expectation: edge → non-GBZ (index (0,0), is_gbz=False).")
    print("With ZM tracking this holds because w2≡0 uniformly (touch, no crossing);")
    print("no clustering-based probe is needed.")


if __name__ == "__main__":
    main()
