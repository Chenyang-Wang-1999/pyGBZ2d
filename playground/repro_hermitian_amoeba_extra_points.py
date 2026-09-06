'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-09-03
Copyright © Department of Physics, Tsinghua University. All rights reserved

Minimal reproduction: amoeba used to return extra PointSubsets for the
Hermitian Haldane test model at real energies E = 3.5, 4.0, 5.0 (up to 70
spurious subsets).  Fixed on 2026-09-05 by the interval-trigger per-pair
tracking change (see log/2026-09-05-interval-trigger-per-pair-tracking.md):
the simultaneous double MR near theta1 = 4.5787 (E = 3.5) had been missed
because the closest-pair identity flickers between the two simultaneously
degenerating pairs and the same-pair guard rejected the derivative flip.

The model is exactly the one built in ``Hermitian_Haldane_sweep`` of
``playground/Haldane-model-gainloss.py`` (M=3, gamma=0, plus four diagonal
InterCell terms), but we do not do any special handling: the same
characteristic polynomial is passed to both amoeba and SGBZ.

Run::

    python playground/repro_hermitian_amoeba_extra_points.py

Expected repro summary (observed on this checkout):
    E = 3.5+0j : SGBZ index=(0,6) | amoeba index=(0,6)
    E = 4.0+0j : SGBZ index=(0,6) | amoeba index=(0,6)
    E = 5.0+0j : SGBZ index=(0,6) | amoeba index=(0,6)
'''

from __future__ import annotations

import sys
from cmath import exp
from math import cos, pi, sin, sqrt
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Repo layout: this file lives in playground/, package lives in src/.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import BerryPy.TightBinding as tb

import pygbz2d.amoeba as bfa
import pygbz2d.sgbz as bfs
from pygbz2d.core import PointSubset, LineSubset, CharPoly
from pygbz2d.continuation import ZeroManager


def non_Hermitian_Haldane_H(u1, u2, v1, v2, phi, M):
    """Identical to the builder in playground/Haldane-model-gainloss.py."""
    dim = 2
    site_num = 2

    lattice_vec = np.array(
        [[-cos(pi / 3), -cos(pi / 3)],
         [-sin(pi / 3), sin(pi / 3)]]
    )

    intra_cell = [
        [0, 0, M],
        [1, 1, -M],
        [1, 0, u1],
        [0, 1, u2]
    ]
    inter_cell = [
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
        [1, 1, v1 * exp(1j * phi), (-1, -1)]
    ]

    site_coord_cart = np.array(
        [[0, 1 / (2 * sqrt(3))],
         [0, -1 / (2 * sqrt(3))]]
    )

    model = tb.TightBindingModel(dim, site_num, lattice_vec, intra_cell, inter_cell)
    model.SiteCoord = model.cart2lattice(site_coord_cart.T).T
    return model


def build_hermitian_haldane_model():
    """Hermitian_Haldane_sweep model, including the four InterCell additions."""
    model = non_Hermitian_Haldane_H(1, 1, 0.5, 0.5, pi / 3, 3)
    model.InterCell += [
        (0, 0, 2, (1, 0)),
        (0, 0, 2, (-1, 0)),
        (1, 1, 2, (1, 0)),
        (1, 1, 2, (-1, 0)),
    ]
    return model


def main():
    model = build_hermitian_haldane_model()
    coeffs, degs = model.get_characteristic_polynomial_data()

    char_poly = CharPoly(coeffs, degs)
    print(f"characteristic polynomial: {len(coeffs)} terms, "
          f"degs shape {degs.shape}, E degree range "
          f"{int(degs[:, 0].min())}..{int(degs[:, 0].max())}")

    for E in (3.5 + 0j, 4.0 + 0j, 5.0 + 0j):
        amoeba_res = bfa.collect_GBZ_subsets(coeffs, degs, E)
        sgbz_res = bfs.collect_GBZ_subsets(coeffs, degs, E)

        a_counts = _count_subsets(amoeba_res.subsets)
        s_counts = _count_subsets(sgbz_res.subsets)
        print(f"\nE = {E.real:g}{E.imag:+.1f}j")
        print(f"  amoeba: success={amoeba_res.success} is_gbz={amoeba_res.is_gbz} "
              f"index={amoeba_res.index} "
              f"subset_counts={a_counts}")
        print(f"  sgbz  : success={sgbz_res.success} is_gbz={sgbz_res.is_gbz} "
              f"index={sgbz_res.index} "
              f"subset_counts={s_counts}")

        color_list = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
        zm = ZeroManager(char_poly, E, 0)
        zm.run()
        for seg_idx, seg in enumerate(zm.segments):
            plt.plot(seg.theta1_arr, np.log(np.abs(seg.tracked_roots)), color_list[seg_idx % len(color_list)])
        for s in amoeba_res.subsets:
            if isinstance(s, PointSubset):
                plt.plot(np.angle(s.beta1) % (2 * pi), np.log(np.abs(s.beta2)), 'x')
            else:
                plt.plot(s.theta1_arr, np.log(np.abs(s.beta2_arr)), ':')
        plt.show()


    # Keep the script importable without side effects.
    return 0


def _count_subsets(subsets):
    n_points = sum(1 for s in subsets if isinstance(s, PointSubset))
    n_lines = sum(1 for s in subsets if isinstance(s, LineSubset))
    return {"PointSubset": n_points, "LineSubset": n_lines}


if __name__ == "__main__":
    raise SystemExit(main())
