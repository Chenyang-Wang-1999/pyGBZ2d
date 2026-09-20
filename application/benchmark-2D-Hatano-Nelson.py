'''
author:        Wang Chenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-09-15
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''


RUN_IN_SRC = True
#### Add source-file location to sys.path #####
if RUN_IN_SRC:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

#### libs ####
import numpy as np
from typing import Literal
import matplotlib.pyplot as plt

#### pyGBZ2d ####
from pygbz2d import core
from pygbz2d import amoeba
from pygbz2d import sgbz
from pygbz2d import continuation

##### Modeling with ChP #####
def get_HN_charpoly(
    Jx1: complex,
    Jx2: complex,
    Jy1: complex,
    Jy2: complex,
    which: Literal["x-y", "y-x", "11-y"]
):
    ''' TODO: add docstring '''

    if which == "x-y" or which == "y-x":
        coeffs = np.array([1, -Jx1, -Jx2, -Jy1, -Jy2])
        degs = np.array([
            # E, betax, betay
            [1, 0, 0],
            [0, -1, 0],
            [0, 1, 0],
            [0, 0, -1],
            [0, 0, 1],
        ])

        if which == "y-x":
            degs = degs[:, [0, 2, 1]]
    elif which == "11-y":
        coeffs = np.array([1, -Jx1, -Jx2, -Jy1, -Jy2])
        degs = np.array([
            # E, beta_([11]), beta_(y)
            [1, 0, 0],
            [0, -1, 1],
            [0, 1, -1],
            [0, 0, -1],
            [0, 0, 1]
        ])
    else:
        raise ValueError(f"Unknown which: {which}")

    return (coeffs, degs)


def calculate_subset_and_print(
    E_ref: complex,
    Jx1: complex,
    Jx2: complex,
    Jy1: complex,
    Jy2: complex,
    which: Literal["amoeba", "x-strip", "y-strip", "11-strip"],
):
    ''' TODO: add docstring '''
    poly_dict = {"x-strip": "x-y", "y-strip": "y-x", "11-strip": "11-y"}
    if which == "amoeba":
        coeffs, degs = get_HN_charpoly(Jx1, Jx2, Jy1, Jy2, "x-y")
        res = (amoeba.collect_GBZ_subsets(coeffs, degs, E_ref, debug_mode=True))
        poly = core.CharPoly(coeffs, degs)
        zm = continuation.ZeroManager(poly, E_ref, 0)
        zm.run()
        print(amoeba.zm_extract.find_crossings(zm, 0.0, -0.09116077839697725, return_refined=True))
        plt.figure()
        for seg in zm.segments:
            plt.plot(seg.theta1_arr, np.log(np.abs(seg.tracked_roots)))
        plt.show()
    else:
        coeffs, degs = get_HN_charpoly(Jx1, Jx2, Jy1, Jy2, poly_dict[which])
        print(sgbz.collect_GBZ_subsets(coeffs, degs, E_ref, debug_mode=True))


def demo_point_subsets():
    ''' Parameters used in arXiv: 2506.22743 '''
    Jx1 = 1 + 1j
    Jx2 = 1.5 + 1.2j
    Jy1 = -1 + 1j
    Jy2 = -1.2 - 0.5j
    E_ref = 1 + 1j
    for which in ["amoeba", "x-strip", "y-strip", "11-strip"]:
        print("="*10, which, "="*10)
        calculate_subset_and_print(E_ref, Jx1, Jx2, Jy1, Jy2, which)


def demo_line_subsets():
    ''' Hermitian model '''
    Jx1 = 1
    Jx2 = 1.5
    Jy1 = -1
    Jy2 = -1.2
    E_ref = 1
    for which in ["amoeba", "x-strip", "y-strip", "11-strip"]:
        print("="*10, which, "="*10)
        calculate_subset_and_print(E_ref, Jx1, Jx2, Jy1, Jy2, which)


def check():
    ''' Hermitian model '''
    Jx1 = 1
    Jx2 = 1.5
    Jy1 = -1
    Jy2 = -1.2
    E_ref = 1
    calculate_subset_and_print(E_ref, Jx1, Jx2, Jy1, Jy2, "amoeba")


if __name__ == "__main__":
    # demo_line_subsets()
    check()
