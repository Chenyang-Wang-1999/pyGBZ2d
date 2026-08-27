'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-22
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

''' Test ZeroManager with some simple examples '''

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import bfgbz2d.core as gbz_types
import bfgbz2d.continuation as continuation

def test_winding_around_0():
    mu1 = 0.2
    c0 = 2 * np.cosh(mu1)
    coeffs = [
        1, -1, -1, -c0
    ]
    degs = np.array(
        [
            [0, 0, 2],
            [0, 1, 0],
            [0, -1, 0],
            [0, 0, 0]
        ],
        dtype=int,
    )
    char_poly = gbz_types.CharPoly(coeffs, degs)

    zm = continuation.ZeroManager(
        char_poly, 0, mu1
    )
    zm.run(verbose=True)

    ## print basic information ##
    print(len(zm.segments))
    print(len(zm.multiple_roots))
    print(zm.has_boundary_mr)
    print(zm.boundary_perm)

    print("Left boundary: ", zm.left_boundary_roots)
    print("segment[0]: ", zm.segments[0].tracked_roots[0, :])
    print("segment[-1]: ", zm.segments[-1].tracked_roots[-1, :])

    ## segments ##
    for seg in zm.segments:
        print(f"seg: from {seg.left_mr} to {seg.right_mr}")
        # plt.plot(seg.theta1_arr, np.abs(seg.tracked_roots))
        plt.plot(seg.tracked_roots.real, seg.tracked_roots.imag, '.-')
    
    ## roots ##
    for mr in zm.multiple_roots:
        plt.plot(mr.roots.real, mr.roots.imag, 'x')
    plt.show()


def test_iterative_solver():
    mu1 = 0.2
    c0 = 2 * np.cosh(mu1)
    coeffs = [
        1, -1, -1, -c0
    ]
    degs = np.array(
        [
            [0, 0, 2],
            [0, 1, 0],
            [0, -1, 0],
            [0, 0, 0]
        ],
        dtype=int,
    )
    char_poly = gbz_types.CharPoly(coeffs, degs)

    continuation.multiple_roots.solve_multiple_roots_iterative(
        char_poly, 0, np.exp(mu1), 0.0003
    )


if __name__ == "__main__":
    test_winding_around_0()
    # test_iterative_solver()
