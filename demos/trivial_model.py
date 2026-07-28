'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-19
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

'''
A simple Hermitian model to test GBZ solving program.
H = tx * cos(kx + phix) + ty * cos(ky + phiy)

or equivalently

H = tx / 2 * exp(1j * phix) * betax + tx / 2 * exp(-1j * phix) / betax 
    + ty / 2 * exp(1j * phiy) * betay + ty / 2 * exp(-1j * phiy) / betay

The expected amoebic GBZ and SGBZ are identical to the BZ: kx, ky \in \mathbb{R}
'''

import BerryPy.TightBinding as tb
from cmath import exp, log, pi
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from continuation import ZeroManager
from gbz_types import CharPoly, LineSubset
import brute_force_amoeba as bfa


DEFAULT_PARAMS = {
    "tx": 1,
    "ty": 1,
    "phix": pi / 12,
    "phiy": 0,
}

def get_model(tx, ty, phix, phiy):
    dim = 2
    site_num = 1
    intra_cell = []
    inter_cell = [
        [0, 0, tx * exp(-1j * phix) / 2, (1, 0)],
        [0, 0, tx * exp(1j * phix) / 2, (-1, 0)],
        [0, 0, ty * exp(-1j * phiy) / 2, (0, 1)],
        [0, 0, ty * exp(1j * phiy) / 2, (0, -1)],
    ]
    return tb.TightBindingModel(
        dim, site_num, np.eye(2), intra_cell, inter_cell
    )


def solve_zeros():
    E_ref = 0.2
    model = get_model(**DEFAULT_PARAMS)
    tx, ty, phix, phiy = DEFAULT_PARAMS["tx"], DEFAULT_PARAMS["ty"], DEFAULT_PARAMS["phix"], DEFAULT_PARAMS["phiy"]
    coeffs, degs = model.get_characteristic_polynomial_data()
    char_poly = CharPoly(coeffs, degs)
    zm = ZeroManager(char_poly, E_ref, 0)
    zm.run(verbose=True)

    print("Multiple roots: ", zm.n_multiple_roots)
    print("Segments: ", zm.n_segments)
    colors = ["r", "g", "b", "m"]
    for j, seg in enumerate(zm.segments):
        theta1_arr = seg.theta1_arr
        tracked_roots = seg.tracked_roots
        plt.plot(theta1_arr, np.log(np.abs(tracked_roots)), '.-', color=colors[j])
    
    for j, r in enumerate(zm.multiple_roots):
        print(r)
        for rj in r.roots:
            plt.plot(r.theta1, np.log(np.abs(rj)), 'x')

    plt.show()


def solve_amoeba_E():
    E_ref = -1
    tx, ty, phix, phiy = DEFAULT_PARAMS["tx"], DEFAULT_PARAMS["ty"], DEFAULT_PARAMS["phix"], DEFAULT_PARAMS["phiy"]
    model = get_model(**DEFAULT_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()
    res = bfa.collect_GBZ_subsets(coeffs, degs, E_ref, 0.0, debug_mode=True)
    print(res)
    for subset in res.subsets:
        if isinstance(subset, LineSubset):
            theta1 = subset.theta1_arr
            theta2 = np.angle(subset.beta2_arr)
            plt.plot(theta1, theta2, '.-')
            print(tx * np.cos(theta1 + phix) + ty * np.cos(theta2 + phiy))
    plt.show()

if __name__ == "__main__":
    # solve_zeros()
    solve_amoeba_E()
