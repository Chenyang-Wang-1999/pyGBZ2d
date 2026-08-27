'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-06-26
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

import numpy as np
import BerryPy.TightBinding as tb
from cmath import exp, log
from math import sin, cos, sqrt, pi
from scipy import linalg as la

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import bfgbz2d.amoeba as bfa
print(bfa.__file__)


def get_model(gamma, u1, u2, v1, v2, u3, v3, w):
    dim = 2
    site_num = 1
    inter_cell = [
        [0, 0, u1 * exp(gamma), (1, 0)],
        [0, 0, u1 * exp(-gamma), (-1, 0)],
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
    return tb.TightBindingModel(
        dim, site_num, np.eye(2), [], inter_cell, [(0, 0)]
    )


def get_parametric_1D_model(
    ky, gamma, u1, u2, v1, v2, u3, v3, w
):
    dim = 1
    site_num = 1
    intra_cell = [
        [0, 0, 2 * v1 * np.cos(ky) + 2 * v2 * np.cos(2 * ky) + 2 * v3 * np.cos(3 * ky)]
    ]
    inter_cell = [
        [0, 0, u1 * exp(gamma), (1,)],
        [0, 0, u1 * exp(-gamma), (-1,)],
        [0, 0, u2, (2,)],
        [0, 0, u2, (-2,)],
        [0, 0, u3, (3,)],
        [0, 0, u3, (-3,)],
        [0, 0, 2 * w * np.cos(ky), (-1,)],
        [0, 0, 2 * w * np.cos(ky), (1,)],
    ]
    return tb.TightBindingModel(
        dim, site_num, [[1]], intra_cell, inter_cell, [(0, 0)]
    )


DEFAULT_PARAMS = {
    "gamma": 0.2,
    "u1": 1,
    "v1": 0.8,
    "w": 0.5,
    "u2": 0.1,
    "v2": 0.1,
    "u3": 0.05,
    "v3": 0.05,
}

MODEL_NAME = "next-nearest-coupling"


def debug_amoeba(E_ref: complex):
    model = get_model(**DEFAULT_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()
    print(bfa.collect_GBZ_subsets(coeffs, degs, E_ref, 0.0, True))

if __name__ == "__main__":
    debug_amoeba(5.916 + 0.1344j)
