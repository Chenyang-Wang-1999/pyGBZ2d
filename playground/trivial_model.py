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
from cmath import exp, pi

import numpy as np


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


