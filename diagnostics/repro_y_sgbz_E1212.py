'''Reproduce the y-SGBZ failure at E_ref = 1.212 recorded in
log/2026-08-16-pairwise-linear-intersection.md ("当前仍存在的问题").

Single energy point from sweep_SGBZ_y's grid: E_real index of 1.212 in
linspace(-3.1, 4.6, 51) does not land exactly on 1.212, so use the same
construction as the log: E_ref = 1.212 (the log reports 1.2120000000000002,
i.e. linspace round-off).  Try the exact grid point nearest 1.212 as well.
'''
import sys
from pathlib import Path
from math import pi

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import importlib.util

import numpy as np
import bfgbz2d.sgbz as bfs

_spec = importlib.util.spec_from_file_location(
    "haldane_gainloss", ROOT / "demos" / "Haldane-model-gainloss.py")
_hg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_hg)
Haldane_non_Hermitian_phase, ALL_PARAMS = _hg.Haldane_non_Hermitian_phase, _hg.ALL_PARAMS

model = Haldane_non_Hermitian_phase(*ALL_PARAMS)
model = model.get_supercell(
    [(0, 0), (1, 0)],
    np.array([[1, 1], [-1, 1]], dtype=int),
)
coeffs, degs = model.get_characteristic_polynomial_data()

# Exact grid point nearest 1.212 (sweep grid: linspace(-3.1, 4.6, 51))
grid = np.linspace(-3.1, 4.6, 51)
E_ref = complex(grid[np.argmin(np.abs(grid - 1.212))])
print(f"E_ref = {E_ref!r}")

res = bfs.collect_GBZ_subsets(coeffs, degs, E_ref, debug_mode=True)
print(res)
