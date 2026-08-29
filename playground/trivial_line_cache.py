'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-26
Copyright © Department of Physics, Tsinghua University. All rights reserved

Shared cache for the trivial_model real-axis scan used by the two
triangulation demos (advancing-front and edge-collapse).

Pickles a list of GBZResults along a real E grid so the two demos don't each
re-run collect_GBZ_subsets (~1.5 s per E point).  Run once to (re)build:

    python playground/trivial_line_cache.py
'''

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pickle
import numpy as np

import pygbz2d.amoeba as bfa
from trivial_model import get_model, DEFAULT_PARAMS

CACHE_PATH = Path(__file__).resolve().parent / "trivial_line_cache.pkl"


def build(E_grid=None, force=False):
    if CACHE_PATH.exists() and not force:
        return load()

    if E_grid is None:
        E_grid = np.arange(-1.9, 1.91, 0.1)

    model = get_model(**DEFAULT_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()

    results = []
    for E in E_grid:
        res = bfa.collect_GBZ_subsets(coeffs, degs, E, 0.0, debug_mode=False)
        results.append(res)

    data = {"E_grid": np.asarray(E_grid, dtype=float), "results": results}
    with open(CACHE_PATH, "wb") as fp:
        pickle.dump(data, fp)
    return data


def load():
    with open(CACHE_PATH, "rb") as fp:
        return pickle.load(fp)


if __name__ == "__main__":
    d = build(force=True)
    print(f"cached {len(d['results'])} results to {CACHE_PATH}")
    n_valid = sum(1 for r in d["results"] if r.success and r.index == (0, 2))
    print(f"  valid (index (0,2)) slices: {n_valid}")
