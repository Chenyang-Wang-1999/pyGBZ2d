'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-27
Copyright © Department of Physics, Tsinghua University. All rights reserved

Sweep + band-clustering test for the Hermitian Haldane model of
``Haldane-model-gainloss.py:Hermitian_Haldane_sweep``.

Physics target: the two bands' E-plane projections OVERLAP (no indirect
gap) while the direct gap never closes.  The radius-graph clustering must
still yield exactly 2 clusters there — separated by theta/mu blocks, NOT
by the E block (which sees both bands at the same E in the overlap).

Sweep spec (as agreed): band bounds from a BZ diagonalization, extended by
0.1 on both sides, 201 real E points via linspace, amoeba
collect_GBZ_subsets.  Output layout matches the other playground sweeps.

Usage:
    python playground/sweep_Hermitian_Haldane.py            # sweep + save
    python playground/demo_band_clustering.py data/Hermitian-Haldane-amoeba.pkl
'''

import sys
import importlib.util
import multiprocessing as mp
from math import pi
from pathlib import Path

import numpy as np
from scipy import linalg as la

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

_spec = importlib.util.spec_from_file_location(
    "haldane_gl", Path(__file__).resolve().parent / "Haldane-model-gainloss.py")
haldane_gl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(haldane_gl)

OUT_FNAME = "data/Hermitian-Haldane-amoeba.pkl"
EXTEND = 0.1
N_POINTS = 201
METHOD = "amoeba"          # or "sgbz" (a1)


def build_model():
    """Exact copy of Hermitian_Haldane_sweep's model construction."""
    model = haldane_gl.Haldane_non_Hermitian_phase(1, 0.5, pi / 3, 3, 0)
    model.InterCell += [
        (0, 0, 2, (1, 0)),
        (0, 0, 2, (-1, 0)),
        (1, 1, 2, (1, 0)),
        (1, 1, 2, (-1, 0)),
    ]
    return model


def bz_band_bounds(model, n_k: int = 101):
    """Band ranges + minimal direct gap from a BZ diagonalization."""
    ks = np.linspace(-pi, pi, n_k)
    k1, k2 = np.meshgrid(ks, ks)
    bands = np.empty((2, k1.size))
    for i in range(k1.size):
        H = model.get_bulk_Hamiltonian((k1.flat[i], k2.flat[i])).todense()
        bands[:, i] = np.sort(la.eigvalsh(H))
    lower, upper = bands[0], bands[1]
    return {
        "lower": (float(lower.min()), float(lower.max())),
        "upper": (float(upper.min()), float(upper.max())),
        "direct_gap_min": float((upper - lower).min()),
        "E_min": float(bands.min()),
        "E_max": float(bands.max()),
    }


def main():
    model = build_model()
    print("building characteristic polynomial (sympy, may take a moment)...")
    coeffs, degs = model.get_characteristic_polynomial_data()
    print(f"n_terms = {len(coeffs)}")

    bounds = bz_band_bounds(model)
    print(f"lower band E in [{bounds['lower'][0]:.4f}, {bounds['lower'][1]:.4f}]")
    print(f"upper band E in [{bounds['upper'][0]:.4f}, {bounds['upper'][1]:.4f}]")
    print(f"projection overlap: "
          f"{bounds['lower'][1] > bounds['upper'][0]} "
          f"(lower.max - upper.min = {bounds['lower'][1] - bounds['upper'][0]:+.4f})")
    print(f"minimal direct gap = {bounds['direct_gap_min']:.4f}")

    E_list = np.linspace(bounds["E_min"] - EXTEND,
                         bounds["E_max"] + EXTEND, N_POINTS)
    print(f"sweep: {len(E_list)} points on the real axis "
          f"[{E_list[0]:.4f}, {E_list[-1]:.4f}], method={METHOD}")

    if METHOD == "amoeba":
        import pygbz2d.amoeba as solver
    else:
        import pygbz2d.sgbz as solver

    # Serial loop: the DSH sandbox forbids the named pipes behind mp.Pool,
    # and ~0.4 s/solve x 201 points is acceptable anyway.
    results = []
    for j, E in enumerate(E_list):
        results.append(solver.collect_GBZ_subsets(coeffs, degs, E + 0j,
                                                  j / len(E_list)))
        if (j + 1) % 25 == 0:
            print(f"  {j + 1}/{len(E_list)} done")

    n_gbz = sum(1 for r in results if r.success and r.is_gbz)
    idxs = {}
    for r in results:
        if r.success and r.is_gbz:
            idxs[r.index] = idxs.get(r.index, 0) + 1
    print(f"in-GBZ: {n_gbz}/{len(results)}, indices={idxs}")

    with open(OUT_FNAME, "wb") as fp:
        import pickle
        pickle.dump({
            "E_real": E_list,
            "E_imag": np.array([0.0]),
            "results": results,
            "bounds": bounds,
            "coeffs": coeffs,
            "degs": degs,
        }, fp)
    print(f"saved -> {OUT_FNAME}")


if __name__ == "__main__":
    main()
