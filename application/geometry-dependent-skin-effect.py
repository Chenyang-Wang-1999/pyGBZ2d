"""Gain-loss Haldane model: directional SGBZs and finite open geometries.

The first positional argument selects a mode, e.g. "sweep"; each mode then takes
its own keyword arguments, and "--help" after the mode lists them, while "--help"
alone lists the modes. OBC, plotting, validation and legacy-data extraction modes
are still to be added.
All paths default to this checkout, so the working directory is immaterial.
BerryPy and SymPy are needed for modeling; plotting additionally uses Matplotlib.
Only load pickle files from trusted research data.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import pickle
import sys
import time

import numpy as np
from math import sin, cos, sqrt, pi
from cmath import log, exp
from typing import Literal
from BerryPy import TightBinding as tb

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pygbz2d import sgbz
from pygbz2d import amoeba

DATA = ROOT / "application/data/geometry-dependent-skin-effect"
FIGURES = ROOT / "application/Figures/geometry-dependent-skin-effect"

#### load parameters from json ####
with open(FIGURES / "param.json", "r", encoding="utf-8") as fp:
    param_txt = json.load(fp)
def _num_parser(num_str: str):
    return eval(num_str.replace("ii", "j").replace("ee", "exp(1)"))
PARAMS = (
    _num_parser(param_txt["t1"]),
    _num_parser(param_txt["t2"]),
    _num_parser(param_txt["phi"]),
    _num_parser(param_txt["M"]),
    _num_parser(param_txt["gamma"])
)
###################################

DIRECTIONS = ("a1", "a2", "x", "y")
COLORS = {"a1": "#305f9e", "a2": "#269c95", "x": "#d29229", "y": "#c24752"}


# def build_model(params=PARAMS):
#     """BerryPy convention: [source, destination, amplitude, cell shift]."""

#     t1, t2, phi, mass, gamma = params
#     lattice = np.array([
#         [-0.5, -0.5], 
#         [-np.sqrt(3) / 2, np.sqrt(3) / 2]
#     ])
#     intracell = [
#         [0, 0, mass], 
#         [1, 1, -mass], 
#         [1, 0, t1], 
#         [0, 1, t1]
#     ]
#     intercell = [
#         [1, 0, t1, (0, -1)], 
#         [1, 0, t1, (1, 0)],
#         [0, 1, t1, (0, 1)], 
#         [0, 1, t1, (-1, 0)]
#     ]
#     for sublattice, sign in ((0, 1), (1, -1)):
#         for shift in ((-1, 0), (0, -1), (1, 1)):
#             intercell.append([sublattice, sublattice,
#                               t2 * np.exp(1j * (gamma + sign * phi)), shift])
#             intercell.append([sublattice, sublattice,
#                               t2 * np.exp(1j * (gamma - sign * phi)),
#                               tuple(-np.array(shift))])
#     model = tb.TightBindingModel(2, 2, lattice, intracell, intercell)
#     sites_cart = np.array([
#         [0, 1 / (2 * np.sqrt(3))], 
#         [0, -1 / (2 * np.sqrt(3))]
#     ])
#     model.SiteCoord = model.cart2lattice(sites_cart.T).T
#     return model

def non_Hermitian_Haldane_H(u1, u2, v1, v2, phi, M):
    dim = 2
    site_num = 2

    lattice_vec = np.array(
        [[-cos(pi/3), -cos(pi/3)],
         [-sin(pi/3), sin(pi/3)]]
    )

    intra_cell = [
        [0, 0, M],
        [1, 1, -M],
        [1, 0, u1],
        [0, 1, u2]
    ]
    inter_cell = [
        [0, 0, v2 * exp(1j*phi), (-1,0)],
        [0, 0, v2 * exp(1j*phi), (0, -1)],
        [0, 0, v2 * exp(1j*phi), (1,1)],
        [0, 0, v1 * exp(-1j*phi), (1,0)],
        [0, 0, v1 * exp(-1j*phi), (0,1)],
        [0, 0, v1 * exp(-1j*phi), (-1,-1)],
        [1, 0, u1, (0, -1)],
        [1, 0, u1, (1, 0)],
        [0, 1, u2, (0, 1)],
        [0, 1, u2, (-1, 0)],
        [1, 1, v2 * exp(-1j * phi), (-1, 0)],
        [1, 1, v2 * exp(-1j*phi), (0, -1)],
        [1, 1, v2 * exp(-1j*phi), (1, 1)],
        [1, 1, v1 * exp(1j * phi), (1, 0)],
        [1, 1, v1 * exp(1j * phi), (0, 1)],
        [1, 1, v1 * exp(1j * phi), (-1,-1)]
    ]

    site_coord_cart = np.array(
        [[0, 1 / (2 * sqrt(3))],
         [0, - 1 / (2 * sqrt(3))]]
    )

    model = tb.TightBindingModel(dim, site_num, lattice_vec, intra_cell, inter_cell)
    model.SiteCoord = model.cart2lattice(site_coord_cart.T).T

    return model


def Haldane_non_Hermitian_phase(t1, t2, phi, M, gamma):
    # v1 = v2 = t2 * exp(i gamma)
    return non_Hermitian_Haldane_H(t1, t1, t2 * exp(1j * gamma), t2 * exp(1j * gamma), phi, M)


def build_model(params=PARAMS):
    return Haldane_non_Hermitian_phase(*params)


def get_transformed_model(which: Literal["a1", "a2", "x", "y"], params=PARAMS):
    model = build_model(params)
    if which == "a1":
        return model
    elif which == "a2":
        # Exchange the two axes
        return model.get_supercell(
            [(0, 0)],
            np.array(
                [[0, -1],
                 [1, 0]],
                dtype=int
            )
        )
    elif which == "x":
        # Supercell with cells (0, 0) and (1, 0)
        return model.get_supercell(
            [(0, 0), (1, 0)],
            np.array(
                [[1, -1],
                 [1, 1]],
                dtype=int
            )
        )
    elif which == "y":
        return model.get_supercell(
            [(0, 0), (1, 0)],
            np.array(
                [[-1, -1],
                 [1, -1]],
                dtype=int
            )
        )
    else:
        raise ValueError(f"Unknown type {which}")


def bz_spectrum(model, nk=101):
    phases = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    # get_bulk_Hamiltonian takes fractional reciprocal coordinates, not radians.
    return np.concatenate([
        np.linalg.eigvals(model.get_bulk_Hamiltonian_complex(np.exp(1j * np.array([a, b]))).toarray())
        for a in phases for b in phases
    ])


##### SGBZ sweep #####
def draw_bar(done, total, elapsed, width=40, prefix="进度"):
    """ Draw process bar in terminal while sweeping """
    ratio = done / total if total else 1.0
    filled = int(width * ratio)
    bar = "█" * filled + "-" * (width - filled)

    speed = done / elapsed if elapsed > 0 else 0.0          # 每秒完成数
    eta = (total - done) / speed if speed > 0 else 0.0      # 预计剩余秒数

    line = (f"\r{prefix}: |{bar}| {done}/{total} "
            f"({ratio:5.1%}) {speed:6.1f}it/s ETA {eta:5.1f}s")
    # 末尾补空格，防止残留上一次更长的内容
    sys.stderr.write(line.ljust(100))
    sys.stderr.flush()


def _sweep_fun(data):
    if data[0] == "sgbz":
        return sgbz.collect_GBZ_subsets(*data[1:])
    else:
        return amoeba.collect_GBZ_subsets(*data[1:])


def sweep_GBZ(
    E_list: np.ndarray[complex, 1],
    model: tb.TightBindingModel,
    n_processes: int = 1,
    which: Literal["amoeba", "sgbz"] = "amoeba"
):
    from multiprocessing import Pool
    coeffs, degs = model.get_characteristic_polynomial_data()
    data_pack = [(which, coeffs, degs, E_ref) for E_ref in E_list]
    results = []

    n_tasks = len(data_pack)
    done = 0
    start = time.time()
    if n_processes == 1:
        for data in data_pack:
            results.append(_sweep_fun(data))
            # Process bar
            done += 1
            draw_bar(done, n_tasks, time.time() - start)

    else:
        with Pool(n_processes) as pool:
            for r in pool.imap_unordered(_sweep_fun, data_pack):
                results.append(r)
                # Process bar
                done += 1
                draw_bar(done, n_tasks, time.time() - start)

    sys.stderr.write("\n")
    return results


def sweep_structured_grid(
    E_re_range: tuple[float, float],
    E_im_range: tuple[float, float],
    n_re: int,
    n_im: int,
    n_processes: int = 1,
    params: tuple = PARAMS,
    which: Literal["amoeba", "x", "y", "a1", "a2"] = "amoeba",
    outdir: str = DATA
):
    ''' Sweep GBZ over a grid '''
    from datetime import datetime, timezone
    os.makedirs(outdir, exist_ok=True)

    #### Build E mesh ####
    E_re = np.linspace(E_re_range[0], E_re_range[1], n_re)
    E_im = np.linspace(E_im_range[0], E_im_range[1], n_im)
    E_re_mesh, E_im_mesh = np.meshgrid(E_re, E_im)
    E_mesh = E_re_mesh + 1j * E_im_mesh
    E_list = E_mesh.flatten(order="C")

    #### Build model ####
    if which == "amoeba":
        model = get_transformed_model("a1", params)
    else:
        model = get_transformed_model(which, params)

    #### Sweep ####
    results = sweep_GBZ(
        E_list,
        model,
        n_processes,
        "amoeba" if which == "amoeba" else "sgbz"
    )

    #### Save data ####
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    fname = outdir / f"sweep_{which}_{n_re}_{n_im}_{run_id}.pkl"
    with open(fname, "wb") as fp:
        pickle.dump(
            {
                "run_id": run_id,
                "E_re_range": E_re_range,
                "E_im_range": E_im_range,
                "which": which,
                "n_re": n_re,
                "n_im": n_im,
                "order": "C",
                "params": params,
                "results": results
            },
            fp
        )
    print(f"Data saved to {fname}")


def add_sweep_arguments(parser: argparse.ArgumentParser) -> None:
    """Register the keyword arguments of the `sweep` mode.

    Each flag is spelled exactly like the sweep_structured_grid parameter it
    feeds, so main_GBZ_sweep can forward vars(kargs) without a translation table.
    """
    parser.add_argument(
        "--E_re_range", type=float, nargs=2, required=True, metavar=("MIN", "MAX"),
        help="Real part of the E_ref grid, given as MIN MAX."
    )
    parser.add_argument(
        "--E_im_range", type=float, nargs=2, required=True, metavar=("MIN", "MAX"),
        help="Imaginary part of the E_ref grid, given as MIN MAX."
    )
    parser.add_argument(
        "--n_re", type=int, required=True,
        help="Number of E_ref samples along the real axis."
    )
    parser.add_argument(
        "--n_im", type=int, required=True,
        help="Number of E_ref samples along the imaginary axis."
    )
    parser.add_argument(
        "--n_processes", type=int, default=1,
        help="Worker processes used by the sweep."
    )
    parser.add_argument(
        "--params", type=_num_parser, nargs=5, default=PARAMS,
        metavar=("t1", "t2", "phi", "M", "gamma"),
        help="Model parameters, overriding param.json; same expression syntax as param.json."
    )
    parser.add_argument(
        "--which", choices=("amoeba", "a1", "a2", "x", "y"), default="amoeba",
        help="Direction and solver used by the sweep."
    )
    parser.add_argument(
        "--outdir", type=Path, default=DATA,
        help="Directory receiving the pickled sweep."
    )


def main_GBZ_sweep(kargs: argparse.Namespace):
    # `mode` only routes the dispatch; every other entry is a sweep parameter.
    kwargs = {key: value for key, value in vars(kargs).items() if key != "mode"}
    sweep_structured_grid(**kwargs)


##### Post processing #####
def plot_GBZ_spectrum(fname):
    import matplotlib.pyplot as plt
    with open(fname, "rb") as fp:
        data = pickle.load(fp)
    results = data["results"]
    E_gbz = np.array([r.E_ref for r in results if r.is_gbz])
    E_failed = np.array([r.E_ref for r in results if not r.success])
    E_non_gbz = np.array([r.E_ref for r in results if not r.is_gbz and r.success])

    ### Print E range ###
    E_re = np.linspace(data["E_re_range"][0], data["E_re_range"][1], data["n_re"])
    E_im = np.linspace(data["E_im_range"][0], data["E_im_range"][1], data["n_im"])
    new_box = (
        max(E_re[E_re < E_gbz.real.min()]), 
        min(E_re[E_re > E_gbz.real.max()]),
        max(E_im[E_im < E_gbz.imag.min()]),
        min(E_im[E_im > E_gbz.imag.max()]),
    )
    print(f"GBZ E range: {new_box[0]} < Re < {new_box[1]}, {new_box[2]} < Im < {new_box[3]}")

    ### Plot ###
    plt.plot(E_gbz.real, E_gbz.imag, '.', label="GBZ")
    plt.plot(E_non_gbz.real, E_non_gbz.imag, '.', label="Non-GBZ")
    plt.plot(E_failed.real, E_failed.imag, '.', label="Failed")
    plt.plot(
        [new_box[0], new_box[1], new_box[1], new_box[0], new_box[0]],
        [new_box[2], new_box[2], new_box[3], new_box[3], new_box[2]],
        'r--'
    )
    plt.legend()
    plt.show()


def main_plot(kargs: argparse.Namespace):
    if kargs.type == "gbz-spectrum":
        plot_GBZ_spectrum(kargs.fname)
    else:
        raise ValueError(f"Unknown plot type {kargs.type}")


def add_plot_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--type", choices=("gbz-spectrum",), default="gbz-spectrum",
        help="Type of plot to make."
    )
    parser.add_argument(
        "--fname", type=Path, required=True,
        help="Pickled sweep file to plot."
    )


def finite_hamiltonian(model, nx, ny):
    """Assemble BerryPy's hopping list in O(number of bonds), discarding exits.

    This avoids constructing a huge periodic supercell just to delete its wrap
    bonds. check_model() compares it with BerryPy's supercell OBC Hamiltonian.
    """
    from scipy.sparse import coo_matrix

    if nx < 1 or ny < 1:
        raise ValueError("Both OBC sizes must be positive")
    cells = [(i, j) for i in range(nx) for j in range(ny)]
    lookup = {cell: k for k, cell in enumerate(cells)}
    n = model.SiteNum
    rows, cols, values = [], [], []
    bonds = [(*b, (0, 0)) for b in model.InCell] + list(model.InterCell)
    for index, (i, j) in enumerate(cells):
        for source, dest, amplitude, shift in bonds:
            target = lookup.get((i + shift[0], j + shift[1]))
            if target is not None:
                rows.append(target * n + dest)
                cols.append(index * n + source)
                values.append(amplitude)
    coords = np.vstack([(np.asarray(model.SiteCoord) + cell) @ model.LatticeVec.T for cell in cells])
    return coo_matrix((values, (rows, cols)), shape=(len(coords), len(coords))).tocsr(), coords


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True, metavar="mode")
    add_sweep_arguments(
        modes.add_parser("sweep", help="Sweep the GBZ over a rectangular E_ref grid.")
    )
    add_plot_arguments(
        modes.add_parser("plot", help="Visualize the data.")
    )
    kargs = parser.parse_args()
    if kargs.mode == "sweep":
        main_GBZ_sweep(kargs)
    elif kargs.mode == "plot":
        main_plot(kargs)
    else:
        raise ValueError(f"Unknown mode {kargs.mode}")
