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


def build_model(params=PARAMS):
    """BerryPy convention: [source, destination, amplitude, cell shift]."""

    t1, t2, phi, mass, gamma = params
    lattice = np.array([
        [-0.5, -0.5], 
        [-np.sqrt(3) / 2, np.sqrt(3) / 2]
    ])
    intracell = [
        [0, 0, mass], 
        [1, 1, -mass], 
        [1, 0, t1], 
        [0, 1, t1]
    ]
    intercell = [
        [1, 0, t1, (0, -1)], 
        [1, 0, t1, (1, 0)],
        [0, 1, t1, (0, 1)], 
        [0, 1, t1, (-1, 0)]
    ]
    for sublattice, sign in ((0, 1), (1, -1)):
        for shift in ((-1, 0), (0, -1), (1, 1)):
            intercell.append([sublattice, sublattice,
                              t2 * np.exp(1j * (gamma + sign * phi)), shift])
            intercell.append([sublattice, sublattice,
                              t2 * np.exp(1j * (gamma - sign * phi)),
                              tuple(-np.array(shift))])
    model = tb.TightBindingModel(2, 2, lattice, intracell, intercell)
    sites_cart = np.array([
        [0, 1 / (2 * np.sqrt(3))], 
        [0, -1 / (2 * np.sqrt(3))]
    ])
    model.SiteCoord = model.cart2lattice(sites_cart.T).T
    return model


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


def _sweep_fun(enum_data):
    idx, data = enum_data
    if data[0] == "sgbz":
        return idx, sgbz.collect_GBZ_subsets(*data[1:])
    else:
        return idx, amoeba.collect_GBZ_subsets(*data[1:])


def sweep_GBZ(
    E_list: np.ndarray[complex, 1],
    model: tb.TightBindingModel,
    n_processes: int = 1,
    which: Literal["amoeba", "sgbz"] = "amoeba"
):
    from multiprocessing import Pool
    coeffs, degs = model.get_characteristic_polynomial_data()
    data_pack = [(which, coeffs, degs, E_ref) for E_ref in E_list]
    results = [0] * len(data_pack)

    n_tasks = len(data_pack)
    done = 0
    start = time.time()
    if n_processes == 1:
        for enum_data in enumerate(data_pack):
            idx, r = _sweep_fun(enum_data)
            results[idx] = r
            # Process bar
            done += 1
            draw_bar(done, n_tasks, time.time() - start)

    else:
        with Pool(n_processes) as pool:
            for idx, r in pool.imap_unordered(_sweep_fun, enumerate(data_pack)):
                results[idx] = r
                # Process bar
                done += 1
                draw_bar(done, n_tasks, time.time() - start)

    sys.stderr.write("\n")
    return results, coeffs, degs


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
    results, coeffs, degs = sweep_GBZ(
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
                "coeffs": coeffs,
                "degs": degs,
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


##### OBC geometries #####
def generate_rhombus_cells(n1: int, n2: int) -> list[tuple[int, int]]:
    """Primitive-cell coordinates with edges along a1 and a2.

    Return n1*n2 cells, ordered by a1 index first, then a2 index.
    """
    if n1 < 1 or n2 < 1:
        raise ValueError("Both rhombus sizes must be positive")
    return [(i, j) for i in range(n1) for j in range(n2)]


def generate_rectangle_cells(nx: int, ny: int) -> list[tuple[int, int]]:
    """Return 2*nx*ny primitive cells forming a Cartesian rectangular cut.

    For build_model's lattice, (1, 1) points along -x and (-1, 1)
    along +y. Each rectangular block contains primitive cells (0, 0)
    and (1, 0); including both preserves the original boundary termination.
    All returned coordinates remain in the original a1-a2 basis.
    """
    if nx < 1 or ny < 1:
        raise ValueError("Both rectangle sizes must be positive")
    ny_2 = ny // 2
    is_odd = ny % 2
    bulk = [
        (i - j + offset, i + j)
        for i in range(nx) for j in range(ny_2) for offset in (0, 1)
    ]
    if is_odd:
        bulk += ([
            (i - ny_2 , i + ny_2)
            for i in range(nx)
        ])
    return bulk


def add_obc_geom_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--geom", choices=("rhombus", "rectangle"),
                        default="rhombus", type=str, help="Geometry types")
    parser.add_argument("--n1", type=int, default=8, help="Parameter n1 for geometry definition")
    parser.add_argument("--n2", type=int, default=8, help="Parameter n2 for geometry definition")
    parser.add_argument("--preview", action="store_true", help="Preview the geometry.")
    parser.add_argument("--outdir", type=Path, default=DATA, help="Directory receiving the pickled OBC geometry.")
    parser.add_argument("--params", type=_num_parser, nargs=5, default=PARAMS,
                        metavar=("t1", "t2", "phi", "M", "gamma"),
                        help="Model parameters, overriding param.json; same expression syntax as param.json.")


def plot_OBC_geom(geom: Literal["rhombus", "rectangle"], n1: int, n2: int):
    model = build_model()
    if geom == "rhombus":
        coords = np.array(generate_rhombus_cells(n1, n2))
    elif geom == "rectangle":
        coords = np.array(generate_rectangle_cells(n1, n2))
    else:
        raise ValueError(f"Unknown geometry type: {geom}")

    site_colors = [COLORS["a1"], COLORS["x"]]

    import matplotlib.pyplot as plt
    for site_idx in range(model.SiteNum):
        site_coords = coords + model.SiteCoord[site_idx, :]
        site_cart = model.lattice2cart(site_coords.T).T
        plt.plot(site_cart[:,0], site_cart[:, 1], '.', color=site_colors[site_idx])
    plt.axis("equal")
    plt.show()


def calculate_OBC_spectrum(
    geom: Literal["rhombus", "rectangle"],
    n1: int,
    n2: int,
    outdir: Path = DATA,
    params: tuple[float] = PARAMS
) -> None:
    ''' Calculate OBC spectrum using scipy.linalg.eig '''
    from scipy import linalg as la
    from datetime import datetime, timezone

    model = build_model(params)
    if geom == "rhombus":
        coords = generate_rhombus_cells(n1, n2)
    elif geom == "rectangle":
        coords = generate_rectangle_cells(n1, n2)
    else:
        raise ValueError(f"Unknown geometry type: {geom}")
    H_mat = model.generate_OBC_bulk(coords).todense()
    eigv, eigvec = la.eig(H_mat)
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    data = {
        "geom": geom,
        "n1": n1,
        "n2": n2,
        "cell_coords": model.lattice2cart(np.array(coords).T).T,
        "site_coords": model.lattice2cart(model.SiteCoord.T).T,
        "eigv": eigv,
        "eigvec": eigvec,
        "run_id": run_id,
    }

    # save
    fname = outdir / f"obc_{geom}_{n1}_{n2}_{run_id}.pkl"
    with open(fname, "wb") as fp:
        pickle.dump(data, fp)
    print(f"Save {fname}")


def main_OBC_geom(kargs: argparse.Namespace):
    if kargs.preview:
        plot_OBC_geom(kargs.geom, kargs.n1, kargs.n2)
    else:
        calculate_OBC_spectrum(kargs.geom, kargs.n1, kargs.n2, kargs.outdir, kargs.params)


##### Post processing #####
def plot_GBZ_spectrum(fname):
    import matplotlib.pyplot as plt
    with open(fname, "rb") as fp:
        data = pickle.load(fp)
    #### Check E sequence ####
    E_re = np.linspace(data["E_re_range"][0], data["E_re_range"][1], data["n_re"])
    E_im = np.linspace(data["E_im_range"][0], data["E_im_range"][1], data["n_im"])
    E_re_mesh, E_im_mesh = np.meshgrid(E_re, E_im)
    E_mesh = E_re_mesh + 1j * E_im_mesh
    E_list = E_mesh.flatten(order=data["order"])
    E_list_res = np.array([r.E_ref for r in data["results"]])
    print("max|E_list - results.E_ref|", max(np.abs(E_list - E_list_res)))

    results = data["results"]
    E_gbz = np.array([r.E_ref for r in results if r.is_gbz])
    E_failed = np.array([r.E_ref for r in results if not r.success])
    E_non_gbz = np.array([r.E_ref for r in results if not r.is_gbz and r.success])
    print("Failed:", E_failed)
    for E_ref in E_failed:
        model = get_transformed_model("y")
        coeffs, degs = model.get_characteristic_polynomial_data()
        print(sgbz.collect_GBZ_subsets(coeffs, degs, E_ref, debug_mode=True))

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


def plot_OBC_spectrum(fname):
    with open(fname, "rb") as fp:
        data = pickle.load(fp)
    import matplotlib.pyplot as plt
    plt.plot(data["eigv"].real, data["eigv"].imag, '.')
    plt.show()


def plot_OBC_amplitude(fname, E_re_range, E_im_range):
    with open(fname, "rb") as fp:
        data = pickle.load(fp)

    # Assemble coordinates
    site_coords = data["site_coords"]
    n_sites = site_coords.shape[0]
    cell_coords = data["cell_coords"]
    all_coords = np.column_stack([cell_coords, cell_coords])
    for site_idx in range(n_sites):
        all_coords[:, site_idx * n_sites:(site_idx + 1) * n_sites] += site_coords[site_idx, :]
    all_coords = all_coords.reshape(-1, site_coords.shape[1])

    # Select eigenstates
    if E_re_range[0] == E_re_range[1]:
        mask = np.argmin(np.abs(data["eigv"] - (E_re_range[0] + 1j * E_im_range[0])))
        mean_amplitude = np.abs(data["eigvec"][:, mask])
    else:
        re_mask = (data["eigv"].real >= E_re_range[0]) & (data["eigv"].real <= E_re_range[1])
        im_mask = (data["eigv"].imag >= E_im_range[0]) & (data["eigv"].imag <= E_im_range[1])
        mask = re_mask & im_mask
        eigvec = data["eigvec"][:, mask]
        mean_amplitude = np.abs(eigvec).mean(axis=1)

    # Plot 
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(data["eigv"].real, data["eigv"].imag, '.', color="#e0e0e0")
    plt.plot(data["eigv"][mask].real, data["eigv"][mask].imag, '.', color=COLORS["x"])
    plt.figure()
    plt.scatter(all_coords[:, 0], all_coords[:, 1], c=mean_amplitude, cmap="viridis")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.show()


def main_plot(kargs: argparse.Namespace):
    if kargs.type == "gbz-spectrum":
        plot_GBZ_spectrum(kargs.fname)
    elif kargs.type == "obc-spectrum":
        plot_OBC_spectrum(kargs.fname)
    elif kargs.type == "obc-mean-amplitude":
        plot_OBC_amplitude(kargs.fname, kargs.E_re_range, kargs.E_im_range)
    else:
        raise ValueError(f"Unknown plot type {kargs.type}")


def add_plot_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--type", 
        choices=("gbz-spectrum", "obc-spectrum", "obc-mean-amplitude"), 
        default="gbz-spectrum",
        help="Type of plot to make."
    )
    parser.add_argument(
        "--fname", type=Path, required=True,
        help="Pickled sweep file to plot."
    )
    parser.add_argument(
        "--E_re_range", type=float, nargs=2, default=(-1e5, 1e5),
        help="Range of eigv.real to plot. If not specified, use the full range. If the lower and upper limits are set equal, use the nearest eigenvalue to (E_re_range[0], E_im_range[0]) instead."
    )
    parser.add_argument(
        "--E_im_range", type=float, nargs=2, default=(-1e5, 1e5),
        help="Range of eigv.imag to plot. If not specified, use the full range."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True, metavar="mode")
    add_sweep_arguments(
        modes.add_parser("sweep", help="Sweep the GBZ over a rectangular E_ref grid.")
    )
    add_plot_arguments(
        modes.add_parser("plot", help="Visualize the data.")
    )
    add_obc_geom_arguments(
        modes.add_parser("obc-geom", help="Build and visualize the OBC geometries.")
    )
    kargs = parser.parse_args()
    if kargs.mode == "sweep":
        main_GBZ_sweep(kargs)
    elif kargs.mode == "plot":
        main_plot(kargs)
    elif kargs.mode == "obc-geom":
        main_OBC_geom(kargs)
    else:
        raise ValueError(f"Unknown mode {kargs.mode}")
