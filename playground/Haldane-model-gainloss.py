'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-01-15
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

import numpy as np
from math import sin, cos, sqrt, pi
from cmath import exp, log
import numpy as np
from BerryPy import TightBinding as tb
import matplotlib.pyplot as plt
import pickle
from scipy import linalg as la
from typing import Literal, Optional
import multiprocessing as mp
import os
from scipy import interpolate
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]) + "/src")
import pygbz2d.sgbz as bfs
import pygbz2d.amoeba as bfa
from pygbz2d.backend import make_laurent
print("Using ", bfs.__file__)
print("Polynomial backend:", type(make_laurent(np.array([1, 1]), np.array([1, 0, 0, 0, 1, 0]))))

ALL_PARAMS = (
    1, #t1,
    0.5, #t2,
    pi / 3, # phi
    0.5j, # M
    0 # gamma
)


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


def plot_BZ_bands():
    t1, t2, phi, M, gamma = ALL_PARAMS

    model = Haldane_non_Hermitian_phase(t1, t2, phi, M, gamma)
 
    k1 = np.linspace(-np.pi, np.pi, 100)
    k2 = np.linspace(-np.pi, np.pi, 100)
    k1_grid, k2_grid = np.meshgrid(k1, k2)
    E1_grid = np.zeros_like(k1_grid, dtype=complex)
    E2_grid = np.zeros_like(k1_grid, dtype=complex)
    for i in range(k1_grid.shape[0]):
        for j in range(k1_grid.shape[1]):
            H = model.get_bulk_Hamiltonian((k1_grid[i, j], k2_grid[i, j])).todense()

            E1_grid[i, j], E2_grid[i, j] = la.eig(H)[0]

    plt.plot(E1_grid.real, E1_grid.imag, 'bx')
    plt.plot(E2_grid.real, E2_grid.imag, 'bx')
    plt.axis("equal")


def sweep_amoeba():
    model = Haldane_non_Hermitian_phase(*ALL_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()

    N_points = 201
    E_real = np.linspace(-3.1, 4.6, N_points)
    E_imag = np.linspace(-0.51, 0.51, N_points)
    E_real_mesh, E_imag_mesh = np.meshgrid(E_real, E_imag)
    E_mesh = E_real_mesh + 1j * E_imag_mesh
    E_list = E_mesh.flatten()

    pool = mp.Pool(12)
    results = pool.starmap(bfa.collect_GBZ_subsets, [(coeffs, degs, E, j / len(E_list), True) for j, E in enumerate(E_list)])
    pool.close()
    pool.join()

    fname_prefix = "data/Haldane-gain-loss-amoeba"

    with open(fname_prefix + ".pkl", "wb") as fp:
        pickle.dump({
            "E_real": E_real,
            "E_imag": E_imag,
            "results": results,
            "params": ALL_PARAMS,
            "coeffs": coeffs,
            "degs": degs
        }, fp)


def sweep_amoeba_multiband():
    model = Haldane_non_Hermitian_phase(*ALL_PARAMS)
    model = model.get_supercell(
        [(0, 0), (1, 0)],
        np.array([
            [1, 1],
            [1, -1]
        ], dtype=int)
    )
 
    coeffs, degs = model.get_characteristic_polynomial_data()

    N_points = 201
    E_real = np.linspace(-3.1, 4.6, N_points)
    E_imag = np.linspace(-0.51, 0.51, N_points)
    E_real_mesh, E_imag_mesh = np.meshgrid(E_real, E_imag)
    E_mesh = E_real_mesh + 1j * E_imag_mesh
    E_list = E_mesh.flatten()

    pool = mp.Pool(12)
    results = pool.starmap(bfa.collect_GBZ_subsets, [(coeffs, degs, E, j / len(E_list), True) for j, E in enumerate(E_list)])
    pool.close()
    pool.join()

    fname_prefix = "data/Haldane-gain-loss-amoeba-xy"

    with open(fname_prefix + ".pkl", "wb") as fp:
        pickle.dump({
            "E_real": E_real,
            "E_imag": E_imag,
            "results": results,
            "params": ALL_PARAMS,
            "coeffs": coeffs,
            "degs": degs
        }, fp)



def sweep_SGBZ_a1():
    model = Haldane_non_Hermitian_phase(*ALL_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()

    N_points = 201
    E_real = np.linspace(-3.1, 4.6, N_points)
    E_imag = np.linspace(-0.51, 0.51, N_points)
    E_real_mesh, E_imag_mesh = np.meshgrid(E_real, E_imag)
    E_mesh = E_real_mesh + 1j * E_imag_mesh
    E_list = E_mesh.flatten()

    pool = mp.Pool(mp.cpu_count())
    results = pool.starmap(bfs.collect_GBZ_subsets, [(coeffs, degs, E, j / len(E_list)) for j, E in enumerate(E_list)])
    pool.close()
    pool.join()

    fname_prefix = "data/Haldane-gain-loss-a1-SGBZ"

    with open(fname_prefix + ".pkl", "wb") as fp:
        pickle.dump({
            "E_real": E_real,
            "E_imag": E_imag,
            "results": results,
            "params": ALL_PARAMS,
            "coeffs": coeffs,
            "degs": degs
        }, fp)


def sweep_SGBZ_a2():
    model = Haldane_non_Hermitian_phase(*ALL_PARAMS)
    model = model.get_supercell(
        [(0,0)],
        np.array([
            [0, 1],
            [1, 0]
        ], dtype=int)
    )
    coeffs, degs = model.get_characteristic_polynomial_data()

    N_points = 201
    E_real = np.linspace(-3.1, 4.6, N_points)
    E_imag = np.linspace(-0.51, 0.51, N_points)
    E_real_mesh, E_imag_mesh = np.meshgrid(E_real, E_imag)
    E_mesh = E_real_mesh + 1j * E_imag_mesh
    E_list = E_mesh.flatten()

    pool = mp.Pool(mp.cpu_count())
    results = pool.starmap(bfs.collect_GBZ_subsets, [(coeffs, degs, E, j / len(E_list)) for j, E in enumerate(E_list)])
    pool.close()
    pool.join()

    fname_prefix = "data/Haldane-gain-loss-a2-SGBZ"

    with open(fname_prefix + ".pkl", "wb") as fp:
        pickle.dump({
            "E_real": E_real,
            "E_imag": E_imag,
            "results": results,
            "params": ALL_PARAMS,
            "coeffs": coeffs,
            "degs": degs
        }, fp)


def sweep_SGBZ_x():
    model = Haldane_non_Hermitian_phase(*ALL_PARAMS)
    model = model.get_supercell(
        [(0, 0), (1, 0)],
        np.array([
            [1, 1],
            [1, -1]
        ], dtype=int)
    )
    coeffs, degs = model.get_characteristic_polynomial_data()

    N_points = 201
    E_real = np.linspace(-3.1, 4.6, N_points)
    E_imag = np.linspace(-0.51, 0.51, N_points)
    E_real_mesh, E_imag_mesh = np.meshgrid(E_real, E_imag)
    E_mesh = E_real_mesh + 1j * E_imag_mesh
    E_list = E_mesh.flatten()

    pool = mp.Pool(mp.cpu_count())
    results = pool.starmap(bfs.collect_GBZ_subsets, [(coeffs, degs, E, j / len(E_list)) for j, E in enumerate(E_list)])
    pool.close()
    pool.join()

    fname_prefix = "data/Haldane-gain-loss-x-SGBZ"

    with open(fname_prefix + ".pkl", "wb") as fp:
        pickle.dump({
            "E_real": E_real,
            "E_imag": E_imag,
            "results": results,
            "params": ALL_PARAMS,
            "coeffs": coeffs,
            "degs": degs
        }, fp)


def sweep_SGBZ_y():
    model = Haldane_non_Hermitian_phase(*ALL_PARAMS)
    model = model.get_supercell(
        [(0, 0), (1, 0)],
        np.array([
            [1, 1],
            [-1, 1]
        ], dtype=int)
    )
    coeffs, degs = model.get_characteristic_polynomial_data()

    N_points = 201
    E_real = np.linspace(-3.1, 4.6, N_points)
    E_imag = np.linspace(-0.51, 0.51, N_points)
    # E_imag = np.array([0])
    # E_real = np.linspace(-3.1, 4.6, 5 * N_points)
    E_real_mesh, E_imag_mesh = np.meshgrid(E_real, E_imag)
    E_mesh = E_real_mesh + 1j * E_imag_mesh
    E_list = E_mesh.flatten()

    pool = mp.Pool(mp.cpu_count())
    results = pool.starmap(bfs.collect_GBZ_subsets, [(coeffs, degs, E, j / len(E_list)) for j, E in enumerate(E_list)])
    pool.close()
    pool.join()

    fname_prefix = "data/Haldane-gain-loss-y-SGBZ"

    with open(fname_prefix + ".pkl", "wb") as fp:
        pickle.dump({
            "E_real": E_real,
            "E_imag": E_imag,
            "results": results,
            "params": ALL_PARAMS,
            "coeffs": coeffs,
            "degs": degs,
        }, fp)


def recompute_failed_SGBZ(
    which: str = "y",
    *,
    fname: Optional[os.PathLike] = None,
    out_fname: Optional[os.PathLike] = None,
    n_procs: Optional[int] = None,
    debug_mode: bool = False,
    failed_predicate=None,
):
    """并行重算已存储 SGBZ sweep 中的失败 E，并写回结果列表。

    从 ``data/Haldane-gain-loss-{which}-SGBZ.pkl`` 中按
    ``failed_predicate``（默认 ``lambda res: not res.success``，与
    :func:`plot_SGBZ` 的 Failed 判定一致）筛选出失败项，取出它们的
    ``E_ref``，再用与 :func:`sweep_SGBZ_y` 相同的
    ``mp.Pool.starmap(bfs.collect_GBZ_subsets, ...)`` 方式并行重算。
    pkl 中已存储的 ``coeffs`` / ``degs`` 被直接复用，因此无需重新构造
    BerryPy 模型，也保证重算对象与存储文件完全一致。

    Parameters
    ----------
    which :
        数据文件后缀，如 ``"y"``、``"x"``、``"a1"``、``"a2"``。
    fname :
        输入 pkl 路径；``None`` 时使用
        ``data/Haldane-gain-loss-{which}-SGBZ.pkl``。
    out_fname :
        输出 pkl 路径；``None`` 时覆盖 *fname*。只想试算时请显式指定
        一个其他路径。
    n_procs :
        并行进程数；``None`` 时使用 ``multiprocessing.cpu_count()``。
    debug_mode :
        传给 ``bfs.collect_GBZ_subsets`` 的 debug 开关。默认 ``False``，
        重算失败会返回 ``success=False`` 的 GBZResult 而不是抛异常。
    failed_predicate :
        失败判定函数 ``(GBZResult) -> bool``；默认 ``not res.success``。

    Returns
    -------
    data : dict
        更新后的 pkl 数据 dict（``results`` 中仅失败项被替换）。
    """
    if failed_predicate is None:
        failed_predicate = lambda res: (not res.is_gbz and not res.success)

    fname = Path(fname) if fname is not None else Path(
        "data/Haldane-gain-loss-%s-SGBZ.pkl" % which)
    if out_fname is None:
        out_fname = fname
    else:
        out_fname = Path(out_fname)

    with open(fname, "rb") as fp:
        data = pickle.load(fp)

    results = list(data["results"])
    coeffs = data["coeffs"]
    degs = data["degs"]

    failed_idx = [
        i for i, res in enumerate(results) if failed_predicate(res)
    ]
    if not failed_idx:
        print(f"{fname}: no failed entries ({len(results)} total)")
        return data

    failed_E = [results[i].E_ref for i in failed_idx]
    n_total = len(results)
    n_failed = len(failed_E)
    print(f"{fname}: recomputing {n_failed}/{n_total} failed entries "
          f"(which={which!r})")

    n_procs = n_procs or mp.cpu_count()
    tasks = [
        (coeffs, degs, E, j / n_failed, debug_mode)
        for j, E in enumerate(failed_E)
    ]

    with mp.Pool(n_procs) as pool:
        new_results = pool.starmap(bfs.collect_GBZ_subsets, tasks)

    if len(new_results) != n_failed:
        raise RuntimeError(
            f"parallel recompute returned {len(new_results)} results, "
            f"expected {n_failed}"
        )

    for idx, new_res in zip(failed_idx, new_results):
        results[idx] = new_res

    still_failed = sum(
        1 for i in failed_idx if failed_predicate(results[i])
    )
    data["results"] = results
    data["recompute_info"] = {
        "which": which,
        "failed_indices": failed_idx,
        "n_recomputed": n_failed,
        "n_still_failed": still_failed,
        "debug_mode": debug_mode,
    }

    with open(out_fname, "wb") as fp:
        pickle.dump(data, fp)

    print(f"saved -> {out_fname}; still failed: {still_failed}/{n_failed}")
    return data


def plot_amoebic_spectrum(suffix=""):
    with open("data/Haldane-gain-loss-amoeba%s.pkl" % (suffix), "rb") as fp:
        data = pickle.load(fp)
    E_real, E_imag, res, params = data["E_real"], data["E_imag"], data["results"], data["params"]
    E_real_mesh, E_imag_mesh = np.meshgrid(E_real, E_imag)
    E_list = (E_real_mesh + 1j * E_imag_mesh).flatten()
    coeffs, degs = data["coeffs"], data["degs"]

    ind_failed = [i for i in range(len(res)) if not res[i].is_gbz and not res[i].success]
    ind_amoeba = [i for i in range(len(res)) if res[i].is_gbz]
    ind_not_amoeba = [i for i in range(len(res)) if not res[i].is_gbz and res[i].success]

    print("Failed:", ind_failed)

    # Plot
    plt.figure()
    plt.plot(E_list[ind_amoeba].real, E_list[ind_amoeba].imag, '.', label="SGBZ")
    plt.plot(E_list[ind_not_amoeba].real, E_list[ind_not_amoeba].imag, '.', label="Non-SGBZ")
    plt.plot(E_list[ind_failed].real, E_list[ind_failed].imag, '.', label="Failed")
    plt.legend()



def plot_SGBZ(which="a1"):
    with open("data/Haldane-gain-loss-%s-SGBZ.pkl" % which, "rb") as fp:
        data = pickle.load(fp)
    E_real, E_imag, res, params = data["E_real"], data["E_imag"], data["results"], data["params"]
    E_real_mesh, E_imag_mesh = np.meshgrid(E_real, E_imag)
    E_list = (E_real_mesh + 1j * E_imag_mesh).flatten()
    coeffs, degs = data["coeffs"], data["degs"]

    ind_failed = [i for i in range(len(res)) if not res[i].is_gbz and not res[i].success]
    ind_amoeba = [i for i in range(len(res)) if res[i].is_gbz]
    ind_not_amoeba = [i for i in range(len(res)) if not res[i].is_gbz and res[i].success]

    print("Failed:", ind_failed)
    if len(ind_failed) > 0:
        print(res[ind_failed[0]].error)

    # Plot
    plt.figure()
    plt.plot(E_list[ind_amoeba].real, E_list[ind_amoeba].imag, '.', label="SGBZ")
    plt.plot(E_list[ind_not_amoeba].real, E_list[ind_not_amoeba].imag, '.', label="Non-SGBZ")
    plt.plot(E_list[ind_failed].real, E_list[ind_failed].imag, '.', label="Failed")
    plt.legend()
    plt.show()


def plot_amoeba_mu(suffix=""):
    with open("data/Haldane-gain-loss-amoeba%s.pkl" % (suffix), "rb") as fp:
        data = pickle.load(fp)
    E_real, E_imag, res, params = data["E_real"], data["E_imag"], data["results"], data["params"]
    E_real_mesh, E_imag_mesh = np.meshgrid(E_real, E_imag)
    E_list = (E_real_mesh + 1j * E_imag_mesh).flatten()
    
    mu1 = []
    mu2 = []
    for item in res:
        if item.success and item.is_gbz:
            for point in item.subsets:
                mu1.append(np.log(abs(point.beta1)))
                mu2.append(np.log(abs(point.beta2)))
    print(len(mu1), len(mu2))
    plt.figure()
    plt.plot(mu1, mu2, '.')


def plot_SGBZ_mu(which="a1"):
    with open("data/Haldane-gain-loss-%s-SGBZ.pkl" % which, "rb") as fp:
        data = pickle.load(fp)
    E_real, E_imag, results, params = data["E_real"], data["E_imag"], data["results"], data["params"]
    all_mu1 = []

    E_real_mesh, E_imag_mesh = np.meshgrid(E_real, E_imag)
    E_list = (E_real_mesh + 1j * E_imag_mesh).flatten()
    E_SGBZ = []
    for ind, res in enumerate(results):
        if not res.is_empty:
            if res["is_PMGBZ"]:
                all_mu1.append(res["mu1"])
                E_SGBZ.append(E_list[ind])
    E_SGBZ = np.asarray(E_SGBZ)
    sca = plt.scatter(E_SGBZ.real, E_SGBZ.imag, c=all_mu1)
    plt.colorbar(sca)


def plot_index_E(which="", kind="amoeba"):
    """子集数-E 图：每个 GBZResult 的 index (n_0D, n_1D) 在能量复平面上着色。

    每种不同的 index（如 (1,0)、(0,2)、(0,0)）映射为一种颜色；(0,0) 表示
    E 在 GBZ 外（或计算失败）。index 是 ``GBZResult.index`` = (PointSubset
    数, LineSubset 数)。

    Parameters
    ----------
    which : str
        数据文件后缀。amoeba 对应 ``Haldane-gain-loss-amoeba<which>.pkl``
        （如 ``""`` 或 ``"-xy"``）；SGBZ 对应
        ``Haldane-gain-loss-<which>-SGBZ.pkl``（如 ``"a1"``、``"a2"``、
        ``"x"``、``"y"``）。
    kind : str
        ``"amoeba"`` 或 ``"SGBZ"``（大小写不敏感），选择数据命名约定。
    """
    kind = kind.lower()
    if kind == "amoeba":
        fname = "data/Haldane-gain-loss-amoeba%s.pkl" % which
    elif kind == "sgbz":
        fname = "data/Haldane-gain-loss-%s-SGBZ.pkl" % which
    else:
        raise ValueError(f"unknown kind={kind!r}; use 'amoeba' or 'SGBZ'")

    with open(fname, "rb") as fp:
        data = pickle.load(fp)
    E_real, E_imag, results = data["E_real"], data["E_imag"], data["results"]
    E_real_mesh, E_imag_mesh = np.meshgrid(E_real, E_imag)
    E_list = (E_real_mesh + 1j * E_imag_mesh).flatten()

    def get_index(res):
        return tuple(res.index)

    index_keys = [get_index(res) for res in results]
    unique = sorted(set(index_keys))
    key_num = {k: n for n, k in enumerate(unique)}
    cvals = np.array([key_num[k] for k in index_keys])
    print(key_num, cvals)

    plt.figure()
    from matplotlib.colors import ListedColormap
    base = plt.cm.get_cmap("tab20" if len(unique) <= 20 else "viridis")
    cmap = ListedColormap([base(i % base.N) for i in range(len(unique))])
    sca = plt.scatter(E_list.real, E_list.imag, c=cvals, cmap=cmap, s=8)
    cbar = plt.colorbar(sca, ticks=range(len(unique)))
    cbar.ax.set_yticklabels([f"({a},{b})" for a, b in unique])
    plt.xlabel("Re E")
    plt.ylabel("Im E")
    plt.title(f"{kind.upper()} GBZ subsets index $(n_{{0D}}, n_{{1D}})$ vs $E$")


def debug_y_SGBZ():
    model = Haldane_non_Hermitian_phase(*ALL_PARAMS)
    model = model.get_supercell(
        [(0, 0), (1, 0)],
        np.array([
            [1, 1],
            [-1, 1]
        ], dtype=int)
    )
    coeffs, degs = model.get_characteristic_polynomial_data()
    E_ref = -1.098 + 0.2445j
    print(bfs.collect_GBZ_subsets(coeffs, degs, E_ref, debug_mode=True))



if __name__ == "__main__":
    sweep_amoeba()
    sweep_amoeba_multiband()
    sweep_SGBZ_a1()
    sweep_SGBZ_a2()
    sweep_SGBZ_x()
    sweep_SGBZ_y()
    # recompute_failed_SGBZ("y")
    # recompute_failed_SGBZ("y", out_fname="data/Haldane-gain-loss-y-SGBZ-recomputed.pkl")
    # recompute_failed_SGBZ("x", out_fname="data/Haldane-gain-loss-x-SGBZ-recomputed.pkl")
    # recompute_failed_SGBZ("a1", out_fname="data/Haldane-gain-loss-a1-SGBZ-recomputed.pkl")
    # recompute_failed_SGBZ("a2", out_fname="data/Haldane-gain-loss-a2-SGBZ-recomputed.pkl")
    # plot_amoebic_spectrum()
    # plot_amoeba_mu()
    # plot_amoebic_spectrum("-xy")
    # plot_SGBZ("a1")
    # plot_SGBZ("a2")
    # plot_SGBZ("x")
    # plot_SGBZ("y")
    # plot_index_E("a1", kind="SGBZ")
    # plot_index_E("a2", kind="SGBZ")
    # plot_index_E("x", kind="SGBZ")
    # plot_index_E("y", kind="SGBZ")
    # plot_index_E("")
    # plot_index_E("-xy")
    # plt.show()
    # debug_y_SGBZ()
