'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-06-01
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''


import BerryPy.TightBinding as tb
import numpy as np
import multiprocessing as mp
import pickle
import matplotlib.pyplot as plt
import poly_tools as pt
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import brute_force_SGBZ as bfs
print("GBZ path: ", bfs.__file__)
import brute_force_amoeba as bfa

DEFAULT_PARAMS = {
    "t": 1,
    "gamma": 0.1,
    "s": 0.3,
}


def get_model(t, gamma, s):
    dim = 2
    site_num = 1
    inter_cell = [
        (0, 0, t - gamma, (-1, 0)),
        (0, 0, t + gamma, (1, 0)),
        (0, 0, t - gamma, (0, -1)),
        (0, 0, t + gamma, (0, 1)),
        (0, 0, s, (1, 1)),
        (0, 0, s, (1, -1)),
        (0, 0, s, (-1, 1)),
        (0, 0, s, (-1, -1)),
    ]

    return tb.TightBindingModel(dim, site_num, np.eye(2), [], inter_cell)


def sweep_general(E_re, E_im, which="amoeba", params=DEFAULT_PARAMS, N_process=10):
    import os
    # Prepare E
    E_re_mesh, E_im_mesh = np.meshgrid(E_re, E_im)
    E_mesh = E_re_mesh + 1j * E_im_mesh
    E_list = E_mesh.flatten()

    # Get characteristic polynomial
    model = get_model(**params)
    coeffs, degs = model.get_characteristic_polynomial_data()

    # Run

    if which == "amoeba":
        fname_prefix = "data/ZWang_amoeba"
        data_pack = [(coeffs, degs, E_list[j], j / len(E_list), True) for j in range(len(E_list))]
        with mp.Pool(N_process) as pool:
            res = pool.starmap(bfa.collect_GBZ_subsets, data_pack)
    elif which == "x-SGBZ":
        fname_prefix = "data/ZWang_x-SGBZ"
        data_pack = [(coeffs, degs, E_list[j], j / len(E_list), True) for j in range(len(E_list))]
        with mp.Pool(N_process) as pool:
            res = pool.starmap(bfs.collect_GBZ_subsets, data_pack)
    elif which == "y-SGBZ":
        fname_prefix = "data/ZWang_y-SGBZ"
        degs = degs[:, [0, 2, 1]]
        data_pack = [(coeffs, degs, E_list[j], j / len(E_list), True) for j in range(len(E_list))]
        with mp.Pool(N_process) as pool:
            res = pool.starmap(bfs.collect_GBZ_subsets, data_pack)
    else:
        raise ValueError(f"Unknown type: {which}")

    fname = f"{fname_prefix}.pkl"
    with open(fname, "wb") as fp:
        pickle.dump({
            "E_re": E_re,
            "E_im": E_im,
            "params": params,
            "res": res,
        }, fp)


def plot_amoeba():
    # Load data
    with open(f"data/ZWang_amoeba.pkl", "rb") as fp:
        data = pickle.load(fp)
    E_re = data["E_re"]
    E_im = data["E_im"]
    E_re_mesh, E_im_mesh = np.meshgrid(E_re, E_im)
    E_mesh = E_re_mesh + 1j * E_im_mesh
    E_list = E_mesh.flatten()
    params = data["params"]
    res = data["res"]

    model = get_model(**params)
    coeffs, degs = model.get_characteristic_polynomial_data()
    ind_failed = [i for i in range(len(res)) if not res[i].success]
    ind_amoeba = [i for i in range(len(res)) if res[i].is_gbz]
    ind_not_amoeba = [i for i in range(len(res)) if not res[i].is_gbz and res[i].success]

    # Plot
    plt.figure()
    plt.plot(E_list[ind_amoeba].real, E_list[ind_amoeba].imag, '.', label="Amoeba")
    plt.plot(E_list[ind_not_amoeba].real, E_list[ind_not_amoeba].imag, '.', label="Non-Amoeba")
    plt.plot(E_list[ind_failed].real, E_list[ind_failed].imag, '.', label="Failed")
    plt.legend()
    plt.show()


def plot_SGBZ(which="x"):
    # Load data
    with open(f"data/ZWang_{which}-SGBZ.pkl", "rb") as fp:
        data = pickle.load(fp)
    E_re = data["E_re"]
    E_im = data["E_im"]
    E_re_mesh, E_im_mesh = np.meshgrid(E_re, E_im)
    E_mesh = E_re_mesh + 1j * E_im_mesh
    E_list = E_mesh.flatten()
    params = data["params"]
    res = data["res"]

    ind_failed = [i for i in range(len(res)) if not res[i].is_gbz and not res[i].success]
    ind_amoeba = [i for i in range(len(res)) if res[i].is_gbz]
    ind_not_amoeba = [i for i in range(len(res)) if not res[i].is_gbz and res[i].success]

    # Plot
    plt.figure()
    plt.plot(E_list[ind_amoeba].real, E_list[ind_amoeba].imag, '.', label="SGBZ")
    plt.plot(E_list[ind_not_amoeba].real, E_list[ind_not_amoeba].imag, '.', label="Non-SGBZ")
    plt.plot(E_list[ind_failed].real, E_list[ind_failed].imag, '.', label="Failed")
    plt.legend()
    plt.show()


def check_roots():
    model = get_model(**DEFAULT_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()
    char_poly = pt.CLaurent(3)
    char_poly.set_Laurent_by_terms(
        pt.CScalarVec(coeffs),
        pt.CLaurentIndexVec(degs.flatten())
    )

    E_ref = -2.4 - 0.1408j
    root_track = bfa._compute_root_tracks(char_poly, E_ref, 0.0)
    # print(root_track)
    theta1_ext = root_track["theta1_ext"]
    tracked_ext = root_track["tracked_ext"]
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    for j in range(2):
        ax.plot(theta1_ext, np.angle(tracked_ext[:,j]), np.log(np.abs(tracked_ext[:,j])))
    plt.show()

    print(bfa.collect_GBZ_subsets(coeffs, degs, E_ref, 0, True))


if __name__ == "__main__":
    N_grid_x = 31
    N_grid_y = 31
    E_re = np.linspace(-4, 6, N_grid_x)
    E_im = np.linspace(-0.16, 0.16, N_grid_y)
    # E_im = np.concatenate([
    #     np.linspace(-0.16, -0.08, N_grid_y),
    #     np.linspace(0.08, 0.16, N_grid_y),
    # ])
    # sweep_general(E_re, E_im, which="amoeba")
    # sweep_general(E_re, E_im, which="x-SGBZ")
    # sweep_general(E_re, E_im, which="y-SGBZ")
    plot_amoeba()
    # plot_SGBZ("x")
    # plot_SGBZ("y")
    plt.show()

    # check_roots()
