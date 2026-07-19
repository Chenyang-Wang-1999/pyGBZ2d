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
from typing import Literal
import multiprocessing as mp
import os
from scipy import interpolate
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import brute_force_SGBZ as bfs
import brute_force_amoeba as bfa
print("Using ", bfs.__file__)

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

    N_points = 51
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

    N_points = 51
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

    N_points = 51
    E_real = np.linspace(-3.1, 4.6, N_points)
    E_imag = np.linspace(-0.51, 0.51, N_points)
    E_real_mesh, E_imag_mesh = np.meshgrid(E_real, E_imag)
    E_mesh = E_real_mesh + 1j * E_imag_mesh
    E_list = E_mesh.flatten()

    pool = mp.Pool(mp.cpu_count())
    results = pool.starmap(bfs.collect_GBZ_subsets, [(coeffs, degs, E, j / len(E_list), True) for j, E in enumerate(E_list)])
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

    N_points = 51
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

    N_points = 51
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

    N_points = 51
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


if __name__ == "__main__":
    # sweep_amoeba()
    # sweep_amoeba_multiband()
    sweep_SGBZ_a1()
    sweep_SGBZ_a2()
    sweep_SGBZ_x()
    sweep_SGBZ_y()
    # plot_amoebic_spectrum()
    # plot_amoeba_mu()
    # plot_amoebic_spectrum("-xy")
    # plot_SGBZ("a1")
    # plot_SGBZ("a2")
    # plot_SGBZ("x")
    # plot_SGBZ("y")
    # plt.show()
