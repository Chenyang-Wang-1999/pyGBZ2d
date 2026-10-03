'''
author:        Wang Chenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-10-03
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''
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



def main():
    model = get_transformed_model("y")
    coeffs, degs = model.get_characteristic_polynomial_data()
    print(sgbz.collect_GBZ_subsets(coeffs, degs, -1.4675, True))


if __name__ == "__main__":
    main()
