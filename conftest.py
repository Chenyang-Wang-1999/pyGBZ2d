import sys
from pathlib import Path

try:
    import pygbz2d  # noqa: F401  — prefer the installed package when present
except ModuleNotFoundError:
    # Fallback for running tests from a source checkout without
    # ``pip install -e .`` (src layout): prepend src/ so ``import pygbz2d``
    # resolves to the same files the installed package would provide.
    sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

import pytest

# NOTE: tests import ``pygbz2d`` (installed package preferred, src/ fallback).
# Shared HN-2D test models below were previously duplicated (with drift)
# across test files.
from cmath import exp as _exp

import numpy as np


def build_HN2D_polynomial(J1, J2, gamma_1, gamma_2, delta_1, delta_2,
                          basis="10"):
    """2D Hatano-Nelson chain-pair characteristic polynomial.

    Returns (coeffs, degs) for f = E − J11 β₁^{p1} − J12 β₁^{p2} −
    J21 β₂⁻¹ − J22 β₂ with anisotropic non-reciprocal hoppings.
    ``basis`` selects the unit cell of the analytic solution the caller
    compares against: "10" puts both β₁ terms at (−1, 0) / (+1, 0);
    "11" couples them to opposite β₂ directions, (−1, +1) / (+1, −1).
    """
    J11 = _exp(gamma_1 + 1j * delta_1) * J1
    J12 = _exp(-gamma_1 + 1j * delta_1) * np.conj(J1)
    J21 = _exp(gamma_2 + 1j * delta_2) * J2
    J22 = _exp(-gamma_2 + 1j * delta_2) * np.conj(J2)

    coeffs = np.array([1, -J11, -J12, -J21, -J22], dtype=complex)

    if basis == "10":
        degs = np.array([
            [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1],
        ], dtype=int)
    elif basis == "11":
        degs = np.array([
            [1, 0, 0], [0, -1, 1], [0, 1, -1], [0, 0, -1], [0, 0, 1],
        ], dtype=int)
    else:
        raise ValueError(f"Unknown basis: {basis}")

    return coeffs, degs


def pytest_addoption(parser):
    parser.addoption(
        "--run-slow", action="store_true", default=False,
        help="run tests marked with @pytest.mark.slow",
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: slow test (skipped by default, use --run-slow)")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-slow"):
        return
    skip_slow = pytest.mark.skip(reason="need --run-slow to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
