"""Data-integrity checks for the application; no optional BerryPy dependency."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from pygbz2d.core import GBZResult


@pytest.fixture
def example():
    path = Path(__file__).resolve().parents[1] / "application/geometry-dependent-skin-effect.py"
    spec = importlib.util.spec_from_file_location("gdse_example", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_density_uses_normalized_right_vectors_before_averaging(example):
    eigenvalues = np.array([-2.2, -1.8, 1.0], complex)
    vectors = np.array([[1, 1, 0], [0, 1j, 0], [0, 0, 1]], complex)
    coords = np.array([[0, 0], [1, 0], [2, 0]])
    summary = example.summarize_eigensystem(eigenvalues, vectors, coords)
    rescaled = example.summarize_eigensystem(eigenvalues, vectors * [3j, .02, -8], coords)
    np.testing.assert_allclose(summary["density_window"], [.75, .25, 0])
    np.testing.assert_allclose(summary["density_all"], [.5, 1/6, 1/3])
    for key in ("density_window", "density_all", "ipr", "selected_density"):
        np.testing.assert_allclose(summary[key], rescaled[key])
    with pytest.raises(ValueError, match="No states"):
        example.summarize_eigensystem(eigenvalues, vectors, coords, window=(4, 5))


def test_sweep_rejects_misordered_energy_cache(example, tmp_path):
    path = tmp_path / "misordered.pkl"
    example.atomic_pickle(path, dict(E_real=[1., 2.], E_imag=[0.],
                                    results=[GBZResult(2), GBZResult(1)]))
    with pytest.raises(ValueError, match="result order"):
        example.load_sweep(path)


def test_resume_preserves_failures_until_retry_is_requested(example, tmp_path, monkeypatch):
    coefficients, degrees = np.array([1.], complex), np.array([[1, 0, 0]])
    monkeypatch.setattr(example, "direction_model", lambda _: None)
    monkeypatch.setattr(example, "polynomial_data", lambda _: (coefficients, degrees))
    calls = []

    def solve(task):
        index, energy = task
        calls.append(index)
        return index, GBZResult(energy)

    monkeypatch.setattr(example, "_solve_energy", solve)
    path = tmp_path / "checkpoint.pkl"
    example.atomic_pickle(path, dict(E_real=np.array([1., 2., 3.]), E_imag=np.array([0.]),
                                    results=[GBZResult(1), GBZResult(2, success=False, error="test failure"), None],
                                    coeffs=coefficients, degs=degrees, params=example.PARAMS, direction="a1"))
    data = example.sweep("a1", np.array([1., 2., 3.]), np.array([0.]), path)
    assert calls == [2]
    assert example.sweep_report(data)["failed"] == 1
    calls.clear()
    example.sweep("a1", np.array([1., 2., 3.]), np.array([0.]), path, retry_failed=True)
    assert calls == [1]
    assert example.sweep_report(example.load_sweep(path))["failed"] == 0
    with pytest.raises(ValueError, match="incompatible checkpoint"):
        example.sweep("x", np.array([1., 2., 3.]), np.array([0.]), path)
