'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-09-07
Copyright © Department of Physics, Tsinghua University. All rights reserved

EXPERIMENTAL playground Chern-number calculation on a triangulated 2D GBZ.

This script consumes the mesh files produced by playground/demo_pipeline.py:

    mesh_cluster{cid}.npz with arrays
      verts      (n, 2)   theta1, theta2 in [0, 2pi)
      triangles  (m, 3)   torus triangle vertex ids
      E          (n,)     complex energy on this band
      mu1, mu2   (n,)     log |beta1|, log |beta2|

For each vertex it reconstructs beta_j = exp(mu_j + i theta_j), evaluates
H(beta1, beta2) through a BerryPy-style model, extracts a right null vector
of E I - H by SVD, then sums the gauge-invariant Wilson phase around every
oriented triangle:

    flux(T) = -arg <u0|u1><u1|u2><u2|u0>.

Biorthogonal mode instead uses Uij = <Li|Rj>, <Li|Ri> = 1, and
Tij = Uij / sqrt(Uij Uji), with one square root per undirected edge and
Tji = 1/Tij. Its complex flux is i Log(Tij Tjk Tki). Raw LR overlap
products contain a second-order complex-metric contribution and cannot be
used as local curvature on shrinking triangles. See doc/experimental.md for the
derivation and references (Fukui-Hatsugai-Suzuki; Shen-Zhen-Fu).

The default model loader matches playground/Haldane-model-gainloss.py and
ALL_PARAMS, because that is the source of the current demo meshes.  For other
models, pass --model-file and --factory NAME; the factory must be a zero-arg
function returning either a BerryPy-like model with get_bulk_Hamiltonian_dense
or a callable H(beta1, beta2).

Usage examples:
    python playground/chern_on_gbz.py \
        playground/band_cluster_out/Haldane-gain-loss-amoeba-enriched/mesh_cluster0.npz

    python playground/chern_on_gbz.py playground/band_cluster_out/.../mesh_cluster*.npz \
        --save-flux

    python playground/chern_on_gbz.py --self-test
'''

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Callable

import numpy as np
from scipy import linalg as la

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pygbz2d.core import TWO_PI


DEFAULT_MODEL_FILE = Path(__file__).resolve().parent / "Haldane-model-gainloss.py"


from pygbz2d.experimental.chern import (
    DEFAULT_RCOND, LINK_EPS, LINK_BRANCH_CUT_TOL,
    ChernResult, VertexDiagnostics, betas_from_mesh, svd_null_vector,
    right_eigenvectors_on_mesh, left_right_eigenvectors_on_mesh,
    make_vertex_diagnostics, triangle_flux_right, biorthogonal_edge_links,
    triangle_flux_biorthogonal_complex, triangle_flux_biorthogonal, integrate_chern,
)
from pygbz2d.experimental.torus_mesh import (
    unwrap_triangle, orient_triangles, regular_torus_mesh,
)


@dataclass
class MeshChernResult(ChernResult):
    """Demo file label attached to the library's array-based result."""
    mesh: str = ""


def calculate_chern(mesh_file: Path, Hfun, *, rcond=None, mode="right",
                    orientation="positive", progress=0, save_flux=False):
    """Demo file I/O wrapper; all numerical integration lives in the package."""
    mesh_file = Path(mesh_file)
    with np.load(mesh_file) as data:
        result, flux = integrate_chern(
            data["verts"], data["triangles"], data["E"], data["mu1"], data["mu2"],
            Hfun, rcond=rcond, mode=mode, orientation=orientation, progress=progress)
    if save_flux:
        out = mesh_file.with_name(mesh_file.stem + "_chern_flux.npz")
        np.savez(out, **flux)
        print(f"    flux saved -> {out}")
    return MeshChernResult(**vars(result), mesh=str(mesh_file))


def import_module_from_file(path: Path):
    """Import a Python file without requiring its filename to be a module id."""
    path = Path(path).resolve()
    spec = importlib.util.spec_from_file_location(path.stem.replace("-", "_"), path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_model(model_file: Path | None, factory: str | None,
               supercell: str):
    """Load a model object for Hamiltonian evaluation."""
    path = Path(model_file) if model_file is not None else DEFAULT_MODEL_FILE
    module = import_module_from_file(path)

    if factory is not None:
        model = getattr(module, factory)()
    elif hasattr(module, "Haldane_non_Hermitian_phase") and hasattr(module, "ALL_PARAMS"):
        model = module.Haldane_non_Hermitian_phase(*module.ALL_PARAMS)
    elif hasattr(module, "make_model"):
        model = module.make_model()
    else:
        raise ValueError(
            f"{path} needs --factory NAME, or a make_model() function, or "
            "Haldane_non_Hermitian_phase(*ALL_PARAMS)"
        )

    if supercell != "none":
        if not hasattr(model, "get_supercell"):
            raise ValueError("--supercell requires a model with get_supercell")
        if supercell == "a2":
            model = model.get_supercell([(0, 0)], np.array([[0, 1], [1, 0]], dtype=int))
        elif supercell in ("x", "xy"):
            model = model.get_supercell(
                [(0, 0), (1, 0)], np.array([[1, 1], [1, -1]], dtype=int)
            )
        elif supercell == "y":
            model = model.get_supercell(
                [(0, 0), (1, 0)], np.array([[1, 1], [-1, 1]], dtype=int)
            )
        else:
            raise ValueError(f"unknown supercell preset {supercell!r}")
    return model


def beta_to_k_cycles(beta1: complex, beta2: complex) -> tuple[complex, complex]:
    """Return k in BerryPy's cycle convention, beta = exp(2 pi i k)."""
    return tuple(np.log(np.array([beta1, beta2], dtype=complex)) / (2j * np.pi))


def make_hamiltonian_function(model) -> Callable[[complex, complex], np.ndarray]:
    """Adapt supported model styles to H(beta1, beta2) -> dense ndarray."""
    if hasattr(model, "get_bulk_Hamiltonian_dense"):
        def H(beta1, beta2):
            return np.asarray(
                model.get_bulk_Hamiltonian_dense(beta_to_k_cycles(beta1, beta2)),
                dtype=complex,
            )
        return H

    if hasattr(model, "get_bulk_Hamiltonian"):
        def H(beta1, beta2):
            return np.asarray(
                model.get_bulk_Hamiltonian(beta_to_k_cycles(beta1, beta2)).todense(),
                dtype=complex,
            )
        return H

    if callable(model):
        def H(beta1, beta2):
            return np.asarray(model(beta1, beta2), dtype=complex)
        return H

    raise TypeError(
        "model must be callable or expose get_bulk_Hamiltonian_dense/get_bulk_Hamiltonian"
    )


class QWZModel:
    """Small self-test model with Chern bands for m in (-2, 0)."""
    def __init__(self, m=-1.0):
        self.m = m

    def __call__(self, beta1, beta2):
        kx = np.angle(beta1)
        ky = np.angle(beta2)
        sx = np.array([[0, 1], [1, 0]], dtype=complex)
        sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
        sz = np.array([[1, 0], [0, -1]], dtype=complex)
        return (np.sin(kx) * sx + np.sin(ky) * sy
                + (self.m + np.cos(kx) + np.cos(ky)) * sz)


def self_test(n: int = 31):
    """Run the Wilson-loop implementation on a known Hermitian Chern model."""
    model = QWZModel(m=-1.0)
    Hfun = make_hamiltonian_function(model)
    verts, tri = regular_torus_mesh(n)
    beta1 = np.exp(1j * verts[:, 0])
    beta2 = np.exp(1j * verts[:, 1])
    E_all = np.array([la.eigvalsh(Hfun(b1, b2)) for b1, b2 in zip(beta1, beta2)])
    tmp = Path(__file__).resolve().parent / "_chern_self_test_mesh.npz"
    try:
        for band, expected in [(0, 1), (1, -1)]:
            np.savez(tmp, verts=verts, triangles=tri, E=E_all[:, band],
                     mu1=np.zeros(len(verts)), mu2=np.zeros(len(verts)))
            res = calculate_chern(tmp, Hfun, rcond=1e-5, progress=0)
            print(
                f"self-test band {band}: C={res.chern:+.6f} "
                f"(nearest {round(res.chern):+d}, expected about {expected:+d})"
            )
    finally:
        try:
            tmp.unlink()
        except OSError:
            pass


def print_result(res: MeshChernResult):
    diag = res.vertex_diagnostics
    print(f"== {res.mesh}")
    print(
        f"C = {res.chern:+.12f}  "
        f"(flux={res.total_flux:+.12f}, nearest={round(res.chern):+d})"
    )
    print(f"links: {res.link_rule}; Im(total flux)={res.total_flux_imag:+.3e}")
    print(
        f"mesh: V={res.n_vertices} F={res.n_triangles}, "
        f"orientation={res.orientation}, flips={res.n_orientation_flips}"
    )
    print(
        f"triangle |flux| p50/p99/max = "
        f"{res.flux_p50_abs:.3e} / {res.flux_p99_abs:.3e} / {res.flux_max_abs:.3e}; "
        f"bad links={res.bad_links}"
    )
    print(
        f"nullity: one={diag.nullity_one}, zero={diag.nullity_zero}, "
        f"multi={diag.nullity_multi}; residual p50/p99/max = "
        f"{diag.residual_p50:.3e} / {diag.residual_p99:.3e} / {diag.residual_max:.3e}"
    )
    print(
        f"sigma_min p50/max = {diag.sigma_min_p50:.3e} / "
        f"{diag.sigma_min_max:.3e}; gap ratio p01/min = "
        f"{diag.gap_ratio_p01:.3e} / {diag.gap_ratio_min:.3e}"
    )


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Chern number on triangulated 2D GBZ meshes (playground)"
    )
    p.add_argument("meshes", nargs="*", help="mesh_cluster*.npz files")
    p.add_argument("--model-file", type=Path, default=None,
                   help="Python file defining a model factory")
    p.add_argument("--factory", default=None,
                   help="zero-arg factory function in --model-file")
    p.add_argument("--supercell", choices=("none", "a2", "x", "xy", "y"),
                   default="none", help="preset supercell applied after loading")
    p.add_argument("--mode", choices=("right", "biorthogonal"), default="right",
                   help="Wilson links from right states or biorthogonal states")
    p.add_argument("--orientation", choices=("positive", "negative", "keep"),
                   default="positive",
                   help="triangle orientation used for Wilson loops")
    p.add_argument("--rcond", type=float, default=DEFAULT_RCOND,
                   help="relative singular-value threshold for nullity diagnostics")
    p.add_argument("--progress", type=int, default=5000,
                   help="print eigenvector progress every N vertices; 0 disables")
    p.add_argument("--save-flux", action="store_true",
                   help="save per-triangle flux next to each mesh")
    p.add_argument("--json", type=Path, default=None,
                   help="write result summary JSON")
    p.add_argument("--self-test", action="store_true",
                   help="run a small QWZ-model smoke test")
    args = p.parse_args(argv)

    if args.self_test:
        self_test()
        if not args.meshes:
            return 0

    if not args.meshes:
        p.error("provide at least one mesh .npz, or use --self-test")

    model = load_model(args.model_file, args.factory, args.supercell)
    Hfun = make_hamiltonian_function(model)

    results = []
    for mesh in args.meshes:
        res = calculate_chern(
            Path(mesh), Hfun, rcond=args.rcond, mode=args.mode,
            orientation=args.orientation, progress=args.progress,
            save_flux=args.save_flux,
        )
        print_result(res)
        results.append(res)

    if args.json is not None:
        payload = [asdict(r) for r in results]
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"summary json -> {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
