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
DEFAULT_RCOND = 1e-6
LINK_EPS = 1e-12


@dataclass
class VertexDiagnostics:
    n_vertices: int
    nullity_one: int
    nullity_zero: int
    nullity_multi: int
    residual_max: float
    residual_p50: float
    residual_p99: float
    sigma_min_max: float
    sigma_min_p50: float
    gap_ratio_min: float
    gap_ratio_p01: float


@dataclass
class MeshChernResult:
    mesh: str
    n_vertices: int
    n_triangles: int
    orientation: str
    n_orientation_flips: int
    chern: float
    total_flux: float
    flux_p50_abs: float
    flux_p99_abs: float
    flux_max_abs: float
    bad_links: int
    vertex_diagnostics: VertexDiagnostics


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


def betas_from_mesh(verts: np.ndarray, mu1: np.ndarray,
                    mu2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    beta1 = np.exp(mu1 + 1j * verts[:, 0])
    beta2 = np.exp(mu2 + 1j * verts[:, 1])
    return beta1, beta2


def svd_null_vector(A: np.ndarray, rcond: float) -> tuple[np.ndarray, int, float, float]:
    """Smallest right-singular vector plus nullity diagnostics."""
    u, s, vh = la.svd(A, full_matrices=True, check_finite=False)
    del u
    scale = float(s[0]) if len(s) else 0.0
    tol = rcond * scale
    nullity = int(np.count_nonzero(s <= tol))
    vec = vh.conj().T[:, -1]
    vec = vec / la.norm(vec)
    sigma_min = float(s[-1]) if len(s) else 0.0
    if len(s) >= 2 and sigma_min > 0:
        gap_ratio = float(s[-2] / sigma_min)
    elif len(s) >= 2:
        gap_ratio = float("inf")
    else:
        gap_ratio = float("nan")
    return vec, nullity, sigma_min, gap_ratio


def right_eigenvectors_on_mesh(Hfun, E: np.ndarray, beta1: np.ndarray,
                               beta2: np.ndarray, rcond: float,
                               progress: int = 0):
    """Compute one right null vector of E I - H at every mesh vertex."""
    vecs = []
    nullities = np.empty(len(E), dtype=int)
    residuals = np.empty(len(E), dtype=float)
    sigma_min = np.empty(len(E), dtype=float)
    gap_ratio = np.empty(len(E), dtype=float)

    for j, (Ej, b1, b2) in enumerate(zip(E, beta1, beta2)):
        H = Hfun(b1, b2)
        A = np.eye(H.shape[0], dtype=complex) * Ej - H
        vec, nullity, smin, gratio = svd_null_vector(A, rcond)
        vecs.append(vec)
        nullities[j] = nullity
        residuals[j] = la.norm(A @ vec)
        sigma_min[j] = smin
        gap_ratio[j] = gratio
        if progress and (j + 1) % progress == 0:
            print(f"    eigenvectors {j + 1}/{len(E)}", flush=True)

    return np.asarray(vecs), make_vertex_diagnostics(
        nullities, residuals, sigma_min, gap_ratio
    )


def left_right_eigenvectors_on_mesh(Hfun, E: np.ndarray, beta1: np.ndarray,
                                    beta2: np.ndarray, rcond: float,
                                    progress: int = 0):
    """Compute biorthogonal left/right null vectors of E I - H."""
    right = []
    left = []
    nullities = np.empty(len(E), dtype=int)
    residuals = np.empty(len(E), dtype=float)
    sigma_min = np.empty(len(E), dtype=float)
    gap_ratio = np.empty(len(E), dtype=float)

    for j, (Ej, b1, b2) in enumerate(zip(E, beta1, beta2)):
        H = Hfun(b1, b2)
        A = np.eye(H.shape[0], dtype=complex) * Ej - H
        rv, rn, smin, gratio = svd_null_vector(A, rcond)
        lv, ln, _, _ = svd_null_vector(A.conj().T, rcond)
        overlap = np.vdot(lv, rv)
        if abs(overlap) > LINK_EPS:
            lv = lv / overlap.conjugate()
        right.append(rv)
        left.append(lv)
        nullities[j] = min(rn, ln)
        residuals[j] = max(la.norm(A @ rv), la.norm(A.conj().T @ lv))
        sigma_min[j] = smin
        gap_ratio[j] = gratio
        if progress and (j + 1) % progress == 0:
            print(f"    eigenvectors {j + 1}/{len(E)}", flush=True)

    return np.asarray(right), np.asarray(left), make_vertex_diagnostics(
        nullities, residuals, sigma_min, gap_ratio
    )


def make_vertex_diagnostics(nullities: np.ndarray, residuals: np.ndarray,
                            sigma_min: np.ndarray,
                            gap_ratio: np.ndarray) -> VertexDiagnostics:
    finite_gap = gap_ratio[np.isfinite(gap_ratio)]
    if len(finite_gap) == 0:
        finite_gap = np.array([np.inf])
    return VertexDiagnostics(
        n_vertices=int(len(nullities)),
        nullity_one=int(np.count_nonzero(nullities == 1)),
        nullity_zero=int(np.count_nonzero(nullities == 0)),
        nullity_multi=int(np.count_nonzero(nullities > 1)),
        residual_max=float(np.max(residuals)),
        residual_p50=float(np.percentile(residuals, 50)),
        residual_p99=float(np.percentile(residuals, 99)),
        sigma_min_max=float(np.max(sigma_min)),
        sigma_min_p50=float(np.percentile(sigma_min, 50)),
        gap_ratio_min=float(np.min(finite_gap)),
        gap_ratio_p01=float(np.percentile(finite_gap, 1)),
    )


def unwrap_triangle(verts: np.ndarray, tri: np.ndarray) -> np.ndarray:
    """Triangle coordinates in one local sheet of the universal cover."""
    p = verts[tri].copy()
    p[:, 1:] -= np.round((p[:, 1:] - p[:, :1]) / TWO_PI) * TWO_PI
    return p


def orient_triangles(verts: np.ndarray, triangles: np.ndarray,
                     positive: bool = True) -> tuple[np.ndarray, int, str]:
    """Orient triangles by the standard theta1,theta2 orientation."""
    tri = np.asarray(triangles, dtype=int).copy()
    p = unwrap_triangle(verts, tri)
    cross = ((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
             - (p[:, 1, 1] - p[:, 0, 1]) * (p[:, 2, 0] - p[:, 0, 0]))
    if positive:
        flip = cross < 0
        orientation = "positive"
    else:
        flip = cross > 0
        orientation = "negative"
    tri[flip, 1], tri[flip, 2] = tri[flip, 2], tri[flip, 1].copy()
    return tri, int(np.count_nonzero(flip)), orientation


def triangle_flux_right(vecs: np.ndarray, tri: np.ndarray) -> tuple[np.ndarray, int]:
    """Wilson-loop Berry flux for right eigenvectors on every triangle."""
    flux = np.empty(len(tri), dtype=float)
    bad = 0
    for n, (i, j, k) in enumerate(tri):
        z01 = np.vdot(vecs[i], vecs[j])
        z12 = np.vdot(vecs[j], vecs[k])
        z20 = np.vdot(vecs[k], vecs[i])
        if min(abs(z01), abs(z12), abs(z20)) <= LINK_EPS:
            bad += 1
        flux[n] = -np.angle(z01 * z12 * z20)
    return flux, bad


def triangle_flux_biorthogonal(right: np.ndarray, left: np.ndarray,
                               tri: np.ndarray) -> tuple[np.ndarray, int]:
    """Biorthogonal Wilson-loop Berry flux on every triangle."""
    flux = np.empty(len(tri), dtype=float)
    bad = 0
    for n, (i, j, k) in enumerate(tri):
        z01 = np.vdot(left[i], right[j])
        z12 = np.vdot(left[j], right[k])
        z20 = np.vdot(left[k], right[i])
        if min(abs(z01), abs(z12), abs(z20)) <= LINK_EPS:
            bad += 1
        flux[n] = -np.angle(z01 * z12 * z20)
    return flux, bad


def calculate_chern(mesh_file: Path, Hfun, *, rcond: float = DEFAULT_RCOND,
                    mode: str = "right", orientation: str = "positive",
                    progress: int = 0, save_flux: bool = False) -> MeshChernResult:
    """Calculate Chern number for one saved 2D-GBZ mesh file."""
    z = np.load(mesh_file)
    verts = np.asarray(z["verts"], dtype=float)
    triangles = np.asarray(z["triangles"], dtype=int)
    E = np.asarray(z["E"], dtype=complex)
    mu1 = np.asarray(z["mu1"], dtype=float)
    mu2 = np.asarray(z["mu2"], dtype=float)
    beta1, beta2 = betas_from_mesh(verts, mu1, mu2)

    if orientation == "keep":
        oriented = triangles
        flips = 0
        orient_label = "kept"
    else:
        oriented, flips, orient_label = orient_triangles(
            verts, triangles, positive=(orientation == "positive")
        )

    if mode == "right":
        vecs, diag = right_eigenvectors_on_mesh(
            Hfun, E, beta1, beta2, rcond, progress=progress
        )
        flux, bad = triangle_flux_right(vecs, oriented)
    elif mode == "biorthogonal":
        right, left, diag = left_right_eigenvectors_on_mesh(
            Hfun, E, beta1, beta2, rcond, progress=progress
        )
        flux, bad = triangle_flux_biorthogonal(right, left, oriented)
    else:
        raise ValueError(f"unknown mode {mode!r}")

    total_flux = float(np.sum(flux))
    result = MeshChernResult(
        mesh=str(mesh_file),
        n_vertices=int(len(verts)),
        n_triangles=int(len(oriented)),
        orientation=orient_label,
        n_orientation_flips=flips,
        chern=total_flux / TWO_PI,
        total_flux=total_flux,
        flux_p50_abs=float(np.percentile(np.abs(flux), 50)),
        flux_p99_abs=float(np.percentile(np.abs(flux), 99)),
        flux_max_abs=float(np.max(np.abs(flux))) if len(flux) else 0.0,
        bad_links=bad,
        vertex_diagnostics=diag,
    )

    if save_flux:
        out = mesh_file.with_name(mesh_file.stem + "_chern_flux.npz")
        np.savez(out, flux=flux, oriented_triangles=oriented,
                 beta1=beta1, beta2=beta2)
        print(f"    flux saved -> {out}")
    return result


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


def regular_torus_mesh(n: int):
    """Periodic n x n square grid split into two oriented triangles/cell."""
    t = np.linspace(0.0, TWO_PI, n, endpoint=False)
    x, y = np.meshgrid(t, t, indexing="ij")
    verts = np.column_stack([x.ravel(), y.ravel()])
    def idx(i, j):
        return (i % n) * n + (j % n)
    tri = []
    for i in range(n):
        for j in range(n):
            tri.append([idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)])
            tri.append([idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)])
    return verts, np.array(tri, dtype=int)


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
