# Copyright © Department of Physics, Tsinghua University. All rights reserved.
# Numerical routines extracted from the playground prototypes by wangchenyang.

"""Experimental single-band Chern integration on triangulated GBZ surfaces.

Array-based numerical API with RR and reciprocal LR links. No file loading,
model discovery, plots, or CLI dependencies. See doc/experimental.md for the
second-order overlap expansion and references (Shen-Zhen-Fu; FHS).
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy import linalg as la
from ..core import TWO_PI
from .torus_mesh import orient_triangles

DEFAULT_RCOND = 1e-6
LINK_EPS = 1e-12
LINK_BRANCH_CUT_TOL = 1e-12

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
class ChernResult:
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
    link_rule: str = "right-overlap"
    total_flux_imag: float = 0.0


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
        if not np.isfinite(overlap) or abs(overlap) <= LINK_EPS:
            raise ValueError(
                f"vertex {j}: left/right self-overlap is zero or ill-conditioned; "
                "an isolated biorthogonal band is required"
            )
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
    # Exactly zero smallest singular values give infinite ratios; percentile
    # interpolation of an all-infinite array otherwise produces NaN.
    gap_min = float(np.min(finite_gap)) if len(finite_gap) else float("inf")
    gap_p01 = float(np.percentile(finite_gap, 1)) if len(finite_gap) else float("inf")
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
        gap_ratio_min=gap_min,
        gap_ratio_p01=gap_p01,
    )


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


def biorthogonal_edge_links(right: np.ndarray, left: np.ndarray,
                            edges: np.ndarray) -> np.ndarray:
    """Reciprocal transports along the supplied edges, with <Li|Ri> = 1.

    The principal square root of Uij*Uji continues the root near +1 on
    sufficiently short edges. This product is invariant under arbitrary
    nonzero complex rescaling of each right vector (with the dual left
    rescaling). Taking sqrt(Uij/Uji) instead would introduce gauge-dependent
    signs. A zero product or a product on the square-root cut is inadmissible;
    the caller must resolve that edge rather than assign it a zero phase.
    """
    i, j = np.asarray(edges, dtype=int).reshape(-1, 2).T
    forward = np.einsum("ij,ij->i", left[i].conj(), right[j])
    backward = np.einsum("ij,ij->i", left[j].conj(), right[i])
    product = forward * backward
    invalid = ~np.isfinite(product) | (np.abs(product) <= LINK_EPS**2)
    on_cut = ((product.real <= 0)
              & (np.abs(product.imag) <= LINK_BRANCH_CUT_TOL * np.abs(product)))
    if np.any(invalid | on_cut):
        n = np.flatnonzero(invalid | on_cut)[0]
        reason = "zero/nonfinite overlap product" if invalid[n] else "square-root branch cut"
        raise ValueError(f"biorthogonal edge ({i[n]}, {j[n]}): {reason}; refine/check the band")
    return forward / np.sqrt(product)


def triangle_flux_biorthogonal_complex(right: np.ndarray, left: np.ndarray,
                                       tri: np.ndarray) -> np.ndarray:
    """Return i Log(W) for each triangle, using reciprocal LR transports.

    Each undirected edge is evaluated once. Signed phase/log-magnitude sums
    implement its exact inverse on backward traversal without multiplying
    potentially very large or small loop factors. The real flux is in
    [-pi, pi]; the imaginary flux is log|W|. Band smoothness and sufficiently
    resolved triangles are still required for continuum accuracy.
    """
    tri = np.asarray(tri, dtype=int).reshape(-1, 3)
    directed = np.stack((tri, np.roll(tri, -1, axis=1)), axis=-1).reshape(-1, 2)
    edges, inverse = np.unique(np.sort(directed, axis=1), axis=0, return_inverse=True)
    links = biorthogonal_edge_links(right, left, edges)
    signs = np.where(directed[:, 0] < directed[:, 1], 1, -1)
    phase = (signs * np.angle(links)[inverse]).reshape(-1, 3).sum(axis=1)
    log_abs = (signs * np.log(np.abs(links))[inverse]).reshape(-1, 3).sum(axis=1)
    return -np.angle(np.exp(1j * phase)) + 1j * log_abs


def triangle_flux_biorthogonal(right: np.ndarray, left: np.ndarray,
                               tri: np.ndarray) -> tuple[np.ndarray, int]:
    """Real LR flux and bad-link count; undefined transports raise ValueError."""
    flux = triangle_flux_biorthogonal_complex(right, left, tri)
    return flux.real, 0


def integrate_chern(verts, triangles, E, mu1, mu2, Hfun, *, rcond=None,
                    mode: str = "right", orientation: str = "positive",
                    progress: int = 0):
    """Return ``(ChernResult, flux_data)`` from in-memory mesh arrays.

    Hfun(beta1, beta2) returns the dense Bloch Hamiltonian. E and mu1/mu2
    are vertex arrays; verts contains real phase angles (n, 2), and triangles
    contains integer vertex ids (m, 3). flux_data holds flux, flux_complex,
    link_rule, oriented_triangles, beta1, and beta2, ready for caller-owned
    persistence or plotting. The Chern sign uses A = i<L|dR>.
    """
    rcond = DEFAULT_RCOND if rcond is None else rcond
    if orientation not in ("positive", "negative", "keep"):
        raise ValueError(f"unknown orientation {orientation!r}")
    verts = np.asarray(verts, dtype=float)
    triangles = np.asarray(triangles, dtype=int)
    E = np.asarray(E, dtype=complex)
    mu1, mu2 = np.asarray(mu1, dtype=float), np.asarray(mu2, dtype=float)
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
        complex_flux = flux.astype(complex)
        link_rule = "right-overlap"
    elif mode == "biorthogonal":
        right, left, diag = left_right_eigenvectors_on_mesh(
            Hfun, E, beta1, beta2, rcond, progress=progress
        )
        complex_flux = triangle_flux_biorthogonal_complex(right, left, oriented)
        flux, bad = complex_flux.real, 0
        link_rule = "biorthogonal-reciprocal-v1"
    else:
        raise ValueError(f"unknown mode {mode!r}")

    total_flux = float(np.sum(flux))
    result = ChernResult(
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
        link_rule=link_rule,
        total_flux_imag=float(np.sum(complex_flux.imag)),
    )

    return result, dict(flux=flux, flux_complex=complex_flux, link_rule=link_rule,
                        oriented_triangles=oriented, beta1=beta1, beta2=beta2)
