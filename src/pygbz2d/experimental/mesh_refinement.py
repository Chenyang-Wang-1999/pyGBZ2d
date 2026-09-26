# Copyright © Department of Physics, Tsinghua University. All rights reserved.
# Numerical routines extracted from the playground prototypes by wangchenyang.

"""Experimental edge-length refinement of a periodic GBZ band mesh.

Dyadic geometry predicts target energies; only matched GBZ solutions become
real vertices. An index-boundary bisection is retained as the fallback.
Algorithms and default tolerances follow the original playground prototype;
reaching the iteration cap does not imply convergence. Uses no plotting or
model-specific code. Multiprocessing workers live in this installed module.
"""
from __future__ import annotations

import time
import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator
from scipy.optimize import linear_sum_assignment
from ..core import TWO_PI
from .band_clustering import flatten_results
from .torus_mesh import (periodic_delaunay, dedup_vertices, torus_midpoint,
                         edge_lengths, edge_len_percentiles, triangle_areas)

MATCH_TOL = 0.5  # Acceptance radius in the (cos,sin) embedding.
BISECT_ITERS = 6
MAX_GEOMETRY_LEVELS = 8
EDGE_THRESH = 0.2  # Radians on the phase torus.
MAX_ITERS = 10

def build_E_interpolator(verts: np.ndarray, E: np.ndarray):
    """Clough-Tocher C1 interpolant of complex E over the torus.

    Built on the replicated cloud so queries near theta = 0 = 2pi see the
    correct neighbours (Re and Im interpolated as separate real fields).
    Anchored on REAL vertices only — the caller never feeds it predicted
    values (self-consistent interpolation drifts; its accuracy is anchored
    at the original data).
    """
    offs = np.array([[di, dj] for di in (-1, 0, 1) for dj in (-1, 0, 1)],
                    dtype=float) * TWO_PI
    pts = (verts[None, :, :] + offs[:, None, :]).reshape(-1, 2)
    val = np.tile(E, 9)
    ct_re = CloughTocher2DInterpolator(pts, val.real)
    ct_im = CloughTocher2DInterpolator(pts, val.imag)
    return lambda q: ct_re(q) + 1j * ct_im(q)


def chord4(pts_theta: np.ndarray) -> np.ndarray:
    """(n, 2) torus angles -> (n, 4) (cos, sin) x 2 embedding."""
    return np.stack([np.cos(pts_theta[:, 0]), np.sin(pts_theta[:, 0]),
                     np.cos(pts_theta[:, 1]), np.sin(pts_theta[:, 1])],
                    axis=1)


def match_result_to_target(res, target: np.ndarray, match_tol: float):
    """Hungarian-match a solver result's points against one target position.

    Returns (theta2-row, E, mu1, mu2, dist) of the best match, or None if
    the result is unusable or the best match exceeds *match_tol*.
    """
    if not (getattr(res, "success", False) and res.is_gbz and res.subsets):
        return None
    bp = flatten_results([res])
    if len(bp.E) == 0:
        return None
    solved = np.stack([bp.theta1, bp.theta2], axis=1)
    cost = np.linalg.norm(chord4(solved) - chord4(target[None, :]), axis=1)
    _, col = linear_sum_assignment(cost[None, :])
    j = int(col[0])
    if cost[j] > match_tol:
        return None
    return (solved[j], complex(res.E_ref), bp.mu1[j], bp.mu2[j],
            float(cost[j]))


def dyadic_geometry_refine(verts: np.ndarray, tri: np.ndarray,
                           edge_thresh: float, interp,
                           max_levels: int | None = None):
    """Insert dyadic midpoints on over-threshold edges until clean.

    *verts* holds the REAL vertices; pending midpoints join the cloud for
    re-triangulation but never feed the predictor.  Each pending remembers
    its two REAL anchors (inherited through nesting) for the boundary
    fallback.

    Returns (positions, predicted E, anchor_left ids, anchor_right ids,
    n_levels used).
    """
    if max_levels is None:
        max_levels = MAX_GEOMETRY_LEVELS
    if max_levels < 1:
        raise ValueError("max_levels must be positive")
    n_real = len(verts)
    real = verts
    cur_tri = tri
    pend_pos: list = []
    pend_E: list = []
    pend_aL: list = []
    pend_aR: list = []
    seen_keys: set = set()

    for level in range(1, max_levels + 1):
        cur = np.vstack([real, np.array(pend_pos)]) if pend_pos else real
        e, L = edge_lengths(cur, cur_tri)
        over = L > edge_thresh
        if not over.any():
            return (np.array(pend_pos), np.array(pend_E, dtype=complex),
                    np.array(pend_aL, dtype=int),
                    np.array(pend_aR, dtype=int), level - 1)
        added = 0
        for (u, v) in e[over]:
            key = (int(min(u, v)), int(max(u, v)))
            if key in seen_keys:
                continue
            seen_keys.add(key)
            pos = torus_midpoint(cur[u], cur[v])
            pend_pos.append(pos)
            pend_E.append(complex(interp(pos[None, :])[0]))
            # anchors: real vertex id, or inherited through pending parents
            pend_aL.append(int(u) if u < n_real else pend_aL[u - n_real])
            pend_aR.append(int(v) if v < n_real else pend_aR[v - n_real])
            added += 1
        if added == 0:
            break
        cur = np.vstack([real, np.array(pend_pos)])
        cur_tri = periodic_delaunay(cur)
    return (np.array(pend_pos), np.array(pend_E, dtype=complex),
            np.array(pend_aL, dtype=int), np.array(pend_aR, dtype=int),
            min(level, max_levels))


def bisect_to_boundary(E_in, idx_in, E_out, solver, coeffs, degs, iters):
    """Bisect [E_in, E_out] between phase labels idx_in and 'other'.

    ``lo`` always carries the anchor phase; a failed probe counts as the
    'other' side.  Returns (probes, lo, lo_res) where ``lo_res`` is the
    last probe RESULT carrying idx_in (None if only E_in itself did).
    """
    probes = []
    lo, hi = complex(E_in), complex(E_out)
    lo_res = None
    for _ in range(iters):
        Em = 0.5 * (lo + hi)
        res = solver(coeffs, degs, complex(Em))
        probes.append((complex(Em), res))
        if res.success and tuple(res.index) == tuple(idx_in):
            lo, lo_res = Em, res
        else:
            hi = Em
    return probes, lo, lo_res


def _solve_one(task):
    """Worker unit: one pending midpoint -> acceptance record.

    Module-level and taking only picklable args so it can run under
    multiprocessing on spawn platforms (the solver is imported inside by
    module-name string).  Used by BOTH the serial and parallel paths of
    :func:`solve_batch` — one code path, no behavioural divergence.
    """
    (pos, E_pred, anc_L, anc_R, coeffs, degs, method,
     match_tol, bisect_iters) = task
    import importlib
    solver = importlib.import_module(
        "pygbz2d." + method).collect_GBZ_subsets

    res = solver(coeffs, degs, complex(E_pred))
    m = match_result_to_target(res, pos, match_tol)
    if m is not None:
        return {"kind": "direct", "theta": m[0], "E": m[1], "mu": m[2:4],
                "idx": tuple(res.index), "dist": m[4], "n_probes": 0}

    idx_pred = tuple(res.index) if res.success else None
    n_probes = 0
    for E_a, idx_a in (anc_L, anc_R):
        # a side whose prediction already carries the anchor phase has
        # no index transition to bisect (pure match failure)
        if idx_pred is not None and idx_pred == tuple(idx_a):
            continue
        probes, lo, lo_res = bisect_to_boundary(
            E_a, idx_a, E_pred, solver, coeffs, degs, bisect_iters)
        n_probes += len(probes)
        if lo_res is None:
            continue
        m = match_result_to_target(lo_res, pos, match_tol)
        if m is not None:
            return {"kind": "boundary", "theta": m[0], "E": lo,
                    "mu": m[2:4], "idx": tuple(lo_res.index),
                    "dist": m[4], "n_probes": n_probes}
    return {"kind": "rejected", "theta": None, "n_probes": n_probes}


def solve_batch(coeffs, degs, pend_pos, pend_E, pend_aL, pend_aR,
                anchor_E: np.ndarray, anchor_idx: np.ndarray,
                method: str, match_tol, bisect_iters=None,
                n_procs: int = 1, *, verbose: bool = False):
    """Solve every pending midpoint; failures trigger the boundary fallback
    toward each real anchor.

    ``anchor_idx`` carries per-real-vertex (n_0D, n_1D) phase labels.
    ``n_procs`` > 1 fans the per-pending work out over a multiprocessing
    Pool (each unit = 1 direct solve + up to 2 x bisect_iters probes).
    Returns (accepted_theta, accepted_E, accepted_mu, accepted_idx, stats).
    """
    if bisect_iters is None:
        bisect_iters = BISECT_ITERS
    tasks = [(pend_pos[i], complex(pend_E[i]),
              (complex(anchor_E[int(pend_aL[i])]),
               tuple(anchor_idx[int(pend_aL[i])])),
              (complex(anchor_E[int(pend_aR[i])]),
               tuple(anchor_idx[int(pend_aR[i])])),
              coeffs, degs, method, match_tol, bisect_iters)
             for i in range(len(pend_E))]

    results = []
    if n_procs and n_procs > 1:
        import multiprocessing as mp
        with mp.Pool(n_procs) as pool:
            for k, rec in enumerate(pool.imap(_solve_one, tasks,
                                              chunksize=4)):
                results.append(rec)
                if verbose and (k + 1) % 25 == 0:
                    print(f"    solved {k + 1}/{len(tasks)}", flush=True)
    else:
        for k, t in enumerate(tasks):
            results.append(_solve_one(t))
            if verbose and (k + 1) % 25 == 0:
                print(f"    solved {k + 1}/{len(tasks)}", flush=True)

    theta_out, E_out, mu_out, idx_out = [], [], [], []
    n_direct = n_boundary = n_rejected = n_probes = 0
    match_dists, boundary_dists = [], []
    for rec in results:
        n_probes += rec["n_probes"]
        if rec["kind"] == "rejected":
            n_rejected += 1
            continue
        theta_out.append(rec["theta"])
        E_out.append(rec["E"])
        mu_out.append(rec["mu"])
        idx_out.append(rec["idx"])
        if rec["kind"] == "direct":
            n_direct += 1
            match_dists.append(rec["dist"])
        else:
            n_boundary += 1
            boundary_dists.append(rec["dist"])

    stats = {
        "n": len(pend_E), "n_direct": n_direct, "n_boundary": n_boundary,
        "n_rejected": n_rejected, "n_probes": n_probes,
        "match_p50": float(np.median(match_dists)) if match_dists else float("nan"),
        "boundary_p50": float(np.median(boundary_dists)) if boundary_dists else float("nan"),
    }
    return (np.array(theta_out), np.array(E_out, dtype=complex),
            np.array(mu_out) if mu_out else np.empty((0, 2)),
            np.array(idx_out, dtype=int).reshape(-1, 2), stats)


def refine_mesh(verts, tri, vdata, coeffs, degs, method, indices, *,
                edge_thresh=None, match_tol=None, max_iters=None,
                n_procs: int = 1, verbose: bool = False):
    """Refine an existing band mesh by solving predicted midpoint energies.

    ``indices`` is an (n_vertices, 2) array of source GBZResult.index labels.
    ``vdata`` contains E, mu1, and mu2 arrays aligned with verts. Returned
    ``(verts, triangles, vdata)`` contains only original or solved samples.
    Arrays supplied by the caller are not modified. Use edge_lengths on the
    result to check whether the requested target was reached; a stalled or
    capped iteration returns its best mesh, as in the original algorithm.
    """
    edge_thresh = EDGE_THRESH if edge_thresh is None else edge_thresh
    match_tol = MATCH_TOL if match_tol is None else match_tol
    max_iters = MAX_ITERS if max_iters is None else max_iters
    verts, tri = np.asarray(verts), np.asarray(tri, dtype=int)
    idx = np.asarray(indices, dtype=int)
    emit = print if verbose else lambda *args, **kwargs: None
    converged_at = None
    for it in range(1, max_iters + 1):
        _, L = edge_lengths(verts, tri)
        if L.max() <= edge_thresh:
            converged_at = it - 1
            break

        t0 = time.perf_counter()
        interp = build_E_interpolator(verts, vdata["E"])
        pend = dyadic_geometry_refine(verts, tri, edge_thresh, interp)
        pend_pos, pend_E, aL, aR, n_levels = pend
        t_geom = time.perf_counter() - t0
        emit(f"  [round {it}] inner geometry: {len(pend_pos)} midpoints "
              f"over {n_levels} nesting level(s) [{t_geom:.0f}s]")
        if len(pend_pos) == 0:
            emit("  stalled: no midpoints available")
            break

        new_theta, new_E, new_mu, new_idx, st = solve_batch(
            coeffs, degs, pend_pos, pend_E, aL, aR,
            vdata["E"], idx, method, match_tol, n_procs=n_procs, verbose=verbose)
        emit(f"    solved: direct={st['n_direct']} "
              f"boundary={st['n_boundary']} rejected={st['n_rejected']} "
              f"(fallback probes={st['n_probes']}, "
              f"match p50={st['match_p50']:.4g}, "
              f"boundary p50={st['boundary_p50']:.4g}) "
              f"[{time.perf_counter() - t0 - t_geom:.0f}s]")

        if len(new_theta) == 0:
            emit("  stalled: no accepted points this round")
            break

        all_theta = np.vstack([verts, new_theta])
        t1, t2, vd2, _ = dedup_vertices(
            all_theta[:, 0], all_theta[:, 1],
            {"E": np.concatenate([vdata["E"], new_E]),
             "mu1": np.concatenate([vdata["mu1"], new_mu[:, 0]]),
             "mu2": np.concatenate([vdata["mu2"], new_mu[:, 1]]),
             "idx0": np.concatenate([idx[:, 0], new_idx[:, 0]]),
             "idx1": np.concatenate([idx[:, 1], new_idx[:, 1]])})
        verts = np.stack([t1, t2], axis=1)
        vdata = vd2
        idx = np.stack([vd2["idx0"], vd2["idx1"]], axis=1)
        tri = periodic_delaunay(verts)
        _, L = edge_lengths(verts, tri)
        emit(f"    mesh: V={len(verts)} F={len(tri)} "
              f"max edge -> {L.max():.4f}")
    else:
        emit(f"  WARNING: hit max_iters={max_iters} before convergence "
              f"(max edge {edge_lengths(verts, tri)[1].max():.4f})")

    if converged_at is not None:
        emit(f"  converged after {converged_at} solving round(s) "
              f"(max edge {edge_lengths(verts, tri)[1].max():.4f} "
              f"<= {edge_thresh})")
    emit(f"  edges p50/p90/p99/max: "
          f"{np.round(edge_len_percentiles(verts, tri), 4)}")
    a = triangle_areas(verts, tri)
    emit(f"  areas p50={np.median(a):.3e} max={a.max():.3e}")
    return verts, tri, vdata
