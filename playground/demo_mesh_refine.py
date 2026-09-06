'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-27
Copyright © Department of Physics, Tsinghua University. All rights reserved

EXPERIMENTAL — adaptive torus-mesh refinement by EDGE LENGTH (prototype).

Architecture (per discussion, 2026-08-27):

* Selection criterion is the torus edge length of the current
  triangulation (triangle AREA is blind to long thin slivers).

* Two-level loop.  The INNER loop is pure geometry + prediction: dyadic
  midpoints are inserted on over-threshold edges (positions = torus short-
  arc midpoints, E = value of a Clough-Tocher interpolant anchored on the
  REAL vertices only — predicted values never feed the predictor), the
  cloud is re-triangulated, and this repeats until every edge satisfies
  the threshold.  Solving happens only AFTER the inner loop converges, in
  one batch.

* Boundary-stall fallback: when a predicted E solves out of spectrum, or
  its subsets are all rejected (an index-transition location — e.g. the
  midpoint of an edge joining two (12,0) points landing in a (6,0)
  region), the prediction is connected to EACH real edge anchor's E and
  the segment is bisected (predicate: index equals the anchor's index)
  until the boundary point is located; the boundary-adjacent probe's
  points are Hungarian-matched against the failing midpoint and merged
  when acceptable.  Sides whose prediction already carries the anchor's
  index are skipped (a pure match-distance failure is not a boundary).
'''

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator
from scipy.optimize import linear_sum_assignment

from pygbz2d.core import TWO_PI
from pygbz2d.experimental import flatten_results
from demo_torus_mesh import periodic_delaunay

#: Default acceptance radius (rad-equivalent chord) for prediction ->
#: solved-point Hungarian matches.
MATCH_TOL = 0.5

#: Bisection iterations for the boundary fallback (precision = initial
#: bracket width / 2**iters).
BISECT_ITERS = 6


# ---------------------------------------------------------------------------
# Torus geometry helpers
# ---------------------------------------------------------------------------

def torus_midpoint(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Midpoint of two torus points along the short way per coordinate."""
    d = q - p
    d = (d + np.pi) % TWO_PI - np.pi
    return (p + 0.5 * d) % TWO_PI


def torus_pair_dist(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Torus arc distance between two (n, 2) point sets."""
    d = np.abs(a - b) % TWO_PI
    d = np.minimum(d, TWO_PI - d)
    return np.hypot(d[:, 0], d[:, 1])


def edge_lengths(verts: np.ndarray, tri: np.ndarray):
    """(edges, lengths) of all triangle edges (with duplicates)."""
    e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    return e, torus_pair_dist(verts[e[:, 0]], verts[e[:, 1]])


def triangle_areas(verts: np.ndarray, tri: np.ndarray) -> np.ndarray:
    """Unsigned areas with per-triangle periodic unwrapping (QC only)."""
    a = verts[tri[:, 0]]
    b = verts[tri[:, 1]] - np.round((verts[tri[:, 1]] - a) / TWO_PI) * TWO_PI
    c = verts[tri[:, 2]] - np.round((verts[tri[:, 2]] - a) / TWO_PI) * TWO_PI
    cross = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) \
        - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    return 0.5 * np.abs(cross)


def edge_len_percentiles(verts, tri):
    e, L = edge_lengths(verts, tri)
    return np.percentile(L, [50, 90, 99, 100])


# ---------------------------------------------------------------------------
# E(theta1, theta2) prediction, seam-correct via 3x3 replication
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Inner loop: pure-geometry dyadic refinement (no solver calls)
# ---------------------------------------------------------------------------

def dyadic_geometry_refine(verts: np.ndarray, tri: np.ndarray,
                           edge_thresh: float, interp,
                           max_levels: int = 8):
    """Insert dyadic midpoints on over-threshold edges until clean.

    *verts* holds the REAL vertices; pending midpoints join the cloud for
    re-triangulation but never feed the predictor.  Each pending remembers
    its two REAL anchors (inherited through nesting) for the boundary
    fallback.

    Returns (positions, predicted E, anchor_left ids, anchor_right ids,
    n_levels used).
    """
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


# ---------------------------------------------------------------------------
# Outer stage: batch solve + boundary-bisection fallback
# ---------------------------------------------------------------------------

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
                method: str, match_tol, bisect_iters=BISECT_ITERS,
                n_procs: int = 1):
    """Solve every pending midpoint; failures trigger the boundary fallback
    toward each real anchor.

    ``anchor_idx`` carries per-real-vertex (n_0D, n_1D) phase labels.
    ``n_procs`` > 1 fans the per-pending work out over a multiprocessing
    Pool (each unit = 1 direct solve + up to 2 x bisect_iters probes).
    Returns (accepted_theta, accepted_E, accepted_mu, accepted_idx, stats).
    """
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
                if (k + 1) % 25 == 0:
                    print(f"    solved {k + 1}/{len(tasks)}", flush=True)
    else:
        for k, t in enumerate(tasks):
            results.append(_solve_one(t))
            if (k + 1) % 25 == 0:
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
