'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-26
Copyright © Department of Physics, Tsinghua University. All rights reserved

Demo 2 (FKU): Fuchs–Kedem–Uselton optimal contour triangulation.

Reference: Fuchs, Kedem, Uselton 1977, "Optimal Surface Reconstruction from
Planar Contours".  Replaces the greedy Hultquist advancing front
(demo_advance_front) with a globally-optimal triangulation of the gap between
two matched LineSubsets, operating purely on the curves' original sample
points — no resampling, no new points, every output vertex is an original
LineSubset sample that exactly satisfies f(E, beta1, beta2) = 0.

Algorithm
---------
A triangulation between two open polylines A (m vertices) and B (n vertices)
corresponds to a monotone path on the (m-1) x (n-1) grid of *quads*:
each quad Q[i,j] = (A_i, A_{i+1}, B_{j+1}, B_j) is split by one of its two
diagonals.  Splitting along A_i—B_{j+1} yields the two triangles
  (A_i, A_{i+1}, B_{j+1})  and  (A_i, B_{j+1}, B_j)
and is reached from Q[i-1,j] or Q[i,j-1].  Choosing a diagonal for every quad
on a monotone path from Q[0,0] to Q[m-2,n-2] tiles the strip with exactly
2*(m+n-2)-2 ... actually 2*( (m-1)+(n-1) ) triangles = 2*(m+n-2) - overlap,
independent of the path — the path only decides *which* triangles.

We solve the optimal split by DP on the quad grid:
  cost[i,j] = cost of best monotone path Q[0,0] -> Q[i,j]
  transition from (i-1,j) or (i,j-1), adding the two triangles of the chosen
  diagonal of Q[i,j].  Arc weight = sum of the two triangle perimeters
  (geometrically neutral; minimising total edge length avoids skinny slivers).

The four anchor vertices (theta1 endpoints, beta2 ~ 1, where the GBZ closes)
are NOT special-cased here: their near-coincident beta2 contributes near-zero
cost to anchor triangles, which is the correct behaviour — the optimal path is
free to route through them cheaply, exactly as the geometry dictates.

Comparison with demo_advance_front:
  - Hultquist: greedy, O(m+n), local choice per step.
  - FKU: globally optimal monotone path, O(m*n) DP.  No leftover-fan needed
    because the quad-grid path covers both curves fully and symmetrically.
Both are interpolation-free; FKU gives the minimum-total-perimeter strip.
'''

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from pygbz2d.core import LineSubset, to_sphere_r3


def _oriented(line: LineSubset):
    """Return (theta1, beta2) with theta1 ascending (consistent front anchor)."""
    th = np.asarray(line.theta1_arr, dtype=float)
    b2 = np.asarray(line.beta2_arr, dtype=complex)
    if len(th) >= 2 and th[0] > th[-1]:
        th = th[::-1]
        b2 = b2[::-1]
    return th, b2


def _info(verts, tris, m, n, base, total_cost, th_a, th_b, D=None):
    """Build the diagnostics dict, recomputing the distance matrix if needed."""
    if D is None:
        # diagnostics-only fallback: beta2 chordal (degenerate fan cases)
        all_b2 = np.array([v[2] for v in verts], dtype=complex)
        r2 = to_sphere_r3(all_b2)
        diff = r2[:, None] - r2[None]
        D = np.sqrt(np.sum(diff * diff, axis=2))
    anchor_idx = {0, m - 1, base, base + n - 1}
    e01 = D[tris[:, 0], tris[:, 1]]
    e12 = D[tris[:, 1], tris[:, 2]]
    e02 = D[tris[:, 0], tris[:, 2]]
    all_edges = np.concatenate([e01, e12, e02])
    tri_anchor = np.array([bool(set(int(x) for x in t) & anchor_idx) for t in tris])
    edge_anchor = np.concatenate([tri_anchor, tri_anchor, tri_anchor])
    interior_edges = all_edges[~edge_anchor]
    anchor_edges = all_edges[edge_anchor]
    return {
        "m": m, "n": n,
        "n_tri": len(tris),
        "total_cost": float(total_cost),
        "min_interior_edge": float(interior_edges.min()) if interior_edges.size else float("nan"),
        "max_edge": float(all_edges.max()) if all_edges.size else float("nan"),
        "min_anchor_edge": float(anchor_edges.min()) if anchor_edges.size else float("nan"),
        "th1_a_range": (float(th_a[0]), float(th_a[-1])),
        "th1_b_range": (float(th_b[0]), float(th_b[-1])),
    }


def _pairwise_chordal_table(b1_a, b2_a, b1_b, b2_b):
    """Full (beta1, beta2) chordal-distance table between every pair of vertices.

    Distance(i, j) = chordal(beta1_i, beta1_j) + chordal(beta2_i, beta2_j),
    i.e. the sum of the per-component Riemann-sphere distances — the same
    metric used by gbz_types.chordal_cost_matrix for PointSubset matching.

    Including beta1 (whose phase is theta1) in the cost prevents FKU from
    running many consecutive advances on one side: a long run moves theta1
    far on that side, making the beta1-chordal term large, so the DP balances
    advances between A and B and the strip stays uniform in theta1, not just
    in beta2.

    Globals are [A_0..A_{m-1}, B_0..B_{n-1}].
    """
    all_b1 = np.concatenate([b1_a, b1_b])
    all_b2 = np.concatenate([b2_a, b2_b])
    r1 = to_sphere_r3(all_b1)
    r2 = to_sphere_r3(all_b2)
    diff1 = r1[:, None, :] - r1[None, :, :]
    diff2 = r2[:, None, :] - r2[None, :, :]
    D = np.sqrt(np.sum(diff1 * diff1, axis=2)) + np.sqrt(np.sum(diff2 * diff2, axis=2))
    return D


def fku_triangulate(La: LineSubset, Lb: LineSubset):
    """Fuchs–Kedem–Uselton optimal triangulation of two matched LineSubsets.

    Models the strip as a shortest-path / DP on the (m x n) vertex-pair grid.
    A monotone path from (0,0) to (m-1, n-1) — advancing A or B one vertex at a
    time — produces exactly m+n-2 triangles, one per step:
      advance A (i -> i+1): triangle (A_i, A_{i+1}, B_j)
      advance B (j -> j+1): triangle (A_i, B_j, B_{j+1})
    This is the same strip topology as the Hultquist advancing front
    (demo_advance_front): exactly m+n-2 triangles, no fan, fully covering both
    curves.  FKU picks the *globally minimal* total-triangle-perimeter path
    rather than a greedy per-step choice.

    Vertices A and B are concatenated; tris reference globals
    verts[:m] = A, verts[m:m+n] = B.  All output vertices are original
    LineSubset samples (interpolation-free).

    Returns:
        verts : (m+n, 3) object ndarray — (slice_tag, theta1, beta2).
        tris  : (m+n-2, 3) int ndarray.
        info  : diagnostics dict.
    """
    th_a, b2_a = _oriented(La)
    th_b, b2_b = _oriented(Lb)
    m, n = len(th_a), len(th_b)
    base = m

    # beta1 = exp(mu1 + i*theta1) carries theta1 as its phase.  Including it in
    # the cost metric makes FKU balance advances in theta1, not just beta2.
    b1_a = np.exp(La.mu1 + 1j * th_a)
    b1_b = np.exp(Lb.mu1 + 1j * th_b)

    verts = np.empty((m + n, 3), dtype=object)
    for i in range(m):
        verts[i] = (0, float(th_a[i]), complex(b2_a[i]))
    for j in range(n):
        verts[base + j] = (1, float(th_b[j]), complex(b2_b[j]))

    # degenerate cases (one curve is a single vertex): fan
    if m == 1 and n >= 1:
        tris = [(0, base + j, base + j + 1) for j in range(n - 1)]
        tris = np.asarray(tris, dtype=int)
        return verts, tris, _info(verts, tris, m, n, base, 0.0, th_a, th_b,
                                  D=_pairwise_chordal_table(b1_a, b2_a, b1_b, b2_b))
    if n == 1 and m >= 1:
        tris = [(i, i + 1, base) for i in range(m - 1)]
        tris = np.asarray(tris, dtype=int)
        return verts, tris, _info(verts, tris, m, n, base, 0.0, th_a, th_b,
                                  D=_pairwise_chordal_table(b1_a, b2_a, b1_b, b2_b))

    # Precompute full (beta1, beta2) chordal distance matrix (vectorised) once.
    D = _pairwise_chordal_table(b1_a, b2_a, b1_b, b2_b)

    # DP on the (m x n) grid of vertex pairs (i, j).
    # cost[i,j] = min total triangle perimeter along a monotone path (0,0)->(i,j).
    # transition into (i,j):
    #   from (i-1,j): advance A, triangle (A_{i-1}, A_i, B_j), cost add D[i-1,i]+D[i-1,bj]+D[i,bj]
    #   from (i,j-1): advance B, triangle (A_i, B_{j-1}, B_j), cost add D[i,bj-1]+D[i,bj]+D[bj-1,bj]
    INF = np.inf
    cost = np.full((m, n), INF)
    back = np.full((m, n), -1, dtype=np.int8)  # -1=seed, 0=came from left (j-1), 1=from above (i-1)
    cost[0, 0] = 0.0
    for i in range(m):
        for j in range(n):
            if i == 0 and j == 0:
                continue
            bj = base + j
            c_up = INF
            c_left = INF
            if i > 0:
                bm1 = i - 1
                tri = D[bm1, i] + D[bm1, bj] + D[i, bj]
                c_up = cost[i - 1, j] + tri
            if j > 0:
                bjm1 = base + j - 1
                tri = D[i, bjm1] + D[i, bj] + D[bjm1, bj]
                c_left = cost[i, j - 1] + tri
            if c_up <= c_left:
                cost[i, j] = c_up
                back[i, j] = 1
            else:
                cost[i, j] = c_left
                back[i, j] = 0

    # backtrack: each step records the triangle that was added entering that node
    tris = []
    i, j = m - 1, n - 1
    while not (i == 0 and j == 0):
        bj = base + j
        if back[i, j] == 1:          # came from (i-1, j): advance A
            tris.append((i - 1, i, bj))
            i -= 1
        else:                         # came from (i, j-1): advance B
            tris.append((i, bj - 1, bj))
            j -= 1
    tris.reverse()
    tris = np.asarray(tris, dtype=int)

    return verts, tris, _info(verts, tris, m, n, base,
                              float(cost[m - 1, n - 1]), th_a, th_b, D=D)


