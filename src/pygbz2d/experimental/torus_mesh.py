# Copyright © Department of Physics, Tsinghua University. All rights reserved.
# Numerical routines extracted from the playground prototypes by wangchenyang.

"""Experimental periodic band meshes and geometry diagnostics (NumPy/SciPy only).

Vertices are real phase coordinates in radians. Delaunay triangulation assumes
one sheet over the phase torus; topology checks alone cannot certify that a
sparse or multivalued GBZ has been sampled adequately. No model, plotting,
filesystem, or command-line policy belongs here.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import Delaunay, cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from ..core import TWO_PI

# Endpoints obtained by bisection can differ by about 1e-6 radians. Merging
# below 1e-5 prevents Qhull from silently excluding near-coincident vertices.
DEDUP_TOL = 1e-5

def dedup_vertices(theta1: np.ndarray, theta2: np.ndarray,
                   data: dict | None = None, tol: float | None = None):
    """Merge near-coincident projected vertices; first occurrence wins.

    Uses a kd-tree tolerance graph + connected components — a rounding-grid
    key would miss pairs straddling a grid-cell boundary.
    """
    if tol is None:
        tol = DEDUP_TOL
    key = np.stack([theta1, theta2], axis=1)
    n = len(key)
    if n > 1:
        pairs = cKDTree(key).query_pairs(tol, output_type="ndarray")
    else:
        pairs = np.empty((0, 2), dtype=int)
    if len(pairs):
        graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])),
                           shape=(n, n))
        _, comp = connected_components(graph, directed=False)
    else:
        comp = np.arange(n)
    _, first_idx = np.unique(comp, return_index=True)
    n_merged = n - len(first_idx)
    t1, t2 = key[first_idx, 0], key[first_idx, 1]
    out_data = None
    if data is not None:
        out_data = {k: v[first_idx] for k, v in data.items()}
    return t1, t2, out_data, n_merged


def periodic_delaunay(vertices: np.ndarray) -> np.ndarray:
    """Triangle vertex-index array of the torus-consistent Delaunay mesh.

    Periodicity is enforced BY CONSTRUCTION (the seam keep-rule bug was
    traced to translate copies differing at ULP level, letting Qhull break
    near-degenerate mirrored configurations differently on the two sides
    of a seam):

    * Equivalence bookkeeping is INDEX-based throughout: with n central
      points, replicated point j lives in tile ``j // n`` and is original
      point ``j % n`` (``owner`` below).  Only the (0,0) tile's decisions
      are kept; triangle dedupe works on sorted original-id triples.
    * Coordinates handed to Qhull are made EXACTLY periodic: rescaled to
      period 1 (exactly representable) and snapped to a 2^-32 grid, so
      replicas are the same float value plus an INTEGER offset — bit-exact
      translates, hence bit-identical predicate inputs in every tile.
    * An index-keyed perturbation (original id only, max 2^-32) breaks
      near-cospherical mirrored configurations the SAME way in every tile
      (translation-covariant), removing the degenerate ties whose
      inconsistent resolution caused over-stuffed seam edges.

    Snap/perturbation are ~2e-10 of the period — three orders below the
    1e-5 dedup tolerance, invisible to the geometry; vertex POSITIONS in
    the returned mesh stay the original thetas.
    """
    n = len(vertices)
    GRID = 2.0 ** 32
    i_orig = np.arange(n)
    # Knuth multiplicative hashes: original-index-keyed, tile-independent
    eps1 = ((i_orig * 2654435761) % 4096) * 2.0 ** -44
    eps2 = ((i_orig * 2246822519) % 4096) * 2.0 ** -44
    u = (np.round(vertices / TWO_PI * GRID) / GRID
         + np.stack([eps1, eps2], axis=1))

    # (0, 0) FIRST; offsets are exact integers so replicas are bit-exact
    # translates of the perturbed canonical cloud
    offsets = np.array(
        [[0, 0]] + [[di, dj] for di in (-1, 0, 1) for dj in (-1, 0, 1)
                    if (di, dj) != (0, 0)], dtype=float)
    reps = (u[None, :, :] + offsets[:, None, :]).reshape(-1, 2)

    j = np.arange(9 * n)                    # replicated point index
    owner = j % n                           # original point id
    tile_of = j // n                        # tile id; (0,0) is tile 0
    central = tile_of == 0

    try:
        tri = Delaunay(reps)
    except Exception:
        # Joggle degenerate (collinear/cospherical) configurations.
        tri = Delaunay(reps, qhull_options="QJ")
    n_dropped_by_qhull = len(getattr(tri, "coplanar", []))
    if n_dropped_by_qhull:
        print(f"    [qhull] {n_dropped_by_qhull} replicated points excluded "
              f"as near-coincident (coplanar) — raise DEDUP_TOL if large")

    simplices = tri.simplices
    keep = central[simplices].any(axis=1)      # touches the central tile
    mapped = owner[simplices[keep]]

    # Degenerate after mod-mapping (two images of the same original vertex)
    # and duplicates (the same torus triangle kept from two tile copies).
    ok = ((mapped[:, 0] != mapped[:, 1]) & (mapped[:, 1] != mapped[:, 2])
          & (mapped[:, 0] != mapped[:, 2]))
    mapped = mapped[ok]
    order = np.sort(mapped, axis=1)
    _, uniq_idx = np.unique(order, axis=0, return_index=True)
    return mapped[np.sort(uniq_idx)]


def mesh_topology(triangles: np.ndarray) -> dict:
    """Euler characteristic, boundary loops and components of the mesh."""
    n_used = int(triangles.max()) + 1
    edges = np.vstack([triangles[:, [0, 1]], triangles[:, [1, 2]],
                       triangles[:, [2, 0]]])
    edges = np.sort(edges, axis=1)
    uniq, counts = np.unique(edges, axis=0, return_counts=True)
    n_edges = len(uniq)

    chi = n_used - n_edges + len(triangles)

    # boundary loops: walk UNDIRECTED edges used by exactly one triangle
    # (stored sorted, so a directed successor lookup would break chains).
    boundary = uniq[counts == 1]
    loops = []
    if len(boundary):
        incident: dict[int, list[int]] = {}
        for k, (a, b) in enumerate(boundary):
            incident.setdefault(int(a), []).append(k)
            incident.setdefault(int(b), []).append(k)
        used = np.zeros(len(boundary), dtype=bool)
        for k0 in range(len(boundary)):
            if used[k0]:
                continue
            a, b = boundary[k0]
            used[k0] = True
            loop = [int(a), int(b)]
            cur = int(b)
            while True:
                nxt = None
                for k in incident.get(cur, []):
                    if not used[k]:
                        nxt = k
                        break
                if nxt is None:
                    break
                used[nxt] = True
                e = boundary[nxt]
                cur = int(e[0]) if e[1] == cur else int(e[1])
                if cur == loop[0]:
                    break
                loop.append(cur)
            loops.append(loop)

    # connected components over the vertex graph
    adj = coo_matrix(
        (np.ones(n_edges), (uniq[:, 0], uniq[:, 1])), shape=(n_used, n_used))
    n_comp, _ = connected_components(adj, directed=False)

    return {
        "n_vertices": n_used,
        "n_edges": n_edges,
        "n_triangles": int(len(triangles)),
        "chi": int(chi),
        "n_boundary_edges": int(len(boundary)),
        "boundary_loops": loops,
        "n_components": int(n_comp),
    }


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


def build_band_mesh(points, mask=None, *, dedup_tol=None):
    """Mesh a BandPoints cloud (or a selected band); retain sample provenance.

    Returns ``(verts, triangles, vertex_data, n_merged)``. Vertex data contains
    E, mu1, mu2, and slice_idx; the first coincident sample wins, as in the
    original prototype. No clustering or energy-band selection is performed.
    """
    selected = slice(None) if mask is None else mask
    data = {name: getattr(points, name)[selected]
            for name in ("E", "mu1", "mu2", "slice_idx")}
    t1, t2, vdata, n_merged = dedup_vertices(
        points.theta1[selected], points.theta2[selected], data, tol=dedup_tol)
    verts = np.stack([t1, t2], axis=1)
    return verts, periodic_delaunay(verts), vdata, n_merged


def build_cluster_mesh(cl, cluster_id: int, *, dedup_tol=None):
    """Adapter from a BandClustering result to :func:`build_band_mesh`."""
    return build_band_mesh(cl.points, cl.labels == cluster_id, dedup_tol=dedup_tol)
