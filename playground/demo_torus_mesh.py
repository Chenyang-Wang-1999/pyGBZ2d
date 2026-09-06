'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-27
Copyright © Department of Physics, Tsinghua University. All rights reserved

Prototype: 2D torus triangle-mesh construction per band cluster — Plan A.

Periodic Delaunay by image replication: copy the cluster's (theta1, theta2)
vertices into the 3x3 neighbouring tiles (+-2pi shifts), run plain
scipy.spatial.Delaunay on the replicated plane, keep every simplex that
touches the central tile, map vertex ids back through mod 2pi, drop
degenerate (repeated-vertex) triangles and deduplicate.

NO alpha filter and NO full-space edge guard (that was Plan B) — this
prototype exists precisely to show where the unfiltered convex-hull
Delaunay is enough (full-torus footprints, e.g. the Hermitian Haldane
model) and where it is not (scattered footprints with holes, e.g. the
gain-loss data).

Exact-duplicate vertices are merged (numerical hygiene, not filtering):
coincident points make Qhull degenerate.  First occurrence wins; the merge
count is reported.

Topology QC per mesh: Euler characteristic chi = V - E + F, boundary-edge
count and boundary loops, connected components.  A mesh that covers the
whole torus has chi = 0 and no boundary; a disk patch has chi = 1 and one
boundary loop; an annulus has two loops.

Usage:
    python playground/demo_torus_mesh.py data/Hermitian-Haldane-amoeba.pkl
    python playground/demo_torus_mesh.py                       # gain-loss default
'''

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import argparse

import numpy as np
from scipy.spatial import Delaunay, cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from pygbz2d.core import TWO_PI
from pygbz2d.experimental import cluster_bands, summarize_clusters

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

DEFAULT_DATA = "data/Haldane-gain-loss-amoeba.pkl"
OUT_ROOT = Path(__file__).resolve().parent / "band_cluster_out"

#: Duplicate-projection merge tolerance (rad).  Spliced arc endpoints agree
#: only to ~1e-6 (bisection precision), and Qhull SILENTLY drops such
#: near-coincident points from the triangulation (they land in
#: ``Delaunay.coplanar``), leaving isolated vertices in the mesh.  1e-5 is
#: still 3.5 orders below any real feature scale (~0.05 rad sampling).
#: Merging is hygiene for Qhull, not shape filtering.
DEDUP_TOL = 1e-5


# ---------------------------------------------------------------------------
# Periodic Delaunay (Plan A)
# ---------------------------------------------------------------------------

def dedup_vertices(theta1: np.ndarray, theta2: np.ndarray,
                   data: dict | None = None, tol: float = DEDUP_TOL):
    """Merge near-coincident projected vertices; first occurrence wins.

    Uses a kd-tree tolerance graph + connected components — a rounding-grid
    key would miss pairs straddling a grid-cell boundary.
    """
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


# ---------------------------------------------------------------------------
# Topology QC
# ---------------------------------------------------------------------------

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


def describe_topology(top: dict) -> str:
    if top["n_boundary_edges"] == 0:
        if top["chi"] == 0:
            kind = "closed torus (genus 1)"
        else:
            kind = f"closed, chi={top['chi']} (NOT a torus — inspect!)"
    else:
        kind = (f"patch with boundary: {len(top['boundary_loops'])} loop(s), "
                f"chi={top['chi']}")
    return (f"V={top['n_vertices']} E={top['n_edges']} "
            f"F={top['n_triangles']} chi={top['chi']} "
            f"boundary_edges={top['n_boundary_edges']} "
            f"components={top['n_components']} -> {kind}")


# ---------------------------------------------------------------------------
# Demo pipeline
# ---------------------------------------------------------------------------

def build_cluster_mesh(cl, cluster_id: int):
    """Dedup + periodic Delaunay for one cluster of a BandClustering."""
    m = cl.labels == cluster_id
    data = {
        "E": cl.points.E[m], "mu1": cl.points.mu1[m], "mu2": cl.points.mu2[m],
        "slice_idx": cl.points.slice_idx[m],
    }
    t1, t2, vdata, n_merged = dedup_vertices(
        cl.points.theta1[m], cl.points.theta2[m], data)
    verts = np.stack([t1, t2], axis=1)
    triangles = periodic_delaunay(verts)
    return verts, triangles, vdata, n_merged


def seam_display_filter(verts: np.ndarray, triangles: np.ndarray):
    """DISPLAY-only: drop triangles whose raw-coordinate edge exceeds pi.

    Seam-crossing triangles have vertices near 0 and near 2pi in the SAME
    coordinate; on a flat [0, 2pi) plot they render as streaks across the
    whole figure.  The MESH keeps them — only the drawing hides them.  The
    3-D torus embedding needs no such filter (cos/sin absorbs periodicity).
    """
    e = verts[triangles]
    long = (np.abs(e[:, [1, 2, 0]] - e[:, [0, 1, 2]]).max(axis=(1, 2))
            > np.pi)
    return triangles[~long] if long.any() else triangles


def make_figure(tag, cluster_id, verts, triangles, top):
    OUT_ROOT.joinpath(tag).mkdir(parents=True, exist_ok=True)
    tri_plot = Triangulation(verts[:, 0], verts[:, 1],
                             seam_display_filter(verts, triangles))

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.triplot(tri_plot, color="0.7", linewidth=0.3)
    ax.scatter(verts[:, 0], verts[:, 1], s=1, c="k")
    for loop in top["boundary_loops"]:
        lv = verts[np.array(loop) % len(verts)]
        ax.plot(lv[:, 0], lv[:, 1], "r-", linewidth=1.5)
    ax.set_xlim(0, TWO_PI)
    ax.set_ylim(0, TWO_PI)
    ax.set_xlabel(r"$\theta_1$")
    ax.set_ylabel(r"$\theta_2$")
    ax.set_title(f"{tag} cluster {cluster_id}: "
                 f"chi={top['chi']}, {len(top['boundary_loops'])} loop(s)")
    fig.tight_layout()
    out = OUT_ROOT / tag / f"torus_mesh_cluster{cluster_id}.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"    figure -> {out}")


def main(argv=None):
    p = argparse.ArgumentParser(description="torus mesh prototype (Plan A)")
    p.add_argument("data", nargs="?", default=DEFAULT_DATA)
    p.add_argument("--eps", type=float, default=None,
                   help="clustering radius (default: auto plateau)")
    p.add_argument("--no-figures", action="store_true")
    args = p.parse_args(argv)

    import pickle
    with open(args.data, "rb") as fp:
        data = pickle.load(fp)
    results = data["results"]
    print(f"data: {args.data} ({len(results)} results)")

    cl = cluster_bands(results, eps=args.eps)
    print(f"clustering: {cl.n_clusters} clusters @ eps={cl.eps:.4g}")
    stats = summarize_clusters(cl.points, cl.labels)
    tag = Path(args.data).stem

    for s in [x for x in stats if x["size"] >= 100][:4]:
        cid = s["label"]
        verts, triangles, vdata, n_merged = build_cluster_mesh(cl, cid)
        top = mesh_topology(triangles)
        print(f"  cluster {cid} (size={s['size']}, merged {n_merged} dup "
              f"projections):")
        print(f"    {describe_topology(top)}")
        if not args.no_figures:
            make_figure(tag, cid, verts, triangles, top)


if __name__ == "__main__":
    main()
