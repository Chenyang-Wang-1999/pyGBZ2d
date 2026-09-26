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

from pygbz2d.core import TWO_PI
from pygbz2d.experimental import cluster_bands, summarize_clusters

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

DEFAULT_DATA = "data/Haldane-gain-loss-amoeba.pkl"
OUT_ROOT = Path(__file__).resolve().parent / "band_cluster_out"

from pygbz2d.experimental.torus_mesh import (
    DEDUP_TOL, dedup_vertices, periodic_delaunay, mesh_topology, build_cluster_mesh,
)


# ---------------------------------------------------------------------------
# Topology QC
# ---------------------------------------------------------------------------


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
