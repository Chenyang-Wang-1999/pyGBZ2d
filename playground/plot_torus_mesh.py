'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-27
Copyright © Department of Physics, Tsinghua University. All rights reserved

Torus-stack visualization of band-cluster meshes (pyvista, EXPERIMENTAL).

Each cluster of the band clustering is meshed on the (theta1, theta2)
torus (Plan-A periodic Delaunay, see demo_torus_mesh) and mapped onto a
3-D torus of radii (R, r) via the standard parametrization

    x = (R + r cos theta2) cos theta1
    y = (R + r cos theta2) sin theta1
    z = r sin theta2 + z_offset(cluster)

The cluster index shifts the torus along z, so several clusters stack
vertically and their MESH QUALITY (edge density, holes, slivers) can be
compared in one figure.  Seam-crossing triangles need NO special casing
here: cos/sin of theta and theta +- 2pi coincide, so they map to proper
small triangles on the torus (unlike a flat theta1-theta2 view, where
they streak across the plot and are display-filtered instead).

Scalars on vertices: Re E / Im E / |E| / mu1 / mu2 (shared color range
across clusters), or flat per-cluster colors with a legend.

Usage:
    python playground/plot_torus_mesh.py                          # default data
    python playground/plot_torus_mesh.py data/Hermitian-Haldane-amoeba.pkl
    python playground/plot_torus_mesh.py --scalar mu2 --flat
    python playground/plot_torus_mesh.py --screenshot out.png --off-screen
'''

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse

import numpy as np

try:
    import pyvista as pv
except ImportError:
    sys.exit("pyvista is required: pip install pyvista")

from pygbz2d.core import TWO_PI
from pygbz2d.experimental import cluster_bands, summarize_clusters
from demo_torus_mesh import build_cluster_mesh

DEFAULT_DATA = "data/Haldane-gain-loss-amoeba.pkl"

TORUS_R = 2.0    # major radius
TORUS_r = 0.6    # minor radius
STACK_SPACING = 3.0   # z-gap between clusters, in units of the minor radius

SCALAR_CHOICES = ("ReE", "ImE", "absE", "mu1", "mu2", "cluster")


def torus_xyz(theta: np.ndarray, R: float, r: float) -> np.ndarray:
    """(n, 2) torus angles -> (n, 3) points on the embedded torus."""
    t1, t2 = theta[:, 0], theta[:, 1]
    return np.stack([(R + r * np.cos(t2)) * np.cos(t1),
                     (R + r * np.cos(t2)) * np.sin(t1),
                     r * np.sin(t2)], axis=1)


def flat_xyz(theta: np.ndarray, z: float) -> np.ndarray:
    """(n, 2) torus angles -> flat development at height z."""
    return np.stack([theta[:, 0], theta[:, 1],
                     np.full(len(theta), z)], axis=1)


def seam_display_filter(verts: np.ndarray, tri: np.ndarray) -> np.ndarray:
    """DISPLAY-only: drop triangles whose raw-coordinate edge exceeds pi.

    Only needed in --flat mode, where seam triangles span the [0, 2pi)
    square; the torus embedding needs no such filter (periodicity is
    absorbed by cos/sin).
    """
    e = verts[tri]
    long = (np.abs(e[:, [1, 2, 0]] - e[:, [0, 1, 2]]).max(axis=(1, 2))
            > np.pi)
    return tri[~long] if long.any() else tri


def cluster_scalar(verts_data: dict, name: str) -> np.ndarray:
    if name == "ReE":
        return verts_data["E"].real
    if name == "ImE":
        return verts_data["E"].imag
    if name == "absE":
        return np.abs(verts_data["E"])
    return verts_data[name]


def plt_get_cmap(n):
    """Distinct solid colors without pulling matplotlib's pyplot state."""
    import matplotlib.colors as mcolors
    base = list(mcolors.TABLEAU_COLORS.values())
    return lambda i: base[i % len(base)]


def main(argv=None):
    p = argparse.ArgumentParser(description="torus-stack mesh viewer")
    p.add_argument("data", nargs="?", default=DEFAULT_DATA,
                   help="sweep pickle (ignored with --mesh)")
    p.add_argument("--mesh", nargs="+", default=None,
                   help="mesh .npz file(s) saved by demo_pipeline.py "
                        "(refined meshes; skips clustering entirely)")
    p.add_argument("--eps", type=float, default=None,
                   help="clustering radius (default: auto plateau)")
    p.add_argument("--scalar", choices=SCALAR_CHOICES, default="ReE")
    p.add_argument("--min-size", type=int, default=100,
                   help="clusters below this size are not plotted")
    p.add_argument("--max-clusters", type=int, default=4)
    p.add_argument("--flat", action="store_true",
                   help="flat (theta1, theta2) development instead of a torus")
    p.add_argument("--R", type=float, default=TORUS_R)
    p.add_argument("--r", type=float, default=TORUS_r)
    p.add_argument("--spacing", type=float, default=STACK_SPACING,
                   help="z-gap between clusters (torus: units of r; flat: rad)")
    p.add_argument("--no-edges", action="store_true",
                   help="hide triangle edges (quality inspection wants them ON)")
    p.add_argument("--off-screen", action="store_true")
    p.add_argument("--screenshot", type=str, default=None)
    args = p.parse_args(argv)

    if args.mesh:
        # refined meshes straight from the pipeline (no clustering rerun)
        entries = []
        for f in args.mesh:
            z = np.load(f)
            vdata = {"E": z["E"], "mu1": z["mu1"], "mu2": z["mu2"]}
            entries.append((str(f), z["verts"], z["triangles"].astype(int),
                            vdata))
        print(f"meshes: {[e[0] for e in entries]}")
    else:
        import pickle
        with open(args.data, "rb") as fp:
            data = pickle.load(fp)
        print(f"data: {args.data}")
        # E-block steps from the ORIGINAL grid axes when available —
        # enriched sweeps have probe energies between grid nodes that
        # would dilute the inferred median step and fragment the
        # clustering (same fix as demo_pipeline.py).
        d_re = d_im = None
        if "E_real" in data and "E_imag" in data:
            re_u = np.unique(data["E_real"])
            im_u = np.unique(data["E_imag"])
            if len(re_u) > 1:
                d_re = float(np.median(np.diff(re_u)))
            if len(im_u) > 1:
                d_im = float(np.median(np.diff(im_u)))
        cl = cluster_bands(data["results"], eps=args.eps,
                           d_re=d_re, d_im=d_im)
        stats = summarize_clusters(cl.points, cl.labels)
        clusters = [s for s in stats if s["size"] >= args.min_size]
        clusters = clusters[:args.max_clusters]
        print(f"{cl.n_clusters} clusters @ eps={cl.eps:.4g}; plotting "
              f"{[s['label'] for s in clusters]}")
        entries = [(f"cluster {s['label']}", None, None, None)
                   for s in clusters]

    meshes = []
    for rank, (label, verts_in, tri_in, vdata_in) in enumerate(entries):
        if args.mesh:
            verts, tri, vdata = verts_in, tri_in, vdata_in
        else:
            verts, tri, vdata, _ = build_cluster_mesh(cl, int(label.split()[-1]))
        scal = (cluster_scalar(vdata, args.scalar)
                if args.scalar != "cluster" else None)
        if args.flat:
            tri = seam_display_filter(verts, tri)
            pts = flat_xyz(verts, rank * args.spacing)
        else:
            pts = torus_xyz(verts, args.R, args.r)
            pts[:, 2] += rank * args.spacing * args.r
        mesh = pv.PolyData(pts)
        mesh.faces = (np.column_stack(
            [np.full(len(tri), 3), tri]).ravel().astype(np.int64))
        if args.scalar != "cluster":
            mesh.point_data[args.scalar] = scal
        meshes.append((label, mesh, scal))
        print(f"  {label}: V={len(verts)} F={len(tri)} "
              f"(z-offset rank {rank})")

    pl = pv.Plotter(off_screen=args.off_screen, window_size=(1400, 1000))
    pl.add_axes()

    if args.scalar == "cluster":
        colors = plt_get_cmap(len(meshes))
        for rank, (label, mesh, _) in enumerate(meshes):
            pl.add_mesh(mesh, color=colors(rank), show_edges=not args.no_edges,
                        edge_color="black", label=label,
                        smooth_shading=False)
        pl.add_legend()
    else:
        vals = np.concatenate([sc for _, _, sc in meshes])
        clim = (float(np.min(vals)), float(np.max(vals)))
        for k, (label, mesh, _) in enumerate(meshes):
            pl.add_mesh(mesh, scalars=args.scalar, clim=clim,
                        cmap="viridis", show_edges=not args.no_edges,
                        edge_color="black", smooth_shading=False,
                        show_scalar_bar=(k == len(meshes) - 1),
                        scalar_bar_args={"title": args.scalar})

    if args.screenshot:
        pl.show(screenshot=args.screenshot)
        print(f"screenshot -> {args.screenshot}")
    else:
        pl.show()


if __name__ == "__main__":
    main()
