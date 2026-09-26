'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-27
Copyright © Department of Physics, Tsinghua University. All rights reserved

Pipeline steps 4-6 (v2): cluster -> per-cluster torus mesh -> EDGE-LENGTH
adaptive refinement iterated to convergence (EXPERIMENTAL prototype).

Refinement architecture (see demo_mesh_refine module docstring):

  outer round:
    1. if every mesh edge <= edge_thresh -> converged;
    2. build the C1 E-predictor on REAL vertices only;
    3. INNER pure-geometry loop: dyadic midpoints on over-threshold edges
       until the (real+pending) triangulation is clean — no solver calls;
    4. batch-solve every pending midpoint; failures fall back to boundary
       bisection toward each real anchor (index-change predicate);
    5. accepted SOLVED points become real vertices, mesh re-triangulated;
  until converged, stalled (zero accepted), or max outer rounds.

Input is the ENRICHED sweep produced by enrich_boundary.py.  E-block
clustering steps come from the ORIGINAL grid axes (probe energies sit
between grid nodes and would dilute the inferred median step).

Usage:
    python playground/demo_pipeline.py data/Haldane-gain-loss-amoeba-enriched.pkl
    python playground/demo_pipeline.py --edge-thresh 0.15 --max-iters 8
'''

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse

import numpy as np

from pygbz2d.experimental import cluster_bands, summarize_clusters
from pygbz2d.experimental.torus_mesh import build_cluster_mesh
from pygbz2d.experimental.mesh_refinement import MATCH_TOL
from demo_torus_mesh import seam_display_filter

DEFAULT_DATA = "data/Haldane-gain-loss-amoeba-enriched.pkl"
EDGE_THRESH = 0.2     # absolute rad threshold on torus edge length
MAX_ITERS = 10        # outer (solving) rounds

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

OUT_ROOT = Path(__file__).resolve().parent / "band_cluster_out"


def refine_to_convergence(cl, cid, coeffs, degs, method, results_idx, *,
                          edge_thresh, match_tol, max_iters, n_procs=1):
    """Demo adapter: select a cluster, then call the installed refinement API."""
    from pygbz2d.experimental.mesh_refinement import refine_mesh
    verts, tri, vdata, _ = build_cluster_mesh(cl, cid)
    indices = np.array([results_idx[s] for s in vdata["slice_idx"]], dtype=int)
    return refine_mesh(verts, tri, vdata, coeffs, degs, method, indices,
                       edge_thresh=edge_thresh, match_tol=match_tol,
                       max_iters=max_iters, n_procs=n_procs, verbose=True)


def make_figure(tag, cid, verts, tri, vdata):
    outdir = OUT_ROOT / tag
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 6))
    # display-only seam filter: hides the 2pi-0 streaks of seam-crossing
    # triangles on the flat development (the mesh keeps them)
    ax.triplot(Triangulation(verts[:, 0], verts[:, 1],
                             seam_display_filter(verts, tri)),
               color="0.7", linewidth=0.3)
    sca = ax.scatter(verts[:, 0], verts[:, 1], s=2, c=vdata["E"].real,
                     cmap="viridis")
    plt.colorbar(sca, ax=ax, label="Re E")
    ax.set_xlim(0, 2 * np.pi)
    ax.set_ylim(0, 2 * np.pi)
    ax.set_xlabel(r"$\theta_1$")
    ax.set_ylabel(r"$\theta_2$")
    ax.set_title(f"{tag} cluster {cid}: refined mesh")
    fig.tight_layout()
    out = outdir / f"pipeline_cluster{cid}.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  figure -> {out}")


def main(argv=None):
    p = argparse.ArgumentParser(description="steps 4-6 pipeline (v2)")
    p.add_argument("data", nargs="?", default=DEFAULT_DATA)
    p.add_argument("--eps", type=float, default=None)
    p.add_argument("--edge-thresh", type=float, default=EDGE_THRESH,
                   help="absolute rad threshold on torus edge length")
    p.add_argument("--max-iters", type=int, default=MAX_ITERS,
                   help="outer (solving) round cap")
    p.add_argument("--match-tol", type=float, default=MATCH_TOL)
    p.add_argument("--min-size", type=int, default=100)
    p.add_argument("--method", choices=("amoeba", "sgbz"), default="amoeba")
    p.add_argument("--n-procs", type=int, default=1,
                   help="worker processes for the solve stage (mp.Pool; "
                        ">1 needs an environment that allows named pipes)")
    p.add_argument("--no-figures", action="store_true")
    args = p.parse_args(argv)

    import pickle
    with open(args.data, "rb") as fp:
        data = pickle.load(fp)
    coeffs, degs = data["coeffs"], data["degs"]

    info = data.get("enrich_info", {})
    print(f"data: {args.data} "
          f"(+{info.get('n_probes', 0)} boundary probes) "
          f"method={args.method}, n_procs={args.n_procs}")

    # E-block steps from the ORIGINAL grid axes (see module docstring)
    d_re = d_im = None
    if "E_real" in data and "E_imag" in data:
        re_u = np.unique(data["E_real"])
        im_u = np.unique(data["E_imag"])
        if len(re_u) > 1:
            d_re = float(np.median(np.diff(re_u)))
        if len(im_u) > 1:
            d_im = float(np.median(np.diff(im_u)))
        print(f"E-grid steps from original axes: d_re={d_re}, d_im={d_im}")

    results_idx = [tuple(r.index) if r.success else None
                   for r in data["results"]]
    cl = cluster_bands(data["results"], eps=args.eps, d_re=d_re, d_im=d_im)
    print(f"clustering: {cl.n_clusters} clusters @ eps={cl.eps:.4g}")
    stats = summarize_clusters(cl.points, cl.labels)
    tag = Path(args.data).stem

    for s in [x for x in stats if x["size"] >= args.min_size][:4]:
        cid = s["label"]
        print(f"\ncluster {cid} (size={s['size']}, "
              f"ReE[{s['E_re'][0]:+.2f},{s['E_re'][1]:+.2f}]):")
        verts, tri, vdata = refine_to_convergence(
            cl, cid, coeffs, degs, args.method, results_idx,
            edge_thresh=args.edge_thresh, match_tol=args.match_tol,
            max_iters=args.max_iters, n_procs=args.n_procs)
        outdir = OUT_ROOT / tag
        outdir.mkdir(parents=True, exist_ok=True)
        mesh_file = outdir / f"mesh_cluster{cid}.npz"
        np.savez(mesh_file, verts=verts, triangles=tri,
                 E=vdata["E"], mu1=vdata["mu1"], mu2=vdata["mu2"])
        print(f"  mesh saved -> {mesh_file}")
        if not args.no_figures:
            make_figure(tag, cid, verts, tri, vdata)


if __name__ == "__main__":
    main()
