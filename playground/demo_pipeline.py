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
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse

import numpy as np

from pygbz2d.experimental import cluster_bands, summarize_clusters
from demo_torus_mesh import (build_cluster_mesh, dedup_vertices,
                             periodic_delaunay, mesh_topology,
                             describe_topology, seam_display_filter)
from demo_mesh_refine import (build_E_interpolator, dyadic_geometry_refine,
                              solve_batch, edge_lengths,
                              edge_len_percentiles, triangle_areas,
                              MATCH_TOL)

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
    verts, tri, vdata, _ = build_cluster_mesh(cl, cid)
    # per-real-vertex phase label (n_0D, n_1D) via the source result
    idx = np.array([results_idx[s] for s in vdata["slice_idx"]], dtype=int)
    print(f"  base mesh: {describe_topology(mesh_topology(tri))}")
    print(f"  base edges: {np.round(edge_len_percentiles(verts, tri), 4)} "
          f"| threshold={edge_thresh}")

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
        print(f"  [round {it}] inner geometry: {len(pend_pos)} midpoints "
              f"over {n_levels} nesting level(s) [{t_geom:.0f}s]")
        if len(pend_pos) == 0:
            print("  stalled: no midpoints available")
            break

        new_theta, new_E, new_mu, new_idx, st = solve_batch(
            coeffs, degs, pend_pos, pend_E, aL, aR,
            vdata["E"], idx, method, match_tol, n_procs=n_procs)
        print(f"    solved: direct={st['n_direct']} "
              f"boundary={st['n_boundary']} rejected={st['n_rejected']} "
              f"(fallback probes={st['n_probes']}, "
              f"match p50={st['match_p50']:.4g}, "
              f"boundary p50={st['boundary_p50']:.4g}) "
              f"[{time.perf_counter() - t0 - t_geom:.0f}s]")

        if len(new_theta) == 0:
            print("  stalled: no accepted points this round")
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
        print(f"    mesh: V={len(verts)} F={len(tri)} "
              f"max edge -> {L.max():.4f}")
    else:
        print(f"  WARNING: hit max_iters={max_iters} before convergence "
              f"(max edge {edge_lengths(verts, tri)[1].max():.4f})")

    if converged_at is not None:
        print(f"  converged after {converged_at} solving round(s) "
              f"(max edge {edge_lengths(verts, tri)[1].max():.4f} "
              f"<= {edge_thresh})")
    print(f"  final: {describe_topology(mesh_topology(tri))}")
    print(f"  edges p50/p90/p99/max: "
          f"{np.round(edge_len_percentiles(verts, tri), 4)}")
    a = triangle_areas(verts, tri)
    print(f"  areas p50={np.median(a):.3e} max={a.max():.3e}")
    return verts, tri, vdata


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
