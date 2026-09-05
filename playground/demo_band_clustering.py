'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-27
Copyright © Department of Physics, Tsinghua University. All rights reserved

Demo CLI for pygbz2d.experimental.band_clustering (EXPERIMENTAL code).

Loads a pickled E-sweep ({'E_real', 'E_imag', 'results', ...} as produced
by the playground sweep scripts), clusters the GBZ point cloud into bands
via radius-graph connected components, prints the diagnostics (eps
stability scan, per-cluster report, inter-cluster margins) and saves
figures under playground/band_cluster_out/<tag>/.

Usage:
    python playground/demo_band_clustering.py                        # default data
    python playground/demo_band_clustering.py data/Haldane-gain-loss-y-SGBZ.pkl
    python playground/demo_band_clustering.py --eps 0.35 --alpha-E 0.25
'''

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import argparse

import numpy as np

from pygbz2d.experimental import cluster_bands, summarize_clusters

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_DATA = "data/Haldane-gain-loss-amoeba.pkl"
OUT_ROOT = Path(__file__).resolve().parent / "band_cluster_out"


def make_figures(cl, outdir: Path):
    bp = cl.points
    outdir.mkdir(parents=True, exist_ok=True)
    cmap = plt.get_cmap("tab10" if cl.n_clusters <= 10 else "tab20")
    colors = cmap(cl.labels % cmap.N)

    # -- eps stability --
    if cl.scan:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot([r["eps"] for r in cl.scan], [r["n_clusters"] for r in cl.scan],
                "o-")
        ax.axvline(cl.eps, color="r", ls="--", label=f"chosen eps={cl.eps:.3g}")
        ax.set_xscale("log")
        ax.set_xlabel("eps")
        ax.set_ylabel("n clusters")
        ax.set_title("cluster count vs eps (plateau = valid window)")
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / "eps_stability.png", dpi=130)
        plt.close(fig)

    # -- cluster views (cap plotted points for speed) --
    n = len(cl.labels)
    sel = np.arange(0, n, max(1, n // 40000))
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    from pygbz2d.core import TWO_PI
    views = [
        (bp.theta1[sel], bp.theta2[sel],
         r"$\theta_1$", r"$\theta_2$", "torus view"),
        (bp.mu1[sel], bp.mu2[sel],
         r"$\mu_1$", r"$\mu_2$", "log-modulus view"),
        (bp.E[sel].real, bp.E[sel].imag,
         r"Re $E$", r"Im $E$", "energy footprint"),
    ]
    for ax, (x, y, xl, yl, title) in zip(axes, views):
        ax.scatter(x, y, s=2, c=colors[sel], alpha=0.5, linewidths=0)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_title(f"{title}  ({cl.n_clusters} clusters)")
        if "torus" in title:
            ax.set_xlim(0, TWO_PI)
            ax.set_ylim(0, TWO_PI)
    fig.tight_layout()
    fig.savefig(outdir / "clusters.png", dpi=130)
    plt.close(fig)
    print(f"  figures -> {outdir}")


def main(argv=None):
    p = argparse.ArgumentParser(
        description="band clustering demo (pygbz2d.experimental)")
    p.add_argument("data", nargs="?", default=DEFAULT_DATA,
                   help="pickle with {'E_real','E_imag','results',...}")
    p.add_argument("--eps", type=float, default=None,
                   help="radius; default = middle of the widest plateau")
    p.add_argument("--alpha-E", type=float, default=None,
                   help="E-block downscale (default: experimental.ALPHA_E)")
    p.add_argument("--w-mu", type=float, default=None,
                   help="mu-block scale (default: auto 1/max(span,1))")
    p.add_argument("--max-scan-eps", type=float, default=None,
                   help="eps ladder cap (default: experimental.EPS_SCAN_MAX)")
    p.add_argument("--subsample", type=int, default=1,
                   help="stride over the results list (quick tests)")
    p.add_argument("--no-figures", action="store_true")
    args = p.parse_args(argv)

    import pickle
    with open(args.data, "rb") as fp:
        data = pickle.load(fp)
    results = data["results"][::args.subsample]
    print(f"data: {args.data}  ({len(results)} results after subsample "
          f"x{args.subsample})")

    cl = cluster_bands(results, eps=args.eps, alpha_E=args.alpha_E,
                       w_mu=args.w_mu, max_scan_eps=args.max_scan_eps)
    print(f"points: {len(cl.points.E)} "
          f"({int(cl.points.is_line.sum())} from LineSubsets, "
          f"{cl.points.n_dropped} dropped as non-finite)")

    knn = cl.knn_percentiles
    print("kNN(k=2) percentiles: "
          + " ".join(f"{k}={v:.4g}" for k, v in knn.items()))

    if cl.scan:
        print("\neps scan (eps : n_clusters / largest / singletons):")
        for r in cl.scan:
            print(f"  {r['eps']:8.4g} : {r['n_clusters']:6d} / "
                  f"{r['largest']:7d} / {r['n_singletons']:5d}")
    if cl.eps_window is not None:
        lo, hi = cl.eps_window
        print(f"\nwidest plateau [{lo:.4g}, {hi:.4g}] "
              f"-> eps={cl.eps:.4g}")
    print(f"\nfinal clustering @ eps={cl.eps:.4g}: "
          f"{cl.n_clusters} clusters")

    stats = summarize_clusters(cl.points, cl.labels)
    big = [s for s in stats if s["size"] >= 5]
    small = [s for s in stats if s["size"] < 5]
    for s in big[:12]:
        print(f"  [{s['label']:>3}] size={s['size']:>7}  "
              f"slices={s['n_slices']:>5}  pts/slice~{s['pts_per_slice_med']:.1f}  "
              f"line={s['frac_line']:.2f}  "
              f"ReE[{s['E_re'][0]:+.2f},{s['E_re'][1]:+.2f}] "
              f"ImE[{s['E_im'][0]:+.2f},{s['E_im'][1]:+.2f}] "
              f"mu1[{s['mu1'][0]:+.2f},{s['mu1'][1]:+.2f}] "
              f"mu2[{s['mu2'][0]:+.2f},{s['mu2'][1]:+.2f}]")
    if small:
        print(f"  (+ {len(small)} small clusters, "
              f"{sum(s['size'] for s in small)} points total)")

    if cl.margins:
        print("inter-cluster min distances (cluster_a, cluster_b, dist):")
        for a, b, d in cl.margins[:8]:
            print(f"  ({a:>3}, {b:>3}) -> {d:.4g}  "
                  f"(safety x{d / cl.eps:.1f} over eps)")

    if not args.no_figures:
        make_figures(cl, OUT_ROOT / Path(args.data).stem)


if __name__ == "__main__":
    main()
