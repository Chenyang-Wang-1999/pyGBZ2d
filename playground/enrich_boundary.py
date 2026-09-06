'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-27
Copyright © Department of Physics, Tsinghua University. All rights reserved

Pipeline steps 2-3: locate index-transition boundaries on the E sweep grid
by bisection, keeping EVERY probe's data points.

The sweep pickle is a mesh over the complex E plane (E_real x E_imag
grid); every node carries a GBZResult whose ``index`` (n_0D, n_1D) is a
phase label ((0, 0) = outside spectrum).  Any 4-neighbour grid edge whose
endpoints carry different labels straddles a transition — spectral
boundary or internal structure boundary (e.g. (6,0) <-> (12,0)).

Each transition edge is bisected along its segment: probe = one
collect_GBZ_subsets call; the bracket keeps two points with DIFFERENT
labels (a third label appearing mid-segment simply tightens the bracket
toward one of the transitions it hides).  ALL probe results are kept —
they are points arbitrarily close to the boundaries, exactly what the
torus mesh needs to close the band edges.

Long-run safety: a checkpoint (probes + edges done) is written every
CHECK_EVERY edges; --resume continues from it.  Output:
data/<stem>-enriched.pkl with the original results + all probes.

Usage:
    python playground/enrich_boundary.py                        # full run
    python playground/enrich_boundary.py --max-edges 3          # smoke test
    python playground/enrich_boundary.py --resume
'''

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import argparse
import pickle

import numpy as np

DEFAULT_DATA = "data/Haldane-gain-loss-amoeba.pkl"
CHECK_EVERY = 25
BISECT_ITERS = 6


def build_grid(data):
    """Reconstruct the (ni, nj) E grid + results grid from the sweep pickle.

    Row index = Im E, column index = Re E (meshgrid 'xy' + flatten order
    of the playground sweeps); verified against ``GBZResult.E_ref``.
    """
    E_real, E_imag = data["E_real"], data["E_imag"]
    results = data["results"]
    E_grid = (np.meshgrid(E_real, E_imag)[0]
              + 1j * np.meshgrid(E_real, E_imag)[1])
    assert len(results) == E_grid.size, \
        f"sweep size mismatch: {len(results)} vs {E_grid.size}"
    res_grid = np.empty(E_grid.shape, dtype=object)
    res_grid.flat[:] = results
    for k in (0, len(results) // 2, len(results) - 1):   # order check
        assert abs(res_grid.flat[k].E_ref - E_grid.flat[k]) < 1e-12, \
            "result/E grid ordering mismatch"
    return E_grid, res_grid


def node_label(res):
    """Phase label of a node: its index tuple, or None if the solve failed
    (edges touching failed nodes are skipped and counted)."""
    return tuple(res.index) if res.success else None


def transition_edges(res_grid):
    """4-neighbour edges whose endpoint labels differ (both known)."""
    ni, nj = res_grid.shape
    labels = np.empty(res_grid.shape, dtype=object)
    for i in range(ni):
        for j in range(nj):
            labels[i, j] = node_label(res_grid[i, j])
    edges = []
    n_failed_nodes = sum(l is None for l in labels.flat)
    for i in range(ni):
        for j in range(nj):
            a = labels[i, j]
            if a is None:
                continue
            for di, dj in ((0, 1), (1, 0)):
                i2, j2 = i + di, j + dj
                if i2 >= ni or j2 >= nj:
                    continue
                b = labels[i2, j2]
                if b is None or b != a:
                    edges.append(((i, j), (i2, j2)))
    return edges, n_failed_nodes


def _get_solver(method: str):
    import importlib
    return importlib.import_module(
        "pygbz2d." + method).collect_GBZ_subsets


def bisect_edge(Ea, idx_a, Eb, idx_b, coeffs, degs, method, iters):
    """Bisect the segment [Ea, Eb] between two different phase labels.

    Returns (probe_results, aborted).  The bracket always holds two
    points with different labels; a probe returning a THIRD label
    replaces the low side (any consistent choice converges to a
    transition hidden on the segment)."""
    solver = _get_solver(method)
    probes = []
    lo, hi, idx_lo = Ea, Eb, idx_a
    for _ in range(iters):
        Em = 0.5 * (lo + hi)
        res = solver(coeffs, degs, complex(Em))
        probes.append(res)
        if not res.success:
            return probes, True
        idx_m = tuple(res.index)
        if idx_m == idx_lo:
            lo = Em
        elif idx_m == tuple(idx_b):
            hi = Em
        else:
            lo, idx_lo = Em, idx_m
    return probes, False


def _enrich_one(task):
    """Pool worker: one transition edge -> (probes, aborted)."""
    Ea, idx_a, Eb, idx_b, coeffs, degs, method, iters = task
    return bisect_edge(Ea, idx_a, Eb, idx_b, coeffs, degs, method, iters)


def main(argv=None):
    p = argparse.ArgumentParser(description="E-grid boundary enrichment")
    p.add_argument("data", nargs="?", default=DEFAULT_DATA)
    p.add_argument("--method", choices=("amoeba", "sgbz"), default="amoeba",
                   help="must match the method that produced the sweep")
    p.add_argument("--iters", type=int, default=BISECT_ITERS)
    p.add_argument("--max-edges", type=int, default=None,
                   help="debug cap; default: all transition edges")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--n-procs", type=int, default=1,
                   help="worker processes for the probe solves (mp.Pool)")
    args = p.parse_args(argv)

    with open(args.data, "rb") as fp:
        data = pickle.load(fp)
    coeffs, degs = data["coeffs"], data["degs"]

    E_grid, res_grid = build_grid(data)
    edges, n_failed = transition_edges(res_grid)
    n_total = len(edges) if args.max_edges is None \
        else min(len(edges), args.max_edges)
    print(f"data: {args.data}  grid {res_grid.shape}, "
          f"{n_failed} failed nodes")
    print(f"transition edges: {len(edges)} total, running {n_total}, "
          f"{args.iters} bisection probes each")

    stem = Path(args.data).stem
    ckpt_file = Path(args.data).parent / f"{stem}-enrich-checkpoint.pkl"
    out_file = Path(args.data).parent / f"{stem}-enriched.pkl"

    probes, n_done, n_aborted = [], 0, 0
    if args.resume and ckpt_file.exists():
        with open(ckpt_file, "rb") as fp:
            ck = pickle.load(fp)
        probes, n_done, n_aborted = ck["probes"], ck["n_done"], ck["n_aborted"]
        print(f"resumed: {n_done} edges already done, "
              f"{len(probes)} probes collected")

    t0 = time.perf_counter()
    tasks = [(E_grid[i1, j1], node_label(res_grid[i1, j1]),
              E_grid[i2, j2], node_label(res_grid[i2, j2]),
              coeffs, degs, args.method, args.iters)
             for (i1, j1), (i2, j2) in edges[n_done:n_total]]

    def _consume(k, new, aborted):
        nonlocal n_done, n_aborted
        probes.extend(new)
        n_aborted += int(aborted)
        n_done += 1
        if n_done % 10 == 0 or n_done == n_total:
            dt = time.perf_counter() - t0
            rate_den = max(1, n_done - (n_total - len(tasks)))
            eta = (dt / rate_den) * (n_total - n_done)
            print(f"  edges {n_done}/{n_total}  probes={len(probes)}  "
                  f"aborted={n_aborted}  elapsed={dt/60:.1f} min  "
                  f"eta={eta/60:.1f} min", flush=True)
        if n_done % CHECK_EVERY == 0:
            with open(ckpt_file, "wb") as fp:
                pickle.dump({"probes": probes, "n_done": n_done,
                             "n_aborted": n_aborted}, fp)

    if args.n_procs and args.n_procs > 1:
        import multiprocessing as mp
        with mp.Pool(args.n_procs) as pool:
            for new, aborted in pool.imap(_enrich_one, tasks, chunksize=4):
                _consume(0, new, aborted)
    else:
        for t in tasks:
            new, aborted = _enrich_one(t)
            _consume(0, new, aborted)

    enriched = dict(data)
    enriched["results"] = list(data["results"]) + probes
    enriched["enrich_info"] = {
        "source": str(args.data), "method": args.method,
        "n_transition_edges": len(edges), "n_edges_done": n_done,
        "bisect_iters": args.iters, "n_probes": len(probes),
        "n_aborted": n_aborted,
    }
    with open(out_file, "wb") as fp:
        pickle.dump(enriched, fp)
    print(f"saved -> {out_file}  "
          f"(+{len(probes)} probe results, total "
          f"{len(enriched['results'])})")


if __name__ == "__main__":
    main()
