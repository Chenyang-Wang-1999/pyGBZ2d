'''cProfile workload for the SGBZ pipeline (diagnostics only — no fixes).

Three representative energies, chosen to cover the distinct code paths of
``brute_force_SGBZ.collect_GBZ_subsets``:

  1. HN2D chain-pair, E inside spectrum   → continuum LineSubset path
     (Mu2MidZM.analyze + extract_continuum_linesubsets + MR joins)
  2. HN2D chain-pair, E outside spectrum  → bracket expansion only,
     early exit (empty result)
  3. Haldane gain/loss 2x2 supercell, E=1.212 grid point → heavier
     polynomial (64 terms, β₂ degree 8), continuum with 6 LineSubsets

Usage:
    python diagnostics/profile_sgbz.py [--out profile_sgbz.prof]

Writes the raw pstats dump and prints hotspots by cumulative time, by
tottime, and per-module aggregates.
'''
import sys
import cProfile
import pstats
import argparse
import importlib.util
from pathlib import Path
from cmath import exp
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import brute_force_SGBZ as bfs


def build_HN2D_polynomial(J1, J2, gamma_1, gamma_2, delta_1, delta_2):
    J11 = exp(gamma_1 + 1j * delta_1) * J1
    J12 = exp(-gamma_1 + 1j * delta_1) * np.conj(J1)
    J21 = exp(gamma_2 + 1j * delta_2) * J2
    J22 = exp(-gamma_2 + 1j * delta_2) * np.conj(J2)
    coeffs = np.array([1, -J11, -J12, -J21, -J22], dtype=complex)
    degs = np.array([
        [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1],
    ], dtype=int)
    return coeffs, degs


def build_haldane_supercell():
    spec = importlib.util.spec_from_file_location(
        "haldane_gainloss", ROOT / "demos" / "Haldane-model-gainloss.py")
    hg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(hg)
    model = hg.Haldane_non_Hermitian_phase(*hg.ALL_PARAMS)
    model = model.get_supercell(
        [(0, 0), (1, 0)],
        np.array([[1, 1], [-1, 1]], dtype=int),
    )
    return model.get_characteristic_polynomial_data()


def workload():
    # 1. HN2D in-spectrum (continuum path)
    coeffs, degs = build_HN2D_polynomial(1.0, 1.0, 0.2, 0.3, 0.0, 0.0)
    r = bfs.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, None, debug_mode=True)
    print(f"[workload 1] HN2D E=1.0      -> index={r.index}")

    # 2. HN2D outside spectrum (fast exit)
    r = bfs.collect_GBZ_subsets(coeffs, degs, 5.0 + 0j, None, debug_mode=True)
    print(f"[workload 2] HN2D E=5.0      -> index={r.index}")

    # 3. Haldane supercell continuum (heavy polynomial)
    hc, hd = build_haldane_supercell()
    grid = np.linspace(-3.1, 4.6, 51)
    E_ref = complex(grid[np.argmin(np.abs(grid - 1.212))])
    r = bfs.collect_GBZ_subsets(hc, hd, E_ref, None, debug_mode=True)
    print(f"[workload 3] Haldane E={E_ref.real:.4f} -> index={r.index}")


def per_module_aggregate(stats: pstats.Stats, top: int = 15) -> None:
    """Aggregate tottime by source file (project files kept separate)."""
    agg = defaultdict(lambda: [0.0, 0])  # file -> [tottime, ncalls]
    for (filename, _, name), (cc, nc, tt, ct, callers) in stats.stats.items():
        key = filename
        agg[key][0] += tt
        agg[key][1] += nc
    rows = sorted(agg.items(), key=lambda kv: -kv[1][0])[:top]
    total = sum(v[0] for v in agg.values())
    print(f"\n=== per-file tottime aggregate (total {total:.1f}s) ===")
    for filename, (tt, nc) in rows:
        rel = str(Path(filename).replace(ROOT, '.'))
        print(f"{tt:9.2f}s {100*tt/total:5.1f}%  {nc:>10,d} calls  {rel}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "diagnostics" / "profile_sgbz.prof"))
    args = ap.parse_args()

    prof = cProfile.Profile()
    prof.enable()
    workload()
    prof.disable()

    prof.dump_stats(args.out)
    stats = pstats.Stats(prof)

    print("\n=== top 35 by cumulative time ===")
    stats.sort_stats("cumulative").print_stats(35)
    print("\n=== top 35 by tottime (self) ===")
    stats.sort_stats("tottime").print_stats(35)
    per_module_aggregate(stats)
    print(f"\nraw pstats dump: {args.out}")


if __name__ == "__main__":
    main()
