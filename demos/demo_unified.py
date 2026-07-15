'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-14
Copyright © Department of Physics, Tsinghua University. All rights reserved

Demonstration of the unified GBZ API using the 2D HN model.

Shows:
  - Building a characteristic polynomial
  - Running check_SGBZ (returns GBZResult directly)
  - Running collect_GBZ_subsets (returns GBZResult directly)
  - Iterating over subsets with match/case
  - Lazy beta2_arr fill for LineSubset
  - Converting to (E, k1, k2) triplets
'''

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from cmath import exp
import poly_tools as pt

import brute_force_SGBZ as bfs
import brute_force_amoeba as bfa
from gbz_types import PointSubset, LineSubset, GBZResult


def build_HN2D_polynomial(J1, J2, gamma_1, gamma_2, delta_1, delta_2):
    """Construct 2D HN model characteristic polynomial."""
    J11 = exp(gamma_1 + 1j * delta_1) * J1
    J12 = exp(-gamma_1 + 1j * delta_1) * np.conj(J1)
    J21 = exp(gamma_2 + 1j * delta_2) * J2
    J22 = exp(-gamma_2 + 1j * delta_2) * np.conj(J2)

    coeffs = np.array([1, -J11, -J12, -J21, -J22], dtype=complex)
    degs = np.array([
        [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1],
    ], dtype=int)
    return coeffs, degs


def main():
    # Parameters
    J1, J2 = 1.0, 1.0
    gamma_1, gamma_2 = 0.2, 0.3
    delta_1, delta_2 = 0.0, 0.0

    coeffs, degs = build_HN2D_polynomial(J1, J2, gamma_1, gamma_2, delta_1, delta_2)

    # Analytic expectations:
    # [10]-SGBZ: mu1 = gamma_1 = 0.2, |beta1| = exp(0.2), |beta2| = exp(0.3)
    # Amoeba GBZ = [10]-SGBZ
    # Spectrum: E = 2cos(theta1) + 2cos(theta2) in [-4, 4]

    E_test = 1.0 + 0j  # Inside the spectrum

    print("=" * 60)
    print(f"2D HN Model: J1={J1}, J2={J2}, gamma=({gamma_1},{gamma_2})")
    print(f"Test energy: E = {E_test}")
    print(f"Analytic SGBZ: mu1 = {gamma_1}")
    print(f"Analytic GBZ: |beta1| = {exp(gamma_1):.6f}, |beta2| = {exp(gamma_2):.6f}")
    print("=" * 60)

    # ---- SGBZ ----
    print("\n--- check_SGBZ ([10]-SGBZ) ---")
    gbz_sgbz = bfs.check_SGBZ(coeffs, degs, E_test, 0.0, N_points=201)
    print(f"  is_gbz: {gbz_sgbz.is_gbz}")
    print(f"  index: {gbz_sgbz.index} (n_0D={gbz_sgbz.index[0]}, n_1D={gbz_sgbz.index[1]})")
    print(f"  n_subsets: {len(gbz_sgbz.subsets)}")

    for i, subset in enumerate(gbz_sgbz.subsets):
        if isinstance(subset, PointSubset):
            print(f"\n  [{i}] PointSubset:")
            print(f"      E     = {subset.E}")
            print(f"      beta1 = {subset.beta1:.6f}  (|beta1| = {abs(subset.beta1):.6f}, mu1 = {subset.mu1:.6f})")
            print(f"      beta2 = {subset.beta2:.6f}  (|beta2| = {abs(subset.beta2):.6f})")
            E_t, k1, k2 = subset.as_triplet()
            print(f"      triplet: (E={E_t:.4f}, k1={k1:.4f}, k2={k2:.4f})")
        elif isinstance(subset, LineSubset):
            print(f"\n  [{i}] LineSubset:")
            print(f"      mu1 = {subset.mu1:.6f}, theta1 ∈ [{subset.theta1_start:.4f}, {subset.theta1_end:.4f}]")
            print(f"      width = {subset.theta1_width:.4f} rad")
            if not subset.is_loaded():
                print("      beta2_arr: not loaded (lazy)")
                subset.fill_beta2()
                print(f"      beta2_arr loaded: shape = {subset.beta2_arr.shape}")

    # ---- Amoeba ----
    print("\n--- collect_GBZ_subsets ---")
    gbz_amoeba = bfa.collect_GBZ_subsets(coeffs, degs, E_test, 0.0)
    print(f"  is_gbz: {gbz_amoeba.is_gbz}")
    print(f"  index: {gbz_amoeba.index}")
    print(f"  n_subsets: {len(gbz_amoeba.subsets)}")

    for i, subset in enumerate(gbz_amoeba.subsets):
        if isinstance(subset, PointSubset):
            print(f"\n  [{i}] PointSubset: |beta1|={abs(subset.beta1):.6f}, |beta2|={abs(subset.beta2):.6f}")
        elif isinstance(subset, LineSubset):
            print(f"\n  [{i}] LineSubset: mu1={subset.mu1:.6f}")

    # ---- Triplet conversion ----
    print("\n--- Triplet conversion ---")
    triplets = bfs.convert_gbz_to_triplets(gbz_sgbz)
    print(f"  n_triplets: {len(triplets)}")
    for i, (E, k1, k2) in enumerate(triplets[:5]):
        print(f"  [{i}] E={E.real:.4f}{E.imag:+.4f}j, "
              f"k1={k1.real:.4f}{k1.imag:+.4f}j, "
              f"k2={k2.real:.4f}{k2.imag:+.4f}j")
    if len(triplets) > 5:
        print(f"  ... ({len(triplets) - 5} more)")

    # ---- Outside spectrum ----
    print("\n--- Outside spectrum (E=5.0) ---")
    gbz_out = bfs.check_SGBZ(coeffs, degs, 5.0 + 0j, 0.0, N_points=101)
    print(f"  is_gbz: {gbz_out.is_gbz}, index: {gbz_out.index}")


if __name__ == "__main__":
    main()
