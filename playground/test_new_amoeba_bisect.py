"""
Smoke test for the redesigned amoeba μ₂ bisection stage.

This script does NOT use collect_GBZ_subsets (the subset-assembly stage has
not been rebuilt yet).  It exercises the new bisection layer directly:

  - _try_fast_mu2          (θ₁=0 gap test)
  - detect_continuum       (std-based, pre-bisection)
  - _find_mu2_for_w2_zero  (continuum-first inner solve)
  - bisect_amoeba_ronkin_min (outer μ₁ bisection)

Run:
    PYTHONPATH=src python playground/test_new_amoeba_bisect.py
"""

import numpy as np
import matplotlib.pyplot as plt

from pygbz2d.core import CharPoly
from pygbz2d.amoeba.bisect import (
    bisect_amoeba_ronkin_min,
    _find_mu2_for_w2_zero,
    _try_fast_mu2,
)
from pygbz2d.amoeba.zm_extract import (
    AmoebaZeroManager,
    calculate_a2_average_winding,
    detect_continuum,
    find_crossings,
)


def build_HN2D_polynomial(J1, J2, gamma_1, gamma_2, delta_1=0.0, delta_2=0.0):
    """2D Hatano-Nelson chain-pair characteristic polynomial (basis='10')."""
    J11 = np.exp(gamma_1 + 1j * delta_1) * J1
    J12 = np.exp(-gamma_1 + 1j * delta_1) * np.conj(J1)
    J21 = np.exp(gamma_2 + 1j * delta_2) * J2
    J22 = np.exp(-gamma_2 + 1j * delta_2) * np.conj(J2)
    coeffs = np.array([1, -J11, -J12, -J21, -J22], dtype=complex)
    degs = np.array([
        [1, 0, 0],
        [0, -1, 0],
        [0, 1, 0],
        [0, 0, -1],
        [0, 0, 1],
    ], dtype=int)
    return coeffs, degs


def main():
    coeffs, degs = build_HN2D_polynomial(1.0, 1.0, 0.2, 0.3)
    poly = CharPoly(coeffs, degs)
    E = 1.0 + 0j

    print("=== 1. _try_fast_mu2 / detect_continuum at fixed μ₁ ===")
    for mu1 in (-0.5, 0.0, 0.2, 0.5):
        zm = AmoebaZeroManager(poly, E, mu1)
        zm.run()
        fast = _try_fast_mu2(zm, poly)
        groups = detect_continuum(zm, 1e-6)
        print(f"mu1={mu1:+.2f}: fast_ok={fast['ok']}, A={fast['A']:.12g}, "
              f"B={fast['B']:.12g}, continuum_groups={groups}")

    print("\n=== 2. inner μ₂ solve (_find_mu2_for_w2_zero) ===")
    for mu1 in (-0.5, 0.0, 0.2, 0.5):
        zm = AmoebaZeroManager(poly, E, mu1)
        zm.run()
        inner = _find_mu2_for_w2_zero(poly, E, mu1, -1.0, 1.0, _zm=zm)
        nz = len(inner.get("zeros") or [])
        print(f"mu1={mu1:+.2f}: mu2={inner['mu2']:.12g}, "
              f"is_continuum={inner['is_continuum']}, n_zeros={nz}, "
              f"winding={inner.get('winding')}")

    print("\n=== 3. calculate_a2_average_winding / find_crossings at μ₁=0.0, μ₂=0.3 ===")
    zm = AmoebaZeroManager(poly, E, 0.0)
    zm.run()
    w = calculate_a2_average_winding(zm, 0.0, 0.3)
    crossings = find_crossings(zm, 0.0, 0.3, return_refined=True)
    print(f"w2={w:.12g}, n_crossings={len(crossings)}")
    for b1, b2 in crossings:
        print(f"  beta1={b1:.6g}, beta2={b2:.6g}")

    print("\n=== 4. outer μ₁ bisection ===")
    for E_ref in (0.0 + 0j, 1.0 + 0j, 5.0 + 0j):
        res = bisect_amoeba_ronkin_min(
            poly, E_ref, mu1_low=-1.0, mu1_high=1.0,
            mu2_low=-1.0, mu2_high=1.0,
        )
        nz = len(res.get("zeros") or [])
        print(f"E={E_ref}: mu1={res['mu1']:.12g}, mu2={res['mu2']:.12g}, "
              f"is_continuum={res['is_continuum']}, n_zeros={nz}, "
              f"exit={res['_exit_reason']}")

    print("\n=== 5. collect_GBZ_subsets (subset assembly) ===")
    import pygbz2d.amoeba as amoeba
    from pygbz2d.core import PointSubset, LineSubset
    for E_ref in (0.0 + 0j, 1.0 + 0j, 5.0 + 0j):
        gbz = amoeba.collect_GBZ_subsets(coeffs, degs, E_ref, 0.0)
        print(f"E={E_ref}: success={gbz.success}, is_gbz={gbz.is_gbz}, "
              f"is_empty={gbz.is_empty}, index={gbz.index}, error={gbz.error}")
        plt.figure()
        plt.title(f"E={E_ref}")
        for s in gbz.subsets:
            if isinstance(s, PointSubset):
                plt.plot(np.angle(s.beta1), np.angle(s.beta2), '.')
                print(f"  PointSubset beta1={s.beta1:.6g}, beta2={s.beta2:.6g}")
            elif isinstance(s, LineSubset):
                plt.plot(np.angle(np.exp(1j * s.theta1_arr)), np.angle(s.beta2_arr))
                print(f"  LineSubset mu1={s.mu1:.12g}, n={len(s.theta1_arr)}")
    plt.show()


if __name__ == "__main__":
    main()
