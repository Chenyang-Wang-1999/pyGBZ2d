"""
Minimal reproduction: Hermitian Haldane amoeba returns PointSubsets that
are NOT zeros of the characteristic polynomial.

Setup
-----
Hermitian Haldane model (gamma = 0, real M):

    Haldane_non_Hermitian_phase(t1=1, t2=0.5, phi=pi/3, M=0.5, gamma=0)

At E_ref = -1.5, ``pygbz2d.amoeba.collect_GBZ_subsets`` returns
``index == (49, 2)``: 49 PointSubsets + 2 LineSubsets.  Substituting the
returned PointSubsets back into the SAME CharPoly gives, e.g.,

    beta1 = 0.5000000001409388 + 0.8660254037030676j
    beta2 = 7520.721655046422  - 4342.087421209515j
    f(E, beta1, beta2) = -0.9946858994321355 + 2.7e-9j   (|f| ~ 0.995)

while other returned PointSubsets give |f| ~ 1e-15.  A PointSubset is
supposed to be an isolated GBZ point, i.e. a solution of f = 0.

Exit status
-----------
1  -> bug reproduced (at least one PointSubset is not a zero)
0  -> bug NOT reproduced
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pygbz2d.amoeba as bfa
import pygbz2d.sgbz as bfs
from pygbz2d.core import CharPoly, PointSubset, TWO_PI
from pygbz2d.amoeba import AmoebaZeroManager
import pygbz2d.continuation as continuation
import matplotlib.pyplot as plt

# Haldane-model-gainloss.py has hyphens in its name -> importlib.
SPEC = importlib.util.spec_from_file_location(
    "haldane_model",
    Path(__file__).parent / "Haldane-model-gainloss.py",
)
HALDANE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HALDANE)

E_REF = -1.5
HERMITIAN_PARAMS = (1.0, 0.5, np.pi / 3, 0.5, 0.0)  # t1, t2, phi, M, gamma
BAD_TOL = 1e-6  # anything above this is definitely not a numerical zero


def debug_zero_manager():
    model = HALDANE.Haldane_non_Hermitian_phase(*HERMITIAN_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()
    poly = CharPoly(coeffs, degs)

    res = bfs.collect_GBZ_subsets(coeffs, degs, E_REF,
                                     debug_mode=True)
    print(res)
    manager = continuation.ZeroManager(poly, -1.5, 0)
    manager.run()

    print(len(manager.segments))
    for seg in manager.segments:
        plt.plot(seg.theta1_arr, np.log(np.abs(seg.tracked_roots)), '-')

    for s in res.subsets:
        if isinstance(s, PointSubset):
            plt.plot(np.angle(s.beta1) % TWO_PI, np.log(np.abs(s.beta2)), 'x')
        else:
            plt.plot(s.theta1_arr, np.log(np.abs(s.beta2_arr)), ':')

    plt.figure()
    for s in res.subsets:
        if isinstance(s, PointSubset):
            plt.plot(np.angle(s.beta1) % TWO_PI, np.angle(s.beta2), 'x')
        else:
            plt.plot(s.theta1_arr, np.angle(s.beta2_arr), '-')
    plt.show()


def main() -> int:
    model = HALDANE.Haldane_non_Hermitian_phase(*HERMITIAN_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()
    coeffs = np.asarray(coeffs, dtype=complex)

    charpoly = CharPoly(coeffs, degs)
    result = bfa.collect_GBZ_subsets(coeffs, degs, E_REF, 0.0,
                                     debug_mode=True)

    points = [s for s in result.subsets if isinstance(s, PointSubset)]
    n_lines = sum(1 for s in result.subsets if not isinstance(s, PointSubset))

    print(f"E_ref          = {E_REF}")
    print(f"result.index   = {result.index}")
    print(f"n_PointSubset  = {len(points)}")
    print(f"n_LineSubset   = {n_lines}")

    residuals = []
    for p in points:
        f = charpoly.eval_val((p.E, p.beta1, p.beta2))
        residuals.append(abs(f))

    residuals = np.array(residuals)
    print(f"|f| at returned PointSubsets: "
          f"min={residuals.min():.3e} "
          f"median={np.median(residuals):.3e} "
          f"max={residuals.max():.3e}")

    bad = np.flatnonzero(residuals > BAD_TOL)
    if len(bad):
        print(f"\n{len(bad)} PointSubset(s) with |f| > {BAD_TOL:.0e}:")
        for k in bad[:5]:
            p = points[k]
            f = charpoly.eval_val((p.E, p.beta1, p.beta2))
            print(f"  [{k}] E={p.E!r}")
            print(f"      beta1 = {p.beta1!r}")
            print(f"      beta2 = {p.beta2!r}")
            print(f"      f     = {f!r}   |f| = {abs(f):.6e}")
        print("\nBUG REPRODUCED: amoeba returned PointSubset(s) that are "
              "not zeros of the characteristic polynomial.")
        return 1

    print("\nBUG NOT REPRODUCED: every returned PointSubset satisfies f = 0.")
    return 0


if __name__ == "__main__":
    # raise SystemExit(main())
    debug_zero_manager()
