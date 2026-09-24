"""Pseudo-arclength vs. uniform-theta1 sampling for beta2 root tracking.

Both methods use exactly the same tooling and the same budget:

  * the same polynomial (the ``playground/demo_zero_manager.py`` case)
  * a FIXED step (no adaptive step-size control anywhere)
  * Hungarian matching (``core.hungarian_match_indices``, chordal distance)
    between consecutive *solved* root sets -- no tangent-based matching

The only difference is how the next theta1 is chosen:

  * arclength : theta1 advances by h / ||V|| with V = [1, d ln beta2_j/dtheta1],
                i.e. a CONSTANT step h in (theta1, ln beta2)-arclength
  * uniform   : theta1 advances by a constant 2*pi/N, with N set to the number
                of solves the arclength run needed for the same period

Because ||V|| -> infinity as two branches approach, the arclength rule spends
its samples where the branches move fast, so the neighbouring root sets stay
close and the Hungarian match stays unambiguous.  The uniform rule keeps its
samples evenly spread in theta1 and jumps straight over the close approach.
The uniform nodes are given a generic phase (GRID_PHASE) so that they do not
accidentally land on theta1 = pi, where the two branches meet.

Model (c0 = 2 cosh(mu1) + DELTA):

    f = beta2**2 - beta1 - beta1**-1 - c0,      beta1 = exp(mu1 + i*theta1)

At DELTA = 0 the two beta2 branches meet in a double root at theta1 = pi.
A small perturbation splits it: the branches then pass within 2*sqrt(|DELTA|)
of each other, which is the regime where a too-coarse theta1 grid makes the
Hungarian match pick the wrong partner.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
PAPER_PATH = str(Path(__file__).resolve().parents[1])
CURR_PATH = str(Path(__file__).resolve().parent)
plt.style.use(PAPER_PATH + "/paper_plot.mplstyle")
CM = 1 / 2.54

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from pygbz2d.core import TWO_PI, CharPoly, hungarian_match_indices  # noqa: E402

MU1 = 0.2
#: Perturbation of c0 = 2 cosh(mu1); splits the double root at theta1 = pi into
#: +-i*sqrt(|DELTA|).  Negative sign puts the split pair across the imaginary
#: axis, which is the branch geometry the uniform grid fails on.
DELTA = -1e-4
#: The FIXED arclength step.  The uniform grid is then handed the same number
#: of root solves, so the only difference between the two runs is the rule that
#: picks the next theta1.
ARC_STEP = 0.2
#: Generic grid phase, so the uniform nodes do not happen to land on the close
#: approach (which is not known a priori in a real problem).
GRID_PHASE = 0.31
N_SUB = 300         # sub-steps per interval for the reference matching


def build_poly() -> CharPoly:
    coeffs = [1, -1, -1, -(2 * np.cosh(MU1) + DELTA)]
    degs = np.array([[0, 0, 2], [0, 1, 0], [0, -1, 0], [0, 0, 0]], dtype=int)
    return CharPoly(coeffs, degs)


def build_poly2() -> CharPoly:
    coeffs = [1, 2, -1, -1, -(2 * np.cosh(MU1) + DELTA) + 1]
    degs = np.array([[0, 0, 2], [0, 0, 1], [0, 1, 0], [0, -1, 0], [0, 0, 0]], dtype=int)
    return CharPoly(coeffs, degs)


def solve_roots(poly: CharPoly, theta1: float) -> np.ndarray:
    """The two beta2 zeros at (E=0, beta1=exp(mu1 + i*theta1))."""
    return poly.solve_roots_1d((0, 1), (0, np.exp(MU1 + 1j * theta1)), (2,))


def track(poly: CharPoly, thetas: np.ndarray):
    """Solve at each theta1 and Hungarian-match consecutive sets.

    ``perm[i]`` is the raw matching between the solved sets i and i+1, so it is
    directly comparable with the reference matching of the same interval.
    """
    raw = [solve_roots(poly, t) for t in thetas]
    order = np.arange(raw[0].size)
    tracked, perms = [raw[0]], []
    for i in range(thetas.size - 1):
        perm = hungarian_match_indices(raw[i], raw[i + 1])
        perms.append(perm)
        order = perm[order]
        tracked.append(raw[i + 1][order])
    return np.array(tracked), perms


def reference_perm(poly: CharPoly, theta0: float, theta1: float) -> np.ndarray:
    """The same Hungarian rule composed over N_SUB sub-steps of the interval.

    This is the operational ground truth: with sub-steps small enough the
    matching just follows the branches continuously.
    """
    raw = [solve_roots(poly, t)
           for t in np.linspace(theta0, theta1, N_SUB + 1)]
    order = np.arange(raw[0].size)
    for i in range(N_SUB):
        order = hungarian_match_indices(raw[i], raw[i + 1])[order]
    return order


def mismatch_intervals(poly: CharPoly, thetas: np.ndarray, perms) -> list:
    """Intervals where one Hungarian match disagrees with the dense matching."""
    return [i for i in range(thetas.size - 1)
            if not np.array_equal(perms[i],
                                  reference_perm(poly, thetas[i], thetas[i + 1]))]


def main() -> None:
    poly = build_poly()
    # Self-adaptive first: its sample count sets the budget for the uniform grid.

    from pygbz2d.continuation import ZeroManager, StepControl
    zm = ZeroManager(poly, 0, MU1)
    zm.run(h0=ARC_STEP, ctrl=StepControl(atol=0.1, rtol=0.1))
    theta_a = []
    roots_a = []
    for seg in zm.segments:
        theta_a.append(seg.theta1_arr[:-1])
        roots_a.append(seg.tracked_roots[:-1,:])
    theta_a.append([zm.segments[-1].theta1_arr[-1]])
    roots_a.append(zm.segments[-1].tracked_roots[-1, :])

    theta_a = np.concatenate(theta_a)
    roots_a = np.vstack(roots_a)

    theta_u = TWO_PI * (np.arange(theta_a.size) + GRID_PHASE) / theta_a.size
    roots_u, perms_u = track(poly, theta_u)


    fig = plt.figure(figsize=(3 * CM / 0.8, 3 * CM / 0.8))
    ax = fig.gca()
    ax.set_position([0.2, 0.2, 0.8, 0.8])
    ax.set_ylim([-0.2, 0.2])
    ax.plot(roots_u.real, roots_u.imag, '.-')
    fig.savefig(CURR_PATH + "/Figures/uniform.pdf", backend="cairo")

    fig = plt.figure(figsize=(3 * CM / 0.8, 3 * CM / 0.8))
    ax = fig.gca()
    ax.set_position([0.2, 0.2, 0.8, 0.8])
    ax.set_ylim([-0.2, 0.2])
    ax.plot(roots_a.real, roots_a.imag, '.-')
    fig.savefig(CURR_PATH + "/Figures/self-adaptive.pdf", backend="cairo")
 
    plt.show()

if __name__ == "__main__":
    main()
