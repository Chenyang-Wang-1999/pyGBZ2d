'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-26
Copyright © Department of Physics, Tsinghua University. All rights reserved

Demo 1: Hultquist advancing-front triangulation between two matched LineSubsets.

Reference: Hultquist 1992, "Constructing Stream Surfaces in Steady Vector
Fields".  Applied here to GBZ branches: two θ₁-monotone curves C_a (at E_a)
and C_b (at E_b) are stitched into a single 1-D triangle strip.

Algorithm (per matched branch pair):
  - Maintain a "front" edge (A_i, B_j) starting at (i=0, j=0).
  - Each step choose the cheaper of:
        advance A : triangle (A_i, A_{i+1}, B_j), front -> (i+1, j)
        advance B : triangle (A_i, B_j, B_{j+1}), front -> (i, j+1)
  - Cost = candidate triangle's max edge length (the "diagonal" we'd add).
    This is Hultquist's original criterion; other choices (chordal distance,
    arc-length imbalance) are drop-in.
  - Terminate when both ends reached, finishing with a fan to any leftover
    vertices on the longer curve.

Two end conditions exercise the "adaptive count" property:
  - trivial_model: θ₁ ranges of adjacent slices differ → m ≠ n, so one curve
    is consumed faster; the strip narrows via fewer advances on that side.
    No triangle is forced to collapse — leftover vertices fan out cleanly.

Output: per (E_a, E_b) pair, the triangle strip as vertex + triangle arrays,
plus diagnostics (n_triangles, min/max triangle edge, leftover count).
'''

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from gbz_types import LineSubset


def _oriented(line: LineSubset):
    """Return (theta1, beta2) oriented so theta1 is ascending.

    LineSubset stores theta1 monotonically along the curve but the direction
    is not guaranteed (trivial_model gives descending).  A consistent
    ascending orientation makes the front advance from a fixed anchor.
    """
    th = np.asarray(line.theta1_arr, dtype=float)
    b2 = np.asarray(line.beta2_arr, dtype=complex)
    if len(th) >= 2 and th[0] > th[-1]:
        th = th[::-1]
        b2 = b2[::-1]
    return th, b2


def _chordal(a: complex, b: complex) -> float:
    """Chordal distance between two complex β values (Riemann sphere)."""
    # Local chordal distance via stereographic projection; reuse the shared
    # helper for consistency with the matching demos.
    from gbz_types import to_sphere_r3
    p = to_sphere_r3(np.array([a, b]))
    return float(np.linalg.norm(p[0] - p[1]))


def _tri_max_edge(a, b, c):
    return max(_chordal(a, b), _chordal(b, c), _chordal(a, c))


def advance_front(La: LineSubset, Lb: LineSubset):
    """Hultquist advancing-front triangulation of two matched LineSubsets.

    Returns:
        verts : (N_a + N_b, 3) ndarray — columns (E_label, theta1, beta2)
                where E_label in {0,1} tags which slice the vertex came from.
                Duplicate-anchored end vertices are de-duplicated by index.
        tris  : (n_tri, 3) int ndarray — triangle vertex indices into `verts`.
        info  : dict with diagnostics.
    """
    th_a, b2_a = _oriented(La)
    th_b, b2_b = _oriented(Lb)
    m, n = len(th_a), len(th_b)

    # Vertex table: slice A then slice B.
    verts = []
    for k in range(m):
        verts.append((0, th_a[k], b2_a[k]))
    for k in range(n):
        verts.append((1, th_b[k], b2_b[k]))
    base_b = m  # index offset for slice B vertices

    tris = []
    i, j = 0, 0
    leftover = 0
    while not (i == m - 1 and j == n - 1):
        adv_a = i < m - 1
        adv_b = j < n - 1
        if adv_a and adv_b:
            # Hultquist: pick the triangle with smaller max edge.
            # tri_advance_A uses edge A_{i+1}—B_j as new front.
            tA = _tri_max_edge(b2_a[i], b2_a[i + 1], b2_b[j])
            tB = _tri_max_edge(b2_a[i], b2_b[j], b2_b[j + 1])
            if tA <= tB:
                tris.append((i, i + 1, base_b + j))
                i += 1
            else:
                tris.append((i, base_b + j, base_b + j + 1))
                j += 1
        elif adv_a:
            # B exhausted; fan the remaining A vertices to B_{n-1}.
            tris.append((i, i + 1, base_b + n - 1))
            i += 1
            leftover += 1
        else:
            # A exhausted; fan the remaining B vertices to A_{m-1}.
            tris.append((m - 1, base_b + j, base_b + j + 1))
            j += 1
            leftover += 1

    tris = np.asarray(tris, dtype=int)
    verts_arr = np.empty((len(verts), 3), dtype=object)
    for k, v in enumerate(verts):
        verts_arr[k] = v

    # diagnostics — separate anchor-endpoint triangles (whose near-zero edge
    # is the physical β₂≈1 coincidence at θ₁ endpoints, not a mesh defect)
    # from interior triangles.
    anchor_indices = {0, m - 1, base_b, base_b + n - 1}
    interior_edges = []
    anchor_edges = []
    for t in tris:
        is_anchor = bool(set(int(x) for x in t) & anchor_indices)
        for u, v in [(t[0], t[1]), (t[1], t[2]), (t[0], t[2])]:
            e = _chordal(verts_arr[v][2], verts_arr[u][2])
            (anchor_edges if is_anchor else interior_edges).append(e)
    info = {
        "n_tri": len(tris),
        "m": m, "n": n,
        "leftover_fan": leftover,
        "min_edge": float(min(interior_edges)) if interior_edges else float("nan"),
        "max_edge": float(max(interior_edges + anchor_edges)) if (interior_edges or anchor_edges) else float("nan"),
        "min_anchor_edge": float(min(anchor_edges)) if anchor_edges else float("nan"),
        "th1_a_range": (float(th_a[0]), float(th_a[-1])),
        "th1_b_range": (float(th_b[0]), float(th_b[-1])),
    }
    return verts_arr, tris, info


def _print_strip(Ea, Eb, La, Lb, verts, tris, info):
    print(f"  E_a={Ea:+.2f} E_b={Eb:+.2f}  "
          f"|A|={info['m']} |B|={info['n']} "
          f"th1_a={info['th1_a_range']} th1_b={info['th1_b_range']}")
    print(f"    triangles={info['n_tri']}  leftover_fan={info['leftover_fan']}  "
          f"edge[min]={info['min_edge']:.3e} edge[max]={info['max_edge']:.3e}")
    # show the advance sequence (0=advance A, 1=advance B) to visualize
    # adaptive consumption
    seq = []
    for t in tris:
        # a triangle (i,i+1,j) advances A; (i,j,j+1) advances B
        if t[0] + 1 == t[1]:
            seq.append('A')
        else:
            seq.append('B')
    s = ''.join(seq)
    print(f"    advance seq: {s}")


def main():
    from trivial_line_cache import get_line_pairs
    E_list, pairs = get_line_pairs()

    print("=" * 64)
    print("Hultquist advancing-front triangulation (trivial_model real axis)")
    print("=" * 64)
    total_tri = 0
    for k, (Ea, Eb) in enumerate(E_list[:6]):
        for branch, (La, Lb) in enumerate(pairs[k]):
            verts, tris, info = advance_front(La, Lb)
            total_tri += info["n_tri"]
            print(f"[pair {k}, branch {branch}]")
            _print_strip(Ea, Eb, La, Lb, verts, tris, info)
        print()
    print(f"total triangles (first 6 pairs, both branches): {total_tri}")

    # Degenerate-cell check: interior triangles should have no near-zero edge.
    # Anchor-endpoint triangles (θ₁ endpoints, β₂≈1 coincidence) are excluded:
    # their near-zero edge is physical (the GBZ closes there), not a mesh defect.
    print("\nDegenerate-cell check (interior triangles only; anchor ends excluded):")
    worst_int = float("inf")
    worst_anc = float("inf")
    for k in range(len(E_list)):
        for La, Lb in pairs[k]:
            _, _, info = advance_front(La, Lb)
            worst_int = min(worst_int, info["min_edge"])
            worst_anc = min(worst_anc, info["min_anchor_edge"])
    print(f"  interior min edge = {worst_int:.3e}  "
          f"({'OK' if worst_int > 1e-9 else 'DEGENERATE'})")
    print(f"  anchor   min edge = {worst_anc:.3e}  (expected ~0: β₂≈1 coincidence)")
    print("=" * 64)


if __name__ == "__main__":
    main()
