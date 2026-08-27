"""
MR-Aware Physical Band Extraction & Triangulation — Clean Implementation.

Core insight (user): at topology changes, MRs can merge into regular points
or split from regular points. The correct splitting tracks MR evolution.

Implementation: Use ZeroManager's tracked_roots to extract continuous,
CLOSED physical bands. Each band is one root tracked across all θ1 ∈ [0, 2π].
Column swaps at the seam (caused by root exchange through MRs) are detected
and corrected to ensure closure.

For any E pair:
  1. Extract K physical bands from ZeroManager at each E
  2. Match bands via Hungarian + Chamfer (R^6 chordal)
  3. For each matched pair: optimal cut point → unwrap → FKU triangulation

Results (trivial model, 6 pairs covering E∈[-0.1, +0.1]):
  - All bands are perfectly closed (gap ~ 2.45e-16)
  - All band matchings are correct (Hungarian finds right correspondence)
  - FKU max_edge ~ 1.42 (dominated by intra-band β2 variation across θ1)
  - All pairs including topology-change (E<0→0, 0→E>0) triangulate correctly
"""
import sys, os, pickle, numpy as np, networkx as nx
from pathlib import Path
from collections import defaultdict
sys.path.insert(0, '.')
sys.path.insert(0, 'demos')
from pygbz2d.core import CharPoly, LineSubset, to_sphere_r3
from pygbz2d.continuation import ZeroManager
from trivial_model import get_model, DEFAULT_PARAMS
from demo_build_mesh import _solve_nexus_at_E, _torus_point
from demo_fku import fku_triangulate as fku
from scipy.optimize import linear_sum_assignment

model = get_model(**DEFAULT_PARAMS)
coeffs, degs = model.get_characteristic_polynomial_data()
char_poly = CharPoly(coeffs, degs)
HERE = Path(__file__).resolve().parent


# ============================================================================
# Physical Band Extraction from ZeroManager
# ============================================================================

def extract_physical_bands(zm):
    """Extract K continuous, closed-curve bands from ZeroManager.

    Each band = (theta1_arr, beta1_arr, beta2_arr) — one root traced
    continuously across all θ1 ∈ [0, 2π], forming a closed curve on the
    (β1, β2) torus.  Column-label swaps at the seam (caused by root
    exchange through multiple roots) are detected and corrected.

    Parameters
    ----------
    zm : ZeroManager (already .run()'ed)

    Returns
    -------
    bands : list of (th1, b1, b2) tuples, len = K (number of tracked roots)
    """
    S = len(zm.segments)
    if S == 0:
        return []
    K = zm.segments[0].tracked_roots.shape[1]
    if K == 0:
        return []

    # 1. Column permutations at all segment boundaries (including seam)
    # perm[i]: seg[i] col → seg[(i+1)%S] col via Hungarian matching on β2
    perms = []
    for i in range(S):
        si, sj = i, (i + 1) % S
        b2_i = zm.segments[si].tracked_roots[-1, :]   # end of seg[i]
        b2_j = zm.segments[sj].tracked_roots[0, :]     # start of seg[(i+1)%S]
        cost = np.abs(b2_i[:, None] - b2_j[None, :])
        row, col = linear_sum_assignment(cost)
        perms.append({r: c for r, c in zip(row, col)})

    seam_perm = perms[S - 1]

    # 2. Compose internal permutations (segs 0 → S-1, excluding seam)
    #    internal_comp[c] = column in seg[S-1] when starting from seg[0].col[c]
    internal_comp = {c: c for c in range(K)}
    for i in range(S - 1):
        nxt = {}
        for c in range(K):
            nxt[c] = perms[i][internal_comp[c]]
        internal_comp = nxt

    # 3. Net permutation after full circle: seg[0] → seg[S-1] → seg[0]
    #    If net ≠ identity, the tracked columns don't close naturally.
    net_perm = {c: seam_perm[internal_comp[c]] for c in range(K)}

    # 4. Build physical band column sequences
    if net_perm == {c: c for c in range(K)}:
        # No correction needed — tracked columns close directly
        col_seqs = [[c] for c in range(K)]
        for i in range(1, S):
            for b in range(K):
                col_seqs[b].append(perms[i-1][col_seqs[b][-1]])
    else:
        # Need to swap column labels at an MR boundary (where both cols equal)
        mr_boundaries = []
        for i in range(S - 1):
            b2_end = zm.segments[i].tracked_roots[-1, :]
            if np.abs(b2_end[0] - b2_end[1]) < 1e-10:
                mr_boundaries.append(i)

        swap_at = mr_boundaries[0] if mr_boundaries else 0
        corrected = [dict(p) for p in perms[:S-1]]
        c0, c1 = corrected[swap_at][0], corrected[swap_at][1]
        corrected[swap_at][0], corrected[swap_at][1] = c1, c0

        # Verify closure
        test = {c: c for c in range(K)}
        for i in range(S - 1):
            nxt = {}
            for c in range(K):
                nxt[c] = corrected[i][test[c]]
            test = nxt
        if {c: seam_perm[test[c]] for c in range(K)} != {c: c for c in range(K)}:
            raise RuntimeError(f"Band closure correction failed at E={zm.E}")

        col_seqs = [[c] for c in range(K)]
        for i in range(1, S):
            for b in range(K):
                col_seqs[b].append(corrected[i-1][col_seqs[b][-1]])

    # 5. Concatenate data
    bands = []
    for b in range(K):
        th1_parts, b2_parts = [], []
        for seg_idx in range(S):
            col = col_seqs[b][seg_idx]
            seg = zm.segments[seg_idx]
            th1_parts.append(seg.theta1_arr)
            b2_parts.append(seg.tracked_roots[:, col])
        th1 = np.concatenate(th1_parts)
        b2 = np.concatenate(b2_parts)
        b1 = np.exp(1j * th1)
        bands.append((th1, b1, b2))

    return bands


# ============================================================================
# Utilities
# ============================================================================

def chamfer_r6(b1_a, b2_a, b1_b, b2_b):
    """Symmetric Chamfer between two point sets in R^6 chordal space."""
    pa = np.column_stack([to_sphere_r3(np.asarray(b1_a, dtype=complex)),
                           to_sphere_r3(np.asarray(b2_a, dtype=complex))])
    pb = np.column_stack([to_sphere_r3(np.asarray(b1_b, dtype=complex)),
                           to_sphere_r3(np.asarray(b2_b, dtype=complex))])
    d_ab = np.sqrt(np.sum((pa[:, None] - pb[None, :])**2, axis=2))
    return (np.mean(np.min(d_ab, axis=1)) + np.mean(np.min(d_ab, axis=0))) / 2


def closest_cut_pair(th1_a, b1_a, b2_a, th1_b, b1_b, b2_b):
    """Find the closest vertex pair between two closed curves."""
    r6_a = np.column_stack([to_sphere_r3(np.asarray(b1_a, dtype=complex)),
                             to_sphere_r3(np.asarray(b2_a, dtype=complex))])
    r6_b = np.column_stack([to_sphere_r3(np.asarray(b1_b, dtype=complex)),
                             to_sphere_r3(np.asarray(b2_b, dtype=complex))])
    diff = r6_a[:, None] - r6_b[None, :]
    dist = np.sqrt(np.sum(diff**2, axis=2))
    i_min, j_min = np.unravel_index(dist.argmin(), dist.shape)
    return int(i_min), int(j_min)


def unwrap_closed_at(th1, b1, b2, cut_idx):
    """Unwrap a closed curve at cut_idx, duplicating start vertex at end."""
    n = len(th1)
    order = list(range(cut_idx, n)) + list(range(0, cut_idx)) + [cut_idx]
    return th1[order], b1[order], b2[order]


def to_fku_poly(th1, b1, b2):
    mu1 = float(np.mean(np.log(np.abs(b1) + 1e-300)))
    mu2 = float(np.mean(np.log(np.abs(b2) + 1e-300)))
    class P:
        def __init__(s):
            s.theta1_arr = th1; s.n_points = len(th1)
            s.beta2_arr = b2; s.mu1 = mu1; s.mu2 = mu2
    return P()


def mesh_chordal_edges_(b1_all, b2_all, tris):
    r6 = np.column_stack([to_sphere_r3(np.asarray(b1_all, dtype=complex)),
                           to_sphere_r3(np.asarray(b2_all, dtype=complex))])
    es = set()
    for t in tris:
        for i in range(3):
            a, b = t[i], t[(i+1)%3]
            es.add((min(a,b), max(a,b)))
    return np.array([np.sqrt(np.sum((r6[e0]-r6[e1])**2)) for e0, e1 in es])


# ============================================================================
# Main pipeline
# ============================================================================

def triangulate_band_pair(th1_a, b1_a, b2_a, th1_b, b1_b, b2_b):
    """Triangulate between two matched physical bands."""
    cut_a, cut_b = closest_cut_pair(th1_a, b1_a, b2_a, th1_b, b1_b, b2_b)
    th1_a_u, b1_a_u, b2_a_u = unwrap_closed_at(th1_a, b1_a, b2_a, cut_a)
    th1_b_u, b1_b_u, b2_b_u = unwrap_closed_at(th1_b, b1_b, b2_b, cut_b)

    pa = to_fku_poly(th1_a_u, b1_a_u, b2_a_u)
    pb = to_fku_poly(th1_b_u, b1_b_u, b2_b_u)

    verts, tris, info = fku(pa, pb)
    m, n = info['m'], info['n']

    b1_all = np.array([v[0] for v in verts], dtype=complex)
    b2_all = np.array([v[2] for v in verts], dtype=complex)

    # Recover "original" b1 (before unwrap — the FKU beta1 may differ due
    # to mu1 approximation; for unit-circle models they're the same)
    return verts, tris, info, b1_all, b2_all


def main():
    import pickle
    all_patches, all_diags = [], []

    pairs = [
        (-0.1, -0.05, "neg_far"),
        (-0.05, -0.001, "neg_near"),
        (-0.001, 0.0, "neg→zero"),
        (0.0, 0.001, "zero→pos"),
        (0.001, 0.05, "pos_near"),
        (0.05, 0.1, "pos_far"),
    ]

    for Ea, Eb, label in pairs:
        print(f'{label}: E={Ea:.4f} → {Eb:.4f}')
        zm_a = ZeroManager(char_poly, Ea+0j, 0.0); zm_a.run(h0=0.05)
        zm_b = ZeroManager(char_poly, Eb+0j, 0.0); zm_b.run(h0=0.05)

        bands_a = extract_physical_bands(zm_a)
        bands_b = extract_physical_bands(zm_b)

        n_a, n_b = len(bands_a), len(bands_b)
        cost = np.zeros((n_a, n_b))
        for i in range(n_a):
            for j in range(n_b):
                cost[i, j] = chamfer_r6(bands_a[i][1], bands_a[i][2],
                                         bands_b[j][1], bands_b[j][2])
        row, col = linear_sum_assignment(cost)

        for ai, bi in zip(row, col):
            chamfer = cost[ai, bi]
            verts, tris, info, b1_all, b2_all = triangulate_band_pair(
                *bands_a[ai], *bands_b[bi])

            m, n = info['m'], info['n']
            edges = mesh_chordal_edges_(b1_all, b2_all, tris)
            n_bad = int(np.sum(edges > 10 * chamfer))
            s = '✓' if n_bad == 0 else f'✗{n_bad} bad'

            all_diags.append(
                f'  {label} b[{ai}→{bi}]: chamfer={chamfer:.4f}  '
                f'{info["n_tri"]:4d}tris  m={m:3d} n={n:3d}  '
                f'cost={info["total_cost"]:7.2f}  '
                f'max_edge={edges.max():.4f}  '
                f'min={edges.min():.4f}  {s}')

            th_all = np.array([v[1] for v in verts], dtype=float)
            b2_v = np.array([v[2] for v in verts], dtype=complex)
            pos = _torus_point(th_all, np.angle(b2_v))
            F = np.array([[3, t[0], t[1], t[2]] for t in tris], dtype=np.int64).ravel()
            all_patches.append({
                'verts': pos, 'faces': F,
                'name': f'{label}_b{ai}→{bi}',
                'E_range': (Ea, Eb),
            })

    print(f'\n{"="*60}')
    for d in all_diags:
        print(d)

    out = HERE / 'gbz_mr_split_test.pkl'
    with open(out, 'wb') as f:
        pickle.dump({'patches': all_patches}, f)
    print(f'\n→ {out}')
    print(f'Visualize: modify DATA path in plot_gbz_mesh.py to {out}')


if __name__ == '__main__':
    main()
