'''
Bisection bridge + MSS for trivial-model E=0 topology change.

Step 1: bisect to narrow the gap.
Step 2: stitch approach segments with FKU.
Step 3: MSS — cut single cycle at merge nodes, match arcs to E=0 sub-cycles.
'''

import sys, pickle
from pathlib import Path
import numpy as np
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from gbz_types import LineSubset, to_sphere_r3, TWO_PI
from demo_nexus import build_nexus, nexus_info
from demo_fku import fku_triangulate
from demo_build_mesh import (
    _solve_nexus_at_E, triangulate_two_nexuses,
    bisect_topology_change, chamfer_distance, mesh_chordal_edges,
    _torus_point, _cycle_vertex_order, _cycle_to_linesubset, mss_split_cycle)
from trivial_line_cache import load
from trivial_model import get_model, DEFAULT_PARAMS


# ---------------------------------------------------------------------------
# FKU stitch helper
# ---------------------------------------------------------------------------

def stitch_slices(slices, name):
    verts_list, faces_list, diag = [], [], []
    offset, prev_B_idx, prev_cut_b = 0, None, None
    for k in range(len(slices) - 1):
        Ea, Ga = slices[k]
        Eb, Gb = slices[k + 1]
        vA, vB, tris, b1, b2, info = triangulate_two_nexuses(
            Ga, Gb, cut_a=prev_cut_b)
        m, n = info['m'], info['n']
        prev_cut_b = info['cut_b_node']

        chamfer = chamfer_distance(b1[:m], b2[:m], b1[m:], b2[m:])
        edges = mesh_chordal_edges(b1, b2, tris)
        n_bad = int(np.sum(edges > 10 * chamfer))

        if prev_B_idx is None:
            a_idx = list(range(offset, offset + m))
            verts_list.append(vA)
            offset += m
        else:
            a_idx = prev_B_idx

        b_idx = list(range(offset, offset + n))
        verts_list.append(vB)
        offset += n
        prev_B_idx = b_idx

        for t in tris:
            gi = [a_idx[x] if x < m else b_idx[x - m] for x in t]
            faces_list.append([3, gi[0], gi[1], gi[2]])

        status = '✓' if n_bad == 0 else f'✗ {n_bad} bad edges'
        diag.append(
            f'  {name} E={Ea:+.4f}→{Eb:+.4f}  '
            f'|A|={m:3d} |B|={n:3d}  tris={info["n_tri"]:3d}  '
            f'cost={info["total_cost"]:7.3f}  chamfer={chamfer:.4f}  '
            f'max_edge={edges.max():.4f}  {status}')

    V = np.concatenate(verts_list, axis=0) if len(verts_list) > 1 else verts_list[0]
    F = np.array([x for f in faces_list for x in f], dtype=np.int64)
    return V, F, diag


# ---------------------------------------------------------------------------
# MSS bridge: single cycle ↔ multi-cycle
# ---------------------------------------------------------------------------

def mss_bridge_patches(G_single, G_multi, E_lo, E_hi, tag):
    """Run MSS on a single-cycle ↔ multi-cycle nexus pair.

    Returns list of patch dicts and diagnostics strings.
    """
    # Build subgraph dict for the multi-cycle nexus.
    multi_subgraphs = {}
    for comp in nx.connected_components(G_multi):
        sub = G_multi.subgraph(comp)
        multi_subgraphs[list(comp)[0]] = sub

    mss_pairs = mss_split_cycle(G_single, multi_subgraphs)

    patches, diag = [], []
    for arc_line, sub_line, comp_id in mss_pairs:
        fku_verts, fku_tris, mss_info = fku_triangulate(arc_line, sub_line)
        m_mss, n_mss = mss_info['m'], mss_info['n']
        th_all = np.array([v[1] for v in fku_verts], dtype=float)
        b2_all = np.array([v[2] for v in fku_verts], dtype=complex)
        pos_torus = _torus_point(th_all, np.angle(b2_all))
        F = np.array([[3, t[0], t[1], t[2]]
                      for t in fku_tris], dtype=np.int64).ravel()
        patches.append({
            'verts': pos_torus, 'faces': F,
            'E_range': (E_lo, E_hi),
            'name': f'{tag}_comp{comp_id}',
        })
        diag.append(
            f'  MSS {tag} comp {comp_id}: |arc|={m_mss} |sub|={n_mss}  '
            f'tris={mss_info["n_tri"]} cost={mss_info["total_cost"]:.3f}')
    return patches, diag


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    data = load()
    E_grid, results = data['E_grid'], data['results']

    print('Building cached nexuses...')
    nexuses = {}
    for i, r in enumerate(results):
        E = float(E_grid[i])
        if r.success and r.index == (0, 2):
            lines = [s for s in r.subsets if isinstance(s, LineSubset)]
            lines.sort(key=lambda L: float(np.median(
                np.angle(np.asarray(L.beta2_arr, dtype=complex)) % (TWO_PI))))
            G = build_nexus(lines, tol=1e-12)
            nexuses[E] = (G, nexus_info(G))

    model = get_model(**DEFAULT_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()
    E_sorted = sorted(nexuses.keys())

    patches, all_diag = [], []

    # ---- E<0 ----
    print('Stitching E<0...')
    slices = [(E, nexuses[E][0]) for E in E_sorted if E < -0.005]
    V, F, diag = stitch_slices(slices, 'E<0')
    if len(V) > 0:
        patches.append({'verts': V, 'faces': F, 'name': 'E<0',
                        'E_range': (slices[0][0], slices[-1][0])})
    all_diag.extend(diag)

    # ---- E>0 ----
    print('Stitching E>0...')
    slices = [(E, nexuses[E][0]) for E in E_sorted if E > 0.005]
    V, F, diag = stitch_slices(slices, 'E>0')
    if len(V) > 0:
        patches.append({'verts': V, 'faces': F, 'name': 'E>0',
                        'E_range': (slices[0][0], slices[-1][0])})
    all_diag.extend(diag)

    # ---- Key E values ----
    E_neg = min(nexuses.keys(), key=lambda x: abs(x - (-0.1)))
    E_zero = min(nexuses.keys(), key=lambda x: abs(x - 0.0))
    E_pos = min(nexuses.keys(), key=lambda x: abs(x - 0.1))

    Ga_neg, _ = nexuses[E_neg]
    G_zero, _ = nexuses[E_zero]
    Gb_pos, _ = nexuses[E_pos]

    # ---- Bridge: E<0 → 0 side ----
    print('\nBisecting E<0 → 0...')
    a_slices, b_slices = bisect_topology_change(
        coeffs, degs, E_neg, Ga_neg, E_zero, G_zero, tol_rel=0.05)

    # Stitch Ea approach: cached E_neg → bisected slices.
    bridge_slices_neg = [(E_neg, Ga_neg)] + a_slices
    V, F, diag = stitch_slices(bridge_slices_neg, 'bridge_neg_approach')
    if len(V) > 0:
        patches.append({'verts': V, 'faces': F, 'name': 'bridge_neg_approach',
                        'E_range': (E_neg, bridge_slices_neg[-1][0])})
    all_diag.extend(diag)

    # MSS: last a_slice (single cycle, E≈-0.0031) → E=0 (two cycles).
    G_tip_neg = bridge_slices_neg[-1][1]
    E_tip_neg = bridge_slices_neg[-1][0]
    mss_p, mss_d = mss_bridge_patches(
        G_tip_neg, G_zero, E_tip_neg, E_zero, 'bridge_neg_MSS')
    patches.extend(mss_p)
    all_diag.extend(mss_d)

    # ---- Bridge: 0 → E>0 side ----
    print('\nBisecting 0 → E>0...')
    a_slices2, b_slices2 = bisect_topology_change(
        coeffs, degs, E_zero, G_zero, E_pos, Gb_pos, tol_rel=0.05)

    # MSS: E=0 (two cycles) → first b_slice (single cycle, E≈+0.0031).
    b_ascending = list(reversed(b_slices2))
    G_tip_pos = b_ascending[0][1]
    E_tip_pos = b_ascending[0][0]
    mss_p2, mss_d2 = mss_bridge_patches(
        G_tip_pos, G_zero, E_zero, E_tip_pos, 'bridge_pos_MSS')
    patches.extend(mss_p2)
    all_diag.extend(mss_d2)

    # Stitch Eb approach: b_slices → cached E_pos.
    bridge_slices_pos = b_ascending + [(E_pos, Gb_pos)]
    V2, F2, diag2 = stitch_slices(bridge_slices_pos, 'bridge_pos_approach')
    if len(V2) > 0:
        patches.append({'verts': V2, 'faces': F2, 'name': 'bridge_pos_approach',
                        'E_range': (E_tip_pos, E_pos)})
    all_diag.extend(diag2)

    # ---- Report ----
    for d in all_diag:
        print(d)

    total_verts = sum(len(p['verts']) for p in patches)
    total_tris = sum(len(p['faces']) // 4 for p in patches)
    print(f'\n{len(patches)} patches:')
    for p in patches:
        er = p['E_range']
        print(f'  {p["name"]}: E=[{er[0]:+.4f}, {er[1]:+.4f}]  '
              f'{len(p["verts"])} verts  {len(p["faces"]) // 4} tris')
    print(f'Total: {total_verts} verts, {total_tris} tris')

    out = Path(__file__).resolve().parent / 'gbz_mesh_full.pkl'
    with open(out, 'wb') as f:
        pickle.dump({'patches': patches}, f)
    print(f'\nExported to {out}')


if __name__ == '__main__':
    main()
