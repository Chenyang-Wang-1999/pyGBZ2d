'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-27
Copyright © Department of Physics, Tsinghua University. All rights reserved

2D GBZ mesh via nexus + FKU, with bisection bridging for topology changes.

Pipeline:
  1. Build nexuses at each E from the cached trivial_model scan.
  2. Same-topology pairs → FKU closed-cycle triangulation → stitch along E.
  3. Topology-change pairs (different component counts) → bisect to narrow
     the interval, collect intermediate slices, stitch them, then bridge the
     narrow gap via MSS (Meyers-Skinner-Sloan branching).
  4. Export mesh for PyVista visualisation.

Acceptance criteria:
  - No new vertices: every triangle vertex is an original LineSubset sample.
  - Boundary = two nexus curves: no extra boundary edges.
  - Max chordal edge ~ Chamfer distance.

Distance metric:
  chordal(beta1_a, beta1_b) + chordal(beta2_a, beta2_b)
  = Euclidean distance on the Riemann sphere (stereographic projection)
    summed over the two beta components.
'''

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import networkx as nx

from gbz_types import LineSubset, to_sphere_r3, TWO_PI
from demo_nexus import build_nexus, nexus_info
from demo_fku import fku_triangulate


# ---------------------------------------------------------------------------
# Torus helpers
# ---------------------------------------------------------------------------

TORUS_R = 3.0
TORUS_r = 1.0


def _torus_point(theta1, theta2):
    t1 = np.asarray(theta1, dtype=float) % (TWO_PI)
    t2 = np.asarray(theta2, dtype=float) % (TWO_PI)
    x = (TORUS_R + TORUS_r * np.cos(t2)) * np.cos(t1)
    y = (TORUS_R + TORUS_r * np.cos(t2)) * np.sin(t1)
    z = TORUS_r * np.sin(t2)
    return np.stack([x, y, z], axis=-1)


# ---------------------------------------------------------------------------
# Nexus re-solving (for bisection midpoints)
# ---------------------------------------------------------------------------

def _solve_nexus_at_E(coeffs, degs, E):
    """Re-solve GBZ at energy E via amoeba, build and return the nexus."""
    import brute_force_amoeba as bfa

    res = bfa.collect_GBZ_subsets(coeffs, degs, E, 0.0, debug_mode=False)
    if not res.success or res.index == (0, 0):
        return None

    lines = [s for s in res.subsets if isinstance(s, LineSubset)]
    if not lines:
        return None

    lines.sort(key=lambda L: float(np.median(
        np.angle(np.asarray(L.beta2_arr, dtype=complex)) % (TWO_PI))))
    return build_nexus(lines, tol=1e-12)


# ---------------------------------------------------------------------------
# Cycle extraction
# ---------------------------------------------------------------------------

def _cycle_vertex_order(G, start_node=None):
    """Walk the cycle: line_id=0 edges first (branch 0), then rest (branch 1 rev)."""
    if start_node is None:
        start_node = list(G.nodes())[0]

    n_nodes = G.number_of_nodes()
    order = [start_node]
    visited = {start_node}
    current = start_node

    while True:
        nxt = None
        for nb in G.neighbors(current):
            if nb not in visited:
                edge = G.get_edge_data(current, nb)
                if edge is not None and edge.get('line_id') == 0:
                    nxt = nb
                    break
        if nxt is None:
            break
        visited.add(nxt)
        order.append(nxt)
        current = nxt

    while len(visited) < n_nodes:
        for nb in G.neighbors(current):
            if nb not in visited:
                visited.add(nb)
                order.append(nb)
                current = nb
                break

    return order


def _cycle_to_linesubset(G, node_order):
    th = np.array([G.nodes[n]['theta1'] for n in node_order], dtype=float)
    b2 = np.array([G.nodes[n]['beta2'] for n in node_order], dtype=complex)
    return LineSubset(E=G.graph['E'], mu1=G.graph['mu1'],
                      theta1_arr=th, beta2_arr=b2)


# ---------------------------------------------------------------------------
# FKU between two closed cycles
# ---------------------------------------------------------------------------

def _closest_vertex_pair(G_a, G_b):
    """Find the pair of nodes (one from each cycle) with minimum chordal
    distance in (beta1, beta2) space.

    This gives a purely distance-based cut-point correspondence without
    any dependence on theta1 wrap or merge-node heuristics.
    """
    nodes_a = list(G_a.nodes())
    nodes_b = list(G_b.nodes())

    # Sample to keep O(N²) tractable (~300² = 90k is fine).
    b1_a = np.array([G_a.nodes[n]['beta1'] for n in nodes_a])
    b2_a = np.array([G_a.nodes[n]['beta2'] for n in nodes_a])
    b1_b = np.array([G_b.nodes[n]['beta1'] for n in nodes_b])
    b2_b = np.array([G_b.nodes[n]['beta2'] for n in nodes_b])

    # Block-wise to avoid O(N²) memory: compute min per A-vertex.
    best_d = float('inf')
    best_pair = (nodes_a[0], nodes_b[0])
    block = 64
    for i_start in range(0, len(nodes_a), block):
        i_end = min(i_start + block, len(nodes_a))
        d_b1 = np.linalg.norm(
            to_sphere_r3(b1_a[i_start:i_end])[:, None, :]
            - to_sphere_r3(b1_b)[None, :, :], axis=2)
        d_b2 = np.linalg.norm(
            to_sphere_r3(b2_a[i_start:i_end])[:, None, :]
            - to_sphere_r3(b2_b)[None, :, :], axis=2)
        d = d_b1 + d_b2
        ij = np.unravel_index(np.argmin(d), d.shape)
        d_min = float(d[ij])
        if d_min < best_d:
            best_d = d_min
            best_pair = (nodes_a[i_start + ij[0]], nodes_b[ij[1]])

    return best_pair[0], best_pair[1], best_d


def _closest_vertex_pair_for_node(G_a, G_b, node_a):
    """Find the node in G_b closest (chordal) to a specific node in G_a."""
    b1_a = G_a.nodes[node_a]['beta1']
    b2_a = G_a.nodes[node_a]['beta2']
    r3_a_b1 = to_sphere_r3(np.array([b1_a]))
    r3_a_b2 = to_sphere_r3(np.array([b2_a]))

    nodes_b = list(G_b.nodes())
    b1_b = np.array([G_b.nodes[n]['beta1'] for n in nodes_b])
    b2_b = np.array([G_b.nodes[n]['beta2'] for n in nodes_b])
    r3_b_b1 = to_sphere_r3(b1_b)
    r3_b_b2 = to_sphere_r3(b2_b)

    d_b1 = np.linalg.norm(r3_a_b1 - r3_b_b1, axis=1)
    d_b2 = np.linalg.norm(r3_a_b2 - r3_b_b2, axis=1)
    d = d_b1 + d_b2
    j = int(np.argmin(d))
    return nodes_b[j], float(d.min()), float(d[j])


# ---------------------------------------------------------------------------
# FKU between two closed cycles
# ---------------------------------------------------------------------------

def triangulate_two_nexuses(G_a, G_b, cut_a=None):
    """Triangulate between two closed-cycle nexuses via FKU.

    Cut points default to the closest vertex pair in chordal distance.
    If ``cut_a`` is provided, it overrides the A-side cut point — used
    during stitching to ensure the A-cycle vertex order matches the
    previous pair's B-cycle order (same nexus, same cut point).

    Cycle order is determined by the deterministic line_id=0-first walk.
    """
    if cut_a is not None and cut_a in G_a:
        cut_b, _, cut_dist = _closest_vertex_pair_for_node(G_a, G_b, cut_a)
    else:
        cut_a, cut_b, cut_dist = _closest_vertex_pair(G_a, G_b)

    order_a = _cycle_vertex_order(G_a, start_node=cut_a)
    order_b = _cycle_vertex_order(G_b, start_node=cut_b)

    m_orig, n_orig = len(order_a), len(order_b)

    La = _cycle_to_linesubset(G_a, order_a + [order_a[0]])
    Lb = _cycle_to_linesubset(G_b, order_b + [order_b[0]])

    fku_verts, fku_tris, fku_info = fku_triangulate(La, Lb)

    m_full, n_full = m_orig + 1, n_orig + 1

    all_th1 = np.array([v[1] for v in fku_verts], dtype=float)
    all_b2 = np.array([v[2] for v in fku_verts], dtype=complex)
    all_b1 = np.exp(La.mu1 + 1j * all_th1)

    # Strip duplicate vertices, remap tri indices.
    remap = {}
    for old in range(m_full + n_full):
        if old == m_orig:
            new = 0
        elif old == m_full + n_orig:
            new = m_orig
        elif old < m_orig:
            new = old
        else:
            new = old - 1
        remap[old] = new

    tris_new = np.array([[remap[int(x)] for x in t] for t in fku_tris], dtype=int)

    pos_all = _torus_point(all_th1, np.angle(all_b2))
    verts_A = pos_all[:m_orig]
    verts_B = pos_all[m_full:m_full + n_orig]
    b1_A, b2_A = all_b1[:m_orig], all_b2[:m_orig]
    b1_B, b2_B = all_b1[m_full:m_full + n_orig], all_b2[m_full:m_full + n_orig]
    b1_all_out = np.concatenate([b1_A, b1_B])
    b2_all_out = np.concatenate([b2_A, b2_B])

    info = {
        'm': m_orig, 'n': n_orig, 'n_tri': len(tris_new),
        'total_cost': fku_info['total_cost'], 'cut_dist': cut_dist,
        'cut_b_node': cut_b,
        'cut_a_th1': float(G_a.nodes[cut_a]['theta1']),
        'cut_b_th1': float(G_b.nodes[cut_b]['theta1']),
    }
    return verts_A, verts_B, tris_new, b1_all_out, b2_all_out, info


# ---------------------------------------------------------------------------
# Chordal diagnostics
# ---------------------------------------------------------------------------

def chamfer_distance(b1_a, b2_a, b1_b, b2_b):
    r3_b1_a, r3_b2_a = to_sphere_r3(b1_a), to_sphere_r3(b2_a)
    r3_b1_b, r3_b2_b = to_sphere_r3(b1_b), to_sphere_r3(b2_b)
    d = (np.linalg.norm(r3_b1_a[:, None] - r3_b1_b[None], axis=2) +
         np.linalg.norm(r3_b2_a[:, None] - r3_b2_b[None], axis=2))
    return 0.5 * (float(np.mean(np.min(d, axis=1))) +
                  float(np.mean(np.min(d, axis=0))))


def mesh_chordal_edges(b1_all, b2_all, tris):
    r3_b1 = to_sphere_r3(b1_all)
    r3_b2 = to_sphere_r3(b2_all)
    edges = []
    for t in tris:
        for a, b in [(0, 1), (1, 2), (0, 2)]:
            i, j = int(t[a]), int(t[b])
            edges.append(float(np.linalg.norm(r3_b1[i] - r3_b1[j]) +
                               np.linalg.norm(r3_b2[i] - r3_b2[j])))
    return np.array(edges)


# ---------------------------------------------------------------------------
# MSS branching — split a single cycle at merge nodes, match to sub-cycles
# ---------------------------------------------------------------------------

def mss_split_cycle(G_single, G_multi_cycles):
    """Split a single-cycle nexus at its merge nodes into arcs,
    match each arc to a sub-cycle of a multi-cycle nexus by mean theta2.

    Parameters:
        G_single: single-cycle nexus.
        G_multi_cycles: dict {comp_id: subgraph} from the multi-cycle nexus.

    Returns:
        List of (arc_line, sub_line, comp_id) tuples ready for FKU.
    """
    merge_nodes = sorted(
        n for n, d in G_single.nodes(data=True)
        if len(d.get('line_ids', set())) > 1
    )
    if len(merge_nodes) < 2:
        raise ValueError(f"Need ≥2 merge nodes to split, got {len(merge_nodes)}")

    n0, n1 = merge_nodes[0], merge_nodes[1]
    full_order = _cycle_vertex_order(G_single)

    pos0, pos1 = full_order.index(n0), full_order.index(n1)
    if pos0 < pos1:
        arc1_order = full_order[pos0:pos1 + 1]
        arc2_order = full_order[pos1:] + full_order[:pos0 + 1]
    else:
        arc1_order = full_order[pos0:] + full_order[:pos1 + 1]
        arc2_order = full_order[pos1:pos0 + 1]

    if len(arc1_order) < 2 or len(arc2_order) < 2:
        raise ValueError(f"Degenerate arcs: {len(arc1_order)}, {len(arc2_order)}")

    arc1_line = _cycle_to_linesubset(G_single, arc1_order)
    arc2_line = _cycle_to_linesubset(G_single, arc2_order)
    arcs = [arc1_line, arc2_line]

    # Match arcs to sub-cycles by chordal Chamfer distance.
    sub_cycles = list(G_multi_cycles.items())
    sub_lines = []
    for _, sub_G in sub_cycles:
        sub_order = _cycle_vertex_order(sub_G)
        sub_lines.append(_cycle_to_linesubset(sub_G, sub_order))

    # Build cost matrix: Chamfer distance between each arc and sub-cycle.
    n_arcs, n_subs = len(arcs), len(sub_lines)
    cost = np.zeros((n_arcs, n_subs))
    for i, arc in enumerate(arcs):
        b1_a = np.exp(arc.mu1 + 1j * np.asarray(arc.theta1_arr, dtype=float))
        b2_a = np.asarray(arc.beta2_arr, dtype=complex)
        r3_b1_a = to_sphere_r3(b1_a)
        r3_b2_a = to_sphere_r3(b2_a)
        for j, sub_line in enumerate(sub_lines):
            b1_s = np.exp(sub_line.mu1 + 1j * np.asarray(sub_line.theta1_arr, dtype=float))
            b2_s = np.asarray(sub_line.beta2_arr, dtype=complex)
            r3_b1_s = to_sphere_r3(b1_s)
            r3_b2_s = to_sphere_r3(b2_s)
            d = (np.linalg.norm(r3_b1_a[:, None] - r3_b1_s[None], axis=2) +
                 np.linalg.norm(r3_b2_a[:, None] - r3_b2_s[None], axis=2))
            cost[i, j] = 0.5 * (float(np.mean(np.min(d, axis=1))) +
                                float(np.mean(np.min(d, axis=0))))

    from scipy.optimize import linear_sum_assignment
    row, col = linear_sum_assignment(cost)

    result = []
    for i, j in zip(row, col):
        comp_id, sub_G = sub_cycles[j]
        result.append((arcs[i], sub_lines[j], comp_id))
    return result


# ---------------------------------------------------------------------------
# Bisection bridge across topology change
# ---------------------------------------------------------------------------

def bisect_topology_change(coeffs, degs, Ea, Ga, Eb, Gb,
                           tol_rel=0.05, max_iter=30):
    """Bisect the E-interval [Ea, Eb] where nexus topology changes.

    At each midpoint the GBZ is re-solved and a nexus built.  Intermediate
    slices on each side are collected.

    Returns:
        a_slices: list of (E, G) on the Ea side, E-ascending toward E*.
        b_slices: list of (E, G) on the Eb side, E-descending toward E*.
    """
    tol = tol_rel * abs(Eb - Ea)
    comp_a = nx.number_connected_components(Ga)
    comp_b = nx.number_connected_components(Gb)

    a_slices = []
    b_slices = []
    Elo, Ehi = Ea, Eb

    for _ in range(max_iter):
        if Ehi - Elo <= tol:
            break

        Emid = 0.5 * (Elo + Ehi)
        Gmid = _solve_nexus_at_E(coeffs, degs, Emid)
        if Gmid is None:
            # Outside spectrum — narrow the other direction.
            Ehi = Emid
            continue

        comp_mid = nx.number_connected_components(Gmid)

        if comp_mid == comp_a:
            a_slices.append((Emid, Gmid))
            Elo = Emid
        elif comp_mid == comp_b:
            b_slices.append((Emid, Gmid))
            Ehi = Emid
        else:
            # Intermediate topology not matching either side —
            # narrow further from the Eb side.
            Ehi = Emid

    a_slices.sort(key=lambda x: x[0])
    b_slices.sort(key=lambda x: x[0], reverse=True)
    return a_slices, b_slices


# ---------------------------------------------------------------------------
# Local FKU helper (add a single pair's triangles to a patch buffer)
# ---------------------------------------------------------------------------

def _emit_fku_pair(verts_list, faces_list, offset, prev_B_idx,
                   Ga, Gb, cut_a=None):
    """Triangulate Ga→Gb, append to verts/faces lists, return updated state.
    Returns (offset, b_idx, fku_info, cut_b) where cut_b propagates to the
    next pair's A-side."""
    vA, vB, tris, _b1, _b2, fku_info = triangulate_two_nexuses(Ga, Gb, cut_a=cut_a)
    m, n = fku_info['m'], fku_info['n']

    if prev_B_idx is None:
        a_idx = list(range(offset, offset + m))
        verts_list.append(vA)
        offset += m
    else:
        a_idx = prev_B_idx
        if len(a_idx) != m:
            raise ValueError(
                f"Vertex count mismatch during stitch: "
                f"prev_B has {len(a_idx)}, new A has {m}")

    b_idx = list(range(offset, offset + n))
    verts_list.append(vB)
    offset += n

    for t in tris:
        gi = [a_idx[x] if x < m else b_idx[x - m] for x in t]
        faces_list.append([3, gi[0], gi[1], gi[2]])

    return offset, b_idx, fku_info, fku_info.get('cut_b_node')


# ---------------------------------------------------------------------------
# Stitching across E
# ---------------------------------------------------------------------------

def stitch_all(nexuses, E_sorted, coeffs, degs):
    """Stitch all adjacent E pairs into mesh patches.

    Same-topology pairs → FKU.
    Topology-change pairs → bisection bridge + MSS.

    Returns:
        patches: list of dicts with 'verts', 'faces', 'E_range', 'name'.
        diagnostics: list of strings.
    """
    patches = []
    diag = []
    verts_list = []
    faces_list = []
    offset = 0
    prev_B_idx = None
    prev_cut_b = None  # propagate A-side cut point for stitching
    seg_E_start = None
    seg_name = None

    def _flush():
        nonlocal verts_list, faces_list, offset, prev_B_idx, prev_cut_b, seg_E_start, seg_name
        if verts_list:
            V = np.concatenate(verts_list, axis=0) if len(verts_list) > 1 else verts_list[0]
            F = np.array([x for f in faces_list for x in f], dtype=np.int64)
            patches.append({
                'verts': V, 'faces': F,
                'E_range': (seg_E_start, seg_E_end),
                'name': seg_name,
            })
        verts_list, faces_list = [], []
        offset = 0
        prev_B_idx = None
        prev_cut_b = None
        seg_E_start = None

    seg_E_end = None
    k = 0
    while k < len(E_sorted) - 1:
        Ea, Eb = E_sorted[k], E_sorted[k + 1]
        Ga, info_a = nexuses[Ea]
        Gb, info_b = nexuses[Eb]

        if info_a['n_components'] == info_b['n_components']:
            # Same topology: plain FKU.
            if seg_E_start is None:
                seg_E_start = Ea
                seg_name = f'E<0' if Ea < 0 else f'E>0'

            offset, prev_B_idx, fku_info, cut_b = _emit_fku_pair(
                verts_list, faces_list, offset, prev_B_idx, Ga, Gb,
                cut_a=prev_cut_b)
            prev_cut_b = cut_b
            seg_E_end = Eb
            diag.append(
                f"  {seg_name} E={Ea:+.2f}→{Eb:+.2f}  "
                f"|A|={fku_info['m']} |B|={fku_info['n']}  "
                f"tris={fku_info['n_tri']} cost={fku_info['total_cost']:.3f}")
            k += 1

        else:
            # Topology change: bisect + bridge.
            diag.append(
                f"\n  [bridge] E={Ea:+.2f} ({info_a['n_components']} comp) "
                f"→ E={Eb:+.2f} ({info_b['n_components']} comp) — bisecting...")

            a_slices, b_slices = bisect_topology_change(
                coeffs, degs, Ea, Ga, Eb, Gb, tol_rel=0.05)

            diag.append(
                f"    bisected: {len(a_slices)} slices on Ea side, "
                f"{len(b_slices)} slices on Eb side")

            # Flush current patch before the bridge.
            _flush()

            # Build bridge patches: Ea side approach.
            bridge_verts = []
            bridge_faces = []
            br_offset = 0
            br_prev_idx = None

            # Ea side: stitch from the cached slice toward E*.
            curr_G = Ga
            curr_E = Ea
            for Emid, Gmid in a_slices:
                br_offset, br_prev_idx, fku_info = _emit_fku_pair(
                    bridge_verts, bridge_faces, br_offset, br_prev_idx,
                    curr_G, Gmid)
                diag.append(
                    f"    bisect_FKU E={curr_E:+.4f}→{Emid:+.4f}  "
                    f"|A|={fku_info['m']} |B|={fku_info['n']}  "
                    f"cost={fku_info['total_cost']:.3f}")
                curr_G, curr_E = Gmid, Emid

            # MSS bridge at the narrowest gap.
            G_lo = curr_G
            G_hi = b_slices[-1][1] if b_slices else Gb
            n_lo_comp = nx.number_connected_components(G_lo)
            n_hi_comp = nx.number_connected_components(G_hi)

            # Determine which side is the single cycle.
            if n_lo_comp == 1 and n_hi_comp > 1:
                G_single, G_multi = G_lo, G_hi
            elif n_lo_comp > 1 and n_hi_comp == 1:
                G_single, G_multi = G_hi, G_lo
            else:
                diag.append(f"    MSS skip: unexpected comp counts "
                           f"({n_lo_comp}, {n_hi_comp})")
                _flush()
                k += 1
                continue

            # Build multi-cycle subgraphs dict.
            multi_subgraphs = {}
            for comp in nx.connected_components(G_multi):
                sub = G_multi.subgraph(comp)
                multi_subgraphs[list(comp)[0]] = sub

            try:
                mss_pairs = mss_split_cycle(G_single, multi_subgraphs)
            except ValueError as e:
                diag.append(f"    MSS failed: {e}")
                _flush()
                k += 1
                continue

            for arc_line, sub_line, comp_id in mss_pairs:
                # Determine order: if G_single is on the lo side, arc→sub;
                # if on the hi side, sub→arc.
                if G_single is G_lo:
                    La, Lb = arc_line, sub_line
                else:
                    La, Lb = sub_line, arc_line

                fku_verts, fku_tris, mss_info = fku_triangulate(La, Lb)
                m_mss, n_mss = mss_info['m'], mss_info['n']
                th_all = np.array([v[1] for v in fku_verts], dtype=float)
                b2_all = np.array([v[2] for v in fku_verts], dtype=complex)
                pos_torus = _torus_point(th_all, np.angle(b2_all))
                F = np.array([[3, t[0], t[1], t[2]]
                              for t in fku_tris], dtype=np.int64).ravel()
                patches.append({
                    'verts': pos_torus, 'faces': F,
                    'E_range': (Ea, Eb),
                    'name': f'bridge_MSS_comp{comp_id}',
                })
                diag.append(
                    f"    MSS comp {comp_id}: |arc|={m_mss} |sub|={n_mss}  "
                    f"tris={mss_info['n_tri']} cost={mss_info['total_cost']:.3f}")

            # Flush any accumulated bridge vertices (from the approach stitches).
            if bridge_verts:
                V = np.concatenate(bridge_verts, axis=0) if len(bridge_verts) > 1 else bridge_verts[0]
                F = np.array([x for f in bridge_faces for x in f], dtype=np.int64)
                patches.append({
                    'verts': V, 'faces': F,
                    'E_range': (Ea, a_slices[-1][0] if a_slices else Ea),
                    'name': f'bridge_approach_Ea_side',
                })

            # Eb side approach (stitch from b_slices back toward Eb).
            if b_slices:
                br2_verts = []
                br2_faces = []
                br2_offset = 0
                br2_prev_idx = None
                # b_slices are sorted E-descending, so first element is closest to E*.
                # We need E-ascending for stitching: reverse b_slices.
                b_ascending = list(reversed(b_slices))  # closest to E* first
                curr_G2 = b_ascending[0][1]
                curr_E2 = b_ascending[0][0]
                for Emid, Gmid in b_ascending[1:]:
                    br2_offset, br2_prev_idx, fku_info = _emit_fku_pair(
                        br2_verts, br2_faces, br2_offset, br2_prev_idx,
                        curr_G2, Gmid)
                    diag.append(
                        f"    bisect_FKU E={curr_E2:+.4f}→{Emid:+.4f}  "
                        f"|A|={fku_info['m']} |B|={fku_info['n']}  "
                        f"cost={fku_info['total_cost']:.3f}")
                    curr_G2, curr_E2 = Gmid, Emid
                # Stitch to Eb.
                br2_offset, br2_prev_idx, fku_info = _emit_fku_pair(
                    br2_verts, br2_faces, br2_offset, br2_prev_idx,
                    curr_G2, Gb)
                diag.append(
                    f"    bisect_FKU E={curr_E2:+.4f}→{Eb:+.4f}  "
                    f"|A|={fku_info['m']} |B|={fku_info['n']}  "
                    f"cost={fku_info['total_cost']:.3f}")
                if br2_verts:
                    V = np.concatenate(br2_verts, axis=0) if len(br2_verts) > 1 else br2_verts[0]
                    F = np.array([x for f in br2_faces for x in f], dtype=np.int64)
                    patches.append({
                        'verts': V, 'faces': F,
                        'E_range': (b_ascending[0][0], Eb),
                        'name': f'bridge_approach_Eb_side',
                    })

            k += 1

    _flush()
    return patches, diag


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    from trivial_line_cache import load

    data = load()
    E_grid = data['E_grid']
    results = data['results']

    # Build all nexuses from cache.
    print("Building nexuses...")
    nexuses = {}
    for i, r in enumerate(results):
        E = float(E_grid[i])
        if r.success and r.index == (0, 2):
            lines = [s for s in r.subsets if isinstance(s, LineSubset)]
            lines.sort(key=lambda L: float(np.median(
                np.angle(np.asarray(L.beta2_arr, dtype=complex)) % (TWO_PI))))
            G = build_nexus(lines, tol=1e-12)
            nexuses[E] = (G, nexus_info(G))

    E_sorted = sorted(nexuses.keys())

    # Get model for re-solving at bisection midpoints.
    from trivial_model import get_model, DEFAULT_PARAMS
    model = get_model(**DEFAULT_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()

    # Stitch everything.
    print("Stitching...")
    patches, diag = stitch_all(nexuses, E_sorted, coeffs, degs)

    for line in diag:
        print(line)

    total_verts = sum(len(p['verts']) for p in patches)
    total_tris = sum(len(p['faces']) // 4 for p in patches)
    print(f"\n{len(patches)} patch(es):")
    for p in patches:
        er = p['E_range']
        print(f"  {p['name']}: E=[{er[0]:+.4f}, {er[1]:+.4f}], "
              f"{len(p['verts'])} verts, {len(p['faces']) // 4} tris")
    print(f"Total: {total_verts} verts, {total_tris} tris")

    # Export.
    import pickle
    out = Path(__file__).resolve().parent / 'gbz_mesh_full.pkl'
    with open(out, 'wb') as f:
        pickle.dump({'patches': patches}, f)
    print(f"\nExported to {out}")


if __name__ == "__main__":
    main()
