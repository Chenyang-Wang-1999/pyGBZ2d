'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-27
Copyright © Department of Physics, Tsinghua University. All rights reserved

LineSubset Nexus Generation Demo
=================================

Step 1: Build a "nexus" — a single 1D graph from all LineSubsets at a given E,
by gluing equivalent (beta1, beta2) points within tolerance ~1e-12.

Step 2: Research summary — survey of algorithms for building 2D triangulations
between adjacent E nexuses (the classic "surface reconstruction from
cross-sections" problem).

----
Step 2 Research: 2D Mesh Between Adjacent E Nexuses
----------------------------------------------------

Acceptance Criteria (User Requirement)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Given two adjacent nexuses G_a (at E_a) and G_b (at E_b), the 2D mesh must:
1. **No new points**: Every vertex is an original LineSubset sample point.
   No interpolation, no resampling.
2. **Boundary = G_a ∪ G_b**: The mesh's 1D boundary is exactly the two nexus
   curves. All G_a edges are incident to exactly one triangle (on the "a-side"),
   all G_b edges to exactly one triangle (on the "b-side"), and all stitching
   edges (connecting a-vertex to b-vertex) are incident to exactly two triangles.

This is the classic **surface reconstruction from cross-sections** problem
(also called "contour stitching" or "tiling between two contours").

Problem Formalization
~~~~~~~~~~~~~~~~~~~~~
Given two 1D graphs G_a, G_b embedded in (beta1, beta2) 4D space, build a 2D
simplicial complex where:
- Vertices = V_a ∪ V_b  (no new points)
- Triangles have exactly one of three forms:
  - (a_i, a_{i+1}, b_j)  — one G_a edge + 2 stitching edges
  - (a_i, b_j, b_{j+1})  — one G_b edge + 2 stitching edges
  - (a_i, a_j, b_k)      — all 3 are stitching edges (only needed at branching)

Case 1: Both Nexuses Are Single Cycles (the trivial model E≠0 case)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Given two closed polylines P (m vertices) and Q (n vertices):

**Classic approach (Keppel 1975, FKU 1977):**
1. Choose a starting vertex pair (p_0, q_0) that genuinely correspond.
2. "Unroll" both cycles by duplicating the start vertex at the end:
   P' = [p_0, p_1, ..., p_{m-1}, p_0']  (m+1 vertices)
   Q' = [q_0, q_1, ..., q_{n-1}, q_0']  (n+1 vertices)
3. Apply FKU (optimal DP on the (m+1)×(n+1) vertex-pair grid).
4. The end-cap edges (p_0–q_0) and (p'_0–q'_0) are the SAME edge (same
   vertices), so they cancel. Boundary = P ∪ Q exactly.

Result: a closed tube (topologically S¹ × [0,1]). No new points.
This satisfies both acceptance criteria.

**Where to cut?** The cut point (p_0, q_0) should be a pair close in
(beta1, beta2). For the trivial model at E≠0, the natural cut points are
the glued endpoint vertices (beta2 = ±1), which are the degree-2 nodes in
the nexus formed by merging the two branches' endpoints.

**Available implementation:** demo_fku.fku_triangulate — already does the
optimal DP for open polylines. Adapt to closed cycles via the
duplicate-start-vertex trick described above.

Case 2: Topology Change (e.g., E crossing 0)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The nexus topology changes: E≠0 has a single cycle, E=0 has two crossing
cycles. This is the **branching problem** in cross-section reconstruction.

**Meyers-Skinner-Sloan (1992) approach:**
1. Detect branching by comparing component counts.
2. For 1-to-N branching: each of the N cycles in G_b corresponds to a
   contiguous arc of the single cycle in G_a.
3. Split the single cycle at the branch points (beta2 = ±1 for the trivial
   model), treat each arc as an open polyline matched to one of the N cycles.
4. Apply FKU to each arc-cycle pair.

Case 3: General Graphs with Branching
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For arbitrary nexus graph topologies (not yet encountered):
- **Correspondence**: Match components across E via Hungarian on centroids.
- **Tiling**: For each matched component pair, apply Case 1 or Case 2.
- **Unmatched components**: Creation (new at E_b) or annihilation (gone at
  E_a) — these are genuine GBZ-edge topological events.

Library Survey
~~~~~~~~~~~~~~
| Library        | Role | Notes |
|----------------|------|-------|
| networkx       | Graph representation, cycle detection | Core for nexus |
| scipy.optimize | Hungarian matching | Already used |
| scipy.spatial  | KDTree for proximity queries | For merge detection |
| scipy.spatial.Delaunay | NOT recommended — creates unconstrained topology | Would violate boundary criteria |
| trimesh        | Post-processing (hole/manifold check) | Optional |
| PyVista        | Visualization on torus | Already used |
| CGAL (Python)  | Production-quality reconstruction | Overkill for now |

**Bottom line:** No new dependencies needed beyond networkx. The core
triangulation algorithm (FKU) is already implemented in demo_fku.py.

Recommended Next Step
~~~~~~~~~~~~~~~~~~~~~
Implement Case 1 (single cycles → FKU + duplicate-start trick) first,
since the trivial model's E≠0 slices are all single cycles. This directly
generalizes the existing demo_build_mesh pipeline but without the
per-E LineSubset matching step — the nexus already unifies the topology.
'''

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pickle
import numpy as np
import networkx as nx
from scipy.spatial import KDTree

from gbz_types import LineSubset, PointSubset, GBZResult, to_sphere_r3, chordal_cost_matrix


# ---------------------------------------------------------------------------
# Nexus construction
# ---------------------------------------------------------------------------

def _all_vertices(lines: list[LineSubset]):
    """Collect all (theta1, beta1, beta2) vertices from a list of LineSubsets.

    Returns:
        theta1_all : (N,) float
        beta1_all  : (N,) complex
        beta2_all  : (N,) complex
        line_ids   : (N,) int — which LineSubset each vertex belongs to
        local_idx  : (N,) int — position within that LineSubset
    """
    parts = []
    for lid, L in enumerate(lines):
        n = len(L.theta1_arr)
        b1 = np.exp(L.mu1 + 1j * np.asarray(L.theta1_arr, dtype=float))
        parts.append((
            np.asarray(L.theta1_arr, dtype=float),
            b1,
            np.asarray(L.beta2_arr, dtype=complex),
            np.full(n, lid, dtype=int),
            np.arange(n, dtype=int),
        ))
    return tuple(np.concatenate([p[i] for p in parts]) for i in range(5))


def _find_merges(beta1_all, beta2_all, tol=1e-12):
    """Find pairs of vertices from DIFFERENT LineSubsets within chordal tol.

    Uses chordal distance in (beta1, beta2) 4D space:
        d(p, q) = chordal(b1_p, b1_q) + chordal(b2_p, b2_q)
    where chordal = Euclidean distance on the Riemann sphere (R^3).

    We build a KDTree in R^6 (concatenated sphere coordinates for beta1 and
    beta2), query all pairs within radius `tol`, then filter by the exact
    chordal-sum distance.  The R^6 Euclidean distance is:
        sqrt(chordal(beta1)^2 + chordal(beta2)^2) <= chordal(beta1) + chordal(beta2)
    so queries with radius tol are guaranteed to cover all true matches,
    and we filter out false positives afterwards.

    Returns:
        merges : list of (i, j, distance) tuples where i < j and distance < tol.
    """
    r3_b1 = to_sphere_r3(beta1_all)
    r3_b2 = to_sphere_r3(beta2_all)
    # Concatenate into R^6
    r6 = np.concatenate([r3_b1, r3_b2], axis=1)

    tree = KDTree(r6)
    # Query all pairs within tol.  p=2 (Euclidean) is the default.
    pairs = tree.query_ball_tree(tree, r=tol)

    merges = []
    for i, neighbors in enumerate(pairs):
        for j in neighbors:
            if j <= i:
                continue
            # Exact chordal-sum distance.
            d = (float(np.linalg.norm(r3_b1[i] - r3_b1[j]))
                 + float(np.linalg.norm(r3_b2[i] - r3_b2[j])))
            if d < tol:
                merges.append((i, j, d))
    return merges


def build_nexus(lines: list[LineSubset], tol=1e-12):
    """Build a nexus graph from a list of LineSubsets at the same E.

    Equivalent (beta1, beta2) points within chordal distance `tol` are
    merged into the same graph node.  LineSubset internal sample sequences
    become graph edges.

    Parameters:
        lines: LineSubset objects at the same (E, mu1).
        tol: Chordal-distance tolerance for point merging (default 1e-12).

    Returns:
        G: networkx.Graph with node attributes:
            - beta1, beta2 : complex
            - theta1       : float (representative)
            - line_ids     : set of LineSubset indices contributing to this node
            - local_indices: list of (line_id, local_idx) for provenance
          Edge attributes:
            - line_id : which LineSubset this edge came from
          Graph attributes:
            - E, mu1
            - merge_count : number of merges performed
    """
    if not lines:
        G = nx.Graph()
        G.graph['E'] = None
        G.graph['mu1'] = None
        G.graph['merge_count'] = 0
        return G

    E = lines[0].E
    mu1 = lines[0].mu1

    th_all, b1_all, b2_all, line_ids, local_idx = _all_vertices(lines)
    n_total = len(th_all)

    # Phase 1: find merges across DIFFERENT LineSubsets only.
    merges = _find_merges(b1_all, b2_all, tol=tol)

    # Phase 2: build union-find over vertices, merging those within tolerance.
    parent = np.arange(n_total)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for i, j, d in merges:
        union(i, j)

    # Phase 3: map each vertex to its canonical node ID (0..n_nodes-1).
    unique_parents = sorted(set(find(i) for i in range(n_total)))
    node_id_map = {p: idx for idx, p in enumerate(unique_parents)}
    vertex_to_node = np.array([node_id_map[find(i)] for i in range(n_total)])

    n_nodes = len(unique_parents)

    # Phase 4: accumulate per-node data.
    node_b1 = np.zeros(n_nodes, dtype=complex)
    node_b2 = np.zeros(n_nodes, dtype=complex)
    node_th = np.zeros(n_nodes, dtype=float)
    node_line_ids = [set() for _ in range(n_nodes)]
    node_count = np.zeros(n_nodes, dtype=int)

    for i in range(n_total):
        nid = vertex_to_node[i]
        node_b1[nid] += b1_all[i]
        node_b2[nid] += b2_all[i]
        node_th[nid] += th_all[i]
        node_line_ids[nid].add(int(line_ids[i]))
        node_count[nid] += 1

    # Canonical values = centroid of merged points (for a single point,
    # centroid = the point itself).
    node_b1 /= node_count
    node_b2 /= node_count
    node_th /= node_count

    # Phase 5: build the graph.
    G = nx.Graph()
    G.graph['E'] = E
    G.graph['mu1'] = mu1
    G.graph['merge_count'] = n_total - n_nodes

    for nid in range(n_nodes):
        G.add_node(nid,
                   beta1=complex(node_b1[nid]),
                   beta2=complex(node_b2[nid]),
                   theta1=float(node_th[nid]),
                   line_ids=node_line_ids[nid])

    # Add edges: for each LineSubset, connect consecutive sample points
    # (their canonical node IDs).
    added_edges = set()
    for lid, L in enumerate(lines):
        n = len(L.theta1_arr)
        # Get the canonical node IDs for this LineSubset's vertices (in order).
        mask = line_ids == lid
        v_in_order = vertex_to_node[mask]

        for k in range(n - 1):
            u, v = int(v_in_order[k]), int(v_in_order[k + 1])
            if u == v:
                continue  # self-loop from merging consecutive identical points
            key = (min(u, v), max(u, v))
            if key not in added_edges:
                G.add_edge(u, v, line_id=lid)
                added_edges.add(key)

    return G


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def nexus_info(G: nx.Graph):
    """Return a dict of diagnostics for a nexus graph."""
    comps = list(nx.connected_components(G))
    n_comps = len(comps)

    degrees = np.array([d for _, d in G.degree()])
    deg_hist = dict(zip(*np.unique(degrees, return_counts=True)))

    try:
        cycle_basis = nx.cycle_basis(G)
    except nx.NetworkXNoCycle:
        cycle_basis = []

    # Per-component analysis
    comp_details = []
    for comp in comps:
        sub = G.subgraph(comp)
        n_v = sub.number_of_nodes()
        n_e = sub.number_of_edges()
        if n_v == 1:
            ctype = "singleton"
        elif n_v == n_e and nx.is_connected(sub):
            # For a connected graph, n_v == n_e means exactly one cycle
            # (a unicyclic graph / simple cycle if all degrees == 2).
            if all(d == 2 for _, d in sub.degree()):
                ctype = "simple_cycle"
            else:
                ctype = "unicyclic"
        elif nx.is_tree(sub):
            ctype = "tree"
        else:
            ctype = "general"
        comp_details.append({
            'n_nodes': n_v, 'n_edges': n_e,
            'type': ctype,
            'theta1_range': (
                float(min(sub.nodes[n]['theta1'] for n in sub.nodes)),
                float(max(sub.nodes[n]['theta1'] for n in sub.nodes)),
            ),
        })

    return {
        'E': G.graph.get('E'),
        'n_raw_vertices': G.number_of_nodes() + G.graph.get('merge_count', 0),
        'n_nodes': G.number_of_nodes(),
        'n_edges': G.number_of_edges(),
        'merge_count': G.graph.get('merge_count', 0),
        'n_components': n_comps,
        'component_types': [c['type'] for c in comp_details],
        'components': comp_details,
        'degree_histogram': deg_hist,
        'n_cycles': len(cycle_basis),
        'cycle_lengths': [len(c) for c in cycle_basis],
    }


# ---------------------------------------------------------------------------
# Torus visualization
# ---------------------------------------------------------------------------

TORUS_R = 3.0
TORUS_r = 1.0


def _torus_point(theta1, theta2):
    """Map (theta1, theta2) → 3D torus coordinates (vectorised)."""
    t1 = np.asarray(theta1, dtype=float) % (2 * np.pi)
    t2 = np.asarray(theta2, dtype=float) % (2 * np.pi)
    x = (TORUS_R + TORUS_r * np.cos(t2)) * np.cos(t1)
    y = (TORUS_R + TORUS_r * np.cos(t2)) * np.sin(t1)
    z = TORUS_r * np.sin(t2)
    return np.stack([x, y, z], axis=-1)


def nexus_torus_plot(G: nx.Graph, plotter=None, color_by='line_id',
                     point_size=6, line_width=3, label=None):
    """Add a nexus to a PyVista plotter as curves on the torus.

    Parameters:
        G: nexus graph.
        plotter: existing pv.Plotter (creates one if None).
        color_by: 'line_id' — one colour per original LineSubset;
                  'component' — one colour per connected component.
        point_size, line_width: visual parameters.
        label: Label for the legend (applied to the first mesh).

    Returns:
        plotter: PyVista Plotter with the nexus added.
    """
    import pyvista as pv

    if plotter is None:
        plotter = pv.Plotter()

    # Gather edge segments.
    th_map = {n: G.nodes[n]['theta1'] for n in G.nodes}
    b2_map = {n: G.nodes[n]['beta2'] for n in G.nodes}

    # One polyline per connected component, arranged to follow edge paths.
    for comp in nx.connected_components(G):
        sub = G.subgraph(comp)
        if sub.number_of_nodes() == 1:
            # Singleton node → just a sphere.
            n = list(comp)[0]
            p = _torus_point(np.array([th_map[n]]),
                             np.array([np.angle(b2_map[n])]))
            plotter.add_mesh(pv.PolyData(p), color='red',
                             point_size=point_size * 2,
                             render_points_as_spheres=True)
            continue

        # For a simple cycle, extract the cycle path.
        # For a general graph, extract edges individually.
        try:
            cyc = nx.find_cycle(sub)
            # networkx < 3 returns (u, v); >= 3 returns (u, v, key).
            verts_th = [th_map[edge[0]] for edge in cyc]
            verts_b2 = [b2_map[edge[0]] for edge in cyc]
            # close the loop
            verts_th.append(verts_th[0])
            verts_b2.append(verts_b2[0])
            pos = _torus_point(np.array(verts_th),
                              np.unwrap(np.angle(np.array(verts_b2, dtype=complex))))
            n_pts = len(pos)
            lines = np.concatenate([[n_pts], np.arange(n_pts), [0]])
            mesh = pv.PolyData(pos, lines=lines)
        except nx.NetworkXNoCycle:
            # Non-cycle: draw each edge individually.
            pos = []
            lines_parts = []
            offset = 0
            for u, v in sub.edges():
                th1_uv = np.array([th_map[u], th_map[v]])
                th2_uv = np.unwrap(np.angle(np.array(
                    [b2_map[u], b2_map[v]], dtype=complex)))
                p_uv = _torus_point(th1_uv, th2_uv)
                pos.append(p_uv)
                lines_parts.append([2, offset, offset + 1])
                offset += 2
            if pos:
                pos = np.concatenate(pos, axis=0)
                lines = np.concatenate([[offset], np.arange(offset)])
                mesh = pv.PolyData(pos, lines=np.array(lines_parts).ravel())
            else:
                continue

        # Color.
        if color_by == 'line_id':
            # Use the first edge's line_id as representative.
            edge = list(sub.edges(data=True))[0] if sub.number_of_edges() > 0 else None
            lid = edge[2].get('line_id', 0) if edge else 0
            colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
                      "#ff7f00", "#a65628", "#f781bf", "#999999"]
            color = colors[lid % len(colors)]
        else:
            # Color by component index.
            colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"]
            color = colors[len(plotter.renderer.actors) % len(colors)]

        plotter.add_mesh(mesh, color=color, line_width=line_width,
                         label=label)

        # Also add node spheres.
        pts_th = np.array([th_map[n] for n in comp])
        pts_b2 = np.array([np.angle(b2_map[n]) for n in comp])
        pts = _torus_point(pts_th, pts_b2)
        plotter.add_mesh(pv.PolyData(pts), color=color,
                         point_size=point_size,
                         render_points_as_spheres=True)

    return plotter


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def main():
    from trivial_line_cache import load

    data = load()
    E_grid = data['E_grid']
    results = data['results']

    # Pick representative E values.
    targets = [-1.0, -0.5, -0.1, 0.0, 0.1, 0.5, 1.0]

    print("=" * 72)
    print("Step 1: Nexus generation for trivial Hermitian model")
    print("  Model: H = cos(kx + pi/12) + cos(ky)")
    print("  Tolerance: 1e-12 (chordal distance in (beta1, beta2) 4D space)")
    print("=" * 72)

    nexuses = {}
    for E_target in targets:
        # Find the cached result closest to E_target.
        best = None
        for i, r in enumerate(results):
            E = float(E_grid[i])
            if r.success and r.index == (0, 2) and abs(E - E_target) < 0.005:
                best = r
                break

        if best is None:
            print(f"\nE={E_target:+.1f}: no cached result found")
            continue

        lines = [s for s in best.subsets if isinstance(s, LineSubset)]
        G = build_nexus(lines, tol=1e-12)
        info = nexus_info(G)
        nexuses[E_target] = G

        print(f"\n--- E = {E_target:+.1f} ---")
        print(f"  LineSubsets: {len(lines)}")
        print(f"  Raw vertices: {info['n_raw_vertices']}")
        print(f"  Merged nodes: {info['n_nodes']}  (merges: {info['merge_count']})")
        print(f"  Edges: {info['n_edges']}")
        print(f"  Connected components: {info['n_components']}")
        print(f"  Component types: {info['component_types']}")
        print(f"  Cycles: {info['n_cycles']}  lengths: {info['cycle_lengths']}")
        print(f"  Degree histogram: {info['degree_histogram']}")

        # Verify merges are at the endpoints.
        if info['merge_count'] > 0:
            # Find nodes with line_ids containing more than one LineSubset.
            merged_nodes = [n for n, d in G.nodes(data=True)
                           if len(d.get('line_ids', set())) > 1]
            for n in merged_nodes:
                d = G.nodes[n]
                print(f"    Merged node {n}: theta1={d['theta1']:.6f} "
                      f"beta2={d['beta2']:.6f} line_ids={d['line_ids']}")

        for comp_d in info['components']:
            print(f"    comp: {comp_d['n_nodes']} nodes, {comp_d['n_edges']} edges, "
                  f"type={comp_d['type']}, "
                  f"th1=[{comp_d['theta1_range'][0]:.4f}, {comp_d['theta1_range'][1]:.4f}]")

    # Show E=0 case in detail.
    print(f"\n--- E = 0.0  detailed analysis ---")
    G0 = nexuses.get(0.0)
    if G0:
        info0 = nexus_info(G0)
        if info0['n_components'] >= 2:
            print("  At E=0: 2 LineSubsets are each independently closed loops.")
            print("  They are complex conjugates, crossing at beta2=+-1 but sample")
            print("  points do not land on the crossings → 2 disconnected cycles.")
            print("  This is correct behavior with tolerance 1e-12.")
            print("  To glue them, one would need to either:")
            print("    a) increase tolerance to ~0.5 to catch the crossing points,")
            print("    b) add the analytic crossing points as extra vertices,")
            print("    c) use the cut-and-match strategy from demo_topo_mismatch.py.")

    # Torus visualization for a few key E values.
    print("\n" + "=" * 72)
    print("Visualizing nexuses on the torus...")
    print("=" * 72)

    # Matplotlib fallback visualization (works headless).
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    plot_Es = [-1.0, -0.5, 0.0, 0.5, 1.0]
    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3']

    for ax, E_target in zip(axes.flat, plot_Es):
        if E_target not in nexuses:
            continue
        G = nexuses[E_target]
        info = nexus_info(G)

        for comp_idx, comp in enumerate(nx.connected_components(G)):
            sub = G.subgraph(comp)
            # Extract ordered path for the cycle.
            try:
                cyc = nx.find_cycle(sub)
                th = [G.nodes[edge[0]]['theta1'] for edge in cyc]
                b2 = [G.nodes[edge[0]]['beta2'] for edge in cyc]
                # close
                th.append(th[0])
                b2.append(b2[0])
            except nx.NetworkXNoCycle:
                # Fallback: use all nodes sorted by theta1.
                nodes = sorted(comp, key=lambda n: G.nodes[n]['theta1'])
                th = [G.nodes[n]['theta1'] for n in nodes]
                b2 = [G.nodes[n]['beta2'] for n in nodes]

            th2 = np.angle(np.array(b2, dtype=complex))
            # Unwrap for smooth curves.
            th2_uw = np.unwrap(th2)
            color = colors[comp_idx % len(colors)]
            ax.plot(th, th2_uw % (2 * np.pi), '.-', color=color,
                    markersize=2, linewidth=1.5)

            # Mark merged nodes (nodes coming from multiple LineSubsets).
            merged = [n for n in comp if len(G.nodes[n].get('line_ids', set())) > 1]
            if merged:
                th_m = [G.nodes[n]['theta1'] for n in merged]
                b2_m = [np.angle(G.nodes[n]['beta2']) for n in merged]
                ax.plot(th_m, b2_m, 'o', color='black', markersize=6,
                        markerfacecolor='none', markeredgewidth=2)

        ax.set_title(f"E = {E_target:+.1f}  ({info['n_components']} comp, "
                     f"{info['n_nodes']} nodes, {info['merge_count']} merges)")
        ax.set_xlabel('theta1')
        ax.set_ylabel('theta2')
        ax.set_xlim(0, 2 * np.pi)
        ax.set_ylim(0, 2 * np.pi)

    # Hide unused subplot if odd number.
    for ax in axes.flat[len(plot_Es):]:
        ax.set_visible(False)

    fig.suptitle("Nexus graphs in (theta1, theta2) plane\n"
                 "Hollow circles = merged nodes (endpoints where branches meet)",
                 fontsize=13)
    fig.tight_layout()
    out_path = Path(__file__).resolve().parent / "nexus_theta12.png"
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Plot saved to {out_path}")


if __name__ == "__main__":
    main()
