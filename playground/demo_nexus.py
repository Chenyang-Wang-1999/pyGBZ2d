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

import numpy as np
import networkx as nx
from scipy.spatial import KDTree

from pygbz2d.core import LineSubset, to_sphere_r3


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


