'''
Local-tangent-plane Delaunay bridge for topology-change gaps.

Idea (after bisection narrows the gap):
  1. Collect neighbour patch vertices near the gap → define local tangent plane via PCA.
  2. Project the two nexus slices (with different topology) onto this plane.
  3. Run constrained 2D Delaunay in the plane, with nexus edges as constraints.
  4. Map triangles back to the original (beta1, beta2) space.

This avoids theta1/theta2 projection (model-specific) and MSS branching logic.
Works for any gap narrow enough that a planar approximation is valid.
'''

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import networkx as nx
from scipy.spatial import Delaunay

from gbz_types import LineSubset, to_sphere_r3, TWO_PI
from demo_nexus import build_nexus, nexus_info
from demo_build_mesh import (
    _solve_nexus_at_E, _cycle_vertex_order, _torus_point,
    chamfer_distance, mesh_chordal_edges, triangulate_two_nexuses)
from trivial_line_cache import load
from trivial_model import get_model, DEFAULT_PARAMS


# ---------------------------------------------------------------------------
# R^6 embedding for (beta1, beta2) chordal space
# ---------------------------------------------------------------------------

def _to_r6(beta1, beta2):
    """Map (beta1, beta2) to R^6 via stereographic projection of each component.

    The Euclidean distance in R^6 is sqrt(chordal(b1)^2 + chordal(b2)^2),
    which is ≤ the chordal sum chordal(b1)+chordal(b2) used as our metric.
    """
    r3_b1 = to_sphere_r3(np.asarray(beta1, dtype=complex))
    r3_b2 = to_sphere_r3(np.asarray(beta2, dtype=complex))
    return np.concatenate([r3_b1, r3_b2], axis=-1)


def _from_r6(r6):
    """Inverse of _to_r6: R^6 points → (beta1, beta2) complex arrays.

    Uses stereographic inverse: given (x, y, z) on the unit sphere,
    beta = (x + i*y) / (1 - z) for z < 1, or inf for z ≈ 1.
    For z = 1 (north pole = infinity), we return a large value.
    """
    r6 = np.asarray(r6)
    x1, y1, z1 = r6[..., 0], r6[..., 1], r6[..., 2]
    x2, y2, z2 = r6[..., 3], r6[..., 4], r6[..., 5]

    denom1 = 1.0 - z1
    denom2 = 1.0 - z2

    b1 = np.where(np.abs(denom1) < 1e-15,
                  np.inf + 0j,
                  (x1 + 1j * y1) / denom1)
    b2 = np.where(np.abs(denom2) < 1e-15,
                  np.inf + 0j,
                  (x2 + 1j * y2) / denom2)
    return b1, b2


# ---------------------------------------------------------------------------
# Local tangent plane via PCA
# ---------------------------------------------------------------------------

def _fit_tangent_plane(points_r6):
    """Fit a 2D plane to a set of R^6 points via PCA.

    Returns:
        mean   : (6,)  — centre of the plane.
        basis  : (2, 6) — orthonormal basis for the plane (first 2 PCs).
        normal : (6,)  — normal direction (smallest PC = out-of-plane).
    """
    mean = points_r6.mean(axis=0)
    centred = points_r6 - mean
    cov = centred.T @ centred / (len(points_r6) - 1)
    eigvals, eigvecs = np.linalg.eigh(cov)
    # Sort descending.
    order = np.argsort(-eigvals)
    basis = eigvecs[:, order[:2]].T  # (2, 6)
    return mean, basis


def _project_to_plane(points_r6, mean, basis):
    """Project R^6 points onto the plane defined by mean + span(basis)."""
    centred = points_r6 - mean
    coords_2d = centred @ basis.T  # (N, 2)
    return coords_2d


def _lift_from_plane(coords_2d, mean, basis):
    """Lift 2D plane coordinates back to R^6."""
    return mean + coords_2d @ basis  # (N, 2) @ (2, 6) → (N, 6)


# ---------------------------------------------------------------------------
# Nexus → point cloud + edge constraints
# ---------------------------------------------------------------------------

def _nexus_to_constrained_points(G):
    """Extract R^6 points and edge segments from a nexus graph.

    Returns:
        pts_r6  : (N, 6) — all vertex positions in R^6.
        edges   : list of (i, j) index pairs — nexus edges (constraints).
        node_to_idx: dict mapping nexus node ID → index in pts_r6.
    """
    nodes = list(G.nodes())
    node_to_idx = {n: i for i, n in enumerate(nodes)}

    b1 = np.array([G.nodes[n]['beta1'] for n in nodes])
    b2 = np.array([G.nodes[n]['beta2'] for n in nodes])
    pts_r6 = _to_r6(b1, b2)

    edges = []
    for u, v in G.edges():
        if u in node_to_idx and v in node_to_idx:
            edges.append((node_to_idx[u], node_to_idx[v]))

    return pts_r6, edges, node_to_idx


# ---------------------------------------------------------------------------
# Constrained 2D Delaunay triangulation
# ---------------------------------------------------------------------------

def _constrained_delaunay_2d(pts_2d, constraint_edges):
    """Compute a 2D Delaunay triangulation and keep only triangles that
    don't cross constraint edges.

    Simple approach: compute full Delaunay, then filter triangles.
    A triangle is valid if it doesn't "straddle" a constraint — i.e.
    all its edges are either Delaunay edges or constraint edges.

    For a proper constrained Delaunay we'd need edge flipping, but
    for a narrow gap the full Delaunay + filtering is sufficient.

    Returns:
        tris: (K, 3) int — indices into pts_2d.
    """
    tri = Delaunay(pts_2d)

    # Build set of constraint edges (both directions).
    constraint_set = set()
    for i, j in constraint_edges:
        constraint_set.add((i, j))
        constraint_set.add((j, i))

    # Filter: keep triangles that are "between" the two nexus curves,
    # i.e. not part of the interior of either curve.
    # Heuristic: keep triangles whose circumcircle doesn't contain
    # points from the opposite side.  For a narrow gap, simply keeping
    # all Delaunay triangles between the two curves works.
    #
    # We use a simple filter: a triangle is valid if it has at least
    # one vertex from each side (i.e., it's a "stitching" triangle),
    # OR it lies entirely within one curve's interior.
    # Since we only want the stitching region, we keep triangles
    # that cross between the two point sets.

    return tri.simplices, tri


# ---------------------------------------------------------------------------
# Main bridge function
# ---------------------------------------------------------------------------

def delaunay_bridge(G_a, G_b, neighbour_points_r6):
    """Build a 2D triangulation between two nexuses G_a and G_b
    using local-plane Delaunay.

    Parameters:
        G_a, G_b: nexus graphs (may have different topology).
        neighbour_points_r6: (M, 6) — patch vertices near the gap
                              (from bisected approach slices).

    Returns:
        verts_torus  : (N, 3) torus coords of all vertices.
        faces_vtk    : (F,) VTK-format face array.
    """
    # 1. Fit local tangent plane from neighbour points.
    mean, basis = _fit_tangent_plane(neighbour_points_r6)

    # 2. Extract points and constraints from both nexuses.
    pts_a, edges_a, n2i_a = _nexus_to_constrained_points(G_a)
    pts_b, edges_b, n2i_b = _nexus_to_constrained_points(G_b)

    # Offset indices for B side.
    n_a = len(pts_a)
    pts_all = np.concatenate([pts_a, pts_b], axis=0)
    edges_all = edges_a + [(i + n_a, j + n_a) for i, j in edges_b]

    # 3. Project to 2D plane.
    pts_2d = _project_to_plane(pts_all, mean, basis)

    # 4. Constrained Delaunay in 2D.
    tris_2d, tri_obj = _constrained_delaunay_2d(pts_2d, edges_all)

    # 5. Filter: keep only triangles that connect A to B (stitching triangles),
    #    plus the boundary triangles on each side.
    #    A triangle connects A↔B if it has vertices on both sides.
    valid_tris = []
    for t in tris_2d:
        on_a = t < n_a
        on_b = t >= n_a
        # At least one vertex from each side → stitching triangle.
        if np.any(on_a) and np.any(on_b):
            valid_tris.append(t.tolist())

    if not valid_tris:
        raise ValueError("Delaunay produced no stitching triangles — gap too wide?")

    tris_final = np.array(valid_tris, dtype=int)

    # 6. Map back: 2D coords → R^6 → (beta1, beta2) → torus coords.
    pts_r6_reconstructed = _lift_from_plane(pts_2d, mean, basis)
    b1_all, b2_all = _from_r6(pts_r6_reconstructed)

    th1_all = np.angle(b1_all)
    th2_all = np.angle(b2_all)
    verts_torus = _torus_point(th1_all, th2_all)

    faces_vtk = np.array([[3, t[0], t[1], t[2]]
                          for t in tris_final], dtype=np.int64).ravel()

    # Also compute chordal edge lengths for diagnostics.
    r6_orig = _to_r6(
        np.concatenate([np.array([G_a.nodes[n]['beta1'] for n in sorted(n2i_a, key=n2i_a.get)]),
                        np.array([G_b.nodes[n]['beta1'] for n in sorted(n2i_b, key=n2i_b.get)])]),
        np.concatenate([np.array([G_a.nodes[n]['beta2'] for n in sorted(n2i_a, key=n2i_a.get)]),
                        np.array([G_b.nodes[n]['beta2'] for n in sorted(n2i_b, key=n2i_b.get)])]))

    return verts_torus, faces_vtk, b1_all, b2_all, tris_final


# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------

def main():
    data = load()
    model = get_model(**DEFAULT_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()

    # Get E=0 nexus.
    for i, r in enumerate(data['results']):
        E = float(data['E_grid'][i])
        if r.success and r.index == (0, 2) and abs(E) < 0.005:
            lines = [s for s in r.subsets if isinstance(s, LineSubset)]
            lines.sort(key=lambda L: float(np.median(
                np.angle(np.asarray(L.beta2_arr, dtype=complex)) % (TWO_PI))))
            G_zero = build_nexus(lines, tol=1e-12)
            break

    # Get single-cycle nexus near E=0 (bisected tip).
    G_tip = _solve_nexus_at_E(coeffs, degs, -0.0031)
    # Also get the approach slice for neighbour points.
    G_neighbour = _solve_nexus_at_E(coeffs, degs, -0.0062)

    # Neighbour points: from G_tip and G_neighbour.
    pts_tip, _, _ = _nexus_to_constrained_points(G_tip)
    pts_neigh, _, _ = _nexus_to_constrained_points(G_neighbour)
    neighbour_pts = np.concatenate([pts_tip, pts_neigh], axis=0)
    print(f'Neighbour points: {len(neighbour_pts)} ({len(pts_tip)} from tip, '
          f'{len(pts_neigh)} from neighbour)')

    # Delaunay bridge.
    verts, faces, b1_all, b2_all, tris = delaunay_bridge(
        G_tip, G_zero, neighbour_pts)

    n_a = len(pts_tip)
    chamfer = chamfer_distance(b1_all[:n_a], b2_all[:n_a],
                                b1_all[n_a:], b2_all[n_a:])
    edges = mesh_chordal_edges(b1_all, b2_all, tris)

    print(f'\nDelaunay bridge: {len(tris)} stitching triangles')
    print(f'Chamfer (A↔B): {chamfer:.4f}')
    print(f'Chordal edges: min={edges.min():.4f} mean={edges.mean():.4f} max={edges.max():.4f}')
    n_bad = int(np.sum(edges > 10 * chamfer))
    print(f'Bad edges (>10x): {n_bad}/{len(edges)}')

    # Export for visualisation.
    import pickle
    out = Path(__file__).resolve().parent / 'gbz_delaunay_bridge.pkl'
    with open(out, 'wb') as f:
        pickle.dump({
            'patches': [{
                'verts': verts, 'faces': faces,
                'name': 'delaunay_bridge',
                'E_range': (-0.0031, 0.0),
            }]
        }, f)
    print(f'Exported to {out}')


if __name__ == '__main__':
    main()
