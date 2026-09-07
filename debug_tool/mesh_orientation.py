'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-09-06
Copyright © Department of Physics, Tsinghua University. All rights reserved

Orientation checker for (theta1, theta2) torus triangle meshes.

The mesh builders (playground periodic Delaunay, demo_pipeline refinement)
inherit Qhull's winding convention and mod-2pi mapping preserves torus
orientation, so a CLEAN mesh should be uniformly oriented — but nothing in
the construction ENFORCES it, and non-manifold pathologies (e.g. the
seam keep-rule over-stuffed fans) destroy orientation silently.  This
tool measures, it does not repair.

Two independent criteria:

1. **Global winding census** (signed areas).  Each triangle is unwrapped
   into the universal cover (shift every vertex to within pi of vertex 0
   per coordinate — faithful iff every torus edge is < pi) and its signed
   area sign is counted.  A uniformly oriented mesh has all nonzero signs
   equal.

2. **Neighbour consistency** (directed-edge balance).  In a consistently
   oriented 2-complex every INTERIOR edge is used once in each direction
   ((u,v) and (v,u) equally often).  Boundary edges (multiplicity 1) are
   asymmetric by design; edges with multiplicity > 2 are non-manifold and
   are flagged separately (orientation is moot there).

CLI: ``python debug_tool/mesh_orientation.py mesh.npz [more.npz ...]``
(each file needs ``verts`` (n, 2) and ``triangles`` (m, 3) arrays, as
saved by playground/demo_pipeline.py).  Exit code 1 if any mesh fails.
'''

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# bootstrap for direct-script execution (python debug_tool/mesh_orientation.py);
# harmless no-op when imported as part of the debug_tool package
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pygbz2d.core import TWO_PI

#: Stored violation samples per category (reports stay printable).
MAX_SAMPLES = 20


@dataclass
class MeshOrientationReport:
    """Orientation diagnostics of one (theta1, theta2) torus mesh."""
    n_vertices: int
    n_triangles: int
    signed_pos: int
    signed_neg: int
    signed_zero: int
    max_edge: float                 # torus metric, rad
    unwrap_valid: bool              # max_edge < pi (signed census faithful)
    n_boundary_edges: int           # undirected multiplicity 1
    n_interior_edges: int           # undirected multiplicity 2
    n_nonmanifold_edges: int        # undirected multiplicity > 2
    balance_violations: list = field(default_factory=list)
    nonmanifold_edges: list = field(default_factory=list)

    @property
    def orientation_uniform(self) -> bool:
        """All nonzero signed areas share one sign (census meaningful only
        when ``unwrap_valid``)."""
        return self.unwrap_valid and (self.signed_neg == 0
                                      or self.signed_pos == 0)

    @property
    def neighbors_consistent(self) -> bool:
        return not self.balance_violations

    @property
    def is_ok(self) -> bool:
        return (self.orientation_uniform and self.neighbors_consistent
                and self.n_nonmanifold_edges == 0 and self.signed_zero == 0)

    def describe(self) -> str:
        lines = [
            f"mesh: V={self.n_vertices} F={self.n_triangles} "
            f"max_edge={self.max_edge:.4f} (unwrap_valid={self.unwrap_valid})",
            f"signed areas: +{self.signed_pos} -{self.signed_neg} "
            f"zero={self.signed_zero} -> "
            + ("uniform" if self.orientation_uniform else "MIXED"),
            f"edges: interior={self.n_interior_edges} "
            f"boundary={self.n_boundary_edges} "
            f"non-manifold={self.n_nonmanifold_edges}",
        ]
        if self.balance_violations:
            lines.append(f"balance violations ({len(self.balance_violations)}"
                         f" shown, edge (u,v): cnt_uv vs cnt_vu): "
                         + "; ".join(f"({u},{v}): {a} vs {b}"
                                     for u, v, a, b
                                     in self.balance_violations))
        if self.nonmanifold_edges:
            lines.append(f"non-manifold edges: "
                         + "; ".join(f"({u},{v}) x{m}" for u, v, m
                                     in self.nonmanifold_edges))
        lines.append("verdict: " + ("OK" if self.is_ok else "INCONSISTENT"))
        return "\n".join(lines)


def _torus_dist(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    d = np.abs(a - b) % TWO_PI
    d = np.minimum(d, TWO_PI - d)
    return np.hypot(d[:, 0], d[:, 1])


def check_mesh_orientation(verts: np.ndarray,
                           triangles: np.ndarray) -> MeshOrientationReport:
    """Full orientation diagnostic of a torus triangle mesh.

    Parameters:
        verts: (n, 2) vertex angles in [0, 2pi).
        triangles: (m, 3) vertex indices.
    """
    verts = np.asarray(verts, dtype=float)
    triangles = np.asarray(triangles, dtype=int)

    # -- 1. signed-area census (unwrap to vertex 0 per coordinate) --------
    a = verts[triangles[:, 0]]
    b = verts[triangles[:, 1]] - np.round(
        (verts[triangles[:, 1]] - a) / TWO_PI) * TWO_PI
    c = verts[triangles[:, 2]] - np.round(
        (verts[triangles[:, 2]] - a) / TWO_PI) * TWO_PI
    cross = ((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
             - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0]))

    dir_edges = np.vstack([triangles[:, [0, 1]], triangles[:, [1, 2]],
                           triangles[:, [2, 0]]])
    max_edge = float(_torus_dist(verts[dir_edges[:, 0]],
                                 verts[dir_edges[:, 1]]).max())
    unwrap_valid = bool(max_edge < np.pi)

    # -- 2. directed-edge balance -----------------------------------------
    dir_keys, dir_cnt = np.unique(dir_edges, axis=0, return_counts=True)
    dir_map = {(int(u), int(v)): int(n)
               for (u, v), n in zip(dir_keys, dir_cnt)}
    und_keys, und_cnt = np.unique(np.sort(dir_edges, axis=1), axis=0,
                                  return_counts=True)

    balance_violations = []
    nonmanifold = []
    n_boundary = n_interior = 0
    for (u, v), m in zip(und_keys, und_cnt):
        u, v, m = int(u), int(v), int(m)
        if m == 1:
            n_boundary += 1
        elif m == 2:
            n_interior += 1
            cuv = dir_map.get((u, v), 0)
            cvu = dir_map.get((v, u), 0)
            if cuv != cvu and len(balance_violations) < MAX_SAMPLES:
                balance_violations.append((u, v, cuv, cvu))
        else:
            nonmanifold.append((u, v, m))

    return MeshOrientationReport(
        n_vertices=int(triangles.max()) + 1,
        n_triangles=int(len(triangles)),
        signed_pos=int((cross > 0).sum()),
        signed_neg=int((cross < 0).sum()),
        signed_zero=int((cross == 0).sum()),
        max_edge=max_edge,
        unwrap_valid=unwrap_valid,
        n_boundary_edges=n_boundary,
        n_interior_edges=n_interior,
        n_nonmanifold_edges=len(nonmanifold),
        balance_violations=balance_violations,
        nonmanifold_edges=nonmanifold[:MAX_SAMPLES],
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description="torus mesh orientation check")
    p.add_argument("files", nargs="+", help="mesh .npz files "
                   "(arrays 'verts' (n,2) and 'triangles' (m,3))")
    args = p.parse_args(argv)

    all_ok = True
    for f in args.files:
        z = np.load(f)
        tri_key = "triangles" if "triangles" in z.files else "tri"
        rep = check_mesh_orientation(z["verts"], z[tri_key])
        print(f"== {f}")
        print(rep.describe())
        print()
        all_ok &= rep.is_ok
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(_main())
