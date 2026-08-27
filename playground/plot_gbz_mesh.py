'''
Standalone PyVista script to visualise the GBZ mesh saved by demo_build_mesh.py.

Usage:
    python playground/plot_gbz_mesh.py
'''

import pickle
from pathlib import Path
import numpy as np

from pygbz2d.core import TWO_PI
import pyvista as pv

HERE = Path(__file__).resolve().parent
DATA = HERE / "gbz_mesh_full.pkl"

TORUS_R = 3.0
TORUS_r = 1.0


def _torus_point(theta1, theta2):
    t1 = np.asarray(theta1, dtype=float) % (TWO_PI)
    t2 = np.asarray(theta2, dtype=float) % (TWO_PI)
    x = (TORUS_R + TORUS_r * np.cos(t2)) * np.cos(t1)
    y = (TORUS_R + TORUS_r * np.cos(t2)) * np.sin(t1)
    z = TORUS_r * np.sin(t2)
    return np.stack([x, y, z], axis=-1)


def plot_mesh(data_path=None):
    if data_path is None:
        data_path = DATA

    with open(data_path, "rb") as f:
        data = pickle.load(f)

    patches = data["patches"]

    plotter = pv.Plotter()
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
              "#ff7f00", "#a65628", "#f781bf", "#999999"]

    for i, p in enumerate(patches):
        verts = np.asarray(p["verts"], dtype=float)
        faces = np.asarray(p["faces"], dtype=np.int64)
        mesh = pv.PolyData(verts, faces)
        color = colors[i % len(colors)]
        name = p.get("name", f"patch_{i}")
        if "MSS" in name:
            continue
        plotter.add_mesh(mesh, color=color, opacity=1, show_edges=True,
                         label=name)

    # # Torus wireframe for reference (major rings + minor rings).
    # n_ring, n_tube = 40, 12
    # for i in range(n_ring):
    #     a = 2 * np.pi * i / n_ring
    #     pts = _torus_point(np.full(n_tube + 1, a),
    #                        np.linspace(0, 2 * np.pi, n_tube + 1))
    #     n = len(pts)
    #     conn = np.concatenate([[n], np.arange(n)])
    #     plotter.add_mesh(pv.PolyData(pts, lines=conn),
    #                      color="gray", line_width=0.5, opacity=0.3)
    # for j in range(n_tube):
    #     b = 2 * np.pi * j / n_tube
    #     pts = _torus_point(np.linspace(0, 2 * np.pi, n_ring + 1),
    #                        np.full(n_ring + 1, b))
    #     n = len(pts)
    #     conn = np.concatenate([[n], np.arange(n)])
    #     plotter.add_mesh(pv.PolyData(pts, lines=conn),
    #                      color="gray", line_width=0.5, opacity=0.3)

    plotter.add_axes()
    plotter.add_legend(bcolor="white")
    plotter.show(title="GBZ mesh on torus (E<0 red, E>0 blue)")


if __name__ == "__main__":
    plot_mesh(HERE / "gbz_mr_split_test.pkl")
