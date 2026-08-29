'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-28
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

"""Torus visualisation of the Morse contour mesh (run WITH a display).

    python playground/morse_torus_view.py               # critical strip debug
    python playground/morse_torus_view.py --mode full   # whole saved mesh
    python playground/morse_torus_view.py --screenshot out.png

The torus embedding maps (theta1, theta2) by
    major angle = theta1 (fixed major radius R)
    minor angle = theta2 (tube radius r)
which is 2pi-periodic in BOTH angles — seam-crossing triangles need no
unwrapping (the flat-chart picture breaks there; the torus does not).

Strip mode shows, on a faint base torus:
  * blue tube  — flank loop at E = -0.05
  * red tube   — flank loop at E = +0.05
  * white thin tubes — the critical fiber at E = 0 (the two crossing
    corridors; their crossings are the saddle critical points)
  * black spheres — inserted critical points
  * orange surface — fan triangles (critical vertex present)
  * teal surface  — remainder-arc FKU strips
so the corridor pairing question ("which blue arc sews to which red
arc") is readable at a glance at each saddle.
"""

import argparse
import pickle
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

TORUS_R = 3.0
TORUS_r = 1.0


def torus_point(theta1, theta2):
    t1 = np.asarray(theta1, dtype=float)
    t2 = np.asarray(theta2, dtype=float)
    x = (TORUS_R + TORUS_r * np.cos(t2)) * np.cos(t1)
    y = (TORUS_R + TORUS_r * np.cos(t2)) * np.sin(t1)
    z = TORUS_r * np.sin(t2)
    return np.stack([x, y, z], axis=-1)


def base_torus(n_major=120, n_minor=60):
    """Faint reference torus surface."""
    import pyvista as pv
    u = np.linspace(0, 2 * np.pi, n_major)
    v = np.linspace(0, 2 * np.pi, n_minor)
    return pv.ParametricTorus() if hasattr(pv, "ParametricTorus") else None


def loop_tube(polyline_2d, radius=0.012):
    """Thin tube following a loop on the torus."""
    import pyvista as pv
    pts = torus_point(polyline_2d[:, 0], polyline_2d[:, 1])
    spline = pv.Spline(pts, int(max(1000, 4 * len(pts))))
    return spline.tube(radius=radius)


def faces_polydata(points_on_torus, tris):
    import pyvista as pv
    if not len(tris):
        return None
    faces = np.column_stack([np.full(len(tris), 3),
                             np.asarray(tris, dtype=int)]).ravel()
    return pv.PolyData(np.asarray(points_on_torus), faces)


def build_strip():
    from trivial_model import get_model, DEFAULT_PARAMS
    import demo_morse_mesh as dmm

    model = get_model(**DEFAULT_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()
    oracle = dmm.make_oracle(coeffs, degs)
    lo, hi, sc = oracle(-0.05), oracle(0.05), oracle(0.0)

    builder = dmm.MeshBuilder()
    builder.register_slice(lo)
    builder.register_slice(hi)
    n_crit = dmm.glue_critical(lo, hi, sc, builder)
    return lo, hi, sc, builder, n_crit


def builder_torus_coords(builder):
    pts = []
    for E, b1, b2 in builder.coords:
        pts.append(torus_point(np.angle(b1) % (2 * np.pi),
                               np.angle(b2) % (2 * np.pi)))
    return np.asarray(pts)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["strip", "full"], default="strip")
    ap.add_argument("--flank-E", type=float, default=0.05,
                    help="flank distance for strip mode")
    ap.add_argument("--screenshot", default=None)
    args = ap.parse_args()

    import pyvista as pv

    pl = pv.Plotter(window_size=(1500, 1100))

    if args.mode == "strip":
        lo, hi, sc, builder, n_crit = build_strip()
        crit_gids = {g for g, _, _ in builder.crit_info}
        print(f"strip: {len(builder.tris)} triangles, "
              f"{n_crit} critical points")

        # flank loops as tubes (through unwrapped polylines for smoothness)
        import demo_morse_mesh as dmm
        pl.add_mesh(loop_tube(dmm._unwrapped_polyline(lo.loops[0])),
                    color="tab:blue", label="E = -|flank|")
        pl.add_mesh(loop_tube(dmm._unwrapped_polyline(hi.loops[0])),
                    color="tab:red", label="E = +|flank|")
        # critical fiber corridors (thin, white-ish)
        for k, loop in enumerate(sc.loops):
            pl.add_mesh(loop_tube(dmm._unwrapped_polyline(loop),
                                  radius=0.005),
                        color="white", opacity=0.8,
                        label="E = 0 corridors" if k == 0 else None)

        pts = builder_torus_coords(builder)
        fan_tris = [t for t in builder.tris
                    if len(set(t)) == 3 and set(t) & crit_gids]
        arc_tris = [t for t in builder.tris
                    if len(set(t)) == 3 and not (set(t) & crit_gids)]
        m_fan = faces_polydata(pts, fan_tris)
        m_arc = faces_polydata(pts, arc_tris)
        if m_fan is not None:
            pl.add_mesh(m_fan, color="orange", opacity=0.55,
                        show_edges=True, edge_color="darkorange",
                        label=f"fan tris ({len(fan_tris)})")
        if m_arc is not None:
            pl.add_mesh(m_arc, color="teal", opacity=0.4,
                        label=f"arc strips ({len(arc_tris)})")
        for g in crit_gids:
            pl.add_mesh(pv.Sphere(0.05, center=pts[g]), color="black")
        pl.add_legend(bcolor="white")

    else:
        src = Path(__file__).resolve().parent / "gbz_morse_mesh.pkl"
        with open(src, "rb") as f:
            m = pickle.load(f)
        E = m["verts"][:, 0]
        th1 = np.angle(m["verts"][:, 2] + 1j * m["verts"][:, 3])
        th2 = np.angle(m["verts"][:, 4] + 1j * m["verts"][:, 5])
        pts = torus_point(th1, th2)
        mesh = faces_polydata(pts, m["tris"])
        pl.add_mesh(mesh, scalars=E, cmap="coolwarm",
                    scalar_bar_args={"title": "E"})
        # pl.add_mesh(mesh, color="blue")
        # crit = np.array([g for g, _, _ in m["crit_info"]], dtype=int)
        # if len(crit):
        #     pl.add_mesh(pv.PolyData(pts[crit]), color="black",
        #                 point_size=15, render_points_as_spheres=True)
        print(f"full mesh: {len(m['tris'])} triangles")

    if args.screenshot:
        pl.screenshot(args.screenshot)
        print("saved ->", args.screenshot)
    else:
        pl.show()


if __name__ == "__main__":
    main()
