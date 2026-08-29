'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-29
Copyright © Department of Physics, Tsinghua University. All rights reserved

Torus view of the critical-slab prototype (demo_critical_slab.py).

    python playground/demo_critical_slab_view.py
    python playground/demo_critical_slab_view.py --screenshot out.png
'''

import argparse
import pickle
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pickle",
                    default=str(Path(__file__).resolve().parent /
                                "gbz_critical_slab.pkl"))
    ap.add_argument("--screenshot", default=None)
    args = ap.parse_args()

    import pyvista as pv

    with open(args.pickle, "rb") as f:
        m = pickle.load(f)

    verts = m["verts"]
    th1 = np.angle(verts[:, 2] + 1j * verts[:, 3])
    th2 = np.angle(verts[:, 4] + 1j * verts[:, 5])
    pts = torus_point(th1, th2)

    # split fan vs FKU strip triangles
    crit_gids = {g for g, _, _ in m["crit_info"]}
    fan = [t for t in m["tris"] if set(t) & crit_gids]
    strip = [t for t in m["tris"] if not (set(t) & crit_gids)]
    print(f"fan tris: {len(fan)}, FKU strip tris: {len(strip)}")

    def polydata(tris):
        faces = np.column_stack([np.full(len(tris), 3),
                                 np.asarray(tris, dtype=int)]).ravel()
        return pv.PolyData(pts, faces)

    pl = pv.Plotter(window_size=(1500, 1100))
    if fan:
        pl.add_mesh(polydata(fan), color="orange", opacity=1,
                    show_edges=True, edge_color="darkorange",
                    label=f"fan ({len(fan)})")
    if strip:
        pl.add_mesh(polydata(strip), color="teal", opacity=1,
                    label=f"FKU strips ({len(strip)})")

    # for g in crit_gids:
    #     pl.add_mesh(pv.Sphere(0.06, center=pts[g]), color="black")

    pl.add_legend(bcolor="white")

    if args.screenshot:
        pl.screenshot(args.screenshot)
        print("saved ->", args.screenshot)
    else:
        pl.show()


if __name__ == "__main__":
    main()
