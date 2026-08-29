"""Measure the angle between hole-fill mesh lines and the original curves.

Suspicion under test: fills add edges that run ALONG the boundary curve
(parallel angle ~0 deg) instead of spanning the hole transversally.

Metric: for an edge (u,v), parallel angle = angle between (p_v - p_u) and
the curve tangent at u (and at v), folded so 0 deg = parallel, 90 deg =
perpendicular.  Tangent at a loop vertex = central difference of its two
loop neighbours, in the same chordal r3 space the mesh builder uses.
"""
import numpy as np

import demo_morse_mesh as D


def replay_until_fill():
    """Copy of demo_morse_mesh.main() up to (not including) fill_holes."""
    from trivial_line_cache import load
    from trivial_model import get_model, DEFAULT_PARAMS
    from pygbz2d.core import LineSubset as LS

    data = load()
    E_grid, results = data["E_grid"], data["results"]

    all_slices = []
    for i, r in enumerate(results):
        lines = [s for s in r.subsets if isinstance(s, LS)]
        if lines:
            all_slices.append(D.nexus_to_slice(
                float(E_grid[i]),
                D.build_nexus(D._sorted_lines(lines), tol=1e-12)))

    from collections import Counter
    far_count = Counter(len(s.loops) for s in all_slices).most_common(1)[0][0]
    slices = [s for s in all_slices
              if not D.is_critical_slice(s, far_count=far_count)]

    for k in range(len(slices) - 1):
        sa, sb = slices[k], slices[k + 1]
        for i, j in D.match_components(sa.loops, sb.loops):
            sb.loops[j] = D.orient_pair(sa.loops[i], sb.loops[j])

    model = get_model(**DEFAULT_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()
    oracle = D.make_oracle(coeffs, degs)

    builder = D.MeshBuilder()
    for sl in slices:
        builder.register_slice(sl)

    for k in range(len(slices) - 1):
        D.mesh_interval(slices[k], slices[k + 1], oracle, builder,
                        log=[], far_count=far_count)
    D.grow_cap(slices[0], -1, oracle, builder)
    D.grow_cap(slices[-1], +1, oracle, builder)
    return slices, builder


def build_tangents(slices, builder):
    """gid -> normalized curve tangent (chordal r3), from loop neighbours."""
    r3 = {g: D._r3_of_coord(c) for g, c in enumerate(builder.coords)}
    tan = {}
    adj = set()          # loop-adjacent gid pairs (unordered)
    for sl in slices:
        for j, loop in enumerate(sl.loops):
            gids = [builder.vid[(sl.E, j, int(loop.perm[p]))]
                    for p in range(loop.n)]
            for p in range(loop.n):
                u, v = gids[p], gids[(p + 1) % loop.n]
                adj.add((min(u, v), max(u, v)))
                if loop.n >= 3:
                    a, c = gids[p - 1], gids[(p + 1) % loop.n]
                    t = r3[c] - r3[a]
                    nn = np.linalg.norm(t)
                    if nn > 1e-15:
                        tan[u] = t / nn
    return r3, tan, adj


def par_angle(vec, t):
    if t is None:
        return None
    nv = np.linalg.norm(vec)
    if nv < 1e-15:
        return None
    c = abs(vec @ t) / nv
    return float(np.degrees(np.arccos(np.clip(c, 0.0, 1.0))))


def hist(angles, label):
    a = [x for x in angles if x is not None]
    if not a:
        print(f"  {label}: (no edges)")
        return
    buckets = [0, 2, 5, 10, 20, 45, 90.0001]
    counts = [sum(1 for x in a if buckets[i] <= x < buckets[i + 1])
              for i in range(len(buckets) - 1)]
    total = len(a)
    q = np.percentile(a, [5, 25, 50, 75, 95])
    bars = " ".join(f"{c:>5}" for c in counts)
    print(f"  {label}: n={total}  p5/25/50/75/95 = "
          + "/".join(f"{x:5.1f}" for x in q) + " deg")
    print(f"      [0-2  2-5  5-10 10-20 20-45 45-90]  {bars}")


def main():
    # record EVERY slice the pipeline registers (grid slices AND the
    # intermediate ladder/growth probe slices) so tangents exist for
    # growth-region vertices too
    seen_slices = []
    _orig_reg = D.MeshBuilder.register_slice

    def _reg(self, sl):
        seen_slices.append(sl)
        _orig_reg(self, sl)
    D.MeshBuilder.register_slice = _reg

    slices, builder = replay_until_fill()
    r3, tan, loop_adj = build_tangents(seen_slices, builder)

    n0 = len(builder.tris)
    print(f"\npre-fill: {n0} triangles")
    D.fill_holes(builder, width_factor=1.3)
    fill_tris = builder.tris[n0:]
    print(f"fill added {len(fill_tris)} triangles")

    old_edges = set()
    for t in builder.tris[:n0]:
        if len(set(t)) < 3:
            continue
        for x, y in ((t[0], t[1]), (t[1], t[2]), (t[0], t[2])):
            old_edges.add((min(x, y), max(x, y)))

    def edge_info(u, v):
        vec = r3[v] - r3[u]
        L = float(np.linalg.norm(vec))
        au, av = par_angle(vec, tan.get(u)), par_angle(vec, tan.get(v))
        return L, au, av

    # ---- per-fill-triangle edge report -------------------------------
    print("\n=== fill triangles ===")
    n_new = 0
    seen = set()
    for ti, tri in enumerate(fill_tris):
        E_t = [builder.coords[g][0] for g in tri]
        deg = " DEGENERATE(repeated gid)" if len(set(tri)) < 3 else ""
        print(f"  tri[{ti}] gids {tri} E {E_t}{deg}")
    print("\n--- NEW edges (not in old mesh) ---")
    for tri in fill_tris:
        for x, y in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[0], tri[2])):
            e = (min(x, y), max(x, y))
            if e in old_edges or e in seen:
                continue
            seen.add(e)
            L, au, av = edge_info(*e)
            n_new += 1
            am = min(a for a in (au, av) if a is not None) \
                if (au is not None or av is not None) else None
            flag = "  <-- PARALLEL" if (am is not None and am < 10) else ""
            print(f"  edge {e}: len {L:.4f}  ang {au}/{av} "
                  f"(min {am:.1f}){flag}" if am is not None
                  else f"  edge {e}: len {L:.4f}  ang None")

    # ---- aggregate comparison ----------------------------------------
    def collect(edges):
        au_l, am_l = [], []
        for u, v in edges:
            vec = r3[v] - r3[u]
            au, av = par_angle(vec, tan.get(u)), par_angle(vec, tan.get(v))
            cands = [a for a in (au, av) if a is not None]
            if cands:
                au_l.extend(cands)
                am_l.append(min(cands))
        return au_l, am_l

    print(f"\n=== angle-to-curve distributions ({n_new} new fill edges) ===")
    print("groups: A=same-loop (must be parallel; metric sanity), "
          "B=old cross-curve (healthy baseline), "
          "C=fill inherited boundary edges, D=fill NEW edges")
    aA, _ = collect([e for e in old_edges if e in loop_adj])
    aB, _ = collect([e for e in old_edges if e not in loop_adj])
    aC, _ = collect([e for e in seen if e in old_edges])  # empty by def
    aD, _ = collect(seen)
    for lab, arr in (("A same-loop", aA), ("B old cross", aB),
                     ("C inherited", aC), ("D fill NEW", aD)):
        hist(arr, lab)

    # ---- non-manifold audit on the final mesh ------------------------
    from collections import Counter
    inc = Counter()
    for t in builder.tris:
        if len(set(t)) < 3:
            continue
        for x, y in ((t[0], t[1]), (t[1], t[2]), (t[0], t[2])):
            inc[(min(x, y), max(x, y))] += 1
    fill_edge_set = set()
    for tri in fill_tris:
        for x, y in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[0], tri[2])):
            fill_edge_set.add((min(x, y), max(x, y)))
    nm = [(e, c) for e, c in inc.items() if c > 2]
    print(f"\nnon-manifold edges (>2 incident): {len(nm)}")
    for e, c in nm[:10]:
        in_fill = e in fill_edge_set
        L, au, av = edge_info(*e)
        print(f"  {e}: {c} tris, in-fill={in_fill}, len {L:.4f}, "
              f"ang {au}/{av}")


if __name__ == "__main__":
    main()
