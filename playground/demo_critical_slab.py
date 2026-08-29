'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-29
Copyright © Department of Physics, Tsinghua University. All rights reserved

Prototype: critical-slab patch for the Hermitian (trivial) model.

Idea (general, not model-specific)
-----------------------------------
Near a Morse saddle P of the E-projection, the GBZ surface is a smooth
graph over the torus chart (theta1, theta2).  The fibre loops of the two
flanking slices (E = -delta and E = +delta) form, in a small chart
neighbourhood of P, TWO U-shaped passages each (the two branches of the
level-set hyperbola).  The patch {|E| <= delta} near P is star-shaped
with respect to P, so it is triangulated exactly by the fan from P over
the boundary points of those four passages sorted by angular order.

Between two saddles the remaining loop arcs (the complement of the four
passages on each flank loop) are topologically consistent cylinders, so
they are stitched with the existing FKU routine.

This replaces the broken ``glue_critical`` fan for the E=0 slice of the
trivial model (which only found ONE passage per flank per crossing and
produced euler=109 / nonmanifold edges).
'''

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import demo_morse_mesh as dmm
from trivial_model import get_model, DEFAULT_PARAMS

TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# Passage detection in the torus chart
# ---------------------------------------------------------------------------

def angular_dist(a, b):
    return np.abs(np.angle(np.exp(1j * (np.asarray(a) - np.asarray(b)))))


def cyclic_runs(mask):
    """Contiguous True runs on a cyclic boolean mask.

    Returns a list of integer arrays, each in cyclic increasing order.
    """
    mask = np.asarray(mask, dtype=bool)
    n = len(mask)
    if not mask.any():
        return []
    if mask.all():
        return [np.arange(n)]
    start = int(np.where(~mask)[0][0])
    m2 = np.concatenate([mask[start:], mask[:start]])
    runs = []
    i = 0
    while i < n:
        if m2[i]:
            j = i
            while j < n and m2[j]:
                j += 1
            runs.append((np.arange(i, j) + start) % n)
            i = j
        else:
            i += 1
    return runs


def lift_run(loop, run, P):
    """Lift a run's (theta1, theta2) to the universal cover near P."""
    t1 = np.asarray(loop.th1)[run].copy()
    t2 = np.asarray(loop.th2)[run].copy()
    t1u = np.unwrap(t1)
    t2u = np.unwrap(t2)
    k1 = int(round((P[0] - t1u.mean()) / TWO_PI))
    k2 = int(round((P[1] - t2u.mean()) / TWO_PI))
    uv = np.column_stack([t1u + k1 * TWO_PI, t2u + k2 * TWO_PI])
    return uv


def passages_on_loop(loop, P, R):
    """All passages of a closed loop within chart-distance R of P.

    Each passage is a dict:
      idx   : loop indices in cyclic increasing order
      uv    : lifted (theta1, theta2) coordinates near P
      s, e  : first / last loop index of the passage (cyclic order)
    """
    th1 = np.asarray(loop.th1) % TWO_PI
    th2 = np.asarray(loop.th2) % TWO_PI
    d = np.sqrt(angular_dist(th1, P[0]) ** 2
                + angular_dist(th2, P[1]) ** 2)
    runs = cyclic_runs(d < R)
    out = []
    for run in runs:
        uv = lift_run(loop, run, P)
        # keep the run in cyclic increasing order; s,e are its endpoints
        out.append({
            "idx": run,
            "uv": uv,
            "s": int(run[0]),
            "e": int(run[-1]),
            "dmin": float(d[run].min()),
        })
    return out


# ---------------------------------------------------------------------------
# Cyclic arcs
# ---------------------------------------------------------------------------

def cyclic_range(n, a, b):
    """Loop indices a -> b inclusive, in increasing cyclic order."""
    if b >= a:
        return np.arange(a, b + 1)
    return np.concatenate([np.arange(a, n), np.arange(0, b + 1)])


def complement_arc_from_endpoint(n, passages, L):
    """Complement arc adjacent to passage endpoint L, oriented away from
    the passage (first index = L)."""
    npass = len(passages)
    for pi, p in enumerate(passages):
        if L == p["s"]:
            prev = passages[(pi - 1) % npass]
            arc = cyclic_range(n, prev["e"], p["s"])
            return arc[::-1]           # starts at L = p['s'], ends at prev['e']
        if L == p["e"]:
            nxt = passages[(pi + 1) % npass]
            return cyclic_range(n, p["e"], nxt["s"])
    raise ValueError(f"index {L} is not a passage endpoint")


def find_complement_arc(n, passages, a, b):
    """Complement arc with endpoints {a,b}, oriented a -> b."""
    npass = len(passages)
    for pi, p in enumerate(passages):
        nxt = passages[(pi + 1) % npass]
        arc = cyclic_range(n, p["e"], nxt["s"])
        if {int(arc[0]), int(arc[-1])} == {int(a), int(b)}:
            if int(arc[0]) == int(a):
                return arc
            return arc[::-1]
    raise ValueError(f"no complement arc with endpoints {a}, {b}")


# ---------------------------------------------------------------------------
# Fan patch at one saddle
# ---------------------------------------------------------------------------

def build_saddle_fan(lo_loop, hi_loop, P, R, builder, key_lo, key_hi):
    """Fan critical point P over the four boundary passages in angular order.

    Returns (rungs, n_tris): rungs is a list of ((side, idx), (side, idx))
    with side 0=lo, 1=hi; n_tris the number of fan triangles added.
    """
    lo_pass = passages_on_loop(lo_loop, P, R)
    hi_pass = passages_on_loop(hi_loop, P, R)
    # keep the two closest passages if more are found
    lo_pass = sorted(lo_pass, key=lambda p: p["dmin"])[:2]
    hi_pass = sorted(hi_pass, key=lambda p: p["dmin"])[:2]
    assert len(lo_pass) == 2 and len(hi_pass) == 2, \
        (len(lo_pass), len(hi_pass))

    # critical vertex at the crossing
    mu1 = float(np.log(abs(lo_loop.b1[0])))   # same mu1 for the whole loop
    b1 = np.exp(mu1 + 1j * P[0])
    b2 = np.exp(1j * P[1])
    P_gid = builder.add_vertex(0.5 * (key_lo[0] + key_hi[0]), b1, b2,
                               note=f"crit E~0")

    # all passage points, each tagged by (side, loop_idx, uv)
    pts = []
    for side, passages, key in ((0, lo_pass, key_lo), (1, hi_pass, key_hi)):
        for p in passages:
            for k in range(len(p["idx"])):
                idx = int(p["idx"][k])
                gid = builder.gid(key[0], key[1], int(idx))
                pts.append((side, idx, gid, p["uv"][k]))

    uv = np.array([q[3] for q in pts])
    ang = np.arctan2(uv[:, 1] - P[1], uv[:, 0] - P[0])
    order = np.argsort(ang)

    rungs = []
    n_tris = 0
    npts = len(order)
    for k in range(npts):
        a = pts[order[k]]
        b = pts[order[(k + 1) % npts]]
        # triangle P - a - b
        tri = (P_gid, a[2], b[2])
        if len({tri[0], tri[1], tri[2]}) == 3:
            builder.tris.append(tri)
            n_tris += 1
        # the edge a-b is a rung when a and b belong to different loops
        if a[0] != b[0]:
            rungs.append(((a[0], a[1]), (b[0], b[1])))
    return rungs, n_tris


# ---------------------------------------------------------------------------
# Critical-point estimation WITHOUT sampling the critical slice
# ---------------------------------------------------------------------------

def estimate_saddles_from_flank(loop, k=8, cyc_min=20, d_max=1.0):
    """Estimate saddle chart positions from ONE clean flank loop.

    Near a saddle the flank loop has two U-shaped passages (the two
    hyperbola branches).  Their tips are the two points of closest
    torus self-approach; the torus midpoint of the two tips is the
    saddle point P.  This works when the critical slice itself was
    never sampled.
    """
    th1 = np.asarray(loop.th1) % TWO_PI
    th2 = np.asarray(loop.th2) % TWO_PI
    n = len(th1)

    # torus self-distance
    D = np.zeros((n, n))
    for i in range(n):
        d1 = angular_dist(th1, th1[i])
        d2 = angular_dist(th2, th2[i])
        D[i] = np.sqrt(d1 ** 2 + d2 ** 2)
    for i in range(n):
        for j in range(n):
            cyc = min((i - j) % n, (j - i) % n)
            if cyc <= k:
                D[i, j] = 1e9

    # local minima with large cyclic separation (= different branches)
    cand = []
    for i in range(n):
        for j in range(n):
            v = D[i, j]
            if v >= d_max:
                continue
            cyc = min((i - j) % n, (j - i) % n)
            if cyc < cyc_min:
                continue
            if v <= D[max(0, i - 4):i + 5, max(0, j - 4):j + 5].min():
                cand.append((i, j, v))

    # deduplicate symmetric (i,j)/(j,i) pairs
    uniq = []
    for i, j, v in sorted(cand, key=lambda x: x[2]):
        if any(abs(u[0] - i) <= 4 and abs(u[1] - j) <= 4 for u in uniq):
            continue
        if any(abs(u[0] - j) <= 4 and abs(u[1] - i) <= 4 for u in uniq):
            continue
        uniq.append((i, j, v))

    # each U-tip = torus midpoint of its two arms (i, j)
    tips = []
    for i, j, v in uniq:
        p1 = float(np.angle(np.exp(1j * th1[i]) + np.exp(1j * th1[j])))
        p2 = float(np.angle(np.exp(1j * th2[i]) + np.exp(1j * th2[j])))
        tips.append((p1 % TWO_PI, p2 % TWO_PI))

    # cluster the two U-tips of each saddle by chart proximity
    saddles = []
    used = set()
    for a in range(len(tips)):
        if a in used:
            continue
        cluster = [a]
        used.add(a)
        for b in range(a + 1, len(tips)):
            if b in used:
                continue
            d = np.sqrt(angular_dist(tips[a][0], tips[b][0]) ** 2
                        + angular_dist(tips[a][1], tips[b][1]) ** 2)
            if d < 1.0:
                cluster.append(b)
                used.add(b)
        x = np.sum([np.exp(1j * tips[p][0]) for p in cluster])
        y = np.sum([np.exp(1j * tips[p][1]) for p in cluster])
        saddles.append((float(np.angle(x) % TWO_PI),
                        float(np.angle(y) % TWO_PI)))
    return saddles


# ---------------------------------------------------------------------------
# Main prototype
# ---------------------------------------------------------------------------

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta", type=float, default=0.05)
    ap.add_argument("--R", type=float, default=None,
                    help="passage window radius in the chart; "
                         "default = 3 * (min passage distance) + 0.15")
    ap.add_argument("--no-critical-slice", action="store_true",
                    help="do NOT sample E=0; estimate saddles from the "
                         "flank loops' self-approach")
    args = ap.parse_args()

    model = get_model(**DEFAULT_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()
    oracle = dmm.make_oracle(coeffs, degs)

    delta = args.delta
    lo = oracle(-delta)
    hi = oracle(delta)

    lo_loop = lo.loops[0]
    hi_loop = hi.loops[0]
    key_lo = (lo.E, 0)
    key_hi = (hi.E, 0)

    if args.no_critical_slice:
        # estimate saddle chart positions from the flanks alone
        s_lo = estimate_saddles_from_flank(lo_loop)
        s_hi = estimate_saddles_from_flank(hi_loop)
        print(f"estimated saddles lo: {s_lo}")
        print(f"estimated saddles hi: {s_hi}")
        # match the two lists by chart proximity and average
        saddles = []
        used_hi = set()
        for p_lo in s_lo:
            best = None
            for q, p_hi in enumerate(s_hi):
                if q in used_hi:
                    continue
                d = np.sqrt(angular_dist(p_lo[0], p_hi[0]) ** 2
                            + angular_dist(p_lo[1], p_hi[1]) ** 2)
                if best is None or d < best[0]:
                    best = (d, q, p_hi)
            if best is not None:
                _, q, p_hi = best
                used_hi.add(q)
                pm = (float(np.angle(np.exp(1j * p_lo[0])
                                     + np.exp(1j * p_hi[0])) % TWO_PI),
                      float(np.angle(np.exp(1j * p_lo[1])
                                     + np.exp(1j * p_hi[1])) % TWO_PI))
                saddles.append(pm)
        if not saddles:
            raise RuntimeError("saddle estimation failed")
        print(f"averaged saddles: {saddles}")
    else:
        sc = oracle(0.0)
        saddles = [tuple(float(x) % TWO_PI for x in P)
                   for P in dmm.find_crossings(sc)]
        print(f"critical-slice crossings: {saddles}")

    # passage window radius: scale with the closest approach distance
    if args.R is not None:
        R = args.R
    else:
        th1 = np.asarray(lo_loop.th1) % TWO_PI
        th2 = np.asarray(lo_loop.th2) % TWO_PI
        dmin = 1e9
        for P in saddles:
            P = (float(P[0]) % TWO_PI, float(P[1]) % TWO_PI)
            d = np.sqrt(angular_dist(th1, P[0]) ** 2
                        + angular_dist(th2, P[1]) ** 2)
            dmin = min(dmin, float(d.min()))
        R = 3.0 * dmin + 0.15
    print(f"delta={delta}, R={R:.4f}")

    builder = dmm.MeshBuilder()
    builder.register_slice(lo)
    builder.register_slice(hi)

    all_rungs = []
    for P in saddles:
        P = (float(P[0]) % TWO_PI, float(P[1]) % TWO_PI)
        rungs, n_tris = build_saddle_fan(lo_loop, hi_loop, P, R,
                                         builder, key_lo, key_hi)
        print(f"  P={P}: fan tris={n_tris}, rungs={len(rungs)}")
        all_rungs.extend(rungs)

    # ---- collect passages on both loops for BOTH crossings ----------------
    lo_passages, hi_passages = [], []
    for P in saddles:
        P = (float(P[0]) % TWO_PI, float(P[1]) % TWO_PI)
        lo_passages.extend(sorted(passages_on_loop(lo_loop, P, R),
                                  key=lambda p: p["dmin"])[:2])
        hi_passages.extend(sorted(passages_on_loop(hi_loop, P, R),
                                  key=lambda p: p["dmin"])[:2])
    print(f"lo passages: {len(lo_passages)}, hi passages: {len(hi_passages)}")

    # sort passages cyclically by their start index on each loop
    lo_passages.sort(key=lambda p: p["s"])
    hi_passages.sort(key=lambda p: p["s"])

    # ---- rung lookup ------------------------------------------------------
    lo_to_hi = {}
    hi_to_lo = {}
    for (side_a, ia), (side_b, ib) in all_rungs:
        if side_a == 0:
            lo_to_hi[ia] = ib
            hi_to_lo[ib] = ia
        else:
            lo_to_hi[ib] = ia
            hi_to_lo[ia] = ib
    print(f"rungs total: {len(all_rungs)}, lo_to_hi {len(lo_to_hi)}")

    # ---- build FKU strips for the complement arcs -------------------------
    n_lo = lo_loop.n
    n_hi = hi_loop.n
    seen_strips = set()
    n_strips = 0
    for L, H in lo_to_hi.items():
        lo_arc = complement_arc_from_endpoint(n_lo, lo_passages, L)
        L2 = int(lo_arc[-1])
        if L2 not in lo_to_hi:
            print(f"  !! lo endpoint {L2} has no rung")
            continue
        H2 = lo_to_hi[L2]
        hi_arc = find_complement_arc(n_hi, hi_passages, H, H2)
        key = (min(L, L2), max(L, L2))
        if key in seen_strips:
            continue
        seen_strips.add(key)
        n_strips += 1
        dmm.arc_strip(lo_loop, lo_arc, hi_loop, hi_arc, builder,
                      key_lo, key_hi)

    print(f"FKU strips: {n_strips}")

    stats = dmm.validate(builder)
    print("\n=== validation (critical slab only) ===")
    for k, v in stats.items():
        print(f"  {k}: {v}")

    out = Path(__file__).resolve().parent / "gbz_critical_slab.pkl"
    import pickle
    with open(out, "wb") as f:
        pickle.dump({"verts": builder.vertex_array(),
                     "tris": np.array(builder.tris, dtype=int),
                     "coords": builder.coords,
                     "crit_info": builder.crit_info,
                     "stats": stats}, f)
    print(f"\nsaved -> {out}")


if __name__ == "__main__":
    main()
