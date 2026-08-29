'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-29
Copyright © Department of Physics, Tsinghua University. All rights reserved

Complete torus reconstruction prototype for the trivial (Hermitian) model.

Pipeline
--------
1.  bulk      : FKU strips between clean slices (certificate-gated).
2.  critical  : adaptive E-bisection until the flank-loop saddle gap
                (torus self-approach distance of the two hyperbola
                branches) is comparable to the local fibre sampling
                spacing.  The ladders are stitched with FKU and the
                final tiny gap is closed by the saddle fan + FKU
                complement strips (see demo_critical_slab.py).
3.  caps      : reuse demo_morse_mesh.grow_cap at the two spectrum
                edges.
4.  orient    : BFS over shared edges flips triangles into a globally
                consistent orientation.
5.  validate  : dmm.validate.

The critical slice E=0 is NOT needed: saddle chart positions are
estimated from ONE flank loop alone, by taking the torus midpoint of
the two U-tip points of closest self-approach.
'''

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from collections import Counter, defaultdict

import numpy as np

import demo_morse_mesh as dmm
from demo_critical_slab import (
    passages_on_loop, cyclic_range,
    complement_arc_from_endpoint, find_complement_arc,
)
from pygbz2d.core import CharPoly, LineSubset
from trivial_model import get_model, DEFAULT_PARAMS

TWO_PI = 2.0 * np.pi

MAX_STEP = 0.1          # bulk E-grid step, as in demo_morse_mesh
TOL_E = 1e-4            # hard floor for interval size (safety net)
MAX_DEPTH = 30

# exact-refinement controls
FTOL = 1e-12            # |f(E, beta1, beta2)| tolerance
DTOL = 1e-14            # |delta beta2| tolerance
MAX_NEWTON = 20

REFINE_STATS = {
    "newton_ok": 0,
    "fallback_ok": 0,
    "failed": 0,
    "max_f_refined": 0.0,
}


# ---------------------------------------------------------------------------
# Exact refinement of resampled points
# ---------------------------------------------------------------------------

def refine_point(E, mu1, theta1, b2_guess, charpoly):
    """Solve f(E, exp(mu1+i*theta1), beta2) = 0 with Newton from b2_guess.

    Returns (beta1, beta2, status) where status is one of
    'newton', 'fallback', 'failed'.
    """
    beta1 = np.exp(mu1 + 1j * theta1)
    b2 = complex(b2_guess)
    var = (E, beta1, b2)
    f = charpoly.eval_val(var)

    if abs(f) < FTOL:
        REFINE_STATS["newton_ok"] += 1
        REFINE_STATS["max_f_refined"] = max(REFINE_STATS["max_f_refined"],
                                            abs(f))
        return beta1, b2, "newton"

    for _ in range(MAX_NEWTON):
        part = charpoly.eval_partials(var)
        df = part[2]
        if abs(df) < 1e-30:
            break
        db2 = f / df
        b2 = b2 - db2
        var = (E, beta1, b2)
        f = charpoly.eval_val(var)
        if abs(db2) < DTOL and abs(f) < FTOL:
            REFINE_STATS["newton_ok"] += 1
            REFINE_STATS["max_f_refined"] = max(REFINE_STATS["max_f_refined"],
                                                abs(f))
            return beta1, b2, "newton"
        if abs(db2) < DTOL:
            break

    # fallback: solve all beta2 roots and pick the one closest to the
    # interpolated guess
    try:
        roots = charpoly.solve_roots_1d((0, 1), (E, beta1), (2,))
        if len(roots):
            roots = np.asarray(roots, dtype=complex)
            k = int(np.argmin(np.abs(roots - complex(b2_guess))))
            b2 = complex(roots[k])
            f = charpoly.eval_val((E, beta1, b2))
            if abs(f) < FTOL:
                REFINE_STATS["fallback_ok"] += 1
                REFINE_STATS["max_f_refined"] = max(
                    REFINE_STATS["max_f_refined"], abs(f))
                return beta1, b2, "fallback"
    except Exception:
        pass

    REFINE_STATS["failed"] += 1
    return beta1, complex(b2_guess), "failed"


# ---------------------------------------------------------------------------
# Uniform spatial resampling of loops
# ---------------------------------------------------------------------------

def loop_r6(loop):
    b1 = np.asarray(loop.b1, dtype=complex)
    b2 = np.asarray(loop.b2, dtype=complex)
    return np.column_stack([b1.real, b1.imag, b2.real, b2.imag])


def loop_spacing_r6(loop):
    p = loop_r6(loop)
    d = np.linalg.norm(np.roll(p, -1, axis=0) - p, axis=1)
    return float(np.median(d))


def _refine_resampled_arrays(loop, b1_new, b2_new, charpoly):
    """Refine interpolated (b1, b2) points to exact GBZ solutions."""
    if charpoly is None:
        return b1_new, b2_new
    mu1 = float(np.log(abs(loop.b1[0])))
    E = loop.E
    out_b1 = np.empty_like(b1_new)
    out_b2 = np.empty_like(b2_new)
    for k in range(len(b1_new)):
        th1 = float(np.angle(b1_new[k]))
        beta1, beta2, status = refine_point(E, mu1, th1,
                                            b2_new[k], charpoly)
        out_b1[k] = beta1
        out_b2[k] = beta2
    return out_b1, out_b2


def resample_loop(loop, spacing, charpoly=None):
    """Resample a closed loop to uniform chordal spacing in R^6, then
    refine every new point to an exact GBZ solution."""
    p = loop_r6(loop)
    n = len(p)
    seg = np.linalg.norm(np.roll(p, -1, axis=0) - p, axis=1)
    total = float(seg.sum())
    if total <= 0 or len(seg) == 0:
        return loop
    n_new = max(3, int(round(total / spacing)))
    s = np.concatenate([[0.0], np.cumsum(seg)])
    s_target = np.linspace(0.0, total, n_new, endpoint=False)
    idx = np.searchsorted(s, s_target, side="right") - 1
    idx = np.clip(idx, 0, n - 1)
    t = (s_target - s[idx]) / np.where(seg[idx] > 0, seg[idx], 1.0)
    j = (idx + 1) % n
    p_new = p[idx] + t[:, None] * (p[j] - p[idx])
    b1_new = p_new[:, 0] + 1j * p_new[:, 1]
    b2_new = p_new[:, 2] + 1j * p_new[:, 3]
    b1_new, b2_new = _refine_resampled_arrays(loop, b1_new, b2_new,
                                              charpoly)
    return dmm.Loop(np.angle(b1_new), np.angle(b2_new),
                    b1_new, b2_new, loop.E,
                    perm=np.arange(n_new))


def resample_loop_to_n(loop, n_target, charpoly=None):
    """Resample a closed loop to exactly n_target vertices, uniformly by
    chordal arc length in R^6, then refine every new point exactly."""
    p = loop_r6(loop)
    n = len(p)
    seg = np.linalg.norm(np.roll(p, -1, axis=0) - p, axis=1)
    total = float(seg.sum())
    n_target = max(3, int(n_target))
    if total <= 0 or n_target == n:
        return loop
    s = np.concatenate([[0.0], np.cumsum(seg)])
    s_target = np.linspace(0.0, total, n_target, endpoint=False)
    idx = np.searchsorted(s, s_target, side="right") - 1
    idx = np.clip(idx, 0, n - 1)
    t = (s_target - s[idx]) / np.where(seg[idx] > 0, seg[idx], 1.0)
    j = (idx + 1) % n
    p_new = p[idx] + t[:, None] * (p[j] - p[idx])
    b1_new = p_new[:, 0] + 1j * p_new[:, 1]
    b2_new = p_new[:, 2] + 1j * p_new[:, 3]
    b1_new, b2_new = _refine_resampled_arrays(loop, b1_new, b2_new,
                                              charpoly)
    return dmm.Loop(np.angle(b1_new), np.angle(b2_new),
                    b1_new, b2_new, loop.E,
                    perm=np.arange(n_target))


def resample_slice_pair(l_orig, r_orig, spacing, charpoly=None):
    """Resample the two flank slices to the SAME number of vertices per
    loop (matched sampling), which is what FKU complement strips need."""
    n_l = sum(l.n for l in l_orig.loops)
    n_r = sum(l.n for l in r_orig.loops)
    if n_l == 0 or n_r == 0:
        return l_orig, r_orig
    # common target N from the mean loop arc length
    lens = []
    for sl in (l_orig, r_orig):
        for l in sl.loops:
            p = loop_r6(l)
            seg = np.linalg.norm(np.roll(p, -1, axis=0) - p, axis=1)
            lens.append(float(seg.sum()))
    n_target = int(round(float(np.mean(lens)) / spacing))
    n_target = max(3, n_target)
    left = dmm.Slice(l_orig.E,
                     [resample_loop_to_n(l, n_target, charpoly)
                      for l in l_orig.loops],
                     is_critical=l_orig.is_critical)
    right = dmm.Slice(r_orig.E,
                      [resample_loop_to_n(l, n_target, charpoly)
                       for l in r_orig.loops],
                      is_critical=r_orig.is_critical)
    return left, right


def resample_slice(sl, spacing, charpoly=None):
    return dmm.Slice(sl.E,
                     [resample_loop(l, spacing, charpoly)
                      for l in sl.loops],
                     is_critical=sl.is_critical)


def global_spacing(slices):
    meds = [loop_spacing_r6(l) for s in slices for l in s.loops]
    return float(np.median(meds))


# ---------------------------------------------------------------------------
# Saddle estimation from one flank loop (no critical slice needed)
# ---------------------------------------------------------------------------

def angular_dist(a, b):
    return np.abs(np.angle(np.exp(1j * (np.asarray(a) - np.asarray(b)))))


def chart_dist(p, q):
    return float(np.sqrt(angular_dist(p[0], q[0]) ** 2
                         + angular_dist(p[1], q[1]) ** 2))


def estimate_saddles_from_flank(loop, k=8, cyc_min=20, d_max=1.0):
    """Return (saddles, gaps) estimated from one clean flank loop.

    saddles : list of (theta1, theta2) torus-chart positions.
    gaps    : self-approach distance between the two U-tips of each
              saddle (~ 2*sqrt(delta) for a Morse saddle).
    """
    th1 = np.asarray(loop.th1) % TWO_PI
    th2 = np.asarray(loop.th2) % TWO_PI
    n = len(th1)

    # torus self-distance matrix
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

    # local minima with large cyclic separation (different hyperbola
    # branches, not the two arms of the same U)
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            v = D[i, j]
            if v >= d_max:
                continue
            cyc = min(j - i, n - (j - i))
            if cyc < cyc_min:
                continue
            win = D[max(0, i - 3):i + 4, max(0, j - 3):j + 4]
            if v <= win.min():
                pairs.append((i, j, v))

    # deduplicate nearby index pairs
    uniq = []
    for i, j, v in sorted(pairs, key=lambda x: x[2]):
        if any(abs(u[0] - i) <= 4 and abs(u[1] - j) <= 4 for u in uniq):
            continue
        uniq.append((i, j, v))

    # each pair -> one saddle candidate = torus midpoint of the two tips
    cand = []
    for i, j, v in uniq:
        p1 = float(np.angle(np.exp(1j * th1[i]) + np.exp(1j * th1[j]))
                   % TWO_PI)
        p2 = float(np.angle(np.exp(1j * th2[i]) + np.exp(1j * th2[j]))
                   % TWO_PI)
        cand.append(((p1, p2), v))

    # merge duplicates (same saddle seen through the chart seam)
    clusters = []
    for p, g in cand:
        for c in clusters:
            if chart_dist(p, c["mean"]) < 0.5:
                c["pts"].append(p)
                c["gaps"].append(g)
                x = np.sum([np.exp(1j * q[0]) for q in c["pts"]])
                y = np.sum([np.exp(1j * q[1]) for q in c["pts"]])
                c["mean"] = (float(np.angle(x) % TWO_PI),
                             float(np.angle(y) % TWO_PI))
                break
        else:
            clusters.append({"pts": [p], "gaps": [g], "mean": p})

    saddles = [c["mean"] for c in clusters]
    gaps = [float(np.mean(c["gaps"])) for c in clusters]
    return saddles, gaps


def estimate_saddles_pair(left, right):
    """Estimate saddles from both flanks and average matched estimates."""
    s_lo, g_lo = estimate_saddles_from_flank(left.loops[0])
    s_hi, g_hi = estimate_saddles_from_flank(right.loops[0])
    out = []
    used = set()
    for a, (p, g) in enumerate(zip(s_lo, g_lo)):
        best = None
        for b, (q, h) in enumerate(zip(s_hi, g_hi)):
            if b in used:
                continue
            d = chart_dist(p, q)
            if best is None or d < best[0]:
                best = (d, b, q, h)
        if best is None:
            continue
        _, b, q, h = best
        used.add(b)
        pm = (float(np.angle(np.exp(1j * p[0]) + np.exp(1j * q[0]))
                    % TWO_PI),
              float(np.angle(np.exp(1j * p[1]) + np.exp(1j * q[1]))
                    % TWO_PI))
        out.append((pm, max(g, h)))
    return out


# ---------------------------------------------------------------------------
# Saddle fan with general critical-point (beta1, beta2) estimation
# ---------------------------------------------------------------------------

def build_saddle_fan(lo_loop, hi_loop, P, R, builder, key_lo, key_hi):
    """Fan P over the four boundary passages in angular order.

    The critical vertex coordinates are estimated from the flank
    passage tips (not from the critical slice): mu1 / mu2 are averaged
    log-moduli of the four U-tips, phases are the estimated chart
    position P.
    """
    lo_pass = sorted(passages_on_loop(lo_loop, P, R),
                     key=lambda p: p["dmin"])[:2]
    hi_pass = sorted(passages_on_loop(hi_loop, P, R),
                     key=lambda p: p["dmin"])[:2]
    if len(lo_pass) < 2 or len(hi_pass) < 2:
        raise RuntimeError(f"passage detection failed at P={P}: "
                           f"{len(lo_pass)} lo, {len(hi_pass)} hi")

    # --- estimate critical vertex (E, beta1, beta2) from the 4 tips ----
    tips = []          # (b1, b2) of each passage tip
    for passages, loop in ((lo_pass, lo_loop), (hi_pass, hi_loop)):
        for p in passages:
            d = np.linalg.norm(p["uv"] - np.array(P), axis=1)
            k = int(np.argmin(d))
            idx = int(p["idx"][k])
            tips.append((complex(loop.b1[idx]), complex(loop.b2[idx])))
    mu1_est = float(np.mean([np.log(abs(b1)) for b1, _ in tips]))
    mu2_est = float(np.mean([np.log(abs(b2)) for _, b2 in tips]))
    E_mid = 0.5 * (key_lo[0] + key_hi[0])
    beta1 = np.exp(mu1_est + 1j * P[0])
    beta2 = np.exp(mu2_est + 1j * P[1])
    P_gid = builder.add_vertex(E_mid, beta1, beta2, note=f"crit E~{E_mid:.5f}")

    # --- fan over each passage in loop order + rung triangles ------------
    # Each passage's boundary order IS its loop order, so fanning P over
    # that order always covers the passage wedge; the four wedges are
    # then connected by rung triangles between loop-order passage
    # ENDPOINTS (not angular extremes, which may be interior points).
    passages = [(0, p, key_lo) for p in lo_pass] + \
               [(1, p, key_hi) for p in hi_pass]

    def passage_angle_info(p):
        uv = p["uv"]
        raw = np.arctan2(uv[:, 1] - P[1], uv[:, 0] - P[0])
        mean = float(np.arctan2(np.sin(raw).mean(), np.cos(raw).mean()))
        unwrap = (raw - mean + np.pi) % TWO_PI - np.pi + mean
        s_pos = int(np.where(p["idx"] == p["s"])[0][0])
        e_pos = int(np.where(p["idx"] == p["e"])[0][0])
        a_s, a_e = float(unwrap[s_pos]), float(unwrap[e_pos])
        if a_s <= a_e:
            lo, hi = p["s"], p["e"]
        else:
            lo, hi = p["e"], p["s"]
        return mean, lo, hi

    info = {id(p): passage_angle_info(p) for _, p, _ in passages}
    passages.sort(key=lambda sp: info[id(sp[1])][0])

    n_tris = 0
    for side, p, key in passages:
        idx = p["idx"]
        for k in range(len(idx) - 1):
            a = builder.gid(key[0], key[1], int(idx[k]))
            b = builder.gid(key[0], key[1], int(idx[k + 1]))
            tri = (P_gid, a, b)
            if len({tri[0], tri[1], tri[2]}) == 3:
                builder.tris.append(tri)
                n_tris += 1

    rungs = []
    for k in range(4):
        side_a, p_a, key_a = passages[k]
        side_b, p_b, key_b = passages[(k + 1) % 4]
        _, _, hi_a = info[id(p_a)]
        _, lo_b, _ = info[id(p_b)]
        ga = builder.gid(key_a[0], key_a[1], int(hi_a))
        gb = builder.gid(key_b[0], key_b[1], int(lo_b))
        tri = (P_gid, ga, gb)
        if len({tri[0], tri[1], tri[2]}) == 3:
            builder.tris.append(tri)
            n_tris += 1
        rungs.append(((side_a, int(hi_a)), (side_b, int(lo_b))))
    return rungs, n_tris


def build_slab_patch(left, right, saddles_gaps, builder):
    """Saddle fans + FKU complement strips for the critical slab."""
    lo_loop = left.loops[0]
    hi_loop = right.loops[0]
    key_lo = (left.E, 0)
    key_hi = (right.E, 0)

    all_rungs = []
    radii = []
    for P, gap in saddles_gaps:
        R = 2.0 * gap + 0.02
        radii.append(R)
        rungs, n_tris = build_saddle_fan(lo_loop, hi_loop, P, R,
                                         builder, key_lo, key_hi)
        print(f"    slab P={P} gap={gap:.4f} R={R:.4f} "
              f"fan_tris={n_tris} rungs={len(rungs)}")
        all_rungs.extend(rungs)

    # passages on both loops for all saddles
    lo_passages, hi_passages = [], []
    for (P, _), R in zip(saddles_gaps, radii):
        lo_passages.extend(sorted(passages_on_loop(lo_loop, P, R),
                                  key=lambda p: p["dmin"])[:2])
        hi_passages.extend(sorted(passages_on_loop(hi_loop, P, R),
                                  key=lambda p: p["dmin"])[:2])
    lo_passages.sort(key=lambda p: p["s"])
    hi_passages.sort(key=lambda p: p["s"])

    lo_to_hi = {}
    for (side_a, ia), (side_b, ib) in all_rungs:
        if side_a == 0:
            lo_to_hi[ia] = ib
        else:
            lo_to_hi[ib] = ia

    n_lo = lo_loop.n
    n_hi = hi_loop.n
    seen = set()
    n_strips = 0
    for L, H in lo_to_hi.items():
        lo_arc = complement_arc_from_endpoint(n_lo, lo_passages, L)
        L2 = int(lo_arc[-1])
        if L2 not in lo_to_hi:
            continue
        H2 = lo_to_hi[L2]
        try:
            hi_arc = find_complement_arc(n_hi, hi_passages, H, H2)
        except ValueError:
            print(f"    !! complement arc not found: H={H} H2={H2}")
            print(f"       hi_passages={[(p['s'], p['e']) for p in hi_passages]}")
            print(f"       lo_passages={[(p['s'], p['e']) for p in lo_passages]}")
            print(f"       lo_to_hi={lo_to_hi}")
            raise
        key = (min(L, L2), max(L, L2))
        if key in seen:
            continue
        seen.add(key)
        n_strips += 1
        dmm.arc_strip(lo_loop, lo_arc, hi_loop, hi_arc, builder,
                      key_lo, key_hi)
    print(f"    slab strips: {n_strips}")


# ---------------------------------------------------------------------------
# Adaptive critical-interval narrowing
# ---------------------------------------------------------------------------

def loop_chart_spacing(loop):
    th1 = np.asarray(loop.th1) % TWO_PI
    th2 = np.asarray(loop.th2) % TWO_PI
    d = np.sqrt(angular_dist(th1, np.roll(th1, 1)) ** 2
                + angular_dist(th2, np.roll(th2, 1)) ** 2)
    return float(np.median(d))


def slab_ready(left, right, factor=2.0):
    """True when every flank loop's saddle gap is comparable to its
    local sampling spacing."""
    for sl in (left, right):
        for loop in sl.loops:
            _, gaps = estimate_saddles_from_flank(loop)
            if not gaps:
                return False
            h = loop_chart_spacing(loop)
            if max(gaps) > factor * h:
                return False
    return True


def stitch_fku(a, b, builder):
    for i in range(len(b.loops)):
        b.loops[i] = dmm.orient_pair(a.loops[0], b.loops[i])
    for i, j in dmm.match_components(a.loops, b.loops):
        dmm.strip_triangles(a.loops[i], b.loops[j], builder,
                            (a.E, i), (b.E, j))


def stitch_ladder(ladder, builder):
    for k in range(len(ladder) - 1):
        a, b = ladder[k], ladder[k + 1]
        stitch_fku(a, b, builder)


def process_interval(sa, sb, oracle, builder, far_count, spacing,
                     charpoly, tol_E=TOL_E, max_depth=MAX_DEPTH, depth=0):
    ok, diag = dmm.certificate(sa, sb)
    if ok and abs(sb.E - sa.E) > MAX_STEP * (1.0 + 1e-9):
        ok = False
        diag["reason"] = f"gap {abs(sb.E - sa.E):.3f} > {MAX_STEP}"
    if ok:
        stitch_fku(sa, sb, builder)
        return

    if abs(sb.E - sa.E) < tol_E or depth >= max_depth:
        # last resort: treat [sa,sb] as a slab if we can estimate saddles
        sg = estimate_saddles_pair(sa, sb)
        if sg:
            build_slab_patch(sa, sb, sg, builder)
        else:
            print(f"    !! fallback grow_patch at "
                  f"[{sa.E:+.5f},{sb.E:+.5f}]")
            dmm.grow_patch(sa, sb, oracle, builder, n_probes=0,
                           extra_pts=[], point_slices=[])
        return

    Em = 0.5 * (sa.E + sb.E)
    sm_orig = oracle(Em)

    if dmm.is_critical_slice(sm_orig, far_count=far_count):
        delta = 0.5 * min(Em - sa.E, sb.E - Em)
        ladder_lo, ladder_hi = [sa], [sb]
        left = right = None
        for _ in range(40):
            l_orig = oracle(Em - delta)
            r_orig = oracle(Em + delta)
            if (dmm.is_critical_slice(l_orig, far_count=far_count)
                    or dmm.is_critical_slice(r_orig, far_count=far_count)
                    or delta <= 1e-7):
                break
            left, right = resample_slice_pair(l_orig, r_orig, spacing,
                                              charpoly)
            builder.register_slice(left)
            builder.register_slice(right)
            ladder_lo.append(left)
            ladder_hi.append(right)
            if slab_ready(left, right):
                break
            delta *= 0.5
        if left is None:
            left, right = resample_slice_pair(oracle(Em - delta),
                                              oracle(Em + delta), spacing,
                                              charpoly)
            builder.register_slice(left)
            builder.register_slice(right)
            ladder_lo.append(left)
            ladder_hi.append(right)

        print(f"    critical [{sa.E:+.3f},{sb.E:+.3f}] narrowed to "
              f"[{left.E:+.6f},{right.E:+.6f}] "
              f"(ladder {len(ladder_lo)}+{len(ladder_hi)})")

        stitch_ladder(ladder_lo, builder)
        stitch_ladder(ladder_hi[::-1], builder)

        sg = estimate_saddles_pair(left, right)
        if sg:
            build_slab_patch(left, right, sg, builder)
        else:
            print(f"    !! slab estimation failed at "
                  f"[{left.E:+.6f},{right.E:+.6f}]; using grow_patch")
            dmm.grow_patch(left, right, oracle, builder, n_probes=0,
                           extra_pts=[], point_slices=[])
        return

    sm = resample_slice(sm_orig, spacing, charpoly)
    builder.register_slice(sm)
    process_interval(sa, sm, oracle, builder, far_count, spacing, charpoly,
                     tol_E, max_depth, depth + 1)
    process_interval(sm, sb, oracle, builder, far_count, spacing, charpoly,
                     tol_E, max_depth, depth + 1)


# ---------------------------------------------------------------------------
# Global orientation fix
# ---------------------------------------------------------------------------

def prune_isolated(builder):
    """Remove vertices that are not used by any triangle (probe-cloud
    leftovers) and remap triangle indices."""
    used = set()
    for t in builder.tris:
        used.update(t)
    old2new = {g: i for i, g in enumerate(sorted(used))}
    builder.coords = [builder.coords[g] for g in sorted(used)]
    builder.tris = [tuple(old2new[g] for g in t) for t in builder.tris]
    builder.crit_info = [(old2new[g], E, note)
                         for g, E, note in builder.crit_info
                         if g in old2new]
    builder.vid = {}       # stale after remap; no longer needed
    return len(old2new)


def orient_mesh(tris):
    """BFS over shared edges; flip neighbours so every interior edge is
    traversed in opposite directions by its two triangles."""
    tris = [list(t) for t in tris if len(set(t)) == 3]
    edge_tris = defaultdict(list)
    for ti, t in enumerate(tris):
        for k in range(3):
            a, b = t[k], t[(k + 1) % 3]
            edge_tris[(min(a, b), max(a, b))].append(ti)

    visited = [False] * len(tris)
    for seed in range(len(tris)):
        if visited[seed]:
            continue
        visited[seed] = True
        stack = [seed]
        while stack:
            ti = stack.pop()
            t = tris[ti]
            for k in range(3):
                a, b = t[k], t[(k + 1) % 3]
                e = (min(a, b), max(a, b))
                for tj in edge_tris[e]:
                    if tj == ti:
                        continue
                    t2 = tris[tj]
                    ia, ib = t2.index(a), t2.index(b)
                    # if t2 traverses a->b in the same direction, flip it
                    if (ia + 1) % 3 == ib:
                        tris[tj] = [t2[0], t2[2], t2[1]]
                    if not visited[tj]:
                        visited[tj] = True
                        stack.append(tj)
    return [tuple(t) for t in tris]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    model = get_model(**DEFAULT_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()
    oracle = dmm.make_oracle(coeffs, degs)
    charpoly = CharPoly(coeffs, degs)

    from trivial_line_cache import load
    data = load()
    E_grid, results = data["E_grid"], data["results"]

    all_slices = []
    for i, r in enumerate(results):
        lines = [s for s in r.subsets if isinstance(s, LineSubset)]
        if lines:
            all_slices.append(dmm.nexus_to_slice(
                float(E_grid[i]),
                dmm.build_nexus(dmm._sorted_lines(lines), tol=1e-12)))

    far_count = Counter(len(s.loops) for s in all_slices) \
        .most_common(1)[0][0]
    slices = [s for s in all_slices
              if not dmm.is_critical_slice(s, far_count=far_count)]
    print(f"{len(all_slices)} slices, {len(all_slices) - len(slices)} "
          f"critical removed, far_count={far_count}")

    # uniform spatial resampling of every clean slice
    spacing = global_spacing(slices)
    print(f"resampling spacing = {spacing:.5f}")
    slices = [resample_slice(s, spacing, charpoly) for s in slices]

    builder = dmm.MeshBuilder()
    for sl in slices:
        builder.register_slice(sl)

    for k in range(len(slices) - 1):
        sa, sb = slices[k], slices[k + 1]
        process_interval(sa, sb, oracle, builder, far_count, spacing,
                         charpoly)
        print(f"  interval [{sa.E:+.2f},{sb.E:+.2f}] "
              f"tris={len(builder.tris)}")

    # orient BEFORE caps (grow_region's outward_ok relies on orientation)
    builder.tris = orient_mesh(builder.tris)

    print("caps ...")
    dmm.grow_cap(slices[0], -1, oracle, builder)
    dmm.grow_cap(slices[-1], +1, oracle, builder)

    # close small cap slits, if any (guarded hole filler)
    dmm.fill_holes(builder)

    builder.tris = orient_mesh(builder.tris)

    # choose the global orientation that matches the (theta1, theta2)
    # chart: flip the whole mesh if the first triangle is negative
    if builder.tris:
        t = builder.tris[0]
        th = np.array([[np.angle(builder.coords[g][1]),
                        np.angle(builder.coords[g][2])] for g in t])
        q = np.angle(np.exp(1j * (th - th[0])))
        area = 0.5 * (q[1, 0] * q[2, 1] - q[2, 0] * q[1, 1])
        if area < 0:
            builder.tris = [(x, z, y) for x, y, z in builder.tris]

    n_before = len(builder.coords)
    n_kept = prune_isolated(builder)
    print(f"pruned {n_before - n_kept} isolated vertices "
          f"({n_kept} kept)")

    stats = dmm.validate(builder)
    print("\n=== validation ===")
    for k, v in stats.items():
        print(f"  {k}: {v}")

    print("\n=== exact refinement ===")
    for k, v in REFINE_STATS.items():
        print(f"  {k}: {v}")

    out = Path(__file__).resolve().parent / "gbz_torus_prototype.pkl"
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
