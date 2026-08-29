'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-28
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

"""Morse-aware contour stitching for LineSubset meshing (trivial model).

GENERALITY CONTRACT
-------------------
The fiber over E is a contour of a smooth function on a closed 2-manifold:

  * away from critical values it is a 1-manifold (union of simple closed
    curves) — bulk stitching only ever sees this case;
  * at critical values it may be an arbitrary graph (saddle, monkey
    saddle, homoclinic tangency, simultaneous events).  No particular
    transition form is assumed anywhere.

Design:
  1. bulk — FKU between matched cycle pairs; a certificate (displacement
     uniformity + stretch bound) decides whether the interval is free of
     critical values;
  2. certificate failure — adaptive E-bisection with the oracle
     (collect_GBZ_subsets) until the certificate passes or |dE| < tol;
  3. residual failing interval — critical points located as "bulges"
     where the flanking fibers diverge, glued by the UNIVERSAL fan rule:
     strand ends sorted by angle around the critical point in the
     intrinsic (theta1, theta2) chart, critical vertex fanned to
     consecutive ends.

Metric: chordal in (beta1, beta2) (Riemann sphere) — no unwrap assumptions.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pickle

import numpy as np
import networkx as nx

from pygbz2d.core import LineSubset, TWO_PI, to_sphere_r3
from demo_nexus import build_nexus
from demo_fku import fku_triangulate


# ---------------------------------------------------------------------------
# Slice representation
# ---------------------------------------------------------------------------

class Loop:
    """One component of a non-critical fiber, in walk order.

    perm maps array position -> ORIGINAL vertex index: orientation fixes
    (orient_pair) reverse loops IN PLACE after they may already be
    registered in the mesh vertex registry, which is keyed by the
    original index — every gid lookup must therefore go through perm.
    (A reversed-without-perm slice once made every later strip fetch the
    WRONG vertices: triangles connected unrelated points across the
    torus, producing hundreds of antipodal edges.)
    """

    def __init__(self, th1, th2, b1, b2, E, perm=None):
        self.th1 = np.asarray(th1, float)
        self.th2 = np.asarray(th2, float)
        self.b1 = np.asarray(b1, complex)
        self.b2 = np.asarray(b2, complex)
        self.E = E
        self.n = len(self.th1)
        self.perm = (np.arange(self.n) if perm is None
                     else np.asarray(perm, dtype=int))

    def reversed(self):
        return Loop(self.th1[::-1], self.th2[::-1],
                    self.b1[::-1], self.b2[::-1], self.E,
                    perm=self.perm[::-1])

    def chordal_r3(self):
        return np.concatenate([to_sphere_r3(self.b1),
                               to_sphere_r3(self.b2)], axis=1)

    def cyclic_arc(self, i, j):
        """Forward arc i->j (cyclic, inclusive).  Returns (Loop, index array
        into self) so triangle vertices can be mapped back exactly."""
        if j >= i:
            idx = np.arange(i, j + 1) % self.n
        else:
            idx = np.concatenate([np.arange(i, self.n), np.arange(0, j + 1)])
        sub = Loop(self.th1[idx], self.th2[idx], self.b1[idx],
                   self.b2[idx], self.E)
        return sub, idx.astype(int)


class Slice:
    def __init__(self, E, loops, is_critical=False):
        self.E = E
        self.loops = loops
        self.is_critical = is_critical



def _sorted_lines(lines):
    """Sort LineSubsets by median arg(beta2) mod 2pi before build_nexus.

    The nexus merge/walk structure depends on line order: unsorted lines
    produce a cycle whose walk order breaks FKU contour matching
    (measured: ~500 edges > 0.5 chordal per strip unsorted vs 0 sorted —
    a habit the original demo_build_mesh.main already had).
    """
    return sorted(lines, key=lambda L: float(np.median(
        np.angle(np.asarray(L.beta2_arr, dtype=complex)) % (2 * np.pi))))


def nexus_to_slice(E, G) -> Slice:
    """Nexus graph -> Slice of ordered Loops (degree != 2 marks critical)."""
    loops = []
    critical = False
    for comp in nx.connected_components(G):
        sub = G.subgraph(comp)
        if any(d != 2 for _, d in sub.degree()) or len(comp) < 3:
            critical = True
        start = min(comp)
        order = [start]
        prev, cur = None, start
        while len(order) < len(comp):
            nbrs = [u for u in sub.neighbors(cur) if u != prev]
            if not nbrs:
                break
            prev, cur = cur, nbrs[0]
            order.append(cur)
        th1 = np.array([G.nodes[n]['theta1'] for n in order])
        b2 = np.array([G.nodes[n]['beta2'] for n in order])
        b1 = np.array([G.nodes[n]['beta1'] for n in order])
        loops.append(Loop(th1, np.angle(b2), b1, b2, E))
    return Slice(E, loops, is_critical=critical)


# ---------------------------------------------------------------------------
# Distances / matching
# ---------------------------------------------------------------------------

def _d_matrix(a: Loop, b: Loop) -> np.ndarray:
    ra, rb = a.chordal_r3(), b.chordal_r3()
    d = ra[:, None, :] - rb[None, :, :]
    return (np.linalg.norm(d[:, :, :3], axis=2)
            + np.linalg.norm(d[:, :, 3:], axis=2))


def chamfer(a: Loop, b: Loop) -> float:
    D = _d_matrix(a, b)
    return 0.5 * (D.min(axis=1).mean() + D.min(axis=0).mean())


def match_components(loops_a, loops_b):
    from scipy.optimize import linear_sum_assignment
    la, lb = len(loops_a), len(loops_b)
    if la == 0 or lb == 0:
        return []
    cost = np.array([[chamfer(a, b) for b in loops_b] for a in loops_a])
    row, col = linear_sum_assignment(cost)
    return [(int(r), int(c)) for r, c in zip(row, col)]


def orient_pair(a: Loop, b: Loop) -> Loop:
    nn = _d_matrix(a, b).argmin(axis=1)
    steps = np.diff(nn)
    if np.sum(steps < 0) > np.sum(steps > 0):
        return b.reversed()
    return b


def dtw_match(a: Loop, b: Loop):
    """Cyclic cut + monotone DTW pairing between two oriented cycles.

    Returns (ia, ib, slope, disp): matched indices and per-step stretch /
    displacement diagnostics.
    """
    D = _d_matrix(a, b)
    m, n = a.n, b.n
    i0, j0 = np.unravel_index(D.argmin(), D.shape)
    A = np.concatenate([np.arange(i0, i0 + m) % m, [i0]])
    B = np.concatenate([np.arange(j0, j0 + n) % n, [j0]])
    M, N = m + 1, n + 1

    cost = np.full((M, N), np.inf)
    back = np.zeros((M, N), dtype=np.int8)
    cost[0, 0] = 0.0
    for i in range(M):
        for j in range(N):
            if i == 0 and j == 0:
                continue
            best, bdir = np.inf, 0
            if i > 0 and j > 0:
                v = cost[i - 1, j - 1] + D[A[i - 1], B[j - 1]]
                if v < best:
                    best, bdir = v, 2
            if i > 0:
                v = cost[i - 1, j] + 0.5 * D[A[i - 1], B[j - 1]]
                if v < best:
                    best, bdir = v, 1
            if j > 0:
                v = cost[i, j - 1] + 0.5 * D[A[i - 1], B[j - 1]]
                if v < best:
                    best, bdir = v, 0
            cost[i, j] = best
            back[i, j] = bdir

    ia, ib = [], []
    i, j = M - 1, N - 1
    while i > 0 or j > 0:
        if back[i, j] == 2 or i == 0:
            ia.append(int(A[i - 1]) if i > 0 else int(A[0]))
            ib.append(int(B[j - 1]))
            i, j = max(i - 1, 0), j - 1
        elif back[i, j] == 1:
            ia.append(int(A[i - 1]))
            ib.append(int(B[j - 1]))
            i -= 1
        else:
            ia.append(int(A[min(i, M - 2)]))
            ib.append(int(B[j - 1]))
            j -= 1
    ia.reverse(); ib.reverse()
    ia, ib = np.array(ia), np.array(ib)
    disp = D[ia, ib]

    slope = np.full(len(ia), np.nan)
    ra, rb = a.chordal_r3(), b.chordal_r3()
    for k in range(1, len(ia)):
        if ia[k] != ia[k - 1] and ib[k] != ib[k - 1]:
            da = np.linalg.norm(ra[ia[k]] - ra[ia[k - 1]])
            db = np.linalg.norm(rb[ib[k]] - rb[ib[k - 1]])
            if da > 1e-15:
                slope[k] = db / da
    return ia, ib, slope, disp


# ---------------------------------------------------------------------------
# Certificate
# ---------------------------------------------------------------------------

def certificate(sa: Slice, sb: Slice, *, disp_factor=4.0, slope_max=6.0):
    diag = {"Ea": sa.E, "Eb": sb.E,
            "n_loops": (len(sa.loops), len(sb.loops))}
    if sa.is_critical or sb.is_critical:
        diag["reason"] = "critical fiber"
        return False, diag
    if len(sa.loops) != len(sb.loops):
        diag["reason"] = "component count changed"
        return False, diag
    pairs = match_components(sa.loops, sb.loops)
    if len(pairs) != len(sa.loops):
        diag["reason"] = "component matching incomplete"
        return False, diag
    disp_ratio, slope_peak = 0.0, 0.0
    for i, j in pairs:
        lb = orient_pair(sa.loops[i], sb.loops[j])
        _, _, slope, disp = dtw_match(sa.loops[i], lb)
        med = max(np.median(disp), 1e-15)
        disp_ratio = max(disp_ratio, disp.max() / med)
        good = slope[np.isfinite(slope)]
        if len(good):
            slope_peak = max(slope_peak, float(np.percentile(good, 99)))
    diag["disp_ratio"] = round(disp_ratio, 2)
    diag["slope_p99"] = round(slope_peak, 2)
    if disp_ratio > disp_factor:
        diag["reason"] = f"displacement x{disp_ratio:.1f}"
        return False, diag
    if slope_peak > slope_max:
        diag["reason"] = f"stretch x{slope_peak:.1f}"
        return False, diag
    diag["reason"] = "ok"
    return True, diag


# ---------------------------------------------------------------------------
# Oracle
# ---------------------------------------------------------------------------

def make_oracle(coeffs, degs):
    import pygbz2d.amoeba as bfa

    def oracle(E):
        res = bfa.collect_GBZ_subsets(coeffs, degs, E, 0.0, debug_mode=False)
        lines = [s for s in res.subsets if isinstance(s, LineSubset)]
        if not res.success or not lines:
            return Slice(E, [], is_critical=True)
        return nexus_to_slice(E, build_nexus(_sorted_lines(lines),
                                          tol=1e-12))
    return oracle


# ---------------------------------------------------------------------------
# Global mesh assembly
# ---------------------------------------------------------------------------

class MeshBuilder:
    """Vertex registry keyed by (E, loop_id, vertex_id) + triangle list."""

    def __init__(self):
        self.vid = {}
        self.coords = []       # (E, b1, b2)
        self.tris = []
        self.crit_info = []    # provenance of inserted critical points

    def register_slice(self, sl):
        for j, loop in enumerate(sl.loops):
            for pos in range(loop.n):
                k = int(loop.perm[pos])   # stable original vertex index
                key = (sl.E, j, k)
                if key not in self.vid:
                    self.vid[key] = len(self.coords)
                    self.coords.append((sl.E, complex(loop.b1[k]),
                                        complex(loop.b2[k])))

    def gid(self, E, j, k):
        return self.vid[(E, j, k)]

    def add_vertex(self, E, b1, b2, note=""):
        gid = len(self.coords)
        self.coords.append((E, complex(b1), complex(b2)))
        self.crit_info.append((gid, E, note))
        return gid

    def vertex_array(self):
        return np.array([[E.real, E.imag, b1.real, b1.imag, b2.real, b2.imag]
                         for E, b1, b2 in self.coords])


# ---------------------------------------------------------------------------
# Bulk strip: FKU between matched closed cycles (duplicate-start trick)
# ---------------------------------------------------------------------------

def strip_triangles(a: Loop, b: Loop, builder, key_a, key_b):
    """FKU strip between matched cycles.  key_* = (E, loop_id)."""
    D = _d_matrix(a, b)
    i0, j0 = np.unravel_index(D.argmin(), D.shape)
    A = np.concatenate([np.arange(i0, i0 + a.n) % a.n, [i0]])
    B = np.concatenate([np.arange(j0, j0 + b.n) % b.n, [j0]])
    La = LineSubset(E=a.E, mu1=float(np.log(abs(a.b1[0]))),
                    theta1_arr=a.th1[A], beta2_arr=a.b2[A])
    Lb = LineSubset(E=b.E, mu1=float(np.log(abs(b.b1[0]))),
                    theta1_arr=b.th1[B], beta2_arr=b.b2[B])
    _v, tris, info = fku_triangulate(La, Lb)

    M, N = a.n + 1, b.n + 1   # FKU sees arrays WITH closing duplicates

    def gmap(g):
        if g < M:
            k = int(A[g]) if g < M - 1 else int(A[0])
            return builder.gid(key_a[0], key_a[1], int(a.perm[k]))
        h = g - M
        k = int(B[h]) if h < N - 1 else int(B[0])
        return builder.gid(key_b[0], key_b[1], int(b.perm[k]))

    for t in tris:
        builder.tris.append(tuple(gmap(int(x)) for x in t))
    return info


def arc_strip(a: Loop, idx_a, b: Loop, idx_b, builder, key_a, key_b):
    """FKU strip between two ARCS (index arrays into their parent loops).

    demo_fku._oriented may reverse a polyline with descending theta1; the
    flip is detected here and the index mapping is reversed accordingly.
    """
    La = LineSubset(E=a.E, mu1=float(np.log(abs(a.b1[0]))),
                    theta1_arr=a.th1[idx_a], beta2_arr=a.b2[idx_a])
    Lb = LineSubset(E=b.E, mu1=float(np.log(abs(b.b1[0]))),
                    theta1_arr=b.th1[idx_b], beta2_arr=b.b2[idx_b])
    verts, tris, info = fku_triangulate(La, Lb)
    m = len(idx_a)
    # demo_fku._oriented silently reverses descending-theta1 polylines;
    # detect by comparing FKU's endpoint thetas against ours
    flip_a = abs(verts[0][1] - a.th1[idx_a][-1]) < \
        abs(verts[0][1] - a.th1[idx_a][0])
    flip_b = abs(verts[m][1] - b.th1[idx_b][-1]) < \
        abs(verts[m][1] - b.th1[idx_b][0])
    ia = idx_a[::-1] if flip_a else idx_a
    ib = idx_b[::-1] if flip_b else idx_b
    for t in tris:
        vs = []
        for x in t:
            x = int(x)
            if x < m:
                vs.append(builder.gid(key_a[0], key_a[1],
                                      int(a.perm[ia[x]])))
            else:
                vs.append(builder.gid(key_b[0], key_b[1],
                                      int(b.perm[ib[x - m]])))
        builder.tris.append(tuple(vs))
    return info


# ---------------------------------------------------------------------------
# Critical-value handling
# ---------------------------------------------------------------------------

def self_touches(sl: Slice, k_local: float = 0.3, adj: int = 4):
    """Genuine near-self-intersections of a fiber.

    Criterion: a non-adjacent sample pair closer than k_local * (LOCAL
    spacing).  Local spacing matters because the solver refines sampling
    at hairpin turning points, where the two legs of a smooth U-turn are
    legitimately close (absolute-distance criteria misfire there); a true
    critical crossing instead brings unrelated strands far closer than
    the local sampling scale.
    """
    from scipy.spatial import cKDTree
    if not sl.loops:
        return []
    pts = np.vstack([l.chordal_r3() for l in sl.loops])
    owner = np.concatenate([[j] * l.n for j, l in enumerate(sl.loops)])
    base = np.concatenate([[sum(l2.n for l2 in sl.loops[:j])]
                           for j in range(len(sl.loops))])
    # per-point local spacing: mean of the two incident segment lengths,
    # computed within each loop (cyclically)
    local_sp = np.zeros(len(pts))
    for j, l in enumerate(sl.loops):
        p = l.chordal_r3()
        d = np.linalg.norm(np.diff(np.vstack([p, p[:1]]), axis=0), axis=1)
        local_sp[base[j]:base[j] + l.n] = 0.5 * (d + np.roll(d, 1))
    global_sp = float(np.median(local_sp))
    tree = cKDTree(pts)
    pairs = tree.query_pairs(10 * global_sp, output_type="ndarray")
    out = []
    for a, b in pairs:
        if owner[a] == owner[b]:
            ia, ib = a - base[owner[a]], b - base[owner[b]]
            n = sl.loops[owner[a]].n
            if min(abs(ia - ib), n - abs(ia - ib)) <= adj:
                continue
        d = float(np.linalg.norm(pts[a] - pts[b]))
        if d < k_local * max(local_sp[a], local_sp[b]):
            out.append((int(a), int(b), d))
    return out


def is_critical_slice(sl: Slice, far_count=None):
    if sl.is_critical or not sl.loops:
        return True
    if far_count is not None and len(sl.loops) != far_count:
        return True
    return len(self_touches(sl)) > 0


def _unwrapped_polyline(loop: Loop) -> np.ndarray:
    """Loop samples lifted to the universal cover (continuity unwrapping).

    Degenerate-fiber LineSubsets can come back in a jumbled sample order
    (observed on the E=0 fiber: theta1 jumps of ~pi between consecutive
    samples).  Consecutive-order continuity is checked; on failure the
    samples are re-sorted by theta1 (valid whenever the loop is a graph
    over theta1 — the generic case for GBZ wall lines) and re-checked.
    """
    th1 = loop.th1 % TWO_PI
    th2 = loop.th2 % TWO_PI

    def order_ok(t1, t2):
        d = np.linalg.norm(np.diff(np.column_stack(
            [np.unwrap(t1), np.unwrap(t2)]), axis=0), axis=1)
        return len(d) == 0 or (np.median(d) < 0.5 and d.max() < 1.0)

    if not order_ok(th1, th2):
        srt = np.argsort(th1, kind="stable")
        if order_ok(th1[srt], th2[srt]):
            th1, th2 = th1[srt], th2[srt]
    # outlier removal: a single off-strand sample (observed on the E=0
    # degenerate fiber — one point per line sitting ~pi off the strand)
    # breaks polyline continuity exactly at the crossing detection; drop
    # samples whose BOTH incident steps exceed 1.5x the median step
    if len(th1) > 4:
        t1u, t2u = np.unwrap(th1), th2.copy()
        d = np.linalg.norm(np.diff(np.column_stack(
            [t1u, np.unwrap(t2u)]), axis=0), axis=1)
        med = np.median(d)
        bad = np.zeros(len(th1), dtype=bool)
        dd = np.concatenate([[np.inf], d, [np.inf]])
        bad = (dd[:-1] > 1.5 * med) & (dd[1:] > 1.5 * med)
        if bad.any() and not bad.all():
            th1, th2 = th1[~bad], th2[~bad]
    # seam handling: rotate to start after the largest theta1 gap so the
    # unwrapped polyline is one continuous strand (incl. past 2pi)
    th1u = np.unwrap(th1)
    n = len(th1u)
    if n > 2:
        gaps = np.diff(th1u)
        wrap_gap = TWO_PI - (th1u[-1] - th1u[0])
        g = np.append(gaps, wrap_gap)
        start = (int(np.argmax(g)) + 1) % n
        rot = np.concatenate([np.arange(start, n), np.arange(0, start)])
        th1, th2 = th1[rot], th2[rot]
    return np.column_stack([np.unwrap(th1), np.unwrap(th2)])


def _polyline_intersections(A: np.ndarray, B: np.ndarray, skip_adjacent):
    """Segment intersections of two open polylines in the plane."""
    out = []
    na, nb = len(A), len(B)
    for i in range(na - 1):
        p1, p2 = A[i], A[i + 1]
        r = p2 - p1
        j0 = i + 2 if skip_adjacent else 0
        for j in range(j0, nb - 1):
            if skip_adjacent and i == 0 and j == nb - 2:
                continue  # wrap-adjacent first/last segments
            q1, q2 = B[j], B[j + 1]
            s = q2 - q1
            denom = r[0] * s[1] - r[1] * s[0]
            if abs(denom) < 1e-14:
                continue
            qp = q1 - p1
            t = (qp[0] * s[1] - qp[1] * s[0]) / denom
            u = (qp[0] * r[1] - qp[1] * r[0]) / denom
            if 0.0 <= t <= 1.0 and 0.0 <= u <= 1.0:
                out.append(p1 + t * r)
    return out


def find_crossings(sl: Slice):
    """Strand crossings of a fiber on the (theta1, theta2) torus.

    Each loop is lifted to the universal cover; intersections are taken
    against all 2π-translate copies (wrapping crossings included), then
    deduplicated modulo the torus.
    """
    polys = [_unwrapped_polyline(l) for l in sl.loops]
    raw = []
    for a in range(len(polys)):
        for b in range(a, len(polys)):
            for s1 in (-TWO_PI, 0.0, TWO_PI):
                for s2 in (-TWO_PI, 0.0, TWO_PI):
                    if a == b and (s1, s2) == (0.0, 0.0):
                        continue  # identical copies never "cross"
                    shifted = polys[b] + np.array([s1, s2])
                    raw.extend(_polyline_intersections(polys[a], shifted,
                                                       skip_adjacent=(a == b)))
    uniq = []
    for c in raw:
        c_mod = c % TWO_PI
        if not any(np.abs(np.angle(np.exp(1j * (c_mod - d)))).max() < 0.2
                   for d in uniq):
            uniq.append(c_mod)
    return uniq


def _passage_ends(loop: Loop, P_th, k_spacing: float = 2.5,
                  min_arc: int = 3):
    """Ends of the loop's passage through P's neighborhood.

    Radius k_spacing * local spacing (small: on a short side of a
    diamond-shaped fiber the passages around two critical points must not
    overlap, or the complementary remainder arcs degenerate to a single
    vertex and FKU collapses to a full-arc fan).  The passage is shrunk
    until it spans at least min_arc vertices.
    """
    P_b1 = np.exp(1j * (P_th[0] % (2 * np.pi)))
    P_b2 = np.exp(1j * (P_th[1] % (2 * np.pi)))
    r3p = np.concatenate([to_sphere_r3(np.array([P_b1]))[0],
                          to_sphere_r3(np.array([P_b2]))[0]])
    r3 = loop.chordal_r3()
    d = np.linalg.norm(r3 - r3p[None, :], axis=1)
    dmin = float(d.min())
    if dmin > 1.0:
        return None
    seg = np.linalg.norm(np.diff(np.vstack([r3, r3[:1]]), axis=0), axis=1)
    spacing = float(np.median(seg))
    near = set(np.where(d <= dmin + k_spacing * spacing)[0].tolist())
    a0 = int(np.argmin(d))
    n = loop.n
    if len(near) < 2:
        near = {a0, (a0 + 1) % n}
    lo_i = hi_i = a0
    while (lo_i - 1) % n in near:
        lo_i = (lo_i - 1) % n
    while (hi_i + 1) % n in near:
        hi_i = (hi_i + 1) % n
    while (hi_i - lo_i) % n + 1 < min_arc:
        lo_i = (lo_i - 1) % n
        hi_i = (hi_i + 1) % n
    return (lo_i, hi_i)


# ---------------------------------------------------------------------------
# Critical interval via GREEDY DELAUNAY GROWTH from the existing mesh
#
# The FKU-grown bulk meshes end at the flank loops; growth continues them
# triangle by triangle over the refined critical cloud: for each boundary
# (front) edge, pick the candidate third point with smallest circumradius
# (computable from the three chordal distances alone — no chart, no
# tangent estimation, torus-periodicity safe), constrained to the outward
# half-space defined by the edge's existing adjacent triangle.  Fronts
# from both flanks advance, meet, and close — the seam problem does not
# exist because there is no separate patch.  Vertices are never moved:
# degenerate configurations are resolved by deterministic tie-breaks.
# ---------------------------------------------------------------------------

def _r3_of_coord(coord):
    from pygbz2d.core import to_sphere_r3
    E, b1, b2 = coord
    return np.concatenate([to_sphere_r3(np.array([b1]))[0],
                           to_sphere_r3(np.array([b2]))[0]])


def grow_region(seed_slices, probe_Es, builder, oracle, extra_pts=(),
                point_slices=(), verbose=True):
    """Greedy Delaunay growth over a cloud built from seed slices, oracle
    probes, and extra vertices.

    Seed front = boundary edges of the CURRENT mesh on the seed slices —
    used both for interior critical patches (two seed slices, probes in
    between) and spectrum-edge caps (one seed slice, probes toward the
    edge, the cap/extremum point in extra_pts).
    """
    import heapq
    from collections import Counter, defaultdict
    from scipy.spatial import cKDTree

    seed_Es = {sl.E for sl in seed_slices}
    slices = list(seed_slices)
    # point-only members: pinched (self-touching) slices are useless as
    # FKU boundaries but perfectly good SURFACE SAMPLES — the growth
    # cloud only cares about positions, so near-critical density joins
    # the cloud without any loop semantics
    for sl in point_slices:
        builder.register_slice(sl)
        slices.append(sl)
    for Em in probe_Es:
        sm = oracle(Em)
        # register ONLY usable probes — registering a filtered slice
        # leaves its vertices isolated in the mesh (measured: 918 such
        # orphans once).  The filter includes the self-touch check: a
        # pinched loop (strands closer than the local sampling) cannot
        # be triangulated across and only poisons the cloud.
        if (not sm.is_critical) and sm.loops and not self_touches(sm):
            builder.register_slice(sm)
            slices.append(sm)
    cloud_gids = []
    for sl in slices:
        for j, loop in enumerate(sl.loops):
            for k in range(loop.n):
                cloud_gids.append(builder.gid(sl.E, j,
                                              int(loop.perm[k])))
    for E, b1, b2 in extra_pts:
        cloud_gids.append(builder.add_vertex(E, b1, b2,
                                             note="cap point"))
    r3_all = np.array([_r3_of_coord(builder.coords[g]) for g in cloud_gids])
    if point_slices:
        # cloud-level dedupe: pinch probes at |E| below the strand-separation
        # scale create near-coincident points across slices; they poison the
        # emptiness test.  SEED vertices are never dropped (the front must
        # keep every flank vertex — a merged-away seed vertex once raised
        # KeyError in the seed-front lookup); non-seed points are dropped
        # only when a survivor already sits within tolerance.
        from scipy.spatial import cKDTree as _KDT
        _t = _KDT(r3_all)
        _nn, _ = _t.query(r3_all, k=2)
        _tol = 0.6 * float(np.median(_nn[:, 1]))
        seed_gids = set()
        for sl in seed_slices:
            for j, loop in enumerate(sl.loops):
                for pos in range(loop.n):
                    seed_gids.add(builder.gid(sl.E, j, int(loop.perm[pos])))
        is_seed = np.array([g in seed_gids for g in cloud_gids])
        keep = list(np.where(is_seed)[0])
        kept_r3 = r3_all[keep]
        _kt = _KDT(kept_r3)
        for idx in range(len(r3_all)):
            if is_seed[idx]:
                continue
            near = _kt.query_ball_point(r3_all[idx], _tol)
            if not near:
                keep.append(idx)
                _kt = _KDT(r3_all[keep])     # rebuild lazily (small clouds)
        keep = sorted(set(keep))
        cloud_gids = [cloud_gids[i] for i in keep]
        r3 = r3_all[keep]
        n_pts = len(cloud_gids)
        if n_pts < 4:
            return 0, 0
        tree = cKDTree(r3)
        nn_dist, _ = tree.query(r3, k=2)
    else:
        r3 = r3_all
        n_pts = len(cloud_gids)
        tree = cKDTree(r3)
        nn_dist, _ = tree.query(r3, k=2)

    # ---- current edge incidence / adjacent triangles of the whole mesh
    inc = Counter()
    adj_tri = defaultdict(list)
    for t in builder.tris:
        if len(set(t)) < 3:
            continue
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[0], t[2])):
            e = (min(a, b), max(a, b))
            inc[e] += 1
            adj_tri[e].append(t)

    # ---- seed front: boundary edges lying on the seed slices
    idx_of_gid = {g: i for i, g in enumerate(cloud_gids)}
    front = set()
    for e, c in inc.items():
        if c == 1 and all(builder.coords[v][0] in seed_Es
                          for v in e):
            front.add((idx_of_gid[e[0]], idx_of_gid[e[1]]))
    if verbose:
        rng = " ".join(f"{sl.E:+.5f}" for sl in seed_slices)
        print(f"    [grow] seeds [{rng}]: {n_pts} cloud pts, "
              f"{len(front)} seed front edges")

    def cdist(i, j):
        return float(np.linalg.norm(r3[i, :3] - r3[j, :3])
                     + np.linalg.norm(r3[i, 3:] - r3[j, 3:]))

    def circumradius(dab, dbc, dca):
        s = 0.5 * (dab + dbc + dca)
        area2 = s * (s - dab) * (s - dbc) * (s - dca)
        if area2 <= 1e-30:
            return np.inf
        return dab * dbc * dca / (4.0 * np.sqrt(area2))

    def tangent_frame(i, j):
        """2D frame at front edge (i,j) from its existing adjacent
        triangle (the mesh itself supplies the local tangent plane)."""
        e = (min(cloud_gids[i], cloud_gids[j]),
             max(cloud_gids[i], cloud_gids[j]))
        tris = adj_tri.get(e)
        if not tris:
            return None
        t = tris[-1]
        l_gid = next(v for v in t if v not in (cloud_gids[i],
                                               cloud_gids[j]))
        rl = (r3[idx_of_gid[l_gid]] if l_gid in idx_of_gid
              else _r3_of_coord(builder.coords[l_gid]))
        ri, rj = r3[i], r3[j]
        e1 = rj - ri
        ne1 = np.linalg.norm(e1)
        if ne1 < 1e-15:
            return None
        e1 = e1 / ne1
        w = rl - ri
        e2 = w - (w @ e1) * e1
        ne2 = np.linalg.norm(e2)
        if ne2 < 1e-15:
            return None
        return ri, e1, e2 / ne2

    def _meb(i, j, k):
        """(center, R) of the minimal enclosing ball of 3 cloud points.

        The chordal metric lives in the flat Euclidean R^6 (Riemann
        sphere x Riemann sphere), so the MEB is exact: half the longest
        side for obtuse triangles, the circumcenter otherwise.  No
        tangent planes, no projections — the projected circle this
        replaces misfired exactly at the high-curvature spots (saddle
        corridors, cap cones) and rejected correct triangles.
        """
        pts = [r3[i], r3[j], r3[k]]
        sides = [(np.linalg.norm(pts[b] - pts[a]), a, b)
                 for a, b in ((0, 1), (1, 2), (0, 2))]
        c_len, ca, cb = max(sides)
        others = [s[0] for s in sides if s is not max(sides)]
        if c_len ** 2 >= others[0] ** 2 + others[1] ** 2 + 1e-30:
            center = 0.5 * (pts[ca] + pts[cb])
            return center, 0.5 * c_len
        p0 = pts[0]
        u = pts[1] - p0
        v = pts[2] - p0
        A = 2.0 * np.array([[u @ u, u @ v], [u @ v, v @ v]])
        b = np.array([pts[1] @ pts[1] - p0 @ p0,
                      pts[2] @ pts[2] - p0 @ p0])
        try:
            alpha, beta = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            center = 0.5 * (pts[ca] + pts[cb])
            return center, 0.5 * c_len
        center = p0 + alpha * u + beta * v
        return center, float(np.linalg.norm(center - p0))

    def empty_ball(i, j, k):
        """Emptiness in the FRONT EDGE's tangent frame (best measured).

        Two ambient variants were tried and were strictly WORSE here:
        the minimal-enclosing-ball rejects every skinny cross-row bridge
        (same-row neighbours intrude into the half-longest-side ball:
        3098 leftover cap edges), and the own-plane circumsphere lost
        the COHERENT comparison frame the front-edge tangent plane
        provides to competing candidates (7637 vs 78 hole edges).  The
        frame comes from the mesh itself (the edge's adjacent triangle).
        """
        frame = tangent_frame(i, j)
        if frame is None:
            return False
        ri, e1, e2 = frame

        def proj(p):
            q = p - ri
            return np.array([q @ e1, q @ e2])

        pi, pj, pk = proj(r3[i]), proj(r3[j]), proj(r3[k])
        ax, ay = pi; bx, by = pj; cx, cy = pk
        d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        if abs(d) < 1e-15:
            return False
        ux = ((ax * ax + ay * ay) * (by - cy)
              + (bx * bx + by * by) * (cy - ay)
              + (cx * cx + cy * cy) * (ay - by)) / d
        uy = ((ax * ax + ay * ay) * (cx - bx)
              + (bx * bx + by * by) * (ax - cx)
              + (cx * cx + cy * cy) * (bx - ax)) / d
        center = np.array([ux, uy])
        rad2 = ((pi - center) ** 2).sum()
        c3 = ri + ux * e1 + uy * e2
        margin = 3.0 * max(nn_dist[i, 1], nn_dist[j, 1])
        for m in tree.query_ball_point(c3, np.sqrt(rad2) + margin):
            m = int(m)
            if m in (i, j, k):
                continue
            pm = proj(r3[m])
            if ((pm - center) ** 2).sum() < rad2 - 1e-12:
                return False
        return True

    def outward_ok(i, j, k):
        """k on the outward side of front edge (i,j).

        The edge's existing adjacent triangle supplies the local tangent
        frame and the inward direction; the candidate must project to the
        opposite side (this is the no-flip / no-back-growth constraint).
        """
        e = (min(cloud_gids[i], cloud_gids[j]),
             max(cloud_gids[i], cloud_gids[j]))
        tris = adj_tri.get(e)
        if not tris:
            return False
        t = tris[0]
        l_gid = next(v for v in t if v not in (cloud_gids[i],
                                               cloud_gids[j]))
        l = idx_of_gid.get(l_gid)
        if l is None:
            # adjacent triangle vertex outside the cloud (the FKU side):
            # use its stored coordinate directly
            rl = _r3_of_coord(builder.coords[l_gid])
        else:
            rl = r3[l]
        ri, rj = r3[i], r3[j]
        e1 = rj - ri
        ne1 = np.linalg.norm(e1)
        if ne1 < 1e-15:
            return False
        e1 = e1 / ne1
        w = rl - ri
        e2 = w - (w @ e1) * e1
        ne2 = np.linalg.norm(e2)
        if ne2 < 1e-15:
            return False
        e2 = e2 / ne2
        y_l = (rl - ri) @ e2
        y_k = (r3[k] - ri) @ e2
        return y_k * y_l < 0.0 and abs(y_k) > 1e-12

    existing_tris = {tuple(sorted(t)) for t in builder.tris
                     if len(set(t)) == 3}

    def local_inc(i, j):
        return inc.get((min(cloud_gids[i], cloud_gids[j]),
                        max(cloud_gids[i], cloud_gids[j])), 0)

    heap = []

    def push_edge(ij):
        i, j = ij
        dij = cdist(i, j)
        # SECTOR-QUOTA candidate selection.  A plain k-nearest pool is
        # density-biased: along-strand samples outnumber cross-strand
        # ones, so bridging candidates never enter the heap and the
        # front stalls into long SLITS along the corridors (measured:
        # the 4 saddle "holes" are 16-20-edge slits).  Quotas per angular
        # sector of the outward half-plane make coverage
        # density-independent.
        frame = tangent_frame(i, j)
        if frame is None:
            return
        ri, e1, e2 = frame
        e_gid = (min(cloud_gids[i], cloud_gids[j]),
                 max(cloud_gids[i], cloud_gids[j]))
        tris_adj = adj_tri.get(e_gid)
        if not tris_adj:
            return
        t0 = tris_adj[-1]
        l_gid = next(v for v in t0 if v not in (cloud_gids[i],
                                                cloud_gids[j]))
        rl = (r3[idx_of_gid[l_gid]] if l_gid in idx_of_gid
              else _r3_of_coord(builder.coords[l_gid]))
        y_l = (rl - ri) @ e2
        sgn = 1.0 if y_l >= 0 else -1.0     # inward side sign

        _, pool = tree.query(0.5 * (r3[i] + r3[j]), k=min(96, n_pts))
        pool = np.atleast_1d(pool)
        mid3 = 0.5 * (r3[i] + r3[j])
        cand = []
        for m in pool:
            m = int(m)
            if m in (i, j):
                continue
            q = r3[m] - ri
            x, y = q @ e1, q @ e2
            y_out = -sgn * y                # outward-positive ordinate
            if y_out <= 1e-12:              # outward half-plane only
                continue
            cand.append((np.hypot(x - (mid3 - ri) @ e1,
                                  y - (mid3 - ri) @ e2),
                         np.arctan2(y_out, x), m))
        n_sec, per_sec = 5, 5
        sec_bins = [[] for _ in range(n_sec)]
        for d2, ang, m in sorted(cand):
            s = min(n_sec - 1, int((ang % np.pi) / np.pi * n_sec))
            if len(sec_bins[s]) < per_sec:
                sec_bins[s].append((d2, m))
        picked = {m for b in sec_bins for _, m in b}
        for m in picked:
            k = int(m)
            R_c = circumradius(dij, cdist(j, k), cdist(k, i))
            if np.isfinite(R_c):
                heapq.heappush(heap, (R_c, (i, j, k)))

    for e in list(front):
        push_edge(e)

    n_added = 0
    while heap:
        R_c, (i, j, k) = heapq.heappop(heap)
        if (i, j) not in front:
            continue
        if local_inc(i, k) > 1 or local_inc(j, k) > 1:
            continue                       # would exceed manifold valence
        tri_gids = tuple(sorted((cloud_gids[i], cloud_gids[j],
                                 cloud_gids[k])))
        if tri_gids in existing_tris:
            continue
        if not empty_ball(i, j, k):
            continue                       # ambient emptiness (no overlap)
        # accept
        new_tri = (cloud_gids[i], cloud_gids[j], cloud_gids[k])
        builder.tris.append(new_tri)
        existing_tris.add(tri_gids)
        n_added += 1
        for a, b in ((i, j), (j, k), (i, k)):
            e = (min(cloud_gids[a], cloud_gids[b]),
                 max(cloud_gids[a], cloud_gids[b]))
            inc[e] += 1
            adj_tri[e].append(new_tri)
            front.discard((a, b))
            front.discard((b, a))
            if inc[e] == 1:
                front.add((a, b))
                push_edge((a, b))
    leftover = len(front)
    if verbose:
        print(f"    [grow] added {n_added} triangles; "
              f"leftover front edges (holes): {leftover}")
    return n_added, leftover


def grow_patch(sa: Slice, sb: Slice, oracle, builder, n_probes: int = 12,
               extra_pts=(), point_slices=(), verbose=True):
    """Critical-interval wrapper: probes linearly between the flanks.

    extra_pts typically carries the CRITICAL POINTS (strand crossings of
    the critical fiber): bridging triangles across a saddle can never be
    Delaunay-empty amid both strands' dense corner samples — a vertex AT
    the crossing turns the local fan into clean small triangles.
    """
    probe_Es = [sa.E + (sb.E - sa.E) * k / (n_probes + 1)
                for k in range(1, n_probes + 1)]
    return grow_region([sa, sb], probe_Es, builder, oracle,
                       extra_pts=extra_pts, point_slices=point_slices,
                       verbose=verbose)


def grow_cap(edge_slice: Slice, direction: int, oracle, builder,
             n_probes: int = 8, tol_edge: float = 1e-4, verbose=True):
    """Close the mesh at a spectrum edge (loop shrinks to an extremum).

    Bisects outward from the last slice until the fiber disappears —
    the index change (0,n)->(0,0) IS the spectrum boundary.  The cap
    point is the chordal centroid of the tightest surviving loop (the
    loop radius shrinks like sqrt(E_edge - E), so at the bisection
    tolerance it is already sub-sampling); growth then closes the cone.
    """
    E0 = edge_slice.E
    lo, hi = E0, E0 + direction * 0.3
    tightest = edge_slice
    for _ in range(40):
        Em = 0.5 * (lo + hi)
        sm = oracle(Em)
        has = (not sm.is_critical) and bool(sm.loops)
        if has:
            tightest = sm
            lo = Em
        else:
            hi = Em
        if abs(hi - lo) < tol_edge:
            break
    E_edge = 0.5 * (lo + hi)
    builder.register_slice(tightest)
    pts = np.concatenate([l.chordal_r3() for l in tightest.loops])

    # cap point = centroid of the tightest loop on the Riemann spheres,
    # stereographically inverted back to beta coordinates
    mean_r3 = pts.mean(axis=0)

    def inv_stereo(v):
        x, y, z = v
        denom = 1.0 - z
        if abs(denom) < 1e-12:
            return None
        return complex(x / denom, y / denom)

    b1 = inv_stereo(mean_r3[:3])
    b2 = inv_stereo(mean_r3[3:])
    if verbose:
        print(f"    [cap] edge at {E_edge:+.5f} (from {E0:+.2f}, "
              f"dir {direction:+d}); cap beta1={b1:.4f} beta2={b2:.4f}")
    # probes clustered toward the edge (t = 1-(1-u)^2): the loops shrink
    # like sqrt(E_edge - E), so the last stretch needs the densest
    # sampling or the final fan to the cap point spans a large gap
    probe_Es = [E0 + (E_edge - E0) * (1.0 - (1.0 - k / (n_probes + 1)) ** 2)
                for k in range(1, n_probes + 1)]
    return grow_region([edge_slice], probe_Es, builder, oracle,
                       extra_pts=[(E_edge, b1, b2)], verbose=verbose)


def glue_critical(sa: Slice, sb: Slice, sc: Slice, builder):
    """Glue flanking clean slices through the critical fiber sc.

    Universal fan rule, per critical point P = a strand crossing of sc:
    each flank loop's passage through P's neighborhood contributes two
    ends; the 2 x 2 ends are ordered by chart angle around P and P is
    fanned over consecutive ends (Morse saddle -> alternating disk;
    other degeneracies -> same rule with more/fewer ends).  The loop
    remainders between passages are stitched pairwise by FKU, paired by
    Chamfer distance.
    """
    crossings = find_crossings(sc)
    if not crossings:
        touches = self_touches(sc)
        th = np.vstack([np.column_stack([l.th1, l.th2]) for l in sc.loops])
        crossings = [0.5 * (th[a] + th[b]) % (2 * np.pi)
                     for a, b, _ in touches]
    if not crossings:
        print(f"    !! no critical structure in [{sa.E:+.5f},{sb.E:+.5f}]"
              f" — plain FKU fallback")
        for i, j in match_components(sa.loops, sb.loops):
            strip_triangles(sa.loops[i], sb.loops[j], builder,
                            (sa.E, i), (sb.E, j))
        return 0

    i, j = match_components(sa.loops, sb.loops)[0]
    lo, hi = sa.loops[i], sb.loops[j]

    cuts = {"lo": [], "hi": []}
    wedge_geo = []   # per crossing: (P_th, sorted ends, their angles)
    for c, P in enumerate(crossings):
        P_th = (float(P[0]), float(P[1]))
        P_b1 = np.exp(1j * (P_th[0] % (2 * np.pi)))
        P_b2 = np.exp(1j * (P_th[1] % (2 * np.pi)))
        cp = builder.add_vertex(0.5 * (sa.E + sb.E), P_b1, P_b2,
                                note=f"crit E~{sc.E:+.5f}")

        ends = []
        for side, loop in (("lo", lo), ("hi", hi)):
            pe = _passage_ends(loop, P_th)
            if pe is None:
                k0 = int(np.argmin(np.abs(loop.th1 - P_th[0])
                                   + np.abs(loop.th2 - P_th[1])))
                pe = (k0, k0)
            ends.append((side, pe[0]))
            ends.append((side, pe[1]))
            cuts[side].append(pe)

        def angle(end):
            side, k = end
            loop = lo if side == "lo" else hi
            a1 = np.angle(np.exp(1j * (loop.th1[k] - P_th[0])))
            a2 = np.angle(np.exp(1j * (loop.th2[k] - P_th[1])))
            return float(np.arctan2(a2, a1))

        ends.sort(key=angle)
        wedge_geo.append((P_th, ends, [angle(e) for e in ends]))
        gids = [builder.gid(sa.E if s == "lo" else sb.E,
                            i if s == "lo" else j, k)
                for s, k in ends]
        for u, v in zip(gids, gids[1:] + gids[:1]):
            builder.tris.append((cp, u, v))

    def chart_angle(loop, k, P_th):
        a1 = np.angle(np.exp(1j * (loop.th1[k] - P_th[0])))
        a2 = np.angle(np.exp(1j * (loop.th2[k] - P_th[1])))
        return float(np.arctan2(a2, a1))

    def corridor_partner(cross, lo_end_k, probe_k):
        """The hi end sharing the wedge that the lo arc leaves INTO.

        probe_k: a loop vertex a step along the arc from lo_end_k — its
        angle around P disambiguates which of the two adjacent wedges is
        the corridor this arc travels (the other wedge belongs to the
        other corridor).
        """
        P_th, ends, angs = wedge_geo[cross]
        a = next(idx for idx, e in enumerate(ends)
                 if e == ("lo", lo_end_k))
        pa = angs[a]
        pp = chart_angle(lo, probe_k, P_th)
        best = None
        for nb in ((a - 1) % len(ends), (a + 1) % len(ends)):
            if ends[nb][0] != "hi":
                continue
            # angular interval from pa to angs[nb] (ccw), does it hold pp?
            lo_a, hi_a = (pa, angs[nb]) if ((angs[nb] - pa) % (2 * np.pi)
                                            < np.pi) else (angs[nb], pa)
            inside = ((pp - lo_a) % (2 * np.pi)) <= ((hi_a - lo_a)
                                                     % (2 * np.pi) + 1e-12)
            if inside:
                best = ends[nb]
        return best

    def remainder_arcs(loop, cut_list):
        arcs = []
        for a in range(len(cut_list)):
            s = cut_list[a][1]
            e = cut_list[(a + 1) % len(cut_list)][0]
            _, idx = loop.cyclic_arc(s, e)
            arcs.append(idx)
        return arcs

    arcs_lo = remainder_arcs(lo, cuts["lo"])
    arcs_hi = remainder_arcs(hi, cuts["hi"])
    n_cross = len(crossings)
    n = lo.n
    for k in range(n_cross):
        exit_k = cuts["lo"][k][1]
        k2 = (k + 1) % n_cross
        entry_k = cuts["lo"][k2][0]
        # probe one vertex ALONG the arc from each end
        probe_exit = (exit_k + 1) % n
        probe_entry = (entry_k - 1) % n
        m1 = corridor_partner(k, exit_k, probe_exit)
        m2 = corridor_partner(k2, entry_k, probe_entry)
        # identify the hi arc whose ends are (m1 at its crossing, m2 at
        # the other) — hi arc kk runs crossing kk (exit) -> kk+1 (entry)
        cand = None
        for kk in range(n_cross):
            ends_kk = {("hi", cuts["hi"][kk][1]), ("hi", cuts["hi"][k2 if
                        kk == k else (kk + 1) % n_cross][0])}
            if (("hi", cuts["hi"][kk][1]) == (m1 if m1 else None) or
                    m1 is None):
                pass
            e1 = m1 if m1 is not None else ("hi", cuts["hi"][kk][1])
            e2 = m2 if m2 is not None else ("hi", cuts["hi"][k2][0])
            if {e1, e2} <= {("hi", cuts["hi"][kk][1]),
                            ("hi", cuts["hi"][(kk + 1) % n_cross][0])}:
                cand = kk
                break
        kk = cand if cand is not None else k
        if len(arcs_lo[k]) < 2 or len(arcs_hi[kk]) < 2:
            continue
        arc_strip(lo, arcs_lo[k], hi, arcs_hi[kk], builder,
                  (sa.E, i), (sb.E, j))
    return len(crossings)


# ---------------------------------------------------------------------------
# Adaptive driver
# ---------------------------------------------------------------------------

def mesh_interval(sa, sb, oracle, builder, tol_E=1e-3, depth=0,
                  max_depth=18, log=None, far_count=None,
                  max_step: float = 0.1):
    """Stitch [sa.E, sb.E]; both slices must be clean (non-critical).

    A certificate PASS is only honoured when |dE| <= max_step: the
    thresholds were calibrated on adjacent grid slices, and after
    critical-value narrowing the recursion can pair qualitatively
    different loops (sharp-cornered near-diamond vs round) across a
    large gap — FKU then produced diametric fans (measured: 186 edges
    above 2.0 chordal, a +0.1 vertex fanned onto a +0.0016 corner).
    """
    ok, diag = certificate(sa, sb)
    if ok and abs(sb.E - sa.E) > max_step * (1.0 + 1e-9):
        # NB the tolerance: the grid step IS max_step, and floating point
        # makes |dE| land a hair above it on half the pairs — a strict >
        # rejects the ENTIRE bulk (measured: every chain pair logged
        # "gap 0.100 > max_step" and fell into bisection cascades)
        ok = False
        diag["reason"] = f"gap {abs(sb.E - sa.E):.3f} > max_step"
    if log is not None:
        log.append(("  " * depth)
                   + f"[{sa.E:+.5f},{sb.E:+.5f}] {diag['reason']}")
    if ok:
        for i, j in match_components(sa.loops, sb.loops):
            strip_triangles(sa.loops[i], sb.loops[j], builder,
                            (sa.E, i), (sb.E, j))
        return
    if abs(sb.E - sa.E) < tol_E or depth >= max_depth:
        # flanks tight, certificate still failing: greedy growth patch
        grow_patch(sa, sb, oracle, builder)
        if log is not None:
            log.append(("  " * depth)
                       + f"  -> grow critical patch "
                         f"[{sa.E:+.5f},{sb.E:+.5f}]")
        return

    Em = 0.5 * (sa.E + sb.E)
    sm = oracle(Em)
    if is_critical_slice(sm, far_count=far_count):
        # Em sits on/near the critical value.  The halving ladder toward
        # E* IS the FKU ladder: every intermediate probe is kept, oriented
        # and stitched to its neighbour (graded shape change per step),
        # so nothing ever pairs a near-diamond with a round loop across a
        # big gap.  Discarding the intermediate probes and recursing on
        # [sa, left] directly produced 186 diametric fan edges (2+ chordal).
        delta = 0.5 * min(Em - sa.E, sb.E - Em)
        ladder_lo, ladder_hi = [sa], [sb]
        left = right = None
        # the narrowing floor is tol_narrow (NOT tol_E): probes stay
        # clean down to |E| ~ 1e-4 (verified: single loop, no self-touch,
        # no degenerate nexus), and the corridor bridge gap shrinks like
        # sqrt(delta) — at +-1e-4 it is ~0.005 << line spacing 0.025,
        # i.e. isotropic sampling, and the bridging triangles become
        # Delaunay-clean by construction.  tol_E=1e-3 once stopped the
        # ladder at +-1.6e-3 where the gap (0.078 ~ 3x spacing) made
        # every bridge fail the emptiness test.
        tol_narrow = 1e-4
        for _ in range(20):
            l = oracle(Em - delta)
            r = oracle(Em + delta)
            if is_critical_slice(l, far_count) or \
                    is_critical_slice(r, far_count) or delta <= tol_narrow:
                break
            left, right = l, r
            builder.register_slice(l)
            builder.register_slice(r)
            ladder_lo.append(l)
            ladder_hi.append(r)
            delta *= 0.5
        if left is None:
            left = oracle(Em - delta)
            right = oracle(Em + delta)
            builder.register_slice(left)
            builder.register_slice(right)
            ladder_lo.append(left)
            ladder_hi.append(right)

        def stitch_ladder(ladder):
            # The ladder lies entirely on ONE side of E* and is graded by
            # construction (halving steps) — no critical value can sit
            # between rungs, so the certificate is unnecessary here.
            # Running it anyway misfires: near a sharp corner the sampling
            # is so dense that the max/median displacement RATIO explodes
            # at every scale, cascading sub-bisections into thousands of
            # pointless growth patches.
            for k in range(len(ladder) - 1):
                a = ladder[k]
                b = ladder[k + 1]
                for i in range(len(b.loops)):
                    b.loops[i] = orient_pair(a.loops[0], b.loops[i])
                for i, j in match_components(a.loops, b.loops):
                    strip_triangles(a.loops[i], b.loops[j], builder,
                                    (a.E, i), (b.E, j))

        stitch_ladder(ladder_lo)
        stitch_ladder(ladder_hi[::-1])   # from sb down to `right`

        # critical points from the critical fiber's strand crossings —
        # inserted into the growth cloud (see grow_patch docstring)
        crit_pts = []
        if sm.loops:
            for P in find_crossings(sm):
                crit_pts.append((sm.E,
                                 np.exp(1j * (float(P[0]) % (2 * np.pi))),
                                 np.exp(1j * (float(P[1]) % (2 * np.pi)))))
        # continue halving INSIDE the pinch zone: these slices cannot
        # be FKU boundaries (self-touching) but their points are exactly
        # the near-critical density the bridging needs — direct
        # flank-to-flank connections are allowed everywhere in growth
        point_slices = []
        d2 = delta
        while d2 > 1e-4:
            d2 *= 0.5
            for s in (oracle(Em - d2), oracle(Em + d2)):
                if (not s.is_critical) and s.loops:
                    point_slices.append(s)
        grow_patch(left, right, oracle, builder, n_probes=0,
                   extra_pts=crit_pts, point_slices=[])
        if log is not None:
            log.append(("  " * depth)
                       + f"[{sa.E:+.5f},{sb.E:+.5f}] critical value "
                           f"~{Em:+.5f}: ladder({len(ladder_lo)}+"
                           f"{len(ladder_hi)}) + grow patch at "
                           f"[{left.E:+.5f},{right.E:+.5f}] "
                           f"({len(crit_pts)} crit pts)")
        return
    builder.register_slice(sm)
    mesh_interval(sa, sm, oracle, builder, tol_E, depth + 1, max_depth,
                  log, far_count)
    mesh_interval(sm, sb, oracle, builder, tol_E, depth + 1, max_depth,
                  log, far_count)


# ---------------------------------------------------------------------------
# Hole filler (guarded mesh repair)
# ---------------------------------------------------------------------------

def _boundary_cycles(builder):
    """Walk boundary cycles of the mesh (undirected boundary graph).

    Every boundary vertex of a manifold-with-boundary has boundary
    degree 2; walks chain unambiguously.  (A directed walker once
    mis-chained at slit corners and silently lost the corridor cycles.)
    """
    from collections import defaultdict
    inc = defaultdict(int)
    for t in builder.tris:
        if len(set(t)) < 3:
            continue
        for k in range(3):
            a, b = t[k], t[(k + 1) % 3]
            inc[(min(a, b), max(a, b))] += 1
    adj = defaultdict(set)
    for (a, b), c in inc.items():
        if c == 1:
            adj[a].add(b)
            adj[b].add(a)
    seen = set()
    cycles = []
    for a, b in sorted(inc):
        if inc[(a, b)] != 1 or (a, b) in seen:
            continue
        cyc = [a]
        visited = {a}
        prev, cur = None, a
        # visited-set + step cap: at figure-eight junctions (boundary
        # degree >= 3) a plain prev-only walk ping-pongs between lobes
        # forever (this once hung the whole pipeline)
        for _step in range(4 * len(adj) + 8):
            nxts = [x for x in adj[cur]
                    if x != prev and x not in visited]
            if not nxts:
                nxts = [x for x in adj[cur]
                        if x != prev and x == a]
            if not nxts:
                break
            prev, cur = cur, nxts[0]
            if cur == a:
                break
            visited.add(cur)
            cyc.append(cur)
        for k in range(len(cyc)):
            seen.add((min(cyc[k], cyc[(k + 1) % len(cyc)]),
                      max(cyc[k], cyc[(k + 1) % len(cyc)])))
        if len(cyc) >= 3:
            cycles.append(cyc)
    return cycles


def fill_holes(builder, max_cycle: int = 32, width_factor: float = 1.2,
               verbose=True):
    """Close small boundary cycles by ear clipping in a local frame.

    Guards (a hole failing ANY guard is left open and reported loudly):
      - cycle length <= max_cycle          (big holes are structural)
      - slit WIDTH <= width_factor * (median cycle-edge length): the
        width of vertex v = distance to the nearest cycle vertex more
        than 2 steps away — a genuinely missing region (e.g. wrong
        pairing) is WIDE and must not be paved over
      - every produced triangle: positive projected area, no edge longer
        than 2x the longest cycle edge, E-span within the cycle's own
        E-span (wrong-corridor bridges would violate this)
    Orientation: the fill fan is flipped as a whole so that its first
    boundary edge runs opposite to the existing adjacent face.
    """
    from collections import defaultdict

    def r3(g):
        return _r3_of_coord(builder.coords[g])

    cycles = _boundary_cycles(builder)
    n_filled, n_skipped = 0, 0
    for cyc in cycles:
        n = len(cyc)
        if n < 3 or n > max_cycle:
            if verbose and n >= 3:
                print(f"    [fill] SKIP cycle of {n} edges (> {max_cycle})")
            n_skipped += 1
            continue
        pts = np.array([r3(g) for g in cyc])
        edges = [float(np.linalg.norm(pts[i] - pts[(i + 1) % n]))
                 for i in range(n)]
        spacing = float(np.median(edges))
        # slit width per vertex; cycles too small to have "far"
        # vertices use the max edge as the width proxy
        if n > 5:
            widths = []
            for i in range(n):
                d = [float(np.linalg.norm(pts[i] - pts[j]))
                     for j in range(n)
                     if min(abs(i - j), n - abs(i - j)) > 2]
                widths.append(min(d) if d else np.inf)
            width_stat = float(np.median(widths))
        else:
            widths = [max(edges)] * n
            width_stat = max(edges)
        if width_stat > width_factor * spacing:
            if verbose:
                print(f"    [fill] SKIP {n}-cycle: width "
                      f"{width_stat:.3f} > "
                      f"{width_factor}*{spacing:.3f}")
            n_skipped += 1
            continue

        # local frame: PCA of cycle + adjacent-face third vertices
        adj3 = []
        inc = defaultdict(list)
        for t in builder.tris:
            if len(set(t)) < 3:
                continue
            for k in range(3):
                a, b = t[k], t[(k + 1) % 3]
                inc[(min(a, b), max(a, b))].append(t)
        extra = []
        for i in range(n):
            e = (min(cyc[i], cyc[(i + 1) % n]),
                 max(cyc[i], cyc[(i + 1) % n]))
            for t in inc.get(e, []):
                v3 = next(v for v in t if v not in e)
                extra.append(r3(v3))
        cloud = np.vstack([pts, np.array(extra)]) if extra else pts
        center = cloud.mean(axis=0)
        _, _, vt = np.linalg.svd(cloud - center, full_matrices=False)
        e1, e2 = vt[0], vt[1]

        def proj(p):
            q = p - center
            return np.array([q @ e1, q @ e2])

        P = np.array([proj(p) for p in pts])

        def area2(a, b, c):
            return 0.5 * ((b[0] - a[0]) * (c[1] - a[1])
                          - (c[0] - a[0]) * (b[1] - a[1]))

        # shoelace sign (the earlier consecutive-triple sum was not the
        # polygon area and mis-signed zigzag slit boundaries)
        poly_sign = np.sign(sum(P[i][0] * P[(i + 1) % n][1]
                                - P[(i + 1) % n][0] * P[i][1]
                                for i in range(n)))
        if poly_sign == 0:
            if verbose:
                print(f"    [fill] SKIP {n}-cycle: degenerate projection")
            n_skipped += 1
            continue

        # Two-chain FKU: a slit is two nearly-parallel chains joined at
        # hairpin ends.  Ear clipping fails on such ribbons (every corner
        # ear contains opposite-side vertices); splitting at the two
        # sharpest turns and FKU-stitching the chains is exactly the
        # problem the bulk solver already handles well.
        # zigzag detection: consecutive pairs on the SAME side of the
        # minor axis (runs along one side) vs all-crossing (zigzag)
        minor_sign = np.sign(P[:, 1] - P[:, 1].mean())
        same_side = sum(1 for i in range(n)
                        if minor_sign[i] == minor_sign[(i + 1) % n])
        zigzag = same_side < n // 4

        if n == 3:
            new_tris = [(cyc[0], cyc[1], cyc[2])]
        elif not zigzag:
            def interior_angle(i):
                a = P[(i - 1) % n]; b = P[i]; c = P[(i + 1) % n]
                u = b - a; v = c - b
                nu, nv = np.linalg.norm(u), np.linalg.norm(v)
                if nu < 1e-15 or nv < 1e-15:
                    return np.pi
                return float(np.arccos(np.clip(u @ v / (nu * nv),
                                               -1.0, 1.0)))
            angles = [interior_angle(i) for i in range(n)]
            ends = sorted(range(n), key=lambda i: angles[i])[:2]
            u0, v0 = sorted(ends)
            chainA = cyc[u0:v0 + 1]
            chainB = cyc[v0:] + cyc[:u0 + 1]
            if len(chainA) < 2 or len(chainB) < 2:
                if verbose:
                    print(f"    [fill] SKIP {n}-cycle: bad chain split")
                n_skipped += 1
                continue

            def chain_arrays(chain):
                th1 = np.array([np.angle(builder.coords[g][1])
                                for g in chain])
                b2 = np.array([builder.coords[g][2] for g in chain])
                mu1 = float(np.log(abs(builder.coords[chain[0]][1])))
                return LineSubset(E=builder.coords[chain[0]][0], mu1=mu1,
                                  theta1_arr=th1, beta2_arr=b2)

            La, Lb = chain_arrays(chainA), chain_arrays(chainB)
            fv, ftris, _info = fku_triangulate(La, Lb)
            m = len(chainA)
            flip_a = abs(fv[0][1] - La.theta1_arr[-1]) < \
                abs(fv[0][1] - La.theta1_arr[0])
            flip_b = abs(fv[m][1] - Lb.theta1_arr[-1]) < \
                abs(fv[m][1] - Lb.theta1_arr[0])
            ia = chainA[::-1] if flip_a else chainA
            ib = chainB[::-1] if flip_b else chainB
            new_tris = []
            for t in ftris:
                vs = []
                for x in t:
                    x = int(x)
                    vs.append(ia[x] if x < m else ib[x - m])
                new_tris.append(tuple(vs))
        else:
            # ZIGZAG slit: every boundary edge crosses the slit; ear
            # clipping on the perimeter (in its own order) is the right
            # tool — the two-chain split would produce garbage chains
            def interior_angle2(i):
                a = P[(i - 1) % n]; b = P[i]; c = P[(i + 1) % n]
                u = b - a; v = c - b
                nu, nv = np.linalg.norm(u), np.linalg.norm(v)
                if nu < 1e-15 or nv < 1e-15:
                    return np.pi
                return float(np.arccos(np.clip(u @ v / (nu * nv), -1, 1)))

            idx = list(range(n))
            new_tris = []
            while len(idx) > 3:
                m2 = len(idx)
                best = None
                for k in range(m2):
                    a, b, c = idx[(k - 1) % m2], idx[k], idx[(k + 1) % m2]
                    A, B, C = P[a], P[b], P[c]
                    if area2(A, B, C) * poly_sign <= 1e-14:
                        continue
                    if any(area2(A, B, P[j]) * poly_sign > 0
                           and area2(B, C, P[j]) * poly_sign > 0
                           and area2(C, A, P[j]) * poly_sign > 0
                           for j in idx if j not in (a, b, c)):
                        continue
                    ang = interior_angle2(b)
                    if best is None or ang < best[0]:
                        best = (ang, k, (a, b, c))
                if best is None:
                    new_tris = None
                    break
                _, k, (a, b, c) = best
                new_tris.append((cyc[a], cyc[b], cyc[c]))
                idx.pop(k)
            if new_tris is not None and len(idx) == 3:
                new_tris.append((cyc[idx[0]], cyc[idx[1]],
                                 cyc[idx[2]]))
            if not new_tris:
                if verbose:
                    print(f"    [fill] SKIP {n}-cycle: no valid ear "
                          f"(zigzag)")
                n_skipped += 1
                continue

        # post-fill assertions
        ok = True
        E_vals = [builder.coords[g][0] for g in cyc]
        e_span = max(E_vals) - min(E_vals) + 1e-9
        max_edge = max(edges)
        # incidence pre-check: a fill edge coinciding with an existing
        # 2-incident edge would make the mesh non-manifold (the corner
        # folds once produced 20 such edges)
        from collections import Counter as _C
        _inc = _C()
        for t in builder.tris:
            if len(set(t)) < 3:
                continue
            for k in range(3):
                a, b = t[k], t[(k + 1) % 3]
                _inc[(min(a, b), max(a, b))] += 1
        bad_inc = False
        for t in new_tris:
            for x, y in ((t[0], t[1]), (t[1], t[2]), (t[0], t[2])):
                if _inc.get((min(x, y), max(x, y)), 0) >= 2:
                    bad_inc = True
        if bad_inc:
            if verbose:
                print(f"    [fill] SKIP {n}-cycle: edge already "
                      f"2-incident")
            n_skipped += 1
            continue

        fail_reason = []
        for t in new_tris:
            tv = np.array([r3(g) for g in t])
            for x, y in ((0, 1), (1, 2), (0, 2)):
                L = float(np.linalg.norm(tv[x] - tv[y]))
                if L > 2.0 * max_edge:
                    ok = False
                    fail_reason.append(f"edge {L:.3f} > {2*max_edge:.3f}")
            if abs(builder.coords[t[0]][0] - builder.coords[t[1]][0]) > e_span:
                ok = False
                fail_reason.append(
                    f"dE {abs(builder.coords[t[0]][0] - builder.coords[t[1]][0]):.6f}"
                    f" > span {e_span:.6f}")
        if not ok:
            if verbose:
                print(f"    [fill] SKIP {n}-cycle: post-check: "
                      f"{fail_reason[:3]}")
            n_skipped += 1
            continue

        # orientation: first fill triangle vs the existing face on its
        # first boundary edge — flip the whole fan if parallel
        a, b, c = new_tris[0]
        e0 = (min(a, b), max(a, b))
        faces = inc.get(e0, [])
        flip = False
        if faces:
            t0 = faces[0]
            k = t0.index(a) if a in t0 else None
            if k is not None and t0[(k + 1) % 3] == b:
                flip = True          # existing face walks a->b like us
        if flip:
            new_tris = [(x, z, y) for x, y, z in new_tris]
        for t in new_tris:
            builder.tris.append(t)
        n_filled += 1
        if verbose:
            print(f"    [fill] closed {n}-cycle with "
                  f"{len(new_tris)} triangles "
                  f"(width/spacing={width_stat/spacing:.2f})")
    if verbose:
        print(f"    [fill] done: {n_filled} filled, {n_skipped} skipped")
    return n_filled, n_skipped


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate(builder):
    from collections import Counter
    verts = builder.vertex_array()
    inc = Counter()
    for t in builder.tris:
        if len(set(t)) < 3:
            continue
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[0], t[2])):
            inc[(min(a, b), max(a, b))] += 1
    n_bnd = sum(1 for c in inc.values() if c == 1)
    n_dang = sum(1 for c in inc.values() if c == 0)
    n_nm = sum(1 for c in inc.values() if c > 2)

    def r3(i):
        E, b1, b2 = builder.coords[i]
        return np.concatenate([to_sphere_r3(np.array([b1]))[0],
                               to_sphere_r3(np.array([b2]))[0]])

    lens = []
    for (a, b) in inc:
        pa, pb = r3(a), r3(b)
        lens.append(float(np.linalg.norm(pa[:3] - pb[:3])
                          + np.linalg.norm(pa[3:] - pb[3:])))
    lens = np.array(lens)

    # chart validity: signed area of each triangle in (theta1, theta2)
    th = np.array([[np.angle(b1), np.angle(b2)]
                   for _, b1, b2 in builder.coords])
    n_flip = 0
    for t in builder.tris:
        if len(set(t)) < 3:
            continue
        p = th[list(t)]
        q = p - p[0]
        q = np.angle(np.exp(1j * q))       # local unwrap on the torus
        area = 0.5 * (q[1, 0] * q[2, 1] - q[2, 0] * q[1, 1])
        if area < 0:
            n_flip += 1
    return {
        "n_vertices": len(verts), "n_edges": len(inc),
        "n_faces": len(builder.tris),
        "euler": len(verts) - len(inc) + len(builder.tris),
        "n_boundary_edges": n_bnd, "n_dangling": n_dang,
        "n_nonmanifold": n_nm, "n_flipped_tris": n_flip,
        "edge_len_median": float(np.median(lens)) if len(lens) else None,
        "edge_len_max": float(lens.max()) if len(lens) else None,
        "n_critical_points": len(builder.crit_info),
    }


