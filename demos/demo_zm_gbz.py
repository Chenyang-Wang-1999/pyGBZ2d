'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-23
Copyright © Department of Physics, Tsinghua University. All rights reserved

Prototype: extract GBZ subsets directly from continuation.ZeroManager output.

Minimal test for replacing the zero-solving layer of brute_force_amoeba /
brute_force_SGBZ with ZeroManager.  Implements the refined design in
TODO/zm-replacement-spec.md.

Validation: 2D Hatano-Nelson model, analytic GBZ at mu1=gamma_1, mu2=gamma_2.
At mu1=gamma_1 the beta1 part collapses to 2cos(theta1), so the continuum
(|b2|=exp(gamma_2) is a root) holds for theta1 with cos(theta1) >= E/2-1.
'''

import sys
import warnings
from pathlib import Path
from cmath import exp, pi
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gbz_types import CharPoly, PointSubset, LineSubset
from continuation import ZeroManager

# A genuine continuum is constant-modulus to ~1e-9; a transversal crossing
# leaves the level after one sample.  1e-6 separates them.
CONTINUUM_TOL = 1e-6
CONTINUUM_FRAC = 0.9
# A crossing whose θ₁ lands within this of a continuum LineSubset endpoint is
# treated as the continuum/MR boundary, not a genuine discrete zero, and
# dropped.  Decoupled from CONTINUUM_TOL so the snap radius is independent of
# band-detection sensitivity.
SNAP_TOL = 1e-3


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

def build_HN2D(J1, J2, gamma_1, gamma_2, delta_1, delta_2):
    J11 = exp(gamma_1 + 1j * delta_1) * J1
    J12 = exp(-gamma_1 + 1j * delta_1) * np.conj(J1)
    J21 = exp(gamma_2 + 1j * delta_2) * J2
    J22 = exp(-gamma_2 + 1j * delta_2) * np.conj(J2)
    coeffs = np.array([1, -J11, -J12, -J21, -J22], dtype=complex)
    degs = np.array([
        [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1],
    ], dtype=int)
    return coeffs, degs


# ---------------------------------------------------------------------------
# Amoeba extraction  (spec 3.1)
# ---------------------------------------------------------------------------

def extract_amoeba_subsets(zm, poly, E, mu1, mu2, tol=CONTINUUM_TOL,
                           frac=CONTINUUM_FRAC, snap_tol=SNAP_TOL):
    """Zeros of f at ln|b1|=mu1, ln|b2|=mu2, from ZeroManager tracks."""
    # Crossings are recorded as (seg_idx, i, j, kind), NOT PointSubsets, so
    # that boundary duplicates between segments sharing an endpoint (an MR)
    # can be resolved after the loop.  kind:
    #   'zero'  — d == 0 at sample i  (exact numerical touch)
    #   'cross' — sign change over [i, i+1]
    hits: list[tuple[int, int, int, str]] = []
    continuum_lines: list[LineSubset] = []
    # Per-segment logabs cache so the post-loop materializer can recompute
    # points without re-taking logs.  Indexed by the same seg_idx (0-based)
    # used in hits.
    seg_logabs: list[np.ndarray] = []

    for seg_idx, seg in enumerate(zm.segments):
        th = seg.theta1_arr
        tr = seg.tracked_roots
        logabs = np.log(np.abs(tr))                       # (N, K)
        seg_logabs.append(logabs)

        # 1. continuum tracks (fraction-based).
        #    A genuine continuum is a whole constant-modulus row, so the
        #    fraction test is applied across all tracks at once and each hit
        #    becomes one LineSubset over the full segment.
        frac_in_band = np.mean(np.abs(logabs - mu2) < tol, axis=0)   # (K,)
        is_cont = frac_in_band > frac
        for j in np.where(is_cont)[0]:
            if frac_in_band[j] < 1.0:
                warnings.warn(
                    f"amoeba continuum track {j}: frac={frac_in_band[j]:.3f} "
                    f"< 1 (possible numerical outlier)")
            continuum_lines.append(LineSubset(
                E=E, mu1=mu1,
                theta1_arr=th.copy(), beta2_arr=tr[:, j].copy()))

        # 2. discrete crossings on the non-continuum tracks.
        #    Drop the continuum rows; detect (a) exact d==0 touches and (b)
        #    sign-change crossings of ln|b2| - mu2, both vectorized across
        #    the remaining tracks.  Only indices are stored here.
        other = np.where(~is_cont)[0]
        if other.size:
            d = logabs[:, other] - mu2                     # (N, J)

            # d == 0 cannot trigger the sign-change product below (0 * x = 0,
            # never < 0), so record it separately.
            i_idx, j_idx = np.where(d == 0)
            for i, jj in zip(i_idx, j_idx):
                hits.append((seg_idx, int(i), int(other[jj]), 'zero'))

            sc = d[:-1] * d[1:] < 0                        # (N-1, J)
            i_idx, j_idx = np.where(sc)
            for i, jj in zip(i_idx, j_idx):
                hits.append((seg_idx, int(i), int(other[jj]), 'cross'))

    # ---- materialize + boundary filtering ----
    # Two boundary rules (segments share endpoints at an MR, so a crossing at
    # a shared sample would otherwise be produced twice):
    #   Rule 1 — drop any crossing whose θ₁ is within snap_tol of a continuum
    #            LineSubset endpoint.  These sit on the continuum/MR edge, not
    #            a genuine discrete zero.
    #   Rule 2 — for 'zero' hits only, drop cross-segment duplicates arising
    #            from the shared boundary sample (same θ₁).  Interior 'zero'
    #            hits (unique θ₁) are kept once.
    cont_endpoints = np.array(
        [t for L in continuum_lines for t in (L.theta1_start, L.theta1_end)]
    ) if continuum_lines else None

    seen_zero_theta1: set[int] = set()
    points: list[tuple[float, complex]] = []
    for seg_idx, i, j, kind in hits:
        seg = zm.segments[seg_idx]
        th, tr = seg.theta1_arr, seg.tracked_roots
        la = seg_logabs[seg_idx]
        if kind == 'zero':
            t1, b2 = float(th[i]), complex(tr[i, j])
        else:
            frac_i = (mu2 - la[i, j]) / (la[i + 1, j] - la[i, j])
            t1 = float(th[i] + frac_i * (th[i + 1] - th[i]))
            b2 = complex(tr[i, j] + frac_i * (tr[i + 1, j] - tr[i, j]))

        # Rule 1: continuum-endpoint snap.
        if cont_endpoints is not None and np.min(np.abs(cont_endpoints - t1)) < snap_tol:
            continue

        # Rule 2: 'zero' boundary-duplicate dedup (shared endpoint → same θ₁).
        if kind == 'zero':
            key = round(t1, 12)
            if key in seen_zero_theta1:
                continue
            seen_zero_theta1.add(key)

        points.append((t1, b2))

    subsets: list = list(continuum_lines)
    for t1, b2 in points:
        subsets.append(PointSubset(E=E, beta1=exp(mu1 + 1j * t1), beta2=b2))
    return subsets


# ---------------------------------------------------------------------------
# SGBZ extraction  (spec 4.1)
# ---------------------------------------------------------------------------

def _runs(mask):
    """Contiguous True runs -> list of (start, end_exclusive)."""
    runs = []
    i, n = 0, len(mask)
    while i < n:
        if mask[i]:
            j = i
            while j < n and mask[j]:
                j += 1
            runs.append((i, j))
            i = j
        else:
            i += 1
    return runs


def extract_sgbz_subsets(zm, poly, E, mu1, tol=CONTINUUM_TOL,
                         frac=CONTINUUM_FRAC):
    """Zeros of f at ln|b2| = mid-gap, from ZeroManager tracks.

    Cluster-aware: 3+ tracks can share one modulus.  Continuum = an
    equal-modulus CLUSTER straddling the M-boundary; each cluster member
    -> one LineSubset.  Accidental = a track crossing the mid-gap level
    mid(theta1) = (m_{M-1}+m_M)/2 (sign change of ln|b2_j| - mid).
    """
    from scipy.optimize import minimize_scalar

    M = poly.M
    K = poly.M + poly.N
    subsets = []
    all_continuum_ranges: list[tuple[float, float]] = []

    for seg in zm.segments:
        th = seg.theta1_arr
        tr = seg.tracked_roots
        absmod = np.abs(tr)
        logabs = np.log(absmod)
        order = seg.abs_argsort
        N = len(th)
        mid_arr = 0.5 * (logabs[np.arange(N), order[:, M - 1]]
                         + logabs[np.arange(N), order[:, M]])

        # 1. equal-modulus clusters (union-find)
        parent = list(range(K))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i in range(K):
            for j in range(i + 1, K):
                f = float(np.mean(np.abs(absmod[:, i] - absmod[:, j]) < tol))
                if f > frac:
                    parent[find(i)] = find(j)
        clusters: dict[int, list[int]] = {}
        for i in range(K):
            clusters.setdefault(find(i), []).append(i)

        # 2. continuum: cluster straddling M-boundary -> LineSubsets
        continuum_members: set[int] = set()
        for members in clusters.values():
            s = len(members)
            if s < 2:
                continue
            c = np.mean(absmod[:, members], axis=1)
            others = np.array([k for k in range(K) if k not in members])
            n_below = (np.sum(absmod[:, others] < c[:, None], axis=1)
                       if others.size else np.zeros(N, dtype=int))
            in_boundary = (n_below <= M - 1) & (n_below + s - 1 >= M)
            for s_idx, e_idx in _runs(in_boundary):
                if e_idx - s_idx < 2:
                    continue
                t_start = float(th[max(0, s_idx - 1)])
                t_end = float(th[min(N - 1, e_idx)])
                for m in members:
                    subsets.append(LineSubset(
                        E=E, mu1=mu1,
                        theta1_arr=th[s_idx:e_idx].copy(),
                        beta2_arr=tr[s_idx:e_idx, m].copy()))
                continuum_members.update(members)
                all_continuum_ranges.append((t_start, t_end))

        # 3. accidental: per-track mid-gap crossing (refined)
        step0 = float(th[1] - th[0]) if N >= 2 else 0.0
        step1 = float(th[-1] - th[-2]) if N >= 2 else 0.0
        for j in range(K):
            if j in continuum_members:
                continue
            d = logabs[:, j] - mid_arr
            for i in range(N - 1):
                if d[i] * d[i + 1] >= 0:
                    continue
                fr = -d[i] / (d[i + 1] - d[i])
                t1_approx = float(th[i] + fr * (th[i + 1] - th[i]))
                if abs(t1_approx - th[0]) < step0 or abs(t1_approx - th[-1]) < step1:
                    continue
                if any(lo_t <= t1_approx <= hi_t
                       for lo_t, hi_t in all_continuum_ranges):
                    continue
                # refine via gap minimization
                t_left = float(th[i])
                t_right = float(th[i + 1])
                if t_right < t_left:
                    t_right += 2 * pi

                def _gap(t):
                    beta1 = exp(mu1 + 1j * t)
                    roots = np.asarray(poly.solve_roots_1d(
                        (0, 1), (E, beta1), (2,)), dtype=complex)
                    o = np.argsort(np.abs(roots))
                    return float(np.log(np.abs(roots[o[M]]))
                                 - np.log(np.abs(roots[o[M - 1]])))
                res = minimize_scalar(
                    _gap, bounds=(t_left, t_right), method='bounded')
                t1 = float(res.x % (2 * pi)) if res.success and res.fun < 1e-8 else t1_approx
                beta1 = exp(mu1 + 1j * t1)
                roots = np.asarray(poly.solve_roots_1d(
                    (0, 1), (E, beta1), (2,)), dtype=complex)
                o = np.argsort(np.abs(roots))
                for r in (roots[o[M - 1]], roots[o[M]]):
                    subsets.append(PointSubset(
                        E=E, beta1=beta1, beta2=complex(r)))

    # 4. MR as candidate GBZ point (isolated MRs only)
    for mr in zm.multiple_roots:
        th_mr = mr.theta1 % (2 * pi)
        near_cont = any(
            min(abs((th_mr - lo) % (2 * pi)),
                abs((th_mr - lo + pi) % (2 * pi)) - pi) < 0.05
            or min(abs((th_mr - hi) % (2 * pi)),
                   abs((th_mr - hi + pi) % (2 * pi)) - pi) < 0.05
            for lo, hi in all_continuum_ranges
        )
        if near_cont:
            continue
        for tup in mr.cluster_indices:
            s = len(tup)
            start = min(tup)
            if start <= M - 1 and start + s - 1 >= M:
                for idx in tup:
                    subsets.append(PointSubset(
                        E=E, beta1=exp(mu1 + 1j * mr.theta1),
                        beta2=complex(mr.roots[idx])))
    return subsets


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def _fmt_th(t):
    return f"{float(t) % (2 * pi):.4f}"


def main():
    coeffs, degs = build_HN2D(1, 1, 0.2, 0.3, 0, 0)
    poly = CharPoly(coeffs, degs)
    mu1, mu2 = 0.2, 0.3

    print("=" * 70)
    print("2D HN model: gamma=(0.2,0.3), analytic GBZ at mu1=0.2, mu2=0.3")
    print("Continuum (|b2|=exp(0.3) is a root) holds when cos(t1) >= E/2-1.")
    print("=" * 70)

    for E in [0.0, 1.0, 2.0]:
        cthr = E / 2 - 1
        if cthr <= -1:
            anal = "full circle"
        elif cthr >= 1:
            anal = "none"
        else:
            t_edge = float(np.arccos(cthr))
            anal = f"[0, {t_edge:.4f}] U [{2*pi - t_edge:.4f}, 2pi]"

        zm = ZeroManager(poly, E + 0j, mu1)
        zm.run(h0=0.1)

        amoeba = extract_amoeba_subsets(zm, poly, E + 0j, mu1, mu2)
        sgbz = extract_sgbz_subsets(zm, poly, E + 0j, mu1)

        print(f"\n--- E = {E}  (analytic continuum: {anal}) ---")
        print(f"  ZeroManager: {zm.n_segments} seg, {zm.n_multiple_roots} MR")
        print(f"  AMOEBA index={_idx(amoeba)}")
        for s in amoeba:
            _print_subset(s, "    ")
        print(f"  SGBZ   index={_idx(sgbz)}")
        for s in sgbz:
            _print_subset(s, "    ")


def cluster_test():
    """3 tracks sharing one modulus (size-3 cluster)."""
    coeffs = np.array([1, -1], dtype=complex)
    degs = np.array([[0, 0, 2], [0, 3, -1]], dtype=int)
    poly = CharPoly(coeffs, degs)
    mu1 = 0.0
    print("=" * 70)
    print(f"Cluster test: f = b2^2 - b1^3/b2  (3 roots at |b2|=exp({mu1}))")
    print("Expect 3 LineSubsets (one per cluster curve), NOT 6.")
    print("=" * 70)

    zm = ZeroManager(poly, 0.0 + 0j, mu1)
    zm.run(h0=0.1)
    print(f"  ZeroManager: {zm.n_segments} seg, {zm.n_multiple_roots} MR, K={zm.K}")

    amoeba = extract_amoeba_subsets(zm, poly, 0.0 + 0j, mu1, mu1)
    print(f"  AMOEBA index={_idx(amoeba)}  (expect (0, 3))")
    for s in amoeba:
        _print_subset(s, "    ")

    sgbz = extract_sgbz_subsets(zm, poly, 0.0 + 0j, mu1)
    print(f"  SGBZ   index={_idx(sgbz)}  (expect (0, 3))")
    for s in sgbz:
        _print_subset(s, "    ")


def accidental_test():
    """MR + accidental PMGBZ."""
    coeffs = np.array([1, -2, -1, 2], dtype=complex)
    degs = np.array([[0, 0, 1], [0, 0, 0], [0, 1, -1], [0, 0, -1]], dtype=int)
    poly = CharPoly(coeffs, degs)
    mu1 = 0.0
    print("=" * 70)
    print("Accidental test: f = b2 - 2 - b1/b2 + 2/b2")
    print("Expect MR@0 (2 pts, b2=1) + accidental@pi (2 pts, |b2|=sqrt3)")
    print("=" * 70)

    zm = ZeroManager(poly, 0.0 + 0j, mu1)
    zm.run(h0=0.1, cluster_tol=1e-4, min_dtheta=1e-6)
    print(f"  ZeroManager: {zm.n_segments} seg, {zm.n_multiple_roots} MR")
    sgbz = extract_sgbz_subsets(zm, poly, 0.0 + 0j, mu1)
    print(f"  SGBZ index={_idx(sgbz)}  (expect (4, 0))")
    for s in sgbz:
        _print_subset(s, "    ")


def _idx(subsets):
    n0 = sum(1 for s in subsets if isinstance(s, PointSubset))
    n1 = sum(1 for s in subsets if isinstance(s, LineSubset))
    return (n0, n1)


def _print_subset(s, indent=""):
    if isinstance(s, PointSubset):
        print(f"{indent}Point: th1={_fmt_th(np.angle(s.beta1))} "
              f"th2={_fmt_th(np.angle(s.beta2))} |b2|={abs(s.beta2):.4f}")
    else:
        print(f"{indent}Line:  th1=[{_fmt_th(s.theta1_start)}, "
              f"{_fmt_th(s.theta1_end)}] w={s.theta1_width:.4f}")


if __name__ == "__main__":
    main()
    cluster_test()
    accidental_test()
