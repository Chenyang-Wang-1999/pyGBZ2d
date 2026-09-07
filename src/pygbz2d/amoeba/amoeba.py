'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-05-19
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

from __future__ import annotations

from cmath import exp
from typing import Optional
import numpy as np

from pygbz2d import core
from pygbz2d.core import live_defaults

#: Plateau-probe non-zero-winding area threshold.
PLATEAU_AREA_THRESHOLD: float = 1e-2
#: Sample snap tolerance in circular θ₁ and chordal β₂ distance; both
#: coordinates must match before a discrete point is removed.
SNAP_TOL: float = 1e-3
from pygbz2d.core import (
    PointSubset, LineSubset, GBZResult, CharPoly,
    JoinableLinePiece, is_mr_cluster_endpoint, TWO_PI,
    check_points_clustered_on_torus, probe_zero_plateau,
    chordal_cost_matrix,
)

from .bisect import (
    bisect_amoeba_ronkin_min, _find_mu2_for_w2_zero,
)
from .ronkin_winding import (
    _get_average_winding_from_zeros,
)
from .zm_extract import find_crossings


# ---- plateau pre-check helpers ----

def _check_zeros_are_clustered(
    zeros: list[tuple[float, float, int]],
    tol_normalized: float,
) -> bool:
    """Whether every amoeba zero has a neighbour within *tol_normalized*.

    Thin amoeba adapter over
    :func:`pygbz2d.core.check_points_clustered_on_torus`: amoeba zeros carry a
    trailing ``jump`` component, so only the first two are forwarded as
    ``(θ₁, θ₂)``.

    At a genuine GBZ point zeros are well-separated (they partition the
    circle into meaningful segments).  At a zero-plateau boundary the
    winding changes sign over a vanishingly narrow angular region, so
    the zeros cluster into nearly degenerate pairs — each zero sits
    within tol_normalized of a neighbour.
    """
    points = [(z[0], z[1]) for z in zeros]
    return check_points_clustered_on_torus(points, tol_normalized)


# ---- plateau detection ----

def _is_zero_plateau_probe(point: dict, winding_tol: float) -> bool:
    return (
        point["success"]
        and (not point["is_continuum"])
        and point["zero_count"] == 0
        and abs(point["w1"]) <= winding_tol
    )


@live_defaults(continuum_tol="core:CONTINUUM_TOL", continuum_perturb="core:CONTINUUM_PERTURB",
    max_iter="amoeba.bisect:BISECT_MAX_ITER", xtol="amoeba.bisect:BISECT_XTOL",
    max_range_expansions="amoeba.bisect:MAX_RANGE_EXPANSIONS",
    range_expand_factor="amoeba.bisect:RANGE_EXPAND_FACTOR")
def _probe_zero_plateau_near_mu1(
    char_poly: CharPoly,
    E_ref: complex,
    mu1: float,
    mu1_bracket: Optional[tuple[float, float]],
    mu2_low: float = -1,
    mu2_high: float = 1,
    continuum_tol: Optional[float] = None,
    continuum_perturb: Optional[float] = None,
    max_iter: Optional[int] = None,
    xtol: Optional[float] = None,
    max_range_expansions: Optional[int] = None,
    range_expand_factor: Optional[float] = None,
    winding_tol: Optional[float] = None,
    probe_radius: Optional[float] = None,
) -> dict:
    """Check whether a nonempty-zero candidate sits next to a zero plateau.

    The decisive plateau signature is strict: after solving a2=0 at a nearby
    mu1, w1 is zero and there are no a2-crossing zeros.
    """
    if winding_tol is None:
        winding_tol = max(xtol, core.WINDING_ZERO_TOL)
    if probe_radius is None:
        probe_radius = continuum_perturb

    bracket_width = 0.0
    if mu1_bracket is not None:
        bracket_width = abs(float(mu1_bracket[1]) - float(mu1_bracket[0]))

    def evaluator(mu1_probe: float) -> dict:
        # The shared ladder protocol (pygbz2d.core.probe_zero_plateau) reads
        # "success" to skip failed probes and keep walking the ladder.  The
        # inner μ₂ bisection legitimately fails on the out-of-spectrum side
        # of a band edge (w2 has no sign change → range-expansion error) —
        # exactly where the probe operates.  Catch and report instead of
        # letting the raise escalate into a failed whole-energy GBZResult.
        try:
            inner = _find_mu2_for_w2_zero(
                char_poly, E_ref, mu1_probe, mu2_low, mu2_high,
                continuum_tol=continuum_tol,
                continuum_perturb=continuum_perturb, max_iter=max_iter,
                xtol=xtol, max_range_expansions=max_range_expansions,
                range_expand_factor=range_expand_factor,
            )
        except Exception as exc:
            return {
                "success": False,
                "error": f"{type(exc).__name__}: {exc}",
                "is_plateau": False,
                "is_continuum": False,
            }
        zeros_probe = inner.get("zeros") or []
        w1, _ = _get_average_winding_from_zeros(
            char_poly, E_ref, mu1_probe, inner["mu2"],
            zeros_probe, direction=1,
        )
        point = {
            "success": True,
            "mu2": inner["mu2"],
            "w1": float(w1),
            "zero_count": len(zeros_probe),
            "is_continuum": bool(inner["is_continuum"]),
        }
        # Collapse the amoeba-specific criterion to the one bool the shared
        # loop reads; keep the raw fields for diagnostics.
        point["is_plateau"] = _is_zero_plateau_probe(point, winding_tol)
        return point

    res = probe_zero_plateau(
        mu1, mu1_bracket,
        zero_tol=winding_tol,
        probe_radius=probe_radius,
        bracket_width=bracket_width,
        evaluator=evaluator,
    )
    # Back-compat alias: amoeba diagnostics historically read ``winding_tol``.
    res["winding_tol"] = winding_tol
    return res


# ---- GBZ subset assembly ----

def _assemble_discrete_subsets(
    E_ref: complex,
    mu1: float,
    mu2: float,
    zeros: list[tuple[float, float]],
) -> list[PointSubset]:
    """Turn discrete (θ₁, θ₂) zeros into PointSubsets."""
    points: list[PointSubset] = []
    for t1, t2 in zeros:
        points.append(PointSubset(
            E=E_ref,
            beta1=exp(mu1 + 1j * float(t1)),
            beta2=exp(mu2 + 1j * float(t2)),
        ))
    return points


def _splice_continuum_pieces(
    zm,
    pieces: list[JoinableLinePiece],
) -> list[JoinableLinePiece]:
    """Splice per-segment continuum pieces that meet end-to-end.

    Interior MR boundaries share one track frame; the θ=0≡2π seam is
    translated through ``zm.boundary_perm``.  A piece whose boundary track is
    in the MR cluster is a genuine terminator and is not spliced.
    """
    n_seg = len(zm.segments)
    if n_seg <= 1:
        return pieces

    def _mod(k: int) -> int:
        return k % n_seg

    def find_by_left(seg_s: int, root: complex) -> int | None:
        for idx, p in enumerate(pieces):
            if p.ml == seg_s and np.abs(p.beta2_arr[0] - root) < 1e-9:
                return idx
        return None

    def find_by_right(seg_s: int, root: complex) -> int | None:
        for idx, p in enumerate(pieces):
            if p.mr == seg_s and np.abs(p.beta2_arr[-1] - root) < 1e-9:
                return idx
        return None

    def _match_column_right(cur_s: int, prev_s: int, col: int) -> int:
        if cur_s == 0 and prev_s == n_seg - 1:
            return int(zm.boundary_perm[int(col)])
        return int(col)

    while True:
        changed = False
        for s in range(n_seg):
            seg = zm.segments[s]
            prev = _mod(s - 1)
            prev_seg = zm.segments[prev]

            left_mr = seg.left_mr
            right_mr = prev_seg.right_mr
            is_circle_seam = (left_mr < 0 and right_mr < 0)
            if not is_circle_seam and left_mr != right_mr:
                continue

            left_b = seg.tracked_roots[0, :]
            right_b = prev_seg.tracked_roots[-1, :]

            for j_l in range(zm.K):
                root = complex(left_b[j_l])
                if is_mr_cluster_endpoint(zm, seg, 'left', int(j_l)):
                    continue
                li_idx = find_by_left(s, root)
                if li_idx is None:
                    continue
                j_prev = _match_column_right(s, prev, int(j_l))
                root_prev = complex(right_b[j_prev])
                if is_mr_cluster_endpoint(zm, prev_seg, 'right', j_prev):
                    raise ValueError(
                        f"Continuum track {j_l} of segment {s} ends at the MR "
                        f"at θ₁={seg.theta1_arr[0]:.4f} as a non-cluster root, "
                        f"but the matched root (track {j_prev}) of segment "
                        f"{prev} is a cluster root there."
                    )
                pi_idx = find_by_right(prev, root_prev)
                if pi_idx is None or pi_idx == li_idx:
                    continue

                _splice_two(pieces, pi_idx, li_idx,
                            cyclic=(s == 0 and prev == n_seg - 1))
                changed = True
                break
            if changed:
                break
        if not changed:
            break
    return pieces


def _splice_two(
    pieces: list[JoinableLinePiece],
    prev_idx: int,
    cur_idx: int,
    *,
    cyclic: bool,
) -> None:
    """Merge pieces[prev_idx] (right side) with pieces[cur_idx] (left side)."""
    rp = pieces[prev_idx]
    cp = pieces[cur_idx]
    th = np.concatenate([rp.theta1_arr, cp.theta1_arr[1:]])
    b2 = np.concatenate([rp.beta2_arr, cp.beta2_arr[1:]])
    if cyclic:
        n_rp = len(rp.theta1_arr)
        th[n_rp:] += TWO_PI
    merged = JoinableLinePiece(
        E=cp.E, mu1=cp.mu1, theta1_arr=th, beta2_arr=b2,
        ml=rp.ml, mr=cp.mr,
    )
    keep = [i for i in range(len(pieces)) if i not in (prev_idx, cur_idx)]
    pieces[:] = [pieces[i] for i in keep] + [merged]


def _assemble_continuum_subsets(
    zm,
    E_ref: complex,
    mu1: float,
    mu2: float,
    amoeba_res: dict,
) -> list:
    """Assemble LineSubsets + PointSubsets for a continuum GBZ point.

    The continuum members from the bisection are materialised as LineSubset
    pieces and spliced end-to-end where they meet.  All other tracks are
    searched for discrete crossings of ``ln|β₂| = μ₂`` (with the continuum
    members passed as ``avoided_segments``) and converted to PointSubsets.
    """
    members = amoeba_res.get("_continuum_members") or []
    if not members:
        raise ValueError("No continuum members found, but is_continuum=True.")

    pieces: list[JoinableLinePiece] = []
    for s, j in members:
        seg = zm.segments[s]
        pieces.append(JoinableLinePiece(
            E=E_ref, mu1=mu1,
            theta1_arr=seg.theta1_arr.copy(),
            beta2_arr=seg.tracked_roots[:, int(j)].copy(),
            ml=s, mr=s,
        ))
    pieces = _splice_continuum_pieces(zm, pieces)

    subsets: list = [
        LineSubset(E=p.E, mu1=p.mu1,
                   theta1_arr=p.theta1_arr, beta2_arr=p.beta2_arr)
        for p in pieces
    ]
    line_subsets = subsets.copy()

    crossings = find_crossings(
        zm, mu1, mu2,
        avoided_segments=members,
        return_refined=True,
    )
    for b1, b2 in crossings:
        point = PointSubset(E=E_ref, beta1=b1, beta2=b2)
        if _point_near_any_line(point, line_subsets, SNAP_TOL):
            continue
        subsets.append(point)

    return subsets


def _point_near_any_line(
    point: PointSubset,
    lines: list[LineSubset],
    snap_tol: float,
) -> bool:
    """Match a point to a sampled line point in both coordinates.

    Assembly supplies points and lines at the same mu1. Only sampled
    coincidences are removed: interpolating across a sparsely sampled curve
    could erase a distinct root track near a multiple root.
    """
    t1 = point.theta1
    for line in lines:
        if abs(point.mu1 - line.mu1) >= snap_tol:
            continue
        d = np.abs((line.theta1_arr - t1 + np.pi) % (TWO_PI) - np.pi)
        nearby = line.beta2_arr[d < snap_tol]
        if len(nearby) and np.any(chordal_cost_matrix(
                np.array([point.beta2]), nearby) < snap_tol):
            return True
    return False


# ---- main entry point ----

@live_defaults(continuum_tol="core:CONTINUUM_TOL", continuum_perturb="core:CONTINUUM_PERTURB",
    max_iter="amoeba.bisect:BISECT_MAX_ITER", xtol="amoeba.bisect:BISECT_XTOL",
    max_range_expansions="amoeba.bisect:MAX_RANGE_EXPANSIONS",
    range_expand_factor="amoeba.bisect:RANGE_EXPAND_FACTOR",
    plateau_cluster_tol="core:PLATEAU_CLUSTER_TOL",
    plateau_area_threshold="amoeba.amoeba:PLATEAU_AREA_THRESHOLD")
def collect_GBZ_subsets(
    coeffs: np.ndarray,
    degs: np.ndarray,
    E_ref: complex,
    perc: float = None,
    debug_mode: bool = False,
    *,
    plateau_check: bool = True,
    plateau_winding_tol: Optional[float] = None,
    plateau_probe_radius: Optional[float] = None,
    plateau_area_threshold: Optional[float] = None,
    plateau_cluster_tol: Optional[float] = None,
    mu1_low: float = -1,
    mu1_high: float = 1,
    mu2_low: float = -1,
    mu2_high: float = 1,
    continuum_tol: Optional[float] = None,
    continuum_perturb: Optional[float] = None,
    max_iter: Optional[int] = None,
    xtol: Optional[float] = None,
    max_range_expansions: Optional[int] = None,
    range_expand_factor: Optional[float] = None,
) -> GBZResult:
    """Check amoeba condition and return GBZ points for a reference energy.

    The zero-solving layer is ``continuation.ZeroManager`` (wrapped as
    ``AmoebaZeroManager``).  The μ₁/μ₂ bisection builds one ZM per μ₁ and
    reuses it across all μ₂ evaluations; the solved ZM is then reused here
    for subset extraction (no second run).

    Returns:
        GBZResult with connected subsets.  ``gbz.is_empty`` means E_ref is
        outside the amoeba GBZ spectrum.
    """
    if perc is not None:
        print("%.2f" % (perc * 100) + r"%")
    char_poly = CharPoly(coeffs, degs)

    # Neighbour threshold note (plateau_cluster_tol vs
    # plateau_area_threshold): deliberately SEPARATE knobs — the former is
    # a torus-clustering radius, the latter a winding-area fraction.
    try:
        amoeba_res = bisect_amoeba_ronkin_min(
            char_poly, E_ref,
            mu1_low=mu1_low, mu1_high=mu1_high,
            mu2_low=mu2_low, mu2_high=mu2_high,
            continuum_tol=continuum_tol,
            continuum_perturb=continuum_perturb,
            max_iter=max_iter, xtol=xtol,
            max_range_expansions=max_range_expansions,
            range_expand_factor=range_expand_factor,
        )

        # Reuse the ZeroManager built inside the bisection at the solved
        # mu1 — subset assembly runs on the same adaptive tracks, no second
        # ZeroManager run.
        zm = amoeba_res["_zm"]
        mu1, mu2 = amoeba_res["mu1"], amoeba_res["mu2"]

        # Determine if this is actually a zero plateau (no GBZ)
        is_amoeba = True
        if (not amoeba_res["is_continuum"]) and (len(amoeba_res["zeros"]) == 0):
            is_amoeba = False
        elif plateau_check and (not amoeba_res["is_continuum"]):
            # Lightweight plateau pre-check: at a zero-plateau boundary
            # three necessary conditions hold simultaneously:
            #
            #   (a) w1 & w2 non-zero winding areas are tiny — most of
            #       the angular circle has u = 0.
            #   (b) zeros cluster into nearly degenerate pairs — at a
            #       genuine GBZ point zeros are well-separated, while a
            #       plateau boundary concentrates all sign changes in a
            #       vanishingly narrow angular region.
            #
            # Skip the expensive plateau probe when any condition fails.
            _should_probe = False
            w1_area = amoeba_res.get("_w1_area")
            if w1_area is None:
                _, w1_area = _get_average_winding_from_zeros(
                    char_poly, E_ref, mu1, mu2,
                    amoeba_res["zeros"], direction=1,
                )
            if w1_area < plateau_area_threshold:
                _, w2_area = _get_average_winding_from_zeros(
                    char_poly, E_ref, mu1, mu2,
                    amoeba_res["zeros"], direction=2,
                )
                if w2_area < plateau_area_threshold:
                    _should_probe = _check_zeros_are_clustered(
                        amoeba_res["zeros"], plateau_cluster_tol,
                    )

            if _should_probe:
                plateau_info = _probe_zero_plateau_near_mu1(
                    char_poly, E_ref, mu1,
                    amoeba_res.get("_mu1_bracket"),
                    mu2_low=mu2_low,
                    mu2_high=mu2_high,
                    continuum_tol=continuum_tol,
                    continuum_perturb=continuum_perturb,
                    max_iter=max_iter, xtol=xtol,
                    max_range_expansions=max_range_expansions,
                    range_expand_factor=range_expand_factor,
                    winding_tol=plateau_winding_tol,
                    probe_radius=plateau_probe_radius,
                )
                if plateau_info["found"]:
                    is_amoeba = False

        if not is_amoeba:
            subsets = []
        elif amoeba_res["is_continuum"]:
            subsets = _assemble_continuum_subsets(
                zm, E_ref, mu1, mu2,
                amoeba_res,
            )
        else:
            subsets = _assemble_discrete_subsets(
                E_ref, mu1, mu2, amoeba_res["zeros"],
            )

        n_0d = sum(1 for s in subsets if isinstance(s, PointSubset))
        n_1d = sum(1 for s in subsets if isinstance(s, LineSubset))
        return GBZResult(E_ref=E_ref, subsets=subsets, index=(n_0d, n_1d))

    except Exception as e:
        if debug_mode:
            raise e
        return GBZResult(E_ref=E_ref, success=False, error=str(e))
