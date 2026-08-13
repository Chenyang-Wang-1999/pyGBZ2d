'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-10
Copyright © Department of Physics, Tsinghua University. All rights reserved

Zero-plateau detection for the SGBZ solver.

At a genuine GBZ point the PMGBZ points (0D ``PointSubset``s) are
well-separated — they partition the θ₁ circle into meaningful segments.
At a zero-plateau boundary all points cluster into nearly degenerate pairs
instead, each within ``tol_normalized`` of a neighbour.  When that
signature is present, ``_probe_zero_plateau_near_mu1`` verifies by probing
``μ₁ ± step``: a plateau shows zero winding with empty GBZ on both sides.
'''

from __future__ import annotations

import math
import cmath
from typing import Optional

import numpy as np

from gbz_types import (
    CharPoly, GBZResult, PointSubset, circ_dist, generate_probe_steps,
)
from continuation import ZeroManager

from .continuum_lines import detect_continuum_simple
from .winding import detect_crossings_and_winding


# ---- clustering pre-check ----

def _check_pmgbz_points_clustered(
    gbz: GBZResult,
    tol_normalized: float = 1e-2,
) -> bool:
    """Check whether every PMGBZ point has a neighbour within tol.

    Port of amoeba's ``_check_zeros_are_clustered`` (amoeba.py:32-75).
    Uses Euclidean distance on the (θ₁, θ₂) torus, normalized by 2π.

    At a genuine GBZ point the PMGBZ points are well-separated (they
    partition the circle into meaningful segments).  At a zero-plateau
    boundary the winding changes sign over a vanishingly narrow angular
    region, so the points cluster into nearly degenerate pairs — each
    point sits within ``tol_normalized`` of a neighbour.

    Returns True when ALL points have a neighbour (suspicious →
    run probe), False when any point is isolated (genuine GBZ →
    skip expensive probe).
    """
    twopi = 2 * math.pi

    # Extract (theta1, theta2) pairs from PointSubsets
    points: list[tuple[float, float]] = []
    for s in gbz.subsets:
        if isinstance(s, PointSubset):
            theta1 = cmath.phase(s.beta1) % twopi
            theta2 = cmath.phase(s.beta2) % twopi
            points.append((theta1, theta2))

    if len(points) < 2:
        return False

    # Tolerance in radians (normalized by 2π)
    tol_rad = tol_normalized * twopi

    for i in range(len(points)):
        t1_i, t2_i = points[i]
        has_neighbor = False
        for j in range(len(points)):
            if i == j:
                continue
            t1_j, t2_j = points[j]
            # Euclidean distance on (θ₁, θ₂) torus
            d1 = circ_dist(t1_i, t1_j)
            d2 = circ_dist(t2_i, t2_j)
            d = math.sqrt(d1 * d1 + d2 * d2)
            if d < tol_rad:
                has_neighbor = True
                break
        if not has_neighbor:
            return False
    return True


# ---- plateau probe ----

def _evaluate_probe(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    zm_run_kwargs: dict,
    *,
    continuum_tol: float,
    dV_tol: float,
    vote_frac: float,
    crossing_tol: float,
    detect_threshold: float,
    max_newton: int,
    dedup_tol: float,
) -> dict:
    """Build a ZeroManager at (E_ref, mu1) and evaluate the winding + subset count.

    Returns a dict with ``success``, ``is_continuum``, ``winding``, ``gbz_count``.
    ``success=False`` if the ZeroManager raised (the probe is skipped then).
    """
    point: dict = {"success": False}
    try:
        zm = ZeroManager(poly, E_ref, mu1)
        zm.run(**zm_run_kwargs)
        if detect_continuum_simple(
            zm, poly,
            continuum_tol=continuum_tol, dV_tol=dV_tol, vote_frac=vote_frac,
        ):
            point.update({
                "success": True, "is_continuum": True,
                "winding": float('nan'), "gbz_count": 0,
            })
            return point
        subsets, W = detect_crossings_and_winding(
            zm, poly,
            crossing_tol=crossing_tol,
            detect_threshold=detect_threshold,
            max_newton=max_newton,
            dedup_tol=dedup_tol,
        )
        point.update({
            "success": True, "is_continuum": False,
            "winding": float(W), "gbz_count": len(subsets),
        })
    except Exception:
        # A failing probe is not fatal — the caller treats it as non-plateau.
        point.update({"success": False, "is_continuum": False,
                      "winding": float('nan'), "gbz_count": -1})
    return point


def _is_zero_plateau_probe(point: dict, zero_tol: float) -> bool:
    """True iff the probe succeeded, is not a continuum point, found an empty
    GBZ, and its winding vanishes within zero_tol."""
    return (
        point["success"]
        and (not point["is_continuum"])
        and point["gbz_count"] == 0
        and abs(point["winding"]) <= zero_tol
    )


def _probe_zero_plateau_near_mu1(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    mu1_bracket: Optional[tuple[float, float]],
    zm_run_kwargs: dict,
    *,
    continuum_tol: float,
    dV_tol: float,
    vote_frac: float,
    crossing_tol: float,
    detect_threshold: float,
    max_newton: int,
    dedup_tol: float,
    zero_tol: float = 1e-10,
    continuum_perturb: float = 1e-2,
    probe_radius: Optional[float] = None,
) -> dict:
    """Check whether a nonempty-PMGBZ candidate sits next to a zero plateau.

    Probes ``mu1 ± step`` for a geometric ladder of steps.  A plateau is
    found when a probe succeeds, is not a continuum, has empty GBZ, and its
    winding vanishes within *zero_tol*.
    """
    if probe_radius is None:
        probe_radius = continuum_perturb

    bracket_width = 0.0
    if mu1_bracket is not None:
        bracket_width = abs(float(mu1_bracket[1]) - float(mu1_bracket[0]))

    steps = generate_probe_steps(bracket_width, probe_radius, zero_tol)

    probe_points = []
    found_plateau = False
    saw_left_nonplateau = False
    saw_right_nonplateau = False

    eval_kwargs = dict(
        continuum_tol=continuum_tol, dV_tol=dV_tol, vote_frac=vote_frac,
        crossing_tol=crossing_tol, detect_threshold=detect_threshold,
        max_newton=max_newton, dedup_tol=dedup_tol,
    )

    for step in steps:
        for side in (-1, 1):
            mu1_probe = mu1 + side * step
            point = {
                "mu1": mu1_probe,
                "side": side,
                "step": step,
                "success": False,
            }
            res = _evaluate_probe(
                poly, E_ref, mu1_probe, zm_run_kwargs, **eval_kwargs,
            )
            point.update(res)
            if _is_zero_plateau_probe(point, zero_tol):
                found_plateau = True
            elif point["success"] and (not point["is_continuum"]):
                if side < 0:
                    saw_left_nonplateau = True
                else:
                    saw_right_nonplateau = True
            probe_points.append(point)
            if found_plateau:
                return {
                    "status": "found",
                    "found": True,
                    "zero_tol": zero_tol,
                    "probe_radius": probe_radius,
                    "bracket_width": bracket_width,
                    "steps": steps,
                    "points": probe_points,
                }

    if found_plateau:
        status = "found"
    elif saw_left_nonplateau and saw_right_nonplateau:
        status = "not_found"
    else:
        status = "inconclusive"

    return {
        "status": status,
        "found": found_plateau,
        "zero_tol": zero_tol,
        "probe_radius": probe_radius,
        "bracket_width": bracket_width,
        "steps": steps,
        "points": probe_points,
    }
