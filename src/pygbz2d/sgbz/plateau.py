'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-10
Copyright © Department of Physics, Tsinghua University. All rights reserved

Zero-plateau detection for the SGBZ solver.

Nearly degenerate pairs of PMGBZ points (0D ``PointSubset``s) can signal a
zero-plateau boundary. The clustering pre-check uses both angular
coordinates. When every point has a nearby neighbour,
``_probe_zero_plateau_near_mu1`` evaluates ``μ₁ ± step``. After evaluating
both sides of a step, a plateau is reported if either probe has zero
winding and an empty GBZ.
'''

from __future__ import annotations

import cmath
from typing import Optional

from ..core import (
    TWO_PI,
    CharPoly, GBZResult, PointSubset,
    check_points_clustered_on_torus, probe_zero_plateau,
    live_defaults
)

from .mu2mid import Mu2MidZM
from .winding import detect_crossings_and_winding


# ---- clustering pre-check ----

@live_defaults(tol_normalized="core:PLATEAU_CLUSTER_TOL")
def _check_pmgbz_points_clustered(
    gbz: GBZResult,
    tol_normalized: Optional[float] = None,
) -> bool:
    """Whether every PMGBZ point has a neighbour within *tol_normalized*.

    Thin SGBZ adapter over
    :func:`pygbz2d.core.check_points_clustered_on_torus`: extracts ``(θ₁, θ₂)``
    from the result's ``PointSubset``s and forwards the angular pairs.

    Nearly degenerate pairs can indicate winding changes confined to a
    narrow angular region near a plateau boundary. Clustering only selects
    candidates for the winding probe; it does not establish a plateau by
    itself. The tolerance is normalized by the angular period ``2π``.
    """
    twopi = TWO_PI
    points: list[tuple[float, float]] = [
        (cmath.phase(s.beta1) % twopi, cmath.phase(s.beta2) % twopi)
        for s in gbz.subsets if isinstance(s, PointSubset)
    ]
    return check_points_clustered_on_torus(points, tol_normalized)


# ---- plateau probe ----

def _evaluate_probe(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    *,
    continuum_tol: float,
    crossing_tol: float,
) -> dict:
    """Build ONE Mu2MidZM at (E_ref, mu1) and evaluate the winding + subset count.

    The single build serves both stages — the continuum gate (inline
    ``has_continuum``) and the crossing + winding — mirroring
    ``sgbz_solver._evaluate_winding`` (the old version built three ZMs per
    probe: a plain ZeroManager, one inside ``detect_continuum_simple`` and
    another inside ``detect_crossings_and_winding``).

    Returns a dict with ``success``, ``is_continuum``, ``winding``, ``gbz_count``.
    Errors propagate directly — a failing probe must not be silently skipped.
    """
    m = Mu2MidZM(poly, E_ref, mu1)
    m.run()
    m.analyze(tie_tol=continuum_tol, crossing_tol=crossing_tol)
    if m.has_continuum:
        return {
            "success": True, "is_continuum": True,
            "winding": float('nan'), "gbz_count": 0,
        }
    subsets, W = detect_crossings_and_winding(
        m, poly, crossing_tol=crossing_tol,
    )
    return {
        "success": True, "is_continuum": False,
        "winding": float(W), "gbz_count": len(subsets),
    }


def _is_zero_plateau_probe(point: dict, zero_tol: float) -> bool:
    """True iff the probe succeeded, is not a continuum point, found an empty
    GBZ, and its winding vanishes within zero_tol."""
    return (
        point["success"]
        and (not point["is_continuum"])
        and point["gbz_count"] == 0
        and abs(point["winding"]) <= zero_tol
    )


@live_defaults(zero_tol="core:WINDING_ZERO_TOL", continuum_perturb="core:CONTINUUM_PERTURB")
def _probe_zero_plateau_near_mu1(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    mu1_bracket: Optional[tuple[float, float]],
    *,
    continuum_tol: float,
    crossing_tol: float,
    zero_tol: Optional[float] = None,
    continuum_perturb: Optional[float] = None,
    probe_radius: Optional[float] = None,
) -> dict:
    """Check whether a nonempty-PMGBZ candidate sits next to a zero plateau.

    Thin SGBZ adapter over :func:`pygbz2d.core.probe_zero_plateau`: the step
    ladder, ``±side`` loop and found/not_found/inconclusive classification
    are shared; the per-probe *evaluation* (build one Mu2MidZM, gate on
    continuum, else run crossing detection + winding) and the plateau
    criterion (empty GBZ + zero winding) are SGBZ-specific, supplied as the
    ``evaluator`` closure.
    """
    if probe_radius is None:
        probe_radius = continuum_perturb

    bracket_width = 0.0
    if mu1_bracket is not None:
        bracket_width = abs(float(mu1_bracket[1]) - float(mu1_bracket[0]))

    def evaluator(mu1_probe: float) -> dict:
        res = _evaluate_probe(
            poly, E_ref, mu1_probe,
            continuum_tol=continuum_tol, crossing_tol=crossing_tol,
        )
        # Collapse the SGBZ-specific criterion to the one bool the shared
        # loop reads; keep the raw fields for diagnostics.
        res["is_plateau"] = _is_zero_plateau_probe(res, zero_tol)
        return res

    return probe_zero_plateau(
        mu1, mu1_bracket,
        zero_tol=zero_tol,
        probe_radius=probe_radius,
        bracket_width=bracket_width,
        evaluator=evaluator,
    )
