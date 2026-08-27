'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-05-19
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

from typing import Optional
import numpy as np

from bfgbz2d import config
from bfgbz2d.core import (
    PointSubset, LineSubset, GBZResult, CharPoly,
    check_points_clustered_on_torus, probe_zero_plateau,
)

from .bisect import (
    bisect_amoeba_ronkin_min, _find_mu2_for_w2_zero,
)
from .ronkin_winding import (
    _get_average_winding_from_zeros,
)
from .zm_extract import extract_amoeba_subsets


# ---- plateau pre-check helpers ----

def _check_zeros_are_clustered(
    zeros: list[tuple[float, float, int]],
    tol_normalized: float,
) -> bool:
    """Whether every amoeba zero has a neighbour within *tol_normalized*.

    Thin amoeba adapter over
    :func:`bfgbz2d.core.check_points_clustered_on_torus`: amoeba zeros carry a
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


def _probe_zero_plateau_near_mu1(
    char_poly: CharPoly,
    E_ref: complex,
    mu1: float,
    mu1_bracket: Optional[tuple[float, float]],
    mu2_low: float = -1,
    mu2_high: float = 1,
    continuum_tol: float = config.CONTINUUM_TOL,
    continuum_perturb: float = 1e-4,
    max_iter: int = 60,
    xtol: float = 1e-10,
    max_range_expansions: int = 10,
    range_expand_factor: float = 2.0,
    winding_tol: Optional[float] = None,
    probe_radius: Optional[float] = None,
) -> dict:
    """Check whether a nonempty-zero candidate sits next to a zero plateau.

    The decisive plateau signature is strict: after solving a2=0 at a nearby
    mu1, w1 is zero and there are no a2-crossing zeros.
    """
    if winding_tol is None:
        winding_tol = max(xtol, config.WINDING_TOL_FLOOR)
    if probe_radius is None:
        probe_radius = continuum_perturb

    bracket_width = 0.0
    if mu1_bracket is not None:
        bracket_width = abs(float(mu1_bracket[1]) - float(mu1_bracket[0]))

    def evaluator(mu1_probe: float) -> dict:
        # The shared ladder protocol (bfgbz2d.core.probe_zero_plateau) reads
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


# ---- main entry point ----

def collect_GBZ_subsets(
    coeffs: np.ndarray,
    degs: np.ndarray,
    E_ref: complex,
    perc: float = None,
    debug_mode: bool = False,
    **options,
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

    solver_options = dict(options)
    plateau_check = solver_options.pop("plateau_check", True)
    plateau_winding_tol = solver_options.pop("plateau_winding_tol", None)
    plateau_probe_radius = solver_options.pop("plateau_probe_radius", None)
    plateau_area_threshold = solver_options.pop("plateau_area_threshold", config.PLATEAU_AREA_THRESHOLD)
    # Neighbour threshold for the torus-clustering pre-check, in units of
    # the 2π torus period.  Deliberately a SEPARATE knob from
    # plateau_area_threshold (a winding-area fraction): the two criteria
    # have different units, and sharing one value couples their tuning.
    # Defaults to the historical shared value for behaviour compatibility.
    plateau_cluster_tol = solver_options.pop("plateau_cluster_tol", config.PLATEAU_CLUSTER_TOL)
    # kwargs forwarded to ZeroManager.run() (h0, ctrl, min_dtheta,
    # cluster_tol, mr_jump, verbose).  Kept separate from the bisection
    # options, which ZeroManager.run does not accept.
    zm_run_kwargs = solver_options.pop("zm_run_kwargs", {})

    try:
        amoeba_res = bisect_amoeba_ronkin_min(
            char_poly, E_ref, zm_run_kwargs=zm_run_kwargs, **solver_options,
        )

        # Reuse the ZeroManager built inside the bisection at the solved
        # mu1 — subset extraction runs on the same adaptive tracks, no
        # second ZeroManager run.
        zm = amoeba_res["_zm"]
        mu1, mu2 = amoeba_res["mu1"], amoeba_res["mu2"]
        subsets = extract_amoeba_subsets(
            zm, char_poly, E_ref, mu1, mu2, mode='solve',
        )

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
                    mu2_low=solver_options.get("mu2_low", -1),
                    mu2_high=solver_options.get("mu2_high", 1),
                    continuum_tol=solver_options.get("continuum_tol", config.CONTINUUM_TOL),
                    continuum_perturb=solver_options.get("continuum_perturb", 1e-4),
                    max_iter=solver_options.get("max_iter", 60),
                    xtol=solver_options.get("xtol", 1e-10),
                    max_range_expansions=solver_options.get("max_range_expansions", 10),
                    range_expand_factor=solver_options.get("range_expand_factor", 2.0),
                    winding_tol=plateau_winding_tol,
                    probe_radius=plateau_probe_radius,
                )
                if plateau_info["found"]:
                    is_amoeba = False

        if not is_amoeba:
            subsets = []

        n_0d = sum(1 for s in subsets if isinstance(s, PointSubset))
        n_1d = sum(1 for s in subsets if isinstance(s, LineSubset))
        return GBZResult(E_ref=E_ref, subsets=subsets, index=(n_0d, n_1d))

    except Exception as e:
        if debug_mode:
            raise e
        return GBZResult(E_ref=E_ref, success=False, error=str(e))
