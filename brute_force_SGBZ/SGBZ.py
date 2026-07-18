'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2025-11-23
Copyright © Department of Physics, Tsinghua University. All rights reserved

Top-level SGBZ solver.

Locates the mu1 = log|beta1| at which the strip winding number vanishes
(via bracket expansion + plain bisection with continuum interception),
classifies the result with an optional zero-plateau probe, and packages the
GBZ points of a reference energy as a GBZResult of connected subsets.
Entry point: collect_GBZ_subsets.
'''

import numpy as np
import poly_tools as pt
from cmath import log
from typing import Optional

from .strip_winding_number import get_strip_winding
from .winding import PolyDiffContext
from .pmgbz_detector import get_roots_and_PMGBZ

from gbz_types import (
    PointSubset, LineSubset, GBZResult, ConnectedSubset,
    generate_probe_steps,
)


# ---- internal helpers ----

def _check_pmgbz_points_clustered(
    gbz: GBZResult,
    tol_normalized: float = 1e-2,
) -> bool:
    """Check whether every PMGBZ point theta1 has a neighbour within tol.

    Port of amoeba's ``_check_zeros_are_clustered`` (amoeba.py:32-75):
    at a genuine GBZ point the PMGBZ zeros are well-separated (they
    partition the circle into meaningful segments).  At a zero-plateau
    boundary all zeros cluster into nearly degenerate pairs — each zero
    sits within ``tol_normalized`` of a neighbour.

    Returns True when ALL points have a neighbour (suspicious →
    run probe), False when any point is isolated (genuine GBZ →
    skip expensive probe).
    """
    from math import pi as _pi
    import cmath
    thetas = []
    for s in gbz.subsets:
        if isinstance(s, PointSubset):
            thetas.append(cmath.phase(s.beta1) % (2 * _pi))
    if len(thetas) < 2:
        return False
    tol_rad = tol_normalized * 2 * _pi
    for i, t in enumerate(thetas):
        has_neighbor = False
        for j, t2 in enumerate(thetas):
            if i == j:
                continue
            d = abs(t - t2)
            d = min(d, 2 * _pi - d)
            if d < tol_rad:
                has_neighbor = True
                break
        if not has_neighbor:
            return False
    return True


def _is_zero_plateau_probe(point: dict, zero_tol: float) -> bool:
    """Check whether a probe point lies on a zero plateau.

    True iff the probe succeeded, is not a continuum point, found an empty
    GBZ, and its winding vanishes within zero_tol.
    """
    return (
        point["success"]
        and (not point["is_continuum"])
        and point["gbz_count"] == 0
        and abs(point["winding"]) <= zero_tol
    )


def _resolve_continuum_winding(
    poly_diff: PolyDiffContext,
    E_ref: complex,
    mu1: float,
    N_points: int = 301,
    continuum_perturb: float = 1e-2,
) -> tuple[Optional[float], Optional[float]]:
    """Compute the left / right winding limits at a continuum-degenerate mu1.

    Mirrors amoeba's _resolve_continuum (bisect.py:215): tries
    perturbation scales (1, 2, 4, 8) to escape the degenerate band.
    Both sides are evaluated in sweep mode (refine_continuum=False).

    Returns:
        (w_left, w_right) — each is float or None (if all scales still
        degenerate on that side).
    """
    for scale in (1.0, 2.0, 4.0, 8.0):
        eps = continuum_perturb * scale
        w_left, _ = get_strip_winding(
            poly_diff, E_ref, mu1 - eps, N_points,
            refine_continuum=False,
        )
        w_right, _ = get_strip_winding(
            poly_diff, E_ref, mu1 + eps, N_points,
            refine_continuum=False,
        )
        if (w_left is not None) and (w_right is not None):
            return w_left, w_right
    return None, None


def _probe_zero_plateau_near_mu1(
    poly_diff: PolyDiffContext,
    E_ref: complex,
    mu1: float,
    mu1_bracket: Optional[tuple[float, float]],
    N_points: int = 301,
    zero_tol: float = 1e-10,
    continuum_perturb: float = 1e-2,
    probe_radius: Optional[float] = None,
) -> dict:
    """Check whether a nonempty-PMGBZ candidate sits next to a zero plateau."""
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

    for step in steps:
        for side in (-1, 1):
            mu1_probe = mu1 + side * step
            point = {
                "mu1": mu1_probe,
                "side": side,
                "step": step,
                "success": False,
            }
            winding, gbz = get_strip_winding(
                poly_diff, E_ref, mu1_probe, N_points,
                continuum_perturb=continuum_perturb,
                zero_tol=zero_tol,
                refine_continuum=False,
            )
            is_continuum = (winding is None)
            point.update({
                "success": True,
                "winding": float(winding) if winding is not None else np.nan,
                "raw_winding": winding,
                "gbz_count": len(gbz.subsets),
                "is_continuum": is_continuum,
            })
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


# ---- SGBZ solver (process-level, mirroring amoeba's bisect_amoeba_ronkin_min) ----

def solve_SGBZ_for_E(
    poly_diff: PolyDiffContext,
    E_ref: complex,
    mu1_guess: tuple[float, float] = (-1, 1),
    zero_tol: float = 1e-10,
    N_points: int = 101,
    continuum_perturb: float = 1e-2,
    refine: bool = True,
    max_iter: int = 60,
    xtol: float = 2e-12,
) -> dict:
    """Locate the strip-winding zero in mu1 and return solve diagnostics.

    Uses bracket expansion + plain bisection (midpoint) that explicitly
    handles continuum-degenerate mu1 values: when the winding evaluates
    to None (continuum), the left/right limits are computed via
    _resolve_continuum_winding.  If they straddle zero the current mu1
    is the SGBZ boundary and the loop returns immediately with
    is_continuum=True.

    Plain bisection is chosen over false-position methods because the
    strip winding has flat plateaus (±1) with a narrow transition zone;
    false position stalls on this shape while midpoint bisection
    guarantees bracket halving every step (mirrors amoeba's
    _find_mu2_for_w2_zero).

    Parameters:
        poly_diff: polynomial evaluation context.
        E_ref: reference energy.
        mu1_guess: initial (left, right) bracket for mu1; each side is
            expanded in unit steps until the winding changes sign.
        zero_tol: tolerance for treating a winding value as zero.
        N_points: theta1 mesh size passed to get_strip_winding.
        continuum_perturb: mu1 offset used for continuum limits.
        refine: if True (default), crisp the gbz for continuum results
            via a full-mode get_roots_and_PMGBZ call before returning.
        max_iter: maximum bisection iterations.
        xtol: minimum bracket width for convergence.

    Returns:
        dict with keys "mu1", "gbz" (GBZResult), "winding" (float or
        None for continuum), "is_continuum", plus debug fields
        "_mu1_bracket", "_winding_bracket", "_w_limits" (continuum
        only), "_exit_reason".
    """

    # --- winding_at: evaluate strip winding in sweep mode ---
    def winding_at(mu1_val: float):
        w, gbz = get_strip_winding(
            poly_diff, E_ref, mu1_val, N_points,
            continuum_perturb=continuum_perturb,
            refine_continuum=False,
        )
        return w, gbz

    # --- handle_continuum: shared logic for bracket/midpoint continuum ---
    def handle_continuum(mu1_val: float, gbz_val: GBZResult):
        """Try to resolve a continuum-degenerate mu1.

        Returns (is_boundary, proxy_w, result_dict).
        is_boundary=True means the winding limits straddle zero →
        this mu1 is the SGBZ point.
        """
        w_l, w_r = _resolve_continuum_winding(
            poly_diff, E_ref, mu1_val, N_points,
            continuum_perturb=continuum_perturb,
        )
        if w_l is None or w_r is None:
            raise ValueError(
                f"Cannot resolve continuum at mu1={mu1_val:.8g}: "
                f"all perturbation scales still degenerate."
            )
        if w_l * w_r <= 0:
            return True, None, {
                "mu1": mu1_val,
                "gbz": gbz_val,
                "winding": None,
                "is_continuum": True,
                "_mu1_bracket": (mu1_low, mu1_high),
                "_winding_bracket": (w_low, w_high),
                "_w_limits": (w_l, w_r),
                "_exit_reason": "continuum_boundary",
            }
        return False, w_l, None  # proxy = w_l (amoeba convention)

    # --- Step 1: bracket expansion ---
    mu1_low = mu1_guess[0]
    mu1_ext_right = None

    while True:
        w_low, gbz_low = winding_at(mu1_low)
        if w_low is None:
            is_boundary, proxy_w, result = handle_continuum(mu1_low, gbz_low)
            if is_boundary:
                return result
            w_low = proxy_w  # use resolved proxy value

        if w_low < zero_tol:
            break
        else:
            mu1_ext_right = mu1_low
            mu1_low -= 1

    if w_low > -zero_tol:
        return {
            "mu1": mu1_low,
            "gbz": gbz_low,
            "winding": w_low,
            "is_continuum": False,
            "_mu1_bracket": (mu1_low, mu1_low),
            "_winding_bracket": (w_low, w_low),
            "_exit_reason": "left_endpoint_zero",
        }

    # Step 1.2: right bracket endpoint
    if mu1_ext_right is None:
        mu1_ext_right = mu1_guess[1]
        while True:
            w_high, gbz_high = winding_at(mu1_ext_right)
            if w_high is None:
                is_boundary, proxy_w, result = handle_continuum(mu1_ext_right, gbz_high)
                if is_boundary:
                    return result
                w_high = proxy_w

            if w_high > -zero_tol:
                break
            else:
                mu1_ext_right += 1

        if w_high < zero_tol:
            return {
                "mu1": mu1_ext_right,
                "gbz": gbz_high,
                "winding": w_high,
                "is_continuum": False,
                "_mu1_bracket": (mu1_ext_right, mu1_ext_right),
                "_winding_bracket": (w_high, w_high),
                "_exit_reason": "right_endpoint_zero",
            }
    else:
        w_high, _ = winding_at(mu1_ext_right)
        if w_high is None:
            is_boundary, proxy_w, _ = handle_continuum(mu1_ext_right, gbz_high)
            w_high = proxy_w

    mu1_high = mu1_ext_right

    # --- Step 2: plain bisection with continuum interception ---
    # The strip winding has flat plateaus (±1) with a narrow transition
    # zone; false-position methods stall on this shape.  Plain midpoint
    # bisection guarantees bracket halving every step (mirrors amoeba's
    # _find_mu2_for_w2_zero and bisect_amoeba_ronkin_min).

    for _ in range(max_iter):
        mu1_mid = 0.5 * (mu1_low + mu1_high)

        w_mid, gbz_mid = winding_at(mu1_mid)

        if w_mid is None:
            is_boundary, proxy_w, result = handle_continuum(mu1_mid, gbz_mid)
            if is_boundary:
                if refine:
                    gbz_refined = get_roots_and_PMGBZ(
                        poly_diff, E_ref, mu1_mid, N_points,
                    )[0]
                    result["gbz"] = gbz_refined
                return result
            f_mid = proxy_w
        else:
            f_mid = w_mid

        # Convergence check
        if abs(f_mid) <= zero_tol or (mu1_high - mu1_low) < xtol:
            exit_reason = "w_zero" if abs(f_mid) <= zero_tol else "bracket_xtol"
            gbz_final = gbz_mid
            winding_final = w_mid
            is_continuum = (w_mid is None)
            if is_continuum and refine:
                gbz_final = get_roots_and_PMGBZ(
                    poly_diff, E_ref, mu1_mid, N_points,
                )[0]
            return {
                "mu1": mu1_mid,
                "gbz": gbz_final,
                "winding": winding_final,
                "is_continuum": is_continuum,
                "_mu1_bracket": (mu1_low, mu1_high),
                "_winding_bracket": (w_low, w_high),
                "_exit_reason": exit_reason,
            }

        # Bracket update
        if w_low * f_mid < 0:
            mu1_high = mu1_mid
        else:
            mu1_low = mu1_mid

    # Max iterations exhausted — last midpoint as fallback
    mu1_final = 0.5 * (mu1_low + mu1_high)
    w_final, gbz_final = winding_at(mu1_final)
    is_continuum = (w_final is None)
    if is_continuum and refine:
        gbz_final = get_roots_and_PMGBZ(
            poly_diff, E_ref, mu1_final, N_points,
        )[0]
    return {
        "mu1": mu1_final,
        "gbz": gbz_final,
        "winding": w_final,
        "is_continuum": is_continuum,
        "_mu1_bracket": (mu1_low, mu1_high),
        "_winding_bracket": (w_low, w_high),
        "_exit_reason": "max_iter",
    }


# ---- main entry point ----

def collect_GBZ_subsets(
    coeffs: np.ndarray,
    degs: np.ndarray,
    E_ref: complex,
    perc: float,
    debug_mode: bool = False,
    **options,
) -> GBZResult:
    """Check the SGBZ condition and return GBZ points for a reference energy.

    Batch entry point: builds the characteristic polynomial from
    (coeffs, degs), solves for the mu1 where the strip winding number
    vanishes using plain bisection with continuum interception, and
    (optionally) reclassifies the candidate as empty when a zero plateau
    exists right next to it.

    For continuum results, the precise beta2 roots are computed once after
    the plateau check, so plateau-detected false positives cost nothing.

    Parameters:
        coeffs: complex coefficients of the characteristic Laurent polynomial
            f(E, beta1, beta2).
        degs: (n_terms, 3) integer exponents of (E, beta1, beta2) per term.
        E_ref: reference energy to test.
        perc: progress fraction in [0, 1], printed as a percentage.
        debug_mode: if True, re-raise solver exceptions instead of returning
            a failed GBZResult.
        **options: solver options — "mu1_guess" (default (-1, 1)),
            "zero_tol" (1e-10), "N_points" (301), "continuum_perturb" (1e-2),
            "plateau_check" (True), "plateau_probe_radius" (None).

    Returns:
        GBZResult with connected subsets.  ``gbz.is_empty`` / ``gbz.index == (0,0)``
        means E_ref is outside the SGBZ spectrum.
    """
    print("%.2f" % (perc * 100) + r"%")
    coeffs_ct = pt.CScalarVec(coeffs)
    degs_ct = pt.CLaurentIndexVec(degs.flatten())
    char_poly = pt.CLaurent(3)
    char_poly.set_Laurent_by_terms(coeffs_ct, degs_ct)
    poly_diff = PolyDiffContext(char_poly)

    solver_options = dict(options)
    plateau_check = solver_options.pop("plateau_check", True)
    plateau_probe_radius = solver_options.pop("plateau_probe_radius", None)
    mu1_guess = solver_options.pop("mu1_guess", (-1, 1))
    zero_tol = solver_options.pop("zero_tol", 1e-10)
    N_points = solver_options.pop("N_points", 301)
    continuum_perturb = solver_options.pop("continuum_perturb", 1e-2)

    try:
        sgbz_res = solve_SGBZ_for_E(
            poly_diff, E_ref, mu1_guess=mu1_guess, zero_tol=zero_tol,
            N_points=N_points, continuum_perturb=continuum_perturb,
            refine=False,
        )
        gbz: GBZResult = sgbz_res["gbz"]
        mu1 = sgbz_res["mu1"]
        is_continuum = sgbz_res["is_continuum"]
    except Exception as e:
        if debug_mode:
            raise e
        print("Error: %s" % str(e))
        return GBZResult(E_ref=E_ref, success=False, error=str(e))

    if not gbz.is_empty:
        # Run plateau check with lightweight pre-filter (mirrors amoeba).
        if plateau_check:
            should_probe = False
            if is_continuum:
                # Continuum boundary — genuine GBZ, no plateau to detect.
                should_probe = False
            else:
                # Point case: only probe when PMGBZ zeros are suspiciously
                # clustered (plateau signature).  Well-separated zeros →
                # genuine GBZ, skip expensive probe.
                should_probe = _check_pmgbz_points_clustered(gbz)

            if should_probe:
                plateau_info = _probe_zero_plateau_near_mu1(
                    poly_diff, E_ref, mu1, sgbz_res.get("_mu1_bracket"),
                    N_points=N_points, zero_tol=zero_tol,
                    continuum_perturb=continuum_perturb,
                    probe_radius=plateau_probe_radius,
                )
                if plateau_info["found"]:
                    gbz = GBZResult(E_ref=E_ref, subsets=[], index=(0, 0))

    # Deferred precise solve for continuum results (after plateau check).
    if is_continuum and gbz.is_gbz:
        gbz = get_roots_and_PMGBZ(
            poly_diff, E_ref, mu1, N_points,
        )[0]

    return gbz
