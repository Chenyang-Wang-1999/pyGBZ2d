'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-10
Copyright © Department of Physics, Tsinghua University. All rights reserved

Top-level SGBZ solver — ZeroManager-based μ₁ bisection.

Locates the ``mu1 = log|beta1|`` at which the average major-axis winding
number ``W(E_ref, mu1)`` vanishes, via bracket expansion + plain midpoint
bisection with continuum interception.  When a continuum-degenerate
``mu1`` is encountered, the left/right winding limits are resolved by
``mu1`` perturbation; if they straddle zero that ``mu1`` IS the SGBZ
boundary (a 1D LineSubset case — materialized by ``extract_continuum_linesubsets``).

Entry point: :func:`collect_GBZ_subsets`.
'''

from __future__ import annotations

from typing import Optional

import numpy as np

from gbz_types import (
    PointSubset, LineSubset, GBZResult, CharPoly,
)

from continuation import ZeroManager

from .continuum_lines import (
    detect_continuum_simple, extract_continuum_linesubsets,
    CONTINUUM_TOL, _DV_TOL, CONTINUUM_FRAC,
)
from .winding import detect_crossings_and_winding
from .plateau import _check_pmgbz_points_clustered, _probe_zero_plateau_near_mu1
from .crossings import _CROSSING_TOL, _DETECT_THRESHOLD, _MAX_NEWTON_ITER, _DEDUP_TOL
from .mu2mid import Mu2MidZM


# ---------------------------------------------------------------------------
# W(E_ref, mu1) evaluation + continuum resolution
# ---------------------------------------------------------------------------

def _evaluate_winding(
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
) -> tuple[Optional[float], Optional[list], Mu2MidZM]:
    """Evaluate W(E_ref, mu1) and the 0D subsets at *mu1*.

    Builds a fresh ``Mu2MidZM`` at *mu1* (mu1 is the bisection variable, so
    each probe needs its own ZM) and runs :meth:`build_mu2_mid` — the build
    both detects the continuum inline (§1/§6.3: ``has_continuum``) and
    provides the piecewise-smooth path the winding integral needs, so the
    root-solving cost is amortised (§6.3).  When a continuum is detected W is
    undefined → returns ``(None, None, zm)``; the caller resolves it via
    left/right limits and, if it is the boundary, materialises the
    LineSubsets from the built *zm*.  Otherwise runs crossing detection +
    winding on the same built *zm* (no second build).

    Returns ``(winding, subsets, zm)``:
      * ``winding`` is ``float``, or ``None`` when a continuum is present.
      * ``subsets`` is the 0D ``PointSubset`` list, or ``None`` when a
        continuum is present (LineSubsets are materialised separately by the
        caller, only when this μ₁ is confirmed as the boundary).
    """
    zm = Mu2MidZM(poly, E_ref, mu1)
    zm.run(**zm_run_kwargs)
    zm.build_mu2_mid(
        tie_tol=continuum_tol,
        # dV_tol / vote_frac are retained for API symmetry; the inline
        # continuum detection uses the whole-segment same-modulus criterion.
    )

    if zm.has_continuum:
        return None, None, zm

    subsets, W = detect_crossings_and_winding(
        zm, poly,
        crossing_tol=crossing_tol,
        detect_threshold=detect_threshold,
        max_newton=max_newton,
        dedup_tol=dedup_tol,
    )
    return W, subsets, zm


def _resolve_continuum_winding(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    zm_run_kwargs: dict,
    *,
    continuum_perturb: float,
    continuum_tol: float,
    dV_tol: float,
    vote_frac: float,
    crossing_tol: float,
    detect_threshold: float,
    max_newton: int,
    dedup_tol: float,
) -> tuple[Optional[float], Optional[float]]:
    """Compute the left / right winding limits at a continuum-degenerate mu1.

    Tries perturbation scales (1, 2, 4, 8) × *continuum_perturb* to escape
    the degenerate band.  Both sides are evaluated in sweep mode (the
    continuum gate short-circuits; the precise 0D subsets are not needed
    here — only the winding sign matters for the bisection).

    Returns ``(w_left, w_right)`` — each is ``float`` or ``None`` (if all
    scales still degenerate on that side).
    """
    eval_kwargs = dict(
        continuum_tol=continuum_tol, dV_tol=dV_tol, vote_frac=vote_frac,
        crossing_tol=crossing_tol, detect_threshold=detect_threshold,
        max_newton=max_newton, dedup_tol=dedup_tol,
    )
    for scale in (1.0, 2.0, 4.0, 8.0):
        eps = continuum_perturb * scale
        w_left, _, _ = _evaluate_winding(
            poly, E_ref, mu1 - eps, zm_run_kwargs, **eval_kwargs,
        )
        w_right, _, _ = _evaluate_winding(
            poly, E_ref, mu1 + eps, zm_run_kwargs, **eval_kwargs,
        )
        if (w_left is not None) and (w_right is not None):
            return w_left, w_right
    return None, None


# ---------------------------------------------------------------------------
# μ₁ bisection (bracket expansion + plain midpoint with continuum interception)
# ---------------------------------------------------------------------------

def solve_SGBZ_for_E(
    poly: CharPoly,
    E_ref: complex,
    mu1_guess: tuple[float, float] = (-1, 1),
    zero_tol: float = 1e-10,
    continuum_perturb: float = 1e-2,
    max_iter: int = 60,
    xtol: float = 2e-12,
    zm_run_kwargs: Optional[dict] = None,
    *,
    continuum_tol: float = CONTINUUM_TOL,
    dV_tol: float = _DV_TOL,
    vote_frac: float = CONTINUUM_FRAC,
    crossing_tol: float = _CROSSING_TOL,
    detect_threshold: float = _DETECT_THRESHOLD,
    max_newton: int = _MAX_NEWTON_ITER,
    dedup_tol: float = _DEDUP_TOL,
) -> dict:
    """Locate the winding-zero mu1 and return solve diagnostics.

    Uses bracket expansion + plain bisection (midpoint).  When a
    continuum-degenerate mu1 is encountered, ``_resolve_continuum_winding``
    resolves the left/right winding limits; if they straddle zero that mu1
    is the SGBZ boundary (``is_continuum=True``) and the 1D LineSubsets are
    materialized (``extract_continuum_linesubsets``) from the built ZM; the
    returned ``subsets`` carries them (the caller signals "in spectrum").

    Plain bisection is chosen over false-position methods because the
    winding has flat plateaus (±1) with a narrow transition zone;
    false position stalls on this shape while midpoint bisection
    guarantees bracket halving every step (mirrors amoeba's
    ``_find_mu2_for_w2_zero``).

    Returns:
        dict with keys "mu1", "subsets" (list[PointSubset] or None for
        continuum), "winding" (float or None), "is_continuum", plus debug
        fields "_mu1_bracket", "_winding_bracket", "_w_limits"
        (continuum only), "_exit_reason".
    """
    if zm_run_kwargs is None:
        zm_run_kwargs = {}

    eval_kwargs = dict(
        continuum_tol=continuum_tol, dV_tol=dV_tol, vote_frac=vote_frac,
        crossing_tol=crossing_tol, detect_threshold=detect_threshold,
        max_newton=max_newton, dedup_tol=dedup_tol,
    )

    # --- winding_at: evaluate average major-axis winding in sweep mode ---
    def winding_at(mu1_val: float):
        return _evaluate_winding(
            poly, E_ref, mu1_val, zm_run_kwargs, **eval_kwargs,
        )

    # --- handle_continuum: resolve a continuum-degenerate mu1 ---
    # When the winding limits straddle zero this mu1 IS the SGBZ boundary
    # (a 1D LineSubset).  The LineSubsets are materialised from the *already
    # built* ``zm`` returned by ``_evaluate_winding`` (no second build) —
    # the continuum tracks are the cluster columns where the boundary pair
    # collapses (j_lo == j_hi), joined across MRs by the amoeba machinery.
    def handle_continuum(
        mu1_val: float,
        zm: Mu2MidZM | None,
        mu1_bracket: tuple[float, float | None],
    ):
        w_l, w_r = _resolve_continuum_winding(
            poly, E_ref, mu1_val, zm_run_kwargs,
            continuum_perturb=continuum_perturb, **eval_kwargs,
        )
        if w_l is None or w_r is None:
            raise ValueError(
                f"Cannot resolve continuum at mu1={mu1_val:.8g}: "
                f"all perturbation scales still degenerate."
            )
        if w_l * w_r < 0:
            line_subsets = (extract_continuum_linesubsets(zm, poly)
                            if zm is not None else [])
            return True, line_subsets, {
                "mu1": mu1_val,
                "subsets": line_subsets,
                "winding": None,
                "is_continuum": True,
                # The mu1 bracket is the best-known enclosing bracket at
                # this call site (the right end may still be None during
                # left-end expansion — see the call sites below).  The
                # winding bracket is the straddling left/right limit pair,
                # which is the decisive winding info at a continuum boundary
                # and is always defined here (unlike the bisection bracket
                # w_low/w_high, which may be incomplete during expansion).
                "_mu1_bracket": mu1_bracket,
                "_winding_bracket": (w_l, w_r),
                "_w_limits": (w_l, w_r),
                "_exit_reason": "continuum_boundary",
            }
        return False, w_l, None  # proxy = w_l (amoeba convention)

    # --- Step 1: bracket expansion ---
    # Initialize bracket variables before defining handle_continuum closure.
    # The closure reads these at call time; during expansion the right bracket
    # may not yet be established, so we initialize to None and set them later.
    mu1_low = mu1_guess[0]
    mu1_high = None
    w_low = None
    w_high = None
    mu1_ext_right = None

    while True:
        w_low, gbz_low, zm_low = winding_at(mu1_low)
        if w_low is None:
            is_boundary, proxy_w, result = handle_continuum(
                mu1_low, zm_low, (mu1_low, mu1_ext_right))
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
            "subsets": gbz_low,
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
            w_high, gbz_high, zm_high = winding_at(mu1_ext_right)
            if w_high is None:
                is_boundary, proxy_w, result = handle_continuum(
                    mu1_ext_right, zm_high, (mu1_low, mu1_ext_right))
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
                "subsets": gbz_high,
                "winding": w_high,
                "is_continuum": False,
                "_mu1_bracket": (mu1_ext_right, mu1_ext_right),
                "_winding_bracket": (w_high, w_high),
                "_exit_reason": "right_endpoint_zero",
            }
    else:
        w_high, _, zm_high = winding_at(mu1_ext_right)
        if w_high is None:
            is_boundary, proxy_w, _ = handle_continuum(
                mu1_ext_right, zm_high, (mu1_low, mu1_ext_right))
            w_high = proxy_w

    mu1_high = mu1_ext_right

    # --- Step 2: plain bisection with continuum interception ---
    # The winding has flat plateaus (±1) with a narrow transition zone;
    # false-position methods stall on this shape.  Plain midpoint bisection
    # guarantees bracket halving every step (mirrors amoeba's
    # _find_mu2_for_w2_zero and bisect_amoeba_ronkin_min).

    for _ in range(max_iter):
        mu1_mid = 0.5 * (mu1_low + mu1_high)

        w_mid, gbz_mid, zm_mid = winding_at(mu1_mid)

        if w_mid is None:
            is_boundary, proxy_w, result = handle_continuum(
                mu1_mid, zm_mid, (mu1_low, mu1_high))
            if is_boundary:
                return result
            # At a non-boundary continuum, W is undefined here.  Use the proxy
            # for bracket update but do NOT converge on a proxy value — the
            # bisection target is the genuine W-zero, not the left/right limit
            # of a continuum point.
            f_mid = proxy_w
        else:
            f_mid = w_mid
            # Convergence check only applies at genuine winding values.
            # At a continuum (w_mid is None) we continue bisection.
            if abs(f_mid) <= zero_tol or (mu1_high - mu1_low) < xtol:
                exit_reason = "w_zero" if abs(f_mid) <= zero_tol else "bracket_xtol"
                return {
                    "mu1": mu1_mid,
                    "subsets": gbz_mid,
                    "winding": w_mid,
                    "is_continuum": False,
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
    w_final, gbz_final, _ = winding_at(mu1_final)
    is_continuum = (w_final is None)
    return {
        "mu1": mu1_final,
        "subsets": gbz_final,
        "winding": w_final,
        "is_continuum": is_continuum,
        "_mu1_bracket": (mu1_low, mu1_high),
        "_winding_bracket": (w_low, w_high),
        "_exit_reason": "max_iter",
    }


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def collect_GBZ_subsets(
    coeffs: np.ndarray,
    degs: np.ndarray,
    E_ref: complex,
    perc: float = None,
    debug_mode: bool = False,
    **options,
) -> GBZResult:
    """Check the SGBZ condition and return GBZ points for a reference energy.

    Builds the characteristic polynomial from (coeffs, degs), solves for the
    mu1 where the average major-axis winding number vanishes using plain
    bisection with continuum interception, and (optionally) reclassifies the
    candidate as empty when a zero plateau exists right next to it.

    For continuum results (a 1D LineSubset case), subset materialization is
    TODO: the returned ``GBZResult`` is marked ``is_continuum=True`` (in
    spectrum) with no subsets.  Spectrum membership still works: the
    bisection's W-zero / left-right-limit straddle decides in vs out.

    Parameters:
        coeffs: complex coefficients of the characteristic Laurent polynomial
            f(E, beta1, beta2).
        degs: (n_terms, 3) integer exponents of (E, beta1, beta2) per term.
        E_ref: reference energy to test.
        perc: progress fraction in [0, 1], printed as a percentage.
        debug_mode: if True, re-raise solver exceptions instead of returning
            a failed GBZResult.
        **options: solver options — "mu1_guess" (default (-1, 1)),
            "zero_tol" (1e-10), "continuum_perturb" (1e-2), "max_iter" (60),
            "xtol" (2e-12), "plateau_check" (True), "plateau_probe_radius"
            (None), "zm_run_kwargs" ({}), plus continuum/crossing tunables.

    Returns:
        GBZResult with connected subsets.  ``gbz.is_empty`` / ``gbz.index
        == (0,0)`` (without ``is_continuum``) means E_ref is outside the
        SGBZ spectrum.  ``gbz.is_continuum`` means in-spectrum; the
        LineSubsets are materialized in ``subsets`` (``index == (0, n_1d)``).
    """
    if perc is not None:
        print("%.2f" % (perc * 100) + r"%")
    poly = CharPoly(coeffs, degs)

    solver_options = dict(options)
    plateau_check = solver_options.pop("plateau_check", True)
    plateau_probe_radius = solver_options.pop("plateau_probe_radius", None)
    mu1_guess = solver_options.pop("mu1_guess", (-1, 1))
    zero_tol = solver_options.pop("zero_tol", 1e-10)
    continuum_perturb = solver_options.pop("continuum_perturb", 1e-2)
    max_iter = solver_options.pop("max_iter", 60)
    xtol = solver_options.pop("xtol", 2e-12)
    zm_run_kwargs = solver_options.pop("zm_run_kwargs", {})

    # Continuum / crossing tunables (rarely overridden).
    continuum_tol = solver_options.pop("continuum_tol", CONTINUUM_TOL)
    dV_tol = solver_options.pop("dV_tol", _DV_TOL)
    vote_frac = solver_options.pop("vote_frac", CONTINUUM_FRAC)
    crossing_tol = solver_options.pop("crossing_tol", _CROSSING_TOL)
    detect_threshold = solver_options.pop("detect_threshold", _DETECT_THRESHOLD)
    max_newton = solver_options.pop("max_newton", _MAX_NEWTON_ITER)
    dedup_tol = solver_options.pop("dedup_tol", _DEDUP_TOL)

    # N_points is obsolete under the adaptive ZeroManager mesh; accept it
    # silently for API compatibility with older callers.
    solver_options.pop("N_points", None)
    if solver_options:
        import warnings
        warnings.warn(
            f"collect_GBZ_subsets: ignoring unrecognized options: "
            f"{sorted(solver_options)}"
        )

    try:
        sgbz_res = solve_SGBZ_for_E(
            poly, E_ref, mu1_guess=mu1_guess, zero_tol=zero_tol,
            continuum_perturb=continuum_perturb, max_iter=max_iter,
            xtol=xtol, zm_run_kwargs=zm_run_kwargs,
            continuum_tol=continuum_tol, dV_tol=dV_tol, vote_frac=vote_frac,
            crossing_tol=crossing_tol, detect_threshold=detect_threshold,
            max_newton=max_newton, dedup_tol=dedup_tol,
        )
        subsets = sgbz_res["subsets"]
        mu1 = sgbz_res["mu1"]
        is_continuum = sgbz_res["is_continuum"]
    except Exception as e:
        if debug_mode:
            raise e
        print("Error: %s" % str(e))
        return GBZResult(E_ref=E_ref, success=False, error=str(e))

    # Continuum boundary → 1D LineSubsets materialised by handle_continuum
    # (the cluster tracks where |β_M| = |β_{M+1}| holds identically, joined
    # across MRs).  is_continuum stays True so existing spectrum-membership
    # assertions still hold; subsets now carry the actual lines.
    if is_continuum:
        line_subsets = subsets or []
        n_1d = len(line_subsets)
        return GBZResult(E_ref=E_ref, success=True, is_continuum=True,
                         subsets=list(line_subsets), index=(0, n_1d))

    # Discrete case: apply the plateau check before trusting the subsets.
    if subsets and plateau_check:
        candidate = GBZResult(E_ref=E_ref, subsets=subsets,
                              index=(len(subsets), 0))
        if _check_pmgbz_points_clustered(candidate):
            plateau_info = _probe_zero_plateau_near_mu1(
                poly, E_ref, mu1, sgbz_res.get("_mu1_bracket"),
                zm_run_kwargs,
                continuum_tol=continuum_tol, dV_tol=dV_tol,
                vote_frac=vote_frac, crossing_tol=crossing_tol,
                detect_threshold=detect_threshold, max_newton=max_newton,
                dedup_tol=dedup_tol, zero_tol=zero_tol,
                continuum_perturb=continuum_perturb,
                probe_radius=plateau_probe_radius,
            )
            if plateau_info["found"]:
                subsets = []

    n_0d = sum(1 for s in subsets if isinstance(s, PointSubset))
    n_1d = sum(1 for s in subsets if isinstance(s, LineSubset))
    return GBZResult(E_ref=E_ref, subsets=list(subsets), index=(n_0d, n_1d))
