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

from ..core import (
    PointSubset, LineSubset, GBZResult, CharPoly, live_defaults
)
from .. import core

from . import pairwise as _pairwise

#: Iteration budget of the μ₁ bisection.
MU1_MAX_ITER: int = 60
#: Cap on μ₁ bracket-expansion steps.
MAX_BRACKET_EXPANSIONS: int = 10
from .continuum_lines import (
    extract_continuum_linesubsets,
)
from .winding import detect_crossings_and_winding
from .plateau import _check_pmgbz_points_clustered, _probe_zero_plateau_near_mu1
from .mu2mid import Mu2MidZM


# Max bracket-expansion steps per side (aligned with amoeba's
# max_range_expansions=10).  The expansion loop only guards against runaway
# cases (unsolvable models, non-monotonic winding); a normal winding crosses
# zero within 1-2 steps of the default guess.
MAX_BRACKET_EXPANSIONS = 10


# ---------------------------------------------------------------------------
# W(E_ref, mu1) evaluation + continuum resolution
# ---------------------------------------------------------------------------

def _evaluate_winding(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    *,
    continuum_tol: float,
    crossing_tol: float,
) -> tuple[Optional[float], Optional[list], Mu2MidZM]:
    """Evaluate W(E_ref, mu1) and the 0D subsets at *mu1*.

    Builds a fresh ``Mu2MidZM`` at *mu1* (mu1 is the bisection variable, so
    each probe needs its own ZM) and runs :meth:`analyze` — the analysis
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
    zm.run()
    zm.analyze(tie_tol=continuum_tol, crossing_tol=crossing_tol)

    if zm.has_continuum:
        return None, None, zm

    subsets, W = detect_crossings_and_winding(
        zm, poly, crossing_tol=crossing_tol,
    )
    return W, subsets, zm


def _resolve_continuum_winding(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    *,
    continuum_perturb: float,
    continuum_tol: float,
    crossing_tol: float,
) -> tuple[Optional[float], Optional[float], Optional[float]]:
    """Compute the left / right winding limits at a continuum-degenerate mu1.

    Tries perturbation scales (1, 2, 4, 8) × *continuum_perturb* to escape
    the degenerate band.  Both sides are evaluated in sweep mode (the
    continuum gate short-circuits; the precise 0D subsets are not needed
    here — only the winding sign matters for the bisection).

    Returns ``(w_left, w_right, eps)`` where *eps* is the first scale at
    which BOTH sides escaped the degenerate band — the caller needs it to
    place the continuation point at ``mu1 ± eps``.  Returns ``(None, None,
    None)`` when all scales are still degenerate on some side.
    """
    for scale in core.ESCAPE_LADDER:
        eps = continuum_perturb * scale
        w_left, _, _ = _evaluate_winding(
            poly, E_ref, mu1 - eps,
            continuum_tol=continuum_tol, crossing_tol=crossing_tol,
        )
        w_right, _, _ = _evaluate_winding(
            poly, E_ref, mu1 + eps,
            continuum_tol=continuum_tol, crossing_tol=crossing_tol,
        )
        if (w_left is not None) and (w_right is not None):
            return w_left, w_right, eps
    return None, None, None


# ---------------------------------------------------------------------------
# μ₁ bisection (bracket expansion + plain midpoint with continuum interception)
# ---------------------------------------------------------------------------

@live_defaults(zero_tol="core:WINDING_ZERO_TOL", continuum_perturb="core:CONTINUUM_PERTURB",
    continuum_tol="core:CONTINUUM_TOL", crossing_tol="sgbz.pairwise:CROSSING_TOL",
    max_iter="sgbz.sgbz_solver:MU1_MAX_ITER")
def solve_SGBZ_for_E(
    poly: CharPoly,
    E_ref: complex,
    mu1_guess: tuple[float, float] = (-1, 1),
    *,
    zero_tol: Optional[float] = None,
    continuum_perturb: Optional[float] = None,
    max_iter: Optional[int] = None,
    continuum_tol: Optional[float] = None,
    crossing_tol: Optional[float] = None,
) -> dict:
    """Locate the winding-zero mu1 and return solve diagnostics.

    Uses bracket expansion + plain bisection (midpoint).  Each side of the
    bracket expansion is capped at ``MAX_BRACKET_EXPANSIONS`` steps (raise
    ``RuntimeError`` instead of looping forever on unsolvable / anomalous
    winding).  When a continuum-degenerate mu1 is encountered,
    ``_resolve_continuum_winding`` resolves the left/right winding limits;
    if they straddle zero that mu1 is the SGBZ boundary
    (``is_continuum=True``) and the 1D LineSubsets are materialized
    (``extract_continuum_linesubsets``) from the built ZM; the returned
    ``subsets`` carries them (the caller signals "in spectrum").

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
    # --- winding_at: evaluate W(E_ref, mu1) + 0D subsets at *mu1_val* ---
    # Partially-applied _evaluate_winding on the fixed (poly, E_ref,
    # tolerances): returns the full (winding, subsets, zm)
    # triple.  The subsets are used by the callers
    # below (gbz_low/gbz_high/gbz_mid/gbz_final), so this is the general
    # evaluation, not a winding-only cheap probe (that lives in
    # _resolve_continuum_winding, which discards the subsets).
    def winding_at(mu1_val: float):
        return _evaluate_winding(
            poly, E_ref, mu1_val,
            continuum_tol=continuum_tol, crossing_tol=crossing_tol,
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
        w_l, w_r, eps = _resolve_continuum_winding(
            poly, E_ref, mu1_val,
            continuum_perturb=continuum_perturb,
            continuum_tol=continuum_tol, crossing_tol=crossing_tol,
        )
        if w_l is None or w_r is None:
            raise ValueError(
                f"Cannot resolve continuum at mu1={mu1_val:.8g}: "
                f"all perturbation scales still degenerate."
            )
        # Boundary: the winding limits STRICTLY straddle zero (opposite
        # signs).  A limit exactly zero is NOT a boundary — it is a
        # zero-plateau signal (the winding vanishes right at the band edge),
        # and a zero plateau is not a GBZ point; that case falls through to
        # the same-sign branch, which returns the =0 edge so the bisection
        # converges there and the discrete path + plateau check handle it.
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
        # Non-boundary (same sign or a zero limit): W is monotonic increasing
        # in mu1, so the winding-zero is OUTSIDE (or at the edge of) the
        # degenerate band.  w_l == 0 means the zero sits at mu1−ε (left edge,
        # a zero plateau); w_r == 0 (with w_l < 0) means it sits at mu1+ε.
        # Return the band edge nearest the zero — INCLUDING the =0 edge —
        # paired with its own winding, so the bisection converges there and
        # the discrete path + plateau check reclassify it.  Do NOT bolt a
        # proxy winding onto the degenerate mu1_val.
        if w_l >= 0:
            return False, (mu1_val - eps, w_l), None   # zero at/left of band
        return False, (mu1_val + eps, w_r), None       # zero right of band

    # --- Step 1: bracket expansion ---
    # Initialize bracket variables before defining handle_continuum closure.
    # The closure reads these at call time; during expansion the right bracket
    # may not yet be established, so we initialize to None and set them later.
    mu1_low = mu1_guess[0]
    mu1_high = None
    w_low = None
    w_high = None
    mu1_ext_right = None

    for _ in range(MAX_BRACKET_EXPANSIONS):
        w_low, gbz_low, zm_low = winding_at(mu1_low)
        if w_low is None:
            is_boundary, proxy, result = handle_continuum(
                mu1_low, zm_low, (mu1_low, mu1_ext_right))
            if is_boundary:
                return result
            mu1_low, w_low = proxy  # corrected (mu1, w) at the band edge
            # Re-materialize the corrected band edge (mirrors the right-end
            # path below): the proxy carries only (mu1, winding), but the
            # left_endpoint_zero return below needs the 0D subsets — without
            # this re-evaluation they stay None (the continuum evaluation
            # returned (None, None, zm)) and collect_GBZ_subsets crashes on
            # them outside its try/except.
            w_low, gbz_low, zm_low = winding_at(mu1_low)
            # A second continuum at the band edge itself (perturbation still
            # inside the degenerate band): keep expanding rather than break
            # on a stale proxy.
            if w_low is None:
                continue

        if w_low < zero_tol:
            break
        mu1_ext_right = mu1_low
        mu1_low -= 1
    else:
        raise RuntimeError(
            f"left bracket expansion exceeded {MAX_BRACKET_EXPANSIONS} steps: "
            f"E_ref={E_ref}, mu1={mu1_low}, winding={w_low}"
        )

    if w_low > -zero_tol:
        return {
            "mu1": mu1_low,
            "subsets": gbz_low if gbz_low is not None else [],
            "winding": w_low,
            "is_continuum": False,
            "_mu1_bracket": (mu1_low, mu1_low),
            "_winding_bracket": (w_low, w_low),
            "_exit_reason": "left_endpoint_zero",
        }

    # Step 1.2: right bracket endpoint (unified for both entry paths: a
    # fresh guess and a mu1_ext_right already established during left-end
    # expansion).  A continuum proxy correction goes through the SAME
    # > -zero_tol check as a plain winding evaluation, so a corrected
    # same-sign / zero endpoint keeps expanding (or exits as
    # right_endpoint_zero) instead of silently entering bisection with a
    # non-straddling bracket.
    if mu1_ext_right is None:
        mu1_ext_right = mu1_guess[1]

    for _ in range(MAX_BRACKET_EXPANSIONS):
        w_high, gbz_high, zm_high = winding_at(mu1_ext_right)
        if w_high is None:
            is_boundary, proxy, result = handle_continuum(
                mu1_ext_right, zm_high, (mu1_low, mu1_ext_right))
            if is_boundary:
                return result
            mu1_ext_right, w_high = proxy
            # Re-evaluate the corrected band-edge point on the next
            # iteration: the proxy is already the nearest candidate, so the
            # generic "+1" step below must not skip past it.
            continue

        if w_high > -zero_tol:
            break
        mu1_ext_right += 1
    else:
        raise RuntimeError(
            f"right bracket expansion exceeded {MAX_BRACKET_EXPANSIONS} steps: "
            f"E_ref={E_ref}, mu1={mu1_ext_right}, winding={w_high}"
        )

    if abs(w_high) <= zero_tol:
        return {
            "mu1": mu1_ext_right,
            "subsets": gbz_high if gbz_high is not None else [],
            "winding": w_high,
            "is_continuum": False,
            "_mu1_bracket": (mu1_ext_right, mu1_ext_right),
            "_winding_bracket": (w_high, w_high),
            "_exit_reason": "right_endpoint_zero",
        }

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
            is_boundary, proxy, result = handle_continuum(
                mu1_mid, zm_mid, (mu1_low, mu1_high))
            if is_boundary:
                return result
            # At a non-boundary continuum, W is undefined at mu1_mid.  The
            # proxy is the (mu1, w) pair at the band edge nearest the zero;
            # use THAT mu1 as the bracket endpoint — do NOT converge on the
            # proxy value (the bisection target is the genuine W-zero, not a
            # continuum edge).
            mu1_adj, f_mid = proxy
            # A proxy winding of EXACTLY zero is the genuine answer, not a
            # bracket update: f_mid == 0 assigned to either side collapses
            # the bracket onto its own endpoint (w_low * 0 == 0 always takes
            # the else-branch) and the bisection then walks that degenerate
            # bracket until max_iter → guaranteed RuntimeError below.  The
            # proxy's band-edge point IS a W-zero; resolve its subsets and
            # return it.
            if f_mid == 0.0:
                w_edge, gbz_edge, _ = winding_at(mu1_adj)
                if w_edge is not None and abs(w_edge) <= zero_tol:
                    return {
                        "mu1": mu1_adj,
                        "subsets": gbz_edge if gbz_edge is not None else [],
                        "winding": w_edge,
                        "is_continuum": False,
                        "_mu1_bracket": (mu1_low, mu1_high),
                        "_winding_bracket": (w_low, w_high),
                        "_exit_reason": "w_zero_continuum_edge",
                    }
        else:
            mu1_adj = mu1_mid
            f_mid = w_mid
            # Convergence check only applies at genuine winding values.
            # At a continuum (w_mid is None) we continue bisection.
            # 2026-08-15: x_tol is no longer an exit reason. Bisection is finished only when:
            #   1. f_mid is zero, corresponding to the discrete case
            #   2. continuum is detected.
            if abs(f_mid) <= zero_tol:
                exit_reason = "w_zero"
                return {
                    "mu1": mu1_mid,
                    "subsets": gbz_mid,
                    "winding": w_mid,
                    "is_continuum": False,
                    "_mu1_bracket": (mu1_low, mu1_high),
                    "_winding_bracket": (w_low, w_high),
                    "_exit_reason": exit_reason,
                }

        # Bracket update (keep the diagnostics winding bracket honest — the
        # returned _winding_bracket should reflect the final bracket; f_mid
        # may be a continuum proxy, whose sign is all the bisection uses).
        if w_low * f_mid < 0:
            mu1_high = mu1_adj
            w_high = f_mid
        else:
            mu1_low = mu1_adj
            w_low = f_mid

    # Max iterations exhausted — bisection failed to converge.  Raise instead
    # of returning a mid-bracket point: the old fallback dressed a
    # non-convergence up as a "success" (a ±1-winding GBZ or an empty set)
    # and could silently misclassify the energy.  Report the failure with the
    # best-known point's full state.
    mu1_final = 0.5 * (mu1_low + mu1_high)
    w_final, gbz_final, _ = winding_at(mu1_final)
    is_continuum = (w_final is None)
    n_subsets = len(gbz_final) if gbz_final is not None else 0
    raise RuntimeError(
        f"bracket bisection failed to converge after {max_iter} iterations: "
        f"E_ref={E_ref}, mu1={mu1_final:.12g}, len(subsets)={n_subsets}, "
        f"is_continuum={is_continuum}, winding={w_final}"
    )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

@live_defaults(zero_tol="core:WINDING_ZERO_TOL", continuum_perturb="core:CONTINUUM_PERTURB",
    continuum_tol="core:CONTINUUM_TOL", crossing_tol="sgbz.pairwise:CROSSING_TOL",
    max_iter="sgbz.sgbz_solver:MU1_MAX_ITER")
def collect_GBZ_subsets(
    coeffs: np.ndarray,
    degs: np.ndarray,
    E_ref: complex,
    debug_mode: bool = False,
    *,
    plateau_check: bool = True,
    plateau_probe_radius: Optional[float] = None,
    mu1_guess: tuple[float, float] = (-1, 1),
    zero_tol: Optional[float] = None,
    continuum_perturb: Optional[float] = None,
    max_iter: Optional[int] = None,
    continuum_tol: Optional[float] = None,
    crossing_tol: Optional[float] = None,
) -> GBZResult:
    """Check the SGBZ condition and return GBZ points for a reference energy.

    Builds the characteristic polynomial from (coeffs, degs), solves for the
    mu1 where the average major-axis winding number vanishes using plain
    bisection with continuum interception, and (optionally) reclassifies the
    candidate as empty when a zero plateau exists right next to it.

    For continuum results (a 1D LineSubset case), the LineSubsets are
    materialised by ``extract_continuum_linesubsets`` from the built ZM and
    returned in ``subsets`` with ``index == (0, n_1d)`` (in spectrum).
    Spectrum membership is decided by the bisection's W-zero /
    left-right-limit straddle.

    Parameters:
        coeffs: complex coefficients of the characteristic Laurent polynomial
            f(E, beta1, beta2).
        degs: (n_terms, 3) integer exponents of (E, beta1, beta2) per term.
        E_ref: reference energy to test.
        perc: progress fraction in [0, 1], printed as a percentage.
        debug_mode: if True, re-raise solver exceptions instead of returning
            a failed GBZResult.
        Explicit keyword arguments only (no catch-all options dict): every
        tunable is named, defaults resolve from the home-module constants
        (see doc/constants.md), and a misspelled keyword raises TypeError
        instead of being silently ignored.

    Returns:
        GBZResult with connected subsets.  ``gbz.index == (0, 0)`` means
        E_ref is outside the SGBZ spectrum.  A 1D continuum result carries
        its materialized LineSubsets in ``subsets`` with
        ``index == (0, n_1d)``.
    """
    poly = CharPoly(coeffs, degs)

    try:
        sgbz_res = solve_SGBZ_for_E(
            poly, E_ref, mu1_guess=mu1_guess, zero_tol=zero_tol,
            continuum_perturb=continuum_perturb, max_iter=max_iter,
            continuum_tol=continuum_tol,
            crossing_tol=crossing_tol,
        )
        subsets = sgbz_res["subsets"]
        mu1 = sgbz_res["mu1"]
        is_continuum = sgbz_res["is_continuum"]

        # Continuum boundary → 1D LineSubsets materialised by handle_continuum
        # (the cluster tracks where |β_M| = |β_{M+1}| holds identically, joined
        # across MRs).  Spectrum membership is carried by index == (0, n_1d).
        if is_continuum:
            line_subsets = subsets or []
            n_1d = len(line_subsets)
            return GBZResult(E_ref=E_ref, success=True,
                             subsets=list(line_subsets), index=(0, n_1d))

        # Discrete case: apply the plateau check before trusting the subsets.
        if subsets and plateau_check:
            candidate = GBZResult(E_ref=E_ref, subsets=subsets,
                                index=(len(subsets), 0))
            if _check_pmgbz_points_clustered(candidate):
                plateau_info = _probe_zero_plateau_near_mu1(
                    poly, E_ref, mu1, sgbz_res.get("_mu1_bracket"),
                    continuum_tol=continuum_tol,
                    crossing_tol=crossing_tol,
                    zero_tol=zero_tol,
                    continuum_perturb=continuum_perturb,
                    probe_radius=plateau_probe_radius,
                )
                if plateau_info["found"]:
                    subsets = []

    except Exception as e:
        if debug_mode:
            raise e
        print("Error: %s" % str(e))
        return GBZResult(E_ref=E_ref, success=False, error=str(e))

    n_0d = sum(1 for s in subsets if isinstance(s, PointSubset))
    n_1d = sum(1 for s in subsets if isinstance(s, LineSubset))
    return GBZResult(E_ref=E_ref, subsets=list(subsets), index=(n_0d, n_1d))
