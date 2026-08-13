"""
Pseudo-arclength continuation for beta2 root tracking along theta1.

Uses analytic derivatives from CharPoly.eval_partials to compute the tangent
vector V = [1, d(ln β₂₁)/dθ₁, ..., d(ln β₂ₙ)/dθ₁]ᵀ, then adapts the step
size Δθ₁ = h / ‖V‖₂ so that the arclength step in (θ₁, ln β₂)-space is
approximately constant.

Adaptive step-size control follows the pattern of scipy's RK45 integrator:
error between tangent-predicted roots and np.roots-actual roots drives the PI
step-size controller with SAFETY, MIN_FACTOR, MAX_FACTOR, and error_exponent.

References
----------
- E. Hairer, S. P. Norsett, G. Wanner, "Solving Ordinary Differential
  Equations I: Nonstiff Problems", Sec. II.4.
"""

from __future__ import annotations

import numpy as np
from math import pi
from cmath import exp
from dataclasses import dataclass
from typing import Optional, NamedTuple

from gbz_types import (
    CharPoly,
    hungarian_match_indices,
    to_sphere_r3,
)

# ---------------------------------------------------------------------------
# Step-size control constants (mirror scipy's RK45)
# ---------------------------------------------------------------------------

SAFETY = 0.9
MIN_FACTOR = 0.2
MAX_FACTOR = 10.0
ERROR_EXPONENT = -0.5  # -1/(p+1) for first-order tangent prediction (p=1)

# Threshold for treating a root as "effectively 0 or ∞"
ZERO_THRESHOLD = 1e-6
INF_THRESHOLD = 1e6


@dataclass(frozen=True)
class StepControl:
    """Tolerances and factors for the adaptive pseudo-arclength step controller.

    Bundles the RK45-style PI controller knobs (SAFETY / MIN_FACTOR /
    MAX_FACTOR / ERROR_EXPONENT) and the arclength step bounds so that
    ``arclength_step``, ``integrate_segment`` and ``ZeroManager.run`` don't
    each re-declare nine parameters.  Add a knob here once; all three layers
    carry the same ``StepControl`` instance.
    """
    max_step: float = 0.5
    min_step: float = 1e-12
    atol: float = 1e-12
    rtol: float = 1e-3
    safety: float = SAFETY
    min_factor: float = MIN_FACTOR
    max_factor: float = MAX_FACTOR
    error_exponent: float = ERROR_EXPONENT
    max_iter: int = 20


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------

class StepResult(NamedTuple):
    """Result of one adaptive pseudo-arclength step."""
    theta1_new: float
    roots_new: np.ndarray       # (n_roots,) solved roots
    h_new: float
    accepted: bool
    V: np.ndarray               # tangent vector at the old point
    norm_V: float               # ‖V‖₂ at the old point

# ---------------------------------------------------------------------------
# Core step functions
# ---------------------------------------------------------------------------

def _is_singular_root(beta2: complex) -> bool:
    """True for the 0 / ∞ padding roots appended by ``CharPoly.solve_roots_1d``
    when the 1-D polynomial is degree-deficient.  Their tangent is undefined,
    so they are held fixed during prediction.
    """
    abs_b2 = np.abs(beta2)
    return (abs_b2 < ZERO_THRESHOLD or abs_b2 > INF_THRESHOLD
            or not np.isfinite(beta2))


def compute_tangent(
    poly: CharPoly,
    E_ref: complex,
    beta1: complex,
    roots: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Compute tangent vector V_j = d(ln β₂ⱼ)/dθ₁ for each root.

    Roots at 0 or ∞ get V_j = 0.  When ∂f/∂β₂ ≈ 0 the implicit-function
    derivative diverges naturally; no artificial cap is applied — the
    caller's min_dtheta / min_step guards catch the divergence.
    """
    n_roots = len(roots)
    V = np.zeros(n_roots, dtype=complex)

    for j in range(n_roots):
        beta2 = roots[j]
        if _is_singular_root(beta2):
            continue

        partials = poly.eval_partials((E_ref, beta1, beta2))
        df_dbeta1 = partials[1]
        df_dbeta2 = partials[2]

        # At a multiple root ∂f/∂β₂ = 0, dβ₂/dθ₁ diverges.  Set V_j = ∞
        # so downstream code can distinguish "undefined" from "truly zero".
        if df_dbeta2 == 0:
            V[j] = np.inf + 0j
            continue

        dbeta2_dtheta1 = -1j * beta1 * df_dbeta1 / df_dbeta2
        V_j = dbeta2_dtheta1 / beta2

        if not np.isfinite(V_j):
            V[j] = np.inf + 0j
            continue

        V[j] = V_j

    norm_V = np.sqrt(1.0 + np.sum(np.abs(V) ** 2))
    return V, norm_V


def predict_roots(
    roots: np.ndarray,
    V: np.ndarray,
    dtheta1: float,
) -> np.ndarray:
    """First-order tangent prediction: β₂ⱼ → β₂ⱼ · exp(Vⱼ · Δθ₁).

    Roots at 0 or ∞ (where Vⱼ = 0) stay unchanged.
    """
    predicted = np.empty_like(roots)

    for j in range(len(roots)):
        beta2 = roots[j]
        if _is_singular_root(beta2):
            predicted[j] = beta2
            continue
        if not np.isfinite(V[j]):
            predicted[j] = beta2
            continue
        predicted[j] = beta2 * np.exp(V[j] * dtheta1)

    return predicted


def predict_roots_hermite(
    theta_target: float,
    theta0: float, roots0: np.ndarray, V0: np.ndarray,
    theta1: Optional[float] = None,
    roots1: Optional[np.ndarray] = None,
    V1: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Predict β₂ roots at *theta_target* from one or two reference endpoints.

    This is the matching anchor used outside the arclength integrator (in
    ``ZeroManager`` and the multiple-root solvers): instead of bare
    Hungarian matching between two solved root sets — which swaps two
    near-degenerate tracks at a closest approach — we predict where each
    track should land and match *predicted → solved*, exactly as
    :func:`arclength_step` does internally.

    Two-endpoint (cubic Hermite): given the value *and* tangent at both
    ends of an interval whose endpoints share a track frame, the per-track
    cubic Hermite polynomial is unique.  Use when *theta_target* lies inside
    such an interval (e.g. a bracketed multiple root, or θ₁ = 2π between two
    segment rows).

    One-endpoint (tangent extrapolation): falls back to the same
    ``β₂·exp(V·Δθ)`` prediction as :func:`predict_roots`, consistent with the
    arclength integrator.  Use when only one reference point is available
    (e.g. ``multiple_root_encountered``).

    Per-track fallback (two-endpoint): Hermite → two-point linear
    interpolation → hold fixed.  Singular roots (the 0/∞ padding roots
    ``solve_roots_1d`` appends, or a divergent tangent) step down the chain.
    One-endpoint singular roots are held fixed by :func:`predict_roots`.
    """
    if theta1 is not None and roots1 is not None and V1 is not None:
        # Degenerate interval (coincident θ): no interpolation possible.
        if abs(theta1 - theta0) < 1e-15:
            return roots0.copy()

        s = (theta_target - theta0) / (theta1 - theta0)
        dt = theta1 - theta0
        h00 = 2 * s ** 3 - 3 * s ** 2 + 1
        h10 = s ** 3 - 2 * s ** 2 + s
        h01 = -2 * s ** 3 + 3 * s ** 2
        h11 = s ** 3 - s ** 2

        predicted = np.empty_like(roots0)
        for j in range(len(roots0)):
            p0 = roots0[j]
            p1 = roots1[j]
            if _is_singular_root(p0) or _is_singular_root(p1):
                predicted[j] = p0
                continue
            m0 = V0[j] * p0
            m1 = V1[j] * p1
            cand = h00 * p0 + h10 * dt * m0 + h01 * p1 + h11 * dt * m1
            if not np.isfinite(cand) or _is_singular_root(cand):
                # Hermite diverged → two-point linear interpolation.
                cand = p0 + s * (p1 - p0)
                if not np.isfinite(cand) or _is_singular_root(cand):
                    cand = p0
            predicted[j] = cand
        return predicted

    # One-endpoint tangent extrapolation.
    return predict_roots(roots0, V0, theta_target - theta0)


def estimate_error(
    roots_predicted: np.ndarray,
    roots_actual: np.ndarray,
    atol: float = 1e-12,
    rtol: float = 1e-3,
) -> float:
    """Compute error norm between tangent-predicted and np.roots-actual roots.

    Hungarian-matches predicted → actual via chordal distance on the Riemann
    sphere, then returns max chordal distance scaled by typical root magnitude.
    < 1 means the step is acceptable.
    """
    matches = hungarian_match_indices(roots_predicted, roots_actual)
    p1_sphere = to_sphere_r3(roots_predicted)
    p2_sphere = to_sphere_r3(roots_actual)

    max_chordal = 0.0
    for from_idx, to_idx in enumerate(matches):
        dist = np.linalg.norm(p1_sphere[from_idx] - p2_sphere[to_idx])
        if dist > max_chordal:
            max_chordal = dist

    scale = atol + rtol * np.median(np.abs(roots_actual))
    return max_chordal / scale


def arclength_step(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    theta1: float,
    roots: np.ndarray,
    h: float,
    ctrl: StepControl = StepControl(),
) -> StepResult:
    """Take one adaptive pseudo-arclength step along θ₁.

    1. Compute tangent, propose dθ₁ = h / ‖V‖₂.
    2. Predict roots via tangent extrapolation.
    3. Solve actual roots via np.roots.
    4. Accept if error_norm < 1; otherwise reduce h and retry.

    *ctrl* carries the RK45-style tolerances and step bounds; see
    :class:`StepControl`.
    """
    V, norm_V = compute_tangent(poly, E_ref,
                                 exp(mu1 + 1j * theta1), roots)

    step_rejected = False

    for _ in range(ctrl.max_iter):
        if h < ctrl.min_step:
            return StepResult(theta1, roots, h, False, V, norm_V)

        dtheta1 = h / norm_V
        theta1_new = theta1 + dtheta1

        predicted = predict_roots(roots, V, dtheta1)

        beta1_new = exp(mu1 + 1j * theta1_new)
        roots_new = poly.solve_roots_1d(
            (0, 1), (E_ref, beta1_new), (2,),
        )

        error_norm = estimate_error(predicted, roots_new, ctrl.atol, ctrl.rtol)

        if error_norm < 1.0:
            if error_norm == 0.0:
                factor = ctrl.max_factor
            else:
                factor = min(ctrl.max_factor,
                             ctrl.safety * error_norm ** ctrl.error_exponent)
            if step_rejected:
                factor = min(1.0, factor)
            h_new = h * factor
            if h_new > ctrl.max_step:
                h_new = ctrl.max_step
            # Return roots_new in track order, anchored on the tangent
            # prediction.  solve_roots_1d returns np.roots order (unstable
            # across theta1), so raw chordal matching old→new would swap two
            # near-degenerate tracks at a closest approach.  
            # So predicted→new is the right anchor for the identity-preserving permutation.
            matches = hungarian_match_indices(predicted, roots_new)
            roots_new = roots_new[matches]
            return StepResult(theta1_new, roots_new, h_new, True, V, norm_V)
        else:
            h *= max(ctrl.min_factor,
                     ctrl.safety * error_norm ** ctrl.error_exponent)
            step_rejected = True

    return StepResult(theta1, roots, h, False, V, norm_V)
