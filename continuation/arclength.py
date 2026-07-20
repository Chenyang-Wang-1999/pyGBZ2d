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
from typing import Optional, NamedTuple

from gbz_types import (
    CharPoly,
    hungarian_match_indices,
    sort_by_root_abs,
    to_sphere_r3,
    cost_from_sphere_r3,
)

# ---------------------------------------------------------------------------
# Step-size control constants (mirror scipy's RK45)
# ---------------------------------------------------------------------------

SAFETY = 0.9
MIN_FACTOR = 0.2
MAX_FACTOR = 10.0
ERROR_EXPONENT = -0.5  # -1/(p+1) for first-order tangent prediction (p=1)

# Threshold for treating a root as "effectively 0 or ∞"
ZERO_THRESHOLD = 1e-14
INF_THRESHOLD = 1e14


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


class MultipleRootInfo(NamedTuple):
    """Information about a detected multiple root."""
    theta1: float
    cluster_indices: np.ndarray  # indices into the modulus-sorted root array


class SegmentResult(NamedTuple):
    """Result of integrating one curve segment until a stop condition."""
    theta1_arr: np.ndarray      # (N,) θ₁ mesh
    tracked_roots: np.ndarray   # (N, K) Hungarian-tracked within segment
    abs_argsort: np.ndarray     # (N, K) per-row |β₂| argsort
    stop_reason: str            # 'completed' | 'multiple_root'
    mr_approx_theta: float      # approximate MR θ₁ (valid if stop='multiple_root')
    n_steps: int
    n_rejected: int


# ---------------------------------------------------------------------------
# Core step functions (unchanged)
# ---------------------------------------------------------------------------

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
        abs_b2 = np.abs(beta2)
        if abs_b2 < ZERO_THRESHOLD or abs_b2 > INF_THRESHOLD or not np.isfinite(beta2):
            continue

        partials = poly.eval_partials((E_ref, beta1, beta2))
        df_dbeta1 = partials[1]
        df_dbeta2 = partials[2]

        dbeta2_dtheta1 = -1j * beta1 * df_dbeta1 / df_dbeta2
        V_j = dbeta2_dtheta1 / beta2

        if not np.isfinite(V_j):
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
        abs_b2 = np.abs(beta2)
        if abs_b2 < ZERO_THRESHOLD or abs_b2 > INF_THRESHOLD or not np.isfinite(beta2):
            predicted[j] = beta2
            continue
        if V[j] == 0:
            predicted[j] = beta2
            continue
        predicted[j] = beta2 * np.exp(V[j] * dtheta1)

    return predicted


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
    for from_idx, to_idx in matches:
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
    direction: float = 1.0,
    max_step: float = 0.5,
    min_step: float = 1e-14,
    atol: float = 1e-12,
    rtol: float = 1e-3,
    safety: float = SAFETY,
    min_factor: float = MIN_FACTOR,
    max_factor: float = MAX_FACTOR,
    error_exponent: float = ERROR_EXPONENT,
    max_iter: int = 20,
) -> StepResult:
    """Take one adaptive pseudo-arclength step along θ₁.

    1. Compute tangent, propose dθ₁ = direction * h / ‖V‖₂.
    2. Predict roots via tangent extrapolation.
    3. Solve actual roots via np.roots.
    4. Accept if error_norm < 1; otherwise reduce h and retry.
    """
    V, norm_V = compute_tangent(poly, E_ref,
                                 exp(mu1 + 1j * theta1), roots)

    step_rejected = False

    for _ in range(max_iter):
        if h < min_step:
            return StepResult(theta1, roots, h, False, V, norm_V)

        dtheta1 = direction * h / norm_V
        theta1_new = theta1 + dtheta1

        predicted = predict_roots(roots, V, dtheta1)

        beta1_new = exp(mu1 + 1j * theta1_new)
        roots_new_list = poly.solve_roots_1d(
            (0, 1), (E_ref, beta1_new), (2,),
        )
        roots_new = np.asarray(roots_new_list, dtype=complex)

        error_norm = estimate_error(predicted, roots_new, atol, rtol)

        if error_norm < 1.0:
            if error_norm == 0.0:
                factor = max_factor
            else:
                factor = min(max_factor, safety * error_norm ** error_exponent)
            if step_rejected:
                factor = min(1.0, factor)
            h_new = h * factor
            if h_new > max_step:
                h_new = max_step
            return StepResult(theta1_new, roots_new, h_new, True, V, norm_V)
        else:
            h *= max(min_factor, safety * error_norm ** error_exponent)
            step_rejected = True

    return StepResult(theta1, roots, h, False, V, norm_V)


# ---------------------------------------------------------------------------
# Multiple-root detection & refinement
# ---------------------------------------------------------------------------

def detect_possible_multiple_root(
    roots: np.ndarray,
    dtheta: float,
    *,
    min_dtheta: float = 1e-10,
    cluster_tol: float = 1e-6,
) -> Optional[np.ndarray]:
    """Check whether a multiple root is likely at the current position.

    Returns cluster indices if ``dtheta < min_dtheta`` AND a root cluster
    is detected, or None otherwise.  This is the single entry point for
    MR detection during integration — extend here to add heuristics.
    """
    if dtheta >= min_dtheta:
        return None
    return detect_cluster(roots, cluster_tol=cluster_tol)


def detect_cluster(
    roots: np.ndarray,
    *,
    cluster_tol: float = 1e-6,
) -> Optional[np.ndarray]:
    """Detect a root cluster by pairwise chordal distance.

    Parameters
    ----------
    roots : np.ndarray, shape (n_roots,)
        β₂ roots (any ordering).
    cluster_tol : float
        Chordal-distance threshold.

    Returns
    -------
    np.ndarray or None
        Indices of roots belonging to a cluster, or None if all roots are
        well-separated.
    """
    r3 = to_sphere_r3(roots)
    pw = cost_from_sphere_r3(r3, r3)
    np.fill_diagonal(pw, np.inf)

    if np.min(pw) > cluster_tol:
        return None

    in_cluster = np.any(pw < cluster_tol, axis=1)
    indices = np.where(in_cluster)[0]
    return indices if len(indices) > 0 else None


def refine_multiple_root_theta(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    theta_guess: float,
    *,
    cluster_tol: float = 1e-6,
    max_iter: int = 40,
) -> float:
    """Refine the θ₁ of a multiple root via ternary search.

    Minimises the minimum pairwise chordal distance of β₂ roots as a
    function of θ₁ near *theta_guess*.

    Returns the refined θ₁ (in [0, 2π)).
    """
    def _min_pairwise_dist(t: float) -> float:
        beta1 = exp(mu1 + 1j * (t % (2 * pi)))
        roots_list = poly.solve_roots_1d((0, 1), (E_ref, beta1), (2,))
        roots = np.asarray(roots_list, dtype=complex)
        r3 = to_sphere_r3(roots)
        pw = cost_from_sphere_r3(r3, r3)
        np.fill_diagonal(pw, np.inf)
        return float(np.min(pw))

    f_mid = _min_pairwise_dist(theta_guess)
    if f_mid < 1e-12:
        return theta_guess % (2 * pi)

    delta = 0.05
    a = theta_guess - delta
    b = theta_guess + delta

    for _ in range(max_iter):
        if b - a < 1e-12:
            break
        m1 = a + (b - a) / 3.0
        m2 = b - (b - a) / 3.0
        if _min_pairwise_dist(m1) < _min_pairwise_dist(m2):
            b = m2
        else:
            a = m1

    result = (0.5 * (a + b)) % (2 * pi)
    # Snap to 0 if the minimum is very close to θ₁=0 ≡ 2π.
    if min(result, 2 * pi - result) < 1e-6:
        result = 0.0
    return result


# ---------------------------------------------------------------------------
# Segment integration
# ---------------------------------------------------------------------------

def integrate_segment(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    theta_start: float,
    roots_start: np.ndarray,
    theta_end: float,
    *,
    h0: float = 0.1,
    max_step: float = 0.5,
    min_step: float = 1e-12,
    min_dtheta: float = 1e-10,
    cluster_tol: float = 1e-6,
    atol: float = 1e-12,
    rtol: float = 1e-3,
    safety: float = SAFETY,
    min_factor: float = MIN_FACTOR,
    max_factor: float = MAX_FACTOR,
    error_exponent: float = ERROR_EXPONENT,
    max_step_iter: int = 20,
) -> SegmentResult:
    """Integrate β₂ roots from *theta_start* toward *theta_end*.

    Stops when *theta_end* is reached (``stop_reason='completed'``) or when
    a multiple root is detected (``stop_reason='multiple_root'``).

    *roots_start* must be track-ordered (column j = physical root j).
    The returned ``tracked_roots`` does **not** include the start point —
    the caller already has it.  Each subsequent row is Hungarian-matched
    to maintain track continuity within the segment.
    """
    n_roots = len(roots_start)

    theta1_list: list[float] = []
    tracked_list: list[np.ndarray] = []

    h = h0
    theta1 = theta_start
    roots = roots_start

    n_accepted = 0
    n_rejected = 0

    while theta1 < theta_end - 1e-15:
        result = arclength_step(
            poly, E_ref, mu1, theta1, roots, h,
            direction=1.0,
            max_step=max_step,
            min_step=min_step,
            atol=atol,
            rtol=rtol,
            safety=safety,
            min_factor=min_factor,
            max_factor=max_factor,
            error_exponent=error_exponent,
            max_iter=max_step_iter,
        )

        if result.accepted:
            dtheta = abs(result.theta1_new - theta1)

            cluster = detect_possible_multiple_root(
                roots, dtheta,
                min_dtheta=min_dtheta, cluster_tol=cluster_tol,
            )
            if cluster is not None:
                theta1_arr = np.array(theta1_list)
                tracked_roots_arr = np.array(tracked_list)
                abs_argsort = np.argsort(np.abs(tracked_roots_arr), axis=1) if len(tracked_list) > 0 else np.empty((0, n_roots), dtype=int)
                return SegmentResult(
                    theta1_arr=theta1_arr,
                    tracked_roots=tracked_roots_arr,
                    abs_argsort=abs_argsort,
                    stop_reason='multiple_root',
                    mr_approx_theta=theta1,
                    n_steps=n_accepted,
                    n_rejected=n_rejected,
                )

            theta1_new = result.theta1_new
            h = result.h_new
            roots_new = result.roots_new

            # Clamp to theta_end.
            if theta1_new > theta_end:
                theta1_new = theta_end
                beta1_end = exp(mu1 + 1j * theta_end)
                roots_end_list = poly.solve_roots_1d(
                    (0, 1), (E_ref, beta1_end), (2,),
                )
                roots_new = np.asarray(roots_end_list, dtype=complex)

            # Hungarian-match to maintain track continuity.
            matches = hungarian_match_indices(roots, roots_new)
            reordered = np.zeros(n_roots, dtype=complex)
            for from_idx, to_idx in matches:
                reordered[from_idx] = roots_new[to_idx]

            theta1_list.append(theta1_new)
            tracked_list.append(reordered)

            theta1 = theta1_new
            roots = reordered
            n_accepted += 1

        else:
            # Step rejected with h < min_step.
            n_rejected += 1
            cluster = detect_cluster(roots, cluster_tol=cluster_tol)
            if cluster is not None:
                theta1_arr = np.array(theta1_list)
                tracked_roots_arr = np.array(tracked_list)
                abs_argsort = np.argsort(np.abs(tracked_roots_arr), axis=1) if len(tracked_list) > 0 else np.empty((0, n_roots), dtype=int)
                return SegmentResult(
                    theta1_arr=theta1_arr,
                    tracked_roots=tracked_roots_arr,
                    abs_argsort=abs_argsort,
                    stop_reason='multiple_root',
                    mr_approx_theta=theta1,
                    n_steps=n_accepted,
                    n_rejected=n_rejected,
                )
            # Steep but no cluster: accept with a tiny fixed step.
            dtheta_jump = min_step * 100
            theta1_new = min(theta1 + dtheta_jump, theta_end)
            beta1_new = exp(mu1 + 1j * theta1_new)
            roots_new_list = poly.solve_roots_1d(
                (0, 1), (E_ref, beta1_new), (2,),
            )
            roots_new = np.asarray(roots_new_list, dtype=complex)
            matches = hungarian_match_indices(roots, roots_new)
            reordered = np.zeros(n_roots, dtype=complex)
            for from_idx, to_idx in matches:
                reordered[from_idx] = roots_new[to_idx]
            theta1_list.append(theta1_new)
            tracked_list.append(reordered)
            theta1 = theta1_new
            roots = reordered
            n_accepted += 1
            h = h0

    # Completed without hitting a multiple root.
    theta1_arr = np.array(theta1_list)
    tracked_roots_arr = np.array(tracked_list)
    abs_argsort = np.argsort(np.abs(tracked_roots_arr), axis=1) if len(tracked_list) > 0 else np.empty((0, n_roots), dtype=int)
    return SegmentResult(
        theta1_arr=theta1_arr,
        tracked_roots=tracked_roots_arr,
        abs_argsort=abs_argsort,
        stop_reason='completed',
        mr_approx_theta=float('nan'),
        n_steps=n_accepted,
        n_rejected=n_rejected,
    )
