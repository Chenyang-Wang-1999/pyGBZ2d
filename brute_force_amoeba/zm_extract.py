'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-24
Copyright © Department of Physics, Tsinghua University. All rights reserved

ZeroManager-based amoeba subset extraction and winding.

Defines ``AmoebaZeroManager`` (a ``ZeroManager`` specialization caching
per-segment ``ln|β₂|``) and the subset-extraction / winding functions that
operate on its segments.

Replaces the old tracks-based ``_compute_winding_from_tracks`` /
``_extract_continuum_intervals`` / ``_detect_crossings_outside_continuum``
layer.  Operates directly on ``AmoebaZeroManager`` segments (adaptive β₂-root
tracks that correctly cross multiple roots) with cached ``ln|β₂|``.

Two consumers:
  - ``extract_amoeba_subsets`` : GBZ subset output (LineSubset / PointSubset).
    Applies boundary dedup (Rule 1 continuum-endpoint snap + Rule 2 d==0
    θ₁-dedup) so shared segment endpoints (an MR) do not produce duplicates.
  - ``amoeba_windings`` : w2 winding (and zeros) for the μ₂-bisection.
    Uses RAW crossing detection — no boundary filtering (the winding integral
    needs every crossing).
'''

from __future__ import annotations

from cmath import exp, pi
from typing import Literal, Optional

import numpy as np

from gbz_types import CharPoly, PointSubset, LineSubset

from .ronkin_winding import _find_exact_crossing, _get_average_winding_from_zeros
from continuation.zero_manager import SegmentData, ZeroManager


# A genuine continuum is constant-modulus to ~1e-9; a transversal crossing
# leaves the level after one sample.  1e-6 separates them.
CONTINUUM_TOL = 1e-6
CONTINUUM_FRAC = 0.9
# A crossing whose θ₁ lands within this of a continuum LineSubset endpoint is
# treated as the continuum/MR boundary, not a genuine discrete zero, and
# dropped.  Decoupled from CONTINUUM_TOL so the snap radius is independent of
# band-detection sensitivity.
SNAP_TOL = 1e-3

ExtractMode = Literal['coarse', 'fine', 'solve']


# ---------------------------------------------------------------------------
# AmoebaZeroManager — ZeroManager with cached per-segment ln|β₂|
# ---------------------------------------------------------------------------

class AmoebaZeroManager(ZeroManager):
    """``ZeroManager`` with per-segment ``ln|β₂|`` cached for amoeba use.

    The amoeba pipeline evaluates the same root tracks at many μ₂ values
    (subset extraction + ~30 μ₂ bisection steps).  Caching ``log|tracked_roots|``
    per segment here lets the extractor and winding code index a precomputed
    array instead of recomputing logs on every call.

    After ``.run()``:
      - ``seg_logabs`` : ``list[np.ndarray]``, ``seg_logabs[s]`` has shape
        ``(N_s, K)`` = ``np.log(np.abs(segments[s].tracked_roots))``.
    """

    seg_logabs: list[np.ndarray]

    def __init__(self, poly: CharPoly, E_ref: complex, mu1: float):
        super().__init__(poly, E_ref, mu1)
        self.seg_logabs = []

    def run(self, **kwargs) -> None:
        super().run(**kwargs)
        # Cache ln|β₂| per segment.  Re-run() resets segments, so recompute
        # the cache unconditionally — do not append to a stale list.
        self.seg_logabs = [
            np.log(np.abs(seg.tracked_roots)) for seg in self.segments
        ]


# ---------------------------------------------------------------------------
# Shared detection
# ---------------------------------------------------------------------------

def _continuum_mask(zm: AmoebaZeroManager, mu2: float, tol: float, frac: float) -> list[np.ndarray]:
    """Per-segment boolean mask of continuum tracks (frac > threshold).

    Returns ``list[np.ndarray]`` of shape ``(K,)`` per segment; ``mask[s][j]``
    is True when track ``j`` of segment ``s`` stays at ``ln|β₂| = mu2`` for a
    fraction ``> frac`` of the segment samples.
    """
    masks: list[np.ndarray] = []
    for la in zm.seg_logabs:
        frac_in_band = np.mean(np.abs(la - mu2) < tol, axis=0)
        masks.append(frac_in_band > frac)
    return masks


# ---------------------------------------------------------------------------
# Subset extraction (3 modes)
# ---------------------------------------------------------------------------

def extract_amoeba_subsets(
    zm: AmoebaZeroManager,
    poly: CharPoly,
    E: complex,
    mu1: float,
    mu2: float,
    *,
    mode: ExtractMode = 'solve',
    tol: float = CONTINUUM_TOL,
    frac: float = CONTINUUM_FRAC,
    snap_tol: float = SNAP_TOL,
) -> list:
    """GBZ subsets of f at ln|b1|=mu1, ln|b2|=mu2, from ZM tracks.

    LineSubset handling is identical across modes: a continuum track (frac >
    ``frac``) yields one ``LineSubset`` over the full segment — detected and
    returned without refinement.

    PointSubset handling depends on ``mode``:
      - ``'coarse'`` : linear-interpolation crossing, no refinement.
      - ``'fine'``   : one ``_find_exact_crossing`` refinement step (fsolve +
        analytic Jacobian), falling back to linear on failure.
      - ``'solve'``  : same as ``fine``; the mode used by ``collect_GBZ_subsets``
        for the full GBZ output.

    Boundary dedup (after snap-to-mean, MR endpoints are exactly degenerate):
      - Rule 1 — drop any crossing within ``snap_tol`` of a continuum LineSubset
        endpoint (continuum/MR edge, not a genuine discrete zero).
      - Rule 2 — ``d == 0`` exact touches are deduplicated by zero-point
        identity, not by θ₁.  A zero's identity is ``(endpoint, track)``:
        an interior touch gets a unique key; a touch on a shared MR endpoint
        collides with the same MR + track on the neighbouring segment; a
        touch on the θ₁=0≡2π circle boundary collides with its counterpart
        on the other side via ``zm.boundary_perm``.  Distinct tracks at the
        same θ₁ (e.g. a conjugate pair at β₁=1) thus stay distinct.
    """
    continuum_masks = _continuum_mask(zm, mu2, tol, frac)
    continuum_lines: list[LineSubset] = []
    # Crossings recorded as (seg_idx, i, j, kind), NOT PointSubsets, so that
    # boundary duplicates between segments sharing an MR endpoint can be
    # resolved after the loop.  kind: 'zero' (d==0 at sample i) | 'cross'
    # (sign change over [i, i+1]).
    hits: list[tuple[int, int, int, str]] = []

    for s, seg in enumerate(zm.segments):
        th = seg.theta1_arr
        tr = seg.tracked_roots
        is_cont = continuum_masks[s]

        # 1. continuum tracks → LineSubsets.
        for j in np.where(is_cont)[0]:
            continuum_lines.append(LineSubset(
                E=E, mu1=mu1,
                theta1_arr=th.copy(), beta2_arr=tr[:, j].copy()))

        # 2. discrete crossings on non-continuum tracks (index only).
        other = np.where(~is_cont)[0]
        if other.size:
            la = zm.seg_logabs[s]
            d = la[:, other] - mu2                      # (N, J)

            # d == 0 cannot trigger the sign-change product below, record
            # separately — it is an exact numerical touch at a sample.
            i_idx, j_idx = np.where(d == 0)
            for i, jj in zip(i_idx, j_idx):
                hits.append((s, int(i), int(other[jj]), 'zero'))

            sc = d[:-1] * d[1:] < 0                     # (N-1, J)
            i_idx, j_idx = np.where(sc)
            for i, jj in zip(i_idx, j_idx):
                hits.append((s, int(i), int(other[jj]), 'cross'))

    # ---- materialize + boundary filtering ----
    cont_endpoints = np.array(
        [t for L in continuum_lines for t in (L.theta1_start, L.theta1_end)]
    ) if continuum_lines else None

    # Inverse of zm.boundary_perm: boundary_perm[inv[k]] == k, so
    # roots_right[j] corresponds to left-boundary track inv[j].
    # Used to fold the θ₁=2π end of the last segment onto the θ₁=0 start of
    # the first segment — they are the same physical circle point.
    boundary_perm_inv = np.empty(zm.K, dtype=int)
    boundary_perm_inv[zm.boundary_perm] = np.arange(zm.K)

    seen_zero: set[tuple] = set()
    points: list[tuple[float, complex]] = []
    for s, i, j, kind in hits:
        seg = zm.segments[s]
        th, tr = seg.theta1_arr, seg.tracked_roots
        la = zm.seg_logabs[s]
        t1, b2 = _finalize_crossing(kind, i, j, th, tr, la, mu2, mode, poly, E, mu1)

        # Rule 1: continuum-endpoint snap.
        if cont_endpoints is not None and np.min(np.abs(cont_endpoints - t1)) < snap_tol:
            continue

        # Rule 2: dedup exact-touch ('zero') hits by zero identity.
        # 'cross' hits are interior sign-changes and are always unique.
        if kind == 'zero':
            key = _zero_identity_key(zm, seg, i, j, boundary_perm_inv)
            if key in seen_zero:
                continue
            seen_zero.add(key)

        points.append((t1, b2))

    subsets: list = list(continuum_lines)
    for t1, b2 in points:
        subsets.append(PointSubset(E=E, beta1=exp(mu1 + 1j * t1), beta2=b2))
    return subsets


def _zero_identity_key(
    zm: AmoebaZeroManager,
    seg: SegmentData,
    i: int,
    j: int,
    boundary_perm_inv: np.ndarray,
) -> tuple:
    """Identity key for an exact-touch (``d == 0``) zero at ``(seg, i, j)``.

    Two hits that are the *same physical zero* must collide here so Rule 2
    keeps only one.  Identity is ``(endpoint, track)``, never θ₁:

    * **Interior sample** — unique, never deduped.
    * **Shared MR endpoint** — a segment's right boundary row equals the
      next segment's left boundary row (both are ``multiple_roots[R].roots``,
      prepended in ``ZeroManager.run``); key ``('mr', R, j)`` collapses them.
    * **Circle boundary θ₁=0≡2π** — the first segment's left end and the
      last segment's right end are the same point, but in different track
      orders; ``boundary_perm_inv`` folds the right end onto the left track,
      key ``('circle', j_left)``.
    """
    n = len(seg.theta1_arr)
    is_left = (i == 0)
    is_right = (i == n - 1)

    # Circle boundary θ₁ = 0 (first segment's left end).
    # seg0.left_mr == -1 (plain start) or == 0 with has_boundary_mr (boundary MR).
    if is_left and seg.left_mr <= 0 and (seg.left_mr < 0 or zm.has_boundary_mr):
        return ('circle', j)
    # Circle boundary θ₁ = 2π (last segment's right end) → fold onto left track.
    if is_right and seg.right_mr <= 0 and (seg.right_mr < 0 or zm.has_boundary_mr):
        return ('circle', int(boundary_perm_inv[j]))

    # Shared interior MR endpoint: same roots array on both sides, same track.
    if is_left and seg.left_mr > 0:
        return ('mr', seg.left_mr, j)
    if is_right and seg.right_mr > 0:
        return ('mr', seg.right_mr, j)

    # Interior sample — a unique exact touch, keep it.
    return ('interior', id(seg), i, j)


def _finalize_crossing(
    kind: str,
    i: int,
    j: int,
    th: np.ndarray,
    tr: np.ndarray,
    la: np.ndarray,
    mu2: float,
    mode: ExtractMode,
    poly: CharPoly,
    E: complex,
    mu1: float,
) -> tuple[float, complex]:
    """Turn a recorded crossing hit into (theta1, beta2) per ``mode``.

    - ``'coarse'`` : linear interpolation only.
    - ``'fine'``/``'solve'`` : ``_find_exact_crossing`` refinement, falling
      back to the linear estimate on solver failure.
    """
    if kind == 'zero':
        t1_lin = float(th[i])
        b2_lin = complex(tr[i, j])
    else:
        denom = la[i + 1, j] - la[i, j]
        frac_i = (mu2 - la[i, j]) / denom
        t1_lin = float(th[i] + frac_i * (th[i + 1] - th[i]))
        b2_lin = complex(tr[i, j] + frac_i * (tr[i + 1, j] - tr[i, j]))

    if mode == 'coarse':
        return t1_lin, b2_lin

    # 'fine' / 'solve' — refine (t1, angle(b2)) via fsolve + analytic Jacobian.
    res = _find_exact_crossing(
        poly, E, mu1, mu2, t1_lin % (2 * pi), float(np.angle(b2_lin)),
    )
    if res is None:
        return t1_lin, b2_lin
    t1, t2 = res
    return float(t1), exp(mu2 + 1j * t2)


# ---------------------------------------------------------------------------
# Winding (for the μ₂-bisection) — RAW crossing detection, no boundary rules
# ---------------------------------------------------------------------------

def amoeba_windings(
    zm: AmoebaZeroManager,
    poly: CharPoly,
    E: complex,
    mu1: float,
    mu2: float,
    *,
    tol: float = CONTINUUM_TOL,
    frac: float = CONTINUUM_FRAC,
    refine: bool = True,
) -> tuple[Optional[float], Optional[list], bool, Optional[float]]:
    """w2 average winding from ZM tracks at ``(E, mu1, mu2)``.

    Replaces ``_compute_winding_from_tracks``.  Same return contract:

    Returns
    -------
    winding : float or None
        w2 average winding number; ``None`` when a continuum is present.
    zeros : list of (theta1, theta2, jump) or None
        Crossing points; ``None`` when a continuum is present.
    has_continuum : bool
        Whether a track is in the mu2 band over a whole segment.
    dW_dmu2 : float or None
        Analytical derivative d(winding)/d(mu2); ``0.0`` when ``refine=False``
        or no crossings, ``None`` when a continuum is present.

    Uses RAW crossing detection — every crossing participates in the winding
    integral.  No Rule 1/Rule 2 boundary filtering (those are a subset-output
    concern, not a winding concern).
    """
    continuum_masks = _continuum_mask(zm, mu2, tol, frac)
    if any(m.any() for m in continuum_masks):
        return None, None, True, None

    zeros: list[tuple[float, float, int]] = []
    for s, seg in enumerate(zm.segments):
        th = seg.theta1_arr
        tr = seg.tracked_roots
        la = zm.seg_logabs[s]

        n_roots = tr.shape[1]
        d = la - mu2                                  # (N, K)

        # d == 0 exact touches: point is exactly on mu2.
        i_idx, j_idx = np.where(d == 0)
        for i, j in zip(i_idx, j_idx):
            t1, t2 = _crossing_thetas('zero', int(i), int(j), th, tr, la, mu2,
                                       refine, poly, E, mu1)
            jump = 1 if la[i, j] < mu2 else -1
            zeros.append((t1 % (2 * pi), t2 % (2 * pi), jump))

        sc = d[:-1] * d[1:] < 0                       # (N-1, K)
        i_idx, j_idx = np.where(sc)
        for i, j in zip(i_idx, j_idx):
            t1, t2 = _crossing_thetas('cross', int(i), int(j), th, tr, la, mu2,
                                       refine, poly, E, mu1)
            jump = 1 if la[i, j] < mu2 else -1
            zeros.append((t1 % (2 * pi), t2 % (2 * pi), jump))

    if not zeros:
        winding, _ = _get_average_winding_from_zeros(
            poly, E, mu1, mu2, [], direction=2,
        )
        return winding, [], False, 0.0

    # dW/dmu2 (only when refined): Σ jump · θ₁_dot / (2π).
    dW_dmu2 = 0.0
    if refine:
        from .ronkin_winding import _compute_zero_dtheta1_dmu2
        for t1, t2, jump in zeros:
            theta1_dot = _compute_zero_dtheta1_dmu2(
                poly, E, mu1, mu2, t1, t2,
            )
            dW_dmu2 += jump * theta1_dot
        dW_dmu2 /= (2 * pi)

    winding, _ = _get_average_winding_from_zeros(
        poly, E, mu1, mu2, zeros, direction=2,
    )
    return winding, zeros, False, dW_dmu2


def _crossing_thetas(
    kind: str,
    i: int,
    j: int,
    th: np.ndarray,
    tr: np.ndarray,
    la: np.ndarray,
    mu2: float,
    refine: bool,
    poly: CharPoly,
    E: complex,
    mu1: float,
) -> tuple[float, float]:
    """(theta1, theta2) for a winding crossing, linear or refined.

    Returns unwrapped θ₁ (caller wraps to [0, 2π)); θ₂ from the (possibly
    refined) β₂ angle.
    """
    if kind == 'zero':
        t1_lin = float(th[i])
        b2_lin = complex(tr[i, j])
    else:
        frac_i = (mu2 - la[i, j]) / (la[i + 1, j] - la[i, j])
        t1_lin = float(th[i] + frac_i * (th[i + 1] - th[i]))
        b2_lin = complex(tr[i, j] + frac_i * (tr[i + 1, j] - tr[i, j]))

    if not refine:
        return t1_lin, float(np.angle(b2_lin))

    res = _find_exact_crossing(
        poly, E, mu1, mu2, t1_lin % (2 * pi), float(np.angle(b2_lin)),
    )
    if res is None:
        return t1_lin, float(np.angle(b2_lin))
    return float(res[0]), float(res[1])
