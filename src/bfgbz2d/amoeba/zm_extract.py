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

from bfgbz2d.core import (
    TWO_PI,
    CharPoly, PointSubset, LineSubset,
    JoinableLinePiece, is_mr_cluster_endpoint,
)

from .ronkin_winding import _find_exact_crossing, _get_average_winding_from_zeros
from bfgbz2d.continuation.zero_manager import SegmentData, ZeroManager


# A genuine continuum is constant-modulus to ~1e-9; a transversal crossing
# leaves the level after one sample.  1e-6 separates them.
CONTINUUM_TOL = 1e-6
CONTINUUM_FRAC = 0.9
# A crossing whose θ₁ lands within this of a continuum LineSubset endpoint is
# treated as the continuum/MR boundary, not a genuine discrete zero, and
# dropped.  Decoupled from CONTINUUM_TOL so the snap radius is independent of
# band-detection sensitivity.
SNAP_TOL = 1e-3
# Root-match tolerance for Rule 1's curve-consistency screen: the shared
# boundary rows carry the SAME root values (an MR's roots are stored once and
# shared by both segments), so the same curve matches to machine precision; a
# different root (same θ₁, unconnected track) differs by O(1).
ROOT_TOL = 1e-9

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
        # ZeroManager.run is single-shot (a second call raises), so this
        # cache is built exactly once over the final segment list.
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
# MR-boundary joining for continuum LineSubsets
# ---------------------------------------------------------------------------
#
# _LinePiece / _is_cluster_endpoint live in bfgbz2d.core (as JoinableLinePiece /
# is_mr_cluster_endpoint) — they are shared with bfgbz2d.sgbz's continuum
# extractor.  The old underscore names remain as aliases for the in-module
# call sites and the tests that construct pieces directly.

_LinePiece = JoinableLinePiece
_is_cluster_endpoint = is_mr_cluster_endpoint


def _join_continuum_across_mrs(
    zm: AmoebaZeroManager,
    continuum_masks: list[np.ndarray],
    line_pieces: list[_LinePiece],
) -> list[_LinePiece]:
    """Join per-segment continuum LineSubsets that pass through MRs.

    A segment boundary is an MR, but only the roots in its ``cluster_indices``
    are genuinely multiple; the rest are regular roots that pass straight
    through.  A continuum track whose endpoint root is NOT in the cluster
    therefore does not terminate at the MR — it continues into the adjacent
    segment on the matched track.  This joins those two pieces, and raises if
    the expected continuation is absent (a topology inconsistency).

    Each original segment's left boundary is matched against the previous
    segment's right boundary by COLUMN IDENTITY.  Interior MR rows share one
    track frame; the first/last-segment seam is translated through
    ``boundary_perm``.  Iterated to a fixpoint so chains and the cyclic seam
    (segment 0 ↔ last segment) both converge.
    """
    n_seg = len(zm.segments)
    if n_seg <= 1:
        return line_pieces

    def _mod(k: int) -> int:
        return k % n_seg

    def find_by_left(seg_s: int, root: complex) -> int | None:
        # Piece whose leftmost spanned segment is seg_s and whose leftmost
        # β₂ ≈ root (continuum tracks on one segment are distinct roots).
        # ROOT_TOL: the same root-value scale used by the Rule 1 screen —
        # continuum tracks are separated by O(1), so 1e-9 cleanly separates
        # "same track" from "adjacent track".
        for idx, p in enumerate(line_pieces):
            if p.ml == seg_s and np.abs(p.beta2_arr[0] - root) < ROOT_TOL:
                return idx
        return None

    def find_by_right(seg_s: int, root: complex) -> int | None:
        for idx, p in enumerate(line_pieces):
            if p.mr == seg_s and np.abs(p.beta2_arr[-1] - root) < ROOT_TOL:
                return idx
        return None

    def _match_column_right(cur_s: int, prev_s: int, col: int) -> int:
        # The only frame change is the cyclic first/last-segment seam.  Its
        # documented monodromy is roots_right[boundary_perm] == roots_left.
        if cur_s == 0 and prev_s == n_seg - 1:
            return int(zm.boundary_perm[int(col)])
        return int(col)

    for _ in range(n_seg):
        changed = False
        for s in range(n_seg):
            seg = zm.segments[s]
            prev = _mod(s - 1 + n_seg)
            prev_seg = zm.segments[prev]

            left_mr = seg.left_mr
            right_mr = prev_seg.right_mr
            # A shared interior MR boundary has the same index on both sides.
            # The circle seam (θ₁=0 ≡ 2π) is the only non-MR shared boundary:
            # it occurs iff segment 0's left_mr and the last segment's
            # right_mr are both -1 (no boundary MR).
            is_circle_seam = (left_mr < 0 and right_mr < 0)
            if not is_circle_seam and left_mr != right_mr:
                continue  # not a shared boundary

            left_b = seg.tracked_roots[0, :]
            right_b = prev_seg.tracked_roots[-1, :]

            for j_l in np.where(continuum_masks[s])[0]:
                root = left_b[j_l]
                if _is_cluster_endpoint(zm, seg, 'left', int(j_l)):
                    continue  # genuine LineSubset terminator at the MR
                li_idx = find_by_left(s, root)
                if li_idx is None:
                    # segment s is not the leftmost of any piece → its left
                    # boundary is already interior to a merged piece (the
                    # continuum was joined here in an earlier pass).
                    continue
                # Match the continuation column on the previous segment's
                # right boundary by track identity (or boundary_perm at the
                # cyclic seam), then read its exact root value.
                j_prev = _match_column_right(s, prev, int(j_l))
                root_prev = right_b[j_prev]
                if _is_cluster_endpoint(zm, prev_seg, 'right', j_prev):
                    # Non-cluster on one side, cluster on the other — the track
                    # ends here in one segment but not the other: inconsistent.
                    raise ValueError(
                        f"Continuum track {j_l} of segment {s} ends at the MR "
                        f"at θ₁={seg.theta1_arr[0]:.4f} as a non-cluster root, "
                        f"but the matched root (track {j_prev}) of segment "
                        f"{prev} is a cluster root there."
                    )
                pi_idx = find_by_right(prev, root_prev)
                if pi_idx is None:
                    raise ValueError(
                        f"Continuum track {j_l} of segment {s} continues "
                        f"through the MR at θ₁={seg.theta1_arr[0]:.4f} as a "
                        f"non-cluster root, but segment {prev} has no "
                        f"continuum on the matched track (track {j_prev})."
                    )
                if pi_idx == li_idx:
                    continue  # already joined (e.g. the cyclic seam)
                _merge_two(line_pieces, pi_idx, li_idx)
                changed = True
                break  # line_pieces changed; restart the segment scan
            if changed:
                break
        if not changed:
            break
    return line_pieces


def _merge_two(
    line_pieces: list[_LinePiece],
    prev_idx: int, cur_idx: int,
) -> None:
    """Merge ``line_pieces[prev_idx]`` (right side) with ``[cur_idx]`` (left).

    The array is ``[prev, cur[1:]]`` for BOTH the interior-MR and the cyclic
    seam case.  Interior MR: ``prev`` (segment s-1) is to the LEFT of ``cur``
    (segment s) in θ₁, so θ₁ stays monotonic and the shared MR row (cur's
    first row) is dropped.  Cyclic seam (θ₁=0 ≡ 2π): ``prev`` is the LAST
    segment (right end at 2π), ``cur`` is segment 0 (left end at 0); putting
    rp first aligns the seam (rp[-1] and cp[0] are the same physical point)
    INSIDE the array — where the track is continuous through θ₁=0/2π — so
    the array's two ends fall on the real terminators (interior cluster
    MRs).  Reversing this — ``[cp, rp[1:]]`` — would splice at the wrong
    physical point and break β₂ continuity.

    The merged endpoints are rp.ml / cp.mr in both cases.  (An earlier
    version swapped them at the seam, which mislabeled the piece and let
    the join loop re-match an already-joined piece — visible only when MR
    snapping makes a far endpoint's root coincide with the seam root.)
    """
    rp = line_pieces[prev_idx]
    cp = line_pieces[cur_idx]
    th = np.concatenate([rp.theta1_arr, cp.theta1_arr[1:]])
    b2 = np.concatenate([rp.beta2_arr, cp.beta2_arr[1:]])
    merged = _LinePiece(
        E=cp.E, mu1=cp.mu1, theta1_arr=th, beta2_arr=b2,
        ml=rp.ml, mr=cp.mr,
    )
    keep = [i for i in range(len(line_pieces)) if i not in (prev_idx, cur_idx)]
    new_list = [line_pieces[i] for i in keep]
    new_list.append(merged)
    line_pieces.clear()
    line_pieces.extend(new_list)


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
        endpoint (continuum/MR edge, not a genuine discrete zero), provided its
        β₂ matches an in-band endpoint root within ``ROOT_TOL`` — the crossing
        is then the line's own endpoint curve.  A same-θ₁ crossing on a
        different, unconnected root is a genuine discrete zero and survives.
      - Rule 2 — ``d == 0`` exact touches are deduplicated by zero-point
        identity, not by θ₁.  A zero's identity is ``(endpoint, track)``:
        an interior touch gets a unique key; a touch on a shared MR endpoint
        collides with the same MR + track on the neighbouring segment; a
        touch on the θ₁=0≡2π circle boundary collides with its counterpart
        on the other side via ``zm.boundary_perm``.  Distinct tracks at the
        same θ₁ (e.g. a conjugate pair at β₁=1) thus stay distinct.
    """
    continuum_masks = _continuum_mask(zm, mu2, tol, frac)
    # One per-segment-per-continuum-track LineSubset, *before* joining across
    # MR boundaries.  ``ml``/``mr`` track the outermost segment indices so that
    # merges can be chained and the join at the cyclic seam detected.
    line_pieces: list[_LinePiece] = []
    # Crossings recorded as (seg_idx, i, j, kind), NOT PointSubsets, so that
    # boundary duplicates between segments sharing an MR endpoint can be
    # resolved after the loop.  kind: 'zero' (d==0 at sample i) | 'cross'
    # (sign change over [i, i+1]).
    hits: list[tuple[int, int, int, str]] = []

    for s, seg in enumerate(zm.segments):
        th = seg.theta1_arr
        tr = seg.tracked_roots
        is_cont = continuum_masks[s]

        # 1. continuum tracks → one LineSubset per segment (joined across MR
        # boundaries below — a track that passes through an MR as a non-cluster
        # root does *not* terminate there).
        for j in np.where(is_cont)[0]:
            line_pieces.append(_LinePiece(
                E=E, mu1=mu1,
                theta1_arr=th.copy(), beta2_arr=tr[:, j].copy(),
                ml=s, mr=s))

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

    # ---- join per-segment LineSubsets across MR boundaries ----
    # A segment ends at an MR only because *some* roots there form a degenerate
    # cluster (``cluster_indices``); the other roots are regular and pass
    # straight through.  A continuum track whose endpoint root is NOT in the
    # cluster therefore does not terminate at the MR — it continues into the
    # adjacent segment on the matched track.  Join such pieces, and raise if the
    # expected continuation is absent (a topology inconsistency).
    line_pieces = _join_continuum_across_mrs(zm, continuum_masks, line_pieces)

    # ---- materialize + boundary filtering ----
    cont_endpoints = np.array(
        [t for p in line_pieces for t in (p.theta1_start, p.theta1_end)]
    ) if line_pieces else None
    # Endpoint roots, element-wise aligned with cont_endpoints (each piece's
    # start / end pair) — feeds Rule 1's curve-consistency screen.
    cont_endpoint_roots = np.array(
        [r for p in line_pieces for r in (p.beta2_arr[0], p.beta2_arr[-1])]
    ) if line_pieces else None

    # Inverse of zm.boundary_perm: boundary_perm[inv[k]] == k, so
    # roots_right[j] corresponds to left-boundary track inv[j].  Used to fold
    # the θ₁=2π end of the last segment onto the θ₁=0 start of the first
    # segment — they are the same physical circle point.
    boundary_perm_inv = np.empty(zm.K, dtype=int)
    boundary_perm_inv[zm.boundary_perm] = np.arange(zm.K)

    seen_zero: set[tuple] = set()
    points: list[tuple[float, complex]] = []
    for s, i, j, kind in hits:
        seg = zm.segments[s]
        th, tr = seg.theta1_arr, seg.tracked_roots
        la = zm.seg_logabs[s]
        t1, b2 = _finalize_crossing(kind, i, j, th, tr, la, mu2, mode, poly, E, mu1)

        # Rule 1: continuum-endpoint snap (θ₁), plus a root-distance screen —
        # only a crossing whose β₂ matches an in-band endpoint root (the SAME
        # curve as the LineSubset) is dropped.  Each hit is judged on its own:
        # b2 is THIS crossing's root; the screen asks whether ANY in-band
        # endpoint root matches it (∃ — a single scalar vs the endpoint-root
        # array).  A boundary point on a different, unconnected track (same
        # θ₁, root differs by O(1)) survives.  θ₁ distance is CIRCULAR:
        # refined crossings wrap to [0, 2π) while a closed piece's right
        # endpoint sits at exactly 2π — a linear distance would read ~2π
        # across the seam and the snap would miss (duplicate PointSubset).
        if cont_endpoints is not None:
            d_circ = np.abs(
                (t1 - cont_endpoints + pi) % (TWO_PI) - pi
            )
            in_band = d_circ < snap_tol
            if np.any(in_band) and np.any(
                np.abs(cont_endpoint_roots[in_band] - b2) < ROOT_TOL
            ):
                continue

        # Rule 2: dedup exact-touch ('zero') hits by zero identity.
        # 'cross' hits are interior sign-changes and are always unique.
        if kind == 'zero':
            key = _zero_identity_key(zm, seg, i, j, boundary_perm_inv)
            if key in seen_zero:
                continue
            seen_zero.add(key)

        points.append((t1, b2))

    # Drop the joining-only ``ml``/``mr`` fields: the public GBZ output is
    # plain ``LineSubset``s, not ``_LinePiece``s.
    subsets: list = [
        LineSubset(E=p.E, mu1=p.mu1,
                   theta1_arr=p.theta1_arr, beta2_arr=p.beta2_arr)
        for p in line_pieces
    ]
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

    # Circle boundary θ₁ = 0 (first segment's left end) ≡ θ₁ = 2π (last
    # segment's right end).  The seam is "no interior MR here": left_mr < 0
    # (plain start) OR the boundary MR at θ₁=0 (index 0, has_boundary_mr=True).
    # When has_boundary_mr is False, MR index 0 is a normal interior MR, NOT the
    # seam — same conflation as _is_cluster_endpoint; do not treat it as circle.
    left_is_circle = seg.left_mr < 0 or (seg.left_mr == 0 and zm.has_boundary_mr)
    right_is_circle = seg.right_mr < 0 or (seg.right_mr == 0 and zm.has_boundary_mr)
    if is_left and left_is_circle:
        return ('circle', j)
    if is_right and right_is_circle:
        return ('circle', int(boundary_perm_inv[j]))

    # Shared interior MR endpoint: same roots array on both sides, same track.
    # Reaches here for any left_mr/right_mr >= 0 that is not the circle seam
    # — including interior MR index 0 when has_boundary_mr is False.
    if is_left:
        return ('mr', seg.left_mr, j)
    if is_right:
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
        poly, E, mu1, mu2, t1_lin % (TWO_PI), float(np.angle(b2_lin)),
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
            zeros.append((t1 % (TWO_PI), t2 % (TWO_PI), jump))

        sc = d[:-1] * d[1:] < 0                       # (N-1, K)
        i_idx, j_idx = np.where(sc)
        for i, j in zip(i_idx, j_idx):
            t1, t2 = _crossing_thetas('cross', int(i), int(j), th, tr, la, mu2,
                                       refine, poly, E, mu1)
            jump = 1 if la[i, j] < mu2 else -1
            zeros.append((t1 % (TWO_PI), t2 % (TWO_PI), jump))

    if not zeros:
        winding, _ = _get_average_winding_from_zeros(
            poly, E, mu1, mu2, [], direction=2,
        )
        return winding, [], False, 0.0

    # dW/dmu2 (only when refined): Σ jump · θ₁_dot / (2π).  The sum runs
    # over UNIQUE zeros: adjacent segments share their boundary MR row, so a
    # d==0 exact touch on that row is recorded once per segment — the same
    # physical zero.  (The winding integral above already folds duplicates
    # via np.unique inside _get_average_winding_from_zeros; the derivative
    # sum must match or the Newton step doubles the shared-touch
    # contribution.)
    dW_dmu2 = 0.0
    if refine:
        from .ronkin_winding import _compute_zero_dtheta1_dmu2
        seen_zero_keys: set[tuple[float, float]] = set()
        for t1, t2, jump in zeros:
            key = (t1 % (TWO_PI), t2 % (TWO_PI))
            if key in seen_zero_keys:
                continue
            seen_zero_keys.add(key)
            theta1_dot = _compute_zero_dtheta1_dmu2(
                poly, E, mu1, mu2, t1, t2,
            )
            dW_dmu2 += jump * theta1_dot
        dW_dmu2 /= (TWO_PI)

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
        poly, E, mu1, mu2, t1_lin % (TWO_PI), float(np.angle(b2_lin)),
    )
    if res is None:
        return t1_lin, float(np.angle(b2_lin))
    return float(res[0]), float(res[1])
