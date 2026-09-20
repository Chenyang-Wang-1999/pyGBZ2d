'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-07-24
Copyright © Department of Physics, Tsinghua University. All rights reserved

ZeroManager-based amoeba helpers.

This module currently contains only the pieces needed by the μ₂ bisection
stage:

  - ``AmoebaZeroManager`` : ``ZeroManager`` specialization caching per-segment
    ``ln|β₂|`` so detection / crossing code can index a precomputed array.
  - ``detect_continuum`` : std-based continuum detection, run BEFORE the μ₂
    bisection.  Returns groups ``[(mu2_c, [(seg_idx, col_idx), ...]), ...]``
    with nearby ``mu2_c`` values merged.
  - ``find_crossings`` : raw crossing detection of ``ln|β₂| = μ₂``, with
    optional persistent mesh refinement and per-(segment, column) avoidance.
'''

from __future__ import annotations

from cmath import exp
from typing import Optional

import numpy as np

from ..core import TWO_PI, CharPoly
from ..continuation.zero_manager import ZeroManager

from .ronkin_winding import (
    _find_exact_crossing,
    _get_average_winding_from_zeros,
)


# ---------------------------------------------------------------------------
# AmoebaZeroManager — ZeroManager with cached per-segment ln|β₂|
# ---------------------------------------------------------------------------

class AmoebaZeroManager(ZeroManager):
    """``ZeroManager`` with per-segment ``ln|β₂|`` cached for amoeba use.

    After ``.run()``:
      - ``seg_logabs`` : ``list[np.ndarray]``, ``seg_logabs[s]`` has shape
        ``(N_s, K)`` = ``np.log(np.abs(segments[s].tracked_roots))``.

    ``insert_solution`` keeps the affected cache synchronized so that samples
    obtained while refining one μ₂ remain usable for subsequent μ₂ levels.
    """

    seg_logabs: list[np.ndarray]

    def __init__(self, poly: CharPoly, E_ref: complex, mu1: float):
        super().__init__(poly, E_ref, mu1)
        self.seg_logabs = []

    def run(self, **kwargs) -> None:
        super().run(**kwargs)
        self.refresh_logabs()

    def insert_solution(self, theta1, seg_idx=None, i=None, *, interp='hermite'):
        if seg_idx is None or i is None:
            seg_idx, i = self.locate(theta1)
        row, changed = super().insert_solution(
            theta1, seg_idx=seg_idx, i=i, interp=interp)
        if changed:
            if len(self.seg_logabs) != len(self.segments):
                self.refresh_logabs()
            else:
                self.seg_logabs[seg_idx] = np.log(
                    np.abs(self.segments[seg_idx].tracked_roots))
        return row, changed

    def refresh_logabs(self) -> None:
        self.seg_logabs = [
            np.log(np.abs(seg.tracked_roots)) for seg in self.segments
        ]


# ---------------------------------------------------------------------------
# Continuum detection (std-based, pre-bisection)
# ---------------------------------------------------------------------------

def detect_continuum(
    zm: AmoebaZeroManager,
    tol: float,
) -> list[tuple[float, list[tuple[int, int]]]]:
    """Detect per-track continuum segments by ``std(ln|β₂|) < tol``.

    Returns groups ``[(mu2_c, [(seg_idx, col_idx), ...]), ...]``:
      - a track ``(s, j)`` is continuum when the segment has at least 2 rows,
        all its ``ln|β₂|`` values are finite, and ``std(la[:, j]) < tol``;
      - ``mu2_c`` is the track's mean ``ln|β₂|``;
      - continuum tracks whose ``mu2_c`` values are closer than ``tol`` are
        merged into one group (numerically they are the same flat level).
    """
    entries: list[tuple[int, int, float]] = []
    for s, seg in enumerate(zm.segments):
        if len(seg.theta1_arr) < 2:
            raise ValueError(f"Segment {s} has less than 2 rows.")
        la = zm.seg_logabs[s]
        for j in range(zm.K):
            col = la[:, j]
            if not np.all(np.isfinite(col)):
                continue
            if float(np.std(col)) < tol:
                entries.append((int(s), int(j), float(np.mean(col))))

    if not entries:
        return []

    entries.sort(key=lambda x: x[2])

    groups: list[list] = []  # [mu2_c, [(s,j), ...], sum_mu, count]
    for s, j, mu in entries:
        if not groups or mu - groups[-1][0] >= tol:
            groups.append([mu, [(s, j)], mu, 1])
        else:
            g = groups[-1]
            g[1].append((s, j))
            g[2] += mu
            g[3] += 1
            g[0] = g[2] / g[3]

    return [(float(g[0]), g[1]) for g in groups]


# ---------------------------------------------------------------------------
# Crossing detection
# ---------------------------------------------------------------------------

def _crossing_intervals(zm, mu2, avoided):
    """Snapshot candidates; callers must restart this scan after insertion."""
    for s, seg in enumerate(zm.segments):
        th = seg.theta1_arr
        la = (zm.seg_logabs[s] if isinstance(zm, AmoebaZeroManager)
              else np.log(np.abs(seg.tracked_roots)))
        d = la - mu2
        for j in range(zm.K):
            if (s, j) in avoided:
                continue
            for i in np.flatnonzero(d[:-1, j] == 0):
                t = float(th[i])
                yield s, j, t, t, None
            left, right = d[:-1, j], d[1:, j]
            cross = ((left < 0) & (right > 0)) | ((left > 0) & (right < 0))
            for i in np.flatnonzero(cross):
                yield s, j, float(th[i]), float(th[i + 1]), bool(left[i] < 0)


def _find_refined_crossings(zm, mu2, avoided):
    # An accepted root is a real mesh row. Keep its crossing orientation as
    # well as its coordinate: two nearby opposite crossings may share a row
    # as a bracket endpoint, but must never stand in for one another.
    accepted = {}
    while True:
        crossings = []
        for s, j, lo, hi, negative_left in _crossing_intervals(zm, mu2, avoided):
            result = accepted.get((s, j, lo, negative_left))
            if result is None:
                result = accepted.get((s, j, hi, negative_left))
            if result is None:
                t1, t2 = _find_exact_crossing(zm, mu2, s, j, lo, hi)
                result = (exp(zm.mu1 + 1j * t1), exp(mu2 + 1j * t2))
                accepted[s, j, t1, negative_left] = result
                # Refinement inserts entire root rows. Rescan every column,
                # including previously visited ones, for newly exposed events.
                break
            crossings.append(result)
        else:
            return crossings


def find_crossings(
    zm: ZeroManager,
    mu1: float,
    mu2: float,
    avoided_segments: Optional[list[tuple[int, int]]] = None,
    return_refined: bool = False,
) -> list[tuple[complex, complex]]:
    """Find all crossings of ``ln|β₂| = μ₂`` over the ZM tracks.

    Traverses every segment and every track column.  A ``(seg_idx, col_idx)``
    present in ``avoided_segments`` is skipped entirely.  Exact touches at
    mesh rows are detected first via ``d[:-1] == 0`` (only the left endpoint,
    so an interval endpoint shared by two intervals is not double-counted),
    then sign changes are detected via ``d[1:] * d[:-1] < 0``.

    Returns a list of ``(beta1, beta2)`` crossing points.  When
    ``return_refined`` is False the coordinates are linearly interpolated;
    otherwise bracketed polynomial solves refine each crossing and persist
    their samples in ``zm``. The updated mesh is rescanned after refinement.
    Failure raises rather than supplying an unrefined winding partition.
    """
    avoided = set(avoided_segments or ())
    if return_refined:
        if mu1 != zm.mu1:
            raise ValueError("crossing refinement requires the ZM's fixed mu1")
        return _find_refined_crossings(zm, mu2, avoided)

    crossings: list[tuple[complex, complex]] = []
    for s, seg in enumerate(zm.segments):
        th = seg.theta1_arr
        tr = seg.tracked_roots
        la = (zm.seg_logabs[s] if isinstance(zm, AmoebaZeroManager)
              else np.log(np.abs(tr)))
        d = la - mu2                                  # (N, K)

        for j in range(zm.K):
            if (s, j) in avoided:
                continue

            # Exact touches at mesh rows (left endpoint only — a zero at
            # d[i+1] would be the right endpoint of interval i, and is also
            # the left endpoint of interval i+1, so this avoids duplicates).
            zero_rows = np.where(d[:-1, j] == 0)[0]
            for i in zero_rows:
                t1 = float(th[i])
                b2 = complex(tr[i, j])
                t2 = float(np.angle(b2))
                crossings.append((exp(mu1 + 1j * t1), exp(mu2 + 1j * t2)))

            # Sign changes inside intervals.
            sign_change = d[:-1, j] * d[1:, j] < 0
            for i in np.where(sign_change)[0]:
                denom = la[i + 1, j] - la[i, j]
                frac = (mu2 - la[i, j]) / denom
                t1 = float(th[i] + frac * (th[i + 1] - th[i]))
                b2 = complex(tr[i, j] + frac * (tr[i + 1, j] - tr[i, j]))
                t2 = float(np.angle(b2))
                crossings.append((exp(mu1 + 1j * t1), exp(mu2 + 1j * t2)))

    return crossings


def _calculate_a2_winding_and_zeros(
    zm: ZeroManager,
    mu1: float,
    mu2: float,
    avoided_segments: Optional[list[tuple[int, int]]] = None,
    return_refined: bool = False,
) -> tuple[float, list[tuple[float, float]]]:
    """Keep the exact partition used by the unchanged midpoint winding solve."""
    crossings = find_crossings(
        zm, mu1, mu2,
        avoided_segments=avoided_segments,
        return_refined=return_refined,
    )
    zeros = [
        (float(np.angle(b1) % TWO_PI), float(np.angle(b2) % TWO_PI))
        for b1, b2 in crossings
    ]
    winding, _ = _get_average_winding_from_zeros(
        zm.poly, zm.E_ref, mu1, mu2, zeros, direction=2,
    )
    return winding, zeros


def calculate_a2_average_winding(
    zm: ZeroManager,
    mu1: float,
    mu2: float,
    avoided_segments: Optional[list[tuple[int, int]]] = None,
    return_refined: bool = False,
) -> float:
    """a2 average winding, with optional persistent crossing refinement."""
    winding, _ = _calculate_a2_winding_and_zeros(
        zm, mu1, mu2, avoided_segments=avoided_segments,
        return_refined=return_refined)
    return winding
