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
    optional fsolve refinement and per-(segment, column) avoidance.
'''

from __future__ import annotations

from cmath import exp
from typing import Optional

import numpy as np

from pygbz2d.core import TWO_PI, CharPoly
from pygbz2d.continuation.zero_manager import ZeroManager

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

    The cache is rebuilt by :meth:`refresh_logabs` after mesh insertion
    (``insert_solution`` mutates the segment rows without touching this cache).
    """

    seg_logabs: list[np.ndarray]

    def __init__(self, poly: CharPoly, E_ref: complex, mu1: float):
        super().__init__(poly, E_ref, mu1)
        self.seg_logabs = []

    def run(self, **kwargs) -> None:
        super().run(**kwargs)
        self.refresh_logabs()

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
    otherwise :func:`_find_exact_crossing` refines each one (falling back to
    the linear estimate when the fsolve fails).
    """
    avoided = set(avoided_segments or ())

    crossings: list[tuple[complex, complex]] = []
    for s, seg in enumerate(zm.segments):
        th = seg.theta1_arr
        tr = seg.tracked_roots
        la = np.log(np.abs(tr))
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
                if return_refined:
                    res = _find_exact_crossing(
                        zm.poly, zm.E_ref, mu1, mu2, t1, t2,
                    )
                    if res is not None:
                        t1, t2 = res
                crossings.append((exp(mu1 + 1j * t1), exp(mu2 + 1j * t2)))

    return crossings


def calculate_a2_average_winding(
    zm: ZeroManager,
    mu1: float,
    mu2: float,
    avoided_segments: Optional[list[tuple[int, int]]] = None,
    return_refined: bool = False,
) -> float:
    """a2 average winding from crossings of ``ln|β₂| = μ₂``.

    Thin consumer of :func:`find_crossings`: it converts the returned
    ``(β₁, β₂)`` points back to angular zeros and feeds them to
    :func:`_get_average_winding_from_zeros` with ``direction=2``.
    """
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
    return winding
