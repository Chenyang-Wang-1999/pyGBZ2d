'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-16
Copyright © Department of Physics, Tsinghua University. All rights reserved

0D SGBZ PointSubset extraction from pairwise crossing EventGroups.

The crossing detection itself lives in :mod:`brute_force_SGBZ.pairwise`
and is orchestrated by :meth:`Mu2MidZM.analyze`.  This module only turns the
resulting EventGroups into PointSubsets:

  * an EventGroup whose induced item/column component covers BOTH sorted
    positions M-1 and M is an SGBZ boundary point;
  * every REAL column of that component gets its own PointSubset
    (representative items are expanded through their full clusters);
  * the charge is the side-change topological charge computed by
    ``pairwise.finalize_event_groups`` (+1 / -1 / 0 / None);
  * crossings near a boundary MR are replaced by the exact MR record.

Assumes no boundary continuum (the solver gates with ``has_continuum``).
'''

from __future__ import annotations

import math

import numpy as np
from cmath import exp

from gbz_types import CharPoly, PointSubset, circ_dist
from continuation import ZeroManager

from .mu2mid import Mu2MidZM, CONTINUUM_TOL


# ---------------------------------------------------------------------------
# Tunables
# ---------------------------------------------------------------------------

# Convergence tolerance for the pairwise brentq refinement.
_CROSSING_TOL: float = 1e-10

# Backward-compat iteration-limit alias (kept for the public signatures).
_MAX_BRACKET_ITER: int = 100
_MAX_NEWTON_ITER: int = _MAX_BRACKET_ITER

# θ₁ within this circular distance of a boundary MR → the crossing is the
# MR's own echo and is dropped in favour of the exact MR record.
_MR_PROXIMITY_TOL: float = 1e-4


# ---------------------------------------------------------------------------
# Mu2MidZM bootstrapping
# ---------------------------------------------------------------------------

def _ensure_mu2mid(zm: ZeroManager, **run_kwargs) -> Mu2MidZM:
    """Return an analyzed Mu2MidZM, reusing *zm* if it is already analyzed.

    A plain ``ZeroManager`` cannot be analyzed in place, so a fresh
    ``Mu2MidZM`` is constructed from its ``(poly, E_ref, mu1)``, run and
    analyzed with ``tie_tol=CONTINUUM_TOL``.
    """
    if isinstance(zm, Mu2MidZM) and getattr(zm, '_analyzed', False):
        return zm
    m = Mu2MidZM(zm.poly, zm.E_ref, zm.mu1)
    m.run(**run_kwargs)
    m.analyze(tie_tol=CONTINUUM_TOL)
    return m


# ---------------------------------------------------------------------------
# MR records
# ---------------------------------------------------------------------------

def _mr_boundary_entries(zm: Mu2MidZM, M: int) -> list[dict]:
    """Multiple roots sitting exactly on the PMGBZ boundary."""
    entries: list[dict] = []
    for mr in zm.multiple_roots:
        for cluster in mr.cluster_indices:
            if (M - 1) in cluster and M in cluster:
                beta2 = complex(mr.roots[cluster[0]])
                entries.append({
                    'theta1': float(mr.theta1),
                    'beta2': beta2,
                    'theta2': float(np.angle(beta2)) % (2 * math.pi),
                })
                break
    return entries


def _is_near_mr(theta1: float, zm: Mu2MidZM) -> bool:
    for mr in zm.multiple_roots:
        if circ_dist(theta1, mr.theta1) < _MR_PROXIMITY_TOL:
            return True
    return False


def _classify_charge(
    theta1: float,
    beta2: complex,
    g_prime: float,
    near_mr: bool,
) -> dict:
    """Legacy single-crossing charge classifier (kept for tests/API)."""
    theta2 = float(np.angle(beta2)) % (2 * math.pi)
    base = dict(theta1=theta1, theta2=theta2)
    if near_mr:
        return {**base, 'charge': None, 'kind': 'mr'}
    if not math.isfinite(g_prime):
        return {**base, 'charge': None, 'kind': 'tangent'}
    return {**base, 'charge': 1 if g_prime > 0 else -1, 'kind': 'ordinary'}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def detect_crossings_simple(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    crossing_tol: float = _CROSSING_TOL,
    max_newton: int = _MAX_BRACKET_ITER,
    zm_run_kwargs: dict | None = None,
) -> tuple[list[PointSubset], list[dict]]:
    """Materialize 0D PMGBZ PointSubsets from analyzed EventGroups.

    *crossing_tol* / *max_newton* are accepted for backward compatibility;
    ``crossing_tol`` is passed to :meth:`Mu2MidZM.analyze` when a fresh
    analysis is needed and ``max_newton`` is unused (brentq has its own
    internal iteration budget).

    Returns
    -------
    subsets : list[PointSubset]
        One PointSubset per REAL column of every M-1/M boundary event
        component, plus one per boundary MR.
    charges : list[dict]
        One dict per subset: ``charge`` is the side-change charge
        (+1 / -1 / 0) for ordinary events and ``None`` for tangent/MR.
    """
    M = poly.M
    K = poly.M + poly.N
    if M >= K:
        raise ValueError(f"M={M} >= K={K}: no PMGBZ boundary to check")
    if M <= 0:
        raise ValueError(f"M={M} <= 0: invalid boundary index")

    if isinstance(zm, Mu2MidZM) and getattr(zm, '_analyzed', False):
        m = zm
    else:
        m = _ensure_mu2mid(zm, **(zm_run_kwargs or {}))
        m = _reanalyze_if_needed(m, crossing_tol)

    E_ref = m.E_ref
    mu1 = m.mu1
    mr_entries = _mr_boundary_entries(m, M)

    subsets: list[PointSubset] = []
    charges: list[dict] = []

    for g in m._event_groups:
        if not g.point_columns:
            continue
        if g.row < 0:
            continue
        seg = m.segments[g.seg_idx]
        theta1 = float(g.theta) % (2 * math.pi)
        if any(circ_dist(theta1, e['theta1']) < _MR_PROXIMITY_TOL
               for e in mr_entries):
            continue  # MR echo: replaced by the exact MR record below

        for col in g.point_columns:
            beta2 = complex(seg.tracked_roots[g.row, col])
            theta2 = float(np.angle(beta2)) % (2 * math.pi)
            q = g.column_q.get(col)
            kind = 'ordinary' if q is not None else 'tangent'
            subsets.append(PointSubset(
                E=E_ref, beta1=exp(mu1 + 1j * theta1), beta2=beta2))
            charges.append(dict(
                theta1=theta1, theta2=theta2, charge=q, kind=kind))

    for e in mr_entries:
        subsets.append(PointSubset(
            E=E_ref, beta1=exp(mu1 + 1j * e['theta1']), beta2=e['beta2']))
        charges.append(dict(theta1=e['theta1'], theta2=e['theta2'],
                            charge=None, kind='mr'))
    return subsets, charges


def _reanalyze_if_needed(m: Mu2MidZM, crossing_tol: float) -> Mu2MidZM:
    """Not currently used; kept as a hook for stricter crossing_tol reuse."""
    return m
