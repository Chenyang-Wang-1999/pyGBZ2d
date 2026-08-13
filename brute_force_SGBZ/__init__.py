'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2025-11-18
Copyright © Department of Physics, Tsinghua University. All rights reserved

SGBZ (strip GBZ) computation via continuation.ZeroManager + μ₂_mid.

The SGBZ subset for a reference energy ``E_ref`` is found by bisecting
``mu1 = log|beta1|`` for the zero of the average major-axis winding number
``W(E_ref, mu1)``.  The piecewise-smooth ``μ₂_mid(θ₁)`` — the boundary-pair
modulus mean — is a first-class object (``Mu2MidZM``): it unifies crossing
detection (per-column vs μ₂_mid, §2) and the winding loop path (quad of
``Im[f'/f]`` over the piecewise-smooth path, §6.4), and inline-detects the
1D continuum (§1).  All faithful to ``log/2026-08-13-SGBZ算法梳理.md``.
'''

from gbz_types import (
    PointSubset, LineSubset, GBZResult, ConnectedSubset, CharPoly,
    get_minor_degrees,
)

from .mu2mid import (
    Mu2MidZM, ItemView, Mu2MidBreakpoint,
    CONTINUUM_TOL, CONTINUUM_FRAC,
)
from .continuum_lines import (
    detect_continuum_simple,
    extract_continuum_linesubsets,
)
from .crossings import detect_crossings_simple
from .winding import (
    detect_crossings_and_winding,
    compute_average_winding,
    WindingFun,
    get_winding_number,
)
from .sgbz_solver import (
    solve_SGBZ_for_E,
    collect_GBZ_subsets,
)

__all__ = [
    "CharPoly",
    "get_minor_degrees",

    # μ₂_mid builder (§1/§2)
    "Mu2MidZM",
    "ItemView",
    "Mu2MidBreakpoint",
    "CONTINUUM_TOL",
    "CONTINUUM_FRAC",

    # continuum detection + materialization (§1)
    "detect_continuum_simple",
    "extract_continuum_linesubsets",

    # crossing detection + winding (§2, §3)
    "detect_crossings_simple",
    "detect_crossings_and_winding",
    "compute_average_winding",
    "WindingFun",
    "get_winding_number",

    # top-level solver
    "solve_SGBZ_for_E",
    "collect_GBZ_subsets",

    "PointSubset",
    "LineSubset",
    "GBZResult",
    "ConnectedSubset",
]
