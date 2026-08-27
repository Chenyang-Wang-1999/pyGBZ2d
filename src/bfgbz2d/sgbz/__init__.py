'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2025-11-18
Copyright © Department of Physics, Tsinghua University. All rights reserved

SGBZ (strip GBZ) computation via continuation.ZeroManager + μ₂_mid.

The SGBZ subset for a reference energy ``E_ref`` is found by bisecting
``mu1 = log|beta1|`` for the zero of the average major-axis winding number
``W(E_ref, mu1)``.  Crossing detection is pairwise on ItemView
representatives (``Mu2MidZM.analyze`` → EventGroups, see ``pairwise.py``);
``winding.detect_crossings_simple`` materializes those groups into
PointSubsets and the loop-winding path ``μ₂_mid(θ₁)`` is the independent
piecewise-smooth ``Mu2Mid`` object.
'''

from bfgbz2d.core import (
    PointSubset, LineSubset, GBZResult, ConnectedSubset, CharPoly,
    get_minor_degrees,
)

from bfgbz2d.config import CONTINUUM_TOL, CONTINUUM_FRAC
from .mu2mid import (
    Mu2MidZM, Mu2Mid, Mu2MidPiece,
    ItemView,
)
from .continuum_lines import (
    detect_continuum_simple,
    extract_continuum_linesubsets,
)
from .winding import (
    detect_crossings_simple,
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
    "Mu2Mid",
    "Mu2MidPiece",
    "ItemView",
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
