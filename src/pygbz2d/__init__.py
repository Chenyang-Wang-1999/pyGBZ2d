'''
pyGBZ2d — Brute-force GBZ computation for 2D non-Hermitian models.

Top-level re-exports the unified data types and the CharPoly entry point;
the two GBZ formulations live in their own subpackages:

  - ``pygbz2d.sgbz``       strip-GBZ / average major-axis winding
  - ``pygbz2d.amoeba``     amoeba / Ronkin-function formulation
  - ``pygbz2d.continuation`` adaptive-step β₂-root tracking (pseudo-arclength
                            parameterization + error-controlled step size)

Both ``pygbz2d.sgbz.collect_GBZ_subsets`` and
``pygbz2d.amoeba.collect_GBZ_subsets`` return the unified ``GBZResult``.
'''

__version__ = "0.1.0"

from .core import (
    CharPoly,
    PointSubset,
    LineSubset,
    ConnectedSubset,
    GBZResult,
    TWO_PI,
)

from . import continuation, sgbz, amoeba

__all__ = [
    "CharPoly",
    "PointSubset",
    "LineSubset",
    "ConnectedSubset",
    "GBZResult",
    "TWO_PI",
    "continuation",
    "sgbz",
    "amoeba",
    "__version__",
]
