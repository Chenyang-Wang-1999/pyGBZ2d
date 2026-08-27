'''
bfGBZ2d — Brute-force GBZ computation for 2D non-Hermitian models.

Top-level re-exports the unified data types and the CharPoly entry point;
the two GBZ formulations live in their own subpackages:

  - ``bfgbz2d.sgbz``       strip-GBZ / average major-axis winding
  - ``bfgbz2d.amoeba``     amoeba / Ronkin-function formulation
  - ``bfgbz2d.continuation`` pseudo-arclength β₂-root tracking along θ₁

Both ``bfgbz2d.sgbz.collect_GBZ_subsets`` and
``bfgbz2d.amoeba.collect_GBZ_subsets`` return the unified ``GBZResult``.
'''

__version__ = "0.1.0"

from bfgbz2d.core import (
    CharPoly,
    PointSubset,
    LineSubset,
    ConnectedSubset,
    GBZResult,
    TWO_PI,
)

from bfgbz2d import continuation, sgbz, amoeba

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
