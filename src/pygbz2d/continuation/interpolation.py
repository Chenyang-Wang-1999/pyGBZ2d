"""
Cubic Hermite interpolation polynomial coefficients.

The continuation machinery turns the discrete β₂-root mesh into a continuous
curve, and every interpolation site in this package uses the same two-point
cubic Hermite recipe.  These two helpers centralise that recipe; all returned
arrays follow the numpy ``poly`` convention (highest power first), so callers
can evaluate with ``np.polyval``, differentiate with ``np.polyder`` and solve
curve intersections with ``np.roots(np.polysub(p1, p2))``.

The two functions are deliberately split by semantics:

* :func:`cubic_hermite_poly` is the pure cubic builder — no input validation.
* :func:`hermite_interp_poly` is the interpolation builder — it checks the
  endpoint derivatives and falls back to the linear polynomial when either
  derivative is divergent / undefined.
"""

from __future__ import annotations

import numpy as np


def cubic_hermite_poly(
    h: float,
    v0,
    dv0,
    v1,
    dv1,
) -> np.ndarray:
    """Pure cubic Hermite coefficients for ``f(0)=v0, f'(0)=dv0, f(h)=v1,
    f'(h)=dv1`` on ``x ∈ [0, h]``.

    Returns ``[a, b, c, d]`` in numpy ``poly`` order (highest power first):
    ``f(x) = a·x³ + b·x² + c·x + d``.

    NO input validation: ``h`` must be positive, and the values / derivatives
    must be finite.  Use :func:`hermite_interp_poly` when endpoint derivatives
    may be divergent.
    """
    h2 = h * h
    h3 = h2 * h
    a = (2.0 * (v0 - v1)) / h3 + (dv0 + dv1) / h2
    b = (3.0 * (v1 - v0)) / h2 - (2.0 * dv0 + dv1) / h
    return np.array([a, b, dv0, v0])


def hermite_interp_poly(
    h: float,
    v0,
    dv0,
    v1,
    dv1,
) -> np.ndarray:
    """Interpolation polynomial coefficients for the same endpoint data.

    Checks ONLY the endpoint derivatives:

    * both finite → :func:`cubic_hermite_poly` (``len(poly) == 4``);
    * either derivative divergent / undefined (inf/nan) →
      the linear polynomial ``[slope, v0]`` with ``slope = (v1 - v0) / h``
      (``len(poly) == 2``, numpy ``poly`` order).

    Callers that need to know which branch was taken can inspect
    ``len(poly)``; no separate flag is returned.  As in
    :func:`cubic_hermite_poly`, ``h`` must be positive (zero-width intervals
    are short-circuited by the callers).
    """
    if np.isfinite(dv0) and np.isfinite(dv1):
        return cubic_hermite_poly(h, v0, dv0, v1, dv1)
    return np.array([(v1 - v0) / h, v0])
