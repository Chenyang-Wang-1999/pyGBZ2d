'''
Pluggable Laurent-polynomial backends for :class:`bfgbz2d.core.CharPoly`.

Two implementations of the same protocol:

  - :class:`PolyToolsLaurent` — wraps the compiled ``poly_tools`` C++
    extension (fast; must be installed separately).  Imported lazily so the
    package itself only needs numpy + scipy.
  - :class:`NumpyLaurent` — pure-numpy fallback with identical semantics
    (validated against the C++ backend to ~1e-14; ``solve_roots_1d`` costs
    ~25% more because ``np.roots`` dominates either way).

Backend selection (first match wins):

  1. explicit ``CharPoly(coeffs, degs, backend=...)`` argument — a name
     (``'poly_tools'`` / ``'numpy'``) or any callable/class satisfying the
     protocol;
  2. the ``POLY_BACKEND`` environment variable (same names);
  3. ``poly_tools`` if importable;
  4. numpy fallback (with a one-time warning).

Protocol contract (what CharPoly consumes):

  - ``dim``: number of variables (3: E, beta1, beta2);
  - ``denom_orders``: per-variable denominator order — the Laurent
    denominator only clears NEGATIVE powers: ``denom_orders[d] =
    max(0, -min_degree_d)``;
  - ``num_max_degrees()``: per-variable max degree of the CLEARED
    numerator (denominator multiplied out) — CharPoly derives the minor
    degrees (M, N) from it;
  - ``eval(var)``: complex value at the variable tuple;
  - ``derivative(var_id)``: partial-derivative polynomial (same protocol);
  - ``partial_terms_1d(param_vals, param_indices, var_indices)``:
    substitute the parameter variables, return ``(coeffs, degs, denom)``
    for the remaining 1-D Laurent polynomial where *coeffs* are the
    numerator term coefficients at (cleared) *degs* and *denom* is the
    1-D denominator order.  Same-degree terms MUST be merged and exact
    zero coefficients dropped — downstream coefficient assembly writes
    one coefficient per degree.
'''

from __future__ import annotations

import os
import warnings
from typing import List, Optional, Protocol, Sequence, Tuple

import numpy as np

__all__ = [
    "LaurentProtocol",
    "PolyToolsLaurent",
    "NumpyLaurent",
    "make_laurent",
]

BACKEND_ENV_VAR = "POLY_BACKEND"

# Tri-state cache for the poly_tools import probe: None = not probed yet.
_POLY_TOOLS_AVAILABLE: Optional[bool] = None

# The numpy fallback warning fires once per process.
_numpy_fallback_warned: bool = False


class LaurentProtocol(Protocol):
    """Structural interface every backend must satisfy (see module docstring)."""

    dim: int
    denom_orders: Sequence[int]

    def eval(self, var) -> complex: ...

    def derivative(self, var_id: int) -> "LaurentProtocol": ...

    def num_max_degrees(self) -> Tuple[int, ...]: ...

    def partial_terms_1d(
        self, param_vals, param_indices, var_indices,
    ) -> Tuple[List[complex], List[int], int]: ...


# ---------------------------------------------------------------------------
# poly_tools (C++) backend
# ---------------------------------------------------------------------------

class PolyToolsLaurent:
    """Laurent polynomial on the compiled ``poly_tools`` extension."""

    def __init__(self, coeffs, degs):
        import poly_tools as pt  # lazy: keeps numpy+scipy the only hard deps

        self._pt = pt
        cl = pt.CLaurent(3)
        cl.set_Laurent_by_terms(
            pt.CScalarVec(np.asarray(coeffs, dtype=complex)),
            pt.CLaurentIndexVec(np.asarray(degs, dtype=int).flatten()),
        )
        self._cl = cl
        self.dim: int = cl.dim
        self.denom_orders: Sequence[int] = list(cl.denom_orders)

        # Per-dim max degree of the cleared numerator (batch_get_data returns
        # a flat per-term degree list; the d-th variable's degrees are every
        # 3rd entry).
        c_ct, d_ct = pt.CScalarVec([]), pt.CIndexVec([])
        cl.num.batch_get_data(c_ct, d_ct)
        d_flat = list(d_ct)
        self._num_max: Tuple[int, ...] = (
            tuple(max(d_flat[d::3]) for d in range(self.dim))
            if d_flat else (0,) * self.dim
        )

    # -- protocol surface --------------------------------------------------

    def num_max_degrees(self) -> Tuple[int, ...]:
        return self._num_max

    def eval(self, var) -> complex:
        return self._cl.eval(self._pt.CScalarVec(var))

    def derivative(self, var_id: int) -> "PolyToolsLaurent":
        return self._wrap(self._cl.derivative(var_id))

    def partial_terms_1d(self, param_vals, param_indices, var_indices):
        pt = self._pt
        poly_1d = self._cl.partial_eval(
            pt.CScalarVec(param_vals),
            pt.CIndexVec(param_indices),
            pt.CIndexVec(var_indices),
        )
        coeffs_ct = pt.CScalarVec([])
        degs_ct = pt.CIndexVec([])
        poly_1d.num.batch_get_data(coeffs_ct, degs_ct)
        degs_list = list(degs_ct)
        if not degs_list:
            # A completely vanished partial polynomial: the actual 1-D
            # denominator degree is zero (mirrors the old inline logic).
            return [], [], 0
        return list(coeffs_ct), degs_list, int(poly_1d.denom_orders[0])

    # -- helpers ------------------------------------------------------------

    def _wrap(self, cl) -> "PolyToolsLaurent":
        other = PolyToolsLaurent.__new__(PolyToolsLaurent)
        other._pt = self._pt
        other._cl = cl
        other.dim = cl.dim
        other.denom_orders = list(cl.denom_orders)
        c_ct, d_ct = self._pt.CScalarVec([]), self._pt.CIndexVec([])
        cl.num.batch_get_data(c_ct, d_ct)
        d_flat = list(d_ct)
        other._num_max = (
            tuple(max(d_flat[d::3]) for d in range(cl.dim))
            if d_flat else (0,) * cl.dim
        )
        return other


# ---------------------------------------------------------------------------
# Pure-numpy backend
# ---------------------------------------------------------------------------

class NumpyLaurent:
    """Sparse multivariate Laurent polynomial: ``sum_i c_i * z^degs[i]``.

    Semantics mirror the C++ backend exactly (validated by
    ``tests/test_backend.py``):

      * exact-zero coefficients are dropped at construction;
      * ``denom_orders[d] = max(0, -min_deg_d)`` — only negative powers are
        cleared;
      * ``partial_terms_1d`` merges same-degree terms (the C++ linked list
        merges them too) and drops terms that cancel to exactly zero.
    """

    def __init__(self, coeffs, degs):
        coeffs = np.asarray(coeffs, dtype=complex)
        degs = np.asarray(degs, dtype=int)
        if degs.ndim == 1:
            degs = degs.reshape(-1, 1)
        keep = coeffs != 0
        self.coeffs = coeffs[keep].copy()
        self.degs = degs[keep].copy()
        self.dim: int = self.degs.shape[1]
        if len(self.degs):
            self.denom_orders: Sequence[int] = np.maximum(
                0, -self.degs.min(axis=0)).astype(int).tolist()
        else:
            self.denom_orders = [0] * self.dim

    # -- protocol surface --------------------------------------------------

    def num_max_degrees(self) -> Tuple[int, ...]:
        if not len(self.degs):
            return (0,) * self.dim
        cleared = self.degs + np.asarray(self.denom_orders, dtype=int)[None, :]
        return tuple(int(v) for v in cleared.max(axis=0))

    def derivative(self, var_id: int) -> "NumpyLaurent":
        d = self.degs[:, var_id]
        keep = d != 0
        nd = self.degs[keep].copy()
        nd[:, var_id] -= 1
        return NumpyLaurent(self.coeffs[keep] * d[keep], nd)

    def eval(self, var) -> complex:
        p = self.coeffs
        for k in range(self.dim):
            p = p * (var[k] ** self.degs[:, k])
        return complex(p.sum())

    def partial_terms_1d(self, param_vals, param_indices, var_indices):
        coeff = self.coeffs.copy()
        for p_id, p_val in zip(param_indices, param_vals):
            coeff = coeff * (p_val ** self.degs[:, p_id])
        free = self.degs[:, list(var_indices)][:, 0]

        # Merge same-degree terms: np.unique sorts, add.at accumulates, then
        # exact zeros are dropped (an all-cancelled polynomial reports no
        # terms at all).
        uniq, inv = np.unique(free, return_inverse=True)
        merged = np.zeros(len(uniq), dtype=complex)
        np.add.at(merged, inv, coeff)
        keep = merged != 0
        uniq, merged = uniq[keep], merged[keep]
        if len(uniq) == 0:
            return [], [], 0
        denom = int(max(0, -uniq.min()))
        return merged.tolist(), (uniq + denom).tolist(), denom


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

def _poly_tools_ok() -> bool:
    """True when ``import poly_tools`` works (result cached)."""
    global _POLY_TOOLS_AVAILABLE
    if _POLY_TOOLS_AVAILABLE is None:
        try:
            import poly_tools  # noqa: F401
            _POLY_TOOLS_AVAILABLE = True
        except Exception:
            _POLY_TOOLS_AVAILABLE = False
    return _POLY_TOOLS_AVAILABLE


def _warn_numpy_fallback_once() -> None:
    global _numpy_fallback_warned
    if not _numpy_fallback_warned:
        _numpy_fallback_warned = True
        warnings.warn(
            "poly_tools is not importable; using the pure-numpy Laurent "
            "backend (numerically equivalent, ~1.25x slower on root "
            "solving). Install poly_tools or set "
            f"{BACKEND_ENV_VAR}=numpy to silence this.",
            stacklevel=3,
        )


def make_laurent(coeffs, degs, backend=None):
    """Build a Laurent-polynomial object on the selected *backend*.

    ``backend``: ``None`` (auto), ``'poly_tools'``, ``'numpy'``, or a
    user-supplied class/callable ``backend(coeffs, degs) -> LaurentProtocol``.
    Auto order: ``POLY_BACKEND`` env var → poly_tools (if importable) →
    numpy fallback with a one-time warning.
    """
    if backend is None:
        backend = os.environ.get(BACKEND_ENV_VAR, "auto")

    if isinstance(backend, str):
        if backend == "numpy":
            return NumpyLaurent(coeffs, degs)
        if backend in ("auto", "poly_tools"):
            if _poly_tools_ok():
                return PolyToolsLaurent(coeffs, degs)
            if backend == "poly_tools":
                raise ImportError(
                    "backend='poly_tools' requested but the poly_tools "
                    "package is not importable"
                )
            _warn_numpy_fallback_once()
            return NumpyLaurent(coeffs, degs)
        raise ValueError(
            f"unknown backend {backend!r}; expected 'poly_tools', 'numpy' "
            f"or a callable"
        )

    # User-supplied backend: any class or factory callable.
    return backend(coeffs, degs)
