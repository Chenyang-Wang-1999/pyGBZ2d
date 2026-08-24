'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-24
Copyright © Department of Physics, Tsinghua University. All rights reserved

Fixed-(E_ref, mu1) GBZ + winding-loop debugging tools.

The solver entry points (:func:`brute_force_SGBZ.collect_GBZ_subsets` /
:func:`brute_force_amoeba.collect_GBZ_subsets`) *solve* for mu1, so when a
result looks wrong there is no way to inspect what the machinery sees at the
candidate mu1.  This module freezes mu1 instead and exposes, per method
("sgbz" / "amoeba"):

  1. :func:`collect_debug_subsets` — the GBZ subsets at the given
     (E_ref, mu1): PointSubsets / LineSubsets for BOTH methods in one call,
     plus the SGBZ crossing charges (``charge`` / ``kind`` dicts aligned with
     the PointSubsets) and the reference windings (SGBZ ``W(E_ref, mu1)``,
     amoeba ``w1`` / ``w2``).

  2. :func:`compute_loop_windings` — the loop winding number at a list of
     user-chosen theta2 values.  A "winding loop" at theta2 is the closed
     loop ``beta1 = exp(mu1 + i*theta1)`` with theta1 sweeping [0, 2pi) and
     beta2 on the method's own path:

       * SGBZ  : ``beta2 = exp(mu2_mid(theta1) + i*theta2)`` — the same
         piecewise-smooth boundary-pair mean the solver integrates over
         (reuses ``brute_force_SGBZ.winding._loop_winding_quad`` /
         ``_loop_min_f``, so the numbers here are bit-comparable with the
         solver's own winding evaluation).
       * amoeba: ``beta2 = exp(mu2 + i*theta2)`` fixed (mu2 = the w2=0
         solution).  The winding is counted EXACTLY by the argument
         principle (number of beta1 roots inside |beta1| < exp(mu1) minus
         the pole order M1) and independently cross-checked by a quadrature
         of Im[f'/f] along the loop.

     The winding profile u(theta2) jumps exactly at the subsets' theta2, by
     the SGBZ charge there — :func:`compute_loop_windings` verifies that on
     every gap between consecutive loops (``gap_checks``), which detects
     missed / duplicated / mis-charged crossings.

  3. :func:`plot_winding_debug` — the theta1-theta2 torus picture: subsets
     (SGBZ PointSubsets charge-coded by marker shape AND a text label),
     LineSubsets as curves, the winding loops as horizontal lines annotated
     with their loop winding numbers, and the MR rows as vertical guides.

Debug philosophy: nothing here fixes or filters the modules' output.  When
the charge-consistency check fails or the amoeba quad/root-count disagree,
that discrepancy is the bug signal this tool exists to surface.
'''

from __future__ import annotations

import math
from cmath import exp
from dataclasses import dataclass, field
from typing import Optional, Sequence, Union

import numpy as np

from gbz_types import (
    CharPoly, PointSubset, LineSubset, TWO_PI, circ_dist,
)

from brute_force_SGBZ import (
    Mu2MidZM, detect_crossings_simple, compute_average_winding,
    extract_continuum_linesubsets, CONTINUUM_TOL,
)
# Private winding helpers imported deliberately: the debug tool must show the
# SAME loop-winding evaluation the solver uses, not a re-implementation that
# could silently diverge from it.
from brute_force_SGBZ.winding import (
    _loop_winding_quad, _loop_min_f, WindingFun, get_winding_number,
)
from brute_force_SGBZ.mu2mid import _CROSSING_TOL
from brute_force_amoeba.bisect import _find_mu2_for_w2_zero
from brute_force_amoeba.zm_extract import AmoebaZeroManager, extract_amoeba_subsets
from brute_force_amoeba.ronkin_winding import _get_average_winding_from_zeros


__all__ = [
    "MethodDebug", "GBZDebugReport", "LoopWinding", "GapCheck", "MethodLoops",
    "collect_debug_subsets", "auto_theta2_grid", "compute_loop_windings",
    "plot_winding_debug",
]


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class MethodDebug:
    """Fixed-(E_ref, mu1) subsets of ONE method ("sgbz" | "amoeba").

    Attributes:
        poly / E_ref / mu1: the frozen evaluation point.
        zm: the built ZeroManager (Mu2MidZM for sgbz, AmoebaZeroManager for
            amoeba) — kept alive so loop windings reuse its tracks / mu2_mid
            path instead of rebuilding them.
        subsets: list[PointSubset | LineSubset] — the requested GBZ subsets.
        charges: SGBZ only — one dict per PointSubset (same order):
            ``theta1``, ``theta2``, ``charge`` (+1/-1 ordinary, None for
            mr/tangent hard boundaries), ``kind``.  Empty for amoeba (it has
            no charge concept) and for the SGBZ continuum case.
        mu2: amoeba only — the w2=0 mu2 the subsets live at.
        is_continuum: whether the method reported a continuum (1D) case.
        W: reference winding — SGBZ average major-axis winding W(E_ref, mu1);
            amoeba w1 at (mu1, mu2).  ``None`` when undefined (continuum).
        w_secondary: amoeba only — the w2 winding at the solved mu2 (≈ 0 by
            construction; a tuple of left/right limits when continuum).
        error: failure message; the report never raises per-method errors.
    """

    method: str
    E_ref: complex
    mu1: float
    poly: CharPoly
    zm: object = None
    subsets: list = field(default_factory=list)
    charges: list = field(default_factory=list)
    mu2: Optional[float] = None
    is_continuum: bool = False
    W: Optional[float] = None
    w_secondary: Optional[float] = None
    error: Optional[str] = None

    @property
    def success(self) -> bool:
        return self.error is None

    @property
    def point_subsets(self) -> list[PointSubset]:
        return [s for s in self.subsets if isinstance(s, PointSubset)]

    @property
    def line_subsets(self) -> list[LineSubset]:
        return [s for s in self.subsets if isinstance(s, LineSubset)]

    def summary(self) -> str:
        if not self.success:
            return f"[{self.method}] FAILED: {self.error}"
        lines = [
            f"[{self.method}] E_ref={self.E_ref}, mu1={self.mu1:.12g}"
            + (f", mu2={self.mu2:.12g}" if self.mu2 is not None else "")
            + f": continuum={self.is_continuum}, "
            f"{len(self.point_subsets)} points + {len(self.line_subsets)} lines"
            + (f", W={self.W:+.6g}" if self.W is not None else ", W=undefined"),
        ]
        if self.method == "amoeba" and self.w_secondary is not None:
            lines.append(f"        w2 at solved mu2 = {self.w_secondary}")
        for c in self.charges:
            q = c["charge"]
            q_txt = "?" if q is None else f"{q:+d}"
            lines.append(
                f"        point (theta1={c['theta1']:.6f}, theta2={c['theta2']:.6f})"
                f"  charge={q_txt}  kind={c['kind']}")
        return "\n".join(lines)


@dataclass
class GBZDebugReport:
    """Both methods' fixed-(E_ref, mu1) results in one container."""

    E_ref: complex
    mu1: float
    poly: CharPoly
    methods: dict = field(default_factory=dict)

    def __getitem__(self, method: str) -> MethodDebug:
        return self.methods[method]

    @property
    def sgbz(self) -> MethodDebug:
        return self.methods["sgbz"]

    @property
    def amoeba(self) -> MethodDebug:
        return self.methods["amoeba"]

    def summary(self) -> str:
        return "\n".join(m.summary() for m in self.methods.values())


@dataclass
class LoopWinding:
    """Loop winding number at one theta2 (one method's loop).

    Attributes:
        winding: the rounded integer loop winding number.
        winding_raw: the unrounded value (deviation from an integer is itself
            a numerical-reliability signal).
        min_abs_f: min |f| along the loop (SGBZ quad reliability metric; a
            mesh evaluation for amoeba).
        margin: amoeba only — smallest |ln|beta1-root| - mu1| at this theta2
            (root-count ambiguity metric; ~0 means a root sits ON the loop).
        winding_quad: amoeba only — the independent Im[f'/f] quadrature
            cross-check of the root-count winding.
        reliable: False when the loop passes too close to a zero of f (or the
            amoeba quad disagrees with the root count) — the winding may be
            off by ±1 there.
    """

    method: str
    theta2: float
    winding: int
    winding_raw: float
    reliable: bool
    min_abs_f: Optional[float] = None
    margin: Optional[float] = None
    winding_quad: Optional[float] = None
    note: str = ""


@dataclass
class GapCheck:
    """Charge-consistency verification on one gap between two loops.

    The winding profile u(theta2) must jump across a gap by exactly the sum
    of the SGBZ charges (or the number of amoeba points) strictly inside it.
    A mismatch means the crossing detection missed / duplicated / mis-charged
    a crossing — the primary bug signal this tool surfaces.
    """

    lo: float
    hi: float
    expected: Optional[int]
    actual: Optional[int]
    status: str  # 'ok' | 'mismatch' | 'skipped-hard' | 'skipped'
    detail: str = ""


@dataclass
class MethodLoops:
    """Loop windings + gap checks for one method.

    ``W_profile`` is the arc-weighted mean of u(theta2) reconstructed from
    the loops (charge-propagated across the crossings inside each gap) —
    comparable with the solver's own ``MethodDebug.W``.  ``None`` when the
    reconstruction is impossible (skipped/mismatched gap checks).
    """

    method: str
    loops: list = field(default_factory=list)
    gap_checks: list = field(default_factory=list)
    W_profile: Optional[float] = None

    def summary(self) -> str:
        out = [f"loops [{self.method}]:"]
        for lw in self.loops:
            flag = "" if lw.reliable else "  (UNRELIABLE)"
            extra = lw.note if lw.note else ""
            out.append(
                f"  theta2={lw.theta2:9.4f}  w={lw.winding:+d}"
                f"  (raw {lw.winding_raw:+.4f}){flag}{extra}")
        n_ok = sum(1 for g in self.gap_checks if g.status == "ok")
        n_bad = sum(1 for g in self.gap_checks if g.status == "mismatch")
        n_skip = len(self.gap_checks) - n_ok - n_bad
        out.append(f"  gap checks: {n_ok} ok, {n_bad} MISMATCH, {n_skip} skipped")
        for g in self.gap_checks:
            if g.status != "ok":
                out.append(f"    ({g.lo:.4f} -> {g.hi:.4f}) {g.status}: {g.detail}")
        if self.W_profile is not None:
            out.append(f"  W from loop profile = {self.W_profile:+.6f}")
        return "\n".join(out)


# ---------------------------------------------------------------------------
# Function 1: subsets (+ SGBZ charges) at fixed (E_ref, mu1)
# ---------------------------------------------------------------------------

def _collect_sgbz(
    poly: CharPoly, E_ref: complex, mu1: float,
    *, zm_run_kwargs: Optional[dict], continuum_tol: float,
    crossing_tol: float,
) -> MethodDebug:
    """SGBZ slice at fixed mu1 — mirrors ``_evaluate_winding`` (sgbz_solver).

    One Mu2MidZM build + analyze feeds everything: the continuum gate, the
    0D crossing materialization and the average winding.  No mu1 bisection
    and no plateau reclassification here — a debug slice must show the raw
    machinery output at the requested mu1.
    """
    zm = Mu2MidZM(poly, E_ref, mu1)
    zm.run(**(zm_run_kwargs or {}))
    zm.analyze(tie_tol=continuum_tol, crossing_tol=crossing_tol)

    if zm.has_continuum:
        subsets = extract_continuum_linesubsets(zm, poly)
        return MethodDebug(
            method="sgbz", E_ref=E_ref, mu1=mu1, poly=poly, zm=zm,
            subsets=list(subsets), charges=[], is_continuum=True, W=None,
        )

    subsets, charges = detect_crossings_simple(
        zm, poly, crossing_tol=crossing_tol,
    )
    W = compute_average_winding(zm, poly, charges)
    return MethodDebug(
        method="sgbz", E_ref=E_ref, mu1=mu1, poly=poly, zm=zm,
        subsets=list(subsets), charges=charges, W=float(W),
    )


def _collect_amoeba(
    poly: CharPoly, E_ref: complex, mu1: float,
    *, zm_run_kwargs: Optional[dict], mu2_guess: tuple,
    options: dict,
) -> MethodDebug:
    """Amoeba slice at fixed mu1 — ``collect_GBZ_subsets`` minus the outer
    mu1 bisection and the plateau reclassification.

    Solves mu2 for the w2=0 condition on the pre-built AmoebaZeroManager
    (the ZM is mu2-independent, so it is built once and reused, exactly like
    ``bisect_amoeba_ronkin_min``), extracts the subsets at (mu1, mu2) and
    reports w1 — the outer-loop quantity whose zero would make this mu1 the
    actual amoeba GBZ radius.
    """
    zm = AmoebaZeroManager(poly, E_ref, mu1)
    zm.run(**(zm_run_kwargs or {}))
    inner = _find_mu2_for_w2_zero(
        poly, E_ref, mu1, float(mu2_guess[0]), float(mu2_guess[1]),
        _zm=zm, **options,
    )
    mu2 = float(inner["mu2"])
    subsets = extract_amoeba_subsets(zm, poly, E_ref, mu1, mu2, mode="solve")

    if inner["is_continuum"]:
        w1 = None
    else:
        w1, _ = _get_average_winding_from_zeros(
            poly, E_ref, mu1, mu2, inner.get("zeros") or [], direction=1,
        )
    return MethodDebug(
        method="amoeba", E_ref=E_ref, mu1=mu1, poly=poly, zm=zm,
        subsets=list(subsets), charges=[], mu2=mu2,
        is_continuum=bool(inner["is_continuum"]),
        W=None if w1 is None else float(w1),
        w_secondary=inner.get("winding"),
    )


def collect_debug_subsets(
    poly: CharPoly,
    E_ref: complex,
    mu1: float,
    *,
    methods: Sequence[str] = ("sgbz", "amoeba"),
    zm_run_kwargs: Optional[dict] = None,
    mu2_guess: tuple = (-1.0, 1.0),
    continuum_tol: float = CONTINUUM_TOL,
    crossing_tol: float = _CROSSING_TOL,
    amoeba_options: Optional[dict] = None,
) -> GBZDebugReport:
    """Compute the GBZ subsets of every requested method at FIXED mu1.

    Parameters:
        poly: ``CharPoly(coeffs, degs)`` of the model under test.
        E_ref, mu1: the frozen evaluation point (mu1 is NOT solved for —
            inspect a solver candidate or probe any radius).
        methods: subset of ``("sgbz", "amoeba")``.
        zm_run_kwargs: forwarded to ``ZeroManager.run`` of both methods.
        mu2_guess: initial mu2 bracket for the amoeba w2=0 bisection.
        continuum_tol / crossing_tol: SGBZ analyze tunables.
        amoeba_options: extra kwargs for the amoeba's ``_find_mu2_for_w2_zero``
            (``continuum_tol``, ``continuum_perturb``, ``max_iter``, ``xtol``).

    Returns:
        GBZDebugReport; ``report[method].subsets`` is the subset list and
        ``report["sgbz"].charges`` carries the PointSubset charges.  A method
        that raises is captured in ``MethodDebug.error`` (never re-raised) so
        one broken module does not hide the other's evidence.
    """
    options = dict(continuum_tol=1e-6, continuum_perturb=1e-4,
                   max_iter=60, xtol=1e-10)
    options.update(amoeba_options or {})

    report = GBZDebugReport(E_ref=E_ref, mu1=mu1, poly=poly, methods={})
    for method in methods:
        try:
            if method == "sgbz":
                report.methods[method] = _collect_sgbz(
                    poly, E_ref, mu1, zm_run_kwargs=zm_run_kwargs,
                    continuum_tol=continuum_tol, crossing_tol=crossing_tol,
                )
            elif method == "amoeba":
                report.methods[method] = _collect_amoeba(
                    poly, E_ref, mu1, zm_run_kwargs=zm_run_kwargs,
                    mu2_guess=mu2_guess, options=options,
                )
            else:
                raise ValueError(f"unknown method {method!r}")
        except Exception as exc:  # noqa: BLE001 — debug capture, by design
            report.methods[method] = MethodDebug(
                method=method, E_ref=E_ref, mu1=mu1, poly=poly,
                error=f"{type(exc).__name__}: {exc}",
            )
    return report


# ---------------------------------------------------------------------------
# Auto theta2 grid ("find the winding loops")
# ---------------------------------------------------------------------------

def point_theta2(sub: PointSubset) -> float:
    """theta2 = arg(beta2) wrapped to [0, 2pi) (PointSubset stores only beta2)."""
    return float(np.angle(sub.beta2) % TWO_PI)


def subset_theta2_values(md: MethodDebug, *, max_line_samples: int = 16) -> list[float]:
    """theta2 coordinates to steer loops away from.

    PointSubset theta2s are the exact jump locations of the winding profile,
    so loops belong in the gaps between them.  LineSubsets contribute a
    subsample of their theta2 track (a loop meeting a continuum track passes
    through / near a root of f, so its winding is unreliable) — the track is
    sampled, not exhaustively excluded; the reliability metrics catch the
    rest.
    """
    values = [c["theta2"] for c in md.charges]
    if not md.charges:  # amoeba / SGBZ-continuum: read the subsets directly
        values = [point_theta2(p) for p in md.point_subsets]
    for line in md.line_subsets:
        t2 = np.angle(line.beta2_arr) % TWO_PI
        step = max(1, len(t2) // max_line_samples)
        values.extend(float(v) for v in t2[::step])
    return sorted(set(round(v, 12) for v in values))


def auto_theta2_grid(
    md: MethodDebug,
    *,
    n_per_gap: int = 1,
    min_gap: float = 0.05,
    max_loops: int = 8,
) -> list[float]:
    """Pick loop theta2s at the midpoints of the widest subset-theta2 gaps.

    Midpoints maximize the distance from every profile jump, which is where
    the loop stays farthest from the zeros of f (the same "safest seed"
    idea as ``brute_force_SGBZ.winding._pick_seed_theta2``).  The *widest*
    gaps are preferred (capped at *max_loops*) so a densely sampled
    continuum line does not flood the plot with dozens of loops.

    With no usable gap at all (no subsets, or a continuum track covering the
    whole circle), the quadrant midpoints are returned.
    """
    occupied = subset_theta2_values(md)
    quadrants = [float(a * np.pi / 4.0) for a in (1, 3, 5, 7)]
    if not occupied:
        return quadrants

    n = len(occupied)
    gaps = []
    for i in range(n):
        a = occupied[i]
        b = occupied[(i + 1) % n]
        width = (b - a) % TWO_PI
        if width > min_gap:
            gaps.append((width, a, b))
    if not gaps:
        return quadrants

    gaps.sort(reverse=True)
    grid: list[float] = []
    for width, a, _b in gaps[:max_loops]:
        for k in range(n_per_gap):
            t2 = (a + (k + 1) * width / (n_per_gap + 1)) % TWO_PI
            grid.append(float(t2))
    return sorted(grid)


# ---------------------------------------------------------------------------
# Function 2a: loop winding numbers
# ---------------------------------------------------------------------------

def _amoeba_loop_root_count(
    poly: CharPoly, E_ref: complex, mu1: float, mu2: float, theta2: float,
) -> tuple[int, float]:
    """Exact loop winding via the argument principle in the beta1 plane.

    The loop ``beta1 = exp(mu1 + i*theta1)``, ``beta2 = exp(mu2 + i*theta2)``
    fixed, has winding ``Z - P`` where Z counts the beta1-roots of
    f(E, ., beta2) inside ``|beta1| < exp(mu1)`` and P = M1 is the pole
    order at beta1 = 0 — the identical computation
    ``_get_average_winding_from_zeros`` performs for its segment samples, so
    the debug numbers stay comparable with the solver's w1.

    Returns ``(winding, margin)`` where *margin* is the smallest
    ``|ln|root| - mu1|`` — the root-count ambiguity measure.
    """
    M1, N1 = poly.get_minor_degrees(1)
    beta2 = exp(mu2 + 1j * theta2)
    roots = poly.solve_roots_1d((0, 2), (E_ref, beta2), (1,), M1, N1)
    la = np.log(np.abs(np.asarray(roots, dtype=complex)))
    winding = int(np.sum(la < mu1)) - M1
    finite = la[np.isfinite(la)]
    margin = float(np.min(np.abs(finite - mu1))) if finite.size else math.inf
    return winding, margin


def _amoeba_loop_quad(
    poly: CharPoly, E_ref: complex, mu1: float, mu2: float, theta2: float,
) -> float:
    """Independent cross-check: quadrature of Im[f'/f] along the loop."""
    beta2 = exp(mu2 + 1j * theta2)

    def loop_fun(t: float):
        beta1 = exp(mu1 + 1j * t)
        return (E_ref, beta1, beta2), (0.0, 1j * beta1, 0.0)

    return get_winding_number(WindingFun(poly, loop_fun, (0.0, TWO_PI)),
                              n_seg=4)


def _amoeba_loop_min_abs_f(
    poly: CharPoly, E_ref: complex, mu1: float, mu2: float, theta2: float,
    n_mesh: int = 256,
) -> float:
    """min |f| on a theta1 mesh along the fixed-beta2 loop."""
    beta2 = exp(mu2 + 1j * theta2)
    t = np.linspace(0.0, TWO_PI, n_mesh, endpoint=False)
    beta1 = np.exp(mu1 + 1j * t)
    vals = [abs(poly.eval_val((E_ref, complex(b1), beta2))) for b1 in beta1]
    return float(min(vals))


def _cyclic_strictly_inside(x: float, a: float, b: float, eps: float = 1e-9) -> bool:
    """Whether *x* lies strictly inside the CCW-open circle interval (a, b)."""
    span = (b - a) % TWO_PI
    dx = (x - a) % TWO_PI
    return bool(eps < dx < span - eps)


def _gap_arc_contribution(
    w_lo: int, lo: float, hi: float,
    inside_marks: list[tuple[float, int]],
) -> float:
    """Exact ∫u dθ₂ over one gap, charge-propagated across its crossings.

    Within the gap the profile starts at *w_lo* and each enclosed crossing
    ``(theta2, charge)`` — in circular order after *lo* — adds its charge to
    the remainder of the gap.  This reconstructs the solver's arc-weighted
    mean from the loop grid even when a gap contains cancelling crossings
    (which a naive ``w * width`` read would miss).
    """
    width = (hi - lo) % TWO_PI
    if width <= 0.0:
        width = TWO_PI
    end = lo + width
    total = float(w_lo) * width
    for t2, q in inside_marks:
        u = lo + ((t2 - lo) % TWO_PI)   # unwrapped position inside the gap
        total += q * (end - u)
    return total


def _gap_checks_sgbz(
    md: MethodDebug, loops: list[LoopWinding],
) -> tuple[list[GapCheck], Optional[float]]:
    """Verify that consecutive loop windings differ by the enclosed charges.

    Skipped on the continuum case (no charges exist) and on gaps containing
    a HARD boundary (charge unknown by design — ``compute_average_winding``
    recomputes the winding independently across it).

    Returns ``(checks, W_profile)`` where *W_profile* is the arc-weighted
    mean of u(θ₂) reconstructed from the loop grid (``None`` when any gap
    check is not "ok", since the reconstruction would build on broken data).
    """
    if len(loops) < 2 or md.is_continuum or not md.charges:
        return [GapCheck(0.0, TWO_PI, None, None, "skipped",
                         "needs >=2 loops and discrete charges")], None

    checks: list[GapCheck] = []
    arc_total = 0.0
    n = len(loops)
    for i in range(n):
        lo = loops[i].theta2
        hi = loops[(i + 1) % n].theta2
        inside = [c for c in md.charges
                  if _cyclic_strictly_inside(c["theta2"], lo, hi)]
        hard = [c for c in inside if c["charge"] is None]
        if hard:
            checks.append(GapCheck(
                lo, hi, None, None, "skipped-hard",
                f"{len(hard)} hard boundary(ies) inside",
            ))
            continue
        expected = int(sum(c["charge"] for c in inside))
        actual = loops[(i + 1) % n].winding - loops[i].winding
        status = "ok" if expected == actual else "mismatch"
        enclosed_txt = ",".join(f"{c['charge']:+d}" for c in inside) or "-"
        detail = (f"enclosed charges [{enclosed_txt}]"
                  if status == "ok" else
                  f"charge sum {expected:+d} but winding jump {actual:+d} "
                  f"({len(inside)} enclosed crossings)")
        checks.append(GapCheck(lo, hi, expected, actual, status, detail))
        if status == "ok":
            arc_total += _gap_arc_contribution(
                loops[i].winding, lo, hi,
                [(c["theta2"], c["charge"]) for c in inside],
            )
    w_profile = (arc_total / TWO_PI
                 if all(g.status == "ok" for g in checks) else None)
    return checks, w_profile


def _amoeba_point_jumps(
    md: MethodDebug,
    occupied_theta2: Sequence[float] = (),
    *,
    probe_eps: float = 1e-3,
) -> dict[float, int]:
    """Signed jump of u(theta2) at each amoeba PointSubset theta2.

    The amoeba module has no charge concept, so the analog of the SGBZ
    charge is obtained numerically: evaluate the exact loop winding
    (:func:`_amoeba_loop_root_count`) just left and just right of every
    point.  Two nearby points with opposite jumps (+1/−1) legitimately
    cancel inside one gap — the gap check must sum SIGNS, not count points.

    *probe_eps* is shrunk when it would reach a neighbouring occupied
    theta2 (another point, or a loop).
    """
    pts = [point_theta2(p) for p in md.point_subsets]
    all_marks = sorted(set([round(v, 12) for v in pts]
                           + [round(v, 12) for v in occupied_theta2]))
    jumps: dict[float, int] = {}
    for t2 in pts:
        eps = probe_eps
        for other in all_marks:
            d = circ_dist(t2, other)
            if 0.0 < d < 2.0 * eps:
                eps = 0.5 * d
        w_lo, _ = _amoeba_loop_root_count(
            md.poly, md.E_ref, md.mu1, md.mu2, (t2 - eps) % TWO_PI)
        w_hi, _ = _amoeba_loop_root_count(
            md.poly, md.E_ref, md.mu1, md.mu2, (t2 + eps) % TWO_PI)
        jumps[round(t2, 12)] = int(w_hi - w_lo)
    return jumps


def _gap_checks_amoeba(
    md: MethodDebug, loops: list[LoopWinding],
) -> tuple[list[GapCheck], Optional[float]]:
    """Amoeba variant: the winding jump must equal the enclosed signed jumps.

    Only meaningful when every subset is a PointSubset: a LineSubset crosses
    horizontal loops at interior theta2s, so the point jumps cannot predict
    the profile — skip rather than produce a false mismatch.

    Returns ``(checks, W_profile)`` (see :func:`_gap_checks_sgbz`).
    """
    if len(loops) < 2 or md.line_subsets:
        return [GapCheck(0.0, TWO_PI, None, None, "skipped",
                         "needs >=2 loops and no LineSubsets")], None

    loop_t2 = [lw.theta2 for lw in loops]
    jumps = _amoeba_point_jumps(md, loop_t2)
    checks: list[GapCheck] = []
    arc_total = 0.0
    n = len(loops)
    for i in range(n):
        lo = loops[i].theta2
        hi = loops[(i + 1) % n].theta2
        inside = {t2: q for t2, q in jumps.items()
                  if _cyclic_strictly_inside(t2, lo, hi)}
        expected = int(sum(inside.values()))
        actual = loops[(i + 1) % n].winding - loops[i].winding
        status = "ok" if expected == actual else "mismatch"
        jumps_txt = ",".join(f"{q:+d}" for q in inside.values()) or "-"
        detail = (f"enclosed point jumps [{jumps_txt}]"
                  if status == "ok" else
                  f"point jumps sum {expected:+d} but winding jump {actual:+d}")
        checks.append(GapCheck(lo, hi, expected, actual, status, detail))
        if status == "ok":
            arc_total += _gap_arc_contribution(
                loops[i].winding, lo, hi, list(inside.items()),
            )
    w_profile = (arc_total / TWO_PI
                 if all(g.status == "ok" for g in checks) else None)
    return checks, w_profile


def compute_loop_windings(
    md: MethodDebug,
    theta2_list: Sequence[float],
    *,
    reliable_f_tol: float = 1e-8,
    reliable_margin_tol: float = 1e-6,
    quad_cross_check: bool = True,
) -> MethodLoops:
    """Compute the loop winding number at each theta2 in *theta2_list*.

    SGBZ loops ride the mu2_mid path (the solver's own quad); amoeba loops
    use the exact beta1 root count with an optional independent quad
    cross-check.  Every loop also reports its reliability metric (min |f|
    along the loop / root-count margin).

    On top of the per-loop numbers, the winding profile u(theta2) is
    verified gap-by-gap between consecutive loops (sorted by theta2): the
    jump across each gap must equal the enclosed SGBZ charges / amoeba
    signed point jumps.  When every gap check passes, the arc-weighted mean
    of u(theta2) is reconstructed charge-propagated into ``W_profile`` —
    directly comparable with the solver's own ``MethodDebug.W``.

    Returns a :class:`MethodLoops` (``.loops`` + ``.gap_checks`` +
    ``.W_profile``).
    """
    if not md.success:
        raise ValueError(f"method {md.method!r} failed: {md.error}")

    loops: list[LoopWinding] = []
    for t2 in theta2_list:
        if md.method == "sgbz":
            raw = _loop_winding_quad(md.zm, md.poly, md.E_ref, md.mu1, float(t2))
            min_f = _loop_min_f(float(t2), md.zm, md.poly)
            reliable = min_f > reliable_f_tol
            note = ""
            if md.is_continuum:
                # The collapsed boundary pair makes the mu2_mid loop touch a
                # root track over whole theta1 intervals — quad results are
                # indicative only.
                note = "  [continuum: mu2_mid touches the tracks]"
                reliable = False
            loops.append(LoopWinding(
                method="sgbz", theta2=float(t2),
                winding=int(round(raw)), winding_raw=float(raw),
                reliable=reliable, min_abs_f=float(min_f), note=note,
            ))
        elif md.method == "amoeba":
            w, margin = _amoeba_loop_root_count(
                md.poly, md.E_ref, md.mu1, md.mu2, float(t2))
            quad = (_amoeba_loop_quad(md.poly, md.E_ref, md.mu1, md.mu2, float(t2))
                    if quad_cross_check else None)
            min_f = _amoeba_loop_min_abs_f(
                md.poly, md.E_ref, md.mu1, md.mu2, float(t2))
            reliable = margin > reliable_margin_tol
            note = ""
            if quad is not None and int(round(quad)) != w:
                reliable = False
                note = f"  [quad cross-check {quad:+.4f} != root count {w:+d}]"
            loops.append(LoopWinding(
                method="amoeba", theta2=float(t2),
                winding=w, winding_raw=float(w),
                reliable=reliable, min_abs_f=min_f, margin=margin,
                winding_quad=None if quad is None else float(quad),
                note=note,
            ))
        else:
            raise ValueError(f"unknown method {md.method!r}")

    # Gap checks read the loops in sorted-theta2 order so "consecutive"
    # means adjacent on the circle.
    loops_sorted = sorted(loops, key=lambda lw: lw.theta2)
    if md.method == "sgbz":
        gap_checks, w_profile = _gap_checks_sgbz(md, loops_sorted)
    else:
        gap_checks, w_profile = _gap_checks_amoeba(md, loops_sorted)
    return MethodLoops(method=md.method, loops=loops, gap_checks=gap_checks,
                       W_profile=w_profile)


# ---------------------------------------------------------------------------
# Function 2b: the theta1-theta2 debug plot
# ---------------------------------------------------------------------------

# Per-method visual identity: SGBZ warm (crimson/orange), amoeba cold
# (teal/green); charge signs share one hue family across methods.
_STYLE = {
    "sgbz": {
        "point": {"+1": ("^", "#d62728"), "-1": ("v", "#1f77b4"),
                  None: ("s", "black"), "+0": ("o", "black"), "-0": ("o", "black")},
        "line": dict(color="darkorange", lw=1.8),
        "loop": dict(color="crimson", ls=(0, (6, 3)), lw=1.4),
        "label": dict(color="crimson", va="bottom"),
        "mr": dict(color="purple", ls="-.", lw=0.9, alpha=0.35),
    },
    "amoeba": {
        "point": {"+1": ("o", "#2ca02c"), "-1": ("o", "#2ca02c"),
                  None: ("o", "#2ca02c")},
        "line": dict(color="teal", lw=1.8),
        "loop": dict(color="teal", ls=(0, (2, 2)), lw=1.4),
        "label": dict(color="teal", va="top"),
        "mr": dict(color="purple", ls="-.", lw=0.9, alpha=0.35),
    },
}


def _split_torus_segments(
    theta1: np.ndarray, theta2: np.ndarray, jump: float = math.pi,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Split a wrapped curve into segments continuous on the [0, 2pi) torus.

    Joined LineSubsets legitimately wrap through the theta1 = 0 ≡ 2pi seam
    (and theta2 wraps as beta2's argument circles zero), so a raw mod-2pi
    plot would draw spurious lines across the whole torus.
    """
    t1 = np.asarray(theta1, dtype=float) % TWO_PI
    t2 = np.asarray(theta2, dtype=float) % TWO_PI
    if len(t1) == 0:
        return []
    breaks = [0] + [k for k in range(1, len(t1))
                    if circ_dist(t1[k], t1[k - 1]) > jump
                    or circ_dist(t2[k], t2[k - 1]) > jump] + [len(t1)]
    return [(t1[a:b], t2[a:b]) for a, b in zip(breaks[:-1], breaks[1:]) if b - a > 1]


def _plot_method(
    ax, md: MethodDebug, ml: MethodLoops, *,
    annotate_charges: bool, mark_mrs: bool,
) -> None:
    style = _STYLE[md.method]

    # --- MR rows (both ZM flavours expose .multiple_roots) ---
    if mark_mrs and md.zm is not None:
        for mr in getattr(md.zm, "multiple_roots", []) or []:
            ax.axvline(float(mr.theta1) % TWO_PI, **style["mr"])

    # --- LineSubsets: theta1 vs arg(beta2) curves on the torus ---
    for line in md.line_subsets:
        segs = _split_torus_segments(line.theta1_arr,
                                     np.angle(line.beta2_arr))
        for t1_seg, t2_seg in segs:
            ax.plot(t1_seg, t2_seg, **style["line"])

    # --- PointSubsets ---
    if md.method == "sgbz":
        # Charge coding: marker SHAPE by sign (triangle up/down, square for
        # the unknown hard charge) + an explicit text label, so the value
        # survives both grayscale printing and marker overlap.
        for sub, ch in zip(md.subsets, md.charges):
            q = ch["charge"]
            key = None if q is None else f"{q:+d}"
            marker, color = style["point"][key]
            ax.plot(sub.theta1, point_theta2(sub), marker=marker, color=color,
                    ms=9, ls="none", mec="black", mew=0.6, zorder=6)
            if annotate_charges:
                label = "?" if q is None else f"{q:+d}"
                ax.annotate(label, (sub.theta1, point_theta2(sub)),
                            xytext=(4, 4), textcoords="offset points",
                            fontsize=8, color=color, zorder=7,
                            fontweight="bold")
    else:
        pts = md.point_subsets
        if pts:
            ax.plot([p.theta1 for p in pts], [point_theta2(p) for p in pts],
                    marker="o", color=style["point"]["+1"][1], ms=8,
                    ls="none", mfc="none", mec="black", mew=1.2, zorder=5)

    # --- winding loops: horizontal lines + winding-number labels ---
    # SGBZ labels sit ABOVE their line, amoeba BELOW, so the two methods'
    # annotations never collide on a shared theta2 grid.
    for lw in ml.loops:
        line_kw = dict(style["loop"])
        if not lw.reliable:
            line_kw.update(ls=":", alpha=0.35)
        ax.plot([0.0, TWO_PI], [lw.theta2, lw.theta2], **line_kw)
        txt = f"{md.method[:2]}: w={lw.winding:+d}" + ("" if lw.reliable else " ?")
        ax.annotate(txt, (0.02, lw.theta2), xycoords="data",
                    xytext=(3, 3 if style["label"]["va"] == "bottom" else -3),
                    textcoords="offset points", fontsize=8,
                    color=style["loop"]["color"],
                    va=style["label"]["va"], ha="left", zorder=7)


def _legend_handles(methods: list) -> list:
    from matplotlib.lines import Line2D
    handles = []
    for md in methods:
        st = _STYLE[md.method]
        if md.method == "sgbz":
            for key, name in (("+1", "SGBZ point, charge +1"),
                              ("-1", "SGBZ point, charge −1"),
                              (None, "SGBZ point, charge ? (hard)")):
                marker, color = st["point"][key]
                handles.append(Line2D([], [], ls="none", marker=marker,
                                      color=color, mec="black", mew=0.6,
                                      ms=8, label=name))
        else:
            handles.append(Line2D([], [], ls="none", marker="o",
                                  color=st["point"]["+1"][1], mfc="none",
                                  mec="black", mew=1.2, ms=8,
                                  label="amoeba point"))
        handles.append(Line2D([], [], **st["line"],
                              label=f"{md.method} LineSubset"))
        handles.append(Line2D([], [], **st["loop"],
                              label=f"{md.method} winding loop (label = w)"))
    return handles


def plot_winding_debug(
    target: Union[MethodDebug, GBZDebugReport, Sequence],
    theta2_list: Optional[Sequence[float]] = None,
    *,
    ax=None,
    annotate_charges: bool = True,
    mark_mrs: bool = True,
    auto_theta2_when_none: bool = True,
    n_per_gap: int = 1,
    quad_cross_check: bool = True,
    verbose: bool = True,
    title: Optional[str] = None,
    figsize: tuple = (9.5, 7.0),
):
    """Compute + plot the winding loops and subsets on the theta1-theta2 torus.

    Parameters:
        target: a :class:`MethodDebug`, a :class:`GBZDebugReport`, or a
            sequence of either — every successful method is plotted on the
            one shared torus.
        theta2_list: the loops' theta2 values (shared by all methods).  When
            ``None`` and *auto_theta2_when_none*, each method gets
            :func:`auto_theta2_grid` midpoints of its subset-theta2 gaps.
        annotate_charges: draw the SGBZ charge next to each PointSubset
            (marker shape/colour codes it regardless).
        mark_mrs: vertical guides at the ZeroManager's multiple-root rows.
        verbose: print the per-method loop / gap-check summaries.

    Returns:
        ``(fig, ax, results)`` where *results* maps method name ->
        :class:`MethodLoops`.
    """
    import matplotlib.pyplot as plt

    if isinstance(target, MethodDebug):
        methods = [target]
    elif isinstance(target, GBZDebugReport):
        methods = list(target.methods.values())
    else:
        methods = [t if isinstance(t, MethodDebug) else t
                   for t in target]

    ok = []
    for m in methods:
        if isinstance(m, MethodDebug) and m.success:
            ok.append(m)
        elif verbose:
            print(f"[{getattr(m, 'method', '?')}] skipped (failed): "
                  f"{getattr(m, 'error', 'not a MethodDebug')}")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    results: dict = {}
    for md in ok:
        t2_list = theta2_list
        if t2_list is None:
            if not auto_theta2_when_none:
                raise ValueError("theta2_list is required when "
                                 "auto_theta2_when_none=False")
            t2_list = auto_theta2_grid(md, n_per_gap=n_per_gap)
        ml = compute_loop_windings(md, t2_list,
                                   quad_cross_check=quad_cross_check)
        results[md.method] = ml
        _plot_method(ax, md, ml, annotate_charges=annotate_charges,
                     mark_mrs=mark_mrs)
        if verbose:
            print(md.summary())
            print(ml.summary())

    ax.set_xlim(-0.15, TWO_PI + 0.15)
    ax.set_ylim(-0.3, TWO_PI + 0.3)
    ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, TWO_PI])
    ax.set_xticklabels(["0", "π/2", "π", "3π/2", "2π"])
    ax.set_yticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, TWO_PI])
    ax.set_yticklabels(["0", "π/2", "π", "3π/2", "2π"])
    ax.set_xlabel(r"$\theta_1$")
    ax.set_ylabel(r"$\theta_2$")
    ax.grid(True, alpha=0.2)
    if title is None:
        if isinstance(target, GBZDebugReport):
            title = (f"GBZ debug: E={target.E_ref}, "
                     f"μ₁={target.mu1:.8g}")
        else:
            m0 = ok[0] if ok else None
            title = (f"GBZ debug [{m0.method}]: E={m0.E_ref}, "
                     f"μ₁={m0.mu1:.8g}") if m0 else "GBZ debug"
    ax.set_title(title, fontsize=11)
    if ok:
        ax.add_artist(ax.legend(handles=_legend_handles(ok),
                                loc="upper left", fontsize=8,
                                bbox_to_anchor=(1.01, 1.0)))
    return fig, ax, results
