# 2026-08-13 — SGBZ μ₂_mid refactor (partial: structure done, continuum-bisection blocked)

Refactored `brute_force_SGBZ/` per `log/2026-08-13-SGBZ算法梳理.md` (μ₂_mid as
a first-class piecewise-smooth object). Plan: `.claude/plans/refactor-sgbz-mu2mid.md`.

## Done

- **`mu2mid.py`** (new): `Mu2MidZM(ZeroManager)` — production port of the
  `demos/demo_zm_gbz.py` prototype. `ItemView` (representative-item abstraction
  collapsing continuum clusters, stable where `abs_argsort` is not),
  `Mu2MidBreakpoint`, cubic-Hermite helpers, the four build stages (walls →
  sort-change detection → cubic-Hermite bracketing iteration → assemble +
  breakpoints), and inline `has_continuum` (§1: j_lo==j_hi ⇒ 1D subset).
- **`crossings.py`**: per-column vs μ₂_mid (§2.1, unified with amoeba structure
  — μ₂ is the curve μ₂_mid(θ₁), not a constant). Cubic-Hermite bracketing on
  `g = ln|β_j| − μ₂_mid` (§2.3, replaces Newton); charge = sign(g') inline
  (§3.1). Deleted the pairwise sweep + `NotImplementedError` dedup (§2.5).
- **`winding.py`**: quad of `Im[f'/f]` over the piecewise-smooth μ₂_mid loop
  path (§6.4): `WindingFun` + `get_winding_number` (split at μ₂_mid breakpoints)
  + `_Mu2MidPath` interpolator (per-row cubic Hermite, left/right derivs at
  breakpoints, linear near MR). §3.2 region/seed/charge-propagation kept.
  **Validated:** matches the old mesh-row winding to machine precision at 0D μ₁.
- **`continuum_lines.py`** (new): LineSubset materialization — the
  `j_lo==j_hi` cluster tracks, joined across MRs via amoeba's
  `_join_continuum_across_mrs` (reused unchanged on Mu2MidZM).
- **`continuum.py`**: `detect_continuum_simple` → delegates to
  `Mu2MidZM.has_continuum` (whole-segment same-modulus, no false negatives, §1);
  the old two-point+dV gate removed.
- **`sgbz_solver.py`**: `_evaluate_winding` builds μ₂_mid once per probe and
  reads `has_continuum` inline (§6.3); `handle_continuum` materializes the
  LineSubsets from the built ZM (no second build); `collect_GBZ_subsets`
  returns them with `is_continuum=True`.
- **`__init__.py`**: exports `Mu2MidZM`, `ItemView`, `Mu2MidBreakpoint`,
  `extract_continuum_linesubsets`, `WindingFun`, `get_winding_number`; API
  otherwise stable.

## Verified

- `pytest tests/test_sgbz.py::TestSGBZ10::test_outside_spectrum`,
  `TestCrossingDetection::test_two_points_one_crossing`,
  `TestCrossingDetection::test_charges_propagate_winding`, `test_all_exports`
  — 4/4 pass (2.67s).
- 0D: 2 PointSubsets, 1 ordinary charge ±1, correct winding sign.
- Exact continuum (μ₁=0.2 direct build): `has_continuum=True`, 2 LineSubsets
  at |β₂|=exp(γ₂)=1.35 (the GBZ boundary radius).

## Blocked — `test_inside_spectrum_is_continuum` FAILS

The bisection can't land on the continuum (μ₁=γ₁=0.2 exact) because the
ZeroManager root solver is **unstable near the degenerate continuum**. The
boundary-pair modulus gap is non-physically discontinuous in eps=μ₁−0.2:

| eps | gap_max |
|-----|---------|
| 0 | ~2.5e-13 (continuum, has_continuum=True) |
| 1e-10 | ~3.8e-7 |
| 1e-9 | ~3.9e-6 |
| 1.4e-8 | ~1.924  ← 5e5× jump over ~14× in eps (non-physical) |
| ≥1.4e-8 | saturates ~1.924 |

So `has_continuum` triggers only at machine precision (eps~1e-12); the
bisection's near-boundary midpoints (eps~1e-8) have gap~1.924 → not clustered →
`has_continuum=False` → per-column detection finds spurious crossings (3, not
1) → garbage W (−1 / 0.816 / 0.667, non-monotonic in eps) → bisection converges
by xtol to μ₁≈0.200000014 (eps=1.4e-8, bracket does NOT contain 0.2) with 6
spurious 0D subsets, `is_continuum=False`.

Root cause = root solver / continuation instability near the degenerate
continuum (the `continuation` module), not the SGBZ μ₂_mid logic — the gap
should be continuous & monotonic (sqrt-branch) in eps but isn't.

The old working-tree two-point gate (gap<1e-6 + dV-match over ≥2 rows) likely
"passed" by catching the tangent touch at θ₁=0 (gap~0 there) as continuum —
arguably a false positive that produced the right `is_continuum` marker.

**Decision deferred to user** (see plan & memory
`sgbz-mu2mid-refactor-continuum-bisection-issue`):
- (A) investigate ZeroManager continuation near the degenerate continuum;
- (B) add W-jump detection in the bisection + refine to the exact continuum;
- (C) accept 0D-near-boundary, update the test expectation.

## Known performance issue (not addressed)

Each `collect_GBZ_subsets` continuum-boundary call is ~230s — `_resolve_continuum_winding`
does up to 8 Mu2MidZM rebuilds × bisection iterations. Separate performance concern.
