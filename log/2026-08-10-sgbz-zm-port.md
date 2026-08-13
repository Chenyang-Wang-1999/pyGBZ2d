# 2026-08-10 — Port SGBZ to continuation.ZeroManager

Replaced the dead/misfit old `brute_force_SGBZ/` modules (uniform-mesh root
solving + scipy-quad winding + a `LineSubset(theta1_start=...)` constructor
call that broke after the `LineSubset` refactor to `(theta1_arr, beta2_arr)`)
with a `continuation.ZeroManager`-based implementation faithful to
`log/SGBZ算法梳理.md` §1–3, ported from the reviewed `demos/demo_zm_gbz.py`.

## Decisions (user, 2026-08-10)
1. **Full replace** — old `pmgbz_detector.py`, `strip_winding_number.py`,
   `root_solver.py`, old `winding.py`, old `SGBZ.py` deleted; tests rewritten.
2. **LineSubset ignored, marked TODO** — when a continuum is detected at the
   solved boundary μ₁, return `GBZResult(is_continuum=True, subsets=[], index=(0,0))`.
   The amoeba `_LinePiece`/`_join_continuum_across_mrs` extractor was NOT ported.
   Spectrum membership still works (decided by the winding-zero / left-right
   bisection, independent of subset materialization).
3. **Zero-plateau check kept**, default ON.

## New module layout
```
brute_force_SGBZ/
├── __init__.py        # whitelist exports (re-exports gbz_types classes)
├── continuum.py       # §1: detect_continuum_simple / detect_continuum_full
├── crossings.py       # §2: detect_crossings_simple (+ charge classification)
├── winding.py         # §3: compute_average_winding, detect_crossings_and_winding
├── plateau.py         # _check_pmgbz_points_clustered, _probe_zero_plateau_near_mu1
└── sgbz_solver.py     # solve_SGBZ_for_E (μ₁ bisection), collect_GBZ_subsets
```

## gbz_types change (backward-compatible)
- `GBZResult.is_continuum: bool = False` added — distinguishes the
  continuum-TODO state (in spectrum, no subsets) from "outside spectrum".
- `is_gbz` now `success and (index != (0,0) or is_continuum)`.
- Amoeba never sets it (it has the extractor) → stays `False`.

## Behavior (2D HN model, γ₁=0.2, γ₂=0.3)
W(μ₁) curve is smooth/monotonic, flat plateau ±1 far from γ₁, crosses zero
at μ₁=γ₁ as a **continuum** (the model's SGBZ boundary is genuinely a 1D
LineSubset, not isolated points — the old tests asserting PointSubsets there
relied on the broken fallback). 0D PointSubsets appear *off* the boundary
(e.g. μ₁=0.1 → 2 points, one crossing, ordinary ±1 charge); exercised by
the new `TestCrossingDetection`.

## Tests
`tests/test_sgbz.py` rewritten: 15 tests, all passing.
- Inside spectrum → `is_gbz and is_continuum`.
- Outside spectrum (E=5) → `not is_gbz and not is_continuum`.
- `solve_SGBZ_for_E` locates μ₁=γ₁ ([10]) / γ₁+γ₂ ([11]) within 2e-3.
- Hermitian limit → continuum on unit circle (cross-checked vs amoeba).
- `TestCrossingDetection` exercises the §2 cubic+Newton+charge path at a
  non-continuum μ₁.

## Pre-existing failures (NOT from this port)
`tests/test_gbz_types.py` has 6 failures (`TestLineSubset::*`,
`TestGBZResult::test_nonempty_line`, `test_mixed`) — they use the obsolete
`LineSubset(theta1_start=..., beta1=..., beta2_mat=...)` constructor from
before the `LineSubset` refactor. Verified pre-existing via `git stash`.

## TODO (future)
- §2.4 full crossing detector (`detect_continuum_full` consumer) — crossing
  detection in the presence of a continuum, for LineSubset endpoint solving.
- LineSubset extraction (port amoeba's `_join_continuum_across_mrs`).
- `_loop_winding_segments` evaluates `poly.eval_val` per-row in a Python loop
  (scalar-only API); vectorization is a future optimization.
- Fix the 6 pre-existing `test_gbz_types.py` LineSubset constructor failures.
