# 2026-08-14 — SGBZ review fixes: zero-mutation crossing phase, clamp policy, full TODO cleanup

Code review of `brute_force_SGBZ/` found 17 findings (recorded in
`TODO/sgbz-review-2026-08-14.md`, now resolved). Two of them were design-level
and are fixed per user decisions; the rest are correctness/dead-code/doc
cleanups. Baseline before the round: 17 passed (~127 s). After: 116 passed.

## 1. Fundamental fix #1 — the crossing phase no longer mutates ZeroManager

**Bug.** `detect_crossings_simple` scanned per-track sign changes (`sc`) and
exact touches (`touch`) on a snapshot of `g = ln|β_j| − μ₂_mid`. The `sc` loop
called `_bracket_crossing`, which inserted one mesh row per probe via
`insert_solution`; the `touch` loop then read `tracked_roots[i, j]` with the
pre-mutation row index — insertions shifted indices, so a touch PointSubset
could pair `beta1` with a β₂ from a different θ₁ row (silent wrong points).

**User decision (architecture, binding):** the mesh may be mutated ONLY during
the μ₂_mid build phase; the crossing-refinement phase after the build must not
insert anything. Implemented as:

- **`continuation/zero_manager.py`** — new `ZeroManager.solve_at(theta1,
  seg_idx, i) -> (roots, V)`: "insert_solution without the insert". Same
  solve + `interpolate_roots` prediction anchor + Hungarian track matching +
  `compute_tangent`; nothing is written into `SegmentData`.
  `insert_solution` is refactored to solve_at + the mesh insert (bit-for-bit
  unchanged behavior; coincidence semantics stay in `insert_solution`).
- **`crossings.py`** — `_bracket_crossing` rewritten with ZERO mutation. New
  `_eval_g`: exact-equality mesh-row read (`_track_g_and_gp`) or a transient
  `solve_at` probe, cached per θ float (a predicted point becomes the next
  iteration's bracket endpoint; the cache key is the float as assigned, so
  dict equality is exact). The bracket lives inside ONE original adjacent
  mesh interval (located once at entry; tightening only shrinks it). Check
  order preserved (touch → width-converged → same-sign → cubic). On width
  convergence it returns the LAST predicted point — its β₂ is the true solved
  root, more accurate than the old mesh-row fallback.
- **Seam behavior preserved bit-for-bit**: the mesh-row hit test is exact
  float equality, NOT a tolerance — inside the 2π seam's `_BOUNDARY_THETA_TOL`
  wrap band `solve_at` raises ValueError (wrapped θ leaves the interval),
  which reproduces the old near-seam abandonment; a tolerance would have
  emitted the crossing at θ≈2π and duplicated the interval-0 detection.
- **`mu2mid.py`** — `_Mu2MidPath` moved here from winding.py (crossings needs
  it to evaluate the built μ₂_mid curve at arbitrary θ; winding imports
  crossings, so the move avoids a cycle). `Mu2MidZM.insert_solution`'s
  post-build refresh branch is now defensive-only (nothing inserts
  post-build); kept as a safety net, docstring updated.
- Tests: `test_solve_at_matches_insert_solution` (roots/V equal what a
  subsequent insert stores; solve_at alone changes no shape),
  `test_solve_at_hint_equals_search`, `test_solve_at_single_row_raises`,
  `test_detection_leaves_mesh_unchanged` (the binding contract, direct).

Side effect: detection no longer re-`_assemble`s per probe, so test_sgbz
dropped from ~127 s to ~66 s (later runs ~130 s under machine load — same
numerics for finite-root models).

## 2. Clamp policy — 0/∞ handling only where μ₂_mid is built

**User decision:** clamp the 0/∞ padding roots' ln|β| to ±14 ONLY when
building μ₂_mid (the boundary-pair mean: "treat the infinity as ±14 and
average with the other"); the rest of the pipeline reads the raw root values.
The roots themselves were never modified — the clamp only ever affected
derived ln|β| views.

- **Clamp points** (μ₂_mid construction only): `mu2mid._assemble` (per-term
  clamp then mean), `_find_breakpoints`'s `value_bp`, `winding._build_root_mu2_mesh`
  (later deleted with #3); `_Mu2MidPath` keeps its band clamp on the
  interpolated value (cubic overshoot guard).
- **Raw everywhere else**: ItemView clustering, near-tie / sort-change
  comparisons, `_cubic_hermite_iterate`, `crossings._track_g_and_gp` /
  `_eval_g` / `detect_crossings_simple`. This is deliberate: a padding
  track's `g = ±∞` never changes sign (contributes no crossing); two
  same-type padding roots give a NaN modulus difference and fail the
  same-modulus cluster test — the old pseudo-continuum risk (#4) disappears
  by construction.
- **M = 0 / N = 0 unsolvable**: when the SGBZ modulus must be 0 or ∞ the code
  cannot solve the model; a prominent warning block added to `doc/SGBZ.md`
  §2.1 telling callers to exclude such models (`CharPoly.get_minor_degrees()`).

## 3. Remaining fixes (#2–#17)

- **#2** `solve_SGBZ_for_E` Step 1.2 discarded `is_boundary` when the right
  bracket endpoint was itself a continuum boundary → bisection returned a
  spurious ±1-winding GBZ at the wrong μ₁. Now returns the continuum result.
- **#3** `_pick_seed_theta2`'s silent fallback to a crossing endpoint: the
  `_loop_min_dist` pre-screen (already removed in a prior edit) and its dead
  machinery (`_loop_min_dist`, `_build_root_mu2_mesh`, dead vars) deleted.
- **#4** pseudo-continuum from clamped padding roots — resolved by §2.
- **#5/#6** dead parameters removed end-to-end: `detect_threshold`,
  `dedup_tol`, `dV_tol`, `vote_frac` (crossings → winding → plateau →
  sgbz_solver), plus dead constants `_DETECT_THRESHOLD`, `_DEDUP_TOL`,
  `_TANGENT_F_TOL`, `_TANGENCY_THRESHOLD`, `_DV_TOL`. The inline continuum
  gate now uses the named `_INLINE_CONTINUUM_FRAC = 0.1` (kept low by design:
  whole-segment clustering already guarantees genuine same-modulus, the frac
  only guards against a spurious row — a missed continuum would break the
  bisection); `CONTINUUM_FRAC = 0.9` documented as a legacy API-compat export.
- **#7/#8** dead `tol` param and `touches_mr` (gone with #1); dead `M` params
  removed from `_loop_winding_quad` and `compute_average_winding` (public
  signature updated).
- **#9** inline clamp duplication in `_track_g_and_gp` — moot after §2.
- **#10/#11** stale docstrings fixed: `detect_crossings_simple` Returns
  ("one PointSubset per detected zero-curve"), `collect_GBZ_subsets` and
  `GBZResult.is_continuum` (materialization no longer TODO).
- **#12** plateau probe built 3 ZMs per probe point (plain ZeroManager +
  one inside `detect_continuum_simple` + another inside
  `detect_crossings_and_winding`); now ONE `Mu2MidZM` build with the inline
  `has_continuum` gate, mirroring `sgbz_solver._evaluate_winding`. Note: the
  probe's clustering tie_tol moves from the old default 1e-8 to
  `continuum_tol` = 1e-6, consistent with the solver.
- **#13** bare `except Exception` in the probe → `warnings.warn`
  (RuntimeWarning); failures surface instead of silently reading as
  non-plateau (CLAUDE.md: unexpected results must be reported).
- **#14** `compute_average_winding` now returns `float` in both branches.
- **#15** `_Mu2MidPath` breakpoint override now applies to ALL rows matching
  the breakpoint θ (segment-boundary MR rows appear twice in the flat array;
  argmin alone silently dropped `deriv_right` past the MR).
- **#16** bisection updates `w_low`/`w_high` with the bracket, so
  `_winding_bracket` diagnostics stay honest (f_mid may be a continuum
  proxy; its sign is all the bisection uses).
- **#17** `get_winding_number` epsabs/epsrel = 1e-3 documented: total error
  ≲ #intervals·1e-3 → ~1.6e-3 winding units, two orders below the 0.5
  rounding margin.

## 4. Docs and API

- `doc/SGBZ.md` synced: §2.3 no-dedup / per-zero semantics, §2.4 (pre-screen
  removed), §3.2–§3.8 signatures, §4 parameter tables, the M=0/N=0 warning.
- `demos/demo_zm_gbz.py` import block updated (`_DV_TOL` removed; constants
  imported from `brute_force_SGBZ.mu2mid`).
- **Breaking API notes**: `collect_GBZ_subsets` no longer accepts
  `dV_tol` / `vote_frac` / `detect_threshold` / `dedup_tol` (passing them
  triggers the unrecognized-options warning); `compute_average_winding`
  loses its `M` parameter.

## Verified

- `python -m pytest tests/ -q` — **116 passed** (was 116 before the round's
  final pass; new tests: 3 × solve_at + 1 × detection-leaves-mesh-unchanged).
- `test_two_points_one_crossing`'s `beta1` rel-1e-9 assertion still holds
  exactly (both boundary tracks' brackets touch the same refined mesh row).

## Known limitation (unchanged)

`_resolve_continuum_winding` perturbation scales cap at 8 × perturb (≈0.08);
a wider continuum band raises ValueError → `success=False`. This is the
2026-08-13 "Blocked" item — ZeroManager root instability near the degenerate
continuum, still awaiting a decision.
