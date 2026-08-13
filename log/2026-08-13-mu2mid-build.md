# 2026-08-13 — Build piecewise-smooth μ₂_mid on ZeroManager

Implemented `Mu2MidZM(ZeroManager)` in `demos/demo_zm_gbz.py` (playground;
**no production code touched**). Builds the piecewise-smooth
`μ₂_mid(θ₁) = (ln|β_{j_lo}| + ln|β_{j_hi}|)/2` (`j_lo = sort pos M-1`,
`j_hi = sort pos M`) as a first-class data object with values, analytic
derivatives, and left/right derivatives at breakpoints. Replaces the old
Newton refinement with a cubic-Hermite bracketing iteration; generalises to
the continuum (full) case via a representative-item abstraction.

Driven by `log/SGBZ算法梳理.md` §2.4 (full version) and the user's design
that "abs_argsort only identifies; column indices do the precise solving".

## Design decisions (user, 2026-08-13)

1. **Subclass ZeroManager** (not a standalone builder) — `insert_solution`,
   `segments`, `abs_argsort`, `tangents`, `locate` reused as `self.*`; wall
   / refinement mesh insertions are the intended side effect (breakpoint
   rows get exact roots + finite tangent). `build_mu2_mid()` is called
   explicitly after `run()` (not overriding `run()`).
2. **Walls cover all three sort-adjacent pairs** `(M-2,M-1)/(M-1,M)/(M,M+1)`
   + MR rows — all three affect μ₂_mid smoothness: the first two change the
   `{j_lo,j_hi}` set (real derivative jump); the (M-1,M) pair is the PMGBZ
   boundary (derivative continuous but sort unstable).
3. **abs_argsort only identifies; column indices do the precise solving.**
   `abs_argsort`/`sort_to_item` identifies the changing pair and `pair_kind`;
   the cubic Hermite iteration, `gap`, `value`, `deriv_left/right` all use
   column indices (rep_col) + the θ* row tangent.
4. **Convergence = derivative-sign + f_pred-sign bracketing + width.** Not
   `|f|<tol` (unreliable near tangency), not pure bisection. The cubic
   Hermite root gives the local derivative sign (increasing/decreasing
   crossing); with the true `f_pred` sign it fixes which side of the root
   θ_pred is on, so the bracket always still brackets the zero. Width
   `<xtol` ⟹ converged (derivative sign guarantees the root is inside).
5. **MR breakpoints taken from `multiple_roots` directly** — MR rows have
   cluster snapping (exact equal modulus), `abs_argsort` fully unstable;
   `_find_breakpoints` skips MR rows, `_find_mr_breakpoints` adds them from
   `multiple_roots.cluster_indices` (cluster covering M-1 or M).
6. **Continuum full version = representative-item set.** Replace the full
   column set by cluster representatives + single roots, each with a
   multiplicity. "M and M-1 → cluster index re-sorted" = `sort_to_item`
   cumulates multiplicities; `j_lo/j_hi` are item indices; `j_lo==j_hi` ⟹
   continuum (1D subset, |β_M|=|β_{M+1}| holds identically).
7. **Continuum detection covers all three pairs, simplified version too.**
   `detect_continuum_full` only checks the (M-1,M) boundary pair; but
   (M-2,M-1)/(M,M+1) continua also destabilise `abs_argsort` (j_lo/j_hi
   jump inside a same-modulus cluster → spurious sort-change). The
   simplified version (`continuum_clusters=None`) runs an internal
   whole-segment same-modulus detection over **all** column pairs (transitive
   closure) and builds the ItemView with representatives — so the simplified
   and full versions are unified (both use representatives; differ only in
   detection source: internal cheap vs external refined).

## Algorithm (4 stages)

### Stage 0 — walls
`_find_near_tie_rows`: mark MR rows + rows where any of the three sort pairs
is near-equal (`|ln|β_p|−ln|β_{p+1}|| < tie_tol`, using item moduli; same-item
pair = continuum, skipped). `_build_walls`: insert at `θ_row ± wall_frac·grid`
(default 0.1) so sort-change detection runs only outside walls; wall points
provide clean (finite-tangent) endpoints.

### Stage 1 — sort-change detection (outside walls)
`_detect_sort_changes`: scan adjacent rows where `j_lo`/`j_hi` (item index)
changes. For each, the changing item pair `(a,b)`; `pair_kind` from
`boundary = min(e_a, e_b)` (last sort position of each item in the left row,
counting multiplicity; `=M-2/M-1/M` ⟹ `M-2_M-1/M-1_M/M_M+1`; mult=1
degenerates to `min(inv_l[a], inv_l[b])`). Outputs `rep_col` pair for the
cubic Hermite step. Dedup by `(a,b,round(t_mid,8))`.

### Stage 2 — cubic-Hermite bracketing iteration (not Newton)
`_cubic_hermite_iterate`: bracket `[θ_lo,θ_hi]`, `f = ln|β_a| − ln|β_b|`.
Each step: build cubic Hermite from `(value, Re(V_a)−Re(V_b))` at both ends
→ `np.roots` predicts θ_pred → `insert_solution` takes true `f_pred` and
true derivative → derivative sign (increasing/decreasing) + `f_pred` sign
fixes which side θ_pred is on → tighten bracket → rebuild cubic. Converge
when width `< xtol`. **Why not Newton:** Newton uses `f/f'` at one point;
near a multiple root `∂f/∂β₂=0` ⟹ `f'` diverges ⟹ huge steps. The bracket
uses `f'` only at the (wall-outside, MR-free) endpoints, which stay finite;
bracketing is monotone. **Why not sign-change check:** sign-change fails when
`f_pred` has the same sign as both ends (cubic approximation error) — the
best_pred fallback could deadloop; derivative sign doesn't depend on the
endpoints and always fixes the side.

### Stage 3 — assemble
`_assemble`: flatten all segments (incl. wall/refinement inserts) into
`mu2_mid_{theta1,values,derivs,jlo,jhi}` (jlo/jhi are item indices; value/deriv
from `item_logabs`/`item_tang_re`). `_find_breakpoints`: item-index changes;
`theta1`/`value` taken at θ* (the smaller-gap row, the refinement point),
`gap = |ln|β_a|−ln|β_b||` at θ*, `deriv_left/right` computed at the θ* row
with left/right-segment item indices + `item_tang_re` (strict left/right
limits; swap ⟹ same set ⟹ strictly continuous). MR rows skipped (Stage MR
handles them). `_find_mr_breakpoints` appends MR breakpoints from
`multiple_roots`.

## ItemView (per segment)
```
rep_cols      (n_items,)  representative column (cluster[0] or single)
mults         (n_items,)  cluster size / 1
item_logabs   (N, n_items) = log|tracked_roots[:, rep_col]|
item_tang_re  (N, n_items) = Re(tangents[:, rep_col])
sort_to_item  (N, K)      sort pos p → item index (item by modulus, repeated mult)
j_lo/j_hi     (N,)        sort pos M-1/M item index
```
Original version = special case `n_items=K, mult=1, sort_to_item=abs_argsort`.

## Continuum detection (simplified, internal)
`_detect_continuum_clusters_internal(tie_tol)`: per segment, all column
pairs `(j,k)` with `max over rows |ln|β_j|−ln|β_k|| < tie_tol` are
same-modulus; transitive closure (BFS) → clusters. Whole-segment same-modulus
= continuum (real-analyticity); accidental (isolated) doesn't pass the
whole-segment max, auto-excluded. Covers all three sort pairs.

## Convergence / correctness (synthetic models)

**Non-continuum** (`build_synthetic_mu2_model`, β₂ odd powers + β₁-dependent
coeffs to break the 2-band self-inversive & even-power degeneracies;
E=2+0.3j, μ₁=−0.1):
- 8 breakpoints, 8/8 refined converged, all gaps ~1e-16 (machine precision).
- 4 PMGBZ swaps (M-1_M): `|Δd| = 0.00` strictly (derivative continuous —
  same `{a,b}` set ⟹ `(V_a+V_b)/2` identical).
- 4 internal (M-2_M-1/M_M+1): `|Δd|` 0.12–0.14 (rel 93–109%, real jump).

**Continuum** (`build_continuum_mu2_model`, β₂² quadratic ⟹ ±-pairs
same-modulus; clusters {(0,1),(2,3)}):
- Internal detection finds `[[(0,1),(2,3)]]`.
- Simplified version (None → internal): **0 spurious breakpoints**
  (vs 126 without representatives), `jlo` adjacent changes 0 (vs 80).
- Matches the full version (external clusters) exactly: mu2_mid range
  `[-9.6150e-02, -4.2806e-02]` both.

## Known limitations / TODO
- **MR left/right derivative** — `_find_mr_breakpoints` sets
  `deriv_left = deriv_right` (cluster snapping leaves one tangent set; the
  ± branch sign of the square-root branch point is not distinguished).
  Marked as a large finite number (tangent → ∞ near the multiple root);
  downstream interpolation must special-case MR.
- **`detect_continuum_full` (production) only checks (M-1,M)** — the demo's
  internal detection covers all three pairs but production `brute_force_SGBZ
  /continuum.py` does not. Extending production to three pairs is a future
  step (the mu2_mid builder itself is ready to consume either).
- **Synthetic models only** — 2-band Haldane has μ₂_mid≡0 (roots pair as
  (β,1/β) ⟹ |β_M·β_{M+1}|=1). A real 3-band / non-self-inversive physical
  model would be a better end-to-end test.
- **`j_lo==j_hi` continuum (1D subset)** — not exercised by the synthetic
  test (boundary pair straddles two clusters there). The algorithm supports
  it (`_find_near_tie_rows` skips same-item pairs; `_assemble` uses item
  indices) but no test hits M-1,M inside one item yet.

## Files
- `demos/demo_zm_gbz.py` only: `Mu2MidZM`, `ItemView`, `Mu2MidBreakpoint`,
  `build_synthetic_mu2_model`, `build_continuum_mu2_model`,
  `test_mu2_mid`, `test_mu2_mid_continuum`. `__main__` runs both.
- Outputs: `demos/mu2_mid_demo.png`, `demos/mu2_mid_continuum_demo.png`.
- Plan: `.claude/plans/build-mu2mid.md`.
