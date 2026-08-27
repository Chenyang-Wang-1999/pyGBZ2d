# pygbz2d.sgbz — SGBZ Spectrum Calculation

Non-Hermitian spectrum computation based on the SGBZ (Strip Generalized Brillouin Zone) formulation. The 2026-08-13 μ₂_mid framework rewrite unifies SGBZ crossing detection with amoeba's per-column-vs-μ₂ structure — the only difference is that μ₂ is the piecewise-smooth boundary-pair mean curve `μ₂_mid(θ₁)` instead of a constant. See `log/2026-08-13-SGBZ算法梳理.md` for the authoritative algorithm writeup.

## 1. Theoretical Background

### 1.1 SGBZ and major-axis winding number

SGBZ integrates along the `θ₁` direction, with the base manifold `μ₂ = μ₂_mid(θ₁; E, μ₁)` — the piecewise-smooth mean of the two boundary-pair moduli `ln|β_{j_lo}|` and `ln|β_{j_hi}|`, where `j_lo`/`j_hi` are the item indices occupying sorted positions `M-1`/`M` (0-based). The average major-axis winding number `W(E, μ₁)` is the winding number around the `θ₁` direction with constant `θ₂`; its zero in `μ₁` defines the SGBZ boundary.

### 1.2 PMGBZ Point Classification

On the `θ₁` axis the `β₂` roots are sorted by modulus. When `|β_M| = |β_{M+1}|` a degeneracy occurs — this is a PMGBZ point:

- **Continuum (1D)**: `|β_M| = |β_{M+1}|` holds identically over a continuous `θ₁` interval. The boundary pair collapses to one continuum item (`j_lo == j_hi`); `W` is undefined there.
- **Accidental (0D)**: `|β_M| = |β_{M+1}|` only at isolated `θ₁` positions. Each accidental point is a transversal crossing of one zero-curve (track) through `μ₂_mid`.

### 1.3 Spectral Inclusion Relation

```
σ_Amoeba ⊇ ∪_j σ_SGBZ,j
```

In general the SGBZ spectrum is a subset of the amoeba spectrum. For uniform bands (geometry-dependent spectral degeneracy) the two are equal.

## 2. Key Components

### 2.1 `mu2mid.py` — μ₂_mid as a First-Class Object

`μ₂_mid(θ₁) = (ln|β_{j_lo}| + ln|β_{j_hi}|) / 2`, piecewise smooth; breakpoints are the `θ₁` where `j_lo`/`j_hi` change, plus the multiple-root (MR) rows where the sort order is undefined. `Mu2MidZM(ZeroManager)` builds this representation as a first-class object with analytic values and derivatives.

**`ItemView`** — a per-segment representative-item view that collapses continuum (whole-segment same-modulus) column clusters into one item with a multiplicity, so `sort_to_item` is stable where `abs_argsort` is not (cluster snapping makes `abs_argsort` jump arbitrarily between equal-modulus columns). `j_lo = sort_to_item[:, M-1]`, `j_hi = sort_to_item[:, M]`; `j_lo == j_hi` iff the M-1/M boundary pair is one continuum item (1D subset). The no-continuum case is the special case `n_items = K, mult = 1, sort_to_item == abs_argsort`.

**`Mu2MidBreakpoint`** — one breakpoint of the curve. `pair_kind` is one of `M-2_M-1` / `M-1_M` (PMGBZ boundary crossing, `is_pmgbz=True`) / `M_M+1` / `multi` (MR). `value` is continuous across a breakpoint; `deriv_left`/`deriv_right` may jump (a swap with the same {a,b} set stays strictly continuous; an internal break changes the set and gives a real jump). MR breakpoints have `deriv_left = deriv_right = inf` (the `dβ/dθ` branch-point divergence).

**`analyze` / pairwise crossing** — the 2026-08-25 build:

| Stage | Action |
|-------|--------|
| 1 — ItemView | Whole-segment same-modulus clustering (vote + BFS) → representative items. Cluster semantics: **columns equal-modulus over the whole segment**; isolated ties are never injected into ItemView. |
| 1b — multi-crossing refinement | Before the pairwise scan, every interval/pair is screened with the cubic-Hermite interpolant of `d = item_logabs[a] − item_logabs[b]` (endpoint values + tangents). Intervals whose same-sign cubic dips near zero, or that already touch / change sign, have their predicted interior roots collected; two or more well-separated roots trigger a sub-mesh (`ρ=4` sub-intervals in the narrowest gap + separation midpoints, capped at 64, up to 3 rounds). This lets the sign-change scan below see even numbers of crossings per interval. |
| 2 — pairwise intersections | For every representative item pair, `d = item_logabs[a] − item_logabs[b]`. `d[:-1]==0` is an exact touch (not refined, left-closed/right-open); `d[:-1]*d[1:]<0` is a transversal crossing → linear prediction + `brentq` refinement. Zero θ deduplication. |
| 3 — EventGroup | θ* closer than `crossing_tol` (including the θ=0≡2π circular seam) merge into one mesh row at the midpoint; connectivity/charges are kept per member event. Every downstream stage consumes EventGroups only — no singleton special case. A merged group's representative row AND every original touch row are event rows. |
| 4 — mesh insert | Each EventGroup inserted once with `insert_solution(interp='linear')`. Adjacent event positions inside a segment are separated by a regular row at their midpoint, directly solved with `insert_solution(interp='hermite')`; rows re-resolved on the final mesh; ItemView is rebuilt. The invariant "no two event rows adjacent inside one segment" is checked explicitly. |
| 5 — Mu2Mid build | Independent `Mu2Mid` object (not a ZeroManager): piecewise-smooth intervals, `hermite_interp_poly` inside each piece, ordinary knots C1, event/MR/±14 knots are breakpoints. Each Hermite piece is split at its ±14 roots; outside-band subpieces become `v=±14, dv=0`. |

**Inline continuum detection**: `has_continuum` is read directly off the ItemView — any row with `j_lo == j_hi` means one continuum item occupies both M-1 and M there, so the boundary pair is a 1D subset; no second scan or fraction gate is needed. So the bisection gate *is* the build itself (§1/§6.3), with no separate two-point detector and its false-negative risk. The `_detect_continuum_clusters_internal` same-modulus criterion is a **column-level vote**: column pair `(j,k)` is clustered when `|ln|β_j| − ln|β_k|| < tie_tol` holds on a fraction of rows strictly above `CONTINUUM_FRAC` (0.9, matching amoeba's `_continuum_mask`), followed by BFS transitive closure. Real-analyticity guarantees same-modulus along the whole zero curve (same column index, until an MR), so a majority vote over sampled rows is the correct detector; the old whole-segment max criterion was over-conservative and could split a genuine cluster on a single noisy row.

**Mesh-mutation contract** (2026-08-25): the mesh is mutated in two batch phases — first by the multi-crossing refinement (stage 1b, `insert_solution(interp='hermite')`), then by the EventGroup insertions (`insert_solution(interp='linear')`, plus directly-solved regular separator rows between adjacent events). The crossing refinement itself solves transiently via `ZeroManager.solve_at` and never inserts; group rows are resolved on the final mesh so later insertions cannot leave stale indices.

**Unified Hermite interpolation** (2026-08-15): cubic-Hermite polynomial construction goes through `continuation.interpolation` — `cubic_hermite_poly` / `hermite_interp_poly` (automatic linear fallback for non-finite endpoint derivatives). The new `Mu2Mid` path uses it per smooth piece; crossing refinement itself is purely linear + `brentq`.

**0/∞ clamp policy**: ItemView sorting and pairwise `d` scanning read the **raw** `ln|β₂|`, so 0/∞ padding roots stay `∓∞` at the outer sort positions. They must not become the M-1/M boundary pair: `_build_item_views` raises explicitly in that case, because the SGBZ modulus would be exactly 0 or ∞. The clamp to `±_LOGABS_CLAMP_L = ±14` applies only to finite boundary-pair knot values while building μ₂_mid, not to padding roots.

> **⚠ Warning — M = 0 or N = 0 is NOT solvable.** When the characteristic polynomial has `M = 0` or `N = 0` (check `CharPoly.get_minor_degrees()`), the SGBZ boundary pair necessarily includes padding roots at `β₂ = 0` or `β₂ = ∞`, i.e. **the SGBZ modulus itself must be 0 or ∞**. This code cannot solve such models: the clamped μ₂_mid (band edge ±14) is only a finite stand-in, not the true boundary, and the crossing/winding on a degenerate 0/∞ boundary is unreliable. Callers must exclude `M = 0` / `N = 0` models before invoking `collect_GBZ_subsets`, or treat any result as unphysical.

### 2.2 `continuum_lines.py` — Continuum Detection and LineSubset Materialization

Continuum detection is folded into the μ₂_mid build (`Mu2MidZM.has_continuum`). This module provides two entry points:

**`detect_continuum_simple`** — a presence-only flag for callers that hold a plain `ZeroManager` and only need "is there a continuum?" (the plateau probe). It builds a fresh `Mu2MidZM` and returns `has_continuum`. The μ₁ bisection itself does NOT call this — it constructs its own `Mu2MidZM` and reads `has_continuum` inline.

**`extract_continuum_linesubsets`** — LIVE materialization of the 1D continuum tracks as `LineSubset`s (not TODO). SGBZ LineSubset semantics differ from amoeba's: an SGBZ LineSubset is the stretch where a specific boundary pair (sorted positions M-1/M) is the degenerate continuum item, terminating at EITHER a modulus-sort change (another root overtakes the boundary pair — a hard terminator *inside* a segment) OR a multiple root (MR) at a segment boundary. So a piece is a **contiguous run** of `j_lo == j_hi == item` rows, NOT the whole segment.

- `_find_boundary_runs`: split each segment's `j_lo == j_hi == item` rows into maximal contiguous runs (multiple runs per segment are possible — the same continuum pair can drop out of the boundary and re-enter later).
- `_runs_to_pieces`: one `_LinePiece` per run per continuum track (each cluster column is a distinct β₂ curve at the same `|β₂|`).
- `_join_runs_across_mrs`: join pieces whose endpoints touch a segment edge (MR / θ₁=0≡2π seam). An endpoint strictly inside a segment is a sort-change terminator — no join. At an MR, if the endpoint track column is in the MR's cluster the track terminates; otherwise it continues into the adjacent segment by column identity (`boundary_perm` at the cyclic seam). Iterated to a fixpoint so chains and the cyclic seam converge.

### 2.3 PointSubset Materialization (in `winding.py`; formerly `crossings.py`)

Crossing DETECTION lives in `pairwise.py` and is orchestrated by
`Mu2MidZM.analyze`. `winding.detect_crossings_simple` only materializes
results (the former `crossings.py` module was dissolved 2026-08-18: its
materialization half merged into `winding.py`, the `Mu2MidZM` bootstrap
`ensure_mu2mid` moved to `mu2mid.py`).

**Detection**: on each segment's ItemView, every representative item pair is
scanned by `d = item_logabs[a] − item_logabs[b]`:
- `d[:-1] == 0` → exact touch, left-closed/right-open, no refinement;
- `d[:-1] * d[1:] < 0` → transversal crossing: linear prediction → `brentq`
  on the true `ln|β_a| − ln|β_b|`.
- No θ deduplication; θ* closer than `crossing_tol` (including across the
  θ=0≡2π circular seam) merge into one `EventGroup` mesh row.
- Event rows = EventGroup representative rows + original touch rows.  Between
  adjacent event rows inside one segment a regular midpoint row is inserted by
  a real polynomial solve (`interp='hermite'`), so `finalize_event_groups`
  never reads another event row as a "regular" side.

**PointSubset rule**: for each EventGroup, find the item/column connected
component whose sort positions cover BOTH `M-1` and `M`; expand it to all
REAL columns (never just representatives) and emit one `PointSubset` per
column. Charge is the side-change rule
`q = (side_right − side_left)/2` where `side = +1` for positions ≥ M and
`-1` for positions < M; q may be 0 for merged events.

**Charge classification** (in the materialized charge dicts):
- Classification is stored per REAL column in `EventGroup.column_kind`, so
  one merged θ can simultaneously contain an MR component and an unrelated
  ordinary component.
- **Ordinary** (charge ±1/0): the side-change charge `q` computed by
  `pairwise.finalize_event_groups`.
- **Hard boundaries — MR / tangent**: the charge is **unknown**. The charge
  dict stores `charge=None` (never a numeric placeholder), so any accidental
  arithmetic on a hard boundary's charge fails loudly with `TypeError`.

MR/tangent/unknown act as **hard** region boundaries (delimit regions, do not propagate winding); ordinary is **soft** (propagates via charge).

### 2.4 `winding.py` — Average Winding Computation

Computes the average major-axis winding number `W(E_ref, μ₁)` from a built `Mu2MidZM` at fixed `(E_ref, μ₁)` together with the charge list from `detect_crossings_simple` (same module, §2.3).

**Loop path** (§6.4): `β₂ = exp(μ₂_mid(θ₁) + iθ₂)` with `θ₁ ∈ [0, 2π)` and fixed `θ₂`. The μ₂_mid path is the independent piecewise-smooth `Mu2Mid` object: each smooth piece is `hermite_interp_poly` from endpoint `(value, derivative)`; ordinary knots share one derivative (C1), while event / MR / ±14-saturation knots are breakpoints with left/right derivatives. Values are bounded to `±14` at build time, and the winding quad is split at every piece boundary.

**Integration**: `WindingFun` computes `Im[f'(t)/f(t)]` along the loop; `get_winding_number` integrates via `scipy.integrate.quad`, **split at the μ₂_mid breakpoints** so every quad segment lies on one smooth piece. The analytic `dβ₂/dθ₁ = β₂ · μ₂_mid'(θ₁)` comes from the path's cubic-Hermite derivative.

**Region/seed/charge propagation** (§3.2): crossings partition `θ₂` into regions (hard boundaries) and intervals (soft boundaries). One seed interval per region, picked by `_pick_seed_theta2` to maximize `min |f(E, β₁(θ₁), β₂_loop(θ₁))|` over `θ₁` (the true safety metric — farthest from char-poly zeros, accounting for β₁ and `|∂f/∂β₂|` that a β₂-distance proxy ignores). Compute `w₀` there via `_loop_winding_quad`, propagate across the region's soft boundaries via charges (ordinary: ±1). Each region contributes exactly one loop-winding evaluation; ±1 numerical noise is confined to the per-region seed.

**Charge conservation**: when every boundary is soft (ordinary, charge ±1), the charges must sum to zero — one full `θ₂` circle must return the winding to itself. A non-zero sum means the crossing detector missed or duplicated a zero and raises `RuntimeError` instead of silently averaging incompatible per-interval windings. Any hard boundary (charge `None`, unknown) disables this check.

### 2.5 `sgbz_solver.py` — μ₁ Bisection Solver

Top-level solver using bracket expansion + plain midpoint bisection with continuum interception.

- **`_evaluate_winding`**: build a fresh `Mu2MidZM` at the probed `μ₁` and run `analyze` — the analysis both detects the continuum inline (`has_continuum`) and provides the `Mu2Mid` path the winding integral needs. When a continuum is detected `W` is undefined → returns `(None, None, zm)`; the caller resolves it via left/right limits and, if it is the boundary, materialises the LineSubsets from the built `zm`. Otherwise runs crossing detection + winding on the same built `zm` (no second analysis).
- **Bracket expansion**: expand left/right until winding signs straddle zero. Each side is capped at `_MAX_BRACKET_EXPANSIONS = 10` steps (aligned with amoeba's `max_range_expansions`); exhausting the cap raises `RuntimeError` instead of looping forever on unsolvable / anomalous winding. The right endpoint is established by one unified loop for both entry paths, and a continuum proxy correction is re-evaluated through the same sign check as a plain winding value (no non-straddling bracket can slip into bisection).
- **Bisection**: plain midpoint (not false-position) because the winding has flat plateaus (±1) with a narrow transition zone; false position stalls on this shape while midpoint guarantees bracket halving.
- **Continuum interception** (`handle_continuum`): when `w_mid` is `None` (continuum), resolve left/right limits via `_resolve_continuum_winding` (perturbation scales `(1, 2, 4, 8) × continuum_perturb`); if they straddle zero, that `μ₁` IS the SGBZ boundary and the 1D LineSubsets are materialised by `extract_continuum_linesubsets` from the built `zm`.
- **Convergence check**: only applied when `w_mid` is a real winding number. When `w_mid` is `None` (continuum proxy), the algorithm continues iteration without convergence testing — converging on a proxy would confuse "left limit is zero" with "W is zero".

### 2.6 `plateau.py` — Zero-Plateau Detection

Two-stage check, unaffected by the μ₂_mid framework:

1. **Clustering pre-check** (`_check_pmgbz_points_clustered`): thin SGBZ adapter over `pygbz2d.core.check_points_clustered_on_torus`. Euclidean distance on the `(θ₁, θ₂)` torus — **must consider both θ₁ and θ₂**, not θ₁ alone: the same `θ₁` can host multiple distinct `β₂` (degenerate pairs), so a θ₁-only check would misclassify. Example: gain-loss Haldane at E=0.5 has 6 points in 3 degenerate pairs; torus minimum distance ≈ 33°, correctly judged non-clustered.
2. **Probe ladder** (`_probe_zero_plateau_near_mu1`): thin adapter over `pygbz2d.core.probe_zero_plateau`. Probe `μ₁ ± step`; a plateau shows zero winding with empty GBZ on both sides.

### 2.7 Continuum Handling Summary

When `|β_M| = |β_{M+1}|` holds identically over a continuous `θ₁` interval:
- `Mu2MidZM.has_continuum` is `True` → winding undefined at that `μ₁`.
- `_resolve_continuum_winding` computes left/right limits via perturbation.
- If limits straddle zero: that `μ₁` is the SGBZ boundary (1D LineSubset case); `extract_continuum_linesubsets` materializes the LineSubsets from the built `zm`. `collect_GBZ_subsets` returns `GBZResult(is_continuum=True, subsets=[LineSubset, ...], index=(0, n_1d))`.
- If limits do not straddle zero: the bisection continues using the left limit as a proxy.

## 3. API Reference

### 3.1 `Mu2MidZM`, `analyze` and `Mu2Mid`

```python
class Mu2MidZM(ZeroManager):
    def analyze(
        self,
        continuum_clusters: list | None = None,
        *,
        tie_tol: float = CONTINUUM_TOL,
        crossing_tol: float = 1e-10,
        min_direction_deriv: float = 1e-12,
        verbose: bool = False,
        refine_multi_crossings: bool = True,
        refine_max_rounds: int = 3,
        refine_safety_factor: float = 4.0,
        refine_max_subintervals: int = 64,
        refine_max_total_inserts: int = 2000,
    ) -> None: ...
```

`ZeroManager` + ItemView analysis + pairwise crossings + μ₂_mid path.
Usage: `zm = Mu2MidZM(poly, E_ref, mu1); zm.run(); zm.analyze()`.
`build_mu2_mid(...)` is retained as a compatibility alias for `analyze`.

After `analyze`:
- `_pair_events` / `_event_groups`: pairwise events and the unified EventGroups;
- `has_continuum`: `any(j_lo == j_hi)`;
- `mu2_mid`: independent `Mu2Mid` object (`Mu2MidPiece` list, `value_deriv`, `breakpoints`);
- `mu2_mid_theta1/values/derivs` and `seg_mu2_*`: compatibility flat arrays.

### 3.1b `refine_mesh_for_multiple_crossings`

```python
def refine_mesh_for_multiple_crossings(
    zm,
    *,
    tie_tol: float = 1e-6,
    crossing_tol: float = 1e-10,
    max_rounds: int = 3,
    safety_factor: float = 4.0,
    max_subintervals: int = 64,
    max_total_inserts: int = 2000,
) -> int
```

Runs between `ZeroManager.run()` and `collect_pair_events` (called
automatically by `analyze` unless `refine_multi_crossings=False`).  It
screens every interval/pair with the cubic-Hermite interpolant of
`d = ln|β_a| − ln|β_b|`, predicts interior roots of suspicious pairs, and
inserts a sub-mesh for intervals with ≥ 2 well-separated predicted roots
(§2.1 stage 1b).  Returns the number of inserted mesh rows.

### 3.2 `detect_continuum_simple`

```python
def detect_continuum_simple(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    continuum_tol: float = 1e-6,
    zm_run_kwargs: dict | None = None,
) -> bool
```

Presence-only continuum detection for plain `ZeroManager` callers. Builds a fresh `Mu2MidZM` from `zm`'s `(poly, E_ref, mu1)` and returns `has_continuum`. The whole-segment same-modulus ItemView criterion replaces the old two-point gate — no false negatives at sub-grid continua. The μ₁ bisection does NOT call this; it reads `has_continuum` inline from its own `Mu2MidZM`.

### 3.3 `extract_continuum_linesubsets`

```python
def extract_continuum_linesubsets(
    zm: ZeroManager, poly: CharPoly, *, zm_run_kwargs: dict | None = None,
) -> list[LineSubset]
```

Materialise the 1D continuum LineSubsets of `zm`. **Precondition**: continuum detection has already run on `zm` and returned `has_continuum == True` — materialization is always a post-detection step. Raises `RuntimeError` when `m.has_continuum` is False (caller skipped the gate, or detection failed), and also when `has_continuum` is True but `_find_boundary_runs` finds no boundary run (the two continuum gates disagree — an invariant violation, never a valid empty result). Otherwise: one `LineSubset` per continuum track per boundary run, joined across MR boundaries where the track passes through as a non-cluster root. Each LineSubset's `θ₁` range is exactly where its track held the M-1/M boundary — sort-change terminators inside a segment end the piece.

### 3.4 `detect_crossings_simple`

```python
def detect_crossings_simple(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    crossing_tol: float = 1e-10,
    zm_run_kwargs: dict | None = None,
) -> tuple[list[PointSubset], list[dict]]
```

0D PMGBZ-boundary crossing detection + charge classification. Detection is done by `Mu2MidZM.analyze` (pairwise ItemView intersections); this function only materializes EventGroups. When a fresh `Mu2MidZM` is built, `_ensure_mu2mid` analyzes with `tie_tol=CONTINUUM_TOL` (1e-6). Assumes no boundary continuum — gate with `detect_continuum_simple` / `has_continuum` first.

`crossing_tol` is accepted for backward compatibility; the brentq refinement tolerance is fixed at `Mu2MidZM.analyze` time.

**Returns**:
- `subsets`: one `PointSubset` per real column of every M-1/M EventGroup component.
- `charges`: one dict per subset — `charge` from the side-change rule (`+1/-1/0`) or `None` for MR/tangent.

### 3.5 `detect_crossings_and_winding`

```python
def detect_crossings_and_winding(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    crossing_tol: float = 1e-10,
    zm_run_kwargs: dict | None = None,
) -> tuple[list[PointSubset], float]
```

Convenience function: crossing detection + average major-axis winding in one call. Builds μ₂_mid once (in `detect_crossings_simple` via `_ensure_mu2mid`); the winding reuses that built `Mu2MidZM`. Returns `(subsets, W_avg)`.

### 3.6 `compute_average_winding`

```python
def compute_average_winding(
    zm: ZeroManager,
    poly: CharPoly,
    charges: list[dict],
) -> float
```

Compute the average major-axis winding number `W(E_ref, μ₁)` from a built `Mu2MidZM` and the charge list from `detect_crossings_simple`. Partitions the `θ₂` circle into regions delimited by hard boundaries (MR/tangent/unknown, charge `None` — unknown), picks one seed interval per region (maximizing `min |f|` along the loop), computes `w₀` via `_loop_winding_quad`, and propagates across the region's soft boundaries (ordinary, charge ±1). When every boundary is soft the charges must sum to zero, otherwise `RuntimeError`; any hard boundary disables the conservation check. Arc-weighted mean of the per-interval windings.

### 3.7 `solve_SGBZ_for_E`

```python
def solve_SGBZ_for_E(
    poly: CharPoly,
    E_ref: complex,
    mu1_guess: tuple[float, float] = (-1, 1),
    zero_tol: float = 1e-10,
    continuum_perturb: float = 1e-2,
    max_iter: int = 60,
    zm_run_kwargs: Optional[dict] = None,
    *,
    continuum_tol: float = 1e-6,
    crossing_tol: float = 1e-10,
) -> dict
```

Locate the winding-zero `μ₁` and return solve diagnostics. Uses bracket expansion + plain bisection (midpoint) with continuum interception.

**Returns**: dict with keys:
- `"mu1"`: the solution `μ₁`.
- `"subsets"`: list of `PointSubset` (0D case) or `LineSubset` (1D continuum case) or `None`.
- `"winding"`: float or `None` (continuum).
- `"is_continuum"`: bool.
- `"_mu1_bracket"`: `(low, high)` bracket at convergence.
- `"_winding_bracket"`: `(w_low, w_high)` at convergence, or the straddling left/right limits at a continuum boundary.
- `"_w_limits"`: `(w_left, w_right)` — present only for continuum boundaries.
- `"_exit_reason"`: one of `"w_zero"`, `"w_zero_continuum_edge"`, `"left_endpoint_zero"`, `"right_endpoint_zero"`, `"continuum_boundary"` (max-iteration exhaustion raises `RuntimeError` instead of returning a reason).

### 3.8 `collect_GBZ_subsets`

```python
def collect_GBZ_subsets(
    coeffs: np.ndarray,
    degs: np.ndarray,
    E_ref: complex,
    perc: float = None,
    debug_mode: bool = False,
    **options,
) -> GBZResult
```

Main entry point. Builds the characteristic polynomial from `(coeffs, degs)`, solves for the `μ₁` where the average major-axis winding number vanishes, and (optionally) reclassifies the candidate as empty when a zero plateau exists next to it.

**Parameters**:
- `coeffs`: complex coefficients of the characteristic Laurent polynomial `f(E, β₁, β₂)`.
- `degs`: `(n_terms, 3)` integer exponents of `(E, β₁, β₂)` per term.
- `E_ref`: reference energy to test.
- `perc`: progress fraction in `[0, 1]`, printed as a percentage.
- `debug_mode`: if `True`, re-raise solver exceptions instead of returning a failed `GBZResult`.
- `**options`: solver options:
  - `"mu1_guess"` (default `(-1, 1)`)
  - `"zero_tol"` (1e-10)
  - `"continuum_perturb"` (1e-2)
  - `"max_iter"` (60)
  - `"plateau_check"` (True)
  - `"plateau_probe_radius"` (None)
  - `"zm_run_kwargs"` ({})
  - `"continuum_tol"` (1e-6)
  - `"crossing_tol"` (1e-10)
  - obsolete `"N_points"` / `"xtol"` / `"max_newton"` are accepted and ignored

**Returns**: `GBZResult` with connected subsets. `gbz.is_empty` means `E_ref` is outside the SGBZ spectrum. `gbz.is_continuum` means in-spectrum with 1D LineSubsets in `subsets` (`index == (0, n_1d)`). Otherwise `index == (n_0d, 0)` with `PointSubset`s.

### 3.9 Other Exports

| Function/Constant | Purpose |
|-------------------|---------|
| `Mu2MidZM` | `ZeroManager` + piecewise-smooth μ₂_mid construction (§2.1) |
| `ItemView`, `Mu2MidBreakpoint` | μ₂_mid data objects (§2.1) |
| `CONTINUUM_TOL` | Default modulus gap threshold for continuum detection (1e-6) |
| `CONTINUUM_FRAC` | Default vote fraction (0.9) |
| `WindingFun`, `get_winding_number` | Loop winding integrand + quad integration (§2.4) |
| `CharPoly` | Characteristic polynomial wrapper (from `pygbz2d.core`) |
| `get_minor_degrees` | Extract (M, N) from polynomial degrees (from `pygbz2d.core`) |
| `PointSubset`, `LineSubset`, `GBZResult`, `ConnectedSubset` | Data types (from `pygbz2d.core`) |

## 4. Key Numerical Parameters

### analyze / Mu2Mid Build

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `tie_tol` | 1e-6 | Whole-segment same-modulus clustering threshold |
| `crossing_tol` | 1e-10 | brentq refinement tolerance AND EventGroup merge gap |
| `min_direction_deriv` | 1e-12 | Minimum \|Re(V_a)−Re(V_b)\|; below → tangent/hard |
| `CONTINUUM_TOL` | 1e-6 | Continuum / ItemView clustering threshold |
| `CONTINUUM_FRAC` | 0.9 | Same-modulus vote fraction |
| `refine_multi_crossings` | True | Enable pre-crossing multi-root mesh refinement |
| `refine_max_rounds` | 3 | Maximum refinement rounds |
| `refine_safety_factor` | 4.0 | Sub-intervals per narrowest predicted root gap (ρ) |
| `refine_max_subintervals` | 64 | Uniform sub-mesh cap per interval per round |
| `refine_max_total_inserts` | 2000 | Total refinement insertion budget |

### Crossing Detection

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `crossing_tol` | 1e-10 | brentq `xtol`; close-event merge gap |
| `MIN_DIRECTION_DERIV` | 1e-12 | Tangent-direction protection floor (hard boundary) |

### Bisection Solver

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `zero_tol` | 1e-10 | Winding zero threshold |
| `continuum_perturb` | 1e-2 | μ₁ perturbation scale for continuum resolution (scales 1/2/4/8 tried) |
| `max_iter` | 60 | Maximum bisection iterations |

### Plateau Detection

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `plateau_check` | True | Enable zero-plateau detection |
| `plateau_probe_radius` | None | Probe radius (default: auto from `continuum_perturb`) |

## 5. Relation to `pygbz2d.amoeba`

| Aspect | pygbz2d.sgbz | pygbz2d.amoeba |
|--------|------------------|---------------------|
| Backend | `continuation.ZeroManager` (adaptive β₂-root tracking) | `continuation.ZeroManager` (same backend) |
| Base manifold | `μ₂_mid(θ₁)` piecewise-smooth boundary-pair mean curve | `μ₂` = constant level surface |
| Root ordering | By `|β₂|` at each `θ₁` mesh row, via `ItemView` (stable under continuum cluster snapping) | Hungarian matching across segments |
| Integration direction | `θ₁` only (average major-axis winding) | Both `θ₁` and `θ₂` (Ronkin winding) |
| Continuum criterion | `j_lo == j_hi` (whole-segment same-modulus ItemView) | `\|ln\|β₂\| − μ₂\| ≈ 0` + Ronkin minimum |
| Crossing detection | Per-column vs μ₂_mid, cubic-Hermite bracketing (not Newton) | Per-column vs μ₂, 2D fsolve |
| Dedup | Rule 1 (continuum-endpoint snap) + Rule 2 (zero-identity-key); no pairwise `NotImplementedError` | Same Rule 1/2 |
| Charge classification | Inline `sign(g')` where `g = ln|β_j| − μ₂_mid` (read atomically at the crossing row) | From Ronkin gradient at crossing |
| Loop path | Exact piecewise-smooth μ₂_mid (cubic Hermite between mesh rows, linear near MR) | Ronkin-image level curve |
| Winding integration | `quad` of `Im[f'/f]` split at μ₂_mid breakpoints | `quad` of `Im[f'/f]` over Ronkin-image loop |
| Solver method | Plain midpoint bisection (flat winding plateaus) | Plain midpoint bisection (same) |
| LineSubset materialization | `_find_boundary_runs` (contiguous `j_lo==j_hi==item` runs) + `_join_runs_across_mrs` | Whole-segment join via `_is_cluster_endpoint` |
| Spectral inclusion | `σ_SGBZ ⊆ σ_Amoeba` | — |
