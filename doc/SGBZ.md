# brute_force_SGBZ — SGBZ Spectrum Calculation

Non-Hermitian spectrum computation based on the SGBZ (Strip Generalized Brillouin Zone) formulation. The 2026-08-13 μ₂_mid framework rewrite unifies SGBZ crossing detection with amoeba's per-column-vs-μ₂ structure — the only difference is that μ₂ is the piecewise-smooth boundary-pair mean curve `μ₂_mid(θ₁)` instead of a constant. See `log/2026-08-13-SGBZ算法梳理.md` for the authoritative algorithm writeup.

## 1. Theoretical Background

### 1.1 SGBZ and Strip Winding Number

SGBZ integrates along the `θ₁` direction, with the base manifold `μ₂ = μ₂_mid(θ₁; E, μ₁)` — the piecewise-smooth mean of the two boundary-pair moduli `ln|β_{j_lo}|` and `ln|β_{j_hi}|`, where `j_lo`/`j_hi` are the item indices occupying sorted positions `M-1`/`M` (0-based). The average major-axis winding number `W(E, μ₁)` is the winding number around the `θ₂` direction integrated over `θ₁`; its zero in `μ₁` defines the SGBZ boundary.

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

**`build_mu2_mid`** — the 4-stage build:

| Stage | Action |
|-------|--------|
| 0 — walls | Insert a point at `wall_frac·grid` on each side of near-tie rows (MR rows + rows where any of the three sort-adjacent pairs `(M-2,M-1)/(M-1,M)/(M,M+1)` is within `tie_tol`). Walls give sort-change detection clean finite-tangent endpoints. |
| 1 — sort-change detection | Scan adjacent rows for `j_lo`/`j_hi` (item-index) changes outside walls → events `(t_mid, a, b, t_lo, t_hi, pair_kind, seg, row)`. Dedup by `(a, b, round(t_mid, 8))`. Same-item pairs (`j_lo == j_hi`, continuum) are skipped. |
| 2 — cubic-Hermite bracketing | Refine each sort-change crossing `f = ln|β_a| − ln|β_b| = 0` by cubic-Hermite bracketing (not Newton): build the cubic from the endpoints' `(value, f' = Re(V_a) − Re(V_b))`, `np.roots` predicts `θ_pred`, `insert_solution` takes the true `f_pred` and true derivative, **derivative sign (direction) + `f_pred` sign** fixes which side of the root `θ_pred` is on, tighten bracket. Converges on bracket width `< xtol`. Robust where Newton diverges (branch points) and where a sign-change check deadlocks (`f_pred` agreeing with both ends). Linear fallback when an endpoint tangent is singular (MR) or absent. |
| 3 — assemble + breakpoints | Flatten the mesh into `mu2_mid_theta1/values/derivs/jlo/jhi` arrays (per-segment arrays `seg_mu2_values`/`seg_mu2_derivs` also stored so the post-build crossing detector treats μ₂_mid as a first-class curve). Breakpoints from sort-change rows (`_find_breakpoints`) and MR rows (`_find_mr_breakpoints`). |

**Inline continuum detection**: after building, `has_continuum` is read directly off the ItemView — `j_lo == j_hi` (mult ≥ 2) over a fraction of a segment means the boundary pair is one continuum item. So the bisection gate *is* the build itself (§1/§6.3), with no separate two-point detector and its false-negative risk. The `_detect_continuum_clusters_internal` whole-segment same-modulus criterion (BFS transitive closure over column pairs with `max over rows |ln|β_j| − ln|β_k|| < tie_tol`) covers all three sort-adjacent pairs, not just the boundary pair — their continua also destabilise `abs_argsort`.

**Post-build insert hook**: `insert_solution` is overridden to refresh the μ₂_mid representation after a post-build insertion (e.g. `_bracket_crossing` refining a multi-index crossing the build's pairwise bracket missed). During `build_mu2_mid` itself the refresh is skipped — Stage 3 `_assemble` rebuilds from the final mesh anyway.

### 2.2 `continuum_lines.py` — Continuum Detection and LineSubset Materialization

Continuum detection is folded into the μ₂_mid build (`Mu2MidZM.has_continuum`). This module provides two entry points:

**`detect_continuum_simple`** — a presence-only flag for callers that hold a plain `ZeroManager` and only need "is there a continuum?" (the plateau probe). It builds a fresh `Mu2MidZM` and returns `has_continuum`. The μ₁ bisection itself does NOT call this — it constructs its own `Mu2MidZM` and reads `has_continuum` inline.

**`extract_continuum_linesubsets`** — LIVE materialization of the 1D continuum tracks as `LineSubset`s (not TODO). SGBZ LineSubset semantics differ from amoeba's: an SGBZ LineSubset is the stretch where a specific boundary pair (sorted positions M-1/M) is the degenerate continuum item, terminating at EITHER a modulus-sort change (another root overtakes the boundary pair — a hard terminator *inside* a segment) OR a multiple root (MR) at a segment boundary. So a piece is a **contiguous run** of `j_lo == j_hi == item` rows, NOT the whole segment.

- `_find_boundary_runs`: split each segment's `j_lo == j_hi == item` rows into maximal contiguous runs (multiple runs per segment are possible — the same continuum pair can drop out of the boundary and re-enter later).
- `_runs_to_pieces`: one `_LinePiece` per run per continuum track (each cluster column is a distinct β₂ curve at the same `|β₂|`).
- `_join_runs_across_mrs`: join pieces whose endpoints touch a segment edge (MR / θ₁=0≡2π seam). An endpoint strictly inside a segment is a sort-change terminator — no join. At an MR, if the endpoint root is in the MR's cluster the track terminates; otherwise it continues into the adjacent segment and is matched by root value. Iterated to a fixpoint so chains and the cyclic seam converge.

### 2.3 `crossings.py` — Crossing Detection

Per-column vs μ₂_mid crossing detection (§2), mirroring amoeba's `extract_amoeba_subsets` structure — the only difference is μ₂ is the curve `μ₂_mid(θ₁)` instead of a constant.

**Why per-column has no false positives** (§2.2): μ₂_mid sits in the gap between the boundary-pair moduli, where no other track lives. A track crossing μ₂_mid must become a boundary column at the crossing instant ⟺ the boundary gap closes ⟺ an SGBZ point. So "track j crosses μ₂_mid" is *equivalent* to the SGBZ point set — including multi-index jumps (a track leaping M+3→M-2 still crosses μ₂_mid once and is caught), which the old sort-change bracket missed.

**Algorithm** (`detect_crossings_simple`):
1. **Build μ₂_mid** (`Mu2MidZM.build_mu2_mid`): walls + sort-change refinement via cubic-Hermite bracketing. The post-build mesh is clean and the μ₂_mid curve is available as the first-class arrays `seg_mu2_values` / `seg_mu2_derivs`.
2. **Per-column detection (post-build)**: for each track `j` compute `g_j = ln|β_j| − μ₂_mid_values` against the BUILT curve (track-vs-curve, numerically distinct from a track-vs-track gap comparison). A genuine **sign change** of `g_j` → bracket & refine via cubic-Hermite bracketing. No near-zero / near-miss shortcut — that used to fire spuriously where `g_j` hovers at machine-ε (the degenerate seam) without an actual crossing. Exact touches (`g == 0` at a mesh row) are read directly off the row.
3. **MR echo drop**: crossings within `_MR_PROXIMITY_TOL` of a boundary MR (an MR whose cluster covers sorted positions M-1 AND M) are dropped and replaced by the exact ZeroManager MR record (`_mr_boundary_entries`).
4. **Charge + PointSubset construction**: charge = `sign(g')` where `g' = Re(V_j) − μ₂_mid_derivs` (the track-vs-curve direction). `β₂`, `g'`, and the charge are read **atomically** from the root row at the instant it is identified — before any later bracket on another track can insert rows and shift indices. Two `PointSubset` per crossing (the boundary β₂ pair) plus one per boundary MR.

**Dedup** (§2.5): the old pairwise `raise NotImplementedError` is gone. Under the per-column-vs-μ₂_mid formulation an interior crossing's identity key `(seg, interval, track j)` is naturally unique — distinct tracks at the same `θ₁` have different `β₂` and do not collide. Dedup is only the two inherited amoeba rules: Rule 1 — a crossing landing within `snap_tol` of a continuum LineSubset endpoint with matching β₂ is dropped (continuum/MR boundary, not a discrete zero); Rule 2 — `d==0` exact touches deduplicated by zero-identity-key. MR-echo replacement (step 3) subsumes Rule 1 in practice.

**Charge classification** (`_classify_charge`):
- **Ordinary** (charge ±1): `charge = sign(g')` where `g = ln|β_j| − μ₂_mid`. At the crossing `j` is a boundary column so `sign(g') = sign(½ gap')`.
- **MR** (charge 0): near a multiple root — branch point, charge unknown (not physically 0).
- **Tangent** (charge 0): `g'` not finite — tangency, physically 0 but treated as hard to be safe against near-tangency transversal pairs.
- **Unknown** (charge 0): `g'` not finite.

MR/tangent/unknown act as **hard** region boundaries (delimit regions, do not propagate winding); ordinary is **soft** (propagates via charge).

### 2.4 `winding.py` — Average Winding Computation

Computes the average major-axis winding number `W(E_ref, μ₁)` from a built `Mu2MidZM` at fixed `(E_ref, μ₁)` together with the charge list from `crossings.py`.

**Loop path** (§6.4): `β₂ = exp(μ₂_mid(θ₁) + iθ₂)` with `θ₁ ∈ [0, 2π)` and fixed `θ₂`. The μ₂_mid path is the *exact* piecewise-smooth curve (not a polyline approximation): `_Mu2MidPath` builds per-row cubic Hermite from `(value, derivative)`, using `deriv_left`/`deriv_right` at breakpoint rows (left cubic uses `deriv_left`, right uses `deriv_right`); MR rows (`deriv == inf`) and rows with unusable tangents fall back to linear interpolation.

**Integration**: `WindingFun` computes `Im[f'(t)/f(t)]` along the loop; `get_winding_number` integrates via `scipy.integrate.quad`, **split at the μ₂_mid breakpoints** so every quad segment lies on one smooth piece. The analytic `dβ₂/dθ₁ = β₂ · μ₂_mid'(θ₁)` comes from the path's cubic-Hermite derivative.

**Region/seed/charge propagation** (§3.2): crossings partition `θ₂` into regions (hard boundaries) and intervals (soft boundaries). One seed interval per region, picked by `_pick_seed_theta2` to maximize `min |f(E, β₁(θ₁), β₂_loop(θ₁))|` over `θ₁` (the true safety metric — farthest from char-poly zeros, accounting for β₁ and `|∂f/∂β₂|` that a β₂-distance proxy ignores; a cheap `min |β₂_loop−β₂_root|` pre-screen filters candidates). Compute `w₀` there via `_loop_winding_quad`, propagate across the region's soft boundaries via charges (ordinary: ±1). Each region contributes exactly one loop-winding evaluation; ±1 numerical noise is confined to the per-region seed.

### 2.5 `sgbz_solver.py` — μ₁ Bisection Solver

Top-level solver using bracket expansion + plain midpoint bisection with continuum interception.

- **`_evaluate_winding`**: build a fresh `Mu2MidZM` at the probed `μ₁` and run `build_mu2_mid` — the build both detects the continuum inline (`has_continuum`) and provides the piecewise-smooth path the winding integral needs, so the root-solving cost is amortised (§6.3). When a continuum is detected `W` is undefined → returns `(None, None, zm)`; the caller resolves it via left/right limits and, if it is the boundary, materialises the LineSubsets from the built `zm`. Otherwise runs crossing detection + winding on the same built `zm` (no second build).
- **Bracket expansion**: expand left/right until winding signs straddle zero.
- **Bisection**: plain midpoint (not false-position) because the winding has flat plateaus (±1) with a narrow transition zone; false position stalls on this shape while midpoint guarantees bracket halving.
- **Continuum interception** (`handle_continuum`): when `w_mid` is `None` (continuum), resolve left/right limits via `_resolve_continuum_winding` (perturbation scales `(1, 2, 4, 8) × continuum_perturb`); if they straddle zero, that `μ₁` IS the SGBZ boundary and the 1D LineSubsets are materialised by `extract_continuum_linesubsets` from the built `zm`.
- **Convergence check**: only applied when `w_mid` is a real winding number. When `w_mid` is `None` (continuum proxy), the algorithm continues iteration without convergence testing — converging on a proxy would confuse "left limit is zero" with "W is zero".

### 2.6 `plateau.py` — Zero-Plateau Detection

Two-stage check, unaffected by the μ₂_mid framework:

1. **Clustering pre-check** (`_check_pmgbz_points_clustered`): thin SGBZ adapter over `gbz_types.check_points_clustered_on_torus`. Euclidean distance on the `(θ₁, θ₂)` torus — **must consider both θ₁ and θ₂**, not θ₁ alone: the same `θ₁` can host multiple distinct `β₂` (degenerate pairs), so a θ₁-only check would misclassify. Example: gain-loss Haldane at E=0.5 has 6 points in 3 degenerate pairs; torus minimum distance ≈ 33°, correctly judged non-clustered.
2. **Probe ladder** (`_probe_zero_plateau_near_mu1`): thin adapter over `gbz_types.probe_zero_plateau`. Probe `μ₁ ± step`; a plateau shows zero winding with empty GBZ on both sides.

### 2.7 Continuum Handling Summary

When `|β_M| = |β_{M+1}|` holds identically over a continuous `θ₁` interval:
- `Mu2MidZM.has_continuum` is `True` → winding undefined at that `μ₁`.
- `_resolve_continuum_winding` computes left/right limits via perturbation.
- If limits straddle zero: that `μ₁` is the SGBZ boundary (1D LineSubset case); `extract_continuum_linesubsets` materializes the LineSubsets from the built `zm`. `collect_GBZ_subsets` returns `GBZResult(is_continuum=True, subsets=[LineSubset, ...], index=(0, n_1d))`.
- If limits do not straddle zero: the bisection continues using the left limit as a proxy.

## 3. API Reference

### 3.1 `Mu2MidZM` and `build_mu2_mid`

```python
class Mu2MidZM(ZeroManager):
    def __init__(self, poly: CharPoly, E_ref: complex, mu1: float): ...

    def build_mu2_mid(
        self,
        continuum_clusters: list | None = None,
        *,
        tie_tol: float = 1e-8,
        wall_frac: float = 0.1,
        xtol: float = 1e-12,
        max_iter: int = 50,
        verbose: bool = False,
    ) -> None
```

`ZeroManager` + piecewise-smooth μ₂_mid construction. Usage: `zm = Mu2MidZM(poly, E_ref, mu1); zm.run(); zm.build_mu2_mid()`.

`continuum_clusters`: per-segment list of column-tuples. `None` (default) runs the cheap internal whole-segment same-modulus detection (`_detect_continuum_clusters_internal`), which covers all three sort-adjacent pairs `(M-2,M-1)/(M-1,M)/(M,M+1)` — not just the boundary pair — because their continua also destabilise `abs_argsort`.

After `build_mu2_mid` the following are available:
- `mu2_mid_theta1` / `mu2_mid_values` / `mu2_mid_derivs` / `mu2_mid_jlo` / `mu2_mid_jhi`: flat arrays over the whole mesh.
- `seg_mu2_values` / `seg_mu2_derivs`: per-segment arrays (used by the post-build crossing detector).
- `mu2_mid_breakpoints`: list of `Mu2MidBreakpoint` (sort-change + MR).
- `has_continuum`: bool — inline continuum detection.
- `_item_views`: per-segment `ItemView`.

`insert_solution` is overridden to refresh the μ₂_mid representation after a post-build insertion (a no-op during `build_mu2_mid` itself).

### 3.2 `detect_continuum_simple`

```python
def detect_continuum_simple(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    continuum_tol: float = 1e-6,
    dV_tol: float = 1e-3,
    vote_frac: float = 0.9,
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

Materialise the 1D continuum LineSubsets of `zm`. Returns `[]` when `zm` is not a continuum. Otherwise: one `LineSubset` per continuum track per boundary run, joined across MR boundaries where the track passes through as a non-cluster root. Each LineSubset's `θ₁` range is exactly where its track held the M-1/M boundary — sort-change terminators inside a segment end the piece.

### 3.4 `detect_crossings_simple`

```python
def detect_crossings_simple(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    crossing_tol: float = 1e-10,
    detect_threshold: float = 1e-2,
    max_newton: int = 100,
    dedup_tol: float = 1e-6,
    zm_run_kwargs: dict | None = None,
) -> tuple[list[PointSubset], list[dict]]
```

0D PMGBZ-boundary crossing detection + charge classification. Per-column vs μ₂_mid (§2): builds μ₂_mid, then detects each zero-curve (track) crossing the μ₂_mid curve. Assumes no continuum — gate with `detect_continuum_simple` (or `Mu2MidZM.has_continuum`) first.

`max_newton` is a backward-compat name for the cubic-Hermite bracketing iteration limit (`_MAX_BRACKET_ITER = 100`). `detect_threshold` is retained for API symmetry but the near-miss shortcut it gated is removed — only genuine sign changes of `g_j = ln|β_j| − μ₂_mid` are refined.

**Returns**:
- `subsets`: list of `PointSubset` (two per crossing plus one per boundary MR, in detection order).
- `charges`: list of charge dicts with keys `theta1`, `theta2`, `charge`, `kind` (`'ordinary'` / `'mr'` / `'tangent'`).

**Warning**: This function MUTATES `zm` via `Mu2MidZM.build_mu2_mid` and the bracket's `insert_solution`. Run continuum detection BEFORE calling this function; do not call it twice on the same `zm`.

### 3.5 `detect_crossings_and_winding`

```python
def detect_crossings_and_winding(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    crossing_tol: float = 1e-10,
    detect_threshold: float = 1e-2,
    max_newton: int = 100,
    dedup_tol: float = 1e-6,
    zm_run_kwargs: dict | None = None,
) -> tuple[list[PointSubset], float]
```

Convenience function: crossing detection + average major-axis winding in one call. Builds μ₂_mid once (in `detect_crossings_simple` via `_ensure_mu2mid`); the winding reuses that built `Mu2MidZM`. Returns `(subsets, W_avg)`.

### 3.6 `compute_average_winding`

```python
def compute_average_winding(
    zm: ZeroManager,
    poly: CharPoly,
    M: int,
    charges: list[dict],
) -> float
```

Compute the average major-axis winding number `W(E_ref, μ₁)` from a built `Mu2MidZM` and the charge list from `crossings.py`. Partitions the `θ₂` circle into regions delimited by hard boundaries (MR/tangent/unknown, charge 0), picks one seed interval per region (maximizing `min |f|` along the loop), computes `w₀` via `_loop_winding_quad`, and propagates across the region's soft boundaries (ordinary, charge ±1). Arc-weighted mean of the per-interval windings.

### 3.7 `solve_SGBZ_for_E`

```python
def solve_SGBZ_for_E(
    poly: CharPoly,
    E_ref: complex,
    mu1_guess: tuple[float, float] = (-1, 1),
    zero_tol: float = 1e-10,
    continuum_perturb: float = 1e-2,
    max_iter: int = 60,
    xtol: float = 2e-12,
    zm_run_kwargs: Optional[dict] = None,
    *,
    continuum_tol: float = 1e-6,
    dV_tol: float = 1e-3,
    vote_frac: float = 0.9,
    crossing_tol: float = 1e-10,
    detect_threshold: float = 1e-2,
    max_newton: int = 100,
    dedup_tol: float = 1e-6,
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
- `"_exit_reason"`: one of `"w_zero"`, `"bracket_xtol"`, `"left_endpoint_zero"`, `"right_endpoint_zero"`, `"continuum_boundary"`, `"max_iter"`.

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
  - `"xtol"` (2e-12)
  - `"plateau_check"` (True)
  - `"plateau_probe_radius"` (None)
  - `"zm_run_kwargs"` ({})
  - `"continuum_tol"` (1e-6)
  - `"dV_tol"` (1e-3)
  - `"vote_frac"` (0.9)
  - `"crossing_tol"` (1e-10)
  - `"detect_threshold"` (1e-2)
  - `"max_newton"` (100)
  - `"dedup_tol"` (1e-6)

**Returns**: `GBZResult` with connected subsets. `gbz.is_empty` means `E_ref` is outside the SGBZ spectrum. `gbz.is_continuum` means in-spectrum with 1D LineSubsets in `subsets` (`index == (0, n_1d)`). Otherwise `index == (n_0d, 0)` with `PointSubset`s.

### 3.9 Other Exports

| Function/Constant | Purpose |
|-------------------|---------|
| `Mu2MidZM` | `ZeroManager` + piecewise-smooth μ₂_mid construction (§2.1) |
| `ItemView`, `Mu2MidBreakpoint` | μ₂_mid data objects (§2.1) |
| `CONTINUUM_TOL` | Default modulus gap threshold for continuum detection (1e-6) |
| `CONTINUUM_FRAC` | Default vote fraction (0.9) |
| `WindingFun`, `get_winding_number` | Loop winding integrand + quad integration (§2.4) |
| `CharPoly` | Characteristic polynomial wrapper (from `gbz_types`) |
| `get_minor_degrees` | Extract (M, N) from polynomial degrees (from `gbz_types`) |
| `PointSubset`, `LineSubset`, `GBZResult`, `ConnectedSubset` | Data types (from `gbz_types`) |

## 4. Key Numerical Parameters

### μ₂_mid Build

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `tie_tol` | 1e-8 | Modulus near-tie threshold for walling (rows this close in any sort-adjacent pair get walled) |
| `wall_frac` | 0.1 | Wall distance as a fraction of local grid spacing |
| `xtol` (build) | 1e-12 | Cubic-Hermite bracketing convergence tolerance for sort-change refinement |
| `max_iter` (build) | 50 | Max cubic-Hermite bracketing iterations per sort-change |
| `CONTINUUM_TOL` | 1e-6 | Modulus gap separating genuine continuum from transversal crossing |
| `CONTINUUM_FRAC` | 0.9 | Fraction of segment rows in-band for a continuum track |

### Continuum Detection (plateau probe path)

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `continuum_tol` | 1e-6 | Same as `tie_tol` when passed to `build_mu2_mid` |
| `dV_tol` | 1e-3 | Retained for API symmetry; inline detection uses whole-segment same-modulus criterion |
| `vote_frac` | 0.9 | Retained for API symmetry; inline detection uses whole-segment same-modulus criterion |

### Crossing Detection

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `crossing_tol` | 1e-10 | Cubic-Hermite bracketing convergence tolerance (`_CROSSING_TOL`) |
| `detect_threshold` | 1e-2 | Near-miss threshold (retained; near-miss shortcut removed — only sign changes refined) |
| `max_newton` | 100 | Max cubic-Hermite bracketing iterations (`_MAX_BRACKET_ITER`; backward-compat alias `_MAX_NEWTON_ITER`) |
| `dedup_tol` | 1e-6 | L²-distance threshold for duplicate-crossing dedup (`_DEDUP_TOL`) |
| `_TANGENT_F_TOL` | 1e-6 | `\|f\|` threshold for tangent touch classification |
| `_MR_PROXIMITY_TOL` | 1e-4 | θ₁ distance threshold for MR echo detection |
| `_TANGENCY_THRESHOLD` | 1e-3 | `\|g'\|` threshold for tangency classification |

### Bisection Solver

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `zero_tol` | 1e-10 | Winding zero threshold |
| `continuum_perturb` | 1e-2 | μ₁ perturbation scale for continuum resolution (scales 1/2/4/8 tried) |
| `max_iter` | 60 | Maximum bisection iterations |
| `xtol` | 2e-12 | Minimum μ₁ bracket width for convergence |

### Plateau Detection

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `plateau_check` | True | Enable zero-plateau detection |
| `plateau_probe_radius` | None | Probe radius (default: auto from `continuum_perturb`) |

## 5. Relation to `brute_force_amoeba`

| Aspect | brute_force_SGBZ | brute_force_amoeba |
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
