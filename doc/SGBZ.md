# brute_force_SGBZ — SGBZ Spectrum Calculation

Non-Hermitian spectrum computation based on the SGBZ (Strip Generalized Brillouin Zone) formulation.

## 1. Theoretical Background

### 1.1 SGBZ and Strip Winding Number

SGBZ integrates along the $\theta_1$ direction, with the base manifold $\mu_2 = \rho_{2,0}(\theta_1; E, \mu_1)$ (the PMGBZ curve determined by $|\beta_M| = |\beta_{M+1}|$). The average major-axis winding number $W(E, \mu_1)$ is the winding number around the $\theta_2$ direction integrated over $\theta_1$.

### 1.2 PMGBZ Point Classification

On the $\theta_1$ axis, $\beta_2$ roots are sorted by modulus. When $|\beta_M| \approx |\beta_{M+1}|$, a degeneracy occurs — this is called a PMGBZ point:

- **Continuum**: $|\beta_M| = |\beta_{M+1}|$ holds over a continuous $\theta_1$ interval
- **Accidental**: $|\beta_M| = |\beta_{M+1}|$ only at isolated $\theta_1$ positions

Each accidental PMGBZ point returns three tuples:
```python
beta2_sols = [pos_tuple, neg_tuple, zero_tuple]
```
Sorted by $d\mu_2/d\theta_1$ in descending order, the top $M - i_{\min}$ are labeled as positive charge.

### 1.3 Spectral Inclusion Relation

$$\sigma_{\text{Amoeba}} \supset \cup_j \sigma_{\text{SGBZ},j}$$

In general, the SGBZ spectrum is a subset of the amoeba spectrum. For uniform bands (geometry-dependent spectral degeneracy), the two are equal.

## 2. Key Components

### 2.1 `continuum.py` — Continuum Detection

Two-stage detection operating directly on `ZeroManager` segments:

**Stage 1**: Candidate scan (`_continuum_runs`)
- Per-row check at sorted positions M-1 and M (0-based)
- Both (a) modulus gap < `continuum_tol` and (b) $\frac{d(\ln|\beta_2|)}{d\theta_1}$ match within `dV_tol`
- Runs are maximal consecutive row intervals with ≥ 2 rows

**Stage 2**: Whole-column voting (`_pair_vote`)
- Fraction of rows across the ENTIRE segment where $|\ln|\beta_{2,j}| - \ln|\beta_{2,k}|| < $ `continuum_tol` must exceed `vote_frac`
- Real-analyticity means genuine pairs vote near 1.0; accidental near-crossings vote near 0

Two versions:
- `detect_continuum_simple`: presence-only flag, used by μ₁ bisection
- `detect_continuum_full`: per-segment degeneracy clusters (for future §2.4 full extractor)

### 2.2 `crossings.py` — Crossing Detection

Simplified 0D PMGBZ-boundary crossing detection with inline charge classification.

**Two-phase batch algorithm**:
1. **Snapshot + Cubic Construction**: Build all cubic Hermite polynomials from endpoint snapshots before any Newton insertions (avoids cross-pair corruption)
2. **Newton Refinement**: Each prediction refined independently; charge classified from frozen $f' = \text{Re}(V_a) - \text{Re}(V_b)$

**Charge classification** (`_classify_charge`):
- **Ordinary** (charge ±1): determined by sign of $f'$
- **MR** (charge 0): near a multiple root
- **Tangent** (charge 0): $|f'|$ too small
- **Unknown** (charge 0): $f'$ not finite

MR boundary entries are detected separately (`_mr_boundary_entries`) and added without Newton refinement.

### 2.3 `winding.py` — Average Winding Computation

Computes the average major-axis winding number $W(E_{\text{ref}}, \mu_1)$ from a `ZeroManager` at fixed $(E_{\text{ref}}, \mu_1)$.

**Design**:
- Loop: $\beta_2 = \exp(\mu_{2,\text{mid}}(\theta_1) + i\theta_2)$ with fixed $\theta_2$
- No spline interpolation — the adaptive mesh IS the polyline through the data
- Crossings partition $\theta_2$ into regions (hard boundaries) and intervals (soft boundaries)
- One seed interval per region: pick $\theta_2$ maximizing distance to all roots, compute $w_0$, propagate via charges

### 2.4 `sgbz_solver.py` — μ₁ Bisection Solver

Top-level solver using bracket expansion + plain midpoint bisection.

**Bracket expansion**: Expand left/right until winding signs straddle zero
**Bisection**: Plain midpoint (not false-position) because winding has flat plateaus
**Continuum interception**: When `w_mid` is `None` (continuum), resolve left/right limits via perturbation; if they straddle zero, that $\mu_1$ IS the SGBZ boundary

### 2.5 `plateau.py` — Zero-Plateau Detection

Two-stage check:
1. **Clustering pre-check** (`_check_pmgbz_points_clustered`): Euclidean distance on $(\theta_1, \theta_2)$ torus
2. **Probe ladder** (`_probe_zero_plateau_near_mu1`): Test $\mu_1 \pm \text{step}$ for zero winding with empty GBZ

**Clustering metric**: `d = sqrt(circ_dist(θ₁_i, θ₁_j)² + circ_dist(θ₂_i, θ₂_j)²)`

Rationale: Degenerate pairs share $\theta_1$ but have different $\beta_2$, so $\theta_1$-only check would misclassify.

### 2.6 Continuum Handling

When $|\beta_M| = |\beta_{M+1}|$ holds over a continuous $\theta_1$ interval:
- `detect_continuum_simple` returns `True` → winding undefined
- `_resolve_continuum_winding` computes left/right limits via perturbation
- If limits straddle zero: that $\mu_1$ is the SGBZ boundary (1D LineSubset case)
- **LineSubset extraction is TODO**: `collect_GBZ_subsets` returns `is_continuum=True` with no subsets

## 3. API Reference

### 3.1 `detect_continuum_simple`

```python
def detect_continuum_simple(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    continuum_tol: float = 1e-6,
    dV_tol: float = 1e-3,
    vote_frac: float = 0.9,
) -> bool:
```

Simplified continuous modulus equality detection — presence only.

Two stages per segment:
1. Candidate scan: ≥ 2 consecutive rows where the M-1/M boundary pair has matching modulus and derivative
2. Whole-column voting: verify that SOME pair of columns is genuinely whole-column modulus-equal

Returns `True` as soon as any run produces one confirmed degenerate pair; false positives are skipped, the scan continues.

Used by the μ₁ bisection for an early-out (a continuum makes W undefined there).

### 3.2 `detect_continuum_full`

```python
def detect_continuum_full(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    continuum_tol: float = 1e-6,
    dV_tol: float = 1e-3,
    vote_frac: float = 0.9,
) -> list[list[tuple[int, ...]]]:
```

Full continuous modulus equality detection — per-segment clusters.

Extends the simplified version: instead of returning at the first confirmed pair, ALL candidate runs of a segment are collected, their column sets are merged, and every pair of involved columns is voted on over the whole segment.

Returns one entry per segment, in segment order: `result[seg_idx] = [(j, ...), ...]` — each cluster is a tuple of tracked_roots column indices sharing one $|\beta_2|(\theta_1)$ curve.

**Note**: Currently no caller consumes this output — the §2.4 full subset extractor is not yet implemented.

### 3.3 `detect_crossings_simple`

```python
def detect_crossings_simple(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    crossing_tol: float = 1e-10,
    detect_threshold: float = 1e-2,
    max_newton: int = 10,
    dedup_tol: float = 1e-6,
) -> tuple[list[PointSubset], list[dict]]:
```

Simplified 0D PMGBZ-boundary crossing detection + charge classification.

Assumes no continuum (caller should gate with `detect_continuum_simple` first).

**Algorithm**:
1. **Detection**: for each sub-interval, check all columns' $\ln|\beta_{2,j}| - \mu_{2,\text{mid}}$ for sign change or near-miss
2. **Solving**: for each candidate column pair (a, b), build a cubic Hermite diff polynomial, find roots via `np.roots`, then Newton-refine
3. **MR echo drop**: crossings near a boundary multiple root are dropped and replaced by the exact ZeroManager record
4. **Dedup check**: raise `NotImplementedError` if two crossings of DIFFERENT column pairs coincide
5. **Charge + PointSubset construction**: each surviving crossing is classified inline from its frozen $f'$

**Returns**:
- `subsets`: list of `PointSubset` (two per crossing plus one per boundary MR)
- `charges`: list of charge dicts with keys `theta1`, `theta2_a`, `theta2_b`, `charge`, `kind`

**Warning**: This function MUTATES `zm` (Newton iterations insert mesh rows). Run continuum detection BEFORE calling this function.

### 3.4 `detect_crossings_and_winding`

```python
def detect_crossings_and_winding(
    zm: ZeroManager,
    poly: CharPoly,
    *,
    crossing_tol: float = 1e-10,
    detect_threshold: float = 1e-2,
    max_newton: int = 10,
    dedup_tol: float = 1e-6,
) -> tuple[list[PointSubset], float]:
```

Convenience function: crossing detection + average major-axis winding in one call.

**Returns**: `(subsets, W_avg)` where `W_avg` is the average major-axis winding number.

### 3.5 `compute_average_winding`

```python
def compute_average_winding(
    zm: ZeroManager,
    poly: CharPoly,
    M: int,
    charges: list[dict],
) -> float:
```

Compute the average major-axis winding number $W(E_{\text{ref}}, \mu_1)$.

**Topology**: crossings come in two kinds:
- **Ordinary** (charge ±1, SOFT): only partitions $\theta_2$ into intervals within a region
- **MR / tangent / unknown** (charge 0, HARD): delimits regions

The algorithm partitions the circle into regions delimited by hard boundaries, picks ONE seed interval per region (farthest from all boundaries), computes $w_0$ there, and propagates across soft boundaries via charges.

### 3.6 `solve_SGBZ_for_E`

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
    max_newton: int = 10,
    dedup_tol: float = 1e-6,
) -> dict:
```

Locate the winding-zero $\mu_1$ and return solve diagnostics.

Uses bracket expansion + plain bisection (midpoint). When a continuum-degenerate $\mu_1$ is encountered, the left/right winding limits are resolved; if they straddle zero that $\mu_1$ is the SGBZ boundary (`is_continuum=True`).

**Plain bisection** is chosen over false-position methods because the winding has flat plateaus (±1) with a narrow transition zone; false position stalls on this shape while midpoint guarantees bracket halving.

**Returns**: dict with keys:
- `"mu1"`: the solution $\mu_1$
- `"subsets"`: list of `PointSubset` or `None` (for continuum)
- `"winding"`: float or `None`
- `"is_continuum"`: bool
- `"_mu1_bracket"`, `"_winding_bracket"`, `"_exit_reason"`: debug fields

**Convergence check**: Only applied when `w_mid` is a real winding number. When `w_mid` is `None` (continuum proxy), the algorithm continues iteration without convergence testing.

### 3.7 `collect_GBZ_subsets`

```python
def collect_GBZ_subsets(
    coeffs: np.ndarray,
    degs: np.ndarray,
    E_ref: complex,
    perc: float = None,
    debug_mode: bool = False,
    **options,
) -> GBZResult:
```

Main entry point. Check the SGBZ condition and return GBZ points for a reference energy.

Builds the characteristic polynomial from `(coeffs, degs)`, solves for the $\mu_1$ where the average major-axis winding number vanishes, and optionally reclassifies the candidate as empty when a zero plateau exists.

**Parameters**:
- `coeffs`: complex coefficients of the characteristic Laurent polynomial
- `degs`: `(n_terms, 3)` integer exponents of `(E, beta1, beta2)` per term
- `E_ref`: reference energy to test
- `perc`: progress fraction in `[0, 1]`, printed as a percentage
- `debug_mode`: if `True`, re-raise solver exceptions
- `**options`: solver options including:
  - `"mu1_guess"`: default `(-1, 1)`
  - `"zero_tol"`: default `1e-10`
  - `"continuum_perturb"`: default `1e-2`
  - `"max_iter"`: default `60`
  - `"xtol"`: default `2e-12`
  - `"plateau_check"`: default `True`
  - `"zm_run_kwargs"`: default `{}`

**Returns**: `GBZResult` with connected subsets. `gbz.is_empty` means $E_{\text{ref}}$ is outside the SGBZ spectrum. `gbz.is_continuum` means in-spectrum but `LineSubset` extraction is TODO.

### 3.8 Other Exports

| Function/Constant | Purpose |
|-------------------|---------|
| `CONTINUUM_TOL` | Default tolerance for continuum detection (1e-6) |
| `CONTINUUM_FRAC` | Default vote fraction for whole-column voting (0.9) |
| `CharPoly` | Characteristic polynomial wrapper (from `gbz_types`) |
| `get_minor_degrees` | Extract (M, N) from polynomial degrees (from `gbz_types`) |
| `PointSubset`, `LineSubset`, `GBZResult` | Data types (from `gbz_types`) |

## 4. Key Numerical Parameters

### Continuum Detection

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `continuum_tol` | 1e-6 | Modulus gap threshold for continuum candidate rows |
| `dV_tol` | 1e-3 | Derivative consistency threshold for continuum detection |
| `vote_frac` | 0.9 | Fraction of rows required for whole-column vote pass |

### Crossing Detection

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `crossing_tol` | 1e-10 | Newton convergence tolerance |
| `detect_threshold` | 1e-2 | Near-miss threshold for suspicious interval detection |
| `max_newton` | 10 | Maximum Newton iterations per crossing |
| `dedup_tol` | 1e-6 | L²-distance threshold for duplicate crossing detection |
| `_TANGENT_F_TOL` | 1e-6 | \|f\| threshold for tangent touch classification |
| `_MR_PROXIMITY_TOL` | 1e-4 | θ₁ distance threshold for MR echo detection |
| `_TANGENCY_THRESHOLD` | 1e-3 | \|f'\| threshold for tangency classification |

### Bisection Solver

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `zero_tol` | 1e-10 | Winding zero threshold |
| `continuum_perturb` | 1e-2 | μ₁ perturbation scale for continuum resolution |
| `max_iter` | 60 | Maximum bisection iterations |
| `xtol` | 2e-12 | Minimum μ₁ bracket width for convergence |

### Plateau Detection

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `plateau_check` | True | Enable zero-plateau detection |
| `plateau_probe_radius` | None | Probe radius (default: auto from bracket) |

**Note**: The `vote_frac` parameter controls the whole-column voting for continuum detection. A genuine degenerate pair should vote near 1.0; an accidental near-crossing votes near 0. The default 0.9 threshold may need tuning for models with numerical noise in the degenerate modulus curves.

## 5. Relation to `brute_force_amoeba`

| Aspect | brute_force_SGBZ | brute_force_amoeba |
|--------|------------------|---------------------|
| Backend | `continuation.ZeroManager` (adaptive β₂-root tracking) | `continuation.ZeroManager` (same backend) |
| Base manifold | $\mu_2 = \rho_{2,0}(\theta_1)$ variable curve | $\mu_2 = \text{const}$ level surface |
| Root ordering | By $|\beta_2|$ at each $\theta_1$ mesh row | Hungarian matching across segments |
| Integration direction | $\theta_1$ only (average major-axis winding) | Both $\theta_1$ and $\theta_2$ (Ronkin winding) |
| Continuum criterion | $|\beta_M| \approx |\beta_{M+1}|$ + derivative match + whole-column vote | $|\ln|\beta_2| - \mu_2| \approx 0$ + Ronkin minimum |
| Crossing detection | Cubic Hermite + Newton refinement, batch processing | Line intersection + Newton refinement |
| Charge classification | Inline from frozen $f'$ at Newton iterate | From Ronkin gradient at crossing |
| Solver method | Plain midpoint bisection (flat winding plateaus) | Plain midpoint bisection (same) |
| Spectral inclusion | $\sigma_{\text{SGBZ}} \subseteq \sigma_{\text{Amoeba}}$ | — |
