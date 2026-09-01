# pygbz2d.amoeba — Amoebic Spectrum Calculation

Non-Hermitian spectrum computation based on the amoeba formulation. Complementary to the SGBZ formulation (`pygbz2d/sgbz/`).

The zero-solving backend is `continuation.ZeroManager` (wrapped as `AmoebaZeroManager`). Subset assembly lives in `amoeba.py`; crossing detection and winding primitives live in `zm_extract.py` and `ronkin_winding.py`.

## 1. Theoretical Background

### 1.1 Amoeba Definition

For an $n$-variable Laurent polynomial $p(\boldsymbol{\beta})$, its **amoeba** is:

$$\mathcal{A}_p \equiv \left\{ \ln|\boldsymbol{\beta}| \;\middle|\; p(\boldsymbol{\beta}) = 0 \right\}$$

For a 2D non-Hermitian lattice at fixed reference energy $E$, the characteristic polynomial $f(E, \boldsymbol{\beta}) = \det[E - h(\boldsymbol{\beta})]$ is a bivariate Laurent polynomial, with amoeba $\mathcal{A}(E) = \mathcal{A}_{f(E,\cdot)}$.

### 1.2 Ronkin Function and Winding Numbers

**Ronkin function**:

$$R_p(\boldsymbol{\mu}) = \int_{\mathbb{T}^n} \left(\frac{d\boldsymbol{\theta}}{2\pi}\right)^n \ln|p(e^{\boldsymbol{\mu}+i\boldsymbol{\theta}})|$$

Its gradient is related to the **average winding numbers** $u_j$. In a hole of the amoeba the $u_j$ are topological invariants.

- **Central hole**: a hole where $\nabla R_E = 0$.
- $E \in \sigma_{\text{Amoeba}}$: $\mathcal{A}(E)$ has **no** central hole.
- $E \notin \sigma_{\text{Amoeba}}$: $\mathcal{A}(E)$ **has** a central hole.

### 1.3 Relation to SGBZ

| Aspect | SGBZ | Amoeba |
|--------|------|--------|
| Base manifold | $X(E,\mu_1)$, $\mu_2$ varies with $\theta_1$ | Level surface $\mu_2 = \text{const}$ |
| Root tracking | Fixed $\theta_1$ mesh + pairwise events | `continuation.ZeroManager` adaptive $\beta_2$-root tracks |
| Continuum detection | ItemView cluster vote | `std(ln|\beta_2|) < tol` per track per segment |
| Subset output | `collect_GBZ_subsets` in `sgbz_solver.py` | `collect_GBZ_subsets` in `amoeba.py` |

## 2. Algorithm Architecture

### 2.1 Zero-Manager Backend

`AmoebaZeroManager(ZeroManager)` caches `seg_logabs[s] = log|segments[s].tracked_roots|` after `run()`. The cache is refreshed after mesh insertion (`insert_solution`). The ZM is $\mu_2$-independent — built once per `(E, \mu_1)` and reused across all $\mu_2$ evaluations.

### 2.2 Continuum Detection (pre-bisection)

`detect_continuum(zm, tol)` detects, for every segment and every track column, whether `std(ln|β₂|[:, j]) < tol`. Detected flat tracks are returned as groups:

```
[(mu2_c, [(seg_idx, col_idx), ...]), ...]
```

where nearby `mu2_c` values (closer than `tol`) are merged. The `(seg_idx, col_idx)` list is the continuum membership consumed later by subset assembly — it is recorded at detection time and passed through the bisection result.

### 2.3 Crossing Detection

`find_crossings(zm, mu1, mu2, avoided_segments=None, return_refined=False)` finds all crossings of `ln|β₂| = μ₂`:

- Traverses every segment and every column.
- Skips `(seg_idx, col_idx)` entries in `avoided_segments`.
- First detects exact touches with `d[:-1] == 0`, then sign changes with `d[1:] * d[:-1] < 0`.
- Returns `list[(β₁, β₂)]`. Linear interpolation unless `return_refined=True`, in which case `_find_exact_crossing` (fsolve + analytic Jacobian) refines each sign-change crossing, falling back to the linear estimate on failure.

### 2.4 Winding

`calculate_a2_average_winding(zm, mu1, mu2, avoided_segments=None, return_refined=False)` is a thin consumer of `find_crossings`. It converts crossings back to angular zeros and computes the a2 average winding with `_get_average_winding_from_zeros(direction=2)`.

`_get_average_winding_from_zeros` partitions the angular circle in the target direction by the zero coordinates and sums `u × width / (2π)`; it also returns the normalized non-zero area used by the plateau pre-check.

### 2.5 Inner μ₂ Solve (`_find_mu2_for_w2_zero`)

Continuum-first dispatcher:

1. Build/reuse `AmoebaZeroManager`.
2. `_try_fast_mu2` — θ₁=0 gap test. Sort roots at θ₁=0 by `|β₂|`, take the first M and last N columns, refine their extrema with Hermite-predicted mesh insertion, and compute `A = max(ln|β₂|` over first M`)`, `B = min(ln|β₂|` over last N`)`. If `A < B`, return `μ₂ = (A+B)/2` with no crossings.
3. `detect_continuum` — if flat tracks exist, probe each merged `μ₂_c` at `μ₂_c ± ε` (`ESCAPE_LADDER`). Opposite w2 signs → return the continuum boundary and the matched `(seg_idx, col_idx)` members. Equal signs → tighten the μ₂ bracket with the signed probe points.
4. `_bisect_mu2_discrete` — pure two-stage bisection on the tightened bracket:
   - coarse stage with `return_refined=False` and `BISECT_COARSE_XTOL`;
   - fine stage on the coarse bracket with `return_refined=True` and `BISECT_XTOL`; the fine bracket is expanded outward when refined endpoint windings no longer straddle zero.

### 2.6 Outer μ₁ Bisection (`bisect_amoeba_ronkin_min`)

- Endpoint expansion and midpoint bisection on μ₁.
- Every inner return is checked for `is_continuum`.
- `_handle_continuum` resolves a continuum inner result: at `μ₁ ± ε` it builds a fresh ZM, runs `find_crossings` at the fixed `μ₂_c`, and computes the a1 average winding from those crossings. Opposite w1 signs end the outer bisection; equal signs update the μ₁ bracket.
- Discrete path computes w1 from the inner `zeros` and bisects normally.

### 2.7 Subset Assembly

`collect_GBZ_subsets` runs the bisection, then:

- **`is_continuum == False`**: each `(θ₁, θ₂)` zero becomes a `PointSubset`.
- **`is_continuum == True`**:
  1. materialize the recorded continuum `(seg_idx, col_idx)` members as per-segment `JoinableLinePiece`s;
  2. splice pieces that meet end-to-end (interior MRs share a track frame; the θ=0≡2π seam is translated through `boundary_perm`; MR-cluster endpoints are genuine terminators);
  3. run `find_crossings` with the continuum members passed as `avoided_segments` to find discrete crossings;
  4. convert those crossings to `PointSubset`s, dropping any point whose θ₁ is within `SNAP_TOL` of a continuum `LineSubset`.

### 2.8 Plateau Detection

Non-continuum results still pass through the plateau pre-check/probe (tiny w1/w2 non-zero area + clustered zeros). Continuum results skip plateau detection.

## 3. API Reference

### 3.1 `collect_GBZ_subsets`

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

Main entry point. Builds `CharPoly`, runs `bisect_amoeba_ronkin_min`, then assembles subsets from the solved ZM. `debug_mode=True` re-raises exceptions instead of returning a failed `GBZResult`.

Options include `mu1_low`, `mu1_high`, `mu2_low`, `mu2_high`, `continuum_tol`, `continuum_perturb`, `max_iter`, `xtol`, `max_range_expansions`, `range_expand_factor`, `plateau_check`, `plateau_winding_tol`, `plateau_probe_radius`, `plateau_area_threshold`, `plateau_cluster_tol`.

### 3.2 `bisect_amoeba_ronkin_min`

```python
def bisect_amoeba_ronkin_min(
    char_poly: CharPoly,
    E_ref: complex,
    mu1_low: float = -1,
    mu1_high: float = 1,
    mu2_low: float = -1,
    mu2_high: float = 1,
    ...
) -> dict:
```

Returns `mu1`, `mu2`, `zeros`, `is_continuum`, `_mu1_bracket`, `_w1_bracket`, `_exit_reason`, `_w1_area` (discrete path), `_continuum_members` (continuum path), `_zm`.

### 3.3 `find_crossings`

```python
def find_crossings(
    zm: ZeroManager,
    mu1: float,
    mu2: float,
    avoided_segments: Optional[list[tuple[int, int]]] = None,
    return_refined: bool = False,
) -> list[tuple[complex, complex]]:
```

### 3.4 `calculate_a2_average_winding`

```python
def calculate_a2_average_winding(
    zm: ZeroManager,
    mu1: float,
    mu2: float,
    avoided_segments: Optional[list[tuple[int, int]]] = None,
    return_refined: bool = False,
) -> float:
```

### 3.5 `AmoebaZeroManager`

```python
class AmoebaZeroManager(ZeroManager):
    seg_logabs: list[np.ndarray]

    def __init__(self, poly: CharPoly, E_ref: complex, mu1: float): ...
    def run(self, **kwargs) -> None: ...
    def refresh_logabs(self) -> None: ...
```

## 4. Internal Functions

### `zm_extract.py`

| Function | Purpose |
|----------|---------|
| `detect_continuum` | std-based continuum detection; returns merged `(μ₂_c, members)` groups. |
| `find_crossings` | All `ln|β₂| = μ₂` crossings with optional refine and segment/column avoidance. |
| `calculate_a2_average_winding` | a2 average winding from crossings. |

### `bisect.py`

| Function | Purpose |
|----------|---------|
| `_try_fast_mu2` | θ₁=0 gap test; returns `{ok, A, B}`. |
| `_screen_extremum_intervals` / `_refine_track_extrema` | Vectorized screen + Hermite-predicted mesh insertion for extremum refinement. |
| `_bisect_mu2_discrete` | Two-stage (coarse/fine) μ₂ bisection. |
| `_find_mu2_for_w2_zero` | Continuum-first inner μ₂ solve. |
| `_handle_continuum` | μ₁ ± ε resolution of a continuum inner result. |
| `bisect_amoeba_ronkin_min` | Outer μ₁ bisection. |

### `amoeba.py`

| Function | Purpose |
|----------|---------|
| `_assemble_discrete_subsets` | zeros → `PointSubset`. |
| `_assemble_continuum_subsets` | continuum members → spliced `LineSubset` + discrete `PointSubset`. |
| `_splice_continuum_pieces` | End-to-end splicing of per-segment continuum pieces. |
| `_point_near_any_line` | θ₁ snap screen for dropping points near a `LineSubset`. |
| `_check_zeros_are_clustered` / `_probe_zero_plateau_near_mu1` | Plateau detection. |

### `ronkin_winding.py`

| Function | Purpose |
|----------|---------|
| `_find_exact_crossing` | fsolve + analytic Jacobian refinement of one crossing. |
| `_get_average_winding_from_zeros` | Average winding and non-zero area from a zero partition. |

## 5. Package Structure

```
pygbz2d/amoeba/
├── __init__.py            # Public API exports
├── amoeba.py              # collect_GBZ_subsets, subset assembly, plateau detection
├── bisect.py              # μ₁/μ₂ bisection, fast μ₂ gap test, extremum refinement
├── ronkin_winding.py      # fsolve crossing refinement, average winding from zeros
└── zm_extract.py          # AmoebaZeroManager, detect_continuum, find_crossings,
                           #   calculate_a2_average_winding
```

Public API: `CharPoly`, `bisect_amoeba_ronkin_min`, `AmoebaZeroManager`, `calculate_a2_average_winding`, `collect_GBZ_subsets`, `PointSubset`, `LineSubset`, `GBZResult`, `ConnectedSubset`.
