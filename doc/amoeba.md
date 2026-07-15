# brute_force_amoeba — Amoebic Spectrum Calculation

Non-Hermitian spectrum computation based on the amoeba formulation. Complementary to the SGBZ formulation (`brute_force_SGBZ/`).

## 1. Theoretical Background

### 1.1 Amoeba Definition

For an $n$-variable Laurent polynomial $p(\boldsymbol{\beta})$, its **amoeba** is defined as:

$$\mathcal{A}_p \equiv \left\{ \ln|\boldsymbol{\beta}| \;\middle|\; p(\boldsymbol{\beta}) = 0 \right\}$$

For a 2D non-Hermitian lattice at fixed reference energy $E$, the characteristic polynomial $f(E, \boldsymbol{\beta}) = \det[E - h(\boldsymbol{\beta})]$ is a bivariate Laurent polynomial, with corresponding amoeba:

$$\mathcal{A}(E) \equiv \mathcal{A}_{f(E,\cdot)}$$

Intuitively, $\mathcal{A}(E)$ is the set of all allowed non-Bloch scaling exponents $\ln|\beta|$ for a given $E$.

### 1.2 Ronkin Function and Winding Numbers

**Ronkin function**:

$$R_p(\boldsymbol{\mu}) = \int_{\mathbb{T}^n} \left(\frac{d\boldsymbol{\theta}}{2\pi}\right)^n \ln|p(e^{\boldsymbol{\mu}+i\boldsymbol{\theta}})|$$

Its gradient is related to the **average winding numbers**:

$$\frac{\partial R_p}{\partial \mu_j} = \int_{\mathbb{T}^{n-1}} \frac{d\theta_1 \cdots \widehat{d\theta_j} \cdots d\theta_n}{(2\pi)^{n-1}} \; u_j(\theta_1, \dots, \widehat{\theta_j}, \dots, \theta_n; \boldsymbol{\mu})$$

where $u_j$ is the winding number with respect to $\theta_j$ with all other variables fixed.

### 1.3 Amoeba Spectrum Criterion

In a **hole** of the amoeba — i.e., a connected open set $\boldsymbol{\mu} \notin \mathcal{A}(E)$ enclosed by $\mathcal{A}(E)$ — the winding numbers $u_j$ are topological invariants and $\nabla R_p$ remains constant.

- **Central hole**: a hole where $\nabla R_E = 0$
- $E \in \sigma_{\text{Amoeba}}$: $\mathcal{A}(E)$ has **no** central hole
- $E \notin \sigma_{\text{Amoeba}}$: $\mathcal{A}(E)$ **has** a central hole

### 1.4 Relation to SGBZ

| Aspect | SGBZ | Amoeba |
|--------|------|--------|
| Base manifold | $X(E,\mu_1)$, $\mu_2 = \rho_{2,0}(\theta_1;E,\mu_1)$ varies with $\theta_1$ | Level surface $\mu_2 = \text{const}$ |
| Root ordering | By $\vert\beta_2\vert$ | Hungarian matching for continuous tracking |
| Continuum detection | $\vert\beta_M\vert \approx \vert\beta_{M+1}\vert$ | $\vert\ln\vert\beta_2\vert - \mu_2\vert < \text{tol}$ |
| Spectral inclusion | $\sigma_{\text{Amoeba}} \supset \cup_j \sigma_{\text{SGBZ},j}$ | — |

## 2. Algorithm Architecture

### 2.1 Hungarian Matching for Root Tracking

Unlike SGBZ which sorts by $\vert\beta_2\vert$, the amoeba method uses Hungarian matching to track each root continuously:

1. Solve for $\beta_2$ roots at each point on a $\theta_1$ mesh
2. Starting from $\theta_1 = 0$, propagate matching relations point by point using the Hungarian algorithm
3. Cost function: chordal distance on the Riemann sphere (avoids $0/\infty$ false swaps)
4. Output: continuous root tracks `tracked[θ₁_idx, root_j]`

The periodic boundary ($\theta_1 = 0 \to 2\pi$) is handled with an extra Hungarian matching step to ensure root identities are preserved across the wrap-around, preventing spurious crossings at the periodic boundary.

### 2.2 Crossing Detection and Lazy Refinement

For each root track, detect where $\ln\vert\beta_2\vert$ crosses $\mu_2$:

- **Discrete crossing**: $d_0 \cdot d_1 < 0$ ($d_i = \ln\vert\beta_2\vert - \mu_2$). Each crossing is labeled with its jump direction: crossing upward through μ₂ → jump=+1, downward → jump=−1.
- **Continuum noise filtering**: if both sides of a crossing lie within the continuum band, it is numerical noise — skip.
- **Lazy refinement**: skip fsolve during bisection iterations (`refine_crossings=False`), using linear interpolation $(\theta_1^{\text{approx}}, \theta_2^{\text{approx}})$ directly as zeros. Only perform one fsolve refinement after bisection converges to the final $\mu_2$.

**Dual continuum+crossing extraction**: When the bisection finds a continuum (`is_continuum=True`), it means the Ronkin minimum is at a continuum boundary — but the θ₁ circle is not entirely continuum. The subset collection step in `collect_GBZ_subsets` does two passes over the root tracks:

1. **Continuum intervals** → `LineSubset` objects: extended θ₁ ranges where $|\ln|\beta_2| - \mu_2| < \mathrm{continuum\_tol}$ spans ≥ `min_continuum_pts` consecutive mesh points.
2. **Discrete crossings in gaps** → `PointSubset` objects: isolated sign changes detected in the non-continuum θ₁ regions, filtered to exclude crossings whose θ₁ falls inside any continuum interval, then refined via `fsolve` with the exact Jacobian.

### 2.3 Analytical Derivative of Average Winding w.r.t. μ₂

After refinement, for each zero $(\theta_1, \theta_2)$ solve a $2\times2$ linear system to obtain $\dot{\theta}_1 = d\theta_1/d\mu_2$:

$$\frac{dW}{d\mu_2} = \frac{1}{2\pi} \sum_i \text{jump}_i \cdot \dot{\theta}_{1,i}$$

This is used for a single Newton correction step on $\mu_2$ after refinement, avoiding re-bisection.

### 2.4 Unified Winding Number Computation

`_get_average_winding_from_zeros` handles both a1 and a2 directions uniformly:

```
Given zero list [(θ₁,θ₂), ...] and direction d ∈ {1,2}:
  1. Extract angular coordinates for direction d → sort → N segments
  2. For each segment midpoint:
     - Fix |β_param| = exp(μ_param), solve for the other variable's roots
     - Count roots with |β_solve| < exp(μ_count)
     - u_d = count_below - M_d
  3. average = Σ(u_d × segment_width) / (2π)
```

The a2 direction has an independent Hungarian matching pipeline; a1 reuses the a2 zero points (different directions use different angular coordinates for partitioning).

Special case — no zeros: the full circle is a single segment; u is evaluated at one midpoint (θ=0 for the appropriate variable). The return includes `non_zero_area`, a lightweight plateau indicator: when the winding is zero but concentrated in a tiny angular region, `non_zero_area` is small.

### 2.5 Bisection Solvers

**Inner loop** (μ₂ bisection, `_find_mu2_for_w2_zero`):

```
Given (E, μ₁), bisect μ₂ so that w2 winding = 0:
  1. Pre-compute _compute_root_tracks (one Hungarian pass) → cached tracks dict
  2. Adaptive range expansion: fast evaluation with unrefined winding (refine_crossings=False)
  3. Bisection: one _compute_winding_from_tracks per step (no fsolve)
  4. Continuum handling inside bisection: μ₂ ± ε perturbation
     - Opposite signs → this μ₂ is the w2=0 boundary → return is_continuum=True
     - Same sign → use perturbed winding to guide bisection
  5. Post-convergence refinement: refine_crossings=True → refined zeros + dW/dμ₂
  6. Newton correction (_refine_and_correct): μ₂ ← μ₂ − W_ref / (dW/dμ₂)
```

Key design: root tracks are μ₂-independent — computed once per `(E, μ₁)` and reused across ~30 μ₂ evaluations.

**Outer loop** (μ₁ bisection, `bisect_amoeba_ronkin_min`):

```
Given E, bisect μ₁ so that w1 winding = 0 (while maintaining w2 = 0):
  1. At each μ₁ endpoint, run inner bisection → (μ₂, zeros)
  2. Evaluate w1 from zeros at each endpoint — adaptive range expansion until w1_low * w1_high ≤ 0
  3. Bisection on μ₁:
     a. At μ₁_mid, inner bisection finds μ₂ with w2 = 0 → (μ₂, zeros, is_continuum)
     b. Continuum path (is_continuum=True):
        - _resolve_continuum: perturb μ₂ ± ε → w2 limits, perturb μ₁ ± ε → w1 limits
        - w2_opposite and w1_opposite → Ronkin minimum (continuum boundary)
        - w2_opposite but not w1_opposite → use w1 sign for bracket update
     c. Normal path (is_continuum=False):
        - Compute w1 from zeros via _get_average_winding_from_zeros (direction=1)
     d. |w1_mid| < xtol or bracket < xtol → converged
  4. Return {"mu1", "mu2", "zeros", "is_continuum"} plus internal fields
```

### 2.6 Plateau Detection

At some energies (particularly near zero-plateau boundaries in next-nearest-neighbor models), the bisection may converge to a false Ronkin minimum where the winding changes sign over a vanishingly narrow angular region. The plateau check in `collect_GBZ_subsets` filters these out:

**Pre-check** (3 lightweight conditions, all must be satisfied to trigger the expensive probe):
1. `w1_area < threshold`: the a1 non-zero winding interval is tiny (most of θ₂ has u₁ = 0)
2. `w2_area < threshold`: same for a2 (most of θ₁ has u₂ = 0)
3. `_check_zeros_are_clustered`: every zero has another zero within `threshold × 2π` distance on the (θ₁, θ₂)-torus

**Probe** (`_probe_zero_plateau_near_mu1`): step away from the candidate μ₁ in both directions; at each step, re-run the inner μ₂ bisection and check for the definitive plateau signature — w1 ≈ 0 and zero w2 crossings. If found, the result is classified as non-amoeba (empty subsets).

This check is applied only in the non-continuum case. Continuum results skip plateau detection since continuum bands are genuinely part of the GBZ.

## 3. API Reference

### 3.1 `get_hungarian_sorted_roots`

```python
def get_hungarian_sorted_roots(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1: float,
    N_points: int = 301,
) -> tuple[np.ndarray, np.ndarray, int, int]:
```

**Returns**: `(theta1_arr, tracked_roots, M, N)`
- `theta1_arr`: (N_points,) θ₁ mesh
- `tracked_roots`: (N_points, M+N) continuous root tracks
- `M`: denominator order in β₂
- `N`: max numerator degree in β₂ − M

### 3.2 `get_a2_average_winding`

```python
def get_a2_average_winding(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1: float,
    mu2: float,
    N_points: int = 301,
    continuum_tol: float = 1e-8,
    min_continuum_pts: int = 3,
) -> float:
```

**Returns**: $\partial R / \partial \mu_2$, the a2-direction average winding number.

Continuum handling: when a root track stays at $|\beta_2| \approx \exp(\mu_2)$ over an extended θ₁ range, spurious crossings from numerical noise are filtered out. Callers should perturb μ₂ to resolve the ambiguity.

### 3.3 `get_a1_average_winding`

```python
def get_a1_average_winding(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1: float,
    mu2: float,
    N_points: int = 301,
    continuum_tol: float = 1e-8,
    min_continuum_pts: int = 3,
) -> float:
```

**Returns**: $\partial R / \partial \mu_1$. Reuses a2 zero points (no independent Hungarian matching needed).

### 3.4 `bisect_amoeba_ronkin_min`

```python
def bisect_amoeba_ronkin_min(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1_low: float = -1,
    mu1_high: float = 1,
    mu2_low: float = -1,
    mu2_high: float = 1,
    N_points: int = 301,
    continuum_tol: float = 1e-8,
    min_continuum_pts: int = 3,
    continuum_perturb: float = 1e-4,
    max_iter: int = 60,
    xtol: float = 1e-10,
    max_range_expansions: int = 10,
    range_expand_factor: float = 2.0,
) -> dict:
```

**Returns**: `{"mu1", "mu2", "zeros", "is_continuum", "_mu1_bracket", "_w1_bracket", "_exit_reason", "_w1_area", "_tracks"}`

Find the Ronkin function minimum — the $(\mu_1, \mu_2)$ satisfying w1 = w2 = 0 simultaneously. Outer bisection on μ₁, inner bisection on μ₂, with adaptive range expansion and continuum handling.

Key return fields:
- `mu1`, `mu2`: critical point where w1 = w2 = 0
- `zeros`: list of `(θ₁, θ₂)` crossing pairs (empty `[]` when `is_continuum=True`)
- `is_continuum`: whether the result is a continuum boundary (vs. discrete zeros)
- `_tracks`: cached root tracks dict (passed through to avoid re-computation in `collect_GBZ_subsets`)
- `_w1_area`: normalized non-zero interval area for w1 (plateau pre-check)

### 3.5 `collect_GBZ_subsets`

```python
def collect_GBZ_subsets(
    coeffs: np.ndarray,
    degs: np.ndarray,
    E_ref: complex,
    perc: float,
    debug_mode: bool = False,
    **options,
) -> GBZResult:
```

**Returns**: `GBZResult` with `subsets` (list of `PointSubset | LineSubset`), `index=(n_0d, n_1d)`, and `is_gbz` / `is_empty` properties.

The main entry point. Runs `bisect_amoeba_ronkin_min`, then converts the result into structured GBZ subsets:

- **Continuum case** (`is_continuum=True`): extracts both continuum intervals (`LineSubset`) and discrete crossings in the gaps (`PointSubset`).
- **Discrete case** (`is_continuum=False`): converts zero pairs into `PointSubset` objects.
- **Plateau check**: optionally probes for false Ronkin minima (zero-plateau boundaries), returning empty subsets when detected.

Key options (passed through `**options`):
| Option | Default | Description |
|--------|---------|-------------|
| `plateau_check` | `True` | Enable zero-plateau detection |
| `plateau_winding_tol` | `None` | w1 tolerance for plateau probe (default: `max(xtol, 1e-10)`) |
| `plateau_probe_radius` | `None` | Step size for plateau probe (default: `continuum_perturb`) |
| `plateau_area_threshold` | `1e-2` | Threshold for pre-check winding area and zero clustering |
| `N_points` | `301` | θ₁ mesh resolution |
| `continuum_tol` | `1e-8` | Tolerance for continuum band detection |
| `min_continuum_pts` | `3` | Minimum mesh points for a valid continuum interval |
| `continuum_perturb` | `1e-4` | Perturbation step for continuum resolution |
| `max_iter` | `60` | Maximum bisection iterations |
| `xtol` | `1e-10` | Convergence tolerance |
| `mu1_low`, `mu1_high` | `-1, 1` | Initial μ₁ bracket |
| `mu2_low`, `mu2_high` | `-1, 1` | Initial μ₂ bracket |

## 4. CLaurent Polynomial Format

The characteristic polynomial $f(E, \beta_1, \beta_2)$ is encoded using triplets:

```
[E_exponent, beta1_exponent, beta2_exponent]
```

**Example**: square lattice nearest-neighbor $f = E - \beta_1 - \beta_1^{-1} - \beta_2 - \beta_2^{-1}$

```python
coeffs, degs = [], []
coeffs.append(1.0);   degs.extend([ 1,  0,  0])   # +E
coeffs.append(-1.0);  degs.extend([ 0,  1,  0])   # -beta1
coeffs.append(-1.0);  degs.extend([ 0, -1,  0])   # -beta1^{-1}
coeffs.append(-1.0);  degs.extend([ 0,  0,  1])   # -beta2
coeffs.append(-1.0);  degs.extend([ 0,  0, -1])   # -beta2^{-1}
```

## 5. Internal Functions

### Root tracking

| Function | File | Purpose |
|----------|------|---------|
| `get_hungarian_sorted_roots` | `tracks.py` | Hungarian-matched root tracks across θ₁ mesh → `(θ₁_arr, tracked, M, N)` |
| `_compute_root_tracks` | `tracks.py` | Extends `get_hungarian_sorted_roots` with periodic boundary handling, pre-computed `ln|β₂|`, and `PolyDiffContext` → tracks dict |

### Winding computation

| Function | File | Purpose |
|----------|------|---------|
| `_find_exact_crossing` | `ronkin_winding.py` | `fsolve` + analytical Jacobian to refine a single (θ₁, θ₂) crossing |
| `_compute_zero_dtheta1_dmu2` | `ronkin_winding.py` | Solve dθ₁/dμ₂ at a refined zero (2×2 Cramer's rule) |
| `_get_average_winding_from_zeros` | `ronkin_winding.py` | Unified winding: partitions angular circle by zeros, sums u×width/(2π) |
| `_compute_winding_from_tracks` | `ronkin_winding.py` | Detect crossings + continuum + w2 + dW/dμ₂ from pre-computed tracks |
| `_compute_crossings_and_winding` | `ronkin_winding.py` | Thin wrapper: `_compute_root_tracks` → `_compute_winding_from_tracks` |
| `get_a2_average_winding` | `ronkin_winding.py` | Public API for w2 (μ₂-derivative of Ronkin function) |
| `get_a1_average_winding` | `ronkin_winding.py` | Public API for w1 (μ₁-derivative, reuses a2 zeros) |

### Bisection

| Function | File | Purpose |
|----------|------|---------|
| `_refine_and_correct` | `bisect.py` | Post-bisection refinement: fsolve zeros + analytical Newton correction on μ₂ |
| `_find_mu2_for_w2_zero` | `bisect.py` | Inner μ₂ bisection: adaptive range expansion, lazy refinement, continuum handling |
| `_resolve_continuum` | `bisect.py` | Resolve continuum point: w2 limits (μ₂ ± ε) and w1 limits (μ₁ ± ε, re-run inner bisection) |
| `bisect_amoeba_ronkin_min` | `bisect.py` | Outer μ₁ bisection for Ronkin minimum (w1 = w2 = 0) |

### Subset collection and plateau detection

| Function | File | Purpose |
|----------|------|---------|
| `collect_GBZ_subsets` | `amoeba.py` | Main entry point: bisection → GBZ subsets (PointSubset + LineSubset) |
| `_extract_continuum_intervals` | `amoeba.py` | Extract continuum θ₁ intervals from root tracks for LineSubset creation |
| `_detect_crossings_outside_continuum` | `amoeba.py` | Detect and refine discrete crossings in non-continuum θ₁ gaps for PointSubset creation |
| `_check_zeros_are_clustered` | `amoeba.py` | Plateau pre-check: test whether all zeros are near-degenerate pairs |
| `_is_zero_plateau_probe` | `amoeba.py` | Plateau probe criterion: w1≈0, no w2 crossings, not continuum |
| `_probe_zero_plateau_near_mu1` | `amoeba.py` | Walk away from candidate μ₁, re-run inner bisection at each step |

### Shared utilities (from `gbz_types`)

| Function | Purpose |
|----------|---------|
| `get_minor_degrees` | Extract (M, N) for a given direction from the characteristic polynomial |
| `find_cyclic_true_intervals` | Find contiguous True intervals in a cyclic boolean array |
| `sort_by_root_abs` | Sort complex roots by modulus |
| `hungarian_match_indices` | Hungarian matching (linear sum assignment) with chordal distance cost |
| `generate_probe_steps` | Generate sorted probe step sizes for plateau detection |

## 6. Package Structure

```
brute_force_amoeba/
├── __init__.py            # Public API exports (8 symbols)
├── amoeba.py              # Main entry: collect_GBZ_subsets, plateau detection, continuum+crossing extraction
├── bisect.py              # μ₂ and μ₁ bisection solvers, continuum resolution
├── ronkin_winding.py      # Winding number computation, zero-crossing detection, fsolve refinement
└── tracks.py              # Hungarian-matched root tracking across θ₁ mesh
```

## 7. Demo

Run `demos/demo_unified.py`:

```bash
python demos/demo_unified.py
```

Demo contents:
1. 2D Hatano-Nelson model: compare SGBZ and Amoeba results at E = 1.0
2. Outside-spectrum test: E = 5.0 returns `is_gbz=False`
3. Triplet conversion: `as_triplet()` transforms `(E, β₁, β₂)` into crystal momenta `(E, k₁, k₂)`
4. LineSubset lazy loading: `fill_beta2()` on-demand for continuum intervals

Additional demos:
- `demos/replication-ZWang.py` — parallel E-mesh sweep with multiprocessing
- `demos/Haldane-model-gainloss.py` — non-Hermitian Haldane model with gain/loss
- `demos/imaginary-degeneracy-splitting.py` — next-nearest-neighbor model with plateau detection
