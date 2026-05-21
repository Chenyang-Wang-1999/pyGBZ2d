# brute_force_amoeba — Amoebic Spectrum Calculation

Non-Hermitian spectrum computation based on the amoeba formulation. Complementary to the SGBZ formulation (`brute_force_solver/`).

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

### 2.2 Crossing Detection and Lazy Refinement

For each track, detect where $\ln\vert\beta_2\vert$ crosses $\mu_2$:

- **Discrete crossing**: $d_0 \cdot d_1 < 0$ ($d_i = \ln\vert\beta_2\vert - \mu_2$). Each crossing is labeled with its jump direction: crossing upward through μ₂ → jump=+1, downward → jump=−1.
- **Continuum noise filtering**: if both sides of a crossing lie within the continuum band, it is numerical noise — skip.
- **Lazy refinement**: skip fsolve during bisection iterations (`refine_crossings=False`), using linear interpolation $(\theta_1^{\text{approx}}, \theta_2^{\text{approx}})$ directly as zeros. Only perform one fsolve refinement after bisection converges to the final $\mu_2$.

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

### 2.5 Bisection Solvers

**Inner loop** ($\mu_2$ bisection, `_find_mu2_for_a2_zero`):

```
Given (E, μ₁), bisect μ₂ so that a2 winding = target:
  1. Pre-compute _compute_root_tracks (one Hungarian pass) → cached tracks dict
  2. Adaptive range expansion: fast evaluation with unrefined winding (refine_crossings=False)
  3. Bisection: one _compute_winding_from_tracks per step (no fsolve)
  4. Post-convergence refinement: refine_crossings=True → refined zeros + dW/dμ₂
  5. Newton correction: μ₂ ← μ₂₀ − W_ref / (dW/dμ₂)
  6. Continuum handling: μ₂ ± ε perturbation, opposite signs → stop, same sign → continue
```

**Outer loop** ($\mu_1$ bisection, `bisect_amoeba_ronkin_min`):

```
Given E, bisect μ₁ so that a1 winding = 0 (while maintaining a2 = 0):
  1. At μ₁_mid, inner bisection finds μ₂ with a2 = 0 → (μ₂, zeros)
  2. Directly compute a1 from zeros (reuse zeros, no extra Hungarian pass)
  3. a1 degenerate (zeros cannot partition θ₂):
     - First check if a1_mid is already zero (zero-plateau early exit)
     - Otherwise μ₁ ± ε perturbation, re-run inner bisection → compare a1 signs
     - Opposite signs → boundary reached, stop
     - Same sign → use sign to guide outer bisection direction
```

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

### 3.4 `bisect_a2_winding`

```python
def bisect_a2_winding(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1: float,
    mu2_low: float,
    mu2_high: float,
    target_winding: float = 0.0,
    N_points: int = 301,
    continuum_tol: float = 1e-8,
    min_continuum_pts: int = 3,
    continuum_perturb: float = 1e-4,
    max_iter: int = 60,
    xtol: float = 1e-10,
) -> dict:
```

**Returns**: `{"mu2", "zeros", "is_continuum", "winding", "success"}`

Bisect μ₂ so that the a2 winding crosses `target_winding`. No adaptive range expansion (caller guarantees opposite signs in the bracket).

### 3.5 `bisect_amoeba_ronkin_min`

```python
def bisect_amoeba_ronkin_min(
    char_poly: pt.CLaurent,
    E_ref: complex,
    mu1_low: float,
    mu1_high: float,
    mu2_low: float,
    mu2_high: float,
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

**Returns**: `{"mu1", "mu2", "zeros", "is_continuum", "success"}`

Find the Ronkin function minimum — the $(\mu_1, \mu_2)$ satisfying a1 = a2 = 0 simultaneously. Outer bisection on μ₁, inner bisection on μ₂, with adaptive range expansion and continuum handling.

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

| Function | Purpose |
|----------|---------|
| `_get_minor_degrees_for_direction` | Extract (M, N) for a given direction |
| `_find_cyclic_true_intervals` | Find contiguous True intervals in a cyclic boolean array |
| `_find_exact_crossing` | `fsolve` + analytical Jacobian to refine a single zero |
| `_compute_zero_dtheta1_dmu2` | Solve dθ₁/dμ₂ at a refined zero (2×2 linear system) |
| `_compute_root_tracks` | Hungarian matching + root trajectory tracking → tracks dict (μ₂-independent) |
| `_compute_winding_from_tracks` | Detect crossings + continuum + winding + dW/dμ₂ from tracks |
| `_compute_crossings_and_winding` | Thin wrapper (calls `_compute_root_tracks` + `_compute_winding_from_tracks`) |
| `_get_average_winding_from_zeros` | Unified winding function (supports both a1/a2 directions) |
| `_refine_and_correct` | Refine zeros + analytical derivative Newton correction on μ₂ |
| `_find_mu2_for_a2_zero` | Inner μ₂ bisection (adaptive range expansion + lazy refinement + Newton correction) |
| `_a1_is_degenerate` | Heuristic to determine whether a1 is degenerate |

## 6. Demo

Run `demos/demo_amoeba.py`:

```bash
python demos/demo_amoeba.py
```

Demo contents:
1. Winding sweep — scan μ₂ at fixed μ₁, compare a1/a2 winding numbers
2. Hungarian root tracking — visualize continuous root tracks
3. μ₂ bisection — find where a2 = 0
4. Ronkin minimum — outer bisection for (μ₁, μ₂) where a1 = a2 = 0
5. Next-nearest-neighbor model — Ronkin minimum under NNN coupling
