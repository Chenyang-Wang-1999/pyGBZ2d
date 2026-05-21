# brute_force_solver — SGBZ Spectrum Calculation

Non-Hermitian spectrum computation based on the SGBZ (Strip Generalized Brillouin Zone) formulation.

## 1. Theoretical Background

### 1.1 SGBZ and Strip Winding Number

SGBZ integrates along the $\theta_1$ direction, with the base manifold $\mu_2 = \rho_{2,0}(\theta_1; E, \mu_1)$ (the PMGBZ curve determined by $|\beta_M| = |\beta_{M+1}|$). The strip winding number $W(E, \mu_1)$ is the winding number around the $\theta_2$ direction integrated over $\theta_1$.

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

### 2.1 `PolyDiffContext` (`winding.py`)

A cached context bundling a polynomial with its precomputed partial derivatives, separated from `WindingFun`:

```python
class PolyDiffContext:
    char_poly: pt.CLaurent          # characteristic polynomial
    dchar_poly: list[pt.CLaurent]   # partial derivatives w.r.t. each variable (precomputed cache)

    def eval_val(var) -> complex             # evaluate polynomial
    def eval_partials(var) -> list[complex]  # evaluate partial derivatives
    def eval_dmu2(var) -> tuple[float, float]  # (dμ₂/dμ₁, dμ₂/dθ₁)
```

### 2.2 `WindingFun` / `MatWindingFun` (`winding.py`)

Winding number computation shells:
- `WindingFun`: polynomial winding, parametrized loop + analytical derivatives
- `MatWindingFun`: matrix winding (sparse solver), for large matrix models

### 2.3 `process_interval()` (`strip_winding_number.py`)

Unified recursive interval processor (core of this refactoring round):

```
process_interval(θ_left, roots_left, θ_right, roots_right):
  1. Analyze boundary cluster matching confidence (_analyze_boundary_matching)
     - Evaluate on the local cost matrix of the M-1/M boundary equimodular cluster
     - Compute cost difference between "current match" and "boundary-swapped alternative"
     - confidence = min_exchange_margin / cost_scale
  2. If not confident → subdivide interval (recurse)
  3. If confident and no cross-boundary exchange → interval complete
  4. If confident with cross-boundary exchange → perform PMGBZ search
     - refine_accidental_point() bisection to locate
     - _validate_pmgbz_point() validate and classify positive/negative/zero charge
  5. Recurse on left and right sub-intervals
```

### 2.4 Riemann Sphere Chordal Distance Matching (`strip_winding_number.py`)

Hungarian matching uses chordal distance on the Riemann sphere as the cost function (avoids $0/\infty$ false swaps in the complex plane):

```python
def _chordal_cost_matrix(roots_from, roots_to):
    # Project complex roots to the Riemann sphere (R³)
    # Use Euclidean distance in R³ as matching cost
```

### 2.5 Continuum Handling

When $|\beta_M| = |\beta_{M+1}|$ holds over a continuous $\theta_1$ interval:
- Strip winding returns a $(W_{\text{left}}, W_{\text{right}})$ tuple
- Left/right limits are resolved by perturbing $\mu_1$ by $\pm\varepsilon$

## 3. API Reference

### 3.1 `calculate_point_roots`

```python
def calculate_point_roots(
    char_poly: pt.CLaurent,
    param_ind_ctype: pt.CIndexVec,   # parameter variable indices (0,1) = (E,β₁)
    param_val: np.ndarray,           # parameter values
    var_ind_ctype: pt.CIndexVec,     # variable to solve for [2] = β₂
    M_max: int,                      # denominator order
    N_max: int,                      # max numerator degree - M
) -> list[complex]:
```

Fix $(E, \beta_1)$ and solve for all $\beta_2$ roots of the resulting univariate Laurent polynomial.

### 3.2 `get_roots_and_PMGBZ`

```python
def get_roots_and_PMGBZ(
    poly_diff: PolyDiffContext,
    E_ref: complex,
    mu1: float,
    N_points: int = 301,
    zero_tol: float = 1e-10,
    GBZ_check_tol: float = 1e-6,
) -> tuple[list[dict], np.ndarray, np.ndarray, dict]:
```

**Returns**: `(PMGBZ_points, theta1_arr, sols_arr, info)`
- `PMGBZ_points`: list of PMGBZ points, each containing `theta1`, `beta2_sols`, `is_continuum`
- `theta1_arr`: extended θ₁ grid (with periodic closure)
- `sols_arr`: corresponding sorted root array
- `info`: `{"M", "N", "continuum_flag"}`

### 3.3 `get_strip_winding`

```python
def get_strip_winding(
    poly_diff: PolyDiffContext,
    E_ref: complex,
    mu1: float,
    N_points: int = 301,
    with_gap_info: bool = False,
    continuum_perturb: float = 1e-2,
    zero_tol: float = 1e-10,
    GBZ_check_tol: float = 1e-6,
) -> tuple:
```

**Returns**: `(W_strip, PMGBZ_points)` or `(W_strip, PMGBZ_points, gap_info)`

When continuum exists, $W_{\text{strip}}$ is a $(W_{\text{left}}, W_{\text{right}})$ tuple.

### 3.4 `get_loop_winding`

```python
def get_loop_winding(
    poly_diff: PolyDiffContext,
    E_ref: complex,
    mu1: float,
    mu2_fun: callable,    # μ₂(θ₁) interpolation function
    theta2: float,        # fixed θ₂
    N_seg: int = 5,
) -> float:
```

Single-loop winding number (fixed θ₂, along the PMGBZ curve).

### 3.5 `SGBZSolver`

```python
class SGBZSolver:
    def __init__(self, char_poly: pt.CLaurent)
    def solve_for_E(
        self, E_ref: complex,
        mu1_guess: tuple[float, float] = (-1, 1),
        zero_tol: float = 1e-10,
        N_points: int = 101,
    ) -> tuple[float, list[dict]]:
```

Bisect μ₁ so that strip winding = 0. Returns (μ₁, PMGBZ_points).

### 3.6 `check_SGBZ`

```python
def check_SGBZ(
    coeffs: np.ndarray, degs: np.ndarray,
    E_ref: complex, perc: float,
    debug_mode: bool = False,
) -> dict:
```

Batch checking entry point — given polynomial coefficients, energy, and progress percentage, returns SGBZ determination results.

### 3.7 Other Exports

| Function | Purpose |
|----------|---------|
| `complex_root` | Complex equation root finding (with analytical Jacobian) |
| `ComplexEqConverter` | Complex → real equation Jacobian conversion |
| `poly_to_np_coefficients` | Laurent polynomial → numpy coefficient array |
| `get_winding_number` | Numerical integration for winding number |
| `get_minor_degrees` | Extract (M, N) |

## 5. Key Numerical Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `zero_tol` | 1e-10 | PMGBZ gap zero-threshold |
| `GBZ_check_tol` | 1e-6 | Equimodular cluster boundary detection tolerance |
| `match_confidence_tol` | 1e-3 | Hungarian matching confidence relative margin |
| `continuum_perturb` | 1e-2 | μ₁ perturbation amount for continuum |
| `double_root_tol` | `GBZ_check_tol` | Double root detection tolerance |

`match_confidence_tol` is an empirical threshold and may need tuning based on actual model scan results.

## 6. Relation to `brute_force_amoeba`

| Aspect | brute_force_solver (SGBZ) | brute_force_amoeba |
|--------|---------------------------|---------------------|
| Base manifold | $\mu_2 = \rho_{2,0}(\theta_1)$ variable curve | $\mu_2 = \text{const}$ level surface |
| Root ordering | By $\vert\beta_2\vert$ | Hungarian matching continuous tracking |
| Integration direction | $\theta_1$ only | Both $\theta_1$ and $\theta_2$ |
| Continuum criterion | $\vert\beta_M\vert \approx \vert\beta_{M+1}\vert$ | $\vert\ln\vert\beta_2\vert - \mu_2\vert \approx 0$ |
| Matching confidence | Exchange cost margin analysis | Not needed (continuous tracking) |
| Spectral inclusion | $\sigma_{\text{SGBZ}} \subseteq \sigma_{\text{Amoeba}}$ | — |
