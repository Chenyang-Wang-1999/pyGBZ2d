# pygbz2d.amoeba — Amoebic Spectrum Calculation

Non-Hermitian spectrum computation based on the amoeba formulation. Complementary to the SGBZ formulation (`pygbz2d/sgbz/`).

The zero-solving backend is `continuation.ZeroManager` (wrapped as `AmoebaZeroManager`); subset extraction and winding computation live in `zm_extract.py`. There is no separate `tracks.py` module — root tracking via Hungarian matching is now handled inside `ZeroManager`.

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
| Root tracking | Fixed $\theta_1$ mesh, Hungarian matching across samples | `continuation.ZeroManager` — adaptive $\beta_2$-root tracks that correctly cross multiple roots, with Hungarian matching handled inside the ZM |
| Continuum detection | $\vert\beta_M\vert \approx \vert\beta_{M+1}\vert$ | Per-segment fraction of samples with $\vert\ln\vert\beta_2\vert - \mu_2\vert < \text{tol}$ exceeding `frac` |
| Spectral inclusion | $\sigma_{\text{Amoeba}} \supset \cup_j \sigma_{\text{SGBZ},j}$ | — |

## 2. Algorithm Architecture

### 2.1 Zero-Manager Backend

The root-tracking layer is `continuation.ZeroManager` (see `continuation/`), wrapped as `AmoebaZeroManager(ZeroManager)` in `zm_extract.py`. Unlike a fixed uniform $\theta_1$ mesh, the ZM produces adaptive $\beta_2$-root tracks that correctly cross through multiple-root (MR) points and segments.

`AmoebaZeroManager` adds one cached field on top of `ZeroManager`:

- `seg_logabs`: `list[np.ndarray]`, where `seg_logabs[s]` has shape `(N_s, K)` = `np.log(np.abs(segments[s].tracked_roots))`.

This caching is central to the amoeba pipeline: the same root tracks are evaluated at many $\mu_2$ values (subset extraction plus ~30 $\mu_2$-bisection steps), and the cached logs let the extractor and winding code index a precomputed array rather than recomputing logs on every call. The ZM is $\mu_2$-independent — built once per `(E, μ₁)` and reused across all $\mu_2$ evaluations.

The periodic boundary ($\theta_1 = 0 \to 2\pi$) is handled inside the ZM via `boundary_perm`, which the extractor folds back through `boundary_perm_inv` (see `_zero_identity_key` in §5).

### 2.2 Continuum Detection and Mask

A track is declared **continuum** for a segment when a fraction `> frac` (default 0.9) of its samples satisfy $|\ln|\beta_2| - \mu_2| < \text{tol}$ (default `CONTINUUM_TOL = 1e-6`). This is computed by `_continuum_mask` (zm_extract.py:93), which returns a `list[np.ndarray]` of shape `(K,)` per segment — `mask[s][j]` is True when track `j` of segment `s` is in the continuum band.

A genuine continuum is constant-modulus to ~1e-9; a transversal crossing leaves the level after one sample. The 1e-6 tolerance separates them.

### 2.3 Crossing Detection and Lazy Refinement

For each non-continuum track, the winding/extraction code detects two kinds of crossings of $\ln|\beta_2| = \mu_2$:

- **Sign-change crossing** (`'cross'`): $d_i \cdot d_{i+1} < 0$ where $d_i = \ln|\beta_2|_i - \mu_2$. Labeled with a jump direction: $d_i < 0 \Rightarrow \text{jump}=+1$ (upward through $\mu_2$), else $-1$.
- **Exact touch** (`'zero'`): $d_i == 0$ at a sample. Recorded separately because the sign-change product cannot trigger when either endpoint is exactly zero.

**Lazy refinement**: during bisection iterations `refine=False`, so crossings use linear interpolation $(\theta_1^{\text{approx}}, \theta_2^{\text{approx}})$ directly. Only at the final $\mu_2$ (or for the `extract_amoeba_subsets` `'solve'` mode) is `_find_exact_crossing` invoked — a `scipy.optimize.fsolve` solve of $f(E, \beta_1, \beta_2) = 0$ with the exact $2\times2$ Jacobian.

### 2.4 Analytical Derivative of Average Winding w.r.t. μ₂

After refinement, for each zero $(\theta_1, \theta_2)$ `_compute_zero_dtheta1_dmu2` solves a $2\times2$ real linear system (Cramer's rule) to obtain $\dot{\theta}_1 = d\theta_1/d\mu_2$:

$$\frac{dW}{d\mu_2} = \frac{1}{2\pi} \sum_i \text{jump}_i \cdot \dot{\theta}_{1,i}$$

This is computed in `amoeba_windings` (zm_extract.py:561) only when `refine=True`. It is used by `_refine_and_correct` (bisect.py:24) for a single Newton correction step on $\mu_2$ after refinement, avoiding a re-bisection:

$$\mu_2 \leftarrow \mu_2 - W_{\text{ref}} / (dW/d\mu_2).$$

### 2.5 Unified Winding Number Computation

`_get_average_winding_from_zeros` (ronkin_winding.py:80) handles both a1 and a2 directions uniformly. It does not itself build root tracks — it consumes the zero list produced by `amoeba_windings` or `extract_amoeba_subsets`:

```
Given zero list [(θ₁,θ₂), ...] and direction d ∈ {1,2}:
  1. Partition the angular circle in the direction-d variable by the zero
     coordinates (θ₁ for d=2, θ₂ for d=1) → sorted cut points → N segments.
  2. For each segment midpoint:
     - Fix |β_param| = exp(μ_param), solve for the other variable's roots
     - count_below = #roots with ln|β_solve| < μ_count
     - u_d = count_below − M
  3. average        = Σ(u_d × segment_width) / (2π)
     non_zero_area = Σ(|u_d| × segment_width) / (2π)
```

- **direction=2** (w2, $\partial R/\partial\mu_2$): partition by $\theta_1$, fix $\beta_1 = e^{\mu_1}$, solve $\beta_2$, count roots with $\ln|\beta_2| < \mu_2$.
- **direction=1** (w1, $\partial R/\partial\mu_1$): partition by $\theta_2$, fix $\beta_2 = e^{\mu_2}$, solve $\beta_1$, count roots with $\ln|\beta_1| < \mu_1$.

**Special case — no zeros**: the full circle is a single segment of width $2\pi$; $u$ is evaluated at one midpoint ($\theta=0$ for the appropriate variable). The returned `non_zero_area` equals $|u|$.

`non_zero_area` is a lightweight plateau indicator: at a zero-plateau boundary both w1 and w2 non-zero areas are tiny (most of the circle has $u = 0$). Callers skip the expensive plateau probe when either area exceeds a safe threshold (default `1e-2`).

### 2.6 Bisection Solvers

**Inner loop** ($\mu_2$ bisection, `_find_mu2_for_w2_zero`, bisect.py:100):

```
Given (E, μ₁), bisect μ₂ so that w2 average winding = 0:
  1. Build one AmoebaZeroManager at (E, μ₁) (or accept a pre-built _zm)
     → reused across ~30 μ₂ evaluations (μ₂-independent).
  2. Adaptive range expansion with unrefined winding (refine=False):
     - If an endpoint falls in a continuum (winding is None), expand and retry.
     - Stop when w_low * w_high ≤ 0.
  3. Bisection with unrefined winding (no fsolve per step):
     - If a midpoint hits a continuum, perturb μ₂ ± continuum_perturb:
       * w_left * w_right < 0 → this μ₂ is the w2=0 boundary → return
         is_continuum=True, zeros=[].
       * Same sign → use the perturbed winding to update the bracket.
  4. On convergence (|w_mid| < xtol or bracket < xtol), refine:
     - _refine_and_correct: refine crossings (refine=True) + analytical
       Newton correction on μ₂ via dW/dμ₂. Falls back to bracket midpoint
       when dW/dμ₂ ≈ 0; clamps the step to the initial bracket.
```

Returns a dict with keys `mu2`, `zeros`, `is_continuum`, `winding`, `_zm`. The `winding` field is `(w_left, w_right)` when `is_continuum=True` (the perturbed limits), else the refined winding scalar.

**Outer loop** ($\mu_1$ bisection, `bisect_amoeba_ronkin_min`, bisect.py:347):

```
Given E, bisect μ₁ so that w1 = 0 while maintaining w2 = 0:
  1. At each μ₁ endpoint, run the inner bisection → (μ₂, zeros, _zm).
     Evaluate w1 via _get_average_winding_from_zeros(direction=1) from the
     inner zeros. Adaptive range expansion until w1_low * w1_high ≤ 0.
  2. Bisection on μ₁:
     a. At μ₁_mid, build a fresh AmoebaZeroManager and run() it.
     b. Inner bisection at μ₁_mid → (μ₂_mid, zeros_mid, is_continuum).
     c. Continuum path (is_continuum=True):
        - _resolve_continuum: perturb μ₂ ± ε → w2 limits; perturb μ₁ ± ε,
          re-run inner bisection at each, compute w1 from the resulting
          (μ₂, zeros) → w1 limits.
        - w2_opposite AND w1_opposite → Ronkin minimum at a continuum
          boundary → return is_continuum=True, zeros=[].
        - w2_opposite but not w1_opposite → use w1 sign for bracket update.
        - Neither → RuntimeError (the inner b2 straddle check disagrees
          with the continuum resolution; the top-level collector turns
          this into a failed GBZResult rather than a fake success).
     d. Normal path (is_continuum=False):
        - Compute w1 (and w1_area) from zeros_mid via
          _get_average_winding_from_zeros(direction=1).
     e. |w1_mid| < xtol or bracket < xtol → converged.
  3. Return dict (see §3.4).
```

### 2.7 Subset Extraction

`extract_amoeba_subsets` (zm_extract.py:317) turns a solved `AmoebaZeroManager` at the critical `(μ₁, μ₂)` into structured GBZ subsets. Three modes: `'coarse'` (linear only), `'fine'`/`'solve'` (fsolve refinement, falling back to linear on failure). `collect_GBZ_subsets` uses `'solve'`.

**Per-segment collection** (before MR joining):
- Continuum tracks (mask True) → one `_LinePiece` per segment per track, carrying the full segment's `theta1_arr` and `beta2_arr` plus `ml`/`mr` = outermost segment indices spanned.
- Non-continuum tracks → record `hits` as `(seg_idx, i, j, kind)` with `kind ∈ {'zero', 'cross'}`, deferred to a second pass for boundary dedup.

**Continuum join across MRs** (`_join_continuum_across_mrs`, zm_extract.py:157): a segment boundary is an MR, but only the roots listed in `multiple_roots[mr].cluster_indices` are genuinely multiple there; every other root passes straight through. A continuum track whose endpoint root is **not** in the cluster continues into the adjacent segment on the matched track. The join:
- Matches continuation by **column identity**: interior MR rows share one track frame, while the cyclic first/last-segment seam is translated through `boundary_perm` (`roots_right[boundary_perm] == roots_left`).
- `_is_cluster_endpoint` (zm_extract.py:127) tests whether an endpoint track column is part of the MR cluster; `mr < 0` marks the $\theta_1=0/2\pi$ circle seam (no MR there). At the θ=2π reuse of boundary MR 0, the right-frame column is translated back through `boundary_perm` before checking MR 0's θ=0-frame `cluster_indices`.
- `_merge_two` (zm_extract.py:263) concatenates as `[prev, cur[1:]]` — drops the shared MR row, keeps $\theta_1$ monotonic. For the cyclic seam ($\theta_1=0 \equiv 2\pi$), `prev` is the last segment and `cur` is segment 0; the seam lands inside the array so the array's two ends fall on the real terminators.
- Iterated to a fixpoint so chains and the cyclic seam converge. Raises `ValueError` on a topology inconsistency (a track ends as a non-cluster root in one segment but the matched track is a cluster root in the other).

**PointSubset finalization** (after MR joining):
- `_finalize_crossing` (zm_extract.py:516) turns each hit into `(θ₁, β₂)` per `mode`.
- **Rule 1** (continuum-endpoint snap): a crossing within `snap_tol` (default `SNAP_TOL = 1e-3`) of a continuum LineSubset endpoint is dropped **only if** its `β₂` matches an in-band endpoint root within `ROOT_TOL = 1e-9` — i.e., the crossing is the line's own endpoint curve. A same-$\theta_1$ crossing on a different, unconnected root survives.
- **Rule 2** (exact-touch dedup): `d == 0` exact touches are deduplicated by zero identity, never by $\theta_1$. `_zero_identity_key` (zm_extract.py:467) returns `(endpoint, track)`:
  - Interior sample → unique key `('interior', id(seg), i, j)`, never deduped.
  - Shared interior MR endpoint → `('mr', R, j)` collapses both sides.
  - Circle boundary $\theta_1=0\equiv2\pi$ → `('circle', j_left)` (or `boundary_perm_inv[j]` on the right) folds the right end onto the left track.
- Sign-change `'cross'` hits are always unique.
- The `ml`/`mr` joining fields are dropped; the public GBZ output is plain `LineSubset`/`PointSubset`.

### 2.8 Winding for the μ₂-Bisection

`amoeba_windings` (zm_extract.py:561) is the winding source used by the $\mu_2$ bisection. It uses **RAW** crossing detection — every crossing participates in the winding integral, with no Rule 1/Rule 2 boundary filtering (those are subset-output concerns).

Returns `(winding, zeros, has_continuum, dW_dmu2)`:
- If any track is in the continuum band on any segment → `(None, None, True, None)`. The caller (`_find_mu2_for_w2_zero`) reacts by perturbing $\mu_2$.
- Otherwise `zeros` is a list of `(θ₁, θ₂, jump)` refined (or linear, if `refine=False`) crossing pairs, `winding` is computed by `_get_average_winding_from_zeros(direction=2)`, and `dW_dmu2` is the analytical derivative when `refine=True` (else `0.0`).

### 2.9 Plateau Detection

At some energies (particularly near zero-plateau boundaries in next-nearest-neighbor models), the bisection may converge to a false Ronkin minimum where the winding changes sign over a vanishingly narrow angular region. The plateau check in `collect_GBZ_subsets` filters these out.

**Pre-check** (three lightweight conditions, all must be satisfied to trigger the expensive probe; non-continuum case only):
1. `w1_area < plateau_area_threshold`: the a1 non-zero winding interval is tiny.
2. `w2_area < plateau_area_threshold`: same for a2.
3. `_check_zeros_are_clustered`: every zero has another zero within `plateau_area_threshold × 2π` distance on the $(\theta_1, \theta_2)$-torus (thin adapter over `pygbz2d.core.check_points_clustered_on_torus`, dropping the trailing `jump`).

**Probe** (`_probe_zero_plateau_near_mu1`, amoeba.py:59): shared `pygbz2d.core.probe_zero_plateau` driver. Steps away from the candidate $\mu_1$ in both directions; at each step, re-runs the inner $\mu_2$ bisection and tests the strict plateau signature — `success`, `not is_continuum`, `zero_count == 0`, and $|w1| \le$ `winding_tol`. If found, the result is classified as non-amoeba (empty subsets).

This check is applied only in the non-continuum case. Continuum results skip plateau detection since continuum bands are genuinely part of the GBZ.

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

**Returns**: `GBZResult` with `subsets` (list of `PointSubset | LineSubset`), `index=(n_0d, n_1d)`, and `is_gbz` / `is_empty` properties.

The main entry point. Builds a `CharPoly(coeffs, degs)`, runs `bisect_amoeba_ronkin_min`, then reuses the solved `_zm` (AmoebaZeroManager built at the converged $\mu_1$) for subset extraction via `extract_amoeba_subsets(..., mode='solve')` — no second ZM run. Plateau detection is optionally applied (non-continuum case). `debug_mode=True` re-raises exceptions instead of returning a failed `GBZResult`.

**Options** (passed through `**options` to `bisect_amoeba_ronkin_min`):

| Option | Default | Description |
|--------|---------|-------------|
| `mu1_low`, `mu1_high` | `-1, 1` | Initial $\mu_1$ bracket |
| `mu2_low`, `mu2_high` | `-1, 1` | Initial $\mu_2$ bracket |
| `N_points` | `301` | Accepted but **ignored** — the mesh is adaptive inside ZeroManager. Kept for API compatibility. |
| `continuum_tol` | `1e-6` | Tolerance for continuum band detection (`CONTINUUM_TOL`) |
| `min_continuum_pts` | `3` | Accepted but **ignored** — continuum is now fraction-based via `frac`. Kept for API compatibility. |
| `continuum_perturb` | `1e-4` | $\mu_2$ / $\mu_1$ perturbation step for continuum resolution |
| `max_iter` | `60` | Maximum bisection iterations (inner and outer) |
| `xtol` | `1e-10` | Convergence tolerance |
| `max_range_expansions` | `10` | Maximum adaptive range expansions |
| `range_expand_factor` | `2.0` | Multiplicative range expansion factor |
| `frac` | `CONTINUUM_FRAC` (0.9) | Continuum fraction threshold |
| `plateau_check` | `True` | Enable zero-plateau detection |
| `plateau_winding_tol` | `None` | w1 tolerance for plateau probe (default: `max(xtol, 1e-10)`) |
| `plateau_probe_radius` | `None` | Step size for plateau probe (default: `continuum_perturb`) |
| `plateau_area_threshold` | `1e-2` | Threshold for pre-check winding area and zero clustering |
| `zm_run_kwargs` | `{}` | Keyword arguments forwarded to `ZeroManager.run()` (e.g. `h0`, `ctrl`, `min_dtheta`, `cluster_tol`, `mr_jump`, `verbose`). Kept separate from the bisection options, which `ZeroManager.run` does not accept. |

### 3.2 `bisect_amoeba_ronkin_min`

```python
def bisect_amoeba_ronkin_min(
    char_poly: CharPoly,
    E_ref: complex,
    mu1_low: float = -1,
    mu1_high: float = 1,
    mu2_low: float = -1,
    mu2_high: float = 1,
    N_points: int = 301,
    continuum_tol: float = 1e-6,
    min_continuum_pts: int = 3,
    continuum_perturb: float = 1e-4,
    max_iter: int = 60,
    xtol: float = 1e-10,
    max_range_expansions: int = 10,
    range_expand_factor: float = 2.0,
    frac: float = CONTINUUM_FRAC,
    zm_run_kwargs: Optional[dict] = None,
) -> dict:
```

Find the Ronkin function minimum — the $(\mu_1, \mu_2)$ satisfying $w_1 = w_2 = 0$ simultaneously. Outer bisection on $\mu_1$, inner bisection on $\mu_2$, with adaptive range expansion and continuum handling.

**Returns** a dict with keys:

| Key | Type | Description |
|-----|------|-------------|
| `mu1` | `float` | Critical $\mu_1$ where $w_1 = 0$ |
| `mu2` | `float` | Critical $\mu_2$ where $w_2 = 0$ |
| `zeros` | `list[(θ₁, θ₂, jump)]` | Crossing pairs at the critical point; `[]` when `is_continuum=True` |
| `is_continuum` | `bool` | Whether the result is a continuum boundary (vs. discrete zeros) |
| `_mu1_bracket` | `(float, float)` | Final $\mu_1$ bisection bracket (consumed by the plateau probe) |
| `_w1_bracket` | `(float, float)` | Final $w_1$ bracket |
| `_exit_reason` | `str` | `"w1_zero"`, `"bracket_xtol"`, or `"continuum_boundary"` |
| `_w1_area` | `float` | Normalized non-zero interval area for $w_1$ (plateau pre-check). Present only in the non-continuum path. |
| `_zm` | `AmoebaZeroManager` | The ZM built at the solved $\mu_1$ — reused by `collect_GBZ_subsets` for subset extraction (no second run). |

Note: `_w2_area` is **not** returned; only `_w1_area` is. (The $w_2$ area is computed on demand in `collect_GBZ_subsets` only when the pre-check fires.)

### 3.3 `extract_amoeba_subsets`

```python
def extract_amoeba_subsets(
    zm: AmoebaZeroManager,
    poly: CharPoly,
    E: complex,
    mu1: float,
    mu2: float,
    *,
    mode: ExtractMode = 'solve',
    tol: float = CONTINUUM_TOL,
    frac: float = CONTINUUM_FRAC,
    snap_tol: float = SNAP_TOL,
) -> list:
```

Build GBZ subsets at $\ln|\beta_1|=\mu_1$, $\ln|\beta_2|=\mu_2$ from a solved ZM's tracks. `mode ∈ {'coarse', 'fine', 'solve'}`; `collect_GBZ_subsets` uses `'solve'`. Returns a `list[PointSubset | LineSubset]`. See §2.7 for the full algorithm.

### 3.4 `amoeba_windings`

```python
def amoeba_windings(
    zm: AmoebaZeroManager,
    poly: CharPoly,
    E: complex,
    mu1: float,
    mu2: float,
    *,
    tol: float = CONTINUUM_TOL,
    frac: float = CONTINUUM_FRAC,
    refine: bool = True,
) -> tuple[Optional[float], Optional[list], bool, Optional[float]]:
```

Returns `(winding, zeros, has_continuum, dW_dmu2)`:
- `winding`: w2 average winding number; `None` when a continuum is present.
- `zeros`: list of `(θ₁, θ₂, jump)`; `None` when a continuum is present.
- `has_continuum`: whether any track is in the $\mu_2$ band over a whole segment.
- `dW_dmu2`: analytical derivative $dW/d\mu_2$; `0.0` when `refine=False` or no crossings; `None` when a continuum is present.

Uses RAW crossing detection — no Rule 1/Rule 2 boundary filtering. See §2.8.

### 3.5 `AmoebaZeroManager`

```python
class AmoebaZeroManager(ZeroManager):
    seg_logabs: list[np.ndarray]

    def __init__(self, poly: CharPoly, E_ref: complex, mu1: float): ...
    def run(self, **kwargs) -> None: ...
```

`ZeroManager` specialization caching per-segment `ln|β₂|`. After `.run()`, `seg_logabs[s]` has shape `(N_s, K) = np.log(np.abs(segments[s].tracked_roots))`. `run()` forwards `**kwargs` to `ZeroManager.run` (e.g. `h0`, `ctrl`, `min_dtheta`, `cluster_tol`, `mr_jump`, `verbose` via `zm_run_kwargs`), then recomputes the cache unconditionally (re-run resets segments).

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

### Ronkin winding primitives (`ronkin_winding.py`)

| Function | Purpose |
|----------|---------|
| `_find_exact_crossing` | `fsolve` + analytical Jacobian to refine a single $(\theta_1, \theta_2)$ crossing; returns `None` on non-convergence (residual ≥ 1e-10) |
| `_get_average_winding_from_zeros` | Unified winding from a zero partition — partitions the angular circle in the direction-d variable, sums $u \times \text{width} / (2\pi)$. Returns `(avg_winding, non_zero_area)`. |
| `_compute_zero_dtheta1_dmu2` | Analytical $d\theta_1/d\mu_2$ at a refined zero (2×2 Cramer's rule). Returns `0.0` when the determinant is near-singular. |

### ZeroManager integration and extraction (`zm_extract.py`)

| Function | Purpose |
|----------|---------|
| `AmoebaZeroManager` | `ZeroManager` subclass caching `seg_logabs` per segment. |
| `_continuum_mask` | Per-segment boolean mask: track is continuum when fraction of in-band samples `> frac`. |
| `extract_amoeba_subsets` | Full GBZ subset output (LineSubset + PointSubset) with Rule 1/Rule 2 boundary dedup. |
| `_finalize_crossing` | Turn a recorded hit into `(θ₁, β₂)` per `mode` (linear / fsolve-refined). |
| `_zero_identity_key` | Identity key for exact-touch dedup: `(endpoint, track)`, never `θ₁`. |
| `_is_cluster_endpoint` | Whether an endpoint track column is part of the boundary MR cluster (column identity; `boundary_perm` translation at the cyclic boundary-MR seam). |
| `_join_continuum_across_mrs` | Join per-segment continuum LineSubsets passing through MRs as non-cluster roots; iterated to fixpoint. |
| `_merge_two` | Concatenate two `_LinePiece`s, dropping the shared MR row; cyclic-seam variant puts the seam inside the array. |
| `_LinePiece` | Per-segment continuum LineSubset being joined; carries `ml`/`mr` (outermost segment indices spanned). |
| `amoeba_windings` | RAW w2 winding + zeros + `dW/dμ₂` for the $\mu_2$ bisection. |
| `_crossing_thetas` | `(θ₁, θ₂)` for a winding crossing, linear or refined. |

### Bisection (`bisect.py`)

| Function | Purpose |
|----------|---------|
| `_find_mu2_for_w2_zero` | Inner $\mu_2$ bisection: builds (or accepts) the ZM, adaptive range expansion, lazy refinement, continuum handling. |
| `_refine_and_correct` | Post-bisection refinement: fsolve zeros + analytical Newton correction on $\mu_2$. |
| `_resolve_continuum` | Resolve a continuum point: w2 limits ($\mu_2 \pm \varepsilon$) and w1 limits ($\mu_1 \pm \varepsilon$, re-run inner bisection). |
| `bisect_amoeba_ronkin_min` | Outer $\mu_1$ bisection for the Ronkin minimum ($w_1 = w_2 = 0$). |

### Plateau detection (`amoeba.py`)

| Function | Purpose |
|----------|---------|
| `_check_zeros_are_clustered` | Thin adapter over `pygbz2d.core.check_points_clustered_on_torus` (drops the trailing `jump`). |
| `_is_zero_plateau_probe` | Plateau probe criterion: `success`, `not is_continuum`, `zero_count == 0`, $|w_1| \le$ tol. |
| `_probe_zero_plateau_near_mu1` | Walk away from candidate $\mu_1$, re-run inner bisection at each step; uses shared `pygbz2d.core.probe_zero_plateau`. |

### Shared utilities (from `pygbz2d.core`)

| Function | Purpose |
|----------|---------|
| `get_minor_degrees` | Extract `(M, N)` for a given direction from the characteristic polynomial |
| `find_cyclic_true_intervals` | Find contiguous True intervals in a cyclic boolean array |
| `sort_by_root_abs` | Sort complex roots by modulus |
| `hungarian_match_indices` | Hungarian matching (linear sum assignment) with chordal distance cost |
| `generate_probe_steps` | Generate sorted probe step sizes for plateau detection |
| `check_points_clustered_on_torus` | Plateau pre-check: test whether all points are near-degenerate pairs |
| `probe_zero_plateau` | Shared plateau-probe driver loop |

## 6. Package Structure

```
pygbz2d/amoeba/
├── __init__.py            # Public API exports (10 symbols)
├── amoeba.py              # 240 lines — collect_GBZ_subsets, plateau detection
├── bisect.py              # 541 lines — μ₂ and μ₁ bisection solvers, continuum resolution
├── ronkin_winding.py      # 188 lines — winding primitives (fsolve refinement, zero partition, dθ₁/dμ₂)
└── zm_extract.py          # 678 lines — AmoebaZeroManager, subset extraction, amoeba_windings
```

Public API (`__init__.py`): `CharPoly`, `bisect_amoeba_ronkin_min`, `AmoebaZeroManager`, `extract_amoeba_subsets`, `amoeba_windings`, `collect_GBZ_subsets`, `PointSubset`, `LineSubset`, `GBZResult`, `ConnectedSubset`.

## 7. Demo

Run `playground/demo_unified.py`:

```bash
python playground/demo_unified.py
```

Demo contents:
1. 2D Hatano-Nelson model: compare SGBZ and Amoeba results at E = 1.0
2. Outside-spectrum test: E = 5.0 returns `is_gbz=False`
3. Triplet conversion: `as_triplet()` transforms `(E, β₁, β₂)` into crystal momenta `(E, k₁, k₂)`
4. LineSubset lazy loading: `fill_beta2()` on-demand for continuum intervals

Additional demos:
- `playground/replication-ZWang.py` — parallel E-mesh sweep with multiprocessing
- `playground/Haldane-model-gainloss.py` — non-Hermitian Haldane model with gain/loss
- `playground/imaginary-degeneracy-splitting.py` — next-nearest-neighbor model with plateau detection
