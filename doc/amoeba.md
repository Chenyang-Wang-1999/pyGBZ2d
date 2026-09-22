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
| Root tracking | Adaptive `ZeroManager` tracks + pairwise event refinement | `continuation.ZeroManager` adaptive $\beta_2$-root tracks |
| Continuum detection | ItemView cluster vote | `std(ln|\beta_2|) < tol` per track per segment |
| Subset output | `collect_GBZ_subsets` in `sgbz_solver.py` | `collect_GBZ_subsets` in `amoeba.py` |

## 2. Algorithm Architecture

### 2.1 Zero-Manager Backend

`AmoebaZeroManager(ZeroManager)` caches `seg_logabs[s] = log|segments[s].tracked_roots|` after `run()`. Its `insert_solution` override automatically refreshes the affected cache after mesh insertion. The ZM is $\mu_2$-independent — built once per `(E, \mu_1)` and reused across all $\mu_2$ evaluations, including the root rows added during crossing refinement.

Before the inner bisection, `_try_fast_mu2` refines bracketed local minima and
maxima on **every** root track, including extrema that do not improve the
global gap bounds. For $g_j(\theta_1)=\ln|\beta_{2,j}(\theta_1)|$, its derivative
is the stored `seg.tangents[:, j].real`. Opposite derivative signs bracket a
stationary point; Brent's method refines the zero of the **actual derivative**,
not the derivative of an interpolating curve. Every new evaluation solves the
polynomial, matches the full root row to the segment tracks, computes its
tangents, and retains the row through `insert_solution`. The converged point
is retained as well. Existing endpoint stationary points already belong to
the mesh, and the $2\pi$ endpoint is read in its own track frame.

These samples separate a local dip or peak into intervals usable by every
subsequent $\mu_2$ level. This avoids repeating an extremum search for each
bisection step. Angular brackets are stored independently of mutable row
indices. Refinement failures raise with the energy, $\mu_1$, segment, track,
and bracket; they are not silently skipped. Intervals next to a non-finite
MR tangent retain one real midpoint sample, and only finite derivative sign
brackets exposed by that sample are refined. The singular endpoint itself
is not passed to Brent's method.

**Current limitation:** the screen uses endpoint derivative signs. It does
**not** search for several interior extrema hidden in a single interval whose
endpoint derivatives have the same sign. In particular, a local maximum and
minimum inside such an interval can still hide multiple level crossings.
This case is intentionally deferred; the preprocessing does not certify that
all intervals are monotone. It also does not recursively resolve arbitrary
oscillations next to MR endpoints. Flat-track continuum detection remains a
separate step. Extrema are located to finite angular accuracy controlled by
`amoeba.bisect.EXTREMUM_INSERT_REL_TOL`.

### 2.2 Continuum Detection (pre-bisection)

`detect_continuum(zm, tol)` detects, for every segment and every track column, whether `std(ln|β₂|[:, j]) < tol`. Detected flat tracks are returned as groups:

```
[(mu2_c, [(seg_idx, col_idx), ...]), ...]
```

where nearby `mu2_c` values (closer than `tol`) are merged. The `(seg_idx, col_idx)` list is the continuum membership consumed later by subset assembly — it is recorded at detection time and passed through the bisection result.

### 2.3 Crossing Detection

`find_crossings(zm, mu1, mu2, avoided_segments=None, return_refined=False)` finds crossings of `ln|β₂| = μ₂` represented by the current mesh:

- Traverses every segment and every column.
- Skips `(seg_idx, col_idx)` entries in `avoided_segments`.
- First detects exact touches with `d[:-1] == 0`, then sign changes with `d[1:] * d[:-1] < 0`.
- Returns `list[(β₁, β₂)]`. Linear interpolation unless `return_refined=True`, in which case `_find_exact_crossing` solves `ln|β₂_j(θ₁)| - μ₂ = 0` with Brent's bracketed method. Each interior evaluation solves the polynomial, matches all roots to the segment tracks, and inserts the complete root row into the ZM. Later μ₂ levels reuse this refined mesh.
- Insertion changes row indices and may expose crossings on other tracks. The refined scan restarts on the updated mesh, reusing accepted crossing rows with their original crossing orientation. Samples always retain actual polynomial roots; projection to the requested torus happens only when returning a crossing. Brent's angular convergence is the stopping criterion, with no additional polynomial-residual or log-modulus acceptance threshold. Non-finite values and failed bracketing/convergence raise with energy, μ values, segment, track, and bracket context; there is no linear fallback in refined mode.

The inner solver prepares extrema before calling this function. A standalone
call on a newly tracked ZM does not perform that preprocessing automatically.
Endpoint-based crossing detection is subject to the limitation in §2.1.

### 2.4 Winding

`calculate_a2_average_winding(zm, mu1, mu2, avoided_segments=None, return_refined=False)` is a thin consumer of `find_crossings`. It converts crossings back to angular zeros and computes the a2 average winding with `_get_average_winding_from_zeros(direction=2)`.

The internal `_calculate_a2_winding_and_zeros` also returns the zero partition. The fine μ₂ bisection retains this partition and returns it on convergence without refining the same crossings again. Integer winding on each angular interval is still obtained by independent midpoint polynomial solves.

`_get_average_winding_from_zeros` partitions the transverse angular circle: `direction=2` partitions theta1 and solves beta2; `direction=1` partitions theta2 and solves beta1. At each interval midpoint, `u` is the number of roots strictly inside the requested radius minus the Laurent denominator order. It sums `u × width / (2π)` and also returns `sum(abs(u) × width)/(2π)` for the plateau pre-check. The latter is weighted by winding magnitude, so it need not be an area fraction bounded by one.

### 2.5 Inner μ₂ Solve (`_find_mu2_for_w2_zero`)

Continuum-first dispatcher:

1. Build/reuse `AmoebaZeroManager`.
2. `_try_fast_mu2` — sort roots at θ₁=0 by `|β₂|` to select the first M and last N columns. Refine the bracketed local minima and maxima of **all** columns with actual derivative solves (§2.1), then compute `A = max(ln|β₂|` over first M`)`, `B = min(ln|β₂|` over last N`)`. If `A < B`, return `μ₂ = (A+B)/2` with no crossings. Otherwise, the same enriched mesh is reused by continuum probes and both bisection stages.
3. `detect_continuum` — if flat tracks exist, probe each merged `μ₂_c` at `μ₂_c ± ε` (`ESCAPE_LADDER`). Opposite w2 signs → return the continuum boundary and the matched `(seg_idx, col_idx)` members. Equal signs → tighten the μ₂ bracket with the signed probe points.
4. `_bisect_mu2_discrete` — pure two-stage bisection on the tightened bracket:
   - coarse stage with `return_refined=False` and `BISECT_COARSE_XTOL`;
   - fine stage on the coarse bracket with `return_refined=True`, stopping when `abs(w2) < wtol`; `wtol` defaults to `core.WINDING_ZERO_TOL` (1e-8). The fine bracket is expanded outward when refined endpoint windings no longer straddle zero.

### 2.6 Outer μ₁ Bisection (`bisect_amoeba_ronkin_min`)

- Endpoint expansion and midpoint bisection on μ₁.
- Every inner return is checked for `is_continuum`.
- `_handle_continuum` resolves a continuum inner result: at `μ₁ ± ε` it builds a fresh ZM, prepares bracketed local extrema, runs `find_crossings` at the fixed `μ₂_c`, and computes the a1 average winding from those crossings. Opposite w1 signs end the outer bisection; equal signs update the μ₁ bracket.
- Discrete path computes w1 from the inner `zeros` and stops when `abs(w1) < wtol`, using the same winding tolerance as the inner solve.

### 2.7 Subset Assembly

`collect_GBZ_subsets` runs the bisection, then:

- **`is_continuum == False`**: each `(θ₁, θ₂)` zero becomes a `PointSubset`.
- **`is_continuum == True`**:
  1. materialize the recorded continuum `(seg_idx, col_idx)` members as per-segment `JoinableLinePiece`s;
  2. splice pieces that meet end-to-end (interior MRs share a track frame; the θ=0≡2π seam is translated through `boundary_perm`; MR-cluster endpoints are genuine terminators);
  3. run `find_crossings` with the continuum members passed as `avoided_segments` to find discrete crossings;
  4. convert those crossings to `PointSubset`s, dropping a point only when its circular θ₁ and chordal β₂ distances to the same continuum line sample are both within `SNAP_TOL` (with matching μ₁). Matching uses samples, not interpolation across tracks.

### 2.8 Plateau Detection

Non-continuum results still pass through the plateau pre-check/probe (tiny w1/w2 non-zero area + clustered zeros). Continuum results skip plateau detection.

Both sides of each mu1 probe step are evaluated. Either successful,
non-continuum probe with no crossing zeros and zero w1 suffices to report
a plateau. Exceptions from an inner mu2 solve are recorded as failed probes
and the ladder continues; they do not count as plateau evidence.

The plateau probe uses the same `wtol` for its inner μ₂ solve and zero-winding predicate. Its diagnostic result exposes this value as `wtol` (and the shared probe's `zero_tol`). There is no separate `winding_tol` or `plateau_winding_tol` parameter.

## 3. API Reference

### 3.1 `collect_GBZ_subsets`

```python
def collect_GBZ_subsets(
    coeffs: np.ndarray,
    degs: np.ndarray,
    E_ref: complex,
    debug_mode: bool = False,
    *,
    plateau_check: bool = True,
    plateau_probe_radius: Optional[float] = None,
    plateau_area_threshold: Optional[float] = None,
    plateau_cluster_tol: Optional[float] = None,
    mu1_low: float = -1,
    mu1_high: float = 1,
    mu2_low: float = -1,
    mu2_high: float = 1,
    continuum_tol: Optional[float] = None,
    continuum_perturb: Optional[float] = None,
    max_iter: Optional[int] = None,
    wtol: Optional[float] = None,
    max_range_expansions: Optional[int] = None,
    range_expand_factor: Optional[float] = None,
) -> GBZResult:
```

Main entry point. Builds `CharPoly`, runs `bisect_amoeba_ronkin_min`, then assembles subsets from the solved ZM. `debug_mode=True` re-raises solver/assembly exceptions instead of returning a failed `GBZResult`. Polynomial construction and signature validation occur outside that handler. Check `success` before interpreting empty subsets as an exterior energy.

Only the explicitly named keyword options above are accepted. `None` defaults resolve from the home-module constants at call time; see [constants.md](constants.md). There is no `perc` parameter or catch-all options dictionary.

The former bisection `xtol` keyword is now `wtol`: it bounds winding, not μ-bracket width. Crossing refinement retains its independent angular `xtol`, and coarse μ₂ bisection retains `coarse_xtol` for bracket width.

### 3.2 `bisect_amoeba_ronkin_min`

```python
def bisect_amoeba_ronkin_min(
    char_poly: CharPoly,
    E_ref: complex,
    mu1_low: float = -1,
    mu1_high: float = 1,
    mu2_low: float = -1,
    mu2_high: float = 1,
    continuum_tol: Optional[float] = None,
    continuum_perturb: Optional[float] = None,
    max_iter: Optional[int] = None,
    wtol: Optional[float] = None,
    max_range_expansions: Optional[int] = None,
    range_expand_factor: Optional[float] = None,
    frac: Optional[float] = None,
) -> dict:
```

Returns `mu1`, `mu2`, `zeros`, `is_continuum`, `_mu1_bracket`, `_w1_bracket`, `_exit_reason`, `_w1_area` (discrete path), `_continuum_members` (continuum path), `_zm`.

`zeros` contains `(theta1, theta2)` pairs. `frac` is a compatibility
argument and is currently unused. Iteration or range-expansion exhaustion
raises; there is no successful fallback at the midpoint of a failed solve.

### 3.3 `find_crossings`

Import this helper from `pygbz2d.amoeba.zm_extract`; it is not re-exported
by `pygbz2d.amoeba`.

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
| `find_crossings` | Mesh-visible `ln|β₂| = μ₂` crossings with optional refinement and segment/column avoidance. |
| `calculate_a2_average_winding` | a2 average winding from crossings. |

### `bisect.py`

| Function | Purpose |
|----------|---------|
| `_try_fast_mu2` | Gap test for groups selected at θ₁=0; returns `{ok, A, B, lo_cols, hi_cols}`. |
| `_screen_extremum_intervals` / `_refine_track_extrema` | Endpoint derivative screen + bracketed actual-derivative solves, retaining all new root rows for later level crossings. |
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
| `_point_near_any_line` | Require matching mu1, circular theta1, and chordal beta2 at the same line sample before dropping a point. |
| `_check_zeros_are_clustered` / `_probe_zero_plateau_near_mu1` | Plateau detection. |

### `ronkin_winding.py`

| Function | Purpose |
|----------|---------|
| `_find_exact_crossing` | Bracketed polynomial-root refinement of one crossing, retaining new ZM rows. |
| `_get_average_winding_from_zeros` | Average winding and non-zero area from a zero partition. |

## 5. Package Structure

```
pygbz2d/amoeba/
├── __init__.py            # Public API exports
├── amoeba.py              # collect_GBZ_subsets, subset assembly, plateau detection
├── bisect.py              # μ₁/μ₂ bisection, fast μ₂ gap test, extremum refinement
├── ronkin_winding.py      # bracketed crossing refinement, average winding from zeros
└── zm_extract.py          # AmoebaZeroManager, detect_continuum, find_crossings,
                           #   calculate_a2_average_winding
```

Public API: `CharPoly`, `bisect_amoeba_ronkin_min`, `AmoebaZeroManager`, `calculate_a2_average_winding`, `collect_GBZ_subsets`, `PointSubset`, `LineSubset`, `GBZResult`, `ConnectedSubset`.
