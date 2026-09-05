# continuation — Pseudo-Arclength Continuation Root Tracking

Adaptive-step root tracking along $\theta_1$ using analytic derivatives, Hungarian matching,
and pseudo-arclength step-size control.  The top-level `ZeroManager` class orchestrates
full-circle integration with automatic multiple-root detection and refinement.

## 1. Motivation

Both SGBZ and amoeba now use `continuation.ZeroManager` as their root-solving backend,
which replaces the prior fixed uniform $\theta_1$ meshes.  When roots change rapidly —
near degeneracies, crossings, or poles — a uniform mesh may miss features or produce
incorrect Hungarian matchings.

The pseudo-arclength approach adapts the step size to the local root dynamics: steps are
small where roots move quickly, large where they are quiescent.

## 2. Algorithm

### 2.1 Pseudo-arclength parameterization

The implicit equation $f(E, \beta_1, \beta_2) = 0$ defines $\beta_2$ as a multi-valued
function of $\theta_1$ (at fixed $E$, $\mu_1$).  Differentiating:

$$\frac{\partial f}{\partial\beta_1} \cdot i\beta_1\, d\theta_1 +
  \frac{\partial f}{\partial\beta_2} \cdot d\beta_2 = 0$$

$$\Rightarrow\quad
\frac{d\beta_2}{d\theta_1} =
-i\beta_1\frac{\partial f/\partial\beta_1}{\partial f/\partial\beta_2}$$

The state variable for continuation is $\ln\beta_2$ (not $\beta_2$), giving the tangent:

$$V_j = \frac{d\ln\beta_{2,j}}{d\theta_1} =
\frac{1}{\beta_{2,j}}\frac{d\beta_{2,j}}{d\theta_1}$$

The full tangent vector in $(\theta_1, \ln\beta_{2,1}, \ldots, \ln\beta_{2,n})$-space is:

$$\mathbf{V} = \begin{pmatrix} 1 \\ V_1 \\ \vdots \\ V_n \end{pmatrix},
\quad \|\mathbf{V}\|_2 = \sqrt{1 + \sum_j |V_j|^2}$$

The arclength step condition: $\Delta\theta_1 = h / \|\mathbf{V}\|_2$, so that one step
covers arclength $h$ in the state space.

### 2.2 Adaptive step-size control

Modelled after scipy's RK45 integrator (`scipy.integrate._ivp.rk.RungeKutta._step_impl`).

Each step:
1. Compute tangent $\mathbf{V}$, propose $\Delta\theta_1 = h / \|\mathbf{V}\|_2$.
2. Predict roots at $\theta_1 + \Delta\theta_1$ via $\beta_{2,j} \leftarrow \beta_{2,j} \cdot \exp(V_j \cdot \Delta\theta_1)$.
3. Solve actual roots at the new $\theta_1$ via `np.roots`.
4. Hungarian-match predicted → actual, compute max chordal distance as error.
5. If `error_norm < 1`: accept step, update $h \leftarrow h \cdot \min(10,\ 0.9 \cdot \text{error}^{-0.5})$.
6. If `error_norm ≥ 1`: reject step, $h \leftarrow h \cdot \max(0.2,\ 0.9 \cdot \text{error}^{-0.5})$, retry.

The error exponent $-0.5 = -1/(p+1)$ corresponds to a first-order method ($p=1$ for
linear tangent extrapolation).

| Constant | Value | Meaning |
|----------|-------|---------|
| `SAFETY` | 0.9 | Safety factor on step update |
| `MIN_FACTOR` | 0.2 | Max step decrease per rejection |
| `MAX_FACTOR` | 10.0 | Max step increase per acceptance |
| `ERROR_EXPONENT` | -0.5 | $-1/(p+1)$ for $p=1$ |

### 2.3 Handling $0$ and $\infty$

Roots at $\beta_2 = 0$ or $\beta_2 = \infty$ have undefined $\ln\beta_2$, so their tangent
components are set to $V_j = \mathrm{nan}$ (undefined — not zero) and their values are kept
unchanged during prediction.  $\|\mathbf{V}\|_2$ ignores the nan components; the inf
components produced by multiple roots are deliberately kept so that an MR drives
$\Delta\theta_1 \to 0$ and triggers the step-collapse detector.

Thresholds: $|\beta_2| < 10^{-6}$ (near 0), $|\beta_2| > 10^{6}$ (near $\infty$).

### 2.4 Multiple root detection

Two complementary triggers detect multiple roots during integration:

| Trigger | Mechanism | Catches |
|---------|-----------|---------|
| **Point trigger** (`multiple_root_point_trigger`) | $\Delta\theta_1 = h/\|\mathbf{V}\|_2$ collapses below `min_dtheta` | $\|\mathbf{V}\|_2$ diverges at the root — the solver is **at** the MR |
| **Interval trigger** (`MultipleRootIntervalTrigger`) | Sign of $d(\min|\beta_i-\beta_j|^2)/d\theta_1$ flips from approaching to separating | "Accidental" MR where tangent stays well-conditioned — the MR lies **between** steps |

The point trigger is essential for **generic double roots** — points where
$\partial f/\partial\beta_2 = 0$ but $\partial f/\partial\beta_1 \neq 0$.  At such
points $d\beta_2/d\theta_1$ diverges, making $\|\mathbf{V}\|_2 \to \infty$ and
$\Delta\theta_1 \to 0$, even though the tangent prediction error may remain small.

The interval trigger catches **accidental multiple roots** where the tangent stays
well-conditioned and the step size never collapses.  It computes the pairwise
distances and their $\theta_1$-derivatives for **every** root pair at once
(vectorized `_pairwise_dist_deriv`) and tracks each pair's derivative sign
across steps; a pair triggers when **its own** sign flips from negative
(approaching) to positive (separating), indicating a local minimum of that
pair's distance.  Per-pair tracking matters for **simultaneous degeneracies**:
when two pairs merge at the same $\theta_1$ (a symmetric double MR) their
distances are tied and the *closest-pair* identity flickers row by row — a
single closest-pair state with a same-pair guard (the retired design) blocks
the flip detection exactly at the sign change.  Each pair is armed only while
its own distance stays below `MIN_DIST_THRESHOLD`.  One stop emits a
`list[MRTriggerRecord]` — one record per flipped pair, all sharing the
trigger bracket $(\theta_{lo}, \theta_{hi})$.

**Non-generic double roots** (where both $\partial f/\partial\beta_1 = 0$ and
$\partial f/\partial\beta_2 = 0$) have a finite tangent via l'Hôpital's rule and
are handled by Hungarian matching without special detection.

Both triggers surface their result through `SegmentResult.stop_reason`:

| `stop_reason` | Meaning |
|---------------|---------|
| `StopReason.completed` | Reached $\theta_1 = 2\pi$ without hitting an MR |
| `StopReason.multiple_root_encountered` | Step-size collapsed → refine at the current $\theta_1$ |
| `StopReason.multiple_root_in_interval` | One or more pairs' derivative sign flipped → per-pair Brent on each record's bracket (`mr_triggers`) |

### 2.5 Multiple root refinement

Two complementary approaches for locating multiple roots:

**Brent-based (1D, fixed μ₁):** `solve_multiple_roots_in_interval` locates the exact
$\theta_1$ of a multiple root via Brent's method on a pairwise distance
derivative.  With `min_pair` given, $g(\theta)$ is **that pair's own**
derivative — continuous across closest-pair identity changes, which an
argmin-based $g$ is not exactly at the flip; without it, the closest pair's
derivative is used.  The root pair is tracked via Hungarian matching to
maintain identity across the bracketing interval.

**Newton-based (4D, free β₁):** `solve_multiple_roots_iterative` solves the
$4 \times 4$ real system

$$\begin{cases}
\operatorname{Re} f(E, \beta_1, \beta_2) = 0 \\
\operatorname{Im} f(E, \beta_1, \beta_2) = 0 \\
\operatorname{Re} \frac{\partial f}{\partial\beta_2}(E, \beta_1, \beta_2) = 0 \\
\operatorname{Im} \frac{\partial f}{\partial\beta_2}(E, \beta_1, \beta_2) = 0
\end{cases}$$

for $(\operatorname{Re}\beta_1, \operatorname{Im}\beta_1, \operatorname{Re}\beta_2,
\operatorname{Im}\beta_2)$ using `scipy.optimize.root` with analytic Jacobian from
`CharPoly.eval_partials_2`.  This directly finds the point in $(\beta_1, \beta_2)$
space where $f = 0$ **and** the root is degenerate ($\partial f/\partial\beta_2 = 0$).
It requires a good initial guess `(beta1_approx, beta2_approx)` — typically obtained
from the Brent-based solver or the ZeroManager's MR detection — and returns the
refined `(beta1_mr, beta2_mr)`.

**Cluster detection:** `detect_cluster` finds all touching root groups via connected
components of the chordal-distance proximity graph (using
`scipy.sparse.csgraph.connected_components`).

### 2.6 ZeroManager orchestration

`ZeroManager.run()` integrates the full circle $[0, 2\pi)$:

1. **Initial check**: solve roots at $\theta_1 = 0$, run `detect_cluster`.  If an MR
   is found, record it and jump past it with a small `mr_jump` step.
2. **Segment loop**: call `integrate_segment` to advance toward $2\pi$.  On each stop:
   - **completed** → Hungarian-match the right boundary to the left boundary
     ($\theta_1 = 0$), store `boundary_perm`, and finish.
   - **MR encountered** → refine via `solve_multiple_roots_iterative`; a single
     event per stop.
   - **MR in interval** → per triggered pair, Brent-refine that pair's own
     derivative zero inside its bracket, then group the candidate $\theta_1$s
     into events by **zero coincidence** (a candidate belongs to an existing
     event when its pair is already within `CLUSTER_TOL` at that event's
     $\theta_1$ — the same physical degeneracy reached through different track
     columns).  Simultaneous degeneracies of different zero pairs merge into
     ONE event with several clusters; genuinely distinct nearby MRs stay
     separate events, bridged by segments of densely sampled regular rows
     (spacing ≤ `MR_DENSE_MAX_STEP`, ≥ `MR_DENSE_MIN_SAMPLES` interior rows).
     If a cluster is confirmed, record the MR and its roots; otherwise treat
     as a false positive and merge the segments.
3. **False-positive MR (`pending`)**: when `_refine_mr` reports `cluster=[]` — the
   trigger fired but `detect_cluster` finds no real cluster at the located $\theta_1$ —
   the segment is held back in a `_PendingSeg` (its data minus the false-MR row, plus
   the current `left_mr`) instead of being appended with a sentinel.  The next
   iteration's segment is vstacked onto the pending one, preserving `left_mr`, and the
   merged segment is then closed at the next real MR or at $2\pi$.  This replaces the
   former `right_mr == -2` post-loop patch.
4. **Fallback**: if `completed` was never reached (e.g. all segments ended at MRs),
   manually match the last segment's right boundary to the left boundary.

Results are stored in:
- `multiple_roots: list[MultipleRootInfo]` — each with `.theta1`, `.cluster_indices`
  (list of tuples, one per cluster), and `.roots` (modulus-sorted).
- `segments: list[SegmentData]` — each with `.theta1_arr`, `.tracked_roots`,
  `.abs_argsort`, `.left_mr`, `.right_mr` (indices into `multiple_roots`, -1 if none).

### 2.7 Root sorting

Roots are tracked via Hungarian matching (chordal distance on the Riemann sphere)
at every accepted step, maintaining continuous identity across $\theta_1$.

Each `SegmentData` stores:
- `tracked_roots`: track-ordered (column $j$ = physical root $j$), shape $(N, K)$.
- `abs_argsort`: per-row `np.argsort(|β₂|)`, shape $(N, K)$ — for compatibility
  with SGBZ's `solve_roots_on_mesh` output.

After the full $[0, 2\pi)$ loop, `ZeroManager.boundary_perm` stores the permutation
from the right boundary ($\theta_1 = 2\pi$) to the left boundary ($\theta_1 = 0$).

### 2.8 Unified Hermite interpolation

The discrete root mesh is turned into a continuous curve by two-point cubic
Hermite interpolation.  All construction sites share `continuation/interpolation.py`:

- `cubic_hermite_poly(h, v0, dv0, v1, dv1)` — the pure cubic builder
  ($f(0)=v0,\ f'(0)=dv0,\ f(h)=v1,\ f'(h)=dv1$); no input validation.
- `hermite_interp_poly(...)` — checks whether `dv0` and `dv1` are finite:
  both finite → the cubic above (`len(poly) == 4`); either divergent/undefined
  → the linear polynomial `[slope, v0]` with `slope = (v1 - v0) / h`
  (`len(poly) == 2`).

Both follow the numpy `poly` convention (highest power first), so evaluation is
`np.polyval`, differentiation `np.polyder`, and curve intersection
`np.roots(np.polysub(p1, p2))`; `len(poly)` distinguishes cubic from linear.
Consumers: `arclength.predict_roots_hermite` (root-track prediction),
`Mu2MidZM` sort-change refinement and μ₂_mid path evaluation, and SGBZ crossing
bracketing.

## 3. API

### 3.1 Top-level entry point

```python
from continuation import ZeroManager

zm = ZeroManager(poly, E_ref, mu1)
zm.run(
    h0=0.1,             # initial arclength step
    ctrl=StepControl(),  # RK45-style tolerances and step bounds
    min_dtheta=1e-10,    # Δθ₁ threshold for point trigger
    cluster_tol=1e-4,    # chordal-distance threshold for clusters
    mr_jump=1e-6,        # θ₁ step to jump past a refined MR
    verbose=False,       # print progress to stdout
)

# Results:
zm.multiple_roots      # list[MultipleRootInfo]
zm.segments            # list[SegmentData]
zm.n_multiple_roots    # int
zm.n_segments          # int
zm.has_boundary_mr     # bool — MR at θ₁=0?
zm.boundary_perm       # np.ndarray (K,) — right→left permutation
```

`MultipleRootInfo` (NamedTuple):
```python
mr.theta1            # float — θ₁ of the multiple root
mr.cluster_indices   # list[tuple[int, ...]] — root indices grouped by cluster
mr.roots             # np.ndarray (K,) — modulus-sorted β₂ roots at this θ₁
```

`SegmentData` (dataclass):
```python
seg.theta1_arr       # np.ndarray (N,) — monotonic θ₁ mesh
seg.tracked_roots    # np.ndarray (N, K) — track-ordered β₂ roots
seg.abs_argsort      # np.ndarray (N, K) — per-row |β₂| argsort
seg.left_mr          # int — MR index at left boundary (-1 = none)
seg.right_mr         # int — MR index at right boundary (-1 = none)
seg.tangents         # np.ndarray (N, K) — per-row analytic V_j, always built by run();
```

`ZeroManager.insert_solution(theta1, seg_idx=None, i=None, *, interp='hermite')` is a **mutating** operation; `interp='linear'` uses the two-point linear matching anchor (used by SGBZ pairwise crossing insertion) and is otherwise identical:
it solves the β₂ roots at an interior θ₁, reorders them onto the segment's track frame via
Hungarian matching against `interpolate_roots` (the prediction anchor), and splices the new
row into `SegmentData.theta1_arr`, `tracked_roots`, `abs_argsort`, and `tangents` *in place*
(at mesh index `i + 1`).  Any external references to those arrays are invalidated by the
splice — callers must re-fetch `segments[seg_idx]` after `insert_solution` returns.  This is
the hook `Mu2MidZM` overrides to densify the mesh at μ₂-refinement sites.

### 3.2 Low-level functions

```python
from continuation import (
    compute_tangent, predict_roots, estimate_error, arclength_step,
    cubic_hermite_poly, hermite_interp_poly,
    multiple_root_point_trigger, MultipleRootIntervalTrigger,
    detect_cluster, solve_multiple_roots_in_interval,
    integrate_segment,
    StepResult, StepControl,
)
```

| Function | Returns | Purpose |
|----------|---------|---------|
| `compute_tangent(poly, E_ref, beta1, roots)` | `(V, norm_V)` | Tangent vector and its norm |
| `predict_roots(roots, V, dtheta1)` | `predicted` | First-order tangent extrapolation |
| `cubic_hermite_poly(h, v0, dv0, v1, dv1)` | `np.ndarray` (len 4) | Pure cubic Hermite poly (numpy poly order) |
| `hermite_interp_poly(h, v0, dv0, v1, dv1)` | `np.ndarray` (len 4 or 2) | Cubic Hermite, or linear fallback when an endpoint derivative is not finite |
| `estimate_error(predicted, actual)` | `error_norm` | Chordal-distance error norm |
| `arclength_step(poly, E_ref, mu1, theta1, roots, h, ctrl=StepControl())` | `StepResult` | One adaptive step |
| `multiple_root_point_trigger(dtheta, *, min_dtheta=1e-10)` | `bool` | Step-size collapse check |
| `MultipleRootIntervalTrigger(min_dist_threshold)` | callable | Per-pair sign-flip detector; returns `list[MRTriggerRecord]` |
| `detect_cluster(roots, *, cluster_tol=1e-6)` | `list[tuple[int,...]]` | Connected components of close roots |
| `solve_multiple_roots_in_interval(poly, E_ref, mu1, theta1_left, theta1_right, roots_ref, min_pair=None)` | `theta1_mr` | Brent refinement (1D, fixed μ₁) |
| `solve_multiple_roots_iterative(poly, E_ref, beta1_approx, beta2_approx)` | `(beta1_mr, beta2_mr)` | Newton refinement (4D, free β₁) |
| `integrate_segment(poly, E_ref, mu1, theta_start, roots_start, theta_end, *, h0, ctrl, min_dtheta, min_dist_threshold)` | `SegmentResult` | Segment integration |

## 4. Synthetic Test Polynomials

The test suite includes analytically constructed polynomials for stress-testing:

| Polynomial | Roots | Tests |
|-----------|-------|-------|
| **A** $(\beta_2-\beta_1-1/\beta_1)(\beta_2-2)$ | $\beta_1+1/\beta_1,\ 2$ | Non-monomial tangent, zero-crossing, non-generic double root |
| **B** $(\beta_2-\beta_1)(\beta_2-1)$ | $\beta_1,\ 1$ | Monomial tangent (constant), exact double root |
| **C** $(\beta_2-\beta_1)(\beta_2-1/\beta_1)(\beta_2-1)$ | $\beta_1,\ 1/\beta_1,\ 1$ | Triple root, double root at $\pi$ |
| **D** $(\beta_2-0.5\beta_1)(\beta_2-1.5\beta_1)(\beta_2-2\beta_1)$ | $0.5\beta_1,\ 1.5\beta_1,\ 2\beta_1$ | Well-separated identical-k roots, modulus ordering |
| **E** $(\beta_2-\beta_1)(\beta_2-(1+\varepsilon)\beta_1)$ | $\beta_1,\ (1+\varepsilon)\beta_1$ | Near-degenerate monomial, Hungarian matching stress |
| **F** $(\beta_2-1)^2 - (\beta_1-1)$ | $1\pm\sqrt{\beta_1-1}$ | **Generic double root** — $\partial f/\partial\beta_2=0$, $\partial f/\partial\beta_1\neq0$, tangent diverges |

## 5. Key Design Decisions

1. **np.roots, not ODE integration**.  The tangent is only used for step-size
   control and prediction; roots are always solved exactly at each step via
   `np.roots`.  This avoids error accumulation from ODE integration.

2. **No artificial tangent cap**.  When $\partial f/\partial\beta_2 \approx 0$,
   $V_j$ grows naturally, $\|\mathbf{V}\|_2 \to \infty$, $\Delta\theta_1 \to 0$.
   The point trigger catches the divergence.  Capping would hide the
   problem and prevent detection.

3. **Riemann sphere chordal distance**.  Hungarian matching uses stereographic
   projection to $\mathbb{S}^2 \subset \mathbb{R}^3$, where $\infty$ maps to
   the north pole.  This gives a proper metric that handles roots at infinity
   without special cases.

4. **Two complementary MR triggers**.  The point trigger catches generic double
   roots where the tangent diverges; the interval trigger catches accidental
   double roots where the tangent stays well-conditioned.  Together they cover
   the full space of multiple-root behaviours.

5. **Compatible output format**.  `SegmentData.tracked_roots` and `.abs_argsort`
   mirror the shape amoeba's per-column extraction and SGBZ's `Mu2MidZM` /
   crossing detector expect, making it straightforward to plug into either
   pipeline.

## 6. File Layout

```
continuation/
├── __init__.py          # Re-exports public API
├── arclength.py         # Low-level step functions (~318 lines):
│                        #   StepControl (RK45-style tolerances bundled so
│                        #     arclength_step / integrate_segment / ZeroManager.run
│                        #     share one knob set), StepResult, compute_tangent,
│                        #     predict_roots, predict_roots_hermite (matching
│                        #     anchor for ZM / MR solvers), estimate_error,
│                        #   arclength_step
├── interpolation.py     # Unified Hermite polynomial builders:
│                        #   cubic_hermite_poly (pure cubic, no checks),
│                        #   hermite_interp_poly (automatic linear fallback when
│                        #     an endpoint derivative is not finite)
├── multiple_roots.py    # MR detection & refinement (~448 lines):
│                        #   MultipleRootInfo (with cluster_stds),
│                        #   snap_clusters_to_mean, multiple_root_point_trigger,
│                        #   MultipleRootIntervalTrigger, detect_cluster,
│                        #   _pair_distance_deriv, _closest_pair_deriv,
│                        #   solve_multiple_roots_in_interval (1D Brent, fixed μ₁),
│                        #   solve_multiple_roots_iterative (4D Newton, free β₁)
└── zero_manager.py      # Orchestration (~1100 lines):
                         #   StopReason, _MREndpoint, SegmentResult, SegmentData,
                         #   _PendingSeg (held-back false-positive-MR segment),
                         #   integrate_segment, ZeroManager (with .locate,
                         #   .interpolate_roots, .insert_solution — mutating,
                         #   overrides hook for Mu2MidZM — ._predict_roots_at_2pi,
                         #   ._boundary_perm_from_right, ._refine_mr)

tests/
├── test_continuation.py    # Low-level + integration tests
└── test_zero_manager.py    # ZeroManager-specific tests
```
