# Numerical Constants Reference

Every numerical constant in pyGBZ2d lives in the module that consumes it
(single-consumer locality); seven shared settings live in `pygbz2d.core`,
alongside the mathematical constant `TWO_PI`. This document maps the named
numerical settings in the core and experimental modules.

## The customization model

There are exactly two ways to tune a computation — no third channel:

1. **Per call** — pass the keyword argument on a public entry point
   (defaults left as `None` resolve from the constants below at call time):

   ```python
   from pygbz2d.sgbz import collect_GBZ_subsets
   gbz = collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, continuum_tol=1e-8)
   ```

2. **Per process** — assign the module constant directly.  The assignment
   takes effect on the next lookup of that module attribute. Decorated
   arguments resolve at function entry; an already resolved local argument
   retains its value for that call. Existing `StepControl` objects likewise
   retain values captured at construction. Direct module reads inside a
   running loop see later assignments:

   ```python
   from pygbz2d.sgbz import pairwise
   pairwise.CROSSING_TOL = 1e-12        # every later call sees 1e-12
   ```

   Scope caveats: "global" means *this Python process*.  On Linux fork,
   child processes inherit values assigned **before** the pool is created;
   later assignments do not propagate to already-running workers. Spawned
   workers import modules afresh and need overrides applied in each worker.

Two rules keep this mechanism sound (enforced by AST lint in
`tests/test_constants.py`):

- Tunable constants are not imported **by value** by consuming modules
  (`from x import CONST` copies the value at import time — the frozen-copy
  trap). Cross-module reads use attribute access: `pairwise.CROSSING_TOL`.
  Package `__init__.py` re-exports are exempt; assign the home-module value,
  rather than the exported copy, to change a default.
- The AST lint rejects numeric literals other than 0/1 outside its input
  whitelist. Public solver entries resolve `None` through `live_defaults`;
  some lower-level helpers resolve it manually. `estimate_error` is an
  exception: its direct-call `atol`/`rtol` defaults are captured at import,
  while `arclength_step` passes values from its `StepControl` explicitly.

### Conventions

- ⚠ **machine-anchored** — the tolerance is tied to float64 resolution.
  Retuning it casually changes discrete decisions; leave it alone unless
  you understand the guard it feeds.
- *Search ranges* (`mu1_guess`, `mu2_low/high`, `mu1_low/high`)
  are model-scale **inputs**, not numerical constants: they stay plain
  keyword arguments and never become module constants.

---

## `pygbz2d.core` — cross-package constants

| Constant | Default | Meaning |
|---|---|---|
| `CONTINUUM_TOL` | 1e-6 | Width of the continuum (degenerate-band) tolerance in μ-space; SGBZ tie detection and amoeba band detection share it. |
| `CONTINUUM_FRAC` | 0.9 | Vote fraction of in-band mesh rows above which an ItemView counts as a continuum cluster. |
| `CONTINUUM_PERTURB` | 1e-4 | Base μ-perturbation for escaping a continuum band when probing the two winding limits (×1 member of `ESCAPE_LADDER`). Unified 2026-08 (SGBZ previously used 1e-2, amoeba 1e-4). |
| `WINDING_ZERO_TOL` | 1e-8 | Live default for SGBZ `zero_tol` in its solve and plateau probe, and amoeba `wtol` in μ₁/μ₂ solves and plateau probes. |
| `PLATEAU_CLUSTER_TOL` | 1e-2 | Torus-clustering radius for plateau probes (PMGBZ points / zeros); shared by the SGBZ and amoeba probes. |
| `ESCAPE_LADDER` | (1.0, 2.0, 4.0, 8.0) | Scale factors applied to `CONTINUUM_PERTURB` when one step fails to escape a degenerate band. |
| `PROBE_XTOL` | 1e-10 | Compatibility default for `generate_probe_steps(xtol=None)`; the current step-generation implementation does not use xtol. |

`generate_probe_steps` includes a geometric ladder from
`max(10*zero_tol, 1e-12)` to
`max(probe_radius, 4*bracket_width, 100*zero_tol, 1e-12)`, together with
positive quarter/half multiples of the input scales. `probe_radius` is a
scale, not a hard upper bound on the probes.

## `pygbz2d.continuation.arclength` — RK45-style stepper

All single-consumer knobs of the pseudo-arclength step controller.
`StepControl` fields left as `None` resolve from these **at construction**,
so assigning `arclength.SAFETY` reaches every controller built afterwards.

| Constant | Default | Meaning |
|---|---|---|
| `SAFETY` | 0.9 | Safety factor in the current-error step-size update. |
| `MIN_FACTOR` / `MAX_FACTOR` | 0.2 / 10.0 | Step-size shrink/grow clamps per rejection/acceptance. |
| `ERROR_EXPONENT` | −0.5 | Error-norm exponent, −1/(p+1) for the first-order tangent predictor. |
| `STEP_ATOL` / `STEP_RTOL` | 1e-12 / 1e-3 | Step-acceptance tolerances (prediction vs. re-solved roots). |
| `STEP_MAX_ITER` | 20 | Rejection iterations allowed inside one step. |
| `MAX_STEP` / `MIN_STEP` | 0.5 / 1e-12 | Hard step bounds. ⚠ `MIN_STEP` is machine-anchored. |
| `ZERO_THRESHOLD` / `INF_THRESHOLD` | 1e-6 / 1e6 | Finite-radius thresholds for treating a root as singular; includes finite roots near 0/∞ as well as exact padding roots. |
| `PREDICT_MAX_ABS_ARG` | 1.0 | \|Vⱼ · Δθ₁\| above which the tangent prediction is held fixed (prevents `exp` overflow). |

## `pygbz2d.continuation.multiple_roots` — MR detection

| Constant | Default | Meaning |
|---|---|---|
| `CLUSTER_TOL` | 1e-4 | Chordal-distance threshold of `detect_cluster` — the single cluster predicate. Unified 2026-08 to the run-side value (the old direct-call default 1e-6 is retired). Also the zero-coincidence tolerance that groups per-pair Brent candidates into MR events. |
| `MIN_DIST_THRESHOLD` | 0.1 | Pair distance below which that pair's interval-trigger tracking arms; above it the pair's state resets (per-pair, not global). |

## `pygbz2d.continuation.zero_manager` — ZM run loop

The MR restart trio (`MR_JUMP`, `MR_RESTART_FACTOR_H0`,
`MR_RESTART_FACTOR_ABS`) is **one formula's** parameters — the effective
restart distance past an MR is
`max(MR_JUMP, H0_FACTOR·min_dtheta/h0, ABS_FACTOR·min_dtheta)`.  Keep them
together.

| Constant | Default | Meaning |
|---|---|---|
| `H0` | 0.1 | Initial θ₁ arclength step. |
| `MIN_DTHETA` | 1e-10 | Step collapse threshold for the MR point trigger (neighbour of `arclength.MIN_STEP`: that is the stepper's floor, this is the MR detector's sensitivity — two roles, deliberately distinct). |
| `MR_JUMP` | 1e-6 | Base restart distance past a refined MR. |
| `MR_DENSE_MAX_STEP` | 1e-4 | Dense sampling between MR events of one trigger bracket: maximum θ₁ spacing of the regular rows sampled between two consecutive events. |
| `MR_DENSE_MIN_SAMPLES` | 8 | Dense sampling between bracket events: minimum number of interior sample rows between two consecutive events. |
| `MR_RESTART_FACTOR_H0` / `MR_RESTART_FACTOR_ABS` | 10.0 / 100.0 | Branch-point-safe floors on the restart distance (too-close restarts re-detect the same MR). |
| `MR_REDETECT_RETRY_FACTOR` / `MR_REDETECT_MAX_RETRIES` | 2.0 / 4 | Multiply the restart jump when the same MR cluster is re-detected without progress; bound the number of retries. |
| `MR_STUCK_TOL` ⚠ | 1e-12 | Two refined MR θ₁ closer than this = no forward progress (error). |
| `BOUNDARY_THETA_TOL` ⚠ | 1e-6 | \|θ − 2π\| below which a boundary MR is pinned to exactly 2π. |
| `MR_GAUGE_TOL` | 0.1 | Warn when the iterative MR solver's θ₁ drifts more than this from the trigger. |
| `MAX_SEGMENTS` ⚠ | 10000 | Hard cap on tracked segments (runaway protection). |

## `pygbz2d.sgbz.pairwise` — the crossing channel

| Constant | Default | Meaning |
|---|---|---|
| `CROSSING_TOL` | 1e-10 | brentq `xtol` when refining a pairwise ln\|β₂\| crossing; also the EventGroup merge distance. |
| `MIN_DIRECTION_DERIV` | 1e-12 | Minimum \|Re(V_a) − Re(V_b)\| for a trustworthy crossing direction; below it the event is a tangent touch. |
| `REFINE_MAX_ROUNDS` | 3 | Multi-crossing mesh-refinement rounds. |
| `REFINE_SAFETY_FACTOR` | 4.0 | Target sub-intervals per narrowest predicted crossing gap, including gaps to interval endpoints. |
| `REFINE_MAX_SUBINTERVALS` | 64 | Sub-mesh cap per interval. |
| `REFINE_MAX_TOTAL_INSERTS` | 2000 | Budget for the regular refinement rounds. Reaching it triggers a warning and one final unbudgeted insertion pass, so it is not a hard total cap. |
| `REFINE_REL_TOL` ⚠ | 1e-12 | Real-root / duplicate-θ filter, relative to max(1, interval length). |
| `THETA_EQ_TOL` ⚠ | 1e-15 | Exact-endpoint float comparison when refinement reads mesh rows. |
| `BRENTQ_MAXITER` | 100 | brentq iteration budget (the bracket is guaranteed by the sign scan). |

## `pygbz2d.sgbz.mu2mid`

| Constant | Default | Meaning |
|---|---|---|
| `LOGABS_CLAMP` | 14.0 | Log-modulus clamp for boundary values and path pieces: radii range from e^-14 ≈ 8.3e-7 to e^14 ≈ 1.2e6. Saturated contributions have zero slope; non-finite unsaturated derivatives cause linear interpolation independently of the clamp magnitude. |

## `pygbz2d.sgbz.winding`

| Constant | Default | Meaning |
|---|---|---|
| `WINDING_QUAD_EPSABS` / `WINDING_QUAD_EPSREL` | 1e-3 / 1e-3 | `scipy.integrate.quad` tolerances; the average-winding solver rounds each seed winding to an integer, while `get_winding_number` itself returns an unrounded value. |
| `WINDING_QUAD_LIMIT` | 200 | quad sub-interval budget. |
| `SEED_N_PER_INTERVAL` | 4 | θ₂ samples per interval when picking the loop-winding seed. |

## `pygbz2d.sgbz.sgbz_solver`

| Constant | Default | Meaning |
|---|---|---|
| `MU1_MAX_ITER` | 60 | Iteration budget of the μ₁ bisection. |
| `MAX_BRACKET_EXPANSIONS` | 10 | Cap on μ₁ bracket-expansion steps (error past this). |

## `pygbz2d.amoeba.bisect`

| Constant | Default | Meaning |
|---|---|---|
| `BISECT_MAX_ITER` | 60 | μ₂ bisection budget. |
| `BISECT_COARSE_XTOL` | 1e-3 | μ₂ bisection coarse-stage bracket-width tolerance (unrefined crossings). |
| `EXTREMUM_INSERT_REL_TOL` | 1e-12 | Angular resolution for actual-derivative extremum refinement: Brent's absolute tolerance is this value times `max(1, abs(theta_lo), abs(theta_hi))`. All new solve samples are retained in the mesh. |
| `MAX_RANGE_EXPANSIONS` | 10 | μ₂ search-range expansion cap. |
| `RANGE_EXPAND_FACTOR` | 2.0 | Range growth per expansion step. |

Fine μ₂ and outer μ₁ bisection use `wtol`, resolved from `core.WINDING_ZERO_TOL`; the former `BISECT_XTOL` constant has been removed.

## `pygbz2d.amoeba.ronkin_winding`

| Constant | Default | Meaning |
|---|---|---|
| `CROSSING_XTOL` | 1e-12 | Absolute θ₁ tolerance for bracketed crossing refinement. Replaces `FSOLVE_XTOL`. |
| `CROSSING_MAXITER` | 500 | Brent iteration budget per crossing. Replaces `FSOLVE_MAXFEV` (formerly a function-evaluation budget). |

## `pygbz2d.amoeba.amoeba`

| Constant | Default | Meaning |
|---|---|---|
| `PLATEAU_AREA_THRESHOLD` | 1e-2 | Threshold on each normalized integral of absolute winding, `sum(abs(u)*width)/(2*pi)`; this differs from the torus radius `PLATEAU_CLUSTER_TOL` and can exceed one when winding magnitudes exceed one. |
| `SNAP_TOL` | 1e-3 | Sample snap tolerance: circular θ₁ and chordal β₂ distances must both be below this value to remove a PointSubset; μ₁ must also match within this tolerance. |

## `pygbz2d.experimental.band_clustering`

These settings belong to the experimental post-processing API. Angular
coordinates are embedded as cosine/sine pairs and energy coordinates are
scaled by the original energy-grid spacing. Line samples are not decimated.

| Constant | Default | Meaning |
|---|---|---|
| `ALPHA_E` | 0.25 | Energy-block weight in the eight-dimensional embedding. |
| `EPS_SCAN_MIN` / `EPS_SCAN_MAX` | 0.02 / 3.0 | Geometric scan bounds for the radius graph; stop early at the first full merge. |
| `N_SCAN_STEPS` | 21 | Number of candidate radii before early stopping. |
| `MIN_LARGEST_FRAC` | 0.1 | Minimum largest-cluster share for an eligible plateau. |
| `MIN_CLUSTERS_PLATEAU` | 2 | Minimum cluster count for an eligible plateau. |
| `MAX_MARGIN_PAIRS` | 12 | Number of largest clusters selected; margins are computed for all pairs among them. |
| `KNN_K` | 2 | Neighbor order, including self; two selects the nearest non-self point. |
| `MIN_CLUSTER_SIZE` | 5 | Size threshold for ordering small clusters after larger ones in summaries; does not discard them. |

---

## Same-pattern families (deliberately NOT merged)

These constants share a *shape* but act on different quantities; merging
them would couple unrelated scales:

- **eps-scale guards**: `THETA_EQ_TOL` (θ endpoints), `REFINE_REL_TOL`
  (relative interval length), `MR_STUCK_TOL` (θ progress),
  `arclength.MIN_STEP` (step size), `CROSSING_XTOL` (θ₁).
- **convergence tolerances**: `CROSSING_TOL` (θ via brentq),
  `WINDING_ZERO_TOL` (winding). `PROBE_XTOL` is retained for compatibility
  but does not currently affect probe distances.
- **step floors**: `arclength.MIN_STEP` (stepper abort) vs
  `zero_manager.MIN_DTHETA` (MR trigger sensitivity) — written at different
  times, both kept; see the zero_manager table.

## Historical merges (2026-08, correlation analysis)

| Former constants | Merged into |
|---|---|
| `zero_tol` (sgbz) + `winding_tol_floor` (amoeba), both 1e-10 | `core.WINDING_ZERO_TOL` |
| sgbz `tol_normalized=1e-2` literal + amoeba `plateau_cluster_tol` | `core.PLATEAU_CLUSTER_TOL` |
| `(1, 2, 4, 8)` ladders ×3 files | `core.ESCAPE_LADDER` |
| `min_dist_threshold` defined in two files | `continuation.multiple_roots.MIN_DIST_THRESHOLD` |
| `detect_cluster` defaults 1e-4 (run side) / 1e-6 (direct) | `continuation.multiple_roots.CLUSTER_TOL = 1e-4` |
| SGBZ `continuum_perturb` 1e-2 / amoeba 1e-4 | `core.CONTINUUM_PERTURB = 1e-4` |
