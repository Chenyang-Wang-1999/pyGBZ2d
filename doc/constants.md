# Numerical Constants Reference

Every numerical constant in bfGBZ2d lives in the module that consumes it
(single-consumer locality); only the seven cross-package constants live in
`bfgbz2d.core`.  This document is the complete map.

## The customization model

There are exactly two ways to tune a computation — no third channel:

1. **Per call** — pass the keyword argument on a public entry point
   (defaults left as `None` resolve from the constants below at call time):

   ```python
   from bfgbz2d.sgbz import collect_GBZ_subsets
   gbz = collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, continuum_tol=1e-8)
   ```

2. **Per process** — assign the module constant directly.  The assignment
   takes effect **immediately and process-wide, on the next read** (module
   attribute lookup at call time), including inside running loops:

   ```python
   from bfgbz2d.sgbz import pairwise
   pairwise.CROSSING_TOL = 1e-12        # every later call sees 1e-12
   ```

   Scope caveats: "global" means *this Python process*.  On Linux fork,
   child processes inherit values assigned **before** the pool is created;
   later assignments do not propagate to already-running workers.

Two rules keep this mechanism sound (enforced by AST lint in
`tests/test_constants.py`):

- Constants are never imported **by value** across modules
  (`from x import CONST` copies the value at import time — the frozen-copy
  trap).  Cross-module reads use attribute access: `pairwise.CROSSING_TOL`.
- Numeric signature defaults are forbidden outside the whitelist below;
  public entries use the `None` sentinel resolved by `live_defaults`.

### Conventions

- ⚠ **machine-anchored** — the tolerance is tied to float64 resolution.
  Retuning it casually changes discrete decisions; leave it alone unless
  you understand the guard it feeds.
- *Search ranges* (`mu1_guess`, `mu2_low/high`, `mu1_low/high`, `perc`)
  are model-scale **inputs**, not numerical constants: they stay plain
  keyword arguments and never become module constants.

---

## `bfgbz2d.core` — cross-package constants

| Constant | Default | Meaning |
|---|---|---|
| `CONTINUUM_TOL` | 1e-6 | Width of the continuum (degenerate-band) tolerance in μ-space; SGBZ tie detection and amoeba band detection share it. |
| `CONTINUUM_FRAC` | 0.9 | Vote fraction of in-band mesh rows above which an ItemView counts as a continuum cluster. |
| `CONTINUUM_PERTURB` | 1e-4 | Base μ-perturbation for escaping a continuum band when probing the two winding limits (×1 member of `ESCAPE_LADDER`). Unified 2026-08 (SGBZ previously used 1e-2, amoeba 1e-4). |
| `WINDING_ZERO_TOL` | 1e-10 | "\|winding\| counts as zero" predicate — SGBZ plateau probe (a₁) and the amoeba winding-tolerance floor (formerly `zero_tol` / `winding_tol_floor`, two names for one value). |
| `PLATEAU_CLUSTER_TOL` | 1e-2 | Torus-clustering radius for plateau probes (PMGBZ points / zeros); shared by the SGBZ and amoeba probes. |
| `ESCAPE_LADDER` | (1.0, 2.0, 4.0, 8.0) | Scale factors applied to `CONTINUUM_PERTURB` when one step fails to escape a degenerate band. |
| `PROBE_XTOL` | 1e-10 | Step resolution of the probe stepper (`core.generate_probe_steps`). |

## `bfgbz2d.continuation.arclength` — RK45-style stepper

All single-consumer knobs of the pseudo-arclength step controller.
`StepControl` fields left as `None` resolve from these **at construction**,
so assigning `arclength.SAFETY` reaches every controller built afterwards.

| Constant | Default | Meaning |
|---|---|---|
| `SAFETY` | 0.9 | PI-controller safety factor. |
| `MIN_FACTOR` / `MAX_FACTOR` | 0.2 / 10.0 | Step-size shrink/grow clamps per rejection/acceptance. |
| `ERROR_EXPONENT` | −0.5 | Error-norm exponent, −1/(p+1) for the first-order tangent predictor. |
| `STEP_ATOL` / `STEP_RTOL` | 1e-12 / 1e-3 | Step-acceptance tolerances (prediction vs. re-solved roots). |
| `STEP_MAX_ITER` | 20 | Rejection iterations allowed inside one step. |
| `MAX_STEP` / `MIN_STEP` | 0.5 / 1e-12 | Hard step bounds. ⚠ `MIN_STEP` is machine-anchored. |
| `ZERO_THRESHOLD` / `INF_THRESHOLD` | 1e-6 / 1e6 | \|β₂\| below/above which a root is a singular 0/∞ padding root. |

## `bfgbz2d.continuation.multiple_roots` — MR detection

| Constant | Default | Meaning |
|---|---|---|
| `CLUSTER_TOL` | 1e-4 | Chordal-distance threshold of `detect_cluster` — the single cluster predicate. Unified 2026-08 to the run-side value (the old direct-call default 1e-6 is retired). |
| `MIN_DIST_THRESHOLD` | 0.1 | Closest-pair distance below which the MR interval trigger arms; above it the trigger state resets. |

## `bfgbz2d.continuation.zero_manager` — ZM run loop

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
| `MR_RESTART_FACTOR_H0` / `MR_RESTART_FACTOR_ABS` | 10.0 / 100.0 | Branch-point-safe floors on the restart distance (too-close restarts re-detect the same MR). |
| `MR_STUCK_TOL` ⚠ | 1e-12 | Two refined MR θ₁ closer than this = no forward progress (error). |
| `BOUNDARY_THETA_TOL` ⚠ | 1e-6 | \|θ − 2π\| below which a boundary MR is pinned to exactly 2π. |
| `MR_GAUGE_TOL` | 0.1 | Warn when the iterative MR solver's θ₁ drifts more than this from the trigger. |
| `MAX_SEGMENTS` ⚠ | 10000 | Hard cap on tracked segments (runaway protection). |

## `bfgbz2d.sgbz.pairwise` — the crossing channel

| Constant | Default | Meaning |
|---|---|---|
| `CROSSING_TOL` | 1e-10 | brentq `xtol` when refining a pairwise ln\|β₂\| crossing; also the EventGroup merge distance. |
| `MIN_DIRECTION_DERIV` | 1e-12 | Minimum \|Re(V_a) − Re(V_b)\| for a trustworthy crossing direction; below it the event is a tangent touch. |
| `REFINE_MAX_ROUNDS` | 3 | Multi-crossing mesh-refinement rounds. |
| `REFINE_SAFETY_FACTOR` | 4.0 | Sub-mesh coarsening factor between rounds. |
| `REFINE_MAX_SUBINTERVALS` | 64 | Sub-mesh cap per interval. |
| `REFINE_MAX_TOTAL_INSERTS` | 2000 | Total inserted-row cap. |
| `REFINE_REL_TOL` ⚠ | 1e-12 | Real-root / duplicate-θ filter, relative to max(1, interval length). |
| `THETA_EQ_TOL` ⚠ | 1e-15 | Exact-endpoint float comparison when refinement reads mesh rows. |
| `BRENTQ_MAXITER` | 100 | brentq iteration budget (the bracket is guaranteed by the sign scan). |

## `bfgbz2d.sgbz.mu2mid`

| Constant | Default | Meaning |
|---|---|---|
| `LOGABS_CLAMP` | 14.0 | Clamp band for ln\|β₂\| in μ₂_mid path pieces (\|β₂\| = e^±14 ≈ 1.2e6); divergent boundary-root slopes exceed it and degrade pieces to linear. |

## `bfgbz2d.sgbz.winding`

| Constant | Default | Meaning |
|---|---|---|
| `WINDING_QUAD_EPSABS` / `EPSREL` | 1e-3 / 1e-3 | `scipy.integrate.quad` tolerances for the loop-winding integral; every caller rounds to an integer, so these stay loose by design. |
| `WINDING_QUAD_LIMIT` | 200 | quad sub-interval budget. |
| `SEED_N_PER_INTERVAL` | 4 | θ₂ samples per interval when picking the loop-winding seed. |

## `bfgbz2d.sgbz.sgbz_solver`

| Constant | Default | Meaning |
|---|---|---|
| `MU1_MAX_ITER` | 60 | Iteration budget of the μ₁ bisection. |
| `MAX_BRACKET_EXPANSIONS` | 10 | Cap on μ₁ bracket-expansion steps (error past this). |

## `bfgbz2d.amoeba.zm_extract`

| Constant | Default | Meaning |
|---|---|---|
| `SNAP_TOL` | 1e-3 | A crossing θ₁ within this of a continuum LineSubset endpoint snaps to the continuum/MR boundary (deliberately decoupled from `CONTINUUM_TOL`: snap radius vs. band sensitivity). |
| `ROOT_TOL` | 1e-9 | Root-match radius of the curve-consistency screen (shared MR rows match to machine precision; unrelated tracks differ by O(1)). |

## `bfgbz2d.amoeba.bisect`

| Constant | Default | Meaning |
|---|---|---|
| `BISECT_MAX_ITER` | 60 | μ₂ bisection budget. |
| `BISECT_XTOL` | 1e-10 | μ₂ bisection x-tolerance. |
| `MAX_RANGE_EXPANSIONS` | 10 | μ₂ search-range expansion cap. |
| `RANGE_EXPAND_FACTOR` | 2.0 | Range growth per expansion step. |

## `bfgbz2d.amoeba.ronkin_winding`

| Constant | Default | Meaning |
|---|---|---|
| `FSOLVE_XTOL` | 1e-12 | fsolve tolerance refining (θ₁, θ₂) crossings. |
| `FSOLVE_MAXFEV` | 500 | fsolve evaluation budget. |
| `CROSSING_RESIDUAL_TOL` | 1e-10 | Residual gate accepting a refined crossing (non-convergence falls back to the unrefined estimate). |

## `bfgbz2d.amoeba.amoeba`

| Constant | Default | Meaning |
|---|---|---|
| `PLATEAU_AREA_THRESHOLD` | 1e-2 | Non-zero-winding-area fraction below which the plateau pre-check arms (a separate quantity from `PLATEAU_CLUSTER_TOL`: area fraction vs. torus radius). |

---

## Same-pattern families (deliberately NOT merged)

These constants share a *shape* but act on different quantities; merging
them would couple unrelated scales:

- **eps-scale guards**: `THETA_EQ_TOL` (θ endpoints), `REFINE_REL_TOL`
  (relative interval length), `MR_STUCK_TOL` (θ progress),
  `arclength.MIN_STEP` (step size), `FSOLVE_XTOL` (residual).
- **1e-10 convergence tolerances**: `CROSSING_TOL` (θ via brentq),
  `BISECT_XTOL` (μ₂), `PROBE_XTOL` (probe steps).
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
