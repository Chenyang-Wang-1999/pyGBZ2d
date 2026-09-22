# pygbz2d.sgbz — Strip GBZ Computation

The SGBZ solver searches for a zero of the average major-axis winding
`W(E_ref, mu1)`. At each fixed `(E_ref, mu1)`, `ZeroManager` tracks the
polynomial's beta2 roots around the theta1 circle. SGBZ uses the mean of
the two boundary log-moduli as a variable transverse radius; amoeba uses
a constant transverse radius instead.

## 1. Boundary condition and data model

Let `M` be the beta2 Laurent denominator order and `K = M + N` the padded
root count. In ascending modulus order, the boundary roots occupy positions
`M-1` and `M` (zero-based). The boundary condition is
`|beta_(M-1)| = |beta_M|`.

- An isolated equality contributes `PointSubset` objects, one for each
  physical root column in the boundary event component. Isolated events can
  be transversal, tangent, or associated with a multiple root (MR).
- Equality over a theta1 interval contributes `LineSubset` objects, one per
  root track. A line may occupy only part of a tracked segment because
  another root can overtake the equal-modulus pair.

The theoretical intermediate curve is
`mu2_mid = (ln|beta_(M-1)| + ln|beta_M|) / 2`. Its numerical representation
is piecewise Hermite interpolation of solved mesh values and analytic
root tangents, with linear fallback for non-finite derivatives. It is an
approximation between mesh rows, with an analytically evaluated polynomial
derivative inside each piece.

A finite boundary pair requires `0 < M < K`; models with `M = 0` or `N = 0`
are outside this solver's supported boundary chart. Zero/infinite padding
roots remain at the outer sort positions. If either occupies the boundary
pair, ItemView construction raises rather than representing that boundary
by a finite clamped radius.

## 2. Analysis pipeline

### 2.1 Representative items and crossings

`Mu2MidZM(ZeroManager)` adds `analyze()` after root tracking:

1. **Build ItemViews.** For each segment, pair columns whose log-moduli
   differ by less than `tie_tol` on a fraction of rows strictly greater
   than `core.CONTINUUM_FRAC` (default 0.9). Graph traversal merges connected
   column pairs. Each cluster becomes one representative item with a
   multiplicity; other columns remain individual items. This avoids
   arbitrary sorting within near-equal continuum columns. The vote is a
   numerical detector for whole-segment equality, not a proof of it.
2. **Refine potential multiple crossings.** Cubic-Hermite predictions of
   pairwise log-modulus differences identify intervals with at least two
   separated predicted intersections. Direct polynomial solves insert a
   finer mesh, subject to round and insertion budgets. Reaching a limit
   triggers a warning and one final unbudgeted insertion pass. This improves
   detection of even numbers of crossings, but cannot certify completeness
   for features absent from the interpolant. Disable with
   `refine_multi_crossings=False`.
3. **Collect pair events.** `pairwise.collect_pair_events` scans every
   representative pair: an exact zero at a left endpoint is a touch;
   opposite endpoint signs bracket a crossing. Touches keep their mesh
   coordinates. Crossings use a linear prediction and Brent refinement of
   the actual log-modulus difference. Trial solves use
   `solve_at(interp='linear')` without mutating the mesh.
4. **Group and insert events.** Within a segment, consecutive events closer
   than `crossing_tol` form an `EventGroup`. A group uses an existing mesh
   row within its span when available, otherwise its midpoint. Circular
   seam merges are placed at zero and their column labels translated into
   the zero-side frame. Member events and their connectivity are retained.
   Group rows are inserted with a linear matching anchor; original touch
   rows also remain event rows. Directly solved regular midpoint rows
   separate adjacent event positions within a segment. Final row indices
   are resolved after insertion.
5. **Finalize on the rebuilt ItemViews.** Read the nearest regular rows
   on both sides of each event to determine boundary-pair changes and
   charges. Tie rows cannot provide a stable sort order. A column's soft
   charge is `(side_right - side_left)/2`, with side -1 below the M cut and
   +1 above it. Ordinary charges can be -1, 0, or +1. Components involving
   MRs or unreliable directions have `charge=None` and form hard boundaries.
6. **Detect continuum and build the path.** In each ItemView,
   `j_lo = sort_to_item[:, M-1]` and `j_hi = sort_to_item[:, M]`.
   `has_continuum` is true if any row has `j_lo == j_hi`. The same analyzed
   mesh supplies the independent `Mu2Mid` path used for winding integration.

`analyze()` is a one-time operation per instance. Its optional supplied
`continuum_clusters` seed the initial ItemViews; later mesh refinement and
finalization recompute internal clusters. Build a fresh instance to rerun
analysis with different settings.

### 2.2 The interpolated mu2_mid path

`Mu2MidPiece` stores interval endpoints, values, derivatives, and a NumPy-order
polynomial in `theta1 - theta0`. `Mu2Mid` holds these pieces independently
of the ZeroManager. It provides `value_deriv`, `value`, `values_at`, and
compatibility arrays `theta1`, `values`, and `derivs`.

Boundary log-moduli are clipped individually to `±mu2mid.LOGABS_CLAMP`
(default 14) before averaging. A saturated contribution has zero derivative;
a non-finite unsaturated contribution forces linear interpolation. Each
interpolant is split at crossings of the clamp band, and outside-band
pieces become constant. Sorting and crossing detection use raw log-moduli.

`Mu2Mid.breakpoints` returns knots whose endpoint values or derivatives do
not agree within its numerical comparison tolerance. Ordinary C1 knots
are omitted, so a quadrature interval can span several Hermite pieces.
Splitting at derivative jumps avoids extending a global interpolant across
a change of boundary pair or an MR singularity.

### 2.3 Materializing subsets

`winding.detect_crossings_simple` consumes analyzed EventGroups. A component
must cover both boundary sort positions on at least one regular side.
Every real root column in that component yields a PointSubset and a charge
dictionary with `theta1`, `theta2`, `charge`, and `kind`. MR events use this
same channel; no separate MR-only materializer or proximity-based echo
removal discards nearby events.

`continuum_lines.extract_continuum_linesubsets` splits segments at endpoint
rows and events with `changes_boundary=True`. A representative regular row
between cuts determines whether that interval belongs to a continuum item;
the returned line includes both cut rows. Events that leave the boundary
pair unchanged do not split a line. This uses refined event locations rather
than terminating at the last regular row before a sort change.

Line pieces join across adjacent segment boundaries only when their track
is outside the MR cluster. Boundary-changing events also block seam joins.
Interior endpoints share column identities; across the circular seam,
`roots_right[boundary_perm] == roots_left` translates those identities.
A seam-spanning line can have theta1 samples above `2*pi` so its stored
parameter remains monotonic.

### 2.4 Average winding and seed selection

For fixed theta2, the integration loop is
`beta1 = exp(mu1 + i*theta1)`,
`beta2 = exp(mu2_mid(theta1) + i*theta2)`.
`WindingFun` evaluates `Im[f'/f]`, using
`d beta1/d theta1 = i*beta1` and
`d beta2/d theta1 = beta2 * mu2_mid'(theta1)`.
`get_winding_number` integrates this quantity and divides by `2*pi`.

Hard boundaries partition the theta2 circle into regions. Soft boundaries
subdivide regions into intervals. For each region with positive width,
`_pick_seed_theta2` samples interval interiors and chooses the candidate
with the largest mesh-sampled minimum of `|f|`. This avoids small
polynomial denominators in `f'/f`; root distance alone omits the local
polynomial scale. The score is a numerical heuristic, not a lower bound
between mesh rows.

One quadrature supplies a rounded integer winding per region. Soft charges
propagate it to adjacent intervals, and their angular widths give the
average. Coincident boundaries create zero-width intervals. With no hard
boundaries, charges must sum to zero; otherwise `RuntimeError` reports
inconsistent crossing data. Hard boundaries require independent winding
seeds because their unknown charges cannot support propagation.

### 2.5 Bisection and plateau detection

`solve_SGBZ_for_E` expands the mu1 bracket, then uses midpoint bisection.
The winding's flat plateaus can stall false-position methods; midpoint
bisection supplies progress on ordinary bracket updates. Each winding
evaluation builds one Mu2MidZM and reuses it for detection and integration.

At a continuum, direct winding is skipped. Probes at
`mu1 ± scale*continuum_perturb`, with scales from `core.ESCAPE_LADDER`,
resolve its two limits. Strictly opposite signs identify a continuum
boundary and return materialized LineSubsets. Otherwise the solver uses
the relevant perturbed coordinate and its winding to update the search.
An exactly zero edge winding is re-evaluated with its subsets before a
discrete return. Convergence requires an actual winding within `zero_tol`
or a verified continuum boundary; bracket width alone is insufficient.
Expansion and iteration exhaustion raise errors.

For nonempty discrete candidates, `collect_GBZ_subsets` optionally checks
for a neighboring zero plateau. Every point must first have a nearby
neighbor in both theta1 and theta2 on the torus. The probe ladder evaluates
both sides of each step; either side finding a successful non-continuum
solve with zero winding and no subsets establishes a plateau. Probe errors
propagate to the entry point's error handler. Continuum results skip this
check. See [constants.md](constants.md) for scales and live defaults.

## 3. API reference

The following signatures use `None` for values resolved from module
constants at call time. Unknown keywords raise `TypeError`.

```python
from pygbz2d import CharPoly
from pygbz2d.sgbz import Mu2MidZM, collect_GBZ_subsets

poly = CharPoly(coeffs, degs)
zm = Mu2MidZM(poly, E_ref, mu1)
zm.run()  # StepControl and root-tracking settings are available here.
zm.analyze(tie_tol=None, crossing_tol=None, min_direction_deriv=None,
           verbose=False, refine_multi_crossings=True,
           refine_max_rounds=None, refine_safety_factor=None,
           refine_max_subintervals=None, refine_max_total_inserts=None)
```

`build_mu2_mid` on Mu2MidZM is a compatibility alias for `analyze`.
`mu2_mid`, `has_continuum`, `seg_mu2_values`, `seg_mu2_derivs`, and the flat
`mu2_mid_*` arrays become available after analysis. Pair events, EventGroups,
and ItemViews are private analysis data.

```python
def solve_SGBZ_for_E(poly, E_ref, mu1_guess=(-1, 1), *, zero_tol=None,
                     continuum_perturb=None, max_iter=None,
                     continuum_tol=None, crossing_tol=None) -> dict: ...

def collect_GBZ_subsets(coeffs, degs, E_ref, debug_mode=False, *,
                        plateau_check=True, plateau_probe_radius=None,
                        mu1_guess=(-1, 1), zero_tol=None,
                        continuum_perturb=None, max_iter=None,
                        continuum_tol=None, crossing_tol=None) -> GBZResult: ...
```

`solve_SGBZ_for_E` returns `mu1`, `subsets`, `winding`, `is_continuum`,
`_mu1_bracket`, `_winding_bracket`, and `_exit_reason`; continuum boundaries
also include `_w_limits`. Exit reasons are `w_zero`,
`w_zero_continuum_edge`, `left_endpoint_zero`, `right_endpoint_zero`, and
`continuum_boundary`. Iteration exhaustion raises instead of returning a
successful diagnostic record.

`collect_GBZ_subsets` returns the unified GBZResult. Check `success` before
interpreting `index=(0, 0)` or `is_empty` as an exterior energy. Solver and
assembly exceptions become failed results unless `debug_mode=True`;
polynomial construction and signature errors occur outside that handler.
Neither entry point accepts `perc`, `zm_run_kwargs`, or legacy catch-all
options. The default `continuum_perturb` is `core.CONTINUUM_PERTURB = 1e-4`.

Other public helpers:

| Signature | Behavior |
|---|---|
| `detect_continuum_simple(zm, poly, *, continuum_tol=None)` | Builds a fresh Mu2MidZM and returns its continuum flag. |
| `extract_continuum_linesubsets(zm, poly)` | Requires continuum; materializes and joins LineSubsets. |
| `detect_crossings_simple(zm, poly, *, crossing_tol=None)` | Materializes points and charges; crossing_tol is compatibility-only here. |
| `detect_crossings_and_winding(zm, poly, *, crossing_tol=None)` | Ensures one analyzed Mu2MidZM and returns `(points, W)`. |
| `compute_average_winding(zm, poly, charges)` | Computes the angularly weighted mean from charges and the path. |
| `get_winding_number(winding_fun, seg_bounds=None, *, n_seg=1)` | Unrounded quadrature winding; n_seg is used if no breakpoints are supplied. |

`ensure_mu2mid` in `mu2mid.py` reuses an analyzed Mu2MidZM; otherwise it
builds a fresh one from `(poly, E_ref, mu1)` with default analysis settings.
To control crossing refinement, pass an already analyzed Mu2MidZM.

Public data exports include `ItemView`, `Mu2Mid`, `Mu2MidPiece`, `CharPoly`,
`PointSubset`, `LineSubset`, `GBZResult`, and `ConnectedSubset`.

## 4. Relation to amoeba

| Aspect | SGBZ | Amoeba |
|---|---|---|
| Tracked roots | ZeroManager at fixed E and mu1 | Same continuation layer |
| Transverse radius | Piecewise interpolated mu2_mid(theta1) | Constant mu2 |
| Crossing refinement | Representative-pair log-modulus difference, transient Brent solves | Track log-modulus minus mu2, Brent solves retained in the mesh |
| Winding | Loop quadrature and side-charge propagation | Polynomial root counts at angular partition midpoints |
| Continuum | Equal-modulus items spanning the M cut | Tracks with nearly constant log-modulus |
| Line extraction | Intervals cut at boundary-changing events | Recorded whole-segment flat tracks |
| Spectrum | Contained in the amoeba spectrum | Contains the union of strip spectra |

For uniform bands the spectra coincide. The two formulations share root
tracking and data types, but use different crossing and winding algorithms.
