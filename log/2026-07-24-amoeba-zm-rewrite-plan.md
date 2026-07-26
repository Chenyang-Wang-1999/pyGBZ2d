# Plan: full §3 rewrite of `brute_force_amoeba` on top of `continuation.ZeroManager`

Confirmed scope: delete old tracks/winding-source, rewrite bisection winding source to consume ZM zeros, add `AmoebaZeroManager` + 3-mode `extract_amoeba_subsets` + `amoeba_windings`, wire `collect_GBZ_subsets` (reuse the bisection's `_zm`), update tests + the one external demo line. Run a *real* ZeroManager. Plateau probe logic kept, fed new windings. **No separate MR-candidate emission** — after snap + boundary dedup, MR points are already produced by the normal crossing/continuum algorithm.

## 0. What stays / what goes (per spec §6)

**Keep untouched**
- `continuation/` (ZeroManager, MultipleRootInfo, snap_clusters_to_mean, …).
- `gbz_types.py` (CharPoly, PointSubset, LineSubset, GBZResult, utils). Format already matches §1.
- `_find_exact_crossing`, `_get_average_winding_from_zeros`, `_compute_zero_dtheta1_dmu2` in `ronkin_winding.py`.
- `_probe_zero_plateau_near_mu1`, `_check_zeros_are_clustered`, `_is_zero_plateau_probe` in `amoeba.py`.

**Delete**
- `brute_force_amoeba/tracks.py` (whole file).
- `ronkin_winding.py`: `_compute_winding_from_tracks`, `_compute_crossings_and_winding`, `get_a2_average_winding`, `get_a1_average_winding`.
- `amoeba.py`: `_extract_continuum_intervals`, `_detect_crossings_outside_continuum`.

**Rewrite**
- `bisect.py`: `_refine_and_correct`, `_find_mu2_for_w2_zero`, `_resolve_continuum`, `bisect_amoeba_ronkin_min` — winding source swapped to ZM zeros.
- `amoeba.py::collect_GBZ_subsets` — uses AmoebaZeroManager + 3-mode extractor; reuses `_zm` from the bisection result.
- `__init__.py` — drop deleted exports, add new ones.

## 1. New file `brute_force_amoeba/zm_manager.py`

### `class AmoebaZeroManager(ZeroManager)`
Subclass that, after `.run()`, caches per-segment `logabs = log|tracked_roots|` so the extractor and winding code never recompute logs.

```python
class AmoebaZeroManager(ZeroManager):
    """ZeroManager + cached ln|beta2| per segment, for amoeba extraction/winding."""
    seg_logabs: list[np.ndarray]   # seg_logabs[s] = log|segments[s].tracked_roots|

    def __init__(self, poly, E_ref, mu1):
        super().__init__(poly, E_ref, mu1)
        self.seg_logabs = []

    def run(self, **kwargs):
        super().run(**kwargs)
        self.seg_logabs = [np.log(np.abs(seg.tracked_roots)) for seg in self.segments]
```

## 2. New file `brute_force_amoeba/zm_extract.py`

### `extract_amoeba_subsets(zm, poly, E, mu1, mu2, *, mode, tol, frac, snap_tol)`
Productionized 3-mode version of the demo extractor. `mode ∈ {'coarse','fine','solve'}`.

Detection (shared, vectorized, index-based — same as the demo):
- continuum: `frac_in_band = mean(|logabs-mu2|<tol, axis=0); is_cont = frac_in_band>frac` → one `LineSubset(theta1_arr=th, beta2_arr=tr[:,j])` per continuum track. In **all** modes a continuum detection is recorded as a `LineSubset` (no refinement; "detect → return").
- discrete crossings on non-continuum tracks: record hits as `(seg_idx, i, j, kind)` with `kind ∈ {'zero','cross'}`:
  - `'zero'`: `d == 0` exact touch at sample i.
  - `'cross'`: sign change `d[:-1]*d[1:] < 0` over `[i,i+1]`.

Post-loop materialization + boundary filtering (rules from the demo):
- **Rule 1 (continuum-endpoint snap):** drop any crossing whose θ₁ is within `snap_tol` of a continuum LineSubset endpoint.
- **Rule 2 (d==0 boundary dedup):** for `'zero'` hits, dedup by rounded θ₁ (shared MR endpoint → same θ₁). Interior `'zero'` kept once.

Mode differences — ONLY how a crossing's `(theta1, theta2)` is finalized for a `PointSubset`:
- **`coarse`** (refine=False): `PointSubset` from linear interpolation only.
  - `'cross'`: linear-interp `(t1, b2)` → `PointSubset(E, exp(mu1+i t1), b2)`.
  - `'zero'`: `(th[i], tr[i,j])` → `PointSubset`.
- **`fine`** (refine=True): one `_find_exact_crossing` refinement step, fall back to linear on failure.
  - `'cross'`: refine `(t1, angle(b2))`; on `None`, use linear.
  - `'zero'`: refine `(th[i], angle(tr[i,j]))`; on `None`, use the sample.
- **`solve`**: same point finalization as `fine` (refined). Used by `collect_GBZ_subsets` to produce the full GBZ output; reuses the bisection's `_zm`. No extra MR-candidate emission (per user: snap + dedup already yields those points).

LineSubsets identical across modes. Return `list[PointSubset|LineSubset]`.

### `amoeba_windings(zm, poly, E, mu1, mu2, *, tol, frac, refine)`
Replaces `_compute_winding_from_tracks`. Returns `(w2, zeros, has_continuum, dW_dmu2)`:
- Detect continuum + crossings from `zm.seg_logabs` (vectorized, same detection as extractor). `has_continuum=True` when any track is continuum at this mu2 → return `(None, None, True, None)` (mirrors old contract).
- `zeros = [(theta1, theta2, jump), ...]`:
  - `jump = +1 if logabs[i,j] < mu2 else -1`.
  - θ₁/θ₂ via linear interp (`refine=False`) for the bisection loop, or `_find_exact_crossing` (`refine=True`) for the final/Newton step.
  - `'zero'` hits: at `d==0`, the point is exactly on mu2 → θ₁=th[i], θ₂=angle(tr[i,j]).
- `w2 = _get_average_winding_from_zeros(poly, E, mu1, mu2, zeros, direction=2)`.
- `dW_dmu2` (only when `refine=True`): `Σ jump·θ₁_dot /(2π)` via `_compute_zero_dtheta1_dmu2`.
- w1 helper: `amoeba_w1_from_zeros(zm, poly, E, mu1, mu2, zeros)` = `_get_average_winding_from_zeros(..., direction=1)`.

Note: winding evaluation must NOT apply Rule 1/Rule 2 filtering — the winding integral needs every crossing (the snap/dedup rules are a subset-extraction concern, not a winding concern). So `amoeba_windings` uses raw crossing detection; only `extract_amoeba_subsets` applies the boundary rules.

## 3. Rewrite `bisect.py`

Structure (outer μ₁ bisection, inner μ₂ bisection, continuum resolution) preserved; only the winding source changes. Each `(E, μ₁)` builds **one** `AmoebaZeroManager` and reuses it across all μ₂ evaluations — the caching role the old `tracks` dict played.

- `_find_mu2_for_w2_zero(..., _zm=None)`: accept a pre-built `AmoebaZeroManager` (replaces `_root_tracks`). `_winding_at(mu2, refine)` → `amoeba_windings(zm, poly, E, mu1, mu2, refine=refine)`.
- `_refine_and_correct`: call `amoeba_windings(..., refine=True)`; reuse its `dW_dmu2` for the Newton step. Continuum short-circuit unchanged.
- `_resolve_continuum`: w2 ±ε limits via `amoeba_windings(zm, ..., mu2±ε, refine=False)`; w1 limits via re-running the inner μ₂ bisection at μ₁±ε (each building its own ZM).
- `bisect_amoeba_ronkin_min`: build `zm = AmoebaZeroManager(poly, E_ref, mu1)` once per μ₁; pass into `_find_mu2_for_w2_zero`. Return dict keeps `mu1, mu2, zeros, is_continuum, _mu1_bracket, _w1_bracket, _w1_area`, and adds `_zm` (the ZM built at the solved μ₁, for `collect_GBZ_subsets` to reuse).

Continuum contract preserved: `amoeba_windings` returns `has_continuum=True` when a track is in the band over the whole segment (frac>frac_thresh).

## 4. Rewrite `amoeba.py::collect_GBZ_subsets`

```python
def collect_GBZ_subsets(coeffs, degs, E_ref, perc, debug_mode=False, **options):
    char_poly = CharPoly(coeffs, degs)
    solver_options = dict(options)
    plateau_check = solver_options.pop("plateau_check", True)
    ...  # same option extraction
    try:
        amoeba_res = bisect_amoeba_ronkin_min(char_poly, E_ref, **solver_options)
        zm = amoeba_res["_zm"]            # reuse the bisection's ZM (no 2nd run)
        mu1, mu2 = amoeba_res["mu1"], amoeba_res["mu2"]
        subsets = extract_amoeba_subsets(zm, char_poly, E_ref, mu1, mu2, mode="solve")
        # plateau pre-check unchanged in structure, fed new windings:
        #   w1_area from amoeba_res["_w1_area"]; w2 via _get_average_winding_from_zeros(direction=2).
        ...
        return GBZResult(E_ref=E_ref, subsets=subsets, index=(n_0d, n_1d))
    except Exception as e:
        ...
```

The non-continuum `else` branch (zeros from `amoeba_res["zeros"]`) is folded into `extract_amoeba_subsets(mode='solve')` — the extractor's crossings *are* the zeros. `_probe_zero_plateau_near_mu1` called unchanged (it internally calls `_find_mu2_for_w2_zero` + `_get_average_winding_from_zeros`, both retained/rewritten compatibly).

## 5. `__init__.py`

Exports after rewrite:
- Keep: `CharPoly, PointSubset, LineSubset, GBZResult, ConnectedSubset, bisect_amoeba_ronkin_min, collect_GBZ_subsets`.
- Add: `AmoebaZeroManager, extract_amoeba_subsets, amoeba_windings`.
- Remove: `get_hungarian_sorted_roots, get_a2_average_winding, get_a1_average_winding, _compute_root_tracks`.

## 6. Tests + demos

- `tests/test_amoeba.py`:
  - `TestAmoeba` (calls `collect_GBZ_subsets`) — should pass once the rewrite is correct.
  - `test_all_exports` — update expected list to new exports.
  - `TestResolveContinuum` — rewrite to build an `AmoebaZeroManager` and call the rewritten `_resolve_continuum(zm, ...)` instead of `(tracks, ...)`. Same three cases / assertions.
  - `TestPlateauEdge` — uses `bisect_amoeba_ronkin_min` + `_get_average_winding_from_zeros` + `_check_zeros_are_clustered`; retained/rewritten compatibly. Verify pass.
- `demos/replication-ZWang.py:145`: replace `bfa._compute_root_tracks(char_poly, E_ref, 0.0)` with `AmoebaZeroManager(char_poly, E_ref, 0.0); zm.run(...)`; use `zm.segments[0].tracked_roots` / `theta1_arr` for the 3D plot.
- `demos/demo_zm_gbz.py`: leave as-is (standalone reference prototype; not required this step).

## 7. Implementation order
1. `zm_manager.py` (AmoebaZeroManager).
2. `zm_extract.py` (extract_amoeba_subsets 3-mode + amoeba_windings) — port from the working demo.
3. `ronkin_winding.py` — delete the 4 mesh functions, keep the 3 retained.
4. `bisect.py` — rewrite winding source.
5. `amoeba.py` — rewrite collect_GBZ_subsets; delete `_extract_continuum_intervals`/`_detect_crossings_outside_continuum`.
6. `__init__.py` — exports.
7. `tests/test_amoeba.py` + `demos/replication-ZWang.py` — adapt.
8. Run tests + demos.

## 8. Verification
- `pytest tests/test_amoeba.py` (all classes) green.
- `pytest tests/test_zero_manager.py tests/test_continuation.py` green (untouched).
- `demos/demo_zm_gbz.py` still runs (unchanged).
- `demos/replication-ZWang.py` runs with the new ZM source.
- HN model E=0/1/2 continuum intervals match analytic `cos θ₁ ≥ E/2−1` (quick check via new `collect_GBZ_subsets`).
