# TODO: warm-start ZeroManager across μ₁ bisection probes

## Motivation (profiled 2026-08-18)

`solve_SGBZ_for_E` evaluates W(E_ref, μ₁) by bisecting over μ₁, and every
probe builds a **fresh** `Mu2MidZM` from scratch: `_evaluate_winding`
(`brute_force_SGBZ/sgbz_solver.py`) constructs the instance, runs
`ZeroManager.run()` (θ₁ = 0 → 2π arclength integration from a cold start),
and then `analyze()`. Nothing is carried between probes, even though
adjacent bisection probes (e.g. μ₁ = 0.195312 and 0.199219) leave the β₂
root topology nearly unchanged.

Measured cost of this (see `diagnostics/profile_sgbz_report.md`):

| workload | ZM rebuilds | avg mesh rows / build | wall time |
|---|---|---|---|
| HN2D chain-pair, E in spectrum | **33** | 202 | ~35 s |
| HN2D, E outside spectrum | 2 | 66 | ~0.8 s |
| Haldane 2×2 supercell, E=1.212 | 5 | 175 | ~22 s |

The K=2 toy polynomial spends ~35 s almost entirely on
`n_probes × cold_rebuild`. After the 2026-08-18 refinements (amplitude
guard + per-segment incremental rescan in `pairwise.py`), the remaining
cost of each probe is dominated by the integration itself — which is
exactly what a warm start would skip.

## Proposed capability

Give the ZeroManager a way to be built **on an inherited mesh** instead of
by cold integration. Pipeline for a new probe at μ₁_new, seeded with the
mesh of μ₁_prev:

1. **Re-solve on the inherited grid** — for every segment row θ₁ᵢ, call
   `CharPoly.solve_roots_1d` at (E_ref, β₁(μ₁_new, θ₁ᵢ)).
2. **Global Hungarian match onto the old track frame** — match
   re-solved roots against the previous probe's `tracked_roots` row by row
   (chordal cost, the existing `hungarian_match_indices`), preserving
   track identity across the probe.
3. **Adaptive refinement** — walk the inherited mesh and refine where the
   prediction error (old-track Hermite prediction vs re-solved roots)
   exceeds tolerance; coarsen or re-integrate intervals that became
   trivial. This reuses the `arclength.py` step controller philosophy.
4. **MR detection / refinement as usual** — the existing point/interval
   triggers run inside step 3, so multiple roots that appear, drift, or
   annihilate between probes are still caught and recorded.

## Design notes / risks

- The output contract must be identical to `run()`: `segments`
  (`left_mr`/`right_mr`, `boundary_perm`, `has_boundary_mr`,
  `tangents` incl. inf/nan sentinels), because `Mu2MidZM.analyze`,
  `pairwise`, and `continuum_lines` consume exactly that. The warm start
  is a *builder alternative*, not a new topology representation.
- MR positions drift with μ₁. Near an inherited MR the old rows are
  badly conditioned (√-branch divergence), so step 3 must be allowed to
  discard and re-integrate those intervals rather than trust the match.
- Probes that jump across a continuum band (the ±ε perturbation probes of
  `_resolve_continuum_winding`) may change topology discontinuously. The
  warm start needs a validity check (e.g. match cost / refined-fraction
  threshold) with a **fallback to the current cold build**, so worst case
  degrades to today's behaviour.
- `run()` is single-shot by design (it appends). Warm start therefore
  means: new ZeroManager instance + inject mesh + build, never a second
  `run()` on the same instance.
- Track-identity corruption is the main correctness risk: a mis-match
  silently injects spurious modulus crossings into the pairwise scan.
  Gate the change behind the existing `test_counterexamples.py` /
  `test_regressions.py` suites and add a probe-to-probe track-identity
  regression test.

## Expected payoff

For HN2D-like cases (many probes, cheap polynomial) the rebuild count is
the dominant multiplier: 33 cold integrations → 1 cold + 32
re-solve+match passes over a ~200-row mesh. The re-solve is one
`solve_roots_1d` per row (the same cost the integrator pays per accepted
step) but without step rejection, MR refinement, and tangent re-computation
overhead — realistically a 2–4× end-to-end speedup on the E-in-spectrum
path, more once combined with caching tangents between probes.
