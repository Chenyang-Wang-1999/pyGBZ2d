# Geometry-dependent skin-effect application

Completed `application/geometry-dependent-skin-effect.py` and the accompanying
eight-page Typst note/PDF. Model and boundary conventions follow the supplied
Haldane scripts. No solver-core changes were made. The pre-existing changes in
`playground/Haldane-model-gainloss.py` and `paper/` were preserved.

## Implementation

- BerryPy two-site gain-loss Haldane model; two-site `a1`, `a2` and four-site
  Cartesian `x`, `y` basis choices. Parameters are `(1, 0.5, pi/3, 0.5j, 0)`.
- CLI commands: `demo` (default), `check`, `sweep`, `obc`, `summarize-obc`, `plot`.
- SGBZ scans use the current `collect_GBZ_subsets` signature. They retain complete
  `GBZResult` objects, all line-subset mesh rows, errors, and unfinished entries.
  Atomic checkpoints resume missing entries; `--retry-failed` retries errors.
- OBC matrices are assembled directly from BerryPy hopping lists. Bonds leaving
  the selected set of cells are discarded. This linear-in-bond-count assembly
  avoids expensive construction of a large periodic supercell solely to cut it.
- Right-eigenvector densities are normalized per column before averaging. The
  common bulk-band window is `-2.5 <= Re(E) <= -1.5`. It is a spectral selection,
  not a proof that every selected finite-system state has bulk character.
- Original pickles are preserved. Full new OBC eigenvectors are optional;
  compact summaries are the default. Only trusted pickles should be loaded.

## Validation actually run

Environment: Python 3.11.5, NumPy 1.26.4, SciPy 1.11.4, Matplotlib 3.7.2,
SymPy 1.11.1, Typst 0.15.1. BerryPy is installed locally without discoverable
distribution-version metadata. SSH extraction used remote NumPy 1.22.3.

1. `python application/geometry-dependent-skin-effect.py check`
   - Eight deterministic complex determinant comparisons in each of four bases.
   - Maximum normalized errors: a1 `4.38e-16`, a2 `6.22e-16`, x `8.19e-16`,
     y `3.24e-15`.
   - Independent 3-by-2 finite samples compared with BerryPy's two-stage
     supercell construction and `get_bulk_Hamiltonian_complex((None, None))`.
     Every matrix difference is exactly zero. Coordinate errors <= `4.45e-16`.
2. `python application/geometry-dependent-skin-effect.py demo --workers 2`
   - Six energies per direction: `Re(E) = [-2, 1.5, 5]`, `Im(E) = [0, 0.3]`.
   - All 24 results completed successfully; each direction has three in-spectrum
     and three outside-spectrum energies. At `E=1.5`, y returns two LineSubsets.
   - The a1/a2/x scans each contain 18 subset samples; y contains 388.
   - Normalized polynomial residual `|sum terms| / max(1, sum |terms|)` over
     every retained point/line sample is <= `4.64e-16`.
   - Two 128-site samples: rhombus 8-by-8, rectangle 8-by-4. Selected right
     eigenpair residuals are `8.54e-15` and `8.36e-15`, respectively.
3. Resume command with the identical 3-by-2 energy grid: all four directions
   reported zero pending energies. Completed results were reused.
4. `python -m pytest tests/test_geometry_dependent_example.py -q`: 3 passed.
   Checks invariance to arbitrary eigenvector normalization, the energy-window
   average, rejection of misordered caches, and preservation/retry of failures.
5. Saved-data plotting, Typst compilation, and visual review of all eight PDF
   pages. Numerical validation and source metadata are also in
   `application/Figures/geometry-dependent-skin-effect/validation.json`.

## Reused data

The four local `data/Haldane-gain-loss-{a1,a2,x,y}-SGBZ.pkl` files each contain
40,401 energies. None is flagged failed. In-spectrum counts are 9,991 for each
zigzag direction and 6,656 for y. No points were suppressed or projected onto
unit circles.

Two remote files were read on `myoffice` in
`~/654/research-data/2D_skin_effect/phcpy-free-algorithm/data`:

- `paper-Haldane-gain-loss-OBC-80-80.pkl` (rhombus, 12,800 sites).
- `paper-Haldane-gain-loss-OBC-square-80-40.pkl` (rectangle, 12,800 sites).

Each original is 2,621,849,888 bytes. A copy of the new script was streamed over
SSH to `python3 - summarize-obc ... --output -`, without installing or writing
anything on the server. Binary stdout was saved locally to `obc-rhombus.npz`
(514,980 bytes) and `obc-rectangle.npz` (556,607 bytes) in the figures directory.
They include the full energy spectrum, coordinates, all-state and window-state
density means, IPRs, and one selected state. They do not contain full eigenvector
matrices. The mean-density normalization sums are within `3e-15` of one.
The window contains 5,178 rhombus states and 5,231 rectangle states.

The original OBC tuples do not store model parameters. Parameter provenance is
the supplied companion script's `ALL_PARAMS`, verified on the server. Source
path, size, modification time and extraction time are retained in NPZ metadata.

## Unexpected results and handling

### Initial smoke-test energy selection

The first trial used `Re(E)=[-2,0,5]`, `Im(E)=[0,0.2]`; all selected energies
were outside all four spectra. This was a poor integration-test selection, not
evidence of a solver defect. Existing data were checked before switching to the
validated grid above. Both spectral interior and exterior, plus a continuum
subset, are now exercised. The initial trial is retained only as ignored local
scratch data in `application/data/geometry-dependent-skin-effect/initial-demo-empty-selection`.
No solver tolerance or physical classification was altered.

### Isolated x-SGBZ radius deviation (unresolved)

At `E=1.674+0.2397j`, six subsets in the existing x cache have
`mu1=1.464784145e-5`; the maximum `abs(mu2)` is `1.154522346e-4`.
Fresh computation with the current solver and the saved polynomial reproduces
the result. A fixed-mu diagnostic gives:

| Fixed mu1 | Average winding | Maximum abs(mu2) |
|---|---|---|
| 0 | -1.545827293e-5 | 9.023851039e-5 |
| 1.464784145e-5 | 4.290513206e-9 | 1.154522319e-4 |

This places the discrepancy in the winding/root-extraction calculation, rather
than the plot or pickle migration. It does not establish its full numerical
cause. The record is not a claim that the solver has been repaired. The x-scan
99th percentile of `abs(mu2)` is `4.58e-11`; every outlier remains in the plots
and maximum statistics. A minimal reproducer is:

```python
import pickle
import numpy as np
from pygbz2d.sgbz import collect_GBZ_subsets
with open("data/Haldane-gain-loss-x-SGBZ.pkl", "rb") as stream:
    data = pickle.load(stream)
result = collect_GBZ_subsets(
    np.asarray(data["coeffs"], complex), np.asarray(data["degs"], int),
    1.674 + 0.2397j,
)
print(result)
```

Run against this checkout's `src` (or its editable installation). Fixed-mu
diagnostics used `_evaluate_winding` with the module's default continuum and
crossing tolerances. These private calls are diagnostic only, not part of the
application's production code.

### Warnings during the small run

BerryPy emits a SymPy deprecation warning about non-Expr objects in a Matrix.
The y real-energy calculation also emits a runtime warning from `pairwise.py`
in the subtraction of root-derived values. All returned samples are finite and
their polynomial residuals are near machine precision. Warnings are not
suppressed, and the solver implementation/tolerances were not changed to hide
them. The `success` flag alone should not be read as proof of all numerical
accuracy; the retained x outlier illustrates this limitation.

## Interpretation

The armchair SGBZ has transverse `abs(mu2)` as large as 0.590, contrasted with
the nearly unit-radius zigzag cases. Its finite rectangular realization shows
strong bulk-window intensity near the armchair boundaries, whereas the rhombus
retains substantial bulk intensity. Corner enhancement in the latter is
retained and described rather than equated with extensive bulk skin modes.
The note distinguishes this numerical geometry comparison from an analytic
proof of exact SGBZ/BZ equality or a general statement about every polygon.
