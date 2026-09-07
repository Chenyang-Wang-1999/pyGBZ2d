# Numerical figures: inputs needed

The main deck contains theory and explicitly labeled vector schematics. It contains no invented data graphs. The analytic Hatano–Nelson equations are a benchmark specification, not a computed result.

To make research figures for your preferred model, please supply or choose:

1. **Model and basis:** `coeffs` and `degs` for `f(E,beta1,beta2)`, hopping parameters, the lattice basis, and the energy unit. Store the arrays in `model.npz`, with complex `coeffs` of shape `(n,)` and integer `degs` of shape `(n,3)` in the order `(E,beta1,beta2)`. Negative Bloch exponents are allowed. The scripts do not silently rotate the basis.
2. **A diagnostic slice:** one complex reference energy and one fixed `mu1` for the root-track plot. Choose a slice that shows the event you want to discuss. The automatically solved radii used in the subset plot are generally different from this fixed diagnostic radius.
3. **An energy grid:** minimum, maximum, and number of points for both Re E and Im E. Both methods must use the same grid for a meaningful spectral comparison. Include a real slice or sufficiently fine boundary sampling when testing a spectrum supported on a line.
4. **A subset energy:** choose a computed grid point with the desired point/line structure. `--energy-index` refers to the unique energies in the exported sweep, with Re E varying fastest.

An additional finite-size OBC comparison requires the actual Hamiltonian builder, shape, boundary termination, and sizes. A Chern-number figure requires left/right eigenvectors, band selection, a gap condition, and a validated oriented mesh. The determinant alone is insufficient for those quantities. These figures are not fabricated or implemented by the current plotting script.

## Included scripts

`collect_slide_data.py` computes data with the checkout's NumPy backend. It exports JSON containing the exact polynomial arrays, command, source revision, library versions, root tracks, subsets, per-energy success/error states, warnings, runtimes, and normalized polynomial residuals. Non-finite values become JSON null rather than unsupported numeric literals. Public solvers use their current repository defaults with the plateau check enabled.

`plot_slide_data.py` consumes that JSON and creates three figure pairs:

| Output | Quantity | Treatment of missing/failed data |
| --- | --- | --- |
| `root-tracks.svg`, `.pdf` | `ln|beta2|` versus theta1 at fixed E and mu1 | Singular samples produce gaps; a failed track run omits this figure |
| `subsets.svg`, `.pdf` | Returned point/line subsets on the angle torus | Empty and failed methods appear explicitly in the legend |
| `spectrum.svg`, `.pdf` | Computed grid points classified by method | Failed runs have red crosses; outside points stay visible; no filled interpolation |

The root-plot colors identify columns within each segment. They do not assert globally fixed labels around the periodic seam. The subset plot splits display segments at wrapped-angle jumps. This does not change the unofficial debug plot excluded by review issue #2.

The plotting script writes `plot-provenance.json` with the input SHA-256 and calculation metadata. It never computes additional model data. Normalized residuals in the JSON use `|sum terms| / sum |terms|`; small residuals alone do not prove correct GBZ classification.

## Runnable examples from the repository root

Real Hatano–Nelson benchmark, J1=J2=1, gamma1=0.2, gamma2=0.3, axis-aligned basis:

```powershell
python paper/slides/collect_slide_data.py --model hn2d --energy 1+0j --mu1 0.25 --sweep-re -4.5 4.5 91 --sweep-im 0 0 1 --output paper/slides/hn2d-results.json
python paper/slides/plot_slide_data.py paper/slides/hn2d-results.json --energy-index 55
python paper/slides/build_slides.py --with-data
```

This is an explicit calculation recipe, not a supplied numerical result. The benchmark parameters are implemented in the collector, so it does not depend on pytest or an unpublished model package. A full sweep can take appreciable time. Begin with a small grid to check your chosen model.

For your research model, first save your actual arrays:

```python
import numpy as np
np.savez("model.npz", coeffs=coeffs, degs=degs)
```

Then select your actual diagnostic slice and grid:

```powershell
python paper/slides/collect_slide_data.py --polynomial model.npz --energy 1.212+0j --mu1 0.135328598265 --sweep-re -3 3 61 --sweep-im -1 1 21 --output paper/slides/model-results.json
python paper/slides/plot_slide_data.py paper/slides/model-results.json --energy-index 640
```

The second command block illustrates CLI syntax with an example grid. Those energies and radii are not automatically appropriate for an arbitrary model. Choose them before launching the calculation. All energies use the unit implied by your polynomial coefficients.

## Convergence and performance figures

The current scripts plot the three observables above. A publishable convergence study additionally needs a specified sequence of solver tolerances, a reference solution or error definition, and identical model/grid settings across runs. A timing comparison needs hardware, threading, backend, warm-up policy, and repeated runs. The deck states these as needed evidence and makes no speedup or convergence-rate claim.

