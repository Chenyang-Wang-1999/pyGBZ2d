# pyGBZ2d — brute-force-non-hermitian

Non-Hermitian skin effect computation for 2D tight-binding models — brute-force polynomial root-solving approaches.

This package implements two complementary formulations for determining the generalized Brillouin zone (GBZ) and energy spectrum of 2D non-Hermitian lattice systems:

| Subpackage | Approach | Key Object |
|--------|----------|------------|
| `pygbz2d.sgbz` | SGBZ / average major-axis winding | PMGBZ points, average winding $W(E, \mu_1)$ |
| `pygbz2d.amoeba` | Amoeba / Ronkin function | Average winding numbers, Ronkin minimum $(\mu_1, \mu_2)$ |

The characteristic polynomial $f(E, \beta_1, \beta_2) = \det[E - h(\beta_1, \beta_2)]$ of a 2D non-Hermitian tight-binding model is a Laurent polynomial in $\beta_j = e^{\mu_j + i\theta_j}$. Both modules solve for the non-Bloch decay factors $(\mu_1, \mu_2)$ that satisfy the GBZ condition, but through different mathematical routes.

## Installation

```bash
pip install .            # from a clone, or: pip install -e . for development
```

Runtime dependencies are **numpy + scipy only** — the package is fully functional out of the box.

### Optional: poly_tools acceleration

Polynomial evaluation runs on a pluggable backend. The default auto-selection uses the compiled [poly_tools](https://atomgit.com/wangchenyang99/PolyTools) C++ extension when it is importable, and otherwise falls back to the built-in pure-numpy backend (numerically equivalent to ~1e-14; ~25% slower on root solving, with a one-time warning):

```bash
git clone https://atomgit.com/wangchenyang99/PolyTools.git
cd PolyTools/pybind
make _poly_tools_cc.cpython-39-x86_64-linux-gnu.so   # adjust suffix to your Python version
cp -r ../python/poly_tools /path/to/site-packages/    # or anywhere on PYTHONPATH
```

Backend selection (first match wins):

```python
CharPoly(coeffs, degs, backend="numpy")        # explicit, per-polynomial
# or the POLY_BACKEND environment variable: "poly_tools" | "numpy"
# or backend=None (default): poly_tools if importable, else numpy fallback
```

A custom backend is any class satisfying the `LaurentProtocol` in `pygbz2d/backend.py` (eval / derivative / partial_terms_1d / num_max_degrees).

- (Optional) **BerryPy** — used by the Haldane playground script and Haldane counterexample tests; those tests skip automatically when BerryPy is absent.

## Quickstart

### Characteristic Polynomial Format

Polynomials use triplet encoding `[E_exponent, beta1_exponent, beta2_exponent]`. For example, the square lattice nearest-neighbor model $f = E - \beta_1 - \beta_1^{-1} - \beta_2 - \beta_2^{-1}$:

```python
import numpy as np
from pygbz2d import CharPoly

coeffs = np.array([1, -1, -1, -1, -1], dtype=complex)
degs = np.array([
    [1, 0, 0],    # +E
    [0, 1, 0],    # -beta1
    [0,-1, 0],    # -beta1^{-1}
    [0, 0, 1],    # -beta2
    [0, 0,-1],    # -beta2^{-1}
], dtype=int)
poly = CharPoly(coeffs, degs)
```

### SGBZ Solver

Find the critical $\mu_1$ where the average major-axis winding number vanishes:

```python
from pygbz2d.sgbz import collect_GBZ_subsets

# Check spectrum membership for a reference energy
gbz = collect_GBZ_subsets(coeffs, degs, E_ref=1.0 + 0j)
print(f"In spectrum: {not gbz.is_empty}, subsets: {len(gbz.subsets)}")
```

### Amoeba Solver

Find the Ronkin function minimum $(\mu_1, \mu_2)$ where both average windings vanish:

```python
from pygbz2d.amoeba import bisect_amoeba_ronkin_min

result = bisect_amoeba_ronkin_min(poly, 1.0 + 0j, -0.5, 0.5, -2.0, 2.0)
print(f"mu1 = {result['mu1']:.6f}, mu2 = {result['mu2']:.6f}")
```

## Customizing Numerical Constants

Every numerical constant lives in the module that consumes it (RK45
stepper knobs in `continuation/arclength.py`, crossing tolerances in
`sgbz/pairwise.py`, ...); the seven cross-package constants live in
`pygbz2d.core`.  The complete per-module reference with defaults and
tuning guidance is [doc/constants.md](doc/constants.md).

The customization model has exactly two layers:

```python
import pygbz2d as bz

# 1) tune ONE call — a plain keyword argument (misspelled names raise
#    TypeError; there is no catch-all options dict)
gbz = bz.sgbz.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, continuum_tol=1e-8)

# 2) tune the WHOLE process — assign the module constant; takes effect
#    immediately, on the next read, everywhere (including running loops)
from pygbz2d.sgbz import pairwise
pairwise.CROSSING_TOL = 1e-12
```

Process-scope note: assignments propagate to forked worker processes only
if made **before** the pool is created.

## API Overview

### `pygbz2d.sgbz`

| Function | Description |
|----------|-------------|
| `collect_GBZ_subsets(coeffs, degs, E_ref)` | Main entry point — check SGBZ condition for reference energy |
| `solve_SGBZ_for_E(poly, E_ref)` | Locate $\mu_1$ where average winding vanishes |
| `Mu2MidZM(poly, E_ref, mu1)` | ZeroManager + ItemView analysis + pairwise crossing detection + μ₂_mid path |
| `detect_continuum_simple(zm, poly)` | Continuum detection (presence only) |
| `detect_crossings_simple(zm, poly)` | Crossing detection + charge classification |
| `compute_average_winding(zm, poly, charges)` | Compute average major-axis winding number |
| `CharPoly(coeffs, degs, backend=None)` | Characteristic polynomial wrapper |

`Mu2MidZM.analyze()` runs a pre-crossing mesh refinement before the pairwise scan: intervals whose cubic-Hermite interpolants predict two or more crossings are sub-divided (disable with `refine_multi_crossings=False`).

### `pygbz2d.amoeba`

| Function | Description |
|----------|-------------|
| `collect_GBZ_subsets(coeffs, degs, E_ref, ...)` | Main entry point — check amoeba condition for reference energy |
| `bisect_amoeba_ronkin_min(char_poly, E_ref, ...)` | Find $(\mu_1, \mu_2)$ where both average windings vanish |
| `AmoebaZeroManager(poly, E_ref, mu1)` | Adaptive $\beta_2$ root tracks via `continuation.ZeroManager` |
| `extract_amoeba_subsets(zm, poly, E, mu1, mu2)` | Extract GBZ subsets from track crossings of $\ln|\beta_2| = \mu_2$ |
| `amoeba_windings(zm, poly, E, mu1, mu2)` | Average windings (a1/a2) and refined zero list |

## Documentation

Detailed documentation is available in the `doc/` directory:

- [doc/SGBZ.md](doc/SGBZ.md) — SGBZ formulation, architecture, API reference
- [doc/amoeba.md](doc/amoeba.md) — Amoeba formulation, algorithm design, API reference
- [doc/continuation.md](doc/continuation.md) — continuation module (arclength / multiple roots / ZeroManager)

## Running Tests

```bash
pip install -e .[dev]      # or: pip install pytest
pytest                     # slow tests are skipped by default
pytest --run-slow          # include BerryPy-dependent slow tests
POLY_BACKEND=numpy pytest   # full suite on the pure-numpy backend
```

## Playground Scripts

```bash
python playground/demo_unified.py               # GBZResult API for both modules
python playground/Haldane-model-gainloss.py     # Haldane sweeps / plots / failed-E recompute
python playground/demo_zero_manager.py          # ZeroManager root tracking
```

See `playground/` for the full set of runnable scripts (unofficial, not part of the package).

## Spectrum Inclusion Relation

$$\sigma_{\text{Amoeba}} \supset \bigcup_j \sigma_{\text{SGBZ}, j}$$

The amoeba spectrum is a superset of the union of SGBZ spectra. For uniform bands the two are equal.

## References

The SGBZ formulation is based on the average major-axis winding number approach for 2D non-Hermitian systems. The amoeba formulation uses the Ronkin function of Laurent polynomials. See [doc/](doc/) for detailed theoretical background.

## License

MIT — see [LICENSE](LICENSE).
