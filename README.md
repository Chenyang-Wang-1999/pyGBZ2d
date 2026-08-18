# brute-force-non-hermitian

Non-Hermitian skin effect computation for 2D tight-binding models — brute-force polynomial root-solving approaches.

This package implements two complementary formulations for determining the generalized Brillouin zone (GBZ) and energy spectrum of 2D non-Hermitian lattice systems:

| Module | Approach | Key Object |
|--------|----------|------------|
| `brute_force_SGBZ` | SGBZ / average major-axis winding | PMGBZ points, average winding $W(E, \mu_1)$ |
| `brute_force_amoeba` | Amoeba / Ronkin function | Average winding numbers, Ronkin minimum $(\mu_1, \mu_2)$ |

The characteristic polynomial $f(E, \beta_1, \beta_2) = \det[E - h(\beta_1, \beta_2)]$ of a 2D non-Hermitian tight-binding model is a Laurent polynomial in $\beta_j = e^{\mu_j + i\theta_j}$. Both modules solve for the non-Bloch decay factors $(\mu_1, \mu_2)$ that satisfy the GBZ condition, but through different mathematical routes.

## Prerequisites

- Python 3.9+
- **[poly_tools](https://atomgit.com/wangchenyang99/PolyTools)** — C-extension library for Laurent polynomial manipulation. Install separately:

```bash
git clone https://atomgit.com/wangchenyang99/PolyTools.git
cd PolyTools/pybind
make _poly_tools_cc.cpython-39-x86_64-linux-gnu.so   # adjust suffix to your Python version
cp -r ../python/poly_tools /path/to/your/workdir/
```

- (Optional) **BerryPy** — only needed for `MatWindingFun` in `brute_force_SGBZ.winding`.

## Installation

```bash
git clone <repo-url>
cd brute-force-non-hermitian
pip install -e .
```

## Quickstart

### Characteristic Polynomial Format

Polynomials use triplet encoding `[E_exponent, beta1_exponent, beta2_exponent]`. For example, the square lattice nearest-neighbor model $f = E - \beta_1 - \beta_1^{-1} - \beta_2 - \beta_2^{-1}$:

```python
import numpy as np
import poly_tools as pt

coeffs = pt.CScalarVec(np.array([1, -1, -1, -1, -1], dtype=complex))
degs = pt.CLaurentIndexVec(np.array([
    1, 0, 0,    # +E
    0, 1, 0,    # -beta1
    0,-1, 0,    # -beta1^{-1}
    0, 0, 1,    # -beta2
    0, 0,-1,    # -beta2^{-1}
], dtype=np.int32))
char_poly = pt.CLaurent(3)
char_poly.set_Laurent_by_terms(coeffs, degs)
```

### SGBZ Solver

Find the critical $\mu_1$ where the average major-axis winding number vanishes:

```python
from brute_force_SGBZ import collect_GBZ_subsets, CharPoly

# Build characteristic polynomial
poly = CharPoly(coeffs, degs)

# Check spectrum membership for a reference energy
gbz = collect_GBZ_subsets(coeffs, degs, E_ref=1.0 + 0j)
print(f"In spectrum: {not gbz.is_empty}, subsets: {len(gbz.subsets)}")
```

### Amoeba Solver

Find the Ronkin function minimum $(\mu_1, \mu_2)$ where both average windings vanish:

```python
from brute_force_amoeba import bisect_amoeba_ronkin_min

result = bisect_amoeba_ronkin_min(char_poly, 1.0 + 0j, -0.5, 0.5, -2.0, 2.0)
print(f"mu1 = {result['mu1']:.6f}, mu2 = {result['mu2']:.6f}")
```

## API Overview

### `brute_force_SGBZ`

| Function | Description |
|----------|-------------|
| `collect_GBZ_subsets(coeffs, degs, E_ref)` | Main entry point — check SGBZ condition for reference energy |
| `solve_SGBZ_for_E(poly, E_ref)` | Locate $\mu_1$ where average winding vanishes |
| `detect_continuum_simple(zm, poly)` | Continuum detection (presence only) |
| `detect_crossings_simple(zm, poly)` | Crossing detection + charge classification |
| `compute_average_winding(zm, poly, M, charges)` | Compute average major-axis winding number |
| `CharPoly(coeffs, degs)` | Characteristic polynomial wrapper |

### `brute_force_amoeba`

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

## Running Tests

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

## Running Demos

```bash
python demos/demo_amoeba.py
python demos/demo_solver.py
```

## Spectrum Inclusion Relation

$$\sigma_{\text{Amoeba}} \supset \bigcup_j \sigma_{\text{SGBZ}, j}$$

The amoeba spectrum is a superset of the union of SGBZ spectra. For uniform bands the two are equal.

## References

The SGBZ formulation is based on the average major-axis winding number approach for 2D non-Hermitian systems. The amoeba formulation uses the Ronkin function of Laurent polynomials. See [doc/](doc/) for detailed theoretical background.

## License

MIT — see [LICENSE](LICENSE).
