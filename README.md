# brute-force-non-hermitian

Non-Hermitian skin effect computation for 2D tight-binding models — brute-force polynomial root-solving approaches.

This package implements two complementary formulations for determining the generalized Brillouin zone (GBZ) and energy spectrum of 2D non-Hermitian lattice systems:

| Module | Approach | Key Object |
|--------|----------|------------|
| `brute_force_SGBZ` | SGBZ / strip winding number | PMGBZ points, strip winding $W(E, \mu_1)$ |
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

Find the critical $\mu_1$ where the strip winding number vanishes:

```python
from brute_force_SGBZ import SGBZSolver

solver = SGBZSolver(char_poly)
mu1, pmgbz_points = solver.solve_for_E(1.0 + 0j)
print(f"mu1 = {mu1:.6f}, PMGBZ points: {len(pmgbz_points)}")
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

| Function / Class | Description |
|-----------------|-------------|
| `SGBZSolver(char_poly)` | Main SGBZ solver. `solve_for_E(E_ref)` returns `(mu1, PMGBZ_points)` |
| `SGBZChecker(char_poly)` | Check if a given $(E, \mu_1)$ satisfies the SGBZ condition |
| `get_strip_winding(poly_diff, E_ref, mu1)` | Compute strip winding number $W(E, \mu_1)$ |
| `get_roots_and_PMGBZ(poly_diff, E_ref, mu1)` | Solve roots and detect PMGBZ points |
| `get_loop_winding(poly_diff, E_ref, mu1, mu2_fun, theta2)` | Single-loop winding number |
| `PolyDiffContext(char_poly)` | Polynomial + precomputed partial derivatives |
| `calculate_point_roots(char_poly, ...)` | Solve 1D polynomial roots at a parameter point |
| `WindingFun(char_poly, loop_fun, loop_range)` | Winding number integrand (polynomial) |
| `MatWindingFun(mat_fun, dmat_fun, ...)` | Winding number integrand (sparse matrix, requires BerryPy) |

### `brute_force_amoeba`

| Function | Description |
|----------|-------------|
| `get_hungarian_sorted_roots(char_poly, E_ref, mu1)` | Continuous $\beta_2$ root tracks via Hungarian matching |
| `get_a2_average_winding(char_poly, E_ref, mu1, mu2)` | $\partial R / \partial \mu_2$ — average winding in a2 direction |
| `get_a1_average_winding(char_poly, E_ref, mu1, mu2)` | $\partial R / \partial \mu_1$ — average winding in a1 direction |
| `bisect_a2_winding(char_poly, E_ref, mu1, mu2_low, mu2_high)` | Bisect $\mu_2$ to find where a2 winding crosses target |
| `bisect_amoeba_ronkin_min(char_poly, E_ref, ...)` | Find $(\mu_1, \mu_2)$ where both a1=a2=0 simultaneously |

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

The SGBZ formulation is based on the strip winding number approach for 2D non-Hermitian systems. The amoeba formulation uses the Ronkin function of Laurent polynomials. See [doc/](doc/) for detailed theoretical background.

## License

MIT — see [LICENSE](LICENSE).
