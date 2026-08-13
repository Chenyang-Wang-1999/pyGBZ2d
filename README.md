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

## Known Issues

### Hungarian matching swap on coarse θ₁ grids in `brute_force_amoeba`

**Symptom**: `get_hungarian_sorted_roots` may produce spurious ln|β₂| sign changes (false crossings of |β₂|=1) on the default coarse grid (N_points=301). This causes `collect_GBZ_subsets` to misclassify a small number of energy points as inside the GBZ spectrum when they are not. The issue is most visible when comparing results from different polynomial representations of the same physical system (e.g., original unit cell vs. supercell).

**Root cause**: When two β₂ roots approach within ~1% of each other in the complex plane between adjacent θ₁ slices, the pure chordal-distance Hungarian matching can swap their identities. The "wrong" matching (crossing in |β₂|) has lower total chordal cost than the "correct" matching (preserving |β₂| ordering), so `linear_sum_assignment` selects it. Finer grids (N≥1001) resolve this by reducing the angular step Δθ₁, making root positions diverge less between slices.

**Affected models**: Supercell / multiband polynomials are more susceptible because the larger total degree compresses root trajectories into the same angular range, increasing the likelihood of near-degeneracies.

**Planned fixes**:

1. Augment the Hungarian cost matrix with a first-order Taylor prediction term using dβ₂/dθ₁ from implicit differentiation of `f(E, β₁, β₂)=0`. The prediction penalizes matches that violate analytic continuity, steering the matcher toward the physically correct assignment. See `brute_force_amoeba/tracks.py:get_hungarian_sorted_roots`.

2. Pseudo arc-length continuation: Instead of matching roots independently at each θ₁ slice, follow each root track along θ₁ by solving an augmented system `[f(E, β₁, β₂), |Δβ₂|² + |Δθ₁|² - ds²]` that parametrizes the root curve by arc length. This naturally handles near-degeneracies because the continuation step is controlled by the local curvature of the root trajectory rather than the θ₁ grid spacing.

**Workaround**: Increase `N_points` from 301 to 1001 in `collect_GBZ_subsets` options, or manually compare results from multiple polynomial representations.

## References

The SGBZ formulation is based on the average major-axis winding number approach for 2D non-Hermitian systems. The amoeba formulation uses the Ronkin function of Laurent polynomials. See [doc/](doc/) for detailed theoretical background.

## License

MIT — see [LICENSE](LICENSE).
