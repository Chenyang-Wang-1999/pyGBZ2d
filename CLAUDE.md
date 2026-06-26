# CLAUDE.md — brute-force-non-hermitian

Non-Hermitian skin effect computation for 2D tight-binding models. Two complementary modules implementing brute-force polynomial root-solving approaches.

## Project Structure

```
brute-force-non-hermitian/
├── brute_force_SGBZ/        # SGBZ / strip winding number formulation
│   ├── __init__.py             # Whitelist exports of all public symbols
│   ├── root_solver.py          # calculate_point_roots, complex_root, poly_to_np_coefficients
│   ├── winding.py              # PolyDiffContext, WindingFun, MatWindingFun, get_winding_number
│   ├── SGBZ.py                 # SGBZSolver, SGBZChecker, check_SGBZ
│   └── strip_winding_number.py # get_roots_and_PMGBZ, get_strip_winding, get_loop_winding
├── brute_force_amoeba/        # Amoeba / Ronkin function formulation
│   ├── __init__.py             # Whitelist exports (6 public functions)
│   └── amoeba.py               # All amoeba logic (~940 lines)
├── doc/                        # Documentation
│   ├── SGBZ.md                 # SGBZ theory, architecture, API
│   └── amoeba.md               # Amoeba theory, algorithm, API
├── tests/                      # pytest tests
├── demos/                      # Runnable demo scripts
├── pyproject.toml
└── README.md
```

## Module: brute_force_SGBZ (SGBZ formulation)

Computes the generalized Brillouin zone (GBZ) using strip winding numbers.

**Core data flow:**
```
Model → CLaurent polynomial → PolyDiffContext (cached derivatives)
  → get_roots_and_PMGBZ() → PMGBZ_points + sorted roots
  → get_strip_winding() → winding number W(E, μ₁)
  → SGBZSolver.solve_for_E() → critical μ₁ where W=0
```

**Key classes:**
- `PolyDiffContext(char_poly)` — bundles polynomial with precomputed partial derivatives. `eval_val(var)`, `eval_partials(var)`, `eval_dmu2(var)`.
- `SGBZSolver(char_poly)` — main entry point. `solve_for_E(E_ref, mu1_guess, zero_tol, N_points)` returns `(mu1, PMGBZ_points)`.
- `WindingFun(char_poly, loop_fun, loop_range)` — callable integrand for `scipy.integrate.quad`.
- `MatWindingFun` — sparse matrix variant (requires BerryPy, optional).

**Key functions:**
- `calculate_point_roots(char_poly, param_ind, param_val, var_ind, M_max, N_max)` — solve 1D polynomial at fixed params using `np.roots`. Returns list of complex roots (may include `0`/`inf` sentinels).
- `get_roots_and_PMGBZ(poly_diff, E_ref, mu1, N_points, ...)` — trace θ₁ around full circle, detect PMGBZ points where |β_M| = |β_{M+1}|. Returns `(PMGBZ_points, theta1_arr, sols_arr, info)`.
- `get_strip_winding(poly_diff, E_ref, mu1, ...)` — compute strip winding number. Continuum case returns `(W_left, W_right)` tuple.
- `get_loop_winding(poly_diff, E_ref, mu1, mu2_fun, theta2, ...)` — single-loop winding at fixed θ₂.

**PMGBZ point structure:** Each point is a dict with keys `theta1`, `beta2_sols` (list of [pos_tuple, neg_tuple] classified by dμ₂/dθ₁ sign), `is_continuum`.

**Key tolerances** (in strip_winding_number.py):
- `zero_tol = 1e-10` — PMGBZ degeneracy detection
- `GBZ_check_tol = 1e-6` — boundary cluster tolerance
- `match_confidence_tol = 1e-3` — Hungarian matching confidence threshold

## Module: brute_force_amoeba (Amoeba formulation)

Computes the amoeba spectrum via the Ronkin function and average winding numbers.

**Algorithm pipeline:**
```
1. get_hungarian_sorted_roots()  → continuous β₂ root tracks (Hungarian matching)
2. _compute_winding_from_tracks() → detect ln|β₂| = μ₂ crossings
3. _get_average_winding_from_zeros() → partition θ circle, compute average winding
4. bisect μ₂ (inner) → find where a₂ = 0
5. bisect μ₁ (outer) → find where a₁ = 0 (with inner loop guaranteeing a₂ = 0)
```

**Public API (6 functions):**
- `get_hungarian_sorted_roots(char_poly, E_ref, mu1, N_points)` → `(theta1_arr, tracked, M, N)`
- `get_a2_average_winding(char_poly, E_ref, mu1, mu2, ...)` → ∂R/∂μ₂ (float)
- `get_a1_average_winding(char_poly, E_ref, mu1, mu2, ...)` → ∂R/∂μ₁ (float), reuses a2 zeros
- `bisect_a2_winding(char_poly, E_ref, mu1, mu2_low, mu2_high, ...)` → dict with mu2, zeros, is_continuum
- `bisect_amoeba_ronkin_min(char_poly, E_ref, mu1_low, mu1_high, mu2_low, mu2_high, ...)` → dict with mu1, mu2, zeros, is_continuum
- `check_amoeba(coeffs, degs, E_ref, perc, ...)` → high-level wrapper for multiprocessing sweeps

**Key design decisions:**
- Hungarian matching with chordal (Riemann sphere) cost metric avoids 0/∞ false swaps.
- Root tracks depend only on (E, μ₁) and are cached — independent of μ₂.
- Lazy refinement: crossings are NOT fsolve-refined during bisection iterations; only the final μ₂ gets Newton-corrected using analytical dW/dμ₂.
- a₁ winding reuses a₂ zero points (no separate Hungarian matching needed).
- Continuum detection uses |ln|β₂| - μ₂| < continuum_tol (simpler than SGBZ's |β_M| ≈ |β_{M+1}|).

**Cross-module dependency:** `amoeba.py` imports `calculate_point_roots` and `PolyDiffContext` from `brute_force_SGBZ`.

## Polynomial Format

All characteristic polynomials use `poly_tools.CLaurent(3)` (3-variable Laurent polynomial). Encoding: each term is a triplet `[E_exponent, beta1_exponent, beta2_exponent]`.

Example — square lattice NN: f = E - β₁ - β₁⁻¹ - β₂ - β₂⁻¹
```python
coeffs = [1, -1, -1, -1, -1]
degs = [1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1]
```

## External Dependencies

- **poly_tools** (required): C-extension Laurent polynomial library. Clone from https://atomgit.com/wangchenyang99/PolyTools, then build via `make` in the `pybind/` directory and copy `python/poly_tools/` to your workdir. Python 3.9+ only (compiled .so).
- **numpy, scipy** (required): standard numerical computing.
- **BerryPy** (optional): only for `MatWindingFun` sparse matrix winding. Not needed for core functionality.

## Build / Test / Run

```bash
pip install -e .                    # install in dev mode
pip install -e ".[dev]"             # with pytest
pytest tests/ -v                    # run tests (skips if poly_tools missing)
python demos/demo_amoeba.py         # full amoeba pipeline demo
python demos/demo_solver.py         # SGBZ solver demo
```

## Code Conventions

- Private/internal functions prefixed with `_`.
- Public API exported via explicit `__all__` in `__init__.py` — no wildcard exports.
- Type annotations used on public functions and class attributes.
- No backwards-compatibility shims — removed interfaces are deleted, not deprecated.
- Comments describe WHY (mathematical motivation), not WHAT (the code is self-documenting).
- Documentation in `doc/` and README are in English.
