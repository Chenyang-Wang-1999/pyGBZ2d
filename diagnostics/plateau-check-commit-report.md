# Commit Report: Zero-Plateau Post-Check for Amoeba and SGBZ

## Suggested Commit Message

```
fix: guard amoeba and SGBZ classification against zero-plateau edge hits
```

## Background

The reference energy `E = 1.648 + 0.0294j` in
`demos/imaginary-degeneracy-splitting.py` exposed a classification error:
the outer bisection can terminate at the edge of a zero-winding plateau.  At
that edge, raw zero/PMGBZ data may be non-empty, so the old final classifier
treated the point as inside the amoeba/SGBZ spectrum even though nearby points
show a finite zero plateau.

The diagnostic sweep found:

- No decreasing interval in the main-direction winding curve.
- A finite zero plateau for the amoeba criterion:
  `mu1 in [0.113282599424, 0.11647195474]`.
- The original `collect_GBZ_subsets` returned `mu1 = 0.11328125`, immediately next to
  the plateau boundary, with non-empty raw zeros but net-zero crossing content.

## Code Changes

### `brute_force_amoeba/amoeba.py`

- Added `_make_ronkin_result(...)` so `bisect_amoeba_ronkin_min(...)` returns
  local bisection metadata:
  - `_mu1_bracket`
  - `_a1_bracket`
  - `_exit_reason`
- Added `_probe_zero_plateau_near_mu1(...)`.
  It uses multi-scale `mu1` perturbations combining:
  - the final bisection bracket width,
  - `xtol`,
  - a fixed probe radius, defaulting to `continuum_perturb`.
- Updated `collect_GBZ_subsets(...)`:
  - `zeros == []` and non-continuum still directly means outside spectrum.
  - `zeros != []` no longer directly means inside spectrum.
  - A nearby point with `a1 == 0`, `zeros == []`, and non-continuum now changes
    the final classification to `is_amoeba=False`.
- Added diagnostic return fields:
  - `_plateau_check`
  - `_plateau_check_found`
  - `_plateau_probe_points`
  - `_classification_reason`

### `brute_force_SGBZ/SGBZ.py`

- Added `SGBZSolver.solve_for_E_info(...)`, preserving the original
  `solve_for_E(...)` return shape for compatibility.
- Added SGBZ zero-plateau probing based on strip winding:
  - A probe point is considered on a plateau when `abs(W) <= zero_tol`,
    `PMGBZ_points == []`, and the point is not a continuum point.
- Updated `check_SGBZ(...)`:
  - Empty `PMGBZ_points` at the returned root is classified as
    `is_PMGBZ=False`.
  - Non-empty `PMGBZ_points` are post-checked with multi-scale perturbations.
  - If a nearby zero plateau is found, the result is changed to
    `is_PMGBZ=False`.
- Added diagnostic return fields analogous to the amoeba path.

### Diagnostics

- Added `diagnostics/diagnose_imaginary_degeneracy_splitting.py`.
- Generated:
  - `diagnostics/imaginary-degeneracy-splitting/report.md`
  - `diagnostics/imaginary-degeneracy-splitting/a1_vs_mu1.csv`

The diagnostic script explicitly prepends the repository root to `sys.path`, so
it tests the working-tree package rather than an installed older package.

## Verification

Syntax checks:

```
python -m py_compile brute_force_amoeba/amoeba.py diagnostics/diagnose_imaginary_degeneracy_splitting.py
python -m py_compile brute_force_SGBZ/SGBZ.py brute_force_SGBZ/strip_winding_number.py brute_force_SGBZ/root_solver.py brute_force_SGBZ/winding.py brute_force_SGBZ/__init__.py
```

Amoeba diagnostic:

```
python diagnostics/diagnose_imaginary_degeneracy_splitting.py
```

Observed result:

```
success=True
is_amoeba=False
classification=nearby_zero_plateau
plateau_check=found
mu1=0.11328125
mu2=0.0
zero_count=4
net_zero_count=0
```

SGBZ smoke check on the same reference energy:

```
is_PMGBZ=False
_classification_reason=empty_PMGBZ_zero_plateau
mu1=0.11427582360136657
winding=-0.0
PMGBZ_count=0
```

## Staging Notes

Include the intended code/report files only.  Do not stage Python cache files.

Suggested relevant files:

```
brute_force_amoeba/amoeba.py
brute_force_SGBZ/SGBZ.py
diagnostics/diagnose_imaginary_degeneracy_splitting.py
diagnostics/imaginary-degeneracy-splitting/report.md
diagnostics/imaginary-degeneracy-splitting/a1_vs_mu1.csv
diagnostics/plateau-check-commit-report.md
```

Current working tree also contains broader unrelated changes, including the
`brute_force_solver` to `brute_force_SGBZ` package transition and documentation
updates.  Review `git status --short` before staging to avoid mixing unrelated
work into this commit.
