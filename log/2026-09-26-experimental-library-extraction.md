# 2026-09-26: remove playground dependencies from the QWZ application

The user requested an audit of the QWZ example's imported functions,
classification into existing library code / reusable algorithms / outer
orchestration, and removal of every playground dependency.

The audit is recorded in `doc/experimental_mesh.md`. Before this change,
`experimental` contained only band clustering and flattening. The reusable
geometry, refinement and Chern functions existed only in four playground
scripts; there was no second package implementation to reuse.

## Boundaries after extraction

- `experimental.torus_mesh`: periodic Delaunay, vertex deduplication, mesh
  topology, torus geometry and orientation. `build_band_mesh` takes an already
  selected `BandPoints` cloud; `build_cluster_mesh` adapts `BandClustering`.
- `experimental.mesh_refinement`: periodic complex-energy interpolation,
  dyadic target construction, matching, index-boundary bisection, batch
  solving, and `refine_mesh` on in-memory arrays. The multiprocessing worker
  is a module-level package function. No Matplotlib or demo imports.
- `experimental.chern`: SVD/eigenvector diagnostics, RR and reciprocal LR
  fluxes, and `integrate_chern` returning `(ChernResult, flux_data)`.
  Numerical integration does not read or write files or discover model files.
- Application: QWZ model, scan windows and workers, band selection, partner
  symmetry, filenames, persistence, plot masks, reports, source provenance,
  and CLI. The NPZ wrapper is local to the example, with no sibling imports.
- Playground: previous demo/CLI entry points and compatibility imports remain;
  their numerical functions now call the package. No duplicate implementation
  is retained in those scripts.

Removed `sys.path` injection and all fixed checkout-source paths from the
application. Default paths are resolved beside the script, so copying it to
another directory is supported with an installed package. Source hashes are
optional metadata collected from the loaded modules. NumPy/SciPy remain the
only package dependencies; example-only plotting/model dependencies stay
outside the package.

Default numerical parameters are unchanged. Optional public numeric defaults
resolve from their home modules at call time. Mesh limitations and the
reciprocal-link Chern convention are unchanged. Package documentation and the
Typst/PDF example now name the installed APIs.

## Verification

- The relevant regression/constant suites report **51 passed**:
  `test_experimental.py`, `test_experimental_mesh.py`,
  `test_chern_biorthogonal.py`, `test_chern_example_independence.py`,
  `test_constants.py`.
- All 16 base/refined, lower/upper, RR/LR Chern results and checked numerical
  diagnostics exactly match the pre-migration JSON. Both cases' two base
  meshes, rebuilt from the saved fine sweep, reproduce all vertex,
  triangle, energy and logarithmic-radius arrays exactly.
- Previous source and result snapshots are retained under
  `application/data/QWZ-benchmark/library-migration-baseline`.
- Built a wheel, installed it into a temporary directory outside the checkout,
  copied only the application script, and launched an isolated Python process
  with checkout paths removed and demo imports explicitly forbidden. Synthetic
  QWZ band meshes passed RR/LR integration, BZ checks, diagnostics, and report
  generation. Two real one-band GBZ solves under spawned worker processes
  reproduced the serial arrays exactly. The loaded package path was verified
  to be the temporary wheel installation, not this checkout or an editable
  package. Evidence is retained in `library-migration-portable.log` and its
  JSON manifest under the default data directory.
- The local interpreter initially had an older non-editable pyGBZ2d install;
  removing path injection exposed that stale version. Installed the current
  checkout with `python -m pip install -e . --no-deps --no-build-isolation`.
  No third-party dependency versions were changed.

No full energy sweep or production mesh refinement was repeated for this
refactoring. Saved mesh files and previous scientific conclusions are retained.
