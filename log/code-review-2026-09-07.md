# Code review — 2026-09-07

## Resolution — 2026-09-07

Implemented findings **1 and 3–7**. Finding **2 is intentionally unchanged at the user's request**. Documentation follows the existing backend-selection code.

| Finding | Resolution |
| --- | --- |
| 1 | Deduplication now requires matching mu1, circular theta1, and chordal beta2 distance to the same line sample. Added assembly regressions for different tracks, different radii, genuine duplicates, and seam duplicates. Matching remains sample-based. |
| 2 | No change requested; the debug seam plotting behavior is retained. |
| 3 | Replaced ndarray `.ptp()` with `np.ptp`. |
| 4 | Checkout tests explicitly load `src` and reject an already-loaded external package. Both direct-script and module imports of the debug demo select checkout source. |
| 5 | Replaced SIGALRM with a child-process run and a 60-second timeout, preserving the solver assertions and timeout protection on Windows and Unix. |
| 6 | An OK mesh verdict now requires zero degenerate faces. Extended the existing zero-area test to assert rejection. |
| 7 | README and CharPoly docstring now document environment override followed by NumPy default, explicit poly_tools selection, and the `auto` alias. Updated affected assembly and test-source documentation as well. |

Validation after revisions:

- `python -m pytest -q`: **278 passed, 16 skipped in 56.12 seconds**, using the checkout with Python 3.11.5, NumPy 1.26.4, and SciPy 1.11.4.
- After the final debug-package import-order adjustment, `python -m pytest tests/test_debug_tool.py -q`: **18 passed, 1 skipped**.
- A fresh interpreter confirmed both the test bootstrap and debug demo module import resolve to this checkout's `src/pygbz2d/__init__.py`.
- NumPy **2.0.2** / SciPy **1.13.1** smoke check: default automatic mu scaling matches explicit scaling; a four-point cloud clusters into the expected two components. The full suite was not run under NumPy 2.
- `git diff --check` passed. Slow and optional skipped tests were not enabled.

The original review below is preserved as a historical record. Checkboxes reflect completed work. Unchecked items identify deliberately excluded work or validation that has not yet been performed. Failure counts and reproducer outputs below describe the pre-fix checkout.

---

Reviewed checkout: `fc65bc9`. Scope: debug CLI and plotting, mesh orientation checker, core types and NumPy backend, amoeba bisection and subset assembly, SGBZ solver, selected continuation paths, experimental band clustering, packaging, and tests. This is a targeted correctness review, not an exhaustive verification of the numerical algorithms. No implementation files were changed.

Seven actionable findings follow. P2 means a correctness or validation issue to address; P3 means a documentation correction. Within P2, address the subset-loss issue first.

## 1. [P2] Amoeba assembly drops distinct points sharing only theta1

Location: `src/pygbz2d/amoeba/amoeba.py:328–339`; caller at line 321.

`_point_near_any_line` compares only the point's first angle to sampled line angles. It never compares beta2. During continuum assembly, a discrete crossing on a different root track is therefore discarded whenever its theta1 matches any line sample within `SNAP_TOL`, even if its second angle is pi away. Mixed line/point results can silently lose components and report an incorrect index.

**Evidence:** a line with theta1 samples `[0, 0.5, 1]` and beta2 identically `+1`, and a point with beta1 `exp(0.5j)` and beta2 `-1`, returns `True` from the duplicate predicate. An assembly-level check with the crossing finder and splicer stubbed returned only `LineSubset`, although the supplied crossing is distinct from that line. This demonstrates the assembly defect; an end-to-end physical model producing that mixed configuration was not established in this review.

- [x] Require proximity in both Bloch coordinates, with appropriate periodic/chordal distance, before removing a crossing. Handle interpolation or endpoint matching explicitly rather than treating a theta1 sample match as point identity.
- [x] Add a regression retaining different-beta2 crossings at the same theta1, alongside a genuine duplicate-removal case.

## 2. [P2] Debug plot never splits seam crossings at its default threshold

**Left unfinished by request:** the user explicitly excluded issue #2 because the debug script is unofficial and its current plotting behavior is acceptable. Both tasks below intentionally remain unchecked.

Location: `debug_tool/gbz_debug.py:799–815`.

`_split_torus_segments` wraps coordinates into `[0, 2pi)` and then checks whether `circ_dist(...) > pi`. Shortest circular distance is at most pi, so neither condition can become true with the default `jump`. Matplotlib consequently draws a long straight segment across the plot when a curve passes through either seam. This can misrepresent the connectivity of the debug output.

**Evidence:** theta1 `[6.0, 6.2, 0.1, 0.3]` and constant theta2 produce one unchanged segment, instead of separate segments on the two sides of the seam.

- [ ] Detect discontinuities using absolute differences of the wrapped display coordinates; optionally insert seam endpoints to preserve visible curve extent.
- [ ] Cover theta1 seams, theta2 seams, and ordinary continuous curves.

## 3. [P2] Default experimental embedding is incompatible with NumPy 2

Location: `src/pygbz2d/experimental/band_clustering.py:213`; dependency declaration in `pyproject.toml` permits NumPy 2.

For any nonempty cloud with the default `w_mu=None`, `embed` calls `bp.mu1.ptp()` and `bp.mu2.ptp()`. These ndarray methods were removed in NumPy 2, so `cluster_bands` fails before constructing the radius graph in an otherwise permitted environment. Explicitly providing `w_mu` avoids this particular branch, but the default API remains broken.

**Evidence level:** source inspection plus the official [NumPy 2 migration guide](https://numpy.org/doc/2.0/numpy_2_0_migration_guide.html#ndarray-and-scalar-methods), which specifies `np.ptp(arr, ...)` as the replacement. The local environment is NumPy 1.26.4; this review did not execute a NumPy 2 environment.

- [x] Replace both method calls with `np.ptp(array)`.
- [ ] Validate the existing experimental tests against NumPy 1 and 2 with compatible SciPy versions.
  **Partially validated:** the full suite passed with NumPy 1.26.4, and a NumPy 2.0.2 / SciPy 1.13.1 embedding-and-clustering smoke test passed. The complete experimental suite has not been run under NumPy 2.

## 4. [P2] Source-checkout tests can silently exercise an unrelated installed package

Location: `conftest.py:4–10`.

The test bootstrap prefers any importable installed `pygbz2d`, adding this checkout's `src` only if that import fails. With a non-editable installation present, edits in the checkout are not necessarily tested. This also affects the debug CLI, which imports `pygbz2d` before its model builder and does not explicitly select the source package.

**Evidence:** the default interpreter resolved `pygbz2d.__file__` to `D:\anaconda3\Lib\site-packages\pygbz2d\__init__.py` in this session. The initial ordinary pytest run therefore cannot be treated as validation of this checkout. The review reran the suite with `src` first on `sys.path` and verified the imported path.

- [x] Make checkout tests explicitly select the local source, or fail early when the imported package is outside the intended checkout. Keep wheel-installation testing as an explicit separate mode.
- [x] Document or enforce the intended package source for the debug CLI.

## 5. [P2] Four regression tests fail before testing the solver on Windows

Location: `tests/test_regressions.py:54–60`.

`TestH1MrRestartPingpong.test_run_terminates_without_zero_mr` unconditionally accesses `signal.SIGALRM` and `signal.alarm`. On this Windows environment the former is absent, so all four parameterized cases raise `AttributeError` before `zm.run()` executes. The suite is red and the intended termination checks are not performed.

**Evidence:** reproduced against this checkout for `(h0, min_dtheta)` values `(0.5, 1e-6)`, `(0.1, 1e-6)`, `(0.2, 1e-4)`, and `(1.0, 1e-6)`.

- [x] Use a portable process-based timeout or a supported timeout test dependency. Preserve a real timeout so a recurrence of the original infinite loop cannot hang the suite.
- [ ] Run the regression on Windows as well as a Unix platform.
  **Partially validated:** all four subprocess regression cases passed on Windows. A Unix run has not been performed in this session.

## 6. [P2] Mesh checker reports zero-area triangles as OK

Location: `debug_tool/mesh_orientation.py:69–82`.

`orientation_uniform` accepts a census with no positive or no negative triangles, and `is_ok` never checks `signed_zero`. A wholly collinear triangle therefore receives `verdict: OK`. The CLI uses that verdict to choose its exit status, so geometrically degenerate meshes can pass the diagnostic.

**Evidence:** vertices `[[0, 0], [0.1, 0], [0.2, 0]]` and triangle `[[0, 1, 2]]` yield `signed_zero=1`, `orientation_uniform=True`, and `is_ok=True`. The existing `test_zero_area_counted` checks the census but does not check rejection.

- [x] Require nondegenerate faces for an OK verdict, or expose an explicit indeterminate/degenerate verdict that causes CLI failure. Retain the separate winding census if its current semantics are useful.
- [x] Extend the zero-area regression to assert the verdict.

## 7. [P3] README promises backend auto-selection that the factory does not perform

Location: `README.md:24–38`; `src/pygbz2d/core.py:117–120`; implementation at `src/pygbz2d/backend.py:267–284`.

The installation documentation and `CharPoly` docstring claim that an importable `poly_tools` is selected automatically, otherwise NumPy is used with a warning. The current factory instead defaults to NumPy unless the argument or `POLY_BACKEND` explicitly selects `poly_tools`; even `"auto"` selects NumPy. Users following the acceleration instructions can therefore install the extension and still run NumPy without realizing it.

- [x] Align the README and constructor docstring with the implemented selection order, including an explicit acceleration example. The backend's comment says the explicit-selection policy is deliberate because of multiple-root issues.

## Validation and limitations

Environment: Windows, Python 3.11.5, NumPy 1.26.4, SciPy 1.11.4.

The checkout-targeted command, run from the repository root, was:

```powershell
python -c "import sys; sys.path.insert(0, 'src'); import pygbz2d, pytest; print('Review target:', pygbz2d.__file__, flush=True); raise SystemExit(pytest.main(['-q']))"
```

Verified import: `D:\research-works\brute-force-SGBZ-amoeba\src\pygbz2d\__init__.py`.

Result: **270 passed, 16 skipped, 4 failed in 57.26 seconds**. All failures were the Windows watchdog issue in finding 5. Slow tests were not enabled; optional-backend and other skipped cases are not validated. The passing suite does not establish correctness for all degeneracies or physical models.

Minimal reproductions for findings 1, 2, and 6, with the checkout selected explicitly:

```python
import sys
sys.path.insert(0, 'src')
import numpy as np
from pygbz2d.core import LineSubset, PointSubset
from pygbz2d.amoeba.amoeba import _point_near_any_line
from debug_tool.gbz_debug import _split_torus_segments
from debug_tool.mesh_orientation import check_mesh_orientation

line = LineSubset(0j, 0.0, np.array([0., .5, 1.]), np.ones(3, complex))
point = PointSubset(0j, np.exp(.5j), -1 + 0j)
print(_point_near_any_line(point, [line], 1e-6))  # True; should be False

segs = _split_torus_segments(np.array([6., 6.2, .1, .3]), np.ones(4))
print([x.tolist() for x, _ in segs])  # One segment crossing the plot

rep = check_mesh_orientation(
    np.array([[0., 0.], [.1, 0.], [.2, 0.]]), np.array([[0, 1, 2]]))
print(rep.describe())  # zero=1, verdict: OK
```

Existing TODOs on exact degeneracy, warm starts, and CI were read and remain separate work. The untracked paper artifacts were outside the code-review scope.
