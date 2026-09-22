# AGENTS.md — brute-force-non-hermitian

Non-Hermitian skin effect computation for 2D tight-binding models. Packaged as **pyGBZ2d** (import name `pygbz2d`, src-layout, `pip install -e .`). Two complementary GBZ modules (SGBZ and amoeba) implementing brute-force polynomial root-solving approaches, plus a pseudo-arclength continuation module for adaptive root tracking along θ₁. All GBZ entry points (`collect_GBZ_subsets`) natively return the unified `GBZResult` type defined in `pygbz2d/core.py`.

## Project Structure

```
brute-force-non-hermitian/
├── pyproject.toml             # pyGBZ2d packaging (deps: numpy, scipy only)
├── src/pygbz2d/               # The installable package
│   ├── __init__.py            # Top-level API: CharPoly, GBZResult, Point/LineSubset, TWO_PI
│   ├── core.py                # Unified data types shared by both modules:
│   │                          #   CharPoly (the single polynomial entry point),
│   │                          #   PointSubset (0D), LineSubset (1D, eager beta2_arr),
│   │                          #   GBZResult, ConnectedSubset = Union[Point, Line],
│   │                          #   TWO_PI (shared 2π constant)
│   │                          # + shared utils: sort_by_root_abs, chordal_cost_matrix,
│   │                          #   hungarian_match_indices, find_cyclic_true_intervals,
│   │                          #   get_minor_degrees, generate_probe_steps, to_sphere_r3,
│   │                          #   circ_dist, check_points_clustered_on_torus,
│   │                          #   probe_zero_plateau, JoinableLinePiece,
│   │                          #   is_mr_cluster_endpoint (cross-module MR join test)
│   ├── backend.py             # Pluggable Laurent backends: LaurentProtocol,
│   │                          #   PolyToolsLaurent (lazy C++ import), NumpyLaurent
│   │                          #   (pure-numpy fallback), make_laurent factory
│   │                          #   (arg > POLY_BACKEND env > numpy default;
│   │                          #    poly_tools only when explicitly requested)
│   ├── (constants)            # Policy, not a directory: constants live in
│   │                          #   their home modules (single-consumer locality);
│   │                          #   the 7 cross-package ones sit in core.py with
│   │                          #   live_defaults.  Full map: doc/constants.md.
│   │                          #   Customization: per-call kwarg OR direct module
│   │                          #   constant assignment (live, process-wide).
│   ├── sgbz/                  # SGBZ / average major-axis winding formulation
│   │   ├── __init__.py        # Whitelist exports (re-exports core classes)
│   │   ├── mu2mid.py          # ItemView + Mu2Mid piecewise-Hermite path + Mu2MidZM
│   │   │                      #   (ZeroManager subclass; analyze() = cluster +
│   │   │                      #   multi-crossing mesh refinement + pairwise + insert
│   │   │                      #   + finalize + μ₂_mid build; has_continuum)
│   │   ├── pairwise.py        # THE crossing channel: refine_mesh_for_multiple_crossings
│   │   │                      #   (pre-scan sub-mesh for intervals with ≥2 predicted roots),
│   │   │                      #   collect_pair_events (touch/cross, brentq refinement,
│   │   │                      #   is_mr marking), EventGroup merge, insert_event_groups
│   │   │                      #   (regular separator rows between adjacent events),
│   │   │                      #   finalize_event_groups (side-change charges)
│   │   ├── continuum_lines.py # Continuum LineSubset extraction (detect_continuum_simple,
│   │   │                      #   extract_continuum_linesubsets, MR/seam joining)
│   │   ├── winding.py         # W(E_ref, mu1) loop-winding + 0D PointSubset materialization
│   │   │                      #   from EventGroups (detect_crossings_simple — pure
│   │   │                      #   materialization, no MR channel / echo drop)
│   │   ├── plateau.py         # Zero-plateau detection (clustering + probe)
│   │   └── sgbz_solver.py     # solve_SGBZ_for_E, collect_GBZ_subsets (returns GBZResult)
│   ├── amoeba/                # Amoeba / Ronkin function formulation
│   │   ├── __init__.py        # Whitelist exports (re-exports core classes)
│   │   ├── ronkin_winding.py  # Brent refinement with retained mesh rows;
│   │   │                      #   average winding from root counts
│   │   ├── bisect.py          # μ₁/μ₂ bisection, fast μ₂ gap test, continuum handling
│   │   ├── zm_extract.py      # AmoebaZeroManager, detect_continuum, find_crossings,
│   │   │                      #   calculate_a2_average_winding
│   │   └── amoeba.py          # collect_GBZ_subsets, subset assembly, plateau check
│   ├── continuation/          # Pseudo-arclength continuation for β₂-root tracking along θ₁
│   │   ├── __init__.py        # Public API re-exports
│   │   ├── interpolation.py   # hermite_interp_poly (cubic Hermite kernel shared by
│   │   │                      #   predict_roots_hermite and Mu2Mid pieces)
│   │   ├── arclength.py       # compute_tangent, predict_roots, estimate_error,
│   │   │                      #   arclength_step, StepControl (fields resolve from
│   │   │                      #   this module's constants at construction)
│   │   ├── multiple_roots.py  # MR detection: point/interval triggers, detect_cluster,
│   │   │                      #   solve_multiple_roots_in_interval, MultipleRootInfo
│   │   └── zero_manager.py    # ZeroManager orchestrator, integrate_segment, SegmentData
│   └── experimental/          # EXPERIMENTAL post-processing, NO stability guarantee
│       ├── __init__.py        # explicit opt-in imports (not in top-level API)
│       └── band_clustering.py # Radius-graph band clustering of GBZ sweep point clouds:
│                              #   flatten_results (no line decimation) → (cos,sin) torus
│                              #   embed → cKDTree + connected components; cluster_bands
│                              #   entry point; eps-window diagnostics (stability scan,
│                              #   widest plateau, inter-cluster margins); scan capped
│                              #   (EPS_SCAN_MAX) + early-stop at full merge
├── tests/                      # pytest: test_backend.py (dual-backend parity),
│                               #   test_constants.py (live-assignment contract + AST lint),
│                               #   test_gbz_types→core, test_sgbz, test_amoeba,
│                               #   test_continuation, test_zero_manager,
│                               #   test_interpolation, test_counterexamples
│                               #   (slow: --run-slow), test_regressions, test_debug_tool
├── conftest.py                 # checkout-src import + shared build_HN2D_polynomial + slow marker
├── application/                # HN benchmark script, Typst guide and compiled PDF;
│                               #   fixed examples, random comparisons, coarse/fine sweeps
├── playground/                 # Unofficial runnable demo scripts (demo_unified.py shows the GBZResult API;
│                               #   Haldane-model-gainloss.py: sweeps, plots, recompute_failed_SGBZ)
├── debug_tool/                 # Fixed-(E_ref, mu1) debugging (see debug_tool/README.md):
│                               #   gbz_debug.py — collect_debug_subsets (GBZDebugReport: subsets +
│                               #   SGBZ charges of BOTH methods at a frozen mu1), compute_loop_windings
│                               #   (loop winding numbers + charge-consistency gap checks),
│                               #   plot_winding_debug (theta1-theta2 torus figure)
│                               #   demo_debug_tool.py — CLI (Haldane E=1.212 debug point, HN2D)
├── diagnostics/                # Debug scripts and analysis reports
├── data/                       # Pickled computation results (demos / replication).
│                               #   NOTE: pre-packaging pickles reference the old
│                               #   gbz_types module path and no longer load.
├── log/                        # Change logs
│   ├── ......
├── paper/                      # Folder for official PRB paper for this project
│   ├── ......
├── doc/                        # Documentation
│   ├── SGBZ.md                 # SGBZ theory, architecture, API
│   ├── amoeba.md               # Amoeba theory, algorithm, API
│   ├── continuation.md         # Continuation module (arclength/MR/ZeroManager/interpolation)
│   └── constants.md           # Per-module numerical-constant reference (user-facing)
├── TODO/                       # Optimization checklists + engineering plans
│                               #   (exact-degeneracy-boundary-behavior.md — open;
│                               #   warm-start ZM — open; resolved items are
│                               #   archived as meeting records in log/)
└── README.md
```

## Development Principles

- **Unexpected output is an opportunity to improve the algorithm.** Investigate the root cause of unexpected behavior, such as non-monotonic winding or a false plateau classification, and repair the underlying algorithm. Do not substitute threshold changes, extra pre-checks, or restricted probe ranges for that investigation. In the plateau regression, commit `3cdb22d` produced the correct result before plateau checking was introduced, locating the regression in the later logic. Comparing commits is a useful diagnostic method.

- Comments should explain why the implementation is needed, rather than restating what it does.

- Ask why before deciding how. Check the facts before drawing conclusions. Do not assert a bug without solid evidence.

- Stop immediately and report any result that differs from expectations.

- Obtain my approval for the modification plan before implementing it.

- You can answer either in Chinese or English. I can understand both.
