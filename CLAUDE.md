# CLAUDE.md — brute-force-non-hermitian

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
│   │                          #   TWO_PI (the single project-wide 2π constant)
│   │                          # + shared utils: sort_by_root_abs, chordal_cost_matrix,
│   │                          #   hungarian_match_indices, find_cyclic_true_intervals,
│   │                          #   get_minor_degrees, generate_probe_steps, to_sphere_r3,
│   │                          #   circ_dist, check_points_clustered_on_torus,
│   │                          #   probe_zero_plateau, JoinableLinePiece,
│   │                          #   is_mr_cluster_endpoint (cross-module MR join test)
│   ├── backend.py             # Pluggable Laurent backends: LaurentProtocol,
│   │                          #   PolyToolsLaurent (lazy C++ import), NumpyLaurent
│   │                          #   (pure-numpy fallback), make_laurent factory
│   │                          #   (arg > POLY_BACKEND env > poly_tools > numpy+warning)
│   ├── (constants)            # NO central config: numerical constants live in
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
│   │   ├── ronkin_winding.py  # Ronkin winding (avg windings from zeros, non-zero area)
│   │   ├── bisect.py          # Bisection + Newton refinement (bisect_amoeba_ronkin_min)
│   │   ├── zm_extract.py      # extract_amoeba_subsets from track crossings of ln|β₂|=μ₂
│   │   └── amoeba.py          # collect_GBZ_subsets (returns GBZResult), plateau check, orchestration
│   └── continuation/          # Pseudo-arclength continuation for β₂-root tracking along θ₁
│       ├── __init__.py        # Public API re-exports
│       ├── interpolation.py   # hermite_interp_poly (cubic Hermite kernel shared by
│       │                      #   predict_roots_hermite and Mu2Mid pieces)
│       ├── arclength.py       # compute_tangent, predict_roots, estimate_error,
│       │                      #   arclength_step, StepControl (fields resolve from
│       │                      #   this module's constants at construction)
│       ├── multiple_roots.py  # MR detection: point/interval triggers, detect_cluster,
│       │                      #   solve_multiple_roots_in_interval, MultipleRootInfo
│       └── zero_manager.py    # ZeroManager orchestrator, integrate_segment, SegmentData
├── tests/                      # pytest: test_backend.py (dual-backend parity),
│                               #   test_constants.py (live-assignment contract + AST lint),
│                               #   test_gbz_types→core, test_sgbz, test_amoeba,
│                               #   test_continuation, test_zero_manager,
│                               #   test_interpolation, test_counterexamples
│                               #   (slow: --run-slow), test_regressions, test_debug_tool
├── conftest.py                 # installed-package-or-src import + shared build_HN2D_polynomial + slow marker
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
│   ├── constants.md           # Per-module numerical-constant reference (user-facing)
│   ├── 拓扑匹配算法说明.md       # Topological matching algorithm
│   ├── sn-main.tex             # Simplified paper for SGBZ, main text
│   └── sn-supp.tex             # Simplified paper for SGBZ, supplementary information. Amoeba GBZ is discussed in section{Comparison with reported frameworks}
├── TODO/                       # Optimization checklists + engineering plans
│                               #   (exact-degeneracy-boundary-behavior.md — open;
│                               #   warm-start ZM — open; resolved items are
│                               #   archived as meeting records in log/)
└── README.md
```

## Development Principles

- **反常输出是改进算法的机会，不是需要绕过的 bug。** 遇到 unexpected behavior（如 winding 非单调、plateau check 误触发）时，优先追查根因并修复底层算法，而不是加 workaround（如调阈值、加预检查、限制探针范围等）。3cdb22d 版本没有 plateau check 时反而结果正确，说明问题出在后来引入的逻辑。对比不同 commit 是定位问题的有效手段。

- 注释应当讲述“为什么”，而不是“做什么”。

- Ask WHY before HOW. Check the FACT before reach the CONCLUSION. Never assert a bug before you get the solid evidence.

- 遇到与预期不符的结果时，要立刻汇报。
