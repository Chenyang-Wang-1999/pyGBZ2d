# CLAUDE.md — brute-force-non-hermitian

Non-Hermitian skin effect computation for 2D tight-binding models. Two complementary GBZ modules (SGBZ and amoeba) implementing brute-force polynomial root-solving approaches, plus a pseudo-arclength continuation module for adaptive root tracking along θ₁. All GBZ entry points (`collect_GBZ_subsets`) natively return the unified `GBZResult` type defined in `gbz_types.py`.

## Project Structure

```
brute-force-non-hermitian/
├── gbz_types.py             # Unified data types shared by both modules:
│                            #   PointSubset (0D), LineSubset (1D, lazy beta2_arr),
│                            #   GBZResult, ConnectedSubset = Union[Point, Line]
│                            # + shared utils: sort_by_root_abs, chordal_cost_matrix,
│                            #   hungarian_match_indices, find_cyclic_true_intervals,
│                            #   get_minor_degrees, generate_probe_steps, to_sphere_r3
├── brute_force_SGBZ/        # SGBZ / average major-axis winding formulation
│   ├── __init__.py             # Whitelist exports (re-exports gbz_types classes)
│   ├── continuum.py            # Continuum detection (detect_continuum_simple/full)
│   ├── crossings.py            # Crossing detection + charge classification
│   ├── winding.py              # Average winding computation (compute_average_winding)
│   ├── plateau.py              # Zero-plateau detection (clustering + probe)
│   └── sgbz_solver.py          # solve_SGBZ_for_E, collect_GBZ_subsets (returns GBZResult)
├── brute_force_amoeba/        # Amoeba / Ronkin function formulation
│   ├── __init__.py             # Whitelist exports (re-exports gbz_types classes)
│   ├── tracks.py               # Root tracks via Hungarian matching (get_hungarian_sorted_roots)
│   ├── ronkin_winding.py       # Ronkin winding (get_a1/a2_average_winding, crossing detection)
│   ├── bisect.py               # Bisection + Newton refinement (bisect_amoeba_ronkin_min)
│   └── amoeba.py               # collect_GBZ_subsets (returns GBZResult), plateau check, orchestration
├── continuation/               # Pseudo-arclength continuation for β₂-root tracking along θ₁
│   ├── __init__.py             # Public API re-exports
│   ├── arclength.py            # compute_tangent, predict_roots, estimate_error, arclength_step
│   ├── multiple_roots.py       # MR detection: point/interval triggers, detect_cluster,
│   │                           #   solve_multiple_roots_in_interval, MultipleRootInfo
│   └── zero_manager.py         # ZeroManager orchestrator, integrate_segment, SegmentData
├── tests/                      # pytest: test_gbz_types.py, test_sgbz.py, test_amoeba.py,
│                               #   test_continuation.py, test_zero_manager.py
├── conftest.py
├── demos/                      # Runnable demo scripts (demo_unified.py shows the GBZResult API)
├── diagnostics/                # Debug scripts and analysis reports
├── data/                       # Pickled computation results (demos / replication)
├── log/                        # Change logs
│   ├── ......
├── doc/                        # Documentation
│   ├── SGBZ.md                 # SGBZ theory, architecture, API
│   ├── amoeba.md               # Amoeba theory, algorithm, API
│   ├── 拓扑匹配算法说明.md       # Topological matching algorithm
│   ├── sn-main.tex             # Simplified paper for SGBZ, main text
│   └── sn-supp.tex             # Simplified paper for SGBZ, supplementary information. Amoeba GBZ is discussed in section{Comparison with reported frameworks}
├── TODO/                       # Optimization checklists
└── README.md
```

## Development Principles

- **反常输出是改进算法的机会，不是需要绕过的 bug。** 遇到 unexpected behavior（如 winding 非单调、plateau check 误触发）时，优先追查根因并修复底层算法，而不是加 workaround（如调阈值、加预检查、限制探针范围等）。3cdb22d 版本没有 plateau check 时反而结果正确，说明问题出在后来引入的逻辑。对比不同 commit 是定位问题的有效手段。

- 注释应当讲述“为什么”，而不是“做什么”。

- Ask WHY before HOW. Check the FACT before reach the CONCLUSION. Never assert a bug before you get the solid evidence.

- 遇到与预期不符的结果时，要立刻汇报。
