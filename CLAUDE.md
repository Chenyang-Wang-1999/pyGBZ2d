# CLAUDE.md — brute-force-non-hermitian

Non-Hermitian skin effect computation for 2D tight-binding models. Two complementary modules implementing brute-force polynomial root-solving approaches. Both entry points (`collect_GBZ_subsets`, `collect_GBZ_subsets`) natively return the unified `GBZResult` type defined in `gbz_types.py`.

## Project Structure

```
brute-force-non-hermitian/
├── gbz_types.py             # Unified data types shared by both modules:
│                            #   PointSubset (0D), LineSubset (1D, lazy beta2_arr),
│                            #   GBZResult, ConnectedSubset = Union[Point, Line]
│                            # + shared utils: sort_by_root_abs, chordal_cost_matrix,
│                            #   hungarian_match_indices, find_cyclic_true_intervals,
│                            #   get_minor_degrees, generate_probe_steps, to_sphere_r3
├── brute_force_SGBZ/        # SGBZ / strip winding number formulation
│   ├── __init__.py             # Whitelist exports (re-exports gbz_types classes)
│   ├── root_solver.py          # calculate_point_roots, complex_root, poly_to_np_coefficients
│   ├── winding.py              # PolyDiffContext, WindingFun, MatWindingFun, get_winding_number
│   ├── pmgbz_detector.py       # get_roots_and_PMGBZ — root solving on theta1 mesh,
│   │                           #   continuum detection, accidental point refinement (~800 lines)
│   ├── strip_winding_number.py # get_strip_winding, get_loop_winding
│   └── SGBZ.py                 # solve_SGBZ_for_E, collect_GBZ_subsets (returns GBZResult)
├── brute_force_amoeba/        # Amoeba / Ronkin function formulation
│   ├── __init__.py             # Whitelist exports (re-exports gbz_types classes)
│   ├── tracks.py               # Root tracks via Hungarian matching (get_hungarian_sorted_roots)
│   ├── ronkin_winding.py       # Ronkin winding (get_a1/a2_average_winding, crossing detection)
│   ├── bisect.py               # Bisection + Newton refinement (bisect_amoeba_ronkin_min)
│   └── amoeba.py               # collect_GBZ_subsets (returns GBZResult), plateau check, orchestration
├── tests/                      # pytest: test_gbz_types.py, test_sgbz.py, test_amoeba.py
├── conftest.py
├── demos/                      # Runnable demo scripts (demo_unified.py shows the GBZResult API)
├── diagnostics/                # Debug scripts and analysis reports
├── data/                       # Pickled computation results (demos / replication)
├── log/                        # Change logs
│   ├── log-old.md              # log before git init
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
