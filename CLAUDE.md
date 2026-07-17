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
│   ├── sn-main.tex             # Simplified paper for SGBZ, main text
|   ├── sn-supp.tex             # Simplified paper for SGBZ, supplementary information. Amoeba GBZ is discussed in section{Comparison with reported frameworks}
│   └── amoeba.md               # Amoeba theory, algorithm, API
├── demos/                      # Runnable demo scripts
├── pyproject.toml
└── README.md
```

## Development Principles

- **反常输出是改进算法的机会，不是需要绕过的 bug。** 遇到 unexpected behavior（如 winding 非单调、plateau check 误触发）时，优先追查根因并修复底层算法，而不是加 workaround（如调阈值、加预检查、限制探针范围等）。3cdb22d 版本没有 plateau check 时反而结果正确，说明问题出在后来引入的逻辑。对比不同 commit 是定位问题的有效手段。

