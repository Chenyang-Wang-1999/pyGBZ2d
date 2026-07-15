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
├── demos/                      # Runnable demo scripts
├── pyproject.toml
└── README.md
```

