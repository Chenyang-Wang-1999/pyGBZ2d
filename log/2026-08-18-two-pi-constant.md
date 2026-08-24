# 2π 常量统一：gbz_types.TWO_PI

日期：2026-08-18

## 动机

代码中 `2 * pi` 有多种拼写（`2 * math.pi`、`2 * np.pi`、`2 * cmath.pi`、
`2.0 * math.pi`、局部 `twopi = ...`）。当前它们恰好全部 bit-identical
（`6.283185307179586`，已核实含 `math.tau`），但分散拼写无法防止将来
分叉——seam 比较（θ % 2π、circ_dist、闭合行 θ=2π、boundary-MR 容差）
必须操作同一个浮点数。

## 修改

- `gbz_types.py`：定义 `TWO_PI: float = 2.0 * math.pi`（带注释说明目的）。
- 全项目代码行 `2 * pi` 类拼写 → `TWO_PI`，共 70 处：
  - 核心库 39 处（gbz_types 5、zero_manager 8、multiple_roots 1、
    winding 6、pairwise 2、plateau 1、zm_extract 10、ronkin_winding 4[+2 局部 twopi]）；
  - tests 19 处（6 个文件）；
  - demos / diagnostics / simple-plot 25 处（8 个文件）。
- 注释与 docstring 中的 Unicode π 未动。
- 顺手清理因此变成死 import 的 `pi`：`multiple_roots.py`、
  `ronkin_winding.py`（其余文件的 pi 仍有 pi/2 等使用者，保留）。

## 验证

- 纯等价重构：替换前核实所有变体 bit-identical。
- 23 个触及文件 `py_compile` + 库 import 通过；代码行残留扫描为零。
- `pytest tests/`：188 passed, 2 skipped（与改动前完全一致，零行为变化）。
