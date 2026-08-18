# 2026-08-15 — Cubic Hermite 插值统一到 continuation.interpolation

用户审查发现全项目至少三处独立实现两节点 cubic Hermite 插值系数构建，且
发散端点的回退逻辑各自为政。本轮将其统一（方案经用户逐点审批）。

## 新模块

`continuation/interpolation.py`：

- `cubic_hermite_poly(h, v0, dv0, v1, dv1)` — 纯三次 Hermite 系数生成，
  不做任何输入检测。返回 numpy poly 降幂格式 `[a, b, c, d]`。
- `hermite_interp_poly(...)` — 插值多项式系数生成：检测 dv0/dv1 是否有限；
  都有限 → 三次（len 4）；任一非有限 → 线性 `[(v1-v0)/h, v0]`（len 2，
  numpy poly 格式）。不返回额外标志，调用方用 `len(poly)` 区分。

统一约定：求值 `np.polyval`，求导 `np.polyder`，求两条曲线交点
`np.roots(np.polysub(p1, p2))`——不再为 linear/cubic 写分支。

`continuation/__init__.py` 导出两个新函数。

## 迁移点

1. `brute_force_SGBZ/mu2mid.py`
   - 删除本地 `_cubic_hermite_coeffs` / `_cubic_roots_in_interval`；
   - `_cubic_hermite_iterate`：`hermite_interp_poly` + `np.roots` /
     `np.polyder`，区间实根筛选与 transversal 选择保留；
   - `_Mu2MidPath.value_deriv`：内联三次公式替换为统一函数 +
     `np.polyval` / `np.polyder`；`h <= 0` 短路和 ±14 clamp 保留；
     `s` 定位变量随之删除；`import math` 移除。
2. `brute_force_SGBZ/crossings.py`
   - `_bracket_crossing`：改为统一函数 + `np.roots` / `np.polyder`；
     原“先手动把 dv 换成 secant”的预处理删除，线性回退由
     `hermite_interp_poly` 承担。
3. `continuation/arclength.py`
   - `predict_roots_hermite` 两端点分支由基函数求和改为
     `hermite_interp_poly(dt, p0, V0·p0, p1, V1·p1)` + `np.polyval`；
     singular-root hold 与求值结果回退链保留。

## 测试

- 新增 `tests/test_interpolation.py`（12 例）：纯三次精确性、端点值/导数
  匹配、有限导数等价性、dv0/dv1 分别为 inf/-inf/nan 时的线性回退、
  复数值支持、`np.polysub` + `np.roots` 交点工作流。
- 回归计划：test_continuation（重点 TestPredictRootsHermiteSynthetic 精度）、
  test_sgbz、test_zero_manager、test_amoeba。

## 文档

- `doc/continuation.md`：新增 §2.8 Unified Hermite interpolation，§3.2 API
  表与 File Layout 同步。
- `doc/SGBZ.md` §2.1：新增 Unified Hermite interpolation 段落。

## 已知风险

`predict_roots_hermite` 从 Hermite 基函数求和改为 poly 求值，数学等价但
浮点舍入路径不同；用户已确认“不回退最好，但回退也没关系”。
