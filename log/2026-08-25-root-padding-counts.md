# CharPoly：按缺失次数补齐 0/∞ padding roots

日期: 2026-08-25

## 背景

`CharPoly.solve_roots_1d` 原来只在 partial polynomial degree 亏缺时各补一个
`β₂ = 0` 和 `β₂ = ∞`：

```python
if len(curr_roots) < M + N:
    if deg_M < M:
        curr_roots = np.append(curr_roots, 0)
    if len(curr_roots) - deg_M < N:
        curr_roots = np.append(curr_roots, np.inf)
```

这只能处理“低阶/高阶各最多消失一次”的情形。参数取特殊值（例如
`β₁ = 0`）可以同时消掉多个分子/分母项，导致返回数组短于 `K = M+N`，
下游 `SegmentData` / tangent / ItemView 的固定 `K` 形状假设随之崩掉。
partial polynomial 完全消失时还会在 `max(degs_list)` 上直接抛
`ValueError: max() arg is an empty sequence`。

## 修复

在 `gbz_types.CharPoly.solve_roots_1d` 中显式计算缺失数量：

- 实际 denominator degree：`deg_M`（partial terms 为空时取 0，而不是
  poly_tools 保留的全局 denominator order）；
- finite polynomial degree：`max_deg`；
- `n_zero = max(0, M - deg_M)`；
- `n_inf = max(0, deg_M + N - max_deg)`；
- 检查 `len(finite_roots) + n_zero + n_inf == M+N`，不满足立即
  `RuntimeError`，不再静默返回短数组；
- 一次性追加全部 `0` / `∞`。

完全消失的 partial polynomial 返回 `M` 个 `0` + `N` 个 `∞`。

## 回归

`tests/test_gbz_types.py::TestCharPolyRootPadding` 覆盖：

- 高阶亏缺 2 → 两个 `∞`；
- 低阶亏缺 2 → 两个 `0`；
- 双侧同时亏缺 3+3 → 三个 `0` + 三个 `∞`；
- partial polynomial 完全消失 → 不再 `max([])` 崩溃，仍返回 K 个
  padding roots；
- 未亏缺时不额外 padding。

相关测试：`test_gbz_types.py` 30 passed；
`test_gbz_types.py + test_continuation.py + test_regressions.py` 93 passed；
全量 `pytest tests/` 216 passed, 3 skipped。
