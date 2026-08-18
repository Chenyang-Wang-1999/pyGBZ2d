# 2026-08-15 — SGBZ review fixes B.1–B.5 + hard-boundary charge=None

本轮回合：`brute_force_SGBZ` 代码审查中发现 B.1–B.5 五类问题，用户逐项
审批后实施；审查过程中用户进一步指出 hard 边界的拓扑荷是**未知**而非
0，因此一并把占位值从 0 改为 `None`（误用时直接 `TypeError`，不静默）。
修改前基线：`tests/test_sgbz.py` 20 passed（约 66 s）。

## B.1 — `_cubic_hermite_iterate` 未收敛不再报 `converged=True`

`brute_force_SGBZ/mu2mid.py`：迭代循环 `max_iter` 耗尽时原来返回
`(midpoint, True)`，`build_mu2_mid` 的 `not converged` 分支永远不触发。
现在返回 `(midpoint, False)`；docstring 明确 `converged=False` 表示
bracket 未达到 `xtol`。数值行为不变（仍返回最佳已知中点，已插入的逼近
行不动），只修正标志语义。

## B.2 — 右端 bracket 统一扩展循环，continuum proxy 不再绕过符号检查

`brute_force_SGBZ/sgbz_solver.py` Step 1.2：原来"左展开已建立
`mu1_ext_right`"的路径在 `handle_continuum` 修正 `(mu1_ext_right, w_high)`
后直接进入 bisection，不检查 `w_high` 符号——同号 bracket 会耗尽
`max_iter` 报错，`w_high == 0` 也不会走 `right_endpoint_zero`。

现在两条路径合并为一个统一循环：

- continuum proxy 修正后 `continue`，下一轮**重新评估修正点**（不会被
  通用 `+1` 步跳过）；
- 修正点与普通 winding 一样过 `> -zero_tol` 检查，不满足则继续向右扩展；
- 右端点 `abs(w_high) <= zero_tol` 时返回 `right_endpoint_zero`。

## B.3 — bracket expansion 上限

`sgbz_solver.py` 新增 `_MAX_BRACKET_EXPANSIONS = 10`（与 amoeba 的
`max_range_expansions` 对齐）。左/右扩展各用 `for ... else raise
RuntimeError` 替代 `while True`，异常模型（M=0/N=0 等不可解情形、W 恒
NaN/恒正）从死循环变为明确报错；`collect_GBZ_subsets` 默认捕获并返回
`success=False`，`debug_mode=True` 原样抛出。

## B.4 — `_ensure_mu2mid` 统一 continuum 阈值（方案 A）

`brute_force_SGBZ/crossings.py`：`_ensure_mu2mid` 构建新 `Mu2MidZM` 时
原来调用 `build_mu2_mid()`（默认 `_TIE_TOL_DEFAULT = 1e-8`），与 solver /
`detect_continuum_simple` 的 `tie_tol=CONTINUUM_TOL`（1e-6）不一致，同一
μ₁ 经不同入口可能得到不同 continuum 判定。现在统一传
`tie_tol=CONTINUUM_TOL`。`build_mu2_mid` 自身默认值不改（避免连带改变
Stage 0 wall 的 near-tie 阈值）。

## B.5 — 纯软边界才做 charge 守恒校验

`brute_force_SGBZ/winding.py::compute_average_winding`：边界列表构造后、
seed/quad 之前，仅当**全部边界都是 soft**（ordinary，charge ±1）时检查
Σcharge == 0；非零说明 crossing 检测漏检/多检零点，raise `RuntimeError`
（不再静默平均不兼容的区间 winding）。存在任何 hard 边界（charge 未知）
则跳过校验。

## Hard-boundary charge 从 0 改为 None

用户指出 hard 边界的 charge 是**未知**，不是 0；原代码/注释多处用
"charge 0" 描述，自相矛盾。现已全面清理：

- `crossings._classify_charge`：mr/tangent 返回 `charge=None`；
  docstring 说明 `None` 是未知，不是占位数值，误用会 `TypeError`。
- `crossings.detect_crossings_simple` 的 MR charge dict：`charge=None`。
- `winding.compute_average_winding`：`boundaries` 类型
  `list[tuple[float, bool, int | None]]`；注释明确 soft 的 dc 是真实
  winding jump，hard 的 dc=None 只作 region 分隔、绝不参与传播。
  传播路径 region 内部必为 soft，正常不会触碰 None；未来误取会直接
  `TypeError`（fail-fast）。
- `doc/SGBZ.md` §2.3/§2.4/§3.4/§3.6 同步，全文清除 "charge 0" 表述。

## Tests

`tests/test_sgbz.py` 新增 `TestReviewFixes`（8 个用例）：

- `test_cubic_hermite_iteration_exhaustion_reports_not_converged`
- `test_left_bracket_expansion_raises_after_cap`
- `test_right_bracket_keeps_expanding_after_same_sign_continuum_proxy`
- `test_right_endpoint_zero_from_continuum_proxy`
- `test_ensure_mu2mid_uses_continuum_tol`
- `test_hard_charge_is_none_sentinel`
- `test_compute_average_winding_rejects_nonconserved_charges`
- `test_compute_average_winding_skips_none_hard_charge`

## Verified

- `python -m pytest tests/test_sgbz.py::TestReviewFixes -q` — 8 passed。
- 回归计划：`tests/test_sgbz.py` 全量、`tests/test_zero_manager.py`。
