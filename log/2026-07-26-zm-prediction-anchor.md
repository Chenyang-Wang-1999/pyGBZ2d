# ZeroManager 匹配改用导数预测锚 — 会议记录

日期：2026-07-26
模块：`continuation/`（zero_manager, arclength, multiple_roots）

## 议题

`arclength_step` 内部用切向预测 `predicted → solved` 做匹配锚，避免最近接近点换轨。但 `ZeroManager` 中 `_to_track_order` 与 `run` 的若干处退化成"两坨已解根之间裸匈牙利"，在换轨点会出错。本轮把 ZeroManager / multiple-root 求解器里所有匹配点统一改成导数预测锚，并修复三个设计缺陷。

## 决策与设计

### 核心原则
匹配锚 = 端点值 + 端点导数，与 `arclength_step` 同源：
- 端点值：已解 `roots`（track-ordered，同 frame）。
- 端点导数：`dβ₂/dθ₁ = V_j · β₂ⱼ`，`V = compute_tangent(poly, E_ref, exp(μ₁+iθ), roots)`。
- 匹配约定沿用 `arclength_step`：`perm = hungarian_match_indices(predicted, roots_target)`，`roots_target[perm]` 落在预测（segment frame）的序上。

### `predict_roots_hermite`（arclength.py 新增）
- 双端（值+导数两端给定，目标在区间内）→ 三次 Hermite，复数算术，唯一解。
- 单端 → `β·exp(V·Δθ)`，复用 `predict_roots`，与 arclength 一致。
- 回退链（逐 track）：Hermite 奇异（`_is_singular_root` 或非有限）→ 两点 lerp `p0+s(p1-p0)` → lerp 退化（同点）→ hold fixed。单端奇异 → hold fixed。与 arclength 对 singular root 一致。

## 三个问题与修复

### 问题 1：MR 锚点不能用 segment 最后两点

`_refine_mr` 原方案用 `new_seg_theta1[-2]` / `[-1]` 做 Hermite，但这俩点不保证跨过重根所在区间。应按 MR 类型分别取锚：

- `multiple_root_in_interval` → 用区间端点 `mr_interval_start` / `mr_interval_end`。两端是连续两步的正则点，同 frame。
- `multiple_root_encountered` → 用 `mr_approx_theta`（步长坍塌处）单端外推。

**实现**：让 `integrate_segment` 在返回时顺手带上 MR 端点的 (θ, roots, V)：
- 新增 `_MREndpoint` NamedTuple（θ + roots + V）。
- `SegmentResult` 加 `mr_ref`（MR-相邻正则端，两 MR 路径都设）和 `mr_ref2`（区间另一端，仅 in_interval）。
- `MultipleRootIntervalTrigger` 现在留存 `_prev_roots` / `_prev_V`，触发时 early-return（不覆盖 prev 状态），保证区间左端点的 roots+V 拿得到。
- `_to_track_order` 加 `ref_V` / `ref_V2` 参数，直接复用 `_MREndpoint` 预存的切向，不重算。

### 问题 2：completed / boundary-MR 的 2π 预测不能假设 [-2] < 2π < [-1]

原方案假设 `completed` 的 `theta1_arr[-2] < 2π < theta1_arr[-1]`，但下标关系不保证。鲁棒做法：

新增 `_predict_roots_at_2pi(theta1_arr, tracked_roots)`：
- 找 `theta1 < 2π` 的最后一行 `idx`。
- 若 `idx+1` 存在（跨 2π）→ 双端 Hermite（`idx` 与 `idx+1` 同 frame）。
- 否则（单行 segment 或 segment 末到 2π）→ 单端外推。

completed 与 boundary-MR 两条路径都用它，消掉重复代码。singular track 自动走 lerp/hold 回退。

### 问题 3：两处 `theta < 2π` 条件重复

`run` 的 `while theta < 2*pi` 和 `integrate_segment` 内的 `while theta1 < theta_end` 是重复判断，且前者需要 post-loop fallback 兜底。

- `run` 改为 `for _ in range(_MAX_SEGMENTS)`（防死循环），循环终止只靠 `completed` / boundary-MR 两个 break 分支。
- 删除 post-loop 的 `boundary_perm` fallback，改用 assert 保证不变量。
- `for` 的 `else` 在耗尽时 raise，明确报错而不是静默产出半成品 topology。
- 两处 `theta < 2π` 去重为一处（`integrate_segment` 内），逻辑统一。

## 顺带修复的潜在 bug

`compute_tangent` 在 `df_dbeta2 == 0`（多重根处，`∂f/∂β₂ = 0`）会 `ZeroDivisionError`。其 docstring 声称"caller 的 min_dtheta/min_step guards 会捕获发散"，但除零直接崩，guards 从未运行。现置 `V_j = 0`（切向未定义），让下游步长坍塌 / 预测 hold-fixed 自然处理。符合 CLAUDE.md "追根因、不加 workaround" 原则。

## 改动文件
- `continuation/arclength.py`：新增 `predict_roots_hermite`；修 `compute_tangent` 除零。
- `continuation/zero_manager.py`：新增 `_MREndpoint`、`SegmentResult.mr_ref/mr_ref2`；`_to_track_order` 加 `ref_V/ref_V2`；新增 `_predict_roots_at_2pi`；`integrate_segment` 收集并返回 MR 端点；`_refine_mr` 改用 `seg.mr_ref/mr_ref2`；completed/boundary 改用 `_predict_roots_at_2pi`；run 循环 `while→for` + 删 fallback。
- `continuation/multiple_roots.py`：`MultipleRootIntervalTrigger` 留存 prev roots/V + early-return；Brent 内部（`solve_multiple_roots_in_interval`）改预测锚（`V_ref` 在 right 端点预算一次复用）。
- `tests/test_continuation.py`：新增 `TestPredictRootsHermiteSynthetic`（3 例：三次 track 精确、singular root hold-fixed、单端等价于 predict_roots）。

## 验证
- `tests/test_continuation.py` + `tests/test_zero_manager.py`：68 passed（含 3 个新 Hermite 单元测试）。
- `tests/test_amoeba.py`：通过。
- `demos/demo_zm_gbz.py`：输出与改前完全一致（4 个用例的 seg 数 / MR 数 / index 不变）。

## 遗留
- `tests/test_sgbz.py` / `tests/test_gbz_types.py` 部分失败，已确认是 session 开始前 amoeba/gbz_types 在制改动所致，与本轮无关。
- 边界 MR 处退化簇的 cluster↔cluster 对应本质模糊，预测锚收益有限，但不恶化。
- Brent 内部改预测锚可能轻微改变 `theta1_mr`（匹配→导数→Brent 路径），demo 未观察到行为变化，仍建议回归对比。
