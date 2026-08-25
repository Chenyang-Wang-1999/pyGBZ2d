# SGBZ pairwise：求交前网格加细（单区间多交点）

日期: 2026-08-25

## 背景

`collect_pair_events` 的扫描判据只有两条：

1. `d[:-1] == 0` —— 精确触零（只在已有 mesh 点上）；
2. `d[:-1] * d[1:] < 0` —— 区间两端异号。

它只能看到**奇数次**穿越。一个网格区间内若同一对代表元曲线有 2 个
交点（两端同号），或 3 个交点（两端异号但中间多一对），扫描会漏掉
偶数个，boundary pair 不在正确 θ₁ 处换边，μ₂_mid 不经过应有交点，
平均绕数随 μ₁ 出现反常跳变。

目标：在 `ZeroManager.run()` 之后、`collect_pair_events` 之前，对
“一个区间内可能含有多个交点”的区间自适应加细，使每个交点都落在
单独的子区间里。

## 设计

### 1. 可疑区间筛选（logabs 距离 + Hermite 探底）

对每个 segment 的每个区间、每对代表元 item，构造
`d = ln|β₂_a| − ln|β₂_b|` 的三次 Hermite 插值（端点值 + 端点导数）。
区间/ pair 满足以下任一条件即视为可疑：

- `d0 == 0` 或 `d1 == 0`（端点触零）；
- `d0 * d1 < 0`（已有奇数次穿越，仍可能还藏着偶次穿越）；
- 两端同号且 Hermite 曲线内部极值探近零（正端点向下探到
  `< tie_tol`，负端点向上探到 `> -tie_tol`）；
- 区间与 MR 边界行相邻（Hermite 锚在 MR 附近不可靠，保守加细）。

非有限 `d`（0/∞ padding 根）的 pair 主动跳过——这是 β₂=0 / β₂=∞
缺口的另一独立问题，不应毒化本流程。

### 2. 交点预测

对每个可疑 pair 求 Hermite 多项式在 `(0,h)` 内的实根。同一区间内
所有 pair 的预测根取并集，间距 `< crossing_tol` 的根合并为一个
cluster（下游 EventGroup 本来也会把它们合并）。**只有合并后根数
≥ 2 时**才制定加细计划；单根交给原有 sign-change + brentq。

### 3. 网格颗粒度

设合并后的根位置为 `p = [0, x₁, ..., xₙ, h]`（局部坐标），最窄间隙

```
Δ = min(p[k+1] − p[k])
m_sub = clamp(ceil(ρ · h / Δ), n + 2, 64),   ρ = 4
```

- 在 `(a,b)` 内插入 `m_sub` 等分网格；
- 每个相邻 `p` 间隙额外插入一个**分离中点**
  `(p[k] + p[k+1]) / 2` —— 位置无关的兜底，保证相邻预测根被隔开；
- 若 `ρ·h/Δ > 64`（根极近），不再提高均匀密度，改为把预测根本身
  也插入 mesh 作为 touch 锚点。

最窄间隙里至少放 4 个子区间，因此 Hermite 预测根即使有偏移，只要
偏移不超过约 `Δ/ρ`，真实根仍落在自己的子区间里。

### 4. 迭代与预算

- 每轮插入后重建 continuum cluster + ItemView，再重扫计划；
- 最多 3 轮；累计插入上限 2000 行；
- 超轮数或超预算后仍有残留计划时，做一次不设预算的
  “根 + 分离点”强制插入并 `warnings.warn`；
- 插入统一 `insert_solution(interp='hermite')`，每 segment 降序；
- 单点插入失败只警告，不打断分析。

## 修改

- `pairwise._RefinementPlan`：一个待加细区间的 seg/interval/局部根。
- `pairwise._find_refinement_plans`：步骤 1 + 2。
- `pairwise._refinement_grid_points`：步骤 3 的 θ 生成。
- `pairwise._insert_refinement_grids`：按 segment 去重、降序插入。
- `pairwise.refine_mesh_for_multiple_crossings`：迭代包装，返回插入行数。
- `pairwise._safe_deriv_diff`：inf/nan 导数差统一折叠为 `inf`，
  既保持 Hermite 线性 fallback 语义，又避免 `inf − inf` 的
  `RuntimeWarning`；`_direction_from_tangents` 同步复用。
- `mu2mid.Mu2MidZM.analyze`：初次 ItemView 之后调用 refine，之后重建
  cluster/ItemView 再进入 `collect_pair_events`。新增开关：
  `refine_multi_crossings=True`、`refine_max_rounds=3`、
  `refine_safety_factor=4.0`、`refine_max_subintervals=64`、
  `refine_max_total_inserts=2000`；`build_mu2_mid` alias 同步透传。

不改动 `collect_pair_events` / `group_events` / `finalize_event_groups`
的语义；它们只看到更细的 mesh。

## 测试

`tests/test_sgbz.py::TestPairwiseAnalysis` 新增：

- `test_multi_root_plan_predicts_two_interior_roots`：
  合成双根区间 → 预测根 `[0.25, 0.55]`；`ρ=4, Δ=0.25 → m_sub=16`；
  分离中点存在。
- `test_multi_root_plan_hits_subinterval_cap`：
  `ρh/Δ ≫ 64` 时均匀网格封顶，预测根本身成为 touch 锚点。
- `test_roots_closer_than_crossing_tol_merge_away`：
  过近预测根合并后不触发加细。
- `test_nonfinite_pair_diff_is_skipped_by_planner`：
  0/∞ padding pair 不崩溃、不产计划。
- `test_refine_isolates_two_real_crossings_before_pair_scan`：
  真实 M=2 Laurent 模型，粗化到单区间含两个真实交点；
  refine 前 `collect_pair_events` 0 个事件，refine 后精确检出
  `θ = arccos(−0.98)` 与 `2π − arccos(−0.98)` 两个事件。

全量 `pytest tests/`：211 passed, 3 skipped。

## 目标场景验证

Gain-loss Haldane `E = 1.289 − 0.4131j`，μ₁ 扫过
`-0.0002 → -0.0001 → 0 → 0.0001 → 0.0002`：
平均绕数从 `-2.27e-4` 平滑过零到 `+2.27e-4`，不再跳变；
`collect_GBZ_subsets(..., debug_mode=True)` 返回
`success=True, is_gbz=True, index=(8, 0)`。
