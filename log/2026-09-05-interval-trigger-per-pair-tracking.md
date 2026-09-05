# 区间触发器逐对追踪与零点重合分组修复

日期: 2026-09-05

## 背景

`playground/repro_hermitian_amoeba_extra_points.py`（厄米 Haldane 模型，
E = 3.5/4.0/5.0）：continuation 的重根检测漏掉 θ₁ ≈ 4.578749 处的**同时
双重简并**（两对根在同一 θ₁ 各自合并），amoeba 因此返回 44/70/46 个多余
PointSubset（SGBZ index (0,4)，amoeba index (44,4)/(70,4)/(46,4)）。

排查过程（worktree A/B 对照 f05d7de 前后版本）：

- f05d7de 的 same-MR 去重守卫**无罪**——修改前后输出逐行一致，两版都漏。
- 逐行回放（对段内每行重算切向、重放触发器）定位机制：两对根距离全程
  打平到 3 位有效数字，*argmin* 最近对身份逐行交替；区间触发器的
  same-pair guard 恰好在导数符号翻转的那一行把翻转拒绝掉。
  点触发也不响：步长在简并点附近只塌缩到 ~5e-9 > MIN_DTHETA = 1e-10，
  （chordal 距离 ~C·√|Δθ|，C≈9，落到 CLUSTER_TOL=1e-4 需要
  |Δθ| ≲ 1e-11，踩不中）。
- 独立逐点解根扫描确认两个 dip 都是真实双重根（非近距擦过）。

## 设计讨论与决策（与王辰阳确认后定稿）

1. **向量化**：距离与 d|βᵢ−βⱼ|²/dθ₁ 广播一次算出所有对
   （`_pairwise_dist_deriv`），`_closest_pair_deriv` 退化为薄封装。
2. **逐对追踪**：触发器状态从"单一 argmin 对 + same-pair guard"改为
   n×n 逐对导数符号矩阵；武装条件逐对判定（该对距离 < 
   `MIN_DIST_THRESHOLD`）。一次停机返回 `list[MRTriggerRecord]`
   （pair + bracket + 诊断量），多对同时翻转返回多条。
3. **逐对二分**：`solve_multiple_roots_in_interval(min_pair=...)` 的
   g(θ) 改为该对自身的导数（连续，不随 argmin 换人跳变）；
   删除原 argmin 身份不符即抛 ValueError 的检查。
4. **事件分组 = 零点重合检验**（替代最初提议的 θ 容差常量
   `MR_SAME_THETA_TOL`）：候选 θ 归入已有事件当且仅当该候选的触发对
   在事件 θ 处 chordal 距离 < `CLUSTER_TOL`。隐含 θ 分辨率
   ~(CLUSTER_TOL/C)² ≈ 1e-11，介于 brentq 噪声与物理间距之间，
   不引入新魔法数。同时简并的不同零点对合并为**一条多 cluster 的
   MultipleRootInfo**（与 5.229 处 (0,1),(2,3) 的既有形态一致）。
5. **相近不同零点事件之间的段**：不做重启重检测，直接在两事件间
   **密集采点**（均匀网格，间距 ≤ `MR_DENSE_MAX_STEP` = 1e-4，内部点数
   ≥ `MR_DENSE_MIN_SAMPLES` = 8），段形如 [左 MR 行 | 常点行 | 右 MR 行]，
   与普通段同构。最后一个事件后重启点仍为 bracket 右端 `ref.theta`。
6. `MIN_DTHETA` 恢复 1e-10（1e-7 实验被本修复取代）。
7. 暂不加端到端回归测试（后续如需，候选方案：把厄米 Haldane 的
   (coeffs, degs) 缓存为 npz 供测试直接构造 CharPoly，绕开 BerryPy 依赖）。

## 修改

`src/pygbz2d/continuation/multiple_roots.py`

- 新增 `_pairwise_dist_deriv(roots, V)`（含 nan/inf 切向剔除，语义同
  旧标量循环；tie-break 顺序一致）。
- 新增 `MRTriggerRecord` NamedTuple；`MultipleRootIntervalTrigger`
  重设计为逐对追踪，`__call__` 返回记录列表。
- `solve_multiple_roots_in_interval`：`min_pair` 给定时以
  `_pair_distance_deriv` 逐对求 g(θ)。

`src/pygbz2d/continuation/zero_manager.py`

- 新常量 `MR_DENSE_MAX_STEP = 1e-4`、`MR_DENSE_MIN_SAMPLES = 8`；
  `MIN_DTHETA` 恢复 1e-10。
- `SegmentResult`：`mr_interval_start/end` 标量对替换为
  `mr_triggers: tuple[MRTriggerRecord, ...]`。
- 新增 `_MREvent` / `_MRRefine`；`_refine_mr` in_interval 分支重写为
  `_events_from_triggers`（逐对 brentq + 零点重合分组 + 假阳性丢弃）
  与 `_dense_rows_between`（事件间密集采样段）。
- run() MR 分支循环处理多事件（逐事件 append MR 与段；边界 2π 钉扎
  适用于落入 BOUNDARY_THETA_TOL 的事件；f05d7de 守卫与 forward-progress
  守卫检查首个事件）。
- 新增 `_chordal_dist`（零点重合检验的度量，与 detect_cluster 一致）。

`src/pygbz2d/continuation/__init__.py`：白名单导出 `MRTriggerRecord`。

`tests/test_continuation.py` / `tests/test_regressions.py`

- 触发器测试更新到 list 返回契约；新增：向量化 vs 标量奇偶校验
  （随机数组含 nan/inf 切向）、双对同时翻转、打平距离下 argmin 抖动
  的回归形状（旧设计不触发、新设计触发）、逐对武装/复位。

`doc/continuation.md`、`doc/constants.md`：触发器机制、事件分组与
两个新常量的文档。

`playground/repro_hermitian_amoeba_extra_points.py`：docstring 的
预期输出更新为修复后观测。

## 验证

- E = 3.5/4.0/5.0 三个能量：ZeroManager 均检出 6 个 MR（含原先漏掉的同时
  双重简并，如 E=3.5 的 4.578749345, clusters [(0,2),(1,3)]），
  其余 MR θ 逐位不变；段以真实积分行缝合。
- 下游：amoeba 与 SGBZ 三个能量全部一致为 index (0,6)、各 6 个
  LineSubset、0 个 PointSubset——repro 的多余点集问题消失。
- 全范围逐点扫描（步长 2e-4）确认每个 chordal 极小 < 0.1 的 dip 均被
  MR 覆盖；扫描中另有两组"最近而未合并"的 dip（chordal ≈ 0.70/0.24，
  θ ≈ 1.2 与 3.2 附近，随 E 缓移）非重根，不检出为正确行为。
- `pytest`：无新增失败；8 个失败经 git stash 对照确认为先例
  （Windows SIGALRM ×4、GBK 读文件解码 ×3、sgbz pairwise ×1）。

## 后续事件记录

- 修复后一度误报"最右侧重根仍未检出"：实为在旧版本代码上复跑；
  在当前树复跑同一脚本确认三能量均正常。
