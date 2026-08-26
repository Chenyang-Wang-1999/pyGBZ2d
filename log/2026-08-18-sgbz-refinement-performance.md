# SGBZ 性能优化：refinement 幅度守卫 + 段级增量 + ±14 clamp 守卫

日期: 2026-08-18

## 背景：profile 结论

对 SGBZ 计算代码做了 cProfile 剖析（负载：HN2D 链对 E=1.0 谱内 /
E=5.0 谱外；Haldane 增益损耗 2×2 超胞 E=1.212），产出
`diagnostics/profile_sgbz.py`（可复跑脚本）、
`diagnostics/profile_sgbz.prof`（原始 dump）、
`diagnostics/profile_sgbz_report.md`（报告）。

求解器净耗时 105.3 s（剔除 Haldane demo 顶层 import matplotlib 的
~25 s 字体缓存噪声）。热点归因：

| 阶段 | cumtime | 占比 |
|---|---|---|
| `refine_mesh_for_multiple_crossings` | 40.0 s | 38% |
| `ZeroManager.run`（积分） | 33.6 s | 32% |
| `build_mu2_mid`（±14 clamp 切分） | 13.4 s | 13% |
| 绕数积分 `detect_crossings_and_winding` | ~13.6 s | ~13% |

内核级：`np.roots` 共 83,478 次调用、30.1 s。按调用方拆解：**二次式
55,896 次（67%）**来自 `_hermite_dips_near_zero`（dip 检验求三次
Hermite 导数式 —— 一个二次式 —— 的根，每次付出伴随矩阵 + LAPACK
eigvals 固定开销）；**三次式 17,460 次（21%）**来自
`_split_piece_at_clamp`（±14 clamp 切分，每区间 2 次）；真根
`solve_roots_1d` 仅 10,122 次（12%）。

即 **88% 的 np.roots 调用在次数 ≤3 的小多项式上**。顶层还有结构性
放大：μ₁ 二分每次探针完整冷重建 `Mu2MidZM`（HN2D E=1.0 重建 33 次 ×
~1.2 s）。

与作者讨论后确定实施范围：① dip 检验前置幅度守卫；③ refine 段级
增量重扫。方向②（小多项式闭式解）不动；方向④（跨 μ₁ 探针暖启动，
需 ZeroManager 大改）写入 TODO（见下）。

## 设计

### 1. dip 检验前置幅度守卫（`pairwise._find_refinement_plans`）

同号对（`d0·d1 > 0`，占绝大多数 —— 两条曲线离得很远）原本全部进入
数值 dip 检验。现利用三次 Hermite 的 **Bernstein 控制点凸包性质**：
控制点

```
B = (d0, d1, d0 + h·m0/3, d1 − h·m1/3)
```

的凸包包含曲线全部值域。控制点全在 `±tie_tol` 带之外的同号侧 ⇒
曲线不可能探近零 ⇒ dip 检验必为 False，直接跳过。**零假阴性**：
与 `_hermite_dips_near_zero` 的触发条件（内部极值 `v < tie_tol` /
`> −tie_tol`）严格一致，被跳过的对不可能满足该条件。

非有限斜率的对（`hermite_interp_poly` 退化为线性 `[slope, d0]`）同号
端点间不可能 dip，同样跳过，与原数值路径返回 False 等价。

变号 / touch / MR 相邻对走原有 `suspicious` 廉价分支，守卫不生效。

### 2. refine 段级增量（`pairwise.refine_mesh_for_multiple_crossings`）

原实现每轮全量重建所有段的 clusters / ItemView 并全量重扫。新成本
模型：clusters / ItemView / plans 是段的纯函数，一次插入只改一个段。
故：

- 首轮复用 `analyze()` 刚建好且与网格同步的 views（行数一致性检查
  `_views_in_sync_with_mesh`），消除与 analyze 的重复重建；
- `_insert_refinement_grids` 返回 `(inserted, modified_segs)`；
- 后续轮只重建、重扫实际收到插入行的段；未重扫段的 plans 原样保留
  （其网格未变，重扫结果恒等）；
- 终止条件（无 plans / 零插入）不变，无死循环风险。

为支持上述复用，`mu2mid.py` 提取单段内核
`_detect_continuum_clusters_for_segment` /
`_build_item_view_for_segment`；`_build_item_views` wrapper 通过类名
显式分发保持 unbound 调用兼容（测试用 SimpleNamespace 传入）。

### 3. ±14 clamp 切分守卫（`mu2mid._split_piece_at_clamp`）

第二批实施（发现过程：优化后 profile 复核显示剩余 np.roots 中
`_split_piece_at_clamp` 占 59.8%，插桩实测两个负载共 8,601 个 piece
**100% 可被守卫跳过，0 次真实切分，0 个 ±14 常值段** —— μ₂_mid 是
边界对 log 模长均值，典型量级 O(0.x)，而 e¹⁴ ≈ 1.2×10⁶，该机制只在
MR 两侧边界根发散的区间才有意义。作者确认后实施）。

与 dip 守卫同源：控制点全在 `[−14, +14]` 带内 ⇒ 两次求根可证扑空，
跳过。两个细节：

- **闭不等式安全**：`_boundary_mean` 已把端点值 clip 到 ±14，端点可
  恰好触界；即使内部切线触界被漏掉也是输出 no-op —— 切分后按端点
  值/导数重拟合 reproduces 同一条三次式（三次 Hermite 唯一性）；
- 提前返回路径与原无切分路径**逐位一致**（同一 `_make_piece` 全区间
  重拟合）。

## 修改

- `brute_force_SGBZ/pairwise.py`
  - `_find_refinement_plans`：同号对 dip 检验前 Bernstein 幅度守卫；
    新增 `seg_filter` 参数（段集合过滤）。
  - `_insert_refinement_grids`：返回值 `int → (int, set[int])`
    （插入行数 + 被修改段集合）。
  - `refine_mesh_for_multiple_crossings`：首轮复用同步 views、后续
    轮段级增量重建/重扫、未重扫段 plans 保留。
  - 新增 `_views_in_sync_with_mesh` 辅助。
- `brute_force_SGBZ/mu2mid.py`
  - 提取 `_detect_continuum_clusters_for_segment` /
    `_build_item_view_for_segment` 单段内核；
    `_detect_continuum_clusters_internal` / `_build_item_views` 变为
    逐段委托的 wrapper。
  - `_split_piece_at_clamp`：±14 求根前 Bernstein 幅度守卫。
- `TODO/warm-start-zeromanager-across-mu1-probes.md`（新）：方向④
  的英文设计文档 —— 继承网格 → 重求根 → 全局 Hungarian 匹配到旧
  track 帧 → 自适应精化 + MR 检测的四步管线、输出契约约束（MR 漂移 /
  continuum 跳变回退冷构建 / track 错配风险）与预期收益。
- `diagnostics/profile_sgbz.py` / `profile_sgbz.prof` /
  `profile_sgbz_report.md`（新）：profile 脚本、dump 与报告（报告含
  两批优化的实测附录与测量噪声警示）。

不动 `collect_pair_events` / `group_events` / `finalize_event_groups` /
`build_mu2_mid` 的语义。

## 测试

- `tests/test_sgbz.py`（48）、`tests/test_zero_manager.py` /
  `test_continuation.py` / `test_interpolation.py`、
  `tests/test_counterexamples.py` / `test_regressions.py` /
  `test_gbz_types.py`：全部通过（65 passed + 2 slow-skipped）。
- 合成双根用例（`d(x) = 4(x−0.25)(x−0.55)`）控制点
  `(0.55, −0.517, −0.25, 1.35)`，min < tie_tol，守卫正确放行 ——
  dip 检验的多交叉保护能力不变。
- 三负载结果不变：HN2D E=1.0 index=(0,2)、E=5.0 index=(0,0)、
  Haldane E=1.212 index=(0,6)。

## 性能验证

**测量噪声警示**：本机负载波动大，同样代码数分钟内 HN2D E=1.0 实测
32.8–39.1 s。单次墙钟只作量级参考；可信口径为负载无关计数
（np.roots / dip 调用数）与同进程交错 A/B。

负载无关计数（profile dump）：

| 指标 | 优化前 | 第一批后 | 第二批后 |
|---|---|---|---|
| `np.roots` 调用 | 83,478 | 29,706 | **12,246** |
| 其中：dip 二次式 | 55,896 | 1,837 | 1,837 |
| 其中：clamp 三次式 | 17,460 | 17,460 | **0** |
| 其中：真根 | 10,122 | 10,122 | 10,122（不可约） |
| `_find_refinement_plans` cum | 35.2 s | 3.4 s | — |
| `_split_piece_at_clamp` cum | 12.45 s | 12.45 s | **5.49 s** |
| `analyze` cum | 25.4 s | 25.4 s | **18.4 s** |
| `solve_SGBZ_for_E` cum | 105.3 s | 70.5 s | **66.6 s** |

（第二批后 12,246 = 真根 10,122 + dip 1,837 + planner 三次式 287，
对账闭合；clamp 剩余 cum 为 polyder/polyval/`_make_piece` 的正当
构建开销。）

同进程交错 A/B（交替执行守卫版/原始版，消除机器负载漂移）：

| 对照 | guard | orig | 加速比 |
|---|---|---|---|
| 第二批 clamp 守卫，HN2D E=1.0 | 32.8 s | 39.1 s | 1.19× |
| 第二批 clamp 守卫，Haldane | 17.8 s | 18.9 s | 1.06× |

（第一批的 39.6→22.0 s Haldane 单次墙钟含机器低谷成分，仅作量级
参考。）

## 遗留

剩余瓶颈排序：μ₁ 冷重建乘数（TODO 暖启动文档）> `_loop_min_f` 绕数
积分（8.8 s cum）> 真根求解（10,122 次，不可约）。
