# SGBZ pairwise 求交重构 — EventGroup + 独立 Mu2Mid

日期: 2026-08-16

## 目标

按讨论确定的方案，取消旧的 μ₂_mid 构建（墙 / sort-change / cubic Hermite
bracketing / `_Mu2MidPath`），改为：

1. `ZeroManager.run()` 后直接用 ItemView 代表元两两求交；
2. `d[:-1]==0` 精确触碰不精化，`d[:-1]*d[1:]<0` 线性预测 + brentq 精化；
3. 零 θ 去重；相近 θ（`< crossing_tol`）合并为一个 `EventGroup` mesh 行；
4. 拓扑荷统一用 side-change 规则，不做单事件/合并事件区分；
5. LineSubset 用 event/MR/seam 定界，`changes_boundary` 的 event 才终止，
   不改变边界对的 event 可 join；
6. μ₂_mid 保留为**独立于 ZeroManager 的简化 `Mu2Mid` 对象**，分段 Hermite
   插值，±14 界在建对象时切分为新的光滑区间。

## 关键实现

- `brute_force_SGBZ/pairwise.py`
  - `PairEvent`：代表元 item 对 + 展开后的真实列集合；
  - `EventGroup`：事件按 `< crossing_tol` 传递合并；item/column 连通分量
    只由事件图生成，**不做任何模长数值判断**；
  - `collect_pair_events`：touch/cross 扫描，零去重；
  - `_refine_pair_crossing`：brentq 精化；内部点 `solve_at(interp='linear')`；
    2π wrap band 内直接按未归一化 θ 求解再 Hungarian 匹配；
  - `finalize_event_groups`：所有排序读两侧最近 regular 行；charge =
    `(side_right − side_left)/2`；seam 通过根值匹配跨 0/2π 映射；
  - `insert_event_groups`：降序插入，seam group 延迟到最后定 row，全部
    group 在最终 mesh 上重新解析 row（避免插入导致的 index 漂移）。
- `continuation/zero_manager.py`
  - `solve_at` / `insert_solution` 增加 `interp='hermite'|'linear'`，
    默认 hermite 完全兼容 continuation/amoeba。
- `brute_force_SGBZ/mu2mid.py`
  - `ItemView` / `_detect_continuum_clusters_internal` / `logabs_clamped`
    原语义保留；event 行不注入 ItemView；
  - `Mu2MidZM.analyze()`：聚类 → pairwise → 插入 → 重建 ItemView →
    `has_continuum = any(j_lo == j_hi)` → 构建 `Mu2Mid`；
  - `Mu2Mid` / `Mu2MidPiece`：不继承 ZeroManager；每段用
    `hermite_interp_poly`；knot 左右导数在 event/MR 断点允许跳变，
    普通 knot C1；
  - `build_mu2_mid`：对每个 Hermite 段求 `poly(t)=±14` 实根并切分，
    界外子段为 `value=±14, dv=0` 常数 piece。
- `brute_force_SGBZ/crossings.py`
  - `detect_crossings_simple` 只物化 EventGroup：覆盖 M-1/M 的连通分量
    展开真实列，每列一个 `PointSubset`；charge 允许 0；MR echo drop 保留。
- `brute_force_SGBZ/continuum_lines.py`
  - run 边界 = 段首/段尾 + `changes_boundary` event；
  - `_join_runs_across_mrs` 只在真 MR / seam 且端点不是
    boundary-changing event 时 join。
- `brute_force_SGBZ/winding.py`
  - `_loop_winding_quad` / `_loop_min_f` 使用 `zm.mu2_mid`；
  - quad 按 `Mu2Mid.breakpoints` 切分。

## 关键修复记录

- `MIN_DIRECTION_DERIV` 最终取 `1e-12`：`1e-8` 会把近 continuum 的
  6.9e-9 导数差误判为切触，导致 W=0、二分提前 w_zero。
- 插入后的 event row index 不能边插边记：后续较小 θ 的插入会整体平移
  已有行；所有 group row 在最终 mesh 上按精确 θ 重新解析。
- 2π seam 附近 θ（进入 ZeroManager 1e-6 wrap band）不能走 `solve_at`，
  必须按未归一化 θ 直接求解后线性预测匹配。
- 单 segment 的 0/2π 两侧 event 必须在圆距 `< crossing_tol` 时合并成
  一个 EventGroup，否则 seam 双计且电荷不守恒。
- seam 跨侧 item pair / side 映射改用 ZeroManager 文档规定的
  `roots_right[boundary_perm] == roots_left` 规范。
- `seg_mu2_values` 兼容数组曾用高级索引生成 `(N, N)` 形状，已改为
  `rows` 对角索引。
- Haldane-y sweep 在 `E=1.212-0.3672j` 失败：MR snap 后 `compute_tangent`
  返回巨大有限导数，μ₂_mid cubic 在 MR 附近冲出 ±14 后被切出伪 ±14 平台，
  loop winding 在 μ₁≈0 出现 ±1 跳变。修复：`SegmentData.tangents` 改为
  run 时必建，MR cluster 列在 snap 后由 `_mark_mr_tangents_inf` 手动置
  `inf`，Hermite 恢复线性 fallback；该点 W 恢复平滑穿过 0。

## 回归

- `pytest tests/`：145 passed。
- HN2D `[10]` / `[11]`：continuum boundary、2 条 LineSubset、μ₁ 解析值一致。
- gain-loss Haldane 历史反例：
  - E=-1.56 → success, is_gbz, index=(12,0)；
  - E=-1.406 → success, is_gbz, index=(6,0)。


## 当前仍存在的问题

Gain-loss Haldane model 的 y-SGBZ 在 E_ref=1.212000 处仍然报错：
``` bash
RuntimeError: bracket bisection failed to converge after 60 iterations: E_ref=(1.2120000000000002+0j), mu1=0.135328598265, len(subsets)=12, is_continuum=False, winding=-0.1028034060294209
```