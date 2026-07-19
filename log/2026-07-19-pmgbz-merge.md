# pmgbz_detector 合并入口 + gap 窗口预筛 + 投影复用

日期：2026-07-19（同日第二轮，接 [2026-07-19-pmgbz-dedup.md](2026-07-19-pmgbz-dedup.md)）

## 背景

上一轮去重后，`_full` / `_sweep` 两个入口的编排结构几乎一致，拆分不再有存在价值。本轮在干净的 `_PmgbzScan` 结构上把入口合并回单函数 `get_roots_and_PMGBZ(..., refine_continuum=True)`——与最初"一个函数里两套算法互相拖累"不同，现在的差异只剩 sweep 模式的 continuum 提前返回一个分支。CLAUDE.md 与 doc/SGBZ.md §3.2 原本就按单函数 + `refine_continuum` 撰写，合并后文档重新变准确。

## 改动

### 1. 合并入口

`get_roots_and_PMGBZ_full` / `get_roots_and_PMGBZ_sweep` → `get_roots_and_PMGBZ(poly, E_ref, mu1, N_points, zero_tol, GBZ_check_tol, refine_continuum=True)`。`refine_continuum=False` 即原 sweep：发现 ≥2 连续退化网格点立即返回 placeholder。同步更新 strip_winding_number.py（删除 if/else 分发）、SGBZ.py（3 处调用）、`__init__.py`、tests。

### 2. Stage 3 wide-gap 窗口预筛（P0.1 的健全版本）

对每对相邻非退化网格点：若相对 gap `(|β_M|-|β_{M-1}|)/scale` 在该对及外侧各一邻点（连续 4 点）都 > `gap_skip_tol=0.5`，跳过全部匹配分析。物理依据：交叉要求 gap 归零，连续函数无法在 4 个大 gap 采样点之间藏一个归零点（除非斜率 > 2·tol/h ≈ 48 @ N=301）。阈值故意取大（50% 根模长尺度）。

- 取代了 sweep 的 Hungarian-only 预筛：不再有"不可信但无交叉的对被跳过"的轻量化取舍，未被 gap 跳过的对一律走完整 `process_interval`（更安全）。
- 与被否决的"模长排序比较"版 P0.1 的区别：后者在两端模长区间重叠时漏真实交叉（反例 0.9/1.0→1.05/0.95）；gap 窗口版在该反例中 gap≈0.1 < 0.5，不会跳过。
- 实测跳过率（HN2D，N=301）：mu1=-1（bracket 扩张区，原白烧场景）**301/301 全跳**；mu1=0（GBZ 附近）44%；mu1=0.2（GBZ 上）31%，交叉附近 gap 收窄、永不误跳。

### 3. P0.2 投影 / cost matrix 复用

`_analyze_boundary_matching` 改为 `_PmgbzScan.analyze_boundary_matching` 方法：每侧 `to_sphere_r3` 投影一次、全 cost 矩阵一次计算；全局 Hungarian 直接用它，局部边界块为其切片（已核实逐位一致），double-root self-cost 由缓存投影的切片导出，簇边界也只算一次。gbz_types 拆出 `cost_from_sphere_r3`（`chordal_cost_matrix` 行为不变，amoeba 不受影响）。

### 4. 顺带统一

单点退化 run 的边界二分从 sweep 的 30 次统一为 60 次（消除无谓的模式差异）。

## 行为变化（相对上一轮）

1. sweep 路径：单点退化 run 的 θ 精修更收敛（30→60 次二分，θ 变化 ~2.5e-11）。
2. sweep Stage 3：不可信但无交叉的网格对现在会被细分探索（原 Hungarian-only 预筛会跳过）——更安全，基线用例中无输出差异。
3. 两种模式 Stage 3 都受 gap 窗口预筛门控（启发式；阈值 0.5 + 4 点窗口，见上）。

## 验证

- 基线比对（对合并前 16 份快照）：14/16 逐位一致；2 处差异全部来自变化 1（θ 差 2.49e-11 ≈ 步长/2³¹），无任何点被预筛漏掉。
- `pytest tests/` — 63/63 通过。
- 耗时（collect_GBZ_subsets，N=301）：point 12.46→12.40s（持平；热点已转移到 np.roots/精修/quad，见 TODO P1/P3），continuum 1.99→**1.49s（-25%）**，谱外 0.74→**0.43s（-42%）**。

## 文件改动

| 文件 | 改动 |
|------|------|
| `pmgbz_detector.py` | 单入口 806 行；`analyze_boundary_matching` 入类 + 投影复用；gap 预筛；统一 60 次边界二分 |
| `gbz_types.py` | 拆出 `cost_from_sphere_r3` |
| `strip_winding_number.py` / `SGBZ.py` / `__init__.py` / `tests/test_sgbz.py` | 引用改为单函数 |
| `TODO/SGBZ-optimization.md` | P0 三项全部关闭；基准表补合并后耗时 |
| `CLAUDE.md` | pmgbz_detector 行数更新 |
