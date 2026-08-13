# SGBZ 提取重构为 cluster 视角 — 会议记录

日期：2026-08-06
模块：`demos/demo_zm_gbz.py`（SGBZ 部分）；关联 amoeba Rule 1 修复见文末

## 议题

把 `extract_sgbz_subsets` 的偶然交点检测从"逐轨道扫描 dⱼ = ln|β₂ⱼ| − μ₂ 的符号变化"重构为 **cluster 视角**：先用连续简并检测得到简并 cluster，每个 cluster 取一个代表元，再查看按模长排序的第 M−1、M 位根属于哪个 cluster，偶然交点通过**边界 cluster 序号是否变化**来检测。同时修正此前 review 发现的若干缺陷。

## 结构事实（用户提供，整个重构的依据）

|β₂ⱼ|² − |β₂ₖ|² 在相邻 MR 之间是 θ₁ 的实解析函数。因此：

- 两条轨道的模长若在某个开区间上重合，则**整段重合** —— 简并曲线一旦出现就持续到 MR，cluster 是**段内静态**的，不需要逐行维护成员。
- 同 cluster 成员共享同一条 |β₂|(θ₁) 曲线，一条代表元足以携带模长信息。
- PMGBZ 边界 = 每行上占据排序位置 M−1、M 的两个 cluster。

## 决策与设计

### 1. `_detect_curve_degeneracy` → 返回段静态 cluster 索引

返回 `list[list[int]]`（如 `[[1,2]]`，空表 = 无简并）。确认标准：

1. 边界 gap（位置 M−1/M 的模长差）< `continuum_tol`；
2. 两条轨道的 d ln|β₂|/dθ₁ 一致（|Δ| < `_DV_TOL`）；
3. 连续 ≥ 2 行。

**不要求轨道身份跨行相同**：相等模长之间排序是任意的，同一物理曲线可以更换占据 M−1/M 的轨道；身份判据会在近相等模长时抖动、且 3+ 重曲线没有唯一的"对"。成员 = 与边界模长在容差内的**全部**轨道（覆盖 3+ 重），多个 run 若成员集相同则合并（一条曲线被第三者插入边界会分裂成两个 run）。

### 2. `_build_cluster_partition`：划分 + 代表元

简并 cluster 取前几个 id，其余轨道各成单元素 cluster；每个 cluster 一个代表元（`representatives[cid]`）。

### 3. 逐行边界 cluster 身份

`c_lo[i]` / `c_hi[i]` = 排序位置 M−1 / M 的根所属 cluster id。

### 4. 连续区（取代原 per-track frac 测试）

`c_lo == c_hi` 的行 = 边界对简并（cluster 模长即 level）。连续行段（≥ 2 行）→ 每个成员每段一条 `LineSubset`（`_true_intervals`，非循环——段内不跨 0/2π）。**修复 F1**：原 frac 测试把"曲线简并 mask"行也算进连续占比，导致整段发射、端点离 level 最多 0.2（ln 尺度）；现在只发射 on-level 行段。

### 5. touch（MR 处的连续区端点）

MR 行（`left_mr`/`right_mr ≥ 0`，0/2π seam 排除）上 `c_lo == c_hi` → 直接出点（代表元的值）。**修复 F3**：原 `i_hit >= N_orig-1` 静默丢弃每段最后一行 `'zero'`，导致连续区**起点**处的 touch 全部丢失；现在两端对称 + 显式去重。

### 6. 偶然交点（change）与细化

相邻行之间 `c_lo`/`c_hi` 变化 → 记录新旧两个 cluster 的代表元 `(rep_old, rep_new)`，细化目标 `d_AB = ln|r_old| − ln|r_new|`（两个代表元的模长差，与 μ₂ 无关）。三道防线：

- **严格变号过滤**（`d_L·d_R < 0`）：一个轨道穿过 level 时 c_lo 和 c_hi 同时变（如 [a,b,c,d]→[b,c,a,d]），其中 c_lo 的变化对应 d_AB 并无变号（b、c 未互穿，是 a 穿过它们）——过滤掉；
- **穿越轨道选择**：交点 β₂ 取 d 相对 level 变号的那条轨道（另一条可能全程在 level 上）；
- **收敛后 level 检查**（`|ln|r| − μ₂| ≤ tol`）：拒绝**次边界交换**——M−1↔M−2（或 M↔M+1）的模长交叉发生在 level 之下/之上，不是 PMGBZ 点（Q2 确认）。

内部网格行的精确落点（交点恰好在网格行上）：`gap == 0`（严格相等）且 `c_lo ≠ c_hi` → 直接出点（Q1 确认）。非零小 gap 由变号路径覆盖，不能用 `gap < zero_tol` 判据。

最终按 (θ₁ mod 2π, |β₂|) 去重（一个交点从两个位置看到、一个 touch 从两段看到）。

## 未决事项

**Q3：MR touch 的跨段处理**（方案 A/B/C 待定）。K4 γ₃=γ₂ 模型 6 个 MR 中 4 个是真端点（边界对根全在 MR 的 `cluster_indices` 里），2 个是**接缝**（0.6435 / 5.6397：边界对混合——一条轨道在 MR cluster 中聚合终止、另一条穿过继续）。候选方案：

- A：保持现状（6 个 touch 全保留）；
- B：丢弃接缝 touch（4 个点）；
- C：完整按 amoeba 的 `_join_continuum_across_mrs` / `_is_cluster_endpoint` 语义：per-track 判终止（在 MR cluster 中 → touch + 线终止；不在 → 与下一段接续，线合并）→ 4 点 + 8 条跨 MR 接续线。

倾向 C（最符合物理：strip winding 需要 touch，接缝处只有聚合的那条轨道有 touch）。

## 验证

四个模型的输出与精度（每个点满足 |ln|β₂| − μ₂(θ)| ≤ ~1e-10，touch 精确）：

| 模型 | 点数 | 线数 | 说明 |
|---|---|---|---|
| HN2D E=0 | 1 | 2 | 边界 MR touch（θ=0，去重后 1 点）|
| HN2D E=1 | 2 | 4 | 两个 MR touch（2.0944 / 4.1888）|
| K4 γ₃=γ₂ | 6 | 12 | 全部 6 个连续区边界 touch（旧代码只出 2 个重复点）|
| K4 γ₃=0.1 | 4 | 4（裁剪）| 2 横向交点 + 2 touch；线从整段裁剪为 on-level 行段 |

K=4 模型（NNN 跳跃、γ₃≠γ₂）首次触发了 `'cross'` 细化路径（demo 的 HN2D K=2 模型 μ₂ ≡ γ₂ 恒等，只有 touch，`'cross'` 是死代码——review 发现 F2），细化残差 ≤ 7e-11。合成测试覆盖：3 重简并曲线、两曲线同段、次边界交换拒绝、精确落点。

## 关联：amoeba Rule 1 根值筛查

同一对话中对 amoeba（demo `extract_amoeba_subsets` 与 production `brute_force_amoeba/zm_extract.py`）的 Rule 1 做了修复：θ 吸附判据**原样保留**，在其筛出的候选中增加**根值距离筛查**——交点 β₂ 与吸附带内端点根同值（`< ROOT_TOL = 1e-9`，同一条曲线）才丢弃；同 θ 但不同根（不相连的离散零点）保留。`d == 0` 的 `'zero'` 路径与 θ 吸附**都不改**（d==0 概率极低但逻辑上存在，不能因测试数据没有而删）。生产代码测试：无新增失败（13 个失败全部预先存在：test_sgbz 7 个 + test_gbz_types 旧 LineSubset API 6 个）。
