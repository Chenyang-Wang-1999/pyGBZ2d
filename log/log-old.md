# 开发日志整理

更新时间：2026-05-02

## 概览

这份日志整理了本轮围绕 `phcpy-free-algorithm/brute_force_solver` 的主要讨论、需求演进、代码修改和测试结果。核心主题集中在以下几类问题：

- `get_roots_and_PMGBZ()` 中 PMGBZ 点的识别、排序和漏解问题
- 匈牙利匹配在 `0 / infinity` 分支附近的误判问题
- `WindingFun` 与多项式求导职责混杂的问题
- 废弃接口和旧兼容层清理
- 新的统一递归逻辑：基于 `M-1 / M` 边界等模簇的交换代价差来判断匹配可信度

## 需求演进时间线

### 1. 早期问题修复

- 对 `brute_force_solver` 做 code review，重点检查已经废弃的接口
- 修复单点简并或周期边界上的“同一点”被误判为 continuum 的问题
- 用户明确要求：废弃接口直接删除，不保留兼容
- 修复 `PMGBZ_pairs` 的有序语义，后来进一步升级为按拓扑荷分类输出

### 2. 匹配与误报问题

- 在 `debug_SGBZ` 中发现：`sols_arr[:, 1]` 与 `sols_arr[:, 2]` 无交点，但程序误报了一个交点
- 排查发现：匈牙利匹配在复平面距离下可能出现 `0 -> infinity` 的伪交换
- 解决方案改为：将根投影到复球面，再嵌入 `R^3` 用欧氏距离做代价
- 用户要求不加“硬校验兜底”，而是直接把匈牙利匹配的距离定义改正确

### 3. 绕数与导数结构重构

- 分析 `WindingFun` 依赖，发现其同时负责：
  - 多项式值和偏导计算
  - 绕数积分参数化
- 用户建议拆分职责
- 最终引入 `PolyDiffContext`，负责：
  - `char_poly`
  - 偏导缓存
  - `eval_val()`
  - `eval_partials()`
  - `eval_dmu2()`
- `WindingFun` 保留为绕数计算壳层
- `get_strip_winding()` 与 `get_roots_and_PMGBZ()` 入参从 `WindingFun` 改为 `PolyDiffContext`
- 兼容层后续被删除，主链条统一切换到 `PolyDiffContext`

### 4. PMGBZ 输出结构升级

- 原先输出为 `beta2_pairs`
- 后来升级为：

```python
beta2_sols = [pos_tuple, neg_tuple]
```

- 分类规则基于 `M-1 / M` 边界等模簇中各根的 `dmu2 / dtheta1`
- 其中：
  - 按 `dmu2 / dtheta1` 从大到小排序
  - 前 `M - i_min` 个记为正荷
  - 其余记为负荷

### 5. 漏解问题讨论

- 用户指出：若两组 PMGBZ 点的 `theta1` 很近，可能出现找到一个却漏掉另一个的问题
- 最早的方案是：给 PMGBZ 点加“局部完整性验证器”
- 后续讨论逐步收敛到更简单的思路：
  - 不需要单独一套“局部完整性验证递归”
  - 局部完整性检查应当融入统一的区间递归流程

### 6. 匹配可信度判据讨论

- 讨论过是否需要“重根模式”
- 最终用户明确：
  - 不需要精确重根求解
  - 漏掉部分重根对总绕数影响不大
  - 真正关键的是：每次做匈牙利匹配时，都要判断这次匹配是否可信
- 又进一步明确：
  - “可信”不是区间的稳定属性
  - 区间变小之后，原先可信的匹配也可能变得不可信
  - 因此“可信”不表示进入某种固定模式，只表示“这个区间当前可以做一次 PMGBZ 搜索”

### 7. 最终采纳的实现方向

- 局部完整性验证不再另写第二套递归
- 统一用一个区间处理器 `process_interval(...)`
- 匹配可信度只看 `|beta_M|` 对应的 `M-1 / M` 边界等模簇
- 判据采用该局部簇中的“交换代价差”
- 其他根即使接近重根，也不作为当前主判据

## 本轮已完成的主要代码修改

### 1. `brute_force_solver/strip_winding_number.py`

这是本轮修改最多、也是当前核心逻辑所在的文件。

已完成的关键改动：

- 复球面匹配代价被抽成独立函数 `_chordal_cost_matrix()`
- 新增 `_get_boundary_cluster_bounds()`
  - 用于围绕 `M-1 / M` 边界提取等模簇
- 新增 `_analyze_boundary_matching()`
  - 只在边界簇局部代价矩阵上分析匹配
  - 计算局部最优 assignment
  - 比较“当前匹配”和“发生边界交换时的替代匹配”的代价差
  - 得到当前区间的匹配可信度
- 原先的 `search_accidental_points_in_interval(...)` 被统一为 `process_interval(...)`
- 新递归语义如下：
  - 每次处理区间，先分析边界簇匹配是否可信
  - 若不可信，则直接细分区间
  - 若可信但无跨边界交换，则当前区间结束
  - 若可信且有跨边界交换，则执行一次 PMGBZ 搜索
  - 找到点之后，左右子区间继续递归调用同一个 `process_interval(...)`
- 这使得原先“局部完整性验证”自然融入统一递归中，而不是额外维护一套递归分支

当前 `match_confidence_tol` 的含义：

- 不是绝对误差
- 是相对于局部匹配代价平均尺度的相对裕量阈值

当前实现形式：

```python
confidence = min_exchange_margin / cost_scale
is_confident = confidence > match_confidence_tol
```

### 2. `brute_force_solver/winding.py`

已完成的主要结构重构：

- 引入 `PolyDiffContext`
- 将 `eval_dmu2()` 从 `WindingFun` 迁移到 `PolyDiffContext`
- `WindingFun` 仅保留绕数计算职责

### 3. `brute_force_solver/SGBZ.py`

- 主调用链切换为持有 `PolyDiffContext`
- 统一通过 `get_strip_winding(self.poly_diff, ...)` 进入 strip winding 计算

### 4. `brute_force_solver/root_solver.py`

已完成的清理与修复：

- 删除旧事件接口和废弃匹配接口
- 修复 `complex_root()` 中 Jacobian 仍按初始值而不是收敛值回算的问题

### 5. `brute_force_solver/__init__.py`

- 由通配导出改为白名单导出
- 废弃接口不再暴露
- 新增 `PolyDiffContext` 导出

### 6. 调试和调用脚本同步

已同步调整以下脚本以适配新接口：

- `Haldane-model-gainloss.py`
- `debug-brute-force-solver.py`
- `test-SGBZ-and-Amoeba.py`
- `test-SGBZ-old-algorithm.py`
- `test-winding.py`
- `test-SGBZ-derivative.py`

## 已解决的问题

### 1. continuum 误判

问题：

- 只有一个简并点
- 或者起点与终点实际上是同一个周期点
- 却被误识别为 continuum interval

修复：

- 增加区间宽度判断
- 零宽或近零宽 continuum 退化为 accidental PMGBZ point

### 2. 匈牙利匹配误报交点

问题：

- `debug_SGBZ` 中出现不存在的交点预测

原因：

- 复平面距离会诱发 `0 / infinity` 分支的伪交换

修复：

- 使用复球面弦距作为匈牙利匹配代价

### 3. 求导职责混杂

问题：

- `WindingFun` 同时承担导数上下文和绕数计算

修复：

- 拆分为 `PolyDiffContext + WindingFun`

### 4. 旧接口残留

问题：

- 已废弃接口仍可导入、仍暴露在外部

修复：

- 直接删除实现与导出，不保留兼容层

## 当前仍需注意的点

### 1. `beta2_sols` 仍是两组结构

当前仍为：

```python
beta2_sols = [pos_tuple, neg_tuple]
```

尚未扩展到包含零荷组的：

```python
beta2_sols = [pos_tuple, neg_tuple, zero_tuple]
```

### 2. 重根处理仍是轻量思路

当前采纳的是用户明确要求的简化版本：

- 不单独进入“重根模式”
- 不追求精确重根定位
- 主流程以匹配可信度和 PMGBZ 搜索为主

### 3. 匹配可信度阈值仍需数值经验校准

当前在 `strip_winding_number.py` 中使用：

```python
match_confidence_tol = 1e-3
```

这是当前版本的经验阈值，后续可能还要结合实际模型扫描结果做微调。

## 本轮测试记录

### 1. 语法和诊断检查

- 对 `brute_force_solver/strip_winding_number.py` 执行了 `GetDiagnostics`
- 结果：无新增 diagnostics

- 运行：

```bash
python -m py_compile /home/wangchenyang/654/research-data/2D_skin_effect/phcpy-free-algorithm/brute_force_solver/strip_winding_number.py
```

- 结果：通过

### 2. 轻量调试入口

运行：

```bash
python /home/wangchenyang/654/research-data/2D_skin_effect/phcpy-free-algorithm/debug-brute-force-solver.py
```

结果：

- 退出码 `0`
- `get_roots_and_PMGBZ()` 正常返回
- 该入口下打印出 `5` 个 accidental PMGBZ 点

### 3. `debug_SGBZ()` 实场景

运行：

```bash
MPLBACKEND=Agg python -c "..."
```

调用：

- 动态加载 `Haldane-model-gainloss.py`
- 将 `plt.show` 改成空函数
- 直接执行 `debug_SGBZ()`

结果：

- 退出码 `0`
- `check_SGBZ(...)` 返回 `success=True, is_PMGBZ=True`
- 一次输出中得到 `6` 个 PMGBZ 点
- 随后 `debug_SGBZ()` 内部另一次 `get_roots_and_PMGBZ()` 调用打印出 `5` 个 PMGBZ 点

运行时还出现了两个外部告警，但与本轮修改无直接关系：

- Matplotlib cache 目录不可写
- `BerryPy` 内部触发的 SymPy deprecation warning

## 本轮最终结论

- 新的统一递归逻辑已经落地
- “局部完整性验证”已经并入统一的区间处理器
- 匹配可信度判据已经从“全局想法”变成实际代码
- 至少在现有两个直接调试入口上没有出现语法错误或运行崩溃

## 建议的下一步

如果后续继续推进，最值得优先做的是以下两项之一：

1. 继续跑更多真实参数例子，校准 `match_confidence_tol`
2. 在 `strip_winding_number.py` 中临时加入调试输出，记录每个区间的：
   - `confidence`
   - 是否细分
   - 是否检测到 crossing
   - 最终是否找到 PMGBZ 点

这样会更容易判断当前判据是否与物理直觉一致。

---

# brute_force_amoeba 开发日志

更新时间：2026-05-20

## 概览

本轮围绕 `phcpy-free-algorithm/brute_force_amoeba` 实现了基于 amoeba formulation 的能谱计算模块，与已有的 SGBZ solver 互补。

## 需求与实现时间线

### 1. 基础 pipeline

- 实现 `get_hungarian_sorted_roots()`：在 θ₁ 网格上求解 β₂ 根，用 Hungarian 匹配（复球面弦距代价）建立连续 root tracks
- 实现 `get_a2_average_winding()`：检测 ln|β₂| = μ₂ 的 crossing，分区计算 u₂ 的加权平均
- 复用了 `brute_force_solver` 的 `calculate_point_roots` 和 `PolyDiffContext`

### 2. Crossing 精化

- 线性插值得到 θ₁ 近似值后，用 `scipy.optimize.fsolve` + 解析 Jacobian 精化
- Jacobian：∂f/∂θⱼ = i·βⱼ·∂f/∂βⱼ（由 `PolyDiffContext.eval_partials` 提供）
- 收敛判据：残差 |f| < 1e-10

### 3. Continuum 检测与过滤

- 问题：根在 |β₂| = exp(μ₂) 边界上时，数值噪声产生大量伪 crossing
- 方案：检测 |ln|β₂| - μ₂| < continuum_tol 的连续 θ₁ 区间（`_find_cyclic_true_intervals`）
- 伪 crossing 过滤：若 crossing 两侧点都在 continuum band 内 → 跳过
- 与 SGBZ 的区别：amoeba 只需判断 ln|β₂| 是否在 μ₂ 附近，不需要比较 |β_M| 和 |β_{M+1}|

### 4. a1 方向绕数

- 实现 `_get_average_winding_from_zeros()`：统一处理 a1/a2 方向
  - a2：θ₁ 分区，解 β₂，统计算子
  - a1：θ₂ 分区，解 β₁，统计算子
- a1 复用 a2 的 (θ₁, θ₂) 零点（不需独立的 Hungarian 匹配）
- `get_a1_average_winding()` 薄封装

### 5. 二分求解

- **内层** `_find_mu2_for_a2_zero()`：二分 μ₂ 使 a2 = 0
  - 自适应区间扩展（winding 单调 → 扩展必达异号）
  - Continuum 处理：μ₂ ± ε 微扰
- **外层** `bisect_amoeba_ronkin_min()`：二分 μ₁ 使 a1 = 0（内层保证 a2 = 0）
  - a1 degenerate 时 μ₁ ± ε 微扰，重做内层二分
  - 异号 → 边界停止；同号 → 指导外层二分方向

### 6. 文档与演示

- `demos/demo-amoeba.py`：完整 pipeline 演示
- `doc/amoeba.md`：理论背景、算法、API 参考
- `doc/solver.md`：SGBZ solver 文档（参考 log.md 和源码）
- 删除 `brute_force_amoeba/readme.md`（迁移到 doc/amoeba.md）

## 关键设计决策

### Hungarian 匹配 vs 模长排序

SGBZ 按模长排序（需要复杂的 boundary matching 分析），amoeba 用 Hungarian 匹配建立连续 track。后者更简单，且自然支持 continuum 检测（只需检查单根 track 而非比较两根）。

### 统一绕数函数

`_get_average_winding_from_zeros` 同时服务 a1 和 a2 方向，避免代码重复。a2 有独立管线（Hungarian matching），a1 复用 a2 零点。

### 函数的返回格式

最终返回值中移除 `winding` 字段：要么是 0，要么是一对正负值，不提供有效信息。返回值聚焦于 `(mu1, mu2, zeros, is_continuum, success)`。

## 测试记录

- Square 模型：a1 ≡ a2（对角线对称），Ronkin 极小值在 (0, 0) ✓
- NNN 模型（g=0.3）：Ronkin 极小值在 (0, 0)，continuum=True ✓
- 自适应区间扩展：窄区间 [1.0, 1.5] → 扩展到异号 → 正确找到 mu2 ≈ 0 ✓
- E=5（带外）：a1 无符号变化 → 无中心孔 ✓
- Crossing 精化残差：~1e-15
- a1/a2 收敛性：非退化情况下两种方法随 N 增大趋于一致（diff ~ 1e-3 @ N=801）
