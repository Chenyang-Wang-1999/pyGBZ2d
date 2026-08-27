# 常数管理重构：归属地原则 + dict 通道处死

日期: 2026-08-26（继 bfgbz2d 打包四阶段之后）

## 背景

Phase 3 的集中式 `config.py` 在代码审阅中被否决，三个问题：
平铺无层级、20 处签名魔数仍存（含 mu2mid.analyze 的 refine_* 四字面量
与 config.REFINE_* 同值失联）、`zm_run_kwargs`/`solver_options` 等
dict 通道构成绕开常量管理的第四条配置路径（拼错键名静默无效）。

用户拍板的三条设计原则：
1. 先做关联性分析——同一逻辑因不同时期所写而分裂的常量要合并；
2. 单一消费者的常量留在所属文件（如 RK45 步长族只有伪弧长步进器
   用），不进中央配置；
3. 风格统一——要么所有参数都能 kwargs 传入，要么不搞 options，
   让用户直接改常量。

## 关联性分析结论（实证）

**真·同逻辑分裂（5 组，全部合并）**：
- 逃逸梯子 `(1,2,4,8)` 裸写 3 处 → `core.ESCAPE_LADDER`
- `min_dist_threshold=0.1` 双定义（trigger + integrate_segment）→
  `multiple_roots.MIN_DIST_THRESHOLD`
- SGBZ `tol_normalized=1e-2` vs amoeba `PLATEAU_CLUSTER_TOL=1e-2`（同谓词
  同值）→ `core.PLATEAU_CLUSTER_TOL`
- SGBZ `zero_tol` vs amoeba `winding_tol_floor`（同为 |winding|≈0 判定，
  同 1e-10）→ `core.WINDING_ZERO_TOL`
- `detect_cluster` 双默认（run 侧 1e-4 / 直调 1e-6）→ 统一 1e-4
  （用户指定用 run 侧值；**本次唯一预期行为变化**）

**同模式不同量纲（不合并，文档化为族）**：eps 级守卫族、1e-10 收敛
容差族、`min_step` vs `min_dtheta`（步进器地板 vs MR 触发灵敏度，
相邻双旋钮）。

## 最终架构

- 常量住在唯一消费模块（core 7 个跨包 + 11 个模块各自落位），
  `config.py` 删除；`live_defaults` 迁至 core，键格式 `"模块:常量"`；
- 定制模型两句话：**调一次 → kwarg；一直调 → 改模块常量**
  （赋值即时全进程生效，下一次读取即见）；
- 防线固化为 AST lint 测试（test_constants.py）：签名数值默认白名单
  外报错；跨模块常量禁止 by-value import（副本陷阱）；
- dict 通道全处死：`solver_options`/`zm_run_kwargs`/`eval_kwargs`/
  `**options`，`collect_GBZ_subsets` 改显式具名 kwargs（拼错即
  TypeError）；sgbz collect 路径 `continuum_perturb` 1e-2→1e-4 与
  solve 侧统一（dict 时代的残留漂移）。

## 新增文档

- `doc/constants.md`：英文、分模块的完整常数参考（含默认值、语义、
  ⚠ 机器精度锚定标记、同模式族说明、历史合并对照表）。

## 验证

- 快套 260 passed（每步护栏）；
- 最终快+慢 × 双后端（见提交信息）。

## 遗留

- SGBZ collect 路径 perturbation 1e-4 化经全套验证无回归，但宽退化带
  模型上 8e-4 梯子顶若不够用，可 per-call `continuum_perturb` 放宽；
- `TODO/exact-degeneracy-boundary-behavior.md` 仍开放（与本次无关）。
