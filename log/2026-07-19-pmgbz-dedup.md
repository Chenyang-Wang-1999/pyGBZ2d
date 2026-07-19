# pmgbz_detector 结构去重

日期：2026-07-19

## 背景

2026-07-18 的拆分（见 [2026-07-18-sgbz-split.md](2026-07-18-sgbz-split.md)）把 `get_roots_and_PMGBZ` 分成 `_full` / `_sweep` 两个入口，解决了 sweep 模式白跑 Stage 3 的问题，但实现是复制粘贴式的：文件 1045 行中 ~300 行逻辑逐字重复（`process_interval` 递归、二分精修、点状退化验证、输出装配各存两份）。本轮做纯结构去重，**行为严格不变**，性能优化留到下一轮。

## 重构方式

新增内部类 `_PmgbzScan` 承载共享机制（网格求根 + cache、连续统边界二分、偶然点二分、统一区间递归 `process_interval`、退化 run 处理、Stage 3 扫描、输出装配）。**编排留在两个公开入口函数里，类不设 run_full/run_sweep**——避免变相回到"单函数 + mode 参数"，full/sweep 的真实差异在入口里直读：

| | `_full` | `_sweep` |
|---|---|---|
| 连续统 run ≥ 2 | `handle_degenerate_interval` 精修边界 → LineSubset | 立即返回 placeholder + `continuum_flag=True` |
| 单点退化 run 边界精修 | `boundary_max_iter=60` | `boundary_max_iter=30` |
| Stage 3 预筛 | 无（`process_interval` 自行决定，含不可信区间的细分探索） | `cross_prefilter=True`（仅 Hungarian 交叉检查，故意跳过"不可信但无交叉"的对） |

其余全部共享。新增模块级小工具 `_normalize_theta`、`_theta_close`（取代两份闭包和 4 处环绕距离比较）。

## 顺带删除的纯冗余

- **full 模式 Stage 3 双重分析**：原先每对网格点先做 `_analyze_boundary_matching` 预筛，通过后 `process_interval` 第一步又以相同参数重算同一分析。预筛结果从未复用，删除后行为逐点一致、full 模式分析调用减半。
- **`_analyze_boundary_matching` 死返回键**：`matches` / `confidence` / `exchange_margin` / `cluster_block` 全库无消费者，从返回 dict 删除（计算逻辑不动，惰性化属下轮性能项）。
- 删除 `_full` docstring 中语义混乱的 FIXME；"Why split" 注释块改述为"共享机制 + 两个编排入口"。

## 重要发现：TODO P0.1 否决

原计划的"模长排序预筛"不可靠：两端模长区间重叠时真实交叉仍可发生（反例：0.9/1.0 → 1.05/0.95），预筛会漏 PMGBZ 点导致绕数错误。已在 [TODO/SGBZ-optimization.md](../TODO/SGBZ-optimization.md) 记录否决理由；P0.3 修订为"double-root 检测惰性化"（结构性修复，无阈值）。

## 验证

- 重构前用临时脚本（/tmp/pmgbz_baseline.py）固化 8 用例 × full/sweep = 16 份快照（覆盖 continuum 提前返回、LineSubset 精修、偶然点、double-root 路径入口、空结果、Haldane 4-band），重构后比对：**16/16 逐位一致**（含 `theta1_arr`、`sols_arr`、`_pmgbz_raw`、subsets、index）。
- `pytest tests/` — 63/63 通过。
- `demos/demo_unified.py` — 输出符合解析预期（SGBZ mu1=0.2、LineSubset、谱外为空）。

## 文件改动

| 文件 | 改动 |
|------|------|
| `pmgbz_detector.py` | 1045 → 830 行；`_PmgbzScan` + 两个编排入口；删双重分析与死返回键 |
| `TODO/SGBZ-optimization.md` | P0.1 否决；P0.3 修订为惰性化；已完成项归档；行号引用改函数名 |

公开 API（函数签名、返回值、`__init__.py` 导出）无任何变化。
