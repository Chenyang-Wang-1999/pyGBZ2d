# brute-force-non-hermitian 全项目 Code Review 报告

日期:2026-08-18
**状态更新(2026-08-18 晚)**:二、三、四、五节所列问题已全部修复(M7 经用户决定保留不修——`_loop_min_f` 是低成本启发式,只需"足够大"无需精确);`demos/demo_zm_gbz.py` 已删除、`mu2_mid_breakpoints` 死代码已清理。回归测试见 `tests/test_regressions.py`(H1/H2/M4/M5)、`tests/test_counterexamples.py`(T4 Haldane 反例/T5 pairwise 深层/T6 probe+MR solver)、`tests/test_gbz_types.py`(T3 重建)。第六、七、八节(文档一致性、工程卫生、疑点)待用户处理。
方法:完整测试套件基线 + 四路并行模块审查(brute_force_SGBZ / brute_force_amoeba / continuation / 测试与文档),所有高严重性发现均经主会话二次实证(复现脚本或代码级数据流推演),非转述。
代码规模:约 10,157 行 Python(源码 + 测试)。
基线:`pytest tests/` 145 项全过(4 分 46 秒)。

---

## 一、总体评价

- **测试基线健康,质量高于科研代码平均水准**:解析值精确断言、历史 bug 回归测试、容差有 docstring 论证(如 test_sgbz 的 2e-3 明确排除 eps~1e-8 非物理区),boundary_perm 3-cycle、合成两段拓扑 join、有限差分交叉验证等均为高质量设计。
- **架构方向正确**:μ₂_mid 重构 + pairwise EventGroup 的设计文档(log/2026-08-16-pairwise-linear-intersection.md)与实现吻合;"predicted→solved 匹配锚"设计贯穿所有匹配点;"关键修复记录"的工程判断(如 MIN_DIRECTION_DERIV=1e-12 的取值依据)质量很高。arclength→multiple_roots→zero_manager 分层清晰。
- **经验证正确的部分**:interpolation.py 三次 Hermite 公式(符号推导 + 数值双重验证)、`_get_average_winding_from_zeros` 与解析值精确一致(非平凡 u∈{0,1} 模型验证)、`_compute_zero_dtheta1_dmu2` Cramer 推导验算无误、`_boundary_perm_from_right` 逆置换推导正确、`find_cyclic_true_intervals` 边界情形正确。
- **但存在 3 个已实证的高严重性 bug**、8 个中严重性问题,且核心源码处于未提交状态。

---

## 二、高严重性问题(3 项,全部二次实证)

### H1. ZeroManager MR 重启无下界钳制 → 死循环,同一 MR 无限重录

位置:`continuation/zero_manager.py:1212`

```python
theta_temp = min(2 * theta1_mr - ref.theta, theta1_mr + mr_jump)
```

镜像公式 `2θ_mr − θ_ref` 在参考点已越过 MR 时落到 MR 之前,`min()` 又偏向更左侧候选。重启点距 MR 仅 `mr_jump`(默认 1e-6)时,平方根分支点上首步 `dθ ≈ 2·h·mr_jump < min_dtheta` → point trigger 立即重新触发 → 再次精化出同一 θ_mr → 无限循环。

**实证(双重复现)**:
- `f = β₂² − (β₁ − i)`(双重根位于 θ₁=π/2、β₂=0),`h0=0.5, min_dtheta=1e-6`:45 秒内 `zm.multiple_roots` 增至 **3572 条完全相同记录**,run() 永不返回;
- 子代理复现:`h0=0.1, min_dtheta=1e-6` → "did not reach 2π within 10000 segments",10000 条相同 θ₁=1.5707963 的 MR 记录。

**触发面**:`min_dtheta ≳ 2·h0·mr_jump`。仓库自己的测试(test_zero_manager.py:153 等)传 `min_dtheta=1e-6`,正处于危险区间边缘,只因其 poly_F 的 MR 恰在 θ=0 边界(由初始 detect_cluster 捕获,不走此重启路径)才未暴露。

**修复建议**:
(a) `theta_temp = max(theta1_mr + mr_jump_eff, min(镜像, θ_mr + mr_jump_eff))`,其中 `mr_jump_eff = max(mr_jump, C·min_dtheta/h0)`(C≈10,由 dθ≈2hΔ 的分支点标度得出);
(b) MR 记录去重(新 θ_mr 与上一 MR 的 θ₁ 差 < ε 时复用索引)作为第二道防线。

### H2. `_join_runs_across_mrs` 的 merge 预算不足 → LineSubset 断裂

位置:`brute_force_SGBZ/continuum_lines.py:171-215`

外层 `for _ in range(n_seg)` 每轮只完成一次 merge(merge 后双重 break 重启),总 merge 数 ≤ n_seg;但所需 merge 数 ≈ 连续体多重数 × 跨越的段边界数,与 n_seg 无关。

**实证**:`diagnostics/review_join_cap_repro.py`——2 段、3-fold continuum,6 pieces 进 **4 pieces 出**(应为 3 条闭合曲线,track C 断裂为两段未闭合 piece)。2 重连续体恰好等于上限,故现有测试全绿未暴露。

**修复建议**:改 `while changed:` 不定轮循环(merge 严格减 piece 数,天然终止),或预算取 `len(pieces)`。

### H3. 左 bracket 扩展的 continuum proxy 使 `subsets=None` 逃出异常处理

位置:`brute_force_SGBZ/sgbz_solver.py:264-291` → 崩溃点 `collect_GBZ_subsets:525`

左端命中 continuum 时 `winding_at` 返回 `(None, None, zm)`,`gbz_low` 保持 None;若 `handle_continuum` 同号分支返回的 proxy 有 `w_l >= 0`(含恰为 0 的零 plateau 情形——`compute_average_winding` 返回精确整数,w=0 可达,代码注释自己承认该情形可达),break 后 `left_endpoint_zero` 分支返回 `"subsets": None`,在 `collect_GBZ_subsets:525` 处 `sum(1 for s in None)` 抛 TypeError,且异常在 try/except **之外**(try 只包住 solve 调用),直接崩溃而非返回 `success=False`。

**证据**:右端点路径有 proxy 后 `continue` 重求值(303-314 行),左端点没有——非对称性确证是遗漏而非设计。

**修复建议**:左端 proxy 分支后重新 `winding_at(mu1_low)` 物化 subsets,与右端对齐。

---

## 三、中严重性问题(8 项)

| # | 位置 | 问题 | 验证 |
|---|------|------|------|
| M1 | `brute_force_amoeba/zm_extract.py:439-443` | Rule 1 seam 去重用线性距离 `np.abs(cont_endpoints - t1)` 而非 `circ_dist`:refined crossing wrap 到 [0,2π) 后与 2π 端点线性距离 ~6.28 > snap_tol,去重失效 → seam 处重复 PointSubset。同文件 `_zero_identity_key`(497-502 行)正确处理了同一 seam 问题,Rule 1 漏配,自洽性缺口。`circ_dist` 在 amoeba 全模块零使用 | 读码 + 数值演示(线性 6.283 vs 圆周 1e-9)确认 |
| M2 | `brute_force_amoeba/amoeba.py:91-114` | plateau 探针 evaluator 硬编码 `"success": True` 且无 try/except;共享协议 `probe_zero_plateau`(gbz_types.py:543-621)设计了 `success=False` 路径但收不到。`bisect.py` 三处 raise(171/191/242)均可穿透——探针恰在能带边缘探测,正是内层 bisection 最易 range-expansion 失败处,单侧 probe 落谱外时必然失败 → 整个能量点变成 `success=False` | 读码确认 |
| M3 | `brute_force_amoeba/zm_extract.py:606-636` | 共享 MR 边界行的 `d==0` 精确触点:相邻两段记同一行,winding 积分有 `np.unique` 保护(ronkin_winding.py:109,数值验证过 dup 与 single 给出相同 winding),dW/dμ₂ 循环(631 行)没有 → Newton 步 `delta = -w_ref/dW_dmu2` 双计偏一倍 | 读码确认;缓解:better-of-two 回退 + d==0 浮点精确相等罕见 |
| M4 | `continuation/multiple_roots.py:287-296` | `_closest_pair_deriv` 不排除奇异根:最近对含 V=nan 哨兵根时 `deriv=nan`,符号比较均 False → approaching→separating 翻转漏检;且 nan 写回 `_prev_deriv` 使后续步骤也永久失明 | 子代理最小复现(4 根数组 + 控制组对比);机制与 compute_tangent nan 哨兵数据流吻合 |
| M5 | `continuation/arclength.py:262` | `estimate_error` 的 `scale = atol + rtol * np.median(...)` 在实际根含 ≥⌈K/2⌉ 个 inf 时 median=inf → error_norm=0,完全错误的步被无条件接受且 h 放大到 max_factor | **本人实证**:`[1,2,3]` vs `[inf,inf,3]` → error_norm = 0.0 |
| M6 | `brute_force_SGBZ/sgbz_solver.py:349-388` | 二分循环中 continuum proxy 的 `f_mid==0` 不走收敛检查(检查只在真实 winding 分支内,368 行),bracket 退化永远走 else 分支 → 收敛到错误端点后必然 RuntimeError;按语义 f_mid==0 的 proxy 点本身就是答案 | 读码确认 |
| M7 | `brute_force_SGBZ/winding.py:169-202` | `_loop_min_f` 种子安全度量用离散 `seg_mu2_values`,而积分实际走连续 Hermite `mu2_mid` 路径(`_loop_winding_quad` 用 `path.value_deriv`)——事件密集处可能选中离零点更近的种子 | 读码确认;影响限于启发式种子质量 |
| M8 | `brute_force_SGBZ/continuum_lines.py:24` | SGBZ 反向 import amoeba 私有符号 `_LinePiece/_is_cluster_endpoint`,分层违规(全项目唯一一处此方向依赖) | grep 确认;建议提升到 gbz_types 或独立模块 |

---

## 四、低严重性 / 代码质量

**死代码与死参数**
- `crossings.py`:`_reanalyze_if_needed`、`_is_near_mr`、`_classify_charge` 无消费者
- `mu2mid.py:424`:`mu2_mid_breakpoints` 永远为空列表,但 `demos/demo_zm_gbz.py:470,583` 仍在读(实证)
- 死参数:`solve_SGBZ_for_E` 的 `xtol`(TODO B1)、`max_newton`、amoeba 的 `N_points`/`min_continuum_pts` 全链(`_refine_and_correct` 形参名 `min_continuum_pts_unused` 自认)
- 未使用 import:`gbz_types.py` 的 `chain/TYPE_CHECKING/Any`、`arclength.py` 的 `pi`、`bisect.py` 的 `CONTINUUM_TOL`、`amoeba.py` 的 `exp`、`zm_extract.py` 的 `_get_average_winding_from_zeros`
- `_merge_two` 的 `cyclic` 形参功能已死(`zm_extract.py:286-301` 两分支逐行相同,实证)
- TODO/dead-code-cleanup.md 部分过时(A1 已修,其余多数仍有效)

**脆弱模式与魔法数**
- `pairwise.py:485` `g in seam_groups` 值相等比较(EventGroup 是 dataclass,应为身份比较 `any(g is s ...)`)——高度对称退化模型下可误跳过,当前非活性 bug
- 魔法数:`bisect.py:56`(1e-15 dW 阈值)、`ronkin_winding.py:75`(1e-10 绝对残差)、`zm_extract.py:189/196`(1e-9 与 ROOT_TOL 重复定义)、`amoeba.py:83`
- `plateau_area_threshold` 一参两义(winding 面积阈值 + torus 聚类距离阈值,量纲不同)
- `continuation/multiple_roots.py:34` NamedTuple 可变默认值 `cluster_stds=[]`(潜伏风险)

**结构与注释**
- `zero_manager.py` 的 `run()` 约 270 行,混合边界初始化/段组装/MR 分派/pending 管理;`_refine_mr` 返回 6 元组,建议 dataclass
- 镜像重启公式 `min(2θ_mr−θ_ref, ...)` 无"为什么取 min"注释——恰是 H1 根源
- 文档失实 ×4:arclength.py:157 "Vⱼ=0" 实为 nan;locate docstring "snap to 0" 无 snap;`insert_solution` 注解 int 实返元组;`_closest_pair_deriv` docstring 与返回不符
- `collect_GBZ_subsets` 的 `perc` 在库入口无条件 print
- 已 analyzed 时 `zm_run_kwargs`/`crossing_tol` 被静默忽略
- 边界 MR 停止分支(zero_manager.py:582-584)留 (θ_mr, 2π) 缝隙,该区间 `locate()` 抛 ValueError(代码级证据,未端到端复现)

---

## 五、测试套件评价

**强项**
- 断言强度高:解析值精确断言为主,容差有论证;boundary_perm 3-cycle、合成两段拓扑 join、有限差分交叉验证设计良好
- 边界覆盖扎实:θ 跨 0/2π、重根、NaN/0/∞ 根、logabs_clamped 的合成 fixture(诚实承认 HN 模型不产生 0/∞ 边界根)
- test_sgbz.py 的 pairwise 章节扎实:单元(FakeZM merge 测试)+ 集成 + TestReviewFixes 回归类

**缺口(按优先级)**
1. **T1** `tests/test_amoeba.py:477` `from BerryPy import ...` 无 `pytest.importorskip` 守卫——无 BerryPy 的机器测试直接 collection 崩溃,与 README "(Optional)" 矛盾。一行修复
2. **T2** H1/H2/H3 均无测试暴露;H1 的危险参数组合就在现有测试传参习惯里
3. **T3** `test_gbz_types.py` 随 a118614 删除后无替代:`find_cyclic_true_intervals`(seam 边界语义核心原语)、`chordal_cost_matrix`、`generate_probe_steps`、`probe_zero_plateau` 共享工具零直接测试
4. **T4** log 里的手工反例(gain-loss Haldane E=-1.56/-1.406)未沉淀为自动化回归,重构时只能靠人肉记忆
5. **T5** pairwise 深层内部(`_refine_pair_crossing` 的 brentq-fallback、`insert_event_groups` 降序插入 + seam 延迟定 row)缺直接单测
6. **T6** SGBZ 侧 `_probe_zero_plateau_near_mu1` 零测试(amoeba 侧有 TestPlateauEdge);`solve_multiple_roots_iterative`(4D Newton)无覆盖
7. **T7** `build_HN2D_polynomial` 在 4 个测试文件近似拷贝(参数还略有差异),应上移 conftest
8. **T8** conftest 的 slow marker 死配置(无测试标记、无 --run-slow 实现);`test_outside_spectrum` 析取弱断言

---

## 六、文档一致性

- **CLAUDE.md 严重过时**:列 `brute_force_SGBZ/continuum.py`(实际 continuum_lines.py)、`brute_force_amoeba/tracks.py`(已删)、`tests/test_gbz_types.py`(已删);漏 `mu2mid.py`、`pairwise.py`、`continuum_lines.py`、`interpolation.py`、`zm_extract.py`、`test_interpolation.py`;称 LineSubset 为 "lazy beta2_arr"(gbz_types.py:213 明确 eager)
- **README**:引用不存在的 `demos/demo_amoeba.py`、`demo_solver.py`;`pip install -e .` 指示必然失败(项目无 pyproject.toml/setup.py)
- **doc/SGBZ.md、amoeba.md、continuation.md 与代码一致性好**(行数标注与 wc -l 吻合,continuation.md 已涵盖新的 interpolation.py)

## 七、工程卫生

- ⚠ **核心源码未提交(最高优先行动)**:`brute_force_SGBZ/pairwise.py`(743 行,测试已依赖)、`continuation/interpolation.py`、`tests/test_interpolation.py`、`TODO/`、`paper/`、6 个 log 文件均 untracked;另有 20 个已跟踪文件 +1339/−1775 行未提交,横跨多个逻辑变更。建议按 log 条目拆 3-4 个提交(μ₂_mid 重构 / pairwise 重构 / interpolation 统一 / paper+TODO)
- **依赖不可复现**:无 requirements.txt/pyproject.toml;poly_tools 是外部 C 扩展(README 有安装说明,尚可)
- **无 CI、无 lint 工具配置**(pyflakes/ruff/mypy 均未安装;本次 review 用 AST 扫描即发现真实死 import)
- `.gitignore` 正确(data、*.pkl 均忽略,git ls-files 确认 0 个 pickle 入库);log/ 维护习惯好

---

## 八、疑点(未确证,值得关注)

1. **切线帧错位(continuation 疑点 1,代码级事实链已确认)**:`_mark_mr_tangents_inf`(zero_manager.py:1099-1121)用 `multiple_roots[mr_idx].cluster_indices`(θ=0 处**模排序帧**索引,441 行在 sort_by_root_abs 之后检测)标注段边界行(**track 帧**)切线为 inf。boundary_perm 非平凡且移动 cluster 成员时(仓库自己的测试证明 perm 可以是 3-cycle)会标错列。`seg.tangents` 被最新的 SGBZ 代码直接消费(pairwise.py:288、mu2mid.py:375)。未端到端复现;建议在 `_append_segment` 加 cluster 列与末行实际重合的断言
2. MR 恰在 β₂=0 时双触发器联合失明(两根同时进 |β₂|<1e-6 窗口 → V 全 nan → point/interval 触发器同时失效)
3. `solve_multiple_roots_iterative` 对 |β₁| 无约束,解可漂离 μ₁ 圆,仅 warning 兜底;漂移超限时可产生 1 行段或破坏段间 θ 单调性
4. `solve_roots_1d` padding 最多各补一个 0/inf,次数亏缺 >1 时返回数组短于 M+N,下游 K 形状校验会崩(gbz_types.py:166-170)
5. SGBZ "W 单调递增于 μ₁" 假设、boundary MR 邻近整组丢弃可能误删合法 crossing、单网格区间双穿越漏检、padding 根进入 ItemView 排序
6. amoeba:`mu2_mid` 配对 `w_left` 的 ε 级偏移、continuum frac=0.9 阈值脆性(borderline 时拓扑不一致直接 raise)、`_is_cluster_endpoint` 值匹配错配风险、Newton clamp 到初始括区间

---

## 九、建议执行顺序

1. 🔴 `git add` untracked 核心源码并分批提交(先于一切——任何误操作都可能丢失大重构成果)
2. 🔴 修 H1(一行钳制 + mr_jump_eff + MR 去重双保险)+ H2(`while changed`)+ H3(左端补重求值),各配回归测试——三个复现脚本/参数组合均已就绪(H1: `f=β₂²−(β₁−i), h0=0.5, min_dtheta=1e-6` + watchdog;H2: `diagnostics/review_join_cap_repro.py`)
3. 🟡 M1-M5 逐项修复(seam 去重换 circ_dist、探针 try/except、dW 去重、nan 过滤、median 有限根)
4. 🟡 测试补强:BerryPy importorskip(T1,一行先做)→ H1/H2/H3 回归(T2)→ 重建 test_gbz_types.py(T3)→ 反例固化(T4)→ pairwise 深层单测(T5)→ 工具去重(T7)
5. 🟢 工程化:pyproject.toml(声明 numpy/scipy/poly_tools)、ruff + CI、CLAUDE.md/README 更新、dead-code 清单执行(TODO 部分项已过时先勾掉)、M8 分层违规修复(`_LinePiece`/`_is_cluster_endpoint` 提升到 gbz_types)
6. 🟢 疑点 1 的 cluster 列重合断言(保护最新的 pairwise 消费链)

## 附:Review 产出文件

- `diagnostics/review_join_cap_repro.py` — H2 证据(修复后可删)
- `diagnostics/code-review-2026-08-18.md` — 本报告
