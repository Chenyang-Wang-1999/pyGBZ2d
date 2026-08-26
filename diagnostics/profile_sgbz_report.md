# SGBZ 性能 profile 报告(2026-08-18,只诊断不修改)

工具:`cProfile` + `pstats`;脚本 `diagnostics/profile_sgbz.py`,原始 dump
`diagnostics/profile_sgbz.prof`(可用 `snakeviz` 交互查看)。

## 工作负载

| 负载 | 多项式 | 结果 | 耗时(未 profile) |
|---|---|---|---|
| HN2D 链对 E=1.0(谱内 continuum) | 5 项, K=2 | index=(0,2) | ~40 s |
| HN2D E=5.0(谱外) | 5 项, K=2 | index=(0,0) | ~0.8 s |
| Haldane 增益损耗超胞 E=1.212 | 64 项, K=8 | index=(0,6) | ~39 s |

profile 总时长 132.5 s,其中 **~25.5 s 是 Haldane demo 模块顶层 `import
matplotlib.pyplot` 的字体缓存噪声**,与求解器无关;**求解器净耗时 105.3 s**。
下文百分比均以 105.3 s 为基数。

## 顶层结构:μ₁ 二分 × 冷重建

`solve_SGBZ_for_E` → 每次 `winding_at(μ₁)` 都新建一个 `Mu2MidZM` 并完整
`run()+analyze()`。三个负载共 **40 次完整重建**:

| 负载 | ZM 重建次数 | 平均网格行数/次 | 每次耗时 |
|---|---|---|---|
| HN2D E=1.0 | **33** | 202 | ~1.2 s |
| HN2D E=5.0 | 2 | 66 | ~0.4 s |
| Haldane E=1.212 | 5 | 175 | ~7.8 s |

- HN2D E=1.0 的 μ₁ 探针序列 −1, 1, 0, 0.5, 0.25, … 一路二分收敛到
  continuum 边界 μ₁=0.2,直到某次探针落进 tie_tol 宽度的简并带才经
  `handle_continuum` 退出。**K=2 的平凡多项式花了 40 s,几乎全部是
  二分迭代次数 × 冷重建的乘积**。
- Haldane 只有 5 次重建,但每次 ~7.8 s(K=8:两两配对筛选平方级、
  8 次 degree-8 `np.roots`/步、段更多)——单次构建成本主导。

## 阶段级占比(40 次重建合计,cProfile 累计时间)

| 阶段 | 累计 | 占比 | 位置 |
|---|---|---|---|
| `refine_mesh_for_multiple_crossings` | 40.0 s | **38%** | `pairwise.py:425` |
| `ZeroManager.run`(积分) | 33.6 s | **32%** | `zero_manager.py:409` |
| `build_mu2_mid`(±14 clamp 切分) | 13.4 s | **13%** | `mu2mid.py:669` → `_split_piece_at_clamp:630` |
| `detect_crossings_and_winding`(绕数积分) | ~13.6 s | ~13% | `winding.py`(`_loop_min_f` 9.1 s cum) |

`analyze()` 合计 58.0 s(refine 40 + build_mu2_mid 13.4 + 配对/分组/收尾 ~4.6)。

## 内核级热点

### 1. `np.roots` 打在小多项式上 —— 30.1 s(29%),83,478 次调用

| 调用源 | 次数 | 说明 |
|---|---|---|
| `_hermite_dips_near_zero`(`pairwise.py:187`) | 55,896 | 求三次 Hermite 导数(=**二次式**)的根,用 `np.roots`→伴随矩阵→LAPACK `eigvals` |
| `_split_piece_at_clamp`(`mu2mid.py:630`) | ~17,460 | 每段 2 次:三次式 −(±14),仍是**三次式** |
| `solve_roots_1d`(真实 β₂ 求根) | 10,122 | arclength 8,207 + `_solve` 1,502 + MR 413 |

即 **88% 的 np.roots 调用在次数 ≤3 的多项式上**,每次付出伴随矩阵构造 +
特征值分解的固定开销(相对闭式解 ~两个数量级)。

### 2. Hungarian 匹配 —— 12.1 s(11%),17,783 次

- 每个 arclength 步 **2 次**:`estimate_error`(`arclength.py:324`)+
  接受后的 track 重排(`arclength.py:342`)。
- 每次经 `chordal_cost_matrix` → `to_sphere_r3` ×2(全程 52,063 次,
  tottime 8.5 s —— 小数组下纯 per-call 开销)+ `linear_sum_assignment`。
  K=2 时本质只是 2 个排列的取优。

### 3. `compute_tangent` —— ~8 s,17,284 次

- 逐根 Python 循环,每根 3 次 `eval_partials`,每次新建 `CScalarVec`:
  全程 290,724 次 poly_tools `eval`(8.2 s)+ 204,482 次 `CScalarVec.__init__`(3.0 s)。
- 其中 7,668 次来自 `_compute_tangents`(`zero_manager.py:1120`)——
  `run()` 收尾时对**整个最终网格重算切向量**,而 `arclength_step` 在每个
  被接受行上其实已经算过 V(注释声称 negligible,实测 ~3.2 s)。

### 4. 多交叉筛本身即使零插入也全量扫描

Haldane μ₁=0 单次构建分相计时(总 3.35 s):

| 相位 | 耗时 | 占比 |
|---|---|---|
| `run`(积分) | 1.87 s | 56% |
| `refine_multi_crossings`(**插入 0 行**) | 1.08 s | 32% |
| `build_mu2_mid` | 0.27 s | 8% |
| 其余 | 0.13 s | 4% |

`_find_refinement_plans` 是 O(段 × 网格区间 × item 对) 的全量 Hermitie
筛选,最多 3 轮,与是否找到需要细化的区间无关;二分的相邻 μ₁ 探针之间
也无任何复用。

## 瓶颈清单(按预期收益排序,待批准后再动)

1. **小多项式求根换闭式解**(`_hermite_dips_near_zero` 的二次式、
   `_split_piece_at_clamp` 的三次式):纯局部替换,不动算法语义,
   预计消掉 ~25 s(≈24%)。
2. **μ₁ 二分的暖启动 / 早期 continuum 带检测**:HN2D 场景 33 次冷重建
   是最大乘数;复用上一个 μ₁ 的网格作为初始网格、或用绕数/边界信息预判
   continuum 带,可把重建次数从 ~33 降到 ~10 以下。涉及算法语义,需讨论。
3. **refine 筛选增量化/提前退出**:零插入时仍全量扫 3 轮;至少可按段
   缓存"无可疑区间"的判定、跳过已排除的区间。
4. **arclength 步内 Hungarian 去重**:`estimate_error` 与接受重排共享同
   一次匹配;K=2 时退化为直接比较。~10 s 里可省一半以上。
5. **`_compute_tangents` 复用积分期切向量**(或只在插入行上增量补算)。
6. **`compute_tangent` 向量化 / 批量 eval**:消 CScalarVec 重建与逐根循环。

---

## 附:2026-08-18 已实施的优化(方向 1 的幅度守卫 + 方向 3 的段级增量)

经批准实施了 dip 检验前置幅度守卫与 refine 段级增量扫描;方向 2(小多项式
闭式解)与方向 4(跨 μ₁ 暖启动,见
`TODO/warm-start-zeromanager-across-mu1-probes.md`)未动。

改动:
- `pairwise.py::_find_refinement_plans` — 同号对在数值 dip 检验前先算三次
  Hermite 的 Bernstein 控制点 (d0, d0+h·m0/3, d1−h·m1/3, d1);凸包性质保证
  控制点全在 ±tie_tol 带之外的同号对不可能触零,直接跳过(零假阴性)。
  非有限斜率(线性退化 poly,同号端点不可能 dip)同样跳过。
- `pairwise.py::refine_mesh_for_multiple_crossings` — 首轮复用 analyze 刚
  建好且与网格同步的 clusters/ItemView(消除重复重建);后续轮只重建、
  重扫上一轮实际收到插入行的段;未重扫段的 plans 原样保留(其网格未变,
  重扫结果恒等)。
- `mu2mid.py` — 提取 `_detect_continuum_clusters_for_segment` /
  `_build_item_view_for_segment` 单段内核供上述复用。

实测(同负载,墙钟):

| 负载 | 优化前 | 优化后 | 变化 |
|---|---|---|---|
| HN2D E=1.0 | 39.6 s | ~35.3 s | −11% |
| HN2D E=5.0 | 0.84 s | 0.79 s | −6% |
| Haldane E=1.212 | 39.0 s | 22.0 s | **−44%** |

> **测量噪声警示(2026-08-18 补)**:本机负载波动大,同样代码数分钟内
> HN2D E=1.0 实测 32.8–39.1 s。上表及任何单次墙钟读数只作量级参考;
> 可信的是负载无关计数(np.roots / dip 调用数,见下)与同进程交错 A/B。

profile 复核:`_find_refinement_plans` cumtime 35.2 s → 3.4 s(−90%);
`np.roots` 调用 83,478 → 29,706(−64%,剩余 = 真根 10.1k + clamp 三次式
17.5k + 少量放行的 dip);求解器净时间 105.3 s → 70.5 s(−33%)。
HN2D 收益有限符合预期:K=2 仅 1 个 item 对,dip 检验基数小,其大头仍是
33 次冷重建(即未实施的方向 2/4)。

### 追加:±14 clamp 切分守卫(2026-08-18 第二批,经批准)

`_split_piece_at_clamp`(`mu2mid.py`)对每个 mesh 区间无条件做两次
`np.roots`(±14 界,复系数伴随矩阵)。插桩实测两个负载 8,601 个 piece
**100% 可被凸包守卫跳过,0 次真实切分,0 个 ±14 常值段**——该机制只在
MR 两侧边界根 ln|β₂| 发散时才有意义(|β₂| = e¹⁴ ≈ 1.2×10⁶)。

改动:同源 Bernstein 凸包守卫(控制点全在带内 ⇒ 跳过两次求根;闭不等式
安全——内部切线触界即使漏掉也是输出 no-op,因为切分后按端点值/导数重
拟合 reproduces 同一三次式;非有限斜率的线性退化 hull 取端点对)。提前
返回与原无切分路径逐位一致(同一 `_make_piece` 重拟合)。

交错 A/B(同进程交替执行守卫版/原始版,消除机器负载漂移):

| 负载 | guard | orig | 加速比 |
|---|---|---|---|
| HN2D E=1.0 | 32.8 s | 39.1 s | **1.19×** |
| Haldane E=1.212 | 17.8 s | 18.9 s | 1.06× |

负载无关计数(profile):`np.roots` 29,706 → **12,246** 次(cum 9.5 →
4.4 s),剩余 = 真根 10,122(不可约)+ dip 二次式 1,837 + planner 三次式
287,对账闭合;`_split_piece_at_clamp` cum 12.45 → 5.49 s(余量为
polyder/polyval/_make_piece 的正当构建开销);`analyze` 25.4 → 18.4 s。

测试:`test_sgbz` + `test_counterexamples` + `test_regressions`
65 passed + 2 slow-skipped,三负载 index 不变。
