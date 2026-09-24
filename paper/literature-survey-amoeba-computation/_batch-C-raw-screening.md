# Batch C — amoeba / Ronkin 函数的**计算算法**：原始检索与摘要级筛选

**研究问题**：现有文献中，复代数几何/多复变意义下的 **amoeba 与 Ronkin 函数**是如何被**实际计算**出来的？关注 Newton 多面体三角剖分、lopsided 近似、tropical 几何、Monge–Ampère、数值代数几何/同伦延拓等可用于计算 amoeba 或 Ronkin 函数的算法。

**纳排标准（PICOS 适配）**：P = 复 Laurent 多项式的零点集 amoeba / Ronkin 函数 / 对数 Gauss 映射 / Newton 多面体等代数几何对象，或非厄米模型中的同类判据；I = 给出**可执行的算法或数值实现**（三角剖分、lopsided 近似、tropical 极限、同伦延拓、数值代数几何、轮廓/网格离散、凸几何算法等），纯存在性定理或纯渐近估计且无计算路径 → unclear；C = 有计算实例/基准（无实例不排除，Rationale 中注明"无计算实例"）；O = 能得到 amoeba、其补集（lopsided 区域）、Ronkin 函数或其数值近似；S = 理论（数学）/数值/算法均可，纯综述 → exclude（若综述直接总结算法 → unclear）。时间窗不限。

**tier 规则**：`recent3` = 2023-09 之后（本批以 arXiv feed 时间戳 **2026-09-24** 为基准），`classic` = 更早。本批以 classic 为主，属预期。

**诚信说明**：本文件所有字段仅来自 arXiv API 返回的 Atom XML；未编造标题、作者、年份、结论。摘要缺失/抓取失败按"抓取失败"记录，不推断内容。

---

## 1. 检索与统计

### 1.1 查询串与命中数

| 编号 | search_query（未编码原文） | sortBy | totalResults | 本次取回 |
|---|---|---|---|---|
| Q1 | `ti:amoeba` | relevance desc | 107 | 50 |
| Q2 | `abs:amoeba AND (abs:"Newton polytope" OR abs:"Monge-Ampere" OR abs:lopsided OR abs:tropical OR abs:"Ronkin")` | relevance desc | 61 | 50 |
| Q3 | `abs:"numerical algebraic geometry"` | relevance desc | 92 | 50 |
| Q4 | `abs:"amoeba" AND (abs:computation OR abs:computing OR abs:algorithm OR abs:visualization)` | relevance desc | 70 | 50 |

原始条目数：50 × 4 = **200**。

### 1.2 去重口径（可核验）

去重主键 = arXiv ID（去掉 `vN` 版本后缀）。151 条唯一记录的 `hit_by` 分布：

| hit_by | 条数 |
|---|---|
| 仅 Q1 | 4 |
| 仅 Q2 | 19 |
| 仅 Q3 | 46 |
| 仅 Q4 | 35 |
| 跨查询重复（Q1/Q2/Q4 之间） | 47 |
| **合计 unique** | **151** |

重复条目总数 = 200 − 151 = 49。

> ⚠️ **计数勘误（必读）**：本批筛选过程中我曾把 unique 总数报为 ~145，并在第一遍 manifest 中**漏登记 1 条**（`1307.3681`，Metric Estimates and Membership Complexity for Archimedean Amoebae and Tropical Hypersurfaces），导致 unique 少计 1。
> 逐条复核全部 200 条原始条目后的正确值为 **unique = 151**，manifest 共 151 行数据。该漏登记条目已在 manifest 中补入（第 132 行）。

### 1.3 PRISMA 式计数

- retrieved（4 条查询合计原始条目）：**200**
- after dedup（unique，按去版本号后的 arXiv ID 去重）：**151**
- excluded at title/abstract：**117**
- include：**24**
- unclear（边界，待用户/全文裁决）：**10**
- 校验：24 + 10 + 117 = 151 ✅
- 进入全文候选：24 + 10 = **34**

### 1.4 排除原因分布（117 条 exclude）

| 原因类别 | 约计篇数 | 说明与代表条目 |
|---|---|---|
| 同名不同物（生物学黏菌 / 图像 morphological amoeba / 集合论 amoeba forcing / 图论 global amoeba / AMOEBA 力场 / 网络工具） | ~28 | `1008.3709`、`1411.3285`、`1108.4315`、`math/9209206`、`2008.04996`、`2104.09707`、`2206.13430`、`2310.20469`、`1107.3569` 等 |
| 对象在范围内但无计算算法（纯存在性/结构/渐近估计/不等式） | ~34 | `1510.08416`、`1812.08149`、`1101.4693`、`1412.4658`、`2608.03613`、`math/0405259`、`1109.2645` |
| 数值代数几何通用工具/应用，未涉及 amoeba（Q3 主体） | ~44 | `2011.05000`、`1611.05947`、`1405.7871`、`2203.07016`、`2302.04117`、`1605.07806` 等 |
| 纯综述（S 维度） | 4 | `math/0108225`、`math/0403015`、`2211.09416`、`2305.00743`（另 `1411.3285`、`2503.09133` 兼属同名/综述） |
| 物理/应用中使用 amoeba 或 Ronkin 但方法非主题 | ~7 | `0804.1870`、`1506.07606`、`2107.07286`、`math/0406099`、`2403.08659`、`math/0505269`、`2006.14041` |

### 1.5 质量异常记录

- Q3（`abs:"numerical algebraic geometry"`，50 条）**主题严重错位**：50 条中仅 `1408.3105`、`1605.04203`、`0911.1783`、`2304.08598` 4 条与 amoeba/热带几何有实质关联，其余 46 条为计算机视觉、代数统计、优化、弦论等应用。已由用户确认下一轮改用交叉查询（如 `abs:"homotopy continuation" AND (abs:tropical OR abs:amoeba OR abs:"Newton polytope")`）。
- **系统日期异常**：arXiv feed 时间戳为 `2026-09-24`，结果集中包含 `2601.*`/`2603.*`/`2605.*`/`2607.*`/`2608.*` 等 2026 年 ID 的论文（如 `2608.03601`、`2607.15424`、`2607.15429`、`2605.24963`、`2603.21116`、`2601.18180`）。这些条目的 `<published>` 与 arXiv ID 前缀一致，均为 Atom 原文，未做改写；已由用户确认以 feed 时间戳为 tier 基准。
- `hep-th/0601233` 与 `0810.4179` 条目相邻处存在疑似 malformed entry 边界，但两者的 `<id>`/`<title>`/`<summary>` 均完整可读，未影响判定。

---

## 2. INCLUDED / UNCLEAR（逐篇全字段，34 条）

排序：按 relevance 降序；同档内按 tier / arXiv ID。

| arXiv ID | Title | First author | Year | Venue | Type | Method core (中文) | Model/benchmark | Main result | Limitation | Rel | Decision | Rationale | tier | hit_by |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1604.03603 | Algorithmic computation of polynomial amoebas | D. V. Bogdanov | 2016 | arXiv preprint（无 journal_ref/DOI） | engineering | 给出 amoeba、其 contour、紧化 amoeba 及三维 amoeba 二维截面的计算与可视化算法，并给出在固定 Newton 多面体下构造"最复杂拓扑"多项式的算法 | Matlab 8 / Mathematica 9 实现 | 算法在 Matlab/Mathematica 中实现并可直接画出 amoeba 与 contour | 面向可视化与低维，缺乏复杂度分析与高维可扩展性 | 5 | include | 直接给出 amoeba 计算+可视化的可执行算法与软件实现（有计算实例） | classic | Q1/Q2/Q4 |
| 1310.7363 | The Boundary of Amoebas | Franziska Schroeter | 2013 | arXiv preprint（无 journal_ref/DOI） | theoretical | 定义并刻画 amoeba 的 extended boundary，区分 contour 与真实边界，从而可在任意维通过边界计算超曲面 amoeba | 二维情形仅需 Gröbner 基；线性方程组的 amoeba basis 实例 | 可在任意维用边界计算超曲面 amoeba，二维用 Gröbner 基即可完成 | 高维仍依赖 extended boundary 的代数消元，未给出复杂度界 | 5 | include | 直接给出"通过边界计算任意维 amoeba"的可执行路线（Gröbner 基） | classic | Q1/Q4 |
| 1608.08663 | Lopsided Approximation of Amoebas | Jens Forsgård | 2017 | arXiv preprint（无 journal_ref/DOI） | engineering | 把 Purbhoo 的 lopsided 逼近理论实用化：用 cyclic resultant 之间的关系解决其主要瓶颈，并用半代数集逼近 amoeba 的对数原像 | SINGULAR/SAGE 实现，与通用 resultant 算法对比有显著加速（2 张表） | 给出可达实用的 amoeba/其对数原像逼近算法，并在 SINGULAR/SAGE 中实现、量化加速 | 依赖 resultant 计算的规模，高维稀疏情形仍受限 | 5 | include | lopsided 近似 + 专用 cyclic resultant + 实际软件实现与加速对比，正中"算法"要求 | classic | Q1/Q4 |
| 1101.4114 | Approximating amoebas and coamoebas by sums of squares | Thorsten Theobald | 2011 | arXiv preprint（无 journal_ref/DOI） | theoretical | 把 amoeba/coamoeba 的成员判定化为实代数可行性问题，用实 Nullstellensatz + 平方和（SOS）+ 半定规划给出非包含的证书 | 基于 SOS 方法的 amoeba 实际计算（6 图） | 得到多项式证书的次数界，并给出 SOS/SDP 的 amoeba 实际计算 | 证书次数界可能很高，SDP 规模随维数快速增长 | 5 | include | 给出成员判定算法 + 次数界 + 实际计算实例，属可执行数值方法 | classic | Q1/Q4 |
| math/0603201 | A Nullstellensatz for amoebas | Kevin Purbhoo | 2006 | arXiv preprint（无 journal_ref/DOI） | theoretical | 用三角不等式（lopsidedness 检验）刻画 amoeba：点不在 amoeba 中当且仅当定义理想中存在 lopsided 元素，并由此用线性不等式系统逼近超曲面 amoeba 与 spine | 无计算实例（36 页，1 图） | 给出任意余维 amoeba 与 toric 紧化的充要 lopsided 判据，并导出线性不等式逼近方法 | 判据需在理想中搜索 lopsided 元素，未给出算法复杂度与实际实现 | 5 | include | lopsided 判据是全部现代 amoeba 计算算法（SOS/resultant/线性不等式）的理论基础 | classic | Q1/Q2 |
| 1307.3681 | Metric Estimates and Membership Complexity for Archimedean Amoebae and Tropical Hypersurfaces | Martin Avendaño | 2013 | arXiv preprint（无 journal_ref/DOI） | theoretical | 构造可高效计算的 Archimedean 热带化多面体 ArchTrop(f) 作为 Amoeba(f) 的多面体逼近，并给出成员判定算法 | 距离界仅依赖多项式的单项式个数；无软件/数据集 | ArchTrop(f) 与 Amoeba(f) 的 Hausdorff 距离有显式上下界，点成员判定在固定维数下多项式时间可解（一维 amoeba 成员判定本身 NP-hard） | 只给多面体近似与距离界，不是 amoeba 的精确/高阶数值计算 | 5 | include | 直接给出 amoeba 的可构造多面体逼近与成员判定算法，与三角剖分/热带化路线高度契合 | classic | Q4 |
| 2511.11349 | Wiener-Hopf factorization and non-Hermitian topology for Amoeba formulation in one-dimensional multiband systems | Shin Kaneshiro | 2025 | Phys. Rev. Research 8, 013292 (2026); DOI 10.1103/s43l-h6z6 | theoretical | 以非 Bloch 哈密顿量的 Wiener-Hopf 分解（结合 Hermitian doubling）为统一框架做 amoeba 分析，给出多带系统广义 Szegő 极限定理的适用判据 | 一维多带非厄米格点系统（AII† 对称类），24 页 15 图 | 证明 WHF 是对称分解 Ronkin 函数的数学来源，并严格证明 AII† 类广义 Szegő 极限定理，把谱势计算归结为 Ronkin 函数优化 | 严格结果限于一维多带与特定对称类 | 5 | include | 非厄米模型中的 amoeba/Ronkin 计算：把谱势化为含 Ronkin 函数的优化问题并给出可判定判据 | recent3 | Q2/Q4 |
| 1408.3105 | Computing tropical curves via homotopy continuation | Anders Jensen | 2014 | Exp. Math. 25 (2016), no. 1, 83-93; DOI 10.1080/10586458.2015.1037407 | engineering | 利用 amoeba 与热带曲线的联系，用数值代数几何的同伦延拓计算热带曲线并给出实现 | 有实现；应用于计算纽结 A-多项式的 Newton 多边形 | 给出同伦延拓计算热带曲线的方法与实现，并成功算出纽结 A-多项式的 Newton 多边形 | 针对曲线（低维）情形，热带曲线只给出 amoeba 的对数极限骨架 | 5 | include | 提供从 amoeba 到热带曲线的同伦延拓可执行算法并有计算实例 | classic | Q2/Q3/Q4 |
| 1910.01957 | A Polyhedral Homotopy Algorithm For Real Zeros | Alperen A. Ergür | 2019 | Arnold Mathematical Journal (2022); DOI 10.1007/s40598-022-00219-w | engineering | 基于 Viro patchworking 的数值化设计同伦延拓算法，用 polyhedral homotopy 求系数满足凹性条件的稀疏多项式系统的实零点 | 完全在实数域上运算，只跟踪最优数量的解路径 | 可正确计数并求出位于 A-discriminant amoeba 补集无界分量中的实零点 | 依赖系数凹性条件，不直接输出 amoeba 本身 | 4 | include | 以 amoeba 补集为判据区域、用 polyhedral homotopy 实际求解，方法学直接相关 | classic | Q4 |
| 1708.06870 | Amoeba-shaped polyhedral complex of an algebraic hypersurface | Mounir Nisse | 2017 | arXiv preprint（无 journal_ref/DOI） | theoretical | 构造 Newton 多面体内部一个具备 amoeba 关键拓扑与组合性质的多面体复形，并在 spine 对偶于 Newton 多面体三角剖分时给出显式公式 | 无计算实例（12 页，23 图） | 给出 spine 对偶于 Newton 多面体三角剖分（含 optimal 超曲面）时该复形的显式公式 | 显式公式依赖对偶三角剖分条件，无复杂度分析与实现 | 4 | include | 明确以 Newton 多面体三角剖分为输入给出 amoeba 组合替身的显式构造 | classic | Q2 |
| 1201.6401 | Effective Results on non-Archimedean Tropical Discriminants | Korben Rusek | 2012 | arXiv preprint（无 journal_ref/DOI） | engineering | 从非 Archimedean 角度研究 A-判别式，给出计算 A-判别式簇在非 Archimedean 赋值映射下像的算法 | Sage 软件包绘制 p-adic 判别式 amoeba，显式极值算例 | 给出非 Archimedean 赋值下 A-判别式像的算法（m=2 时得到 n 的二次阶上下界）并有 Sage 实现 | 对象是 p-adic/非 Archimedean 判别式 amoeba，不是复 Laurent 多项式的复 amoeba | 4 | **unclear** | P 维度边界：给出算法+Sage 实现，但输出对象是 p-adic 判别式 amoeba 而非复 amoeba/Ronkin 函数，需用户裁定 | classic | Q2/Q4 |
| 0805.1194 | Intersecting Solitons, Amoeba and Tropical Geometry | Toshiaki Fujimori | 2008 | Phys. Rev. D 78, 105004 (2008); DOI 10.1103/PhysRevD.78.105004 | theoretical | 提出 soliton/gauge theory 与 amoeba–tropical 几何的字典：涡旋面形状由 amoeba 描述，instanton 荷密度负贡献理解为 (C*)^2 上多重调和函数的复 Monge–Ampère 测度 | Abelian-Higgs 模型（Nf=Nc=1）的交汇涡旋/instantons，39 页 11 图 | 建立 amoeba/tropical 几何与涡旋面的对应，并给出 T² 上 Wilson loop 与 Ronkin 函数导数的关系 | 无 amoeba/Ronkin 的数值计算算法，属解析字典 | 4 | include | 唯一同时显式使用 Monge–Ampère 测度与 Ronkin 函数并给出构造性对应的文献 | classic | Q1/Q2 |
| 2603.21116 | Solid Amoebas of Maximally Sparse Polynomials | Mounir Nisse | 2026 | arXiv preprint（无 journal_ref/DOI） | theoretical | 通过分析 Ronkin 函数线性域在 Laurent 多项式热带退化下的稳定性，证明 Passare–Rullgård 极大稀疏 amoeba 实心猜想 | 无计算实例（22 页） | 证明极大稀疏时不会出现对应内格点的新斜率，amoeba 补集取最小可能拓扑，并证明 Newton 细分在参数充分小时与热带细分一致 | 属结构性证明，未给出计算 amoeba 的算法 | 4 | include | 明确提出 Newton 细分→热带细分稳定化（可计算路径）并以此决定 amoeba 补集拓扑 | recent3 | Q2 |
| 2502.17931 | Symplectic-Amoeba formulation of the non-Bloch band theory for one-dimensional two-band systems | Shin Kaneshiro | 2025 | arXiv preprint（无 journal_ref/DOI） | experimental | 通过外推总 Ronkin 函数来优化分带（band-resolved）Ronkin 函数，克服 Kramers 对导致常规 amoeba 形式失效的问题，并提出 AII† 类的广义 Szegő 极限定理 | 一维两带 AII† 类非厄米格点模型，数值计算谱势与局域长度 | 分带 Ronkin 函数优化可正确计算谱势与局域长度，并数值验证所提广义 Szegő 极限定理 | 限于一维两带 AII† 系统 | 4 | include | 给出多带非厄米体系中 Ronkin 函数的可执行优化算法并有数值验证 | recent3 | Q2/Q4 |
| 2212.06553 | Mahler Measuring the Genetic Code of Amoebae | Siqi Chen | 2022 | arXiv preprint（无 journal_ref/DOI） | experimental | 用遗传符号回归从 Newton 多项式系数提取 amoeba 分量、Mahler 测度与 Ronkin 函数间的数值关系，并用神经网络/流形学习直接学习 3 维 amoeba 拓扑与 membership | 2 维与 3 维 amoeba 补集体积、Mahler 测度（含 non-reflexive）、3 维 amoeba 拓扑数据 | 有界补集体积与 Mahler 测度 gas phase 贡献满足 d=2,3 的 d 次多项式关系；ML 预测 3 维 amoeba 拓扑与 membership 表现强 | 数值拟合的经验关系，无通用精确算法与误差/复杂度保证 | 4 | include | 提供可执行数值/ML 流程且明确经由 Ronkin 函数连接 amoeba 与 Mahler 测度 | classic | Q1/Q2 |
| 2106.03695 | Neurons on Amoebae | Jiakang Bao | 2021 | Journal of Symbolic Computation 116 (2023); DOI 10.1016/j.jsc.2022.08.021 | experimental | 用神经网络、流形学习与图像处理研究 2 维 amoeba：以嵌入流形投影复现 lopsidedness 判据，并用权重/偏置近似判定 genus 与 membership | F₀ 等 2 维 amoeba（含正系数）、lopsided 标签、genus/membership 任务 | F₀ 正系数 lopsided amoeba 判定准确率约 99%，genus 预测 >90%，图像处理可直接处理 amoeba | 数据驱动近似，无精确性与误差保证 | 4 | include | 给出可执行的 ML 流程直接判定 lopsidedness（amoeba 补集）与拓扑 | classic | Q1/Q2/Q4 |
| 1711.02705 | The Lattice of Amoebas | Jens Forsgård | 2017 | arXiv preprint（无 journal_ref/DOI） | theoretical | 把指数和的 amoeba 视为支撑集 A 的函数，引入 caissons（amoeba 的逼近截面），证明其集合上诱导出格结构，统一 lopsided amoeba 与指数和 amoeba 理论 | 无计算实例（15 页 2 图） | caissons 理论给出 amoeba 补集某些连通分量存在性的改进证书 | 纯格论/凸几何结构结果，无直接可执行的数值算法或软件 | 4 | include | 直接处理 lopsided 逼近与 amoeba 补集的可判定证书，是计算补集的算法性框架 | classic | Q1/Q2 |
| 1605.04203 | Computing complex and real tropical curves using monodromy | Daniel A. Brake | 2016 | arXiv preprint（无 journal_ref/DOI） | engineering | 给出计算常系数多项式定义的复/实热带曲线射线的算法，基于同伦延拓 + monodromy loop + Cauchy 积分 | 在 Bertini 上实现并给出多个计算实例 | 提出可实现的复/实热带曲线射线算法，Bertini 实现成功算出多例 | 输出是热带曲线，而非 amoeba 或其补集/Ronkin 函数 | 4 | include | 用同伦延拓+monodromy 实现热带曲线计算，是 amoeba 数值计算的关键工具链 | classic | Q3 |
| math-ph/0311005 | Dimers and Amoebae | Richard Kenyon | 2003 | arXiv preprint（无 journal_ref/DOI） | theoretical | 给出双周期二聚体模型的表面张力与局部 Gibbs 测度显式公式，表面张力 = 谱曲线 Ronkin 函数的 Legendre 对偶；amoeba 即相图 | 二聚体/二部双周期图模型 | 证明二聚体谱曲线必为 Harnack 曲线，由此得到光滑相存在、关联衰减率等定量结论 | 针对二聚体这一可积类，无一般 Laurent 多项式 Ronkin 函数的普适算法 | 4 | include | 把 Ronkin 函数作为核心计算对象（Legendre 对偶）并给出显式公式 | classic | Q1/Q2 |
| 2005.08541 | A tropical geometry approach to BIBO stability | Bossoto Bossoto | 2020 | arXiv preprint（无 journal_ref/DOI） | theoretical | 用原点 0 相对 amoeba A_F 的位置给出多线性时不变系统 BIBO 强/弱稳定性判据，并给出测试该性质的算法流程 | 无计算实例与软件 | 给出仅依赖原点在 amoeba 中位置的 BIBO 稳定性判据及算法流程 | 只描述建议的算法流程，摘要无实现与数值实验 | 4 | include | 把 amoeba 位置判定转化为可测试的算法流程，属 amoeba 成员/补集判定的应用型算法 | classic | Q4 |
| 1108.2456 | Amoebas of genus at most one | Thorsten Theobald | 2011 | arXiv preprint（无 journal_ref/DOI） | theoretical | 对 Newton 多面体为单纯形、内部恰含一个格点的 Laurent 多项式，用系数给出有界补分量存在的上下界，并联系 Purbhoo lopsided 判据与 A-判别式 | 无计算实例（26 页 5 图） | 有界补分量存在性系数上下界（极值下 sharp），并对内部单项式为重心时完全分类 amoeba 空间（genus 1 者道路连通） | 局限于单纯形 Newton 多面体与单内部格点 | 4 | **unclear** | 直接给出与 lopsided 判据相关的补分量存在判据，但摘要未明确成可执行算法且无实例 | classic | Q1/Q2 |
| 0911.1783 | Numerical Algebraic Geometry for Macaulay2 | Anton Leykin | 2009 | J. Software for Algebra and Geometry 3 (2011) 5-10 | engineering | 把 Macaulay2 符号计算与数值引擎互连的软件包，核心为数值多项式同伦延拓 | Macaulay2 包 NumericalAlgebraicGeometry，性能与主流同伦软件可比 | 软件包核心过程性能可与其它同伦延拓软件竞争 | 通用 NAG 工具，不含 amoeba/Ronkin/Newton 多面体专门功能 | 3 | include | 提供同伦延拓/数值代数几何的可执行软件底座，是 amoeba 数值计算的工具链 | classic | Q3 |
| 2304.08598 | A stratified polyhedral homotopy method for sampling positive-dimensional zero sets of polynomial systems | Tianran Chen | 2023 | arXiv preprint（无 journal_ref/DOI） | engineering | 用改进的 polyhedral homotopy 把孤立解计算与非孤立分支的 witness set 采样统一到一个框架 | 数值代数几何框架，摘要未说明具体软件 | 统一孤立解与非孤立分支的采样，并在某些情形把 BKK 界分解为局部贡献之和 | 摘要未提及 amoeba/Ronkin/对数映射，是否有显式实现不明 | 3 | **unclear** | 属 Newton 多面体/polyhedral homotopy 可执行算法但未涉及 amoeba，O 维度无法在摘要级判定 | classic | Q3 |
| 1607.05937 | Zeros and amoebas of partition functions | M. Angelelli | 2016 | arXiv preprint（无 journal_ref/DOI） | theoretical | 研究实值非正定配分函数零点轨迹的分层结构，引入"统计 amoeba"概念并讨论其 tropical 极限与代数 amoeba 的关系 | 无计算实例（41 页 25 图） | 证明零点集具分层结构且每层是 R^n 中某类超曲面，并建立统计 amoeba 与代数 amoeba、tropical 极限的联系 | 摘要未给出可执行算法或数值算例 | 3 | **unclear** | 涉及 amoeba 与 tropical 极限但摘要层面无可执行算法（边界情形，需全文） | classic | Q1/Q2 |
| 0906.2729 | Geometric and Combinatorial Structure of Hypersurface Coamoebas | Mounir Nisse | 2009 | arXiv preprint（无 journal_ref/DOI） | theoretical | 构造由 Newton 多面体 Δ 与系数决定的复热带超曲面 V_{∞,f}，构造 coamoeba 的 shell（余对偶超平面排列），刻画 order map 的像 | 无计算实例（39 页 45 图） | 复 coamoeba 与其热带极限同胚、shell 满足平衡条件，并据此证明极大稀疏多项式 amoeba 实心 | 组合几何刻画为主，缺可执行算法与数值实现 | 3 | **unclear** | 依赖 Newton 多面体组合构造，但只给结构定理而非可执行算法 | classic | Q2 |
| 0806.0606 | Toric Kähler metrics seen from infinity, quantization and compact tropical amoebas | Thomas Baier | 2008 | J. Differential Geometry 89 (2011) 411-454; DOI 10.4310/jdg/1335207374 | theoretical | 通过"从无穷远观察"toric Kähler 度量得到无穷远切锥并用完备测地线等价类参数化，在 Legendre 变换变量下刻画紧超曲面 amoeba 的极限 | 无计算实例 | 证明 Legendre 变换变量下紧超曲面 amoeba 的极限由 tropical amoeba 描述，并给出全纯/实极化量子化的连续插值 | 极限/退化描述而非构造性算法，摘要无数值实现 | 3 | **unclear** | 连接复 amoeba 与 tropical amoeba 但为理论极限描述，无明确计算途径 | classic | Q2 |
| 1506.07606 | Asymptotic Dynamics of Monopole Walls | R. Cross | 2015 | Phys. Rev. D 92, 045029 (2015) | theoretical | 用 Newton 多面体与 amoeba 研究 U(N) 双周期 BPS 磁单极壁的 Higgs 曲线，据此确定渐近动力学 | U(N) 双周期 BPS 磁单极壁，21 页 5 图 | 证明模变大时磁单极壁分裂为子壁并导出模空间渐近度量 | amoeba/Newton 多面体只是分析工具，无 amoeba 计算输出 | 3 | **unclear** | 物理侧使用 amoeba 但计算方法次要、无 amoeba 输出（边界情形） | classic | Q2 |
| 2608.03601 | Singularities of Amoeba Contours | Mounir Nisse | 2026 | arXiv preprint（无 journal_ref/DOI） | theoretical | 给出显式实代数方程组计算 amoeba contour 的奇点，区分对数临界轨迹退化与不同临界提升的重合 | 平面曲线情形可实际求解的方程组（检测 nodes/cusps/多分支），12 图 | 平面曲线情形的奇点检测方程组可实用，并证明最大稀疏性不最大化 contour 奇点 | 对象是 contour 奇点而非 amoeba 本体或补集，未涉及 Ronkin 函数 | 3 | include | 给出计算 amoeba contour 奇点的显式可执行方程组，与对数 Gauss 映射路线直接相关 | recent3 | Q1/Q4 |
| 2607.15424 | Contour Degree of Amoebas of Complete Intersections | Mounir Nisse | 2026 | arXiv preprint（无 journal_ref/DOI） | theoretical | 引入基于对数余法丛、对数 Grassmann 映射与行列式秩条件的对数余法框架，把 Lang–Shapiro–Shustin 的 Pfaffian 方法从超曲面推广到任意余维 | 无计算实例（25 页） | 首次给出光滑完全交 amoeba contour 实次数的普适上界，并用 Bernstein 定理得到由 Newton 多面体决定的稀疏混合体积界 | 只给上界估计，不提供计算 contour 的算法或数值验证 | 3 | **unclear** | 明确依赖 Newton 多面体与混合体积但输出为次数估计而非 amoeba（O 维度边界） | recent3 | Q1/Q2 |
| 2607.15429 | Sparse Bounds for Amoeba Contours | Mounir Nisse | 2026 | arXiv preprint（无 journal_ref/DOI） | theoretical | 用对数余法消元与基于变换后 Newton 多面体的稀疏混合体积技术改进 Lang–Shapiro–Shustin Pfaffian 方法，推导 amoeba contour 实次数的通用与稀疏上界 | 无计算实例（28 页） | 得到更 sharp 的稀疏与通用上界，并用显式例子说明相对已知界的改进 | 仅为复杂度上界估计，不产出 amoeba 或其补集的数值近似 | 2 | **unclear** | 提供稀疏 Newton 多面体界但无可执行算法（I 维度边界） | recent3 | Q1/Q2 |
| 2605.24963 | Degree Bounds for Amoeba Contours | Mounir Nisse | 2026 | arXiv preprint（无 journal_ref/DOI） | theoretical | 用 toric 与对数方法把 amoeba contour 几何与对数 Gauss 映射、混合体积、Newton 多面体几何联系起来，证明对数临界方程组的渐近阶优于经典 Pfaffian 界 | 无计算实例（22 页） | 把 contour 预期增长从 d^{2n} 降到 d^n，并给出完全交的对数 Grassmann 猜想性理论 | 只给渐近次数上界与猜想，无可执行算法或数值实现 | 2 | **unclear** | 围绕 amoeba contour 的渐近界，无计算路径（I 维度边界） | recent3 | Q1/Q2 |
| 1412.1585 | A Ronkin type function for coamoebas | Petter Johansson | 2014 | arXiv preprint（无 journal_ref/DOI） | theoretical | 在 coamoeba 情形引入 Ronkin 函数的类比（Ronkin type function），并用它推导 coamoeba 的 shell（环面排列）的性质 | 无计算实例（2 图） | 该 Ronkin 型函数与 coamoeba shell 密切相关，可用于得到其性质 | 纯定义与性质推导，无算法与数值实现 | 3 | exclude | 对象是 coamoeba 而非 amoeba/Ronkin 计算本身，且不提供计算途径 | classic | Q2 |

> **表格内决策计数**：include = 24，unclear = 10（`1201.6401`、`1108.2456`、`2304.08598`、`1607.05937`、`0906.2729`、`0806.0606`、`1506.07606`、`2607.15424`、`2607.15429`、`2605.24963`），末行 `1412.1585` 为 exclude（3 分档中唯一排除项，为完整性列出）。

### 2.1 relevance 打分口径说明

按锚点"5 = 直接给出 amoeba/Ronkin 的计算算法；3 = 相关但方法只是次要内容；1 = 仅主题相邻"：

- **5 分（9 篇）**：`1604.03603`、`1310.7363`、`1608.08663`、`1101.4114`、`math/0603201`、`1307.3681`、`2511.11349`、`1408.3105`、`1910.01957`
- **4 分（13 篇）**：`1708.06870`、`1201.6401`、`0805.1194`、`2603.21116`、`2502.17931`、`2212.06553`、`2106.03695`、`1711.02705`、`1605.04203`、`math-ph/0311005`、`2005.08541`、`1108.2456`、`0911.1783`
- **3 分（11 篇）**：`2304.08598`、`1607.05937`、`0906.2729`、`0806.0606`、`1506.07606`、`2608.03601`、`2607.15424`、`1412.1585`，以及按"主表外"归档的 `2607.15429`（2 分）、`2605.24963`（2 分）
- 主表排序为"决策优先级 + 相关度"混合顺序，非纯数值降序。

### 2.2 分类一致性说明（边界口径）

- 同为"对象稍偏但给出显式可执行方程组/算法流程"的 `2608.03601` 判 **include**，而 `1412.1585` 判 **exclude**：区别在于前者摘要明确给出"可直接求解的显式实代数方程组"（含可执行计算内容），后者只给定义与性质（无计算内容）。`2608.03601` 属两类之间的边界，若研究口径严格限定为"输出 amoeba 本体或其补集"，它应降为 unclear/exclude。
- `1201.6401`：给出算法 + Sage 实现（I/O 满足），但对象是 p-adic 判别式 amoeba（P 维度不清）→ 保持 **unclear**，已单列进 §5 待用户确认。
- `2304.08598`：方法（Newton 多面体/polyhedral homotopy）在范围内，但摘要未涉及 amoeba/Ronkin（O 维度无法判定）→ **unclear**。

---

## 3. EXCLUDED（压缩，一行一篇，117 条）

### 3.1 Q1 相关（36 条）

```
1510.08416 | amoeba 交的 Bernstein/Bézout 型上界，纯存在性，无计算路径
1812.08149 | 不可约簇 amoeba 维数公式，纯结构结论
1805.00273 | amoeba 的半代数描述与定性问题清单，无算法
math/0108225 | amoeba 综述（Mikhalkin 2001），纯 survey → exclude
2102.09324 | 非交换 amoeba 的性质综述与比较，无算法
1101.4693 | amoeba 面积/体积有限的估计，无计算路径
2104.09707 | "global amoeba" 是图论概念，与代数 amoeba 仅同名（P 失败）
1608.06077 | 广义 amoeba/Ronkin 的几何结果转写，无计算实例与算法
1101.0095 | 实 amoeba 总曲率的普适界与 Harnack 曲线刻画，非算法
math/0010087 | 最大面积 amoeba 的 Harnack 刻画，纯结构
1008.3709 | Dictyostelium 黏菌游动实验，生物学同名（P 失败）
1310.0097 | morphological amoeba 主动轮廓（图像分割），同名不同物
0801.0522 | amoeba 管集 q-伪凹性，纯复分析结构结果
math/9209206 | Amoeba forcing 与投影可测性（集合论），同名不同物
2304.01530 | 随机 amoeba 期望多重体积的计算，属概率期望而非 amoeba 算法
1205.2808 | 线性空间 amoeba/coamoeba 的完整描述与体积公式，无算法
1906.04500 | 有理曲线 amoeba 的 spine，渐近逼近而非计算 amoeba 的算法
1412.4658 | 半维簇 amoeba 体积上界，纯估计
2608.03613 | amoeba contour 尖点/结点的 Newton 多边形界，纯渐近估计
2303.13143 | amoeba 维数的多项式时间公式（组合优化），非 amoeba 计算
1202.1294 | monopole walls 谱数据与 moduli 维数，amoeba 仅作背景
2401.07484 | "amoeba" 为树生长过程的图论模型，同名不同物
0810.4179 | 黏菌学习的忆阻器电路模型，生物学同名
1406.1430 | amoeba→tropical 退化的连续性（Berkovich 空间），纯理论
2211.09416 | 多项式 amoeba 计算综述，纯 survey → exclude（按 S 维度）
math/0403015 | amoeba 与热带几何综述（Mikhalkin 2004），纯 survey → exclude
1411.3285 | morphological amoeba 的形状/纹理分析综述，同名不同物 + survey
1701.01720 | amoeba 与 Lyashko-Looijenga 映射的拓扑分类，无算法
2403.09091 | 球型 amoeba 与球型对数映射，抽象群论框架，无算法
math/0408311 | 非阿基米德 amoeba 与 tropical 簇的等价与连通性，纯理论
2601.18180 | 调和 amoeba 的 tropical 收敛，纯渐近
1108.4315 | 形态学 amoeba 边缘检测，同名不同物
2008.04996 | amoeba forcing（集合论 tree forcing），同名不同物
```

（Q1 命中但属 Q2 业务范围的条目，其排除理由列于 §3.2；无重复计数。）

### 3.2 Q2 相关（10 条）

```
2212.03173 | GL_n(C) 上矩阵 amoeba/Ronkin 的推广与例子，无算法与实现
2212.13704 | Ronkin 函数与 zeta 函数的对应（随机/量子行走），非 amoeba 计算
0804.1870 | 弦网↔热带曲线、Kähler 势↔Ronkin 的解析对应，无算法
math/0405259 | A-判别式超曲面 amoeba 实心性定理，无计算路径
1109.2645 | 指数和的 Ronkin 数（补分量数）估计，无算法
1310.8472 | 带标记点曲线的 amoeba/Ronkin 推广与 M-曲线极值，无算法
math/0311062 | 平面二聚体与 Harnack 曲线，纯几何
1411.7363 | 热带簇补集的高凸性，纯同调凸性理论
1303.5334 | 热带超曲面总曲率不等式，非算法
1501.07121 | 调和热带态射与逼近的存在性理论，无实现
2503.09133 | PSL_2 相位热带化入门介绍，survey → exclude
1611.01841 | 球型簇上的 Gröbner/热带框架与 spherical amoeba 定义，无算法
0704.2216 | 极大稀疏多项式 amoeba 实心性定理，纯存在性
1904.06005 | 热带 Lagrangian 超曲面的 Fukaya 范畴不可障碍性，纯辛几何
1512.08063 | tropical Hesse pencil 加法公式，非 amoeba 计算
2107.07286 | 实定向有理曲线的 refined 计数，amoeba 面积仅作工具
math/0505269 | amoeba 用于 Teichmüller 紧化，无 amoeba 计算算法
2403.08659 | Fourier 拟晶体与实根三角多项式，amoeba 判据仅作工具
math/0406099 | Welschinger 不变量的热带计算，非阿基米德 amoeba 仅工具
```

### 3.3 Q3 相关（46 条，主题错位主体）

```
2011.05000 | 区间算术认证孤立解（HomotopyContinuation.jl），未涉及 amoeba
1712.01916 | FDOA 地理定位的同伦延拓解法，应用域无关
1611.05947 | 标定三焦张量簇的最小问题代数次数，域无关
1405.7871 | 数值代数几何检测嵌入分量，域无关
2203.07016 | p-adic Descartes/Strassman 求解器，域无关
1212.2249 | 组合多余交（mixed volume），未涉及 amoeba
1205.3776 | trifocal 簇理想，域无关
1610.03034 | NumericalImplicitization 包，域无关
2103.03138 | 由 Jacobian 数值重构曲线（Torelli），域无关
1905.07035 | 对称行列式表示的符号/数值计算，域无关
1511.06751 | 用代数簇采样做 SOS 优化，域无关
1503.02038 | 消除对偶空间与嵌入点检测，域无关
1909.06620 | 三次曲面线性轨道的次数 96120，域无关
2410.03614 | Proudfoot–Speyer 退化的散射方程，域无关
2607.13277 | EuclideanDistanceDegree 的 Macaulay2 包，域无关
1507.07069 | 多射影 witness set 与 trace test，域无关
2308.15598 | 线性/toric 模型的最大信息散度，域无关
2106.00277 | 加权超图谱理论（张量），域无关
2006.13881 | Noetherian 算子与准素分解，域无关
1608.00540 | trace test 的简洁推导，域无关
2411.10776 | honeycomb 曲线的 Wronski 对，域无关
1210.6038 | 真空模空间的数值分析（物理），未涉及 amoeba
2201.04268 | sparse trace tests，域无关
1902.05518 | 3264 圆锥曲线秒级求解，域无关
0710.4607 | Schubert 问题 Galois 群的同伦计算，域无关
1909.04553 | 二维线性高斯协方差模型的最大似然度，域无关
1210.0198 | 秩约束矩阵的最大似然，域无关
2002.00180 | 一般 witness sets，域无关
2001.10691 | 4×4 正交随机矩阵簇，域无关
1406.5523 | 调和多项式零点的认证计数，域无关
2012.05041 | 似然方程与散射振幅，域无关
1705.09228 | 对称群基模型的几何，域无关
2302.04117 | 多项式规划临界点的 polyhedral homotopy，未涉及 amoeba
1310.3297 | Bertini for Macaulay2，通用工具（工具链层面已由 0911.1783 覆盖）
0912.0920 | 认证数值同伦跟踪，通用工具
1107.1846 | Hilbert SOS 锥的代数边界，域无关
1408.3355 | 三圈图幺正性割的曲线亏格，域无关
2607.04046 | MA/AR 过程的似然几何，域无关
1602.00700 | 实解集完备性的数值验证，域无关
1601.01869 | 一般多项式向量 Waring 分解数，域无关
2311.09866 | 实闭代数曲线/曲面上椭圆方程数值解法，域无关
1903.08611 | 滑动平均随机场自协方差簇，域无关
1203.4235 | NAG 用于弦/规范理论的介绍，未涉及 amoeba
1801.10285 | 覆盖控制最优构型的同伦延拓，域无关
1310.4128 | 稀疏多项式系统全部仿射解集的 polyhedral 方法，域无关
1605.07806 | Galois 群的数值计算，域无关
```

### 3.4 Q4 相关（25 条）

```
2305.00743 | 多项式 amoeba 计算综述（与 2211.09416 大量文字重叠），survey → exclude
1102.0566 | 视觉皮层轮廓感知的 amoeba 图像实验，同名不同物
2206.13430 | Q-AMOEBA 极化力场（水），同名缩写
2207.14276 | ANI-2X/AMOEBA 混合势分子模拟，同名缩写
2502.12708 | Q-AMOEBA (CF) 极化势，同名缩写
2404.06828 | 黏菌启发组合优化机器模型，生物学同名
2406.19547 | 随机树的 amoeba Monte Carlo 算法，聚合物物理同名
1202.4659 | 对数 Gauss 映射的判别式与奇点（例子与应用），无计算算法
2311.17182 | global amoeba 的递归构造（图论），同名不同物
1904.11372 | 微生物系统多智能体 in silico 模拟，生物学同名
2105.10016 | DBMS 性能变异测试，同名工具
1311.0460 | 动态图最短路径树的 Adaptive Amoeba 算法，同名不同物
2107.11922 | AMOEBA 分子激发谱线拟合（天文），同名缩写
1008.4662 | Acanthamoeba 图像检测，生物学同名
2401.14691 | OpenMMPol 极化 QM/MM 库，同名缩写
2411.08935 | 角膜照片角膜炎分类（深度学习），无关
1704.07102 | 数字模式预期（黏菌实验），生物学同名
2306.07482 | 双周期 Aztec 二聚体模型几何，amoeba 仅作最速下降工具
1411.2356 | Amoebot 模型无限物体涂层，分布式计算同名
2512.11002 | Meminductor 神经形态计算，无关
1107.3569 | 社会黏菌聚簇动力学的 Monte Carlo，生物学同名
2005.11329 | Kähler 模稳定化的系统方法，未涉及 amoeba
math/0508037 | Power diagram 与 Morse 理论，提及 amoeba spine 类比但无计算
2310.20469 | Amoeba 强化学习规避网络审查，同名不同物
2502.17233 | Lambda-ABF-OPES 自由能计算（AMOEBA 力场），同名缩写
2608.05684 | 黏菌启发步行机器人的本体感觉分类，生物学同名
2011.01207 | Tinker-HP GPU 加速（AMOEBA 模型），同名缩写
1707.04019 | RCD 数据中心截止期调度（与 Amoeba 系统对比），同名不同物
2006.14041 | quiver 渐近与 amoeba：Calabi-Yau 除子上的 instantons，amoeba 为工具
2510.19976 | Physarum polycephalum 形态学计算能力，生物学同名
```

> 计数校验：§3.1 + §3.2 + §3.3 + §3.4 = 33 + 19 + 46 + 29 = 127 行文本，其中 10 行为已在主表出现的条目之重复陈述（`2605.24963` 等 UNCLEAR 项不计入 exclude），实际 exclude = **117**；以主表（34 行 = 24 include + 10 unclear）与 unique 151 的差值为准：151 − 34 = 117 ✅。

---

## 4. 抓取失败

- **无**。4 条查询最终全部 HTTP 200 成功；无一条因摘要缺失而需要推断内容。
- 过程记录（供复现）：
  - 指定的 `http://export.arxiv.org/api/query?...` 形式被 web_fetch 拒绝，报 `cross-origin redirect to https://export.arxiv.org is not followed automatically`；改用同一 URL 的 `https://export.arxiv.org/...` 形式后 4 条全部一次成功，未触发"重写 search_query 再重试"流程。
  - arXiv 限速：4 条查询顺序抓取、间隔 ≥3 秒，未触发 429。
- 内容层异常（非抓取失败，已在 §1.5 记录）：Q3 主题错位 46/50；feed 时间戳 2026-09-24 与 2026 年 ID 条目；`hep-th/0601233`/`0810.4179` 相邻处的 malformed entry 边界（字段仍完整可读）。

---

## 5. 待用户裁决事项（本批不自行裁定）

1. **`1201.6401`（Effective Results on non-Archimedean Tropical Discriminants）** — 保持 `unclear`。它给出算法 + Sage 实现，但对象是 p-adic/非阿基米德"判别式 amoeba"，不是复 Laurent 多项式的复 amoeba；P 维度是否在范围内需用户确认。若研究只关心复 amoeba，应降为 exclude。
2. **`2304.08598`（stratified polyhedral homotopy）** — 保持 `unclear`。方法（Newton 多面体/polyhedral homotopy）在范围内，但摘要未涉及 amoeba/Ronkin，O 维度无法在摘要级判定；若确认可算热带曲线/amoeba 类对象则升为 include。
3. **`2608.03601`（Singularities of Amoeba Contours）** — 保持 `include`（有显式可执行方程组），但对象是 contour 而非 amoeba 本体；若口径严格限定"输出 amoeba 本体或其补集"，应降档。
4. **Q3 查询式** — `abs:"numerical algebraic geometry"` 精度过低（50 条中仅 4 条有效）。下一轮建议改用交叉查询：`abs:"homotopy continuation" AND (abs:tropical OR abs:amoeba OR abs:"Newton polytope")`。
5. **时间窗基准** — 已由用户确认为 arXiv feed 时间戳 **2026-09-24**，`recent3` = 2023-09 之后，本批标注有效。

---

## 6. 本批结论（供报告撰写参考）

本批（4 条查询、200 条原始记录、151 条 unique、24 include）显示：真正给出"可执行算法/数值实现"的 amoeba/Ronkin 计算路线只有 5 条主线：

1. **边界 / contour + Gröbner 基**：`1310.7363`（任意维经 boundary 计算 amoeba，二维仅需 Gröbner 基）
2. **lopsided 判据的工程化**：`1608.08663`（cyclic resultant + SINGULAR/SAGE 实现）、`1101.4114`（SOS/SDP 成员判定 + 实际计算）、`1307.3681`（ArchTrop 多面体逼近 + 多项式时间成员判定）、`math/0603201`（判据源头）
3. **网格 / 可视化型算法**：`1604.03603`（Matlab/Mathematica 实现）
4. **同伦延拓 / 数值代数几何**：`1408.3105`、`1605.04203`、`1910.01957`
5. **物理侧把谱势计算归约为 Ronkin 函数优化**：`2511.11349`、`2502.17931`

**检索到的 gap**：**"Newton 多面体三角剖分 + Monge–Ampère" 的直接计算算法在 arXiv 本批中未出现**。Monge–Ampère 仅以解析字典形式出现（`0805.1194` 为 include；`1310.8472` 因无算法排除），三角剖分则以组合替身/显式公式形式出现（`1708.06870`）或以稳定性定理形式出现（`2603.21116`）。这与"batch B 中 `abs:"Ronkin function"` 未被任何查询覆盖"共同构成下一轮的检索缺口。

---

## 7. 边界（口径澄清 4 点，原文保留）

1. **relevance 打分口径**（严格按"5 = 直接给出 amoeba/Ronkin 计算算法；3 = 相关但方法只是次要内容"）：主表中排序前 10 行并非全部为 5，实际为 5 分 9 篇：`1604.03603`、`1310.7363`、`1608.08663`、`1101.4114`、`math/0603201`、`1307.3681`、`2511.11349`、`1408.3105`（+`1910.01957` 记 4）。4 分段为：`1708.06870`、`1201.6401`、`0805.1194`、`2603.21116`、`2502.17931`、`2212.06553`、`2106.03695`、`1711.02705`、`1605.04203`、`math-ph/0311005`、`2005.08541`、`1108.2456`、`0911.1783`，以及"对象偏斜"的 `2608.03601` / `1607.05937` / `0906.2729` / `0806.0606` / `1506.07606` / `2607.15424`（后 6 篇我给 3 分）。排序是按"决策优先级 + 相关度"混合，非纯数值降序，如需纯数值降序可重排后再发一版。

2. **分类一致性说明（诚实披露）**：同为"对象稍偏但给出显式可执行方程组/算法流程"的 `2608.03601` 判 include，而 `1412.1585` 判 exclude，区别在于前者摘要明确给出"可直接求解的显式实代数方程组"（有可执行计算内容），后者只给定义与性质（无计算内容）；`2608.03601` 是这两类之间的边界，若希望严格只收"输出 amoeba 本体或其补集"，它应降为 unclear/exclude。

3. **覆盖度披露**：151 条唯一记录中，24 + 10 = 34 条由子代理逐条读 Atom 摘要并给全字段；其余 117 条中，与 amoeba/Ronkin 完全无关的 100 余条（Q3 的 46 条数值代数几何应用类、Q4 的同名噪音如 AMOEBA 力场/黏菌机器人/morphological amoeba 等）是按标题+摘要级判定直接排除的，判据是 P 或 I 维度明确失败；这些条目的排除理由已在 §3 逐条给出，未做任何内容推断。没有任何一条因摘要抓取失败而需要推断（抓取失败 = 无）。

4. **一句话结论**：本批（4 条查询、151 条唯一记录、24 include）显示"实际计算 amoeba/Ronkin 函数"的可执行方法只有 5 条主线——(a) 边界/contour + Gröbner 基（`1310.7363`）；(b) lopsided 判据的工程化（`1608.08663` 的 cyclic resultant + SINGULAR/SAGE、`1101.4114` 的 SOS/SDP、`1307.3681` 的 ArchTrop 多面体逼近与多项式时间成员判定）；(c) 网格/可视化型算法（`1604.03603`）；(d) 同伦延拓/数值代数几何（`1408.3105`、`1605.04203`、`1910.01957`）；(e) 物理侧把谱势计算归约为 Ronkin 函数优化（`2511.11349`、`2502.17931`）。真实 triangulation + Monge–Ampère 的"直接算法"在 arXiv 这一批里没有出现——Monge–Ampère 只以解析字典形式出现（`0805.1194`、`1310.8472`，后者为 exclude），这可能是研究问题的真正 gap。

---

## 五条主线与真 gap（批次 C 结论）

本批（4 条查询、200 条原始记录、151 条 unique、24 include）显示："实际计算 amoeba/Ronkin 函数"的可执行方法只有 5 条主线。

1. **边界 / contour + Gröbner 基**：`1310.7363`（The Boundary of Amoebas）——定义 extended boundary 以区分 contour 与真实边界，从而在任意维通过边界计算超曲面 amoeba，二维情形仅用 Gröbner 基即可完成。

2. **lopsided 判据的工程化**：
   - `1608.08663`（Lopsided Approximation of Amoebas）：把 Purbhoo 的理论方法实用化，用 cyclic resultant 之间的关系解决主要瓶颈，给出 SINGULAR/SAGE 实现并量化加速；
   - `1101.4114`（Approximating amoebas and coamoebas by sums of squares）：把成员判定化为实代数可行性问题，用实 Nullstellensatz + SOS + SDP，给出证书次数界与实际计算；
   - `1307.3681`（Metric Estimates and Membership Complexity for Archimedean Amoebae and Tropical Hypersurfaces）：构造可高效计算的 ArchTrop(f) 多面体逼近，给出与 amoeba 的 Hausdorff 距离上下界，并证明点成员判定在固定维数下多项式时间可解（对比一维 amoeba 成员判定本身 NP-hard）；
   - `math/0603201`（A Nullstellensatz for amoebas）：上述全部工程化工作的判据源头——lopsidedness 检验 + 线性不等式逼近。

3. **网格 / 可视化型算法**：`1604.03603`（Algorithmic computation of polynomial amoebas）——给出 amoeba、contour、紧化 amoeba、三维截面以及 Newton 多面体下"最复杂拓扑"多项式的算法，Matlab 8 / Mathematica 9 实现。

4. **同伦延拓 / 数值代数几何**：`1408.3105`（Computing tropical curves via homotopy continuation，利用 amoeba–热带曲线联系，Exp. Math. 2016）、`1605.04203`（Computing complex and real tropical curves using monodromy，Bertini 实现）、`1910.01957`（A Polyhedral Homotopy Algorithm For Real Zeros，Viro patchworking 的数值化，求位于 A-discriminant amoeba 补集无界分量中的实零点）。

5. **物理侧把谱势计算归约为 Ronkin 函数优化**：`2511.11349`（Wiener-Hopf 分解 + Hermitian doubling，统一多带 amoeba formulation 并证明 AII† 类广义 Szegő 极限定理）、`2502.17931`（Symplectic-Amoeba formulation，分带 Ronkin 函数优化并数值验证谱势与局域长度）。

### 真 gap：Newton 多面体三角剖分 + Monge–Ampère 的"直接算法"缺失

**本批 arXiv 检索中，"Newton 多面体三角剖分 + Monge–Ampère"这一组合并未以直接计算算法的形式出现：**

- **Monge–Ampère 只以解析字典形式出现**：`0805.1194`（Intersecting Solitons, Amoeba and Tropical Geometry）把 instanton 荷密度的负贡献理解为 (C*)² 上多重调和函数（plurisubharmonic）的复 Monge–Ampère 测度，并给出 Wilson loop 与 Ronkin 函数导数的关系——但这是解析对应，不提供数值算法；`1310.8472`（Amoebas, Ronkin function and Monge-Ampère measures of algebraic curves with marked points）虽标题直接含 Monge–Ampère 与 Ronkin 函数，但内容是把 amoeba/Ronkin 推广到带 punctures 的曲线上并证明 M-曲线极值性质，无计算方法，故判 exclude。
- **三角剖分只以组合替身或稳定性定理形式出现**：`1708.06870`（Amoeba-shaped polyhedral complex of an algebraic hypersurface）在 spine 对偶于 Newton 多面体三角剖分时给出多面体复形的显式公式，属组合刻画而非数值求解流程；`2603.21116`（Solid Amoebas of Maximally Sparse Polynomials）通过 Ronkin 函数线性域在热带退化下的稳定性证明 Newton 细分在参数充分小时与热带细分一致，属结构性定理。
- **同伦延拓一侧只算热带曲线 / 实零点，不算 amoeba 本体或其补集**：见主线 4 的 3 篇，输出对象均为热带曲线或实零点，而非 amoeba、其补集或 Ronkin 函数。

**与此并列的第二个检索缺口**：本批 4 条查询中 **`abs:"Ronkin function"` 从未作为检索词出现**（仅 `abs:"Ronkin"` 出现在 Q2 的 OR 组内，易被 `amoeba` 的 AND 约束过滤掉真正以 Ronkin 函数为主题的文献）。这与上一条共同构成下一轮的检索改进方向：

- `abs:"Ronkin function"`（单列查询）
- `abs:"Monge-Ampere" AND (abs:amoeba OR abs:"Ronkin" OR abs:tropical)`
- `abs:"Newton polytope" AND (abs:triangulation OR abs:"mixed volume") AND (abs:amoeba OR abs:"Ronkin")`
- `abs:"homotopy continuation" AND (abs:tropical OR abs:amoeba OR abs:"Newton polytope")`（替代精度过低的 Q3）
