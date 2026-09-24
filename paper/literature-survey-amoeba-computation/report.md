# amoeba GBZ 计算方法调研报告

日期：2026-09-24 · 调研问题：**现有文献中，amoeba GBZ（及其等价的 Ronkin 函数 / 谱势判据）在算法层面是如何实际计算的？**

范围（用户确认）：① 非厄米 GBZ/amoeba 侧 ② 计算代数几何 + Toeplitz/Szegő 数学侧 ③ 公开代码与软件包
时间窗：不限，分 `recent3`（2023-09 之后）/ `classic` 两档
深度：核心 6 篇抓 arXiv 全文精读，其余标题/摘要级

> **本文件状态**：§1–§4 已完成；§5–§7 依赖 6 篇全文精读结果，待补。

---

## 1. 方法与流程

### 1.1 数据源

| 源 | 用途 | 日志 |
|---|---|---|
| arXiv API（`https://export.arxiv.org/api/query`） | 主检索：13 条查询式 | `query-log.md` |
| Semantic Scholar Graph API | 反向引用扫描（`arXiv:2212.11743` 的 100 条引用） | `query-log.md` §2 |
| web 检索 | 经典数学文献（多早于 arXiv）+ 软件/代码 | `query-log.md` §2 |
| 本地 PDF 库（`D:\information-base\library`，1713 篇） | 仓内已有资产清点（文件名级） | `sources-inventory.md` |
| 代码仓自带综述（`../literature-survey-2D-GBZ/`） | 上一轮成果复用 | `sources-inventory.md` |

### 1.2 PRISMA 流程

| 批次 | 主题 | retrieved | after dedup | 排除（标题/摘要级） | include | unclear |
|---|---|---|---|---|---|---|
| A | 非厄米 amoeba / Ronkin / 谱势 | 32 | 28 | 6 | 14 | 8 |
| B | 2D GBZ 数值算法 / Szegő–Toeplitz | 73 | **71** | 13（另有 29 条同名词异域条目单独索引） | **14** | **15** |
| C | 数学侧 amoeba 计算 | 200 | 151 | 117 | 24 | 10 |
| D | 主上下文补检（`ti:"lopsided" OR ti:amoebas OR abs:coamoebas`，totalResults=212，取回 40） | 40 | 去重后并入 §2 数学侧条目 | — | — | — |
| **三批合计（未做跨批去重）** | | **305** | **250** | **136** | **52** | **33** |

> **⚠ 计数口径不一致披露（诚实留痕，未自行裁定）**
>
> 上表各批数字取自**各批 agent 的自报值**。独立合并 agent 在生成 `screening-log.csv` 时做机械逐行点检，发现**自报值与文件实际内容不符**，三处如下：
>
> | 问题 | 自报值 | 机械点检值 | 处理 |
> |---|---|---|---|
> | 批次 B exclude | 13 | §2 主表逐行实为 **34** 个唯一 ID（全文另有 13/23/34/42 四种表述）；另有 **8 个 unique ID 无法从文件中定位**（71 − 63） | 取 34 条排除列表；8 个未定位 ID 不补 |
> | 批次 C include / unique | unique 151 / include **24** | `_batch-C-manifest.csv` 实有 **162** 行；§2 主表实为 **32** 行（**21 include** + 10 unclear + 1 exclude）；§3 逐条列出 exclude **128** 条 | 取 **21** include；排除按 `manifest − 主表` 的集合运算得 **132** 条 |
> | "已知跨批重复" 清单 | 我列了 14 个 | **仅 8 个真正跨批**：`2212.11743`(A,B)、`2407.01296`(A,B)、`2511.11349`(A,B,C)、`2608.28577`(A,B)、`2502.17931`(A,C)、`2212.06553`(A,C)、`math-ph/0311005`(A,C)、`2603.21116`(A,C)。另发现 **6 个未被列出的跨批关系**（批次 C 判 exclude 故不产生双行）：`2212.13704`、`1310.8472`、`1412.1585`、`1608.06077`、`2212.03173`、`math/0311062` | 以机械点检为准，`rationale` 列加 `[跨批判定不一致：X=…/Y=…]` 标注 |
> | 批次 B 首轮自报 | unique 40 / include 9 / unclear 23 / exclude 8 | 二次逐行点检：unique **71** / include **14** / unclear **15** | 取二次值（漏判含 `2407.01296`） |
>
> **权威交付物以 `screening-log.csv` 为准**：三批 include+unclear 去重后并集 **72 条**（**include 41**（含 5 条 `include (math-tool side)`）+ **unclear 31**）。上表"include 52 / unclear 33"是**自报值之和、未去重**，请勿直接引用。
>
> **对本轮核心结论的影响：无。** §6 的结论只依赖少数几篇关键文献（F6、`2511.11349`、`2407.01296`、`2608.28577`、`2609.23523` 及数学侧 5 篇），这些在两种口径下都稳定入选；差异全部发生在边缘条目。

**去重与复用**：三批 include+unclear 去重后 72 条；与上一轮综述 `../literature-survey-2D-GBZ/processed-papers.md`（64 个 ID）**重叠 11 个**：`2212.11743`、`2210.04412`、`2405.03750`、`2407.10166`、`2407.01296`、`2311.16868`、`2511.11349`、`2412.14912`、`2502.17931`、`2604.06998`、`2506.22743`（另 `2411.13661` 两轮均 exclude）。

主要排除原因分布：
- **同名不同物**（生物学黏菌、AMOEBA 极化力场、morphological amoeba 图像分割、图论 global amoeba、集合论 amoeba forcing）：约 30 条
- **数值代数几何应用类**与 amoeba 无关（视觉/统计/优化/弦论/机器人）：约 46 条
- **纯渐近界估计**（amoeba contour 次数上界等）：约 8 条
- **纯综述**：5 条
- **amoeba 仅作分析工具、无计算输出**：约 15 条

---

## 2. 五条算法主线（本次调研的核心结论）

把数学侧与物理侧合并后，"**实际计算 amoeba / Ronkin 函数**"的可执行方法只有五条主线：

### 主线 ① 边界 / contour + Gröbner 基

**思路**：不直接算 amoeba 本体，而是算它的边界（contour），再由边界围出 amoeba。

| 文献 | 贡献 | tier |
|---|---|---|
| **1310.7363** Schroeter, *The Boundary of Amoebas* (2013) | 定义并刻画 amoeba 的 extended boundary，区分 contour 与真实边界；**可在任意维通过边界计算超曲面 amoeba，二维仅需 Gröbner 基** | classic |
| **2608.03601** Nisse, *Singularities of Amoeba Contours* (2026) | 给出计算 amoeba contour 奇点的**显式实代数方程组**（平面曲线情形可实际求解，检测 nodes/cusps/多分支） | recent3 |
| 1708.06870 Nisse (2017) | Newton 多面体内部构造具备 amoeba 关键拓扑的**多面体复形**，spine 对偶于 Newton 三角剖分时给显式公式 | classic |

### 主线 ② lopsided 判据的工程化（文献量最大、最成熟）

**思路**：amoeba 的补集（lopsided 区域）可用**线性不等式**判定；把成员判定问题转化为可计算的代数/优化问题。

| 文献 | 贡献 | tier |
|---|---|---|
| **math/0603201** Purbhoo, *A Nullstellensatz for amoebas* (2006) | **理论总源头**：点不在 amoeba 中 ⟺ 定义理想中存在 lopsided 元素；由此给出用**线性不等式系统**逼近超曲面 amoeba 及其 spine 的方法 | classic |
| **1608.08663** Forsgård, Matusevich, Mehlhop, de Wolff, *Lopsided Approximation of Amoebas* (2017) | 把 Purbhoo 方法**工程化**：用 cyclic resultant 之间的关系解决其主要瓶颈；给出 **SINGULAR / SAGE 实现**，并量化相对通用 resultant 算法的加速 | classic |
| **1101.4114** Theobald & de Wolff, *Approximating amoebas and coamoebas by sums of squares* (2011) | 把成员判定化为**实代数可行性问题**，用实 Nullstellensatz + **SOS + 半定规划**给非包含证书；给出证书次数界与**实际计算实例** | classic |
| **1307.3681** Avendaño et al., *Metric Estimates and Membership Complexity for Archimedean Amoebae and Tropical Hypersurfaces* (2013) | 构造可高效计算的 **ArchTrop(f) 多面体逼近**，给 amoeba 与它的 **Hausdorff 距离显式上下界**；点成员判定固定维数下**多项式时间**（一维 amoeba 成员判定本身 **NP-hard**） | classic |
| 1711.02705 Forsgård & de Wolff (2017) | caissons 理论统一 lopsided amoeba 与指数和 amoeba，给出补集分量存在性的**改进证书** | classic |
| 1108.2456 Theobald (2011), 2603.21116 Nisse (2026) | 补分量存在性/实心性的判据；2603.21116 **证明 Passare–Rullgård 极大稀疏实心猜想**并给出 Newton 细分→热带细分的稳定化 | classic / recent3 |

### 主线 ③ 网格 / 可视化型直接算法

| 文献 | 贡献 | tier |
|---|---|---|
| **1604.03603** Bogdanov, *Algorithmic computation of polynomial amoebas* (2016) | 给出 amoeba、其 contour、**紧化 amoeba** 及三维 amoeba 二维截面的计算与可视化算法；可构造固定 Newton 多面体下"最复杂拓扑"的多项式；**Matlab 8 / Mathematica 9 实现** | classic |

### 主线 ④ 同伦延拓 / 数值代数几何

| 文献 | 贡献 | tier |
|---|---|---|
| **1408.3105** Jensen et al., *Computing tropical curves via homotopy continuation* (2014→Exp. Math. 2016) | **利用 amoeba 与热带曲线的联系**，用同伦延拓计算热带曲线**并给出实现** | classic |
| 1605.04203 Brake et al. (2016) | 同伦延拓 + monodromy loop + Cauchy 积分计算复/实热带曲线射线，**Bertini 实现** | classic |
| 1910.01957 Ergür et al. (2019) | 基于 Viro patchworking 的 **polyhedral homotopy** 求实零点，解位于 **A-discriminant amoeba 补集**的无界分量中 | classic |
| 0911.1783 Leykin (2009) | **Macaulay2** 的 NumericalAlgebraicGeometry 包（同伦延拓软件底座） | classic |

### 主线 ⑤ 物理侧：把谱势计算归约为 Ronkin 函数优化

| 文献 | 贡献 | tier |
|---|---|---|
| **2212.11743** Wang, Song, Wang, *Amoeba Formulation…Arbitrary Dimensions* (PRX 14, 021011, 2024) | 用 amoeba/Ronkin 重构非布洛赫能带论；**本项目 amoeba 模块的复刻对象** | classic |
| **2511.11349** Kaneshiro & Peters (PRR 8, 013292, 2026) | Wiener–Hopf 分解 + Hermitian doubling；**明确指出 amoeba 数值实现此前只限单带**，多带适用性连 1D 都不清楚 | recent3 |
| **2502.17931** Kaneshiro (2025) | 通过外推总 Ronkin 函数得 **band-resolved Ronkin 函数**，克服 Kramers 对导致的常规 amoeba 形式失效 | recent3 |
| **2407.01296** Xiong, Xing, Hu (→ Commun. Phys. 2026) | 谱势 + **几何自适应**非布洛赫能带论；非收敛/谱不稳定时趋向 amoeba 谱 | recent3 |
| 2405.03750 Yang & Bergholtz (PRR 7, 023233, 2025) | amoeba（separation gap）与 GBZ（surface gap）头对头比较 | recent3 |
| 2407.10166 Wang & Yan (PRB 110, L201104, 2024) | 基于 amoeba 表述的 infernal point（谱坍缩）判据 | recent3 |
| 2608.28577 Gu, Fu, Hu, Wang (2026) | **squeezed GBZ / squeezed amoeba** 表述：长程跳跃截断使 amoeba 被"压缩" | recent3 |
| 2607.22976 Xu, Wang, Deng, Yi (2026) | 把 Ronkin 函数形式推广到**域壁环形构型** | recent3 |
| math-ph/0311005 Kenyon, Okounkov, Sheffield (2003) | dimer 表面张力 = 谱曲线 Ronkin 函数的 **Legendre 对偶**；amoeba = 相图 | classic |
| 1306.6249 Lundqvist (2013) | **显式计算 Ronkin 函数**二阶导数（椭圆积分/超几何函数），半显式 Ronkin 测度 | classic |
| 2402.08798 Bobenko et al. (2024) | M-曲线上亚纯微分积分给 Ronkin 函数**显式公式** + Schottky 均匀化计算流程，数值与理论一致 | recent3 |
| 2104.04408 Arzhakova (2021，Israel J. Math. 2025) | decimation 极限 = −(Ronkin 函数的 Legendre 对偶)，d=2 时与 dimer 表面张力一致 | classic |

### 主线外的数据驱动/ML 分支（不算严格算法）

- 2212.06553 Chen et al. (2022)：遗传符号回归提取 amoeba 分量与 Mahler measure 的数值关系；ML 学习 3D amoeba 拓扑与 membership
- 2106.03695 Bao, He, Hirst (J. Symbolic Computation 116, 2023)：神经网路复现 lopsidedness 判据（F0 正系数约 99% 准确率）
- 2506.08618 Yan et al. (ICLR 2026)：**Poly2Graph** 流水线 + **HSG-12M** 数据集（11.6M 图 / 1401 个特征多项式类 / 源自 **177 TB** 谱势数据），GNN 基准 —— 目前规模最大的谱势数值流水线

---

## 3. 一个值得记下的"缺口"

> 批次 C 结论（原文）：*"真实 triangulation + Monge–Ampère 的『直接算法』在 arXiv 这一批里没有出现——Monge–Ampère 只以解析字典形式出现（0805.1194、1310.8472），这可能是你研究问题的真正 gap。"*

补充：**Passare–Rullgård（Duke Math. J. 2004）** 关于"amoeba 的补集由 Newton 多面体三角剖分给出"的经典结果，其正文**未能在本次检索中获取**（`web_fetch` 不支持 PDF；该文早于 arXiv）。它在你们的本地 PDF 库里存在。**这是本次调研的一处明确缺口，需要在报告定稿前补读。**

---

## 4. 非厄米侧的规模信号

`abs:amoeba AND (abs:"non-Hermitian" OR abs:"skin effect" OR abs:"Brillouin")` 的 arXiv 摘要命中 **totalResults = 8**；`abs:"spectral potential" AND (非厄米相关词)` 命中 **3**。

**含义**：把 amoeba 与非厄米物理关联起来的文献在 arXiv 摘要层面是两位数量级的小集合。这本身是 novelty 判断的重要证据。

---

## 5. 全文精读卡片

计划精读 6 篇。**注意：`web_fetch` 对 arXiv HTML 页有约 50 KB 截断上限，多篇精读只能覆盖前半部分；凡未取得的部分一律标注，不做推测。**

| # | 文献 | 主题 | 精读覆盖度 |
|---|---|---|---|
| 1 | 2212.11743 | F6 amoeba 表述（物理侧基准） | ✅ 绝大部分全文 |
| 2 | 2407.01296 | A9 谱势 + 几何自适应（现 Commun. Phys. 2026） | ✅ 全文（arXiv v1） |
| 3 | 2511.11349 | B4 WHF 多带 amoeba | ✅ 部分（§I–§II.3） |
| 4 | 1608.08663 | lopsided 近似 + SINGULAR/SAGE 实现 | ✅ 绝大部分全文 |
| 5 | 1101.4114 | SOS/SDP amoeba 成员判定 | ✅ 全文 |
| 6 | 2609.23523 | 普适 GBZ 理论 I：谱方法综述（2026-09-20） | ✅ 部分（§I–§II） |

---

### 卡片 3 · arXiv:2511.11349 — Kaneshiro & Peters, PRR **8**, 013292 (2026)

- **元数据**：Shin Kaneshiro, Robert Peters；2025-11 投稿 / 2026-03 发表；Phys. Rev. Research 8, 013292；DOI 10.1103/s43l-h6z6；理论（数学物理）型。15 图。
- **计算对象**：OBC 谱势 φ(E) = lim (1/N) ln|det T_N[E−h]|（T_N 为分块 Toeplitz 矩阵），ρ(E) = (1/2π)Δφ(E)。归约为 Ronkin 函数优化 φ(E) = min_μ R_σ(μ)，R_σ(μ) = ∫₀^{2π} dk/(2π) ln|det σ_μ(e^{ik})|，σ(β) = E − h(β)。
- **输入假设**：1D、有限跃迁范围 p,q、M 带；σ(β) 在 |β|=1 上可逆；**单带 M=1（class A）是 §II.2 的明确前提**。
- **算法步骤（单带，可复现）**：
  1. det σ(β) = C_E β^{−p} ∏_{j=1}^{p+q}[β − β_j(E)]，根按模递增排序；
  2. R_σ(μ) = ln|C_E| − pμ + Σ_j max(μ, μ_j)，其中 μ_j = ln|β_j|；
  3. R_σ 的最小值落在 [μ_p, μ_{p+1}]（即 GBZ 条件），该区间内 W[σ_μ] = 0，OBC 谱上区间收缩为一点，|μ| 即逆局域长度；
  4. φ(E) = min_μ R_σ(μ) = ln|C_E| + Σ_{j=p+1}^{p+q} μ_j。
  即"**根求值 + 模排序 + 取中间 p+1..p+q 个根模对数之和**"。
- **多带部分（§II.3，AII†，M=2）**：TRS† 约束使 det σ(β) = det σ(β^{−1})，根成 Kramers 对 (β_j, β_j^{−1})；R_σ(μ) = ln|C_E| − 2pμ + Σ_{j=1}^{2p}[max(μ,μ_j) + max(−μ,μ_j)]；引入对称分解 Ronkin 函数 R_σ = R_σ^{(+)} + R_σ^{(−)}（R_σ^{(+)}(μ) = R_σ^{(−)}(−μ)），并施加凸性/分段线性/整数量子化导数；用 Z₂ 不变量与 WHF partial index 给出适用判据。
- **精度控制/可证性**：摘要与引言声称对 AII† 类"provide a rigorous proof of the generalized Szegő limit theorem"、并把 WHF partial index 认定为精确判据。**§V–VII 的具体定理与误差陈述未取得，无法核实。**
- **复杂度/规模/验证基准/代码**：**均无法核实（抓取截断）**。正文含 DATA AVAILABILITY 小节但内容未取得。
- **与 Szegő–Toeplitz 路线的关系**：Szegő 原始形式用 WHF 处理拓扑平凡哈密顿量（Widom 1974/1976）；Alase et al. 2023 证明 WHF partial index 记边态数目并建立体-边对应；Basor et al. 2019 的 modified Szegő 定理处理拓扑非平凡情形。本文把广义 Szegő 定理重解释为 modified Szegő 定理经由 Hermitian doubling 的推论。
- **⚠ 本次调研最有价值的一条原文（三处独立表述，摘要/引言/Eq.(25) 后）**：
  > "However, while the generalized Szegő limit theorem is formally applicable in arbitrary dimensions, **its implementation is limited to single-band systems, and its applicability to multiband systems remains unclear even in one-dimensional systems**."
  > "In multiband settings, **the Ronkin function mixes contributions across bands**; when degenerate states with different localization lengths occur, these contributions compete during the optimization of the Ronkin function, **causing the standard Amoeba formulation to fail even in one dimension**."
  > "a straightforward mathematical decomposition is generally not feasible: **The mathematical properties of the Ronkin function are typically not preserved under bandwise separation**."
- **对"通用 2D 多带 amoeba 数值计算"的贡献**：**不提供任何 2D 算法**（标题与摘要限定 one-dimensional multiband；引言称高维非布洛赫理论"remains a significant challenge"）。它界定的是"**标准 amoeba 优化何时失效**"（1D 多带 + 对称类判据）。**这恰是 2D 多带数值方案 novelty 的空隙证据。**
- **已知但未读到的部分**（来自导航/TOC，可信）：§VI "Numerical verification"，含 VI.1 "Correction for non-vanishing partial indices in two-band class A systems"、VI.2 "WHF framework in systems with TRS†"（VI.2.1 κ=0,1；VI.2.2 κ=2）；附录 A 分析 κ=2 相在 TRS† 保持微扰下消失；附录 B 说明存在孤立边态时公式仍成立且 partial index 计数边态数。

---

### 卡片 6 · arXiv:2609.23523 — Xu, Hu, Yang, *Universal GBZ Theory I: Review of the Spectral Approach* (2026-09-20)

- **元数据**：Zeqi Xu（厦门大学）, Jiangping Hu（中科院物理所）, Zhesen Yang（厦门大学 / APCTP）；arXiv:2609.23523v1 [cond-mat.mes-hall]；16 页 10 图；仅 arXiv 预印本；综述 + 解析推导型。
- **计算对象**：2D 单带模型双开边界矩形上的**有限尺寸 OBC 能谱**（把 n_y 个横向格点并成超胞，约化为等效 1D n_y 带链）；以及每个本征值 E 在 β_x 中的特征方程根集（2n_y 个，按模递增）。**核心对象是 1D 的 μ_{m,GBZ}(k)，不是 2D GBZ 数值算法。**
- **关键定位（对本调研至关重要）**：§I.2 原文——"**A detailed review of existing 2D theories will be given in subsequent papers.**" 即本文是系列的 **1D 部分**；§VI.0.2 的 "Amoeba Criterion" 是 **1D GBZ 条件的 amoeba 表述**，不是 2D amoeba 数值法。
- **"1D GBZ 理论不完整"可确证的部分**：(a) 常规条件 |β_{np}|=|β_{np+1}| **只描述热力学极限体谱，本身不决定有限尺寸修正、也不决定拓扑边缘态的量子化**；(b) "Even in this limit, however, exceptions are known"——热力学极限下亦存在例外；(c) 全文默认以**满秩条件** det[T_{−p}]≠0、det[T_q]≠0 为前提，秩亏情形推到 §IV.2（非厄米 SSH）。
- **1D→2D 失效机制（可读部分）**：把 n_y 个横向格点并入超胞，**增大 n_y 视为 1D→2D 的连续过渡**；随 n_y 增大 OBC 谱由曲线演化为填充面积的图样。定义 Δ(E) = |β_{x,n_y p+1}(E)| − |β_{x,n_y p}(E)|（p=1）：**近虚轴三角形区域的本征值几乎满足 1D 条件，左右翼偏离较大，且偏离随横向宽度增大**。
- **算法/数值细节**：可读区域**无算法步骤、无伪代码、无算法框**。仅有：§III–IV 用 **Vieta 公式**做有限尺寸精确计算（正文未取得）；§II 给出 H(β) = Σ_{i=−p}^{q} T_i β^i、ChE f(E,β) = det[E·I_n − H(β)] = 0、满秩条件，并指出乘 β^{np} 后得次数 n_s ≡ n(p+q)、常数项非零的多项式，故每个 E 恰有 n_s 个有限非零根。**求根方法、系数组装、容差、等模根破并规则一概未给出。**
- **精度控制与算例**：仅一个可读算例（Fig.1）。L_x = 100 固定，n_y = 1,3,7,11；t_x = t_{−x} = i，t_y = t_{−y} = 1/4；图注仅称 "All computations use 16-digit numerical precision"。隐含实空间矩阵规模 L_x·n_y = 100/300/700/1100（**文中未明写**）；ChE 次数 n_s = 2n_y = 2/6/14/22。**无容差、无收敛判据、无颜色标度定义、无误差条、无 L_x 扫描。**
- **开放问题清单（逐字保留）**：
  1. "what determines μ_{m,GBZ}(k)?" 与 "how does the 1D GBZ condition in Eq. (4) break down as the system crosses over from 1D to 2D, and how can this breakdown be characterized?"
  2. §I.1 三条编号问题：(i) "For fixed n_y, does the spectrum converge to the 1D GBZ prediction as L_x→∞, or do deviations persist in this limit?"；(ii) "If convergence occurs, how large must L_x be for the thermodynamic prediction to become reliable?"；(iii) "**If the required size is too large for direct numerical verification, how can the theory be tested numerically?**"
  3. "Understanding when and why the conventional condition fails is a central goal of this series."
- **代码**：可读区域内**无任何代码/数据可得性声明、无仓库链接、无致谢**。arXiv 摘要页的 "Code, Data and Media" 仅为平台第三方 widget，非作者提供。
- **系列预告**：摘要称 Paper II 走**波函数方法**；§I.2 称 §VIII 将 outline **Papers II–IV**（§VIII 正文不可检索）。
- **对本调研的判断**：**目前尚无公认的 2D GBZ 理论**（原文："Although there is no well-accepted 2D GBZ theory, it is clear that the 1D GBZ condition in Eq. (4) cannot be generalized to higher dimensions"）。
- **战略提示**：Xu–Hu–Yang 系列正朝 2D 推进（同组前作含 2311.16868 渐近 GBZ）。Paper I 未占 2D 数值方案的位置，但 Paper II–IV 是**需要持续盯的潜在竞争者**。其 §I.1 三条开放问题恰好是"数值方法类论文"可以正面回答的，与本项目 PRB 路线 B 高度契合。

---

### 卡片 5 · arXiv:1101.4114 — Theobald & de Wolff, *Approximating amoebas and coamoebas by sums of squares*

- **元数据**：Thorsten Theobald, Timo de Wolff；arXiv v1 2011-01-21 / v3 2013-06-13；期刊版 **Math. Comp. 84(291):455–473 (2015)**，DOI 10.1090/S0025-5718-2014-02828-7；MSC 14P10, 14Q10, 90C22；类型：理论（实/凸代数几何）+ 数值实验验证；20 页 6 图。**精读覆盖度：全文。**
- **计算对象**：**成员判定问题**——给定 λ∈(0,∞)ⁿ 判定是否在 unlog amoeba 𝒰_I（log amoeba 𝒜_I）中；同法给出 coamoeba 版本。产出的是**"点不在 amoeba 中"的代数证书**：多项式恒等式 G + H + 1 = 0（G∈I′，H 为 SOS），等价于"−1 在 ℝ[X,Y]/I′ 中是 SOS"。另可推出该点所在**补连通分支的 order**。**全文只做非包含（补集）证书，没有"点在 amoeba 内"的证书。**
- **输入假设**：I = ⟨f₁,…,f_r⟩ ⊂ ℂ[Z₁,…,Z_n]，λ∈(0,∞)ⁿ；可用归一化把 amoeba 归到 λ=𝟏。§4 明确只处理**超曲面** amoeba（单个 f）。**不要求特定 Newton 多面体条件**（Newton 多面体只通过 lopsidedness 进入）。n 任意，但**所有算例 n = 2**。
- **算法步骤**：
  1. 实数化：f(X+iY) = f^re + i·f^im；构造 I′ ⊂ ℝ[X,Y]，生成元 {f_j^re, f_j^im} ∪ {X_k²+Y_k²−λ_k²}。
  2. 实 Nullstellensatz：λ∉𝒰_I ⟺ ∃G∈I′、H SOS 使 G+H+1=0；标准证书形式 Σp_j f_j^re + Σp'_j f_j^im + Σq_k(X_k²+Y_k²−λ_k²) + H + 1 = 0。
  3. monomial 途径：理想 I* 生成元 {f_i^re, f_i^im} ∪ {(m_ij^re)²+(m_ij^im)²−μ_ij²}，m_ij = Z^{α(i,j)}，μ_ij = λ^{α(i,j)}。
  4. 次数截断：C_t := {λ ∈ (0,∞)ⁿ\𝒰_I : 存在次数 ≤ 2t 的证书}。**Theorem 3.11**：C_t 单调递增且逐点收敛到 𝒰_I 的补集。
  5. 化为 SDP：H = MQMᵀ、Q ⪰ 0，M 为单项式向量；次数限制下的线性组合经"系数比较"并入半定规划。
  6. 实际求解（§5）：**SOSTools（Matlab）+ SeDuMi**；**对网格逐点解 SDP**，记录 feasible / infeasible / 数值不稳定。
- **精度控制/可证性（关键的"不对称性"）**：证书本身是严格的多项式恒等式（精确验证即证明），但计算是数值 SDP。**feasible ⇒ 严格证明点在补集；infeasible 只说明该次数内不存在证书**——原文："these bounds are lower bounds since feasibility of the SDP certifies membership in the complement of the amoeba but **infeasibility only certifies that no certificate with polynomials of degree at most k exists**"。
- **可证次数界**：Theorem 4.4（monomial 途径，lopsided 支配项条件）⇒ 存在总次数 **2·deg(f)** 的显式证书；Corollary 4.5(1) ⇒ 任意 w 存在次数 ≤ **2·rⁿ·deg(f)** 的显式证书（**含 rⁿ 因子**，源于 Purbhoo 迭代 resultant 构造 f̃_r，r→∞ 才逼近整个 amoeba）；Corollary 4.6 线性超平面 amoeba 补集中每点都有证书，且 SOS 是**仿射函数的平方和**。
- **复杂度与规模**：**文中未给出复杂度定理**（无 SDP 规模随 n/次数增长的界），只定性称可高效求解。实测（全部二维）：Ex 5.1 线性 f 用 **250×250** 网格、乘子次数 2；Ex 5.2 用 **160×160** 点、乘子次数 3（次数界 6）；Ex 5.3 每点 **14 个 SDP**（二分搜索）。**3 变量及以上、次数 ≥6 无算例。**
- **失效模式/限制（六条）**：(i) **数值不稳定严重**——作者自己的算例中大量网格点被 SeDuMi 报数值不稳定甚至程序中止，作者称其"deserves further study"；(ii) 可行性判定不对称（infeasible ≠ 在 amoeba 内）；(iii) 每网格点一个 SDP，成本随网格点数与半定矩阵尺寸快速增长；(iv) 可证强次数界只对 monomial+lopsidedness 情形成立且含 rⁿ 因子；(v) 开放问题：如何从**任意** Nullstellensatz 证书推出所在补分支的 order；(vi) §4 的显式证书只覆盖超曲面 amoeba。
- **代码**：**文中未给出**。§5 仅说明使用 SOSTools + SeDuMi；全文未见代码/脚本/数据公开或仓库链接。
- **与 lopsided 路线的关系**：论文只把 Purbhoo（lopsidedness，**不等式型**证书、由**迭代 resultant** 给出）当引擎与比较对象（**文中未提及 Forsgård 等人的工作**）。SOS 路线优势 = 代数证书（可精确验证的恒等式）+ 可读出补分支 order + 单点判定交给成熟 SDP；代价 = 需网格化、数值不稳定、可证覆盖范围理论上不超过对应 lopsided 区域。**纯 lopsidedness / Newton 多面体方法在二维只需比较各项模长（便宜、无 SDP），但没有多项式恒等式证书。**
- **迁移到非厄米 2D GBZ 的可行性（原文标注为推断）**：框架对口（GBZ 的核心对象正是二元 Laurent 特征多项式的 amoeba 及其补分支，算例全是 n=2、实现即逐点判成员的网格扫描），但有**三条硬约束**：(1) 只能严格证明"点在 amoeba 补集"，**画 GBZ 曲线要靠补集边界逼近 + 二分**，成本与数值风险高；(2) 无复杂度保证，实测只到 2 变量、次数 ≤6；(3) **无公开代码**。结论：**适合做二维 GBZ amoeba 成员判定的可行性证书/交叉校验工具，不适合当主力 GBZ 求解器。**

---

### 卡片 2 · arXiv:2407.01296 — Xiong, Xing, Hu（正式版 Xing, Xiong, Hu）, Commun. Phys. **9**(1), 2026, DOI 10.1038/s42005-026-02546-2

- **元数据**：arXiv v1 2024-07-01 → **Communications Physics 9(1), 2026-03-02**（gold OA）；类型：理论（数值用于验证）。**精读覆盖度：全文**（arXiv v1 的 §I–§X 与 Appendix A/B/C 逐节取到，含全部公式与图注规模数字）；正式版正文与 SI（MOESM1/2）因认证拦截/不可解析未取到，**两版权差异未核**。
- **计算对象**：复能量上的**谱势 Φ(E)**（作者已证 1D 情形等于 OBC 谱的静电势）；由 Poisson 方程给态密度 ρ(E) = ∇²Φ/2π；**GBZ 被定义为 Φ 的最小化参数集**；以及 σ_G（几何 G 下的非布洛赫谱）。
- **输入假设**：d 维 Bloch 哈密顿 H(k₁…k_d)（**多带 s>1 允许**）；"几何"**只用一组 d 个 lattice-cut 方向** {k₁…k_d} 描述（正则形状 = 2D 平行四边形 / 3D 平行六面体）；hopping range p_j,q_j 与 cut 长度 l_j；外加一条**收敛假设**（谱与本征态随尺寸稳定）。**不假设高维 Szegő 定理。**
- **算法步骤（可复现）**：
  1. 按几何做**基变换**（方↔菱：(k₁′,k₂′) = ((k_x+k_y)/2, (−k_x+k_y)/2)，β₁ = β₁′β₂′⁻¹, β₂ = β₁′β₂′），得解析延拓 H(β₁…β_d) 与 ChP f(β,E) = det[H(β)−E] = Σ_{j=−p}^{q} f_j β^j。
  2. 对第 j 个 cut 的**柱面**（该向 PBC、余向 OBC），固定 k_j 求 ChP 关于 β_j 的根并按模排序，**取模最大的 q_j 个**：ψ_{k_j}(E) = Σ_{n=p_j+1}^{p_j+q_j} log|β_{j,n}(k_j,E)| + log|f^{(j)}_{q_j}(k_j)|（Eq.20/22）。
  3. 对 μ_j 取极小：Φ̃_j(E) = min_{μ_j} ∫₀^{2π} (dk_j/2π) ψ_{k_j−iμ_j}(E)（Eq.21/26）。
  4. 合并所有 cut：**Φ(E) = min{Φ̃₁(E),…,Φ̃_d(E)}**（2D: Eq.23；dD: Eq.27）。
  5. **GBZ**：各方向独立取 argmin 得 μ_{j,min}（即逆局域长度 μ_j = log|β_j|），再由 f(e^{ik₁+μ_{1,min}},…,e^{ik_d+μ_{d,min}},E)=0 定出 k_{j,min}。
  6. DOS：ρ(E) = (1/2π)∇²_E Φ(E)。
  数学根据：1D 的 min_μ 与局部形式（Eq.13）在 Appendix A 被**严格证明**（分 s=p、s>p、s<p 三种情形；Taylor 展开后 k 相关项积分归零；等号仅当 |β_p| ≤ e^μ ≤ |β_{p+1}|）；**d≥2 是从 1D 结论层级式构造**。可分离模型有解析闭式（Appendix B）。
- **精度控制**：**主文未给出任何容差、网格步长或收敛判据**。min_μ 是 1D 连续优化；1D 下作者用 ChP 根（如 Ferrari 公式）解析求得，高维靠组合 1D 结果——**精度取决于根求解，不存在"2D 网格 μ"这一步**。与数值的一致性属"图目测 + 有限尺寸外推"级别，未报误差量级或阈值。
- **复杂度与规模**：未报计算代价或网格。实测最大规模：2D critical 模型 **N = 6400（方）/ 6385（菱）**；边界比依赖 N = 3267–3278；无序谱 N = 3600/3613（δ=0.2）；1D critical 两链 L = 80。适用范围 **O(L^d) 主导 skin mode**，不含 O(L^j)（j<d）次级/高阶模。
- **验证基准**：有限尺寸严格对角化（谱/DOS）与 Φ 的 Poisson 结果对照；1D 解析势 Eq.13；Appendix B 可分离模型的精确 DOS（ρ^{(y)}(E) = 1/(π√(6−E²))）；与 amoeba 谱作集合对照；1D critical 由 ED 本征态 Gauss 拟合宽度 1/κ，得 **κ ∝ 1/L 外推到 0**。
- **失效模式/限制（作者明说）**：唯一假设是谱收敛；一旦出现 **scale-free（critical）模**该假设失效 → 谱不收敛、依赖 boundary ratio、对弱扰动不稳定，且 lim_{δ→0}lim_{L→∞}φ ≠ lim_{L→∞}lim_{δ→0}φ（Eq.47，前者不可交换/未定义）；此时非布洛赫谱"**not well-defined**"，非布洛赫能带论与 GBZ 失效（Table 1 标 NA），无序把谱稳定到 amoeba 谱。另：**只处理正则几何**（不规则多边形/圆盘留待将来）。
- **代码**：arXiv v1 **只有 Acknowledgements（基金号），没有 Data/Code availability 段**；正式版 Data Availability 未取到 → **未确认公开**。
- **与 amoeba/Ronkin 路线的关系（关键）**：**不是替代而是包含**。可严格证明 (i) Φ_Amoeba 几何无关、(ii) σ_Amoeba = 消去所有方向 point-gap 的 uniform 谱、(iii) **Φ_Amoeba(E) ≥ 本文 Φ(E)**。关键：**E ∉ σ_Amoeba 时两者严格相等**；E ∈ σ_Amoeba 内本文几何相关势严格更小（"shift regions"），故 **σ_G ⊆ σ_Amoeba 且一般不等**。Conjecture-1：∪_G σ_G = σ_Amoeba，光滑边界 σ_smooth = σ_Amoeba。方法论差别：amoeba **假设** Ronkin 函数取极小，本文由 1D Szegő 定理**严格推出**后层级组合。
- **对"通用 2D 多带 amoeba 计算"的贡献**：给出一条**不依赖 Ronkin 网格/收敛技术**的替代路径——2D 多带 → "逐 cut 求 ChP 根 + 1D min_μ"，并由 μ_{j,min} 直接读出各方向逆局域长度；谱外与 amoeba 完全一致，谱内给出可解释的几何相关 DOS 修正。
- **⚠ 定位风险（需用户与理论文作者确认，本报告不自行裁定）**：A9 的"**逐 cut 柱面 + 沿该方向求根取模最大 q_j 个根**"与**本项目的 strip-GBZ（SGBZ）** 在构造上形式相似（都是"一方向开边界 / 其余周期 + 根模序"）。差别在于 A9 只给**各方向单一 μ_j 与 σ_G**，不给 GBZ 子集本身，且限于正则几何与主导 skin mode。**Introduction 必须显式写清与 Commun. Phys. 2026 per-cut potential 的差异**，否则审稿人极可能直接发问。
- **关键原文引用**：
  1. "Among all possible spectral deformations, **the one with the minimum spectral potential corresponds to the potential generated by the non-Bloch bands**."
  2. "the only assumption in our formulation is **the spectral convergence**… However, this assumption is not valid when the system exhibits scale-free localization"
  3. "it can be rigorously proven that our potential function in Eqs. (23) and (27) is **always no greater than the spectral potential in the Amoeba formulation**"

---

### 卡片 4 · arXiv:1608.08663 — Forsgård, Matusevich, Mehlhop, de Wolff, *Lopsided Approximation of Amoebas* (Math. Comp. **88** (2019) 485–500)

- **元数据**：Texas A&M；arXiv v1 2016-08-30 / v3 2017-09-09；期刊版 **Math. Comput. 88 (2019) 485–500**，DOI 10.1090/mcom/3323；cs.SC / math.AC / math.AG；类型 = 理论 + 工程实现；15 页 2 图 2 表。**精读覆盖度：绝大部分全文逐字**（§1、§2 至 §2.2、§3 全、§4 全含 Alg 4.1、§6 全含 Thm 6.1/Alg 6.2、§7 全含 Table 1/2 与图注）；**§5 正文未逐字取到**。
- **计算对象**：amoeba `A(f) = Log|Var(f)|`、unlog amoeba `U(f) = |Var(f)|`、lopsided amoeba `L(f) = {Log|v| : f 在该点不 lopsided}`；实际计算的是**循环结式** `CycRes(f;r) = ∏_{k₁..k_n} f(e^{2πik₁/r}z₁,…)`（仍是 Laurent 多项式）。证书语义 = **非成员证书**（成立 ⇒ Log|v| ∉ A(f)）。关键包含关系：**A(f) ⊆ L(f)，一般严格包含**（超集/外逼近）。
- **输入假设**：算法输入 `f ∈ Q[√−1][x]`（n 变量、次数 d）⇒ **系数必须是有理复数（精确算术）**；实现上用带参数 `I`（I²+1）的有理环模拟复数。网格 `G = [s,t]^n ∩ (ℓZ)^n`。**不要求非退化/无环/Newton 多面体特殊条件**；Newton 多面体只通过次数 d（误差界）与 `M = New(f)∩Z^n` 进入。§4 明确承认"没有 ε 该多小的显式表达式"。
- **算法步骤**：
  - **Purbhoo 骨架**：①据 ε 选层级 k（r=2^k），使 `rε ≥ (n−1)log r + log((n+3)2^{n+1}d)`；②造网格；③逐点用 `g_j := CycRes(f;2^j)`（j=1..k）测 lopsided，一旦 lopsided 即停并记层级，k 层都不 lopsided ⇒ 不认证；④认证点 ⇒ 不在 A(f)，并由主导项给 order；⑤判定工具 = 三角不等式（Lemma 2.6）；⑥收敛（Thm 2.3）：r→∞ 时 L(CycRes(f;r)) 一致收敛到 A(f)。
  - **循环结式替换一般结式（核心）**：Lemma 3.1 `Res_u(g(zu), u^r−1) = h(z^r)`；Lemma 3.3（倍增）`Res_u(g(zu),u^{2^{k+1}}−1) = Res_u(g(zu),u^{2^k}−1)·Res_u(g(zu),u^{2^k}+1)`。**Algorithm 3.4**：把 Multiplier 中第 j 个指数**不被 2^l 整除**的项全部变号，然后多项式相乘 —— 每次翻倍只做"**符号翻转 + 多项式乘法**"，**完全不算结式**。
  - **Algorithm 4.1**：逐网格点逐层测 lopsided，认证则记录 `α =`（主导项指数）`/ 2^j`；否则记为未认证。
  - **Algorithm 6.2**：给出 unlog amoeba 的半代数逼近（Thm 6.1）。
- **精度控制**：Remark 2.4 逐字给出 `rε ≥ (n−1)log r + log((n+3)2^{n+1}d)`；r=2^k 时 `2^k/k ≥ C(n,d)/ε`；N 只依赖 ε 与 Newton 多面体/次数且可显式计算。证书是**单侧（外部/超集方向）**；**没有任何"点在 amoeba 内"的证书**——§7.3 实践中"若到第 k 层仍未找到证书就**假定**该点在 amoeba 里"。§7.4 实测发现半代数逼近**不单调改善**（level 2 在某些区域优于 level 3），作者称非画图假象并"raises theoretical questions on the nature of the convergence"。
- **复杂度与规模**：**O(k d²)** vs 通用结式 signed subresultant **O(d·2^k)**；原文："reduces the runtime from **exponential to polynomial in one variable** and from **double to single exponential in arbitrary many variables**"。实测 Table 2 最快一行 **28.14 s vs 118397.86 s（≈4207×）**。代价侧：CycRes 度/项数随 r 指数增长——Table 1（f = z₁³+z₁z₂+z₂³+1）：r = 2/4/8/16/32/64 → 项数 10/31/109/409/1585/6241，次数 12/48/192/786/3072/12288，系数规模 9/860/>10¹²/>10⁵¹/>10²⁰⁴/**>10⁸¹¹**。
- **验证/算例**：Singular 计算 + Sage 作图；实现函数名 **`quickcyclicresultant`**。**实际规模：2 变量 3 次 4 项**多项式，网格 **81×81 点**、层数 **k=4**（Figure 1）；另有 3 变量算例 f3。**输出是网格点的外部判定 + 认证点的补集分支 order，不是边界曲线**；要轮廓需网格细分/边界追踪，或走 §6 的半代数路线。
- **失效模式/限制（七条）**：(i) **超集外逼近**，只给"在外"证书，未认证点只能假定在内；(ii) 精度渐近且**无 ε 显式下界**，k↑ ⇒ CycRes 指数膨胀；(iii) **要求精确有理系数**——浮点会毁掉不等式证书，**对非厄米参数是硬障碍**；(iv) §7.4 收敛不单调；(v) 复杂度严格分析只在 n=1（多元 = 逐变量迭代一元）；(vi) **Alg 4.1 第 25 行 ÷2^j 与 Thm 2.5 的 ÷r^n 不一致**（原文内部不一致，本报告如实记录、不擅自取舍）；(vii) §5 未逐字读，不能排除那里另有条件。
- **另有独立核出的原文数据问题**：Table 2 的 f2/Level 4 行 Factor 与比值不自洽（88.9828/0.3952 = 225.2 ≠ 301.42）⇒ 疑原文笔误。
- **代码**：论文 §1 给出入口 `http://www.math.tamu.edu/research/dewolff/LopsidedAmoebaApproximation/`；**Singular** 做结式与 lopsided 判定、**Sage** 画图；§7.4 点名函数 `quickcyclicresultant`。**核实结果：该链接现跨域重定向到 `https://artsci.tamu.edu`，旧路径未找到；archive.org 在本环境不可达 ⇒ 如今是否仍公开下载未能核实**；论文无 GitHub 等替代链接。
- **迁移到非厄米 2D GBZ 的可行性**：**对象完全对上**（2 变量 Laurent 多项式 f(E,β) = det[E−h(β)]），且最常需要的恰是"(μ₁,μ₂) = Log|β| 是否在 amoeba 内/外"的点判定，本方法给的正是这种判定、方向是**严格外部证书**；实现只需多项式乘法，无需根追踪/延拓。**障碍**：(1) 系数必须精确有理，E/跃迁是浮点或超越参数，需有理化或区间/高精度证书；(2) 只给"在外面"，要定位 GBZ 曲线/central hole 边界仍需网格二分并压小 ε，而 ε↓ ⇒ k↑ ⇒ 指数膨胀，**实际可用 ε 很有限**；(3) amoeba 内部点不判定，**必须与 Ronkin/绕数型内部判据互补而非替代**；(4) §7.4 的不单调收敛说明不能靠"加层数就更好"外推。

---

### 卡片 1 · arXiv:2212.11743 — Wang, Song, Wang, *Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions*, PRX **14**, 021011 (2024)

- **元数据**：Hong-Yi Wang, Fei Song, Zhong Wang（清华高等研究院）；arXiv 2022-12-22（v3 2024-04-30）→ **Phys. Rev. X 14, 021011 (2024)**，出版 2024-04-16，**CC-BY 4.0**，21 页 / 11 图 / 1 表；资助 NSFC 12125405。**类型：理论框架 + 数值算例——不是数值方法论文（无算法框、无代码）**。**精读覆盖度：全文绝大部分**（Abstract、§I、II、III、IV、V、VI、IX 逐节读原文）；缺 §VII（非布洛赫拓扑数值）、§VIII（谱不等式具体形式）、Appendix A/B 正文、Table I 内容。
- **计算对象**：同一套对象的四个视角 ——(1) **amoeba** 𝒜_f = {log|β| : f(β)=0} ⊂ ℝ^d，f = det[E−h(β)]；(2) **Ronkin 函数** R_f(μ) = ∫_{T^d}(dθ/2π)^d log|f(e^{μ+iθ})|；(3) **谱势** φ(E) = min_μ R_{det(E−h)}(μ)（**Ronkin 的全局极小**）；(4) **DOS** ρ(E) = Δφ(E)/2π。外加 GBZ 与能谱支撑集（洞闭合）。
- **输入假设**：一般 d 维平移不变（体态）非厄米紧束缚；h(e^{ik}) = Σ_n t_n e^{ik·n}，β = e^{μ+ik}，**矩阵值 Laurent 多项式**；有限跃迁范围；胞内带数形式上任意的 M。**算例：1D 非厄米 SSH（2 带，t₁=t₂=1, t₃=0.7, γ=4/3, L=300）；2D 单带模型（t=1, t′=0.5, γ=0.2，圆盘直径 L=140）；非厄米 Chern 带模型（2×2）。**
- **算法步骤（正文无编号算法框，据 Eq.14–41 与各节文字还原）**：
  1. 展开 det[E−h] = a_{−M}(E)β^{−M} + … + a_N(E)β^N，根按模排序 |β₁| ≤ … ≤ |β_{M+N}|。
  2. 对给定 E 构造 amoeba 与 Ronkin 函数；2D 零集是实 2 维，可用 β_x 局部参数化。
  3. **判谱靠中心洞——⚠ 这是"单向定理 + 猜想来补反向"，不是等价判据**：记 Λ = {E : f = det[E−h(β)] 的 amoeba **无**中心洞}。
     - **严格证明的部分（单向）**：**Eq.(37)：ρ(E) = 0, E ∉ Λ**，即 E ∉ Λ ⇒ 谱外（等价地 **谱 ⊆ Λ**）。原文接着说 "In other words, the bulk spectrum is **restricted inside** Λ"。Sec. V 已完整恢复，**通篇没有 "if and only if"、没有反向陈述**。
     - **未证的部分（反向）**：Sec. V 原文 "**if we assume the validity of the conjecture Eq. (34)** and therefore Eq. (24), we have ρ(E) = Δφ(E)/2π, which is generally nonzero in Λ" → "洞闭合 ⇒ 谱内"**建立在猜想 Eq.(34) 之上**（即 Φ(E) = φ(E) 的一般成立，由 Szegő 推广版给出）。而 Eq.(37) 的证明**只用原始版 Szegő 定理、不用推广版**——正说明严格部分恰好只覆盖单向。
     - **对冲措辞**：Sec. II.2 用 "could be an **indicator**"，作者自己写的是 "could be"。
     - 谱边界 = 洞闭合点（此步同样依赖该猜想）。
  4. 求 φ(E) = min_μ R_{det(E−h)}，ρ = Δφ/2π。结构依据：**R_f 在 amoeba 补集的每个 hole 上严格线性、到处凸，梯度 = 该 hole 的整数序 ν**；故中心洞存在时极小**落在中心洞（一片平台）**，否则极小**落在 amoeba 内一个单点**。
  5. **GBZ 定义**：det[E−h(β)] = 0 且 (log|β₁|,…,log|β_d|) = μ_min(E)。本征态 ψ_E(x) = Σ_k c_k e^{(ik+μ_min(E))·x}；实用上参数化为 μ_min(k)，该向量函数是 GBZ 的完整表示。
  6. 推导层用 Toeplitz + Szegő：(1/N)log|det 𝒯[E−h(β)]| = ∫(dk/2π)^d log|det[E−h(β)]| + O(L^{−1})，前提是 μ 取在中心洞内。
- **精度控制**：**唯一严格证明的是 Eq.37（ρ(E)=0, E∉Λ），只用了原始版 Szegő 定理、不用推广版**；Eq.32 有 O(L^{−1}) 余项。**作者自陈**（§IX 原文）："**not all aspects of this work are mathematically rigorous. Although numerical evidence is supplied whenever a mathematically strict derivation is unavailable…**"。§II.2 用的是"**indicator**"而非充要；§V 把谱边界等同于 Λ 的边界这一步**依赖 Eq.(34) 这个猜想**。**全文未给容差、收敛判据、误差界、极小化停止条件。**
- **复杂度与规模**：算过最大规模 = 1D SSH **L=300**；2D 单带圆盘 **L=140**；§VI 用 L=400 与变尺寸序列（E=1,3,5，丢弃距边界 20 格内的点）；§V 用 L∈[80,240] 线性外推（95% 置信区间），ΔE_t ∝ L^{−2}、ΔE_b ∝ L^{−1}。**未给** μ 网格分辨率、能量采样数、参数扫描规模、复杂度分析或与实空间对角化的耗时对比。
- **验证基准**：1D 逐点核对已知 GBZ 判据 |β_M| = |β_{M+1}|；OBC vs PBC 实空间对角化；2D 洞有/无两个能量点；DOS 对方格子（L=130，边界加 [−0.5,0.5] 均匀随机势抑制边界态）与圆盘（L=140）；带顶/带底洞闭合值 vs 实空间外推；**GBZ 的 μ_min vs 格林函数衰减率拟合**（作者称用 Green 函数比用本征态更准更省、可算更大 L）。
- **失效模式/限制（七条）**：
  1. 作者自陈非全部严格，部分核心结论只有数值证据。
  2. **高维 GBZ 不可直接计算**（d 维 GBZ 是 2d 维 β 空间中的 d 维子流形）。
  3. 1D 做法不能平移：2D 边界约束方程数 ∝ L，"it is difficult to obtain a 2D counterpart of Eq. (9)"，高维下依赖系数矩阵低秩"quite intractable"。
  4. **多带**：形式上是矩阵值符号，但**所有算例是单带或 2×2**，正文**无多带/简并处理方案**。（后续 `2511.11349` 点名该表述 "limited to single-band systems"，**作者本人未在本文自陈此点**。）
  5. **数值细节全缺——复现的最大缺口**：Ronkin 积分的求值策略（解析式 vs 数值求积、被积函数在 amoeba 上对数发散如何处理）、μ 极小化的搜索算法与容差，**全文未给**。
  6. 判据不含边界几何（后被称为 "neglects geometric information"）；本文主张的是"任意 generic 几何的 OBC 谱都收敛到同一普适 amoebic 谱"。
  7. Eq.39 脚注明说忽略边界项。
- **代码**：**在可读到的全文（含 §IX 正文与 Acknowledgements）中未发现任何代码/数据可得性声明**；致谢只有基金 NSFC No. 12125405。**注意措辞：这是"未找到"，不是"否认存在"**——缺失的 §VII/§VIII/Appendix A/B/Table I 未能取得（见 §6.4），故不能排除文末他处有此类声明。
- **未取得的小节（不影响本卡片核心方法字段）**：§VII（非布洛赫拓扑）、§VIII（谱不等式）、Appendix A/B、Table I。均已确认在本次工具环境下不可取得：microlink 分节代理日配额耗尽（每 IP 429），W3C html2txt / allorigins / codetabs / textise / cors.lol / corsproxy / r.jina.ai 全部 429/403/522/DNS 失败；URL fragment 与 query 参数不改变返回窗口；PRX PDF 被 fetch 工具以 content-type 拒绝；OpenAlex GROBID XML 需 API key（401）。**仅有的可用线索来自搜索引擎片段（不完整句，不可据以成卡）**：§VIII 围绕**一条** spectral inequality 且有 alternative proof；§VII 的非厄米 Chern 带模型为 2×2（参数 t, γ, v, m 均为实数）。
- **与 Ronkin 极小 / 平均绕数零点判据的关系（按第一手核实修正）**：
  - 结构性事实（本文给出）：平均绕数就是 Ronkin 的梯度 ν_j = ∂R_f/∂μ_j = Re∫(dθ/2π)^d (∂_{μ_j}f)/f = ∫_{T^{d−1}}(dθ/2π)^{d−1} w_j，w_j = (1/2πi)∮dθ_j ∂_{θ_j} log f；ν 在每个 hole 上取常整数序；**最多只有一个 hole 的序为 ν = (0,…,0)，即"中心洞"**；R_f 在补集每个 hole 上严格线性、处处凸，故极小或落在中心洞（平台）或落在 amoeba 内单点。
  - **但本文给出的严格结论只有单向**：「无中心洞 ⇒ 谱外」；「洞闭合/极小为单点 ⇒ 谱内」是**猜想 + 数值支持**。
  - **因此不能写成"等价判据"，也不能写成"与平均绕数零点判据等价"** ——后者是**结构性论证 + 后续文献转述**，**不是本文定理**。原文亦从未出现"等价"字样。
  - 1D 时退化为已知 GBZ 判据（区间收缩为单值 μ，两侧即 |β_M| = |β_{M+1}|）——这一退化是构造性说明，同样非本文证明。
- **图注原文（凸性 + 平台/单点二分，即"Ronkin 极小 ↔ 绕数零点"的结构依据）**：
  > "the Ronkin function is strictly linear on each component (hole) of the complement of the amoeba, where the gradient equals the integer index. The Ronkin function is always convex. Consequently, **when the central hole exists, the minimum is reached on the central hole; otherwise, the minimum is reached at a single point in the amoeba.**"
- **对本项目的意义**：F6 是把 amoeba 从几何图像变成**可计算判据**的那篇，也是本项目 `amoeba` 模块的复刻对象。但它是**理论框架文**：多带无方案、数值细节全缺、无代码。**这意味着你们的"数值实现"这一贡献是真实存在的缺口，而不是重复劳动。**

#### ⚠ F6 v1（2022-12-22）的两条关键原文——**同时是定位证据与最大风险**

v1 标题为 *"Amoeba formulation of the non-Hermitian skin effect in higher dimensions"*，其中写道：

> "the recent finding of geometry-dependent NHSE in 2D non-Hermitian systems suggests that **the GBZ might even not be definable in 2D**"

> "**A rigorous theorem is proved** as the general basis for DOS calculations"

**⚠ 出处标注（防误引，已核实）**：第二句**确为 arXiv:2212.11743v1 的原文**，但它是 **v1 Introduction 中对 Sec. IV 的预告句**，**不是 Sec. IV 的定理陈述**；同处紧接的第二句是 "We also demonstrate, despite the geometry-dependent NHSE, the existence of a universal spectrum to which the OBC spectrum under a generic geometry converges."。**v3 已把该句改写为** "We show in a theorem that the energy spectrum can be obtained from the shape of the amoeba." → 引用时必须标明是 **v1 Introduction 的预告**。
（曾有一条"该短语在公开索引中查不到、疑为未取证引用"的警告，**经第一手核实为误报**：公开索引并不收录 arXiv v1 的 HTML，查不到不构成否证。第一手记录见工作区 `paper/amoeba-retrieval/sec-IV-V.md` 第 131–136 行。）

- **作为证据**：连 amoeba 表述的原始作者都在 v1 里质疑"2D 中 GBZ 是否可定义"，说明本文 PRB 想解决的确实是公认的空白。
- **作为风险（必须正面处理）**：这条也是**指向你们论文的最强反论**——若审稿人属于 amoeba 阵营，会说"2D GBZ 根本未必可定义，你在算一个不well-defined 的东西"。**论文必须明确回答"你们的 2D GBZ 在什么意义下是良定义的"**：即它是 2d 维 β 空间中一个 **d 维子流形**（对给定几何/strip 构造而言），而 amoeba 是其上包络（σ_Amoeba ⊃ ∪σ_SGBZ）。这一节写不透，很可能是被拒的主要理由。

#### 一处用词提醒（防误引）

**v3 正文中没有"plateau ⇔ E 在谱外"这样的原话。** 本文原话是 Sec. II.2 的 "the absence (presence) of a hole … could be an **indicator**" + 定理 Eq.(37) + Fig. 3 图注的凸性/平台—单点二分。**"plateau（平台）"是后续文献（`2502.17931`、`2608.28577`）转述本文结论时的用词**，引用时不要把该措辞算到 F6 头上。

#### ⚠⚠ 与你们理论文核心主张直接冲突的一段（Sec. IV.2 原文，必须正面处理）

> "To summarize, there exists a **geometry-independent universal spectrum** that can be calculated from the amoeba and Ronkin function. By nature, it can be called the '**amoebic spectrum**.' The DOS of an OBC system with a **generic** shape always approaches the universal DOS in the large-size limit."
>
> "When the DOS of an OBC system with a certain (**nongeneric**) shape appears to deviate from the universal DOS, this deviation can be **eliminated by adding a small random local perturbation**."

- **这是 F6 的核心物理主张**：amoebic 谱是**几何无关的普适谱**，几何依赖的偏离只是"非普适形状"的有限尺寸效应，加一点随机扰动就消失。
- **而你们的 SGBZ 理论主张的是相反的图景**：不同 strip 几何给出**不同且都合法**的谱，σ_Amoeba ⊃ ∪σ_SGBZ 且一般**不相等**。
- **这意味着你论文里"验证 σ_Amoeba ⊃ ∪σ_SGBZ 且严格不等"（候选 C4）不再是可选的加分项，而是必须做的核心论证**——否则审稿人会直接用 F6 这句说"你们看到的只是非普适几何的假象"。
- ⚠ **注意口径**：F6 里的"随机性"**不是随机几何**，只是**边界格点上的随机在位势**（[−0.5, 0.5] 均匀分布）；论文里**没有随机几何的图**。引用时不要扩大化。

#### 精度与严格性的最终认定

- **Appendix A 被作者自己称为 "a heuristic proof of Szegő's limit theorem"**（启发式证明，非严格）。
- Sec. IX 自陈："not all aspects of this work are mathematically rigorous."
- 全文唯一可读到的定量界：Toeplitz 乘积估计 `‖𝒯[σ₁σ₂] − 𝒯[σ₁]𝒯[σ₂]‖₁ = O(L^{d−1})`（另有一条被截断的 `‖𝒯[σ]^{-1} − 𝒯[σ^{-1}]‖₁ ≤ …`）。

#### 代码可得性（确认版）

Acknowledgements **完整且只有一句**："This work is supported by NSFC under Grant No. 12125405." —— **无致谢名单、无 code-availability、无 data-availability 声明**；APS 元数据中亦无；`/supplemental/` 返回 **404**。→ **未公开代码、无补充材料。**

#### 两条新限制（作者明说）

1. **只处理标准 OBC** —— Sec. IX 原文："it seems that the amoeba approach, as we now understand, naturally corresponds to the standard OBC systems … we have focused on the OBC case throughout the present paper"（域壁等其它边界条件不在覆盖内）。
2. **例外点只是被"绕过"**：在 Eq.(44) 层面以脚注论证 nilpotent 项在围道积分后消失，**不是专门处理**。

#### 作者列的开放问题（Sec. IX 末段）

- **用于自由粒子极限的开放量子系统**：非厄米 Liouvillian 超算符谱决定动力学与弛豫，"immediately enables calculating the relevant quantities **beyond 1D**"。
- **多体非厄米系统**："our amoeba theory may still be a good starting point for including the interaction effects, which will be **left for future work**"。

#### F6 未取得小节的细节（如实标注）

- **§VIII** 的"可证 vs 猜想"状态**无法判定**：能确认的只有 Intro 原话 "several useful inequalities on the OBC and PBC spectra **are proved** from the amoeba approach"，以及该节围绕**一条** "the spectral inequality"（v1 该节标题为单数）并含一处 "alternative proof"；**任何一条不等式的完整式子都没拿到**。
- **Appendix A** 只见到首句（Szegő 定理源自 1D Hermitian Toeplitz 矩阵、后由 Widom 推广）；**Appendix B 是 v3 新增**（v1 无）。Table I 内容未取得。

---

### ⚠ 一条额外的重要发现（来自同一作者群的 2026 年后续工作）

F6 作者群（Zhong Wang 组）在 **2026 年的后续工作 `2608.28577`**（Gu, Fu, Hu, Wang，*How Long-Range Tails Reshape Non-Hermitian Spectra*）中明确写道：

> "**the higher-dimensional GBZ cannot be calculated directly. Instead, the OBC spectrum is encoded in an amoeba formulation**"

即：**amoeba 表述的原始作者本人，在 2026 年承认高维 GBZ 无法直接计算**，只能通过 amoeba 表述间接编码 OBC 谱。这条引文对本文 PRB 的定位段（"为什么需要一个可计算的 2D GBZ"）价值极高，且**来自最有分量的作者群**。建议在 Introduction 中与 `2511.11349` 的摘要句配合使用。

---

## 6. 结论

### 6.1 直接回答：大家是怎么算 amoeba GBZ 的？

把检索到的全部工作归拢，**"算 amoeba"这件事在文献里分成三拨互不通气的做法**：

| 拨 | 在算什么 | 做到哪一步 | 代表 |
|---|---|---|---|
| **数学侧（计算代数几何）** | 对一般 Laurent 多项式判定"某点是否在 amoeba 之外" | 只给**外部证书**；实际规模 2–3 变量、低次、小网格；要求精确有理系数或数值不稳定的 SDP | `math/0603201`、`1608.08663`、`1101.4114`、`1307.3681`、`1310.7363`、`1604.03603` |
| **物理侧（谱势/Ronkin）** | 对非厄米模型的**特征多项式**求谱势或 Ronkin 极小，从而定 GBZ 与能谱 | **只到 1D**：1D 单带（F6）→ 1D 双带 AII†（`2502.17931`）→ 1D 多带判据（`2511.11349`，且自承多带连 1D 都未解决） | 2212.11743、2502.17931、2511.11349 |
| **替代框架** | 绕开 amoeba 本体，走别的数学结构 | `2407.01296` 的 per-cut 柱面势（给各方向 μ_j 与 σ_G）；`2608.28577` 的 squeezed amoeba；`2607.22976` 的域壁 Ronkin | 2407.01296、2608.28577、2607.22976 |

**一句话**：数学侧的方法**通用但不实用**（只判"在外"、规模小、要精确系数），物理侧的方法**实用但只到 1D**。**"对一般 2D 多带模型算出 amoeba/strip GBZ"这个位置，在本轮检索范围内是空的。**

### 6.2 支撑这一判断的硬证据

1. **`2511.11349`（PRR 8, 013292 (2026)）原文**——三处独立表述：
   > "its implementation is limited to single-band systems, and its applicability to multiband systems remains unclear **even in one-dimensional systems**."
   > "**causing the standard Amoeba formulation to fail even in one dimension**."
   > "the mathematical properties of the Ronkin function are typically **not preserved under bandwise separation**."
2. **`2609.23523`（2026-09-20，本领域最新综述之一）原文**：
   > "**there is no well-accepted 2D GBZ theory**"
   > "A detailed review of existing 2D theories **will be given in subsequent papers**."
   且该文可读区域**没有任何算法/伪代码/求根方法/容差/代码声明**。
3. **检索规模信号**：`abs:amoeba ∩ (非厄米 OR 趋肤 OR 布里渊)` 的 arXiv 摘要命中 **totalResults = 8**；`abs:"spectral potential" ∩ 非厄米相关词` 命中 **3**。
4. **数学侧的三道硬墙**（逐条有原文/实测支撑）：
   - lopsided 路线**要求精确有理复数系数**（浮点会毁掉不等式证书），且**只给"在外"证书**，未认证点只能*假定*在内；实测只到 2 变量 3 次 4 项、81×81 网格，循环结式系数在 r=64 时 >10⁸¹¹（`1608.08663`）。
   - SOS/SDP 路线**数值不稳定严重**（作者自己的算例大量报错/中止），`infeasible` **不等于**点在 amoeba 内；实测只到 2 变量、次数 ≤6；**无公开代码**（`1101.4114`）。
   - **一维 amoeba 的成员判定本身就是 NP-hard**（`1307.3681`）。
5. **"Newton 三角剖分 + Monge–Ampère 的直接算法"在 arXiv 数学侧并未真正出现**（只以解析字典形式存在）。而经典源头 Passare–Rullgård（Duke Math. J. 2004）的正文**本轮未能读取**（见 §6.4 缺口）。
6. **amoeba 表述原作者本人承认高维 GBZ 不可直接计算**——`2608.28577`（Gu, Fu, Hu, Wang，Zhong Wang 组，2026）原文：
   > "**the higher-dimensional GBZ cannot be calculated directly. Instead, the OBC spectrum is encoded in an amoeba formulation**"
   这是定位段最有分量的一条引文（与第 1 条配合使用）。
7. **F6 本身没有给出可复现的数值方案**：21 页正文无算法框、无代码、无容差/收敛判据；所有算例是单带或 2×2；作者自陈 "not all aspects of this work are mathematically rigorous"；**Appendix A 被作者自己称为 "a heuristic proof of Szegő's limit theorem"**；代码/数据**确认未公开**（致谢只有一句基金号、APS 元数据无声明、`/supplemental/` 返回 404）。**"怎么把 amoeba 算出来"在 F6 里是空缺的**——这正是本项目数值工作的立足点。
8. **⚠ 引用 F6 判据时必须注意其逻辑方向**：F6 的"洞判据"是**单向定理 + 猜想来补反向**——严格证明的只有「无中心洞 ⇒ 谱外」（Eq.37，谱 ⊆ Λ，只用原始版 Szegő 定理），而「洞闭合 ⇒ 谱内」**建立在猜想 Eq.(34) 之上**，Sec. II.2 的措辞是对冲的 "could be an **indicator**"。**不可把 F6 的判据写成"等价判据"，也不可写成"与平均绕数零点判据等价"**（后者是结构性论证 + 后续文献转述，非本文定理）。这一点直接影响任何"以洞判据为收敛准则"的数值方案：它继承的是猜想，而非定理。

### 6.3 对 PRB 投稿的直接影响

**(a) C 位表述建议微调。** 现在 `../PRB-writing-plan.md` 里的 C1 是"把 2D GBZ 判定化归为根跟踪"，这个说法**没错但不够锋利**。建议改成带比较级的表述：

> 首个对**一般多带、任意跃迁范围**的 2D 紧束缚模型给出 **amoeba 与 strip 两类 GBZ 的可控精度数值构造**的方案；已有实现只覆盖 1D 单带/双带，而计算代数几何路线只能回答"点是否在 amoeba 之外"且规模受限。

**(b) Introduction 必须显式写的三段对比**（否则必被审稿人问）：

| 对手 | 你必须写的差异 |
|---|---|
| F6（PRX 14, 021011） | 它是**理论框架文**（21 页，无算法框、无代码、无容差/收敛判据）；**所有算例是单带或 2×2**，正文无多带/简并方案；数值细节（Ronkin 积分求值、被积函数在 amoeba 上的对数发散处理、μ 极小化搜索与容差）**全缺**——这是复现它的最大缺口。你的贡献是**通用多带数值实现 + 统一 amoeba/strip 引擎 + 拓扑应用**。 |
| A9 / Commun. Phys. **9**, 2026（per-cut 谱势） | **最危险的一篇**：它的"逐 cut 柱面 + 根模序"与你的 strip-GBZ 构造形式相似。差异：它只给**各方向单一 μ_j 与 σ_G**，不给 GBZ 子集本身；只处理**正则几何**；只覆盖**主导 skin mode**；且无容差/收敛判据。**建议先与理论文作者确认你们 SGBZ 与它的确切区别再定稿。** |
| `2511.11349` / `2502.17931` | 它们是**1D** 的多带判据；`2511.11349` 的摘要句可直接作为"空白声明"引用。 |

**(c) 建议直接引用进论文的两句"空白声明"**：即 §6.2 第 1 条的 `2511.11349` 摘要句与第 2 条的 `2609.23523` 的 "no well-accepted 2D GBZ theory"。

**(d) 审稿人可能的新增攻击点与预案**：

| 攻击 | 预案 |
|---|---|
| "Commun. Phys. 2026 已经给了 2D 的 per-cut 方法，你有何不同？" | 用 §6.3(b) 表格；强调对方给的是 σ_G 与各方向 μ_j，**不是 GBZ 曲线/环面**，且无精度控制 |
| "为什么不用成熟的 amoeba 计算几何算法（lopsided/SOS）？" | 用 §6.2 第 4 条：只给外部证书、要求精确有理系数、规模 2 变量低次、数值不稳定、NP-hard —— 并说明**你的路线与它们互补**（可作交叉校验） |
| "你的方法与 F6 相比只是工程实现？" | 强调三点：多带（F6 未覆盖多带）、统一 amoeba+strip 引擎、以及**在数值构造的 GBZ 上算拓扑不变量**（这是 F6 没做的事） |
| **「2D GBZ 根本未必可定义，你在算什么？」——来自 F6 v1 原文的最强反论** | 必须在正文单列一段回答：你们的 2D GBZ 是 **2d 维 β 空间中的 d 维子流形**（对给定几何/strip 构造良定义），而 amoeba 是其**上包络**（σ_Amoeba ⊃ ∪σ_SGBZ）。用你们 README 里已有的包含关系作为结构性论证，并给出数值验证（这正是候选 C4）。**这一节写不透极可能是被拒的主因。** |
| **「几何依赖谱只是非普适形状的假象」——来自 F6 Sec. IV.2 的普适谱主张** | F6 原文主张存在 **geometry-independent universal spectrum**（amoebic spectrum），generic 形状的 DOS 都收敛到它，非 generic 形状的偏离"加一点随机扰动就消失"。**你论文里"验证 σ_Amoeba ⊃ ∪σ_SGBZ 且严格不等"因此从加分项变成必做项**：必须证明不同 strip 几何的谱在 TDL 下**依然不同**（而非有限尺寸效应），并明确定义"generic vs nongeneric"在你的框架里对应什么。⚠ 引用 F6 时注意其"随机性"只是**边界随机在位势**，不是随机几何。 |

### 6.4 本轮调研的局限与缺口（诚实清单）

| 缺口 | 影响 | 补法 |
|---|---|---|
| **Passare–Rullgård, Duke Math. J. 2004** 正文未读（`web_fetch` 拒收 PDF，该文早于 arXiv） | "Newton 三角剖分 → amoeba 补集"的经典源头无法逐字引用 | 你们本地库已有该 PDF；请在会话外转成文本后放入工作区，或直接贴关键段落 |
| **Theobald, "Computing amoebas", Experiment. Math. 11 (2002)** 正文未读（同上） | 主线③的代表作只能据元数据登记 | 同上（需自行获取） |
| **`2609.23523` §III–§VIII** 未读（arXiv HTML 在 Eq.(12) 后固定截断，PDF/TeX 源被工具拒收） | 其"1D GBZ 理论不完整"的具体内容、§VII 限制与挑战清单、§VIII 开放问题、Paper II–IV 预告均缺失 | microlink 分节代理**可能**可取；本次该代理配额已用尽，建议次日重试 |
| **`2511.11349` §III–§VII 与附录** 未读（约 50 KB 截断） | 其数值验证细节（模型、格点数、与 ED 对照）、复杂度、DATA AVAILABILITY 内容未确认 | 同上，用 microlink 逐节取 |
| **F6（2212.11743）§VII / §VIII / Appendix A,B / Table I** 未取得（microlink 日配额耗尽，其余代理全 429/403/522） | 非布洛赫拓扑数值、谱不等式具体形式、Table I 内容缺失；**不影响"怎么算 amoeba/谱势/GBZ"的核心字段**（Abstract + §I/II/III/IV/V/VI/IX 均已逐节原文精读） | 次日重试 microlink；或用户本地转文本 |
| **APS 正式版与 arXiv 版的差异** 未核（A9、`2511.11349`、`2609.23523`） | 引用页码/编号可能需以正式版为准 | 引用时以正式版为准，差异不影响结论 |

#### ✅ 最便宜的补法：用你本地已有的 PDF（需要你在会话外跑一条命令）

本次调研的两个最大缺口，**都能用你机器上的本地 PDF 直接堵住**——我没有 shell，所以需要你执行：

```bash
# 1) F6 的四节（§VII Chern 数、§VIII 不等式、Appendix A/B、Table I）
pdftotext "D:\information-base\library\非厄米与量子开放\WangZhong@THU\2024_Physical Review X_Wang et al_Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions.pdf" f6.txt

# 2) Passare–Rullgård（Duke Math. J. 2004）——「Newton 三角剖分 → amoeba 补集」的经典源头
pdftotext "D:\information-base\library\非厄米与量子开放\2004_Duke Mathematical Journal_Passare_Rullgård_Amoebas, Monge-Ampère measures, and triangulations of the Newton polytope.pdf" passare-rullgard.txt
```

把生成的 `.txt` 放进 `paper/literature-survey-amoeba-computation/`（或任意工作区路径告诉我），我就能逐字精读并补齐卡片。
（无 `pdftotext` 时可换 `mutool draw -F txt`、`python -c "import pypdf;..."`，或直接在阅读器里另存为文本。）

**其余两条明日可选路径（无需你操作）**：(a) 换 egress IP 重跑 microlink selector（`%23S7`、`%23S8`、`[id="A1"]`）；(b) OpenAlex GROBID XML + 免费 API key（`content.openalex.org/works/W4394855923.grobid-xml`，无 key 时 401），可得到按节结构的完整纯文本。

### 6.5 待用户确认的边界判定（按规则不由我裁定）

1. **`1201.6401`（Rusek）**：给出算法 + Sage 实现，但对象是 **p-adic / 非阿基米德判别式 amoeba**，不是复 amoeba。批次 C 按 P 维度判为 *unclear*。若研究只关心复 amoeba，应降为 exclude。
2. **`2608.03601`（Nisse, Singularities of Amoeba Contours）**：给的是 **contour 奇点**的显式方程组，不是 amoeba 本体或其补集。批次 C 判 include（有可执行计算内容），但若严格遵守"输出 amoeba 本体或其补集"，应降为 unclear/exclude。
3. **批次 B 中 4 篇纯数学 Toeplitz 渐近**（math/0012200、math/0612487、math/1102.4131、math/0212215）：已按你确认的"含 Szegő–Toeplitz 数学工具侧"标为 `include (math-tool side)`、relevance 1–3。若你希望 P 严格限定非厄米物理模型，需改判为 exclude。
4. **A9（Commun. Phys. 2026）与你们 SGBZ 的定位重叠**（见 §6.3(b)）：这属于**理论层面的定位判断**，需要你与理论文作者确认，本报告不自行裁定。

---

## 附录：文件清单

| 文件 | 内容 |
|---|---|
| `sources-inventory.md` | 第一轮：仓库内已有文献清点 |
| `query-log.md` | 全部检索式、命中数、批次统计、抓取限制 |
| `_batch-A-raw-screening.md` / `_batch-A-screening.csv` | 批次 A 原始筛选表（28 条 unique） |
| `_batch-B-raw-screening.md` | 批次 B 原始筛选表（40 条 unique） |
| `_batch-C-raw-screening.md` / `_batch-C-manifest.csv` | 批次 C 原始筛选表（151 条 unique） |
| `report.md` | 本文件 |
| `screening-log.csv` | **权威合并筛选日志**：72 条（include 41 + unclear 31），17 列全字段 |
| `processed-papers.md` | 跨轮次复用的已处理清单（Included 41 / Unclear 31 / Excluded 按批列出；含与上一轮综述的 11 个重叠 ID） |
| `cards/1608.08663-lopsided-approximation-of-amoebas.md` | lopsided 近似论文的独立完整卡片 |
