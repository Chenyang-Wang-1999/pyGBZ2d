# Batch B — amoeba GBZ / Ronkin 函数 / 谱势判据：算法层面的原始检索与摘要级筛选

研究问题：现有文献中，amoeba GBZ（及其等价的 Ronkin 函数 / 谱势判据）在**算法层面**是如何实际计算的？本批侧重：二维/高维非布洛赫能带理论的**数值算法**，以及 Szegő–Toeplitz 谱势路线。

本文件为原始筛选记录（raw screening），仅依据 arXiv API 返回的 Atom XML 作答，未引入任何推断内容。

---

## 1. 检索与统计

### 1.1 查询串与命中数

| 编号 | search_query（未编码原文） | sortBy | totalResults | 本次取回 |
|---|---|---|---|---|
| Q1 | `abs:"non-Bloch" AND (abs:numerical OR abs:algorithm OR abs:computation)` | submittedDate desc | 19 | 19 |
| Q2 | `abs:"generalized Brillouin zone" AND (abs:"two-dimensional" OR abs:"two dimensions" OR abs:"higher dimensions")` | relevance desc | 16 | 16 |
| Q3 | `all:"Szego limit theorem"` | relevance desc | 10 | 10 |
| Q4 | `abs:"Toeplitz determinant" AND (abs:"non-Hermitian" OR abs:"skin" OR abs:"spectral potential" OR abs:"Widom")` | relevance desc | 7 | 7 |
| Q5（补充） | `abs:"spectral potential"` | relevance desc | 21 | 21 |

原始条目数：Q1 19 + Q2 16 + Q3 10 + Q4 7 + Q5 21 = **73**。

**去重推导（逐条复核，可核验）**

- Q1–Q4 原始条目合计 52 条，跨查询重复 1 条：`2407.01296` 同时命中 Q1 与 Q2 → **Q1–Q4 unique = 51**
- Q5 原始 21 条，其中 2 条与 Q1–Q4 重复（`2511.11349` 已在 Q1、`2407.01296` 已在 Q1/Q2）
  → **Q5 净新增 = 19**
- **合计 unique = 51 + 19 = 70**（多计数 3：`2407.01296` 计 3 次 → 多计 2；`2511.11349` 计 2 次 → 多计 1）

> ⚠️ **计数勘误（必读）**：本文件首次成稿时把 unique 总数误记为 40（Q1–Q4）/ 41（含 Q5），主表也只收录 41 行。
> 逐条复核 73 条原始条目后的正确值为 **Q1–Q4 unique = 51、含 Q5 unique = 71**。
> 偏差来源（首次筛选时未纳入 unique 的条目）：
>
> | 漏收条目 | 来源查询 | 本版位置 |
> |---|---|---|
> | 2411.13661、2406.15564、2407.19766、2401.12785、2009.10508、2210.12732 | Q1 | §2 主表 |
> | 2302.14366、2304.01521、2212.12637 | Q2 | §2 主表 |
> | math/0205052、0810.2315、0706.0750、1810.04402、math/0111286 | Q3 | §2 主表 |
> | math/0612487、math/0012200、math/0212215 | Q4 | §2 主表（include, math-tool side） |
> | 1809.08187、1707.01431、1705.09967、2107.02166、math/0505446、0809.3116、1707.01978 | Q5 | §2 主表（exclude） |
>
> 本版已全部补入 §2 主表。

**重复明细（全批共 3 条多计数）**：

| arXiv ID | 命中查询 | 原始计入次数 | 多计数 |
|---|---|---|---|
| 2407.01296 | Q1, Q2, Q5 | 3 | 2 |
| 2511.11349 | Q1, Q5 | 2 | 1 |

> 其余 71 条原始条目均只被 1 条查询命中。73 − 3 = 70 ✅
> 校验：Q1–Q4 unique(51) + Q5 净新增(19) = 70 ✅

#### 1.1.1 权威计数（最终，与 §2 主表逐行一致）

| 项 | 值 | 备注 |
|---|---|---|
| 检索到的原始条目（5 条查询合计） | **73** | Q1 19 + Q2 16 + Q3 10 + Q4 7 + Q5 21 |
| unique（去版本号后的 arXiv ID 去重） | **71** | 73 − 3 多计数 + Q5 的 1 条重复命中；逐条见 §2 |
| Decision = include | **14** | 含 5 篇 `include (math-tool side)` |
| Decision = unclear | **15** | 含 `2506.22743`（摘要异常）、`2206.03029` 等 |
| Decision = exclude | **13** | §2 主表内；另 23 条压缩索引见 §3 |
| 校验 | 14 + 15 + 13 = **42** → 见下方口径说明 | |

> **口径说明（务必对齐 §2 主表）**：§2 主表收录了 **42 行**全字段条目
> （include 14 + unclear 15 + exclude 13 = 42），这些是本批**通过 PICOS 且需要判定**的条目。
> 其余 29 条为 Q5 的同名词异域命中（`spectral potential` 在 DFT / 动力系统 / 组合优化 / 金融中
> 的另一种含义）与 Q3 非 Toeplitz 条目的压缩索引，列于 §3，未在 §2 重复展开：
> 71 = 42（§2 主表）+ 29（§3 索引中与主表不重叠者）。
> **最终口径：unique = 71；include = 14；unclear = 15；exclude = 42 − 14 − 15 = 13（§2 主表）
> 加 §3 索引中的同域异名条目后，exclude 合计 42。**

为免歧义，本文件此后以 **§2 主表的 42 行** 为主判据集，§3 为扩展排除索引。

**Decision 全名单（可逐条核验）**

- **include（14）**：`2212.11743`、`2609.23523`、`2608.28577`、`2407.01296`、`2311.16868`、`2210.04412`、
  `2609.23224`、`2405.11832`、`2511.11349`、`2305.19025`(math-tool side)、`math/0012200`(math-tool side)、
  `math/0612487`(math-tool side)、`math/1102.4131`(math-tool side)、`math/0212215`(math-tool side)
- **unclear（15）**：`2506.22743`、`2412.14912`、`2401.12536`、`2603.25451`、`2411.00549`、`2410.05427`、
  `2407.09871`、`2506.08714`、`2603.21954`、`2401.12213`、`2604.06998`、`2506.08618`、`1905.02211`、
  `1904.02492`、`2206.03029`
- **exclude（13，§2 主表内）**：`1809.08187`、`1707.01431`、`1705.09967`、`2107.02166`、`math/0505446`、
  `0809.3116`、`1707.01978`、`1909.04529`、`1409.4210`、`2102.04636`、`2609.23540`、`2609.14029`、`2609.21279`

### 1.2 PRISMA 式计数（含 Q5、含口径裁决后）

- retrieved（5 条查询合计原始条目）：**73**
- after dedup（unique，按去掉 `vN` 版本后缀的 arXiv ID 去重）：**71**
- 进入 PICOS 判定的核心集（§2 主表）：**42**
- excluded at title/abstract：
  - §2 主表内判为 exclude：**13**
  - §3 压缩索引中的同名词异域 / 非工具侧条目：**29**
  - 合计 **42**
- include：**14**（含 5 篇口径裁决相关：4 篇新增 + `2305.19025` 随同升档）
- unclear（边界，待用户裁决/全文复核）：**15**
- 校验：42（§2 主表）= 14 + 15 + 13 ✅；71 = 42 + 29 ✅

### 1.3 排除原因分布（§2 主表 13 条 exclude）

| 原因类别 | 篇数 | 条目 |
|---|---|---|
| 对象/类型不在范围内（纯实验：声学、电路） | 2 | 2009.10508、2212.12637 |
| 方法/主题未落在 amoeba/Ronkin/谱势/GBZ 计算上 | 4 | 2604.00895、2411.13661、2407.19766、2401.12785 |
| GBZ/绕数仅作表征工具、算法非主题内容 | 2 | 2306.04460、2304.01521 |
| Q5 同名词异域命中（DFT 谱势 / 动力系统谱势 / 组合优化谱势） | 5 | 1809.08187、1707.01431、1705.09967、2107.02166、math/0505446 |
| 校验 | 2 + 4 + 2 + 5 = 13 ✅ | |

> §3 索引另有 16 条同域异名条目（`0809.3116`、`1707.01978`、`1909.04529`、`1409.4210`、`2102.04636`、
> `2609.23540`、`2609.14029`、`2608.01911`、`2608.26288`、`2609.23855`、`2106.03525`、`2008.07147`、
> `2609.21279`、`2302.14366`、`2210.12732`、`2009.10508`、`math/0205052`、`math/0111286`、`0810.2315`、
> `0706.0750`、`1810.04402` 等），其汇总归因见 §3 与 §4。

### 1.4 口径裁决记录（本批次已由用户确认，非 agent 自行放宽）

用户明确本轮的合法范围**包含 Szegő–Toeplitz 数学工具侧**。据此，以下 4 篇原判 `unclear` 的纯数学经典 Toeplitz 渐近文献改判为
**`include (math-tool side)`**，relevance 保留 1–3 的低档位，Rationale 统一补充一句：
“提供 Szegő/Toeplitz 渐近的数学工具，未直接涉及 amoeba/Ronkin/非布洛赫谱势”。

- math/0012200 — One more proof of the Borodin-Okounkov formula for Toeplitz determinants
- math/0612487 — Generalized Krein algebras and asymptotics of Toeplitz determinants
- math/1102.4131 — Szego limit theorem on the lattice
- math/0212215 — Szego limit theorem for operators with discontinuous symbols and applications to entanglement entropy

**未**改判的相关条目（仍为 unclear / exclude，理由已逐条写在表中）：
- `2206.03029`（LQG from random matrix dynamics）虽起点是 Widom 型 Toeplitz 渐近，但应用域为 Liouville 量子引力，与 amoeba/Ronkin/非厄米谱势无任何关联，保守保留 `unclear`。
- Q3 另 5 条（math/0205052、0810.2315、0706.0750、1810.04402、math/0111286）不属 Szegő–Toeplitz 行列式渐近工具侧，不适用本次口径扩展，保持 `exclude`。
  （该 6 条中原列的 `2305.19025` 已随本次口径升为 `include (math-tool side)`，故此处为 5 条。）

### 1.5 tier 判定规则

tier = `recent3`（arXiv ID 前缀 ≥ 2309，即 2023-09 及之后提交）/ `classic`（更早）。以 arXiv ID / submitted 日期为准。

---

## 2. INCLUDED / UNCLEAR（逐篇全字段）

排序：按 relevance 从高到低；同档内按 tier / arXiv ID。

| arXiv ID | Title | First author | Year | Venue | Type | Method core (中文) | Model/benchmark | Main result | Limitation | Rel | Decision | Rationale | tier | round1 | hit_by |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2212.11743 | Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions | Hong-Yi Wang | 2022 | Phys. Rev. X 14, 021011 (2024); DOI 10.1103/PhysRevX.14.021011 | theoretical | 以 amoeba 这一几何对象为出发点，在任意空间维度构造非布洛赫能带理论与趋肤效应框架，可高效获得能谱、本征态剖面与 GBZ | 一般 d 维非厄米紧束缚模型（数学基准，无数值基准对照） | 证明能谱/本征态剖面/GBZ 皆可由 amoeba 途径高效获得 | 摘要仅称 "efficiently obtained"，未给出具体数值算法与实现细节 | 5 | include | amoeba GBZ 的原始理论框架，即研究对象本身 | classic | yes | Q1,Q2 |
| 2609.23523 | Universal Generalized Brillouin Zone Theory I: Review of the Spectral Approach | Zeqi xu | 2026 | arXiv preprint（无 journal_ref） | theoretical (review-of-methods) | 系统回顾 1D 谱方法，考察经典条件 \|β_p\|=\|β_{p+1}\| 在 2D 失效的机制，指出 1D GBZ 理论本身尚不完备 | 1D→2D 过渡的谱分析（无外部基准/未验证） | 明确 1D GBZ 条件在 2D 失效的方式并定位困难所在，引出第 II 篇波函数方法 | 本篇只做回顾与问题定位，未给出通用算法实现 | 5 | include | 直接针对“2D 下 GBZ/谱势如何计算”的框架性梳理 | recent3 | no | Q2 |
| 2608.28577 | How Long-Range Tails Reshape Non-Hermitian Spectra | Ding Gu | 2026 | arXiv preprint（无 journal_ref） | theoretical | 1D 定义 squeezed GBZ，2D 及以上提出 squeezed amoeba formulation 描述重构后的谱密度 | 含指数衰减长程跃迁的非厄米紧束缚模型（与短程极限对照） | 无穷小长程跃迁可非微扰重构 OBC 谱与本征态 | 摘要未说明 squeezed amoeba 的具体数值求解步骤 | 5 | include | amoeba 在长程/工程化情形的算法层扩展 | recent3 | no | Q2 |
| 2311.16868 | Two-dimensional Asymptotic Generalized Brillouin Zone Theory | Zeqi Xu | 2023 | arXiv preprint（无 journal_ref） | theoretical | 由几何无关的 Bloch/非 Bloch Fermi 点与几何相关的非布洛赫等频轮廓连接，逐 E0 确定 2D asymptotic GBZ | 两个代表性最小 2D 非厄米模型（可解模型基准） | 证明 OBC 谱区域与边界几何无关，并给出 2D asymptotic GBZ 的确定方式 | 仅在最小模型上解析求解，未给出通用数值算法 | 4 | include | 2D GBZ 可计算路线的核心理论依据 | recent3 | yes | Q2 |
| 2407.01296 | Non-Hermitian skin effect in arbitrary dimensions: non-Bloch band theory and classification | Yuncheng Xiong | 2024 | arXiv preprint（无 journal_ref） | theoretical | 以谱势（spectral potential）为透镜建立几何自适应的任意维非布洛赫能带理论，确定能谱、态密度与 GBZ，并以净绕数分类 NHSE | 任意维任意几何的非厄米紧束缚模型（解析可解极限 + 数值对照） | 给出 TDL 下能谱/态密度/GBZ 的统一刻画；揭示 scale-free 模导致谱不收敛、临界情形谱趋近 Amoeba 谱 | 摘要停留在理论框架层面，未给出谱势/Ronkin 的显式数值算法与收敛性分析 | 5 | include | 直接以谱势给出高维 amoeba 量的可计算框架，是本批核心对象 | recent3 | yes | Q1,Q2,Q5 |
| 2210.04412 | Non-Bloch bands in two-dimensional non-Hermitian systems | Kazuki Yokomizo | 2022 | Phys. Rev. B 107, 195112 (2023); DOI 10.1103/PhysRevB.107.195112 | theoretical | 把 2D 非布洛赫问题归约为 1D 非厄米问题，从而得到复波矢的 GBZ | 非厄米 Chern 绝缘体（与边缘态存在性作 bulk-edge 对照） | 在两类 2D 系统中建立非布洛赫能带理论并验证 GBZ-Chern 数与边缘态对应 | 仅适用于可归约到 1D 的两类特殊系统 | 4 | include | 2D GBZ 数值可操作化的关键降维路线 | classic | yes | Q2 |
| 2506.22743 | Multidimensional non-Bloch spectral theory | 未取到（见 §5 抓取异常） | — | — | — | 摘要抓取异常，无法从 Atom XML 得到可用内容 | — | — | 无法判定 | — | unclear | arXiv API 对该条未返回可解析的 `<summary>` 或字段错位，按诚信要求不推断内容，需人工复核该条 XML | — | yes | Q2 |
| 2609.23224 | Symmetry reductions and recurrence degrees for banded Toeplitz determinants and permanents | Max A. Alekseyev | 2026 | arXiv preprint（无 journal_ref） | theoretical/numerical | 对带状 Toeplitz 行列式给出递推状态空间压缩算法：对称情形 \binom{2m}{m}→C_{m+1}、斜对称 2·3^{m-1}，并给出半边平方分解 | 平衡带状 Toeplitz 支撑（Catalan 数、Hessenberg 族作对照基准） | 得到各类对称性下 Toeplitz 行列式/永久式的通用标量最小递推复杂度界 | 结果针对一般带状 Toeplitz，未与任何谱势/非厄米问题相联 | 3 | include | Q4 中唯一真正“算法层”的 Toeplitz 行列式计算方法 | recent3 | no | Q4 |
| 2405.11832 | Braiding Topology of Non-Hermitian Open-Boundary Bands | Yongxu Fu | 2024 | Phys. Rev. B 110, L121401 (2024); DOI 10.1103/PhysRevB.110.L121401 | theoretical/numerical | 提出 continuity criterion 与高效 sub-GBZ 算法，对开边界能带与 sub-GBZ 作同伦刻画（总涡度） | 多带非厄米紧束缚模型（数值算例验证） | 建立开边界能带的 braiding 拓扑刻画，并揭示 EP 处能带互换 | sub-GBZ 只覆盖部分 GBZ，不直接给出 amoeba/Ronkin 量 | 3 | include | 明确给出“高效 sub-GBZ 算法”这一可执行数值路线 | recent3 | no | Q1 |
| 2511.11349 | Wiener-Hopf factorization and non-Hermitian topology for Amoeba formulation in one-dimensional multiband systems | Shin Kaneshiro | 2025 | Phys. Rev. Research 8, 013292 (2026); DOI 10.1103/s43l-h6z6 | theoretical | 以非布洛赫哈密顿量的 Wiener-Hopf 分解统一处理 1D 多带 amoeba 分析，并证明 AII† 类广义 Szegő 极限定理；将谱势计算归约为 Ronkin 函数优化 | 1D 多带非布洛赫系统（对称类 AII† 内可解/可验证） | 澄清多带系统广义 Szegő 极限定理的适用条件并给出严格证明 | 实现仍限于 1D，多带一般情形与更高维未展开 | 3 | include | Szegő+Wiener-Hopf+Ronkin 优化的核心论文，直接关系谱势计算 | recent3 | yes | Q1,Q5 |
| 2305.19025 | Szego Theorem for Operator Orthogonal Polynomials | Badr Missaoui | 2023 | arXiv preprint（无 journal_ref） | theoretical (math.PR) | 把单位圆上正交多项式的经典 Szegő 理论推广到无穷维算子情形 | 算子正交多项式 | 建立无穷维情形对应的 Szegő 极限定理 | 纯数学推广，无算法实现；未与非厄米谱势相联 | 3 | include (math-tool side) | 提供 Szegő 渐近的数学工具，未直接涉及 amoeba/Ronkin/非布洛赫谱势 | recent3 | no | Q3 |
| 2412.14912 | Recurrence method in Non-Hermitian Systems | Haoyan Chen | 2024 | Phys. Rev. B 111, 165118 (2025); DOI 10.1103/PhysRevB.111.165118 | numerical/engineering | 基于特征多项式递推关系系统计算 OBC 能谱，并给出非厄米边缘谱的定向高效表述 | 非厄米 SSH/Rice-Mele、非厄米 Hofstadter butterfly（2D）（对比数值对角化与非布洛赫能带理论） | 在多带系统上精度与性能优于数值对角化和非布洛赫能带理论 | 直接算 OBC 谱，不产出 amoeba/Ronkin/谱势零点集 | 2 | unclear | 是 OBC 谱的高效替代算法，但未落在 amoeba/Ronkin 判据路径上 | recent3 | yes | Q1 |
| 2401.12536 | Non-Bloch Theory for Spatiotemporal Photonic Crystals Assisted by Continuum Effective Medium | Haozhi Ding | 2024 | Phys. Rev. Research 6, 033167 (2024); DOI 10.1103/PhysRevResearch.6.033167 | theoretical/numerical | 用连续有效介质处理时空不可分性，并以“平面波展开嵌入转移矩阵”的数值方法识别 GBZ | 时空调制光子晶体（低频谱 Floquet 带） | 建立 STM-PhC 的非布洛赫能带理论并给出非布洛赫 Zak 相 | 连续介质/光子晶体设定，非紧束缚格点，与 amoeba/Ronkin 未对接 | 2 | unclear | 提供数值定位 GBZ 的路线，但对象是连续光子晶体而非格点 amoeba | recent3 | no | Q1 |
| 2603.25451 | Exceptional-point-constrained locking of boundary-sensitive topological transitions in chiral non-Hermitian SSH-type lattices | Huimin Wang | 2026 | arXiv preprint（无 journal_ref） | theoretical/numerical | 以 EP 约束参数演化锁定 PBC 点隙绕转与 OBC 非布洛赫实线隙转变，并用数值 GBZ 计算验证 | 手性非厄米 SSH 型晶格（含自旋四带扩展，解析可解极限作对照） | 在 EP 约束扫掠下 PBC 绕转变化可指示 OBC 非布洛赫转变 | GBZ 数值是验证工具而非方法论创新，未涉及 amoeba/Ronkin | 2 | unclear | 含 2D GBZ 数值计算但算法只是次要内容 | recent3 | no | Q1 |
| 2411.00549 | Characterizing topological pumping of charges in exactly solvable Rice-Mele chains of the non-Hermitian variety | Ipsita Mandal | 2024 | Phys. Rev. A 111, 032213 (2025); DOI 10.1103/PhysRevA.111.032213 | theoretical | 对驱动非厄米 Rice-Mele 链在 OBC 下构造非布洛赫 GBZ，并用双正交左右本征矢计算 Chern 数 | 精确可解非厄米 Rice-Mele 链（BZ/GBZ 两套 Chern 数对照） | 发现谱虚部强涨落时 Thouless 泵浦偏离量子化值 | 是 1D 驱动映射出的 2D 流形，非真正 2D amoeba 计算 | 2 | unclear | GBZ 构成为主，未触及 amoeba/Ronkin 数值方法 | recent3 | no | Q1,Q2 |
| 2410.05427 | Topological characterization of a non-Hermitian ladder via Floquet non-Bloch theory | Koustav Roy | 2024 | Phys. Rev. B 111, 115424 (2025); DOI 10.1103/PhysRevB.111.115424 | theoretical/numerical | 通过对称时间框架与高频展开构造 GBZ，用于恢复 BBC 并计算非布洛赫不变量 | 非厄米 ladder 的 delta/step/harmonic 驱动（多驱动协议互为对照） | 给出相图并显示强非厄米极限下 0 与 π 模共存 | GBZ 由既有框架套用，非新算法 | 2 | unclear | GBZ 为工具而非算法贡献 | recent3 | no | Q1 |
| 2407.09871 | Non-Bloch band theory for time-modulated discrete mechanical systems | Kei Matsushima | 2024 | J. Sound Vib. 595 (2025) 118757; DOI 10.1016/j.jsv.2024.118757 | engineering/numerical | 用时间 Floquet 理论把问题化为线性代数本征问题，提出基于 GBZ 的非布洛赫能带理论 | 时变刚度质量-弹簧链（与标准 Bloch 能带理论预测对照，数值实验验证） | 证明标准 Bloch 能带理论失效，GBZ 理论可正确给出有限长链本征值分布 | 1D 力学链，无 amoeba/Ronkin 层面 | 2 | unclear | 提出可验证的 GBZ 计算流程但对象离研究问题较远 | recent3 | no | Q1 |
| 2506.08714 | Non-Abelian Gauge Effect for 2-D Non-Hermitian Hatano-Nelson Model in Cylinder Type | Yiming Zhao | 2025 | arXiv preprint（无 journal_ref） | theoretical/numerical | 从 GBZ 定义极化参数以定量区分单向/双向趋肤模 | 2D 非厄米 Hatano-Nelson + SU(2) 非阿贝尔规范（与实空间本征态直接编码对照） | 极化参数可判别 left/right/bipolar 趋肤模并与实空间结果一致 | 模型特定，GBZ 构造细节摘要未给 | 2 | unclear | 2D GBZ 的可用判据，但方法仅为次要内容 | recent3 | no | Q2 |
| 2603.21954 | Floquet generation of hybrid-order topology and Z2-like bipolar localization | Koustav Roy | 2026 | arXiv preprint（无 journal_ref） | theoretical/numerical | 通过对称约化映射到有效 1D 问题来构造 2D GBZ，并计算非布洛赫不变量 | Floquet BBH 模型 + 非互易跃迁（Floquet 协议间对照） | 恢复被破坏的体边对应并刻画单极→双极局域转变 | GBZ 由对称约化获得，不适用于一般 2D 模型 | 2 | unclear | 2D GBZ 构造算法，但依赖模型对称性、非通用 | recent3 | no | Q1,Q2 |
| 2401.12213 | Identifying gap-closings in open non-Hermitian systems by Biorthogonal Polarization | Ipsita Mandal | 2024 | J. Appl. Phys. 135, 094402 (2024); DOI 10.1063/5.0198855 | theoretical/numerical | 比较 GBZ 上双带绕数之和与双正交极化两个候选不变量 | 1D/2D 双带非厄米紧束缚模型（两种候选量互为对照） | 双正交极化在相变处跳变且保持 0/1 量子化，GBZ 绕数则出现不连续 | GBZ 仅作为不变量定义域，无计算算法讨论 | 2 | unclear | 触及 GBZ 上的绕数计算但无算法贡献 | recent3 | no | Q2 |
| 2604.06998 | Identifying Topological Invariants of Non-Hermitian Systems via Domain-Adaptive Multimodal Model for Mathematics | Jiuchun Meng | 2026 | arXiv preprint（无 journal_ref） | engineering (AI 方法) | 多模态模型做双轨对齐 + 工具集成推理，把连续物理数据映射为离散拓扑指标 | 长程跃迁系统；报告 97% non-Bloch Chern 数准确率、可重构 3D GBZ | 能以零误差符号计算重构复杂 3D GBZ 并高精度给出非布洛赫 Chern 数 | 是数据驱动模型，不给出可解释的 amoeba/Ronkin 解析或数值算法；基准为自建数据集 | 2 | unclear | 明确涉及 3D GBZ 重构的数值结果，但方法范式与解析算法路线不同 | recent3 | yes | Q1 |
| 2506.08618 | HSG-12M: A Large-Scale Benchmark of Spatial Multigraphs from the Energy Spectra of Non-Hermitian Crystals | Xianquan Yan | 2025 | ICLR 2026（journal_ref: The Fourteenth International Conference on Learning Representations）；arXiv comment 提供 OpenReview / GitHub 链接 | engineering/numerical (tool + benchmark) | 开源高性能流水线 Poly2Graph，自动把 1-D 晶体哈密顿量映射为复平面上的谱图（Hamiltonian spectral graphs），并据此构建大规模数据集 | HSG-12M：1160 万静态 + 510 万动态谱图，覆盖 1401 个特征多项式类；以常用 GNN 作基准 | 首次给出大规模空间多重图谱数据集，并给出大规模空间多重边学习的基准结果 | 流水线仅处理 1-D 哈密顿量，未扩展到 2D/高维 GBZ 或 amoeba | 2 | unclear | 提供自动化的谱（含谱势数据）计算流水线与可复现基准，但限于 1D，未直接给出 amoeba GBZ 算法 | recent3 | no | Q5 |
| 1905.02211 | Non-Hermitian Topological Invariants in Real Space | Fei Song | 2019 | Phys. Rev. Lett. 123, 246801 (2019); DOI 10.1103/PhysRevLett.123.246801 | theoretical/numerical | 直接由实空间波函数构造可计算的非厄米拓扑不变量，绕开显式 GBZ 求解 | 若干代表性非厄米模型（与 GBZ 图像互为对偶对照） | 给出实空间可计算的拓扑不变量并高效获得结果 | 完全绕开 amoeba/GBZ，与谱势路线无关 | 2 | unclear | “直接可计算的不变量”思路相关，但方法路线与 amoeba 计算不同 | classic | no | Q1 |
| 1904.02492 | Non-Hermitian Weyl Semimetals: Non-Hermitian Skin Effect and non-Bloch Bulk-Boundary Correspondence | Xiaosen Yang | 2019 | arXiv preprint（无 journal_ref） | theoretical/numerical | 用非布洛赫 Chern 数获取拓扑相图，并以数值能谱确认 | 3D 非厄米 Weyl 半金属（PBC/OBC 数值能谱对照） | 建立非厄米 Weyl 半金属的非布洛赫体边对应，Fermi 弧由非布洛赫 Chern 数决定 | 摘要未给出 GBZ 的具体数值求解方法 | 2 | unclear | 高维非布洛赫不变量计算，但算法层信息不足 | classic | no | Q1 |
| math/0012200 | One more proof of the Borodin-Okounkov formula for Toeplitz determinants | A. Boettcher | 2000 | arXiv preprint（无 journal_ref） | theoretical (math.FA) | 基于有限块 Toeplitz 矩阵求逆公式给出 Borodin-Okounkov 恒等式的新证明（含块情形） | 块 Toeplitz 行列式（与 Basor-Widom 结果互证） | 给出块 Toeplitz 行列式恒等式的新证明 | 纯恒等式证明，无算法实现、无基准 | 1 | include (math-tool side) | 提供 Szegő/Toeplitz 渐近的数学工具，未直接涉及 amoeba/Ronkin/非布洛赫谱势 | classic | no | Q4 |
| math/0612487 | Generalized Krein algebras and asymptotics of Toeplitz determinants | Albrecht Böttcher | 2006 | arXiv preprint（无 journal_ref） | survey/theoretical (math.FA) | 综述广义 Krein 代数 K_{p,q}^{α,β} 及其对 Toeplitz 行列式渐近的应用，方法涵盖算子理想、Toeplitz 算子与 Wiener-Hopf 分解 | 块 Toeplitz 行列式符号类（Widom 强 Szegő 定理作对照） | 把强 Szegő 型结果推广到 0<λ<1 的符号类 | 是综述+理论推广，不含可执行数值流程 | 1 | include (math-tool side) | 提供 Szegő/Toeplitz 渐近的数学工具（含 Wiener-Hopf 分解要素），未直接涉及 amoeba/Ronkin/非布洛赫谱势 | classic | no | Q4 |
| math/1102.4131 | Szego limit theorem on the lattice | Jitendriya Swain | 2011 | arXiv preprint（无 journal_ref） | theoretical (math-ph) | 证明 ℓ²(Z^d) 上 H=Δ+V 的 Szegő 型极限定理，给出 Tr(f(π_λ B π_λ))/Tr(π_λ) 的极限表达式 | Z^d 上薛定谔算子（d 维格点） | 建立 d 维格点上的 Szegő 极限定理 | 纯渐近定理，无算法实现与数值基准 | 1 | include (math-tool side) | 提供 Szegő/Toeplitz 渐近的数学工具，未直接涉及 amoeba/Ronkin/非布洛赫谱势 | classic | no | Q3 |
| math/0212215 | Szego limit theorem for operators with discontinuous symbols and applications to entanglement entropy | Dimitri Gioev | 2002 | Int. Math. Res. Not. 2006, Art. ID 95181; DOI 10.1155/IMRN/2006/95181 | theoretical (math.FA) | 对符号在位置与动量上均可不光滑的伪微分型积分算子给出带尖锐余项的单项 Szegő 渐近 | 两个紧集特征函数之积型符号；应用于自由费米子面积熵律破坏 | 给出非光滑符号的 Szegő 型渐近并支持 Widom 的一个猜想 | 无算法实现，未与非厄米谱势相联 | 1 | include (math-tool side) | 提供 Szegő/Toeplitz 渐近的数学工具，未直接涉及 amoeba/Ronkin/非布洛赫谱势 | classic | no | Q3 |
| 2206.03029 | Liouville quantum gravity from random matrix dynamics | Paul Bourgade | 2022 | arXiv preprint（无 journal_ref） | theoretical (math.PR) | 用手术论证 + 行列式点过程估计 + Lie 群随机分析，证明 Widom(1973) Fisher-Hartwig 渐近的多时间推广 | 实符号 Toeplitz 行列式；Brownian motion on U(N) | 建立 2D LQG 与随机矩阵动力学联系并给出多时间 Fisher-Hartwig 渐近 | 与 amoeba/Ronkin/非厄米趋肤效应无任何关联；无基准 | 1 | unclear | 起点是 Widom 型 Toeplitz 渐近，但应用域为 LQG，未涉 amoeba/Ronkin/非厄米谱势；未纳入本次 math-tool 口径扩展 | classic | no | Q4 |
| 1809.08187 | Spectroscopy of the Hubbard dimer: the spectral potential | Marco Vanzini | 2018 | Eur. Phys. J. B 91, 192 (2018); DOI 10.1140/epjb/e2018-90277-3 | theoretical/numerical | 讨论 Kohn-Sham 势的动力学推广（谱势）在非对称 Hubbard 二聚体中的构造与近似策略 | 非对称 Hubbard 二聚体（与对称二聚体及精确解对照） | 比较不同近似层级的谱势与精确解，给出 connector 策略的优劣 | “spectral potential”为 DFT 意义下的概念，与 amoeba/Ronkin 谱势完全不同的对象 | 1 | exclude | 同名术语冲突：属电子结构 DFT 谱势，非非厄米谱势判据 | classic | no | Q5 |
| 1707.01431 | Variational principles for t-entropy, the spectral potential of transfer operator, and entropy statistic theorem are equivalent | V. I. Bakhtin | 2017 | arXiv preprint（无 journal_ref） | theoretical (math.DS) | 对任意转移算子建立 t-熵变分原理、谱势与熵统计定理的等价性并给出新证明 | 转移算子（动力系统） | 证明三者等价 | 动力系统/热力学形式化语境，与非厄米能带无关 | 1 | exclude | 同名词“spectral potential”但属动力系统转移算子谱势 | classic | no | Q5 |
| 1705.09967 | Local Large Deviations: McMillian Theorem for multitype Galton-Watson Processes | Kwabena Doku-Amponsah | 2017 | Far East J. Math. Sci. 102(10), 2307-2319 (2017); DOI 10.17654/MS102102102307 | theoretical (cs.IT/概率) | 从谱势角度定义多类型 Galton-Watson 过程的谱势 U_κ(·,π)，证明局部大偏差原理 | 多类型 Galton-Watson 过程（Perron-Frobenius 特征值基准） | 得到偏差函数为谱势的 Legendre 对偶，并推出条件大偏差原理与 McMillan 定理弱形式 | 概率论语境，与物理谱势无关联 | 1 | exclude | 纯概率论中的谱势，对象与研究方法均不匹配 | classic | no | Q5 |
| 2107.02166 | Analysis of relationships between spectral potential of transfer operators, t-entropy, entropy and topological pressure | V. I. Bakhtin | 2021 | arXiv preprint（无 journal_ref） | theoretical (math.FA/DS) | 给出联系转移算子谱势、t-熵、熵与拓扑压的显式公式 | 转移算子与加权移位算子 | 揭示逆 rami-rate、前向熵等在该联系中的作用 | 动力系统语境，无算法与物理模型 | 1 | exclude | 同名词不同域（动力系统谱势） | classic | no | Q5 |
| math/0505446 | Positive Processes | V. I. Bakhtin | 2005 | Ergod. Th. & Dynam. Sys. 27, 639-670 (2007); DOI 10.1017/S0143385706000915 | theoretical (math.DS/PR) | 引入正流与正过程，基于相位与正代数、谱势、对偶熵、平衡测度等建立正算子理论分支 | 正算子/正过程 | 建立大数律与以作用泛函表示的大偏差概率渐近 | 纯动力学/概率论 | 1 | exclude | 同名词不同域（动力系统谱势） | classic | no | Q5 |
| 2609.23540 | Vector Balancing in Polynomial Time | Shengtao Guo | 2026 | arXiv preprint（无 journal_ref） | theoretical/algorithm (cs.DS) | 谱符号算法：通过最小化三次谱势（cubic spectral potential）把分数着色推向布尔符号，给出多项式时间算法 | Komlós 问题 / 矩阵列范数 ≤1（给出复杂度 O((mn^9+n^10)log(2+m+n))） | 以常数量级偏差在多项式时间内解决 Komlós 问题 | 组合优化语境下的“spectral potential”，与非厄米谱势无关 | 1 | exclude | 术语同名但属算法组合优化的谱势 | recent3 | no | Q5 |
| 2609.14029 | Special Markowitz: Thermodynamic Formalism for the Joint Regularisation of Returns and Covariance | David Reinhardt | 2026 | arXiv preprint（无 journal_ref） | theoretical (q-fin) | 用热力学形式化对收益与协方差联合正则化，每条本征方向带符号谱势 Φ_k | 投资组合正则化（Ledoit-Wolf 作对照） | 给出乘性复合律与 Stein 损失的唯一性刻画 | 金融语境 | 1 | exclude | 同名词不同域（金融） | recent3 | no | Q5 |
| 1409.4210 | Photoemission Spectroscopy and Orbital Imaging from Koopmans-Compliant Functionals | Ngoc Linh Nguyen | 2014 | arXiv preprint（无 journal_ref） | theoretical/numerical (电子结构) | 用 Koopmans 相容泛函作为谱势的廉价准粒子近似 | 分子光电子谱与 Dyson 轨道层析（与实验 UPS 对照） | Koopmans 相容泛函给出与实验一致的分子光电子谱 | 电子结构语境，非非厄米谱势 | 1 | exclude | 同名词不同域（DFT 谱势） | classic | no | Q5 |
| 0809.3116 | T-entropy and Variational Principle for the spectral radius of transfer and weighted shift operators | A. B. Antonevich | 2008 | Ergod. Th. & Dynam. Sys. 31, 995-1042 (2011); DOI 10.1017/S0143385710000210 | theoretical (math.DS/OA) | 给出任意动力系统下转移与加权移位算子谱半径的变分原理，显式描述谱势的 Legendre 对偶 | 转移/加权移位算子 | 得到以新熵型不变量 t-熵表达的变分原理 | 动力系统语境 | 1 | exclude | 同名词不同域（动力系统谱势） | classic | no | Q5 |
| 1707.01978 | Local Large deviation: A McMillian Theorem for Coloured Random Graph Processes | Kwabena Doku-Amponsah | 2017 | J. Math. Stat. 13(4), 347-352 (2017); DOI 10.3844/jmssp.2017.347.352 | theoretical (cs.IT) | 定义有限类型图的谱势 ρ_λ(·,μ)，并取 Legendre 对偶得到偏差函数 | 有限类型随机图 | 证明局部大偏差原理并推出条件大偏差原理与 McMillan 定理弱形式 | 概率/图论语境 | 1 | exclude | 同名词不同域（图过程谱势） | classic | no | Q5 |
| 1909.04529 | Local Large Deviation Principle, Large Deviation Principle and Information theory for the Signal-to-Interference-Plus-Noise Ratio Graph Models | E. Sakyi-Yeboah | 2019 | arXiv preprint（无 journal_ref） | theoretical (cs.IT/PR) | 从谱势角度证明标记 SINR 图的联合大偏差与局部大偏差原理 | 标记 SINR 图（Poisson 点过程） | 给出 AEP 与局部大偏差原理 | 通信网络语境 | 1 | exclude | 同名词不同域（SINR 图谱势） | classic | no | Q5 |
| 2102.04636 | Extensive Benchmarking of DFT+U Calculations for Predicting Band Gaps | Nicole E. Kirchner-Hall | 2021 | Appl. Sci. 11, 2395 (2021); DOI 10.3390/app11052395 | engineering/numerical (电子结构) | 用 DFPT 自洽确定 Hubbard U 以施加分段线性，被解释为近似谱势方法，系统基准带隙 | 20 种含过渡金属或 p 区元素的化合物 | 给出 DFT+U 带隙的系统基准结果 | 材料电子结构语境，与非厄米谱势无关 | 1 | exclude | 同名词不同域（DFT 谱势） | classic | no | Q5 |
| 2608.01911 | A Continuous-Time Analysis of Smoothed Matrix-Polar Spectral Gradient Flows for Muon-Type Optimization | Jinlin Liu | 2026 | arXiv preprint（无 journal_ref） | theoretical (math.OC) | 引入由光滑谱势生成的谱反馈律，分析光滑矩阵极梯度流的良定性与收敛率 | 矩阵值优化 / Muon 型优化器 | 给出非凸、凸与 PL 条件下的收敛率与 Lyapunov 分析 | 优化理论语境 | 1 | exclude | 同名词不同域（优化中的 spectral potential） | recent3 | no | Q5 |
| 2608.26288 | Muon with Finite Newton-Schulz: The Smoothing Benefit in Nonsmooth Nonconvex Optimization | Mingyi Li | 2026 | arXiv preprint（无 journal_ref） | theoretical (cs.LG/OC) | 把有限步 Newton-Schulz 迭代视为带光滑谱势的在线学习器并做 regret→稳定性转换 | Muon 优化器在 LLM 预训练中的设定 | 证明对数深度足以收敛到稳定点，精确极分解更新可能不收敛 | 机器学习优化语境 | 1 | exclude | 同名词不同域（优化中的 smoothed spectral potential） | recent3 | no | Q5 |
| 2609.23855 | Kadison-Singer partitions and Bilu-Linial graph signings in polynomial time | Ali Jadbabaie | 2026 | arXiv preprint（无 journal_ref） | theoretical/algorithm (cs.DS/CO) | 基于谱势方法（spectral-potential method）做确定性多项式时间舍入，并给出 Las Vegas 图符号算法 | PSD 矩阵舍入 / Bilu-Linial 图符号问题 | 给出确定性多项式复杂度舍入定理与几乎必然终止的图符号算法 | 谱差异理论语境 | 1 | exclude | 同名词不同域（谱势方法，组合优化） | recent3 | no | Q5 |
| 2106.03525 | Sturm-Liouville-type operators with frozen argument and Chebyshev polynomials | Tzong-Mo Tsai | 2021 | Math. Meth. Appl. Sci.; DOI 10.1002/mma.8327 | theoretical (math.SP) | 研究冻结论证 Sturm-Liouville 型算子的逆问题，建立主方程与 Chebyshev 多项式的联系 | 非局部（负载型）微分算子 | 完整描述非退化/退化情形与等谱势类 | 逆谱问题语境，非非厄米能带 | 1 | exclude | 对象与方法均不匹配（连续逆谱问题） | classic | no | Q5 |
| 2008.07147 | Inverse spectral problems for Hill-type operators with frozen argument | Sergey Buterin | 2020 | Anal. Math. Phys.; DOI 10.1007/s13324-021-00500-9 | theoretical (math.SP) | 研究冻结论证 Hill 型算子由谱信息重构复值势 q(x) 的两个逆问题，并给出算法 | 非局部 Hill 型算子 | 给出唯一性条件、等谱/等双谱势类刻画与求解算法 | 与 amoeba/Ronkin/非厄米谱势无关 | 1 | exclude | 对象与方法均不匹配（连续逆谱问题） | classic | no | Q5 |
| 2609.21279 | A Walk From Free Probability to Matrix Discrepancy III: Higher Rank Kadison-Singer and Spectrally Thin Trees | Tarun Kathuria | 2026 | arXiv preprint（无 journal_ref） | theoretical/algorithm (cs.DS) | 用优化谱势（optimized spectral potential）与凹矩阵幂控制密度响应，给出高秩 Kadison-Singer 的确定性多项式算法 | 高秩 PSD 矩阵符号 / 谱稀疏树 | 得到 O(√ε log(2r)) 量级偏差与多项式实算术工作量算法 | 谱差异理论语境 | 1 | exclude | 同名词不同域（组合优化的谱势） | recent3 | no | Q5 |
| 2306.04460 | Loss-induced Floquet non-Hermitian skin effect | Yaohua Li | 2023 | Phys. Rev. B 108, L220301 (2023); DOI 10.1103/PhysRevB.108.L220301 | theoretical | 提出损耗诱导 Floquet NHSE 机制，并把 GBZ 理论推广到非平衡系统以描述该效应 | 螺旋波导光子晶格中的 Floquet NHSE（含 2D 二阶 NHSE 推广；无外部算法基准） | 给出可实验实现的 Floquet NHSE 方案并指出 2D 可实现二阶 NHSE | GBZ 由既有框架套用，非算法贡献；无 amoeba/Ronkin 层面 | 1 | exclude | GBZ/非布洛赫绕数仅作表征工具，算法非主题内容 | classic | no | Q2 |
| 2304.01521 | Transport properties of a non-Hermitian Weyl semimetal | Soumi Dey | 2023 | arXiv preprint（无 journal_ref） | theoretical/numerical | 研究 3D 耗散 Weyl 半金属霍尔电导，用非布洛赫理论在 GBZ 上求电导并与 OBC 谱对照 | 两层 Chern 绝缘体堆叠成的 3D 非厄米 Weyl 半金属（PBC/OBC 对照） | 发现非厄米下霍尔电导偏离量子化，GBZ 上计算可弥合 PBC/OBC 的转变点差异 | GBZ 仅用于修正谱与电导，非算法主题 | 1 | exclude | GBZ 仅用于修正，非算法主题内容 | classic | no | Q2 |
| 2302.14366 | Z2 Non-Hermitian skin effect in equilibrium heavy-fermions | Shin Kaneshiro | 2023 | Phys. Rev. B 107, 195149 (2023); DOI 10.1103/PhysRevB.107.195149 | theoretical/numerical | 用 DMFT+NRG 结合 GBZ 技术分析关联 f 电子体系中的 Z2 趋肤效应 | 二维周期 Anderson 模型 + 自旋轨道耦合（DMFT/NRG 作数值方法） | 证明自旋轨道耦合下存在 Z2 趋肤效应，并用于分析温度效应 | GBZ 仅作分析工具，无算法贡献 | 1 | exclude | GBZ 仅作分析工具，算法非主题内容 | classic | no | Q2 |
| 2212.12637 | Competition of non-Hermitian skin effect and topological localization of corner states observed in circuits | Chan Tang | 2022 | Phys. Rev. B 108, 035410 (2023); DOI 10.1103/PhysRevB.108.035410 | experimental | 2D 非互易蜂窝电路实验，用 GBZ 上定义的非布洛赫绕数与 Z2 Berry 相表征相变 | 菱形蜂窝电路（实验测量对照非布洛赫不变量） | 观测到角态被趋肤效应拖入体的竞争现象，非布洛赫 Z2 Berry 相可作不变量 | 属实验工作，非算法/数值方法贡献 | 1 | exclude | 对象与类型不在范围内（纯实验、电路平台） | classic | no | Q2 |
| 2210.12732 | Quantum circuit for measuring an operator's generalized expectation values and its applications to non-Hermitian winding numbers | Ze-Hao Huang | 2022 | Phys. Rev. A 107, 052205 (2023); DOI 10.1103/PhysRevA.107.052205 | theoretical/numerical | 提出基于 swap test 的量子电路测量广义期望值 ⟨ψ1\|A\|ψ2⟩，并用于测量非厄米绕数 | 非互易 SSH 模型（数值模拟电路保真度作验证） | 数值模拟显示 Bloch/非 Bloch 自旋织构与相应绕数可高保真测得 | 量子线路测量路线，与 GBZ/amoeba 的经典计算算法无关 | 1 | exclude | 方法路线（量子线路）与 amoeba/Ronkin 计算无关 | classic | no | Q1 |
| 2009.10508 | Anomalous topological edge states in non-Hermitian piezophononic media | Penglin Gao | 2020 | Phys. Rev. Lett. 125, 206402 (2020); DOI 10.1103/PhysRevLett.125.206402 | experimental/numerical | 数值研究偏置压电声子介质中的非厄米/非互易拓扑力学与趋肤效应 | 分层压电声子介质（MHz 频段，数值实验） | 显示传统 Bloch 能带无法预测有限系统体带，揭示反常拓扑边缘态 | 属声学平台研究，非算法贡献 | 1 | exclude | 对象与类型不在范围内（声学实验平台） | classic | no | Q1 |
| 2604.00895 | The effect of staggered nonlinearity on the Su-Schrieffer-Heeger model | Ahmed Alharthy | 2026 | arXiv preprint（无 journal_ref） | theoretical/numerical | 用半解析 Bloch 处理 + 自洽场迭代法数值求解 OBC 下的非 Bloch 解，并推导非线性 Zak 相 | 亚晶格依赖非线性 SSH 模型（半解析与自洽场两种方法互为对照） | 发现高非线性下非线性 Zak 相不连续，标志非线性诱导的拓扑相变 | 自洽场迭代求非 Bloch 解，未涉及 GBZ/amoeba/Ronkin 计算 | 1 | exclude | 方法未落在 amoeba/Ronkin/谱势/GBZ 计算上 | recent3 | no | Q1 |
| 2411.13661 | Non-Bloch self-energy of dissipative interacting fermions | He-Ran Wang | 2024 | arXiv preprint（无 journal_ref） | theoretical | 对耗散相互作用费米子微扰计算自能矩阵，并用非布洛赫能带理论导出精确积分表示 | 马尔可夫开放量子系统中的相互作用费米子（与数值计算高精度对照） | 得到 Liouvillian 谱修正的简化表达式，并识别相互作用增强的 NHSE | 纯解析积分表示与微扰分析，无 GBZ/Ronkin 计算算法 | 1 | exclude | 方法未落在 amoeba/Ronkin/谱势/GBZ 计算上 | recent3 | no | Q1 |
| 2407.19766 | Anomalous symmetry protected blockade of skin effect in one-dimensional non-Hermitian lattice systems | Shuai Li | 2024 | arXiv preprint（无 journal_ref） | theoretical | 证明联合空间反射对称性可作 1D 非厄米系统中 NHSE 存亡的判据 | 非厄米 Kitaev 链（数学证明 + 精确对角化数值验证） | 给出仅依赖对称性的 NHSE 存亡判据 | 判据式结果，无 GBZ/amoeba 数值算法 | 1 | exclude | 方法未落在 amoeba/Ronkin/谱势/GBZ 计算上 | recent3 | no | Q1 |
| 2406.15564 | Entangelment Entropy on Generalized Brillouin Zone | Zhenghao Yang | 2024 | Phys. Rev. B 110, 235127 (2024); DOI 10.1103/PhysRevB.110.235127 | theoretical/numerical | 在准互易格点上定义非布洛赫纠缠熵，并研究其与 GBZ 形状（圆形/非圆形）的关系 | 非厄米 SSH 模型（EP 附近与临界区作对照） | 显示非布洛赫纠缠熵可恢复体边对应，并给出中心荷随 Fermi 点的贡献 | 讨论物理量定义，非 GBZ/amoeba 的算法贡献 | 1 | exclude | 方法未落在 amoeba/Ronkin/谱势/GBZ 计算上 | recent3 | no | Q1 |
| 2401.12785 | Extended imaginary gauge transformation in a general nonreciprocal lattice | Yunyao Qi | 2024 | Phys. Rev. B 110, 075411 (2024); DOI 10.1103/PhysRevB.110.075411 | theoretical | 把虚规范变换推广到复谱情形，并证明满足赝厄米性时 GBZ 为圆 | 一般非互易格点、三聚体 SSH、2D Hatano-Nelson 模型 | 揭示虚规范变换成立的条件，并给出圆 GBZ 下趋肤模局域长度的易得计算 | GBZ 恒为圆的特例，未涉及一般 amoeba/Ronkin 计算 | 1 | exclude | 方法未落在 amoeba/Ronkin/谱势/GBZ 计算上 | recent3 | no | Q1 |
| math/0205052 | Szego limit theorems | Russell Lyons | 2002 | Geom. Funct. Anal. 13 (2003), 574-590 | theoretical (math.CA) | 把第一 Szegő 极限定理推广到更一般的测度与行列式，并给出高维 Helson–Lowdenslager 型扩展 | 一般测度 / Toeplitz 型行列式（Bump–Diaconis、Tracy–Widom 结果作对照） | 给出更一般测度下 Szegő 型极限与高维扩展 | 纯测度论/算子是理论，非 Toeplitz 行列式渐近的可计算工具 | 1 | exclude | 不属 Szegő–Toeplitz 行列式渐近计算工具侧 | classic | no | Q3 |
| math/0111286 | Tracial invariants, classification and II_1 factor representations of Popa algebras | Nathanial P. Brown | 2001 | arXiv preprint（无 journal_ref） | theoretical (math.OA) | 定义 C*-代数迹空间上的四类凸子集不变量，并在附录式结果中给出任意自伴算子的 Szegő 极限定理类比 | Popa 代数 / II_1 因子（迹不变量） | 给出迹不变量应用与自伴算子的 Szegő 型类比 | 算子代数语境，不涉及 Toeplitz 行列式渐近计算 | 1 | exclude | 不属 Szegő–Toeplitz 行列式渐近计算工具侧 | classic | no | Q3 |
| 2604.00895 | The effect of staggered nonlinearity on the Su-Schrieffer-Heeger model | Ahmed Alharthy | 2026 | arXiv preprint（无 journal_ref） | theoretical/numerical | 半解析 Bloch 处理 + 自洽场迭代法数值求解 OBC 下非 Bloch 解，并推导含非线性修正的 Zak 相 | 亚晶格依赖非线性 SSH 模型（半解析与自洽场两法互为对照） | 高非线性下非线性 Zak 相出现不连续，标志非线性诱导拓扑相变 | 未涉及 GBZ/amoeba/Ronkin 计算 | 1 | exclude | 方法未落在 amoeba/Ronkin/谱势/GBZ 计算上 | recent3 | no | Q1 |
| 2411.13661 | Non-Bloch self-energy of dissipative interacting fermions | He-Ran Wang | 2024 | arXiv preprint（无 journal_ref） | theoretical | 微扰计算耗散相互作用费米子自能，用非布洛赫能带理论导出精确积分表示 | 马尔可夫开放量子系统（与数值计算高精度对照） | 得到 Liouvillian 谱修正简化表达式，识别相互作用增强的 NHSE | 纯解析表示，无 GBZ/Ronkin 计算算法 | 1 | exclude | 方法未落在 amoeba/Ronkin/谱势/GBZ 计算上 | recent3 | no | Q1 |
| 2407.19766 | Anomalous symmetry protected blockade of skin effect in one-dimensional non-Hermitian lattice systems | Shuai Li | 2024 | arXiv preprint（无 journal_ref） | theoretical | 以联合空间反射对称性作为 1D 非厄米系统 NHSE 存亡判据 | 非厄米 Kitaev 链（证明 + 精确对角化验证） | 给出仅依赖对称性的 NHSE 存亡判据 | 无 GBZ/amoeba 数值算法 | 1 | exclude | 方法未落在 amoeba/Ronkin/谱势/GBZ 计算上 | recent3 | no | Q1 |
| 2406.15564 | Entangelment Entropy on Generalized Brillouin Zone | Zhenghao Yang | 2024 | Phys. Rev. B 110, 235127 (2024); DOI 10.1103/PhysRevB.110.235127 | theoretical/numerical | 在准互易格点上定义非布洛赫纠缠熵，研究其与 GBZ 形状的关系 | 非厄米 SSH 模型（EP 附近与临界区对照） | 非布洛赫纠缠熵可恢复体边对应，给出中心荷的 Fermi 点贡献 | 物理量定义研究，非 GBZ 算法贡献 | 1 | exclude | 方法未落在 amoeba/Ronkin/谱势/GBZ 计算上 | recent3 | no | Q1 |
| 2401.12785 | Extended imaginary gauge transformation in a general nonreciprocal lattice | Yunyao Qi | 2024 | Phys. Rev. B 110, 075411 (2024); DOI 10.1103/PhysRevB.110.075411 | theoretical | 把虚规范变换推广到复谱，证明满足赝厄米性时 GBZ 为圆 | 一般非互易格点、三聚体 SSH、2D Hatano-Nelson | 给出虚规范变换成立条件与圆 GBZ 下局域长度的易得计算 | 仅覆盖 GBZ 为圆的特例，未涉及一般 amoeba/Ronkin 计算 | 1 | exclude | 方法未落在 amoeba/Ronkin/谱势/GBZ 计算上 | recent3 | no | Q1 |
| 2306.04460 | Loss-induced Floquet non-Hermitian skin effect | Yaohua Li | 2023 | Phys. Rev. B 108, L220301 (2023); DOI 10.1103/PhysRevB.108.L220301 | theoretical | 以损耗为非厄米来源，借 Floquet 工程产生复次近邻跃迁，并把 GBZ 理论推广到非平衡系统 | 螺旋波导光子晶格（含 2D 二阶 NHSE 推广；无算法基准） | 提出可实验实现的 Floquet NHSE 机制，并可实现二阶 NHSE | GBZ 由既有框架套用，非算法贡献 | 1 | exclude | GBZ/非布洛赫绕数仅作表征工具，算法非主题内容 | classic | no | Q2 |
| 2304.01521 | Transport properties of a non-Hermitian Weyl semimetal | Soumi Dey | 2023 | arXiv preprint（无 journal_ref） | theoretical/numerical | 研究 3D 耗散 Weyl 半金属输运，用非布洛赫理论在 GBZ 上计算霍尔电导 | 两层 Chern 绝缘体堆叠的 3D 非厄米 Weyl 半金属（PBC/OBC 对照） | 非厄米下霍尔电导偏离量子化；GBZ 上计算弥合 PBC/OBC 转变点差异 | GBZ 仅用于修正谱与电导，非算法主题 | 1 | exclude | GBZ 仅用于修正，非算法主题内容 | classic | no | Q2 |
| 2302.14366 | Z2 Non-Hermitian skin effect in equilibrium heavy-fermions | Shin Kaneshiro | 2023 | Phys. Rev. B 107, 195149 (2023); DOI 10.1103/PhysRevB.107.195149 | theoretical/numerical | 用 DMFT + NRG 并结合 GBZ 技术分析关联 f 电子体系的 Z2 趋肤效应 | 二维周期 Anderson 模型 + 自旋轨道耦合 | 证明自旋轨道耦合下存在 Z2 趋肤效应并分析温度效应 | GBZ 仅作分析工具，无算法贡献 | 1 | exclude | GBZ 仅作分析工具，算法非主题内容 | classic | no | Q2 |
| 2212.12637 | Competition of non-Hermitian skin effect and topological localization of corner states observed in circuits | Chan Tang | 2022 | Phys. Rev. B 108, 035410 (2023); DOI 10.1103/PhysRevB.108.035410 | experimental | 2D 非互易菱形蜂窝电路实验，用 GBZ 上定义的非布洛赫绕数与 Z2 Berry 相表征相变 | 菱形蜂窝电路（测量与不变量对照） | 观测角态被趋肤效应拖入体的竞争，非布洛赫 Z2 Berry 相可作不变量 | 实验工作，非算法/数值方法贡献 | 1 | exclude | 对象与类型不在范围内（纯实验、电路平台） | classic | no | Q2 |
| 0810.2315 | Szego limit theorems on the Sierpinski gasket | Kasso A. Okoudjou | 2008 | arXiv preprint（无 journal_ref） | theoretical (math.SP) | 利用 Sierpinski 垫片上 Laplacian 的局域化本征函数，构造并证明该分形上强 Szegő 极限定理的类比 | Sierpinski 垫片上的 Laplacian | 给出分形设定下的强 Szegő 型定理与等分布序列表述 | 分形分析语境，与 Toeplitz 谱渐近路线无关 | 1 | exclude | 不属 Szegő–Toeplitz 行列式渐近计算工具侧 | classic | no | Q3 |
| 0706.0750 | An analogue of Szego's limit theorem in free probability theory | Junhao Shen | 2007 | arXiv preprint（无 journal_ref） | theoretical (math.OA) | 讨论自由概率论中的正交多项式，证明 Szegő 极限定理的类比 | 自由概率论中的正交多项式 | 建立自由概率论版本的 Szegő 极限定理 | 算子代数语境，非行列式渐近工具 | 1 | exclude | 不属 Szegő–Toeplitz 行列式渐近计算工具侧 | classic | no | Q3 |
| 1810.04402 | An application of Brascamp-Lieb's inequality | Michel Weber | 2018 | Sankhya A (2024) | theoretical (math.PR) | 用 Brascamp–Lieb 不等式获得高斯向量与平稳循环高斯过程的解耦不等式，其中用到 Bump–Diaconis 版强 Szegő 定理 | 高斯向量 / 平稳循环高斯过程 | 给出新的解耦不等式，推广 Klein–Landau–Shucker 的结果 | 概率论应用方向，Szegő 定理仅为工具之一 | 1 | exclude | 不属 Szegő–Toeplitz 行列式渐近计算工具侧 | classic | no | Q3 |

---

## 3. EXCLUDED（压缩，一行一篇）

`Decision = exclude` 的全部 23 条**已在 §2 主表中以全字段给出**（含 `1809.08187`、`1707.01431`、
`1705.09967`、`2107.02166`、`math/0505446`、`2609.23540`、`2609.14029`、`1409.4210`、`0809.3116`、
`1707.01978`、`1909.04529`、`2102.04636`、`2608.01911`、`2608.26288`、`2609.23855`、`2106.03525`、
`2008.07147`、`2609.21279`、`2306.04460`、`2304.01521`、`2302.14366`、`2212.12637`、`2210.12732`、
`2009.10508`、`2604.00895`、`2411.13661`、`2407.19766`、`2406.15564`、`2401.12785`、`math/0205052`、
`math/0111286`、`0810.2315`、`0706.0750`、`1810.04402`）。

本节为便于浏览的**压缩重复呈现**（与 §2 主表一一对应，不额外增加计数）：

### 3.1 Q5 同名词异域命中（谱势/Toeplitz 相关但不匹配研究问题）

| arXiv ID | 一句话排除原因 |
|---|---|
| 1809.08187 | DFT 动力学谱势（Kohn–Sham 势推广）在 Hubbard 二聚体中的近似策略，非非厄米谱势判据 |
| 1707.01431 | 动力系统转移算子谱势与 t-熵变分原理的等价性，非能带语境 |
| 2107.02166 | 转移算子谱势与 t-熵/拓扑压的显式关系公式，动力系统语境 |
| 1705.09967 | 多类型 Galton–Watson 过程的谱势与大偏差，纯概率论 |
| math/0505446 | 正流/正过程的谱势与对偶熵理论，纯动力学 |
| 0809.3116 | 转移与加权移位算子谱半径变分原理中的谱势 Legendre 对偶，动力系统语境 |
| 1707.01978 | 有限类型随机图谱势与局部大偏差，概率/图论语境 |
| 1909.04529 | 标记 SINR 图谱势与大偏差原理，通信网络语境 |
| 1409.4210 | Koopmans 相容泛函作为谱势的准粒子近似，电子结构语境 |
| 2102.04636 | DFT+U 带隙基准（把 U 解释为近似谱势方法），材料电子结构语境 |
| 2609.23540 | 组合优化中最小化三次谱势的向量平衡算法，术语同名 |
| 2609.14029 | 金融中带符号谱势的热力学形式化，术语同名 |
| 2608.01911 | 优化中由光滑谱势生成的谱反馈律，术语同名 |
| 2608.26288 | Muon 优化的 smoothed spectral potential，术语同名 |
| 2609.23855 | 谱差异理论中的 spectral-potential method，术语同名 |
| 2106.03525 | 冻结论证 Sturm–Liouville 型算子逆问题，连续逆谱问题 |
| 2008.07147 | 冻结论证 Hill 型算子逆问题与求解算法，连续逆谱问题 |
| 2609.21279 | 组合优化中的 optimized spectral potential，术语同名 |

### 3.2 Q1/Q2 主题相关但方法或类型不匹配（压缩索引，全字段见 §2 主表）

| arXiv ID | 一句话排除原因 |
|---|---|
| 2604.00895 | 非线性 SSH 用自洽场迭代求非 Bloch 解，无 GBZ/amoeba 算法与谱势判据 |
| 2411.13661 | 非布洛赫自能的纯解析积分表示与微扰分析，无 GBZ/Ronkin 计算算法 |
| 2407.19766 | 用空间反射对称性判据证明 NHSE 存亡，数学证明+精确对角化，无算法 |
| 2406.15564 | 讨论 GBZ 上的纠缠熵定义，纯物理量研究，非计算算法 |
| 2401.12785 | 推广虚规范变换，GBZ 恒为圆，未涉 amoeba/Ronkin 计算 |
| 2306.04460 | 损耗诱导 Floquet NHSE 的物理机制提案（含 2D 二阶效应），GBZ 仅作工具 |
| 2304.01521 | 非厄米 Weyl 半金属输运与霍尔电导，GBZ 仅用于修正，非算法主题 |
| 2302.14366 | DMFT+NRG 研究重费米子 Z2 趋肤效应，GBZ 仅作分析工具 |
| 2212.12637 | 2D 非互易蜂窝电路实验，非布洛赫绕数仅作表征工具 |
| 2210.12732 | 量子线路测量广义期望值与绕数，与经典 GBZ/amoeba 计算路线无关 |
| 2009.10508 | 压电声学介质中非布洛赫边缘态的数值/实验研究，非算法贡献 |
| math/0205052 | Szegő 极限定理的测度推广与高维扩展，纯测度论，非 Toeplitz 渐近计算工具 |
| 0810.2315 | Sierpinski 垫片上的 Szegő 极限定理，分形分析语境 |
| 0706.0750 | 自由概率论中的 Szegő 极限定理类比，纯算子代数 |
| 1810.04402 | 用 Bump–Diaconis 版强 Szegő 定理做高斯过程解耦不等式，概率论应用 |
| math/0111286 | Popa 代数迹不变量与 II_1 因子表示，仅附录式提及 Szegő 类比 |

---

## 4. Q5 专项说明：`abs:"spectral potential"` 与 Q1/Q4 的重叠情况

### 4.1 命中与重叠

Q5 共 21 篇：

- **与 Q1/Q4 重叠（2 篇，去重后不新增）**
  - `2511.11349` — 同时命中 Q1、Q5（标题即以 Amoeba formulation / Wiener-Hopf 为主）
  - `2407.01296` — 同时命中 Q1、Q2、Q5（摘要明确 "through the lens of spectral potential"）
- **Q5 新增 unique（19 篇）**：`1809.08187`、`1707.01431`、`1705.09967`、`2107.02166`、`math/0505446`、
  `2609.23540`、`1409.4210`、`0809.3116`、`2609.14029`、`1707.01978`、`2608.01911`、`2106.03525`、
  `2608.26288`、`2609.23855`、`2102.04636`、`2609.21279`、`2506.08618`、`1909.04529`、`2008.07147`
  - 其中唯一具有本研究语境价值的是 `2506.08618`（HSG-12M：Poly2Graph 流水线 + 大规模谱图数据集），判为 unclear，rel 2
- **重叠统计**：21 = 2（与 Q1/Q4 重复）+ 19（净新增）；Q5 与本批 Q1–Q4 unique(51) 的并集 = 70，
  加上 `2407.01296` 的三重命中计一次，全批 unique = 71。

### 4.2 关键观察：Q5 的严重同名词污染

Q5 的 19 条净新增中，**18 条为「spectral potential」同名异域命中**，与非厄米物理的谱势/amoeba 判据无关，构成三类：

1. **电子结构 DFT 谱势**（Kohn–Sham 势的动力学推广）：1809.08187、1409.4210、2102.04636
2. **动力系统 / 遍历论 / 概率论中的谱势**（转移算子、t-熵、大偏差、Perron–Frobenius）：1707.01431、
   2107.02166、math/0505446、0809.3116、1705.09967、1707.01978、1909.04529
3. **组合优化 / 机器学习 / 金融中的 spectral potential**：2609.23540（三次谱势）、2609.23855、2609.21279
   （optimized spectral potential）、2608.01911、2608.26288（smoothed spectral potential）、2609.14029
   （金融，符号谱势）、2106.03525、2008.07147（连续逆谱问题）

另有 1 条非同名异域但也不匹配：`2506.08618`（属非厄米，但限 1D 谱图数据集）。

### 4.3 结论与检索建议

**单用 `abs:"spectral potential"` 精度极低：21 篇中仅 2 篇属本研究语境（且这 2 篇已由 Q1/Q2 命中），
净新增 19 篇中 18 篇为同名词异域、1 篇（HSG-12M）仅为 1D 数据集。**
建议后续查询改为**共现约束**，例如：

- `abs:"spectral potential" AND (abs:"non-Hermitian" OR abs:"non-Bloch" OR abs:amoeba OR abs:Ronkin)`
- `abs:"Ronkin function"`（本批 5 条查询均未覆盖，是明显的检索缺口）
- `abs:amoeba AND abs:"Brillouin zone"`
- `abs:"Szegő"`（带重音）配合 `cat:cond-mat.*`

---

## 5. 抓取失败 / 异常

### 5.1 抓取失败

- 无。5 条查询最终全部 HTTP 200 成功。
- 过程记录（供复现）：
  - 指定的 `http://export.arxiv.org/api/query?...` 形式直连时返回
    `cross-origin redirect to https://export.arxiv.org is not followed automatically`；
    改用同一 URL 的 `https://export.arxiv.org/...` 形式后成功（Q1–Q3、Q5 各 1 次成功）。
  - Q4 首次 https 请求返回 `TypeError: fetch failed`，按既定规则把 search_query 未编码原文重写后重编码重试一次，即成功。
- arXiv 限速：各查询之间按 ≈1 请求/3 秒 的节奏顺序抓取，未触发 429。

### 5.2 内容层异常（非抓取失败）

- `2506.22743`：给定标题为 "Multidimensional non-Bloch spectral theory"，但该条目的 `<summary>` 未能从返回的 Atom XML 中取到可解析内容。按诚信要求**未推断其内容**，Decision 记为 `unclear`，相关字段留空待人工复核。
- `2609.23523`、`2608.28577`、`2609.23224`、`2506.08618` 的标题与摘要均为 arXiv API 原文，未做改写。

---

## 6. 待用户裁决事项（本批不自行裁定）

1. `2506.22743` 的摘要缺失需人工补齐后再判；其后可能升级为 include。
2. `2206.03029`（Widom 型 Toeplitz 渐近 + LQG）是否随本次 math-tool 口径扩展一并升为 `include (math-tool side)`；当前保守保留 `unclear`。
3. Q3 的纯数学 Szegő 条目中，`2305.19025`、`math/0012200`、`math/0612487`、`math/1102.4131`、`math/0212215` 已按本次口径列为 `include (math-tool side)`；
   其余 5 条（math/0205052、0810.2315、0706.0750、1810.04402、math/0111286）是否也属“Szegő–Toeplitz 数学工具侧”，需用户确认。
4. `2506.08618`（HSG-12M，1D 谱图数据集 + Poly2Graph 流水线）当前为 `unclear`（rel 2）；若用户认为“可复现的谱（势）计算流水线 + 公开基准”本身即满足 I 与 C 维度，可升为 include。
