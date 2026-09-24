# Round-2 arXiv screening log — algorithmic computation of amoeba GBZ / Ronkin function

## Queries (verbatim, unencoded)
- Q1: `abs:amoeba AND (abs:"non-Hermitian" OR abs:"skin effect" OR abs:"Brillouin")` — totalResults=8, retrieved=8
- Q2: `abs:"Ronkin function"` — totalResults=21, retrieved=21
- Q3: `abs:"spectral potential" AND (abs:"non-Hermitian" OR abs:"skin" OR abs:"Toeplitz" OR abs:"Brillouin" OR abs:"lattice")` — totalResults=3, retrieved=3

Raw retrieved = 32; duplicates merged = 4 (2511.11349 x3, 2502.17931 x2, 2407.01296 x2);
unique_after_dedup = 28; include = 14; unclear = 8; exclude = 6.

Retrieval note: Q1 via `http://export.arxiv.org/...` failed with "cross-origin redirect to https://export.arxiv.org
is not followed automatically"; the same encoded query re-issued against `https://export.arxiv.org/...` returned HTTP 200.

## INCLUDED / UNCLEAR (sorted by relevance)

| arXiv ID | Title | First author | Year | Venue | Type | Method core | Model/benchmark | Main result | Limitation | Rel | Decision | Rationale | tier | round1 | hit_by |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2212.11743 | Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions | Hong-Yi Wang | 2022 | PRX 14, 021011 (2024), DOI 10.1103/PhysRevX.14.021011 | theoretical | 用 amoeba（代数几何对象）重构非布洛赫能带论：由 amoeba/Ronkin 函数确定 GBZ 与能谱 | 任意维度非厄米紧束缚模型 | 给出任意维度 amoeba 公式，可高效获得能谱、本征态轮廓与 GBZ | 摘要未给出数值算法细节、精度与复杂度分析 | 5 | include | 直接给出 amoeba 定义 GBZ 的构造，是本问题的核心文献 | classic | yes | Q1 |
| 2502.17931 | Symplectic-Amoeba formulation of the non-Bloch band theory for one-dimensional two-band systems | Shin Kaneshiro | 2025 | 无期刊信息 (v2 2025-06-27) | theoretical+numerical | 谱势计算归结为 Ronkin 函数优化；对多带/对称退化用"外推总 Ronkin 函数"得 band-resolved Ronkin 函数再逐一优化 | 1D 双带 AII† 类非厄米模型 | 提出 AII† 类广义 Szegő 极限定理，并数值验证谱势与局域长度计算正确 | 仅 1D 双带且依赖对称类；未报告复杂度/精度控制 | 5 | include | 明确讨论 Ronkin 函数的实际优化/外推实现并附数值验证 | recent3 | yes | Q1,Q2 |
| 2407.01296 | Non-Hermitian skin effect in arbitrary dimensions: non-Bloch band theory and classification | Yuncheng Xiong | 2024 | 无期刊信息 | theoretical | 以 spectral potential 构造几何自适应非布洛赫能带论，用净绕数分类 NHSE | 任意维度非厄米紧束缚（热力学极限） | 精确给出给定几何下的能谱、态密度与 GBZ，并揭示标度无关模导致的谱非收敛/不稳定（谱趋向 Amoeba 谱） | 摘要未给出数值实现细节与复杂度 | 5 | include | 直接针对谱势/GBZ 的数值收敛与稳定性，属算法层面核心 | recent3 | yes | Q1,Q3 |
| 2511.11349 | Wiener-Hopf factorization and non-Hermitian topology for Amoeba formulation in one-dimensional multiband systems | Shin Kaneshiro | 2025 | Phys. Rev. Research 8, 013292 (2026), DOI 10.1103/s43l-h6z6 | theoretical | 用非布洛赫哈密顿量的 Wiener-Hopf 分解 + Hermitian doubling 给出多带体系广义 Szegő 定理适用性判据，并解释 AII† 类对称分解 Ronkin 函数的来源 | 1D 多带非厄米模型 | 严格证明 AII† 类广义 Szegő 极限定理，为多带 amoeba 分析提供统一框架 | 摘要未给出具体数值算法/基准；推广到其他对称类仍开放 | 5 | include | 明确指出 amoeba 实现以往只限单带并给出多带判据（适用边界/失效模式） | recent3 | yes | Q1,Q2,Q3 |
| 1306.6249 | An explicit calculation of the Ronkin function | Johannes Lundqvist | 2013 | math.CV 预印本 | theoretical | 对三变量仿射线性多项式显式计算 Ronkin 函数二阶导数，用完全椭圆积分与超几何函数表示，半显式给出 Ronkin 测度 | 三变量仿射线性多项式 | 给出 Ronkin 函数二阶导数的椭圆积分/超几何函数闭式 | 仅仿射线性多项式情形，未涉及非厄米谱 | 5 | include | 标题即"显式计算 Ronkin 函数"，直接对应 Ronkin 函数求值问题 | classic | no | Q2 |
| 2608.28577 | How Long-Range Tails Reshape Non-Hermitian Spectra | Ding Gu | 2026 | 无期刊信息 | theoretical | 提出 squeezed GBZ / squeezed amoeba 表述，把无穷长程跳跃截断效应编码进被压缩的 amoeba 以重构开边界谱与谱密度 | 1D 与 ≥2D 含指数衰减长程跳跃的非厄米紧束缚模型 | 长程跳跃使 GBZ 被压缩，高维谱密度由 squeezed amoeba 给出并可改变格林函数 | 摘要未给出数值实现与精度控制 | 4 | include | 给出 amoeba 表述的可执行变体，是算法层面的扩展 | recent3 | no | Q1 |
| 2607.22976 | Spectral Topology and Non-Bloch Band Theory for Domain-Wall Systems | Mingtao Xu | 2026 | 无期刊信息 | theoretical | 把 Ronkin 函数形式推广到域壁环形构型，得到各类模的 GBZ 条件与谱绕数 | 1D 非厄米域壁/环几何模型 | 得到基于 Ronkin 函数形式化的 GBZ 条件，区分驻波型与行波型皮肤模 | 摘要未给出数值算法细节/基准 | 4 | include | 直接扩展 Ronkin 函数形式以解 GBZ，属算法层面推广 | recent3 | no | Q2 |
| 2402.08798 | Dimers and M-Curves | Alexander I. Bobenko | 2024 | math-ph 预印本 | theoretical+numerical | 用 M-曲线上的亚纯微分积分给出 Ronkin 函数与表面张力显式公式，并基于 Schottky 均匀化给出具体计算流程（算权重、采样 dimer 构型） | 双周期二分图 dimer 模型 / M-曲线族 | 给出 Ronkin 函数显式积分公式，计算结果与理论预言完全一致 | 对象为 dimer 表面张力而非厄米 GBZ；摘要未报告复杂度 | 4 | include | 少数给出 Ronkin 函数显式公式并附数值验证的文献，可作算法与基准参考 | recent3 | no | Q2 |
| 2212.06553 | Mahler Measuring the Genetic Code of Amoebae | Siqi Chen | 2022 | hep-th 预印本 | engineering/numerical | 用符号回归/机器学习数值提取 amoeba 分量（有界补集体积等）与 Mahler measure 的数值关系，并用 ML 直接学习 3d amoeba 拓扑；另给部分 amoeba 边界解析式 | 2d/3d amoeba（Newton 多项式，quiver/dimer 背景） | 有界补集体积与 Mahler measure 的 gas 相贡献满足 d 次多项式关系（d=2,3），ML 学习 3d 拓扑表现良好 | 数据驱动拟合而非严格算法；未报告复杂度与精度界 | 4 | include | 直接给出 amoeba 数值计算/学习的实现与对照，可作算法与基准参考 | classic | no | Q2 |
| 2506.08618 | HSG-12M: A Large-Scale Benchmark of Spatial Multigraphs from the Energy Spectra of Non-Hermitian Crystals | Xianquan Yan | 2025 | ICLR 2026 | engineering | 开源高性能流水线 Poly2Graph，自动把 1D 晶体哈密顿量映射为复平面谱图并构建数据集与 GNN 基准 | HSG-12M（11.6M 静态 + 5.1M 动态谱图，1401 个特征多项式类，源自 177 TB 谱势数据） | 首个大规模空间多重图数据集，暴露 GNN 学习空间多重边的新挑战 | 仅 1D 哈密顿量；关注图学习而非 GBZ/Ronkin 数值精度 | 4 | include | 提供谱势/谱图的规模化数值流水线与基准，是"数值实现与验证基准"的直接素材 | recent3 | no | Q3 |
| 2405.03750 | Anatomy of Higher-Order Non-Hermitian Skin and Boundary Modes | Fan Yang | 2024 | Phys. Rev. Research 7, 023233 (2025), DOI 10.1103/PhysRevResearch.7.023233 | theoretical | 对任意维度一大类模型的精确解析解，结合考虑 separation gap 的 amoeba 理论及双正交极化/GBZ 的高维推广，追踪体模与边界模的拓扑起源 | d 维非厄米模型（d_c ≤ d 个方向开边界） | 给出 amoeba 理论 + 高维 GBZ 推广的统一图景，可分辨体模与边界模 | 以解析解为主，摘要未给出通用数值算法与复杂度 | 4 | include | 在高维实际使用 amoeba 理论区分体/边界谱，是算法适用范围的参考 | recent3 | yes | Q1 |
| 2212.13704 | Ronkin/Zeta Correspondence | Takashi Komatsu | 2022 | math-ph 预印本 | theoretical | 建立 Ronkin 函数与随机游走/量子游走 zeta 函数之间的新关系（1D、高维随机游走；1D 量子游走） | 1D/高维随机游走与 1D 量子游走 | 首次通过量子游走把 Ronkin 函数与 zeta 函数联系起来 | 未涉及非厄米 GBZ；摘要未给出数值精度/复杂度 | 3 | include | 给出 Ronkin 函数的一种可计算对应（zeta 函数），可作求值途径参考 | classic | no | Q2 |
| math-ph/0311005 | Dimers and Amoebae | Richard Kenyon | 2003 | math-ph 预印本 | theoretical | 由 Kasteleyn 算子谱曲线导出表面张力与局部 Gibbs 测度显式公式：表面张力 = 谱曲线 Ronkin 函数的 Legendre 对偶，amoeba = 相图 | 平面双周期二分图 dimer 模型 | 证明 dimer 谱曲线必为 Harnack 曲线，并给出表面张力与 Ronkin 函数的显式关系 | 无算法实现与数值基准；对象为 dimer 统计力学 | 3 | include | amoeba/Ronkin 与谱曲线联系的开创性文献，提供 Ronkin 计算的对偶途径 | classic | no | Q2 |
| 2407.10166 | General theory for infernal points in non-Hermitian systems | Shu-Xuan Wang | 2024 | Phys. Rev. B 110, L201104, DOI 10.1103/PhysRevB.110.L201104 | theoretical | 基于非布洛赫能带论与 amoeba 表述建立 1D 与高维开边界体系中 infernal point（宏观本征态聚并）的判据，并解释极端局域化 | 1D 及高维开边界非厄米模型 | 给出任意维度 infernal point 存在判据与波函数极端局域化机制 | 摘要未给出具体数值实现/精度 | 3 | include | 给出基于 amoeba 表述的判据，属 amoeba 判据的极端退化（失效）情形 | recent3 | yes | Q1 |
| 2310.18627 | Constraints of internal symmetry on the non-Hermitian skin effect and bidirectional skin effect under the action of the Hermitian conjugate of time-reversal symmetry | Shu-Xuan Wang | 2023 | Phys. Rev. B 109, L081108, DOI 10.1103/PhysRevB.109.L081108 | theoretical | 采用 Amoeba 表述建立内部对称性与 NHSE 行为的一般对应 | 任意维度非厄米体系（含时间反演† 对称） | 发现时间反演† 对称下本征态可同时局域在相对边界（双向趋肤效应），超出 Amoeba 表述 | 摘要未给出 amoeba/Ronkin 具体计算步骤；无基准 | 3 | unclear | Amoeba 表述是核心但摘要只见对称性对应，未说明"如何计算"；其 beyond-Amoeba 结论是有价值的失效模式线索 | recent3 | no | Q1 |
| 1310.8472 | Amoebas, Ronkin function and Monge-Ampère measures of algebraic curves with marked points | I. Krichever | 2013 | math.AG 预印本 | theoretical | 对带 punctures 代数曲线上的调和函数对推广 amoeba 与 Ronkin 函数，并证明 M-曲线的极值性质 | 平面代数曲线 / 差分算子谱理论 | 建立 M-曲线极值性质并与正系数差分算子谱理论联系 | 纯形式化结果，摘要未给出可执行的 Ronkin 函数计算方式 | 3 | unclear | 极值性质与差分算子谱理论相关，但摘要无法判断是否含可执行算法 | classic | no | Q2 |
| 2104.04408 | Decimation limits of principal algebraic Z^d-actions | Elizaveta Arzhakova | 2021 | Israel J. Math 265 (2025) 255-299, DOI 10.1007/s11856-024-2676-z | theoretical | 证明 principal Z^d-作用的 decimation 重正化极限存在，且等于 Ronkin 函数 Legendre 对偶的负值 | 整数系数 Laurent 多项式的 principal 代数 Z^d-作用 | 识别 decimation 极限 = −(Ronkin 函数的 Legendre 对偶)，d=2 时与 dimer 表面张力一致 | 纯定理；摘要未报告数值实现、精度或复杂度 | 3 | unclear | 给出 Ronkin 函数（对偶）的一种极限构造，是否可作实用算法摘要层面不明 | classic | no | Q2 |
| 1412.1585 | A Ronkin type function for coamoebas | Petter Johansson | 2014 | math.CV 预印本 | theoretical | 在 coamoeba 框架引入 Ronkin 型函数，并用以刻画 coamoeba 的 shell（toric arrangement）性质 | coamoeba / 复代数簇 | 建立 coamoeba 的 Ronkin 型函数并导出 shell 的性质 | 定义性/形式化工作，摘要未涉及数值计算 | 2 | unclear | 属 Ronkin 函数的正式推广，但摘要未说明"如何计算"，需全文判断 | classic | no | Q2 |
| 1608.06077 | Geometry of generalized amoebas | Yury Eliyashev | 2016 | math.AG 预印本 | theoretical | 把 Krichever 的广义 amoeba/Ronkin 函数推广到高维，并平移已有几何结果 | 广义 amoeba 与 Ronkin 函数（高维） | 给出广义 amoeba 的高维几何结果 | 纯几何定理，摘要未给出计算途径 | 2 | unclear | 广义 Ronkin 函数的形式化几何结果，是否含可执行计算摘要不明 | classic | no | Q2 |
| 2212.03173 | Generalized amoebas for subvarieties of GL_n(C) | Rémi Delloque | 2022 | math.AG 预印本 | theoretical | 把 amoeba 概念推广到 GL_n(C) 子簇（矩阵 amoeba），推广 Ronkin 函数并用 Newton 多面体刻画渐近方向 | 矩阵 amoeba（GL_n(C) 子簇） | 证明矩阵 amoeba 闭、超曲面补集分支凸，并部分推广 amoeba 收敛到 tropical variety | 摘要未涉及数值实现 | 2 | unclear | 涉及 Newton 多面体/tropical 方向（与算法相关），但摘要无计算方法 | classic | no | Q2 |
| 2603.21116 | Solid Amoebas of Maximally Sparse Polynomials | Mounir Nisse | 2026 | math.AG 预印本 | theoretical | 通过分析 Ronkin 函数线性域在 tropical 退化下的稳定性证明 Passare–Rullgård 猜想，并给出 Newton 细分与 tropical 细分的一致性稳定性结果 | 极大稀疏 Laurent 多项式（Newton 多面体顶点支撑） | 证明极大稀疏多项式的 amoeba 必为 solid（补集分支数 = Newton 多面体顶点数） | 纯证明，摘要未提供算法或数值实现 | 2 | unclear | 涉及 Ronkin 线性域/Newton 细分等算法相关结构，但为形式化证明 | recent3 | no | Q2 |
| math/0311062 | Planar dimers and Harnack curves | Richard Kenyon | 2003 | math.AG 预印本 | theoretical | 证明 Harnack 曲线均为某 dimer 模型的谱曲线，用 amoeba 洞面积/触须间距作全局坐标，并刻划 Ronkin 函数下体积的极小化 | 平面 dimer 模型 / Harnack 曲线 | 给出 Harnack 曲线空间与闭八分体同胚，并证明属零 Harnack 曲线是 Ronkin 函数下体积的极小者 | 形式化，无算法/数值基准 | 2 | unclear | 提供 Ronkin 函数的变分极值刻划（与 Ronkin 极小相关），但摘要无计算方法 | classic | no | Q2 |

## EXCLUDED
- 0804.1870 | 弦网络/超引力解中 amoeba 仅作 tropical 曲线比喻、Ronkin 函数识别为 Kähler 势，无任何计算方法
- hep-th/0601233 | 平面分划/5D SYM 中 amoeba 与 Ronkin 函数仅作 limit shape 与超椭圆曲线的几何中介，无算法
- 0902.3996 | 晶体融化/CY 几何中 Ronkin 函数与全纯 3-form 的对应，纯几何识别，无计算
- 1711.00710 | 环簇超曲面高度公式，Ronkin 函数仅作积分项出现，无算法/数值实现
- 2412.16308 | 完全交高度极限公式，Ronkin 函数仅作积分项，无算法
- 0805.1194 | 涡旋/瞬子与 amoeba-tropical 几何的对应字典，无计算方法

## BibTeX — included (14)
```bibtex
@misc{wang2022amoeba,
  title={Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions},
  author={Wang, Hong-Yi and Song, Fei and Wang, Zhong},
  year={2022}, eprint={2212.11743}, archivePrefix={arXiv}, primaryClass={cond-mat.mes-hall},
  note={Phys. Rev. X 14, 021011 (2024); DOI 10.1103/PhysRevX.14.021011}
}
@misc{kaneshiro2025symplectic,
  title={Symplectic-Amoeba formulation of the non-Bloch band theory for one-dimensional two-band systems},
  author={Kaneshiro, Shin and Peters, Robert},
  year={2025}, eprint={2502.17931}, archivePrefix={arXiv}, primaryClass={cond-mat.mes-hall}
}
@misc{xiong2024nhse,
  title={Non-Hermitian skin effect in arbitrary dimensions: non-Bloch band theory and classification},
  author={Xiong, Yuncheng and Xing, Ze-Yu and Hu, Haiping},
  year={2024}, eprint={2407.01296}, archivePrefix={arXiv}, primaryClass={cond-mat.mes-hall}
}
@misc{kaneshiro2025wienerhopf,
  title={Wiener-Hopf factorization and non-Hermitian topology for Amoeba formulation in one-dimensional multiband systems},
  author={Kaneshiro, Shin and Peters, Robert},
  year={2025}, eprint={2511.11349}, archivePrefix={arXiv}, primaryClass={cond-mat.mes-hall},
  note={Phys. Rev. Research 8, 013292 (2026); DOI 10.1103/s43l-h6z6}
}
@misc{lundqvist2013explicit,
  title={An explicit calculation of the Ronkin function},
  author={Lundqvist, Johannes},
  year={2013}, eprint={1306.6249}, archivePrefix={arXiv}, primaryClass={math.CV}
}
@misc{gu2026longrange,
  title={How Long-Range Tails Reshape Non-Hermitian Spectra},
  author={Gu, Ding and Fu, Zhanpeng and Hu, Yu-Min and Wang, Zhong},
  year={2026}, eprint={2608.28577}, archivePrefix={arXiv}, primaryClass={quant-ph}
}
@misc{xu2026domainwall,
  title={Spectral Topology and Non-Bloch Band Theory for Domain-Wall Systems},
  author={Xu, Mingtao and Wang, Rui and Deng, Tian-Shu and Yi, Wei},
  year={2026}, eprint={2607.22976}, archivePrefix={arXiv}, primaryClass={cond-mat.mes-hall}
}
@misc{bobenko2024dimers,
  title={Dimers and M-Curves},
  author={Bobenko, Alexander I. and Bobenko, Nikolai and Suris, Yuri B.},
  year={2024}, eprint={2402.08798}, archivePrefix={arXiv}, primaryClass={math-ph}
}
@misc{chen2022mahler,
  title={Mahler Measuring the Genetic Code of Amoebae},
  author={Chen, Siqi and He, Yang-Hui and Hirst, Edward and Nestor, Andrew and Zahabi, Ali},
  year={2022}, eprint={2212.06553}, archivePrefix={arXiv}, primaryClass={hep-th}
}
@misc{yan2025hsg12m,
  title={HSG-12M: A Large-Scale Benchmark of Spatial Multigraphs from the Energy Spectra of Non-Hermitian Crystals},
  author={Yan, Xianquan and Akg{\"u}n, Hakan and Kawaguchi, Kenji and Loh, N. Duane and Lee, Ching Hua},
  year={2025}, eprint={2506.08618}, archivePrefix={arXiv}, primaryClass={cs.LG},
  note={The Fourteenth International Conference on Learning Representations (ICLR 2026)}
}
@misc{yang2024anatomy,
  title={Anatomy of Higher-Order Non-Hermitian Skin and Boundary Modes},
  author={Yang, Fan and Bergholtz, Emil J.},
  year={2024}, eprint={2405.03750}, archivePrefix={arXiv}, primaryClass={cond-mat.mes-hall},
  note={Phys. Rev. Research 7, 023233 (2025); DOI 10.1103/PhysRevResearch.7.023233}
}
@misc{komatsu2022ronkinzeta,
  title={Ronkin/Zeta Correspondence},
  author={Komatsu, Takashi and Konno, Norio and Sato, Iwao and Sato, Kohei},
  year={2022}, eprint={2212.13704}, archivePrefix={arXiv}, primaryClass={math-ph}
}
@misc{kenyon2003dimersamoebae,
  title={Dimers and Amoebae},
  author={Kenyon, Richard and Okounkov, Andrei and Sheffield, Scott},
  year={2003}, eprint={math-ph/0311005}, archivePrefix={arXiv}, primaryClass={math-ph}
}
@misc{wang2024infernal,
  title={General theory for infernal points in non-Hermitian systems},
  author={Wang, Shu-Xuan and Yan, Zhongbo},
  year={2024}, eprint={2407.10166}, archivePrefix={arXiv}, primaryClass={cond-mat.mes-hall},
  note={Phys. Rev. B 110, L201104; DOI 10.1103/PhysRevB.110.L201104}
}
```
