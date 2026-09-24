# amoeba GBZ 计算方法调研 · 第一轮：仓库内已有文献清点

日期：2026-09-24 · 调研问题：**现有文献中，amoeba GBZ（及其等价的 Ronkin 函数 / 谱势判据）在算法层面是如何实际计算的？**

本文件只记录**仓库内已存在的**文献资产（不含本轮新检索），供后续合成报告引用。

---

## A. 代码仓自带综述的命中（`paper/literature-survey-2D-GBZ/`，2026-08-18，title/abstract 级）

| 编号 | arXiv / 出处 | 对"怎么算"的贡献 | 仓内位置 |
|---|---|---|---|
| **F6** | 2212.11743 → PRX **14**, 021011 (2024) | **方法学奠基**：谱势 + Ronkin 函数 + amoeba；本项目 `amoeba` 模块复刻对象 | `processed-papers.md`；`playground/replication-ZWang.py` |
| **B4** | 2511.11349 → PRR **8**, 013292 (2026) | Wiener–Hopf 分解给多带 amoeba 严格基础；明确指出广义 Szegő 定理形式上任意维适用、**实现只限单带**，多带连 1D 都不清楚 | `report.md` §3B |
| — | 2502.17931（Kaneshiro，2025） | Symplectic-Amoeba：1D 双带的**能带分辨 Ronkin 优化**（简并多带、外推总 Ronkin 函数） | `report.md` 附录 A |
| **A9** | 2407.01296 → **Commun. Phys. (2026)**, DOI 10.1038/s42005-026-02546-2 | 谱势 + 几何自适应非布洛赫能带理论；指出非收敛时谱趋向 amoeba 谱 | `report.md` A9 |
| **A3** | 2405.03750 → PRR **7**, 023233 (2025) | amoeba 理论（bulk separation gap）与 GBZ 路线（surface gap）头对头比较 | `report.md` A3 |
| **B5** | 2407.10166 → PRB **110**, L201104 (2024) | 用 amoeba 表述给任意维 infernal point（谱坍缩）判据 | `report.md` B5 |
| **F8** | 2102.05059 → Nat. Commun. **13**, 2496 (2022) | 面积定理 + 几何依赖趋肤——amoeba 谱的物理前身 | `processed-papers.md` |
| **B10** | 2501.15209 → PRB **111**, 214305 (2025) | Wasserstein 度量提取 (auxiliary) GBZ——替代数值路线 | `report.md` B10 |
| **B11** | 2604.06998 | ML 重建 3D GBZ | `report.md` B11 |
| F7 / A10 | 2210.04412 / 2311.16868 | strip / 渐近 GBZ（对比项，非 amoeba） | `report.md` F7、A10 |

---

## B. 本地 PDF 库命中（`D:\information-base\library`，共 1713 篇）

按文件名检索，与"怎么算 amoeba"直接相关者：

### B1. 非厄米物理侧（amoeba / 谱势路线）

| 文件 | 为什么关键 |
|---|---|
| `非厄米与量子开放\WangZhong@THU\2024_PRX_Wang et al_Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions.pdf` | **F6 全文**（仓内综述只做了摘要级） |
| `非厄米与量子开放\2025-PRB-kaneshiro_and_peters-symplectic-amoeba-formulation...1D two-band...pdf` | 2502.17931 全文 |
| `非厄米与量子开放\HuHaiping@IOP\2024-xiong_et_al.-...arbitrary_dimensions_non-bloch_band_theory_and_classification.pdf` | A9 全文（谱势路线） |
| `非厄米与量子开放\YangZhesen@IOP\2024_Xu et al_Two-dimensional Asymptotic GBZ Theory.pdf` | A10 全文 |
| `非厄米与量子开放\LiuYongchun@THU\2025-wang_et_al.-non-bloch_band_theory_for_2d_geometry-dependent...pdf` | 本组理论文（SGBZ） |
| `非厄米与量子开放\2023_PRB_Sarkar et al_...Haldane model.pdf`、`JinghuiPi@THU\2024_Lin et al_Imaginary Stark Skin Effect.pdf`、`LiuYongchun@THU\2022_PRL_Li et al_Gain-Loss-Induced Hybrid Skin-Topological Effect.pdf` | 候选测试模型 |
| `非厄米与量子开放\2022_Nature Communications_Zhang et al_Universal non-Hermitian skin effect in two and higher dimensions.pdf` | F8 全文 |

### B2. 数学侧（计算 amoeba / Ronkin 的经典基础）——**仓内综述完全未覆盖**

| 文件 | 为什么关键 |
|---|---|
| `非厄米与量子开放\2004_Duke Mathematical Journal_Passare_Rullgård_Amoebas, Monge-Ampère measures, and triangulations of the Newton polytope.pdf` | 计算 amoeba 的经典源头：Newton 多面体三角剖分 + lopsided 近似 |
| `非厄米与量子开放\2019_Complex Analysis and Operator Theory_Ossete Ingoba_A New Insight on Ronkin Functions or Currents.pdf` | Ronkin 函数/流的一般理论 |
| `非厄米与量子开放\` 下 **Toeplitz / truncated Toeplitz 算子文献群**（Berger–Coburn 1986、Reichel–Trefethen 1992、Garcia–Ross 2012、O'Loughlin 2021/2022、de Monvel 1978、Upmeier 1991、Hayashi 2022、Venugopalkrishna 1972 等） | Szegő/谱势路线的数学底座 |
| `mathematics\1996-NumerAlgor-bini-numerical_computation_of_polynomial_zeros_by_means_of_aberth's_method.pdf` | brute-force 求根的直接算法邻居 |
| `电子书\2006_Decker_Computing in algebraic geometry.pdf` | 代数几何计算入门 |

---

## C. 能力边界（诚实说明）

本轮清点对本地 PDF 只做到**文件名清单级检索**。当前工具链只有 UTF-8 文本读取与图片读取，**没有 PDF 正文提取能力**，因此上述 PDF 的正文内容未被阅读，也不据此下任何结论。要拿正文，可行路径是抓同篇的 arXiv 全文（`arxiv.org/html/<id>` 或 `ar5iv.labs.arxiv.org/html/<id>`）。

---

## D. 从已有文献初步看出的路线骨架（待第二轮验证/补全）

1. **Ronkin 极小 / 双绕数零点路线** — F6 及本项目：在 μ 空间找 ∇R_E = 0（或 u₁ = u₂ = 0），用平均绕数判谱成员。**本项目现用路线。**
2. **Szegő / Toeplitz 行列式渐近（谱势）路线** — B4、2502.17931：把 σ_Amoeba 归结为谱势零点集；严格性强，但实现目前只到 1D 单带 / 1D 双带。
3. **几何自适应谱势路线** — A9（现 Commun. Phys. 2026）：对给定几何求谱、态密度与 GBZ，含收敛性/稳定性讨论。
4. **计算代数几何路线** — Passare–Rullgård / lopsided / tropical / Newton 多面体三角剖分：这是"算 amoeba"的数学原义，本地库已有奠基文献，仓内综述未覆盖。
5. **间接/替代路线** — 面积定理（F8）、渐近 GBZ（A10）、Wasserstein（B10）、ML 重建（B11）。

---

## 相关文件

- 本轮检索与筛选：`report.md`、`screening-log.csv`、`processed-papers.md`
- 上一轮（2D GBZ 总览）：`../literature-survey-2D-GBZ/report.md`
