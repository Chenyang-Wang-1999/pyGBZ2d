# amoeba GBZ 计算调研 · 检索日志（可复现）

日期：2026-09-24 · 工具：arXiv API（`export.arxiv.org`，注意 `http://` 会 302 到 `https://`，**必须直接请求 https**）+ Semantic Scholar Graph API + web 检索
限速：arXiv 要求 ≤1 请求/3 秒

研究问题：现有文献中，amoeba GBZ（及其等价的 Ronkin 函数 / 谱势判据）在**算法层面**是如何实际计算的？

范围（用户 2026-09-24 确认）：① 非厄米 GBZ/amoeba 侧 ② 计算代数几何 + Toeplitz/Szegő 数学侧 ③ 公开代码与软件包
时间窗：不限，分 `recent3`（2023-09 之后）/ `classic` 两档标注
深度：核心 6 篇抓 arXiv 全文精读，其余标题/摘要级

---

## 1. arXiv API 查询

| # | 批次 | search_query（原文） | sortBy | totalResults | 取回 |
|---|---|---|---|---|---|
| A1 | A | `abs:amoeba AND (abs:"non-Hermitian" OR abs:"skin effect" OR abs:"Brillouin")` | relevance | 8 | 8 |
| A2 | A | `abs:"Ronkin function"` | relevance | 21 | 21 |
| A3 | A | `abs:"spectral potential" AND (abs:"non-Hermitian" OR abs:"skin" OR abs:"Toeplitz" OR abs:"Brillouin" OR abs:"lattice")` | relevance | 3 | 3 |
| B1 | B | `abs:"non-Bloch" AND (abs:numerical OR abs:algorithm OR abs:computation)` | submittedDate↓ | 19 | 19 |
| B2 | B | `abs:"generalized Brillouin zone" AND (abs:"two-dimensional" OR abs:"two dimensions" OR abs:"higher dimensions")` | relevance | 16 | 16 |
| B3 | B | `all:"Szego limit theorem"` | relevance | 10 | 10 |
| B4 | B | `abs:"Toeplitz determinant" AND (abs:"non-Hermitian" OR abs:"skin" OR abs:"spectral potential" OR abs:"Widom")` | relevance | 7 | 7 |
| B5 | B | `abs:"spectral potential"` | relevance | 21 | 21 |
| C1 | C | `ti:amoeba` | relevance | 107 | 50 |
| C2 | C | `abs:amoeba AND (abs:"Newton polytope" OR abs:"Monge-Ampere" OR abs:lopsided OR abs:tropical OR abs:"Ronkin")` | relevance | 61 | 50 |
| C3 | C | `abs:"numerical algebraic geometry"` | relevance | 92 | 50 |
| C4 | C | `abs:"amoeba" AND (abs:computation OR abs:computing OR abs:algorithm OR abs:visualization)` | relevance | 70 | 50 |
| D1 | D（主上下文） | `ti:"lopsided" OR ti:amoebas OR abs:coamoebas` | relevance | 212 | 40 |

### 批次统计

| 批次 | raw | unique（去版本号） | include | unclear | exclude |
|---|---|---|---|---|---|
| A（非厄米 amoeba/Ronkin/谱势） | 32 | 28 | 14 | 8 | 6 |
| B（2D GBZ 数值/Szegő–Toeplitz） | 73 | **71** | **14** | **15** | **13**（+29 条同名词异域索引） |
| C（数学侧 amoeba 计算） | 200 | 151 | 24 | 10 | 117 |
| D（主上下文） | 40 | 待合并 | — | — | — |
| **三批合计（未跨批去重）** | **305** | **250** | **52** | **33** | **136** |

**计数勘误留痕**：批次 B 首轮曾报 `unique=40 / include=9 / unclear=23 / exclude=8`，系心算去重漏计 24 条（漏收 Q1 6 条、Q2 3 条、Q3 5 条、Q4 3 条、Q5 7 条），其中 **`2407.01296` 被写入 hit_by 却未做判定**（该篇恰是最该精读的对象之一）。上表为二次逐行点检后的权威值；该错误由批次 B 的 agent 主动披露，此处留痕。

**⚠ 上表为各批 agent 的自报值，与文件实际内容存在三处不一致**（由独立合并 agent 机械逐行点检发现）：批次 B 自报 exclude 13，§2 主表实为 34 个唯一 ID，另有 8 个 unique ID 无法定位；批次 C 自报 unique 151 / include 24，而 manifest 实有 162 行、§2 主表实为 21 include + 10 unclear + 1 exclude；我先前列的"跨批重复 14 个"中**仅 8 个真正跨批**，另发现 6 个未列出的跨批关系（批次 C 判 exclude 故不产生双行）。**权威合并结果以 `screening-log.csv` 为准（72 条 = include 41 + unclear 31）**，详见 `report.md` §1.2 的口径披露表。此差异不影响 §6 结论。

**真正的跨批重复（8 条，出现在 include/unclear 集内）**：`2212.11743`(A,B)、`2407.01296`(A,B)、`2511.11349`(A,B,C)、`2608.28577`(A,B)、`2502.17931`(A,C)、`2212.06553`(A,C)、`math-ph/0311005`(A,C)、`2603.21116`(A,C)。
**另 6 个跨批关系因批次 C 判 exclude 故未产生双行**：`2212.13704`、`1310.8472`、`1412.1585`、`1608.06077`、`2212.03173`、`math/0311062`（已在 `screening-log.csv` 的 rationale 列标注 `[跨批判定不一致：…]`）。

**查询精度观察**
- `abs:"spectral potential"`（B5，21 条）：与 Q1/Q4 仅重叠 2 篇，净新增 19 篇中 **18 篇是同名词异域**（DFT 谱势、动力系统/概率谱势、组合优化/金融谱势），仅 `2506.08618`（HSG-12M）属非厄米但限 1D。**该查询精度极低，建议改用共现约束。**
- **`abs:"Ronkin function"` 是批次 B 五条查询的检索缺口**（一条都没覆盖）——但它在批次 A 的 Q2 中被单独覆盖（21 条），故整体未漏。
- 主体信号：`abs:amoeba ∩ (非厄米 OR 趋肤 OR 布里渊)` = **8 条**；`abs:"spectral potential" ∩ 非厄米相关词` = **3 条**。该交叉区间在 arXiv 摘要层是两位数量级的小集合。

---

## 2. 其他数据源

| 源 | 查询 | 命中 | 用途 |
|---|---|---|---|
| Semantic Scholar Graph API | `arXiv:2212.11743` 的引用列表（fields=title,year,abstract,externalIds,venue, limit=100，`next=100`） | 100（还有下一页） | 反向定位"还有谁在用 amoeba 表述" |
| web 检索 | Theobald「Computing amoebas」、Purbhoo Nullstellensatz、Bogdanov lopsided、SINGULAR/SAGE 实现、phcpy/PHCpack、GitHub/PyPI 代码 | 见 `report.md` §数学侧 | 经典数学文献（多在 arXiv 之前）+ 软件 |
| 本地 PDF 库文件名检索 | `D:\information-base\library`（1713 篇） | 见 `sources-inventory.md` | 仓内已有资产清点 |

**已知抓取限制**
- `web_fetch` 不支持 `application/pdf`（返回 `unsupported content type`），因此 Theobald 2002 等无 arXiv 版本的经典文献**正文未读**，仅据检索元数据登记。
- 本地 PDF 库无正文提取能力，仅文件名级检索。
- arXiv `http://` 端点会 302 到 `https://` 且不被自动跟随，需直接使用 https。
- **arXiv HTML 长页会在约 45–65 KB 处硬截断**（fragment / query 参数均无效），导致多篇论文只能读到前半部分。

**绕行方法（本次验证有效，建议复用）**
- **microlink 分节代理**：可逐节取 arXiv HTML 原文，绕过整页截断。
  ```
  https://api.microlink.io/?url=<URL编码后的 arxiv html 地址>&meta=false&data.x.selector=<CSS选择器>&data.x.type=text
  ```
  selector 用分节 id（`%23S4` 即 `#S4`）；**含点号的 id 必须用属性选择器**，如 `%5Bid%3D%22S7.T2%22%5D`（即 `[id="S7.T2"]`，Table 2 就是这样取到的）。
  ⚠ **microlink 免费额度按日限流**，用尽后返回 `{"code":"ERATE"}` —— 请把要取的节一次性列好再取。
- 备选服务端代理：W3C `https://www.w3.org/services/html2txt?url=<encoded>` 可把截断点从 §2 推后到 §4 末，但自身也受 ~65 KB 上限，且对 AMS/zbMATH（Cloudflare）与 archive.org 无效。
- 效果：据此完成了 `2407.01296` 全文（§I–§X + Appendix A/B/C）、`1101.4114` 全文、`1608.08663` 绝大部分全文的逐字精读。

---

## 3. 未覆盖 / 建议下一轮扩词

批次 A/B 的查询饱和度偏低（A1=8、A3=3、B4=7），说明该方向在 arXiv 摘要层的文献量本身很小（这本身是有意义的信号）。若要进一步扩覆盖，建议：

- `abs:"amoeba"` 单独全量检索（配合 cat:cond-mat.*）
- `abs:"generalized Brillouin zone" AND abs:(algorithm OR numerical)`
- `abs:"Szegő" OR abs:"Szego"` 配合 `cat:cond-mat.*`（注意 A/B 批次已验证：`all:"Szego limit theorem"` 只返回纯数学，物理侧应用不被该式覆盖）
- `abs:"winding number" AND abs:"non-Bloch"`
- 反向引用检索：2412.14912（递推法）、2311.16868（渐近 GBZ）的引用列表
