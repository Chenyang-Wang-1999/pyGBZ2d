# 算法实质卡片 — 2212.11743（定稿）

精读覆盖度：Abstract + Sec. I, II, III, IV(IV.1/IV.2 图注/IV.3), V, VI, IX = 全文正文主体。
用 microlink 分节代理（CSS selector 按节取文本）绕过 fetch 工具 ~40k 字符截断。
未覆盖：Sec. VII（非布洛赫拓扑）、Sec. VIII（谱不等式）、Appendix A/B，以及 Table I 的内容（表在哪个小节未能确认；
已确认存在的只有 arXiv 元数据里的 "1 table" 与章节标题）。

---

## 2212.11743 Amoeba Formulation of Non-Bloch Band Theory in Arbitrary Dimensions

- **元数据**：Hong-Yi Wang, Fei Song, Zhong Wang（清华大学高等研究院）；arXiv 2022-12-22（v3 2024-04-30）；
  Phys. Rev. X **14**, 021011 (2024)；21 页 / 11 图 / 1 表；类型：理论（含数值算例，非数值方法论文）

- **计算对象**：同时计算四个对象，且它们是**同一套对象的不同视角**：
  1. **amoeba**（log|β| 空间中的零集）：$\mathcal A_f=\{\log|\bm\beta|: f(\bm\beta)=0\}\subset\mathbb R^d$，取 $f=\det[E-h(\bm\beta)]$（Eq. 14）。
  2. **Ronkin 函数**：$R_f(\bm\mu)=\int_{T^d}\left(\frac{d\bm\theta}{2\pi}\right)^d\log|f(e^{\bm\mu+i\bm\theta})|$，$T^d=[0,2\pi]^d$（Eq. 15）。
  3. **谱势 / 库仑势**：$\Phi(E)=\phi(E)=\min_{\bm\mu}R_{\det(E-h)}(\bm\mu)$（Eq. 24–25），即 Ronkin 函数的**极小值**。
  4. **态密度 DOS**：$\rho(E)=\frac{1}{2\pi}\Delta\phi(E)$（Eq. 26）。
  另外给出**GBZ**（本征态指数行为）与**能谱支撑集**（洞闭合判据，Eq. 37）。

  **普适性（amoebic spectrum，Sec. IV.2 原文）**："To summarize, there exists a geometry-independent universal spectrum that
  can be calculated from the amoeba and Ronkin function. By nature, it can be called the 'amoebic spectrum.' The DOS of an OBC
  system with a generic shape always approaches the universal DOS in the large-size limit." 并补充："When the DOS of an OBC
  system with a certain (nongeneric) shape appears to deviate from the universal DOS, this deviation can be eliminated by
  adding a small random local perturbation."（注意：**没有**"随机几何"图；随机性只以"边界格点上的随机在位势"形式出现。）

- **输入假设**：一般 $d$ 维平移不变（体态）非厄米紧束缚模型；Bloch 哈密顿量 $h(e^{i\bm k})=\sum_{\bm n}t_{\bm n}e^{i\bm k\cdot\bm n}$，
  $\bm\beta=e^{\bm\mu+i\bm k}$ 为**矩阵值 Laurent 多项式**；**有限跃迁范围**（1D 设 $t_n=0,\,|n|>n_c$）；
  胞内自由度（带数）任意 $M$。数值算例用的模型：
  (a) 1D 非厄米 SSH（2 带，四阶特征方程），$t_1=t_2=1,\ t_3=0.7,\ \gamma=4/3$，链长 $L=300$；
  (b) **单带** 2D 模型 Eq. (12)：$h=t(\beta_x+\beta_x^{-1}+\beta_y+\beta_y^{-1})+t'(\beta_x+\beta_x^{-1})(\beta_y+\beta_y^{-1})+\gamma(\beta_x-\beta_x^{-1}+\beta_y-\beta_y^{-1})$，
  $t=1,\ t'=0.5,\ \gamma=0.2$，(ii) 圆盘 OBC 直径 $L=140$、插图环面 PBC；
  (c) Sec. VII 的非厄米 Chern 带模型（$2\times2$ 矩阵，Fig. 2(b)，Eq. (VII.13)）。

- **算法步骤**（按"可复现"程度整理；原文并未给编号算法框，以下是我从 Eq. 14–41 与各节文字还原）：
  1. 写出特征多项式并展开：$\det[E-h(\bm\beta)]=a_{-M}(E)\beta^{-M}+\cdots+a_N(E)\beta^N$，根按模排序
     $|\beta_1(E)|\le\cdots\le|\beta_{M+N}(E)|$（1D 与 2D 皆然）。
  2. 对给定 $E$，构造 amoeba（Eq. 14）与 Ronkin 函数（Eq. 15）。$d=2$ 时零集是实 2 维，可用 $\beta_x$（或 $\beta_y$）局部参数化。
  3. 判定 $E$ 是否在 OBC 谱内——**看中心洞**：
     - 中心洞**存在** → $E$ **不在** OBC 谱内（$\rho(E)=0$）；
     - 中心洞**不存在/闭合** → $E$ **在** OBC 谱内。
     谱边界 = 中心洞闭合点（Fig. 6：沿实轴降低能量，洞在 $E_t\approx5.95995$ 闭合 = 带顶）。
  4. 求 $\phi(E)=\min_{\bm\mu}R_{\det(E-h)}(\bm\mu)$，再由 $\rho(E)=\Delta\phi(E)/2\pi$ 得 DOS。
     关键结构事实（使极小可判）：**$R_f$ 在 amoeba 补集的每个连通分支（hole）上严格线性、处处凸**，
     梯度即该 hole 的整数序 $\nu_j$；故中心洞存在时极小值**落在中心洞上（一片平台）**，
     否则极小值**落在 amoeba 内部一个单点**。
  5. GBZ：原文定义（Eq. 41）
     $$\det[E-h(\bm\beta)]=0,\qquad (\log|\beta_1|,\dots,\log|\beta_d|)=\bm\mu_{\min}(E).$$
     即"**位置由 Ronkin 极小点决定的那些 $\bm\beta$**"。本征态渐近形式
     $\psi_E(\bm x)=\sum_{\bm k}c_{\bm k}\exp[(i\bm k+\bm\mu_{\min}(E))\cdot\bm x]$（Eq. 39），复波矢 $\tilde{\bm k}=\bm k-i\bm\mu_{\min}$（Eq. 40）。
     实用上把 GBZ 参数化为 $\bm k$ 的函数 $\bm\mu_{\min}(\bm k)$，"该向量函数是 GBZ 的一个完整表示"。
     因为中心洞不存在时极小点唯一（凸性保证），此定义无歧义。
  6. （推导层）用 Toeplitz 矩阵 + Szegő 极限定理把 $\frac1N\log|\det\mathcal T[E-h(\bm\beta)]|$ 化为 torus 积分 $+O(L^{-1})$
     （Eq. 29、31、32）；当 $\bm\mu$ 取在 $\det(E-h)$ 的中心洞内时 $\det[E-h(\bm\beta)]\neq0$、
     沿任一 $k_j$ 圆的绕数为 0，于是可套用 Szegő（Eq. 32 的推导）。

- **精度控制**：
  - **可证部分**：定理（Eq. 37）$\rho(E)=0,\ E\notin\Lambda$（$\Lambda$ = amoeba 无中心洞的 $E$ 集合）给出**严格证明**，
    证明只用"原始版 Szegő 定理"（Eq. 29），不用推广版（Eq. 34）；DOS 推导给出 $O(L^{-1})$ 余项（Eq. 32）。
    Sec. VIII 声称"从 amoeba 途径证明若干 OBC/PBC 谱不等式"；Appendix A 给 Szegő 极限定理的简短证明。
  - **关键单向性警告（务必标注）**：Sec. V 的定理是**单向**的——只证明了
    "$E\notin\Lambda\Rightarrow\rho(E)=0$"（即**谱 ⊆ $\Lambda$**），**没有** if-and-only-if，
    也没有反向陈述；Sec. II.2 的原始表述亦明确对冲："could be an **indicator**"。
    反向方向（谱边界 = $\Lambda$ 边界）在正文中依赖 **Eq. (34) 这一猜想**：
    "if we assume the validity of the conjecture Eq. (34) and therefore Eq. (24), we have $\rho(E)=\Delta\phi(E)/2\pi$,
    which is generally nonzero in $\Lambda$." → 即"洞闭合 ⇔ 在谱内"的充分性未证，是**猜想**。
  - **v1 序言里对 Sec. IV 的预告（原文）**："**A rigorous theorem is proved as the general basis for DOS calculations.**
    We also demonstrate, despite the geometry-dependent NHSE, the existence of a universal spectrum to which the OBC
    spectrum under a generic geometry converges."（v3 把同一句改成"We show in a theorem that the energy spectrum can be
    obtained from the shape of the amoeba"。**注意**：这句是 v1 **Introduction 的预告**，不是 Sec. IV 的定理陈述；
    曾有子代理因"公开搜索索引查不到该短语"而误判为不可靠——实为索引未收录 arXiv v1 HTML 所致，已核实为原文。）
  - **作者自陈的范围**（Sec. IX 原文）："not all aspects of this work are mathematically rigorous. Although numerical
    evidence is supplied whenever a mathematically strict derivation is unavailable, a fully rigorous proof of all our
    main results is of course desirable." → 即部分核心命题（尤其极小化的严格性）依赖数值证据。
    **Appendix A 被作者自己描述为 "a heuristic proof of Szegő's limit theorem"**（启发式证明，非严格）。
  - **推导中出现的定量界（Toeplitz 乘积估计）**：$\|\mathcal T[\sigma_1\sigma_2]-\mathcal T[\sigma_1]\mathcal T[\sigma_2]\|_1=O(L^{d-1})$，
    以及被截断的第二条 $\|\mathcal T[\sigma]^{-1}-\mathcal T[\sigma^{-1}]\|_1\le\|\mathcal T[\sigma]\|\cdots$。
    这是全文唯一可读到的"阶"型误差陈述（Sec. VIII 的谱不等式原文本次未取得）。
  - 全文**未给出**数值容差、收敛判据、误差界或极小化停止条件（极小化在文中是解析/构造性的，未描述数值求解器细节）。

- **复杂度与规模**：
  - 实际算过的最大模型：1D SSH $L=300$；2D 单带模型圆盘 $L=140$（+ 环面 PBC 插图）。
    Sec. VI 的 GBZ 验证用圆盘 $L=400$（固定 $L$ 比对衰减率）与**变尺寸序列**，$E=1,3,5$；
    边界拟合时丢弃距边界 20 格以内的点。
  - Sec. V 的带顶/带底外推：用 $L\in[80,240]$ 的数据线性外推 $L\to\infty$，误差棒为 95% 置信区间；
    有限尺寸修正 $\Delta E_t\propto L^{-2}$（带顶）、$\Delta E_b\propto L^{-1}$（带底）。
  - **文中未给出**：$\bm\mu$ 网格分辨率、能量采样点数、参数扫描规模、复杂度分析、或与实空间对角化的耗时对比。
    成本只有定性表述："不需要对角化随尺寸增长的大矩阵"、"hot 上避开有限尺寸误差"。

- **验证基准**：  - 1D：已知 GBZ 判据 $|\beta_M(E)|=|\beta_{M+1}(E)|$ 的逐点核对（$E_1=-1+0.3i$ 不在谱：$\mu_2\neq\mu_3$；
    $E_2=-1$ 在谱：$\mu_2=\mu_3$），Fig. 1(a)(b)(c)；OBC 与 PBC 实空间对角化对照。
  - 2D：OBC 圆盘 vs PBC 环面对角化（Fig. 1(d)）；洞存在/不存在两个能量点 $E_3=-2.5$（有洞，谱外）/ $E_4=-1$（无洞，谱内）
    （Fig. 1(e)(f)）。
  - DOS：Ronkin 函数算出的 DOS（Fig. 4(a)，由 Eq. 26）对比实空间对角化 DOS——**方格子 $L=130$**（边界加 $[-0.5,0.5]$ 均匀
    随机势以抑制边界态）与**圆盘 $L=140$**（不加随机势）；两组 DOS 都从库仑势经 $\rho=\Delta\Phi/2\pi$ 得到以利比较（Fig. 4）。
  - 带顶/带底：amoeba 洞闭合值 vs 实空间对角化随 $L$ 外推值（Fig. 7），"strikingly close"/"excellent agreement"。
  - GBZ：$\mu_{\min}$ vs 格林函数指数衰减率 $\log\langle\bm x|(E-H)^{-1}|\bm 0\rangle\sim\mu_x x+\mu_y y$ 的线性拟合
    （Fig. 8(c)(d)）；作者说明用格林函数而非本征态本身来验证是"计算上更准、更省"，能算更大的 $L$。
  - **注意**：边界随机势是**数值工具**（防止边界态污染 DOS），不是"对照物"。

- **失效模式/限制**（作者明说或全文结构上明示）：
  1. "not all aspects of this work are mathematically rigorous"（Sec. IX）——数值证据替代严格推导。
  2. 高维 GBZ 的"不可直接计算"：$d$ 维 GBZ 是 $\bm\beta$ 空间（实维 $2d$）中的 $d$ 维子流形；
     对 $d\ge2$ 它比 1D 曲线难算得多——同组 2026 年后续工作（arXiv:2608.28577，引 [3] 即本文）明确写
     "the higher-dimensional GBZ cannot be calculated directly. Instead, the OBC spectrum is encoded in an amoeba
     formulation"。
  3. 1D 那套做法不能平移到 2D：2D OBC 的边界约束方程数 $\propto L$，"it is difficult to obtain a 2D counterpart of
     Eq. (9)"；依赖系数矩阵低秩的做法在高维"quite intractable"。
  4. **多带**：本文形式上是矩阵值符号（任意胞内自由度），但**所有算例是单带或 $2\times2$**，正文没有多带/简并的处理方案；
     后续文献（arXiv:2511.11349, PRR 8, 013292 (2026)）明确指出该 amoeba 表述"limited to single-band systems"，
     多带时 Ronkin 函数会混合各带贡献、在简并且局域长度不同的情形下"fail even in one dimension"（作者未在本文读到相关自陈）。
  5. **数值细节缺失**：Ronkin 积分的求值策略（解析式？数值 quadrature？被积函数在 amoeba 上对数发散如何处理）、
     $\bm\mu$ 极小化的搜索算法与容差、收敛判据——全文未给（这是"复现"的最大缺口）。
  6. 边界的几何形状未被纳入判据（后来者批评本文"neglects geometric information"，谱与几何无关）；
     本文自己主张的是"amoebic spectrum：任意 generic 几何下 OBC 谱都收敛到的普适谱"。
     值得一提：**v1（arXiv:2212.11743v1，标题为 "Amoeba formulation of the non-Hermitian skin effect in higher dimensions"）
     的原文更直白**："the recent finding of geometry-dependent NHSE in 2D non-Hermitian systems suggests that
     **the GBZ might even not be definable in 2D**"（v1 第一手引文）。
  7. **只处理标准 OBC**（Sec. IX 原文）："it seems that the amoeba approach, as we now understand, naturally corresponds to
     the standard OBC systems … we have focused on the OBC case throughout the present paper"——域壁（domain-wall）等其它
     边界条件不在覆盖范围内（1D GBZ 可推广过去，本文的 amoeba 途径不行）。
  8. 例外点/简并只在 Eq. (44) 这一层被"绕过"（脚注给出 nilpotent 项在围道积分后消失的论证），并非专门处理。
  9. **有限尺寸误差不对称**（Sec. V 原文）："the finite-size correction $\Delta E_t\propto L^{-2}$ for the former …, and
     $\Delta E_b\propto L^{-1}$ for the latter … Therefore, the error of extrapolation is larger for the latter."

- **作者列出的开放问题**（Sec. IX 末段原文）：把该表述用于**自由粒子极限下的开放量子系统**（非厄米 Liouvillian 超算符谱
  决定动力学与弛豫，"Our theory immediately enables calculating the relevant quantities beyond 1D"）；
  以及**多体非厄米系统**——"our amoeba theory may still be a good starting point for including the interaction effects,
  which will be left for future work"。
  7. 边界项被忽略：GBZ 波函数式（Eq. 39）的脚注明确"possible boundary terms are omitted, in the same spirit as in 1D"。

- **代码**：**未公开**。Acknowledgements 完整且只有一句："This work is supported by NSFC under Grant No. 12125405."——
  无致谢名单、**无 code-availability、无 data-availability 声明**；APS 元数据亦无；`/supplemental/` 返回 404。

- **与 Ronkin 极小 / 平均绕数零点判据的关系**：
  本文不是把两条判据并列，而是**用 Ronkin 极小把两者统一起来**：
  - 平均绕数就是 Ronkin 的梯度：$\nu_j=\partial R_f/\partial\mu_j=\mathrm{Re}\int(\frac{d\theta}{2\pi})^d\frac{\partial_{\mu_j}f}{f}$，
    并且 $=\int_{T^{d-1}}\frac{d\bm\theta}{(2\pi)^{d-1}}w_j$，其中 $w_j=\frac{1}{2\pi i}\oint d\theta_j\,\partial_{\theta_j}\log f$
    是沿 $\theta_j$ 圆的绕数（Eq. 16–18）。$\nu$ 在每个 hole 上取常整数"序"；**最多只有一个 hole 的序为 $\nu=(0,\dots,0)$，
    即"中心洞"**。
  - 判据链：中心洞存在 ⟺ Ronkin 极小落在该洞上（平台，$\nu=0$ 区域）⟺ $E$ 在谱外（$\rho=0$，严格证明 Eq. 37）；
    中心洞不存在 ⟺ 极小在 amoeba 内单点 $\bm\mu_{\min}(E)$ ⟺ $E$ 在谱内且本征态按 $e^{\bm\mu_{\min}\cdot\bm x}$ 指数行为。
    1D 时该套退化为已知 GBZ 判据：区间收缩为单值 $\mu$，区间两侧即 $|\beta_M|=|\beta_{M+1}|$。
  - 因此：**与"Ronkin 极小"是同一件事**（$\phi=\min R$）；**与"平均绕数零点"是等价的实现方式**
    （极小点处 $\nabla R=\nu=0$；洞 = 常绕数区域，中心洞 = $\nu=0$）。原文没有写"等价"二字的显式断言，
    而是给出凸性 + 平台/单点二分（Fig. 3 caption）作为两者等价的结构性依据。

- **对"通用 2D 多带 amoeba 数值计算"的贡献**（1–2 句）：
  它把 amoeba 从几何图像变成**可计算判据**：中心洞（= 绕数 $\nu=0$ 的洞）的有无直接判定 $E$ 是否在 OBC 谱内，
  并把谱势归约为 Ronkin 函数 $R_E(\mu_x,\mu_y)$ 的**全局极小**、DOS 归约为其 Laplacian；同时给出 2D GBZ 的可写公式
  （Eq. 41）与本征态指数行为。对多带，真正可复用的是"$R_E$ 凸 + 极小为平台或单点 ⟺ 谱内/谱外"这一二分，
  但论文本身只演示单带/$2\times2$ 情形，且**未公开可复现的数值求值与极小化细节，也未公开代码**。

- **精读覆盖度**：**全文的绝大部分**——Abstract、Sec. I、II、III、IV（IV.1、IV.2 图注、IV.3）、V、VI、IX 均已逐节读取原文；
  缺 Sec. VII（非布洛赫拓扑数值）、Sec. VIII（谱不等式具体形式）、Appendix A/B 正文、Table I 内容（这些不影响上面关于
  "怎么算"的核心字段）。
  另经检索确认：**Sec. VIII "Spectral inequalities" 的每条不等式原文与"已验证/属猜想"的判定状态未能取得**——仅能确认
  Intro 的表述"several useful inequalities on the OBC and PBC spectra **are proved** from the amoeba approach"，
  以及该节围绕**一条**"the spectral inequality"（v1 该节标题为单数）、并有一处 "alternative proof"（§ 提示性片段）。
  Appendix A 仅见到开头一句（"Szegő's limit theorem was originally established for 1D Hermitian Toeplitz matrices [41],
  but thereafter generalized by Widom…"）；Appendix B 只知是 v3 新增（v1 无 Appendix B）。

- **关键原文引用**（≤3 条，便于核对）：
  1. （Sec. VI 开头，GBZ 的定义方式）"In this section, we establish another proposal mentioned at the end of Sec. III, namely,
     **the location of the Ronkin minimum determines the complex momenta and hence the GBZ.**"
  2. （Fig. 3 图注，凸性 + 平台/单点二分，这是"Ronkin 极小 ↔ 绕数零点"的结构性依据）"the Ronkin function is strictly linear
     on each component (hole) of the complement of the amoeba, where the gradient equals the integer index. The Ronkin function
     is always convex. Consequently, **when the central hole exists, the minimum is reached on the central hole; otherwise, the
     minimum is reached at a single point in the amoeba.**"
  3. （Sec. V 定理 Eq. 37 及其证明思路）"$\rho(E)=0,\ E\notin\Lambda$ … Note that this proof makes use only of the original
     version of Szegő's theorem Eq. (29), without invoking the generalized version Eq. (34)."

---

## 附：与"怎么算"相关的旁证（非本文原文，已标注来源）

用于交叉核对"复现本文需要什么"，不作为本文事实：
- **Ronkin 积分的实际求值**：本文未给数值求值细节；同一路线/同组后续工作对 2D 的处理是**直接在环面 $T^2$ 上做数值积分**
  ——$R_E(\mu_x,\mu_y)=\int_{T^2}\frac{d\theta_x d\theta_y}{(2\pi)^2}\log|f(e^{\mu_x+i\theta_x},e^{\mu_y+i\theta_y})|$，
  再对 $(\mu_x,\mu_y)$ 求全局极小（arXiv:2608.28577，2026，Tsinghua 同组；其 2D 基准为方格子 $L=130$/边界无序、
  以及 $g=10^{-4},L=256$ 的扰动算例）。即：**最小二重积分 + 极小化**是这条路线可复现的核心操作，
  而本文对"被积函数在 amoeba 上对数发散时如何处理"没有给出方案。
- **v3 的"显式 Ronkin 公式"线索（未读全）**：由搜索引擎片段可见 v3 有一处
  "We now apply the explicit formula of the Ronkin function to $\det[E-h(\beta)]=a_{-M}(E)\beta^{-M}+\dots$"，
  即 1D 情形下 $R$ 可用系数 $a_{-M}(E)$ 与根显式写出；但**含 $\log|a_{-M}(E)|$ 的完整式子未读到**，
  本卡片不作补全（上述 Eq. (19)/Jensen 显式式子是同一件事的已读版本）。
- **勿混淆**：v3 正文中**没有**"plateau ⇔ $E$ 在谱外"这样的原话；本文原话是
  "the absence (presence) of a hole … could be an **indicator**"（Sec. II.2）+ 定理 Eq. (37)（$\rho(E)=0,\ E\notin\Lambda$）
  + Fig. 3 图注的凸性/平台-单点二分。"plateau" 一词是后续文献（arXiv:2502.17931、2608.28577）对本文结论的转述用词。
- **1D 的闭式 Ronkin 公式**（可直接当算法用）：$R_\sigma(\mu)=\ln|C_E|-p\mu+\sum_{j=1}^{p+q}\max(\mu,\mu_j)$，
  $\mu_j=\ln|\beta_j|$；极小落在 $[\mu_p,\mu_{p+1}]$（即 1D GBZ 条件），在 OBC 谱上该区间收缩为单值 $\mu$
  （arXiv:2511.11349v2，Kaneshiro & Peters, PRR 8, 013292 (2026)，其 1D 小节）。本文 Sec. III 的 Eq. (19)（Jensen 形式
  $\frac1{2\pi}\int_0^{2\pi}d\theta\log|g(Re^{i\theta})|=\log|g(0)|+\sum_k\log|R/z_k|$）是同一件事的 1D 原版。
- **多带的适用性争议**：arXiv:2511.11349 明确写该 amoeba 表述"limited to single-band systems"，多带时 Ronkin 混合各带
  贡献、简并且局域长度不同时"fail even in one dimension"；arXiv:2407.01296（Xiong–Xing–Hu）则把本文表征为
  "neglects geometric information and yields geometry-irrelevant non-Bloch spectra"，并称其 GBZ 是
  "$d$D object embedded in a $2d$D space"。这些是**他人评述**，本文正文未读到相应自陈。
