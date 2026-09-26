// Standalone source; no external Typst packages or image assets.
// From the repository root:
// typst compile paper/algorithm-analysis/algorithm-analysis.typ
// Analysis baseline: 2026-09-24, src/ at Git HEAD af595fa.

#set document(title: "pyGBZ2d 算法分析：误差、稳定性、复杂度与验证计划")
#set page(
  paper: "a4",
  margin: (x: 22mm, top: 21mm, bottom: 21mm),
  header: align(right, text(size: 8.5pt, fill: rgb("64748b"))[pyGBZ2d · 算法分析工作记录]),
  footer: context align(center, text(size: 9pt)[#counter(page).display("1 / 1", both: true)]),
)
#set text(font: ("New Computer Modern", "Noto Serif SC", "SimSun"), size: 10.5pt, lang: "zh")
#set par(justify: true, leading: 0.55em, spacing: 0.75em)
#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)")
#show math.equation: set text(font: "New Computer Modern Math")
#show heading: set text(font: ("Noto Sans SC", "Microsoft YaHei"), fill: rgb("213b53"))
#show raw: set text(font: ("Cascadia Mono", "Noto Sans SC"), size: 8.5pt)
#set table(inset: 5pt, stroke: 0.35pt + rgb("cbd5e1"), align: left + horizon)
#show table.cell.where(y: 0): set text(weight: "bold")
#show link: set text(fill: rgb("285c83"))

#let note(title, body) = block(
  width: 100%, breakable: true, inset: 9pt,
  fill: rgb("f1f5f9"), stroke: (left: 2pt + rgb("446b89")),
)[*#title* #body]
#let source(path, body) = link("../../" + path, body)
#let smalltable(..args) = text(size: 9pt, table(..args))
#let lesssim = $≲$

#align(center)[
  #text(size: 23pt, font: ("Noto Sans SC", "Microsoft YaHei"), weight: "bold")[pyGBZ2d 算法分析]
  #v(5pt)
  #text(size: 14pt)[误差、稳定性、复杂度与数值验证计划]
  #v(8pt)
  源码分析工作记录 · 2026 年 9 月 24 日
]
#v(8pt)

#note("文档状态。", [
  本文整理已有源码分析、局部数学推导、文档中已记录的 benchmark，以及尚未实施的数值验证计划。
  未据此修改求解器，也未开展新的数值实验。下文的条件估计不是当前程序已提供的严格误差证书；待验证风险也不等同于已复现的程序错误。
])

分析以 #source("paper/slides/pygbz2d-details.typ")[算法说明的 Typst 源码] 与当前求解器源码为依据。
源码基线为 Git 提交 #raw("af595fa")；算法说明按分析时工作区内容读取。
本文未使用图片或 PDF 作为分析依据。适用对象为固定参考能量下的 SGBZ、amoeba 和共享 continuation 主流程，不包含 experimental band clustering。

核心区分是：*在拓扑结构判定正确时分析连续误差传播；在判定余量不足时分析离散结构改变。*
误差主链为
$
  "系数与求根误差" arrow.r "交点角度误差" arrow.r "平均绕数误差" arrow.r "μ 与 GBZ 误差".
$
轨迹换序、漏交点、错误重根或 continuum 分类、错误谱内外判定，需要另行检查，不能全部归入上述连续误差链。

#outline(title: [目录], depth: 2)
#pagebreak()

= 算法结构与分析对象

== 两种求解路径

两种方法都固定参考能量 $E$，搜索衰减参数，并通过多项式零点重建绕数。差别在于横向半径和绕数求值方式。

#smalltable(
  columns: (0.85fr, 2fr, 2fr),
  table.header([环节], [SGBZ], [amoeba]),
  [外层变量], [$mu_1$，求平均主轴绕数 $W=0$], [$mu_1$，沿内层极小值求 $w_1=0$],
  [横向半径], [随 $theta_1$ 变化的边界两根对数模均值], [内层搜索得到的常数 $mu_2$],
  [零点事件], [两条根轨迹等模], [根轨迹与水平线 $ln abs(beta_2)=mu_2$ 相交],
  [整数绕数], [选取环路，积分 $Im(f'/f)$ 后取整；普通电荷传播], [在角度区间中点重新求根，按圆内根数计数],
  [退化处理], [continuum 双侧探针、MR 事件、线段提取], [continuum 内外层探针、平坦轨迹与线段提取],
)

共享 #raw("ZeroManager") 在 $theta_1$ 方向追踪全部 $beta_2$ 根。一次普通步先作切向预测，再实际求根，以弦距离指派匹配根，并用预测误差调节下一步长度。
普通节点并不是将微分方程数值积分得到的近似根，因此不存在通常 ODE 积分方式的根值累计截断误差；不过错误轨迹身份可以传播到后续节点。
MR 聚类吸附后的节点是例外，应单独记账。

== 输出精度的不同含义

需要分别检查：

+ 特征方程残差：输出是否接近 $f(E,beta_1,beta_2)=0$。
+ GBZ 约束误差：$mu$、等模条件、角度与目标 GBZ 的几何偏差。
+ 拓扑正确性：点数、线数、轨迹闭合置换、MR 成员、电荷、连通关系与谱内外分类。

特征方程的零集通常大于目标 GBZ。很小的方程残差不能单独证明 GBZ 约束满足，也不能证明所有连通分支已被发现。
成功返回空集与失败返回空集必须区分；后者不能用来判断谱外。

== 主要容差控制什么

#smalltable(
  columns: (1.5fr, 0.7fr, 2.2fr),
  table.header([参数], [默认值], [直接控制的对象]),
  [#raw("STEP_ATOL / STEP_RTOL")], [$10^(-12)$ / $10^(-3)$], [预测根与实际根的归一化弦距离；影响网格与匹配],
  [amoeba #raw("CROSSING_XTOL")], [$10^(-12)$], [交点 Brent 求解的角度停止容差],
  [SGBZ #raw("CROSSING_TOL")], [$10^(-10)$], [交点角度精化，以及相近事件分组尺度],
  [#raw("WINDING_ZERO_TOL")], [$10^(-8)$], [平均绕数接近零的停止判据],
  [#raw("CONTINUUM_TOL")], [$10^(-6)$], [对数模的平坦或等模检测],
  [#raw("CONTINUUM_FRAC")], [$0.9$], [SGBZ 等模聚类的网格行投票阈值],
  [MR #raw("CLUSTER_TOL")], [$10^(-4)$], [有限根的弦距离聚类阈值],
  [#raw("CONTINUUM_PERTURB")], [$10^(-4)$], [continuum 双侧探针的基础参数偏移],
  [SGBZ 积分绝对/相对容差], [$10^(-3)$ / $10^(-3)$], [未取整环路积分的数值精度目标],
)

这些量的单位和用途不同，不能比较数值大小后选一个作为“总误差”。
完整参数定义与生效规则见 #source("doc/constants.md")[数值常数说明]。
下文首先固定参数与拓扑结构，建立局部误差估计，再讨论这些前提何时失效。

= 条件性误差传播

== 多项式求根、机器误差与条件数

固定 $E,beta_1$，清除 Laurent 分母后写成
$ p(z)=sum_(k=0)^K a_k z^k, quad z=beta_2. $
对单根 $z_j$，一阶扰动关系为
$ delta z_j approx - frac(delta p(z_j),p'(z_j)). $
若系数扰动满足 $abs(delta a_k)<=eta abs(a_k)$，则
$
  frac(abs(delta z_j),abs(z_j)) lesssim kappa_j eta,
  quad kappa_j := frac(sum_(k=0)^K abs(a_k) abs(z_j)^k, abs(z_j p'(z_j))).
$ <eq-root-condition>
该式要求根为有限非零单根，且扰动足够小，尚未进入重根分裂的非线性尺度。
近重根处 $p'(z_j)$ 变小，前向误差可能被显著放大。

当前 #source("src/pygbz2d/core.py")[#raw("CharPoly.solve_roots_1d")] 调用 #raw("np.roots")。
NumPy 通过伴随矩阵特征值求根 [N1]；不能直接将矩阵特征值算法的稳定性替换成原多项式逐系数相对误差为机器精度的保证。
双精度单位舍入误差约为 $2^(-53)$，但它只描述底层算术尺度，实际前向误差还取决于条件数、系数构造和求值抵消。

建议将误差来源拆成两层：

+ 原始模型参数到双精度系数的误差，包括部分代入、合并同次幂及抵消。
+ 将给定双精度系数视作精确输入后，求根与后续运算产生的误差。

单根可同时记录归一化残差与参考根位置误差：
$
  eta_("res")(z)=frac(abs(p(z)),sum_k abs(a_k) abs(z)^k),
  quad epsilon_(log z)=abs(log z_("num")-log z_("ref")).
$
对数虚部应在匹配后选择一致分支。残差最好用更高精度重新计算，避免残差求值本身的抵消掩盖问题。
该残差描述一个候选根与多项式的相容程度，不自动认证整组根的完整性及共同后向误差。

对有限非零根，
$
  delta mu_2 approx Re frac(delta z,z),
  quad delta theta_2 approx Im frac(delta z,z).
$
根的相对误差由此进入对数模和相位。弦距离靠近球面两极时不能直接替代对数坐标误差；靠近 $0$ 或 $infinity$ 要分别报告度量。

== 交点位置精度与横截性

amoeba 与 SGBZ 分别求解
$
  g_j(theta_1)=ln abs(beta_(2,j)(theta_1))-mu_2=0,
  quad
  g_(a b)(theta_1)=ln abs(beta_(2,a)(theta_1))-ln abs(beta_(2,b)(theta_1))=0.
$
共享解析切向量为
$
  V_j := frac(d log beta_(2,j),d theta_1)
  =-i frac(beta_1 f_(beta_1), beta_(2,j) f_(beta_2)),
  quad g_j'=Re V_j,
  quad g_(a b)'=Re(V_a-V_b).
$
固定 $mu$，在事件已被发现、轨迹身份正确、交点为横截交点的条件下，有局部估计
$
  abs(delta theta_(1,*)) lesssim epsilon_("Brent")
  +frac(epsilon_g,abs(g'(theta_(1,*)))).
$ <eq-crossing-error>
其中 $epsilon_g$ 包含多项式求根和对数模求值误差。
Brent 的停止尺度包括绝对角度容差及相对项；其收敛前提是连续函数与有效异号括号 [S1]。
角度容差不包括函数值误差通过 $1/abs(g')$ 的放大，也不保证发现所有交点。

相切时该式失效。例如 $g(theta) approx a(theta-theta_*)^2$ 受到 $epsilon_g$ 量级的扰动后，角度变化可能为 $O(sqrt(epsilon_g/abs(a)))$，并可能产生两个交点或消去交点。

主轴绕数的分区通常位于 $theta_2$，因此还要传播一次：
$
  abs(delta theta_2) lesssim abs(Im V_j) abs(delta theta_1)
  +abs(Im frac(delta beta_(2,j),beta_(2,j))).
$
第二项为固定角度上的求根误差；外层 $mu$ 误差将在后文单独加入。
仅测量 $theta_1$ 精度会遗漏相位变化很快时的放大。

当前实现还存在两个应单独记录的路径：SGBZ Brent 失败时会警告并退回线性预测；EventGroup 会把相近事件放到代表角度。
前者不再具有 Brent 精度，后者应记录实际组跨度。连续相邻事件按阈值合并时，整个组跨度也可能大于一次合并阈值。
详见 #source("src/pygbz2d/sgbz/pairwise.py")[SGBZ 事件精化与分组] 和 #source("src/pygbz2d/amoeba/ronkin_winding.py")[amoeba 精化]。

== 从交点角度到平均绕数

在正确的循环分区上，两种方法最终均可写成
$
  W=frac(1,2 pi)sum_(k=1)^C u_k(alpha_(k+1)-alpha_k),
$
其中 $u_k$ 是区间内整数绕数，循环展开满足 $alpha_(C+1)=alpha_1+2 pi$。
令事件跳变量为 $q_k=u_k-u_(k-1)$。若事件对应关系、循环顺序与整数绕数保持正确，则
$
  delta W=-frac(1,2 pi)sum_k q_k delta alpha_k,
  quad abs(delta W)<=frac(1,2 pi)sum_k abs(q_k) abs(delta alpha_k).
$ <eq-winding-error>
这一步对分区角度是精确线性的；非线性主要在底层根误差到事件角度的传播中。
它是后续数值实验最适合逐项验证的公式之一。

若某区间的整数绕数错了 $Delta u$，则该区间贡献约为
$ Delta W=frac(Delta u Delta alpha,2 pi). $
该误差可能很大，也可能由于区间很窄而很小。
漏掉一对相反电荷可能仍满足总电荷守恒，因此小绕数误差或电荷总和为零均不认证完整拓扑。

== 环路整数的恢复条件

正确的整数绕数一旦确定，就不再贡献连续小误差；确定这个整数的数值过程仍可能出错。

amoeba 的区间中点根计数需要足够的圆内外判定余量，例如
$
  min_j abs(ln abs(beta_j)-mu_("count"))
  >epsilon_(log abs(beta))+epsilon_(mu_("count")).
$
还需要根列表按正确重数完整给出，并且该区间没有漏掉的奇异角度。

SGBZ 对数值环路积分取整。设真整数绕数为 $w$，积分结果为 $hat(w)$；只有在 $abs(hat(w)-w)<1/2$ 时，取整才恢复正确整数。
当前 #source("src/pygbz2d/sgbz/winding.py")[#raw("get_winding_number")] 仅使用 #raw("quad") 返回的积分值。
给定积分容差本身不认证此不等式；#raw("quad") 的误差输出也是估计量 [S2]。
数值结果很接近某个整数，也不能排除漏积分尖峰后接近了错误整数。

== 从平均绕数到二分变量

SGBZ 离散分支以实际绕数满足容差停止，不以区间宽度单独停止。
设准确解满足 $W(mu_1^*)=0$，且局部 $W'(mu_1)>=m>0$。
如果 $abs(hat(W)-W)<=epsilon_W$，而停止条件为 $abs(hat(W))<=tau_W$，则
$
  abs(hat(mu)_1-mu_1^*)<=frac(tau_W+epsilon_W,m).
$ <eq-mu-sgbz>
因此绕数容差不是 $mu_1$ 的位置容差。平台上解可能不唯一，应改用到解集的距离；跳跃或不可微点则需使用单侧极限与局部结构。

amoeba 还包含内外层误差耦合。
设 Ronkin 函数在严格局部极小值附近二阶可微，Hessian 为 $H$，并且 $H_(22)>0$。
沿内层极小值曲线，外层有效曲率为
$
  kappa_("out")=H_(11)-frac(H_(12)H_(21),H_(22)).
$
若该量为正，记两方向绕数计算误差为 $epsilon_1,epsilon_2$，停止容差为 $tau_1,tau_2$，则局部一阶估计为
$
  abs(delta mu_1) lesssim
  frac(tau_1+epsilon_1+abs(H_(12)/H_(22))(tau_2+epsilon_2),kappa_("out")),
$ <eq-mu-amoeba>
$
  abs(delta mu_2) lesssim
  abs(H_(21)/H_(22)) abs(delta mu_1)+frac(tau_2+epsilon_2,H_(22)).
$
内层误差可以通过耦合项主导最终 $mu_1$ 误差。
这些式子不适用于 continuum、平坦方向或 Hessian 退化情形。
实际停止与复用逻辑见 #source("src/pygbz2d/amoeba/bisect.py")[amoeba 二分] 与 #source("src/pygbz2d/sgbz/sgbz_solver.py")[SGBZ 二分]。

== 最终 GBZ 坐标与插值误差

普通根的总误差还满足
$
  delta log beta_2 approx
  -frac(beta_1 f_(beta_1),beta_2 f_(beta_2))delta mu_1
  +V delta theta_1+frac(delta beta_(2,"solve"),beta_2).
$
式中最后一项为固定参数上的求根误差。它说明即便二分变量足够精确，近重根处的输出坐标仍可能较差。

amoeba 返回离散交点时将根相位投影到指定半径 $exp(mu_2)$；SGBZ 普通点通常使用重新求解的根。
因此应同时检查返回坐标的原方程残差与约束残差，不能仅检查内部求根残差。

在光滑区间，节点与导数精确时，三次 Hermite 插值有 $O(Delta theta^4)$ 的函数误差。
实际插值还包含节点误差与 $O(Delta theta)$ 乘以导数误差的贡献。
非有限导数触发的线性退化，以及 MR 附近低光滑度，都不能直接沿用四阶结论。
线集节点精确也不等于整条曲线表示精确；节点之间需要独立抽样验证。

= 重根、拓扑判定与稳定性

== 重根分裂与增广系统的不同条件数

若 $p(z) approx A(z-z_*)^r$ 为 $r$ 重根局部模型，则一般扰动的分裂尺度为
$ abs(delta z) ∼ abs(delta p/A)^(1/r). $
这是问题本身的病态性，不能对所有根要求机器精度级前向误差。

但定位重根事件未必同样病态。
当前 #source("src/pygbz2d/continuation/multiple_roots.py")[迭代 MR 求解器] 解增广系统
$
  F(beta_1,beta_2)=mat(f; f_(beta_2))=0,
  quad J=mat(f_(beta_1),f_(beta_2); f_(beta_1 beta_2),f_(beta_2 beta_2)).
$
在普通二重分支点，若 $f_(beta_1)!=0$ 且 $f_(beta_2 beta_2)!=0$，则该 Jacobian 可非奇异。
固定 $beta_1$ 求分裂根与联合定位分支点是不同条件问题。
局部定位误差应结合 $norm(J^(-1))$、最终方程残差及求值误差估计；高重数或更强退化需要另外分析。

实现调用默认 #raw("scipy.optimize.root(method='hybr')") 并提供解析 Jacobian，检查成功标志，但未单独提供尺度归一化的残差保证。
#raw("hybr") 的迭代停止条件不等同于真实位置误差界 [S3, S4]。

完整 MR 流程还包含候选触发、固定 $mu_1$ 后重新求根、弦距离聚类，以及吸附到几何均值。
#raw("cluster_stds") 记录原始聚类散布，可作为诊断指标，但不是严格误差上界。
求得的复 $beta_1$ 分支点与目标固定半径是否相容，也属于需要检查的条件。

== 判定余量与可能的离散变化

#smalltable(
  columns: (0.85fr, 1.65fr, 1.5fr),
  table.header([数值判定], [应比较的余量], [余量不足的后果]),
  [轨迹匹配], [预测与求根误差相对于根间距], [换轨、闭合置换改变],
  [圆内外计数], [根到计数圆的距离], [整数绕数改变],
  [交点与排序], [事件间距、交点斜率、网格分辨率], [漏交点、不同事件合并],
  [SGBZ 电荷], [两侧常规行的模排序间隙], [电荷或边界根对改变],
  [MR 检测], [最小根间距相对于聚类阈值], [近重根误分类、MR 漏检],
  [continuum], [模差或起伏相对于检测阈值], [点集与线集互相误判],
  [线段拼接], [轨迹身份、MR 成员、接缝置换], [连通关系改变],
  [谱内外分类], [小绕数区间与零平台的完整发现], [错误空集或错误非空集],
)

这张表给出待验证的风险机制，不意味着上述错误已经在当前版本复现。

== 轨迹匹配与步长接受

若在弦距离下，每条预测根到其对应数值根的误差界为 $e$，并且数值根之间最小距离大于 $2e$，则正确对应具有直接的最近邻安全余量。
这个充分条件不是当前步长接受判据本身。
当前归一化尺度为
$ s="atol"+"rtol" dot op("median")_(j:"finite") abs(beta_(2,j)), $
并未直接除以局部根间距。
因此接受一步不构成对任意近简并问题均不换轨的证明；现有 polynomial G 测试支持的是对应模型和参数范围。

== SGBZ 环路的拓扑保持条件

令准确环路上的函数为 $f(t)$，数值路径、系数和求值共同产生扰动 $delta f(t)$。
如果
$ norm(delta f)_infinity < min_t abs(f(t)), $
则直线同伦不会穿过零，绕数保持不变。
当前 seed 搜索最大化网格采样得到的最小 $abs(f)$，改善了这个余量，但采样最小值不是整个环路的严格下界。

数值 $mu_(2,"mid")$ 还应保持与理论路径相同的拓扑关系。
其 Hermite 插值误差与 #raw("LOGABS_CLAMP=14") 触发后的路径变化需要单独检查。
在未触发 clamp 的光滑区间可研究通常插值误差；触发时不能将其统称为浮点舍入误差。

== continuum 与有限分辨率

amoeba 用每段轨迹 $ln abs(beta_2)$ 的标准差小于容差检测平坦轨迹；SGBZ 用超过 90% 的网格行满足等模容差建立连通聚类。
有限精度不能无条件区分严格重合与起伏低于分辨率的曲线。
此外，两个统计量依赖自适应网格行权重：改变采样密度可能改变检测结果，即使几何曲线不变。

需要分别研究：保持简并结构的扰动、打破简并的扰动，以及算法无法区分两者的参数尺度。
固定 continuum 与 MR 阈值，仅将二分容差趋于零，不构成完整算法趋于精确问题的收敛过程。
双侧探针距离是识别尺度，也不能直接解释成返回点的严格位置误差界。

== 事件完整性与零平台

当前算法有极值插入、Hermite 多交点预扫描、事件之间的常规分隔行、未知电荷的独立积分，以及周期接缝置换等保护。
它们改善实际稳定性，但还不是全局事件完整性证明。
amoeba 的极值筛选不覆盖同号端点导数间任意隐藏的多个极值；SGBZ 插值预扫描不能发现插值模型没有表现出来的特征。

零平台判定会把候选非空结果改判为谱外，必须同时检查绕数与零点集合完整性。
预筛选的点聚类、小非零绕数面积，以及邻近探针，均应与解析谱边界对照。
失败探针、未决探针与成功空集应分别统计，不能把全部异常合并为谱外结果。

= 复杂度与收敛性

== 计数变量与基本内核

复杂度应依赖实际网格与事件数。固定多项式次数后，逼近退化仍可能显著增加计算量，因而仅给出 $O(K^3)$ 不足以描述完整算法。

#smalltable(
  columns: (0.8fr, 3fr),
  table.header([记号], [含义]),
  [$T$], [Laurent 多项式项数],
  [$K_j=M_j+N_j$], [沿方向 $j$ 的根数；下文 $K=K_2$],
  [$n$], [实际保留的总网格行数],
  [$A$], [continuation 尝试步数，包含拒绝步],
  [$C$], [待处理交点或事件数；需区分原始候选数与最终边界数],
  [$U$], [额外插入网格行数],
  [$B_1,B_2$], [外层、内层绕数求值次数],
  [$Q$], [SGBZ 积分的函数求值次数],
  [$L$], [单个标量精化的函数求值次数],
)

以下采用浮点运算成本模型，忽略整数指数的位复杂度；稠密求根与指派匹配使用通常的三次成本估计。
这不是对底层线性代数每一种实现路径的紧确界。

#smalltable(
  columns: (2.5fr, 1.4fr),
  table.header([子过程], [成本估计]),
  [一次 $f$ 或固定阶偏导求值], [$O(T)$],
  [当前 NumPy 后端部分代入与合并同次幂], [$O(T log T)$],
  [一次 $K$ 次多项式求根], [$O(K^3)$],
  [一组根的解析切向量], [$O(T K)$],
  [弦距离矩阵、全部根对检查], [$O(K^2)$],
  [一次根指派匹配], [$O(K^3)$],
  [每行模排序], [$O(K log K)$],
  [SGBZ 全网格根对模差扫描], [$O(n K^2)$],
  [amoeba 一次水平线交点扫描], [$O(n K)$],
  [对 $C$ 个事件角度排序], [$O(C log C)$],
)

定义一次求根、匹配与切向量计算的成本
$ L_K=O(T log T+K^3+T K). $
普通 continuation 主体约为 $O(A L_K)$，还需加 MR 精化、MR 间密集桥接以及整理成本。
单个交点 Brent 精化若调用 $L$ 次函数，则成本约为 $O(L L_K)$。
MR 增广系统只有四个实未知量，但每次方程与 Jacobian 求值依赖 $T$；其迭代次数必须实测，不能由未知量个数视作固定常数。

== 两种算法的组合成本

SGBZ 每个 $mu_1$ 求值的主要组成是
$
  T_("SGBZ oracle") approx T_("ZM")+O(n K^2)+T_("crossing")
  +T_("seed")+O(Q(T+log n))+T_("bookkeeping").
$
这里的 $log n$ 来自数值路径区间定位。
若 seed 候选数为 $N_("seed")$，逐候选扫描所有节点约需 $O(N_("seed") n T)$。
最终总成本再乘实际 $mu_1$ 求值次数，并计入 continuum 与 plateau 的额外探针。
预扫描的多轮细化及事件管理成本不能遗漏在一次扫描的 $O(n K^2)$ 中。

amoeba 的结构性成本可写为
$
  T_("amoeba") approx B_1(T_("ZM")+T_("extrema")+B_2 T_(w_2)+T_(w_1))+T_("probes").
$
该式以各层典型成本表示；不同参数点成本差异很大时应改写为逐次求和。
每个 $mu_1$ 的 ZeroManager 在内层 $mu_2$ 搜索中复用，所以 $T_("ZM")$ 不再乘 $B_2$。
对由 $C$ 个角度分隔的区间进行中点根计数，沿方向 $j$ 的部分主要成本为
$ O(C(T log T+K_j^3)+C log C). $
因此 amoeba 的 $w_1$ 求值必须保留 $K_1$，不能只使用 minor-axis 次数 $K_2$。
精化交点、缓存刷新和重新扫描另计。

== 网格插入、重新扫描与存储

#source("src/pygbz2d/continuation/zero_manager.py")[#raw("insert_solution")] 使用数组插入，并重新计算受影响段全部行的模排序。
amoeba 还刷新该段 logabs 缓存。
若 $U$ 次插入集中在初始大小为 $n_0$ 的一段，其相关成本可达到
$ O(U(n_0+U)K log K). $
精化交点后重新扫描所有轨迹也会增加成本；若一次水平线求解为 $C$ 个接受事件反复扫描，保守扫描项可按 $O((C+1) n K)$ 记录，实际应使用每轮网格大小。
SGBZ 事件分组、常规行定位和线段拼接同样有数据管理开销；上述结构公式保留了该项，并未将其假定为线性或常数。

主要轨迹存储为 $O(n K)$，根对矩阵和稠密求根工作空间为 $O(K^2)$，再加事件、输出和局部缓存。
线段反复拼接与数组复制可能提高峰值内存，应在成本实验中测量。

== 二分、交点精化与 continuation 的收敛条件

准确连续标量函数与有效异号括号下，普通二分区间宽度为 $D_0 2^(-B)$。
若目标位置误差为 $tau_mu$，通常需要 $O(log(D_0/tau_mu))$ 次区间缩减。
若停止条件是绕数残差，则还需要函数的局部连续性尺度以及计算误差界，不能仅由区间缩小推出停止成功。
例如局部有 $abs(W(mu)-W(mu^*))<=L_W abs(mu-mu^*)$，且误差界 $epsilon_W<tau_W$，才可由足够小的区间推出满足残差目标；位置误差则需要前述下斜率条件。

Brent 在连续函数与有效括号下具有稳健收敛性质，光滑简单零点通常快于纯二分；有噪声或换轨时前提需要重新检查。
它只定位已发现事件，不解决全局事件完备性。

本程序使用一次切向预测加重新求根，步长控制只是沿用 RK 风格，并不是五阶 RK 根轨迹解法。
光滑区间内预测差通常为 $O(Delta theta^2)$，因此在尺度与条件数固定、误差控制主导且尚未触及步长上下限时，可预期步数约随有效容差的 $-1/2$ 次幂增长。
这属于待验证的局部渐近预期，不能推广为 MR 附近的统一收敛阶。

MR 间距平方的导数在普通平方根型分支点附近还可能失去通常光滑性。
不能将普通横截交点的 Brent 分析直接用于全部 MR 触发分支；应分别分析区间触发、点触发和增广系统回退。

== continuum 外层判据的待证条件

光滑情形下，沿准确内层极小值计算外层导数有通常的包络关系。
不可微情形则需要次梯度论证。
当前 amoeba 的 #raw("_handle_continuum") 固定 $mu_(2,c)$，在 $mu_1 plus.minus epsilon$ 两侧求主轴绕数。
一般不可微凸函数的两个坐标方向分别夹住零，并不足以证明联合极小。

一个纯数学例子是
$ R(x,y)=max(x-2y,y-2x). $
原点沿任一坐标轴均为极小值，但沿 $x=y=t>0$ 有 $R(t,t)=-t$，仍可下降。
其完整次梯度不包含原点，尽管两个坐标投影分别包含零。

#note("待证明，不是已复现错误。", [
  上例说明一般凸性不足以支撑该推理，并未证明当前 Ronkin 输入模型会出现同样的错误结果。
  需要证明目标 continuum 情形的额外结构，或直接论证 $0 in partial R(mu_1,mu_2)$。
  在此之前，不应将坐标方向探针称作普适的联合极小值认证。
])

= 已记录的数值结果及证据范围

== HN 模型与误差定义

以下数据抄录自 #source("paper/slides/pygbz2d-details.typ")[算法说明的 benchmark 段落]，相应程序为
#source("application/benchmark-2D-Hatano-Nelson.py")[HN benchmark 脚本]。
说明文档记录计算日期为 2026-09-22，环境为 NumPy 1.26.4、SciPy 1.11.4，使用 NumPy backend 和默认设置。
本次整理没有重新运行这些计算，日期与环境是原记录的属性。

模型为
$ h(beta_x,beta_y)=J_(x 1)/beta_x+J_(x 2)beta_x+J_(y 1)/beta_y+J_(y 2)beta_y. $
两组固定参数分别为：

#smalltable(
  columns: (0.8fr, 2.8fr, 0.6fr),
  table.header([案例], [$(J_(x 1),J_(x 2),J_(y 1),J_(y 2))$], [$E$]),
  [点集], [$(1+i,1.5+1.2i,-1+i,-1.2-0.5i)$], [$1+i$],
  [线集], [$(1,1.5,-1,-1.2)$], [$1$],
)

对每个返回点及每条线的全部存储样本，原记录使用
$
  epsilon_mu=max_(j=1,2) max_n abs(ln abs(beta_(j,n))-mu_j^("exact")(theta_(1,n))),
$
$
  epsilon_E=max_n frac(abs(E-h(beta_(x,n),beta_(y,n))),
  abs(E)+abs(J_(x 1)/beta_(x,n))+abs(J_(x 2)beta_(x,n))+abs(J_(y 1)/beta_(y,n))+abs(J_(y 2)beta_(y,n))).
$
在 $[11]$ 几何下按 benchmark 的坐标映射与角度依赖解析半径比较。
其中 $epsilon_E$ 是归一化方程残差，不是谱能量位置误差。
原验收阈值为 $epsilon_mu<=10^(-5)$、$epsilon_E<=10^(-7)$。

== 八组固定案例

#smalltable(
  columns: (0.6fr, 1fr, 0.6fr, 1fr, 1fr),
  table.header([案例], [GBZ], [样本数], [$epsilon_mu$], [$epsilon_E$]),
  [点集], [Amoeba], [4], [$9.52 times 10^(-9)$], [$2.01 times 10^(-16)$],
  [点集], [$x$ strip], [4], [$2.00 times 10^(-9)$], [$1.00 times 10^(-16)$],
  [点集], [$y$ strip], [4], [$1.44 times 10^(-9)$], [$1.00 times 10^(-16)$],
  [点集], [$[11]$ strip], [4], [$1.44 times 10^(-9)$], [$2.10 times 10^(-16)$],
  [线集], [Amoeba], [354], [$5.97 times 10^(-7)$], [$1.32 times 10^(-10)$],
  [线集], [$x$ strip], [350], [$4.22 times 10^(-7)$], [$1.32 times 10^(-10)$],
  [线集], [$y$ strip], [346], [$5.80 times 10^(-7)$], [$1.28 times 10^(-10)$],
  [线集], [$[11]$ strip], [380], [$3.98 times 10^(-7)$], [$5.51 times 10^(-11)$],
)

原记录中八例均通过验收；每个点集案例有四个点，每个线集案例有两条线。
样本数是输出点或线的存储采样总数，不是 continuation 尝试步数，也不是精度等级。
数据支持“方程残差与 GBZ 约束误差需要分别报告”，但不能据此确定误差来源、收敛阶或其他模型的统一精度。
线集的误差高于点集，也不能未经实验就全部归因于某一个 continuum 或 MR 阈值。

== 已记录的能量扫描

说明文档另记录四种方法使用相同的 $100 times 100$ 能量网格，失败数均为零。
Amoeba、$x$ strip、$y$ strip 共有 4,832 个谱内网格点，$[11]$ strip 为 1,108 个。
这些是有限网格上的结果；零失败率不等同于零误分类率，也不认证网格之间的谱边界位置。

== 已有回归与诊断材料

下列源码为后续实验提供模型和检查方式，本轮只阅读，未重新执行：

+ #source("tests/test_continuation.py")[continuation 测试]：解析根、多重根及 polynomial G 的近简并步长检查。
+ #source("tests/test_amoeba_refinement.py")[amoeba 精化测试]：解析 affine tracks、隐藏双交点、极值插入与持久化求根行。
+ #source("tests/test_counterexamples.py")[反例回归测试]：Haldane、plateau、MR 以及精化回退相关案例。
+ #source("debug_tool/gbz_debug.py")[固定参数绕数调试工具]：区间环路绕数与电荷一致性检查。

历史 #source("diagnostics/profile_sgbz_report.md")[性能报告] 和
#source("diagnostics/imaginary-degeneracy-splitting/report.md")[网格敏感性报告]
可用于设计实验，但涉及历史实现、环境或采样策略，不作为当前版本重新实测的性能与错误结论。

= 数值验证计划（尚未实施）

== 分层实验与检验目标

#smalltable(
  columns: (0.8fr, 1.8fr, 1.8fr),
  table.header([实验], [控制变量与参考], [要检验的结论]),
  [A · 求根], [固定参数；解析或独立高精度根；扫描间距与系数尺度], [条件数放大与重根 $1/r$ 次幂规律],
  [B · 交点], [解析轨迹；扫描水平线、斜率与角度容差], [交点误差公式及相切过渡],
  [C · 绕数], [固定 $mu$；分别用参考事件和独立整数绕数替换数值部件], [角度权重公式与整数误判分离],
  [D · 二分], [准确参考绕数与完整数值绕数分别驱动搜索], [斜率和 Hessian 对 $mu$ 误差的放大],
  [E · 拓扑], [近 MR、相切、continuum、谱边界、接缝参数族], [结构稳定范围与分类转变尺度],
  [F · 成本], [改变次数、项数、容差及事件间距], [内核成本、复用收益与数据维护开销],
)

实验 A 区分“对给定双精度系数的高精度参考”与“从原模型高精度构造的参考”。
当前 CharPoly 会转换为普通复数数组，不能只更换入口数组 dtype 就宣称得到高精度版本。
参考求值与求根需要独立路径，并检查参考精度再次提高后结果不变到目标位数。

实验 B 可复用 $beta_2=a+b exp(i theta_1)$ 的解析轨迹，直接求水平线交点；也应加入 SGBZ 的成对等模交点。
先固定完整拓扑验证普通误差，再靠近相切测量事件数变化。

实验 C 对 amoeba 可比较根计数与独立环路积分；对 SGBZ 可使用更高精度或独立自适应相位追踪作交叉检查。
独立实现一致只能增加可信度，不能自动构成严格认证；必须记录环路离零余量与采样收敛。

实验 D 先隔离二分误差，再接入数值绕数；用多个差分尺度或高精度参考估计局部斜率及曲率。
在平台、非光滑点和近奇异 Hessian 处，不拟合普通局部公式。

实验 E 应优先检查 continuum 外层判据的适用条件，分别考虑保持简并与打破简并的扰动。
现有参数族若不足以区分坐标极小与联合极小，需要另构造具有明确 Ronkin 结构的诊断模型，并先验证模型适用假设。

== 统一记录的指标

+ 连续误差：归一化残差、根与对数根误差、交点两角误差、平均绕数误差、$mu$ 误差、点集匹配距离与线集几何误差。
+ 离散结构：点数、线数、MR 聚类、电荷与事件顺序、闭合置换、连通关系、谱内外分类。
+ 成本计数：多项式求根、尝试与拒绝步、网格行、插入、Brent 求值、积分求值、二分与探针次数，以及各阶段耗时和峰值内存。

失败需拆分为：显式失败、未决、错误谱外、错误谱内、错误连接，以及结构正确但坐标偏差。
不同拓扑之间不应强行逐点配对后给出看似很小的平均坐标误差。
线集除了已存储节点，还应在节点之间增加独立参考采样；离散样本的 Hausdorff 距离也要检查采样收敛。

== 容差扫描与等价变换检查

一次只改变一组参数：步长控制、交点精化、绕数停止、continuum 检测、MR 聚类及探针尺度分别扫描。
其余误差源先固定在足够低且经过检查的水平，再拟合收敛规律。
若发现结构变化，应先定位触发的判定，再决定如何扩展实验，不能直接调整阈值来恢复预期输出。

可加入保持理论问题等价的检查：$f ↦ c f$（$c!=0$）、角度原点平移；amoeba 可检查交换两轴后的对应性。
改变次数或系数尺度的成本实验应同时监测条件数和事件密度，避免将病态性变化误认为单纯次数效应。
SGBZ 与 amoeba 仅在理论上应一致的模型中互作参照，不能要求任意模型两者一致。

== 建议实施顺序与产物

+ 建立解析根、解析交点与独立高精度参考，完成 A、B 的普通区间实验。
+ 固定 $mu$ 分解绕数误差，验证事件角度加权公式，完成 C。
+ 用 HN 等解析模型验证二分误差传播，完成 D。
+ 核实 continuum 联合极小值判据，随后展开 E 的近退化参数扫描。
+ 在已知正确与可解释的案例上完成 F，拟合成本与局部收敛趋势。

建议产物为可复现的参数清单、逐例结构与误差记录、内部计数、独立参考说明，以及误差—容差、误差—间距和耗时—规模的静态图。
实验脚本可先采用独立诊断入口，是否需要增加求解器诊断接口，应根据首轮实验另行提出具体修改计划。
本文只规划这些工作，不将其写作已完成的验证结论。

= 来源索引与后续维护

== 仓库源码入口

文档链接相对于本文件所在目录指向仓库，方便从源码追溯；函数名比行号更适合后续版本维护。

#smalltable(
  columns: (1.2fr, 2.8fr),
  table.header([分析主题], [主要来源]),
  [总体算法], [#source("paper/slides/pygbz2d-details.typ")[pygbz2d-details.typ]；#source("doc/SGBZ.md")[SGBZ.md]；#source("doc/amoeba.md")[amoeba.md]],
  [求根与后端], [#source("src/pygbz2d/core.py")[core.py]；#source("src/pygbz2d/backend.py")[backend.py]],
  [追踪与插值], [#source("src/pygbz2d/continuation/arclength.py")[arclength.py]；#source("src/pygbz2d/continuation/zero_manager.py")[zero_manager.py]；#source("doc/continuation.md")[continuation.md]],
  [重根], [#source("src/pygbz2d/continuation/multiple_roots.py")[multiple_roots.py]],
  [SGBZ 事件与路径], [#source("src/pygbz2d/sgbz/pairwise.py")[pairwise.py]；#source("src/pygbz2d/sgbz/mu2mid.py")[mu2mid.py]],
  [SGBZ 绕数与搜索], [#source("src/pygbz2d/sgbz/winding.py")[winding.py]；#source("src/pygbz2d/sgbz/sgbz_solver.py")[sgbz_solver.py]；#source("src/pygbz2d/sgbz/plateau.py")[plateau.py]],
  [amoeba 事件与搜索], [#source("src/pygbz2d/amoeba/ronkin_winding.py")[ronkin_winding.py]；#source("src/pygbz2d/amoeba/zm_extract.py")[zm_extract.py]；#source("src/pygbz2d/amoeba/bisect.py")[bisect.py]；#source("src/pygbz2d/amoeba/amoeba.py")[amoeba.py]],
  [数值参数], [#source("doc/constants.md")[constants.md]],
)

== 外部算法接口说明

以下链接用于核实接口与算法前提。网站当前版本可能随时间更新，具体重现实验仍须保存实际 NumPy、SciPy 版本。

- N1：#link("https://numpy.org/doc/stable/reference/generated/numpy.roots.html")[NumPy · numpy.roots]。伴随矩阵特征值求根。
- S1：#link("https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html")[SciPy · brentq]。连续性、异号括号与角度停止容差。
- S2：#link("https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html")[SciPy · quad]。积分结果、误差估计及求值诊断。
- S3：#link("https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html")[SciPy · root]。默认 hybrid 方法与解析 Jacobian 接口。
- S4：#link("https://docs.scipy.org/doc/scipy/reference/optimize.root-hybr.html")[SciPy · root(method='hybr')]。迭代停止与选项说明。

== 更新约定

新增实验时应记录代码提交、参数、环境、参考解来源及失败分类。
文档中的三类内容应继续明确区分：源码事实、条件性推导、数值证据。
将待证条件写成已验证结论之前，需要补入相应证明或实验及其适用范围。
