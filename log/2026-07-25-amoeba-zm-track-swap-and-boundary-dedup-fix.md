# amoeba 根追踪身份交换 + 边界零点 dedup 修复 — 消除 BASE/XY 谱分歧

日期：2026-07-25

## 背景

gain-loss Haldane 模型，同一物理参数 `ALL_PARAMS=(1, 0.5, π/3, 0.5j, 0)` 用两种坐标系算 amoeba 谱：原胞 BASE（2-site，多项式 26 项）与超胞 XY（2-site supercell，基 `[[1,1],[1,-1]]`，64 项）。两者描述同一物理模型，amoeba 谱（E 是否在 GBZ 内）必须一致。

101×101 网格上发现 **87 个 PointSubset 计数不一致**，其中 **5 个是真正的能谱错误**：BASE 判 E 在 GBZ 外（`index=(0,0)`），XY 判 E 在 GBZ 内（`index=(2,0)`）。复现点 **E = 0.365 + 0j**，两模型都收敛到 μ₁=μ₂=0。调查记录见 [TODO/amoeba-coord-mismatch-summary.md](../TODO/amoeba-coord-mismatch-summary.md)。

该 summary 把根因候选定为"疑点 A：BASE 与 XY 的 θ₁ 参数化是否对应同一物理 k 点"。本轮排查推翻此判断，定位到两个独立的 ZeroManager/提取层 bug，修复后 BASE 与 XY 在全部 10201 点上完全一致。

## 推翻 summary 的判断

**疑点 A 作为根因候选站不住——amoeba 谱是胞不变量，不存在"合法不同"。** BASE 与 XY 的特征多项式对应同一代数簇（超胞是原胞的 2 重覆盖重参数化），Ronkin 函数、amoeba、所有 winding/零计数都是该簇的几何对象，与用哪个胞参数化无关。`is_gbz`、w1/w2、孤立零点计数原则上必须胞一致；出现差异只可能是某一侧算错，不可能是"扫的 k 方向不同所以根结构不同属正常"。

**真正的线索是"全场统一 +2"**——87 个不一致点的 PointSubset 计数全部差 2，是系统性过计数签名，不是单点采样精度问题能解释的。summary 把精力放在 E=0.365 单点的"touch-but-don't-cross gap ~0.002"故事上，该故事哪怕成立也只能解释 1 个点。

## Bug 1：根追踪身份交换（integrate_segment）

### E=0.365 真相

密网格（20001 点，θ₁∈[3.0,3.3]）连续性追踪（nearest-neighbor 跨采样匹配）XY 的中间两个 β₂ 根：

- track1：ln|β₂| ∈ [−0.0759, −0.0021]，**始终 < 0**
- track2：ln|β₂| ∈ [+0.0021, +0.0759]，**始终 > 0**
- 两支在 θ₁=π 处最近，gap |track1−track2| ≈ 4.1e-3，**永不相交**
- 复数序列光滑（Im 光滑过零，无 π=−π 跳变），无身份交换

BASE 同样密网格结果结构一致（track1∈[−0.038,−0.001]、track2∈[+0.001,+0.038]，π 处擦肩 gap ~2e-3）。**两侧密网格下都不跨 |β₂|=1。**

但 ZeroManager（ZM，自适应粗网格步长 ~0.016）跑出来的 track1/track2 **各自** ln|β₂| 范围都是 [−2.72, +2.72]——两端跨过 0。定位到 i=119→120 这一步（θ₁ 跨过 π）：

- i=119：track1=−0.9974（<1）、track2=−1.0024（>1）—— 正确
- i=120：track1=−1.0021（>1）、track2=−0.9979（<1）—— **标签对调**

summary 记录的"事实 1：XY 确有 2 个 crossing"是错的——那不是两条支各自跨 1，而是同一步里两支互换标签后伪饰出的假象。

### 根因：裸 chordal 匹配选错锚

`integrate_segment`（continuation/zero_manager.py）每个接受步对 `roots_old → roots_new` 做裸 chordal Hungarian 匹配：

```python
matches = hungarian_match_indices(roots, roots_new)
reordered = roots_new[matches]
```

在最近接触点附近两支几乎共轭对称时，交叉配对的总 chordal 距离比保持身份更短。i=119→120 这步的 cost 矩阵中间 2×2 块：

```
              r2[1](−0.9979)  r2[2](−1.0021)
r1[1](−0.9974)    0.0158         0.0091
r1[2](−1.0024)    0.0091         0.0158
```

- 保持身份（对角）：0.0158+0.0158 = 0.0316
- 交换身份（反对角）：0.0091+0.0091 = 0.0183 ← Hungarian 选了这个

交换光滑且不可逆（翻面后永久不回），整条后半段 track1/track2 标签对调，每条 track 看起来都跨过 |β₂|=1，凭空多 2 个零点 → 谱判定错误。密网格也救不了——只要两支共轭对称，交叉配恒更短，与 gap 大小无关。

**关键对照**：`arclength_step` 内部 `estimate_error`（arclength.py）早已用 `hungarian_match_indices(predicted, actual)` 做过一次正确匹配——tangent 预测点在"正确延续"那一侧，交叉配反而更贵（同步 keep=0.0016 < swap=0.0172）。但那次 `matches` 算完被丢弃，只返回标量 error。`integrate_segment` 反而用错误锚重算一次。每步做了两次 Hungarian：一次对、扔了，一次错、用了。

`roots_new` 来自 `solve_roots_1d`，返回是 `np.roots` 顺序（unsorted，跨 θ₁ 不稳定）——这正是"必须做一次 Hungarian"的理由，而 bug 在于用了错的锚。

### 修复

按"在 arclength_step 层把根排好序"的分层：

- **`arclength_step`**（continuation/arclength.py）：接受步时用预测锚 `hungarian_match_indices(predicted, roots_new)` 把 `roots_new` 就地排成 track 顺序再返回。`StepResult.roots_new` 现在是 track-ordered。
- **`integrate_segment`**（continuation/zero_manager.py）：删掉冗余的第二次裸 chordal 匹配，直接用 `result.roots_new`。
- `estimate_error` 签名不变（测试依赖其返回标量）。

padding 根（0/inf）对齐安全：`predict_roots` 把它们原样保留，`to_sphere_r3` 映射到固定锚点（0→南极、inf→北极），Hungarian 自然配对。

### 验证

**单点 E=0.365**：修复后 XY 与 BASE 的 track1/track2 都单侧不跨 0（track1 全负、track2 全正），假 crossing 消失，与密网格物理事实一致。

**全网格复跑**（101×101 = 10201 点，BASE 1373s + XY 1425s，12 进程并行）：

| | 修复前 | Bug1 修复后 |
|---|---|---|
| PointSubset 计数不一致 | 87 | **2** |
| 真正能谱错误（is_gbz 分歧） | 5 | **0** |

5 个真谱错误全部消除；系统性 ±2 过计数签名消失，印证诊断（87 点绝大多数源于同一身份交换机制）。残余 2 个 PointSubset 计数差（模式 BASE=6 vs XY=5，计数差 1）is_gbz 一致、对谱判定无害——见 Bug 2。

## Bug 2：边界零点 dedup 用错指纹（extract_amoeba_subsets）

### 残余 2 点的真相

Bug 1 修复后剩余 2 个不一致点（i=5107 E=1.289、i=5128 E=2.906），模式都是 BASE=6 vs XY=5，is_gbz 一致。两点的共同特征：在 β₁=1（θ₁=0≡2π 圆边界）处有一对共轭单位根（E=1.289 的 `−0.1248±0.9922`），BASE 保留共轭对（6 点），XY 只输出一个（5 点）。

**单点重跑稳定 6、12-worker mp 网格偶发 5**——典型的浮点抖动签名，但根因不是 mp 噪声。

### 根因：Rule 2 dedup 用 θ₁ 而非零点身份

`extract_amoeba_subsets`（brute_force_amoeba/zm_extract.py）的 Rule 2：

```python
if kind == 'zero':
    key = round(t1, 12)
    if key in seen_zero_theta1:
        continue        # ← 误丢共轭根
    seen_zero_theta1.add(key)
```

两个问题：

1. **误丢共轭根**：θ₁=0≡2π 边界上两个共轭单位根是**两个不同零点**（不同 track j），但 θ₁ 相同 → `round(t1,12)` 同 key → 第二个被 dedup 丢弃。本该 6 点变 5 点。
2. **依赖 fsolve 浮点收敛**：`_finalize_crossing` 调 `_find_exact_crossing`（fsolve）精修 θ₁，返回值未归一化到 [0,2π)。两个共轭根的精修 t1 一个收敛到 0、一个到 2π（同一物理点的两种表示）。`round(0,12)≠round(6.283…,12)` 时 → key 不同 → 都保留（6）；若 fsolve 把两者都收敛到同侧 → key 相同 → dedup 丢一个（5）。哪种发生取决于 `np.roots`（伴随矩阵特征值）的浮点微差，mp 子进程下 BLAS 线程数/累加顺序不同就会变 → "单点稳定 6、网格偶发 5"。

**根本错误**：零点身份本该由 `(端点, track j)` 唯一标识，而 Rule 2 退而用 θ₁ 浮点值做指纹——既丢共轭根，又对浮点收敛敏感。`zm_extract` 完全没用 `boundary_perm`，不知边界处左右 track 的对应关系。

### 修复

dedup 改用零点身份键，新增 `_zero_identity_key(zm, seg, i, j, boundary_perm_inv)`：

- **segment 内部**（`0 < i < N−1`）→ `('interior', id(seg), i, j)`，唯一，永不 dedup。
- **共享 MR 端点**：相邻 segment 在 MR 处共享同一行（`multiple_roots[R].roots`，`ZeroManager.run` 把它 prepend 到下一段 i=0）；右端 `('mr', seg.right_mr, j)` 与下一段左端 `('mr', seg.left_mr, j)` 同键 → 折叠。
- **圆边界 θ₁=0≡2π**：首段左端 `('circle', j)`；末段右端经 `boundary_perm_inv[j]` 折叠到左端 track（`roots_right[boundary_perm]==roots_left` → 右端 track j ↔ 左端 track `boundary_perm[j]`）→ `('circle', boundary_perm_inv[j])` 同键。

这样共轭根（不同 track）键不同 → 都保留；同一零点的左右/MR 双重表示键相同 → dedup。与 fsolve 收敛到 0 还是 2π 完全无关。`boundary_perm` 在无 MR、单 segment 时也有效（`completed` 分支用 `hungarian_match_indices(right_end, left_boundary_roots)` 设置）。

### 验证

- E=1.289、E=2.906 单点 4 次 + **12-worker mp 100 次**全部稳定 6（β₁=1 共轭对都保留）——mp 偶发 5 彻底消除。
- 相关测试 78 项全绿。

## 文件改动

| 文件 | 改动 |
|------|------|
| continuation/arclength.py | `arclength_step` 接受分支增加预测锚排序 `roots_new = roots_new[matches]`，StepResult 返回 track-ordered roots |
| continuation/zero_manager.py | `integrate_segment` 删除冗余 `hungarian_match_indices(roots, roots_new)`，直接用 `result.roots_new` |
| brute_force_amoeba/zm_extract.py | Rule 2 dedup 从 θ₁ 改为零点身份键；新增 `_zero_identity_key`；导入 `SegmentData` |

## 反思（呼应 CLAUDE.md）

- 两个 bug 都属"反常输出改进算法"而非 workaround：Bug 1 是匹配判据选错锚（密网格也救不了，只要共轭对称交叉配恒更短）；Bug 2 是 dedup 用错指纹（既丢共轭根又对浮点敏感），非阈值/采样问题。
- 一行错锚 → 全场 87 点连锁错误。底层判据选错会沿调用链放大。
- `estimate_error` 内部已有正确匹配却被丢弃、外部用错锚重算——分层时的隐藏冗余值得作为代码审查的模式警惕。
- summary 把根因 frame 成"θ₁ 参数化对应关系"违反胞不变量约束；正确方向是用胞不变量排除"合法不同"，再找系统性过计数的来源。
- "单点稳定、mp 偶发"不一定是 mp 噪声——本案是浮点收敛路径分歧经 dedup 放大而成，根因仍在算法层的指纹选择。
