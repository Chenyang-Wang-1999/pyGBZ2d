# 会议记录：GBZ torus 网格构建原型成功 + 失败方案清理

日期: 2026-08-29

## 一、目前成功的原型

文件：`playground/demo_torus_prototype.py`（输出 `playground/gbz_torus_prototype.pkl`）。

### 流水线

1. **bulk**：干净 slice 之间用 FKU（`demo_fku.fku_triangulate`），certificate 门控。
2. **临界区间**：`[sa,sb]` certificate 失败 → E 二分；二分出的非临界 slice 全部重采样后注册。收紧到 flank loop 的 saddle gap（两个 hyperbola 分支的自接近距离，约 2√δ）与局部采样间距 `h` 同量级（`gap <= 2h`）才停。本轮自动停在 `E ≈ ±0.000391`（不是硬编码 0.0001）。
3. **临界 slab**：不再依赖 E=0 采样。saddle 位置由单条 flank loop 的自接近局部极小估计（torus 中点）；每个 saddle 的 4 个 passage 在展开的 (θ₁,θ₂) 图里 fan 到临界顶点；四个互补弧用 FKU 缝合。
4. **cap**：谱边缘复用 `demo_morse_mesh.grow_cap` + `fill_holes` 补小洞。
5. **定向统一**：BFS 沿共享边翻转，最后按 (θ₁,θ₂) 图定向翻转整体符号。
6. **均匀重采样**：所有 clean slice 按 R⁶ 弦长弧长均匀重采样；slab 左右 flank 用 `resample_slice_pair` 重采样到相同点数。
7. **精确求解**：重采样点先用 Newton 解 `f(E, exp(mu1+iθ₁), β₂)=0` 精确化（初值=插值 β₂）；失败回退到 `CharPoly.solve_roots_1d` 全根选最近根。

### 验证（trivial model 全 torus）

```
n_vertices: 19958
n_edges:    59874
n_faces:    39916
euler:      0          <- 闭 torus
boundary:   0
nonmanifold:0
flipped:    6         <- 局部 seam 展开符号，非定向问题
edge_len_median: 0.0499
edge_len_max:   0.2667

exact refinement:
  newton_ok:   17987
  fallback_ok: 23
  failed:      0
  max |f|:     8.88e-16
```

### 仍然存在的已知问题

- slab 互补 FKU 条带几乎全是钝角三角形（最大角 ~178°）。原因是左右 flank 在互补弧区域间距 ~δ（极小），薄带三角剖分固有；当前边长（median ~0.057，max ~0.104）与沿 loop 点距相当，暂时接受。
- fan 三角形也是钝角（fan 几何固有），fan 半径 `2*gap+0.02` 有界。
- cap 仍靠 `fill_holes` 补 3-圈小洞。

---

## 二、失败方案（已删除）

### 1. Gudhi Tangential Complex（`src/pygbz2d/unsupported.py`，`demo_morse_mesh.tc_patch`）

- 想法：R⁶ 点云上跑 Gudhi TangentialComplex（intrinsic dim=2）。
- 失败原因：
  - 局部 PCA 盲估切平面，没有利用已知纤维结构；
  - 云边界处 star 不可靠，需要 padding slice 并丢弃边界三角形，仍然无法保证 fiber 边被尊重；
  - 临界 patch 上 gudhi 对孤立点/小簇会崩。
- 结论：全局 TC 不适合这个问题的精度要求。

### 2. Hultquist advancing front（`playground/demo_advance_front.py`）

- 想法：贪心最短对角线逐三角形推进。
- 失败/弃用原因：局部贪心，整体不是最优；FKU 的全局 DP 更优且同样 O(nm)，bulk 直接采用 FKU。

### 3. nexus + FKU + MSS 分支桥接（`playground/demo_build_mesh.py`）

- 想法：相邻 E 的 nexus 图，同拓扑 FKU，拓扑变化用 bisection + Meyers-Skinner-Sloan 分支。
- 失败/弃用原因：拓扑变化处 MSS 的 1-to-N 分支处理在通用 nexus 图上脆弱；后续被 Morse 思路替代。依赖它的调试脚本一并删除：`demo_bisect_only.py`、`demo_delaunay_bridge.py`、`_study_mr.py`、`_test_mr_split.py`、`plot_gbz_mesh.py`。

### 4. 局部切平面 Delaunay 桥接（`playground/demo_delaunay_bridge.py`）

- 想法：临界 gap 附近 PCA 切平面投影 + 2D constrained Delaunay。
- 失败/弃用原因：切平面由邻域 PCA 估计，约束边恢复不可靠；没有用纤维结构。

### 5. Morse-aware contour stitching（`playground/demo_morse_mesh.py` 中的 `glue_critical` / `grow_patch` / `tc_patch`）

- 想法：临界 fiber 交叉点 fan + FKU 补弧；或 greedy Delaunay growth。
- 失败原因：
  - `glue_critical` 每个 flank 每个 crossing 只找到一个 passage（实际有两个），euler=109、nm=1、flipped 234；
  - `grow_patch`/`grow_region` 的 greedy growth 留洞、错误配对，全 torus euler=24、nm=1、flipped 16195。
- 处置：`tc_patch` 已从文件中删除；`demo_morse_mesh.py` 文件本身保留，仅作为 helper 库（oracle、nexus、MeshBuilder、validate、FKU 封装、grow_cap 等），其中 `glue_critical`/`grow_patch` 仍作为 legacy 保留（viewer strip 模式和 fallback 还引用）。

### 6. NN-bridging v1/v2（`playground/nn_bridging.py`、`playground/nn_bridging_v2.py`、`playground/nn_bridging_v2_view.py`）

- 想法：两遍最近邻 rung（A→B 全部，B→A 补剩余）+ 在包含这些边的前提下建三角形网格。
- 失败/弃用原因：拓扑一致时 rung 单调、效果尚可；拓扑不一致处 rung 交叉，流形条带三角剖分在数学上无法同时包含所有 rung。v2 用 LNDS 保留最大不交叉子集并报告 excluded，但跨拓扑区域仍桥接错误（大量 many-to-one fan、长边），最终被“FKU bulk + 鞍点 fan slab”替代。

---

## 三、本次文件清理

删除：

- `src/pygbz2d/unsupported.py`
- `playground/demo_advance_front.py`
- `playground/demo_build_mesh.py`
- `playground/demo_bisect_only.py`
- `playground/demo_delaunay_bridge.py`
- `playground/_study_mr.py`
- `playground/_test_mr_split.py`
- `playground/plot_gbz_mesh.py`
- `playground/gbz-mesh-demo.py`
- `playground/nn_bridging.py`
- `playground/nn_bridging_v2.py`
- `playground/nn_bridging_v2_view.py`

保留：

- `playground/demo_torus_prototype.py` — 成功原型
- `playground/demo_critical_slab.py` — slab helper（passages / fan / FKU strips）
- `playground/demo_critical_slab_view.py`、`playground/morse_torus_view.py` — 可视化
- `playground/demo_morse_mesh.py` — helper 库（oracle、nexus、MeshBuilder、validate、grow_cap 等）
- `playground/demo_fku.py` — FKU 实现
- `playground/demo_nexus.py`、`trivial_model.py`、`trivial_line_cache.py` — 依赖

未提交的旧数据/图片（`gbz_*.pkl`、`*.png`）保留在 workspace，未纳入 git。
