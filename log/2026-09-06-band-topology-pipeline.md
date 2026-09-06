# GBZ 能带拓扑构建流水线：band 聚类 → 环面三角网格 → 自适应加细

日期: 2026-09-06

## 背景与目标

输入是 `collect_GBZ_subsets` 的 E 扫描数据（复 E 网格，每节点一个
`GBZResult`）。目标：把离散 GBZ 数据点组织成每个 band 的
θ₁–θ₂ 环面三角网格（npz 落盘），供拓扑分析使用。

约定假设：band 之间 well-separated、带隙不闭合；超胞折叠带
（y-SGBZ 案例）不支持，遇到时聚类/网格质量诊断会如实暴露。

整体流水线六步（本次定型）：

1. 在 E 网格上求解（现有 sweep pkl）；
2. 找 index 变化边（含 ↔`(0,0)` 出谱边）；
3. 二分定位边界并**保留全部探测点**（enriched pkl）；
4. 聚类；
5. 每 cluster 建三角网格；
6. 按边长自适应加细，迭代到收敛。

## 一、band 聚类（`pygbz2d.experimental.band_clustering`）

- **选型**：半径图连通分量（DBSCAN 取 min_samples=1）。对比过
  HDBSCAN/谱聚类/ToMATo（需新依赖，项目仅 numpy+scipy）、k-means/GMM
  （凸簇模型不符）、切片匈牙利跟踪（E 为二维网格无全序）后选定。
- **实现**：`flatten_results`（线采样**不抽稀**）→ R⁸ 嵌入
  （角度 (cos,sin)² 解决 θ=0≡2π 缝；E 块按网格步降缩放
  `ALPHA_E=0.25`；μ 块 1/max(span,1)）→ cKDTree + 连通分量。
  `cluster_bands` 主入口；eps 自动取"簇数-eps 曲线最宽平台"的几何中点，
  扫描上限 `EPS_SCAN_MAX=3.0`，全并拢早停防密集云的点对爆炸。
- **两个几何事实**（设计依据）：固定 β 处不同 band 的能量不同，片与片
  只在含 E 的全空间分离，故 E 必须进度量；μ 分离能力经合成云
  （仅 μ 不同的两 sheet + `w_mu=0` 对照）验证。
- **两个教训**：按 θ₁ 宽度抽稀线采样会摧毁近竖直弧（Hermitian 模型
  等能线的常态形状）→ 已移除，性能由扫描上限+早停承担；enriched 数据
  的 E 步长必须取自**原始网格轴**（探测点位于网格点之间，会稀释中位
  步长、放大 E 块、碎片化聚类）→ `cluster_bands` 增加 `d_re`/`d_im` 参数。
- **验证**：gain-loss amoeba 2 簇（margin 11.5）；Hermitian Haldane
  （E 投影重叠 +4.46、直接带隙 1.52 不闭合）2 簇（margin 0.44，
  分离完全由 θ 块完成）；合成云 μ 分离、缝穿越、早停等 11 个单测。

## 二、环面三角网格（playground `demo_torus_mesh.py`，Plan A）

- 3×3 复制周期 Delaunay + 中心 tile 筛选 + mod 映射 + 去重；
  拓扑 QC：χ、边界环、连通分量（χ=0 无边界 = 闭合 genus-1）。
- **数值卫生**：1e-5 rad 容差去重（拼接弧共享端点只收敛到 ~1e-6，
  Qhull 会把这类近重合点静默扔出三角化 → 孤立顶点、χ 虚高）。
- **seam keep-rule 超填（已修复，2026-09-06 晚）**：根因是复制平面的
  三角化不保证周期性——同一 torus 边在 θ≈0 和 θ≈2π 各有一个平面实现，
  浮点平移副本相差 1 ULP，Qhull 在近退化的镜像对称构型上两侧判定不一
  致，keep 规则取并集 → 边重数 3/4（y-SGBZ 触发：χ=19/41）。修复
  （按"等价性只看 index"的原则）：等价/规范判定全程用 `j // n`（tile）
  与 `j % n`（原始点）；坐标在喂给 Qhull 前定标到周期 1 + 2⁻³² 网格
  量化（副本 = 量化值 + 整数偏移 = 比特级平移）；叠加 index 键控微扰
  （仅依赖 `j % n`，≤2⁻³² 周期）以平移协变方式打破退化平局。修复后
  y-SGBZ 两 cluster 均 χ=0、E=3F/2 精确、定向一致；其余数据集回归
  不变。
- 结果：Hermitian / gain-loss 全部 cluster 为 χ=0 闭合环面，单分量。

## 三、边界富集（playground `enrich_boundary.py`，步 2–3）

- E 网格 4-邻接边，两端 `index` 不同即为转换边；沿连线二分
  （默认 6 iters，精度 = 网格步/64），谓词为 index 比较（`(0,0)` 与
  `(6,0)`/`(12,0)` 同权）；**全部探测结果并入数据集**。
- checkpoint 每 25 条边，`--resume` 续跑；`--n-procs` 并行。
- gain-loss amoeba：1192 条转换边 × 6 = 7152 探测点，零失败，
  串行 6.8 h（用户机并行可显著缩短）。

## 四、自适应加细（playground `demo_mesh_refine.py` + `demo_pipeline.py`）

- **v1（废弃）**：三角形面积标准 + 三边全插中点 + 交错求解。诊断出
  sliver 盲区：cluster0 的 452 条 >0.25 rad 长边 100% 属于面积达标的
  细长三角形，"0 轮收敛"但 max edge 0.541。
- **v2（现行）**：
  - 选择标准 = torus 边长（绝对阈值 `--edge-thresh`，默认 0.2 rad）；
  - 两级循环：内层纯几何（超限边上 dyadic 中点嵌套 + 重三角化，直到
    所有边 ≤ 阈值；预测器 = 只锚定真实顶点的 Clough-Tocher C1 插值，
    建在 3×3 复制云上保证缝两侧邻居正确；预测值永不反馈进预测器）；
    外层批量求解 → Hungarian 归并（(cos,sin)² 嵌入，容差 0.5）→
                以求解位置合并 → 重三角化 → 复查；
  - 边界 stall 回退：预测出谱或全部子集被拒（index 变化位置，如两个
    (12,0) 连线中点落在 (6,0) 区）时，与两个真实锚点的 E 连线二分
    （谓词 = index 相对锚点是否变化）定位边界点并归并；预测侧 index
    与锚点相同者跳过（纯匹配失败不是边界）。
- **结果（用户并行跑）**：cluster1 max edge 0.197 ≤ 0.2 收敛，
  χ=0 闭合环面保持，定向一致。
- **并行**：`_solve_one` worker（一个待定点 = 1 直接求解 + 至多
  2×iters 探测），串行/并行同一条代码路径。

## 五、基础设施与修复

- **`GBZResult`**：`is_gbz` 改为 `success and index != (0,0)`；
  `is_continuum` 字段删除（审计：amoeba 从不设置；sgbz continuum 走
  `index == (0, n_1d)`，escape 条款实为死代码）。同名内部 dict key、
  `GBZDebugReport.is_continuum` 不受影响。
- `tests/test_constants.py`：`read_text` 加 `encoding="utf-8"`
  （GBK locale 下解码 backend.py 失败的预存 bug）。
- 预存失败 `test_mu01_has_single_seam_event_group`：seam touch 的
  `PairEvent.direction` 两侧单侧导数天然反号（d 在缝处 ~−|θ| 尖点），
  符号零下游消费；测试已删。
- `playground/plot_torus_mesh.py`：pyvista 环面堆叠视图（cluster index
  作 z 偏移；`--mesh` 直接读 npz；跨缝三角形在环面视图零特判）。
- `demo_pipeline.make_figure` 平面展开图的 2π−0 长条过滤
  （`seam_display_filter`，显示专用，网格不动）。
- `debug_tool/mesh_orientation.py`：网格定向检查（带符号面积普查 +
  有向边平衡 + 非流形边），CLI 退出码可挂 CI；现有网格全部 OK。
- `demo_pipeline.py` 每 cluster 落盘
  `mesh_cluster{cid}.npz`（`verts`(n,2)、`triangles`(m,3)、
  `E`/`mu1`/`mu2`）。

## 六、使用指南：从初次扫描数据到最终 npz 网格

前提：已有初次扫描 pkl（如 `playground/Haldane-model-gainloss.py` 的
`sweep_amoeba()` 产物），内含 `E_real` / `E_imag` / `results` /
`coeffs` / `degs`。以下命令均在**仓库根目录**运行；
`--n-procs` 按机器核数调整；sgbz 方法的扫描在两个脚本都加
`--method sgbz`（必须与初扫方法一致）。

```bash
# 步骤 2-3：E 网格边界富集（数小时量级；断点续存可 --resume）
python playground/enrich_boundary.py data/Haldane-gain-loss-amoeba.pkl --n-procs 12
#   → data/Haldane-gain-loss-amoeba-enriched.pkl

# 步骤 4-6：聚类 + 环面网格 + 边长加细收敛（主要耗时在求解，同样可并行）
python playground/demo_pipeline.py data/Haldane-gain-loss-amoeba-enriched.pkl --n-procs 12
#   → playground/band_cluster_out/Haldane-gain-loss-amoeba-enriched/mesh_cluster{cid}.npz

# 可选：定向/拓扑检查（退出码 0 = 全部合格）
python debug_tool/mesh_orientation.py \
    playground/band_cluster_out/Haldane-gain-loss-amoeba-enriched/mesh_cluster*.npz

# 可选：可视化（3D 环面堆叠，多个 npz 依次为多层）
python playground/plot_torus_mesh.py --mesh \
    playground/band_cluster_out/Haldane-gain-loss-amoeba-enriched/mesh_cluster1.npz \
    playground/band_cluster_out/Haldane-gain-loss-amoeba-enriched/mesh_cluster0.npz
```

常用旋钮：`--edge-thresh`（绝对 rad，默认 0.2）、`--max-iters`
（外层求解轮上限，默认 10）、`--match-tol`（归并容差，默认 0.5）、
`--eps`（手动指定聚类半径，默认自动平台）。
