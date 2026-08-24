# SGBZ pairwise：MR 单通道化 + seam 合并列帧重标

日期：2026-08-18

## 背景

y-SGBZ Haldane（gain-loss，supercell [[1,1],[-1,1]]）在 `E_ref = 1.212`
处 `solve_SGBZ_for_E` 抛出
`bracket bisection failed to converge after 60 iterations`
（μ₁=0.135328598265，W=−0.1028034060294209 非整数 → W(μ₁) 在该处
拓扑跳变，二分永不收敛）。2026-08-16 日志遗留问题，本轮定位并修复。

## 根因（两处，均在事件合并，不在检测）

用户人工审查定位方向；实证补充：

1. **MR 双通道漏检**：`collect_pair_events` 把落在段 MR 边界行的 touch
   事件 skip 掉（"由 MR 记录通道负责"），但该通道只物化跨骑 M-1/M 的
   cluster——其他模长重合被两个通道同时丢弃。
2. **seam 合并列帧错标（E=1.212 直接根因）**：整圆段的事件在 2π 侧
   检出时携带段 track 帧列号；θ=0 帧与 2π 端 track 帧相差单项式置换
   `boundary_perm`（本例 [3,6,5,7,0,2,1,4]，非平凡 monodromy）。
   `group_events` 的 seam merge 把 2π 侧事件直接拼进 θ=0 锚定的
   group 而不翻译列号 → `finalize_event_groups` 按 frame-0 解释
   track 帧标号 → `{M-1,M}` 覆盖测试失败 → `point_columns=()`
   静默丢弃。E=1.212 数据：等模对 {0,7}@2π → frame-0 应为 {3,4}；
   三对 seam 等模对 {5,6}/{1,2}/{3,4} 逐一验证。

检测/精化层经逐 pair 原始扫描确认无漏检（26 变号 = 26 事件）；
brentq 三层防护确认不会把 2π 侧根拉回 0。

## 修改

- `brute_force_SGBZ/pairwise.py`
  - `collect_pair_events`：删除 MR 边界行 skip；touch 统一收集，
    落在 MR 边界行的标 `is_mr=True`（荷 → None，hard）。
  - `PairEvent.is_mr` / `EventGroup.is_mr` 字段（`_fill_group_components`
    中 OR 汇聚）。
  - seam merge：`last.theta > π` 时对 `last.events` 调
    `_relabel_event_to_theta0_frame`——`inv_perm` 翻译
    cols/rep，ia/ib 按 frame-0 rep 在 ItemView 重解析；
    `theta_star` 保留检出位置。0 侧事件不动。
    （经审查删去"段覆盖整圆"冗余守卫：merge_tol=1e-10 下
    circ_gap 触发本身已蕴含该条件。）
- `brute_force_SGBZ/winding.py`
  - 删除 `_mr_boundary_entries`、`_MR_PROXIMITY_TOL` echo drop、
    MR 物化块。`detect_crossings_simple` 成为纯
    EventGroup→PointSubset 物化；`kind='mr'` 读 group 的 `is_mr`。
- 全项目 `2 * pi` → `gbz_types.TWO_PI`（70 处，纯等价重构；
  替换前已核实 `2*math.pi == 2*np.pi == 2*cmath.pi == math.tau`）。

## 回归

- `pytest tests/`：189 passed, 2 skipped。
- `--run-slow` Haldane 反例（E=-1.56 → (12,0)；E=-1.406 → (6,0)）：不变，7 passed。
- E=1.212 复现（`diagnostics/repro_y_sgbz_E1212.py`）：不再
  RuntimeError；返回 `index=(0, 6), is_continuum=True`（6 条
  LineSubset 的 continuum 边界，物理合理）。失败 μ₁ 处 seam group
  `point_columns=(3,4)`、`column_q={3:-1, 4:+1}`，soft 电荷总和为 0
  （守恒检查通过），W=0.192（整数 0 的种子区 + 边界层贡献）。
- 新测试（均 A/B 验证：旧逻辑 fail / 新逻辑 pass）：
  - `test_mr_boundary_touch_is_collected_not_skipped`（poly_F 边界 MR）
  - `test_seam_merge_relabels_2pi_side_events_to_theta0_frame`
    （非平凡 boundary_perm=[1,0] 合成 ZM）
  - `test_transient_inf_row_does_not_discard_pair_events`（瞬态 0/∞）

## 残留

- MR 单通道化后，非 cluster 根在 MR 边界行的模长重合现在按普通
  touch 事件处理（soft/possible charge）；`is_mr` group 整体标
  hard（charge None）——若物理上 MR 处存在可定荷的横向穿越，
  需要后续按事件粒度（而非 group 粒度）细分 hard 标注。
- review 报告第八节疑点 1（`_mark_mr_tangents_inf` 帧错位）仍在。
