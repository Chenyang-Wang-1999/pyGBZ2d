# SGBZ 性能分析与函数拆分

日期：2026-07-18 — 2026-07-19

## 背景

`collect_GBZ_subsets`（SGBZ 入口）单次调用 point case ~16s（N_points=301），bisection 31 次迭代累计 ~3258 次 `_analyze_boundary_matching` 调用。需要做性能分析并寻找优化点。

## 性能分析结论

### Profile 数据（N_points=101）

| 排名 | 函数 | cumtime | 占比 | ncalls |
|------|------|---------|------|--------|
| 1 | `_analyze_boundary_matching` | 7.36s | 56.8% | 3,411 |
| 2 | `chordal_cost_matrix` | 5.20s | 40.1% | 13,644 |
| 3 | `to_sphere_r3` | 4.02s | 31.0% | 27,288 |
| 4 | `_solve_sorted_roots_at` | 3.86s | 29.8% | 10,301 |

### 根因定位

sweep mode（`refine_continuum=False`）在 bisection 的 31 次迭代中，即使没有连续统（continuum），Stage 3 仍对全部 ~99 对相邻网格点做完整的 `_analyze_boundary_matching`。**2972 次调用全部返回 `is_confident and not cross_matches` → 跳过，纯浪费。**

关键诊断数据：
```
mu1=-1.000000 (bracket 扩张, 远在 GBZ 外): 101 calls (stage3=101, recurse=0)
mu1=+0.000000 (bisection 初期):            99 calls  (stage3=99,  recurse=0)
mu1=+0.200012 (bisection 近收敛):         105 calls (stage3=99,  recurse=6)
...
31 次合计: 3258 calls (stage3=2972, recurse=286)
```

## 架构决策：拆分 `get_roots_and_PMGBZ`

### 问题

原函数用 `refine_continuum` 参数在一个函数里混合了两种本质不同的算法：

- **LineSubset（连续统）检测** — 扫描 degenerate mask 找连续 True 区间，不需要 Hungarian matching
- **PointSubset（偶然点）检测** — 对每对相邻网格点做 Hungarian matching + double-root + exchange margin，O(k³)

混合导致 sweep 模式即使无 continuum 也跑完整 Stage 3。

### 方案

拆分为两个独立函数：

- `get_roots_and_PMGBZ_full` — 完整流水线（连续统边界精修 + Stage 3）
- `get_roots_and_PMGBZ_sweep` — 轻量扫描（先检 continuum，若无则 Hungarian-only 预筛）

### Sweep 函数的 PointSubset 路径

Stage 3 用 `_cross_boundary_matches`（纯 Hungarian + 跨边界检查）替代 `_analyze_boundary_matching`，省去 double-root detection、cluster bounds、exchange margin 分析。对无交叉的网格对，speedup ~2.3x（226ms → 96ms at mu1=-1）。

### Continuum 数据流修正

发现 sweep 函数的 continuum 返回值（假的 LineSubset）会被传到 `collect_GBZ_subsets` → 用 `is_gbz` 判断是否触发精确求解。这是错误的数据流。

修正：`handle_continuum` 不再接收 sweep 的 GBZResult，直接调用 `get_roots_and_PMGBZ_full` 计算精确结果。删除 `solve_SGBZ_for_E` 的 `refine` 参数。

修正后的 continuum 路径：
```
winding_at → sweep 检测 continuum → winding=None
  → handle_continuum(mu1) → resolve limits
    → straddle zero → get_roots_and_PMGBZ_full → 精确 GBZResult
    → 直接返回给 collect_GBZ_subsets
```
sweep 的 placeholder 只在内部传递，不会被用户看到。

## 文件改动

| 文件 | 改动 |
|------|------|
| `pmgbz_detector.py` | 拆出 `get_roots_and_PMGBZ_full` + `get_roots_and_PMGBZ_sweep`，sweep 用 Hungarian-only 预筛 |
| `strip_winding_number.py` | `get_strip_winding` 根据 `refine_continuum` 分发到 `_full` 或 `_sweep` |
| `SGBZ.py` | `handle_continuum` 改为自己算 GBZResult；删除 `refine` 参数；简化 `collect_GBZ_subsets` |
| `__init__.py` | 导出 `get_roots_and_PMGBZ_full` + `get_roots_and_PMGBZ_sweep` |
| `gbz_types.py` | 移除 `_SWEEP_CONTINUUM_PLACEHOLDER` sentinel（不再需要） |
| `tests/test_sgbz.py` | 更新 export 检查 + 直接调用处 |

## 待优化

- `TODO/SGBZ-optimization.md` 记录了后续优化方向（P0: Stage 3 轻量预筛、缓存 to_sphere_r3、跳过 gap 大时的 double-root 检测）
- `get_roots_and_PMGBZ_full` 仍混合 LineSubset + PointSubset，可进一步拆分
- Sweep 函数的 PointSubset 路径可以加更激进的预筛（模长 gap + 欧氏距离）

## 测试

63/63 通过（SGBZ 21 + Amoeba 13 + gbz_types 29）
