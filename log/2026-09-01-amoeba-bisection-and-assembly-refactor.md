# Amoeba 交点探测、二分与 subset 组装重构

日期: 2026-09-01

## 背景

旧的 amoeba 管线在 μ₂ 二分过程中同时做 continuum 检测与 crossing
detection，且 subset 提取逻辑（`extract_amoeba_subsets` 及配套去重 /
MR 拼接）与绕数计算耦合在一起。本轮按 continuum-first 的思路重构：

- continuum 检测前移，二分只处理离散路径；
- crossing 检测独立成 `find_crossings`；
- 绕数计算改成薄封装 `calculate_a2_average_winding`；
- subset 组装放到二分结束后；
- 默认 Laurent backend 改为 numpy。

## 主要改动

### 1. `amoeba/zm_extract.py`

- 删除旧 subset 提取层：`extract_amoeba_subsets`、
  `_continuum_mask`、`_finalize_crossing`、`_zero_identity_key`、
  `_join_continuum_across_mrs`、`_merge_two`、`_LinePiece` 别名等。
- 保留并更新：
  - `AmoebaZeroManager`：新增 `refresh_logabs()`，供插点后刷新缓存。
  - `detect_continuum(zm, tol)`：std 版 continuum 检测，返回
    `[(mu2_c, [(seg_idx, col_idx), ...]), ...]`，近邻 `mu2_c` 合并。
  - `find_crossings(zm, mu1, mu2, avoided_segments=None,
    return_refined=False)`：统一 crossing 检测；支持跳过指定
    `(seg_idx, col_idx)`；`return_refined=True` 时用
    `_find_exact_crossing` 精修 sign-change crossing，exact touch 不精修。
  - `calculate_a2_average_winding(zm, mu1, mu2, ...)`：由
    `find_crossings` 结果计算 a2 平均绕数。

### 2. `amoeba/bisect.py`

- 删除 `_resolve_continuum`、`_refine_and_correct`（Newton 精修）。
- 新增：
  - `_try_fast_mu2`：θ₁=0 处按 |β₂| 排序取前 M / 后 N 列，插点精化
    极值后计算 A/B gap，有 gap 直接返回 μ₂=(A+B)/2。
  - `_screen_extremum_intervals` / `_refine_track_extrema`：向量化初筛 +
    Hermite 预测 + `insert_solution` 插点，借鉴 SGBZ 网格加细思路。
  - `_bisect_mu2_discrete`：纯离散 μ₂ 两阶段二分（粗：不 refine +
    `BISECT_COARSE_XTOL`；细：refine + `BISECT_XTOL`，细 bracket 自动
    扩展直到 refined 端点反号）。
  - `_handle_continuum`：μ₁ ± ε 处各建 ZM、用 `find_crossings` 找零点
    后算 w1 极限，判断是否 continuum boundary。
- `_find_mu2_for_w2_zero` 重写为 continuum-first 调度器：
  `_try_fast_mu2` → `detect_continuum` → continuum ε 探针 →
  `_bisect_mu2_discrete`。
- `bisect_amoeba_ronkin_min`：所有 inner 返回均显式处理
  `is_continuum == True`（端点扩展与主循环中点）。

### 3. `amoeba/amoeba.py`

- 新增 subset 组装：
  - `_assemble_discrete_subsets`：zeros → PointSubset。
  - `_assemble_continuum_subsets`：continuum `(seg_idx, col_idx)` →
    per-segment `JoinableLinePiece` → `_splice_continuum_pieces` 首尾
    拼接检查 → 其余 track 用 `find_crossings(avoided_segments=...)`
    找离散点 → PointSubset；距 LineSubset θ₁ 小于 `SNAP_TOL` 的点剔除。
  - `_splice_continuum_pieces` / `_splice_two`：MR 列身份 + boundary_perm
    接缝拼接。
- `collect_GBZ_subsets` 恢复完整流程（bisect → plateau check → 组装）。

### 4. `amoeba/ronkin_winding.py`

- 删除未使用的 `_compute_zero_dtheta1_dmu2`（原 Newton 精修用）。

### 5. `continuation/arclength.py`

- `predict_roots` 在 `np.exp` 前先判断
  `|Vⱼ·Δθ₁| > PREDICT_MAX_ABS_ARG`，避免溢出；新增常量
  `PREDICT_MAX_ABS_ARG = 1.0`。

### 6. `backend.py`

- 默认 backend 改为 `NumpyLaurent`；`PolyToolsLaurent` 仅显式请求时
  使用，并在类旁备注 known issues（可能漏掉部分重根）。

### 7. 文档 / 测试 / 引用清理

- 重写 `doc/amoeba.md`；更新 `doc/constants.md`、`doc/SGBZ.md` 中
  过时引用。
- 更新 `README.md` 的 amoeba API 表与 backend 说明。
- 更新 `CLAUDE.md` 的项目结构描述。
- 重写 `tests/test_amoeba.py` 适配新 API；更新 `tests/test_constants.py`
  与 `tests/test_backend.py`。
- 更新 `debug_tool/gbz_debug.py` 与旧 playground 脚本中的
  `amoeba_windings` / `extract_amoeba_subsets` 引用。

## 验证

- `pytest -q tests/test_amoeba.py tests/test_constants.py tests/test_backend.py`
  全部通过。
- `playground/test_new_amoeba_bisect.py` 冒烟通过：
  E=0 / E=1 返回 continuum LineSubset，E=5 返回空 GBZ。
