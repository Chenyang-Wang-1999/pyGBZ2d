# SGBZ / amoeba：列身份 join、padding 边界拒绝与 component 粒度 MR 标注

日期: 2026-08-25

## 背景

本批修复三个 review 残留：

1. padding β₂ root 可以参与 ItemView 全局排序，但不应成为 M-1/M
   boundary pair；
2. `_is_cluster_endpoint` 用最近复数值判断 MR cluster membership，在高简并
   /浮点误差下可能错配；
3. EventGroup 的 `is_mr` 是 group 粒度，同一 θ 上无关的 ordinary crossing
   会被 MR component 连带标成 hard。

## 1. padding root 进入 M-1/M 时显式失败

`Mu2MidZM._build_item_views()` 仍用 raw `ln|β₂|` 排序：

- `β₂ = 0` → `-∞`，位于最外侧；
- `β₂ = ∞` → `+∞`，位于最外侧。

这保持不变。新增检查：若任一 mesh row 的 `j_lo/j_hi`
（排序位置 M-1/M）对应的 `item_logabs` 非有限，则直接
`ValueError`，错误信息包含 segment / row / item 列。

理由：此时 SGBZ 边界模长本身是 0 或 ∞，有限 μ₂_mid 路径无法表示；
继续用 ±14 clamp 只会得到物理上不可靠的有限代理。批量
`collect_GBZ_subsets` 会把该异常转换为 `success=False` 的 `GBZResult`。

## 2. MR cluster membership / continuation 改用列身份

### `_is_cluster_endpoint`

API 从：

```python
is_mr_cluster_endpoint(zm, seg, side, root)
```

改为：

```python
is_mr_cluster_endpoint(zm, seg, side, col)
```

语义：

- interior MR：`cluster_indices` 与两侧 boundary row 共享 track frame，
  直接检查列号；
- θ=0 左边界 MR0：直接检查列号；
- θ=2π 右边界复用 MR0：MR0 的 `cluster_indices` 在 θ=0 模排序帧中，
  用 `boundary_perm` 翻译后检查。

约定仍为：

```text
roots_right[boundary_perm] == roots_left
```

因此左列 `k` 对应右列 `boundary_perm[k]`。

### SGBZ / amoeba join 调用方

`_join_runs_across_mrs()` 与 `_join_continuum_across_mrs()` 不再用
`argmin |root_i - root_j|` 找上一段右边界 continuation：

- interior MR：列号相同；
- first/last segment cyclic seam：左列 `j_l` 对应右列
  `boundary_perm[j_l]`。

LinePiece 的端点定位仍用根值筛选，但 MR membership 与跨边界 continuation
不再依赖最近值。

## 3. MR hard 标注细化到 component / column

`EventGroup` 新增：

```python
column_kind: dict[int, str]
```

`finalize_event_groups()` 对每个 connected component 判断：

```python
component_mr = any(e.is_mr for events in this component)
```

然后逐 real column 写入：

- `'mr'`：该 component 含 MR event；
- `'tangent'`：非 MR 但 direction/charge unknown；
- `'ordinary'`：可定荷事件。

`detect_crossings_simple()` 改为读取 `g.column_kind[col]`，不再用
`g.is_mr` 覆盖整组。缺失 `column_kind` 时 raise，暴露未 finalize 的对象。

保留的保守行为：若 MR event 与 ordinary event 在 item 图中连通，它们属于
同一 component，仍整体 hard。这种情形表示拓扑上确实耦合，宁可多独立计算
一条 loop winding。

## 测试

新增 / 更新：

- `tests/test_sgbz.py::TestItemViewBoundaryPadding`
  - 外侧 padding root 正常排序；
  - padding root 进入 M-1/M 显式 raise。
- `tests/test_sgbz.py::TestEventGroupColumnKind`
  - `finalize_event_groups()` 对同一 group 内两个 disconnected components
    分别写出 MR / ordinary kind 和 charge；
  - materialization 对应输出 `kind='mr'` / `kind='ordinary'`。
- `tests/test_gbz_types.py::TestIsMrClusterEndpoint`
  - 左 / interior 右边界列号直接判定；
  - θ=2π boundary MR 右列经 `boundary_perm` 翻译；
  - interior MR index 0 不翻译；
  - 缺少 `boundary_perm` 时 raise。
- `tests/test_regressions.py` 的合成两段 topology 补充 seam monodromy。
- 原最近值匹配测试改为列身份测试。

回归：

- `test_gbz_types.py + test_continuation.py + test_zero_manager.py +
  test_regressions.py`：132 passed；
- `test_sgbz.py`：48 passed；
- `test_amoeba.py`：17 passed；
- 全量 `pytest tests/`：227 passed, 3 skipped。
