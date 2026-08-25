# ZeroManager：修复 2π boundary MR tangent 的列帧错位

日期: 2026-08-25

## 背景

当 θ₁ = 0 处存在 boundary MR 时，`ZeroManager.run()` 的记录顺序是：

1. 在 θ=0 解 β₂ 根；
2. 按 `|β₂|` 模长排序；
3. 对模排序后的根做 `detect_cluster` / `snap_clusters_to_mean`；
4. 记录 `multiple_roots[0]`。

因此 MR0 的 `cluster_indices` 是 **θ=0 模排序帧**（即
`left_boundary_roots` 帧）的列号。

随后积分进入段内 track frame。若在 θ≈2π 处再次遇到同一个 boundary MR，
最后一段 closing row 是 **该段 track frame**，但 `right_mr` 仍指向 MR0。
原 `_mark_mr_tangents_inf()` 直接把

```python
multiple_roots[0].cluster_indices
```

当作 closing-row 列号使用。`boundary_perm` 非平凡时，这会标错列：
下游 `pairwise` 的 direction 判定与 `Mu2Mid` 的 event-knot 导数都直接消费
`seg.tangents`。

`boundary_perm` 的类约定是：

```text
roots_right[boundary_perm] == roots_left
```

所以左帧 / MR0 记录列 `k` 对应右 closing-row track 列：

```text
boundary_perm[k]
```

## 修改

### 1. 2π boundary MR 分支先建立 monodromy

`run()` 中 θ≈2π 的 boundary-MR 分支原来先 `_append_segment()`、后计算
`boundary_perm`。现在改为：

1. `_predict_roots_at_2π()`；
2. `_boundary_perm_from_right()`；
3. `_append_segment(..., right_mr=0)`。

原因是 `_append_segment()` 内部会立即调用 `_mark_mr_tangents_inf()`；
没有 `boundary_perm` 就无法翻译 MR0 的 cluster 列。

附近注释明确记录了 θ=0 模排序帧与 2π track frame 的区别，以及为什么
此处必须先做 `boundary_perm`。

### 2. `_mark_mr_tangents_inf()` 增加 seam 列翻译

仅对以下情形翻译：

```python
row == n - 1
and mr_idx == 0
and has_boundary_mr
```

即最后一段右边界复用 θ=0 boundary MR。此时：

```python
right_col = boundary_perm[left_col]
```

再置 `tangents[-1, right_col] = inf`。

不翻译的情形：

- 第一段左边界 MR0：row 0 本身就是 `multiple_roots[0].roots` 的左帧；
- 普通 interior MR：其 `cluster_indices` 与对应边界 row 共享 track frame；
- `has_boundary_mr=False` 时 MR index 0 是普通 interior MR，不是 seam
  boundary MR。

若需要 seam 翻译但 `boundary_perm` 尚未设置，立即 `RuntimeError`，防止
append 顺序 regression。

## 测试

新增 `tests/test_zero_manager.py::TestBoundaryMrTangentFrame`：

- 非平凡 `boundary_perm=[2,0,3,1]`、左帧 cluster `(1,2)`：
  - 右边界应标 inf 的列是 `{0,3}`，不是 `{1,2}`；
  - 左边界仍直接标 `{1,2}`；
- `has_boundary_mr=False` 时，interior MR index 0 不做 seam 翻译；
- 右 boundary MR 需要翻译但缺少 `boundary_perm` 时 raise。

回归：

- `pytest tests/test_zero_manager.py`：36 passed；
- `pytest tests/test_sgbz.py tests/test_regressions.py`：56 passed。
