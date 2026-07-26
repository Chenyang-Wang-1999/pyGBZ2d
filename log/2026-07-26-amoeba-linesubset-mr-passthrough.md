# Amoeba LineSubset 端点修正：continuum 穿过 MR 的拼接

日期：2026-07-26
模块：`brute_force_amoeba/zm_extract.py`、`tests/test_amoeba.py`

## 问题

`extract_amoeba_subsets` 对每个 segment 的每条 continuum track 直接生成一个 `LineSubset`，端点就取 segment 端点（即 MR 所在位置）。但 MR 处只有 `cluster_indices` 列出的根才是真简并根，其余根是穿过 MR 的正则根。因此：

- 若一条 continuum track 以**非 cluster 根**穿过 MR，它并不该在那儿终止——目前却被截断，端点是假端点。
- 只有端点根**落在 cluster 里**时，MR 才是真端点。

这是 CLAUDE.md "反常输出是改进算法的机会" 类问题：根因在于把"segment 边界 = MR = LineSubset 端点"当成恒等关系，忽略了 MR 只局部简并。

## 决策与设计

### 判断真端点的方式

按用户给的判断方式：
1. 看 LineSubsets 的左右端点根是否真落在 multiple root 的 cluster 里——是则 MR 是真端点。
2. 否则看相邻 segment 是否也是 LineSubset（同 track 的 continuum）——是则拼起来；否则报错（拓扑不一致）。

### 实现（`zm_extract.py`）

**`_LinePiece`**：`LineSubset` 子类，多带 `ml`/`mr`（拼接后跨过的最左/最右 segment 原始下标），用于链式拼接与环形 seam 检测。joining 完成后行为上就是一个 `LineSubset`。

**`_is_cluster_endpoint(zm, seg, side, root)`**：判断端点根 *值* 是否在 boundary MR 的 cluster 里。
- `mr = left_mr` / `right_mr`；`-1` 或（`mr==0 and not has_boundary_mr`）→ 无 MR，非 cluster。
- 按 *值* 匹配 `multiple_roots[mr].roots`（`argmin |roots - root|`），再查该下标是否在某个 `cluster_indices` 元组里。
- **frame-independent**：boundary MR 的 roots 是 modulus-sorted，interior MR 的 roots 是 track-ordered，但 `cluster_indices` 始终是对同一行 roots 的下标，所以按值匹配在两种 frame 下都对。

**`_join_continuum_across_mrs(zm, continuum_masks, line_pieces)`**：fixpoint 迭代，每个 segment 的左端点（= 上一 segment 的右端点）逐 continuum track 处理：
- **circle seam 识别**：`left_mr < 0 and right_mr < 0`（θ₁=0≡2π，唯一非 MR 共享边界）。注意不是 `mr==0 and not has_boundary_mr`——`has_boundary_mr=False` 时 MR index 0 是普通 interior MR，仍可拼接。
- 端点根在 cluster → 跳过（真终止）。
- 非 cluster → 在相邻 segment 右端点按值找匹配根：
  - 也非 cluster 且有同 track continuum → `_merge_two` 拼接（共享边界行去重；环形 seam 让 segment 0 在前，丢 θ₁=0/2π 重复行）。
  - 相邻无 continuum → `ValueError`（"continuum continues through MR as non-cluster root but no matching continuum"）。
  - 一边 cluster 一边非 cluster → `ValueError`（track 一边真终止一边穿过，不一致）。
- **避免重复/死循环**：已拼好的内部边界，`find_by_left` 返回 None（该 segment 不再是任一 piece 的左端）或 `pi_idx == li_idx`（两端已属同一 piece，如环形 seam 已拼）→ 跳过。
- fixpoint：每次拼接后 `break` 重扫，直到一轮无变化。

**`extract_amoeba_subsets`** 改动：
- continuum LineSubsets 先建成 `line_pieces`（per-segment），joining 后再用于 `cont_endpoints`（Rule 1 snap）与最终 `subsets`。
- `_merge_two` 后 `boundary_perm_inv` 仍需为 Rule 2 dedup 计算（一度被误删，已补回）。

## 改动文件

- `brute_force_amoeba/zm_extract.py`：新增 `_LinePiece`、`_is_cluster_endpoint`、`_join_continuum_across_mrs`、`_merge_two`；`extract_amoeba_subsets` 的 continuum 路径改用 `line_pieces` + joining；返回前 `_LinePiece` 转回 `LineSubset`（不泄漏内部 `ml`/`mr`）。
- `tests/test_amoeba.py`：新增 `TestContinuumThroughMR`（4 例）+ `_build_two_segment_zm` 合成拓扑构造器。

## 环形 seam 拼接方向修正

`_merge_two` 的 cyclic 分支方向一度写反。seam 处 `prev` 是最后一段（右端 θ₁=2π），`cur` 是 segment 0（左端 θ₁=0），共享点是 seam 本身——即 `rp[-1]`（θ=2π）和 `cp[0]`（θ=0）。要让 track 在 seam 处连续，必须让这两点在数组中相邻，即 `[rp, cp[1:]]`（seam 落数组内部，两端落在真正的 interior cluster 终止子上）。

写成 `[cp, rp[1:]]` 会把 segment 0 的右端（某 interior MR）接到最后一段的第二个点——不同的物理点，β₂ 在那里断；且 seam（θ₁=0/2π）被推到数组两端，把 track 连续穿过的地方伪造成端点。已改为 `[rp, cp[1:]]`。

附带：当 continuum 同时穿过 interior MR 和 seam（形成闭环、无真端点）时，正确行为是在 interior MR 处切开（seam 落内部），而非在 seam 处切开。原 `test_continuum_passing_through_mr_joins` 断言"结果端点在 [0,2π]"是错的，已改为"两端在 interior MR θ₁=π、seam 在内部、β₂ 跨 wrap 步长与段内一致"。

## 验证

- `tests/test_amoeba.py`：17 passed（含 4 个新测试）。
- `tests/test_zero_manager.py` + `tests/test_continuation.py`：68 passed（未触 ZeroManager，仅回归确认）。
- 合成验证：
  - 2-seg、continuum 穿过 interior MR（cluster={0,1}，track 2 非集群）→ 闭环，切开在 π，1 条 n=21，seam 处 β₂ 步长与段内一致。
  - seam-only（track 2 在 interior cluster、无 boundary MR）→ seam 拼 1 条，两端在 π，seam 内部。
  - 缺延续（seg1 track 2 离开 continuum）→ raise。
  - track 2 在 cluster + boundary MR（seam 也是 cluster 终止）→ interior MR 拼或不拼视 track 2 是否在 cluster。
  - frame 不匹配（seg0 开口 modulus-sorted col 2、seg1 闭口 track-ordered col 0，同值不同列）→ 按值匹配仍拼。

## 设计要点

- **按值匹配而非 track 下标**：segment 边界行 frame 不统一（boundary MR modulus-sorted、interior MR track-ordered、`completed` 分支闭口行 track-ordered），按值匹配一概适用，无需 `boundary_perm_inv` 做 frame 折叠（该 inv 仅 Rule 2 dedup 还在用）。
- **circle seam 判据**：`left_mr < 0 and right_mr < 0`。`mr==0 and not has_boundary_mr` 是错的——会把 `has_boundary_mr=False` 下的 interior MR index 0 误判为 seam。
- **seam 拼接方向**：`[rp, cp[1:]]`（rp=最后段在右端、cp=seg0 在左端），让 seam（rp[-1]↔cp[0]）落数组内部；反方向会把 interior MR 误当接缝。
- **fixpoint 而非单遍**：链式拼接（>2 segment）与环形 seam 都需多轮收敛；每轮拼一条就 `break` 重扫，避免在变化的 list 上继续遍历。
- **闭环表示**：continuum 无真端点（穿过所有 MR 与 seam）时是闭环，在 interior MR 处切开、seam 落内部；θ₁ 数组会跨 0/2π wrap（非单调），这是表示穿过 seam 的弧所必需。

## 遗留

- 合成测试用 `SimpleNamespace` 代替真 `AmoebaZeroManager`（只测 joining 逻辑，不跑多项式求根）；真物理模型的 multi-segment + continuum-through-MR 用例目前依赖 Haldane demo 的 sweep 间接覆盖，未加针对性单测。
- 若一条 continuum 在 interior MR 处一边 cluster、一边非 cluster，当前直接 raise。理论上是拓扑不一致，但未在实际模型中遇到验证过。
- 闭环 LineSubset 的 θ₁ 非单调，下游消费者（若假设单调）需留意。
