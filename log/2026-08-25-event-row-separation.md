# SGBZ pairwise：event 行全集 + 相邻 event 间 regular 行

日期: 2026-08-25

## 背景

`finalize_event_groups` 的拓扑荷 side-change 规则要求读每个 EventGroup
左右两侧**最近的 regular 行**。原实现只有两个机制：

1. `group_events` 合并间距 `< crossing_tol` 的事件；
2. `insert_event_groups` 每个 EventGroup 只插一行。

但没有保证两个相邻 event 行之间一定存在 regular 行，也没有把合并后
group 内原始的 `d==0` touch 行保留在 event 集合里。事件密集时，左右
regular 搜索会跨越整簇事件，导致 side/charge 读到错误侧。

## 修改

- `pairwise._event_thetas_for_segment`：一个 segment 的 event θ 全集 =
  EventGroup 代表行 θ + 每个 `kind == 'touch'` 事件的原始 `theta_star`。
  deferred seam group 用 `TWO_PI` 代表。
- `pairwise._insert_regular_rows_between_events`：相邻 event θ 之间若无
  任何 mesh 行，则在**中点直接求解多项式**插入一个 regular 行；
  `insert_solution(interp='hermite')` 的 `interp` 只选择 track 匹配锚，
  根值本身总是真实求解。
- `pairwise._resolve_event_rows`：event 行集合统一由 group 行 + touch
  原始行构成；`finalize_event_groups` 使用同一集合。
- `pairwise._check_event_row_separation`：显式不变量——同一 segment 内
  任意两个 event 行不得相邻，违反即 `RuntimeError`。
- 交叉事件合并进 group 后**不单独插成员 θ***；group 行仍是唯一代表。
- 插入 regular 行后不重扫交点（一轮完成）。

## 测试

- `test_touch_member_rows_stay_in_event_set_after_grouping`
- `test_regular_row_is_inserted_between_adjacent_events`
- `test_event_row_separation_invariant_raises`
- `test_analyze_has_a_regular_row_next_to_every_event`
- 全量 `pytest tests/`：206 passed, 3 skipped。
