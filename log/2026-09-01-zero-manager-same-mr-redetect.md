# ZeroManager 同簇 MR 重检修复

日期: 2026-09-01

## 背景

`TODO/exact-degeneracy-boundary-behavior.md` 现象 2：非通有二重根模型
（poly A）在重启行仍处于同一 MR 简并邻域内时，会把同一个 cluster 再次
append 成新的 MR 记录，导致 `n_multiple_roots` 由末位 ULP 噪声决定。

## 修改

`src/pygbz2d/continuation/zero_manager.py`

- 新增常量：
  - `MR_REDETECT_RETRY_FACTOR = 2.0`
  - `MR_REDETECT_MAX_RETRIES = 4`
- 新增 `_same_mr_cluster(cluster_a, cluster_b)`：对 cluster 列表做
  顺序无关比较（`frozenset(frozenset(c) ...)`）。
- run 循环 MR 分支新增 same-MR re-detection 守卫：
  - `left_mr >= 0`；
  - 新 MR 的 θ 在 `prev.θ + mr_jump_eff` 内；
  - cluster 成员与 `prev_mr.cluster_indices` 相同。
  满足时**不 append MR、不 append segment**，而是按
  `mr_jump_eff * (MR_REDETECT_RETRY_FACTOR ** retry_count)` 推远
  重启点，重新求解并匹配到 prev_mr 的 track frame 后重试；超过
  `MR_REDETECT_MAX_RETRIES` 则抛 RuntimeError。
- 原 forward-progress 守卫保留，负责不同 cluster / 无前进的异常。

`tests/test_zero_manager.py`

- `test_poly_A_no_generic_mr` 恢复严格断言：`n_multiple_roots == 1`。

`TODO/exact-degeneracy-boundary-behavior.md`

- 验收标准更新：现象 2 勾选完成，现象 1 仍未修复。

## 验证

- `pytest -q tests/test_zero_manager.py` — 36 passed。
- `pytest -q tests/test_sgbz.py` — 48 passed。
- `pytest -q tests/test_amoeba.py tests/test_constants.py tests/test_backend.py`
  — 48 passed。
- `POLY_BACKEND=poly_tools pytest -q tests/test_zero_manager.py::TestMultipleRootDetection::test_poly_A_no_generic_mr`
  — passed（与 numpy backend 一致）。
