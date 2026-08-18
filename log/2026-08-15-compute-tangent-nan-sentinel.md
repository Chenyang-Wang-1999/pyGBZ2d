# 2026-08-15 — compute_tangent: 0/∞ padding root 的切线从 0 改为 nan

用户审查指出：`continuation/arclength.py:compute_tangent` 把 0/∞ padding
root 的切线 V 设为 0，而 0 表示“根不随 θ₁ 移动”，但 padding root 的
d(ln β₂)/dθ₁ 实际**无定义**。本轮将 sentinel 改为 nan，并同步修正
ZeroManager 中的过时注释。

## compute_tangent 行为澄清

- 0/∞ padding roots（`_is_singular_root`）→ `V = nan`（undefined）；
- multiple root（∂f/∂β₂ == 0）→ `V = inf`（divergent，保持不变）；
- 其他非有限结果 → `inf`。

`norm_V` 改为 `sqrt(1 + Σ|V[~isnan(V)]|²)`：只忽略 nan，保留 inf，因此
MR 仍使 `dθ₁ = h/norm_V → 0` 并触发步长塌缩 MR 检测，而 padding 的
undefined 切线不会污染步长控制器。

## 注释/文档同步

- `continuation/arclength.py:compute_tangent` docstring。
- `continuation/zero_manager.py`：
  - `SegmentData.tangents` 字段注释（旧注释错误声称 MR cluster tracks
    V=0，现改为 nan=padding undefined / inf=MR divergent / finite=真实切线）；
  - `_compute_tangents` docstring；
  - `run()` 中 boundary-MR 预测处 "zeroed tangent" → "divergent (inf)
    tangent"。
- `doc/continuation.md` §2.3：V_j = nan 语义、norm 对 nan/inf 的处理，
  并顺手修正阈值笔误（10^-14/10^14 → 1e-6/1e6，与代码一致）。

## 测试

- `tests/test_continuation.py::test_finite_difference_agreement`：padding
  root 断言从 `abs(V) < 1e-10` 改为 `np.isnan(V)`。
- 新增 `test_padding_roots_get_nan_and_norm_ignores_nan`：roots 末尾拼
  0/∞ 后 V 对应位置为 nan，norm_V 与仅含有限 roots 的 norm 一致。

## 回归

test_continuation 51 passed；其余回归（test_sgbz / test_zero_manager /
test_amoeba）待跑。
