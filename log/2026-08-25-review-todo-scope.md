# Code Review TODO 范围调整

日期: 2026-08-25

## 决定

按当前维护范围决定，以下 review 疑点**不再作为待办事项跟踪**。这不是
“已修复”，而是经讨论后暂不处理，从
`TODO/review-2026-08-18-remaining.md` 中移除：

1. 原 TODO 疑点 2 —— MR 恰在 β₂ = 0 时双触发器联合失明。
2. 原 TODO 疑点 3 —— `solve_multiple_roots_iterative` 对
   `|β₁|` 无 μ₁ 圆约束。
3. 原 TODO 疑点 5a —— SGBZ 求解器中“W 随 μ₁ 单调递增”的理论假设。
4. 原 TODO 疑点 6b —— amoeba continuum 检测中
   `CONTINUUM_FRAC = 0.9` 的阈值脆性。

TODO 中保留的编号仍沿用原始 review 编号（1、4、5、6），便于与
`diagnostics/code-review-2026-08-18.md` 对照；编号缺口表示上述条目已按
维护范围移除。

## TODO 中的保留部分

- 疑点 1：MR tangent 列帧 / boundary MR 复用索引的潜在错标。
- 疑点 5 中仍保留：
  - boundary MR 邻近事件是否还有误删风险；
  - padding root 进入 ItemView 排序。
- 疑点 6 中仍保留：
  - amoeba `mu2_mid` 与 `w_left` 的 ε 级不一致；
  - `_is_cluster_endpoint` 最近值匹配的错配风险；
  - Newton 修正 clamp 到初始 bracket。
- 2026-08-18 pairwise 重构后的 group 粒度 MR hard 标注问题。

## 说明

这些移除不改变代码行为，也不表示对相应机制作出理论结论。若后续需要
重新审查，可从本 log 恢复上下文，再按原始 review 编号追查。

## 后续更新（同日）

本节记录的是当时 TODO 的保留范围，其中部分项已在后续提交中处理：

- 疑点 1 / boundary MR tangent 列帧：见
  `2026-08-25-boundary-mr-tangent-frame.md`；
- 疑点 5 / padding root 进入 M-1/M、group 粒度 MR hard 标注，以及
  `_is_cluster_endpoint` 最近值匹配风险：见
  `2026-08-25-column-identity-and-hardening.md`。

截至该后续更新，amoeba 侧仍待处理的是 `mu2_mid/w_left` 的 ε 级不一致
与 Newton 修正 clamp 到初始 bracket。
