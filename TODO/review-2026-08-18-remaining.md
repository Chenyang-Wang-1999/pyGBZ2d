# Code Review 未竟事项

来源：`diagnostics/code-review-2026-08-18.md` 第六、七、八节
（二~五节已全部修复，2026-08-18 晚；本文件承接剩余项）。

## 文档一致性（第六节）

- [x] CLAUDE.md 更新（2026-08-18）：文件结构与实际对齐——补 mu2mid.py、
      pairwise.py、continuum_lines.py、zm_extract.py、interpolation.py、
      doc/continuation.md、conftest 描述、tests 完整清单；删已不存在的
      continuum.py / crossings.py / tracks.py；修正 "lazy beta2_arr"
      → eager；gbz_types 工具列表补 CharPoly/TWO_PI/circ_dist/
      probe_zero_plateau/is_mr_cluster_endpoint 等（引用文件存在性已
      脚本核验，20/20）。
- [ ] README：引用不存在的 `demos/demo_amoeba.py`、`demo_solver.py`；
      `pip install -e .` 指示必然失败（无 pyproject.toml/setup.py）。

## 工程卫生（第七节）

- [ ] **提交未入库的核心源码**（最高优先）：`brute_force_SGBZ/pairwise.py`、
      `continuation/interpolation.py`、`tests/test_interpolation.py`、
      `TODO/`、`paper/`、6 个 log 文件 untracked；另有 20+ 已跟踪文件
      +1339/−1775 未提交。建议按 log 条目拆分提交。
- [ ] 依赖不可复现：无 requirements.txt / pyproject.toml。
- [ ] 无 CI、无 lint 配置（ruff/pyflakes/mypy 均未装）。

## 疑点（第八节，未确证，需逐项实证）

- [ ] 1. `_mark_mr_tangents_inf`（zero_manager.py:1099-1121）用 θ=0 处
      **模排序帧**的 cluster_indices 标注段边界行 **track 帧**切线为
      inf；boundary_perm 非平凡且 cluster 成员移动时可能标错列。
      `seg.tangents` 被 pairwise.py / mu2mid.py 直接消费。建议在
      `_append_segment` 加 cluster 列与末行实际重合的断言。
- [ ] 2. MR 恰在 β₂=0 时双触发器联合失明（两根同进 |β₂|<1e-6 窗口 →
      V 全 nan → point/interval 触发器同时失效）。
- [ ] 3. `solve_multiple_roots_iterative` 对 |β₁| 无约束，解可漂离 μ₁ 圆，
      仅 warning 兜底。
- [ ] 4. `solve_roots_1d` padding 最多各补一个 0/inf，次数亏缺 >1 时
      返回数组短于 M+N，下游 K 形状校验会崩（gbz_types.py:166-170）。
- [ ] 5. SGBZ "W 单调递增于 μ₁" 假设、boundary MR 邻近整组丢弃可能误删
      合法 crossing、单网格区间双穿越漏检、padding 根进入 ItemView 排序。
- [ ] 6. amoeba：`mu2_mid` 配对 `w_left` 的 ε 级偏移、continuum frac=0.9
      阈值脆性、`_is_cluster_endpoint` 值匹配错配风险、Newton clamp 到
      初始括区间。

## 本轮新增残留（2026-08-18 pairwise 重构后）

- [ ] MR 单通道化后 `is_mr` 是 **group 粒度**的 hard 标注（charge=None）；
      若 MR 边界行上存在可定荷的横向穿越（非 cluster 根的模长重合），
      需要按事件粒度细分 hard 标注，否则该穿越荷被吞为 None。


## 新发现的 bug
- [ ] Gain-loss Haldane model 在 E = 1.2890000000000001-0.4131j 时，在 mu1 从 -0.0002 变化到 0.0002 时，SGBZ 平均绕数有反常跳变。同时查看 mu2_mid，并没有经过应该经过的交点。
