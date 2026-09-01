# TODO: 精确简并输入下的判定边界行为（touch/cross 通道与 MR 重启簿记）

来源：Phase 2 backend 整合时 `POLY_BACKEND=numpy` 全套测试暴露（2026-08）。
这不是 numpy 后端的 bug——两个后端**最终 GBZ 输出完全一致**（μ₁★ 逐位相同、
同样的 LineSubset）；暴露的是算法在"精确简并恰好落在判定边界上"时，行为由
末位 ULP 噪声决定。更高精度算术（quad）会把同类问题放大，此项是前置排查。

## 现象 1：seam 精确交叉被分类为 touch（pairwise.py）

- 模型：HN2D 实跃迁（δ₁=δ₂=0, basis "10"），E=1.0, μ₁=0.1。对称性使两条
  representative 轨道的 ln|β₂| 在 seam（θ=0≡2π）**严格相等**。
- `collect_pair_events` 的 touch 通道由 `d[i] == 0.0` **逐位判定**触发；
  cross 通道要求严格变号 `d[i]·d[i+1] < 0`。
- 后端系数末位舍入不同 → 一个走 touch（θ*=0），一个走 cross（θ*=2π）。
  方向、电荷、最终结果完全一致；只有 kind 标签和 θ* 簿记不同。

### 根因与改进方向

touch 的语义本应是"方向不可信"（MIN_DIRECTION_DERIV 保护），但触发条件是
"逐位相等"，与方向可信度无关。**精确落在网格行上的横向交叉**应判 cross：

- 行 i 处 d[i]==0 时看两侧符号：`d[i-1]·d[i+1] < 0` → 横向交叉（cross，
  用切向方向）；同号 → 真切向 touch。seam 行 i=0 需跨段环绕取邻居。
- 改动面：`pairwise.py` touch 通道分类（约 ±15 行）+ 针对两侧情形的单测。
- 风险：pairwise 是最敏感通道；默认（C++噪声）路径几乎不会命中逐位相等，
  改动对其惰性，全套回归护栏在场。

## 现象 2：非通有二重根的 MR 重复记录（zero_manager.py）

- 模型：poly A（μ₁=0，θ=0 处精确二重根，随 θ 线性慢分裂）。
- 重启行 θ=mr_jump 处真实分裂量低于双精度可分辨度（系数噪声 1e-16 →
  根分裂噪声 ~1e-8 > 真实分裂 ~1e-9）。
- 重启行 ressolve 恰好返回相等对（切向 inf → 点触发）则同一簇被记第二次；
  返回噪声分裂对（切向大而有限）则不会。计数由噪声决定。

### 根因与改进方向

重启行落在 MR 简并邻域内时，"再次检测到同一簇"应识别为**同一 MR 的重检**
而非新记录：

- 新 MR 的 θ 在 `[prev.θ, prev.θ + mr_jump_eff]` 内且 cluster 成员相同 →
  视为重检测：加大重启距离重试，而不是 append（现有 `_MR_STUCK_TOL` 只防
  "无前进"，不防"近距离重检"）。
- 改动面：`zero_manager.py` run 循环 MR 分支（约 +20 行）+ 非通有模型单测。

## 验收标准

- [x] 现象 2 已修复（2026-09-01）：`test_zero_manager.
      test_poly_A_no_generic_mr` 恢复 `== 1`，且 numpy / poly_tools 两个
      backend 下均通过。`zero_manager.py` 新增 same-cluster re-detection
      守卫（`MR_REDETECT_RETRY_FACTOR` / `MR_REDETECT_MAX_RETRIES`）。
- [ ] 现象 1 仍未修复：`test_sgbz.test_mu01_has_single_seam_event_group`
      仍为 `ev.kind in ('cross', 'touch')`，待改 pairwise touch 通道。
- [ ] 默认全套 + `POLY_BACKEND=poly_tools` 全套绿（未重新全量跑）。
