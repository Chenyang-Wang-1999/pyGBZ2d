# 2026-08-15 — SGBZ continuum cluster 判据改为列级投票

用户审查指出：zero curve 的复解析性保证的是**相同 column index** 的两条
曲线若在一段同模则整条同模（至 MR）；但数据中 M-1/M 是模长排序位，不是
列身份。因此 continuum_clusters 的构建应在列对上用投票判据，而不是
whole-segment max。

## 改动

### `brute_force_SGBZ/mu2mid.py`

- 删除 `_INLINE_CONTINUUM_FRAC = 0.1` 常量。
- `_detect_continuum_clusters_internal` 的列对判据从
  `max over rows |ln|β_j| − ln|β_k|| < tie_tol` 改为：
  ```python
  frac_in_band = np.mean(np.abs(logabs[:, j] - logabs[:, k]) < tie_tol)
  if frac_in_band > CONTINUUM_FRAC:   # 0.9，与 amoeba _continuum_mask 一致
  ```
  旧 max 判据过保守：near-MR / seam 附近单个噪声行会拆散真实 continuum
  cluster。投票判据允许少数异常行，但仍排除只在一个孤立 θ 点同模的
  transversal 0D 点（它们 frac ≈ 0）。
- 删除 `_detect_continuum_inline` 方法；`build_mu2_mid` 中
  `has_continuum` 直接由 ItemView 判断：
  ```python
  self.has_continuum = any(
      bool(np.any(view.j_lo == view.j_hi)) for view in self._item_views
  )
  ```
  j_lo == j_hi 即 boundary pair 落入同一 continuum item，无需二次扫描。
- `CONTINUUM_FRAC` 注释更新：不再是 legacy-only，而是 SGBZ 列级同模投票
  分数（demo 仍继续使用）。

### 文档

- `doc/SGBZ.md`：inline continuum detection 段落与 §4 `CONTINUUM_FRAC`
  参数表同步（删除 `_INLINE_CONTINUUM_FRAC` 描述）。

## 验证

- `tests/test_sgbz.py` 全量回归通过（28 passed）；其他测试无改动。
