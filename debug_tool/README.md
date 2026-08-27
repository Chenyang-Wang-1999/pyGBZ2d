# debug_tool — 固定 (E_ref, mu1) 的 GBZ / winding 调试工具

求解入口 `collect_GBZ_subsets` 会自行解出 mu1，结果异常时无法查看"求解器在候选 mu1 处到底看到了什么"。
本工具把 **mu1 冻结**，在给定 (E_ref, mu1) 处同时计算 amoeba GBZ 与 SGBZ 两套子集、SGBZ 拓扑荷，
以及各 theta2 上的 winding loop 绕数，并画在 θ₁–θ₂ 图上。

## API（两个核心函数）

### 1. `collect_debug_subsets(poly, E_ref, mu1) -> GBZDebugReport`

一次调用返回两种方法在固定 (E_ref, mu1) 的全部 GBZ subsets：

```python
from pygbz2d import CharPoly
from debug_tool import collect_debug_subsets, plot_winding_debug

poly = CharPoly(coeffs, degs)
report = collect_debug_subsets(poly, E_ref=1.212 + 0j, mu1=0.135328598265)
print(report.summary())

report["sgbz"]    # MethodDebug: subsets + charges（PointSubset 的拓扑荷）
report["amoeba"]  # MethodDebug: subsets + mu2 + w1/w2
```

- `report[method].subsets` — `PointSubset` / `LineSubset` 列表（离散 0D / continuum 1D）。
- `report["sgbz"].charges` — 与 PointSubset 一一对应的字典
  `{"theta1", "theta2", "charge", "kind"}`；`charge` 为 `+1/-1`（ordinary）、
  `None`（mr/tangent 硬边界，荷未知）。
- `report["sgbz"].W` — 平均主轴 winding `W(E_ref, mu1)`（continuum 时为 `None`）。
- `report["amoeba"].mu2 / .W / .w_secondary` — amoeba 在该 mu1 解得的 w2=0 的 mu2、
  w1 以及 w2。
- 单方法失败不会抛出，错误记录在 `report[method].error`，另一方法的证据保留。

### 2. `compute_loop_windings(md, theta2_list) -> MethodLoops`

对给定 theta2 列表计算 winding loop 绕数（loop = β₁ 绕整圆 θ₁∈[0,2π)、
β₂ 取该方法自己的路径），并自动做 **荷一致性检查**：

- SGBZ：loop 沿 μ₂_mid(θ₁) 逐段积分 `Im[f'/f]`（与求解器同一代码路径）。
  相邻两 loop 的绕数跳变必须等于区间内 SGBZ crossing 的荷之和
  （含硬边界的区间自动跳过——其荷未知是设计行为）。
- amoeba：β₂ 固定在 `exp(mu2 + iθ₂)`，绕数用辐角原理精确计数
  （|β₁|<e^{μ₁} 内的 β₁ 根数 − M₁），并用独立的 `Im[f'/f]` 数值积分交叉验证；
  离散情形下 |绕数跳变| 必须等于区间内的 PointSubset 个数。

每个 loop 附带可靠性指标：`min_abs_f`（loop 上 min|f|）、amoeba 的 `margin`
（根模与 μ₁ 的最小距离）与 `winding_quad`；不可靠的 loop 标记
`reliable=False`（例如 loop 恰好穿过 crossing 的 θ₂）。

当所有 gap check 通过时，还会给出 `ml.W_profile` —— 用电荷传播重建的
u(θ₂) 弧长加权平均，可直接与求解器自己的 `W`（`md.W`）对照；
二者不一致说明 crossing 检测或荷分配有问题。

```python
from debug_tool import compute_loop_windings, auto_theta2_grid

grid = auto_theta2_grid(report["sgbz"])        # 自动取子集 θ₂ 间隙的中点
ml = compute_loop_windings(report["sgbz"], grid)
print(ml.summary())                            # 各 loop 绕数 + gap check 汇总
```

### 3. `plot_winding_debug(target, theta2_list=None)`

把一切画在 θ₁–θ₂ 环面图上（`target` 为 `GBZDebugReport` / `MethodDebug` / 二者序列）：

- SGBZ PointSubset：marker 形状编码电荷（▲=+1，▼=−1，■=硬边界/未知），并加文字标注 `+1/−1/?`；
- amoeba PointSubset：空心圆；LineSubset：按方法着色的曲线（自动处理 θ₁/θ₂ 环面接缝拆段）；
- winding loop：水平虚线，线上标注 `w=±n`（SGBZ 标在线上方、amoeba 在下方，避免重叠）；
  不可靠 loop 退化为点线并加 `?`；
- MR（multiple-root）θ₁ 行：紫色点划竖线。
- `theta2_list=None` 时自动用 `auto_theta2_grid` 取各方法最宽间隙的中点。

## 命令行 demo

```bash
# Haldane gain-loss 超胞，E=1.212、mu1=0.1353...（log/2026-08-16 的调试点）
python debug_tool/demo_debug_tool.py

# 2D Hatano-Nelson 链对（快；解析答案 mu1 = gamma1 = 0.2 为 continuum 情形）
python debug_tool/demo_debug_tool.py --model hn2d --mu1 0.2     # continuum
python debug_tool/demo_debug_tool.py --model hn2d --mu1 0.25    # 离散点

# 手动指定 loop theta2、只跑单一方法、保存路径
python debug_tool/demo_debug_tool.py --theta2 0.5 2.5 4.5 --methods sgbz --save fig.png
```

图默认存到 `debug_tool/out/`；加 `--show` 弹窗显示。

## 测试

`tests/test_debug_tool.py`：HN-2D 模型（快）覆盖离散荷一致性、continuum 情形、
失败捕获、绘图产物；`--run-slow` 额外跑 Haldane E=1.212 调试点
（14 个 crossing、8 项 gap check 全部通过）。

## 设计说明

- **不修补模块输出**：gap check 出现 `mismatch`、amoeba 根计数与积分不一致、
  loop 不可靠等，本身就是该暴露的 bug 信号，本工具只呈现不掩盖。
- SGBZ loop 复用 `pygbz2d.sgbz.winding._loop_winding_quad / _loop_min_f`，
  amoeba 计数与 `_get_average_winding_from_zeros` 同一公式 —— 调试数字与求解器
  内部使用的严格可比。
