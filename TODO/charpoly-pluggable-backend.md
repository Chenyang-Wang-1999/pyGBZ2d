# TODO: CharPoly 可插拔 backend（用户免装 poly_tools）

状态：方案已验证（原型对拍通过），**未实施**。2026-08-26 方案评审时定：
先做 bfGBZ2d 打包工程化（见 TODO-2026-08-26.md 与 packaging 方案），backend
机制随打包一并整合进包内。本文件保留实施依据与验收标准。

## Motivation

- `poly_tools` 是需单独编译的 C++ 扩展（README 要求用户手动 `make` 并拷贝
  `.so`），是"开箱即用"的最大障碍；它还在包顶层 import 了 sympy。
- 好消息：全项目仅 `gbz_types.py` 一处 `import poly_tools`，没有任何外部
  代码触碰 `CharPoly` 的私有属性；C++ API 实际用到只有 6 个
  （`set_Laurent_by_terms` / `derivative` / `eval` / `partial_eval` /
  `num.batch_get_data` / `denom_orders`）。

## 原型实测结论（2026-08，已验证可行性）

纯 numpy 后端 vs C++ poly_tools（HN2D + 随机 24 项 Laurent）：

| 指标 | 结果 |
|---|---|
| M/N 次数、eval、一阶/二阶偏导 | 全部一致，rel-err ~1e-14 |
| `solve_roots_1d` 根 | abs-err ~1e-13，计数一致 |
| `solve_roots_1d` 速度 | 仅慢 ~1.25×（np.roots/LAPACK 主导） |
| `eval_partials` 速度 | 慢 ~6.5×（低频路径，可接受） |
| `__init__` 速度 | 慢 ~4×（一次性构造，可接受） |

## 实施要点（实施时不得丢的三个语义细节）

1. `denom_orders[d] = max(0, -min_deg_d)` —— Laurent 分母只清负幂；
2. **同次项必须合并**（C++ 链表会合并；`solve_roots_1d` 的系数组装循环按
   唯一次数逐项赋值，不合并会静默丢项）；
3. 系数恰为 0 的项要丢弃（否则破坏"多项式完全消失"空分支）。

## 验收标准

- [ ] `poly_backend.py`：`LaurentProtocol` + `PolyToolsLaurent`（延迟
      import）+ `NumpyLaurent` + `make_laurent`（backend=None →
      `GBZ_BACKEND` 环境变量 → poly_tools → numpy 回退并警告一次）；
- [ ] `CharPoly(coeffs, degs, backend=None)` 接受自定义 backend；
- [ ] 双后端对拍测试进 `tests/`（eval/偏导/二阶偏导/solve_roots_1d，
      含 HN2D-10、HN2D-11、随机 Laurent）；
- [ ] `GBZ_BACKEND=numpy pytest` 全绿；
- [ ] 下游 SGBZ / amoeba / continuation / debug_tool 零改动；
- [ ] README 更新依赖说明（numpy+scipy 即可运行，poly_tools 可选加速）。

## 已知代价（接受）

- numpy 后端结果与 C++ 不逐位相同（求和顺序，~1e-14）；
- `np.roots` 在简并处的根序可能不同 —— 下游全部 Hungarian 匹配，天然容忍。
