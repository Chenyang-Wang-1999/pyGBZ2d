# 会议记录：已解决 TODO 归档（backend 整合 · 打包工程化）

日期: 2026-08-26
性质: 工作总结会（本日多轮方案评审的决议汇总；执行细节见同日两份
工程日志，本记录只保留决议与关键结论）

---

## 议题一：CharPoly 可插拔 backend（原 TODO/charpoly-pluggable-backend.md，已整件完成）

### 动机

`poly_tools` 为需手动编译并拷贝 `.so` 的 C++ 扩展，且包顶层 import
sympy——是"开箱即用"的最大障碍。事前审计确认耦合面极小：全项目仅
`gbz_types.py`（今 `bfgbz2d/core.py`）一处 import，实际用到 C++ API
仅 6 个。

### 决议与执行

- `src/bfgbz2d/backend.py`：`LaurentProtocol` + `PolyToolsLaurent`
  （延迟 import）+ `NumpyLaurent`（纯 numpy 回退）+ `make_laurent`
  工厂，选择顺序：显式参数 > `POLY_BACKEND` 环境变量 > poly_tools
  > numpy 回退并一次性警告（环境变量 2026-08-26 由 GBZ_BACKEND
  改名，名实相符）。
- `CharPoly(coeffs, degs, backend=None)`；`solve_roots_1d` 的 0/∞
  padding 逻辑原样保留，仅容器提取改走协议方法
  `partial_terms_1d`。下游 SGBZ / amoeba / continuation 零改动。
- 原型对拍数据（HN2D-10/11 + 随机 24 项 Laurent）：eval/偏导/二阶
  偏导 ~1e-14 一致，根 ~1e-13 一致；`solve_roots_1d` 热路径仅慢
  ~25%（np.roots/LAPACK 主导）。

### 必须传承的三条语义细节（回退后端实现铁律）

1. `denom_orders[d] = max(0, -min_deg_d)` —— Laurent 分母只清负幂；
2. **同次项必须合并**——`solve_roots_1d` 系数装配按唯一次数逐项赋值，
   HN2D 三个 β₂⁰ 次项若不合并会被静默覆盖（C++ 链表会合并，纯
   numpy 实现用 `np.unique + np.add.at` 对齐）；
3. 系数恰为 0 的项必须丢弃——否则破坏"多项式完全消失"空分支
   （空 `degs_list` → deg_M=0 → 全 padding 根）。

### 验收结论

- 双后端对拍测试进 `tests/test_backend.py`（21 项）；
- `POLY_BACKEND=numpy pytest --run-slow` 全绿；
- README 依赖说明完成（numpy+scipy 即可运行，poly_tools 可选加速）。
  **全部验收项通过，原 TODO 文件删除。**

### 衍生问题（另立新 TODO，未关闭）

numpy 后端全套测试暴露"精确简并落在判定边界"的噪声敏感行为
（对称模型 seam 逐位相等触发 touch 通道；非通有二重根重启行 MR
重复记录）。两个后端最终 GBZ 输出逐位一致，分歧仅在内部簿记；
算法层修复方向见 `TODO/exact-degeneracy-boundary-behavior.md`。

---

## 议题二：依赖不可复现（原 TODO/TODO-2026-08-26.md 工程卫生第 1 项，已解决）

### 决议与执行

pyproject.toml + src-layout 打包为 bfGBZ2d（导入名 `bfgbz2d`），
运行时依赖仅 numpy+scipy；无兼容 shim，49 个消费方文件 import 重写，
conftest 改"已装包优先、src 回退"。

### 验收结论

- `pip install .` wheel 构建成功；全套测试绿。
- 原条目从 TODO 文件移除（CI/lint 与 amoeba 两条仍开放）。
- 附带决议：`data/*.pkl` 旧 pickle 因引用 `gbz_types` 旧模块路径失效；
  经实测 11 个文件可用路径重定向 Unpickler 全部无损找回，用户决定
  不迁移（历史数据用打包前 git tag 环境读取）。

---

## 同日关联工作（详见各自日志，未列入 TODO 归档）

- 四阶段打包工程化全记录：`log/2026-08-26-bfgbz2d-packaging.md`
- 常数管理重构（关联性分析 → 归属地原则 → dict 通道处死）：
  `log/2026-08-26-constants-locality.md`；常数总参考 `doc/constants.md`

## 遗留清单（TODO/ 现存）

1. `exact-degeneracy-boundary-behavior.md` —— 精确简并判定边界的
   算法层修复（touch/cross 两侧符号判据、MR 近距重检去重）；
2. `TODO-2026-08-26.md` —— CI/lint 配置；amoeba `mu2_mid` 配对
   `w_left` 的 ε 级偏移与 Newton clamp；
3. `warm-start-zeromanager-across-mu1-probes.md` —— ZM 暖启动
   （μ₁ 探针间继承网格，性能优化，未动工）。
