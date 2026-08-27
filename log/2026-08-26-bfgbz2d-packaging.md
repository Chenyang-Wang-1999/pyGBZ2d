# bfGBZ2d 工程化：打包、可插拔 backend、常量集中

日期: 2026-08-26

## 背景

项目要"开箱即用"：用户不应被强制安装需手工编译的 poly_tools C++ 扩展，
数值常量散落在 11 个模块（约 45 个，含两处重复定义），仓库无打包元数据。
方案经四阶段审批执行（记录于会话评审；TODO/charpoly-pluggable-backend.md
保存了 backend 的原型验证数据）。

## 变更总览

### Phase 1 — 打包（70d993d）

- `pyproject.toml`（分发名 bfGBZ2d，导入名 `bfgbz2d`，src-layout，
  依赖仅 numpy+scipy）；运行时代码 git mv 入 `src/bfgbz2d/`：
  gbz_types→core、brute_force_SGBZ→sgbz、brute_force_amoeba→amoeba、
  continuation→continuation。
- 无兼容 shim（决策）：49 个消费方文件 import 全部重写；`data/*.pkl`
  旧 pickle 引用 `gbz_types` 模块路径，已标记失效（决策：不迁移）。
- conftest 改为"已装包优先、src 回退"，克隆后免安装可跑测试。

### Phase 2 — 可插拔 backend（834755b）

- `bfgbz2d/backend.py`：LaurentProtocol + PolyToolsLaurent（延迟 import）
  + NumpyLaurent（纯 numpy 回退）+ make_laurent 工厂
  （显式参数 > GBZ_BACKEND 环境变量 > poly_tools > numpy 回退+一次性警告）。
- `CharPoly(coeffs, degs, backend=None)`；solve_roots_1d 的 0/∞ padding
  逻辑不变，仅容器提取改走协议方法 partial_terms_1d。
- 对拍验收：eval/偏导/二阶偏导 ~1e-14 一致，solve_roots_1d 根 ~1e-13
  一致；numpy 后端热路径仅慢 ~25%（np.roots 主导）。
- 语义细节（NumpyLaurent 必须遵守）：denom 只清负幂；**同次项必须合并**
  （HN2D 三个 β₂⁰ 次项若不合并会被系数装配循环静默覆盖）；零系数项丢弃。

### Phase 3a — 常量集中（1e6b0a5）

- `bfgbz2d/config.py`：45 个常量唯一定义点，按调参安全性分三组
  （模型尺度 / 步长与预算 / 机器精度锚定），override() 上下文管理器。
- 全部使用点切 `config.X` 活值读取；重复定义合并（CONTINUUM_TOL、
  CONTINUUM_FRAC 双份合一）；用户决策：CONTINUUM_PERTURB 合并为单一
  1e-4（原 SGBZ 1e-2 / amoeba 1e-4），全套测试（含 --run-slow）验证。

### Phase 3b — 活值默认（d91b958）

- `config.live_defaults` 装饰器：公共入口 kwargs 改 None-sentinel，
  调用时从 config 解析（per-call > 全局 config）；StepControl 字段构造时
  解析；三处 import 时求值的 StepControl() 可变默认实例改 None。
- 优先级契约、装饰器结构不变量进 `tests/test_config.py`（AST 审计防
  装饰器漂移——评审中真实抓到一例错挂，见下）。

## 过程中的重要发现

1. **精确简并的噪声敏感行为**（详见
   TODO/exact-degeneracy-boundary-behavior.md）：对称模型 seam 处
   ln|β₂| 严格相等（对称性，非噪声），touch 通道按逐位零触发；非通有
   二重根重启行的 MR 重复记录。两个后端**最终 GBZ 输出逐位一致**，
   分歧仅在内部簿记（kind 标签、MR 计数）。算法层修复方向已记录，
   两处测试放宽为后端不变断言。
2. **装饰器错挂**：批量改造脚本用 `find("):")` 截取签名，遇 `) -> bool:`
   （带空格）越界，把 live_defaults 挂到无此参数的函数上。AST 审计
   （已固化为测试）抓到两处，已修复。
3. C++ poly_tools 的算术与 Python/cmath 在末位 ULP 上不同（如
   3.9999999999979994 vs 3.999999999998），且 C++ 侧更"脏"——逐位复现
   不可行，也不必要。

## 验证

- 默认后端：`pytest --run-slow` 265 passed；
- 纯 numpy 后端：`GBZ_BACKEND=numpy pytest --run-slow` 265 passed；
- `pip install .` 构建成功（沙箱内仅安装落盘被拦，真实环境无碍）。

## 遗留

- TODO/exact-degeneracy-boundary-behavior.md（算法层 touch/cross 与
  MR 重启去重修复，含恢复严格断言的验收标准）；
- data/ 旧 pickle 失效（用打包前 git tag 环境读取）；
- CI / lint 未配置（原 TODO-2026-08-26 工程卫生项保留）。
