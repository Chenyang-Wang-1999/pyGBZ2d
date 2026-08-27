# 2D Hatano-Nelson 模型解析解

## 模型定义

2D HN 模型的 Hamiltonian：

$$h(\beta_1, \beta_2) = J_{11} \beta_1^{-1} + J_{12} \beta_1 + J_{21} \beta_2^{-1} + J_{22} \beta_2$$

其中 $\beta_1$, $\beta_2$ 是 $(1,0)$ 和 $(0,1)$ 方向的 Bloch 因子。

参数化（将耦合系数分解为 Hermitian 和非 Hermitian 部分）：

$$J_{\alpha 1} = e^{\gamma_\alpha + i\delta_\alpha} J_\alpha, \quad J_{\alpha 2} = e^{-\gamma_\alpha + i\delta_\alpha} J_\alpha^*, \quad \alpha = 1,2$$

- $J_\alpha \in \mathbb{C}$：Hermitian 部分
- $\gamma_\alpha \in \mathbb{R}$：非 Hermitian 部分
- $\delta_\alpha \in \mathbb{R}$：相位因子

[10]-basis 特征多项式（三个变量：$E$, $\beta_1$, $\beta_2$）：

$$f_{xy}(E, \beta_1, \beta_2) = E - J_{11}\beta_1^{-1} - J_{12}\beta_1 - J_{21}\beta_2^{-1} - J_{22}\beta_2$$

---

## [10]-SGBZ（主方向 β₁）

### SGBZ 条件

$$\mu_1 = \gamma_1$$

### SGBZ 点集

$$\beta_1 = e^{\gamma_1 + i\theta_1}, \quad \beta_2 = e^{\gamma_2 + i\theta_2}, \quad \theta_1, \theta_2 \in [-\pi, \pi]$$

即 $|\beta_1| = e^{\gamma_1}$, $|\beta_2| = e^{\gamma_2}$——两个方向的半径均为常数。

### 能谱

$$E(\theta_1, \theta_2) = 2e^{i\delta_1}\,\text{Re}\!\left(J_1^* e^{i\theta_1}\right) + 2e^{i\delta_2}\,\text{Re}\!\left(J_2^* e^{i\theta_2}\right)$$

在复平面上，能谱是由 $\pm 2|J_1|e^{i\delta_1}$ 和 $\pm 2|J_2|e^{i\delta_2}$ 张成的平行四边形。

**特例**：若 $J_1, J_2$ 为实数且 $\delta_1 = \delta_2 = 0$：
$$E(\theta_1, \theta_2) = 2J_1\cos\theta_1 + 2J_2\cos\theta_2 \in [-2|J_1|-2|J_2|,\; 2|J_1|+2|J_2|] \subset \mathbb{R}$$

### GBZ 拓扑类型

- **离散点**：$E$ 在平行四边形内部，有有限个 $(\theta_1, \theta_2)$ 解
- **连续域**：仅在能谱边界（van Hove 奇点附近），多个解简并

---

## [11]-SGBZ（主方向 β_{[11]}）

### 基底变换

[11]-strip 使用基底 $\mathbf{a}_{[11]} = (1,1)$, $\mathbf{a}_y = (0,1)$，对应的 Bloch 因子为 $\beta_{[11]}$ 和 $\tilde{\beta}_y$。

变量代换关系：
$$\beta_1 = \beta_{[11]} \tilde{\beta}_y^{-1}, \quad \beta_2 = \tilde{\beta}_y$$

即在 [10]-basis 多项式中做代换 $\beta_1^{-1} \to \beta_{[11]}^{-1}\tilde{\beta}_y$, $\beta_1 \to \beta_{[11]}\tilde{\beta}_y^{-1}$，得到 [11]-basis 多项式：

$$f_{[11]}(E, \beta_{[11]}, \tilde{\beta}_y) = E - J_{11}\beta_{[11]}^{-1}\tilde{\beta}_y - J_{12}\beta_{[11]}\tilde{\beta}_y^{-1} - J_{21}\tilde{\beta}_y^{-1} - J_{22}\tilde{\beta}_y$$

### SGBZ 条件

$$\mu_{[11]} = \gamma_1 + \gamma_2$$

### SGBZ 点集

$$\tilde{\beta}_{[11]} = e^{\gamma_1+\gamma_2 + i\theta_{[11]}}, \quad \theta_{[11]} \in [-\pi, \pi]$$

$$\tilde{\beta}_y = e^{\gamma_2 + i\theta_y} \sqrt{\left|\frac{J_1^* e^{i\Delta_{12} + i\theta_{[11]}} + J_2}{J_1 e^{i\Delta_{12} - i\theta_{[11]}} + J_2^*}\right|}, \quad \theta_y \in [-\pi, \pi]$$

其中 $\Delta_{12} \equiv \delta_1 - \delta_2$。

即 $|\tilde{\beta}_{[11]}| = e^{\gamma_1+\gamma_2}$ 为常数，而 $|\tilde{\beta}_y|$ 随 $\theta_{[11]}$ 变化。

### 能谱

定义辅助函数：

$$v_{[11]}(\theta_{[11]}; \mu_{[11]}) = J_{12}J_{22}\, e^{\mu_{[11]}+i\theta_{[11]}} + J_{11}J_{21}\, e^{-\mu_{[11]}-i\theta_{[11]}} + J_{11}J_{12} + J_{21}J_{22}$$

以及相位：
$$\varphi_1(\theta_{[11]}) = \text{Arg}\!\left(J_{11}e^{-\mu_{[11]}-i\theta_{[11]}} + J_{22}\right)$$
$$\varphi_2(\theta_{[11]}) = \text{Arg}\!\left(J_{12}e^{\mu_{[11]}+i\theta_{[11]}} + J_{21}\right)$$
$$\bar{\varphi} = \frac{\varphi_1 + \varphi_2}{2}, \quad \Delta\varphi = \varphi_2 - \varphi_1$$

则能谱为：
$$E(\theta_{[11]}, \theta_y) = 2e^{i\bar{\varphi}(\theta_{[11]})} \sqrt{|v_{[11]}(\theta_{[11]}; \gamma_1+\gamma_2)|}\; \cos\!\left(\theta_y - \frac{\Delta\varphi(\theta_{[11]})}{2}\right)$$

**特例**：若 $J_1=J_2=1$, $\delta_1=\delta_2=0$（参数集 A），则：
$$v_{[11]} = e^{i\theta_{[11]}} + e^{-i\theta_{[11]}} + 2 = 4\cos^2(\theta_{[11]}/2)$$
$$\bar{\varphi} = 0, \quad \Delta\varphi = 0$$
$$E(\theta_{[11]}, \theta_y) = 4|\cos(\theta_{[11]}/2)| \cos\theta_y \in [-4, 4]$$

且 $|\tilde{\beta}_y| = \sqrt{|e^{i\theta_{[11]}}+1| / |e^{-i\theta_{[11]}}+1|} = 1$，即 $|\tilde{\beta}_y| = 1 = e^0$（常数）。

---

## Amoeba GBZ

Amoeba GBZ = [10]-SGBZ（或 [01]-SGBZ）。参数：

$$\mu_1 = \gamma_1, \quad \mu_2 = \gamma_2$$

GBZ 点集与 [10]-SGBZ 一致：

$$\beta_1 = e^{\gamma_1 + i\theta_1}, \quad \beta_2 = e^{\gamma_2 + i\theta_2}$$

---

## 解析解验证规则汇总

| GBZ 类型 | μ 条件 | \|β₁\| 条件 | \|β₂\| 条件 |
|---|---|---|---|---|
| [10]-SGBZ | μ₁ = γ₁ | e^{γ₁} | e^{γ₂} |
| [11]-SGBZ | μ₁ = γ₁+γ₂ | e^{γ₁+γ₂} | 变化（由解析公式给出） |
| Amoeba GBZ | μ₁ = γ₁, μ₂ = γ₂ | e^{γ₁} | e^{γ₂} |

---

## 测试参数

### 参数集 A（实数，简单）

```python
J1 = 1.0; J2 = 1.0
gamma_1 = 0.2; gamma_2 = 0.3
delta_1 = 0.0; delta_2 = 0.0
```

- [10]-SGBZ：$\mu_1 = 0.2$, $|\beta_1| = e^{0.2}$, $|\beta_2| = e^{0.3}$
- [11]-SGBZ：$\mu_{[11]} = 0.5$, $|\tilde{\beta}_{[11]}| = e^{0.5}$, $|\tilde{\beta}_y| = 1$
- 能谱：$[-4, 4]$（纯实数，三种 SGBZ 相同）

### 参数集 B（复数，一般情况）

```python
J1 = 1.0 + 0.5j; J2 = 0.8 - 0.3j
gamma_1 = 0.2; gamma_2 = 0.3
delta_1 = 0.1; delta_2 = -0.05
```

### 测试 E_ref 点

| E_ref | [10]-GBZ? | [11]-GBZ? | 说明 |
|---|---|---|---|
| $1.0 + 0i$ | ✅ | ✅ | 参数 A 内部 |
| $0 + 0i$ | ✅ | ✅ | 参数 A 中心 |
| $5.0 + 0i$ | ❌ | ❌ | 参数 A 外部 |
| $2.0 + 0.5i$ | ✅ | ✅ | 参数 B 内部 |

---

## 构造特征多项式

```python
import numpy as np
import poly_tools as pt
from cmath import exp

def build_HN2D_polynomial(J1, J2, gamma_1, gamma_2, delta_1, delta_2, basis="10"):
    """构造 2D HN 模型的特征多项式。

    basis="10": f(E, β₁, β₂) = E - J₁₁β₁^{-1} - J₁₂β₁ - J₂₁β₂^{-1} - J₂₂β₂
    basis="11": f(E, β_{[11]}, β̃_y) = E - J₁₁β_{[11]}^{-1}β̃_y - J₁₂β_{[11]}β̃_y^{-1} - J₂₁β̃_y^{-1} - J₂₂β̃_y
    """
    J11 = exp(gamma_1 + 1j*delta_1) * J1
    J12 = exp(-gamma_1 + 1j*delta_1) * np.conj(J1)
    J21 = exp(gamma_2 + 1j*delta_2) * J2
    J22 = exp(-gamma_2 + 1j*delta_2) * np.conj(J2)

    coeffs = np.array([1, -J11, -J12, -J21, -J22], dtype=complex)

    if basis == "10":
        # f(E, β₁, β₂)
        degs = np.array([
            [1, 0, 0],    # E
            [0, -1, 0],   # -J₁₁ β₁^{-1}
            [0, 1, 0],    # -J₁₂ β₁
            [0, 0, -1],   # -J₂₁ β₂^{-1}
            [0, 0, 1],    # -J₂₂ β₂
        ], dtype=int)
    elif basis == "11":
        # f(E, β_{[11]}, β̃_y)
        degs = np.array([
            [1, 0, 0],    # E
            [0, -1, 1],   # -J₁₁ β_{[11]}^{-1} β̃_y
            [0, 1, -1],   # -J₁₂ β_{[11]} β̃_y^{-1}
            [0, 0, -1],   # -J₂₁ β̃_y^{-1}
            [0, 0, 1],    # -J₂₂ β̃_y
        ], dtype=int)
    else:
        raise ValueError(f"Unknown basis: {basis}")

    return coeffs, degs

def make_char_poly(coeffs, degs):
    """从 coeffs/degs 构造 pt.CLaurent。"""
    char_poly = pt.CLaurent(3)
    char_poly.set_Laurent_by_terms(
        pt.CScalarVec(coeffs),
        pt.CLaurentIndexVec(degs.flatten())
    )
    return char_poly
```

## 验证函数

```python
def verify_point_subset_sgbz10(subset, E_ref, gamma_1, gamma_2, tol=1e-6):
    """验证 [10]-SGBZ 的 PointSubset 满足解析解条件。"""
    assert subset.E == E_ref
    assert abs(abs(subset.beta1) - exp(gamma_1)) < tol
    assert abs(abs(subset.beta2) - exp(gamma_2)) < tol

def verify_point_subset_sgbz11(subset, E_ref, gamma_1, gamma_2, tol=1e-6):
    """验证 [11]-SGBZ 的 PointSubset 满足解析解条件。"""
    assert subset.E == E_ref
    assert abs(abs(subset.beta1) - exp(gamma_1 + gamma_2)) < tol
    # |β₂| 随 θ_{[11]} 变化，不做常数验证

def verify_gbz_result(gbz, E_ref, expect_nonempty=True):
    """验证 GBZResult 结构。"""
    assert gbz.E_ref == E_ref
    if expect_nonempty:
        assert gbz.index != (0, 0)
        assert len(gbz.subsets) > 0
        assert all(isinstance(s, (PointSubset, LineSubset)) for s in gbz.subsets)
        n_0d = sum(1 for s in gbz.subsets if isinstance(s, PointSubset))
        n_1d = sum(1 for s in gbz.subsets if isinstance(s, LineSubset))
        assert gbz.index == (n_0d, n_1d)
    else:
        assert gbz.is_empty
        assert gbz.index == (0, 0)
```
