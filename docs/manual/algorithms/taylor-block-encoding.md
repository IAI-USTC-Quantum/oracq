# 截断 Taylor 块编码（Truncated Taylor Block Encoding）

> 类别 C2 · 模块 `pyqecclang.algorithms.hamiltonian` · 阶段 V2

## 概述

把矩阵指数的截断 Taylor 级数组装成块编码：给定生成元的块编码与阶数 $d$，

$$
e^{-iHt} \;\approx\; \sum_{k=0}^{d} \frac{(-it)^k}{k!}\, H^k,
$$

对 $k$ 次幂 $H^k$（输入 BE 的 $k$ 重乘积）与解析系数做线性组合。模块 docstring 将其定位为"可闭合的普通 Hamiltonian-function BE"，供 QODE 求解器作为可注入的 `hamiltonian_function` 使用；它不要求 Hermitian，也不宣称是高效最优的 HamSim 算法，可替换为 QSP/HamSim 协议（见[哈密顿量模拟协议](hamiltonian-simulation.md)）。

## 接口与输入模型

```python
taylor_hamiltonian(hamiltonian, time, *, degree=2)
```

- `hamiltonian`：`BlockEncoding`，演化生成元的块编码（input model 为 BE）；数学语义（如厄米性）由调用方负责。
- `time`：演化时间，有限实数。
- `degree`：截断阶数 $d$，非负整数，默认 2。

返回 `BlockEncoding`，零信号块约等于 $\sum_k (-it)^k H^k / k!$。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"truncated_taylor_hamiltonian_function"` |
| `degree` / `time` | 调用参数回显 |
| `be_alpha` | LCU 归一化 $\sum_{k=0}^{d} (\alpha t)^k / k!$（$\alpha$ 为输入 BE 归一化） |
| `correctness` / `success_condition` | `"pending"` / `"signal == 0"` |

## 实现要点

幂序列迭代生成：从 $(1, I)$ 出发，`product(hamiltonian, current)` 逐次右乘得到 $H^k$，系数解析取 $(-it)^k / k!$，随后整体经 `lcu(powers)` 组合。寄存器布局由 LCU 决定：`target = Bits(width)`，`signal` 为 selector（$\lceil\log_2(d+1)\rceil$ 位）与各项信号空间的拼接；复数系数的相位经全局相位门实现。

适用边界：门数随 $d$ 线性放大（深幂链），本函数不做截断误差界分析，误差由调用方按 $\lVert H\rVert t$ 与 $d$ 自行评估；生成期校验 `degree` 为非负整数、`time` 为有限实数。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）。与验证覆盖矩阵一致的说法是：`hamiltonian.py` 行未给本入口单列数值见证，本函数的证据来自被 ODE 求解器消费——

- 结构：`tests/core/test_differential.py:DifferentialStructureTests.test_four_methods_keep_input_oracles` 以 `degree=1` 注入 lchs / cbmd / schrodingerization / carleman 四条路线，检查输入 oracle 槽位保留与 IR JSON 往返（登记于验证覆盖矩阵 ode 行）。
- 数值：`tests/core/test_sde.py:SolverContractTests.test_qode_problem_accepted_by_lchs` 同样以 `degree=1` 注入 LCHS 契约并求解（登记于 sde 行）。截断精度本身的"误差 vs 阶数"直接对拍缺失。
- 绑定：无独立绑定见证；LCU 组合子的绑定语义由 `tests/core/test_language.py:BlockEncodingTests.test_alpha_survives_ir_serialization` 覆盖（见[块编码组合代数](block-encoding-algebra.md)）。

## 已知缺口与计划阶段

与验证覆盖矩阵 `hamiltonian.py` 行一致：模块登记缺口为 Trotter 阶数误差率扫描（阶段 V2 收敛性扫描框架）；本入口自身的 Taylor 截断误差对拍属于同类收敛性证据，待同一框架接入后补齐。

## 相关链接

- 源码：`src/pyqecclang/algorithms/hamiltonian.py`
- 同族页面：[哈密顿量模拟协议](hamiltonian-simulation.md)、[Trotter 乘积公式模拟](trotter.md)、[块编码组合代数](block-encoding-algebra.md)
- API 参考：[Hamiltonian 演化](../../api/algorithms/hamiltonian.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
