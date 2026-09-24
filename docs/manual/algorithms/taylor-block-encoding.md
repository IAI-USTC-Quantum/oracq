# 截断 Taylor 块编码（Truncated Taylor Block Encoding）

> 类别 C2 · 模块 [`oracq.algorithms.common.hamiltonian`](../../api/algorithms/common/hamiltonian.rst) · 阶段 V2

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

API 入口：{obj}`taylor_hamiltonian <oracq.algorithms.common.hamiltonian.taylor_hamiltonian>`

- `hamiltonian`：{obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`，演化生成元的块编码（input model 为 BE）；数学语义（如厄米性）由调用方负责。
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

幂序列迭代生成：从 $(1, I)$ 出发，{obj}`product(hamiltonian, current) <oracq.algorithms.input_model.operators.product>` 逐次右乘得到 $H^k$，系数解析取 $(-it)^k / k!$，随后整体经 {obj}`lcu(powers) <oracq.algorithms.input_model.block_encoding.lcu>` 组合。寄存器布局由 LCU 决定：`target = Bits(width)`，`signal` 为 selector（$\lceil\log_2(d+1)\rceil$ 位）与各项信号空间的拼接；复数系数的相位经全局相位门实现。

适用边界：门数随 $d$ 线性放大（深幂链），本函数不做截断误差界分析，误差由调用方按 $\lVert H\rVert t$ 与 $d$ 自行评估；生成期校验 `degree` 为非负整数、`time` 为有限实数。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）。与验证覆盖矩阵一致的说法是：`hamiltonian.py` 行未给本入口单列数值见证，本函数的证据来自被 ODE 求解器消费——

- 结构：`tests/core/test_differential.py:DifferentialStructureTests.test_four_methods_keep_input_oracles` 以 `degree=1` 注入 lchs / cbmd / schrodingerization / carleman 四条路线，检查输入 oracle 槽位保留与 IR JSON 往返（登记于验证覆盖矩阵 ode 行）。
- 数值：`tests/core/test_sde.py:SolverContractTests.test_qode_problem_accepted_by_lchs` 同样以 `degree=1` 注入 LCHS 契约并求解（登记于 sde 行）。截断精度本身的"误差 vs 阶数"直接对拍缺失。
- 绑定：无独立绑定见证；LCU 组合子的绑定语义由 `tests/core/test_language.py:BlockEncodingTests.test_alpha_survives_ir_serialization` 覆盖（见[块编码组合代数](block-encoding-algebra.md)）。

## 已知缺口与计划阶段

与验证覆盖矩阵 `hamiltonian.py` 行一致：模块登记缺口为 Trotter 阶数误差率扫描（阶段 V2 收敛性扫描框架）；本入口自身的 Taylor 截断误差对拍属于同类收敛性证据，待同一框架接入后补齐。

## 相关链接

- 源码：`src/oracq/algorithms/common/hamiltonian.py`
- 同族页面：[哈密顿量模拟协议](hamiltonian-simulation.md)、[Trotter 乘积公式模拟](trotter.md)、[块编码组合代数](block-encoding-algebra.md)
- API 参考：[Hamiltonian 演化](../../api/algorithms/common/hamiltonian.rst)
- 概念：[Oracle 与算子表示](../operators.md)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值验证脚本：`tests/verification/verify_blockencoding.py`（真实后端执行，无 mock、无 skip；2026-09-16 共 73 个案例全部通过），产物 `out/verification/blockencoding.json`。本页对应 `taylor-block-encoding-d{1,2,3,4}` 共 4 个案例，补上"截断精度本身的误差 vs 阶数直接对拍缺失"的登记缺口。

实验设计：输入生成元为 2×2 厄米矩阵 $H=\begin{pmatrix}1.0&0.4\\0.4&-0.6\end{pmatrix}$（经 {obj}`matrix_pauli_encoding <oracq.algorithms.input_model.block_encoding.matrix_pauli_encoding>` 编码，$\alpha_{\rm in}=1.4$），演化时间 $t=0.7$，阶数 $d=1..4$。块语义为 $\big(\sum_{k\le d}(-itH)^k/k!\big)/\alpha$，$\alpha=\sum_{k\le d}(\alpha_{\rm in}t)^k/k!$。两类误差分离报告：**实现误差**——reference 逐列提取的零信号块对照同阶截断级数（numpy 独立计算）；**方法误差**——块对照 $e^{-iHt}/\alpha$（numpy 特征分解独立计算）；另对照 `be_alpha` 与级数闭式，并以 rir-pysparq / adapter-pysparq 逐列交叉。

| 案例 | 阶数 | 实现误差 | 方法误差（信息性） | $\alpha$ 偏差 | 后端交叉 |
|---|---|---|---|---|---|
| `taylor-block-encoding-d1` | 1 | 8.3e-17 | 1.41e-1 | 0 | 0 |
| `taylor-block-encoding-d2` | 2 | 1.8e-16 | 2.81e-2 | 0 | 0 |
| `taylor-block-encoding-d3` | 3 | 1.7e-16 | 5.20e-3 | 0 | 0 |
| `taylor-block-encoding-d4` | 4 | 1.9e-16 | 7.75e-4 | 0 | 0 |

方法误差按阶数单调下降（约 $(\alpha_{\rm in}t)^{d+1}/(d+1)!$ 量级，$\alpha_{\rm in}t=0.98$），实现误差始终处于机器精度——截断误差确实由调用方按 $d$ 控制，组装本身不引入额外误差。

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_blockencoding.py
```

产物：`out/verification/blockencoding.json`。
