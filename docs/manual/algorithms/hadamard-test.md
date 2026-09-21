# Hadamard 检验（Hadamard Test）

> 类别 C3 · 模块 `pyqecclang.algorithms.common.estimation` · 阶段 V2

## 概述

测量复期望值 $\langle\psi|U|\psi\rangle$ 的实部或虚部：单探针比特经受控 $U$ 干涉后，`probe` 的 Z 期望等于所选分量。该操作不执行测量；采样与概率差的经典计算由宿主完成。

## 接口与输入模型

```python
hadamard_test(unitary, preparation=None, *, component="real")
```

- `unitary`：完整酉 `Operation`，或无信号位且 `alpha=1` 的块编码（UO，经 `as_block_encoding` 适配）；带信号位或 $\alpha\ne 1$ 的输入在生成期抛 `ValidationError`。
- `preparation`：初态制备（SP，零输入 + 复净工作区）；省略时用目标空间的零基态（`basis_state(width)`）。
- `component`：`"real"` 或 `"imag"`。

返回 `Operation`，寄存器为 `target`、`work`、`probe`。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"hadamard_test"` |
| `readout_register` | `"probe"` |
| `component` | `"real"` / `"imag"` |

制备宽度必须等于酉宽度，`component` 取值非法时生成期报错。

## 实现要点

电路为制备 $\to$ H(probe) $\to$ 受控 $U$ $\to$（虚部分量时对 probe 施加 $-\pi/2$ 相位）$\to$ H(probe)。虚部读出把 $S^\dagger$ 型相位放在受控调用之后，将正交分量旋进 Z 基。受控调用以零宽信号视图（`b["work"][:0]`）适配 BE 接口，因此 BE 形式的输入不引入信号位。

## 验证方案

类别 C3（判定准则见 `docs/development/validation-plan.md` §2：输出分布等于闭式期望）。三层证据与 `docs/development/validation-coverage.md` 的 `estimation.py` 行一致：

- 结构：`tests/core/test_algorithm_expansion.py:AlgorithmExpansionTests` 的构造与属性断言。
- 数值：`AlgorithmExpansionTests.test_hadamard_test_real_and_imaginary`——单比特相位门（角度 0.6，本征值 $e^{0.6i}$）作用在 $|1\rangle$ 上，`real`/`imag` 两个分量分别对拍 $p(0)-p(1)=\cos 0.6$ 与 $\sin 0.6$，places=10。
- 绑定：模块行登记为"—"。

见证技术为精确态矢量模拟取 probe 边际分布后对拍闭式，无采样断言。

## 已知缺口与计划阶段

本算法无单独缺口；`estimation.py` 登记的唯一缺口属于 QAE 的置信区间声明（见[振幅估计](qae.md)），模块整体处于阶段 V2。

## 相关链接

- 源码：`src/pyqecclang/algorithms/estimation.py`
- API 参考：[相位、振幅与重叠估计](../../api/algorithms/common/estimation.rst)
- 同组页面：[SWAP 检验](swap-test.md)（两态重叠的对应读出）
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

`tests/verification/verify_estimation.py` 在真实后端上对本接口做期望值级数值验证（共 6 个案例，全部通过）。实验设计：

- 单比特相位门 $U=\mathrm{diag}(1,e^{i\theta})$ 作用于 $|1\rangle$（$\theta=0.6$ 与 $-1.1$）：闭式期望 $e^{i\theta}$。
- 双比特对角幺正 $U|x\rangle=e^{i\theta_x}|x\rangle$（角度表 $(0.35,-0.9,1.7,0.55)$，受控全局相位实现）配 `gate_state_prep` 复幅度制备 $\psi=(0.5,\,0.5i,\,0.5,\,-0.5)$：期望值 $\sum_x|\psi_x|^2 e^{i\theta_x}$ 由 math 库独立计算。
- 后端路径：reference、rir-pysparq、adapter-pysparq、originir-ext 四条；probe 的 Z 期望 $p(0)-p(1)$ 对照期望值实/虚部分量。

| 案例 | 规模 | 路径 | 指标值 |
|---|---|---|---|
| 相位门 $\theta=0.6$ | target 1 位 | 全部四条 | real $=0.8253356$（$\cos\theta$），imag $=0.5646425$（$\sin\theta$），误差 $\le 5.6\times10^{-16}$ |
| 相位门 $\theta=-1.1$ | target 1 位 | 全部四条 | real $=0.4535961$，imag $=-0.8912074$，误差 $\le 4.4\times10^{-16}$ |
| 对角幺正 + 复制备 | target 2 位 | 全部四条 | real $=0.5711657$，imag $=0.2684807$，误差 $\le 1.1\times10^{-16}$ |

复现命令：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_estimation.py
```

产物路径：`out/verification/estimation.json`（案例名前缀 `hadamard-`）。
