# Schrödingerization 非酉演化模拟（Schrödingerization）

[English](../../../index.html) · **简体中文**

> 类别 C2 · 模块 [`oracq.algorithms.qode.schrodingerization`](../../api/algorithms/qode/schrodingerization.rst) · 阶段 V2

## 概述

对 $u'=Gu$ 写 $G=H_1+iH_2$（$H_1,H_2$ 均为 Hermitian）。引入辅助变量 $p$，在有效恢复区域使用 $v(t,p)=e^{-p}u(t)$，初值对两侧延拓为 $v(0,p)=e^{-|p|}u_0$，则提升方程为

$$
\partial_t v=-H_1\,\partial_p v+iH_2 v ,
$$

Fourier 化后成为 Hermitian Hamiltonian $K=P\otimes H_1-I\otimes H_2$ 的演化（依据 [Schrödingerization 技术论文 §3.1](https://arxiv.org/html/2212.14703v1)，另见[微分方程章 §5](../differential-equations.md)）。本实现把辅助坐标保留为量子寄存器，最终选择一个 $p$ 通道读出物理解。

## 接口与输入模型

```python
schrodinger_qode(generator, initial, time, *, plan=None,
                 hamiltonian_function=taylor_hamiltonian)
SchrodingerPlan(auxiliary_width=2, period=8.0, selected_index=1)
fourier_momentum(width, period)
```

API 入口：{obj}`schrodinger_qode <oracq.algorithms.qode.schrodingerization.schrodinger_qode>`、{obj}`SchrodingerPlan <oracq.algorithms.qode.schrodingerization.SchrodingerPlan>`、{obj}`fourier_momentum <oracq.algorithms.qode.schrodingerization.fourier_momentum>`

- `generator`：$G$ 的 {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`（input model 为 BE），$H_1=(G+G^\dagger)/2$、$H_2=(G-G^\dagger)/(2i)$ 由 `HermitianParts.from_operator` 生成。
- `initial`：{obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>`（SP），作用在物理寄存器。
- {obj}`SchrodingerPlan <oracq.algorithms.qode.schrodingerization.SchrodingerPlan>`：辅助寄存器位数 `auxiliary_width`（1..63）、周期 `period` 与选中通道 `selected_index`（$0..2^p-1$）。
- `hamiltonian_function`：`(K, t) -> BlockEncoding`，默认 {obj}`taylor_hamiltonian <oracq.algorithms.common.hamiltonian.taylor_hamiltonian>`。
- {obj}`fourier_momentum(width, period) <oracq.algorithms.qode.schrodingerization.fourier_momentum>`：频率对角算子 $P$ 的独立构造入口。

返回 {obj}`StateOracle <oracq.algorithms.input_model.oracles.StateOracle>`（`target`/`signal`），模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"schrodingerization_qode"` |
| `auxiliary_grid` | 辅助网格坐标的 JSON 表 |
| `selected_p` / `recovery_scale` | 选中通道的 $p_j$ 与 $e^{p_j}$ |
| `recovery_assumption` | `"selected p in valid warped region; periodic truncation pending"` |

## 实现要点

$P$ 用 `fourier_momentum` 逐位构造：每位一个投影 BE（alpha 为 1），系数 $(2^b)\cdot 2\pi/\text{period}$，最高位取负以实现二补码负频率，再 LCU 求和，避免稠密矩阵。寄存器布局为 `target = 物理 n 位（低）+ 辅助 p 位（高）`，`work` 承载初态制备的工作位。

生成流程：初态作用在 `target[:n]`；辅助位制备归一化的 $e^{-|p_j|}$（{obj}`gate_state_prep <oracq.algorithms.input_model.oracles.gate_state_prep>`）后做 QFT；对 $K=P\otimes H_1-I\otimes H_2$ 调用 `hamiltonian_function`；逆 QFT 后用 {obj}`select_subspace <oracq.algorithms.common.state_preparation.select_subspace>` 选辅助寄存器等于 `selected_index` 的通道，选择条件并入 signal。辅助网格按编码顺序 $p_j=(j$ 若 $j<2^{p-1}$ 否则 $j-2^p)\cdot\text{period}/2^p$；默认配置 $p=2$、period 8、index 1 对应网格 $[0,2,-4,-2]$、选中 $p=2$。

适用边界：实现不自动检验所选通道位于恢复区，也不保证周期窗口足够大（需按 $H_1$ 的传播速度、时间与周期边界选择窗口与通道，不能仅以 $p>0$ 为充分条件）。`recovery_scale` 只记录 warp 的一个因子；理想恢复关系下成功块约为 $e^{-p_j}|u(t)\rangle/(rZ\alpha_E)$（$r$ 为初值范数、$Z$ 为离散 warp 范数、$\alpha_E$ 为演化 BE 归一化），完整范数恢复接口尚未提供，这些量需宿主自行保存。入口处先执行 {obj}`operator_state_contract("schrodingerization") <oracq.algorithms.input_model.interfaces.operator_state_contract>` 的能力检查。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）。三层证据：

- 结构：`tests/core/test_differential.py:DifferentialStructureTests.test_four_methods_keep_input_oracles` 的 schrodingerization 子测试——抽象 BE/SP 输入下开放槽（`input_A`/`input_b`）保留、`algorithm` 属性正确、RIR 序列化往返。
- 数值：`tests/core/test_differential.py` 定位为结构测试；Fourier 正负号约定、有限网格与恢复区的整体数值验证仍待完成（见已知缺口）。
- 绑定：同类注册与契约检查——入口契约经 `operator_state_contract` 强制，协议级假设（辅助窗口、Fourier 约定、恢复区域需应用层验证）由 `QODEProtocol.contract` 记录（见 [QODE 问题对象与协议](qode-problem.md)）。

## 已知缺口与计划阶段

统一的"解析可解 ODE 族"收敛基准缺失，归入阶段 V2 的收敛性扫描框架；与 `validation-coverage.md` 的 schrodingerization.py 行一致。辅助窗口与恢复区域的数值验证按协议 assumptions 属应用层责任，当前未自动化。

## 数值验证

论文级数值实验见 `tests/verification/verify_ode.py`（ode 组，本文件覆盖 `schrodingerization.py` 的部分）。经典参考全部独立：numpy 全堆叠仿真（warp→正号 DFT→截断 Taylor 或精确 `expm` 演化→逆 DFT→通道选择）、解析解（$e^{-t}$ 衰减、解析旋转矩阵）；量子程序在真实后端 reference 与 OriginIR-ext 上执行（嵌套 LCU 的 $8^{\text{degree}}$ 分支标度使 pysparq 系路径超出预算，案例中注明）。

**实验设计**：(a) 动量块——`fourier_momentum(2, 8.0)` 的零信号块对二补码有符号频率对角矩阵（幺正 `to_matrix` + reference 交叉）；(b) 标量衰减 $G=-I$、plan(auxiliary_width=2, period=1.2, selected_index=1)（$\Delta p=t$ 使网格平移逐点精确）、Taylor degree 4、$t=0.3$：组装保真度、恢复幅值与误差分解；(c) 独立重组装——验证脚本内用公开组合子组装 $K'=-P\otimes H_1-I\otimes H_2$（修复后的库约定）的同构真实量子程序，作为符号约定的回归钉；(d) 旋转生成元 $G=J$（$H_1=0$，纯 $H_2$）对解析旋转的恢复。

**关键指标**：

| 案例 | 规模 | 路径 | 指标 | 数值 |
|---|---|---|---|---|
| fourier-momentum-block | 4 量子位 | originir+to_matrix、reference | 动量块最大偏差 | 2.2e-16 |
| schrodingerization-scalar-decay-grid | 21 量子位 | reference、originir | 组装保真度（对 $K'$ 全堆叠仿真） | 2.0e-18 |
| 同上 | | | 恢复幅值对 $e^{-t}u_0$（正向流）的精确演化网格误差 | 1.3e-16 |
| 同上 | | | 对 $e^{+t}u_0$（时间反演）的失配（方向性证据） | 0.43 |
| 同上 | | | 端到端恢复误差（Nyquist 模 Taylor 余项） | 3.3e-2 |
| schrodingerization-sign-flipped-degree2/4 | 同上 | reference | 独立重组装实现误差 / 网格误差 | 8.8e-18 / 1.3e-16 |
| 同上 | | | 量子恢复误差（Nyquist 模 Taylor 余项，随阶数下降） | 9.1e-2 → 3.3e-2 |
| schrodingerization-rotation-recovery | 同上 | reference | 实现误差 / 恢复误差 | 9.3e-19 / 1.5e-5 |

**符号约定的发现与修复**：历史版本装配 $K=P\otimes H_1-I\otimes H_2$，配合正号 DFT 约定给出的 warp 传输方向与恢复关系相反——恢复幅值精确拟合时间反演解（对 $e^{+t}u_0$ 拟合 4.8e-16、对 $e^{-t}u_0$ 失配 0.43，互为决定性证据）。验证轮据此把库内生成元的动量项符号修复为 $K'=-P\otimes H_1-I\otimes H_2$：同一网格上精确演化的恢复误差降到 1.3e-16，证明符号是唯一结构性失配；库现以 $K'$ 为约定，并由"独立重组装与库程序逐振幅一致"的回归案例钉住（防回退）。纯 $H_2$ 生成元（旋转案例）不受该项符号影响，端到端恢复误差 1.5e-5。端到端残余误差是 Nyquist 动量模（相位 $\pi$）的 Taylor 截断，随阶数下降（degree 2→4：9.1e-2→3.3e-2），属于可替换 `hamiltonian_function` 协议的方法误差。

**复现**：

```bash
PYTHONPATH=src /home/agony/projects/qcfd-dev/quantum-cfd-software/.venv/bin/python tests/verification/verify_ode.py
```

产物：`out/verification/ode.json`（`schrodingerization-*` 与 `fourier-momentum-block` 共 6 个案例）。

## 相关链接

- 源码：`src/oracq/algorithms/qode/schrodingerization.py`
- 教程：[为同一个线性问题替换 QODE 方法](../../tutorials/differential-equations.md)
- API 参考：[Schrödingerization](../../api/algorithms/qode/schrodingerization.rst)
- 相关页：[QODE 问题对象与协议](qode-problem.md) · [LCHS](lchs.md)（另一条线性路线）· [Carleman 线性化](carleman.md)（提升后接本方法）
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
