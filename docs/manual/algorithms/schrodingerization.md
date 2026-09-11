# Schrödingerization 非酉演化模拟（Schrödingerization）

> 类别 C2 · 模块 `pyqecclang.algorithms.schrodingerization` · 阶段 V2

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

- `generator`：$G$ 的 `BlockEncoding`（input model 为 BE），$H_1=(G+G^\dagger)/2$、$H_2=(G-G^\dagger)/(2i)$ 由 `HermitianParts.from_operator` 生成。
- `initial`：`StatePreparation`（SP），作用在物理寄存器。
- `SchrodingerPlan`：辅助寄存器位数 `auxiliary_width`（1..63）、周期 `period` 与选中通道 `selected_index`（$0..2^p-1$）。
- `hamiltonian_function`：`(K, t) -> BlockEncoding`，默认 `taylor_hamiltonian`。
- `fourier_momentum(width, period)`：频率对角算子 $P$ 的独立构造入口。

返回 `StateOracle`（`target`/`signal`），模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"schrodingerization_qode"` |
| `auxiliary_grid` | 辅助网格坐标的 JSON 表 |
| `selected_p` / `recovery_scale` | 选中通道的 $p_j$ 与 $e^{p_j}$ |
| `recovery_assumption` | `"selected p in valid warped region; periodic truncation pending"` |

## 实现要点

$P$ 用 `fourier_momentum` 逐位构造：每位一个投影 BE（alpha 为 1），系数 $(2^b)\cdot 2\pi/\text{period}$，最高位取负以实现二补码负频率，再 LCU 求和，避免稠密矩阵。寄存器布局为 `target = 物理 n 位（低）+ 辅助 p 位（高）`，`work` 承载初态制备的工作位。

生成流程：初态作用在 `target[:n]`；辅助位制备归一化的 $e^{-|p_j|}$（`gate_state_prep`）后做 QFT；对 $K=P\otimes H_1-I\otimes H_2$ 调用 `hamiltonian_function`；逆 QFT 后用 `select_subspace` 选辅助寄存器等于 `selected_index` 的通道，选择条件并入 signal。辅助网格按编码顺序 $p_j=(j$ 若 $j<2^{p-1}$ 否则 $j-2^p)\cdot\text{period}/2^p$；默认配置 $p=2$、period 8、index 1 对应网格 $[0,2,-4,-2]$、选中 $p=2$。

适用边界：实现不自动检验所选通道位于恢复区，也不保证周期窗口足够大（需按 $H_1$ 的传播速度、时间与周期边界选择窗口与通道，不能仅以 $p>0$ 为充分条件）。`recovery_scale` 只记录 warp 的一个因子；理想恢复关系下成功块约为 $e^{-p_j}|u(t)\rangle/(rZ\alpha_E)$（$r$ 为初值范数、$Z$ 为离散 warp 范数、$\alpha_E$ 为演化 BE 归一化），完整范数恢复接口尚未提供，这些量需宿主自行保存。入口处先执行 `operator_state_contract("schrodingerization")` 的能力检查。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）。三层证据：

- 结构：`tests/core/test_differential.py:DifferentialStructureTests.test_four_methods_keep_input_oracles` 的 schrodingerization 子测试——抽象 BE/SP 输入下开放槽（`input_A`/`input_b`）保留、`algorithm` 属性正确、RIR 序列化往返。
- 数值：`tests/core/test_differential.py` 定位为结构测试；Fourier 正负号约定、有限网格与恢复区的整体数值验证仍待完成（见已知缺口）。
- 绑定：同类注册与契约检查——入口契约经 `operator_state_contract` 强制，协议级假设（辅助窗口、Fourier 约定、恢复区域需应用层验证）由 `QODEProtocol.contract` 记录（见 [QODE 问题对象与协议](qode-problem.md)）。

## 已知缺口与计划阶段

统一的"解析可解 ODE 族"收敛基准缺失，归入阶段 V2 的收敛性扫描框架；与 `validation-coverage.md` 的 schrodingerization.py 行一致。辅助窗口与恢复区域的数值验证按协议 assumptions 属应用层责任，当前未自动化。

## 相关链接

- 源码：`src/pyqecclang/algorithms/schrodingerization.py`
- API 参考：[Schrödingerization](../../api/algorithms/schrodingerization.rst)
- 相关页：[QODE 问题对象与协议](qode-problem.md) · [LCHS](lchs.md)（另一条线性路线）· [Carleman 线性化](carleman.md)（提升后接本方法）
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
