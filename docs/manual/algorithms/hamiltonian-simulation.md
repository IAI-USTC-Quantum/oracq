# 哈密顿量模拟协议（Hamiltonian Simulation Protocol）

> 类别 C2 · 模块 [`oracq.algorithms.common.hamiltonian`](../../api/algorithms/common/hamiltonian.rst) · 阶段 V2

## 概述

按输入算子实际具备的访问方式选择哈密顿量模拟路线：可分解的算子走 Trotter 乘积公式（单项演化由 [Trotter 乘积公式模拟](trotter.md) 生成），只有块编码的算子走调用方注入的 QSP 内核。数学目标都是逼近 $e^{-iHt}$；两条路线的精度语义不同——Trotter 是乘积公式近似，QSP 是多项式近似——因此返回对象携带 `correctness = "product formula approximation pending"` 一类的显式标注，不冒称已验证。

## 接口与输入模型

```python
hamiltonian_simulation(operator, time, *, method="auto", steps=2, qsp=None)
```

API 入口：{obj}`hamiltonian_simulation <oracq.algorithms.common.hamiltonian.hamiltonian_simulation>`

- `operator`：厄米算子，需满足 {obj}`HermitianProtocol <oracq.algorithms.common.hamiltonian.HermitianProtocol>`（`hermitian` 属性为 `True`）；input model 为 HAM（Pauli/可分解表示）或 BE（块编码）。
- `time`：演化时间，有限实数。
- `method`：`"auto"`（默认；Trotterizable 走 trotter，否则走 qsp）、`"trotter"`、`"qsp"` 之一，其他值拒绝。
- `steps`：Trotter 步数，正整数，默认 2。
- `qsp`：`qsp(BE, time) -> BlockEncoding` 的可调用实现，`method="qsp"` 时必须注入。

Trotter 路径要求 `operator` 满足 {obj}`TrotterizableProtocol <oracq.algorithms.common.hamiltonian.TrotterizableProtocol>`：`trotter_list()` 返回非空 {obj}`TrotterTerm(coefficient, operator) <oracq.algorithms.common.hamiltonian.TrotterTerm>` 序列，各项算子满足 {obj}`EvolvableProtocol <oracq.algorithms.common.hamiltonian.EvolvableProtocol>`（`evolution(time)` 返回 {obj}`Operation <oracq.infrastructure.builder.Operation>`）。返回 {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`（`alpha = 1.0`），模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"trotter_hamiltonian_protocol"` |
| `steps` | 重复步数 |
| `correctness` | `"product formula approximation pending"` |

内置算子类型：{obj}`PauliOperator(word) <oracq.algorithms.common.hamiltonian.PauliOperator>`（演化即 {obj}`trotter_hamsim <oracq.algorithms.common.hamiltonian.trotter_hamsim>` 单项单步）、{obj}`PauliHamiltonian(terms) <oracq.algorithms.common.hamiltonian.PauliHamiltonian>`（同时提供 BE 与 `trotter_list()` 两种访问）、{obj}`EncodedOperator(encoding, hermitian) <oracq.algorithms.common.hamiltonian.EncodedOperator>`（只持有 BE 与厄米性声明，不带分解）。

## 实现要点

QSP 路径是显式注入点而非已实现内核：库内没有通用 QSP-HamSim 内核，缺少 `qsp` 可调用时直接报错，不静默降级；注入实现的返回值须是 `BlockEncoding` 且目标宽度与输入一致。Trotter 路径逐项调用 `term.operator.evolution(term.coefficient * time / steps)`，要求每项演化是无需后选择的酉——除 `target` 外公开寄存器须为零宽、各项同宽；带信号寄存器、需要后选择的一般多项式 BE 不能冒充单项精确演化。协议组装层用 {obj}`Repeat(steps) <oracq.infrastructure.ir.Repeat>` 加模块调用（资源名 `term0`、`term1`……）引用各项演化，不在生成阶段展开。

适用边界：非厄米动力学被拒绝并明确指向 QODE 路线；"非 Hermitian 输入"与"缺 qsp 实现"是两种独立报错，互不混淆。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）。与验证覆盖矩阵 `hamiltonian.py` 行一致：

- 结构：`tests/core/test_algorithm_protocols.py:AlgorithmProtocolTests.test_trotter_protocol_keeps_phase_and_repeat`（Repeat `count == steps`、`alpha == 1`）；`test_trotter_validates_each_term`（`TrotterTerm` / `PauliHamiltonian` 的项级校验）。
- 数值：`test_trotter_protocol_keeps_phase_and_repeat` 对 $H = 0.3I + 0.7X$、$t = 0.4$、`steps = 3` 逐幅度对拍解析值（places = 11）；`test_trotter_only_input_does_not_need_block_encoding` 见证只有 Trotter 访问、没有 BE 的自定义算子可被求解（单 Z 项系数 0.5、$t = 0.2$ 得 $e^{-0.1i}$）；`test_nonhermitian_and_missing_qsp_are_distinct_errors` 钉死两种失败模式分别抛出。
- 绑定：无独立绑定见证（输入要求具体协议实现，不含抽象槽位）。

## 已知缺口与计划阶段

Trotter 阶数误差率扫描缺失（阶段 V2 收敛性扫描框架，validation-plan §5）；QSP 路径只有契约检查与注入点，内核本身不属于本模块。与验证覆盖矩阵 `hamiltonian.py` 行的缺口列一致。

## 相关链接

- 源码：`src/oracq/algorithms/common/hamiltonian.py`
- 教程：[提供自己的 Hamiltonian 分解](../../tutorials/hamiltonian.md)
- 同族页面：[Trotter 乘积公式模拟](trotter.md)、[截断 Taylor 块编码](taylor-block-encoding.md)
- API 参考：[Hamiltonian 演化](../../api/algorithms/common/hamiltonian.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_hamiltonian.py`（真实后端执行，无模拟替身；27 案例全 PASS）。实验设计覆盖协议两条路线：Trotter 路由用 `PauliHamiltonian`（$H = 0.9\,ZZ + 0.6\,XI + 0.45\,IX$）在 steps ∈ {2, 16} 下经 OriginIR-ext + UniQC `to_matrix` 提取全幺正，与乘积公式经典矩阵逐元素对拍，并以 `scipy.linalg.expm` 为独立 oracle 报告 Trotter 间隙随步数的下降；QSP 注入路径用 `EncodedOperator`（$H = 0.6\,Z + 0.4\,X$，$\alpha = 1.0$，`method="auto"` 自动路由）注入库内 QSVT 内核 {obj}`qsvt_hamiltonian_simulation <oracq.algorithms.common.qsvt.qsvt_hamiltonian_simulation>`（取 $\tau = -t\cdot\alpha$ 实现 $e^{-iHt}$，`error = 0.02`），在 reference / rir-pysparq / originir-ext 三条路径上读出零信号块并与 $e^{-itH}/\text{sim\_scale}$ 对拍。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `hamsim-protocol-trotter-2q` | 3 项，steps ∈ {2, 16} | originir-ext + UniQC to_matrix | 与乘积公式 max_error | 4.5e-15 |
| 同上 | — | — | expm 间隙 r=2 → r=16 | 0.3475 → 0.0420 |
| `hamsim-protocol-qsp-injection` | 1 量子位，$t = 0.5$ | reference、rir-pysparq、originir-ext | 零信号块 max_error | 1.07e-4 |
| 同上 | — | — | 成功概率（signal = 0） | 0.1113 |

Trotter 收敛阶的系统扫描（拟合阶 1.008 / 1.013，与一阶理论一致）见 [Trotter 乘积公式模拟](trotter.md) 的数值验证节。

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_hamiltonian.py
```

产物：`out/verification/hamiltonian.json`。
