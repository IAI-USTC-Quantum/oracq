# Trotter 乘积公式模拟（Trotter Product-Formula Simulation）

> 类别 C2 · 模块 [`oracq.algorithms.common.hamiltonian`](../../api/algorithms/common/hamiltonian.rst) · 阶段 V2

## 概述

给定 Pauli 分解 $H = \sum_j c_j P_j$（各 $P_j$ 为同宽 Pauli 词），用一阶乘积公式模拟 $e^{-iHt}$：

$$
e^{-iHt} \approx \left(\prod_j e^{-i c_j P_j\, t / r}\right)^{r},
$$

重复步数 $r$ 越大，乘积公式误差越小。单项 $e^{-i\theta P}$ 经基变换化为计算基上的相位旋转：对 $P$ 的非恒等位做基变换（X 位加 H，Y 位先加相位门 $-\pi/2$ 再加 H），用 CNOT 链把激活位的奇偶折叠到末位，对末位施加相位旋转后逆序复原。

## 接口与输入模型

```python
trotter_hamsim(terms, final_time, *, steps=2)
```

API 入口：{obj}`trotter_hamsim <oracq.algorithms.common.hamiltonian.trotter_hamsim>`

- `terms`：`((coefficient, word), ...)` 序列，`coefficient` 为实数（内部经 `float` 转换），`word` 为同宽 I/X/Y/Z 字符串；input model 为 HAM（Pauli 项分解）。
- `final_time`：总演化时间 $t$（有限实数）。
- `steps`：重复步数 $r$，正整数，默认 2。

返回 {obj}`Operation <oracq.infrastructure.builder.Operation>`（裸操作而非 {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>`），寄存器为 `target = Bits(width)`、`signal = Bits(0)`——演化是无需后选择的确切酉。模块属性：

| 属性 | 含义 |
|---|---|
| `algorithm` | `"trotter_hamsim"` |
| `validation_stage` | `"paradigm"` |

生成期校验：`terms` 为空或 `steps < 1`、各 Pauli 词宽度不一致时抛出 {obj}`ValidationError <oracq.infrastructure.ir.ValidationError>`。全 I 的项不消耗量子比特，退化为全局相位 $e^{-i c_j t / r}$。

## 实现要点

电路主体是单个 {obj}`Repeat(steps) <oracq.infrastructure.ir.Repeat>` 块，块内按输入顺序逐项演化；不在生成阶段展开，文本序列化与后端降低保留 Repeat 结构。每项按"基变换 → CNOT 链 → 末位相位旋转 → 复原"生成：CNOT 链只覆盖非 I 位，末位旋转角为 $2 c_j t / r$（`rz(θ)` 按 $e^{-i\theta Z/2}$ 约定，乘 2 后恰好实现 $e^{-i\theta Z^{\otimes k}}$ 型的相位）。

适用边界：只接受 Pauli 词分解，不做高阶 Suzuki 编排（更高阶公式需调用方自行重排 `terms` 并配合更大 `steps`）；系数经 `float()` 转换，复系数直接报错。协议层入口见[哈密顿量模拟协议](hamiltonian-simulation.md)，它把任意可 Trotter 化算子的单项演化路由回本函数（`PauliOperator.evolution` 即以 `steps=1` 调用本函数）。

## 验证方案

类别 C2（近似连续语义，判定准则见 `../development/validation-plan.md` §2）。与验证覆盖矩阵 `hamiltonian.py` 行一致，三层证据：

- 结构：`tests/core/test_algorithm_protocols.py:AlgorithmProtocolTests.test_trotter_protocol_keeps_phase_and_repeat` 钉死 Repeat 结构（模块 body 首块 `count == steps`）；`test_trotter_validates_each_term` 覆盖项级校验（{obj}`TrotterTerm <oracq.algorithms.common.hamiltonian.TrotterTerm>` 拒绝非 Evolvable 算子，宽度不一致的 {obj}`PauliHamiltonian <oracq.algorithms.common.hamiltonian.PauliHamiltonian>` 拒绝构造）。
- 数值：`test_trotter_protocol_keeps_phase_and_repeat` 对 $H = 0.3I + 0.7X$、$t = 0.4$、`steps = 3` 逐幅度对拍 $e^{-0.12i}\bigl(\cos 0.28\,|0\rangle - i\sin 0.28\,|1\rangle\bigr)$（places = 11），间接覆盖本函数生成的逐项阶段顺序；`test_trotter_only_input_does_not_need_block_encoding` 见证单项 Z 演化（系数 0.5、$t = 0.2$）输出 $e^{-0.1i}$。
- 绑定：无独立绑定见证（输入是显式 Pauli 分解，不含抽象槽位）。

## 已知缺口与计划阶段

Trotter 阶数误差率扫描缺失：误差随 `steps` / $\Delta t$ 变化的批量曲线未自动化，待阶段 V2 的收敛性扫描框架（validation-plan §5）落地后接入；当前数值见证只在单一参数点上对拍解析值。

## 相关链接

- 源码：`src/oracq/algorithms/common/hamiltonian.py`
- 教程：[提供自己的 Hamiltonian 分解](../../tutorials/hamiltonian.md)
- 同族页面：[哈密顿量模拟协议](hamiltonian-simulation.md)、[截断 Taylor 块编码](taylor-block-encoding.md)
- API 参考：[Hamiltonian 演化](../../api/algorithms/common/hamiltonian.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值实验见 `tests/verification/verify_hamiltonian.py`（真实后端执行，无模拟替身；27 案例全 PASS，总耗时约 28 s）。实验设计：单项 Pauli 演化（10 个实例，含 Y 基变换、2–4 量子位 CNOT 链、恒等项全局相位、负时间）经 OriginIR-ext 导出并由 UniQC `Circuit.to_matrix` 提取全幺正，与独立闭式 $\cos(ct)I - i\sin(ct)P$ 逐元素对拍；多项非对易分解（1–4 量子位，steps ∈ {1,2,3,5}，t ∈ {0.4, −0.7, 0.9, 1.3}）与乘积公式经典矩阵 $\bigl(\prod_j e^{-ic_jP_jt/r}\bigr)^r$ 对拍；3 量子位 TFIM 在均匀叠加输入下做四路径交叉对拍；收敛阶用 2 与 3 量子位 TFIM 在 $t = 1.0$ 下扫 steps ∈ {1,…,64}，以 `scipy.linalg.expm` 为独立 oracle 取谱范数误差并做 log–log 拟合。UniQC `to_matrix` 会裁掉尾部无门量子位（Pauli 词尾随 I 的情形），验证脚本按恒等因子嵌回后对拍（见脚本 `embed_unitary` 注释）。

| 案例 | 规模 | 后端路径 | 指标 | 数值 |
|---|---|---|---|---|
| `trotter-single-term-exact` | 10 实例，1–4 量子位 | originir-ext + UniQC to_matrix | max_error | 4.4e-16 |
| `trotter-product-formula-1q` | 3 项（含恒等项），steps×t 网格 | originir-ext + UniQC to_matrix | max_error | 5.7e-16 |
| `trotter-product-formula-2q` | 4 项，steps ∈ {1,2,5} | originir-ext + UniQC to_matrix | max_error | 1.5e-15 |
| `trotter-product-formula-3q` | 5 项 TFIM，steps ∈ {2,3} | originir-ext + UniQC to_matrix | max_error | 1.4e-15 |
| `trotter-product-formula-4q` | 5 项，steps=3 | originir-ext + UniQC to_matrix | max_error | 1.7e-15 |
| `trotter-superposition-cross-3q` | 5 项 TFIM，steps=3，叠加穷举 | reference、rir-pysparq、adapter-pysparq、originir-ext | max_deviation | 1.2e-15 |
| `trotter-convergence-2q` | $t = 1.0$，$r$ = 1…64 | originir-ext + UniQC to_matrix | 拟合收敛阶（理论 1） | 1.008（误差 0.7627 → 0.01049） |
| `trotter-convergence-3q` | $t = 1.0$，$r$ = 1…32 | originir-ext + UniQC to_matrix | 拟合收敛阶（理论 1） | 1.013（误差 1.0680 → 0.02809） |

乘积公式语义验证到机器精度（≤ 1.7e-15），一阶 Lie–Trotter 的收敛阶拟合值 1.008 / 1.013 与理论值 1 一致，填补了本节原先"Trotter 阶数误差率扫描缺失"的缺口。

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_hamiltonian.py
```

产物：`out/verification/hamiltonian.json`（含全部步数的误差序列与逐词误差）。
