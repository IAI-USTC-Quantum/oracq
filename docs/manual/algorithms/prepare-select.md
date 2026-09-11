# PREPARE–SELECT 分解（PREPARE–SELECT Decomposition）

> 类别 C5/C1 · 模块 `pyqecclang.algorithms.prepare_select` · 阶段 V4

## 概述

LCU 分解 $H = \sum_i c_i P_i$ 的块编码标准组装（Low & Chuang 2019；Babbush et al. 2018，模块 docstring 引用；input model 词汇中的 BE / HAM）：

$$
\mathrm{BE} = (\mathrm{PREPARE}^\dagger \otimes I)\cdot \mathrm{SELECT}\cdot(\mathrm{PREPARE} \otimes I), \qquad \alpha = \sum_i |c_i| .
$$

PREPARE 在 selector 寄存器上制备振幅 $\propto \sqrt{|c_i|}$ 的态，SELECT 按 selector 值对 target 施加第 $i$ 个 Pauli 字。系数相位按仓库惯例进 SELECT（经 `global_phase` 施加），PREPARE 只编码模长。本页覆盖通用分解与 gate / QRAM / abstract 三层 PREPARE；alias 采样路线单独成页（见 [Alias 采样制备](alias-preparation.md)）。

## 接口与输入模型

```python
abstract_prepare(coefficients, *, work_width=0, name=None)
gate_prepare(coefficients, *, name=None)
qram_prepare(coefficients, *, angle_width=8)
select_pauli(terms)
lcu_prepare_select(terms, *, prepare=None)
```

- `coefficients`：至少两个有限复系数、不全为零（单项请直接用 `scale`）。
- `abstract_prepare`：开放声明，系数写入声明属性，`work_width` 须与后续绑定实现的 work 宽度一致（gate 为 0，QRAM 为 selector 位宽加角度位宽）。
- `qram_prepare`：返回 `QramPreparation(preparation, memory)`，角表 `memory = {"angles": ...}` 由调用方作为 QRAM 数据绑定（`angle_width` 2..32）。
- `select_pauli` / `lcu_prepare_select` 的 `terms`：`(系数, Pauli 字)` 序列或 `PauliHamiltonian`；字仅允许 `IXYZ` 且各项等宽，SELECT 至少两项。
- `lcu_prepare_select`：`prepare` 缺省 `gate_prepare`，也可传 abstract / QRAM / alias 句柄（经 `as_state_preparation` 适配）；单项退化为 `scale(c, pauli_word)`。

返回对象属性：

| 对象 | 属性 |
|---|---|
| PREPARE 各版（`StatePreparation`） | `width` = selector 位宽；模块属性 `prepare_terms` / `prepare_alpha` / `selector_width` |
| `select_pauli`（`Operation`） | 寄存器 `selector` / `target`；属性 `algorithm="select_pauli"`、`select_terms`、`selector_width` |
| `lcu_prepare_select`（`BlockEncoding`） | `be_alpha` = α；属性 `be_form="prepare_select"`、`lcu_terms`、`selector_width`、`prepare_work` |

## 实现要点

gate 版直接复用 `gate_state_prep` 的多路复用 Ry 旋转树（幅度 $\sqrt{|c_i|/\alpha}$ 补零到 $2^{\text{selector}}$）；QRAM 版复用 `qram_state_prep` 旋转树，角表由 `qram_state_angles` 生成。块编码的 signal 寄存器布局：低 `selector_width` 位是 selector，高位是 PREPARE 的 work（alias 版为纯化 junk）；PREPARE 需满足 `zero_input` 契约且目标宽度与 selector 位宽一致，否则生成期拒绝。

SELECT 的控制条件按 selector 二进制值用 RIR Control 原语表达；需要一元迭代时由后端把多比特控制降级实现，生成阶段不展开。产物可直接交给 `transforms.qubitization_walk`。适用边界：gate PREPARE 门数随项数指数增长；QRAM 版有角表量化误差（分布对拍容差 0.02）。

## 验证方案

类别 C5/C1：PREPARE / SELECT oracle 属数据访问层（C5），块编码组装的 (0,0) 块按 C1 与稠密矩阵**逐点对拍**（判定准则见 `../development/validation-plan.md` §2）。三层证据位于 `tests/core/test_prepare_select.py:PrepareSelectTests`（4 项测试哈密顿量 `TERMS`，α = 2.2）：

- 结构：`test_gate_and_qram_prepare_bindings_agree` 的声明与槽位断言；`test_pauli_hamiltonian_input_and_single_term_shortcut`（`PauliHamiltonian` 输入、单项 scale 捷径）；`test_bad_inputs_fail_at_generation`（单项 / 全零 / 非有限系数、SELECT 项数与字宽与字母违例、PREPARE 宽度不匹配等生成期拒绝）。
- 数值：`test_select_applies_indexed_pauli_word_with_phase`——SELECT 逐 selector 值与稠密 Pauli 矩阵 × 系数相位逐元对拍（places = 11）；`test_block_encoding_recovers_hamiltonian_over_alpha`——BE 的 (0,0) 块乘 α 后与稠密哈密顿量对拍（places = 10，经 `witness.assert_block_equals`）；`test_qram_prepare_block_encoding_approaches_hamiltonian`（QRAM PREPARE 版 delta = 0.02）；`test_qubitization_walk_composition`（直接接 `qubitization_walk`，signal = 0 块 = H/α，places = 10）。
- 绑定：`test_gate_and_qram_prepare_bindings_agree`——gate 版 selector 幅度与 $\sqrt{|c|/\alpha}$ 对拍（places = 11），abstract 槽经 `bind` 绑定 gate 实现后分布精确不变，QRAM 绑定经 `Binding` 挂角表后幅度对拍（delta = 0.02）；abstract → gate 在块编码内部的绑定由 `test_abstract_prepare_binds_inside_block_encoding` 覆盖（绑定后块不变，places = 10）。

## 已知缺口与计划阶段

模块级缺口集中在 alias 路线（见 [Alias 采样制备](alias-preparation.md)）；此外统一的三绑定一致性参数化（同一声明逐绑定对拍，oracles / prepare_select / graph_walks / data_loading 铺开）待 V4，现有 gate / QRAM / abstract 各自独立见证。与 `validation-coverage.md` 的 prepare_select.py 行一致。

## 相关链接

- 源码：`src/pyqecclang/algorithms/prepare_select.py`
- 同组页面：[Alias 采样制备](alias-preparation.md)、[态制备 Oracle](state-preparation.md)、[XOR 数据库](xor-database.md)
- API 参考：[PREPARE-SELECT 分解](../../api/algorithms/prepare_select.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
