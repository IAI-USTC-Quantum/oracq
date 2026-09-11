# 态制备 Oracle（State Preparation）

> 类别 C5 · 模块 `pyqecclang.algorithms.oracles` · 阶段 V4

## 概述

幅值编码输入模型（input model 词汇中的 SP）的标准接口：范式 `state_prep_isometry` 把制备操作 $V$ 定义为从零态子空间出发的等距——$V|0\rangle|0\rangle = |\psi\rangle|0\rangle$，work 寄存器承诺复净；逆与受控操作要求存在可逆扩张 $U$。幅度向量按二叉树分解：每个节点把左右子块的权重写成

$$
\theta = 2\,\arctan\frac{\sqrt{w_{\mathrm{right}}}}{\sqrt{w_{\mathrm{left}}}} ,
$$

由多路复用 Ry 旋转树或 QRAM 角度表两种绑定实现。

## 接口与输入模型

```python
basis_state(width, value=0, *, work_width=0)
uniform_state(width, *, work_width=0)
gate_state_prep(amplitudes, *, work_width=0, name=None)
qram_state_prep(width, angle_width=8)
qram_state_angles(amplitudes, angle_width=8)
abstract_state_prep(name, width, work_width=0, *, reversible=True)
StatePreparation.from_unitary(operation, *, target=None, work=None, clean_work=False)
```

- `gate_state_prep`：任意复幅度向量（长度为二的幂、范数非零），返回门级旋转树制备。
- `qram_state_prep` / `qram_state_angles`：QRAM 资源版。前者只给位宽与角度位宽，后者由幅度向量生成角表字典（键为树节点下标），执行时作为 memory 提供；当前仅接受非负实幅度。
- `abstract_state_prep`：开放声明；`reversible=False` 时声明不支持 adjoint/controlled。
- `from_unitary`：把普通酉操作显式赋予 $U|0\rangle$ 初态角色；指定非零宽 work 时必须承诺 `clean_work=True`。

返回 `StatePreparation`（`OracleView`，`state_preparation()` 返回自身）。属性：

| 属性 | 含义 |
|---|---|
| `operation` | 底层 `Operation`（寄存器恰为 `target`、`work`） |
| `width` / `work_width` | target / work 位宽 |
| `capabilities` | adjoint / controlled 能力（来自声明或标注） |

模块属性：`zero_input=True`、`clean_work`；gate 版 `implementation="multiplexed_rotations"`，QRAM 版 `implementation="qram_rotation_tree"` 且附 `qram_queries=2*width`。

## 实现要点

gate 版：根节点施加无条件 Ry，其余节点按 target 高位前缀受控；复相位对每个非零相位基态经受控 `global_phase` 逐点补偿，因此幅度可为复数。至少需要一个目标位（单幅度向量拒绝）。

QRAM 版：寄存器 `target(width)`、`work(address_width + angle_width)`（`address_width = max(1, width)`）。逐层把 target 高位前缀 XOR 进地址、加树偏移后查询角表，按角度字的各个位受控施加 $R_y(2\pi k / 2^{\text{angle\_width}})$，再反向重放复净地址；每层两次查询。角度分辨率 $2\pi/2^{\text{angle\_width}}$ 是量化误差来源，与 gate 版的分布一致性在 prepare-select 侧以 delta = 0.02 对拍（见下节）。

`from_unitary` 的 `target=None` 分支把整个寄存器空间视为 target（work 宽 0）；`annotate` 对该范式缺省授予 `zero_input` / `clean_work` 承诺。适用边界：gate 版门数随维度指数增长；QRAM 版角表须由调用方与电路一同提供。

## 验证方案

类别 C5（数据访问层，判定准则见 `../development/validation-plan.md` §2）。三层证据：

- 结构：`tests/core/test_contracts.py:OracleContractTests.test_explicit_unitary_state_prep_adapter`（Hadamard 经 `from_unitary` 适配后满足 `StatePreparationProtocol` 契约，酉本体同样可作输入）；`test_unitary_adapter_requires_work_promise`（未承诺 `clean_work` 抛 `ContractError`）；`tests/core/test_open_ir.py:OpenIRTests.test_isometry_requires_declared_adjoint_capability`（`reversible=False` 声明在 adjoint 上下文中被拒）。
- 数值：`OracleContractTests.test_dirty_work_is_rejected_by_algorithm`——`clean_work=False` 的制备被 `linear_qode` 契约在运行内核之前拒绝（validation-coverage 的 oracles.py 行数值锚点）；QRAM 旋转树的数值行为由 `tests/core/test_qlss_input_models.py:QLSSInputTests.test_tree_preparation_queries_each_layer_coherently` 见证（8 次 Load 查询、归一化幅度对拍 delta = 0.07）。
- 绑定：`OracleContractTests.test_capability_restrictions_survive_annotation_and_wrapping`（能力限制在标注与包装后保留）；`test_legacy_state_prep_clean_work_in_bind`（缺失 `clean_work` 属性的旧制备仍可绑定闭合）；gate / QRAM 两实现的分布一致性由 `tests/core/test_prepare_select.py:PrepareSelectTests.test_gate_and_qram_prepare_bindings_agree` 经 `gate_prepare` / `qram_prepare` 组装间接见证（gate 版 places = 11，QRAM 版 delta = 0.02）。

## 已知缺口与计划阶段

三绑定一致性参数化（V4 铺开）：abstract / gate / QRAM 三层的统一参数化对拍未铺开，现有各绑定独立见证。与 `validation-coverage.md` 的 oracles.py 行一致。

## 相关链接

- 源码：`src/pyqecclang/algorithms/oracles.py`
- 同组页面：[PREPARE–SELECT 分解](prepare-select.md)、[Alias 采样制备](alias-preparation.md)、[XOR 数据库](xor-database.md)
- API 参考：[Oracle 声明与实现](../../api/algorithms/oracles.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)
