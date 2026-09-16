# 布尔网络（Boolean Networks）

> 类别 C1 · 模块 `pyqecclang.algorithms.arithmetic` · 阶段 V1

## 概述

布尔网络是可逆算术的中间表示：一个顺序有向无环的 SSA 图，节点为常量（编号 0/1）、输入位与不可变的布尔值（`not` / `xor` / `and`），在图上以经典组合逻辑（逐位加法、移位累加乘法、带恢复余数的除法、逐位开方等）描述算法，再统一编译为 compute/XOR/uncompute 形式的可逆线路，或生成同一张图的原生 C++ 算子。[可逆定点算术](fixed-point-arithmetic.md)、mathfunc 定点降级以及 Roe / QFVM / QHAM 等应用都以它为算术底座。

## 接口与输入模型

```python
BooleanNetwork()
net.operation(name=None, *, attributes=None)
BooleanNetwork.from_payload(value)
arithmetic_native_registry(program, *, cache_dir="out/native-cache")
```

`BooleanNetwork` 的构图层 API（节选，均返回节点编号或位列表）：

| 方法 | 含义 |
|---|---|
| `input(name, width)` / `const(value, width)` / `resize(bits, width)` | 声明输入位、常量与位宽调整 |
| `inv` / `xor` / `and_` / `or_` / `any` | 基本布尔代数（含 $a \oplus a = 0$ 等化简与规范化去重） |
| `mux(select, yes, no)` | 按位选择 |
| `add` / `neg` / `sub` / `lt` / `abs` / `mul` / `div` / `sqrt` | 逐位组合算术 |
| `evaluate(**inputs)` | 图的经典求值（输入为整数，输出各端口的整数值） |
| `payload()` / `from_payload(value)` | JSON 序列化与严格校验的反序列化 |

`net.operation()` 返回 `Operation`，寄存器为全部输入与输出端口（位宽即端口宽），并写入属性 `arithmetic_network`（payload JSON）、`correctness = "pending"`、`workspace_contract = "zero_in_zero_out"`；模块名缺省为 payload 哈希前缀。`arithmetic_native_registry` 扫描程序中带 `arithmetic_network` 属性的模块，端口与 payload 不一致时抛 `ValidationError`。

## 实现要点

编译策略是"计算—拷贝—复净"：每个非常量非常输入节点分配到宽度 64 的私有 bank（`ssa_0`, `ssa_1`, …），`not`/`xor` 用 X 门与 XOR 拷贝实现、`and` 用受控 X 实现，逐节点写入 bank；随后把输出位 XOR 拷贝到公开输出寄存器；最后对前向已发出的帧整体取 `Adjoint` 复净全部 bank，满足零进零出契约，调用方无需管理脏工作区。图构造侧以元组键缓存节点实现 hash-consing：`xor`/`and_` 对操作数排序规范化，配合代数化简（常量吸收、`not(not a) = a`）让相同子表达式共享同一节点。

`from_payload` 把网络当不可信输入校验：常量前缀、端口命名、位宽 1..64、位引用范围、顺序 DAG 性质与输入端口映射一致性，任一违例抛 `ValidationError`。原生路线 `BooleanCppFactory` 为同一张图生成真实 PySparQ C++ 算子：按调用点的实际布局（跨寄存器视图的 `start/width`）缓存编译产物，代码把节点值算进 `v[]` 数组后对输出寄存器做 XOR 写回，与门级线路语义一致。

## 验证方案

类别 C1（精确离散语义，判定准则见 `../development/validation-plan.md` §2）：编译出的线路须与图的经典求值逐点相等。三层证据：

- 结构：`tests/core/test_stage2.py:Stage2StructureTests`——`test_arithmetic_construction_and_basis` 对全部 14 种算术 `kind` 生成的网络模块检查私有工作区存在、程序 JSON 往返相等、Toffoli/U3/CZ 基导出只含三种基门且无 `controlled_by`。
- 数值：`Stage2StructureTests.test_gate_network_is_executable_small_example`——`add` 模块的 `arithmetic_network` payload 经 `from_payload` 复原后 `evaluate(a=1, b=2) == {"out": 3, "status": 0}`，且模拟器上同一程序的幅度精确集中在对应基态（线路与经典求值对拍）。
- 绑定：本机制无独立绑定见证（矩阵口径为 —）。原生执行路径在 `tests/integration/test_stage2_native.py:Stage2NativeTests.test_compiled_arithmetic_views_controls_and_adjoint`（真实 PySparQ 上原生算子与参考模拟器幅度相等，跨寄存器视图 + 受控 + 伴随场景，`native_calls = 2`）与 `test_native_module_private_workspace`（私有工作区经伴随调用复净）覆盖，属 L4 冒烟。

## 已知缺口与计划阶段

无已知缺口（验证矩阵缺口列为 —），阶段 V1。

## 相关链接

- 同模块：[可逆定点算术](fixed-point-arithmetic.md)
- 源码：`src/pyqecclang/algorithms/arithmetic.py`
- API 参考：[可逆算术](../../api/algorithms/arithmetic.rst)
- 验证矩阵：[验证覆盖矩阵](../../development/validation-coverage.md)

## 数值验证

论文级数值验证脚本：`tests/verification/verify_oracles.py`（`dj-boolean-network-*` 案例），产物 `out/verification/oracles.json`。

实验设计：把 BooleanNetwork 编译产物当作查询算法的 oracle 使用，端到端检验"编译线路 ≡ 图经典求值"。两个 3 bit 网络——parity3（`xor` 节点构成的平衡函数）与 const3（常量 1）——经 `operation()` 编译后包装为 XOR database 接口（address/data，输出为 XOR 拷贝语义，私有 bank 由伴随帧复净）：先在 address/data 全叠加下用 rir-pysparq 与 originir-ext 各跑一次，与 `net.evaluate` 逐点生成的真值表做全部 16 个分支的逐振幅对拍；再接入 [Deutsch–Jozsa](deutsch-jozsa.md) 电路做常量/平衡判定。

| 案例 | 规模 | 路径 | 指标值 |
|---|---|---|---|
| dj-boolean-network-parity3 | 3 bit 平衡，16 分支 | rir-pysparq, originir-ext | max_error = 8.3e-17，p_zero_error = 0.0 |
| dj-boolean-network-const3 | 3 bit 常量，16 分支 | rir-pysparq, originir-ext | max_error = 8.3e-17，p_zero_error = 1.1e-15 |

复现命令：

```bash
PYTHONPATH=src <含 pysparq+uniqc 的解释器> tests/verification/verify_oracles.py
```

产物路径：`out/verification/oracles.json`。
