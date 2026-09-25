# OriginIR-ext 与 PySparQ 审阅

<a href="../../en/reference/backend-compatibility.html">English</a> · **简体中文**

审阅日期为 2026-09-08。实际提交记录在 backend-revisions.json。以下结论来自本地实现、测试和真实兼容性实验，没有修改两个上游仓库。

## OriginIR-ext 的优点与边界

OriginIR-ext 支持 QRAMDECL、DEF/ENDDEF、标量参数、named registers、受控和伴随结构。它适合作为 oracq 的第一个线路交换后端。oracq 侧的适配器见 API 参考[OriginIR-ext 后端](../api/infrastructure/backends/originir.rst)。

然而，DEF 和命名寄存器在当前实现中主要属于文本书写面。OriginIR_BaseParser._expand_def_call 会把调用内联，Circuit 的常用内部形式是 opcode_list。命名寄存器也在解析后映射为全局物理量子位。因此它不能直接承担本语言的模块化 register-level IR。

| 检查项 | 实际行为 | oracq 的处理 |
|---|---|---|
| DEF 参数 | 可以声明多个固定宽度寄存器。 | 每个 RIR 模块生成一个或多个 DEF。 |
| 嵌套 DEF | 可以展开，但内部调用必须正确重映射实参。 | 所有调用使用显式逐比特实参。 |
| DEF 中的 whole-register 实参 | 已记录的验证中最小实验中，嵌套调用 load(a,d) 被解析器拒绝。 | 导出 load(a[0],d[0],...)，保留模块调用。 |
| CONTROL 块 | 文档要求 ENDCONTROL 携带列表，但当前 parse_line 对它返回空量子位，实际解析仍失败。 | 使用逐门 controlled_by；受控模块通过附加控制参数特化，保持 DEF。 |
| QRAMDECL | 声明地址宽度与数据宽度；内存值由运行时提供。 | 类型声明进入导出文本，数据单独传给执行适配器。 |
| QRAM 查询 | 地址不变，数据目标执行 XOR。最低位按列表首元素解释。 | 接口和对拍覆盖非零目标及叠加态地址。 |
| QRAM 的模块形参 | DEF 没有资源句柄形式参数。 | 根据实际全局 QRAM 名称特化模块定义。 |
| 重复结构 | 当前静态子集没有保留任意 Repeat 的执行 IR。 | 导出对数大小的辅助 DEF 调用图。 |
| 回读与再导出 | Circuit 再导出会输出扁平线路。 | oracq 序列化文本（YAML/JSON）才是结构往返的依据。 |

审阅位置包括：
- UnifiedQuantum/uniqc/circuit_builder/originir_ext_spec.py。
- UnifiedQuantum/uniqc/compile/originir/originir_base_parser.py 中的 _expand_def_call、_process_statement 和 _apply_op。
- UnifiedQuantum/uniqc/test/core/test_originir_def.py、test_originir_named_registers.py 和 test_qram.py。
- UnifiedQuantum/uniqc/simulator/base_simulator.py 中的 _register_qrams。
- UnifiedQuantum/docs/source/1_basic_usage/originir.md。

### 数据注入和执行预算

Simulator.simulate_preprocess 创建 QRAM 对象。适配器填入 qram_objects 的数据，再调用 simulate_statevector。当前实现会保留相同声明对应的数据。

适配器关闭 least_qubit_remapping，保证状态向量索引与导出 artifact 的寄存器映射一致。默认最多模拟 24 位，实际项目通常应选择更小实例。QRAM 容器额外限制地址与数据位宽总和不超过 30；这不是 RIR 的位宽限制。

在将文本交给解析器之前，适配器计算保守展开预算。默认最多一百万个估计步骤。因而一个 Repeat(2^40) 可以紧凑存储和导出，但不会意外在 UnifiedQuantum 中展开执行。

已记录的验证中遇到的解析错误会被 Simulator 的格式自动回退包装成 OpenQASM 寄存器错误。调试时应直接检查 OriginIR_BaseParser，而不能只根据最终异常判断源文本类型。

## PySparQ 的寄存器模型

PySparQ 使用命名寄存器与全局注册表，每个基态保存寄存器的整数值。AddRegister 的 C++ 构造器拒绝超过 64 位的寄存器。这个限制不约束全系统总位数，也不能被解释为最多 64 个寄存器。oracq 侧的适配器见 API 参考[PySparQ 后端](../api/infrastructure/backends/pysparq.rst)。

RIR 的 Bits、UInt、SInt、Rational 分别对应 General、UnsignedInteger、SignedInteger、Rational 存储解释。Rational 表示无符号字除以 2^width，并不等于任意精度的 QFixed 类型。当前位操作处理原始位模式，add_const 只作用于 UInt，语义为模 2^width 加法。

PySparQ 的 SplitRegister 会把父寄存器的低位取出，并改变父寄存器的存储与宽度。CombineRegister 按对应低位规则拼回。RIR 的切片只是不可变视图，不直接调用这些有状态的结构变换。两者在位序上相容，但生命周期语义不同。

原生控制 API 的单个同类条件调用可能覆盖前一次条件。适配器将所有正控制位一次传给 conditioned_by_bit 的列表重载；零控制通过调用前后翻转相关位实现。

### 全局注册表和 QRAM 视图

System 注册表是进程全局的。适配器要求接管前注册表为空，拒绝清除已有模拟的寄存器。适配器使用自身互斥锁，并在 finally 中清理自己创建的状态；这不意味着可以与其他不使用此锁的 PySparQ 代码并发修改注册表。

对任意 QRAM 视图，适配器将地址按低位顺序 XOR 复制到临时整数寄存器，执行原生 QRAMLoad，把结果按控制条件 XOR 到原目标，再执行逆查询和逆复制。临时空间复净后才移除。它支持非零目标和跨寄存器视图，不依赖未知量子态的复制。

本版适配器物化 QRAM 数据时限制最多 2^20 项。RIR 本身只保存声明，并可以使用稀疏字典表示外部内存；这个执行限制不反向改变 IR。

审阅位置包括：
- QRAM-Simulator/PySparQ/pysparq/_core.pyi。
- QRAM-Simulator/PySparQ/pysparq/operators/condition_mixin.py。
- QRAM-Simulator/SparQ/src/system_operations.cpp。
- QRAM-Simulator/PySparQ/test/test_register_split.py、test_register_capacity.py 和 test_controlled_dagger.py。

## 已记录的验证中兼容性验证

真实测试覆盖嵌套 DEF、资源绑定特化、叠加态地址、非零数据目标、QRAM 视图、零控制、寄存器模加、Repeat、Adjoint 和复系数 LCU。比较对象为完整复幅度，不仅是测量概率。

语法导出支持模块化并不代表下游执行器具有模块级资源复用。后续若需大规模原生模块解释，应扩展下游执行架构或使用直接的 PySparQ 解释，不通过一次完整扁平化绕过该问题。

### 已复现的控制解析缺陷

在审阅提交中，OriginIR_LineParser.parse_line("ENDCONTROL q[0]") 返回的量子位字段为 None，而 OriginIR_BaseParser._apply_op 随后会遍历该字段。这会触发 TypeError。适配器没有修改上游，而是生成 inline controlled_by。对受控模块调用，适配器生成带额外控制形式参数的 DEF，并将这些参数传到内部逐门控制。Repeat 的辅助模块采用同样方式。

PySparQ 执行适配器默认限制最多 65536 个稀疏基态，并在原生门调用后检查。调用结束或异常退出都会清理适配器拥有的全局寄存器。


## 2026-09-09 应用描述验收补充

已记录的验证中 33 个绑定案例通过真实 OriginIR-ext 解析。没有以数值结果作为现有实现验收条件。

BaseParser 不识别 RESET，动态解析器虽支持 RESET 却不接收 DEF。因此末端测量/重置由独立 ReadoutAction 描述；正常 program.originir 保留 DEF。显式宿主读出适配才把已闭合模块展平，生成 execution-with-readout.originir 并交给真实动态解析器。这个局部适配不改变 RIR 或普通应用的模块化导出。

## 2026-09-09 自定义算子路径

补充审阅 dynamic_operator/{compiler,operator_wrapper}.py 与 SparQ/include/basic_components.h。动态算子的控制通过 split_systems/combine_systems 接线；生成的 C++ 直接访问实例 registers 存储，避免动态共享库的静态名字注册表与 Python 核心分离。不是 CACHED_REGISTER_SIZE 的固定容量问题。真实受控调用、逆、视图、私有工作区和 Roe 算术已运行。QRAM 逻辑 patch 由 oracq 管理，当前原生 QRAM bank 更新采用重新物化。细节见 [第二阶段说明](../manual/backends.md)。
