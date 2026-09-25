# 开放声明、能力与绑定

<a href="../../reference/open-ir.html">English</a> · **简体中文**

本章定义 RIR 0.3 中开放声明的行为。开放主体自 0.2 引入，0.3 保留同一语义。算法协议与契约视角的背景见手册[算法自己的约定：从一个 gate 开始](../manual/contracts.md)。

## 开放声明不是恒等操作

Module.body 有两种形式：

| body | 含义 |
|---|---|
| null | 有签名但尚未实现的 oracle 声明。 |
| 指令数组，包括空数组 | 已定义实现。空数组明确表示恒等操作。 |

开放声明必须有字符串属性 oracle_paradigm。寄存器、资源与标量属性依然必须明确。声明可以出现在模块调用、Control、Adjoint 和 Repeat 中，并参与正常的类型与别名检查。

Program 可以保存部分实现。生成、验证和文本序列化不要求所有声明已绑定。执行与后端导出要求入口可达的 oracle 已闭合，否则报告槽名、范式与调用路径。未被入口调用的开放声明不阻止该入口导出。

所有调用目标仍然必须存在。因此开放调用引用的是明确声明，不是拼写错误或未知符号。

## 能力与等距制备

属性 supports_adjoint 和 supports_controlled 是布尔值。缺省值保持与 0.1 一致。具体模块的能力由自身声明及被调用模块的能力保守合取。

如果调用位于控制或伴随语境，验证器检查相应能力。等距态制备可以声明仅提供制备，也可以要求提供完整 unitary 扩张以及它的 inverse/controlled 版本。数学上的等距制备不再被定义成“物理不可逆制备”。

zero_input 是库契约。当前阶段不证明调用位置上的量子态为零，也不证明辅助空间已复净。所有实际量子位依然显式保留在接口中。

## 绑定

{obj}`bind(program, mapping) <oracq.infrastructure.linking.bind>` 返回新的 Program，不修改输入。mapping 的键必须是开放声明的名字，值为 {obj}`Operation <oracq.infrastructure.builder.Operation>` 或 {obj}`Binding <oracq.infrastructure.linking.Binding>`。

实现的寄存器参数按位置匹配声明的 kind 和 width。实现可以使用不同的局部参数名。明确声明的范式、能力和 be_alpha 必须兼容。绑定的实现自身可以包含新的开放依赖，因此绑定不等于全程序立即闭合。

实现通过一个同名包装模块接到原来的槽位。原始调用模块和接口用途仍然可见，实际实现不会被直接内联复制。

`oracle_paradigm` 是可扩展的接口角色标识；`PARADIGMS` 列出内置约定，应用
可用自定义标识符声明新角色。带有角色标记的具体实现须兼容声明；历史上没有
角色标记的操作仍按 ABI 和能力检查。角色名不会自动证明矩阵或函数的数学性质。

声明包含 `fixed_width`、`fixed_fraction`、`fixed_signed`、`rounding` 时，绑定
要求对应实现属性一致；`be_alpha` 同样必须一致。实现显式否定声明所需的
`zero_input` 或 `clean_work` 时拒绝绑定。缺少这些承诺的历史实现仍可绑定，
因此绑定成功表示结构兼容，不表示已经验证零输入或复净。

{obj}`bind_with_report(program, mapping) <oracq.infrastructure.linking.bind_with_report>` 只执行一次链接，返回 {obj}`BindingResult <oracq.infrastructure.linking.BindingResult>`。
`result.require()` 取得程序；`result.report.to_dict()` 输出输入和输出程序的
SHA-256、显式资源映射、剩余依赖和结构化错误。错误含槽位调用路径、期望和
实际值。部分绑定可以成功且仍有剩余依赖。报告不保存 Python 回调，也不包含
QRAM 数据内容；运行实验须另行记录内存快照的指纹。

同名但不同定义、循环依赖、错误参数数量或不兼容类型都会被拒绝。这些属于组装结构错误，不是算法数值正确性的判定。端到端绑定操作示例见教程[给算法替换 oracle](../tutorials/oracle-binding.md)。

## QRAM 捕获

抽象槽可以不声明 QRAM 资源，而它的具体实现可以需要资源。Binding.resources 将实现的资源参数映射为入口的逻辑资源名：

```python
bound = bind(open_program, {
    "A": Binding(implementation, {"table": "matrix_values"})
})
```

链接器会显式增加中间模块的捕获参数，更新沿调用链的资源实参，并把 matrix_values 放到入口资源表。相同逻辑资源可以被多个只读 oracle 共享。类型不一致或入口名字冲突会报错。

资源捕获只绑定句柄身份，不把数据表嵌入 IR。运行时仍然单独提供内存内容。分批绑定时，已提升的资源会保留。

新生成的捕获使用 `binding_captures` 属性保存来源；独立绑定顺序变化时，
新增资源及对应调用实参按逻辑名排列。序列化往返后继续绑定同一 bank 会复用
已有捕获。属性格式与跨节点检查见 [RIR 规范](rir.md)。

## 开放程度的边界

开放的是实现主体，签名中的宽度和用于组装的常量仍然具体。改变 BE 的 M 或 alpha，需要重新运行相同的 Python 生成器，得到对应的新 IR。链接器不会在已经生成的 LCU 角度里偷偷替换 alpha。

这符合当前的两阶段模型：Python 决定结构和配置，RIR 保留未提供的量子实现。任意符号宽度、延迟常量表达式求解和普遍依赖类型推断不属于本版。

RIR 允许多个寄存器和多个 QRAM bank。banked_database 演示了 96 位逻辑数据字在 64/32 位寄存器中的表示。部分算法便捷接口仍将信号打包为一个不超过 64 位的寄存器；更大实例应拆分信号接口，这属于库接口扩展，不应扩大 PySparQ 单寄存器位宽。

## 旧规则迁移

coverage.md 对 61 个旧正例和 16 个旧负例逐项说明。它记录范式表达路径，不宣称旧 .qec 逐字迁移或旧 CLIR golden 等价。Python 高阶函数、闭包与显式资源绑定取代旧 require 和禁止部分应用等规则。

末端测量/重置采用 ReadoutAction 宿主计划。正常 RIR 与 OriginIR 导出仍保留模块。当前 UnifiedQuantum 的动态解析器不接受 DEF；只有显式宿主读出适配会在最后一步展平已闭合模块，再交给动态解析器。目录同时保存模块化程序和这份执行适配文本。
