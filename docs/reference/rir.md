# RIR：模块化寄存器级中间表示

本规范定义 RIR 0.3。开放声明、能力与绑定的详细规则见[开放 IR](open-ir.md)。

规范日期：2026-09-09。Python API、序列化器、参考执行器与后端适配器必须遵守本文。JSON Schema 只覆盖对象形状；跨节点约束由 validate 检查。

## 1. 设计范围

RIR 表示已经完成编译期参数求值的量子操作。所有寄存器宽度、整数常量和旋转角度均已具体化。RIR 保留模块定义、模块调用、静态重复、控制及伴随块，不要求展开成量子位级线路。

RIR 的已实现主体由酉操作构成；未实现模块以开放声明表示。它不包含测量、重置、运行期经典反馈或不透明 Python callback。对固定的外部 QRAM 内容，所有合法指令在完整量子空间上具有下文规定的酉语义。

公开工作区和信号寄存器出现在调用接口中；模块还可以通过 locals 声明私有零输入、零输出工作区。结构验证不提供任意辅助位复净证明。

## 2. 对象模型

Program 具有 entry、modules 和 version 三个字段。当前 version 为字符串 "0.3"，读取器也接受 "0.1" 和 "0.2"。entry 必须引用一个已定义模块。

Module 具有 name、registers、resources、body、attributes 和 locals。registers 是有序量子参数列表，resources 是有序 QRAM 参数列表，body 是有序指令体，或表示开放声明的 null。attributes 是由键和值组成的有序二元组列表。属性值只允许字符串、整数、有限浮点数或布尔值。

模块名在 Program 中唯一。寄存器名和资源名在各自模块内唯一且不能相互冲突。标识符匹配 [A-Za-z_][A-Za-z0-9_]*。模块至少具有一个非空量子接口，寄存器参数本身可以含零宽度项。

属性不改变指令语义。库可以基于属性约定数学解释。例如 be_alpha 与库规定的零投影布局共同定义块编码的尺度，但执行器不根据属性额外缩放量子态。

## 3. 存储类型与位序

Register 由 name 和 RegType 构成。RegType 由 kind 和 width 构成。width 必须是严格的 Python/JSON 整数，布尔值不能冒充宽度。

| kind | 位模式解释 |
|---|---|
| bits | 不指定数值意义的位串，对应 PySparQ General。 |
| uint | 无符号整数，范围为 0 到 2^width−1。 |
| sint | 二补码整数，位操作仍作用于原始字。 |
| rational | 无符号字除以 2^width，对应 PySparQ Rational。 |

一个寄存器或视图的宽度必须处于 0..64。零宽度表示空接口；原生适配器不为它创建实际寄存器。总量子位数和寄存器数量不以 64 为上限。

寄存器下标零表示最低位。假设入口签名依次声明 r0、r1 等寄存器，那么状态向量导出索引满足：

```text
index = value(r0) + (value(r1) << width(r0)) + ...
```

这个约定只定义向量化和后端物理映射。RIR 本身不保存全局物理位号。

## 4. 寄存器引用与视图

Ref 由 parts 和 type 构成。parts 是 Span 的有序列表。Span 具有 register、start 和 width，表示当前模块中一个根寄存器的连续范围。

Ref.parts 按低位到高位拼接。各段宽度的和必须等于 Ref.type.width。每个段都必须处于根寄存器范围内。同一 Ref 中的段不能引用重复量子位。

切片、fuse 和 reinterpret 在 Python 生成阶段形成 Ref，不生成量子指令。切片产生 bits 解释；reinterpret 明确改变解释但不改变位宽或量子状态。fuse 可以跨根寄存器，但合并后的视图仍然不能超过 64 位。

如果 f 是某模块的形式寄存器，调用实参是跨根寄存器的 Ref，那么 f 的局部切片必须映射为该 Ref 相同低位偏移的子视图。映射不能假定实际参数是连续物理量子位。

RIR 引用采用词法名字绑定。它们是声明式的引用，不是对 Python 变量生命周期或对象唯一引用的证明。

## 5. 指令集合

所有指令均有确定的参数字段。未知指令不是合法扩展点，必须先修订规范并提供后端行为。

### 5.1 Primitive

Primitive 具有 op、operands、angle 和 value。operands 是 Ref 列表。没有使用的 angle 或 value 必须为 null。

| op | 操作数 | 参数 | 语义 |
|---|---|---|---|
| h、x、y、z、s、t | 一个寄存器视图 | 无 | 对视图的每一位按低位到高位广播标准门。 |
| rx、ry、rz | 一个寄存器视图 | angle | 对每位施加 exp(−i angle P/2)。 |
| phase | 一个寄存器视图 | angle | 对每位施加 diag(1, exp(i angle))。 |
| gphase | 无 | angle | 对当前量子状态乘以 exp(i angle)，必须保留受控时的相对相位。 |
| xor | 两个同宽视图 | 无 | 保持第一个输入，第二个输入按位异或第一个输入。 |
| swap | 两个同宽视图 | 无 | 交换两个字的原始位模式。 |
| add_const | 一个 uint 视图 | value | 对该字执行模 2^width 加法。 |

所有角度均为有限实数，以弧度为单位。add_const 的 value 必须为 0..2^width−1 内的整数。二元操作的两个视图必须完全不重叠。

广播门在 IR 中仍然是单条寄存器级操作，不因为内部有多个物理门就变成多个 RIR 节点。零宽度上的广播、异或、交换和加零均为空操作。

### 5.2 Load

Load 具有 resource、address 和 data。resource 引用当前模块的 QRAM 形式参数。address 和 data 的宽度必须与声明一致，两个视图不能重叠。

QRAM 类型具有 address_width 与 data_width，两者均为 1..64。内存是由地址到无符号数据字的固定映射 M，未指定的单元为零。内存数据不写入程序 JSON，而是作为执行输入另行绑定。

```text
|address>|data>  ↦  |address>|data XOR M[address]>
```

Load 对任意数据目标成立，不要求目标初始为零。Load 自逆，且不修改地址或经典内存。

### 5.3 Call

Call 具有 module、arguments 和 resources。module 引用被调 Module。量子实参和资源实参均按被调模块的签名顺序绑定。

每个量子实参的 kind 和 width 必须与形式参数完全相同。需要改变位解释时，调用方必须显式 reinterpret。所有量子实参之间必须不重叠。

资源实参引用当前模块的资源名字，其 QRAM 类型必须与形式参数完全一致。多个只读资源形式参数允许绑定到同一个实际 QRAM。资源绑定不绑定量子地址或数据寄存器。

调用语义是将被调模块的所有局部引用替换为实际视图，在同一量子状态上执行主体。该定义不要求存储时内联主体。

### 5.4 Repeat

Repeat 具有 count 和 body。count 为 0..2^63−1 的整数。语义是顺序施加 body 共 count 次。count 为零时是恒等操作。

即使 count 为零，body 也必须结构合法。序列化和生成阶段不得根据 count 无条件复制指令体。

### 5.5 Control

Control 具有 register、value 和 body。register 必须是非空视图，value 必须处于其无符号位模式范围。

当控制视图等于 value 时施加 body，否则施加恒等。控制寄存器作为量子条件参与相干控制，不被测量或转换成 Python 条件。

控制位在整个 body 内受到保护。基元、QRAM 和模块调用的量子实参不能与它们重叠。模块调用采用保守规则，即使被调模块实际上没有修改某个形式参数，也不能将控制位作为该参数传入。

嵌套控制的视图不能相互重叠。不同控制条件按逻辑合取组合。gphase 不具有操作数，因此允许控制覆盖模块全部量子位；这时语义仍然是受控相位。

### 5.6 Adjoint

Adjoint 具有 body。其语义是将 body 的指令顺序反转，并对每条操作取伴随。rx、ry、rz、phase 和 gphase 的角度取负；add_const 转为模减法；Load、xor、swap、h、x、y、z 自逆；s 和 t 采用相应的逆相位。

Call 的伴随指向被调模块的逆操作，Repeat 的伴随重复其逆主体，Control 的伴随保持控制条件并对其主体取逆。双重伴随恢复原操作。

这些规则是指令的语义，不要求 RIR 在创建 Adjoint 时立即改写或展开其主体。

## 6. 调用图与模块复用

Program 的调用图必须无环。所有声明的模块都要检查，包括不可达模块。所有调用目标都必须存在。当前实现将模块调用深度和结构块嵌套深度限制为 127。

模块在各个调用点共享定义；调用不会复制签名或主体到 Program.modules。Python 生成器使用同一名字给出两个不同定义时必须报错。

每个调用目标都必须有明确的 Module 记录。该记录可以是开放声明；未定义的名字不等同于未完成的实现。

## 7. JSON 编码

每个数据类使用具有 tag 字段的 JSON 对象。tag 的取值对应以下记录名：

```text
Program, Module, Register, RegType, Span, Ref, QRAM, Resource,
Primitive, Load, Call, Repeat, Control, Adjoint
```

所有字段都必须写出，包括 null、空列表和空 attributes。解码器拒绝多余字段、缺失字段、未知 tag、重复 JSON 键、非有限浮点数以及未知版本。

Python 中的不可变元组编码为 JSON 数组。反序列化将其恢复为元组，不允许将任意 Python 对象反序列化成可调用代码。

以下示例是一条作用于两位整数寄存器的 H 广播指令：

```json
{
  "tag": "Primitive",
  "op": "h",
  "operands": [{
    "tag": "Ref",
    "parts": [{"tag": "Span", "register": "address", "start": 0, "width": 2}],
    "type": {"tag": "RegType", "kind": "uint", "width": 2}
  }],
  "angle": null,
  "value": null
}
```

规范输出采用 UTF-8、两空格缩进、按键排序和末尾一个换行。模块定义按模块名排序，签名参数和指令的列表顺序保留。浮点数采用 Python JSON 的有限浮点表示，禁止 NaN 与 Infinity。

模块 attributes 的顺序也会保留。Builder 按键排序属性；外部直接构造的 IR 应采用相同次序，以获得相同的规范输出。当前规范保证同一 IR 的确定性输出，不要求所有语义等价线路具有相同 JSON。

Schema 文件见 [rir.schema.json](schemas/rir.schema.json)。它描述 JSON 结构和局部范围；引用解析、别名、元数、跨节点类型、控制保护及调用图规则仍须运行语义验证器。对象形状的形式产生式汇总见第 8 节。

## 8. 形式文法

本节以产生式汇总第 2 节至第 7 节定义的对象形状，供独立实现对照。文法只覆盖结构与字段；别名、控制保护、调用图和数值范围等语义规则以正文和语义验证器为准。

词法约定：name 匹配 `[A-Za-z_][A-Za-z0-9_]*`；integer 是严格整数，布尔值不能冒充；float 是有限浮点数；string 是任意 JSON 字符串。`∅` 表示字段缺失（开放声明的空体，或未使用的 angle/value）。`X*` 表示有序不可变元组，允许为空。

```text
program     = Program { entry: name;
                        modules: module*;
                        version: "0.1" | "0.2" | "0.3" } .

module      = Module { name: name;
                       registers: register*;
                       resources: resource*;
                       body: instruction* | ∅;
                       attributes: attribute*;
                       locals: register* } .

register    = Register { name: name; type: regtype } .
regtype     = RegType { kind: "bits" | "uint" | "sint" | "rational";
                        width: integer } .
resource    = Resource { name: name; type: qram } .
qram        = QRAM { address_width: integer; data_width: integer } .
attribute   = name "×" (string | integer | float | boolean) .
span        = Span { register: name; start: integer; width: integer } .
ref         = Ref { parts: span*; type: regtype } .

instruction = primitive | load | call | repeat | control | adjoint .

primitive   = Primitive { op: gate-op;
                          operands: ref*;
                          angle: float | ∅;
                          value: integer | ∅ } .
gate-op     = "h" | "x" | "y" | "z" | "s" | "t"
            | "rx" | "ry" | "rz" | "phase" | "gphase"
            | "xor" | "swap" | "add_const" .
load        = Load { resource: name; address: ref; data: ref } .
call        = Call { module: name; arguments: ref*; resources: name* } .
repeat      = Repeat { count: integer; body: instruction* } .
control     = Control { register: ref; value: integer; body: instruction* } .
adjoint     = Adjoint { body: instruction* } .
```

angle 与 value 是否出现由 op 决定（见第 5.1 节），未使用者序列化为 null。开放模块的 body 为 ∅ 且不得声明 locals。版本 "0.1" 与 "0.2" 的 Module 不携带 locals 字段。

主要数值范围：RegType.width 处于 0..64；QRAM 两个宽度处于 1..64；Repeat.count 处于 0..2^63−1；Control.value 处于控制视图的无符号范围；add_const 的 value 处于 0..2^width−1；模块调用深度和结构块嵌套深度不超过 127。

JSON 编码把每条记录映射为携带 "tag" 的对象，字段名与记录字段一致；元组映射为 JSON 数组，标量与 null 原样传递：

```text
json(T, f1: v1, …, fn: vn) = { "tag": T, "f1": enc(v1), …, "fn": enc(vn) }
enc((e1, …, ek)) = [ enc(e1), …, enc(ek) ]
enc(∅)           = null
enc(标量)        = 标量
```

解码要求字段集合与 tag 记录完全一致，并拒绝未知 tag、多余或缺失字段、重复键、非有限数和未知版本（见第 7 节）。规范输出为 UTF-8、按键排序、两空格缩进、末尾一个换行。

## 9. 后端降低规则

OriginIR-ext 后端可以将寄存器操作降低为物理位操作，但必须保留模块定义与调用。QRAM 资源按实际绑定进行模块特化。Repeat 可以转换成共享的辅助模块图。

只有运行具体下游执行器时，才允许按其执行能力遍历或展开模块。下游自身可能展平，但这种处理不能反向改变 RIR 或替代 RIR 的结构序列化。

PySparQ 后端将非空入口寄存器映射到原生命名整数寄存器。对视图执行的临时重排或复制必须恢复其临时空间，不得修改 RIR 的根寄存器解释。

后端可以施加比 IR 更严格的执行预算，例如状态向量位数、QRAM 物化长度和展开次数。遇到限制时必须报错，不能静默截断 Repeat、量子位或内存内容。

## 10. 模块私有工作寄存器

Module.locals 是有序 Register 数组，不属于公开调用签名。每次调用从零态借入，必须在返回前复净；IR 只检查宽度和引用，复净是实现义务，模拟器提供运行期检查。开放模块不得声明 locals。Adjoint 和 Control 包括完整模块行为，工作区不能跨调用逃逸。OriginIR 导出为模块工作参数，顺序调用复用一段物理工作区；PySparQ 可在模块边界截获 native 实现，跳过其内部工作区与分解。原生注册表不是 IR 的一部分，不能将原生可执行误报为门级闭合。0.1/0.2 旧 JSON 仍可读写，其 Module 不含 locals；0.3 显式携带该字段。
