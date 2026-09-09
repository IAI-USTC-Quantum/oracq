# pyqecclang 完整规范

| 项 | 值 |
| --- | --- |
| 包版本 | 0.6.0（`pyproject.toml`） |
| RIR 版本 | 0.3（`pyqecclang.infrastructure.ir.VERSION`，兼容读写 0.1/0.2） |
| MIR 版本 | 0.1（数学函数图） |
| PDE / QCL plan 版本 | 0.1 / 0.1（QHAM 生成层） |
| 规范日期 | 2026-09-09 |
| 适用对象 | Python API、序列化器、参考执行器、后端适配器、CLI |

本文档是 pyqecclang 的统一总规范，把分散在各专题文档中的规范性内容合并为一份完整描述，覆盖语言核心、四层生成性中间表示、算法组装层与全部后端契约。各专题文档（[rir-spec.md](rir.md)、[open-ir.md](open-ir.md)、[math-ir-spec.md](math-ir.md)、[oracle-paradigms.md](../archive/oracle-paradigms.md)、[backend-review.md](backend-compatibility.md) 等）仍然是各自领域的详细规范；当本文档与专题文档出现出入时，以更接近代码现状的专题文档为准，并应当修订本文档。使用导向的教程与案例集见 [pyqecclang-guide.md](../archive/guide-0.6.md)。

## 1. 语言定位与设计原则

pyqecclang 是一门基于 Python 的量子操作生成语言。开发者编写普通的 Python 函数，这些函数在生成阶段执行，产出由模块定义、模块调用、静态重复、相干控制和伴随块组成的**寄存器级中间表示（Register-level IR，下称 RIR）**。RIR 保存模块结构而不展开成量子位门列表，这一点是整个设计的中心约束：宽度在生成期具体化，但被调模块的主体不会被复制进调用方。

本语言独立于 QECC.Lang，不继承其文本语法，也不受其扁平 CLIR 约束。它与后者的关系仅在案例层面：[coverage.md](../archive/coverage.md) 记录了 QECC.Lang spec-tests 的 61 个正例与 16 个负例如何映射到本语言的范式（该映射证明存在表达路径，不证明源码或 golden 等价）。

设计原则可以概括为五条，后续章节的每条规则都能追溯到其中之一。

**生成阶段与运行阶段严格分离。** Python 负责普通参数求值、函数复用、闭包和实现选择；生成结束后的 RIR 中不含任何 Python 可调用对象，不含包路径加载钩子，也没有运行期经典分支。普通的 Builder 流程不解析 Python AST；唯一的例外是 `compile_function`（第 13 章），它读取受限纯数学函数的源码并翻译成数学图，这一入口同样不引入新的文本语法。

**模块结构是一等公民。** RIR 中的模块调用、`Repeat`、`Control`、`Adjoint` 都是显式节点。生成与 JSON 序列化阶段不得按重复次数或调用图无条件展开指令；寄存器名、位宽和视图必须原样保留到后端降低阶段。只有当某个具体下游执行器需要遍历时，才允许按其自身能力展开，且展开不能反向改写 RIR。

**允许未完成的实现。** 一个只有签名、没有主体的模块是合法的 RIR（第 8 章）。开放的是实现主体，不是签名形状：宽度和组装常量在生成期已经固定，链接器不会在已生成的角度或系数里偷换参数。空的主体与"空线路冒充实现"有本质区别，后者被明确禁止。

**核心不依赖量子后端。** 语言核心零第三方运行依赖（Python ≥ 3.11，仅标准库）。后端分为描述侧（OriginIR-ext 文本导出）与执行侧（参考执行器、UnifiedQuantum、PySparQ），执行侧的原生依赖只在调用时导入。测量、重置与读出后选择属于宿主层，不进入语言核心；近似精度界（eps）、资源估计与算法收敛性同样不在语言核心之内。

**正确性分级声明。** 语言核心的结构规则、序列化确定性和后端描述可达性是已验证的；算法层面的数学正确性、数值精度、成功概率与复杂度优势一律标记为待核验（`correctness: pending` / `not_assessed`），并且不作为组装与导出的前置条件。第 19 章给出完整的一致性分级。

## 2. 总体架构与五层表示

一个 pyqecclang 程序从书写到执行经过五种表示形态，各层职责如下表。

| 层 | 表示 | 是否可序列化 | 所在模块 |
| --- | --- | --- | --- |
| 宿主生成层 | Python 函数、闭包、常量 | 否（生成期即消失） | 用户代码 + `builder.py`、`oracles.py` |
| 生成性中间表示 | MIR 0.1（数学函数图）、PDE 0.1、QCL plan 0.1 | 是（独立 JSON） | `mathfunc/`、`qham/` |
| 核心中间表示 | RIR 0.3 | 是（`dumps`/`loads`） | `ir.py`、`serialization.py` |
| 描述后端 | OriginIR-ext 文本；Toffoli/U3/CZ 严格门集文本 | 文本文件 | `backends/originir.py`、`backends/basis.py` |
| 执行后端 | 参考执行器（纯标准库）；uniqc 模拟器；PySparQ | 运行期 | `execution.py`、`backends/pysparq.py` |

生成性中间表示是可选的：直接用 Builder 手写模块时不存在 MIR 或 PDE；只有调用 `compile_function`（产生 MIR）或 QHAM 生成器（产生 PDE 与 QCL plan）时才会出现。它们的共同特征是不含量子门，只描述数学结构，且必须在降低到 RIR 之前完成自身的一致性检查。

后端的锁定关系记录在仓库根的 [backend-revisions.json](../../backend-revisions.json)：UnifiedQuantum 审阅提交 `5555b2c`（OriginIR-ext 解析与 uniqc 模拟器），QRAM-Simulator 审阅提交 `111630a`（pysparq System、动态算子编译）。pyqecclang 不通过包索引安装或更换这两个后端；集成测试要求解释器环境真实提供它们，不以缺少后端为理由跳过。

[architecture.md](../archive/architecture.md) 提供同一架构的图示索引；本章与该文档冲突时以本章为准并修订该文档。

## 3. 寄存器与类型系统

### 3.1 存储类型

每个量子寄存器（或寄存器视图）携带一个 `RegType`，由存储类 `kind` 与位宽 `width` 组成。存储类共有四种：`bits`（原始位）、`uint`（无符号整数）、`sint`（二补码有符号整数）、`rational`（无符号字除以 2^width 的定点解释，对应 PySparQ 的 Rational 存储）。构造函数 `Bits(w)`、`UInt(w)`、`SInt(w)`、`Rational(w)` 生成这四种类型。

存储类是 ABI 元数据。基元门一律按原始位语义作用：在 `SInt(3)` 或 `Rational(3)` 上发一个 `x` 门得到全 1 的位模式，与在 `Bits(3)` 上完全相同。存储类影响两件事：后端存储映射（PySparQ 侧 bits/uint/sint/rational 分别映射到 General、UnsignedInteger、SignedInteger、Rational），以及个别指令的类型约束（`add_const` 只接受 `uint` 目标）。语言不为定点解释内建数值精度证明。

### 3.2 宽度与位序

单个寄存器或合并视图的宽度范围是 0 到 64。零宽度寄存器是合法的空接口，不为它分配任何原生寄存器。寄存器数量与总量子位数不受 64 限制：一个模块可以持有多个 64 位寄存器。

位序约定与 PySparQ 一致：下标 0 是最低位。状态向量导出时按入口签名顺序打包，索引公式为 `value(r0) + (value(r1) << width(r0)) + …`，即第一个寄存器占据最低位。RIR 本身不保存全局物理量子位编号；这一映射只存在于向量化与后端物理布局阶段。

### 3.3 引用与视图

指令不直接引用寄存器，而是引用 `Ref`——一个由若干 `Span`（某根寄存器的连续位段 `[start, start+width)`）按低位到高位拼接成的视图，附带自身的 `RegType`。视图操作发生在生成阶段，不产生任何量子指令：

- `builder["name"]` 取整个寄存器的引用；
- `ref[i]` 与 `ref[a:b]`（步长必须为 1）取下标与连续切片，切片结果按 `bits` 解释，零宽切片合法；
- `ref.reinterpret(kind)` 改变存储类解释而不改变任何量子位；
- `fuse(ref_a, ref_b, …)` 跨根寄存器拼接，低位参数在前，结果宽度为各段之和且不得超过 64。

一个 `Ref` 内部的各 `Span` 不得占据重复的量子位。形参寄存器的局部切片在调用时映射到实参 `Ref` 的相同低位偏移子视图，因此不能假定实参是连续的物理量子位。引用是词法名字绑定，不是 Python 变量生命周期证明。

### 3.4 别名规则

同一条指令的多个操作数所占据的物理位集合必须两两不相交。违反时验证器报"操作数存在别名或重叠"。控制位受更强的保护：`Control` 块内的基元操作数、QRAM 地址与数据、以及被调模块的全部实参都不得与控制视图重叠——即使被调模块实际上不修改某个形参，也不能把控制位作为实参传入，这是保守规则。

### 3.5 QRAM 资源

`QRAM(address_width, data_width)` 描述一种量子可查询内存资源，两个维度各自取值 1..64。模块通过 `resources` 声明 QRAM 形参（`Resource(name, type)`），调用点通过名字把形参绑定到调用方的具体资源，绑定要求 `QRAM` 类型逐字段相等。QRAM 内存内容不写入程序 JSON：它是执行期输入，按入口资源名提供（密集字序列，或地址到字的稀疏字典），未指定的单元读作零。

## 4. RIR 对象模型

RIR 的全部节点都是不可变 dataclass（`frozen=True`），定义在 `ir.py`。下表是完整节点清单，字段名即 JSON 键名。

| 节点 | 字段 | 说明 |
| --- | --- | --- |
| `Program` | `entry, modules, version` | 入口模块名、模块集合、RIR 版本 |
| `Module` | `name, registers, resources, body, attributes, locals` | 见 4.1 |
| `Register` | `name, type` | 模块接口或局部寄存器声明 |
| `RegType` | `kind, width` | 存储类 × 位宽 |
| `Span` | `register, start, width` | 根寄存器内的连续位段 |
| `Ref` | `parts, type` | 视图：parts 按低位到高位拼接 |
| `QRAM` | `address_width, data_width` | QRAM 资源规格 |
| `Resource` | `name, type` | 具名 QRAM 形参 |
| `Primitive` | `op, operands, angle, value` | 基元门；operands 为 Ref 元组 |
| `Load` | `resource, address, data` | QRAM 异或载入 |
| `Call` | `module, arguments, resources` | 模块调用；resources 为资源名字符串元组 |
| `Repeat` | `count, body` | 静态重复 |
| `Control` | `register, value, body` | 相干控制 |
| `Adjoint` | `body` | 伴随 |

### 4.1 Module 与 Program

`Module` 的 `registers` 是有序的公开量子形参表，`resources` 是有序的 QRAM 形参表，两者内部名字唯一且互不冲突。公开寄存器的总宽必须大于零（单个寄存器允许零宽）。`body` 有两种形态：指令元组表示已定义的实现，空元组是明确的恒等模块；`None` 表示开放声明（第 8 章），此时必须携带字符串属性 `oracle_paradigm`，且不得声明 `locals`。一个主体为空元组的普通模块不是开放操作——`Builder(...).finish()` 产出的就是合法的空实现。

`attributes` 是 `(key, 标量)` 二元组的有序元组。属性值只允许字符串、整数、有限浮点数与布尔。属性不改变指令语义：`be_alpha` 与库约定共同定义块编码尺度，但执行器不会据此缩放量子态。

`locals`（RIR 0.3 新增）是模块私有工作区寄存器的有序表，不属于公开签名。每次调用时局部寄存器从零态借入，模块必须在返回前把它们复净为零。IR 只检查宽度与引用合法性，复净是实现义务；参考执行器与 OriginIR 后端在模块边界做运行期检查。工作区不能跨调用逃逸信息。

`Program.entry` 指向入口模块；`module_map` 与 `main` 提供按名访问。模块名在 Program 内唯一，同名不同定义会在 `Operation.program()` 合并依赖时被拒绝。

### 4.2 Builder 与 Operation

`Builder(name, registers, resources=None, *, attributes=None)` 是唯一的常规构造入口。常用方法：`local(name, type)` 声明私有工作区；`__getitem__` 取视图；`gate/h/x/z/ry/rz` 发基元门；`xor/swap` 发二元操作；`add_const(target, value)` 发常量加法（目标必须是 `uint`，value 取 0..2^w−1，按模 2^w 回绕）；`global_phase(angle)`；`qram(resource, address, data)` 发 `Load`；`call(operation, *, resources=None, **arguments)` 发 `Call`，关键字名集合必须恰好等于目标的寄存器名集合与资源名集合；`repeat(count)`、`control(register, value=None)`（缺省值为全 1）、`adjoint()` 是三个上下文管理器，各自产生显式 IR 节点；`finish()` 冻结并返回 `Operation`。

块机制带有回滚语义：`with` 体内抛出异常时弹帧且不发射任何指令。全部验证延迟到 `finish()`（经 `validate`），构造期的形状错误在那里统一报告。`finish()` 之后构造器冻结，再次调用或继续发射都会报错。

`Operation(module, dependencies)` 携带一个模块及其依赖集合；`program()` 把依赖按名去重合并、排序、验证后产出 `Program`。普通 Python 函数、闭包或可调用对象都可以作为 protocol 的实现（返回 `Operation`），但生成结果必须是可验证的 RIR 模块——不能把任意 Python callback 当作未定义的 IR 指令保存在程序里。

## 5. 指令语义

RIR 表示已完成编译期参数求值的量子操作：寄存器宽度、整数常量与旋转角度均已具体化。0.1 版语义范围是封闭酉操作子集——没有测量、重置、运行期经典反馈与不透明回调。对固定的外部 QRAM 内容，所有合法指令在完整量子空间上具有规定的酉语义。模块不隐式申请或释放量子寄存器；语言不提供任意辅助位复净的证明。

未知指令不是合法扩展点。新增指令必须先修订规范，再进入实现。

### 5.1 Primitive

| op | 操作数 | 语义 |
| --- | --- | --- |
| `h, x, y, z, s, t` | 1 | 单量子位门，在寄存器视图上逐位广播 |
| `rx, ry, rz` | 1 | 旋转门 exp(−i·angle·P/2)，angle 为有限实数弧度 |
| `phase` | 1 | diag(1, e^{i·angle})，逐位广播 |
| `gphase` | 0 | 全局相位 e^{i·angle}；受控时必须保留为受控相对相位 |
| `xor` | 2 | 第二操作数按位异或第一操作数（CNOT 链），两视图完全不重叠 |
| `swap` | 2 | 等宽交换 |
| `add_const` | 1 | `uint` 目标模 2^width 加常数，value ∈ 0..2^width−1 |

广播门是单条寄存器级指令：64 位寄存器上的一个 `h` 是一个 IR 节点。在零宽度视图上，广播、异或、交换与加零都是空操作。旋转门与 `gphase` 要求有限角度，其余基元的 `angle` 必须为 `None`。

### 5.2 Load

`Load` 实现 `|address⟩|data⟩ ↦ |address⟩|data ⊕ M[address]⟩`。内存是地址到无符号数据字的固定映射。该操作对任意数据目标成立（不要求初始为零）、自逆（连续两次载入同一地址为恒等）、不改地址寄存器与经典内存。每个相干计算段内 QRAM 内容必须固定。

### 5.3 Call

调用语义是把被调模块的形参引用替换为实参视图，不要求存储时内联。量子实参的 kind 与 width 必须逐字段相等（改解释须显式 `reinterpret`），实参互不重叠；资源实参解析到调用方声明的资源且类型完全相等。多个只读资源形参可以绑定同一个 QRAM。调用深度与结构块嵌套深度均限制在 127 以内。

### 5.4 Repeat

`count` 取值 0..2^63−1，count 为 0 时是恒等。即使 count 为 0，body 也必须结构合法。生成与序列化阶段不得按 count 复制指令体；后端可以把它编译为对数规模的辅助模块图（17.1）。

### 5.5 Control

`Control(register, value, body)` 在视图的无符号位模式等于 value 时施加 body。这是量子相干控制，不测量、不转为 Python 条件。控制位在整个 body 内受保护（3.4 节）。嵌套控制的视图不得重叠；不同控制条件按合取组合。`gphase` 没有操作数，允许控制覆盖全部量子位，此时它仍是受控相位。

### 5.6 Adjoint

`Adjoint(body)` 表示 body 的逆：逆序遍历并逐条取逆——旋转角取负、`add_const` 变为模减、`Load/xor/swap/h/x/y/z` 自逆、`s/t/phase` 逆为相位取负。`Call` 的伴随指向被调模块的逆，`Repeat` 的伴随重复逆主体，`Control` 保持条件取逆主体，双伴随恢复原操作。这些是语义规则，不要求在创建 `Adjoint` 时立即改写或展开主体。

## 6. 验证规则与硬限制

`validate(program, *, require_closed=False)` 是唯一的验证入口，定义在 `validation.py`。验证失败抛出 `ValidationError`（`ValueError` 子类），错误消息为中文并包含具体规则短语（如"重叠"、"递归"、"冲突"、"未绑定 oracle：…"），测试以消息子串断言。验证涵盖以下层次。

结构规则包括：标识符匹配 `[A-Za-z_][A-Za-z0-9_]*`；版本属于 {0.1, 0.2, 0.3} 且 0.1 不允许开放声明、`locals` 仅 0.3 携带；寄存器与资源命名唯一；公开接口总宽大于零；QRAM 两维各 1..64；属性键为标识符、值为标量、浮点有限、能力属性为布尔。

指令级规则包括：3.3 与 3.4 节的视图与别名约束；基元操作数个数（二元恰两个、`gphase` 零个、其余一个）与二元等宽；`add_const` 的类型与取值范围；`Load` 的资源声明与宽度匹配；`Call` 的 ABI 与资源解析；`Repeat` 计数范围；`Control` 寄存器非空、value 在位宽范围内、嵌套不重叠；嵌套深度小于 128。

图级规则包括：模块名唯一、入口存在、调用目标存在、调用图无环（对包括不可达模块在内的全部声明检查）、调用深度小于 128，以及能力需求核对——`Control` 语境中的调用要求被调者 `supports_controlled`，`Adjoint` 语境中的调用要求 `supports_adjoint`（详见第 9 章）。

`require_closed=True` 在以上基础上叠加开放缺口检查，所有导出与执行入口都使用它。非 `ValidationError` 的常见异常会被包装成 `ValidationError` 再抛出。

全量硬限制汇总见附录 B。

## 7. JSON 序列化与版本兼容

序列化定义在 `serialization.py`。编码规则为：每个 dataclass 编成 `{"tag": 类名, 字段…}`；元组成 JSON 数组；标量原样；其余类型报错。所有字段必须写出（包括 null 与空表）；解码器拒绝多余字段、缺失字段、未知 tag、重复 JSON 键与非有限浮点，也不允许把任意 Python 对象反序列化成可调用代码。根必须是 `Program`。

规范输出（`dumps`）是确定性的：UTF-8、两空格缩进、键按字母排序、末尾一个换行、禁 NaN/Infinity、`ensure_ascii=False`；模块定义按名排序，签名参数与指令顺序保留。同一 IR 的输出逐字节相等；语义等价的不同线路不保证同 JSON。`loads` 之后 `dumps` 得到相同文本，对象层面也相等——模块与调用边原样往返。

版本处理：0.1 与 0.2 的 JSON 没有 `locals` 字段，读入时回填空表，写出时删去该键；0.3 显式携带。`Builder` 的属性按键排序后冻结，保证确定性。

配套的结构校验是 [rir.schema.json](schemas/rir.schema.json)（Draft 2020-12，`Program.version` 枚举三种版本，`Module.locals` 为可选属性）；跨节点规则不在 schema 内，由 `validate` 检查。

## 8. 开放 IR 与链接

开放 IR（RIR 0.2 引入）允许程序携带尚未提供的实现。设计立场是：开放的是实现主体，签名宽度与组装常量仍然具体；任意符号宽度、延迟常量求解与普遍依赖类型推断不属于本版。

### 8.1 开放声明与缺口报告

`Module.body is None` 即开放声明，必须携带字符串属性 `oracle_paradigm`（取值见第 10 章），不得声明 `locals`。`declare(name, registers, *, paradigm, resources=None, attributes=None, supports_adjoint=True, supports_controlled=True)`（`oracles.py`）是标准构造入口，它额外写入 `implementation_status="unresolved"`。

`unresolved(program)` 先验证，再从入口做广度优先遍历（只走可达模块），对每个可达的开放声明产出一条 `OracleRequirement(name, paradigm, registers, attributes, path)`，其中 path 是从入口出发的最短调用路径。结果按名字排序。不可达的开放声明不阻塞导出。

### 8.2 绑定

`bind(program, bindings)` 返回新的 `Program`，不修改输入。键必须是当前开放声明的名字，值是 `Operation` 或 `Binding(operation, resources)`；实现模块名不得与槽同名。绑定前做四项核对：实现的寄存器类型元组必须与声明逐位相等（形状变化需要回到 Python 生成器重新生成）；声明携带 `be_alpha` 时实现必须相等；`oracle_paradigm` 必须一致（声明为 `unitary` 时豁免）；声明要求的能力（缺省为 True）必须是实现推断能力的子集。绑定的实现自身可以包含新的开放依赖——分批绑定是正常用法，部分绑定后其余槽位保持开放。

槽被改写为一条对实现模块的调用（同名包装，不内联复制），属性追加 `oracle_bound_to` 与 `implementation_status="bound"`。

### 8.3 QRAM 捕获提升

实现需要的 QRAM 资源若在抽象槽上未声明，链接器做显式捕获：全局资源名取 `Binding.resources[formal]`，缺省为 `{槽名}__{形参名}`；随后计算每模块的资源传递闭包（同时检测绑定引入的循环与未知模块），入口直接采用全局名（与已有资源或寄存器冲突时报错），中间模块以 `capture_<sha256(资源名)[:16]>` 为名新增形参并改写全部调用边。共享同一捕获资源的多个槽在入口只产生一个资源。捕获只绑定句柄身份，不把数据表嵌入 IR。

## 9. 能力系统

`supports_adjoint` 与 `supports_controlled` 是模块的布尔属性，缺省为 True。模块能力由自身声明与全部被调模块的能力保守合取推导（`capabilities(program)` 记忆化递归，穿透 Repeat/Control/Adjoint）。验证器在控制或伴随语境中发现调用时核对相应能力。

能力边界声明：当前指令集均为固定 QRAM 数据下的酉操作，因此默认可控制可取逆。今后引入测量或其他效果时必须扩展能力检查，不能沿用这一结论。等距态制备是典型例子：可以只声明制备能力（`abstract_state_prep(..., reversible=False)` 会同时关掉两种能力），也可以提供完整酉扩张；只有声明了逆能力的等距制备才能进入 Costa/Grover 一类的反射结构。`zero_input`（制备假设输入为零）与 `clean_work` 是库契约，不由编译器证明。

原生注册表（17.3）不是 IR 的一部分：一个模块可以在 PySparQ 里有原生可执行实现，但这不使它获得门级闭合的身份，`require_closed` 检查的是 RIR 主体。

## 10. Oracle 范式目录

`oracle_paradigm` 的合法取值共九种（`oracles.py` 的 `PARADIGMS`）。范式是接口约定，编译器检查的是结构；数学契约由库文档与验证记录支持。下表给出每种范式的接口与仓库已提供的实现。

| 范式 | 接口约定 | 已提供实现 |
| --- | --- | --- |
| `unitary` | 任意固定签名酉 | 用户自建模块 |
| `block_encoding` | target、signal，属性 `be_alpha` | Pauli/LCU 门级、对角角表查询、稀疏适配（第 11、14 章） |
| `database_xor` | address 保留，data 按位异或查表结果 | `gate_database`（真值表 X 网络）、`qram_database`（单条 Load）、`banked_database`（宽字多 bank） |
| `state_prep_isometry` | 从零输入子空间制备 target，work 显式保留 | `basis_state`、`uniform_state`、`gate_state_prep`（多路旋转树）、`qram_state_prep`（QRAM 角表树） |
| `sparse_location_inplace` | CKS 位置接口：给定 column，把第 index 个非零位置原地写入 index 寄存器 | `sparse_location_gate`（完整置换扩张的换位网络）、`sparse_location_qram`（正反向双 QRAM） |
| `sparse_entry_xor` | CKS 元素接口：任意 row/column，data 异或写入矩阵条目 | `sparse_entry`（包装 XOR 数据库） |
| `phase_oracle` | 谓词匹配基态加相位 | `phase_marks`、`phase_from_database` |
| `reversible_function` | 输入保留、输出可逆更新 | 算术模块、`reversible_lookup` 查表 |
| `algorithm_stage` | 自定义但固定的量子签名阶段 | 领域组装保留为声明 |

CKS 位置接口（arXiv:1511.02306 §1.1）值得单独说明：它是把给定列的第 l 个非零位置**原地写入** index 寄存器，需要有效逆映射；保留 l 的 XOR 查表不满足该接口，不能把 XOR 位置表直接当原地置换实现。仓库实现采用完整 n 位索引寄存器与完整置换扩张，稀疏度是编译期参数。

`oracles.py` 同时提供包装类型 `XorDatabase`、`StatePreparation`、`StateOracle`、`SparseAccess`（在构造时校验接口形状），以及 `annotate`（把已构造操作标记为某范式）、`invoke`（带前缀资源映射的调用发射）与 `resources_for`。`access.py` 提供适配层：`reversible_lookup`（查表函数的可逆调用）、`word_rotation`（查询字到旋转角的转导）、`batch_lookup`（分批查询）；其中 `sparse_block_encoding` 是历史候选，正式稀疏适配见 14.4 的 `real_symmetric_sparse_encoding`，它带有 `legacy_input_model` 标记。

## 11. 块编码代数

块编码库（`library.py`，公共 API）维护"信号位视角"的编码代数。`BlockEncoding(operation)` 要求模块接口恰好是 `{target, signal}` 且均为 bits、target 宽度大于零、`be_alpha` 为有限正数；暴露 `alpha`、`width`、`signal_qubits`。alpha 保存在模块属性 `be_alpha` 中，JSON 往返存活。

公共组合子的语义：

- `identity(width)`：恒等，signal 宽 0，α = 1。
- `pauli_x(width)`：target 上的 x，α = 1。
- `zero(width)`：零算子编码（signal 位上 x 门）。
- `scale(c, a)`：c 为有限复数；c 为 0 时退化为 `zero`；否则 `gphase(arg c)` 后调用 a，α = |c|·αa。
- `product(a, b)`：等宽；signal 拼接（a 高 b 低）；执行顺序为先 b 后 a；α = αa·αb；资源形参加 `a__`/`b__` 前缀。
- `linear_combination(ca, a, cb, b)`：等宽；α = |ca|·αa + |cb|·αb；制备角 θ = 2·acos(sqrt(|ca|·αa/α))；结构为 `ry(select, θ)`、按 select 分支调用并施加系数相位、`ry(select, −θ)`。矩阵语义满足 ⟨0|U†(αH)|0⟩ = c·H（测试以 α = 7、含负系数的 2×2 实例断言）。

语言不内建 eps：目标矩阵近似误差、相位求解误差与求解器收敛性不进入语言核心，alpha 也不是可以随意修改的成本标签。

`combinators.py` 提供未导出的算法层组合子：`reflect_zero`（绕 |0⟩ 反射）、`pad_signal`、`tensor`、`adjoint_be`、`lcu`（N 项 LCU，选择器宽 ⌈log2 N⌉，权重经态制备、逐项受控调用、逆制备，属性携带 `lcu_terms`）、`kronecker_sum`、`projector`、`direct_sum`、`truncated_shift`、`pauli_word`、`matrix_pauli_encoding`（显式 Pauli 分解，仅支持不超过 5 个目标位的小实例，容差 1e-12；它服务于后端描述可达性，不宣称量子加速）。

## 12. 可逆定点算术层

`arithmetic.py` 提供从定点算术到可逆电路的三条路径：RIR 模块、严格门集文本与 PySparQ 原生算子。三条路径共用同一张 Boolean 网络，没有用 Python 原函数求值冒充量子执行。

`FixedFormat(width=8, fraction=3, signed=True)` 描述二补码定点格式（宽度 2..64），`encode/decode` 在浮点与位模式之间转换。`fixed_arithmetic(kind, fmt)` 按种类生成模块，共十四种：`add, sub, neg, abs, mul, div, reciprocal, sqrt, lt, eq, select, and, or, xor`。二元算术的接口是 `a, b, out, status`（各为 bits 存储），输出按 XOR 写入 out，`status[0]` 表示定义域失效（除零、负数开方等），`status[1]` 表示越出可表示范围；这两个位是状态标志，不是误差界。舍入策略记录在属性 `rounding="toward_zero; modular_wrap"` 中。

内部表示是 `BooleanNetwork`：顺序 DAG，节点为 `not/xor/and` 原语及加法进位网络、乘法移位累加、恢复除法、逐双位恢复开方等组合；支持公共子表达式复用、经典求值（`evaluate`）与独立 JSON 往返（`payload/from_payload`）。降低到 RIR 时（`operation()`），SSA 节点打包进 64 位的 `ssa_i` 私有 bank（模块 `locals`），not 对应 `x`、xor 对应 `xor`、and 对应受控 x；计算结果 XOR 到输出后，整段前向计算被 `Adjoint` 反算，满足属性 `workspace_contract="zero_in_zero_out"`。复杂度随位宽多项式增长，不用真值表替代分解；可逆 pebbling 与最优算术综合不在当前范围。

模块名默认 `arith_<sha256[:20]>`，属性携带完整网络 JSON（`arithmetic_network`），供原生注册表识别。`BooleanCppFactory` 从同一网络生成 C++ 求值与输出异或代码，经 PySparQ 的 `compile_operator` 编译为真实自定义算子（17.3）；`arithmetic_native_registry(program)` 扫描属性并完成登记。严格门集降低（17.2）把全部算术电路降到 Toffoli/U3/CZ。

## 13. 数学函数编译器与 MIR 0.1

`compile_function` 是唯一读取 Python AST 的入口。它把受限的纯数学函数编译为可逆模块，链路是 Python 源码 → MIR 0.1 → RIR 0.3 → OriginIR-ext / PySparQ。签名：

```python
compile_function(function, *, fmt=None, inputs=None, constants=None,
                 helpers=None, output_names=None, config=None,
                 max_unroll=128, entry=None)
```

`function` 可以是函数对象或 `def` 源码字符串（多函数时用 `entry` 指定入口）；`fmt` 缺省 `FixedFormat(12, 6)`；`config` 缺省 `MathConfig()`。返回 `CompiledFunction(operation, math_ir, fmt, output_names, input_layout, output_layout)`，入口模块属性携带完整 `math_ir` JSON。

### 13.1 前端可接受子集

前端读取源码解释受限 AST，不执行待编译函数、不用 eval。接受的内容包括：数值常量与常量折叠；局部赋值与元组解包；算术、比较链、布尔运算（短路编译为 select）；条件表达式与结构化 `if`（两侧环境合并为 select 节点，返回路径必须一致）；静态有界的 `for k in range(...)`（上限 `max_unroll` = 128，循环体内禁止 return）；标量或一层元组的返回；纯 helper 调用（保留为独立模块）；白名单的 `math`/`cmath` 函数与常量。参数类型经注解或 `inputs` 映射声明为 real/complex/bool，`Index(width)`（1..64）为无符号整数输入提供较短的公开寄存器。默认参数与 `constants` 在生成期固化。

拒绝的内容包括：I/O 与对象突变、任意方法调用、动态循环、递归、异常处理、生成器、lambda、源码不可读的可调用。不支持的构造抛出 `FunctionCompileError`（`ValidationError` 子类）——它不是可以编译任意 Python 程序的工具。

### 13.2 MIR 0.1 对象模型

`MathProgram{version="0.1", entry, functions}`；每个 `MathFunction` 含唯一名字、源码标签、参数表、顺序节点表与返回索引表。参数类型为 real/complex/bool/index（index 带位宽，其余 width=0；具体字长由降低时的 FixedFormat 决定，格式必须容纳 index 的完整范围）。节点 `MathNode{op, kind, args, data}` 的 op 集合是 `input, const, neg, not, real, imag, conj, abs, add, sub, mul, div, pow, lt, eq, and, or, complex, select, intrinsic, call`；args 只能引用编号更小的节点；函数图无递归，helper 调用保留为 call 节点，多返回值的不同索引在降低时合并为一次 RIR 调用。MIR 可独立 JSON 往返（[math-ir.schema.json](schemas/math-ir.schema.json)），同一数学图在不同配置下得到不同 RIR 模块符号，可以同时组装。

### 13.3 降低与数值核

降低复用第 12 章的定点 Boolean 电路；复数降低为两个独立实数寄存器。整数常量幂用平方-乘（指数绝对值上限 128），一般实数幂按 exp(b·log a) 处理并要求正底数。非多项式实函数（sqrt/exp/log/log10/三角/双曲/反三角等十五类核）用配置区间的 Chebyshev 采样系数加 Clenshaw 递推实现：采样数由阶数决定，独立于输入位模式总数，不生成真值表。`MathConfig(degree=6, intervals=())` 控制阶数（1..32）与逐核区间覆写；配置必须覆盖每个中间数学核的输入区间，当前没有自动区间推导。复数函数以实数核与代数关系分解（如复对数 = log|z| + i·atan2，复指数用欧拉式）；atan2 与 hypot 由分解合成。

### 13.4 输出契约与状态位

入口模块满足 `|inputs, outputs, status, 0_private⟩ ↦ |inputs, outputs ⊕ F̃(inputs), status ⊕ flags(inputs), 0_private⟩`：输入端口保持不变，输出按 XOR 写入（非零初始输出合法），私有工作区复净。F̃ 是电路精确实现的有限布尔函数；它与理想数学函数 F 的接近程度不由语言证明。`status[0]` 表示定义域失效，`status[1]` 表示数学核超出配置区间或字长溢出；它们不是精度界，也不模拟 IEEE NaN/Inf 与有符号零。select 的状态等于条件与被选中分支状态的 OR——未选中分支虽然计算并反算，其定义域标志不污染最终状态。

边界声明：Python `cmath` 在分支切线上依赖有符号零，定点编码没有这一信息，当前约定零虚部取非负一侧；复杂分支切线行为与近似精度待核验。Roe face 公式（`roe_formulas.frozen_roe_face`，六场量加行列索引）由 `roe.roe_face()` 编译并供第 15 章的 QFVM 物理核调用。

## 14. QLSS 输入模型与协议

本章与 [qfvm-qlss-input-model-review.md](../manual/qfvm.md)（该路径的权威文档）共同规定量子线性系统求解器的输入契约。核心决策是四层分离：QRAM 资源（可逆数据查询）→ 经典可更新数据结构 → Oracle 操作（具体寄存器更新契约）→ Protocol（普通 Python 生成器，在 RIR 中留下模块调用与资源绑定）。

### 14.1 问题对象与谱声明

`LinearSystem(sparse=None, block=None, physical_width, physical_high_value=0, rhs_norm=None, data_assumptions=())` 是问题侧入口，只能指定一种源输入模型。`SparseSystem(access, value_format, entry_bound, rhs, spectrum, diagonal_nonnegative=False, hermitian=False)` 携带 CKS 形式的稀疏访问；`BlockSystem(encoding, rhs, spectrum)` 携带块编码输入。`block_input()` 在需要时从稀疏访问显式适配出块编码并返回适配轨迹；反方向（从任意 BE 恢复稀疏 oracle）一般没有保证，语言不提供。

`SpectralPromise(norm_upper, sigma_min_lower, evidence="caller_declared_unverified")` 记录范数上界与最小奇异值下界及证据来源。它是调用者声明，不由语言证明。库 BE 的零信号角块是 D/alpha，因此 Costa 侧看到的编码逆谱界是 `encoded_inverse_bound = alpha / sigma_min_lower`（下限 1）——它一般不同于条件数 cond(D)。Costa 配置中的 `kappa` 指的正是这个编码矩阵逆谱界。

### 14.2 协议对象与求解结果

`QLSSProtocol(name, input_model, kernel)` 的 `input_model` 取 `sparse`（CKS 路线，消费 SparseSystem）或 `block_encoding`（Costa 路线，请求 BlockSystem 并在必要时触发稀疏→BE 适配）。`solve(problem)` 按输入模型分派内核、执行能力检查、选择物理输出通道（`select_subspace`），并把全部谱声明、适配轨迹与数据假设写入模块属性，同时构造独立的矩阵范数探针。

`SolveResult(state, norm_probe, input_alpha, encoded_inverse_bound, rhs_norm, adapter_trace, kernel_status)` 统一两条路线的输出。范数恢复公式是 `‖x‖ = ‖r‖ / (alpha · sqrt(p_joint / p_solver))`，其中 `p_joint / p_solver = ‖D|x̂⟩‖² / alpha²` 由范数探针给出。Costa filtering 的成功率不能直接套用 CKS 逆算子 LCU 的归一化因子，这是范数探针独立存在的原因。已知 RHS 范数为零时入口拒绝制备归一化 RHS，要求经典侧直接处理零更新。

### 14.3 CKS 路线

`make_cks_qlss(CKSConfig(order=2))` 生成基于 Chebyshev 截断逆多项式的基础 LCU 求解器（arXiv:1511.02306 §4），没有 VTAA。它要求 Hermitian 稀疏系统；经 `real_symmetric_sparse_encoding`（14.4）把稀疏访问转为块编码后组装 LCU。实现范围明确标注为基础构造。

### 14.4 稀疏→块编码适配

`real_symmetric_sparse_encoding`（`sparse_models.py`）是当前正式适配，专门支持实 Hermitian、非负对角输入（QFVM 的 Hermitian 扩张满足这一条件）；一般稀疏矩阵不能仅靠相同位宽套用。它生成 CKS 型 T：PREP（前缀态制备、位置查询、元素 XOR、幅度旋转）之后交换 target 与 neighbor 坐标并交换两侧失败旗标，形成自伴酉扩张 T†ST。每行每列 s 个结构位置且 amax 约束条目幅值时，角块为 D/alpha 且 `alpha = s · amax`。字宽不超过 12 位时幅度转导有普通门实现，更宽时保留为显式待绑定的 transducer。见证检查（2×2 负非对角角块与 walk 三次幂投影）检验符号与归一化约定，不构成对所有矩阵与字长的证明。

### 14.5 Costa 路线

`costa_qlss(a, bprep, config)` 组装 general walk（按插值调度 s → f(s) 的受控块编码调用、RHS 零输入反射、信号位反射）多步链，末端可接相干 filtering（unary 权重制备、按 clock 位受控的 walk 幂、反制备；Dolph–Chebyshev 权重经 Laurent 递推生成，负幂用 Adjoint 表达）。`CostaConfig(steps, kappa, schedule_power, filter_degree, filter_attenuation)` 中 `kappa` 的含义见 14.1。内核状态标注为原型，解态正确性与调度误差不在本阶段认证范围。

## 15. QFVM 路径

QFVM 路径（`qfvm.py`、`qfvm_sparse.py`、`roe.py`、`flow_data.py`）实现了从 QRAM 原始流场到 QLSS 问题的完整管线，专门化范围是周期一维 Euler、三个守恒量与 frozen-Roe Jacobian（arXiv:2102.03557 的一个子集，不宣称覆盖原文全部维度、网格与边界）。

`roe_qfvm_inputs(cell_width=2, angle_width=8)` 构造输入集合：密度、动量、能量三个场数据库，几何数据库，PTheta 角表数据库，以及 RHS 制备。几何表按 (row, slot) 打包邻居、反向槽、源单元、行列分量、谱带与有效位，是 O(Ns) 的纯结构表，不含任何矩阵元素——数据库不预存按行列计算的 Jacobian 条目。

`roe_entry` 是可逆 Roe 矩阵元模块：对给定源单元取西、中、东三个单元的原始守恒量，调用两次 `roe_face`（由 `frozen_roe_face` 纯函数经第 13 章自动编译，保留熵修正参数），加中心质量项，按谱带三路选择。全部中间量在私有工作区计算并反算。算术失效或溢出时，当前有限实现把该条目总化为零——谱声明应当适用于实际量化后的矩阵，而不是只针对理想连续 Jacobian。

`qfvm_sparse_access` 产出 CKS 形式的 O_F（九个结构位置的相干换位链，补全为完整置换，逆由 adjoint 给出）与 O_A（row/column/data 接口，结构域外返回零，三分量补齐对角为 `padding_value`）。Hermitian 扩张采用 D = [[0, M], [Mᵀ, 0]]、b_D = [r, 0]，统一选择扩张标志为 1 的物理通道；物理 M 的可逆性与量化后的谱界由问题方声明。

`roe_qfvm_problem(inputs, *, spectrum, rhs_norm=None, amax, padding_value, **entry_options)` 构造 `LinearSystem`，`SpectralPromise` 是必填参数（padding 并入谱界扩展）；`data_assumptions` 记录六条数据侧假设（相干 XOR QRAM、快照跨相干调用固定、无矩阵元表、符号残差 oracle 加范数树、干净可逆 RHS 制备、原生 bank 当前重物化）。同一个问题对象可以交给 `make_cks_qlss` 或 `make_costa_qlss`——替换发生在问题与 protocol 之间并保留显式适配，不能把不同输入接口当成相同签名。替换 protocol 后应重新生成上层寄存器布局；两种完整求解器的"无缝等价替换"尚未认证。

经典侧 `RoeFlowData`（`flow_data.py`）维护原始流场、界面通量、残差符号与平方范数树：长度 2^n 向量的新制备线路做 2n 次角度 bank 查询（含反算），单点修改只重算相邻界面、三个残差单元与树路径。PySparQ 当前没有 QRAM 局部写接口，变更 bank 需要重新物化——不能算作论文假设的常数时间物理写。

## 16. 一般 QHAM 生成层

`qham` 包把自治有限多项式演化 PDE 自动推导为量子线性系统，依据是量子适配线性化 QHAM（arXiv:2411.06759v2，DOI 10.1007/s11433-024-2584-2）。数学推导的权威记录在 [qham-general-derivation.md](qham-derivation.md)，实现说明在 [qham-general-implementation.md](../manual/qham.md)。适用范围是自治一阶时间演化 PDE：任意有限分量与空间维、未知场及其空间导数的有限多项式、已知空间系数与强迫。规则之外的内容包括未知场除法、超越非线性、未消去的约束、一般隐式时间方程与时间依赖绑定——数学位置保留，但当前组装要求自治。

### 16.1 PDE 0.1 与 QCL plan 0.1

`Field`/`Known`/`Expr` 构成表达式代数（乘、数除、`d(axis, order)` 外导数、归一化）；`PolynomialPDE.from_equations(dict, axes, label)` 从方程组构造 PDE 并按 [pde.schema.json](schemas/pde.schema.json) 序列化。`ports()` 把表达式拆成线性算子 L、强迫 F 与每个非线性单项式一个多线性映射 B_tau——端口包含 PDE 系数，绑定裸乘法是不够的。

`QHAMPlan(pde, order)` 是惰性闭包对象：只存 PDE 与阶数 m，按需生成各阶同伦递推（Ui' = L·Ui − η·Σ(1+η)^(i−1−l)·C_l）、有序张量字 Y_a、物理输出块 u_sum、强迫常量分量 one，以及维数、偏移与索引映射。闭包权重条件是 weight(a) = grade·Σa + len(a) ≤ grade·m + 1（grade = max(1, D−1)），对给定 HAM 截断是精确闭包而非二次截断。二次无强迫 m = 20 的计划可以描述 2^21 个函数块而不写任何 JSON——显式导出超预算时计划仍然有效，不伪装成已展开线路。QCL plan 按 [qcl-plan.schema.json](schemas/qcl-plan.schema.json) 序列化，两层表示都不含 Python callback。

### 16.2 结构化差分端口与量子降低

`structured_fd_bindings(discretization, initial)` 生成结构化端口绑定：周期二次幂网格上的中心差分导数降低为少量可逆移位的 LCU，分量选择是输入投影加输出编码，同点乘积是坐标 XOR、对角条件与矩形收缩的组合，已知系数是对角乘法块编码。全程不物化 N^r × N^r 的端口矩阵，更不物化整个提升矩阵。

降低层（`quantum.py`）把每个端口作为矩形线性映射放入张量字：`place_port` 做槽位置换，`embed_rectangular` 同时检查输入与输出窗口并用两个失败旗标防止零填充溢出到相邻块。生成元逐行逐耦合组装为 LCU（`generator_encoding`，块与项数预算超限时报错并建议保留开放生成元）。`lifted_initial` 按对数范数稳定的分支权重（r, r, r², …, r^K，有强迫时补 1）制备提升初态——调用的是可重复可受控的制备操作，不是复制未知量子态；每次制备复净工作区。初态整体范数以 `log_initial_norm` 记录，程序不会为了成功率偷偷归一化原 PDE。

### 16.3 QODE 输入模型与两级开放

线性 QODE、Carleman 与空间离散化的具体调用签名，以及多种 given-oracle 的适配方式，见 [QPDE/QODE 实作指南](../manual/differential-equations.md)。其中区分 `make_qpde` 的三参数线性接口与 `qpde_solver` 的 model/time 接口；这些组合不增加新的 RIR 节点。

`qham_input_model(plan, bindings, eta)` 返回 `QHAMInputModel`：携带生成元块编码、提升初态制备、状态宽度、η 与初始范数，`solve(qode, time)` 调用任意 QODE 生成器并选择 `qham_physical_sum` 输出通道。`dissipative_shift()` 施加显式整体移位 G − μI（μ ≥ α_G ≥ ‖G‖）供 LCHS/CBMD 使用，保存 growth_shift·t 作为对数幅值恢复因子——它保持归一化方向但可能显著影响成功概率与成本。有限 Taylor QODE 是打通门级执行与对照的普通候选，不要求 Hermitian 或耗散前提，也不代表高效最优算法。

实现边界有两级：`QHAMBindings.declare` 逐个声明 L/F/B_tau/初态等基础端口再生成 G 与 QODE 调用，可以逐个绑定；`open_qham_input` 把整个生成元保留为未完成模块（属性携带完整 QCL plan），适合没有高效全局端口或超出显式块预算的场合。更换端口算法若改变 alpha 或辅助位宽度，应从同一数学计划重新生成上层 RIR；只有保持实例化签名与 alpha 时才适合在已有 IR 上晚绑定。稀疏输入不是从块编码自动恢复的：`qcl_row/qcl_entry` 是可审阅的经典参考，选择稀疏输入型求解器需要另行提供满足约定的稀疏访问（强迫注入可产生稠密列，行稀疏不保证 CKS 双边稀疏假设成立）。

CLI `python -m pyqecclang.applications.qham` 与 `report.export_derivation`（pde.json、qcl-plan.json、rows.json、qode-input.json、derivation.md 五件套）见第 18 章。验证边界：链式法则见证残差约 1e-16（检验代数生成的正确性，不是原 PDE 的收敛误差）；HAM 收敛、量子求解精度、时间依赖适配、IQHAM 外层重启与完整范数估计仍是待办。

## 17. 后端契约

### 17.1 OriginIR-ext 导出

`export_originir(program)` 要求程序闭合，产出 `OriginIRArtifact(text, registers, resources, workspace_qubits)`。文本结构是：每个 QRAM 一行 `QRAMDECL ram_名 aw,dw`、`QINIT 总量子位数`、`CREG 0`、按缓存去重的 `DEF` 定义序列、末尾一条对入口模块的调用。每个模块（及重复体）编译为 `DEF m_模块名_哈希 … ENDDEF`，形参是 `v_寄存器[宽]`（跳过零宽）、`pw_work[工作区宽]` 与 `pc_control[n]`（模块级控制经附加形参传递，不展开）。入口寄存器连续映射到 `q[i]`，私有工作区追加其后。

降低规则：`xor` 降为逐位 CNOT，`swap` 降为逐位 SWAP，`add_const` 降为按置位位的 X 级联（高位 X 受低位区间控制），无控制 `gphase` 用一对 U1 与 X 门精确表达全局相位。控制用逐门 `controlled_by` 后缀表达并合并去重，零值控制位用 X 共轭。`Repeat` 采用对数规模的二分策略（count 为 1 时内联，否则半量 DEF 调用两次、奇数补一次），2^40 次重复产出不到 45 个 DEF。同一模块按不同的 QRAM 资源绑定特化为不同 DEF，同绑定复用同一定义。QRAM 按全局名特化是应对 OriginIR DEF 无资源形参的既定策略；UnifiedQuantum 解析器会把 DEF 内联展开，这一限制影响它的执行阶段，不影响 pyqecclang 保存与导出的模块化结构。

### 17.2 严格门集

`export_toffoli_u3_cz(program)` 在 DEF 边界内把门集降到 {TOFFOLI, U3, CZ}（X 用 U3(π,0,π) 精确表达）：CNOT 与 SWAP 进入多控制 X 梯子，相位门族进入 U3 参数表，多控制经 Toffoli 梯子借池位（`pb_work`）完成，受控 U3 用相位分解。未知门报错。`emit --basis toffoli-u3-cz` 走同一路径。

### 17.3 PySparQ 适配与原生算子

`run_pysparq(program, memory, *, max_steps, max_states, native_registry, report)` 把闭合程序交给 PySparQ 执行。适配器要求接管前全局寄存器表为空并在结束时清理；入口寄存器改名 `pyqec_i` 并按存储类映射；QRAM 物化为稠密表（地址宽不超过 20）；跨寄存器视图与切片用临时寄存器做可逆 XOR 复制后完整反算，不依赖未知量子态的复制。每次原生门后检查稀疏基态数不超过 `max_states`。

`NativeRegistry.register(module, factory)` 按"同名模块逐字段相等"登记原生实现；`matching` 核对程序中的模块与注册版本一致，`missing` 报告可达但未覆盖的开放声明。原生执行时控制条件经真实的 split_systems/combine_systems 分区满足。`DynamicCppFactory` 惰性调用 PySparQ 的 `compile_operator`（需要 C++17 编译器），生成代码直接访问实例寄存器存储；C++ 源码不进入 RIR。原生路径跳过模块内部的算术门与工作区，但 native-only 模块不因此获得门级导出能力。report 记录 native_calls、gate_events 与 native_labels，初始 `correctness="pending"`。

### 17.4 参考执行器

`simulate(program, memory=None, *, initial=None, max_steps=1_000_000, max_states=65536)` 是纯标准库实现。初态是各寄存器取给定整数值（缺省全零）的单基矢；态以寄存器整数值元组为键的稀疏字典保存；幅度绝对值不超过 1e-15 的分量被截去；基态数超限报错。事件遍历惰性进行、不修改原始 IR：`Repeat` 只在体成本非零时迭代（空重复即使 2^63 也安全），`Adjoint` 逆序翻转，`locals` 每次调用分配唯一合成名并在退出时检查全部分支复净为零。`RegisterState.statevector(max_qubits=20)` 按寄存器顺序 LSB-first 打包。执行预算（`expanded_steps`）按指令成本模型饱和累加，超限报"执行超过展开预算"。参考执行器服务于小规模交叉验证，不提供精度保证。

`run_originir` 走真实 uniqc 模拟器：总宽度加工作区不超过 24 量子位、每 QRAM 地址加数据宽不超过 30 且地址宽不超过 20、展开步数预算 10^6（默认）。工作区量子位幅度大于 1e-12 判"未复净"。后端可以施加更严格的预算，但遇限必须报错，不能静默截断。

### 17.5 读出与重置

末端测量与重置不在 RIR 内。`ReadoutAction(kind, register)`（kind ∈ measure/reset）以宿主计划形式独立保存；`export_with_readout` 在需要动态解析器时先把已闭合模块展平，再逐量子位追加 MEASURE/RESET 并把 CREG 扩为相应数量。正常导出的 OriginIR 始终保留模块结构。

## 18. CLI 规范

入口 `pyqecclang`（`python -m pyqecclang` 同义），错误统一以退出码 2 与 `pyqecclang: 消息` 报告。

- `pyqecclang validate input.rir.json`：打印版本、模块数与开放槽计数。
- `pyqecclang canonicalize input`：输出规范 JSON（`dumps`）。
- `pyqecclang emit input [-o out.originir] [--basis {default,toffoli-u3-cz}]`：导出 OriginIR-ext 文本，或严格门集文本。
- `pyqecclang run input --memory memory.json [--backend {reference,originir,pysparq}] [--native-arithmetic] [--native-cache dir]`：执行。memory 是资源名到字序列或稀疏字典的 JSON；originir 后端按量子位索引输出幅度，reference/pysparq 按寄存器基矢输出；`--native-arithmetic` 仅限 pysparq，构造算术原生注册表并输出报告。
- `pyqecclang requirements input`：输出 `unresolved` 的 JSON 数组（槽名、范式、路径、寄存器、属性）。
- `pyqecclang bind input --bindings bindings.json -o out.rir.json`：bindings 是 `{槽名: {"program": 路径, "resources": {…}}}` 的清单，实现程序以 entry 为模块。
- `pyqecclang compile-function source.py --function 名 [--width 12] [--fraction 6] [--degree 6] [--constants json] [--inputs json] [--mir-output mir.json] -o out.rir.json`：把源文件中的纯数学函数编译为 RIR。`--inputs` 接受 `{名: 宽}`（Index）或 `{名: {"index": 宽}}` 形式的类型映射。

QHAM 子命令 `python -m pyqecclang.applications.qham [input.json] [--example {burgers,kdv,reaction,coupled,vector_burgers_2d}] [--order 2] [--eta -0.4] [--state-width 2] [--max-blocks 256] [--row physical|one|0,1] [-o 目录]`：从 PDE JSON 或内置案例生成推导五件套并输出 manifest。

## 19. 一致性分级与非保证

已验证的部分：语言核心的结构规则与错误诊断、JSON 序列化的确定性与往返保等、参考执行器与 OriginIR/PySparQ 后端在小规模实例上的复幅度对拍、全部目录案例的原生解析可达性（详见 [validation.md](../archive/validation.md)、[workboard.md](../archive/workboard.md)、[stage2-board.md](../archive/stage2-board.md) 与各阶段 validation JSON）。

待核验的部分：全部算法组装的数学正确性、数值精度、成功概率与复杂度优势。具体地：块编码的目标矩阵误差、态制备的分支切线行为、定点算术的溢出与符号边界、数学函数近似的区间覆盖、Costa/CKS 的解态正确性、QFVM 的数值求解器等价与完整经典-量子闭环、QHAM 的收敛性与范数恢复。这些以 `correctness: pending` / `not_assessed` 属性记录在模块与报告中，且不作为组装、导出与执行的前置条件。

明确不属于语言承诺的内容：QRAM 查询的物理时间复杂度（论文资源模型假设，写出 QRAMDECL 不等于证明硬件性质）、辅助位复净的静态证明、宿主读出与后选择、eps 与资源估计、以及任何"量子加速"声明。

## 20. 版本与兼容

包版本 0.7.0。从各阶段验证记录可得包版本与功能里程碑的对应：0.3.0 对应第二阶段交付（私有工作区、原生算子、算术分解、QODE/QPDE 组装），0.4.0 对应纯函数编译，0.5.0 对应 QFVM/QLSS 输入模型审查，0.6.0 对应一般 QHAM，0.7.0 对应算法自有的 Python 协议、可查询契约与工程加固；更早的 0.1.x/0.2.x 覆盖语言核心与范式阶段（开放 IR、绑定、目录案例）。包版本与 RIR 版本编号相互独立。RIR 版本当前 0.3；0.1 的 JSON 不含开放声明，0.2 不含 `locals`，读取时自动回填，写出旧版本时删去 `locals` 键。MIR、PDE 与 QCL plan 均为 0.1。

## 附录 A：公共 API 总表

顶层 `pyqecclang`（`__init__.py`）导出以下名字。

| 分组 | 名字 |
| --- | --- |
| IR 对象 | `VERSION, ValidationError, Bits, UInt, SInt, Rational, RegType, Register, Span, Ref, fuse, QRAM, Resource, Primitive, Load, Call, Repeat, Control, Adjoint, Module, Program` |
| 构造 | `Builder, Operation` |
| 验证与序列化 | `validate, dumps, loads` |
| 开放 IR | `OracleRequirement, Binding, unresolved, capabilities, bind` |
| Oracle 声明 | `declare` |
| 块编码 | `BlockEncoding, Generator, block_encoding, identity, pauli_x, zero, product, scale, linear_combination` |
| 执行 | `RegisterState, simulate, run_originir, run_pysparq` |
| 后端导出 | `OriginIRArtifact, export_originir, export_toffoli_u3_cz` |
| 原生注册 | `NativeRegistry, DynamicCppFactory` |
| 算术 | `FixedFormat, fixed_arithmetic, arithmetic_native_registry` |
| 数学函数编译 | `CompiledFunction, FunctionCompileError, Index, MathConfig, MathProgram, compile_function, lower_math_ir` |
| QLSS 契约 | `SpectralPromise, SparseSystem, BlockSystem, LinearSystem, SolveResult, QLSSProtocol` |

主要子模块入口：`pyqecclang.algorithms.oracles`（范式工厂与包装类型）、`pyqecclang.algorithms.block_encoding`（算法层组合子）、`pyqecclang.algorithms.sparse`（适配层）、`pyqecclang.algorithms.elementary/cks/costa/differential/solvers`（算法原型）、`pyqecclang.applications.legacy`（QFVM 玩具负载与 m=1 QHAM，阶段历史）、`pyqecclang.applications.catalog`（目录案例）、`pyqecclang.applications.qfvm/qfvm_sparse/roe/flow_data/sparse_models`（QFVM 路径）、`pyqecclang.applications.qham`（一般 QHAM：`Field, Known, PolynomialPDE, QHAMPlan, Grid, Discretization, structured_fd_bindings, qham_input_model, open_qham_input` 等）、`pyqecclang.infrastructure.readout`（ReadoutAction）。

## 附录 B：硬限制总表

| 限制 | 值 |
| --- | --- |
| 单寄存器/合并视图宽度 | 0..64（总量子位数与寄存器数不受限） |
| QRAM 地址宽 / 数据宽 | 各 1..64；执行侧物化另见 17.4 |
| Repeat 计数 | 0..2^63−1 |
| 结构嵌套深度 / 调用深度 | < 128 |
| 标识符 | `[A-Za-z_][A-Za-z0-9_]*` |
| 参考执行 | 基态数 ≤ 65536（默认）、幅度截断阈值 1e-15、statevector ≤ 20 位、展开步数预算 10^6（默认） |
| OriginIR 执行 | 总宽 + 工作区 ≤ 24 量子位；QRAM aw+dw ≤ 30 且 aw ≤ 20；工作区复净阈值 1e-12 |
| PySparQ | QRAM 地址宽 ≤ 20；稀疏基态数 ≤ max_states（默认 65536） |
| 幂指数（数学函数） | 绝对值 ≤ 128；循环展开 ≤ 128 步 |
| Chebyshev 阶数 | 1..32（默认 6） |
| 幅度转导门实现 | 定点字宽 ≤ 12 位，更宽保留开放槽 |
| `matrix_pauli_encoding` | ≤ 5 个目标位，容差 1e-12 |
| 公开接口 | 寄存器总宽 > 0（单个寄存器可为 0 宽） |

## 附录 C：术语表

| 术语 | 含义 |
| --- | --- |
| RIR | 寄存器级中间表示，pyqecclang 的架构中心 |
| 开放声明 / 槽 | `body=None` 的模块；带签名、待实现 |
| 绑定 | 用具体实现模块替换开放槽的链接操作 |
| QRAM 捕获 | 绑定时把实现私需资源提升为入口资源并贯通调用链 |
| 范式 | `oracle_paradigm` 的取值，接口约定的分类 |
| 块编码 / alpha | 信号位视角的算子编码及其归一化常数，存于 `be_alpha` |
| 能力 | `supports_adjoint` / `supports_controlled`，控制与伴随语境的核对项 |
| locals | 模块私有工作区寄存器，零态借入、返回前复净 |
| MIR | 数学函数图（0.1），`compile_function` 的中间表示 |
| QCL plan | QHAM 的量子适配线性化计划（0.1），PDE 加 HAM 阶数的惰性闭包 |
| SpectralPromise | 调用者声明的范数上界与最小奇异值下界 |
| 编码逆谱界 | alpha / sigma_min_lower，Costa 侧 kappa 的真实含义 |
| 范数探针 | SolveResult 中独立于求解器的矩阵范数估计通道 |
| 严格门集 | Toffoli/U3/CZ 三门字母表 |
| 原生算子 | 经 PySparQ 动态编译的模块级 C++ 实现，不改变 RIR 身份 |
| 读出计划 | 宿主层 ReadoutAction，RIR 之外的测量/重置描述 |

## 0.7 算法库约定补充

见 [算法自己的约定](../manual/contracts.md)。结构协议、requires、输入/输出报告和算法选择属于宿主算法库；对象可满足多个协议。新增算法无需修改 RIR，也不存在穷尽算法的全局类型枚举。Operation 的 unitary/state_preparation/block_encoding 方法生成保留模块边界的普通操作视图。
