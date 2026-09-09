# pyqecclang

pyqecclang 是基于 Python 的量子操作生成语言。生成函数返回可组合的操作对象；它们最终形成 **模块化、register-level 的 RIR**，而不是必须展开的量子位门列表。

首个导出后端是 UnifiedQuantum 的 **OriginIR-ext**。导出结果保留 `DEF`、嵌套调用和 `QRAMDECL`。仓库另外提供独立的寄存器参考执行器和可选的 PySparQ 执行适配器，用于交叉验证。

本仓库独立于 QECC.Lang，不继承其文本语法或扁平 CLIR 约束。当前版本为 0.6.0。

新增一般 QHAM 自动生成：有限多项式 PDE → 任意阶 HAM → 惰性张量闭包 → 模块化 QODE 输入。先看 [数学推导](docs/qham-general-derivation.md)，再看 [实现与运行说明](docs/qham-general-implementation.md)。支持强迫、耦合分量和高次项；数学恒等式与小型真实后端已验证，收敛与量子精度仍单独保留。

QFVM / QLSS 的输入模型已重新审查并修正。请看独立文档 [QFVM 中替换 QLSS](docs/qfvm-qlss-input-model-review.md)：QFVM 提供稀疏 oracle，CKS 和 Costa 通过显式 protocol/适配器接入；alpha、谱声明、物理输出通道和范数探针均独立记录。算法数值等价和完整求解正确性仍待核验。

新增 [普通数学函数自动量子编译](docs/function-compiler.md)：compile_function 将受限的纯 Python math/cmath 函数编译为可逆 XOR 模块，支持实数、复数、helper 调用和分支；自动生成工作区与反算。QFVM 的 Roe face 已改用普通 Python 公式自动编译。MIR 0.1 和 RIR 0.3 都可独立序列化，近似精度仍由生成配置决定且待核验。

第二阶段已实现 PySparQ 自定义算子、14 类算术的 Toffoli/U3/CZ 分解、QRAM + Roe 算术 QFVM，以及 Carleman / Schrödingerization / LCHS / CBMD 的开放 QODE/QPDE 组装。请先看 [实施说明与使用方式](docs/stage2-implementation.md)、[工作面板](docs/stage2-board.md) 和 [验收记录](docs/stage2-validation.json)。算法正确性仍待核验。

## 本轮范式面板

- [实现工作面板](docs/workboard.md)列出 P0–P8 计划和每个案例的开放、部分绑定、闭合 IR 与 OriginIR 产物。
- [开放 IR 0.2](docs/open-ir.md)定义 null 主体、分批绑定、缺口报告和资源提升。
- [Oracle 范式目录](docs/oracle-paradigms.md)说明 BE、XOR database、state-prep isometry 和 CKS sparse。
- [旧用例覆盖矩阵](docs/coverage.md)逐项对应 61 个正例与 16 个负例。
- [应用组装说明](docs/application-assemblies.md)说明 Costa/filter、QFVM 和 QHAM 的当前范围。

本阶段评价语言范式和组装可用性，不认证算法精度、成功概率或流体数值结果。未完成的 oracle 是有效 IR，不会用空线路代替。

```bash
uv run python tools/build_catalog.py
uv run pyqecclang requirements out/catalog/costa_qram/open.rir.json
uv run pyqecclang bind out/catalog/costa_qram/open.rir.json --bindings out/catalog/costa_qram/bindings.json -o out/costa-bound.rir.json
uv run pyqecclang emit out/costa-bound.rir.json -o out/costa.originir
```

如需使用真实后端验收全部描述，请在已安装 uniqc 的环境运行 tools/build_catalog.py --native-parse。读出/重置计划独立保存，只有显式宿主读出适配会使用下游要求的末端展平。

## 安装和验证

语言核心没有第三方运行依赖，要求 Python 3.11 或更新版本。

```bash
uv sync --all-extras
uv run python -m unittest discover -s tests/core -v
uv run python -m unittest discover -s tests/schema -v
uv run ruff check src tests examples tools
uv run python examples/qram_modules.py
```

也可以不安装包，直接用标准库运行：

```bash
PYTHONPATH=src python3 -m unittest discover -s tests/core -v
PYTHONPATH=src python3 examples/qram_modules.py
```

真实后端集成测试要求解释器已安装 `uniqc` 与 `pysparq`。请使用已有后端环境运行：

```bash
PYTHONPATH=src /path/to/backend/python -m unittest discover -s tests/integration -v
```

该集成测试不以缺少后端为理由跳过。测试环境必须真正提供两个后端。本工作区审阅和验证的后端提交见 [backend-revisions.json](backend-revisions.json)。

## 一个模块化 QRAM 程序

下面的 Python 函数在生成阶段执行。它没有执行量子查询，也没有把生成函数保存在 IR 中。

```python
from pyqecclang import Builder, QRAM, UInt, dumps, export_originir, simulate

def make_lookup():
    b = Builder("lookup",
                {"address": UInt(2), "data": UInt(3)},
                {"table": QRAM(2, 3)})
    b.qram("table", b["address"], b["data"])
    return b.finish()

lookup = make_lookup()

b = Builder("main",
            {"address": UInt(2), "data": UInt(3)},
            {"values": QRAM(2, 3)})
b.h(b["address"])
b.x(b["data"][0])
b.call(lookup, address=b["address"], data=b["data"],
       resources={"table": "values"})
program = b.finish().program()

rir_json = dumps(program)
artifact = export_originir(program)
state = simulate(program, {"values": [1, 2, 4, 7]})

print(artifact.text)
print(state.amplitudes)
```

结果中的寄存器元组按入口签名的顺序排列。此例产生 `(0,0)`、`(1,3)`、`(2,5)`、`(3,6)` 四个分量，每个幅度为 1/2。数据目标初始为 1，查询采用 XOR，不是覆盖赋值。

RIR 中一个 `h(address)` 仍然是一条作用于整个寄存器的指令。只有导出到 OriginIR-ext 时，它才降低为逐比特 H 门。模块调用在这一步仍然保留。

## 模块、视图和控制

```python
from pyqecclang import Builder, Bits, UInt, fuse

b = Builder("arithmetic", {"word": UInt(8), "flag": Bits(1)})
with b.control(b["flag"], 0):
    with b.repeat(1000):
        b.add_const(b["word"], 7)

with b.adjoint():
    b.rz(b["word"][:2], 0.3)

operation = b.finish()
```

`Repeat`、`Control`、`Adjoint` 都是显式 IR 节点。大重复导出为对数数量的复用 DEF 定义。执行器按需遍历调用，不修改原始 IR。

寄存器存储模型对齐 PySparQ：

- 每个整数寄存器或合并视图的宽度为 0..64。零宽度只表达空接口，不分配原生寄存器。
- 寄存器数量及总量子位数不受 64 限制。
- 支持 `Bits`、`UInt`、`SInt` 和 `Rational` 的存储解释。
- 下标零表示最低位。连续切片与 `fuse` 是零成本的逻辑视图。
- 视图不能重叠，调用实参不能相互别名，受保护控制位不能被调用修改。
- QRAM 具有独立的地址、数据位宽和模块资源参数。

## 块编码生成器

```python
from pyqecclang import identity, pauli_x, product, linear_combination, scale

a = scale(2, identity(1))
b = scale(3, pauli_x(1))
c = linear_combination(-0.5, a, 2, b)
assert c.alpha == 7

d = product(a, b)
assert d.alpha == 6
```

`alpha` 保存在模块属性 `be_alpha` 中，经过 JSON 往返仍然存在。乘积保留独立信号空间；加权和按归一化常数选择分支幅度，并保留复系数的相对相位。语言不内建 `eps`，不保证目标矩阵近似误差或求解器收敛性。

Python 的普通函数、闭包或可调用对象即可作为 protocol 的实现。生成结果必须是可验证的 `Operation`，不能把任意 Python callback 当作未定义的 IR 指令。本版还提供 Costa general walk 与相干 filtering、QODE/QPDE、QFVM 和 m=1 QHAM 的组装原型。它们可以保留开放 oracle，也提供普通 gate/QRAM 小绑定；数值正确性未在本阶段认证。

## CLI

```bash
pyqecclang validate examples/qram.rir.json
pyqecclang emit examples/qram.rir.json -o out.originir
pyqecclang run examples/qram.rir.json --memory examples/memory.json
pyqecclang run examples/qram.rir.json --memory examples/memory.json --backend originir
pyqecclang run examples/qram.rir.json --memory examples/memory.json --backend pysparq
```

编译导出不需要安装任何量子模拟器。执行后端在调用时才导入其依赖。

## 文档与边界

先读两份总览文档：

- [pyqecclang 完整规范](docs/pyqecclang-spec.md)统一描述语言核心、四层生成性表示、算法组装层与后端契约。
- [pyqecclang 指南与案例集](docs/pyqecclang-guide.md)用九个可运行案例（QRAM 模块、开放绑定、块编码、Grover、定点算术、数学函数编译、QFVM 替换 QLSS、一般 QHAM）串起全部用法。

专题文档：

- [RIR v0.1 规范](docs/rir-spec.md) 定义对象模型、JSON 编码和每条指令的语义。
- [JSON Schema](docs/rir.schema.json) 提供结构校验，跨节点规则由 `validate` 检查。
- [后端审阅与映射](docs/backend-review.md) 记录实际 OriginIR-ext 和 PySparQ 行为。
- [架构](docs/architecture.md) 解释生成阶段、模块依赖和后端分工。

0.3 保留开放 oracle，并新增模块私有 locals、PySparQ 原生算子注册以及模块内的算术分解。私有工作区由后端分配并在模拟返回时检查复净；语言不证明复净或算法精度。测量反馈和符号形状尚未纳入，读出与后选择属于宿主层。

当前 UnifiedQuantum 解析器会展开 DEF。该限制影响它的执行阶段，不影响 pyqecclang 保存和导出的模块化结构。适配器在执行前检查展开预算，避免把巨大模块图意外交给扁平解析器。

本地参考执行器会截去幅度绝对值小于等于 `1e-15` 的分量，并限制稀疏态数量；它服务于小规模验证，不提供精度保证。
