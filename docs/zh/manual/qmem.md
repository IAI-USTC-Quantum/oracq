# QRAM 数据结构：指针、偏移与随机写

[English](../../manual/qmem.html) · **简体中文**

{obj}`QMem <oracq.infrastructure.qmem.QMem>` 把 QRAM 资源抽象成 C 风格的数组访问：基地址、常量与量子偏移、多维视图，以及随机读写。所有寻址都是 Python 生成阶段的糖衣——落到 RIR 里的只有寄存器算术、{obj}`Load <oracq.infrastructure.ir.Load>` 和 {obj}`Store <oracq.infrastructure.ir.Store>`（见[RIR 规范](../reference/rir.md) 3.2 节）。API 见 [QRAM 指针式读写](../api/infrastructure/qmem.rst)。

## 指针与读

{obj}`QMem <oracq.infrastructure.qmem.QMem>` 绑定构造器里已声明的 QRAM 资源。`ptr()` 返回指向 0 号单元的指针，整数加减产生常量偏移；传入寄存器视图则得到持有量子地址的指针。解引用 `load` 是现有 XOR-Load：`|addr⟩|d⟩ ↦ |addr⟩|d ⊕ M[addr]⟩`，叠加地址天然支持。

```{doctest}
>>> from oracq import Builder, QRAM, QMem, UInt, simulate
>>> b = Builder("demo", {"idx": UInt(2), "out": UInt(4)}, {"rom": QRAM(2, 4)})
>>> mem = QMem(b, "rom")
>>> (mem.ptr() + 1).load(b["out"])
>>> state = simulate(b.finish().program(), {"rom": [10, 11, 12, 13]})
>>> sorted(state.amplitudes.items())
[((0, 11), (1+0j))]
```

量子偏移（`p + b["idx"]`）和多维展平需要寄存器加法：`QMem` 自动合成逐位进位加法器与移位拼接，地址临时寄存器在解引用后由整体 Adjoint 复净，模块保持酉。地址算术的门成本如实计入[资源估计](resource-estimation.md)。当地址表达式恰为单个全宽寄存器且无常量分量时，不引入任何寻址算术。

## 多维视图

`shape` 声明 row-major 展平；下标可以是整数、寄存器视图或切片，量子下标的宽度必须满足 2^width 不超过该维长度。切片返回重新定址的子数组视图。

```{doctest}
>>> b = Builder("grid", {"row": UInt(2), "out": UInt(4)}, {"rom": QRAM(4, 4)})
>>> grid = QMem(b, "rom", shape=(4, 4))
>>> grid[b["row"], 2].load(b["out"])
>>> state = simulate(b.finish().program(), {"rom": list(range(16))}, initial={"row": 3})
>>> sorted(state.amplitudes.items())
[((3, 14), (1+0j))]
```

## 随机写

`store` 是 RIR 的 `Store` 指令：`M[addr] := data`。存储单元按经典单元建模，随机写不计入门成本（资源估计中单独计入 `qram_writes`）。执行时地址与数据寄存器必须处于确定基矢；叠加地址下的写没有线性语义，执行器会直接报错。结构上 `Store` 不能出现在 {obj}`Control <oracq.infrastructure.ir.Control>` 或 {obj}`Adjoint <oracq.infrastructure.ir.Adjoint>` 体内，含 `Store` 的模块不具备受控与伴随能力。

```{doctest}
>>> b = Builder("write", {"addr": UInt(2), "val": UInt(4), "out": UInt(4)}, {"ram": QRAM(2, 4)})
>>> mem = QMem(b, "ram")
>>> mem[b["addr"]].store(b["val"])
>>> mem.ptr(3).load(b["out"])
>>> state = simulate(b.finish().program(), {"ram": [0, 0, 0, 0]}, initial={"addr": 1, "val": 9})
>>> sorted(state.amplitudes.items())
[((1, 9, 0), (1+0j))]
```

后续 `Load` 读到演化后的内存：上例中单元 1 已写入 9，单元 3 仍为 0。[OriginIR-ext](backends.md#originir-ext) 导出把 `Store` 降低为 `QRAMWRITE` 扩展行；文本执行器（UnifiedQuantum、PySparQ）暂不接受运行期写，遇到含 `Store` 的程序会在执行入口报错，文本导出不受影响。

## 能力一览

| 能力 | 形式 | 降低产物 |
|---|---|---|
| 常量偏移 | `mem.ptr(2) + 3`、`p - 1` | `add_const`（模 2^地址宽度） |
| 量子指针 | `mem.ptr(b["addr"])` | 直接寻址或移位拼接 |
| 量子偏移 | `p + b["idx"]` | 与指针区间重叠时合成进位加法器（Toffoli/CNOT） |
| 多维索引 | `grid[b["row"], 3]` | 移位拼接；各维占据不相交位段时不需加法器 |
| 子数组视图 | `grid[2:]` | 常量重定址，无指令 |

完整语义测试见 `tests/core/test_qmem.py`。
