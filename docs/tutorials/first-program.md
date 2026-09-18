# 第一个寄存器程序

这篇教程生成一个 Bell 态。你会用到一个寄存器、两条门操作和参考执行器。完成后，再把同一程序保存为 RIR 和 OriginIR-ext。

## 定义操作

`Bits(2)` 定义一个两位寄存器。我们先让低位进入叠加态，再以它为控制，对高位做 XOR。

```{testcode}
from pyqecclang import Bits, Builder, simulate

# 声明一个名为 bell_pair 的模块，公开接口是一个名为 pair 的两位 bits 寄存器。
b = Builder("bell_pair", {"pair": Bits(2)})
# 对 pair 的第 0 位（最低位）广播 H 门：|0> -> (|0>+|1>)/sqrt(2)。
b.h(b["pair"][0])
# 以第 0 位为源、第 1 位为目标做 XOR（逐位 CNOT）：得到 (|00>+|11>)/sqrt(2)。
b.xor(b["pair"][0], b["pair"][1])
# 结束构建，得到不可变的 Operation（含模块与指令体）。
operation = b.finish()
# 取出其中的 RIR Program：entry 指向 bell_pair，modules 里只有这一个模块。
program = operation.program()
# 用无依赖的参考执行器模拟，返回稀疏振幅。
state = simulate(program)

# 振幅键是入口寄存器按声明顺序拼接出的整数元组；(0,) 即 00，(3,) 即 11。
assert set(state.amplitudes) == {(0,), (3,)}
# 两个基态的幅度都应为 1/sqrt(2)。
assert abs(state.amplitudes[(0,)] - 2**-0.5) < 1e-12
assert abs(state.amplitudes[(3,)] - 2**-0.5) < 1e-12
```

结果用入口寄存器的整数值组成元组作为键。`0` 表示 `00`，`3` 表示 `11`；它们的幅度都是 `1/√2`。

## 保存和导出

```{testcode}
from pyqecclang import dumps, loads, export_originir

# 序列化为规范 JSON：按键排序、两空格缩进，模块与指令结构原样保留。
text = dumps(program)
# 反序列化回 Program 对象；解码是严格的，并会重新跑语义校验。
restored = loads(text)
# 往返应逐字节等价：寄存器名、视图和模块结构都不丢失。
assert restored == program
# 导出为模块化 OriginIR-ext 文本。
artifact = export_originir(restored)
# 文本里应出现 DEF 模块定义：导出不会内联调用、不会展开成平坦门列表。
assert "DEF" in artifact.text
```

JSON 保存的是寄存器和模块结构。OriginIR-ext 输出包含模块定义；导出不会先把所有调用复制成一个平坦门列表。

`simulate` 适合这种小规模检查。需要实际后端时，仍使用同一个 `Program`，改为调用 `run_pysparq` 或 `run_originir`。
