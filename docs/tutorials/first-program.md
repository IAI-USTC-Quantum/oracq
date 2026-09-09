# 第一个寄存器程序

这篇教程生成一个 Bell 态。你会用到一个寄存器、两条门操作和参考执行器。完成后，再把同一程序保存为 RIR 和 OriginIR-ext。

## 定义操作

`Bits(2)` 定义一个两位寄存器。我们先让低位进入叠加态，再以它为控制，对高位做 XOR。

```{testcode}
from pyqecclang import Bits, Builder, simulate

b = Builder("bell_pair", {"pair": Bits(2)})
b.h(b["pair"][0])
b.xor(b["pair"][0], b["pair"][1])
operation = b.finish()
program = operation.program()
state = simulate(program)

assert set(state.amplitudes) == {(0,), (3,)}
assert abs(state.amplitudes[(0,)] - 2**-0.5) < 1e-12
assert abs(state.amplitudes[(3,)] - 2**-0.5) < 1e-12
```

结果用入口寄存器的整数值组成元组作为键。`0` 表示 `00`，`3` 表示 `11`；它们的幅度都是 `1/√2`。

## 保存和导出

```{testcode}
from pyqecclang import dumps, loads, export_originir

text = dumps(program)
restored = loads(text)
assert restored == program
artifact = export_originir(restored)
assert "DEF" in artifact.text
```

JSON 保存的是寄存器和模块结构。OriginIR-ext 输出包含模块定义；导出不会先把所有调用复制成一个平坦门列表。

`simulate` 适合这种小规模检查。需要实际后端时，仍使用同一个 `Program`，改为调用 `run_pysparq` 或 `run_originir`。
