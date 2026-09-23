# 搜索一个元素，并估计成功概率

Grover 搜索使用相位 oracle 标记好状态。下面在四个基态中标记 `3`，一次迭代后读取该状态。

```{testcode}
from pyqecclang.algorithms.input_model.oracles import phase_marks
from pyqecclang.algorithms.common.search import grover
from pyqecclang import simulate

result = grover(phase_marks(2, [3]), 2, iterations=1)
state = simulate(result.operation.program())
print(state.amplitudes)
assert abs(state.amplitudes[(3, 0)] - 1) < 1e-12
```

```{testoutput}
{(3, 0): (0.9999999999999996-1.8369701987210287e-16j)}
```

打印出的字典只有一个基态：键 `(3, 0)` 表示 target 读数为 `3`、signal 为零，幅度的模方约等于 `1`（虚部是浮点噪声）。输出含 target 和 signal。这里 signal 为空接口，结果字为零；更复杂的输入制备可能有工作寄存器。

## 估计概率

标准振幅估计对 Grover iterate 做相位估计。它接收初态制备和好状态集合，输出 phase 寄存器。下面的初态在 `0` 和 `1` 上均匀分布，故好状态 `1` 的概率为 `1/2`。

```{testcode}
from pyqecclang.algorithms.input_model.oracles import uniform_state
from pyqecclang.algorithms.common.estimation import amplitude_estimation, amplitude_from_phase

operation = amplitude_estimation(uniform_state(1), [1], precision=3)
state = simulate(operation.program())
phases = {key[2] for key in state.amplitudes}
print({value: amplitude_from_phase(value, 3) for value in sorted(phases)})
assert phases == {2, 6}
assert all(abs(amplitude_from_phase(value, 3) - 0.5) < 1e-12 for value in phases)
```

```{testoutput}
{2: 0.4999999999999999, 6: 0.5000000000000001}
```

打印出的字典把两个相位字解码回概率：`2` 和 `6` 都给出约 `1/2`，末位的偏差是有限精度读出的舍入。一般情况下，有限精度的相位读出会产生近似估计；多次采样和统计处理在宿主侧进行。
