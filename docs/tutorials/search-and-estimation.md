# 搜索一个元素，并估计成功概率

Grover 搜索使用相位 oracle 标记好状态。下面在四个基态中标记 `3`，一次迭代后读取该状态。

```{testcode}
from pyqecclang.algorithms.oracles import phase_marks
from pyqecclang.algorithms.search import grover
from pyqecclang import simulate

result = grover(phase_marks(2, [3]), 2, iterations=1)
state = simulate(result.operation.program())
assert abs(state.amplitudes[(3, 0)] - 1) < 1e-12
```

输出含 target 和 signal。这里 signal 为空接口，结果字为零；更复杂的输入制备可能有工作寄存器。

## 估计概率

标准振幅估计对 Grover iterate 做相位估计。它接收初态制备和好状态集合，输出 phase 寄存器。下面的初态在 `0` 和 `1` 上均匀分布，故好状态 `1` 的概率为 `1/2`。

```{testcode}
from pyqecclang.algorithms.oracles import uniform_state
from pyqecclang.algorithms.estimation import amplitude_estimation, amplitude_from_phase

operation = amplitude_estimation(uniform_state(1), [1], precision=3)
state = simulate(operation.program())
phases = {key[2] for key in state.amplitudes}
assert phases == {2, 6}
assert all(abs(amplitude_from_phase(value, 3) - 0.5) < 1e-12 for value in phases)
```

两个相位都对应同一个概率。一般情况下，有限精度的相位读出会产生近似估计；多次采样和统计处理在宿主侧进行。
