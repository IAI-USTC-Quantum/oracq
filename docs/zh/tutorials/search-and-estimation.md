# 搜索一个元素，并估计成功概率

<a href="../../tutorials/search-and-estimation.html">English</a> · **简体中文**

[Grover 搜索](../manual/algorithms/grover.md)使用相位 oracle 标记好状态。下面在四个基态中标记 `3`，一次迭代后读取该状态。

```{testcode}
from oracq.algorithms.input_model.oracles import phase_marks
from oracq.algorithms.common.search import grover
from oracq import simulate

result = grover(phase_marks(2, [3]), 2, iterations=1)
state = simulate(result.operation.program())
print(state.amplitudes)
assert abs(state.amplitudes[(3, 0)] - 1) < 1e-12
```

```{testoutput}
{(3, 0): (0.9999999999999996-1.8369701987210287e-16j)}
```

{obj}`grover <oracq.algorithms.common.search.grover>` 接收 {obj}`phase_marks <oracq.algorithms.input_model.oracles.phase_marks>` 声明的相位 oracle。打印出的字典只有一个基态：键 `(3, 0)` 表示 target 读数为 `3`、signal 为零，幅度的模方约等于 `1`（虚部是浮点噪声）。输出含 target 和 signal。这里 signal 为空接口，结果字为零；更复杂的输入制备可能有工作寄存器。

## 估计概率

标准[振幅估计](../manual/algorithms/qae.md)对 Grover iterate 做相位估计，入口是 {obj}`amplitude_estimation <oracq.algorithms.common.estimation.amplitude_estimation>`。它接收初态制备和好状态集合，输出 phase 寄存器。下面的初态在 `0` 和 `1` 上均匀分布，故好状态 `1` 的概率为 `1/2`。

```{testcode}
from oracq.algorithms.input_model.oracles import uniform_state
from oracq.algorithms.common.estimation import amplitude_estimation, amplitude_from_phase

operation = amplitude_estimation(uniform_state(1), [1], precision=3)
state = simulate(operation.program())
phases = {key[2] for key in state.amplitudes}
# 解码值与 1/2 只差浮点末位，不同平台的三角函数库在末位上可能不同；
# 打印前舍入到 12 位，输出在各平台逐字一致。
print({value: round(amplitude_from_phase(value, 3), 12) for value in sorted(phases)})
assert phases == {2, 6}
assert all(abs(amplitude_from_phase(value, 3) - 0.5) < 1e-12 for value in phases)
```

```{testoutput}
{2: 0.5, 6: 0.5}
```

初态由 {obj}`uniform_state <oracq.algorithms.input_model.oracles.uniform_state>` 制备。打印出的字典把两个相位字经 {obj}`amplitude_from_phase <oracq.algorithms.common.estimation.amplitude_from_phase>` 解码回概率：`2` 和 `6` 都给出约 `1/2`，未舍入的解码值与 `1/2` 的差异在浮点末位，来自有限精度读出的舍入。一般情况下，有限精度的相位读出会产生近似估计；多次采样和统计处理在宿主侧进行。

## 相关页面

- 算法页：[Grover 搜索](../manual/algorithms/grover.md)、[振幅放大](../manual/algorithms/amplitude-amplification.md)、[振幅估计](../manual/algorithms/qae.md)、[相位估计](../manual/algorithms/qpe.md)
- API 参考：[搜索与振幅放大](../api/algorithms/common/search.rst)、[相位、振幅与重叠估计](../api/algorithms/common/estimation.rst)、[Oracle 声明与实现](../api/algorithms/input_model/oracles.rst)
- 继续教程：[运行与修改算法展示目录](gallery.md)（覆盖本例两个算法的更多案例）
