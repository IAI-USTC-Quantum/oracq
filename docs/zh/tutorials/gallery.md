# 运行与修改算法展示目录

[English](../../tutorials/gallery.html) · **简体中文**

展示目录把不同类别的算法放在统一的运行流程中，便于比较它们的输入、寄存器和读出方式。全部案例由 {obj}`algorithm_gallery <oracq.applications.gallery.algorithm_gallery>` 统一构造。

```bash
uv run python examples/algorithm_gallery.py
```

每个子目录包含 `closed.rir.yaml`、`modular.originir` 和 `toffoli_u3_cz.originir`。Bernstein–Vazirani 还保存绑定前的开放描述。`index.json` 记录类别、模块数量和读出说明。

## 选择要检查的结果

不要对所有算法使用同一种读出方式。BV 读取秘密字符串；QPE 和振幅估计读取 phase；Hadamard test 估计 probe 的 Z 期望；重复码恢复则检查逻辑态，并保留 syndrome。

用以下代码查找案例：

```{testcode}
from oracq.applications.gallery import algorithm_gallery

# 构造全部展示案例，并按名字建索引方便查询。
cases = {case.name: case for case in algorithm_gallery()}
# 打印当前展示库覆盖的全部案例名。
print("\n".join(sorted(cases)))
# 确认 QAOA 最大割案例存在（展示库覆盖变分类）。
assert "qaoa_maxcut" in cases
# 每个案例带族标签：序数查找属于 number_theory 族。
assert cases["order_finding"].family == "number_theory"
# 展示库共 22 个案例，增删案例时这行断言会提醒同步更新。
assert len(cases) == 22
```

```{testoutput}
amplitude_amplification
amplitude_estimation
ansatz
bernstein_vazirani
cycle_walk
deutsch_jozsa
fourier_add
grover
hadamard_imag
hadamard_real
hamiltonian_trotter
modular_multiply
order_finding
phase_estimation
qaoa_maxcut
qft
quantum_counting
repetition_bit
repetition_phase
simon
swap_test
vqe_pauli_measurement
```

打印出的 22 个名字就是展示库当前覆盖的全部案例，横跨搜索、估计、算术、行走、纠错和变分类；增删案例时这行输出要与 `len(cases) == 22` 断言一起同步更新。

需要添加应用例子时，可以在自己的脚本里调用对应算法文件。若要扩展公开展示目录，则为 {obj}`GalleryCase <oracq.applications.gallery.GalleryCase>` 提供操作和具体读出说明，并增加独立的数学见证。

## 真实后端对拍

```bash
PYTHONPATH=src /path/to/backend/python examples/algorithm_gallery.py --native
```

脚本逐项比较参考执行器、PySparQ 和 OriginIR 的复幅度，超出容差立即退出。它验证同一线路在各后端中的语义一致；某个参数是否适合真实问题，仍需要应用层分析。

## 相关页面

- 手册：[算法手册](../manual/algorithms/index.md)（展示库覆盖的多数算法都有词条页）、[后端与导出](../manual/backends.md)
- API 参考：[算法展示目录](../api/applications/gallery.rst)
- 继续教程：[搜索一个元素，并估计成功概率](search-and-estimation.md)、[提供自己的 Hamiltonian 分解](hamiltonian.md)
