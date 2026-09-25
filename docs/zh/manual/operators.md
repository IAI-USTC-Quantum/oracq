# Oracle 与算子表示

[English](../../manual/operators.html) · **简体中文**

选择输入表示时，应先问上层算法需要怎样访问数据。一个矩阵的 BE、稀疏位置查询和数值 XOR 查询，提供的是不同能力。

| 表示 | 可读取的信息 | 主要用途 |
|---|---|---|
| {obj}`BlockEncoding <oracq.algorithms.input_model.operators.BlockEncoding>` | width、signal_qubits、alpha、capabilities | 线性算子组合、QLSS 和演化 |
| {obj}`StatePreparation <oracq.algorithms.input_model.oracles.StatePreparation>` | width、work_width、零输入和复净约定 | 初态、右端态与反射 |
| {obj}`XorDatabase <oracq.algorithms.input_model.oracles.XorDatabase>` | address_width、data_width | 可逆数据访问 |
| {obj}`SparseAccess <oracq.algorithms.input_model.oracles.SparseAccess>` | 位置/元素两个操作、sparsity、value_width | CKS 输入及稀疏适配 |
| {obj}`Operation <oracq.infrastructure.builder.Operation>` | 公开寄存器、资源与模块依赖 | 完整量子操作 |

同一个对象可以提供多种视图。例如完整 unitary `U` 可以制备 `U|0>`，也可以作为 alpha=1 的 BE 参加 LCU。算法通过方法和 Python 协议检查这些能力，约定详见[算法自己的约定](contracts.md)。

## 归一化常数

若零信号角块满足 `⟨0|U_A|0⟩=A/alpha_A`，那么 alpha 是组合所必需的量。乘积的 alpha 相乘；LCU 的 alpha 为各项 `abs(coefficient)*alpha` 之和。

这个常数不会让执行器额外缩放量子态。它描述实际线路的成功角块。修改 alpha 时必须重新生成依赖它的分支权重和角度。

## 数据字与幅度

XOR database 返回位模式，不直接提供幅度访问。{obj}`diagonal_block_encoding <oracq.algorithms.input_model.oracles.diagonal_block_encoding>` 的输入是角度字，所编码的对角元为 `alpha*cos(angle/2)`。一般矩阵元需要相应的数值到幅度转换。

稀疏输入还需要明确位置 oracle 的可逆形式。当前 CKS 位置接口原地置换 index；保留 index 的 XOR 表不能直接替代它。

## 开放 oracle

使用 {obj}`abstract_block_encoding <oracq.algorithms.input_model.oracles.abstract_block_encoding>`、{obj}`abstract_state_prep <oracq.algorithms.input_model.oracles.abstract_state_prep>`、{obj}`abstract_database <oracq.algorithms.input_model.oracles.abstract_database>` 或 {obj}`abstract_sparse_access <oracq.algorithms.input_model.oracles.abstract_sparse_access>` 可以先写算法，再提供实现。开放状态只影响是否能执行或导出，结构检查仍然有效。

完整的绑定例子见[教程：替换 oracle](../tutorials/oracle-binding.md)，API 见[oracle 目录](../api/algorithms/input_model/oracles.rst) 与 [BE 组合](../api/algorithms/input_model/block_encoding.rst)。
