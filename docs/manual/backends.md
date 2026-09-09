# 导出与执行后端

pyqecclang 将导出和执行分开。`export_originir` 与 `export_toffoli_u3_cz` 只生成描述，不需要安装量子模拟器。

## OriginIR-ext

导出器保留 `DEF`、模块调用和 `QRAMDECL`。寄存器操作在导出时降低为具体门；Repeat 使用可复用的辅助定义表达。QRAM 内存以外部资源形式传入，不嵌入 RIR。

`export_toffoli_u3_cz` 将普通门降低到 Toffoli、U3 和 CZ，QRAM 仍保留为独立资源指令。它没有展开物理 QRAM 的器件网络。

`run_originir` 使用 UnifiedQuantum 执行。该后端的解析器会展开 DEF，因此执行前有展开预算；这不改变保存的模块化 RIR。

## PySparQ

`run_pysparq` 按寄存器事件执行模块。`NativeRegistry` 可以在模块边界调用自定义原生算子，适合大型可逆算术。原生实现与门级主体分别管理：只有原生实现的开放 oracle 可以用于相应模拟，但不会因此获得门级导出能力。

动态算子需要匹配的 PySparQ ABI 和 C++17 编译器。具体接口行为与已审阅的版本见[后端兼容说明](../reference/backend-compatibility.md)。

## 参考执行器

`simulate` 是只依赖标准库的小型寄存器执行器。它便于检查位序、相位、XOR 语义和辅助位状态。状态数量、步数及小幅度截断都有限制，不能把它的可运行规模视为真实硬件资源估计。

## 读出

RIR 核心不含测量和重置。最终测量、后选择和统计处理在宿主层完成；需要动态 OriginIR 时，可以显式使用读出适配。带 syndrome 的纠错恢复电路会保留错误信息，重新使用这些寄存器前需要适当的宿主处理。
