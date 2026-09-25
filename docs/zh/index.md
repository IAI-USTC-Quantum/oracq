[English](../index.html) · **简体中文**

# oracq

oracq 用 Python 组织量子算法，并生成可保存、可组合的寄存器级中间表示（RIR）。算法可以依赖尚未实现的 oracle；在提供具体实现后，同一份描述可以导出为 OriginIR-ext，或交给 PySparQ 执行（导出与执行细节见[后端与导出](manual/backends.md)）。

如果你第一次使用这个项目，请先读[核心概念](manual/concepts.md)，再从教程开始。若你正在实现或审查算法，完整文档会说明接口、输入前提、组合规则和后端限制。API 参考直接来自当前源码。

```{toctree}
:maxdepth: 1
:caption: 完整文档

manual/index
reference/index
api/index
```

```{toctree}
:maxdepth: 2
:caption: 教程

tutorials/index
```

```{toctree}
:maxdepth: 2
:caption: 项目开发

development/index
```

项目目前提供从 oracle 查询、搜索、估计和 Hamiltonian 演化，到 QLSS、QODE、QFVM 与 QHAM 的组装路径。各实现的成熟度不同；请在选择算法前查看[适用范围与验证状态](manual/limits.md)。
