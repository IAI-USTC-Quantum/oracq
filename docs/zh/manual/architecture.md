# 源码结构

[English](../../manual/architecture.html) · **简体中文**

源码分为基础设施、算法库和领域应用。目录名对应代码的职责，算法实现按类别放在独立文件中。

```text
src/oracq/
├── infrastructure/       RIR、Builder、验证、序列化、执行和后端
│   ├── backends/         OriginIR-ext、严格门集和 PySparQ 适配
│   └── mathfunc/         纯数学函数前端、MIR 与可逆 lowering
├── algorithms/           量子算法、oracle 接口和组合工具
│   ├── fourier.py        Fourier 变换与 Fourier 加法
│   ├── search.py         Grover 和振幅放大
│   ├── estimation.py     QPE、振幅估计、Hadamard/Swap test
│   ├── hamiltonian.py    Hamiltonian 演化与可替换实现
│   ├── qlss.py           线性系统契约、Costa 和 CKS
│   ├── lchs.py           LCHS
│   ├── schrodingerization.py
│   ├── carleman.py       Carleman 线性化
│   └── ...
└── applications/         QFVM、Roe 数据、QHAM 数学支持和案例目录
```

## 基础设施

`infrastructure` 负责描述的结构与执行机制，不根据算法名字选择线路。数学函数前端可以调用算法库的算术生成器；{obj}`Operation <oracq.infrastructure.builder.Operation>` 的便捷视图方法也会在需要时调用算法适配器。这些调用发生在 Python 生成阶段，不增加 RIR 指令种类。

可选后端在执行入口导入（见[导出与执行后端](backends.md)）。生成与序列化 RIR 不要求安装 PySparQ 或 UnifiedQuantum。

## 算法库

`algorithms` 中每个类别有自己的文件和输入约定。共享的 [`interfaces.py`](../api/algorithms/input_model/interfaces.rst) 提供常见 Python 结构协议；[`contracts.py`](../api/algorithms/input_model/contracts.rst) 提供检查报告。应用仍可定义自己的协议。

[`operators.py`](../api/algorithms/input_model/operators.rst) 保存 BE 包装类与基本缩放操作，[`block_encoding.py`](../api/algorithms/input_model/block_encoding.rst) 提供 LCU、张量和小矩阵构造。它们是算法的组合工具。[`ode.py`](../api/algorithms/qode/ode.rst) 提供可替换线性求解接口；LCHS、Schrödingerization、CBMD 和 Carleman 各自维护实现。

## 领域应用

[QFVM](qfvm.md) 的物理参数、Roe 公式和经典数据更新放在 `applications`。[QHAM](qham.md) 的量子组装入口是 [`algorithms/qham.py`](../api/algorithms/input_model/qham.rst)；PDE 表达式、同伦推导和空间离散化支持放在 `applications/qham/`。这样可以阅读量子生成步骤，而不必同时展开全部数学推导代码。

早期演示实现集中在 `applications/legacy.py` 和 `algorithms/legacy.py`。它们用于兼容旧案例，当前文档不会将其作为新应用的默认入口。

## 兼容与迁移

旧导入路径由根目录的 `_compat.py` 集中转发。兼容层引用同一个函数或类，不保留另一份实现。仓内源码、测试与教程使用规范路径；外部代码可逐步迁移，见[导入路径迁移](compatibility.md)。
