# 目录整理与算法扩展面板

本轮先整理代码和文档，再扩展一批可实际生成线路的算法。所有工作在 pyqecclang 仓库内完成，保留已有未提交改动。

| 阶段 | 内容 | 状态 |
|---|---|---|
| R1 | infrastructure / algorithms / applications 分层；拆开 elementary、solvers、differential，保留旧导入兼容 | 已完成 |
| R2 | Sphinx + MyST + API 文档；完整文档与教程分开，历史记录归档 | 已完成 |
| R3 | oracle 查询、Fourier 算术、搜索与振幅放大、估计、变分算法、量子行走、数论、简单纠错 | 已完成 |
| R4 | 新算法数学见证、旧案例回归、真实后端验证、严格文档构建、打包与迁移验收 | 已完成 |

算法文件按用途命名：oracle_algorithms.py、fourier.py、search.py、estimation.py、variational.py、walks.py、number_theory.py、error_correction.py，以及现有的 hamiltonian.py、qlss.py、lchs.py、schrodingerization.py、cbmd.py、carleman.py、qham.py。

首批新增实现包含 Bernstein–Vazirani、Simon 采样与 GF(2) 后处理、QFT 加法、通用振幅放大、Hadamard/Swap test、标准振幅估计、MaxCut QAOA、VQE 测量电路、参数化 ansatz、周期 coined walk、模乘求阶及因子后处理、三位 bit/phase flip 编解码。门级输出和算法适用范围分别记录，不把小型模乘表的实现当作可扩展的 Shor 算术。

文档使用中文完整句子说明问题、输入、步骤和结果。以用户任务组织教程；规范给出明确契约；API 从源码生成。构建要求无 Sphinx warning，教程可运行。现有阶段报告移入 archive，保留引用但不混入主阅读路径。
