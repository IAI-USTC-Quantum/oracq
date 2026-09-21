# pyqecclang 开发约定

本仓库是独立的 Python 包。所有工作只在本仓库进行；不要为兼容性测试修改相邻的 UnifiedQuantum、QRAM-Simulator 或 QECC.Lang。

RIR 是架构中心。修改 IR 时，必须同步 docs/reference/rir.md、JSON Schema、序列化器和语义测试。寄存器名、位宽和视图必须保留到后端降低阶段。模块调用和 Repeat 不能在生成或 JSON 序列化阶段无条件展开。

核心不依赖量子后端。后端导出与执行分离，原生依赖只在执行入口导入。代码标识符使用英文，文档和注释使用中文。

核心验证命令是 python -m unittest discover -s tests/core -v。真实后端验证使用具有 uniqc 和 pysparq 的解释器运行 tests/integration，不使用模拟替身或 skip 代替真实对拍。静态检查使用 ruff check src tests examples，类型检查使用 mypy src examples。

不要提交虚拟环境、缓存、构建产物和模拟输出。没有明确要求时不执行 git commit 或 push。

源码按 infrastructure、algorithms、applications 分层。新增算法放在相应类别文件，不写入旧导入兼容层。文档使用 Sphinx，分为 manual、tutorials、reference 和 api；历史资料放 archive。文档改动运行 HTML 与 doctest builder，并将 warning 视为错误。
