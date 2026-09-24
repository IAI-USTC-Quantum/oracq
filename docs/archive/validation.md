# 首版验证记录

验证日期：2026-09-08。

本轮共有 52 个用例通过，未跳过任何用例：

| 范围 | 用例数 | 结果 |
|---|---:|---|
| 标准库语言核心和 BE 代数 | 37 | 全部通过。 |
| Draft 2020-12 JSON Schema 与序列化 | 4 | 全部通过。 |
| 真实 UnifiedQuantum 与 PySparQ 集成 | 11 | 全部通过。 |

核心与 Schema 测试在本仓库自己的 uv 环境中运行。原生集成测试使用本工作区已有的 QECC.Lang/.venv 解释器，只借用其安装的后端依赖；导入的 oracq 来自新仓库 src。两个后端源码提交见 backend-revisions.json。

真实集成覆盖寄存器 64 位及总位数超过 64、signed/rational 存储、嵌套 DEF、不同 QRAM 绑定、任意非零数据目标、切片查询、零控制、复系数相对相位、受控重复模块、伴随、原生注册表所有权和稀疏态预算。小规模实例对拍完整复幅度。

此外，Ruff 检查通过；安装后的 CLI 已完成 JSON 验证、QRAM 数据绑定执行和 OriginIR-ext 导出。源码包与 wheel 均可构建。文档中的相对链接与代码块闭合检查通过。

## 重现命令

```bash
uv sync --all-extras
uv run python -m unittest discover -s tests/core -v
uv run python -m unittest discover -s tests/schema -v
uv run ruff check src tests examples
PYTHONPATH=src ../QECC.Lang/.venv/bin/python -m unittest discover -s tests/integration -v
uv build
```

原生后端可以换成其他已安装对应依赖的解释器，源码包本身不依赖相邻仓库路径。

这些结果只针对首版语言核心和已列出的后端能力，不表示 QLSS、QODE 或 QHAM 算法已经实现或获得精度与复杂度保证。
