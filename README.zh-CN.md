# oracq

[English](README.md) · **简体中文**

[![CI](https://github.com/IAI-USTC-Quantum/oracq/actions/workflows/python-ci.yml/badge.svg)](https://github.com/IAI-USTC-Quantum/oracq/actions/workflows/python-ci.yml)
[![Docs](https://github.com/IAI-USTC-Quantum/oracq/actions/workflows/docs.yml/badge.svg)](https://github.com/IAI-USTC-Quantum/oracq/actions/workflows/docs.yml)
[![PyPI](https://img.shields.io/pypi/v/oracq.svg)](https://pypi.org/project/oracq/)
[![Python](https://img.shields.io/pypi/pyversions/oracq.svg)](https://pypi.org/project/oracq/)
[![License: MIT](https://img.shields.io/badge/License-MIT-informational.svg)](LICENSE)

oracq 是面向量子算法研究者的科学计算算法实现框架。用 Python 按访问模型编写可组合算法，生成保留模块结构的寄存器级中间表示（RIR）；保存尚未实现的 oracle，比较门网络、QRAM 和可逆算术等实现，完成绑定后进行数值验证和资源分析。

当前包版本为 **0.1.0**，RIR 格式为 **0.1**（[规范全文](docs/zh/reference/rir.md)）。算法通过普通 Python 协议定义输入和输出，不要求扩展语言类型系统。

## 核心概念速览

| 概念 | 一句话 | 手册 / 规范 | API 参考 |
|---|---|---|---|
| RIR | 寄存器级中间表示；模块调用与 Repeat 结构保留到文本形式 | [RIR 0.1 规范](docs/zh/reference/rir.md) | [ir](docs/api/infrastructure/ir.rst) |
| 构造器 | `Builder`/`Operation` 三段式生成模块与程序 | [核心概念](docs/zh/manual/concepts.md#模块与未完成的实现) | [builder](docs/api/infrastructure/builder.rst) |
| 寄存器与视图 | bits/uint/qubit 解释、切片与重解释，下标 0 是最低位 | [核心概念](docs/zh/manual/concepts.md#寄存器与视图) | [ir](docs/api/infrastructure/ir.rst) |
| 开放 oracle 与绑定 | 先声明后绑定：能力合取、候选实现比较、QRAM 捕获 | [绑定教程](docs/zh/tutorials/oracle-binding.md) | [linking](docs/api/infrastructure/linking.rst) |
| 算法契约 | 可检查的输入协议、能力规格与验收报告 | [契约](docs/zh/manual/contracts.md) | [contracts](docs/api/algorithms/input_model/contracts.rst) |
| 结构验证 | 生成期校验位宽、重叠与结构约束 | [核心概念](docs/zh/manual/concepts.md) | [validation](docs/api/infrastructure/validation.rst) |
| 序列化 | RIR 默认 YAML，JSON 可选，两种文本逐字段一致 | [RIR 规范](docs/zh/reference/rir.md) | [serialization](docs/api/infrastructure/serialization.rst) |
| 执行与读出 | 寄存器参考执行器与宿主读出 | [后端](docs/zh/manual/backends.md) | [execution](docs/api/infrastructure/execution.rst)、[readout](docs/api/infrastructure/readout.rst) |
| QMem | QRAM 资源的指针式读写与多维视图 | [QMem](docs/zh/manual/qmem.md) | [qmem](docs/api/infrastructure/qmem.rst) |
| QRAM 数据文件 | `*.qram.yaml` 内存定义的加载与写出 | [QRAM 内存](docs/zh/reference/qram-memory.md) | [qram_schema](docs/api/infrastructure/qram_schema.rst) |
| 数学函数前端 | Python 可调用编译为可逆线路 | [数学函数](docs/zh/manual/math-functions.md) | [mathfunc](docs/api/infrastructure/mathfunc.rst) |
| 资源估计 | T / 旋转 / QRAM 访问的计数模型与开放分析 | [资源估计](docs/zh/manual/resource-estimation.md) | [estimate](docs/api/infrastructure/estimate.rst) |
| 后端导出 | OriginIR-ext、严格网表、PySparQ、quantikz | [后端](docs/zh/manual/backends.md) | [backends](docs/api/infrastructure/backends/originir.rst) |

## 一个最小程序

```python
from oracq import Bits, Builder, export_originir, simulate

b = Builder("bell_pair", {"pair": Bits(2)})
b.h(b["pair"][0])
b.xor(b["pair"][0], b["pair"][1])
program = b.finish().program()

print(simulate(program).amplitudes)
print(export_originir(program).text)
```

结果在 `00` 与 `11` 上具有相等幅度。模块调用和 Repeat 在 RIR 与 YAML/JSON 文本中保留，不在生成阶段展开。

本例用到的 API：[`Builder`](docs/api/infrastructure/builder.rst) · [`Bits`](docs/api/infrastructure/ir.rst) · [`simulate`](docs/api/infrastructure/execution.rst) · [`export_originir`](docs/api/infrastructure/backends/originir.rst)。逐步讲解见教程[第一个寄存器程序](docs/zh/tutorials/first-program.md)。

## 典型工作流

| 步骤 | 入口 | 深入阅读 |
|---|---|---|
| 1. 声明开放 oracle | [`declare`](docs/api/algorithms/input_model/oracles.rst)、能力与规格 | [算法契约](docs/zh/manual/contracts.md) |
| 2. 构造程序 | [`Builder`](docs/api/infrastructure/builder.rst) 组合模块调用 | [教程：第一个程序](docs/zh/tutorials/first-program.md) |
| 3. 结构验证 | [`validate`](docs/api/infrastructure/validation.rst) | [核心概念](docs/zh/manual/concepts.md) |
| 4. 绑定实现 | [`bind`](docs/api/infrastructure/linking.rst)、`bind_with_report` 比较候选 | [教程：替换 oracle](docs/zh/tutorials/oracle-binding.md)、[算法研究教程](docs/zh/tutorials/algorithm-research.md) |
| 5. 导出 / 执行 / 估计 | [`export_originir`](docs/api/infrastructure/backends/originir.rst)、[`simulate`](docs/api/infrastructure/execution.rst)、[`estimate_resources`](docs/api/infrastructure/estimate.rst) | [后端](docs/zh/manual/backends.md)、[资源估计](docs/zh/manual/resource-estimation.md) |

## 算法库

算法库按用途组织为十个子包，完整清单与选型指引见[算法目录](docs/zh/manual/algorithms/index.md)，每个算法一页手册（接口、实现要点、验证方案与已知缺口）。各类代表：

| 类别 | 代表算法页 |
|---|---|
| 查询与搜索 | [Grover 搜索](docs/zh/manual/algorithms/grover.md)、[振幅估计](docs/zh/manual/algorithms/qae.md)、[量子计数](docs/zh/manual/algorithms/quantum-counting.md) |
| 基础查询算法 | [Deutsch–Jozsa](docs/zh/manual/algorithms/deutsch-jozsa.md)、[Simon](docs/zh/manual/algorithms/simon.md)、[Bernstein–Vazirani](docs/zh/manual/algorithms/bernstein-vazirani.md) |
| Fourier 与算术 | [QFT](docs/zh/manual/algorithms/qft.md)、[Fourier 加法](docs/zh/manual/algorithms/fourier-addition.md)、[求阶](docs/zh/manual/algorithms/order-finding.md) |
| Hamiltonian 演化 | [Trotter](docs/zh/manual/algorithms/trotter.md)、[QSP 相位合成](docs/zh/manual/algorithms/qsp-phase-synthesis.md)、[QSVT HamSim](docs/zh/manual/algorithms/qsvt-hamiltonian-simulation.md) |
| 估计与测试 | [QPE](docs/zh/manual/algorithms/qpe.md)、[Hadamard test](docs/zh/manual/algorithms/hadamard-test.md)、[Swap test](docs/zh/manual/algorithms/swap-test.md) |
| 量子线性系统 | [Costa 行走](docs/zh/manual/algorithms/costa-walk.md)、[CKS](docs/zh/manual/algorithms/cks.md)、[VTAA-CKS](docs/zh/manual/algorithms/vtaa-cks.md) |
| QODE / QPDE | [Schrödingerization](docs/zh/manual/algorithms/schrodingerization.md)、[LCHS](docs/zh/manual/algorithms/lchs.md)、[Carleman](docs/zh/manual/algorithms/carleman.md) |
| 变分与优化 | [VQE](docs/zh/manual/algorithms/vqe.md)、[QAOA MaxCut](docs/zh/manual/algorithms/qaoa-maxcut.md)、[DQI](docs/zh/manual/algorithms/dqi.md) |
| 量子行走 | [coined walk](docs/zh/manual/algorithms/coined-cycle-walk.md)、[Szegedy](docs/zh/manual/algorithms/szegedy-walk.md)、[MNRS](docs/zh/manual/algorithms/mnrs-search.md) |
| 量子机器学习 | [QPCA](docs/zh/manual/algorithms/qpca.md)、[QCNN](docs/zh/manual/algorithms/qcnn.md)、KP 推荐系统（[API](docs/api/algorithms/qml/recommendation.rst)，手册页待补） |
| 数据加载与输入模型 | [态制备](docs/zh/manual/algorithms/state-preparation.md)、[XOR 数据库](docs/zh/manual/algorithms/xor-database.md)、[Select-Swap QROM](docs/zh/manual/algorithms/select-swap.md) |
| 量子纠错 | [重复码](docs/zh/manual/algorithms/repetition-codes.md) |

## 输入模型与算子

算法通过五类访问模型消费输入，彼此可组合：

- 算子包装与基本组合（`identity`/`product`/`scale`/LCU）：[手册](docs/zh/manual/operators.md)、[API](docs/api/algorithms/input_model/operators.rst)；块编码代数见 [block_encoding](docs/api/algorithms/input_model/block_encoding.rst)。
- Oracle 范式（XorDatabase、StatePreparation、StateOracle、SparseAccess）：[API](docs/api/algorithms/input_model/oracles.rst)。
- 量子数据结构 QVector/QMatrix：[手册](docs/zh/manual/qdata.md)、[API](docs/api/algorithms/input_model/qdata.rst)。
- 密度矩阵与 Gibbs 态、谱与低秩分解：[density](docs/api/algorithms/input_model/density.rst)、[spectral](docs/api/algorithms/input_model/spectral.rst)、[lowrank](docs/api/algorithms/input_model/lowrank.rst)。

## 应用层

- **QFVM**（量子流体求解）：[手册](docs/zh/manual/qfvm.md)、[API](docs/api/applications/qfvm.rst)、输入模型审阅[规范](docs/zh/reference/qfvm-input-models.md)。
- **QHAM**（PDE → HAM → QHAM 流程）：[手册](docs/zh/manual/qham.md)、[推导规范](docs/zh/reference/qham-derivation.md)、[API](docs/api/applications/qham/linearization.rst)、[教程](docs/zh/tutorials/qham.md)。
- **Roe 矩阵元**：[roe](docs/api/applications/roe.rst)、[roe_formulas](docs/api/applications/roe_formulas.rst)。
- **案例目录**：22 个参考工作负载（[catalog](docs/api/applications/catalog.rst)）与展示生成器（[gallery](docs/api/applications/gallery.rst)、[教程](docs/zh/tutorials/gallery.md)）。

## 文档地图

按角色选择阅读路径：

- **上手**：[第一个程序](docs/zh/tutorials/first-program.md) → [核心概念](docs/zh/manual/concepts.md) → [替换 oracle](docs/zh/tutorials/oracle-binding.md)。
- **算法研究**：[算法研究教程](docs/zh/tutorials/algorithm-research.md) → [算法目录](docs/zh/manual/algorithms/index.md) → [验证覆盖矩阵](docs/zh/development/validation-coverage.md)，适用边界见[验证与适用范围](docs/zh/manual/limits.md)。
- **微分方程应用**：[教程](docs/zh/tutorials/differential-equations.md) → [手册](docs/zh/manual/differential-equations.md) → [QHAM 手册](docs/zh/manual/qham.md)。
- **后端工程**：[后端手册](docs/zh/manual/backends.md) → [兼容性审阅](docs/zh/reference/backend-compatibility.md) → [基础设施 API](docs/api/infrastructure/index.rst)。
- **语言与规范**：[RIR 规范](docs/zh/reference/rir.md)、[生成层边界](docs/zh/reference/language.md)、[开放 IR](docs/zh/reference/open-ir.md)、[QRAM 内存格式](docs/zh/reference/qram-memory.md)、[数学 IR](docs/zh/reference/math-ir.md)。

## 安装 · 构建 · 检查

```bash
uv sync --locked --extra dev --extra docs
uv run python tools/build_docs.py --lang all
```

构建双语两个语言树（warning 视为错误）：英文输出到 `out/docs/en/html`，中文输出到 `out/docs/zh/html`，并执行 doctest。打开 `out/docs/zh/html/index.html` 浏览中文站点。算法展示与工程检查：

```bash
uv run python examples/algorithm_gallery.py
uv run python tools/check_project.py --docs
```

算法展示目录生成 22 个小实例的 RIR、OriginIR-ext 和读出说明。完整原生检查要求另一个已安装 `pysparq` 与 `uniqc` 的环境：

```bash
PYTHONPATH=src /path/to/backend/python examples/algorithm_gallery.py --native
```

语言核心只依赖 PyYAML 做文本序列化。可选后端在执行入口导入；生成产物、环境和构建文件均不提交。开发流程见 [CONTRIBUTING.zh-CN.md](CONTRIBUTING.zh-CN.md) 与[文档编写规范](docs/zh/development/writing-docs.md)；源码分类与迁移说明见[架构](docs/zh/manual/architecture.md)与[导入路径](docs/zh/manual/compatibility.md)。

QLSS/QODE/QHAM 等高级算法仍有数值精度、成功通道或收敛性待核验项（逐项状态见[验证覆盖矩阵](docs/zh/development/validation-coverage.md)）。通用 QSP-HamSim 内核尚需提供；当前模乘采用有限规模置换合成，VQE/QAOA 的经典优化器由应用选择。
