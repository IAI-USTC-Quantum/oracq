# pyqecclang

pyqecclang 是面向量子算法研究者的科学计算算法实现框架。用 Python 按访问模型编写可组合算法，生成保留模块结构的寄存器级中间表示（RIR）；保存尚未实现的 oracle，比较门网络、QRAM 和可逆算术等实现，完成绑定后进行数值验证和资源分析。

[算法研究教程](docs/tutorials/algorithm-research.md) 展示同一开放程序的三种实现及成本比较。

当前包版本为 **0.8.0**，RIR 格式为 **0.3**。算法通过普通 Python 协议定义输入和输出，不要求扩展语言类型系统。

## 文档

文档使用 Sphinx 生成，包含两条阅读路径和自动生成的 API 参考：

- [完整文档](docs/manual/index.md)：操作、oracle、算法约定、微分方程、应用与后端。
- [教程](docs/tutorials/index.md)：从第一个寄存器程序开始，逐步完成绑定、搜索、估计、Hamiltonian 与 QHAM 任务。
- [API 参考](docs/api/index.rst)：按照规范源码目录生成。
- [验证与适用范围](docs/manual/limits.md)：区分接口可用、线路见证和完整数值认证。

```bash
uv sync --locked --extra dev --extra docs
uv run sphinx-build -W --keep-going -b html docs out/docs/html
uv run sphinx-build -W --keep-going -b doctest docs out/docs/doctest
```

打开 `out/docs/html/index.html` 浏览生成站点。

## 一个最小程序

```python
from pyqecclang import Bits, Builder, export_originir, simulate

b = Builder("bell_pair", {"pair": Bits(2)})
b.h(b["pair"][0])
b.xor(b["pair"][0], b["pair"][1])
program = b.finish().program()

print(simulate(program).amplitudes)
print(export_originir(program).text)
```

结果在 `00` 与 `11` 上具有相等幅度。寄存器下标 0 是最低位；模块调用和 Repeat 在 RIR 与 JSON 中保留。

## 源码分类

```text
src/pyqecclang/
├── infrastructure/   RIR、构造器、验证、序列化、数学前端和后端
├── algorithms/       按类别组织的量子算法与组合工具
└── applications/     QFVM、Roe、QHAM 数学支持和案例目录
```

算法库包括查询、Fourier 算术、搜索与振幅放大、估计、变分电路、量子行走、求阶、简单纠错、Hamiltonian 演化、QLSS 和多种 QODE 方法。各类别的实现位于独立文件，详见[算法目录](docs/manual/algorithms/index.md)。

旧导入路径集中转发到同一份实现，新代码使用规范路径。迁移说明见[导入路径](docs/manual/compatibility.md)。

## 运行与检查

```bash
uv run python examples/algorithm_gallery.py
uv run python tools/check_project.py --docs
```

算法展示目录生成 22 个小实例的 RIR、OriginIR-ext 和读出说明。完整原生检查要求另一个已安装 `pysparq` 与 `uniqc` 的环境：

```bash
PYTHONPATH=src /path/to/backend/python examples/algorithm_gallery.py --native
```

语言核心没有第三方运行时依赖。可选后端在执行入口导入；生成产物、环境和构建文件均不提交。开发流程见 [CONTRIBUTING.md](CONTRIBUTING.md)。

QLSS/QODE/QHAM 等高级算法仍有数值精度、成功通道或收敛性待核验项。通用 QSP-HamSim 内核尚需提供；当前模乘采用有限规模置换合成，VQE/QAOA 的经典优化器由应用选择。
