# 算法目录

算法库按用途组织。下表给出入口文件、已实现的内容以及使用时需要留意的边界。API 参考列出了完整签名。

| 类别 | 文件 | 实现与范围 |
|---|---|---|
| Oracle 查询 | `oracle_algorithms.py` | D-J、Bernstein–Vazirani、Simon 采样；Simon 的 GF(2) 消元在经典侧进行 |
| Fourier 变换与算术 | `fourier.py` | 正/逆 QFT、零宽 work 适配、模 2^n 的 Fourier 加法 |
| 搜索与放大 | `search.py` | Grover、Grover iterate、成功子空间振幅放大 |
| 估计与重叠 | `estimation.py` | QPE、标准振幅估计、Hadamard test、Swap test |
| 变分算法 | `variational.py` | 参数化 ansatz、MaxCut QAOA、Pauli 测量和 VQE 测量电路集合 |
| 量子行走 | `walks.py` | 周期格点上的 Hadamard coined walk |
| 数论 | `number_theory.py` | 有限规模模乘置换、QPE 求阶、连分数因子候选后处理 |
| 简单纠错 | `error_correction.py` | 三位 bit/phase flip 重复码的编码与相干恢复 |
| Hamiltonian 演化 | `hamiltonian.py` | Pauli 项演化、Trotter 组合、Taylor BE、可注入的 QSP 接口 |
| 矩阵变换 | `transforms.py` | qubitization、显式相位序列 QSVT、oblivious amplification 的组装 |
| 线性系统 | `qlss.py` | 问题契约、Costa walk/filter、CKS 基础 Chebyshev/LCU 路线 |
| 线性演化 | `ode.py` | QODE 协议和 Euler history 组装；具体方法在独立文件中 |
| 微分方程方法 | `lchs.py`、`schrodingerization.py`、`cbmd.py`、`carleman.py` | 各自维护输入模型、配置和生成步骤 |
| QHAM | `qham.py` | 从有限 HAM 闭包构造 QODE 输入及物理输出通道 |

## 如何选择起点

如果输入是布尔函数，先看查询与搜索类别。如果需要估计概率或期望值，先看 `estimation`。已有 Hamiltonian 项分解时可以选择 Trotter；只有算子 BE 时，需要选择能够消费该访问模型的实现。

PDE 算法先构造空间离散化的输入，再选择 QODE 或线性化路线。数据在 QRAM 中并不意味着矩阵已经被 block encoded，输入适配步骤仍需明确。

## 运行展示目录

```bash
uv run python examples/algorithm_gallery.py
PYTHONPATH=src /path/to/backend/python examples/algorithm_gallery.py --native
```

展示目录包含 22 个小实例，也包含同一算法的不同应用，例如振幅估计与量子计数。它们用于说明接口和读出方法，不能视为 22 种互不相关的算法。

## 实现依据

Fourier 加法采用 [Draper 的 QFT 加法构造](https://arxiv.org/abs/quant-ph/0008033)。振幅放大和估计采用 [Brassard 等人的框架](https://arxiv.org/abs/quant-ph/0005055)。QAOA 的 cost/mixer 分层依据 [Farhi 等人的原始算法](https://arxiv.org/abs/1411.4028)。求阶和经典因子后处理依据 [Shor 的构造](https://arxiv.org/abs/quant-ph/9508027)。

当前模乘使用有上限的置换合成，默认最多 8 位；它用于检查求阶接口与线路，不代表已实现可扩展的 Shor 模算术。VQE 和 QAOA 提供量子电路，经典优化器由应用选择。QSVT 接收调用方提供的相位序列，不包含通用相位求解器。

## 算法页面

每个算法一页，说明接口、输入模型、实现要点与验证证据的位置。页面按文件名排序；每页首行下方注明类别（C1–C6 或应用层）与所属模块，类别定义见[验证计划](../../development/validation-plan.md)。

```{toctree}
:maxdepth: 1
:caption: 算法页面
:glob:

*
```
