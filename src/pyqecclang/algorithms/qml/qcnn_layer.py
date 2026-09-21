"""QCNN 的量子构件（arXiv:1911.01117 §5.1）。

``qcnn_vector_prep`` 实现 Eq. (13)-(15) 的 QRAM 行/列装载：索引寄存器 p
选择一行，目标寄存器制备归一化行向量态 ``|A_p⟩``。旋转角与符号存为运行时 QRAM
bank——角度树（RY 半角约定，编址 (1<<depth)-1+prefix）加符号反冲，机制
与 QFVM 符号残差树一致；换数据只换内存表，线路不变。

``qcnn_inner_product`` 实现 Eq. (16)-(19)：p、q 均匀叠加，flag 上的
Hadamard 测试把内积编码进幅度。测量 (p, q, flag) 的概率满足
P0(p,q) = (1+<A_p|F_q>)/(2 H'W'D')（Eq. 20）。

两个构件的语义验证对照带角度量化的经典镜像，见 tests/core/test_qcnn.py。
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from pyqecclang.algorithms.input_model.operators import _name
from pyqecclang.algorithms.input_model.oracles import (
    XorDatabase,
    invoke,
    qram_database,
    resources_for,
)
from pyqecclang.algorithms.qml.qcnn import ConvSpec
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.ir import Bits, Ref, ValidationError, fuse


def _log2_exact(value: int) -> int:
    """校验 ``value`` 为 2 的幂，返回其以 2 为底的对数。"""
    if value < 1 or value & (value - 1):
        raise ValidationError("QCNN 量子构件要求行/列数为 2 的幂")
    return value.bit_length() - 1


def vector_angle_tables(
    rows: Sequence[Sequence[float]], angle_width: int
) -> tuple[dict[int, int], dict[int, int]]:
    """把一组等长行向量编译成角度树 bank 与符号 bank（运行时内存表）。

    角度按 RY 半角约定从平方部分和导出：RY 参数为
    2*atan2(sqrt(S_right), sqrt(S_left))；树节点编址 (1<<depth)-1+prefix
    （与 qram_state_prep 的编址一致）。返回 (angles, signs)：angles 键为
    fuse(index, node)，signs 键为 fuse(index, leaf)。

    Args:
        rows: 等长非零范数行向量序列，长度须为 2 的幂。
        angle_width: 角度字的量化位宽，弧角按 2π/2^angle_width 量化。

    Returns:
        tuple[dict[int, int], dict[int, int]]: 角度树 bank 与符号 bank 的
        运行时内存表。
    """
    length = len(rows[0])
    if length < 1 or length & (length - 1):
        raise ValidationError("向量长度必须是 2 的幂")
    width = length.bit_length() - 1
    angles: dict[int, int] = {}
    signs: dict[int, int] = {}
    scale = (1 << angle_width) / (2 * math.pi)
    for index, row in enumerate(rows):
        if any(len(r) != length for r in rows):
            raise ValidationError("所有行必须等长")
        squared = [value * value for value in row]
        if sum(squared) == 0:
            raise ValidationError("行向量的范数为零")
        for depth in range(width):
            bit = width - depth - 1
            for prefix in range(1 << depth):
                start = prefix << (bit + 1)
                left = sum(squared[start : start + (1 << bit)])
                right = sum(squared[start + (1 << bit) : start + (1 << (bit + 1))])
                node = (1 << depth) - 1 + prefix
                angle = 2 * math.atan2(math.sqrt(right), math.sqrt(left))
                angles[(index << width) | node] = round(angle * scale) % (1 << angle_width)
        for leaf, value in enumerate(row):
            signs[(index << width) | leaf] = int(value < 0)
    return angles, signs


def _constant_register(
    builder: Builder, name: str, width: int, value: int
) -> Ref:
    """X 制备的常量寄存器（从 |0> 出发，调用方负责对称恢复）。"""
    register = builder.local(name, Bits(width))
    for bit in range(width):
        if (value >> bit) & 1:
            builder.x(register[bit])
    return register


def _prep_node_rotation(
    builder: Builder,
    angles: XorDatabase,
    node_register: Ref,
    bit: int,
    angle_width: int,
) -> None:
    """一个树节点上的"查角-受控 RY 合成-反查"三步。

    bank 键为 (index<<width)|node：index 占高位，因此地址为
    fuse(node_register, index)。
    """
    invoke(
        builder,
        angles.operation,
        "angles",
        address=fuse(node_register, builder["index"]),
        data=builder["work"],
    )
    for k in range(angle_width):
        with builder.control(builder["work"][k], 1):
            builder.ry(
                builder["target"][bit], 2 * math.pi * (1 << k) / (1 << angle_width)
            )
    invoke(
        builder,
        angles.operation,
        "angles",
        address=fuse(node_register, builder["index"]),
        data=builder["work"],
    )


def qcnn_vector_prep(
    count: int, vector_width: int, angle_width: int, *, name: str | None = None
) -> Operation:
    """Eq. (15)：按索引把归一化向量态写入目标寄存器（QRAM 数据驱动）。

    寄存器：index（log2 count 位）、target（vector_width 位）、
    work（angle_width 位）。资源：angles（地址 fuse(index, node)，数据宽
    angle_width）与 signs（地址 fuse(index, leaf)，数据宽 1）。角度字按位
    加权合成为受控 RY：2*pi*2^k/2^angle_width。符号经 Z 反冲写入并在查询
    复净后只剩分支相位。index >= count 的行为未定义，调用方须保证叠加只
    覆盖有效行。

    Args:
        count: 向量个数，须为 2 的幂，决定 index 寄存器位宽。
        vector_width: 单个向量的量子位宽（向量长度为 2 的幂）。
        angle_width: 角度 bank 的字宽，旋转按其逐位加权合成。
        name: 生成的模块名；缺省按参数自动生成。

    Returns:
        Operation: 以 index 选取行、把归一化向量态写入 target 的制备操作。
    """
    index_width = _log2_exact(count)
    angles = qram_database(index_width + vector_width, angle_width)
    signs = qram_database(index_width + vector_width, 1)
    b = Builder(
        name or _name("qcnn_vector_prep", count, vector_width, angle_width),
        {
            "index": Bits(index_width),
            "target": Bits(vector_width),
            "work": Bits(angle_width),
        },
        resources_for(("angles", angles.operation), ("signs", signs.operation)),
        attributes={
            "algorithm": "qcnn_qram_row_preparation",
            "implementation": "angle_tree_with_sign_kickback",
            "qram_queries": 2 * vector_width,
            "correctness": "pending",
        },
    )
    for depth in range(vector_width):
        bit = vector_width - depth - 1
        for prefix in range(1 << depth):
            node = (1 << depth) - 1 + prefix
            node_register = _constant_register(b, f"node_{depth}_{prefix}", vector_width, node)
            if depth:
                with b.control(b["target"][bit + 1 :], prefix):
                    _prep_node_rotation(b, angles, node_register, bit, angle_width)
            else:
                _prep_node_rotation(b, angles, node_register, bit, angle_width)
            for bit_ in range(vector_width):
                if (node >> bit_) & 1:
                    b.x(node_register[bit_])
    invoke(
        b,
        signs.operation,
        "signs",
        address=fuse(b["target"], b["index"]),
        data=b["work"][:1],
    )
    b.z(b["work"][0])
    invoke(
        b,
        signs.operation,
        "signs",
        address=fuse(b["target"], b["index"]),
        data=b["work"][:1],
    )
    return b.finish()


def qcnn_inner_product(
    row_count: int,
    column_count: int,
    vector_width: int,
    angle_width: int,
    *,
    name: str | None = None,
) -> Operation:
    """Eq. (16)-(19)：Hadamard 型内积电路。

    寄存器：p、q（行/列索引）、flag（1 位）、vec（vector_width 位）、
    work（angle_width + 1 位）。p、q 均匀叠加后，flag=0 分支装载行向量、
    flag=1 分支装载列向量，最后的 Hadamard 把 (1+<A_p|F_q>)/2 编码进
    flag=0 的幅度。资源：行/列各一套角度与符号 bank（资源名前缀 row 与 col）。

    Args:
        row_count: 行向量个数，须为 2 的幂。
        column_count: 列向量个数，须为 2 的幂。
        vector_width: 行/列向量的量子位宽。
        angle_width: 角度 bank 的字宽，与制备操作一致。
        name: 生成的模块名；缺省按参数自动生成。

    Returns:
        Operation: 把 (1+<A_p|F_q>)/2 编码进 flag 幅度的内积估计操作。
    """
    row_prep = qcnn_vector_prep(row_count, vector_width, angle_width, name="qcnn_row_prep")
    column_prep = qcnn_vector_prep(
        column_count, vector_width, angle_width, name="qcnn_col_prep"
    )
    b = Builder(
        name or _name("qcnn_inner_product", row_count, column_count, vector_width, angle_width),
        {
            "p": Bits(_log2_exact(row_count)),
            "q": Bits(_log2_exact(column_count)),
            "flag": Bits(1),
            "vec": Bits(vector_width),
            "work": Bits(angle_width),
        },
        resources_for(
            ("row", row_prep),
            ("col", column_prep),
        ),
        attributes={
            "algorithm": "qcnn_hadamard_inner_product",
            "equation": "P(p,q,flag=0) = (1+<A_p|F_q>)/(2*rows*cols)",
            "correctness": "pending",
        },
    )
    for register in (b["p"], b["q"], b["flag"]):
        b.h(register)
    with b.control(b["flag"], 0):
        invoke(
            b,
            row_prep,
            "row",
            index=b["p"],
            target=b["vec"],
            work=b["work"],
        )
    with b.control(b["flag"], 1):
        invoke(
            b,
            column_prep,
            "col",
            index=b["q"],
            target=b["vec"],
            work=b["work"],
        )
    b.h(b["flag"])
    return b.finish()

def quantized_prepared_state(vector: Sequence[float], angle_width: int) -> list[float]:
    """角度量化镜像：行/列制备出的（量化）归一化向量。

    与电路语义一致：逐层旋转，每层角度取 bank 中的量化角字按位加权
    合成；符号按符号 bank 施加。返回实数列表。

    Args:
        vector: 待制备的实数向量，长度为 2 的幂且范数非零。
        angle_width: 角度字的量化位宽，与电路 bank 一致。

    Returns:
        list[float]: 与电路输出一致的量化归一化向量。
    """
    length = len(vector)
    if length & (length - 1):
        raise ValidationError("向量长度必须是 2 的幂")
    width = length.bit_length() - 1
    step = 2 * math.pi / (1 << angle_width)
    squared = [value * value for value in vector]
    if sum(squared) == 0:
        raise ValidationError("行向量的范数为零")
    states = {0: 1.0}
    for depth in range(width):
        bit = width - depth - 1
        nxt: dict[int, float] = {}
        for prefix, amplitude in states.items():
            start = prefix << (bit + 1)
            left = sum(squared[start : start + (1 << bit)])
            right = sum(squared[start + (1 << bit) : start + (1 << (bit + 1))])
            angle = 2 * math.atan2(math.sqrt(right), math.sqrt(left))
            word = round(angle / step) % (1 << angle_width)
            theta = sum(
                2 * math.pi * (1 << k) / (1 << angle_width)
                for k in range(angle_width)
                if (word >> k) & 1
            )
            low, high = prefix * 2, prefix * 2 + 1
            nxt[low] = nxt.get(low, 0.0) + amplitude * math.cos(theta / 2)
            nxt[high] = nxt.get(high, 0.0) + amplitude * math.sin(theta / 2)
        states = nxt
    magnitude = [states[index] for index in range(length)]
    norm = math.sqrt(sum(value * value for value in magnitude))
    return [m / norm * (1.0 if vector[index] >= 0 else -1.0) for index, m in enumerate(magnitude)]


def qcnn_sampled_layer(
    x: Sequence[float],
    kernel: Sequence[float],
    spec: ConvSpec,
    *,
    samples: int,
    eta: float,
    seed: int = 0,
) -> tuple[list[list[float]], dict[str, int]]:
    """论文 Algorithm 1 的采样驱动（小规模端到端）。

    量子语义分两层：行/列制备与 Hadamard 内积已由电路验证
    （qcnn_inner_product，Eq. 16-20）；条件旋转、幅度放大与 l_inf 层析
    （Eq. 28-38）在此按其精确分布建模——采样概率取 f(Y_pq)^2/sum f^2
    （Eq. 34），每次采样得到三元组 (p, q, f(Y_pq))，低于 eta 的像素视为
    未采样置零（Eq. 36）。采样值按 Eq. 39 映射进池化区域，经 QRAM 在线
    覆写规则（max 保高 / average 均摊）得到输出张量。返回
    (pooled_rows, stats)。

    Args:
        x: 按行主序扁平存储的输入张量数据。
        kernel: 卷积核系数的扁平序列。
        spec: 卷积层规格，决定展开、激活上限与池化方式。
        samples: 采样次数，取正整数。
        eta: 像素保留阈值，低于该值的采样视为未采样置零。
        seed: 随机数生成器种子。

    Returns:
        tuple[list[list[float]], dict[str, int]]: 池化输出矩阵与采样统计
        （sampled/kept 计数）。
    """
    import random

    from pyqecclang.algorithms.qml.qcnn import (
        QCNNQRAM,
        cap_relu,
        im2col,
        kernel_columns,
    )

    rows = im2col(x, spec)
    columns = kernel_columns(kernel, spec)
    if len(rows) & (len(rows) - 1) or len(columns) & (len(columns) - 1):
        raise ValidationError("采样驱动要求行数/列数为 2 的幂")
    qram = QCNNQRAM(rows)
    prepared_rows = [quantized_prepared_state(row, 10) for row in rows]
    prepared_columns = [quantized_prepared_state(column, 10) for column in columns]
    values: dict[tuple[int, int], float] = {}
    for p, a in enumerate(prepared_rows):
        for q, f in enumerate(prepared_columns):
            overlap = sum(u * v for u, v in zip(a, f, strict=True))
            inner = overlap * qram.norm(p) * _column_norm(columns[q])
            values[(p, q)] = cap_relu(inner, spec.cap)
    weight = {key: value * value for key, value in values.items()}
    total = sum(weight.values())
    if total == 0:
        return _zero_pooled(spec), {"sampled": 0, "kept": 0}
    rng = random.Random(seed)
    oh, ow, _ = spec.output_shape
    ph, pw, d = spec.pooled_shape
    pooled = [[0.0] * d for _ in range(ph * pw)]
    counts = [[0] * d for _ in range(ph * pw)]
    kept = 0
    for _ in range(samples):
        pick, acc = rng.random() * total, 0.0
        for key, w in weight.items():
            acc += w
            if acc >= pick:
                p, q = key
                break
        value = values[(p, q)]
        if value < eta:
            continue
        kept += 1
        i, j = p // ow, p % ow
        index = (i // spec.pool) * pw + (j // spec.pool)
        if spec.pool_kind == "max":
            pooled[index][q] = max(pooled[index][q], value)
        else:
            pooled[index][q] = (pooled[index][q] * counts[index][q] + value) / (counts[index][q] + 1)
        counts[index][q] += 1
    stats = {"sampled": samples, "kept": kept}
    return pooled, stats


def _column_norm(column: Sequence[float]) -> float:
    """返回列向量的欧几里得范数。"""
    return math.sqrt(sum(value * value for value in column))


def _zero_pooled(spec: ConvSpec) -> list[list[float]]:
    """与 ``spec.pooled_shape`` 一致的全零池化输出张量。"""
    ph, pw, d = spec.pooled_shape
    return [[0.0] * d for _ in range(ph * pw)]
