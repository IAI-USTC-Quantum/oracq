"模块化 BE、投影、LCU 与小矩阵表示组合。"

from __future__ import annotations

import cmath
import itertools
import math

from pyqecclang.algorithms.input_model.operators import BlockEncoding, _name, identity, scale
from pyqecclang.algorithms.input_model.oracles import (
    annotate,
    gate_state_prep,
    invoke,
    resources_for,
)
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError


def reflect_zero(builder, register, *, positive=False):
    """在给定寄存器的零态上追加符号翻转反射 ``I-2|0><0|``。

    Args:
        builder: 追加门的目标 Builder。
        register: 参与反射的寄存器，可传入多个寄存器拼接后的视图。
        positive: 为 True 时附加全局相位 pi，得到 ``2|0><0|-I``，零态分量取正号。

    行走与振幅放大类算法共享该构件。"""
    if positive:
        builder.global_phase(math.pi)
    if register.width:
        with builder.control(register, 0):
            builder.global_phase(math.pi)
    else:
        builder.global_phase(math.pi)


def pad_signal(a: BlockEncoding, width):
    """把 BE 的信号寄存器扩张到指定位宽。

    Args:
        a: 输入 block encoding。
        width: 新的信号位宽，不得小于原信号位宽。

    Returns:
        BlockEncoding: 角块语义与 alpha 不变，新增的高位信号恒为零，供同签名晚绑定。

    Raises:
        ValidationError: 试图缩小信号空间。
    """
    if width < a.signal_qubits:
        raise ValidationError("不能缩小 BE 信号空间")
    b = Builder(
        _name("pad_be", a.operation, width),
        {"target": Bits(a.width), "signal": Bits(width)},
        resources_for(("a", a.operation)),
    )
    invoke(b, a.operation, "a", target=b["target"], signal=b["signal"][: a.signal_qubits])
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=a.alpha))


def tensor(a: BlockEncoding, b: BlockEncoding):
    """构造两个 BE 的张量积。

    Args:
        a: 作用于目标高位的 BE。
        b: 作用于目标低位的 BE。

    Returns:
        BlockEncoding: 编码 A⊗B，alpha 为两者之积；signal 同样按 a 在低、b 在高拼接。
    """
    out = Builder(
        _name("tensor", a.operation, b.operation),
        {"target": Bits(a.width + b.width), "signal": Bits(a.signal_qubits + b.signal_qubits)},
        resources_for(("a", a.operation), ("b", b.operation)),
    )
    invoke(
        out,
        b.operation,
        "b",
        target=out["target"][: b.width],
        signal=out["signal"][a.signal_qubits :],
    )
    invoke(
        out,
        a.operation,
        "a",
        target=out["target"][b.width :],
        signal=out["signal"][: a.signal_qubits],
    )
    return BlockEncoding(annotate(out.finish(), "block_encoding", be_alpha=a.alpha * b.alpha))


def adjoint_be(a):
    """返回编码伴随矩阵 A† 的 BE。

    Args:
        a: 输入 block encoding。

    Returns:
        BlockEncoding: 在伴随上下文中调用原操作，alpha 保持不变。
    """
    out = Builder(
        _name("adjoint_be", a.operation),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits)},
        resources_for(("a", a.operation)),
    )
    with out.adjoint():
        invoke(out, a.operation, "a", target=out["target"], signal=out["signal"])
    return BlockEncoding(annotate(out.finish(), "block_encoding", be_alpha=a.alpha))


def lcu(terms):
    """以 PREPARE/SELECT 结构组装若干 BE 的线性组合。

    Args:
        terms: (系数, BE) 二元组序列；零系数项被剔除。

    Returns:
        BlockEncoding: 编码 Σ c_j A_j，alpha 为 ``Σ |c_j|*alpha_j``；仅一项时退化为 ``scale``。

    Raises:
        ValidationError: 没有非零项或各项目标宽度不一致。

    复系数的相位经选择器控制下的全局相位实现；signal 为选择位与各分支信号位的拼接，末尾逆制备恢复选择器。"""
    from pyqecclang.algorithms.input_model.interfaces import as_block_encoding

    terms = tuple((complex(c), as_block_encoding(a)) for c, a in terms if c != 0)
    if not terms:
        raise ValidationError("LCU 至少需要一个非零项")
    if len(terms) == 1:
        return scale(terms[0][0], terms[0][1])
    width = terms[0][1].width
    if any(a.width != width for _, a in terms):
        raise ValidationError("LCU 目标宽度不匹配")
    alpha = sum(abs(c) * a.alpha for c, a in terms)
    selector_width = (len(terms) - 1).bit_length()
    work_width = max(a.signal_qubits for _, a in terms)
    weights = [math.sqrt(abs(c) * a.alpha / alpha) for c, a in terms]
    weights += [0] * ((1 << selector_width) - len(weights))
    prep = gate_state_prep(weights)
    resources = resources_for(*[(f"term{i}", a.operation) for i, (_, a) in enumerate(terms)])
    b = Builder(
        _name("lcu_many", *(a.operation for _, a in terms), tuple(c for c, _ in terms)),
        {"target": Bits(width), "signal": Bits(selector_width + work_width)},
        resources,
    )
    selector, signal = b["signal"][:selector_width], b["signal"][selector_width:]
    invoke(b, prep.operation, target=selector, work=selector[:0])
    for i, (coefficient, a) in enumerate(terms):
        with b.control(selector, i):
            b.global_phase(cmath.phase(coefficient))
            invoke(b, a.operation, f"term{i}", target=b["target"], signal=signal[: a.signal_qubits])
    with b.adjoint():
        invoke(b, prep.operation, target=selector, work=selector[:0])
    return BlockEncoding(
        annotate(b.finish(), "block_encoding", be_alpha=alpha, lcu_terms=len(terms))
    )


def kronecker_sum(a, b=None):
    """构造两个 BE 的 Kronecker 和 A⊗I+I⊗B。

    Args:
        a: 第一个 BE。
        b: 第二个 BE；省略时取 a 自身。

    Returns:
        BlockEncoding: 两个张量项的 LCU，alpha 为两者 alpha 之和。
    """
    b = a if b is None else b
    return lcu([(1, tensor(a, identity(b.width))), (1, tensor(identity(a.width), b))])


def projector(width, accepted):
    """构造到指定基态子空间的投影 BE。

    Args:
        width: target 位宽。
        accepted: 被接受的基态整数值集合；重复值会合并并排序。

    Returns:
        BlockEncoding: 零信号角块为对角投影，集合外的基态被打入 signal 分支，alpha 为 1。
    """
    accepted = tuple(sorted(set(accepted)))
    out = Builder(_name("projector", width, accepted), {"target": Bits(width), "signal": Bits(1)})
    for value in range(1 << width):
        if value not in accepted:
            with out.control(out["target"], value):
                out.x(out["signal"])
    return BlockEncoding(annotate(out.finish(), "block_encoding", be_alpha=1.0))


def direct_sum(a, b):
    """构造同宽矩阵的直和。

    Args:
        a: 选择位为 0 时生效的 BE。
        b: 选择位为 1 时生效的 BE，宽度必须与 a 相同。

    Returns:
        BlockEncoding: 目标最高位作为选择位的块对角组合，alpha 为两者之和。

    Raises:
        ValidationError: 两个 BE 宽度不同。
    """
    if a.width != b.width:
        raise ValidationError("当前 direct_sum 需要同宽矩阵")
    return lcu([(1, tensor(projector(1, [0]), a)), (1, tensor(projector(1, [1]), b))])


def truncated_shift(width, last):
    """构造截断上移位 BE。

    Args:
        width: target 位宽。
        last: 截断阈值，范围为 1..2**width-1。

    Returns:
        BlockEncoding: 零信号角块把基态 v 映射到 v+1（v 小于 last），其余基态被打入 signal 分支，alpha 为 1。

    Raises:
        ValidationError: last 不在允许范围内。
    """
    if not 1 <= last < 1 << width:
        raise ValidationError("截断移位范围无效")
    out = Builder(_name("shift", width, last), {"target": Bits(width), "signal": Bits(1)})
    for value in range(last, 1 << width):
        with out.control(out["target"], value):
            out.x(out["signal"])
    out.add_const(out["target"].reinterpret("uint"), 1)
    return BlockEncoding(annotate(out.finish(), "block_encoding", be_alpha=1.0))


def pauli_word(word):
    """把 Pauli 字符串编码为无信号位的 BE。

    Args:
        word: I/X/Y/Z 字符串，第一个字符作用在最低位。

    Returns:
        BlockEncoding: 逐位显式单量子位门序列，alpha 为 1。

    Raises:
        ValidationError: 出现 I/X/Y/Z 之外的字符。
    """
    b = Builder("pauli_" + word, {"target": Bits(len(word)), "signal": Bits(0)})
    for bit, letter in enumerate(word):
        if letter not in "IXYZ":
            raise ValidationError("Pauli 字只允许 I/X/Y/Z")
        if letter != "I":
            b.gate(letter.lower(), b["target"][bit])
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))


def matrix_pauli_encoding(matrix, *, drop_tolerance=1e-12):
    """小型应用的显式门实现；不宣称矩阵输入或经典展开具有量子加速。"""
    matrix = tuple(tuple(complex(v) for v in row) for row in matrix)
    d = len(matrix)
    if d < 2 or d & (d - 1) or any(len(row) != d for row in matrix):
        raise ValidationError("矩阵必须是二的幂维方阵")
    n = (d - 1).bit_length()
    if n > 5:
        raise ValidationError("显式 Pauli 展开仅用于最多 5 位的小实例；大实例使用访问 oracle")
    terms = []
    for letters in itertools.product("IXYZ", repeat=n):
        coefficient = 0j
        for column in range(d):
            row, phase = column, 1 + 0j
            for bit, letter in enumerate(letters):
                value = (column >> bit) & 1
                if letter in "XY":
                    row ^= 1 << bit
                if letter == "Y":
                    phase *= -1j if value else 1j
                elif letter == "Z":
                    phase *= -1 if value else 1
            coefficient += phase.conjugate() * matrix[row][column]
        coefficient /= d
        if abs(coefficient) > drop_tolerance:
            terms.append((coefficient, pauli_word("".join(letters))))
    if not terms:
        from pyqecclang.algorithms.input_model.operators import zero

        return zero(n)
    return lcu(terms)
