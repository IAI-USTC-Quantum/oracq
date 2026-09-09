"""正逆 Fourier 变换、接口适配与不使用进位寄存器的模加法。"""

from __future__ import annotations

import math

from pyqecclang.algorithms.contracts import positive_integer
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits


def qft(width):
    """生成正号离散 Fourier 变换。

    Args:
        width: target 位宽，范围为 1..64。

    Returns:
        Operation: 只有 target 寄存器，包含末尾交换，位序为 little endian。

    矩阵元素为 exp(2πi*x*y/2**width)/sqrt(2**width)。"""
    positive_integer(width, "qft.width", maximum=64)
    b = Builder(f"qft_{width}", {"target": Bits(width)})
    for high in reversed(range(width)):
        b.h(b["target"][high])
        for low in reversed(range(high)):
            with b.control(b["target"][low]):
                b.gate("phase", b["target"][high], math.pi / (1 << (high - low)))
    for bit in range(width // 2):
        b.swap(b["target"][bit], b["target"][width - bit - 1])
    return b.finish()


def qft_with_work(width):
    """保留早期零宽 work 接口的 QFT 适配。"""
    b = Builder("qft_with_work_" + str(width), {"target": Bits(width), "work": Bits(0)})
    b.call(qft(width), target=b["target"])
    return b.finish()


def inverse_qft(width):
    """生成 QFT 的伴随操作，保持模块调用。"""
    b = Builder("inverse_qft_" + str(width), {"target": Bits(width)})
    with b.adjoint():
        b.call(qft(width), target=b["target"])
    return b.finish()


def fourier_add(width):
    """用 QFT 实现无进位寄存器的模加法。

    Args:
        width: a 和 b 的共同位宽。

    Returns:
        Operation: ``|a,b> -> |a,(a+b) mod 2**width>``。a 保留，无额外公开工作区。

    Fourier 变换和逆变换保留为模块调用。"""
    positive_integer(width, "fourier_add.width", maximum=64)
    b = Builder("fourier_add_" + str(width), {"a": Bits(width), "b": Bits(width)})
    b.call(qft(width), target=b["b"])
    for i in range(width):
        for j in range(width - i):
            with b.control(b["a"][i]):
                b.gate("phase", b["b"][j], 2 * math.pi / (1 << (width - i - j)))
    b.call(inverse_qft(width), target=b["b"])
    return b.finish()
