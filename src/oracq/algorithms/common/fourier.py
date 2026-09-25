"""Forward and inverse Fourier transforms, interface adapters, and modular addition
without a carry register."""

from __future__ import annotations

import math

from oracq.algorithms.input_model.contracts import positive_integer
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits


def qft(width: int) -> Operation:
    """Generate the positive-sign discrete Fourier transform.

    Args:
        width: Target bit width, in the range 1..64.

    Returns:
        Operation: Has only a target register, includes the trailing swaps, and the bit
        order is little endian.

    The matrix elements are exp(2πi*x*y/2**width)/sqrt(2**width)."""
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


def qft_with_work(width: int) -> Operation:
    """QFT adapter preserving the earlier zero-width work interface.

    Args:
        width: Target bit width, in the range 1..64.

    Returns:
        Operation: QFT operation that additionally carries a zero-width work register
        besides target.
    """
    b = Builder("qft_with_work_" + str(width), {"target": Bits(width), "work": Bits(0)})
    b.call(qft(width), target=b["target"])
    return b.finish()


def inverse_qft(width: int) -> Operation:
    """Generate the adjoint of the QFT, keeping the module call.

    Args:
        width: Target bit width, in the range 1..64.

    Returns:
        Operation: Inverse QFT operation wrapped in an adjoint module call.
    """
    b = Builder("inverse_qft_" + str(width), {"target": Bits(width)})
    with b.adjoint():
        b.call(qft(width), target=b["target"])
    return b.finish()


def fourier_add(width: int) -> Operation:
    """Modular addition without a carry register, implemented with the QFT.

    Args:
        width: Common bit width of a and b.

    Returns:
        Operation: ``|a,b> -> |a,(a+b) mod 2**width>``. a is preserved and there is no
        additional public workspace.

    The Fourier transform and its inverse are kept as module calls."""
    positive_integer(width, "fourier_add.width", maximum=64)
    b = Builder("fourier_add_" + str(width), {"a": Bits(width), "b": Bits(width)})
    b.call(qft(width), target=b["b"])
    for i in range(width):
        for j in range(width - i):
            with b.control(b["a"][i]):
                b.gate("phase", b["b"][j], 2 * math.pi / (1 << (width - i - j)))
    b.call(inverse_qft(width), target=b["b"])
    return b.finish()
