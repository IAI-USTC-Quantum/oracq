"""小规模模乘、量子求阶与经典因子后处理。"""

import math
from fractions import Fraction

from pyqecclang.algorithms.contracts import positive_integer
from pyqecclang.algorithms.estimation import phase_estimation
from pyqecclang.algorithms.operators import _name
from pyqecclang.algorithms.oracles import _transposition, invoke, resources_for
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError


def modular_multiply(multiplier, modulus, *, width=None, max_width=8):
    """生成可逆的有限规模模乘置换。

    Args:
        multiplier: 与 modulus 互素的正整数。
        modulus: 大于等于二的模数。
        width: 目标位宽；省略时取能够容纳模数的最小位宽。
        max_width: 置换合成预算，默认 8，最多允许 12。

    Returns:
        Operation: x<modulus 时映射到 multiplier*x mod modulus，其余基态保持不变。

    实现枚举有限置换，不代表可扩展的 Shor 模算术。"""
    positive_integer(modulus, "modular_multiply.modulus", minimum=2)
    positive_integer(multiplier, "modular_multiply.multiplier", minimum=1)
    positive_integer(max_width, "modular_multiply.max_width", maximum=12)
    width = (modulus - 1).bit_length() if width is None else width
    positive_integer(width, "modular_multiply.width", maximum=max_width)
    if modulus > (1 << width) or math.gcd(multiplier, modulus) != 1:
        raise ValidationError("模数必须适合寄存器，且乘数必须与模数互素")
    multiplier %= modulus
    permutation = [multiplier * x % modulus if x < modulus else x for x in range(1 << width)]
    b = Builder(
        _name("modular_multiply", multiplier, modulus, width),
        {"target": Bits(width)},
        attributes={
            "algorithm": "modular_multiply",
            "implementation_scope": "bounded permutation synthesis",
            "modulus": modulus,
        },
    )
    visited = set()
    for start in range(1 << width):
        if start in visited:
            continue
        cycle, current = [], start
        while current not in visited:
            visited.add(current)
            cycle.append(current)
            current = permutation[current]
        for other in cycle[1:]:
            _transposition(b, b["target"], start, other)
    return b.finish()


def order_finding(multiplier, modulus, *, precision=3, max_width=8):
    """从整数一开始，对模乘 unitary 做量子求阶。

    Args:
        multiplier: 与 modulus 互素的乘数。
        modulus: 模数。
        precision: QPE 相位位宽。
        max_width: 底层模乘置换的位宽预算。

    Returns:
        Operation: target/phase 接口。phase 提供阶的分数信息，需要经典后处理。"""
    operation = modular_multiply(multiplier, modulus, max_width=max_width)
    n = operation.module.registers[0].type.width
    qpe = phase_estimation(operation, precision=precision)
    b = Builder(
        _name("order_finding", multiplier, modulus, precision),
        {"target": Bits(n), "phase": Bits(precision)},
        resources_for(("qpe", qpe)),
        attributes={
            "algorithm": "order_finding",
            "modulus": modulus,
            "multiplier": multiplier,
            "readout_register": "phase",
        },
    )
    b.x(b["target"][0])
    invoke(b, qpe, "qpe", target=b["target"], phase=b["phase"])
    return b.finish()


def factors_from_phase(value, precision, multiplier, modulus):
    """用相位样本的连分数候选阶尝试得到非平凡因子。

    Args:
        value: 相位寄存器的整数读出值。
        precision: 相位寄存器位宽。
        multiplier: 量子求阶使用的乘数。
        modulus: 待处理整数。

    Returns:
        tuple or None: 已验证的因子对，或表示该样本未能给出因子的 None。"""
    positive_integer(precision, "factors.precision", maximum=63)
    positive_integer(value, "factors.phase", minimum=0, maximum=(1 << precision) - 1)
    positive_integer(modulus, "factors.modulus", minimum=3)
    positive_integer(multiplier, "factors.multiplier")
    common = math.gcd(multiplier, modulus)
    if 1 < common < modulus:
        return tuple(sorted((common, modulus // common)))
    if value == 0:
        return None
    order = Fraction(value, 1 << precision).limit_denominator(modulus).denominator
    if order % 2 or pow(multiplier, order, modulus) != 1:
        return None
    half = pow(multiplier, order // 2, modulus)
    for candidate in (math.gcd(half - 1, modulus), math.gcd(half + 1, modulus)):
        if 1 < candidate < modulus:
            return tuple(sorted((candidate, modulus // candidate)))
    return None
