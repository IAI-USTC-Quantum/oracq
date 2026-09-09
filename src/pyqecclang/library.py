"""块编码组合库；alpha 作为 RIR 模块属性保存，精度不属于语言核心。"""

from __future__ import annotations

import cmath
import hashlib
import math
from dataclasses import dataclass, replace
from typing import Protocol

from .builder import Builder, Operation
from .ir import Bits, ValidationError
from .serialization import dumps


class Generator(Protocol):
    def __call__(self, *args, **kwargs) -> Operation: ...


def _name(kind, *values):
    serialized = [dumps(v.program()) if isinstance(v, Operation) else repr(v) for v in values]
    return kind + "_" + hashlib.sha256("\n".join(serialized).encode()).hexdigest()[:20]


@dataclass(frozen=True)
class BlockEncoding:
    operation: Operation

    def __post_init__(self):
        self.operation.program()
        registers = {r.name: r.type for r in self.operation.module.registers}
        if set(registers) != {"target", "signal"} or any(
            t.kind != "bits" for t in registers.values()
        ):
            raise ValidationError("BE 必须具有 target 和 signal 两个 bits 接口")
        if not registers["target"].width:
            raise ValidationError("BE 目标不能为空")
        alpha = dict(self.operation.module.attributes).get("be_alpha")
        if type(alpha) not in (int, float) or not math.isfinite(alpha) or alpha <= 0:
            raise ValidationError("BE 必须声明有限正数 be_alpha")

    @property
    def alpha(self):
        return dict(self.operation.module.attributes)["be_alpha"]

    @property
    def width(self):
        return next(r.type.width for r in self.operation.module.registers if r.name == "target")

    @property
    def signal_qubits(self):
        return next(r.type.width for r in self.operation.module.registers if r.name == "signal")


def block_encoding(operation: Operation, alpha: float = 1.0) -> BlockEncoding:
    attributes = dict(operation.module.attributes)
    attributes["be_alpha"] = alpha
    module = replace(
        operation.module,
        name=_name("be", operation, alpha),
        attributes=tuple(sorted(attributes.items())),
    )
    return BlockEncoding(Operation(module, operation.dependencies))


def identity(width: int) -> BlockEncoding:
    b = Builder(f"identity_{width}", {"target": Bits(width), "signal": Bits(0)})
    return block_encoding(b.finish())


def pauli_x(width: int) -> BlockEncoding:
    b = Builder(f"pauli_x_{width}", {"target": Bits(width), "signal": Bits(0)})
    b.x(b["target"])
    return block_encoding(b.finish())


def zero(width: int) -> BlockEncoding:
    b = Builder(f"zero_{width}", {"target": Bits(width), "signal": Bits(1)})
    b.x(b["signal"])
    return block_encoding(b.finish())


def _resources(a, b=None):
    result = {}
    for prefix, operand in (("a", a), ("b", b)):
        if operand is not None:
            for resource in operand.operation.module.resources:
                result[prefix + "__" + resource.name] = resource.type
    return result


def _call(builder, operand, target, signal, prefix):
    resources = {r.name: prefix + "__" + r.name for r in operand.operation.module.resources}
    builder.call(operand.operation, target=target, signal=signal, resources=resources)


def product(a: BlockEncoding, b: BlockEncoding) -> BlockEncoding:
    if a.width != b.width:
        raise ValidationError("BE 乘积的目标宽度不同")
    builder = Builder(
        _name("product", a.operation, b.operation),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits + b.signal_qubits)},
        _resources(a, b),
    )
    signal = builder["signal"]
    _call(builder, b, builder["target"], signal[a.signal_qubits :], "b")
    _call(builder, a, builder["target"], signal[: a.signal_qubits], "a")
    return block_encoding(builder.finish(), a.alpha * b.alpha)


def scale(coefficient: complex, a: BlockEncoding) -> BlockEncoding:
    if not (math.isfinite(coefficient.real) and math.isfinite(coefficient.imag)):
        raise ValidationError("BE 系数必须有限")
    if coefficient == 0:
        return zero(a.width)
    builder = Builder(
        _name("scale", a.operation, coefficient),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits)},
        _resources(a),
    )
    builder.global_phase(cmath.phase(coefficient))
    _call(builder, a, builder["target"], builder["signal"], "a")
    return block_encoding(builder.finish(), abs(coefficient) * a.alpha)


def linear_combination(
    ca: complex, a: BlockEncoding, cb: complex, b: BlockEncoding
) -> BlockEncoding:
    if a.width != b.width:
        raise ValidationError("BE 求和的目标宽度不同")
    if not all(math.isfinite(x) for c in (ca, cb) for x in (c.real, c.imag)):
        raise ValidationError("BE 系数必须有限")
    if ca == 0:
        return scale(cb, b)
    if cb == 0:
        return scale(ca, a)
    alpha = abs(ca) * a.alpha + abs(cb) * b.alpha
    theta = 2 * math.acos(math.sqrt(abs(ca) * a.alpha / alpha))
    builder = Builder(
        _name("lcu", a.operation, b.operation, ca, cb),
        {"target": Bits(a.width), "signal": Bits(1 + a.signal_qubits + b.signal_qubits)},
        _resources(a, b),
    )
    signal = builder["signal"]
    select, sa, sb = signal[:1], signal[1 : 1 + a.signal_qubits], signal[1 + a.signal_qubits :]
    builder.ry(select, theta)
    with builder.control(select, 0):
        builder.global_phase(cmath.phase(ca))
        _call(builder, a, builder["target"], sa, "a")
    with builder.control(select, 1):
        builder.global_phase(cmath.phase(cb))
        _call(builder, b, builder["target"], sb, "b")
    builder.ry(select, -theta)
    return block_encoding(builder.finish(), alpha)
