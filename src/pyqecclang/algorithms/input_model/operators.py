"块编码组合库；alpha 作为 RIR 模块属性保存，精度不属于语言核心。"

from __future__ import annotations

import cmath
import hashlib
import math
from dataclasses import dataclass, replace
from typing import Protocol, cast

from pyqecclang.algorithms.input_model.contracts import OracleView, validate_signature
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.ir import QRAM, Bits, Ref, ValidationError
from pyqecclang.infrastructure.serialization import dumps


class Generator(Protocol):
    """量子操作生成器的结构协议。

    任何以任意位置与关键字参数调用并返回 ``Operation`` 的可调用对象都在结构上
    满足本协议；算法侧据此接受操作工厂而不绑定具体的生成签名。
    """

    def __call__(self, *args: object, **kwargs: object) -> Operation:
        """调用生成器并返回其产出的操作。"""
        ...


def _name(kind: str, *values: object) -> str:
    """按各值的序列化内容生成 ``kind_`` 前缀的确定性模块名。"""
    serialized = [dumps(v.program()) if isinstance(v, Operation) else repr(v) for v in values]
    return kind + "_" + hashlib.sha256("\n".join(serialized).encode()).hexdigest()[:20]


@dataclass(frozen=True)
class BlockEncoding(OracleView):
    """以零信号投影角块表示线性算子的视图。

    包装的 ``Operation`` 恰含 ``target`` 与 ``signal`` 两个 bits 寄存器，并以
    模块属性 ``be_alpha`` 声明有限正的归一化常数：若酉 ``U`` 的零信号角块满足
    ``<0|U|0> = A/alpha``，则本视图把 ``U`` 当作算子 ``A`` 的块编码。alpha 只在
    生成期参与组合代数，执行器不会据此缩放量子态。

    Attributes:
        operation: 被包装的 ``Operation``。
    """

    oracle_kind = "block_encoding"
    operation: Operation

    def block_encoding(self) -> BlockEncoding:
        """返回自身；实现 ``BlockEncodingProtocol`` 的视图适配方法。

        Returns:
            BlockEncoding: 该视图自身。
        """
        return self

    def __post_init__(self) -> None:
        """校验 target/signal 签名与有限正 ``be_alpha`` 属性约束。"""
        validate_signature(self.operation, ("target", "signal"), "BlockEncoding")
        registers = {r.name: r.type for r in self.operation.module.registers}
        if set(registers) != {"target", "signal"} or any(
            t.kind != "bits" for t in registers.values()
        ):
            raise ValidationError("BE 必须具有 target 和 signal 两个 bits 接口")
        if not registers["target"].width:
            raise ValidationError("BE 目标不能为空")
        alpha = dict(self.operation.module.attributes).get("be_alpha")
        if (
            type(alpha) not in (int, float)
            or not math.isfinite(cast("int | float", alpha))
            or cast("int | float", alpha) <= 0
        ):
            raise ValidationError("BE 必须声明有限正数 be_alpha")

    @property
    def alpha(self) -> float:
        """模块属性 ``be_alpha`` 中声明的归一化常数。"""
        return cast("float", dict(self.operation.module.attributes)["be_alpha"])

    @property
    def width(self) -> int:
        """``target`` 寄存器的位宽。"""
        return next(r.type.width for r in self.operation.module.registers if r.name == "target")

    @property
    def signal_qubits(self) -> int:
        """``signal`` 寄存器的位宽。"""
        return next(r.type.width for r in self.operation.module.registers if r.name == "signal")


def block_encoding(operation: Operation, alpha: float = 1.0) -> BlockEncoding:
    """把 ``target``/``signal`` 签名的操作包装为块编码。

    在模块属性中写入 ``be_alpha`` 与 ``oracle_paradigm``，并按内容确定性重命名
    模块；签名与 alpha 约束由 ``BlockEncoding`` 的构造检查完成。

    Args:
        operation: 恰含 ``target`` 与 ``signal`` bits 寄存器的 ``Operation``。
        alpha: 有限正的归一化常数，缺省为一。

    Returns:
        BlockEncoding: 补全属性后的块编码视图。

    Raises:
        ValidationError: 操作不符合块编码的签名或 alpha 约束。
    """
    attributes = dict(operation.module.attributes)
    attributes["be_alpha"] = alpha
    attributes["oracle_paradigm"] = "block_encoding"
    module = replace(
        operation.module,
        name=_name("be", operation, alpha),
        attributes=tuple(sorted(attributes.items())),
    )
    return BlockEncoding(Operation(module, operation.dependencies))


def identity(width: int) -> BlockEncoding:
    """构造单位算子的块编码；alpha 为一且不需要信号位。

    Args:
        width: target 寄存器位宽。

    Returns:
        BlockEncoding: 单位算子的块编码。
    """
    b = Builder(f"identity_{width}", {"target": Bits(width), "signal": Bits(0)})
    return block_encoding(b.finish())


def pauli_x(width: int) -> BlockEncoding:
    """构造 ``width`` 个 X 门张量幂的块编码；alpha 为一且无信号位。

    Args:
        width: target 寄存器位宽。

    Returns:
        BlockEncoding: ``X^⊗width`` 的块编码。
    """
    b = Builder(f"pauli_x_{width}", {"target": Bits(width), "signal": Bits(0)})
    b.x(b["target"])
    return block_encoding(b.finish())


def zero(width: int) -> BlockEncoding:
    """构造零算子的块编码；用一个被翻转的信号位使零信号角块恒为零，alpha 为一。

    Args:
        width: target 寄存器位宽。

    Returns:
        BlockEncoding: 零算子的块编码。
    """
    b = Builder(f"zero_{width}", {"target": Bits(width), "signal": Bits(1)})
    b.x(b["signal"])
    return block_encoding(b.finish())


def _resources(a: BlockEncoding, b: BlockEncoding | None = None) -> dict[str, QRAM]:
    """汇总各操作声明的资源并按 ``a__``/``b__`` 前缀重命名。"""
    result: dict[str, QRAM] = {}
    for prefix, operand in (("a", a), ("b", b)):
        if operand is not None:
            for resource in operand.operation.module.resources:
                result[prefix + "__" + resource.name] = resource.type
    return result


def _call(
    builder: Builder, operand: BlockEncoding, target: Ref, signal: Ref, prefix: str
) -> None:
    """按前缀映射资源后在 ``builder`` 中调用块编码操作。"""
    resources = {r.name: prefix + "__" + r.name for r in operand.operation.module.resources}
    builder.call(operand.operation, target=target, signal=signal, resources=resources)


def product(a: BlockEncoding, b: BlockEncoding) -> BlockEncoding:
    """组合两个块编码的矩阵乘积 ``A·B``。

    两个操作依次作用于共享的 ``target``（先 ``b`` 后 ``a``），信号位按 ``a`` 在
    高位拼接；返回块编码的 alpha 为 ``a.alpha * b.alpha``，信号位数为两者之和。

    Args:
        a: 左因子块编码。
        b: 右因子块编码。

    Returns:
        BlockEncoding: 编码 ``A·B`` 的块编码，两个操作以模块调用保留。

    Raises:
        ValidationError: 两个目标宽度不同。
    """
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
    """用复系数缩放块编码所表示的算子。

    系数相位以 ``global_phase`` 记账，返回块编码的 alpha 为
    ``abs(coefficient) * a.alpha``；系数为零时直接返回 ``zero(a.width)``。

    Args:
        coefficient: 有限复系数。
        a: 被缩放的块编码。

    Returns:
        BlockEncoding: 编码 ``coefficient * A`` 的块编码。

    Raises:
        ValidationError: 系数的实部或虚部不有限。
    """
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
    """组合两个块编码的线性组合 ``ca*A + cb*B``。

    单个选择位按 ``abs(ca)*a.alpha`` 与 ``abs(cb)*b.alpha`` 的权重比例分支，
    两路各经 ``global_phase`` 补偿系数相位后作用于共享的 ``target``；返回块编码
    的 alpha 为两个权重之和，信号位在两操作信号之外多出一个选择位。某项系数
    为零时退化为 ``scale``。

    Args:
        ca: ``a`` 的有限复系数。
        a: 第一个块编码。
        cb: ``b`` 的有限复系数。
        b: 第二个块编码。

    Returns:
        BlockEncoding: 编码 ``ca*A + cb*B`` 的块编码。

    Raises:
        ValidationError: 目标宽度不同，或任一系数含非有限分量。
    """
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
