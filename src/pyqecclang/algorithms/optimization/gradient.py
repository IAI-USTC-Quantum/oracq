"""Jordan 量子梯度估计：相位 oracle 输入模型与单查询梯度读出。

实现 Jordan 2005（PRL 95, 050501）的梯度估计。网格约定：第 i 个坐标寄存器占
位段 [i*grid_bits, (i+1)*grid_bits)，整数值 k 对应定点网格点 x = k/N，N = 2**grid_bits。
缩放约定与 Gilyén–Arunachalam–Wiebe 2019 的推广一致：相位 oracle 实现
``O|x> = exp(2πi·N·f(x))|x>`` ，即 phase_scale 必须等于 N。f 在网格上近似线性时，
单次 oracle 调用加各坐标逆 QFT 即把 N·∂f/∂x_i 写入第 i 个坐标寄存器；
一次查询得到全部 d 个分量，经典确定性评估同一梯度需要 O(d) 次函数查询。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from pyqecclang.algorithms.common.arithmetic import FixedFormat
from pyqecclang.algorithms.common.fourier import inverse_qft
from pyqecclang.algorithms.input_model.contracts import (
    OracleView,
    fail,
    finite_real,
    positive_integer,
    validate_signature,
)
from pyqecclang.algorithms.input_model.operators import _name
from pyqecclang.algorithms.input_model.oracles import annotate, declare, invoke, resources_for
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.ir import Bits, ValidationError


@dataclass(frozen=True)
class PhaseOracle(OracleView):
    """网格相位 oracle 视图：``O|x> = exp(2πi·phase_scale·f(x))|x>`` 的对角相位作用。

    target 的整数值编码 d 维定点网格点；缩放因子记入 phase_scale 属性，
    供算法侧在生成阶段核对缩放约定。"""

    oracle_kind = "phase_oracle"
    operation: Operation

    def phase_oracle(self):
        """返回相位 oracle 视图；本类自身即包装相位 oracle，直接返回自身。"""
        return self

    def __post_init__(self):
        validate_signature(self.operation, ("target",), "PhaseOracle")
        scale = dict(self.operation.module.attributes).get("phase_scale")
        if scale is None:
            raise ValidationError("相位 oracle 必须声明 phase_scale 属性")
        finite_real(scale, "PhaseOracle.phase_scale", minimum=0, strict=True)

    @property
    def width(self):
        """target 寄存器的位宽，即相位 oracle 作用的网格寄存器总宽度。"""
        return next(r.type.width for r in self.operation.module.registers if r.name == "target")

    @property
    def phase_scale(self):
        """oracle 声明的相位缩放因子，即 ``O|x> = exp(2πi·phase_scale·f(x))|x>`` 中的缩放。"""
        return dict(self.operation.module.attributes)["phase_scale"]


def abstract_phase_oracle(name, width, *, phase_scale):
    """开放声明一个相位 oracle；相位缩放写入 phase_scale 属性，实现留待 bind。"""
    positive_integer(width, "abstract_phase_oracle.width", maximum=64)
    finite_real(phase_scale, "abstract_phase_oracle.phase_scale", minimum=0, strict=True)
    return PhaseOracle(
        declare(
            name,
            {"target": Bits(width)},
            paradigm="phase_oracle",
            attributes={"phase_scale": phase_scale},
        )
    )


def gate_phase_oracle(width, angles, *, phase_scale, name=None):
    """由显式相位表构造对角相位 oracle，angles[x] 是基态 ``|x>`` 获得的相位弧度。

    小尺度见证可直接枚举全部 2**width 个相位；大网格应改用
    function_phase_oracle 或绑定其他实现。"""
    positive_integer(width, "gate_phase_oracle.width", maximum=64)
    finite_real(phase_scale, "gate_phase_oracle.phase_scale", minimum=0, strict=True)
    angles = tuple(angles)
    if len(angles) != 1 << width:
        raise ValidationError("相位表长度必须等于 2**width")
    for angle in angles:
        finite_real(angle, "gate_phase_oracle.angles")
    b = Builder(
        name or _name("phase_table", width, angles, phase_scale), {"target": Bits(width)}
    )
    for value, angle in enumerate(angles):
        reduced = math.remainder(angle, 2 * math.pi)
        if reduced:
            with b.control(b["target"], value):
                b.global_phase(reduced)
    return PhaseOracle(
        annotate(
            b.finish(),
            "phase_oracle",
            phase_scale=phase_scale,
            implementation="diagonal_phase_table",
        )
    )


def function_phase_oracle(
    source,
    *,
    dimension,
    grid_bits,
    fmt=None,
    scale=None,
    name=None,
    constants=None,
    helpers=None,
    config=None,
    max_unroll=128,
    entry=None,
):
    """由 mathfunc 算术构造相位 oracle：计算定点 f(x) 后对输出做相位踢回并复原。

    Args:
        source: 纯 Python 函数源码，签名 f(x0, ..., x{d-1})，返回单个实数。
        dimension: 网格维数 d，范围为 1..16。
        grid_bits: 每个坐标的位数 m，网格点 x = k/2**m 必须可被 fmt 精确表示。
        fmt: 定点格式，默认 FixedFormat(grid_bits+12, grid_bits+8)；必须带符号且 fraction >= grid_bits。
        scale: 相位缩放，默认 2**grid_bits，即 Jordan 缩放约定。
        name: 覆盖自动生成的模块名。
        constants: 传给 mathfunc 前端的常量绑定。
        helpers: 传给 mathfunc 前端的辅助函数。
        config: 传给 mathfunc 前端的 MathConfig。
        max_unroll: 传给 mathfunc 前端的展开上限。
        entry: 传给 mathfunc 前端的入口函数名。

    Returns:
        PhaseOracle: 相位为 ``exp(2πi·scale·f(x))`` ；f 超出 fmt 值域时由 status 标记，数值不保证。

    坐标寄存器经零态工作位写入函数输入，相位踢回按二进制补码解码输出，
    随后逆调用复原全部工作位。"""
    positive_integer(dimension, "function_phase_oracle.dimension", maximum=16)
    positive_integer(grid_bits, "function_phase_oracle.grid_bits", maximum=32)
    if dimension * grid_bits > 64:
        raise ValidationError("网格总宽度不能超过 64 位")
    fmt = fmt or FixedFormat(grid_bits + 12, grid_bits + 8)
    if not isinstance(fmt, FixedFormat):
        raise ValidationError("fmt 必须是 FixedFormat")
    if not fmt.signed or fmt.fraction < grid_bits:
        raise ValidationError("定点格式需要符号位，且 fraction >= grid_bits 以精确表示网格坐标")
    scale = (1 << grid_bits) if scale is None else scale
    finite_real(scale, "function_phase_oracle.scale", minimum=0, strict=True)
    from pyqecclang.infrastructure.mathfunc import compile_function

    compiled = compile_function(
        source,
        fmt=fmt,
        inputs={f"x{i}": "real" for i in range(dimension)},
        constants=constants,
        helpers=helpers,
        config=config,
        max_unroll=max_unroll,
        entry=entry,
    )
    operation = compiled.operation
    expected = {f"x{i}" for i in range(dimension)} | {"out", "status"}
    if {r.name for r in operation.module.registers} != expected:
        raise ValidationError("函数必须有 dimension 个实数参数并返回单个实数")
    b = Builder(
        name or _name("function_phase", operation, grid_bits, scale),
        {"target": Bits(dimension * grid_bits)},
        resources_for(("f", operation)),
    )
    xin = [b.local(f"x{i}", Bits(fmt.width)) for i in range(dimension)]
    out = b.local("out", Bits(fmt.width))
    status = b.local("status", Bits(2))
    shift = fmt.fraction - grid_bits
    for i in range(dimension):
        b.xor(b["target"][i * grid_bits : (i + 1) * grid_bits], xin[i][shift : shift + grid_bits])
    arguments = {**{f"x{i}": xin[i] for i in range(dimension)}, "out": out, "status": status}
    invoke(b, operation, "f", **arguments)
    for bit in range(fmt.width):
        weight = (
            -math.ldexp(1.0, fmt.width - 1 - fmt.fraction)
            if bit == fmt.width - 1
            else math.ldexp(1.0, bit - fmt.fraction)
        )
        angle = math.remainder(2 * math.pi * scale * weight, 2 * math.pi)
        if angle:
            with b.control(out[bit]):
                b.global_phase(angle)
    with b.adjoint():
        invoke(b, operation, "f", **arguments)
    for i in range(dimension):
        b.xor(b["target"][i * grid_bits : (i + 1) * grid_bits], xin[i][shift : shift + grid_bits])
    return PhaseOracle(
        annotate(
            b.finish(),
            "phase_oracle",
            phase_scale=scale,
            implementation="mathfunc_kickback",
            math_function=dict(operation.module.attributes).get("math_function"),
        )
    )


def gradient_estimation(oracle, *, dimension, grid_bits):
    """生成 Jordan 梯度估计电路。

    Args:
        oracle: PhaseOracle 或具 target 签名的相位 oracle 操作，实现 ``O|x> = exp(2πi·N·f(x))|x>`` 。
        dimension: 网格维数 d，范围为 1..16。
        grid_bits: 每个坐标寄存器的位数 m，N = 2**m，范围为 1..32。

    Returns:
        Operation: 寄存器 target，宽 d*m。读出后用 gradient_from_readout 解码各梯度分量。

    Raises:
        ValidationError: 网格参数无效、oracle 宽度不符或 phase_scale 不等于 2**grid_bits。

    各坐标寄存器制备均匀叠加，单次调用相位 oracle，再逐坐标施加逆 QFT。
    f 近似线性时相位把 N·∂f/∂x_i 写入坐标寄存器 i 的 Fourier 基，
    逆 QFT 后读出即各分量的定点近似。一次查询得到全部 d 个分量，
    经典确定性梯度评估需要 O(d) 次函数查询（Jordan 2005, PRL 95, 050501）。"""
    positive_integer(dimension, "gradient.dimension", maximum=16)
    positive_integer(grid_bits, "gradient.grid_bits", maximum=32)
    width = dimension * grid_bits
    if width > 64:
        raise ValidationError("网格总宽度不能超过 64 位")
    view = oracle if isinstance(oracle, PhaseOracle) else PhaseOracle(oracle)
    if view.width != width:
        raise ValidationError("相位 oracle 宽度必须等于 dimension*grid_bits")
    grid_points = 1 << grid_bits
    if view.phase_scale != grid_points:
        fail(
            "INPUT_PROMISE",
            "gradient.oracle.phase_scale",
            grid_points,
            view.phase_scale,
            "缩放约定要求 phase_scale == 2**grid_bits",
        )
    b = Builder(
        _name("jordan_gradient", view.operation, dimension, grid_bits),
        {"target": Bits(width)},
        resources_for(("oracle", view.operation)),
        attributes={
            "algorithm": "jordan_gradient",
            "readout_register": "target",
            "dimension": dimension,
            "grid_bits": grid_bits,
            "oracle_queries": 1,
            "classical_queries": "O(dimension)",
            "decoder": "gradient_from_readout",
            "reference": "Jordan 2005, PRL 95, 050501",
        },
    )
    b.h(b["target"])
    invoke(b, view.operation, "oracle", target=b["target"])
    for i in range(dimension):
        invoke(
            b,
            inverse_qft(grid_bits),
            "qft",
            target=b["target"][i * grid_bits : (i + 1) * grid_bits],
        )
    return b.finish()


def gradient_from_readout(value, *, dimension, grid_bits):
    """把 target 的整数读出解码为梯度各分量估计。

    Args:
        value: target 寄存器的整数读出。
        dimension: 网格维数。
        grid_bits: 每个坐标寄存器的位数 m。

    Returns:
        tuple: 第 i 个分量把位段 [i*m, (i+1)*m) 按二进制补码解释后除以 2**m。"""
    positive_integer(dimension, "gradient.dimension", maximum=16)
    positive_integer(grid_bits, "gradient.grid_bits", maximum=32)
    positive_integer(
        value, "gradient.value", minimum=0, maximum=(1 << (dimension * grid_bits)) - 1
    )
    result = []
    for i in range(dimension):
        chunk = (value >> (i * grid_bits)) & ((1 << grid_bits) - 1)
        if chunk >> (grid_bits - 1):
            chunk -= 1 << grid_bits
        result.append(chunk / (1 << grid_bits))
    return tuple(result)
