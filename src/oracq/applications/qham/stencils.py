"周期网格上的结构化 PDE 端口：移位 LCU、分量选择和同点收缩。"

from __future__ import annotations

import cmath
import math
from collections.abc import Callable, Sequence

from oracq.algorithms.input_model.block_encoding import lcu
from oracq.algorithms.input_model.operators import (
    BlockEncoding,
    _name,
    identity,
    product,
    scale,
    zero,
)
from oracq.algorithms.input_model.oracles import (
    abstract_database,
    annotate,
    diagonal_block_encoding,
    gate_state_prep,
    invoke,
    resources_for,
)
from oracq.algorithms.input_model.qham import PortBinding, QHAMBindings
from oracq.applications.qham.pde import EquationTerm, Monomial
from oracq.applications.qham.reference import Discretization, Grid, centered_weights
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits, Ref, ValidationError, fuse


def derivative_encoding(grid: Grid, derivative: tuple[tuple[str, int], ...]) -> BlockEncoding:
    """构造周期网格上空间导数的移位 LCU 块编码。

    每个轴的导数按中心差分模板分解为若干循环移位项，项系数为模板权重除以
    ``spacing**order``；各轴算子复合到整个空间地址寄存器上，各算子只作用
    于自己轴的地址位。无导数的轴贡献恒等，模板系数全部相消的轴贡献零算子。

    Args:
        grid: 周期边界、各轴长度均为 2 的幂的 ``Grid``。
        derivative: 形如 ``((axis, order), ...)`` 的导数说明。

    Returns:
        BlockEncoding: target 为空间寄存器（宽 ``grid.spatial_width``）的
        差分算子块编码。

    Raises:
        ValidationError: 网格边界不是周期，或存在非 2 的幂的轴长。
    """
    if grid.boundary != "periodic" or any(n & (n - 1) for n in grid.shape):
        raise ValidationError("结构化移位端口需要各轴为二次幂的周期网格；其他边界可提供自己的 BE")
    width = grid.spatial_width
    result = identity(width)
    cursor = 0
    for axis, length, spacing in zip(grid.axes, grid.shape, grid.spacing, strict=True):
        axis_width = (length - 1).bit_length()
        order = dict(derivative).get(axis, 0)
        if order:
            grouped: dict[int, float] = {}
            for offset, coefficient in centered_weights(order):
                step = (-offset) % length
                grouped[step] = grouped.get(step, 0) + coefficient / spacing**order
            terms: list[tuple[float, BlockEncoding]] = []
            for step, coefficient in grouped.items():
                if not coefficient:
                    continue
                b = Builder(
                    _name("fd_shift", grid, axis, step), {"target": Bits(width), "signal": Bits(0)}
                )
                if axis_width:
                    b.add_const(b["target"][cursor : cursor + axis_width].reinterpret("uint"), step)
                encoded = BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))
                terms.append((coefficient, encoded))
            axis_op = lcu(terms) if terms else zero(width)
            result = product(axis_op, result)
        cursor += axis_width
    return result


def coefficient_encoding(
    discretization: Discretization, monomial: Monomial, *, max_words: int = 4096
) -> BlockEncoding:
    """构造已知系数对角乘子的门实现块编码。

    逐地址求出已知场（含其空间导数）与单项式系数的乘积值，按地址受控旋转
    单个 signal 量子比特并补偿复相位，使 signal 投影回零时实现对角值除以
    ``alpha`` 的作用；值取遍整个空间地址寄存器，网格外的填充地址为 0。

    Args:
        discretization: 提供网格与已知场数据的 ``Discretization``。
        monomial: 待编码的单项式；无已知场且网格无填充位时退化为常数缩放。
        max_words: 门系数表允许的最大地址数。

    Returns:
        BlockEncoding: target 为空间寄存器、``alpha`` 为对角值最大模的
        块编码；对角值全为零时返回零算子。

    Raises:
        ValidationError: 系数表超过 ``max_words`` 预算。
    """
    grid = discretization.grid
    width = grid.spatial_width
    if not monomial.known and grid.size == 1 << width:
        return scale(monomial.coefficient, identity(width))
    if 1 << width > max_words:
        raise ValidationError("已知系数的门实现超过预算；请绑定 QRAM/自定义系数 BE")
    values = [
        discretization.known_product(monomial, row) if row < grid.size else 0j
        for row in range(1 << width)
    ]
    alpha = max((abs(v) for v in values), default=0)
    if not alpha:
        return zero(width)
    b = Builder(_name("known_multiplier", values), {"target": Bits(width), "signal": Bits(1)})
    for address, value in enumerate(values):
        with b.control(b["target"], address):
            b.ry(b["signal"], 2 * math.acos(min(1, abs(value) / alpha)))
            if value:
                b.global_phase(cmath.phase(value))
    return BlockEncoding(
        annotate(
            b.finish(), "block_encoding", be_alpha=alpha, implementation="known_diagonal_multiplier"
        )
    )


def qram_coefficient_encoding(
    discretization: Discretization,
    monomial: Monomial,
    *,
    angle_width: int = 8,
    max_words: int = 4096,
) -> BlockEncoding:
    """已知系数对角线的开放角数据库编码；同一份程序可绑定 gate 或 QRAM 数据库。

    与 coefficient_encoding 的合同一致（target 为空间位、alpha 相同），但系数数据
    不烧进门里：对角值 alpha*cos(theta_a/2) 的角度字留在开放的 XOR 数据库槽位中，
    运行时表由 qram_coefficient_memory 单独计算。仅接受实系数数据；
    角度量化引入不超过 alpha*pi/2**angle_width 的幅值误差。

    Args:
        discretization: 提供网格与已知场数据的 ``Discretization``。
        monomial: 待编码的单项式；无已知场且网格无填充位时退化为常数缩放。
        angle_width: 角度字的位宽，决定对角值的量化精度。
        max_words: 角数据库允许的最大地址数。

    Returns:
        BlockEncoding: 对角值留在开放 XOR 数据库槽位中的块编码，``alpha``
        为对角值的最大模，对角值全为零时返回零算子。
    """
    grid = discretization.grid
    width = grid.spatial_width
    if not monomial.known and grid.size == 1 << width:
        return scale(monomial.coefficient, identity(width))
    if 1 << width > max_words:
        raise ValidationError("已知系数的 QRAM 角表超过预算；请提供自定义系数 BE")
    values = [
        discretization.known_product(monomial, row) if row < grid.size else 0j
        for row in range(1 << width)
    ]
    if any(v.imag for v in values):
        raise ValidationError("QRAM 角编码的系数数据必须为实数")
    alpha = max((abs(v.real) for v in values), default=0)
    if not alpha:
        return zero(width)
    db = abstract_database(_name("coefficient_angles", values), width, angle_width)
    return diagonal_block_encoding(db, alpha=alpha)


def qram_coefficient_memory(
    discretization: Discretization, monomial: Monomial, *, angle_width: int = 8
) -> dict[int, int]:
    """与 qram_coefficient_encoding 对应的运行时角表（地址 -> 角度字）。

    Args:
        discretization: 提供网格与已知场数据的 ``Discretization``。
        monomial: 编码时使用的同一单项式，须与编码调用保持一致。
        angle_width: 角度字的位宽，须与编码调用保持一致。

    Returns:
        dict[int, int]: 空间地址到角度字的映射；可退化为常数缩放或对角值
        全为零时返回空表。
    """
    grid = discretization.grid
    width = grid.spatial_width
    if not monomial.known and grid.size == 1 << width:
        return {}
    values = [
        discretization.known_product(monomial, row) if row < grid.size else 0j
        for row in range(1 << width)
    ]
    alpha = max((abs(v.real) for v in values), default=0)
    if not alpha:
        return {}
    step = 2 * math.pi / (1 << angle_width)
    return {
        address: round(2 * math.acos(min(1, max(-1, v.real / alpha))) / step)
        % (1 << angle_width)
        for address, v in enumerate(values)
    }


def term_encoding(
    discretization: Discretization,
    term: EquationTerm,
    *,
    max_coefficient_words: int = 4096,
    coefficient_encoder: Callable[..., BlockEncoding] = coefficient_encoding,
) -> BlockEncoding:
    """为单个 PDE 方程项构造结构化差分端口的多线性块编码。

    端口由各因子字段的中心差分导数、已知系数对角乘子、分量选择与外导数
    组装而成：因子各自在空间位上差分，多因子情形先收缩到同一点再对角
    合并；signal 中的两个拒绝位分别标记分量不符与（多因子时）各因子不在
    同一点的输入坐标，零阶项则标记非零地址并在空间位上取均匀叠加。
    不物化端口矩阵。

    Args:
        discretization: 提供网格、分量布局与已知数据的 ``Discretization``。
        term: 含输出分量与一个单项式的 ``EquationTerm``。
        max_coefficient_words: 传给系数编码器的地址数预算。
        coefficient_encoder: 系数对角乘子的编码函数，合同同
            ``coefficient_encoding``（如 ``qram_coefficient_encoding``）。

    Returns:
        BlockEncoding: target 宽 ``max(1, arity)*discretization.width`` 的
        矩形块编码，``rectangular_arity`` 属性记录元数；``alpha`` 为系数、
        外导数与各因子导数的 alpha 之积，零阶项另乘 ``sqrt(2**spatial_width)``。

    Raises:
        ValidationError: 端口 target 宽超过 64，或系数编码器报告超预算。
    """
    n = discretization.width
    ns = discretization.grid.spatial_width
    nc = discretization.component_width
    monomial = term.monomial
    arity = len(monomial.fields)
    width = max(1, arity) * n
    if width > 64:
        raise ValidationError("单个多线性端口超过当前 BE target 包装宽度")
    derivatives = [derivative_encoding(discretization.grid, a.derivative) for a in monomial.fields]
    outer = derivative_encoding(discretization.grid, monomial.outer_derivative)
    coefficient = coefficient_encoder(discretization, monomial, max_words=max_coefficient_words)
    operands = [(f"d{i}", op.operation) for i, op in enumerate(derivatives)]
    operands += [("coefficient", coefficient.operation), ("outer", outer.operation)]
    signal_width = (
        sum(op.signal_qubits for op in derivatives)
        + coefficient.signal_qubits
        + outer.signal_qubits
        + 2
    )
    b = Builder(
        _name(
            "pde_multilinear_term",
            term,
            discretization.grid,
            discretization.pde.fields,
            coefficient.operation,
            outer.operation,
        ),
        {"target": Bits(width), "signal": Bits(signal_width)},
        resources_for(*operands),
        attributes={
            "pde_arity": arity,
            "input_coordinate_order": "factor 0 low",
            "structured_stencil": True,
        },
    )
    cursor = 0
    signal_views: list[Ref] = []
    for op in derivatives:
        signal_views.append(b["signal"][cursor : cursor + op.signal_qubits])
        cursor += op.signal_qubits
    coeff_signal = b["signal"][cursor : cursor + coefficient.signal_qubits]
    cursor += coefficient.signal_qubits
    outer_signal = b["signal"][cursor : cursor + outer.signal_qubits]
    cursor += outer.signal_qubits
    reject_input, reject_diagonal = b["signal"][cursor], b["signal"][cursor + 1]
    if arity == 0:
        b.x(reject_input)
        with b.control(b["target"], 0):
            b.x(reject_input)
        b.h(b["target"][:ns])
    else:
        component_refs: list[Ref] = []
        for i, atom in enumerate(monomial.fields):
            field = discretization.pde.fields.index(atom.name)
            component = b["target"][i * n + ns : (i + 1) * n]
            for bit in range(nc):
                if (field >> bit) & 1:
                    b.x(component[bit])
            component_refs.append(component)
        if nc:
            b.x(reject_input)
            with b.control(fuse(*component_refs), 0):
                b.x(reject_input)
        for i, op in enumerate(derivatives):
            invoke(
                b,
                op.operation,
                f"d{i}",
                target=b["target"][i * n : i * n + ns],
                signal=signal_views[i],
            )
    invoke(b, coefficient.operation, "coefficient", target=b["target"][:ns], signal=coeff_signal)
    if arity > 1:
        for i in range(1, arity):
            b.xor(b["target"][:ns], b["target"][i * n : i * n + ns])
        b.x(reject_diagonal)
        with b.control(b["target"][n:], 0):
            b.x(reject_diagonal)
    output = discretization.pde.fields.index(term.output)
    for bit in range(nc):
        if (output >> bit) & 1:
            b.x(b["target"][ns + bit])
    invoke(b, outer.operation, "outer", target=b["target"][:ns], signal=outer_signal)
    alpha = coefficient.alpha * outer.alpha * math.prod(op.alpha for op in derivatives)
    if arity == 0:
        alpha *= math.sqrt(1 << ns)
    return BlockEncoding(
        annotate(
            b.finish(),
            "block_encoding",
            be_alpha=alpha,
            rectangular_arity=arity,
            implementation="shifts_components_diagonal_contraction",
        )
    )


def structured_fd_bindings(
    discretization: Discretization,
    initial: Sequence[complex],
    *,
    max_coefficient_words: int = 4096,
    coefficient_encoder: Callable[..., BlockEncoding] = coefficient_encoding,
) -> QHAMBindings:
    """基本矩阵从移位与收缩生成，不物化 N^r × N^r 的端口矩阵。

    Args:
        discretization: 提供网格、分量布局与已知数据的 ``Discretization``。
        initial: 长度为 ``discretization.dimension`` 的初值向量；范数为零时
            改用第一个基矢制备。
        max_coefficient_words: 传给系数编码器的地址数预算。
        coefficient_encoder: 系数对角乘子的编码函数，合同同
            ``coefficient_encoding``（如 ``qram_coefficient_encoding``）。

    Returns:
        QHAMBindings: 各端口绑定到移位 LCU 块编码、并含初值制备与范数的
        QHAM 绑定集合。
    """
    if len(initial) != discretization.dimension:
        raise ValidationError("初值需要完整寄存器布局")
    ports: list[tuple[str, PortBinding]] = []
    for port in discretization.pde.ports:
        encoded = lcu(
            [
                (
                    1,
                    term_encoding(
                        discretization,
                        term,
                        max_coefficient_words=max_coefficient_words,
                        coefficient_encoder=coefficient_encoder,
                    ),
                )
                for term in port.terms
            ]
        )
        ports.append((port.name, PortBinding(encoded, port.arity)))
    norm = math.sqrt(sum(abs(v) ** 2 for v in initial))
    prep = (
        gate_state_prep(initial)
        if norm
        else gate_state_prep([1.0] + [0.0] * (discretization.dimension - 1))
    )
    return QHAMBindings(discretization.width, tuple(ports), prep, norm)
