"周期网格上的结构化 PDE 端口：移位 LCU、分量选择和同点收缩。"

import cmath
import math

from pyqecclang.algorithms.block_encoding import lcu
from pyqecclang.algorithms.operators import BlockEncoding, _name, identity, product, scale, zero
from pyqecclang.algorithms.oracles import (
    abstract_database,
    annotate,
    diagonal_block_encoding,
    gate_state_prep,
    invoke,
    resources_for,
)
from pyqecclang.algorithms.qham import PortBinding, QHAMBindings
from pyqecclang.applications.qham.reference import centered_weights
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError, fuse


def derivative_encoding(grid, derivative):
    if grid.boundary != "periodic" or any(n & (n - 1) for n in grid.shape):
        raise ValidationError("结构化移位端口需要各轴为二次幂的周期网格；其他边界可提供自己的 BE")
    width = grid.spatial_width
    result = identity(width)
    cursor = 0
    for axis, length, spacing in zip(grid.axes, grid.shape, grid.spacing, strict=True):
        axis_width = (length - 1).bit_length()
        order = dict(derivative).get(axis, 0)
        if order:
            grouped = {}
            for offset, coefficient in centered_weights(order):
                step = (-offset) % length
                grouped[step] = grouped.get(step, 0) + coefficient / spacing**order
            terms = []
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


def coefficient_encoding(discretization, monomial, *, max_words=4096):
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


def qram_coefficient_encoding(discretization, monomial, *, angle_width=8, max_words=4096):
    """已知系数对角线的开放角数据库编码；同一份程序可绑定 gate 或 QRAM 数据库。

    与 coefficient_encoding 的合同一致（target 为空间位、alpha 相同），但系数数据
    不烧进门里：对角值 alpha*cos(theta_a/2) 的角度字留在开放的 XOR 数据库槽位中，
    运行时表由 qram_coefficient_memory 单独计算。仅接受实系数数据；
    角度量化引入不超过 alpha*pi/2**angle_width 的幅值误差。
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


def qram_coefficient_memory(discretization, monomial, *, angle_width=8):
    """与 qram_coefficient_encoding 对应的运行时角表（地址 -> 角度字）。"""
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


def term_encoding(discretization, term, *, max_coefficient_words=4096, coefficient_encoder=coefficient_encoding):
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
    signal_views = []
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
        component_refs = []
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


def structured_fd_bindings(discretization, initial, *, max_coefficient_words=4096, coefficient_encoder=coefficient_encoding):
    """基本矩阵从移位与收缩生成，不物化 N^r × N^r 的端口矩阵。"""
    if len(initial) != discretization.dimension:
        raise ValidationError("初值需要完整寄存器布局")
    ports = []
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
