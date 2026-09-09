"QCL 行规则到模块化寄存器 BE，以及保留尺度的 QODE 输入。"

from __future__ import annotations

import math
from dataclasses import dataclass, replace

from pyqecclang.algorithms.arithmetic import BooleanNetwork
from pyqecclang.algorithms.block_encoding import lcu, matrix_pauli_encoding
from pyqecclang.algorithms.operators import BlockEncoding, _name, identity, product, zero
from pyqecclang.algorithms.oracles import (
    StateOracle,
    StatePreparation,
    abstract_block_encoding,
    abstract_state_prep,
    annotate,
    gate_state_prep,
    invoke,
    resources_for,
)
from pyqecclang.applications.qham.linearization import Block
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError


@dataclass(frozen=True)
class PortBinding:
    encoding: BlockEncoding
    arity: int


@dataclass(frozen=True)
class QHAMBindings:
    state_width: int
    ports: tuple[tuple[str, PortBinding], ...]
    initial: StatePreparation
    initial_norm: float

    def validate(self, plan):
        if not 1 <= self.state_width <= 64 or self.initial.width != self.state_width:
            raise ValidationError("QHAM 状态空间与初值寄存器布局不匹配")
        if not math.isfinite(self.initial_norm) or self.initial_norm < 0:
            raise ValidationError("初始物理向量范数必须为有限非负数")
        supplied = dict(self.ports)
        if len(supplied) != len(self.ports) or set(supplied) != {p.name for p in plan.pde.ports}:
            raise ValidationError("QHAM 输入端口集合不匹配")
        for p in plan.pde.ports:
            binding = supplied[p.name]
            if (
                binding.arity != p.arity
                or binding.encoding.width != max(1, p.arity) * self.state_width
            ):
                raise ValidationError("QHAM 矩形端口的输入/输出宽度不匹配：" + p.name)
        return self

    @classmethod
    def declare(cls, plan, state_width, port_specs, *, initial_norm, initial_work=0, prefix="Qham"):
        if set(port_specs) != {p.name for p in plan.pde.ports}:
            raise ValidationError("开放 QHAM 端口的 alpha/信号位规格不完整")
        ports = []
        for port in plan.pde.ports:
            alpha, signal = port_specs[port.name]
            be = abstract_block_encoding(
                prefix + "_" + port.name, max(1, port.arity) * state_width, signal, alpha
            )
            attrs = dict(be.operation.module.attributes)
            attrs.update(
                qham_pde_contract=plan.pde.dumps(), qham_port=port.name, input_rank=port.arity
            )
            op = replace(
                be.operation,
                module=replace(be.operation.module, attributes=tuple(sorted(attrs.items()))),
            )
            ports.append((port.name, PortBinding(BlockEncoding(op), port.arity)))
        return cls(
            state_width,
            tuple(ports),
            abstract_state_prep(prefix + "_initial", state_width, initial_work),
            initial_norm,
        ).validate(plan)


def gate_bindings(discretization, initial, *, max_port_qubits=5):
    """仅物化小型基础端口，不构造整个提升矩阵 G。"""
    width = discretization.width
    if len(initial) != discretization.dimension or any(
        not math.isfinite(complex(v).real) or not math.isfinite(complex(v).imag) for v in initial
    ):
        raise ValidationError("初值必须是有限的、按完整寄存器布局编码的向量")
    ports = []
    for port in discretization.pde.ports:
        local_width = max(1, port.arity) * width
        if local_width > max_port_qubits:
            raise ValidationError("基础端口门实现超过小规模物化预算；可提供开放/QRAM/结构化 BE")
        size = 1 << local_width
        matrix = [
            [
                discretization.entry(port.name, row, col)
                if row < discretization.dimension and col < discretization.dimension**port.arity
                else 0j
                for col in range(size)
            ]
            for row in range(size)
        ]
        ports.append((port.name, PortBinding(matrix_pauli_encoding(matrix), port.arity)))
    norm = math.sqrt(sum(abs(v) ** 2 for v in initial))
    prep = (
        gate_state_prep(initial)
        if norm
        else gate_state_prep([1.0] + [0.0] * (discretization.dimension - 1))
    )
    return QHAMBindings(width, tuple(ports), prep, norm)


def _permute_slots(builder, target, width, desired):
    current = list(range(len(desired)))
    for position, token in enumerate(desired):
        other = current.index(token)
        if position != other:
            builder.swap(
                target[position * width : (position + 1) * width],
                target[other * width : (other + 1) * width],
            )
            current[position], current[other] = current[other], current[position]


def place_port(binding, state_width, source_rank, position):
    arity = binding.arity
    output_rank = source_rank - arity + 1
    rank = max(source_rank, output_rank)
    if source_rank < arity or not 0 <= position < output_rank:
        raise ValidationError("矩形张量位置无效")
    base = binding.encoding
    b = Builder(
        _name("qham_place", base.operation, state_width, source_rank, position),
        {"target": Bits(rank * state_width), "signal": Bits(base.signal_qubits)},
        resources_for(("port", base.operation)),
        attributes={
            "qham_input_rank": source_rank,
            "qham_output_rank": output_rank,
            "qham_tensor_position": position,
            "rectangular_arity": arity,
        },
    )
    if arity == 0:
        desired = [*range(position), source_rank, *range(position, source_rank)]
        _permute_slots(b, b["target"], state_width, desired)
    start = position * state_width
    invoke(
        b,
        base.operation,
        "port",
        target=b["target"][start : start + max(1, arity) * state_width],
        signal=b["signal"],
    )
    if arity > 1:
        desired = [
            *range(position + 1),
            *range(position + arity, source_rank),
            *range(position + 1, position + arity),
        ]
        _permute_slots(b, b["target"], state_width, desired)
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=base.alpha))


def embed_rectangular(local, global_width, row_offset, column_offset, input_width, output_width):
    """同时约束输入和输出窗口，避免矩形零填充污染相邻块。"""
    if (
        local.width > global_width
        or max(row_offset + (1 << output_width), column_offset + (1 << input_width))
        > 1 << global_width
    ):
        raise ValidationError("矩形块超出全局布局")
    b = Builder(
        _name(
            "qham_embed",
            local.operation,
            global_width,
            row_offset,
            column_offset,
            input_width,
            output_width,
        ),
        {"target": Bits(global_width), "signal": Bits(local.signal_qubits + 2)},
        resources_for(("local", local.operation)),
        attributes={
            "qham_row_offset": row_offset,
            "qham_column_offset": column_offset,
            "qham_input_width": input_width,
            "qham_output_width": output_width,
        },
    )
    b.add_const(b["target"].reinterpret("uint"), (-column_offset) % (1 << global_width))

    def reject_nonzero(ref, flag):
        if ref.width:
            b.x(flag)
            with b.control(ref, 0):
                b.x(flag)

    reject_nonzero(b["target"][input_width:], b["signal"][local.signal_qubits])
    invoke(
        b,
        local.operation,
        "local",
        target=b["target"][: local.width],
        signal=b["signal"][: local.signal_qubits],
    )
    reject_nonzero(b["target"][output_width:], b["signal"][local.signal_qubits + 1])
    b.add_const(b["target"].reinterpret("uint"), row_offset)
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=local.alpha))


def generator_encoding(plan, bindings, eta, *, max_blocks=256, max_terms=4096):
    bindings.validate(plan)
    n = bindings.state_width
    if plan.max_rank * n > 64:
        raise ValidationError("当前 BE target 包装限制为 64 位；惰性 QCL plan 不受此物化限制")
    dimension = 1 << n
    width = (plan.raw_dimension(dimension) - 1).bit_length()
    blocks = tuple(plan.blocks(max_blocks=max_blocks))
    terms = []
    placements = {}
    for block in blocks:
        for coupling in plan.row_terms(block):
            coefficient = coupling.weight.evaluate(eta)
            if not coefficient:
                continue
            if len(terms) >= max_terms:
                raise ValidationError("显式 QCL BE 耦合数超过预算；请保留开放生成元")
            key = (coupling.operator, coupling.column.rank, coupling.position)
            if key not in placements:
                placements[key] = place_port(
                    dict(bindings.ports)[coupling.operator],
                    n,
                    coupling.column.rank,
                    coupling.position,
                )
            embedded = embed_rectangular(
                placements[key],
                width,
                plan.offset(block, dimension),
                plan.offset(coupling.column, dimension),
                coupling.column.rank * n,
                block.rank * n,
            )
            terms.append((coefficient, embedded))
    result = lcu(terms) if terms else zero(width)
    result = BlockEncoding(
        replace(
            result.operation,
            module=replace(
                result.operation.module,
                name=_name("qham_generator", result.operation, plan.dumps(), eta),
            ),
        )
    )
    return BlockEncoding(
        annotate(
            result.operation,
            "block_encoding",
            be_alpha=result.alpha,
            algorithm="qham_generated_qcl",
            ham_order=plan.order,
            pde_degree=plan.pde.degree,
            tensor_rank_limit=plan.max_rank,
            raw_dimension=plan.raw_dimension(dimension),
            explicit_couplings=len(terms),
            qham_plan=plan.dumps(),
            correctness="linearization_witnessed; solver pending",
        )
    )


def lifted_initial(plan, bindings):
    n = bindings.state_width
    if plan.max_rank * n > 64:
        raise ValidationError("当前初值 target 包装超过 64 位；请保留开放初值模型")
    dimension = 1 << n
    width = (plan.raw_dimension(dimension) - 1).bit_length()
    candidates = [(Block("physical"), 1.0)]
    candidates += [(Block("tensor", (0,) * k), float(k)) for k in range(1, plan.max_rank + 1)]
    if plan.has_forcing:
        candidates.append((Block("tensor"), 0.0))
    if bindings.initial_norm == 0:
        candidates = [item for item in candidates if item[1] == 0]
        if not candidates:
            raise ValidationError("提升初值为零；无强迫时应在经典侧返回零解")
    log_r = math.log(bindings.initial_norm) if bindings.initial_norm else 0.0
    logs = [k * log_r for _, k in candidates]
    maximum = max(logs)
    scaled = [math.exp(value - maximum) for value in logs]
    norm = math.sqrt(sum(x * x for x in scaled))
    log_norm = maximum + math.log(norm)
    branches = len(candidates)
    selector_width = max(1, (branches - 1).bit_length())
    weights = scaled + [0.0] * ((1 << selector_width) - branches)
    selector_prep = gate_state_prep(weights)
    # 每次已知零输入的制备都复净 work，可顺序复用同一段工作寄存器。
    work_width = bindings.initial.work_width + selector_width
    b = Builder(
        _name("qham_initial", bindings.initial.operation, plan.dumps(), bindings.initial_norm),
        {"target": Bits(width), "work": Bits(work_width)},
        resources_for(("initial", bindings.initial.operation)),
        attributes={
            "algorithm": "qham_tensor_initial",
            "log_initial_norm": log_norm,
            "tensor_preparation": "independent repeated oracle invocations; no cloning",
        },
    )
    selector = b["work"][work_width - selector_width :]
    invoke(b, selector_prep.operation, target=selector, work=selector[:0])
    intervals = []
    for branch, (block, _) in enumerate(candidates):
        rank = block.rank
        with b.control(selector, branch):
            for slot in range(rank):
                invoke(
                    b,
                    bindings.initial.operation,
                    "initial",
                    target=b["target"][slot * n : (slot + 1) * n],
                    work=b["work"][: bindings.initial.work_width],
                )
            offset = plan.offset(block, dimension)
            b.add_const(b["target"].reinterpret("uint"), offset)
        intervals.append((offset, offset + dimension**rank, branch))
    # 块的支持互不重叠，利用地址区间反算选择标签。
    net = BooleanNetwork()
    address = net.input("address", width)
    outputs = [0] * selector_width
    for lo, hi, branch in intervals:
        inside = net.inv(net.lt(address, net.const(lo, width)))
        if hi < (1 << width):
            inside = net.and_(inside, net.lt(address, net.const(hi, width)))
        for bit in range(selector_width):
            if (branch >> bit) & 1:
                outputs[bit] = net.xor(outputs[bit], inside)
    net.outputs = {"selector": outputs}
    b.call(net.operation(), address=b["target"], selector=selector)
    return StatePreparation(
        annotate(b.finish(), "state_prep_isometry", zero_input=True, clean_work=True)
    ), log_norm


@dataclass(frozen=True)
class QHAMInputModel:
    plan: object
    generator: BlockEncoding
    initial: StatePreparation
    state_width: int
    eta: complex
    log_initial_norm: float
    growth_shift: float = 0.0

    def solve(self, qode, time):
        from pyqecclang.algorithms.state_preparation import select_subspace

        state = qode(self.generator, self.initial, time)
        physical = select_subspace(state, self.state_width, 0, label="qham_physical_sum")
        return StateOracle(
            annotate(
                physical.operation,
                "unitary",
                algorithm="general_qham",
                ham_order=self.plan.order,
                pde_degree=self.plan.pde.degree,
                qham_plan=self.plan.dumps(),
                qham_log_initial_norm=self.log_initial_norm,
                qham_growth_rescale_log=self.growth_shift * time,
                output_semantics="normalized truncated HAM sum; magnitude follows chosen QODE normalization",
                correctness="HAM algebra checked; convergence and quantum solver accuracy pending",
            )
        )

    def dissipative_shift(self, shift=None):
        """给 LCHS/CBMD 显式适配；整个提升向量统一乘 exp(-shift*t)。"""
        shift = self.generator.alpha if shift is None else float(shift)
        if not math.isfinite(shift) or shift < self.generator.alpha:
            raise ValidationError("自动耗散适配要求 shift >= 已声明的 BE 范数界 alpha")
        g = lcu([(1, self.generator), (-shift, identity(self.generator.width))])
        return QHAMInputModel(
            self.plan,
            g,
            self.initial,
            self.state_width,
            self.eta,
            self.log_initial_norm,
            self.growth_shift + shift,
        )


def qham_input_model(plan, bindings, *, eta=-1.0, max_blocks=256, max_terms=4096):
    if not math.isfinite(complex(eta).real) or not math.isfinite(complex(eta).imag):
        raise ValidationError("eta 必须有限")
    g = generator_encoding(plan, bindings, eta, max_blocks=max_blocks, max_terms=max_terms)
    initial, log_norm = lifted_initial(plan, bindings)
    return QHAMInputModel(plan, g, initial, bindings.state_width, eta, log_norm)


def open_qham_input(
    plan, bindings, *, generator_alpha, generator_signal, eta=-1.0, name="QhamGenerator"
):
    """显式保留整个提升生成元未实现；不给它伪造一个空主体。"""
    bindings.validate(plan)
    width = (plan.raw_dimension(1 << bindings.state_width) - 1).bit_length()
    abstract = abstract_block_encoding(name, width, generator_signal, generator_alpha)
    attrs = dict(abstract.operation.module.attributes)
    attrs.update(
        qham_plan=plan.dumps(),
        eta_real=complex(eta).real,
        eta_imag=complex(eta).imag,
        input_model="QCL generator defined by row rules; implementation unresolved",
    )
    operation = replace(
        abstract.operation,
        module=replace(abstract.operation.module, attributes=tuple(sorted(attrs.items()))),
    )
    initial, log_norm = lifted_initial(plan, bindings)
    return QHAMInputModel(
        plan, BlockEncoding(operation), initial, bindings.state_width, eta, log_norm
    )


def taylor_qode(generator, initial, time, *, degree=2):
    """普通可闭合的有限 Taylor 候选；无 Hermitian/耗散前提，不承诺效率。"""
    from pyqecclang.algorithms.state_preparation import apply_be_to_state

    if type(degree) is not int or degree < 0 or not math.isfinite(time):
        raise ValidationError("Taylor QODE 参数无效")
    powers = [(1, identity(generator.width))]
    current = identity(generator.width)
    for order in range(1, degree + 1):
        current = product(generator, current)
        powers.append((time**order / math.factorial(order), current))
    evolution = lcu(powers)
    state = apply_be_to_state(evolution, initial)
    return StateOracle(
        annotate(
            state.operation,
            "unitary",
            algorithm="finite_taylor_qode",
            degree=degree,
            evolution_alpha=evolution.alpha,
            correctness="finite polynomial candidate; error pending",
        )
    )
