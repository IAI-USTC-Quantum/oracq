"QCL 行规则到模块化寄存器 BE，以及保留尺度的 QODE 输入。"

from __future__ import annotations

import math
from dataclasses import dataclass, replace

from pyqecclang.algorithms.common.arithmetic import BooleanNetwork
from pyqecclang.algorithms.input_model.block_encoding import lcu, matrix_pauli_encoding
from pyqecclang.algorithms.input_model.operators import (
    BlockEncoding,
    _name,
    identity,
    product,
    zero,
)
from pyqecclang.algorithms.input_model.oracles import (
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
    """QHAM 基础端口的绑定结果：一个矩形算子编码及其张量元数。

    Attributes:
        encoding: 端口算子的块编码，target 宽度为 ``max(1, arity) * state_width``。
        arity: 端口消耗的输入张量因子数；常量注入为 0，线性映射为 1，r 次非线性端口为 r。
    """

    encoding: BlockEncoding
    arity: int


@dataclass(frozen=True)
class QHAMBindings:
    """一个 QHAM 问题的全部端口与初值绑定。

    Attributes:
        state_width: 物理状态寄存器位宽，对应 ``2**state_width`` 维离散空间。
        ports: ``(端口名, PortBinding)`` 元组，端口名集合须与 plan 的 PDE 端口一致。
        initial: 物理初态的制备 oracle。
        initial_norm: 经典初值向量的 Euclidean 范数。
    """

    state_width: int
    ports: tuple[tuple[str, PortBinding], ...]
    initial: StatePreparation
    initial_norm: float

    def validate(self, plan):
        """校验绑定与 QCL plan 的端口集合和寄存器布局是否匹配。

        Args:
            plan: QCL 计划对象，提供 ``pde.ports`` 端口规格。

        Returns:
            QHAMBindings: 校验通过的 ``self``，便于链式书写。

        Raises:
            ValidationError: 状态位宽或初值宽度非法、范数不是有限非负数、端口名集合不匹配，或某端口的元数/编码宽度与规格不符。
        """
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
        """以抽象声明方式构造开放的端口与初值绑定。

        为每个 PDE 端口生成一个未实现的块编码声明，并把 QCL plan 序列、端口名与输入秩写入模块属性；初值同样保留为抽象制备，构造完成后整体执行 ``validate``。

        Args:
            plan: QCL 计划对象，提供 ``pde.ports`` 端口规格。
            state_width: 物理状态寄存器位宽。
            port_specs: 从端口名到 ``(alpha, signal)`` 的映射，给出各端口的 BE 范数界与信号位数。
            initial_norm: 经典初值向量的 Euclidean 范数。
            initial_work: 初值制备声明预留的工作位宽度。
            prefix: 抽象声明名的前缀；各端口声明名与初值声明名由该前缀派生。

        Returns:
            QHAMBindings: 通过 ``validate`` 校验的开放绑定集合。

        Raises:
            ValidationError: 端口规格集合与 plan 不一致，或随后的 ``validate`` 校验失败。
        """
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
    """把基础端口放置到张量字的指定槽位，返回秩收缩后的块编码。

    端口作用在 ``source_rank`` 重输入张量从 ``position`` 开始的连续 ``arity`` 个槽位上，其余槽位作恒等映射；元数为 0 时先把一个空槽换入 ``position`` 供常量注入，元数大于 1 时把多余的被消耗槽位置换到末尾，输出张量秩为 ``source_rank - arity + 1``。

    Args:
        binding: 端口的 ``PortBinding``。
        state_width: 每个张量槽位的位宽。
        source_rank: 输入张量的秩。
        position: 端口输入窗口的起始槽位。

    Returns:
        BlockEncoding: 目标宽度为 ``max(source_rank, source_rank - arity + 1) * state_width`` 的提升编码，alpha 继承自基础端口，并附带输入/输出秩与张量位置属性。

    Raises:
        ValidationError: ``source_rank`` 小于端口元数，或 ``position`` 不在输出秩范围内。
    """
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
    """按 QCL plan 的行规则显式装配提升生成元的块编码。

    枚举各块的行耦合，把每个非零耦合表示为端口放置加矩形嵌入后的 LCU 项；相同 ``(算子, 列秩, 张量位置)`` 的放置被缓存复用。

    Args:
        plan: QCL 计划对象。
        bindings: 与 plan 匹配的端口/初值绑定。
        eta: HAM 同伦参数，代入各耦合权重后系数为零的项被剔除。
        max_blocks: 显式枚举的块数上限。
        max_terms: 显式 LCU 耦合项数上限。

    Returns:
        BlockEncoding: 整个提升生成元 G 的编码，无耦合时退化为零编码；模块属性记录 plan、HAM 阶数与显式耦合数等。

    Raises:
        ValidationError: 绑定校验失败、提升寄存器超过 64 位，或显式耦合数超过预算。
    """
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
    """构造提升张量空间的初态制备及其对数范数。

    分支权重保持各提升块的相对范数 r, r, r**2, ..., r**K（有强迫时常量分量为 1）：先按缩放后的权重制备选择标签，各分支内按块秩重复调用初值 oracle 并把结果平移到块偏移，最后用地址区间反算选择标签；对已知零输入的制备均复净工作位。

    Args:
        plan: QCL 计划对象。
        bindings: 端口/初值绑定。

    Returns:
        tuple: ``(StatePreparation, log_initial_norm)``，后者为提升初态整体范数的自然对数。

    Raises:
        ValidationError: 提升寄存器超过 64 位，或初值为零且无强迫因而不存在可制备的非零向量。
    """
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
    """QHAM 装配完成的 QODE 输入模型。

    求解对象是提升后的线性系统 ``Y' = G Y``，非齐次强迫已通过常量分量齐次化；G 不是需要求逆的 QLSS 矩阵。

    Attributes:
        plan: 生成该模型的 QCL 计划对象。
        generator: 提升生成元 G 的块编码。
        initial: 提升初态的制备。
        state_width: 物理状态寄存器位宽。
        eta: HAM 同伦参数。
        log_initial_norm: 提升初态范数的自然对数。
        growth_shift: 已累积的显式耗散移位量，作为幅值恢复的对数因子随时间线性增长。
    """

    plan: object
    generator: BlockEncoding
    initial: StatePreparation
    state_width: int
    eta: complex
    log_initial_norm: float
    growth_shift: float = 0.0

    def solve(self, qode, time):
        """用给定的 QODE 协议演化提升系统并选出物理输出块。

        Args:
            qode: 三参数线性求解协议 ``(generator, initial, time) -> StateOracle``。
            time: 演化时长。

        Returns:
            StateOracle: 第一个物理块的归一化截断 HAM 和；物理幅值仍需结合所选协议的范数恢复信息，模块属性记录了 plan、对数初值范数与移位因子。
        """
        from pyqecclang.algorithms.common.state_preparation import select_subspace

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
    """装配显式 QHAM 输入模型：提升生成元块编码加提升初态。

    Args:
        plan: QCL 计划对象。
        bindings: 与 plan 匹配的端口/初值绑定。
        eta: HAM 同伦参数，须为有限复数。
        max_blocks: 显式枚举的块数上限。
        max_terms: 显式 LCU 耦合项数上限。

    Returns:
        QHAMInputModel: 可直接交给 QODE 协议求解的输入模型。

    Raises:
        ValidationError: eta 非有限，或生成元/初态装配过程中的校验失败。
    """
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
    from pyqecclang.algorithms.common.state_preparation import apply_be_to_state

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
