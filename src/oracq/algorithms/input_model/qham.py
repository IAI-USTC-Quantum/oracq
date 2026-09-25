"QCL row rules to modular register BEs, plus scale-preserving QODE inputs."

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, replace

from oracq.algorithms.common.arithmetic import BooleanNetwork
from oracq.algorithms.input_model.block_encoding import lcu, matrix_pauli_encoding
from oracq.algorithms.input_model.operators import (
    BlockEncoding,
    _name,
    identity,
    product,
    zero,
)
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    StatePreparation,
    abstract_block_encoding,
    abstract_state_prep,
    annotate,
    gate_state_prep,
    invoke,
    resources_for,
)
from oracq.applications.qham.linearization import Block, QHAMPlan
from oracq.applications.qham.reference import Discretization
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits, Ref, ValidationError


@dataclass(frozen=True)
class PortBinding:
    """Binding result of a QHAM elementary port: a rectangular operator encoding and its tensor arity.

    Attributes:
        encoding: Block encoding of the port operator, with target width ``max(1, arity) * state_width``.
        arity: Number of input tensor factors consumed by the port; 0 for constant injection, 1 for a linear map, r for a degree-r nonlinear port.
    """

    encoding: BlockEncoding
    arity: int


@dataclass(frozen=True)
class QHAMBindings:
    """All port and initial-value bindings of one QHAM problem.

    Attributes:
        state_width: Bit width of the physical state register, corresponding to the ``2**state_width``-dimensional discrete space.
        ports: ``(port name, PortBinding)`` tuples; the set of port names must match the PDE ports of the plan.
        initial: Preparation oracle of the physical initial state.
        initial_norm: Euclidean norm of the classical initial vector.
    """

    state_width: int
    ports: tuple[tuple[str, PortBinding], ...]
    initial: StatePreparation
    initial_norm: float

    def validate(self, plan: QHAMPlan) -> QHAMBindings:
        """Validate that the bindings match the port set and register layout of the QCL plan.

        Args:
            plan: QCL plan object providing the ``pde.ports`` port specifications.

        Returns:
            QHAMBindings: ``self`` once validation passes, for chaining.

        Raises:
            ValidationError: The state bit width or initial-value width is illegal, the norm is not a finite nonnegative number, the port name set does not match, or some port's arity or encoding width disagrees with the specification.
        """
        if not 1 <= self.state_width <= 64 or self.initial.width != self.state_width:
            raise ValidationError("QHAM state space and initial-value register layout do not match")
        if not math.isfinite(self.initial_norm) or self.initial_norm < 0:
            raise ValidationError("the initial physical vector norm must be a finite nonnegative number")
        supplied = dict(self.ports)
        if len(supplied) != len(self.ports) or set(supplied) != {p.name for p in plan.pde.ports}:
            raise ValidationError("QHAM input port set does not match")
        for p in plan.pde.ports:
            binding = supplied[p.name]
            if (
                binding.arity != p.arity
                or binding.encoding.width != max(1, p.arity) * self.state_width
            ):
                raise ValidationError("QHAM rectangular port input/output width mismatch for port " + p.name)
        return self

    @classmethod
    def declare(
        cls,
        plan: QHAMPlan,
        state_width: int,
        port_specs: Mapping[str, tuple[float, int]],
        *,
        initial_norm: float,
        initial_work: int = 0,
        prefix: str = "Qham",
    ) -> QHAMBindings:
        """Construct open port and initial-value bindings via abstract declarations.

        Generates an unimplemented block encoding declaration for each PDE port, writes the QCL plan serialization, port name, and input rank into the module attributes; the initial value likewise stays an abstract preparation, and ``validate`` runs over the whole after construction.

        Args:
            plan: QCL plan object providing the ``pde.ports`` port specifications.
            state_width: Bit width of the physical state register.
            port_specs: Mapping from port name to ``(alpha, signal)``, giving each port's BE norm bound and signal bit count.
            initial_norm: Euclidean norm of the classical initial vector.
            initial_work: Work bit width reserved by the initial-value preparation declaration.
            prefix: Prefix of the abstract declaration names; each port declaration name and the initial-value declaration name derive from this prefix.

        Returns:
            QHAMBindings: The open binding set that passed ``validate``.

        Raises:
            ValidationError: The port specification set disagrees with the plan, or the subsequent ``validate`` check fails.
        """
        if set(port_specs) != {p.name for p in plan.pde.ports}:
            raise ValidationError("the open QHAM port alpha or signal-bit specifications are incomplete")
        ports: list[tuple[str, PortBinding]] = []
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


def gate_bindings(
    discretization: Discretization,
    initial: Sequence[complex],
    *,
    max_port_qubits: int = 5,
) -> QHAMBindings:
    """Materialize only small elementary ports, without building the whole lifted matrix G.

    Args:
        discretization: Discretization object carrying port and dimension information.
        initial: Finite initial vector encoded under the full register layout.
        max_port_qubits: Maximum bit count allowed for gate implementations of elementary ports.

    Returns:
        QHAMBindings: Bindings made of small-port gate implementations and the initial state preparation.
    """
    width = discretization.width
    if len(initial) != discretization.dimension or any(
        not math.isfinite(complex(v).real) or not math.isfinite(complex(v).imag) for v in initial
    ):
        raise ValidationError("the initial value must be a finite vector encoded under the full register layout")
    ports: list[tuple[str, PortBinding]] = []
    for port in discretization.pde.ports:
        local_width = max(1, port.arity) * width
        if local_width > max_port_qubits:
            raise ValidationError("the elementary port gate implementation exceeds the small-scale materialization budget; provide an open, QRAM, or structured BE instead")
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


def _permute_slots(builder: Builder, target: Ref, width: int, desired: Sequence[int]) -> None:
    """Rearrange ``target`` into the slot order given by ``desired`` using equal-width slot swaps."""
    current = list(range(len(desired)))
    for position, token in enumerate(desired):
        other = current.index(token)
        if position != other:
            builder.swap(
                target[position * width : (position + 1) * width],
                target[other * width : (other + 1) * width],
            )
            current[position], current[other] = current[other], current[position]


def place_port(
    binding: PortBinding, state_width: int, source_rank: int, position: int
) -> BlockEncoding:
    """Place an elementary port onto the given slots of the tensor word, returning the rank-contracted block encoding.

    The port acts on the ``arity`` consecutive slots of the ``source_rank``-fold
    input tensor starting at ``position``, with identity on the remaining
    slots; for arity 0 an empty slot is first swapped into ``position`` for
    constant injection, and for arity greater than 1 the surplus consumed
    slots are permuted to the end, giving output tensor rank
    ``source_rank - arity + 1``.

    Args:
        binding: The ``PortBinding`` of the port.
        state_width: Bit width of each tensor slot.
        source_rank: Rank of the input tensor.
        position: Starting slot of the port's input window.

    Returns:
        BlockEncoding: The lifted encoding with target width ``max(source_rank, source_rank - arity + 1) * state_width``; alpha is inherited from the elementary port, with input/output rank and tensor position attributes attached.

    Raises:
        ValidationError: ``source_rank`` is smaller than the port arity, or ``position`` is outside the output rank range.
    """
    arity = binding.arity
    output_rank = source_rank - arity + 1
    rank = max(source_rank, output_rank)
    if source_rank < arity or not 0 <= position < output_rank:
        raise ValidationError("invalid rectangular tensor position")
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


def embed_rectangular(
    local: BlockEncoding,
    global_width: int,
    row_offset: int,
    column_offset: int,
    input_width: int,
    output_width: int,
) -> BlockEncoding:
    """Constrain both the input and output windows, preventing rectangular zero padding from polluting neighboring blocks.

    Args:
        local: Block encoding of the local block.
        global_width: Bit width of the global target register.
        row_offset: Row offset of the output window in the global layout.
        column_offset: Column offset of the input window in the global layout.
        input_width: Bit count of the input window.
        output_width: Bit count of the output window.

    Returns:
        BlockEncoding: The window-constrained block encoding in the global bit layout, with alpha inherited from local.
    """
    if (
        local.width > global_width
        or max(row_offset + (1 << output_width), column_offset + (1 << input_width))
        > 1 << global_width
    ):
        raise ValidationError("the rectangular block exceeds the global layout")
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

    def reject_nonzero(ref: Ref, flag: Ref) -> None:
        """Write whether ``ref`` is all zeros into ``flag``: ``flag`` is 1 when all zeros, 0 otherwise."""
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


def generator_encoding(
    plan: QHAMPlan,
    bindings: QHAMBindings,
    eta: complex,
    *,
    max_blocks: int = 256,
    max_terms: int = 4096,
) -> BlockEncoding:
    """Explicitly assemble the block encoding of the lifted generator from the row rules of the QCL plan.

    Enumerates the row couplings of each block, expressing every nonzero
    coupling as an LCU term after port placement plus rectangular embedding;
    placements with the same ``(operator, column rank, tensor position)`` are
    cached and reused.

    Args:
        plan: QCL plan object.
        bindings: Port/initial-value bindings matching the plan.
        eta: HAM homotopy parameter; terms whose coefficient evaluates to zero after substitution into the coupling weights are dropped.
        max_blocks: Maximum number of blocks enumerated explicitly.
        max_terms: Maximum number of explicit LCU coupling terms.

    Returns:
        BlockEncoding: Encoding of the whole lifted generator G, degenerating to the zero encoding with no couplings; module attributes record the plan, HAM order, explicit coupling count, and more.

    Raises:
        ValidationError: Binding validation fails, the lifted register exceeds 64 bits, or the explicit coupling count exceeds the budget.
    """
    bindings.validate(plan)
    n = bindings.state_width
    if plan.max_rank * n > 64:
        raise ValidationError("the current BE target wrapper is limited to 64 bits; a lazy QCL plan is not subject to this materialization limit")
    dimension = 1 << n
    width = (plan.raw_dimension(dimension) - 1).bit_length()
    blocks = tuple(plan.blocks(max_blocks=max_blocks))
    terms: list[tuple[complex, BlockEncoding]] = []
    placements: dict[tuple[str, int, int], BlockEncoding] = {}
    for block in blocks:
        for coupling in plan.row_terms(block):
            coefficient = coupling.weight.evaluate(eta)
            if not coefficient:
                continue
            if len(terms) >= max_terms:
                raise ValidationError("the explicit QCL BE coupling count exceeds the budget; keep the generator open instead")
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


def lifted_initial(plan: QHAMPlan, bindings: QHAMBindings) -> tuple[StatePreparation, float]:
    """Build the initial state preparation on the lifted tensor space and its logarithmic norm.

    Branch weights preserve the relative norms r, r, r**2, ..., r**K of the
    lifted blocks (with forcing, the constant component is 1): first prepare
    the selector label from the scaled weights, then within each branch invoke
    the initial-value oracle repeatedly by block rank and shift the result to
    the block offset, and finally recompute the selector label from the
    address intervals; preparations with known zero input all restore the work
    bits to zero.

    Args:
        plan: QCL plan object.
        bindings: Port/initial-value bindings.

    Returns:
        tuple: ``(StatePreparation, log_initial_norm)``, the latter being the natural logarithm of the overall norm of the lifted initial state.

    Raises:
        ValidationError: The lifted register exceeds 64 bits, or the initial value is zero with no forcing so no nonzero vector can be prepared.
    """
    n = bindings.state_width
    if plan.max_rank * n > 64:
        raise ValidationError("the current initial-value target wrapper exceeds 64 bits; keep an open initial-value model instead")
    dimension = 1 << n
    width = (plan.raw_dimension(dimension) - 1).bit_length()
    candidates = [(Block("physical"), 1.0)]
    candidates += [(Block("tensor", (0,) * k), float(k)) for k in range(1, plan.max_rank + 1)]
    if plan.has_forcing:
        candidates.append((Block("tensor"), 0.0))
    if bindings.initial_norm == 0:
        candidates = [item for item in candidates if item[1] == 0]
        if not candidates:
            raise ValidationError("the lifted initial value is zero; without forcing the classical side should return the zero solution")
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
    # Each preparation with known zero input restores work to zero, so the same work register can be reused sequentially.
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
    intervals: list[tuple[int, int, int]] = []
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
    # The blocks have disjoint support; recompute the selector label from the address intervals.
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
    """The QODE input model assembled by QHAM.

    The system being solved is the lifted linear system ``Y' = G Y``; the
    inhomogeneous forcing has been homogenized through the constant component,
    and G is not the QLSS matrix to be inverted.

    Attributes:
        plan: QCL plan object that generated this model.
        generator: Block encoding of the lifted generator G.
        initial: Preparation of the lifted initial state.
        state_width: Bit width of the physical state register.
        eta: HAM homotopy parameter.
        log_initial_norm: Natural logarithm of the lifted initial state norm.
        growth_shift: Accumulated explicit dissipation shift, growing linearly in time as the logarithmic factor for magnitude recovery.
    """

    plan: object
    generator: BlockEncoding
    initial: StatePreparation
    state_width: int
    eta: complex
    log_initial_norm: float
    growth_shift: float = 0.0

    def solve(
        self,
        qode: Callable[[BlockEncoding, StatePreparation, float], StateOracle],
        time: float,
    ) -> StateOracle:
        """Evolve the lifted system with the given QODE protocol and select the physical output block.

        Args:
            qode: Three-argument linear solver protocol ``(generator, initial, time) -> StateOracle``.
            time: Evolution duration.

        Returns:
            StateOracle: The normalized truncated HAM sum over the first physical block; the physical magnitude still requires the norm recovery information of the chosen protocol, and the module attributes record the plan, logarithmic initial norm, and shift factor.
        """
        from oracq.algorithms.common.state_preparation import select_subspace

        state = qode(self.generator, self.initial, time)
        physical = select_subspace(state, self.state_width, 0, label="qham_physical_sum")
        return StateOracle(
            annotate(
                physical.operation,
                "unitary",
                algorithm="general_qham",
                ham_order=self.plan.order,  # type: ignore[attr-defined]
                pde_degree=self.plan.pde.degree,  # type: ignore[attr-defined]
                qham_plan=self.plan.dumps(),  # type: ignore[attr-defined]
                qham_log_initial_norm=self.log_initial_norm,
                qham_growth_rescale_log=self.growth_shift * time,
                output_semantics="normalized truncated HAM sum; magnitude follows chosen QODE normalization",
                correctness="HAM algebra checked; convergence and quantum solver accuracy pending",
            )
        )

    def dissipative_shift(self, shift: float | None = None) -> QHAMInputModel:
        """Explicit adaptation for LCHS/CBMD; the whole lifted vector is uniformly multiplied by exp(-shift*t).

        Args:
            shift: Explicit dissipation shift; defaults to the generator's alpha and must not be smaller than it.

        Returns:
            QHAMInputModel: New input model whose generator becomes G − shift·I.
        """
        shift = self.generator.alpha if shift is None else float(shift)
        if not math.isfinite(shift) or shift < self.generator.alpha:
            raise ValidationError("automatic dissipation adaptation requires shift >= the declared BE norm bound alpha")
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


def qham_input_model(
    plan: QHAMPlan,
    bindings: QHAMBindings,
    *,
    eta: complex = -1.0,
    max_blocks: int = 256,
    max_terms: int = 4096,
) -> QHAMInputModel:
    """Assemble an explicit QHAM input model: the lifted generator block encoding plus the lifted initial state.

    Args:
        plan: QCL plan object.
        bindings: Port/initial-value bindings matching the plan.
        eta: HAM homotopy parameter, which must be a finite complex number.
        max_blocks: Maximum number of blocks enumerated explicitly.
        max_terms: Maximum number of explicit LCU coupling terms.

    Returns:
        QHAMInputModel: The input model ready to be solved by a QODE protocol.

    Raises:
        ValidationError: eta is non-finite, or a check fails while assembling the generator or initial state.
    """
    if not math.isfinite(complex(eta).real) or not math.isfinite(complex(eta).imag):
        raise ValidationError("eta must be finite")
    g = generator_encoding(plan, bindings, eta, max_blocks=max_blocks, max_terms=max_terms)
    initial, log_norm = lifted_initial(plan, bindings)
    return QHAMInputModel(plan, g, initial, bindings.state_width, eta, log_norm)


def open_qham_input(
    plan: QHAMPlan,
    bindings: QHAMBindings,
    *,
    generator_alpha: float,
    generator_signal: int,
    eta: complex = -1.0,
    name: str = "QhamGenerator",
) -> QHAMInputModel:
    """Explicitly keep the whole lifted generator unimplemented; do not fake an empty body for it.

    Args:
        plan: QCL plan object.
        bindings: Port/initial-value bindings matching the plan.
        generator_alpha: Norm bound of the open generator declaration.
        generator_signal: Signal bit width of the open generator declaration.
        eta: HAM homotopy parameter, which must be a finite complex number.
        name: Name of the open generator declaration, ``QhamGenerator`` by default.

    Returns:
        QHAMInputModel: Input model whose generator stays an open declaration awaiting later binding.
    """
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


def taylor_qode(
    generator: BlockEncoding,
    initial: StatePreparation,
    time: float,
    *,
    degree: int = 2,
) -> StateOracle:
    """A plain closable finite Taylor candidate; no Hermitian or dissipative assumption, and no efficiency is promised.

    Args:
        generator: Block encoding of the generator.
        initial: Initial state preparation.
        time: Evolution duration, which must be a finite real number.
        degree: Taylor truncation order, a nonnegative integer.

    Returns:
        StateOracle: Output state oracle after evolution by the truncated Taylor polynomial.
    """
    from oracq.algorithms.common.state_preparation import apply_be_to_state

    if type(degree) is not int or degree < 0 or not math.isfinite(time):
        raise ValidationError("invalid Taylor QODE parameters")
    powers: list[tuple[float, BlockEncoding]] = [(1, identity(generator.width))]
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
