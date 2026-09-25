"""Jordan quantum gradient estimation: a phase oracle input model and single-query gradient readout.

Implements the gradient estimation of Jordan 2005 (PRL 95, 050501). Grid
convention: the i-th coordinate register occupies the bit segment
[i*grid_bits, (i+1)*grid_bits), and the integer value k corresponds to the
fixed-point grid point x = k/N with N = 2**grid_bits. The scaling convention
matches the generalization of Gilyén–Arunachalam–Wiebe 2019: the phase
oracle implements ``O|x> = exp(2πi·N·f(x))|x>``, i.e. phase_scale must equal
N. When f is approximately linear on the grid, a single oracle call plus an
inverse QFT per coordinate writes N·∂f/∂x_i into the i-th coordinate
register; one query yields all d components, while classical deterministic
evaluation of the same gradient takes O(d) function queries.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from typing import cast

from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.algorithms.common.fourier import inverse_qft
from oracq.algorithms.input_model.contracts import (
    OracleView,
    fail,
    finite_real,
    positive_integer,
    validate_signature,
)
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import annotate, declare, invoke, resources_for
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError
from oracq.infrastructure.mathfunc import MathConfig


@dataclass(frozen=True)
class PhaseOracle(OracleView):
    """Grid phase oracle view: the diagonal phase action of ``O|x> = exp(2πi·phase_scale·f(x))|x>``.

    Integer values of target encode d-dimensional fixed-point grid points;
    the scaling factor is recorded in the phase_scale attribute so the
    algorithm side can verify the scaling convention at generation time."""

    oracle_kind = "phase_oracle"
    operation: Operation

    def phase_oracle(self) -> PhaseOracle:
        """Return the phase oracle view; this class itself wraps a phase oracle, so it returns itself.

        Returns:
            PhaseOracle: A reference to itself, keeping the role accessor
            interface uniform.
        """
        return self

    def __post_init__(self) -> None:
        """Validate the target signature and the phase_scale attribute declaration."""
        validate_signature(self.operation, ("target",), "PhaseOracle")
        scale = dict(self.operation.module.attributes).get("phase_scale")
        if scale is None:
            raise ValidationError("a phase oracle must declare the phase_scale attribute")
        finite_real(scale, "PhaseOracle.phase_scale", minimum=0, strict=True)

    @property
    def width(self) -> int:
        """Bit width of the target register, i.e. the total width of the grid register the phase oracle acts on."""
        return next(r.type.width for r in self.operation.module.registers if r.name == "target")

    @property
    def phase_scale(self) -> float:
        """The phase scaling factor declared by the oracle, i.e. the scale in ``O|x> = exp(2πi·phase_scale·f(x))|x>``."""
        return cast("float", dict(self.operation.module.attributes)["phase_scale"])


def abstract_phase_oracle(name: str, width: int, *, phase_scale: float) -> PhaseOracle:
    """Open-declare a phase oracle; the phase scale is written into the phase_scale attribute and the implementation is left for bind.

    Args:
        name: Declared module name of the phase oracle slot.
        width: Bit width of the target register, in 1..64.
        phase_scale: Phase scaling factor, a positive real, registered along
            with the attribute declaration.

    Returns:
        PhaseOracle: Handle of the phase oracle slot with an empty body,
        whose implementation is bound later by bind.
    """
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


def gate_phase_oracle(
    width: int,
    angles: Iterable[float],
    *,
    phase_scale: float,
    name: str | None = None,
) -> PhaseOracle:
    """Build a diagonal phase oracle from an explicit phase table; angles[x] is the phase in radians gained by basis state ``|x>``.

    Small-scale witnesses can enumerate all 2**width phases directly; larger
    grids should use function_phase_oracle instead or bind another
    implementation.

    Args:
        width: Bit width of the target register, in 1..64.
        angles: Table of phases in radians per basis state, of length exactly
            2**width.
        phase_scale: Phase scaling factor, a positive real, registered into
            the attributes.
        name: Name of the generated module; generated from the parameters by
            default.

    Returns:
        PhaseOracle: Oracle handle that applies the diagonal phase per basis
        state according to the phase table.
    """
    positive_integer(width, "gate_phase_oracle.width", maximum=64)
    finite_real(phase_scale, "gate_phase_oracle.phase_scale", minimum=0, strict=True)
    angles = tuple(angles)
    if len(angles) != 1 << width:
        raise ValidationError("the phase table length must equal 2 to the power of width")
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
    source: str | Callable[..., float],
    *,
    dimension: int,
    grid_bits: int,
    fmt: FixedFormat | None = None,
    scale: float | None = None,
    name: str | None = None,
    constants: Mapping[str, bool | int | float | complex] | None = None,
    helpers: Mapping[str, Callable[..., object]] | None = None,
    config: MathConfig | None = None,
    max_unroll: int = 128,
    entry: str | None = None,
) -> PhaseOracle:
    """Build a phase oracle from mathfunc arithmetic: compute fixed-point f(x), kick the phase back from the output, then restore.

    Args:
        source: Source of a pure Python function with signature
            f(x0, ..., x{d-1}) returning a single real.
        dimension: Grid dimension d, in 1..16.
        grid_bits: Bits m per coordinate; the grid points x = k/2**m must be
            exactly representable by fmt.
        fmt: Fixed-point format, default FixedFormat(grid_bits+12,
            grid_bits+8); must be signed with fraction >= grid_bits.
        scale: Phase scale, default 2**grid_bits, i.e. the Jordan scaling
            convention.
        name: Overrides the automatically generated module name.
        constants: Constant bindings passed to the mathfunc frontend.
        helpers: Helper functions passed to the mathfunc frontend.
        config: MathConfig passed to the mathfunc frontend.
        max_unroll: Unrolling cap passed to the mathfunc frontend.
        entry: Entry function name passed to the mathfunc frontend.

    Returns:
        PhaseOracle: Phase ``exp(2πi·scale·f(x))``; when f leaves the value
        range of fmt this is flagged by status and the value is not
        guaranteed.

    The coordinate registers write the function inputs through zero-state
    work bits, the phase kickback decodes the output in two's complement,
    and a subsequent inverse invocation restores all work bits."""
    positive_integer(dimension, "function_phase_oracle.dimension", maximum=16)
    positive_integer(grid_bits, "function_phase_oracle.grid_bits", maximum=32)
    if dimension * grid_bits > 64:
        raise ValidationError("the total grid width must not exceed 64 bits")
    fmt = fmt or FixedFormat(grid_bits + 12, grid_bits + 8)
    if not isinstance(fmt, FixedFormat):
        raise ValidationError("fmt must be a FixedFormat")
    if not fmt.signed or fmt.fraction < grid_bits:
        raise ValidationError("the fixed-point format needs a sign bit and a fraction of at least grid_bits to represent grid coordinates exactly")
    scale = (1 << grid_bits) if scale is None else scale
    finite_real(scale, "function_phase_oracle.scale", minimum=0, strict=True)
    from oracq.infrastructure.mathfunc import compile_function

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
    operation = cast("Operation", compiled.operation)
    expected = {f"x{i}" for i in range(dimension)} | {"out", "status"}
    if {r.name for r in operation.module.registers} != expected:
        raise ValidationError("the function must take dimension real arguments and return a single real")
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
            math_function=cast("str", dict(operation.module.attributes).get("math_function")),
        )
    )


def gradient_estimation(
    oracle: PhaseOracle | Operation, *, dimension: int, grid_bits: int
) -> Operation:
    """Generate the Jordan gradient estimation circuit.

    Args:
        oracle: A PhaseOracle or a phase oracle operation with the target
            signature implementing ``O|x> = exp(2πi·N·f(x))|x>``.
        dimension: Grid dimension d, in 1..16.
        grid_bits: Bits m per coordinate register, N = 2**m, in 1..32.

    Returns:
        Operation: Register target of width d*m. Decode the gradient
        components with gradient_from_readout after readout.

    Raises:
        ValidationError: Invalid grid parameters, a mismatched oracle width,
        or a phase_scale not equal to 2**grid_bits.

    Each coordinate register is prepared in uniform superposition, the phase
    oracle is invoked once, then an inverse QFT is applied per coordinate.
    When f is approximately linear the phase writes N·∂f/∂x_i into the
    Fourier basis of coordinate register i, and after the inverse QFT the
    readout is the fixed-point approximation of each component. One query
    yields all d components, while classical deterministic gradient
    evaluation takes O(d) function queries (Jordan 2005, PRL 95, 050501)."""
    positive_integer(dimension, "gradient.dimension", maximum=16)
    positive_integer(grid_bits, "gradient.grid_bits", maximum=32)
    width = dimension * grid_bits
    if width > 64:
        raise ValidationError("the total grid width must not exceed 64 bits")
    view = oracle if isinstance(oracle, PhaseOracle) else PhaseOracle(oracle)
    if view.width != width:
        raise ValidationError("the phase oracle width must equal dimension times grid_bits")
    grid_points = 1 << grid_bits
    if view.phase_scale != grid_points:
        fail(
            "INPUT_PROMISE",
            "gradient.oracle.phase_scale",
            grid_points,
            view.phase_scale,
            "the scaling convention requires phase_scale to equal 2 to the power of grid_bits",
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


def gradient_from_readout(value: int, *, dimension: int, grid_bits: int) -> tuple[float, ...]:
    """Decode the integer readout of target into per-component gradient estimates.

    Args:
        value: Integer readout of the target register.
        dimension: Grid dimension.
        grid_bits: Bits m per coordinate register.

    Returns:
        tuple: The i-th component interprets the bit segment [i*m, (i+1)*m)
        in two's complement and divides by 2**m."""
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
