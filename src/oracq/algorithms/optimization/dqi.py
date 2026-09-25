"""GF(2) max-XORSAT generator for DQI (Decoded Quantum Interferometry).

Implements the circuit skeleton of Jordan et al. 2024 (arXiv:2408.08292,
figure 4): prepare a Dicke state of weight l on the m-bit error register,
apply the right-hand-side phase ``(-1)**(v.y)``, reversibly compute ``B^T y``
into the n-bit syndrome register, then unload the error register back to
``|0>`` with a reversible classical decoder, and finally apply the Hadamard
transform to the syndrome register and measure. The decoder is part of the
input model and plugs in through the three-layer paradigm of open
declaration (abstract) plus batched binding (bind).

Only GF(2) is supported for now: the GF(q) case needs a q-ary discrete
Fourier transform and generalized Dicke states, and the register would
additionally have to be organized into subregisters of ``log2(q)`` bits,
left as future work. Keeping the branches where the error register is zero
after measurement (postselection) is up to the caller.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from math import comb

from oracq.algorithms.input_model.contracts import (
    OracleView,
    fail,
    positive_integer,
    require_instance,
    validate_signature,
)
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
    annotate,
    declare,
    gate_state_prep,
    invoke,
    resources_for,
)
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError

# Explicit-amplitude and bruteforce decoding only serve small-instance witnesses; beyond that scale, switch to dedicated preparation circuits or efficient decoders.
_EXPLICIT_LIMIT = 16


@dataclass(frozen=True)
class XorSatInstance:
    """GF(2) max-XORSAT instance: a sparse row representation of the constraint system ``Bx = v``.

    rows holds one constraint per row, listing the variable indices that take
    part in that constraint, with no repetition inside a row; rhs is the
    right-hand side ``v``, of the same length as the number of rows and valued
    0 or 1; num_variables is the number of variables n, also the bit width of
    the syndrome register.
    """

    rows: tuple[tuple[int, ...], ...]
    rhs: tuple[int, ...]
    num_variables: int

    def __post_init__(self) -> None:
        """Normalize and validate the construction-time constraints on constraint rows, the right-hand side, and variable indices."""
        positive_integer(self.num_variables, "XorSatInstance.num_variables", maximum=64)
        rows = tuple(tuple(row) for row in self.rows)
        rhs = tuple(self.rhs)
        if not rows or len(rows) != len(rhs):
            fail(
                "CONFIG_VALUE",
                "XorSatInstance.rows",
                "non-empty and of the same length as rhs",
                (len(rows), len(rhs)),
                "constraint rows must be non-empty and of the same length as the right-hand side",
            )
        for row in rows:
            if not row or len(set(row)) != len(row):
                fail(
                    "CONFIG_VALUE",
                    "XorSatInstance.rows",
                    "non-empty with no duplicate indices",
                    row,
                    "constraint rows must be non-empty with no duplicated variable indices inside a row",
                )
            for index in row:
                positive_integer(
                    index,
                    "XorSatInstance.index",
                    minimum=0,
                    maximum=self.num_variables - 1,
                )
        if any(type(bit) is not int or bit not in (0, 1) for bit in rhs):
            fail(
                "CONFIG_VALUE",
                "XorSatInstance.rhs",
                "0 or 1",
                rhs,
                "right-hand side entries must be 0 or 1 and must not be bool",
            )
        object.__setattr__(self, "rows", rows)
        object.__setattr__(self, "rhs", rhs)

    @property
    def num_constraints(self) -> int:
        """Number of constraints m, i.e. the bit width of the error register."""
        return len(self.rows)

    def satisfied_count(self, assignment: int) -> int:
        """Count the constraints satisfied by a given integer assignment.

        Args:
            assignment: Integer encoding of the assignment, in 0..2^number of
                variables−1, with bit j holding the value of variable j.

        Returns:
            int: Number of constraints satisfied by that assignment.
        """
        positive_integer(
            assignment,
            "XorSatInstance.assignment",
            minimum=0,
            maximum=(1 << self.num_variables) - 1,
        )
        return sum(
            (sum((assignment >> j) & 1 for j in row) & 1) == bit
            for row, bit in zip(self.rows, self.rhs, strict=True)
        )


@dataclass(frozen=True)
class DecoderOracle(OracleView):
    """Reversible form of a syndrome decoder, the key open input of DQI.

    The semantics are ``|syndrome, error> -> |syndrome, error XOR
    D(syndrome)>``, where ``D`` is the classical decoding function (for
    example a reversible implementation of belief propagation). Undecodable
    syndromes may map to arbitrary error patterns; the corresponding branches
    are eliminated by postselection.
    """

    oracle_kind = "reversible_function"
    operation: Operation

    def decoder(self) -> DecoderOracle:
        """Decoder role accessor, returns itself; consistent with the role methods of other ``OracleView`` classes.

        Returns:
            DecoderOracle: A reference to itself, keeping the role accessor
            interface uniform.
        """
        return self

    def __post_init__(self) -> None:
        """Validate that the wrapped operation has the syndrome and error signature."""
        validate_signature(self.operation, ("syndrome", "error"), "DecoderOracle")

    @property
    def syndrome_width(self) -> int:
        """Bit width of the syndrome register, read directly from the RIR register signature."""
        return next(r.type.width for r in self.operation.module.registers if r.name == "syndrome")

    @property
    def error_width(self) -> int:
        """Bit width of the error register, read directly from the RIR register signature."""
        return next(r.type.width for r in self.operation.module.registers if r.name == "error")


def abstract_decoder(name: str, syndrome_width: int, error_width: int) -> DecoderOracle:
    """Declare a decoder slot; efficient classical decoding algorithms are bound in batches via bind.

    Args:
        name: Declared module name of the decoder slot.
        syndrome_width: Bit width of the syndrome register, in 1..64.
        error_width: Bit width of the error register, in 1..64.

    Returns:
        DecoderOracle: Handle of the decoder slot with an empty body, whose
        implementation is bound later by bind.
    """
    positive_integer(syndrome_width, "abstract_decoder.syndrome_width", maximum=64)
    positive_integer(error_width, "abstract_decoder.error_width", maximum=64)
    return DecoderOracle(
        declare(
            name,
            {"syndrome": Bits(syndrome_width), "error": Bits(error_width)},
            paradigm="reversible_function",
            attributes={"decoder_role": "syndrome_decode"},
        )
    )


def table_decoder(
    syndrome_width: int, error_width: int, table: Mapping[int, int], *, name: str | None = None
) -> DecoderOracle:
    """A gate witness implementing the decoder with an explicit lookup table; only the syndromes listed in the table are decoded.

    Args:
        syndrome_width: Bit width of the syndrome register, i.e. the
            instance's variable count n.
        error_width: Bit width of the error register, i.e. the instance's
            constraint count m.
        table: Map from syndrome integer values to error pattern integer
            values.
        name: Optional module name.

    Returns:
        DecoderOracle: Syndromes not listed map to the zero error pattern."""
    positive_integer(syndrome_width, "table_decoder.syndrome_width", maximum=64)
    positive_integer(error_width, "table_decoder.error_width", maximum=64)
    items = tuple(sorted(table.items()))
    b = Builder(
        name or _name("decoder_table", syndrome_width, error_width, items),
        {"syndrome": Bits(syndrome_width), "error": Bits(error_width)},
    )
    for syndrome, error in items:
        if (
            type(syndrome) is not int
            or type(error) is not int
            or not (0 <= syndrome < 1 << syndrome_width and 0 <= error < 1 << error_width)
        ):
            raise ValidationError("a syndrome or error pattern in the decoder lookup table is out of range")
        if error:
            with b.control(b["syndrome"], syndrome):
                for bit in range(error_width):
                    if (error >> bit) & 1:
                        b.x(b["error"][bit])
    return DecoderOracle(
        annotate(b.finish(), "reversible_function", implementation="gate_truth_table")
    )


def bruteforce_decoder(
    instance: XorSatInstance, *, max_weight: int | None = None, name: str | None = None
) -> DecoderOracle:
    """Bruteforce minimum-weight decoder: gives each syndrome the lightest error of weight at most max_weight.

    Enumerates all 2**m error patterns and only accepts small instances with
    m <= 16; larger instances should bind a reversible implementation of an
    efficient decoder. On equal-weight ties the error pattern that comes
    first in enumeration order is kept.

    Args:
        instance: The XOR instance being decoded, with at most 16 constraints.
        max_weight: Weight cap for error patterns, in 0..constraint count;
            defaults to the constraint count.
        name: Optional name of the generated module.

    Returns:
        DecoderOracle: An explicit lookup-table decoder mapping each syndrome
        to the lightest error pattern.
    """
    require_instance(instance, XorSatInstance, "bruteforce_decoder.instance")
    m = instance.num_constraints
    if m > _EXPLICIT_LIMIT:
        fail(
            "CONFIG_VALUE",
            "bruteforce_decoder.m",
            f"<= {_EXPLICIT_LIMIT}",
            m,
            "the bruteforce decoder only accepts small instances",
        )
    if max_weight is None:
        max_weight = m
    positive_integer(max_weight, "bruteforce_decoder.max_weight", minimum=0, maximum=m)
    row_masks = [sum(1 << j for j in row) for row in instance.rows]
    best: dict[int, int] = {}
    for error in range(1 << m):
        weight = error.bit_count()
        if weight > max_weight:
            continue
        syndrome = 0
        for i, mask in enumerate(row_masks):
            if (error >> i) & 1:
                syndrome ^= mask
        if syndrome not in best or weight < best[syndrome].bit_count():
            best[syndrome] = error
    return table_decoder(
        instance.num_variables,
        m,
        best,
        name=name
        or _name(
            "decoder_bruteforce",
            instance.rows,
            instance.rhs,
            instance.num_variables,
            max_weight,
        ),
    )


def dicke_state(m: int, weight: int) -> StatePreparation:
    """Prepare the m-bit Dicke state ``|D_l^m>`` of weight weight.

    Built by multi-way rotations from explicit amplitudes and only accepts
    small sizes with m <= 16; large-scale preparation needs dedicated Dicke
    state circuits (such as the O(l*m) construction of Bartschi-Eidenbenz),
    which can be attached separately as a state_prep_isometry open
    declaration.

    Args:
        m: Number of target bits, in 1..16.
        weight: Hamming weight, in 0..m.

    Returns:
        StatePreparation: A zero-input preparation whose target width is m."""
    positive_integer(m, "dicke_state.m", maximum=_EXPLICIT_LIMIT)
    positive_integer(weight, "dicke_state.weight", minimum=0, maximum=m)
    amplitude = 1 / math.sqrt(comb(m, weight))
    amplitudes = [amplitude if value.bit_count() == weight else 0 for value in range(1 << m)]
    return gate_state_prep(amplitudes, name=f"dicke_{m}_{weight}")


def dqi(instance: XorSatInstance, decoder: DecoderOracle, *, weight: int) -> Operation:
    """Assemble the main DQI circuit (GF(2), the single-weight version of figure 4 in the paper).

    Args:
        instance: XorSatInstance; the constraint count m and variable count n
            determine the register widths.
        decoder: DecoderOracle whose interface widths must match the
            instance's n and m; it may be an abstract declaration.
        weight: Dicke state weight l, in 0..m; meaningful only when the
            decoding radius covers this weight.

    Returns:
        Operation: Registers error(m) and syndrome(n). Starting from the
        all-zero state, error returns to ``|0>`` on the branches where
        decoding succeeds, and after the Hadamard transform syndrome measures
        assignment x with probability proportional to ``K_l(u(x))**2``, where
        u(x) is the number of unsatisfied constraints and ``K_l`` is the
        Krawtchouk polynomial. Postselecting on zero error is up to the
        caller."""
    require_instance(instance, XorSatInstance, "dqi.instance")
    require_instance(decoder, DecoderOracle, "dqi.decoder")
    m, n = instance.num_constraints, instance.num_variables
    positive_integer(weight, "dqi.weight", minimum=0, maximum=m)
    if decoder.syndrome_width != n or decoder.error_width != m:
        fail(
            "INPUT_WIDTH",
            "dqi.decoder",
            (n, m),
            (decoder.syndrome_width, decoder.error_width),
            "the decoder widths must match the n and m of the instance",
        )
    preparation = dicke_state(m, weight)
    b = Builder(
        _name("dqi", instance, decoder.operation, weight),
        {"error": Bits(m), "syndrome": Bits(n)},
        resources_for(("prep", preparation.operation), ("decoder", decoder.operation)),
        attributes={
            "algorithm": "dqi",
            "field": "GF(2)",
            "num_constraints": m,
            "num_variables": n,
            "dicke_weight": weight,
            "readout_register": "syndrome",
            "postselection": "error_zero",
        },
    )
    invoke(b, preparation.operation, "prep", target=b["error"], work=b["syndrome"][:0])
    for i, bit in enumerate(instance.rhs):
        if bit:
            b.z(b["error"][i])
    for i, row in enumerate(instance.rows):
        for j in row:
            b.xor(b["error"][i], b["syndrome"][j])
    invoke(b, decoder.operation, "decoder", syndrome=b["syndrome"], error=b["error"])
    b.h(b["syndrome"])
    return b.finish()
