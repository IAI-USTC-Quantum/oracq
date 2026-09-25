"""Binding experiments pairing one diagonal-operator algorithm with a gate table, QRAM and reversible integer arithmetic."""

from __future__ import annotations

import math

from oracq.algorithms.input_model.oracles import (
    abstract_database,
    annotate,
    diagonal_block_encoding,
    gate_database,
    qram_database,
)
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import Bits, Program, UInt
from oracq.infrastructure.linking import Binding


def oracle_study(
    width: int = 2, repetitions: int = 1,
) -> tuple[Program, dict[str, tuple[Binding, dict[str, list[int]]]]]:
    """Build an open program and three implementations; the data word is ``(address + 1) mod 2**width``.

    The program prepares target uniformly and repeatedly calls the angle-word
    driven diagonal rotation. The zero-signal block diagonal element is
    ``cos(repetitions * pi * word / 2**width)``. This identity is used only for
    the rotation dilation in this example and is not claimed as a general
    property of powers of arbitrary block encodings.
    """
    slot = abstract_database("AngleWord", width, width)
    encoding = diagonal_block_encoding(slot)
    b = Builder(
        f"OracleStudy_{width}_{repetitions}",
        {"target": Bits(width), "signal": Bits(width + 1)},
    )
    b.h(b["target"])
    with b.repeat(repetitions):
        b.call(encoding.operation, target=b["target"], signal=b["signal"])

    arithmetic = Builder(f"IncrementWord_{width}", {"address": Bits(width), "data": Bits(width)})
    work = arithmetic.local("word", UInt(width))
    arithmetic.xor(arithmetic["address"], work)
    arithmetic.add_const(work, 1)
    arithmetic.xor(work, arithmetic["data"])
    with arithmetic.adjoint():
        arithmetic.xor(arithmetic["address"], work)
        arithmetic.add_const(work, 1)
    table = [(i + 1) % (1 << width) for i in range(1 << width)]
    variants = {
        "gate_table": (Binding(gate_database(width, width, table).operation), {}),
        "qram_table": (
            Binding(qram_database(width, width).operation, {"table": "angle_words"}),
            {"angle_words": table},
        ),
        "arithmetic": (Binding(annotate(arithmetic.finish(), "database_xor")), {}),
    }
    return b.finish().program(), variants


def oracle_study_reference(width: int, repetitions: int) -> dict[tuple[int, ...], complex]:
    """Independent trigonometric reference, including the complex amplitudes of the success and failure signals."""
    size = 1 << width
    state: dict[tuple[int, ...], complex] = {}
    for address in range(size):
        angle = repetitions * math.pi * ((address + 1) % size) / size
        state[address, 0] = complex(math.cos(angle) / math.sqrt(size))
        state[address, 1 << width] = complex(math.sin(angle) / math.sqrt(size))
    return state
