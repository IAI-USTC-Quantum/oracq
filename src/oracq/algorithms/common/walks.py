"""Discrete-time coined quantum walk."""

from __future__ import annotations

from oracq.algorithms.input_model.contracts import positive_integer
from oracq.algorithms.input_model.operators import _name
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits


def cycle_walk(width: int, *, steps: int = 1) -> Operation:
    """Generate a Hadamard coined walk on a periodic lattice.

    Args:
        width: Position bit width; the period length is 2**width.
        steps: Nonnegative number of steps.

    Returns:
        Operation: Exposes position and a one-bit coin. Each step updates the coin
        first, then moves by +1 or -1 according to the 0/1 value.

    The input state is prepared by the caller; the zero input corresponds to starting
    at position zero with coin zero."""
    positive_integer(width, "cycle_walk.width", maximum=64)
    positive_integer(steps, "cycle_walk.steps", minimum=0)
    b = Builder(
        _name("cycle_walk", width, steps),
        {"position": Bits(width), "coin": Bits(1)},
        attributes={"algorithm": "coined_cycle_walk", "steps": steps},
    )
    with b.repeat(steps):
        b.h(b["coin"])
        with b.control(b["coin"], 0):
            b.add_const(b["position"].reinterpret("uint"), 1)
        with b.control(b["coin"], 1):
            b.add_const(b["position"].reinterpret("uint"), (1 << width) - 1)
    return b.finish()
