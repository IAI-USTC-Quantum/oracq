"""Coherent encoding and single-error recovery circuits of the three-bit repetition code."""

from __future__ import annotations

from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError, fuse


def repetition_encode(*, error: str = "bit") -> Operation:
    """Encode one logical bit into the three-bit repetition code.

    Args:
        error: ``bit`` for a single X error, ``phase`` for a single Z error.

    Returns:
        Operation: target is the one-bit logical input; syndrome is two bits and must be
        zero on input.

    The three encoded physical bits are laid out as target, syndrome[0], syndrome[1]."""
    if error not in {"bit", "phase"}:
        raise ValidationError("Repetition code error must be bit or phase")
    b = Builder(
        "repetition_encode_" + error,
        {"target": Bits(1), "syndrome": Bits(2)},
        attributes={"algorithm": "repetition_encode", "error_kind": error},
    )
    b.xor(b["target"], b["syndrome"][0])
    b.xor(b["target"], b["syndrome"][1])
    if error == "phase":
        b.h(fuse(b["target"], b["syndrome"]))
    return b.finish()


def repetition_recover(*, error: str = "bit") -> Operation:
    """Coherently recover a single error of the given kind in the three-bit repetition code.

    Args:
        error: ``bit`` or ``phase``, matching the encoder.

    Returns:
        Operation: target restores the logical state; the error information stays in syndrome.

    No measurement or reset is included; a nonzero syndrome must not be treated as
    already-clean work space."""
    if error not in {"bit", "phase"}:
        raise ValidationError("Repetition code error must be bit or phase")
    b = Builder(
        "repetition_recover_" + error,
        {"target": Bits(1), "syndrome": Bits(2)},
        attributes={
            "algorithm": "repetition_recover",
            "error_kind": error,
            "syndrome_policy": "retained; host reset required before reuse",
        },
    )
    if error == "phase":
        b.h(fuse(b["target"], b["syndrome"]))
    b.xor(b["target"], b["syndrome"][1])
    b.xor(b["target"], b["syndrome"][0])
    with b.control(b["syndrome"], 3):
        b.x(b["target"])
    return b.finish()
