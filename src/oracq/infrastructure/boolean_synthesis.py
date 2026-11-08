"""Exact affine synthesis of Boolean networks with clean nonlinear workspace.

XOR and NOT nodes are represented by linear forms instead of automatically
owning private qubits. Shared nontrivial affine values are materialized once
to bound recomputation. At most two scratch bits materialize AND controls,
and are cleaned immediately.
The complete compute/copy/uncompute circuit preserves arbitrary output words.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from collections.abc import Iterator, Mapping
from typing import TYPE_CHECKING, cast

from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Adjoint, Bits, Ref, fuse

if TYPE_CHECKING:
    from oracq.algorithms.common.arithmetic import BooleanNetwork


def _terms(form: int) -> Iterator[int]:
    """Yield nonconstant node indices of a GF(2) form in ascending order."""
    form &= ~1
    while form:
        bit = form & -form
        yield bit.bit_length() - 1
        form ^= bit


def lower_affine(
    network: BooleanNetwork,
    name: str | None = None,
    *,
    attributes: Mapping[str, str | int | float | bool] | None = None,
) -> Operation:
    """Synthesize the exact Boolean function, without allocating XOR/NOT nodes.

    A form's bit zero denotes its constant term. Other set bits name input
    bits or previously computed values. ANDs with zero, one, equal or
    complementary operands are folded. Backward reachability is checked again
    after folding, so dead nonlinear work is not emitted either.
    """
    needed: set[int] = set()
    pending = [node for bits in network.outputs.values() for node in bits]
    while pending:
        index = pending.pop()
        if index < 2 or index in needed:
            continue
        needed.add(index)
        record = network.nodes[index]
        if record[0] != "input":
            pending.extend(cast("tuple[int, ...]", record[1:]))

    uses: Counter[int] = Counter(node for bits in network.outputs.values() for node in bits)
    for index in needed:
        if network.nodes[index][0] != "input":
            uses.update(cast("tuple[int, ...]", network.nodes[index][1:]))
    forms = {0: 0, 1: 1}
    plans: dict[int, tuple[int, ...]] = {}
    for index in sorted(needed):
        op, *args = network.nodes[index]
        if op == "input":
            forms[index] = 1 << index
            continue
        left = forms[cast("int", args[0])]
        if op == "not":
            forms[index] = left ^ 1
        elif op == "xor":
            forms[index] = left ^ forms[cast("int", args[1])]
            if uses[index] > 1 and (forms[index] >> 1).bit_count() > 1:
                plans[index] = (forms[index],)
                forms[index] = 1 << index
        else:
            right = forms[cast("int", args[1])]
            if left == 0 or right == 0 or left == right ^ 1:
                forms[index] = 0
            elif left == 1 or left == right:
                forms[index] = right
            elif right == 1:
                forms[index] = left
            else:
                plans[index] = (left, right)
                forms[index] = 1 << index

    live: set[int] = set()
    pending = [forms[node] for bits in network.outputs.values() for node in bits]
    while pending:
        for index in _terms(pending.pop()):
            if index in plans and index not in live:
                live.add(index)
                pending.extend(plans[index])
    consumers: Counter[int] = Counter()
    control_inputs: set[int] = set()
    for index in live:
        for form in plans[index]:
            consumers.update(_terms(form))
            if len(plans[index]) == 2:
                control_inputs.update(_terms(form))
    for bits in network.outputs.values():
        for node in bits:
            consumers.update(_terms(forms[node]))
    # An AND used only by one affine copy can XOR directly into its target.
    # It must remain stored if another AND needs it as a control.
    inline = {
        index for index in live
        if len(plans[index]) == 2 and consumers[index] == 1 and index not in control_inputs
    }
    ordered = sorted(live - inline)
    scratch_count = max(
        (sum((form >> 1).bit_count() > 1 for form in plans[i]) for i in live if len(plans[i]) == 2),
        default=0,
    )
    payload = network.payload()
    symbol = name or "arith_" + hashlib.sha256((payload + ":affine-v1").encode()).hexdigest()[:20]
    b = Builder(
        symbol,
        {key: Bits(len(bits)) for key, bits in {**network.inputs, **network.outputs}.items()},
        attributes={
            "arithmetic_network": payload,
            "arithmetic_lowering": "affine-v1",
            "correctness": "pending",
            "workspace_contract": "zero_in_zero_out",
            **(attributes or {}),
        },
    )
    refs = {node: b[key][i] for key, bits in network.inputs.items() for i, node in enumerate(bits)}
    banks = [
        b.local("value_" + str(i // 64), Bits(min(64, len(ordered) - i)))
        for i in range(0, len(ordered), 64)
    ]
    refs.update({node: banks[i // 64][i % 64] for i, node in enumerate(ordered)})
    scratch = b.local("affine_scratch", Bits(scratch_count)) if scratch_count else None

    def copy(form: int, target: Ref) -> None:
        if form & 1:
            b.x(target)
        for index in _terms(form):
            if index in inline:
                product(plans[index], target)
            else:
                b.xor(refs[index], target)

    def product(factors: tuple[int, ...], target: Ref) -> None:
        controls: list[Ref] = []
        values = 0
        temporaries: list[tuple[int, Ref]] = []
        for form in factors:
            terms = list(_terms(form))
            if len(terms) == 1:
                control, value = refs[terms[0]], 1 ^ (form & 1)
            else:
                assert scratch is not None
                control, value = scratch[len(temporaries)], 1
                copy(form, control)
                temporaries.append((form, control))
            values |= value << len(controls)
            controls.append(control)
        with b.control(fuse(*controls), values):
            b.x(target)
        for form, control in reversed(temporaries):
            copy(form, control)

    for index in ordered:
        if len(plans[index]) == 1:
            copy(plans[index][0], refs[index])
        else:
            product(plans[index], refs[index])

    forward = tuple(b._frames[0])
    for key, bits in network.outputs.items():
        for bit, node in enumerate(bits):
            copy(forms[node], b[key][bit])
    if forward:
        b.emit(Adjoint(forward))
    return b.finish()
