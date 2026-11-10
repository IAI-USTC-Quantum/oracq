"""Optional read-only QRAM assisted Boolean and phase-gradient compilation.

QRAM is an explicit external resource, never a free gate hidden in the ledger.
The returned extra memory contains a single four-word AND truth table. A shared
phase gradient is prepared and unprepared in the entry circuit, so its synthesis
and workspace are counted. Calls and Repeat nodes remain compressed.
"""

from __future__ import annotations

import math
from dataclasses import replace

from oracq.infrastructure.builder import Builder
from oracq.infrastructure.ir import (
    QRAM,
    Adjoint,
    Bits,
    Call,
    Control,
    Instruction,
    Load,
    Module,
    Primitive,
    Program,
    Ref,
    Register,
    Repeat,
    Resource,
    Span,
    ValidationError,
    fuse,
)
from oracq.infrastructure.linking import uses_store
from oracq.infrastructure.validation import validate

_PREFIX = "__qr_"
_BANK = _PREFIX + "truth"
_GRAD = _PREFIX + "gradient"
_CTRL = _PREFIX + "controls"


def _view(name: str, width: int) -> Ref:
    return Ref((Span(name, 0, width),), Bits(width))


def _bits(ref: Ref) -> tuple[Ref, ...]:
    return tuple(ref[i] for i in range(ref.width))


def _overlaps(left: Ref, right: Ref) -> bool:
    return any(
        a.register == b.register and a.start < b.start + b.width and b.start < a.start + a.width
        for a in left.parts
        for b in right.parts
    )


def _group(nodes: tuple[Instruction, ...], start: int):
    first = nodes[start]
    if not (
        isinstance(first, Control)
        and first.register.width == 1
        and first.value == 1
        and len(first.body) == 1
        and isinstance(first.body[0], Primitive)
    ):
        return None
    rotation = first.body[0]
    if rotation.op != "ry" or rotation.operands[0].width != 1 or not rotation.angle:
        return None
    ratio = 2 * math.pi / abs(rotation.angle)
    width = round(math.log2(ratio)) if math.isfinite(ratio) and ratio > 0 else 0
    if not 2 <= width <= 32 or abs(rotation.angle) != 2 * math.pi / (1 << width):
        return None
    if start + width > len(nodes):
        return None
    refs = []
    for bit, node in enumerate(nodes[start : start + width]):
        if not (
            isinstance(node, Control)
            and node.register.width == 1
            and node.value == 1
            and len(node.body) == 1
            and isinstance(node.body[0], Primitive)
        ):
            return None
        gate = node.body[0]
        if (
            gate.op != "ry"
            or gate.operands != rotation.operands
            or gate.angle != rotation.angle * (1 << bit)
        ):
            return None
        refs.append(node.register)
    angle = fuse(*refs)
    if len(set(refs)) != width or _overlaps(angle, rotation.operands[0]):
        return None
    return width, angle, rotation.operands[0], 1 if rotation.angle > 0 else -1


def _groups(nodes: tuple[Instruction, ...]):
    i = 0
    while i < len(nodes):
        group = _group(nodes, i)
        if group:
            yield group
            i += group[0]
        else:
            node = nodes[i]
            if isinstance(node, (Control, Repeat, Adjoint)):
                yield from _groups(node.body)
            i += 1


def _adder(width: int) -> Module:
    """Cuccaro modular addition, preserving the left word and carry input."""
    b = Builder(_PREFIX + "adder", {"a": Bits(width), "b": Bits(width), "carry": Bits(1)})

    def majority(a, target, carry):
        b.xor(a, target)
        b.xor(a, carry)
        with b.control(carry), b.control(target):
            b.x(a)

    def unmajority(a, target, carry):
        with b.control(carry), b.control(target):
            b.x(a)
        b.xor(a, carry)
        b.xor(carry, target)

    for bit in range(width - 1):
        majority(b["a"][bit], b["b"][bit], b["carry"] if bit == 0 else b["a"][bit - 1])
    b.xor(b["a"][width - 1], b["b"][width - 1])
    b.xor(b["carry"] if width == 1 else b["a"][width - 2], b["b"][width - 1])
    for bit in reversed(range(width - 1)):
        unmajority(b["a"][bit], b["b"][bit], b["carry"] if bit == 0 else b["a"][bit - 1])
    return b.finish().module


def compile_qram(program: Program, *, phase_gradient: bool = True, boolean: bool = True):
    """Return ``(compiled_program, extra_memory, report)`` with the same public interface.

    The gradient is reusable throughout one coherent entry execution. Arbitrary
    non-dyadic rotations are left unchanged. The four-word truth table implements
    XOR, including nonzero targets, superposed addresses and negative controls.
    Write-bearing or open programs are rejected, not partially optimized.
    """
    validate(program, require_closed=True)
    if uses_store(program):
        raise ValidationError("QRAM assisted compilation requires read-only memory")
    if dict(program.main.attributes).get("qram_assisted_compilation"):
        return program, {}, {"applied": False, "reason": "already_compiled"}
    for module in program.modules:
        names = ([module.name] + [r.name for r in module.registers]
                 + [r.name for r in module.locals] + [r.name for r in module.resources])
        if any(name.startswith(_PREFIX) for name in names):
            raise ValidationError("reserved QRAM compiler namespace already in use")
    source = program.module_map.copy()

    def calls(nodes):
        for node in nodes:
            if isinstance(node, Call):
                yield node.module
            elif isinstance(node, (Control, Repeat, Adjoint)):
                yield from calls(node.body)

    reachable = {program.entry}
    pending = [program.entry]
    while pending:
        for callee in calls(source[pending.pop()].body or ()):
            if callee not in reachable:
                reachable.add(callee)
                pending.append(callee)
    found = (
        [g for name in sorted(reachable) for g in _groups(source[name].body or ())]
        if phase_gradient
        else []
    )
    gradient_width = max((g[0] + 1 for g in found), default=0)
    if gradient_width:
        source[_PREFIX + "adder"] = _adder(gradient_width)
    gradient_users = {
        name for name in reachable if phase_gradient and any(_groups(source[name].body or ()))
    }
    while True:
        expanded = gradient_users | {
            name
            for name, m in source.items()
            if any(c in gradient_users for c in calls(m.body or ()))
        }
        if expanded == gradient_users:
            break
        gradient_users = expanded
    generated: dict[tuple[str, tuple[int, ...]], Module] = {}
    metrics = {
        "phase_gradient_bits": gradient_width,
        "rotation_word_groups": 0,
        "boolean_xor_replacements": 0,
        "controlled_call_alias_fallbacks": 0,
    }

    def compile_module(name: str, values: tuple[int, ...]) -> Module:
        key = name, values
        if key in generated:
            return generated[key]
        original = source[name]
        new_name = _PREFIX + str(len(generated)) + "_" + name
        # Reserve a stable name before visiting callees (the source graph is a DAG).
        generated[key] = replace(original, name=new_name)
        controls = tuple((_view(_CTRL, len(values))[i], value) for i, value in enumerate(values))
        gradient = _view(_GRAD, gradient_width)
        scratch_width = 0
        padding_width = 0
        need_carry = False

        def normalize(active):
            result = dict()
            for ref, value in active:
                if ref in result and result[ref] != value:
                    return None
                result[ref] = value
            return tuple(result.items())

        def controlled(nodes, active):
            nonlocal scratch_width
            if not active:
                return tuple(nodes)
            if boolean and len(active) >= 2 and not all(isinstance(n, Load) for n in nodes):
                refs_ = [r for r, _ in active]
                negative = tuple(Primitive("x", (r,)) for r, v in active if not v)
                scratch_width = max(scratch_width, len(refs_) - 1)
                flags = _view(_PREFIX + "flags", len(refs_) - 1)
                query = tuple(
                    Load(_BANK, fuse(refs_[0] if i == 0 else flags[i - 1], refs_[i + 1]), flags[i])
                    for i in range(len(refs_) - 1)
                )
                return (
                    *negative,
                    *query,
                    Control(flags[len(refs_) - 2], 1, tuple(nodes)),
                    *reversed(query),
                    *reversed(negative),
                )
            refs = fuse(*(r for r, _ in active))
            mask = sum(v << k for k, (_, v) in enumerate(active))
            return (Control(refs, mask, tuple(nodes)),)

        def flip(target, active):
            nonlocal scratch_width
            active = normalize(active)
            if active is None:
                return ()
            if not boolean or len(active) < 2:
                return controlled((Primitive("x", (target,)),), active)
            metrics["boolean_xor_replacements"] += 1
            negative = tuple(Primitive("x", (r,)) for r, v in active if not v)
            refs = [r for r, _ in active]
            scratch_width = max(scratch_width, len(refs) - 2)
            flags = _view(_PREFIX + "flags", max(0, len(refs) - 2))
            prefix = []
            for i in range(len(refs) - 2):
                prefix.append(
                    Load(_BANK, fuse(refs[0] if i == 0 else flags[i - 1], refs[i + 1]), flags[i])
                )
            final = Load(
                _BANK, fuse(refs[0] if len(refs) == 2 else flags[len(refs) - 3], refs[-1]), target
            )
            return (*negative, *prefix, final, *reversed(prefix), *reversed(negative))

        def call(node, active):
            active = normalize(active)
            if active is None:
                return ()
            # A legal controlled call may read its control through an existing
            # argument. Keep that original alias contract instead of adding an
            # overlapping formal argument forbidden by RIR call validation.
            alias = any(_overlaps(r, a) for r, _ in active for a in node.arguments)
            passing = () if alias else active
            callee = compile_module(node.module, tuple(v for _, v in passing))
            args = node.arguments + ((gradient,) if node.module in gradient_users else ())
            if passing:
                args += (fuse(*(r for r, _ in passing)),)
            resources = node.resources + ((_BANK,) if boolean else ())
            new = Call(callee.name, args, resources)
            if alias:
                metrics["controlled_call_alias_fallbacks"] += 1
                return controlled((new,), active)
            return (new,)

        def word_rotation(group, active):
            nonlocal padding_width, need_carry
            width, angle, target, sign = group
            padding_width = max(padding_width, gradient_width - width)
            padding = _view(_PREFIX + "padding", gradient_width - width)
            shift = gradient_width - width - 1
            word = fuse(padding[:shift], angle, padding[shift:])
            carry = _view(_PREFIX + "carry", 1)
            need_carry = True
            add = Call(_PREFIX + "adder", (word, gradient, carry))
            plus = call(add, (*active, (target, 0)))
            minus = call(add, (*active, (target, 1)))
            if sign < 0:
                plus, minus = minus, plus
            metrics["rotation_word_groups"] += 1
            return (
                Primitive("phase", (target,), angle=-math.pi / 2),
                Primitive("h", (target,)),
                *plus,
                Adjoint(minus),
                Primitive("h", (target,)),
                Primitive("phase", (target,), angle=math.pi / 2),
            )

        def body(nodes, active):
            result = []
            i = 0
            while i < len(nodes):
                group = _group(nodes, i) if gradient_width else None
                if group:
                    result.extend(word_rotation(group, active))
                    i += group[0]
                    continue
                node = nodes[i]
                i += 1
                if isinstance(node, Control):
                    extra = tuple(
                        (r, (node.value >> j) & 1) for j, r in enumerate(_bits(node.register))
                    )
                    resolved = normalize((*active, *extra))
                    if resolved is not None:
                        result.extend(body(node.body, resolved))
                elif isinstance(node, Call):
                    result.extend(call(node, active))
                elif isinstance(node, (Repeat, Adjoint)):
                    result.append(replace(node, body=tuple(body(node.body, active))))
                elif (
                    isinstance(node, Primitive)
                    and boolean
                    and node.op in {"x", "xor", "swap", "add_const"}
                ):
                    if node.op == "x":
                        for target in _bits(node.operands[0]):
                            result.extend(flip(target, active))
                    elif node.op in {"xor", "swap"}:
                        for left, right in zip(*(_bits(r) for r in node.operands), strict=True):
                            result.extend(flip(right, (*active, (left, 1))))
                            if node.op == "swap":
                                result.extend(flip(left, (*active, (right, 1))))
                                result.extend(flip(right, (*active, (left, 1))))
                    else:
                        target = node.operands[0]
                        value = node.value % (1 << target.width)
                        for offset in range(target.width):
                            if value >> offset & 1:
                                for high in reversed(range(offset + 1, target.width)):
                                    result.extend(
                                        flip(
                                            target[high],
                                            (
                                                *active,
                                                *((target[j], 1) for j in range(offset, high)),
                                            ),
                                        )
                                    )
                                result.extend(flip(target[offset], active))
                else:
                    result.extend(controlled((node,), active))
            return result

        compiled = tuple(body(original.body or (), controls))
        registers = original.registers + (
            (Register(_GRAD, Bits(gradient_width)),) if name in gradient_users else ()
        )
        if values:
            registers += (Register(_CTRL, Bits(len(values))),)
        resources = original.resources + ((Resource(_BANK, QRAM(2, 1)),) if boolean else ())
        locals_ = original.locals
        if scratch_width:
            locals_ += (Register(_PREFIX + "flags", Bits(scratch_width)),)
        if padding_width:
            locals_ += (Register(_PREFIX + "padding", Bits(padding_width)),)
        if need_carry:
            locals_ += (Register(_PREFIX + "carry", Bits(1)),)
        result = replace(
            original,
            name=new_name,
            registers=registers,
            resources=resources,
            locals=locals_,
            body=compiled,
        )
        generated[key] = result
        return result

    entry = compile_module(program.entry, ())
    wrapper = Builder(
        _PREFIX + "entry",
        {r.name: r.type for r in program.main.registers},
        {r.name: r.type for r in entry.resources},
        attributes={
            **dict(program.main.attributes),
            "qram_assisted_compilation": True,
            "gradient_bits": gradient_width,
        },
    )
    gradient = (
        wrapper.local(_GRAD, Bits(gradient_width)) if program.entry in gradient_users else None
    )
    preparation: list[Instruction] = []
    if gradient is not None:
        for bit in range(gradient_width):
            preparation.extend(
                (
                    Primitive("h", (gradient[bit],)),
                    Primitive(
                        "phase",
                        (gradient[bit],),
                        angle=2 * math.pi * (1 << bit) / (1 << gradient_width),
                    ),
                )
            )
    args = tuple(wrapper[r.name] for r in program.main.registers) + (
        (gradient,) if gradient is not None else ()
    )
    shell_body = (*preparation, Call(entry.name, args, tuple(r.name for r in entry.resources)))
    if preparation:
        shell_body += (Adjoint(tuple(preparation)),)
    shell = Module(
        wrapper.name,
        wrapper.registers,
        wrapper.resources,
        shell_body,
        wrapper.attributes,
        tuple(wrapper.locals),
    )
    compiled = replace(
        program,
        entry=shell.name,
        modules=tuple(sorted((shell, *generated.values()), key=lambda m: m.name)),
    )
    validate(compiled)
    report = {
        "applied": True,
        **metrics,
        "source_modules": len(program.modules),
        "compiled_modules": len(compiled.modules),
        "truth_table_words": 4 if boolean else 0,
        "truth_table_nonzero_words": 1 if boolean else 0,
        "gradient_preparation_and_cleanup_in_circuit": bool(gradient_width),
    }
    return compiled, {_BANK: {3: 1}} if boolean else {}, report
