"按需遍历模块图，并提供基于寄存器整数元组的小规模参考执行器。"

from __future__ import annotations

import cmath
import math
from collections.abc import Mapping
from dataclasses import dataclass

from pyqecclang.infrastructure.ir import (
    Adjoint,
    Call,
    Control,
    Load,
    Primitive,
    Ref,
    Repeat,
    ValidationError,
)
from pyqecclang.infrastructure.validation import locations, validate


def check_memory(program, memory):
    memory = {} if memory is None else memory
    if not isinstance(memory, Mapping):
        raise ValidationError("QRAM 数据必须按资源名提供映射")
    specs = {r.name: r.type for r in program.main.resources}
    if set(memory) != set(specs):
        raise ValidationError("必须为入口的每个 QRAM 提供数据，且不能有多余资源")
    result = {}
    for name, spec in specs.items():
        raw = memory[name]
        try:
            items = raw.items() if isinstance(raw, Mapping) else enumerate(raw)
        except TypeError as exc:
            raise ValidationError("每个 QRAM 必须提供字序列或稀疏字典") from exc
        cells = {}
        for address, value in items:
            if type(address) is not int or not 0 <= address < 1 << spec.address_width:
                raise ValidationError("QRAM 地址越界")
            if type(value) is not int or not 0 <= value < 1 << spec.data_width:
                raise ValidationError("QRAM 字必须是位宽内的无符号整数")
            if value:
                cells[address] = value
        result[name] = cells
    return result


def _counter(program, limit, native_modules=frozenset()):
    modules, cached = program.module_map, {}

    def count(nodes):
        total = 0
        for node in nodes:
            if isinstance(node, Primitive):
                cost = max(4, sum(r.width for r in node.operands) ** 2)
            elif isinstance(node, Load):
                cost = 1
            elif isinstance(node, Call):
                if node.module in native_modules:
                    cost = 1
                    total = min(limit + 1, total + cost)
                    continue
                if node.module not in cached:
                    cached[node.module] = count(modules[node.module].body)
                cost = 1 + cached[node.module]
            elif isinstance(node, Repeat):
                cost = node.count * (1 + count(node.body)) if node.body else 0
            elif isinstance(node, Control):
                cost = 1 + (1 + 2 * node.register.width) * count(node.body)
            else:
                cost = 1 + count(node.body)
            total = min(limit + 1, total + cost)
        return total

    return count


def expanded_steps(program, limit=1_000_000, native_modules=frozenset()):
    if program.entry in native_modules:
        return 1
    return _counter(program, limit, native_modules)(program.main.body)


@dataclass(frozen=True)
class LocalEnter:
    name: str
    type: object


@dataclass(frozen=True)
class LocalExit:
    name: str
    width: int


def remap(ref, mapping):
    parts = []
    for span in ref.parts:
        parts.extend(mapping[span.register][span.start : span.start + span.width].parts)
    return Ref(tuple(parts), ref.type)


def events(program, *, max_steps=1_000_000, native_modules=frozenset()):
    """调用逐层展开为迭代器，原始 RIR 保持不变。"""
    from pyqecclang.infrastructure.ir import Span

    program = validate(program, require_closed=not bool(native_modules))
    counter = _counter(program, max_steps, native_modules)
    if expanded_steps(program, max_steps, native_modules) > max_steps:
        raise ValidationError("执行超过展开预算")
    modules = program.module_map

    def walk(nodes, mapping, resources, controls=(), inverse=False):
        for node in reversed(nodes) if inverse else nodes:
            if isinstance(node, Call):
                target = modules[node.module]
                actuals = {
                    r.name: remap(a, mapping)
                    for r, a in zip(target.registers, node.arguments, strict=True)
                }
                bound = {
                    r.name: resources[a]
                    for r, a in zip(target.resources, node.resources, strict=True)
                }
                yield from enter_module(target, actuals, bound, controls, inverse)
            elif isinstance(node, Repeat):
                if counter(node.body):
                    for _ in range(node.count):
                        yield from walk(node.body, mapping, resources, controls, inverse)
            elif isinstance(node, Control):
                yield from walk(
                    node.body,
                    mapping,
                    resources,
                    controls + ((remap(node.register, mapping), node.value),),
                    inverse,
                )
            elif isinstance(node, Adjoint):
                yield from walk(node.body, mapping, resources, controls, not inverse)
            elif isinstance(node, Primitive):
                yield (
                    Primitive(
                        node.op,
                        tuple(remap(r, mapping) for r in node.operands),
                        node.angle,
                        node.value,
                    ),
                    controls,
                    inverse,
                )
            elif isinstance(node, Load):
                yield (
                    Load(
                        resources[node.resource],
                        remap(node.address, mapping),
                        remap(node.data, mapping),
                    ),
                    controls,
                    inverse,
                )

    local_counter = 0

    def enter_module(module, mapping, resources, controls=(), inverse=False):
        nonlocal local_counter
        if module.name in native_modules:
            from pyqecclang.infrastructure.native import NativeSite

            yield (
                NativeSite(
                    module,
                    tuple(mapping[r.name] for r in module.registers),
                    tuple(resources[r.name] for r in module.resources),
                ),
                controls,
                inverse,
            )
            return
        if module.body is None:
            raise ValidationError(f"未实现模块：{module.name}")
        mapping = dict(mapping)
        allocated = []
        for register in module.locals:
            local_counter += 1
            name = f"__local_{local_counter}"
            mapping[register.name] = Ref((Span(name, 0, register.type.width),), register.type)
            if register.type.width:
                allocated.append((name, register.type.width))
                yield LocalEnter(name, register.type), (), False
        yield from walk(module.body, mapping, resources, controls, inverse)
        for name, width in reversed(allocated):
            yield LocalExit(name, width), (), False

    roots = {r.name: Ref((Span(r.name, 0, r.type.width),), r.type) for r in program.main.registers}
    resources = {r.name: r.name for r in program.main.resources}
    yield from enter_module(program.main, roots, resources)


def gate_matrix(op, angle=None, inverse=False):
    if op == "h":
        a = 1 / math.sqrt(2)
        matrix = ((a, a), (a, -a))
    elif op == "x":
        matrix = ((0, 1), (1, 0))
    elif op == "y":
        matrix = ((0, -1j), (1j, 0))
    elif op == "z":
        matrix = ((1, 0), (0, -1))
    elif op in {"s", "t", "phase"}:
        theta = {"s": math.pi / 2, "t": math.pi / 4}.get(op, angle)
        matrix = ((1, 0), (0, cmath.exp(1j * theta)))
    elif op == "rz":
        matrix = ((cmath.exp(-0.5j * angle), 0), (0, cmath.exp(0.5j * angle)))
    elif op in {"rx", "ry"}:
        c, s = math.cos(angle / 2), math.sin(angle / 2)
        matrix = ((c, -1j * s), (-1j * s, c)) if op == "rx" else ((c, -s), (s, c))
    else:
        raise ValidationError(f"未知单比特矩阵：{op}")
    if inverse:
        return tuple(tuple(complex(matrix[j][i]).conjugate() for j in range(2)) for i in range(2))
    return matrix


@dataclass(frozen=True)
class RegisterState:
    registers: tuple
    amplitudes: dict[tuple[int, ...], complex]

    def statevector(self, max_qubits=20):
        width = sum(r.type.width for r in self.registers)
        if width > max_qubits:
            raise ValidationError("密集状态转换超过位数预算")
        vector = [0j] * (1 << width)
        for values, amplitude in self.amplitudes.items():
            index, offset = 0, 0
            for reg, value in zip(self.registers, values, strict=True):
                index |= value << offset
                offset += reg.type.width
            vector[index] = amplitude
        return vector


def simulate(program, memory=None, *, initial=None, max_steps=1_000_000, max_states=65536):
    program = validate(program, require_closed=True)
    memories = check_memory(program, memory)
    registers = program.main.registers
    index = {r.name: i for i, r in enumerate(registers)}
    initial = initial or {}
    if set(initial) - set(index):
        raise ValidationError("初态包含未知寄存器")
    values = []
    for reg in registers:
        value = initial.get(reg.name, 0)
        if type(value) is not int or not 0 <= value < (1 << reg.type.width):
            raise ValidationError("初态寄存器值越界")
        values.append(value)
    state = {tuple(values): 1 + 0j}

    def read(ref, values):
        result, offset = 0, 0
        for span in ref.parts:
            result |= (
                (values[index[span.register]] >> span.start) & ((1 << span.width) - 1)
            ) << offset
            offset += span.width
        return result

    def write(ref, values, value):
        result, offset = list(values), 0
        for span in ref.parts:
            pos = index[span.register]
            mask = ((1 << span.width) - 1) << span.start
            result[pos] = (result[pos] & ~mask) | (
                ((value >> offset) & ((1 << span.width) - 1)) << span.start
            )
            offset += span.width
        return tuple(result)

    for node, controls, inverse in events(program, max_steps=max_steps):
        if isinstance(node, LocalEnter):
            index[node.name] = len(next(iter(state)))
            state = {values + (0,): amplitude for values, amplitude in state.items()}
            continue
        if isinstance(node, LocalExit):
            position = index.pop(node.name)
            if any(values[position] != 0 for values in state):
                raise ValidationError(f"局部寄存器未复净：{node.name}")
            state = {
                values[:position] + values[position + 1 :]: amplitude
                for values, amplitude in state.items()
            }
            continue

        def active(values, controls=controls):
            return all(read(ref, values) == expected for ref, expected in controls)

        if isinstance(node, Primitive) and node.op not in {"xor", "swap", "add_const", "gphase"}:
            matrix = gate_matrix(node.op, node.angle, inverse)
            for register, bit in locations(node.operands[0]):
                result = {}
                pos = index[register]
                for values, amplitude in state.items():
                    if not active(values):
                        result[values] = result.get(values, 0j) + amplitude
                        continue
                    incoming = (values[pos] >> bit) & 1
                    for outgoing in range(2):
                        target = list(values)
                        target[pos] = (target[pos] & ~(1 << bit)) | (outgoing << bit)
                        target = tuple(target)
                        result[target] = (
                            result.get(target, 0j) + matrix[outgoing][incoming] * amplitude
                        )
                state = {k: v for k, v in result.items() if abs(v) > 1e-15}
                if len(state) > max_states:
                    raise ValidationError("参考模拟超过稀疏态数量预算")
            continue
        result = {}
        for values, amplitude in state.items():
            target = values
            if active(values):
                if isinstance(node, Load):
                    value = memories[node.resource].get(read(node.address, values), 0)
                    target = write(node.data, values, read(node.data, values) ^ value)
                elif node.op == "gphase":
                    amplitude *= cmath.exp(1j * node.angle * (-1 if inverse else 1))
                elif node.op == "xor":
                    source, dest = node.operands
                    target = write(dest, values, read(source, values) ^ read(dest, values))
                elif node.op == "swap":
                    a, b = node.operands
                    av, bv = read(a, values), read(b, values)
                    target = write(b, write(a, values, bv), av)
                elif node.op == "add_const":
                    ref = node.operands[0]
                    value = read(ref, values) + node.value * (-1 if inverse else 1)
                    target = write(ref, values, value & ((1 << ref.width) - 1))
            result[target] = result.get(target, 0j) + amplitude
        state = {k: v for k, v in result.items() if abs(v) > 1e-15}
    return RegisterState(registers, state)
