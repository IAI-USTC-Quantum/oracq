"按需遍历模块图，并提供基于寄存器整数元组的小规模参考执行器。"

from __future__ import annotations

import cmath
import math
from collections.abc import Callable, Collection, Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, cast

from pyqecclang.infrastructure.ir import (
    Adjoint,
    Call,
    Control,
    Instruction,
    Load,
    Module,
    Primitive,
    Program,
    Ref,
    RegType,
    Repeat,
    Span,
    Store,
    ValidationError,
)
from pyqecclang.infrastructure.validation import locations, validate

if TYPE_CHECKING:
    from pyqecclang.infrastructure.native import NativeSite


def check_memory(
    program: Program,
    memory: Mapping[str, Sequence[int] | Mapping[int, int]] | None,
) -> dict[str, dict[int, int]]:
    """校验并规整入口的 QRAM 绑定数据。

    Args:
        program: 封闭 RIR 程序。
        memory: 按资源名提供的数据；每项为字序列或 ``地址 -> 字`` 的稀疏字典，``None`` 视为全零。

    Returns:
        dict: 资源名到稀疏单元映射的字典；值为零的单元被丢弃，缺省单元按零解释。

    Raises:
        ValidationError: 数据不是按资源名的映射、资源集合与入口声明不符、地址越界或字不是位宽内的无符号整数。
    """
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


def _counter(
    program: Program, limit: int, native_modules: Collection[str] = frozenset()
) -> Callable[[tuple[Instruction, ...]], int]:
    """构造按模块记忆化的递归计步函数，供 ``events`` 与预算检查复用。"""
    modules, cached = program.module_map, {}

    def count(nodes: tuple[Instruction, ...]) -> int:
        """逐节点按加权规则累计步数，总数达到 ``limit + 1`` 后钳制。"""
        total = 0
        for node in nodes:
            if isinstance(node, Primitive):
                cost = max(4, sum(r.width for r in node.operands) ** 2)
            elif isinstance(node, (Load, Store)):
                cost = 1
            elif isinstance(node, Call):
                if node.module in native_modules:
                    cost = 1
                    total = min(limit + 1, total + cost)
                    continue
                if node.module not in cached:
                    cached[node.module] = count(
                        cast("tuple[Instruction, ...]", modules[node.module].body)
                    )
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


def expanded_steps(
    program: Program,
    limit: int = 1_000_000,
    native_modules: Collection[str] = frozenset(),
) -> int:
    """估计程序完全展开后的执行步数，用于执行前的预算检查。

    原语按操作数总位宽加权，Load/Store 计 1 步；模块调用递归计入被调体，Repeat 按次数相乘。
    计数一旦达到上限即被钳制在 ``limit + 1``，不再继续增长。

    Args:
        program: RIR 程序。
        limit: 计数钳制上限。
        native_modules: 按单步计的原生模块名集合；入口本身为原生模块时直接返回 1。

    Returns:
        int: 展开步数的饱和估计值。
    """
    if program.entry in native_modules:
        return 1
    return _counter(program, limit, native_modules)(
        cast("tuple[Instruction, ...]", program.main.body)
    )


@dataclass(frozen=True)
class LocalEnter:
    """局部寄存器进入作用域的事件，寄存器以零值加入执行状态。

    Attributes:
        name: 分配给局部寄存器的唯一合成名。
        type: 局部寄存器的类型。
    """

    name: str
    type: RegType


@dataclass(frozen=True)
class LocalExit:
    """局部寄存器离开作用域的事件，退出时要求全部分支复净为零。

    Attributes:
        name: 被释放的局部寄存器名。
        width: 局部寄存器位宽。
    """

    name: str
    width: int


def remap(ref: Ref, mapping: Mapping[str, Ref]) -> Ref:
    """按寄存器名替换表重写引用，用于模块调用时把形参绑定到实参视图。

    Args:
        ref: 待重写的 ``Ref``。
        mapping: 寄存器名到 ``Ref`` 的映射；原引用的每段截取替换结果的对应区间后顺序拼接。

    Returns:
        Ref: 重写后的新引用，类型与原引用一致。
    """
    parts: list[Span] = []
    for span in ref.parts:
        parts.extend(mapping[span.register][span.start : span.start + span.width].parts)
    return Ref(tuple(parts), ref.type)


def events(
    program: Program,
    *,
    max_steps: int = 1_000_000,
    native_modules: Collection[str] = frozenset(),
) -> Iterator[
    tuple[
        Primitive | Load | Store | LocalEnter | LocalExit | NativeSite,
        tuple[tuple[Ref, int], ...],
        bool,
    ]
]:
    """调用逐层展开为迭代器，原始 RIR 保持不变。

    Args:
        program: 待展开执行的 RIR 程序；无原生模块时须闭合。
        max_steps: 展开后的指令步数预算上限。
        native_modules: 视为原生执行、不再内联展开的模块名集合。

    Returns:
        Iterator: 惰性产出 ``(事件, 控制链, 逆序标志)`` 三元组；事件为
        基元、QRAM 读写、局部寄存器进出或原生调用点，控制链为
        ``(Ref, int)`` 元组，逆序标志表示处于 ``Adjoint`` 语境。
    """
    from pyqecclang.infrastructure.ir import Span

    program = validate(program, require_closed=not bool(native_modules))
    counter = _counter(program, max_steps, native_modules)
    if expanded_steps(program, max_steps, native_modules) > max_steps:
        raise ValidationError("执行超过展开预算")
    modules = program.module_map

    def walk(
        nodes: tuple[Instruction, ...],
        mapping: Mapping[str, Ref],
        resources: Mapping[str, str],
        controls: tuple[tuple[Ref, int], ...] = (),
        inverse: bool = False,
    ) -> Iterator[
        tuple[
            Primitive | Load | Store | LocalEnter | LocalExit | NativeSite,
            tuple[tuple[Ref, int], ...],
            bool,
        ]
    ]:
        """按序重写指令节点并惰性产出事件三元组，控制链与逆序标志沿递归传递。"""
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
            elif isinstance(node, Store):
                yield (
                    Store(
                        resources[node.resource],
                        remap(node.address, mapping),
                        remap(node.data, mapping),
                    ),
                    controls,
                    inverse,
                )

    local_counter = 0

    def enter_module(
        module: Module,
        mapping: Mapping[str, Ref],
        resources: Mapping[str, str],
        controls: tuple[tuple[Ref, int], ...] = (),
        inverse: bool = False,
    ) -> Iterator[
        tuple[
            Primitive | Load | Store | LocalEnter | LocalExit | NativeSite,
            tuple[tuple[Ref, int], ...],
            bool,
        ]
    ]:
        """进入单个模块：绑定形参映射、按需合成局部寄存器名并转发体事件流。"""
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


def gate_matrix(
    op: str, angle: float | None = None, inverse: bool = False
) -> tuple[tuple[complex, complex], tuple[complex, complex]]:
    """返回单比特门的标准 2x2 西矩阵。

    Args:
        op: 门名，取 ``h``、``x``、``y``、``z``、``s``、``t``、``phase``、``rx``、``ry``、``rz`` 之一。
        angle: 旋转与相位门的角度（弧度）；``s`` 与 ``t`` 使用固定角度，忽略该参数。
        inverse: 为真时返回共轭转置，即门的逆矩阵。

    Returns:
        tuple: 2x2 复数矩阵，以行元组的形式给出。

    Raises:
        ValidationError: 门名不在支持列表中。
    """
    if op == "h":
        a = 1 / math.sqrt(2)
        matrix: tuple[tuple[complex, complex], tuple[complex, complex]] = ((a, a), (a, -a))
    elif op == "x":
        matrix = ((0, 1), (1, 0))
    elif op == "y":
        matrix = ((0, -1j), (1j, 0))
    elif op == "z":
        matrix = ((1, 0), (0, -1))
    elif op in {"s", "t", "phase"}:
        theta = {"s": math.pi / 2, "t": math.pi / 4}.get(op, angle)
        matrix = ((1, 0), (0, cmath.exp(1j * cast(float, theta))))
    elif op == "rz":
        matrix = (
            (cmath.exp(-0.5j * cast(float, angle)), 0),
            (0, cmath.exp(0.5j * cast(float, angle))),
        )
    elif op in {"rx", "ry"}:
        c, s = math.cos(cast(float, angle) / 2), math.sin(cast(float, angle) / 2)
        matrix = ((c, -1j * s), (-1j * s, c)) if op == "rx" else ((c, -s), (s, c))
    else:
        raise ValidationError(f"未知单比特矩阵：{op}")
    if inverse:
        return cast(
            "tuple[tuple[complex, complex], tuple[complex, complex]]",
            tuple(tuple(complex(matrix[j][i]).conjugate() for j in range(2)) for i in range(2)),
        )
    return matrix


@dataclass(frozen=True)
class RegisterState:
    """参考执行的终态：入口寄存器布局与稀疏振幅。

    Attributes:
        registers: 入口模块的寄存器元组，规定取值元组的顺序。
        amplitudes: 寄存器取值元组到复振幅的映射。
    """

    registers: tuple
    amplitudes: dict[tuple[int, ...], complex]

    def statevector(self, max_qubits: int = 20) -> list[complex]:
        """把稀疏态打包为密集态向量，各寄存器按声明顺序 LSB-first 占据下标位段。

        Args:
            max_qubits: 允许的最大总位数。

        Returns:
            list: 长度为 ``2**总位数`` 的复振幅列表。

        Raises:
            ValidationError: 寄存器总位数超过 ``max_qubits``。
        """
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


def simulate(
    program: Program,
    memory: Mapping[str, Sequence[int] | Mapping[int, int]] | None = None,
    *,
    initial: Mapping[str, int] | None = None,
    max_steps: int = 1_000_000,
    max_states: int = 65536,
) -> RegisterState:
    """在稀疏寄存器态上参考执行封闭程序，只依赖标准库。

    初态是各寄存器取给定整数值的单基矢（缺省全零），态保存为寄存器整数值元组到复振幅的稀疏
    字典。指令流由 ``events`` 惰性展开：单比特门对首个操作数逐位作用，``xor``、``swap``、
    ``add_const`` 与 ``gphase`` 按 RIR 语义更新寄存器取值或幅值；``Control`` 条件按逻辑合取
    限定作用分支，``Adjoint`` 逆序取逆，``Repeat`` 按次数重复。QRAM 的 ``Load`` 是 XOR 读，
    叠加地址天然支持；``Store`` 是随机写，要求地址与数据处于确定基矢。局部寄存器进入时置零，
    退出时检查全部分支复净。

    Args:
        program: 待执行的封闭 RIR 程序。
        memory: 按资源名提供的 QRAM 数据，格式同 ``check_memory``。
        initial: 寄存器名到初始整数值的映射，缺省为全零。
        max_steps: 展开步数预算，超限报错。
        max_states: 稀疏基矢数预算，超限报错。

    Returns:
        RegisterState: 入口寄存器布局与终态稀疏振幅；幅值不超过 1e-15 的分量被截去。

    Raises:
        ValidationError: 程序或 QRAM 数据校验失败、初态寄存器未知或越界、Store 不在确定基矢、
            局部寄存器未复净，或超出步数与稀疏态数量预算。
    """
    program = validate(program, require_closed=True)
    memories = check_memory(program, memory)
    registers = program.main.registers
    index = {r.name: i for i, r in enumerate(registers)}
    initial = initial or {}
    if set(initial) - set(index):
        raise ValidationError("初态包含未知寄存器")
    values: list[int] | tuple[int, ...] = []
    for reg in registers:
        value = initial.get(reg.name, 0)
        if type(value) is not int or not 0 <= value < (1 << reg.type.width):
            raise ValidationError("初态寄存器值越界")
        cast("list[int]", values).append(value)
    state = {tuple(values): 1 + 0j}

    def read(ref: Ref, values: tuple[int, ...]) -> int:
        """按视图分段低位到高位拼接，从取值元组读出无符号整数。"""
        result, offset = 0, 0
        for span in ref.parts:
            result |= (
                (values[index[span.register]] >> span.start) & ((1 << span.width) - 1)
            ) << offset
            offset += span.width
        return result

    def write(ref: Ref, values: tuple[int, ...], value: int) -> tuple[int, ...]:
        """把整数值按分段掩码写入取值元组，返回新的取值元组。"""
        result, offset = list(values), 0
        for span in ref.parts:
            pos = index[span.register]
            mask = ((1 << span.width) - 1) << span.start
            result[pos] = (result[pos] & ~mask) | (
                ((value >> offset) & ((1 << span.width) - 1)) << span.start
            )
            offset += span.width
        return tuple(result)

    for node, controls, inverse in cast(
        "Iterator[tuple[Primitive | Load | Store | LocalEnter | LocalExit, tuple[tuple[Ref, int], ...], bool]]",
        events(program, max_steps=max_steps),
    ):
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

        def active(
            values: tuple[int, ...], controls: tuple[tuple[Ref, int], ...] = controls
        ) -> bool:
            """判断基矢是否满足全部控制条件的比较值。"""
            return all(read(ref, values) == expected for ref, expected in controls)

        if isinstance(node, Store):
            observed = {
                (read(node.address, values), read(node.data, values)) for values in state
            }
            if len(observed) != 1:
                raise ValidationError("Store 要求地址与数据寄存器在执行时处于确定基矢")
            address, value = observed.pop()
            if value:
                memories[node.resource][address] = value
            else:
                memories[node.resource].pop(address, None)
            continue

        if isinstance(node, Primitive) and node.op not in {"xor", "swap", "add_const", "gphase"}:
            matrix = gate_matrix(node.op, node.angle, inverse)
            for register, bit in locations(node.operands[0]):
                result: dict[tuple[int, ...], complex] = {}
                pos = index[register]
                for values, amplitude in state.items():
                    if not active(values):
                        result[values] = result.get(values, 0j) + amplitude
                        continue
                    incoming = (values[pos] >> bit) & 1
                    for outgoing in range(2):
                        target: list[int] | tuple[int, ...] = list(values)
                        cast("list[int]", target)[pos] = (
                            (target[pos] & ~(1 << bit)) | (outgoing << bit)
                        )
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
                    amplitude *= cmath.exp(1j * cast(float, node.angle) * (-1 if inverse else 1))
                elif node.op == "xor":
                    source, dest = node.operands
                    target = write(dest, values, read(source, values) ^ read(dest, values))
                elif node.op == "swap":
                    a, b = node.operands
                    av, bv = read(a, values), read(b, values)
                    target = write(b, write(a, values, bv), av)
                elif node.op == "add_const":
                    ref = node.operands[0]
                    value = read(ref, values) + cast(int, node.value) * (-1 if inverse else 1)
                    target = write(ref, values, value & ((1 << ref.width) - 1))
            result[target] = result.get(target, 0j) + amplitude
        state = {k: v for k, v in result.items() if abs(v) > 1e-15}
    return RegisterState(registers, state)
