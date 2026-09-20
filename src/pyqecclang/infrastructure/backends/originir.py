"直接生成模块化 OriginIR-ext，控制以 inline controlled_by 表达。"

from __future__ import annotations

import hashlib
from dataclasses import dataclass

from pyqecclang.infrastructure.ir import (
    Adjoint,
    Call,
    Control,
    Load,
    Primitive,
    Program,
    Repeat,
    Store,
    ValidationError,
)
from pyqecclang.infrastructure.validation import validate


@dataclass(frozen=True)
class OriginIRArtifact:
    """模块化 OriginIR-ext 导出的文本结果及配套映射信息。

    Attributes:
        text: 完整的 OriginIR-ext 源文本，含 QRAMDECL、全部 DEF 定义和入口调用。
        registers: 入口寄存器名到其占用全局量子位下标元组的映射。
        resources: RIR 资源名到导出文本中 QRAM 名称的映射。
        workspace_qubits: 入口私有工作区占用的全局量子位下标元组，默认为空。
    """

    text: str
    registers: dict[str, tuple[int, ...]]
    resources: dict[str, str]
    workspace_qubits: tuple[int, ...] = ()


def export_originir(program: Program) -> OriginIRArtifact:
    """把闭合程序导出为模块化 OriginIR-ext 文本。

    纯文本导出，不需要安装任何量子后端。模块调用与 Repeat 保留为可复用的
    DEF 定义及调用，不在导出阶段展开；控制以逐门 ``controlled_by`` 和附加
    控制形式参数表达。

    Args:
        program: 待导出的闭合程序。

    Returns:
        OriginIRArtifact: 导出文本及寄存器、资源与工作区的映射信息。

    Raises:
        ValidationError: 程序结构非法或存在未绑定的 oracle。
    """
    return _Exporter(validate(program, require_closed=True)).run()


class _Exporter:
    def __init__(self, program):
        self.program = program
        self.modules = program.module_map
        self.definitions = []
        self.cache = {}
        from pyqecclang.infrastructure.layout import workspace_table

        self.workspace = workspace_table(program)

    def _symbol(self, key):
        digest = hashlib.sha256(repr(key).encode()).hexdigest()[:24]
        return f"m_{key[1]}_{digest}"

    def _mapping(self, module):
        result = {
            reg.name: tuple(f"v_{reg.name}[{i}]" for i in range(reg.type.width))
            for reg in module.registers
        }
        cursor = 0
        for reg in module.locals:
            result[reg.name] = tuple(
                f"pw_work[{i}]" for i in range(cursor, cursor + reg.type.width)
            )
            cursor += reg.type.width
        return result

    def _workspace_args(self, module):
        return [f"pw_work[{i}]" for i in range(self.workspace[module.name])]

    def _bits(self, ref, mapping):
        return [
            bit
            for span in ref.parts
            for bit in mapping[span.register][span.start : span.start + span.width]
        ]

    def _args(self, module, mapping):
        return [bit for reg in module.registers for bit in mapping[reg.name]]

    def _control_formals(self, count):
        return tuple((f"pc_control[{i}]", 1) for i in range(count))

    def _header(self, symbol, module, control_count):
        args = [f"v_{r.name}[{r.type.width}]" for r in module.registers if r.type.width]
        if self.workspace[module.name]:
            args.append(f"pw_work[{self.workspace[module.name]}]")
        if control_count:
            args.append(f"pc_control[{control_count}]")
        return f"DEF {symbol}({', '.join(args)})"

    def _wrap(self, lines, controls):
        """仅包装普通门或 QRAM 调用；合并已有的内层逐门控制。"""
        if not lines or not controls:
            return lines
        zeros = [bit for bit, value in controls if value == 0]
        result = [f"X {bit}" for bit in zeros]
        for line in lines:
            head, separator, tail = line.partition(" controlled_by (")
            existing = tail[:-1].split(", ") if separator else []
            combined = list(dict.fromkeys([bit for bit, _ in controls] + existing))
            result.append(head + " controlled_by (" + ", ".join(combined) + ")")
        result.extend(f"X {bit}" for bit in reversed(zeros))
        return result

    def _call(self, symbol, args, controls):
        """模块控制通过附加形式参数传递，定义和调用都不展开。"""
        zeros = [bit for bit, value in controls if value == 0]
        full_args = args + [bit for bit, _ in controls]
        return (
            [f"X {bit}" for bit in zeros]
            + [f"{symbol}({', '.join(full_args)})"]
            + [f"X {bit}" for bit in reversed(zeros)]
        )

    def module(self, module, resources, control_count=0):
        key = ("module", module.name, tuple(sorted(resources.items())), control_count)
        if key in self.cache:
            return self.cache[key]
        symbol = self._symbol(key)
        self.cache[key] = symbol
        mapping = self._mapping(module)
        body = self.body(
            module.body, module, mapping, resources, self._control_formals(control_count)
        )
        self.definitions.append([self._header(symbol, module, control_count), *body, "ENDDEF"])
        return symbol

    def repeat(self, nodes, count, module, resources, control_count):
        key = ("repeat", module.name, tuple(sorted(resources.items())), nodes, count, control_count)
        if key in self.cache:
            return self.cache[key]
        symbol = self._symbol(key)
        self.cache[key] = symbol
        mapping = self._mapping(module)
        controls = self._control_formals(control_count)
        args = self._args(module, mapping) + self._workspace_args(module)
        if count == 1:
            lines = self.body(nodes, module, mapping, resources, controls)
        else:
            half = self.repeat(nodes, count // 2, module, resources, control_count)
            lines = self._call(half, args, controls) + self._call(half, args, controls)
            if count % 2:
                base = self.repeat(nodes, 1, module, resources, control_count)
                lines.extend(self._call(base, args, controls))
        self.definitions.append([self._header(symbol, module, control_count), *lines, "ENDDEF"])
        return symbol

    def body(self, nodes, module, mapping, resources, controls=()):
        lines = []
        for node in nodes:
            if isinstance(node, Control):
                bits = self._bits(node.register, mapping)
                extra = tuple((bit, (node.value >> i) & 1) for i, bit in enumerate(bits))
                lines.extend(self.body(node.body, module, mapping, resources, controls + extra))
            elif isinstance(node, Adjoint):
                inner = self.body(node.body, module, mapping, resources, controls)
                if inner:
                    lines.extend(["DAGGER", *inner, "ENDDAGGER"])
            elif isinstance(node, Repeat):
                if node.count and node.body:
                    symbol = self.repeat(node.body, node.count, module, resources, len(controls))
                    lines.extend(
                        self._call(
                            symbol,
                            self._args(module, mapping) + self._workspace_args(module),
                            controls,
                        )
                    )
            elif isinstance(node, Call):
                target = self.modules[node.module]
                bound = {
                    formal.name: resources[actual]
                    for formal, actual in zip(target.resources, node.resources, strict=True)
                }
                symbol = self.module(target, bound, len(controls))
                args = [bit for ref in node.arguments for bit in self._bits(ref, mapping)]
                offset = sum(r.type.width for r in module.locals)
                args += [
                    f"pw_work[{i}]" for i in range(offset, offset + self.workspace[target.name])
                ]
                lines.extend(self._call(symbol, args, controls))
            elif isinstance(node, Load):
                args = self._bits(node.address, mapping) + self._bits(node.data, mapping)
                lines.extend(
                    self._wrap([f"{resources[node.resource]} {', '.join(args)}"], controls)
                )
            elif isinstance(node, Store):
                args = self._bits(node.address, mapping) + self._bits(node.data, mapping)
                lines.append(f"QRAMWRITE {resources[node.resource]} {', '.join(args)}")
            elif isinstance(node, Primitive):
                lines.extend(self.primitive(node, module, mapping, controls))
        return lines

    def primitive(self, node, module, mapping, controls):
        operands = [self._bits(ref, mapping) for ref in node.operands]
        if node.op == "gphase":
            angle = repr(float(node.angle))
            if controls:
                zeros = [bit for bit, value in controls if not value]
                active = tuple((bit, 1) for bit, _ in controls[:-1])
                target = controls[-1][0]
                return (
                    [f"X {bit}" for bit in zeros]
                    + self._wrap([f"U1 {target}, ({angle})"], active)
                    + [f"X {bit}" for bit in reversed(zeros)]
                )
            anchor = self._args(module, mapping)[0]
            return [
                f"U1 {anchor}, ({angle})",
                f"X {anchor}",
                f"U1 {anchor}, ({angle})",
                f"X {anchor}",
            ]
        if node.op in {"xor", "swap"}:
            gate = "CNOT" if node.op == "xor" else "SWAP"
            raw = [f"{gate} {a}, {b}" for a, b in zip(operands[0], operands[1], strict=True)]
        elif node.op == "add_const":
            bits, raw = operands[0], []
            for offset in range(len(bits)):
                if (node.value >> offset) & 1:
                    for i in reversed(range(offset + 1, len(bits))):
                        raw.extend(
                            self._wrap([f"X {bits[i]}"], tuple((b, 1) for b in bits[offset:i]))
                        )
                    raw.append(f"X {bits[offset]}")
        else:
            gate = "U1" if node.op == "phase" else node.op.upper()
            suffix = "" if node.angle is None else f", ({repr(float(node.angle))})"
            raw = [f"{gate} {bit}{suffix}" for bit in operands[0]]
        return self._wrap(raw, controls)

    def run(self):
        entry = self.program.main
        resources = {res.name: "ram_" + res.name for res in entry.resources}
        registers, cursor = {}, 0
        for reg in entry.registers:
            registers[reg.name] = tuple(range(cursor, cursor + reg.type.width))
            cursor += reg.type.width
        private = tuple(range(cursor, cursor + self.workspace[entry.name]))
        cursor += len(private)
        main = self.module(entry, resources)
        text = [
            f"QRAMDECL {resources[r.name]} {r.type.address_width},{r.type.data_width}"
            for r in entry.resources
        ]
        text.extend([f"QINIT {cursor}", "CREG 0"])
        for definition in self.definitions:
            text.extend(definition)
        text.append(f"{main}({', '.join(f'q[{i}]' for i in range(cursor))})")
        return OriginIRArtifact("\n".join(text) + "\n", registers, resources, private)


def run_originir(program: Program, memory=None, *, max_qubits=24, max_steps=1_000_000):
    """通过实际 UnifiedQuantum 后端执行小规模实例；依赖为可选安装。"""
    from pyqecclang.infrastructure.execution import check_memory, expanded_steps
    from pyqecclang.infrastructure.linking import uses_store

    program = validate(program, require_closed=True)
    if uses_store(program):
        raise ValidationError("UnifiedQuantum 执行暂不支持运行期 QRAM 写（Store）；文本导出不受影响")
    width = sum(r.type.width for r in program.main.registers)
    from pyqecclang.infrastructure.layout import workspace_table

    if width + workspace_table(program)[program.entry] > max_qubits:
        raise ValidationError("超过 OriginIR 状态向量的本次量子位预算")
    if expanded_steps(program, max_steps) > max_steps:
        raise ValidationError("超过 UnifiedQuantum 解析展开预算")
    memories = check_memory(program, memory)
    for resource in program.main.resources:
        if resource.type.address_width + resource.type.data_width > 30:
            raise ValidationError("当前 UnifiedQuantum QRAM 容器要求地址位和数据位总数不超过 30")
        if resource.type.address_width > 20:
            raise ValidationError("当前适配器拒绝物化超过 2^20 项的 QRAM")
    try:
        from uniqc.simulator import Simulator
    except ImportError as exc:
        raise ValidationError(
            "OriginIR 执行需要已安装 UnifiedQuantum 的环境；文本导出不需要该依赖"
        ) from exc
    artifact = export_originir(program)
    sim = Simulator(least_qubit_remapping=False)
    sim.simulate_preprocess(artifact.text)
    for resource, cells in memories.items():
        ram = sim.qram_objects[artifact.resources[resource]]
        for address, value in cells.items():
            ram.write(address, value)
    vector = sim.simulate_statevector(artifact.text)
    if artifact.workspace_qubits:
        if any(abs(value) > 1e-12 for value in vector[1 << width :]):
            raise ValidationError("OriginIR 局部工作区未复净")
        return vector[: 1 << width]
    return vector
