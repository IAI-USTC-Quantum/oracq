"Directly emit modular OriginIR-ext, with controls expressed via inline controlled_by."

from __future__ import annotations

import hashlib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TypeAlias, cast

from oracq.infrastructure.ir import (
    Adjoint,
    Call,
    Control,
    Instruction,
    Load,
    Module,
    Primitive,
    Program,
    Ref,
    Repeat,
    Store,
    ValidationError,
)
from oracq.infrastructure.validation import validate


@dataclass(frozen=True)
class OriginIRArtifact:
    """Text result of a modular OriginIR-ext export, with companion mapping information.

    Attributes:
        text: The complete OriginIR-ext source text, including QRAMDECL, all DEF definitions and the entry call.
        registers: Mapping from entry register names to their tuples of occupied global qubit indices.
        resources: Mapping from RIR resource names to QRAM names in the exported text.
        workspace_qubits: Tuple of global qubit indices occupied by the entry's private workspace, empty by default.
    """

    text: str
    registers: dict[str, tuple[int, ...]]
    resources: dict[str, str]
    workspace_qubits: tuple[int, ...] = ()


def export_originir(program: Program) -> OriginIRArtifact:
    """Export a closed program as modular OriginIR-ext text.

    Pure text export, requiring no quantum backend installation. Module
    calls and Repeats are kept as reusable DEF definitions and calls, not
    expanded at export time; controls are expressed via per-gate
    ``controlled_by`` and additional control formal parameters.

    Args:
        program: The closed program to export.

    Returns:
        OriginIRArtifact: The exported text plus register, resource and workspace mapping information.

    Raises:
        ValidationError: The program structure is illegal, or an unbound oracle exists.
    """
    return _Exporter(validate(program, require_closed=True)).run()


_DefinitionKey: TypeAlias = (
    tuple[str, str, tuple[tuple[str, str], ...], int]
    | tuple[str, str, tuple[tuple[str, str], ...], tuple[Instruction, ...], int, int]
)
"""Key of the definition cache: a ``("module", ...)`` module key or a ``("repeat", ...)`` repeat-body key."""


class _Exporter:
    """Translate a closed program into an OriginIR-ext definition table with per-module memoization."""

    def __init__(self, program: Program) -> None:
        """Bind the program, module table and definition cache, and precompute workspace quotas."""
        self.program: Program = program
        self.modules: dict[str, Module] = program.module_map
        self.definitions: list[list[str]] = []
        self.cache: dict[_DefinitionKey, str] = {}
        from oracq.infrastructure.layout import workspace_table

        self.workspace: dict[str, int] = workspace_table(program)

    def _symbol(self, key: _DefinitionKey) -> str:
        """Generate a deterministic ``DEF`` symbol name from the cache key."""
        digest = hashlib.sha256(repr(key).encode()).hexdigest()[:24]
        return f"m_{key[1]}_{digest}"

    def _mapping(self, module: Module) -> dict[str, tuple[str, ...]]:
        """Build the mapping of all module register names to tuples of wire names."""
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

    def _workspace_args(self, module: Module) -> list[str]:
        """Generate the list of workspace wire actual arguments passed at module call sites."""
        return [f"pw_work[{i}]" for i in range(self.workspace[module.name])]

    def _bits(self, ref: Ref, mapping: Mapping[str, tuple[str, ...]]) -> list[str]:
        """Expand a view into a list of OriginIR wire names."""
        return [
            bit
            for span in ref.parts
            for bit in mapping[span.register][span.start : span.start + span.width]
        ]

    def _args(self, module: Module, mapping: Mapping[str, tuple[str, ...]]) -> list[str]:
        """Concatenate the module's actual argument wires in public register declaration order."""
        return [bit for reg in module.registers for bit in mapping[reg.name]]

    def _control_formals(self, count: int) -> tuple[tuple[str, int], ...]:
        """Generate the tuple of ``pc_control`` control formal parameters whose value is always 1."""
        return tuple((f"pc_control[{i}]", 1) for i in range(count))

    def _header(self, symbol: str, module: Module, control_count: int) -> str:
        """Generate the ``DEF`` line, appending workspace and control formal parameters as needed."""
        args = [f"v_{r.name}[{r.type.width}]" for r in module.registers if r.type.width]
        if self.workspace[module.name]:
            args.append(f"pw_work[{self.workspace[module.name]}]")
        if control_count:
            args.append(f"pc_control[{control_count}]")
        return f"DEF {symbol}({', '.join(args)})"

    def _wrap(self, lines: list[str], controls: tuple[tuple[str, int], ...]) -> list[str]:
        """Wrap only plain gates or QRAM calls; merge any existing inner per-gate controls."""
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

    def _call(
        self, symbol: str, args: list[str], controls: tuple[tuple[str, int], ...]
    ) -> list[str]:
        """Module controls are passed through additional formal parameters; neither definitions nor calls expand them."""
        zeros = [bit for bit, value in controls if value == 0]
        full_args = args + [bit for bit, _ in controls]
        return (
            [f"X {bit}" for bit in zeros]
            + [f"{symbol}({', '.join(full_args)})"]
            + [f"X {bit}" for bit in reversed(zeros)]
        )

    def module(
        self, module: Module, resources: Mapping[str, str], control_count: int = 0
    ) -> str:
        """Generate or reuse the module's ``DEF`` definition for the given control count, returning the symbol name.

        Args:
            module: The module definition to translate.
            resources: Mapping from module resource names to exported QRAM names.
            control_count: Number of extra control bits the definition must accommodate; no controls when omitted.

        Returns:
            str: The ``DEF`` definition's symbol name; a same-key module is generated once and reused directly on repeated calls.
        """
        key = ("module", module.name, tuple(sorted(resources.items())), control_count)
        if key in self.cache:
            return self.cache[key]
        symbol = self._symbol(key)
        self.cache[key] = symbol
        mapping = self._mapping(module)
        # Exporting the entry requires a closed program; module.body is verified non-None by require_closed.
        body = self.body(
            cast("tuple[Instruction, ...]", module.body),
            module,
            mapping,
            resources,
            self._control_formals(control_count),
        )
        self.definitions.append([self._header(symbol, module, control_count), *body, "ENDDEF"])
        return symbol

    def repeat(
        self,
        nodes: tuple[Instruction, ...],
        count: int,
        module: Module,
        resources: Mapping[str, str],
        control_count: int,
    ) -> str:
        """Generate or reuse a ``DEF`` definition for a repeat body by the binary method, returning the symbol name.

        Args:
            nodes: The repeat body's instruction sequence.
            count: The repetition count, a positive integer.
            module: The module definition the repeat body belongs to.
            resources: Mapping from module resource names to exported QRAM names.
            control_count: Number of extra control bits the definition must accommodate.

        Returns:
            str: The repeat body ``DEF`` definition's symbol name; generated once per body and count.
        """
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

    def body(
        self,
        nodes: tuple[Instruction, ...],
        module: Module,
        mapping: Mapping[str, tuple[str, ...]],
        resources: Mapping[str, str],
        controls: tuple[tuple[str, int], ...] = (),
    ) -> list[str]:
        """Translate an instruction body node by node into OriginIR-ext text lines.

        Args:
            nodes: The instruction sequence to translate.
            module: The module the instruction body belongs to, providing local registers and the call context.
            mapping: Mapping from register names to tuples of wire names.
            resources: Mapping from module resource names to exported QRAM names.
            controls: (wire, active value) control tuples inherited from outer ``Control`` nodes.

        Returns:
            list[str]: The sequence of translated OriginIR-ext text lines.
        """
        lines: list[str] = []
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

    def primitive(
        self,
        node: Primitive,
        module: Module,
        mapping: Mapping[str, tuple[str, ...]],
        controls: tuple[tuple[str, int], ...],
    ) -> list[str]:
        """Lower a single primitive instruction into OriginIR-ext gate text lines.

        Args:
            node: The primitive instruction to lower.
            module: The module the instruction belongs to; an uncontrolled gphase anchors on its first actual argument wire.
            mapping: Mapping from register names to tuples of wire names.
            controls: (wire, active value) control tuples inherited from outer ``Control`` nodes.

        Returns:
            list[str]: The gate text lines of the primitive, with controls merged into ``controlled_by``.
        """
        operands = [self._bits(ref, mapping) for ref in node.operands]
        if node.op == "gphase":
            # Validation guarantees the gphase angle is a finite real number.
            angle = repr(float(cast(float, node.angle)))
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
                # Validation guarantees the integer value of add_const exists.
                if (cast(int, node.value) >> offset) & 1:
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

    def run(self) -> OriginIRArtifact:
        """Export the entry module and assemble the full OriginIR-ext text and mapping information.

        Returns:
            OriginIRArtifact: The full exported text, plus the mapping of registers,
            resource names and the workspace to global qubits.
        """
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


def run_originir(
    program: Program,
    memory: Mapping[str, Sequence[int] | Mapping[int, int]] | None = None,
    *,
    max_qubits: int = 24,
    max_steps: int = 1_000_000,
) -> list[complex]:
    """Execute a small instance through a real UnifiedQuantum backend; the dependency is an optional install.

    Args:
        program: The closed RIR program to execute; it must not contain runtime QRAM writes via Store.
        memory: Binding from resource names to QRAM contents; sequences give data words by index, mappings by address.
        max_qubits: The qubit budget cap for state-vector simulation, workspace included.
        max_steps: The instruction-step budget cap after expansion.

    Returns:
        list[complex]: The final-state amplitudes over the entry's public register space, ordered from least to most significant bit.
    """
    from oracq.infrastructure.execution import check_memory, expanded_steps
    from oracq.infrastructure.linking import uses_store

    program = validate(program, require_closed=True)
    if uses_store(program):
        raise ValidationError("UnifiedQuantum execution does not yet support runtime QRAM writes via Store; text export is unaffected")
    width = sum(r.type.width for r in program.main.registers)
    from oracq.infrastructure.layout import workspace_table

    if width + workspace_table(program)[program.entry] > max_qubits:
        raise ValidationError("the qubit budget for this OriginIR state-vector run was exceeded")
    if expanded_steps(program, max_steps) > max_steps:
        raise ValidationError("the UnifiedQuantum parse expansion budget was exceeded")
    memories = check_memory(program, memory)
    for resource in program.main.resources:
        if resource.type.address_width + resource.type.data_width > 30:
            raise ValidationError("the current UnifiedQuantum QRAM container requires at most 30 total address and data bits")
        if resource.type.address_width > 20:
            raise ValidationError("the current adapter refuses to materialize QRAMs with more than 2^20 entries")
    try:
        from uniqc.simulator import Simulator
    except ImportError as exc:
        raise ValidationError(
            "OriginIR execution requires an environment with UnifiedQuantum installed; text export does not need that dependency"
        ) from exc
    artifact = export_originir(program)
    sim = Simulator(least_qubit_remapping=False)
    sim.simulate_preprocess(artifact.text)
    for resource, cells in memories.items():  # type: ignore[assignment]
        # resource is bound as a Resource object in the loop above; here it carries the resource name string.
        ram = sim.qram_objects[artifact.resources[cast(str, resource)]]
        for address, value in cells.items():
            ram.write(address, value)
    vector = sim.simulate_statevector(artifact.text)
    if artifact.workspace_qubits:
        if any(abs(value) > 1e-12 for value in vector[1 << width :]):
            raise ValidationError("the OriginIR local workspace was not restored to zero")
        return vector[: 1 << width]
    return vector
