"""Shared modular solve circuits, memory binding, and optional native execution."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal

from oracq.infrastructure.backends.originir import OriginIRArtifact, export_originir, run_originir
from oracq.infrastructure.backends.pysparq import run_pysparq
from oracq.infrastructure.estimate import ResourceEstimate, estimate_resources
from oracq.infrastructure.execution import RegisterState, check_memory
from oracq.infrastructure.ir import Call, Instruction, Load, Program, ValidationError
from oracq.infrastructure.native import NativeRegistry
from oracq.infrastructure.qram_schema import RegisteredQRAM
from oracq.infrastructure.serialization import dumps


def resolve_memory(program: Program, banks: tuple[RegisteredQRAM, ...]) -> dict[str, dict[int, int]]:
    """Trace resource arguments through calls without expanding gates or repeats."""
    registered: dict[str, RegisteredQRAM] = {}
    for bank in banks:
        if bank.name in registered and bank != registered[bank.name]:
            raise ValidationError("conflicting registered QRAM snapshots for resource " + bank.name)
        registered[bank.name] = bank
    memories: dict[str, dict[int, int]] = {}
    seen: set[tuple[str, tuple[tuple[str, str], ...]]] = set()
    modules = program.module_map

    def visit(module_name: str, mapping: Mapping[str, str]) -> None:
        """Walk one module instance under a resource mapping, once per mapping."""
        key = (module_name, tuple(sorted(mapping.items())))
        if key in seen:
            return
        seen.add(key)
        module = modules[module_name]
        types = {r.name: r.type for r in module.resources}

        def walk(body: tuple[Instruction, ...]) -> None:
            """Inspect one instruction body, recursing into nested structures."""
            for node in body:
                if isinstance(node, Load):
                    if node.resource not in registered:
                        raise ValidationError("no registered memory snapshot for leaf QRAM resource " + node.resource)
                    bank = registered[node.resource]
                    if bank.type != types[node.resource]:
                        raise ValidationError("registered QRAM widths disagree with the leaf resource " + node.resource)
                    target, cells = mapping[node.resource], bank.snapshot()
                    if target in memories and memories[target] != cells:
                        raise ValidationError("multiple QRAM snapshots resolve to the same entry resource " + target)
                    memories[target] = cells
                elif isinstance(node, Call):
                    visit(node.module, {r.name: mapping[actual] for r, actual in zip(modules[node.module].resources, node.resources, strict=True)})
                elif hasattr(node, "body"):
                    walk(node.body)

        if module.body is not None:
            walk(module.body)

    visit(program.entry, {r.name: r.name for r in program.main.resources})
    # Resources may be declared by an oracle whose query is unreachable in this
    # particular specialization, e.g. a degree-zero evolution polynomial.
    for resource in program.main.resources:
        if resource.name not in memories:
            raise ValidationError("an entry QRAM resource has no reachable registered query: " + resource.name)
    return check_memory(program, memories)


@dataclass(frozen=True)
class SolveCircuit:
    """A modular RIR circuit with separately snapshotted execution memory."""

    program: Program
    _memory: tuple[tuple[str, tuple[tuple[int, int], ...]], ...] = ()

    @property
    def memory(self) -> dict[str, dict[int, int]]:
        """Copy entry bindings, including all resource argument renamings."""
        return {name: dict(cells) for name, cells in self._memory}

    def dumps(self, *, format: Literal["yaml", "json"] = "yaml") -> str:
        """Serialize RIR without memory contents or module expansion."""
        return dumps(self.program, format=format)

    def originir_ext(self) -> OriginIRArtifact:
        """Export modular OriginIR-ext without importing a native backend."""
        return export_originir(self.program)

    def resource_estimate(self, basic_gates: str = "clifford+t+qram") -> ResourceEstimate:
        """Estimate circuit resources compositionally in the requested basis.

        Clifford+T uses the exact seven-T Toffoli decomposition: four T,
        three T-dagger, six CNOT, and two H gates. Arbitrary rotations remain
        synthesis estimates in ``ResourceEstimate.t_total(epsilon)``. QRAM
        queries remain independent oracle costs.
        """
        if basic_gates not in {"clifford+t+qram", "toffoli+clifford+t+qram"}:
            raise ValidationError("supported resource bases are clifford+t+qram and toffoli+clifford+t+qram")
        result = estimate_resources(self.program)
        if basic_gates == "clifford+t+qram":
            count = result.atoms.pop("toffoli", 0)
            if count:
                for atom, increment in {"t": 4 * count, "tdg": 3 * count, "cnot": 6 * count, "h": 2 * count}.items():
                    result.atoms[atom] = result.atoms.get(atom, 0) + increment
        return result

    def run_pysparq(
        self, *, max_steps: int = 1_000_000, max_states: int = 65536,
        native_registry: NativeRegistry | None = None,
    ) -> RegisterState:
        """Execute the prepared RIR and memory with the optional PySparQ backend."""
        return run_pysparq(self.program, self.memory, max_steps=max_steps, max_states=max_states, native_registry=native_registry)

    def run_originir_ext(self, *, max_qubits: int = 24, max_steps: int = 1_000_000) -> RegisterState:
        """Execute exported OriginIR and decode its original register boundaries."""
        vector = run_originir(self.program, self.memory, max_qubits=max_qubits, max_steps=max_steps)
        amplitudes = {}
        for index, amplitude in enumerate(vector):
            if abs(amplitude) <= 1e-15:
                continue
            remaining, values = index, []
            for register in self.program.main.registers:
                values.append(remaining & ((1 << register.type.width) - 1))
                remaining >>= register.type.width
            amplitudes[tuple(values)] = amplitude
        return RegisterState(self.program.main.registers, amplitudes)
