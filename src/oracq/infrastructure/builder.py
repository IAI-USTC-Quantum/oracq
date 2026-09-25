"Python generation-stage builder. Generated results contain no Python callables."

from __future__ import annotations

from collections.abc import Callable, Iterator
from contextlib import AbstractContextManager, contextmanager
from dataclasses import dataclass
from typing import TYPE_CHECKING, TypeAlias

from oracq.infrastructure.ir import (
    QRAM,
    Adjoint,
    Call,
    Control,
    Instruction,
    Load,
    Module,
    Primitive,
    Program,
    Ref,
    Register,
    RegType,
    Repeat,
    Resource,
    Span,
    Store,
    ValidationError,
)
from oracq.infrastructure.validation import validate

if TYPE_CHECKING:
    from oracq.algorithms.input_model.operators import BlockEncoding
    from oracq.algorithms.input_model.oracles import StatePreparation
    from oracq.infrastructure.estimate import ResourceEstimate

RegisterSpec: TypeAlias = dict[str, RegType]
"""Register specification for the Builder constructor: a mapping from register names to storage types."""


@dataclass(frozen=True)
class Operation:
    """A complete RIR operation consisting of an entry module and its dependency modules.

    Attributes:
        module: entry module definition.
        dependencies: all dependency module definitions other than the entry.
    """

    module: Module
    dependencies: tuple[Module, ...] = ()

    def unitary(self) -> Operation:
        """Full-register unitary action of the RIR operation; does not claim any target matrix is itself unitary.

        Returns:
            Operation: returns itself unchanged; an RIR operation is itself unitary evolution on the full register space.
        """
        return self

    def state_preparation(self) -> StatePreparation:
        """Prepare U|0> with the full public register space as target.

        Returns:
            StatePreparation: state preparation object obtained by wrapping this operation via ``from_unitary``.
        """
        from oracq.algorithms.input_model.oracles import StatePreparation

        return StatePreparation.from_unitary(self)

    def block_encoding(self) -> BlockEncoding:
        """Trivial alpha=1, zero-signal block encoding of the full unitary itself; keeps the original module calls.

        Returns:
            BlockEncoding: trivial block encoding with scaling factor 1.0 and a zero-width signal register.
        """
        from oracq.algorithms.input_model.operators import _name, block_encoding
        from oracq.algorithms.input_model.oracles import invoke, resources_for
        from oracq.infrastructure.ir import Bits

        prep = self.state_preparation()
        b = Builder(
            _name("unitary_as_be", self),
            {"target": Bits(prep.width), "signal": Bits(0)},
            resources_for(("unitary", prep.operation)),
        )
        invoke(b, prep.operation, "unitary", target=b["target"], work=b["signal"])
        return block_encoding(b.finish(), alpha=1.0)

    @classmethod
    def from_program(cls, program: Program) -> Operation:
        """Recover the entry operation and module dependencies from a serialized program, keeping open declarations.

        Args:
            program: complete RIR program to recover; modules other than the entry become dependencies.

        Returns:
            Operation: operation whose entry is ``program.main`` and whose dependencies are all remaining modules.
        """
        validate(program)
        return cls(program.main, tuple(m for m in program.modules if m.name != program.entry))

    def program(self) -> Program:
        """Merge the entry module and dependencies into a validated ``Program`` with modules sorted by name.

        Returns:
            Program: program with ``module.name`` as its entry.

        Raises:
            ValidationError: conflicting definitions exist for a dependency module name, or the merged result violates RIR rules.
        """
        modules: dict[str, Module] = {}
        for item in (*self.dependencies, self.module):
            if item.name in modules and modules[item.name] != item:
                raise ValidationError(f"conflicting definitions for module name: {item.name}")
            modules[item.name] = item
        return validate(Program(self.module.name, tuple(modules[k] for k in sorted(modules))))

    def quantikz(self, **kwargs: bool | str | None) -> str:
        """Render the entry module as quantikz code; structure stays RIR with no call expansion.

        Returns:
            str: quantikz LaTeX source produced by rendering the entry module.
        """
        from oracq.infrastructure.backends.quantikz import quantikz

        # **kwargs forwards into keyword-only parameters; mypy cannot verify heterogeneous keyword packs.
        return quantikz(self, **kwargs)  # type: ignore[arg-type]

    def estimate(self, **kwargs: bool) -> ResourceEstimate:
        """Composable Toffoli+Clifford+T+QRAM level resource estimation with Repeat counts multiplied symbolically.

        Returns:
            ResourceEstimate: estimate summarized by Toffoli, Clifford, T and QRAM query counts.
        """
        from oracq.infrastructure.estimate import estimate_resources

        return estimate_resources(self.program(), **kwargs)


class Builder:
    """Generation-stage builder that constructs an RIR module instruction by instruction; generated results contain no Python callables.

    Args:
        name: name of the module to generate.
        registers: mapping from public register names to bit width types.
        resources: mapping from public QRAM resource names to specifications; omitted means the module has no resources.
        attributes: attribute key-value pairs appended to the module; omitted means empty.

    Public registers are accessed as ``Ref`` views via ``builder[name]``, ``local``
    declares module-private work registers, and ``finish`` closes the builder
    and returns a validated ``Operation``.
    """

    def __init__(
        self,
        name: str,
        registers: RegisterSpec,
        resources: dict[str, QRAM] | None = None,
        *,
        attributes: dict | None = None,
    ) -> None:
        """Declare the module signature and initialize the instruction frames, reference table and dependency table."""
        self.name: str = name
        self.attributes: tuple[tuple[str, str | int | float | bool], ...] = tuple(
            sorted((attributes or {}).items())
        )
        self.registers: tuple[Register, ...] = tuple(
            Register(k, v) for k, v in registers.items()
        )
        self.resources: tuple[Resource, ...] = tuple(
            Resource(k, v) for k, v in (resources or {}).items()
        )
        self._refs: dict[str, Ref] = {
            r.name: Ref((Span(r.name, 0, r.type.width),), r.type) for r in self.registers
        }
        self.locals: list[Register] = []
        self._frames: list[list[Instruction]] = [[]]
        self._dependencies: dict[str, Module] = {}
        self._closed: bool = False

    def local(self, name: str, type: RegType) -> Ref:
        """Declare a module-private work register and return the full view covering it.

        Args:
            name: register name; must not collide with a public register, a resource or another local register.
            type: type and bit width of the register.

        Returns:
            Ref: reference covering the entire new register.

        Raises:
            ValidationError: name conflict, or the builder is already closed.
        """
        if self._closed or name in self._refs or name in {r.name for r in self.resources}:
            raise ValidationError("local register name conflict or the builder is already closed")
        self.locals.append(Register(name, type))
        ref = Ref((Span(name, 0, type.width),), type)
        self._refs[name] = ref
        return ref

    def __getitem__(self, name: str) -> Ref:
        """Return the full-register ``Ref`` view for a public or local register name."""
        return self._refs[name]

    def emit(self, instruction: Instruction) -> None:
        """Append one RIR instruction to the current innermost scope.

        Args:
            instruction: RIR ``Instruction`` node.

        Raises:
            ValidationError: the builder is already closed.
        """
        if self._closed:
            raise ValidationError("the builder is already closed")
        self._frames[-1].append(instruction)

    def gate(self, name: str, target: Ref, angle: float | None = None) -> None:
        """Emit a single-target ``Primitive`` gate.

        Args:
            name: gate name.
            target: target view.
            angle: angle in radians for rotation gates; omitted for non-rotation gates.
        """
        self.emit(Primitive(name, (target,), angle=angle))

    def h(self, target: Ref) -> None:
        """Emit a Hadamard gate.

        Args:
            target: target view the gate acts on.
        """
        self.gate("h", target)

    def x(self, target: Ref) -> None:
        """Emit a Pauli-X gate.

        Args:
            target: target view the gate acts on.
        """
        self.gate("x", target)

    def z(self, target: Ref) -> None:
        """Emit a Pauli-Z gate.

        Args:
            target: target view the gate acts on.
        """
        self.gate("z", target)

    def ry(self, target: Ref, angle: float) -> None:
        """Emit a gate rotating by ``angle`` radians about the Y axis.

        Args:
            target: target view the rotation acts on.
            angle: rotation angle in radians.
        """
        self.gate("ry", target, angle)

    def rz(self, target: Ref, angle: float) -> None:
        """Emit a gate rotating by ``angle`` radians about the Z axis.

        Args:
            target: target view the rotation acts on.
            angle: rotation angle in radians.
        """
        self.gate("rz", target, angle)

    def xor(self, source: Ref, target: Ref) -> None:
        """Emit a bitwise XOR gate that XORs the value of ``source`` into the equal-width ``target``.

        Args:
            source: source view providing the value to XOR in.
            target: target view of the same width as ``source`` that accumulates the XOR result.
        """
        self.emit(Primitive("xor", (source, target)))

    def swap(self, first: Ref, second: Ref) -> None:
        """Emit a swap gate between the two equal-width views ``first`` and ``second``.

        Args:
            first: first view participating in the swap.
            second: second view of the same width as ``first``.
        """
        self.emit(Primitive("swap", (first, second)))

    def add_const(self, target: Ref, value: int) -> None:
        """Emit a constant-addition gate that adds ``value`` to ``target`` with wraparound by bit width.

        Args:
            target: target view that undergoes the addition and wraps by its bit width.
            value: integer constant to add.
        """
        self.emit(Primitive("add_const", (target,), value=value))

    def global_phase(self, angle: float) -> None:
        """Emit a global phase gate that multiplies the overall amplitude by ``exp(i*angle)``.

        Args:
            angle: overall phase angle in radians.
        """
        self.emit(Primitive("gphase", (), angle=angle))

    def qram(self, resource: str, address: Ref, data: Ref) -> None:
        """Emit a QRAM read: query ``resource`` at ``address`` and XOR the hit word into ``data``.

        Args:
            resource: name of the QRAM resource to query.
            address: view providing the address of the cell to query.
            data: equal-width target view into which the hit data word is XORed.
        """
        self.emit(Load(resource, address, data))

    def store(self, resource: str, address: Ref, data: Ref) -> None:
        """Emit a QRAM write: store ``data`` into cell ``address`` of ``resource``.

        Args:
            resource: name of the target QRAM resource to write.
            address: address view specifying the cell to write.
            data: source view providing the data to write.
        """
        self.emit(Store(resource, address, data))

    def call(
        self, operation: Operation, *, resources: dict[str, str] | None = None, **arguments: Ref
    ) -> None:
        """Call another ``Operation`` as-is; the callee stays a ``Call`` node and is not expanded at generation time.

        Args:
            operation: operation to call; its module and dependencies are registered as dependencies of this module.
            resources: mapping from callee resource formal parameter names to this module's resource names.
            **arguments: argument views keyed by callee register names.

        Raises:
            ValidationError: quantum or resource argument names do not match, or conflicting definitions exist for a dependency module name.
        """
        target = operation.module
        if set(arguments) != {r.name for r in target.registers}:
            raise ValidationError("module call quantum argument names do not match")
        resources = resources or {}
        if set(resources) != {r.name for r in target.resources}:
            raise ValidationError("module call resource argument names do not match")
        for module in (*operation.dependencies, target):
            existing = self._dependencies.get(module.name)
            if existing is not None and existing != module:
                raise ValidationError(f"conflicting definitions for module name: {module.name}")
            self._dependencies[module.name] = module
        self.emit(
            Call(
                target.name,
                tuple(arguments[r.name] for r in target.registers),
                tuple(resources[r.name] for r in target.resources),
            )
        )

    @contextmanager
    def _block(
        self, factory: Callable[[tuple[Instruction, ...]], Instruction]
    ) -> Iterator[Builder]:
        """Open a new instruction frame; on normal exit aggregate the frame's instructions through ``factory`` and emit them."""
        if self._closed:
            raise ValidationError("the builder is already closed")
        self._frames.append([])
        try:
            yield self
        except BaseException:
            self._frames.pop()
            raise
        else:
            nodes = tuple(self._frames.pop())
            self.emit(factory(nodes))

    def repeat(self, count: int) -> AbstractContextManager[Builder]:
        """Open a scope repeated ``count`` times; body instructions aggregate into a ``Repeat`` node whose count stays symbolic and is not expanded at generation time.

        Args:
            count: number of times the body scope repeats; stored symbolically in the ``Repeat`` node.

        Returns:
            scope context for use with a ``with`` statement; all blocks must be closed per their ``with`` nesting before ``finish``.
        """
        return self._block(lambda body: Repeat(count, body))

    def control(self, register: Ref, value: int | None = None) -> AbstractContextManager[Builder]:
        """Open a controlled scope; body instructions aggregate into a ``Control`` node.

        Args:
            register: control view.
            value: integer mask that activates the scope; omitted means the all-ones mask of the register.

        Returns:
            scope context for use with a ``with`` statement; all blocks must be closed per their ``with`` nesting before ``finish``.
        """
        if value is None:
            value = (1 << register.width) - 1
        return self._block(lambda body: Control(register, value, body))

    def adjoint(self) -> AbstractContextManager[Builder]:
        """Open an inverse scope; body instructions aggregate into an ``Adjoint`` block executed in reverse order with each operation inverted.

        Returns:
            scope context for use with a ``with`` statement; all blocks must be closed per their ``with`` nesting before ``finish``.
        """
        return self._block(Adjoint)

    def finish(self) -> Operation:
        """Finish construction: assemble the module with all registered dependencies, validate and return an ``Operation``, then close the builder.

        Returns:
            Operation: operation whose entry is the module built here and whose dependencies are the registered callees.

        Raises:
            ValidationError: the builder is already closed, scopes are still open, or the assembled result violates RIR rules.
        """
        if self._closed or len(self._frames) != 1:
            raise ValidationError("the builder is already closed or control blocks are still open")
        operation = Operation(
            Module(
                self.name,
                self.registers,
                self.resources,
                tuple(self._frames[0]),
                self.attributes,
                tuple(self.locals),
            ),
            tuple(self._dependencies[k] for k in sorted(self._dependencies)),
        )
        operation.program()
        self._closed = True
        return operation
