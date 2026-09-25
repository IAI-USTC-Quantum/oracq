"Registration of PySparQ module-level implementations; runtime capabilities do not masquerade as gate-level closure."

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from types import ModuleType
from typing import cast

from oracq.infrastructure.builder import Operation
from oracq.infrastructure.ir import Module, Program, Ref, ValidationError
from oracq.infrastructure.linking import calls
from oracq.infrastructure.validation import validate


@dataclass(frozen=True)
class NativeSite:
    """A module call site in the execution flow taken over by a native implementation.

    Attributes:
        module: declaration of the called module.
        arguments: argument views corresponding one-to-one with the module registers.
        resources: resource binding names corresponding one-to-one with the module resources.
    """

    module: Module
    arguments: tuple[Ref, ...]
    resources: tuple[str, ...]


@dataclass(frozen=True)
class NativeEntry:
    """A native implementation registration entry under one module name.

    Attributes:
        module: module description declared at registration, used to compare against same-named modules in programs.
        factory: implementation factory only called with a ``NativeContext`` argument at execution time.
        label: implementation name displayed in execution reports.
        kind: implementation category tag, ``pysparq_custom`` by default.
    """

    module: Module
    factory: Callable
    label: str
    kind: str


class NativeRegistry:
    """Register PySparQ module-level native implementations by module name and verify coverage against programs."""

    def __init__(self) -> None:
        """Initialize an empty registration table from module names to native implementations."""
        self.entries: dict[str, NativeEntry] = {}

    def register(
        self,
        operation: Operation | Module,
        factory: Callable[[NativeContext], object],
        *,
        label: str | None = None,
        kind: str = "pysparq_custom",
    ) -> NativeRegistry:
        """Register the native implementation of one module.

        Args:
            operation: ``Operation`` or ``Module`` providing the module declaration.
            factory: factory called with a ``NativeContext`` argument at execution time; heavy work such as
                C++ compilation is deferred inside the factory as needed, not done at registration.
            label: display name of the implementation in execution reports; the module name is used when omitted.
            kind: implementation category tag.

        Returns:
            NativeRegistry: returns itself to support chained registration.

        Raises:
            ValidationError: the module name is already registered with different content.
        """
        module = operation.module if isinstance(operation, Operation) else operation
        entry = NativeEntry(module, factory, label or module.name, kind)
        if module.name in self.entries and self.entries[module.name] != entry:
            raise ValidationError("duplicate native implementation name")
        self.entries[module.name] = entry
        return self

    def matching(self, program: Program) -> frozenset[str]:
        """Compute the set of module names in a program taken over by native implementations.

        Args:
            program: program to verify.

        Returns:
            frozenset[str]: module names that match a registered entry by name and module declaration.

        Raises:
            ValidationError: a same-named module does not match the declaration at registration.
        """
        result = set()
        for module in program.modules:
            if module.name in self.entries:
                if self.entries[module.name].module != module:
                    raise ValidationError(f"native implementation does not match the current module description: {module.name}")
                result.add(module.name)
        return frozenset(result)

    def missing(self, program: Program) -> tuple[str, ...]:
        """List entry-reachable open declaration module names not yet covered by native implementations.

        Args:
            program: program to check.

        Returns:
            tuple[str, ...]: open declaration module names sorted by name; an empty tuple when there is no gap.

        Raises:
            ValidationError: the program failed RIR structural validation.
        """
        validate(program)
        supported = self.matching(program)
        result, visited = [], set()

        def visit(key: str) -> None:
            """Skip taken-over modules and recursively collect entry-reachable open declaration names."""
            if key in visited or key in supported:
                return
            visited.add(key)
            module = program.module_map[key]
            if module.body is None:
                result.append(key)
            else:
                for call in calls(module.body):
                    visit(call.module)

        visit(program.entry)
        return tuple(sorted(result))


@dataclass
class NativeContext:
    """Execution environment passed to native factories, bridging RIR formal parameters and PySparQ runtime objects.

    Attributes:
        ps: currently imported ``pysparq`` module.
        site: native call site being executed.
        names: mapping from RIR register names to PySparQ register names.
        qrams: mapping from resource names to PySparQ QRAM objects.
        memories: mapping from resource names to sparse word data dicts.
    """

    ps: object
    site: NativeSite
    names: dict
    qrams: dict
    memories: dict

    def register_id(self, parameter: str) -> int:
        """Resolve a quantum formal parameter of the module to a PySparQ integer register id.

        Args:
            parameter: module register formal parameter name.

        Returns:
            id given by PySparQ ``System.get_id`` for the corresponding register.

        Raises:
            ValidationError: the argument is not a single view covering exactly one whole register, such as a slice or concatenation.
        """
        index = next(i for i, r in enumerate(self.site.module.registers) if r.name == parameter)
        ref = self.site.arguments[index]
        if len(ref.parts) != 1:
            raise ValidationError("the current native factory requires whole-register arguments")
        span = ref.parts[0]
        name = self.names[span.register]
        # ps is the imported pysparq module object; the System attribute comes from the dynamic module.
        if span.start != 0 or span.width != cast(ModuleType, self.ps).System.size_of(name):
            raise ValidationError("the current native factory does not accept slices; use whole registers or gate-level implementations")
        return cast(ModuleType, self.ps).System.get_id(name)

    def qram(self, parameter: str) -> object:
        """Resolve a resource formal parameter of the module to the QRAM object bound at this call site.

        Args:
            parameter: module resource formal parameter name.

        Returns:
            PySparQ QRAM object bound to this formal parameter.
        """
        index = next(i for i, r in enumerate(self.site.module.resources) if r.name == parameter)
        return self.qrams[self.site.resources[index]]


class DynamicCppFactory:
    """Calls the real compile_operator on demand; C++ code is not written into RIR."""

    def __init__(
        self,
        name: str,
        source: str,
        parameters: Iterable[str],
        *,
        cache_dir: str,
        base_class: str = "SelfAdjointOperator",
    ) -> None:
        """Record the operator name, C++ source and parameters, deferring compilation to the first call."""
        self.name: str = name
        self.source: str = source
        self.parameters: tuple[str, ...] = tuple(parameters)
        self.cache_dir: str = str(cache_dir)
        self.base_class: str = base_class
        self._class: type | None = None

    def __call__(self, context: NativeContext) -> object:
        """Lazily compile the C++ operator source on first call, then construct an instance from formal register ids."""
        if self._class is None:
            from pysparq.dynamic_operator import compile_operator

            self._class = compile_operator(
                name=self.name,
                cpp_code=self.source,
                base_class=self.base_class,
                constructor_args=[("size_t", name) for name in self.parameters],
                cache_dir=self.cache_dir,
            )
        return self._class(**{name: context.register_id(name) for name in self.parameters})
