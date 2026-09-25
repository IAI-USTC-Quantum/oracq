"Checkable input contracts for the generation stage; structure/capabilities and mathematical truth are handled separately."

from __future__ import annotations

import math
from collections.abc import Callable, Iterable
from dataclasses import asdict, dataclass
from typing import NoReturn, TypeVar, cast

from oracq.infrastructure.builder import Operation
from oracq.infrastructure.ir import Register, Resource, ValidationError

_T = TypeVar("_T")


@dataclass(frozen=True)
class OracleCapabilities:
    """Adjoint and controlled invocation capabilities supported by an oracle; values are derived by conservative conjunction over the dependency graph.

    Attributes:
        adjoint: Whether adjoint invocation is supported.
        controlled: Whether controlled invocation is supported.
    """

    adjoint: bool
    controlled: bool


@dataclass(frozen=True)
class OracleSpec:
    """Snapshot read from the actual operation; CKS is a collection of operations with no unified total anc width."""

    type: str
    main_qubit: int | None
    anc_qubit: int | None
    capabilities: OracleCapabilities
    implementation: str
    registers: tuple[Register, ...] = ()
    resources: tuple[Resource, ...] = ()
    alpha: float | None = None
    parameters: tuple[tuple[str, str | int | float | bool], ...] = ()
    components: tuple[tuple[str, OracleSpec], ...] = ()

    def to_dict(self) -> dict[str, object]:
        """Convert the snapshot into a JSON-storable plain dictionary.

        Returns:
            dict[str, object]: Nested plain dictionary with fields of the same names.
        """
        return asdict(self)


class OracleView:
    """Wrapper classes share read-only queries; widths always follow the RIR register signature.

    Host view classes must declare the wrapped ``operation`` and the
    ``oracle_kind`` class marker.
    """

    operation: Operation
    oracle_kind: str

    @property
    def spec(self) -> OracleSpec:
        """Current ``OracleSpec`` snapshot of this view; re-derived on each access."""
        return self.describe()

    def describe(self) -> OracleSpec:
        """Derive the ``OracleSpec`` of this view from the RIR module signature, attributes, and dependency graph.

        Returns:
            OracleSpec: The description snapshot of this view.
        """
        return _operation_spec(self.operation, self.oracle_kind)

    @property
    def type(self) -> str:
        """Descriptive paradigm role name such as ``block_encoding``; not used as a unique type or dispatch basis."""
        return self.spec.type

    @property
    def main_qubit(self) -> int | None:
        """Number of bits of the main register (``target``/``address``/``column``); ``None`` when it cannot be recognized."""
        return self.spec.main_qubit

    @property
    def anc_qubit(self) -> int | None:
        """Number of bits of the public auxiliary register (``signal``/``work``); ``None`` when the operation collection has no unified width."""
        return self.spec.anc_qubit

    @property
    def capabilities(self) -> OracleCapabilities:
        """Effective invocation capabilities derived from the whole dependency graph."""
        return self.spec.capabilities


@dataclass(frozen=True)
class ContractIssue:
    """A concrete issue found by one contract check.

    Attributes:
        code: Issue category code, e.g. ``INPUT_PROTOCOL``, ``INPUT_WIDTH``.
        path: Input path where the issue is located, e.g. ``protocol.input.field``.
        expected: The value expected by the contract.
        actual: The value actually received.
        message: User-facing description of the issue.
    """

    code: str
    path: str
    expected: object
    actual: object
    message: str


class ContractError(ValidationError):
    """Exception that raises all discovered contract issues at once; a ``ValidationError`` subclass.

    Args:
        issues: Iterable of ``ContractIssue`` objects joined as ``code path: message``.

    Attributes:
        issues: Immutable tuple of structured issues.
    """

    def __init__(self, issues: Iterable[ContractIssue]) -> None:
        """Construct the exception with all discovered issues and join their messages one by one."""
        self.issues: tuple[ContractIssue, ...] = tuple(issues)
        super().__init__("; ".join(f"{i.code} {i.path}: {i.message}" for i in self.issues))


def fail(code: str, path: str, expected: object, actual: object, message: str) -> NoReturn:
    """Raise a ``ContractError`` immediately with a single issue.

    Args:
        code: Issue category code.
        path: Input path where the issue is located.
        expected: The value expected by the contract.
        actual: The value actually received.
        message: User-facing description of the issue.

    Returns:
        NoReturn: Never returns normally, always exits by raising.

    Raises:
        ContractError: Always raised, carrying only this one issue.
    """
    raise ContractError((ContractIssue(code, path, expected, actual, message),))


def require_instance(value: object, cls: type[_T], path: str) -> _T:
    """Require ``value`` to be an instance of ``cls``, otherwise raise ``ContractError``.

    Args:
        value: The value to check.
        cls: The expected type.
        path: Input path used when reporting the issue.

    Returns:
        _T: The ``value`` itself once the instance check passes, narrowed to ``cls``.

    Raises:
        ContractError: Reported as ``INPUT_TYPE`` when the instance check fails.
    """
    if not isinstance(value, cls):
        fail("INPUT_TYPE", path, cls.__name__, type(value).__name__, f"requires an instance of {cls.__name__}")
    return value


def finite_real(
    value: object,
    path: str,
    *,
    minimum: float | None = None,
    strict: bool = False,
) -> int | float:
    """Require ``value`` to be a finite real number with an optional lower bound; ``bool`` is not accepted.

    Args:
        value: The numeric value to check.
        path: Input path used when reporting the issue.
        minimum: The minimum value allowed; when omitted no lower bound is checked.
        strict: When ``True``, require strictly greater than ``minimum``; otherwise equality is allowed.

    Returns:
        int | float: The ``value`` itself once the check passes, narrowed to a finite real number.

    Raises:
        ContractError: Reported as ``CONFIG_VALUE`` when the value is not a finite real number or is below the bound.
    """
    if type(value) not in (int, float) or not math.isfinite(cast("int | float", value)):
        fail("CONFIG_VALUE", path, "finite real", repr(value), "requires a finite real number, not a bool")
    if minimum is not None and (
        cast("int | float", value) <= minimum if strict else cast("int | float", value) < minimum
    ):
        fail("CONFIG_VALUE", path, f"{'>' if strict else '>='}{minimum}", value, "out of the allowed range")
    return cast("int | float", value)


def positive_integer(
    value: object,
    path: str,
    *,
    minimum: int = 1,
    maximum: int | None = None,
) -> int:
    """Require ``value`` to be an integer within the range; ``bool`` is not accepted.

    Args:
        value: The integer to check.
        path: Input path used when reporting the issue.
        minimum: The minimum value allowed, default 1.
        maximum: The maximum value allowed; when omitted no upper bound is checked.

    Returns:
        int: The ``value`` itself once the check passes, narrowed to ``int``.

    Raises:
        ContractError: Reported as ``CONFIG_VALUE`` when the value is not an integer or is out of range.
    """
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        fail("CONFIG_VALUE", path, f"integer {minimum}..{maximum}", repr(value), "integer configuration out of range")
    return cast(int, value)


def validate_signature(operation: object, fields: Iterable[str], path: str) -> None:
    """Validate that the public registers of the operation are exactly the agreed ``bits`` register set.

    Args:
        operation: The ``Operation`` to check.
        fields: The agreed set of register names.
        path: Input path used when reporting the issue.

    Raises:
        ContractError: Raised when the input is not an ``Operation``, when the
        register names disagree with ``fields``, or when a non-``bits``
        register exists.
        ValidationError: Passed through when program merging validation fails for ``operation``.
    """
    from oracq.infrastructure.builder import Operation

    require_instance(operation, Operation, path)
    cast("Operation", operation).program()
    actual = {r.name: r.type for r in cast("Operation", operation).module.registers}
    if set(actual) != set(fields) or any(t.kind != "bits" for t in actual.values()):
        fail("ORACLE_SIGNATURE", path, list(fields), list(actual), "requires the agreed bits register signature")


def _operation_spec(operation: Operation, kind: str | None = None) -> OracleSpec:
    """Derive the ``OracleSpec`` snapshot of an operation from module attributes and the dependency graph."""
    from oracq.infrastructure.linking import capabilities, unresolved

    program = operation.program()
    module = operation.module
    attrs = dict(module.attributes)
    caps = capabilities(program)
    regs = {r.name: r.type.width for r in module.registers}
    kind = kind or cast("str", attrs.get("oracle_paradigm", "unitary"))
    main = regs.get("target", regs.get("address", regs.get("column")))
    anc = regs.get("signal", regs.get("work"))
    return OracleSpec(
        str(kind),
        main,
        anc,
        OracleCapabilities(caps["supports_adjoint"], caps["supports_controlled"]),
        "open" if unresolved(program) else "closed",
        module.registers,
        module.resources,
        cast("float | None", attrs.get("be_alpha")),
        module.attributes,
    )


def describe_oracle(value: object) -> OracleSpec:
    """Obtain the ``OracleSpec`` description snapshot of an arbitrary input object.

    An ``Operation`` is derived directly from module attributes and the
    dependency graph; an object with a callable ``describe`` method has that
    method invoked and must return an ``OracleSpec``; any other host object
    only reports the type name and optional ``main_qubit``/``anc_qubit``
    attributes, without consulting any global type table.

    Args:
        value: An ``Operation``, an object implementing ``describe``, or any host object.

    Returns:
        OracleSpec: The read-only description snapshot of the input.

    Raises:
        ContractError: Reported as ``INPUT_TYPE`` when the ``describe`` return
        value is not an ``OracleSpec``.
        ValidationError: Passed through when program merging validation fails in the ``Operation`` branch.
    """
    from oracq.infrastructure.builder import Operation

    if isinstance(value, Operation):
        return _operation_spec(value)
    describe = getattr(value, "describe", None)
    if callable(describe):
        result = describe()
        require_instance(result, OracleSpec, "describe.output")
        return result
    # Custom algorithm objects need not register with the language; without a describe method only the host type is reported.
    return OracleSpec(
        type(value).__name__,
        getattr(value, "main_qubit", None),
        getattr(value, "anc_qubit", None),
        OracleCapabilities(False, False),
        "host",
    )


def requires(value: _T, protocol: type, *, path: str = "input") -> _T:
    """Algorithm-side structural check; the protocol is defined by the caller, and no global type table is consulted.

    Args:
        value: The input value to check.
        protocol: The protocol class required to be satisfied.
        path: Input path used when reporting the issue, ``input`` by default.

    Returns:
        _T: The ``value`` itself once the protocol check passes.
    """
    if not isinstance(value, protocol):
        fail(
            "INPUT_PROTOCOL",
            path,
            protocol.__name__,
            type(value).__name__,
            "requires satisfying " + protocol.__name__,
        )
    return value


@dataclass(frozen=True)
class InputRequirement:
    """Declaration of a protocol's requirement on a single named input.

    ``protocols`` lists the acceptable protocol classes (satisfying any one is
    enough); the optional ``adapter`` converts the raw value into a concrete
    view before checking layout, capabilities, and promises.

    Attributes:
        name: Input name, unique within the owning ``AlgorithmContract``.
        protocols: Tuple of acceptable protocol classes.
        adapter: Optional adapter function such as ``as_block_encoding``; when omitted the raw value is checked directly.
        main_qubit: Expected number of bits of the main register; ``None`` means no check.
        anc_qubit: Expected number of bits of the auxiliary register; ``None`` means no check.
        alpha: Expected BE normalization constant; ``None`` means no check.
        adjoint: Require the input to support adjoint invocation.
        controlled: Require the input to support controlled invocation.
        zero_input: Require the input to explicitly promise zero-input preparation.
        clean_work: Require the input to explicitly promise the work space is restored to zero.

    Raises:
        ValidationError: Raised at construction when the name is empty, the
        protocol tuple is invalid, the adapter is not callable, or a field
        value is illegal.
    """

    name: str
    protocols: tuple[type, ...]
    adapter: Callable[..., object] | None = None
    main_qubit: int | None = None
    anc_qubit: int | None = None
    alpha: float | None = None
    adjoint: bool = False
    controlled: bool = False
    zero_input: bool = False
    clean_work: bool = False

    def __post_init__(self) -> None:
        """Validate the name, protocol tuple, adapter callability, and values of the constraint fields."""
        if (
            not self.name
            or type(self.protocols) is not tuple
            or not self.protocols
            or any(not isinstance(t, type) for t in self.protocols)
        ):
            raise ValidationError("InputRequirement requires a name and a non-empty tuple of protocol classes")
        if self.adapter is not None and not callable(self.adapter):
            raise ValidationError("InputRequirement.adapter must be callable")
        for key in ("main_qubit", "anc_qubit"):
            value = getattr(self, key)
            if value is not None:
                positive_integer(value, key, minimum=0)
        if self.alpha is not None:
            finite_real(self.alpha, "alpha", minimum=0, strict=True)
        for key in ("adjoint", "controlled", "zero_input", "clean_work"):
            if type(getattr(self, key)) is not bool:
                raise ValidationError(f"InputRequirement.{key} must be a bool")

    def to_dict(self) -> dict[str, object]:
        """Convert into a JSON-storable dictionary; protocols record class names, the adapter records the function name.

        Returns:
            dict[str, object]: Plain dictionary of the requirement fields and constraint values.
        """
        return {
            "name": self.name,
            "protocols": [p.__name__ for p in self.protocols],
            "adapter": getattr(self.adapter, "__name__", None),
            **{
                key: getattr(self, key)
                for key in (
                    "main_qubit",
                    "anc_qubit",
                    "alpha",
                    "adjoint",
                    "controlled",
                    "zero_input",
                    "clean_work",
                )
            },
        }

    def inspect(
        self,
        value: object,
        *,
        prefix: str = "",
    ) -> tuple[OracleSpec | None, tuple[ContractIssue, ...]]:
        """Check whether a candidate value satisfies this requirement, returning the description and issues instead of raising.

        Args:
            value: The input value to check.
            prefix: Path prefix for reporting, usually the name of the owning protocol.

        Returns:
            tuple: ``(spec, issues)``. ``spec`` is ``None`` when the protocol
        does not match or obtaining the description fails; layout,
        capability, and promise issues are appended one by one to ``issues``
        without blocking later checks.

        Raises:
            Exception: Exceptions other than ValidationError from inside the provider are passed through unchanged, preserving the cause.
        """
        _, spec, issues = self.resolve(value, prefix=prefix)
        return spec, issues

    def resolve(
        self, value: object, *, prefix: str = ""
    ) -> tuple[object, OracleSpec | None, tuple[ContractIssue, ...]]:
        """Adapt and check once, also keeping the actual view for use by the algorithm kernel."""
        path = f"{prefix}.{self.name}" if prefix else self.name
        expected_names = tuple(p.__name__ for p in self.protocols)
        if not any(isinstance(value, p) for p in self.protocols):
            return value, None, (
                ContractIssue(
                    "INPUT_PROTOCOL",
                    path,
                    expected_names,
                    type(value).__name__,
                    f"requires satisfying at least one protocol: {expected_names}",
                ),
            )
        try:
            value = self.adapter(value) if self.adapter is not None else value
            spec = describe_oracle(value)
        except ContractError as exc:
            return value, None, tuple(
                ContractIssue(
                    "INPUT_ADAPTER" if issue.code == "INPUT_TYPE" else issue.code,
                    path + "." + issue.path, issue.expected, issue.actual, issue.message,
                ) for issue in exc.issues
            )
        except ValidationError as exc:
            return value, None, (
                ContractIssue(
                    "INPUT_ADAPTER", path, expected_names, type(value).__name__, str(exc)
                ),
            )
        issues: list[ContractIssue] = []

        def issue(code: str, field: str, expected: object, actual: object, message: str) -> None:
            """Construct an issue and append it with the path ``path.field``."""
            issues.append(ContractIssue(code, path + "." + field, expected, actual, message))

        for field in ("main_qubit", "anc_qubit", "alpha"):
            expected, actual = getattr(self, field), getattr(spec, field)
            if expected is not None and expected != actual:
                issue(
                    "INPUT_LAYOUT" if field != "alpha" else "INPUT_ALPHA",
                    field,
                    expected,
                    actual,
                    f"expected {expected}, received {actual}",
                )
        for field in ("adjoint", "controlled"):
            if getattr(self, field) and not getattr(spec.capabilities, field):
                issue("INPUT_CAPABILITY", field, True, False, f"missing supports_{field}")
        attrs = dict(spec.parameters)
        for field in ("zero_input", "clean_work"):
            if getattr(self, field) and attrs.get(field) is not True:
                issue("INPUT_PROMISE", field, True, attrs.get(field), f"requires an explicit {field} promise")
        return value, spec, tuple(issues)


@dataclass(frozen=True)
class ContractReport:
    """Result report of one contract check, covering only structure and declarations, proving no mathematical property.

    Attributes:
        protocol: Name of the checked protocol.
        inputs: ``(input name, OracleSpec)`` tuples recording successfully obtained input descriptions.
        issues: All issues found by the check; empty means passed.
        assumptions: Assumptions declared by the protocol, kept as-is.
        adapters: Adapter names involved in the check; empty by default.
    """

    protocol: str
    inputs: tuple[tuple[str, OracleSpec], ...]
    issues: tuple[ContractIssue, ...] = ()
    assumptions: tuple[str, ...] = ()
    adapters: tuple[str, ...] = ()

    @property
    def ok(self) -> bool:
        """``True`` when no issue was found."""
        return not self.issues

    def require(self) -> ContractReport:
        """Return self when passed, otherwise raise all issues at once.

        Returns:
            ContractReport: Self, for chaining.

        Raises:
            ContractError: Raised when any unresolved issue exists.
        """
        if self.issues:
            raise ContractError(self.issues)
        return self

    def to_dict(self) -> dict[str, object]:
        """Convert into a JSON-storable dictionary, adding the top-level ``ok`` field.

        Returns:
            dict[str, object]: Plain dictionary containing ``ok`` and all report fields.
        """
        return {"ok": self.ok, **asdict(self)}


@dataclass(frozen=True)
class AlgorithmContract:
    """The algorithm's contract declaration over all inputs, checkable as a whole with an aggregated report.

    Attributes:
        name: Protocol name.
        inputs: Tuple of ``InputRequirement`` distinguished by name.
        same_width: Pairs of input names required to share the same ``main_qubit``.
        assumptions: Assumption declarations passed through to the report as-is.

    Raises:
        ValidationError: Raised at construction when the name is empty, input
        names are duplicated, an inputs element is illegal, or a same-width
        constraint references a nonexistent input.
    """

    name: str
    inputs: tuple[InputRequirement, ...]
    same_width: tuple[tuple[str, str], ...] = ()
    assumptions: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Validate the inputs type, name uniqueness, and same-width constraint references."""
        if type(self.inputs) is not tuple or any(
            not isinstance(r, InputRequirement) for r in self.inputs
        ):
            raise ValidationError("AlgorithmContract.inputs requires an immutable tuple of InputRequirement")
        names = [r.name for r in self.inputs]
        if not self.name or len(names) != len(set(names)):
            raise ValidationError("invalid protocol name or duplicated input names")
        if any(len(pair) != 2 or any(n not in names for n in pair) for pair in self.same_width):
            raise ValidationError("a same-width constraint references a nonexistent input")

    def check(self, **values: object) -> ContractReport:
        """Check the inputs one by one by name and generate an aggregated report.

        Only methods of the input objects are invoked to obtain descriptions;
        no algorithm kernel or quantum simulation is run.

        Args:
            **values: Input values keyed by requirement name; extras and missing values are both recorded as issues.

        Returns:
            ContractReport: Aggregated check result over all inputs.
        """
        return self.resolve(**values).report

    def resolve(self, **values: object) -> ResolvedInputs:
        """Obtain and check one concrete set of inputs; the report can be saved, and the views are used only in the host generation stage."""
        resolved: dict[str, object] = {}
        specs: dict[str, OracleSpec] = {}
        issues: list[ContractIssue] = []
        expected = {r.name for r in self.inputs}
        for key in sorted(values.keys() - expected):
            issues.append(
                ContractIssue("INPUT_UNEXPECTED", f"{self.name}.{key}", None, key, "unexpected input")
            )
        for req in self.inputs:
            if req.name not in values:
                issues.append(
                    ContractIssue(
                        "INPUT_MISSING",
                        f"{self.name}.{req.name}",
                        tuple(p.__name__ for p in req.protocols),
                        None,
                        "missing input",
                    )
                )
                continue
            value, spec, found = req.resolve(values[req.name], prefix=self.name)
            if spec is not None:
                resolved[req.name] = value
                specs[req.name] = spec
            issues.extend(found)
        for left, right in self.same_width:
            if (
                left in specs
                and right in specs
                and specs[left].main_qubit != specs[right].main_qubit
            ):
                issues.append(
                    ContractIssue(
                        "INPUT_WIDTH",
                        f"{self.name}.{right}.main_qubit",
                        specs[left].main_qubit,
                        specs[right].main_qubit,
                        f"must match the target width of {left}",
                    )
                )
        report = ContractReport(self.name, tuple(specs.items()), tuple(issues), self.assumptions)
        return ResolvedInputs(tuple(resolved.items()), report)

    def to_dict(self) -> dict[str, object]:
        """Convert into a JSON-storable dictionary; each input recursively uses its own ``to_dict``.

        Returns:
            dict[str, object]: Plain dictionary of the contract fields and per-input requirements.
        """
        return {
            "name": self.name,
            "inputs": [r.to_dict() for r in self.inputs],
            "same_width": self.same_width,
            "assumptions": self.assumptions,
        }


@dataclass(frozen=True)
class ResolvedInputs:
    """Views and check report of one boundary adaptation; host objects are not cached, and this is not serialized as executable IR."""

    values: tuple[tuple[str, object], ...]
    report: ContractReport

    def get(self, name: str, cls: type[_T]) -> _T:
        """After the contract passes, obtain the checked concrete view by name."""
        self.report.require()
        return require_instance(dict(self.values)[name], cls, self.report.protocol + "." + name)


# Compatibility for old imports; new algorithms use names that express the responsibility precisely.
ProtocolContract = AlgorithmContract
