"Catalog of oracle paradigms and plain gate/QRAM implementations. Mathematical correctness certification is out of scope for this module."

from __future__ import annotations

import cmath
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from typing import cast

from oracq.algorithms.input_model.contracts import (
    OracleSpec,
    OracleView,
    fail,
    positive_integer,
    require_instance,
    validate_signature,
)
from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import (
    QRAM,
    Bits,
    Module,
    Ref,
    Register,
    RegType,
    Resource,
    ValidationError,
    fuse,
)

PARADIGMS = {
    "unitary": "quantum operation interface over the full register space",
    "block_encoding": "combined convention of the zero-signal projected-angle block and be_alpha",
    "database_xor": "|address,data> -> |address,data XOR memory[address]>",
    "state_prep_isometry": "prepares the target state from the zero-state subspace; adjoint and controlled use requires a reversible extension",
    "sparse_location_inplace": "CKS P_A(1): |column,index> -> |column,nu(column,index)>",
    "sparse_entry_xor": "CKS P_A(2): |row,column,data> -> |row,column,data XOR A[row,column]>",
    "phase_oracle": "applies a phase to target basis states per a defined predicate",
    "reversible_function": "domain computation interface with inputs preserved and outputs reversibly updated",
    "algorithm_stage": "an unfinished algorithm stage with an explicit register interface",
}
"""Built-in oracle role conventions; applications may use custom identifiers without modifying this catalog."""


def declare(
    name: str,
    registers: Mapping[str, RegType],
    *,
    paradigm: str = "unitary",
    resources: Mapping[str, QRAM] | None = None,
    attributes: Mapping[str, str | int | float | bool] | None = None,
    supports_adjoint: bool = True,
    supports_controlled: bool = True,
) -> Operation:
    """Declare an open oracle operation.

    Generates an ``Operation`` with an empty body and ``implementation_status``
    set to ``unresolved``; the structure is validated at generation time, and
    the implementation can be closed via ``annotate`` or binding.

    Args:
        name: Module name.
        registers: Ordered mapping from register names to types such as ``Bits``.
        paradigm: Built-in or application-defined role identifier; the role name itself proves no mathematical property.
        resources: Mapping from resource names to resource types such as ``QRAM``.
        attributes: Extra attributes merged into the module.
        supports_adjoint: Whether to declare support for the adjoint operation.
        supports_controlled: Whether to declare support for the controlled operation.

    Returns:
        Operation: The unresolved open declaration.

    Raises:
        ValidationError: The paradigm name is invalid, or the generated module fails structural validation.
    """
    from oracq.infrastructure.validation import name as validate_name

    validate_name(paradigm)
    attrs = dict(attributes or {})
    attrs.update(
        oracle_paradigm=paradigm,
        supports_adjoint=supports_adjoint,
        supports_controlled=supports_controlled,
        implementation_status="unresolved",
    )
    operation = Operation(
        Module(
            name,
            tuple(Register(k, v) for k, v in registers.items()),
            tuple(Resource(k, v) for k, v in (resources or {}).items()),
            None,
            tuple(sorted(attrs.items())),
        )
    )
    operation.program()
    return operation


def annotate(
    operation: Operation,
    paradigm: str,
    **attributes: str | int | float | bool,
) -> Operation:
    """Register the paradigm and attributes of a generated implementation operation.

    Sets ``implementation_status`` to ``constructed`` and writes ``oracle_paradigm``;
    adjoint/controlled capabilities explicitly restricted to ``False`` cannot be
    re-granted through annotation. When the paradigm is ``state_prep_isometry``,
    the ``zero_input`` and ``clean_work`` promises are added by default.

    Args:
        operation: The generated ``Operation``.
        paradigm: Built-in or application-defined role identifier.
        **attributes: Additional key-value pairs merged into the module attributes.

    Returns:
        Operation: A new ``Operation`` with the attributes replaced.

    Raises:
        ValidationError: The input type or paradigm name is invalid, or a restricted capability is being restored.
    """
    require_instance(operation, Operation, "annotate.operation")
    from oracq.infrastructure.validation import name as validate_name

    validate_name(paradigm)
    attrs = dict(operation.module.attributes)
    for cap in ("supports_adjoint", "supports_controlled"):
        if attrs.get(cap) is False and attributes.get(cap) is True:
            fail("INPUT_CAPABILITY", "annotate." + cap, False, True, "annotation cannot grant a capability that was already restricted")
    attrs.update(
        oracle_paradigm=paradigm,
        implementation_status="constructed",
        **attributes,
    )
    # An annotation never grants transformation capabilities that were already restricted.
    attrs.setdefault("supports_adjoint", True)
    attrs.setdefault("supports_controlled", True)
    if paradigm == "state_prep_isometry":
        attrs.setdefault("zero_input", True)
        attrs.setdefault("clean_work", True)
    return Operation(
        replace(operation.module, attributes=tuple(sorted(attrs.items()))), operation.dependencies
    )


def resources_for(*items: tuple[str, Operation]) -> dict[str, QRAM]:
    """Aggregate the resources of several operations and rename them as ``prefix__resource_name``.

    Args:
        *items: Pairs of the form ``(prefix, operation)``, where ``operation`` is
            an ``Operation`` carrying resource declarations.

    Returns:
        dict: Mapping from ``prefix__resource_name`` keys to resource types, for
        the nested ``Builder`` to declare resources.
    """
    return {f"{prefix}__{r.name}": r.type for prefix, op in items for r in op.module.resources}


def invoke(builder: Builder, operation: Operation, prefix: str = "", **arguments: Ref) -> None:
    """Invoke an operation in a builder with a prefixed resource mapping.

    Args:
        builder: The current ``Builder``, whose resources must already be declared under the same prefix.
        operation: The ``Operation`` to invoke.
        prefix: Resource name prefix; an empty string passes resource names through unchanged.
        **arguments: Bindings from the invoked module's register names to argument views.

    Raises:
        ValidationError: A register or resource argument name does not match the invoked module, or a dependency definition conflicts.
    """
    mapping = {
        r.name: f"{prefix}__{r.name}" if prefix else r.name for r in operation.module.resources
    }
    builder.call(operation, resources=mapping, **arguments)


@dataclass(frozen=True)
class XorDatabase(OracleView):
    """XOR database view: ``|address,data> -> |address,data XOR memory[address]>``.

    Reversible quantum access interface for a classical lookup table; the
    wrapped ``Operation`` contains exactly the two bits registers ``address``
    and ``data``, and the XOR semantics is self-inverse under any initial
    value.

    Attributes:
        operation: The wrapped ``Operation``.
    """

    oracle_kind = "database_xor"
    operation: Operation

    def xor_database(self) -> XorDatabase:
        """Return self; the view adapter method implementing ``XorDatabaseProtocol``.

        Returns:
            XorDatabase: This view itself.
        """
        return self

    def __post_init__(self) -> None:
        """Validate that the wrapped operation has exactly the two bits registers ``address`` and ``data``."""
        validate_signature(self.operation, ("address", "data"), "XorDatabase")
        regs = {r.name: r.type for r in self.operation.module.registers}
        if set(regs) != {"address", "data"} or any(r.kind != "bits" for r in regs.values()):
            raise ValidationError("XOR database requires an address/data bits interface")

    @property
    def address_width(self) -> int:
        """Bit width of the ``address`` register."""
        return next(r.type.width for r in self.operation.module.registers if r.name == "address")

    @property
    def data_width(self) -> int:
        """Bit width of the ``data`` register."""
        return next(r.type.width for r in self.operation.module.registers if r.name == "data")


@dataclass(frozen=True)
class StatePreparation(OracleView):
    """Zero-input state preparation view: an isometry acting from the zero-state subspace.

    The convention is ``V|0,0> = |psi,0>``, i.e. the ``work`` register is
    restored to zero; the adjoint and controlled operations require a
    reversible extension to exist. The wrapped ``Operation`` contains exactly
    the two bits registers ``target`` and ``work``.

    Attributes:
        operation: The wrapped ``Operation``.
    """

    oracle_kind = "state_prep_isometry"
    operation: Operation

    def state_preparation(self) -> StatePreparation:
        """Return self; the view adapter method implementing ``StatePreparationProtocol``.

        Returns:
            StatePreparation: This view itself.
        """
        return self

    def __post_init__(self) -> None:
        """Validate that the wrapped operation has exactly the two bits registers ``target`` and ``work``."""
        validate_signature(self.operation, ("target", "work"), "StatePreparation")
        regs = {r.name: r.type for r in self.operation.module.registers}
        if set(regs) != {"target", "work"} or any(r.kind != "bits" for r in regs.values()):
            raise ValidationError("StatePreparation requires a target/work bits interface")

    @classmethod
    def from_unitary(
        cls,
        operation: Operation,
        *,
        target: str | None = None,
        work: str | None = None,
        clean_work: bool = False,
    ) -> StatePreparation:
        """Explicitly assign the U|0> zero-input role; with a work register, the caller promises it is restored to zero.

        Args:
            operation: The full-unitary ``Operation``.
            target: Name of the register carrying the state preparation; by default all registers are concatenated into target.
            work: Work register name; may only be given when target is given explicitly.
            clean_work: Must be True when work has nonzero width, promising that work is restored to zero after preparation.

        Returns:
            StatePreparation: The preparation view wrapping the operation in the zero-input role.
        """
        require_instance(operation, Operation, "StatePreparation.from_unitary")
        if target is None:
            if work is not None:
                raise ValidationError("target must also be specified when work is given")
            operation.program()
            width = sum(r.type.width for r in operation.module.registers)
            b = Builder(
                _name("unitary_zero_input", operation),
                {"target": Bits(width), "work": Bits(0)},
                resources_for(("unitary", operation)),
            )
            cursor = 0
            arguments: dict[str, Ref] = {}
            for reg in operation.module.registers:
                arguments[reg.name] = b["target"][cursor : cursor + reg.type.width].reinterpret(
                    reg.type.kind
                )
                cursor += reg.type.width
            invoke(b, operation, "unitary", **arguments)
            return cls(
                annotate(b.finish(), "state_prep_isometry", zero_input=True, clean_work=True)
            )
        expected = (target,) if work is None else (target, work)
        if len(set(expected)) != len(expected):
            raise ValidationError("target and work cannot refer to the same register")
        validate_signature(operation, expected, "StatePreparation.from_unitary")
        widths = {r.name: r.type.width for r in operation.module.registers}
        ww = 0 if work is None else widths[work]
        if ww and clean_work is not True:
            fail(
                "INPUT_PROMISE",
                "StatePreparation.from_unitary.clean_work",
                True,
                clean_work,
                "must promise that work is restored to zero after applying U to the zero input",
            )
        b = Builder(
            _name("as_state_preparation", operation, target, work),
            {"target": Bits(widths[target]), "work": Bits(ww)},
            resources_for(("unitary", operation)),
        )
        arguments = {target: b["target"]}
        if work is not None:
            arguments[work] = b["work"]
        invoke(b, operation, "unitary", **arguments)
        return cls(annotate(b.finish(), "state_prep_isometry", zero_input=True, clean_work=True))

    @property
    def width(self) -> int:
        """Bit width of the ``target`` register."""
        return next(r.type.width for r in self.operation.module.registers if r.name == "target")

    @property
    def work_width(self) -> int:
        """Bit width of the ``work`` register."""
        return next(r.type.width for r in self.operation.module.registers if r.name == "work")


@dataclass(frozen=True)
class StateOracle(OracleView):
    """State output view: the success subspace is defined by ``signal`` being all zeros.

    The output carrier for algorithm stages such as solver kernels; the
    wrapped ``Operation`` contains exactly the two bits registers ``target``
    and ``signal``, and reading the ``signal == 0`` branch yields the target
    state.

    Attributes:
        operation: The wrapped ``Operation``.
    """

    oracle_kind = "state_oracle"
    operation: Operation

    def state_oracle(self) -> StateOracle:
        """Return self; the view adapter method implementing ``StateOracleProtocol``.

        Returns:
            StateOracle: This view itself.
        """
        return self

    def __post_init__(self) -> None:
        """Validate that the wrapped operation has exactly the ``target`` and ``signal`` registers."""
        validate_signature(self.operation, ("target", "signal"), "StateOracle")

    @property
    def width(self) -> int:
        """Bit width of the ``target`` register."""
        return next(r.type.width for r in self.operation.module.registers if r.name == "target")

    @property
    def signal_qubits(self) -> int:
        """Bit width of the ``signal`` register."""
        return next(r.type.width for r in self.operation.module.registers if r.name == "signal")


@dataclass(frozen=True)
class SparseAccess(OracleView):
    """Bundle view of CKS sparse access: a composition of the position and entry query operations, not disguised as a total unitary.

    Attributes:
        location: Position operation ``P_A(1)``, permuting ``index`` in place,
            with signature ``column``, ``index``, ``work``.
        entry: Entry operation ``P_A(2)``, queryable at any row and column,
            with signature ``row``, ``column``, ``data``.
        width: Number of bits of the matrix dimension, range 1..64.
        value_width: Word width of an entry value, range 1..64.
        sparsity: Upper bound on nonzeros per column, range 1..2**width.
    """

    location: Operation
    entry: Operation
    width: int
    value_width: int
    sparsity: int

    def sparse_access(self) -> SparseAccess:
        """Return self; the view adapter method implementing ``CKSSparseProtocol``.

        Returns:
            SparseAccess: This view itself.
        """
        return self

    def describe(self) -> OracleSpec:
        """Aggregate the two components into a read-only description of the whole bundle.

        The type is ``cks_sparse`` and ``anc_qubit`` is always ``None``; the
        adjoint and controlled capabilities are the conjunction of the two
        components, the bundle is marked ``closed`` only when every component
        is closed, and the ``sparsity`` and ``value_width`` parameters are
        attached.

        Returns:
            OracleSpec: Description of the whole bundle, with components named ``position`` and ``entry``.
        """
        from oracq.algorithms.input_model.contracts import (
            OracleCapabilities,
            OracleSpec,
            describe_oracle,
        )

        components = (
            ("position", describe_oracle(self.location)),
            ("entry", describe_oracle(self.entry)),
        )
        return OracleSpec(
            "cks_sparse",
            self.width,
            None,
            OracleCapabilities(
                all(s.capabilities.adjoint for _, s in components),
                all(s.capabilities.controlled for _, s in components),
            ),
            "closed" if all(s.implementation == "closed" for _, s in components) else "open",
            parameters=(("sparsity", self.sparsity), ("value_width", self.value_width)),
            components=components,
        )

    def __post_init__(self) -> None:
        """Validate the bit widths, the sparsity value, and the register signatures of the two component operations."""
        positive_integer(self.width, "SparseAccess.width", maximum=64)
        positive_integer(self.value_width, "SparseAccess.value_width", maximum=64)
        positive_integer(self.sparsity, "SparseAccess.sparsity", maximum=1 << self.width)
        validate_signature(self.location, ("column", "index", "work"), "SparseAccess.position")
        validate_signature(self.entry, ("row", "column", "data"), "SparseAccess.entry")
        if (
            not 1 <= self.width <= 64
            or not 1 <= self.value_width <= 64
            or not 1 <= self.sparsity <= 1 << self.width
        ):
            raise ValidationError("SparseAccess bit width or sparsity is invalid")
        location = {r.name: r.type for r in self.location.module.registers}
        entry = {r.name: r.type for r in self.entry.module.registers}
        if (
            set(location) != {"column", "index", "work"}
            or location["column"] != Bits(self.width)
            or location["index"] != Bits(self.width)
        ):
            raise ValidationError("SparseAccess requires an in-place column/index/work position interface")
        if entry != {
            "row": Bits(self.width),
            "column": Bits(self.width),
            "data": Bits(self.value_width),
        }:
            raise ValidationError("SparseAccess requires an entry XOR interface for any row and column")


def abstract_database(name: str, address_width: int, data_width: int) -> XorDatabase:
    """Return an open XOR database declaration with no bound implementation.

    Args:
        name: Declaration name.
        address_width: Bit width of the ``address`` register.
        data_width: Bit width of the ``data`` register.

    Returns:
        XorDatabase: An open declaration with an empty body, closed by binding an implementation.
    """
    return XorDatabase(
        declare(
            name,
            {"address": Bits(address_width), "data": Bits(data_width)},
            paradigm="database_xor",
        )
    )


def abstract_state_prep(
    name: str,
    width: int,
    work_width: int = 0,
    *,
    reversible: bool = True,
) -> StatePreparation:
    """Return an open state preparation declaration with no bound implementation.

    Args:
        name: Declaration name.
        width: Bit width of the ``target`` register.
        work_width: Bit width of the ``work`` register, zero by default.
        reversible: When ``False``, declare no support for adjoint or controlled.

    Returns:
        StatePreparation: An open declaration carrying the ``zero_input`` and ``clean_work`` promises.
    """
    return StatePreparation(
        declare(
            name,
            {"target": Bits(width), "work": Bits(work_width)},
            paradigm="state_prep_isometry",
            attributes={"zero_input": True, "clean_work": True},
            supports_adjoint=reversible,
            supports_controlled=reversible,
        )
    )


def abstract_block_encoding(
    name: str,
    width: int,
    signal_width: int,
    alpha: float,
) -> BlockEncoding:
    """Return an open block encoding declaration with no bound implementation.

    Args:
        name: Declaration name.
        width: Bit width of the ``target`` register.
        signal_width: Bit width of the ``signal`` register.
        alpha: Normalization constant, written into the ``be_alpha`` attribute.

    Returns:
        BlockEncoding: An open declaration with an empty body.
    """
    return BlockEncoding(
        declare(
            name,
            {"target": Bits(width), "signal": Bits(signal_width)},
            paradigm="block_encoding",
            attributes={"be_alpha": alpha},
        )
    )


def abstract_sparse_access(
    name: str,
    width: int,
    value_width: int,
    sparsity: int,
    work_width: int | None = None,
) -> SparseAccess:
    """Return an open sparse-access bundle with no bound implementation.

    Generates the two open operations ``name_position`` and ``name_entry``; the
    position operation carries the ``sparsity`` and full permutation extension
    attributes.

    Args:
        name: Declaration name prefix.
        width: Number of bits of the matrix dimension.
        value_width: Word width of an entry value.
        sparsity: Upper bound on nonzeros per column, at most ``2**width``.
        work_width: ``work`` bit width of the position operation, defaults to ``width``.

    Returns:
        SparseAccess: A bundle of the two open declarations for position and entry.

    Raises:
        ValidationError: ``sparsity`` exceeds the dimension.
    """
    if not 1 <= sparsity <= 1 << width:
        raise ValidationError("sparsity exceeds the dimension")
    work_width = width if work_width is None else work_width
    location = declare(
        name + "_position",
        {"column": Bits(width), "index": Bits(width), "work": Bits(work_width)},
        paradigm="sparse_location_inplace",
        attributes={"sparsity": sparsity, "full_permutation_extension": True},
    )
    entry = declare(
        name + "_entry",
        {"row": Bits(width), "column": Bits(width), "data": Bits(value_width)},
        paradigm="sparse_entry_xor",
    )
    return SparseAccess(location, entry, width, value_width, sparsity)


def gate_database(
    address_width: int,
    data_width: int,
    table: Mapping[int, int] | Sequence[int],
    *,
    name: str | None = None,
) -> XorDatabase:
    """Implement an XOR database with a gate-level truth table.

    Nonzero table words apply X gates bit by bit under address control; the
    gate count grows with the table size, suitable as a small-instance
    witness.

    Args:
        address_width: Bit width of the ``address`` register.
        data_width: Bit width of the ``data`` register.
        table: Mapping from addresses to words, or a word sequence enumerated by address.
        name: Module name; generated deterministically from the content by default.

    Returns:
        XorDatabase: Gate-level implementation annotated ``implementation="gate_truth_table"``.

    Raises:
        ValidationError: A table entry is not an integer, or an address or word exceeds its bit-width range.
    """
    items = tuple(sorted(table.items())) if isinstance(table, dict) else tuple(enumerate(table))
    b = Builder(
        name or _name("xor_table", address_width, data_width, items),
        {"address": Bits(address_width), "data": Bits(data_width)},
    )
    for address, value in items:
        if (
            type(address) is not int
            or type(value) is not int
            or not (0 <= address < 1 << address_width and 0 <= value < 1 << data_width)
        ):
            raise ValidationError("lookup table address or word is out of range")
        if value:
            with b.control(b["address"], address):
                for bit in range(data_width):
                    if (value >> bit) & 1:
                        b.x(b["data"][bit])
    return XorDatabase(annotate(b.finish(), "database_xor", implementation="gate_truth_table"))


def qram_database(address_width: int, data_width: int, *, name: str | None = None) -> XorDatabase:
    """Implement an XOR database with a QRAM resource.

    Declares a ``QRAM(address_width, data_width)`` resource named ``table`` and
    emits one query primitive; the data table is supplied at execution time as
    a memory and never enters the IR.

    Args:
        address_width: Bit width of the ``address`` register.
        data_width: Bit width of the ``data`` register.
        name: Module name; generated from the bit widths by default.

    Returns:
        XorDatabase: Resource-based implementation annotated ``implementation="qram"``.
    """
    b = Builder(
        name or f"qram_xor_{address_width}_{data_width}",
        {"address": Bits(address_width), "data": Bits(data_width)},
        {"table": QRAM(address_width, data_width)},
    )
    b.qram("table", b["address"], b["data"])
    return XorDatabase(annotate(b.finish(), "database_xor", implementation="qram"))


def basis_state(width: int, value: int = 0, *, work_width: int = 0) -> StatePreparation:
    """Prepare the computational basis state ``|value>`` with bit-by-bit X gates.

    Args:
        width: Bit width of the ``target`` register.
        value: Integer index of the target basis state, range ``0 .. 2**width - 1``.
        work_width: Bit width of the ``work`` register, zero by default.

    Returns:
        StatePreparation: Gate-level preparation carrying the ``zero_input`` promise.

    Raises:
        ValidationError: ``value`` is not an integer or is out of range.
    """
    if type(value) is not int or not 0 <= value < 1 << width:
        raise ValidationError("basis state value is out of range")
    b = Builder(
        _name("basis", width, value, work_width), {"target": Bits(width), "work": Bits(work_width)}
    )
    for bit in range(width):
        if (value >> bit) & 1:
            b.x(b["target"][bit])
    return StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True))


def uniform_state(width: int, *, work_width: int = 0) -> StatePreparation:
    """Apply Hadamard to the whole ``target`` to prepare the uniform superposition state.

    Args:
        width: Bit width of the ``target`` register.
        work_width: Bit width of the ``work`` register, zero by default.

    Returns:
        StatePreparation: Gate-level preparation carrying the ``zero_input`` promise.
    """
    b = Builder(
        _name("uniform", width, work_width), {"target": Bits(width), "work": Bits(work_width)}
    )
    b.h(b["target"])
    return StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True))


def _state_angles(
    amplitudes: Iterable[complex],
) -> tuple[tuple[complex, ...], int, list[tuple[int, int, int, float]]]:
    """Validate the amplitude vector and derive the rotation angle of each node in the Ry rotation tree."""
    values = tuple(complex(v) for v in amplitudes)
    if not values or len(values) & (len(values) - 1):
        raise ValidationError("state vector length must be a power of two")
    if not all(math.isfinite(v.real) and math.isfinite(v.imag) for v in values):
        raise ValidationError("state amplitudes must be finite")
    if sum(abs(v) ** 2 for v in values) == 0:
        raise ValidationError("the zero vector cannot be prepared as a normalized state")
    n = (len(values) - 1).bit_length()
    nodes: list[tuple[int, int, int, float]] = []
    for depth in range(n):
        bit = n - depth - 1
        for prefix in range(1 << depth):
            start, size = prefix << (bit + 1), 1 << bit
            left = sum(abs(v) ** 2 for v in values[start : start + size])
            right = sum(abs(v) ** 2 for v in values[start + size : start + 2 * size])
            nodes.append((depth, prefix, bit, 2 * math.atan2(math.sqrt(right), math.sqrt(left))))
    return values, n, nodes


def gate_state_prep(
    amplitudes: Iterable[complex],
    *,
    work_width: int = 0,
    name: str | None = None,
) -> StatePreparation:
    """Prepare an arbitrary complex-amplitude state with a multiplexed Ry rotation tree.

    The amplitudes decompose over a binary tree into controlled ``Ry`` angles,
    and nonzero phase components are compensated pointwise by a controlled
    ``global_phase``, so amplitudes may be complex.

    Args:
        amplitudes: Finite complex amplitude sequence of power-of-two length with nonzero norm.
        work_width: Bit width of the ``work`` register, zero by default.
        name: Module name; generated deterministically from the amplitude content by default.

    Returns:
        StatePreparation: Gate-level preparation annotated
        ``implementation="multiplexed_rotations"``.

    Raises:
        ValidationError: The amplitude vector is invalid, or there is a single amplitude and no target bit.
    """
    values, n, nodes = _state_angles(amplitudes)
    if n == 0:
        raise ValidationError("this state preparation needs at least one target bit")
    b = Builder(
        name or _name("state", values, work_width), {"target": Bits(n), "work": Bits(work_width)}
    )
    for depth, prefix, bit, angle in nodes:
        if not depth:
            b.ry(b["target"][bit], angle)
        else:
            with b.control(b["target"][bit + 1 :], prefix):
                b.ry(b["target"][bit], angle)
    for index, value in enumerate(values):
        if value and cmath.phase(value):
            with b.control(b["target"], index):
                b.global_phase(cmath.phase(value))
    return StatePreparation(
        annotate(
            b.finish(),
            "state_prep_isometry",
            zero_input=True,
            implementation="multiplexed_rotations",
        )
    )


def qram_state_prep(width: int, angle_width: int = 8) -> StatePreparation:
    """Implement rotation-tree state preparation with a QRAM angle table.

    Layer by layer, the high-bit prefix of ``target`` is swapped into the
    address and the angle table is queried; controlled Ry rotations are
    applied according to the bits of the angle word and the address is then
    restored to zero, with two queries per layer; the angle resolution
    ``2*pi/2**angle_width`` is the source of the quantization error.

    Args:
        width: Bit width of the ``target`` register.
        angle_width: Bit width of the angle word.

    Returns:
        StatePreparation: Resource-based preparation annotated
        ``implementation="qram_rotation_tree"``; ``work`` is
        ``max(1, width) + angle_width`` bits, with ``qram_queries`` attached
        as an attribute.
    """
    address_width = max(1, width)
    b = Builder(
        f"qram_state_{width}_{angle_width}",
        {"target": Bits(width), "work": Bits(address_width + angle_width)},
        {"angles": QRAM(address_width, angle_width)},
    )
    addr, angle = b["work"][:address_width], b["work"][address_width:]
    for depth in range(width):
        bit = width - depth - 1
        if depth:
            b.xor(b["target"][bit + 1 :], addr[:depth])
        offset = (1 << depth) - 1
        b.add_const(addr.reinterpret("uint"), offset)
        b.qram("angles", addr, angle)
        for k in range(angle_width):
            with b.control(angle[k]):
                b.ry(b["target"][bit], 2 * math.pi * (1 << k) / (1 << angle_width))
        b.qram("angles", addr, angle)
        b.add_const(addr.reinterpret("uint"), (-offset) % (1 << address_width))
        if depth:
            b.xor(b["target"][bit + 1 :], addr[:depth])
    return StatePreparation(
        annotate(
            b.finish(),
            "state_prep_isometry",
            zero_input=True,
            implementation="qram_rotation_tree",
            clean_work=True,
            qram_queries=2 * width,
        )
    )


def qram_state_angles(amplitudes: Iterable[complex], angle_width: int = 8) -> dict[int, int]:
    """Generate the QRAM rotation-tree angle table from an amplitude vector.

    Keys are tree node indices of the form ``(1 << depth) - 1 + prefix`` and
    values are Ry angle words quantized to the ``2*pi/2**angle_width`` grid,
    supplied at execution time as the angle-table memory.

    Args:
        amplitudes: Nonnegative real amplitude sequence of power-of-two length with nonzero norm.
        angle_width: Bit width of the angle word.

    Returns:
        dict: Mapping from tree node indices to angle words.

    Raises:
        ValidationError: The amplitude vector is invalid, or contains amplitudes with negative real part or nonzero imaginary part.
    """
    values, _, nodes = _state_angles(amplitudes)
    if any(v.imag or v.real < 0 for v in values):
        raise ValidationError("the current QRAM angle-table implementation accepts only nonnegative real amplitudes")
    scale = (1 << angle_width) / (2 * math.pi)
    return {
        ((1 << depth) - 1 + prefix): round(angle * scale) % (1 << angle_width)
        for depth, prefix, _, angle in nodes
    }


def phase_marks(width: int, marked: Iterable[int]) -> Operation:
    """Generate a phase oracle applying phase ``pi`` to the given computational basis states.

    Args:
        width: Bit width of the ``target`` register.
        marked: Iterable of integer indices of basis states to mark; duplicates are removed automatically.

    Returns:
        Operation: Operation annotated with the ``phase_oracle`` paradigm.
    """
    marked = tuple(sorted(set(marked)))
    b = Builder(_name("phase_marks", width, marked), {"target": Bits(width)})
    for value in marked:
        with b.control(b["target"], value):
            b.global_phase(math.pi)
    return annotate(b.finish(), "phase_oracle")


def _transposition(b: Builder, ref: Ref, first: int, second: int) -> None:
    """Compile the transposition of ``first`` and ``second`` into a round-trip gate sequence of controlled X."""
    path = [first]
    for bit in range(ref.width):
        if ((first ^ second) >> bit) & 1:
            path.append(path[-1] ^ (1 << bit))
    edges = list(zip(path, path[1:], strict=False))
    for left, right in edges + list(reversed(edges[:-1])):
        bit = (left ^ right).bit_length() - 1
        controls = fuse(ref[:bit], ref[bit + 1 :])
        value = (left & ((1 << bit) - 1)) | ((left >> (bit + 1)) << bit)
        if controls.width:
            with b.control(controls, value):
                b.x(ref[bit])
        else:
            b.x(ref[bit])


def sparse_location_gate(
    width: int,
    permutations: Iterable[Iterable[int]],
    *,
    work_width: int | None = None,
    name: str | None = None,
) -> Operation:
    """Implement the CKS in-place position query with a gate-level permutation network.

    Each column decomposes its permutation into transpositions by cycle under
    column control, and each transposition becomes a round-trip path of
    controlled X; the gate count grows with ``2**width``, suitable only as a
    small-instance witness.

    Args:
        width: Bit width of the ``column`` and ``index`` registers.
        permutations: ``2**width`` permutations given per column, each a complete
            arrangement of ``0 .. 2**width - 1``.
        work_width: Bit width of the ``work`` register, defaults to ``width``.
        name: Module name; generated deterministically from the content by default.

    Returns:
        Operation: Operation annotated with the ``sparse_location_inplace`` paradigm and the full permutation extension attribute.

    Raises:
        ValidationError: The number of columns is not ``2**width``, or some column is not a complete permutation.
    """
    work_width = width if work_width is None else work_width
    permutations = tuple(tuple(row) for row in permutations)
    dimension = 1 << width
    if len(permutations) != dimension or any(
        sorted(row) != list(range(dimension)) for row in permutations
    ):
        raise ValidationError("the CKS in-place position implementation requires a full permutation extension for every column")
    b = Builder(
        name or _name("sparse_position", width, permutations, work_width),
        {"column": Bits(width), "index": Bits(width), "work": Bits(work_width)},
    )
    for column, permutation in enumerate(permutations):
        visited: set[int] = set()
        with b.control(b["column"], column):
            for first in range(dimension):
                if first in visited:
                    continue
                cycle: list[int] = []
                current = first
                while current not in visited:
                    visited.add(current)
                    cycle.append(current)
                    current = cast("tuple[int, ...]", permutation)[current]
                for second in cycle[1:]:
                    _transposition(b, b["index"], first, second)
    return annotate(b.finish(), "sparse_location_inplace", full_permutation_extension=True)


def sparse_location_qram(width: int) -> Operation:
    """Implement the CKS in-place position query with a forward and an inverse QRAM table.

    Done in three steps: ``work ^= forward[column, index]``,
    ``swap(index, work)``, ``work ^= inverse[column, index]``; the last step
    restores ``work`` to zero while ``index`` holds the new value.

    Args:
        width: Bit width of the ``column``, ``index``, and ``work`` registers.

    Returns:
        Operation: Operation annotated with the ``sparse_location_inplace``
    paradigm, declaring the two ``QRAM(2*width, width)`` resources ``forward``
    and ``inverse``.
    """
    b = Builder(
        f"sparse_position_qram_{width}",
        {"column": Bits(width), "index": Bits(width), "work": Bits(width)},
        {"forward": QRAM(2 * width, width), "inverse": QRAM(2 * width, width)},
    )
    b.qram("forward", fuse(b["column"], b["index"]), b["work"])
    b.swap(b["index"], b["work"])
    b.qram("inverse", fuse(b["column"], b["index"]), b["work"])
    return annotate(b.finish(), "sparse_location_inplace", full_permutation_extension=True)


def sparse_entry(database: XorDatabase, width: int) -> Operation:
    """Adapt an XOR database into a matrix entry query at any row and column.

    ``row`` and ``column`` concatenate into a ``2*width``-bit address driving
    the database, and ``data`` holds the XOR result.

    Args:
        database: An XOR database whose ``address_width`` is exactly ``2*width``.
        width: Bit width of rows and columns respectively.

    Returns:
        Operation: Entry query operation annotated with the ``sparse_entry_xor`` paradigm.

    Raises:
        ValidationError: The database address width is not ``2*width``.
    """
    if database.address_width != 2 * width:
        raise ValidationError("matrix entry query requires both row and column address bits")
    b = Builder(
        _name("matrix_entry", database.operation),
        {"row": Bits(width), "column": Bits(width), "data": Bits(database.data_width)},
        resources_for(("db", database.operation)),
    )
    invoke(b, database.operation, "db", address=fuse(b["row"], b["column"]), data=b["data"])
    return annotate(b.finish(), "sparse_entry_xor")


def diagonal_block_encoding(
    database: XorDatabase,
    *,
    alpha: float = 1.0,
    angle_scale: float | None = None,
) -> BlockEncoding:
    """Control a signal Ry by the queried word, encoding the diagonal matrix determined by that angle table.

    Args:
        database: XOR database storing rotation angle words, with address driving target.
        alpha: Block encoding normalization constant.
        angle_scale: Angle step of bit k of the data word; defaults to ``2π/2**data_width``.

    Returns:
        BlockEncoding: Block encoding whose diagonal entries are ``alpha*cos(encoded_angle/2)``.
    """
    angle_scale = 2 * math.pi / (1 << database.data_width) if angle_scale is None else angle_scale
    n, w = database.address_width, database.data_width
    b = Builder(
        _name("diagonal_be", database.operation, alpha, angle_scale),
        {"target": Bits(n), "signal": Bits(w + 1)},
        resources_for(("db", database.operation)),
    )
    word, flag = b["signal"][:w], b["signal"][w:]
    invoke(b, database.operation, "db", address=b["target"], data=word)
    for k in range(w):
        with b.control(word[k]):
            b.ry(flag, angle_scale * (1 << k))
    with b.adjoint():
        invoke(b, database.operation, "db", address=b["target"], data=word)
    return BlockEncoding(
        annotate(
            b.finish(),
            "block_encoding",
            be_alpha=alpha,
            matrix_interpretation="alpha*cos(encoded_angle/2) on diagonal",
        )
    )


def banked_database(
    address_width: int,
    data_width: int,
    *,
    name: str = "BankedData",
    abstract: bool = False,
) -> Operation:
    """Split a wide data word into registers and QRAM banks of at most 64 bits each.

    Args:
        address_width: Bit width of the address word.
        data_width: Total bit width of the logical data word, which must be a positive integer.
        name: Name of the generated operation, ``BankedData`` by default.
        abstract: When True, generate an open declaration with an empty body instead of a QRAM implementation.

    Returns:
        Operation: Database operation split little-endian into banks of at most 64 bits each.
    """
    if type(data_width) is not int or data_width < 1:
        raise ValidationError("banked data width must be a positive integer")
    chunks = tuple(min(64, data_width - offset) for offset in range(0, data_width, 64))
    registers = {
        "address": Bits(address_width),
        **{f"data{i}": Bits(w) for i, w in enumerate(chunks)},
    }
    if abstract:
        return declare(
            name,
            registers,
            paradigm="database_xor",
            attributes={"logical_data_width": data_width, "word_order": "little_endian_banks"},
        )
    b = Builder(
        _name("banked_qram", address_width, data_width),
        registers,
        {f"bank{i}": QRAM(address_width, w) for i, w in enumerate(chunks)},
    )
    for i, _ in enumerate(chunks):
        b.qram(f"bank{i}", b["address"], b[f"data{i}"])
    return annotate(
        b.finish(), "database_xor", logical_data_width=data_width, word_order="little_endian_banks"
    )
