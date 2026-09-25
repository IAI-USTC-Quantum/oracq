"""Toffoli+Clifford+T+QRAM level resource estimation.

Counting descends the OriginIR-ext → basis lowering chain following the real
emission recipes (``basis.mcx`` and ``basis.controlled_u3``): the module graph
aggregates compositionally, ``Repeat`` multiplies symbolic counts, and no
exponential expansion happens. Rotation gates are classified by angle: integer
multiples of π/4 fall into exact Clifford+T (T/SDG/S/Z counts), the rest are
recorded as rotation atoms pending synthesis, with the T cost estimated by the
Ross–Selinger style model ``synthesis_t_per_rotation`` (default
``ceil(3*log2(1/eps))``). Global phase is not counted (consistent with strict
export).

Control costs align with the exporter line by line: zero-valued control bits
are flipped once on each side of every emitted line group (primitive gate /
call / Repeat call lines), for 2 X gates per zero bit; controlled single-qubit
gates follow the ``controlled_u3`` recipe (a 2(c−1) Toffoli ladder + an
intermediate network when c≥2), and X-type multi-controls follow the ``mcx``
recipe (2c−3 Toffoli when c≥3).

QRAM query counts are tallied per resource over ``Load`` nodes (each Load = 1
query) as a first-class metric independent of gate-level cost. QRAM random
writes are tallied per resource over ``Store`` nodes into ``qram_writes``:
storage cells are modeled as classical cells, and random writes do not enter
gate cost.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import asdict, dataclass, field, replace
from typing import cast, overload

from oracq.infrastructure.ir import (
    Adjoint,
    Call,
    Control,
    Instruction,
    Load,
    Primitive,
    Program,
    Repeat,
    Store,
    ValidationError,
)
from oracq.infrastructure.layout import workspace_table
from oracq.infrastructure.linking import unresolved
from oracq.infrastructure.validation import validate

_PI = math.pi

CLIFFORD_ATOMS = ("h", "x", "y", "z", "s", "sdg", "cnot", "cz")
"""Set of Clifford atom names; ``ResourceEstimate.clifford`` summarizes ``atoms`` counts over these names."""


def _rz_exact_counts(k: int) -> Counter[str]:
    """Exact Clifford+T counts of rz(k·π/4), with k normalized mod 8."""
    k %= 8
    if k > 4:
        counts = _rz_exact_counts(8 - k)
        result: Counter[str] = Counter()
        for atom, n in counts.items():
            result[{"t": "tdg", "s": "sdg", "tdg": "t", "sdg": "s"}.get(atom, atom)] = n
        return result
    return Counter(
        {
            0: {},
            1: {"t": 1},
            2: {"s": 1},
            3: {"s": 1, "t": 1},
            4: {"z": 1},
        }[k]
    )


def _quarter_turns(angle: float) -> int | None:
    """Whether angle is an integer multiple of π/4; returns k mod 8 if so, otherwise None."""
    k = round(angle / (_PI / 4))
    if abs(angle - k * _PI / 4) < 1e-12:
        return k % 8
    return None


def classify_rz(angle: float) -> tuple[Counter[str], list[tuple[str, float]]]:
    """Classify an rz/phase atom: exact Clifford+T or a rotation pending synthesis.

    Args:
        angle: rotation angle in radians.

    Returns:
        tuple[Counter[str], list[tuple[str, float]]]: pair ``(exact atom counts,
        rotations pending synthesis)``; an angle that is an integer multiple of
        π/4 goes into the former, otherwise the whole entry is recorded as the
        rotation ``("rz", angle)`` pending synthesis.
    """
    k = _quarter_turns(angle)
    if k is None:
        return Counter(), [("rz", angle)]
    return _rz_exact_counts(k), []


def classify_ry(angle: float) -> tuple[Counter[str], list[tuple[str, float]]]:
    """Classify an ry atom: RY(θ) = S·H·RZ(θ)·H·S†.

    Args:
        angle: rotation angle in radians.

    Returns:
        tuple[Counter[str], list[tuple[str, float]]]: pair ``(exact atom counts,
        rotations pending synthesis)``; an angle that is not an integer multiple
        of π/4 is recorded whole as the rotation ``("ry", angle)`` pending
        synthesis, and a zero angle returns two empty containers.
    """
    k = _quarter_turns(angle)
    if k is None:
        return Counter(), [("ry", angle)]
    if k == 0:
        return Counter(), []
    counts = Counter({"s": 1, "h": 2, "sdg": 1})
    counts += _rz_exact_counts(k)
    return counts, []


def classify_u3(
    theta: float, phi: float, lam: float
) -> tuple[Counter[str], list[tuple[str, float]]]:
    """U3(θ,φ,λ) = RZ(φ)·RY(θ)·RZ(λ), dropping global phase, consistent with strict export.

    Common gates H/X/Y are first matched as exact single atoms, then fall back
    to the Euler decomposition.

    Args:
        theta: Y-axis Euler angle in radians.
        phi: left Z-axis Euler angle in radians.
        lam: right Z-axis Euler angle in radians.

    Returns:
        tuple[Counter[str], list[tuple[str, float]]]: pair ``(exact atom counts,
        rotations pending synthesis)``; the H/X/Y special cases match a single
        atom, the rest are summarized item by item via the Euler decomposition.
    """

    def close(a: float, b: float) -> bool:
        """Decide whether two floats are equal within tolerance."""
        return abs(a - b) < 1e-12

    if close(theta, _PI / 2) and close(phi, 0.0) and close(lam, _PI):
        return Counter({"h": 1}), []
    if close(theta, _PI) and close(phi, 0.0) and close(lam, _PI):
        return Counter({"x": 1}), []
    if close(theta, _PI) and close(phi, _PI / 2) and close(lam, _PI / 2):
        return Counter({"y": 1}), []
    counts: Counter[str]
    rotations: list[tuple[str, float]]
    counts, rotations = Counter(), []
    for counts_part, rotations_part in (
        classify_rz(lam),
        classify_ry(theta),
        classify_rz(phi),
    ):
        counts += counts_part
        rotations += rotations_part
    return counts, rotations


def mcx_counts(n_controls: int) -> Counter[str]:
    """basis.mcx recipe, the X-type multi-control.

    Args:
        n_controls: number of effective control bits, a nonnegative integer.

    Returns:
        Counter[str]: atom counts of the multi-control X; c of 0/1/2 gives a
        single X, H-CZ-H and a single Toffoli respectively, and c≥3 gives
        2c−3 Toffoli gates.
    """
    if n_controls <= 0:
        return Counter({"x": 1})
    if n_controls == 1:
        return Counter({"h": 2, "cz": 1})
    if n_controls == 2:
        return Counter({"toffoli": 1})
    return Counter({"toffoli": 2 * n_controls - 3})


def _controlled_u3_counts(
    n_controls: int, theta: float, phi: float, lam: float, global_angle: float = 0.0
) -> tuple[Counter[str], list[tuple[str, float]]]:
    """Counts of the basis.controlled_u3 recipe; includes a 2(c−1) Toffoli ladder when c≥2."""
    counts: Counter[str]
    rotations: list[tuple[str, float]]
    counts, rotations = Counter(), []
    if n_controls == 0:
        part, rot = classify_u3(theta, phi, lam)
        counts += part
        rotations += rot
        if global_angle:
            for _ in range(2):
                part, rot = classify_rz(global_angle)
                counts += part
                rotations += rot
            counts["x"] += 2
        return counts, rotations
    if n_controls >= 2:
        counts["toffoli"] += 2 * (n_controls - 1)
    counts += Counter({"h": 4, "cz": 2})
    for angle in ((lam + phi) / 2 + global_angle, (lam - phi) / 2):
        part, rot = classify_rz(angle)
        counts += part
        rotations += rot
    for theta_part, phi_part, lam_part in (
        (-theta / 2, 0.0, -(phi + lam) / 2),
        (theta / 2, phi, 0.0),
    ):
        part, rot = classify_u3(theta_part, phi_part, lam_part)
        counts += part
        rotations += rot
    return counts, rotations


# Table of (θ, φ, λ, global phase) parameters for single-qubit gates in basis.py
_GATE_PARAMETERS = {
    "h": (_PI / 2, 0.0, _PI, 0.0),
    "y": (_PI, _PI / 2, _PI / 2, 0.0),
    "z": (0.0, 0.0, _PI, 0.0),
    "s": (0.0, 0.0, _PI / 2, 0.0),
    "t": (0.0, 0.0, _PI / 4, 0.0),
}


def _gate_counts(
    op: str, angle: float | None, n_controls: int
) -> tuple[Counter[str], list[tuple[str, float]]]:
    """Counts of a single 1q gate under n_controls effective controls."""
    if op == "x":
        return mcx_counts(n_controls), []
    if op in _GATE_PARAMETERS:
        theta, phi, lam, glob = _GATE_PARAMETERS[op]
        return _controlled_u3_counts(n_controls, theta, phi, lam, glob)
    if op == "phase":
        return _controlled_u3_counts(n_controls, 0.0, 0.0, cast(float, angle), 0.0)
    if op == "ry":
        return _controlled_u3_counts(n_controls, cast(float, angle), 0.0, 0.0, 0.0)
    if op == "rx":
        return _controlled_u3_counts(n_controls, cast(float, angle), -_PI / 2, _PI / 2, 0.0)
    if op == "rz":
        return _controlled_u3_counts(n_controls, 0.0, 0.0, cast(float, angle), -cast(float, angle) / 2)
    raise ValidationError(f"primitive gate not supported by resource estimation: {op}")


class RotationCounts(Sequence[tuple[str, float]]):
    """Store rotation multiplicities by axis and angle; supports read-only length, indexing and lazy iteration."""

    def __init__(self, values: Iterable[tuple[str, float]] = ()) -> None:
        self.counts: Counter[tuple[str, float]] = values.counts.copy() if isinstance(values, RotationCounts) else Counter(values)

    @property
    def total(self) -> int:
        """Total number of rotations; not limited by the machine integer bound on Python sequence lengths."""
        return sum(self.counts.values())

    def __len__(self) -> int:
        return self.total

    def __iter__(self) -> Iterator[tuple[str, float]]:
        for value, count in self.counts.items():
            for _ in range(count):
                yield value

    @overload
    def __getitem__(self, index: int) -> tuple[str, float]: ...

    @overload
    def __getitem__(self, index: slice) -> list[tuple[str, float]]: ...

    def __getitem__(self, index: int | slice) -> tuple[str, float] | list[tuple[str, float]]:
        if isinstance(index, slice):
            indices = range(*index.indices(self.total))
            if len(indices) > 100_000:
                raise ValidationError("rotation slice too large; read the compact multiplicities from counts")
            return [self[i] for i in indices]
        if index < 0:
            index += self.total
        if index >= 0:
            for value, count in self.counts.items():
                if index < count:
                    return value
                index -= count
        raise IndexError("rotation index out of range")

    def add_scaled(self, other: RotationCounts, factor: int = 1) -> None:
        """Combine multiplicities without copying the rotation list."""
        for value, count in other.counts.items():
            self.counts[value] += count * factor


@dataclass(frozen=True, order=True)
class OracleCall:
    """Call form of an open oracle; resource names are relative to the current module, with argument substitution already done in the entry report."""

    module: str
    controls: int = 0
    adjoint: bool = False
    resources: tuple[str, ...] = ()


@dataclass
class ResourceEstimate:
    "Toffoli+Clifford+T+QRAM level resource ledger."

    qubits: int | None
    atoms: Counter = field(default_factory=Counter)
    rotations: RotationCounts = field(default_factory=RotationCounts)
    qram_queries: Counter = field(default_factory=Counter)
    qram_writes: Counter = field(default_factory=Counter)
    mcx_ancilla: int = 0
    oracle_calls: Counter[OracleCall] = field(default_factory=Counter)
    unknown_workspace: tuple[str, ...] = ()
    qubits_lower_bound: int = 0

    @property
    def complete(self) -> bool:
        """Whether known implementations cover all entry-reachable slots."""
        return not self.unknown_workspace

    @property
    def toffoli(self) -> int:
        """Total number of Toffoli gates in ``atoms``."""
        return self.atoms.get("toffoli", 0)

    @property
    def t_exact(self) -> int:
        """Total number of T and TDG gates falling exactly into Clifford+T, excluding the T cost of rotations pending synthesis."""
        return self.atoms.get("t", 0) + self.atoms.get("tdg", 0)

    @property
    def clifford(self) -> int:
        """Sum of counts of the atoms listed in ``CLIFFORD_ATOMS`` within ``atoms``."""
        return sum(self.atoms.get(name, 0) for name in CLIFFORD_ATOMS)

    @property
    def qram_total(self) -> int:
        """Sum of query counts over all resources in ``qram_queries``."""
        return sum(self.qram_queries.values())

    @property
    def qram_write_total(self) -> int:
        """Sum of random-write counts over all resources in ``qram_writes``."""
        return sum(self.qram_writes.values())

    @property
    def gate_total(self) -> int:
        """Sum of all atom counts in ``atoms``, excluding QRAM queries and writes."""
        return sum(self.atoms.values())

    def synthesis_t_per_rotation(self, epsilon: float = 1e-10) -> int:
        """Ross–Selinger style leading-term model; replaceable as a whole.

        Args:
            epsilon: synthesis precision of a single rotation, in (0, 1).

        Returns:
            int: T gate cost per rotation pending synthesis, i.e. ``ceil(3*log2(1/epsilon))``, at least 1.
        """
        return max(1, math.ceil(3 * math.log2(1 / epsilon)))

    def t_synthesis(self, epsilon: float = 1e-10) -> int:
        """Total T cost of rotations pending synthesis: rotation count times ``synthesis_t_per_rotation(epsilon)``.

        Args:
            epsilon: synthesis precision of a single rotation.

        Returns:
            int: estimated T gate cost of all rotations pending synthesis.
        """
        return self.rotations.total * self.synthesis_t_per_rotation(epsilon)

    def t_total(self, epsilon: float = 1e-10) -> int:
        """Sum of the exact T count and the T cost of rotations pending synthesis, i.e. ``t_exact + t_synthesis``.

        Args:
            epsilon: synthesis precision of a single rotation.

        Returns:
            int: total T gate count: exact T/TDG counts plus the synthesis T cost.
        """
        return self.t_exact + self.t_synthesis(epsilon)

    def to_dict(self, epsilon: float = 1e-10) -> dict[str, object]:
        """Export a JSON-friendly flat resource ledger.

        Args:
            epsilon: rotation synthesis precision, also recorded in the result under the ``epsilon`` key.

        Returns:
            dict: includes the qubit count and ``mcx_ancilla``, Toffoli/Clifford/exact T
            counts, the total of rotations pending synthesis with per-axis statistics,
            T cost estimates, per-atom ``atoms`` details, per-resource QRAM
            query/write counts and totals, plus ``gate_total``.
        """
        axes: Counter[str] = Counter()
        for (axis, _), count in self.rotations.counts.items():
            axes[axis] += count
        return {
            "complete": self.complete,
            "qubits": self.qubits,
            "qubits_lower_bound": self.qubits_lower_bound,
            "unknown_workspace": self.unknown_workspace,
            "oracle_calls": [
                {**asdict(call), "count": count}
                for call, count in sorted(self.oracle_calls.items()) if count
            ],
            "mcx_ancilla": self.mcx_ancilla,
            "toffoli": self.toffoli,
            "clifford": self.clifford,
            "t_exact": self.t_exact,
            "rotations": self.rotations.total,
            "rotation_axes": dict(axes),
            "rotation_counts": [
                {"axis": axis, "angle": angle, "count": count}
                for (axis, angle), count in sorted(self.rotations.counts.items()) if count
            ],
            "t_synthesis": self.t_synthesis(epsilon),
            "t_total": self.t_total(epsilon),
            "atoms": dict(sorted(self.atoms.items())),
            "qram_queries": dict(sorted(self.qram_queries.items())),
            "qram_total": self.qram_total,
            "qram_writes": dict(sorted(self.qram_writes.items())),
            "qram_write_total": self.qram_write_total,
            "gate_total": self.gate_total,
            "epsilon": epsilon,
        }


@dataclass
class _Cost:
    """Module-relative cost; resource arguments are substituted along each call edge."""

    atoms: Counter[str] = field(default_factory=Counter)
    rotations: RotationCounts = field(default_factory=RotationCounts)
    queries: Counter[str] = field(default_factory=Counter)
    writes: Counter[str] = field(default_factory=Counter)
    calls: Counter[OracleCall] = field(default_factory=Counter)

    def add(self, other: _Cost, factor: int = 1, resources: dict[str, str] | None = None) -> None:
        mapping = resources or {}
        for atom, count in other.atoms.items():
            self.atoms[atom] += count * factor
        self.rotations.add_scaled(other.rotations, factor)
        for target, source in ((self.queries, other.queries), (self.writes, other.writes)):
            for resource, count in source.items():
                target[mapping.get(resource, resource)] += count * factor
        for call, count in other.calls.items():
            actual = replace(call, resources=tuple(mapping.get(r, r) for r in call.resources))
            self.calls[actual] += count * factor


def estimate_resources(program: Program, *, require_closed: bool = True) -> ResourceEstimate:
    """Composable resource estimation: memoized by module and control count, with Repeat counts multiplied symbolically.

    Args:
        program: RIR program to estimate; module implementations must be closed unless validation is relaxed.
        require_closed: when false, open declaration modules may remain; by default all must be closed.

    Returns:
        ResourceEstimate: aggregated resource ledger; the qubit count includes the workspace, and the maximum control count
        determines ``mcx_ancilla``.
    """
    program = validate(program, require_closed=require_closed)
    modules = program.module_map
    workspace = workspace_table(program)
    memo: dict[tuple[str, int, bool], _Cost] = {}
    max_controls = 0

    def primitive_cost(
        node: Primitive, n_controls: int
    ) -> _Cost:
        """Compute the gate-level cost of a single primitive instruction under the given control count."""
        nonlocal max_controls
        counts: Counter[str]
        rotations: list[tuple[str, float]]
        queries: Counter[str]
        writes: Counter[str]
        counts, rotations, queries, writes = Counter(), [], Counter(), Counter()
        widths = [ref.width for ref in node.operands]
        if node.op in {"xor", "swap"}:
            for _ in range(widths[0]):
                counts += mcx_counts(n_controls + 1)
                if node.op == "swap":
                    counts += mcx_counts(n_controls + 1)
                    counts += mcx_counts(n_controls + 1)
            max_controls = max(max_controls, n_controls + 1)
        elif node.op == "add_const":
            value = cast(int, node.value) % (1 << widths[0])
            for offset in range(widths[0]):
                if (value >> offset) & 1:
                    for i in reversed(range(offset + 1, widths[0])):
                        counts += mcx_counts(n_controls + (i - offset))
                        max_controls = max(max_controls, n_controls + (i - offset))
                    counts += mcx_counts(n_controls)
        elif node.op == "gphase":
            theta = cast(float, node.angle)
            if n_controls:
                part, rot = _controlled_u3_counts(n_controls - 1, 0.0, 0.0, theta, 0.0)
                counts += part
                rotations += rot
                max_controls = max(max_controls, n_controls - 1)
            else:
                for _ in range(2):
                    part, rot = classify_rz(theta)
                    counts += part
                    rotations += rot
                counts["x"] += 2
        else:
            for _ in range(widths[0]):
                part, rot = _gate_counts(node.op, node.angle, n_controls)
                counts += part
                rotations += rot
            max_controls = max(max_controls, n_controls)
        return _Cost(counts, RotationCounts(rotations), queries, writes)

    def body_cost(
        nodes: tuple[Instruction, ...], n_controls: int, n_zeros: int, inverse: bool
    ) -> _Cost:
        """Aggregate instruction body cost, handling zero-valued control bit flips, repetition and module calls."""
        total = _Cost()
        for node in nodes:
            if isinstance(node, Primitive):
                if n_zeros and not (
                    node.op == "add_const"
                    and not cast(int, node.value) % (1 << node.operands[0].width)
                ):
                    total.atoms["x"] += 2 * n_zeros
                total.add(primitive_cost(node, n_controls))
            elif isinstance(node, Load):
                total.queries[node.resource] += 1
            elif isinstance(node, Store):
                total.writes[node.resource] += 1
            elif isinstance(node, Call):
                if n_zeros:
                    total.atoms["x"] += 2 * n_zeros
                resource_map = dict(zip(
                    (r.name for r in modules[node.module].resources), node.resources, strict=True,
                ))
                total.add(module_cost(node.module, n_controls, inverse), resources=resource_map)
            elif isinstance(node, Repeat):
                if node.count and node.body:
                    if n_zeros:
                        total.atoms["x"] += 2 * n_zeros
                    total.add(body_cost(node.body, n_controls, 0, inverse), node.count)
            elif isinstance(node, Control):
                total.add(
                    body_cost(
                        node.body,
                        n_controls + node.register.width,
                        n_zeros + node.register.width - bin(node.value).count("1"),
                        inverse,
                    ),
                )
            elif isinstance(node, Adjoint):
                total.add(body_cost(node.body, n_controls, n_zeros, not inverse))
        return total

    def module_cost(
        name: str, n_controls: int, inverse: bool
    ) -> _Cost:
        """Memoize by module, control count and adjoint context; repetitions are not expanded."""
        key = (name, n_controls, inverse)
        if key in memo:
            return memo[key]
        module = modules[name]
        if module.body is None:
            call = OracleCall(name, n_controls, inverse, tuple(r.name for r in module.resources))
            memo[key] = _Cost(calls=Counter({call: 1}))
        else:
            memo[key] = body_cost(module.body, n_controls, 0, inverse)
        return memo[key]

    total = module_cost(program.entry, 0, False)
    missing = tuple(item.name for item in unresolved(program))
    known_qubits = sum(r.type.width for r in program.main.registers) + workspace[program.entry]
    return ResourceEstimate(
        qubits=None if missing else known_qubits,
        atoms=total.atoms,
        rotations=total.rotations,
        qram_queries=total.queries,
        qram_writes=total.writes,
        mcx_ancilla=max(0, max_controls - 1),
        oracle_calls=total.calls,
        unknown_workspace=missing,
        qubits_lower_bound=known_qubits,
    )
