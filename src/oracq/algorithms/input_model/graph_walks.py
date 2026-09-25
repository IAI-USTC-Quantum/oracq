"""Graph oracle input model and the Szegedy/MNRS quantum walk search framework.

The graph is supplied as a neighbor-table oracle:
``|v>|j>|c> -> |v>|j>|c XOR N(v,j)>``, with the out-edge table uniformly padded
to D=2**degree_bits columns (usually padded with self-loops). The Szegedy
walk step is defined on the (current, peer, index) bipartite walk space and
consists of two reflections built from neighbor superposition preparations,
W = R_B R_A; the MNRS search skeleton is setup initial state preparation
followed by Repeat{ walk; marked phase flip }.

References: Szegedy 2004 (the spectrum of the walk operator);
Magniez–Nayak–Roland–Santha 2011 (Search via quantum walk, hitting time and
step scaling).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import cast

from oracq.algorithms.input_model.block_encoding import reflect_zero
from oracq.algorithms.input_model.contracts import (
    OracleView,
    positive_integer,
    require_instance,
    validate_signature,
)
from oracq.algorithms.input_model.interfaces import (
    StatePreparationProtocol,
    checked_state_preparation,
)
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
    XorDatabase,
    annotate,
    declare,
    invoke,
    resources_for,
)
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import QRAM, Bits, Ref, ValidationError, fuse


@dataclass(frozen=True)
class AdjacencyOracle(OracleView):
    """Neighbor-table query: ``|v>|j>|c> -> |v>|j>|c XOR N(v,j)>``.

    The out-edge table is arranged by vertex number, with exactly
    D=2**degree_bits entries per row; the caller pads the shortfall with
    self-loops or the like, to keep the neighbor superposition normalized as
    1/sqrt(D). The paradigm reuses database_xor, with address interpreted as
    the low-bit concatenation of vertex and index."""

    oracle_kind = "database_xor"
    operation: Operation

    def adjacency_oracle(self) -> AdjacencyOracle:
        """Return the ``AdjacencyOracle`` view of itself, as the adapter entry point of the adjacency oracle protocol.

        Returns:
            AdjacencyOracle: This view itself.
        """
        return self

    def __post_init__(self) -> None:
        """Validate the operation signature and structural constraints such as equal vertex/neighbor widths."""
        validate_signature(self.operation, ("vertex", "index", "neighbor"), "AdjacencyOracle")
        regs = {r.name: r.type for r in self.operation.module.registers}
        if regs["vertex"].width != regs["neighbor"].width:
            raise ValidationError("vertex and neighbor bit widths must agree")
        if regs["index"].width < 0:
            raise ValidationError("invalid index bit width")

    @property
    def vertex_bits(self) -> int:
        """Bit width of the vertex register ``vertex``."""
        return next(r.type.width for r in self.operation.module.registers if r.name == "vertex")

    @property
    def degree_bits(self) -> int:
        """Bit width of the out-edge index register ``index``."""
        return next(r.type.width for r in self.operation.module.registers if r.name == "index")

    @property
    def vertices(self) -> int:
        """Total number of addressable vertices, i.e. ``2**vertex_bits``."""
        return 1 << self.vertex_bits

    @property
    def degree(self) -> int:
        """Number ``D`` of out-edge table columns per vertex after padding, i.e. ``2**degree_bits``."""
        return 1 << self.degree_bits

    def xor_database(self) -> XorDatabase:
        """Treat the neighbor table as an XOR database with address=``vertex||index`` and data=neighbor.

        Returns:
            XorDatabase: XOR database view wrapping the same neighbor-table operation.
        """
        from oracq.algorithms.input_model.oracles import XorDatabase

        v, g = self.vertex_bits, self.degree_bits
        b = Builder(
            _name("adjacency_as_database", self.operation),
            {"address": Bits(v + g), "data": Bits(v)},
            resources_for(("adj", self.operation)),
        )
        invoke(
            b,
            self.operation,
            "adj",
            vertex=b["address"][:v],
            index=b["address"][v:],
            neighbor=b["data"],
        )
        return XorDatabase(annotate(b.finish(), "database_xor"))


def as_adjacency(value: AdjacencyOracle | Operation) -> AdjacencyOracle:
    """Adapt a graph adjacency input into an ``AdjacencyOracle`` view.

    Args:
        value: Returned as-is when already an ``AdjacencyOracle``; otherwise it
        must be a complete ``Operation`` with the ``(vertex, index, neighbor)``
        signature.

    Returns:
        AdjacencyOracle: Adjacency oracle view wrapping the given operation.

    Raises:
        ValidationError: value is neither an ``AdjacencyOracle`` nor an ``Operation``.
    """
    if isinstance(value, AdjacencyOracle):
        return value
    return AdjacencyOracle(require_instance(value, Operation, "adjacency"))


def _as_operation(value: object, path: str) -> Operation:
    """Unify a view or a bare ``Operation`` into a program-validated ``Operation``."""
    operation = getattr(value, "operation", None)
    if operation is None:
        operation = value
    require_instance(operation, Operation, path)
    cast("Operation", operation).program()
    return cast("Operation", operation)


def abstract_adjacency(
    vertex_bits: int, degree_bits: int, *, name: str | None = None
) -> AdjacencyOracle:
    """Open declaration of graph adjacency; the body is empty and can be bound to a gate/QRAM implementation via bind.

    Args:
        vertex_bits: Bit width of the vertex register, range 1..32.
        degree_bits: Bit width of the out-edge index register, range 0..32.
        name: Name of the declaration slot; ``GraphAdjacency_{vertex_bits}_{degree_bits}`` by default.

    Returns:
        AdjacencyOracle: Adjacency declaration with an empty body awaiting a bound implementation.
    """
    positive_integer(vertex_bits, "abstract_adjacency.vertex_bits", maximum=32)
    positive_integer(degree_bits, "abstract_adjacency.degree_bits", minimum=0, maximum=32)
    return AdjacencyOracle(
        declare(
            name or f"GraphAdjacency_{vertex_bits}_{degree_bits}",
            {
                "vertex": Bits(vertex_bits),
                "index": Bits(degree_bits),
                "neighbor": Bits(vertex_bits),
            },
            paradigm="database_xor",
            attributes={"oracle_role": "graph_adjacency"},
        )
    )


def _check_table(
    neighbors: Sequence[Sequence[int]],
) -> tuple[tuple[tuple[int, ...], ...], int, int]:
    """Validate that the neighbor table is a nonempty rectangle whose entries are legal vertex numbers, returning the normalized row table."""
    rows = tuple(tuple(row) for row in neighbors)
    if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValidationError("the neighbor table must be a nonempty rectangle")
    n = len(rows)
    for row in rows:
        for target in row:
            if type(target) is not int or not 0 <= target < n:
                raise ValidationError("neighbor table entries must be legal vertex numbers")
    return rows, n, len(rows[0])


def _widths(n: int, d: int) -> tuple[int, int]:
    """Encoding bit widths needed for the vertex count ``n`` and the out-edge column count ``d``."""
    return max(1, (n - 1).bit_length()), (d - 1).bit_length()


def gate_adjacency(
    neighbors: Sequence[Sequence[int]], *, name: str | None = None
) -> AdjacencyOracle:
    """Gate-level implementation for a small neighbor table: a controlled X network per (v, j) branch.

    Args:
        neighbors: Nonempty rectangular neighbor table with rows for vertices
        0..V-1 and D entries per row; entries are legal vertex numbers.
        name: Name of the generated operation; derived from the table content by default.

    Returns:
        AdjacencyOracle: Adjacency oracle implemented as a controlled X network.
    """
    rows, n, d = _check_table(neighbors)
    v, g = _widths(n, d)
    b = Builder(
        name or _name("gate_adjacency", rows),
        {"vertex": Bits(v), "index": Bits(g), "neighbor": Bits(v)},
    )
    for vertex, row in enumerate(rows):
        for index, target in enumerate(row):
            if target:
                with b.control(fuse(b["vertex"], b["index"]), vertex | (index << v)):
                    for bit in range(v):
                        if (target >> bit) & 1:
                            b.x(b["neighbor"][bit])
    return AdjacencyOracle(
        annotate(b.finish(), "database_xor", implementation="gate_neighbor_table")
    )


def qram_adjacency(
    vertex_bits: int, degree_bits: int, *, name: str | None = None
) -> AdjacencyOracle:
    """QRAM implementation of the neighbor table; the address is ``vertex | (index << vertex_bits)`` and the word is the neighbor number.

    Args:
        vertex_bits: Bit width of the vertex register, range 1..32.
        degree_bits: Bit width of the out-edge index register, range 0..32.
        name: Name of the generated operation; ``qram_adjacency_{vertex_bits}_{degree_bits}`` by default.

    Returns:
        AdjacencyOracle: Adjacency oracle backed by the table bank.
    """
    positive_integer(vertex_bits, "qram_adjacency.vertex_bits", maximum=32)
    positive_integer(degree_bits, "qram_adjacency.degree_bits", minimum=0, maximum=32)
    b = Builder(
        name or f"qram_adjacency_{vertex_bits}_{degree_bits}",
        {
            "vertex": Bits(vertex_bits),
            "index": Bits(degree_bits),
            "neighbor": Bits(vertex_bits),
        },
        {"table": QRAM(vertex_bits + degree_bits, vertex_bits)},
    )
    b.qram("table", fuse(b["vertex"], b["index"]), b["neighbor"])
    return AdjacencyOracle(
        annotate(b.finish(), "database_xor", implementation="qram_neighbor_table")
    )


def _neighbor_reflection(
    b: Builder, adjacency: AdjacencyOracle, prefix: str, vertex: Ref, peer: Ref, index: Ref
) -> None:
    """Reflection about ``span_v{ A|v>|0,0> }``, with ``A|v>|0,0> = 1/sqrt(D) sum_j |v>|N(v,j)>|j>``.

    The neighbor state preparation A = Adj · H_index is unitary as a whole
    (Adj is self-inverse), so ``R = A(2|0><0|-I)A†`` can be assembled in the
    order Adj, H, Phi, H, Adj; the index register stays inside the walk
    space."""
    invoke(b, adjacency.operation, prefix, vertex=vertex, index=index, neighbor=peer)
    b.h(index)
    reflect_zero(b, fuse(peer, index), positive=True)
    b.h(index)
    invoke(b, adjacency.operation, prefix, vertex=vertex, index=index, neighbor=peer)


def szegedy_walk(
    adjacency: AdjacencyOracle | Operation, *, name: str | None = None
) -> Operation:
    """Szegedy walk step W = R_B R_A, acting on the (current, peer, index) walk space.

    R_A is the reflection about the subspace of neighbor superpositions on the
    current end, and R_B swaps the roles of the two ends (equivalent to SWAP
    conjugation). The graph need not be regular; when irregular, the padding
    scheme determines the effective transition.

    Args:
        adjacency: Adjacency oracle, or an operation with the ``(vertex, index, neighbor)`` signature.
        name: Name of the generated operation; derived from the adjacency operation by default.

    Returns:
        Operation: The unitary operation of one walk step.
    """
    adjacency = as_adjacency(adjacency)
    v, g = adjacency.vertex_bits, adjacency.degree_bits
    b = Builder(
        name or _name("szegedy_walk", adjacency.operation),
        {"current": Bits(v), "peer": Bits(v), "index": Bits(g)},
        resources_for(("adj", adjacency.operation)),
        attributes={
            "algorithm": "szegedy_walk",
            "vertex_bits": v,
            "degree_bits": g,
            "composition": "R_B.R_A",
        },
    )
    _neighbor_reflection(b, adjacency, "adj", b["current"], b["peer"], b["index"])
    _neighbor_reflection(b, adjacency, "adj", b["peer"], b["current"], b["index"])
    return b.finish()


def szegedy_setup(
    adjacency: AdjacencyOracle | Operation, *, name: str | None = None
) -> StatePreparation:
    """Prepare the MNRS initial state ``1/sqrt(N) sum_v A|v>``: a uniform superposition on current followed by neighbor state preparation.

    The target layout matches the register declaration order of szegedy_walk:
    the low vertex_bits bits are current, the next vertex_bits bits are peer,
    and the high degree_bits bits are index.

    Args:
        adjacency: Adjacency oracle, or an operation with the ``(vertex, index, neighbor)`` signature.
        name: Name of the generated operation; derived from the adjacency operation by default.

    Returns:
        StatePreparation: Preparation view preparing the MNRS initial state.
    """
    adjacency = as_adjacency(adjacency)
    v, g = adjacency.vertex_bits, adjacency.degree_bits
    b = Builder(
        name or _name("szegedy_setup", adjacency.operation),
        {"target": Bits(2 * v + g), "work": Bits(0)},
        resources_for(("adj", adjacency.operation)),
    )
    current, peer, index = b["target"][:v], b["target"][v : 2 * v], b["target"][2 * v :]
    b.h(current)
    b.h(index)
    invoke(b, adjacency.operation, "adj", vertex=current, index=index, neighbor=peer)
    return StatePreparation(
        annotate(
            b.finish(),
            "state_prep_isometry",
            zero_input=True,
            clean_work=True,
            implementation="uniform_vertex_plus_neighbor",
        )
    )


def quantum_walk_search(
    setup: StatePreparationProtocol,
    walk: Operation,
    marked: Operation,
    steps: int,
    *,
    name: str | None = None,
) -> Operation:
    """MNRS marked-vertex search skeleton: setup followed by Repeat{ walk; marked phase flip }.

    Args:
        setup: Initial state preparation (StatePreparation protocol), with target being the whole walk space.
        walk: Walk step handle; its registers map in declaration order to consecutive slices of setup.target.
        marked: Phase oracle (target with optional work) acting on the first target.width bits of the walk space, which is the current register for szegedy_walk.
        steps: Nonnegative iteration count; choose with reference to suggest_steps.

    All three handles may be open declarations; the generated program becomes
    executable after bind closes the implementations in batches.

    Returns:
        Operation: The search program of repeated walking and phase flipping.
    """
    positive_integer(steps, "quantum_walk_search.steps", minimum=0)
    prep = checked_state_preparation(setup, path="quantum_walk_search.setup")
    walk_op = _as_operation(walk, "quantum_walk_search.walk")
    marked_op = _as_operation(marked, "quantum_walk_search.marked")
    width = sum(r.type.width for r in walk_op.module.registers)
    if width != prep.width:
        raise ValidationError("walk and setup disagree on the walk space width")
    marked_regs = {r.name: r.type for r in marked_op.module.registers}
    if "target" not in marked_regs or set(marked_regs) - {"target", "work"}:
        raise ValidationError("the marked phase oracle requires a target interface with optional work")
    marked_width = marked_regs["target"].width
    marked_work = marked_regs.get("work", Bits(0)).width
    if marked_width > width:
        raise ValidationError("the marked target exceeds the walk space width")
    b = Builder(
        name or _name("quantum_walk_search", prep.operation, walk_op, marked_op, steps),
        {"target": Bits(width), "work": Bits(prep.work_width + marked_work)},
        resources_for(("setup", prep.operation), ("walk", walk_op), ("marked", marked_op)),
        attributes={
            "algorithm": "quantum_walk_search",
            "framework": "MNRS",
            "steps": steps,
            "walk_width": width,
        },
    )
    setup_work, marked_ref = b["work"][: prep.work_width], b["work"][prep.work_width :]
    invoke(b, prep.operation, "setup", target=b["target"], work=setup_work)
    with b.repeat(steps):
        cursor = 0
        arguments: dict[str, Ref] = {}
        for reg in walk_op.module.registers:
            arguments[reg.name] = (
                b["target"][cursor : cursor + reg.type.width].reinterpret(reg.type.kind)
            )
            cursor += reg.type.width
        invoke(b, walk_op, "walk", **arguments)
        marked_arguments = {"target": b["target"][:marked_width]}
        if marked_work:
            marked_arguments["work"] = marked_ref
        invoke(b, marked_op, "marked", **marked_arguments)
    return b.finish()


def transition_matrix(neighbors: Sequence[Sequence[int]]) -> tuple[tuple[float, ...], ...]:
    """Build the classical random-walk transition matrix from the neighbor table, ``P[v][u] = |{j: N(v,j)=u}| / D``.

    Args:
        neighbors: Nonempty rectangular neighbor table whose entries are legal vertex numbers.

    Returns:
        tuple[tuple[float, ...], ...]: Row-major transition probability matrix, with each row summing to 1.
    """
    rows, n, d = _check_table(neighbors)
    result = []
    for row in rows:
        counts = [0.0] * n
        for target in row:
            counts[target] += 1.0 / d
        result.append(tuple(counts))
    return tuple(result)


def hitting_times(
    transition: Sequence[Sequence[float]], marked: Iterable[int]
) -> tuple[float, ...]:
    """Expected number of steps for each vertex to first reach the marked set (zero for marked vertices).

    Solves the linear system (I - P_free) h = 1, where P_free is the
    sub-transition matrix with the marked rows and columns removed.

    Args:
        transition: Nonnegative square transition matrix with each row summing to 1.
        marked: Set of marked vertex numbers, which must be nonempty and in range.

    Returns:
        tuple[float, ...]: Expected hitting times ordered by vertex number.
    """
    n = len(transition)
    if not n or any(len(row) != n for row in transition):
        raise ValidationError("the transition matrix must be a nonempty square matrix")
    for row in transition:
        if any(type(p) not in (int, float) or p < 0 for p in row):
            raise ValidationError("transition probabilities must be nonnegative real numbers")
        if abs(sum(row) - 1.0) > 1e-9:
            raise ValidationError("each row of the transition matrix must sum to 1")
    marked_set = set(marked)
    if not marked_set or any(type(m) is not int or not 0 <= m < n for m in marked_set):
        raise ValidationError("marked must be a nonempty set of legal vertex numbers")
    free = [v for v in range(n) if v not in marked_set]
    if not free:
        return tuple(0.0 for _ in range(n))
    size = len(free)
    system = [
        [
            (1.0 if i == j else 0.0) - float(transition[free[i]][free[j]])
            for j in range(size)
        ]
        + [1.0]
        for i in range(size)
    ]
    for col in range(size):
        pivot = max(range(col, size), key=lambda r: abs(system[r][col]))
        if abs(system[pivot][col]) < 1e-12:
            raise ValidationError("the marked set is unreachable from some vertices and the hitting time diverges")
        system[col], system[pivot] = system[pivot], system[col]
        for row in range(size):  # type: ignore[assignment]
            if row != col:
                factor = system[cast("int", row)][col] / system[col][col]
                for k in range(col, size + 1):
                    system[cast("int", row)][k] -= factor * system[col][k]
    hits = {free[i]: system[i][size] / system[i][i] for i in range(size)}
    return tuple(0.0 if v in marked_set else hits[v] for v in range(n))


def suggest_steps(transition: Sequence[Sequence[float]], marked: Iterable[int]) -> int:
    """MNRS suggested walk step count: ceil(pi/4 * sqrt(H_avg)), where H_avg is the average hitting time.

    Matches the Grover limit: on the complete graph ``H_avg = N/|M|`` and the
    step count recovers ``pi/4 * sqrt(N/|M|)``.

    Args:
        transition: Nonnegative square transition matrix with each row summing to 1.
        marked: Set of marked vertex numbers, which must be nonempty and in range.

    Returns:
        int: Suggested iteration count, at least 1; zero when there are no free vertices.
    """
    hits = hitting_times(transition, marked)
    n = len(transition)
    free = n - len(set(marked))
    if not free:
        return 0
    average = sum(hits) / free
    return max(1, math.ceil(math.pi / 4 * math.sqrt(average)))
