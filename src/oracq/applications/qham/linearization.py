"Lazily generate the QCL by HAM order and ordered tensor words; spatial matrices are not materialized."

from __future__ import annotations

import json
import math
from collections.abc import Iterator, Sequence
from dataclasses import dataclass

from oracq.applications.qham.pde import PolynomialPDE
from oracq.infrastructure.ir import ValidationError


def compositions(total: int, length: int) -> Iterator[tuple[int, ...]]:
    """Lazily enumerate all ordered splits of ``length`` non-negative integers summing to ``total``, in increasing first-component order.

    Args:
        total: Target sum, a non-negative integer.
        length: Number of components in each split.

    Returns:
        Iterator[tuple[int, ...]]: Lazy generator yielding one split at a
        time without materializing all splits.

    Yields:
        tuple[int, ...]: Ordered tuples of length ``length`` whose elements
            sum to ``total``; when ``length`` is 0 the empty tuple is yielded
            only if ``total`` is 0.
    """
    if length == 0:
        if total == 0:
            yield ()
    elif length == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, length - 1):
                yield (first, *rest)


@dataclass(frozen=True)
class HomotopyWeight:
    """Coefficient on a QCL linear edge that varies with the homotopy parameter ``eta``.

    Attributes:
        kind: Weight kind: ``one`` (identically 1), ``correction``, or
            ``physical``.
        power: Exponent used by the ``correction`` and ``physical`` kinds.
    """

    kind: str = "one"
    power: int = 0

    def evaluate(self, eta: complex) -> complex:
        """Evaluate the weight at the homotopy parameter ``eta``.

        Args:
            eta: Homotopy parameter value.

        Returns:
            ``one`` is identically 1, ``correction`` is
            ``-eta*(1+eta)**power``, and ``physical`` is
            ``1-(1+eta)**power``.

        Raises:
            ValidationError: ``kind`` is not a supported weight kind.
        """
        if self.kind == "one":
            return 1.0
        if self.kind == "correction":
            return -eta * (1 + eta) ** self.power
        if self.kind == "physical":
            return 1 - (1 + eta) ** self.power
        raise ValidationError("unknown homotopy weight")

    def formula(self) -> str:
        """Return the formula text of the weight, with the homotopy parameter written ``eta``.

        Returns:
            str: Infix expression matching ``evaluate``, e.g.
            ``-eta*(1+eta)^3``.
        """
        if self.kind == "one":
            return "1"
        if self.kind == "correction":
            return f"-eta*(1+eta)^{self.power}"
        return f"1-(1+eta)^{self.power}"


@dataclass(frozen=True)
class Block:
    """One block variable of the QCL finite closure.

    Attributes:
        kind: Block kind; ``physical`` is the physical output block, the sum
            of the HAM order components, and ``tensor`` is a tensor block on
            independent coordinates.
        orders: Tensor word, i.e. the ordered tuple of HAM orders of the
            factors; always the empty tuple for physical blocks, and the empty
            tuple with ``tensor`` denotes the constant block introduced by
            forcing homogenization.
    """

    kind: str
    orders: tuple[int, ...] = ()

    @property
    def rank(self) -> int:
        """Number of tensor factors of the block; 1 for physical blocks and 0 for the empty-word constant block."""
        return 1 if self.kind == "physical" else len(self.orders)

    @property
    def label(self) -> str:
        """Return the human-readable label of the block.

        Physical blocks are ``u_sum``, the empty-word constant block is
        ``one``, and other tensor blocks are ``Y_`` followed by the
        underscore-joined orders.
        """
        if self.kind == "physical":
            return "u_sum"
        if not self.orders:
            return "one"
        return "Y_" + "_".join(map(str, self.orders))


@dataclass(frozen=True)
class Coupling:
    """One linear edge of the QCL closure, delivering a column block's contribution into the row block's equation.

    Attributes:
        row: Target block of the edge (the owning equation row).
        column: Source block of the edge (the input variable).
        operator: Name of the applied multilinear port, e.g. ``L``, ``F``, or
            ``B_i``.
        position: Which tensor factor of the row block the port acts on,
            replacing it; counted from 0.
        arity: Number of port inputs.
        weight: Homotopy weight of this edge, identically 1 by default.
    """

    row: Block
    column: Block
    operator: str
    position: int
    arity: int
    weight: HomotopyWeight = HomotopyWeight()


@dataclass(frozen=True)
class QHAMPlan:
    """QCL finite closure plan of a polynomial PDE at a HAM truncation order.

    The closure is bounded by the weight criterion
    ``p*sum(orders)+len(orders) <= p*m+1``, where ``p = max(1, D-1)``, ``D``
    is the PDE polynomial degree, and ``m`` is the truncation order. The block
    set and each row rule are enumerated lazily; the closure matrix is never
    materialized.

    Attributes:
        pde: The ``PolynomialPDE``, validated at construction.
        order: HAM truncation order ``m``, must be a non-negative integer.
        version: Serialization format version, currently ``0.1``.

    Raises:
        ValidationError: ``pde`` fails validation, or ``order`` is not a
            non-negative integer, or ``version`` is not ``0.1``.
    """

    pde: PolynomialPDE
    order: int
    version: str = "0.1"

    def __post_init__(self) -> None:
        """Validate the PDE, truncation order, and version, raising ``ValidationError`` when invalid."""
        self.pde.validate()
        if self.version != "0.1" or type(self.order) is not int or self.order < 0:
            raise ValidationError("the QHAM truncation order must be a non-negative integer")

    @property
    def grade(self) -> int:
        """Degree scaling of the closure weight criterion, ``p = max(1, D-1)``, with ``D`` the PDE polynomial degree."""
        return max(1, self.pde.degree - 1)

    @property
    def max_rank(self) -> int:
        """Closure weight bound ``p*m+1``, the maximum tensor rank allowed for tensor blocks."""
        return self.grade * self.order + 1

    @property
    def has_forcing(self) -> bool:
        """Whether the PDE has a constant forcing term, i.e. port ``F`` exists and the closure needs the empty-word constant block."""
        return any(port.name == "F" for port in self.pde.ports)

    def weight(self, word: Sequence[int]) -> int:
        """Compute the closure weight ``p*sum(word)+len(word)`` of a tensor word.

        Args:
            word: Tuple of the HAM orders of the factors.

        Returns:
            int: The closure weight of this tensor word.
        """
        return self.grade * sum(word) + len(word)

    def contains(self, block: Block) -> bool:
        """Test whether a block belongs to this plan's finite closure.

        Physical blocks must carry no orders; the empty-word tensor block
        belongs only when forcing exists; other tensor blocks require every
        order to be a non-negative integer with closure weight within the
        bound; block kinds other than these never belong to the closure.

        Args:
            block: The block to check.

        Returns:
            bool: True when the block belongs to the closure.
        """
        if block.kind == "physical":
            return block.orders == ()
        if block.kind != "tensor":
            return False
        if not block.orders:
            return self.has_forcing
        return (
            all(type(a) is int and a >= 0 for a in block.orders)
            and self.weight(block.orders) <= self.max_rank
        )

    def rank_count(self, rank: int) -> int:
        """Count the closure tensor blocks of a given tensor rank.

        Args:
            rank: Number of tensor factors, counted from 1.

        Returns:
            int: Number of tensor words of this rank satisfying the weight
            constraint; 0 when the rank is below 1 or above the bound.
        """
        limit = (self.max_rank - rank) // self.grade
        return math.comb(limit + rank, rank) if limit >= 0 and rank >= 1 else 0

    @property
    def block_count(self) -> int:
        """Total closure block count: the physical block, the optional empty-word constant block, plus the tensor blocks of each rank."""
        return (
            1 + int(self.has_forcing) + sum(self.rank_count(k) for k in range(1, self.max_rank + 1))
        )

    def blocks(self, *, max_blocks: int | None = None) -> Iterator[Block]:
        """Lazily enumerate all blocks of the closure in canonical order.

        The physical block comes first, tensor blocks follow by ascending rank
        and, within a rank, by ascending word sum (ties follow the
        ``compositions`` enumeration order), and the empty-word constant
        block, when forcing exists, comes last.

        Args:
            max_blocks: Budget for explicitly enumerated blocks; when the
                total block count exceeds it an error is raised immediately,
                pointing to lazy queries instead.

        Returns:
            Iterator[Block]: Lazy generator yielding closure blocks one at a
            time without materializing the full block list.

        Yields:
            Block: A block of the closure.

        Raises:
            ValidationError: ``max_blocks`` is given and the closure block
                count exceeds the budget.
        """
        if max_blocks is not None and self.block_count > max_blocks:
            raise ValidationError("explicit block enumeration exceeds the budget; QHAMPlan remains lazily queryable")
        yield Block("physical")
        for rank in range(1, self.max_rank + 1):
            limit = (self.max_rank - rank) // self.grade
            for total in range(limit + 1):
                for word in compositions(total, rank):
                    yield Block("tensor", word)
        if self.has_forcing:
            yield Block("tensor")

    def row_terms(self, row: Block) -> tuple[Coupling, ...]:
        """Enumerate all linear edges of a row block within the closure.

        Three edge classes are covered: the same-word linear edge ``L``;
        forcing insertion edges ``F`` (physical rows take the empty-word
        constant block, tensor rows take the word with that zero-order factor
        removed); and edges replacing a higher-order factor with a lower-order
        tuple and applying a nonlinear port, carrying a ``physical`` or
        ``correction`` homotopy weight. The row of the empty-word constant
        block has no edges.

        Args:
            row: Row block, must belong to the closure.

        Returns:
            tuple[Coupling, ...]: All edges of this row; each edge's ``row``
            field equals the input.

        Raises:
            ValidationError: The row block is outside the finite closure.
            AssertionError: A generated column block escapes the closure,
                i.e. the closure construction invariant is broken.
        """
        if not self.contains(row):
            raise ValidationError("the QCL row does not belong to the finite closure")
        if row.kind == "tensor" and not row.orders:
            return ()
        ports = self.pde.ports
        linear = any(p.name == "L" for p in ports)
        nonlinear = [p for p in ports if p.arity >= 2]
        result: list[Coupling] = []
        if row.kind == "physical":
            if linear:
                result.append(Coupling(row, row, "L", 0, 1))
            if self.has_forcing:
                result.append(Coupling(row, Block("tensor"), "F", 0, 0))
            for total in range(self.order):
                for port in nonlinear:
                    for word in compositions(total, port.arity):
                        result.append(
                            Coupling(
                                row,
                                Block("tensor", word),
                                port.name,
                                0,
                                port.arity,
                                HomotopyWeight("physical", self.order - total),
                            )
                        )
        else:
            for position, order in enumerate(row.orders):
                if linear:
                    result.append(Coupling(row, row, "L", position, 1))
                if order == 0:
                    if self.has_forcing:
                        word = row.orders[:position] + row.orders[position + 1 :]
                        result.append(Coupling(row, Block("tensor", word), "F", position, 0))
                else:
                    for total in range(order):
                        for port in nonlinear:
                            for word in compositions(total, port.arity):
                                source = row.orders[:position] + word + row.orders[position + 1 :]
                                result.append(
                                    Coupling(
                                        row,
                                        Block("tensor", source),
                                        port.name,
                                        position,
                                        port.arity,
                                        HomotopyWeight("correction", order - 1 - total),
                                    )
                                )
        if any(not self.contains(term.column) for term in result):
            raise AssertionError("QCL closure construction error")
        return tuple(result)

    def raw_dimension(self, dimension: int) -> int:
        """Compute the total dimension of the closure linear system.

        Args:
            dimension: Dimension of a single block's local state space.

        Returns:
            int: Sum of the local dimensions of the physical block and each
            rank's tensor blocks (``dimension**rank``), plus the single
            coordinate of the empty-word constant block when forcing exists.
        """
        return (
            dimension
            + sum(self.rank_count(k) * dimension**k for k in range(1, self.max_rank + 1))
            + int(self.has_forcing)
        )

    def offset(self, block: Block, dimension: int) -> int:
        """Return the global start index of local coordinate 0 of a block within the closure system.

        Blocks are laid out in the canonical order of ``blocks``: the physical
        block first, each rank's tensor blocks occupying ``dimension**rank``
        consecutive coordinates, and the empty-word constant block last with a
        single coordinate.

        Args:
            block: A block within the closure.
            dimension: Dimension of a single block's local state space.

        Returns:
            int: The starting global index of this block.

        Raises:
            ValidationError: The block is outside the closure.
        """
        if not self.contains(block):
            raise ValidationError("unknown QCL block")
        if block.kind == "physical":
            return 0
        if not block.orders:
            return self.raw_dimension(dimension) - 1
        rank, total = len(block.orders), sum(block.orders)
        offset = dimension + sum(self.rank_count(k) * dimension**k for k in range(1, rank))
        before = math.comb(total + rank - 1, rank) if total else 0
        remainder = total
        for pos, value in enumerate(block.orders[:-1]):
            slots = rank - pos - 1
            for first in range(value):
                before += math.comb(remainder - first + slots - 1, slots - 1)
            remainder -= value
        return offset + before * dimension**rank

    def locate(self, index: int, dimension: int) -> tuple[Block, int]:
        """Decompose a global index into its owning block and local index within the block; inverse of ``offset``.

        Args:
            index: Global index within the closure system.
            dimension: Dimension of a single block's local state space.

        Returns:
            tuple[Block, int]: The block containing the index and the local
            index within it.

        Raises:
            ValidationError: The index is outside the closure system's total
                dimension.
        """
        if not 0 <= index < self.raw_dimension(dimension):
            raise ValidationError("QCL raw index out of range")
        if index < dimension:
            return Block("physical"), index
        index -= dimension
        for rank in range(1, self.max_rank + 1):
            size = dimension**rank
            group = self.rank_count(rank) * size
            if index >= group:
                index -= group
                continue
            ordinal, local = divmod(index, size)
            total = 0
            while ordinal >= math.comb(total + rank - 1, rank - 1):
                ordinal -= math.comb(total + rank - 1, rank - 1)
                total += 1
            word: list[int] = []
            remaining = total
            for pos in range(rank - 1):
                slots = rank - pos - 1
                first = 0
                while ordinal >= math.comb(remaining - first + slots - 1, slots - 1):
                    ordinal -= math.comb(remaining - first + slots - 1, slots - 1)
                    first += 1
                word.append(first)
                remaining -= first
            word.append(remaining)
            return Block("tensor", tuple(word)), local
        return Block("tensor"), 0

    def summary(self) -> dict[str, object]:
        """Summarize the plan's key quantities for reports and diagnostics.

        Returns:
            dict: Dictionary with the version, truncation order, PDE degree,
            maximum tensor rank, total block count, whether forcing exists, a
            description of the closure criterion, and each port's arity.
        """
        return {
            "version": self.version,
            "order": self.order,
            "pde_degree": self.pde.degree,
            "max_tensor_rank": self.max_rank,
            "block_count": self.block_count,
            "homogeneous_constant": self.has_forcing,
            "closure": "(D-1)*sum(orders)+rank <= (D-1)*m+1; D<=1 uses grade 1",
            "operators": [
                {"name": p.name, "input_rank": p.arity, "output_rank": 1} for p in self.pde.ports
            ],
        }

    def dumps(self) -> str:
        """Serialize the plan to canonical JSON text.

        Returns:
            str: JSON text with ``version``, ``order``, and ``pde`` fields,
            keys sorted, indentation 2, ending with a newline.
        """
        return (
            json.dumps(
                {"version": self.version, "order": self.order, "pde": json.loads(self.pde.dumps())},
                ensure_ascii=False,
                sort_keys=True,
                indent=2,
            )
            + "\n"
        )

    @classmethod
    def loads(cls, text: str) -> QHAMPlan:
        """Rebuild a plan from the JSON text produced by ``dumps``.

        Args:
            text: JSON text.

        Returns:
            QHAMPlan: The rebuilt plan; validation runs again at construction.

        Raises:
            ValidationError: The JSON structure, field set, or data is
                invalid.
        """
        try:
            raw = json.loads(text)
            if set(raw) != {"version", "order", "pde"}:
                raise ValidationError("unknown QCL plan field")
            return cls(PolynomialPDE.loads(json.dumps(raw["pde"])), raw["order"], raw["version"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValidationError("invalid QCL plan: " + str(exc)) from exc
