"Small-scale mathematical witnesses: finite-difference ports and independent HAM/chain-rule evaluation, not quantum simulation."

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from fractions import Fraction
from functools import lru_cache

from oracq.applications.qham.linearization import QHAMPlan, compositions
from oracq.applications.qham.pde import Atom, Monomial, OperatorPort, PolynomialPDE
from oracq.infrastructure.ir import ValidationError


@lru_cache(maxsize=32)
def centered_weights(order: int) -> tuple[tuple[int, float], ...]:
    """Compute the weights of a unit-step one-dimensional centered difference stencil.

    The derivative interpolation conditions are solved by exact rational
    elimination on integer lattice points centered at the origin; the result is
    exact for polynomials of degree below the point count. Zero weights are
    dropped.

    Args:
        order: Derivative order, a non-negative integer.

    Returns:
        tuple: ``(lattice offset, weight)`` pairs where the weight is a float
        for unit step size; usage still requires division by
        ``spacing**order``. The result is cached per order.
    """
    if order == 0:
        return ((0, 1.0),)
    radius = max(1, (order + 1) // 2)
    points = list(range(-radius, radius + 1))
    size = len(points)
    rows = [
        [Fraction(x) ** k for x in points] + [Fraction(math.factorial(k) if k == order else 0)]
        for k in range(size)
    ]
    for col in range(size):
        pivot = next(i for i in range(col, size) if rows[i][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        scale = rows[col][col]
        rows[col] = [x / scale for x in rows[col]]
        for i in range(size):
            if i != col:
                scale = rows[i][col]
                rows[i] = [x - scale * y for x, y in zip(rows[i], rows[col], strict=True)]
    return tuple((x, float(rows[i][-1])) for i, x in enumerate(points) if rows[i][-1])


@dataclass(frozen=True)
class Grid:
    """Spatial discretization grid: axes, shape, spacing, and boundary type, with indexing and difference stencils.

    Linear grid indices use a mixed-radix layout with axis 0 in the least
    significant position; ``coordinates`` and ``address`` are inverse mappings
    of each other.

    Attributes:
        axes: Tuple of spatial axis names, must match the PDE axes.
        shape: Number of grid points per axis, all positive integers.
        spacing: Grid point spacing per axis, all positive finite values.
        boundary: Boundary type: ``periodic`` wraps out-of-range coordinates by
            the axis length; ``dirichlet_zero`` drops out-of-range stencil
            points by zero extension.

    Raises:
        ValidationError: Axis counts disagree, the shape or spacing is invalid,
            or the boundary type is unsupported.
    """
    axes: tuple[str, ...]
    shape: tuple[int, ...]
    spacing: tuple[float, ...]
    boundary: str = "periodic"
    _row_cache: dict[tuple[tuple[tuple[str, int], ...], int], tuple[tuple[int, float], ...]] = field(
        default_factory=dict, init=False, compare=False, repr=False, hash=False
    )

    def __post_init__(self) -> None:
        """Validate the grid shape, spacing, and boundary type, raising ``ValidationError`` when invalid."""
        if (
            not len(self.axes) == len(self.shape) == len(self.spacing)
            or any(type(n) is not int or n < 1 for n in self.shape)
            or any(x <= 0 or not math.isfinite(x) for x in self.spacing)
        ):
            raise ValidationError("invalid spatial grid shape or spacing")
        if self.boundary not in {"periodic", "dirichlet_zero"}:
            raise ValidationError("only periodic or zero-extension Dirichlet differences are currently provided")

    @property
    def size(self) -> int:
        """Total number of grid points, i.e. the product of all ``shape`` entries."""
        return math.prod(self.shape)

    @property
    def spatial_width(self) -> int:
        """Bit width required to index a spatial coordinate, at least 1.

        Spatial storage is aligned to ``1 << spatial_width``, which may exceed
        ``size``."""
        return max(1, (self.size - 1).bit_length())

    def coordinates(self, index: int) -> list[int]:
        """Decompose a linear grid index into per-axis coordinates; inverse of ``address``.

        Args:
            index: Linear grid index.

        Returns:
            list[int]: Per-axis coordinates, axis 0 in the least significant
            position, with the same length as ``shape``.
        """
        result: list[int] = []
        for n in self.shape:
            index, value = divmod(index, n)
            result.append(value)
        return result

    def address(self, coords: Sequence[int]) -> int:
        """Combine per-axis coordinates into a linear grid index.

        Args:
            coords: Sequence of per-axis coordinates; its length must match
                ``shape``.

        Returns:
            int: Linear grid index under the mixed-radix layout.
        """
        stride, result = 1, 0
        for n, x in zip(self.shape, coords, strict=True):
            result += stride * x
            stride *= n
        return result

    def derivative_row(
        self, derivative: tuple[tuple[str, int], ...], row: int
    ) -> tuple[tuple[int, float], ...]:
        """Compute the sparse difference row of a mixed spatial derivative at a given grid point.

        One-dimensional centered difference stencils are composed axis by axis,
        with weights already divided by the matching power of ``spacing``; the
        periodic boundary wraps out-of-range points while ``dirichlet_zero``
        drops them.

        Args:
            derivative: Derivative specification built from ``(axis name,
                order)`` pairs.
            row: Linear grid index.

        Returns:
            tuple: ``(column address, weight)`` pairs sorted by column address,
            zero weights dropped; an empty tuple when ``row`` is out of range.
            The result is cached by ``(derivative, row)``.

        Raises:
            ValueError: ``derivative`` references an axis not present in
                ``axes``.
        """
        key = (derivative, row)
        if key in self._row_cache:
            return self._row_cache[key]
        if not 0 <= row < self.size:
            return ()
        result = {row: 1.0}
        for axis, order in derivative:
            position = self.axes.index(axis)
            step = self.spacing[position]
            following: dict[int, float] = {}
            for source, value in result.items():
                origin = self.coordinates(source)
                for offset, weight in centered_weights(order):
                    coords = list(origin)
                    coords[position] += offset
                    if self.boundary == "periodic":
                        coords[position] %= self.shape[position]
                    elif not 0 <= coords[position] < self.shape[position]:
                        continue
                    address = self.address(coords)
                    following[address] = following.get(address, 0) + value * weight / (step**order)
            result = following
        self._row_cache[key] = tuple((i, v) for i, v in sorted(result.items()) if v)
        return self._row_cache[key]


def tensor_values(vectors: Iterable[Sequence[complex]]) -> list[complex]:
    """Compute the tensor product of several vectors flattened into a one-dimensional complex list.

    Factor 0 occupies the least significant position: the index stride of the
    k-th factor is the product of the lengths of the preceding factors; with no
    factors a single-element ``[1+0j]`` is returned.

    Args:
        vectors: Sequence of factor vectors.

    Returns:
        list[complex]: Flattened tensor product, with length equal to the
        product of the factor lengths.
    """
    result = [1.0 + 0j]
    for vector in vectors:
        result = [a * b for b in vector for a in result]
    return result


def digits(value: int, dimension: int, count: int) -> tuple[int, ...]:
    """Decompose a non-negative integer into ``count`` digits in a fixed radix, least significant digit first.

    Inverse of ``pack``; the part of ``value`` above ``count`` digits is
    discarded.

    Args:
        value: Non-negative integer to decompose.
        dimension: Radix base.
        count: Number of output digits.

    Returns:
        tuple[int, ...]: Tuple of digits with length ``count``.
    """
    result: list[int] = []
    for _ in range(count):
        value, bit = divmod(value, dimension)
        result.append(bit)
    return tuple(result)


def pack(values: Iterable[int], dimension: int) -> int:
    """Combine a least-significant-first digit sequence into an integer in a fixed radix; inverse of ``digits``.

    Args:
        values: Digit sequence; the i-th entry carries weight
            ``dimension**i``.
        dimension: Radix base.

    Returns:
        int: Weighted sum of the digits.
    """
    return sum(value * dimension**i for i, value in enumerate(values))


class Discretization:
    """Classical spatial discretization reference for a polynomial PDE on a grid.

    The single-component state is addressed with ``width`` bits: the low
    ``grid.spatial_width`` bits are the padded spatial coordinates and the high
    ``component_width`` bits select the field component; the row rules for port
    matrix elements and closure generators are both evaluated from difference
    stencils and cached.

    Args:
        pde: The ``PolynomialPDE`` to discretize; validated again at
            construction.
        grid: Grid whose spatial axes match the PDE.
        known: Known coefficient and forcing data, keyed by known field name
            with numeric sequences of length equal to the grid point count as
            values; must cover every known field referenced by the PDE.

    Attributes:
        component_width: Bit width of the field component index.
        width: Total bit width of a single-component state.
        dimension: Dimension of a single-component state, ``2**width``.
        spatial_storage: Storage length occupied by each field component,
            ``2**spatial_width``.
        ports: Mapping from port name to ``OperatorPort``.

    Raises:
        ValidationError: The PDE fails validation, the PDE and grid spatial
            axes differ, known data is missing, or a known data length differs
            from the grid point count.
    """
    def __init__(
        self,
        pde: PolynomialPDE,
        grid: Grid,
        known: Mapping[str, Sequence[complex]] | None = None,
    ) -> None:
        """Bind the PDE, grid, and known data, and initialize the derived widths and caches."""
        pde.validate()
        if tuple(pde.axes) != tuple(grid.axes):
            raise ValidationError("the PDE and the grid have different spatial axes")
        self.pde: PolynomialPDE = pde
        self.grid: Grid = grid
        self.component_width: int = (len(pde.fields) - 1).bit_length()
        self.width: int = grid.spatial_width + self.component_width
        self.dimension: int = 1 << self.width
        self.spatial_storage: int = 1 << grid.spatial_width
        self.known: dict[str, tuple[complex, ...]] = {
            k: tuple(map(complex, v)) for k, v in (known or {}).items()
        }
        required = {a.name for t in pde.terms for a in t.monomial.known}
        if not required <= self.known.keys():
            raise ValidationError(
                "missing known coefficient or forcing data: " + ", ".join(sorted(required - self.known.keys()))
            )
        if any(len(v) != grid.size for v in self.known.values()):
            raise ValidationError("known field data must cover the physical grid")
        self.ports: dict[str, OperatorPort] = {port.name: port for port in pde.ports}
        self._known_cache: dict[tuple[Atom, int], complex] = {}
        self._entry_cache: dict[tuple[str, int, int], complex] = {}
        self._row_entries_cache: dict[tuple[str, int], tuple[tuple[int, complex], ...]] = {}

    def known_value(self, atom: Atom, row: int) -> complex:
        """Evaluate a known-field atom at a physical grid point.

        The difference stencil of ``atom.derivative`` is applied to the known
        data for ``atom.name`` and summed with weights.

        Args:
            atom: Known-field atom, possibly carrying a spatial derivative
                specification.
            row: Physical grid point index.

        Returns:
            complex: Value of the known field or its derivative at the grid
            point, cached by ``(atom, row)``.
        """
        key = (atom, row)
        if key in self._known_cache:
            return self._known_cache[key]
        value = sum(
            weight * self.known[atom.name][col]
            for col, weight in self.grid.derivative_row(atom.derivative, row)
        )
        self._known_cache[key] = value
        return value

    def known_product(self, monomial: Monomial, row: int) -> complex:
        """Evaluate the product of the known part of a monomial at a physical grid point.

        Args:
            monomial: PDE monomial.
            row: Physical grid point index.

        Returns:
            complex: Running product of ``monomial.coefficient`` and the values
            of the known-field atoms, each carrying its own difference.
        """
        return monomial.coefficient * math.prod(
            self.known_value(atom, row) for atom in monomial.known
        )

    def entry(self, key: str, row: int, column: int) -> complex:
        """Evaluate a single matrix element of a multilinear port.

        The row is a single-component state index (component and spatial
        coordinates) and the column packs ``arity`` input state indices in a
        least-significant-first mixed radix; after matching output and input
        components term by term, the contributions of the outer derivative,
        known coefficients, and each field derivative are computed with
        difference stencils.

        Args:
            key: Port name, e.g. ``L``, ``F``, or ``B_i``.
            row: Output row index.
            column: Packed input column index.

        Returns:
            complex: The matrix element; 0 when the row or column is out of
            range. The result is cached by ``(key, row, column)``.
        """
        cache_key = (key, row, column)
        if cache_key in self._entry_cache:
            return self._entry_cache[cache_key]
        port = self.ports[key]
        if not 0 <= row < self.dimension or not 0 <= column < self.dimension**port.arity:
            return 0j
        component, x = divmod(row, self.spatial_storage)
        coordinates = digits(column, self.dimension, port.arity)
        result = 0j
        for term in port.terms:
            monomial = term.monomial
            if component != self.pde.fields.index(term.output):
                continue
            local: list[int] = []
            for atom, index in zip(monomial.fields, coordinates, strict=True):
                field, coordinate = divmod(index, self.spatial_storage)
                if field != self.pde.fields.index(atom.name):
                    break
                local.append(coordinate)
            else:
                for center, outer in self.grid.derivative_row(monomial.outer_derivative, x):
                    value = outer * self.known_product(monomial, center)
                    for atom, coordinate in zip(monomial.fields, local, strict=True):
                        value *= dict(self.grid.derivative_row(atom.derivative, center)).get(
                            coordinate, 0
                        )
                    result += value
        self._entry_cache[cache_key] = result
        return result

    def row_entries(self, key: str, row: int) -> tuple[tuple[int, complex], ...]:
        """Enumerate a row from the base difference stencils without scanning N^arity input coordinates.

        Args:
            key: Port name, e.g. ``L``, ``F``, or ``B_i``.
            row: Output row index.

        Returns:
            tuple: ``(packed column index, weight)`` pairs sorted by column,
            zero weights dropped; an empty tuple when ``row`` is out of range.
            The result is cached by ``(key, row)``.
        """
        import itertools

        cache_key = (key, row)
        if cache_key in self._row_entries_cache:
            return self._row_entries_cache[cache_key]
        if not 0 <= row < self.dimension:
            return ()
        component, x = divmod(row, self.spatial_storage)
        result: dict[int, complex] = {}
        for term in self.ports[key].terms:
            if component != self.pde.fields.index(term.output):
                continue
            monomial = term.monomial
            for center, outer in self.grid.derivative_row(monomial.outer_derivative, x):
                prefactor = outer * self.known_product(monomial, center)
                rows: list[tuple[tuple[int, float], ...]] = []
                for atom in monomial.fields:
                    field = self.pde.fields.index(atom.name) * self.spatial_storage
                    rows.append(
                        tuple(
                            (field + col, weight)
                            for col, weight in self.grid.derivative_row(atom.derivative, center)
                        )
                    )
                for choice in itertools.product(*rows):
                    column = pack((col for col, _ in choice), self.dimension)
                    value = prefactor * math.prod(weight for _, weight in choice)
                    result[column] = result.get(column, 0) + value
        self._row_entries_cache[cache_key] = tuple(
            (column, value) for column, value in sorted(result.items()) if value
        )
        return self._row_entries_cache[cache_key]

    def qcl_row(self, plan: QHAMPlan, eta: complex, row: int) -> tuple[tuple[int, complex], ...]:
        """Enumerate one row of the QCL closure generator under global lifted indices.

        ``plan.locate`` first locates the owning block of the row and its local
        index, then each linear edge embeds the port ``row_entries`` into the
        source block coordinates and multiplies the homotopy weight evaluated
        at ``eta``.

        Args:
            plan: ``QHAMPlan`` closure plan.
            eta: Homotopy parameter.
            row: Global row index within the closure system.

        Returns:
            tuple: ``(global column index, weight)`` pairs sorted by column,
            zero weights dropped; an empty tuple when ``row`` is out of range.
        """
        if not 0 <= row < plan.raw_dimension(self.dimension):
            return ()
        block, local = plan.locate(row, self.dimension)
        coordinates = digits(local, self.dimension, block.rank)
        result: dict[int, complex] = {}
        for term in plan.row_terms(block):
            weight = term.weight.evaluate(eta)
            for column, value in self.row_entries(term.operator, coordinates[term.position]):
                source = (
                    coordinates[: term.position]
                    + digits(column, self.dimension, term.arity)
                    + coordinates[term.position + 1 :]
                )
                index = plan.offset(term.column, self.dimension) + pack(source, self.dimension)
                result[index] = result.get(index, 0) + weight * value
        return tuple((column, value) for column, value in sorted(result.items()) if value)

    def qcl_entry(self, plan: QHAMPlan, eta: complex, row: int, column: int) -> complex:
        """Classical entry reference usable for review; does not masquerade as a reversible quantum oracle.

        Args:
            plan: ``QHAMPlan`` closure plan.
            eta: Homotopy parameter.
            row: Global row index within the closure system.
            column: Global column index within the closure system.

        Returns:
            complex: Matrix element of the closure generator at the given row
            and column; 0 when either index is out of range.
        """
        return dict(self.qcl_row(plan, eta, row)).get(column, 0j)

    def apply_port(self, key: str, vectors: Sequence[Sequence[complex]]) -> list[complex]:
        """Differentiate and multiply the fields directly, for independent validation of the port and closure rules.

        Args:
            key: Port name, e.g. ``L``, ``F``, or ``B_i``.
            vectors: Sequence of input state vectors; its length must equal the
                port arity and each vector must have length ``dimension``.

        Returns:
            list[complex]: Multilinear evaluation of the port with length
            ``dimension``.
        """
        port = self.ports[key]
        if len(vectors) != port.arity or any(len(v) != self.dimension for v in vectors):
            raise ValidationError("multilinear port input count or dimension mismatch")
        result = [0j] * self.dimension
        for term in port.terms:
            m = term.monomial
            values: list[complex] = []
            for row in range(self.grid.size):
                value = self.known_product(m, row)
                for atom, vector in zip(m.fields, vectors, strict=True):
                    offset = self.pde.fields.index(atom.name) * self.spatial_storage
                    value *= sum(
                        weight * vector[offset + col]
                        for col, weight in self.grid.derivative_row(atom.derivative, row)
                    )
                values.append(value)
            offset = self.pde.fields.index(term.output) * self.spatial_storage
            for row in range(self.grid.size):
                result[offset + row] += sum(
                    weight * values[col]
                    for col, weight in self.grid.derivative_row(m.outer_derivative, row)
                )
        return result

    def ham_rhs(self, values: Sequence[Sequence[complex]], eta: complex) -> list[list[complex]]:
        """Compute the classical right-hand sides of the order components ``Ui'`` by the HAM recursion.

        ``U0' = L U0 + f``; ``Ui' = L Ui - eta*sum_l (1+eta)^(i-1-l) C_l``,
        where ``C_l`` is the multilinear evaluation of the nonlinear ports on
        the order split ``l``.

        Args:
            values: Order components ``[U0, ..., Um]``, each a vector of length
                ``dimension``.
            eta: Homotopy parameter.

        Returns:
            list: List of vectors with the same length as ``values``; the i-th
            entry is ``Ui'``.

        Raises:
            ValidationError: Some vector length differs from ``dimension``.
        """
        size = self.dimension

        def linear(v: Sequence[complex]) -> list[complex]:
            """Action of the linear port ``L`` on a single component; returns the zero vector when no ``L`` port exists."""
            return self.apply_port("L", [v]) if "L" in self.ports else [0j] * size

        forcing = self.apply_port("F", []) if "F" in self.ports else [0j] * size
        result = [[a + b for a, b in zip(linear(values[0]), forcing, strict=True)]]
        previous = [0j] * size
        for order in range(1, len(values)):
            nonlinear = [0j] * size
            for port in self.pde.ports:
                if port.arity < 2:
                    continue
                for indices in compositions(order - 1, port.arity):
                    image = self.apply_port(port.name, [values[i] for i in indices])
                    nonlinear = [a + b for a, b in zip(nonlinear, image, strict=True)]
            previous = [(1 + eta) * a - eta * b for a, b in zip(previous, nonlinear, strict=True)]
            result.append([a + b for a, b in zip(linear(values[order]), previous, strict=True)])
        return result

    def lift(self, plan: QHAMPlan, values: Sequence[Sequence[complex]]) -> list[complex]:
        """Lift the order components into the block layout of the closure system.

        The physical block holds the sum of the order components; tensor blocks
        hold the tensor products of their factors via ``tensor_values``; the
        empty-word constant block is identically 1.

        Args:
            plan: ``QHAMPlan`` closure plan.
            values: Order-component vectors ``[U0, ..., Um]``.

        Returns:
            list[complex]: Lifted vector with length
            ``plan.raw_dimension(dimension)``.
        """
        result = [0j] * plan.raw_dimension(self.dimension)
        for block in plan.blocks():
            if block.kind == "physical":
                data = [sum(v[i] for v in values) for i in range(self.dimension)]
            else:
                data = tensor_values([values[i] for i in block.orders])
            offset = plan.offset(block, self.dimension)
            result[offset : offset + len(data)] = data
        return result

    def chain_rule(
        self, plan: QHAMPlan, values: Sequence[Sequence[complex]], eta: complex
    ) -> list[complex]:
        """Compute the time derivative of the lifted vector directly with the product rule.

        Each tensor-factor position of the ``lift`` output is differentiated
        and the order-wise ``Ui'`` from ``ham_rhs`` substituted; the result
        should agree component by component with the generator matrix acting on
        the lifted vector, forming an independent mathematical witness.

        Args:
            plan: ``QHAMPlan`` closure plan.
            values: Order-component vectors ``[U0, ..., Um]``.
            eta: Homotopy parameter.

        Returns:
            list[complex]: Derivative vector with length
            ``plan.raw_dimension(dimension)``.

        Raises:
            ValidationError: Some vector length differs from ``dimension``.
        """
        rhs = self.ham_rhs(values, eta)
        result = [0j] * plan.raw_dimension(self.dimension)
        for block in plan.blocks():
            if block.kind == "physical":
                data = [sum(v[i] for v in rhs) for i in range(self.dimension)]
            else:
                data = [0j] * self.dimension**block.rank
                for pos, _index in enumerate(block.orders):
                    factors = [
                        rhs[a] if k == pos else values[a] for k, a in enumerate(block.orders)
                    ]
                    image = tensor_values(factors)
                    data = [a + b for a, b in zip(data, image, strict=True)]
            offset = plan.offset(block, self.dimension)
            result[offset : offset + len(data)] = data
        return result

    def linear_action(
        self, plan: QHAMPlan, eta: complex, values: Sequence[complex]
    ) -> list[complex]:
        """Act block row by block row, avoiding storage of raw_dimension^2 matrix elements.

        Args:
            plan: ``QHAMPlan`` closure plan.
            eta: Homotopy parameter.
            values: Lifted vector of length ``plan.raw_dimension(dimension)``.

        Returns:
            list[complex]: Vector after the closure generator acts, with the
            same length as the input vector.
        """
        size = plan.raw_dimension(self.dimension)
        if len(values) != size:
            raise ValidationError("wrong lifted vector dimension")
        result = [0j] * size
        for block in plan.blocks():
            offset = plan.offset(block, self.dimension)
            for term in plan.row_terms(block):
                coefficient = term.weight.evaluate(eta)
                if not coefficient:
                    continue
                source = plan.offset(term.column, self.dimension)
                for local_row in range(self.dimension**block.rank):
                    coordinates = digits(local_row, self.dimension, block.rank)
                    for local_column, value in self.row_entries(
                        term.operator, coordinates[term.position]
                    ):
                        inputs = (
                            coordinates[: term.position]
                            + digits(local_column, self.dimension, term.arity)
                            + coordinates[term.position + 1 :]
                        )
                        result[offset + local_row] += (
                            coefficient * value * values[source + pack(inputs, self.dimension)]
                        )
        return result

    def matrix(
        self, plan: QHAMPlan, eta: complex, *, max_dimension: int = 512
    ) -> list[list[complex]]:
        """Materialize the full square matrix of the closure generator.

        Args:
            plan: ``QHAMPlan`` closure plan.
            eta: Homotopy parameter.
            max_dimension: Maximum total dimension allowed for
                materialization.

        Returns:
            list: Square complex matrix of order
            ``plan.raw_dimension(dimension)``, organized by row.

        Raises:
            ValidationError: The total dimension exceeds ``max_dimension``.
        """
        size = plan.raw_dimension(self.dimension)
        if size > max_dimension:
            raise ValidationError("materializing G is only allowed for small-scale mathematical witnesses; use the lazy row rules instead")
        result = [[0j] * size for _ in range(size)]
        for block in plan.blocks():
            offset = plan.offset(block, self.dimension)
            for term in plan.row_terms(block):
                coefficient = term.weight.evaluate(eta)
                if coefficient == 0:
                    continue
                source = plan.offset(term.column, self.dimension)
                for local_row in range(self.dimension**block.rank):
                    row_coordinates = digits(local_row, self.dimension, block.rank)
                    for local_column in range(self.dimension**term.arity):
                        value = coefficient * self.entry(
                            term.operator, row_coordinates[term.position], local_column
                        )
                        if not value:
                            continue
                        source_coordinates = (
                            row_coordinates[: term.position]
                            + digits(local_column, self.dimension, term.arity)
                            + row_coordinates[term.position + 1 :]
                        )
                        column = source + pack(source_coordinates, self.dimension)
                        result[offset + local_row][column] += value
        return result

    def encode_fields(self, values: Mapping[str, Sequence[complex]]) -> tuple[complex, ...]:
        """Pack initial data given per field component into a single-component state vector.

        Args:
            values: Mapping keyed by field component name with data sequences
                of length equal to the grid point count as values; must cover
                all components.

        Returns:
            tuple[complex]: Vector of length ``dimension`` whose c-th component
            occupies the coordinates ``[c*spatial_storage,
            c*spatial_storage+size)``.

        Raises:
            ValidationError: The component set is incomplete or a data length
                mismatches.
        """
        if set(values) != set(self.pde.fields):
            raise ValidationError("initial fields must cover all components")
        output = [0j] * self.dimension
        for component, name in enumerate(self.pde.fields):
            if len(values[name]) != self.grid.size:
                raise ValidationError("initial field length mismatch")
            output[
                component * self.spatial_storage : component * self.spatial_storage + self.grid.size
            ] = map(complex, values[name])
        return tuple(output)
