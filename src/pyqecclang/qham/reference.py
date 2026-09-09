"""小规模数学见证：有限差分端口与独立 HAM/链式法则求值，不是量子模拟。"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction
from functools import lru_cache

from ..ir import ValidationError
from .linearization import compositions


@lru_cache(maxsize=32)
def centered_weights(order):
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
    axes: tuple[str, ...]
    shape: tuple[int, ...]
    spacing: tuple[float, ...]
    boundary: str = "periodic"
    _row_cache: dict = field(
        default_factory=dict, init=False, compare=False, repr=False, hash=False
    )

    def __post_init__(self):
        if (
            not len(self.axes) == len(self.shape) == len(self.spacing)
            or any(type(n) is not int or n < 1 for n in self.shape)
            or any(x <= 0 or not math.isfinite(x) for x in self.spacing)
        ):
            raise ValidationError("空间网格形状/间距无效")
        if self.boundary not in {"periodic", "dirichlet_zero"}:
            raise ValidationError("目前提供周期或零延拓 Dirichlet 差分")

    @property
    def size(self):
        return math.prod(self.shape)

    @property
    def spatial_width(self):
        return max(1, (self.size - 1).bit_length())

    def coordinates(self, index):
        result = []
        for n in self.shape:
            index, value = divmod(index, n)
            result.append(value)
        return result

    def address(self, coords):
        stride, result = 1, 0
        for n, x in zip(self.shape, coords, strict=True):
            result += stride * x
            stride *= n
        return result

    def derivative_row(self, derivative, row):
        key = (derivative, row)
        if key in self._row_cache:
            return self._row_cache[key]
        if not 0 <= row < self.size:
            return ()
        result = {row: 1.0}
        for axis, order in derivative:
            position = self.axes.index(axis)
            step = self.spacing[position]
            following = {}
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


def tensor_values(vectors):
    result = [1.0 + 0j]
    for vector in vectors:
        result = [a * b for b in vector for a in result]
    return result


def digits(value, dimension, count):
    result = []
    for _ in range(count):
        value, bit = divmod(value, dimension)
        result.append(bit)
    return tuple(result)


def pack(values, dimension):
    return sum(value * dimension**i for i, value in enumerate(values))


class Discretization:
    def __init__(self, pde, grid, known=None):
        pde.validate()
        if tuple(pde.axes) != tuple(grid.axes):
            raise ValidationError("PDE 与网格空间轴不同")
        self.pde, self.grid = pde, grid
        self.component_width = (len(pde.fields) - 1).bit_length()
        self.width = grid.spatial_width + self.component_width
        self.dimension = 1 << self.width
        self.spatial_storage = 1 << grid.spatial_width
        self.known = {k: tuple(map(complex, v)) for k, v in (known or {}).items()}
        required = {a.name for t in pde.terms for a in t.monomial.known}
        if not required <= self.known.keys():
            raise ValidationError(
                "缺少已知系数/强迫数据：" + ", ".join(sorted(required - self.known.keys()))
            )
        if any(len(v) != grid.size for v in self.known.values()):
            raise ValidationError("已知场数据必须覆盖物理网格")
        self.ports = {port.name: port for port in pde.ports}
        self._known_cache, self._entry_cache, self._row_entries_cache = {}, {}, {}

    def known_value(self, atom, row):
        key = (atom, row)
        if key in self._known_cache:
            return self._known_cache[key]
        value = sum(
            weight * self.known[atom.name][col]
            for col, weight in self.grid.derivative_row(atom.derivative, row)
        )
        self._known_cache[key] = value
        return value

    def known_product(self, monomial, row):
        return monomial.coefficient * math.prod(
            self.known_value(atom, row) for atom in monomial.known
        )

    def entry(self, key, row, column):
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
            local = []
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

    def row_entries(self, key, row):
        """利用基础差分模板枚举行，不扫描 N^arity 个输入坐标。"""
        import itertools

        cache_key = (key, row)
        if cache_key in self._row_entries_cache:
            return self._row_entries_cache[cache_key]
        if not 0 <= row < self.dimension:
            return ()
        component, x = divmod(row, self.spatial_storage)
        result = {}
        for term in self.ports[key].terms:
            if component != self.pde.fields.index(term.output):
                continue
            monomial = term.monomial
            for center, outer in self.grid.derivative_row(monomial.outer_derivative, x):
                prefactor = outer * self.known_product(monomial, center)
                rows = []
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

    def qcl_row(self, plan, eta, row):
        if not 0 <= row < plan.raw_dimension(self.dimension):
            return ()
        block, local = plan.locate(row, self.dimension)
        coordinates = digits(local, self.dimension, block.rank)
        result = {}
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

    def qcl_entry(self, plan, eta, row, column):
        """可用于审阅的经典 entry reference，不冒充可逆量子 oracle。"""
        return dict(self.qcl_row(plan, eta, row)).get(column, 0j)

    def apply_port(self, key, vectors):
        """直接对场求导并相乘，用于独立验证矩阵端口/闭包规则。"""
        port = self.ports[key]
        if len(vectors) != port.arity or any(len(v) != self.dimension for v in vectors):
            raise ValidationError("多线性端口的输入数目/维度不符")
        result = [0j] * self.dimension
        for term in port.terms:
            m = term.monomial
            values = []
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

    def ham_rhs(self, values, eta):
        size = self.dimension

        def linear(v):
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

    def lift(self, plan, values):
        result = [0j] * plan.raw_dimension(self.dimension)
        for block in plan.blocks():
            if block.kind == "physical":
                data = [sum(v[i] for v in values) for i in range(self.dimension)]
            else:
                data = tensor_values([values[i] for i in block.orders])
            offset = plan.offset(block, self.dimension)
            result[offset : offset + len(data)] = data
        return result

    def chain_rule(self, plan, values, eta):
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

    def linear_action(self, plan, eta, values):
        """逐个块行作用，避免存储 raw_dimension^2 个矩阵元。"""
        size = plan.raw_dimension(self.dimension)
        if len(values) != size:
            raise ValidationError("提升向量维数错误")
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

    def matrix(self, plan, eta, *, max_dimension=512):
        size = plan.raw_dimension(self.dimension)
        if size > max_dimension:
            raise ValidationError("仅小规模数学见证允许物化 G；请使用惰性行规则")
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

    def encode_fields(self, values):
        if set(values) != set(self.pde.fields):
            raise ValidationError("初始场必须覆盖全部分量")
        output = [0j] * self.dimension
        for component, name in enumerate(self.pde.fields):
            if len(values[name]) != self.grid.size:
                raise ValidationError("初始场长度不匹配")
            output[
                component * self.spatial_storage : component * self.spatial_storage + self.grid.size
            ] = map(complex, values[name])
        return tuple(output)
