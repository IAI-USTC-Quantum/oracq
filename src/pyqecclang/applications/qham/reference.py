"小规模数学见证：有限差分端口与独立 HAM/链式法则求值，不是量子模拟。"

from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction
from functools import lru_cache

from pyqecclang.applications.qham.linearization import compositions
from pyqecclang.infrastructure.ir import ValidationError


@lru_cache(maxsize=32)
def centered_weights(order):
    """计算单位步长一维中心差分模板的权重。

    在以原点为中心的整数格点上用精确有理数消元求解导数插值条件，
    对次数小于格点数的多项式精确；零权重被剔除。

    Args:
        order: 导数阶数，非负整数。

    Returns:
        tuple: ``(格点偏移, 权重)`` 对，权重为 float，对应单位步长；
        使用时还需除以 ``spacing**order``。结果按阶数缓存。
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
    """空间离散网格：轴、形状、间距与边界类型，并给出索引和差分模板。

    线性格点编号采用混合基布局，轴 0 为最低位；``coordinates`` 与
    ``address`` 互为逆映射。

    Attributes:
        axes: 空间轴名元组，须与 PDE 的轴一致。
        shape: 各轴格点数，均为正整数。
        spacing: 各轴格点间距，均为正有限值。
        boundary: 边界类型：``periodic`` 把越界坐标按轴长回绕；
            ``dirichlet_zero`` 按零延拓丢弃越界模板点。

    Raises:
        ValidationError: 轴数不一致、形状/间距非法或边界类型不受支持。
    """
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
        """网格总格点数，即 ``shape`` 各分量之积。"""
        return math.prod(self.shape)

    @property
    def spatial_width(self):
        """索引空间坐标所需的位宽，至少为 1。

        空间存储按 ``1 << spatial_width`` 对齐，可能大于 ``size``。"""
        return max(1, (self.size - 1).bit_length())

    def coordinates(self, index):
        """把线性格点编号分解为各轴坐标，为 ``address`` 的逆。

        Args:
            index: 线性格点编号。

        Returns:
            list[int]: 各轴坐标，轴 0 为最低位，与 ``shape`` 同长。
        """
        result = []
        for n in self.shape:
            index, value = divmod(index, n)
            result.append(value)
        return result

    def address(self, coords):
        """把各轴坐标合成为线性格点编号。

        Args:
            coords: 各轴坐标序列，长度须与 ``shape`` 一致。

        Returns:
            int: 混合基布局下的线性格点编号。
        """
        stride, result = 1, 0
        for n, x in zip(self.shape, coords, strict=True):
            result += stride * x
            stride *= n
        return result

    def derivative_row(self, derivative, row):
        """计算混合空间导数在给定格点处的稀疏差分行。

        逐轴组合一维中心差分模板，权重已除以相应 ``spacing`` 的幂次；
        周期边界把越界格点回绕，``dirichlet_zero`` 丢弃越界格点。

        Args:
            derivative: ``(轴名, 阶数)`` 对组成的导数说明。
            row: 线性格点编号。

        Returns:
            tuple: 按列地址排序的 ``(列地址, 权重)`` 对，零权重被剔除；
            ``row`` 越界时为空元组。结果按 ``(derivative, row)`` 缓存。

        Raises:
            ValueError: ``derivative`` 引用了不在 ``axes`` 中的轴。
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
    """计算多个向量的张量积并展平为一维复数列表。

    因子 0 占据最低位：第 k 个因子的下标跨度为其前各因子长度之积；
    无因子时返回单元素 ``[1+0j]``。

    Args:
        vectors: 因子向量序列。

    Returns:
        list[complex]: 展平后的张量积，长度为各因子长度之积。
    """
    result = [1.0 + 0j]
    for vector in vectors:
        result = [a * b for b in vector for a in result]
    return result


def digits(value, dimension, count):
    """把非负整数按固定进制分解为 ``count`` 个数位，最低位在前。

    与 ``pack`` 互为逆运算；``value`` 高于 ``count`` 位的部分被丢弃。

    Args:
        value: 待分解的非负整数。
        dimension: 进制基数。
        count: 输出数位个数。

    Returns:
        tuple[int, ...]: 长度为 ``count`` 的数位元组。
    """
    result = []
    for _ in range(count):
        value, bit = divmod(value, dimension)
        result.append(bit)
    return tuple(result)


def pack(values, dimension):
    """把低位在前的数位序列按固定进制合成为整数，为 ``digits`` 的逆。

    Args:
        values: 数位序列，第 i 项的权重为 ``dimension**i``。
        dimension: 进制基数。

    Returns:
        int: 各数位的加权和。
    """
    return sum(value * dimension**i for i, value in enumerate(values))


class Discretization:
    """多项式 PDE 在网格上的经典空间离散参考。

    单分量状态按 ``width`` 位编址：低 ``grid.spatial_width`` 位是补齐的
    空间坐标，高 ``component_width`` 位选择场分量；端口矩阵元与闭包生成
    元的行规则均基于差分模板求值并缓存。

    Args:
        pde: 待离散的 ``PolynomialPDE``，构造时会再次验证。
        grid: 与 PDE 空间轴一致的网格。
        known: 已知系数/强迫数据，键为已知场名，值为长度等于网格格点数
            的数值序列，须覆盖 PDE 引用的全部已知场。

    Attributes:
        component_width: 场分量索引位宽。
        width: 单分量状态总位宽。
        dimension: 单分量状态维数 ``2**width``。
        spatial_storage: 每个场分量占用的存储长度 ``2**spatial_width``。
        ports: 端口名到 ``OperatorPort`` 的映射。

    Raises:
        ValidationError: PDE 未通过验证、PDE 与网格空间轴不同、缺少已知
            数据或已知数据长度不等于网格格点数。
    """
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
        """求已知场原子在物理格点处的取值。

        对 ``atom.name`` 对应的已知数据施加 ``atom.derivative`` 的差分
        模板并加权求和。

        Args:
            atom: 已知场原子，可带空间导数说明。
            row: 物理格点编号。

        Returns:
            complex: 该已知场（导数）在格点处的值，按 ``(atom, row)``
            缓存。
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

    def known_product(self, monomial, row):
        """求单项式已知部分在物理格点处的乘积。

        Args:
            monomial: PDE 单项式。
            row: 物理格点编号。

        Returns:
            complex: ``monomial.coefficient`` 与各已知场原子取值（各自
            带差分）的连乘积。
        """
        return monomial.coefficient * math.prod(
            self.known_value(atom, row) for atom in monomial.known
        )

    def entry(self, key, row, column):
        """求多线性端口的单个矩阵元。

        行是单分量状态索引（分量与空间坐标），列是 ``arity`` 个输入状态
        索引按低位在前的混合基打包；逐项匹配输出分量与输入分量后，用
        差分模板计算外导数、已知系数与各场导数的贡献。

        Args:
            key: 端口名，如 ``L``、``F`` 或 ``B_i``。
            row: 输出行索引。
            column: 打包后的输入列索引。

        Returns:
            complex: 矩阵元；行或列越界时为 0。结果按
            ``(key, row, column)`` 缓存。
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
        """枚举 QCL 闭包生成元在全局提升索引下的一行。

        先用 ``plan.locate`` 定位行所属块与块内局部索引，再逐条线性边把
        端口 ``row_entries`` 嵌入到源块坐标，并乘上同伦权重在 ``eta``
        处的取值。

        Args:
            plan: ``QHAMPlan`` 闭包计划。
            eta: 同伦参数。
            row: 闭包系统内的全局行索引。

        Returns:
            tuple: 按列排序的 ``(全局列索引, 权重)`` 对，零权重被剔除；
            ``row`` 越界时为空元组。
        """
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
        """按 HAM 递推计算各阶分量 ``Ui'`` 的经典右端。

        ``U0' = L U0 + f``；``Ui' = L Ui - eta*sum_l (1+eta)^(i-1-l) C_l``，
        其中 ``C_l`` 是非线性端口在阶数拆分 ``l`` 上的多线性求值。

        Args:
            values: 各阶分量 ``[U0, ..., Um]``，每个均为长度
                ``dimension`` 的向量。
            eta: 同伦参数。

        Returns:
            list: 与 ``values`` 等长的向量列表，第 i 项为 ``Ui'``。

        Raises:
            ValidationError: 某个向量长度不等于 ``dimension``。
        """
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
        """把各阶分量提升到闭包系统的块布局。

        物理块放各阶分量之和；张量块放各因子按 ``tensor_values`` 的
        张量积；空字常量块恒为 1。

        Args:
            plan: ``QHAMPlan`` 闭包计划。
            values: 各阶分量 ``[U0, ..., Um]`` 向量。

        Returns:
            list[complex]: 长度为 ``plan.raw_dimension(dimension)`` 的
            提升向量。
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

    def chain_rule(self, plan, values, eta):
        """用乘积法则直接计算提升向量对时间的导数。

        对 ``lift`` 输出的每个张量因子位置求导并代入 ``ham_rhs`` 得到的
        各阶 ``Ui'``；结果应与生成元矩阵作用于提升向量逐分量一致，构成
        独立的数学见证。

        Args:
            plan: ``QHAMPlan`` 闭包计划。
            values: 各阶分量 ``[U0, ..., Um]`` 向量。
            eta: 同伦参数。

        Returns:
            list[complex]: 长度为 ``plan.raw_dimension(dimension)`` 的
            导数向量。

        Raises:
            ValidationError: 某个向量长度不等于 ``dimension``。
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
        """物化闭包生成元的完整方阵。

        Args:
            plan: ``QHAMPlan`` 闭包计划。
            eta: 同伦参数。
            max_dimension: 允许物化的最大总维数。

        Returns:
            list: ``plan.raw_dimension(dimension)`` 阶的复数方阵，按行
            组织。

        Raises:
            ValidationError: 总维数超过 ``max_dimension``。
        """
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
        """把按场分量给出的初始数据打包成单分量状态向量。

        Args:
            values: 键为场分量名、值为长度等于网格格点数的数据序列，
                须覆盖全部分量。

        Returns:
            tuple[complex]: 长度为 ``dimension`` 的向量，第 c 个分量占据
            ``[c*spatial_storage, c*spatial_storage+size)`` 的坐标。

        Raises:
            ValidationError: 分量集合不全或数据长度不匹配。
        """
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
