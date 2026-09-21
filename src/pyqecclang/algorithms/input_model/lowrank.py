"""量子化学低秩分解（DF/THC）哈密顿量的 LCU/块编码组装。

实现 Berry et al. 2019、von Burg et al. 2021 的双因子分解（DF）与
Lee et al. 2021 的张量超收缩（THC）的"低秩张量 → LCU → BE"管道，
对应算法覆盖工作板的 C1 项。积分张量作为经典输入数据给出，本模块不做
真实的化学积分计算；显式小矩阵路径面向验证用的小实例，大规模时
``U_r``/``G_r`` 的谱与旋转角度表应改走 QRAM 数据绑定（参见 prepare_select
的 qram_prepare/alias_prepare 三层范式）。
"""

from __future__ import annotations

import cmath
import contextlib
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from pyqecclang.algorithms.input_model.block_encoding import adjoint_be, lcu, matrix_pauli_encoding
from pyqecclang.algorithms.input_model.contracts import finite_real, require_instance
from pyqecclang.algorithms.input_model.operators import BlockEncoding, _name, identity, product
from pyqecclang.algorithms.input_model.oracles import annotate, invoke, resources_for
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.ir import Bits, Ref, ValidationError, fuse

_SYNTHESIS_TOLERANCE = 1e-15
_MATRIX_LIMIT_QUBITS = 5


def _as_complex_matrix(
    value: Iterable[Iterable[complex]], path: str
) -> tuple[tuple[complex, ...], ...]:
    """把输入规范化为元素有限的复方阵元组，维度限制为 2..32 的二的幂。"""
    try:
        matrix = tuple(tuple(complex(v) for v in row) for row in value)
    except TypeError as exc:
        raise ValidationError(f"{path} 的元素必须是数值") from exc
    if not matrix or any(len(row) != len(matrix) for row in matrix):
        raise ValidationError(f"{path} 必须是非空方阵")
    d = len(matrix)
    if d < 2 or d & (d - 1) or d > (1 << _MATRIX_LIMIT_QUBITS):
        raise ValidationError(f"{path} 的维度必须是 2..32 内二的幂")
    if not all(math.isfinite(v.real) and math.isfinite(v.imag) for row in matrix for v in row):
        raise ValidationError(f"{path} 的元素必须有限")
    return matrix


def _matmul(
    a: Sequence[Sequence[complex]], b: Sequence[Sequence[complex]]
) -> tuple[tuple[complex, ...], ...]:
    """计算两个复方阵的矩阵乘积，返回行主序嵌套元组。"""
    return tuple(
        tuple(sum(ar[k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))) for ar in a
    )


def _check_unitary(
    matrix: Sequence[Sequence[complex]], path: str, *, tolerance: float = 1e-9
) -> None:
    """校验方阵的列正交归一性，不满足时抛 ``ValidationError``。"""
    d = len(matrix)
    for i in range(d):
        for j in range(d):
            value = sum(matrix[k][i].conjugate() * matrix[k][j] for k in range(d))
            if abs(value - (1 if i == j else 0)) > tolerance:
                raise ValidationError(f"{path} 必须是酉矩阵（列正交归一）")


def diagonalize_symmetric(
    matrix: Iterable[Iterable[float]], *, tolerance: float = 1e-12, max_sweeps: int = 100
) -> tuple[tuple[float, ...], tuple[tuple[float, ...], ...]]:
    """实对称矩阵的 Jacobi 特征分解：返回 (特征值, 特征向量矩阵)，满足 ``G = V diag(λ) Vᵀ``。

    这是 DF 输入模型的经典预处理：对称阵 ``G_r`` 对角化后，特征向量矩阵按
    ``from_symmetric`` 折叠进旋转 ``U_r``。特征向量矩阵的列是特征向量。

    Args:
        matrix: 实对称方阵，元素须为有限实数。
        tolerance: 对称性检查与收敛判定共用的容差。
        max_sweeps: Jacobi 扫描轮数上限，超限未收敛时报错。

    Returns:
        tuple[tuple[float, ...], tuple[tuple[float, ...], ...]]:
        ``(特征值元组, 列为特征向量的矩阵)``，满足 ``G = V diag(λ) Vᵀ``。
    """
    raw = tuple(tuple(row) for row in matrix)
    if not raw or any(len(row) != len(raw) for row in raw):
        raise ValidationError("diagonalize_symmetric.matrix 必须是非空方阵")
    d = len(raw)
    for row in raw:
        for value in row:
            finite_real(value, "diagonalize_symmetric.matrix")
    a = [list(map(float, row)) for row in raw]
    for p in range(d):
        for q in range(p + 1, d):
            if abs(a[p][q] - a[q][p]) > tolerance:
                raise ValidationError("diagonalize_symmetric.matrix 必须是实对称矩阵")
    v = [[1.0 if i == j else 0.0 for j in range(d)] for i in range(d)]
    for _ in range(max_sweeps):
        off = max((abs(a[p][q]) for p in range(d) for q in range(p + 1, d)), default=0.0)
        if off <= tolerance:
            break
        for p in range(d - 1):
            for q in range(p + 1, d):
                if abs(a[p][q]) <= tolerance:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = (1 if theta >= 0 else -1) / (abs(theta) + math.hypot(1.0, theta))
                c = 1.0 / math.sqrt(1.0 + t * t)
                s = t * c
                for k in range(d):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(d):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(d):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq
    else:
        raise ValidationError("Jacobi 特征分解未在 max_sweeps 内收敛")
    return tuple(a[k][k] for k in range(d)), tuple(tuple(row) for row in v)


@dataclass(frozen=True)
class DoubleFactorization:
    """DF 输入模型（CP input model）：``H = scalar·I + Σ_r U_r diag(g_r) U_r†`` 的小实例表示。

    rotations 是显式小酉矩阵，spectra 是对称阵 ``G_r`` 经经典对角化后的实特征值；
    物理 DF 中 ``G_r`` 由 ``from_symmetric`` 接收并折叠进旋转。显式矩阵路径仅面向
    验证用小实例；大规模时谱范数表与旋转角度走 QRAM 数据绑定。
    """

    scalar: float
    rotations: tuple[tuple[tuple[complex, ...], ...], ...]
    spectra: tuple[tuple[float, ...], ...]

    def __post_init__(self) -> None:
        """校验标量、各秩旋转的酉性与谱长度，并规范化存储格式。"""
        finite_real(self.scalar, "DoubleFactorization.scalar")
        rotations = tuple(
            _as_complex_matrix(r, f"DoubleFactorization.rotations[{i}]")
            for i, r in enumerate(self.rotations)
        )
        if not rotations:
            raise ValidationError("DoubleFactorization 至少需要一个秩项")
        if len(rotations) != len(self.spectra):
            raise ValidationError("DoubleFactorization 的 rotations 与 spectra 数量不一致")
        d = len(rotations[0])
        if any(len(r) != d for r in rotations):
            raise ValidationError("DoubleFactorization 的各秩旋转维度不一致")
        for i, rotation in enumerate(rotations):
            _check_unitary(rotation, f"DoubleFactorization.rotations[{i}]")
        spectra: list[tuple[float, ...]] = []
        for i, spectrum in enumerate(self.spectra):
            values = tuple(spectrum)
            if len(values) != d:
                raise ValidationError(f"DoubleFactorization.spectra[{i}] 的长度必须等于维度 {d}")
            for value in values:
                finite_real(value, f"DoubleFactorization.spectra[{i}]")
            spectra.append(tuple(float(v) for v in values))
        object.__setattr__(self, "scalar", float(self.scalar))
        object.__setattr__(self, "rotations", rotations)
        object.__setattr__(self, "spectra", tuple(spectra))

    @classmethod
    def from_symmetric(
        cls,
        scalar: float,
        rotations: Iterable[Iterable[Iterable[complex]]],
        factors: Iterable[Iterable[Iterable[float]]],
    ) -> DoubleFactorization:
        """从显式酉 ``U_r`` 与实对称 ``G_r`` 构造：``G_r = V_r diag(g_r) V_rᵀ`` 折叠为 ``U_r V_r``。

        Args:
            scalar: 恒等项系数。
            rotations: 各秩的显式酉矩阵 ``U_r``。
            factors: 与旋转一一对应、同维的实对称阵 ``G_r``。

        Returns:
            DoubleFactorization: 折叠 ``U_r V_r`` 后谱不变的 DF 输入模型。
        """
        rotations = tuple(rotations)
        factors = tuple(factors)
        if len(rotations) != len(factors):
            raise ValidationError("from_symmetric 的 rotations 与 factors 数量不一致")
        combined: list[tuple[tuple[complex, ...], ...]] = []
        spectra: list[tuple[float, ...]] = []
        for rotation, factor in zip(rotations, factors, strict=True):
            unitary = _as_complex_matrix(rotation, "from_symmetric.rotations")
            eigenvalues, vectors = diagonalize_symmetric(factor)
            if len(vectors) != len(unitary):
                raise ValidationError("from_symmetric 的旋转与对称阵维度不一致")
            combined.append(_matmul(unitary, vectors))
            spectra.append(eigenvalues)
        return cls(scalar, tuple(combined), tuple(spectra))

    @property
    def width(self) -> int:
        """目标量子位数。"""
        return (len(self.rotations[0]) - 1).bit_length()

    @property
    def rank(self) -> int:
        """DF 秩项数。"""
        return len(self.rotations)


@dataclass(frozen=True)
class THCDecomposition:
    """THC 输入模型：``H = Σ_{μν} ζ_{μν} L_μ L_ν†``，``L_μ`` 为显式小矩阵叶算符。

    coefficients 是实对称的 ``ζ`` 矩阵，leaves 是同维显式小矩阵；叶算符无需酉或
    Hermitian，乘积 ``L_μ L_ν†`` 的块编码由 matrix_pauli_encoding 与 BE 乘积组装。
    """

    coefficients: tuple[tuple[float, ...], ...]
    leaves: tuple[tuple[tuple[complex, ...], ...], ...]

    def __post_init__(self) -> None:
        """校验 ζ 矩阵的实对称性、维度一致性与叶算符格式。"""
        zeta = tuple(tuple(row) for row in self.coefficients)
        if not zeta or any(len(row) != len(zeta) for row in zeta):
            raise ValidationError("THCDecomposition.coefficients 必须是非空方阵")
        for row in zeta:
            for value in row:
                finite_real(value, "THCDecomposition.coefficients")
        size = len(zeta)
        for mu in range(size):
            for nu in range(mu + 1, size):
                if abs(zeta[mu][nu] - zeta[nu][mu]) > 1e-12:
                    raise ValidationError("THCDecomposition.coefficients 必须是实对称矩阵")
        if len(self.leaves) != size:
            raise ValidationError("THCDecomposition 的 leaves 数量必须与 ζ 的维度一致")
        leaves = tuple(
            _as_complex_matrix(leaf, f"THCDecomposition.leaves[{i}]")
            for i, leaf in enumerate(self.leaves)
        )
        if len({len(leaf) for leaf in leaves}) != 1:
            raise ValidationError("THCDecomposition 的叶算符维度不一致")
        object.__setattr__(self, "coefficients", zeta)
        object.__setattr__(self, "leaves", leaves)

    @property
    def width(self) -> int:
        """目标量子位数。"""
        return (len(self.leaves[0]) - 1).bit_length()

    @property
    def leaf_count(self) -> int:
        """THC 叶算符个数。"""
        return len(self.leaves)


def _emit_transposition(builder: Builder, ref: Ref, first: int, second: int) -> None:
    """沿 Gray 路径用多控 X 交换 ``|first>`` 与 ``|second>``，其余基态不动。"""
    path = [first]
    for bit in range(ref.width):
        if ((first ^ second) >> bit) & 1:
            path.append(path[-1] ^ (1 << bit))
    edges = list(zip(path, path[1:], strict=False))
    for left, _right in edges + list(reversed(edges[:-1])):
        bit = (left ^ _right).bit_length() - 1
        controls = fuse(ref[:bit], ref[bit + 1 :])
        value = (left & ((1 << bit) - 1)) | ((left >> (bit + 1)) << bit)
        with builder.control(controls, value) if controls.width else contextlib.nullcontext():
            builder.x(ref[bit])


def _zyz(matrix: Sequence[Sequence[complex]]) -> tuple[float, float, float, float]:
    """二阶酉的 ``e^{iφ} Rz(α) Ry(β) Rz(γ)`` 精确分解（模拟器约定见 execution.gate_matrix）。"""
    a00, a01, a10, a11 = matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1]
    beta = 2 * math.atan2(abs(a10), abs(a00))
    if abs(a10) <= _SYNTHESIS_TOLERANCE:
        phi = (cmath.phase(a00) + cmath.phase(a11)) / 2
        return phi, cmath.phase(a11) - cmath.phase(a00), 0.0, 0.0
    if abs(a00) <= _SYNTHESIS_TOLERANCE:
        lower, upper = cmath.phase(a10), cmath.phase(-a01)
        return (lower + upper) / 2, lower - upper, beta, 0.0
    phi = (cmath.phase(a00) + cmath.phase(a11)) / 2
    alpha = cmath.phase(a11) + cmath.phase(a10) - 2 * phi
    gamma = cmath.phase(a11) - cmath.phase(a10)
    return phi, alpha, beta, gamma


def _emit_controlled_two_by_two(
    builder: Builder, ref: Ref, bit: int, anchor: int, matrix: Sequence[Sequence[complex]]
) -> None:
    """在其余位固定为 ``anchor`` 的子空间上施加 ``matrix``（基序为 bit=0,1）。"""
    phi, alpha, beta, gamma = _zyz(matrix)
    controls = fuse(ref[:bit], ref[bit + 1 :])
    value = (anchor & ((1 << bit) - 1)) | ((anchor >> (bit + 1)) << bit)
    with builder.control(controls, value) if controls.width else contextlib.nullcontext():
        if abs(gamma) > _SYNTHESIS_TOLERANCE:
            builder.rz(ref[bit], gamma)
        if abs(beta) > _SYNTHESIS_TOLERANCE:
            builder.ry(ref[bit], beta)
        if abs(alpha) > _SYNTHESIS_TOLERANCE:
            builder.rz(ref[bit], alpha)
        if abs(phi) > _SYNTHESIS_TOLERANCE:
            builder.global_phase(phi)


def _emit_two_level(
    builder: Builder, ref: Ref, p: int, q: int, matrix: Sequence[Sequence[complex]]
) -> None:
    """施加只作用于 ``span{|p>, |q>}`` 的两能级酉，基序为 (``|p>``, ``|q>``)。"""
    bit = ((p ^ q) & -(p ^ q)).bit_length() - 1
    moved = p ^ (1 << bit)
    if moved != q:
        _emit_transposition(builder, ref, q, moved)
    if (p >> bit) & 1:
        swapped = ((matrix[1][1], matrix[1][0]), (matrix[0][1], matrix[0][0]))
        _emit_controlled_two_by_two(builder, ref, bit, p, swapped)
    else:
        _emit_controlled_two_by_two(builder, ref, bit, p, matrix)
    if moved != q:
        _emit_transposition(builder, ref, q, moved)


def _emit_unitary(builder: Builder, ref: Ref, matrix: Sequence[Sequence[complex]]) -> None:
    """两能级分解合成显式小酉矩阵：逐列消元为对角相位后按逆序回放。"""
    d = len(matrix)
    u = [list(row) for row in matrix]
    steps: list[tuple[int, int, tuple[tuple[complex, complex], tuple[complex, complex]]]] = []
    for column in range(d - 1):
        for row in range(d - 1, column, -1):
            value = u[row][column]
            if abs(value) <= _SYNTHESIS_TOLERANCE:
                continue
            pivot = row - 1
            head = u[pivot][column]
            norm = math.sqrt(abs(head) ** 2 + abs(value) ** 2)
            block = (
                (head.conjugate() / norm, value.conjugate() / norm),
                (-value / norm, head / norm),
            )
            steps.append((pivot, row, block))
            for k in range(d):
                top, bottom = u[pivot][k], u[row][k]
                u[pivot][k] = block[0][0] * top + block[0][1] * bottom
                u[row][k] = block[1][0] * top + block[1][1] * bottom
    for basis in range(d):
        phase = cmath.phase(u[basis][basis])
        if abs(phase) > _SYNTHESIS_TOLERANCE:
            with builder.control(ref, basis):
                builder.global_phase(phase)
    for pivot, row, block in reversed(steps):
        inverse = (
            (block[0][0].conjugate(), block[1][0].conjugate()),
            (block[0][1].conjugate(), block[1][1].conjugate()),
        )
        _emit_two_level(builder, ref, pivot, row, inverse)


def _rotation_operation(matrix: Sequence[Sequence[complex]]) -> Operation:
    """把显式小酉矩阵经两能级分解合成为量子操作。"""
    n = (len(matrix) - 1).bit_length()
    b = Builder(_name("df_rotation", matrix), {"target": Bits(n)})
    _emit_unitary(b, b["target"], matrix)
    return annotate(b.finish(), "unitary", algorithm="two_level_synthesis")


def _diagonal_encoding(spectrum: Sequence[float]) -> BlockEncoding:
    """对角阵的受控旋转块编码：单比特信号，``cos(θ_t/2) = g_t/α``，``α = Σ_p |g_p|``。

    对每个基态 ``|t>`` 施加受控 ``Ry(θ_t)``，(0,0) 块恰为 ``diag(g)/α``；
    ``α`` 取谱的 1-范数，与外层 DF/THC 组装的报告口径一致。
    """
    d = len(spectrum)
    n = (d - 1).bit_length()
    alpha = sum(abs(g) for g in spectrum)
    if not alpha:
        raise ValidationError("对角谱不能全部为零")
    b = Builder(
        _name("df_diagonal", spectrum),
        {"target": Bits(n), "signal": Bits(1)},
    )
    for t, g in enumerate(spectrum):
        theta = 2 * math.acos(max(-1.0, min(1.0, g / alpha)))
        if theta:
            with b.control(b["target"], t):
                b.ry(b["signal"], theta)
    return BlockEncoding(
        annotate(b.finish(), "block_encoding", be_alpha=alpha, be_form="diagonal_controlled_rotation")
    )


def _conjugated_encoding(
    rotation: Sequence[Sequence[complex]], spectrum: Sequence[float]
) -> BlockEncoding:
    """``U diag(g) U†`` 的块编码：对角块编码两侧共轭施加两能级合成的 ``U``。"""
    n = (len(rotation) - 1).bit_length()
    rotation_op = _rotation_operation(rotation)
    diagonal = _diagonal_encoding(spectrum)
    b = Builder(
        _name("df_term", rotation, spectrum),
        {"target": Bits(n), "signal": Bits(diagonal.signal_qubits)},
        resources_for(("rot", rotation_op), ("diag", diagonal.operation)),
    )
    with b.adjoint():
        invoke(b, rotation_op, "rot", target=b["target"])
    invoke(b, diagonal.operation, "diag", target=b["target"], signal=b["signal"])
    invoke(b, rotation_op, "rot", target=b["target"])
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=diagonal.alpha))


def double_factorized_encoding(
    df: DoubleFactorization, *, name: str | None = None
) -> BlockEncoding:
    """DF 哈密顿量的 LCU 块编码：外层 PREPARE 在秩指标 ``r`` 上，权重 ``∝ ‖g_r‖₁``。

    SELECT 受控施加 ``U_r diag(g_r)/‖g_r‖₁ U_r†`` 的项块编码，其中对角 ``G_r``
    演化用显式 PREPARE–SELECT 实现，``U_r`` 用两能级分解合成；``scalar`` 项并入
    同一外层 LCU。alpha 取 ``|scalar| + Σ_r ‖g_r‖₁`` 并在 ``df_lambda`` 属性报告。
    产物可直接交给 transforms.qubitization_walk。

    Args:
        df: DF 输入模型实例。
        name: 当前实现未使用，仅为接口兼容保留。

    Returns:
        BlockEncoding: ``H = scalar·I + Σ_r U_r diag(g_r) U_r†`` 的 LCU 块编码。
    """
    require_instance(df, DoubleFactorization, "double_factorized_encoding.df")
    del name
    n = df.width
    terms: list[tuple[float, BlockEncoding]] = []
    if df.scalar:
        terms.append((df.scalar, identity(n)))
    for rotation, spectrum in zip(df.rotations, df.spectra, strict=True):
        if sum(abs(g) for g in spectrum):
            terms.append((1.0, _conjugated_encoding(rotation, spectrum)))
    if not terms:
        raise ValidationError("DF 哈密顿量的所有系数为零")
    out = lcu(terms)
    return BlockEncoding(
        annotate(
            out.operation,
            "block_encoding",
            be_alpha=out.alpha,
            be_form="double_factorization",
            df_rank=df.rank,
            df_lambda=out.alpha,
            lcu_terms=len(terms),
        )
    )


def thc_encoding(thc: THCDecomposition) -> BlockEncoding:
    """THC 哈密顿量的 LCU 块编码：外层 PREPARE 在 ``(μ, ν)`` 对上，权重 ``∝ |ζ_{μν}|·α_μ α_ν``。

    每个 ``(μ, ν)`` 项是 ``L_μ L_ν†`` 的块编码乘积（叶算符经 matrix_pauli_encoding
    编码，``α_μ`` 为其 Pauli l1 上界）。alpha 取 ``Σ_{μν} |ζ_{μν}| α_μ α_ν`` 并在
    ``thc_lambda`` 属性报告。产物可直接交给 transforms.qubitization_walk。

    Args:
        thc: THC 输入模型实例。

    Returns:
        BlockEncoding: ``Σ_{μν} ζ_{μν} L_μ L_ν†`` 的 LCU 块编码。
    """
    require_instance(thc, THCDecomposition, "thc_encoding.thc")
    encodings = [matrix_pauli_encoding(leaf) for leaf in thc.leaves]
    adjoints = [adjoint_be(be) for be in encodings]
    terms: list[tuple[float, BlockEncoding]] = []
    for mu, row in enumerate(thc.coefficients):
        for nu, zeta in enumerate(row):
            if zeta:
                terms.append((zeta, product(encodings[mu], adjoints[nu])))
    if not terms:
        raise ValidationError("THC 哈密顿量的所有系数为零")
    out = lcu(terms)
    return BlockEncoding(
        annotate(
            out.operation,
            "block_encoding",
            be_alpha=out.alpha,
            be_form="thc",
            thc_leaves=thc.leaf_count,
            thc_lambda=out.alpha,
            lcu_terms=len(terms),
        )
    )
