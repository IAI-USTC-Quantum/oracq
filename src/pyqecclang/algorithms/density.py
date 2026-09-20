"""密度矩阵（DM）input model 与 Gibbs 态制备。

DM input model 的核心抽象是纯化访问（purification access）：对密度矩阵 ρ 的访问
定义为制备其纯化态的量子操作 U，U ``|0>`` = ``|ψ_ρ⟩`` 作用在 system 与 environment
两个寄存器上，且对环境取偏迹后 Tr_env ``|ψ⟩⟨ψ|`` = ρ。本模块给出三层范式：
abstract_purification 开放声明、gate_purification 显式小矩阵见证（特征分解后由
多重旋转制备）、以及 PurificationAccess.from_state_preparation 纯态适配（纯态即
环境复净的平凡纯化）。该视图是 B2（量子 SDP）依赖的稳定接口。

Gibbs 态制备走 QSVT 纯化路线（Chowdhury–Somma 2017、van Apeldoorn–Gilyén 2019、
Gilyén et al. 2019, arXiv:1806.01838）：先制备 system 与 environment 上最大混合态
的纯化（n 对 Bell 对），再对 system 作用零信号块正比于 g(H/α) 的 QSVT 块编码，
其中 g(x) = exp(−βα(x+1)/2) 在 [−1,1] 上取值于 (0,1]。g 按奇偶分解为
e^{−c}·cosh(cx) 与 −e^{−c}·sinh(cx)（c = βα/2）两支，各自以修正 Bessel 截断逼近、
经虚部补全合成相位，实部由 (U_Φ + U_{−Φ})/2 提取，最后经 LCU 相加。后置选择
signal == 0 后，system 的约化密度矩阵正比于 g(H/α)² = e^{−βH}（至多相差归一化）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

from pyqecclang.algorithms.contracts import (
    OracleView,
    fail,
    finite_real,
    positive_integer,
    require_instance,
    validate_signature,
)
from pyqecclang.algorithms.operators import BlockEncoding, _name, linear_combination
from pyqecclang.algorithms.oracles import (
    StatePreparation,
    annotate,
    declare,
    gate_state_prep,
    invoke,
    resources_for,
)
from pyqecclang.algorithms.qsvt import (
    _MAX_DEGREE,
    _chebyshev_t,
    _eval,
    _real_qsvt_be,
    _scale,
    _sup_norm,
    _synthesize_with_imag,
    _trim,
)
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError, fuse

__all__ = [
    "ApproximatePurification",
    "PurificationAccess",
    "abstract_purification",
    "gate_purification",
    "gibbs_purification",
    "gibbs_state",
    "maximally_mixed_purification",
    "partial_trace",
    "trace_distance",
]


# ---------------------------------------------------------------------------
# 经典小矩阵工具（纯 Python 复数矩阵，供见证与下游算法的经典侧复用）。
# ---------------------------------------------------------------------------


def _as_complex_matrix(matrix, path):
    try:
        result = tuple(tuple(complex(v) for v in row) for row in matrix)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{path} 需要方阵形式的数值矩阵") from exc
    d = len(result)
    if d < 1 or any(len(row) != d for row in result):
        raise ValidationError(f"{path} 需要非空方阵")
    if not all(math.isfinite(v.real) and math.isfinite(v.imag) for row in result for v in row):
        raise ValidationError(f"{path} 的矩阵元必须有限")
    return result


def _check_hermitian(matrix, path, *, tol=1e-9):
    matrix = _as_complex_matrix(matrix, path)
    d = len(matrix)
    scale = max(1.0, max(abs(v) for row in matrix for v in row))
    for i in range(d):
        for j in range(i + 1, d):
            if abs(matrix[i][j] - matrix[j][i].conjugate()) > tol * scale:
                raise ValidationError(f"{path} 必须是 Hermitian 矩阵")
    return matrix


def _check_density_matrix(rho, *, tol=1e-7):
    matrix = _check_hermitian(rho, "gate_purification.rho")
    d = len(matrix)
    if d < 2 or d & (d - 1):
        raise ValidationError("密度矩阵维度必须是至少为 2 的二的幂")
    if d > 16:
        raise ValidationError("显式纯化见证仅支持不超过 16 维的小密度矩阵")
    if abs(sum(matrix[i][i] for i in range(d)) - 1.0) > tol:
        raise ValidationError("密度矩阵的迹必须为 1")
    return matrix


def _hermitian_eigendecomposition(matrix, *, tol=1e-13, max_sweeps=64):
    """循环 Jacobi 方法求小 Hermitian 矩阵的特征分解。

    返回 (特征值元组, 特征向量表)；特征值按降序排列，特征向量表的第 i 行第 j 列
    是第 j 个特征向量的第 i 个分量。
    """
    d = len(matrix)
    a = [[complex(matrix[i][j]) for j in range(d)] for i in range(d)]
    vectors = [[1.0 + 0j if i == j else 0.0 + 0j for j in range(d)] for i in range(d)]
    scale = max(1.0, max(abs(a[i][j]) for i in range(d) for j in range(d)))
    for _ in range(max_sweeps):
        off = max((abs(a[i][j]) for i in range(d) for j in range(i + 1, d)), default=0.0)
        if off <= tol * scale:
            break
        for p in range(d):
            for q in range(p + 1, d):
                if abs(a[p][q]) <= tol * scale:
                    continue
                # 先做对角相位旋转使非对角元变为正实数，再做实 Jacobi 旋转消元。
                # 相似变换 D†AD 不改变对角元，列缩放后必须恢复 a[q][q]。
                phase = a[p][q] / abs(a[p][q])
                diagonal_qq = a[q][q].real
                for k in range(d):
                    a[k][q] *= phase
                    vectors[k][q] *= phase
                a[q][q] = diagonal_qq + 0j
                for k in range(d):
                    if k != q:
                        a[q][k] = a[k][q].conjugate()
                app, aqq, b = a[p][p].real, a[q][q].real, a[p][q].real
                tau = (aqq - app) / (2 * b)
                t = (1.0 if tau >= 0 else -1.0) / (abs(tau) + math.sqrt(1.0 + tau * tau))
                cos, sin = 1.0 / math.sqrt(1.0 + t * t), t / math.sqrt(1.0 + t * t)
                for k in range(d):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = cos * akp - sin * akq
                    a[k][q] = sin * akp + cos * akq
                for j in range(d):
                    apj, aqj = a[p][j], a[q][j]
                    a[p][j] = cos * apj - sin * aqj
                    a[q][j] = sin * apj + cos * aqj
                for k in range(d):
                    vkp, vkq = vectors[k][p], vectors[k][q]
                    vectors[k][p] = cos * vkp - sin * vkq
                    vectors[k][q] = sin * vkp + cos * vkq
    order = sorted(range(d), key=lambda k: -a[k][k].real)
    values = tuple(a[k][k].real for k in order)
    return values, [[vectors[i][k] for k in order] for i in range(d)]


def partial_trace(amplitudes, system_width, environment_width):
    """对 environment 取偏迹，返回 system 上的约化密度矩阵（行主序嵌套元组）。

    amplitudes 是长度 2^(system_width + environment_width) 的稠密态向量，
    基态下标约定为 system | (environment << system_width)。
    """
    positive_integer(system_width, "partial_trace.system_width", minimum=0, maximum=20)
    positive_integer(environment_width, "partial_trace.environment_width", minimum=0, maximum=20)
    values = [complex(v) for v in amplitudes]
    dim_s, dim_e = 1 << system_width, 1 << environment_width
    if len(values) != dim_s * dim_e:
        raise ValidationError("态向量长度与寄存器宽度不符")
    return tuple(
        tuple(
            sum(
                values[i + (e << system_width)] * values[j + (e << system_width)].conjugate()
                for e in range(dim_e)
            )
            for j in range(dim_s)
        )
        for i in range(dim_s)
    )


def gibbs_state(hamiltonian, beta):
    """经典参考 Gibbs 态 e^{−βH}/Tr(e^{−βH})；仅用于小矩阵的经典见证。"""
    matrix = _check_hermitian(hamiltonian, "gibbs_state.hamiltonian")
    finite_real(beta, "gibbs_state.beta")
    values, vectors = _hermitian_eigendecomposition(matrix)
    d = len(matrix)
    weights = [math.exp(-beta * v) for v in values]
    partition = sum(weights)
    return tuple(
        tuple(
            sum(weights[k] * vectors[i][k] * vectors[j][k].conjugate() for k in range(d))
            / partition
            for j in range(d)
        )
        for i in range(d)
    )


def trace_distance(rho, sigma):
    """迹距离 T(ρ,σ) = ‖ρ−σ‖₁/2，经差矩阵的 Hermitian 特征分解计算。"""
    a = _check_hermitian(rho, "trace_distance.rho")
    b = _check_hermitian(sigma, "trace_distance.sigma")
    if len(a) != len(b):
        raise ValidationError("迹距离要求两个矩阵维度一致")
    d = len(a)
    values, _ = _hermitian_eigendecomposition(
        [[a[i][j] - b[i][j] for j in range(d)] for i in range(d)]
    )
    return 0.5 * sum(abs(v) for v in values)


# ---------------------------------------------------------------------------
# 纯化访问视图：DM input model 的稳定接口（B2 量子 SDP 复用）。
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PurificationAccess(OracleView):
    """密度矩阵的纯化访问视图：制备 ``|ψ_ρ⟩`` 且 Tr_env ``|ψ⟩⟨ψ|`` = ρ。

    操作从全零态出发，在 system 与 environment 两个寄存器上制备 ρ 的纯化态；
    environment 宽度不小于 system 的秩所需位数（一般取等于 system 宽度即可）。
    """

    oracle_kind = "purification_access"
    operation: object

    def purification_access(self):
        """纯化访问的结构化协议访问器，返回 ``self``。

        宿主对象实现同名方法并返回 ``PurificationAccess`` 即可被按协议适配，
        与 ``state_preparation``、``block_encoding`` 等访问器同构。
        """
        return self

    def __post_init__(self):
        validate_signature(self.operation, ("system", "environment"), "PurificationAccess")
        if self.width < 1:
            raise ValidationError("PurificationAccess 的 system 寄存器不能为空")

    @property
    def width(self):
        """system 寄存器位宽，以 RIR 寄存器签名为准。"""
        return next(
            r.type.width for r in self.operation.module.registers if r.name == "system"
        )

    @property
    def environment_width(self):
        """environment 寄存器位宽，以 RIR 寄存器签名为准。"""
        return next(
            r.type.width for r in self.operation.module.registers if r.name == "environment"
        )

    def describe(self):
        """返回类型为 ``purification_access`` 的 ``OracleSpec`` 快照。

        system 宽度记入 main_qubit，environment 宽度记入 anc_qubit。
        """
        from pyqecclang.algorithms.contracts import describe_oracle

        base = describe_oracle(self.operation)
        return replace(
            base,
            type="purification_access",
            main_qubit=self.width,
            anc_qubit=self.environment_width,
        )

    @classmethod
    def from_state_preparation(cls, preparation):
        """纯态即平凡纯化：work 复净的 StatePreparation 适配为 PurificationAccess。

        work 寄存器扮演 environment 角色；制备契约承诺 work 复净，因此偏迹环境后
        system 上仍是原来的纯态。
        """
        require_instance(
            preparation, StatePreparation, "PurificationAccess.from_state_preparation"
        )
        attributes = dict(preparation.operation.module.attributes)
        if preparation.work_width and attributes.get("clean_work") is not True:
            fail(
                "INPUT_PROMISE",
                "PurificationAccess.from_state_preparation.clean_work",
                True,
                attributes.get("clean_work"),
                "需要承诺制备后 work 复净，纯态才能作为平凡纯化",
            )
        n, m = preparation.width, preparation.work_width
        b = Builder(
            _name("pure_state_purification", preparation.operation),
            {"system": Bits(n), "environment": Bits(m)},
            resources_for(("prep", preparation.operation)),
        )
        invoke(b, preparation.operation, "prep", target=b["system"], work=b["environment"])
        return cls(
            annotate(
                b.finish(),
                "unitary",
                density_model="purification_access",
                implementation="pure_state_adapter",
                zero_input=True,
            )
        )

    def as_state_preparation(self):
        """把纯化操作整体视为 system⊕environment 上的 StatePreparation（供 B2 组合）。"""
        n, m = self.width, self.environment_width
        b = Builder(
            _name("purification_as_state_preparation", self.operation),
            {"target": Bits(n + m), "work": Bits(0)},
            resources_for(("purification", self.operation)),
        )
        invoke(
            b,
            self.operation,
            "purification",
            system=b["target"][:n],
            environment=b["target"][n:],
        )
        return StatePreparation(
            annotate(
                b.finish(),
                "state_prep_isometry",
                zero_input=True,
                clean_work=True,
                density_model="purification_access",
            )
        )


@dataclass(frozen=True)
class ApproximatePurification(OracleView):
    """后置选择的近似纯化：signal == 0 分支上 system 与 environment 承载近似纯化态。

    与 PurificationAccess 的区别在于存在信号寄存器与可控近似误差；算法参数
    （如 Gibbs 制备的 beta、error）保存在模块属性中，经 attributes 读取。
    """

    oracle_kind = "approximate_purification"
    operation: object

    def approximate_purification(self):
        """近似纯化的结构化协议访问器，返回 ``self``。

        宿主对象实现同名方法并返回 ``ApproximatePurification`` 即可被按协议适配，
        与 ``state_preparation``、``block_encoding`` 等访问器同构。
        """
        return self

    def __post_init__(self):
        validate_signature(
            self.operation, ("system", "environment", "signal"), "ApproximatePurification"
        )
        if self.width < 1:
            raise ValidationError("ApproximatePurification 的 system 寄存器不能为空")

    @property
    def width(self):
        """system 寄存器位宽，以 RIR 寄存器签名为准。"""
        return next(
            r.type.width for r in self.operation.module.registers if r.name == "system"
        )

    @property
    def environment_width(self):
        """environment 寄存器位宽，以 RIR 寄存器签名为准。"""
        return next(
            r.type.width for r in self.operation.module.registers if r.name == "environment"
        )

    @property
    def signal_qubits(self):
        """signal 寄存器位宽；后置选择要求其读出全 0。"""
        return next(
            r.type.width for r in self.operation.module.registers if r.name == "signal"
        )

    @property
    def attributes(self):
        """模块属性字典的副本；保存 beta、error 等算法参数。"""
        return dict(self.operation.module.attributes)

    @property
    def beta(self):
        """逆温度 β，取自模块属性；未记录时为 ``None``。"""
        return self.attributes.get("beta")

    @property
    def error(self):
        """多项式一致逼近误差参数，取自模块属性；未记录时为 ``None``。"""
        return self.attributes.get("error")

    def describe(self):
        """返回类型为 ``approximate_purification`` 的 ``OracleSpec`` 快照。

        system 宽度记入 main_qubit，environment 与 signal 宽度之和记入 anc_qubit。
        """
        from pyqecclang.algorithms.contracts import describe_oracle

        base = describe_oracle(self.operation)
        return replace(
            base,
            type="approximate_purification",
            main_qubit=self.width,
            anc_qubit=self.environment_width + self.signal_qubits,
        )


# ---------------------------------------------------------------------------
# 三层范式：abstract 开放声明与 gate 见证实现。
# ---------------------------------------------------------------------------


def abstract_purification(name, width, environment_width=None, *, reversible=True):
    """DM input model 的开放声明：制备 ``|ψ_ρ⟩`` 的纯化访问槽，供分批绑定。

    environment_width 缺省取 width（任何密度矩阵都有等宽环境的纯化）；
    声明经 linking.bind 绑定 gate_purification 等见证实现后程序闭合。
    """
    positive_integer(width, "abstract_purification.width", maximum=63)
    environment_width = width if environment_width is None else environment_width
    positive_integer(
        environment_width, "abstract_purification.environment_width", minimum=0, maximum=63
    )
    return PurificationAccess(
        declare(
            name,
            {"system": Bits(width), "environment": Bits(environment_width)},
            paradigm="unitary",
            attributes={"density_model": "purification_access", "zero_input": True},
            supports_adjoint=reversible,
            supports_controlled=reversible,
        )
    )


def gate_purification(rho, *, name=None):
    """显式小密度矩阵的纯化见证：特征分解 ρ = Σ_j p_j ``|v_j⟩⟨v_j|`` 后受控制备。

    纯化态取 ``|ψ_ρ⟩ = Σ_j √p_j |v_j⟩_s |j⟩_e``，其幅度向量经
    gate_state_prep 的多重旋转树在 system 与 environment 的拼接寄存器上制备；
    对环境取偏迹恰好回到 ρ。
    """
    matrix = _check_density_matrix(rho)
    d = len(matrix)
    n = (d - 1).bit_length()
    values, vectors = _hermitian_eigendecomposition(matrix)
    if values[-1] < -1e-7:
        fail(
            "INPUT_PROMISE",
            "gate_purification.rho",
            "positive semidefinite",
            values[-1],
            "密度矩阵必须半正定",
        )
    amplitudes = [0j] * (d * d)
    for j, p in enumerate(values):
        root = math.sqrt(max(p, 0.0))
        for i in range(d):
            amplitudes[i + (j << n)] = root * vectors[i][j]
    preparation = gate_state_prep(amplitudes)
    b = Builder(
        name or _name("purification", matrix),
        {"system": Bits(n), "environment": Bits(n)},
        resources_for(("prep", preparation.operation)),
    )
    invoke(
        b,
        preparation.operation,
        "prep",
        target=fuse(b["system"], b["environment"]),
        work=b["system"][:0],
    )
    return PurificationAccess(
        annotate(
            b.finish(),
            "unitary",
            density_model="purification_access",
            implementation="gate_eigendecomposition",
            zero_input=True,
        )
    )


def maximally_mixed_purification(width, *, name=None):
    """最大混合态 I/2^n 的纯化生成器：n 对 Bell 对 ``|Φ+⟩`` 的张量积。"""
    positive_integer(width, "maximally_mixed_purification.width", maximum=32)
    b = Builder(
        name or f"bell_purification_{width}",
        {"system": Bits(width), "environment": Bits(width)},
    )
    for bit in range(width):
        b.h(b["system"][bit])
        b.xor(b["system"][bit], b["environment"][bit])
    return PurificationAccess(
        annotate(
            b.finish(),
            "unitary",
            density_model="purification_access",
            implementation="bell_pairs",
            zero_input=True,
        )
    )


# ---------------------------------------------------------------------------
# Gibbs 态制备：QSVT 纯化路线。
# ---------------------------------------------------------------------------


def _bessel_i(n, x):
    """第一类修正 Bessel 函数 I_n(x)，幂级数纯 Python 实现。"""
    term = (x / 2) ** n / math.factorial(n)
    total = term
    m = 0
    while term > 1e-18 * max(1.0, total) and m < 100000:
        term *= ((x / 2) ** 2) / ((m + 1) * (m + n + 1))
        total += term
        m += 1
    return total


def _gibbs_branches(c, error):
    """g(x) = e^{−c(x+1)} 的偶/奇 Chebyshev 截断：e^{−c}cosh(cx) 与 −e^{−c}sinh(cx)。

    每支截断尾部按 2e^{−c}·Σ_{k>d} I_k(c) ≤ error/8 控制，合计一致误差不超过
    error/4；返回 (偶支升幂系数, 奇支升幂系数, 偶支度数, 奇支度数)。
    """
    kmax = min(_MAX_DEGREE, int(math.ceil(c)) + 8 * int(math.ceil(math.log10(8 / error))) + 8)
    ivals = [_bessel_i(k, c) for k in range(kmax + 2)]
    suffix = [0.0] * (kmax + 3)
    for j in range(kmax + 1, -1, -1):
        suffix[j] = suffix[j + 1] + ivals[j]

    def tail_ok(k):
        return 2.0 * math.exp(-c) * suffix[k + 1] <= error / 8

    d_even = next((k for k in range(2, kmax + 1, 2) if tail_ok(k)), None)
    d_odd = next((k for k in range(1, kmax + 1, 2) if tail_ok(k)), None)
    if d_even is None or d_odd is None:
        raise ValidationError(
            "β·α 过大：Gibbs 多项式度数超过合成上限；请减小 β 或先缩小 H 的谱尺度"
        )
    shift = math.exp(-c)
    even = [0.0] * (d_even + 1)
    for k in range(d_even // 2 + 1):
        coef = shift * (1.0 if k == 0 else 2.0) * ivals[2 * k]
        for i, v in enumerate(_chebyshev_t(2 * k)):
            even[i] += coef * v
    odd = [0.0] * (d_odd + 1)
    for k in range((d_odd + 1) // 2):
        coef = -2.0 * shift * ivals[2 * k + 1]
        for i, v in enumerate(_chebyshev_t(2 * k + 1)):
            odd[i] += coef * v
    return _trim(even), _trim(odd), d_even, d_odd


def _even_imag_candidates(f):
    d = len(f) - 1
    a0 = math.sqrt(max(0.0, 1.0 - _eval(f, 0.0).real ** 2))
    a1 = math.sqrt(max(0.0, 1.0 - _eval(f, 1.0).real ** 2))
    for m in range(1, d // 2 + 1):
        yield tuple(
            (a0 if i == 0 else 0.0) + ((a1 - a0) if i == 2 * m else 0.0)
            for i in range(d + 1)
        )


def _odd_imag_candidates(f):
    d = len(f) - 1
    a1 = math.sqrt(max(0.0, 1.0 - _eval(f, 1.0).real ** 2))
    for m in range(0, (d - 1) // 2 + 1):
        yield tuple(a1 * (1.0 if i == 2 * m + 1 else 0.0) for i in range(d + 1))


def gibbs_purification(hamiltonian, beta, *, error=0.01):
    """Gibbs 态 ρ = e^{−βH}/Z 的近似纯化制备（QSVT 纯化路线）。

    hamiltonian 是 H 的 BlockEncoding，约定谱含于 [−α,α]（α = be_alpha）；
    谱变量 x = λ/α ∈ [−1,1] 上目标函数 g(x) = exp(−βα(x+1)/2) ∈ (0,1]。
    构造分两步：(a) 在 system 与 environment 上制备 n 对 Bell 对（最大混合态
    的纯化 ``|Φ⟩``）；(b) 对 system 作用零信号块为 g(H/α)/(2s) 的 QSVT 块编码。
    由于 g 的偶支 e^{−c}cosh(cx) 与奇支 −e^{−c}sinh(cx) 均为凸函数，端点匹配的
    虚部补全恒满足单位圆盘约束（f²(x) 不超过端点连线），相位合成必然可行。
    后置选择 signal == 0 后态正比于 (g(H/α) ⊗ I)``|Φ⟩``，system 的约化密度矩阵
    即为 e^{−βH}/Z；归一化因子 2s 与配分函数无关，不影响约化态。

    error 控制多项式一致逼近误差（每支截断尾部 ≤ error/8）；返回
    ApproximatePurification，模块属性含 algorithm="gibbs_purification"、beta、
    error、qsp_degree 与 gibbs_scale = 2s。β = 0 时退化为最大混合态纯化。
    """
    require_instance(hamiltonian, BlockEncoding, "gibbs_purification.hamiltonian")
    finite_real(beta, "gibbs_purification.beta", minimum=0)
    finite_real(error, "gibbs_purification.error", minimum=0, strict=True)
    if error >= 1:
        raise ValidationError("近似误差 error 必须在 (0,1) 内")
    n = hamiltonian.width
    if beta == 0:
        bells = maximally_mixed_purification(n)
        b = Builder(
            _name("gibbs_purification", hamiltonian.operation, 0.0, error),
            {"system": Bits(n), "environment": Bits(n), "signal": Bits(0)},
            resources_for(("bells", bells.operation)),
        )
        invoke(b, bells.operation, "bells", system=b["system"], environment=b["environment"])
        return ApproximatePurification(
            annotate(
                b.finish(),
                "unitary",
                density_model="approximate_purification",
                algorithm="gibbs_purification",
                beta=0.0,
                error=float(error),
                qsp_degree=0,
                gibbs_scale=1.0,
                success_condition="signal == 0",
            )
        )
    c = float(beta) * hamiltonian.alpha / 2.0
    g_even, g_odd, d_even, d_odd = _gibbs_branches(c, error)
    s = 1.5 * max(_sup_norm(g_even), _sup_norm(g_odd), 1e-3)
    f_even, f_odd = _scale(1.0 / s, g_even), _scale(1.0 / s, g_odd)
    phases_even = _synthesize_with_imag(f_even, _even_imag_candidates(f_even), "Gibbs 偶支")
    phases_odd = _synthesize_with_imag(f_odd, _odd_imag_candidates(f_odd), "Gibbs 奇支")
    be_even = _real_qsvt_be(hamiltonian, phases_even)
    be_odd = _real_qsvt_be(hamiltonian, phases_odd)
    gibbs_be = linear_combination(1.0, be_even, 1.0, be_odd)
    b = Builder(
        _name("gibbs_purification", hamiltonian.operation, beta, error),
        {"system": Bits(n), "environment": Bits(n), "signal": Bits(gibbs_be.signal_qubits)},
        resources_for(("gibbs_be", gibbs_be.operation)),
    )
    for bit in range(n):
        b.h(b["system"][bit])
        b.xor(b["system"][bit], b["environment"][bit])
    invoke(b, gibbs_be.operation, "gibbs_be", target=b["system"], signal=b["signal"])
    return ApproximatePurification(
        annotate(
            b.finish(),
            "unitary",
            density_model="approximate_purification",
            algorithm="gibbs_purification",
            beta=float(beta),
            error=float(error),
            qsp_degree=max(d_even, d_odd),
            gibbs_scale=2.0 * s,
            success_condition="signal == 0",
            spectral_variable="x = eigenvalue/alpha in [-1,1]",
            target_function="g(x) = exp(-beta*alpha*(x+1)/2)",
        )
    )
