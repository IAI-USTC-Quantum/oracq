"""线性系统输入契约、Costa walk 与 CKS Chebyshev 求解器。"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, replace

from pyqecclang.algorithms.arithmetic import FixedFormat
from pyqecclang.algorithms.block_encoding import lcu, reflect_zero
from pyqecclang.algorithms.contracts import (
    ContractIssue,
    ContractReport,
    InputRequirement,
    ProtocolContract,
    finite_real,
    positive_integer,
    require_instance,
)
from pyqecclang.algorithms.interfaces import (
    BlockEncodingProtocol,
    CKSSparseProtocol,
    StateOracleProtocol,
    StatePreparationProtocol,
    as_block_encoding,
    as_qlss_matrix,
    as_sparse_access,
    as_state_preparation,
    operator_state_contract,
)
from pyqecclang.algorithms.operators import BlockEncoding, _name
from pyqecclang.algorithms.oracles import (
    SparseAccess,
    StateOracle,
    StatePreparation,
    annotate,
    invoke,
    resources_for,
)
from pyqecclang.algorithms.sparse import chebyshev_block, real_symmetric_sparse_encoding
from pyqecclang.algorithms.state_preparation import apply_be_to_state, select_subspace
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError, fuse


@dataclass(frozen=True)
class SpectralPromise:
    norm_upper: float
    sigma_min_lower: float
    evidence: str = "caller_declared_unverified"

    def __post_init__(self):
        finite_real(self.norm_upper, "spectrum.norm_upper", minimum=0, strict=True)
        finite_real(self.sigma_min_lower, "spectrum.sigma_min_lower", minimum=0, strict=True)
        if not 0 < self.sigma_min_lower <= self.norm_upper or not math.isfinite(self.norm_upper):
            raise ValidationError("需要有限正的范数上界及最小奇异值下界")

    def inverse_bound(self, alpha):
        if alpha < self.sigma_min_lower:
            raise ValidationError("BE alpha 与奇异值下界声明冲突")
        return max(1.0, alpha / self.sigma_min_lower)


@dataclass(frozen=True)
class SparseSystem:
    """有限元素字与同一 RHS 的稀疏 Hermitian 问题；矩阵性质由调用者声明。"""

    access: SparseAccess
    value_format: FixedFormat
    entry_bound: float
    rhs: StatePreparation
    spectrum: SpectralPromise
    diagonal_nonnegative: bool = False
    hermitian: bool = False

    def __post_init__(self):
        object.__setattr__(self, "access", as_sparse_access(self.access))
        object.__setattr__(self, "rhs", as_state_preparation(self.rhs))
        require_instance(self.value_format, FixedFormat, "SparseSystem.value_format")
        require_instance(self.rhs, StatePreparation, "SparseSystem.rhs")
        require_instance(self.spectrum, SpectralPromise, "SparseSystem.spectrum")
        finite_real(self.entry_bound, "SparseSystem.entry_bound", minimum=0, strict=True)
        if type(self.hermitian) is not bool or type(self.diagonal_nonnegative) is not bool:
            raise ValidationError("稀疏矩阵性质声明必须是 bool")
        if (
            self.value_format.width != self.access.value_width
            or self.rhs.width != self.access.width
        ):
            raise ValidationError("稀疏矩阵、值格式与 RHS 的布局不匹配")
        if not math.isfinite(self.entry_bound) or self.entry_bound <= 0:
            raise ValidationError("元素幅值上界必须为有限正数")


@dataclass(frozen=True)
class BlockSystem:
    encoding: BlockEncoding
    rhs: StatePreparation
    spectrum: SpectralPromise

    def __post_init__(self):
        object.__setattr__(self, "encoding", as_block_encoding(self.encoding))
        object.__setattr__(self, "rhs", as_state_preparation(self.rhs))
        require_instance(self.rhs, StatePreparation, "BlockSystem.rhs")
        require_instance(self.spectrum, SpectralPromise, "BlockSystem.spectrum")
        if self.encoding.width != self.rhs.width:
            raise ValidationError("BE 与 RHS 宽度不匹配")

    @property
    def inverse_norm_bound(self):
        return self.spectrum.inverse_bound(self.encoding.alpha)


@dataclass(frozen=True)
class LinearSystem:
    sparse: SparseSystem | None = None
    block: BlockSystem | None = None
    physical_width: int | None = None
    physical_high_value: int = 0
    rhs_norm: float | None = None
    data_assumptions: tuple[str, ...] = ()

    def __post_init__(self):
        if (self.sparse is None) == (self.block is None):
            raise ValidationError("线性问题需要且只能指定一种源输入模型，其他模型由显式适配生成")
        source = self.block if self.block is not None else self.sparse
        require_instance(
            source, BlockSystem if self.block is not None else SparseSystem, "LinearSystem.source"
        )
        width = source.encoding.width if self.block is not None else source.access.width
        if self.physical_width is not None:
            positive_integer(self.physical_width, "LinearSystem.physical_width", maximum=width)
        positive_integer(
            self.physical_high_value,
            "LinearSystem.physical_high_value",
            minimum=0,
            maximum=(1 << (width - (self.physical_width or width))) - 1,
        )
        if self.rhs_norm is not None:
            finite_real(self.rhs_norm, "LinearSystem.rhs_norm", minimum=0)
        if self.rhs_norm is not None and (self.rhs_norm < 0 or not math.isfinite(self.rhs_norm)):
            raise ValidationError("右端范数声明无效")

    def block_input(self):
        if self.block is not None:
            return self.block, ("block input supplied",)
        from pyqecclang.algorithms.sparse import real_symmetric_sparse_encoding

        s = self.sparse
        if not s.hermitian:
            raise ValidationError(
                "当前稀疏适配需要显式 Hermitian 声明；非 Hermitian 系统需先提供行列访问并扩张"
            )
        encoded = real_symmetric_sparse_encoding(
            s.access, s.value_format, s.entry_bound, diagonal_nonnegative=s.diagonal_nonnegative
        )
        return BlockSystem(encoded, s.rhs, s.spectrum), (
            "sparse location + entry",
            "CKS Tdag S T to BE",
        )


@dataclass(frozen=True)
class SolveResult:
    state: StateOracle
    norm_probe: StateOracle
    input_alpha: float
    encoded_inverse_bound: float
    rhs_norm: float | None
    adapter_trace: tuple[str, ...]
    kernel_status: str = "prototype; solver accuracy pending"

    def state_oracle(self):
        return self.state

    @property
    def operation(self):
        return self.state.operation

    def recover_norm(self, solver_success, joint_matrix_success):
        if self.rhs_norm is None:
            raise ValidationError("缺少经典右端范数，不能恢复物理更新幅值")
        if not 0 < joint_matrix_success <= solver_success <= 1:
            raise ValidationError("需要有效的求解成功概率和求解+矩阵探针联合成功概率")
        return self.rhs_norm / (self.input_alpha * math.sqrt(joint_matrix_success / solver_success))


@dataclass(frozen=True)
class QLSSProtocol:
    provides = (StateOracleProtocol,)
    name: str
    input_model: str
    kernel: object
    legacy: object = None

    def __post_init__(self):
        if not self.name or self.input_model not in {"sparse", "block_encoding"}:
            raise ValidationError("QLSSProtocol 的名称或 input_model 无效")
        if not callable(self.kernel) or (self.legacy is not None and not callable(self.legacy)):
            raise ValidationError("QLSSProtocol 内核必须可调用")

    @property
    def contract(self):
        types = (
            (CKSSparseProtocol,)
            if self.input_model == "sparse"
            else (BlockEncodingProtocol, CKSSparseProtocol)
        )
        return ProtocolContract(
            self.name,
            (
                InputRequirement("A", types, adapter=as_qlss_matrix, adjoint=True, controlled=True),
                InputRequirement(
                    "b",
                    (StatePreparationProtocol,),
                    adapter=as_state_preparation,
                    adjoint=True,
                    controlled=True,
                    zero_input=True,
                    clean_work=True,
                ),
            ),
            same_width=(("A", "b"),),
            assumptions=("谱界与矩阵解释是调用者声明；算法精度与成功通道待核验",),
        )

    def check(self, problem):
        """不运行内核或适配器；给出结构、能力及所需适配的完整报告。"""
        if not isinstance(problem, LinearSystem):
            return ContractReport(
                self.name,
                (),
                (
                    ContractIssue(
                        "INPUT_TYPE",
                        self.name,
                        "LinearSystem",
                        type(problem).__name__,
                        "需要 LinearSystem",
                    ),
                ),
            )
        source = problem.block if problem.block is not None else problem.sparse
        a = source.encoding if problem.block is not None else source.access
        report = self.contract.check(A=a, b=source.rhs)
        issues, adapters = list(report.issues), []
        if problem.rhs_norm == 0:
            issues.append(
                ContractIssue(
                    "INPUT_ZERO_RHS", self.name + ".b", "nonzero vector", 0, "零右端应在经典侧处理"
                )
            )
        if self.input_model == "sparse" and problem.sparse is None:
            issues.append(
                ContractIssue(
                    "INPUT_ADAPTER",
                    self.name + ".A",
                    "cks_sparse",
                    "block_encoding",
                    "不能从一般 BE 自动恢复高效稀疏 oracle",
                )
            )
        if problem.sparse is not None:
            for key in ("hermitian", "diagonal_nonnegative"):
                if getattr(source, key) is not True:
                    issues.append(
                        ContractIssue(
                            "INPUT_PROMISE",
                            self.name + ".A." + key,
                            True,
                            False,
                            "当前 CKS 稀疏适配需要 Hermitian 与非负对角声明",
                        )
                    )
            adapters.append("CKS real symmetric sparse -> Tdag S T block encoding")
            alpha = source.access.sparsity * source.entry_bound
        else:
            alpha = source.encoding.alpha
        if alpha < source.spectrum.sigma_min_lower:
            issues.append(
                ContractIssue(
                    "INPUT_SPECTRUM",
                    self.name + ".A.alpha",
                    f">={source.spectrum.sigma_min_lower}",
                    alpha,
                    "BE alpha 与奇异值下界声明冲突",
                )
            )
        return replace(report, issues=tuple(issues), adapters=tuple(adapters))

    def __call__(self, *args):
        if len(args) == 1 and isinstance(args[0], LinearSystem):
            return self.solve(args[0])
        if len(args) == 2 and self.legacy is not None:
            operator_state_contract(self.name, matrix="A", state="b").check(
                A=args[0], b=args[1]
            ).require()
            return self.legacy(as_block_encoding(args[0]), as_state_preparation(args[1]))
        raise ValidationError("QLSS protocol 需要 LinearSystem；不能根据调用形状猜测输入模型")

    def solve(self, problem):
        self.check(problem).require()
        if problem.rhs_norm == 0:
            raise ValidationError("零右端应在经典侧返回零更新，不存在归一化 RHS 态")
        if self.input_model == "sparse":
            if problem.sparse is None:
                raise ValidationError("不能从一般 BE 自动恢复高效稀疏 oracle")
            source = problem.sparse
            state = self.kernel(source)
            block, trace = problem.block_input()
            trace = ("sparse input consumed by " + self.name,) + trace
        elif self.input_model == "block_encoding":
            block, trace = problem.block_input()
            state = self.kernel(block)
            trace = trace + (self.name + " consumes BE",)
        else:
            raise ValidationError("未知 QLSS input_model")
        require_instance(state, StateOracle, self.name + ".output")
        from pyqecclang.infrastructure.linking import capabilities

        capability = capabilities(state.operation.program())
        if not capability["supports_adjoint"] or not capability["supports_controlled"]:
            raise ValidationError("QFVM 的可组合 QLSS 输出需要 adjoint 与 controlled 能力")
        if state.width != block.encoding.width:
            raise ValidationError("QLSS 输出布局与所求矩阵不匹配")
        width = problem.physical_width or state.width
        selected = select_subspace(
            state, width, problem.physical_high_value, label="qlss_physical_solution"
        )
        attributes = dict(selected.operation.module.attributes)
        attributes.update(
            qlss_protocol=self.name,
            qlss_input_model=self.input_model,
            input_alpha=block.encoding.alpha,
            encoded_inverse_bound=block.inverse_norm_bound,
            spectral_norm_upper=block.spectrum.norm_upper,
            sigma_min_lower=block.spectrum.sigma_min_lower,
            spectral_evidence=block.spectrum.evidence,
            adapter_trace=json.dumps(trace),
            data_assumptions=json.dumps(problem.data_assumptions),
            kernel_status="prototype; success channel and accuracy require validation",
            norm_recovery="rhs_norm / (input_alpha * sqrt(p_joint/p_solver))",
        )
        selected = StateOracle(
            replace(
                selected.operation,
                module=replace(
                    selected.operation.module, attributes=tuple(sorted(attributes.items()))
                ),
            )
        )
        # 独立矩阵范数探针：不能拿 Costa filtering 的成功率套用 CKS 的 LCU 因子。
        a = block.encoding
        probe = Builder(
            _name("solution_matrix_norm_probe", selected.operation, a.operation),
            {"target": Bits(a.width), "signal": Bits(selected.signal_qubits + a.signal_qubits)},
            resources_for(("solution", selected.operation), ("matrix", a.operation)),
        )
        high = probe["target"][width:]
        for bit in range(high.width):
            if (problem.physical_high_value >> bit) & 1:
                probe.x(high[bit])
        invoke(
            probe,
            selected.operation,
            "solution",
            target=probe["target"][:width],
            signal=probe["signal"][: selected.signal_qubits],
        )
        invoke(
            probe,
            a.operation,
            "matrix",
            target=probe["target"],
            signal=probe["signal"][selected.signal_qubits :],
        )
        norm_probe = StateOracle(
            annotate(
                probe.finish(),
                "unitary",
                algorithm="solution_matrix_norm_probe",
                probability_contract="joint solver success and fresh matrix-ancilla success",
            )
        )
        return SolveResult(
            selected, norm_probe, a.alpha, block.inverse_norm_bound, problem.rhs_norm, trace
        )


@dataclass(frozen=True)
class CostaConfig:
    """kappa 为实际编码矩阵 A/alpha 的逆谱界；并非任意尺度下的 cond(A)。"""

    steps: int = 2
    kappa: float = 4.0
    schedule_power: float = 1.5
    filter_degree: int = 2
    filter_attenuation: float = 0.2

    def __post_init__(self):
        positive_integer(self.steps, "CostaConfig.steps")
        finite_real(self.kappa, "CostaConfig.kappa", minimum=1)
        finite_real(self.schedule_power, "CostaConfig.schedule_power")
        if self.schedule_power == 1:
            raise ValidationError("Costa schedule_power 不能等于 1")
        positive_integer(self.filter_degree, "CostaConfig.filter_degree", minimum=2)
        if self.filter_degree % 2:
            raise ValidationError("Costa filter_degree 需要偶数")
        finite_real(
            self.filter_attenuation, "CostaConfig.filter_attenuation", minimum=0, strict=True
        )
        if self.filter_attenuation >= 1:
            raise ValidationError("Costa filter_attenuation 必须小于 1")


@dataclass(frozen=True)
class FilterPlan:
    weights: tuple[float, ...]
    stride: int = 1
    offset: int = 0
    method: str = "explicit"


def dolph_chebyshev_plan(degree=2, attenuation=0.2):
    """以 Laurent 多项式递推构造一个偶数阶 Dolph–Chebyshev LCU。"""
    if type(degree) is not int or degree < 2 or degree % 2:
        raise ValidationError("Dolph–Chebyshev 原型需要正偶数阶")
    if not 0 < attenuation < 1:
        raise ValidationError("attenuation 必须介于 0 和 1 之间")
    beta = math.cosh(math.acosh(1 / attenuation) / degree)
    previous, current = {0: 1.0}, {-1: beta / 2, 1: beta / 2}
    for _ in range(2, degree + 1):
        following = {}
        for power, value in current.items():
            following[power - 1] = following.get(power - 1, 0.0) + beta * value
            following[power + 1] = following.get(power + 1, 0.0) + beta * value
        for power, value in previous.items():
            following[power] = following.get(power, 0.0) - value
        previous, current = current, following
    weights = tuple(
        max(0.0, attenuation * current.get(k, 0.0)) for k in range(-degree, degree + 1, 2)
    )
    return FilterPlan(weights, stride=2, offset=-degree, method="dolph_chebyshev")


def schedule(s, kappa, power=1.5):
    if kappa < 1 or power == 1 or not 0 <= s <= 1:
        raise ValidationError("Costa 调度参数无效")
    if kappa == 1:
        return s
    return kappa / (kappa - 1) * (1 - (1 + s * (kappa ** (power - 1) - 1)) ** (1 / (1 - power)))


def costa_walk(a: BlockEncoding, bprep: StatePreparation, fs: float):
    operator_state_contract("costa_walk", matrix="A", state="b").check(A=a, b=bprep).require()
    finite_real(fs, "costa_walk.fs", minimum=0)
    a, bprep = as_block_encoding(a), as_state_preparation(bprep)
    if a.width != bprep.width or not 0 <= fs <= 1:
        raise ValidationError("Costa 输入形状或调度点不匹配")
    total = a.signal_qubits + bprep.work_width + 4
    b = Builder(
        _name("costa_walk", a.operation, bprep.operation, fs),
        {"target": Bits(a.width), "signal": Bits(total)},
        resources_for(("a", a.operation), ("b", bprep.operation)),
        attributes={"algorithm": "costa_general_walk", "validation_stage": "paradigm", "fs": fs},
    )
    enc = b["signal"][: a.signal_qubits]
    bw = b["signal"][a.signal_qubits : a.signal_qubits + bprep.work_width]
    start = a.signal_qubits + bprep.work_width
    a1, a2, a3, a4 = (b["signal"][start + i] for i in range(4))

    def prep(inverse=False):
        if inverse:
            with b.adjoint():
                invoke(b, bprep.operation, "b", target=b["target"], work=bw)
        else:
            invoke(b, bprep.operation, "b", target=b["target"], work=bw)

    def reflect_rhs_input_zero():
        # U_b 是 target+work 上的酉扩张，投影必须同时要求工作位为零。
        with b.control(b["target"], 0):
            if bw.width:
                with b.control(bw, 0):
                    b.global_phase(math.pi)
            else:
                b.global_phase(math.pi)

    def rotation():
        # R(s) 是反射矩阵，写成 Ry(2 atan2(f,1-f)) Z。
        b.z(a2)
        b.ry(a2, 2 * math.atan2(fs, 1 - fs))

    b.h(a3)
    prep(True)
    b.x(a1)
    with b.control(fuse(a1, a3, a4)):
        reflect_rhs_input_zero()
    b.x(a1)
    prep()
    b.x(a4)
    with b.control(a4):
        rotation()
    b.x(a4)
    with b.control(a4):
        b.h(a2)
    with b.control(fuse(a1, a2)):
        invoke(b, a.operation, "a", target=b["target"], signal=enc)
    with b.control(a2):
        b.x(a1)
    with b.control(a1):
        reflect_zero(b, a2)
    with b.control(fuse(a1, a2)):
        with b.adjoint():
            invoke(b, a.operation, "a", target=b["target"], signal=enc)
    b.x(a4)
    with b.control(a4):
        b.h(a2)
    b.x(a4)
    with b.control(a4):
        rotation()
    b.x(a4)
    prep(True)
    b.x(a1)
    with b.control(fuse(a1, a3, a4)):
        reflect_rhs_input_zero()
    b.x(a1)
    prep()
    b.h(a3)
    reflect_zero(b, fuse(enc, a2, a3), positive=True)
    b.global_phase(math.pi / 2)
    return b.finish()


def unary_weight_preparation(weights):
    weights = tuple(float(v) for v in weights)
    if not weights or any(not math.isfinite(v) or v < 0 for v in weights) or sum(weights) <= 0:
        raise ValidationError("unary PREP 需要非负、非零总和的权重")
    width = len(weights) - 1
    if width < 1:
        raise ValidationError("unary PREP 至少需要两个权重")
    b = Builder(_name("unary_prepare", weights), {"target": Bits(width)})
    for bit in range(width):
        remaining = sum(weights[bit:])
        tail = sum(weights[bit + 1 :])
        angle = 0.0 if remaining == 0 else 2 * math.asin(math.sqrt(tail / remaining))
        if bit:
            with b.control(b["target"][bit - 1]):
                b.ry(b["target"][bit], angle)
        else:
            b.ry(b["target"][bit], angle)
    return annotate(b.finish(), "state_prep_isometry", zero_input=True, encoding="unary_prefix")


def lcu_filter(walk, plan: FilterPlan):
    widths = {r.name: r.type.width for r in walk.module.registers}
    if set(widths) != {"target", "signal"}:
        raise ValidationError("filtering 接收 target/signal walk 接口")
    if type(plan.stride) is not int or plan.stride < 1 or type(plan.offset) is not int:
        raise ValidationError("filter 幂配置无效")
    prep = unary_weight_preparation(plan.weights)
    clock_width = len(plan.weights) - 1
    b = Builder(
        _name("filter", walk, plan.weights, plan.stride, plan.offset),
        {"target": Bits(widths["target"]), "signal": Bits(widths["signal"] + clock_width)},
        resources_for(("walk", walk)),
        attributes={
            "algorithm": "coherent_lcu_filter",
            "filter_method": plan.method,
            "filter_terms": len(plan.weights),
            "stride": plan.stride,
            "offset": plan.offset,
            "validation_stage": "paradigm",
        },
    )
    work, clock = b["signal"][: widths["signal"]], b["signal"][widths["signal"] :]

    def repeat_walk(count):
        with b.repeat(count):
            invoke(b, walk, "walk", target=b["target"], signal=work)

    if plan.offset < 0:
        with b.adjoint():
            repeat_walk(-plan.offset)
    else:
        repeat_walk(plan.offset)
    invoke(b, prep, target=clock)
    for bit in range(clock_width):
        with b.control(clock[bit]):
            repeat_walk(plan.stride)
    with b.adjoint():
        invoke(b, prep, target=clock)
    return b.finish()


def costa_qlss(a: BlockEncoding, bprep: StatePreparation, config=None, *, filtering=None):
    a, bprep = as_block_encoding(a), as_state_preparation(bprep)
    config = CostaConfig() if config is None else config
    if type(config.steps) is not int or config.steps < 1:
        raise ValidationError("Costa steps 必须为正整数")
    walks = tuple(
        costa_walk(a, bprep, schedule((i + 1) / config.steps, config.kappa, config.schedule_power))
        for i in range(config.steps)
    )
    plan = filtering or dolph_chebyshev_plan(config.filter_degree, config.filter_attenuation)
    final_filter = lcu_filter(walks[-1], plan)
    ws = a.signal_qubits + bprep.work_width + 4
    signal_width = ws + len(plan.weights) - 1
    operands = [("prep", bprep.operation), ("filter", final_filter)]
    operands += [(f"walk{i}", walk) for i, walk in enumerate(walks)]
    b = Builder(
        _name("costa_qlss", a.operation, bprep.operation, config, plan),
        {"target": Bits(a.width), "signal": Bits(signal_width)},
        resources_for(*operands),
        attributes={
            "algorithm": "costa_qlss",
            "input_alpha": a.alpha,
            "encoded_inverse_bound": config.kappa,
            "normalization_assumption": "sigma_min(A / alpha) >= 1 / kappa",
            "kernel_status": "prototype; initial walk eigenstate and readout channel unverified",
            "validation_stage": "paradigm",
            "steps": config.steps,
            "filtering": plan.method,
            "success_condition": "signal == 0; probability and solution accuracy unverified",
        },
    )
    bw = b["signal"][a.signal_qubits : a.signal_qubits + bprep.work_width]
    invoke(b, bprep.operation, "prep", target=b["target"], work=bw)
    for i, walk in enumerate(walks):
        invoke(b, walk, f"walk{i}", target=b["target"], signal=b["signal"][:ws])
    invoke(b, final_filter, "filter", target=b["target"], signal=b["signal"])
    return StateOracle(b.finish())


def make_costa_qlss(config=None):
    """声明 BE 输入；问题层自动按 alpha/sigma_min 推导实际调度参数。"""
    from dataclasses import replace

    config = CostaConfig() if config is None else config
    require_instance(config, CostaConfig, "make_costa_qlss.config")

    def kernel(system):
        effective = replace(config, kappa=system.inverse_norm_bound)
        return costa_qlss(system.encoding, system.rhs, effective)

    return QLSSProtocol(
        "costa_general_walk", "block_encoding", kernel, legacy=lambda a, b: costa_qlss(a, b, config)
    )


@dataclass(frozen=True)
class CKSConfig:
    order: int = 2
    terms: int | None = None

    def __post_init__(self):
        positive_integer(self.order, "CKSConfig.order", maximum=128)
        if self.terms is not None:
            positive_integer(self.terms, "CKSConfig.terms", maximum=self.order)

    def coefficients(self):
        if type(self.order) is not int or not 1 <= self.order <= 128:
            raise ValidationError("CKS 基础原型的 order 范围为 1..128")
        terms = self.order if self.terms is None else self.terms
        if not 1 <= terms <= self.order:
            raise ValidationError("CKS 截断项数无效")
        return tuple(
            4
            * (-1) ** j
            * sum(math.comb(2 * self.order, self.order + i) for i in range(j + 1, self.order + 1))
            / 2 ** (2 * self.order)
            for j in range(terms)
        )


def cks_chebyshev(system, config=None):
    config = config or CKSConfig()
    if not system.hermitian:
        raise ValidationError("CKS 稀疏输入需要 Hermitian 声明或显式 Hermitian dilation")
    a = real_symmetric_sparse_encoding(
        system.access,
        system.value_format,
        system.entry_bound,
        diagonal_nonnegative=system.diagonal_nonnegative,
    )
    coefficients = config.coefficients()
    inverse = lcu(
        [(weight, chebyshev_block(a, 2 * j + 1)) for j, weight in enumerate(coefficients)]
    )
    state = apply_be_to_state(inverse, system.rhs)
    return StateOracle(
        annotate(
            state.operation,
            "unitary",
            algorithm="cks_chebyshev_basic",
            input_model="sparse_location_inplace_and_entry_xor",
            input_alpha=a.alpha,
            encoded_inverse_bound=system.spectrum.inverse_bound(a.alpha),
            polynomial_order=config.order,
            polynomial_terms=len(coefficients),
            inverse_lcu_normalization=inverse.alpha,
            implementation_scope="CKS section 4 basic LCU; no VTAA",
            correctness="pending",
        )
    )


def make_cks_qlss(config=None):
    config = CKSConfig() if config is None else config
    require_instance(config, CKSConfig, "make_cks_qlss.config")
    return QLSSProtocol(
        "cks_chebyshev_basic", "sparse", lambda system: cks_chebyshev(system, config)
    )
