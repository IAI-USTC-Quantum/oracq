"""QLSS 的问题契约、显式输入适配及统一物理输出/范数探针。"""

import json
import math
from dataclasses import dataclass, replace

from .algorithms.solvers import select_subspace
from .builder import Builder
from .ir import Bits, ValidationError
from .library import _name
from .oracles import StateOracle, annotate, invoke, resources_for


@dataclass(frozen=True)
class SpectralPromise:
    norm_upper: float
    sigma_min_lower: float
    evidence: str = "caller_declared_unverified"

    def __post_init__(self):
        if not 0 < self.sigma_min_lower <= self.norm_upper or not math.isfinite(self.norm_upper):
            raise ValidationError("需要有限正的范数上界及最小奇异值下界")

    def inverse_bound(self, alpha):
        if alpha < self.sigma_min_lower:
            raise ValidationError("BE alpha 与奇异值下界声明冲突")
        return max(1.0, alpha / self.sigma_min_lower)


@dataclass(frozen=True)
class SparseSystem:
    """有限元素字与同一 RHS 的稀疏 Hermitian 问题；矩阵性质由调用者声明。"""

    access: object
    value_format: object
    entry_bound: float
    rhs: object
    spectrum: SpectralPromise
    diagonal_nonnegative: bool = False
    hermitian: bool = False

    def __post_init__(self):
        if (
            self.value_format.width != self.access.value_width
            or self.rhs.width != self.access.width
        ):
            raise ValidationError("稀疏矩阵、值格式与 RHS 的布局不匹配")
        if not math.isfinite(self.entry_bound) or self.entry_bound <= 0:
            raise ValidationError("元素幅值上界必须为有限正数")


@dataclass(frozen=True)
class BlockSystem:
    encoding: object
    rhs: object
    spectrum: SpectralPromise

    def __post_init__(self):
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
        if self.rhs_norm is not None and (self.rhs_norm < 0 or not math.isfinite(self.rhs_norm)):
            raise ValidationError("右端范数声明无效")

    def block_input(self):
        if self.block is not None:
            return self.block, ("block input supplied",)
        from .sparse_models import real_symmetric_sparse_encoding

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
    name: str
    input_model: str
    kernel: object
    legacy: object = None

    def __call__(self, *args):
        if len(args) == 1 and isinstance(args[0], LinearSystem):
            return self.solve(args[0])
        if len(args) == 2 and self.legacy is not None:
            return self.legacy(*args)
        raise ValidationError("QLSS protocol 需要 LinearSystem；不能根据调用形状猜测输入模型")

    def solve(self, problem):
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
        from .linking import capabilities

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
