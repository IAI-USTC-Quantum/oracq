"""线性 QODE 的问题对象、可替换协议与 history 组装。"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace
from functools import partial

from oracq.algorithms.common.hamiltonian import taylor_hamiltonian
from oracq.algorithms.common.state_preparation import extend_initial, select_subspace
from oracq.algorithms.input_model.block_encoding import lcu, projector, tensor, truncated_shift
from oracq.algorithms.input_model.contracts import (
    AlgorithmContract,
    ContractIssue,
    ContractReport,
    InputRequirement,
    ResolvedInputs,
    finite_real,
    require_instance,
)
from oracq.algorithms.input_model.interfaces import (
    BlockEncodingProtocol,
    StateOracleProtocol,
    StatePreparationProtocol,
    as_block_encoding,
    as_state_preparation,
    operator_state_contract,
)
from oracq.algorithms.input_model.operators import BlockEncoding, identity, scale
from oracq.algorithms.input_model.oracles import StateOracle, StatePreparation
from oracq.algorithms.qode.cbmd import ContourPlan, cbmd_qode
from oracq.algorithms.qode.lchs import QuadraturePlan, lchs_qode
from oracq.algorithms.qode.ode_models import HermitianParts, LinearODE
from oracq.algorithms.qode.schrodingerization import SchrodingerPlan, schrodinger_qode
from oracq.infrastructure.ir import ValidationError


@dataclass(frozen=True)
class QODEProblem:
    """u'=Gu 的问题；dissipative 是数学声明，不是语言证明。"""

    generator: BlockEncoding
    initial: StatePreparation
    dissipative: bool | None = None
    initial_norm: float | None = None
    evidence: str = "caller_declared_unverified"

    def __init__(
        self, generator: BlockEncodingProtocol, initial: StatePreparationProtocol,
        dissipative: bool | None = None, initial_norm: float | None = None,
        evidence: str = "caller_declared_unverified",
    ) -> None:
        """提供方在问题构造边界归一，后续算法读取具体视图。"""
        for key, value in (
            ("generator", generator), ("initial", initial), ("dissipative", dissipative),
            ("initial_norm", initial_norm), ("evidence", evidence),
        ):
            object.__setattr__(self, key, value)
        self.__post_init__()

    def __post_init__(self) -> None:
        """把生成元与初态规范化为对应视图，并校验耗散声明、初值范数与说明字段。"""
        object.__setattr__(self, "generator", as_block_encoding(self.generator))
        object.__setattr__(self, "initial", as_state_preparation(self.initial))
        require_instance(self.initial, StatePreparation, "QODEProblem.initial")
        if self.dissipative is not None and type(self.dissipative) is not bool:
            raise ValidationError("dissipative 必须是 bool 或 None")
        if self.initial_norm is not None:
            finite_real(self.initial_norm, "QODEProblem.initial_norm", minimum=0, strict=True)
        if not isinstance(self.evidence, str) or not self.evidence:
            raise ValidationError("QODEProblem.evidence 需要非空说明")


@dataclass(frozen=True)
class QODESolver:
    """带契约检查的线性 QODE 求解协议包装。

    Attributes:
        provides: 输出声明满足的协议集合。
        name: 协议名，用于契约、属性与错误信息。
        kernel: 实际内核，形如 ``(generator, initial, time) -> StateOracle`` 的可调用。
        requires_dissipative: 为 True 时，问题级 ``solve`` 必须收到 ``dissipative=True`` 的声明。
    """

    provides = (StateOracleProtocol,)
    name: str
    kernel: Callable[..., StateOracle]
    requires_dissipative: bool = False

    def __post_init__(self) -> None:
        """校验协议名非空、内核可调用与耗散前提字段为 bool。"""
        if (
            not self.name
            or not callable(self.kernel)
            or type(self.requires_dissipative) is not bool
        ):
            raise ValidationError("QODESolver 需要名称、可调用内核与 bool 前提字段")

    @property
    def contract(self) -> AlgorithmContract:
        """返回该协议的算子-态契约。

        前提包含自治齐次 ``u'=Gu``；``requires_dissipative`` 为 True 时追加 ``Hermitian(G)<=0`` 的要求，Schrödingerization 协议还追加辅助窗口与恢复区域的验证责任。

        Returns:
            AlgorithmContract: 协议契约，``assumptions`` 为按协议定制的前提元组。
        """
        assumptions = ["自治齐次 u'=Gu；矩阵契约与数值近似待算法层核验"]
        if self.requires_dissipative:
            assumptions.append("Hermitian(G)<=0；问题级 solve 必须声明 dissipative=True")
        if self.name == "schrodingerization":
            assumptions.append("辅助窗口、Fourier 约定与所选恢复区域需要应用层验证")
        return replace(operator_state_contract(self.name), assumptions=tuple(assumptions))

    def check(
        self,
        generator: BlockEncodingProtocol | QODEProblem,
        initial: StatePreparationProtocol | None = None,
        time: float | None = None,
    ) -> ContractReport:
        """检查输入是否满足该协议的契约，并把未满足项汇集成报告。

        Args:
            generator: 生成元块编码；也可直接传入 ``QODEProblem``，此时改用问题对象内的生成元与初态。
            initial: 初态制备；传入 ``QODEProblem`` 时被忽略。
            time: 演化时长；给出时一并检查其为有限非负实数。

        Returns:
            ContractReport: 检查报告，``issues`` 汇集布局、调用能力与前提声明方面的全部未满足项。
        """
        return self._resolve(generator, initial, time).report

    def _resolve(
        self, generator: BlockEncodingProtocol | QODEProblem,
        initial: StatePreparationProtocol | None = None, time: float | None = None,
    ) -> ResolvedInputs:
        """保留本次取得的视图，并把问题级前提检查合并到同一报告。"""
        problem = generator if isinstance(generator, QODEProblem) else None
        if problem is not None:
            generator, initial = problem.generator, problem.initial
        resolved = self.contract.resolve(generator=generator, initial=initial)
        report = resolved.report
        issues = list(report.issues)
        if time is not None:
            try:
                finite_real(time, self.name + ".time", minimum=0)
            except ValidationError as exc:
                issues.append(
                    ContractIssue(
                        "CONFIG_VALUE",
                        self.name + ".time",
                        "finite time >= 0",
                        repr(time),
                        str(exc),
                    )
                )
        if problem is not None and self.requires_dissipative and problem.dissipative is not True:
            issues.append(
                ContractIssue(
                    "INPUT_PROMISE",
                    self.name + ".generator.dissipative",
                    True,
                    problem.dissipative,
                    "需要显式声明生成元耗散，或先做可追踪的整体移位",
                )
            )
        return replace(resolved, report=replace(report, issues=tuple(issues)))

    def _generate(
        self, generator: BlockEncoding, initial: StatePreparation, time: float
    ) -> StateOracle:
        """运行内核生成输出，并核验输出宽度与伴随、受控能力契约。"""
        result = self.kernel(generator, initial, time)
        require_instance(result, StateOracle, self.name + ".output")
        if result.width != generator.width:
            from oracq.algorithms.input_model.contracts import fail

            fail(
                "OUTPUT_LAYOUT",
                self.name + ".output",
                generator.width,
                result.width,
                "输出目标宽度必须与输入相同",
            )
        _, issues = InputRequirement(
            "output", (StateOracleProtocol,), adjoint=True, controlled=True
        ).inspect(result, prefix=self.name)
        if issues:
            from oracq.algorithms.input_model.contracts import ContractError

            raise ContractError(issues)
        return result

    def __call__(
        self, generator: BlockEncodingProtocol, initial: StatePreparationProtocol, time: float
    ) -> StateOracle:
        """兼容入口：数学前提仍由调用者承担；推荐问题级 solve。"""
        resolved = self._resolve(generator, initial, time)
        resolved.report.require()
        finite_real(time, self.name + ".time", minimum=0)
        return self._generate(
            resolved.get("generator", BlockEncoding), resolved.get("initial", StatePreparation), time
        )

    def solve(self, problem: QODEProblem, time: float) -> StateOracle:
        """问题级求解入口：检查契约后演化，并把声明来源写入输出属性。

        Args:
            problem: ``QODEProblem`` 问题对象。
            time: 演化时长，须为有限非负实数。

        Returns:
            StateOracle: 输出态 oracle；模块属性保存 ``evidence``、耗散声明与初值范数（若提供），不添加虚构的范数恢复能力。

        Raises:
            ContractError: 契约检查未通过（含耗散前提未声明）。
            ValidationError: ``problem`` 不是 ``QODEProblem``，或 ``time`` 不是有限非负实数。
        """
        require_instance(problem, QODEProblem, self.name + ".problem")
        resolved = self._resolve(problem, time=time)
        resolved.report.require()
        finite_real(time, self.name + ".time", minimum=0)
        result = self._generate(
            resolved.get("generator", BlockEncoding), resolved.get("initial", StatePreparation), time
        )
        # 保存声明来源；不添加虚构的物理范数恢复能力。
        attrs = dict(result.operation.module.attributes)
        attrs["qode_input_evidence"] = problem.evidence
        if problem.dissipative is not None:
            attrs["qode_dissipative_promise"] = problem.dissipative
        if problem.initial_norm is not None:
            attrs["qode_initial_norm"] = problem.initial_norm
        return StateOracle(
            replace(
                result.operation,
                module=replace(result.operation.module, attributes=tuple(sorted(attrs.items()))),
            )
        )


# 旧名称保留为同一类型的别名。
QODEProtocol = QODESolver


def make_euler_history_qode(
    qlss: Callable[[BlockEncoding, StatePreparation], StateOracle],
    *,
    steps: int = 2,
) -> Callable[[BlockEncoding, StatePreparation, float], StateOracle]:
    """把 QLSS 求解器包装成隐式 Euler 时间离散的 QODE 生成函数。

    在历史寄存器上装配一次求解全部时间步的线性系统：矩阵为 ``I - dt*(Q⊗G) - S``（Q 选择全部非零时刻，S 为截断的步进移位），右端把初态放在零号时刻，随后调用 ``qlss`` 求解并选取末时刻子空间。

    Args:
        qlss: 线性系统求解协议，接收矩阵与右端制备并返回 ``StateOracle``。
        steps: 时间步数，须为正整数。

    Returns:
        callable: 形如 ``(generator, initial, final_time) -> StateOracle`` 的生成函数，输出为末时刻的态。

    Raises:
        ValidationError: ``steps`` 不是正整数；生成阶段中输入布局或时间无效。
    """
    if type(steps) is not int or steps < 1:
        raise ValidationError("时间步数必须为正整数")

    def generate(
        generator: BlockEncoding, initial: StatePreparation, final_time: float
    ) -> StateOracle:
        """在历史寄存器上求解全部时间步并读出末时刻的态。"""
        if generator.width != initial.width or final_time <= 0:
            raise ValidationError("QODE 输入布局或时间无效")
        nt = steps.bit_length()
        dt = final_time / steps
        q = projector(nt, range(1, steps + 1))
        shift = truncated_shift(nt, steps)
        c = lcu(
            [
                (1, identity(nt + generator.width)),
                (-dt, tensor(q, generator)),
                (-1, tensor(shift, identity(generator.width))),
            ]
        )
        rhs = extend_initial(initial, nt)
        history = qlss(c, rhs)
        return select_subspace(history, generator.width, steps, label="qode_final_time")

    return generate


def linear_qode(
    method: str,
    *,
    hamiltonian_function: Callable[[BlockEncoding, float], BlockEncoding] = taylor_hamiltonian,
    **options: object,
) -> QODESolver:
    """通用 u'=Gu 接口，可直接注入既有 QHAM / make_qpde。

    Args:
        method: 求解方法，取 ``schrodingerization``、``lchs`` 或 ``cbmd``。
        hamiltonian_function: 形如 (BE, time) 返回 BlockEncoding 的哈密顿量模拟实现，
            缺省为截断 Taylor。
        **options: 方法配置；仅接受 ``plan``，类型须与所选方法匹配。

    Returns:
        QODESolver: 包装所选方法内核的求解协议，lchs/cbmd 路径要求耗散声明。
    """

    if not callable(hamiltonian_function):
        raise ValidationError("hamiltonian_function 必须可调用")
    if set(options) - {"plan"}:
        raise ValidationError("未知 QODE 配置：" + ", ".join(sorted(set(options) - {"plan"})))
    if method == "schrodingerization":
        if options.get("plan") is not None:
            require_instance(options["plan"], SchrodingerPlan, "schrodingerization.plan")
        return QODESolver(
            method,
            partial(
                schrodinger_qode,
                hamiltonian_function=hamiltonian_function,
                **options,  # type: ignore[arg-type]
            ),
        )
    if method not in {"lchs", "cbmd"}:
        raise ValidationError("未知线性 QODE 方法")
    algorithm = lchs_qode if method == "lchs" else cbmd_qode
    if options.get("plan") is not None:
        require_instance(
            options["plan"], QuadraturePlan if method == "lchs" else ContourPlan, method + ".plan"
        )

    def generate(generator: BlockEncoding, initial: StatePreparation, time: float) -> StateOracle:
        """构造 LinearODE 输入模型并调用对应算法内核。"""
        model = LinearODE(HermitianParts.from_operator(scale(-1, generator)), initial)
        return algorithm(
            model,
            time,
            hamiltonian_function=hamiltonian_function,
            **options,  # type: ignore[arg-type]
        )

    return QODESolver(method, generate, requires_dissipative=True)
