"""Linear-system input contracts, the Costa walk, and the CKS Chebyshev solver."""

from __future__ import annotations

import json
import math
from collections.abc import Callable, Iterable
from dataclasses import dataclass, replace
from typing import cast

from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.algorithms.common.state_preparation import apply_be_to_state, select_subspace
from oracq.algorithms.input_model.block_encoding import lcu, reflect_zero
from oracq.algorithms.input_model.contracts import (
    AlgorithmContract,
    ContractIssue,
    ContractReport,
    InputRequirement,
    finite_real,
    positive_integer,
    require_instance,
)
from oracq.algorithms.input_model.interfaces import (
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
from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.algorithms.input_model.oracles import (
    SparseAccess,
    StateOracle,
    StatePreparation,
    annotate,
    invoke,
    resources_for,
)
from oracq.algorithms.input_model.sparse import chebyshev_block, real_symmetric_sparse_encoding
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError, fuse


@dataclass(frozen=True)
class SpectralPromise:
    """Matrix spectral bounds declared by the caller; the language does not verify the declarations.

    Attributes:
        norm_upper: Declared upper bound on the matrix norm.
        sigma_min_lower: Declared lower bound on the minimum singular value; must be positive and not exceed norm_upper.
        evidence: Marker for the source of the declaration; the default denotes an unverified caller declaration.
    """

    norm_upper: float
    sigma_min_lower: float
    evidence: str = "caller_declared_unverified"

    def __post_init__(self) -> None:
        """Validate the declared values of the norm upper bound and singular value lower bound."""
        finite_real(self.norm_upper, "spectrum.norm_upper", minimum=0, strict=True)
        finite_real(self.sigma_min_lower, "spectrum.sigma_min_lower", minimum=0, strict=True)
        if not 0 < self.sigma_min_lower <= self.norm_upper or not math.isfinite(self.norm_upper):
            raise ValidationError("Finite positive norm upper bound and minimum singular value lower bound are required")

    def inverse_bound(self, alpha: float) -> float:
        """Derive an upper spectral bound for the encoded inverse operator from the declared singular value lower bound.

        Args:
            alpha: Block-encoding normalization factor actually consumed.

        Returns:
            float: ``max(1, alpha / sigma_min_lower)``.

        Raises:
            ValidationError: ``alpha`` is below the declared singular value lower bound.
        """
        if alpha < self.sigma_min_lower:
            raise ValidationError("BE alpha conflicts with the declared minimum singular value lower bound")
        return max(1.0, alpha / self.sigma_min_lower)


@dataclass(frozen=True)
class SparseSystem:
    """Sparse Hermitian problem over finite element words and a shared RHS; matrix properties are caller-declared."""

    access: SparseAccess
    value_format: FixedFormat
    entry_bound: float
    rhs: StatePreparation
    spectrum: SpectralPromise
    diagonal_nonnegative: bool = False
    hermitian: bool = False

    def __init__(
        self, access: CKSSparseProtocol, value_format: FixedFormat, entry_bound: float,
        rhs: StatePreparationProtocol, spectrum: SpectralPromise,
        diagonal_nonnegative: bool = False, hermitian: bool = False,
    ) -> None:
        """Accept provider inputs; after construction the fields always hold concrete views."""
        for key, value in (
            ("access", access), ("value_format", value_format), ("entry_bound", entry_bound),
            ("rhs", rhs), ("spectrum", spectrum),
            ("diagonal_nonnegative", diagonal_nonnegative), ("hermitian", hermitian),
        ):
            object.__setattr__(self, key, value)
        self.__post_init__()

    def __post_init__(self) -> None:
        """Normalize the sparse access and RHS views, and validate the layout declarations and the entry magnitude upper bound."""
        object.__setattr__(self, "access", as_sparse_access(self.access))
        object.__setattr__(self, "rhs", as_state_preparation(self.rhs))
        require_instance(self.value_format, FixedFormat, "SparseSystem.value_format")
        require_instance(self.rhs, StatePreparation, "SparseSystem.rhs")
        require_instance(self.spectrum, SpectralPromise, "SparseSystem.spectrum")
        finite_real(self.entry_bound, "SparseSystem.entry_bound", minimum=0, strict=True)
        if type(self.hermitian) is not bool or type(self.diagonal_nonnegative) is not bool:
            raise ValidationError("Sparse matrix property declarations must be bool")
        if (
            self.value_format.width != self.access.value_width
            or self.rhs.width != self.access.width
        ):
            raise ValidationError("Layouts of the sparse matrix, value format, and RHS do not match")
        if not math.isfinite(self.entry_bound) or self.entry_bound <= 0:
            raise ValidationError("The entry magnitude upper bound must be a finite positive number")


@dataclass(frozen=True)
class BlockSystem:
    """Linear problem under the block-encoding input model.

    Attributes:
        encoding: ``BlockEncoding`` of the matrix A.
        rhs: ``StatePreparation`` of the right-hand side b, with the same width as encoding.
        spectrum: Caller-declared spectral bounds.
    """

    encoding: BlockEncoding
    rhs: StatePreparation
    spectrum: SpectralPromise

    def __init__(
        self, encoding: BlockEncodingProtocol, rhs: StatePreparationProtocol,
        spectrum: SpectralPromise,
    ) -> None:
        """Inputs accept structural protocols; the public fields hold normalized views."""
        object.__setattr__(self, "encoding", encoding)
        object.__setattr__(self, "rhs", rhs)
        object.__setattr__(self, "spectrum", spectrum)
        self.__post_init__()

    def __post_init__(self) -> None:
        """Normalize the block encoding and RHS views, and check that both target widths agree."""
        object.__setattr__(self, "encoding", as_block_encoding(self.encoding))
        object.__setattr__(self, "rhs", as_state_preparation(self.rhs))
        require_instance(self.rhs, StatePreparation, "BlockSystem.rhs")
        require_instance(self.spectrum, SpectralPromise, "BlockSystem.spectrum")
        if self.encoding.width != self.rhs.width:
            raise ValidationError("BE and RHS widths do not match")

    @property
    def inverse_norm_bound(self) -> float:
        """Inverse spectral bound ``max(1, alpha / sigma_min_lower)`` of the encoded matrix A/alpha."""
        return self.spectrum.inverse_bound(self.encoding.alpha)


@dataclass(frozen=True)
class LinearSystem:
    """Problem-layer linear system: exactly one source input model plus readout and norm declarations.

    Attributes:
        sparse: Sparse source input model; exactly one of this and block.
        block: Block-encoding source input model; exactly one of this and sparse.
        physical_width: Width of the low target bits occupied by the physical solution; omitted means full width.
        physical_high_value: Integer value that the discarded high bits must equal when selecting the physical subspace.
        rhs_norm: Classical right-hand side norm; the physical magnitude of the solution can only be recovered when provided.
        data_assumptions: Data assumptions recorded alongside the solve result.
    """

    sparse: SparseSystem | None = None
    block: BlockSystem | None = None
    physical_width: int | None = None
    physical_high_value: int = 0
    rhs_norm: float | None = None
    data_assumptions: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Validate the exactly-one source input model plus the physical subspace and right-hand side norm declarations."""
        if (self.sparse is None) == (self.block is None):
            raise ValidationError("A linear problem requires exactly one source input model; other models are produced by explicit adaptation")
        source = cast(
            "BlockSystem | SparseSystem",
            self.block if self.block is not None else self.sparse,
        )
        require_instance(
            source, BlockSystem if self.block is not None else SparseSystem, "LinearSystem.source"
        )
        width = (
            cast("BlockSystem", source).encoding.width
            if self.block is not None
            else cast("SparseSystem", source).access.width
        )
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
            raise ValidationError("Invalid right-hand side norm declaration")

    def block_input(self) -> tuple[BlockSystem, tuple[str, ...]]:
        """Normalize the source input to a ``BlockSystem`` and report the adaptation trace.

        Returns:
            tuple[BlockSystem, tuple[str, ...]]: The block-encoding problem and a step-by-step adaptation description.

        Raises:
            ValidationError: The sparse source lacks a Hermitian declaration, so the CKS sparse adaptation cannot proceed.
        """
        if self.block is not None:
            return self.block, ("block input supplied",)
        from oracq.algorithms.input_model.sparse import real_symmetric_sparse_encoding

        s = cast("SparseSystem", self.sparse)
        if not s.hermitian:
            raise ValidationError(
                "The current sparse adaptation requires an explicit Hermitian declaration; non-Hermitian systems must first provide row and column access and dilation"
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
    """Protocol solve output: the physical-subspace solution state, an independent norm probe, and provenance.

    Attributes:
        state: Solution state oracle after selecting the physical subspace.
        norm_probe: Norm probe oracle on the branch where the solve and matrix probes jointly succeed.
        input_alpha: Block-encoding normalization factor actually consumed by the solve.
        encoded_inverse_bound: Upper spectral bound of the encoded inverse operator.
        rhs_norm: Classical right-hand side norm; when missing, the physical magnitude of the solution cannot be recovered.
        adapter_trace: Record of adaptation steps from the source input to the kernel.
        kernel_status: Kernel maturity declaration.
    """

    state: StateOracle
    norm_probe: StateOracle
    input_alpha: float
    encoded_inverse_bound: float
    rhs_norm: float | None
    adapter_trace: tuple[str, ...]
    kernel_status: str = "prototype; solver accuracy pending"

    def state_oracle(self) -> StateOracle:
        """Return the ``StateOracle`` of the physical-subspace solution state.

        Returns:
            StateOracle: Readout handle of the physical-subspace solution state.
        """
        return self.state

    @property
    def operation(self) -> Operation:
        """The RIR ``Operation`` corresponding to the solution state."""
        return self.state.operation

    def recover_norm(self, solver_success: float, joint_matrix_success: float) -> float:
        """Recover the norm of the solution vector from two success probabilities.

        Args:
            solver_success: Success probability of the solver probe alone.
            joint_matrix_success: Joint success probability of the solver and matrix probes.

        Returns:
            float: The right-hand side norm divided by the product of ``input_alpha`` and the
            square root of the joint-to-standalone success probability ratio.

        Raises:
            ValidationError: The classical right-hand side norm is missing, or the probability values or ordering are invalid.
        """
        if self.rhs_norm is None:
            raise ValidationError("Missing classical right-hand side norm; the physical update magnitude cannot be recovered")
        if not 0 < joint_matrix_success <= solver_success <= 1:
            raise ValidationError("Valid solver success probability and joint solver plus matrix probe success probability are required")
        return self.rhs_norm / (self.input_alpha * math.sqrt(joint_matrix_success / solver_success))


@dataclass(frozen=True)
class QLSSSolver:
    """Quantum linear-system protocol: contract validation, kernel dispatch, and assembly of a composable result.

    Attributes:
        name: Protocol name.
        input_model: Source input model, ``sparse`` or ``block_encoding``.
        kernel: Solve kernel; receives the source input model and returns a ``StateOracle``.
        legacy: Optional two-argument legacy entry; when omitted that call shape is unsupported.
    """

    provides = (StateOracleProtocol,)
    name: str
    input_model: str
    kernel: Callable[..., StateOracle]
    legacy: Callable[..., StateOracle] | None = None

    def __post_init__(self) -> None:
        """Validate the protocol name, the input model enumeration, and the callability of the kernel entries."""
        if not self.name or self.input_model not in {"sparse", "block_encoding"}:
            raise ValidationError("Invalid name or input_model for QLSSSolver")
        if not callable(self.kernel) or (self.legacy is not None and not callable(self.legacy)):
            raise ValidationError("QLSSSolver kernel must be callable")

    @property
    def contract(self) -> AlgorithmContract:
        """Return the ``AlgorithmContract`` for A and b according to the input model.

        The sparse model accepts only sparse oracles; the block-encoding model accepts both
        block-encoding and sparse interfaces.
        A and b must have the same width; spectral bounds and the matrix interpretation
        are caller-declared.
        """
        types = (
            (CKSSparseProtocol,)
            if self.input_model == "sparse"
            else (BlockEncodingProtocol, CKSSparseProtocol)
        )
        return AlgorithmContract(
            self.name,
            (
                InputRequirement(
                    "A", types,
                    adapter=as_sparse_access if self.input_model == "sparse" else as_qlss_matrix,
                    adjoint=True, controlled=True,
                ),
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
            assumptions=("spectral bounds and matrix interpretation are caller-declared; algorithm accuracy and success channel pending validation",),
        )

    def check(self, problem: object) -> ContractReport:
        """Run neither the kernel nor the adapters; produce a full report on structure, capabilities, and required adaptations.

        Args:
            problem: The linear-system problem to check; the valid type is ``LinearSystem``.

        Returns:
            ContractReport: A contract report summarizing the input type, a zero right-hand side,
            spectral declaration conflicts, and the required adaptation steps.
        """
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
                        "A LinearSystem is required",
                    ),
                ),
            )
        source = cast(
            "BlockSystem | SparseSystem",
            problem.block if problem.block is not None else problem.sparse,
        )
        a = (
            cast("BlockSystem", source).encoding
            if problem.block is not None
            else cast("SparseSystem", source).access
        )
        report = self.contract.check(A=a, b=source.rhs)
        issues: list[ContractIssue] = list(report.issues)
        adapters: list[str] = []
        if problem.rhs_norm == 0:
            issues.append(
                ContractIssue(
                    "INPUT_ZERO_RHS", self.name + ".b", "nonzero vector", 0, "A zero right-hand side should be handled on the classical side"
                )
            )
        if self.input_model == "sparse" and problem.sparse is None:
            issues.append(
                ContractIssue(
                    "INPUT_ADAPTER",
                    self.name + ".A",
                    "cks_sparse",
                    "block_encoding",
                    "Cannot automatically recover an efficient sparse oracle from a general BE",
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
                            "The current CKS sparse adaptation requires Hermitian and non-negative diagonal declarations",
                        )
                    )
            adapters.append("CKS real symmetric sparse -> Tdag S T block encoding")
            alpha = (
                cast("SparseSystem", source).access.sparsity
                * cast("SparseSystem", source).entry_bound
            )
        else:
            alpha = cast("BlockSystem", source).encoding.alpha
        if alpha < source.spectrum.sigma_min_lower:
            issues.append(
                ContractIssue(
                    "INPUT_SPECTRUM",
                    self.name + ".A.alpha",
                    f">={source.spectrum.sigma_min_lower}",
                    alpha,
                    "BE alpha conflicts with the declared minimum singular value lower bound",
                )
            )
        return replace(report, issues=tuple(issues), adapters=tuple(adapters))

    def __call__(self, *args: object) -> SolveResult | StateOracle:
        """Dispatch the protocol entry according to the call shape.

        A single ``LinearSystem`` argument goes to ``solve``; the two-argument shape is first
        checked against the operator-state contract of ``legacy``, then handed to the legacy
        entry.

        Args:
            *args: A single ``LinearSystem``, or two arguments giving the block-encoding
                matrix and the RHS preparation (the latter requires ``legacy`` to be set).

        Returns:
            SolveResult | StateOracle: A ``SolveResult`` for the single-argument shape, or the
                solution state output by the kernel for the legacy two-argument shape.

        Raises:
            ValidationError: The call shape is unsupported, or the legacy entry contract is not satisfied.
        """
        if len(args) == 1 and isinstance(args[0], LinearSystem):
            return self.solve(args[0])
        if len(args) == 2 and self.legacy is not None:
            resolved = operator_state_contract(self.name, matrix="A", state="b").resolve(
                A=args[0], b=args[1]
            )
            resolved.report.require()
            return self.legacy(
                resolved.get("A", BlockEncoding), resolved.get("b", StatePreparation),
            )
        raise ValidationError("QLSS protocol requires a LinearSystem; the input model cannot be guessed from the call shape")

    def solve(self, problem: LinearSystem) -> SolveResult:
        """Solve a ``LinearSystem`` and assemble the ``SolveResult``.

        The contract check runs first, then the source input goes to the kernel according to
        input_model; afterwards the physical subspace is selected, protocol provenance
        attributes are written, and the independent matrix norm probe is assembled.

        Args:
            problem: The linear system to solve.

        Returns:
            SolveResult: The solution state, norm probe, and provenance.

        Raises:
            ValidationError: The contract is not satisfied, the right-hand side is zero, the kernel output is invalid, or widths mismatch.
        """
        self.check(problem).require()
        if problem.rhs_norm == 0:
            raise ValidationError("A zero right-hand side should return a zero update on the classical side; no normalized RHS state exists")
        if self.input_model == "sparse":
            if problem.sparse is None:
                raise ValidationError("Cannot automatically recover an efficient sparse oracle from a general BE")
            source = problem.sparse
            state = self.kernel(source)
            block, trace = problem.block_input()
            trace = ("sparse input consumed by " + self.name,) + trace
        elif self.input_model == "block_encoding":
            block, trace = problem.block_input()
            state = self.kernel(block)
            trace = trace + (self.name + " consumes BE",)
        else:
            raise ValidationError("Unknown QLSS input_model")
        require_instance(state, StateOracle, self.name + ".output")
        from oracq.infrastructure.linking import capabilities

        capability = capabilities(state.operation.program())
        if not capability["supports_adjoint"] or not capability["supports_controlled"]:
            raise ValidationError("The composable QLSS output of QFVM requires adjoint and controlled capabilities")
        if state.width != block.encoding.width:
            raise ValidationError("QLSS output layout does not match the solved matrix")
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
        # Independent matrix norm probe: the Costa filtering success rate cannot be reused for the CKS LCU factor.
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


# The old name is retained as an alias of the same type.
QLSSProtocol = QLSSSolver


@dataclass(frozen=True)
class CostaConfig:
    """Schedule and Dolph–Chebyshev filter configuration for the Costa parameterized quantum walk.

    kappa is the inverse spectral bound of the actually encoded matrix A/alpha, not
    cond(A) at an arbitrary scale; the problem layer overrides it automatically using
    the system's ``inverse_norm_bound``.

    Attributes:
        steps: Number of chained walk steps; step i takes the schedule point (i+1)/steps.
        kappa: Inverse spectral bound of the encoded matrix A/alpha, i.e. the assumption
            sigma_min(A/alpha) >= 1/kappa; at least 1.
        schedule_power: Exponent of the interpolation schedule; cannot equal 1.
        filter_degree: Degree of the Dolph–Chebyshev filter; an even number of at least 2.
        filter_attenuation: Sidelobe level of the filter; smaller means stronger stopband
            suppression; lies strictly between 0 and 1.
    """

    steps: int = 2
    kappa: float = 4.0
    schedule_power: float = 1.5
    filter_degree: int = 2
    filter_attenuation: float = 0.2

    def __post_init__(self) -> None:
        """Validate the walk step count, schedule exponent, and the filter degree and sidelobe ranges."""
        positive_integer(self.steps, "CostaConfig.steps")
        finite_real(self.kappa, "CostaConfig.kappa", minimum=1)
        finite_real(self.schedule_power, "CostaConfig.schedule_power")
        if self.schedule_power == 1:
            raise ValidationError("Costa schedule_power cannot equal 1")
        positive_integer(self.filter_degree, "CostaConfig.filter_degree", minimum=2)
        if self.filter_degree % 2:
            raise ValidationError("Costa filter_degree must be even")
        finite_real(
            self.filter_attenuation, "CostaConfig.filter_attenuation", minimum=0, strict=True
        )
        if self.filter_attenuation >= 1:
            raise ValidationError("Costa filter_attenuation must be less than 1")


@dataclass(frozen=True)
class FilterPlan:
    """LCU description of a Laurent polynomial filter over walk powers.

    Attributes:
        weights: LCU weights for each walk power.
        stride: Power spacing between adjacent terms; at least 1.
        offset: Starting power; when negative, inverse powers are realized via the adjoint.
        method: Construction method marker.
    """

    weights: tuple[float, ...]
    stride: int = 1
    offset: int = 0
    method: str = "explicit"


def dolph_chebyshev_plan(degree: int = 2, attenuation: float = 0.2) -> FilterPlan:
    """Construct an even-degree Dolph–Chebyshev LCU via Laurent polynomial recursion.

    Args:
        degree: Filter degree; an even number of at least 2.
        attenuation: Stopband attenuation, in (0,1).

    Returns:
        FilterPlan: A walk-power filter plan starting at power −degree with power spacing 2.
    """
    if type(degree) is not int or degree < 2 or degree % 2:
        raise ValidationError("The Dolph–Chebyshev prototype requires a positive even degree")
    if not 0 < attenuation < 1:
        raise ValidationError("attenuation must lie between 0 and 1")
    beta = math.cosh(math.acosh(1 / attenuation) / degree)
    previous, current = {0: 1.0}, {-1: beta / 2, 1: beta / 2}
    for _ in range(2, degree + 1):
        following: dict[int, float] = {}
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


def schedule(s: float, kappa: float, power: float = 1.5) -> float:
    """Evaluate the Costa interpolation schedule at schedule point s.

    Args:
        s: Schedule point, between 0 and 1.
        kappa: Inverse spectral bound of the encoded matrix; at least 1.
        power: Schedule exponent; cannot equal 1.

    Returns:
        float: The interpolation parameter f(s); degenerates to s when kappa is 1.

    Raises:
        ValidationError: kappa is below 1, power equals 1, or s is out of range.
    """
    if kappa < 1 or power == 1 or not 0 <= s <= 1:
        raise ValidationError("Invalid Costa schedule parameters")
    if kappa == 1:
        return s
    return kappa / (kappa - 1) * (1 - (1 + s * (kappa ** (power - 1) - 1)) ** (1 / (1 - power)))


def costa_walk(a: BlockEncoding, bprep: StatePreparation, fs: float) -> Operation:
    """Assemble the single-step operator of the Costa parameterized quantum walk.

    It builds in order the RHS zero-state reflection, the schedule rotation, and the
    controlled forward and inverse block-encoding invocations, and finishes with a
    positive reflection on all signal bits.

    Args:
        a: Block encoding of the matrix A.
        bprep: State preparation of the right-hand side b, with the same width as a.
        fs: Interpolation parameter at the schedule point, between 0 and 1.

    Returns:
        Operation: The walk operator with a target/signal interface; signal consists of the
        BE signal bits, the RHS work bits, and four walk ancilla bits, in that order.

    Raises:
        ValidationError: The input contract is not satisfied, widths mismatch, or the schedule point is out of range.
    """
    operator_state_contract("costa_walk", matrix="A", state="b").check(A=a, b=bprep).require()
    finite_real(fs, "costa_walk.fs", minimum=0)
    a, bprep = as_block_encoding(a), as_state_preparation(bprep)
    if a.width != bprep.width or not 0 <= fs <= 1:
        raise ValidationError("Costa input shape or schedule point mismatch")
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

    def prep(inverse: bool = False) -> None:
        """Invoke the RHS state preparation; apply its adjoint when ``inverse`` is true."""
        if inverse:
            with b.adjoint():
                invoke(b, bprep.operation, "b", target=b["target"], work=bw)
        else:
            invoke(b, bprep.operation, "b", target=b["target"], work=bw)

    def reflect_rhs_input_zero() -> None:
        """Flip the phase of the branch where target and work bits are all zero, realizing the RHS input zero-state reflection."""
        # U_b is a unitary extension over target+work; the projection must also require the work bits to be zero.
        with b.control(b["target"], 0):
            if bw.width:
                with b.control(bw, 0):
                    b.global_phase(math.pi)
            else:
                b.global_phase(math.pi)

    def rotation() -> None:
        """Apply the rotation R(s) at schedule point ``fs`` on the walk ancilla bit."""
        # R(s) is a reflection matrix, written as Ry(2 atan2(f,1-f)) Z.
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


def unary_weight_preparation(weights: Iterable[float]) -> Operation:
    """Prepare an LCU weight superposition over a unary prefix chain.

    Args:
        weights: Weight sequence; entries are non-negative, the sum is positive, and there are at least two entries.

    Returns:
        Operation: A zero-input state preparation on the target register; the prefix state
        with the first k bits one and the rest zero has probability proportional to
        ``weights[k]``.

    Raises:
        ValidationError: The weights are empty, contain negative or non-finite values, sum to zero, or number fewer than two.
    """
    weights = tuple(float(v) for v in weights)
    if not weights or any(not math.isfinite(v) or v < 0 for v in weights) or sum(weights) <= 0:
        raise ValidationError("unary PREP requires non-negative weights with a nonzero sum")
    width = len(weights) - 1
    if width < 1:
        raise ValidationError("unary PREP requires at least two weights")
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


def lcu_filter(walk: Operation, plan: FilterPlan) -> Operation:
    """Coherently superpose a Laurent polynomial filter over walk powers.

    The weight ``plan.weights[k]`` acts on walk power ``plan.offset + k * plan.stride``;
    negative powers are realized via the adjoint of the walk operator.

    Args:
        walk: Walk operator exposing only a target/signal interface.
        plan: Weight and power configuration of the polynomial.

    Returns:
        Operation: A target/signal interface; signal appends clock bits beyond the walk signal bits.

    Raises:
        ValidationError: The walk register interface does not match, or the filter power configuration is invalid.
    """
    widths = {r.name: r.type.width for r in walk.module.registers}
    if set(widths) != {"target", "signal"}:
        raise ValidationError("filtering expects a target/signal walk interface")
    if type(plan.stride) is not int or plan.stride < 1 or type(plan.offset) is not int:
        raise ValidationError("Invalid filter power configuration")
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

    def repeat_walk(count: int) -> None:
        """Invoke the walk operator ``count`` times in a row."""
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


def costa_qlss(
    a: BlockEncoding,
    bprep: StatePreparation,
    config: CostaConfig | None = None,
    *,
    filtering: FilterPlan | None = None,
) -> StateOracle:
    """Costa walk solve kernel: prepare the RHS, chain the walks at each schedule point, and apply the filter.

    Args:
        a: Block encoding of the matrix A.
        bprep: State preparation of the right-hand side b, with the same width as a.
        config: Walk schedule and filter configuration; defaults to the ``CostaConfig`` defaults when omitted.
        filtering: Explicit filter plan; when omitted, a Dolph–Chebyshev plan is generated from the configuration.

    Returns:
        StateOracle: The solve circuit; the success condition is all-zero signal, with probability and accuracy unverified.

    Raises:
        ValidationError: steps is not a positive integer, or validation of an underlying building block failed.
    """
    a, bprep = as_block_encoding(a), as_state_preparation(bprep)
    config = CostaConfig() if config is None else config
    if type(config.steps) is not int or config.steps < 1:
        raise ValidationError("Costa steps must be a positive integer")
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


def make_costa_qlss(config: CostaConfig | None = None) -> QLSSSolver:
    """Declare BE input; the problem layer automatically derives the actual schedule parameters from alpha/sigma_min.

    Args:
        config: Costa walk solve configuration; defaults to the default configuration.

    Returns:
        QLSSSolver: The Costa solve protocol under the block-encoding input model, with a two-argument legacy entry.
    """
    from dataclasses import replace

    config = CostaConfig() if config is None else config
    require_instance(config, CostaConfig, "make_costa_qlss.config")

    def kernel(system: BlockSystem) -> StateOracle:
        """Override kappa with the system's inverse spectral bound and invoke the Costa kernel."""
        effective = replace(config, kappa=system.inverse_norm_bound)
        return costa_qlss(system.encoding, system.rhs, effective)

    return QLSSSolver(
        "costa_general_walk", "block_encoding", kernel, legacy=lambda a, b: costa_qlss(a, b, config)
    )


@dataclass(frozen=True)
class CKSConfig:
    """Truncation configuration of the CKS basic Chebyshev inverse-operator expansion.

    Attributes:
        order: Order of the Chebyshev expansion, in the range 1..128.
        terms: Number of truncated terms actually kept; defaults to order when omitted.
    """

    order: int = 2
    terms: int | None = None

    def __post_init__(self) -> None:
        """Validate the expansion order and the optional truncation term count."""
        positive_integer(self.order, "CKSConfig.order", maximum=128)
        if self.terms is not None:
            positive_integer(self.terms, "CKSConfig.terms", maximum=self.order)

    def coefficients(self) -> tuple[float, ...]:
        """Compute the LCU weights of the truncated Chebyshev expansion of 1/x.

        Returns:
            tuple[float, ...]: The truncated coefficient sequence; term j pairs with the
            odd power ``T_(2j+1)``.

        Raises:
            ValidationError: order is outside 1..128, or the truncation term count is invalid.
        """
        if type(self.order) is not int or not 1 <= self.order <= 128:
            raise ValidationError("The order of the CKS basic prototype must lie in 1..128")
        terms = self.order if self.terms is None else self.terms
        if not 1 <= terms <= self.order:
            raise ValidationError("Invalid CKS truncation term count")
        return tuple(
            4
            * (-1) ** j
            * sum(math.comb(2 * self.order, self.order + i) for i in range(j + 1, self.order + 1))
            / 2 ** (2 * self.order)
            for j in range(terms)
        )


def cks_chebyshev(system: SparseSystem, config: CKSConfig | None = None) -> StateOracle:
    """CKS basic Chebyshev solve kernel, accepting the sparse Hermitian input model.

    It applies the truncated Chebyshev expansion of 1/x to the sparse block encoding,
    then applies the inverse-operator LCU to the RHS preparation.

    Args:
        system: The sparse Hermitian linear system described by a ``SparseSystem``.
        config: Expansion order and truncation configuration; defaults to the ``CKSConfig`` defaults when omitted.

    Returns:
        StateOracle: The solve circuit; basic LCU only, without VTAA, with prototype-level accuracy claims.

    Raises:
        ValidationError: The input lacks a Hermitian declaration.
    """
    config = config or CKSConfig()
    if not system.hermitian:
        raise ValidationError("CKS sparse input requires a Hermitian declaration or an explicit Hermitian dilation")
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


def make_cks_qlss(config: CKSConfig | None = None) -> QLSSSolver:
    """Construct the CKS solve protocol for sparse input.

    Args:
        config: The ``CKSConfig`` passed to the kernel; defaults when omitted.

    Returns:
        QLSSSolver: A solve protocol with input_model ``sparse``.
    """
    config = CKSConfig() if config is None else config
    require_instance(config, CKSConfig, "make_cks_qlss.config")
    return QLSSSolver(
        "cks_chebyshev_basic", "sparse", lambda system: cks_chebyshev(system, config)
    )
