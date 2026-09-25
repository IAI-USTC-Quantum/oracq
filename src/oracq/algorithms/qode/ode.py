"""Problem objects, replaceable protocols, and history assembly for linear QODEs."""

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
    """The u'=Gu problem; dissipative is a mathematical declaration, not a language proof."""

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
        """Providers are normalized at the problem construction boundary; downstream algorithms read concrete views."""
        for key, value in (
            ("generator", generator), ("initial", initial), ("dissipative", dissipative),
            ("initial_norm", initial_norm), ("evidence", evidence),
        ):
            object.__setattr__(self, key, value)
        self.__post_init__()

    def __post_init__(self) -> None:
        """Normalize the generator and initial state into the corresponding views, and validate the dissipativity declaration, the initial norm, and the evidence field."""
        object.__setattr__(self, "generator", as_block_encoding(self.generator))
        object.__setattr__(self, "initial", as_state_preparation(self.initial))
        require_instance(self.initial, StatePreparation, "QODEProblem.initial")
        if self.dissipative is not None and type(self.dissipative) is not bool:
            raise ValidationError("dissipative must be bool or None")
        if self.initial_norm is not None:
            finite_real(self.initial_norm, "QODEProblem.initial_norm", minimum=0, strict=True)
        if not isinstance(self.evidence, str) or not self.evidence:
            raise ValidationError("QODEProblem.evidence requires a nonempty description")


@dataclass(frozen=True)
class QODESolver:
    """Linear QODE solve-protocol wrapper with contract checking.

    Attributes:
        provides: Set of protocols the output is declared to satisfy.
        name: Protocol name, used in contracts, attributes, and error messages.
        kernel: The actual kernel; a callable of the form ``(generator, initial, time) -> StateOracle``.
        requires_dissipative: When True, the problem-level ``solve`` must receive a ``dissipative=True`` declaration.
    """

    provides = (StateOracleProtocol,)
    name: str
    kernel: Callable[..., StateOracle]
    requires_dissipative: bool = False

    def __post_init__(self) -> None:
        """Validate that the protocol name is nonempty, the kernel is callable, and the dissipativity precondition field is bool."""
        if (
            not self.name
            or not callable(self.kernel)
            or type(self.requires_dissipative) is not bool
        ):
            raise ValidationError("QODESolver requires a name, a callable kernel, and a bool precondition field")

    @property
    def contract(self) -> AlgorithmContract:
        """Return the operator-state contract of this protocol.

        The preconditions include the autonomous homogeneous ``u'=Gu``; when ``requires_dissipative`` is True the requirement ``Hermitian(G)<=0`` is appended, and the Schrödingerization protocol additionally appends the validation responsibility for the auxiliary window and the recovery region.

        Returns:
            AlgorithmContract: The protocol contract; ``assumptions`` is the tuple of protocol-specific preconditions.
        """
        assumptions = ["autonomous homogeneous u'=Gu; the matrix contract and numerical approximations pending algorithm-layer validation"]
        if self.requires_dissipative:
            assumptions.append("Hermitian(G)<=0; the problem-level solve must declare dissipative=True")
        if self.name == "schrodingerization":
            assumptions.append("the auxiliary window, Fourier convention, and the chosen recovery region require application-layer validation")
        return replace(operator_state_contract(self.name), assumptions=tuple(assumptions))

    def check(
        self,
        generator: BlockEncodingProtocol | QODEProblem,
        initial: StatePreparationProtocol | None = None,
        time: float | None = None,
    ) -> ContractReport:
        """Check whether the inputs satisfy this protocol's contract, collecting unsatisfied items into a report.

        Args:
            generator: Generator block encoding; a ``QODEProblem`` may also be passed directly, in which case the generator and initial state inside the problem object are used instead.
            initial: Initial state preparation; ignored when a ``QODEProblem`` is passed.
            time: Evolution duration; when given, it is also checked to be a finite non-negative real number.

        Returns:
            ContractReport: The check report; ``issues`` collects all unsatisfied items regarding layout, invocation capabilities, and precondition declarations.
        """
        return self._resolve(generator, initial, time).report

    def _resolve(
        self, generator: BlockEncodingProtocol | QODEProblem,
        initial: StatePreparationProtocol | None = None, time: float | None = None,
    ) -> ResolvedInputs:
        """Keep the views obtained this time and merge the problem-level precondition checks into the same report."""
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
                    "An explicit generator dissipativity declaration is required, or a traceable global shift must be applied first",
                )
            )
        return replace(resolved, report=replace(report, issues=tuple(issues)))

    def _generate(
        self, generator: BlockEncoding, initial: StatePreparation, time: float
    ) -> StateOracle:
        """Run the kernel to generate the output, and verify the output width and the adjoint and controlled capability contracts."""
        result = self.kernel(generator, initial, time)
        require_instance(result, StateOracle, self.name + ".output")
        if result.width != generator.width:
            from oracq.algorithms.input_model.contracts import fail

            fail(
                "OUTPUT_LAYOUT",
                self.name + ".output",
                generator.width,
                result.width,
                "The output target width must equal the input width",
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
        """Compatibility entry: mathematical preconditions remain the caller's responsibility; the problem-level solve is recommended."""
        resolved = self._resolve(generator, initial, time)
        resolved.report.require()
        finite_real(time, self.name + ".time", minimum=0)
        return self._generate(
            resolved.get("generator", BlockEncoding), resolved.get("initial", StatePreparation), time
        )

    def solve(self, problem: QODEProblem, time: float) -> StateOracle:
        """Problem-level solve entry: check the contract, evolve, and write the declaration provenance into the output attributes.

        Args:
            problem: A ``QODEProblem`` problem object.
            time: Evolution duration; must be a finite non-negative real number.

        Returns:
            StateOracle: The output state oracle; module attributes store ``evidence``, the dissipativity declaration, and the initial norm when provided, without adding a fabricated norm recovery capability.

        Raises:
            ContractError: The contract check failed (including an undeclared dissipativity precondition).
            ValidationError: ``problem`` is not a ``QODEProblem``, or ``time`` is not a finite non-negative real number.
        """
        require_instance(problem, QODEProblem, self.name + ".problem")
        resolved = self._resolve(problem, time=time)
        resolved.report.require()
        finite_real(time, self.name + ".time", minimum=0)
        result = self._generate(
            resolved.get("generator", BlockEncoding), resolved.get("initial", StatePreparation), time
        )
        # Store declaration provenance; no fabricated physical norm recovery capability is added.
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


# The old name is retained as an alias of the same type.
QODEProtocol = QODESolver


def make_euler_history_qode(
    qlss: Callable[[BlockEncoding, StatePreparation], StateOracle],
    *,
    steps: int = 2,
) -> Callable[[BlockEncoding, StatePreparation, float], StateOracle]:
    """Wrap a QLSS solver into an implicit-Euler time-discretized QODE generation function.

    It assembles, on the history register, a single linear system solving all time steps at once: the matrix is ``I - dt*(Q⊗G) - S`` (Q selects all nonzero times, S is the truncated step shift), the right-hand side places the initial state at time zero, then ``qlss`` is invoked to solve and the final-time subspace is selected.

    Args:
        qlss: A linear-system solve protocol that takes a matrix and a right-hand-side preparation and returns a ``StateOracle``.
        steps: Number of time steps; must be a positive integer.

    Returns:
        callable: A generation function of the form ``(generator, initial, final_time) -> StateOracle`` whose output is the final-time state.

    Raises:
        ValidationError: ``steps`` is not a positive integer; or the input layout or time is invalid at generation time.
    """
    if type(steps) is not int or steps < 1:
        raise ValidationError("The number of time steps must be a positive integer")

    def generate(
        generator: BlockEncoding, initial: StatePreparation, final_time: float
    ) -> StateOracle:
        """Solve all time steps on the history register and read out the final-time state."""
        if generator.width != initial.width or final_time <= 0:
            raise ValidationError("Invalid QODE input layout or time")
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
    """Generic u'=Gu interface that can be injected directly into existing QHAM / make_qpde.

    Args:
        method: Solve method; one of ``schrodingerization``, ``lchs``, or ``cbmd``.
        hamiltonian_function: Hamiltonian simulation implementation of the form (BE, time)
            returning a BlockEncoding; defaults to truncated Taylor.
        **options: Method configuration; only ``plan`` is accepted, and its type must match the
            chosen method.

    Returns:
        QODESolver: The solve protocol wrapping the chosen method's kernel; the lchs/cbmd paths
        require a dissipativity declaration.
    """

    if not callable(hamiltonian_function):
        raise ValidationError("hamiltonian_function must be callable")
    if set(options) - {"plan"}:
        raise ValidationError("Unknown QODE options: " + ", ".join(sorted(set(options) - {"plan"})))
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
        raise ValidationError("Unknown linear QODE method")
    algorithm = lchs_qode if method == "lchs" else cbmd_qode
    if options.get("plan") is not None:
        require_instance(
            options["plan"], QuadraturePlan if method == "lchs" else ContourPlan, method + ".plan"
        )

    def generate(generator: BlockEncoding, initial: StatePreparation, time: float) -> StateOracle:
        """Build the LinearODE input model and invoke the corresponding algorithm kernel."""
        model = LinearODE(HermitianParts.from_operator(scale(-1, generator)), initial)
        return algorithm(
            model,
            time,
            hamiltonian_function=hamiltonian_function,
            **options,  # type: ignore[arg-type]
        )

    return QODESolver(method, generate, requires_dissipative=True)
