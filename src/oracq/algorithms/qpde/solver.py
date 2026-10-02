"""Prepared quantum PDE instances, backend execution, and physical readout."""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from functools import partial
from typing import Protocol

from oracq.algorithms.common.solve import SolveCircuit, resolve_memory
from oracq.algorithms.input_model.contracts import (
    AlgorithmContract,
    ContractReport,
    InputRequirement,
    positive_integer,
    require_instance,
)
from oracq.algorithms.input_model.interfaces import (
    BlockEncodingProtocol,
    StatePreparationProtocol,
    as_block_encoding,
    as_state_preparation,
)
from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import StateOracle, StatePreparation, annotate
from oracq.algorithms.input_model.pde import InitialLayout, PDEProblem, PreparedInitial
from oracq.algorithms.input_model.qham import qham_input_model, taylor_qode
from oracq.applications.qham.linearization import QHAMPlan
from oracq.applications.qham.reference import Discretization
from oracq.applications.qham.stencils import structured_fd_bindings
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.execution import RegisterState
from oracq.infrastructure.ir import Bits, ValidationError
from oracq.infrastructure.native import NativeRegistry
from oracq.infrastructure.qram_schema import RegisteredQRAM
from oracq.infrastructure.validation import validate


@dataclass(frozen=True)
class QHAMConfig:
    """HAM truncation and independent linear evolution configuration.

    The default linear solver is a finite Taylor polynomial of the lifted
    generator. HAM convergence, spatial error, angle quantization, and Taylor
    error remain separate approximation questions. Supply ``linear_solver``
    to use an existing QODE protocol with its own input requirements.
    """

    order: int = 1
    eta: complex = -1.0
    taylor_degree: int = 1
    linear_solver: Callable[[BlockEncoding, StatePreparation, float], StateOracle] | None = None
    max_blocks: int = 256
    max_terms: int = 4096
    max_initial_words: int = 4096
    max_coefficient_words: int = 4096

    def __post_init__(self) -> None:
        positive_integer(self.order, "QHAM.order", minimum=0)
        positive_integer(self.taylor_degree, "QHAM.taylor_degree", minimum=0)
        for key in ("max_blocks", "max_terms", "max_initial_words", "max_coefficient_words"):
            positive_integer(getattr(self, key), "QHAM." + key)
        if not math.isfinite(complex(self.eta).real) or not math.isfinite(complex(self.eta).imag):
            raise ValidationError("QHAM eta must be finite")
        if self.linear_solver is not None and not callable(self.linear_solver):
            raise ValidationError("QHAM linear_solver must be a three-argument QODE protocol")


@dataclass(frozen=True)
class PreparedPDE:
    """A method's generated state oracle and its validated input contract.

    ``amplitude_scale`` recovers physical values from unnormalized success
    amplitudes when the chosen method supplies this normalization information.
    Custom methods may leave it unresolved and still support state readout.
    """

    state: StateOracle
    report: ContractReport
    qrams: tuple[RegisteredQRAM, ...] = ()
    model: object | None = None
    amplitude_scale: float | None = None


class PDEMethod(Protocol):
    """An algorithm that owns the contracts and adaptations for its PDE inputs.

    Implementations check their own grid, oracle capabilities, and mathematical
    promises, returning a physical ``target/signal`` state oracle. There is no
    implicit conversion from a raw value query into amplitude preparation.
    """

    def prepare(self, problem: PDEProblem) -> PreparedPDE:
        """Adapt the problem's inputs and return the validated state oracle."""


@dataclass(frozen=True)
class QHAMMethod:
    """QHAM adaptation from independently chosen grid and initial inputs."""

    config: QHAMConfig = QHAMConfig()

    def prepare(self, problem: PDEProblem) -> PreparedPDE:
        """Check oracle contracts, lift HAM, and generate the selected QODE circuit."""
        cfg = self.config
        layout = InitialLayout(problem.grid, problem.equation.fields, cfg.max_initial_words)
        initial = require_instance(problem.initial.prepare_initial(layout), PreparedInitial, "QHAM.initial")
        initial_requirement = InputRequirement(
            "initial", (StatePreparationProtocol,), adapter=as_state_preparation,
            main_qubit=layout.width, adjoint=True, controlled=True,
            zero_input=True, clean_work=True,
        )
        # Reject a bad initial view before spending the spatial generation budget.
        AlgorithmContract("QHAM", (initial_requirement,)).check(initial=initial.preparation).require()
        plan = QHAMPlan(problem.equation, cfg.order)
        if initial.norm == 0 and not plan.has_forcing:
            b = Builder("zero_pde_solution", {"target": Bits(layout.width), "signal": Bits(1)})
            b.x(b["signal"])
            state = StateOracle(annotate(b.finish(), "unitary", algorithm="zero_pde_solution"))
            report = AlgorithmContract("QHAM", (initial_requirement,)).check(initial=initial.preparation)
            return PreparedPDE(state, report, amplitude_scale=1.0)
        space = Discretization(problem.equation, problem.grid, problem.known)
        bindings = structured_fd_bindings(
            space, initial.preparation, initial_norm=initial.norm,
            derivative_encoder=problem.grid.derivative_encoding,
            max_coefficient_words=cfg.max_coefficient_words,
        ).validate(plan)
        requirements = (initial_requirement,) + tuple(InputRequirement(
            name, (BlockEncodingProtocol,), adapter=as_block_encoding,
            main_qubit=max(1, port.arity) * layout.width, adjoint=True, controlled=True,
        ) for name, port in bindings.ports)
        report = AlgorithmContract(
            "QHAM", requirements,
            assumptions=("The grid encodings implement the declared discrete PDE ports.",
                         "HAM convergence and the linear solver approximation need separate error validation."),
        ).check(initial=initial.preparation, **{name: port.encoding for name, port in bindings.ports}).require()
        model = qham_input_model(plan, bindings, eta=cfg.eta, max_blocks=cfg.max_blocks, max_terms=cfg.max_terms)
        linear_solver = cfg.linear_solver or partial(taylor_qode, degree=cfg.taylor_degree)
        state = model.solve(linear_solver, problem.final_time)
        scale = None
        if cfg.linear_solver is None:
            # This factor is the exact LCU normalization of the selected finite
            # polynomial, not a convergence guarantee for the underlying PDE.
            x = problem.final_time * model.generator.alpha
            term, alpha = 1.0, 1.0
            for k in range(1, cfg.taylor_degree + 1):
                term *= x / k
                alpha += term
            try:
                scale = math.exp(model.log_initial_norm) * alpha
            except OverflowError:
                scale = None
            if scale is not None and not math.isfinite(scale):
                scale = None
        return PreparedPDE(state, report, (*problem.grid.qrams, *initial.qrams), model, scale)


@dataclass(frozen=True)
class PDECircuit(SolveCircuit):
    """A prepared PDE circuit with shared export, estimation, and execution."""


@dataclass(frozen=True)
class PDEResult:
    """Backend state and physical success-channel readout of a prepared PDE."""

    state: RegisterState
    layout: InitialLayout
    backend: str
    amplitude_scale: float | None = None

    def _physical(self) -> dict[str, tuple[complex, ...]]:
        names = [r.name for r in self.state.registers]
        target, signal = names.index("target"), names.index("signal")
        amplitudes = {values[target]: amplitude for values, amplitude in self.state.amplitudes.items() if values[signal] == 0}
        storage = 1 << self.layout.grid.spatial_width
        return {name: tuple(amplitudes.get(i * storage + node, 0j) for node in range(self.layout.grid.size)) for i, name in enumerate(self.layout.fields)}

    @property
    def success_probability(self) -> float:
        """Joint probability of zero signal and a physical output address."""
        return sum(abs(value) ** 2 for values in self._physical().values() for value in values)

    def physical_amplitudes(self, *, normalize: bool = True) -> dict[str, tuple[complex, ...]]:
        """Read physical samples; optionally normalize the joint success channel."""
        values = self._physical()
        if not normalize:
            return values
        probability = self.success_probability
        if probability <= 0:
            raise ValidationError("the physical success channel has zero probability")
        return {name: tuple(v / math.sqrt(probability) for v in row) for name, row in values.items()}

    def physical_values(self) -> dict[str, tuple[complex, ...]]:
        """Recover physical values of the generated approximation when its scale is known."""
        if self.amplitude_scale is None:
            raise ValidationError("this linear solver did not supply physical magnitude recovery; read conditional amplitudes instead")
        return {name: tuple(v * self.amplitude_scale for v in row) for name, row in self._physical().items()}


class PDESolveInstance:
    """Explicit prepare/inspect/execute lifecycle for a quantum PDE algorithm.

    Preparation is cached and imports no optional quantum backend. Execution
    requires successful preparation and never silently regenerates a circuit.
    """

    def __init__(self, problem: PDEProblem, method: PDEMethod) -> None:
        self._problem = problem
        self._method = method
        self._prepared: PreparedPDE | None = None
        self._circuit: PDECircuit | None = None

    @property
    def problem(self) -> PDEProblem:
        """The problem whose configuration is captured by this instance."""
        return self._problem

    @property
    def method(self) -> PDEMethod:
        """The algorithm implementation selected for this instance."""
        return self._method

    @property
    def prepared(self) -> PreparedPDE:
        """Inspect the method's model and checked contract after preparation."""
        if self._prepared is None:
            raise ValidationError("call prepare() before inspecting or running the PDE instance")
        return self._prepared

    def prepare(self) -> PDESolveInstance:
        """Generate and validate once, then resolve every entry QRAM snapshot."""
        if self._prepared is None:
            prepared = require_instance(self.method.prepare(self.problem), PreparedPDE, "PDE.method")
            prepared.report.require()
            layout = InitialLayout(self.problem.grid, self.problem.equation.fields)
            if prepared.state.width != layout.width:
                raise ValidationError("a PDE method must return the physical target register width")
            program = validate(prepared.state.operation.program(), require_closed=True)
            memory = resolve_memory(program, prepared.qrams)
            circuit = PDECircuit(program, tuple((name, tuple(sorted(cells.items()))) for name, cells in sorted(memory.items())))
            self._prepared, self._circuit = prepared, circuit
        return self

    def circuit(self) -> PDECircuit:
        """Inspect or export the generated circuit after preparation."""
        if self._circuit is None:
            raise ValidationError("call prepare() before inspecting or running the PDE instance")
        return self._circuit

    def _result(self, state: RegisterState, backend: str) -> PDEResult:
        return PDEResult(state, InitialLayout(self.problem.grid, self.problem.equation.fields), backend, self.prepared.amplitude_scale)

    def run_pysparq(
        self, *, max_steps: int = 1_000_000, max_states: int = 65536,
        native_registry: NativeRegistry | None = None,
    ) -> PDEResult:
        """Execute with the real optional PySparQ backend and registered memory."""
        circuit = self.circuit()
        state = circuit.run_pysparq( max_steps=max_steps, max_states=max_states, native_registry=native_registry)
        return self._result(state, "pysparq")

    def run_originir_ext(self, *, max_qubits: int = 24, max_steps: int = 1_000_000) -> PDEResult:
        """Export and execute through the real optional UnifiedQuantum backend."""
        circuit = self.circuit()
        return self._result(circuit.run_originir_ext(max_qubits=max_qubits, max_steps=max_steps), "originir_ext")


def qpde_solve(
    problem: PDEProblem, *, method: str | PDEMethod = "QHAM", config: QHAMConfig | None = None,
) -> PDESolveInstance:
    """Create an unprepared solve instance with a named or application-defined method.

    Only QHAM is a built-in named method. Other methods implement ``PDEMethod``
    and own their input contracts; the host lifecycle and backends are shared.
    """
    require_instance(problem, PDEProblem, "PDE.problem")
    if isinstance(method, str):
        if method.upper() != "QHAM":
            raise ValidationError("unknown PDE method; use QHAM or supply a PDEMethod implementation")
        selected: PDEMethod = QHAMMethod(config or QHAMConfig())
    else:
        if config is not None:
            raise ValidationError("configure a custom PDE method directly instead of passing QHAM config")
        if not callable(getattr(method, "prepare", None)):
            raise ValidationError("a PDE method must implement prepare(problem)")
        selected = method
    return PDESolveInstance(problem, selected)
