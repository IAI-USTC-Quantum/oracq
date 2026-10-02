"""Prepared ODE solve instances with replaceable methods and physical vector readout."""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from functools import partial
from typing import Protocol

from oracq.algorithms.common.hamiltonian import taylor_hamiltonian
from oracq.algorithms.common.solve import SolveCircuit, resolve_memory
from oracq.algorithms.input_model.contracts import (
    ContractReport,
    finite_real,
    positive_integer,
    require_instance,
)
from oracq.algorithms.input_model.initial import PreparedInitial
from oracq.algorithms.input_model.interfaces import operator_state_contract
from oracq.algorithms.input_model.ode import (
    ODELayout,
    ODEProblem,
    PreparedGenerator,
    physical_generator,
)
from oracq.algorithms.input_model.operators import BlockEncoding
from oracq.algorithms.input_model.oracles import StateOracle, annotate
from oracq.algorithms.input_model.qham import taylor_qode
from oracq.algorithms.qode.cbmd import ContourPlan
from oracq.algorithms.qode.lchs import QuadraturePlan
from oracq.algorithms.qode.ode import QODEProblem, QODESolver, linear_qode
from oracq.algorithms.qode.schrodingerization import SchrodingerPlan
from oracq.infrastructure.builder import Builder
from oracq.infrastructure.execution import RegisterState
from oracq.infrastructure.ir import Bits, ValidationError
from oracq.infrastructure.native import NativeRegistry
from oracq.infrastructure.qram_schema import RegisteredQRAM
from oracq.infrastructure.validation import validate


@dataclass(frozen=True)
class ODEConfig:
    """Linear evolution settings; approximation choices remain method-specific."""

    taylor_degree: int = 2
    plan: QuadraturePlan | ContourPlan | SchrodingerPlan | None = None
    hamiltonian_function: Callable[[BlockEncoding, float], BlockEncoding] = taylor_hamiltonian
    max_initial_words: int = 4096

    def __post_init__(self) -> None:
        positive_integer(self.taylor_degree, "ODE.taylor_degree", minimum=0)
        positive_integer(self.max_initial_words, "ODE.max_initial_words")
        if not callable(self.hamiltonian_function):
            raise ValidationError("ODE hamiltonian_function must be callable")


@dataclass(frozen=True)
class PreparedODE:
    """Method output, checked contracts, memory, and optional physical magnitude scale."""

    state: StateOracle
    report: ContractReport
    qrams: tuple[RegisteredQRAM, ...] = ()
    model: object | None = None
    amplitude_scale: float | None = None

    def __post_init__(self) -> None:
        if self.amplitude_scale is not None:
            finite_real(self.amplitude_scale, "ODE.amplitude_scale", minimum=0)


class ODEMethod(Protocol):
    """A method owns its input adaptations, mathematical requirements, and readout scale."""

    def prepare(self, problem: ODEProblem) -> PreparedODE:
        """Adapt the problem to this method and return the checked executable output."""


@dataclass(frozen=True)
class LinearODEMethod:
    """Adapt shared ODE inputs to an existing checked QODESolver protocol."""

    solver: QODESolver
    max_initial_words: int = 4096
    recover_taylor_scale: bool = False

    def __post_init__(self) -> None:
        require_instance(self.solver, QODESolver, "ODE.solver")
        positive_integer(self.max_initial_words, "ODE.max_initial_words")

    def prepare(self, problem: ODEProblem) -> PreparedODE:
        """Resolve inputs, enforce method promises, and generate one physical state oracle."""
        layout = ODELayout(problem.size, problem.components, self.max_initial_words)
        initial = require_instance(problem.initial.prepare_initial(layout), PreparedInitial, "ODE.initial")
        generator = require_instance(problem.generator.prepare_generator(layout), PreparedGenerator, "ODE.generator")
        # Validate the raw provider before composing padding projectors.
        operator_state_contract(self.solver.name).check(generator=generator.encoding, initial=initial.preparation).require()
        model = QODEProblem(
            physical_generator(generator.encoding, layout), initial.preparation,
            dissipative=problem.dissipative, initial_norm=initial.norm or None,
            evidence=problem.evidence,
        )
        report = self.solver.check(model, time=problem.final_time).require()
        if initial.norm == 0:
            b = Builder("zero_ode_solution", {"target": Bits(layout.width), "signal": Bits(1)})
            b.x(b["signal"])
            return PreparedODE(StateOracle(annotate(b.finish(), "unitary", algorithm="zero_ode_solution")), report, model=model, amplitude_scale=1.0)
        state = self.solver.solve(model, problem.final_time)
        scale = None
        if self.recover_taylor_scale:
            alpha = dict(state.operation.module.attributes).get("evolution_alpha")
            if not isinstance(alpha, (int, float)):
                raise ValidationError("the Taylor solver did not supply its evolution normalization")
            scale = initial.norm * alpha
            if not math.isfinite(scale):
                scale = None
        return PreparedODE(state, report, (*generator.qrams, *initial.qrams), model, scale)


@dataclass(frozen=True)
class ODECircuit(SolveCircuit):
    """A modular ODE circuit sharing export, resource estimates, and backend execution."""


@dataclass(frozen=True)
class ODEResult:
    """Backend state and final-time vector readout, excluding register padding."""

    state: RegisterState
    layout: ODELayout
    backend: str
    amplitude_scale: float | None = None

    def _physical(self) -> tuple[complex, ...]:
        names = [r.name for r in self.state.registers]
        target, signal = names.index("target"), names.index("signal")
        amplitudes = {values[target]: amplitude for values, amplitude in self.state.amplitudes.items() if values[signal] == 0}
        return tuple(amplitudes.get(i, 0j) for i in range(self.layout.size))

    @property
    def success_probability(self) -> float:
        """Joint probability of zero signal and a physical vector address."""
        return sum(abs(v) ** 2 for v in self._physical())

    def physical_amplitudes(self, *, normalize: bool = True) -> tuple[complex, ...]:
        """Return physical vector amplitudes, optionally conditioned on success."""
        values = self._physical()
        if not normalize:
            return values
        probability = self.success_probability
        if probability <= 0:
            raise ValidationError("the physical success channel has zero probability")
        return tuple(v / math.sqrt(probability) for v in values)

    def physical_values(self) -> tuple[complex, ...]:
        """Recover the generated approximation's magnitude when the method supplies its scale."""
        if self.amplitude_scale is None:
            raise ValidationError("this ODE method did not supply physical magnitude recovery; read conditional amplitudes instead")
        return tuple(v * self.amplitude_scale for v in self._physical())


class ODESolveInstance:
    """Cached prepare/inspect/execute lifecycle, independent of native quantum backends."""

    def __init__(self, problem: ODEProblem, method: ODEMethod) -> None:
        self._problem, self._method = problem, method
        self._prepared: PreparedODE | None = None
        self._circuit: ODECircuit | None = None

    @property
    def problem(self) -> ODEProblem:
        """Captured mathematical problem and input choices."""
        return self._problem

    @property
    def method(self) -> ODEMethod:
        """Selected method and its configuration."""
        return self._method

    @property
    def prepared(self) -> PreparedODE:
        """Inspect the method's checked model and readout scale after preparation."""
        if self._prepared is None:
            raise ValidationError("call prepare() before inspecting or running the ODE instance")
        return self._prepared

    def prepare(self) -> ODESolveInstance:
        """Generate once, check closure and output layout, and resolve QRAM snapshots."""
        if self._prepared is None:
            prepared = require_instance(self.method.prepare(self.problem), PreparedODE, "ODE.method")
            prepared.report.require()
            layout = ODELayout(self.problem.size, self.problem.components)
            if prepared.state.width != layout.width:
                raise ValidationError("an ODE method must return the physical target register width")
            program = validate(prepared.state.operation.program(), require_closed=True)
            memory = resolve_memory(program, prepared.qrams)
            circuit = ODECircuit(program, tuple((name, tuple(sorted(cells.items()))) for name, cells in sorted(memory.items())))
            self._prepared, self._circuit = prepared, circuit
        return self

    def circuit(self) -> ODECircuit:
        """Inspect the exact prepared circuit used for export and execution."""
        if self._circuit is None:
            raise ValidationError("call prepare() before inspecting or running the ODE instance")
        return self._circuit

    def _result(self, state: RegisterState, backend: str) -> ODEResult:
        return ODEResult(state, ODELayout(self.problem.size, self.problem.components), backend, self.prepared.amplitude_scale)

    def run_pysparq(
        self, *, max_steps: int = 1_000_000, max_states: int = 65536,
        native_registry: NativeRegistry | None = None,
    ) -> ODEResult:
        """Execute the prepared circuit through the real optional PySparQ backend."""
        return self._result(self.circuit().run_pysparq(max_steps=max_steps, max_states=max_states, native_registry=native_registry), "pysparq")

    def run_originir_ext(self, *, max_qubits: int = 24, max_steps: int = 1_000_000) -> ODEResult:
        """Execute the prepared circuit through the real optional UnifiedQuantum backend."""
        return self._result(self.circuit().run_originir_ext(max_qubits=max_qubits, max_steps=max_steps), "originir_ext")


def qode_solve(
    problem: ODEProblem, *, method: str | QODESolver | ODEMethod = "taylor", config: ODEConfig | None = None,
) -> ODESolveInstance:
    """Select Taylor, LCHS, CBMD, Schrödingerization, or an application-defined method.

    Named methods use ODEConfig. A QODESolver keeps its own configuration and
    preconditions; an ODEMethod owns preparation and physical scale recovery.
    """
    require_instance(problem, ODEProblem, "ODE.problem")
    selected: ODEMethod
    if isinstance(method, str):
        cfg = require_instance(config, ODEConfig, "ODE.config") if config is not None else ODEConfig()
        name = method.lower()
        if name == "taylor":
            if cfg.plan is not None:
                raise ValidationError("Taylor evolution does not accept a quadrature or auxiliary-grid plan")
            solver = QODESolver("taylor", partial(taylor_qode, degree=cfg.taylor_degree))
        else:
            solver = linear_qode(name, plan=cfg.plan, hamiltonian_function=cfg.hamiltonian_function)
        selected = LinearODEMethod(solver, cfg.max_initial_words, recover_taylor_scale=name == "taylor")
    else:
        if config is not None:
            raise ValidationError("configure a custom ODE method or QODESolver directly instead of passing ODEConfig")
        if isinstance(method, QODESolver):
            selected = LinearODEMethod(method)
        elif callable(getattr(method, "prepare", None)):
            selected = method
        else:
            raise ValidationError("an ODE method must implement prepare(problem) or be a QODESolver")
    return ODESolveInstance(problem, selected)
