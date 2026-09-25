"""CKS section 5 VTAA linear-system solver: GPE clock, banded inverse LCU, and variable-time amplitude amplification.

Based on Childs–Kothari–Somma (arXiv:1511.02306 §5, SIAM J. Comput. 2017). Variable-time layer structure:
the clock bits C_j are written by gapped phase estimation (Lemma 22), which controlled-triggers the banded
Chebyshev inverse LCU (Lemma 23, including the alpha_max uniformization rotation of eq. (98)); the outer layer
applies the clock-aware nested amplification of Ambainis (arXiv:1010.4458) (operator form taken from
Low–Su arXiv:2410.18178 eq. (47)–(53)); finally the GPE garbage is erased with the inverse of A'
(eq. (99)–(112)).

The GPE decision follows the deterministic route of Low–Su Prop 23: a Yoder–Low–Chuang fixed-point decision
polynomial (reusing the verified synthesis of qsvt.fixed_point_search_phases) is applied to the qubitization
walk, and the decision bit is flipped on the signal==0 branch: when the magnitude of lambda is at least
theta_j the decision is 1 (fire, hard bound epsilon) with amplitude near 1, while at x=0 the decision is 0;
the transition band response is deterministic and exactly computable. The query count
O((alpha/theta_j) log(1/eps)) is of the same order as the PEA plus majority vote of CKS Lemma 22; the P_j
garbage is kept as in the paper and erased by the inverse of A'. No public implementation of VTAA or
CKS §5 was available before this one.
"""

from __future__ import annotations

import cmath
import json
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import lru_cache

from oracq.algorithms.common.arithmetic import BooleanNetwork
from oracq.algorithms.common.qsvt import fixed_point_search_phases, qsp_response
from oracq.algorithms.common.transforms import qsvt_sequence
from oracq.algorithms.input_model.block_encoding import reflect_zero
from oracq.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
)
from oracq.algorithms.input_model.interfaces import as_block_encoding
from oracq.algorithms.input_model.operators import BlockEncoding, _name
from oracq.algorithms.input_model.oracles import (
    StateOracle,
    annotate,
    gate_state_prep,
    invoke,
    resources_for,
)
from oracq.algorithms.input_model.sparse import chebyshev_block, real_symmetric_sparse_encoding
from oracq.algorithms.qlss.qlss import CKSConfig, QLSSProtocol, SparseSystem
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError, fuse


@lru_cache(maxsize=64)
def gpe_fire_phases(
    threshold: float, x_edge: float, epsilon: float, degree_cap: int = 40
) -> tuple[tuple[float, ...], int]:
    """GPE fire decision phases: YLC fixed-point polynomial with ``|P| >= 1-epsilon`` when ``|x| >= threshold``.

    The threshold decreases monotonically with degree; odd degrees from 3 upward are searched until
    the [threshold, x_edge] grid validation passes; at x=0, P=0, and the response over the transition
    band (0, threshold) is deterministic and exactly computable (a band CKS makes no promise about).

    Args:
        threshold: Fire decision threshold (encoded spectral units), in (0, x_edge].
        x_edge: Right end of the validation grid, in [threshold, 1].
        epsilon: Hard bound on the decision response, in (0,1).
        degree_cap: Search cap for the decision polynomial degree; odd degrees from 3 upward are searched.

    Returns:
        tuple[tuple[float, ...], int]: The grid-validated QSP phase sequence and the actual degree.
    """
    finite_real(threshold, "gpe_fire.threshold", minimum=0, strict=True)
    finite_real(x_edge, "gpe_fire.x_edge", minimum=0, strict=True)
    finite_real(epsilon, "gpe_fire.epsilon", minimum=0, strict=True)
    if not 0 < threshold <= x_edge <= 1 or not 0 < epsilon < 1:
        raise ValidationError("Invalid GPE decision geometry or precision")
    fire_grid = [threshold + (x_edge - threshold) * i / 16.0 for i in range(17)]
    for degree in range(3, degree_cap + 1, 2):
        try:
            phases = fixed_point_search_phases(epsilon, degree)
        except ValidationError:
            continue
        if min(abs(qsp_response(x, phases)) for x in fire_grid) >= 1 - epsilon:
            return phases, degree
    raise ValidationError("The GPE decision polynomial degree exceeded the cap; increase phi or lower the kappa declaration")


@lru_cache(maxsize=64)
def clock_or_operation(width: int) -> Operation:
    """Clock-prefix OR predicate: the coherent criterion for stopped<=j, used by the VTAA reflections.

    Args:
        width: Number of bits in the clock prefix; a positive integer.

    Returns:
        Operation: A Boolean-circuit operation writing the OR of the prefix bits into the 1-bit stopped register.
    """
    net = BooleanNetwork()
    bits = net.input("prefix", width)
    net.outputs = {"stopped": [net.any(bits)]}
    return net.operation(attributes={"algorithm": "clock_prefix_or", "width": width})


def gapped_phase_estimation(
    a: BlockEncoding,
    threshold: float,
    x_edge: float,
    *,
    epsilon: float = 0.02,
    degree_cap: int = 40,
) -> Operation:
    """GPE of CKS Lemma 22 (deterministic QSP route, Low–Su arXiv:2410.18178 Prop 23).

    Applies the YLC fixed-point decision polynomial P(lambda/alpha) to the qubitization walk of
    the BE and flips the decision bit on the signal==0 branch: when ``|lambda/alpha| >= threshold``
    (fire hard bound epsilon) the decision is 1 with amplitude near 1 (stop in this band and
    invert), while at the low end of the spectrum the decision tends to 0 (deferred to a finer
    band). The signal register keeps gamma garbage (P_j of the paper §5.2), erased at the end by
    the inverse of A'.

    Args:
        a: Hermitian BE of the self-adjoint unitary dilation.
        threshold: Fire decision threshold (encoded spectral units, i.e. relative to alpha).
        x_edge: Spectral upper bound divided by alpha; right end of the validation grid.
        epsilon: Hard bound on the fire-side decision response.
        degree_cap: Cap on the decision polynomial degree.

    Returns:
        Operation: Registers target/signal/decision; decision=1 means stopping in this band.
    """
    a = as_block_encoding(a)
    phases, degree = gpe_fire_phases(threshold, x_edge, epsilon, degree_cap)
    marker = qsvt_sequence(a, phases)
    b = Builder(
        _name("vtaa_gpe", a.operation, threshold, x_edge, epsilon, degree),
        {
            "target": Bits(a.width),
            "signal": Bits(a.signal_qubits),
            "decision": Bits(1),
        },
        resources_for(("marker", marker)),
        attributes={
            "algorithm": "gapped_phase_estimation",
            "fire_threshold": threshold,
            "marker_degree": degree,
            "marker_epsilon": epsilon,
            "walk_queries_max": degree,
            "references": "CKS Lemma 22 via Low-Su Prop 23 deterministic QSP decision",
            "validation_stage": "paradigm",
        },
    )
    invoke(b, marker, "marker", target=b["target"], signal=b["signal"])
    with b.control(b["signal"], 0):
        b.x(b["decision"])
    return b.finish()


def band_inverse_step(
    a: BlockEncoding, coefficients: Sequence[float], alpha_max: float
) -> Operation:
    """W(lambda, delta) of CKS Lemma 23: banded Chebyshev inverse LCU plus uniformization rotation.

    Args:
        a: BE of the self-adjoint unitary dilation; chebyshev_block provides the odd Chebyshev powers.
        coefficients: Truncated inverse polynomial coefficients c_k of this band (T_{2k+1} basis).
        alpha_max: Maximum LCU 1-norm over all bands; eq. (98) uniformly scales the success amplitude down to 1/alpha_max.

    Returns:
        Operation: Registers target/signal/flag; the branch with signal==0 and flag==1 carries
        h(A)psi/alpha_max (psi being the input state), where h is this band's polynomial.
    """
    a = as_block_encoding(a)
    terms = tuple(
        (complex(c), chebyshev_block(a, 2 * k + 1))
        for k, c in enumerate(coefficients)
        if c != 0
    )
    if not terms:
        raise ValidationError("All band inverse polynomial coefficients are zero")
    beta = sum(abs(c) for c, _ in terms)
    if alpha_max <= 0 or beta > alpha_max * (1 + 1e-9):
        raise ValidationError("The band LCU normalization conflicts with the alpha_max declaration")
    selector_width = (len(terms) - 1).bit_length()
    b = Builder(
        _name("vtaa_band_inverse", a.operation, tuple(coefficients), alpha_max),
        {
            "target": Bits(a.width),
            "signal": Bits(selector_width + a.signal_qubits),
            "flag": Bits(1),
        },
        resources_for(*[(f"term{i}", t.operation) for i, (_, t) in enumerate(terms)]),
        attributes={
            "algorithm": "vtaa_band_inverse",
            "band_lcu_normalization": beta,
            "uniform_normalization": alpha_max,
            "polynomial_terms": len(terms),
            "implementation_scope": "CKS eq. 96-98 piecewise inverse with fixed alpha_max",
            "validation_stage": "paradigm",
        },
    )
    selector, work = b["signal"][:selector_width], b["signal"][selector_width:]
    if selector_width:
        weights = [math.sqrt(abs(c)) for c, _ in terms]
        weights += [0.0] * ((1 << selector_width) - len(weights))
        prep = gate_state_prep(weights)
        invoke(b, prep.operation, target=selector, work=selector[:0])
        for i, (coefficient, term) in enumerate(terms):
            with b.control(selector, i):
                b.global_phase(cmath.phase(coefficient))
                invoke(b, term.operation, f"term{i}", target=b["target"], signal=work)
        with b.adjoint():
            invoke(b, prep.operation, target=selector, work=selector[:0])
    else:
        coefficient, term = terms[0]
        b.global_phase(cmath.phase(coefficient))
        invoke(b, term.operation, "term0", target=b["target"], signal=work)
    # Eq. (98): on the LCU success branch, rotate the flag to amplitude beta/alpha_max, uniformizing the success amplitudes across bands.
    with b.control(b["signal"], 0):
        b.ry(b["flag"], 2 * math.asin(beta / alpha_max))
    return b.finish()


@dataclass(frozen=True)
class VTAAConfig:
    """order is the base Chebyshev degree of band 1; band j scales it geometrically by 2^(j-1).

    rounds holds the VTAA amplification rounds r_j per stage (None means all zero, i.e. the pure
    variable-time layer plus postselection); correctness does not depend on the schedule, which
    only affects the success rate. Production deployments should choose r_j via the amplitude
    estimation of Ambainis Algorithm 2 or the deterministic schedule of Low–Su (see tunable_rounds).
    """

    order: int = 2
    terms: int | None = None
    clock_steps: int | None = None
    marker_epsilon: float = 0.02
    degree_cap: int = 40
    rounds: tuple[int, ...] | None = None

    def __post_init__(self) -> None:
        """Validate the degree, truncation, clock steps, decision hard bound, and amplification rounds."""
        positive_integer(self.order, "VTAAConfig.order", maximum=128)
        if self.terms is not None:
            positive_integer(self.terms, "VTAAConfig.terms", maximum=self.order)
        if self.clock_steps is not None:
            positive_integer(self.clock_steps, "VTAAConfig.clock_steps", maximum=16)
        finite_real(self.marker_epsilon, "VTAAConfig.marker_epsilon", minimum=0, strict=True)
        if self.marker_epsilon >= 0.2:
            raise ValidationError("The GPE decision hard bound must be below 0.2")
        positive_integer(self.degree_cap, "VTAAConfig.degree_cap", minimum=8, maximum=40)
        if self.rounds is not None:
            rounds = tuple(int(r) for r in self.rounds)
            if not rounds or any(not 0 <= r <= 64 for r in rounds):
                raise ValidationError("VTAA rounds must be non-negative and not exceed 64")
            object.__setattr__(self, "rounds", rounds)

    def band_order(self, step: int) -> int:
        """Return the Chebyshev inverse polynomial degree of band ``step``.

        Band 1 takes ``order`` and each subsequent band doubles it (geometric growth in
        powers of 2), capped at 128; ``step`` below 1 is treated as band 1.

        Args:
            step: Band index, counting from 1.

        Returns:
            int: The Chebyshev degree of that band, at most 128.
        """
        return min(128, self.order * (1 << max(0, step - 1)))

    def band_coefficients(self, step: int) -> tuple[float, ...]:
        """Return the truncated inverse polynomial coefficients of band ``step`` (T_{2k+1} basis).

        Computed via ``CKSConfig.coefficients`` with degree ``band_order(step)`` and truncation
        term count ``terms`` (None means the full degree).

        Args:
            step: Band index, counting from 1.

        Returns:
            tuple[float, ...]: The truncated inverse polynomial coefficient sequence on the odd Chebyshev basis.
        """
        return CKSConfig(self.band_order(step), self.terms).coefficients()


def tunable_rounds(
    stage_amplitudes: Iterable[float], thresholds: Iterable[float] | None = None
) -> tuple[int, ...]:
    """Low–Su tunable VTAA schedule (arXiv:2410.18178 eq. (52)–(53)).

    Args:
        stage_amplitudes: Estimates of the "not yet failed" amplitude norms x_j per stage (obtainable via the amplitude estimation channel).
        thresholds: Thresholds alpha_j, uniformly split by default; their sum is constant to guarantee constant loss.

    Returns:
        A tuple of per-stage rounds r_j = max(ceil(sqrt(alpha_j)/(6 x_j) - 1/2), 0),
        satisfying the no-overshoot condition (2 r_j + 1) x_j <= sqrt(alpha_j).
    """
    values = tuple(float(v) for v in stage_amplitudes)
    if not values or any(not math.isfinite(v) or not 0 < v <= 1 for v in values):
        raise ValidationError("Stage norms must be finite values between 0 exclusive and 1 inclusive")
    thresholds = (
        tuple(1.0 / len(values) for _ in values)
        if thresholds is None
        else tuple(float(t) for t in thresholds)
    )
    if len(thresholds) != len(values) or any(
        not math.isfinite(t) or t <= 0 for t in thresholds
    ):
        raise ValidationError("Invalid VTAA threshold count or values")
    return tuple(
        max(0, math.ceil(math.sqrt(alpha) / (6 * x) - 0.5))
        for x, alpha in zip(values, thresholds, strict=True)
    )


def vtaa_cks(system: SparseSystem, config: VTAAConfig | None = None) -> StateOracle:
    """CKS §5 VTAA solve kernel; takes a SparseSystem and returns the solution StateOracle.

    Args:
        system: Sparse Hermitian linear system; must carry Hermitian and spectral bound declarations.
        config: VTAA configuration; defaults to the default configuration.

    Returns:
        StateOracle: The solution state oracle output by the variable-time amplification cascade, with correctness marked pending.
    """
    config = config or VTAAConfig()
    require_instance(config, VTAAConfig, "vtaa_cks.config")
    if not system.hermitian:
        raise ValidationError("VTAA sparse input requires a Hermitian declaration or an explicit Hermitian dilation")
    a = real_symmetric_sparse_encoding(
        system.access,
        system.value_format,
        system.entry_bound,
        diagonal_nonnegative=system.diagonal_nonnegative,
    )
    x_edge = system.spectrum.norm_upper / a.alpha
    physical_kappa = system.spectrum.norm_upper / system.spectrum.sigma_min_lower
    # Coverage condition: the fire threshold of the finest band, 2^(1-steps)*x_edge, must not exceed sigma_min/alpha.
    if config.clock_steps is not None:
        steps = config.clock_steps
    else:
        steps = math.ceil(math.log2(physical_kappa)) + 1 if physical_kappa > 1 else 1
    if physical_kappa > 2 ** (steps - 1):
        raise ValidationError("clock_steps cannot cover the declared minimum singular value band")
    if config.rounds is not None and len(config.rounds) != steps:
        raise ValidationError("VTAA rounds length must equal clock_steps")

    bands: list[tuple[float, tuple[float, ...], float]] = []
    for step in range(1, steps + 1):
        threshold = x_edge * 2.0 ** (1 - step)
        coefficients = config.band_coefficients(step)
        bands.append((threshold, coefficients, sum(abs(c) for c in coefficients)))
    alpha_max = max(beta for _, _, beta in bands)
    gpes = tuple(
        gapped_phase_estimation(
            a,
            threshold,
            x_edge,
            epsilon=config.marker_epsilon,
            degree_cap=config.degree_cap,
        )
        for threshold, _, _ in bands
    )
    inverses = tuple(
        band_inverse_step(a, coefficients, alpha_max) for _, coefficients, _ in bands
    )
    rounds = config.rounds or (0,) * steps
    selector_width = max(
        (len([c for c in coefficients if c != 0]) - 1).bit_length()
        for _, coefficients, _ in bands
    )
    # Signal layout: [shared inverse-LCU region selw+s][per-step independent GPE garbage m*s]; P_j garbage is kept for later erasure.
    signal_width = selector_width + a.signal_qubits * (steps + 1)
    registers = {
        "target": Bits(a.width),
        "clock": Bits(steps),
        "flag": Bits(1),
        "signal": Bits(signal_width),
    }
    gpe_base = selector_width + a.signal_qubits

    def build_step(step: int, *, uncompute: bool) -> Operation:
        """A_j: controlled GPE writes C_j; when C_j=1 (fire) W_j is applied (A'_j only flips the flag, eq. (99))."""
        threshold, coefficients, beta = bands[step - 1]
        kind = "vtaa_uncompute_step" if uncompute else "vtaa_variable_step"
        resources = resources_for(("gpe", gpes[step - 1]))
        if not uncompute:
            resources.update(resources_for(("inverse", inverses[step - 1])))
        b = Builder(
            _name(kind, a.operation, system.rhs.operation, step, coefficients, config),
            registers,
            resources,
            attributes={
                "algorithm": kind,
                "band_index": step,
                "fire_threshold": threshold,
                "band_lcu_normalization": beta,
                "validation_stage": "paradigm",
            },
        )
        prefix = b["clock"][: step - 1]
        gpe_slice = b["signal"][gpe_base + (step - 1) * a.signal_qubits :][
            : a.signal_qubits
        ]
        inverse_width = (
            (len([c for c in coefficients if c != 0]) - 1).bit_length() + a.signal_qubits
        )

        def run_gpe() -> None:
            """Invoke the GPE on this step's dedicated signal slot, writing the decision into the clock bit."""
            invoke(
                b,
                gpes[step - 1],
                "gpe",
                target=b["target"],
                signal=gpe_slice,
                decision=b["clock"][step - 1],
            )

        # All-zero prefix = not yet stopped; C_j=1 means fire in this band.
        if prefix.width:
            with b.control(prefix, 0):
                run_gpe()
        else:
            run_gpe()
        with b.control(b["clock"][step - 1]):
            if uncompute:
                b.x(b["flag"])
            else:
                invoke(
                    b,
                    inverses[step - 1],
                    "inverse",
                    target=b["target"],
                    signal=b["signal"][:inverse_width],
                    flag=b["flag"],
                )
        return b.finish()

    initial = Builder(
        _name("vtaa_initial", system.rhs.operation),
        registers,
        resources_for(("rhs", system.rhs.operation)),
    )
    rhs_work = initial.local("rhs_work", Bits(system.rhs.work_width))
    invoke(initial, system.rhs.operation, "rhs", target=initial["target"], work=rhs_work)
    initial_op = initial.finish()

    def run_chain(b: Builder, items: Sequence[tuple[str, Operation]]) -> None:
        """Wire the operations in ``items`` into ``b`` one by one under their own prefixes."""
        for prefix, op in items:
            invoke(
                b,
                op,
                prefix,
                target=b["target"],
                clock=b["clock"],
                flag=b["flag"],
                signal=b["signal"],
            )

    def prefix_module(step: int, items: Sequence[tuple[str, Operation]]) -> Operation:
        """Wrap the operation chain up to step ``step`` into a single module."""
        b = Builder(
            _name("vtaa_prefix", a.operation, system.rhs.operation, config, step),
            registers,
            resources_for(*items),
            attributes={
                "algorithm": "vtaa_prefix",
                "stage": step,
                "validation_stage": "paradigm",
            },
        )
        run_chain(b, items)
        return b.finish()

    def amplification(prefix_op: Operation, step: int, count: int) -> Operation:
        """M_j = (R_s R_f)^{r_j} P_j; R_f flips the phase of stopped<=j with flag 0."""
        b = Builder(
            _name("vtaa_amplified", prefix_op, step, count),
            registers,
            resources_for(("prefix", prefix_op)),
            attributes={
                "algorithm": "vtaa_amplified_stage",
                "stage": step,
                "rounds": count,
                "implementation_scope": "Ambainis VTAA cascade; operator form Low-Su eq. 47",
                "validation_stage": "paradigm",
            },
        )
        run_chain(b, (("prefix", prefix_op),))

        def run_prefix() -> None:
            """Replay the amplified prefix chain once."""
            run_chain(b, (("prefix", prefix_op),))

        if count:
            orop = clock_or_operation(step)
            stopped = b.local("stopped", Bits(1))

            def mark_stopped() -> None:
                """Write the OR decision of the clock prefix into the ``stopped`` ancilla bit."""
                b.call(orop, prefix=b["clock"][:step], stopped=stopped)

            with b.repeat(count):
                mark_stopped()
                # Branches with stopped=1 (already fired) and flag=0 (failed) acquire a pi phase.
                with b.control(fuse(stopped, b["flag"]), 1):
                    b.global_phase(math.pi)
                mark_stopped()
                with b.adjoint():
                    run_prefix()
                reflect_zero(
                    b, fuse(b["target"], b["clock"], b["flag"], b["signal"]), positive=True
                )
                run_prefix()
        return b.finish()

    chain = [("start", initial_op)]
    amplified = initial_op
    for step in range(1, steps + 1):
        chain.append((f"step{step}", build_step(step, uncompute=False)))
        amplified = amplification(prefix_module(step, chain), step, rounds[step - 1])
        chain = [("run", amplified)]
    erase_ops = tuple(build_step(step, uncompute=True) for step in range(1, steps + 1))

    top = Builder(
        _name("vtaa_cks", a.operation, system.rhs.operation, config),
        {"target": Bits(a.width), "signal": Bits(steps + 1 + signal_width)},
        resources_for(
            ("run", amplified), *((f"erase{s}", op) for s, op in enumerate(erase_ops, 1))
        ),
        attributes={
            "algorithm": "vtaa_cks",
            "input_model": "sparse_location_inplace_and_entry_xor",
            "input_alpha": a.alpha,
            "encoded_inverse_bound": system.spectrum.inverse_bound(a.alpha),
            "clock_steps": steps,
            "fire_thresholds": json.dumps([t for t, _, _ in bands]),
            "band_orders": json.dumps([config.band_order(s) for s in range(1, steps + 1)]),
            "band_lcu_normalizations": json.dumps([beta for _, _, beta in bands]),
            "alpha_max": alpha_max,
            "marker_degrees": json.dumps(
                [dict(op.module.attributes)["marker_degree"] for op in gpes]
            ),
            "marker_epsilon": config.marker_epsilon,
            "rounds": json.dumps(list(rounds)),
            "implementation_scope": (
                "CKS arXiv:1511.02306 section 5 VTAA; GPE per Lemma 22 realized as the"
                " deterministic QSP decision of Low-Su arXiv:2410.18178 Prop 23 on the"
                " qubitization walk with YLC fixed-point markers; clock-aware cascade"
                " per Ambainis arXiv:1010.4458 with operator form from Low-Su eq. 47-53"
            ),
            "kernel_status": "prototype; band polynomial accuracy and VTAA schedule pending",
            "correctness": "pending",
            "success_condition": "signal == 0; erased branch carries h(A)|b> mixture",
            "validation_stage": "paradigm",
        },
    )
    clock, flag, rest = (
        top["signal"][:steps],
        top["signal"][steps],
        top["signal"][steps + 1 :],
    )

    def call_top(op: Operation, prefix: str) -> None:
        """Wire ``op`` into the top-level circuit's clock, flag, and signal layout under the given prefix."""
        invoke(top, op, prefix, target=top["target"], clock=clock, flag=flag, signal=rest)

    call_top(amplified, "run")
    # (A')^dagger: erase the GPE clock and flag in reverse order (eq. (110)–(112)); the solution state stays in target.
    with top.adjoint():
        for step in range(1, steps + 1):
            call_top(erase_ops[step - 1], f"erase{step}")
    return StateOracle(
        annotate(
            top.finish(),
            "unitary",
            algorithm="vtaa_cks",
            implementation_scope="CKS section 5 variable-time amplitude amplification",
            correctness="pending",
        )
    )


def make_vtaa_cks_qlss(config: VTAAConfig | None = None) -> QLSSProtocol:
    """CKS §5 VTAA solve protocol; the input model is the same sparse access as cks_chebyshev.

    Args:
        config: VTAA configuration; defaults to the default configuration.

    Returns:
        QLSSProtocol: The VTAA solve protocol under the sparse input model.
    """
    config = VTAAConfig() if config is None else config
    require_instance(config, VTAAConfig, "make_vtaa_cks_qlss.config")
    return QLSSProtocol("vtaa_cks", "sparse", lambda system: vtaa_cks(system, config))
