"""Publication-grade numerical validation of the oracle catalog and query algorithms.

Validation targets: src/oracq/algorithms/oracles.py and oracle_algorithms.py.

- XOR database: the gate truth table exhausts the whole input domain in both
  per-basis-state and full-superposition modes; QRAM binding and gate/QRAM
  dual-binding consistency of abstract declarations; four backend paths
  cross-checked;
- Bernstein-Vazirani (an algorithm named in the paper): 3-8 bit multiple
  scales x multiple secret strings x both biases, two primary paths
  (originir-ext and rir-pysparq, plus reference and adapter cross-checks),
  with the probability of recovering the secret string exactly 1; includes
  the open-declaration -> binding end-to-end route;
- Deutsch-Jozsa: constant/balanced function-family decisions (including the
  end-to-end route of a BooleanNetwork-compiled oracle, which also
  cross-checks the circuit semantics of boolean-networks against classical
  evaluation);
- Simon: phase-sensitive pointwise comparison of the sampling distribution +
  GF(2) elimination to recover the period (including the wide-width CNOT
  linear-oracle route).

Run: PYTHONPATH=src <python with pysparq+uniqc> tests/verification/verify_oracles.py
"""

from __future__ import annotations

import math

from harness import (
    Report,
    adapter_pysparq,
    amplitude_error,
    amplitudes_to_statevector,
    originir_ext,
    probabilities,
    reference,
    rir_pysparq,
    statevector_error,
    tvd,
)

from oracq import Binding, Bits, Builder, bind
from oracq.algorithms.basics.oracle_algorithms import (
    affine_boolean_oracle,
    bernstein_vazirani,
    deutsch_jozsa,
    simon_nullspace,
    simon_sample,
)
from oracq.algorithms.common.arithmetic import BooleanNetwork
from oracq.algorithms.input_model.oracles import (
    XorDatabase,
    abstract_database,
    annotate,
    gate_database,
    qram_database,
)

# ---------------------------------------------------------------------------
# Classical reference tools independent of the library implementation
# ---------------------------------------------------------------------------


def _dot(a, b):
    """GF(2) dot product."""
    return (a & b).bit_count() & 1


def _pseudo_table(address_width, data_width, seed):
    """Deterministic xorshift64* pseudo-random truth table (independently constructed, no library code reused)."""
    state = seed | 1
    mask64 = (1 << 64) - 1
    words = []
    for _ in range(1 << address_width):
        state ^= state >> 12
        state ^= (state << 25) & mask64
        state ^= state >> 27
        words.append((state * 0x2545F4914F6CDD1D & mask64) & ((1 << data_width) - 1))
    return words


def _gf2_independent(vectors, count):
    """Deterministic selection of count linearly independent vectors by highest-bit pivot elimination."""
    pivots = {}
    chosen = []
    for vector in vectors:
        w = vector
        while w:
            pivot = w.bit_length() - 1
            if pivot in pivots:
                w ^= pivots[pivot]
            else:
                pivots[pivot] = w
                chosen.append(vector)
                break
        if len(chosen) == count:
            break
    return chosen


def _db_program(db_operation, *, superpose=(), initial=None, name="drive"):
    """XOR database driver.

    The harness basis/superposition helpers do not thread QRAM resources, so
    they are mapped explicitly by name here; combinations of basis-state
    initial values and partial-register superposition are supported.
    """
    module = db_operation.module
    b = Builder(
        name,
        {r.name: r.type for r in module.registers},
        {r.name: r.type for r in module.resources},
    )
    for key, value in (initial or {}).items():
        for bit in range(b[key].width):
            if (value >> bit) & 1:
                b.x(b[key][bit])
    for key in superpose:
        b.h(b[key])
    b.call(
        db_operation,
        resources={r.name: r.name for r in module.resources},
        **{r.name: b[r.name] for r in module.registers},
    )
    return b.finish().program()


def _xor_expected(address_width, data_width, table):
    """Exact reference state of the XOR database under full superposition: |a,d> -> |a, d XOR table[a]> uniform superposition."""
    uniform = 1 / math.sqrt(1 << (address_width + data_width))
    return {
        (a, d ^ table[a]): uniform
        for a in range(1 << address_width)
        for d in range(1 << data_width)
    }


def _bv_expected(secret, bias):
    """Closed-form final state of the BV/DJ phase-kickback circuit: input concentrated on s, answer left in |->."""
    sign = 1.0 if bias == 0 else -1.0
    root = 1 / math.sqrt(2)
    return {(secret, 0): sign * root, (secret, 1): -sign * root}


def _input_marginal(state_probs, input_index=0):
    marginal = {}
    for key, p in state_probs.items():
        marginal[key[input_index]] = marginal.get(key[input_index], 0.0) + p
    return marginal


def _origin_input_probability(vector, input_width, target):
    """Sum the probability of input == target from the OriginIR full-amplitude state vector (input in the low bits)."""
    mask = (1 << input_width) - 1
    return sum(abs(v) ** 2 for i, v in enumerate(vector) if (i & mask) == target)


# ---------------------------------------------------------------------------
# Group A: XOR database truth-table semantics
# ---------------------------------------------------------------------------


def verify_xor_gate_basis(report):
    """Basis-state mode: pointwise exhaustive sweep over all (address, data) initial states; the output must be a single basis state."""
    # (2,3) on all paths; (3,2) on reference + rir-pysparq to control total runtime
    configs = [
        (2, 3, [5, 7, 0, 3], ("reference", "rir-pysparq", "originir-ext")),
        (3, 2, [1, 0, 3, 2, 0, 1, 2, 3], ("reference", "rir-pysparq")),
    ]
    for aw, dw, table, paths in configs:
        db = gate_database(aw, dw, table)
        failures = 0
        for a in range(1 << aw):
            for d in range(1 << dw):
                program = _db_program(
                    db.operation, initial={"address": a, "data": d}, name=f"basis_{a}_{d}"
                )
                expected_key = (a, d ^ table[a])
                states = {}
                if "reference" in paths:
                    states["reference"] = reference(program)
                if "rir-pysparq" in paths:
                    states["rir"] = rir_pysparq(program)
                for state in states.values():
                    if set(state) != {expected_key} or abs(state[expected_key] - 1) > 1e-12:
                        failures += 1
                if "originir-ext" in paths:
                    vector = originir_ext(program)
                    index = a | ((d ^ table[a]) << aw)
                    if any(
                        abs(v - (1.0 if i == index else 0.0)) > 1e-12
                        for i, v in enumerate(vector)
                    ):
                        failures += 1
        report.case(
            f"xor-gate-basis-{aw}x{dw}",
            paths=list(paths),
            parameters={"address_width": aw, "data_width": dw, "inputs": 1 << (aw + dw)},
            metrics={"failures": failures},
            criterion="all 2^(aw+dw) basis-state inputs output exactly |a, d XOR table[a]> (failures == 0)",
            passed=failures == 0,
        )


def verify_xor_gate_superposition(report):
    """Superposition mode: one run exhausts the whole input domain; four paths cross-checked amplitude by amplitude against the classical permutation."""
    for aw, dw, seed in ((2, 3, 11), (3, 4, 17), (4, 3, 23)):
        table = _pseudo_table(aw, dw, seed)
        db = gate_database(aw, dw, table)
        program = _db_program(db.operation, superpose=("address", "data"), name=f"sup_{aw}_{dw}")
        expected = _xor_expected(aw, dw, table)
        ref = reference(program)
        worst = amplitude_error(ref, expected)
        for runner in (rir_pysparq, adapter_pysparq):
            worst = max(worst, amplitude_error(runner(program), expected))
        vector = originir_ext(program)
        worst = max(
            worst,
            statevector_error(vector, amplitudes_to_statevector(expected, [aw, dw])),
        )
        report.case(
            f"xor-gate-superposition-{aw}x{dw}",
            paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
            parameters={"address_width": aw, "data_width": dw, "branches": len(expected)},
            metrics={"max_error": worst},
            criterion="superposition branches match the classical permutation amplitude by amplitude (max_error < 1e-9)",
            passed=worst < 1e-9,
        )


def verify_xor_qram(report):
    """QRAM binding: superposition sweep + basis-state spot checks with non-zero initial data."""
    aw, dw = 3, 4
    table = _pseudo_table(aw, dw, 29)
    db = qram_database(aw, dw)
    memory = {"table": table}
    program = _db_program(db.operation, superpose=("address", "data"), name="qram_sup")
    expected = _xor_expected(aw, dw, table)
    worst = amplitude_error(reference(program, memory), expected)
    for runner in (rir_pysparq, adapter_pysparq):
        worst = max(worst, amplitude_error(runner(program, memory), expected))
    vector = originir_ext(program, memory)
    worst = max(
        worst, statevector_error(vector, amplitudes_to_statevector(expected, [aw, dw]))
    )
    failures = 0
    for a, d in ((0, 9), (3, 15), (7, 1)):
        single = _db_program(db.operation, initial={"address": a, "data": d}, name="qram_b")
        state = rir_pysparq(single, memory)
        if set(state) != {(a, d ^ table[a])}:
            failures += 1
    report.case(
        "xor-qram-superposition-3x4",
        paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
        parameters={"address_width": aw, "data_width": dw, "branches": len(expected)},
        metrics={"max_error": worst, "basis_failures": failures},
        criterion="QRAM superposition matches amplitude by amplitude (max_error < 1e-9) and all basis-state spot checks hit",
        passed=worst < 1e-9 and failures == 0,
    )


def verify_xor_binding_consistency(report):
    """gate / QRAM dual-binding consistency of an abstract declaration (the parameterized cross-check of the validation-matrix V4 gap)."""
    aw, dw = 3, 2
    table = _pseudo_table(aw, dw, 31)
    slot = abstract_database("Mem", aw, dw)
    program = _db_program(slot.operation, superpose=("address", "data"), name="open_sup")
    gate_bound = bind(program, {"Mem": gate_database(aw, dw, table).operation})
    qram_bound = bind(
        program, {"Mem": Binding(qram_database(aw, dw).operation, {"table": "mem"})}
    )
    memory = {"mem": table}
    expected = _xor_expected(aw, dw, table)
    gate_ref = reference(gate_bound)
    qram_ref = reference(qram_bound, memory)
    worst = max(
        amplitude_error(gate_ref, expected),
        amplitude_error(qram_ref, expected),
        amplitude_error(rir_pysparq(gate_bound), expected),
        amplitude_error(rir_pysparq(qram_bound, memory), expected),
    )
    gate_vector = originir_ext(gate_bound)
    qram_vector = originir_ext(qram_bound, memory)
    expected_vector = amplitudes_to_statevector(expected, [aw, dw])
    worst = max(
        worst,
        statevector_error(gate_vector, expected_vector),
        statevector_error(qram_vector, expected_vector),
        statevector_error(
            [complex(a - b) for a, b in zip(gate_vector, qram_vector, strict=True)],
            [0j] * len(expected_vector),
        ),
    )
    report.case(
        "xor-abstract-binding-consistency-3x2",
        paths=["reference", "rir-pysparq", "originir-ext"],
        parameters={"address_width": aw, "data_width": dw, "bindings": ["gate", "qram"]},
        metrics={"max_error": worst},
        criterion="both bindings of the same open declaration pairwise match the classical permutation (max_error < 1e-9)",
        passed=worst < 1e-9,
    )


# ---------------------------------------------------------------------------
# Group B: Bernstein-Vazirani (named in the paper, prioritized)
# ---------------------------------------------------------------------------


def _bv_secrets(width):
    mask = (1 << width) - 1
    alternating = sum(1 << i for i in range(0, width, 2))
    hashed = 0x9E3779B97F4A7C15 & mask
    return sorted({1, mask, alternating, hashed})


def verify_bv_recovery(report):
    """3-8 bit multiple scales: all secret strings x bias instances; the four-path recovery probability is exactly 1."""
    for width in range(3, 9):
        secrets = _bv_secrets(width)
        worst = {"p_success": 1.0, "tvd": 0.0, "max_error": 0.0}
        instances = 0
        for secret in secrets:
            for bias in (0, 1):
                instances += 1
                program = bernstein_vazirani(
                    affine_boolean_oracle(width, secret, bias=bias)
                ).program()
                expected = _bv_expected(secret, bias)
                states = {
                    "reference": reference(program),
                    "rir-pysparq": rir_pysparq(program),
                    "adapter-pysparq": adapter_pysparq(program),
                }
                for state in states.values():
                    worst["max_error"] = max(
                        worst["max_error"], amplitude_error(state, expected)
                    )
                    marginal = _input_marginal(probabilities(state))
                    worst["p_success"] = min(worst["p_success"], marginal.get(secret, 0.0))
                    worst["tvd"] = max(worst["tvd"], tvd(marginal, {secret: 1.0}))
                vector = originir_ext(program)
                worst["max_error"] = max(
                    worst["max_error"],
                    statevector_error(
                        vector, amplitudes_to_statevector(expected, [width, 1])
                    ),
                )
                worst["p_success"] = min(
                    worst["p_success"],
                    _origin_input_probability(vector, width, secret),
                )
        report.case(
            f"bv-recovery-w{width}",
            paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
            parameters={"width": width, "secrets": secrets, "instances": instances},
            metrics={
                "success_probability": worst["p_success"],
                "tvd": worst["tvd"],
                "max_error": worst["max_error"],
            },
            criterion="all instances read out the secret string (success_probability > 1 - 1e-9, tvd < 1e-9)",
            passed=worst["p_success"] > 1 - 1e-9 and worst["tvd"] < 1e-9,
        )


def verify_bv_oracle_truth_table(report):
    """Truth-table exhaustive check of affine_boolean_oracle itself: f(x) = s.x XOR c compared branch by branch."""
    for width, secret, bias in ((4, 0b1011, 1), (8, 0xA5, 0)):
        oracle = affine_boolean_oracle(width, secret, bias=bias)
        table = [_dot(secret, x) ^ bias for x in range(1 << width)]
        program = _db_program(
            oracle.operation, superpose=("address", "data"), name=f"bv_tab_{width}"
        )
        expected = _xor_expected(width, 1, table)
        worst = amplitude_error(reference(program), expected)
        worst = max(worst, amplitude_error(rir_pysparq(program), expected))
        vector = originir_ext(program)
        worst = max(
            worst,
            statevector_error(vector, amplitudes_to_statevector(expected, [width, 1])),
        )
        report.case(
            f"bv-oracle-truth-table-w{width}",
            paths=["reference", "rir-pysparq", "originir-ext"],
            parameters={"width": width, "secret": secret, "bias": bias},
            metrics={"max_error": worst},
            criterion="oracle superposition branches match s.x XOR c amplitude by amplitude (max_error < 1e-9)",
            passed=worst < 1e-9,
        )


def verify_bv_open_binding(report):
    """Open declaration -> gate/QRAM binding end-to-end recovery (the open-oracle route described in the docs)."""
    width = 4
    for secret, bias in ((0b1011, 1), (0b0110, 0)):
        slot = abstract_database("BVFunction", width, 1)
        program = bernstein_vazirani(slot).program()
        table = [_dot(secret, a) ^ bias for a in range(1 << width)]
        bound = {
            "gate": bind(program, {"BVFunction": gate_database(width, 1, table).operation}),
            "qram": bind(
                program,
                {"BVFunction": Binding(qram_database(width, 1).operation, {"table": "f"})},
            ),
        }
        expected = _bv_expected(secret, bias)
        p_min, deviation = 1.0, 0.0
        states = {}
        for label, prog in bound.items():
            memory = {"f": table} if label == "qram" else None
            for tag, runner in (
                ("reference", reference),
                ("rir-pysparq", rir_pysparq),
            ):
                state = runner(prog, memory)
                states[(label, tag)] = state
                p_min = min(
                    p_min, _input_marginal(probabilities(state)).get(secret, 0.0)
                )
                deviation = max(deviation, amplitude_error(state, expected))
            vector = originir_ext(prog, memory)
            p_min = min(p_min, _origin_input_probability(vector, width, secret))
            deviation = max(
                deviation,
                statevector_error(
                    vector, amplitudes_to_statevector(expected, [width, 1])
                ),
            )
        deviation = max(
            deviation, amplitude_error(states[("gate", "rir-pysparq")], states[("qram", "rir-pysparq")])
        )
        report.case(
            f"bv-open-bind-w{width}-s{secret}b{bias}",
            paths=["reference", "rir-pysparq", "originir-ext"],
            parameters={
                "width": width,
                "secret": secret,
                "bias": bias,
                "bindings": ["gate", "qram"],
            },
            metrics={"success_probability": p_min, "max_error": deviation},
            criterion="both bindings recover with probability 1 and match the closed-form state (p > 1 - 1e-9, max_error < 1e-9)",
            passed=p_min > 1 - 1e-9 and deviation < 1e-9,
        )


# ---------------------------------------------------------------------------
# Group C: Deutsch-Jozsa
# ---------------------------------------------------------------------------


def _dj_functions(width):
    size = 1 << width
    half = size // 2
    return {
        "const0": [0] * size,
        "const1": [1] * size,
        "balanced_parity": [x.bit_count() & 1 for x in range(size)],
        "balanced_msb": [(x >> (width - 1)) & 1 for x in range(size)],
        "balanced_threshold": [1 if x < half else 0 for x in range(size)],
    }


def verify_dj_decision(report):
    """Constant/balanced function families: P(input==0) must be exactly 1 (constant) or 0 (balanced)."""
    for width in (2, 3, 4, 5):
        worst_p0, decision_errors, worst_state = 0.0, 0, 0.0
        for label, table in _dj_functions(width).items():
            constant = label.startswith("const")
            program = deutsch_jozsa(gate_database(width, 1, table)).program()
            expected_p0 = 1.0 if constant else 0.0
            states = {
                "reference": reference(program),
                "rir-pysparq": rir_pysparq(program),
            }
            for state in states.values():
                p0 = _input_marginal(probabilities(state)).get(0, 0.0)
                worst_p0 = max(worst_p0, abs(p0 - expected_p0))
                if (p0 > 0.5) != constant:
                    decision_errors += 1
            vector = originir_ext(program)
            p0 = _origin_input_probability(vector, width, 0)
            worst_p0 = max(worst_p0, abs(p0 - expected_p0))
            if (p0 > 0.5) != constant:
                decision_errors += 1
            # Constant and linearly balanced (parity = s.x with s all ones) have closed-form final states; phase-sensitive comparison
            if label == "const0":
                closed = _bv_expected(0, 0)
            elif label == "const1":
                closed = _bv_expected(0, 1)
            elif label == "balanced_parity":
                closed = _bv_expected((1 << width) - 1, 0)
            else:
                continue
            for state in states.values():
                worst_state = max(worst_state, amplitude_error(state, closed))
            worst_state = max(
                worst_state,
                statevector_error(vector, amplitudes_to_statevector(closed, [width, 1])),
            )
        report.case(
            f"dj-decision-w{width}",
            paths=["reference", "rir-pysparq", "originir-ext"],
            parameters={"width": width, "functions": sorted(_dj_functions(width))},
            metrics={
                "p_zero_max_error": worst_p0,
                "decision_errors": decision_errors,
                "closed_form_max_error": worst_state,
            },
            criterion="all decisions correct with P(input=0) deviation < 1e-9 (closed-form cases amplitude by amplitude < 1e-9)",
            passed=decision_errors == 0 and worst_p0 < 1e-9 and worst_state < 1e-9,
        )


def _network_database(net, input_name, output_name, width):
    """Wrap a BooleanNetwork compilation product into the XOR database interface (output has XOR-copy semantics)."""
    net_op = net.operation()
    b = Builder("net_db", {"address": Bits(width), "data": Bits(1)})
    b.call(net_op, **{input_name: b["address"], output_name: b["data"]})
    return XorDatabase(annotate(b.finish(), "database_xor", implementation="boolean_network"))


def verify_dj_boolean_network(report):
    """DJ end to end with a BooleanNetwork-compiled oracle + circuit semantics exhaustively vs classical evaluation."""
    parity_net = BooleanNetwork()
    xs = parity_net.input("x", 3)
    parity_net.outputs["out"] = [parity_net.xor(parity_net.xor(xs[0], xs[1]), xs[2])]
    const_net = BooleanNetwork()
    const_net.input("x", 3)
    const_net.outputs["out"] = BooleanNetwork.const(1, 1)
    for label, net, constant in (
        ("parity3", parity_net, False),
        ("const3", const_net, True),
    ):
        database = _network_database(net, "x", "out", 3)
        table = [net.evaluate(x=a)["out"] for a in range(8)]
        # Circuit-semantics exhaustive sweep: one superposition run against all 16 branches of net.evaluate
        sweep = _db_program(database.operation, superpose=("address", "data"), name="net_sw")
        expected = _xor_expected(3, 1, table)
        worst = amplitude_error(rir_pysparq(sweep), expected)
        vector = originir_ext(sweep)
        worst = max(
            worst, statevector_error(vector, amplitudes_to_statevector(expected, [3, 1]))
        )
        # DJ decision
        program = deutsch_jozsa(database).program()
        p0 = _input_marginal(probabilities(rir_pysparq(program))).get(0, 0.0)
        expected_p0 = 1.0 if constant else 0.0
        p0_error = abs(p0 - expected_p0)
        p0 = _origin_input_probability(originir_ext(program), 3, 0)
        p0_error = max(p0_error, abs(p0 - expected_p0))
        report.case(
            f"dj-boolean-network-{label}",
            paths=["rir-pysparq", "originir-ext"],
            parameters={"width": 3, "constant": constant, "branches": 16},
            metrics={"max_error": worst, "p_zero_error": p0_error},
            criterion="network circuit matches classical evaluation branch by branch and the DJ decision is correct (error < 1e-9)",
            passed=worst < 1e-9 and p0_error < 1e-9,
        )


# ---------------------------------------------------------------------------
# Group D: Simon
# ---------------------------------------------------------------------------


def _simon_table(n, s):
    """Truth table of a two-to-one linear function with period exactly s (independently constructed + self-checked by classical exhaustion)."""
    pivot = (s & -s).bit_length() - 1

    def f(x):
        z = x ^ (((x >> pivot) & 1) * s)
        return (z & ((1 << pivot) - 1)) | ((z >> (pivot + 1)) << pivot)

    table = [f(x) for x in range(1 << n)]
    assert all(table[x] == table[x ^ s] for x in range(1 << n))
    assert len(set(table)) == 1 << (n - 1)
    return table


def _simon_cnot_oracle(n, s):
    """CNOT realization of a linear Simon function (replaces the truth table at wide widths to reduce the gate count)."""
    pivot = (s & -s).bit_length() - 1
    b = Builder(f"simon_linear_{n}_{s}", {"address": Bits(n), "data": Bits(n - 1)})
    for j in range(n - 1):
        m = j if j < pivot else j + 1
        b.xor(b["address"][m], b["data"][j])
        if (s >> m) & 1:
            b.xor(b["address"][pivot], b["data"][j])
    return XorDatabase(annotate(b.finish(), "database_xor", implementation="cnot_linear"))


def _simon_expected(n, s, table):
    """Exact joint distribution (phases included) of the Simon sampling circuit: 2^{2n-2} equal-amplitude basis states in the support."""
    pivot = (s & -s).bit_length() - 1
    representative = {}
    for x, value in enumerate(table):
        if not (x >> pivot) & 1:
            representative[value] = x
    amplitude = 1.0 / (1 << (n - 1))
    expected = {}
    for y in range(1 << n):
        if _dot(y, s):
            continue
        for value, x0 in representative.items():
            expected[(y, value)] = amplitude * (1 if _dot(x0, y) == 0 else -1)
    return expected


def verify_simon(report):
    """Phase-sensitive comparison of the sampling distribution + elimination recovering the period (3-8 bit)."""
    for n in range(3, 9):
        use_cnot = n >= 7  # the truth table's gate count grows as 2^n; wide widths use the CNOT linear realization
        secrets = sorted({1, (1 << n) - 1, (0x9E3779B97F4A7C15 & ((1 << n) - 1)) | 1})
        worst_amp, worst_tvd, recoveries, attempts = 0.0, 0.0, 0, 0
        case_paths = ["rir-pysparq"] if use_cnot else ["reference", "rir-pysparq"]
        if not use_cnot and n <= 5:
            case_paths.append("originir-ext")  # dense small instances go through UniQC full amplitude
        for s in secrets:
            table = _simon_table(n, s)
            oracle = (
                _simon_cnot_oracle(n, s) if use_cnot else gate_database(n, n - 1, table)
            )
            program = simon_sample(oracle).program()
            expected = _simon_expected(n, s, table)
            expected_probs = probabilities(expected)
            states = {"rir-pysparq": rir_pysparq(program, max_states=1 << 16)}
            if not use_cnot:
                states["reference"] = reference(program, max_states=1 << 16)
            for state in states.values():
                worst_amp = max(worst_amp, amplitude_error(state, expected))
                probs = probabilities(state)
                worst_tvd = max(worst_tvd, tvd(probs, expected_probs))
                if any(_dot(y, s) for (y, _value) in probs):
                    worst_tvd = math.inf  # support escapes the subspace orthogonal to s
                # End-to-end elimination recovery using this path's own sampling support
                support = sorted({y for (y, _value), p in probs.items() if p > 1e-15})
                samples = _gf2_independent(support, n - 1)
                attempts += 1
                if simon_nullspace(samples, n) == (s,):
                    recoveries += 1
            if "originir-ext" in case_paths:
                vector = originir_ext(program)
                expected_vector = amplitudes_to_statevector(expected, [n, n - 1])
                worst_amp = max(worst_amp, statevector_error(vector, expected_vector))
        report.case(
            f"simon-sampling-recovery-w{n}",
            paths=case_paths,
            parameters={
                "width": n,
                "secrets": secrets,
                "oracle": "cnot_linear" if use_cnot else "gate_truth_table",
                "support_states": 1 << (2 * n - 2),
            },
            metrics={
                "max_error": worst_amp,
                "tvd": worst_tvd,
                "recoveries": recoveries,
                "recovery_attempts": attempts,
            },
            criterion="joint distribution matches amplitude by amplitude and every sampled elimination recovers the period (recoveries == attempts)",
            passed=worst_amp < 1e-9 and worst_tvd < 1e-9 and recoveries == attempts,
        )


def verify_simon_rank_deficiency(report):
    """With insufficient samples the null space stays multidimensional and the period lies in its span (informative metric)."""
    n, s = 4, 0b1011
    table = _simon_table(n, s)
    program = simon_sample(gate_database(n, n - 1, table)).program()
    state = rir_pysparq(program)
    support = sorted({y for (y, _v), p in probabilities(state).items() if p > 1e-15})
    samples = _gf2_independent(support, n - 2)  # one independent sample short
    basis = simon_nullspace(samples, n)
    span = {0}
    for vector in basis:
        span |= {x ^ vector for x in list(span)}
    in_span = s in span
    report.case(
        "simon-rank-deficiency-w4",
        paths=["rir-pysparq"],
        parameters={"width": n, "secret": s, "independent_samples": len(samples)},
        metrics={"nullspace_dim": len(basis), "secret_in_span": in_span},
        criterion="the null space has exactly one extra dimension and s lies in its span",
        passed=len(basis) == 2 and in_span,
    )


def run():
    report = Report(
        "oracles",
        "Exhaustive and end-to-end numerical validation of the three XOR database bindings and the BV/DJ/Simon query algorithms "
        "on four real backend paths (BV covers 3-8 bit multiple scales on dual paths).",
    )
    verify_xor_gate_basis(report)
    verify_xor_gate_superposition(report)
    verify_xor_qram(report)
    verify_xor_binding_consistency(report)
    verify_bv_recovery(report)
    verify_bv_oracle_truth_table(report)
    verify_bv_open_binding(report)
    verify_dj_decision(report)
    verify_dj_boolean_network(report)
    verify_simon(report)
    verify_simon_rank_deficiency(report)
    report.write()
    return report


if __name__ == "__main__":
    run()
