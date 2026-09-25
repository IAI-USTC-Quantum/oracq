"""Publication-grade numerical validation of search and quantum-walk algorithms.

Covers search.py (Grover success-rate curves, amplitude amplification),
walks.py (Hadamard-coin walks on a cycle), graph_walks.py (adjacency-oracle
end to end, Szegedy walk, MNRS search, classical Markov-chain reference
routines), and provides end-to-end numerics for adjacent entries touched by
the same documentation group: transforms.qubitization_walk (Chebyshev
rotation), qsvt.fixed_point_search (YLC fixed-point guarantee),
estimation.amplitude_estimation (quantum-counting readout), and
qlss.costa_walk (operator unitarity).

All classical references are constructed independently: closed-form formulas
sin^2((2k+1)theta) / YLC P_S / the Dirichlet kernel, and numpy-assembled
walk-space reflection operators and Markov-chain linear systems, without
reusing internal helpers of the modules under test.

Run: PYTHONPATH=src <python with pysparq+uniqc> tests/verification/verify_search_walks.py
"""

from __future__ import annotations

import math

from harness import (
    Report,
    adapter_pysparq,
    amplitude_error,
    originir_ext,
    originir_unitary,
    reference,
    rir_pysparq,
    superposition_program,
    tvd,
)

from oracq import Bits, Builder
from oracq.algorithms.input_model.oracles import invoke, resources_for

ALL_PATHS = ("reference", "rir-pysparq", "adapter-pysparq", "originir-ext")


def _register_amplitudes(vector, widths):
    """OriginIR state vector -> register-tuple sparse dictionary (key order matching the reference path)."""
    result = {}
    for index, amplitude in enumerate(vector):
        if amplitude:
            key = ()
            for w in widths:
                key += (index & ((1 << w) - 1),)
                index >>= w
            result[key] = amplitude
    return result


def _execute(path, program, memory=None):
    """Execute per path name and unify the output into a register-tuple -> amplitude dictionary."""
    if path == "reference":
        return reference(program, memory)
    if path == "rir-pysparq":
        return rir_pysparq(program, memory)
    if path == "adapter-pysparq":
        return adapter_pysparq(program, memory)
    if path == "originir-ext":
        if memory is not None:
            raise ValueError("the OriginIR-ext backend does not accept QRAM data")
        widths = [r.type.width for r in program.main.registers]
        return _register_amplitudes(originir_ext(program), widths)
    raise ValueError(path)


def _marginal(state, index, value):
    """Marginal probability of the target register taking value."""
    return sum(abs(a) ** 2 for key, a in state.items() if key[index] == value)


# ---------------------------------------------------------------------------
# Grover and amplitude amplification (search.py)
# ---------------------------------------------------------------------------


def _grover_case(report, name, phase_oracle, width, marked, iterations, note):
    """Success rates vs sin^2((2k+1)theta) over several iteration counts, with the two-level distribution checked per basis state."""
    from oracq.algorithms.common.search import grover

    size, t = 1 << width, len(marked)
    theta = math.asin(math.sqrt(t / size))
    per_iteration = []
    prob_error = dist_error = 0.0
    for k in iterations:
        program = grover(phase_oracle, width, iterations=k).operation.program()
        theory = math.sin((2 * k + 1) * theta) ** 2
        deviation = 0.0
        for path in ALL_PATHS:
            state = _execute(path, program)
            deviation = max(deviation, abs(_marked_prob(state, marked) - theory))
            for value in range(size):
                expect = theory / t if value in marked else (1 - theory) / (size - t)
                dist_error = max(dist_error, abs(_marginal(state, 0, value) - expect))
        prob_error = max(prob_error, deviation)
        per_iteration.append({"k": k, "theory": theory, "max_deviation": deviation})
    report.case(
        name,
        paths=list(ALL_PATHS),
        parameters={"width": width, "marked": list(marked), "iterations": list(iterations), "oracle": note},
        metrics={
            "max_prob_error": prob_error,
            "max_distribution_error": dist_error,
            "per_iteration": per_iteration,
        },
        criterion="each path's marked probability and per-basis-state distribution match the two-level sin^2((2k+1)theta) formulas (error < 1e-9)",
        passed=prob_error < 1e-9 and dist_error < 1e-9,
    )


def _marked_prob(state, marked):
    return sum(abs(a) ** 2 for key, a in state.items() if key[0] in marked)


def verify_grover_phase_marks(report):
    from oracq.algorithms.input_model.oracles import phase_marks

    _grover_case(
        report,
        "grover-phase-marks-n3-t1",
        phase_marks(3, (5,)),
        3,
        (5,),
        (0, 1, 2, 3, 4),  # includes the oscillating descent past the optimal iteration count
        "phase_marks",
    )


def verify_grover_xor_database(report):
    from oracq.algorithms.common.search import phase_from_database
    from oracq.algorithms.input_model.oracles import gate_database

    # t=2/8 gives theta=pi/6: k=1 amplifies exactly to 1, k=2 falls back to 1/4, testing the overshoot segment
    database = gate_database(3, 1, {5: 1, 6: 1})
    _grover_case(
        report,
        "grover-xor-database-n3-t2",
        phase_from_database(database),
        3,
        (5, 6),
        (0, 1, 2, 3),
        "phase_from_database(gate_database)",
    )


def verify_amplify_success(report):
    from oracq.algorithms.common.search import amplify_success
    from oracq.algorithms.input_model.oracles import StateOracle, annotate

    # Initial probability cos^2 theta = sin^2(pi/8) of the signal==0 success subspace, theta_a = pi/8
    theta = 3 * math.pi / 8
    b = Builder("amplify_source", {"target": Bits(2), "signal": Bits(1)})
    b.h(b["target"])
    b.ry(b["signal"][0], 2 * theta)
    state_oracle = StateOracle(annotate(b.finish(), "unitary"))
    theta_a = math.pi / 8
    per_iteration = []
    worst = 0.0
    for k in (0, 1, 2, 3):
        program = amplify_success(state_oracle, iterations=k).operation.program()
        theory = math.sin((2 * k + 1) * theta_a) ** 2
        deviation = 0.0
        for path in ALL_PATHS:
            state = _execute(path, program)
            deviation = max(deviation, abs(_marginal(state, 1, 0) - theory))
        worst = max(worst, deviation)
        per_iteration.append({"k": k, "theory": theory, "max_deviation": deviation})
    report.case(
        "amplify-success-curve",
        paths=list(ALL_PATHS),
        parameters={"initial_probability": math.cos(theta) ** 2, "iterations": [0, 1, 2, 3]},
        metrics={"max_prob_error": worst, "per_iteration": per_iteration},
        criterion="zero-signal success probability vs sin^2((2k+1)theta_a) (error < 1e-9, overshoot and fallback included)",
        passed=worst < 1e-9,
    )


# ---------------------------------------------------------------------------
# Coined walks on a cycle (walks.py)
# ---------------------------------------------------------------------------


def _coined_cycle_reference(width, steps):
    """Independent numpy semantics reference: H coin then a conditional +/-1 shift on coin 0/1, starting from |0,0>."""
    size = 1 << width
    state = {(0, 0): 1.0 + 0j}
    for _ in range(steps):
        coined = {}
        for (pos, coin), amp in state.items():
            sign = 1.0 if coin == 0 else -1.0
            coined[(pos, 0)] = coined.get((pos, 0), 0j) + amp / math.sqrt(2)
            coined[(pos, 1)] = coined.get((pos, 1), 0j) + sign * amp / math.sqrt(2)
        shifted = {}
        for (pos, coin), amp in coined.items():
            new = (pos + 1) % size if coin == 0 else (pos - 1) % size
            shifted[(new, coin)] = shifted.get((new, coin), 0j) + amp
        state = shifted
    return state


def _classical_cycle_distribution(width, steps):
    """Position distribution of a classical symmetric random walk (p=1/2) on the cycle: binomial mod N."""
    size = 1 << width
    distribution = [0.0] * size
    for right in range(steps + 1):
        prob = math.comb(steps, right) / 2**steps
        distribution[(2 * right - steps) % size] += prob
    return distribution


def _ring_distance(position, size):
    return min(position, size - position)


def _cycle_walk_case(report, width, steps_range):
    from oracq.algorithms.common.walks import cycle_walk

    size = 1 << width
    amp_error = tv_distance = 0.0
    for steps in steps_range:
        program = cycle_walk(width, steps=steps).program()
        expected = _coined_cycle_reference(width, steps)
        # With steps=0 the program contains no gates and UniQC has no qubit mapping for zero-gate circuits (max of an
        # empty sequence), so the zero-step case runs only the register-level paths; this is an OriginIR-backend
        # limitation on empty programs, not an issue of the module under test.
        paths = ALL_PATHS if steps else tuple(p for p in ALL_PATHS if p != "originir-ext")
        for path in paths:
            state = _execute(path, program)
            amp_error = max(amp_error, amplitude_error(state, expected))
            actual_probs = {pos: _marginal(state, 0, pos) for pos in range(size)}
            expected_probs = {pos: 0.0 for pos in range(size)}
            for (pos, _coin), amp in expected.items():
                expected_probs[pos] += abs(amp) ** 2
            tv_distance = max(tv_distance, tvd(actual_probs, expected_probs))
    # Informative metric: quantum ballistic transport vs classical diffusion (mean ring distance at the final step)
    final = max(steps_range)
    quantum_state = _execute("reference", cycle_walk(width, steps=final).program())
    quantum_mean = sum(
        _ring_distance(pos, size) * _marginal(quantum_state, 0, pos) for pos in range(size)
    )
    classical = _classical_cycle_distribution(width, final)
    classical_mean = sum(_ring_distance(pos, size) * classical[pos] for pos in range(size))
    report.case(
        f"cycle-walk-w{width}",
        paths=list(ALL_PATHS),
        parameters={"width": width, "steps": list(steps_range)},
        metrics={
            "max_amplitude_error": amp_error,
            "max_position_tvd": tv_distance,
            f"mean_ring_distance_quantum_s{final}": quantum_mean,
            f"mean_ring_distance_classical_s{final}": classical_mean,
        },
        criterion="full amplitudes vs the independent coin-walk reference (error < 1e-9); mean ring distance is an informative metric",
        passed=amp_error < 1e-9 and tv_distance < 1e-9,
    )


# ---------------------------------------------------------------------------
# Graph adjacency and Szegedy/MNRS (graph_walks.py)
# ---------------------------------------------------------------------------


def _cycle_table(n):
    """Alternating edge-coloring neighbor table of an even cycle (satisfying the involution N(N(v,j),j)=v)."""
    return [
        [(v + 1) % n if v % 2 == 0 else (v - 1) % n, (v - 1) % n if v % 2 == 0 else (v + 1) % n]
        for v in range(n)
    ]


def _complete_table(n):
    """Complete graph padded with self-loops: N(v,j)=j, D=N."""
    return [list(range(n)) for _ in range(n)]


HYPERCUBE_Q3 = [[v ^ 1, v ^ 2, v ^ 4, v] for v in range(8)]


def _szegedy_reference(neighbors):
    """Independently assemble the walk-space operators: returns (dimension layout, W, initial state, marked phase diagonal)."""
    import numpy as np

    n, d = len(neighbors), len(neighbors[0])
    v = max(1, (n - 1).bit_length())
    g = (d - 1).bit_length()
    dim = n * n * d

    def index(current, peer, j):
        return current | (peer << v) | (j << (2 * v))

    a_states = np.zeros((dim, n))
    b_states = np.zeros((dim, n))
    for vertex in range(n):
        for j in range(d):
            a_states[index(vertex, neighbors[vertex][j], j), vertex] = 1.0 / math.sqrt(d)
            b_states[index(neighbors[vertex][j], vertex, j), vertex] = 1.0 / math.sqrt(d)
    identity = np.eye(dim)
    walk = (2.0 * b_states @ b_states.T - identity) @ (2.0 * a_states @ a_states.T - identity)
    setup = a_states.sum(axis=1) / math.sqrt(n)
    return {"v": v, "g": g, "dim": dim, "index": index, "walk": walk, "setup": setup}


def _mnrs_reference(neighbors, marked, steps):
    """Independent numpy reference of (M*W)^steps acting on the initial state, MNRS style."""
    import numpy as np

    ref = _szegedy_reference(neighbors)
    n = len(neighbors)
    phase = np.ones(ref["dim"])
    for current in marked:
        for peer in range(n):
            for j in range(len(neighbors[0])):
                phase[ref["index"](current, peer, j)] = -1.0
    state = ref["setup"].copy()
    for _ in range(steps):
        state = phase * (ref["walk"] @ state)
    return ref, state


def _walk_space_dict(ref, vector):
    """numpy walk-space vector -> (current, peer, index) sparse amplitude dictionary."""
    n = 1 << ref["v"]
    d = 1 << ref["g"] if ref["g"] else 1
    result = {}
    for current in range(n):
        for peer in range(n):
            for j in range(d):
                amplitude = vector[ref["index"](current, peer, j)]
                if amplitude:
                    result[(current, peer, j)] = complex(amplitude)
    return result


def _embed_walk_dict(ref, amplitudes):
    """Fold the (current, peer, index) keys into the target integer value, aligning with the driver's (target, work) keys."""
    return {
        (ref["index"](current, peer, j), 0): amplitude
        for (current, peer, j), amplitude in amplitudes.items()
    }


def verify_adjacency_superposition(report):
    """Adjacency oracle end to end: vertex/index full superposition exhausts all 32 queries in one run."""
    from oracq.algorithms.input_model.graph_walks import gate_adjacency

    oracle = gate_adjacency(HYPERCUBE_Q3)
    program = superposition_program(oracle.operation, ["vertex", "index"])
    uniform = 1.0 / math.sqrt(8 * 4)
    expected = {(v, j, HYPERCUBE_Q3[v][j]): uniform for v in range(8) for j in range(4)}
    worst = foreign = 0.0
    for path in ALL_PATHS:
        state = _execute(path, program)
        worst = max(worst, amplitude_error(state, expected))
        foreign = max(
            foreign,
            sum(abs(a) ** 2 for key, a in state.items() if key[2] != HYPERCUBE_Q3[key[0]][key[1]]),
        )
    report.case(
        "adjacency-superposition-hypercube",
        paths=list(ALL_PATHS),
        parameters={"vertices": 8, "degree": 4, "queries": 32},
        metrics={"max_amplitude_error": worst, "foreign_branch_weight": foreign},
        criterion="superposed queries match the XOR neighbor table branch by branch (amplitude error < 1e-9, zero weight off the table)",
        passed=worst < 1e-9 and foreign < 1e-18,
    )


def verify_szegedy_walk(report):
    """Szegedy walk step: both the unitary matrix and the evolution on the initial state vs the independent reflection-operator assembly."""
    import numpy as np

    from oracq.algorithms.input_model.graph_walks import (
        gate_adjacency,
        szegedy_setup,
        szegedy_walk,
    )

    neighbors = _cycle_table(8)
    ref = _szegedy_reference(neighbors)
    v, g = ref["v"], ref["g"]
    adjacency = gate_adjacency(neighbors)
    walk = szegedy_walk(adjacency)
    # Unitary level: to_matrix vs (2*Pi_B - I)(2*Pi_A - I)
    driver = Builder(
        "szegedy_unitary",
        {"current": Bits(v), "peer": Bits(v), "index": Bits(g)},
        resources_for(("walk", walk)),
    )
    invoke(driver, walk, "walk", current=driver["current"], peer=driver["peer"], index=driver["index"])
    unitary = originir_unitary(driver.finish().program())
    matrix_error = float(np.abs(unitary - ref["walk"]).max())
    # State level: 1/3 steps of evolution after setup, compared on all four paths
    setup = szegedy_setup(adjacency)
    state_error = 0.0
    for steps in (1, 3):
        b = Builder(
            f"szegedy_evolve_{steps}",
            {"target": Bits(2 * v + g), "work": Bits(0)},
            resources_for(("setup", setup.operation), ("walk", walk)),
        )
        invoke(b, setup.operation, "setup", target=b["target"], work=b["work"])
        with b.repeat(steps):
            invoke(
                b,
                walk,
                "walk",
                current=b["target"][:v],
                peer=b["target"][v : 2 * v],
                index=b["target"][2 * v :],
            )
        program = b.finish().program()
        vector = np.linalg.matrix_power(ref["walk"], steps) @ ref["setup"]
        expected = _embed_walk_dict(ref, _walk_space_dict(ref, vector))
        for path in ALL_PATHS:
            state_error = max(state_error, amplitude_error(_execute(path, program), expected))
    report.case(
        "szegedy-walk-cycle8",
        paths=[*ALL_PATHS, "originir-ext+to_matrix"],
        parameters={"vertices": 8, "degree": 2, "walk_qubits": 2 * v + g, "steps": [1, 3]},
        metrics={"matrix_max_error": matrix_error, "state_max_error": state_error},
        criterion="walk unitary and evolved states vs the independent reflection-operator assembly (error < 1e-9)",
        passed=matrix_error < 1e-9 and state_error < 1e-9,
    )


def _hitting_reference(transition, marked):
    """Independently solve (I - P_free) h = 1 (marked vertices fixed at 0)."""
    import numpy as np

    n = len(transition)
    free = [v for v in range(n) if v not in set(marked)]
    if not free:
        return [0.0] * n
    sub = np.array([[transition[a][b] for b in free] for a in free])
    solved = np.linalg.solve(np.eye(len(free)) - sub, np.ones(len(free)))
    hits = {v: 0.0 for v in marked}
    hits.update({free[i]: float(solved[i]) for i in range(len(free))})
    return [hits[v] for v in range(n)]


def _transition_reference(neighbors):
    d = len(neighbors[0])
    return [[row.count(u) / d for u in range(len(neighbors))] for row in neighbors]


def _mnrs_case(report, name, neighbors, marked, paths, memory=None):
    from oracq.algorithms.input_model.graph_walks import (
        gate_adjacency,
        quantum_walk_search,
        szegedy_setup,
        szegedy_walk,
    )

    n = len(neighbors)
    v = max(1, (n - 1).bit_length())
    from oracq.algorithms.input_model.oracles import phase_marks

    # Independent classical reference: hitting times determine the MNRS step count ceil(pi/4 * sqrt(H_avg))
    transition = _transition_reference(neighbors)
    hits = _hitting_reference(transition, marked)
    average_h = sum(hits) / (n - len(marked))
    steps = max(1, math.ceil(math.pi / 4 * math.sqrt(average_h)))
    adjacency = gate_adjacency(neighbors)
    program = quantum_walk_search(
        szegedy_setup(adjacency), szegedy_walk(adjacency), phase_marks(v, marked), steps
    ).program()
    ref, vector = _mnrs_reference(neighbors, marked, steps)
    expected = _embed_walk_dict(ref, _walk_space_dict(ref, vector))
    mask = (1 << v) - 1
    expected_prob = sum(
        abs(a) ** 2 for key, a in expected.items() if (key[0] & mask) in set(marked)
    )
    state_error = prob_error = 0.0
    for path in paths:
        state = _execute(path, program, memory=memory)
        state_error = max(state_error, amplitude_error(state, expected))
        measured_prob = sum(
            abs(a) ** 2 for key, a in state.items() if (key[0] & mask) in set(marked)
        )
        prob_error = max(prob_error, abs(measured_prob - expected_prob))
    report.case(
        name,
        paths=list(paths),
        parameters={
            "vertices": n,
            "degree": len(neighbors[0]),
            "marked": list(marked),
            "steps": steps,
            "memory": "QRAM neighbor table" if memory else None,
        },
        metrics={
            "marked_probability": expected_prob,
            "max_state_error": state_error,
            "marked_prob_error": prob_error,
            "classical_avg_hitting_time": average_h,
            "steps_over_sqrt_h": steps / math.sqrt(average_h),
        },
        criterion="search final state's full amplitudes vs the independent MNRS reference (error < 1e-9); steps/sqrt(H) vs pi/4 is an informative metric",
        passed=state_error < 1e-9 and prob_error < 1e-9,
    )


def verify_mnrs_search(report):
    _mnrs_case(
        report,
        "mnrs-search-complete-k4",
        _complete_table(4),
        (0,),
        ALL_PATHS,
    )
    _mnrs_case(
        report,
        "mnrs-search-hypercube-q3",
        HYPERCUBE_Q3,
        (0,),
        ALL_PATHS,
    )


def verify_mnrs_search_qram(report):
    """End-to-end search with QRAM adjacency; OriginIR-ext carries no QRAM resources, so only reference and PySparQ are used."""
    from oracq.algorithms.input_model.graph_walks import qram_adjacency

    neighbors = _complete_table(4)
    table = {v | (j << 2): neighbors[v][j] for v in range(4) for j in range(4)}
    memory = {"setup__adj__table": table, "walk__adj__table": table}
    from oracq.algorithms.input_model.graph_walks import (
        quantum_walk_search,
        szegedy_setup,
        szegedy_walk,
    )
    from oracq.algorithms.input_model.oracles import phase_marks

    transition = _transition_reference(neighbors)
    average_h = sum(_hitting_reference(transition, (0,))) / 3
    steps = max(1, math.ceil(math.pi / 4 * math.sqrt(average_h)))
    adjacency = qram_adjacency(2, 2)
    program = quantum_walk_search(
        szegedy_setup(adjacency), szegedy_walk(adjacency), phase_marks(2, (0,)), steps
    ).program()
    ref, vector = _mnrs_reference(neighbors, (0,), steps)
    expected = _embed_walk_dict(ref, _walk_space_dict(ref, vector))
    state_error = 0.0
    paths = ("reference", "rir-pysparq")
    for path in paths:
        state_error = max(
            state_error, amplitude_error(_execute(path, program, memory=memory), expected)
        )
    marked_probability = sum(
        abs(a) ** 2 for key, a in expected.items() if (key[0] & 3) == 0
    )
    report.case(
        "mnrs-search-qram-k4",
        paths=list(paths),
        parameters={
            "vertices": 4,
            "degree": 4,
            "steps": steps,
            "originir_excluded": "the OriginIR-ext circuit carries no QRAM resources and cannot execute QRAM programs",
        },
        metrics={
            "marked_probability": marked_probability,
            "max_state_error": state_error,
        },
        criterion="after QRAM binding, the search final state matches the independent MNRS reference (error < 1e-9)",
        passed=state_error < 1e-9,
    )


def verify_markov_helpers(report):
    """transition_matrix / hitting_times / suggest_steps vs independent Markov-chain references."""
    from oracq.algorithms.input_model.graph_walks import (
        hitting_times,
        suggest_steps,
        transition_matrix,
    )

    tables = {
        "complete_k4": _complete_table(4),
        "hypercube_q3": HYPERCUBE_Q3,
        "cycle8": _cycle_table(8),
        # Non-regular graph: the path 0-1-2-3 padded with self-loops to D=2
        "padded_path4": [[1, 0], [0, 2], [1, 3], [2, 3]],
    }
    transition_error = hitting_error = 0.0
    suggest_mismatch = 0
    for name, neighbors in tables.items():
        marked = (0,) if name != "padded_path4" else (3,)
        expected_transition = _transition_reference(neighbors)
        actual = transition_matrix(neighbors)
        transition_error = max(
            transition_error,
            max(abs(actual[v][u] - expected_transition[v][u]) for v in range(len(neighbors)) for u in range(len(neighbors))),
        )
        expected_hits = _hitting_reference(expected_transition, marked)
        actual_hits = hitting_times(actual, set(marked))
        hitting_error = max(
            hitting_error, max(abs(a - e) for a, e in zip(actual_hits, expected_hits, strict=True))
        )
        average_h = sum(expected_hits) / (len(neighbors) - len(marked))
        expected_steps = max(1, math.ceil(math.pi / 4 * math.sqrt(average_h)))
        if suggest_steps(actual, set(marked)) != expected_steps:
            suggest_mismatch += 1
    # Closed form on the cycle h(v) = d(N-d): an independent closed-form comparison
    cycle_hits = hitting_times(transition_matrix(_cycle_table(8)), {0})
    closed_form_error = max(
        abs(cycle_hits[k] - min(k, 8 - k) * (8 - min(k, 8 - k))) for k in range(8)
    )
    # Complete-graph Grover limit: H_avg = N/|M|, steps recover ceil(pi/4 sqrt(N/|M|))
    grover_limit = []
    for n in (4, 8, 16):
        transition = transition_matrix(_complete_table(n))
        for m in (1, 2):
            expected = max(1, math.ceil(math.pi / 4 * math.sqrt(n / m)))
            actual = suggest_steps(transition, set(range(m)))
            grover_limit.append({"n": n, "m": m, "expected": expected, "actual": actual})
            if actual != expected:
                suggest_mismatch += 1
    report.case(
        "markov-chain-helpers",
        paths=["classical"],
        parameters={"tables": sorted(tables), "grover_limit": grover_limit},
        metrics={
            "transition_max_error": transition_error,
            "hitting_max_error": hitting_error,
            "cycle_closed_form_error": closed_form_error,
            "suggest_mismatch": suggest_mismatch,
        },
        criterion="transition matrix / hitting times vs independent numpy solves (error < 1e-9), step formula agrees case by case",
        passed=transition_error < 1e-9
        and hitting_error < 1e-9
        and closed_form_error < 1e-9
        and suggest_mismatch == 0,
    )


# ---------------------------------------------------------------------------
# Adjacent entries touched by the same documentation group (transforms/qsvt/estimation/qlss)
# ---------------------------------------------------------------------------


def _householder_be(theta):
    """Hermitian (Householder-style) block encoding whose zero-signal block is cos theta: U = Ry(2*theta)*Z tensor I_target."""
    from oracq.algorithms.input_model.operators import BlockEncoding
    from oracq.algorithms.input_model.oracles import annotate

    b = Builder("householder_be", {"target": Bits(1), "signal": Bits(1)})
    b.z(b["signal"][0])
    b.ry(b["signal"][0], 2 * theta)
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))


def verify_qubitization_walk(report):
    """Zero-signal probability of k walk applications on a Hermitian BE vs Chebyshev T_k(x)^2, plus the unitary decomposition and spectrum."""
    import numpy as np

    from oracq.algorithms.common.transforms import qubitization_walk

    x = math.cos(math.pi / 6)
    walk = qubitization_walk(_householder_be(math.pi / 6))
    prob_error = 0.0
    per_step = []
    for k in (1, 2, 3, 4):
        b = Builder(
            f"qubitization_repeat_{k}",
            {"target": Bits(1), "signal": Bits(1)},
            resources_for(("walk", walk)),
        )
        with b.repeat(k):
            invoke(b, walk, "walk", target=b["target"], signal=b["signal"])
        program = b.finish().program()
        theory = math.cos(k * math.acos(x)) ** 2
        deviation = 0.0
        for path in ALL_PATHS:
            deviation = max(deviation, abs(_marginal(_execute(path, program), 1, 0) - theory))
        prob_error = max(prob_error, deviation)
        per_step.append({"k": k, "theory": theory, "max_deviation": deviation})
    # Unitary relation W == (2*Pi - I)U and the spectrum: each of the two 2D subspaces contributes e^{+/-i*pi/6}
    unitary_walk = originir_unitary(walk.program())
    unitary_be = originir_unitary(_householder_be(math.pi / 6).operation.program())
    projector = np.diag([1.0 if (i >> 1) == 0 else 0.0 for i in range(4)])
    relation_error = float(
        np.abs(unitary_walk - (2 * projector - np.eye(4)) @ unitary_be).max()
    )
    angles = sorted(float(np.angle(value)) for value in np.linalg.eigvals(unitary_walk))
    expected_angles = sorted([math.pi / 6, math.pi / 6, -math.pi / 6, -math.pi / 6])
    eigenphase_error = max(
        abs(a - e) for a, e in zip(angles, expected_angles, strict=True)
    )
    report.case(
        "qubitization-walk-chebyshev",
        paths=[*ALL_PATHS, "originir-ext+to_matrix"],
        parameters={"x": x, "steps": [1, 2, 3, 4]},
        metrics={
            "max_prob_error": prob_error,
            "unitary_relation_error": relation_error,
            "eigenphase_error": eigenphase_error,
            "per_step": per_step,
        },
        criterion="zero-signal probability vs T_k(x)^2, W=(2*Pi-I)U with rotation angles +/-arccos(x) (error < 1e-9)",
        passed=prob_error < 1e-9 and relation_error < 1e-9 and eigenphase_error < 1e-9,
    )


def _fixed_point_case(report, name, x, delta, degree):
    from oracq.algorithms.common.qsvt import fixed_point_search
    from oracq.algorithms.input_model.oracles import diagonal_block_encoding, gate_database

    # Scalar BE: diagonal block cos(angle/2) == x (both basis states identical)
    database = gate_database(1, 1, {0: 1, 1: 1})
    be = diagonal_block_encoding(database, angle_scale=2 * math.acos(x))
    program = fixed_point_search(be, delta, degree).operation.program()
    # Independent closed form: P_S(x) = 1 - delta^2 T_L^2(c*sqrt(1-x^2)), c = T_{1/L}(1/delta)
    c = math.cosh(math.acosh(1.0 / delta) / degree)
    u = c * math.sqrt(1.0 - x * x)
    chebyshev = math.cosh(degree * math.acosh(u)) if u > 1 else math.cos(degree * math.acos(u))
    theory = 1.0 - delta * delta * chebyshev * chebyshev
    measured = []
    for path in ALL_PATHS:
        state = _execute(path, program)
        measured.append(_marginal(state, 1, 0))
    error = max(abs(value - theory) for value in measured)
    threshold = math.sqrt(1.0 - 1.0 / (c * c))
    report.case(
        name,
        paths=list(ALL_PATHS),
        parameters={"x": x, "delta": delta, "degree": degree, "threshold": threshold},
        metrics={
            "success_probability": sum(measured) / len(measured),
            "ylc_closed_form": theory,
            "max_error": error,
            "bound_1_minus_delta2": 1.0 - delta * delta,
        },
        criterion="zero-signal success probability vs the YLC closed form (error < 1e-9); outside the threshold one should have P_S >= 1-delta^2",
        passed=error < 1e-9 and (x < threshold or theory >= 1.0 - delta * delta),
    )


def verify_fixed_point_search(report):
    _fixed_point_case(report, "fixed-point-search-above-threshold", math.cos(math.pi / 6), 0.3, 5)
    _fixed_point_case(report, "fixed-point-search-below-threshold", 0.3, 0.3, 5)


def _qae_reference_distribution(a, precision):
    """Independent reference of the QAE phase readout: Dirichlet-kernel superposition at eigenphases +/-theta/pi."""
    grid = 1 << precision
    phi = math.asin(math.sqrt(a)) / math.pi

    def kernel(delta):
        if abs(delta) < 1e-15:
            return 1.0
        return math.sin(math.pi * grid * delta) / (grid * math.sin(math.pi * delta))

    return {
        y: 0.5 * (kernel(y / grid - phi) ** 2 + kernel(y / grid + phi) ** 2)
        for y in range(grid)
    }


def verify_quantum_counting(report):
    """Quantum counting with n=3 and 3/8 marked: phase distribution vs the Dirichlet kernel; readout of t-hat."""
    from oracq.algorithms.common.estimation import amplitude_estimation
    from oracq.algorithms.input_model.oracles import uniform_state

    n, marked, precision = 3, (1, 5, 7), 4
    size, t, grid = 1 << n, len(marked), 1 << precision
    program = amplitude_estimation(uniform_state(n), marked, precision=precision).program()
    expected = _qae_reference_distribution(t / size, precision)
    tv_distance = 0.0
    measured_dist = None
    for path in ALL_PATHS:
        state = _execute(path, program)
        actual = {y: _marginal(state, 2, y) for y in range(grid)}
        tv_distance = max(tv_distance, tvd(actual, expected))
        if path == "reference":
            measured_dist = actual
    peak = max(measured_dist, key=measured_dist.get)
    estimate = size * math.sin(math.pi * peak / grid) ** 2
    # QAE guarantee: probability of landing within one grid step of the truth >= 8/pi^2
    central_mass = sum(
        prob
        for y, prob in measured_dist.items()
        if abs(size * math.sin(math.pi * y / grid) ** 2 - t) <= 1.0 + 1e-12
    )
    report.case(
        "quantum-counting-n3-t3",
        paths=list(ALL_PATHS),
        parameters={"width": n, "marked": list(marked), "precision": precision},
        metrics={
            "phase_distribution_tvd": tv_distance,
            "estimated_count": estimate,
            "estimate_error": abs(estimate - t),
            "central_mass": central_mass,
            "qae_mass_bound": 8 / math.pi**2,
        },
        criterion="phase distribution TVD < 1e-9; |t-hat - t| <= 1 and central grid mass >= 8/pi^2",
        passed=tv_distance < 1e-9 and abs(estimate - t) <= 1.0 and central_mass >= 8 / math.pi**2,
    )


def verify_costa_walk_unitarity(report):
    """Unitarity and cross-backend agreement of the operator assembled by costa_walk (the kernel's physical channel upstream is labeled an unverified prototype)."""
    import numpy as np

    from oracq.algorithms.input_model.oracles import (
        basis_state,
        diagonal_block_encoding,
        gate_database,
    )
    from oracq.algorithms.qlss.qlss import costa_walk

    database = gate_database(1, 1, {0: 1, 1: 1})
    be = diagonal_block_encoding(database, angle_scale=math.pi / 3)
    walk = costa_walk(be, basis_state(1, 0), 0.5)
    program = walk.program()
    unitary = originir_unitary(program)
    deviation = float(
        np.abs(unitary.conj().T @ unitary - np.eye(unitary.shape[0])).max()
    )
    ref = _execute("reference", program)
    state_error = 0.0
    for path in ("rir-pysparq", "adapter-pysparq", "originir-ext"):
        state_error = max(state_error, amplitude_error(_execute(path, program), ref))
    report.case(
        "costa-walk-unitarity",
        paths=["originir-ext+to_matrix", "reference", "rir-pysparq", "adapter-pysparq"],
        parameters={"target_width": 1, "signal_width": 6, "fs": 0.5},
        metrics={"unitarity_deviation": deviation, "cross_path_state_error": state_error},
        criterion="max deviation of W^dagger W - I < 1e-9 and the four paths agree on the zero-input state",
        passed=deviation < 1e-9 and state_error < 1e-9,
    )


def run():
    report = Report(
        "search_walks",
        "Grover/amplitude-amplification success curves, coined cycle walks, adjacency-oracle and Szegedy/MNRS end-to-end walk numerics, "
        "plus the qubitization/fixed-point search/quantum counting/Costa walk entries touched by the same documentation group.",
    )
    verify_grover_phase_marks(report)
    verify_grover_xor_database(report)
    verify_amplify_success(report)
    _cycle_walk_case(report, 2, range(0, 9))
    _cycle_walk_case(report, 3, range(0, 9))
    verify_adjacency_superposition(report)
    verify_szegedy_walk(report)
    verify_mnrs_search(report)
    verify_mnrs_search_qram(report)
    verify_markov_helpers(report)
    verify_qubitization_walk(report)
    verify_fixed_point_search(report)
    verify_quantum_counting(report)
    verify_costa_walk_unitarity(report)
    report.write()
    return report


if __name__ == "__main__":
    run()
