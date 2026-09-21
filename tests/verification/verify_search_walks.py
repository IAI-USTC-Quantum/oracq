"""搜索与量子行走算法的论文级数值验证。

覆盖 search.py（Grover 成功率曲线、振幅放大）、walks.py（Hadamard coin 环上
行走）、graph_walks.py（邻接 oracle 端到端、Szegedy 行走、MNRS 搜索、经典
Markov 链参考例程），并为同组文档页涉及的相邻入口提供端到端数值：
transforms.qubitization_walk（Chebyshev 旋转）、qsvt.fixed_point_search
（YLC 不动点保证）、estimation.amplitude_estimation（量子计数读出）、
qlss.costa_walk（算符幺正性）。

经典参考全部独立构造：闭式公式 sin²((2k+1)θ) / YLC P_S / Dirichlet 核、
numpy 组装的行走空间反射算子与 Markov 链线性方程组，不复用被测模块内部
辅助函数。

运行：PYTHONPATH=src <含 pysparq+uniqc 的 python> tests/verification/verify_search_walks.py
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

from pyqecclang import Bits, Builder
from pyqecclang.algorithms.input_model.oracles import invoke, resources_for

ALL_PATHS = ("reference", "rir-pysparq", "adapter-pysparq", "originir-ext")


def _register_amplitudes(vector, widths):
    """OriginIR 态向量 → 寄存器元组稀疏字典（键序与 reference 路径一致）。"""
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
    """按路径名执行并把输出统一为寄存器元组 → 振幅的字典。"""
    if path == "reference":
        return reference(program, memory)
    if path == "rir-pysparq":
        return rir_pysparq(program, memory)
    if path == "adapter-pysparq":
        return adapter_pysparq(program, memory)
    if path == "originir-ext":
        if memory is not None:
            raise ValueError("OriginIR-ext 后端不接收 QRAM 数据")
        widths = [r.type.width for r in program.main.registers]
        return _register_amplitudes(originir_ext(program), widths)
    raise ValueError(path)


def _marginal(state, index, value):
    """目标寄存器取 value 的边缘概率。"""
    return sum(abs(a) ** 2 for key, a in state.items() if key[index] == value)


# ---------------------------------------------------------------------------
# Grover 与振幅放大（search.py）
# ---------------------------------------------------------------------------


def _grover_case(report, name, phase_oracle, width, marked, iterations, note):
    """多个迭代数上的成功率对照 sin²((2k+1)θ)，并逐基态核对两级分布。"""
    from pyqecclang.algorithms.common.search import grover

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
        criterion="各路径 marked 概率与逐基态分布对照 sin²((2k+1)θ) 两级公式（误差 < 1e-9）",
        passed=prob_error < 1e-9 and dist_error < 1e-9,
    )


def _marked_prob(state, marked):
    return sum(abs(a) ** 2 for key, a in state.items() if key[0] in marked)


def verify_grover_phase_marks(report):
    from pyqecclang.algorithms.input_model.oracles import phase_marks

    _grover_case(
        report,
        "grover-phase-marks-n3-t1",
        phase_marks(3, (5,)),
        3,
        (5,),
        (0, 1, 2, 3, 4),  # 含越过最优迭代数后的振荡下行段
        "phase_marks",
    )


def verify_grover_xor_database(report):
    from pyqecclang.algorithms.common.search import phase_from_database
    from pyqecclang.algorithms.input_model.oracles import gate_database

    # t=2/8 时 θ=π/6：k=1 精确放大到 1，k=2 回落 1/4，检验过冲段
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
    from pyqecclang.algorithms.common.search import amplify_success
    from pyqecclang.algorithms.input_model.oracles import StateOracle, annotate

    # 成功子空间 signal==0 的初态概率 cos²θ = sin²(π/8)，θ_a = π/8
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
        criterion="零信号成功概率对照 sin²((2k+1)θ_a)（误差 < 1e-9，含过冲回落）",
        passed=worst < 1e-9,
    )


# ---------------------------------------------------------------------------
# Coined 环上行走（walks.py）
# ---------------------------------------------------------------------------


def _coined_cycle_reference(width, steps):
    """独立的 numpy 语义参考：H coin 后按 0/1 条件 ±1 移位，从 |0,0> 出发。"""
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
    """经典对称随机游走（p=1/2）在环上的位置分布，二项分布 mod N。"""
    size = 1 << width
    distribution = [0.0] * size
    for right in range(steps + 1):
        prob = math.comb(steps, right) / 2**steps
        distribution[(2 * right - steps) % size] += prob
    return distribution


def _ring_distance(position, size):
    return min(position, size - position)


def _cycle_walk_case(report, width, steps_range):
    from pyqecclang.algorithms.common.walks import cycle_walk

    size = 1 << width
    amp_error = tv_distance = 0.0
    for steps in steps_range:
        program = cycle_walk(width, steps=steps).program()
        expected = _coined_cycle_reference(width, steps)
        # steps=0 时程序不含任何门，UniQC 对零门线路没有 qubit mapping（空序列 max），
        # 故零步情形只走寄存器级路径；这是 OriginIR 后端对空程序的限制，非被测模块问题。
        paths = ALL_PATHS if steps else tuple(p for p in ALL_PATHS if p != "originir-ext")
        for path in paths:
            state = _execute(path, program)
            amp_error = max(amp_error, amplitude_error(state, expected))
            actual_probs = {pos: _marginal(state, 0, pos) for pos in range(size)}
            expected_probs = {pos: 0.0 for pos in range(size)}
            for (pos, _coin), amp in expected.items():
                expected_probs[pos] += abs(amp) ** 2
            tv_distance = max(tv_distance, tvd(actual_probs, expected_probs))
    # 信息性指标：量子弹道输运 vs 经典扩散（末步平均环距离）
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
        criterion="全振幅对照独立 coin-walk 参考（误差 < 1e-9）；平均环距离为信息性指标",
        passed=amp_error < 1e-9 and tv_distance < 1e-9,
    )


# ---------------------------------------------------------------------------
# 图邻接与 Szegedy/MNRS（graph_walks.py）
# ---------------------------------------------------------------------------


def _cycle_table(n):
    """偶环的交替边染色邻居表（满足对合性 N(N(v,j),j)=v）。"""
    return [
        [(v + 1) % n if v % 2 == 0 else (v - 1) % n, (v - 1) % n if v % 2 == 0 else (v + 1) % n]
        for v in range(n)
    ]


def _complete_table(n):
    """带自环补齐的完全图：N(v,j)=j，D=N。"""
    return [list(range(n)) for _ in range(n)]


HYPERCUBE_Q3 = [[v ^ 1, v ^ 2, v ^ 4, v] for v in range(8)]


def _szegedy_reference(neighbors):
    """独立组装行走空间算子：返回 (维度布局, W, 初态, marked 相位对角)。"""
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
    """MNRS 迭代 (M·W)^steps 作用在初态上的独立 numpy 参考。"""
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
    """numpy 行走空间向量 → (current, peer, index) 稀疏振幅字典。"""
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
    """把 (current, peer, index) 键折成 target 整值，与驱动程序的 (target, work) 键对齐。"""
    return {
        (ref["index"](current, peer, j), 0): amplitude
        for (current, peer, j), amplitude in amplitudes.items()
    }


def verify_adjacency_superposition(report):
    """邻接 oracle 端到端：vertex/index 全叠加一次穷举全部 32 个查询。"""
    from pyqecclang.algorithms.input_model.graph_walks import gate_adjacency

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
        criterion="叠加查询逐分支 XOR 邻居表（振幅误差 < 1e-9，表外分支权重为 0）",
        passed=worst < 1e-9 and foreign < 1e-18,
    )


def verify_szegedy_walk(report):
    """Szegedy 行走步：幺正矩阵与初态上演化均对照独立反射算子组装。"""
    import numpy as np

    from pyqecclang.algorithms.input_model.graph_walks import (
        gate_adjacency,
        szegedy_setup,
        szegedy_walk,
    )

    neighbors = _cycle_table(8)
    ref = _szegedy_reference(neighbors)
    v, g = ref["v"], ref["g"]
    adjacency = gate_adjacency(neighbors)
    walk = szegedy_walk(adjacency)
    # 幺正层：to_matrix 对照 (2Π_B-I)(2Π_A-I)
    driver = Builder(
        "szegedy_unitary",
        {"current": Bits(v), "peer": Bits(v), "index": Bits(g)},
        resources_for(("walk", walk)),
    )
    invoke(driver, walk, "walk", current=driver["current"], peer=driver["peer"], index=driver["index"])
    unitary = originir_unitary(driver.finish().program())
    matrix_error = float(np.abs(unitary - ref["walk"]).max())
    # 态层：setup 后 1/3 步演化，四路径对照
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
        criterion="行走幺正与演化态对照独立反射算子组装（误差 < 1e-9）",
        passed=matrix_error < 1e-9 and state_error < 1e-9,
    )


def _hitting_reference(transition, marked):
    """独立求解 (I - P_free) h = 1（marked 顶点为 0）。"""
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
    from pyqecclang.algorithms.input_model.graph_walks import (
        gate_adjacency,
        quantum_walk_search,
        szegedy_setup,
        szegedy_walk,
    )

    n = len(neighbors)
    v = max(1, (n - 1).bit_length())
    from pyqecclang.algorithms.input_model.oracles import phase_marks

    # 独立经典参考：首达时间决定 MNRS 步数 ceil(pi/4 * sqrt(H_avg))
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
            "memory": "QRAM 邻居表" if memory else None,
        },
        metrics={
            "marked_probability": expected_prob,
            "max_state_error": state_error,
            "marked_prob_error": prob_error,
            "classical_avg_hitting_time": average_h,
            "steps_over_sqrt_h": steps / math.sqrt(average_h),
        },
        criterion="搜索末态全振幅对照独立 MNRS 参考（误差 < 1e-9）；步数/√H 与 π/4 为信息性指标",
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
    """QRAM 邻接的端到端搜索；OriginIR-ext 不含 QRAM 资源，仅走参考与 PySparQ。"""
    from pyqecclang.algorithms.input_model.graph_walks import qram_adjacency

    neighbors = _complete_table(4)
    table = {v | (j << 2): neighbors[v][j] for v in range(4) for j in range(4)}
    memory = {"setup__adj__table": table, "walk__adj__table": table}
    from pyqecclang.algorithms.input_model.graph_walks import (
        quantum_walk_search,
        szegedy_setup,
        szegedy_walk,
    )
    from pyqecclang.algorithms.input_model.oracles import phase_marks

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
            "originir_excluded": "OriginIR-ext 线路不含 QRAM 资源，无法执行 QRAM 程序",
        },
        metrics={
            "marked_probability": marked_probability,
            "max_state_error": state_error,
        },
        criterion="QRAM 绑定后搜索末态对照独立 MNRS 参考（误差 < 1e-9）",
        passed=state_error < 1e-9,
    )


def verify_markov_helpers(report):
    """transition_matrix / hitting_times / suggest_steps 对照独立 Markov 链参考。"""
    from pyqecclang.algorithms.input_model.graph_walks import (
        hitting_times,
        suggest_steps,
        transition_matrix,
    )

    tables = {
        "complete_k4": _complete_table(4),
        "hypercube_q3": HYPERCUBE_Q3,
        "cycle8": _cycle_table(8),
        # 非正则图：路径 0-1-2-3 以自环补齐到 D=2
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
    # 环上闭式解 h(v) = d(N-d)：独立闭式对照
    cycle_hits = hitting_times(transition_matrix(_cycle_table(8)), {0})
    closed_form_error = max(
        abs(cycle_hits[k] - min(k, 8 - k) * (8 - min(k, 8 - k))) for k in range(8)
    )
    # 完全图 Grover 极限：H_avg = N/|M|，步数恢复 ceil(pi/4 sqrt(N/|M|))
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
        criterion="转移矩阵/首达时间对照独立 numpy 解（误差 < 1e-9），步数公式逐例一致",
        passed=transition_error < 1e-9
        and hitting_error < 1e-9
        and closed_form_error < 1e-9
        and suggest_mismatch == 0,
    )


# ---------------------------------------------------------------------------
# 同组文档页涉及的相邻入口（transforms/qsvt/estimation/qlss）
# ---------------------------------------------------------------------------


def _householder_be(theta):
    """零信号块为 cos θ 的厄米（Householder 型）块编码：U = Ry(2θ)·Z ⊗ I_target。"""
    from pyqecclang.algorithms.input_model.operators import BlockEncoding
    from pyqecclang.algorithms.input_model.oracles import annotate

    b = Builder("householder_be", {"target": Bits(1), "signal": Bits(1)})
    b.z(b["signal"][0])
    b.ry(b["signal"][0], 2 * theta)
    return BlockEncoding(annotate(b.finish(), "block_encoding", be_alpha=1.0))


def verify_qubitization_walk(report):
    """厄米 BE 上行走 k 次的零信号概率对照 Chebyshev T_k(x)²，并核对幺正分解与谱。"""
    import numpy as np

    from pyqecclang.algorithms.common.transforms import qubitization_walk

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
    # 幺正关系 W == (2Π-I)U 与谱：两个二维子空间各贡献 e^{±i·π/6}
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
        criterion="零信号概率对照 T_k(x)²，W=(2Π-I)U 与转角 ±arccos(x)（误差 < 1e-9）",
        passed=prob_error < 1e-9 and relation_error < 1e-9 and eigenphase_error < 1e-9,
    )


def _fixed_point_case(report, name, x, delta, degree):
    from pyqecclang.algorithms.common.qsvt import fixed_point_search
    from pyqecclang.algorithms.input_model.oracles import diagonal_block_encoding, gate_database

    # 标量 BE：对角块 cos(angle/2) ≡ x（两基态相同）
    database = gate_database(1, 1, {0: 1, 1: 1})
    be = diagonal_block_encoding(database, angle_scale=2 * math.acos(x))
    program = fixed_point_search(be, delta, degree).operation.program()
    # 独立闭式：P_S(x) = 1 - δ² T_L²(c√(1-x²))，c = T_{1/L}(1/δ)
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
        criterion="零信号成功概率对照 YLC 闭式（误差 < 1e-9）；阈值外应有 P_S ≥ 1-δ²",
        passed=error < 1e-9 and (x < threshold or theory >= 1.0 - delta * delta),
    )


def verify_fixed_point_search(report):
    _fixed_point_case(report, "fixed-point-search-above-threshold", math.cos(math.pi / 6), 0.3, 5)
    _fixed_point_case(report, "fixed-point-search-below-threshold", 0.3, 0.3, 5)


def _qae_reference_distribution(a, precision):
    """QAE 相位读出的独立参考：本征相位 ±θ/π 上的 Dirichlet 核叠加。"""
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
    """n=3、3/8 标记的量子计数：相位分布对照 Dirichlet 核，读出 t̂。"""
    from pyqecclang.algorithms.common.estimation import amplitude_estimation
    from pyqecclang.algorithms.input_model.oracles import uniform_state

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
    # QAE 保证：落在真值一个栅格步内的概率 ≥ 8/π²
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
        criterion="相位分布 TVD < 1e-9；|t̂ - t| ≤ 1 且中心栅格质量 ≥ 8/π²",
        passed=tv_distance < 1e-9 and abs(estimate - t) <= 1.0 and central_mass >= 8 / math.pi**2,
    )


def verify_costa_walk_unitarity(report):
    """costa_walk 组装算子的幺正性与后端一致性（kernel 物理通道上游标注为未验证原型）。"""
    import numpy as np

    from pyqecclang.algorithms.input_model.oracles import (
        basis_state,
        diagonal_block_encoding,
        gate_database,
    )
    from pyqecclang.algorithms.qlss.qlss import costa_walk

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
        criterion="W†W - I 最大偏差 < 1e-9 且四路径零输入态一致",
        passed=deviation < 1e-9 and state_error < 1e-9,
    )


def run():
    report = Report(
        "search_walks",
        "Grover/振幅放大成功率曲线、coined 环行走、邻接 oracle 与 Szegedy/MNRS "
        "端到端行走数值，并覆盖同组文档页涉及的 qubitization/定点搜索/量子计数/Costa 行走。",
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
