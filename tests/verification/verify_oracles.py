"""Oracle 目录与查询算法的论文级数值验证。

验证对象：src/oracq/algorithms/oracles.py 与 oracle_algorithms.py。

- XOR database：gate 真值表在逐基态与全叠加两种模式下穷举全部输入域，
  QRAM 绑定、abstract 声明的 gate/QRAM 双绑定一致性，四条后端路径对拍；
- Bernstein–Vazirani（论文点名算法）：3–8 bit 多规模 × 多秘密串 × 两种偏置，
  originir-ext 与 rir-pysparq 两条主路径（外加 reference 与 adapter 对拍），
  恢复串等于秘密串的概率精确为 1；含开放声明 → 绑定的端到端路线；
- Deutsch–Jozsa：常量/平衡函数族判定（含 BooleanNetwork 编译 oracle 的
  端到端路线，顺带对拍 boolean-networks 的线路语义与经典求值）；
- Simon：采样分布的相位敏感逐点对拍 + GF(2) 消元恢复周期（含宽位 CNOT
  线性 oracle 路线）。

运行：PYTHONPATH=src <含 pysparq+uniqc 的 python> tests/verification/verify_oracles.py
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
# 与库实现无关的经典参考工具
# ---------------------------------------------------------------------------


def _dot(a, b):
    """GF(2) 点积。"""
    return (a & b).bit_count() & 1


def _pseudo_table(address_width, data_width, seed):
    """确定性 xorshift64* 伪随机真值表（独立构造，不复用库代码）。"""
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
    """按最高位主元消元，确定性选取 count 个线性无关向量。"""
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
    """XOR database 驱动程序。

    harness 的 basis/superposition 助手不穿线 QRAM 资源，这里显式按名映射，
    同时支持基态初值与部分寄存器叠加的组合。
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
    """全叠加下 XOR database 的精确参考态：|a,d> → |a, d XOR table[a]> 均匀叠加。"""
    uniform = 1 / math.sqrt(1 << (address_width + data_width))
    return {
        (a, d ^ table[a]): uniform
        for a in range(1 << address_width)
        for d in range(1 << data_width)
    }


def _bv_expected(secret, bias):
    """BV/DJ 相位反冲电路的闭式末态：input 集中于 s，answer 留在 |->。"""
    sign = 1.0 if bias == 0 else -1.0
    root = 1 / math.sqrt(2)
    return {(secret, 0): sign * root, (secret, 1): -sign * root}


def _input_marginal(state_probs, input_index=0):
    marginal = {}
    for key, p in state_probs.items():
        marginal[key[input_index]] = marginal.get(key[input_index], 0.0) + p
    return marginal


def _origin_input_probability(vector, input_width, target):
    """从 OriginIR 全振幅态向量汇总 input == target 的概率（input 在低位）。"""
    mask = (1 << input_width) - 1
    return sum(abs(v) ** 2 for i, v in enumerate(vector) if (i & mask) == target)


# ---------------------------------------------------------------------------
# A 组：XOR database 真值表语义
# ---------------------------------------------------------------------------


def verify_xor_gate_basis(report):
    """基态模式：逐点穷举全部 (address, data) 初态，输出必须恰为单一基态。"""
    # (2,3) 全路径；(3,2) 走 reference + rir-pysparq 控制运行总量
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
            criterion="全部 2^(aw+dw) 个基态初态输出恰为 |a, d XOR table[a]>（failures == 0）",
            passed=failures == 0,
        )


def verify_xor_gate_superposition(report):
    """叠加模式：一次运行穷举全部输入域，四路径与经典置换逐振幅对拍。"""
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
            criterion="叠加分支与经典置换逐振幅一致（max_error < 1e-9）",
            passed=worst < 1e-9,
        )


def verify_xor_qram(report):
    """QRAM 绑定：叠加穷举 + 非零 data 初值的基态抽点。"""
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
        criterion="QRAM 叠加逐振幅一致（max_error < 1e-9）且基态抽点全部命中",
        passed=worst < 1e-9 and failures == 0,
    )


def verify_xor_binding_consistency(report):
    """abstract 声明的 gate / QRAM 双绑定一致性（验证矩阵 V4 缺口的参数化对拍）。"""
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
        criterion="同一开放声明的两种绑定与经典置换两两一致（max_error < 1e-9）",
        passed=worst < 1e-9,
    )


# ---------------------------------------------------------------------------
# B 组：Bernstein–Vazirani（论文点名，优先）
# ---------------------------------------------------------------------------


def _bv_secrets(width):
    mask = (1 << width) - 1
    alternating = sum(1 << i for i in range(0, width, 2))
    hashed = 0x9E3779B97F4A7C15 & mask
    return sorted({1, mask, alternating, hashed})


def verify_bv_recovery(report):
    """3–8 bit 多规模：全秘密串 × 偏置实例，四路径恢复概率精确为 1。"""
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
            criterion="全部实例读出串 == 秘密串（success_probability > 1 - 1e-9，tvd < 1e-9）",
            passed=worst["p_success"] > 1 - 1e-9 and worst["tvd"] < 1e-9,
        )


def verify_bv_oracle_truth_table(report):
    """affine_boolean_oracle 自身的真值表穷举：f(x) = s·x XOR c 逐分支对拍。"""
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
            criterion="oracle 叠加分支与 s·x XOR c 逐振幅一致（max_error < 1e-9）",
            passed=worst < 1e-9,
        )


def verify_bv_open_binding(report):
    """开放声明 → gate/QRAM 绑定的端到端恢复（docs 描述的开放 oracle 路线）。"""
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
            criterion="两种绑定恢复概率均为 1 且与闭式态一致（p > 1 - 1e-9，max_error < 1e-9）",
            passed=p_min > 1 - 1e-9 and deviation < 1e-9,
        )


# ---------------------------------------------------------------------------
# C 组：Deutsch–Jozsa
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
    """常量/平衡函数族：P(input==0) 须精确等于 1（常量）或 0（平衡）。"""
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
            # 常量与线性平衡（parity = s·x，s 全 1）有闭式末态，做相位敏感对拍
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
            criterion="判定全部正确且 P(input=0) 偏差 < 1e-9（闭式情形逐振幅 < 1e-9）",
            passed=decision_errors == 0 and worst_p0 < 1e-9 and worst_state < 1e-9,
        )


def _network_database(net, input_name, output_name, width):
    """把 BooleanNetwork 编译产物包装成 XOR database 接口（输出为 XOR 拷贝语义）。"""
    net_op = net.operation()
    b = Builder("net_db", {"address": Bits(width), "data": Bits(1)})
    b.call(net_op, **{input_name: b["address"], output_name: b["data"]})
    return XorDatabase(annotate(b.finish(), "database_xor", implementation="boolean_network"))


def verify_dj_boolean_network(report):
    """BooleanNetwork 编译 oracle 的 DJ 端到端 + 线路语义对经典求值穷举。"""
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
        # 线路语义穷举：叠加一次运行对照 net.evaluate 的全部 16 个分支
        sweep = _db_program(database.operation, superpose=("address", "data"), name="net_sw")
        expected = _xor_expected(3, 1, table)
        worst = amplitude_error(rir_pysparq(sweep), expected)
        vector = originir_ext(sweep)
        worst = max(
            worst, statevector_error(vector, amplitudes_to_statevector(expected, [3, 1]))
        )
        # DJ 判定
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
            criterion="网络线路与经典求值逐分支一致且 DJ 判定正确（误差 < 1e-9）",
            passed=worst < 1e-9 and p0_error < 1e-9,
        )


# ---------------------------------------------------------------------------
# D 组：Simon
# ---------------------------------------------------------------------------


def _simon_table(n, s):
    """周期恰为 s 的二对一线性函数真值表（独立构造 + 经典穷举自证）。"""
    pivot = (s & -s).bit_length() - 1

    def f(x):
        z = x ^ (((x >> pivot) & 1) * s)
        return (z & ((1 << pivot) - 1)) | ((z >> (pivot + 1)) << pivot)

    table = [f(x) for x in range(1 << n)]
    assert all(table[x] == table[x ^ s] for x in range(1 << n))
    assert len(set(table)) == 1 << (n - 1)
    return table


def _simon_cnot_oracle(n, s):
    """线性 Simon 函数的 CNOT 实现（宽位时替代真值表以降低门数）。"""
    pivot = (s & -s).bit_length() - 1
    b = Builder(f"simon_linear_{n}_{s}", {"address": Bits(n), "data": Bits(n - 1)})
    for j in range(n - 1):
        m = j if j < pivot else j + 1
        b.xor(b["address"][m], b["data"][j])
        if (s >> m) & 1:
            b.xor(b["address"][pivot], b["data"][j])
    return XorDatabase(annotate(b.finish(), "database_xor", implementation="cnot_linear"))


def _simon_expected(n, s, table):
    """Simon 采样电路的精确联合分布（含相位）：支持集 2^{2n-2} 个等幅基态。"""
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
    """采样分布相位敏感对拍 + 消元恢复周期（3–8 bit）。"""
    for n in range(3, 9):
        use_cnot = n >= 7  # 真值表门数随 2^n 增长，宽位走 CNOT 线性实现
        secrets = sorted({1, (1 << n) - 1, (0x9E3779B97F4A7C15 & ((1 << n) - 1)) | 1})
        worst_amp, worst_tvd, recoveries, attempts = 0.0, 0.0, 0, 0
        case_paths = ["rir-pysparq"] if use_cnot else ["reference", "rir-pysparq"]
        if not use_cnot and n <= 5:
            case_paths.append("originir-ext")  # 稠密小实例走 UniQC 全振幅
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
                    worst_tvd = math.inf  # 支持集越出 s 的正交子空间
                # 用该路径自身的采样支持集做端到端消元恢复
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
            criterion="联合分布逐振幅一致且每次采样消元都恢复周期（recoveries == attempts）",
            passed=worst_amp < 1e-9 and worst_tvd < 1e-9 and recoveries == attempts,
        )


def verify_simon_rank_deficiency(report):
    """样本不足时零空间保持多维且周期落在其张成内（信息性指标）。"""
    n, s = 4, 0b1011
    table = _simon_table(n, s)
    program = simon_sample(gate_database(n, n - 1, table)).program()
    state = rir_pysparq(program)
    support = sorted({y for (y, _v), p in probabilities(state).items() if p > 1e-15})
    samples = _gf2_independent(support, n - 2)  # 少一个独立样本
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
        criterion="零空间恰好多出一维且 s 在其张成内",
        passed=len(basis) == 2 and in_span,
    )


def run():
    report = Report(
        "oracles",
        "XOR database 三种绑定、BV/DJ/Simon 查询算法在四条真实后端路径上的"
        "穷举式与端到端数值验证（BV 覆盖 3–8 bit 多规模双路径）。",
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
