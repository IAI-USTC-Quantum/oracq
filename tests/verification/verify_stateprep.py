"""态制备与数据加载组的论文级数值验证。

覆盖内容（按 docs/manual/algorithms 页面对齐）：
- 态制备 oracle（oracles.gate_state_prep / qram_state_prep）：制备振幅 vs 目标向量
  （max_error / fidelity），OriginIR-ext 小宽度与 rir-pysparq 中宽度；QRAM 角度树
  把"实现误差（对照量化角表的经典树展开）"与"方法误差（量化树 vs 精确向量）"分开。
- state_preparation.py 组合子：extend_initial / select_subspace / apply_be_to_state
  的全振幅字典级 oracle（解析 Pauli 作用与指标置换由 numpy 独立实现）。
- data_loading.py：qrom_lookup / select_swap_qrom 的真值表穷举（地址叠加一次运行
  读出全表；地址+数据叠加穷举 XOR 语义），多分区 λ 权衡；qrom_cost 与论文闭式公式
  对拍；qram_database（QRAM 资源绑定）真值表穷举。
- prepare_select.alias_prepare：selector 边缘分布与目标分布 / 量化分布的 TVD。
- density.gate_purification / maximally_mixed_purification / from_state_preparation：
  偏迹（numpy 独立实现）与目标密度矩阵的迹距离。

已知后端问题（绕行说明）：PySparQ RIR 解释器（pysparq 0.1.2.dev16）对 add_const 作用
于寄存器切片 reinterpret 视图（RIR 文本中 Span 操作数）执行结果错误——最小复现见
backend-rir-sliced-add-const 案例；oracq 的序列化、参考执行器、PySparQ 事件适配器
与 OriginIR-ext 均给出正确结果。因此 qram_state_prep 系列案例的正确性判据建立在
reference / adapter-pysparq / originir-ext 三条相互独立的路径上，rir-pysparq 的偏差
作为信息性指标记录（metrics.rir_deviation），不参与判据。

运行：PYTHONPATH=src <含 pysparq+uniqc 的 python> tests/verification/verify_stateprep.py
"""

from __future__ import annotations

import math
import random

from harness import (
    Report,
    adapter_pysparq,
    amplitude_error,
    amplitudes_to_statevector,
    basis_program,
    fidelity,
    originir_ext,
    originir_unitary,
    reference,
    rir_pysparq,
    statevector_error,
    superposition_program,
    tvd,
)

from oracq import Bits, Builder
from oracq.algorithms.common.prepare_select import alias_prepare, alias_table
from oracq.algorithms.common.state_preparation import (
    apply_be_to_state,
    extend_initial,
    select_subspace,
)
from oracq.algorithms.input_model.block_encoding import pauli_word
from oracq.algorithms.input_model.data_loading import qrom_cost, qrom_lookup, select_swap_qrom
from oracq.algorithms.input_model.density import (
    PurificationAccess,
    gate_purification,
    maximally_mixed_purification,
)
from oracq.algorithms.input_model.oracles import (
    diagonal_block_encoding,
    gate_database,
    gate_state_prep,
    qram_database,
    qram_state_angles,
    qram_state_prep,
)
from oracq.infrastructure.layout import workspace_table

ORIGINIR_QUBIT_BUDGET = 24
TABLE16 = (3, 0, 5, 2, 7, 1, 6, 4, 0, 2, 1, 7, 5, 3, 6, 4)


# ---------------------------------------------------------------------------
# 独立经典 oracle（不依赖被测实现的辅助函数）
# ---------------------------------------------------------------------------


def _complex_vector(seed, dim, nonzero=None):
    """确定性伪随机归一化复向量；nonzero 给定非零分量个数（稀疏态）。"""
    rng = random.Random(seed)
    indices = list(range(dim)) if nonzero is None else sorted(rng.sample(range(dim), nonzero))
    values = {i: complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) for i in indices}
    norm = math.sqrt(sum(abs(v) ** 2 for v in values.values()))
    return [values.get(i, 0j) / norm for i in range(dim)]


def _prep_expected(vector):
    """StatePreparation（work 宽 0）从零态出发的期望幅度字典。"""
    return {(i, 0): amplitude for i, amplitude in enumerate(vector) if amplitude}


def _dict_fidelity(actual, expected):
    keys = set(actual) | set(expected)
    overlap = sum(actual.get(k, 0j).conjugate() * expected.get(k, 0j) for k in keys)
    return abs(overlap) ** 2


def _register_widths(program):
    return [r.type.width for r in program.main.registers]


def _originir_within_budget(program):
    """预判 OriginIR-ext 全振幅态向量是否落在 24 qubit 预算内。"""
    width = sum(_register_widths(program))
    return width + workspace_table(program)[program.entry] <= ORIGINIR_QUBIT_BUDGET


def _qram_driver(operation, hadamard=(), name="verify_qram"):
    """携带 QRAM 资源的驱动程序：资源按同名映射，memory 键即资源名。"""
    module = operation.module
    b = Builder(
        name,
        {r.name: r.type for r in module.registers},
        {r.name: r.type for r in module.resources},
    )
    for key in hadamard:
        b.h(b[key])
    b.call(
        operation,
        resources={r.name: r.name for r in module.resources},
        **{r.name: b[r.name] for r in module.registers},
    )
    return b.finish().program()


def _tree_state(probs, angle_width=None):
    """多路复用旋转树的独立经典展开：逐层 cos/sin 分裂。

    与文档公式 θ = 2·atan2(√w_right, √w_left) 一致；angle_width 给定时角度按
    q = round(θ·2^aw/2π) 量化（QRAM 角表语义），None 表示不量化（gate 版语义）。
    """
    dim = len(probs)
    result = [0.0] * dim

    def visit(lo, hi, amplitude):
        if hi - lo == 1:
            result[lo] = amplitude
            return
        mid = (lo + hi) // 2
        left = sum(probs[lo:mid])
        right = sum(probs[mid:hi])
        theta = 2.0 * math.atan2(math.sqrt(right), math.sqrt(left))
        if angle_width is not None:
            quantized = round(theta * (1 << angle_width) / (2 * math.pi)) % (1 << angle_width)
            theta = quantized * (2 * math.pi) / (1 << angle_width)
        visit(lo, mid, amplitude * math.cos(theta / 2))
        visit(mid, hi, amplitude * math.sin(theta / 2))

    visit(0, dim, 1.0)
    return result


def _alias_distribution(table):
    """alias 采样分布的独立闭式：q_i ∝ quantized_i + Σ_{j: alt_j=i} (2^p − quantized_j)。"""
    scale = 1 << table.precision
    size = len(table.keep)
    counts = [0] * size
    for j in range(size):
        counts[j] += table.quantized[j]
        counts[table.alt[j]] += scale - table.quantized[j]
    return [c / (size * scale) for c in counts]


def _marginal(amplitudes, index):
    """幅度字典第 index 个寄存器的边缘分布。"""
    result = {}
    for key, amplitude in amplitudes.items():
        result[key[index]] = result.get(key[index], 0.0) + abs(amplitude) ** 2
    return result


def _originir_padded(program, memory=None):
    """OriginIR 态向量按全部入口寄存器宽度补零到完整长度。

    UniQC 的 simulate_statevector 只覆盖被使用的低位 qubit（高位空闲 qubit 隐含
    |0> 且不返回值）；补零后即可与完整稠密期望向量逐点比较。
    """
    total = 1 << sum(_register_widths(program))
    origin = list(originir_ext(program, memory))
    if len(origin) > total or len(origin) & (len(origin) - 1):
        raise AssertionError(f"OriginIR 态向量长度异常：{len(origin)}（预算 {total}）")
    return origin + [0j] * (total - len(origin))


def _reduced_density(amplitudes, system_width, environment_width):
    """numpy 独立偏迹：基态下标 system | (environment << system_width)。"""
    import numpy as np

    dim_s, dim_e = 1 << system_width, 1 << environment_width
    view = np.zeros((dim_e, dim_s), dtype=complex)
    for key, amplitude in amplitudes.items():
        system, environment = key[0], key[1]
        view[environment, system] += amplitude
    return view.T @ view.conj()


def _trace_distance(rho, sigma):
    """T(ρ,σ) = ‖ρ−σ‖₁/2，经 Hermitian 特征值（numpy）计算。"""
    import numpy as np

    values = np.linalg.eigvalsh(np.asarray(rho) - np.asarray(sigma))
    return float(0.5 * np.abs(values).sum())


def _matrix_error(rho, sigma):
    import numpy as np

    return float(np.abs(np.asarray(rho) - np.asarray(sigma)).max())


def _pauli_action(word, vector):
    """按仓库小端约定（第 i 个字母作用于 target[i]）的 Pauli 字稠密作用。"""
    import numpy as np

    gates = {
        "I": np.eye(2),
        "X": np.array([[0, 1], [1, 0]]),
        "Y": np.array([[0, -1j], [1j, 0]]),
        "Z": np.array([[1, 0], [0, -1]]),
    }
    operator = np.array([[1.0 + 0j]])
    for letter in reversed(word):  # 高位字母在 kron 左侧
        operator = np.kron(operator, gates[letter])
    return operator @ np.asarray(vector, dtype=complex)


# ---------------------------------------------------------------------------
# 态制备 oracle：gate_state_prep
# ---------------------------------------------------------------------------


def _check_state_case(report, name, vector, paths_extra_note=None):
    """单振幅向量的四路径对拍：max_error / fidelity vs 目标向量。"""
    prep = gate_state_prep(vector)
    program = basis_program(prep.operation, {})
    expected = _prep_expected(vector)
    actual = {"reference": reference(program)}
    actual["rir-pysparq"] = rir_pysparq(program, max_states=1 << (len(vector).bit_length() + 6))
    actual["adapter-pysparq"] = adapter_pysparq(program)
    paths = ["reference", "rir-pysparq", "adapter-pysparq"]
    max_error = max(amplitude_error(actual[path], expected) for path in paths)
    fid = min(_dict_fidelity(actual[path], expected) for path in paths)
    if _originir_within_budget(program):
        vector_actual = _originir_padded(program)
        max_error = max(max_error, statevector_error(vector_actual, vector))
        fid = min(fid, fidelity(vector_actual, vector))
        paths.append("originir-ext")
    parameters = {"width": prep.width, "dimension": len(vector), "nonzero": len(expected)}
    if paths_extra_note:
        parameters["note"] = paths_extra_note
    report.case(
        name,
        paths=paths,
        parameters=parameters,
        metrics={"max_error": max_error, "fidelity": fid},
        criterion="全振幅与目标向量一致（max_error < 1e-9 且 fidelity > 1 - 1e-12）",
        passed=max_error < 1e-9 and fid > 1 - 1e-12,
    )


def verify_gate_state_prep(report):
    for width in range(1, 5):
        _check_state_case(
            report,
            f"gate-state-prep-complex-w{width}",
            _complex_vector(20260900 + width, 1 << width),
        )
    _check_state_case(
        report, "gate-state-prep-dense-w6", _complex_vector(20260906, 64)
    )
    _check_state_case(
        report,
        "gate-state-prep-sparse-w8",
        _complex_vector(20260908, 256, nonzero=5),
        paths_extra_note="稀疏态（5 个非零振幅）；宽寄存器走 rir-pysparq 稀疏路径",
    )


# ---------------------------------------------------------------------------
# state_preparation.py 组合子
# ---------------------------------------------------------------------------


def verify_extend_initial(report):
    vector = _complex_vector(20260911, 8)
    prep = gate_state_prep(vector)
    extended = extend_initial(prep, 2)
    program = basis_program(extended.operation, {})
    expected = {(i, 0): amplitude for i, amplitude in enumerate(vector)}
    ref = reference(program)
    deviation = amplitude_error(ref, expected)
    deviation = max(deviation, amplitude_error(rir_pysparq(program), expected))
    deviation = max(deviation, amplitude_error(adapter_pysparq(program), expected))
    paths = ["reference", "rir-pysparq", "adapter-pysparq"]
    if _originir_within_budget(program):
        origin = _originir_padded(program)
        dense = [0j] * 32
        for i, amplitude in enumerate(vector):
            dense[i] = amplitude
        deviation = max(deviation, statevector_error(origin, dense))
        paths.append("originir-ext")
    report.case(
        "extend-initial-w3e2",
        paths=paths,
        parameters={"prep_width": 3, "extra_width": 2},
        metrics={"max_error": deviation},
        criterion="高 2 位保持 |0>，低 3 位振幅与目标向量一致（max_error < 1e-9）",
        passed=deviation < 1e-9,
    )


def verify_apply_be_to_state(report):
    """对角块编码作用于制备态：signal==0 分支应等于 cos(π·T/4)⊙ψ。"""
    vector = _complex_vector(20260912, 4)
    table = {0: 1, 1: 2, 2: 3, 3: 0}
    prep = gate_state_prep(vector)
    block = diagonal_block_encoding(gate_database(2, 2, table))
    applied = apply_be_to_state(block, prep)
    program = basis_program(applied.operation, {})
    # 独立 oracle：D·ψ，D[x,x] = cos(π·T[x]/4)（alpha = 1）
    diagonal = [math.cos(math.pi * table[x] / 4) for x in range(4)]
    expected_block = [diagonal[x] * vector[x] for x in range(4)]
    expected_success = sum(abs(v) ** 2 for v in expected_block)
    ref = reference(program)
    block_error = 0.0
    success = 0.0
    for (target, signal), amplitude in ref.items():
        if signal == 0:
            block_error = max(block_error, abs(amplitude - expected_block[target]))
            success += abs(amplitude) ** 2
    # 全态跨路径一致性
    rir = rir_pysparq(program)
    adapter = adapter_pysparq(program)
    pairwise = max(amplitude_error(ref, rir), amplitude_error(ref, adapter))
    paths = ["reference", "rir-pysparq", "adapter-pysparq"]
    if _originir_within_budget(program):
        origin = _originir_padded(program)
        pairwise = max(
            pairwise,
            statevector_error(origin, amplitudes_to_statevector(ref, _register_widths(program))),
        )
        paths.append("originir-ext")
    report.case(
        "apply-be-diagonal-w2",
        paths=paths,
        parameters={"table": table, "alpha": 1.0},
        metrics={
            "block_max_error": block_error,
            "success_probability": success,
            "success_expected": expected_success,
            "pairwise_deviation": pairwise,
        },
        criterion=(
            "signal==0 块等于对角 oracle（< 1e-9）、成功概率一致（< 1e-12）、"
            "各路径全态两两偏差 < 1e-9"
        ),
        passed=(
            block_error < 1e-9
            and abs(success - expected_success) < 1e-12
            and pairwise < 1e-9
        ),
    )


def _verify_select_subspace(report, name, width, output_width, high_value, word, seed):
    """Pauli 字态 oracle + 子空间选择：全幅度字典 oracle 与 signal==0 子块。"""
    vector = _complex_vector(seed, 1 << width)
    prep = gate_state_prep(vector)
    state = apply_be_to_state(pauli_word(word), prep)
    selected = select_subspace(state, output_width, high_value)
    program = basis_program(selected.operation, {})
    # 独立 oracle：ψ' = P·ψ；终态振幅 (t, h) ↦ ψ'[t | ((h ^ high_value) << k)]
    acted = _pauli_action(word, vector)
    extra = width - output_width
    expected = {}
    for h in range(1 << extra):
        for t in range(1 << output_width):
            expected[(t, h)] = acted[t | ((h ^ high_value) << output_width)]
    ref = reference(program)
    max_error = amplitude_error(ref, expected)
    max_error = max(max_error, amplitude_error(rir_pysparq(program), expected))
    paths = ["reference", "rir-pysparq"]
    if _originir_within_budget(program):
        origin = _originir_padded(program)
        dense = amplitudes_to_statevector(expected, _register_widths(program))
        max_error = max(max_error, statevector_error(origin, dense))
        paths.append("originir-ext")
    # signal==0 分支即"高位 == high_value"的后选子向量
    selected_block_error = max(
        (
            abs(amplitude - acted[t | (high_value << output_width)])
            for (t, h), amplitude in ref.items()
            if h == 0
        ),
        default=0.0,
    )
    report.case(
        name,
        paths=paths,
        parameters={
            "width": width,
            "output_width": output_width,
            "high_value": high_value,
            "pauli_word": word,
        },
        metrics={"max_error": float(max_error), "selected_block_error": float(selected_block_error)},
        criterion="全幅度字典与 Pauli oracle 一致，且 signal==0 块为后选子向量（均 < 1e-9）",
        passed=max_error < 1e-9 and selected_block_error < 1e-9,
    )


def verify_select_subspace(report):
    _verify_select_subspace(report, "select-subspace-w3-k2-hv1", 3, 2, 1, "YXI", 20260913)
    _verify_select_subspace(report, "select-subspace-w3-k1-hv3", 3, 1, 3, "XZY", 20260914)


# ---------------------------------------------------------------------------
# QRAM 角度树态制备（含已知后端偏差的记录）
# ---------------------------------------------------------------------------


def _verify_qram_state_prep(report, name, width, angle_width, seed):
    vector = [abs(v) for v in _complex_vector(seed, 1 << width)]  # 非负实幅度
    norm = math.sqrt(sum(v * v for v in vector))
    vector = [v / norm for v in vector]
    prep = qram_state_prep(width, angle_width)
    memory = {"angles": qram_state_angles(vector, angle_width)}
    program = prep.operation.program()
    # 独立 oracle：量化角表的经典树展开
    quantized = _tree_state([v * v for v in vector], angle_width)
    expected = _prep_expected(quantized)
    exact = _prep_expected(vector)
    ref = reference(program, memory)
    adapter = adapter_pysparq(program, memory)
    origin = _originir_padded(program, memory)
    widths = _register_widths(program)
    dense_quantized = amplitudes_to_statevector(expected, widths)
    dense_exact = amplitudes_to_statevector(exact, widths)
    impl_error = max(amplitude_error(ref, expected), amplitude_error(adapter, expected))
    impl_error = max(impl_error, statevector_error(origin, dense_quantized))
    method_error = statevector_error(quantized, vector)
    fid = fidelity(origin, dense_exact)
    rir = rir_pysparq(program, memory)
    rir_deviation = amplitude_error(rir, expected)
    # 自洽检查：同一 oracle 不量化时应精确复现 gate 版线路
    gate = gate_state_prep(vector)
    gate_ref = reference(basis_program(gate.operation, {}))
    fine = _prep_expected(_tree_state([v * v for v in vector]))
    oracle_self_check = amplitude_error(gate_ref, fine)
    report.case(
        name,
        paths=["reference", "adapter-pysparq", "originir-ext"],
        parameters={
            "width": width,
            "angle_width": angle_width,
            "angle_resolution": 2 * math.pi / (1 << angle_width),
        },
        metrics={
            "impl_error": impl_error,
            "method_error": method_error,
            "fidelity_vs_exact": fid,
            "oracle_self_check": oracle_self_check,
            "rir_deviation": rir_deviation,
        },
        criterion=(
            "实现误差（对照量化角树 oracle）< 1e-9 且三路径一致；"
            "rir-pysparq 因 Span add_const 后端问题仅记录 rir_deviation"
        ),
        passed=impl_error < 1e-9 and oracle_self_check < 1e-9,
    )


def verify_qram_state_prep(report):
    _verify_qram_state_prep(report, "qram-state-prep-w2-a8", 2, 8, 20260915)
    _verify_qram_state_prep(report, "qram-state-prep-w3-a10", 3, 10, 20260916)


def verify_backend_rir_sliced_add_const(report):
    """已知后端问题的最小复现：add_const 作用于切片 reinterpret 视图。

    RIR 文本中操作数为 Span(register=work, start=0, width=2)、值 1 与 3；
    参考执行器、PySparQ 事件适配器均正确复净到 |0>，而 PySparQ RIR 解释器把
    常量加在错误的位偏移上（work = 4）。该问题影响 qram_state_prep 的地址
    簿记，故 qram-state-prep 案例的判据不含 rir-pysparq 路径。
    """
    b = Builder("stateprep_diag_sliced_add", {"work": Bits(4)})
    b.add_const(b["work"][:2].reinterpret("uint"), 1)
    b.add_const(b["work"][:2].reinterpret("uint"), 3)
    program = b.finish().program()
    ref = reference(program)
    adapter = adapter_pysparq(program)
    rir = rir_pysparq(program)
    clean = set(ref) == {(0,)} and set(adapter) == {(0,)}
    deviation = amplitude_error(rir, ref)
    report.case(
        "backend-rir-sliced-add-const",
        paths=["reference", "adapter-pysparq", "rir-pysparq"],
        parameters={"pattern": "add_const on Span(work, 0, 2) reinterpreted as uint"},
        metrics={"reference_and_adapter_clean": clean, "rir_deviation": deviation},
        criterion="参考与适配器路径复净（判据）；rir-pysparq 的 Span 偏差仅记录",
        passed=clean,
    )


# ---------------------------------------------------------------------------
# QROM 数据加载
# ---------------------------------------------------------------------------


def _truth_table_expected(table, address_bits, data_bits, data_superposition):
    """叠加查询的期望幅度字典：{(a, d0 ^ T[a])} 均匀分布。"""
    n_a, n_d = 1 << address_bits, 1 << data_bits
    uniform = 1.0 / math.sqrt(n_a * (n_d if data_superposition else 1))
    expected = {}
    for address in range(n_a):
        word = table[address] if address < len(table) else 0
        if data_superposition:
            for data in range(n_d):
                expected[(address, data ^ word)] = uniform
        else:
            expected[(address, word)] = uniform
    return expected


def _check_database_case(report, name, database, table, data_superposition, note=None):
    """单数据库程序的叠加穷举：期望字典 vs 各后端路径。"""
    registers = ["address", "data"] if data_superposition else ["address"]
    program = superposition_program(database.operation, registers)
    address_bits = database.address_width
    data_bits = database.data_width
    expected = _truth_table_expected(table, address_bits, data_bits, data_superposition)
    ref = reference(program)
    max_error = amplitude_error(ref, expected)
    max_error = max(max_error, amplitude_error(rir_pysparq(program), expected))
    max_error = max(max_error, amplitude_error(adapter_pysparq(program), expected))
    paths = ["reference", "rir-pysparq", "adapter-pysparq"]
    parameters = {
        "addresses": 1 << address_bits,
        "data_bits": data_bits,
        "data_superposition": data_superposition,
        "branches": len(expected),
    }
    if _originir_within_budget(program):
        origin = _originir_padded(program)
        dense = amplitudes_to_statevector(expected, _register_widths(program))
        max_error = max(max_error, statevector_error(origin, dense))
        paths.append("originir-ext")
    else:
        total = sum(_register_widths(program)) + workspace_table(program)[program.entry]
        parameters["originir_note"] = f"workspace 共 {total} qubit，超 24 预算，仅走 pysparq 路径"
    if note:
        parameters["note"] = note
    report.case(
        name,
        paths=paths,
        parameters=parameters,
        metrics={"max_error": max_error},
        criterion="叠加穷举的查询结果与真值表逐振幅一致（max_error < 1e-9）",
        passed=max_error < 1e-9,
    )


def verify_qrom_lookup(report):
    database = qrom_lookup(TABLE16)
    _check_database_case(report, "qrom-lookup-truth-table-n16", database, TABLE16, False)
    sparse = {0: 5, 3: 2, 6: 7}
    database = qrom_lookup(sparse)
    padded = tuple(sparse.get(i, 0) for i in range(8))
    _check_database_case(
        report,
        "qrom-lookup-sparse-dict",
        database,
        padded,
        False,
        note="稀疏字典表，缺失地址按 0 处理",
    )
    wide = tuple((7 * i + 3) % 16 for i in range(64))
    database = qrom_lookup(wide, data_bits=4)
    _check_database_case(report, "qrom-lookup-wide-n64", database, wide, False)


def verify_select_swap(report):
    for partitions in (1, 2, 4, 8, 16):
        database = select_swap_qrom(TABLE16, partitions=partitions)
        _check_database_case(
            report,
            f"select-swap-truth-table-l{partitions}",
            database,
            TABLE16,
            False,
            note=f"λ = {partitions} 个窗口分区",
        )
    database = select_swap_qrom(TABLE16, partitions=4)
    _check_database_case(
        report,
        "select-swap-xor-superposition-l4",
        database,
        TABLE16,
        True,
        note="地址+数据全叠加，穷举 XOR 语义 d ↦ d ⊕ T[a]",
    )
    wide = tuple((7 * i + 3) % 16 for i in range(64))
    for partitions in (2, 8):
        database = select_swap_qrom(wide, partitions=partitions, data_bits=4)
        _check_database_case(
            report,
            f"select-swap-wide-n64-l{partitions}",
            database,
            wide,
            False,
        )


def verify_qrom_cost(report):
    """qrom_cost 与论文闭式 4(N/λ − 1 + b(λ − 1)) 及生成属性的双向对拍。"""
    failures = []
    grid = [(16, 3, p) for p in (1, 2, 4, 8, 16)]
    grid += [(64, 4, p) for p in (1, 2, 4, 8)]
    grid.append((1024, 32, 8))
    for n_addresses, data_bits, partitions in grid:
        cost = qrom_cost(n_addresses, data_bits, partitions)
        address_bits = max(1, (n_addresses - 1).bit_length())
        low_bits = partitions.bit_length() - 1
        high_bits = address_bits - low_bits
        checks = {
            "select_toffoli": (cost.select_toffoli, (1 << high_bits) - 1),
            "swap_toffoli": (cost.swap_toffoli, data_bits * (partitions - 1)),
            "t_count": (
                cost.t_count,
                4 * ((n_addresses >> low_bits) - 1 + data_bits * (partitions - 1)),
            ),
            "t_depth": (cost.t_depth, (1 << high_bits) - 1 + max(0, partitions - 1)),
            "round_trip": (cost.round_trip_t_count, 2 * cost.t_count),
            "work_qubits": (cost.work_qubits, partitions * data_bits),
            "fanout_qubits": (cost.fanout_qubits, (partitions - 1) * low_bits),
        }
        for key, (actual, expected) in checks.items():
            if actual != expected:
                failures.append((n_addresses, data_bits, partitions, key, actual, expected))
    # 生成操作的模块属性与 qrom_cost 输出一致
    for partitions in (1, 2, 4, 8, 16):
        database = select_swap_qrom(TABLE16, partitions=partitions)
        attrs = dict(database.operation.module.attributes)
        cost = qrom_cost(16, 3, partitions)
        if attrs["t_count"] != cost.t_count or attrs["work_qubits"] != cost.work_qubits:
            failures.append((16, 3, partitions, "attributes"))
    baseline = qrom_lookup(TABLE16)
    attrs = dict(baseline.operation.module.attributes)
    if attrs["t_count"] != qrom_cost(16, 3, 1).t_count:
        failures.append((16, 3, 1, "qrom_lookup.attributes"))
    report.case(
        "qrom-cost-formula",
        paths=["classical"],
        parameters={"grid": [list(item) for item in grid]},
        metrics={"failures": len(failures)},
        criterion="成本模型与闭式公式及生成属性完全一致（failures == 0）",
        passed=not failures,
    )


def verify_qram_database(report):
    """QRAM 资源绑定的 XOR 数据库：真值表穷举（memory 提供表数据）。"""
    database = qram_database(2, 3)
    table = (5, 2, 7, 1)
    memory = {"table": dict(enumerate(table))}
    for data_superposition in (False, True):
        registers = ["address", "data"] if data_superposition else ["address"]
        program = _qram_driver(database.operation, registers)
        expected = _truth_table_expected(table, 2, 3, data_superposition)
        ref = reference(program, memory)
        max_error = amplitude_error(ref, expected)
        max_error = max(max_error, amplitude_error(rir_pysparq(program, memory), expected))
        max_error = max(max_error, amplitude_error(adapter_pysparq(program, memory), expected))
        paths = ["reference", "rir-pysparq", "adapter-pysparq"]
        origin = _originir_padded(program, memory)
        dense = amplitudes_to_statevector(expected, _register_widths(program))
        max_error = max(max_error, statevector_error(origin, dense))
        paths.append("originir-ext")
        report.case(
            "qram-database-xor-superposition" if data_superposition else "qram-database-truth-table",
            paths=paths,
            parameters={
                "address_bits": 2,
                "data_bits": 3,
                "data_superposition": data_superposition,
                "branches": len(expected),
            },
            metrics={"max_error": max_error},
            criterion="QRAM 查询与真值表逐振幅一致（max_error < 1e-9）",
            passed=max_error < 1e-9,
        )


# ---------------------------------------------------------------------------
# alias 采样制备
# ---------------------------------------------------------------------------


def _check_alias_case(report, name, coefficients, precision, use_qram):
    values = [complex(c) for c in coefficients]
    alpha = sum(abs(c) for c in values)
    width = (len(values) - 1).bit_length()
    size = 1 << width
    exact = {i: abs(c) / alpha for i, c in enumerate(values)}
    table = alias_table(coefficients, precision=precision)
    if use_qram:
        alias = alias_prepare(coefficients, precision=precision)
        program = alias.preparation.operation.program()
        memory = alias.memory
        binding = "qram"
    else:
        database = gate_database(width, precision + width, dict(table.table))
        alias = alias_prepare(coefficients, precision=precision, database=database)
        program = basis_program(alias.preparation.operation, {})
        memory = None
        binding = "gate_database"
    quantized = {i: p for i, p in enumerate(_alias_distribution(table))}
    # 模块自带 distribution() 与独立闭式的一致性（信息性指标）
    self_consistency = tvd(dict(enumerate(table.distribution())), quantized)
    max_states = 1 << (width + precision + 4)
    results = {}
    for path, runner in (("reference", reference), ("rir-pysparq", rir_pysparq)):
        amplitudes = runner(program, memory, max_states=max_states)
        results[path] = _marginal(amplitudes, 0)
    tvd_exact = max(tvd(results[path], exact) for path in results)
    tvd_quantized = max(tvd(results[path], quantized) for path in results)
    cross = tvd(results["reference"], results["rir-pysparq"])
    bound = size / (1 << precision)
    report.case(
        name,
        paths=["reference", "rir-pysparq"],
        parameters={
            "coefficients": len(values),
            "selector_width": width,
            "precision": precision,
            "binding": binding,
            "tvd_bound": bound,
        },
        metrics={
            "tvd_vs_exact": tvd_exact,
            "tvd_vs_quantized": tvd_quantized,
            "cross_path_tvd": cross,
            "distribution_self_consistency": self_consistency,
        },
        criterion=(
            f"边缘分布与量化 alias 分布一致（TVD < 1e-12），与精确分布 TVD ≤ 2^w·2^-p "
            f"= {bound:.6f}，两路径互差 < 1e-12"
        ),
        passed=(
            tvd_quantized < 1e-12 and tvd_exact <= bound + 1e-12 and cross < 1e-12
        ),
    )


def verify_alias_preparation(report):
    coefficients = (0.6, -0.8, 0.3j, -0.5)
    _check_alias_case(report, "alias-preparation-gate-w2-p8", coefficients, 8, False)
    _check_alias_case(report, "alias-preparation-qram-w2-p8", coefficients, 8, True)
    wide = (1.1, -0.4, 0.7j, 0.0, 2.3, -0.9)
    _check_alias_case(report, "alias-preparation-gate-w3-p10", wide, 10, False)
    _check_alias_case(report, "alias-preparation-uniform-p8", (1.0, 1.0, 1.0, 1.0), 8, False)


# ---------------------------------------------------------------------------
# 纯化访问
# ---------------------------------------------------------------------------


def _check_purification_case(report, name, operation, system_width, environment_width, target):
    program = basis_program(operation, {})
    paths = ["reference", "rir-pysparq", "adapter-pysparq"]
    distance = 0.0
    element_error = 0.0
    amplitudes = {}
    for path, runner in (
        ("reference", reference),
        ("rir-pysparq", rir_pysparq),
        ("adapter-pysparq", adapter_pysparq),
    ):
        amplitudes[path] = runner(program)
        reduced = _reduced_density(amplitudes[path], system_width, environment_width)
        distance = max(distance, _trace_distance(reduced, target))
        element_error = max(element_error, _matrix_error(reduced, target))
    if _originir_within_budget(program):
        origin = _originir_padded(program)
        as_dict = {
            ((index & ((1 << system_width) - 1)), index >> system_width): value
            for index, value in enumerate(origin)
        }
        reduced = _reduced_density(as_dict, system_width, environment_width)
        distance = max(distance, _trace_distance(reduced, target))
        element_error = max(element_error, _matrix_error(reduced, target))
        paths.append("originir-ext")
    # Schmidt 结构指标：非零幅度中 system != environment 的最大概率（信息性）
    off_diagonal = max(
        (
            abs(amplitude) ** 2
            for key, amplitude in amplitudes["reference"].items()
            if key[0] != key[1]
        ),
        default=0.0,
    )
    report.case(
        name,
        paths=paths,
        parameters={
            "system_width": system_width,
            "environment_width": environment_width,
            "dimension": 1 << system_width,
        },
        metrics={
            "trace_distance": distance,
            "max_element_error": element_error,
            "off_schmidt_probability": off_diagonal,
        },
        criterion="偏迹后的约化密度矩阵与目标一致（trace_distance < 1e-9）",
        passed=distance < 1e-9,
    )


def verify_purification(report):
    import numpy as np

    rho2 = ((0.7 + 0j, 0.1 - 0.05j), (0.1 + 0.05j, 0.3 + 0j))
    _check_purification_case(
        report, "purification-gate-rho2", gate_purification(rho2).operation, 1, 1, rho2
    )
    for dimension, seed in ((4, 7), (16, 11)):
        n = dimension.bit_length() - 1
        rng = np.random.default_rng(seed)
        gram = rng.standard_normal((dimension, dimension)) + 1j * rng.standard_normal(
            (dimension, dimension)
        )
        rho = gram.conj().T @ gram
        rho /= np.trace(rho)
        rho = tuple(tuple(complex(v) for v in row) for row in rho)
        _check_purification_case(
            report,
            f"purification-gate-rho{dimension}",
            gate_purification(rho).operation,
            n,
            n,
            rho,
        )
    width = 3
    dimension = 1 << width
    identity = tuple(
        tuple((1.0 / dimension if i == j else 0.0) + 0j for j in range(dimension))
        for i in range(dimension)
    )
    _check_purification_case(
        report,
        "purification-bell-w3",
        maximally_mixed_purification(width).operation,
        width,
        width,
        identity,
    )
    prep = gate_state_prep([math.sqrt(0.3), math.sqrt(0.7)])
    access = PurificationAccess.from_state_preparation(prep)
    pure = (
        (0.3 + 0j, math.sqrt(0.21) + 0j),
        (math.sqrt(0.21) + 0j, 0.7 + 0j),
    )
    _check_purification_case(
        report, "purification-pure-adapter", access.operation, 1, 0, pure
    )


# ---------------------------------------------------------------------------
# 块编码幺正层面的交叉验证（pysparq 自带块编码模块之外的独立路径）
# ---------------------------------------------------------------------------


def verify_gate_state_prep_unitary(report):
    """gate_state_prep 的幺正第一列等于目标向量（UniQC Circuit.to_matrix 路径）。"""
    vector = _complex_vector(20260917, 8)
    prep = gate_state_prep(vector)
    unitary = originir_unitary(basis_program(prep.operation, {}))
    column = [unitary[row, 0] for row in range(8)]
    error = statevector_error(column, vector)
    fid = fidelity(column, vector)
    report.case(
        "gate-state-prep-unitary-column-w3",
        paths=["originir-ext+to_matrix"],
        parameters={"width": 3},
        metrics={"max_error": float(error), "fidelity": float(fid)},
        criterion="线路幺正第一列与目标向量一致（max_error < 1e-12）",
        passed=error < 1e-12 and fid > 1 - 1e-15,
    )


def run():
    report = Report(
        "stateprep",
        "态制备 oracle、state_preparation 组合子、QROM/QRAM 数据加载、alias 采样"
        "与纯化访问的数值验证：振幅/真值表逐点穷举、分布 TVD、偏迹迹距离。",
    )
    verify_gate_state_prep(report)
    verify_extend_initial(report)
    verify_apply_be_to_state(report)
    verify_select_subspace(report)
    verify_qram_state_prep(report)
    verify_backend_rir_sliced_add_const(report)
    verify_qrom_lookup(report)
    verify_select_swap(report)
    verify_qrom_cost(report)
    verify_qram_database(report)
    verify_alias_preparation(report)
    verify_purification(report)
    verify_gate_state_prep_unitary(report)
    report.write()
    return report


if __name__ == "__main__":
    run()
