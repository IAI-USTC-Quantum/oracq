"""块编码组的论文级数值验证。

覆盖模块：
- ``pyqecclang.algorithms.block_encoding``（BE 组合代数）
- ``pyqecclang.algorithms.prepare_select``（PREPARE–SELECT 分解）
- ``pyqecclang.algorithms.sparse``（稀疏访问辅助与 CKS 稀疏 BE）
- ``pyqecclang.algorithms.lowrank``（DF/THC 低秩块编码）

正确性 oracle：块编码的零信号角块 == A/α。小规模用 OriginIR-ext 经 UniQC
``Circuit.to_matrix`` 取全幺正，再用 ``harness.effective_block`` 提取有效块并
报告失败分支（信号≠0）泄漏上界；对角/稀疏块编码跑多矩阵、多规模；三对角与
一般稀疏矩阵另用 pysparq 自带块编码模块（``BlockEncodingTridiagonal`` /
``BlockEncodingViaQRAM``）独立编码同一矩阵，双向对拍有效块。

后端备注：早期 pysparq.rir 解释器对切片视图上的 ``add_const`` 不按视图位宽
回卷（进位写入相邻位），曾影响 ``qram_state_prep`` 旋转树；该缺陷已于
2026-09 修复，本脚本对切片 ``add_const`` 往返与 ``qram_state_prep``
width=1..3 的最小复测全部为 0 偏差，QRAM PREPARE 案例恢复 reference 与
rir_pysparq 双路径对拍。

运行：PYTHONPATH=src <含 pysparq+uniqc 的 python> tests/verification/verify_blockencoding.py
"""

from __future__ import annotations

import math

import numpy as np
from harness import (
    Report,
    adapter_pysparq,
    amplitude_error,
    amplitudes_to_statevector,
    basis_program,
    effective_block,
    originir_ext,
    originir_unitary,
    reference,
    rir_pysparq,
    statevector_error,
)

from pyqecclang import Builder, bind, simulate, unresolved
from pyqecclang.algorithms.arithmetic import FixedFormat
from pyqecclang.algorithms.block_encoding import (
    adjoint_be,
    direct_sum,
    kronecker_sum,
    lcu,
    matrix_pauli_encoding,
    pad_signal,
    pauli_word,
    projector,
    tensor,
    truncated_shift,
)
from pyqecclang.algorithms.data_loading import select_swap_qrom
from pyqecclang.algorithms.lowrank import (
    DoubleFactorization,
    THCDecomposition,
    _diagonal_encoding,
    diagonalize_symmetric,
    double_factorized_encoding,
    thc_encoding,
)
from pyqecclang.algorithms.operators import product
from pyqecclang.algorithms.oracles import (
    SparseAccess,
    gate_database,
    qram_database,
    sparse_entry,
    sparse_location_gate,
    sparse_location_qram,
)
from pyqecclang.algorithms.prepare_select import (
    abstract_prepare,
    alias_prepare,
    gate_prepare,
    lcu_prepare_select,
    qram_prepare,
    select_pauli,
)
from pyqecclang.algorithms.sparse import (
    batch_lookup,
    chebyshev_block,
    compare_words,
    magnitude_rotation,
    prefix_state,
    real_symmetric_sparse_encoding,
    reversible_lookup,
    value_transposition,
    word_rotation,
)

EXACT = 1e-9  # 精确路径（门级、无量化）的块对拍容差

PAULI = {
    "I": np.eye(2),
    "X": np.array([[0, 1], [1, 0]]),
    "Y": np.array([[0, -1j], [1j, 0]]),
    "Z": np.array([[1, 0], [0, -1]]),
}


def pauli_matrix(word):
    """Pauli 字的稠密矩阵；word[0] 作用于最低位（与 pauli_word 的比特约定一致）。"""
    result = PAULI[word[0]]
    for letter in word[1:]:
        result = np.kron(PAULI[letter], result)
    return result


def pauli_l1(matrix):
    """独立计算 Pauli 展开的 l1 上界：c_P = Tr(P†M)/2^n（numpy 直算，不经被测实现）。"""
    n = (len(matrix) - 1).bit_length()
    total = 0.0
    for letters in __import__("itertools").product("IXYZ", repeat=n):
        word = "".join(letters)
        total += abs(np.trace(pauli_matrix(word).conj().T @ matrix) / len(matrix))
    return total


def block_via_reference(operation, alpha, width):
    """reference 路径逐列提取 (0,0) 块（乘回 alpha 前）与失败分支泄漏上界。"""
    program = operation.program() if hasattr(operation, "program") else operation
    dim = 1 << width
    block = np.zeros((dim, dim), dtype=complex)
    leakage = 0.0
    for column in range(dim):
        state = simulate(program, initial={"target": column})
        for (t, sig), amplitude in state.amplitudes.items():
            if sig == 0:
                block[t, column] = amplitude
            else:
                leakage = max(leakage, abs(amplitude))
    return block, leakage


def cross_columns(operation, columns, paths=("rir", "adapter")):
    """basis_program 驱动下各后端路径与 reference 的逐列最大偏差。"""
    program = operation.program()
    widths = [r.type.width for r in program.main.registers]
    runners = {"rir": rir_pysparq, "adapter": adapter_pysparq, "origin": None}
    worst = 0.0
    for column in columns:
        p = basis_program(operation, {"target": column})
        ref = reference(p)
        for name in paths:
            if name == "origin":
                vec = originir_ext(p)
                deviation = statevector_error(vec, amplitudes_to_statevector(ref, widths))
            else:
                deviation = amplitude_error(runners[name](p), ref)
            worst = max(worst, deviation)
    return worst


def sampled_columns(width, count=8):
    """宽寄存器的列抽样：边界列 + 确定性伪随机列（与 harness.sampled_inputs 同种子）。"""
    from harness import sampled_inputs

    values, _ = sampled_inputs(width, samples=count)
    return values


def unitary_block(be):
    """originir_unitary + effective_block：全幺正的有效块与泄漏（小规模）。"""
    unitary = originir_unitary(be.operation.program())
    return effective_block(unitary, be.width)


def driver(operation, initial=None):
    """带 QRAM 资源的驱动程序：初态 X 置位 + 资源直通（harness.basis_program 无资源版）。"""
    b = Builder(
        "verify_driver",
        {r.name: r.type for r in operation.module.registers},
        {r.name: r.type for r in operation.module.resources},
    )
    for key, value in (initial or {}).items():
        for bit in range(b[key].width):
            if (value >> bit) & 1:
                b.x(b[key][bit])
    b.call(
        operation,
        **{r.name: b[r.name] for r in operation.module.registers},
        resources={r.name: r.name for r in operation.module.resources},
    )
    return b.finish().program()


def tridiagonal(dim, alpha, beta):
    """经典三对角矩阵 alpha*I + beta*T（numpy 独立参考）。"""
    a = np.zeros((dim, dim))
    for i in range(dim):
        a[i, i] = alpha
        if i > 0:
            a[i - 1, i] = beta
        if i < dim - 1:
            a[i + 1, i] = beta
    return a


def structural_permutations(matrix, sparsity_rows):
    """按每列结构位置构造 CKS 完整置换扩张（前 s 项为结构行，其余任意补全）。"""
    dim = len(matrix)
    permutations = []
    for column in range(dim):
        rows = list(sparsity_rows[column])
        permutations.append(rows + [r for r in range(dim) if r not in rows])
    return permutations


def tridiagonal_rows(dim):
    """三对角矩阵每列的结构行；边界列用条目为零的非结构行补足 s 个（位置置换须互异）。"""
    rows = []
    for j in range(dim):
        column_rows = [r for r in (j - 1, j, j + 1) if 0 <= r < dim]
        padding = (r for r in range(dim) if r not in column_rows)
        while len(column_rows) < min(3, dim):
            column_rows.append(next(padding))  # 补足行的矩阵条目为 0，幅度转导后置零
        rows.append(column_rows)
    return rows


def sparse_be_for(matrix, fmt, amax, sparsity_rows, *, signed, entry_database=None):
    """由稠密对称矩阵组装 gate 版 CKS 稀疏块编码（access 层逐列结构位置）。"""
    dim = len(matrix)
    n = (dim - 1).bit_length()
    sparsity = len(sparsity_rows[0])
    permutations = structural_permutations(matrix, sparsity_rows)
    if entry_database is None:
        table = {}
        for c in range(dim):
            for r in range(dim):
                word = fmt.encode(float(matrix[r, c]))
                if word:
                    table[r | (c << n)] = word
        entry_database = gate_database(2 * n, fmt.width, table)
    access = SparseAccess(
        sparse_location_gate(n, permutations),
        sparse_entry(entry_database, n),
        n,
        fmt.width,
        sparsity,
    )
    return real_symmetric_sparse_encoding(
        access, fmt, amax, diagonal_nonnegative=True
    )


def ps_tridiagonal_block(alpha, beta, n_bits):
    """pysparq BlockEncodingTridiagonal 的有效块：逐列 |j>|0> 演化后读 anc==0 振幅。"""
    import pysparq as ps
    from pysparq.algorithms.block_encoding import BlockEncodingTridiagonal

    ps.System.clear()
    try:
        ps.System.add_register("main_reg", ps.UnsignedInteger, n_bits)
        ps.System.add_register("anc_UA", ps.UnsignedInteger, 4)
        dim = 1 << n_bits
        main_id = ps.System.get_id("main_reg")
        anc_id = ps.System.get_id("anc_UA")
        block = np.zeros((dim, dim))
        for column in range(dim):
            state = ps.SparseState()
            ps.Init_Unsafe("main_reg", column)(state)
            ps.Init_Unsafe("anc_UA", 0)(state)
            BlockEncodingTridiagonal("main_reg", "anc_UA", alpha, beta)(state)
            for basis in state.basis_states:
                if basis.get(anc_id).value == 0:
                    block[basis.get(main_id).value, column] = basis.amplitude.real
        return block
    finally:
        ps.System.clear()


def ps_qram_block(matrix, n_bits, *, data_size=50, rational_size=51, exponent=15):
    """pysparq BlockEncodingViaQRAM 的有效块（配置同 C++ CorrectnessTest）。"""
    import pysparq as ps
    from pysparq.algorithms.block_encoding import BlockEncodingViaQRAM
    from pysparq.algorithms.qram_utils import make_vector_tree, scale_and_convert_vector

    dim = 1 << n_bits
    converted = scale_and_convert_vector(
        matrix.flatten().tolist(), exponent=exponent, data_size=data_size, from_matrix=True
    )
    tree = make_vector_tree(converted, data_size)
    qram = ps.QRAMCircuit_qutrit(2 * n_bits + 1, data_size, tree)
    ps.System.clear()
    try:
        ps.System.add_register("main_reg", ps.UnsignedInteger, n_bits)
        ps.System.add_register("anc_UA", ps.UnsignedInteger, n_bits)
        main_id = ps.System.get_id("main_reg")
        anc_id = ps.System.get_id("anc_UA")
        block = np.zeros((dim, dim))
        for column in range(dim):
            state = ps.SparseState()
            ps.Init_Unsafe("main_reg", column)(state)
            ps.Init_Unsafe("anc_UA", 0)(state)
            BlockEncodingViaQRAM(qram, "main_reg", "anc_UA", data_size, rational_size)(state)
            for basis in state.basis_states:
                if basis.get(anc_id).value == 0:
                    block[basis.get(main_id).value, column] = basis.amplitude.real
        return block
    finally:
        ps.System.clear()


# ---------------------------------------------------------------------------
# 对角块编码（论文点名 diagonal BE；lowrank._diagonal_encoding 与 DF 公开组装）
# ---------------------------------------------------------------------------


def verify_diagonal_be_unitary(report):
    """小规模对角 BE：originir_unitary 提取块，三后端逐列交叉。"""
    spectra = [
        (2, (0.7, -1.3)),
        (4, (0.5, -1.0, 0.25, 1.5)),
        (8, (0.6, -0.4, 0.0, 1.1, -0.9, 0.3, 0.2, -0.5)),  # 含零元与混合符号
    ]
    for dim, spectrum in spectra:
        be = _diagonal_encoding(spectrum)
        alpha = sum(abs(g) for g in spectrum)
        expected = np.diag(spectrum) / alpha
        block, leakage = unitary_block(be)
        error = float(np.abs(block - expected).max())
        cross = cross_columns(be.operation, range(dim))
        report.case(
            f"diagonal-be-unitary-d{dim}",
            paths=["originir-ext+to_matrix", "reference", "rir-pysparq", "adapter-pysparq"],
            parameters={"dim": dim, "spectrum": list(spectrum), "alpha": alpha},
            metrics={
                "max_error": error,
                "leakage": leakage,
                "cross_deviation": cross,
                "alpha": alpha,
            },
            criterion="有效块 == diag(g)/α（max_error < 1e-9），三后端逐列一致",
            passed=error < EXACT and cross < EXACT and abs(be.alpha - alpha) < 1e-12,
        )


def verify_diagonal_be_wide(report):
    """多规模对角 BE（d=16/32/64）：reference/rir/adapter 逐列对拍（超 OriginIR 预算）。"""
    for n in (4, 5, 6):
        dim = 1 << n
        spectrum = tuple(math.sin(0.7 * t + 0.3) + 0.2 * math.cos(1.3 * t) for t in range(dim))
        be = _diagonal_encoding(spectrum)
        alpha = sum(abs(g) for g in spectrum)
        expected = np.diag(spectrum) / alpha
        block, leakage = block_via_reference(be.operation, alpha, n)
        error = float(np.abs(block - expected).max())
        cross = cross_columns(be.operation, sampled_columns(n), paths=("rir",))
        report.case(
            f"diagonal-be-wide-d{dim}",
            paths=["reference", "rir-pysparq"],
            parameters={
                "dim": dim,
                "alpha": alpha,
                "note": "寄存器规模按位预算走 pysparq 路径；交叉列抽样 8 列（全列由 reference 覆盖）",
            },
            metrics={"max_error": error, "leakage": leakage, "cross_deviation": cross},
            criterion="逐列块 == diag(g)/α（max_error < 1e-9），后端间一致",
            passed=error < EXACT and cross < EXACT,
        )


def verify_diagonal_be_public_df(report):
    """同一对角谱的两条编码路径：_diagonal_encoding 与 double_factorized_encoding(U=I)。"""
    spectrum = (0.5, -1.0, 0.25, 1.5)
    direct = _diagonal_encoding(spectrum)
    identity4 = [[1.0 if i == j else 0.0 for j in range(4)] for i in range(4)]
    df = DoubleFactorization(0.0, [identity4, identity4], [spectrum, spectrum])
    assembled = double_factorized_encoding(df)
    expected = np.diag(spectrum) / direct.alpha
    block_a, _ = block_via_reference(direct.operation, direct.alpha, 2)
    block_b, _ = block_via_reference(assembled.operation, assembled.alpha, 2)
    error = max(
        float(np.abs(block_a - expected).max()),
        float(np.abs(block_b - expected).max()),
    )
    report.case(
        "diagonal-be-public-df-d4",
        paths=["reference"],
        parameters={"dim": 4, "alpha_direct": direct.alpha, "alpha_df": assembled.alpha},
        metrics={"max_error": error, "alpha_direct": direct.alpha, "alpha_df": assembled.alpha},
        criterion="公开 DF 组装（恒等旋转、双秩同谱）与直连对角 BE 块一致且 == diag(g)/α",
        passed=error < EXACT and abs(assembled.alpha - 2 * direct.alpha) < 1e-12,
    )


# ---------------------------------------------------------------------------
# BE 组合代数（block_encoding.py 与 operators.py 基元）
# ---------------------------------------------------------------------------


def verify_pauli_word_and_embeddings(report):
    """pauli_word（无信号、泄漏为零）与布尔嵌入 projector/truncated_shift。"""
    cases = [
        ("pauli-word-YZX", pauli_word("YZX"), pauli_matrix("YZX"), 3),
        ("projector-w2-03", projector(2, [0, 3]), np.diag([1, 0, 0, 1]).astype(complex), 2),
    ]
    shift = np.zeros((4, 4))
    for v in range(3):
        shift[v + 1, v] = 1.0
    cases.append(("truncated-shift-w2-l3", truncated_shift(2, 3), shift, 2))
    for name, be, expected, width in cases:
        block, leakage = unitary_block(be)
        error = float(np.abs(block - expected).max())
        cross = cross_columns(be.operation, range(1 << width))
        report.case(
            name,
            paths=["originir-ext+to_matrix", "reference", "rir-pysparq", "adapter-pysparq"],
            parameters={"width": width, "alpha": be.alpha},
            metrics={"max_error": error, "leakage": leakage, "cross_deviation": cross},
            criterion="有效块逐元等于经典矩阵（max_error < 1e-9）",
            passed=error < EXACT and cross < EXACT,
        )


def verify_matrix_pauli_encoding(report):
    """小矩阵显式 Pauli LCU：Hermitian 与非 Hermitian，alpha 与独立 l1 上界对拍。"""
    rng = np.random.default_rng(20260916)
    hermitian = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    hermitian = (hermitian + hermitian.conj().T) / 2
    general = np.array([[0.5, 0.2 - 0.1j], [-0.3j, 0.4]], dtype=complex)
    for name, matrix in (("2x2-general", general), ("4x4-hermitian", hermitian)):
        be = matrix_pauli_encoding(matrix)
        block, leakage = unitary_block(be)
        error = float(np.abs(block * be.alpha - matrix).max())
        l1 = pauli_l1(matrix)
        cross = cross_columns(be.operation, range(len(matrix)))
        report.case(
            f"matrix-pauli-encoding-{name}",
            paths=["originir-ext+to_matrix", "reference", "rir-pysparq", "adapter-pysparq"],
            parameters={"dim": len(matrix), "alpha": be.alpha, "pauli_l1_independent": l1},
            metrics={
                "max_error": error,
                "leakage": leakage,
                "cross_deviation": cross,
                "alpha_minus_l1": be.alpha - l1,
            },
            criterion="块*α == M（max_error < 1e-9）且 α 不超过独立 Pauli l1 上界（+1e-9 容差）",
            passed=error < EXACT and cross < EXACT and be.alpha <= l1 + 1e-9,
        )


def verify_be_algebra_combinators(report):
    """组合子逐项：tensor/adjoint/pad_signal/lcu/kronecker_sum/direct_sum/product。"""
    rng = np.random.default_rng(7)
    a_mat = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    a_mat = (a_mat + a_mat.conj().T) / 2
    b_mat = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    b_mat = (b_mat + b_mat.conj().T) / 2
    c_mat = np.array([[0.3, 0.5 + 0.2j], [0.1 - 0.4j, -0.6]])  # 非 Hermitian
    ea, eb, ec = (matrix_pauli_encoding(m) for m in (a_mat, b_mat, c_mat))
    ca, cb = 0.6 + 0.3j, -0.8
    entries = [
        ("tensor", tensor(ea, eb), np.kron(a_mat, b_mat)),
        ("adjoint", adjoint_be(ec), c_mat.conj().T),
        ("pad-signal", pad_signal(ea, ea.signal_qubits + 2), a_mat),
        ("lcu-complex", lcu([(ca, ea), (cb, eb)]), ca * a_mat + cb * b_mat),
        (
            "kronecker-sum",
            kronecker_sum(ea, eb),
            np.kron(a_mat, np.eye(2)) + np.kron(np.eye(2), b_mat),
        ),
        ("direct-sum", direct_sum(ea, eb), np.block([[a_mat, np.zeros((2, 2))], [np.zeros((2, 2)), b_mat]])),
        ("product", product(ea, ec), a_mat @ c_mat),
    ]
    for name, be, expected in entries:
        block, leakage = unitary_block(be)
        error = float(np.abs(block * be.alpha - expected).max())
        report.case(
            f"be-algebra-{name}",
            paths=["originir-ext+to_matrix", "reference"],
            parameters={"alpha": be.alpha, "signal_qubits": be.signal_qubits},
            metrics={"max_error": error, "leakage": leakage, "alpha": be.alpha},
            criterion="块*α 逐元等于组合语义矩阵（max_error < 1e-9）",
            passed=error < EXACT,
        )


# ---------------------------------------------------------------------------
# PREPARE–SELECT（prepare_select.py）
# ---------------------------------------------------------------------------


def verify_gate_prepare_distribution(report):
    """gate PREPARE 振幅 == sqrt(|c|/α)，含补零槽位；四路径对拍。"""
    for name, coefficients in (
        ("terms4", [0.7, -0.4, 1.1, 0.2]),
        ("terms5-padded", [0.5, 0.25, -0.75, 1.0, 0.1]),
    ):
        alpha = sum(abs(c) for c in coefficients)
        prep = gate_prepare(coefficients)
        program = prep.operation.program()
        expected = {i: math.sqrt(abs(c) / alpha) for i, c in enumerate(coefficients)}
        state = reference(program)
        error = max(
            abs(abs(state.get((i, 0), 0j)) - amp) for i, amp in expected.items()
        )
        padded = (1 << prep.width) - len(coefficients)
        for i in range(len(coefficients), 1 << prep.width):
            error = max(error, abs(state.get((i, 0), 0j)))
        cross = max(
            amplitude_error(rir_pysparq(program), state),
            amplitude_error(adapter_pysparq(program), state),
        )
        report.case(
            f"gate-prepare-distribution-{name}",
            paths=["reference", "rir-pysparq", "adapter-pysparq"],
            parameters={"terms": len(coefficients), "selector_width": prep.width, "padded_slots": padded},
            metrics={"max_amplitude_error": error, "cross_deviation": cross},
            criterion="selector 振幅逐点等于 sqrt(|c|/α)，补零槽位振幅为零（< 1e-9）",
            passed=error < EXACT and cross < EXACT,
        )


def verify_select_pauli(report):
    """SELECT 逐 selector 值施加对应 Pauli 字与系数相位（稠密对拍）。"""
    terms = [(0.7, "XI"), (-0.4, "ZZ"), (1.1, "IY"), (0.2j, "YX")]
    operation = select_pauli(terms)
    program = operation.program()
    worst = 0.0
    for index, (coefficient, word) in enumerate(terms):
        expected_op = np.exp(1j * np.angle(coefficient)) * pauli_matrix(word)
        for t in range(4):
            state = simulate(program, initial={"selector": index, "target": t})
            (out_sel, out_t), amplitude = next(iter(state.amplitudes.items()))
            got = np.zeros(4, dtype=complex)
            got[out_t] = amplitude
            worst = max(worst, float(np.abs(got - expected_op[:, t]).max()))
            if out_sel != index:
                worst = max(worst, 1.0)
    report.case(
        "select-pauli-indexed-words",
        paths=["reference"],
        parameters={"terms": len(terms), "selector_width": 2},
        metrics={"max_error": worst},
        criterion="selector==i 时 target 上恰为 e^{iφ_i}·P_i（max_error < 1e-9）",
        passed=worst < EXACT,
    )


def verify_lcu_prepare_select_block(report):
    """PREPARE–SELECT 块编码：多规模、多系数（含复相位），(0,0) 块 == H/α。"""
    families = [
        ("w1-balanced", [(0.6, "X"), (-0.8, "Z")], 1, True),
        ("w2-mixed", [(0.7, "XI"), (-0.4, "ZZ"), (1.1, "IY")], 2, True),
        ("w2-complex", [(0.5 + 0.2j, "XY"), (0.3, "II"), (-0.4j, "ZY"), (0.9, "XZ")], 2, False),
        ("w3-six-terms", [(0.4, "XII"), (0.3, "IZI"), (-0.5, "ZXZ"), (0.2j, "YYY"), (0.6, "IXX"), (-0.35, "ZZI")], 3, False),
    ]
    for name, terms, width, use_unitary in families:
        be = lcu_prepare_select(terms)
        hamiltonian = sum(c * pauli_matrix(w) for c, w in terms)
        alpha = sum(abs(c) for c, _ in terms)
        paths = ["reference", "rir-pysparq", "adapter-pysparq"]
        block, leakage = block_via_reference(be.operation, alpha, width)
        error = float(np.abs(block * alpha - hamiltonian).max())
        cross = cross_columns(be.operation, range(1 << width))
        unitary_error = None
        if use_unitary:
            ublock, uleak = unitary_block(be)
            unitary_error = float(np.abs(ublock * alpha - hamiltonian).max())
            leakage = max(leakage, uleak)
            paths = ["originir-ext+to_matrix", *paths]
        metrics = {"max_error": error, "leakage": leakage, "cross_deviation": cross, "alpha": alpha}
        if unitary_error is not None:
            metrics["unitary_error"] = unitary_error
        report.case(
            f"lcu-prepare-select-block-{name}",
            paths=paths,
            parameters={"width": width, "terms": len(terms), "alpha": alpha},
            metrics=metrics,
            criterion="块*α == Σc_i P_i（max_error < 1e-9），各后端一致",
            passed=error < EXACT and cross < EXACT and (unitary_error is None or unitary_error < EXACT),
        )


def verify_qram_prepare_quantization(report):
    """QRAM PREPARE：角表量化误差随 angle_width 收敛；reference 与 rir 双路径。"""
    coefficients = [0.7, -0.4, 1.1, 0.2]
    alpha = sum(abs(c) for c in coefficients)
    expected = {i: abs(c) / alpha for i, c in enumerate(coefficients)}
    for angle_width, bound in ((8, 0.02), (12, 0.002)):
        qp = qram_prepare(coefficients, angle_width=angle_width)
        program = driver(qp.state_preparation().operation)
        state = reference(program, qp.memory)
        probs = {k[0]: abs(v) ** 2 for k, v in state.items()}
        tvd = 0.5 * sum(abs(probs.get(i, 0.0) - p) for i, p in expected.items())
        cross = amplitude_error(rir_pysparq(program, qp.memory), state)
        report.case(
            f"qram-prepare-quantization-a{angle_width}",
            paths=["reference", "rir-pysparq"],
            parameters={"angle_width": angle_width},
            metrics={"tvd": tvd, "bound": bound, "cross_deviation": cross},
            criterion=f"量化分布与精确分布的 TVD ≤ {bound}（仓库 QRAM 约定容差）且两后端逐振幅一致",
            passed=tvd <= bound and cross < EXACT,
        )


def verify_alias_prepare(report):
    """alias 采样 PREPARE：分布 TVD 不超过 2^selector·2^-precision，rir 精确一致。"""
    coefficients = [0.7, -0.4, 1.1, 0.2]
    alpha = sum(abs(c) for c in coefficients)
    expected = {i: abs(c) / alpha for i, c in enumerate(coefficients)}
    for precision in (6, 10):
        ap = alias_prepare(coefficients, precision=precision)
        program = driver(ap.state_preparation().operation)
        state = reference(program, ap.memory)
        marginal = {}
        for (t, _w), amplitude in state.items():
            marginal[t] = marginal.get(t, 0.0) + abs(amplitude) ** 2
        tvd = 0.5 * sum(abs(marginal.get(i, 0.0) - p) for i, p in expected.items())
        bound = 4 * 2.0**-precision  # 文档解析界：TVD ≤ 2^selector·2^-precision
        cross = amplitude_error(rir_pysparq(program, ap.memory), state)
        report.case(
            f"alias-prepare-distribution-p{precision}",
            paths=["reference", "rir-pysparq"],
            parameters={"precision": precision, "selector_width": 2},
            metrics={"tvd": tvd, "bound": bound, "cross_deviation": cross},
            criterion="target 边缘分布 TVD ≤ 2^2·2^-precision 且 rir 与 reference 逐振幅一致",
            passed=tvd <= bound and cross < EXACT,
        )


def verify_alias_prepare_select_block(report):
    """alias PREPARE 组装的 PREPARE–SELECT BE：块误差受量化界控制，两后端一致。"""
    terms = [(0.7, "XI"), (-0.4, "ZZ"), (1.1, "IY")]
    coefficients = [c for c, _ in terms]
    alpha = sum(abs(c) for c in coefficients)
    ap = alias_prepare(coefficients, precision=10)
    be = lcu_prepare_select(terms, prepare=ap)
    resource = be.operation.module.resources[0].name
    memory = {resource: next(iter(ap.memory.values()))}
    hamiltonian = sum(c * pauli_matrix(w) for c, w in terms)
    bound = 4 * 2.0**-10  # TVD 界 ×（≤2 的振幅因子）宽松取
    error, cross = 0.0, 0.0
    for column in range(4):
        program = driver(be.operation, {"target": column})
        ref = reference(program, memory)
        cross = max(cross, amplitude_error(rir_pysparq(program, memory), ref))
        for (t, sig), amplitude in ref.items():
            if sig == 0:
                error = max(error, abs(amplitude * alpha - hamiltonian[t, column]))
    report.case(
        "alias-prepare-select-block",
        paths=["reference", "rir-pysparq"],
        parameters={"precision": 10, "alpha": alpha},
        metrics={"max_error": error, "bound": bound, "cross_deviation": cross},
        criterion="块*α 与 H 的偏差 ≤ 4·2^-10（量化界）且后端一致",
        passed=error <= bound and cross < EXACT,
    )


def verify_abstract_prepare_bind(report):
    """开放 PREPARE 声明：在 BE 内 bind gate 实现后块与直接组装逐列一致。"""
    terms = [(0.7, "XI"), (-0.4, "ZZ"), (1.1, "IY")]
    coefficients = [c for c, _ in terms]
    alpha = sum(abs(c) for c in coefficients)
    be_abstract = lcu_prepare_select(terms, prepare=abstract_prepare(coefficients, work_width=0))
    slots = [r.name for r in unresolved(be_abstract.operation.program())]
    if len(slots) != 1:
        raise AssertionError(f"抽象 BE 应恰有一个未绑定槽位：{slots}")
    bound = bind(be_abstract.operation.program(), {slots[0]: gate_prepare(coefficients).operation})
    hamiltonian = sum(c * pauli_matrix(w) for c, w in terms)
    be_gate = lcu_prepare_select(terms)
    error, deviation = 0.0, 0.0
    for column in range(4):
        bound_state = simulate(bound, initial={"target": column})
        gate_state = simulate(be_gate.operation.program(), initial={"target": column})
        deviation = max(deviation, amplitude_error(dict(bound_state.amplitudes), dict(gate_state.amplitudes)))
        for (t, sig), amplitude in bound_state.amplitudes.items():
            if sig == 0:
                error = max(error, abs(amplitude * alpha - hamiltonian[t, column]))
    report.case(
        "abstract-prepare-bind-in-be",
        paths=["reference"],
        parameters={"slot": slots[0], "alpha": alpha},
        metrics={"max_error": error, "bind_deviation": deviation},
        criterion="绑定后块*α == H 且与 gate 直接组装逐振幅一致（< 1e-9）",
        passed=error < EXACT and deviation < EXACT,
    )


# ---------------------------------------------------------------------------
# 稀疏访问辅助（sparse.py）与稀疏块编码
# ---------------------------------------------------------------------------


def verify_sparse_rotation_helpers(report):
    """word_rotation / magnitude_rotation 的概率语义穷举（reference + rir 交叉）。"""
    worst_word = 0.0
    for width in (2, 3, 4):
        operation = word_rotation(width)
        scale = 2 * math.pi / (1 << width)
        for value in range(1 << width):
            state = simulate(operation.program(), initial={"value": value})
            p1 = sum(abs(a) ** 2 for k, a in state.amplitudes.items() if k[1] == 1)
            worst_word = max(worst_word, abs(p1 - math.sin(scale * value / 2) ** 2))
    fmt = FixedFormat(4, 1, signed=True)
    amax = 3.0
    rotation = magnitude_rotation(fmt, amax)
    worst_mag = 0.0
    for raw in range(1 << fmt.width):
        state = simulate(rotation.program(), initial={"value": raw})
        p0 = sum(abs(a) ** 2 for k, a in state.amplitudes.items() if k[1] == 0)
        worst_mag = max(worst_mag, abs(p0 - min(1.0, abs(fmt.decode(raw)) / amax)))
    cross = 0.0
    for raw in range(1 << fmt.width):
        p = basis_program(rotation, {"value": raw})
        cross = max(cross, amplitude_error(reference(p), rir_pysparq(p)))
    report.case(
        "sparse-rotation-helpers",
        paths=["reference", "rir-pysparq"],
        parameters={"word_widths": [2, 3, 4], "fmt": "4.1s", "amax": amax},
        metrics={
            "word_rotation_max_error": worst_word,
            "magnitude_rotation_max_error": worst_mag,
            "cross_deviation": cross,
        },
        criterion="P(flag=1)==sin²(θ_v/2)、P(成功)==min(1,|v|/amax) 逐点成立（< 1e-9）",
        passed=worst_word < EXACT and worst_mag < EXACT and cross < EXACT,
    )


def verify_sparse_boolean_helpers(report):
    """compare_words / value_transposition 全输入穷举；prefix_state 均匀性。"""
    mismatches = 0
    checked = 0
    for width in (1, 2, 3, 4):
        for kind in ("eq", "lt"):
            operation = compare_words(width, kind)
            for a in range(1 << width):
                for b in range(1 << width):
                    state = simulate(operation.program(), initial={"a": a, "b": b})
                    flag = next(iter(state.amplitudes))[2]
                    expect = int(a == b) if kind == "eq" else int(a < b)
                    mismatches += flag != expect
                    checked += 1
    trans = value_transposition(3)
    for index in range(8):
        for a in range(8):
            for b in range(8):
                state = simulate(trans.program(), initial={"index": index, "a": a, "b": b})
                (oi, oa, ob), _ = next(iter(state.amplitudes.items()))
                expect_index = b if index == a else (a if index == b else index)
                mismatches += (oi, oa, ob) != (expect_index, a, b)
                checked += 1
    prefix_error = 0.0
    for width, count in ((1, 1), (2, 3), (3, 5), (3, 8), (4, 13)):
        state = simulate(prefix_state(width, count).program())
        if len(state.amplitudes) != count:
            prefix_error = max(prefix_error, 1.0)
        for amplitude in state.amplitudes.values():
            prefix_error = max(prefix_error, abs(abs(amplitude) ** 2 - 1 / count))
    report.case(
        "sparse-boolean-helpers",
        paths=["reference"],
        parameters={"boolean_checked": checked, "prefix_cases": [(1, 1), (2, 3), (3, 5), (3, 8), (4, 13)]},
        metrics={"mismatches": mismatches, "prefix_max_error": prefix_error},
        criterion="布尔网络全输入语义正确（mismatches == 0）且前缀叠加均匀（< 1e-12）",
        passed=mismatches == 0 and prefix_error < 1e-12,
    )


def verify_sparse_lookup_helpers(report):
    """reversible_lookup 融合视图与 batch_lookup 的 XOR 语义（含非零 data 初值）。"""
    table = {i: (i * i + 3) % 16 for i in range(8)}
    database = gate_database(3, 4, table)
    lookup = reversible_lookup({"x": 2, "y": 1}, {"lo": 2, "hi": 2}, database)
    mismatches = 0
    for x in range(4):
        for y in range(2):
            for lo in range(4):
                for hi in range(4):
                    state = simulate(
                        lookup.program(), initial={"x": x, "y": y, "lo": lo, "hi": hi}
                    )
                    (ox, oy, olo, ohi), _ = next(iter(state.amplitudes.items()))
                    word = table[(y << 2) | x]
                    mismatches += (ox, oy, olo, ohi) != (
                        x,
                        y,
                        lo ^ (word & 3),
                        hi ^ (word >> 2),
                    )
    batch = batch_lookup(database, 3)
    state = simulate(
        batch.program(),
        initial={"address0": 1, "address1": 2, "address2": 3, "data0": 5, "data1": 6, "data2": 7},
    )
    key = next(iter(state.amplitudes))
    batch_ok = key == (1, 5 ^ table[1], 2, 6 ^ table[2], 3, 7 ^ table[3])
    report.case(
        "sparse-lookup-helpers",
        paths=["reference"],
        parameters={"table_size": 8, "lookup_inputs": 128},
        metrics={"mismatches": mismatches, "batch_ok": int(batch_ok)},
        criterion="融合地址/数据视图 XOR 语义逐点正确（mismatches == 0）",
        passed=mismatches == 0 and batch_ok,
    )


def verify_sparse_access_layer(report):
    """CKS 访问层：位置置换逐列对拍、元素 XOR（非零 data）、QRAM 位置复净。"""
    matrix = tridiagonal(4, 1.5, -0.5)
    rows = tridiagonal_rows(4)
    permutations = structural_permutations(matrix, rows)
    location = sparse_location_gate(2, permutations)
    mismatches = 0
    for column in range(4):
        for index in range(4):
            state = simulate(location.program(), initial={"column": column, "index": index})
            (oc, oi, ow), _ = next(iter(state.amplitudes.items()))
            mismatches += (oc, oi, ow) != (column, permutations[column][index], 0)
    fmt = FixedFormat(3, 1, signed=True)
    table = {}
    for c in range(4):
        for r in range(4):
            word = fmt.encode(float(matrix[r, c]))
            if word:
                table[r | (c << 2)] = word
    entry = sparse_entry(gate_database(4, fmt.width, table), 2)
    for row in range(4):
        for column in range(4):
            for data in (0, 5):
                state = simulate(
                    entry.program(), initial={"row": row, "column": column, "data": data}
                )
                (or_, oc, od), _ = next(iter(state.amplitudes.items()))
                mismatches += (or_, oc, od) != (
                    row,
                    column,
                    data ^ fmt.encode(float(matrix[row, column])),
                )
    # QRAM 位置实现：正反思表绑定后 work 复净且 index 持新值（reference + rir）
    forward = {}
    inverse = {}
    for column in range(4):
        for index in range(4):
            forward[column | (index << 2)] = permutations[column][index]
            inverse[column | (permutations[column][index] << 2)] = index
    qlocation = sparse_location_qram(2)
    memory = {"forward": forward, "inverse": inverse}
    cross = 0.0
    qmismatch = 0
    for column in range(4):
        for index in range(4):
            program = driver(qlocation, {"column": column, "index": index})
            ref = reference(program, memory)
            cross = max(cross, amplitude_error(rir_pysparq(program, memory), ref))
            (oc, oi, ow), _ = next(iter(ref.items()))
            qmismatch += (oc, oi, ow) != (column, permutations[column][index], 0)
    report.case(
        "sparse-access-layer",
        paths=["reference", "rir-pysparq"],
        parameters={"dim": 4, "gate_checks": 16 + 32, "qram_checks": 16},
        metrics={
            "location_mismatches": mismatches,
            "qram_mismatches": qmismatch,
            "qram_cross_deviation": cross,
        },
        criterion="位置/元素 oracle 逐基态语义正确，QRAM 位置 work 复净且两后端一致",
        passed=mismatches == 0 and qmismatch == 0 and cross < EXACT,
    )


def verify_sparse_be_signed(report):
    """带符号稀疏 BE（CKS T†ST）：三对角多规模；d2 走 OriginIR 态向量提取。"""
    # d=2：15 qubits，OriginIR 态向量预算内
    fmt = FixedFormat(3, 1, signed=True)
    matrix2 = tridiagonal(2, 1.5, -0.5)
    be2 = sparse_be_for(matrix2, fmt, 1.5, [[0, 1], [0, 1]], signed=True)
    block2, leakage2 = block_via_reference(be2.operation, be2.alpha, 1)
    error2 = float(np.abs(block2 * be2.alpha - matrix2).max())
    cross2 = cross_columns(be2.operation, range(2), paths=("rir", "adapter", "origin"))
    report.case(
        "sparse-be-signed-tridiagonal-d2",
        paths=["originir-ext", "reference", "rir-pysparq", "adapter-pysparq"],
        parameters={"dim": 2, "sparsity": 2, "alpha": be2.alpha, "originir_qubits": 15},
        metrics={"max_error": error2, "leakage": leakage2, "cross_deviation": cross2, "alpha": be2.alpha},
        criterion="块*α == A（max_error < 1e-9），OriginIR 态向量与稀疏路径一致",
        passed=error2 < EXACT and cross2 < EXACT and abs(be2.alpha - 2 * 1.5) < 1e-12,
    )
    # d=4：24 qubits 达到 OriginIR 预算，但门数使单列态向量模拟 ~68 s，按预算只走 pysparq 路径
    matrix4 = tridiagonal(4, 1.5, -0.5)
    be4 = sparse_be_for(matrix4, fmt, 1.5, tridiagonal_rows(4), signed=True)
    block4, leakage4 = block_via_reference(be4.operation, be4.alpha, 2)
    error4 = float(np.abs(block4 * be4.alpha - matrix4).max())
    cross4 = cross_columns(be4.operation, range(4))
    report.case(
        "sparse-be-signed-tridiagonal-d4",
        paths=["reference", "rir-pysparq", "adapter-pysparq"],
        parameters={
            "dim": 4,
            "sparsity": 3,
            "alpha": be4.alpha,
            "note": "OriginIR 导出 24 qubits 且门数大（单列态向量 ~68 s），按运行时预算只走 pysparq/reference 路径",
        },
        metrics={"max_error": error4, "leakage": leakage4, "cross_deviation": cross4, "alpha": be4.alpha},
        criterion="块*α == A（max_error < 1e-9），三后端一致",
        passed=error4 < EXACT and cross4 < EXACT and abs(be4.alpha - 3 * 1.5) < 1e-12,
    )
    # d=8：34 qubits 超 OriginIR 预算
    matrix8 = tridiagonal(8, 1.5, -0.5)
    be8 = sparse_be_for(matrix8, fmt, 1.5, tridiagonal_rows(8), signed=True)
    block8, leakage8 = block_via_reference(be8.operation, be8.alpha, 3)
    error8 = float(np.abs(block8 * be8.alpha - matrix8).max())
    cross8 = cross_columns(be8.operation, [0, 1, 3, 7], paths=("rir",))
    report.case(
        "sparse-be-signed-tridiagonal-d8",
        paths=["reference", "rir-pysparq"],
        parameters={
            "dim": 8,
            "sparsity": 3,
            "alpha": be8.alpha,
            "note": "34 qubits 超 OriginIR 24 位预算；交叉列抽样 4 列（全列由 reference 覆盖）",
        },
        metrics={"max_error": error8, "leakage": leakage8, "cross_deviation": cross8, "alpha": be8.alpha},
        criterion="块*α == A（max_error < 1e-9），两后端一致",
        passed=error8 < EXACT and cross8 < EXACT,
    )


def verify_sparse_be_unsigned(report):
    """无符号稀疏 BE：d4 用 originir_unitary 提取全幺正块；d16 多规模 pysparq 路径。"""
    fmt = FixedFormat(3, 1, signed=False)
    matrix4 = tridiagonal(4, 1.5, 0.5)
    be4 = sparse_be_for(matrix4, fmt, 1.5, tridiagonal_rows(4), signed=False)
    block4, leakage4 = block_via_reference(be4.operation, be4.alpha, 2)
    error4 = float(np.abs(block4 * be4.alpha - matrix4).max())
    # 12 qubits 在 OriginIR 态向量预算内逐列提取；to_matrix 因门数多（~65 s/次）超出运行时预算
    cross4 = cross_columns(be4.operation, range(4), paths=("rir", "adapter", "origin"))
    report.case(
        "sparse-be-unsigned-unitary-d4",
        paths=["originir-ext", "reference", "rir-pysparq", "adapter-pysparq"],
        parameters={"dim": 4, "sparsity": 3, "alpha": be4.alpha, "originir_qubits": 12},
        metrics={"max_error": error4, "leakage": leakage4, "cross_deviation": cross4, "alpha": be4.alpha},
        criterion="块*α == A（max_error < 1e-9），OriginIR 态向量与稀疏路径逐列一致",
        passed=error4 < EXACT and cross4 < EXACT,
    )
    matrix16 = tridiagonal(16, 1.5, 0.5)
    be16 = sparse_be_for(matrix16, fmt, 1.5, tridiagonal_rows(16), signed=False)
    block16, leakage16 = block_via_reference(be16.operation, be16.alpha, 4)
    error16 = float(np.abs(block16 * be16.alpha - matrix16).max())
    cross16 = cross_columns(be16.operation, sampled_columns(4), paths=("rir",))
    report.case(
        "sparse-be-unsigned-d16",
        paths=["reference", "rir-pysparq"],
        parameters={
            "dim": 16,
            "sparsity": 3,
            "alpha": be16.alpha,
            "note": "超 OriginIR 位预算；交叉列抽样 8 列（全列由 reference 覆盖）",
        },
        metrics={"max_error": error16, "leakage": leakage16, "cross_deviation": cross16, "alpha": be16.alpha},
        criterion="块*α == A（max_error < 1e-9），两后端一致",
        passed=error16 < EXACT and cross16 < EXACT,
    )


def verify_sparse_be_qram_access(report):
    """QRAM 数据绑定的稀疏 BE：位置正反表 + 元素表全在 memory，两后端对拍。"""
    fmt = FixedFormat(3, 1, signed=True)
    matrix = tridiagonal(4, 1.5, -0.5)
    rows = tridiagonal_rows(4)
    permutations = structural_permutations(matrix, rows)
    entry_database = qram_database(4, fmt.width)
    access = SparseAccess(
        sparse_location_qram(2), sparse_entry(entry_database, 2), 2, fmt.width, 3
    )
    be = real_symmetric_sparse_encoding(access, fmt, 1.5, diagonal_nonnegative=True)
    forward, inverse, entries = {}, {}, {}
    for column in range(4):
        for index in range(4):
            forward[column | (index << 2)] = permutations[column][index]
            inverse[column | (permutations[column][index] << 2)] = index
        for row in range(4):
            word = fmt.encode(float(matrix[row, column]))
            if word:
                entries[row | (column << 2)] = word
    logical = {"forward": forward, "inverse": inverse, "table": entries}
    memory = {}
    for resource in be.operation.module.resources:
        for key, value in logical.items():
            if resource.name.endswith(key):
                memory[resource.name] = value
    if len(memory) != 3:
        raise AssertionError(f"QRAM 资源映射不全：{list(memory)}")
    error, cross = 0.0, 0.0
    leakage = 0.0
    for column in range(4):
        program = driver(be.operation, {"target": column})
        ref = reference(program, memory)
        cross = max(cross, amplitude_error(rir_pysparq(program, memory), ref))
        for (t, sig), amplitude in ref.items():
            if sig == 0:
                error = max(error, abs(amplitude * be.alpha - matrix[t, column]))
            else:
                leakage = max(leakage, abs(amplitude))
    report.case(
        "sparse-be-qram-access-d4",
        paths=["reference", "rir-pysparq"],
        parameters={"dim": 4, "sparsity": 3, "alpha": be.alpha, "resources": sorted(memory)},
        metrics={"max_error": error, "leakage": leakage, "cross_deviation": cross},
        criterion="QRAM 数据绑定下块*α == A（max_error < 1e-9）且两后端逐振幅一致",
        passed=error < EXACT and cross < EXACT,
    )


def verify_chebyshev_walk(report):
    """chebyshev_block：自伴酉扩张行走幂的零信号块 == T_k(A/α)。"""
    fmt = FixedFormat(3, 1, signed=True)
    matrix = np.array([[1.5, -0.5], [-0.5, 1.5]])
    be = sparse_be_for(matrix, fmt, 1.5, [[0, 1], [0, 1]], signed=True)
    h = matrix / be.alpha
    polynomials = {
        1: h,
        2: 2 * h @ h - np.eye(2),
        3: 4 * h @ h @ h - 3 * h,
        4: 8 * h @ h @ h @ h - 8 * h @ h + np.eye(2),
    }
    for degree, expected in polynomials.items():
        walk = chebyshev_block(be, degree)
        block, leakage = block_via_reference(walk.operation, walk.alpha, 1)
        error = float(np.abs(block - expected).max())
        cross = cross_columns(walk.operation, range(2)) if degree <= 2 else 0.0
        report.case(
            f"chebyshev-walk-k{degree}",
            paths=["reference", "rir-pysparq", "adapter-pysparq"],
            parameters={"degree": degree, "argument_scale": be.alpha},
            metrics={"max_error": error, "leakage": leakage, "cross_deviation": cross},
            criterion=f"行走幂零信号块 == T_{degree}(A/α)（max_error < 1e-9）",
            passed=error < EXACT and cross < EXACT,
        )


def verify_select_swap(report):
    """Select-Swap QROM：λ 全扫描对拍 gate_database；并作为稀疏 BE 的 entry 端到端。"""
    rng = np.random.default_rng(3)
    words = [int(v) for v in rng.integers(0, 1 << 3, size=16)]
    baseline = gate_database(4, 3, {i: v for i, v in enumerate(words) if v})
    for partitions in (1, 2, 4, 8, 16):
        database = select_swap_qrom(words, partitions=partitions, data_bits=3)
        worst = 0.0
        for address in range(16):
            for initial_data in (0, 5):
                a = simulate(
                    driver(database.operation, {"address": address, "data": initial_data})
                )
                b = simulate(
                    driver(baseline.operation, {"address": address, "data": initial_data})
                )
                worst = max(worst, amplitude_error(dict(a.amplitudes), dict(b.amplitudes)))
        report.case(
            f"select-swap-lambda{partitions}",
            paths=["reference"],
            parameters={"partitions": partitions, "addresses": 16, "data_bits": 3},
            metrics={"max_deviation": worst},
            criterion="全部地址 × data 初值读出与 gate_database 基线逐振幅一致（< 1e-9）",
            passed=worst < EXACT,
        )
    # 端到端：select_swap 作为稀疏块编码的元素数据库
    fmt = FixedFormat(3, 1, signed=True)
    matrix = np.array([[1.5, -0.5], [-0.5, 1.5]])
    table = {}
    for c in range(2):
        for r in range(2):
            word = fmt.encode(float(matrix[r, c]))
            if word:
                table[r | (c << 1)] = word
    database = select_swap_qrom([table.get(i, 0) for i in range(4)], partitions=2, data_bits=fmt.width)
    be = sparse_be_for(matrix, fmt, 1.5, [[0, 1], [0, 1]], signed=True, entry_database=database)
    block, leakage = block_via_reference(be.operation, be.alpha, 1)
    error = float(np.abs(block * be.alpha - matrix).max())
    report.case(
        "select-swap-as-sparse-entry",
        paths=["reference"],
        parameters={"dim": 2, "partitions": 2, "alpha": be.alpha},
        metrics={"max_error": error, "leakage": leakage},
        criterion="select_swap 作 entry 的稀疏 BE 块*α == A（max_error < 1e-9）",
        passed=error < EXACT,
    )


# ---------------------------------------------------------------------------
# 低秩 DF/THC（lowrank.py）
# ---------------------------------------------------------------------------


def verify_diagonalize_symmetric(report):
    """Jacobi 特征分解：与 numpy.linalg.eigh 独立对拍，重构 G == V diag(λ) Vᵀ。"""
    rng = np.random.default_rng(11)
    for dim in (2, 4, 8, 16, 32):
        raw = rng.normal(size=(dim, dim))
        matrix = (raw + raw.T) / 2
        eigenvalues, vectors = diagonalize_symmetric(matrix.tolist())
        v = np.array(vectors)
        lam = np.array(eigenvalues)
        reconstruct_error = float(np.abs(v @ np.diag(lam) @ v.T - matrix).max())
        reference_values = np.linalg.eigvalsh(matrix)
        eigen_error = float(np.abs(np.sort(lam) - reference_values).max())
        orthogonality = float(np.abs(v.T @ v - np.eye(dim)).max())
        report.case(
            f"diagonalize-symmetric-d{dim}",
            paths=["classical(numpy-independent)"],
            parameters={"dim": dim},
            metrics={
                "reconstruct_error": reconstruct_error,
                "eigenvalue_error": eigen_error,
                "orthogonality_error": orthogonality,
            },
            criterion="重构/特征值/正交性误差均 < 1e-9（对照 numpy.linalg.eigh）",
            passed=max(reconstruct_error, eigen_error, orthogonality) < EXACT,
        )


def _random_orthonormal(dim, rng):
    """numpy QR 独立生成酉/正交旋转矩阵。"""
    q, r = np.linalg.qr(rng.normal(size=(dim, dim)))
    signs = np.sign(np.diag(r))
    return (q * signs).tolist()


def verify_double_factorization(report):
    """DF 块编码：多秩、多规模；originir_unitary 提取块，α 与独立口径对拍。"""
    rng = np.random.default_rng(23)
    # 2x2 多秩（含 scalar 与 from_symmetric 预处理）
    g1 = [[2.0, 0.5], [0.5, 1.0]]
    g2 = [[1.0, -0.25], [-0.25, 0.75]]
    had = (np.ones((2, 2)) / math.sqrt(2)).tolist()
    had[1][1] *= -1
    df = DoubleFactorization.from_symmetric(0.3, [[[1.0, 0.0], [0.0, 1.0]], had], [g1, g2])
    be = double_factorized_encoding(df)
    hamiltonian = 0.3 * np.eye(2) + np.array(g1) + np.array(had) @ np.array(g2) @ np.array(had).T
    expected_alpha = 0.3 + sum(abs(v) for v in np.linalg.eigvalsh(np.array(g1))) + sum(
        abs(v) for v in np.linalg.eigvalsh(np.array(g2))
    )
    block, leakage = unitary_block(be)
    error = float(np.abs(block * be.alpha - hamiltonian).max())
    cross = cross_columns(be.operation, range(2))
    report.case(
        "df-encoding-rank2-scalar-d2",
        paths=["originir-ext+to_matrix", "reference", "rir-pysparq", "adapter-pysparq"],
        parameters={"rank": 2, "scalar": 0.3, "alpha": be.alpha},
        metrics={
            "max_error": error,
            "leakage": leakage,
            "cross_deviation": cross,
            "alpha_minus_closed_form": be.alpha - expected_alpha,
        },
        criterion="块*α == scalar·I + Σ U_r G_r U_rᵀ（max_error < 1e-9），α == |scalar|+Σ‖g_r‖₁",
        passed=error < EXACT and cross < EXACT and abs(be.alpha - expected_alpha) < 1e-9,
    )
    # 4x4 三秩随机旋转（直接谱构造）
    dim, rank = 4, 3
    rotations = [_random_orthonormal(dim, rng) for _ in range(rank)]
    spectra = [tuple(rng.uniform(-1.5, 1.5, size=dim).tolist()) for _ in range(rank)]
    df4 = DoubleFactorization(0.0, rotations, spectra)
    be4 = double_factorized_encoding(df4)
    hamiltonian4 = sum(
        np.array(u) @ np.diag(g) @ np.array(u).T for u, g in zip(rotations, spectra, strict=True)
    )
    alpha4 = sum(sum(abs(v) for v in g) for g in spectra)
    block4, leakage4 = unitary_block(be4)
    error4 = float(np.abs(block4 * be4.alpha - hamiltonian4).max())
    report.case(
        "df-encoding-rank3-d4",
        paths=["originir-ext+to_matrix", "reference"],
        parameters={"rank": rank, "dim": dim, "alpha": be4.alpha},
        metrics={"max_error": error4, "leakage": leakage4, "alpha": be4.alpha},
        criterion="块*α == Σ U_r diag(g_r) U_rᵀ（max_error < 1e-9）",
        passed=error4 < EXACT and abs(be4.alpha - alpha4) < 1e-9,
    )
    # 8x8 双秩：reference + rir 逐列
    dim8 = 8
    rotations8 = [_random_orthonormal(dim8, rng) for _ in range(2)]
    spectra8 = [tuple(rng.uniform(-1.0, 1.0, size=dim8).tolist()) for _ in range(2)]
    df8 = DoubleFactorization(0.7, rotations8, spectra8)
    be8 = double_factorized_encoding(df8)
    hamiltonian8 = 0.7 * np.eye(dim8) + sum(
        np.array(u) @ np.diag(g) @ np.array(u).T for u, g in zip(rotations8, spectra8, strict=True)
    )
    block8, leakage8 = block_via_reference(be8.operation, be8.alpha, 3)
    error8 = float(np.abs(block8 * be8.alpha - hamiltonian8).max())
    cross8 = cross_columns(be8.operation, sampled_columns(3), paths=("rir",))
    report.case(
        "df-encoding-rank2-d8",
        paths=["reference", "rir-pysparq"],
        parameters={
            "rank": 2,
            "dim": dim8,
            "alpha": be8.alpha,
            "note": "寄存器规模按预算走 pysparq 路径；交叉列抽样 8 列（全列由 reference 覆盖）",
        },
        metrics={"max_error": error8, "leakage": leakage8, "cross_deviation": cross8},
        criterion="块*α == H（max_error < 1e-9），两后端一致",
        passed=error8 < EXACT and cross8 < EXACT,
    )


def verify_thc(report):
    """THC 块编码：块*α == Σ ζ_{μν} L_μ L_ν†；α 与独立 Pauli-l1 手算口径对拍。"""
    leaf0 = [[0.6, 0.2], [0.1, -0.5]]
    leaf1 = [[0.3, -0.4], [0.2, 0.7]]
    zeta = [[0.8, 0.15], [0.15, -0.5]]
    thc = THCDecomposition(zeta, [leaf0, leaf1])
    be = thc_encoding(thc)
    leaves = [np.array(leaf0), np.array(leaf1)]
    hamiltonian = sum(
        zeta[mu][nu] * leaves[mu] @ leaves[nu].conj().T for mu in range(2) for nu in range(2)
    )
    expected_alpha = sum(
        abs(zeta[mu][nu]) * pauli_l1(leaves[mu]) * pauli_l1(leaves[nu])
        for mu in range(2)
        for nu in range(2)
    )
    block, leakage = unitary_block(be)
    error = float(np.abs(block * be.alpha - hamiltonian).max())
    cross = cross_columns(be.operation, range(2))
    report.case(
        "thc-encoding-leaves2-d2",
        paths=["originir-ext+to_matrix", "reference", "rir-pysparq", "adapter-pysparq"],
        parameters={"leaves": 2, "alpha": be.alpha},
        metrics={
            "max_error": error,
            "leakage": leakage,
            "cross_deviation": cross,
            "alpha_minus_independent": be.alpha - expected_alpha,
        },
        criterion="块*α == ΣζL_μL_ν†（max_error < 1e-9），α == Σ|ζ|α_μα_ν（独立 Pauli-l1）",
        passed=error < EXACT and cross < EXACT and abs(be.alpha - expected_alpha) < 1e-9,
    )
    # 4x4 三叶（对角叶 + 全耦合 ζ）：reference + rir。对角叶的 Pauli 项少、信号位宽小，
    # 电路规模可控；稠密叶的一般性由 leaves2-d2 的全幺正案例覆盖。
    rng = np.random.default_rng(29)
    leaves4 = [np.diag(rng.uniform(-0.9, 0.9, size=4)).tolist() for _ in range(3)]
    zeta4 = [[0.5, 0.1, -0.2], [0.1, 0.7, 0.05], [-0.2, 0.05, 0.6]]
    thc4 = THCDecomposition(zeta4, leaves4)
    be4 = thc_encoding(thc4)
    la = [np.array(leaf) for leaf in leaves4]
    hamiltonian4 = sum(
        zeta4[mu][nu] * la[mu] @ la[nu].conj().T for mu in range(3) for nu in range(3)
    )
    block4, leakage4 = block_via_reference(be4.operation, be4.alpha, 2)
    error4 = float(np.abs(block4 * be4.alpha - hamiltonian4).max())
    cross4 = cross_columns(be4.operation, range(4))
    report.case(
        "thc-encoding-leaves3-d4",
        paths=["reference", "rir-pysparq", "adapter-pysparq"],
        parameters={"leaves": 3, "dim": 4, "alpha": be4.alpha, "leaf_form": "diagonal"},
        metrics={"max_error": error4, "leakage": leakage4, "cross_deviation": cross4},
        criterion="块*α == ΣζL_μL_ν†（max_error < 1e-9），三后端一致",
        passed=error4 < EXACT and cross4 < EXACT,
    )


# ---------------------------------------------------------------------------
# 截断 Taylor 块编码（taylor-block-encoding 页面；由本组 lcu/product 组合而成）
# ---------------------------------------------------------------------------


def verify_taylor_block_encoding(report):
    """taylor_hamiltonian：块 == 截断级数/α（实现误差）与 e^{-iHt}/α（方法误差）分离报告。"""
    from pyqecclang.algorithms.hamiltonian import taylor_hamiltonian  # 页面归属本组

    hamiltonian = np.array([[1.0, 0.4], [0.4, -0.6]], dtype=complex)
    source = matrix_pauli_encoding(hamiltonian)
    time = 0.7
    eigenvalues, vectors = np.linalg.eigh(hamiltonian)
    exact = vectors @ np.diag(np.exp(-1j * time * eigenvalues)) @ vectors.conj().T
    for degree in (1, 2, 3, 4):
        be = taylor_hamiltonian(source, time, degree=degree)
        # 块的语义：截断级数 Σ(-itH)^k/k! 除以 LCU 归一化 α = Σ(α_source·t)^k/k!
        series = sum(
            (-1j * time) ** k / math.factorial(k) * np.linalg.matrix_power(hamiltonian, k)
            for k in range(degree + 1)
        )
        expected_alpha = sum((source.alpha * time) ** k / math.factorial(k) for k in range(degree + 1))
        block, leakage = block_via_reference(be.operation, be.alpha, 1)
        impl_error = float(np.abs(block - series / expected_alpha).max())
        method_error = float(np.abs(block - exact / expected_alpha).max())
        cross = cross_columns(be.operation, range(2))
        report.case(
            f"taylor-block-encoding-d{degree}",
            paths=["reference", "rir-pysparq", "adapter-pysparq"],
            parameters={"degree": degree, "time": time, "source_alpha": source.alpha},
            metrics={
                "impl_error": impl_error,
                "method_error": method_error,
                "cross_deviation": cross,
                "alpha_minus_series": be.alpha - expected_alpha,
                "leakage": leakage,
            },
            criterion="块 == Σ(-itH)^k/k!/α（impl_error < 1e-9）；方法误差（vs e^{-iHt}/α）随阶数下降（信息性）",
            passed=impl_error < EXACT and cross < EXACT and abs(be.alpha - expected_alpha) < 1e-12,
        )


# ---------------------------------------------------------------------------
# 与 pysparq 块编码模块的独立交叉验证
# ---------------------------------------------------------------------------


def verify_cross_tridiagonal(report):
    """同一三对角矩阵：pyqecclang 稀疏/Pauli 两路 × pysparq BlockEncodingTridiagonal。

    pysparq 侧归一化：dim≥4 时为 Frobenius 范数 norm_f（其 C++ 正确性测试域
    randint(2,5) 覆盖的规模）；dim=2 时 (0,0) 块实测退化为 A/(|α|+2|β|)
    （1 位主寄存器上加一/减一都触发溢出分支，anc==0 角块失去 Frobenius 归一，
    β=0 时无移位分支、仍等于 A/norm_f）。该边界行为作为实证语义记录，
    不是断言 pysparq 构造文档的标称归一化。
    """
    fmt = FixedFormat(3, 1, signed=True)
    for dim, alpha, beta, fmt_case in (
        (2, 1.5, -0.5, fmt),
        (4, 1.5, 0.5, fmt),
        (4, 2.0, -1.0, FixedFormat(4, 1, signed=True)),  # amax=2.0 超出 3 位有符号值域
        (8, 1.5, -0.5, fmt),
        (16, 1.25, 0.75, FixedFormat(4, 2, signed=True)),  # 1.25/0.75 需 0.25 步长精确表示
    ):
        n_bits = (dim - 1).bit_length()
        amax = max(abs(alpha), abs(beta))
        matrix = tridiagonal(dim, alpha, beta)
        rows = [[0, 1], [0, 1]] if dim == 2 else tridiagonal_rows(dim)
        be = sparse_be_for(matrix, fmt_case, amax, rows, signed=True)
        qecc_block, _ = block_via_reference(be.operation, be.alpha, n_bits)
        qecc_matrix = qecc_block * be.alpha
        norm_f = math.sqrt(dim * alpha**2 + 2 * (dim - 1) * beta**2)
        # pysparq 侧有效归一化（dim=2 且 β≠0 时的实测退化，见函数 docstring）
        ps_norm = norm_f if dim > 2 or beta == 0 else abs(alpha) + 2 * abs(beta)
        ps_matrix = ps_tridiagonal_block(alpha, beta, n_bits) * ps_norm
        qecc_err = float(np.abs(qecc_matrix - matrix).max())
        ps_err = float(np.abs(ps_matrix - matrix).max())
        cross_err = float(np.abs(qecc_matrix - ps_matrix).max())
        metrics = {
            "qecc_sparse_error": qecc_err,
            "psparq_error": ps_err,
            "cross_deviation": cross_err,
            "qecc_alpha": be.alpha,
            "psparq_norm_f": norm_f,
            "psparq_alpha_effective": ps_norm,
        }
        parameters = {"dim": dim, "alpha_diag": alpha, "beta_offdiag": beta}
        if ps_norm != norm_f:
            parameters["note"] = (
                "pysparq 侧 dim=2 边界行为：(0,0) 块归一化实测为 |α|+2|β| 而非 Frobenius 范数"
            )
        passed = qecc_err < EXACT and ps_err < EXACT and cross_err < EXACT
        if dim <= 8:
            # Pauli 路线：同一矩阵的第二条 pyqecclang 编码路径
            pauli_be = matrix_pauli_encoding(matrix)
            pauli_block, _ = block_via_reference(pauli_be.operation, pauli_be.alpha, n_bits)
            pauli_err = float(np.abs(pauli_block * pauli_be.alpha - matrix).max())
            cross_pauli = float(np.abs(pauli_block * pauli_be.alpha - ps_matrix).max())
            metrics["qecc_pauli_error"] = pauli_err
            metrics["cross_pauli_psparq"] = cross_pauli
            metrics["qecc_pauli_alpha"] = pauli_be.alpha
            passed = passed and pauli_err < EXACT and cross_pauli < EXACT
        report.case(
            f"cross-tridiagonal-d{dim}-a{alpha}-b{beta}",
            paths=["reference", "pysparq.BlockEncodingTridiagonal"],
            parameters=parameters,
            metrics=metrics,
            criterion="两实现有效块各乘自身归一化后均 == A 且互相一致（< 1e-9）",
            passed=passed,
        )


def verify_cross_qram_block_encoding(report):
    """pysparq BlockEncodingViaQRAM × pyqecclang 稀疏 BE：三对角与非三对角稀疏矩阵。"""
    fmt = FixedFormat(3, 1, signed=True)
    families = []
    families.append(("tridiagonal-d4", tridiagonal(4, 1.5, -0.5), tridiagonal_rows(4)))
    # 非三对角：对匹配稀疏图（每列一个非对角伙伴），s=2
    matched = np.eye(4) * 1.5
    for i, j in ((0, 2), (1, 3)):
        matched[i, j] = matched[j, i] = -0.5
    matched_rows = [[0, 2], [1, 3], [0, 2], [1, 3]]
    families.append(("matched-pairs-d4", matched, matched_rows))
    for name, matrix, rows in families:
        dim = len(matrix)
        n_bits = (dim - 1).bit_length()
        sparsity = len(rows[0])
        amax = float(np.abs(matrix).max())
        scaled = matrix / np.linalg.norm(matrix, "fro")  # pysparq QRAM 侧按 Frobenius 归一化入表
        be = sparse_be_for(matrix, fmt, amax, rows, signed=True)
        qecc_block, _ = block_via_reference(be.operation, be.alpha, n_bits)
        qecc_matrix = qecc_block * be.alpha
        ps_matrix = ps_qram_block(scaled, n_bits) * np.linalg.norm(matrix, "fro")
        qecc_err = float(np.abs(qecc_matrix - matrix).max())
        ps_err = float(np.abs(ps_matrix - matrix).max())
        cross_err = float(np.abs(qecc_matrix - ps_matrix).max())
        report.case(
            f"cross-qram-be-{name}",
            paths=["reference", "pysparq.BlockEncodingViaQRAM"],
            parameters={
                "dim": dim,
                "sparsity": sparsity,
                "qecc_alpha": be.alpha,
                "psparq_config": "data_size=50,rational=51,exponent=15（C++ 判据配置）",
            },
            metrics={
                "qecc_sparse_error": qecc_err,
                "psparq_error": ps_err,
                "cross_deviation": cross_err,
                "psparq_tolerance": 5e-3,
            },
            criterion="pyqecclang 侧 < 1e-9；pysparq 侧定点量化 ≤ 5e-3（C++ 容差 2^-15 量级）",
            passed=qecc_err < EXACT and ps_err < 5e-3 and cross_err < 5e-3,
        )


def run():
    report = Report(
        "blockencoding",
        "块编码组合代数、PREPARE–SELECT、稀疏与低秩块编码的角块对拍（块 == A/α）"
        "与 pysparq 自带块编码的独立交叉验证；对角/稀疏多矩阵多规模双路径覆盖。",
    )
    verify_diagonal_be_unitary(report)
    verify_diagonal_be_wide(report)
    verify_diagonal_be_public_df(report)
    verify_pauli_word_and_embeddings(report)
    verify_matrix_pauli_encoding(report)
    verify_be_algebra_combinators(report)
    verify_gate_prepare_distribution(report)
    verify_select_pauli(report)
    verify_lcu_prepare_select_block(report)
    verify_qram_prepare_quantization(report)
    verify_alias_prepare(report)
    verify_alias_prepare_select_block(report)
    verify_abstract_prepare_bind(report)
    verify_sparse_rotation_helpers(report)
    verify_sparse_boolean_helpers(report)
    verify_sparse_lookup_helpers(report)
    verify_sparse_access_layer(report)
    verify_sparse_be_signed(report)
    verify_sparse_be_unsigned(report)
    verify_sparse_be_qram_access(report)
    verify_chebyshev_walk(report)
    verify_select_swap(report)
    verify_diagonalize_symmetric(report)
    verify_double_factorization(report)
    verify_thc(report)
    verify_taylor_block_encoding(report)
    verify_cross_tridiagonal(report)
    verify_cross_qram_block_encoding(report)
    report.write()
    return report


if __name__ == "__main__":
    run()
