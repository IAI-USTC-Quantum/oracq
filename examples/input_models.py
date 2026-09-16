"""输入模型变体教程：四个求解器家族 × 结构化 / 谱嵌入 / QRAM 数据路径。

运行：PYTHONPATH=src python examples/input_models.py
产物只记录结构证据（开放槽位、资源、闭合与导出），求解精度属于后续验证项。
每个家族在同一问题实例下展示不同的输入模型：结构化门端口（基线）、
谱嵌入（Pauli 精确展开）、以及数据路径走开放角数据库/QRAM 资源的变体。
"""

import argparse
import json
from collections import defaultdict
from functools import partial
from pathlib import Path

from pyqecclang import (
    Binding,
    FixedFormat,
    bind,
    dumps,
    export_originir,
    export_toffoli_u3_cz,
    loads,
    scale,
    unresolved,
)
from pyqecclang.algorithms.block_encoding import lcu, matrix_pauli_encoding
from pyqecclang.algorithms.carleman import PolynomialODE, carleman_qode
from pyqecclang.algorithms.hamiltonian import taylor_hamiltonian
from pyqecclang.algorithms.lchs import QuadraturePlan
from pyqecclang.algorithms.ode import linear_qode
from pyqecclang.algorithms.oracles import (
    abstract_block_encoding,
    abstract_sparse_access,
    abstract_state_prep,
    gate_database,
    gate_state_prep,
    qram_database,
    qram_state_angles,
    qram_state_prep,
    sparse_entry,
    sparse_location_qram,
)
from pyqecclang.algorithms.pde import DiscretePDE, make_qpde
from pyqecclang.algorithms.prepare_select import lcu_prepare_select, qram_prepare
from pyqecclang.algorithms.sparse import real_symmetric_sparse_encoding
from pyqecclang.applications.qham import (
    Discretization,
    Field,
    Grid,
    Known,
    PolynomialPDE,
    QHAMBindings,
    QHAMPlan,
    gate_bindings,
    qham_input_model,
    qram_coefficient_encoding,
    qram_coefficient_memory,
    structured_fd_bindings,
)


def save_case(root, name, state, bindings=None, memory=None, *, notes=None):
    """同时留下开放/部分绑定/闭合产物；不内联 oracle 主体。"""
    folder = root / name
    folder.mkdir(parents=True, exist_ok=True)
    opened = state.operation.program()
    assert loads(dumps(opened)) == opened
    (folder / "open.rir.json").write_text(dumps(opened), encoding="utf-8")
    closed = opened
    for index, (slot, implementation) in enumerate((bindings or {}).items()):
        closed = bind(closed, {slot: implementation})
        if index == 0:
            (folder / "partial.rir.json").write_text(dumps(closed), encoding="utf-8")
    assert not unresolved(closed), unresolved(closed)
    assert loads(dumps(closed)) == closed
    (folder / "closed.rir.json").write_text(dumps(closed), encoding="utf-8")
    (folder / "modular.originir").write_text(export_originir(closed).text, encoding="utf-8")
    basis = export_toffoli_u3_cz(closed).text
    (folder / "toffoli_u3_cz.originir").write_text(basis, encoding="utf-8")
    # 教程的小表统一存为数组，避免 JSON 对象的字符串键被误当作整数地址。
    memory = dict(memory or {})
    for resource in closed.main.resources:
        bank = memory[resource.name]
        if isinstance(bank, dict):
            memory[resource.name] = [
                bank.get(address, 0) for address in range(1 << resource.type.address_width)
            ]
    (folder / "memory.json").write_text(json.dumps(memory, indent=2), encoding="utf-8")
    record = {
        "case": name,
        "open_slots": [r.name for r in unresolved(opened)],
        "modules": len(closed.modules),
        "resources": [r.name for r in closed.main.resources],
        "target_width": state.width,
        "signal_width": state.signal_qubits,
        "notes": notes or {},
        "correctness": "solver accuracy and success probability pending",
    }
    (folder / "report.json").write_text(
        json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(name, "modules", len(closed.modules), "open slots", len(record["open_slots"]), flush=True)
    return record


def polynomial_from_bindings(bindings):
    """教程宿主适配器：按次数合并 PDE 多线性端口，交给 Carleman。"""
    grouped = defaultdict(list)
    for _, port in bindings.ports:
        grouped[port.arity].append((1, port.encoding))
    return PolynomialODE(
        bindings.state_width,
        tuple((order, lcu(terms)) for order, terms in sorted(grouped.items())),
        bindings.initial,
        bindings.initial_norm,
    )


def reopen_coefficients(problem, label):
    """把具体系数 BE 重新声明为开放槽位；返回 (PolynomialODE, 绑定)。"""
    coefficient_slots, bindings = [], {}
    for order, coefficient in problem.coefficients:
        name = label + "F" + str(order)
        slot = abstract_block_encoding(
            name, coefficient.width, coefficient.signal_qubits, coefficient.alpha
        )
        coefficient_slots.append((order, slot))
        bindings[name] = coefficient.operation
    initial_slot = abstract_state_prep(
        label + "Initial", problem.width, problem.initial.work_width
    )
    bindings[label + "Initial"] = problem.initial.operation
    return (
        PolynomialODE(problem.width, tuple(coefficient_slots), initial_slot, problem.initial_norm),
        bindings,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", type=Path, default=Path("out/input-models"))
    args = parser.parse_args()
    write = partial(save_case, args.output)
    records = []
    time = 0.05
    kernel = partial(taylor_hamiltonian, degree=1)
    quadrature = QuadraturePlan.cauchy(cutoff=1, spacing=1.0)
    lchs = linear_qode("lchs", plan=quadrature, hamiltonian_function=kernel)
    schrodinger = linear_qode("schrodingerization", hamiltonian_function=kernel)
    cbmd = linear_qode("cbmd", hamiltonian_function=kernel)

    # ---- QHAM 家族：强迫 Burgers（order 2, eta=-0.4），输入模型三种 ----------
    u = Field("u")
    pde = PolynomialPDE.from_equations(
        {"u": 0.1 * u.d("x", 2) - u * u.d("x") + Known("forcing")},
        label="forced_burgers",
    )
    plan = QHAMPlan(pde, order=2)
    grid = Grid(("x",), (4,), (1.0,), boundary="periodic")
    space = Discretization(pde, grid, {"forcing": [0.05, 0, -0.05, 0]})

    # 1. 结构化门端口 + 门初态：基线（与 examples/general_qham.py 相同输入）。
    initial = space.encode_fields({"u": [0.1, 0.2, 0, -0.1]})
    stencil_bindings = structured_fd_bindings(space, initial)
    model = qham_input_model(plan, stencil_bindings, eta=-0.4)
    records.append(write("qham_stencil_gate", model.solve(schrodinger, 0.01)))

    # 2. 谱嵌入端口：基础端口矩阵的精确 Pauli-LCU 展开 + 门初态。
    spectral_bindings = gate_bindings(space, initial)
    model = qham_input_model(plan, spectral_bindings, eta=-0.4)
    records.append(write("qham_spectral_ports", model.solve(schrodinger, 0.01)))

    # 3+4. 系数数据留在开放角数据库槽位、初态声明为槽位：同一份开放程序
    # 分别绑定 gate 数据库与 QRAM 数据库（运行时角表单独提供）。
    # QRAM 角表当前仅接受非负幅度，故这一对案例使用非负剖面。
    profile = [0.1, 0.2, 0.15, 0.05]
    qinitial = space.encode_fields({"u": profile})
    encoder = partial(qram_coefficient_encoding, angle_width=8)
    open_bindings = structured_fd_bindings(space, qinitial, coefficient_encoder=encoder)
    open_bindings = QHAMBindings(
        open_bindings.state_width,
        open_bindings.ports,
        abstract_state_prep("QhamInitial", open_bindings.state_width, 10),
        open_bindings.initial_norm,
    )
    model = qham_input_model(plan, open_bindings, eta=-0.4)
    state = model.solve(schrodinger, 0.01)
    coeff_slots = [
        r.name for r in unresolved(state.operation.program()) if r.name != "QhamInitial"
    ]
    assert len(coeff_slots) == 1, coeff_slots
    (db_slot,) = coeff_slots
    (forcing_monomial,) = (
        term.monomial for port in space.pde.ports for term in port.terms if term.monomial.known
    )
    angle_words = qram_coefficient_memory(space, forcing_monomial, angle_width=8)
    word_table = [angle_words.get(address, 0) for address in range(4)]
    records.append(
        write(
            "qham_coeff_gate",
            state,
            {
                db_slot: gate_database(2, 8, word_table).operation,
                "QhamInitial": gate_state_prep(profile, work_width=10).operation,
            },
        )
    )
    qram_memory = {
        "coeff_angles": word_table,
        "initial_angles": qram_state_angles(profile, 8),
    }
    qram_dict = {
        db_slot: Binding(qram_database(2, 8).operation, {"table": "coeff_angles"}),
        "QhamInitial": Binding(qram_state_prep(2, 8).operation, {"angles": "initial_angles"}),
    }
    records.append(write("qham_coeff_qram", state, qram_dict, memory=qram_memory))

    # ---- Carleman 家族：无粘性 Burgers，cutoff=2，输入模型三种 -------------
    pde_c = PolynomialPDE.from_equations({"u": 0.1 * u.d("x", 2) - u * u.d("x")})
    space_c = Discretization(pde_c, grid)
    initial_c = space_c.encode_fields({"u": profile})

    # 5. 结构化门端口（槽位重新抽象，同 ode_input_models 案例 6 的模式）。
    concrete = polynomial_from_bindings(structured_fd_bindings(space_c, initial_c))
    problem, rebind = reopen_coefficients(concrete, "Burgers")
    records.append(
        write(
            "carleman_stencil_gate",
            carleman_qode(problem, time, schrodinger, cutoff=2),
            rebind,
        )
    )

    # 6. 谱嵌入端口（Pauli-LCU 基础端口矩阵）。
    concrete = polynomial_from_bindings(gate_bindings(space_c, initial_c))
    problem, rebind = reopen_coefficients(concrete, "Spectral")
    records.append(
        write(
            "carleman_spectral_ports",
            carleman_qode(problem, time, schrodinger, cutoff=2),
            rebind,
        )
    )

    # 7. 空间变化粘性：系数走开放角数据库 + QRAM 初态。
    pde_nu = PolynomialPDE.from_equations(
        {"u": 0.1 * Known("nu") * u.d("x", 2) - u * u.d("x")}, label="burgers_nu"
    )
    space_nu = Discretization(pde_nu, grid, {"nu": [1.0, 0.8, 0.6, 1.2]})
    concrete = polynomial_from_bindings(
        structured_fd_bindings(space_nu, initial_c, coefficient_encoder=encoder)
    )
    problem, rebind = reopen_coefficients(concrete, "QramBurgers")
    # 初态槽位按 QRAM 制备的寄存器布局声明（gate 实现的 work 宽度为 0）。
    del rebind["QramBurgersInitial"]
    qinit_slot = abstract_state_prep("QramBurgersInitial", concrete.width, 10)
    problem = PolynomialODE(
        concrete.width, problem.coefficients, qinit_slot, concrete.initial_norm
    )
    extra = [
        r.name
        for _, coefficient in concrete.coefficients
        for r in unresolved(coefficient.operation.program())
    ]
    assert len(extra) == 1, extra
    (nu_slot,) = extra
    (nu_monomial,) = (
        term.monomial for port in space_nu.pde.ports for term in port.terms if term.monomial.known
    )
    nu_words = qram_coefficient_memory(space_nu, nu_monomial, angle_width=8)
    rebind[nu_slot] = Binding(qram_database(2, 8).operation, {"table": "nu_angles"})
    rebind["QramBurgersInitial"] = Binding(
        qram_state_prep(2, 8).operation, {"angles": "initial_angles"}
    )
    state = carleman_qode(problem, time, schrodinger, cutoff=2)
    records.append(
        write(
            "carleman_coeff_qram",
            state,
            rebind,
            memory={
                "nu_angles": [nu_words.get(address, 0) for address in range(4)],
                "initial_angles": qram_state_angles(profile, 8),
            },
        )
    )

    # ---- LCHS / CBMD：给定生成元的谱嵌入（与 given_sparse 案例同矩阵） ------
    a_matrix = [[1, -0.5], [-0.5, 1]]
    spectral_be = matrix_pauli_encoding(a_matrix)
    initial_slot = abstract_state_prep("SpectralInitial", 1, work_width=3)
    initial_gate = gate_state_prep([1, 1], work_width=3)
    gen_slot = abstract_block_encoding(
        "SpectralGenerator", 1, spectral_be.signal_qubits, spectral_be.alpha
    )
    for family, solver in (("lchs", lchs), ("cbmd", cbmd)):
        state = solver(scale(-1, gen_slot), initial_slot, time)
        records.append(
            write(
                family + "_given_spectral",
                state,
                {
                    "SpectralGenerator": scale(-1, spectral_be).operation,
                    "SpectralInitial": initial_gate.operation,
                },
            )
        )

    # 谱嵌入 + QRAM PREPARE：Pauli 系数从 QRAM 角度表加载。A = I - 0.5 X。
    terms = [(1.0, "I"), (-0.5, "X")]
    qprep = qram_prepare([1.0, -0.5], angle_width=8)
    qbe = lcu_prepare_select(terms, prepare=qprep.preparation)
    qgen_slot = abstract_block_encoding(
        "QramSpectralGenerator", 1, qbe.signal_qubits, qbe.alpha
    )
    qgen_op = scale(-1, qbe).operation
    (angles_formal,) = [r.name for r in qgen_op.program().main.resources]
    for family, solver in (("lchs", lchs), ("cbmd", cbmd)):
        state = solver(scale(-1, qgen_slot), initial_slot, time)
        records.append(
            write(
                family + "_spectral_qram_prepare",
                state,
                {
                    "QramSpectralGenerator": Binding(
                        qgen_op, {angles_formal: "spectral_angles"}
                    ),
                    "SpectralInitial": initial_gate.operation,
                },
                memory={"spectral_angles": qprep.memory["angles"]},
            )
        )

    # ---- LCHS / CBMD：QRAM 网格离散化 ---------------------------------------
    # 四点周期环 Laplacian 的稀疏访问（位置 + 条目表来自 QRAM）+ QRAM 初态。
    fmt = FixedFormat(4, 1)
    access = abstract_sparse_access("RingStencil", 2, fmt.width, sparsity=3)
    laplacian = real_symmetric_sparse_encoding(access, fmt, 2.0, diagonal_nonnegative=True)
    heat = DiscretePDE(
        scale(-0.1, laplacian),
        abstract_state_prep("GridInitial", 2, 10),
        label="ring_heat",
    )
    ring = [
        [2 if row == col else (-1 if (row - col) % 4 in (1, 3) else 0) for col in range(4)]
        for row in range(4)
    ]
    slots = {col: [row for row in range(4) if ring[row][col]] for col in range(4)}
    forward, inverse = [0] * 16, [0] * 16
    for col in range(4):
        for slot, row in enumerate(slots[col]):
            forward[col + slot * 4] = row
            inverse[col + row * 4] = slot
    entries = [fmt.encode(ring[row][col]) for col in range(4) for row in range(4)]
    grid_bindings = {
        "RingStencil_position": Binding(
            sparse_location_qram(2),
            {"forward": "positions", "inverse": "inverse_positions"},
        ),
        "RingStencil_entry": Binding(
            sparse_entry(qram_database(4, fmt.width), 2),
            {"db__table": "entries"},
        ),
        "GridInitial": Binding(qram_state_prep(2, 8).operation, {"angles": "grid_angles"}),
    }
    grid_memory = {
        "positions": forward,
        "inverse_positions": inverse,
        "entries": entries,
        "grid_angles": qram_state_angles([1, 0, 0, 0], 8),
    }
    for family, solver in (("lchs", lchs), ("cbmd", cbmd)):
        state = make_qpde(solver)(heat, time)
        records.append(
            write(family + "_qram_grid", state, grid_bindings, memory=grid_memory)
        )

    # ---- CBMD：QHAM 开放系数输入经耗散移位（输入模型跨家族复用） ------------
    model = qham_input_model(plan, open_bindings, eta=-0.4)
    state = model.dissipative_shift().solve(cbmd, 0.01)
    records.append(
        write("cbmd_qham_coeff_qram", state, qram_dict, memory=qram_memory)
    )

    (args.output / "index.json").write_text(
        json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
