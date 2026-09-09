"""QPDE/QODE 教程：不同输入范式、分批绑定和可替换的生成函数。

运行：PYTHONPATH=src python examples/ode_input_models.py
可选 --native-parse 使用真实 uniqc 解析导出；--native-bindings 对拍小型 QRAM 绑定。
这些检查不认证求解精度。
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
    identity,
    loads,
    run_pysparq,
    scale,
    unresolved,
    zero,
)
from pyqecclang.algorithms.block_encoding import lcu
from pyqecclang.algorithms.carleman import PolynomialODE, carleman_qode
from pyqecclang.algorithms.hamiltonian import taylor_hamiltonian
from pyqecclang.algorithms.lchs import QuadraturePlan, lchs_qode
from pyqecclang.algorithms.ode import linear_qode
from pyqecclang.algorithms.ode_models import HermitianParts, LinearODE
from pyqecclang.algorithms.oracles import (
    abstract_block_encoding,
    abstract_database,
    abstract_sparse_access,
    abstract_state_prep,
    diagonal_block_encoding,
    gate_database,
    gate_state_prep,
    qram_database,
    qram_state_angles,
    qram_state_prep,
    sparse_entry,
    sparse_location_gate,
    sparse_location_qram,
)
from pyqecclang.algorithms.pde import DiscretePDE, PDEInput, make_qpde, qpde_solver
from pyqecclang.algorithms.schrodingerization import SchrodingerPlan
from pyqecclang.algorithms.sparse import real_symmetric_sparse_encoding
from pyqecclang.applications.qham import (
    Discretization,
    Field,
    Grid,
    PolynomialPDE,
    structured_fd_bindings,
)
from pyqecclang.applications.qham.stencils import derivative_encoding


def save_case(root, name, state, bindings=None, memory=None, *, notes=None, native_parse=False):
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
    if native_parse:
        from uniqc.compile.originir.originir_base_parser import OriginIR_BaseParser

        OriginIR_BaseParser().parse(basis)
    record = {
        "case": name,
        "open_slots": [r.name for r in unresolved(opened)],
        "modules": len(closed.modules),
        "resources": [r.name for r in closed.main.resources],
        "target_width": state.width,
        "signal_width": state.signal_qubits,
        "native_parse": "passed" if native_parse else "not_requested",
        "notes": notes or {},
        "correctness": "solver accuracy and success probability pending",
    }
    (folder / "report.json").write_text(
        json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(name, "modules", len(closed.modules), "open slots", len(record["open_slots"]), flush=True)
    return record


def verify_bindings(root):
    """真实 PySparQ 比较角数据库案例两种绑定的完整复幅度。"""
    records = []
    for stem in ("given_xor",):
        states = []
        for implementation in ("gate", "qram"):
            directory = root / f"{stem}_{implementation}_lchs"
            program = loads((directory / "closed.rir.json").read_text(encoding="utf-8"))
            memory = json.loads((directory / "memory.json").read_text(encoding="utf-8"))
            states.append(run_pysparq(program, memory).amplitudes)
        keys = states[0].keys() | states[1].keys()
        error = max(abs(states[0].get(k, 0) - states[1].get(k, 0)) for k in keys)
        assert error < 1e-10, (stem, error)
        records.append(
            {
                "case": stem,
                "gate_vs_qram_max_amplitude_difference": error,
                "backend": "real pysparq",
                "scope": "same finite circuit construction, not ODE convergence",
            }
        )
        print(stem, "gate/QRAM amplitude difference", error, flush=True)
    (root / "binding-validation.json").write_text(json.dumps(records, indent=2), encoding="utf-8")


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


def shifted_solver(qode, recovery):
    """给提升后的线性系统施加移位；范数恢复记录放在宿主报告中。"""

    def solve(generator, initial, time):
        mu = generator.alpha
        shifted = lcu([(1, generator), (-mu, identity(generator.width))])
        recovery.update(growth_shift=mu, log_amplitude_rescale=mu * time)
        return qode(shifted, initial, time)

    return solve


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", type=Path, default=Path("out/ode-input-models"))
    parser.add_argument("--native-parse", action="store_true", help="用真实 uniqc 解析全部导出")
    parser.add_argument("--native-bindings", action="store_true", help="用 PySparQ 对拍角数据库绑定")
    args = parser.parse_args()
    write = partial(save_case, args.output, native_parse=args.native_parse)
    records = []
    time = 0.05
    hamiltonian_function = partial(taylor_hamiltonian, degree=1)
    quadrature = QuadraturePlan.cauchy(cutoff=1, spacing=1.0)
    lchs = linear_qode("lchs", plan=quadrature, hamiltonian_function=hamiltonian_function)
    schrodinger = linear_qode(
        "schrodingerization",
        plan=SchrodingerPlan(auxiliary_width=2, period=8.0, selected_index=1),
        hamiltonian_function=hamiltonian_function,
    )

    # 1. 直接给 G 的 BE；G=-I。work=3 预留给后面的 QRAM 初态实现。
    initial = abstract_state_prep("Initial", 1, work_width=3)
    initial_gate = gate_state_prep([1, 1], work_width=3)
    generator = abstract_block_encoding("Generator", 1, 0, 1.0)
    state = lchs(generator, initial, time)
    records.append(
        write(
            "given_be_lchs",
            state,
            {
                "Generator": scale(-1, identity(1)).operation,
                "Initial": initial_gate.operation,
            },
        )
    )

    # 2. 已经给 A=L+iH 时直接用 parts 入口，避免先合并再拆开。
    dissipation = abstract_block_encoding("Dissipation", 1, 0, 1.0)
    hamiltonian = abstract_block_encoding("Hamiltonian", 1, 1, 1.0)
    state = lchs_qode(
        LinearODE(HermitianParts(dissipation, hamiltonian), initial),
        time,
        plan=quadrature,
        hamiltonian_function=hamiltonian_function,
    )
    records.append(
        write(
            "given_parts_lchs",
            state,
            {
                "Dissipation": identity(1).operation,
                "Hamiltonian": zero(1).operation,
                "Initial": initial_gate.operation,
            },
        )
    )

    # 3. 相同开放图，分别绑定 gate / QRAM 数据库及态制备。
    # 数据是旋转角的整数编码；A=diag(1, cos(pi/4))，G=-A。
    angles = abstract_database("DiagonalAngles", 1, 2)
    a = diagonal_block_encoding(angles, alpha=1.0)
    state = lchs(scale(-1, a), initial, time)
    records.append(
        write(
            "given_xor_gate_lchs",
            state,
            {
                "DiagonalAngles": gate_database(1, 2, [0, 1]).operation,
                "Initial": initial_gate.operation,
            },
        )
    )
    records.append(
        write(
            "given_xor_qram_lchs",
            state,
            {
                "DiagonalAngles": Binding(
                    qram_database(1, 2).operation, {"table": "diagonal_angles"}
                ),
                "Initial": Binding(qram_state_prep(1, 2).operation, {"angles": "initial_angles"}),
            },
            memory={
                "diagonal_angles": [0, 1],
                "initial_angles": qram_state_angles([1, 1], 2),
            },
        )
    )

    # 4. A=[[1,-.5],[-.5,1]]：实对称、非负对角且 A>=0。
    fmt = FixedFormat(4, 1)
    access = abstract_sparse_access("SparseA", 1, fmt.width, sparsity=2)
    a = real_symmetric_sparse_encoding(access, fmt, 1.0, diagonal_nonnegative=True)
    state = lchs(scale(-1, a), initial, time)
    values = [fmt.encode(v) for v in (1, -0.5, -0.5, 1)]
    records.append(
        write(
            "given_sparse_gate_lchs",
            state,
            {
                "SparseA_position": sparse_location_gate(1, [(0, 1), (0, 1)]),
                "SparseA_entry": sparse_entry(gate_database(2, fmt.width, values), 1),
                "Initial": initial_gate.operation,
            },
            notes={"sparsity": 2, "amax": 1.0, "alpha_A": a.alpha},
        )
    )
    records.append(
        write(
            "given_sparse_qram_lchs",
            state,
            {
                "SparseA_position": Binding(
                    sparse_location_qram(1),
                    {
                        "forward": "positions",
                        "inverse": "inverse_positions",
                    },
                ),
                "SparseA_entry": Binding(
                    sparse_entry(qram_database(2, fmt.width), 1),
                    {
                        "db__table": "entries",
                    },
                ),
                "Initial": initial_gate.operation,
            },
            memory={
                "positions": [0, 0, 1, 1],
                "inverse_positions": [0, 0, 1, 1],
                "entries": values,
            },
        )
    )

    # 5. 周期四点热方程，结构化移位 BE 直接作为空间离散化结果。
    grid = Grid(("x",), (4,), (1.0,), boundary="periodic")
    heat = DiscretePDE(
        scale(0.1, derivative_encoding(grid, (("x", 2),))),
        gate_state_prep([1, 0, 0, 0]),
        label="periodic_heat",
    )
    for method, solver in (("lchs", lchs), ("schrodingerization", schrodinger)):
        records.append(write("heat_" + method, make_qpde(solver)(heat, time)))

    # 6. 普通 PDE -> F_p 端口 -> Carleman -> 可替换线性 solver。
    u = Field("u")
    pde = PolynomialPDE.from_equations({"u": 0.1 * u.d("x", 2) - u * u.d("x")})
    space = Discretization(pde, grid)
    ports = structured_fd_bindings(space, space.encode_fields({"u": [0.1, 0.2, 0, -0.1]}))
    concrete = polynomial_from_bindings(ports)
    coefficient_slots, bindings = [], {}
    for order, coefficient in concrete.coefficients:
        name = "BurgersF" + str(order)
        slot = abstract_block_encoding(
            name, coefficient.width, coefficient.signal_qubits, coefficient.alpha
        )
        coefficient_slots.append((order, slot))
        bindings[name] = coefficient.operation
    initial_slot = abstract_state_prep(
        "BurgersInitial", concrete.width, concrete.initial.work_width
    )
    bindings["BurgersInitial"] = concrete.initial.operation
    problem = PolynomialODE(
        concrete.width, tuple(coefficient_slots), initial_slot, concrete.initial_norm
    )
    recovery = {}
    for method, solver in (
        ("schrodingerization", schrodinger),
        ("shifted_lchs", shifted_solver(lchs, recovery)),
    ):
        nonlinear_solver = partial(carleman_qode, linear_solver=solver, cutoff=2)
        state = qpde_solver(nonlinear_solver)(PDEInput(problem, "burgers"), time)
        records.append(
            write(
                "burgers_carleman_" + method,
                state,
                bindings,
                notes={
                    "initial_norm": concrete.initial_norm,
                    "cutoff": 2,
                    **(recovery if method == "shifted_lchs" else {}),
                },
            )
        )

    # 7. 同一 QPDE，独立替换内部 Hamiltonian-function protocol 的配置。
    degree_two = linear_qode(
        "lchs", plan=quadrature, hamiltonian_function=partial(taylor_hamiltonian, degree=2)
    )
    records.append(write("heat_lchs_taylor2", make_qpde(degree_two)(heat, time)))
    (args.output / "index.json").write_text(
        json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    if args.native_bindings:
        verify_bindings(args.output)


if __name__ == "__main__":
    main()
