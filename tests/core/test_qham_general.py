"""一般 QHAM 的数学身份、结构预算、矩形映射和初值尺度。"""

import math
import random
import unittest
from functools import partial

from oracq import ValidationError, bind, dumps, loads, simulate, unresolved
from oracq.algorithms.input_model.oracles import gate_database
from oracq.algorithms.input_model.qham import embed_rectangular, place_port
from oracq.applications.qham import (
    Block,
    Discretization,
    Field,
    Grid,
    Known,
    PolynomialPDE,
    QHAMBindings,
    QHAMPlan,
    gate_bindings,
    open_qham_input,
    qham_input_model,
    qram_coefficient_encoding,
    qram_coefficient_memory,
    structured_fd_bindings,
    taylor_qode,
)


class QhamGeneralTests(unittest.TestCase):
    def test_rules_roundtrip_and_unsupported_nonlinearity(self):
        u, v = Field("u"), Field("v")
        pde = PolynomialPDE.from_equations(
            {"u": -0.1 * (u * u).d("x") + Known("a") * v, "v": v.d("x", 3) + u * v * v},
            label="coupled",
        )
        self.assertEqual(PolynomialPDE.loads(pde.dumps()), pde)
        self.assertEqual(pde.degree, 3)
        with self.assertRaises(ValidationError):
            _ = 1 / u
        with self.assertRaises(ValidationError):
            _ = u**0.5

    def test_lazy_closure_and_rank_unrank(self):
        u = Field("u")
        for degree in (2, 3, 4):
            pde = PolynomialPDE.from_equations({"u": u + u**degree + Known("f")})
            for order in range(4):
                plan = QHAMPlan(pde, order)
                self.assertEqual(QHAMPlan.loads(plan.dumps()), plan)
                for block in plan.blocks():
                    for term in plan.row_terms(block):
                        self.assertTrue(plan.contains(term.column))
                        if term.arity >= 2:
                            before = sum(block.orders) if block.kind == "tensor" else order
                            self.assertLess(sum(term.column.orders), before)
                    offset = plan.offset(block, 2)
                    self.assertEqual(plan.locate(offset, 2), (block, 0))
                    self.assertEqual(
                        plan.locate(offset + 2**block.rank - 1, 2), (block, 2**block.rank - 1)
                    )
        huge = QHAMPlan(PolynomialPDE.from_equations({"u": u * u}), 20)
        self.assertEqual(huge.block_count, 1 << 21)
        self.assertLess(len(huge.dumps()), 3000)
        with self.assertRaises(ValidationError):
            tuple(huge.blocks(max_blocks=16))

    def test_quadratic_dimension_and_nontrivial_eta_weight(self):
        u = Field("u")
        pde = PolynomialPDE.from_equations({"u": u * u})
        for order in range(5):
            plan = QHAMPlan(pde, order)
            self.assertEqual(plan.block_count, 2 ** (order + 1))
            self.assertEqual(plan.raw_dimension(4), 5 ** (order + 1) + 3)
        row = QHAMPlan(pde, 2).row_terms(Block("tensor", (2,)))
        weights = {term.column.orders: term.weight.evaluate(-0.4) for term in row}
        self.assertAlmostEqual(weights[(0, 0)], 0.24)
        self.assertAlmostEqual(weights[(0, 1)], 0.4)
        self.assertAlmostEqual(weights[(1, 0)], 0.4)

    def test_chain_rule_identity_general_cases(self):
        u, v = Field("u"), Field("v")
        cases = [
            ({"u": 0.1 * u.d("x", 2) - u * u.d("x") + Known("f")}, 4, 3),
            ({"u": -u.d("x", 3) - 6 * u * u.d("x")}, 4, 2),
            ({"u": -0.2 * u + 0.1 * u**3 + Known("f")}, 2, 2),
            (
                {"u": -0.1 * u + 0.2 * u * v + Known("f"), "v": 0.3 * u - 0.2 * v + 0.05 * v * v},
                2,
                2,
            ),
            ({"u": -(Known("a") * u * u).d("x") + 0.1 * u.d("x", 2)}, 4, 2),
        ]
        for equations, size, order in cases:
            with self.subTest(equations=tuple(equations), size=size, order=order):
                pde = PolynomialPDE.from_equations(equations)
                disc = Discretization(
                    pde,
                    Grid(("x",), (size,), (1.0,)),
                    {"f": [0.03] * size, "a": [1 + 0.1 * i for i in range(size)]},
                )
                plan = QHAMPlan(pde, order)
                rng = random.Random(11)
                values = [
                    [rng.uniform(-0.2, 0.2) for _ in range(disc.dimension)]
                    for _ in range(order + 1)
                ]
                for eta in (-1.0, -0.4, 0.2):
                    lifted = disc.lift(plan, values)
                    matrix = disc.matrix(plan, eta, max_dimension=2048)
                    direct = disc.chain_rule(plan, values, eta)
                    actual = [
                        sum(a * b for a, b in zip(row, lifted, strict=True)) for row in matrix
                    ]
                    self.assertLess(
                        max(abs(a - b) for a, b in zip(actual, direct, strict=True)), 1e-10
                    )

    def test_structured_initial_preserves_relative_tensor_norms(self):
        u = Field("u")
        pde = PolynomialPDE.from_equations({"u": -0.2 * u + 0.1 * u * u + Known("f")})
        plan = QHAMPlan(pde, 2)
        disc = Discretization(pde, Grid(("x",), (2,), (1.0,)), {"f": [0.1, 0.2]})
        initial = [0.2, 0.1]
        bindings = gate_bindings(disc, initial)
        model = qham_input_model(plan, bindings, eta=-0.4)
        actual = simulate(model.initial.operation.program())
        expected = disc.lift(plan, [initial, [0, 0], [0, 0]])
        norm = math.sqrt(sum(abs(v) ** 2 for v in expected))
        self.assertAlmostEqual(math.exp(model.log_initial_norm), norm)
        for row, value in enumerate(expected):
            self.assertAlmostEqual(actual.amplitudes.get((row, 0), 0), value / norm, places=10)
        self.assertTrue(all(work == 0 for _, work in actual.amplitudes))

    def test_rectangular_forcing_and_contraction_do_not_spill(self):
        u = Field("u")
        pde = PolynomialPDE.from_equations({"u": 0.1 * u * u + Known("f")})
        disc = Discretization(pde, Grid(("x",), (2,), (1.0,)), {"f": [0.2, -0.1]})
        bindings = gate_bindings(disc, [0.2, 0])
        for key in ("F", next(k for k, _ in bindings.ports if k.startswith("B_"))):
            port = dict(bindings.ports)[key]
            source_rank = port.arity
            placed = place_port(port, 1, source_rank, 0)
            # 非对齐、窄输出窗口，其他行/列必须为零。
            encoded = embed_rectangular(placed, 3, 1, 4, source_rank, 1)
            for column in range(8):
                actual = simulate(encoded.operation.program(), initial={"target": column})
                for row in range(8):
                    expected = (
                        disc.entry(key, row - 1, column - 4)
                        if 1 <= row < 3 and 4 <= column < 4 + 2**source_rank
                        else 0
                    )
                    self.assertAlmostEqual(
                        actual.amplitudes.get((row, 0), 0) * encoded.alpha, expected, places=10
                    )

    def test_open_qode_inputs_keep_operator_and_initial_oracles(self):
        u = Field("u")
        pde = PolynomialPDE.from_equations({"u": -0.2 * u + 0.1 * u * u})
        plan = QHAMPlan(pde, 2)
        disc = Discretization(pde, Grid(("x",), (2,), (1.0,)))
        concrete = gate_bindings(disc, [0.2, 0.1])
        specs = {
            key: (port.encoding.alpha, port.encoding.signal_qubits) for key, port in concrete.ports
        }
        abstract = QHAMBindings.declare(plan, 1, specs, initial_norm=concrete.initial_norm)
        model = qham_input_model(plan, abstract, eta=-0.4)
        solution = model.solve(partial(taylor_qode, degree=1), 0.01)
        names = {r.name for r in unresolved(solution.operation.program())}
        self.assertEqual(names, {"Qham_" + p.name for p in plan.pde.ports} | {"Qham_initial"})
        mapping = {
            dict(abstract.ports)[k].encoding.operation.module.name: v.encoding.operation
            for k, v in concrete.ports
        }
        mapping[abstract.initial.operation.module.name] = concrete.initial.operation
        closed = bind(solution.operation.program(), mapping)
        self.assertFalse(unresolved(closed))
        self.assertEqual(loads(dumps(closed)), closed)
        opaque = open_qham_input(
            plan,
            abstract,
            generator_alpha=model.generator.alpha,
            generator_signal=model.generator.signal_qubits,
        )
        holes = {
            r.name
            for r in unresolved(
                opaque.solve(partial(taylor_qode, degree=1), 0.01).operation.program()
            )
        }
        self.assertEqual(holes, {"QhamGenerator", "Qham_initial"})

    def test_dissipative_adaptation_is_explicit(self):
        u = Field("u")
        pde = PolynomialPDE.from_equations({"u": u + 0.1 * u * u})
        plan = QHAMPlan(pde, 1)
        disc = Discretization(pde, Grid(("x",), (2,), (1.0,)))
        model = qham_input_model(plan, gate_bindings(disc, [0.2, 0]))
        shifted = model.dissipative_shift()
        self.assertEqual(shifted.growth_shift, model.generator.alpha)
        self.assertAlmostEqual(shifted.generator.alpha, 2 * model.generator.alpha)
        with self.assertRaises(ValidationError):
            model.dissipative_shift(0.0)

    def test_multidimensional_coupled_identity(self):
        u, v = Field("u"), Field("v")
        pde = PolynomialPDE.from_equations(
            {"u": 0.1 * u.d("x", 2) - u * v.d("y", 2), "v": 0.2 * v.d("y", 2) + v * u.d("x")},
            axes=("x", "y"),
        )
        plan = QHAMPlan(pde, 1)
        disc = Discretization(pde, Grid(("x", "y"), (3, 2), (1.0, 1.0)))
        rng = random.Random(3)
        values = [[rng.random() * 0.1 for _ in range(disc.dimension)] for _ in range(2)]
        lifted = disc.lift(plan, values)
        matrix = disc.matrix(plan, -0.6)
        actual = [sum(a * b for a, b in zip(row, lifted, strict=True)) for row in matrix]
        expected = disc.chain_rule(plan, values, -0.6)
        self.assertLess(max(abs(a - b) for a, b in zip(actual, expected, strict=True)), 1e-10)

    def test_m1_reduces_to_previous_special_case(self):
        from oracq.applications.legacy import qham_lift_m1
        from oracq.applications.qham import structured_fd_bindings

        u = Field("u")
        pde = PolynomialPDE.from_equations({"u": -0.2 * u + 0.1 * u * u})
        disc = Discretization(pde, Grid(("x",), (2,), (1.0,)))
        bindings = structured_fd_bindings(disc, [0.2, 0])
        ports = dict(bindings.ports)
        nonlinear = next(v.encoding for k, v in ports.items() if k.startswith("B_"))
        old = qham_lift_m1(ports["L"].encoding, nonlinear, h=-0.4)
        new = qham_input_model(QHAMPlan(pde, 1), bindings, eta=-0.4).generator
        for column in range(16):
            a = simulate(old.operation.program(), initial={"target": column}).amplitudes
            b = simulate(new.operation.program(), initial={"target": column}).amplitudes
            for row in range(16):
                self.assertAlmostEqual(
                    a.get((row, 0), 0) * old.alpha, b.get((row, 0), 0) * new.alpha, places=10
                )

    def test_zero_initial_forcing_and_reused_work(self):
        from dataclasses import replace

        from oracq.algorithms.input_model.oracles import gate_state_prep
        from oracq.algorithms.input_model.qham import lifted_initial
        from oracq.applications.qham import structured_fd_bindings

        u = Field("u")
        pde = PolynomialPDE.from_equations({"u": u * u + Known("f")})
        plan = QHAMPlan(pde, 2)
        disc = Discretization(pde, Grid(("x",), (2,), (1.0,)), {"f": [0.1, 0.1]})
        bindings = structured_fd_bindings(disc, [0, 0])
        state, log_norm = lifted_initial(plan, bindings)
        value = simulate(state.operation.program()).amplitudes
        self.assertEqual(value, {(plan.raw_dimension(2) - 1, 0): 1 + 0j})
        self.assertEqual(log_norm, 0)
        bindings = replace(
            bindings,
            initial=gate_state_prep([0.2, 0.1], work_width=3),
            initial_norm=math.sqrt(0.05),
        )
        state, _ = lifted_initial(plan, bindings)
        self.assertEqual(state.work_width, 6)  # 3 个共享工作位 + 5 个分支所需的 3 位标签。
        self.assertTrue(
            all(work == 0 for _, work in simulate(state.operation.program()).amplitudes)
        )

    def test_lazy_stencil_rows_match_explicit_linearization(self):
        u = Field("u")
        pde = PolynomialPDE.from_equations({"u": 0.1 * u.d("x", 2) - u * u.d("x") + Known("f")})
        plan = QHAMPlan(pde, 2)
        disc = Discretization(pde, Grid(("x",), (4,), (1.0,)), {"f": [0.1, 0.2, 0, -0.1]})
        matrix = disc.matrix(plan, -0.4)
        for row in range(len(matrix)):
            entries = dict(disc.qcl_row(plan, -0.4, row))
            for column, value in enumerate(matrix[row]):
                self.assertAlmostEqual(entries.get(column, 0), value, places=11)

    def test_qram_coefficient_encoding_and_bindings(self):
        u = Field("u")
        pde = PolynomialPDE.from_equations({"u": u.d("x", 2) + 0.1 * Known("f")})
        grid = Grid(("x",), (4,), (1.0,), boundary="periodic")
        disc = Discretization(pde, grid, {"f": [0.5, -0.25, 0.0, 1.0]})
        monomial = next(
            term.monomial
            for port in disc.pde.ports
            for term in port.terms
            if term.monomial.known
        )
        encoding = qram_coefficient_encoding(disc, monomial, angle_width=8)
        self.assertAlmostEqual(encoding.alpha, 0.1)
        (db,) = [r.name for r in unresolved(encoding.operation.program())]
        words = qram_coefficient_memory(disc, monomial, angle_width=8)
        step = 2 * math.pi / (1 << 8)
        for row, value in enumerate([0.05, -0.025, 0.0, 0.1]):
            decoded = encoding.alpha * math.cos(words[row] * step / 2)
            self.assertAlmostEqual(decoded, value, delta=0.1 * step / 2)
        closed = bind(
            encoding.operation.program(),
            {db: gate_database(2, 8, [words.get(address, 0) for address in range(4)]).operation},
        )
        self.assertFalse(unresolved(closed))
        bindings = structured_fd_bindings(
            disc, [0.1, 0.2, 0.15, 0.05], coefficient_encoder=qram_coefficient_encoding
        )
        plan = QHAMPlan(pde, 2)
        self.assertIs(bindings.validate(plan), bindings)
        complex_disc = Discretization(pde, grid, {"f": [1j, 0, 0, 0]})
        with self.assertRaises(ValidationError):
            qram_coefficient_encoding(complex_disc, monomial)
