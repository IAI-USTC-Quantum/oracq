"""一般 QHAM 的真实寄存器执行和模块化后端描述见证。"""

import math
import unittest
from functools import partial

from pyqecclang import arithmetic_native_registry, export_toffoli_u3_cz, run_pysparq
from pyqecclang.qham import (
    Block,
    Discretization,
    Field,
    Grid,
    Known,
    PolynomialPDE,
    QHAMPlan,
    qham_input_model,
    structured_fd_bindings,
    taylor_qode,
)


class QhamNativeTests(unittest.TestCase):
    def test_structured_nonlinear_ports_beyond_dense_binding_budget(self):
        u, v = Field("u"), Field("v")
        pde = PolynomialPDE.from_equations({"u": -u * u.d("x"), "v": u * v})
        disc = Discretization(pde, Grid(("x",), (4,), (1.0,)))
        ports = structured_fd_bindings(disc, [0.1] + [0.0] * 7)
        for name, port in ports.ports:
            column = (
                8
                if pde.ports[[p.name for p in pde.ports].index(name)].terms[0].output == "u"
                else 32
            )
            from pyqecclang import Builder

            b = Builder(
                "port_witness", {r.name: r.type for r in port.encoding.operation.module.registers}
            )
            for bit in range(port.encoding.width):
                if (column >> bit) & 1:
                    b.x(b["target"][bit])
            b.call(port.encoding.operation, target=b["target"], signal=b["signal"])
            result = run_pysparq(b.finish().program())
            self.assertGreater(port.encoding.width, 5)
            for row in range(disc.dimension):
                self.assertAlmostEqual(
                    result.amplitudes.get((row, 0), 0) * port.encoding.alpha,
                    disc.entry(name, row, column),
                    places=10,
                )

    def test_generated_lift_and_taylor_solution_match_finite_polynomial(self):
        u = Field("u")
        pde = PolynomialPDE.from_equations({"u": -0.2 * u + 0.1 * u * u})
        plan = QHAMPlan(pde, 2)
        disc = Discretization(pde, Grid(("x",), (2,), (1.0,)))
        bindings = structured_fd_bindings(disc, [0.2, 0])
        model = qham_input_model(plan, bindings, eta=-0.4)
        state = model.solve(partial(taylor_qode, degree=1), 0.01)
        program = state.operation.program()
        result = run_pysparq(program, native_registry=arithmetic_native_registry(program))
        initial = disc.lift(plan, [[0.2, 0], [0, 0], [0, 0]])
        matrix = disc.matrix(plan, -0.4)
        expected = [
            initial[row] + 0.01 * sum(a * b for a, b in zip(matrix[row], initial, strict=True))
            for row in range(2)
        ]
        normalization = math.exp(model.log_initial_norm) * (1 + 0.01 * model.generator.alpha)
        for row, value in enumerate(expected):
            self.assertAlmostEqual(
                result.amplitudes.get((row, 0), 0), value / normalization, places=10
            )

    def test_forcing_generator_and_actual_originir_parser(self):
        from uniqc.compile.originir.originir_base_parser import OriginIR_BaseParser

        from pyqecclang import Builder

        u = Field("u")
        pde = PolynomialPDE.from_equations({"u": -0.2 * u + 0.1 * u * u + Known("f")})
        plan = QHAMPlan(pde, 2)
        disc = Discretization(pde, Grid(("x",), (2,), (1.0,)), {"f": [0.03, 0.02]})
        model = qham_input_model(plan, structured_fd_bindings(disc, [0.2, 0]), eta=-0.4)
        matrix = disc.matrix(plan, -0.4)
        for block in (Block("tensor", (0,)), Block("tensor", (1,)), Block("tensor")):
            column = plan.offset(block, disc.dimension)
            b = Builder(
                "forcing_generator_column",
                {r.name: r.type for r in model.generator.operation.module.registers},
            )
            for bit in range(model.generator.width):
                if (column >> bit) & 1:
                    b.x(b["target"][bit])
            b.call(model.generator.operation, target=b["target"], signal=b["signal"])
            result = run_pysparq(b.finish().program())
            for row in range(len(matrix)):
                self.assertAlmostEqual(
                    result.amplitudes.get((row, 0), 0) * model.generator.alpha,
                    matrix[row][column],
                    places=10,
                )
        solution = model.solve(partial(taylor_qode, degree=1), 0.01)
        artifact = export_toffoli_u3_cz(solution.operation.program())
        parser = OriginIR_BaseParser()
        parser.parse(artifact.text)
        self.assertGreater(parser.n_qubit, model.generator.width)
