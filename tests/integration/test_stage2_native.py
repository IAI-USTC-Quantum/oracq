"""第二阶段真实编译/执行冒烟；不验证 Roe 或微分方程的数值精度。"""

import math
import unittest

from oracq import Bits, Builder, ValidationError, export_originir, run_pysparq, simulate
from oracq.algorithms.common.arithmetic import (
    FixedFormat,
    arithmetic_native_registry,
    fixed_arithmetic,
)
from oracq.algorithms.input_model.oracles import declare
from oracq.infrastructure.backends.basis import export_toffoli_u3_cz
from oracq.infrastructure.native import DynamicCppFactory, NativeRegistry


class Stage2NativeTests(unittest.TestCase):
    def test_compiled_arithmetic_views_controls_and_adjoint(self):
        op = fixed_arithmetic("add", FixedFormat(3, 0, False))
        b = Builder("compiled_views", {"words": Bits(11), "control": Bits(1)})
        b.x(b["words"][0])
        b.x(b["words"][4])
        b.h(b["control"])
        kwargs = dict(
            a=b["words"][:3], b=b["words"][3:6], out=b["words"][6:9], status=b["words"][9:]
        )
        with b.control(b["control"], 0):
            b.call(op, **kwargs)
        with b.control(b["control"], 1):
            with b.adjoint():
                b.call(op, **kwargs)
        program = b.finish().program()
        report = {}
        actual = run_pysparq(
            program, native_registry=arithmetic_native_registry(program), report=report
        )
        self.assertEqual(actual.amplitudes, simulate(program).amplitudes)
        self.assertEqual(report["native_calls"], 2)
        self.assertLess(report["gate_events"], 10)

    def test_native_only_declaration_is_not_gate_closed(self):
        op = declare("NativeFlip", {"q": Bits(1)})
        source = """
class NativeFlipImpl : public SelfAdjointOperator {
 size_t q;
 public: NativeFlipImpl(size_t q_) : q(q_) {}
 void operator()(std::vector<System>& state) const override {
   for(auto &s:state) s.registers.at(q).value ^= 1ULL;
 }
};
"""
        registry = NativeRegistry().register(
            op, DynamicCppFactory("NativeFlipImpl", source, ["q"], cache_dir="out/native-cache")
        )
        self.assertEqual(
            run_pysparq(op.program(), native_registry=registry).amplitudes, {(1,): 1 + 0j}
        )
        with self.assertRaises(ValidationError):
            export_originir(op.program())

    def test_native_module_private_workspace(self):
        op = fixed_arithmetic("add", FixedFormat(2, 0, False))
        b = Builder("native_private_parent", {"a": Bits(2)})
        tmp, flag = b.local("result", Bits(2)), b.local("status", Bits(2))
        rhs = b.local("rhs", Bits(2))
        b.x(rhs[0])
        b.call(op, a=b["a"], b=rhs, out=tmp, status=flag)
        with b.adjoint():
            b.call(op, a=b["a"], b=rhs, out=tmp, status=flag)
        b.x(rhs[0])
        p = b.finish().program()
        self.assertEqual(
            run_pysparq(p, native_registry=arithmetic_native_registry(p)).amplitudes, {(0,): 1 + 0j}
        )

    def test_strict_basis_real_originir_execution(self):
        from uniqc.simulator import Simulator

        leaf = Builder("strict_leaf", {"q": Bits(1)})
        leaf.ry(leaf["q"], 0.37)
        leaf.rz(leaf["q"], 0.49)
        leaf.global_phase(0.19)
        leaf = leaf.finish()
        b = Builder("strict_main", {"q": Bits(4)})
        b.h(b["q"])
        with b.control(b["q"][:3], 5):
            with b.repeat(3):
                b.call(leaf, q=b["q"][3])
            with b.adjoint():
                b.call(leaf, q=b["q"][3])
        p = b.finish().program()
        artifact = export_toffoli_u3_cz(p)
        vector = Simulator(least_qubit_remapping=False).simulate_statevector(artifact.text)
        expected = simulate(p).statevector()
        for actual, ref in zip(vector[:16], expected, strict=True):
            self.assertAlmostEqual(complex(actual), ref, places=10)
        self.assertLess(sum(abs(x) ** 2 for x in vector[16:]), 1e-20)

    def test_roe_face_compiled_arithmetic_smoke(self):
        from oracq.applications.roe import roe_face

        fmt = FixedFormat(4, 1)
        op = roe_face(fmt=fmt)
        b = Builder("roe_native_smoke", {r.name: r.type for r in op.module.registers})
        for name in ("rho_l", "rho_r", "e_l", "e_r"):
            value = fmt.encode(1 if name.startswith("rho") else 2)
            for bit in range(fmt.width):
                if (value >> bit) & 1:
                    b.x(b[name][bit])
        b.call(op, **{r.name: b[r.name] for r in op.module.registers})
        p = b.finish().program()
        report = {}
        state = run_pysparq(p, native_registry=arithmetic_native_registry(p), report=report)
        self.assertEqual(len(state.amplitudes), 1)
        self.assertTrue(all(math.isfinite(abs(v)) for v in state.amplitudes.values()))
        self.assertGreater(report["native_calls"], 100)
