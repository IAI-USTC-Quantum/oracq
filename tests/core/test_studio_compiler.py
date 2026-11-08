"""Actual compiler, backend layout and standalone artifact contracts for the UI."""
import hashlib
import tempfile
import unittest
from pathlib import Path

from oracq import estimate_resources, loads, simulate
from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.applications.studio.compiler import CompilationRequest, compile_comparison
from oracq.infrastructure.mathfunc import MathConfig, compile_function


class StudioCompilerTests(unittest.TestCase):
    def test_resource_comparison_uses_real_backend_workspace(self):
        request = CompilationRequest(width=4, fraction=1)
        report = compile_comparison(request)
        selected = report["variants"]["selected"]
        program = loads(report["artifacts"]["program.rir.json"])
        ledger = estimate_resources(program)
        netlist = report["artifacts"]["program.strict.originir"]
        qinit = next(line for line in netlist.splitlines() if line.startswith("QINIT "))
        self.assertEqual(selected["logical_qubits"], int(qinit.split()[1]))
        self.assertEqual(selected["toffoli"], ledger.toffoli)
        self.assertEqual(selected["rir_qubits"], ledger.qubits)
        self.assertGreaterEqual(selected["backend_scratch_qubits"], 0)
        self.assertEqual(report["source_sha256"], hashlib.sha256(request.source.encode()).hexdigest())
        self.assertFalse(report["scope"]["fault_tolerant_schedule"])

    def test_selected_artifact_matches_reference_on_all_small_inputs(self):
        source = "def f(x):\n a=x+1\n b=a*a\n c=b+1\n return c if x<0 else -b"
        report = compile_comparison(CompilationRequest(source, 3, 0, 3, 2))
        selected = loads(report["artifacts"]["program.rir.json"])
        reference = compile_function(source, fmt=FixedFormat(3, 0),
                                     config=MathConfig(degree=3, checkpoint_size=0)).program()
        for value in range(8):
            initial = {"x": value, "out": 2, "status": 1}
            self.assertEqual(simulate(selected, initial=initial).amplitudes,
                             simulate(reference, initial=initial).amplitudes)

    def test_invalid_bounded_requests(self):
        for kwargs in ({"width": 1}, {"width": True}, {"width": 17}, {"fraction": 7},
                       {"degree": 9}, {"checkpoint_size": -1}, {"source": ""},
                       {"source": "def f(x):\n for i in range(2):\n  x=x+1\n return x"},
                       {"source": "def f(:"}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                CompilationRequest(**kwargs).validate()

    def test_source_is_compiled_and_never_executed_as_python(self):
        with tempfile.TemporaryDirectory() as folder:
            marker = Path(folder)/"side-effect"
            source = f"def f(x):\n open({str(marker)!r}, 'w').write('bad')\n return x"
            with self.assertRaises(ValueError):
                compile_comparison(CompilationRequest(source))
            self.assertFalse(marker.exists())
