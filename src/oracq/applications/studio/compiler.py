"""Compile bounded pure-function examples and compare checkpoint trade-offs.

This is the native compiler pipeline, not a fitted resource model. Logical
space is read from the actual strict backend layout, including its scratch.
"""
from __future__ import annotations

import ast
import hashlib
import time
from dataclasses import asdict, dataclass
from importlib.metadata import version

from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.infrastructure.backends.strict import export_strict
from oracq.infrastructure.estimate import estimate_resources
from oracq.infrastructure.mathfunc import MathConfig, compile_function
from oracq.infrastructure.serialization import dumps

DEFAULT_SOURCE = """def f(x):
    a = x + 1
    b = a * a
    c = b + 1
    d = c * c
    return d if x < 0 else -c
"""


@dataclass(frozen=True)
class CompilationRequest:
    """One numeric configuration shared by both resource comparison variants."""

    source: str = DEFAULT_SOURCE
    width: int = 8
    fraction: int = 3
    degree: int = 4
    checkpoint_size: int = 2

    def validate(self) -> None:
        """Keep GUI workloads bounded; the ordinary compiler CLI remains general."""
        if not isinstance(self.source, str) or not 1 <= len(self.source.encode()) <= 6000:
            raise ValueError("source must contain 1..6000 UTF-8 bytes")
        limits = {"width": (2, 16), "degree": (1, 8), "checkpoint_size": (0, 64)}
        for name, (lo, hi) in limits.items():
            value = getattr(self, name)
            if type(value) is not int or not lo <= value <= hi:
                raise ValueError(f"{name} must be an integer in [{lo}, {hi}]")
        if type(self.fraction) is not int or not 0 <= self.fraction < self.width-1:
            raise ValueError("fraction must be an integer in [0,width-1)")
        try:
            tree = ast.parse(self.source)
        except SyntaxError as exc:
            raise ValueError(f"source syntax error: {exc.msg}") from exc
        nodes = list(ast.walk(tree))
        if len(nodes) > 128:
            raise ValueError("studio supports at most 128 AST nodes; use compile-function for larger sources")
        if any(isinstance(node, (ast.For, ast.While, ast.AsyncFor, ast.ListComp,
                                 ast.DictComp, ast.SetComp, ast.GeneratorExp)) for node in nodes):
            raise ValueError("studio examples cannot contain loops; use compile-function for static loops")


def compile_comparison(request: CompilationRequest) -> dict[str, object]:
    """Compare retained intermediates with selected checkpointing at fixed precision.

    Returns resource summaries and the selected program's standalone JSON RIR,
    MIR and strict modular netlist. No source function or hardware job is run.
    """
    request.validate()
    started = time.perf_counter()
    summaries: dict[str, dict[str, object]] = {}
    artifacts: dict[str, str] = {}
    for label, checkpoint in (("baseline", 0), ("selected", request.checkpoint_size)):
        tick = time.perf_counter()
        compiled = compile_function(
            request.source, fmt=FixedFormat(request.width, request.fraction),
            config=MathConfig(degree=request.degree, checkpoint_size=checkpoint), max_unroll=16,
        )
        program = compiled.program()
        ledger = estimate_resources(program)
        strict = export_strict(program)
        wires = [wire for register in strict.registers.values() for wire in register]
        wires.extend(strict.workspace_qubits)
        logical_qubits = max(wires, default=-1)+1
        summaries[label] = {
            "checkpoint_size": checkpoint, "complete": ledger.complete,
            "logical_qubits": logical_qubits, "rir_qubits": ledger.qubits,
            "backend_scratch_qubits": logical_qubits-ledger.qubits_lower_bound,
            "toffoli": ledger.toffoli, "clifford": ledger.clifford,
            "t_exact": ledger.t_exact, "rotations_pending_synthesis": ledger.rotations.total,
            "qram_queries": ledger.qram_total, "compile_and_export_s": time.perf_counter()-tick,
        }
        if label == "selected":
            artifacts = {"program.rir.json": dumps(program, format="json"),
                         "program.strict.originir": strict.text,
                         "program.mir.json": compiled.math_ir.dumps(), "source.py": request.source}
    if sum(len(text.encode()) for text in artifacts.values()) > 8*1024*1024:
        raise ValueError("studio export exceeds 8 MiB; use compile-function for larger programs")
    return {
        "schema": "oracq.studio.v1", "version": version("oracq"),
        "configuration": asdict(request),
        "source_sha256": hashlib.sha256(request.source.encode()).hexdigest(),
        "numeric_configuration_shared": True, "variants": summaries,
        "artifacts": artifacts, "total_s": time.perf_counter()-started,
        "scope": {"level": "logical Toffoli/Clifford/exact-T/unsynthesized rotations",
                  "numerical_error_certified": False, "fault_tolerant_schedule": False,
                  "hardware_run": False},
    }
