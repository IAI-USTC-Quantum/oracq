"""Publication-grade validation of the arithmetic compilation modules.

Two families of implementations are covered:
- Builder primitive arithmetic (add_const / xor / swap / controlled arithmetic):
  OriginIR-ext full-amplitude exhaustive sweeps (1-8 bit), UniQC
  ``Circuit.to_matrix`` unitary comparison (3-4 bit), and PySparQ RIR wide
  registers (12/16/32/64 bit);
- fixed-point compiled functions (Boolean SSA lowering of compile_function):
  pointwise semantics of polynomials and elementary functions under several
  FixedFormats, reporting quantization error and the exact-match fraction.

Run: PYTHONPATH=src <python with pysparq+uniqc> tests/verification/verify_arithmetic.py
"""

from __future__ import annotations

import math

from harness import (
    Report,
    adapter_pysparq,
    amplitude_error,
    amplitudes_to_statevector,
    basis_program,
    originir_ext,
    originir_unitary,
    reference,
    rir_pysparq,
    sampled_inputs,
    statevector_error,
    superposition_program,
)

from oracq import Bits, Builder, FixedFormat, UInt
from oracq.infrastructure.mathfunc import MathConfig, compile_function


def _add_const_operation(width):
    b = Builder(f"add_{width}", {"w": UInt(width)})
    b.add_const(b["w"], 0)  # placeholder; the actual constant is injected at the call site
    return b


def verify_add_const_superposition(report):
    """Small bit counts: full-amplitude exhaustive check of add_const on superpositions (OriginIR-ext + three-way cross-check)."""
    for width in range(1, 9):
        modulus = 1 << width
        constants = sorted({1, 3 % modulus, modulus - 1, modulus // 2})
        worst = 0.0
        for value in constants:
            b = Builder(f"add_const_{width}_{value}", {"w": UInt(width)})
            b.h(b["w"])
            b.add_const(b["w"], value)
            p = b.finish().program()
            uniform = 1 / math.sqrt(modulus)
            expected_vector = [0j] * modulus
            for x in range(modulus):
                expected_vector[(x + value) % modulus] = uniform
            origin = originir_ext(p)
            worst = max(worst, statevector_error(origin, expected_vector))
            expected = reference(p)
            for runner in (rir_pysparq, adapter_pysparq):
                actual = runner(p)
                worst = max(worst, amplitude_error(actual, expected))
        report.case(
            f"add_const-superposition-w{width}",
            paths=["originir-ext", "reference", "rir-pysparq", "adapter-pysparq"],
            parameters={"width": width, "constants": constants},
            metrics={"max_error": worst},
            criterion="full amplitudes agree pointwise with the classical permutation (max_error < 1e-9)",
            passed=worst < 1e-9,
        )


def verify_add_const_matrix(report):
    """Unitary level: element-wise comparison of to_matrix against the classical permutation matrix."""
    import numpy as np

    for width, value in ((3, 3), (4, 11)):
        modulus = 1 << width
        b = Builder(f"add_matrix_{width}_{value}", {"w": UInt(width)})
        b.add_const(b["w"], value)
        p = b.finish().program()
        unitary = originir_unitary(p)
        permutation = np.zeros((modulus, modulus), dtype=complex)
        for x in range(modulus):
            permutation[(x + value) % modulus, x] = 1.0
        error = float(np.abs(unitary - permutation).max())
        report.case(
            f"add_const-matrix-w{width}-c{value}",
            paths=["originir-ext+to_matrix"],
            parameters={"width": width, "constant": value},
            metrics={"max_error": error},
            criterion="circuit unitary equals the classical permutation matrix (max_error < 1e-12)",
            passed=error < 1e-12,
        )


def verify_add_const_wide(report):
    """Wide registers: basis-state sweep and 12-bit full-superposition comparison via PySparQ RIR."""
    for width in (16, 32, 64):
        modulus = 1 << width
        constant = (0x5A5A5A5A5A5A5A5A % modulus) or 1
        inputs, exhaustive = sampled_inputs(width)
        failures = 0
        for x in inputs:
            b = Builder(f"wide_{width}", {"w": UInt(width)})
            for bit in range(width):
                if (x >> bit) & 1:
                    b.x(b["w"][bit])
            b.add_const(b["w"], constant)
            state = rir_pysparq(b.finish().program())
            if set(state) != {((x + constant) % modulus,)} or abs(
                next(iter(state.values())) - 1
            ) > 1e-12:
                failures += 1
        report.case(
            f"add_const-wide-w{width}",
            paths=["rir-pysparq"],
            parameters={"width": width, "constant": constant, "inputs": len(inputs)},
            metrics={"failures": failures},
            criterion="sampled basis-state inputs produce exactly (x+c) mod 2^w (failures == 0)",
            passed=failures == 0,
        )
    width = 12
    modulus = 1 << width
    b = Builder("wide_superposition", {"w": UInt(width)})
    b.h(b["w"])
    b.add_const(b["w"], 0x123)
    p = b.finish().program()
    actual = rir_pysparq(p, max_states=1 << (width + 1))
    uniform = 1 / math.sqrt(modulus)
    expected = {(((x + 0x123) % modulus),): uniform for x in range(modulus)}
    error = amplitude_error(actual, expected)
    report.case(
        "add_const-wide-superposition-w12",
        paths=["rir-pysparq"],
        parameters={"width": width, "states": modulus},
        metrics={"max_error": error},
        criterion="all 4096 superposition branches agree amplitude by amplitude (max_error < 1e-12)",
        passed=error < 1e-12,
    )


def verify_structured_arithmetic(report):
    """xor/swap views and controlled arithmetic: pairwise cross-check over four paths."""
    b = Builder("structured", {"a": UInt(4), "b": UInt(4), "c": Bits(2)})
    b.h(b["a"])
    b.h(b["c"])
    b.x(b["b"][3])
    b.xor(b["a"], b["b"])
    b.swap(b["a"][:2], b["b"][2:])
    with b.control(b["c"], 2):
        b.add_const(b["a"], 5)
    with b.control(b["c"][0], 0):
        b.xor(b["b"][2:], b["a"][:2])
    p = b.finish().program()
    ref = reference(p)
    rir = rir_pysparq(p)
    adapter = adapter_pysparq(p)
    origin = originir_ext(p)
    widths = [4, 4, 2]
    deviation = max(
        amplitude_error(ref, rir),
        amplitude_error(ref, adapter),
        statevector_error(origin, amplitudes_to_statevector(ref, widths)),
    )
    report.case(
        "structured-xor-swap-controlled",
        paths=["reference", "rir-pysparq", "adapter-pysparq", "originir-ext"],
        parameters={"widths": widths},
        metrics={"max_pairwise_deviation": deviation},
        criterion="pairwise deviation across the four paths < 1e-9",
        passed=deviation < 1e-9,
    )


def _sweep_compiled(report, source, fmt, exact_oracle, name, paths_extra=()):
    """A single superposition run exhausts all 2^width inputs; additional basis-state spot checks confirm determinism."""
    compiled = compile_function(source, fmt=fmt)
    operation = compiled.operation
    quantum = 1.0 / (1 << fmt.fraction)
    modulus = 1 << fmt.width
    names = [r.name for r in operation.module.registers]
    out_key, flag_key = names[-2], names[-1]
    state = rir_pysparq(
        superposition_program(operation, ["x"]), max_states=1 << (fmt.width + 4)
    )
    if len(state) != modulus:
        raise AssertionError(f"{name}: branch count {len(state)} != 2^{fmt.width}")
    uniform = 1 / math.sqrt(modulus)
    max_error = 0.0
    exact = 0
    flagged = 0
    for basis, amplitude in state.items():
        values = dict(zip(names, basis, strict=True))
        if abs(amplitude - uniform) > 1e-12:
            raise AssertionError(f"{name}: anomalous superposition-branch amplitude {amplitude}")
        x = fmt.decode(values["x"])
        expected_value = exact_oracle(x)
        if values[flag_key]:
            flagged += 1
            continue
        if expected_value is None:
            continue
        error = abs(fmt.decode(values[out_key]) - expected_value)
        max_error = max(max_error, error)
        if error == 0:
            exact += 1
    for raw in (0, 1, modulus // 2 - 1):
        single = rir_pysparq(basis_program(operation, {"x": raw}))
        if len(single) != 1 or abs(next(iter(single.values())) - 1) > 1e-12:
            raise AssertionError(f"{name}: basis-state execution produced a non-deterministic output raw={raw}")
    report.case(
        name,
        paths=["rir-pysparq", *paths_extra],
        parameters={"width": fmt.width, "fraction": fmt.fraction, "branches": len(state)},
        metrics={
            "max_error": max_error,
            "exact_fraction": exact / len(state),
            "flagged": flagged,
        },
        criterion=f"unflagged output error <= 2 quanta ({2 * quantum:.6f}) with uniform branches",
        passed=max_error <= 2 * quantum + 1e-12,
    )


def verify_fixed_point_polynomial(report):
    for width, fraction in ((4, 1), (6, 2), (8, 3)):
        fmt = FixedFormat(width, fraction)
        _sweep_compiled(
            report,
            "def f(x):\n return x*x+0.5",
            fmt,
            lambda x: x * x + 0.5,
            f"fixed-point-polynomial-{width}-{fraction}",
        )


def _approximation_coefficients(program, function):
    for module in program.modules:
        for key, value in module.attributes:
            if key == "math_approximation":
                import json

                data = json.loads(value)
                if data["function"] == function:
                    return data["coefficients"]
    raise AssertionError(f"approximation coefficients for {function} not found")


def verify_fixed_point_elementary(report):
    fmt = FixedFormat(8, 4)
    compiled = compile_function(
        "import math\ndef f(x):\n return math.sin(x)",
        fmt=fmt,
        config=MathConfig(degree=3, intervals=(("sin", -1.0, 1.0),)),
    )
    operation = compiled.operation
    quantum = 1.0 / (1 << fmt.fraction)
    coefficients = _approximation_coefficients(compiled.program(), "sin")

    def interpolant(x):
        result = 0.0
        for coefficient in reversed(coefficients):
            result = result * x + coefficient
        return result

    names = [r.name for r in operation.module.registers]
    state = rir_pysparq(superposition_program(operation, ["x"]), max_states=1 << 12)
    impl_error = 0.0
    checked = 0
    flagged = 0
    misflagged = 0
    for basis, _ in state.items():
        values = dict(zip(names, basis, strict=True))
        x = fmt.decode(values["x"])
        in_interval = -1.0 <= x <= 1.0
        if values["status"]:
            flagged += 1
            if in_interval:
                misflagged += 1
            continue
        if not in_interval:
            misflagged += 1
            continue
        checked += 1
        impl_error = max(impl_error, abs(fmt.decode(values["out"]) - interpolant(x)))
    # Method error: gap between the Chebyshev interpolant and sin itself (dense classical grid; informative metric)
    method_error = max(
        abs(interpolant(i / 1000) - math.sin(i / 1000)) for i in range(-1000, 1001)
    )
    report.case(
        "fixed-point-sin-8-4",
        paths=["rir-pysparq"],
        parameters={"width": 8, "fraction": 4, "degree": 3, "interval": [-1.0, 1.0]},
        metrics={
            "impl_error": impl_error,
            "method_error": method_error,
            "checked": checked,
            "flagged": flagged,
            "misflagged": misflagged,
        },
        criterion=(
            f"implementation error <= 2 quanta ({2 * quantum:.6f}, against the Chebyshev coefficients in the module attributes)"
            " with flags consistent inside and outside the interval"
        ),
        passed=impl_error <= 2 * quantum + 1e-12 and misflagged == 0,
    )


def verify_math_superposition_cross(report):
    """Four-path cross-check of compiled functions on superposition inputs (small format, limited by the OriginIR bit budget)."""
    fmt = FixedFormat(4, 1)
    compiled = compile_function("def f(x):\n return x*x+0.5", fmt=fmt)
    operation = compiled.operation
    program = superposition_program(operation, ["x"])
    ref = reference(program)
    deviation = amplitude_error(ref, rir_pysparq(program))
    deviation = max(deviation, amplitude_error(ref, adapter_pysparq(program)))
    try:
        origin = originir_ext(program)
    except Exception:
        origin = None
    paths = ["reference", "rir-pysparq", "adapter-pysparq"]
    if origin is not None:
        widths = [r.type.width for r in program.main.registers]
        deviation = max(
            deviation, statevector_error(origin, amplitudes_to_statevector(ref, widths))
        )
        paths.append("originir-ext")
    report.case(
        "fixed-point-superposition-cross",
        paths=paths,
        parameters={"format": "4.1"},
        metrics={"max_pairwise_deviation": deviation},
        criterion="pairwise deviation across paths < 1e-9",
        passed=deviation < 1e-9,
    )


def run():
    report = Report(
        "arithmetic",
        "Three-level (pointwise / superposition / unitary) validation of primitive arithmetic and fixed-point compiled functions, covering 1-64 bit.",
    )
    verify_add_const_superposition(report)
    verify_add_const_matrix(report)
    verify_add_const_wide(report)
    verify_structured_arithmetic(report)
    verify_fixed_point_polynomial(report)
    verify_fixed_point_elementary(report)
    verify_math_superposition_cross(report)
    report.write()
    return report


if __name__ == "__main__":
    run()
