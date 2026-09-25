"Compilation entry point for pure Python math functions."

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence

from oracq.algorithms.common.arithmetic import FixedFormat
from oracq.infrastructure.mathfunc.frontend import Frontend, FunctionCompileError, Source
from oracq.infrastructure.mathfunc.graph import Index, MathProgram
from oracq.infrastructure.mathfunc.lowering import CompiledFunction, Lowerer
from oracq.infrastructure.mathfunc.numeric import MathConfig


def compile_function(
    function: str | Callable[..., object] | Source,
    *,
    fmt: FixedFormat | None = None,
    inputs: Mapping[str, Index | type | str] | None = None,
    constants: Mapping[str, bool | int | float | complex] | None = None,
    helpers: Mapping[str, Callable[..., object]] | None = None,
    output_names: Sequence[str] | None = None,
    config: MathConfig | None = None,
    max_unroll: int = 128,
    entry: str | None = None,
) -> CompiledFunction:
    """Compile an ordinary pure Python function into a reversible quantum module in one step.

    The frontend first interprets the restricted AST into MIR, which is then
    lowered to modular RIR under the fixed-point format and math kernel
    configuration.

    Args:
        function: An ordinary Python function object, or a source string containing several defs and math/cmath imports; the module docstring and ``__future__`` imports are ignored.
        fmt: ``FixedFormat`` fixed-point format; ``FixedFormat(12, 6)`` is used when omitted.
        inputs: Mapping from parameter names to real/complex/bool or ``Index(width)``; inferred from annotations when omitted.
        constants: Mapping from parameter names to finite numeric values, bound as generation-time constants.
        helpers: Helper functions appended to the entry namespace.
        output_names: Output register names of the entry's return values; when omitted a single return is out and multiple returns are out_0, out_1 and so on.
        config: ``MathConfig`` math kernel degree and approximation interval configuration; the default configuration is used when omitted.
        max_unroll: Maximum number of iterations a static range loop may unroll.
        entry: Entry function name for the source-string form; the last def is used when omitted.

    Returns:
        CompiledFunction: The compilation result, carrying the reversible module, the original MIR and the input/output layouts.

    Raises:
        FunctionCompileError: The source has a syntax error or contains statements the frontend does not support.
        ValidationError: Validation of input kinds, the fixed-point format or the lowering stage failed.
    """
    fmt = fmt or FixedFormat(12, 6)
    config = config or MathConfig()
    program = Frontend(
        function,
        inputs=inputs,
        constants=constants,
        helpers=helpers,
        max_unroll=max_unroll,
        entry=entry,
    ).program()
    return lower_math_ir(program, fmt=fmt, config=config, output_names=output_names)


def lower_math_ir(
    program: MathProgram,
    *,
    fmt: FixedFormat | None = None,
    config: MathConfig | None = None,
    output_names: Sequence[str] | None = None,
) -> CompiledFunction:
    """Lower an already constructed MIR program into a reversible RIR module under the given generation configuration.

    Args:
        program: The ``MathProgram`` to lower; fully validated when the lowerer is constructed.
        fmt: ``FixedFormat`` fixed-point format; ``FixedFormat(12, 6)`` is used when omitted.
        config: ``MathConfig`` math kernel degree and approximation interval configuration; the default configuration is used when omitted.
        output_names: Output register names of the entry's return values; named automatically by the number of return values when omitted.

    Returns:
        CompiledFunction: The compilation result of the entry function.

    Raises:
        ValidationError: MIR validation failed, or lowering-stage layout and fixed-point format validation failed.
    """
    compiler = Lowerer(program, fmt or FixedFormat(12, 6), config or MathConfig(), output_names)
    return compiler.lower(program.entry)


__all__ = [
    "compile_function",
    "lower_math_ir",
    "CompiledFunction",
    "FunctionCompileError",
    "MathProgram",
    "MathConfig",
    "Index",
]
