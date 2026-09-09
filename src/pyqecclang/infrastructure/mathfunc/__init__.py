"纯 Python 数学函数编译入口。"

from pyqecclang.algorithms.arithmetic import FixedFormat
from pyqecclang.infrastructure.mathfunc.frontend import Frontend, FunctionCompileError
from pyqecclang.infrastructure.mathfunc.graph import Index, MathProgram
from pyqecclang.infrastructure.mathfunc.lowering import CompiledFunction, Lowerer
from pyqecclang.infrastructure.mathfunc.numeric import MathConfig


def compile_function(
    function,
    *,
    fmt=None,
    inputs=None,
    constants=None,
    helpers=None,
    output_names=None,
    config=None,
    max_unroll=128,
    entry=None,
):
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


def lower_math_ir(program, *, fmt=None, config=None, output_names=None):
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
