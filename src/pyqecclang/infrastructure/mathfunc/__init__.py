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
    """把普通 Python 纯函数一步编译为可逆量子模块。

    先由前端解释受限 AST 生成 MIR，再按定点格式与数学核配置降低为
    模块化 RIR。

    Args:
        function: 普通 Python 函数对象，或含若干 def 与 math/cmath 导入的源码字符串。
        fmt: ``FixedFormat`` 定点格式；省略时使用 ``FixedFormat(12, 6)``。
        inputs: 参数名到 real/complex/bool 或 ``Index(width)`` 的映射；省略时按注解推导。
        constants: 参数名到有限数值的映射，作为生成期常量绑定。
        helpers: 追加到入口命名空间的辅助函数。
        output_names: 入口各返回值的输出寄存器名；省略时单返回值为 out，多返回值依次为 out_0、out_1 等。
        config: ``MathConfig`` 数学核阶数与近似区间配置；省略时使用默认配置。
        max_unroll: 静态 range 循环允许展开的最大迭代数。
        entry: 源码字符串形式下的入口函数名；省略时使用最后一个 def。

    Returns:
        CompiledFunction: 编译结果，携带可逆模块、原 MIR 与输入/输出布局。

    Raises:
        FunctionCompileError: 源码有语法错误或包含前端不支持的语句。
        ValidationError: 输入类型、定点格式或降低阶段校验失败。
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


def lower_math_ir(program, *, fmt=None, config=None, output_names=None):
    """把已构造的 MIR 程序按给定生成配置降低为可逆 RIR 模块。

    Args:
        program: 待降低的 ``MathProgram``；构造降低器时先做完整校验。
        fmt: ``FixedFormat`` 定点格式；省略时使用 ``FixedFormat(12, 6)``。
        config: ``MathConfig`` 数学核阶数与近似区间配置；省略时使用默认配置。
        output_names: 入口各返回值的输出寄存器名；省略时按返回值个数自动命名。

    Returns:
        CompiledFunction: 入口函数的编译结果。

    Raises:
        ValidationError: MIR 校验失败，或降低阶段布局与定点格式校验失败。
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
