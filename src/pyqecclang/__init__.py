"""pyqecclang：Python 生成器与模块化寄存器级 IR。"""

from .arithmetic import FixedFormat, arithmetic_native_registry, fixed_arithmetic
from .backends import OriginIRArtifact, export_originir, run_originir
from .backends.basis import export_toffoli_u3_cz
from .backends.pysparq import run_pysparq
from .builder import Builder, Operation
from .execution import RegisterState, simulate
from .ir import (
    QRAM,
    VERSION,
    Adjoint,
    Bits,
    Call,
    Control,
    Load,
    Module,
    Primitive,
    Program,
    Rational,
    Ref,
    Register,
    RegType,
    Repeat,
    Resource,
    SInt,
    Span,
    UInt,
    ValidationError,
    fuse,
)
from .library import (
    BlockEncoding,
    Generator,
    block_encoding,
    identity,
    linear_combination,
    pauli_x,
    product,
    scale,
    zero,
)
from .linking import Binding, OracleRequirement, bind, capabilities, unresolved
from .mathfunc import (
    CompiledFunction,
    FunctionCompileError,
    Index,
    MathConfig,
    MathProgram,
    compile_function,
    lower_math_ir,
)
from .native import DynamicCppFactory, NativeRegistry
from .oracles import declare
from .qlss import (
    BlockSystem,
    LinearSystem,
    QLSSProtocol,
    SolveResult,
    SparseSystem,
    SpectralPromise,
)
from .serialization import dumps, loads
from .validation import validate

__all__ = [
    "VERSION",
    "Adjoint",
    "Bits",
    "BlockEncoding",
    "Builder",
    "Call",
    "Control",
    "Generator",
    "Load",
    "Module",
    "Operation",
    "OriginIRArtifact",
    "Primitive",
    "Program",
    "QRAM",
    "Rational",
    "Ref",
    "Register",
    "RegisterState",
    "RegType",
    "Repeat",
    "Resource",
    "SInt",
    "Span",
    "UInt",
    "ValidationError",
    "block_encoding",
    "dumps",
    "export_originir",
    "fuse",
    "identity",
    "linear_combination",
    "loads",
    "pauli_x",
    "product",
    "run_originir",
    "run_pysparq",
    "scale",
    "simulate",
    "validate",
    "zero",
]


__all__ += ["Binding", "OracleRequirement", "bind", "capabilities", "declare", "unresolved"]


__all__ += [
    "FixedFormat",
    "arithmetic_native_registry",
    "fixed_arithmetic",
    "export_toffoli_u3_cz",
    "DynamicCppFactory",
    "NativeRegistry",
]

__all__ += [
    "CompiledFunction",
    "FunctionCompileError",
    "Index",
    "MathConfig",
    "MathProgram",
    "compile_function",
    "lower_math_ir",
]

__all__ += [
    "BlockSystem",
    "LinearSystem",
    "QLSSProtocol",
    "SolveResult",
    "SparseSystem",
    "SpectralPromise",
]
