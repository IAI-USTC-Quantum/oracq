"pyqecclang：Python 生成器与模块化寄存器级 IR。"

from pyqecclang.algorithms.arithmetic import (
    FixedFormat,
    arithmetic_native_registry,
    fixed_arithmetic,
)
from pyqecclang.algorithms.contracts import (
    ContractError,
    ContractIssue,
    ContractReport,
    InputRequirement,
    OracleCapabilities,
    OracleSpec,
    ProtocolContract,
    describe_oracle,
    requires,
)
from pyqecclang.algorithms.ode import QODEProblem, QODEProtocol
from pyqecclang.algorithms.operators import (
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
from pyqecclang.algorithms.oracles import declare
from pyqecclang.algorithms.qlss import (
    BlockSystem,
    LinearSystem,
    QLSSProtocol,
    SolveResult,
    SparseSystem,
    SpectralPromise,
)
from pyqecclang.infrastructure.backends import OriginIRArtifact, export_originir, run_originir
from pyqecclang.infrastructure.backends.basis import export_toffoli_u3_cz
from pyqecclang.infrastructure.backends.pysparq import run_pysparq
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.execution import RegisterState, simulate
from pyqecclang.infrastructure.ir import (
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
from pyqecclang.infrastructure.linking import (
    Binding,
    OracleRequirement,
    bind,
    capabilities,
    unresolved,
)
from pyqecclang.infrastructure.mathfunc import (
    CompiledFunction,
    FunctionCompileError,
    Index,
    MathConfig,
    MathProgram,
    compile_function,
    lower_math_ir,
)
from pyqecclang.infrastructure.native import DynamicCppFactory, NativeRegistry
from pyqecclang.infrastructure.serialization import dumps, loads
from pyqecclang.infrastructure.validation import validate

from ._compat import install as _install_compatibility

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

__all__ += [
    "ContractError",
    "ContractIssue",
    "ContractReport",
    "InputRequirement",
    "OracleCapabilities",
    "OracleSpec",
    "ProtocolContract",
    "describe_oracle",
    "QODEProblem",
    "QODEProtocol",
]

__all__ += ["requires"]


_install_compatibility()
del _install_compatibility
