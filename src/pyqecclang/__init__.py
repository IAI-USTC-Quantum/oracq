"pyqecclang：Python 生成器与模块化寄存器级 IR。"

from pyqecclang.algorithms.common.arithmetic import (
    FixedFormat,
    arithmetic_native_registry,
    fixed_arithmetic,
)
from pyqecclang.algorithms.input_model.contracts import (
    AlgorithmContract,
    ContractError,
    ContractIssue,
    ContractReport,
    InputRequirement,
    OracleCapabilities,
    OracleSpec,
    ProtocolContract,
    ResolvedInputs,
    describe_oracle,
    requires,
)
from pyqecclang.algorithms.input_model.operators import (
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
from pyqecclang.algorithms.input_model.oracles import declare
from pyqecclang.algorithms.input_model.qdata import QMatrix, QVector
from pyqecclang.algorithms.qlss.qlss import (
    BlockSystem,
    LinearSystem,
    QLSSProtocol,
    QLSSSolver,
    SolveResult,
    SparseSystem,
    SpectralPromise,
)
from pyqecclang.algorithms.qml.recommendation import (
    KPRecommendationConfig,
    RecommendationResult,
    kp_recommendation,
    sigma_from_phase,
)
from pyqecclang.algorithms.qode.ode import QODEProblem, QODEProtocol, QODESolver
from pyqecclang.infrastructure.backends import (
    OriginIRArtifact,
    StrictArtifact,
    export_originir,
    export_strict,
    quantikz,
    run_originir,
)
from pyqecclang.infrastructure.backends.basis import export_toffoli_u3_cz
from pyqecclang.infrastructure.backends.pysparq import run_pysparq, run_pysparq_rir
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.estimate import (
    OracleCall,
    ResourceEstimate,
    RotationCounts,
    estimate_resources,
)
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
    Store,
    UInt,
    ValidationError,
    fuse,
)
from pyqecclang.infrastructure.linking import (
    Binding,
    BindingError,
    BindingReport,
    BindingResult,
    OracleRequirement,
    bind,
    bind_with_report,
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
from pyqecclang.infrastructure.qmem import QMem, QPtr
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
    "QMem",
    "QPtr",
    "ResourceEstimate",
    "StrictArtifact",
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
    "Store",
    "UInt",
    "ValidationError",
    "block_encoding",
    "dumps",
    "estimate_resources",
    "export_originir",
    "export_strict",
    "fuse",
    "identity",
    "linear_combination",
    "loads",
    "pauli_x",
    "product",
    "quantikz",
    "run_originir",
    "run_pysparq",
    "run_pysparq_rir",
    "scale",
    "simulate",
    "validate",
    "zero",
]


__all__ += ["Binding", "OracleRequirement", "bind", "capabilities", "declare", "unresolved"]
__all__ += ["BindingError", "BindingReport", "BindingResult", "bind_with_report", "OracleCall", "RotationCounts"]


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
    "AlgorithmContract",
    "ResolvedInputs",
    "QLSSSolver",
    "QODESolver",
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


__all__ += ["QMatrix", "QVector"]


__all__ += [
    "KPRecommendationConfig",
    "RecommendationResult",
    "kp_recommendation",
    "sigma_from_phase",
]

__all__ += ["requires"]


_install_compatibility()
del _install_compatibility
