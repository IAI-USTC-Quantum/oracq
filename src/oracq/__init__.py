"oracq：Python 生成器与模块化寄存器级 IR。"

from oracq.algorithms.common.arithmetic import (
    FixedFormat,
    arithmetic_native_registry,
    fixed_arithmetic,
)
from oracq.algorithms.input_model.contracts import (
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
from oracq.algorithms.input_model.operators import (
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
from oracq.algorithms.input_model.oracles import declare
from oracq.algorithms.input_model.qdata import QMatrix, QVector
from oracq.algorithms.qlss.qlss import (
    BlockSystem,
    LinearSystem,
    QLSSProtocol,
    QLSSSolver,
    SolveResult,
    SparseSystem,
    SpectralPromise,
)
from oracq.algorithms.qml.recommendation import (
    KPRecommendationConfig,
    RecommendationResult,
    kp_recommendation,
    sigma_from_phase,
)
from oracq.algorithms.qode.ode import QODEProblem, QODEProtocol, QODESolver
from oracq.infrastructure.backends import (
    OriginIRArtifact,
    StrictArtifact,
    export_originir,
    export_strict,
    quantikz,
    run_originir,
)
from oracq.infrastructure.backends.basis import export_toffoli_u3_cz
from oracq.infrastructure.backends.pysparq import run_pysparq, run_pysparq_rir
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.estimate import (
    OracleCall,
    ResourceEstimate,
    RotationCounts,
    estimate_resources,
)
from oracq.infrastructure.execution import RegisterState, simulate
from oracq.infrastructure.ir import (
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
from oracq.infrastructure.linking import (
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
from oracq.infrastructure.mathfunc import (
    CompiledFunction,
    FunctionCompileError,
    Index,
    MathConfig,
    MathProgram,
    compile_function,
    lower_math_ir,
)
from oracq.infrastructure.native import DynamicCppFactory, NativeRegistry
from oracq.infrastructure.qmem import QMem, QPtr
from oracq.infrastructure.qram_schema import dump_qram_yaml, load_qram_yaml
from oracq.infrastructure.serialization import dumps, loads
from oracq.infrastructure.validation import validate

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
    "dump_qram_yaml",
    "dumps",
    "estimate_resources",
    "export_originir",
    "export_strict",
    "fuse",
    "identity",
    "linear_combination",
    "load_qram_yaml",
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
