"0.7 导入路径的集中兼容表；仓内代码只使用规范路径。"

import importlib
import sys
from types import ModuleType

ALIASES = {
    "pyqecclang.ir": "pyqecclang.infrastructure.ir",
    "pyqecclang.builder": "pyqecclang.infrastructure.builder",
    "pyqecclang.validation": "pyqecclang.infrastructure.validation",
    "pyqecclang.serialization": "pyqecclang.infrastructure.serialization",
    "pyqecclang.linking": "pyqecclang.infrastructure.linking",
    "pyqecclang.execution": "pyqecclang.infrastructure.execution",
    "pyqecclang.native": "pyqecclang.infrastructure.native",
    "pyqecclang.layout": "pyqecclang.infrastructure.layout",
    "pyqecclang.readout": "pyqecclang.infrastructure.readout",
    "pyqecclang.backends.originir": "pyqecclang.infrastructure.backends.originir",
    "pyqecclang.backends": "pyqecclang.infrastructure.backends",
    "pyqecclang.backends.pysparq": "pyqecclang.infrastructure.backends.pysparq",
    "pyqecclang.backends.basis": "pyqecclang.infrastructure.backends.basis",
    "pyqecclang.mathfunc.graph": "pyqecclang.infrastructure.mathfunc.graph",
    "pyqecclang.mathfunc.frontend": "pyqecclang.infrastructure.mathfunc.frontend",
    "pyqecclang.mathfunc.numeric": "pyqecclang.infrastructure.mathfunc.numeric",
    "pyqecclang.mathfunc.lowering": "pyqecclang.infrastructure.mathfunc.lowering",
    "pyqecclang.mathfunc": "pyqecclang.infrastructure.mathfunc",
    "pyqecclang.mathfunc.roe_formulas": "pyqecclang.applications.roe_formulas",
    "pyqecclang.library": "pyqecclang.algorithms.operators",
    "pyqecclang.combinators": "pyqecclang.algorithms.block_encoding",
    "pyqecclang.contracts": "pyqecclang.algorithms.contracts",
    "pyqecclang.oracles": "pyqecclang.algorithms.oracles",
    "pyqecclang.arithmetic": "pyqecclang.algorithms.arithmetic",
    "pyqecclang.sparse_models": "pyqecclang.algorithms.sparse",
    "pyqecclang.access": "pyqecclang.algorithms.sparse",
    "pyqecclang.qlss": "pyqecclang.algorithms.qlss",
    "pyqecclang.algorithms.costa": "pyqecclang.algorithms.qlss",
    "pyqecclang.algorithms.cks": "pyqecclang.algorithms.qlss",
    "pyqecclang.qode": "pyqecclang.algorithms.ode",
    "pyqecclang.qfvm": "pyqecclang.applications.qfvm",
    "pyqecclang.qfvm_sparse": "pyqecclang.applications.qfvm",
    "pyqecclang.flow_data": "pyqecclang.applications.flow_data",
    "pyqecclang.roe": "pyqecclang.applications.roe",
    "pyqecclang.workloads": "pyqecclang.applications.catalog",
    "pyqecclang.qham.pde": "pyqecclang.applications.qham.pde",
    "pyqecclang.qham.linearization": "pyqecclang.applications.qham.linearization",
    "pyqecclang.qham.reference": "pyqecclang.applications.qham.reference",
    "pyqecclang.qham.quantum": "pyqecclang.algorithms.qham",
    "pyqecclang.qham.examples": "pyqecclang.applications.qham.examples",
    "pyqecclang.qham.report": "pyqecclang.applications.qham.report",
    "pyqecclang.qham.stencils": "pyqecclang.applications.qham.stencils",
}
SPLIT_EXPORTS = {
    "pyqecclang.algorithms.elementary": {
        "deutsch_jozsa": ("pyqecclang.algorithms.oracle_algorithms", "deutsch_jozsa"),
        "grover": ("pyqecclang.algorithms.search", "grover"),
        "phase_from_database": ("pyqecclang.algorithms.search", "phase_from_database"),
        "qft": ("pyqecclang.algorithms.fourier", "qft"),
        "phase_estimation": ("pyqecclang.algorithms.estimation", "phase_estimation"),
        "qubitization_walk": ("pyqecclang.algorithms.transforms", "qubitization_walk"),
        "qsvt_sequence": ("pyqecclang.algorithms.transforms", "qsvt_sequence"),
        "oblivious_amplification": ("pyqecclang.algorithms.transforms", "oblivious_amplification"),
    },
    "pyqecclang.algorithms.differential": {
        "tagged": ("pyqecclang.algorithms._dynamics", "tagged"),
        "_lcu_dynamics": ("pyqecclang.algorithms._dynamics", "_lcu_dynamics"),
        "HermitianParts": ("pyqecclang.algorithms.ode_models", "HermitianParts"),
        "LinearODE": ("pyqecclang.algorithms.ode_models", "LinearODE"),
        "QuadraturePlan": ("pyqecclang.algorithms.lchs", "QuadraturePlan"),
        "lchs_qode": ("pyqecclang.algorithms.lchs", "lchs_qode"),
        "ContourPlan": ("pyqecclang.algorithms.cbmd", "ContourPlan"),
        "cbmd_qode": ("pyqecclang.algorithms.cbmd", "cbmd_qode"),
        "cbmd_function": ("pyqecclang.algorithms.cbmd", "cbmd_function"),
        "taylor_hamiltonian": ("pyqecclang.algorithms.hamiltonian", "taylor_hamiltonian"),
        "qft": ("pyqecclang.algorithms.fourier", "qft_with_work"),
        "fourier_momentum": ("pyqecclang.algorithms.schrodingerization", "fourier_momentum"),
        "SchrodingerPlan": ("pyqecclang.algorithms.schrodingerization", "SchrodingerPlan"),
        "schrodinger_qode": ("pyqecclang.algorithms.schrodingerization", "schrodinger_qode"),
        "PolynomialODE": ("pyqecclang.algorithms.carleman", "PolynomialODE"),
        "_carleman_term": ("pyqecclang.algorithms.carleman", "_carleman_term"),
        "carleman_lift": ("pyqecclang.algorithms.carleman", "carleman_lift"),
        "carleman_initial": ("pyqecclang.algorithms.carleman", "carleman_initial"),
        "carleman_qode": ("pyqecclang.algorithms.carleman", "carleman_qode"),
        "linear_qode": ("pyqecclang.algorithms.ode", "linear_qode"),
        "PDEInput": ("pyqecclang.algorithms.pde", "PDEInput"),
        "qpde_solver": ("pyqecclang.algorithms.pde", "qpde_solver"),
    },
    "pyqecclang.algorithms.solvers": {
        "extend_initial": ("pyqecclang.algorithms.state_preparation", "extend_initial"),
        "select_subspace": ("pyqecclang.algorithms.state_preparation", "select_subspace"),
        "apply_be_to_state": ("pyqecclang.algorithms.state_preparation", "apply_be_to_state"),
        "make_euler_history_qode": ("pyqecclang.algorithms.ode", "make_euler_history_qode"),
        "DiscretePDE": ("pyqecclang.algorithms.pde", "DiscretePDE"),
        "make_qpde": ("pyqecclang.algorithms.pde", "make_qpde"),
        "make_lchs_qode": ("pyqecclang.algorithms.legacy", "make_lchs_qode"),
        "make_schrodingerisation_qode": (
            "pyqecclang.algorithms.legacy",
            "make_schrodingerisation_qode",
        ),
        "trotter_hamsim": ("pyqecclang.algorithms.hamiltonian", "trotter_hamsim"),
    },
}


def install():
    for old, new in ALIASES.items():
        module = importlib.import_module(new)
        sys.modules[old] = module
    for old, exports in SPLIT_EXPORTS.items():
        module = ModuleType(old)
        module.__package__ = old.rpartition(".")[0]
        for name, (owner, symbol) in exports.items():
            setattr(module, name, getattr(importlib.import_module(owner), symbol))
        module.__all__ = [n for n in exports if not n.startswith("_")]
        sys.modules[old] = module
    for old in (*ALIASES, *SPLIT_EXPORTS):
        parent, _, leaf = old.rpartition(".")
        owner = sys.modules.get(parent) or importlib.import_module(parent)
        setattr(owner, leaf, sys.modules[old])
