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
    "pyqecclang.library": "pyqecclang.algorithms.input_model.operators",
    "pyqecclang.combinators": "pyqecclang.algorithms.input_model.block_encoding",
    "pyqecclang.contracts": "pyqecclang.algorithms.input_model.contracts",
    "pyqecclang.oracles": "pyqecclang.algorithms.input_model.oracles",
    "pyqecclang.arithmetic": "pyqecclang.algorithms.common.arithmetic",
    "pyqecclang.sparse_models": "pyqecclang.algorithms.input_model.sparse",
    "pyqecclang.access": "pyqecclang.algorithms.input_model.sparse",
    "pyqecclang.qlss": "pyqecclang.algorithms.qlss.qlss",
    "pyqecclang.algorithms.costa": "pyqecclang.algorithms.qlss.qlss",
    "pyqecclang.algorithms.cks": "pyqecclang.algorithms.qlss.qlss",
    "pyqecclang.qode": "pyqecclang.algorithms.qode.ode",
    "pyqecclang.qfvm": "pyqecclang.applications.qfvm",
    "pyqecclang.qfvm_sparse": "pyqecclang.applications.qfvm",
    "pyqecclang.flow_data": "pyqecclang.applications.flow_data",
    "pyqecclang.roe": "pyqecclang.applications.roe",
    "pyqecclang.workloads": "pyqecclang.applications.catalog",
    "pyqecclang.qham.pde": "pyqecclang.applications.qham.pde",
    "pyqecclang.qham.linearization": "pyqecclang.applications.qham.linearization",
    "pyqecclang.qham.reference": "pyqecclang.applications.qham.reference",
    "pyqecclang.qham.quantum": "pyqecclang.algorithms.input_model.qham",
    "pyqecclang.qham.examples": "pyqecclang.applications.qham.examples",
    "pyqecclang.qham.report": "pyqecclang.applications.qham.report",
    "pyqecclang.qham.stencils": "pyqecclang.applications.qham.stencils",
}
SPLIT_EXPORTS = {
    "pyqecclang.algorithms.elementary": {
        "deutsch_jozsa": ("pyqecclang.algorithms.basics.oracle_algorithms", "deutsch_jozsa"),
        "grover": ("pyqecclang.algorithms.common.search", "grover"),
        "phase_from_database": ("pyqecclang.algorithms.common.search", "phase_from_database"),
        "qft": ("pyqecclang.algorithms.common.fourier", "qft"),
        "phase_estimation": ("pyqecclang.algorithms.common.estimation", "phase_estimation"),
        "qubitization_walk": ("pyqecclang.algorithms.common.transforms", "qubitization_walk"),
        "qsvt_sequence": ("pyqecclang.algorithms.common.transforms", "qsvt_sequence"),
        "oblivious_amplification": ("pyqecclang.algorithms.common.transforms", "oblivious_amplification"),
    },
    "pyqecclang.algorithms.differential": {
        "tagged": ("pyqecclang.algorithms.qode._dynamics", "tagged"),
        "_lcu_dynamics": ("pyqecclang.algorithms.qode._dynamics", "_lcu_dynamics"),
        "HermitianParts": ("pyqecclang.algorithms.qode.ode_models", "HermitianParts"),
        "LinearODE": ("pyqecclang.algorithms.qode.ode_models", "LinearODE"),
        "QuadraturePlan": ("pyqecclang.algorithms.qode.lchs", "QuadraturePlan"),
        "lchs_qode": ("pyqecclang.algorithms.qode.lchs", "lchs_qode"),
        "ContourPlan": ("pyqecclang.algorithms.qode.cbmd", "ContourPlan"),
        "cbmd_qode": ("pyqecclang.algorithms.qode.cbmd", "cbmd_qode"),
        "cbmd_function": ("pyqecclang.algorithms.qode.cbmd", "cbmd_function"),
        "taylor_hamiltonian": ("pyqecclang.algorithms.common.hamiltonian", "taylor_hamiltonian"),
        "qft": ("pyqecclang.algorithms.common.fourier", "qft_with_work"),
        "fourier_momentum": ("pyqecclang.algorithms.qode.schrodingerization", "fourier_momentum"),
        "SchrodingerPlan": ("pyqecclang.algorithms.qode.schrodingerization", "SchrodingerPlan"),
        "schrodinger_qode": ("pyqecclang.algorithms.qode.schrodingerization", "schrodinger_qode"),
        "PolynomialODE": ("pyqecclang.algorithms.qnlss.carleman", "PolynomialODE"),
        "_carleman_term": ("pyqecclang.algorithms.qnlss.carleman", "_carleman_term"),
        "carleman_lift": ("pyqecclang.algorithms.qnlss.carleman", "carleman_lift"),
        "carleman_initial": ("pyqecclang.algorithms.qnlss.carleman", "carleman_initial"),
        "carleman_qode": ("pyqecclang.algorithms.qnlss.carleman", "carleman_qode"),
        "linear_qode": ("pyqecclang.algorithms.qode.ode", "linear_qode"),
        "PDEInput": ("pyqecclang.algorithms.qpde.pde", "PDEInput"),
        "qpde_solver": ("pyqecclang.algorithms.qpde.pde", "qpde_solver"),
    },
    "pyqecclang.algorithms.solvers": {
        "extend_initial": ("pyqecclang.algorithms.common.state_preparation", "extend_initial"),
        "select_subspace": ("pyqecclang.algorithms.common.state_preparation", "select_subspace"),
        "apply_be_to_state": ("pyqecclang.algorithms.common.state_preparation", "apply_be_to_state"),
        "make_euler_history_qode": ("pyqecclang.algorithms.qode.ode", "make_euler_history_qode"),
        "DiscretePDE": ("pyqecclang.algorithms.qpde.pde", "DiscretePDE"),
        "make_qpde": ("pyqecclang.algorithms.qpde.pde", "make_qpde"),
        "make_lchs_qode": ("pyqecclang.algorithms.qode.legacy", "make_lchs_qode"),
        "make_schrodingerisation_qode": (
            "pyqecclang.algorithms.qode.legacy",
            "make_schrodingerisation_qode",
        ),
        "trotter_hamsim": ("pyqecclang.algorithms.common.hamiltonian", "trotter_hamsim"),
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
