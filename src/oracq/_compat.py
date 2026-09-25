"Central compatibility table for 0.7 import paths; in-repo code uses only the canonical paths."

from __future__ import annotations

import importlib
import sys
from types import ModuleType

ALIASES = {
    "oracq.ir": "oracq.infrastructure.ir",
    "oracq.builder": "oracq.infrastructure.builder",
    "oracq.validation": "oracq.infrastructure.validation",
    "oracq.serialization": "oracq.infrastructure.serialization",
    "oracq.linking": "oracq.infrastructure.linking",
    "oracq.execution": "oracq.infrastructure.execution",
    "oracq.native": "oracq.infrastructure.native",
    "oracq.layout": "oracq.infrastructure.layout",
    "oracq.readout": "oracq.infrastructure.readout",
    "oracq.backends.originir": "oracq.infrastructure.backends.originir",
    "oracq.backends": "oracq.infrastructure.backends",
    "oracq.backends.pysparq": "oracq.infrastructure.backends.pysparq",
    "oracq.backends.basis": "oracq.infrastructure.backends.basis",
    "oracq.mathfunc.graph": "oracq.infrastructure.mathfunc.graph",
    "oracq.mathfunc.frontend": "oracq.infrastructure.mathfunc.frontend",
    "oracq.mathfunc.numeric": "oracq.infrastructure.mathfunc.numeric",
    "oracq.mathfunc.lowering": "oracq.infrastructure.mathfunc.lowering",
    "oracq.mathfunc": "oracq.infrastructure.mathfunc",
    "oracq.mathfunc.roe_formulas": "oracq.applications.roe_formulas",
    "oracq.library": "oracq.algorithms.input_model.operators",
    "oracq.combinators": "oracq.algorithms.input_model.block_encoding",
    "oracq.contracts": "oracq.algorithms.input_model.contracts",
    "oracq.oracles": "oracq.algorithms.input_model.oracles",
    "oracq.arithmetic": "oracq.algorithms.common.arithmetic",
    "oracq.sparse_models": "oracq.algorithms.input_model.sparse",
    "oracq.access": "oracq.algorithms.input_model.sparse",
    "oracq.qlss": "oracq.algorithms.qlss.qlss",
    "oracq.algorithms.costa": "oracq.algorithms.qlss.qlss",
    "oracq.algorithms.cks": "oracq.algorithms.qlss.qlss",
    "oracq.qode": "oracq.algorithms.qode.ode",
    "oracq.qfvm": "oracq.applications.qfvm",
    "oracq.qfvm_sparse": "oracq.applications.qfvm",
    "oracq.flow_data": "oracq.applications.flow_data",
    "oracq.roe": "oracq.applications.roe",
    "oracq.workloads": "oracq.applications.catalog",
    "oracq.qham.pde": "oracq.applications.qham.pde",
    "oracq.qham.linearization": "oracq.applications.qham.linearization",
    "oracq.qham.reference": "oracq.applications.qham.reference",
    "oracq.qham.quantum": "oracq.algorithms.input_model.qham",
    "oracq.qham.examples": "oracq.applications.qham.examples",
    "oracq.qham.report": "oracq.applications.qham.report",
    "oracq.qham.stencils": "oracq.applications.qham.stencils",
}
SPLIT_EXPORTS = {
    "oracq.algorithms.elementary": {
        "deutsch_jozsa": ("oracq.algorithms.basics.oracle_algorithms", "deutsch_jozsa"),
        "grover": ("oracq.algorithms.common.search", "grover"),
        "phase_from_database": ("oracq.algorithms.common.search", "phase_from_database"),
        "qft": ("oracq.algorithms.common.fourier", "qft"),
        "phase_estimation": ("oracq.algorithms.common.estimation", "phase_estimation"),
        "qubitization_walk": ("oracq.algorithms.common.transforms", "qubitization_walk"),
        "qsvt_sequence": ("oracq.algorithms.common.transforms", "qsvt_sequence"),
        "oblivious_amplification": ("oracq.algorithms.common.transforms", "oblivious_amplification"),
    },
    "oracq.algorithms.differential": {
        "tagged": ("oracq.algorithms.qode._dynamics", "tagged"),
        "_lcu_dynamics": ("oracq.algorithms.qode._dynamics", "_lcu_dynamics"),
        "HermitianParts": ("oracq.algorithms.qode.ode_models", "HermitianParts"),
        "LinearODE": ("oracq.algorithms.qode.ode_models", "LinearODE"),
        "QuadraturePlan": ("oracq.algorithms.qode.lchs", "QuadraturePlan"),
        "lchs_qode": ("oracq.algorithms.qode.lchs", "lchs_qode"),
        "ContourPlan": ("oracq.algorithms.qode.cbmd", "ContourPlan"),
        "cbmd_qode": ("oracq.algorithms.qode.cbmd", "cbmd_qode"),
        "cbmd_function": ("oracq.algorithms.qode.cbmd", "cbmd_function"),
        "taylor_hamiltonian": ("oracq.algorithms.common.hamiltonian", "taylor_hamiltonian"),
        "qft": ("oracq.algorithms.common.fourier", "qft_with_work"),
        "fourier_momentum": ("oracq.algorithms.qode.schrodingerization", "fourier_momentum"),
        "SchrodingerPlan": ("oracq.algorithms.qode.schrodingerization", "SchrodingerPlan"),
        "schrodinger_qode": ("oracq.algorithms.qode.schrodingerization", "schrodinger_qode"),
        "PolynomialODE": ("oracq.algorithms.qnlss.carleman", "PolynomialODE"),
        "_carleman_term": ("oracq.algorithms.qnlss.carleman", "_carleman_term"),
        "carleman_lift": ("oracq.algorithms.qnlss.carleman", "carleman_lift"),
        "carleman_initial": ("oracq.algorithms.qnlss.carleman", "carleman_initial"),
        "carleman_qode": ("oracq.algorithms.qnlss.carleman", "carleman_qode"),
        "linear_qode": ("oracq.algorithms.qode.ode", "linear_qode"),
        "PDEInput": ("oracq.algorithms.qpde.pde", "PDEInput"),
        "qpde_solver": ("oracq.algorithms.qpde.pde", "qpde_solver"),
    },
    "oracq.algorithms.solvers": {
        "extend_initial": ("oracq.algorithms.common.state_preparation", "extend_initial"),
        "select_subspace": ("oracq.algorithms.common.state_preparation", "select_subspace"),
        "apply_be_to_state": ("oracq.algorithms.common.state_preparation", "apply_be_to_state"),
        "make_euler_history_qode": ("oracq.algorithms.qode.ode", "make_euler_history_qode"),
        "DiscretePDE": ("oracq.algorithms.qpde.pde", "DiscretePDE"),
        "make_qpde": ("oracq.algorithms.qpde.pde", "make_qpde"),
        "make_lchs_qode": ("oracq.algorithms.qode.legacy", "make_lchs_qode"),
        "make_schrodingerisation_qode": (
            "oracq.algorithms.qode.legacy",
            "make_schrodingerisation_qode",
        ),
        "trotter_hamsim": ("oracq.algorithms.common.hamiltonian", "trotter_hamsim"),
    },
}


def install() -> None:
    """Install the legacy 0.7 import paths as canonical modules and attach them to their parent packages.

    Legacy paths in ``ALIASES`` directly reuse the corresponding canonical module
    objects; legacy packages in ``SPLIT_EXPORTS`` are synthesized as forwarding
    modules from their canonical owner modules per the export manifest.
    """
    for old, new in ALIASES.items():
        module = importlib.import_module(new)
        sys.modules[old] = module
    for old, exports in SPLIT_EXPORTS.items():
        module = ModuleType(old)
        module.__package__ = old.rpartition(".")[0]
        for name, (owner, symbol) in exports.items():
            setattr(module, name, getattr(importlib.import_module(owner), symbol))
        # typeshed does not declare __all__ as assignable on ModuleType; dynamically synthesized modules need it.
        module.__all__ = [n for n in exports if not n.startswith("_")]  # type: ignore[attr-defined]
        sys.modules[old] = module
    for old in (*ALIASES, *SPLIT_EXPORTS):
        parent, _, leaf = old.rpartition(".")
        # owner is bound to module path strings in the unpacking loop above and carries module objects here.
        owner = sys.modules.get(parent) or importlib.import_module(parent)  # type: ignore[assignment]
        setattr(owner, leaf, sys.modules[old])
