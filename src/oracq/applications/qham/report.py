"Export reviewable automated derivations and the QODE input manifest."

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import cast

from oracq.applications.qham.linearization import Block, QHAMPlan


def input_manifest(plan: QHAMPlan, state_width: int) -> dict[str, object]:
    """Build the input manifest handed to the QODE solver side.

    Args:
        plan: ``QHAMPlan`` closure plan.
        state_width: Bit width of a single-component state.

    Returns:
        dict: ``plan.summary()`` extended with the state width and dimension,
        the generator target width, the total lifted dimension, per-port input
        and output widths and block encoding requirements, plus constraint
        notes on initial-state preparation, time dependence, and sparse
        access.
    """
    dimension = 1 << state_width
    raw = plan.raw_dimension(dimension)
    return {
        **plan.summary(),
        "state_width": state_width,
        "state_dimension": dimension,
        "generator_target_width": (raw - 1).bit_length(),
        "raw_lifted_dimension": raw,
        "dynamics": "dY/dt=G*Y; homogeneous constant included iff forcing exists",
        "physical_window": [0, dimension],
        "basis_order": "factor 0 occupies least significant coordinate bits",
        "ports": [
            {
                "name": p.name,
                "input_width": p.arity * state_width,
                "output_width": state_width,
                "padded_be_width": max(1, p.arity) * state_width,
                "required": "block encoding operation, alpha and ancilla width",
            }
            for p in plan.pde.ports
        ],
        "initial_input": "reversible preparation of u_in/||u_in|| and classical ||u_in||",
        "time_dependence": "autonomous bindings; time-dependent families require a different QODE adapter",
        "solver_requirements": "arbitrary linear generator, or explicit dissipative shift before LCHS/CBMD",
        "sparse_model": "requires base-operator sparse access; never inferred from a generic BE",
    }


def row_description(
    plan: QHAMPlan, block: Block, state_width: int, eta: complex
) -> dict[str, object]:
    """Build the derivation description of a single row block.

    Args:
        plan: ``QHAMPlan`` closure plan.
        block: Row block within the closure.
        state_width: Bit width of a single-component state.
        eta: Homotopy parameter.

    Returns:
        dict: Contains the block label, tensor word, block offset, and block
        dimension, plus for each linear edge the source block, port, action
        position, weight formula, and value at ``eta`` (real and imaginary
        parts).
    """
    dimension = 1 << state_width
    return {
        "row": block.label,
        "orders": block.orders,
        "offset": plan.offset(block, dimension),
        "dimension": dimension**block.rank,
        "terms": [
            {
                "column": term.column.label,
                "column_orders": term.column.orders,
                "column_offset": plan.offset(term.column, dimension),
                "operator": term.operator,
                "arity": term.arity,
                "position": term.position,
                "weight": asdict(term.weight),
                "formula": term.weight.formula(),
                "evaluated_weight": [
                    complex(term.weight.evaluate(eta)).real,
                    complex(term.weight.evaluate(eta)).imag,
                ],
            }
            for term in plan.row_terms(block)
        ],
    }


def export_derivation(
    plan: QHAMPlan,
    directory: str | Path,
    *,
    state_width: int = 2,
    eta: complex = -1.0,
    max_blocks: int = 256,
    row: Block | None = None,
) -> dict[str, object]:
    """Export the automated derivation results as reviewable JSON, Markdown, and a QODE manifest.

    Writes ``pde.json``, ``qcl-plan.json``, ``rows.json``,
    ``qode-input.json``, and ``derivation.md`` under ``directory``; the
    directory is created recursively when absent.

    Args:
        plan: ``QHAMPlan`` closure plan.
        directory: Export directory.
        state_width: Bit width of a single-component state.
        eta: Homotopy parameter; must be a finite complex number.
        max_blocks: Budget for explicitly enumerating row blocks; when the
            closure block count exceeds it, rows are not enumerated and the
            manifest keeps the plan summary with a note that the budget was
            exceeded.
        row: When given, export only this row block; otherwise export all row
            blocks.

    Returns:
        dict: The manifest written to ``qode-input.json`` (including ``eta``).

    Raises:
        ValidationError: ``eta`` is not finite, or ``row`` is not in the QCL
            closure.
    """
    import math
    from pathlib import Path

    from oracq.infrastructure.ir import ValidationError

    if not math.isfinite(complex(eta).real) or not math.isfinite(complex(eta).imag):
        raise ValidationError("eta must be finite")
    path = Path(directory)
    path.mkdir(parents=True, exist_ok=True)
    (path / "pde.json").write_text(plan.pde.dumps())
    (path / "qcl-plan.json").write_text(plan.dumps())
    manifest = input_manifest(plan, state_width)
    manifest["eta"] = [complex(eta).real, complex(eta).imag]
    if row is not None:
        blocks: tuple[Block, ...] = (row,)
        manifest["rows_materialized"] = "single requested row"
    elif plan.block_count <= max_blocks:
        blocks = tuple(plan.blocks())
        manifest["rows_materialized"] = True
    else:
        blocks = ()
        manifest["rows_materialized"] = False
        manifest["note"] = "plan remains valid and queryable; explicit enumeration budget exceeded"
    if row is not None and not plan.contains(row):
        raise ValidationError("the requested row is not in the QCL closure")
    rows = [row_description(plan, b, state_width, eta) for b in blocks]
    (path / "rows.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2) + "\n")
    (path / "qode-input.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    text = [
        f"# Automated QHAM derivation: {plan.pde.label}",
        "",
        f"HAM order m={plan.order}, highest nonlinear degree D={plan.pde.degree}, maximum tensor rank {plan.max_rank}.",
        f"Block count {plan.block_count}, raw lifted dimension {manifest['raw_lifted_dimension']}.",
        "",
        "Homotopy recursion: Ui' = L Ui - eta sum_l (1+eta)^(i-1-l) C_l.",
        "The physical output is u_sum; the one component, when present, satisfies one'=0, one(0)=1.",
        "",
        "| Row block | Source block | Linear port | Position | Coefficient |",
        "|---|---|---|---|---|",
    ]
    for item in rows:
        # rows is built by row_description; "terms" is always a list of per-row term dictionaries.
        for term in cast(list[dict[str, object]], item["terms"]):
            text.append(
                f"| {item['row']} | {term['column']} | {term['operator']} ({term['arity']}→1) | {term['position']} | {term['formula']} |"
            )
    text += ["", "This is the automated linearization of the truncated HAM; convergence and quantum solver accuracy are validated separately."]
    (path / "derivation.md").write_text("\n".join(text) + "\n")
    return manifest
