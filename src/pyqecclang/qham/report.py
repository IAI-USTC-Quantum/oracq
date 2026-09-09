"""导出可审阅的自动推导与 QODE input manifest。"""

import json
from dataclasses import asdict


def input_manifest(plan, state_width):
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


def row_description(plan, block, state_width, eta):
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


def export_derivation(plan, directory, *, state_width=2, eta=-1.0, max_blocks=256, row=None):
    import math
    from pathlib import Path

    from ..ir import ValidationError

    if not math.isfinite(complex(eta).real) or not math.isfinite(complex(eta).imag):
        raise ValidationError("eta 必须有限")
    path = Path(directory)
    path.mkdir(parents=True, exist_ok=True)
    (path / "pde.json").write_text(plan.pde.dumps())
    (path / "qcl-plan.json").write_text(plan.dumps())
    manifest = input_manifest(plan, state_width)
    manifest["eta"] = [complex(eta).real, complex(eta).imag]
    if row is not None:
        blocks = (row,)
        manifest["rows_materialized"] = "single requested row"
    elif plan.block_count <= max_blocks:
        blocks = tuple(plan.blocks())
        manifest["rows_materialized"] = True
    else:
        blocks = ()
        manifest["rows_materialized"] = False
        manifest["note"] = "plan remains valid and queryable; explicit enumeration budget exceeded"
    if row is not None and not plan.contains(row):
        raise ValidationError("请求的行不在 QCL 闭包")
    rows = [row_description(plan, b, state_width, eta) for b in blocks]
    (path / "rows.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2) + "\n")
    (path / "qode-input.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    text = [
        f"# 自动 QHAM 推导：{plan.pde.label}",
        "",
        f"HAM 阶数 m={plan.order}，非线性最高次数 D={plan.pde.degree}，最大张量秩 {plan.max_rank}。",
        f"函数块数 {plan.block_count}，原始提升维数 {manifest['raw_lifted_dimension']}。",
        "",
        "同伦递推：Ui' = L Ui - eta sum_l (1+eta)^(i-1-l) C_l。",
        "物理输出是 u_sum；one 分量若存在则满足 one'=0，one(0)=1。",
        "",
        "| 行块 | 源块 | 线性端口 | 作用位置 | 系数 |",
        "|---|---|---|---|---|",
    ]
    for item in rows:
        for term in item["terms"]:
            text.append(
                f"| {item['row']} | {term['column']} | {term['operator']} ({term['arity']}→1) | {term['position']} | {term['formula']} |"
            )
    text += ["", "以上是对截断 HAM 的自动线性化，收敛与量子求解精度另行验证。"]
    (path / "derivation.md").write_text("\n".join(text) + "\n")
    return manifest
