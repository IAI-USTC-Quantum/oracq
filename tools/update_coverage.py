"""把旧语言用例逐项映射到新范式案例；不声称逐字翻译或 golden 等价。"""

import json
from pathlib import Path

from oracq.applications.catalog import CASES

ROOT = Path(__file__).resolve().parents[1]


def positive(case):
    category, filename = case.split("/", 1)
    name = filename.split("/")[-1].removesuffix(".qec")
    if category == "00-primitives":
        return {
            "bell-pair": ["bell"],
            "ghz": ["ghz"],
            "measure-reset": ["measure_reset"],
            "rot-angle": ["qsvt"],
            "adjoint-block": ["oaa"],
            "control-nested": ["costa_gate"],
            "dagger-self-inverse": ["oaa"],
        }.get(name, ["register_views"])
    if category == "01-registers":
        if "qram" in name or "array" in name:
            return ["batch_qram", "banked_qram"]
        if "fixed" in name:
            return ["arithmetic"]
        if "ancilla" in name:
            return ["costa_gate"]
        return ["register_views"]
    if category == "02-generics-const":
        return ["python_generators", "be_algebra"]
    if category == "03-oracle":
        if "sparse" in name:
            return ["sparse_gate", "sparse_qram"]
        if "state-prep" in name:
            return ["stateprep_gate", "stateprep_qram"]
        if "qram" in name:
            return ["dj_qram"]
        return ["grover_gate", "grover_qram"]
    if category == "04-block-encoding":
        if "sparse" in name:
            return ["sparse_qram"]
        if "walk" in name:
            return ["qpe"]
        return ["be_algebra"]
    if category == "05-qsvt-qpe":
        if "filter" in name:
            return ["costa_gate", "costa_qram"]
        if "qpe" in name or "walk" in name:
            return ["qpe"]
        if "oaa" in name:
            return ["oaa"]
        return ["qsvt", "python_generators"]
    if category == "06-protocols":
        if "trotter" in name:
            return ["trotter_hamsim"]
        if "order" in name:
            return ["dj_gate", "dj_qram"]
        return ["costa_qram", "qham_qode"]
    if category == "07-scientific":
        return {
            "heat-equation": ["heat_qode"],
            "schrodingerization": ["schrodingerisation"],
            "carleman-step": ["carleman_step"],
            "lchs-ode": ["lchs"],
            "poisson-qlss": ["poisson_qlss"],
            "cfd-implicit-step": ["qfvm_qram"],
            "stateprep-vs-oracle": ["stateprep_gate", "costa_gate"],
        }[name]
    if category == "08-modules":
        return ["dj_gate", "dj_qram", "qham_qpde"]
    raise ValueError(case)


def main():
    inventory = json.loads((ROOT / "docs/archive/case-inventory.json").read_text())
    changes = {
        "ancilla-clean": "工作区显式进入接口；不继承旧版 clean_ancilla 静态证明。",
        "measure-reset": "由显式宿主 ReadoutAction 生成末端测量和重置；不混入酉 oracle 主体。",
        "reinterpret-fixed": "定点解释作为存储位模式和库元数据，不内建数值精度证明。",
        "fixed-point-add": "可逆算术接口与小型查表实现；RNE/溢出算法语义待下一阶段。",
    }
    negative = {
        "unbound-require-program": "现在允许保存开放 IR，在后端导出时报告缺口。",
        "oracle-takes-operation": "由 Python 高阶生成函数承担，旧语法限制不继承。",
        "partial-application": "允许 Python 绑定和闭包；生成后 IR 不保留 Python callback。",
        "require-in-program": "require 文本语法不迁移，依赖通过 Python 参数和显式绑定表达。",
        "require-mid-body": "require 文本语法不迁移。",
        "isometry-nonzero-input": "零输入是接口契约，本阶段不证明零态。",
        "isometry-adjoint-cap": "等距角色可以提供 unitary 扩张；无声明逆能力时拒绝逆调用。",
        "ancilla-not-clean": "所有工作区显式持有；不继承未实现的自动复净判定。",
        "fused-source-use": "逻辑视图不进行持久冻结；同次调用的重叠和控制修改仍拒绝。",
        "measure-in-if": "量子结果读出在宿主层，RIR 没有动态经典分支。",
    }
    rows = []
    for original in inventory["cases"]:
        row = dict(original)
        key = row["case"].split("/")[-1].removesuffix(".qec")
        if row["kind"] == "positive":
            row["examples"] = positive(row["case"])
            row["status"] = "paradigm_mapped"
            row["note"] = changes.get(
                key, "映射到可生成范式案例；未执行旧 .qec 或比较旧 CLIR golden。"
            )
        else:
            row["examples"] = []
            row["status"] = (
                "intentional_design_change" if key in negative else "structural_rejection"
            )
            row["note"] = negative.get(key, "对应 RIR 的签名、类型、别名、常量或调用图验证。")
            row["evidence"] = ["tests/core/test_language.py", "tests/core/test_open_ir.py"]
        assert all(example in CASES for example in row["examples"])
        rows.append(row)
    data = {
        "mode": "paradigm_coverage_not_golden_equivalence",
        "source": inventory["source"],
        "positive": sum(r["kind"] == "positive" for r in rows),
        "negative": sum(r["kind"] == "negative" for r in rows),
        "cases": rows,
    }
    (ROOT / "docs/archive/coverage.json").write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
    lines = [
        "# 旧案例与新范式覆盖矩阵",
        "",
        "本表覆盖 61 个正例和 16 个负例。映射证明有对应的表达和组装路径，不表示旧源码逐字迁移、旧 golden 等价或算法正确性认证。",
        "",
        "| 旧用例 | 状态 | 新例子/证据 | 说明 |",
        "|---|---|---|---|",
    ]
    for row in rows:
        evidence = ", ".join(row["examples"] or row["evidence"])
        lines.append(f"| {row['case']} | {row['status']} | {evidence} | {row['note']} |")
    lines += [
        "",
        "六组参考负载分别由 qfvm_gate/qram、be_algebra、arithmetic、costa/sparse、oracle 目录与 qham_qode/qpde 覆盖。"
        "Roe 物理核、一般高阶 QHAM、严格 QSVT 相位与 PDE 收敛证明保留为下一阶段工作，已定义对应的开放接口。",
        "",
    ]
    (ROOT / "docs/archive/coverage.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
