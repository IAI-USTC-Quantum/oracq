"""Generate Sphinx API pages from the canonical source tree and a top-level
overview page from the root package __all__.

Two language trees are generated: ``docs/api`` (English titles) and
``docs/zh/api`` (Chinese titles). Both document the same modules; page bodies
come from the English docstrings in the source.

Legacy import paths are kept for compatibility only and never appear in the
API documentation.
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TITLES = {
    "oracle_algorithms": "Oracle query algorithms",
    "fourier": "Fourier transforms and arithmetic",
    "search": "Search and amplitude amplification",
    "estimation": "Phase, amplitude, and overlap estimation",
    "variational": "Variational algorithm circuits",
    "walks": "Quantum walks",
    "graph_walks": "Graph-walk search",
    "number_theory": "Modular multiplication and order finding",
    "error_correction": "Repetition codes and error recovery",
    "hamiltonian": "Hamiltonian evolution",
    "qlss": "Quantum linear systems",
    "vtaa_cks": "VTAA-CKS variable-time linear-system solver",
    "newton": "Quantum Newton method",
    "ode": "QODE assembly interface",
    "lchs": "LCHS",
    "schrodingerization": "Schrödingerization",
    "carleman": "Carleman linearization",
    "cbmd": "CBMD",
    "sde": "SDE/Fokker–Planck input models",
    "qham": "QHAM",
    "qsvt": "QSVT standard transforms",
    "spectral": "Spectral input/output primitives and operator-level optimization",
    "spectral_synthesis": "Operator-level synthesis optimization of spectral circuits",
    "lowrank": "Chemistry low-rank decomposition block encodings",
    "dqi": "DQI decoding quantum interference optimization",
    "recommendation": "KP quantum recommendation system",
    "gradient": "Quantum gradient estimation",
    "integration": "Quantum summation and integration",
    "qpca": "QPCA quantum principal component analysis",
    "qsdp": "Quantum semidefinite programming framework",
    "qcnn": "Quantum convolutional neural networks",
    "qcnn_layer": "Quantum building blocks of QCNNs",
    "density": "Density-matrix input models and Gibbs states",
    "data_loading": "Select-Swap QROM data loading",
    "qdata": "Quantum data structures (qsample and sample-and-query)",
    "contracts": "Algorithm contracts and reports",
    "interfaces": "Algorithm structure protocols",
    "operators": "Operator wrappers and basic composition",
    "block_encoding": "Block-encoding composition",
    "oracles": "Oracle declarations and implementations",
    "arithmetic": "Reversible arithmetic",
    "prepare_select": "PREPARE-SELECT decomposition",
    "sparse": "Sparse-access adapters",
    "pde": "PDE models and adapters",
    "ode_models": "Linear ODE input models",
    "state_preparation": "State-preparation composition",
    "transforms": "Matrix transform sequences",
    "ir": "RIR objects",
    "builder": "Module builder",
    "serialization": "RIR serialization",
    "validation": "Structural validation",
    "linking": "Binding and capability analysis",
    "execution": "Register reference executor",
    "native": "Native implementation registry",
    "layout": "Register layout",
    "readout": "Host readout",
    "qmem": "QRAM pointer-style read/write",
    "qram_schema": "QRAM YAML memory format",
    "estimate": "Resource estimation",
    "originir": "OriginIR-ext backend",
    "pysparq": "PySparQ backend",
    "basis": "Toffoli / U3 / CZ lowering",
    "quantikz": "Quantikz circuit export",
    "strict": "Strict netlist export",
    "mathfunc": "Math-function compilation entry",
    "frontend": "Python math-function frontend",
    "graph": "MIR objects",
    "lowering": "Math-function lowering",
    "numeric": "Math kernel generation",
    "qfvm": "QFVM application",
    "qfvm_qmem": "QFVM direct QMem data path",
    "flow_data": "Flow fields and QRAM data",
    "roe": "Roe matrix-element generation",
    "roe_formulas": "Roe classical math formulas",
    "gallery": "Algorithm gallery",
    "catalog": "Reference workload catalog",
    "linearization": "QHAM finite closure",
    "reference": "Spatial discretization and classical reference",
    "stencils": "Structured finite-difference ports",
    "report": "QHAM derivation report",
    "examples": "QHAM PDE examples",
}
TITLES_ZH = {
    "oracle_algorithms": "Oracle 查询算法",
    "fourier": "Fourier 变换与算术",
    "search": "搜索与振幅放大",
    "estimation": "相位、振幅与重叠估计",
    "variational": "变分算法电路",
    "walks": "量子行走",
    "graph_walks": "图行走搜索",
    "number_theory": "模乘与求阶",
    "error_correction": "重复码与错误恢复",
    "hamiltonian": "Hamiltonian 演化",
    "qlss": "量子线性系统",
    "vtaa_cks": "VTAA-CKS 变时线性系统求解器",
    "newton": "量子牛顿法",
    "ode": "QODE 组装接口",
    "lchs": "LCHS",
    "schrodingerization": "Schrödingerization",
    "carleman": "Carleman 线性化",
    "cbmd": "CBMD",
    "sde": "SDE/Fokker–Planck 输入模型",
    "qham": "QHAM",
    "qsvt": "QSVT 标准变换",
    "spectral": "谱输入输出原语与算子级优化",
    "spectral_synthesis": "谱线路的算子级合成优化",
    "lowrank": "化学低秩分解块编码",
    "dqi": "DQI 解码量子干涉优化",
    "recommendation": "KP 量子推荐系统",
    "gradient": "量子梯度估计",
    "integration": "量子求和与积分",
    "qpca": "QPCA 量子主成分分析",
    "qsdp": "量子半定规划框架",
    "qcnn": "量子卷积神经网络",
    "qcnn_layer": "量子卷积神经网络的量子构件",
    "density": "密度矩阵输入模型与 Gibbs 态",
    "data_loading": "Select-Swap QROM 数据加载",
    "qdata": "量子数据结构（qsample 与 sample-and-query）",
    "contracts": "算法契约与报告",
    "interfaces": "算法结构协议",
    "operators": "算子包装与基本组合",
    "block_encoding": "Block encoding 组合",
    "oracles": "Oracle 声明与实现",
    "arithmetic": "可逆算术",
    "prepare_select": "PREPARE-SELECT 分解",
    "sparse": "稀疏访问适配",
    "pde": "PDE 模型与适配",
    "ode_models": "线性 ODE 输入模型",
    "state_preparation": "态制备组合",
    "transforms": "矩阵变换序列",
    "ir": "RIR 对象",
    "builder": "模块构造器",
    "serialization": "RIR 序列化",
    "validation": "结构验证",
    "linking": "绑定与能力分析",
    "execution": "寄存器参考执行器",
    "native": "原生实现注册",
    "layout": "寄存器布局",
    "readout": "宿主读出",
    "qmem": "QRAM 指针式读写",
    "qram_schema": "QRAM YAML 内存格式",
    "estimate": "资源估计",
    "originir": "OriginIR-ext 后端",
    "pysparq": "PySparQ 后端",
    "basis": "Toffoli / U3 / CZ 降低",
    "quantikz": "Quantikz 线路导出",
    "strict": "严格网表导出",
    "mathfunc": "数学函数编译入口",
    "frontend": "Python 数学函数前端",
    "graph": "MIR 对象",
    "lowering": "数学函数降低",
    "numeric": "数学核生成",
    "qfvm": "QFVM 应用",
    "qfvm_qmem": "QFVM 的 QMem 直连数据路径",
    "flow_data": "流场与 QRAM 数据",
    "roe": "Roe 矩阵元生成",
    "roe_formulas": "Roe 经典数学公式",
    "gallery": "算法展示目录",
    "catalog": "参考工作负载目录",
    "linearization": "QHAM 有限闭包",
    "reference": "空间离散与经典参考",
    "stencils": "结构化差分端口",
    "report": "QHAM 推导报告",
    "examples": "QHAM PDE 示例",
}
GROUPS = ("infrastructure", "algorithms", "applications")
GROUP_TITLES = {
    "infrastructure": "Infrastructure API",
    "algorithms": "Algorithms API",
    "applications": "Domain Applications API",
}
GROUP_TITLES_ZH = {
    "infrastructure": "基础设施 API",
    "algorithms": "算法 API",
    "applications": "领域应用 API",
}
STRINGS = {
    "en": {
        "titles": TITLES,
        "group_titles": GROUP_TITLES,
        "toplevel_title": "oracq package overview",
        "toplevel_intro": [
            "The root package ``oracq`` re-exports the public API ({count} names).",
            "Names are grouped by their defining module; module titles link to the "
            "corresponding API page and members link to their full descriptions there.",
        ],
        "joiner": ", ",
        "dash": "—",
        "api_title": "API Reference",
        "api_intro": (
            "API pages are generated from the canonical source paths. Legacy paths "
            "are kept for import compatibility only and are not listed twice."
        ),
    },
    "zh": {
        "titles": TITLES_ZH,
        "group_titles": GROUP_TITLES_ZH,
        "toplevel_title": "oracq 包总览",
        "toplevel_intro": [
            "根包 ``oracq`` 汇总导出公开 API（共 {count} 个名字）。",
            "名字按定义模块分组；模块标题链接到对应 API 页，成员链接到模块页内的完整说明。",
        ],
        "joiner": "、",
        "dash": "——",
        "api_title": "API 参考",
        "api_intro": "API 从规范源码路径生成。旧路径只保留导入兼容，不重复列出。",
    },
}


def _underline(title: str, char: str) -> str:
    return char * max(12, len(title) * 2)


def write_toplevel(docs: Path, modules: dict[str, dict], strings: dict) -> None:
    """Generate the top-level overview page from the root package __all__,
    grouping names by their defining module."""
    sys.path.insert(0, str(ROOT / "src"))
    import importlib

    import oracq

    grouped: dict[str, dict[str, list[str]]] = {group: {} for group in GROUPS}
    for name in oracq.__all__:
        obj = getattr(oracq, name)
        module = getattr(obj, "__module__", None)
        if module not in modules:
            module = None
        if module is None:
            for candidate in modules:
                attribute = getattr(importlib.import_module(candidate), name, None)
                if attribute is obj:
                    module = candidate
                    break
        if module is None:
            raise SystemExit(f"{name} from the root __all__ was not located in any canonical module")
        grouped[module.split(".")[1]].setdefault(module, []).append(name)

    lines = [
        strings["toplevel_title"],
        _underline(strings["toplevel_title"], "="),
        "",
    ]
    lines += [line.format(count=len(oracq.__all__)) for line in strings["toplevel_intro"]]
    lines.append("")
    for group in GROUPS:
        modules_in_group = grouped[group]
        if not modules_in_group:
            continue
        title = strings["group_titles"][group]
        lines += [title, _underline(title, "-"), ""]
        for module in sorted(modules_in_group):
            info = modules[module]
            names = modules_in_group[module]
            refs = strings["joiner"].join(
                f":obj:`{name} <{module}.{name}>`" for name in names
            )
            lines += [
                f":doc:`{info['title']} <{info['doc']}>`",
                f"    ``{module}`` {strings['dash']} {refs}",
                "",
            ]
    (docs / "toplevel.rst").write_text("\n".join(lines), encoding="utf-8")


def generate_language(docs: Path, lang: str) -> int:
    strings = STRINGS[lang]
    titles = strings["titles"]
    docs.mkdir(parents=True, exist_ok=True)
    groups = {}
    modules: dict[str, dict] = {}
    for group in GROUPS:
        pages = []
        for source in sorted((ROOT / "src/oracq" / group).rglob("*.py")):
            if source.name.startswith("_") or source.stem == "legacy":
                if source.name != "__init__.py" or source.parent.name != "mathfunc":
                    continue
            relative = source.relative_to(ROOT / "src/oracq").with_suffix("")
            parts = relative.parts[:-1] if source.name == "__init__.py" else relative.parts
            module = "oracq." + ".".join(parts)
            target = docs.joinpath(*parts).with_suffix(".rst")
            target.parent.mkdir(parents=True, exist_ok=True)
            title = titles.get(parts[-1], parts[-1])
            modules[module] = {"title": title, "doc": "/".join(parts)}
            target.write_text(
                title
                + "\n"
                + _underline(title, "=")
                + "\n\n"
                + "``"
                + module
                + "``\n\n"
                + ".. automodule:: "
                + module
                + "\n   :members:\n   :undoc-members:\n   :show-inheritance:\n",
                encoding="utf-8",
            )
            pages.append(target.relative_to(docs / group).with_suffix("").as_posix())
        groups[group] = pages
        title = strings["group_titles"][group]
        (docs / group / "index.rst").write_text(
            title
            + "\n"
            + _underline(title, "=")
            + "\n\n.. toctree::\n   :maxdepth: 1\n\n"
            + "".join("   " + page + "\n" for page in pages),
            encoding="utf-8",
        )
    write_toplevel(docs, modules, strings)
    (docs / "index.rst").write_text(
        strings["api_title"]
        + "\n"
        + "=" * max(12, len(strings["api_title"]) * 2)
        + "\n\n"
        + strings["api_intro"]
        + "\n\n.. toctree::\n   :maxdepth: 2\n\n"
        "   toplevel\n   infrastructure/index\n   algorithms/index\n   applications/index\n",
        encoding="utf-8",
    )
    return sum(map(len, groups.values()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lang", choices=["en", "zh", "all"], default="all")
    args = parser.parse_args()
    langs = ["en", "zh"] if args.lang == "all" else [args.lang]
    total = 0
    for lang in langs:
        docs = ROOT / "docs/api" if lang == "en" else ROOT / "docs/zh/api"
        count = generate_language(docs, lang)
        print(f"Generated {count} API pages + toplevel ({lang}) -> {docs.relative_to(ROOT)}")
        total += count


if __name__ == "__main__":
    main()
