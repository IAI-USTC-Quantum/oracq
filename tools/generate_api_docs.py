"""按规范源码目录生成 Sphinx API 页面，并从根包 __all__ 生成顶层总览页。

旧导入路径只保留兼容，不出现在 API 文档中。
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TITLES = {
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
    "infrastructure": "基础设施 API",
    "algorithms": "算法 API",
    "applications": "领域应用 API",
}


def _underline(title: str, char: str) -> str:
    return char * max(12, len(title) * 2)


def write_toplevel(docs: Path, modules: dict[str, dict]) -> None:
    """从根包 __all__ 生成顶层总览页，名字按定义模块归组。"""
    sys.path.insert(0, str(ROOT / "src"))
    import importlib

    import pyqecclang

    grouped: dict[str, dict[str, list[str]]] = {group: {} for group in GROUPS}
    for name in pyqecclang.__all__:
        obj = getattr(pyqecclang, name)
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
            raise SystemExit(f"根 __all__ 中的 {name} 未定位到任何规范模块")
        grouped[module.split(".")[1]].setdefault(module, []).append(name)

    lines = [
        "pyqecclang 包总览",
        _underline("pyqecclang 包总览", "="),
        "",
        f"根包 ``pyqecclang`` 汇总导出公开 API（共 {len(pyqecclang.__all__)} 个名字）。",
        "名字按定义模块分组；模块标题链接到对应 API 页，成员链接到模块页内的完整说明。",
        "",
    ]
    for group in GROUPS:
        modules_in_group = grouped[group]
        if not modules_in_group:
            continue
        title = GROUP_TITLES[group]
        lines += [title, _underline(title, "-"), ""]
        for module in sorted(modules_in_group):
            info = modules[module]
            names = modules_in_group[module]
            refs = "、".join(
                f":obj:`{name} <{module}.{name}>`" for name in names
            )
            lines += [
                f":doc:`{info['title']} <{info['doc']}>`",
                f"    ``{module}`` —— {refs}",
                "",
            ]
    (docs / "toplevel.rst").write_text("\n".join(lines), encoding="utf-8")


def main():
    docs = ROOT / "docs/api"
    docs.mkdir(parents=True, exist_ok=True)
    groups = {}
    modules: dict[str, dict] = {}
    for group in GROUPS:
        pages = []
        for source in sorted((ROOT / "src/pyqecclang" / group).rglob("*.py")):
            if source.name.startswith("_") or source.stem == "legacy":
                if source.name != "__init__.py" or source.parent.name != "mathfunc":
                    continue
            relative = source.relative_to(ROOT / "src/pyqecclang").with_suffix("")
            parts = relative.parts[:-1] if source.name == "__init__.py" else relative.parts
            module = "pyqecclang." + ".".join(parts)
            target = docs.joinpath(*parts).with_suffix(".rst")
            target.parent.mkdir(parents=True, exist_ok=True)
            title = TITLES.get(parts[-1], parts[-1])
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
        title = GROUP_TITLES[group]
        (docs / group / "index.rst").write_text(
            title
            + "\n"
            + _underline(title, "=")
            + "\n\n.. toctree::\n   :maxdepth: 1\n\n"
            + "".join("   " + page + "\n" for page in pages),
            encoding="utf-8",
        )
    write_toplevel(docs, modules)
    (docs / "index.rst").write_text(
        "API 参考\n========\n\n"
        "API 从规范源码路径生成。旧路径只保留导入兼容，不重复列出。\n\n"
        ".. toctree::\n   :maxdepth: 2\n\n"
        "   toplevel\n   infrastructure/index\n   algorithms/index\n   applications/index\n",
        encoding="utf-8",
    )
    print("Generated", sum(map(len, groups.values())), "API pages + toplevel")


if __name__ == "__main__":
    main()
