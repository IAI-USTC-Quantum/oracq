"""按规范源码目录生成 Sphinx API 页面，不导入兼容模块。"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TITLES = {
    "oracle_algorithms": "Oracle 查询算法",
    "fourier": "Fourier 变换与算术",
    "search": "搜索与振幅放大",
    "estimation": "相位、振幅与重叠估计",
    "variational": "变分算法电路",
    "walks": "量子行走",
    "number_theory": "模乘与求阶",
    "error_correction": "重复码与错误恢复",
    "hamiltonian": "Hamiltonian 演化",
    "qlss": "量子线性系统",
    "ode": "QODE 组装接口",
    "lchs": "LCHS",
    "schrodingerization": "Schrödingerization",
    "carleman": "Carleman 线性化",
    "cbmd": "CBMD",
    "qham": "QHAM",
    "contracts": "算法契约与报告",
    "interfaces": "算法结构协议",
    "operators": "算子包装与基本组合",
    "block_encoding": "Block encoding 组合",
    "oracles": "Oracle 声明与实现",
    "arithmetic": "可逆算术",
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
    "originir": "OriginIR-ext 后端",
    "pysparq": "PySparQ 后端",
    "basis": "Toffoli / U3 / CZ 降低",
    "mathfunc": "数学函数编译入口",
    "frontend": "Python 数学函数前端",
    "graph": "MIR 对象",
    "lowering": "数学函数降低",
    "numeric": "数学核生成",
    "qfvm": "QFVM 应用",
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


def main():
    docs = ROOT / "docs/api"
    docs.mkdir(parents=True, exist_ok=True)
    groups = {}
    for group in ("infrastructure", "algorithms", "applications"):
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
            target.write_text(
                title
                + "\n"
                + "=" * max(12, len(title) * 2)
                + "\n\n"
                + "``"
                + module
                + "``\n\n"
                + ".. automodule:: "
                + module
                + "\n   :members:\n   :undoc-members:\n   :show-inheritance:\n",
                encoding="utf-8",
            )
            pages.append(str(target.relative_to(docs / group).with_suffix("")))
        groups[group] = pages
        title = {
            "infrastructure": "基础设施 API",
            "algorithms": "算法 API",
            "applications": "领域应用 API",
        }[group]
        (docs / group / "index.rst").write_text(
            title
            + "\n"
            + "=" * 30
            + "\n\n.. toctree::\n   :maxdepth: 1\n\n"
            + "".join("   " + page + "\n" for page in pages),
            encoding="utf-8",
        )
    (docs / "index.rst").write_text(
        "API 参考\n========\n\nAPI 从规范源码路径生成。旧路径只保留导入兼容，不重复列出。\n\n.. toctree::\n   :maxdepth: 2\n\n   infrastructure/index\n   algorithms/index\n   applications/index\n",
        encoding="utf-8",
    )
    print("Generated", sum(map(len, groups.values())), "API pages")


if __name__ == "__main__":
    main()
