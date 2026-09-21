"""论文级数值验证的共享设施。

verify_*.py 脚本使用这里封装的后端路径与比较器，产出 out/verification/<group>.json。
运行解释器必须同时安装 pysparq 与 uniqc（真实后端，不使用替身）：

    PYTHONPATH=src /path/to/python tests/verification/verify_<group>.py
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
import platform
import time
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ARTIFACTS = ROOT / "out" / "verification"

# ---------------------------------------------------------------------------
# 后端路径
# ---------------------------------------------------------------------------


def reference(program, memory=None, **kwargs):
    """pyqecclang 内置参考执行器（字典稀疏态）。"""
    from pyqecclang import simulate

    return dict(simulate(program, memory, **kwargs).amplitudes)


def rir_pysparq(program, memory=None, **kwargs):
    """PySparQ 原生 RIR 解释器（QRAM-Simulator 仓 pysparq.rir）。"""
    from pyqecclang import run_pysparq_rir

    return dict(run_pysparq_rir(program, memory, **kwargs).amplitudes)


def adapter_pysparq(program, memory=None, **kwargs):
    """pyqecclang 的 PySparQ 事件适配器。"""
    from pyqecclang import run_pysparq

    return dict(run_pysparq(program, memory, **kwargs).amplitudes)


def originir_ext(program, memory=None, **kwargs):
    """OriginIR-ext + UnifiedQuantum 全振幅态向量（返回 list[complex]）。"""
    from pyqecclang import run_originir

    return [complex(value) for value in run_originir(program, memory, **kwargs)]


def originir_unitary(program):
    """OriginIR-ext 线路经 UniQC ``Circuit.to_matrix`` 的全幺正矩阵。

    只适用于不含 QRAM 资源的程序；矩阵覆盖全部量子位（含工作区），
    量子位 0 为最低位。
    """
    import numpy as np
    from uniqc.circuit_builder import Circuit

    from pyqecclang import export_originir

    text = export_originir(program).text
    circuit = Circuit.from_originir_ext(text)
    return np.asarray(circuit.to_matrix(), dtype=complex)


def amplitudes_to_statevector(amplitudes, widths):
    vector = [0j] * (1 << sum(widths))
    for values, amplitude in amplitudes.items():
        index = 0
        offset = 0
        for value, w in zip(values, widths, strict=True):
            index |= value << offset
            offset += w
        vector[index] = amplitude
    return vector


# ---------------------------------------------------------------------------
# 比较器
# ---------------------------------------------------------------------------


def amplitude_error(actual, expected):
    """两个字典稀疏态（寄存器元组→振幅）的最大绝对偏差。"""
    keys = set(actual) | set(expected)
    if not keys:
        return 0.0
    return max(abs(actual.get(k, 0j) - expected.get(k, 0j)) for k in keys)


def statevector_error(actual, expected):
    if len(actual) != len(expected):
        return math.inf
    if not actual:
        return 0.0
    return max(abs(a - b) for a, b in zip(actual, expected, strict=True))


def fidelity(actual, expected):
    """两态向量的保真度 |<a|b>|^2。"""
    overlap = sum(a.conjugate() * b for a, b in zip(actual, expected, strict=True))
    return abs(overlap) ** 2


def tvd(actual_probs, expected_probs):
    """两个概率字典的全变差距离。"""
    keys = set(actual_probs) | set(expected_probs)
    return 0.5 * sum(
        abs(actual_probs.get(k, 0.0) - expected_probs.get(k, 0.0)) for k in keys
    )


def probabilities(amplitudes):
    return {key: abs(value) ** 2 for key, value in amplitudes.items()}


def effective_block(unitary, data_qubits):
    """提取 data 低位、工作区为 0 的有效算子块与泄漏上界。

    OriginIR 导出中入口寄存器占据低位量子位、工作区紧随其后；幺正矩阵
    的低 2^data_qubits 维块即工作区复净时的有效算子。
    """
    import numpy as np

    dim = 1 << data_qubits
    block = unitary[:dim, :dim]
    leakage = float(np.abs(unitary[dim:, :dim]).max()) if unitary.shape[0] > dim else 0.0
    return block, leakage


# ---------------------------------------------------------------------------
# 驱动程序构造
# ---------------------------------------------------------------------------


def _driver(operation, name):
    from pyqecclang import Builder

    module = operation.module
    return Builder(name, {r.name: r.type for r in module.registers})


def basis_program(operation, initial, *, name="verify_basis"):
    """把 entry 寄存器置为 classical 初态后调用 operation。"""
    b = _driver(operation, name)
    for key, value in initial.items():
        for bit in range(b[key].width):
            if (value >> bit) & 1:
                b.x(b[key][bit])
    b.call(operation, **{r.name: b[r.name] for r in operation.module.registers})
    return b.finish().program()


def superposition_program(operation, registers, *, name="verify_superposition"):
    """对指定寄存器加 H 后调用 operation（其余寄存器保持 |0>）。"""
    b = _driver(operation, name)
    for key in registers:
        b.h(b[key])
    b.call(operation, **{r.name: b[r.name] for r in operation.module.registers})
    return b.finish().program()


# ---------------------------------------------------------------------------
# 报告与产物
# ---------------------------------------------------------------------------


@dataclass
class Case:
    name: str
    paths: list[str]
    parameters: dict
    metrics: dict
    criterion: str
    passed: bool


@dataclass
class Report:
    group: str
    summary: str
    cases: list = field(default_factory=list)
    started: float = field(default_factory=time.time)

    def case(self, name, *, paths, parameters, metrics, criterion, passed):
        self.cases.append(
            Case(
                name=name,
                paths=list(paths),
                parameters=dict(parameters),
                metrics=dict(metrics),
                criterion=criterion,
                passed=bool(passed),
            )
        )
        status = "PASS" if passed else "FAIL"
        print(f"[{status}] {self.group}/{name} {metrics}", flush=True)
        if not passed:
            raise AssertionError(f"验证失败：{self.group}/{name}（判据：{criterion}）")

    def write(self):
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        path = ARTIFACTS / f"{self.group}.json"
        payload = {
            "python": platform.python_version(),
            "source_sha256": hashlib.sha256(b"".join(
                str(p.relative_to(ROOT)).encode() + b"\0" + p.read_bytes()
                for p in sorted((ROOT / "src").rglob("*.py"))
            )).hexdigest(),
            "verification_sha256": hashlib.sha256(b"".join(
                p.name.encode() + b"\0" + p.read_bytes()
                for p in sorted((ROOT / "tests/verification").glob("*.py"))
            )).hexdigest(),
            "backend_versions": {
                dist.metadata["Name"]: dist.version
                for dist in importlib.metadata.distributions()
                if dist.metadata["Name"].lower() in {"pysparq", "unified-quantum", "uniqc-cppsimulator", "numpy", "scipy"}
            },
            "group": self.group,
            "summary": self.summary,
            "elapsed_seconds": round(time.time() - self.started, 3),
            "cases": [
                {
                    "name": case.name,
                    "paths": case.paths,
                    "parameters": case.parameters,
                    "metrics": case.metrics,
                    "criterion": case.criterion,
                    "passed": case.passed,
                }
                for case in self.cases
            ],
        }
        path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
        print(f"[artifact] {path}", flush=True)
        return path


def sampled_inputs(width, *, exhaustive_below=8, samples=17):
    """比特数扫描：小宽度穷举，大宽度取边界值与确定性伪随机样本。"""
    if width <= exhaustive_below:
        return list(range(1 << width)), True
    limit = 1 << width
    values = {0, 1, limit // 2, limit - 2, limit - 1}
    state = 0x9E3779B97F4A7C15
    while len(values) < samples:
        state = (state * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
        values.add(state % limit)
    return sorted(values), False
