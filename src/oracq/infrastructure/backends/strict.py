"""Toffoli+Clifford+T+QRAM 严格网表导出。

在 ``basis.lower_toffoli_u3_cz`` 的产物上再做一遍 U3 分类：π/4 整数倍
角度的旋转落入精确 Clifford+T 原子（H/S/SDG/T/TDG/Z/X/Y），其余保留为
``RZ``/``RY`` 待合成旋转原子；模块调用、Repeat 符号结构与 QRAM 行保持
不展开。原子集合：TOFFOLI、CZ、H、S、SDG、T、TDG、X、Y、Z、RZ、RY、
QRAM。角度分类与 ``estimate`` 模块共享同一实现，保证网表计数与资源
估计逐项一致（测试中有逐行对拍）。"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from oracq.infrastructure.backends.basis import lower_toffoli_u3_cz
from oracq.infrastructure.backends.originir import OriginIRArtifact, export_originir
from oracq.infrastructure.estimate import classify_ry, classify_rz, classify_u3
from oracq.infrastructure.ir import Program, ValidationError

_ATOM_NAMES = {
    "t": "T",
    "tdg": "TDG",
    "s": "S",
    "sdg": "SDG",
    "h": "H",
    "x": "X",
    "y": "Y",
    "z": "Z",
}


@dataclass(frozen=True)
class StrictArtifact:
    "严格网表：文本 + 逐原子计数 + 布局信息。"

    text: str
    counts: Counter[str]
    registers: dict[str, tuple[int, ...]]
    resources: dict[str, str]
    workspace_qubits: tuple[int, ...]
    qram_queries: Counter[str]
    qram_writes: Counter[str] = field(default_factory=Counter)


def _emit_rz(lines: list[str], counts: Counter[str], target: str, angle: float) -> None:
    """按 ``classify_rz`` 的分类向网表发射 RZ 角度的原子序列。"""
    atoms, rotations = classify_rz(angle)
    for atom, n in sorted(atoms.items()):
        lines.extend(f"{_ATOM_NAMES[atom]} {target}" for _ in range(n))
    for axis, value in rotations:
        lines.append(f"{axis.upper()} {target}, ({value!r})")
    counts.update(atoms)
    counts.update(axis for axis, _ in rotations)


def _emit_ry(lines: list[str], counts: Counter[str], target: str, angle: float) -> None:
    """按 ``classify_ry`` 的分类向网表发射 RY 角度的原子序列。"""
    atoms, rotations = classify_ry(angle)
    if rotations:
        lines.append(f"RY {target}, ({rotations[0][1]!r})")
        counts["ry"] += 1
        return
    if not atoms:
        return  # 恒等：RY(0) 不发射
    # RY(θ) = S·H·RZ(θ)·H·S†：S/H 环绕 + rz 精确序列
    lines.append(f"S {target}")
    lines.append(f"H {target}")
    counts.update({"s": 1, "h": 2, "sdg": 1})
    _emit_rz(lines, counts, target, angle)
    lines.append(f"H {target}")
    lines.append(f"SDG {target}")


def _emit_u3(
    lines: list[str], counts: Counter[str], target: str, theta: float, phi: float, lam: float
) -> None:
    """按 ``classify_u3`` 的分类向网表发射 U3 角度的原子序列。"""
    atoms, rotations = classify_u3(theta, phi, lam)
    if not rotations and set(atoms) <= {"h", "x", "y"} and sum(atoms.values()) == 1:
        atom = next(iter(atoms))
        lines.append(f"{_ATOM_NAMES[atom]} {target}")
        counts[atom] += 1
        return
    # 与 classify_u3 相同的 Euler 顺序：RZ(λ) → RY(θ) → RZ(φ)
    _emit_rz(lines, counts, target, lam)
    _emit_ry(lines, counts, target, theta)
    _emit_rz(lines, counts, target, phi)


def lower_strict(artifact: OriginIRArtifact) -> StrictArtifact:
    """把 Toffoli+U3+CZ 网表进一步降低为严格原子网表。

    ``U3`` 行按与 ``estimate`` 共享的角度分类改写：π/4 整数倍角度发射
    精确 Clifford+T 原子，其余保留为 ``RZ``/``RY`` 待合成旋转原子；其余
    行原样保留，同时逐行累计原子计数、QRAM 查询（``ram_`` 行）与随机写
    （``QRAMWRITE`` 行）。

    Args:
        artifact: ``lower_toffoli_u3_cz`` 产出的网表。

    Returns:
        StrictArtifact: 严格网表文本、逐原子计数与原布局信息。

    Raises:
        ValidationError: ``U3`` 行的角度参数不是三个。
    """
    lines: list[str] = []
    counts: Counter[str] = Counter()
    qram: Counter[str] = Counter()
    writes: Counter[str] = Counter()
    for line in artifact.text.splitlines():
        if line.startswith("U3 "):
            head, _, tail = line.partition("(")
            target = head.split(None, 1)[1].strip().rstrip(",")
            angles = [float(part) for part in tail.rstrip(")").split(",")]
            if len(angles) != 3:
                raise ValidationError("strict  lowering 遇到非法 U3 行：" + line)
            _emit_u3(lines, counts, target, *angles)
            continue
        keyword = line.split(" ", 1)[0]
        lowered = keyword.lower()
        if lowered in {"toffoli", "cz"}:
            counts[lowered] += 1
        elif keyword.startswith("ram_"):
            qram[keyword] += 1
        elif line.startswith("QRAMWRITE "):
            writes[line.split()[1]] += 1
        lines.append(line)
    return StrictArtifact(
        "\n".join(lines) + "\n",
        counts,
        artifact.registers,
        artifact.resources,
        artifact.workspace_qubits,
        qram,
        writes,
    )


def export_strict(program: Program) -> StrictArtifact:
    """编译到 Toffoli+Clifford+T(+待合成旋转)+QRAM 级别的模块保持网表。

    Args:
        program: 待导出的闭合 RIR 程序。

    Returns:
        StrictArtifact: 网表文本，连同 Toffoli/CZ 门计数、QRAM 查询与
        写入统计以及寄存器和资源映射。
    """
    return lower_strict(lower_toffoli_u3_cz(export_originir(program)))
