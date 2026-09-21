"""Toffoli+Clifford+T+QRAM 级别的资源估计。

沿 OriginIR-ext → basis 降低链路的真实发射配方（``basis.mcx`` 与
``basis.controlled_u3``）逐级计数：模块图组合式聚合，``Repeat`` 以符号
计数相乘，不做指数展开。旋转门按角度分类：π/4 整数倍落入精确
Clifford+T（T/SDG/S/Z 计数），其余记为待合成旋转原子，T 开销按
Ross–Selinger 型模型 ``synthesis_t_per_rotation`` 估算（默认
``ceil(3*log2(1/eps))``）。全局相位不计（与 strict 导出一致）。

控制成本与导出器逐行对齐：零值控制位在每个发射行组（原始门 / 调用 /
Repeat 调用行）两侧各翻转一次（2×零位数个 X）；受控单比特门走
``controlled_u3`` 配方（c≥2 时 2(c−1) Toffoli 梯子 + 中间网络），X 型
多控走 ``mcx`` 配方（c≥3 时 2c−3 Toffoli）。

QRAM 查询数按 ``Load`` 节点逐资源计数（每次 Load = 1 次查询），是
独立于门级成本的第一类指标。QRAM 随机写按 ``Store`` 节点逐资源计入
``qram_writes``：存储单元按经典单元建模，随机写不计入门成本。
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field
from typing import cast

from pyqecclang.infrastructure.ir import (
    Adjoint,
    Call,
    Control,
    Instruction,
    Load,
    Primitive,
    Program,
    Repeat,
    Store,
    ValidationError,
)
from pyqecclang.infrastructure.layout import workspace_table
from pyqecclang.infrastructure.validation import validate

_PI = math.pi

CLIFFORD_ATOMS = ("h", "x", "y", "z", "s", "sdg", "cnot", "cz")
"""Clifford 原子名集合；``ResourceEstimate.clifford`` 按这些名字汇总 ``atoms`` 计数。"""


def _rz_exact_counts(k: int) -> Counter[str]:
    """rz(k·π/4) 的精确 Clifford+T 计数（k 按 mod 8 归一）。"""
    k %= 8
    if k > 4:
        counts = _rz_exact_counts(8 - k)
        result: Counter[str] = Counter()
        for atom, n in counts.items():
            result[{"t": "tdg", "s": "sdg", "tdg": "t", "sdg": "s"}.get(atom, atom)] = n
        return result
    return Counter(
        {
            0: {},
            1: {"t": 1},
            2: {"s": 1},
            3: {"s": 1, "t": 1},
            4: {"z": 1},
        }[k]
    )


def _quarter_turns(angle: float) -> int | None:
    """angle 是否为 π/4 的整数倍；是则返回 k（mod 8），否则 None。"""
    k = round(angle / (_PI / 4))
    if abs(angle - k * _PI / 4) < 1e-12:
        return k % 8
    return None


def classify_rz(angle: float) -> tuple[Counter[str], list[tuple[str, float]]]:
    """rz/phase 原子分类：精确 Clifford+T 或待合成旋转。

    Args:
        angle: 旋转角，单位为弧度。

    Returns:
        tuple[Counter[str], list[tuple[str, float]]]: 二元组 ``(精确原子计数,
        待合成旋转列表)``；角度为 π/4 整数倍时计入前者，否则整笔记为
        ``("rz", angle)`` 待合成旋转。
    """
    k = _quarter_turns(angle)
    if k is None:
        return Counter(), [("rz", angle)]
    return _rz_exact_counts(k), []


def classify_ry(angle: float) -> tuple[Counter[str], list[tuple[str, float]]]:
    """ry 原子分类：RY(θ) = S·H·RZ(θ)·H·S†。

    Args:
        angle: 旋转角，单位为弧度。

    Returns:
        tuple[Counter[str], list[tuple[str, float]]]: 二元组 ``(精确原子计数,
        待合成旋转列表)``；角度非 π/4 整数倍时整笔记为 ``("ry", angle)`` 待
        合成旋转，零角度返回两个空容器。
    """
    k = _quarter_turns(angle)
    if k is None:
        return Counter(), [("ry", angle)]
    if k == 0:
        return Counter(), []
    counts = Counter({"s": 1, "h": 2, "sdg": 1})
    counts += _rz_exact_counts(k)
    return counts, []


def classify_u3(
    theta: float, phi: float, lam: float
) -> tuple[Counter[str], list[tuple[str, float]]]:
    """U3(θ,φ,λ) = RZ(φ)·RY(θ)·RZ(λ)（丢弃全局相位，与 strict 导出一致）。

    常见门 H/X/Y 先按精确单原子匹配，再走 Euler 分解。

    Args:
        theta: Y 轴欧拉角，单位为弧度。
        phi: 左侧 Z 轴欧拉角，单位为弧度。
        lam: 右侧 Z 轴欧拉角，单位为弧度。

    Returns:
        tuple[Counter[str], list[tuple[str, float]]]: 二元组 ``(精确原子计数,
        待合成旋转列表)``；H/X/Y 特例匹配单原子，其余按 Euler 分解逐项汇总。
    """

    def close(a: float, b: float) -> bool:
        """判断两个浮点数在容差内相等。"""
        return abs(a - b) < 1e-12

    if close(theta, _PI / 2) and close(phi, 0.0) and close(lam, _PI):
        return Counter({"h": 1}), []
    if close(theta, _PI) and close(phi, 0.0) and close(lam, _PI):
        return Counter({"x": 1}), []
    if close(theta, _PI) and close(phi, _PI / 2) and close(lam, _PI / 2):
        return Counter({"y": 1}), []
    counts: Counter[str]
    rotations: list[tuple[str, float]]
    counts, rotations = Counter(), []
    for counts_part, rotations_part in (
        classify_rz(lam),
        classify_ry(theta),
        classify_rz(phi),
    ):
        counts += counts_part
        rotations += rotations_part
    return counts, rotations


def mcx_counts(n_controls: int) -> Counter[str]:
    """basis.mcx 配方（X 型多控）。

    Args:
        n_controls: 有效控制位个数，取非负整数。

    Returns:
        Counter[str]: 多控 X 的原子计数；c 为 0/1/2 时分别为单 X、
        H-CZ-H 与单个 Toffoli，c≥3 时为 2c−3 个 Toffoli。
    """
    if n_controls <= 0:
        return Counter({"x": 1})
    if n_controls == 1:
        return Counter({"h": 2, "cz": 1})
    if n_controls == 2:
        return Counter({"toffoli": 1})
    return Counter({"toffoli": 2 * n_controls - 3})


def _controlled_u3_counts(
    n_controls: int, theta: float, phi: float, lam: float, global_angle: float = 0.0
) -> tuple[Counter[str], list[tuple[str, float]]]:
    """basis.controlled_u3 配方的计数；c≥2 时含 2(c−1) Toffoli 梯子。"""
    counts: Counter[str]
    rotations: list[tuple[str, float]]
    counts, rotations = Counter(), []
    if n_controls == 0:
        part, rot = classify_u3(theta, phi, lam)
        counts += part
        rotations += rot
        if global_angle:
            for _ in range(2):
                part, rot = classify_rz(global_angle)
                counts += part
                rotations += rot
            counts["x"] += 2
        return counts, rotations
    if n_controls >= 2:
        counts["toffoli"] += 2 * (n_controls - 1)
    counts += Counter({"h": 4, "cz": 2})
    for angle in ((lam + phi) / 2 + global_angle, (lam - phi) / 2):
        part, rot = classify_rz(angle)
        counts += part
        rotations += rot
    for theta_part, phi_part, lam_part in (
        (-theta / 2, 0.0, -(phi + lam) / 2),
        (theta / 2, phi, 0.0),
    ):
        part, rot = classify_u3(theta_part, phi_part, lam_part)
        counts += part
        rotations += rot
    return counts, rotations


# basis.py 中单比特门的 (θ, φ, λ, 全局相位) 参数表
_GATE_PARAMETERS = {
    "h": (_PI / 2, 0.0, _PI, 0.0),
    "y": (_PI, _PI / 2, _PI / 2, 0.0),
    "z": (0.0, 0.0, _PI, 0.0),
    "s": (0.0, 0.0, _PI / 2, 0.0),
    "t": (0.0, 0.0, _PI / 4, 0.0),
}


def _gate_counts(
    op: str, angle: float | None, n_controls: int
) -> tuple[Counter[str], list[tuple[str, float]]]:
    """单个 1q 门在 n_controls 个有效控制下的计数。"""
    if op == "x":
        return mcx_counts(n_controls), []
    if op in _GATE_PARAMETERS:
        theta, phi, lam, glob = _GATE_PARAMETERS[op]
        return _controlled_u3_counts(n_controls, theta, phi, lam, glob)
    if op == "phase":
        return _controlled_u3_counts(n_controls, 0.0, 0.0, cast(float, angle), 0.0)
    if op == "ry":
        return _controlled_u3_counts(n_controls, cast(float, angle), 0.0, 0.0, 0.0)
    if op == "rx":
        return _controlled_u3_counts(n_controls, cast(float, angle), -_PI / 2, _PI / 2, 0.0)
    if op == "rz":
        return _controlled_u3_counts(n_controls, 0.0, 0.0, cast(float, angle), -cast(float, angle) / 2)
    raise ValidationError(f"资源估计不支持的原始门：{op}")


@dataclass
class ResourceEstimate:
    "Toffoli+Clifford+T+QRAM 级别的资源台账。"

    qubits: int
    atoms: Counter = field(default_factory=Counter)
    rotations: list = field(default_factory=list)
    qram_queries: Counter = field(default_factory=Counter)
    qram_writes: Counter = field(default_factory=Counter)
    mcx_ancilla: int = 0

    @property
    def toffoli(self) -> int:
        """``atoms`` 中的 Toffoli 门总数。"""
        return self.atoms.get("toffoli", 0)

    @property
    def t_exact(self) -> int:
        """精确落入 Clifford+T 的 T 与 TDG 门总数（不含待合成旋转的 T 开销）。"""
        return self.atoms.get("t", 0) + self.atoms.get("tdg", 0)

    @property
    def clifford(self) -> int:
        """``atoms`` 中 ``CLIFFORD_ATOMS`` 所列原子的计数总和。"""
        return sum(self.atoms.get(name, 0) for name in CLIFFORD_ATOMS)

    @property
    def qram_total(self) -> int:
        """``qram_queries`` 中全部资源的查询次数总和。"""
        return sum(self.qram_queries.values())

    @property
    def qram_write_total(self) -> int:
        """``qram_writes`` 中全部资源的随机写次数总和。"""
        return sum(self.qram_writes.values())

    @property
    def gate_total(self) -> int:
        """``atoms`` 中全部原子计数的总和（不含 QRAM 查询与写入）。"""
        return sum(self.atoms.values())

    def synthesis_t_per_rotation(self, epsilon: float = 1e-10) -> int:
        """Ross–Selinger 型前导项模型；可整体替换。

        Args:
            epsilon: 单个旋转的合成精度，取值范围为 (0, 1)。

        Returns:
            int: 每个待合成旋转的 T 门开销，即 ``ceil(3*log2(1/epsilon))``，至少为 1。
        """
        return max(1, math.ceil(3 * math.log2(1 / epsilon)))

    def t_synthesis(self, epsilon: float = 1e-10) -> int:
        """待合成旋转的总 T 开销：旋转数乘以 ``synthesis_t_per_rotation(epsilon)``。

        Args:
            epsilon: 单个旋转的合成精度。

        Returns:
            int: 全部待合成旋转的 T 门开销估计。
        """
        return len(self.rotations) * self.synthesis_t_per_rotation(epsilon)

    def t_total(self, epsilon: float = 1e-10) -> int:
        """精确 T 数与待合成旋转的 T 开销之和（``t_exact + t_synthesis``）。

        Args:
            epsilon: 单个旋转的合成精度。

        Returns:
            int: 精确 T/TDG 计数加上合成 T 开销的总 T 门数。
        """
        return self.t_exact + self.t_synthesis(epsilon)

    def to_dict(self, epsilon: float = 1e-10) -> dict[str, object]:
        """导出 JSON 友好的扁平资源台账。

        Args:
            epsilon: 旋转合成精度，同时作为 ``epsilon`` 键记入结果。

        Returns:
            dict: 包含量子位数与 ``mcx_ancilla``、Toffoli/Clifford/精确 T
            计数、待合成旋转总数及按轴统计、T 开销估计、逐原子 ``atoms``
            明细、逐资源 QRAM 查询/写入计数与总和，以及 ``gate_total``。
        """
        return {
            "qubits": self.qubits,
            "mcx_ancilla": self.mcx_ancilla,
            "toffoli": self.toffoli,
            "clifford": self.clifford,
            "t_exact": self.t_exact,
            "rotations": len(self.rotations),
            "rotation_axes": dict(Counter(axis for axis, _ in self.rotations)),
            "t_synthesis": self.t_synthesis(epsilon),
            "t_total": self.t_total(epsilon),
            "atoms": dict(sorted(self.atoms.items())),
            "qram_queries": dict(sorted(self.qram_queries.items())),
            "qram_total": self.qram_total,
            "qram_writes": dict(sorted(self.qram_writes.items())),
            "qram_write_total": self.qram_write_total,
            "gate_total": self.gate_total,
            "epsilon": epsilon,
        }


def estimate_resources(program: Program, *, require_closed: bool = True) -> ResourceEstimate:
    """组合式资源估计：按 (module, 控制数) 记忆化，Repeat 符号相乘。

    Args:
        program: 待估计的 RIR 程序；模块实现须封闭，除非放宽校验。
        require_closed: 为假时允许保留开放声明模块，默认要求全部封闭。

    Returns:
        ResourceEstimate: 汇总后的资源台账；量子位数含工作空间，最大控制数
        决定 ``mcx_ancilla``。
    """
    program = validate(program, require_closed=require_closed)
    modules = program.module_map
    workspace = workspace_table(program)
    memo: dict[
        tuple[str, int],
        tuple[Counter[str], list[tuple[str, float]], Counter[str], Counter[str]],
    ] = {}
    max_controls = 0

    def merge(
        total: tuple[Counter[str], list[tuple[str, float]], Counter[str], Counter[str]],
        part: tuple[Counter[str], list[tuple[str, float]], Counter[str], Counter[str]],
        factor: int = 1,
    ) -> None:
        """按重复倍数把分段成本累加进总成本台账。"""
        counts, rotations, queries, writes = part
        for atom, n in counts.items():
            total[0][atom] += n * factor
        total[1].extend(rotations * factor)
        for resource, n in queries.items():
            total[2][resource] += n * factor
        for resource, n in writes.items():
            total[3][resource] += n * factor

    def primitive_cost(
        node: Primitive, n_controls: int
    ) -> tuple[Counter[str], list[tuple[str, float]], Counter[str], Counter[str]]:
        """计算单条基元指令在给定控制数下的门级成本。"""
        nonlocal max_controls
        counts: Counter[str]
        rotations: list[tuple[str, float]]
        queries: Counter[str]
        writes: Counter[str]
        counts, rotations, queries, writes = Counter(), [], Counter(), Counter()
        widths = [ref.width for ref in node.operands]
        if node.op in {"xor", "swap"}:
            for _ in range(widths[0]):
                counts += mcx_counts(n_controls + 1)
                if node.op == "swap":
                    counts += mcx_counts(n_controls + 1)
                    counts += mcx_counts(n_controls + 1)
            max_controls = max(max_controls, n_controls + 1)
        elif node.op == "add_const":
            value = cast(int, node.value) % (1 << widths[0])
            for offset in range(widths[0]):
                if (value >> offset) & 1:
                    for i in reversed(range(offset + 1, widths[0])):
                        counts += mcx_counts(n_controls + (i - offset))
                        max_controls = max(max_controls, n_controls + (i - offset))
                    counts += mcx_counts(n_controls)
        elif node.op == "gphase":
            theta = cast(float, node.angle)
            if n_controls:
                part, rot = _controlled_u3_counts(n_controls - 1, 0.0, 0.0, theta, 0.0)
                counts += part
                rotations += rot
                max_controls = max(max_controls, n_controls - 1)
            else:
                for _ in range(2):
                    part, rot = classify_rz(theta)
                    counts += part
                    rotations += rot
                counts["x"] += 2
        else:
            for _ in range(widths[0]):
                part, rot = _gate_counts(node.op, node.angle, n_controls)
                counts += part
                rotations += rot
            max_controls = max(max_controls, n_controls)
        return counts, rotations, queries, writes

    def body_cost(
        nodes: tuple[Instruction, ...], n_controls: int, n_zeros: int
    ) -> tuple[Counter[str], list[tuple[str, float]], Counter[str], Counter[str]]:
        """聚合指令体成本，处理零值控制位翻转、重复与模块调用。"""
        counts: Counter[str]
        rotations: list[tuple[str, float]]
        queries: Counter[str]
        writes: Counter[str]
        counts, rotations, queries, writes = Counter(), [], Counter(), Counter()
        for node in nodes:
            if isinstance(node, Primitive):
                if n_zeros and not (
                    node.op == "add_const"
                    and not cast(int, node.value) % (1 << node.operands[0].width)
                ):
                    counts["x"] += 2 * n_zeros
                merge((counts, rotations, queries, writes), primitive_cost(node, n_controls))
            elif isinstance(node, Load):
                queries[node.resource] += 1
            elif isinstance(node, Store):
                writes[node.resource] += 1
            elif isinstance(node, Call):
                if n_zeros:
                    counts["x"] += 2 * n_zeros
                merge((counts, rotations, queries, writes), module_cost(node.module, n_controls))
            elif isinstance(node, Repeat):
                if node.count and node.body:
                    if n_zeros:
                        counts["x"] += 2 * n_zeros
                    merge(
                        (counts, rotations, queries, writes),
                        body_cost(node.body, n_controls, 0),
                        node.count,
                    )
            elif isinstance(node, Control):
                merge(
                    (counts, rotations, queries, writes),
                    body_cost(
                        node.body,
                        n_controls + node.register.width,
                        n_zeros + node.register.width - bin(node.value).count("1"),
                    ),
                )
            elif isinstance(node, Adjoint):
                merge((counts, rotations, queries, writes), body_cost(node.body, n_controls, n_zeros))
        return counts, rotations, queries, writes

    def module_cost(
        name: str, n_controls: int
    ) -> tuple[Counter[str], list[tuple[str, float]], Counter[str], Counter[str]]:
        """按 (模块名, 控制数) 记忆化地计算模块体成本。"""
        key = (name, n_controls)
        if key in memo:
            return memo[key]
        module = modules[name]
        if module.body is None:
            raise ValidationError(f"资源估计需要封闭程序：模块 {name} 无实现")
        memo[key] = body_cost(module.body, n_controls, 0)
        return memo[key]

    counts, rotations, queries, writes = module_cost(program.entry, 0)
    return ResourceEstimate(
        qubits=sum(r.type.width for r in program.main.registers) + workspace[program.entry],
        atoms=counts,
        rotations=rotations,
        qram_queries=queries,
        qram_writes=writes,
        mcx_ancilla=max(0, max_controls - 1),
    )
