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
独立于门级成本的第一类指标。
"""""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field

from pyqecclang.infrastructure.ir import (
    Adjoint,
    Call,
    Control,
    Load,
    Primitive,
    Repeat,
    ValidationError,
)
from pyqecclang.infrastructure.layout import workspace_table
from pyqecclang.infrastructure.validation import validate

_PI = math.pi

CLIFFORD_ATOMS = ("h", "x", "y", "z", "s", "sdg", "cnot", "cz")


def _rz_exact_counts(k):
    """rz(k·π/4) 的精确 Clifford+T 计数（k 按 mod 8 归一）。"""
    k %= 8
    if k > 4:
        counts = _rz_exact_counts(8 - k)
        result = Counter()
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


def _quarter_turns(angle):
    """angle 是否为 π/4 的整数倍；是则返回 k（mod 8），否则 None。"""
    k = round(angle / (_PI / 4))
    if abs(angle - k * _PI / 4) < 1e-12:
        return k % 8
    return None


def classify_rz(angle):
    """rz/phase 原子分类：精确 Clifford+T 或待合成旋转。"""
    k = _quarter_turns(angle)
    if k is None:
        return Counter(), [("rz", angle)]
    return _rz_exact_counts(k), []


def classify_ry(angle):
    """ry 原子分类：RY(θ) = S·H·RZ(θ)·H·S†。"""
    k = _quarter_turns(angle)
    if k is None:
        return Counter(), [("ry", angle)]
    if k == 0:
        return Counter(), []
    counts = Counter({"s": 1, "h": 2, "sdg": 1})
    counts += _rz_exact_counts(k)
    return counts, []


def classify_u3(theta, phi, lam):
    """U3(θ,φ,λ) = RZ(φ)·RY(θ)·RZ(λ)（丢弃全局相位，与 strict 导出一致）。

    常见门 H/X/Y 先按精确单原子匹配，再走 Euler 分解。
    """

    def close(a, b):
        return abs(a - b) < 1e-12

    if close(theta, _PI / 2) and close(phi, 0.0) and close(lam, _PI):
        return Counter({"h": 1}), []
    if close(theta, _PI) and close(phi, 0.0) and close(lam, _PI):
        return Counter({"x": 1}), []
    if close(theta, _PI) and close(phi, _PI / 2) and close(lam, _PI / 2):
        return Counter({"y": 1}), []
    counts, rotations = Counter(), []
    for counts_part, rotations_part in (
        classify_rz(lam),
        classify_ry(theta),
        classify_rz(phi),
    ):
        counts += counts_part
        rotations += rotations_part
    return counts, rotations


def mcx_counts(n_controls):
    """basis.mcx 配方（X 型多控）。"""
    if n_controls <= 0:
        return Counter({"x": 1})
    if n_controls == 1:
        return Counter({"h": 2, "cz": 1})
    if n_controls == 2:
        return Counter({"toffoli": 1})
    return Counter({"toffoli": 2 * n_controls - 3})


def _controlled_u3_counts(n_controls, theta, phi, lam, global_angle=0.0):
    """basis.controlled_u3 配方的计数；c≥2 时含 2(c−1) Toffoli 梯子。"""
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


def _gate_counts(op, angle, n_controls):
    """单个 1q 门在 n_controls 个有效控制下的计数。"""
    if op == "x":
        return mcx_counts(n_controls), []
    if op in _GATE_PARAMETERS:
        theta, phi, lam, glob = _GATE_PARAMETERS[op]
        return _controlled_u3_counts(n_controls, theta, phi, lam, glob)
    if op == "phase":
        return _controlled_u3_counts(n_controls, 0.0, 0.0, angle, 0.0)
    if op == "ry":
        return _controlled_u3_counts(n_controls, angle, 0.0, 0.0, 0.0)
    if op == "rx":
        return _controlled_u3_counts(n_controls, angle, -_PI / 2, _PI / 2, 0.0)
    if op == "rz":
        return _controlled_u3_counts(n_controls, 0.0, 0.0, angle, -angle / 2)
    raise ValidationError(f"资源估计不支持的原始门：{op}")


@dataclass
class ResourceEstimate:
    "Toffoli+Clifford+T+QRAM 级别的资源台账。"

    qubits: int
    atoms: Counter = field(default_factory=Counter)
    rotations: list = field(default_factory=list)
    qram_queries: Counter = field(default_factory=Counter)
    mcx_ancilla: int = 0

    @property
    def toffoli(self):
        return self.atoms.get("toffoli", 0)

    @property
    def t_exact(self):
        return self.atoms.get("t", 0) + self.atoms.get("tdg", 0)

    @property
    def clifford(self):
        return sum(self.atoms.get(name, 0) for name in CLIFFORD_ATOMS)

    @property
    def qram_total(self):
        return sum(self.qram_queries.values())

    @property
    def gate_total(self):
        return sum(self.atoms.values())

    def synthesis_t_per_rotation(self, epsilon=1e-10):
        "Ross–Selinger 型前导项模型；可整体替换。"
        return max(1, math.ceil(3 * math.log2(1 / epsilon)))

    def t_synthesis(self, epsilon=1e-10):
        return len(self.rotations) * self.synthesis_t_per_rotation(epsilon)

    def t_total(self, epsilon=1e-10):
        return self.t_exact + self.t_synthesis(epsilon)

    def to_dict(self, epsilon=1e-10):
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
            "gate_total": self.gate_total,
            "epsilon": epsilon,
        }


def estimate_resources(program, *, require_closed=True):
    "组合式资源估计：按 (module, 控制数) 记忆化，Repeat 符号相乘。"
    program = validate(program, require_closed=require_closed)
    modules = program.module_map
    workspace = workspace_table(program)
    memo = {}
    max_controls = 0

    def merge(total, part, factor=1):
        counts, rotations, queries = part
        for atom, n in counts.items():
            total[0][atom] += n * factor
        total[1].extend(rotations * factor)
        for resource, n in queries.items():
            total[2][resource] += n * factor

    def primitive_cost(node, n_controls):
        nonlocal max_controls
        counts, rotations, queries = Counter(), [], Counter()
        widths = [ref.width for ref in node.operands]
        if node.op in {"xor", "swap"}:
            for _ in range(widths[0]):
                counts += mcx_counts(n_controls + 1)
                if node.op == "swap":
                    counts += mcx_counts(n_controls + 1)
                    counts += mcx_counts(n_controls + 1)
            max_controls = max(max_controls, n_controls + 1)
        elif node.op == "add_const":
            value = node.value % (1 << widths[0])
            for offset in range(widths[0]):
                if (value >> offset) & 1:
                    for i in reversed(range(offset + 1, widths[0])):
                        counts += mcx_counts(n_controls + (i - offset))
                        max_controls = max(max_controls, n_controls + (i - offset))
                    counts += mcx_counts(n_controls)
        elif node.op == "gphase":
            theta = node.angle
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
        return counts, rotations, queries

    def body_cost(nodes, n_controls, n_zeros):
        counts, rotations, queries = Counter(), [], Counter()
        for node in nodes:
            if isinstance(node, Primitive):
                if n_zeros and not (node.op == "add_const" and not node.value % (1 << node.operands[0].width)):
                    counts["x"] += 2 * n_zeros
                merge((counts, rotations, queries), primitive_cost(node, n_controls))
            elif isinstance(node, Load):
                queries[node.resource] += 1
            elif isinstance(node, Call):
                if n_zeros:
                    counts["x"] += 2 * n_zeros
                merge((counts, rotations, queries), module_cost(node.module, n_controls))
            elif isinstance(node, Repeat):
                if node.count and node.body:
                    if n_zeros:
                        counts["x"] += 2 * n_zeros
                    merge(
                        (counts, rotations, queries),
                        body_cost(node.body, n_controls, 0),
                        node.count,
                    )
            elif isinstance(node, Control):
                merge(
                    (counts, rotations, queries),
                    body_cost(
                        node.body,
                        n_controls + node.register.width,
                        n_zeros + node.register.width - bin(node.value).count("1"),
                    ),
                )
            elif isinstance(node, Adjoint):
                merge((counts, rotations, queries), body_cost(node.body, n_controls, n_zeros))
        return counts, rotations, queries

    def module_cost(name, n_controls):
        key = (name, n_controls)
        if key in memo:
            return memo[key]
        module = modules[name]
        if module.body is None:
            raise ValidationError(f"资源估计需要封闭程序：模块 {name} 无实现")
        memo[key] = body_cost(module.body, n_controls, 0)
        return memo[key]

    counts, rotations, queries = module_cost(program.entry, 0)
    return ResourceEstimate(
        qubits=sum(r.type.width for r in program.main.registers) + workspace[program.entry],
        atoms=counts,
        rotations=rotations,
        qram_queries=queries,
        mcx_ancilla=max(0, max_controls - 1),
    )
